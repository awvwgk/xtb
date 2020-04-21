! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

!> Implementation of a generalized Born solvation model with surface area
!> contributions.
module xtb_solv_gbsa
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : symv, dot
   use xtb_mctc_convert, only : aatoau
   use xtb_mctc_la, only : contract
   use xtb_param_vdwradd3, only : vanDerWaalsRadD3
   use xtb_solv_born, only : TBornIntegrator, init_ => init
   use xtb_solv_gbsaparam, only : TGBSAData
   use xtb_solv_sasa, only : TSurfaceIntegrator, init_ => init
   use xtb_solv_state, only : getStateShift
   use xtb_type_environment, only : TEnvironment
   use xtb_type_neighbourlist, only : TNeighbourList
   use xtb_type_solvation, only : TSolvation
   implicit none
   private

   public :: TGeneralizedBorn, gbKernel, init


   !> Possible kernels for the generalized Born model
   type :: TGBKernelEnum

      !> Classical Still kernel
      integer :: still = 1

      !> P16 kernel by Lange (JCTC 2012, 8, 1999-2011)
      integer :: p16 = 2

   end type TGBKernelEnum

   !> Actual enumerator for the generalized Born kernels
   type(TGBKernelEnum), parameter :: gbKernel = TGBKernelEnum()


   !> Generalized Born solvation model with surface area contributions
   type, extends(TSolvation) :: TGeneralizedBorn

      !> Number of atoms
      integer :: nAtom

      !> Generalized Born kernel
      integer :: gbKernel

      !> Evaluator for Born radii
      type(TBornIntegrator) :: born

      !> Evaluator for Surface area
      type(TSurfaceIntegrator) :: surface

      !> Dielectric scaling
      real(wp) :: kEps

      !> Shift for free energy calculations
      real(wp) :: freeEnergyShift

      !> Surface tension scaling for each element
      real(wp), allocatable :: surfaceTension(:)

      !> Hydrogen bond strength for each element
      real(wp), allocatable :: hBondStrength(:)

      !> Energy contribution from the SASA
      real(wp) :: gsasa

      !> Born radii
      real(wp), allocatable :: bornRad(:)

      !> Derivative of Born radii w.r.t. Cartesian coordinates
      real(wp), allocatable :: dbrdr(:, :, :)

      !> Derivative of Born radii w.r.t. strain deformations
      real(wp), allocatable :: dbrdL(:, :, :)

      !> Solvent accessible surface area
      real(wp), allocatable :: sasa(:)

      !> Derivative of SASA w.r.t. Cartesian coordinates
      real(wp), allocatable :: dsdr(:, :, :)

      !> Derivative of SASA w.r.t. strain deformations
      real(wp), allocatable :: dsdL(:, :, :)

      !> Born matrix
      real(wp), allocatable :: bornMat(:, :)

      !> Scratch memory for storing the shifts
      real(wp), allocatable :: shift(:)

   contains

      !> Update geometry dependent data
      procedure :: update

      !> Add potential shift
      procedure :: addShift

      !> Get energy contributions
      procedure :: getEnergy

      !> Add contributions to the gradient
      procedure :: getGradient

   end type TGeneralizedBorn


   !> Initialize GBSA model
   interface init
      module procedure :: initGeneralizedBorn
   end interface init


   !> P16 zeta parameter
   real(wp), parameter :: zetaP16 = 1.028_wp

   !> P16 zeta parameter over 16
   real(wp), parameter :: zetaP16o16 = zetaP16 / 16.0_wp

   !> exp zeta parameter
   real(wp), parameter :: zetaExp = 1.021_wp


contains


!> Initialize GBSA model from parametrisation data
subroutine initGeneralizedBorn(self, env, input, nAtom, rad, state, temperature, &
      & cutoff, kernel, num)

   !> Source for the error generation
   character(len=*), parameter :: source = 'solv_gbsa_initGeneralizedBorn'

   !> Instance of the solvation model
   type(TGeneralizedBorn), intent(out) :: self

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Parametrisation data
   type(TGBSAData), intent(in) :: input

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Van-der-Waals radii
   real(wp), intent(in) :: rad(:)

   !> Reference state for solution
   integer, intent(in), optional :: state

   !> Temperature of solution
   real(wp), intent(in), optional :: temperature

   !> Real space cutoff for the Born integrator
   real(wp), intent(in), optional :: cutoff

   !> Generalized Born kernel
   integer, intent(in), optional :: kernel

   !> Atomic numbers for each species
   integer, intent(in), optional :: num(:)

   integer :: ii
   real(wp) :: intCutoff

   self%nAtom = nAtom
   allocate(self%bornRad(nAtom))
   allocate(self%dbrdr(3, nAtom, nAtom))
   allocate(self%dbrdL(3, 3, nAtom))
   allocate(self%sasa(nAtom))
   allocate(self%dsdr(3, nAtom, nAtom))
   !allocate(self%dsdL(3, 3, nAtom))
   allocate(self%bornMat(nAtom, nAtom))
   allocate(self%shift(nAtom))

   if (present(kernel)) then
      self%gbKernel = kernel
   else
      self%gbKernel = gbKernel%still
   end if

   if (present(cutoff)) then
      intCutoff = cutoff
   else
      intCutoff = 35.0_wp*aatoau
   end if

   call init_(self%born, rad, input%descreening, input%bornScale, &
      & input%bornOffset, cutoff=intCutoff, num=num)

   call init_(self%surface, env, rad, input%probeRad, num=num)

   self%kEps = 1.0_wp / input%dielectricConst - 1.0_wp
   if (present(state) .and. present(temperature)) then
      self%freeEnergyShift = input%freeEnergyShift + getStateShift(state, &
         & temperature, input%density, input%molecularMass)
   else
      self%freeEnergyShift = input%freeEnergyShift
   end if

   if (present(num)) then
      allocate(self%surfaceTension(size(num)))
      do ii = 1, size(num)
         self%surfaceTension(ii) = input%surfaceTension(num(ii))
      end do
   else
      self%surfaceTension = input%surfaceTension
   end if

   if (allocated(input%hBondStrength)) then
      if (present(num)) then
         allocate(self%hBondStrength(size(num)))
         do ii = 1, size(num)
            self%hBondStrength(ii) = input%hBondStrength(num(ii))
         end do
      else
         self%hBondStrength = input%hBondStrength
      end if
   end if

end subroutine initGeneralizedBorn


!> Update Born radii and surface area from current geometry
subroutine update(self, neighList, num, bornRad)

   !> Instance of the solvation model
   class(TGeneralizedBorn), intent(inout) :: self

   !> Neighbourlist
   class(TNeighbourList), intent(in) :: neighList

   !> Atomic identifiers
   integer, intent(in) :: num(:)

   !> Born radii
   real(wp), intent(in), optional :: bornRad(:)

   !> Number of neighbours for each atom
   integer, allocatable :: neighs(:)

   integer :: iat, izp
   real(wp) :: hbterm

   allocate(neighs(self%nAtom))

   call neighlist%getNeighs(neighs, self%surface%cutoff)
   call self%surface%getSASA(neighs, neighList, num, self%sasa, self%dsdr)

   self%gsasa = 0.0_wp
   do iat = 1, self%nAtom
      izp = num(iat)
      self%gsasa = self%gsasa + self%surfaceTension(izp) * self%sasa(iat)
   end do

   if (present(bornRad)) then
      self%bornRad(:) = bornRad
      self%dbrdr(:, :, :) = 0.0_wp
      self%dbrdL(:, :, :) = 0.0_wp
   else
      call neighlist%getNeighs(neighs, self%born%cutoff)
      call self%born%getBornRad(neighs, neighList, num, self%bornRad, self%dbrdr, &
         & self%dbrdL)
   end if

   select case(self%gbKernel)
   case default
      call getBornMatrixStill(self%nAtom, neighList%coords, self%kEps, &
         & self%bornRad, self%bornMat)
   case(gbKernel%p16)
      call getBornMatrixP16(self%nAtom, neighList%coords, self%kEps, &
         & self%bornRad, self%bornMat)
   end select

   if (allocated(self%hBondStrength)) then
      do iat = 1, self%nAtom
         izp = num(iat)
         !hbTerm = self%sasa(iat) / (fourpi*self%surface%probeRad(izp)**2) &
         hbTerm = self%sasa(iat) / self%surface%probeRad(izp)**2 &
            & * self%hBondStrength(izp)
         self%bornMat(iat, iat) = self%bornMat(iat, iat) + hbTerm
      end do
   end if

end subroutine update


subroutine addShift(self, qat, qsh, atomicShift, shellShift)

   !> Instance of the solvation model
   class(TGeneralizedBorn), intent(inout) :: self

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Atomic potential shift
   real(wp), intent(inout) :: atomicShift(:)

   !> Shell-resolved potential shift
   real(wp), intent(inout) :: shellShift(:)

   call symv('l', self%nAtom, 1.0_wp, self%bornMat, self%nAtom, qat, 1, &
      & 1.0_wp, atomicShift, 1)

end subroutine addShift


!> Get energy from solvation model
pure subroutine getEnergy(self, qat, qsh, energy)

   !> Instance of the solvation model
   class(TGeneralizedBorn), intent(inout) :: self

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Third order contribution to the energy
   real(wp), intent(out) :: energy

   call symv('l', self%nAtom, 1.0_wp, self%bornMat, self%nAtom, qat, 1, &
      & 0.0_wp, self%shift, 1)
   energy = 0.5_wp * dot(self%nAtom, self%shift, 1, qat, 1) &
      & + self%gsasa + self%freeEnergyShift

end subroutine getEnergy


subroutine getGradient(self, neighList, num, qat, qsh, gradient, sigma)

   !> Instance of the solvation model
   class(TGeneralizedBorn), intent(inout) :: self

   !> Neighbourlist
   class(TNeighbourList), intent(in) :: neighList

   !> Atomic identifiers
   integer, intent(in) :: num(:)

   !> Partial charges
   real(wp), intent(in) :: qat(:)

   !> Partial charges
   real(wp), intent(in) :: qsh(:)

   !> Number of neighbours for each atom
   integer, allocatable :: neighs(:)

   !> Molecular Gradient
   real(wp), intent(inout) :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, izp
   real(wp) :: hbTerm
   real(wp), allocatable :: djdr(:, :, :)
   real(wp), allocatable :: djdtr(:, :)
   real(wp), allocatable :: djdL(:, :, :)

   allocate(djdr(3, self%nAtom, self%nAtom))
   allocate(djdtr(3, self%nAtom))
   allocate(djdL(3, 3, self%nAtom))

   select case(self%gbKernel)
   case default
      call getBornDerivStill(self%nAtom, neighList%coords, qat, self%kEps, &
         & self%bornRad, self%dbrdr, self%dbrdL, djdr, djdtr, djdL)
   case(gbKernel%p16)
      call getBornDerivP16(self%nAtom, neighList%coords, qat, self%kEps, &
         & self%bornRad, self%dbrdr, self%dbrdL, djdr, djdtr, djdL)
   end select

   call contract(djdr, qat, gradient, beta=1.0_wp)

   do iat = 1, self%nAtom
      izp = num(iat)
      gradient(:, :) = gradient + self%dsdr(:, :, iat) * self%surfaceTension(izp)
   end do

   if (allocated(self%hBondStrength)) then
      do iat = 1, self%nAtom
         izp = num(iat)
         hbTerm = qsh(iat)**2 / self%surface%probeRad(izp)**2 &
            & * self%hBondStrength(izp)
         gradient(:, :) = gradient + self%dsdr(:, :, iat) * hbTerm
      end do
   end if

end subroutine getGradient


!> Calculate Born matrix from Born radii. Still kernel is used for evaluation.
subroutine getBornMatrixStill(nAtom, xyz, kEps, bornRad, jmat)

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Dielectric screening
   real(wp), intent(in) :: kEps

   !> Born radii
   real(wp), intent(in) :: bornRad(:)

   !> Born matrix
   real(wp), intent(out) :: jmat(:, :)

   integer :: iat, jat
   real(wp) :: r2, ab, arg, eab, fgb2, dfgb

   jmat(:, :) = 0.0_wp

   !$omp parallel do default(none) shared(jmat, nAtom, xyz, bornRad, kEps) &
   !$omp private(iat, jat, r2, ab, arg, eab, fgb2, dfgb)
   do iat = 1, nAtom
      do jat = 1, iat-1
         r2 = sum((xyz(:, iat) - xyz(:, jat))**2)
         ab = bornRad(iat) * bornRad(jat)
         arg = 0.25_wp * r2 / ab
         eab = exp(-arg)
         fgb2 = r2 + ab*eab
         dfgb = 1.0_wp / sqrt(fgb2)
         jmat(jat, iat) = dfgb * kEps
         jmat(iat, jat) = dfgb * kEps
      end do
      jmat(iat, iat) = kEps / bornRad(iat)
   end do
   !$omp end parallel do

end subroutine getBornMatrixStill


!> Calculate Born matrix from Born radii. Still kernel is used for evaluation.
subroutine getBornMatrixP16(nAtom, xyz, kEps, bornRad, jmat)

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Dielectric screening
   real(wp), intent(in) :: kEps

   !> Born radii
   real(wp), intent(in) :: bornRad(:)

   !> Born matrix
   real(wp), intent(out) :: jmat(:, :)

   integer :: iat, jat
   real(wp) :: r1, ab, arg, eab, fgb, dfgb

   jmat(:, :) = 0.0_wp

   !$omp parallel do default(none) shared(jmat, nAtom, xyz, bornRad, kEps) &
   !$omp private(iat, jat, r1, ab, arg, fgb, dfgb)
   do iat = 1, nAtom
      do jat = 1, iat-1
         r1 = sqrt(sum((xyz(:, iat) - xyz(:, jat))**2))
         ab = sqrt(bornRad(iat) * bornRad(jat))
         arg = ab / (ab + zetaP16o16*r1) ! ab / (1 + ζR/(16·ab))
         arg = arg * arg ! ab / (1 + ζR/(16·ab))²
         arg = arg * arg ! ab / (1 + ζR/(16·ab))⁴
         arg = arg * arg ! ab / (1 + ζR/(16·ab))⁸
         arg = arg * arg ! ab / (1 + ζR/(16·ab))¹⁶
         fgb = r1 + ab*arg
         dfgb = 1.0_wp / fgb
         jmat(jat, iat) = dfgb * kEps
         jmat(iat, jat) = dfgb * kEps
      end do
      jmat(iat, iat) = kEps / bornRad(iat)
   end do
   !$omp end parallel do

end subroutine getBornMatrixP16


!> Calculate unfolded derivatives of the Born matrix.
subroutine getBornDerivStill(nAtom, xyz, qVec, kEps, bornRad, dbrdr, dbrdL, &
      & djdr, djdtr, djdL)

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Partial charges
   real(wp), intent(in) :: qVec(:)

   !> Dielectric screening
   real(wp), intent(in) :: kEps

   !> Born radii
   real(wp), intent(in) :: bornRad(:)

   !> Derivative of Born radii w.r.t. Cartesian coordinates
   real(wp), intent(in) :: dbrdr(:, :, :)

   !> Derivative of Born radii w.r.t. strain deformations
   real(wp), intent(in) :: dbrdL(:, :, :)

   !> Derivative of Born matrix w.r.t. Cartesian coordinates
   real(wp), intent(out) :: djdr(:, :, :)

   !> Trace derivative of Born matrix w.r.t. Cartesian coordinates
   real(wp), intent(out) :: djdtr(:, :)

   !> Derivative of Born radii w.r.t. strain deformations
   real(wp), intent(out) :: djdL(:, :, :)

   integer :: iat, jat
   real(wp) :: vec(3), r2, ab, arg, eab, fgb2, dfgb, dfgb2, dfgb3, dG(3), ap, bp
   real(wp) :: dEdbri, dEdbrj

   djdr(:, :, :) = 0.0_wp
   djdtr(:, :) = 0.0_wp
   djdL(:, :, :) = 0.0_wp

   !$omp parallel do default(none) reduction(+:djdr, djdtr, djdL) &
   !$omp shared(nAtom, xyz, bornRad, qvec, kEps, dbrdr, dbrdL) &
   !$omp private(iat, jat, vec, r2, ab, arg, eab, fgb2, dfgb, dfgb2, dfgb3, ap, &
   !$omp& bp, dEdbri, dEdbrj, dG)
   do iat = 1, nAtom
      do jat = 1, iat-1
         vec(:) = xyz(:, iat) - xyz(:, jat)
         r2 = sum(vec**2)

         ab = bornRad(iat) * bornRad(jat)
         arg = 0.25_wp * r2 / ab
         eab = exp(-arg)
         fgb2 = r2 + ab*eab
         dfgb2 = 1.0_wp / fgb2
         dfgb = sqrt(dfgb2)
         dfgb3 = dfgb*dfgb2

         ap = (1.0_wp-0.25_wp*eab)*dfgb3
         dG(:) = ap * vec * kEps
         djdr(:, iat, jat) = djdr(:, iat, jat) - dG*qvec(iat)
         djdr(:, jat, iat) = djdr(:, jat, iat) + dG*qvec(jat)

         bp = -0.5_wp*eab*(1.0_wp+arg)*dfgb3
         dEdbri = bornRad(jat) * bp * kEps
         dEdbrj = bornRad(iat) * bp * kEps
         djdr(:, :, jat) = djdr(:, :, jat) + dEdbrj * dbrdr(:, :, jat) * qvec(iat)
         djdr(:, :, iat) = djdr(:, :, iat) + dEdbri * dbrdr(:, :, iat) * qvec(jat)

      end do

      bp = 1.0_wp/bornRad(iat)
      dEdbri = -kEps*bp*bp*0.5_wp
      djdr(:, :, iat) = djdr(:, :, iat) + dEdbri * dbrdr(:, :, iat) * qvec(iat)
   end do
   !$omp end parallel do

end subroutine getBornDerivStill


!> Calculate unfolded derivatives of the Born matrix.
subroutine getBornDerivP16(nAtom, xyz, qVec, kEps, bornRad, dbrdr, dbrdL, &
      & djdr, djdtr, djdL)

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Partial charges
   real(wp), intent(in) :: qVec(:)

   !> Dielectric screening
   real(wp), intent(in) :: kEps

   !> Born radii
   real(wp), intent(in) :: bornRad(:)

   !> Derivative of Born radii w.r.t. Cartesian coordinates
   real(wp), intent(in) :: dbrdr(:, :, :)

   !> Derivative of Born radii w.r.t. strain deformations
   real(wp), intent(in) :: dbrdL(:, :, :)

   !> Derivative of Born matrix w.r.t. Cartesian coordinates
   real(wp), intent(out) :: djdr(:, :, :)

   !> Trace derivative of Born matrix w.r.t. Cartesian coordinates
   real(wp), intent(out) :: djdtr(:, :)

   !> Derivative of Born radii w.r.t. strain deformations
   real(wp), intent(out) :: djdL(:, :, :)

   integer :: iat, jat
   real(wp) :: vec(3), r2, r1, ab, arg1, arg16, fgb, fgb2, dfgb, dfgb2
   real(wp) :: dEdbri, dEdbrj, dG(3), ap, bp

   djdr(:, :, :) = 0.0_wp
   djdtr(:, :) = 0.0_wp
   djdL(:, :, :) = 0.0_wp

   !$omp parallel do default(none) reduction(+:djdr, djdtr, djdL) &
   !$omp shared(nAtom, xyz, bornRad, qvec, kEps, dbrdr, dbrdL) &
   !$omp private(iat, jat, vec, r1, r2, ab, arg1, arg16, fgb, dfgb, dfgb2, ap, &
   !$omp& bp, dEdbri, dEdbrj, dG)
   do iat = 1, nAtom
      do jat = 1, iat-1
         vec(:) = xyz(:, iat) - xyz(:, jat)
         r2 = sum(vec**2)
         r1 = sqrt(r2)

         ab = sqrt(bornRad(iat) * bornRad(jat))
         arg1 = ab / (ab + zetaP16o16*r1) ! 1 / (1 + ζR/(16·ab))
         arg16 = arg1 * arg1 ! 1 / (1 + ζR/(16·ab))²
         arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁴
         arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁸
         arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))¹⁶

         fgb = r1 + ab*arg16
         dfgb = 1.0_wp / fgb
         dfgb2 = dfgb * dfgb

         ! (1 - ζ/(1 + Rζ/(16 ab))^17)/(R + ab/(1 + Rζ/(16 ab))¹⁶)²
         ap = (1.0_wp - zetaP16 * arg1 * arg16) * dfgb2
         dG(:) = ap * vec * kEps / r1
         djdr(:, iat, jat) = djdr(:, iat, jat) - dG*qvec(iat)
         djdr(:, jat, iat) = djdr(:, jat, iat) + dG*qvec(jat)

         ! -(Rζ/(2·ab²·(1 + Rζ/(16·ab))¹⁷) + 1/(2·ab·(1 + Rζ/(16·ab))¹⁶))/(R + ab/(1 + Rζ/(16·ab))¹⁶)²
         bp = -0.5_wp*(r1 * zetaP16 / ab * arg1 + 1.0_wp) / ab * arg16 * dfgb2
         dEdbri = bornRad(jat) * bp * kEps
         dEdbrj = bornRad(iat) * bp * kEps
         djdr(:, :, jat) = djdr(:, :, jat) + dEdbrj * dbrdr(:, :, jat) * qvec(iat)
         djdr(:, :, iat) = djdr(:, :, iat) + dEdbri * dbrdr(:, :, iat) * qvec(jat)

      end do

      bp = 1.0_wp/bornRad(iat)
      dEdbri = -kEps*bp*bp*0.5_wp
      djdr(:, :, iat) = djdr(:, :, iat) + dEdbri * dbrdr(:, :, iat) * qvec(iat)
   end do
   !$omp end parallel do

end subroutine getBornDerivP16


end module xtb_solv_gbsa
