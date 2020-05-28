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
   use xtb_mctc_blas, only : mctc_symv, mctc_dot, mctc_gemv
   use xtb_mctc_convert, only : aatoau
   use xtb_mctc_math, only : matDet3x3, derivDet3x3
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

      !> Analytical linearized Possion Boltzmann factor
      real(wp) :: alpbet

      !> Dielectric scaling
      real(wp) :: kEps

      !> Shift for free energy calculations
      real(wp) :: freeEnergyShift

      !> Surface tension scaling for each element
      real(wp), allocatable :: surfaceTension(:)

      !> Hydrogen bond strength for each element
      real(wp), allocatable :: hBondStrength(:)

      !> Van-der-Waals radii
      real(wp), allocatable :: vdwRad(:)

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

      !> Get realspace cutoff for neighbour lists
      procedure :: getCutoff

   end type TGeneralizedBorn


   !> Initialize GBSA model
   interface init
      module procedure :: initGeneralizedBorn
   end interface init


   !> Parameter for the analytical linearized Possion Boltzmann
   real(wp), parameter :: alpha = 0.571412_wp

   !> P16 zeta parameter
   real(wp), parameter :: zetaP16 = 1.028_wp

   !> P16 zeta parameter over 16
   real(wp), parameter :: zetaP16o16 = zetaP16 / 16.0_wp

   !> exp zeta parameter
   real(wp), parameter :: zetaExp = 1.021_wp


contains


!> Initialize GBSA model from parametrisation data
subroutine initGeneralizedBorn(self, env, input, nAtom, rad, state, temperature, &
      & cutoff, kernel, alpb, num)

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

   !> Construct an analytical linearized Possion-Boltzmann model instead
   logical, intent(in), optional :: alpb

   !> Atomic numbers for each species
   integer, intent(in), optional :: num(:)

   integer :: ii
   logical :: isALPB
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

   if (present(alpb)) then
      isALPB = alpb
   else
      isALPB = .false.
   end if

   call init_(self%born, rad, input%descreening, input%bornScale, &
      & input%bornOffset, cutoff=intCutoff, num=num)

   call init_(self%surface, env, rad, input%probeRad, num=num)

   if (isALPB) then
      if (present(num)) then
         allocate(self%vdwRad(size(num)))
         do ii = 1, size(num)
            self%vdwRad(ii) = rad(num(ii))
         end do
      else
         self%vdwRad = rad
      end if
      ! β = 1 / ε, we save αβ because they are always used together
      self%alpbet = alpha / input%dielectricConst
   else
      self%alpbet = 0.0_wp
   end if
   self%kEps = (1.0_wp / input%dielectricConst - 1.0_wp) / (1.0_wp + self%alpbet)
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
subroutine update(self, neighList, num)

   !> Instance of the solvation model
   class(TGeneralizedBorn), intent(inout) :: self

   !> Neighbourlist
   class(TNeighbourList), intent(in) :: neighList

   !> Atomic identifiers
   integer, intent(in) :: num(:)

   !> Number of neighbours for each atom
   integer, allocatable :: neighs(:)

   integer :: iat, izp
   real(wp) :: hbterm
   real(wp) :: aDet

   allocate(neighs(self%nAtom))

   call neighlist%getNeighs(neighs, self%surface%cutoff)
   call self%surface%getSASA(neighs, neighList, num, self%sasa, self%dsdr)

   self%gsasa = 0.0_wp
   do iat = 1, self%nAtom
      izp = num(iat)
      self%gsasa = self%gsasa + self%surfaceTension(izp) * self%sasa(iat)
   end do

   call neighlist%getNeighs(neighs, self%born%cutoff)
   call self%born%getBornRad(neighs, neighList, num, self%bornRad, self%dbrdr, &
      & self%dbrdL)

   select case(self%gbKernel)
   case default
      call getBornMatrixStill(self%nAtom, neighList%coords, self%kEps, &
         & self%bornRad, self%bornMat)
   case(gbKernel%p16)
      call getBornMatrixP16(self%nAtom, neighList%coords, self%kEps, &
         & self%bornRad, self%bornMat)
   end select

   if (self%alpbet > 0.0_wp) then
      call getADet(self%nAtom, neighList%coords, num, self%vdwRad, aDet)
      self%bornMat = self%bornMat + self%alpbet / aDet
   end if

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


!> Return realspace cutoff for the generation of neighbour lists
pure function getCutoff(self) result(cutoff)

   !> Instance of the solvation model
   class(TGeneralizedBorn), intent(in) :: self

   !> Maximal needed real space cutoff
   real(wp) :: cutoff

   cutoff = max(self%surface%cutoff, self%born%cutoff)

end function getCutoff


!> Add potential shift from solvation model
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

   call mctc_symv(self%bornMat, qat, atomicShift, beta=1.0_wp)

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

   call mctc_symv(self%bornMat, qat, self%shift)
   energy = 0.5_wp * mctc_dot(self%shift, qat) &
      & + self%gsasa + self%freeEnergyShift

end subroutine getEnergy


!> Get gradient contributions from solvation model
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

   call mctc_gemv(djdr, qat, gradient, beta=1.0_wp)

   do iat = 1, self%nAtom
      izp = num(iat)
      gradient(:, :) = gradient + self%dsdr(:, :, iat) * self%surfaceTension(izp)
   end do

   if (self%alpbet > 0.0_wp) then
      call getADetDeriv(self%nAtom, neighList%coords, num, self%vdwRad, &
         & self%kEps*self%alpbet, qat, gradient, sigma)
   end if

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


subroutine getADet(nAtom, xyz, num, rad, aDet)

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Species identifiers for each atom
   integer, intent(in) :: num(:)

   !> Atomic radii
   real(wp), intent(in) :: rad(:)

   !> Shape descriptor of the structure
   real(wp), intent(out) :: aDet

   integer :: iat, izp
   real(wp) :: r2, rad2, rad3, totRad3, vec(3), center(3), inertia(3, 3)
   real(wp), parameter :: tof = 2.0_wp/5.0_wp, unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   totRad3 = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nAtom
      izp = num(iat)
      rad2 = rad(izp) * rad(izp)
      rad3 = rad2 * rad(izp)
      totRad3 = totRad3 + rad3
      center(:) = center + xyz(:, iat) * rad3
   end do
   center = center / totRad3

   inertia(:, :) = 0.0_wp
   do iat = 1, nAtom
      izp = num(iat)
      rad2 = rad(izp) * rad(izp)
      rad3 = rad2 * rad(izp)
      vec(:) = xyz(:, iat) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
         & - spread(vec, 1, 3) * spread(vec, 2, 3))
   end do

   aDet = sqrt(matDet3x3(inertia)**(1.0_wp/3.0_wp)/(tof*totRad3))

end subroutine getADet


subroutine getADetDeriv(nAtom, xyz, num, rad, kEps, qvec, gradient, sigma)

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Species identifiers for each atom
   integer, intent(in) :: num(:)

   !> Atomic radii
   real(wp), intent(in) :: rad(:)

   real(wp), intent(in) :: kEps
   real(wp), intent(in) :: qvec(:)

   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(inout) :: sigma(:, :)

   integer :: iat, izp
   real(wp) :: r2, rad2, rad3, totRad3, vec(3), center(3), inertia(3, 3), aDet
   real(wp) :: aDeriv(3, 3), qtotal
   real(wp), parameter :: tof = 2.0_wp/5.0_wp, unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   qtotal = 0.0_wp
   totRad3 = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nAtom
      izp = num(iat)
      rad2 = rad(izp) * rad(izp)
      rad3 = rad2 * rad(izp)
      totRad3 = totRad3 + rad3
      center(:) = center + xyz(:, iat) * rad3
      qtotal = qtotal + qvec(iat)
   end do
   center = center / totRad3

   inertia(:, :) = 0.0_wp
   do iat = 1, nAtom
      izp = num(iat)
      rad2 = rad(izp) * rad(izp)
      rad3 = rad2 * rad(izp)
      vec(:) = xyz(:, iat) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
         & - spread(vec, 1, 3) * spread(vec, 2, 3))
   end do

   aDet = sqrt(matDet3x3(inertia)**(1.0_wp/3.0_wp)/(tof*totRad3))
   aDeriv(:, :) = reshape([&
      & inertia(1,1)*(inertia(2,2)+inertia(3,3))-inertia(1,2)**2-inertia(1,3)**2, &
      & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), &
      & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), &
      & inertia(1,2)*inertia(3,3)-inertia(1,3)-inertia(2,3), &
      & inertia(2,2)*(inertia(1,1)+inertia(3,3))-inertia(1,2)**2-inertia(2,3)**2, &
      & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), &
      & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), &
      & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), &
      & inertia(3,3)*(inertia(1,1)+inertia(2,2))-inertia(1,3)**2-inertia(2,3)**2],&
      & shape=[3, 3], order=[2, 1]) &
      & * 125.0_wp / (48.0_wp * aDet**5 * totRad3**3) &
      & * (-0.5_wp * kEps * qtotal**2 / aDet**2)

   do iat = 1, nAtom
      izp = num(iat)
      rad2 = rad(izp) * rad(izp)
      rad3 = rad2 * rad(izp)
      vec(:) = xyz(:, iat) - center
      gradient(:, iat) = gradient(:, iat) + rad3 * matmul(aDeriv, vec)
   end do

end subroutine getADetDeriv


end module xtb_solv_gbsa
