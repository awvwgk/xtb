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

!> TODO
module xtb_solv_sasa
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_solv_lebedev, only : getAngGrid, gridSize
   use xtb_type_environment, only : TEnvironment
   use xtb_type_neighbourlist, only : TNeighbourList
   implicit none
   private


   public :: TSurfaceIntegrator, init


   !> Default angular Lebedev-Laikov grid, corresponds to 230 grid points
   integer, parameter :: angDefault = 12

   !> Default smoothing parameter for numerical integration
   real(wp), parameter :: smoothingDefault = 0.3_wp*aatoau

   !> Default offset for generation of neighbourlist
   real(wp), parameter :: offsetDefault = 2.0_wp*aatoau

   !> Default tolerance for including contributions to surface area
   real(wp), parameter :: tolDefault = 1.0e-6_wp


   !> Data for the solvent accessible surface area model
   type :: TSurfaceIntegrator

      !> Real space cutoff
      real(wp) :: cutoff

      !> Real space cut-offs for surface area contribution
      real(wp) :: tolerance

      !> Smoothing dielectric function parameters
      real(wp) :: smoothingPar(3)

      !> Van-der-Waals radii + probe radius
      real(wp), allocatable :: probeRad(:)

      !> Surface tension for each species
      real(wp), allocatable :: surfaceTension(:)

      !> Thresholds for smooth numerical integration
      real(wp), allocatable :: thresholds(:, :)

      !> Radial weight
      real(wp), allocatable :: radWeight(:)

      !> Angular grid for surface integration
      real(wp), allocatable :: angGrid(:, :)

      !> Weights of grid points for surface integration
      real(wp), allocatable :: angWeight(:)

   contains

      !> Calculate solvent accessible surface area for every atom
      procedure :: getSASA

   end type TSurfaceIntegrator


   !> Initialize solvent accessible surface area model from input data
   interface init
      module procedure :: initSurfaceIntegrator
   end interface init


contains


!> Initialize solvent accessible surface area model from input data
subroutine initSurfaceIntegrator(self, env, vdwRad, probeRad, offset, angSize, &
      & smoothingPar, tolerance, num)

   !> Source for error creation
   character(len=*), parameter :: source = 'solv_sasa_initSurfaceIntegrator'

   !> Data structure
   type(TSurfaceIntegrator), intent(out) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Real space cutoff
   real(wp), intent(in), optional :: offset

   !> Grid for numerical integration of atomic surfaces
   integer, intent(in), optional :: angSize

   !> Probe radius for solvent molecules in atomic surfaces calculations
   real(wp), intent(in) :: probeRad

   !> Smoothing dielectric function parameters
   real(wp), intent(in), optional :: smoothingPar

   !> Real space cut-offs for surface area contribution
   real(wp), intent(in), optional :: tolerance

   !> Van-der-Waals radii
   real(wp), intent(in) :: vdwRad(:)

   !> Mapping from identity to atomic numbers
   integer, intent(in), optional :: num(:)

   integer :: ii, nId, stat, angGrid
   real(wp) :: wPar

   if (present(angSize)) then
      angGrid = angSize
   else
      angGrid = angDefault
   end if
   allocate(self%angGrid(3, gridSize(angGrid)))
   allocate(self%angWeight(gridSize(angGrid)))
   call getAngGrid(angGrid, self%angGrid, self%angWeight, stat)
   if (stat /= 0) then
      call env%error("Cannot generate angular grid from invalid size", source)
      return
   end if

   if (present(num)) then
      allocate(self%probeRad(size(num)))
      do ii = 1, size(num)
         self%probeRad(ii) = vdwRad(num(ii)) + probeRad
      end do
   else
      allocate(self%probeRad(size(vdwRad)))
      self%probeRad(:) = probeRad + vdwRad
   end if

   if (present(smoothingPar)) then
      wPar = smoothingPar
   else
      wPar = smoothingDefault
   end if
   nId = size(self%probeRad)
   allocate(self%thresholds(2, nId))
   allocate(self%radWeight(nId))
   call getIntegrationParam(nId, wPar, self%smoothingPar, &
      & self%probeRad, self%thresholds, self%radWeight)

   self%cutoff = 2.0_wp * (maxval(self%probeRad) + wPar)
   if (present(offset)) then
      self%cutoff = self%cutoff + offset
   else
      self%cutoff = self%cutoff + offsetDefault
   end if
   if (present(tolerance)) then
      self%tolerance = tolerance
   else
      self%tolerance = tolDefault
   end if

end subroutine initSurfaceIntegrator


!> Get parameters for smooth numerical integration
subroutine getIntegrationParam(nSpecies, w, ah, rad, thr, wrp)

   !> Number of species
   integer, intent(in) :: nSpecies

   !> Smoothing parameter
   real(wp), intent(in) :: w

   !> Smoothing parameters
   real(wp), intent(out) :: ah(3)

   !> Probe radii
   real(wp), intent(in) :: rad(:)

   !> Numerical thresholds
   real(wp), intent(out) :: thr(:, :)

   !> Radial weight
   real(wp), intent(out) :: wrp(:)

   integer :: iSp1
   real(wp) :: w3, rm, rp

   ah(1) = 0.5_wp
   ah(2) = 3.0_wp/(4.0_wp*w)
   ah(3) = -1.0_wp/(4.0_wp*w*w*w)

   do iSp1 = 1, nSpecies
      rm = rad(iSp1) - w
      rp = rad(iSp1) + w
      thr(1, iSp1) = rm**2
      thr(2, iSp1) = rp**2
      wrp(iSp1) = (0.25_wp/w + 3.0_wp*ah(3)*(0.2_wp*rp*rp-0.5_wp*rp*rad(iSp1) &
         & + rad(iSp1)*rad(iSp1)/3.0_wp))*rp*rp*rp - (0.25_wp/w &
         & + 3.0_wp*ah(3)*(0.2_wp*rm*rm - 0.5_wp*rm*rad(iSp1) &
         & + rad(iSp1)*rad(iSp1)/3.0_wp))*rm*rm*rm
   end do

end subroutine getIntegrationParam


!> Calculate solvent accessible surface area for every atom
pure subroutine getSASA(self, neighs, neighList, num, sasa, dsdr)

   !> Data structure
   class(TSurfaceIntegrator), intent(inout) :: self

   !> Nr. of neighbours for each atom
   integer, intent(in) :: neighs(:)

   !> Neighbourlist
   type(TNeighbourList), intent(in) :: neighList

   !> Species, shape: [nAtom]
   integer, intent(in) :: num(:)

   !> Solvent accessible surface area
   real(wp), intent(out) :: sasa(:)

   !> Derivative of solvent accessible surface area w.r.t. coordinates
   real(wp), intent(out) :: dsdr(:, :, :)

   integer iat, izp, img, jat, jzp, ij, ip
   integer :: nAtom, mNeighbour, nNeighs, nEval
   real(wp) :: vec(3), dist2, dist, weight
   real(wp) :: uj, uj3, ah3uj2
   real(wp) :: sasaij, dsasaij
   real(wp) :: rsas, sasai, point(3), sasap, wsa, dGr(3)
   real(wp), allocatable :: grds(:,:), derivs(:,:)
   integer, allocatable :: grdi(:)

   nAtom = size(neighs)
   mNeighbour = maxval(neighs)

   allocate(derivs(3,nAtom))
   allocate(grds(3,mNeighbour))
   allocate(grdi(mNeighbour))

   do iat = 1, nAtom
      izp = num(iat)

      rsas = self%probeRad(izp)

      ! initialize storage
      derivs(:, :) = 0.0_wp
      sasai = 0.0_wp

      ! loop over grid points
      do ip = 1, size(self%angGrid, dim=2)
         ! grid point position
         point(:) = neighList%coords(:,iat) + rsas * self%angGrid(:, ip)

         ! atomic surface function at the grid point
         nEval = 0
         grds(:, :) = 0.0_wp
         grdi(:) = 0
         sasap = 1.0_wp
         do ij = 1, neighs(iat)
            img = neighList%ineigh(ij, iat)
            jat = neighList%image(img)
            jzp = num(jat)
            ! compute the distance to the atom
            vec(:) = point(:) - neighList%coords(:, img)
            dist2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            ! if within the outer cut-off compute
            if (dist2 < self%thresholds(2,jzp)) then
               if (dist2 < self%thresholds(1,jzp)) then
                  sasap = 0.0_wp
                  exit
               else
                  dist = sqrt(dist2)
                  uj = dist - self%probeRad(jzp)
                  ah3uj2 = self%smoothingPar(3)*uj*uj
                  dsasaij = self%smoothingPar(2)+3.0_wp*ah3uj2
                  sasaij =  self%smoothingPar(1)+(self%smoothingPar(2)+ah3uj2)*uj

                  ! accumulate the molecular surface
                  sasap = sasap*sasaij
                  ! compute the gradient wrt the neighbor
                  dsasaij = dsasaij/(sasaij*dist)
                  nEval = nEval+1
                  grdi(nEval) = jat
                  grds(:, nEval) = dsasaij*vec(:)
               end if
            end if
         end do

         if (sasap > self%tolerance) then
            ! numerical quadrature weight
            wsa = self%angWeight(ip)*self%radWeight(izp)*sasap
            ! accumulate the surface area
            sasai = sasai + wsa
            ! accumulate the surface gradient
            do ij = 1, nEval
               jat = grdi(ij)
               dGr(:) = wsa * grds(:, ij)
               derivs(:, iat) = derivs(:, iat) + dGr(:)
               derivs(:, jat) = derivs(:, jat) - dGr(:)
            end do
         end if
      end do

      sasa(iat) = sasai
      dsdr(:,:,iat) = derivs

   end do

end subroutine getSASA


end module xtb_solv_sasa
