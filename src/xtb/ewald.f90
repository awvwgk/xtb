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

!> Implementation of Ewald summation specific tasks
module xtb_xtb_ewald
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   implicit none
   private

   public :: ewaldMatPBC3D, ewaldDerivPBC3D

contains


pure function ewaldMatPBC3D(vec, recPoint, qpc, volume, alpha, scale) result(Amat)

   !> Distance from i to WSC atom
   real(wp), intent(in) :: vec(:)

   !> Reciprocal lattice
   real(wp), intent(in) :: recPoint(:, :)

   !> Pseudo-quadrupole charge
   real(wp), intent(in) :: qpc

   !> Direct cell volume
   real(wp), intent(in) :: volume

   !> Convergence factor
   real(wp), intent(in) :: alpha

   !> Additional scaling factor
   real(wp), intent(in) :: scale

   !> Element of interaction matrix
   real(wp) :: Amat

   integer :: iG
   real(wp) :: rik2, rik(3), expterm

   Amat = 0.0_wp
   do iG = 1, size(recPoint, dim=2)
      rik(:) = recPoint(:, iG)
      rik2 = dot_product(rik, rik)
      expterm = exp(-rik2/(4.0_wp*alpha**2))/rik2
      Amat = Amat + cos(dot_product(rik, vec)) * expterm * (1.0_wp + 2*rik2*qpc**2)
   end do
   Amat = Amat * 4.0_wp*pi/volume * scale

end function ewaldMatPBC3D


pure subroutine ewaldDerivPBC3D(vec, recPoint, qpc, volume, alpha, scale, &
      & dAmat, sigma)

   !> Distance from i to WSC atom
   real(wp),intent(in) :: vec(:)

   !> Reciprocal lattice
   real(wp),intent(in) :: recPoint(:, :)

   !> Pseudo-quadrupole charge
   real(wp), intent(in) :: qpc

   !> Direct cell volume
   real(wp),intent(in) :: volume

   !> Convergence factor
   real(wp),intent(in) :: alpha

   !> Additional scaling factor
   real(wp), intent(in) :: scale

   !> Derivative of interaction matrix
   real(wp),intent(out) :: dAmat(:)

   !> Strain of interaction matrix
   real(wp),intent(out) :: sigma(:, :)

   integer  :: iG
   real(wp) :: rik2, rik(3), dtmp, expterm, arg
   real(wp) :: fqpc, falp, dS(3, 3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   dAmat = 0.0_wp
   sigma = 0.0_wp
   fqpc = 2*qpc**2
   falp = 2.0_wp/3.0_wp/(4.0_wp*alpha**2)
   do iG = 1, size(recPoint, dim=2)
      rik = recPoint(:, iG)
      rik2 = dot_product(rik,rik)
      expterm = exp(-rik2/(4.0_wp*alpha**2))/rik2
      arg = dot_product(rik,vec)
      dtmp = -sin(arg) * expterm
      dAmat = dAmat + rik*dtmp
      dS = spread(rik,1,3)*spread(rik,2,3)
      sigma = sigma + expterm * cos(arg) * ( &
         & - unity * (1.0_wp + rik2*falp + rik2*fqpc) &
         & + (2.0_wp/rik2 + 0.5_wp/alpha**2 + 0.5_wp*fqpc) * dS)
   end do
   dAmat = dAmat * 4.0_wp*pi/volume * scale
   sigma = sigma * 4.0_wp*pi/volume * scale

end subroutine ewaldDerivPBC3D


end module xtb_xtb_ewald
