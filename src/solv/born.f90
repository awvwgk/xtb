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

!> Implementation for Born radii evaluation
module xtb_solv_born
   use xtb_mctc_accuracy, only : wp
   use xtb_type_neighbourlist, only : TNeighbourList
   implicit none
   private

   public :: TBornIntegrator, init


   !> Default parameters for GB^OBC II
   real(wp), parameter :: obcDefault(3) = [1.00_wp, 0.80_wp, 4.85_wp]


   !> Integrator for Born radii
   type :: TBornIntegrator

      !> Real space cutoff
      real(wp) :: cutoff

      !> Scaling factor for Born radii
      real(wp) :: bornScale

      !> Offset parameter for Born radii
      real(wp) :: bornOffset

      !> Onufriev--Bashford--Case correction to Born radii
      real(wp) :: obcParam(3)

      !> Van-der-Waals radii
      real(wp), allocatable :: rad(:)

      !> Pair descreening approximation radii
      real(wp), allocatable :: rho(:)

   contains

      !> Calculate Born radii for a given geometry
      procedure :: getBornRad

   end type TBornIntegrator


   !> Initialize integrator
   interface init
      module procedure :: initBornIntegrator
   end interface init


contains


!> Initialize integrator for Born radii
subroutine initBornIntegrator(self, rad, descreening, bornScale, bornOffset, &
      & obcParam, cutoff, num)

   !> Data structure
   type(TBornIntegrator), intent(out) :: self

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Scaling factor for Born radii
   real(wp), intent(in) :: bornScale

   !> Offset parameter for Born radii
   real(wp), intent(in) :: bornOffset

   !> Onufriev--Bashford--Case correction to Born radii
   real(wp), intent(in), optional :: obcParam(3)

   !> Van-der-Waals radii
   real(wp), intent(in) :: rad(:)

   !> Pair descreening approximation radii
   real(wp), intent(in) :: descreening(:)

   !> Mapping from identity to atomic numbers
   integer, intent(in), optional :: num(:)

   integer :: ii

   self%cutoff = cutoff
   self%bornScale = bornScale
   self%bornOffset = bornOffset
   if (present(obcParam)) then
      self%obcParam = obcParam
   else
      self%obcParam = obcDefault
   end if

   if (present(num)) then
      allocate(self%rad(size(num)))
      do ii = 1, size(num)
         self%rad(ii) = rad(num(ii))
      end do
   else
      self%rad = rad
   end if

   if (present(num)) then
      allocate(self%rho(size(num)))
      do ii = 1, size(num)
         self%rho(ii) = rad(num(ii)) * descreening(num(ii))
      end do
   else
      self%rho = rad * descreening
   end if

end subroutine initBornIntegrator


!> Calculate Born radii for a given geometry
pure subroutine getBornRad(self, neighs, neighList, num, bornRad, dbrdr, dbrdL)

   !> Data structure
   class(TBornIntegrator), intent(in) :: self

   !> Nr. of neighbours for each atom
   integer, intent(in) :: neighs(:)

   !> Neighbourlist
   class(TNeighbourList), intent(in) :: neighList

   !> Central cell chemical species
   integer, intent(in) :: num(:)

   !> Born radii
   real(wp), intent(out) :: bornRad(:)

   !> Gradient of the Born radii
   real(wp), intent(out) :: dbrdr(:, :, :)

   !> Strain derivative of the Born radii
   real(wp), intent(out) :: dbrdL(:, :, :)

   integer :: iat, izp
   real(wp) :: br, dpsi, svdwi,vdwri
   real(wp) :: s1, v1, s2, arg, arg2, th, ch
   real(wp), allocatable :: psi(:), dpsidr, dpsidL

   call getPsi(neighs, neighList, num, self%rad, self%rho, &
      & bornRad, dbrdr, dbrdL)

   do iat = 1, size(bornRad)
      izp = num(iat)

      br = bornRad(iat)

      svdwi = self%rad(izp) - self%bornOffset
      vdwri = self%rad(izp)
      s1 = 1.0_wp/svdwi
      v1 = 1.0_wp/vdwri
      s2 = 0.5_wp*svdwi

      br = br*s2

      arg2 = br*(self%obcParam(3)*br-self%obcParam(2))
      arg = br*(self%obcParam(1)+arg2)
      arg2 = 2.0_wp*arg2+self%obcParam(1)+self%obcParam(3)*br*br

      th = tanh(arg)
      ch = cosh(arg)

      br = 1.0_wp/(s1-v1*th)
      ! Include GBMV2-like scaling
      br = self%bornScale*br

      dpsi = ch*(s1-v1*th)
      dpsi = s2*v1*arg2/(dpsi*dpsi)
      dpsi = self%bornScale*dpsi

      bornRad(iat) = br
      dbrdr(:, :, iat) = dbrdr(:, :, iat) * dpsi
      dbrdL(:, :, iat) = dbrdL(:, :, iat) * dpsi

   end do

end subroutine getBornRad


!> Evaluate volume integrals, intermediate values are stored in Born radii fields
pure subroutine getPsi(neighs, neighList, num, rad, rho, psi, dpsidr, dpsidL)

   !> Nr. of neighbours for each atom
   integer, intent(in) :: neighs(:)

   !> Neighbourlist
   class(TNeighbourList), intent(in) :: neighList

   !> Central cell chemical species
   integer, intent(in) :: num(:)

   !> Van-der-Waals radii
   real(wp), intent(in) :: rad(:)

   !> Descreened radii
   real(wp), intent(in) :: rho(:)

   !> Volume integral
   real(wp), intent(out) :: psi(:)

   !> Derivative of volume integral w.r.t. cartesian coordinates
   real(wp), intent(out) :: dpsidr(:, :, :)

   !> Derivative of volume integral w.r.t. strain deformations
   real(wp), intent(out) :: dpsidL(:, :, :)

   integer :: iat, ij, img, jat, izp, jzp
   logical :: tOvij, tOvji
   real(wp) :: vec(3), dist, rhoi, rhoj, rvdwi, rvdwj, weight
   real(wp) :: gi, gj, ap, am, lnab, rhab, ab, dgi, dgj
   real(wp) :: dG(3), dS(3, 3)
   real(wp) :: rh1, rhr1, r24, rh2, r1, aprh1, r12

   psi(:) = 0.0_wp
   dpsidr(:, :, :) = 0.0_wp
   dpsidL(:, :, :) = 0.0_wp

   do iat = 1, size(psi)
      izp = num(iat)
      do ij = 1, neighs(iat)
         img = neighList%ineigh(ij, iat)
         jat = neighList%image(img)
         jzp = num(jat)
         vec(:) = neighList%coords(:, iat) - neighList%coords(:, img)
         dist = sqrt(neighList%dist2(ij, iat))
         weight = neighList%weight(ij, iat)

         rhoi = rho(izp)
         rhoj = rho(jzp)
         rvdwi = rad(izp)
         rvdwj = rad(jzp)

         tOvij = dist < (rvdwi + rhoj)
         tOvji = dist < (rhoi + rvdwj)

         tOverlap: if (.not. tOvij .and. .not. tOvji) then ! ij do not overlap; ji do not overlap
            ! nonoverlaping spheres
            if (abs(rhoi-rhoj) < 1.e-8_wp) then
               ! equal reduced radii
               r1 = 1.0_wp/dist
               ap = dist+rhoj
               am = dist-rhoj
               ab = ap*am
               rhab = rhoj/ab
               lnab = 0.5_wp*log(am/ap)*r1
               gi = (rhab+lnab)*weight
               dgi = (-2.0_wp*rhab/ab+(rhab-lnab)*r1*r1)*weight
               ! accumulate psi
               psi(iat) = psi(iat) + gi
               psi(jat) = psi(jat) + gi
               ! accumulate psi gradient
               dG(:) = dgi*vec(:)
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
               dpsidr(:,iat,iat) = dpsidr(:,iat,iat) + dG(:)
               dpsidr(:,jat,iat) = dpsidr(:,jat,iat) - dG(:)
               dpsidr(:,jat,jat) = dpsidr(:,jat,jat) - dG(:)
               dpsidr(:,iat,jat) = dpsidr(:,iat,jat) + dG(:)
               dpsidL(:,:,iat) = dpsidL(:,:,iat) + dS
               if (iat /= jat) then
                  dpsidL(:,:,jat) = dpsidL(:,:,jat) + dS
               end if
            else
               ! unequal reduced radii
               ! ij contribution
               r1 = 1.0_wp/dist
               ap = dist+rhoj
               am = dist-rhoj
               ab = ap*am
               rhab = rhoj/ab
               lnab = 0.5_wp*log(am/ap)*r1
               gi = (rhab+lnab)*weight
               dgi = (-2.0_wp*rhab/ab+(rhab-lnab)*r1*r1)*weight
               ! ji contribution
               ap = dist+rhoi
               am = dist-rhoi
               ab = ap*am
               rhab = rhoi/ab
               lnab = 0.5_wp*log(am/ap)*r1
               gj = (rhab+lnab)*weight
               dgj = (-2.0_wp*rhab/ab+(rhab-lnab)*r1*r1)*weight
               ! accumulate psi
               psi(iat) = psi(iat) + gi
               psi(jat) = psi(jat) + gj
               ! accumulate psi gradient
               dG(:) = dgi*vec(:)
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
               dpsidr(:,iat,iat) = dpsidr(:,iat,iat) + dG(:)
               dpsidr(:,jat,iat) = dpsidr(:,jat,iat) - dG(:)
               dpsidL(:,:,iat) = dpsidL(:,:,iat) + dS

               dG(:) = dgj*vec(:)
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
               dpsidr(:,jat,jat) = dpsidr(:,jat,jat) - dG(:)
               dpsidr(:,iat,jat) = dpsidr(:,iat,jat) + dG(:)
               dpsidL(:,:,jat) = dpsidL(:,:,jat) + dS
            end if

         else if (.not. tOvij .and. tOvji) then tOverlap ! ij do not overlap; ji overlap

            ! ij contribution
            r1 = 1.0_wp/dist
            ap = dist+rhoj
            am = dist-rhoj
            ab = ap*am
            rhab = rhoj/ab
            lnab = 0.5_wp*log(am/ap)*r1
            gi = (rhab+lnab)*weight
            dgi = (-2.0_wp*rhab/ab+(rhab-lnab)*r1*r1)*weight
            ! accumulate psi
            psi(iat) = psi(iat) + gi
            ! accumulate psi gradient
            dG(:) = dgi*vec(:)
            dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
            dpsidr(:,iat,iat) = dpsidr(:,iat,iat) + dG(:)
            dpsidr(:,jat,iat) = dpsidr(:,jat,iat) - dG(:)
            dpsidL(:,:,iat) = dpsidL(:,:,iat) + dS

            if((dist+rhoi) > rvdwj) then
               ! ji contribution
               r1 = 1.0_wp/dist
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = dist+rhoi
               am = dist-rhoi
               rh1 = 1.0_wp/rvdwj
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gj = (rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab))*weight

               dgj = (rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  & rhoi*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1))*r1*weight
               ! accumulate psi
               psi(jat) = psi(jat) + gj
               ! accumulate psi gradient
               dG(:) = dgj*vec(:)
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
               dpsidr(:,jat,jat) = dpsidr(:,jat,jat) - dG(:)
               dpsidr(:,iat,jat) = dpsidr(:,iat,jat) + dG(:)
               dpsidL(:,:,jat) = dpsidL(:,:,jat) + dS
            end if

         else if (tOvij .and. .not. tOvji) then ! ij overlap; ji do not overlap

            if((dist+rhoj) > rvdwi) then
               ! ij contribution
               r1 = 1.0_wp/dist
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = dist+rhoj
               am = dist-rhoj
               rh1 = 1.0_wp/rvdwi
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gi = (rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab))*weight

               dgi = (rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  & rhoj*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1))*r1*weight
               ! accumulate psi
               psi(iat) = psi(iat) + gi
               ! accumulate psi gradient
               dG(:) = dgi*vec(:)
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
               dpsidr(:,iat,iat) = dpsidr(:,iat,iat) + dG(:)
               dpsidr(:,jat,iat) = dpsidr(:,jat,iat) - dG(:)
               dpsidL(:,:,iat) = dpsidL(:,:,iat) + dS
            end if

            ! ji contribution
            ap = dist+rhoi
            am = dist-rhoi
            ab = ap*am
            rhab = rhoi/ab
            lnab = 0.5_wp*log(am/ap)*r1
            gj = (rhab+lnab)*weight
            dgj = (-2.0_wp*rhab/ab+(rhab-lnab)*r1*r1)*weight
            ! accumulate psi
            psi(jat) = psi(jat) + gj
            ! accumulate psi gradient
            dG(:) = dgj*vec(:)
            dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
            dpsidr(:,jat,jat) = dpsidr(:,jat,jat) - dG(:)
            dpsidr(:,iat,jat) = dpsidr(:,iat,jat) + dG(:)
            dpsidL(:,:,jat) = dpsidL(:,:,jat) + dS

         else if (tOvij .and. tOvji) then tOverlap ! ij and ji overlap
            ! overlaping spheres
            if((dist+rhoj) > rvdwi) then
               ! ij contribution
               r1 = 1.0_wp/dist
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = dist+rhoj
               am = dist-rhoj
               rh1 = 1.0_wp/rvdwi
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gi = (rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab))*weight

               dgi = (rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  & rhoj*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1))*r1*weight
               ! accumulate psi
               psi(iat) = psi(iat) + gi
               ! accumulate psi gradient
               dG(:) = dgi*vec(:)
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
               dpsidr(:,iat,iat) = dpsidr(:,iat,iat) + dG(:)
               dpsidr(:,jat,iat) = dpsidr(:,jat,iat) - dG(:)
               dpsidL(:,:,iat) = dpsidL(:,:,iat) + dS
            end if

            if((dist+rhoi) > rvdwj) then
               ! ji contribution
               r1 = 1.0_wp/dist
               r12 = 0.5_wp*r1
               r24 = r12*r12

               ap = dist+rhoi
               am = dist-rhoi
               rh1 = 1.0_wp/rvdwj
               rhr1 = 1.0_wp/ap
               aprh1 = ap*rh1
               lnab = log(aprh1)

               gj = (rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab))*weight

               dgj = (rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
                  & rhoi*r24*(rhr1-rh1*aprh1)+r12*(r1*lnab-rhr1))*r1*weight
               ! accumulate psi
               psi(jat) = psi(jat) + gj
               ! accumulate psi gradient
               dG(:) = dgj*vec(:)
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
               dpsidr(:,jat,jat) = dpsidr(:,jat,jat) - dG(:)
               dpsidr(:,iat,jat) = dpsidr(:,iat,jat) + dG(:)
               dpsidL(:,:,jat) = dpsidL(:,:,jat) + dS
            end if

         end if tOverlap

      end do
   end do

end subroutine getPsi


end module xtb_solv_born
