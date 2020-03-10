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

!> Evalulation of Coulomb like terms in xTB
module xtb_xtb_coulomb
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi, sqrtpi
   use xtb_mctc_boundaryconditions, only : boundaryCondition
   use xtb_type_molecule, only : TMolecule, len
   implicit none
   private

   public :: get_gfn_coulomb_matrix, get_gfn_coulomb_derivs


   type :: enum_gam_average
      integer :: gfn1 = 1
      integer :: gfn2 = 2
   end type enum_gam_average

   type(enum_gam_average), parameter :: tb_gam_type = enum_gam_average()


   abstract interface
      real(wp) pure function gam_average(gi, gj) result(xij)
         import wp
         real(wp), intent(in) :: gi, gj
      end function gam_average
   end interface

contains


!> Wrapper for GFN-xTB electrostatics using the Mataga--Nishimoto--Ohno--Klopman
!  potential shape for the gamma-function.
!
!  This wrapper decides based on the systems periodicity which implementation
!  of the electrostatics to choose (molecular or Ewald summation), based on
!  the GFN-method the averaging function for the gamma-function is chosen.
pure subroutine get_gfn_coulomb_matrix(mol, nshell, ash, gam, gtype, cf, lqpc, jab)
   !> Molecular structure information.
   type(TMolecule), intent(in) :: mol

   !> Number of shells in the system.
   integer, intent(in) :: nshell

   !> Mapping from shells to atoms.
   integer, intent(in) :: ash(:)

   !> GFN-method identifying the averaging function.
   integer, intent(in) :: gtype

   !> Chemical hardness of every shell.
   real(wp), intent(in) :: gam(:)

   !> Convergence for the Ewald summation (only used under PBC).
   real(wp), intent(in) :: cf

   !> Quadrupole correction for Ewald summation (only used under PBC).
   logical, intent(in) :: lqpc

   !> Final Coulomb matrix over all shells.
   real(wp), intent(inout) :: jab(:, :)

   select case(gtype)
   case(tb_gam_type%gfn1)
      if (mol%npbc > 0) then
         call coulomb_matrix_3d_impl(mol, nshell, ash, gfn1_gam_average, gam, &
            &                        cf, lqpc, jab)
      else
         call coulomb_matrix_0d_impl(mol, nshell, ash, gfn1_gam_average, gam, jab)
      endif
   case(tb_gam_type%gfn2)
      if (mol%npbc > 0) then
         call coulomb_matrix_3d_impl(mol, nshell, ash, gfn2_gam_average, gam, &
            &                        cf, lqpc, jab)
      else
         call coulomb_matrix_0d_impl(mol, nshell, ash, gfn2_gam_average, gam, jab)
      endif
   end select

end subroutine get_gfn_coulomb_matrix

!> Average of chemical hardnesses used in GFN1-xTB gamma-function, returns 1/eta12.
real(wp) pure function gfn1_gam_average(gi, gj) result(xij)
   real(wp), intent(in) :: gi, gj
   xij = 0.5_wp*(1.0_wp/gi+1.0_wp/gj)
end function gfn1_gam_average

!> Average of chemical hardnesses used in GFN2-xTB gamma-function, returns 1/eta12.
real(wp) pure function gfn2_gam_average(gi, gj) result(xij)
   real(wp), intent(in) :: gi, gj
   xij = 2.0_wp/(gi+gj)
end function gfn2_gam_average

!> Implementation of Mataga--Nishimoto--Ohno--Klopman potential for molecular
!  systems.
pure subroutine coulomb_matrix_0d_impl(mol, nshell, ash, gav, gam, jab)

   !> Molecular structure information.
   type(TMolecule), intent(in) :: mol

   !> Number of shells in the system.
   integer, intent(in) :: nshell

   !> Mapping from shells to atoms.
   integer, intent(in) :: ash(:)

   !> Averaging function for chemical hardnesses.
   procedure(gam_average) :: gav

   !> Chemical hardness of every shell.
   real(wp), intent(in) :: gam(:)

   !> Final Coulomb matrix over all shells.
   real(wp), intent(inout) :: jab(:, :)

   integer :: is, iat, js, jat
   real(wp) :: gi, gj, xij, rij(3), r2

   do is = 1, nshell
      iat = ash(is)
      gi = gam(is)
      do js = 1, is-1
         jat = ash(js)
         gj = gam(js)
         xij = gav(gi, gj)
         rij = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = sum(rij**2)
         jab(js, is) = 1.0_wp/sqrt(r2 + xij**2)
         jab(is, js) = jab(js, is)
      enddo
      jab(is, is) = gi
   enddo

end subroutine coulomb_matrix_0d_impl

!> Implementation of Mataga--Nishimoto--Ohno--Klopman potential for systems
!  under 3D periodic boundary conditions. This Ewald summation explicitly
!  compensates higher multipole moments arising from the chosen potential shape.
pure subroutine coulomb_matrix_3d_impl(mol, nshell, ash, gav, gam, cf, lqpc, jab)

   !> Molecular structure information.
   type(TMolecule), intent(in) :: mol

   !> Number of shells in the system.
   integer, intent(in) :: nshell

   !> Mapping from shells to atoms.
   integer, intent(in) :: ash(:)

   !> Averaging function for chemical hardnesses.
   procedure(gam_average) :: gav

   !> Chemical hardness of every shell.
   real(wp), intent(in) :: gam(:)

   !> Convergence for the Ewald summation (only used under PBC).
   real(wp), intent(in) :: cf

   !> Quadrupole correction for Ewald summation (only used under PBC).
   logical, intent(in) :: lqpc

   !> Final Coulomb matrix over all shells.
   real(wp), intent(inout) :: jab(:, :)

   integer, parameter :: ewaldCutD(3) = 2
   integer, parameter :: ewaldCutR(3) = 2
   real(wp), parameter :: zero(3) = 0.0_wp

   integer :: is, iat, js, jat, img
   real(wp) :: qpc, gi, gj, xij, jc, jself
   real(wp) :: ri(3), rj(3), rw(3), riw(3)

   qpc = 0.0_wp
   jab = 0.0_wp
   do is = 1, nshell
      iat = ash(is)
      ri = mol%xyz(:, iat)
      gi = gam(is)
      do js = 1, is
         jat = ash(js)
         rj = mol%xyz(:, jat)
         gj = gam(js)
         xij = gav(gi, gj)
         if (lqpc) qpc = xij
         jc = 0.0_wp
         do img = 1, mol%wsc%itbl(jat, iat)
            rw = rj + matmul(mol%lattice, mol%wsc%lattr(:, img, jat, iat))
            riw = ri - rw
            jc = jc + mol%wsc%w(jat, iat) * (&
               & + gfn_ewald_3d_rec(riw,ewaldCutR,mol%rec_lat,qpc,mol%volume,cf) &
               & + gfn_ewald_3d_dir(riw,ewaldCutD,mol%lattice,xij,qpc,cf))
         enddo
         if (iat == jat) then
            jself = xij - 2.0_wp*cf/sqrtpi + 2.0_wp/3.0_wp*cf**3/sqrtpi*qpc**2
            jab(is, js) = jc + jself
            jab(js, is) = jc + jself
         else
            jab(is, js) = jc
            jab(js, is) = jc
         endif
      enddo
   enddo

end subroutine coulomb_matrix_3d_impl

pure function gfn_ewald_3d_dir(riw,rep,dlat,xij,qpc,cf) result(Amat)
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: dlat(3,3) !< direct lattice
   real(wp),intent(in) :: xij       !< interaction radius
   real(wp),intent(in) :: qpc       !< pseudo-quadrupole charge
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp) :: Amat                 !< element of interaction matrix
   real(wp),parameter :: eps = 1.0e-9_wp
   integer  :: dx,dy,dz
   real(wp) :: r2,r1,rij(3)
   real(wp) :: t(3),arg2,eij
   Amat = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), dy = -rep(2):rep(2), dz = -rep(3):rep(3))
      rij = riw + matmul(dlat, [dx,dy,dz])
      r2 = sum(rij**2)
      if(r2 < eps) cycle
      r1 = sqrt(r2)
      arg2 = cf**2*r2
      eij = erf(cf*r1)
      Amat = Amat + 1.0_wp/sqrt(r1**2 + xij**2) &
         & - (eij/r1 + 0.5_wp*qpc**2*(2*cf/sqrtpi*exp(-arg2)/r2 - eij/(r2*r1)))
   end do
end function gfn_ewald_3d_dir

pure function gfn_ewald_3d_rec(riw,rep,rlat,qpc,vol,cf) result(Amat)
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: rlat(3,3) !< reciprocal lattice
   real(wp),intent(in) :: qpc       !< pseudo-quadrupole charge
   real(wp),intent(in) :: vol       !< direct cell volume
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp) :: Amat                 !< element of interaction matrix
   integer  :: dx,dy,dz
   real(wp) :: rik2,rik(3)
   real(wp) :: t(3)
   real(wp) :: fpivol
   Amat = 0.0_wp
   fpivol = 4.0_wp*pi/vol
   do concurrent(dx = -rep(1):rep(1), dy = -rep(2):rep(2), dz = -rep(3):rep(3), &
         &       dx/=0 .or. dy/=0 .or. dz/=0)
      rik = matmul(rlat, [dx,dy,dz])
      rik2 = sum(rik**2)
      Amat = Amat + cos(dot_product(rik,riw)) &
         * exp(-rik2/(4.0_wp*cf**2))/rik2 * (1.0_wp + 2*rik2*qpc**2)
   end do
   Amat = Amat * fpivol
end function gfn_ewald_3d_rec

!! ========================================================================
!  shellwise electrostatic gradient for GFN1
!! ========================================================================
subroutine get_gfn_coulomb_derivs(mol, nshell, ash, gam, gtype, cf, lqpc, qsh, &
      &                           gradient, sigma)

   !> Molecular structure information.
   type(TMolecule), intent(in) :: mol

   !> Number of shells in the system.
   integer, intent(in) :: nshell

   !> Mapping from shells to atoms.
   integer, intent(in) :: ash(:)

   !> GFN-method identifying the averaging function.
   integer, intent(in) :: gtype

   !> Chemical hardness of every shell.
   real(wp), intent(in) :: gam(:)

   !> Convergence for the Ewald summation (only used under PBC).
   real(wp), intent(in) :: cf

   !> Quadrupole correction for Ewald summation (only used under PBC).
   logical, intent(in) :: lqpc

   real(wp), intent(in) :: qsh(:)

   real(wp), intent(inout) :: gradient(:, :)

   real(wp), intent(inout) :: sigma(:, :)

   select case(gtype)
   case(tb_gam_type%gfn1)
      if (mol%npbc > 0) then
         call coulomb_derivs_3d_impl(mol, nshell, ash, gam, gfn1_gam_average, &
            &                        cf, lqpc, qsh, gradient, sigma)
      else
         call coulomb_derivs_0d_impl(mol, nshell, ash, gam, gfn1_gam_average, &
            &                        qsh, gradient)
      endif
   case(tb_gam_type%gfn2)
      if (mol%npbc > 0) then
         call coulomb_derivs_3d_impl(mol, nshell, ash, gam, gfn2_gam_average, &
            &                        cf, lqpc, qsh, gradient, sigma)
      else
         call coulomb_derivs_0d_impl(mol, nshell, ash, gam, gfn2_gam_average, &
            &                        qsh, gradient)
      endif
   end select

end subroutine get_gfn_coulomb_derivs

subroutine coulomb_derivs_0d_impl(mol, nshell, ash, gam, gav, qsh, gradient)

   !> Molecular structure information.
   type(TMolecule), intent(in) :: mol

   !> Number of shells in the system.
   integer, intent(in) :: nshell

   !> Mapping from shells to atoms.
   integer, intent(in) :: ash(:)

   !> Averaging function for chemical hardnesses.
   procedure(gam_average) :: gav

   !> Chemical hardness of every shell.
   real(wp),intent(in) :: gam(:)

   real(wp),intent(in) :: qsh(:)

   real(wp),intent(inout) :: gradient(:,:)

   integer  :: is, js, iat, jat
   real(wp) :: gi, gj, r2, rij(3), xij, dG(3)

   do is = 1, nshell
      iat = ash(is)
      gi = gam(is)
      do js = 1, is-1
         jat = ash(js)
         if(jat >= iat) cycle
         rij = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = sum(rij**2)
         gj = gam(js)
         xij = gav(gi, gj)
         dG = -qsh(is)*qsh(js)/sqrt(r2 + xij**2)**3 * rij
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
      enddo
   enddo

end subroutine coulomb_derivs_0d_impl

subroutine coulomb_derivs_3d_impl(mol, nshell, ash, gam, gav, cf, lqpc, qsh, &
      &                           gradient, sigma)

   !> Molecular structure information.
   type(TMolecule), intent(in) :: mol

   !> Number of shells in the system.
   integer, intent(in) :: nshell

   !> Mapping from shells to atoms.
   integer, intent(in) :: ash(:)

   !> Averaging function for chemical hardnesses.
   procedure(gam_average) :: gav

   !> Chemical hardness of every shell.
   real(wp),intent(in) :: gam(:)

   !> Convergence for the Ewald summation (only used under PBC).
   real(wp),intent(in) :: cf

   !> Quadrupole correction for Ewald summation (only used under PBC).
   logical, intent(in) :: lqpc

   real(wp),intent(in) :: qsh(:)

   real(wp),intent(inout) :: gradient(:, :)

   real(wp),intent(inout) :: sigma(:, :)

   integer, parameter :: ewaldCutD(3) = 2
   integer, parameter :: ewaldCutR(3) = 2
   real(wp), parameter :: zero(3) = 0.0_wp

   integer :: is, js, iat, jat, img
   real(wp) :: gi, gj, r2, ri(3), rj(3), rw(3), riw(3), qpc, xij, wqq
   real(wp) :: dG(3), dS(3, 3), tG(3), tS(3, 3)

   qpc = 0.0_wp

   do is = 1, nshell
      iat = ash(is)
      ri = mol%xyz(:, iat)
      gi = gam(is)
      do js = 1, is
         jat = ash(js)
         rj = mol%xyz(:, jat)
         gj = gam(js)
         xij = gav(gi, gj)
         if (lqpc) qpc = xij
         dG = 0.0_wp
         dS = 0.0_wp
         do img = 1, mol%wsc%itbl(jat, iat)
            rw = rj + matmul(mol%lattice, mol%wsc%lattr(:, img, jat, iat))
            riw = ri - rw
            wqq = mol%wsc%w(jat,iat)*qsh(is)*qsh(js)
            call gfn_ewald_dx_3d_rec(riw, ewaldCutR, mol%rec_lat, qpc, &
               &                     mol%volume, cf, tG, tS)
            dG = dG + tG * wqq
            dS = dS + tS * wqq
            call gfn_ewald_dx_3d_dir(riw, ewaldCutD, mol%lattice, xij, qpc, cf, &
               &                     tG, tS)
            dG = dG + tG * wqq
            dS = dS + tS * wqq
         enddo
         if (is == js) then
            sigma = sigma + dS*0.5_wp
         else
            if (iat /= jat) then
               gradient(:, iat) = gradient(:, iat) + dG
               gradient(:, jat) = gradient(:, jat) - dG
            endif
            sigma = sigma + dS
         endif
      enddo
   enddo

end subroutine coulomb_derivs_3d_impl

pure subroutine gfn_ewald_dx_3d_rec(riw,rep,rlat,qpc,vol,cf,dG,dS)
   use mctc_constants
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: rlat(3,3) !< reciprocal lattice
   real(wp),intent(in) :: vol       !< direct cell volume
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp),intent(in) :: qpc       !< pseudo-quadrupole charge
   real(wp),intent(out) :: dG(3) !< element of interaction matrix
   real(wp),intent(out) :: dS(3,3)
   integer  :: dx,dy,dz
   real(wp) :: rik2,rik(3)
   real(wp) :: t(3),dtmp,fpivol
   real(wp) :: expterm,arg
   real(wp), parameter :: unity(3, 3) = reshape(&
      &[1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & shape(unity))
   dG = 0.0_wp
   dS = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), dy = -rep(2):rep(2), dz = -rep(3):rep(3), &
         &       dx/=0 .or. dy/=0 .or. dz/=0)
      rik = matmul(rlat, [dx,dy,dz])
      rik2 = sum(rik**2)
      expterm = exp(-rik2/(4.0_wp*cf**2))/rik2
      arg = dot_product(rik,riw)
      ! d/dx (sin**2 + cos**2) = -2*sin*cos - 2*cos*sin
      dtmp = -sin(arg)*expterm
      dG = dG + rik*dtmp
      dS = dS + cos(arg)*expterm * (spread(rik, 1, 3)*spread(rik, 2, 3) &
         & * (-2.0_wp/rik2 - 0.5_wp/cf**2 + 4*qpc**2) &
         & - unity*(1.0_wp + 2*rik2*qpc**2))
   end do
   fpivol = 4.0_wp*pi/vol
   dG = dG * fpivol
   dS = dS * fpivol
end subroutine gfn_ewald_dx_3d_rec

pure subroutine gfn_ewald_dx_3d_dir(riw,rep,dlat,xij,qpc,cf,dG,dS)
   use mctc_constants
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: dlat(3,3) !< direct lattice
   real(wp),intent(in) :: xij       !< interaction radius
   real(wp),intent(in) :: qpc       !< pseudo-quadrupole charge
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp),intent(out) :: dG(3) !< element of interaction matrix
   real(wp),intent(out) :: dS(3,3)
   real(wp),parameter :: eps = 1.0e-9_wp
   integer  :: dx,dy,dz
   real(wp) :: r2,r1,rij(3),arg2
   real(wp) :: t(3),dtmp,stmp(3)
   dG = 0.0_wp
   dS = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), dy = -rep(2):rep(2), dz = -rep(3):rep(3))
      ! real contributions
      t = [dx,dy,dz]
      rij = riw + matmul(dlat,t)
      r1 = norm2(rij)
      r2 = r1*r1
      if(r1 < eps) cycle ! self-interaction handled elsewhere
      arg2 = cf**2 * r2
      dtmp = - 2*cf*exp(-arg2)/(sqrtpi*r2) + erf(cf*r1)/(r2*r1) &
         &   + 0.5_wp*qpc**2 * ((4*cf**3*r2 + 6*cf)*exp(-arg2)/(sqrtpi*r2*r2) &
         &                      - 3*erf(cf*r1)/(r2*r2*r1)) &
         &   - 1.0_wp/sqrt(r1**2 + xij**2)**3
      dG = dG + dtmp*rij
      dS = dS + dtmp*spread(rij, 1, 3)*spread(rij, 2, 3)
   enddo

end subroutine gfn_ewald_dx_3d_dir

end module xtb_xtb_coulomb
