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

!> Implementation of repulsion functions in xTB
module xtb_xtb_repulsion
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_neighbourlist, only : TNeighbourList
   use xtb_type_param, only : scc_parameter
   use xtb_aoparam, only : rep, en
   implicit none
   private

   public :: repulsionEnGrad0, repulsionEnGrad1, repulsionEnGrad2


contains


!> Repulsion gradient of GFN0-xTB
subroutine repulsionEnGrad0(mol, neighs, neighList, param, energy, gradient, sigma)
   type(scc_parameter), intent(in) :: param
   type(TMolecule), intent(in) :: mol
   integer, intent(in) :: neighs(:)
   class(TNeighbourList), intent(in) :: neighList
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(out) :: energy

   integer :: iat, jat, ati, atj, ij, img
   real(wp) :: r2, r1, r2top34, r2top34top2, rij(3)
   real(wp) :: den2, den4, alpha, repab, expterm, dtmp
   real(wp), allocatable :: energies(:)

   energy = 0.0_wp
   allocate(energies(len(mol)), source=0.0_wp)

   do iat = 1, mol%n
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         img = neighList%iNeigh(ij, iat)
         jat = neighList%image(img)
         rij(:) = neighList%coords(:, iat) - neighList%coords(:, img)
         r2 = neighList%dist2(ij, iat)
         atj = mol%at(jat)
         r1 = sqrt(r2)
         den2 = (en(ati) - en(atj))**2
         den4 = den2**2
         alpha = sqrt(rep(1,ati)*rep(1,atj))&
            & * (1.0_wp+(0.01_wp*den2+0.01_wp*den4)*param%xbrad)
         repab = rep(2,ati)*rep(2,atj)
         r2top34 = r2**0.75_wp
         r2top34top2 = r2top34**2
         expterm = exp(-alpha*r2top34)*repab
         dtmp = expterm*(1.5_wp*alpha*r2top34 + 1)/r2top34top2
         energies(iat) = energies(iat) + 0.5_wp*expterm/r1
         if (iat == jat) then
            sigma = sigma - 0.5_wp*dtmp*spread(rij,1,3)*spread(rij,2,3)
         else
            energies(jat) = energies(jat) + 0.5_wp*expterm/r1
            gradient(:, iat) = gradient(:, iat) - dtmp*rij
            gradient(:, jat) = gradient(:, jat) + dtmp*rij
            sigma = sigma - dtmp*spread(rij,1,3)*spread(rij,2,3)
         end if
      enddo
   enddo

   energy = sum(energies)

end subroutine repulsionEnGrad0


!> Repulsion gradient of GFN1-xTB
subroutine repulsionEnGrad1(mol, neighs, neighlist, kexp, rexp, &
      &                  energy, gradient, sigma)
   type(TMolecule), intent(in) :: mol
   class(TNeighbourList), intent(in) :: neighlist
   integer, intent(in) :: neighs(:)
   real(wp), intent(inout) :: gradient(:, :)
   real(wp), intent(inout) :: sigma(:, :)
   real(wp), intent(out) :: energy
   real(wp), intent(in) :: kexp
   real(wp), intent(in) :: rexp

   integer  :: iat, jat, ati, atj, ij, img
   real(wp) :: t16, t26, t27
   real(wp) :: alpha, repab
   real(wp) :: r1, r2, rij(3), dS(3, 3), dG(3), dE
   real(wp), allocatable :: energies(:)

   energy = 0.0_wp
   allocate(energies(len(mol)), source=0.0_wp)

   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, sigma) &
   !$omp shared(mol, neighs, neighlist, rep, kexp, rexp)&
   !$omp private(ij, img, jat, ati, atj, r2, rij, r1, alpha, repab, &
   !$omp&        t16, t26, t27, dE, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dist2(ij, iat)
         rij = mol%xyz(:, iat) - neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         r1 = sqrt(r2)

         alpha = sqrt(rep(1, ati)*rep(1, atj))
         repab = rep(2, ati)*rep(2, atj)
         t16 = r1**kexp
         t26 = exp(-alpha*t16)
         t27 = r1**rexp
         dE = repab * t26/t27
         dG = -(alpha*t16*kexp + rexp) * dE * rij/r2
         dS = spread(dG, 1, 3) * spread(rij, 2, 3)
         energies(iat) = energies(iat) + 0.5_wp * dE
         sigma = sigma + 0.5_wp * dS
         if (iat /= jat) then
            energies(jat) = energies(jat) + 0.5_wp * dE
            sigma = sigma + 0.5_wp * dS
         endif
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
      enddo
   enddo
   !$omp end parallel do

   energy = sum(energies)

end subroutine repulsionEnGrad1


!> Repulsion gradient of GFN2-xTB
subroutine repulsionEnGrad2(mol, neighs, neighlist, rexp, energy, gradient, sigma)
   use xtb_aoparam, only : rep
   use xtb_type_molecule
   use xtb_type_neighbourlist
   type(TMolecule), intent(in) :: mol
   class(TNeighbourlist), intent(in) :: neighlist
   integer, intent(in) :: neighs(:)
   real(wp),intent(inout) :: gradient(:, :)
   real(wp),intent(inout) :: sigma(:, :)
   real(wp),intent(out) :: energy
   real(wp),intent(in) :: rexp

   integer  :: iat, jat, ati, atj, ij, img
   real(wp) :: t16, t26, t27
   real(wp) :: alpha, repab
   real(wp) :: r1, r2, rij(3), dE, dG(3), dS(3, 3)
   real(wp), allocatable :: energies(:)

   energy = 0.0_wp
   allocate(energies(len(mol)), source=0.0_wp)

   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, sigma) &
   !$omp shared(mol, neighs, neighlist, rep, rexp)&
   !$omp private(ij, img, jat, ati, atj, r2, rij, r1, alpha, repab, &
   !$omp&        t16, t26, t27, dE, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dist2(ij, iat)
         rij = mol%xyz(:, iat) - neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         r1 = sqrt(r2)

         alpha = sqrt(rep(1,ati)*rep(1,atj))
         repab = rep(2,ati)*rep(2,atj)
         t16 = r1**kexp(ati, atj)
         t26 = exp(-alpha*t16)
         t27 = r1**rexp
         dE = repab * t26/t27
         dG = -(alpha*t16*kexp(ati, atj) + rexp) * dE * rij / r2
         dS = spread(dG, 1, 3) * spread(rij, 2, 3)
         energies(iat) = energies(iat) + 0.5_wp * dE
         sigma = sigma + 0.5_wp * dS
         if (iat /= jat) then
            energies(jat) = energies(jat) + 0.5_wp * dE
            sigma = sigma + 0.5_wp * dS
         endif
         gradient(:, iat) = gradient(:, iat) + dG
         gradient(:, jat) = gradient(:, jat) - dG
      enddo
   enddo
   !$omp end parallel do

   energy = sum(energies)

contains

real(wp) pure elemental function kexp(ati, atj)
   integer, intent(in) :: ati, atj
   if(ati <= 2 .and. atj <= 2) then
      kexp = 1.0_wp
   else
      kexp = 1.5_wp
   endif
end function kexp

end subroutine repulsionEnGrad2


end module xtb_xtb_repulsion
