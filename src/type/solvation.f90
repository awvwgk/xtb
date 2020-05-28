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

!> Interface definition for solvation models
module xtb_type_solvation
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_type_neighbourlist, only : TNeighbourList
   implicit none
   private

   public :: TSolvation


   !> Generic wrapper of for a solvation model
   type, abstract :: TSolvation
   contains

      !> Update geometry dependent data
      procedure(update), deferred :: update

      !> Add potential shift
      procedure(addShift), deferred :: addShift

      !> Get energy contributions
      procedure(getEnergy), deferred :: getEnergy

      !> Add contributions to the gradient
      procedure(getGradient), deferred :: getGradient

      !> Return realspace cutoff for the generation of neighbour lists
      procedure(getCutoff), deferred :: getCutoff

   end type TSolvation


   abstract interface
      !> Update solvation model from current geometry
      subroutine update(self, neighList, num)
         import :: TSolvation, TNeighbourList

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Neighbourlist
         class(TNeighbourList), intent(in) :: neighList

         !> Atomic identifiers
         integer, intent(in) :: num(:)
      end subroutine update

      !> Add potential shift from solvation model
      subroutine addShift(self, qat, qsh, atomicShift, shellShift)
         import :: TSolvation, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)

         !> Shell-resolved partial charges
         real(wp), intent(in) :: qsh(:)

         !> Atomic potential shift
         real(wp), intent(inout) :: atomicShift(:)

         !> Shell-resolved potential shift
         real(wp), intent(inout) :: shellShift(:)

      end subroutine addShift

      !> Get energy from solvation model
      pure subroutine getEnergy(self, qat, qsh, energy)
         import :: TSolvation, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)

         !> Shell-resolved partial charges
         real(wp), intent(in) :: qsh(:)

         !> Third order contribution to the energy
         real(wp), intent(out) :: energy

      end subroutine getEnergy

      !> Get gradient contributions from solvation model
      subroutine getGradient(self, neighList, num, qat, qsh, gradient, sigma)
         import :: TSolvation, TNeighbourList, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Neighbourlist
         class(TNeighbourList), intent(in) :: neighList

         !> Atomic identifiers
         integer, intent(in) :: num(:)

         !> Partial charges
         real(wp), intent(in) :: qat(:)

         !> Partial charges
         real(wp), intent(in) :: qsh(:)

         !> Molecular Gradient
         real(wp), intent(inout) :: gradient(:, :)

         !> Strain derivatives
         real(wp), intent(inout) :: sigma(:, :)

      end subroutine getGradient

      !> Return realspace cutoff for the generation of neighbour lists
      pure function getCutoff(self) result(cutoff)
         import :: TSolvation, wp

         !> Instance of the solvation model
         class(TSolvation), intent(in) :: self

         !> Maximal needed real space cutoff
         real(wp) :: cutoff

      end function getCutoff
   end interface


end module xtb_type_solvation
