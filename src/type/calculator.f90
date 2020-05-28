! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

!> abstract calculator that hides implementation details from calling codes
module xtb_type_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_neighbourlist, only : TNeighbourList
   use xtb_type_latticepoint, only : TLatticePoint
   use xtb_type_solvation, only : TSolvation
   use xtb_type_solvent, only : TSolvent
   use xtb_type_wavefunction, only : TWavefunction
   use xtb_type_wignerseitzcell, only : TWignerSeitzCell
   implicit none

   public :: TCalculator
   private


   !> Base calculator
   type, abstract :: TCalculator

      real(wp) :: accuracy
      type(TLatticePoint) :: latp
      type(TNeighbourList) :: neighList
      type(TWignerSeitzCell) :: wsCell
      class(TSolvation), allocatable :: solv_
      type(TSolvent), allocatable :: solv

   contains

      !> Perform single point calculation
      procedure(singlepoint), deferred :: singlepoint

      !> Write informative printout
      procedure(writeInfo), deferred :: writeInfo

      !> Update internal storage of coordinates
      procedure :: update

      !> Get realspace cutoff for users of neighbourlists
      procedure :: getCutoff

   end type TCalculator


   abstract interface
      subroutine singlepoint(self, env, mol, wfn, printlevel, restart, &
            & energy, gradient, sigma, hlgap, results)
         import :: TCalculator, TEnvironment, TMolecule, TWavefunction, wp
         import :: scc_results

         !> Calculator instance
         class(TCalculator), intent(inout) :: self

         !> Computational environment
         type(TEnvironment), intent(inout) :: env

         !> Molecular structure data
         type(TMolecule), intent(inout) :: mol

         !> Wavefunction data
         type(TWavefunction), intent(inout) :: wfn

         !> Print level for IO
         integer, intent(in) :: printlevel

         !> Restart from previous results
         logical, intent(in) :: restart

         !> Total energy
         real(wp), intent(out) :: energy

         !> Molecular gradient
         real(wp), intent(out) :: gradient(:, :)

         !> Strain derivatives
         real(wp), intent(out) :: sigma(:, :)

         !> HOMO-LUMO gap
         real(wp), intent(out) :: hlgap

         !> Detailed results
         type(scc_results), intent(out) :: results

      end subroutine singlepoint


      subroutine writeInfo(self, unit, mol)
         import :: TCalculator, TMolecule

         !> Calculator instance
         class(TCalculator), intent(in) :: self

         !> Unit for I/O
         integer, intent(in) :: unit

         !> Molecular structure data
         type(TMolecule), intent(in) :: mol

      end subroutine writeInfo
   end interface


contains


!> Update internal copy of real neighbour lists
subroutine update(self, env, mol, reset)

   !> Calculator instance
   class(TCalculator), intent(inout) :: self

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Force update
   logical, intent(in), optional :: reset

   real(wp), allocatable :: trans(:, :)
   real(wp) :: cutoff
   integer :: nTrans
   logical :: forced_update

   if (present(reset)) then
      forced_update = reset
   else
      forced_update = .false.
   end if

   cutoff = self%getCutoff()
   nTrans = self%latp%nTrans
   if (cutoff > 0.0_wp) then
      call self%latp%update(env, mol%lattice)
      call self%latp%getLatticePoints(trans, cutoff)
      if (nTrans == self%latp%nTrans .or. forced_update) then
         call self%neighList%update(mol%xyz, trans)
         call self%wsCell%update(mol%xyz, trans)
      else
         call self%neighList%generate(env, mol%xyz, cutoff, trans, .false.)
         call self%wsCell%generate(env, mol%xyz, cutoff, trans, .false.)
      end if
   end if

end subroutine update


!> Return realspace cutoff for the generation of neighbour lists
pure function getCutoff(self) result(cutoff)

   !> Calculator instance
   class(TCalculator), intent(in) :: self

   !> Maximal needed real space cutoff
   real(wp) :: cutoff

   if (allocated(self%solv_)) then
      cutoff = self%solv_%getCutoff()
   else
      cutoff = 0.0_wp
   end if

end function getCutoff


end module xtb_type_calculator
