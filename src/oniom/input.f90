! This file is part of xtb.
!
! Copyright (C) 2020 Sebastian Ehlert
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

!> Input data to setup an ONIOM calculator stack
module xtb_oniom_input
   use xtb_mctc_accuracy, only : wp
   use xtb_main_methods, only : calcMethod
   use xtb_type_atomlist, only : TAtomlist
   implicit none
   private

   public :: TOniomInput, validOuterMethod


   !> Valid method for the outer region, method should be able to run parallel
   integer, parameter :: validOuterMethod(*) = &
      & [calcMethod%xtb, calcMethod%eht, calcMethod%gfnff]


   !> Information to setup an ONIOM calculator stack
   type :: TOniomInput

      !> Identifier for inner method
      integer :: inner = calcMethod%xtb

      !> Identifier for outer method
      integer :: outer = calcMethod%gfnff

      !> List of atoms for inner region
      type(TAtomlist) :: list

   end type TOniomInput


contains


end module xtb_oniom_input
