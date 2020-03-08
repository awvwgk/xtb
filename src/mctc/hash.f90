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
module xtb_mctc_hash
   use iso_c_binding
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: hash


   interface hash
      module procedure :: hashCharacter
   end interface hash


   !> Interface to C-implementation of hashing function
   interface hashFunction
      pure subroutine cMurmurHash3_x86_32(key, len, seed, out) &
            & bind(C, name="MurmurHash3_x86_32")
         import :: c_ptr, c_int
         type(c_ptr), value, intent(in) :: key
         integer(c_int), value, intent(in) :: len
         integer(c_int), value, intent(in) :: seed
         integer(c_int), intent(out) :: out
      end subroutine cMurmurHash3_x86_32
   end interface hashFunction


contains


pure function hashCharacter(var, seed) result(out)

   !> Variable to hash
   character(len=*), target, intent(in) :: var

   !> Seed for the hashing function
   integer, intent(in) :: seed

   !> Hash value
   integer :: out

   call hashFunction(c_loc(var), len(var), seed, out)

end function hashCharacter


end module xtb_mctc_hash
