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

!> Provides access to the solvation model class
module xtb_api_solvation
   use, intrinsic :: iso_c_binding
   use xtb_mctc_accuracy, only : wp
   use xtb_api_environment
   use xtb_api_utils
   use xtb_solv_input
   use xtb_solv_model
   use xtb_type_environment
   implicit none
   private

   public :: VSolvation
   public :: newSolvation_api, delSolvation_api


   !> Void pointer to solvation model blueprint
   type :: VSolvation
      type(TSolvModel) :: ptr
   end type VSolvation


contains


function newSolvation_api() result(vsolv) &
      & bind(C, name="xtb_newSolvation")
   type(VSolvation), pointer :: solv
   type(c_ptr) :: vsolv

   call checkGlobalEnv

   allocate(solv)
   vsolv = c_loc(solv)

end function newSolvation_api


subroutine delSolvation_api(vsolv) &
      & bind(C, name="xtb_delSolvation")
   type(c_ptr), intent(inout) :: vsolv
   type(VSolvation), pointer :: solv

   call checkGlobalEnv

   if (c_associated(vsolv)) then
      call c_f_pointer(vsolv, solv)
      deallocate(solv)
      vsolv = c_null_ptr
   end if

end subroutine delSolvation_api


end module xtb_api_solvation
