! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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

!> Docking implementation for xtb
module xtb_prog_dock
   use xtb_type_environment, only : TEnvironment, init
   use xtb_prog_argparser, only : TArgParser
   implicit none

contains


!> Entry point for performing non-covalent docking
subroutine xtbDock(env, argParser)

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "prog_dock"

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Command line arguments
   type(TArgParser), intent(inout) :: argParser

   call env%error("Not implemented", source)
   call env%checkpoint("Docking submodule")
end subroutine xtbDock


end module xtb_prog_dock
