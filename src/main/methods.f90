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

!> Available methods for computations
module xtb_main_methods
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: calcMethod, parseCalcMethod, writeCalcMethod


   !> Possible calculation methods
   type :: TCalcMethodEnum

      !> Invalid method
      integer :: invalid = -1

      !> Standard extended tight binding
      integer :: xtb = 1

      !> Extended Hückel
      integer :: eht = 0

      !> Quantum mechanical derived force field
      integer :: qmdff = 2

      !> Driving Turbomole
      integer :: turbomole = 4

      !> Driving Orca
      integer :: orca = 5

      !> Driving Mopac
      integer :: mopac = 12

      !> Geometry, frequency, non-covalent force field
      integer :: gfnff = 13

      !> Oniom calculator, meaning several of the above
      integer :: oniom = 10

   end type TCalcMethodEnum

   type(TCalcMethodEnum), parameter :: calcMethod = TCalcMethodEnum()


contains


!> Deserialize method identifier from a character representation
subroutine parseCalcMethod(method, string)

   !> Method identifier
   integer, intent(out) :: method

   !> String representation of the method identifier
   character(len=*), intent(in) :: string

   select case(string)
   case default
      method = calcMethod%invalid
   case('eht')
      method = calcMethod%eht
   case('xtb')
      method = calcMethod%xtb
   case('qmdff')
      method = calcMethod%qmdff
   case('orca')
      method = calcMethod%orca
   case('turbomole')
      method = calcMethod%turbomole
   case('mopac')
      method = calcMethod%mopac
   case('ff')
      method = calcMethod%gfnff
   case('oniom')
      method = calcMethod%oniom
   end select

end subroutine parseCalcMethod


!> Serialize method identifier to a character representation
subroutine writeCalcMethod(string, method)

   !> String representation of the method identifier
   character(len=:), allocatable, intent(out) :: string

   !> Method identifier
   integer, intent(in) :: method

   select case(method)
   case(calcMethod%eht)
      string = 'eht'
   case(calcMethod%xtb)
      string = 'xtb'
   case(calcMethod%qmdff)
      string = 'qmdff'
   case(calcMethod%orca)
      string = 'orca'
   case(calcMethod%turbomole)
      string = 'turbomole'
   case(calcMethod%mopac)
      string = 'mopac'
   case(calcMethod%gfnff)
      string = 'gfnff'
   case(calcMethod%oniom)
      string = 'oniom'
   end select

end subroutine writeCalcMethod


end module xtb_main_methods
