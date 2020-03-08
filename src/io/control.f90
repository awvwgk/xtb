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

!> Reading and writing of Turbomole control like data structures
module xtb_io_control
   use xtb_mctc_accuracy, only : wp
   use xtb_type_abstractdata
   use xtb_type_cuckooroost
   use xtb_type_environment
   use xtb_type_linkedlist
   use xtb_type_reader
   implicit none
   private


   type :: TControlTokenEnum

      integer :: invalid = 0

      integer :: flag = 1

      integer :: equal = 2

      integer :: colon = 3

      integer :: comma = 4

      integer :: newline = 5

      integer :: whitespace = 6

      integer :: string = 7

   end type TControlTokenEnum

   type(TControlTokenEnum), parameter :: controlToken = TControlTokenEnum()

   character(len=*), parameter :: controlLetters = &
      & 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(len=*), parameter :: controlWhitespace = ' '//char(9)
   character(len=*), parameter :: controlDigits = '0123456789'
   character(len=*), parameter :: controlPunctuation = &
      & '!"%&'//achar(39)//'()*+-./<=>?@['//achar(92)//']^_`{|}~'
   character(len=*), parameter :: controlLiterals = &
      controlLetters//controlDigits//controlPunctuation


   type :: TControlToken

      integer :: tokenType = controlToken%invalid

      character(len=:), pointer :: ptr => null()

      integer :: len = 0

   end type TControlToken


   type :: TControlDeserializer

      logical :: finished

      type(TReader) :: reader

      type(TControlToken) :: token

      integer :: pos

      character(len=:), allocatable :: buffer

   end type TControlDeserializer


contains


subroutine parseRoot(de, env, ctrl)

   character(len=*), parameter :: source = 'io_control_parseRoot'

   type(TControlDeserializer), intent(inout) :: de

   type(TEnvironment), intent(inout) :: env

   type(TValue), intent(inout) :: ctrl

   class(TDict), pointer :: root

   character(len=:), allocatable :: message

   logical :: exitRun

   call ctrl%get(env, root)
   call env%check(exitRun)
   if (exitRun) then
      call env%error("Data structure is not a dictionary", source)
      return
   end if

   de%pos = 0
   de%buffer = ''
   de%token = newToken(controlToken%newline, de%buffer, 0)

   do while(.not.de%finished)
      call nextToken(de)
      select case(de%token%tokenType)
      case(controlToken%newline, controlToken%whitespace)
         cycle
      case(controlToken%flag)
         call parseGroup(de, env, root)
         call env%check(exitRun)
         if (exitRun) call env%error("Failed to parse group", source)
      case default
         call syntaxError(de, "Invalid characters, expected '$'", message)
         call env%error(message, source)
         exit
      end select
   end do

end subroutine parseRoot


subroutine parseGroup(de, env, root)

   character(len=*), parameter :: source = 'io_control_parseGroup'

   type(TControlDeserializer), intent(inout) :: de

   type(TEnvironment), intent(inout) :: env

   class(TDict), pointer, intent(in) :: root

   type(TValue), pointer :: group

   class(TDict), pointer :: dict

   character(len=:), allocatable :: key

   character(len=:), allocatable :: message

   integer :: groupType

   logical :: exitRun

   call nextToken(de)

   if (de%token%tokenType /= controlToken%string) then
      call syntaxError(de, "Invalid characters, expected string", message)
      call env%error(message, source)
      return
   end if

   call getKey(de%token, env, key)
   call root%insert(key, group)

   groupType = group%getType()
   if (any(groupType == [valueType%list, valueType%rArr, valueType%iArr])) then
      call env%error("Incompatible value already present", source)
      return
   end if

   do
      call nextToken(de)
      select case(de%token%tokenType)
      case(controlToken%whitespace); cycle
      case(controlToken%string)
         if (groupType == valueType%dict) then
            call env%error("Syntax error", source)
            return
         end if

      case(controlToken%newline)
         if (groupType == valueType%unknown) then
            call createDict(env, group, dict)
         end if
         if (groupType /= valueType%dict) then
            call env%error("Syntax error", source)
            return
         end if

         call parseKeyVal(de, env, dict)

      case default
         call env%error("Syntax error", source)
         return
      end select
   end do

end subroutine parseGroup


subroutine parseKeyVal(de, env, dict)

   character(len=*), parameter :: source = 'io_control_parseKeyVal'

   type(TControlDeserializer), intent(inout) :: de

   type(TEnvironment), intent(inout) :: env

   class(TDict), pointer, intent(in) :: dict

   character(len=:), allocatable :: key

   type(TValue), pointer :: val

   type(TLinkedList), pointer :: list

   logical :: exist, isArray

   lpKeyVal: do

      lpKey: do
         call nextToken(de)
         select case(de%token%tokenType)
         case(controlToken%flag)
            exit lpKeyVal

         case(controlToken%whitespace)
            cycle lpKey

         case(controlToken%string)
            call getKey(de%token, env, key)
            call dict%get(key, val)
            exist = associated(val)
            exit lpKey

         case default
            call env%error("Syntax error", source)
            return
         end select
      end do lpKey

      lpDelim: do
         call nextToken(de)
         select case(de%token%tokenType)
         case(controlToken%whitespace)
            cycle lpDelim

         case(controlToken%colon)
            isArray = .true.
            exit lpDelim

         case(controlToken%equal)
            isArray = .false.
            exit lpDelim

         case default
            call env%error("Syntax error", source)
            exit lpKeyVal

         end select
      end do lpDelim

      lpVal: do
         call nextToken(de)
         select case(de%token%tokenType)
         case(controlToken%whitespace)
            cycle lpVal

         case(controlToken%comma)
            ! TODO
            cycle lpVal

         case(controlToken%string)
            ! TODO
            cycle lpVal

         case(controlToken%newline)
            ! TODO
            exit lpVal

         case default
            call env%error("Syntax error", source)
            exit lpKeyVal
         end select
      end do lpVal

   end do lpKeyVal

end subroutine parseKeyval


subroutine getKey(token, env, key)

   type(TControlToken), intent(in) :: token

   type(TEnvironment), intent(inout) :: env

   character(len=:), allocatable, intent(out) :: key

   key = token%ptr(:token%len)

end subroutine getKey


subroutine createDict(env, ctrl, dict)

   character(len=*), parameter :: source = 'io_control_createDict'

   type(TEnvironment), intent(inout) :: env

   type(TValue), intent(inout) :: ctrl

   class(TDict), pointer, intent(out) :: dict

   type(TCuckooSkein), pointer :: skein

   logical :: exitRun

   allocate(skein)
   call init(skein, 20)

   call ctrl%set(env, skein)
   call env%check(exitRun)
   if (exitRun) then
      call env%error("Creation of dictionary failed", source)
      deallocate(skein)
      dict => null()
   else
      dict => skein
   end if


end subroutine createDict


subroutine nextToken(de)

   type(TControlDeserializer), intent(inout) :: de

   integer :: ii, err, last

   if (de%finished) return

   de%pos = de%pos + de%token%len
   if (de%pos>=len(de%buffer) .and. de%token%tokenType==controlToken%newline) then
      de%pos = 1
      call de%reader%read(de%buffer, err)
      de%finished = is_iostat_end(err)
   end if
   last = len(de%buffer)

   do while(de%pos < last)
      select case(de%buffer(de%pos:de%pos))

      case('#')
         de%pos = last
         cycle

      case('$')
         de%token = newToken(controlToken%flag, de%buffer(de%pos:), 1)
         return

      case('=')
         de%token = newToken(controlToken%equal, de%buffer(de%pos:), 1)
         return

      case(':')
         de%token = newToken(controlToken%colon, de%buffer(de%pos:), 1)
         return

      case(',', ';')
         de%token = newToken(controlToken%comma, de%buffer(de%pos:), 1)
         return

      case(' ', char(9))
         ii = verify(de%buffer(de%pos:), controlWhitespace)
         if (ii == 0) ii = last - de%pos + 1
         de%token = newToken(controlToken%whitespace, de%buffer(de%pos:), ii)
         return

      end select

      ii = verify(de%buffer(de%pos:), controlLiterals)
      if (ii == 0) ii = last - de%pos + 1
      de%token = newToken(controlToken%string, de%buffer(de%pos:), ii)
      return

   end do

   de%token = newToken(controlToken%newline, de%buffer(len(de%buffer):), 0)

end subroutine nextToken


function newToken(token, conf, len)

   integer, intent(in) :: token

   character(len=*), target, intent(in) :: conf

   integer, intent(in) :: len

   type(TControlToken) :: newToken

   newToken%tokenType = token
   newToken%ptr => conf(:len)
   newToken%len = len

end function newToken


subroutine syntaxError(de, message, error)

   type(TControlDeserializer), intent(inout) :: de

   character(len=*), intent(in) :: message

   character(len=:), allocatable, intent(out) :: error

   error = "Syntax Error:" // message // new_line('a') // de%buffer // &
      & new_line('a') // repeat('-', de%pos-1) // repeat('^', de%token%len)

end subroutine syntaxError


end module xtb_io_control
