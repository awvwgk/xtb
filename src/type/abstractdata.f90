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

! --------------------------------------------------------------------------------
!> Derived types for working with abstract data structures.
!>
!> This module implements an additional abstraction layer around Fortran standard
!> types, standard types are assumed to be
!>
!>  - deferred-length characters
!>  - reals (wanted precision)
!>  - one-dimensional arrays of reals (wanted precision)
!>  - integers (default precision)
!>  - one-dimensional arrays of integers (default precision)
!>  - logicals (default precision)
!>
!> Additionally this module provides two abstract data types, which are
!>
!>  - dictionaries (either hash maps or binary search trees)
!>  - lists (either linked lists or arrays)
!>
!> Together with the empty type they form the nine base types of the abstract
!> value container implemented here. The value container should provide sufficient
!> getter and setter functions to manipulate its contained data.
!> The dictionary and list implementation have to provided externally by inheriting
!> from the abstract base class and providing all the deferred procedures.
!>
!> Since all the datastructures are opaque you cannot directly interact with them,
!> An abstract visitor is provided, any class inheriting from this visitor can
!> be accepted from the abstract data structures and work on them in a read-only
!> fashion.
!> Note that the visitor must be able to visit all the nine data types defined
!> in the abstract value container, the acceptor function of the value container
!> will redirect you to the appropiate data in the container. Therefore, you cannot
!> explicitly visit a value container, the visit function on the whold value
!> container is only called for empty values.
module xtb_type_abstractdata
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: TValue, TDict, TList, TVisitor, valueType, len


   !> Possible types of a cuckoo value
   type :: TValueTypeEnum

      !> Uninitialized value
      integer :: unknown = 0

      !> Character value
      integer :: cVal = 1

      !> Real value
      integer :: rVal = 2

      !> Real array value
      integer :: rArr = 3

      !> Integer value
      integer :: iVal = 4

      !> Integer array value
      integer :: iArr = 5

      !> Logical value
      integer :: lVal = 6

      !> Dictionary value
      integer :: dict = 7

      !> List value
      integer :: list = 8

   end type TValueTypeEnum

   !> Actual value type enumerator
   type(TValueTypeEnum), parameter :: valueType = TValueTypeEnum()


   !> Dictionary mapping character keys to values
   type, abstract :: TDict
   contains

      !> Get length of the data type
      procedure(lenD), deferred :: getLength

      !> Insert a new key into the dictionary
      procedure(insert), deferred :: insert

      !> Get an existing key from the dictionary
      procedure(getD), deferred :: get

      !> Drop a key from the dictionary
      procedure(drop), deferred :: drop

      !> Export keys as list
      procedure(export), deferred :: export

      !> Accept a visitor to transverse the data structure
      procedure :: accept => acceptDict

      !> Clear the dictionary
      procedure(clearDict), deferred :: clear

   end type TDict


   !> List of values
   type, abstract :: TList
   contains

      !> Get length of the data type
      procedure(lenL), deferred :: getLength

      !> Get pointer value at this index
      procedure(getL), deferred :: get

      !> Find a index of a value
      procedure(find), deferred :: find

      !> Create a new value and return pointer to this entry
      procedure(append), deferred :: append

      !> Accept a visitor to transverse the data structure
      procedure :: accept => acceptList

      !> Clear the list
      procedure(clearList), deferred :: clear

   end type TList


   !> Abstract value container
   type :: TValue
      private

      !> Value for character data structure
      character(len=:), allocatable :: cVal

      !> Value for real data structure
      real(wp), allocatable :: rVal

      !> Array for real data structure
      real(wp), allocatable :: rArr(:)

      !> Value for integer data structure
      integer, allocatable :: iVal

      !> Array for integer data structure
      integer, allocatable :: iArr(:)

      !> Value for logical data structure
      logical, allocatable :: lVal

      !> Dictionary for the data structure
      class(TDict), pointer :: dict => null()

      !> List for the data structure
      class(TList), pointer :: list => null()

      !> Type of this data structure
      integer :: valueType = valueType%unknown

   contains

      !> Overloaded equals operator
      generic :: operator(==) => equalsValue

      !> Compare to values
      procedure, private :: equalsValue

      !> Get value information
      procedure :: getType

      !> Setter interface for the data structure
      generic :: set => setCharacter, setReal, setRealArray, setInteger, &
         & setIntegerArray, setLogical, setDict, setList

      !> Set a character variable in the data structure
      procedure, private :: setCharacter

      !> Set a real variable in the data structure
      procedure, private :: setReal

      !> Set a real array in the data structure
      procedure, private :: setRealArray

      !> Set an integer variable in the data structure
      procedure, private :: setInteger

      !> Set an integer array in the data structure
      procedure, private :: setIntegerArray

      !> Set a logical variable in the data structure
      procedure, private :: setLogical

      !> Set a dictionary in the data structure
      procedure, private :: setDict

      !> Set a list in the data structure
      procedure, private :: setList

      !> Getter interface for the data structure
      generic :: get => getCharacter, getReal, getRealArray, getInteger, &
         & getIntegerArray, getLogical, getDict, getList

      !> Get a character variable from the data structure
      procedure, private :: getCharacter

      !> Get a real variable from the data structure
      procedure, private :: getReal

      !> Get a real array from the data structure
      procedure, private :: getRealArray

      !> Get an integer variable from the data structure
      procedure, private :: getInteger

      !> Get an integer arrays from the data structure
      procedure, private :: getIntegerArray

      !> Get a logical variable from the data structure
      procedure, private :: getLogical

      !> Get reference to a dictionary from the data structure
      procedure, private :: getDict

      !> Get reference to a list from the data structure
      procedure, private :: getList

      !> Accept a visitor to transverse the data structure
      procedure :: accept

      !> Clear value form data structure
      procedure :: clear => clearValue

   end type


   !> Abstract visitor for abstract data structures
   type, abstract :: TVisitor
   contains

      !> Visitor entering a node of the data structure
      generic :: visit => visitValue, visitCharacter, visitReal, visitRealArray, &
         & visitInteger, visitIntegerArray, visitLogical, visitDict, visitList

      !> Visitor entering a node of the data structure
      procedure(visitValue), deferred :: visitValue

      !> Visitor entering a character variable of the data structure
      procedure(visitCharacter), deferred :: visitCharacter

      !> Visitor entering a real variable of the data structure
      procedure(visitReal), deferred :: visitReal

      !> Visitor entering a real array of the data structure
      procedure(visitRealArray), deferred :: visitRealArray

      !> Visitor entering an integer variable of the data structure
      procedure(visitInteger), deferred :: visitInteger

      !> Visitor entering an integer array of the data structure
      procedure(visitIntegerArray), deferred :: visitIntegerArray

      !> Visitor entering a logical variable of the data structure
      procedure(visitLogical), deferred :: visitLogical

      !> Visitor entering a dictionary of the data structure
      procedure(visitDict), deferred :: visitDict

      !> Visitor entering a list of the data structure
      procedure(visitList), deferred :: visitList

   end type TVisitor


   !> Get length of a data structure
   interface len
      module procedure :: lenDict
      module procedure :: lenList
   end interface


   abstract interface
   !> Return length of a hash map
   pure function lenD(self) result(length)
      import :: TDict

      !> Instance of the data structure
      class(TDict), intent(in) :: self

      !> Length of the data structure
      integer :: length

   end function lenD

   !> Insert a new key into the hash map
   recursive subroutine insert(self, key, ptr)
      import TDict, TValue

      !> Instance of the hash map
      class(TDict), target, intent(inout) :: self

      !> Character key to insert
      character(len=*), intent(in) :: key

      !> Pointer for the newly occupied nest
      type(TValue), pointer, intent(out) :: ptr

   end subroutine insert

   !> Get a pointer to the value with the corresponding character key
   subroutine getD(self, key, ptr)
      import :: TDict, TValue

      !> Instance of the hash map
      class(TDict), target, intent(inout) :: self

      !> Character key
      character(len=*), intent(in) :: key

      !> Pointer to the value
      type(TValue), pointer, intent(out) :: ptr

   end subroutine getD

   !> Drop a value from the hash map
   subroutine drop(self, key)
      import :: TDict

      !> Instance of the hash map
      class(TDict), intent(inout) :: self

      !> Character key
      character(len=*), intent(in) :: key

   end subroutine drop

   !> Export keys as list
   subroutine export(self, env, list)
      import :: TDict, TEnvironment, TList

      !> Instance of the hash map
      class(TDict), intent(inout) :: self

      !> Computation environment
      type(TEnvironment), intent(inout) :: env

      !> Instance of the key list
      class(TList), intent(inout) :: list

   end subroutine export

   !> Clear the hash map
   subroutine clearDict(self)
      import :: TDict

      !> Instance of the hash map
      class(TDict), intent(inout) :: self

   end subroutine clearDict

   !> Return length of a list
   pure function lenL(self) result(length)
      import :: TList

      !> Instance of the data structure
      class(TList), intent(in) :: self

      !> Length of the data structure
      integer :: length

   end function lenL

   !> Get pointer value at this index
   function getL(self, ind) result(val)
      import :: TList, TValue

      !> Instance of the linked list
      class(TList), intent(inout) :: self

      !> Index to lookup
      integer, intent(in) :: ind

      !> Pointer to the content of the node
      type(TValue), pointer :: val

   end function getL

   !> Find a index of a value
   function find(self, val) result(ind)
      import :: TList, TValue

      !> Instance of the linked list
      class(TList), intent(inout) :: self

      !> Value to search for
      type(TValue), intent(in) :: val

      !> Index of the value
      integer :: ind

   end function find

   !> Create a new value and return pointer to this entry
   subroutine append(self, ptr)
      import :: TList, TValue

      !> Instance of the linked list
      class(TList), intent(inout) :: self

      !> Pointer to newly created element
      type(TValue), pointer, intent(out) :: ptr

   end subroutine append

   !> Clear the list
   subroutine clearList(self)
      import :: TList

      !> Instance of the list
      class(TList), intent(inout) :: self

   end subroutine clearList

   !> Visitor entering a node of the data structure
   subroutine visitValue(visitor, val)
      import TVisitor, TValue

      !> Visitor for this value
      class(TVisitor), intent(inout) :: visitor

      !> Instance of the data structure.
      class(TValue), intent(in) :: val

   end subroutine visitValue

   !> Visitor entering a node of the data structure
   subroutine visitCharacter(visitor, val)
      import TVisitor

      !> Visitor for this value
      class(TVisitor), intent(inout) :: visitor

      !> Value to vist
      character(len=*), intent(in) :: val

   end subroutine visitCharacter

   !> Visitor entering a node of the data structure
   subroutine visitReal(visitor, val)
      import TVisitor, wp

      !> Visitor for this value
      class(TVisitor), intent(inout) :: visitor

      !> Value to visit
      real(wp), intent(in) :: val

   end subroutine visitReal

   !> Visitor entering a node of the data structure
   subroutine visitRealArray(visitor, val)
      import TVisitor, wp

      !> Visitor for this value
      class(TVisitor), intent(inout) :: visitor

      !> Value to visit
      real(wp), intent(in) :: val(:)

   end subroutine visitRealArray

   !> Visitor entering a node of the data structure
   subroutine visitInteger(visitor, val)
      import TVisitor

      !> Visitor for this value
      class(TVisitor), intent(inout) :: visitor

      !> Value to visit
      integer, intent(in) :: val

   end subroutine visitInteger

   !> Visitor entering a node of the data structure
   subroutine visitIntegerArray(visitor, val)
      import TVisitor

      !> Visitor for this value
      class(TVisitor), intent(inout) :: visitor

      !> Value to visit
      integer, intent(in) :: val(:)

   end subroutine visitIntegerArray

   !> Visitor entering a node of the data structure
   subroutine visitLogical(visitor, val)
      import TVisitor

      !> Visitor for this value
      class(TVisitor), intent(inout) :: visitor

      !> Value to visit
      logical, intent(in) :: val

   end subroutine visitLogical

   !> Visitor entering a node of the data structure
   subroutine visitDict(visitor, val)
      import TVisitor, TDict

      !> Visitor for this value
      class(TVisitor), intent(inout) :: visitor

      !> Value to visit
      class(TDict), intent(in) :: val

   end subroutine visitDict

   !> Visitor entering a node of the data structure
   subroutine visitList(visitor, val)
      import TVisitor, TList

      !> Visitor for this value
      class(TVisitor), intent(inout) :: visitor

      !> Value to visit
      class(TList), intent(in) :: val

   end subroutine visitList
   end interface


contains


!> Return length of a hash map
function lenDict(self) result(length)

   !> Instance of the data structure
   class(TDict), intent(in) :: self

   !> Length of the data structure
   integer :: length

   length = self%getLength()

end function lenDict


!> Return length of a list
function lenList(self) result(length)

   !> Instance of the data structure
   class(TList), intent(in) :: self

   !> Length of the data structure
   integer :: length

   length = self%getLength()

end function lenList


!> Clear a cuckoo nest
subroutine clearValue(self)

   !> Cuckoo nest instance
   class(TValue), intent(inout) :: self

   select case(self%valueType)
   case(valueType%cVal); deallocate(self%cVal)
   case(valueType%rVal); deallocate(self%rVal)
   case(valueType%rArr); deallocate(self%rArr)
   case(valueType%iVal); deallocate(self%iVal)
   case(valueType%iArr); deallocate(self%iArr)
   case(valueType%lVal); deallocate(self%lVal)
   case(valueType%dict); nullify(self%dict)
   case(valueType%list); nullify(self%list)
   end select
   self%valueType = valueType%unknown

end subroutine clearValue


!> Set a character variable in the data structure
subroutine setCharacter(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_setCharacter"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to set
   character(len=*), intent(in) :: val

   if (self%valueType /= valueType%cVal) then
      if (self%valueType /= valueType%unknown) then
         call self%clear
         call env%warning("Value type changed", source)
      end if
      self%valueType = valueType%cVal
   end if
   self%cVal = val

end subroutine setCharacter


!> Set a real variable in the data structure
subroutine setReal(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_setReal"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to set
   real(wp), intent(in) :: val

   if (self%valueType /= valueType%rVal) then
      if (self%valueType /= valueType%unknown) then
         call self%clear
         call env%warning("Value type changed", source)
      end if
      self%valueType = valueType%rVal
   end if
   self%rVal = val

end subroutine setReal


!> Set a real array in the data structure
subroutine setRealArray(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_setRealArray"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to set
   real(wp), intent(in) :: val(:)

   if (self%valueType /= valueType%rArr) then
      if (self%valueType /= valueType%unknown) then
         call self%clear
         call env%warning("Value type changed", source)
      end if
      self%valueType = valueType%rArr
   end if
   self%rArr = val

end subroutine setRealArray


!> Set an integer variable in the data structure
subroutine setInteger(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_setInteger"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to set
   integer, intent(in) :: val

   if (self%valueType /= valueType%iVal) then
      if (self%valueType /= valueType%unknown) then
         call self%clear
         call env%warning("Value type changed", source)
      end if
      self%valueType = valueType%iVal
   end if
   self%iVal = val

end subroutine setInteger


!> Set an integer array in the data structure
subroutine setIntegerArray(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_setIntegerArray"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to set
   integer, intent(in) :: val(:)

   if (self%valueType /= valueType%iArr) then
      if (self%valueType /= valueType%unknown) then
         call self%clear
         call env%warning("Value type changed", source)
      end if
      self%valueType = valueType%iArr
   end if
   self%iArr = val

end subroutine setIntegerArray


!> Set a logical variable in the data structure
subroutine setLogical(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_setLogical"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to set
   logical, intent(in) :: val

   if (self%valueType /= valueType%lVal) then
      if (self%valueType /= valueType%unknown) then
         call self%clear
         call env%warning("Value type changed", source)
      end if
      self%valueType = valueType%lVal
   end if
   self%lVal = val

end subroutine setLogical


!> Set a dictionary in the data structure
subroutine setDict(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_setDict"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to set
   class(TDict), target, intent(in) :: val

   if (self%valueType /= valueType%dict) then
      if (self%valueType /= valueType%unknown) then
         call self%clear
         call env%warning("Value type changed", source)
      end if
      self%valueType = valueType%dict
   end if
   self%dict => val

end subroutine setDict


!> Set a list in the data structure
subroutine setList(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_setList"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to set
   class(TList), target, intent(in) :: val

   if (self%valueType /= valueType%list) then
      if (self%valueType /= valueType%unknown) then
         call self%clear
         call env%warning("Value type changed", source)
      end if
      self%valueType = valueType%list
   end if
   self%list => val

end subroutine setList


!> Get a character variable from the data structure
subroutine getCharacter(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_getCharacter"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to get
   character(len=:), allocatable, intent(out) :: val

   if (self%valueType == valueType%cVal) then
      val = self%cVal
   else
      call env%error("Data structure does not contain character variable", source)
   end if

end subroutine getCharacter


!> Get a real variable from the data structure
subroutine getReal(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_getReal"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to get
   real(wp), intent(out) :: val

   if (self%valueType == valueType%rVal) then
      val = self%rVal
   else
      call env%error("Data structure does not contain real variable", source)
   end if

end subroutine getReal


!> Get a real array from the data structure
subroutine getRealArray(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_getRealArray"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to get
   real(wp), intent(out) :: val(:)

   if (self%valueType == valueType%rArr) then
      if (size(val) >= size(self%rArr)) then
         val(:) = self%rArr
      else
         call env%error("Insufficient size for real array", source)
      end if
   else
      call env%error("Data structure does not contain real array", source)
   end if

end subroutine getRealArray


!> Get an integer variable from the data structure
subroutine getInteger(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_getInteger"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to get
   integer, intent(out) :: val

   if (self%valueType == valueType%iVal) then
      val = self%iVal
   else
      call env%error("Data structure does not contain integer variable", source)
   end if

end subroutine getInteger


!> Get an integer array from the data structure
subroutine getIntegerArray(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_getIntegerArray"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to get
   integer, intent(out) :: val(:)

   if (self%valueType == valueType%iArr) then
      if (size(val) >= size(self%iArr)) then
         val(:) = self%iArr
      else
         call env%error("Insufficient size for integer array", source)
      end if
   else
      call env%error("Data structure does not contain integer array", source)
   end if

end subroutine getIntegerArray


!> Get a logical variable from the data structure
subroutine getLogical(self, env, val)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_getLogical"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Value to get
   logical, intent(out) :: val

   if (self%valueType == valueType%lVal) then
      val = self%lVal
   else
      call env%error("Data structure does not contain logical variable", source)
   end if

end subroutine getLogical


!> Get reference to a dictionary from the data structure
subroutine getDict(self, env, ptr)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_getDict"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Reference to get dictionary
   class(TDict), pointer, intent(out) :: ptr

   if (self%valueType == valueType%dict) then
      ptr => self%dict
   else
      call env%error("Data structure does not contain dictionary", source)
   end if

end subroutine getDict


!> Get reference to a list from the data structure
subroutine getList(self, env, ptr)

   !> Source for error creation
   character(len=*), parameter :: source = "type_abstractdata_getList"

   !> Instance of data structure
   class(TValue), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Reference to get dictionary
   class(TList), pointer, intent(out) :: ptr

   if (self%valueType == valueType%list) then
      ptr => self%list
   else
      call env%error("Data structure does not contain list", source)
   end if

end subroutine getList


!> Accept a visitor to transverse the data structure
recursive subroutine accept(self, visitor)

   !> Instance of the data structure
   class(TValue), intent(in) :: self

   !> Visitor for this value.
   class(TVisitor), intent(inout) :: visitor

   select case(self%valueType)
   case(valueType%list); call self%list%accept(visitor)
   case(valueType%dict); call self%dict%accept(visitor)
   case(valueType%cVal); call visitor%visit(self%cVal)
   case(valueType%rVal); call visitor%visit(self%rVal)
   case(valueType%rArr); call visitor%visit(self%rArr)
   case(valueType%iVal); call visitor%visit(self%iVal)
   case(valueType%iArr); call visitor%visit(self%iArr)
   case(valueType%lVal); call visitor%visit(self%lVal)
   case default; call visitor%visit(self)
   end select

end subroutine accept


!> Accept a visitor to transverse the hash map
recursive subroutine acceptDict(self, visitor)

   !> Instance of the data structure
   class(TDict), intent(in) :: self

   !> Visitor for this value.
   class(TVisitor), intent(inout) :: visitor

   call visitor%visit(self)

end subroutine acceptDict


!> Accept a visitor to transverse the list
recursive subroutine acceptList(self, visitor)

   !> Instance of the data structure
   class(TList), intent(in) :: self

   !> Visitor for this value.
   class(TVisitor), intent(inout) :: visitor

   call visitor%visit(self)

end subroutine acceptList


!> Return information on value type
function getType(self) result(vtype)

   !> Instance of the data structure
   class(TValue), intent(in) :: self

   !> Value type
   integer :: vtype

   vtype = self%valueType

end function getType


!> Compare to values
function equalsValue(lhs, rhs) result(equals)

   !> Left hand side value
   class(TValue), intent(in) :: lhs

   !> Right hand side value
   class(TValue), intent(in) :: rhs

   !> Values are equal
   logical :: equals

   equals = lhs%valueType == rhs%valueType
   if (equals) then
      select case(lhs%valueType)
      case(valueType%list); equals = associated(lhs%list, rhs%list)
      case(valueType%dict); equals = associated(lhs%dict, rhs%dict)
      case(valueType%cVal); equals = lhs%cVal == rhs%cVal
      case(valueType%rVal); equals = lhs%rVal == rhs%rVal
      case(valueType%rArr); equals = size(lhs%rArr) == size(rhs%rArr)
         if (equals) equals = all(lhs%rArr == rhs%rArr)
      case(valueType%iVal); equals = lhs%iVal == rhs%iVal
      case(valueType%iArr); equals = size(lhs%iArr) == size(rhs%iArr)
         if (equals) equals = all(lhs%iArr == rhs%iArr)
      case(valueType%lVal); equals = lhs%lVal .eqv. rhs%lVal
      end select
   end if

end function equalsValue


end module xtb_type_abstractdata
