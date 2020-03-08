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

!> Linked list implementation
module xtb_type_linkedlist
   use xtb_mctc_accuracy, only : wp
   use xtb_type_abstractdata, only : TList, TValue
   implicit none
   private

   public :: TLinkedList


   !> Node of the linked list
   type :: TNode

      !> Value stored in the list
      type(TValue), allocatable :: val

      !> Reference to the next node
      type(TNode), pointer :: next

   end type TNode


   !> Factory to manipulate linked values
   type, extends(TList) :: TLinkedList

      !> Length of the list
      integer :: length

      !> Start of the linked list
      type(TNode), pointer :: head

      !> Last element of the linked list
      type(TNode), pointer :: last

      !> Cached index from last lookup
      integer :: cacheInd

      !> Cached pointer from last lookup
      type(TNode), pointer :: cachePtr

   contains

      !> Return length of list
      procedure :: getLength

      !> Get pointer value at this index
      procedure :: get

      !> Find a index of a value
      procedure :: find

      !> Create a new value and return pointer to this entry
      procedure :: append

      !> Clear the complete list
      procedure :: clear

   end type TLinkedList


   !> Initialize instance of linked list
   interface init
      module procedure :: initLinkedList
   end interface init


contains


subroutine initLinkedList(self)

   !> Instance of the linked list
   type(TLinkedList), intent(out) :: self

   self%length = 0
   self%head => null()
   self%last => null()

   self%cacheInd = 0
   self%cachePtr => null()

end subroutine initLinkedList


!> Clear the complete list
subroutine clear(self)

   !> Instance of the linked list
   class(TLinkedList), intent(inout) :: self

   type(TNode), pointer :: ptr, next

   ptr => self%head
   do while(associated(ptr))
      next => ptr%next
      call ptr%val%clear
      deallocate(ptr%val)
      deallocate(ptr)
      ptr => next
   end do

end subroutine clear


!> Create a new value and return pointer to this entry
subroutine append(self, ptr)

   !> Instance of the linked list
   class(TLinkedList), intent(inout) :: self

   type(TValue), pointer, intent(out) :: ptr

   if (associated(self%last)) then
      allocate(self%last%next)
      self%last => self%last%next
   else
      allocate(self%head)
      self%last => self%head
   end if

   self%length = self%length + 1
   nullify(self%last%next)
   allocate(self%last%val)

   ptr => self%last%val

end subroutine append


pure function getLength(self) result(length)

   !> Instance of the linked list
   class(TLinkedList), intent(in) :: self

   integer :: length

   length = self%length

end function getLength


!> Find a index of a value
function find(self, val) result(ind)

   !> Instance of the linked list
   class(TLinkedList), intent(inout) :: self

   type(TValue), intent(in) :: val

   integer :: ind

   type(TNode), pointer :: ptr

   integer :: ii

   ptr => self%head
   ii = 1
   do while(associated(ptr))
      if (ptr%val == val) then
         exit
      end if
      ptr => ptr%next
      ii = ii + 1
   end do

   if (associated(ptr)) then
      ind = ii
      self%cacheInd = ii
      self%cachePtr => ptr
   else
      ind = 0
   end if

end function find


!> Get pointer value at this index
function get(self, ind) result(val)

   !> Instance of the linked list
   class(TLinkedList), intent(inout) :: self

   !> Index to lookup
   integer, intent(in) :: ind

   !> Pointer to the content of the node
   type(TValue), pointer :: val

   type(TNode), pointer :: ptr

   integer :: ii, iStart

   if (self%cacheInd == ind) then
      val => self%cachePtr%val
      return
   end if

   if (self%cacheInd > 0 .and. self%cacheInd < ind) then
      iStart = self%cacheInd
      ptr => self%cachePtr
   else
      iStart = 1
      ptr => self%head
   end if

   do ii = iStart + 1, ind
      ptr => ptr%next
   end do
   self%cachePtr => ptr
   self%cacheInd = ind
   val => ptr%val

end function get


end module xtb_type_linkedlist
