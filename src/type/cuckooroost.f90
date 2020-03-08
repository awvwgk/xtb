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

!> Hash map implementation
module xtb_type_cuckooroost
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_hash, only : hash
   use xtb_mctc_io, only : stderr
   use xtb_type_abstractdata, only : TValue, TDict, TList, len
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: TCuckooSkein, len, init


   !> Wrapper for the key-value type
   type :: TCuckooNest

      !> Hashed key of the nest
      character(len=:), allocatable :: key

      type(TValue), allocatable :: val

   end type TCuckooNest


   !> Cuckoo hash map
   type, extends(TDict) :: TCuckooSkein

      !> Capacity of the flog
      integer :: capacity

      !> Seed for the hash functions
      integer :: seed(3)

      !> Maximum tries for claiming a nest before rehashing the flog
      integer :: maxTry

      !> Allocated values in the dictionary
      integer :: nAlloc

      !> Possible nests for cuckoos
      type(TCuckooNest), allocatable :: nest(:, :)

   contains

      !> Get length of the flog
      procedure :: getLength

      !> Insert a new key into the flog
      procedure :: insert

      !> Insert a new key into the flog
      procedure, private :: insertNest

      !> Get an existing key
      procedure :: get

      !> Get an existing nest
      procedure, private :: getNest

      !> Drop a key from the flog
      procedure :: drop

      !> Export key map as linked list
      procedure :: export

      !> Clear the flog
      procedure :: clear

   end type TCuckooSkein


   !> Initialize the cuckoo hash map
   interface init
      module procedure :: initSkein
   end interface init


   !> Overload operator for cuckoo nest
   interface operator(==)
      module procedure :: equals
   end interface


contains


!> Compare a cuckoo nest against a key
elemental function equals(lhs, rhs)

   !> Cuckoo nest
   type(TCuckooNest), intent(in) :: lhs

   !> Character key
   character(len=*), intent(in) :: rhs

   !> Cuckoo key matches key
   logical :: equals

   if (allocated(lhs%key)) then
      equals = lhs%key == rhs
   else
      equals = .false.
   end if

end function equals


!> Initialize a new flog
subroutine initSkein(self, capacity, seed)

   !> Instance of the flog
   type(TCuckooSkein), intent(out) :: self

   !> Capacity of this flog
   integer, intent(in) :: capacity

   !> Seeds for the hashing function, all most be unique
   integer, intent(in), optional :: seed(3)

   self%nAlloc = 0
   self%capacity = capacity/3+1
   self%maxTry = ceiling(3*log(real(3*self%capacity)))
   allocate(self%nest(1:self%capacity, 3))

   self%seed(:) = [2946901, 5039, 433494437]
   if (present(seed)) then
      if (seed(1)/=seed(2).and.seed(2)/=seed(3).and.seed(1)/=seed(3)) then
         self%seed(:) = seed
      end if
   end if

end subroutine initSkein


!> Get current number of occupied nest in the flog
pure function getLength(self) result(length)

   !> Instance of the flog
   class(TCuckooSkein), intent(in) :: self

   !> Number of allocated values
   integer :: length

   length = self%nAlloc

end function getLength


!> First hash function of the flog
pure function hash1(self, key) result(pos)

   !> Instance of the flog
   class(TCuckooSkein), intent(in) :: self

   !> Character key to be hashed
   character(len=*), intent(in) :: key

   !> Position in the flog
   integer :: pos

   pos = modulo(hash(key, self%seed(1)), self%capacity)+1

end function hash1


!> Second hash function of the flog
pure function hash2(self, key) result(pos)

   !> Instance of the flog
   class(TCuckooSkein), intent(in) :: self

   !> Character key to be hashed
   character(len=*), intent(in) :: key

   !> Position in the flog
   integer :: pos

   pos = modulo(hash(key, self%seed(2)), self%capacity)+1

end function hash2


!> Third hash function of the flog
pure function hash3(self, key) result(pos)

   !> Instance of the flog
   class(TCuckooSkein), intent(in) :: self

   !> Character key to be hashed
   character(len=*), intent(in) :: key

   !> Position in the flog
   integer :: pos

   pos = modulo(hash(key, self%seed(3)), self%capacity)+1

end function hash3


!> Get a pointer to the nest with the corresponding character key
subroutine getNest(self, key, ptr)

   !> Instance of the flog
   class(TCuckooSkein), target, intent(in) :: self

   !> Character key
   character(len=*), intent(in) :: key

   !> Pointer to the occupied nest
   type(TCuckooNest), pointer, intent(out) :: ptr

   integer :: ii

   ii = hash1(self, key)
   if (self%nest(ii, 1) == key) then
      ptr => self%nest(ii, 1)
   else
      ii = hash2(self, key)
      if (self%nest(ii, 2) == key) then
         ptr => self%nest(ii, 2)
      else
         ii = hash3(self, key)
         if (self%nest(ii, 3) == key) then
            ptr => self%nest(ii, 3)
         else
            ptr => null()
         end if
      end if
   end if

end subroutine getNest


!> Get a pointer to the nest with the corresponding character key
subroutine get(self, key, ptr)

   !> Source of the error
   character(len=*), parameter :: source = "type_cuckooroost_get"

   !> Instance of the flog
   class(TCuckooSkein), target, intent(inout) :: self

   !> Character key
   character(len=*), intent(in) :: key

   !> Pointer to the occupied nest
   type(TValue), pointer, intent(out) :: ptr

   type(TCuckooNest), pointer :: nest

   integer :: ii

   call self%getNest(key, nest)
   if (associated(nest)) then
      if (.not.allocated(nest%val)) then
         write(stderr, '("[FATAL]",1x,a,":",1x,a)') source, &
            & "Key '"//key//"' present, but value not allocated"
         ptr => null()
      end if
      ptr => nest%val
   else
      ptr => null()
   end if

end subroutine get


!> Insert a new key into the flog
recursive subroutine insert(self, key, ptr)

   !> Source of the error
   character(len=*), parameter :: source = "type_cuckooroost_get"

   !> Instance of the flog
   class(TCuckooSkein), target, intent(inout) :: self

   !> Character key to insert
   character(len=*), intent(in) :: key

   !> Pointer for the newly occupied nest
   type(TValue), pointer, intent(out) :: ptr

   !> Pointer for the newly occupied nest
   type(TCuckooNest), pointer :: nest

   call self%insertNest(key, nest)
   if (allocated(nest%val)) then
      ptr => nest%val
   else
      if (allocated(nest%val)) then
         ptr => nest%val
      else
         write(stderr, '("[FATAL]",1x,a,":",1x,a)') source, &
            & "Key '"//key//"' present, but value not allocated"
         ptr => null()
      end if
   end if

end subroutine insert


!> Insert a new key into the flog
recursive subroutine insertNest(self, key, ptr)

   !> Instance of the flog
   class(TCuckooSkein), intent(inout) :: self

   !> Character key to insert
   character(len=*), intent(in) :: key

   !> Pointer for the newly occupied nest
   type(TCuckooNest), pointer, intent(out) :: ptr

   !> Pointer for the newly occupied nest
   type(TCuckooNest), pointer :: tmp

   !> Temporay cuckoo searching for a new nest
   type(TCuckooNest) :: cuckoo

   integer :: ii, iTry

   call self%getNest(key, ptr)
   if (.not.associated(ptr)) then
      cuckoo%key = key
      allocate(cuckoo%val)
      do iTry = 1, self%maxTry
         ii = hash1(self, cuckoo%key)
         if (allocated(self%nest(ii, 1)%key)) then
            call swap(self%nest(ii, 1), cuckoo)
         else
            call move(cuckoo, self%nest(ii, 1))
            self%nAlloc = self%nAlloc + 1
            exit
         end if
         ii = hash2(self, cuckoo%key)
         if (allocated(self%nest(ii, 2)%key)) then
            call swap(self%nest(ii, 2), cuckoo)
         else
            call move(cuckoo, self%nest(ii, 2))
            self%nAlloc = self%nAlloc + 1
            exit
         end if
         ii = hash3(self, cuckoo%key)
         if (allocated(self%nest(ii, 3)%key)) then
            call swap(self%nest(ii, 3), cuckoo)
         else
            call move(cuckoo, self%nest(ii, 3))
            self%nAlloc = self%nAlloc + 1
            exit
         end if
      end do
      if (allocated(cuckoo%key)) then
         call rehash(self)
         call self%insertNest(cuckoo%key, tmp)
      end if
      call self%getNest(key, ptr)
   end if

end subroutine insertNest


!> Rehash the flog
recursive subroutine rehash(self)

   !> Instance of the flog
   class(TCuckooSkein), intent(inout) :: self

   !> Old cuckoo nests
   type(TCuckooNest), allocatable :: nest(:, :)

   !> Cuckoo for reassigning the value pointer
   type(TCuckooNest), pointer :: cuckoo

   character(len=:), allocatable :: key

   integer :: ii, jj
   integer :: oldCapacity
   integer :: nAlloc

   call move_alloc(self%nest, nest)

   oldCapacity = self%capacity
   self%capacity = oldCapacity + oldCapacity/2 + 1
   allocate(self%nest(1:self%capacity, 3))
   self%maxTry = ceiling(3*log(real(2*self%capacity)))

   do ii = 1, 3
      do jj = 1, oldCapacity
         if (allocated(nest(jj, ii)%key)) then
            self%nAlloc = self%nAlloc - 1
            call move_alloc(nest(jj, ii)%key, key)
            call self%insertNest(key, cuckoo)
            call move_alloc(nest(jj, ii)%val, cuckoo%val)
         end if
      end do
   end do

   deallocate(nest)

end subroutine rehash


!> Clear the flog
subroutine clear(self)

   !> Source of the error
   character(len=*), parameter :: source = "type_cuckooroost_clear"

   !> Instance of the flog
   class(TCuckooSkein), intent(inout) :: self

   integer :: ii, jj

   do ii = 1, 3
      do jj = 1, self%capacity
         if (allocated(self%nest(jj, ii)%key)) then
            call self%drop(self%nest(jj, ii)%key)
         end if
      end do
   end do

   if (self%nAlloc /= 0) then
      write(stderr, '("[FATAL]",1x,a,":",1x,a,1x,"(",i0,"/",i0,")")') source, &
         & "Memory leak in cuckoo hash map implementation", self%nAlloc, len(self)
   end if

end subroutine clear


!> Swap two cuckoos
elemental subroutine swap(lhs, rhs)

   !> Cuckoo currently occupying a nest
   type(TCuckooNest), intent(inout) :: lhs

   !> New cuckoo resident
   type(TCuckooNest), intent(inout) :: rhs

   type(TCuckooNest) :: tmp

   call move(lhs, tmp)
   call move(rhs, lhs)
   call move(tmp, rhs)

end subroutine swap


!> Move content from one nest to another
elemental subroutine move(from, to)

   !> New cuckoo resident
   type(TCuckooNest), intent(inout) :: from

   !> Empty cuckoo nest
   type(TCuckooNest), intent(inout) :: to

   call move_alloc(from%key, to%key)
   call move_alloc(from%val, to%val)

end subroutine move


!> Drop a cuckoo from the flog
subroutine drop(self, key)

   !> Instance of the flog
   class(TCuckooSkein), intent(inout) :: self

   !> Character key
   character(len=*), intent(in) :: key

   !> Pointer to the nest to be cleared
   type(TCuckooNest), pointer :: nest

   call self%getNest(key, nest)

   if (associated(nest)) then
      self%nAlloc = self%nAlloc - 1
      if (allocated(nest%val)) then
         call nest%val%clear
         deallocate(nest%val)
      end if
      if (allocated(nest%key)) then
         deallocate(nest%key)
      end if
   end if

end subroutine drop


!> Export keys as list
subroutine export(self, env, list)

   !> Instance of the flog
   class(TCuckooSkein), intent(inout) :: self

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Instance of the key list
   class(TList), intent(inout) :: list

   type(TValue), pointer :: key

   integer :: ii, jj

   do ii = 1, 3
      do jj = 1, self%capacity
         if (allocated(self%nest(jj, ii)%key)) then
            call list%append(key)
            call key%set(env, self%nest(jj, ii)%key)
         end if
      end do
   end do

end subroutine export


end module xtb_type_cuckooroost
