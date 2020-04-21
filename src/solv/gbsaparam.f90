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
module xtb_solv_gbsaparam
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : fourpi
   use xtb_mctc_convert, only : kcaltoau, aatoau
   use xtb_mctc_systools, only : getline
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: TGBSAData, readGBSAParameter


   !> Parametrisation data for the GBSA model
   type :: TGBSAData

      !> Dielectric constant of the solvent
      real(wp) :: dielectricConst

      !> Molecular mass of the solvent
      real(wp) :: molecularMass

      !> Density of the solvent
      real(wp) :: density

      !> Born radius scaling parameter
      real(wp) :: bornScale

      !> Shift for free energy calculations
      real(wp) :: freeEnergyShift

      !> Offset parameter for Born radii
      real(wp) :: bornOffset

      !> Probe radius for the numerical integration
      real(wp) :: probeRad

      !> Dielectric descreening parameter for each element
      real(wp), allocatable :: descreening(:)

      !> Surface tension scaling for each element
      real(wp), allocatable :: surfaceTension(:)

      !> Hydrogen bond strength for each element
      real(wp), allocatable :: hBondStrength(:)

   end type TGBSAData


   interface readGBSAParameter
      module procedure :: readGBSAParameterFile
      module procedure :: readGBSAParameterUnit
   end interface readGBSAParameter


contains


!> Read GBSA parametrisation data from parameter file
subroutine readGBSAParameterFile(file, env, input)

   !> Source for the error generation
   character(len=*), parameter :: source = 'solv_gbsaparam_readGBSAParameterFile'

   !> Parameter file name
   character(len=*), intent(in) :: file

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Parametrisation data
   type(TGBSAData), intent(out) :: input

   integer :: unit

   call open_file(unit, file, 'r')
   if (unit == -1) then
      call env%error("Could not open parameter file '"//file//"'", source)
      return
   end if

   call readGBSAParameter(unit, env, input)

   call close_file(unit)

end subroutine readGBSAParameterFile


!> Read GBSA parametrisation data from unit
subroutine readGBSAParameterUnit(unit, env, input)

   !> Source for the error generation
   character(len=*), parameter :: source = 'solv_gbsaparam_readGBSAParameterUnit'

   !> Unit bound to parameter file
   integer, intent(in) :: unit

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Parametrisation data
   type(TGBSAData), intent(out) :: input

   integer :: ii, error
   real(wp) :: param(8), hBondStrength, surfaceTension
   character(len=:), allocatable :: line
   character(len=1) :: tmp1
   character(len=3) :: tmp2

   allocate(input%descreening(94))
   allocate(input%surfaceTension(94))
   allocate(input%hBondStrength(94))

   do ii = 1, 8
      call getline(unit, line, error)
      if (error /= 0) then
         write(tmp1, '(i1)') ii
         call env%error("Unexpected end of data in line "//tmp1, source)
         exit
      end if
      read(line, *, iostat=error) param(8)
      if (error /= 0) then
         write(tmp1, '(i1)') ii
         call env%error("Could not interpret data in line "//tmp1//" as real", &
            & source)
         exit
      end if
   end do
   if (error /= 0) then
      return
   end if

   input%dielectricConst = param(1)
   input%molecularMass = param(2)
   input%density = param(3)
   input%bornScale = param(4)
   input%probeRad = param(5) * aatoau
   input%freeEnergyShift = param(6) * kcaltoau
   input%bornOffset = param(7) * aatoau * 0.1_wp

   do ii = 1, 94
      call getline(unit, line, error)
      if (error /= 0) then
         write(tmp2, '(i3)') ii + 8
         call env%error("Unexpected end of data in line "//trim(tmp2), source)
         exit
      end if
      read(line, *, iostat=error) surfaceTension, input%descreening(ii), &
         & hBondStrength
      if (error /= 0) then
         write(tmp2, '(i3)') ii + 8
         call env%error("Could not interpret data in line "//trim(tmp2)// &
            & " as reals", source)
         exit
      end if
      input%surfaceTension(ii) = fourpi * surfaceTension * 1.0e-5_wp
      input%hBondStrength(ii) = -hBondStrength**2 * kcaltoau
   end do

end subroutine readGBSAParameterUnit


end module xtb_solv_gbsaparam
