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

!> Implementation of an ONIOM calculator stack
module xtb_oniom_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_gbsa, only : TBorn
   use xtb_solv_model, only : info, newSolvationModel, newBornModel
   use xtb_type_atomlist, only : TAtomlist
   use xtb_type_basisset, only : TBasisset
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, len, init
   use xtb_type_param, only : scc_parameter
   use xtb_type_pcem
   use xtb_type_solvation, only : TSolvation
   use xtb_type_restart, only : TRestart
   use xtb_type_wsc, only : tb_wsc
   use xtb_xtb_data, only : TxTBData
   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_gfnff_calculator, only : TGFFCalculator
   use xtb_setparam
   use xtb_fixparam
   use xtb_scanparam
   use xtb_sphereparam
   use xtb_scf, only : scf
   use xtb_qmdff, only : ff_eg,ff_nonb,ff_hb
   use xtb_peeq, only : peeq
   use xtb_embedding, only : read_pcem
   use xtb_metadynamic
   use xtb_constrainpot
   implicit none
   private

   public :: TOniomCalculator


   !> Calculator stack to evaluate the ONIOM scheme for an inner:outer method
   !
   !  E(i:o) = lambda*(E(total@outer) - E(inner@outer)) + E(inner@inner)
   type, extends(TCalculator) :: TOniomCalculator

      !> Inner region, inner method
      class(TCalculator), allocatable :: inner

      !> Inner region, outer method
      class(TCalculator), allocatable :: subtr

      !> Outer region, outer method
      class(TCalculator), allocatable :: outer

      !> List of atoms for inner region
      type(TAtomlist) :: list

      !> Interaction scaling parameter for outer region
      real(wp), allocatable :: lambda

   contains

      !> Perform single point calculation
      procedure :: singlepoint

      !> Write informative printout
      procedure :: writeInfo

   end type TOniomCalculator


   character(len=*),private,parameter :: outfmt = &
      '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'


contains


!> Perform single point calculation
subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)

   !> Source of the generated errors
   character(len=*), parameter :: source = 'oniom_calculator_singlepoint'

   !> Calculator instance
   class(TOniomCalculator), intent(inout) :: self

   !> Computational environment
   type(TEnvironment), intent(inout) :: env

   !> Molecular structure data
   type(TMolecule), intent(inout) :: mol

   !> Wavefunction data
   type(TRestart), intent(inout) :: chk

   !> Print level for IO
   integer, intent(in) :: printlevel

   !> Restart from previous results
   logical, intent(in) :: restart

   !> Total energy
   real(wp), intent(out) :: energy

   !> Molecular gradient
   real(wp), intent(out) :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(out) :: sigma(:, :)

   !> HOMO-LUMO gap
   real(wp), intent(out) :: hlgap

   !> Detailed results
   type(scc_results), intent(out) :: results

   type(TMolecule) :: inner
   type(scc_results) :: resInner
   logical :: exitRun
   real(wp) :: energyInner, energyOuter, energySubtr, efix
   real(wp) :: gradDiff1, gradDiff2
   real(wp) :: sigmaInner(3, 3), sigmaOuter(3, 3), sigmaSubtr(3, 3)
   real(wp) :: gapInner, gapOuter, gapSubtr
   integer, allocatable :: at(:)
   real(wp), allocatable :: xyz(:, :), gradOuter(:, :), gradInner(:, :), charges(:)

   call mol%update
   if (mol%npbc > 0) call generate_wsc(mol,mol%wsc)

   if (.not.(allocated(self%outer) .and. allocated(self%subtr) &
      & .and. allocated(self%inner))) then
      call env%error("ONIOM calculation is not setup correctly", source)
      return
   end if

   allocate(gradOuter(3, len(mol)))
   energy = 0.0_wp
   gradOuter(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   hlgap = 0.0_wp
   efix = 0.0_wp

   if (allocated(self%outer)) then
      call self%outer%singlepoint(env, mol, chk, printlevel-1, restart, &
         & energyOuter, gradOuter, sigmaOuter, gapOuter, results)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Calculation of outer region failed", source)
         return
      end if

      energy = energy + energyOuter
      sigma = sigma + sigmaOuter
      select type(calc => self%outer)
      type is (TxTBCalculator)
         call self%list%gather(chk%wfn%q, charges)
      type is (TGFFCalculator)
         call self%list%gather(calc%topo%q, charges)
      end select
   end if

   ! Collect the inner region
   call self%list%gather(mol%at, at)
   call self%list%gather(mol%xyz, xyz)
   call init(inner, at, xyz, chrg=mol%chrg, uhf=mol%uhf)

   ! Allocate space to hold inner region results
   allocate(gradInner(3, len(inner)))
   gradInner(:, :) = 0.0_wp

   if (allocated(self%subtr)) then
      ! Advance pointer in restart list
      !chkp => chkp%next

      call self%subtr%singlepoint(env, inner, chk%next, printlevel-1, restart, &
         & energySubtr, gradInner, sigmaSubtr, gapSubtr, resInner)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Calculation of inner region failed (-)", source)
         return
      end if

      call self%list%scatter(gradInner, gradient, scale=+0.0_wp)
      energy = energy - energySubtr
      sigma = sigma - sigmaSubtr
      hlgap = hlgap - gapSubtr
      results%dipole(:) = results%dipole - resInner%dipole
   end if

   ! backup gradient of inner region
   xyz(:, :) = gradInner
   ! Reset inner region results
   gradInner(:, :) = 0.0_wp

   if (allocated(self%inner)) then
      ! Advance pointer in restart list
      !chkp => chkp%next

      call self%inner%singlepoint(env, inner, chk%next%next, printlevel-1, restart, &
         & energyInner, gradInner, sigmaInner, gapInner, resInner)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Calculation of inner region failed (+)", source)
         return
      end if

      call self%list%scatter(gradInner, gradient, scale=-1.0_wp)
      energy = energy + energyInner
      sigma = sigma + sigmaInner
      hlgap = hlgap + gapInner
      results%dipole(:) = results%dipole + resInner%dipole
   end if

   ! Difference between inner method on inner region vs outer method on inner region
   gradDiff1 = norm2(gradInner - xyz)
   ! Difference between outer method on outer region vs outer method on outer region
   call self%list%gather(gradOuter, gradInner)
   gradDiff2 = norm2(gradInner - xyz)

   ! ------------------------------------------------------------------------
   !  various external potentials
   call constrain_pot(potset,mol%n,mol%at,mol%xyz,gradient,efix)
   call constrpot   (mol%n,mol%at,mol%xyz,gradient,efix)
   call cavity_egrad(mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,gradient)

   energy = energy + efix
   gradient(:, :) = gradOuter + gradient
   results%gnorm = norm2(gradient)
   results%e_total = energy
   results%hl_gap = hlgap

   if (printlevel >= 2) then
      ! start with summary header
      write(env%unit,'(9x,53(":"))')
      write(env%unit,'(9x,"::",21x,a,21x,"::")') "SUMMARY"
      write(env%unit,'(9x,53(":"))')
      write(env%unit,outfmt) "total energy      ", results%e_total,"Eh   "
      write(env%unit,outfmt) "gradient norm     ", results%gnorm,  "Eh/a0"
      write(env%unit,outfmt) "HOMO-LUMO gap     ", results%hl_gap, "eV   "
      write(env%unit,'(9x,"::",49("."),"::")')
      write(env%unit,outfmt) "total energy (+)  ", energyInner,    "Eh   "
      write(env%unit,outfmt) "total energy (-)  ", energySubtr,    "Eh   "
      write(env%unit,outfmt) "total energy (o)  ", energyOuter,    "Eh   "
      write(env%unit,outfmt) "grad. diff. (+:-) ", gradDiff1,      "Eh/a0"
      write(env%unit,outfmt) "grad. diff. (o:-) ", gradDiff2,      "Eh/a0"
      if (allocated(charges)) then
         write(env%unit,outfmt) &
            & "total charge (+)  ", inner%chrg,   "e    "
         write(env%unit,outfmt) &
            & "total charge (o)  ", sum(charges), "e    "
      end if
      write(env%unit,'(9x,"::",49("."),"::")')
      write(env%unit,outfmt) "add. restraining  ", efix,       "Eh   "
      write(env%unit,'(9x,53(":"))')
      write(env%unit,'(a)')
   end if

end subroutine singlepoint


!> Write informative printout
subroutine writeInfo(self, unit, mol, printlevel)

   !> Calculator instance
   class(TOniomCalculator), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Print level for IO
   integer, intent(in) :: printlevel

   type(TMolecule) :: inner
   integer, allocatable :: at(:)
   real(wp), allocatable :: xyz(:, :)
   character(len=:), allocatable :: region

   ! Collect the inner region
   call self%list%gather(mol%at, at)
   call self%list%gather(mol%xyz, xyz)
   call init(inner, at, xyz, chrg=mol%chrg, uhf=mol%uhf)

   call self%list%to_string(region)

   write(unit, '(/,*(1x,a))') "***", "ONIOM CALCULATION", "***"
   write(unit, '(a,1x,a)') "inner region", region
   write(unit, '(/,*(1x,a))') "***", "CALCULATOR FOR OUTER REGION", "***"
   call self%outer%writeInfo(unit, mol, printlevel-1)
   write(unit, '(/,*(1x,a))') "***", "CALCULATOR FOR INNER REGION", "***"
   call self%inner%writeInfo(unit, inner, printlevel-1)

end subroutine writeInfo


end module xtb_oniom_calculator
