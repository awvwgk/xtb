! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

!! ========================================================================
!  DRIVER FOR SINGLEPOINT CALCULATIONS
!  => call on (mini)scf to perform GFN-xTB single point calculations
!  => call on qmdff to run the QMDFF on a solvent file
!  => call on external driver library, interfaces
!     -> Turbomole, ORCA, driver (which is by itself an interface)
!  use the singlepoint function to get the appropiate call on the necessary
!  functions
!! ========================================================================
module xtb_single
   use xtb_mctc_accuracy, only : wp
   implicit none

   character(len=*),private,parameter :: outfmt = &
      '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'

contains

subroutine singlepoint&
      & (env, mol, wfn, calc, egap, et, maxiter, prlevel, restart, lgrad, &
      &  acc, etot, g, sigma, res, forceUpdate)
   use xtb_mctc_convert

!! ========================================================================
!  type definitions
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_calculator
   use xtb_type_data
   use xtb_type_pcem
   use xtb_type_latticepoint

!! ========================================================================
   use xtb_aoparam
   use xtb_setparam
   use xtb_fixparam
   use xtb_scanparam
   use xtb_sphereparam

!! ========================================================================
   use xtb_solv_gbobc, only : lgbsa
   use xtb_scf, only : scf
   use xtb_qmdff, only : ff_eg,ff_nonb,ff_hb
   use xtb_extern_mopac, only : runMopac
   use xtb_extern_orca, only : runOrca
   use xtb_peeq, only : peeq
   use xtb_embedding, only : read_pcem
   implicit none

   character(len=*), parameter :: source = 'single'

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(inout) :: mol
   type(TWavefunction),intent(inout) :: wfn
   type(TCalculator), intent(inout) :: calc
   type(tb_pcem) :: pcem
   real(wp),intent(inout) :: egap
   real(wp),intent(in)    :: et
   integer, intent(in)    :: maxiter
   integer, intent(in)    :: prlevel
   logical, intent(in)    :: restart
   logical, intent(in)    :: lgrad
   real(wp),intent(in)    :: acc
   logical, intent(in), optional :: forceUpdate
   real(wp),intent(out)   :: etot
   real(wp),intent(out)   :: g(3,mol%n)
   type(scc_results),intent(out) :: res
   real(wp),intent(out)   :: sigma(3,3)
   real(wp), allocatable :: latticePoint(:, :)
   integer  :: i,ich
   integer  :: mode_sp_run = 1
   real(wp) :: efix
   logical  :: inmol
   logical, parameter :: ccm = .true.
   logical :: exitRun, update
!  real(wp) :: efix1,efix2
!  real(wp),dimension(3,n) :: gfix1,gfix2

   call mol%update
   if (mol%npbc > 0) call generate_wsc(mol,mol%wsc,[1,1,1])

   call calc%latp%update(env, mol%lattice, update)
   call env%check(exitRun)
   if (exitRun) then
      call env%error("Generation of lattice points failed", source)
      return
   end if

   if (present(forceUpdate)) then
      update = update .or. forceUpdate
   end if
   call calc%latp%getLatticePoints(latticePoint, 40.0_wp)

   if (update) then
      call calc%neighList%generate(env, mol%xyz, 40.0_wp, latticePoint, .false.)
      call calc%wsCell%generate(env, mol%xyz, 40.0_wp, latticePoint, .false.)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Generation of neighbour lists failed", source)
         return
      end if
   else
      call calc%neighList%update(mol%xyz, latticePoint)
      call calc%wsCell%update(mol%xyz, latticePoint)
   end if

   etot = 0.0_wp
   efix = 0.0_wp
   g = 0.0_wp
   sigma = 0.0_wp

!! ========================================================================
!  external constrains which can be applied beforehand

   ! check for external point charge field
   call open_file(ich,pcem_file,'r')
   if (ich.ne.-1) then
      call read_pcem(ich,pcem)
      call close_file(ich)
   endif

!! ========================================================================
!  actual calculation
   select case(mode_extrun)
   case default
      call scf &
         & (env, mol, wfn, calc%basis, calc%param, pcem, calc%latp, &
         &  calc%neighList, calc%wsCell, egap, et, maxiter, prlevel, restart, &
         &  lgrad, acc, etot, g, res)

   case(p_ext_eht)
      call peeq &
         & (env, mol, wfn, calc%basis, calc%param, calc%latp, calc%neighList, &
         &  calc%wsCell, egap, et, prlevel, lgrad, ccm, acc, etot, g, sigma, res)

   case(p_ext_qmdff)
      call ff_eg  (mol%n,mol%at,mol%xyz,etot,g)
      call ff_nonb(mol%n,mol%at,mol%xyz,etot,g)
      call ff_hb  (mol%n,mol%at,mol%xyz,etot,g)

   case(p_ext_orca)
      call runOrca(env,mol,etot,g)

   case(p_ext_turbomole)
      call external_turbomole(mol%n,mol%at,mol%xyz,wfn%nel,wfn%nopen, &
         &                    lgrad,etot,g,res%dipole,lgbsa)

   case(p_ext_mopac)
      call runMopac(env,mol%n,mol%at,mol%xyz,etot,g)

   end select

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Electronic structure method terminated", source)
      return
   end if

!! ========================================================================
!  post processing of gradient and energy

   ! point charge embedding gradient file
   if (pcem%n > 0) then
      call open_file(ich,pcem_grad,'w')
      do i=1,pcem%n
         write(ich,'(3f12.8)')pcem%grd(1:3,i)
      enddo
      call close_file(ich)
   endif

!! ------------------------------------------------------------------------
!  various external potentials
   call constrain_pot(potset,mol%n,mol%at,mol%xyz,g,efix)
   call constrpot   (mol%n,mol%at,mol%xyz,g,efix)
   call cavity_egrad(mol%n,mol%at,mol%xyz,efix,g)
   call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,g)

!! ------------------------------------------------------------------------
!  fixing of certain atoms
!  print*,abs(efix/etot)
   etot = etot + efix
   res%e_total = etot
   res%gnorm = norm2(g)
   if (fixset%n.gt.0) then
      do i=1, fixset%n
         !print*,i,fixset%atoms(i)
         g(1:3,fixset%atoms(i))=0
      enddo
   endif

   if (prlevel.ge.2) then
      ! start with summary header
      if (.not.silent) then
         write(env%unit,'(9x,53(":"))')
         write(env%unit,'(9x,"::",21x,a,21x,"::")') "SUMMARY"
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,outfmt) "total energy      ", res%e_total,"Eh   "
      if (.not.silent.and.lgbsa) then
         write(env%unit,outfmt) "total w/o Gsasa/hb", &
            &  res%e_total-res%g_sasa-res%g_hb-res%g_shift, "Eh   "
      endif
      write(env%unit,outfmt) "gradient norm     ", res%gnorm,  "Eh/a0"
      write(env%unit,outfmt) "HOMO-LUMO gap     ", res%hl_gap, "eV   "
      if (.not.silent) then
         if (verbose) then
            write(env%unit,'(9x,"::",49("."),"::")')
            write(env%unit,outfmt) "HOMO orbital eigv.", wfn%emo(wfn%ihomo),  "eV   "
            write(env%unit,outfmt) "LUMO orbital eigv.", wfn%emo(wfn%ihomo+1),"eV   "
         endif
         write(env%unit,'(9x,"::",49("."),"::")')
         if (gfn_method.eq.2) call print_gfn2_results(env%unit,res,verbose,lgbsa)
         if (gfn_method.eq.1) call print_gfn1_results(env%unit,res,verbose,lgbsa)
         if (gfn_method.eq.0) call print_gfn0_results(env%unit,res,verbose,lgbsa)
         write(env%unit,outfmt) "add. restraining  ", efix,       "Eh   "
         if (verbose) then
            write(env%unit,'(9x,"::",49("."),"::")')
            write(env%unit,outfmt) "atomisation energy", res%e_atom, "Eh   "
         endif
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,'(a)')
   endif

end subroutine singlepoint

subroutine print_gfn0_results(iunit,res,verbose,lgbsa)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lgbsa
   write(iunit,outfmt) "H0 energy         ", res%e_elec, "Eh   "
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
   write(iunit,outfmt) "electrostat energy", res%e_es,   "Eh   "
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   !write(iunit,outfmt) "   -> Gborn       ", res%g_born, "Eh   " ! not saved
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   !write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   " ! not saved
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   write(iunit,outfmt) "dispersion energy ", res%e_disp, "Eh   "
   write(iunit,outfmt) "short-range corr. ", res%e_xb,   "Eh   "
end subroutine print_gfn0_results

subroutine print_gfn1_results(iunit,res,verbose,lgbsa)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lgbsa
   write(iunit,outfmt) "SCC energy        ", res%e_elec, "Eh   "
   write(iunit,outfmt) "-> electrostatic  ", res%e_es,   "Eh   "
   if (lgbsa) then
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   write(iunit,outfmt) "   -> Gborn       ", res%g_born, "Eh   "
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   "
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   endif
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
   write(iunit,outfmt) "dispersion energy ", res%e_disp, "Eh   "
   write(iunit,outfmt) "halogen bond corr.", res%e_xb,   "Eh   "
end subroutine print_gfn1_results

subroutine print_gfn2_results(iunit,res,verbose,lgbsa)
   use xtb_type_data
   integer, intent(in) :: iunit ! file handle (usually output_unit=6)
   type(scc_results),    intent(in) :: res
   logical,intent(in) :: verbose,lgbsa
   write(iunit,outfmt) "SCC energy        ", res%e_elec, "Eh   "
   write(iunit,outfmt) "-> isotropic ES   ", res%e_es,   "Eh   "
   write(iunit,outfmt) "-> anisotropic ES ", res%e_aes,  "Eh   "
   write(iunit,outfmt) "-> anisotropic XC ", res%e_axc,  "Eh   "
   write(iunit,outfmt) "-> dispersion     ", res%e_disp, "Eh   "
   if (lgbsa) then
   write(iunit,outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
   write(iunit,outfmt) "   -> Gborn       ", res%g_born, "Eh   "
   write(iunit,outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
   write(iunit,outfmt) "   -> Ghb         ", res%g_hb,   "Eh   "
   write(iunit,outfmt) "   -> Gshift      ", res%g_shift,"Eh   "
   endif
   write(iunit,outfmt) "repulsion energy  ", res%e_rep,  "Eh   "
end subroutine print_gfn2_results


end module xtb_single
