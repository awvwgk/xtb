subroutine test_latticepoint_pbc3d
   use assertion
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_boundaryconditions, only : boundaryCondition
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_latticepoint, only : TLatticePoint, init
   implicit none

   real(wp), parameter :: lattice_SiO2(3,3) = reshape(&
      &[ 8.7413053236641_wp,  0.0000000000000_wp,  0.0000000000000_wp,  &
      &  0.0000000000000_wp,  8.7413053236641_wp,  0.0000000000000_wp,  &
      &  0.0000000000000_wp,  0.0000000000000_wp,  8.7413053236641_wp], &
      & shape(lattice_SiO2))

   real(wp), parameter :: lattice_CaF2(3,3) = reshape(&
      &[ 5.9598811567890_wp,  2.1071361905157_wp,  3.6496669404404_wp,  &
      &  0.0000000000000_wp,  6.3214085715472_wp,  3.6496669404404_wp,  &
      &  0.0000000000000_wp,  0.0000000000000_wp,  7.2993338808807_wp], &
      & shape(lattice_CaF2))

   real(wp), parameter :: lattice_ammonia(3,3) = reshape(&
      &[ 6.4411018522600_wp,  0.0492571261505_wp,  0.2192046129910_wp,  &
      &  0.0462076831749_wp,  6.6435057067500_wp,  0.1670513770770_wp,  &
      &  0.2262248220170_wp, -0.9573234940220_wp,  6.7608039126200_wp], &
      & shape(lattice_ammonia)) * aatoau

   type(TEnvironment) :: env
   type(TLatticePoint) :: latp
   logical :: fail
   real(wp), allocatable :: latticePoint(:, :)

   call init(env)

   call init(latp, env, lattice_SiO2, boundaryCondition%pbc3d, 40.0_wp)
   call env%check(fail)
   call assert(.not.fail)
   call assert(allocated(latp%trans))
   call assert(allocated(latp%dist2))
   call assert_eq(latp%nTrans, 739)

   call latp%getLatticePoints(latticePoint)
   call assert_eq(size(latticePoint, dim=2), 739)

   !> Simply switch lattice without reinitalization
   call latp%update(env, lattice_CaF2)
   call env%check(fail)
   call assert(.not.fail)
   call assert(allocated(latp%trans))
   call assert(allocated(latp%dist2))
   call assert_eq(latp%nTrans, 1553)

   call latp%getLatticePoints(latticePoint, 30.0_wp)
   call assert_eq(size(latticePoint, dim=2), 791)

   call latp%getLatticePoints(latticePoint, 50.0_wp)
   call assert_eq(size(latticePoint, dim=2), 1553)

   !> Reinitialize with new lattice and new cutoff
   call init(latp, env, lattice_ammonia, boundaryCondition%pbc3d, 60.0_wp)
   call env%check(fail)
   call assert(.not.fail)
   call assert(allocated(latp%trans))
   call assert(allocated(latp%dist2))
   call assert_eq(latp%nTrans, 833)

   call latp%getLatticePoints(latticePoint)
   call assert_eq(size(latticePoint, dim=2), 833)

   call latp%getLatticePoints(latticePoint, 40.0_wp)
   call assert_eq(size(latticePoint, dim=2), 317)

   !> Reinitialize generator to exclude inversion symmetry
   call init(latp, env, lattice_ammonia, boundaryCondition%pbc3d, 50.0_wp, &
      & excludeInversion=.true.)
   call env%check(fail)
   call assert(.not.fail)
   call assert(allocated(latp%trans))
   call assert(allocated(latp%dist2))
   call assert_eq(latp%nTrans, 264)

   call latp%getLatticePoints(latticePoint, 40.0_wp)
   call assert_eq(size(latticePoint, dim=2), 159) ! 317 / 2 + 1

   call latp%getLatticePoints(latticePoint, 33.0_wp)
   call assert_eq(size(latticePoint, dim=2), 105)

   call terminate(afail)

end subroutine test_latticepoint_pbc3d
