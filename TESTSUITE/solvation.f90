subroutine test_solv_bornrad
   use assertion
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_solv_born
   implicit none
   real(wp), parameter :: thr = 1.0e-10_wp
   real(wp), parameter :: thr2 = 1.0e-8_wp
   integer, parameter :: nat = 24
   integer, parameter :: at(nat) = [6,7,6,7,6,6,6,8,7,6,8,7,6,6, &
      &                             1,1,1,1,1,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[ 2.02799738646442_wp,  0.09231312124713_wp, -0.14310895950963_wp, &
      &  4.75011007621000_wp,  0.02373496014051_wp, -0.14324124033844_wp, &
      &  6.33434307654413_wp,  2.07098865582721_wp, -0.14235306905930_wp, &
      &  8.72860718071825_wp,  1.38002919517619_wp, -0.14265542523943_wp, &
      &  8.65318821103610_wp, -1.19324866489847_wp, -0.14231527453678_wp, &
      &  6.23857175648671_wp, -2.08353643730276_wp, -0.14218299370797_wp, &
      &  5.63266886875962_wp, -4.69950321056008_wp, -0.13940509630299_wp, &
      &  3.44931709749015_wp, -5.48092386085491_wp, -0.14318454855466_wp, &
      &  7.77508917214346_wp, -6.24427872938674_wp, -0.13107140408805_wp, &
      & 10.30229550927022_wp, -5.39739796609292_wp, -0.13672168520430_wp, &
      & 12.07410272485492_wp, -6.91573621641911_wp, -0.13666499342053_wp, &
      & 10.70038521493902_wp, -2.79078533715849_wp, -0.14148379504141_wp, &
      & 13.24597858727017_wp, -1.76969072232377_wp, -0.14218299370797_wp, &
      &  7.40891694074004_wp, -8.95905928176407_wp, -0.11636933482904_wp, &
      &  1.38702118184179_wp,  2.05575746325296_wp, -0.14178615122154_wp, &
      &  1.34622199478497_wp, -0.86356704498496_wp,  1.55590600570783_wp, &
      &  1.34624089204623_wp, -0.86133716815647_wp, -1.84340893849267_wp, &
      &  5.65596919189118_wp,  4.00172183859480_wp, -0.14131371969009_wp, &
      & 14.67430918222276_wp, -3.26230980007732_wp, -0.14344911021228_wp, &
      & 13.50897177220290_wp, -0.60815166181684_wp,  1.54898960808727_wp, &
      & 13.50780014200488_wp, -0.60614855212345_wp, -1.83214617078268_wp, &
      &  5.41408424778406_wp, -9.49239668625902_wp, -0.11022772492007_wp, &
      &  8.31919801555568_wp, -9.74947502841788_wp,  1.56539243085954_wp, &
      &  8.31511620712388_wp, -9.76854236502758_wp, -1.79108242206824_wp],&
      &  shape(xyz))
   real(wp), parameter :: vdwRad(4) = aatoau * [&
      & 1.45515_wp, 1.31125_wp, 1.24085_wp, 1.09155_wp]
   real(wp), parameter :: descreening1(4) = [&
      & 0.76509182_wp, 1.03142132_wp, 0.30000000_wp, 0.59799831_wp]
   real(wp), parameter :: descreening2(4) = [&
      & 0.81869146_wp, 1.16315604_wp, 0.30000000_wp, 0.30000000_wp]
   real(wp), parameter :: obc1(3) = [0.80_wp, 0.00_wp, 2.91_wp]
   real(wp), parameter :: obc2(3) = [1.00_wp, 0.80_wp, 4.85_wp]

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TBornIntegrator) :: gbobc
   type(TNeighbourlist) :: neighList

   integer :: ii, jj, kk
   logical :: fail
   real(wp) :: eps(3, 3)
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   real(wp), allocatable :: bornRad(:)
   real(wp), allocatable :: dbrdr(:, :, :)
   real(wp), allocatable :: dbrdL(:, :, :)
   real(wp), allocatable :: brl(:), brr(:), drdum(:, :, :), dLdum(:, :, :)

   call init(env)
   call init(mol, at, xyz)
   call init(neighList, len(mol))

   allocate(bornRad(nat))
   allocate(dbrdr(3, nat, nat))
   allocate(dbrdL(3, 3, nat))

   call neighList%generate(env, mol%xyz, 30.0_wp, trans, .false.)

   call init(gbobc, vdwRad, descreening1, 1.54870752_wp, 0.008071897_wp * aatoau, &
      & obc1, cutoff=30.0_wp)

   call gbobc%getBornRad(neighList%neighs, neighList, mol%id, bornRad, dbrdr, dbrdL)

   call assert_close(bornRad( 1), 4.8657451881819_wp, thr)
   call assert_close(bornRad( 3), 5.5380217411303_wp, thr)
   call assert_close(bornRad( 6), 7.0662782350137_wp, thr)
   call assert_close(bornRad(10), 5.7163033485686_wp, thr)
   call assert_close(bornRad(12), 5.6996090188203_wp, thr)
   call assert_close(bornRad(17), 3.8864765131601_wp, thr)
   call assert_close(bornRad(22), 3.8540074597995_wp, thr)

   call init(gbobc, vdwRad, descreening2, 1.65675010_wp, 0.0_wp, obc2, &
      & cutoff=30.0_wp)

   call gbobc%getBornRad(neighList%neighs, neighList, mol%id, bornRad, dbrdr, dbrdL)

   call assert_close(bornRad( 1), 5.5593224957980_wp, thr)
   call assert_close(bornRad( 3), 7.1245650398770_wp, thr)
   call assert_close(bornRad( 6), 12.497674302988_wp, thr)
   call assert_close(bornRad(10), 7.5659071458152_wp, thr)
   call assert_close(bornRad(12), 7.8618127813664_wp, thr)
   call assert_close(bornRad(17), 4.5297343895426_wp, thr)
   call assert_close(bornRad(22), 4.4662958502376_wp, thr)

   if (afail > 0) call terminate(afail)

   allocate(brr(nat))
   allocate(brl(nat))
   allocate(drdum(3, nat, nat))
   allocate(dLdum(3, 3, nat))

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbobc%getBornRad(neighList%neighs, neighList, mol%id, &
            & brr, drdum, dLdum)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbobc%getBornRad(neighList%neighs, neighList, mol%id, &
            & brl, drdum, dLdum)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step

         do kk = 1, nat
            call assert_close(dbrdr(jj, ii, kk), (brr(kk) - brl(kk))*step2, thr2)
         end do
      end do
   end do

   if (afail > 0) call terminate(afail)

   ! check numerical strain derivatives
   eps = unity
   do ii = 1, 3
      do jj = 1, ii
         eps(jj, ii) = eps(jj, ii) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbobc%getBornRad(neighList%neighs, neighList, mol%id, &
            & brr, drdum, dLdum)

         eps(jj, ii) = eps(jj, ii) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbobc%getBornRad(neighList%neighs, neighList, mol%id, &
            & brl, drdum, dLdum)

         eps(jj, ii) = eps(jj, ii) + step
         do kk = 1, nat
            call assert_close(dbrdL(jj, ii, kk), (brr(kk) - brl(kk))*step2, thr2)
         end do
      end do
   end do

   call terminate(afail)
end subroutine test_solv_bornrad

subroutine test_solv_sasaint
   use assertion
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_solv_sasa
   implicit none
   real(wp), parameter :: thr = 1.0e-10_wp
   real(wp), parameter :: thr2 = 1.0e-8_wp
   integer, parameter :: nat = 24
   integer, parameter :: at(nat) = [6,7,6,7,6,6,6,8,7,6,8,7,6,6, &
      &                             1,1,1,1,1,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[ 2.02799738646442_wp,  0.09231312124713_wp, -0.14310895950963_wp, &
      &  4.75011007621000_wp,  0.02373496014051_wp, -0.14324124033844_wp, &
      &  6.33434307654413_wp,  2.07098865582721_wp, -0.14235306905930_wp, &
      &  8.72860718071825_wp,  1.38002919517619_wp, -0.14265542523943_wp, &
      &  8.65318821103610_wp, -1.19324866489847_wp, -0.14231527453678_wp, &
      &  6.23857175648671_wp, -2.08353643730276_wp, -0.14218299370797_wp, &
      &  5.63266886875962_wp, -4.69950321056008_wp, -0.13940509630299_wp, &
      &  3.44931709749015_wp, -5.48092386085491_wp, -0.14318454855466_wp, &
      &  7.77508917214346_wp, -6.24427872938674_wp, -0.13107140408805_wp, &
      & 10.30229550927022_wp, -5.39739796609292_wp, -0.13672168520430_wp, &
      & 12.07410272485492_wp, -6.91573621641911_wp, -0.13666499342053_wp, &
      & 10.70038521493902_wp, -2.79078533715849_wp, -0.14148379504141_wp, &
      & 13.24597858727017_wp, -1.76969072232377_wp, -0.14218299370797_wp, &
      &  7.40891694074004_wp, -8.95905928176407_wp, -0.11636933482904_wp, &
      &  1.38702118184179_wp,  2.05575746325296_wp, -0.14178615122154_wp, &
      &  1.34622199478497_wp, -0.86356704498496_wp,  1.55590600570783_wp, &
      &  1.34624089204623_wp, -0.86133716815647_wp, -1.84340893849267_wp, &
      &  5.65596919189118_wp,  4.00172183859480_wp, -0.14131371969009_wp, &
      & 14.67430918222276_wp, -3.26230980007732_wp, -0.14344911021228_wp, &
      & 13.50897177220290_wp, -0.60815166181684_wp,  1.54898960808727_wp, &
      & 13.50780014200488_wp, -0.60614855212345_wp, -1.83214617078268_wp, &
      &  5.41408424778406_wp, -9.49239668625902_wp, -0.11022772492007_wp, &
      &  8.31919801555568_wp, -9.74947502841788_wp,  1.56539243085954_wp, &
      &  8.31511620712388_wp, -9.76854236502758_wp, -1.79108242206824_wp],&
      &  shape(xyz))
   real(wp), parameter :: vdwRad(4) = aatoau * [&
      & 1.45515_wp, 1.31125_wp, 1.24085_wp, 1.09155_wp]

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TSurfaceIntegrator) :: numint
   type(TNeighbourlist) :: neighList

   integer :: ii, jj, kk
   logical :: fail
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), allocatable :: sasa(:)
   real(wp), allocatable :: dsdr(:, :, :)
   real(wp), allocatable :: sl(:), sr(:), drdum(:, :, :)

   call init(env)
   call init(mol, at, xyz)
   call init(neighList, len(mol))

   allocate(sasa(nat))
   allocate(dsdr(3, nat, nat))

   call init(numint, env, vdwRad, 0.71032070_wp*aatoau)

   call env%check(fail)
   call assert(.not.fail)

   call assert_close(numint%cutoff, 13.097579816638_wp, thr)

   call neighList%generate(env, mol%xyz, numint%cutoff, trans, .false.)

   call env%check(fail)
   call assert(.not.fail)

   call numint%getSASA(neighList%neighs, neighList, mol%id, sasa, dsdr)

   call assert_close(sasa( 1), 1.2447730920594_wp, thr)
   call assert_close(sasa( 3), 3.7052686142047_wp, thr)
   call assert_close(sasa( 6), 4.8211340465422_wp, thr)
   call assert_close(sasa(10), 4.0746959304882_wp, thr)
   call assert_close(sasa(12), 7.8351317349355_wp, thr)
   call assert_close(sasa(17), 11.517101987449_wp, thr)
   call assert_close(sasa(22), 7.0235464411875_wp, thr)

   if (afail > 0) call terminate(afail)

   allocate(sr(nat))
   allocate(sl(nat))
   allocate(drdum(3, nat, nat))

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call numint%getSASA(neighList%neighs, neighList, mol%id, &
            & sr, drdum)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call numint%getSASA(neighList%neighs, neighList, mol%id, &
            & sl, drdum)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step

         do kk = 1, nat
            call assert_close(dsdr(jj, ii, kk), (sr(kk) - sl(kk))*step2, thr2)
         end do
      end do
   end do

   call terminate(afail)
end subroutine test_solv_sasaint

subroutine test_solv_gbsa
   use assertion
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : fourpi
   use xtb_mctc_convert, only : aatoau, kcaltoau
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_solv_gbsaparam
   use xtb_solv_gbsa
   implicit none
   real(wp), parameter :: thr = 1.0e-10_wp
   real(wp), parameter :: thr2 = 1.0e-8_wp
   integer, parameter :: nat = 24
   integer, parameter :: at(nat) = [6,7,6,7,6,6,6,8,7,6,8,7,6,6, &
      &                             1,1,1,1,1,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[ 2.02799738646442_wp,  0.09231312124713_wp, -0.14310895950963_wp, &
      &  4.75011007621000_wp,  0.02373496014051_wp, -0.14324124033844_wp, &
      &  6.33434307654413_wp,  2.07098865582721_wp, -0.14235306905930_wp, &
      &  8.72860718071825_wp,  1.38002919517619_wp, -0.14265542523943_wp, &
      &  8.65318821103610_wp, -1.19324866489847_wp, -0.14231527453678_wp, &
      &  6.23857175648671_wp, -2.08353643730276_wp, -0.14218299370797_wp, &
      &  5.63266886875962_wp, -4.69950321056008_wp, -0.13940509630299_wp, &
      &  3.44931709749015_wp, -5.48092386085491_wp, -0.14318454855466_wp, &
      &  7.77508917214346_wp, -6.24427872938674_wp, -0.13107140408805_wp, &
      & 10.30229550927022_wp, -5.39739796609292_wp, -0.13672168520430_wp, &
      & 12.07410272485492_wp, -6.91573621641911_wp, -0.13666499342053_wp, &
      & 10.70038521493902_wp, -2.79078533715849_wp, -0.14148379504141_wp, &
      & 13.24597858727017_wp, -1.76969072232377_wp, -0.14218299370797_wp, &
      &  7.40891694074004_wp, -8.95905928176407_wp, -0.11636933482904_wp, &
      &  1.38702118184179_wp,  2.05575746325296_wp, -0.14178615122154_wp, &
      &  1.34622199478497_wp, -0.86356704498496_wp,  1.55590600570783_wp, &
      &  1.34624089204623_wp, -0.86133716815647_wp, -1.84340893849267_wp, &
      &  5.65596919189118_wp,  4.00172183859480_wp, -0.14131371969009_wp, &
      & 14.67430918222276_wp, -3.26230980007732_wp, -0.14344911021228_wp, &
      & 13.50897177220290_wp, -0.60815166181684_wp,  1.54898960808727_wp, &
      & 13.50780014200488_wp, -0.60614855212345_wp, -1.83214617078268_wp, &
      &  5.41408424778406_wp, -9.49239668625902_wp, -0.11022772492007_wp, &
      &  8.31919801555568_wp, -9.74947502841788_wp,  1.56539243085954_wp, &
      &  8.31511620712388_wp, -9.76854236502758_wp, -1.79108242206824_wp],&
      &  shape(xyz))
   real(wp), parameter :: charges(nat) = [&
      &-0.05476196_wp, 0.00248420_wp, 0.08809027_wp,-0.28313223_wp, 0.12260575_wp, &
      &-0.04041664_wp, 0.26374283_wp,-0.43589031_wp,-0.11704515_wp, 0.31376962_wp, &
      &-0.44064154_wp,-0.07607650_wp,-0.05087676_wp,-0.04160170_wp, 0.06582772_wp, &
      & 0.08200852_wp, 0.08188349_wp, 0.05475583_wp, 0.09364446_wp, 0.07149007_wp, &
      & 0.07154523_wp, 0.09011030_wp, 0.06924195_wp, 0.06924253_wp]
   real(wp), parameter :: vdwRad(4) = aatoau * [&
      & 1.45515_wp, 1.31125_wp, 1.24085_wp, 1.09155_wp]

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TGeneralizedBorn) :: gbsa
   type(TNeighbourlist) :: neighList
   type(TGBSAData) :: input

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3), eps(3, 3), er, el
   real(wp), allocatable :: gradient(:, :), bornRad(:)
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   input = TGBSAData(&
      & dielectricConst = 7.6_wp, &
      & molecularMass = 72.1061_wp, density = 0.883_wp, &
      & bornScale = 1.51283258_wp, bornOffset = 0.000321234_wp * aatoau, &
      & freeEnergyShift = 2.13313378_wp * kcaltoau, &
      & probeRad = 0.97384049_wp * aatoau,&
      & descreening = &
      & [ 0.75304099_wp, 0.83443070_wp, 0.30000000_wp, 1.12338649_wp], &
      & surfaceTension = fourpi * 1.0e-5_wp * &
      & [-2.77665627_wp,-1.85054315_wp, 0.44495052_wp,-1.26534834_wp])

   allocate(gradient(3, nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call init(env)
   call init(mol, at, xyz)
   call init(neighList, len(mol))
   call init(gbsa, env, input, len(mol), vdwRad, cutoff=30.0_wp)

   call neighList%generate(env, mol%xyz, 30.0_wp, trans, .false.)

   call env%check(fail)
   call assert(.not.fail)

   call gbsa%update(neighList, mol%id)
   bornRad = gbsa%bornRad
   call gbsa%getEnergy(charges, charges, energy)
   call gbsa%getGradient(neighList, mol%id, charges, charges, gradient, sigma)

   call assert_close(energy, -0.36992594997156E-01_wp, thr)
   call assert_close(norm2(gradient), 0.35406332631899E-02_wp, thr)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbsa%update(neighList, mol%id)
         gbsa%bornRad = bornRad
         call gbsa%getEnergy(charges, charges, er)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbsa%update(neighList, mol%id)
         call gbsa%getEnergy(charges, charges, el)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step

         call assert_close(gradient(jj, ii), (er - el)*step2, thr2)
      end do
   end do

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call init(gbsa, env, input, len(mol), vdwRad, kernel=gbKernel%p16, &
      & cutoff=30.0_wp)
   call mol%update
   call neighList%update(mol%xyz, trans)
   call gbsa%update(neighList, mol%id)

   call gbsa%getEnergy(charges, charges, energy)
   call gbsa%getGradient(neighList, mol%id, charges, charges, gradient, sigma)

   call assert_close(energy, -0.37622153331423E-01_wp, thr)
   call assert_close(norm2(gradient), 0.37925825515919E-02_wp, thr)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbsa%update(neighList, mol%id)
         call gbsa%getEnergy(charges, charges, er)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbsa%update(neighList, mol%id)
         call gbsa%getEnergy(charges, charges, el)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step

         call assert_close(gradient(jj, ii), (er - el)*step2, thr2)
      end do
   end do

   call terminate(afail)
end subroutine test_solv_gbsa

subroutine test_solv_alpb
   use assertion
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : fourpi
   use xtb_mctc_convert, only : aatoau, kcaltoau
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init, len
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_solv_gbsaparam
   use xtb_solv_gbsa
   implicit none
   real(wp), parameter :: thr = 1.0e-10_wp
   real(wp), parameter :: thr2 = 1.0e-8_wp
   integer, parameter :: nat = 24
   integer, parameter :: at(nat) = [6,7,6,7,6,6,6,8,7,6,8,7,6,6, &
      &                             1,1,1,1,1,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[ 2.02799738646442_wp,  0.09231312124713_wp, -0.14310895950963_wp, &
      &  4.75011007621000_wp,  0.02373496014051_wp, -0.14324124033844_wp, &
      &  6.33434307654413_wp,  2.07098865582721_wp, -0.14235306905930_wp, &
      &  8.72860718071825_wp,  1.38002919517619_wp, -0.14265542523943_wp, &
      &  8.65318821103610_wp, -1.19324866489847_wp, -0.14231527453678_wp, &
      &  6.23857175648671_wp, -2.08353643730276_wp, -0.14218299370797_wp, &
      &  5.63266886875962_wp, -4.69950321056008_wp, -0.13940509630299_wp, &
      &  3.44931709749015_wp, -5.48092386085491_wp, -0.14318454855466_wp, &
      &  7.77508917214346_wp, -6.24427872938674_wp, -0.13107140408805_wp, &
      & 10.30229550927022_wp, -5.39739796609292_wp, -0.13672168520430_wp, &
      & 12.07410272485492_wp, -6.91573621641911_wp, -0.13666499342053_wp, &
      & 10.70038521493902_wp, -2.79078533715849_wp, -0.14148379504141_wp, &
      & 13.24597858727017_wp, -1.76969072232377_wp, -0.14218299370797_wp, &
      &  7.40891694074004_wp, -8.95905928176407_wp, -0.11636933482904_wp, &
      &  1.38702118184179_wp,  2.05575746325296_wp, -0.14178615122154_wp, &
      &  1.34622199478497_wp, -0.86356704498496_wp,  1.55590600570783_wp, &
      &  1.34624089204623_wp, -0.86133716815647_wp, -1.84340893849267_wp, &
      &  5.65596919189118_wp,  4.00172183859480_wp, -0.14131371969009_wp, &
      & 14.67430918222276_wp, -3.26230980007732_wp, -0.14344911021228_wp, &
      & 13.50897177220290_wp, -0.60815166181684_wp,  1.54898960808727_wp, &
      & 13.50780014200488_wp, -0.60614855212345_wp, -1.83214617078268_wp, &
      &  5.41408424778406_wp, -9.49239668625902_wp, -0.11022772492007_wp, &
      &  8.31919801555568_wp, -9.74947502841788_wp,  1.56539243085954_wp, &
      &  8.31511620712388_wp, -9.76854236502758_wp, -1.79108242206824_wp],&
      &  shape(xyz))
   real(wp), parameter :: charges(nat) = 1.0_wp/real(nat, wp) + [&
      &-0.05476196_wp, 0.00248420_wp, 0.08809027_wp,-0.28313223_wp, 0.12260575_wp, &
      &-0.04041664_wp, 0.26374283_wp,-0.43589031_wp,-0.11704515_wp, 0.31376962_wp, &
      &-0.44064154_wp,-0.07607650_wp,-0.05087676_wp,-0.04160170_wp, 0.06582772_wp, &
      & 0.08200852_wp, 0.08188349_wp, 0.05475583_wp, 0.09364446_wp, 0.07149007_wp, &
      & 0.07154523_wp, 0.09011030_wp, 0.06924195_wp, 0.06924253_wp]
   real(wp), parameter :: vdwRad(4) = aatoau * [&
      & 1.45515_wp, 1.31125_wp, 1.24085_wp, 1.09155_wp]

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TGeneralizedBorn) :: gbsa
   type(TNeighbourlist) :: neighList
   type(TGBSAData) :: input

   integer :: ii, jj
   logical :: fail
   real(wp) :: energy, sigma(3, 3), eps(3, 3), er, el
   real(wp), allocatable :: gradient(:, :), bornRad(:)
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-5_wp, step2 = 0.5_wp/step
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   input = TGBSAData(&
      & dielectricConst = 7.6_wp, &
      & molecularMass = 72.1061_wp, density = 0.883_wp, &
      & bornScale = 1.51283258_wp, bornOffset = 0.000321234_wp * aatoau, &
      & freeEnergyShift = 2.13313378_wp * kcaltoau, &
      & probeRad = 0.97384049_wp * aatoau,&
      & descreening = &
      & [ 0.75304099_wp, 0.83443070_wp, 0.30000000_wp, 1.12338649_wp], &
      & surfaceTension = fourpi * 1.0e-5_wp * &
      & [-2.77665627_wp,-1.85054315_wp, 0.44495052_wp,-1.26534834_wp])

   allocate(gradient(3, nat))
   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call init(env)
   call init(mol, at, xyz)
   call init(neighList, len(mol))
   call init(gbsa, env, input, len(mol), vdwRad, cutoff=30.0_wp, alpb=.true.)

   call neighList%generate(env, mol%xyz, 30.0_wp, trans, .false.)

   call env%check(fail)
   call assert(.not.fail)

   call gbsa%update(neighList, mol%id)
   bornRad = gbsa%bornRad
   call gbsa%getEnergy(charges, charges, energy)
   call gbsa%getGradient(neighList, mol%id, charges, charges, gradient, sigma)

   call assert_close(energy, -0.36992594997156E-01_wp, thr)
   call assert_close(norm2(gradient), 0.35406332631899E-02_wp, thr)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbsa%update(neighList, mol%id)
         gbsa%bornRad = bornRad
         call gbsa%getEnergy(charges, charges, er)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbsa%update(neighList, mol%id)
         call gbsa%getEnergy(charges, charges, el)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step

         call assert_close(gradient(jj, ii), (er - el)*step2, thr2)
      end do
   end do

   energy = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   call init(gbsa, env, input, len(mol), vdwRad, kernel=gbKernel%p16, &
      & cutoff=30.0_wp, alpb=.true.)
   call mol%update
   call neighList%update(mol%xyz, trans)
   call gbsa%update(neighList, mol%id)

   call gbsa%getEnergy(charges, charges, energy)
   call gbsa%getGradient(neighList, mol%id, charges, charges, gradient, sigma)

   call assert_close(energy, -0.37622153331423E-01_wp, thr)
   call assert_close(norm2(gradient), 0.37925825515919E-02_wp, thr)

   ! check numerical gradient
   do ii = 1, nat
      do jj = 1, 3
         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbsa%update(neighList, mol%id)
         call gbsa%getEnergy(charges, charges, er)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) - 2*step
         call mol%update
         call neighList%update(mol%xyz, trans)
         call gbsa%update(neighList, mol%id)
         call gbsa%getEnergy(charges, charges, el)

         mol%xyz(jj, ii) = mol%xyz(jj, ii) + step

         call assert_close(gradient(jj, ii), (er - el)*step2, thr2)
      end do
   end do

   call terminate(afail)
end subroutine test_solv_alpb
