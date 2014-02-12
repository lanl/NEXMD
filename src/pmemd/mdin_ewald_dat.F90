!*******************************************************************************
!
! Module: mdin_ewald_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module mdin_ewald_dat_mod

use file_io_dat_mod

  implicit none

! ew_type = 0 for pme, no other value supported.
! netfrc is 1 if remove average force (due to analytic forces in pme).
! use_pme = 0 basically implies a vacuum simulation and no switch function;
! vdwmeth = 0 for cutoff of vdw, 1 for analytic fix.
! verbose controls level of output. Look in force_info.
! frameon, chngmask - extra points options; read the manual &ewald section.

! Performance-related switches, newly added for PMEMD 10.  Most of these should
! default to reasonable values, but they are exposed to allow fine tuning and
! experimentation.  Some of the defaults can also be changed with defined
! constants in config.h.
!
! block_fft = -1 - auto select
!           =  0 - use slab fft
!           =  1 - use block fft; requires at least 4 processors, and not
!                  permitted for minimizations or if nrespa > 1.
! 
! When using block fft's, we essentially start with a collection of fft
! x-y slabs, and further divide each slab into fft_blk_y_divisor "blocks" in
! the y dimension.  Thus each "block is a collection contiguous x-runs. The
! y divisor must be a value between 2 and nfft2, where nfft2 is the nfft2 value
! AFTER the axis optimization operation has been performed if it is allowed
! (thus it is really only safe to either use min(nfft1,nfft2,nfft3), or turn
! off axis optimization, or think carefully about what axis optimization will
! actually do...  In all instances tested so far, relatively small values work
! best for fft_blk_y_divisor.  This value is used only if block_fft .eq. 1.
!
! fft_blk_y_divisor = 2 .. nfft2 (reoriented);
!                     default=2 or 4 depending on numtasks. 
!
! excl_recip = 0..1 - Exclusive reciprocal tasks flag.  This flag, when 1,
!                     specifies that tasks that do reciprocal force calcs will
!                     not also do direct force calculations.  This has some
!                     benefits at higher task count.  At lower task count,
!                     setting this flag can result in significant
!                     underutilization of reciprocal tasks.  This flag will
!                     automatically be cleared if block fft's are not in use.
!
! excl_master = 0..1 - Exclusive master task flag.  This flag, when 1,
!                      specifies that the master task will not do force and
!                      energy calculations.  At high scaling, what this does
!                      is insure that no tasks are waiting for the master to
!                      initiate collective communications events.  The master
!                      is thus basically dedicated to handling loadbalancing and
!                      output.  At lower task count, this is obviously
!                      wasteful.  This flag will automatically be cleared if
!                      block fft's are not in use or if excl_recip .ne. 1.
!
!                      AND NOTE - when block fft's are in use, that implies that
!                      you are not doing a minimization and not using
!                      nrespa > 1.
!
! atm_redist_freq = 16..1280 - The frequency (in pairlist build events) for
!                              reassigning atom ownership to tasks.  As a run
!                              progresses, diffusion causes the atoms originally
!                              collocated and assigned to one task to occupy a
!                              larger volume.  With time, this starts to cause
!                              a higher communications load, though the impact
!                              is lower than one might expect.  Currently, by
!                              default we reassign atoms to tasks every 320
!                              pairlist builds at low to medium task count and
!                              we reassign atoms to tasks every 32 pairlist
!                              builds at higher task counts (currently defined
!                              as >= 96 tasks, redefinable in config.h).  The
!                              user can however specify the specific value he
!                              desires.  At low task count, frequent atom
!                              redistribution tends to have a noticeable cost
!                              and little benefit. At higher task count, the
!                              cost is lower and the benefit is higher.

  integer, parameter    :: mdin_ewald_int_cnt = 16

  integer                       nfft1, nfft2, nfft3, bspl_order, netfrc, &
                                use_axis_opt, vdwmeth, verbose, frameon, &
                                chngmask, block_fft, fft_blk_y_divisor, &
                                excl_recip, excl_master, atm_redist_freq, &
                                use_pme

  common / mdin_ewald_int /     nfft1, nfft2, nfft3, bspl_order, netfrc, &
                                use_axis_opt, vdwmeth, verbose, frameon, &
                                chngmask, block_fft, fft_blk_y_divisor, &
                                excl_recip, excl_master, atm_redist_freq, &
                                use_pme
                                

  integer, parameter    :: mdin_ewald_dbl_cnt = 5

  save  :: /mdin_ewald_int /

  double precision              dsum_tol, ew_coeff, skinnb, eedtbdns, &
                                fft_grids_per_ang

  common / mdin_ewald_dbl /     dsum_tol, ew_coeff, skinnb, eedtbdns, &
                                fft_grids_per_ang

  save  :: /mdin_ewald_dbl /

  integer, parameter            :: gridlo = 6
  integer, parameter            :: gridhi = 512

  double precision, parameter   :: ew_coefflo = 0.1d0
  double precision, parameter   :: ew_coeffhi = 0.7d0

  integer, parameter            :: orderlo = 3
  integer, parameter            :: orderhi = 25
  
  double precision, parameter   :: skinlo = 0.d0
  double precision, parameter   :: skinhi = 5.d0
  
  double precision, parameter   :: denslo = 100.d0
  double precision, parameter   :: denshi = 25000.d0

  double precision, parameter   :: boxlo = 1.d0
  double precision, parameter   :: boxhi = 1000.d0

  double precision, parameter   :: anglo = 30.d0
  double precision, parameter   :: anghi = 150.d0

  double precision, parameter   :: tollo = 1.d-12
  double precision, parameter   :: tolhi = 1.d-2

  double precision, parameter   :: fft_grids_per_ang_lo = 0.25d0
  double precision, parameter   :: fft_grids_per_ang_hi = 10.d0

! Namelist values that are local because they only have one legal value, or are
! not legitimately used.

  integer, private              :: ew_type, mlimit(3), nbflag, nbtell, &
                                   eedmeth, ee_type, ischrgd, indmeth, &
                                   maxiter, irstdip, scaldip, &
                                   column_fft

  double precision, private     :: rsum_tol, diptol, dipmass, diptau

! Namelist entries that are returned with disambiguated names for clarity:

  double precision, private     :: a, b, c
  double precision, private     :: alpha, beta, gamma
  integer, private              :: order

  private       :: ewald

  namelist /ewald/      nfft1, nfft2, nfft3, order, verbose, ew_type, &
                        dsum_tol, rsum_tol, mlimit, ew_coeff, nbflag, &
                        skinnb, nbtell, netfrc, use_pme, vdwmeth, &
                        eedmeth, eedtbdns, ee_type, ischrgd, frameon, &
                        chngmask, indmeth, diptol, maxiter, dipmass, &
                        diptau, irstdip, scaldip, alpha, beta, gamma, &
                        a, b, c, use_axis_opt, fft_grids_per_ang, column_fft, &
                        block_fft, fft_blk_y_divisor, excl_recip, excl_master, &
                        atm_redist_freq

! Private subroutines:

  private       compute_nfft, compute_even_nfft, compute_ew_coeff

contains

!*******************************************************************************
!
! Subroutine:  init_mdin_ewald_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_mdin_ewald_dat(ew_box_alpha, ew_box_beta, ew_box_gamma, &
                               ew_box)

! CODING NOTE: Initialization of the common blocks is NOT complete after this
!              subroutine executes.  There is additional initialization in
!              subroutine TBS which depends on external values not defined at
!              the time of execution of this subroutine.  So don't broadcast
!              the common areas until that initialization is complete!

  use axis_optimize_mod
  use file_io_mod
  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  double precision, intent(in out)      :: ew_box_alpha
  double precision, intent(in out)      :: ew_box_beta
  double precision, intent(in out)      :: ew_box_gamma
  double precision, intent(in out)      :: ew_box(3)

! Local variables:

  integer               :: ifind
  integer               :: inerr
  double precision      :: erfc_val, y
  integer               :: excl_recip_switch_val
  integer               :: excl_master_switch_val
  integer               :: atm_redist_switch_val ! from lo freq to hi freq

  inerr = 0

  nfft1 = 0
  nfft2 = 0
  nfft3 = 0
  order = 4
  verbose = 0
  ew_type = 0
  dsum_tol = 1.d-5
  rsum_tol = 5.d-5
  mlimit(:) = 0
  ew_coeff = 0.d0
  nbflag = 1
! Optimal values for skinnb selected below (or best guess anyway).
#ifdef MPI
  skinnb = 2.d0
#else
  skinnb = 1.d0
#endif

! Override default optimal values if SKINNB_DEF is defined.  It
! must have a real value.

#ifdef SKINNB_DEF
  skinnb = SKINNB_DEF
#endif

#ifdef EXCL_RECIP_SWITCH
  excl_recip_switch_val = EXCL_RECIP_SWITCH
#else
  excl_recip_switch_val = 96
#endif

#ifdef EXCL_MASTER_SWITCH
  excl_master_switch_val = EXCL_MASTER_SWITCH
#else
  excl_master_switch_val = 96
#endif

#ifdef ATM_REDIST_SWITCH
  atm_redist_switch_val = ATM_REDIST_SWITCH
#else
  atm_redist_switch_val = 96
#endif

  nbtell = 0
  netfrc = -1                   ! -1 means "set to default value", once you
                                ! know imin.
  use_pme = 1
  vdwmeth = 1
  eedmeth = 1
  !If IPS is turned on, set in cntrl namelist then we need to turn off PME
  !so we will override the default here.
  !vdwmeth and eedmeth are set to 2 and 6, respectively, to match sander
  if (ips .gt. 0)then
    use_pme = 0
    vdwmeth = 2
    eedmeth = 6
  endif
  eedtbdns = 5000.d0
  ee_type = 1
  ischrgd = RETIRED_INPUT_OPTION
  frameon = 1
  chngmask = 1
  indmeth = 3                   ! Not really supported...
  diptol = 1.d-4                ! Not really supported...
  maxiter = 20                  ! Not really supported...
  dipmass = 0.33d0              ! Not really supported...
  diptau = 11.d0                ! Not really supported...
  irstdip = 0                   ! Not really supported...
  scaldip = 1                   ! Not really supported...
  alpha = 0.d0
  beta = 0.d0
  gamma = 0.d0
  a = 0.d0
  b = 0.d0
  c = 0.d0
  use_axis_opt = -1             ! -1 means "set to default value, once you have
                                ! final box info"
  fft_grids_per_ang = 1.d0

  ! PMEMD-specific performance options:

  block_fft = -1                ! auto-select
  fft_blk_y_divisor = -1        ! auto-select
  excl_recip = -1               ! auto-select
  excl_master = -1              ! auto-select
  atm_redist_freq = -1          ! auto-select

  column_fft = 0                ! Unsupported sander performance option.

  rewind(mdin)                  ! Insurance against maintenance mods.

  call nmlsrc('ewald', mdin, ifind)

  if (ifind .ne. 0) read(mdin, nml = ewald)   ! Namelist found. Read it:

  bspl_order = order    ! Store order by unambiguous name:

  ! If there was box length or angle input, use it to override the global
  ! values unless ntx is 6 or 7.

  if (ntx .ne. 6 .and. ntx .ne. 7) then
    if (alpha .ne. 0.d0) ew_box_alpha = alpha
    if (beta .ne. 0.d0)  ew_box_beta = beta
    if (gamma .ne. 0.d0) ew_box_gamma = gamma
    if (a .ne. 0.d0)     ew_box(1) = a
    if (b .ne. 0.d0)     ew_box(2) = b
    if (c .ne. 0.d0)     ew_box(3) = c
  endif

  ! For nonorthogonal unit cells, we force axis flipping off; we don't handle
  ! changing angles in the axis flipping code.

  if (ew_box_alpha .ne. 90.d0 .or. &
      ew_box_beta .ne. 90.d0 .or. &
      ew_box_gamma .ne. 90.d0) then
    if (use_axis_opt .gt. 0) then
      write(mdout, '(a,a,/)') info_hdr, &
        'Axis order optimization not used - unit cell is nonorthogonal.'
    end if
    use_axis_opt = 0
  end if

  ! Now we know enough to set a default use_axis_opt, assuming the user did
  ! not set a value.  We don't select axis optimizations if the likelihood of
  ! hotspots is high (minimizations ) because results will not be as comparable
  ! (probably due to slight differences in mpi fft on different axes that get
  ! magnified in shake; axis flipping is nicely reproducible in sane systems).
  ! We also don't select it unless the box has at least a 3:2 aspect ratio, as
  ! it won't make much difference in relatively cubic systems.

#ifdef MPI
  ! For minimizations, anything using random number generation, and boxes
  ! that are not particularly oblong we turn off axis optimization.  This is
  ! a conservative approach, insuring better consistency for simulations
  ! likely to have hotspots, avoiding problems with the sequence of random
  ! number generations (and random number techniques may also generate
  ! hotspots, I would think), and basically only turning it on where there
  ! is a clear benefit.
    
  if (use_axis_opt .eq. -1) then
    if (imin .ne. 0 .or. &
        ntx .eq. 1 .or. ntx .eq. 2 .or. &
        ntt .eq. 2 .or. ntt .eq. 3 .or. &
        maxval(ew_box(:))/minval(ew_box(:)) .lt. 1.5d0) then
      use_axis_opt = 0
    else
      use_axis_opt = 1
    end if
  end if
#else
  if (use_axis_opt .eq. -1)  use_axis_opt = 0
#endif

  if (use_axis_opt .ne. 0) then

    ! Okay, one last check.  If axis flipping produces an odd value for
    ! nfft1, then it can't be used. We have to set up the axis flips, and
    ! revert them if nfft1 would be odd.

    call setup_axis_opt(ew_box(1), ew_box(2), ew_box(3))

    if (.not. flipped_nfft1_even()) then
      call setup_axis_opt(1.d0, 1.d0, 1.d0)        ! Basically no flipping done.
      use_axis_opt = 0
      write(mdout, '(a,a,/)') info_hdr, &
        'Axis order optimization not used - nfft1 would be odd.'
    else
      write(mdout, '(a,a,/)') info_hdr, 'Axis order optimization will be used.'
      call axes_flip(ew_box(1), ew_box(2), ew_box(3))
      ! After this point, ew_box has been flipped.  Other unflipped copies of 
      ! it are used for output in the setup printout.
      call int_axes_flip(nfft1, nfft2, nfft3)
    end if

  else

    call setup_axis_opt(1.d0, 1.d0, 1.d0)        ! Basically no flipping done.

  end if

  ! The above flipping needed to occur before computing nfft1..3


#ifdef CUDA
  ! CUDA needs all dimensions divisible by 4 (but supports powers of 7) ...
  if (nfft1 .eq. 0) &
    call compute_cuda_nfft(ew_box(1), fft_grids_per_ang, nfft1)
  if (nfft2 .eq. 0) &
    call compute_cuda_nfft(ew_box(2), fft_grids_per_ang, nfft2)
  if (nfft3 .eq. 0) &
    call compute_cuda_nfft(ew_box(3), fft_grids_per_ang, nfft3)
#else
  ! RCFFT needs an even dimension in the x direction (nfft1) ...

  if (nfft1 .eq. 0) &
    call compute_even_nfft(ew_box(1), fft_grids_per_ang, nfft1)
  if (nfft2 .eq. 0) &
    call compute_nfft(ew_box(2), fft_grids_per_ang, nfft2)
  if (nfft3 .eq. 0) &
    call compute_nfft(ew_box(3), fft_grids_per_ang, nfft3)
#endif

  ! Assign ewald coefficient.  The input values are checked here to insure
  ! we don't have a fp error...

  if (ew_coeff .lt. 1.d-6) then
    call float_legal_range('dsum_tol', dsum_tol, tollo, tolhi, inerr)
    if (inerr .eq. 1) call mexit(6, 1)
    call compute_ew_coeff(es_cutoff) ! depends on dsum_tol
  else
    call float_legal_range('ew_coeff', ew_coeff, ew_coefflo, ew_coeffhi, inerr)
    if (inerr .eq. 1) call mexit(6, 1)
    y = ew_coeff * es_cutoff
    call derfcfun(y, erfc_val)
    dsum_tol = erfc_val / es_cutoff
  end if

  ! Select performance defaults if they were not made explicit...

#ifdef MPI

  if (block_fft .lt. 0) then
    if (nfft3 / numtasks .ge. 2) then
      block_fft = 0
    else
      block_fft = 1
    end if
  end if

  ! Hammer block_fft to allowable values if there is a problem with processor
  ! count...

  if (block_fft .eq. 1 .and. numtasks .lt. 4) block_fft = 0

  ! Reset block_fft for unallowed conditions:

  if (block_fft .eq. 1 .and. imin .eq. 1) block_fft = 0
  if (block_fft .eq. 1 .and. nrespa .gt. 1) block_fft = 0

  if (fft_blk_y_divisor .le. 0) fft_blk_y_divisor = 4

  if (numtasks .le. 16) fft_blk_y_divisor = 2

  if (excl_recip .lt. 0) then
    if (numtasks .lt. excl_recip_switch_val) then
      excl_recip = 0
    else
      excl_recip = 1
    end if
  end if

  ! Reset excl_recip if not using block fft's:

  if (block_fft .eq. 0) excl_recip = 0

  if (excl_master .lt. 0) then
    if (numtasks .lt. excl_master_switch_val) then
      excl_master = 0
    else
      excl_master = 1
    end if
  end if

  ! Reset excl_master if not using block fft's:

  if (block_fft .eq. 0) excl_master = 0
  if (excl_recip .eq. 0) excl_master = 0

  if (atm_redist_freq .le. 0) then
    if (numtasks .lt. atm_redist_switch_val) then
      atm_redist_freq = 320
    else
      atm_redist_freq = 32
    end if
  end if

#else
  ! These are not used in the uniprocessor implementation; set to 0 to 
  ! insure no grief...
  block_fft = 0         ! legal value; meaningless for uniprocessor code...
  fft_blk_y_divisor = 2 ! legal value; meaningless for uniprocessor code...
  excl_recip = 0        ! legal value, meaningless for uniprocessor code...
  excl_master = 0       ! legal value, meaningless for uniprocessor code...
  atm_redist_freq = 0   ! legal value; meaningless for uniprocessor code...
#endif

  ! Determine netfrc value to use.

  if (netfrc .lt. 0) then
    if (imin .eq. 0) then
      netfrc = 1
    else
      netfrc = 0
    end if
  end if

  return

contains

!*******************************************************************************
!
! Internal Function:  flipped_nfft1_even
!
! Description: Check to see if flipping the nfft's will produce an odd nfft1.
!              
!*******************************************************************************

function flipped_nfft1_even()

  implicit none

! Return value:

  logical       :: flipped_nfft1_even

! Local variables:

  integer       :: test_nfft1, test_nfft2, test_nfft3

  test_nfft1 = nfft1
  test_nfft2 = nfft2
  test_nfft3 = nfft3

  call int_axes_flip(test_nfft1, test_nfft2, test_nfft3)

  flipped_nfft1_even = (mod(test_nfft1, 2) .eq. 0)

end function flipped_nfft1_even

end subroutine init_mdin_ewald_dat

!*******************************************************************************
!
! Subroutine:  validate_mdin_ewald_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine validate_mdin_ewald_dat(ew_box_alpha, ew_box_beta, ew_box_gamma, &
                                   ew_box)

  use fft1d_mod
  use gbl_constants_mod
  use mdin_ctrl_dat_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use charmm_mod, only : charmm_active

  implicit none

! Formal arguments:

  double precision, intent(in)  :: ew_box_alpha
  double precision, intent(in)  :: ew_box_beta
  double precision, intent(in)  :: ew_box_gamma
  double precision, intent(in)  :: ew_box(3)

! Local variables

  integer       :: inerr

! Check on bogus data and unsupported options:

  inerr = 0

  call int_legal_range('nfft1', nfft1, gridlo, gridhi, inerr)
  call int_legal_range('nfft2', nfft2, gridlo, gridhi, inerr)
  call int_legal_range('nfft3', nfft3, gridlo, gridhi, inerr)
  
  call test_prime_factors('nfft1', nfft1, inerr)
  call test_prime_factors('nfft2', nfft2, inerr)
  call test_prime_factors('nfft3', nfft3, inerr)

  call int_legal_range('order (bspline interpolation order) ', order, &
                       orderlo, orderhi, inerr)

  call int_legal_range('verbose', verbose, 0, 3, inerr)

  call int_legal_range('use_pme', use_pme, 0, 1, inerr)

  if (ew_type .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'ew_type must be 0 for PMEMD!'
    inerr = 1
  end if

  if (rsum_tol .ne. 5.d-5) then
    write(mdout, '(a,a)') error_hdr, 'rsum_tol is only used if ew_type == 1!'
    inerr = 1
  end if

  if (mlimit(1) .ne. 0 .or. mlimit(2) .ne. 0 .or. mlimit(3) .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'mlimit is only used if ew_type == 1!'
    inerr = 1
  end if

  ! Check for user input of non-default values of nbflag and nbtell.  These
  ! switches support old list processing/debugging output, and are no longer
  ! really used.  No need to store values.

  if (nbflag .ne. 1) then
    write(mdout, '(a,a)') warn_hdr, &
      'nbflag is ignored; skin checks are always done.'
  end if

  call float_legal_range('skinnb', skinnb, skinlo, skinhi, inerr)

  if (nbtell .ne. 0) then
    write(mdout, '(a,a)') info_hdr, 'PMEMD ignores nbtell != 0.'
  end if

  call int_legal_range('netfrc', netfrc, 0, 1, inerr)

  call int_legal_range('vdwmeth', vdwmeth, 0, 2, inerr)

  if (eedmeth .ne. 1) then
    if (eedmeth .ne. 6) then
      write(mdout, '(a,a)') error_hdr, 'eedmeth must be 1 for PME or 6 for IPS!'
      inerr = 1
    end if
  end if

  call float_legal_range('eedtbdns', eedtbdns, denslo, denshi, inerr)

  if (ee_type .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'ee_type must be 1 for PMEMD!'
    inerr = 1
  end if

  if (imin .eq. 1 .and. netfrc .eq. 1) then
    write(mdout, '(a,a)') warn_hdr, &
      'Use of netfrc == 1 in minimizations is not recommended!'
  end if

  if (ischrgd .ne. RETIRED_INPUT_OPTION) then
    write(mdout, '(a,a)') warn_hdr, &
      'The ischrgd ewald option is deprecated and ignored.'
    write(mdout, '(a,a)') extra_line_hdr, &
      'A nonzero net charge due to roundoff error is always neutralized.'
  end if

  call int_legal_range('frameon', frameon, 0, 1, inerr)

  call int_legal_range('chngmask', chngmask, 0, 1, inerr)

  ! If CHARMM force field is active we have to set chngmask to 0 to
  ! stop the code rebuilding the 1-4 list.
  ! If we are using the charmm force field we want the exlcuded atom list in
  ! the prmtop file to be used and NOT rebuilt. This is forced
  ! by setting the ewald namelist variable chngmask to 0.
  ! RCW: Not sure we strictly need this but we will do it this way for the time being
  !      to avoid confusion.
  if (charmm_active) then
     write(mdout,'(a)') '|CHARMM: Overriding default value of chngmask.'
     write(mdout,'(a)') '|CHARMM: Setting chngmask = 0.'

     chngmask = 0
  end if

  if (indmeth .ne. 3) then
    write(mdout, '(a,a)') error_hdr, 'indmeth is only used if ipol != 0!'
    inerr = 1
  end if

  if (diptol .ne. 1.d-4) then
    write(mdout, '(a,a)') error_hdr, 'diptol is only used if ipol != 0!'
    inerr = 1
  end if

  if (maxiter .ne. 20) then
    write(mdout, '(a,a)') error_hdr, 'maxiter is only used if ipol != 0!'
    inerr = 1
  end if

  if (dipmass .ne. 0.33d0) then
    write(mdout, '(a,a)') error_hdr, 'dipmass is only used if ipol != 0!'
    inerr = 1
  end if

  if (diptau .ne. 11.d0) then
    write(mdout, '(a,a)') error_hdr, 'diptau is only used if ipol != 0!'
    inerr = 1
  end if

  if (irstdip .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'irstdip is only used if ipol != 0!'
    inerr = 1
  end if

  if (scaldip .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'scaldip is only used if ipol != 0!'
    inerr = 1
  end if

  call float_legal_range('a', ew_box(1), boxlo, boxhi, inerr)
  call float_legal_range('b', ew_box(2), boxlo, boxhi, inerr)
  call float_legal_range('c', ew_box(3), boxlo, boxhi, inerr)
  call float_legal_range('alpha', ew_box_alpha, anglo, anghi, inerr)
  call float_legal_range('beta',  ew_box_beta, anglo, anghi, inerr)
  call float_legal_range('gamma', ew_box_gamma, anglo, anghi, inerr)
  call float_legal_range('fft_grids_per_ang', fft_grids_per_ang, &
                         fft_grids_per_ang_lo, fft_grids_per_ang_hi, inerr)

#ifdef MPI
  call int_legal_range('block_fft', block_fft, 0, 1, inerr)
  call int_legal_range('fft_blk_y_divisor', fft_blk_y_divisor, 2, nfft2, inerr)
  call int_legal_range('excl_recip', excl_recip, 0, 1, inerr)
  call int_legal_range('excl_master', excl_master, 0, 1, inerr)
  call int_legal_range('atm_redist_freq', atm_redist_freq, 16, 1280, inerr)
#endif

  call int_legal_range('use_axis_opt', use_axis_opt, 0, 1, inerr)

  if (column_fft .ne. 0) then
    write(mdout, '(a,a)') info_hdr, &
      'The column_fft ewald option is not supported and ignored.'
    write(mdout, '(a,a)') extra_line_hdr, &
      'This is a sander-specific optimization option.'
  end if

#ifdef CUDA

  ! Trap things in the ewald namelist that the CUDA implementation does not
  ! support.

  if (order > 4) then
    write(mdout, '(a)') 'CUDA (GPU): Implementation does not support interpolation orders of greater than 4.'
    write(mdout, '(a)') '            Require order <= 4.'
    inerr = 1
  end if

#endif

! Field any errors and bag out.

  if (inerr .eq. 1) then
    write(mdout, '(/,a)') ' Input errors occurred. Terminating execution.'
    call mexit(6, 1)
  else
    write(mdout, '(a)') ' '
  end if

  return

end subroutine validate_mdin_ewald_dat

!*******************************************************************************
!
! Subroutine:  int_legal_range
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine int_legal_range(string, param, lo, hi, inerr)

  use gbl_constants_mod

  implicit none

  character(*)  :: string
  integer       :: param
  integer       :: lo
  integer       :: hi
  integer       :: inerr

  if (param .lt. lo .or. param .gt. hi) then
    write(mdout, '(a,a,a,i8,a,i8,a)') error_hdr, &
      string, ' must be in the range of ', lo, ' to ', hi, '!'
    inerr = 1
  end if

  return

end subroutine int_legal_range

!*******************************************************************************
!
! Subroutine:  float_legal_range
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine float_legal_range(string, param, lo, hi, inerr)

  use gbl_constants_mod

  implicit none

  character(*)          :: string
  double precision      :: param
  double precision      :: lo
  double precision      :: hi
  integer               :: inerr

  if (param .lt. lo .or. param .gt. hi) then
    write(mdout, '(a,a,a,e12.5,a,e12.5,a)') error_hdr, &
      string, ' must be in the range of ', lo, ' to ', hi, '!'
    inerr = 1
  end if

  return

end subroutine float_legal_range

!*******************************************************************************
!
! Subroutine:  print_mdin_ewald_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine print_mdin_ewald_dat(alpha, beta, gamma, ew_box, es_cut)

! NOTE: The box data is not really part of the ewald data, but is passed in
!       because we must print it here for consistency with sander.  It should
!       already be axis-flipped if necessary, and otherwise represents the
!       values stored in pbc_mod.

  use axis_optimize_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: alpha
  double precision, intent(in)  :: beta
  double precision, intent(in)  :: gamma
  double precision, intent(in)  :: ew_box(3)
  double precision, intent(in)  :: es_cut

! Local variables:

  integer               :: ord1, ord2, ord3
  integer               :: nffts(3)

  if (numextra .gt. 0) then
    write(mdout,'(/a)') 'Extra-points options:'
    write(mdout,'(5x,4(a,i8))') 'frameon =', frameon, &
            ', chngmask=', chngmask
  end if

! If the axis flipping optimization is in effect, we want to restore the
! box lengths and nfft's to original values for printout...

  ord1 = axis_flipback_ords(1)
  ord2 = axis_flipback_ords(2)
  ord3 = axis_flipback_ords(3)

  nffts(1) = nfft1
  nffts(2) = nfft2
  nffts(3) = nfft3

  write(mdout,'(/a)') 'Ewald parameters:'
  write(mdout,'(5x,4(a,i8))') 'verbose =', verbose, &
          ', ew_type =', ew_type, ', nbflag  =', nbflag, &
          ', use_pme =', use_pme
  write(mdout,'(5x,4(a,i8))') 'vdwmeth =', vdwmeth, &
          ', eedmeth =', eedmeth,', netfrc  =', netfrc
  write(mdout, 9002) ew_box(ord1), ew_box(ord2), ew_box(ord3)
  write(mdout, 9003) alpha, beta, gamma
  write(mdout, 9004) nffts(ord1), nffts(ord2), nffts(ord3)
  write(mdout, 9006) es_cut, dsum_tol
  write(mdout, 9007) ew_coeff
  write(mdout, 9005) order

#ifdef MPI
  write(mdout,'(/a)') '| PMEMD ewald parallel performance parameters:'
  write(mdout, 9008) block_fft
  write(mdout, 9009) fft_blk_y_divisor
  write(mdout, 9010) excl_recip
  write(mdout, 9011) excl_master
  write(mdout, 9012) atm_redist_freq
#endif /* MPI */

  return

9002 format (5x,'Box X =',f9.3,3x,'Box Y =',f9.3,3x,'Box Z =',f9.3)
9003 format (5x,'Alpha =',f9.3,3x,'Beta  =',f9.3,3x,'Gamma =',f9.3)
9004 format (5x,'NFFT1 =',i5  ,7x,'NFFT2 =',i5  ,7x,'NFFT3 =',i5)
9005 format (5x,'Interpolation order =',i5)
9006 format (5x,'Cutoff=',f9.3,3x,'Tol   =',e9.3)
9007 format (5x,'Ewald Coefficient =',f9.5)
#ifdef MPI
9008 format ('|     block_fft =',i5)
9009 format ('|     fft_blk_y_divisor =',i5)
9010 format ('|     excl_recip =',i5)
9011 format ('|     excl_master =',i5)
9012 format ('|     atm_redist_freq =',i5)
#endif /* MPI */

end subroutine print_mdin_ewald_dat

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_mdin_ewald_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_mdin_ewald_dat

  use parallel_dat_mod

  implicit none

  call mpi_bcast(nfft1, mdin_ewald_int_cnt, mpi_integer, &
                 0, pmemd_comm, err_code_mpi)

  call mpi_bcast(dsum_tol, mdin_ewald_dbl_cnt, mpi_double_precision, &
                 0, pmemd_comm, err_code_mpi)

  return

end subroutine bcast_mdin_ewald_dat
#endif

!*******************************************************************************
!
! Subroutine:  compute_nfft
!
! Description:  This routine uses check_prime_factors to get the smallest
!               integer that will give at least the indicated grids per angstrom
!               given the specified box length (one axis, of course), with the
!               restriction that the integer have only prime factors of
!               2, 3, 5, and for fftw only, 7 (basically, check_prime_factors
!               is provided by the fft implementation being used, and that
!               determines the acceptable prime factors).
!*******************************************************************************

subroutine compute_nfft(box_len, grids_per_ang, nfft)

  use fft1d_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: box_len
  double precision, intent(in)  :: grids_per_ang
  integer, intent(out)          :: nfft

! Local variables:

  integer               :: candidate_grid_cnt
  integer               :: ret_code

  candidate_grid_cnt = ceiling(box_len * grids_per_ang)

  do
    call check_prime_factors(candidate_grid_cnt, ret_code)
    if (ret_code .eq. 1) then
      nfft = candidate_grid_cnt
      exit
    end if
    candidate_grid_cnt = candidate_grid_cnt + 1
  end do

  return

end subroutine compute_nfft

!*******************************************************************************
!
! Subroutine:  compute_even_nfft
!
! Description:  This routine uses check_prime_factors to get the smallest
!               integer that will give at least the indicated grids per angstrom
!               given the specified box length (one axis, of course), with the
!               restriction that the integer have only prime factors of
!               2, 3, 5, and for fftw only, 7 (basically, check_prime_factors
!               is provided by the fft implementation being used, and that
!               determines the acceptable prime factors).  This routine also
!               REQUIRES that one of the prime factors is two.
!
!*******************************************************************************

subroutine compute_even_nfft(box_len, grids_per_ang, nfft)

  use fft1d_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: box_len
  double precision, intent(in)  :: grids_per_ang
  integer, intent(out)          :: nfft

! Local variables:

  integer               :: candidate_grid_cnt
  integer               :: ret_code

  candidate_grid_cnt = ceiling(box_len * grids_per_ang)
  if (mod(candidate_grid_cnt, 2) .ne. 0) &
    candidate_grid_cnt = candidate_grid_cnt + 1
  do
    call check_prime_factors(candidate_grid_cnt, ret_code)
    if (ret_code .eq. 1) then
      nfft = candidate_grid_cnt
      exit
    end if
    candidate_grid_cnt = candidate_grid_cnt + 2
  end do

  return

end subroutine compute_even_nfft

#ifdef CUDA
!*******************************************************************************
!
! Subroutine:  compute_cuda_nfft
!
! Description:  This routine uses check_prime_factors to get the smallest
!               integer that will give at least the indicated grids per angstrom
!               given the specified box length (one axis, of course), with the
!               restriction that the integer is a multiple of 4 which only has prime 
!               factors of 2, 3, 5, and 7 (basically, check_prime_factors
!               is provided by the fft implementation being used, and that
!               determines the acceptable prime factors).
!
!*******************************************************************************

subroutine compute_cuda_nfft(box_len, grids_per_ang, nfft)

  use fft1d_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: box_len
  double precision, intent(in)  :: grids_per_ang
  integer, intent(out)          :: nfft

! Local variables:

  integer               :: candidate_grid_cnt
  integer               :: ret_code

  candidate_grid_cnt = ((ceiling(box_len * grids_per_ang) + 3) / 4) * 4
  
! Clamp to powers of 2 if close to one of them
  if ((candidate_grid_cnt .ge. 60) .and. (candidate_grid_cnt .le. 68)) then
    candidate_grid_cnt = 64  
  else if ((candidate_grid_cnt .ge. 120) .and. (candidate_grid_cnt .le. 136)) then
    candidate_grid_cnt = 128  
  else if ((candidate_grid_cnt .ge. 240) .and. (candidate_grid_cnt .le. 272)) then
    candidate_grid_cnt = 256    
  else if ((candidate_grid_cnt .ge. 480) .and. (candidate_grid_cnt .le. 544)) then
    candidate_grid_cnt = 512    
  else if ((candidate_grid_cnt .ge. 964) .and. (candidate_grid_cnt .le. 1088)) then
    candidate_grid_cnt = 1024
  end if    
  
  do
    call check_prime_factors(candidate_grid_cnt, ret_code)
    if (ret_code .eq. 1) then
      nfft = candidate_grid_cnt
      exit
    end if
    candidate_grid_cnt = candidate_grid_cnt + 4
  end do

  return

end subroutine compute_cuda_nfft
#endif

!*******************************************************************************
!
! Subroutine:  compute_ew_coeff
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine compute_ew_coeff(cutoff)

  implicit none

! Formal arguments:

  double precision      :: cutoff

! Local variables:

  integer               :: i, n

  double precision      :: erfc_val
  double precision      :: term, x, xlo, xhi, y

! First get direct sum tolerance. How big must ew_coeff be to get
! terms outside the cutoff below tol?

  x = 0.5d0
  i = 0

  do
    x = 2.d0 * x
    i = i + 1
    y = x * cutoff
    call derfcfun(y, erfc_val)
    term = erfc_val/cutoff
    if (term .lt. dsum_tol) exit
  end do

! binary search tolerance is 2 to the -50th

  n = i + 50
  xlo = 0.d0
  xhi = x

  do i = 1, n
    x = (xlo + xhi)/2
    y = x * cutoff
    call derfcfun(y, erfc_val)
    term = erfc_val/cutoff
    if (term .ge. dsum_tol) then
      xlo = x
    else 
      xhi = x
    end if
  end do

  ew_coeff = x

  return

end subroutine compute_ew_coeff

end module mdin_ewald_dat_mod
