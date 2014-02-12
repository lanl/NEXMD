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
!             pmemd does not support this...
! vdwmeth = 0 for cutoff of vdw, 1 for analytic fix.
! verbose controls level of output. Look in force_info.


  integer, parameter    :: mdin_ewald_int_cnt = 8

  integer                       nfft1, nfft2, nfft3, bspl_order, netfrc, &
                                use_axis_opt, vdwmeth, verbose

  common / mdin_ewald_int /     nfft1, nfft2, nfft3, bspl_order, netfrc, &
                                use_axis_opt, vdwmeth, verbose

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
                                   eedmeth, ee_type, ischrgd, frameon, &
                                   chngmask, indmeth, maxiter, irstdip, &
                                   scaldip, use_cit, use_pme, column_fft

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
                        a, b, c, use_axis_opt, fft_grids_per_ang, column_fft

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

  nbtell = 0
  netfrc = -1                   ! -1 means "set to default value", once you
                                ! know imin.
  use_pme = 1
  vdwmeth = 1
  eedmeth = 1
  eedtbdns = 5000.d0
  ee_type = 1
  ischrgd = RETIRED_INPUT_OPTION
  frameon = UNSUPPORTED_INPUT_OPTION
  chngmask = UNSUPPORTED_INPUT_OPTION
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

  use_cit = RETIRED_INPUT_OPTION

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

#ifdef AMOEBA
  ! For now, until we understand all requirements, we turn axis optimization
  ! off under amoeba.
  ! BUGBUG - Need to fix this at some point...

  if (iamoeba .eq. 1) then
    if (use_axis_opt .gt. 0) then
      write(mdout, '(a,a,/)') info_hdr, &
        'Axis order optimization not used - not yet supported in amoeba.'
    end if
    use_axis_opt = 0
  end if
#endif /* AMOEBA */

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

  ! RCFFT needs an even dimension in the x direction (nfft1) ...

  if (nfft1 .eq. 0) &
    call compute_even_nfft(ew_box(1), fft_grids_per_ang, nfft1)
  if (nfft2 .eq. 0) &
    call compute_nfft(ew_box(2), fft_grids_per_ang, nfft2)
  if (nfft3 .eq. 0) &
    call compute_nfft(ew_box(3), fft_grids_per_ang, nfft3)

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

  call int_legal_range('vdwmeth', vdwmeth, 0, 1, inerr)

  if (eedmeth .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'eedmeth must be 1 for PMEMD!'
    inerr = 1
  end if

  if (use_pme .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'use_pme must be 1 for PMEMD!'
    inerr = 1
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

  if (frameon .ne. UNSUPPORTED_INPUT_OPTION) then
    write(mdout, '(a,a)') error_hdr, &
      'Extra points options (frameon) are not supported by PMEMD!'
    inerr = 1
  end if

  if (chngmask .ne. UNSUPPORTED_INPUT_OPTION) then
    write(mdout, '(a,a)') error_hdr, &
      'Extra points options (chngmask) are not supported by PMEMD!'
    inerr = 1
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

  call int_legal_range('use_axis_opt', use_axis_opt, 0, 1, inerr)

  if (use_cit .ne. RETIRED_INPUT_OPTION) then
    write(mdout, '(a,a)') info_hdr, &
      'The use_cit ewald option is deprecated and ignored.'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' always uses the CIT-based code paths.'
  end if

  if (column_fft .ne. 0) then
    write(mdout, '(a,a)') info_hdr, &
      'The column_fft ewald option is not supported and ignored.'
    write(mdout, '(a,a)') extra_line_hdr, &
      'This is a sander-specific optimization option.'
  end if

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

  return

9002 format (5x,'Box X =',f9.3,3x,'Box Y =',f9.3,3x,'Box Z =',f9.3)
9003 format (5x,'Alpha =',f9.3,3x,'Beta  =',f9.3,3x,'Gamma =',f9.3)
9004 format (5x,'NFFT1 =',i5  ,7x,'NFFT2 =',i5  ,7x,'NFFT3 =',i5)
9005 format (5x,'Interpolation order =',i5)
9006 format (5x,'Cutoff=',f9.3,3x,'Tol   =',e9.3)
9007 format (5x,'Ewald Coefficient =',f9.5)

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
                 0, mpi_comm_world, err_code_mpi)

  call mpi_bcast(dsum_tol, mdin_ewald_dbl_cnt, mpi_double_precision, &
                 0, mpi_comm_world, err_code_mpi)

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
