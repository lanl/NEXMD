!*******************************************************************************
!
! Module: mdin_ctrl_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module mdin_ctrl_dat_mod

use file_io_dat_mod

  implicit none

! The mdin ctrl input options are grouped in the order of input, with the
! common blocks containing valid options and the private data containing
! either deprecated options or sander 10 options not supported by pmemd.
! Some variables here actually are not &ctrl namelist parameters, but values
! derived from them that are needed during mdin &ctrl printout, etc.

  integer                       imin, nmropt, ntx, irest, ntrx, ntxo, &
                                ntpr, ntave, ntwr, iwrap, ntwx, ntwv, &
                                ntwe, ioutfm, ntwprt, ntf, ntb, nsnb, &
                                ipol, ibelly, ntr, maxcyc, ncyc, ntmin, &
                                nstlim, nscm, nrespa, ntt, ig, vrand, ntp, &
                                ntc, jfastw, ivcap, igb, alpb, rbornstat, &
                                gbsa, nrespai, ndfmin, jar, &
                                no_intermolecular_bonds, ene_avg_sampling, &
                                loadbal_verbose, ips, numexchg, exch_ntpr, &
                                csurften, ninterface,iamd,iamdlag  
  character(4)                  ihwtnm(2), iowtnm, iwtnm

  common / mdin_ctrl_int /      imin, nmropt, ntx, irest, ntrx, ntxo, &     !  6
                                ntpr, ntave, ntwr, iwrap, ntwx, ntwv, &     ! 12
                                ntwe, ioutfm, ntwprt, ntf, ntb, nsnb, &     ! 18
                                ipol, ibelly, ntr, maxcyc, ncyc, ntmin, &   ! 24
                                nstlim, nscm, nrespa, ntt, ig, vrand, &     ! 30
                                ntp, ntc, jfastw, ivcap, igb, alpb, &       ! 36
                                rbornstat, gbsa, nrespai, ndfmin, jar, &    ! 41
                                no_intermolecular_bonds, ene_avg_sampling, &! 43
                                loadbal_verbose, numexchg, exch_ntpr, &     ! 46
                                ihwtnm, iowtnm, iwtnm, ips, &               ! 51
                                csurften, ninterface,iamd, iamdlag        ! 55

  integer, parameter    :: mdin_ctrl_int_cnt = 55

  save  :: / mdin_ctrl_int /

  ! Note - ndfmin is not intended to be settable by the user.  We permit it to
  !        show up in the namelist but issue a warning and ignore it.

  ! New cut variables for pmemd.  These are set based on actual input of the
  ! values, or alternatively on the value of the 'cut' input variable, or based
  ! on default.  They are what actually gets used to determine cutoffs for pme:
  ! es_cutoff - electrostatic cutoff.
  ! vdw_cutoff - vdw cutoff.

  ! no_intermolecular_bonds -
  !
  ! New pmemd variable controlling intermolecular bonding.  If 1, any
  ! molecules joined by a covalent bond are fused; if 0 then the old
  ! behaviour (use prmtop molecule definitions) pertains.  The default is 1;
  ! a value of 0 is not supported with forcefields using extra points.

  ! ene_avg_sampling -
  !
  ! Number of steps between energy samples used in averages.  If not specified
  ! then ntpr is used (default).  To match the behaviour of earlier releases,
  ! this variable should be set to 1.  This variable is only used for MD,
  ! not minimization and will also effectively be turned off if ntave is in use
  ! or RESPA is in use.

  double precision              dielc, es_cutoff, vdw_cutoff, &
                                dx0, drms, t, dt, temp0, tempi, &
                                tautp, gamma_ln, vlimit, pres0, comp, taup, &
                                tol, fcap, pencut, intdiel, extdiel, saltcon, &
                                rgbmax, offset, surften, cut_inner, gb_cutoff, &
                                gb_alpha, gb_beta, gb_gamma, gb_fs_max, &
                                gb_kappa, gb_neckscale, arad,ips_cut,dvbips, &
                                restraint_wt, gb_alpha_h, gb_beta_h, gb_gamma_h, &
                                gb_alpha_c, gb_beta_c, gb_gamma_c, gb_alpha_n, &
                                gb_beta_n, gb_gamma_n, gb_alpha_os, gb_beta_os, &
                                gb_gamma_os, gb_alpha_p, gb_beta_p, gb_gamma_p, &
                                screen_h, screen_c, screen_n, screen_o, &
                                screen_s, screen_p, gamma_ten,EthreshD,alphaD,EthreshP,alphaP

  common / mdin_ctrl_dbl /      dielc, es_cutoff, vdw_cutoff, &                  !3
                                dx0, drms, t, dt, temp0, tempi, &                !9
                                tautp, gamma_ln, vlimit, pres0, comp, taup, &    !15
                                tol, fcap, pencut, intdiel, extdiel, saltcon, &  !21
                                rgbmax, offset, surften, cut_inner, gb_cutoff, & !26
                                gb_alpha, gb_beta, gb_gamma, gb_fs_max, &        !30
                                gb_kappa, gb_neckscale, arad, ips_cut, &         !34
                                dvbips, restraint_wt, gb_alpha_h, gb_beta_h, &   !38
                                gb_gamma_h, gb_alpha_c, gb_beta_c, gb_gamma_c, & !42
                                gb_alpha_n, gb_beta_n, gb_gamma_n, gb_alpha_os,& !46
                                gb_beta_os, gb_gamma_os, gb_alpha_p, gb_beta_p,& !50
                                gb_gamma_p, screen_h, screen_c, screen_n, &      !54
                                screen_o, screen_s, screen_p, gamma_ten,EthreshD, & !59
                                alphaD,EthreshP, alphaP                           !62

  integer, parameter    :: mdin_ctrl_dbl_cnt = 62


  save  :: /mdin_ctrl_dbl /

  ! Restraint and belly masks

  character(256), public   :: restraintmask, bellymask

  ! Note - gb_alpha, gb_beta, gb_gamma, gb_fs_max, gb_kappa and gb_neckscale
  !        would be better factored elsewhere, but as usual the sander order
  !        of operations makes this difficult...

  ! Master-only PMEMD-specific control variables:

  integer, save :: mdinfo_flush_interval        ! controls mdinfo flushing.
  integer, save :: mdout_flush_interval         ! controls mdout flushing.

  ! Undocumented control variables, for dev use:

  integer, save :: dbg_atom_redistribution      ! force more frequent atom
                                                ! redistribution; actually get
                                                ! better testing now in assoc
                                                ! with fft slab redistribution.

  ! Logical variables used to indicate the type of potential to use.
  ! Currently the options include PME and Generalized Born.  These are
  ! set in the broadcast module for mpi, and are based on the igb values.
  ! Direct use of igb is undesirable because it designates multiple gb
  ! implementations.

  logical, save :: using_gb_potential
  logical, save :: using_pme_potential
  logical, save :: TEIPS,TVIPS

  !no_ntt3_sync - This tells the code to not synchronize the random number stream between
  !               MPI tasks when running with ntt=3. This improves scaling but means that
  !               the random number streams and therefore results are dependent on the 
  !               number of MPI tasks. It is therefore only turned on when ig=-1.

  logical, save :: no_ntt3_sync

  ! The t variable gets stomped on by inpcrd code before printout so...

  double precision, save, private     :: original_t

  ! Namelist entries that are ambiguous and replaced with more meaningful
  ! names in the code:

  double precision, save, private       :: cut  ! local use only, for clarity.

  ! Namelist entries used in water name processing, but not exported:

  character(4), save, private           :: hwtnm1, hwtnm2, owtnm, watnam

  ! Namelist entries that are used for amber 10 or backward compatibility,
  ! but that are not actually used.  These are here so we can provide
  ! meaningful error messages to the user instead of just falling over
  ! in namelist processing.  These are visible in this module only, in order
  ! that they can be seen for validation/printing.

  integer, private              :: idecomp, itgtmd, ifqnt, icnstph, ntcnstph, &
                                   iscale, noeskp, ipnlty, mxsub, icfe, &
                                   klambda, n3b, nion, iesp, npscal, &
                                   lastist, lastrst, isgld, isgsta, isgend, &
                                   ievb, ifcr

  integer, private              :: amber7_compat, amber8_mdout_format

  double precision, private     :: tgtrmsd, tgtmdfrc, solvph, temp0les, &  
                                   scalm, tausw, clambda, timlim, dtemp, &
                                   dxm, heat, rdt, tsgavg, tempsg, sgft


  character(256), private       :: tgtrmsmask, tgtfitmask

  ! The active namelist:

  private       :: cntrl

  namelist /cntrl/      imin, nmropt, ntx, irest, ntrx, ntxo, &
                        ntpr, ntave, ntwr, iwrap, ntwx, ntwv, &
                        ntwe, ioutfm, ntwprt, idecomp, ntf, ntb, &
                        dielc, cut, nsnb, ipol, igb, &
                        intdiel, extdiel, saltcon, rgbmax, rbornstat, &
                        offset, gbsa, surften, rdt, nrespai, cut_inner, &
                        ibelly, ntr, restraint_wt, restraintmask, bellymask, &
                        itgtmd, tgtrmsd, tgtmdfrc, tgtrmsmask, tgtfitmask, &
                        maxcyc, ncyc, ntmin, dx0, drms, nstlim, nscm, t, dt, &
                        nrespa, ntt, temp0, temp0les, tempi, ig, &
                        tautp, gamma_ln, vrand, vlimit, ntp, &
                        pres0, comp, taup, ntc, tol, jfastw, &
                        hwtnm1, hwtnm2, owtnm, watnam, &
                        ivcap, fcap, iscale, noeskp, ipnlty, mxsub, &
                        scalm, pencut, tausw, icfe, clambda, klambda, &
                        timlim, ndfmin, jar, n3b, nion, iesp, npscal, &
                        lastist, lastrst, no_intermolecular_bonds, &
                        ene_avg_sampling, &
                        amber7_compat, amber8_mdout_format, &
                        mdinfo_flush_interval, &
                        mdout_flush_interval, &
                        dbg_atom_redistribution, &
                        loadbal_verbose, &
                        es_cutoff, vdw_cutoff, &
                        dtemp, dxm, heat, alpb, arad, &
                        ifqnt, icnstph, ntcnstph, solvph, &
                        isgld, tsgavg, tempsg, sgft, isgsta, isgend, &
                        ievb, ips, dvbips, numexchg, exch_ntpr, &
                        csurften, gamma_ten, ninterface, ifcr, iamd,iamdlag, &
                        EthreshD,alphaD,EthreshP,alphaP

contains

!*******************************************************************************
!
! Subroutine:  init_mdin_ctrl_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine init_mdin_ctrl_dat

  use file_io_mod
  use gbl_constants_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Local variables:

  integer               ifind, sec_temp
  character(4)          watdef(4)

! Define default water residue name and the names of water oxygen & hydrogens

  data watdef/'WAT ', 'O   ', 'H1  ', 'H2  '/

  hwtnm1 = '    '
  hwtnm2 = '    '
  owtnm =  '    '
  watnam = '    '

!AMD initial values
!iamd=0 no boost is used, 1 boost on the total energy, 2 boost on the dihedrals, 3 boost on dihedrals and total energy
  iamd = 0
  iamdlag = 0 !frequency of boosting in steps
  EthreshD = 0.d0
  alphaD = 0.d0
  EthreshP = 0.d0
  alphaP = 0.d0
  
!IPS initial values
  ips = 0    ! no isotropic periodic sum
  dvbips = 1.0d-8 ! Volume change tolerance. aips done when change > dvbips

  imin = 0
  nmropt = 0
  jar = 0                       ! Jarzynski Relationship - is supported.
  no_intermolecular_bonds = 1
  ene_avg_sampling = -1        ! -1 means "set to default value".
  ntx = 1
  irest = 0
  ntrx = 1
  ntxo = 1
  ntpr = 50
  ntave = 0
  ntwr = 500

  ! pmemd adjusts ntwr upward for high processor count...

#ifdef MPI
  if (numtasks .gt. 10) ntwr = 50 * numtasks
#endif

  iwrap = 0
  ntwx = 0
  ntwv = 0
  ntwe = 0
  ioutfm = 0
  ntwprt = 0
  idecomp = 0           ! not supported...
  ntf = 1
  ntb = NO_INPUT_VALUE  ! default depends on igb and ntp
  dielc = 1.d0
  cut = 0.d0            ! indicates value not set...
  es_cutoff = 0.d0      ! indicates value not set...
  vdw_cutoff = 0.d0     ! indicates value not set...
  nsnb = 25
  ipol = 0
  igb = 0
  alpb = 0              ! used for Analytical Linearized Poisson-Boltzmann
  arad = 15.d0          ! used for Analytical Linearized Poisson-Boltzmann
  intdiel = 1.d0        ! used for generalized Born
  extdiel = 78.5d0      ! used for generalized Born
  saltcon = 0.d0        ! used for generalized Born
  rgbmax = 25.d0        ! used for generalized Born
  rbornstat = 0         ! used for generalized Born
  offset = 0.09d0       ! used for generalized Born. Overwritten for igb=8
  gbsa = 0              ! used for generalized Born
  surften = 0.005d0     ! used for generalized Born
  nrespai = 1           ! used for generalized Born
  cut_inner = 8.d0      ! used for generalized Born

! Generalized Born variables, perhaps better factored elsewhere, but...

  gb_cutoff = 0.d0
  gb_alpha = 0.d0
  gb_beta = 0.d0
  gb_gamma = 0.d0
  gb_fs_max = 0.d0
  gb_kappa = 0.d0
  gb_neckscale = 1.d0 / 3.d0
! Used for igb = 8, the defaults are set below
  gb_alpha_h = 0.d0
  gb_beta_h = 0.d0
  gb_gamma_h = 0.d0
  gb_alpha_c = 0.d0
  gb_beta_c = 0.d0
  gb_gamma_c = 0.d0
  gb_alpha_n = 0.d0
  gb_beta_n = 0.d0
  gb_gamma_n = 0.d0
  gb_alpha_os = 0.d0
  gb_beta_os = 0.d0
  gb_gamma_os = 0.d0
  gb_alpha_p = 0.d0
  gb_beta_p = 0.d0
  gb_gamma_p = 0.d0
  screen_h = 0.d0
  screen_c = 0.d0
  screen_n = 0.d0
  screen_o = 0.d0
  screen_s = 0.d0
  screen_p = 0.d0

  ibelly = 0
  ntr = 0
  restraint_wt = 0.d0
  restraintmask = ''
  bellymask = ''
  itgtmd = 0            ! not supported...
  tgtrmsd = 0.d0        ! not actually used...
  tgtmdfrc = 0.d0       ! not actually used...
  tgtrmsmask = ''       ! not actually used...
  tgtfitmask = ''       ! not actually used...
  maxcyc = 1
  ncyc = 10
  ntmin = 1
  dx0 = 0.01d0
  drms = 1.0d-4
  nstlim = 1
  nscm = 1000
  t = 0.d0
  dt = 0.001d0
  nrespa = 1
  ntt = 0
  temp0 = 300.d0
  temp0les = -1.d0      ! not actually used...
  rdt = 0.0             ! not actually used...
  tempi = 0.d0
  ig = 71277; no_ntt3_sync = .false.
  tautp = 1.d0
  gamma_ln = 0.d0
  vrand = 1000
#ifdef CUDA
  vlimit = -1.0d0 !Vlimit not supported on GPUs
#else
  vlimit = 20.0d0
#endif
  ntp = 0
  pres0 = 1.d0
  taup = 1.d0
  comp = 44.6d0
  ntc = 1
  tol = 0.00001d0
  jfastw = 0
  ivcap = 0
  fcap = 1.5d0
  numexchg = 0
  exch_ntpr = 1
  iscale = 0            ! not supported...
  noeskp = 1            ! not supported...
  ipnlty = 1            ! not supported...
  mxsub = 1             ! not supported...
  scalm = 100.d0        ! not supported...
  pencut = 0.1d0
  tausw = 0.1d0         ! not supported...
  icfe = 0              ! not supported...
  clambda = 0.d0        ! not supported...
  klambda = 1           ! not supported...
  ifqnt = 0             ! not supported...
  icnstph = 0           ! not supported...
  ntcnstph = 10         ! not supported...
  solvph = 7.d0         ! not supported...

  ! Self-guided Langevin dynamics (not currently supported)

  isgld = 0             ! not supported...
  isgsta = 1            ! not supported...
  isgend = 0            ! not supported...
  tsgavg = 0.2d0        ! not supported...
  sgft = 0.d0           ! not supported...
  tempsg = 1.d0         ! not supported...

  ! Empirical Valence Bond dynamics (not currently supported)

  ievb = 0
  ndfmin = RETIRED_INPUT_OPTION
  dtemp = RETIRED_INPUT_OPTION
  dxm = RETIRED_INPUT_OPTION
  heat = RETIRED_INPUT_OPTION

  timlim =  999999.d0   ! this option has not really done anything since sander6

  n3b = 0               ! not supported...
  nion = 0              ! not supported...
  iesp = 0              ! not supported...
  ifcr = 0              ! not supported...
  npscal = 1            ! only supported value.
  lastist = 0           ! ignored
  lastrst = 0           ! ignored

! Constant Surface Tension
  csurften = 0      !constant surface tension off (valid options are 0,1,2,3)
  gamma_ten = 0.0d0 !0.0 dyne/cm - default used in charmm. Ignored for csurften=0
  ninterface = 2   !Number of interfaces in the surface tension (Must be greater than 2)

! PMEMD-specific options:

  mdinfo_flush_interval = 60    ! Flush mdinfo every 60 seconds by default
  mdout_flush_interval = 300    ! Flush mdout every 300 seconds by default

! Retired options that we will complain about, but not fail on.

  ! PMEMD 11 always runs in Amber 11 compatibility mode.

  amber7_compat = RETIRED_INPUT_OPTION

  ! Amber 11 mdout format used in pmemd 3.1. Now is default and only option.

  amber8_mdout_format = RETIRED_INPUT_OPTION

  dbg_atom_redistribution = 0   ! Don't do atom redivision at every list
                                ! build (this is really only good for debugging
                                ! the atom redivision code).
  
  loadbal_verbose = 0           ! Loadbalancing verbosity.  Default is minimal
                                ! loadbalancing info in log file.  Values of
                                ! 1,2,3... increase amount of info dumped.
                                ! This variable only has an effect for
                                ! parallel runs.

! Read input in namelist format:

  rewind(mdin)                  ! Insurance against maintenance mods.

  call nmlsrc('cntrl', mdin, ifind)

  if (ifind .ne. 0) then        ! Namelist found. Read it:
    read(mdin, nml = cntrl)
  else                          ! No namelist was present, 
    write(mdout, '(a,a)') error_hdr, 'Could not find cntrl namelist'
    call mexit(6, 1)
  end if

  ! Set the proper value of ntb based on igb and ntp
  if (ntb .eq. NO_INPUT_VALUE) then
    if (ntp .gt. 0) then
      ntb = 2
    else if (igb .gt. 0) then
      ntb = 0
    else
      ntb = 1
    end if
  end if

  ! Set the definition of the water molecule. The default def is in watdef(4).
  ! Reset definitions if there are new definitions.

  read(watdef(1), '(A4)') iwtnm
  read(watdef(2), '(A4)') iowtnm
  read(watdef(3), '(A4)') ihwtnm(1)
  read(watdef(4), '(A4)') ihwtnm(2)

  if (watnam .ne. '    ') read(watnam, '(A4)') iwtnm
  if (owtnm  .ne. '    ') read(owtnm,  '(A4)') iowtnm
  if (hwtnm1 .ne. '    ') read(hwtnm1, '(A4)') ihwtnm(1)
  if (hwtnm2 .ne. '    ') read(hwtnm2, '(A4)') ihwtnm(2)

  ! The variable pair es_cutoff, vdw_cutoff currently only applies to pme
  ! simulations.  We use a variable gb_cutoff to capture the generalized Born
  ! cutoff (specified as cut in &ctrl).

  if (igb .eq. 0) then

    ! If none of the cutoff variables were set, set the default:

    if (cut .eq. 0.d0 .and. es_cutoff .eq. 0.d0 .and. vdw_cutoff .eq. 0.d0) then
      es_cutoff = 8.d0
      vdw_cutoff = 8.d0
    end if

    ! If only cut was set, propagate it to es_cutoff and vdw_cutoff.

    if (cut .ne. 0.d0 .and. es_cutoff .eq. 0.d0 .and. vdw_cutoff .eq. 0.d0) then
      es_cutoff = cut
      vdw_cutoff = cut
    end if

  else if (igb .eq. 1 .or. igb .eq. 2 .or. igb .eq. 5 .or. &
           igb .eq. 7 .or. igb .eq. 8) then

    ! We will catch incorrect input of es_cutoff/vdw_cutoff later.  Here we just
    ! set gb_cutoff to cut, or the default value if cut was not set.

    if (cut .eq. 0.d0) then
      gb_cutoff = 9999.d0
    else
      gb_cutoff = cut
    end if

    if (igb .eq. 1) then
      gb_alpha = 1.d0
      gb_beta = 0.d0
      gb_gamma = 0.d0
    else if (igb .eq. 2) then
      ! Use our best guesses for Onufriev/Case GB  (GB^OBC I):
      gb_alpha = 0.8d0
      gb_beta = 0.d0
      gb_gamma = 2.909125d0
    else if (igb .eq. 5) then
      ! Use our second best guesses for Onufriev/Case GB (GB^OBC II):
      gb_alpha = 1.d0
      gb_beta = 0.8d0
      gb_gamma = 4.851d0
    else if (igb .eq. 7) then
      ! Use parameters for Mongan et al. CFA GBNECK:
      gb_alpha = 1.09511284d0
      gb_beta = 1.90792938d0
      gb_gamma = 2.50798245d0
      gb_neckscale = 0.361825d0
    else if (igb .eq. 8) then
      gb_alpha_h   = 0.788440d0
      gb_beta_h    = 0.798699d0
      gb_gamma_h   = 0.437334d0
      gb_alpha_c   = 0.733756d0
      gb_beta_c    = 0.506378d0
      gb_gamma_c   = 0.205844d0
      gb_alpha_n   = 0.503364d0
      gb_beta_n    = 0.316828d0
      gb_gamma_n   = 0.192915d0
      gb_alpha_os  = 0.867814d0
      gb_beta_os   = 0.876635d0
      gb_gamma_os  = 0.387882d0
      gb_alpha_p   = 1.0d0
      gb_beta_p    = 0.8d0
      gb_gamma_p   = 4.851d0  
      gb_neckscale = 0.826836d0
      screen_h = 1.425952d0
      screen_c = 1.058554d0
      screen_n = 0.733599d0
      screen_o = 1.061039d0
      screen_s = -0.703469d0
      screen_p = 0.5d0
      offset  = 0.195141d0
    end if

    ! Set gb_kappa as long as saltcon is .ge. 0.d0 (and catch bug below if not).

    if (saltcon .ge. 0.d0) then

      ! Get Debye-Huckel kappa (A**-1) from salt concentration (M), assuming:
      !   T = 298.15, epsext=78.5,
      gb_kappa = sqrt(0.10806d0 * saltcon)

      ! Scale kappa by 0.73 to account(?) for lack of ion exclusions:

      gb_kappa = 0.73d0 * gb_kappa

    end if

  end if

  ! If ene_avg_sampling was not set, set to default:

  if (ene_avg_sampling .eq. -1) then
    if (ntpr .gt. 0) then
      ene_avg_sampling = ntpr
    else
      ene_avg_sampling = 1
    end if
  end if

  ! If running averages have been invoked, 
  ! or we are doing minimization or we are using nrespa .ne. 1,
  ! we silently override any ene_avg_sampling setting.

  if (imin .eq. 1 .or. ntave .gt. 0) ene_avg_sampling = 1

  if (nrespa .ne. 1) ene_avg_sampling = nrespa

! Here we reset some possibly errant variables without killing the run.

  if (dielc .le. 0.0d0) dielc = 1.0d0
  if (taup .le. 0.d0) taup = 0.2d0
  if (tautp <= 0.d0 ) tautp = 0.2d0
  if (nstlim .lt. 0) nstlim = 0
  if (ntpr .eq. 0) ntpr = nstlim + 1

! If Jarzynski Relationship is being used, force nmropt on.

  if (jar .ne. 0) nmropt = 1

  if (ntxo .eq. 0) then
    ntxo = 1
    write(mdout, '(a,a)') warn_hdr, &
    'NTXO was been reset to 1 (formatted restrt output); because binary'
    write(mdout, '(a,a)') extra_line_hdr, &
    'restrt files are not portable, writing them is no longer supported.'
  end if

! Save original value of t for later printing.

  original_t = t

! Set the potential flags:

  if (igb .eq. 0) then
    using_pme_potential = .true.
    using_gb_potential = .false.
  else if (igb .eq. 1 .or. igb .eq. 2 .or. igb .eq. 5 .or. &
           igb .eq. 7 .or. igb .eq. 8) then
    using_pme_potential = .false.
    using_gb_potential = .true.
  else                                  ! in a world of trouble...
    using_pme_potential = .false.
    using_gb_potential = .false.
  end if

! Set IPS potential flags

   !--------------------------------------------------------------------
   ! parameters for IPS and for SGLD:
   ! To enhance performance in pmemd we have decided to allow only 
   ! IPS calculations with ips=1 so if ips=1 teips and tvips are set to TRUE 
   ! ips=1  3D IPS for electrostatic and Lennard-Jones potentials
   !--------------------------------------------------------------------
   
   teips=.false.
   tvips=.false.
   if( (ips -1) == 0 )tvips =.true.
   if( (ips -1) == 0 )teips =.true.
   if( teips ) then
    using_pme_potential = .true.
    using_gb_potential = .false.
   end if
   if( tvips ) then
    using_pme_potential = .true.
    using_gb_potential = .false.
   end if
 
  !Support for random seed using time of day in microsecs
  if( ig==-1 ) then
    !Turn on NO_NTT3_SYNC when ig=-1. This means the code no
    !longer synchronized the random numbers between streams when
    !running in parallel giving better scaling.
    no_ntt3_sync = .true.
#ifdef MPI
    write (mdout, '(a)') "Note: ig = -1. Setting random seed based on wallclock &
                               &time in microseconds"
    write (mdout, '(a)') "      and disabling the synchronization of random &
                               &numbers between tasks"
    write (mdout, '(a)') "      to improve performance."
#else 
    write (mdout, '(a)') "Note: ig = -1. Setting random seed based on wallclock &
                                &time in microseconds."
#endif
    call get_wall_time(sec_temp,ig)
  end if

  return

end subroutine init_mdin_ctrl_dat

!*******************************************************************************
!
! Subroutine:  validate_mdin_ctrl_dat
!
! Description: <TBS>
!              
!*******************************************************************************
#ifdef MPI
subroutine validate_mdin_ctrl_dat(remd_method)
#else
subroutine validate_mdin_ctrl_dat
#endif

  use gbl_constants_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

#ifdef MPI
! Passed variables

  integer, intent(in) :: remd_method
#endif

! Local variables

  integer       :: inerr

! Check on bogus data and unsupported options:

  inerr = 0

  if (imin .ne. 0 .and. imin .ne. 1) then
    write(mdout, '(a,a)') error_hdr, &
      'imin must be set to 0 for (molecular dynamics) or 1 for (minimization)!'
    if (imin .eq. 5) then
      write(mdout, '(a,a,a)') error_hdr, prog_name, &
        ' does not support imin == 5 (trajectory analysis).'
      write(mdout, '(a)') use_sander
    end if
    inerr = 1
  end if

  if (nmropt .ne. 0 .and. nmropt .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'nmropt must be set to 0 or 1!'
    if (nmropt .eq. 2) then
      write(mdout, '(a,a,a)') error_hdr, prog_name, &
        ' does not support nmropt == 2 (NOESY vol or chem shift restraints).'
      write(mdout, '(a)') use_sander
    end if
    inerr = 1
  end if

  if (ntx .lt. 1 .or. ntx .gt. 7 .or. ntx .eq. 3) then
    write(mdout, '(a,a)') error_hdr, 'ntx must be set to 1, 2, 4, 5, 6, or 7!'
    inerr = 1
  end if

  if (irest .ne. 0 .and. irest .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'irest must be set to 0 or 1!'
    inerr = 1
  end if

  if  ((ntx .eq. 1 .or. ntx .eq. 2) .and. irest .ne. 0)  then
    write(mdout, '(a,a)') error_hdr, 'ntx and irest are inconsistent!'
    inerr = 1
  end if

  if (ntrx .ne. 0 .and. ntrx .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'ntrx must be set to 0, or 1!'
    inerr = 1
  end if

  if (ntxo .ne. 0 .and. ntxo .ne. 1 .and. ntxo .ne. 2) then
    write(mdout, '(a,a)') error_hdr, 'ntxo must be set to 0, 1, or 2!'
    inerr = 1
  end if

  if (ntpr .lt. 0) then
    write(mdout, '(a,a)') error_hdr, 'ntpr must be >= 0!'
    inerr = 1
  end if

  if (ntave .lt. 0) then
    write(mdout, '(a,a)') error_hdr, 'ntave must be >= 0!'
    inerr = 1
  end if

  if (iwrap .ne. 0 .and. iwrap .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'iwrap must be set to 0 or 1!'
    inerr = 1
  end if

  if (iwrap .eq. 1 .and. ntb .eq. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'iwrap = 1 must be used with a periodic box!'
    inerr = 1
  end if

  if (ntwx .lt. 0) then
    write(mdout, '(a,a)') error_hdr, 'ntwx must be >= 0!'
    inerr = 1
  end if

  if (ntwv .lt. -2) then
    write(mdout, '(a,a)') error_hdr, 'ntwv must be >= -1!'
    inerr = 1
  end if

  ! We only allow velocity output in the mdcrd file at the same frequency of
  ! the coords (designated by ntwv .eq. -1) if we are writing an netCDF
  ! binary file.

  if (ntwv .eq. -1 .and. ioutfm .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'ntwv may be -1 only if ioutfm == 1!'
    inerr = 1
  end if

  if (ntwv .eq. -1 .and. ntwx .eq. 0) then
    write(mdout, '(a,a)') error_hdr, 'ntwv may be -1 only if ntwx > 0!'
    inerr = 1
  end if

  if (ntwe .lt. 0) then
    write(mdout, '(a,a)') error_hdr, 'ntwe must be >= 0!'
    inerr = 1
  end if

  if (ioutfm .lt. 0 .or. ioutfm .gt. 1) then
    write(mdout, '(a,a)') error_hdr, 'ioutfm must be 0 or 1!'
    inerr = 1
  end if

  if (ntwprt .lt. 0) then
    write(mdout, '(a,a)') error_hdr, 'ntwprt must be >= 0!'
    inerr = 1
  end if

  if (idecomp .ne. 0) then
    write(mdout, '(a,a,a)') error_hdr, prog_name, &
      ' does not support idecomp != 0!'
    inerr = 1
  end if

  if (ntf .lt. 1 .or. ntf .gt. 8) then
    write(mdout, '(a,a)') error_hdr, 'ntf must be in the range 1..8!'
    inerr = 1
  end if

  if (ntb .ne. 0 .and. ntb .ne. 1 .and. ntb .ne. 2) then
    write(mdout, '(a,a)') error_hdr, 'ntb must be 0, 1 or 2!'
    inerr = 1
  end if

  if (nsnb .ne. 25) then
    write(mdout, '(a,a)') info_hdr, &
      'The nsnb ctrl option does not affect nonbonded list update frequency.'
    write(mdout, '(a,a)') extra_line_hdr, &
      'It does affect steepest descent minimization freq if ntmin == 0'
  end if

  if (ipol .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'ipol must == 0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support this type of polarizable force field calculations.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if


  if (jar .ne. 0 .and. jar .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'jar must be 0 or 1!'
  end if

  if (no_intermolecular_bonds .ne. 0 .and. no_intermolecular_bonds .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'no_intermolecular_bonds must be 0 or 1!'
  end if

  if (ene_avg_sampling .le. 0) then
    write(mdout, '(a,a)') error_hdr, 'ene_avg_sampling must be .gt. 0!'
  end if

  if (igb .ne. 0 .and. &
      igb .ne. 1 .and. &
      igb .ne. 2 .and. &
      igb .ne. 5 .and. &
      igb .ne. 7 .and. &
      igb .ne. 8) then
    if (igb .eq. 10) then
      write(mdout, '(a,a,a)') error_hdr, prog_name, &
        ' does not support Poisson-Boltzmann simulations (igb .eq. 10)!'
      write(mdout, '(a)') use_sander
    else
      write(mdout, '(a,a)') error_hdr, 'igb must be 0,1,2,5,7, or 8!'
    end if
    inerr = 1
  end if

  if (alpb .ne. 0 .and. alpb .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'alpb must be 0 or 1!'
  end if

  if (alpb .eq. 1 .and. .not. using_gb_potential) then
    write(mdout, '(a,a)') error_hdr, 'igb must be 1,2,5,7, or 8 if alpb == 1!'
  end if

  if (saltcon .lt. 0.d0) then
    write(mdout, '(a,a)') error_hdr, 'saltcon must be a positive value!'
    inerr = 1
  end if

  if (rgbmax .lt. 5.d0 * gb_fs_max) then
    write(mdout, '(a,a,f8.2)') error_hdr, 'rgbmax must be at least ', &
      5.d0 * gb_fs_max
    inerr = 1
  end if

  if (rbornstat .ne. 0 .and. rbornstat .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'rbornstat must be 0 or 1!'
    inerr = 1
  end if

  if (intdiel .ne. 1.d0) then
    write(mdout, '(a,a)') error_hdr, 'Non-default values of intdiel are not supported!'
    inerr = 1
  end if

  if (gbsa .ne. 0 .and. gbsa .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'gbsa must be set to 0 or 1!'
    inerr = 1
  else ! this is a valid choice of gbsa... but only for GB sims.
    if ( gbsa .eq. 1 .and. .not. using_gb_potential) then
      write(mdout, '(a,a)') error_hdr, 'gbsa can only be used with a Generalized &
        &Born solvent model.'
      inerr = 1
    end if
  end if

  ! This is a valid option in sander -- print out a warning here

  if (gbsa .eq. 2) then
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support gbsa == 2'
    write(mdout, '(a)') use_sander
  end if

  if (ibelly .ne. 0 .and. ibelly .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'ibelly must be set to 0 or 1!'
    inerr = 1
  end if

  if (ibelly == 1 .and. ntt == 3) then
    write(mdout, '(a,a)') error_hdr, 'ibelly = 1 with ntt = 3 is not a valid option.'
    inerr = 1
  end if

  if (ntr .ne. 0 .and. ntr .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'ntr must be set to 0 or 1!'
    inerr = 1
  end if

  if (restraint_wt .lt. 0) then
    write(mdout, '(a,a)') error_hdr, 'restraint_wt must be >= 0!'
    inerr = 1
  end if

  if (itgtmd .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'itgtmd must == 0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support targeted MD.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (tgtrmsd .ne. 0.d0) then
    write(mdout, '(a,a)') error_hdr, 'tgtrmsd is only used if itgtmd != 0!'
    inerr = 1
  end if

  if (tgtmdfrc .ne. 0.d0) then
    write(mdout, '(a,a)') error_hdr, 'tgtmdfrc is only used if itgtmd != 0!'
    inerr = 1
  end if

  if (tgtrmsmask .ne. '') then
    write(mdout, '(a,a)') error_hdr, 'tgtrmsmask is only used if itgtmd != 0!'
    inerr = 1
  end if

  if (tgtfitmask .ne. '') then
    write(mdout, '(a,a)') error_hdr, 'tgtfitmask is only used if itgtmd != 0!'
    inerr = 1
  end if

  if (ifqnt .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'ifqnt must == 0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support QM/MM calculations.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (icnstph .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'icnstph must == 0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support constant pH calculations.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (ntcnstph .ne. 10) then
    write(mdout, '(a,a)') error_hdr, 'ntcnstph is only used if icnstph != 0!'
    inerr = 1
  end if

  if (solvph .ne. 7.d0) then
    write(mdout, '(a,a)') error_hdr, 'solvph is only used if icnstph != 0!'
    inerr = 1
  end if

  if (isgld .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'isgld must == 0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support Self-guided Langevin dynamics.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (isgsta .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'isgsta is only used if isgld != 0!'
    inerr = 1
  end if

  if (isgend .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'isgend is only used if isgld != 0!'
    inerr = 1
  end if

  if (tsgavg .ne. 0.2d0) then
    write(mdout, '(a,a)') error_hdr, 'tsgavg is only used if isgld != 0!'
    inerr = 1
  end if

  if (sgft .ne. 0.d0) then
    write(mdout, '(a,a)') error_hdr, 'sgft is only used if isgld != 0!'
    inerr = 1
  end if

  if (tempsg .ne. 1.d0) then
    write(mdout, '(a,a)') error_hdr, 'tempsg is only used if isgld != 0!'
    inerr = 1
  end if

  if (ievb .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'ievb must == 0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support Empirical Valence Bond dynamics.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (ntmin .lt. 0 .or. ntmin .gt. 2) then
    write(mdout, '(a,a)') error_hdr, 'ntmin must be set to 0, 1 or 2!'
    if (ntmin .eq. 3) then
      write(mdout, '(a,a,a)') error_hdr, prog_name, &
        ' does not support ntmin == 3 (LMOD XMIN minimization method).'
      write(mdout, '(a)') use_sander
    else if (ntmin .eq. 4) then
      write(mdout, '(a,a,a)') error_hdr, prog_name, &
        ' does not support ntmin == 4 (LMOD LMOD minimization method).'
      write(mdout, '(a)') use_sander
    end if
    inerr = 1
  end if

  if (nrespa .lt. 1) then
    write(mdout, '(a,a)') error_hdr, 'nrespa must be >= 1!'
    inerr = 1
  end if

  if (nrespai .lt. 1) then
    write(mdout, '(a,a)') error_hdr, 'nrespai must be >= 1!'
    inerr = 1
  end if

  if (ntt .lt. 0 .or. ntt .gt. 3) then
    write(mdout, '(a,a)') error_hdr, 'ntt must be set to 0, 1, 2, or 3!'
    inerr = 1
  end if

  if (temp0les .ne. -1.d0) then
    write(mdout, '(a,a)') error_hdr, 'temp0les must == -1.d0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support LES.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (rdt .ne. 0.d0) then
    write(mdout, '(a,a)') error_hdr, 'rdt must == 0.d0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support LES.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (ntp .lt. 0 .or. ntp .gt. 3) then
    write(mdout, '(a,a)') error_hdr, 'ntp must be set to 0, 1, 2 or 3!'
    inerr = 1
  end if

  if (ntp == 3 .and. csurften < 1) then
    write(mdout, '(a,a)') error_hdr, 'ntp=3 may only be used with csurften>0!'
    inerr = 1
  end if

  if (ntc .lt. 1 .or. ntc .gt. 3) then
    write(mdout, '(a,a)') error_hdr, 'ntc must be set to 1, 2 or 3!'
    inerr = 1
  end if

  if (jfastw .lt. 0 .or. jfastw .gt. 4) then
    ! This is very sloppy, but sander-compliant.  Basically, values 0, 1, 2
    ! and 3 all specify to use fast water shake, with names as provided by
    ! other variables.  If jfastw .eq. 4, then don't use fast water shake.
    write(mdout, '(a,a)') error_hdr, 'jfastw must be between 0 and 4!'
    inerr = 1
  end if

  if (ivcap .ne. 0 .and. ivcap .ne. 2) then
    write(mdout, '(a,a)') error_hdr, 'ivcap must be set to 0 or 2!'
    inerr = 1
  end if

  if (iscale .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'iscale must == 0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support this NMR refinement option.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (noeskp .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'noeskp must == 1!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support this NMR refinement option.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (ipnlty .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'ipnlty must == 1!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support this NMR refinement option.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (mxsub .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'mxsub must == 1!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support this NMR refinement option.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (scalm .ne. 100.d0) then
    write(mdout, '(a,a)') error_hdr, 'scalm must == 100.d0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support this NMR refinement option.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (tausw .ne. 0.1d0) then
    write(mdout, '(a,a)') error_hdr, 'tausw must == 0.1d0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support this NMR refinement option.'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (icfe .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'icfe must == 0!'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' does not support free energies via thermo integration calcs!'
    write(mdout, '(a)') use_sander
    inerr = 1
  end if

  if (clambda .ne. 0.d0) then
    write(mdout, '(a,a)') error_hdr, 'clambda is only used if icfe != 0!'
    inerr = 1
  end if

  if (klambda .ne. 1) then
    write(mdout, '(a,a)') error_hdr, 'klambda is only used if icfe != 0!'
    inerr = 1
  end if

  if (ndfmin .ne. RETIRED_INPUT_OPTION) then
    write(mdout, '(a,a)') info_hdr, &
      'The ndfmin ctrl option is deprecated and ignored.'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' always sets this value based on other information.'
  end if

  if (dtemp .ne. RETIRED_INPUT_OPTION) then
    write(mdout, '(a,a)') info_hdr, &
      'The dtemp ctrl option is deprecated and ignored.'
  end if

  if (dxm .ne. RETIRED_INPUT_OPTION) then
    write(mdout, '(a,a)') info_hdr, &
      'The dxm ctrl option is deprecated and ignored.'
  end if

  if (heat .ne. RETIRED_INPUT_OPTION) then
    write(mdout, '(a,a)') info_hdr, &
      'The heat ctrl option is deprecated and ignored.'
  end if

  if (timlim .ne.  999999.d0) then
    write(mdout, '(a,a)') info_hdr, &
      'The timlim ctrl option is deprecated and ignored.'
  end if

  if (n3b .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'n3b may only be used if ipol != 0!'
    inerr = 1
  end if

  if (nion .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'nion may only be used if ipol != 0!'
    inerr = 1
  end if

  if (iesp .ne. 0) then
    write(mdout, '(a,a,a)') error_hdr, prog_name, &
      ' does not support iesp != 0!'
    inerr = 1
  end if

  if (ifcr .ne. 0) then
    write(mdout, '(a,a,a)') error_hdr, prog_name, &
      ' does not support ifcr != 0!'
    inerr = 1
  end if

  if (npscal .ne. 1) then
    write(mdout, '(a,a,a)') error_hdr, prog_name, &
      ' does not support npscal != 1!'
    write(mdout, '(a,a)') extra_line_hdr, &
      'Atom-based pressure scaling is no longer available.'
    inerr = 1
  end if

! Let the user know that lastist/lastrst serve no purpose in PMEMD:

  if (lastist .ne. 0) then
    write(mdout, '(a,a)') warn_hdr, &
    'The sander lastist option is not needed and is ignored.'
  end if

  if (lastrst .ne. 0) then
    write(mdout, '(a,a)') warn_hdr, &
    'The sander lastrst option is not needed and is ignored.'
  end if

  if (mdinfo_flush_interval .lt. 0) then
    write(mdout, '(a,a)') error_hdr, 'mdinfo_flush_interval must be >= 0!'
    inerr = 1
  else if (mdinfo_flush_interval .gt. 60 * 60) then
    write(mdout, '(a,a)') info_hdr, &
      'Excessive mdinfo_flush_interval reset to 1 hour!'
  end if

  if (mdout_flush_interval .lt. 0) then
    write(mdout, '(a,a)') error_hdr, 'mdout_flush_interval must be >= 0!'
    inerr = 1
  else if (mdout_flush_interval .gt. 60 * 60) then
    write(mdout, '(a,a)') info_hdr, &
      'Excessive mdout_flush_interval reset to 1 hour!'
  end if

  if (dbg_atom_redistribution .ne. 0 .and. dbg_atom_redistribution .ne. 1) then
    write(mdout, '(a,a)') error_hdr, &
      'dbg_atom_redistribution must be set to 0 or 1!'
    inerr = 1
  end if

! Report retired pmemd 3.1 options, but don't fail.

  if (amber7_compat .ne. RETIRED_INPUT_OPTION) then
    write(mdout, '(a,a)') info_hdr, &
      'The amber7_compat ctrl option is deprecated and ignored.'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' always runs in Amber 11 compatibility mode.'
  end if

  if (amber8_mdout_format .ne. RETIRED_INPUT_OPTION) then
    write(mdout, '(a,a)') info_hdr, &
      'The amber8_mdout_format ctrl option is deprecated and ignored.'
    write(mdout, '(a,a,a)') extra_line_hdr, prog_name, &
      ' always runs in Amber 11 mdout format mode.'
  end if

! Consistency checks:

  if ((nrespa .gt. 1 .or. nrespai .gt. 1) .and. imin .ne. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'For minimization, nrespa and nrespai must be 1 (default)!'
    inerr = 1
  end if

  if (nrespa .gt. 1 .and. ntp .ne. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'For constant pressure MD, nrespa must be 1!'
    inerr = 1
  end if

  if (tautp .lt. dt .and. ntt .eq. 1) then
    write(mdout, '(a,a)') error_hdr, 'tautp must be >= dt (step size)!'
    inerr = 1
  end if

  if (gamma_ln .gt. 0.d0 .and. ntt .ne. 3) then
    write(mdout, '(a,a)') error_hdr, 'ntt must be 3 if gamma_ln > 0!'
    inerr = 1
  end if

  if (ntp .eq. 0 .and. ntb .eq. 2) then
    write(mdout, '(a,a)') error_hdr, 'ntp must be 1 or 2 if ntb == 2!'
    inerr = 1
  end if

  if (ntp .ne. 0 .and. ntb .ne. 2) then
    write(mdout, '(a,a)') error_hdr, 'ntp must be 0 if ntb != 2!'
    inerr = 1
  end if 

  if (taup .lt. dt .and. ntp .ne. 0) then
    write(mdout, '(a,a)') error_hdr, 'taup must be >= dt (step size)!'
    inerr = 1
  end if

  if (ntb .ne. 0 .and. igb .gt. 0) then
    write(mdout, '(a,a)') error_hdr, &
                          ' igb > 0 is only compatible with ntb == 0!'
    inerr = 1
  end if

  if (igb .eq. 0) then

    ! User must either specify no cutoffs, or cut only, or both es_cutoff and
    ! vdw_cutoff.  If the user specified no cutoffs, then all values were set
    ! to defaults in the input subroutine.  If es_cutoff and vdw_cutoff are
    ! specified, es_cutoff must not be greater than vdw_cutoff.  We actually do
    ! allow specification of cut, es_cutoff and vdw_cutoff if they are all
    ! equal.

    if (cut .eq. 0.d0) then
      if (es_cutoff .eq. 0.d0 .or. vdw_cutoff .eq. 0.d0) then
        write(mdout, '(a,a)') error_hdr, &
          'Both es_cutoff and vdw_cutoff must be specified!'
        inerr = 1
      else if (es_cutoff .gt. vdw_cutoff) then
        write(mdout, '(a,a)') error_hdr, &
          'vdw_cutoff must be greater than es_cutoff!'
        inerr = 1
      end if
    else
      if (es_cutoff .ne. cut .or. vdw_cutoff .ne. cut) then
        write(mdout, '(a,a)') error_hdr, &
          'If cut is used then es_cutoff and vdw_cutoff must be consistent!'
        inerr = 1
      end if
    end if

    if (es_cutoff .le. 0.d0 .or. vdw_cutoff .le. 0.d0) then
      write(mdout, '(a,a)') error_hdr, &
            'Cutoffs (cut, es_cutoff, vdw_cutoff) must be positive!'
      inerr = 1
    end if

  else if (igb .eq. 1 .or. igb .eq. 2 .or. igb .eq. 5 .or. &
           igb .eq. 7 .or. igb .eq. 8) then

    if (gb_cutoff .lt. 8.05) then
      write(mdout, '(a,a)') error_hdr, &
        'Cut for Generalized Born simulation too small!'
      inerr = 1
    end if

  end if

!AMD initialization
!iamd=0 no boost is used, 1 boost on the total energy, 2 boost on the dohedrals, 3 boost on dihedrals and total energy
  if(iamd.gt.0)then
    if(iamd.eq.1)then !only total potential energy will be boosted
      EthreshD=0.0
      alphaD=0.0
    endif
    if(iamd.eq.2)then !only dihedral energy will be boosted
      EthreshP=0.0
      alphaP=0.0
    endif
    if (master) then
      write(mdout,'(a,i3)')'| Using Accelerated MD (AMD) to enhance sampling iamd =',iamd
      write(mdout,'(a,2f22.12)')'| AMD boost to total energy: EthreshP,alphaP',EthreshP,alphaP
      write(mdout,'(a,2f22.12)')'| AMD boost to dihedrals: EthreshD,alphaD',EthreshD,alphaD
    endif
  endif

!IPS validation
!the only valid ips option in pmemd is ips=1
   if (ips < 0 .or. ips > 1) then
      write(mdout,'(/2x,a,i3,a)') 'IPS (',ips,') must be between 0 and 1'
      inerr = 1
   end if
   if (ips /= 0 .and. ipol /= 0 ) then
      write(6,'(/2x,a)') 'IPS and IPOL are inconsistent options'
      inerr = 1
   endif
   if ( igb > 0 .and. ips > 0 ) then
      write(mdout,'(/2x,a,i3,a,i3,a)') 'IGB (',igb,') and ips (',ips, &
          ') cannot both be turned on'
      inerr = 1
   end if
   ! For IPS es and vdw cutoffs must be the same and are read through ips_cut 
   if (ips /= 0) then 
    if (cut .eq. 0.d0) then !if cut=0 set it
      if (es_cutoff .ne. vdw_cutoff) then
        write(mdout, '(a,a)') error_hdr, &
          'IPS only supports es_cutoff and vdw_cutoff being equal at the time!'
        inerr = 1
      else
        ips_cut = es_cutoff
      endif
    else
      ips_cut = cut
    endif
  endif

! Constant surface tension valid options
  if (csurften > 0) then
    if (csurften < 0 .or. csurften > 3) then
      write(mdout,'(/2x,a)') &
      'Invalid csurften value. csurften must be between 0 and 3'
      inerr = 1
    end if
    if (ntb /= 2) then
      write(mdout,'(/2x,a)') &
      'ntb (periodic boundary) invalid. ntb must be 2 for constant surface tension.'
      inerr = 1
    end if
    if (ntp < 2 .or. ntp > 3) then
      write(mdout,'(/2x,a)') &
      'ntp (constant pressure) invalid. ntp must be 2 or 3 for constant surface tension.'
      inerr = 1
    end if
    if (ninterface < 2) then
      write(mdout,'(/2x,a)') &
      'ninterface (number of interfaces) must be greater than 2 for constant surface tension.'
      inerr = 1
    end if
  end if

#ifdef MPI
  if (numexchg .gt. 0 .and. remd_method .eq. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'NUMEXCHG is a REMD variable and requires a non-zero value for -rem!'
    inerr = 1
  else if (numexchg .lt. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'NUMEXCHG must be non-negative!'
    inerr = 1
  else if (numexchg .eq. 0 .and. remd_method .ne. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'You are running a REMD simulation with no exchanges!'
    inerr = 1
  end if

  if (remd_method .ne. 0 .and. remd_method .ne. 1 .and. remd_method .ne. 3) then
    write(mdout, '(a,a)') error_hdr, &
      'The only -rem values currently supported are 0, 1, or 3!'
    inerr = 1
  end if

  if (remd_method .gt. 0) then

    if (exch_ntpr .lt. 0) then

      write(mdout, '(a,a)') error_hdr, 'exch_ntpr must be non-negative!'
      inerr = 1

    else if (exch_ntpr .eq. 0) then

      exch_ntpr = numexchg ! turn off printing except on last step

    end if

    if (ntp .gt. 0) then

      write(mdout, '(a,a)') error_hdr, 'REMD cannot be run with ntp > 0!'
      inerr = 1
    
    end if

    if (mod(numgroups, 2) .ne. 0) then

      write(mdout, '(a,a)') error_hdr, 'REMD requires an even number of &
         &replicas!'
      inerr = 1

    end if

  end if
#else
  if (numexchg .ne. 0) then
    write(mdout, '(a,a)') error_hdr, &
      'NUMEXCHG is a REMD variable that requires a parallel build of pmemd!'
    inerr = 1
  end if
#endif

#ifdef CUDA
!Trap specific limitations for NVIDIA GPU CUDA implementations.
!Current limitations are:
!
! ibelly /= 0                            !No support for belly.
! if (igb/=0) and cut < systemsize       !GPU sims do not use a cutoff in GB.
! gbsa /= 0                              !The GBSA method is not currently supported.
! nrespa /= 1                            !No multiple time stepping supported. 
! vlimit /= -1                           !vlimit is not supported on GPUs.
! es_cutoff /= vdw_cutoff                !Independent cutoffs for electrostatics
!                                        !and van der Waals are not supported on GPUs.
! order <= 4                             !Ewald interpolation order must be <= 4.
!
! Ewald Error Estimate is not printed    !This information is NOT calculated when running on GPU.
!
! Parallel CUDA limitations
!
! imin /= 0                              !Minimization is not supported in parallel on GPUs.
!

  if (ibelly /= 0) then
    write(mdout, '(a)') 'CUDA (GPU): Implementation does not support Belly options.'
    write(mdout, '(a)') '            Require ibelly == 0.'
    inerr = 1
  end if

! nmropt is currently calculated on the CPU
! if (nmropt /= 0) then
!   write(mdout, '(a)') 'CUDA (GPU): Implementation does not support nmropt options.'
!   write(mdout, '(a)') '            Require nmropt == 0.'
!   inerr = 1
! end if

  if (nrespa /= 1) then
    write(mdout, '(a)') 'CUDA (GPU): Implementation does not support nrespa.'
    write(mdout, '(a)') '            Require nrespa == 1.'
    inerr = 1
  end if
 
  if ( gbsa /= 0 ) then
    write(mdout, '(a)') 'CUDA (GPU): Implementation does not support GBSA method.'
    write(mdout, '(a)') '            Require gbsa == 0'
    inerr = 1
  end if
 
  !Check the size of the cutoff with respect to the system size.
  if (ntb == 0) then
    !Check not fully implemented. Ideally this check should loop over the current
    !coordinates and find the largest distance between any two atoms and make sure
    !that cut if larger than that value. For the moment though we shall just post
    !an error if cut looks to be too small in GB sims.
    if ( cut < 999.0d0 ) then
      write(mdout, '(a)') 'CUDA (GPU): Implementation does not support the use of a cutoff in GB simulations.'
      write(mdout, '(a)') '            Require cut > 999.0d0.'
      inerr = 1
    end if
  end if

  if (vlimit /= -1) then
      write(mdout, '(a)') 'CUDA (GPU): Implementation does not support the use of vlimit.'
      write(mdout, '(a)') '            Require vlimit == -1.0d0 (default).'
      inerr = 1
  end if

  if ( es_cutoff /= vdw_cutoff ) then
      write(mdout, '(a)') 'CUDA (GPU): Implementation does not support the use of multiple cutoffs.'
      write(mdout, '(a)') '            Require es_cutoff /= vdw_cutoff.'
      inerr = 1
  end if

#ifdef MPI
  if ( imin /= 0 ) then
      write(mdout, '(a)') 'CUDA (GPU): Minimization is NOT supported in parallel on GPUs. Please use'
      write(mdout, '(a)') '            the single GPU code for minimizations.'
      inerr = 1
  end if
#endif
#endif

! Field any errors and bag out.

  if (inerr .eq. 1) then
    write(mdout, '(/,a)') ' Input errors occurred. Terminating execution.'
    call mexit(6, 1)
  else
    write(mdout, '(a)') ' '
  end if

  return

end subroutine validate_mdin_ctrl_dat

!*******************************************************************************
!
! Subroutine:  print_mdin_ctrl_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine print_mdin_ctrl_dat(remd_method)

  implicit none

  integer, intent(in) :: remd_method

  write(mdout,'(/a)') 'General flags:'
  write(mdout,'(5x,2(a,i8))') 'imin    =', imin,', nmropt  =', nmropt

#ifdef MPI
  if (remd_method .ne. 0) then
    write(6,'(/a)') 'Replica exchange'
    write(6,'(5x,4(a,i8))') 'numexchg=',numexchg,', rem=',remd_method
  end if
#endif

  write(mdout,'(/a)') 'Nature and format of input:'
  write(mdout,'(5x,4(a,i8))') 'ntx     =', ntx, ', irest   =', irest, &
          ', ntrx    =', ntrx

  write(mdout,'(/a)') 'Nature and format of output:'
  write(mdout,'(5x,4(a,i8))') 'ntxo    =', ntxo,', ntpr    =', ntpr, &
          ', ntrx    =', ntrx, ', ntwr    =', ntwr
  write(mdout,'(5x,4(a,i8))') 'iwrap   =', iwrap,', ntwx    =', ntwx, &
          ', ntwv    =', ntwv,', ntwe    =', ntwe
  write(mdout,'(5x,3(a,i8),a,i7)') 'ioutfm  =', ioutfm, ', ntwprt  =', ntwprt, &
          ', idecomp =', idecomp, ', rbornstat=', rbornstat

  write(mdout,'(/a)') 'Potential function:'

  write(mdout,'(5x,5(a,i8))') 'ntf     =', ntf, ', ntb     =', ntb, &
          ', igb     =', igb, ', nsnb    =', nsnb
  write(mdout,'(5x,3(a,i8))') 'ipol    =', ipol, ', gbsa    =', gbsa, &
          ', iesp    =', iesp
  if (igb .eq. 0) then

    ! If there is only one cutoff, we use the old format for printout, just
    ! to create fewer output deltas; otherwise a new format shows the different
    ! values for es_cutoff and vdw_cutoff:

    if (es_cutoff .eq. vdw_cutoff) then
      write(mdout,'(5x,3(a,f10.5))') 'dielc   =', dielc, &
            ', cut     =', vdw_cutoff, ', intdiel =', intdiel
    else
      write(mdout,'(5x,2(a,f10.5))') 'es_cutoff   =', es_cutoff, &
            ', vdw_cutoff     =', vdw_cutoff
      write(mdout,'(5x,2(a,f10.5))') 'dielc   =', dielc, ', intdiel =', intdiel
    end if

  else if (igb .eq. 1 .or. igb .eq. 2 .or. igb .eq. 5 .or. igb .eq. 7) then

      write(mdout,'(5x,3(a,f10.5))') 'dielc   =', dielc, &
            ', cut     =', gb_cutoff, ', intdiel =', intdiel
      write(mdout,'(5x,3(a,f10.5))') 'saltcon =', saltcon, &
            ', offset  =', offset, ', gbalpha= ', gb_alpha
      write(mdout,'(5x,3(a,f10.5))') 'gbbeta  =', gb_beta, &
            ', gbgamma =', gb_gamma, ', surften =', surften
      write(mdout,'(5x,3(a,f10.5))') 'rdt     =', rdt, &
            ', rgbmax  =', rgbmax,'  extdiel =',extdiel
      write(mdout,'(5x,3(a,i8))') 'alpb  = ', alpb

      if (alpb .ne. 0) then
        write(mdout,'(5x,3(a,f10.5))') 'Arad =', Arad
      end if

  else if (igb .eq. 8) then

      write(mdout, '(5x,3(a,f10.5))') 'dielc   =', dielc, &
        ', cut     =', gb_cutoff, ', intdiel =', intdiel
      write(mdout, '(5x,3(a,f10.5))') 'saltcon =', saltcon, &
        ', offset  =', offset, ', surften =', surften
      write(mdout,'(5x,3(a,f10.5))') 'rdt     =', rdt, ', rgbmax  =', rgbmax, &
        '  extdiel =', extdiel
      write(mdout,'(5x,3(a,i8))') 'alpb  = ', alpb
      write(mdout,'(5x,3(a,f10.5))') 'gbalphaH  =', gb_alpha_h, &
        ', gbbetaH   =', gb_beta_h, ',  gbgammaH  = ', gb_gamma_h
      write(mdout,'(5x,3(a,f10.5))') 'gbalphaC  =', gb_alpha_c, &
        ', gbbetaC   =', gb_beta_c, ',  gbgammaC  = ', gb_gamma_c
      write(mdout,'(5x,3(a,f10.5))') 'gbalphaN  =', gb_alpha_n, &
        ', gbbetaN   =', gb_beta_n, ',  gbgammaN  = ', gb_gamma_n
      write(mdout,'(5x,3(a,f10.5))') 'gbalphaOS =', gb_alpha_os, &
        ', gbbetaOS  =', gb_beta_os, ',  gbgammaOS = ', gb_gamma_os
      write(mdout,'(5x,3(a,f10.5))') 'gbalphaP  =', gb_alpha_p, &
        ', gbbetaP   =', gb_beta_p, ',  gbgammaP  = ', gb_gamma_p
      
  end if
  
  write(mdout,'(/a)') 'Frozen or restrained atoms:'

  write(mdout,'(5x,4(a,i8))') 'ibelly  =', ibelly,', ntr     =', ntr
  if( ntr == 1 ) write(6,'(5x,a,f10.5)') 'restraint_wt =', restraint_wt

  if (imin .ne. 0) then
    write(mdout,'(/a)') 'Energy minimization:'
    ! print inputable variables applicable to all minimization methods.
    write(mdout,'(5x,4(a,i8))') 'maxcyc  =', maxcyc, ', ncyc    =', ncyc, &
            ', ntmin   =', ntmin
    write(mdout,'(5x,2(a,f10.5))') 'dx0     =', dx0, ', drms    =', drms
  else
    write(mdout,'(/a)') 'Molecular dynamics:'
    write(mdout,'(5x,3(a,i10))') 'nstlim  =', nstlim, ', nscm    =', nscm, &
            ', nrespa  =', nrespa
    write(mdout,'(5x,3(a,f10.5))') 't       =', original_t, &
            ', dt      =', dt,', vlimit  =', vlimit
      
    if (ntt .eq. 1) then
      write(mdout,'(/a)') 'Berendsen (weak-coupling) temperature regulation:'
      write(mdout,'(5x,3(a,f10.5))') 'temp0   =', temp0, &
              ', tempi   =', tempi, ', tautp   =', tautp
    else if (ntt .eq. 2) then
      write(mdout,'(/a)') 'Anderson (strong collision) temperature regulation:'
      write(mdout,'(5x,4(a,i8))') 'ig      =', ig, ', vrand   =', vrand
      write(mdout,'(5x,3(a,f10.5))') 'temp0   =', temp0, ', tempi   =', tempi
    else if (ntt .eq. 3) then
      write(mdout,'(/a)') 'Langevin dynamics temperature regulation:'
      write(mdout,'(5x,4(a,i8))') 'ig      =', ig
      write(mdout,'(5x,3(a,f10.5))') 'temp0   =', temp0, &
               ', tempi   =', tempi,', gamma_ln=',  gamma_ln
    end if

    if (ntp .ne. 0) then
      write(mdout,'(/a)') 'Pressure regulation:'
      write(mdout,'(5x,a,i8)') 'ntp     =', ntp
      write(mdout,'(5x,3(a,f10.5))') 'pres0   =', pres0, &
               ', comp    =', comp, ', taup    =', taup
    end if

    if (csurften /= 0) then
      write(mdout,'(/a)') 'Constant surface tension:'
      write(mdout,'(5x,a,i8)') 'csurften  =', csurften
      write(mdout,'(5x,a,f10.5,a,i8)') 'gamma_ten =', gamma_ten, ' ninterface =', ninterface
    end if

  end if

  if (ntc .ne. 1) then
    write(mdout,'(/a)') 'SHAKE:'
    write(mdout,'(5x,4(a,i8))') 'ntc     =', ntc, ', jfastw  =', jfastw
    write(mdout,'(5x,3(a,f10.5))') 'tol     =', tol
  end if

  if (nmropt .gt. 0) then
    write(mdout,'(/a)') 'NMR refinement options:'
    write(mdout,'(5x,4(a,i8))')'iscale  =', iscale, ', noeskp  =', noeskp, &
            ', ipnlty  =', ipnlty, ', mxsub   =', mxsub
    write(mdout,'(5x,3(a,f10.5))') 'scalm   =', scalm, &
            ', pencut  =', pencut, ', tausw   =', tausw
  end if

! Polarization not currently supported...
!
!if (ipol .ne. 0) then
!  write(mdout,'(/a)') 'Polarizable options:'
!  write(mdout,'(5x,4(a,i8))') 'indmeth =', indmeth, &
!          ', maxiter =', maxiter,', irstdip =', irstdip, &
!          ', scaldip =', scaldip
!  write(mdout,'(5x,3(a,f10.5))') &
!          'diptau  =', diptau,', dipmass =', dipmass
!end if

! Free energies via thermodynamic integration not currently supported...
!
! if (icfe .ne. 0) then
!   write(mdout,'(/a)') 'Free energy options:'
!   write(mdout,'(5x,4(a,i8))') 'klambda =', klambda
!   write(mdout,'(5x,3(a,f10.5))') 'clambda =', clambda
! end if

! Targetted MD not currently supported...
!
! if (itgtmd .ne. 0) then
!   write(mdout,'(/a)') 'Targeted molecular dynamics:'
!   write(mdout,'(5x,3(a,f10.5))') 'tgtrmsd =', tgtrmsd, &
!           ', tgtmdfrc=', tgtmdfrc
! end if

! Intermolecular bond treatment:

  write(mdout,'(/a)') '| Intermolecular bonds treatment:'
  write(mdout,'(a,i8)') '|     no_intermolecular_bonds =', &
    no_intermolecular_bonds

! Frequency of averaging energies:

  write(mdout,'(/a)') '| Energy averages sample interval:'
  write(mdout,'(a,i8)') '|     ene_avg_sampling =', &
    ene_avg_sampling

  return

end subroutine print_mdin_ctrl_dat

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_mdin_ctrl_dat
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_mdin_ctrl_dat

  use parallel_dat_mod

  implicit none

  call mpi_bcast(imin, mdin_ctrl_int_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)
  call mpi_bcast(dielc, mdin_ctrl_dbl_cnt, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

! Broadcast the NO_NTT3_SYNC Flag
  call mpi_bcast(no_ntt3_sync, 1, mpi_integer, 0, pmemd_comm, err_code_mpi)

! Set the potential flags:

  if (igb .eq. 0) then
    using_gb_potential = .false.
    using_pme_potential = .true.
  else if (igb .eq. 1 .or. igb .eq. 2 .or. igb .eq. 5 .or. &
           igb .eq. 7 .or. igb .eq. 8) then
    using_gb_potential = .true.
    using_pme_potential = .false.
  else                                  ! in a world of trouble...
    using_pme_potential = .false.
    using_gb_potential = .false.
  end if

  if ( no_ntt3_sync ) then
    !Here we are not synching the random number generator across threads for
    !NTT=3 but we don't want every thread to have the same random seed. So for
    !the moment just add mytask id to ig.
#ifndef CUDA
    !In the CUDA case this seems to cause massive problems with NTT3 in parallel so
    !only do it for regular CPU runs. NO_NTT3_SYNC shouldn't help much with GPUs
    !anyway since they use their own random number system.
    ig = ig+mytaskid
#endif
  end if

  return

end subroutine bcast_mdin_ctrl_dat
#endif

end module mdin_ctrl_dat_mod
