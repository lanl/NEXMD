#include "dprec.fh"
module qm2_extern_tc_module
! ----------------------------------------------------------------
! Interface for TeraChem based QM and QM/MM MD 
!
! Currently supports:
! pure QM
! QM/MM with cutoff for QM-MM electrostatics under periodic
! boundary conditions
!
! Author: Andreas Goetz (agoetz@sdsc.edu)
!
! Date: November 2010
!
! ----------------------------------------------------------------

  implicit none

  private
  public :: get_tc_forces, tc_finalize
  logical, save :: do_mpi = .false.  ! Used in finalize subroutine

  type tc_nml_type
     character(len=20) :: basis
     character(len=20) :: method
     character(len=20) :: precision
     character(len=20) :: executable
     character(len=20) :: dftd
     character(len=20) :: guess
     character(len=20) :: cis
     _REAL_ :: threall
     _REAL_ :: convthre
     integer :: charge
     integer :: spinmult
     integer :: maxit
     integer :: dftgrid
     integer :: ngpus
     integer, dimension(:), pointer :: gpuids => null()
     integer :: cisnumstates
     integer :: cistarget
     integer :: mpi
     integer :: ntpr
     integer :: verbosity
     logical :: dipole
     logical :: use_template
  end type tc_nml_type

  integer, save         :: newcomm ! Initialized in mpi_init subroutine

contains

  ! --------------------------------------
  ! Get QM energy and forces from TeraChem
  ! --------------------------------------
  subroutine get_tc_forces( do_grad, nstep, ntpr_default, id, nqmatoms, qmcoords,&
    nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge )

    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO

    logical, intent(in) :: do_grad              ! Return gradient/not
    integer, intent(in) :: nstep                ! MD step number
    integer, intent(in) :: ntpr_default         ! frequency of printing
    character(len=3), intent(in) :: id          ! ID number for PIMD or REMD
    integer, intent(in) :: nqmatoms             ! Number of QM atoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM coordinates
    integer, intent(in) :: nclatoms             ! Number of MM atoms
    _REAL_,  intent(in) :: clcoords(4,nclatoms) ! MM coordinates and charges in au
    _REAL_, intent(out) :: escf                 ! SCF energy
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM force
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM force
    integer, intent(out) :: charge
    _REAL_              :: dipmom(4,3)          ! Dipole moment {x, y, z, |D|}, {QM, MM, TOT}

    type(tc_nml_type), save :: tc_nml
    logical, save :: first_call = .true.
    integer :: i
    integer :: printed =-1 ! Used to tell if we have printed this step yet 
                           ! since the same step may be called multiple times
    character(len=150) :: call_buffer
    character(len=*),  parameter :: basename = 'tc_job'
    character(len=*),  parameter :: runext = '.inp'
    character(len=*),  parameter :: datext = '.dat'
    character(len=*),  parameter :: dipext = '.dip'
    character(len=*),  parameter :: tplext = '.tpl'
    character(len=14) :: inpfile, datfile, dipfile, crdfile, chrgfile, tplfile
    character(len=25)            :: subdir ! May need to prepend subdirectory
    
    ! assemble input - / output data filenames
    inpfile = basename//trim(id)//runext 
    datfile = basename//trim(id)//datext
    dipfile = basename//trim(id)//dipext 
    crdfile = 'inpfile'//trim(id)//'.xyz' 
    chrgfile  = 'ptchrg'//trim(id)//'.xyz'
    tplfile = basename//tplext 
    
    ! Setup on first program call
    if ( first_call ) then
      first_call = .false.
      write (6,'(/,a,/)') '   >>> Running QM calculation with TeraChem <<<'
      call get_namelist( ntpr_default, tc_nml )
      call check_installation(trim(tc_nml%executable))
      call print_namelist(tc_nml)
      charge = tc_nml%charge
      write (6,'(80a)') ('-', i=1,80)
      write (6,'(a)') '   4.  RESULTS'
      write (6,'(80a)') ('-', i=1,80)
      call system('rm -f '//dipfile)
    end if

#ifdef MPI
# ifndef MPI_1
    if(tc_nml%mpi==1) then ! Do mpi (forced to 0 ifndef MPI)
       call mpi_hook( trim(tplfile), nqmatoms, qmcoords, nclatoms, clcoords,&
            tc_nml, escf, dxyzqm, dxyzcl, dipmom, do_grad, id )
    else
# else
    ! If we are using MPI 1.x the code will not compile since
    ! MPI_LOOKUP_NAME is part of the MPI 2 standard, so  just quit
    if(tc_nml%mpi==1) then 
    call sander_bomb('(qm2_extern_tc_module)', &
        '&unsupported MPI version', &
        'Will quit now.')
    else
# endif
#endif
      
       call system('rm -f '//inpfile)
       call write_inpfile( trim(inpfile), trim(crdfile), trim(chrgfile), trim(tplfile), &
            nqmatoms, qmcoords, nclatoms, clcoords, tc_nml, do_grad )

       call_buffer=''
       subdir=''
       if(trim(id)/='') then
          subdir='./'//trim(id)//'/'
          call_buffer=' mkdir -p '//trim(subdir)//&
                    '; cd '//trim(subdir)//&
                    '; mv ../'//inpfile//' .;'//&
                     ' mv ../'//crdfile//' .;'
          if(nclatoms > 0) then
             call_buffer=trim(call_buffer)//' mv ../'//chrgfile//' .;'
          end if
       end if

      ! Run TeraChem with file inpfile
       call system('rm -f '//datfile)
       write(call_buffer,'(2a)') &
            trim(call_buffer),'$TeraChem/'//trim(tc_nml%executable)//' '//trim(inpfile)//' > '//datfile
      
       call system(trim(call_buffer))

       ! If working in a subdirectory, move datfile back for reading
       if(trim(id)/='') then
          call system('mv '//trim(subdir)//datfile//' .;')
       end if

       ! Read TeraChem results
       call read_results(trim(datfile), nqmatoms, nclatoms, escf, dxyzqm, dxyzcl, dipmom, do_grad)
 
       ! Write dipole moment to file
       if ( tc_nml%ntpr > 0 .and. mod(nstep, tc_nml%ntpr) == 0 ) then
          if( printed /= nstep .and. tc_nml%dipole ) then
             call write_results(trim(dipfile), dipmom)
             printed = nstep
          end if
       end if

       call system('mv '//trim(subdir)//trim(inpfile)//' '//trim(subdir)//'old.'//inpfile)
       call system('mv '//trim(datfile)//' '//trim(subdir)//'old.'//datfile)
#ifdef MPI
    end if
#endif

    if(tc_nml%verbosity > 0)then
       write(6,'(/,a,/)') 'get_tc_forces - gradient from terachem in au:'
       write(6,'(a)') 'Quantum region:'
       do i = 1, nqmatoms
          write(6,'(3(x,f16.10))') dxyzqm(:, i)
       end do
       if (nclatoms > 0 ) then
          write(6,'(a)') 'MM region:'
          do i = 1, nclatoms
             write(6,'(3(x,f16.10))') dxyzcl(:, i)
          end do
       end if
    end if

    if ( do_grad ) then
       ! Convert Hartree/Bohr -> kcal/(mol*A)
       dxyzqm(:,:) = dxyzqm(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
       if ( nclatoms > 0 ) then
          dxyzcl(:,:) = dxyzcl(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
       end if
    else
       dxyzqm = ZERO
       if ( nclatoms > 0 ) dxyzcl = ZERO
    end if

    if(tc_nml%verbosity > 0)then
       write(6,'(/,a,/)') 'get_tc_forces - gradient in kcal/mol/bohr:'
       write(6,'(a)') 'Quantum region:'
       do i = 1, nqmatoms
          write(6,'(3(x,f16.10))') dxyzqm(:, i)
       end do
       if (nclatoms > 0 ) then
          write(6,'(a)') 'MM region:'
          do i = 1, nclatoms
             write(6,'(3(x,f16.10))') dxyzcl(:, i)
          end do
       end if
    end if
    
    escf = escf * CODATA08_AU_TO_KCAL

  end subroutine get_tc_forces

  ! -----------------------------------------------
  ! Read TeraChem tc namelist values from file mdin,
  ! use default values if none are present.
  ! -----------------------------------------------
  subroutine get_namelist( ntpr_default, tc_nml)

    use UtilitiesModule, only: Upcase
    implicit none

    integer, intent(in) :: ntpr_default
    type(tc_nml_type), intent(out) :: tc_nml
    character(len=20):: basis, method, dftd, precision, executable, guess, cis
    _REAL_ :: threall, convthre
    integer, parameter :: maxgpus = 64
    integer :: charge, spinmult, maxit, dftgrid, ngpus, gpuids(maxgpus),  &
         cisnumstates, cistarget, mpi, ntpr, verbosity, dipole, use_template
    namelist /tc/ basis, method, dftd, precision, executable, guess, cis, &
         threall, convthre, &
         charge, spinmult, maxit, dftgrid, ngpus, gpuids, cisnumstates, cistarget, &
           mpi, ntpr, verbosity, dipole, use_template
    integer :: i, ierr

    ! Default values
    basis        = '6-31g'
    method       = 'blyp'
    dftd         = 'no'
    precision    = 'mixed'
    executable   = 'terachem'
    guess        = 'scr/c0'
    cis          = 'no'
    threall      = 1.0d-11
    convthre     = 3.0d-05
    charge       = 0
    spinmult     = 1
    maxit        = 100
    dftgrid      = 1
    ngpus        = 0 ! Use all available GPUs
    do i = 1, maxgpus
       gpuids(i) = i-1
    end do
    cisnumstates = 1
    cistarget    = 1
    mpi          = 1 ! Default to using MPI if available
    ntpr         = ntpr_default
    verbosity    = 0
    dipole       = 0
    use_template = 0

    ! Read namelist
    rewind 5
    read(5,nml=tc,iostat=ierr)

    if ( ierr > 0 ) then
       call sander_bomb('get_namelist (qm2_extern_tc_module)', &
            '&tc namelist read error', &
            'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a/a)') '&tc namelist read encountered end of file', &
            'Please check your input if the calculation encounters a problem'
    end if

    ! Assign namelist values to tc_nml data type
    tc_nml%basis        = Upcase(basis)
    tc_nml%method       = method
    tc_nml%dftd         = dftd
    tc_nml%precision    = precision
    tc_nml%executable   = executable
    tc_nml%guess        = guess
    tc_nml%cis          = cis
    tc_nml%threall      = threall
    tc_nml%convthre     = convthre
    tc_nml%charge       = charge
    tc_nml%spinmult     = spinmult
    tc_nml%maxit        = maxit
    tc_nml%dftgrid      = dftgrid
    tc_nml%ngpus        = ngpus
    tc_nml%cisnumstates = cisnumstates
    tc_nml%cistarget    = cistarget
    if ( ngpus > 0 ) then
       allocate ( tc_nml%gpuids(ngpus), stat = ierr )
       if ( ierr /= 0 ) then
          call sander_bomb('get_namelist (qm2_extern_tc_module)', &
               'Allocation error for gpuids(:)', &
               'Will quit now')
       end if
       tc_nml%gpuids(:) = gpuids(:ngpus)
    end if
#ifndef MPI
        if ( tc_nml%mpi == 1 ) then
          write(6,'(a)') '| Warning: mpi=1 selected but sander was not compiled with MPI support.'
          write(6,'(a)') '| Continuing with mpi=0'
        end if
        tc_nml%mpi         = 0 ! Can't pick MPI if not available 
#else
        tc_nml%mpi         = mpi
#endif

    ! Need this variable so we don't call MPI_Send in the finalize subroutine
    if(mpi==1) then
      do_mpi=.true.
    end if

    tc_nml%ntpr      = ntpr
    tc_nml%verbosity = verbosity

    if ( dipole == 0 ) then
       tc_nml%dipole = .false.
    else if ( dipole == 1 ) then
       tc_nml%dipole = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_tc_module)', &
            '&tc dipole value not allowed', &
            'Please check your input. dipole can only be 0 or 1.')
    end if

    if ( use_template == 0 ) then
       tc_nml%use_template = .false.
    else if ( use_template == 1 ) then
       tc_nml%use_template = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_tc_module)', &
            '&tc use_template value not allowed', &
            'Please check your input. use_template can only be 0 or 1.')
    end if


  end subroutine get_namelist

  ! --------------------------------
  ! Print TeraChem namelist settings
  ! --------------------------------
  subroutine print_namelist( tc_nml )

    implicit none
    type(tc_nml_type), intent(in) :: tc_nml

    integer :: i, j, jstart, jend
    integer, parameter :: jstep = 10
    character(len=30) :: tmpstr

    write(6, '(/,a)')      '| &tc'
    write(6, '(2a)')       '|   basis        = ', tc_nml%basis
    write(6, '(2a)')       '|   method       = ', tc_nml%method
    write(6, '(2a)')       '|   dftd         = ', tc_nml%dftd
    write(6, '(2a)')       '|   precision    = ', tc_nml%precision
    write(6, '(2a)')       '|   executable   = ', tc_nml%executable
    write(6, '(2a)')       '|   guess        = ', tc_nml%guess
    write(6, '(2a)')       '|   cis          = ', tc_nml%cis
    write(6, '(a,es10.2)') '|   threall      = ', tc_nml%threall
    write(6, '(a,es10.2)') '|   convthre     = ', tc_nml%convthre
    write(6, '(a,i4)')     '|   charge       = ', tc_nml%charge
    write(6, '(a,i4)')     '|   spinmult     = ', tc_nml%spinmult
    write(6, '(a,i4)')     '|   maxit        = ', tc_nml%maxit
    write(6, '(a,i4)')     '|   dftgrid      = ', tc_nml%dftgrid
    write(6, '(a,i4)')     '|   ngpus        = ', tc_nml%ngpus
    write(6, '(a,i4)')     '|   cisnumstates = ', tc_nml%cisnumstates
    write(6, '(a,i4)')     '|   cistarget    = ', tc_nml%cistarget
    if ( tc_nml%ngpus > 0 ) then
       jstart = 1
       do i = 1, tc_nml%ngpus / jstep + 1
          if ( i == 1 ) then
             tmpstr =      '|   gpuids       = '
          else
             tmpstr =      '                  '
          end if
          jend = min ( (jstart + jstep - 1), tc_nml%ngpus )
          write(6,'(a,9999(i5))') tmpstr, (tc_nml%gpuids(j), j = jstart, jend)
          jstart = jstart + jstep
       end do
    end if
    write(6, '(a,i1)')     '|   mpi          = ', tc_nml%mpi
    write(6, '(a,i0)')     '|   ntpr         = ', tc_nml%ntpr
    write(6, '(a,i2)')     '|   verbosity    = ', tc_nml%verbosity
    write(6, '(a,l)')      '|   dipole       = ', tc_nml%dipole
    write(6, '(a,l)')      '|   use_template = ', tc_nml%use_template
    write(6,'(a)')         '| /'

  end subroutine print_namelist

  ! --------------------------------------------------------
  ! Check whether ADF/GAMESS/TeraChem are properly installed
  ! --------------------------------------------------------
  subroutine check_installation( program )

    implicit none

    character(len=*), intent(in) :: program

    logical :: exist
    character(len=80) :: env_path
    integer :: len, stat
    
    ! Search for $TeraChem
    call get_environment_variable('TeraChem', env_path, len, stat)
    if ( stat /= 0 ) then
       call sander_bomb('check_installation (qm2_extern_tc_module)', &
            'Environment variable "$TeraChem" not set', &
            'Please check you TeraChem installation')
    end if

    env_path = trim(env_path)//'/'//program
    write(6,'(2a)') '| Searching for ', env_path
    inquire (FILE=env_path, EXIST=exist)

    if(exist .eqv. .true.) then
       write(6,'(3a)') '| Program ', program, ' found!'
    else
       call sander_bomb('check_installation (qm2_extern_tc_module)', &
            'Program '//trim(program)//' not found', &
            'External program "'//trim(program)//'" required to use the extern_tc module.')
    end if

  end subroutine check_installation


#if defined(MPI) && !defined(MPI_1)
  ! Perform MPI communications with terachem. Requires MPI 2.0 or above to use
  subroutine mpi_hook( tplfile, nqmatoms, qmcoords, nclatoms, clcoords,&
       tc_nml, escf, dxyzqm, dxyzcl, dipmom, do_grad, id )
    
    use qmmm_module, only : qmmm_struct
    use ElementOrbitalIndex, only : elementSymbol
    
    implicit none
    include 'mpif.h'

    character(len=*), intent(in)  :: tplfile
    integer, intent(in) :: nqmatoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) 
    integer, intent(in) :: nclatoms
    _REAL_,  intent(in) :: clcoords(4,nqmatoms)
    type(tc_nml_type), intent(in) :: tc_nml
    _REAL_, intent(out) :: escf
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)
    _REAL_, intent(out) :: dipmom(4,3)
    logical, intent(in) :: do_grad
    character(len=3), intent(in) :: id

    character(len=2)    :: atom_types(nqmatoms)
    _REAL_              :: coords(3,nqmatoms+nclatoms)
    _REAL_              :: charges(nclatoms)
    _REAL_              :: dxyz_all(3,nclatoms+nqmatoms)

    _REAL_              :: buf(3)
    logical,save        :: first_call=.true.
    integer             :: i, status(MPI_STATUS_SIZE)
    integer             :: ierr

    character*(255)       :: port_name
    character(len=128) :: dbuffer(2,32)

    ! AWG FIXME: write charges if requested
    ! for now just receive them locally to make the interface work
    _REAL_ :: qmcharges(nqmatoms)

    if (tc_nml%verbosity > 0) then
       write (6, '(a)') '>>>>> Entered mpi_hook() (qm2_extern_tc_module)'
       call flush(6)
    end if

    ! Determine atom types
    ! TeraChem needs those both for initialization and later during the MD run
    do i = 1, nqmatoms
      atom_types(i)=elementSymbol(qmmm_struct%iqm_atomic_numbers(i))
   end do

    ! ---------------------------------------------------
    ! Initialization: Connect to "terachem_port", set    
    ! newcomm (global), send relevant namelist variables.
    ! ---------------------------------------------------
    if(first_call) then 
      first_call=.false.
      call connect_to_terachem( tplfile, tc_nml, nqmatoms, atom_types, do_grad, id )
    end if

    ! -----------------------------------------
    ! Begin sending data each step to terachem
    ! -----------------------------------------
    if (tc_nml%verbosity > 0) then
       write(6,'(a)') 'Sending data to TeraChem'
       call flush(6)
    end if

    ! Send nqmatoms and the type of each qmatom
    if (tc_nml%verbosity > 1) then
       write(6,'(/, a, i0)') 'Sending nqmatoms = ', nqmatoms
       call flush(6)
    end if
    call MPI_Send( nqmatoms, 1, MPI_INTEGER, 0, 2, newcomm, ierr )

    if (tc_nml%verbosity > 1) then
       write(6,'(/,a)') 'Sending QM atom types: '
       do i = 1, nqmatoms
          write(6,'(a)') atom_types(i)
          call flush(6)
       end do
    end if
    call MPI_Send( atom_types, 2*size(atom_types), MPI_CHARACTER, 0, 2, newcomm, ierr )

    ! Send QM coordinate array
    if (tc_nml%verbosity > 1) then
       write(6,'(a)') 'Sending QM coords: '
    end if
    do i=1, nqmatoms
       if (tc_nml%verbosity > 1) then
          write(6,*) 'Atom ',i,': ',qmcoords(:,i)
          call flush(6)
       end if
    end do 
    call MPI_Send( qmcoords, 3*nqmatoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

    ! Send nclatoms and the charge of each atom
    if (tc_nml%verbosity > 1) then
       write(6,'(a, i0)') 'Sending nclatoms = ', nclatoms
       call flush(6)
    end if
    call MPI_Send( nclatoms, 1, MPI_INTEGER, 0, 2, newcomm, ierr ) 

    if (tc_nml%verbosity > 1) then
       write(6,'(a)') 'Sending charges: '
    end if
    do i=1, nclatoms
      charges(i)=clcoords(4,i)
      if (tc_nml%verbosity > 1) then
         write(6,*) 'Charge ',i,':',charges(i)
         call flush(6)
      end if
    end do
    call MPI_Send( charges, nclatoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

    ! Send MM point charge coordinate array
    if (tc_nml%verbosity > 1) then
       write(6,'(a)') 'Sending CL coords: '
    end if
    do i=1, nclatoms
      coords(:,i)=clcoords(:3,i)
      if (tc_nml%verbosity > 1) then
         write(6,*) 'Atom ',i,': ',coords(:,i)
         call flush(6)
      end if
    end do 
    call MPI_Send( coords, 3*nclatoms, MPI_DOUBLE_PRECISION, 0, 2, newcomm, ierr ) 

    ! -----------------------------------
    ! Begin receiving data from terachem
    ! -----------------------------------

    ! Energy
    if (tc_nml%verbosity > 1) then
       write(6,'(a)') 'Waiting to receive scf energy from TeraChem...'
       call flush(6)
    end if
    call MPI_Recv( escf, 1, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if (tc_nml%verbosity > 0) then
       write(6,'(a,es15.6)') 'Received scf energy from server:', escf
       call flush(6)
    end if

    ! Charges (Mulliken or other)
    if (tc_nml%verbosity > 1) then
       write(6,'(a)') 'Waiting to receive charges...'
    end if
    call MPI_Recv( qmcharges(:), nqmatoms, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if (tc_nml%verbosity > 1) then
       write(6,'(a)') 'Received the following charges from server:'
       do i=1, nqmatoms
          write(6,*) 'Atom ',i, ': ', qmcharges(i)
       end do
       call flush(6)
    end if

    ! Dipole moment
    if (tc_nml%verbosity > 1) then
       write(6,'(a)') 'Waiting to receive dipole moment...'
    end if
    ! QM dipole moment
    call MPI_Recv( dipmom(:,1), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if (tc_nml%verbosity > 0) then
       write(6,'(a,4es15.6)') 'Received QM  dipole moment from server:', dipmom(:,1)
       call flush(6)
    end if
    ! MM dipole moment
    call MPI_Recv( dipmom(:,2), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if (tc_nml%verbosity > 0) then
       write(6,'(a,4es15.6)') 'Received MM  dipole moment from server:', dipmom(:,2)
       call flush(6)
    end if
    ! TOT dipole moment
    call MPI_Recv( dipmom(:,3), 4, MPI_DOUBLE_PRECISION, MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
    if (tc_nml%verbosity > 0) then
       write(6,'(a,4es15.6)') 'Received TOT dipole moment from server:', dipmom(:,3)
       call flush(6)
    end if
    
    ! QM gradients
    if ( do_grad ) then
       if (tc_nml%verbosity > 1) then
          write(6,'(a)') 'Waiting to receive gradients...'
       end if
       call MPI_Recv( dxyz_all, 3*(nqmatoms+nclatoms), MPI_DOUBLE_PRECISION, &
            MPI_ANY_SOURCE, MPI_ANY_TAG, newcomm, status, ierr )
       if (tc_nml%verbosity > 1) then
          write(6,'(a)') 'Received the following gradients from server:'
          do i=1, nqmatoms+nclatoms
             write(6,*) 'Atom ',i, ': ',dxyz_all(:,i)
          end do
          call flush(6)
       end if
       
       ! Poplulate our output arrays with gradients from terachem
       do i=1, nqmatoms
          dxyzqm(:,i)=dxyz_all(:,i)
       end do
       do i=1, nclatoms
          dxyzcl(:,i)=dxyz_all(:,i+nqmatoms)
       end do

    end if

    if (tc_nml%verbosity > 0) then
       write (6, '(a)') '<<<<< Leaving mpi_hook()'
       call flush(6)
    end if

  end subroutine mpi_hook

  ! -------------------------------------------------
  ! Search for name published by TeraChem and connect
  ! (this step initializes newcomm)
  ! Send relevant namelist variables to terachem
  ! -------------------------------------------------
  subroutine connect_to_terachem( tplfile, tc_nml, nqmatoms, atom_types, do_grad, id )

    implicit none
    include 'mpif.h'

    character(len=*), intent(in)  :: tplfile
    type(tc_nml_type), intent(in) :: tc_nml
    integer, intent(in) :: nqmatoms
    character(len=2), intent(in) :: atom_types(nqmatoms)
    logical, intent(in) :: do_grad
    character(len=3), intent(in) :: id

    character(len=17) :: server_name="terachem_port"
    integer, parameter  :: clen=128 ! Length of character strings we are using
    character(255) :: port_name
    character(len=clen) :: dbuffer(2,32)
    _REAL_          :: timer
    integer         :: ierr, i, j, irow
    logical         :: done=.false.

    _REAL_ :: escf

    if (tc_nml%verbosity > 0) then
       write(6,'(a)') '>>>>> Entered connect_to_terachem() (qm2_extern_tc_module)'
       call flush(6)
    end if

    ! -----------------------------------
    ! Look for server_name, get port name
    ! After 60 seconds, exit if not found
    ! -----------------------------------
    if ( trim(id) /= '' ) then
       server_name = trim(server_name)//'.'//trim(id)
    end if
    if (tc_nml%verbosity > 0) then
       write(6,'(2a)') 'Looking up server under name:', trim(server_name)
       call flush(6)
    end if
    timer = MPI_WTIME(ierr)
    do while (done .eqv. .false.)

       call MPI_LOOKUP_NAME(trim(server_name), MPI_INFO_NULL, port_name, ierr)
       if (ierr == MPI_SUCCESS) then
          if ( tc_nml%verbosity > 0 ) then
             write(6,'(2a)') 'Found port: ', trim(port_name)
             call flush(6)
          end if
          done=.true.

       end if

       if ( (MPI_WTIME(ierr)-timer) > 60 ) then ! Time out after 60 seconds
          call sander_bomb('connect_to_terachem() (qm2_extern_tc_module)', &
               '"'//trim(server_name)//'" not found. Timed out after 60 seconds.', &
               'Will quit now')
       end if

    end do

    ! ----------------------------------------
    ! Establish new communicator via port name
    ! ----------------------------------------
    call MPI_COMM_CONNECT(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, newcomm, ierr)
    if (tc_nml%verbosity > 0) then
       write(6,'(a,i0)') 'Established new communicator:', newcomm
       call flush(6)
    end if

    ! --------------------------------
    ! Send job information to terachem
    ! --------------------------------
    dbuffer(:,:) = ''
    if ( .not. tc_nml%use_template ) then
       ! Send namelist data
       irow = 1
       write(dbuffer(:,irow),'(a,/,a)') 'basis',      tc_nml%basis
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,a)') 'method',     tc_nml%method
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,a)') 'dftd',       tc_nml%dftd
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,a)') 'guess',       tc_nml%guess
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,a)') 'precision',  tc_nml%precision
       irow = irow + 1
       if( do_grad ) then
          write(dbuffer(:,irow),'(a,/,a)') 'run', 'gradient'
       else
          write(dbuffer(:,irow),'(a,/,a)') 'run', 'energy'
       end if
       ! This is set to 1 because this instructs terachem to skip 
       ! calculating the self-energy of the charges
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,a)') 'amber', 'yes'
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,E22.16)') 'threall', tc_nml%threall
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,E22.16)') 'convthre', tc_nml%convthre
       
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,i0)') 'charge', tc_nml%charge
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,i0)') 'spinmult', tc_nml%spinmult
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,i0)') 'maxit', tc_nml%maxit
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,i0)') 'dftgrid', tc_nml%dftgrid
       irow = irow + 1
       write(dbuffer(:,irow),'(a,/,a)') 'cis',        tc_nml%cis
       if ( trim(tc_nml%cis) == 'yes' ) then
          irow = irow + 1
          write(dbuffer(:,irow),'(a,/,i0)') 'cisnumstates', tc_nml%cisnumstates
          irow = irow + 1
          write(dbuffer(:,irow),'(a,/,i0)') 'cistarget', tc_nml%cistarget
       end if

       ! Write gpus
       if ( tc_nml%ngpus > 0 ) then
          irow = irow + 1
          write(dbuffer(:,irow), '(a,/,9999(i3))') 'gpus      ', tc_nml%ngpus, (tc_nml%gpuids(i), i = 1, tc_nml%ngpus)
       end if
       ! Finish writing - send 'end'
       irow = irow + 1
       write(dbuffer(:,irow),'(a)') 'end', ''
    else
       call read_template(tplfile, dbuffer)
    end if

    if (tc_nml%verbosity > 1) then
       write(6,'(a)') '(debug) sending namelist data:'
       do j=1, 32
          write(6,*) trim(dbuffer(1,j)), ' = ', trim(dbuffer(2,j))
       end do
       call flush(6)
    end if

    call MPI_Send( dbuffer, 2*clen*size(dbuffer,2), MPI_CHARACTER, 0, 2, newcomm, ierr )

    ! -----------------------------------------
    ! Send nqmatoms and the type of each qmatom
    ! TeraChem needs this information to correctly initialize the GPUs
    ! (depending on whether d orbitals are in use or not)
    ! -----------------------------------------
    if (tc_nml%verbosity > 1) then
       write(6,'(/, a, i0)') 'Sending nqmatoms = ', nqmatoms
       call flush(6)
    end if
    call MPI_Send( nqmatoms, 1, MPI_INTEGER, 0, 2, newcomm, ierr )

    if (tc_nml%verbosity > 1) then
       write(6,'(/,a)') 'Sending QM atom types: '
       do i = 1, nqmatoms
          write(6,'(a)') atom_types(i)
          call flush(6)
       end do
    end if
    call MPI_Send( atom_types, 2*size(atom_types), MPI_CHARACTER, 0, 2, newcomm, ierr )


    if (tc_nml%verbosity > 0) then
       write(6,'(a)') '<<<<< Leaving connect_to_terachem()'
       call flush(6)
    end if

  end subroutine connect_to_terachem

#endif

  ! ---------------------------------
  ! Write the input file for TeraChem
  ! ---------------------------------
  subroutine write_inpfile( inpfile, crdfile, chrgfile, tplfile, nqmatoms, qmcoords,&
    nclatoms, clcoords, tc_nml, do_grad )

    use qmmm_module, only : qmmm_struct
    use ElementOrbitalIndex, only : elementSymbol

    implicit none

    character (len=*) , intent(in) :: inpfile
    character (len=*) , intent(in) :: crdfile
    character (len=*) , intent(in) :: chrgfile
    character (len=*) , intent(in) :: tplfile
    integer           , intent(in) :: nqmatoms
    _REAL_            , intent(in) :: qmcoords(:,:)
    integer           , intent(in) :: nclatoms
    _REAL_            , intent(in) :: clcoords(:,:)
    type(tc_nml_type) , intent(in) :: tc_nml
    logical, intent(in) :: do_grad

    character(len=80)  :: read_buffer
    character(len=20)  :: bakfile
    integer            :: i, j, ios
    integer, parameter :: iurun = 10
    logical, save :: first_call = .true.

    if(tc_nml%verbosity > 0) then
      write(6,*)'Writing TeraChem inpfile...'//inpfile
    end if

    if ( tc_nml%use_template ) then
      bakfile = tplfile//'.bak'
      if ( first_call ) then
        call copy_template(tplfile, trim(bakfile))
      end if
      call system('cp '//trim(bakfile)//' '//inpfile)
    end if

    open(iurun, file=inpfile, position='append', iostat=ios)
    if ( ios > 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_tc_module)', &
           'Error opening TeraChem input file '//inpfile//' for writing', &
           'Will quit now')
    end if

    if ( .not. tc_nml%use_template ) then ! We only need the last keys if using template file
      write(iurun, '(a,/,a)')'# Run using SANDER file-based interface for TeraChem','#'
      write(iurun, '(2a)')       'basis        ', trim(tc_nml%basis)
      write(iurun, '(2a)')       'method       ', trim(tc_nml%method)
      write(iurun, '(2a)')       'precision    ', trim(tc_nml%precision)
      write(iurun, '(a,E22.16)') 'threall      ', tc_nml%threall
      write(iurun, '(a,E22.16)') 'convthre     ', tc_nml%convthre
      write(iurun, '(2a)')       'dftd         ', trim(tc_nml%dftd)
      write(iurun, '(a,i3)')     'charge       ', tc_nml%charge
      write(iurun, '(a,i3)')     'spinmult     ', tc_nml%spinmult
      write(iurun, '(a,i4)')     'maxit        ', tc_nml%maxit
      write(iurun, '(a,i3)')     'dftgrid      ', tc_nml%dftgrid
      write(iurun, '(2a)')       'cis          ', trim(tc_nml%cis)
      if ( trim(tc_nml%cis) == 'yes' ) then
         write(iurun, '(a,i3)')     'cisnumstates ', tc_nml%cisnumstates
         write(iurun, '(a,i3)')     'cistarget    ', tc_nml%cistarget
      end if
      
      if ( tc_nml%ngpus > 0 ) then
         write(iurun, '(a,65(i3))')    'gpus      ', tc_nml%ngpus, (tc_nml%gpuids(i), i = 1, tc_nml%ngpus)
      end if
    end if

    if( .not. first_call ) then
       write(iurun, '(2a)') 'guess       ', tc_nml%guess
    end if

    if( do_grad ) then
      write(iurun, '(a)') 'run         gradient'
    else
      write(iurun, '(a)') 'run         energy'
    end if
    write(iurun, '(a)') 'coordinates '//crdfile

    if ( nclatoms > 0 ) then
       write(iurun, '(a)') 'pointcharges '//chrgfile
       write(iurun, '(a)') 'amber       yes'
    end if

    write(iurun, '(a)') 'end'

    close(iurun)
    
    ! Now write xyz file with geometry
    open(iurun, file=crdfile, iostat=ios)
    if ( ios > 0 ) then
      call sander_bomb('write_inpfile (qm2_extern_tc_module)', &
           'Error opening coordinate file '//crdfile//' for writing', &
           'Will quit now')
    end if
    write(iurun,'(i5,/)') nqmatoms
    do i = 1, nqmatoms
       write(iurun,'(a2,1x,3f21.16)') &
          elementSymbol(qmmm_struct%iqm_atomic_numbers(i)), &
          qmcoords(:,i)
    enddo
    close(iurun)

    ! Now write xyz file with point charges
    if ( nclatoms > 0 ) then
       open(iurun, file=chrgfile, iostat=ios)
       if ( ios > 0 ) then
          call sander_bomb('write_inpfile (qm2_extern_tc_module)', &
               'Error opening point charge coordinate file '//chrgfile//' for writing', &
               'Will quit now')
       end if
       write(iurun,'(i6,/)') nclatoms
       do i = 1, nclatoms
          write(iurun,'(4f21.16)') clcoords(4,i), clcoords(:3,i)
       enddo
       close(iurun)
    end if

    if (first_call) then
       first_call = .false.
    end if

  end subroutine write_inpfile

  ! ---------------------
  ! Read TeraChem results
  ! ---------------------
  subroutine read_results( datfile, nqmatoms, nclatoms, escf, dxyzqm, dxyzcl,&
       dipmom, do_grad )

    implicit none

    character(len=*), intent(in) :: datfile
    integer, intent(in) :: nqmatoms, nclatoms
    _REAL_, intent(out) :: escf, dxyzqm(3,nqmatoms), dxyzcl(3,nclatoms)
    logical, intent(in) :: do_grad

    integer :: ios, i, itmp
    integer, parameter :: iunit = 351
    _REAL_, intent(out) :: dipmom(4,3)
    character(len=100) :: read_buffer

    open(iunit, file=datfile, iostat=ios)
    if ( ios /= 0 ) then
       call sander_bomb('read_results (qm2_extern_tc_module)', &
            'Error opening TeraChem output file '//datfile//' for reading', &
            'Will quit now')
    end if

    do

       read (iunit, '(a)', iostat = ios) read_buffer
       ! End of file; data not found
       if (ios < 0) then
          call sander_bomb('read_results (qm2_extern_tc_module)', &
               'Error reading TeraChem output from file '//datfile, &
               '("FINAL ENERGY" or "dE/dX" or "MM / Point charge part" not found)')
       end if

       if ( read_buffer(1:14) == 'FINAL ENERGY: ' ) then
          ! Read energy
          itmp = len_trim(read_buffer)
          read (read_buffer(15:itmp-5),*) escf
          if( .not. do_grad ) exit
       else if ( read_buffer(9:13) == 'dE/dX' ) then
          ! Read QM gradient data
          do i = 1, nqmatoms
             read (iunit, '(3(f16.10,1x))') dxyzqm(:,i)
          end do
          if ( nclatoms == 0 ) exit
       else if ( read_buffer(9:31) == 'MM / Point charge part' ) then
          ! Read MM atoms gradient data
          do i = 1, nclatoms
             read (iunit, '(3(f16.10,1x))') dxyzcl(:,i)
          end do
          exit
       end if

       if (index(read_buffer,'QM  DIPOLE MOMENT: ') > 0) then
         read(read_buffer(index(read_buffer,'{')+1:index(read_buffer,'}')-1),*) dipmom(1:3,1)
         read(read_buffer(index(read_buffer,'|D| = ')+5:index(read_buffer,') DEBYE')-1),*) dipmom(4,1)
       end if
       if (index(read_buffer,'MM  DIPOLE MOMENT: ') > 0) then
         read(read_buffer(index(read_buffer,'{')+1:index(read_buffer,'}')-1),*) dipmom(1:3,2)
         read(read_buffer(index(read_buffer,'|D| = ')+5:index(read_buffer,') DEBYE')-1),*) dipmom(4,2)
       end if
       if (index(read_buffer,'TOT DIPOLE MOMENT: ') > 0) then
         read(read_buffer(index(read_buffer,'{')+1:index(read_buffer,'}')-1),*) dipmom(1:3,3)
         read(read_buffer(index(read_buffer,'|D| = ')+5:index(read_buffer,') DEBYE')-1),*) dipmom(4,3)
       end if

    end do
    close(iunit)

  end subroutine read_results

  ! ----------------------
  ! Write TeraChem results
  ! ----------------------
  subroutine write_results( dipfile, dipmom )

    use constants, only: CODATA08_AU_TO_DEBYE
    implicit none

    character(len=*), intent(in) :: dipfile
    _REAL_, intent(in)  :: dipmom(4,3)
    integer :: iunit = 351, ios
    logical, save :: first_call = .true.

    ! write dipole moment to dipole moment property file
    open (iunit, file=dipfile, position='append', iostat=ios)
    if ( ios /= 0 ) then
       call sander_bomb('write_results (qm2_extern_tc_module)', &
            'Error opening file '//dipfile//' for appending.', &
            'Will quit now')
    end if
    if(first_call) then
      first_call=.false.
      write(iunit,'(a, f22.12)') "Using DEBYE_TO_AU = ", 1.0d0 / CODATA08_AU_TO_DEBYE      
      write(iunit,'(a)') "| Dipole moment (a.u.-Å): QM {x, y, z}, |D|; MM {x y z}, |D|; TOT { x y z}, |D|"
    end if
    write(iunit,'(12f10.6)') dipmom(:,1)/CODATA08_AU_TO_DEBYE, &
         dipmom(:,2)/CODATA08_AU_TO_DEBYE, &
         dipmom(:,3)/CODATA08_AU_TO_DEBYE
    close(iunit)

  end subroutine write_results

  subroutine tc_finalize()

    implicit none
#ifdef MPI
    include 'mpif.h'

    integer :: ierr
    _REAL_  :: empty
    if(do_mpi) then
      call MPI_Send( empty, 1, MPI_DOUBLE_PRECISION, 0, 0, newcomm, ierr )
    end if
#endif

  end subroutine tc_finalize
      

  subroutine copy_template( tplfile, bakfile )

    use UtilitiesModule, only: Upcase
    implicit none
    character(len=*), intent(in) :: tplfile, bakfile
    integer, parameter :: tplunit = 351, bakunit=352
    character(len=100) :: read_buffer
    integer :: tplerr, bakerr, ios
    logical :: end_key = .false.

    open(tplunit, file=tplfile, iostat=tplerr )
    open(bakunit, file=bakfile, iostat=bakerr )
    if ( tplerr /= 0 ) then
      call sander_bomb('copy_template (qm2_extern_adf_module)', &
        'Error opening ADF template file '//tplfile//' for reading', &
        'Will quit now.')
    end if
    if ( bakerr /= 0 ) then
      call sander_bomb('copy_template (qm2_extern_adf_module)', &
        'Error opening ADF template backup file '//bakfile//' for writing', &
        'Will quit now.')
    end if

    ! Write tplfile (without end key) to bakfile
    do
       read (tplunit, '(a)', iostat = ios) read_buffer
       ! End of file; stop writing
       if (ios < 0) then
         exit
       end if
       if ( index(Upcase(read_buffer), 'END') > 0 ) then
          exit
       end if
       write(bakunit, '(a)') trim(read_buffer)
    end do

    close(tplunit)
    close(bakunit)

  end subroutine copy_template

  subroutine read_template( tplfile, dbuffer )

    use UtilitiesModule, only: Upcase
    implicit none
    character(len=*), intent(in) :: tplfile
    character(len=*), intent(out):: dbuffer(:,:)
    integer, parameter :: tplunit = 351
    character(len=100) :: read_buffer
    integer :: tplerr, ios, i, ispace, j

    dbuffer=''
    open(tplunit, file=tplfile, iostat=tplerr )
    if ( tplerr /= 0 ) then
      call sander_bomb('copy_template (qm2_extern_adf_module)', &
        'Error opening ADF template file '//tplfile//' for reading', &
        'Will quit now.')
    end if
    ! Write tplfile (without end key) to dbuffer
    i=1
    do
       read (tplunit, '(a)', iostat = ios) read_buffer
       ! End of file; stop reading
       if (ios < 0) then
         exit
       end if
       if ( index(Upcase(read_buffer), 'END') > 0 ) then
          exit
       end if
       ispace=index(read_buffer,' ')
       write(dbuffer(:,i), '(a,/,a)') read_buffer(1:ispace-1), adjustl(read_buffer(ispace:))
       i=i+1
    end do
    write(dbuffer(:,i),   '(a,/,a)') 'amber','yes'
    write(dbuffer(:,i+1), '(a,/,a)') 'run', 'gradient'
    write(dbuffer(:,i+2), '(a,/,a)') 'end',''

    close(tplunit)

  end subroutine read_template

end module qm2_extern_tc_module
