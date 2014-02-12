#include "dprec.fh"
module qm2_extern_gau_module
! ----------------------------------------------------------------
! Interface for Gaussian based QM MD 
!
! Currently supports:
! pure QM
!
! Initial implementation by
! Matthew Clark
! under supervision of
! Andreas Goetz and Ross Walker (SDSC)
! 
! Date: February 2011
!
! ----------------------------------------------------------------

  implicit none

  private
  public :: get_gau_forces

  type gau_nml_type
     character(len=20) :: method
     character(len=20) :: basis
     integer :: scf_conv
     integer :: charge
     integer :: spinmult
     integer :: ntpr
     integer :: num_threads 
     integer :: verbosity
     logical :: dipole
     logical :: use_template
  end type gau_nml_type

contains

  ! --------------------------------------------
  ! Get QM energy and forces from Gaussian
  ! --------------------------------------------
  subroutine get_gau_forces( do_grad, nstep, ntpr_default, id, nqmatoms, qmcoords,&
    nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge)

    use constants, only: CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO
    use file_io_dat

    implicit none

    logical, intent(in) :: do_grad              ! Return gradient/not
    integer, intent(in) :: nstep                ! MD step number
    integer, intent(in) :: ntpr_default         ! frequency of printing
    character(len=3), intent(in) :: id          ! ID number for PIMD or REMD
    integer, intent(in) :: nqmatoms             ! Number of QM atoms
    _REAL_,  intent(in) :: qmcoords(3,nqmatoms) ! QM atom coordinates
    integer, intent(in) :: nclatoms             ! Number of MM atoms
    _REAL_,  intent(in) :: clcoords(4,nclatoms) ! MM atom coordinates and charges in au
    _REAL_, intent(out) :: escf                 ! SCF energy
    _REAL_, intent(out) :: dxyzqm(3,nqmatoms)   ! SCF QM force
    _REAL_, intent(out) :: dxyzcl(3,nclatoms)   ! SCF MM force
    _REAL_              :: dipxyz(3), dipole    ! Dipole moment
    integer, intent(out) :: charge

    type(gau_nml_type), save     :: gau_nml
    logical, save                :: first_call = .true.
    integer                      :: i
    integer                      :: printed =-1 ! Used to tell if we have printed this step yet 
                                                ! since the same step may be called multiple times
    character (len=150)          :: call_buffer
    character(len=3), save       :: program
    character(len=*), parameter  :: basename = 'gau_job'
    character(len=*), parameter  :: inpext = '.inp'
    character(len=*), parameter  :: logext = '.log'
    character(len=*), parameter  :: rstext = '.chk' ! Restart from checkpoint files
    character(len=*), parameter  :: dipext = '.dip'
    character(len=*), parameter  :: tplext = '.tpl'
    character(len=14)            :: inpfile, rstfile, logfile, fortfile, dipfile, tplfile
    ! Need to prepend subdirectory if doing REMD, PIMD
    character(len=25)            :: subdir 

    ! for system call
    integer :: system
    integer :: stat

    ! assemble input - / output data filenames
    inpfile = basename//trim(id)//inpext
    logfile = basename//trim(id)//logext
    rstfile = basename//trim(id)//rstext
    dipfile = basename//trim(id)//dipext
    tplfile = basename//tplext
    fortfile='fort.7'

    ! Setup on first call
    if ( first_call ) then
       first_call = .false.
       write (6,'(/,a,/)') '  >>> Running calculations with Gaussian <<<'
       call get_namelist( ntpr_default, gau_nml )
       call print_namelist( gau_nml ) 
       ! Check for version of Gaussian to use; store as 'program'
       call check_installation( program, id )
       charge = gau_nml%charge
       write (6,'(80a)') ('-', i=1,80)
       write (6,'(a)') '   4.  RESULTS'
       write (6,'(80a)') ('-', i=1,80)
       ! Remove old inpfile, logfile, fort.7, and dipfile at the 
       ! beginning of a run so only the latest run is stored.
       stat = system('rm -f '//inpfile//' '//logfile//' '//rstfile//' '//fortfile//' '//dipfile)
       if ( stat /= 0 ) then
          call sander_bomb('get_gau_forces (qm2_extern_gau_module)', & 
               'Error with system call (removing files)', &
               'Will quit now.')
       end if
    end if
    
    call write_inpfile( trim(inpfile), trim(rstfile), trim(tplfile), &
      nqmatoms, qmcoords, nclatoms, clcoords, gau_nml, do_grad )

    if ( gau_nml%verbosity > 0 ) then
      write(6,'(a)') 'Runfile has written successfully; Calling Gaussian...'
    end if

    ! Run g09/g03
    ! Separate runs into different directories if we are doing PIMD
    subdir=''
    call_buffer=''
    if(trim(id)/='') then 
      subdir='./'//trim(id)//'/'
      call_buffer=' mkdir -p '//trim(subdir)//'; cd '//trim(subdir)//'; mv ../'//inpfile//' .;'
    end if
    stat = system(trim(call_buffer)//' '//program//' '//inpfile) 
    if ( stat /= 0 ) then
       call sander_bomb('get_gau_forces (qm2_extern_gau_module)', & 
            'Error with system call (executing Gaussian)', &
            'Will quit now.')
    end if

    if ( gau_nml%verbosity > 0 ) then    
      write(6,'(a)') 'Gaussian execution success; Processing Gaussian results...'
    end if    

    ! Call read_results - retrieve data from Gaussian .log and .fort7 files
    ! Will output data to escf and dxyqm for pure QM or MM runs
    ! For QM/MM runs will also return dxyzcl containing electric field strength
    
    ! Search in subdir for logfile and fortfile if doing PIMD
    ! Otherwise, search current directory
    call read_results( trim(subdir)//trim(logfile), trim(subdir)//trim(fortfile),&
      nqmatoms, escf, dxyzqm, nclatoms, dxyzcl, dipxyz, dipole, do_grad )

    ! Call write_results with dipfile, dipxyz, and magnitude of dipxyz;
    ! will write output to .dip and .chg files
    if ( gau_nml%ntpr > 0 .and. mod(nstep, gau_nml%ntpr) == 0 ) then
      if ( printed /= nstep .and. gau_nml%dipole ) then
        call write_results( trim(dipfile), dipxyz, dipole )
        printed = nstep
      end if
    end if

    ! Save copy of last input and log files and delete old fort.7 file
    stat = system('mv '//trim(subdir)//inpfile//' '//trim(subdir)//'old.'//inpfile)
    stat = stat + system('mv '//trim(subdir)//logfile//' '//trim(subdir)//'old.'//logfile)
    stat = stat + system('rm -f '//trim(subdir)//fortfile)
    if ( stat /= 0 ) then
       call sander_bomb('get_gau_forces (qm2_extern_gau_module)', & 
            'Error with system call (moving / removing files)', &
            'Will quit now.')
    end if
    

    ! F = E*q to get gradients
    ! Note dxyzcl is currently holding the electric field strength
    do i = 1, nclatoms
        dxyzcl(:,i)= -dxyzcl(:,i) * clcoords(4,i)
    end do

    ! Convert gradient from au to kcal/(mol*A)
    if ( do_grad ) then
      dxyzqm(:,:) = dxyzqm(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
      if ( nclatoms > 0 ) then
          dxyzcl(:,:) = dxyzcl(:,:) * CODATA08_AU_TO_KCAL * CODATA08_A_TO_BOHRS
      end if
    else
      dxyzqm = ZERO
      if ( nclatoms > 0 ) dxyzcl = ZERO
    end if
    
    if( gau_nml%verbosity > 0 ) then
      write (6,'(a)') 'get_gau_forces - final gradient:'
      do i = 1, nqmatoms
        write(6,*) dxyzqm(1:3,i)
      end do
      write(6,*)
      if ( nclatoms > 0 ) then
         write(6,'(a)') 'MM region:'
         do i = 1, nclatoms
            write(6,'(3(x,f16.10))') dxyzcl(1:3, i)
         end do
      end if
    end if

    escf = escf * CODATA08_AU_TO_KCAL

  end subroutine get_gau_forces

  ! ---------------------------------------------
  ! Read Gaussian namelist values from file mdin,
  ! use default values if none are present.
  ! ---------------------------------------------
    
 
  subroutine get_namelist(ntpr_default, gau_nml)

    implicit none
    integer, intent(in) :: ntpr_default
    type(gau_nml_type), intent(out) :: gau_nml

    character(len=20) :: method, basis
    integer :: charge, spinmult, verbosity
    integer :: scf_conv, ntpr, num_threads, dipole, use_template
    namelist /gau/ method, basis, scf_conv, ntpr, charge, spinmult, &
      num_threads, verbosity, dipole, use_template

    integer :: ifind, ierr

    ! Set default values for gau namelist values
    method       = 'BLYP'
    basis        = '6-31G*'
    scf_conv     = 8 
    charge       = 0
    spinmult     = 1
    ntpr         = ntpr_default
    num_threads  = 1
    verbosity    = 0
    dipole       = 0
    use_template = 0

    ! Read namelist
    rewind 5
    read(5,nml=gau,iostat=ierr)

    if ( ierr > 0 ) then
       call sander_bomb('get_namelist (qm2_extern_gau_module)', &
            '&gau namelist read error', &
            'Please check your input.')
    else if ( ierr < 0 ) then
       write(6,'(a/a)') '&gau namelist read encountered end of file', &
            'Please check your input if the calculation encounters a problem'
    end if

    ! Assign namelist values to gau_nml data type
    gau_nml%method       = method
    gau_nml%basis        = basis
    gau_nml%scf_conv     = scf_conv
    gau_nml%charge       = charge
    gau_nml%spinmult     = spinmult
    gau_nml%ntpr         = ntpr
    gau_nml%num_threads  = num_threads
    gau_nml%verbosity    = verbosity
    if ( dipole == 0 ) then
       gau_nml%dipole = .false.
    else if ( dipole == 1 ) then
       gau_nml%dipole = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_gau_module)', &
            '&gau dipole value not allowed', &
            'Please check your input. dipole can only be 0 or 1.')
    end if

    if ( use_template == 0 ) then
       gau_nml%use_template = .false.
    else if ( use_template == 1 ) then
       gau_nml%use_template = .true.
    else
       call sander_bomb('get_namelist (qm2_extern_gau_module)', &
            '&gau use_template value not allowed', &
            'Please check your input. use_template can only be 0 or 1.')
    end if

  end subroutine get_namelist

  ! --------------------------------
  ! Print Gaussian namelist settings
  ! --------------------------------
  subroutine print_namelist(gau_nml)

    implicit none
    type(gau_nml_type), intent(in) :: gau_nml

    write(6, '(a)')       '| &gau'
    write(6, '(2a)')      '|   method       = ', gau_nml%method
    write(6, '(2a)')      '|   basis        = ', gau_nml%basis
    write(6, '(a,i2)')    '|   scf_conv     = ', gau_nml%scf_conv
    write(6, '(a,i2)')    '|   charge       = ', gau_nml%charge
    write(6, '(a,i2)')    '|   spinmult     = ', gau_nml%spinmult
    write(6, '(a,i0)')    '|   ntpr         = ', gau_nml%ntpr
    write(6, '(a,i2)')    '|   num_threads  = ', gau_nml%num_threads
    write(6, '(a,i2)')    '|   verbosity    = ', gau_nml%verbosity
    write(6, '(a,l)')     '|   dipole       = ', gau_nml%dipole
    write(6, '(a,l)')     '|   use_templte  = ', gau_nml%use_template
    write(6,'(a)')        '| /'

  end subroutine print_namelist


  ! --------------------------------------------
  ! Check whether Gaussian is properly installed
  ! Also, return version of Gaussian to use
  ! --------------------------------------------
  subroutine check_installation( program, id )

    implicit none

    character(len=3), intent(out) :: program
    character(len=3), intent(in) :: id

    character(len=80) :: read_buffer
    character(len=80) :: call_buffer
    character(len=80) :: filename
    integer :: iunit = 77
    integer :: stat
    integer :: system

    filename = 'extern_location'//trim(id)
    
    ! Search for g03 or g09 executable
    call_buffer = 'which g09 > '//trim(filename)
    stat = system(trim(call_buffer))
    if ( stat == 0 ) then
       !using g09:      
       program='g09'
    else
       call_buffer = 'which g03 > '//trim(filename)
       stat = system(trim(call_buffer))
       if ( stat == 0 ) then
          !using g03:
          program='g03'
       else
          call sander_bomb('check_installation (qm2_extern_gau_module)', & 
               'Executable g03 or g09 not found', &
               'Please check your Gaussian installation')
       end if
    end if

    write(6,'(3a,/)') '| Program ', program, ' found!'

    ! Get complete executable path
    open (unit=iunit, file=trim(filename), form='formatted', iostat=stat)
    if ( stat /= 0 ) then
       call sander_bomb('check_installation (qm2_extern_gau_module)', & 
            'Internal error opening file with path location of executable.', &
            'Quitting now.')
    end if
    read (iunit, '(a)', iostat=stat) read_buffer
    if ( stat /= 0 ) then
       call sander_bomb('check_installation (qm2_extern_gau_module)', & 
            'Internal error reading from file with path location of executable.', &
            'Quitting now.')
    end if
    close (unit=iunit, status='delete', iostat=stat)
    if ( stat /= 0 ) then
       call sander_bomb('check_installation (qm2_extern_gau_module)', & 
            'Internal error closing and deleting file with path location of executable.', &
            'Quitting now.')
    end if

    write(6,'(2a,/)') '| Executable location: ', trim(read_buffer)

  end subroutine check_installation

  ! -----------------------------
  ! Write input file for Gaussian
  ! -----------------------------

  subroutine write_inpfile( inpfile, rstfile, tplfile, nqmatoms, qmcoords,&
    nclatoms, clcoords, gau_nml, do_grad )

    use qmmm_module, only : qmmm_struct
    use ElementOrbitalIndex, only : elementSymbol

    implicit none

    character(len=*), intent(in)   :: inpfile, rstfile, tplfile
    integer, intent(in)            :: nqmatoms
    _REAL_,  intent(in)            :: qmcoords(:,:)
    integer, intent(in)            :: nclatoms
    _REAL_,  intent(in)            :: clcoords(:,:)
    type(gau_nml_type), intent(in) :: gau_nml
    logical, intent(in)            :: do_grad

    integer, parameter :: iunit = 351, tplunit = 352
    integer            :: i, j, ierr
    integer            :: tplerr, ios
    character(len=100) :: route
    character(len=100) :: read_buffer
    logical, save      :: first_call = .true.

    if( gau_nml%verbosity > 0 ) then
      write(6,'(a)') 'Writing Gaussian inpfile '//inpfile
    end if

    open(iunit, file=inpfile, iostat=ierr)
    if ( ierr /= 0 ) then
       call sander_bomb('write_inpfile (qm2_extern_gau_module)', &
            'Error opening Gaussian inpfile '//inpfile//' for writing', &
            'Will quit now.')
    end if

    ! Write link option
    ! chkfile so we can restart from previous runs
    write (iunit,'(a)') '%chk='//rstfile

    ! Shared memory parallelism
    write (iunit, '(a,i0)') '%NProcShared=',gau_nml%num_threads


    ! Assemble route
    ! Force calculation, punch gradient, never use symmetry
    ! If using template, write route information in template file to route card
    route=''
    if( gau_nml%use_template) then
      open(tplunit, file=tplfile, iostat=tplerr)
      if ( tplerr /= 0 ) then
        call sander_bomb('write_inpfile (qm2_extern_gau_module)', &
          'Error opening Gaussian template file '//tplfile//' for reading', &
          'Will quit now.')
      end if
      read (tplunit, '(a)', iostat = ios) read_buffer
      if (ios < 0) then
        call sander_bomb('write_inpfile (qm2_extern_gau_module)', &
          'Error reading Gaussian template file '//tplfile, &
          'Will quit now.')
      end if
      write(route,'(a)') trim(route)//' '//trim(read_buffer)//' '

    else
      ! Method/Basis
      route = '#P '//trim(gau_nml%method)//'/'//trim(gau_nml%basis)
      ! SCF convergence setting
      write(route,'(a,i0,a)') trim(route)//' SCF=(Conver=', gau_nml%scf_conv, ')'
    end if
    ! Always need these options
    if ( do_grad ) then
      route = trim(route)//' Force Punch=Derivatives NoSymm'
    else
      ! Don't calculate gradient anymore (punch needed for energy) 
      route = trim(route)//' Punch=Derivatives NoSymm'
    end if
    ! If we are not on our first run, restart from .chk file 
    if( .not. first_call ) then
      route = trim(route)//' Guess=Read'
    end if
    ! If doing electrostatic embadding QM/MM, 
    ! read external point charges, print to .log
    if( nclatoms > 0 ) then
      route = trim(route)//' Charge Prop=(Field,Read)'
    end if

    write(iunit,'(a)') trim(route)
    write(iunit,'(a)')


    ! Write optional comment, net charge and spin multiplicity
    if( gau_nml%use_template) then
      do i=1, 2
        read (tplunit, '(a)', iostat = ios) read_buffer
        if (ios < 0) then
          call sander_bomb('write_inpfile (qm2_extern_gau_module)', &
            'Error reading Gaussian template file '//tplfile, &
            'Will quit now.')
        end if
        write(iunit,'(a)') trim(read_buffer)
        if(i==1) then
          write(iunit,'(a)')
        end if
      end do
      close(tplunit)
    else
      write (iunit,'(a,/)') 'Gaussian run using SANDER external interface.'
      write (iunit,'(i0,a,i0)') gau_nml%charge,' ', gau_nml%spinmult
    end if


    ! Write QM atoms and coordinates
    do i = 1, nqmatoms
      write(iunit,'(a2,1x,3f25.16)') elementSymbol(qmmm_struct%iqm_atomic_numbers(i)), qmcoords(1:3,i)
    end do
    write(iunit,'()')

    ! When electrostatic embadding QM/MM is in use 
    ! write MM coordinates with point charges
    if ( nclatoms > 0 ) then
      do i = 1, nclatoms
        write(iunit,'(4f21.16)') clcoords(:,i)
      end do
      write(iunit,'()')

      ! Write a second time without charges
      do i = 1, nclatoms
        write(iunit,'(4f21.16)') clcoords(:3,i)
      end do

      close(iunit)
    end if

   write(iunit,'()')
   close(iunit, iostat=ierr)

   if ( ierr /= 0 ) then
     call sander_bomb('write_inpfile (qm2_extern_gau_module)', &
       'Error closing Gaussian runfile after writing', &
       'Will quit now.')
   end if
   first_call = .false.

  end subroutine write_inpfile

  subroutine read_results( datfile, fortfile, nqmatoms, escf, dxyzqm,&
    nclatoms, dxyzcl, dipxyz, dipole, do_grad )

    implicit none

    character(len=*), intent(in) :: datfile, fortfile
    integer, intent(in)          :: nqmatoms, nclatoms
    _REAL_, intent(out)          :: escf, dxyzqm(3,nqmatoms), & 
                                    dxyzcl(3,nclatoms) ! dxyzcl will return containing the electric field at x,y,z
    _REAL_, intent(out)          :: dipxyz(3), dipole
    logical, intent(in)          :: do_grad


    _REAL_ :: self_energy = 0 ! Temporary variable to hold self energy of point charges
    integer :: ios, i
    integer, parameter :: iunit = 351
    character(len=120) :: read_buffer

    open(iunit, file=datfile, status='old', iostat=ios)
    if ( ios /= 0 ) then
      call sander_bomb('read_results (qm2_extern_gau_module)', &
        'Error opening Gaussian log file '//datfile//' (expected in same dir as input file).', &
        'Will quit now')
    end if

    do
      read (iunit, '(a)', iostat = ios) read_buffer
      ! End of file; nothing left to read
      if (ios < 0) then
        exit
      end if 
    
      ! Store SCF energy to escf 
      if ( read_buffer(1:10) == ' SCF Done:' ) then
        ! Read value after the first "=" sign
        read(read_buffer(index(read_buffer,'=')+1:),*) escf
      end if
       
      ! Energy is stored under EUMP2 when running MP2 calculations
      ! If found, write this value over our old escf
      if ( read_buffer(1:5) == ' E2 =' ) then
        ! Read first value after the second "=" sign
        read(read_buffer(index(read_buffer,'=',.true.)+1:),*) escf 
      end if

      ! If doing QM/MM...
      if ( nclatoms > 1) then
        ! Retrieve self energy of charges (must add to escf)
        if (read_buffer(1:30) == ' Self energy of the charges =') then
          read(read_buffer(index(read_buffer,'=')+1:),*) self_energy
        end if

        if( do_grad ) then
          ! Read efield at CL atoms into dxyzcl
          if( read_buffer(1:10) == '    1 Atom') then
            ! Skip QM atoms
            do  i = 1, nqmatoms-1
              read(iunit, '(a)') read_buffer
            end do
            ! Read into dxyzcl
            do i = 1, nclatoms
              read (iunit, '(a)', iostat = ios) read_buffer
              read(read_buffer(25:),*) dxyzcl(:,i ) 
            end do
          end if
        end if
      end if

      ! Read dipole moment to vars dipxyz and dipole
      if (read_buffer(1:15) == ' Dipole moment ') then
        ! Read next line in
        read(iunit, '(a)') read_buffer
        ! Store dipole moment and magnitude
        read(read_buffer(index(read_buffer,'X=')+2:index(read_buffer,'Y=')-2),*),   dipxyz(1)
        read(read_buffer(index(read_buffer,'Y=')+2:index(read_buffer,'Z=')-2),*),   dipxyz(2) 
        read(read_buffer(index(read_buffer,'Z=')+2:index(read_buffer,'Tot=')-2),*), dipxyz(3) 
        read(read_buffer(index(read_buffer,'Tot=')+5:),*) dipole
      end if

    end do
    
    escf = escf - self_energy
    close(iunit)
     
    ! Read gradients from the fort.7 file
    open(iunit, file=fortfile) !Will always exist
    do i = 1, nqmatoms
      read (iunit, '(3f20.16)', iostat=ios) dxyzqm(:,i)

      if (ios < 0) then
        call sander_bomb('read_results (qm2_extern_gau_module)', &
          'Error reading file '//fortfile, &
          'Will quit now')
      end if
    end do
 
  end subroutine read_results

  ! ----------------------
  ! Write Gaussian results
  ! ----------------------
  subroutine write_results( dipfile, dipxyz, dipole )

    use constants, only: CODATA08_AU_TO_DEBYE
    implicit none

    character(len=*), intent(in) :: dipfile
    _REAL_  :: dipxyz(3), dipole
    integer :: iunit = 351, ios
    logical, save :: first_call = .true.

    ! write dipole moment to dipole moment property file
    open (iunit, file=dipfile, position='append', iostat=ios)
    if ( ios /= 0 ) then
       call sander_bomb('write_results (qm2_extern_gau_module)', &
            'Error opening file '//dipfile//' for appending.', &
            'Will quit now')
    end if
    if(first_call) then
      first_call=.false.
      write(iunit,'(a, f22.12)') "Using DEBYE_TO_AU = ", 1.0d0 / CODATA08_AU_TO_DEBYE 
      write(iunit,'(a)') "| Dipole moment (a.u.-Å): {x, y, z}, |D|"
    end if
    write(iunit,'(4f15.6)') dipxyz(:)/CODATA08_AU_TO_DEBYE, dipole/CODATA08_AU_TO_DEBYE
    close(iunit)

  end subroutine write_results

end module qm2_extern_gau_module
