#include "dprec.fh"
module qm2_extern_nw_module
! ----------------------------------------------------------------
! Interface for NWChem based QM MD
!
! Currently supports:
! pure QM
!
! Mark J. Williamson (Unilever Centre, Cambridge)
!
! Date: February 2012
! Still TODO:
!   Templates
!   Dipole support
!   PIMD
!   Deprecate use of null parameter holders null{i,r}
!   Check energy conservation and possibly let users adjust
!    integral neglect thresholds / XC quadrature grid params etc
!
! Current issues:
!   NWchem crashes when a large number of BQ points (~1000) are
!   defined within the geometry module.
!
! ----------------------------------------------------------------

    implicit none

    private
    public :: get_nw_forces

    type nw_nml_type
        character(len=20) :: method
        character(len=20) :: basis
        integer :: scf_conv
        integer :: charge
        integer :: spinmult
        integer :: ntpr
        integer :: num_threads
        integer :: verbosity
    end type nw_nml_type

contains

    ! --------------------------------------------
    ! Get QM energy and forces from NWChem
    ! --------------------------------------------
    subroutine get_nw_forces(do_grad, nstep, ntpr_default, id, nqmatoms, qmcoords, &
        nclatoms, clcoords, escf, dxyzqm, dxyzcl, charge)

        use constants, only : CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO
        use file_io_dat

        implicit none

        logical, intent(in) :: do_grad              ! Return gradient/not
        integer, intent(in) :: nstep                ! MD step number
        integer, intent(in) :: ntpr_default         ! frequency of printing
        character(len=3), intent(in) :: id          ! ID number for PIMD or REMD
        integer, intent(in) :: nqmatoms             ! Number of QM atoms
        _REAL_, intent(in) :: qmcoords(3, nqmatoms) ! QM atom coordinates
        integer, intent(in) :: nclatoms             ! Number of MM atoms
        _REAL_, intent(in) :: clcoords(4, nclatoms) ! MM atom coordinates and charges in au
        _REAL_, intent(out) :: escf                 ! SCF energy
        _REAL_, intent(out) :: dxyzqm(3, nqmatoms)   ! SCF QM force
        _REAL_, intent(out) :: dxyzcl(3, nclatoms)   ! SCF MM force
        integer, intent(out) :: charge

        type(nw_nml_type), save      :: nw_nml
        logical, save                :: first_call = .true.
        integer                      :: i
        integer                      :: printed = -1 ! Used to tell if we have printed this step yet
        ! since the same step may be called multiple times
        character(len=150)          :: call_buffer
        character(len=6), save       :: program
        character(len=*), parameter  :: basename = 'nwchem'
        character(len=*), parameter  :: inpext = '.nw'
        character(len=*), parameter  :: logext = '.log'
        character(len=*), parameter  :: movecext = '.movecs' ! MO vector file
        character(len=*), parameter  :: dbext = '.db'     ! Database file; essentially current state of calculation
        character(len=14)            :: inpfile, movecfile, logfile, dbfile
        ! Need to prepend subdirectory if doing REMD, PIMD
        character(len=25)            :: subdir

        ! for system calls
        integer :: system
        integer :: stat

        ! assemble input - / output data filenames
        inpfile = basename//trim(id)//inpext
        logfile = basename//trim(id)//logext
        movecfile = basename//trim(id)//movecext
        dbfile = basename//trim(id)//dbext

        ! Setup on first call
        if (first_call) then
            first_call = .false.
            write (6, '(/,a,/)') '  >>> Running calculations with NWChem <<<'
            call get_namelist(ntpr_default, nw_nml)
            call print_namelist(nw_nml)
            charge = nw_nml%charge
            ! Check for version of NWChem to use; store as 'program'
            call check_installation(program, id)
            write (6, '(80a)') ('-', i=1, 80)
            write (6, '(a)') '   4.  RESULTS'
            write (6, '(80a)') ('-', i=1, 80)
            ! Remove old inpfile, logfile, database file and movecs file at the
            ! beginning of a run so only the latest run is stored.
            stat = system('rm -f '//inpfile//' '//dbfile//' '//logfile//' '//movecfile)

            if (stat /= 0) then
                call sander_bomb('get_nw_forces (qm2_extern_nw_module)', &
                    'Error with system call (removing files)', &
                    'Will quit now.')
            end if
        end if

        call write_inpfile(trim(inpfile), &
            nqmatoms, qmcoords, nclatoms, clcoords, nw_nml, do_grad)

        if (nw_nml%verbosity > 0) then
            write (6, '(a)') 'Runfile has written successfully; Calling NWChem...'
        end if

        ! Run NWChem
        ! Separate runs into different directories if we are doing PIMD
        subdir = ''
        call_buffer = ''
        if (trim(id) /= '') then
            subdir = './'//trim(id)//'/'
            call_buffer = ' mkdir -p '//trim(subdir)//'; cd '//trim(subdir)//'; mv ../'//inpfile//' .;'
        end if
        stat = system(program//' '//inpfile//'>'//logfile)
        if (stat /= 0) then
            call sander_bomb('get_nw_forces (qm2_extern_nw_module)', &
                'Error with system call (executing NWChem)', &
                'Will quit now.')
        end if

        if (nw_nml%verbosity > 0) then
            write (6, '(a)') 'NWChem execution success; Processing NWChem results...'
        end if

        ! Call read_results - retrieve data from NWChem .log file
        ! Will output data to escf and dxyqm for pure QM or MM runs
        ! For QM/MM runs will also return dxyzcl containing electric field strength

        ! Search in subdir for logfile if doing PIMD
        ! Otherwise, search current directory
        call read_results(nw_nml, trim(subdir)//trim(logfile), &
            nqmatoms, escf, dxyzqm, nclatoms, dxyzcl, do_grad)

        ! Save copy of last input and log files
        stat = system('mv '//trim(subdir)//inpfile//' '//trim(subdir)//'old.'//inpfile)
        stat = stat + system('mv '//trim(subdir)//logfile//' '//trim(subdir)//'old.'//logfile)
        if (stat /= 0) then
            call sander_bomb('get_nw_forces (qm2_extern_nw_module)', &
                'Error with system call (moving / removing files)', &
                'Will quit now.')
        end if

        ! F = E*q to get gradients
        ! Note dxyzcl is currently holding the electric field strength
        do i = 1, nclatoms
            dxyzcl(:, i) = -dxyzcl(:, i)*clcoords(4, i)
        end do

        ! Convert gradient from au to kcal/(mol*A)
        if (do_grad) then
            dxyzqm(:, :) = dxyzqm(:, :)*CODATA08_AU_TO_KCAL*CODATA08_A_TO_BOHRS
            if (nclatoms > 0) then
                dxyzcl(:, :) = dxyzcl(:, :)*CODATA08_AU_TO_KCAL*CODATA08_A_TO_BOHRS
            end if
        else
            dxyzqm = ZERO
            if (nclatoms > 0) dxyzcl = ZERO
        end if

        if (nw_nml%verbosity > 0) then
            write (6, '(a)') 'get_nw_forces - final gradient(s):'
            write (6, '(a)') 'QM region:'
            do i = 1, nqmatoms
                write (6, *) dxyzqm(1:3, i)
            end do
            write (6, *)
            if (nclatoms > 0) then
                write (6, '(a)') 'MM region:'
                do i = 1, nclatoms
                    write (6, '(3(x,f16.10))') dxyzcl(1:3, i)
                end do
            end if
        end if

        escf = escf*CODATA08_AU_TO_KCAL

    end subroutine get_nw_forces

    ! ---------------------------------------------
    ! Read NWChem namelist values from file mdin,
    ! use default values if none are present.
    ! ---------------------------------------------

    subroutine get_namelist(ntpr_default, nw_nml)

        implicit none
        integer, intent(in) :: ntpr_default
        type(nw_nml_type), intent(out) :: nw_nml

        character(len=20) :: method, basis
        integer :: charge, spinmult, verbosity
        integer :: scf_conv, ntpr, num_threads
        namelist /nw/ method, basis, scf_conv, ntpr, charge, spinmult, &
            num_threads, verbosity

        integer :: ierr

        ! Set default values for nw namelist values
        method = 'BLYP'
        basis = '6-31G*'
        scf_conv = 8
        charge = 0
        spinmult = 1
        ntpr = ntpr_default
        num_threads = 1
        verbosity = 0

        ! Read namelist
        rewind 5
        read (5, nml=nw, iostat=ierr)

        if (ierr > 0) then
            call sander_bomb('get_namelist (qm2_extern_nw_module)', &
                '&nw namelist read error', &
                'Please check your input.')
        else if (ierr < 0) then
            write (6, '(a/a)') '&nw namelist read encountered end of file', &
                'Please check your input if the calculation encounters a problem'
        end if

        ! Assign namelist values to nw_nml data type
        nw_nml%method = method
        nw_nml%basis = basis
        nw_nml%scf_conv = scf_conv
        nw_nml%charge = charge
        nw_nml%spinmult = spinmult
        nw_nml%ntpr = ntpr
        nw_nml%num_threads = num_threads
        nw_nml%verbosity = verbosity

    end subroutine get_namelist

    ! --------------------------------
    ! Print NWChem namelist settings
    ! --------------------------------
    subroutine print_namelist(nw_nml)

        implicit none
        type(nw_nml_type), intent(in) :: nw_nml

        write (6, '(a)') '| &nw'
        write (6, '(2a)') '|   method       = ', nw_nml%method
        write (6, '(2a)') '|   basis        = ', nw_nml%basis
        write (6, '(a,i2)') '|   scf_conv     = ', nw_nml%scf_conv
        write (6, '(a,i2)') '|   charge       = ', nw_nml%charge
        write (6, '(a,i2)') '|   spinmult     = ', nw_nml%spinmult
        write (6, '(a,i0)') '|   ntpr         = ', nw_nml%ntpr
        write (6, '(a,i2)') '|   num_threads  = ', nw_nml%num_threads
        write (6, '(a,i2)') '|   verbosity    = ', nw_nml%verbosity
        write (6, '(a)') '| /'

    end subroutine print_namelist

    ! --------------------------------------------
    ! Check whether NWChem is properly installed
    ! Also, return version of NWchem to use
    ! --------------------------------------------
    subroutine check_installation(program, id)

        implicit none

        character(len=6), intent(out) :: program
        character(len=3), intent(in) :: id

        character(len=80) :: read_buffer
        character(len=80) :: call_buffer
        character(len=80) :: filename
        integer :: iunit = 77
        integer :: stat
        integer :: system

        filename = 'extern_location'//trim(id)

        ! Search for NWChem executable
        call_buffer = 'which nwchem > '//trim(filename)
        stat = system(trim(call_buffer))
        if (stat == 0) then
            program = 'nwchem'
        else
            call sander_bomb('check_installation (qm2_extern_nw_module)', &
                'Executable NWChem not found', &
                'Please check your NWChem installation')
        end if

        write (6, '(3a,/)') '| Program ', program, ' found!'

        ! Get complete executable path
        open (unit=iunit, file=trim(filename), form='formatted', iostat=stat)
        if (stat /= 0) then
            call sander_bomb('check_installation (qm2_extern_nw_module)', &
                'Internal error opening file with path location of executable.', &
                'Quitting now.')
        end if
        read (iunit, '(a)', iostat=stat) read_buffer
        if (stat /= 0) then
            call sander_bomb('check_installation (qm2_extern_nw_module)', &
                'Internal error reading from file with path location of executable.', &
                'Quitting now.')
        end if
        close (unit=iunit, status='delete', iostat=stat)
        if (stat /= 0) then
            call sander_bomb('check_installation (qm2_extern_nw_module)', &
                'Internal error closing and deleting file with path location of executable.', &
                'Quitting now.')
        end if

        write (6, '(2a,/)') '| Executable location: ', trim(read_buffer)

    end subroutine check_installation

    ! -----------------------------
    ! Write input file for NWChem
    ! -----------------------------

    subroutine write_inpfile(inpfile, nqmatoms, qmcoords, &
        nclatoms, clcoords, nw_nml, do_grad)

        use qmmm_module, only : qmmm_struct
        use ElementOrbitalIndex, only : elementSymbol

        implicit none

        character(len=*), intent(in)   :: inpfile
        integer, intent(in)            :: nqmatoms
        _REAL_, intent(in)            :: qmcoords(:, :)
        integer, intent(in)            :: nclatoms
        _REAL_, intent(in)            :: clcoords(:, :)
        type(nw_nml_type), intent(in)  :: nw_nml
        logical, intent(in)            :: do_grad

        integer, parameter :: iunit = 351, tplunit = 352
        integer            :: i, ierr
        logical, save      :: first_call = .true.

        if (nw_nml%verbosity > 0) then
            write (6, '(a)') 'Writing NWChem inpfile '//inpfile
        end if

        open (iunit, file=inpfile, iostat=ierr)
        if (ierr /= 0) then
            call sander_bomb('write_inpfile (qm2_extern_nw_module)', &
                'Error opening NWChem inpfile '//inpfile//' for writing', &
                'Will quit now.')
        end if

        ! Write title
        write (iunit, '(a)') 'title "NWChem run using SANDER external interface."'

        ! Charge
        write (iunit, '(a,i0)'), 'charge ', nw_nml%charge
        write (iunit, '(a)')

        ! Write QM atoms and coordinates
        ! More info here: http://www.nwchem-sw.org/index.php/Release61:Geometry
        ! TODO, check if bqbq is needed or not
        write (iunit, '(a)'), 'geometry units an noautoz nocenter noautosym noprint'
        do i = 1, nqmatoms
            write (iunit, '(a2,1x,3f25.16)') elementSymbol(qmmm_struct%iqm_atomic_numbers(i)), qmcoords(1:3, i)
        end do

        ! If QM/MM, write the external point charges as a set of dummy points here
        ! so that the corresponding Electric field is calculated for them in addition
        ! to the QM atoms. This is fragile since NWChem can only deal with a finite
        ! number of dummy bq points within the geometry module.

        if (nclatoms > 0) then
            do i = 1, nclatoms
                write (iunit, '(a2,1x,3f25.16)') "bq", clcoords(1:3, i)
            end do
        end if

        write (iunit, '(a)'), 'end'
        write (iunit, '(a)')

        ! When electrostatic embedding (QM/MM) is in use write MM coordinates with point charges.
        ! This is written to the bq module.
        if (nclatoms > 0) then
            write (iunit, '(a)'), 'bq'
            do i = 1, nclatoms
                ! Note the reordering here:
                !  clcoords(1:3,i) are coordinates and clcoords(4,i) is charge
                !  However, NWChem expects the format charge,x,y,z in the bq module
                write (iunit, '(f21.16,3f21.16)') clcoords(4, i), clcoords(1:3, i)
            end do
            write (iunit, '(a)'), 'end'
            write (iunit, '(a)')
        end if

        ! Basis
        write (iunit, '(a)'), 'basis'
        write (iunit, '(a,a)'), '  * library ', trim(nw_nml%basis)
        write (iunit, '(a)'), 'end'
        write (iunit, '(a)')

        ! Method
        ! TODO; this is far too fagile

        ! HF Method
        ! More info here: http://www.nwchem-sw.org/index.php/Release61:Hartree-Fock_Theory_for_Molecules
        if (nw_nml%method == "hf") then
            write (iunit, '(a)'), 'scf'
            ! Convergence
            write (iunit, '(a,i1)'), '  thresh 1e-', nw_nml%scf_conv
            ! Calculate all integrals "on-the-fly"
            write (iunit, '(a)'), '  direct '
            ! Do not print the MO vector coefficients; just too much data.
            write (iunit, '(a)'), '  noprint "final vectors analysis"'
            write (iunit, '(a)'), 'end'
        end if
        write (iunit, '(a)')

        ! DFT Method
        ! More info here: http://www.nwchem-sw.org/index.php/Density_Functional_Theory_for_Molecules
        if (nw_nml%method == "BLYP") then
            write (iunit, '(a)'), 'dft'
            ! Convergence
            write (iunit, '(a,i1)'), '  convergence energy 1e-', nw_nml%scf_conv
            ! Grid
            write (iunit, '(a)'), '  grid fine'
            write (iunit, '(a,a,a)'), '  xc', ' becke88', ' lyp'
            write (iunit, '(a,i1)'), '  mult ', nw_nml%spinmult
            write (iunit, '(a)'), '  noio '
            ! Calculate all integrals "on-the-fly"
            write (iunit, '(a)'), '  direct '
            ! Do not print the MO vector coefficients; just too much data.
            write (iunit, '(a)'), '  noprint "final vectors analysis"'
            write (iunit, '(a)'), 'end'
        end if
        write (iunit, '(a)')

        ! QM/MM
        !
        ! This outputs the electric field on all atoms + point charges in the
        ! system, which is harvested into dxyzcl(3,nclatoms).
        ! This needs the property task (below) to be executed to actually
        ! generate output.

        if (nclatoms > 0) then
            write (iunit, '(a)') 'property'
            write (iunit, '(a)') '  efield'
            write (iunit, '(a)') 'end'
            write (iunit, '(a)')
        end if

        ! One will need at least one of these task directives

        if (nw_nml%method == "hf") then
            if (do_grad) then
                if (nclatoms > 0) then
                    ! QM/MM
                    ! Note, 'property' is needed in the task directive to trigger the
                    ! evaluation of the property directive which has the 'efield' keyword:
                    ! http://www.nwchem-sw.org/index.php/Release61:Properties
                    !
                    write (iunit, '(a)') 'task scf gradient'
                    write (iunit, '(a)') 'task scf property'
                else
                    ! QM only; 'efield' output not needed
                    write (iunit, '(a)') 'task scf gradient'
                end if
            else
                ! Don't calculate gradient anymore
                if (nclatoms > 0) then
                    ! QM/MM
                    write (iunit, '(a)') 'task scf energy'
                    write (iunit, '(a)') 'task scf property'
                else
                    ! QM
                    write (iunit, '(a)') 'task scf energy'
                end if
            end if
            write (iunit, '(a)')
        end if

        ! TODO, this is just not good enough.
        ! We  needs a generic list of DFT methods here
        ! i.e. if (nw_nml%method == "BLYP" | "B3LYP" | ) then
        ! or perhaps even change the "nw_nml_type" to include
        ! a functional entry?

        if (nw_nml%method == "BLYP") then
            if (do_grad) then
                if (nclatoms > 0) then
                    ! QM/MM
                    ! Note, 'property' is needed in the task directive to trigger the
                    ! evaluation of the property directive which has the 'efield' keyword.
                    ! This outputs the electric field on all atoms + point charges in the
                    ! system, which is harvested into dxyzcl(3,nclatoms)
                    write (iunit, '(a)') 'task dft gradient'
                    write (iunit, '(a)') 'task dft property'
                else
                    ! QM
                    write (iunit, '(a)') 'task dft gradient'
                end if
            else
                ! Don't calculate gradient anymore
                if (nclatoms > 0) then
                    write (iunit, '(a)') 'task dft energy property'
                    write (iunit, '(a)') 'task dft property'
                else
                    write (iunit, '(a)') 'task dft energy'
                end if
            end if
            write (iunit, '(a)')
        end if

        close (iunit, iostat=ierr)

        if (ierr /= 0) then
            call sander_bomb('write_inpfile (qm2_extern_nw_module)', &
                'Error closing NWChem runfile after writing', &
                'Will quit now.')
        end if
        first_call = .false.

    end subroutine write_inpfile

    ! Parse the output of the nwchem.log file It is aiming to extract the following:
    !
    !   * QM Energy
    !   * Gradient on QM atoms
    !   * Field on MM atoms

    subroutine read_results(nw_nml, datfile, nqmatoms, escf, dxyzqm, &
        nclatoms, dxyzcl, do_grad)

        implicit none

        type(nw_nml_type), intent(in) :: nw_nml ! Needed to determine type of calculation for parsing
        character(len=*), intent(in)  :: datfile
        integer, intent(in)           :: nqmatoms, nclatoms
        _REAL_, intent(out)           :: escf, dxyzqm(3, nqmatoms), &
            dxyzcl(3, nclatoms) ! dxyzcl will return containing the electric field at x,y,z
        logical, intent(in)           :: do_grad

        integer :: ios, i
        integer :: nulli
        _REAL_  :: nullr
        integer, parameter :: iunit = 351
        character(len=120) :: read_buffer

        open (iunit, file=datfile, status='old', iostat=ios)
        if (ios /= 0) then
            call sander_bomb('read_results (qm2_extern_nw_module)', &
                'Error opening NWChem log file '//datfile//' (expected in same dir as input file).', &
                'Will quit now')
        end if

        do
            read (iunit, '(a)', iostat=ios) read_buffer
            ! End of file; nothing left to read
            if (ios < 0) then
                exit
            end if

            if (nw_nml%method == "hf") then
                !         Total SCF energy =   -243.826232163838
                !      One-electron energy =   -684.067100443940
                !      Two-electron energy =    263.613349226485
                ! Nuclear repulsion energy =    176.627519053618

                if (read_buffer(1:25) == '         Total SCF energy') then
                    ! Read value after "=" sign
                    read (read_buffer(index(read_buffer, '=') + 1:), *) escf
                end if

            end if

            if (nw_nml%method == "BLYP") then
                !         Total DFT energy =      -74.827214054279
                !      One electron energy =     -115.784049745746
                !           Coulomb energy =       43.899915115726
                !    Exchange-Corr. energy =       -8.837093294400
                ! Nuclear repulsion energy =        5.894013870142

                if (read_buffer(1:25) == '         Total DFT energy') then
                    ! Read value after "=" sign
                    read (read_buffer(index(read_buffer, '=') + 1:), *) escf
                end if

            end if

            if (do_grad) then
                ! Gradients
                !                         DFT ENERGY GRADIENTS
                !
                !    atom               coordinates                        gradient
                !                 x          y          z           x          y          z
                !   1 O       0.000000   0.000000   0.400871    0.000000   0.000000   0.151975
                !   2 H       2.004357   0.000000  -1.603486    0.073745   0.000000  -0.075987
                !   3 H      -2.004357   0.000000  -1.603486   -0.073745   0.000000  -0.075987

                if (read_buffer(1:80) == '    atom               coordinates                        gradient') then
                    ! Skip over this line
                    read (iunit, '(a)') read_buffer
                    ! Skip over next line (x,y,z heading)
                    read (iunit, '(a)') read_buffer

                    if (nw_nml%verbosity > 0) then
                        write (6, '(a)') 'read_results() - read in gradients:'
                        write (6, '(a)') 'QM region:'
                    end if

                    do i = 1, nqmatoms

                        ! MJW TODO; this is too dirty
                        read (read_buffer, '(1X,I3,1X,A4,2(1X,3(1X,F10.6)))', iostat=ios) &
                            nulli, nulli, nullr, nullr, nullr, dxyzqm(1:3, i)
                        if (nw_nml%verbosity > 0) then
                            write (6, *) dxyzqm(1:3, i)
                        end if

                        ! Next line
                        read (iunit, '(a)') read_buffer
                    end do
                end if

                ! QM/MM
                if (nclatoms > 0) then
                    !  Read Efield of MM charges
                    !   Atom       X         Y         Z                        Electric field (a.u.)
                    !                                              X              Y              Z           Field
                    !  ------------------------------------------------------------------------------------------------
                    !    1 N    7.10169   7.39654  -0.33737       0.003260      -0.023523       0.036570       0.043604
                    !    2 H    5.61573   8.49194  -0.45940      -0.072568       0.062168      -0.005779       0.095730
                    !    3 C    9.60125   8.63375  -0.50305       0.006860      -0.001083      -0.002420       0.007354
                    !   11 H    6.84332   5.35446  -0.26037       0.010656       0.063760      -0.000629       0.064647
                    !   12 H   10.76253  12.61212   0.19451       0.074180       0.055897       0.023027       0.095694
                    !   ........
                    !   13 bq   3.92568   1.96181   0.26579      -0.006272      -0.004868       0.001409       0.008064
                    !   14 bq   3.99289   4.04591   0.22189      -0.018218      -0.002074       0.004034       0.018774

                    if (read_buffer(1:80) == '   Atom       X         Y         Z                        Electric field (a.u.)') then
                        ! Skip over this line
                        read (iunit, '(a)') read_buffer
                        ! Skip over next line (x,y,z,Field heading)
                        read (iunit, '(a)') read_buffer
                        ! Skip over next line ------------------
                        read (iunit, '(a)') read_buffer

                        ! Skip over all QM atoms
                        do i = 1, nqmatoms
                            read (iunit, '(a)') read_buffer
                        end do

                        if (nw_nml%verbosity > 0) then
                            write (6, '(a)') 'MM region:'
                        end if

                        do i = 1, nclatoms
                            read (read_buffer, '(i5,1x,a2,3f10.5,4f15.6)', iostat=ios) &
                                nulli, nulli, nullr, nullr, nullr, dxyzcl(1:3, i), nullr

                            if (nw_nml%verbosity > 0) then
                                write (6, *) dxyzcl(1:3, i)
                            end if

                            ! Next line
                            read (iunit, '(a)') read_buffer
                        end do

                    end if

                end if

            end if
        end do

        close (iunit)

    end subroutine read_results

end module qm2_extern_nw_module
