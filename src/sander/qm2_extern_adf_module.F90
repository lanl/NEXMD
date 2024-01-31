#include "dprec.fh"
module qm2_extern_adf_module
! ----------------------------------------------------------------
! Interface for ADF based QM MD
!
! Currently supports:
! pure QM
!
! Initial implementation by
! Matthew Clark and Prithvi Undavalli
! (SDSC LSSI summer highschool students)
! under supervision of
! Andreas Goetz and Ross Walker (SDSC)
!
! Date: August 2010
!
! Extensions by Andreas Goetz (agoetz@sdsc.edu)
!
! Date: November 2010
!
! ----------------------------------------------------------------

    implicit none

    private
    public :: get_adf_forces

    type adf_nml_type
        character(len=20) :: xc
        character(len=20) :: basis
        character(len=20) :: core
        character(len=20) :: fit_type
        _REAL_ :: integration
        _REAL_ :: scf_conv
        integer :: charge
        integer :: spin
        integer :: scf_iter
        integer :: ntpr
        integer :: num_threads
        integer :: linear_scaling
        integer :: verbosity
        logical :: use_dftb
        logical :: oldgradients
        logical :: dipole
        logical :: exactdensity
        logical :: use_template
    end type adf_nml_type

contains

    ! -------------------------------------------------------
    ! Get QM energy and forces from ADF (adf or dftb program)
    ! -------------------------------------------------------
    subroutine get_adf_forces(do_grad, nstep, ntpr_default, id, natoms, coords, escf, dxyzqm, charge)

        use constants, only : CODATA08_AU_TO_KCAL, CODATA08_A_TO_BOHRS, ZERO
        use file_io_dat

        implicit none

        logical, intent(in) :: do_grad            ! Return gradient/not
        integer, intent(in) :: nstep              ! MD step number
        integer, intent(in) :: ntpr_default       ! frequency of printing
        character(len=3), intent(in) :: id        ! ID number for PIMD or REMD
        integer, intent(in) :: natoms             ! Total number of atoms
        _REAL_, intent(in) :: coords(3, natoms)   ! QM atom coordinates
        _REAL_, intent(out) :: escf               ! SCF energy
        _REAL_, intent(out) :: dxyzqm(3, natoms)   ! SCF QM force
        integer, intent(out) :: charge

        _REAL_              :: dipxyz(3)          ! Dipole Moment

        type(adf_nml_type), save :: adf_nml
        logical, save :: first_call = .true.
        integer :: i
        integer :: printed = -1 ! Used to tell if we have printed this step yet
        ! since the same step may be called multiple times
        character(len=150)          :: call_buffer
        character(len=*), parameter :: basename = 'adf_job'
        character(len=*), parameter :: inpext = '.inp'
        character(len=*), parameter :: outext = '.out'
        character(len=*), parameter :: dipext = '.dip'
        character(len=*), parameter :: tplext = '.tpl'
        character(len=14)            :: inpfile, outfile, dipfile, keyfile, tplfile
        ! Need to prepend subdirectory if doing REMD, PIMD
        character(len=25)            :: subdir

        ! assemble input - / output data filenames
        inpfile = basename//trim(id)//inpext
        outfile = basename//trim(id)//outext
        dipfile = basename//trim(id)//dipext
        tplfile = basename//tplext

        ! Setup on first call
        if (first_call) then
            first_call = .false.
            write (6, '(/,a,/)') '  >>> Running QM calculation with ADF <<<'
            call get_namelist(ntpr_default, adf_nml)
            call print_namelist(adf_nml)
            if (adf_nml%use_dftb) then
                call check_installation('dftb')
            else
                call check_installation('adf')
            end if
            charge = adf_nml%charge
            write (6, '(80a)') ('-', i=1, 80)
            write (6, '(a)') '   4.  RESULTS'
            write (6, '(80a)') ('-', i=1, 80)
        end if
        if (adf_nml%use_dftb) then
            keyfile = 'DFTB.kf'
        else
            keyfile = 'TAPE21'
        end if
        ! Remove the logfile at the beginning of a run so only the latest run is
        ! stored. ADF will also throw a false 'error detected' if it finds a TAPE13
        ! file in the working directory; this is to prevent confusion.
        call system('rm -f logfile TAPE13')

        call system('rm -f '//inpfile)
        call write_inpfile(trim(inpfile), trim(tplfile), coords, natoms, &
            adf_nml, do_grad)

        if (adf_nml%verbosity > 0) then
            write (6, *) 'Runfile has written successfully; Calling ADF...'
        end if

        ! Run ADF
        ! Separate runs into different directories if we are doing PIMD
        subdir = ''
        call_buffer = ''
        if (trim(id) /= '') then
            subdir = './'//trim(id)//'/'
            call_buffer = ' mkdir -p '//trim(subdir)//'; cd '//trim(subdir)//'; mv ../'//inpfile//' .;'
        end if

        ! Inputting 0 will cause ADF to use default number of threads
        if (adf_nml%num_threads == 0) then
            ! Enable the use of DFTB if use_dftb flag is set
            if (adf_nml%use_dftb) then
                write (call_buffer, '(a)') trim(call_buffer)//'$ADFBIN/dftb'
            else
                write (call_buffer, '(a)') trim(call_buffer)//'$ADFBIN/adf'
            end if
            ! Otherwise, we write what the user specified
        else
            if (adf_nml%use_dftb) then
                write (call_buffer, '(a,i0)') trim(call_buffer)//'$ADFBIN/dftb -n ', adf_nml%num_threads
            else
                write (call_buffer, '(a,i0)') trim(call_buffer)//'$ADFBIN/adf -n ', adf_nml%num_threads
            end if
        end if

        write (call_buffer, '(a)') trim(call_buffer)//' < ./'//trim(inpfile)//' > '//outfile
        call system(trim(call_buffer))

        if (adf_nml%verbosity > 0) then
            write (6, *) 'ADF execution success; Processing ADF results...'
        end if

        ! Call read_adf_results - a function in C program (qm2_read_adf_results.c)
        ! to read TAPE21 file; will output the data to escf and dxyzqm
        call read_adf_results(natoms, escf, dxyzqm, adf_nml%use_dftb, dipxyz, &
            trim(subdir)//trim(keyfile), do_grad)

        ! Call write_results - dipfile, dipxyz, and magnitude of dipxyz;
        ! will write output to .dip and .chg files
        if (adf_nml%ntpr > 0 .and. mod(nstep, adf_nml%ntpr) == 0) then
            if (printed /= nstep .and. adf_nml%dipole) then
                call write_results(dipfile, dipxyz, (dipxyz(1)**2 + dipxyz(2)**2 + dipxyz(3)**2)**.5)
                printed = nstep
            end if
        end if

        if (adf_nml%use_dftb) then
            call system('mv '//trim(subdir)//keyfile//' '//trim(subdir)//'kf.DFTB')
        else
            ! Save old TAPE21 file to restart from in next calculation
            call system('mv '//trim(subdir)//keyfile//' '//trim(subdir)//'adf.t21')
        end if

        ! Convert gradient from au to kcal/(mol*A)
        if (do_grad) then
            dxyzqm(:, :) = dxyzqm(:, :)*CODATA08_AU_TO_KCAL*CODATA08_A_TO_BOHRS
        else
            dxyzqm = ZERO
        end if

        if (adf_nml%verbosity > 0) then
            write (6, '(a)') 'get_adf_forces - final gradient:'
            do i = 1, natoms
                write (6, *) dxyzqm(1:3, i)
            end do
        end if

        escf = escf*CODATA08_AU_TO_KCAL

    end subroutine get_adf_forces

    ! ----------------------------------------
    ! Read ADF namelist values from file mdin,
    ! use default values if none are present.
    ! ----------------------------------------
    subroutine get_namelist(ntpr_default, adf_nml)

        implicit none
        integer, intent(in) :: ntpr_default
        type(adf_nml_type), intent(out) :: adf_nml

        character(len=20) :: xc, basis, core, fit_type
        integer :: charge, spin, scf_iter, ntpr, num_threads, linear_scaling, verbosity, use_dftb, &
            oldgradients, dipole, exactdensity, use_template
        _REAL_ :: integration, scf_conv, scf_conv_dftb, scf_conv_adf
        namelist /adf/ xc, basis, core, fit_type, integration, scf_conv, charge, spin, scf_iter, ntpr, &
            num_threads, linear_scaling, verbosity, use_dftb, oldgradients, dipole, exactdensity, use_template

        integer :: ifind, ierr

        ! Set default values for adf namelist values
        xc = 'GGA BLYP'
        basis = 'DZP'
        core = 'None'
        fit_type = 'ZORA/QZ4P'
        integration = 5.0D0
        scf_conv = 1.0D30    ! dummy value
        scf_conv_dftb = 1.0d-12   ! dftb default value
        scf_conv_adf = 1.0d-06   ! adf default value
        charge = 0
        spin = 0
        ntpr = ntpr_default
        scf_iter = 50
        num_threads = 0         ! 0 = ADF default (all available)
        linear_scaling = -1
        verbosity = 0
        use_dftb = 0
        oldgradients = 0
        dipole = 0
        exactdensity = 0
        use_template = 0

        ! Read namelist
        rewind 5
        read (5, nml=adf, iostat=ierr)

        if (ierr > 0) then
            call sander_bomb('get_namelist (qm2_extern_adf_module)', &
                '&adf namelist read error', &
                'Please check your input.')
        else if (ierr < 0) then
            write (6, '(a/a)') '&adf namelist read encountered end of file', &
                'Please check your input if the calculation encounters a problem'
        end if

        ! Assign namelist values to adf_nml data type
        adf_nml%xc = xc
        adf_nml%basis = basis
        adf_nml%core = core
        adf_nml%fit_type = fit_type
        adf_nml%integration = integration
        if (use_dftb > 0) then
            adf_nml%use_dftb = .true.
        else
            adf_nml%use_dftb = .false.
        end if
        if (scf_conv > 1.0D29) then
            adf_nml%scf_conv = scf_conv_dftb ! use default value
        else
            adf_nml%scf_conv = scf_conv
        end if
        if (scf_conv > 1.0D29) then
            adf_nml%scf_conv = scf_conv_adf  ! use default value
        else
            adf_nml%scf_conv = scf_conv
        end if
        adf_nml%charge = charge
        adf_nml%spin = spin
        adf_nml%scf_iter = scf_iter
        adf_nml%ntpr = ntpr
        adf_nml%num_threads = num_threads
        adf_nml%linear_scaling = linear_scaling
        adf_nml%verbosity = verbosity

        if (oldgradients == 1) then
            adf_nml%oldgradients = .true.
        else if (oldgradients == 0) then
            adf_nml%oldgradients = .false.
        else
            call sander_bomb('get_namelist (qm2_extern_adf_module)', &
                '&adf oldgradients value not allowed', &
                'Please check your input. oldgradients can only be 0 or 1.')
        end if

        if (dipole == 1) then
            adf_nml%dipole = .true.
        else if (dipole == 0) then
            adf_nml%dipole = .false.
        else
            call sander_bomb('get_namelist (qm2_extern_adf_module)', &
                '&adf dipole value not allowed', &
                'Please check your input. dipole can only be 0 or 1.')
        end if

        if (exactdensity == 1) then
            adf_nml%exactdensity = .true.
        else if (exactdensity == 0) then
            adf_nml%exactdensity = .false.
        else
            call sander_bomb('get_namelist (qm2_extern_adf_module)', &
                '&adf exactdensity value not allowed', &
                'Please check your input. exactdensity can only be 0 or 1.')
        end if

        if (use_template == 0) then
            adf_nml%use_template = .false.
        else if (use_template == 1) then
            adf_nml%use_template = .true.
        else
            call sander_bomb('get_namelist (qm2_extern_adf_module)', &
                '&adf use_template value not allowed', &
                'Please check your input. use_template can only be 0 or 1.')
        end if

    end subroutine get_namelist

    ! ---------------------------
    ! Print ADF namelist settings
    ! ---------------------------
    subroutine print_namelist(adf_nml)

        implicit none
        type(adf_nml_type), intent(in) :: adf_nml

        write (6, '(a)') '| &adf'
        if (.not. adf_nml%use_dftb) then
            write (6, '(2a)') '|   xc             = ', adf_nml%xc
            write (6, '(2a)') '|   basis          = ', adf_nml%basis
            write (6, '(2a)') '|   core           = ', adf_nml%core
            write (6, '(2a)') '|   fit_type       = ', adf_nml%fit_type
            write (6, '(a,E22.16)') '|   integration    = ', adf_nml%integration
        end if
        write (6, '(a,es10.2)') '|   scf_conv       = ', adf_nml%scf_conv
        write (6, '(a,i2)') '|   charge         = ', adf_nml%charge
        write (6, '(a,i2)') '|   spin           = ', adf_nml%spin
        write (6, '(a,i5)') '|   scf_iter       = ', adf_nml%scf_iter
        write (6, '(a,i0)') '|   ntpr           = ', adf_nml%ntpr
        write (6, '(a,i3)') '|   num_threads    = ', adf_nml%num_threads
        write (6, '(a,i3)') '|   linear_scaling = ', adf_nml%linear_scaling
        write (6, '(a,l)') '|   use_dftb       = ', adf_nml%use_dftb
        if (.not. adf_nml%use_dftb) then
            write (6, '(a,l)') '|   oldgradients   = ', adf_nml%oldgradients
            write (6, '(a,l)') '|   dipole         = ', adf_nml%dipole
            write (6, '(a,l)') '|   exactdensity   = ', adf_nml%exactdensity
        end if
        write (6, '(a,l)') '|   use_template   = ', adf_nml%use_template
        write (6, '(a)') '| /'

    end subroutine print_namelist

    ! ---------------------------------------
    ! Check whether ADF is properly installed
    ! ---------------------------------------
    subroutine check_installation(program)

        implicit none

        character(len=*), intent(in) :: program

        logical :: exist
        character(len=80) :: env_path
        integer :: len, stat

        ! Search for $ADFBIN
        call get_environment_variable('ADFBIN', env_path, len, stat)
        if (stat /= 0) then
            call sander_bomb('check_installation (qm2_extern_adf_module)', &
                'Environment variable "$ADFBIN" not set', &
                'Please check your ADF installation')
        end if

        env_path = trim(env_path)//'/'//program
        write (6, '(2a)') '| Searching for ', env_path
        inquire (FILE=env_path, EXIST=exist)

        if (exist .eqv. .true.) then
            write (6, '(3a)') '| Program ', program, ' found!'
        else
            call sander_bomb('check_installation (qm2_extern_adf_module)', &
                'Program '//program//' not found', &
                'External program "'//program//'" required to use the extern_adf module.')
        end if

    end subroutine check_installation

    ! ------------------------
    ! Write input file for ADF
    ! ------------------------
    subroutine write_inpfile(inpfile, tplfile, coords, natoms, adf_nml, do_grad)

        use qmmm_module, only : qmmm_struct
        use ElementOrbitalIndex, only : elementSymbol

        implicit none

        character(len=*), intent(in) :: inpfile, tplfile
        _REAL_, intent(in) :: coords(:, :)
        integer, intent(in) :: natoms
        type(adf_nml_type), intent(in) :: adf_nml
        logical, intent(in) :: do_grad

        character(len=20) :: bakfile

        integer :: temp_atoms(natoms), i, j, ierr, tempvar
        integer, parameter :: iunit = 351
        character(len=20) :: tmp_buffer = ''
        logical, save :: first_call = .true.

        ! 'temp_atoms' is used to write only the unique atoms in the fragments keyword
        temp_atoms = qmmm_struct%iqm_atomic_numbers

        if (adf_nml%verbosity > 0) then
            write (6, *) 'Writing ADF inpfile '//inpfile
        end if

        if (adf_nml%use_template) then
            bakfile = tplfile//'.bak'
            if (first_call) then
                call copy_template(tplfile, trim(bakfile), do_grad)
                call system('cp '//tplfile//' '//inpfile)
            else
                call system('cp '//trim(bakfile)//' '//inpfile)
            end if
        end if

        open (iunit, file=inpfile, iostat=ierr, position='append')
        if (ierr /= 0) then
            call sander_bomb('write_inpfile (qm2_extern_adf_module)', &
                'Error opening ADF inpfile '//inpfile//' for writing', &
                'Will quit now.')
        end if

        ! ATOMS keyword
        write (iunit, '(a)') 'Atoms'
        do i = 1, natoms
            write (iunit, '(a2,1x,3f25.16)') elementSymbol(qmmm_struct%iqm_atomic_numbers(i)), coords(1:3, i)
        end do
        write (iunit, '(a,/)') 'End'

        ! Don't include BASIS/FRAGMENTS etc if we are on the first call
        if (first_call .and. adf_nml%use_template) then
            first_call = .false.
            if (do_grad) write (iunit, '(a)') 'GRADIENT'
            return
        end if

        if (.not. adf_nml%use_dftb) then

            ! ---------------------------
            ! ADF PROGRAM INPUT SPECIFICS
            ! ---------------------------

            ! BASIS/RESTART/FRAGMENTS keyword
            if (first_call) then

                first_call = .false.
                if (adf_nml%fit_type /= 'standard') then
                    tmp_buffer = ' FitType '//adf_nml%fit_type
                end if
                ! write this only in the first call, not needed for restarts
                write (iunit, '(a,/,2(a,a,/),3(a,/))') &
                    'Basis', &
                    ' type ', trim(adf_nml%basis), &
                    ' core ', trim(adf_nml%core), &
                    trim(tmp_buffer), &
                    ' createoutput None', &
                    'End'
            else
                ! Use last t21 file as restart file
                write (iunit, '(a,/,a,/,a,/)') 'Restart adf.t21 &', 'nogeo', 'END'
                ! Set any duplicate atomic numbers in our array to zero
                do i = 1, natoms
                    tempvar = qmmm_struct%iqm_atomic_numbers(i)
                    do j = i + 1, natoms
                        if (temp_atoms(j) == tempvar) then
                            temp_atoms(j) = 0
                        end if
                    end do
                end do
                ! Begin writing fragments
                write (iunit, '(a)') 'Fragments'
                do i = 1, natoms
                    if (temp_atoms(i) /= 0) then
                        write (iunit, '(a,a,a)') elementSymbol(temp_atoms(i)), ' t21.', &
                            elementSymbol(temp_atoms(i))
                    end if
                end do
                write (iunit, '(a,/)') 'End'
            end if

            if (adf_nml%use_template) then
                return
            end if

            ! ELECTRIC field (for external point charges for QM/MM)
            ! AWG: This is disabled at the moment
            ! AWG: We need to add QM region extraction and link atom setup first

            ! XC keyword
            write (iunit, '(a,/,a,/,a,/)') &
                'XC', &
                trim(adf_nml%xc), &
                'END'

            ! SCF keyword
            write (iunit, '(a,/,a,i0,/,a,E22.16,/,a,/)') &
                'SCF ', &
                'iterations ', adf_nml%scf_iter, &
                'converge ', adf_nml%scf_conv, &
                'END'

            ! INTEGRATION keyword
            write (iunit, '(a,E22.16,/)') 'INTEGRATION ', adf_nml%integration

            ! GRADIENT, SYMMETRY, EXACTDENSITY keywords!
            if (adf_nml%oldgradients) then
                write (iunit, '(a)') 'OLDGRADIENTS'
            end if
            ! Need to retreive force data
            if (do_grad) write (iunit, '(a)') 'GRADIENT'
            ! Run without symmetry, we won't need it for MD
            write (iunit, '(a)') 'SYMMETRY NOSYM'
            ! Use exact density for XC calculations if specified
            if (adf_nml%exactdensity) then
                write (iunit, '(a)') 'EXACTDENSITY'
            end if

            ! LINEARSCALING Keyword
            if (adf_nml%linear_scaling /= -1) then
                write (iunit, '(a,i3)') 'LINEARSCALING ', adf_nml%linear_scaling
            end if

            ! CHARGE keyword
            write (iunit, '(a,i2,i2,/)') 'CHARGE ', adf_nml%charge, adf_nml%spin
            if (adf_nml%spin > 0) then
                write (iunit, '(a,/)') 'UNRESTRICTED'
            end if

            ! SAVE files (Need TAPE21 to extract sander data and restart)
            write (iunit, '(a,/)') 'SAVE TAPE21'

        else

            ! ----------------------------
            ! DFTB PROGRAM INPUT SPECIFICS
            ! ----------------------------

            write (iunit, '(a,i2,/)') 'CHARGE ', adf_nml%charge

            ! SCF Keywords
            write (iunit, '(a,/,a,E22.16,/,a,/)') &
                'SCF ', &
                'converge ', adf_nml%scf_conv, &
                'END'

        end if

        ! End writing inpfile for ADF
        close (iunit, iostat=ierr)
        if (ierr /= 0) then
            call sander_bomb('write_inpfile (qm2_extern_adf_module)', &
                'Error closing ADF inpfile for writing', &
                'Will quit now.')
        end if

    end subroutine write_inpfile

    ! ----------------------
    ! Write ADF results
    ! ----------------------
    ! This subroutine will write the dipole moment
    ! to a dipole moment property file
    subroutine write_results(dipfile, dipxyz, dipole)

        use constants, only : CODATA08_AU_TO_DEBYE
        implicit none

        character(len=*), intent(in) :: dipfile
        _REAL_  :: dipxyz(3), dipole
        integer :: iunit = 351, ios
        logical, save :: first_call = .true.

        ! Remove any existing dipole file on first run
        if (first_call) then ! This is set to false later
            call system('rm -f '//dipfile)
        end if
        open (iunit, file=dipfile, position='append', iostat=ios)
        if (ios /= 0) then
            call sander_bomb('write_results (qm2_extern_adf_module)', &
                'Error opening file '//dipfile//' for appending.', &
                'Will quit now')
        end if
        if (first_call) then
            first_call = .false.
            write (iunit, '(a, f22.12)') "Using DEBYE_TO_AU = ", 1.0d0/CODATA08_AU_TO_DEBYE
            write (iunit, '(a)') "| Dipole moment (a.u.-Å): {x, y, z}, |D|"
        end if
        write (iunit, '(4f15.6)') dipxyz(:)/CODATA08_AU_TO_DEBYE, dipole/CODATA08_AU_TO_DEBYE
        close (iunit)

    end subroutine write_results

    subroutine copy_template(tplfile, bakfile, do_grad)

        use UtilitiesModule, only : Upcase

        implicit none
        character(len=*), intent(in) :: tplfile, bakfile
        logical, intent(in) :: do_grad

        integer, parameter :: tplunit = 351, bakunit = 352
        character(len=100) :: read_buffer
        integer :: tplerr, bakerr, ios
        logical :: in_basis = .false.

        open (tplunit, file=tplfile, iostat=tplerr)
        open (bakunit, file=bakfile, iostat=bakerr)
        if (tplerr /= 0) then
            call sander_bomb('copy_template (qm2_extern_adf_module)', &
                'Error opening ADF template file '//tplfile//' for reading', &
                'Will quit now.')
        end if
        if (bakerr /= 0) then
            call sander_bomb('copy_template (qm2_extern_adf_module)', &
                'Error opening ADF template backup file '//bakfile//' for writing', &
                'Will quit now.')
        end if

        ! Write tplfile (without basis key) to bakfile
        do
            read (tplunit, '(a)', iostat=ios) read_buffer
            ! End of file; stop writing
            if (ios < 0) then
                exit
            end if
            if (index(Upcase(read_buffer), 'BASIS') > 0) then
                ! Stop writing until past basis
                in_basis = .true.
            end if
            if (.not. in_basis) then
                write (bakunit, '(a)') read_buffer
            else if (in_basis .and. index(Upcase(read_buffer), 'END') > 0) then
                in_basis = .false.
            end if
        end do

        if (do_grad) write (bakunit, '(a)') 'GRADIENT'

        close (tplunit)
        close (bakunit)

    end subroutine copy_template

end module qm2_extern_adf_module
