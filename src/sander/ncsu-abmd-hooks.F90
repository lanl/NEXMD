!<compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_abmd_hooks

    use ncsu_constants, only : SL => STRING_LENGTH, ABMD_MONITOR_UNIT

    use ncsu_umbrella, only : umbrella_t, &
        MAX_NUMBER_OF_COLVARS => UMBRELLA_MAX_NEXTENTS

    use ncsu_colvar_type, only : colvar_t

#ifdef MPI
    use ncsu_constants, only : ABMD_REMLOG_UNIT
#endif /* MPI */

    implicit none

    private

!-----------------------------------------------------------------------------
!                          T H E    H O O K S
!-----------------------------------------------------------------------------

! multisander.f
#ifdef MPI
    public :: on_delta
    public :: on_exchange
#endif /* MPI */
    public :: on_multisander_exit

!-----------------------------------------------------------------------------

! sander.f
    public :: on_sander_init
    public :: on_sander_exit

! force.f
    public :: on_force

! mdwrit.f
    public :: on_mdwrit

!-----------------------------------------------------------------------------
!                             * P R I V A T E *
!-----------------------------------------------------------------------------

    integer, private, parameter :: MONITOR_UNIT = ABMD_MONITOR_UNIT

    character(*), private, parameter :: SECTION = 'ncsu_abmd'

    character(*), private, parameter :: &
        DEFAULT_MONITOR_FILE = 'ncsu-abmd-monitor', &
        DEFAULT_UMBRELLA_FILE = 'ncsu-abmd-umbrella', &
        DEFAULT_SNAPSHOTS_BASENAME = 'ncsu-abmd-umbrella-snapshot'

    integer, private, parameter :: MODE_NONE = 0

    integer, private, parameter :: MODE_ANALYSIS = 123
    integer, private, parameter :: MODE_UMBRELLA = 234
    integer, private, parameter :: MODE_FLOODING = 345

    integer, private, save :: mode = MODE_NONE

!-----------------------------------------------------------------------------
    integer, private, save :: ncolvars = 0

    type(colvar_t), private, save :: colvars(MAX_NUMBER_OF_COLVARS)
    type(umbrella_t), private, save :: umbrella ! master only
!-----------------------------------------------------------------------------
#ifndef MPI
    NCSU_REAL, private, save :: instantaneous(MAX_NUMBER_OF_COLVARS) ! master only
#else
    NCSU_REAL, private, save :: instantaneous(MAX_NUMBER_OF_COLVARS + 1)
#endif /* MPI */
!-----------------------------------------------------------------------------
    character(len=SL), private, save :: monitor_file ! master only
    character(len=SL), private, save :: monitor_fmt ! master only
    integer, private, save :: monitor_freq
!-----------------------------------------------------------------------------
    character(len=SL), private, save :: snapshots_basename ! master only
    integer, private, save :: snapshots_freq ! master only
!-----------------------------------------------------------------------------
    character(len=SL), private, save :: umbrella_file ! master only
!-----------------------------------------------------------------------------
    NCSU_REAL, private, save :: timescale ! master only
!-----------------------------------------------------------------------------
    integer, private, save :: mdstep ! = runmd.f::nstep + 1
!-----------------------------------------------------------------------------
    integer, private, save :: nhills ! not zeroed on exchange
!-----------------------------------------------------------------------------

#ifdef MPI

!
! multiwalk/rem
!

    NCSU_REAL, private, allocatable, save :: all_hills(:)

!
! rem-specific
!

    character(*), private, parameter :: REMLOG_FILE = 'ncsu-abmd.log'
    integer, private, parameter :: REMLOG_UNIT = ABMD_REMLOG_UNIT

! rem_* arrays live on sander masters; indexed by masterrank

    character(len=SL), private, allocatable, save :: rem_monitor_files(:)
    integer, private, allocatable, save :: rem_monitor_freqs(:)

    character(len=SL), private, allocatable, save :: rem_umbrella_files(:)
    type(umbrella_t), private, allocatable, save :: rem_umbrellas(:)

    character(len=SL), private, allocatable, save :: rem_snapshots_basenames(:)
    integer, private, allocatable, save :: rem_snapshots_freqs(:)

    private :: rem_postinit
    private :: rem_cleanup

    integer, private, save :: LOG_UNIT
    integer, private, save :: exchange_partner

#else
#  define LOG_UNIT OUT_UNIT
#endif /* MPI */

    NCSU_REAL, parameter, private :: TINY = 0.00001d0 ! NCSU_TO_REAL(0.00001)

!-----------------------------------------------------------------------------

contains

!-----------------------------------------------------------------------------

#ifdef MPI
    subroutine on_delta(o_masterrank, need_U_xx, U_mm, U_mo, U_om, U_oo)

        NCSU_USE_AFAILED

        use ncsu_umbrella
        use ncsu_sander_proxy

        implicit none

        integer, intent(in) :: o_masterrank
        logical, intent(in) :: need_U_xx

        NCSU_REAL, intent(inout) :: U_mm, U_mo, U_om, U_oo

#  include "ncsu-mpi.h"

        NCSU_REAL :: o_instantaneous(MAX_NUMBER_OF_COLVARS), U_o(2), U_m(2)

        integer :: error

        if (mode .eq. MODE_NONE .or. mode .eq. MODE_ANALYSIS) &
            return

        ncsu_assert(multisander_rem() .eq. 1)
        ncsu_assert(sanderrank .eq. 0) ! master
        ncsu_assert(commmaster .ne. mpi_comm_null)

        ! exchange instantaneous(:) with the partner
        call mpi_sendrecv &
            (instantaneous, ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
            o_instantaneous, ncolvars, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
            commmaster, MPI_STATUS_IGNORE, error)
        ncsu_assert(error .eq. 0)

        ! evaluate 'my' values
        U_m(1) = umbrella_eval_v(umbrella, instantaneous)   ! U_mm = U_m(x_m)
        U_m(2) = umbrella_eval_v(umbrella, o_instantaneous) ! U_mo = U_m(x_o)

        ! get partner's U_m? (i.e., U_o? in this replica)
        call mpi_sendrecv &
            (U_m, 2, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
            U_o, 2, MPI_DOUBLE_PRECISION, o_masterrank, mdstep, &
            commmaster, MPI_STATUS_IGNORE, error)
        ncsu_assert(error .eq. 0)

        if (need_U_xx) then
            U_mm = U_mm + U_m(1)
            U_mo = U_mo + U_m(2)
            U_om = U_om + U_o(2)
            U_oo = U_oo + U_o(1)
        end if

    end subroutine on_delta

!-----------------------------------------------------------------------------

    subroutine on_exchange(o_masterrank)

#ifndef NCSU_DISABLE_ASSERT
        use ncsu_utils
        use ncsu_sander_proxy
#endif /* NCSU_DISABLE_ASSERT */

        implicit none

        integer, intent(in) :: o_masterrank

#  include "ncsu-mpi.h"

        if (mode .eq. MODE_NONE) &
            return

        ncsu_assert(multisander_rem() .eq. 1)
        ncsu_assert(sanderrank .eq. 0) ! master
        ncsu_assert(commmaster .ne. mpi_comm_null)

        exchange_partner = o_masterrank

    end subroutine on_exchange
#endif /* MPI */

!-----------------------------------------------------------------------------

    subroutine on_multisander_exit()

        NCSU_USE_AFAILED

        use ncsu_colvar
        use ncsu_umbrella
        use ncsu_sander_proxy

        implicit none

#  include "ncsu-mpi.h"

        integer :: n

        if (mode .ne. MODE_NONE) then
            ncsu_assert(ncolvars .gt. 0)
            do n = 1, ncolvars
                call colvar_cleanup(colvars(n))
            end do
        end if

        NCSU_MASTER_ONLY_BEGIN

#ifdef MPI
        if (allocated(all_hills)) &
            deallocate (all_hills)

        if (multisander_rem() .eq. 1) &
            call rem_cleanup()
#endif /* MPI */

        if (mode .eq. MODE_FLOODING .or. mode .eq. MODE_UMBRELLA) &
            call umbrella_fini(umbrella)

        NCSU_MASTER_ONLY_END

        mode = MODE_NONE

    end subroutine on_multisander_exit

!-----------------------------------------------------------------------------

!
! 'on_sander_init()' is called at the point when
! MPI is initialized and MDIN/PRMTOP/INPCRD are loaded
!

#ifdef MPI
    subroutine on_sander_init(mdin_name, mdin_unit, amass, rem_idx)
#else
    subroutine on_sander_init(mdin_name, mdin_unit, amass)
#endif /* MPI */

        use ncsu_utils
        use ncsu_cftree
        use ncsu_colvar
        use ncsu_parser
        use ncsu_constants
        use ncsu_umbrella
        use ncsu_sander_proxy

        implicit none

        character(SL), intent(in) :: mdin_name
        integer, intent(in) :: mdin_unit

        NCSU_REAL, intent(in) :: amass(*)

#ifdef MPI
        integer, intent(in) :: rem_idx
#endif /* MPI */

        type(node_t), pointer :: root
        type(child_t), pointer :: child

        logical :: umbrella_file_exists
        type(umbrella_t) :: umbrella_from_file

        logical :: found, do_transfer
        integer :: n, error

        integer :: cv_extents(UMBRELLA_MAX_NEXTENTS)
        logical :: cv_periodicity(UMBRELLA_MAX_NEXTENTS)

        NCSU_REAL :: cv_origin(UMBRELLA_MAX_NEXTENTS)
        NCSU_REAL :: cv_spacing(UMBRELLA_MAX_NEXTENTS)

        NCSU_REAL :: tmp

        character(len=SL) :: astring

#ifdef MPI
#  include "ncsu-mpi.h"
        integer, allocatable :: perm(:)
        integer, save :: counter = 0
#endif /* MPI */

        mdstep = 0 ! used in on_force() below

#ifdef MPI
        if (counter .gt. 0) then
            ncsu_assert(multisander_rem() .gt. 0)

            NCSU_MASTER_ONLY_BEGIN
            if (mode .ne. MODE_NONE) then

                ! re-arrange rem_* arrays

                allocate (perm(mastersize), stat=error)
                if (error .ne. 0) &
                    NCSU_OUT_OF_MEMORY

                call mpi_allgather(exchange_partner, 1, MPI_INTEGER, &
                    perm, 1, MPI_INTEGER, commmaster, error)
                ncsu_assert(error .eq. 0)

                if (mode .ne. MODE_ANALYSIS) &
                    call umbrella_swap(umbrella, rem_umbrellas(masterrank))

                do n = 0, mastersize - 1
                    if (n .gt. perm(n + 1)) then

                        call swap(rem_monitor_files(n), rem_monitor_files(perm(n + 1)))
                        call swap(rem_monitor_freqs(n), rem_monitor_freqs(perm(n + 1)))

                        if (mode .eq. MODE_FLOODING) then
                            call swap(rem_snapshots_freqs(n), &
                                rem_snapshots_freqs(perm(n + 1)))
                            call swap(rem_snapshots_basenames(n), &
                                rem_snapshots_basenames(perm(n + 1)))
                            call swap(rem_umbrella_files(n), &
                                rem_umbrella_files(perm(n + 1)))
                        end if ! mode.eq.MODE_FLOODING

                        if (mode .ne. MODE_ANALYSIS) &
                            call umbrella_swap(rem_umbrellas(n), &
                            rem_umbrellas(perm(n + 1)))
                    end if ! n.gt.perm(n + 1)
                    exchange_partner = masterrank
                end do

                deallocate (perm)

                ! adjust filenames/freqs/umbrellas

                monitor_file = rem_monitor_files(masterrank)
                monitor_freq = rem_monitor_freqs(masterrank)

                if (mode .eq. MODE_FLOODING) then
                    snapshots_freq = rem_snapshots_freqs(masterrank)
                    snapshots_basename = rem_snapshots_basenames(masterrank)
                    umbrella_file = rem_umbrella_files(masterrank)
                end if

                if (mode .ne. MODE_ANALYSIS) &
                    call umbrella_swap(umbrella, rem_umbrellas(masterrank))

                ! re-open MONITOR_UNIT after exchange (closed in on_sander_exit())
                open (unit=MONITOR_UNIT, file=monitor_file, &
                    iostat=error, form='FORMATTED', action='WRITE', &
                    position='APPEND', status='OLD')

                if (error .ne. 0) then
                    write (unit=ERR_UNIT, fmt='(/a,a,a,a/)') NCSU_ERROR, &
                        'could not open ''', trim(monitor_file), ''' file for writing'
                    call terminate()
                end if

            end if ! mode.ne.MODE_NONE
            NCSU_MASTER_ONLY_END

            if (mode .eq. MODE_ANALYSIS) then
                call mpi_bcast(monitor_freq, 1, MPI_INTEGER, 0, commsander, error)
                ncsu_assert(error .eq. 0)
            end if

            return

        end if ! counter.gt.0

        counter = counter + 1

        NCSU_MASTER_ONLY_BEGIN

        if (multisander_rem() .ne. 0) then
            ncsu_assert(multisander_numgroup() .gt. 1)
            exchange_partner = masterrank
            n = rem_idx
        else if (multisander_numgroup() .gt. 1) then
            n = masterrank + 1
        end if

        if (multisander_numgroup() .gt. 1) then
            write (unit=monitor_file, fmt='(a,a,i3.3)') &
                DEFAULT_MONITOR_FILE, '-', n
            write (unit=umbrella_file, fmt='(a,a,i3.3,a)') &
                DEFAULT_UMBRELLA_FILE, '-', n, '.nc'
            write (unit=snapshots_basename, fmt='(a,a,i3.3)') &
                DEFAULT_SNAPSHOTS_BASENAME, '-', n
        else
            monitor_file = DEFAULT_MONITOR_FILE
            umbrella_file = DEFAULT_UMBRELLA_FILE//'.nc'
            snapshots_basename = DEFAULT_SNAPSHOTS_BASENAME
        end if
#else
        monitor_file = DEFAULT_MONITOR_FILE
        umbrella_file = DEFAULT_UMBRELLA_FILE
        snapshots_basename = DEFAULT_SNAPSHOTS_BASENAME
#endif /* MPI */

        ! parse MDIN (on master)
        ncsu_assert(ncolvars .eq. 0)
        root => parse_cf(mdin_name, SECTION, mdin_unit)

        ! no root section
        if (.not. associated(root)) then
            mode = MODE_NONE
            goto 2
        end if

        ! discover the run-mode
        found = node_lookup_string(root, 'mode', astring)
        if (.not. found) &
            call fatal('could not find ''mode'' in the '''//SECTION//''' section')

        if (astring == 'NONE') then
            mode = MODE_NONE
            goto 1
        else if (astring == 'ANALYSIS') then
            mode = MODE_ANALYSIS
        else if (astring == 'UMBRELLA') then
            mode = MODE_UMBRELLA
        else if (astring == 'FLOODING') then
            mode = MODE_FLOODING
        else
            write (unit=ERR_UNIT, fmt='(/a,a,a,a/)') &
                NCSU_ERROR, 'unknown mode ''', trim(astring), ''''
            call terminate()
        end if

        found = node_lookup_string(root, 'umbrella_file', astring)
        if (found) &
            umbrella_file = astring

#ifdef NCSU_NO_NETCDF
        umbrella_file_exists = .false.
        write (unit=ERR_UNIT, fmt='(a,a)') NCSU_WARNING, &
            'netCDF is not available (try ''-bintraj'' configure option)'
#else
        inquire (file=umbrella_file, exist=umbrella_file_exists)
#endif /* NCSU_NO_NETCDF */

        if (.not. umbrella_file_exists .and. mode .eq. MODE_UMBRELLA) then
            write (unit=ERR_UNIT, fmt='(/a,a,a,a/)') NCSU_ERROR, '''', &
                trim(umbrella_file), ''' does not exist (required for UMBRELLA mode)'
            call terminate()
        end if

        if (mode .eq. MODE_ANALYSIS) &
            umbrella_file_exists = .false.

#ifndef NCSU_NO_NETCDF
        if (umbrella_file_exists) &
            call umbrella_load(umbrella_from_file, umbrella_file)
#endif /* NCSU_NO_NETCDF */

        ! collective variables
        ncsu_assert(ncolvars .eq. 0)

        child => node_children(root)
        do while (associated(child))
            if (node_title(child%node) == 'variable') &
                ncolvars = ncolvars + 1
            child => child%next
        end do

        if (ncolvars .eq. 0) &
            call fatal('no variable(s) in the '''//SECTION//''' section')

        if (ncolvars .gt. MAX_NUMBER_OF_COLVARS) &
            call fatal('too many variables in the '''//SECTION//''' section')

        if (umbrella_file_exists) then
            if (umbrella_nextents(umbrella_from_file) .ne. ncolvars) &
                call fatal('number of variables in the '''//SECTION//''' does not &
            &match with the number of extents found in the umbrella_file')
        end if ! umbrella_file_exists

        n = 1
        child => node_children(root)

        do while (associated(child))
            if (node_title(child%node) == 'variable') then

                ncsu_assert(n <= ncolvars)
                call colvar_mdread(colvars(n), child%node, n)

                if (mode .eq. MODE_FLOODING) then
                    found = node_lookup_positive_real(child%node, &
                        'resolution', cv_spacing(n))
                    cv_spacing(n) = cv_spacing(n)/4
                    if (.not. found .and. .not. umbrella_file_exists) then
                        write (unit=ERR_UNIT, fmt='(/a,a,i1/)') NCSU_ERROR, &
                            'could not determine ''resolution'' for CV #', n
                        call terminate()
                    end if

                    if (.not. found) &
                        cv_spacing(n) = umbrella_spacing(umbrella_from_file, n)

                    cv_periodicity(n) = colvar_is_periodic(colvars(n))

                    if (cv_periodicity(n)) then

                        ncsu_assert(colvar_has_min(colvars(n)))
                        ncsu_assert(colvar_has_max(colvars(n)))

                        cv_origin(n) = colvar_min(colvars(n))

                        ncsu_assert(cv_spacing(n) .gt. ZERO)
                        ncsu_assert(colvar_max(colvars(n)) .gt. cv_origin(n))

                        cv_extents(n) = &
                            int((colvar_max(colvars(n)) - cv_origin(n))/cv_spacing(n))

                        if (cv_extents(n) .lt. UMBRELLA_MIN_EXTENT) then
                            write (unit=ERR_UNIT, fmt='(/a,a,i1,a/)') NCSU_ERROR, &
                                'CV #', n, ' : ''resolution'' is too big'
                            call terminate()
                        end if

                        cv_spacing(n) = &
                            (colvar_max(colvars(n)) - cv_origin(n))/cv_extents(n)

                    else ! .not.periodic

                        found = node_lookup_real(child%node, 'min', cv_origin(n))
                        if (.not. found) then
                            if (umbrella_file_exists) then
                                cv_origin(n) = umbrella_origin(umbrella_from_file, n)
                            else if (colvar_has_min(colvars(n))) then
                                cv_origin(n) = colvar_min(colvars(n))
                            else
                                write (unit=ERR_UNIT, fmt='(/a,a,i1/)') NCSU_ERROR, &
                                    'could not determine ''min'' for CV #', n
                                call terminate()
                            end if
                        end if

                        found = node_lookup_real(child%node, 'max', tmp)
                        if (.not. found) then
                            if (umbrella_file_exists) then
                                tmp = umbrella_origin(umbrella_from_file, n) &
                                    + umbrella_spacing(umbrella_from_file, n) &
                                    *(umbrella_extent(umbrella_from_file, n) - 1)
                            else if (colvar_has_max(colvars(n))) then
                                tmp = colvar_max(colvars(n))
                            else
                                write (unit=ERR_UNIT, fmt='(/a,a,i1/)') NCSU_ERROR, &
                                    'could not determine ''max'' for CV #', n
                                call terminate()
                            end if
                        end if

                        if (cv_origin(n) .ge. tmp) then
                            write (unit=ERR_UNIT, fmt='(/a,a,i1/)') NCSU_ERROR, &
                                'min.ge.max for CV #', n
                            call terminate()
                        end if

                        ncsu_assert(cv_spacing(n) .gt. ZERO)
                        cv_extents(n) = 1 &
                            + int((tmp - cv_origin(n))/cv_spacing(n))

                        if (cv_extents(n) .lt. UMBRELLA_MIN_EXTENT) then
                            write (unit=ERR_UNIT, fmt='(/a,a,i1,a/)') NCSU_ERROR, &
                                'CV #', n, ' : the ''resolution'' is too big'
                            call terminate()
                        end if

                        cv_spacing(n) = (tmp - cv_origin(n))/(cv_extents(n) - 1)

                    end if ! cv_periodicity(n)
                end if ! mode.eq.MODE_FLOODING
                n = n + 1
            end if ! title == 'variable'
            child => child%next
        end do

        ! monitor
        found = node_lookup_string(root, 'monitor_file', astring)
        if (found) &
            monitor_file = astring

        found = node_lookup_positive_integer(root, 'monitor_freq', monitor_freq)
        if (.not. found) &
            monitor_freq = 50

        monitor_freq = min(monitor_freq, sander_nstlim())
        monitor_freq = max(1, monitor_freq)

        ! umbrella snapshots
        found = node_lookup_string(root, 'snapshots_basename', astring)
        if (found) &
            snapshots_basename = astring

        found = node_lookup_integer(root, 'snapshots_freq', snapshots_freq)
        if (.not. found) &
            snapshots_freq = -1 ! no snapshots

        if (mode .eq. MODE_FLOODING) then
            found = node_lookup_positive_real(root, 'timescale', timescale)
            if (.not. found) &
                call fatal('could not find ''timescale'' &
            &in the '''//SECTION//''' section')

        end if ! mode.eq.MODE_FLOODING

1       call node_cleanup(root)
        deallocate (root)

2       continue ! done with parsing

        NCSU_MASTER_ONLY_END

#ifdef MPI
        call mpi_bcast(mode, 1, MPI_INTEGER, 0, commsander, error)
        ncsu_assert(error .eq. 0)
#endif /* MPI */

        if (mode .eq. MODE_NONE) &
            return

#ifdef MPI
        ncsu_assert(.not. is_master() .or. ncolvars .gt. 0)

        call mpi_bcast(ncolvars, 1, MPI_INTEGER, 0, commsander, error)
        ncsu_assert(error .eq. 0)

        call mpi_bcast(monitor_freq, 1, MPI_INTEGER, 0, commsander, error)
        ncsu_assert(error .eq. 0)
#endif /* MPI */

        ncsu_assert(ncolvars .gt. 0)
        ncsu_assert(ncolvars .le. MAX_NUMBER_OF_COLVARS)

        if (multisander_numwatkeep() .gt. 0) &
            call fatal('numwatkeep.gt.0 is not supported')

        if (sander_imin() .ne. 0) &
            call fatal('imin.ne.0 is not supported')

        do n = 1, ncolvars
            call colvar_bootstrap(colvars(n), n, amass)
        end do

#ifdef MPI
        if (multisander_numgroup() .gt. 1) &
            call rem_checks()
#endif /* MPI */

        NCSU_MASTER_ONLY_BEGIN

        if (mode .eq. MODE_UMBRELLA) then
            ncsu_assert(umbrella_file_exists)
            do_transfer = .false.
            call umbrella_swap(umbrella, umbrella_from_file)
        else if (mode .eq. MODE_FLOODING) then
            if (umbrella_file_exists) then
                do_transfer = .false.
                do n = 1, ncolvars
                    do_transfer = do_transfer &
                        .or. (cv_extents(n) .ne. umbrella_extent(umbrella_from_file, n))
                    do_transfer = do_transfer &
                        .or. (cv_periodicity(n) .neqv. &
                        umbrella_periodicity(umbrella_from_file, n))
                    do_transfer = do_transfer &
                        .or. (abs(cv_origin(n) - umbrella_origin(umbrella_from_file, n)) &
                        .gt. TINY)
                    do_transfer = do_transfer &
                        .or. (abs(cv_spacing(n) - umbrella_spacing(umbrella_from_file, &
                        n)) .gt. TINY)
                    if (do_transfer) &
                        exit
                end do
                if (do_transfer) then
                    call umbrella_init(umbrella, ncolvars, cv_extents, &
                        cv_origin, cv_spacing, cv_periodicity)
                    call umbrella_transfer(umbrella, umbrella_from_file)
                    call umbrella_fini(umbrella_from_file)
                else
                    call umbrella_swap(umbrella, umbrella_from_file)
                end if ! do_transfer
            else
                call umbrella_init(umbrella, ncolvars, cv_extents, &
                    cv_origin, cv_spacing, cv_periodicity)
            end if ! umbrella_file_exits
#ifdef MPI
            if (multisander_numgroup() .gt. 1) then
                allocate (all_hills((ncolvars + 1)*multisander_numgroup()), &
                    stat=error)
                if (error .ne. 0) &
                    NCSU_OUT_OF_MEMORY
                if (multisander_rem() .eq. 0) then
                    if (masterrank .gt. 0) &
                        call umbrella_fini(umbrella)

                    call umbrella_bcast(umbrella, commmaster, 0)

                    call mpi_bcast(timescale, 1, &
                        MPI_DOUBLE_PRECISION, 0, commmaster, error)
                    ncsu_assert(error .eq. 0)
                end if ! .not.REM
            end if ! ng.gt.1
#endif /* MPI */
        end if

        ! prepare monitor_fmt & open MONITOR_UNIT

        open (unit=MONITOR_UNIT, file=monitor_file, iostat=error, &
            form='FORMATTED', action='WRITE', status='REPLACE')

        if (error .ne. 0) then
            write (unit=ERR_UNIT, fmt='(/a,a,a,a/)') &
                NCSU_ERROR, 'failed to open ''', trim(monitor_file), ''' for writing'
            call terminate()
        end if

        write (unit=MONITOR_UNIT, fmt='(a,/a)', advance='NO') &
            '#', '# MD time (ps), '
        do n = 1, ncolvars - 1
            write (unit=MONITOR_UNIT, fmt='(a,i1,a)', advance='NO') &
                'CV #', n, ', '
        end do

        if (mode == MODE_FLOODING) then
            write (unit=MONITOR_UNIT, fmt='(a,i1,a,/a)') &
                'CV #', ncolvars, ', E_{bias} (kcal/mol)', '#'
            write (unit=monitor_fmt, fmt='(a,i1,a)') &
                '(f12.4,', ncolvars, '(1x,f16.10),1x,f16.10)'
        else
            write (unit=MONITOR_UNIT, fmt='(a,i1,/a)') &
                'CV #', ncolvars, '#'
            write (unit=monitor_fmt, fmt='(a,i1,a)') &
                '(f12.4,', ncolvars, '(1x,f16.10))'
        end if

        call flush_UNIT(MONITOR_UNIT)

        ! print summary & return

#ifdef MPI
        if (multisander_rem() .eq. 0) then
            LOG_UNIT = OUT_UNIT ! write to MDOUT
        else
            LOG_UNIT = REMLOG_UNIT

            if (masterrank .eq. 0) then
                open (unit=LOG_UNIT, file=REMLOG_FILE, iostat=error, &
                    form='FORMATTED', action='WRITE', status='REPLACE')

                if (error .ne. 0) then
                    write (unit=ERR_UNIT, fmt='(/a,a,a,a/)') NCSU_ERROR, &
                        'failed to open ''', trim(REMLOG_FILE), ''' for writing'
                    call terminate()
                end if
                write (unit=LOG_UNIT, fmt='(80(''*''))')
                close (unit=LOG_UNIT)
            end if ! masterrank.eq.0

            call cpus_enter(commmaster, 12345)

            open (unit=LOG_UNIT, file=REMLOG_FILE, &
                iostat=error, form='FORMATTED', action='WRITE', &
                position='APPEND', status='OLD')

            if (error .ne. 0) then
                write (unit=ERR_UNIT, fmt='(/a,a,a,a/)') NCSU_ERROR, &
                    'could not open ''', trim(REMLOG_FILE), ''' file for writing'
                call terminate()
            end if

            write (unit=LOG_UNIT, fmt= &
                '(/23x,a,'//pfmt(masterrank)//',a,'//pfmt(sander_temp0(), 3)//',a/)') &
                'REPLICA #', masterrank, ' (temp0 = ', sander_temp0(), ')'
            write (unit=LOG_UNIT, fmt='(1x,a,a,a)') &
                'MDIN = ''', trim(mdin_name), ''''

        end if
#else
#  define LOG_UNIT OUT_UNIT
#endif
        write (unit=LOG_UNIT, fmt='(/a,a)') NCSU_INFO, &
            '() () () () () () () () () ()   A. B. M. D.  () () () () () () () () ()'
        write (unit=LOG_UNIT, fmt='(a,/a,a)', advance='NO') &
            NCSU_INFO, NCSU_INFO, 'mode = '

        select case (mode)
          case (MODE_ANALYSIS)
            write (unit=LOG_UNIT, fmt='(a)') 'ANALYSIS'
          case (MODE_UMBRELLA)
            write (unit=LOG_UNIT, fmt='(a)') 'UMBRELLA'
          case (MODE_FLOODING)
            write (unit=LOG_UNIT, fmt='(a)') 'FLOODING'
          case default
            ncsu_assert_not_reached()
            continue
        end select

        write (unit=LOG_UNIT, fmt='(a)') NCSU_INFO

        do n = 1, ncolvars
            write (unit=LOG_UNIT, fmt='(a,a,i1)') NCSU_INFO, 'CV #', n
            call colvar_print(colvars(n), LOG_UNIT)
            write (unit=LOG_UNIT, fmt='(a)') NCSU_INFO
        end do

        write (unit=LOG_UNIT, fmt='(a,a,a)') NCSU_INFO, &
            'monitor_file = ', trim(monitor_file)
        write (unit=LOG_UNIT, &
            fmt='(a,a,'//pfmt(monitor_freq)//',a,'//pfmt( &
            monitor_freq*sander_timestep(), 4)//',a)') NCSU_INFO, &
            'monitor_freq = ', monitor_freq, ' (', &
            monitor_freq*sander_timestep(), ' ps)'

        if (mode .eq. MODE_ANALYSIS) &
            goto 3

#  ifdef MPI
        if (multisander_numgroup() .gt. 1 .and. multisander_rem() .eq. 0) then
            if (masterrank .gt. 0) &
                write (unit=LOG_UNIT, fmt='(a,/a,a,/a)') NCSU_INFO, NCSU_INFO, &
                'ng.gt.1.and.rem.eq.0 => using umbrella from replica #1', NCSU_INFO
        end if
#  endif /* MPI */

        write (unit=LOG_UNIT, fmt='(a,a,a,a)', advance='NO') NCSU_INFO, &
            'umbrella_file = ', trim(umbrella_file), ' ('

        if (umbrella_file_exists) then
            write (unit=LOG_UNIT, fmt='(a)') 'loaded)'
        else
            write (unit=LOG_UNIT, fmt='(a)') 'not found)'
        end if

        write (unit=LOG_UNIT, fmt='(a)') NCSU_INFO
        write (unit=LOG_UNIT, fmt='(a,a)', advance='NO') NCSU_INFO, &
            'umbrella discretization '

        if (umbrella_file_exists) then
            if (do_transfer) then
                write (unit=LOG_UNIT, fmt='(a)') '(modified) :'
            else
                write (unit=LOG_UNIT, fmt='(a)') '(unchanged) :'
            end if
        else
            write (unit=LOG_UNIT, fmt='(a)') '(new) :'
        end if

        do n = 1, ncolvars
            write (unit=LOG_UNIT, fmt='(a,a,i1)', advance='NO') &
                NCSU_INFO, 'CV #', n
            if (umbrella_periodicity(umbrella, n)) then
                write (unit=LOG_UNIT, fmt='(a)', advance='NO') ' periodic, '
                tmp = umbrella_origin(umbrella, n) &
                    + umbrella_spacing(umbrella, n)*umbrella_extent(umbrella, n)
            else
                write (unit=LOG_UNIT, fmt='(a)', advance='NO') ' not periodic, '
                tmp = umbrella_origin(umbrella, n) &
                    + umbrella_spacing(umbrella, n)*(umbrella_extent(umbrella, n) - 1)
            end if

            write (unit=LOG_UNIT, &
                fmt='('//pfmt(umbrella_extent(umbrella, &
                n))//',a,'//pfmt(umbrella_origin(umbrella, n), &
                6)//',a,'//pfmt(tmp, 6)//')') &
                umbrella_extent(umbrella, n), ' points, min/max = ', &
                umbrella_origin(umbrella, n), '/', tmp
        end do

        if (mode .eq. MODE_UMBRELLA) &
            goto 3

        write (unit=LOG_UNIT, fmt='(a/,a,a,'//pfmt(timescale, 3)//',a)') &
            NCSU_INFO, NCSU_INFO, 'flooding timescale = ', timescale, ' ps'

        if (snapshots_freq .gt. 0) then
            write (unit=LOG_UNIT, fmt='(a,a,a)') NCSU_INFO, &
                'snapshots_basename = ', trim(snapshots_basename)
            write (unit=LOG_UNIT, &
                fmt='(a,a,'//pfmt(snapshots_freq)//',a,'//pfmt(snapshots_freq*sander_timestep(), 4)//',a)') &
                NCSU_INFO, 'snapshots_freq = ', snapshots_freq, ' (', &
                snapshots_freq*sander_timestep(), ' ps)'
        end if

        nhills = 0

3       write (unit=LOG_UNIT, fmt='(a)') NCSU_INFO
        write (unit=LOG_UNIT, fmt='(a,a/)') NCSU_INFO, &
            '() () () () () () () () () () () () () () () () () () () () () () () ()'
        call flush_UNIT(LOG_UNIT)

#ifdef MPI

        if (multisander_rem() .ne. 0) then
            write (unit=LOG_UNIT, fmt='(/80(''*''))')
            close (unit=LOG_UNIT)
            call cpus_leave(commmaster, 12345)
            call rem_postinit() ! allocates/populates rem_* arrays
        end if

        NCSU_MASTER_ONLY_END

    contains

!.............................................................................

        subroutine rem_checks()

            implicit none

            integer, allocatable :: int_recv(:)

            if (commmaster .eq. mpi_comm_null) &
                return

            ncsu_assert(mastersize .gt. 0)
            ncsu_assert(masterrank .lt. mastersize)

            allocate (int_recv(mastersize), stat=error)
            if (error .ne. 0) &
                NCSU_OUT_OF_MEMORY

            call mpi_allgather(mode, 1, MPI_INTEGER, int_recv, 1, MPI_INTEGER, &
                commmaster, error)
            ncsu_assert(error .eq. 0)

            do n = 1, mastersize
                if (int_recv(n) .ne. mode) &
                    call fatal('''mode'' has different values in different replicas')
            end do

            call mpi_allgather(ncolvars, 1, MPI_INTEGER, int_recv, 1, MPI_INTEGER, &
                commmaster, error)
            ncsu_assert(error .eq. 0)

            do n = 1, mastersize
                if (int_recv(n) .ne. ncolvars) &
                    call fatal('number of collective variables is &
                &different in different replicas')
            end do

            ! FIXME: more here

            deallocate (int_recv)

        end subroutine rem_checks
#endif /* MPI */

    end subroutine on_sander_init

!-----------------------------------------------------------------------------

    subroutine on_sander_exit()

        use ncsu_utils, only : close_UNIT

        implicit none

#  include "ncsu-mpi.h"

        NCSU_MASTER_ONLY_BEGIN
        call close_UNIT(MONITOR_UNIT)
        NCSU_MASTER_ONLY_END

    end subroutine on_sander_exit

!-----------------------------------------------------------------------------

    subroutine on_force(x, f, virial)

        use ncsu_utils
        use ncsu_colvar
        use ncsu_umbrella
        use ncsu_constants
        use ncsu_sander_proxy

        implicit none

        NCSU_REAL, intent(in) :: x(*)

        NCSU_REAL, intent(inout) :: f(*)
        NCSU_REAL, intent(inout) :: virial(4)

#ifdef MPI
#  include "ncsu-mpi.h"
        integer :: error
#endif /* MPI */

        NCSU_REAL :: u_value, u_derivative(UMBRELLA_MAX_NEXTENTS), alt

        character(len=SL + 16) :: snapshot

        integer :: n
        logical :: real_mdstep

        if (mode .eq. MODE_NONE) &
            return

        real_mdstep = (.not. multisander_initremd()) .and. (sander_init() .eq. 4)

        ncsu_assert(ncolvars .gt. 0)

        if (mode .eq. MODE_ANALYSIS) then
            if (real_mdstep .and. mod(mdstep - 1, monitor_freq) .eq. 0) then
                do n = 1, ncolvars
                    instantaneous(n) = colvar_value(colvars(n), x)
                end do

                NCSU_MASTER_ONLY_BEGIN
                write (unit=MONITOR_UNIT, fmt=monitor_fmt) &
                    sander_mdtime(), instantaneous(1:ncolvars)
                call flush_UNIT(MONITOR_UNIT)
                NCSU_MASTER_ONLY_END
            end if
            goto 1
        end if ! mode.eq.MODE_ANALYSIS

        !
        ! either UMBRELLA or FLOODING
        !

        do n = 1, ncolvars
            instantaneous(n) = colvar_value(colvars(n), x)
        end do

        NCSU_MASTER_ONLY_BEGIN
        call umbrella_eval_vdv(umbrella, instantaneous, u_value, u_derivative)
        NCSU_MASTER_ONLY_END

#ifdef MPI
        call mpi_bcast(u_derivative, ncolvars, &
            MPI_DOUBLE_PRECISION, 0, commsander, error)
        ncsu_assert(error .eq. 0)
#endif /* MPI */

        ! FIXME: virial
        do n = 1, ncolvars
            call colvar_force(colvars(n), x, -u_derivative(n), f)
        end do

        if (real_mdstep) then
            NCSU_MASTER_ONLY_BEGIN

            ! mdstep.ge.monitor_freq - 1 is to avoid the sample right after exchange
            if ((mdstep + 1) .ge. monitor_freq &
                .and. mod(mdstep + 1, monitor_freq) .eq. 0) then
                if (mode .eq. MODE_FLOODING) then
                    write (unit=MONITOR_UNIT, fmt=monitor_fmt) &
                        sander_mdtime(), instantaneous(1:ncolvars), u_value
                else
                    write (unit=MONITOR_UNIT, fmt=monitor_fmt) &
                        sander_mdtime(), instantaneous(1:ncolvars)
                end if
                call flush_UNIT(MONITOR_UNIT)
            end if

            if (mode .eq. MODE_FLOODING) then
#ifndef NCSU_NO_NETCDF
                if (snapshots_freq .gt. 0 .and. mod(nhills, snapshots_freq) .eq. 0) then
                    write (unit=snapshot, fmt='(a,a,i10.10,a)') &
                        trim(snapshots_basename), '.', nhills, '.nc'
                    call umbrella_save(umbrella, snapshot)
                    write (unit=OUT_UNIT, fmt='(/a,a,f16.4,a,/a,a,a,a)') &
                        NCSU_INFO, 'biasing potential snapshot at t = ', &
                        sander_mdtime(), ' ps', NCSU_INFO, 'saved as ''', &
                        trim(snapshot), ''''
                end if
#endif /* NCSU_NO_NETCDF */
                alt = sander_timestep()/timescale
#ifdef MPI
                if (multisander_numgroup() .gt. 1) then
                    ncsu_assert(commmaster .ne. mpi_comm_null)
                    ! get all instantaneous/altitudes
                    ncsu_assert(allocated(all_hills))
                    instantaneous(ncolvars + 1) = alt
                    call mpi_allgather(instantaneous, ncolvars + 1, &
                        MPI_DOUBLE_PRECISION, all_hills, ncolvars + 1, &
                        MPI_DOUBLE_PRECISION, commmaster, error)
                    ncsu_assert(error .eq. 0)

                    if (multisander_rem() .eq. 0) then
                        do n = 0, multisander_numgroup() - 1
                            call umbrella_hill(umbrella, &
                                all_hills(n*(ncolvars + 1) + 1:), &
                                all_hills((n + 1)*(ncolvars + 1)))
                        end do
                    else if (multisander_rem() .eq. 1) then
                        do n = 0, multisander_numgroup() - 1
                            if (n .eq. masterrank) then
                                call umbrella_hill(umbrella, &
                                    all_hills(n*(ncolvars + 1) + 1:), &
                                    all_hills((n + 1)*(ncolvars + 1)))
                            else
                                call umbrella_hill(rem_umbrellas(n), &
                                    all_hills(n*(ncolvars + 1) + 1:), &
                                    all_hills((n + 1)*(ncolvars + 1)))
                            end if ! n.eq.masterrank
                        end do
                    else
                        continue
                        ncsu_assert_not_reached()
                    end if
                else
#endif /* MPI */
                    call umbrella_hill(umbrella, instantaneous, alt)
#ifdef MPI
                end if ! multisander_numgroup().gt.1
#endif /* MPI */
                nhills = nhills + 1
            end if
            NCSU_MASTER_ONLY_END

        end if ! real_mdstep

1       if (real_mdstep) &
            mdstep = mdstep + 1

    end subroutine on_force

!-----------------------------------------------------------------------------

    subroutine on_mdwrit()

        use ncsu_utils
        use ncsu_umbrella
        use ncsu_sander_proxy

        implicit none

#ifndef NCSU_NO_NETCDF
        if (mode .eq. MODE_FLOODING) then
            ncsu_assert(is_master())
            call umbrella_save(umbrella, umbrella_file)
        end if
#endif /* NCSU_NO_NETCDF */

    end subroutine on_mdwrit

!><><><><><><><><><><><><><><><>< R E M D ><><><><><><><><><><><><><><><><><><

#ifdef MPI

!
! populates rem_* arrays (indexed by masterrank)
!

    subroutine rem_postinit()

        use ncsu_utils
        use ncsu_umbrella
        use ncsu_sander_proxy

        implicit none

#  include "ncsu-mpi.h"

        integer :: ng, error, i

        ncsu_assert(multisander_rem() .ne. 0)
        ncsu_assert(commmaster .ne. mpi_comm_null)

        ng = multisander_numgroup() - 1

        select case (mode)
          case (MODE_ANALYSIS)
            allocate (rem_monitor_files(0:ng), &
                rem_monitor_freqs(0:ng), stat=error)
          case (MODE_UMBRELLA)
            allocate (rem_monitor_files(0:ng), &
                rem_monitor_freqs(0:ng), rem_umbrellas(0:ng), stat=error)
          case (MODE_FLOODING)
            allocate (rem_monitor_files(0:ng), &
                rem_monitor_freqs(0:ng), rem_umbrella_files(0:ng), &
                rem_snapshots_basenames(0:ng), rem_snapshots_freqs(0:ng), &
                rem_umbrellas(0:ng), stat=error)
          case default
            ncsu_assert_not_reached()
            continue
        end select

        if (error .ne. 0) &
            NCSU_OUT_OF_MEMORY

        do i = 0, ng

            if (masterrank .eq. i) then
                rem_monitor_files(i) = monitor_file
                rem_monitor_freqs(i) = monitor_freq
            end if ! masterrank.eq.i

            call mpi_bcast(rem_monitor_files(i), SL, &
                MPI_CHARACTER, i, commmaster, error)
            ncsu_assert(error .eq. 0)

            call mpi_bcast(rem_monitor_freqs(i), 1, &
                MPI_INTEGER, i, commmaster, error)
            ncsu_assert(error .eq. 0)

            if (mode .eq. MODE_FLOODING) then

                if (masterrank .eq. i) then
                    rem_umbrella_files(i) = umbrella_file
                    rem_snapshots_basenames(i) = snapshots_basename
                    rem_snapshots_freqs(i) = snapshots_freq
                end if ! masterrank.eq.i

                call mpi_bcast(rem_umbrella_files(i), SL, &
                    MPI_CHARACTER, i, commmaster, error)
                ncsu_assert(error .eq. 0)

                call mpi_bcast(rem_snapshots_basenames(i), SL, &
                    MPI_CHARACTER, i, commmaster, error)
                ncsu_assert(error .eq. 0)

                call mpi_bcast(rem_snapshots_freqs(i), 1, &
                    MPI_INTEGER, i, commmaster, error)
                ncsu_assert(error .eq. 0)

            end if ! FLOODING

            if (mode .eq. MODE_UMBRELLA .or. mode .eq. MODE_FLOODING) then
                if (masterrank .eq. i) &
                    call umbrella_swap(rem_umbrellas(i), umbrella)

                call umbrella_bcast(rem_umbrellas(i), commmaster, i)

                if (masterrank .eq. i) &
                    call umbrella_swap(rem_umbrellas(i), umbrella)
            end if

        end do

    end subroutine rem_postinit

!-----------------------------------------------------------------------------

    subroutine rem_cleanup()

        NCSU_USE_AFAILED

        implicit none

#  include "ncsu-mpi.h"

        ncsu_assert(sanderrank .eq. 0)

        if (mode .eq. MODE_NONE) then
            continue
            ncsu_assert(.not. allocated(rem_umbrellas))
            ncsu_assert(.not. allocated(rem_monitor_files))
        else if (mode .eq. MODE_ANALYSIS) then
            ncsu_assert(allocated(rem_monitor_files))
            ncsu_assert(allocated(rem_monitor_freqs))
            deallocate (rem_monitor_files, rem_monitor_freqs)
        else if (mode .eq. MODE_UMBRELLA) then
            ncsu_assert(allocated(rem_monitor_files))
            ncsu_assert(allocated(rem_monitor_freqs))
            ncsu_assert(allocated(rem_umbrellas))
            call finalize_umbrellas()
            deallocate (rem_monitor_files, rem_monitor_freqs, rem_umbrellas)
        else if (mode .eq. MODE_FLOODING) then
            ncsu_assert(allocated(rem_monitor_files))
            ncsu_assert(allocated(rem_monitor_freqs))
            ncsu_assert(allocated(rem_umbrella_files))
            ncsu_assert(allocated(rem_snapshots_basenames))
            ncsu_assert(allocated(rem_snapshots_freqs))
            ncsu_assert(allocated(rem_umbrellas))
            call finalize_umbrellas()
            deallocate (rem_monitor_files, rem_monitor_freqs, &
                rem_snapshots_basenames, rem_snapshots_freqs, &
                rem_umbrella_files, rem_umbrellas)
        else
            continue
            ncsu_assert_not_reached()
        end if

    contains

        subroutine finalize_umbrellas

            use ncsu_umbrella, only : umbrella_fini
            use ncsu_sander_proxy, only : multisander_numgroup

            implicit none

            integer :: i

            do i = 0, multisander_numgroup() - 1
                if (masterrank .ne. i) &
                    call umbrella_fini(rem_umbrellas(i))
            end do

        end subroutine finalize_umbrellas

    end subroutine rem_cleanup

#endif /* MPI */

end module ncsu_abmd_hooks
