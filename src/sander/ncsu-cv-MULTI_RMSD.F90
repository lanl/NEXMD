! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! input:
!
! cv%i = (a1, a2, a3, 0, b1, b2, b3, b4, 0, ..., c1, c2, c3, 0)
!
!     (a[1-3] - 1st group, b[1-4] - 2nd group, etc; an atom may
!      enter a few groups simultaneously; last zero is optional;
!      empty groups [e.g., 2+ zeros in a row] are not allowed)
!
! cv%r = (a1x, a1y, a1z, a2x, a2y, a2z, a3x, a3y, a3z, b1x, ...)
!
!        (reference coordinates without '0' sentinel(s))
!
! value = sqrt[(M_1*rmsd1^2 + ... + M_N*rmsdN^2)/(M_1 + ... + M_N)]
!         M_1 - mass of group 1, ..., M_N - mass of group N
!

module ncsu_cv_MULTI_RMSD

!=============================================================================

    implicit none

    private

!=============================================================================

!
! V-table
!

    public :: colvar_value
    public :: colvar_force

    public :: colvar_bootstrap
    public :: colvar_cleanup

    public :: print_details

!=============================================================================

    type, private :: group_t

        integer :: i0, i1, r0

        NCSU_REAL, pointer :: mass(:) => null()

        NCSU_REAL, pointer :: cm_crd(:) => null()
        NCSU_REAL, pointer :: ref_crd(:) => null()

        NCSU_REAL :: total_mass
        NCSU_REAL :: ref_nrm
        NCSU_REAL :: quaternion(4)
        NCSU_REAL :: rmsd2

    end type group_t

    private :: group_bootstrap
    private :: group_finalize
    private :: group_print
    private :: group_evaluate

    type, private :: priv_t
        NCSU_REAL :: value
        NCSU_REAL :: total_mass
        type(group_t), pointer :: groups(:) => null()
#ifdef MPI
        integer :: first_cpu
#endif /* MPI */
#  include "ncsu-cv-priv.type"
    end type priv_t

#include "ncsu-cv-priv.decl"

!=============================================================================

contains

!=============================================================================

#include "ncsu-cv-priv.impl"

!=============================================================================

    subroutine group_evaluate(cv, grp, x)

        use ncsu_rmsd
        use ncsu_utils
        use ncsu_constants
        use ncsu_colvar_type

        implicit none

        type(colvar_t), intent(in)    :: cv
        type(group_t), intent(inout) :: grp

        NCSU_REAL, intent(in) :: x(*)

        NCSU_REAL :: cm(3), cm_nrm, lambda
        integer :: i, n, a, a3

        ! compute the center of mass (of the moving atoms)
        cm = ZERO
        n = 1
        do a = grp%i0, grp%i1
            a3 = 3*cv%i(a)
            cm = cm + grp%mass(n)*x(a3 - 2:a3)
            n = n + 1
        end do

        cm = cm/grp%total_mass

        ! populate cm_crd && compute the "norm"
        cm_nrm = ZERO
        n = 1
        do a = grp%i0, grp%i1
            a3 = 3*(n - 1)
            do i = 1, 3
                grp%cm_crd(a3 + i) = x(3*(cv%i(a) - 1) + i) - cm(i)
                cm_nrm = cm_nrm + grp%mass(n)*grp%cm_crd(a3 + i)**2
            end do
            n = n + 1
        end do

        call rmsd_q(grp%i1 - grp%i0 + 1, grp%mass, &
            grp%cm_crd, grp%ref_crd, lambda, grp%quaternion)

        grp%rmsd2 = max(ZERO, ((grp%ref_nrm + cm_nrm) - 2*lambda))

    end subroutine group_evaluate

!=============================================================================

    function colvar_value(cv, x) result(value)

        NCSU_USE_AFAILED

        use ncsu_constants
        use ncsu_colvar_type

        implicit none

        NCSU_REAL :: value

        type(colvar_t), intent(inout) :: cv

        NCSU_REAL, intent(in) :: x(*)

#  include "ncsu-mpi.h"

        type(priv_t), pointer :: priv
        integer :: g

#ifdef MPI
        integer :: error
#endif /* MPI */

        ncsu_assert(cv%type == COLVAR_MULTI_RMSD)

        priv => get_priv(cv)
        ncsu_assert(associated(priv%groups))
        ncsu_assert(size(priv%groups) .gt. 0)

        priv%value = ZERO

        do g = 1, size(priv%groups)
#ifdef MPI
            if (mod(g + priv%first_cpu - 1, sandersize) .eq. sanderrank) then
#endif /* MPI */
                call group_evaluate(cv, priv%groups(g), x)
                priv%value = priv%value + priv%groups(g)%rmsd2
#ifdef MPI
            end if
#endif /* MPI */
        end do

#ifdef MPI
        call mpi_reduce(priv%value, value, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, 0, commsander, error)
        ncsu_assert(error .eq. 0)
        priv%value = sqrt(value/priv%total_mass)
#else
        priv%value = sqrt(priv%value/priv%total_mass)
#endif /* MPI */

        value = priv%value

    end function colvar_value

!=============================================================================

!
! assumes that atom positions have not been
!   changed since last call to 'value()'
!

    subroutine colvar_force(cv, x, fcv, f)

        use ncsu_rmsd
        use ncsu_utils
        use ncsu_constants
        use ncsu_colvar_type

        implicit none

        type(colvar_t), intent(in) :: cv

        NCSU_REAL, intent(in) :: x(*), fcv
        NCSU_REAL, intent(inout) :: f(*)

#  include "ncsu-mpi.h"

        integer :: n, g, a, a3, f3
        NCSU_REAL :: U(3, 3), tmp
        type(priv_t), pointer :: priv

#ifdef MPI
        integer :: error
#endif /* MPI */

        ncsu_assert(cv%type == COLVAR_MULTI_RMSD)
        ncsu_assert(associated(cv%i))

        priv => get_priv(cv)

        ncsu_assert(associated(priv))
        ncsu_assert(associated(priv%groups))
        ncsu_assert(size(priv%groups) .gt. 0)

#ifdef MPI
        call mpi_bcast(priv%value, 1, MPI_DOUBLE_PRECISION, 0, commsander, error)
        ncsu_assert(error .eq. 0)
#endif /* MPI */

        do g = 1, size(priv%groups)

#ifdef MPI
            if (mod(g + priv%first_cpu - 1, sandersize) .ne. sanderrank) &
                cycle
#endif /* MPI */

            call rmsd_q2u(priv%groups(g)%quaternion, U)
            tmp = fcv/max(priv%value, NCSU_TO_REAL(0.000001))/priv%total_mass

            n = 1
            do a = priv%groups(g)%i0, priv%groups(g)%i1
                a3 = 3*n
                f3 = 3*cv%i(a)

                f(f3 - 2:f3) = f(f3 - 2:f3) &
                    + tmp*priv%groups(g)%mass(n)*(priv%groups(g)%cm_crd(a3 - 2:a3) &
                    - matmul(U, priv%groups(g)%ref_crd(a3 - 2:a3)))
                n = n + 1
            end do
        end do

#ifdef MPI
        priv%first_cpu = mod(priv%first_cpu + 1, sandersize)
#endif /* MPI */

    end subroutine colvar_force

!=============================================================================

    subroutine group_bootstrap(grp, cv, cvno, amass, i0, i1, r0)

        use ncsu_utils
        use ncsu_constants
        use ncsu_colvar_type
        use ncsu_sander_proxy

        implicit none

        type(group_t), intent(inout) :: grp

        type(colvar_t), intent(in) :: cv
        integer, intent(in) :: cvno
        NCSU_REAL, intent(in) :: amass(*)

        integer, intent(in) :: i0, i1, r0

#  include "ncsu-mpi.h"

        integer :: a, b, n_atoms, error
        NCSU_REAL :: cm(3)

        ncsu_assert(cv%type == COLVAR_MULTI_RMSD)
        ncsu_assert(associated(cv%i) .and. associated(cv%r))
        ncsu_assert(i0 .gt. 0 .and. i0 .le. size(cv%i))
        ncsu_assert(i1 .gt. 0 .and. i1 .le. size(cv%i))
        ncsu_assert(r0 .gt. 0 .and. r0 .lt. size(cv%r))

        ! basic checks
        if (i0 + 2 > i1) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, &
                fmt='(/a,a,'//pfmt(cvno)//',a,'//pfmt &
                (i0)//',a,'//pfmt(i1)//',a/)') &
                NCSU_ERROR, 'CV #', cvno, &
                ' (MULTI_RMSD) : too few integers in group (', i0, ':', i1, ')'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if

        do a = i0, i1
            ncsu_assert(a .le. size(cv%i))
            if (cv%i(a) < 1 .or. cv%i(a) > sander_natoms()) then
                NCSU_MASTER_ONLY_BEGIN
                write (unit=ERR_UNIT, &
                    fmt='(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                    (cv%i(a))//',a,'//pfmt(sander_natoms())//',a/)') &
                    NCSU_ERROR, 'CV #', cvno, &
                    ' (MULTI_RMSD) : integer #', a, ' (', cv%i(a), &
                    ') is out of range [1, ', sander_natoms(), ']'
                NCSU_MASTER_ONLY_END
                call terminate()
            end if
        end do

        ! check for duplicates
        do a = i0, i1
            do b = a + 1, i1
                if (cv%i(a) == cv%i(b)) then
                    NCSU_MASTER_ONLY_BEGIN
                    write (unit=ERR_UNIT, &
                        fmt='(/a,a,'//pfmt(cvno)//',a,'//pfmt(a)//',a,'//pfmt &
                        (b)//',a,'//pfmt(cv%i(a))//',a/)') &
                        NCSU_ERROR, 'CV #', cvno, &
                        ' (MULTI_RMSD) : integers #', a, ' and #', b, &
                        ' are equal (', cv%i(a), ')'
                    NCSU_MASTER_ONLY_END
                    call terminate()
                end if ! cv%i(a) == cv%i(b)
            end do
        end do

        ! allocate/setup

        n_atoms = i1 - i0 + 1

        grp%i0 = i0
        grp%i1 = i1
        grp%r0 = r0

        allocate (grp%mass(n_atoms), grp%cm_crd(3*n_atoms), &
            grp%ref_crd(3*n_atoms), stat=error)

        if (error /= 0) &
            NCSU_OUT_OF_MEMORY

        cm = ZERO
        grp%total_mass = ZERO

        do a = 1, n_atoms
            grp%mass(a) = amass(cv%i(a + grp%i0 - 1))
            grp%total_mass = grp%total_mass + grp%mass(a)
            cm = cm + grp%mass(a)*cv%r(r0 + 3*(a - 1):r0 + 3*a - 1)
        end do

        cm = cm/grp%total_mass

        ! translate reference coordinates to CM frame
        grp%ref_nrm = ZERO
        do a = 1, n_atoms
            do b = 1, 3
                grp%ref_crd(3*(a - 1) + b) = cv%r(r0 + 3*(a - 1) + b - 1) - cm(b)
                grp%ref_nrm = grp%ref_nrm + grp%mass(a)*grp%ref_crd(3*(a - 1) + b)**2
            end do
        end do

    end subroutine group_bootstrap

!=============================================================================

    subroutine group_finalize(grp)

        NCSU_USE_AFAILED

        implicit none

        type(group_t), intent(inout) :: grp

        ncsu_assert(associated(grp%mass))
        ncsu_assert(associated(grp%cm_crd))
        ncsu_assert(associated(grp%ref_crd))

        deallocate (grp%mass, grp%cm_crd, grp%ref_crd)

    end subroutine group_finalize

!=============================================================================

    subroutine group_print(cv, grp, lun)

        use ncsu_utils
        use ncsu_colvar_type
        use ncsu_sander_proxy

        implicit none

        type(colvar_t), intent(in) :: cv
        type(group_t), intent(in) :: grp
        integer, intent(in) :: lun

        integer :: a, c
        character(4) :: aname

        ncsu_assert(associated(cv%i) .and. associated(cv%r))
        ncsu_assert(grp%i0 .gt. 0 .and. grp%i0 .le. size(cv%i))
        ncsu_assert(grp%i1 .gt. 0 .and. grp%i1 .le. size(cv%i))

        write (unit=lun, fmt='(a,a)', advance='NO') NCSU_INFO, 'atoms = ('

        c = 1
        do a = grp%i0, grp%i1

            ncsu_assert(cv%i(a) > 0 .and. cv%i(a) <= sander_natoms())
            aname = sander_atom_name(cv%i(a))

            write (unit=lun, fmt='('//pfmt(cv%i(a))//',a,a,a)', advance='NO') &
                cv%i(a), ' [', trim(aname), ']'

            if (a == grp%i1) then
                write (unit=lun, fmt='(a)') ')'
            else if (mod(c, 5) == 0) then
                write (unit=lun, fmt='(a,/a,3x)', advance='NO') ',', NCSU_INFO
            else
                write (unit=lun, fmt='(a)', advance='NO') ', '
            end if

            c = c + 1
        end do

        write (unit=lun, fmt='(a,a)') &
            NCSU_INFO, 'reference coordinates :'

        c = 0
        do a = grp%i0, grp%i1
            ncsu_assert(grp%r0 + c + 2 <= size(cv%r))
            write (unit=lun, fmt='(a,3x,i5,a,f8.3,a,f8.3,a,f8.3)') &
                NCSU_INFO, cv%i(a), ' : ', &
                cv%r(grp%r0 + c), ', ', cv%r(grp%r0 + c + 1), &
                ', ', cv%r(grp%r0 + c + 2)
            c = c + 3
        end do

    end subroutine group_print

!=============================================================================

    subroutine colvar_bootstrap(cv, cvno, amass)

        use ncsu_utils
        use ncsu_constants
        use ncsu_colvar_type
        use ncsu_sander_proxy

        implicit none

        type(colvar_t), intent(inout) :: cv
        integer, intent(in)    :: cvno
        NCSU_REAL, intent(in)    :: amass(*)

        integer :: i, i0, n_atoms, n_groups, error

        type(priv_t), pointer :: priv

#  include "ncsu-mpi.h"

        ncsu_assert(cv%type == COLVAR_MULTI_RMSD)

        ! very basic checks
        if (.not. associated(cv%i)) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                NCSU_ERROR, 'CV #', cvno, &
                ' (MULTI_RMSD) : no integers found'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if ! .not. associated(cv%i)

        if (.not. associated(cv%r)) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                NCSU_INFO, 'CV #', cvno, &
                ' (MULTI_RMSD) : no reals (reference coordinates) found'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if

        ! count the groups (number of zeros in the cv%i array)
        n_atoms = 0
        n_groups = 0
        i0 = 1

        do i = 1, size(cv%i)
            if (cv%i(i) == 0) then
                n_groups = n_groups + 1
                if (i .eq. i0) then
                    NCSU_MASTER_ONLY_BEGIN
                    write (unit=ERR_UNIT, &
                        fmt='(/a,a,'//pfmt(cvno)//',a,'//pfmt(i)//',a/)') &
                        NCSU_ERROR, 'CV #', cvno, &
                        ' (MULTI_RMSD) : unexpected zero (integer #', i, ')'
                    NCSU_MASTER_ONLY_END
                    call terminate()
                end if ! i .eq. i0
                i0 = i + 1
            else if (cv%i(i) .lt. 0) then
                NCSU_MASTER_ONLY_BEGIN
                write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                    NCSU_ERROR, 'CV #', cvno, &
                    ' (MULTI_RMSD) : negative integer'
                NCSU_MASTER_ONLY_END
                call terminate()
            else
                n_atoms = n_atoms + 1
            end if
        end do

        if (size(cv%i) .gt. 0 .and. cv%i(size(cv%i)) .gt. 0) &
            n_groups = n_groups + 1

        ncsu_assert(n_groups .gt. 0)

        if (size(cv%r) /= 3*n_atoms) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                NCSU_ERROR, 'CV #', cvno, &
                ' (MULTI_RMSD) : wrong number of reals'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if ! size(cv%r) /= 3*n_atoms

        ! allocate priv_t instance for this variable
        priv => new_priv(cv)

        ! allocate/setup groups
        allocate (priv%groups(n_groups), stat=error)
        if (error /= 0) &
            NCSU_OUT_OF_MEMORY

        i0 = 1 ! first atom
        n_atoms = 1
        n_groups = 0

        priv%total_mass = ZERO

        do i = 1, size(cv%i)
            if (cv%i(i) .eq. 0) then
                n_groups = n_groups + 1
                ncsu_assert(n_groups .le. size(priv%groups))
                call group_bootstrap(priv%groups(n_groups), &
                    cv, cvno, amass, i0, i - 1, n_atoms)
                n_atoms = n_atoms + 3*(i - i0)
                i0 = i + 1
                priv%total_mass = priv%total_mass &
                    + priv%groups(n_groups)%total_mass
            end if
        end do

        if (size(cv%i) .gt. 0 .and. cv%i(size(cv%i)) .gt. 0) then
            n_groups = n_groups + 1
            ncsu_assert(n_groups .le. size(priv%groups))
            call group_bootstrap(priv%groups(n_groups), &
                cv, cvno, amass, i0, size(cv%i), n_atoms)
            priv%total_mass = priv%total_mass &
                + priv%groups(n_groups)%total_mass
        end if

#ifdef MPI
        priv%first_cpu = 0
#endif /* MPI */

    end subroutine colvar_bootstrap

!=============================================================================

    subroutine colvar_cleanup(cv)

        NCSU_USE_AFAILED

        use ncsu_colvar_type

        implicit none

        type(colvar_t), intent(inout) :: cv

        integer :: g
        type(priv_t), pointer :: priv

        ncsu_assert(cv%type == COLVAR_MULTI_RMSD)
        ncsu_assert(associated(cv%i) .and. size(cv%i) .gt. 0)

        priv => get_priv(cv)
        ncsu_assert(associated(priv%groups))

        do g = 1, size(priv%groups)
            call group_finalize(priv%groups(g))
        end do

        deallocate (priv%groups)
        call del_priv(cv)

    end subroutine colvar_cleanup

!=============================================================================

    subroutine print_details(cv, lun)

        use ncsu_utils
        use ncsu_colvar_type
        use ncsu_colvar_utils
        use ncsu_sander_proxy

        implicit none

        type(colvar_t), intent(in) :: cv
        integer, intent(in) :: lun

        integer :: g
        type(priv_t), pointer :: priv

        ncsu_assert(is_master())
        ncsu_assert(cv%type == COLVAR_MULTI_RMSD)
        ncsu_assert(associated(cv%i))
        ncsu_assert(associated(cv%r))

        priv => get_priv(cv)
        ncsu_assert(associated(priv%groups))

        do g = 1, size(priv%groups)
            write (unit=lun, fmt='(a,a,'//pfmt(g)//',a)') &
                NCSU_INFO, '<> group <> #', g, ':'
            call group_print(cv, priv%groups(g), lun)
        end do

    end subroutine print_details

!=============================================================================

end module ncsu_cv_MULTI_RMSD
