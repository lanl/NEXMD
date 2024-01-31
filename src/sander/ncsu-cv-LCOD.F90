! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! LCOD (Linear Combination Of Distances)
!
! cv%i = (a11, a12, a21, a22, ..., aN1, aN2)
!
!     indexes of the participating atoms
!
! cv%r = (r1, r2, ..., rN)
!
!     non-zero weights
!
! value = r1*d1 + r2*d2 + ... + rN*dN
!
!     d1 -- distance between atoms #a11 and #a12
!     d2 -- distance between atoms #a21 and #a22
!                       . . .
!     dN -- distance between atoms #aN1 and #aN2
!

module ncsu_cv_LCOD

!=============================================================================

    implicit none

    private

!=============================================================================

    public :: colvar_value
    public :: colvar_force

    public :: colvar_bootstrap
    public :: print_details

!=============================================================================

contains

!=============================================================================

    function colvar_value(cv, x) result(value)

#  ifndef NCSU_DISABLE_ASSERT
        use ncsu_utils
        use ncsu_sander_proxy
#  endif /* NCSU_DISABLE_ASSERT */

        use ncsu_constants
        use ncsu_colvar_type
        use ncsu_colvar_math

        implicit none

        NCSU_REAL :: value

        type(colvar_t), intent(in) :: cv

        NCSU_REAL, intent(in) :: x(*)

#  include "ncsu-mpi.h"

        integer :: a1, a2, n

        ncsu_assert(cv%type .eq. COLVAR_LCOD)

        ncsu_assert(associated(cv%i))
        ncsu_assert(associated(cv%r))

        ncsu_assert(size(cv%i) .gt. 0)
        ncsu_assert(size(cv%r) .gt. 0)

        ncsu_assert(mod(size(cv%i), 2) .eq. 0)
        ncsu_assert((size(cv%i)/2) .eq. size(cv%r))

        value = ZERO

        NCSU_MASTER_ONLY_BEGIN

        do n = 1, size(cv%r)

            ncsu_assert(cv%i(2*n - 1) .gt. 0 .and. cv%i(2*n - 1) .le. sander_natoms())
            a1 = 3*cv%i(2*n - 1) - 2

            ncsu_assert(cv%i(2*n) .gt. 0 .and. cv%i(2*n) .le. sander_natoms())
            a2 = 3*cv%i(2*n) - 2

            ncsu_assert(a1 .ne. a2)

            value = value + cv%r(n)*distance(x(a1:a1 + 2), x(a2:a2 + 2))

        end do

        NCSU_MASTER_ONLY_END

    end function colvar_value

!=============================================================================

    subroutine colvar_force(cv, x, fcv, f)

#  ifndef NCSU_DISABLE_ASSERT
        use ncsu_utils
        use ncsu_sander_proxy
#  endif /* NCSU_DISABLE_ASSERT */

        use ncsu_colvar_type
        use ncsu_colvar_math

        implicit none

        type(colvar_t), intent(in) :: cv

        NCSU_REAL, intent(in) :: x(*), fcv

        NCSU_REAL, intent(inout) :: f(*)

#  include "ncsu-mpi.h"

        NCSU_REAL :: d1(3), d2(3)

        integer :: a1, a2, n

        ncsu_assert(cv%type .eq. COLVAR_LCOD)

        ncsu_assert(associated(cv%i))
        ncsu_assert(associated(cv%r))

        ncsu_assert(size(cv%i) .gt. 0)
        ncsu_assert(size(cv%r) .gt. 0)

        ncsu_assert(mod(size(cv%i), 2) .eq. 0)
        ncsu_assert((size(cv%i)/2) .eq. size(cv%r))

        NCSU_MASTER_ONLY_BEGIN

        do n = 1, size(cv%r)

            ncsu_assert(cv%i(2*n - 1) .gt. 0 .and. cv%i(2*n - 1) .le. sander_natoms())
            a1 = 3*cv%i(2*n - 1) - 2

            ncsu_assert(cv%i(2*n) .gt. 0 .and. cv%i(2*n) .le. sander_natoms())
            a2 = 3*cv%i(2*n) - 2

            ncsu_assert(a1 .ne. a2)

            call distance_d(x(a1:a1 + 2), x(a2:a2 + 2), d1, d2)

            f(a1:a1 + 2) = f(a1:a1 + 2) + cv%r(n)*fcv*d1
            f(a2:a2 + 2) = f(a2:a2 + 2) + cv%r(n)*fcv*d2

        end do

        NCSU_MASTER_ONLY_END

    end subroutine colvar_force

!=============================================================================

    subroutine colvar_bootstrap(cv, cvno, amass)

        use ncsu_utils
        use ncsu_constants
        use ncsu_colvar_type
        use ncsu_colvar_utils
        use ncsu_sander_proxy

        implicit none

        type(colvar_t), intent(inout) :: cv
        integer, intent(in)    :: cvno
        NCSU_REAL, intent(in)    :: amass(*)

#  include "ncsu-mpi.h"

        integer :: n, j

        ncsu_assert(cv%type .eq. COLVAR_LCOD)

        if (.not. associated(cv%i)) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                NCSU_ERROR, 'CV #', cvno, ' (LCOD) : no integers found'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if ! .not.associated(cv%i)

        if (.not. associated(cv%r)) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                NCSU_ERROR, 'CV #', cvno, ' (LCOD) : no reals found'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if ! .not.associated(cv%r)

        if (mod(size(cv%i), 2) .ne. 0) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                NCSU_ERROR, 'CV #', cvno, ' (LCOD) : odd number of integers'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if ! mod(size(cv%i), 2).ne.0

        if (size(cv%i) .lt. 2) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                NCSU_ERROR, 'CV #', cvno, ' (LCOD) : too few integers'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if ! size(cv%i).lt.2

        if (size(cv%r) .ne. size(cv%i)/2) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                NCSU_ERROR, 'CV #', cvno, ' (LCOD) : size(cv%r).ne.size(cv%i)/2'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if ! size(cv%r).ne.size(cv%i)/2

        do n = 1, size(cv%r)
            if (abs(cv%r(n)) .lt. 1.0D-8) then
                NCSU_MASTER_ONLY_BEGIN
                write (unit=ERR_UNIT, &
                    fmt='(/a,a,'//pfmt(cvno)//',a,'//pfmt(n)//',a/)') &
                    NCSU_ERROR, 'CV #', cvno, ' (LCOD) : real number #', n, &
                    ' is too small'
                NCSU_MASTER_ONLY_END
                call terminate()
            end if ! abs(cv%r(n)).lt.1.0D-8

            do j = 0, 1
                if (cv%i(2*n - j) .lt. 1 .or. cv%i(2*n - j) .gt. sander_natoms()) then
                    NCSU_MASTER_ONLY_BEGIN
                    write (unit=ERR_UNIT, &
                        fmt='(/a,a,'//pfmt(cvno)//',a,'//pfmt(2*n - j)//',a,'//pfmt &
                        (cv%i(2*n - j))//',a,'//pfmt(sander_natoms())//',a/)') &
                        NCSU_ERROR, 'CV #', cvno, ' (LCOD) : integer #', &
                        (2*n - j), ' (', cv%i(2*n - j), ') is out of range [1, ', &
                        sander_natoms(), ']'
                    NCSU_MASTER_ONLY_END
                    call terminate()
                end if
            end do

            if (cv%i(2*n - 1) .eq. cv%i(2*n)) then
                NCSU_MASTER_ONLY_BEGIN
                write (unit=ERR_UNIT, &
                    fmt='(/a,a,'//pfmt(cvno)//',a,'//pfmt(2*n - 1)//',a,'//pfmt &
                    (2*n)//',a,'//pfmt(cv%i(2*n))//',a/)') NCSU_ERROR, 'CV #', &
                    cvno, ' (LCOD) : integers #', (2*n - 1), ' and #', (2*n), &
                    ' are equal (', cv%i(2*n), ')'
                NCSU_MASTER_ONLY_END
                call terminate()
            end if
        end do

    end subroutine colvar_bootstrap

!=============================================================================

    subroutine print_details(cv, lun)

        use ncsu_utils
        use ncsu_colvar_type
        use ncsu_colvar_utils
        use ncsu_sander_proxy

        implicit none

        type(colvar_t), intent(in) :: cv
        integer, intent(in) :: lun

        integer :: n, a1, a2
        character(4) :: aname1, aname2

        ncsu_assert(is_master())
        ncsu_assert(cv%type .eq. COLVAR_LCOD)

        ncsu_assert(associated(cv%i))
        ncsu_assert(associated(cv%r))

        ncsu_assert(size(cv%i) .gt. 0)
        ncsu_assert(size(cv%r) .gt. 0)

        ncsu_assert(mod(size(cv%i), 2) .eq. 0)
        ncsu_assert((size(cv%i)/2) .eq. size(cv%r))

        do n = 1, size(cv%r)

            a1 = cv%i(2*n - 1)
            a2 = cv%i(2*n)

            aname1 = sander_atom_name(a1)
            aname2 = sander_atom_name(a2)

            write (unit=lun, fmt='(a,4x,f8.3,a,'//pfmt &
                (a1)//',a,a,a,'//pfmt(a2)//',a,a,a)') NCSU_INFO, cv%r(n), ' * (', &
                a1, ' [', trim(aname1), '] <=> ', a2, ' [', trim(aname2), '])'

        end do

    end subroutine print_details

!=============================================================================

end module ncsu_cv_LCOD
