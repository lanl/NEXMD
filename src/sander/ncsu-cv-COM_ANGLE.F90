! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

! vbabin-at-ncsu-dot-edu, 09/08/2010
!
! COM_ANGLE -- angle formed by the centers of mass of three atom groups
!
! input: i = (a1, ..., aN, 0, b1, ..., bM, 0, c1, ..., cK, 0)
! last zero is optional; a?/b?/c? -- indices of the participating atoms

module ncsu_cv_COM_ANGLE

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

        NCSU_USE_AFAILED

        use ncsu_constants, only : ZERO
        use ncsu_colvar_type
        use ncsu_colvar_math
        use ncsu_colvar_utils

        implicit none

        NCSU_REAL :: value

        type(colvar_t), intent(in) :: cv
        NCSU_REAL, intent(in) :: x(*)

#  include "ncsu-mpi.h"

        integer :: n

        NCSU_REAL :: cm1(3), cm2(3), cm3(3)

        ncsu_assert(cv%type == COLVAR_COM_ANGLE)

        ncsu_assert(associated(cv%i))
        ncsu_assert(size(cv%i) .ge. 5)
        ncsu_assert(associated(cv%r))
        ncsu_assert(size(cv%r) .eq. size(cv%i))

#ifdef MPI
        if (sanderrank .eq. 0) then
#endif /* MPI */
            n = 1
            call group_com(cv, x, n, cm1)
            ncsu_assert(n .le. size(cv%i))
            call group_com(cv, x, n, cm2)
            ncsu_assert(n .le. size(cv%i))
            call group_com(cv, x, n, cm3)

            value = angle(cm1, cm2, cm3)
#ifdef MPI
        else
            value = ZERO
        end if
#endif /* MPI */

    end function colvar_value

!=============================================================================

    subroutine colvar_force(cv, x, fcv, f)

        NCSU_USE_AFAILED

        use ncsu_constants, only : ERR_UNIT
        use ncsu_colvar_type
        use ncsu_colvar_math
        use ncsu_colvar_utils
        use ncsu_sander_proxy

        implicit none

        type(colvar_t), intent(in) :: cv
        NCSU_REAL, intent(in) :: x(*), fcv
        NCSU_REAL, intent(inout) :: f(*)

#  include "ncsu-mpi.h"

        NCSU_REAL :: d1(3), d2(3), d3(3)
        NCSU_REAL :: cm1(3), cm2(3), cm3(3)

        integer :: n

        ncsu_assert(cv%type == COLVAR_COM_ANGLE)

        ncsu_assert(associated(cv%i))
        ncsu_assert(size(cv%i) .ge. 5)
        ncsu_assert(associated(cv%r))
        ncsu_assert(size(cv%r) .eq. size(cv%i))

        NCSU_MASTER_ONLY_BEGIN

        n = 1
        call group_com(cv, x, n, cm1)
        call group_com(cv, x, n, cm2)
        call group_com(cv, x, n, cm3)

        if (distance(cm1, cm2) .lt. 1.0d-8 .or. distance(cm2, cm3) .lt. 1.0d-8) then
            write (unit=ERR_UNIT, fmt='(/a,a/)') NCSU_ERROR, &
                'COM_ANGLE : centers of mass got too close to each other'
            call terminate()
        end if ! distance too small

        call angle_d(cm1, cm2, cm3, d1, d2, d3)

        d1 = fcv*d1
        d2 = fcv*d2
        d3 = fcv*d3

        n = 1
        call group_com_d(cv, f, d1, n)
        call group_com_d(cv, f, d2, n)
        call group_com_d(cv, f, d3, n)

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

        integer :: error, ngroups

        ncsu_assert(cv%type == COLVAR_COM_ANGLE)

        call com_check_i(cv%i, cvno, 'COM_ANGLE', ngroups)
        if (ngroups .ne. 3) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(/a,a,'//pfmt(cvno)//',a/)') &
                NCSU_ERROR, 'CV #', cvno, &
                ' (COM_ANGLE) : number of atom groups is not three'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if ! ngroups.ne.3

        if (associated(cv%r)) &
            deallocate (cv%r)

        allocate (cv%r(size(cv%i)), stat=error)
        if (error .ne. 0) &
            NCSU_OUT_OF_MEMORY

        cv%r = ZERO
        do error = 1, size(cv%i)
            if (cv%i(error) .gt. 0) &
                cv%r(error) = amass(cv%i(error))
        end do

        print *, cv%r
        call com_init_weights(cv%r)
        print *, cv%r

    end subroutine colvar_bootstrap

!=============================================================================

    subroutine print_details(cv, lun)

#ifndef NCSU_DISABLE_ASSERT
        use ncsu_utils
        use ncsu_sander_proxy
#endif /* NCSU_DISABLE_ASSERT */

        use ncsu_colvar_type
        use ncsu_colvar_utils

        implicit none

        type(colvar_t), intent(in) :: cv
        integer, intent(in) :: lun

        ncsu_assert(is_master())
        ncsu_assert(cv%type == COLVAR_COM_ANGLE)
        ncsu_assert(associated(cv%i))

        call com_print_i(cv%i, lun)

    end subroutine print_details

!=============================================================================

end module ncsu_cv_COM_ANGLE
