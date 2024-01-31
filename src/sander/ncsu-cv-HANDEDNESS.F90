! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

!
! cv%i = (i1, ..., iN) -- list of participating atoms
!

module ncsu_cv_HANDEDNESS

!=============================================================================

    implicit none

    private

!=============================================================================

    public :: colvar_value
    public :: colvar_force

    public :: colvar_bootstrap
    public :: print_details

    NCSU_REAL, parameter, private :: TINY = 0.00000100000000000000D0 ! NCSU_TO_REAL(0.000001)
    private :: value4, derivative4

!=============================================================================

contains

!=============================================================================

    function colvar_value(cv, x) result(value)

        use ncsu_utils
        use ncsu_constants
        use ncsu_colvar_type

        implicit none

        NCSU_REAL :: value

        type(colvar_t), intent(inout) :: cv

        NCSU_REAL, intent(in) :: x(*)

        integer :: natoms, a, j1, j2, j3, j4

        ncsu_assert(cv%type .eq. COLVAR_HANDEDNESS)

        ncsu_assert(associated(cv%i))
        ncsu_assert(associated(cv%r))

        natoms = size(cv%i)
        ncsu_assert(natoms .gt. 3)

        value = ZERO

        do a = 1, natoms - 3

            j1 = 3*cv%i(a + 0)
            j2 = 3*cv%i(a + 1)
            j3 = 3*cv%i(a + 2)
            j4 = 3*cv%i(a + 3)

            value = value + value4(cv%r(1), &
                x(j1 - 2:j1), x(j2 - 2:j2), x(j3 - 2:j3), x(j4 - 2:j4))

        end do

    end function colvar_value

!=============================================================================

    subroutine colvar_force(cv, x, fcv, f)

        use ncsu_utils
        use ncsu_colvar_type

        implicit none

        type(colvar_t), intent(in) :: cv

        NCSU_REAL, intent(in) :: x(*), fcv
        NCSU_REAL, intent(inout) :: f(*)

#  include "ncsu-mpi.h"

        integer :: natoms, a, j1, j2, j3, j4
        NCSU_REAL :: d1(3), d2(3), d3(3), d4(3)

        ncsu_assert(cv%type .eq. COLVAR_HANDEDNESS)

        ncsu_assert(associated(cv%i))
        ncsu_assert(associated(cv%r))

        natoms = size(cv%i)
        ncsu_assert(natoms .gt. 3)

        NCSU_MASTER_ONLY_BEGIN

        do a = 1, natoms - 3

            j1 = 3*cv%i(a + 0)
            j2 = 3*cv%i(a + 1)
            j3 = 3*cv%i(a + 2)
            j4 = 3*cv%i(a + 3)

            call derivative4(cv%r(1), &
                x(j1 - 2:j1), x(j2 - 2:j2), x(j3 - 2:j3), x(j4 - 2:j4), &
                d1, d2, d3, d4)

            f(j1 - 2:j1) = f(j1 - 2:j1) + fcv*d1
            f(j2 - 2:j2) = f(j2 - 2:j2) + fcv*d2
            f(j3 - 2:j3) = f(j3 - 2:j3) + fcv*d3
            f(j4 - 2:j4) = f(j4 - 2:j4) + fcv*d4

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

        integer :: natoms, error

#  include "ncsu-mpi.h"

        ncsu_assert(cv%type == COLVAR_HANDEDNESS)

        natoms = size(cv%i)

        call check_i(cv%i, cvno, 'HANDEDNESS')
        if (.not. natoms .gt. 3) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(a,a,'//pfmt(cvno)//',a)') &
                NCSU_ERROR, 'CV #', cvno, &
                ' (HANDEDNESS) : too few integers'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if

        if (.not. associated(cv%r)) then
            allocate (cv%r(1), stat=error)
            if (error .ne. 0) &
                NCSU_OUT_OF_MEMORY
            cv%r(1) = ZERO
        end if

        if (size(cv%r) .ne. 1) then
            NCSU_MASTER_ONLY_BEGIN
            write (unit=ERR_UNIT, fmt='(a,a,'//pfmt(cvno)//',a)') &
                NCSU_ERROR, 'CV #', cvno, &
                ' (HANDEDNESS) : unexpected number of reals (just 1 is needed)'
            NCSU_MASTER_ONLY_END
            call terminate()
        end if

        cv%r(1) = max(ZERO, cv%r(1))
        cv%r(1) = min(ONE, cv%r(1))

    end subroutine colvar_bootstrap

!=============================================================================

    subroutine print_details(cv, lun)

        use ncsu_utils
        use ncsu_colvar_type
        use ncsu_colvar_utils

        implicit none

        type(colvar_t), intent(in) :: cv
        integer, intent(in) :: lun

        ncsu_assert(cv%type .eq. COLVAR_HANDEDNESS)

        ncsu_assert(associated(cv%i))
        ncsu_assert(associated(cv%r))

        write (unit=lun, fmt='(a,a,'//pfmt(cv%r(1), 3)//')') &
            NCSU_INFO, '      w = ', cv%r(1)
        call print_i(cv%i, lun)

    end subroutine print_details

!=============================================================================

!
! ** value4() & derivative4() are maple-generated & postprocessed by hands **
!
! dot3 := proc(u, v)
!     u[1]*v[1] + u[2]*v[2] + u[3]*v[3]
! end proc;
!
! cross3 := proc(u, v)
!     [
!         u[2]*v[3] - u[3]*v[2],
!         u[3]*v[1] - u[1]*v[3],
!         u[1]*v[2] - v[1]*u[2]
!     ]
! end proc;
!
!
! r1 := array(1..3, [r1x, r1y, r1z]);
! r2 := array(1..3, [r2x, r2y, r2z]);
! r3 := array(1..3, [r3x, r3y, r3z]);
! r4 := array(1..3, [r4x, r4y, r4z]);
!
! u1 := evalm(r2 - r1);
! u2 := evalm(r4 - r3);
! u3 := evalm((1 - w)*(r3 - r2) + w*(r4 - r1));
!
! u1n := dot3(u1, u1) + TINY;
! u2n := dot3(u2, u2) + TINY;
! u3n := dot3(u3, u3) + TINY;
!
! dofs := [r1x, r1y, r1z, r2x, r2y, r2z, r3x, r3y, r3z, r4x, r4y, r4z];
!
! handedness := dot3(u3, cross3(u1, u2))/sqrt(u1n*u2n*u3n);
! handedness := simplify(handedness);
!
! handedness_v := codegen[makeproc](handedness, dofs);
! handedness_v := codegen[optimize](handedness_v, tryhard);
!
! handedness_d := codegen[GRADIENT]
!    (codegen[split](handedness_v, dofs), dofs,
!     function_value = false, mode = reverse);
!
! handedness_v := codegen[packargs](handedness_v, [r1x, r1y, r1z], `x1`);
! handedness_v := codegen[packargs](handedness_v, [r2x, r2y, r2z], `x2`);
! handedness_v := codegen[packargs](handedness_v, [r3x, r3y, r3z], `x3`);
! handedness_v := codegen[packargs](handedness_v, [r4x, r4y, r4z], `x4`);
!
! handedness_d := codegen[optimize](handedness_d, tryhard);
! handedness_d := codegen[packargs](handedness_d, [r1x, r1y, r1z], `x1`);
! handedness_d := codegen[packargs](handedness_d, [r2x, r2y, r2z], `x2`);
! handedness_d := codegen[packargs](handedness_d, [r3x, r3y, r3z], `x3`);
! handedness_d := codegen[packargs](handedness_d, [r4x, r4y, r4z], `x4`);
!
! CodeGeneration[Fortran](handedness_v, deducetypes = false,
!                         optimize, defaulttype = numeric,
!                                                 output = `ncsu-handedness.maple.f`);
!
! CodeGeneration[Fortran](handedness_d, deducetypes = false,
!                         optimize, defaulttype = numeric,
!                                                 output = `ncsu-handedness.maple.f`);
!

    NCSU_REAL function value4(w, x1, x2, x3, x4)

    implicit none

    NCSU_REAL, intent(in) :: w
    NCSU_REAL, intent(in) :: x1(3)
    NCSU_REAL, intent(in) :: x2(3)
    NCSU_REAL, intent(in) :: x3(3)
    NCSU_REAL, intent(in) :: x4(3)

    NCSU_REAL ::  t13
    NCSU_REAL ::  t14
    NCSU_REAL ::  t119
    NCSU_REAL ::  t1
    NCSU_REAL ::  t23
    NCSU_REAL ::  t22
    NCSU_REAL ::  t34
    NCSU_REAL ::  t5
    NCSU_REAL ::  t6
    NCSU_REAL ::  t12
    NCSU_REAL ::  t40
    NCSU_REAL ::  t18
    NCSU_REAL ::  t7
    NCSU_REAL ::  t55
    NCSU_REAL ::  t32
    NCSU_REAL ::  t8
    NCSU_REAL ::  t25
    NCSU_REAL ::  t30
    NCSU_REAL ::  t29
    NCSU_REAL ::  t43
    NCSU_REAL ::  t35
    NCSU_REAL ::  t21
    NCSU_REAL ::  t2
    NCSU_REAL ::  t9
    NCSU_REAL ::  t27
    NCSU_REAL ::  t36
    NCSU_REAL ::  t24
    NCSU_REAL ::  t28
    NCSU_REAL ::  t26
    NCSU_REAL ::  t3
    NCSU_REAL ::  t15
    NCSU_REAL ::  t31
    NCSU_REAL ::  t33
    NCSU_REAL ::  t17
    NCSU_REAL ::  t16
    NCSU_REAL ::  t39
    NCSU_REAL ::  t52
    NCSU_REAL ::  t41
    NCSU_REAL ::  t42
    NCSU_REAL ::  t10
    NCSU_REAL ::  t4
    NCSU_REAL ::  t11
    NCSU_REAL ::  t20
    NCSU_REAL ::  t19
    NCSU_REAL ::  t37
    NCSU_REAL ::  t38

    t1 = x2(3)
    t7 = t1**2
    t2 = x2(2)
    t9 = t2**2
    t3 = x2(1)
    t11 = t3**2
    t32 = t7 + t9 + t11
    t4 = x3(3)
    t8 = t4**2
    t5 = x3(2)
    t10 = t5**2
    t6 = x3(1)
    t12 = t6**2
    t31 = t8 + t10 + t12
    t13 = 0.2D1*t2
    t30 = -t13
    t33 = 0.2D1*t3
    t29 = -t33
    t34 = 0.2D1*t4
    t28 = -t34
    t35 = x4(1)
    t36 = x1(1)
    t27 = t35 - t36
    t37 = x4(2)
    t38 = x1(2)
    t26 = t37 - t38
    t39 = x4(3)
    t40 = x1(3)
    t25 = t39 - t40
    t24 = -t35 + t6
    t23 = -t37 + t5
    t22 = t4 - t39
    t21 = t33
    t20 = t13
    t19 = 0.2D1*t1
    t41 = t40**2
    t42 = t38**2
    t43 = t36**2
    t18 = t41 + t42 + t43 + t32
    t17 = -t6*t37 + t5*t35
    t16 = -t5*t39 + t4*t37
    t15 = t6*t39 - t4*t35
    t52 = 0.2D1*t5
    t55 = 0.2D1*t6
    t14 = t31 + (t28 + t39)*t39 + (-t52 + t37)*t37 + (-t55 + t35)*t35
    t119 = sqrt((t36*t29 + t38*t30 - 0.2D1*t1*t40 + TINY + t18) &
        *(TINY + t14)*(t1*t28 + TINY + t6*t29 + t5*t30 + (-0.2D1 &
        *t12 - 0.2D1*t11 - 0.2D1*t7 - 0.2D1*t9 - 0.2D1*t10 - &
        0.2D1*t8 + (t52 - t26)*t20 + (t55 - t27)*t21 + (t34 - t25)* &
        t19 + (t14 - t24*t21 - t23*t20 - t22*t19 + t18)*w)*w + &
        0.2D1*(t27*t6 + t25*t4 + t26*t5 + ((-t3 + t24)*t36 + (-t1 &
        + t22)*t40 + (-t2 + t23)*t38)*w)*w + t31 + t32))

    value4 = (t17*t1 + t15*t2 + t16*t3 + (-t24*t2 + t23*t3 - t17) &
        *t40 + (t24*t1 - t22*t3 - t15)*t38 + (-t23*t1 + t22* &
        t2 - t16)*t36)/t119

end module ncsu_cv_HANDEDNESS

subroutine derivative4(w, x1, x2, x3, x4, d1, d2, d3, d4)

    implicit none

    NCSU_REAL, intent(in) :: w
    NCSU_REAL, intent(in) :: x1(3)
    NCSU_REAL, intent(in) :: x2(3)
    NCSU_REAL, intent(in) :: x3(3)
    NCSU_REAL, intent(in) :: x4(3)

    NCSU_REAL, intent(out) :: d1(3)
    NCSU_REAL, intent(out) :: d2(3)
    NCSU_REAL, intent(out) :: d3(3)
    NCSU_REAL, intent(out) :: d4(3)

    NCSU_REAL ::  t216
    NCSU_REAL ::  t41
    NCSU_REAL ::  t42
    NCSU_REAL ::  t43
    NCSU_REAL ::  t40
    NCSU_REAL ::  t128
    NCSU_REAL ::  t184
    NCSU_REAL ::  t191
    NCSU_REAL ::  t20
    NCSU_REAL ::  t93
    NCSU_REAL ::  t307
    NCSU_REAL ::  t19
    NCSU_REAL ::  t18
    NCSU_REAL ::  s1
    NCSU_REAL ::  t17
    NCSU_REAL ::  t106
    NCSU_REAL ::  t83
    NCSU_REAL ::  t85
    NCSU_REAL ::  t124
    NCSU_REAL ::  df(27)
    NCSU_REAL ::  t55
    NCSU_REAL ::  t16
    NCSU_REAL ::  t15
    NCSU_REAL ::  t38
    NCSU_REAL ::  t39
    NCSU_REAL ::  t14
    NCSU_REAL ::  t110
    NCSU_REAL ::  t100
    NCSU_REAL ::  t59
    NCSU_REAL ::  t60
    NCSU_REAL ::  t87
    NCSU_REAL ::  t52
    NCSU_REAL ::  t137
    NCSU_REAL ::  t1
    NCSU_REAL ::  t2
    NCSU_REAL ::  t114
    NCSU_REAL ::  t117
    NCSU_REAL ::  t125
    NCSU_REAL ::  t122
    NCSU_REAL ::  t146
    NCSU_REAL ::  t295
    NCSU_REAL ::  t13
    NCSU_REAL ::  t33
    NCSU_REAL ::  t22
    NCSU_REAL ::  t21
    NCSU_REAL ::  t34
    NCSU_REAL ::  t35
    NCSU_REAL ::  t321
    NCSU_REAL ::  t3
    NCSU_REAL ::  t174
    NCSU_REAL ::  t36
    NCSU_REAL ::  t37
    NCSU_REAL ::  t95
    NCSU_REAL ::  s0
    NCSU_REAL ::  t6
    NCSU_REAL ::  t7
    NCSU_REAL ::  t123
    NCSU_REAL ::  t9
    NCSU_REAL ::  t11
    NCSU_REAL ::  t4
    NCSU_REAL ::  t5
    NCSU_REAL ::  t32
    NCSU_REAL ::  t8
    NCSU_REAL ::  t10
    NCSU_REAL ::  t12
    NCSU_REAL ::  t31
    NCSU_REAL ::  t30
    NCSU_REAL ::  t29
    NCSU_REAL ::  t28
    NCSU_REAL ::  t179
    NCSU_REAL ::  t27
    NCSU_REAL ::  t237
    NCSU_REAL ::  t201
    NCSU_REAL ::  t204
    NCSU_REAL ::  t206
    NCSU_REAL ::  t79
    NCSU_REAL ::  t291
    NCSU_REAL ::  t293
    NCSU_REAL ::  t26
    NCSU_REAL ::  t25
    NCSU_REAL ::  t58
    NCSU_REAL ::  t24
    NCSU_REAL ::  t23
    NCSU_REAL ::  t226
    NCSU_REAL ::  t310

    t1 = x2(3)
    t7 = t1**2
    t2 = x2(2)
    t9 = t2**2
    t3 = x2(1)
    t11 = t3**2
    t32 = t7 + t9 + t11
    t4 = x3(3)
    t8 = t4**2
    t5 = x3(2)
    t10 = t5**2
    t6 = x3(1)
    t12 = t6**2
    t31 = t8 + t10 + t12
    t13 = 0.2D1*t2
    t30 = -t13
    t33 = 0.2D1*t3
    t29 = -t33
    t34 = 0.2D1*t4
    t28 = -t34
    t35 = x4(1)
    t36 = x1(1)
    t27 = t35 - t36
    t37 = x4(2)
    t38 = x1(2)
    t26 = t37 - t38
    t39 = x4(3)
    t40 = x1(3)
    t25 = t39 - t40
    t24 = -t35 + t6
    t23 = -t37 + t5
    t22 = t4 - t39
    t21 = t33
    t20 = t13
    t19 = 0.2D1*t1
    t41 = t40**2
    t42 = t38**2
    t43 = t36**2
    t18 = t41 + t42 + t43 + t32
    t17 = -t6*t37 + t5*t35
    t16 = -t5*t39 + t4*t37
    t15 = t6*t39 - t4*t35
    t52 = 0.2D1*t5
    t55 = 0.2D1*t6
    t14 = t31 + (t28 + t39)*t39 + (-t52 + t37)*t37 + (-t55 + t35)*t35
    t58 = t1*t28
    t59 = t6*t29
    t60 = t5*t30
    t79 = (-0.2D1*t12 - 0.2D1*t11 - 0.2D1*t7 - 0.2D1*t9 - 0.2D1*t10 &
        - 0.2D1*t8 + (t52 - t26)*t20 + (t55 - t27)*t21 + &
        (t34 - t25)*t19 + (t14 - t24*t21 - t23*t20 - t22*t19 + t18)*w)*w
    t83 = -t3 + t24
    t85 = -t1 + t22
    t87 = -t2 + t23
    t93 = 0.2D1*(t27*t6 + t25*t4 + t26*t5 + &
        (t83*t36 + t85*t40 + t87*t38)*w)*w
    t95 = TINY + t14
    s0 = (t58 + TINY + t59 + t60 + t79 + t93 + t31 + t32)*t95
    t100 = t36*t29 + t38*t30 - 0.2D1*t1*t40 + TINY + t18
    s1 = t100*s0
    t106 = -t24*t2 + t23*t3 - t17
    t110 = t24*t1 - t22*t3 - t15
    t114 = -t23*t1 + t22*t2 - t16
    t117 = sqrt(s1)
    df(27) = -(t17*t1 + t15*t2 + t16*t3 + t106*t40 + t110*t38 &
        + t114*t36)/t117/s1/0.2D1
    t122 = df(27)
    df(26) = t122*t100
    t123 = df(26)
    t124 = w**2
    t125 = t124*t95
    df(25) = t123*(t125 + t58 + TINY + t59 + t60 + t79 + t93 + t31 + t32)
    t128 = 0.1D1/t117
    df(24) = (t2 - t38)*t128
    df(23) = (t3 - t36)*t128
    df(22) = (t1 - t40)*t128
    df(21) = t122*s0 + t123*t124*t95
    t137 = w*t95
    df(20) = t123*(t34 - t25 - t22*w)*t137
    df(19) = t123*(t52 - t26 - t23*w)*t137
    df(18) = t123*(t55 - t27 - t24*w)*t137
    t146 = 0.2D1*t40*t124
    df(17) = t123*(-t19*t124 + t146)*t95 + (-t3*t38 + t2*t36)*t128
    df(16) = t123*(-t20*t124 + 0.2D1*t38*t124)*t95 + &
        (t3*t40 - t1*t36)*t128
    df(15) = t123*(-t21*t124 + 0.2D1*t36*t124)*t95 + &
        (-t2*t40 + t1*t38)*t128
    t174 = t19*w
    df(14) = t123*(-t174 + 0.2D1*t4*w)*t95
    t179 = t20*w
    df(13) = t123*(-t179 + 0.2D1*t5*w)*t95
    t184 = t21*w
    df(12) = t123*(-t184 + 0.2D1*t6*w)*t95
    t191 = df(25)
    df(11) = t123*t1*t95 + t191*t39
    df(10) = t122*t36*s0 + t123*t6*t95
    df(9) = t122*t38*s0 + t123*t5*t95
    t201 = t123*t95
    df(8) = t201 + t191
    t204 = 0.2D1*t123*w*t95
    df(7) = -t204 + df(8)
    df(6) = df(7)
    df(5) = df(6)
    t206 = df(21)
    df(4) = t201 + t206
    df(3) = -t204 + df(4)
    df(2) = df(3)
    df(1) = df(2)
    t216 = df(12)
    t226 = df(13)
    t237 = df(14)
    t291 = df(24)
    t293 = df(22)
    t295 = df(15)
    t307 = df(23)
    t310 = df(16)
    t321 = df(17)
    d1(1) = t114*t128 + t122*t29*s0 + 0.2D1*t123*t83*t125 &
        + 0.2D1*t206*t36 - t216
    d1(2) = t110*t128 + t122*t30*s0 + 0.2D1*t123*t87*t125 &
        + 0.2D1*t206*t38 - t226
    d1(3) = t106*t128 - 0.2D1*t122*t1*s0 + 0.2D1*t123*t85*t125 &
        + 0.2D1*t206*t40 - t237
    d2(1) = (t16 + t23*t40 - t22*t38)*t128 - 0.2D1*t123*t36*t125 &
        + 0.2D1*df(18) - 0.2D1*df(10) + 0.2D1*df(2)*t3
    d2(2) = (t15 - t24*t40 + t22*t36)*t128 - 0.2D1*t123*t38*t125 &
        + 0.2D1*df(19) - 0.2D1*df(9) + 0.2D1*df(1)*t2
    d2(3) = (t17 + t24*t38 - t23*t36)*t128 - 0.2D1*t122* &
        t40*s0 + t123*(t28 - t146)*t95 + 0.2D1*df(20) + 0.2D1*df(1)*t1
    d3(1) = t123*(t29 + 0.2D1*t184 + 0.2D1*t27*w)*t95 - &
        0.2D1*t191*t35 + t291*t39 - t293*t37 + t295 + 0.2D1*df(6)*t6
    d3(2) = t123*(t30 + 0.2D1*t179 + 0.2D1*t26*w)*t95 - &
        0.2D1*t191*t37 - t307*t39 + t293*t35 + t310 + 0.2D1*df(5)*t5
    d3(3) = t123*(0.2D1*t174 + 0.2D1*t25*w)*t95 - t291 &
        *t35 + t307*t37 + t321 - 0.2D1*df(11) + 0.2D1*df(5)*t4
    d4(1) = -0.2D1*t191*t24 - t291*t4 + t293*t5 - t295 + t216
    d4(2) = -0.2D1*t191*t23 + t307*t4 - t293*t6 - t310 + t226
    d4(3) = t191*(0.2D1*t39 + t28) + t291*t6 - t307*t5 - t321 + t237

end subroutine derivative4

end module ncsu_cv_HANDEDNESS
