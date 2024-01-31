! <compile=optimized>

#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_colvar_math

    implicit none

    private

    public :: distance
    public :: distance_d

    public :: angle
    public :: angle_d

    public :: torsion
    public :: torsion_d

    public :: dot3
    public :: norm3
    public :: cross3

    public :: cosine
    public :: cosine_d

    private :: torsion3
    private :: torsion3_d

! leaving the pucker stuff here -- may be useful in the future
#if 0
    public :: pucker_z
    public :: pucker_dz
#endif

!=============================================================================

    NCSU_REAL, public, parameter :: epsilon = 1.0d-8

!=============================================================================

contains

!=============================================================================

    pure NCSU_REAL function distance(r1, r2)

    implicit none

    NCSU_REAL, intent(in) :: r1(3), r2(3)

    distance = norm3(r1 - r2)

end module ncsu_colvar_math

!=============================================================================

subroutine distance_d(r1, r2, d1, d2)

    implicit none

    NCSU_REAL, intent(in)  :: r1(3), r2(3)
    NCSU_REAL, intent(out) :: d1(3), d2(3)

    NCSU_REAL :: dr(3)

    dr = r1 - r2

    d1 = dr/norm3(dr)
    d2 = -d1

end subroutine distance_d

!=============================================================================

pure NCSU_REAL function angle(r1, r2, r3)

implicit none

NCSU_REAL, intent(in) :: r1(3), r2(3), r3(3)

angle = acos(cosine(r1 - r2, r3 - r2))

end function angle

!=============================================================================

subroutine angle_d(r1, r2, r3, d1, d2, d3)

    use ncsu_constants, only : ONE

    implicit none

    NCSU_REAL, intent(in)  :: r1(3), r2(3), r3(3)
    NCSU_REAL, intent(out) :: d1(3), d2(3), d3(3)

    NCSU_REAL :: c, d

    c = cosine_d(r1 - r2, r3 - r2, d1, d3)
    d = -ONE/sqrt(ONE - c*c)

    d1 = d*d1
    d3 = d*d3

    d2 = -(d1 + d3)

end subroutine angle_d

!=============================================================================

!
! for A == B == C == D, a torsion along B == C is arccos([ABxBC]*[CDx(-BC)])
!        and its sign is given by sign(BC*[[ABxBC]x[CDx(-BC)]])
!

pure NCSU_REAL function torsion(r1, r2, r3, r4)

implicit none

NCSU_REAL, intent(in) :: r1(3), r2(3), r3(3), r4(3)

torsion = torsion3(r2 - r1, r3 - r2, r4 - r3)

end function torsion

!=============================================================================

subroutine torsion_d(r1, r2, r3, r4, d1, d2, d3, d4)

    implicit none

    NCSU_REAL, intent(in)  :: r1(3), r2(3), r3(3), r4(3)
    NCSU_REAL, intent(out) :: d1(3), d2(3), d3(3), d4(3)

    NCSU_REAL :: t(3)

    call torsion3_d(r2 - r1, r3 - r2, r4 - r3, d1, t, d4)

    d2 = d1 - t
    d1 = -d1
    d3 = t - d4

end subroutine torsion_d

!=============================================================================

pure NCSU_REAL function dot3(v1, v2)

implicit none

NCSU_REAL, intent(in) :: v1(3), v2(3)

dot3 = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)

end function dot3

!=============================================================================

pure NCSU_REAL function norm3(v)

implicit none

NCSU_REAL, intent(in) :: v(3)

norm3 = sqrt(v(1)**2 + v(2)**2 + v(3)**2)

end function norm3

!=============================================================================

pure function cross3(v1, v2) result(cross)

    implicit none

    NCSU_REAL             :: cross(3)
    NCSU_REAL, intent(in) :: v1(3), v2(3)

    cross(1) = v1(2)*v2(3) - v1(3)*v2(2)
    cross(2) = v1(3)*v2(1) - v1(1)*v2(3)
    cross(3) = v1(1)*v2(2) - v1(2)*v2(1)

end function cross3

!=============================================================================

pure NCSU_REAL function cosine(v1, v2)

use ncsu_constants, only : ONE

implicit none

NCSU_REAL, intent(in) :: v1(3), v2(3)

cosine = dot3(v1, v2)/(norm3(v1)*norm3(v2))

if (cosine > ONE) cosine = ONE
if (cosine < -ONE) cosine = -ONE

end function cosine

!=============================================================================

NCSU_REAL function cosine_d(v1, v2, d1, d2)

use ncsu_constants, only : ONE

implicit none

NCSU_REAL, intent(in)  :: v1(3), v2(3)
NCSU_REAL, intent(out) :: d1(3), d2(3)

NCSU_REAL :: n1, n2, p12

n1 = norm3(v1)
n2 = norm3(v2)
p12 = dot3(v1, v2)

cosine_d = p12/(n1*n2)

d1 = (v2 - v1*p12/(n1**2))/(n1*n2)
d2 = (v1 - v2*p12/(n2**2))/(n1*n2)

if (cosine_d > ONE) cosine_d = ONE
if (cosine_d < -ONE) cosine_d = -ONE

end function cosine_d

!=============================================================================

pure NCSU_REAL function torsion3(v1, v2, v3)

implicit none

NCSU_REAL, intent(in) :: v1(3), v2(3), v3(3)

NCSU_REAL :: n1(3), n2(3)

n1 = cross3(v1, v2)
n2 = cross3(v2, v3)

torsion3 = sign(acos(cosine(n1, n2)), dot3(v2, cross3(n1, n2)))

end function torsion3

!=============================================================================

subroutine torsion3_d(v1, v2, v3, d1, d2, d3)

    use ncsu_constants, only : ONE

    implicit none

    NCSU_REAL, intent(in)  :: v1(3), v2(3), v3(3)
    NCSU_REAL, intent(out) :: d1(3), d2(3), d3(3)

    NCSU_REAL :: n1(3), n2(3), dc1(3), dc2(3), c, c2, s, d

    n1 = cross3(v1, v2)
    n2 = cross3(v2, v3)

    c = cosine_d(n1, n2, dc1, dc2)
    s = dot3(v2, cross3(n1, n2))

    c2 = c*c - epsilon
    d = sign(ONE/sqrt(ONE - c2), s)

    d1 = d*cross3(dc1, v2)
    d2 = d*(cross3(v1, dc1) + cross3(dc2, v3))
    d3 = d*cross3(v2, dc2)

end subroutine torsion3_d

!=============================================================================
#if 0
!
! computes out-of-plane (puckering) coordinates as described in
!        http://dx.doi.org/10.1021/ja00839a011
!

subroutine pucker_z(R, z)

    NCSU_USE_AFAILED

    use ncsu_constants

    implicit none

    NCSU_REAL, intent(in) :: R(:, :)
    NCSU_REAL, intent(out) :: z(:)

    NCSU_REAL :: R0(3), Rt(3), Rc(3), Rs(3), phase, N(3), pi
    integer :: nring, k

    nring = size(R, 2)
    ncsu_assert(nring .gt. 0)

    pi = 4*atan(NCSU_TO_REAL(1))

    R0 = ZERO
    do k = 1, nring
        R0 = R0 + R(:, k)
    end do
    R0 = R0/nring

    Rc = ZERO
    Rs = ZERO

    do k = 1, nring
        Rt = R(:, k) - R0
        phase = 2*pi*(k - 1)/NCSU_TO_REAL(nring)
        Rc = Rc + Rt*cos(phase)
        Rs = Rs + Rt*sin(phase)
    end do

    N = cross3(Rs, Rc)
    N = N/norm3(N)

    do k = 1, nring
        Rt = R(:, k) - R0
        z(k) = dot3(Rt, N)
    end do

end subroutine pucker_z

!-----------------------------------------------------------------------------

subroutine pucker_dz(R, z, dz)

#ifndef NCSU_DISABLE_ASSERT
    use ncsu_utils
#endif /* NCSU_DISABLE_ASSERT */

    use ncsu_constants

    implicit none

    NCSU_REAL, intent(in) :: R(:, :)
    NCSU_REAL, intent(out) :: z(:), dz(:, :, :)

    NCSU_REAL :: R0(3), Rt(3), Rc(3), Rs(3), pi
    NCSU_REAL :: N(3), phi(3), z0(3), NN, phase

    integer :: nring, k, m

    nring = size(R, 2)
    ncsu_assert(nring .gt. 0)

    pi = 4*atan(NCSU_TO_REAL(1))

    R0 = ZERO
    do k = 1, nring
        R0 = R0 + R(:, k)
    end do
    R0 = R0/nring

    Rc = ZERO
    Rs = ZERO

    do k = 1, nring
        Rt = R(:, k) - R0
        phase = 2*pi*(k - 1)/NCSU_TO_REAL(nring)
        Rc = Rc + Rt*cos(phase)
        Rs = Rs + Rt*sin(phase)
    end do

    N = cross3(Rs, Rc)
    NN = norm3(N)
    N = N/NN

    do m = 1, nring
        Rt = R(:, m) - R0
        z(m) = dot3(Rt, N)
        Rt = Rt - z(m)*N
        z0 = ZERO
        do k = 1, nring
            phase = 2*pi*(k - 1)/NCSU_TO_REAL(nring)
            phi = (Rs*cos(phase) - Rc*sin(phase))/NN
            dz(:, k, m) = cross3(Rt, phi)
            if (m .eq. k) &
                dz(:, k, m) = dz(:, k, m) + N
            z0 = z0 + dz(:, k, m)
        end do
        z0 = z0/nring
        do k = 1, nring
            dz(:, k, m) = dz(:, k, m) - z0
        end do
    end do

end subroutine pucker_dz
#endif

!=============================================================================

end module ncsu_colvar_math
