! <compile=optimized>
!----------------------------------------------------------------------
! Copyright (C) 2004, 2005 Chaok Seok, Evangelos Coutsias and Ken Dill
!      UCSF, Univeristy of New Mexico, Seoul National University
! Witten by Chaok Seok and Evangelos Coutsias 2004.
! modified by Mahmoud Moradi and Volodymyr Babin at NCSU, 2007
! The original version is at http://www.dillgroup.ucsf.edu/rmsd

! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!

! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!

! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
! MA  02110-1301  USA
!-----------------------------------------------------------------------

!
! make up by Mahmoud Moradi and Volodymyr Babin at NCSU, Cox 308
!

#include "ncsu-utils.h"
#include "ncsu-config.h"

module ncsu_rmsd

!=============================================================================

    implicit none

    private

!=============================================================================

#ifdef NCSU_ENABLE_RMSD_CANNED
    public :: rmsd_canned
#endif /* NCSU_ENABLE_RMSD_CANNED */

    public :: rmsd_q
    public :: rmsd_q1
    public :: rmsd_q2u
    public :: rmsd_q3u
!=============================================================================

contains

!=============================================================================

!
! a "driver" in FORTRAN slang (for testing only)
!

#ifdef NCSU_RMSD_CANNED
    NCSU_REAL function rmsd_canned(m, x1, x2, g)

    use ncsu_utils
    use ncsu_constants

    implicit none

    NCSU_REAL, intent(in)    :: m(:)
    NCSU_REAL, intent(inout) :: x1(:), x2(:) ! destroyed upon return

    NCSU_REAL, intent(out), optional :: g(:)

    integer   :: n, a, a3, i
    NCSU_REAL :: mass, c1(3), c2(3), q(4), lambda, x1n, x2n, tmp, U(3, 3)

    n = size(m)

    ncsu_assert(n > 1)

    ncsu_assert(3*n == size(x1))
    ncsu_assert(3*n == size(x2))

    !
    ! find centers of mass
    !

    c1 = ZERO
    c2 = ZERO

    mass = ZERO

    do a = 1, n
        a3 = 3*(a - 1)
        mass = mass + m(a)
        do i = 1, 3
            c1(i) = c1(i) + m(a)*x1(a3 + i)
            c2(i) = c2(i) + m(a)*x2(a3 + i)
        end do
    end do

    ncsu_assert(mass > ZERO)

    c1 = c1/mass
    c2 = c2/mass

    ! center x1, x2 && find "norms"

    x1n = ZERO
    x2n = ZERO

    do a = 1, n
        a3 = 3*(a - 1)
        do i = 1, 3
            x1(a3 + i) = x1(a3 + i) - c1(i)
            x1n = x1n + m(a)*x1(a3 + i)**2
            x2(a3 + i) = x2(a3 + i) - c2(i)
            x2n = x2n + m(a)*x2(a3 + i)**2
        end do
    end do

    call rmsd_q(n, m, x1, x2, lambda, q)

    rmsd_canned = sqrt(max(ZERO, ((x1n + x2n) - 2*lambda))/mass)

    ! g is w.r.t x1

    if (present(g)) then
        ncsu_assert(3*n == size(g))
        call rmsd_q2u(q, U)

        tmp = ONE/mass/max(rmsd_canned, NCSU_TO_REAL(0.000001))
        do a = 1, n
            a3 = 3*a
            g(a3 - 2:a3) = m(a)*tmp*(x1(a3 - 2:a3) - matmul(U, x2(a3 - 2:a3)))
        end do
    end if ! present g

end module ncsu_rmsd
#endif /* NCSU_RMSD_CANNED */

!=============================================================================

!
! finds optimal rotation; size(w) == n; size(x?) == 3*n
!

subroutine rmsd_q(n, w, x1, x2, lambda, q)

    use ncsu_constants

    implicit none

    integer, intent(in) :: n
    NCSU_REAL, intent(in) :: w(*), x1(*), x2(*)

    NCSU_REAL, intent(out) :: lambda, q(4)

    integer   :: a, a3, i, j
    NCSU_REAL :: R(4, 4), S(4, 4)

    ! calculate the R matrix

    R = ZERO

    do a = 0, n - 1
        a3 = 3*a
        do i = 1, 3
            do j = 1, 3
                R(i, j) = R(i, j) + w(a + 1)*x1(a3 + i)*x2(a3 + j)
            end do
        end do
    end do

    ! S matrix

    S(1, 1) = R(1, 1) + R(2, 2) + R(3, 3)
    S(2, 1) = R(2, 3) - R(3, 2)
    S(3, 1) = R(3, 1) - R(1, 3)
    S(4, 1) = R(1, 2) - R(2, 1)

    S(1, 2) = S(2, 1)
    S(2, 2) = R(1, 1) - R(2, 2) - R(3, 3)
    S(3, 2) = R(1, 2) + R(2, 1)
    S(4, 2) = R(1, 3) + R(3, 1)

    S(1, 3) = S(3, 1)
    S(2, 3) = S(3, 2)
    S(3, 3) = -R(1, 1) + R(2, 2) - R(3, 3)
    S(4, 3) = R(2, 3) + R(3, 2)

    S(1, 4) = S(4, 1)
    S(2, 4) = S(4, 2)
    S(3, 4) = S(4, 3)
    S(4, 4) = -R(1, 1) - R(2, 2) + R(3, 3)

    call dstmev(S, lambda, q)

end subroutine rmsd_q

!=============================================================================

! add variable control

subroutine rmsd_q1(n, state, w, x1, x2, lambda, q)

    use ncsu_constants

    implicit none

    integer, intent(in) :: n, state(*)
    NCSU_REAL, intent(in) :: w(*), x1(*), x2(*)

    NCSU_REAL, intent(out) :: lambda, q(4)

    integer   :: a, a3, i, j
    NCSU_REAL :: R(4, 4), S(4, 4)

    ! calculate the R matrix

    R = ZERO

    do a = 0, n - 1
        if (state(a + 1) == 0) cycle
        a3 = 3*a
        do i = 1, 3
            do j = 1, 3
                R(i, j) = R(i, j) + w(a + 1)*x1(a3 + i)*x2(a3 + j)
            end do
        end do
    end do

    ! S matrix

    S(1, 1) = R(1, 1) + R(2, 2) + R(3, 3)
    S(2, 1) = R(2, 3) - R(3, 2)
    S(3, 1) = R(3, 1) - R(1, 3)
    S(4, 1) = R(1, 2) - R(2, 1)

    S(1, 2) = S(2, 1)
    S(2, 2) = R(1, 1) - R(2, 2) - R(3, 3)
    S(3, 2) = R(1, 2) + R(2, 1)
    S(4, 2) = R(1, 3) + R(3, 1)

    S(1, 3) = S(3, 1)
    S(2, 3) = S(3, 2)
    S(3, 3) = -R(1, 1) + R(2, 2) - R(3, 3)
    S(4, 3) = R(2, 3) + R(3, 2)

    S(1, 4) = S(4, 1)
    S(2, 4) = S(4, 2)
    S(3, 4) = S(4, 3)
    S(4, 4) = -R(1, 1) - R(2, 2) + R(3, 3)

    call dstmev(S, lambda, q)

end subroutine rmsd_q1

!=============================================================================

!
! This subroutine constructs (transposed) rotation matrix U from quaternion q.
!

subroutine rmsd_q2u(q, U)

    use ncsu_constants

    implicit none

    NCSU_REAL, intent(in)  :: q(*) ! 4
    NCSU_REAL, intent(out) :: U(3, 3)

    NCSU_REAL :: b0, b1, b2, b3
    NCSU_REAL :: q00, q01, q02, q03
    NCSU_REAL :: q11, q12, q13, q22, q23, q33

    b0 = q(1) + q(1)
    b1 = q(2) + q(2)
    b2 = q(3) + q(3)
    b3 = q(4) + q(4)

    q00 = b0*q(1) - ONE
    q01 = b0*q(2)
    q02 = b0*q(3)
    q03 = b0*q(4)

    q11 = b1*q(2)
    q12 = b1*q(3)
    q13 = b1*q(4)

    q22 = b2*q(3)
    q23 = b2*q(4)

    q33 = b3*q(4)

    U(1, 1) = q00 + q11
    U(2, 1) = q12 - q03
    U(3, 1) = q13 + q02

    U(1, 2) = q12 + q03
    U(2, 2) = q00 + q22
    U(3, 2) = q23 - q01

    U(1, 3) = q13 - q02
    U(2, 3) = q23 + q01
    U(3, 3) = q00 + q33

end subroutine rmsd_q2u

!=============================================================================
!
! This subroutine constructs (UNtransposed) rotation matrix U from quaternion q.
!

subroutine rmsd_q3u(q, U)

    use ncsu_constants

    implicit none

    NCSU_REAL, intent(in)  :: q(*) ! 4
    NCSU_REAL, intent(out) :: U(3, 3)

    NCSU_REAL :: b0, b1, b2, b3
    NCSU_REAL :: q00, q01, q02, q03
    NCSU_REAL :: q11, q12, q13, q22, q23, q33

    b0 = q(1) + q(1)
    b1 = q(2) + q(2)
    b2 = q(3) + q(3)
    b3 = q(4) + q(4)

    q00 = b0*q(1) - ONE
    q01 = b0*q(2)
    q02 = b0*q(3)
    q03 = b0*q(4)

    q11 = b1*q(2)
    q12 = b1*q(3)
    q13 = b1*q(4)

    q22 = b2*q(3)
    q23 = b2*q(4)

    q33 = b3*q(4)

    U(1, 1) = q00 + q11
    U(1, 2) = q12 - q03
    U(1, 3) = q13 + q02

    U(2, 1) = q12 + q03
    U(2, 2) = q00 + q22
    U(2, 3) = q23 - q01

    U(3, 1) = q13 - q02
    U(3, 2) = q23 + q01
    U(3, 3) = q00 + q33

end subroutine rmsd_q3u

!=============================================================================

!
! a simple subroutine to compute the leading eigenvalue and eigenvector
! of a symmetric, traceless 4x4 matrix A by an inverse power iteration:
! (1) the matrix is converted to tridiagonal form by 3 Givens
! rotations;  V*A*V' = T
! (2) Gershgorin's theorem is used to estimate a lower
! bound for the leading negative eigenvalue:
! lambda_1 > g=min(T11-t12,-t21+T22-t23,-t32+T33-t34,-t43+T44)
!          =
! where tij=abs(Tij)
! (3) Form the positive definite matrix
!     B = T-gI
! (4) Use svd (algorithm svdcmp from "Numerical Recipes")
!     to compute eigenvalues and eigenvectors for SPD matrix B
! (5) Shift spectrum back and keep leading singular vector
!     and largest eigenvalue.
! (6) Convert eigenvector to original matrix A, through
!     multiplication by V'.
!

subroutine dstmev(A, lambda, evec)

    implicit none

    NCSU_REAL, intent(in)  :: A(4, 4)
    NCSU_REAL, intent(out) :: lambda, evec(4)

    NCSU_REAL :: T(4, 4), V(4, 4), SV(4, 4), SW(4)

    integer :: i, max_loc(1)

    ! (I) Convert to tridiagonal form, keeping similarity transform
    !            (a product of 3 Givens rotations)
    call givens4(A, T, V)

    ! (II) Estimate lower bound of smallest eigenvalue by Gershgorin's theorem
    lambda = min(T(1, 1) - abs(T(1, 2)), &
        T(2, 2) - abs(T(2, 1)) - abs(T(2, 3)), &
        T(3, 3) - abs(T(3, 2)) - abs(T(3, 4)), &
        T(4, 4) - abs(T(4, 3)))

    ! (III) Form positive definite matrix  T = lambda*I - T
    do i = 1, 4
        T(i, i) = T(i, i) - lambda
    end do

    ! (IV) Compute singular values/vectors of SPD matrix B
    call svdcmp(T, SW, SV)

    ! (V) Shift spectrum back
    max_loc = maxloc(SW)
    lambda = lambda + SW(max_loc(1))

    ! (VI) Convert eigenvector to original coordinates: (V is transposed!)
    evec = matmul(V, SV(:, max_loc(1)))

end subroutine dstmev

!=============================================================================

!
! performs givens rotations to reduce symmetric 4x4 matrix to tridiagonal
!

subroutine givens4(S, T, V)

    use ncsu_constants

    implicit none

    NCSU_REAL, dimension(4, 4), intent(in)  :: S
    NCSU_REAL, dimension(4, 4), intent(out) :: T, V

    NCSU_REAL :: c1, c2, c3, s1, s2, s3, r1, r2, r3, c1c2, s1c2

    T = S
    V = ZERO

    ! zero out entries T(4,1) and T(1,4)
    ! compute cos and sin of rotation angle in the 3-4 plane

    r1 = pythag(T(3, 1), T(4, 1))

    if (r1 .ne. ZERO) then
        c1 = T(3, 1)/r1
        s1 = T(4, 1)/r1

        V(3, 3) = c1
        V(3, 4) = s1
        V(4, 3) = -s1
        V(4, 4) = c1

        T(3, 1) = r1
        T(4, 1) = ZERO

        T(3:4, 2:4) = matmul(V(3:4, 3:4), T(3:4, 2:4))
        T(1:2, 3:4) = transpose(T(3:4, 1:2))
        T(3:4, 3:4) = matmul(T(3:4, 3:4), transpose(V(3:4, 3:4)))
    else
        c1 = ONE
        s1 = ZERO
    end if

    ! zero out entries T(3,1) and T(1,3)
    ! compute cos and sin of rotation angle in the 2-3 plane

    r2 = pythag(T(3, 1), T(2, 1))

    if (r2 .ne. ZERO) then
        c2 = T(2, 1)/r2
        s2 = T(3, 1)/r2

        V(2, 2) = c2
        V(2, 3) = s2
        V(3, 2) = -s2
        V(3, 3) = c2

        T(2, 1) = r2
        T(3, 1) = ZERO

        T(2:3, 2:4) = matmul(V(2:3, 2:3), T(2:3, 2:4))
        T(1, 2:3) = T(2:3, 1)
        T(4, 2:3) = T(2:3, 4)
        T(2:3, 2:3) = matmul(T(2:3, 2:3), transpose(V(2:3, 2:3)))
    else
        c2 = ONE
        s2 = ZERO
    end if

    ! zero out entries T(4,2) and T(2,4)
    ! compute cos and sin of rotation angle in the 3-4 plane

    r3 = pythag(T(4, 2), T(3, 2))

    if (r3 .ne. ZERO) then
        c3 = T(3, 2)/r3
        s3 = T(4, 2)/r3

        V(3, 3) = c3
        V(3, 4) = s3
        V(4, 3) = -s3
        V(4, 4) = c3

        T(3, 2) = r3
        T(4, 2) = ZERO

        T(3:4, 3:4) = matmul(V(3:4, 3:4), T(3:4, 3:4))
        T(1:2, 3:4) = transpose(T(3:4, 1:2))
        T(3:4, 3:4) = matmul(T(3:4, 3:4), transpose(V(3:4, 3:4)))
    else
        c3 = ONE
        s3 = ZERO
    end if

    ! compute net rotation matrix (accumulate similarity for evec. computation)
    ! To save transposing later, This is the transpose!

    V(1, 1) = ONE
    V(1, 2:4) = ZERO
    V(2:4, 1) = ZERO

    V(2, 2) = c2
    V(3, 2) = c1*s2
    V(4, 2) = s1*s2

    c1c2 = c1*c2
    s1c2 = s1*c2

    V(2, 3) = -s2*c3
    V(3, 3) = c1c2*c3 - s1*s3
    V(4, 3) = s1c2*c3 + c1*s3
    V(2, 4) = s2*s3
    V(3, 4) = -c1c2*s3 - s1*c3
    V(4, 4) = -s1c2*s3 + c1*c3

end subroutine givens4

!============================================================================

subroutine svdcmp(a, w, v)

    use ncsu_utils
    use ncsu_constants

    implicit none

    integer, parameter :: N = 4

    NCSU_REAL, intent(inout) :: a(N, *)
    NCSU_REAL, intent(out)   :: v(N, *), w(*)

    integer :: i, its, j, jj, k, l, nm

    NCSU_REAL :: anorm, c, f, g, h, s, scale, x, y, z, rv1(2*N)

    g = ZERO
    scale = ZERO
    anorm = ZERO

    nm = 0 ! for g95

    do i = 1, N

        l = i + 1
        rv1(i) = scale*g

        g = ZERO
        s = ZERO
        scale = ZERO

        do k = i, N
            scale = scale + abs(a(k, i))
        end do

        if (scale .ne. ZERO) then
            do k = i, N
                a(k, i) = a(k, i)/scale
                s = s + a(k, i)*a(k, i)
            end do

            f = a(i, i)
            g = -sign(sqrt(s), f)
            h = f*g - s
            a(i, i) = f - g

            do j = l, N
                s = ZERO
                do k = i, N
                    s = s + a(k, i)*a(k, j)
                end do

                f = s/h
                do k = i, N
                    a(k, j) = a(k, j) + f*a(k, i)
                end do
            end do

            do k = i, N
                a(k, i) = scale*a(k, i)
            end do
        end if ! scale .ne. ZERO

        w(i) = scale*g
        g = ZERO
        s = ZERO
        scale = ZERO

        if (i .ne. N) then
            do k = l, N
                scale = scale + abs(a(i, k))
            end do
            if (scale .ne. ZERO) then
                do k = l, N
                    a(i, k) = a(i, k)/scale
                    s = s + a(i, k)*a(i, k)
                end do
                f = a(i, l)
                g = -sign(sqrt(s), f)
                h = f*g - s
                a(i, l) = f - g
                do k = l, N
                    rv1(k) = a(i, k)/h
                end do
                do j = l, N
                    s = ZERO
                    do k = l, N
                        s = s + a(j, k)*a(i, k)
                    end do
                    do k = l, N
                        a(j, k) = a(j, k) + s*rv1(k)
                    end do
                end do
                do k = l, N
                    a(i, k) = scale*a(i, k)
                end do
            end if
        end if
        anorm = max(anorm, (abs(w(i)) + abs(rv1(i))))
    end do

    do i = N, 1, -1
        if (i .lt. N) then
            if (g .ne. ZERO) then
                do j = l, N
                    v(j, i) = (a(i, j)/a(i, l))/g
                end do
                do j = l, N
                    s = ZERO
                    do k = l, N
                        s = s + a(i, k)*v(k, j)
                    end do
                    do k = l, N
                        v(k, j) = v(k, j) + s*v(k, i)
                    end do
                end do
            end if ! g .ne. ZERO
            do j = l, N
                v(i, j) = ZERO
                v(j, i) = ZERO
            end do
        end if
        v(i, i) = ONE
        g = rv1(i)
        l = i
    end do

    do i = N, 1, -1
        l = i + 1
        g = w(i)
        do j = l, N
            a(i, j) = ZERO
        end do
        if (g .ne. ZERO) then
            g = ONE/g
            do j = l, N
                s = ZERO
                do k = l, N
                    s = s + a(k, i)*a(k, j)
                end do
                f = (s/a(i, i))*g
                do k = i, N
                    a(k, j) = a(k, j) + f*a(k, i)
                end do
            end do
            do j = i, N
                a(j, i) = a(j, i)*g
            end do
        else
            do j = i, N
                a(j, i) = ZERO
            end do
        end if ! g .ne. ZERO
        a(i, i) = a(i, i) + ONE
    end do

    do k = N, 1, -1
        do its = 1, 30
            do l = k, 1, -1
                nm = l - 1
                if ((abs(rv1(l)) + anorm) .eq. anorm) &
                    goto 2
                if ((abs(w(nm)) + anorm) .eq. anorm) &
                    goto 1
            end do
1           c = ZERO
            s = ONE
            do i = l, k
                f = s*rv1(i)
                rv1(i) = c*rv1(i)
                if ((abs(f) + anorm) .eq. anorm) &
                    goto 2
                g = w(i)
                h = pythag(f, g)
                w(i) = h
                h = ONE/h
                c = (g*h)
                s = -(f*h)
                do j = 1, N
                    y = a(j, nm)
                    z = a(j, i)
                    a(j, nm) = (y*c) + (z*s)
                    a(j, i) = -(y*s) + (z*c)
                end do
            end do
2           z = w(k)
            if (l .eq. k) then
                if (z .lt. ZERO) then
                    w(k) = -z
                    do j = 1, N
                        v(j, k) = -v(j, k)
                    end do
                end if
                goto 3
            end if
            if (its .eq. 30) &
                call fatal('ncsu-rmsd: no convergence in svdcmp()')
            x = w(l)
            nm = k - 1
            y = w(nm)
            g = rv1(nm)
            h = rv1(k)
            f = ((y - z)*(y + z) + (g - h)*(g + h))/(2*h*y)
            g = pythag(f, ONE)
            f = ((x - z)*(x + z) + h*((y/(f + sign(g, f))) - h))/x
            c = ONE
            s = ONE
            do j = l, nm
                i = j + 1
                g = rv1(i)
                y = w(i)
                h = s*g
                g = c*g
                z = pythag(f, h)
                rv1(j) = z
                c = f/z
                s = h/z
                f = (x*c) + (g*s)
                g = -(x*s) + (g*c)
                h = y*s
                y = y*c
                do jj = 1, N
                    x = v(jj, j)
                    z = v(jj, i)
                    v(jj, j) = (x*c) + (z*s)
                    v(jj, i) = -(x*s) + (z*c)
                end do
                z = pythag(f, h)
                w(j) = z
                if (z .ne. ZERO) then
                    z = ONE/z
                    c = f*z
                    s = h*z
                end if
                f = (c*g) + (s*y)
                x = -(s*g) + (c*y)
                do jj = 1, N
                    y = a(jj, j)
                    z = a(jj, i)
                    a(jj, j) = (y*c) + (z*s)
                    a(jj, i) = -(y*s) + (z*c)
                end do
            end do
            rv1(l) = ZERO
            rv1(k) = f
            w(k) = x
        end do
3       continue
    end do

end subroutine svdcmp

!============================================================================

!
! computes sqrt(a**2 + b**2) carefully
!

pure NCSU_REAL function pythag(a, b)

use ncsu_constants

implicit none

NCSU_REAL, intent(in) :: a, b

NCSU_REAL :: absa, absb

absa = abs(a)
absb = abs(b)

if (absa .gt. absb) then
    pythag = absa*sqrt(ONE + (absb/absa)**2)
else
    if (absb .eq. ZERO) then
        pythag = ZERO
    else
        pythag = absb*sqrt(ONE + (absa/absb)**2)
    end if
end if

end function pythag

!============================================================================

end module ncsu_rmsd
