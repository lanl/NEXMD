! <compile=optimized>
!  -*- mode: f90; coding: iso-8859-15; -*-

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#include "dprec.fh"
#include "copyright.h"
subroutine broyden(niter, natoms, qmold, qmat)
    ! Calculates the next iteration qmat vector using
    ! Broyden mixing as described in:
    !
    !     D. D. Johnson, Phys. Rev. B, 38, 12807 (1988)
    !
    ! The original version of this routine is from
    ! Elstner, Porezag and Hajnal and was
    ! distributed with the SCC-DFTB code.
    !
    ! Adapted for Fortran90/Amber by:
    !  Gustavo Seabra
    !  Quantum Theory Project
    !  University of Florida

    use constants, only : zero, one
    use qm2_dftb_module, only : MAX_BRD_ITER, MAXSIZ, brd_struct
    implicit none

!! Broyden Mixing Parameter
    _REAL_, parameter :: amix = 0.2d0

!! Passed in:
    integer, intent(in)   :: niter       ! iteration number
    integer, intent(in)   :: natoms      ! Length of the vectors qmold and qmat
    _REAL_, intent(in)   :: qmat(*)     ! The most recent charges vector
    _REAL_, intent(inout) :: qmold(*)    ! In:  Old charges vector
    ! Out: New charges vector
!! Locals
    _REAL_  :: w0
    _REAL_  :: wtmp
    _REAL_  :: fnorm
    _REAL_  :: dfnorm
    _REAL_  :: fac1
    _REAL_  :: fac2
    _REAL_  :: aij, cmj, gamma

    integer :: lastit
    integer :: lastm1
    integer :: lastm2
    integer :: ilastit
    integer :: iter
    integer :: i, j, k, lm, ln, ip

!! Pointers
    ! These are continually updated:
    !   ui  -> Johnson's u(i)
    !   vti -> DeltaF(transpose)
    !
    ! The results from all iterations are stored in unit32
    ! this was originally done on tape, to avoid the (then)
    ! prohibitive storage of the whole jacobian.
    !
    !   t1  -> vt from earlier iteration
    !   f   -> F = qmat - qmold
    !   df  -> DeltaF = F(i) - F(i-1)
    !   dumvi -> Predicted vector, returned in qmold.

    _REAL_, pointer :: f(:)           !(maxsiz) qmat - qmold
    _REAL_, pointer :: ui(:)          !(maxsiz)
    _REAL_, pointer :: vti(:)         !(maxsiz)
    _REAL_, pointer :: t1(:)          !(maxsiz)
    _REAL_, pointer :: dumvi(:)       !(maxsiz)
    _REAL_, pointer :: df(:)          !(maxsiz)
    _REAL_, pointer :: cm(:)          !(imatsz)
    _REAL_, pointer :: w(:)           !(imatsz)
    _REAL_, pointer :: a(:, :)         !(imatsz,imatsz)
    _REAL_, pointer :: beta(:, :)      !(imatsz,imatsz)
    _REAL_, pointer :: d(:, :)         !(imatsz,imatsz)
    _REAL_, pointer :: unit31(:, :)    !(maxsiz,2)
    _REAL_, pointer :: unit32(:, :, :)  !(maxsiz,2,MAX_BRD_ITER)

    save !Saves every variable from one iteration to the next!

!! Pointers

    f => brd_struct%f
    ui => brd_struct%ui
    vti => brd_struct%vti
    t1 => brd_struct%t1
    dumvi => brd_struct%dumvi
    df => brd_struct%df
    cm => brd_struct%cm
    w => brd_struct%w
    a => brd_struct%a
    beta => brd_struct%b
    d => brd_struct%d
    unit31 => brd_struct%unit31
    unit32 => brd_struct%unit32

!!--

    iter = niter
    if (niter > MAX_BRD_ITER) iter = mod(iter, MAX_BRD_ITER) + 1

    if (ITER == 0) return

    !*******************  begin broyden's method  **********************

    ! weighting factor for the zeroth iteration
    w0 = 0.01d0

    ! f:      the difference of previous output and input vectors
    ! dumvi:  a dummy vector, here it is the previous input vector
    if (iter /= 1) then
        ! ==================================================
        !         ALL ITERATIONS EXCEPT THE FIRST
        ! ==================================================

        lastit = ilastit

        f(1:natoms) = unit31(1:natoms, 1) ! F(n-1) = F from the previous iteration
        dumvi(1:natoms) = unit31(1:natoms, 2) ! qmold(n-1)

        dumvi(1:natoms) = qmold(1:natoms) - dumvi(1:natoms)    ! qmold(n) - qmold(n-1)
        df(1:natoms) = qmat(1:natoms) - qmold(1:natoms) - f(1:natoms)

        f(1:natoms) = qmat(1:natoms) - qmold(1:natoms) ! F for this iteration.

        ! Norm of F and DF, used for normalization:
        ! for i-th iter.,
        !      dfnorm = | df(i) - df(i-1) |
        !       fnorm = |  f(i) -  f(i-1) |

        dfnorm = zero
        fnorm = zero
        do k = 1, natoms
            dfnorm = dfnorm + df(k)*df(k)
            fnorm = fnorm + f(k)*f(k)
        end do
        dfnorm = sqrt(dfnorm)
        fnorm = sqrt(fnorm)

        ! Normalization Factors
        fac2 = one/dfnorm
        fac1 = amix*fac2

        ! Normalized vectors
        ui(1:natoms) = fac1*df(1:natoms) + fac2*dumvi(1:natoms)
        vti(1:natoms) = fac2*df(1:natoms)

        !*********** calculation of coefficient matrices *************
        !***********    and the sum for corrections      *************

        lastit = lastit + 1
        lastm1 = lastit - 1
        lastm2 = lastit - 2

        ! dumvi is the u(of i) and t1 is the vt(of i)
        ! from the previous iterations
        if (lastit > 2) then
            do j = 1, lastm2

                dumvi(1:natoms) = unit32(1:natoms, 1, j)
                t1(1:natoms) = unit32(1:natoms, 2, j)

                aij = zero
                cmj = zero

                do k = 1, natoms
                    cmj = cmj + t1(k)*f(k)
                    aij = aij + t1(k)*vti(k)
                end do

                a(lastm1, j) = aij
                a(j, lastm1) = aij
                cm(j) = cmj

            end do
        end if

        aij = zero
        cmj = zero

        do k = 1, natoms
            cmj = cmj + vti(k)*f(k)
            aij = aij + vti(k)*vti(k)
        end do

        a(lastm1, lastm1) = aij
        cm(lastm1) = cmj

        unit32(1:natoms, 1, lastm1) = ui(1:natoms)
        unit32(1:natoms, 2, lastm1) = vti(1:natoms)

        ! ===================
        !  WEIGHTING FACTORS
        ! ===================
        ! the weighting factors for each iteration have been chosen
        ! equal to one over the r.m.s. error. this need not be the case.
        if (fnorm > 1.0d-7) then
            wtmp = 0.010d0/fnorm
        else
            wtmp = 1.0d5
        end if

        if (wtmp < one) wtmp = one

        w(lastm1) = wtmp
        ilastit = lastit

        unit31(1:natoms, 1) = f(1:natoms)
        unit31(1:natoms, 2) = qmold(1:natoms)

        ! *** set up and calculate beta matrix
        do lm = 1, lastm1
            do ln = 1, lastm1
                d(ln, lm) = a(ln, lm)*w(ln)*w(lm)
                beta(ln, lm) = zero
            end do
            beta(lm, lm) = one
            d(lm, lm) = w0**2 + a(lm, lm)*w(lm)*w(lm)
        end do

        call inverse(d, beta, lastm1)

        ! *** calculate the vector for the new iteration
        dumvi(1:natoms) = qmold(1:natoms) + amix*f(1:natoms)
        do i = 1, lastm1

            gamma = zero
            do ip = 1, lastm1
                gamma = gamma + cm(ip)*beta(ip, i)*w(ip)
            end do

            dumvi(1:natoms) = dumvi(1:natoms) - gamma*unit32(1:natoms, 1, i)*w(i)
        end do
        ! *** end of the calculation of dumvi, the new vector

    else ! if (iter /= 1) then

        ! ==================================================
        !            THIS IS THE FIRST ITERATION.
        ! ==================================================
        !
        ! In the first iteration, we just do a simple mixing
        ! of: F (=qmat-qmold) and qmold.

        ! Sets some variables (for the next iteration?)
        lastit = 1
        ilastit = lastit

        f(1:natoms) = qmat(1:natoms) - qmold(1:natoms)

        ! Store F(i) into unit31(i,1)
        unit31(1:natoms, 1) = f(1:natoms)

        ! Store qmold(i) into unit32(i,2)
        unit31(1:natoms, 2) = qmold(1:natoms)

        ! since we are on the first iteration, simple mix the vectors.
        dumvi(1:natoms) = qmold(1:natoms) + amix*f(1:natoms)

    end if ! if (iter /= 1) then

    !
    !*************  the end of the broyden method **************

    ! Finally, load the new vector into the appropriate arrays.
    qmold(1:natoms) = dumvi(1:natoms)

    return

end subroutine broyden
! =============================================================

! =============================================================
subroutine inverse(a, b, m)
    !
    ! Subroutine to perform Gaussian Elimination.
    !
    ! Gaussian Elimination:
    ! ---------------------
    !
    ! A method to solve matrix equations of the form A . x = b,
    ! If b is the identity matrix, the result will be the inverse
    ! of A.
    !
    ! The result is returned into B.
    !
    ! NO ZEROS ALONG THE DIAGONAL !!
    use constants, only : zero, one
    use qm2_dftb_module, only : IMATSZ, brd_struct
    implicit none

!!   ! Static memory
!!   integer, parameter :: imatsz = 80   ! Maximum matrix dimension

!! Passed in
    _REAL_, intent(inout) :: a(imatsz, imatsz)
    _REAL_, intent(inout) :: b(imatsz, imatsz)
    integer, intent(in) :: m

!! Locals
    _REAL_ :: atmp

    integer :: i, j, k, n

!! Pointers
    _REAL_, pointer :: td(:)    ! (imatsz)
    _REAL_, pointer :: ad(:)    ! (imatsz)
    _REAL_, pointer :: bd(:)    ! (imatsz)

    save !What is this for?

!!Pointers
    td => brd_struct%td
    ad => brd_struct%ad
    bd => brd_struct%bd

    n = m

    if (n > imatsz) then
        call sander_bomb("inverse <qm2_dftb_broyden.f> : ", "Matrix is too large.", "Exiting")
    end if

    ! Check if any diagonal element is zero. If so, bomb.
    do i = 1, n
        atmp = a(i, i)
        if (abs(atmp) < 1.0d-08) then
            write (6, '('' inverse: matrix to be inverted has zero diagonal element in the '',i4,'' row'')') i
            call sander_bomb("inverse <qm2_dftb_broyden.f> : ", "Zero diagonal element.", "Exiting")
        end if
    end do

    if (n /= 1) then
        do i = 1, n

            do j = 1, n
                td(j) = a(j, i)/a(i, i)
            end do

            td(i) = zero

            do k = 1, n
                bd(k) = b(i, k)
                ad(k) = a(i, k)
            end do

            do k = 1, n
                do j = 1, n
                    b(j, k) = b(j, k) - (td(j)*bd(k))
                    a(j, k) = a(j, k) - (td(j)*ad(k))
                end do
            end do

        end do

        do i = 1, n
            do j = 1, n
                b(j, i) = b(j, i)/a(j, j)
            end do
        end do

        return
    else
        b(1, 1) = one/a(1, 1)
    end if
    return
end subroutine inverse
