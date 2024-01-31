! <compile=optimized>

#include "dprec.fh"

! INTEGER FUNCTION SKSPAR
! =======================
!
! It is actually used to get the SK parameters,
! which are returned in 'dd' (_REAL_) after the call.
!
! That should be a subroutine, not a function!!
!
integer function skspar(i, j, r2, dd)

    use qm2_dftb_module, only : sktab, lmax

    implicit none
    integer:: in, ind, inu, maxmax, minmax
    _REAL_ :: cubicspline, grdr, r, spline5th

    integer :: i, j
    integer :: mxind
    _REAL_  :: r2, dd(13)
    _REAL_  :: x0, x1, x2, f0, f1, f2
    _REAL_  :: xh, hl

    skspar = 0.0d0

    ! Orbital limits
    maxmax = max(lmax(i), lmax(j))
    minmax = min(lmax(i), lmax(j))

    if (maxmax <= 1) then
        inu = 10
    else
        if (maxmax <= 2) then
            if (minmax <= 1) then
                inu = 9
            else
                inu = 6
            end if
        else
            if (minmax <= 1) then
                inu = 8
            else
                if (minmax <= 2) then
                    inu = 4
                else
                    inu = 1
                end if
            end if
        end if
    end if

    ! Maximum index
    ! In case a distance is larger than the
    ! distances stored, the program tries a
    ! 5th order spline extrapolation, as
    ! long as it is not too much larger than
    ! the largest distance stored, basically
    ! seting a confidence limit on the extrapolations.
    mxind = sktab%dimens(i, j) + (0.3/sktab%sr(i, j) - 1.0)
    r = sqrt(r2)

    ! Index for THE integral needed
    ind = r/sktab%sr(i, j) + 1.0

    if (r2 < 1e-8) then

        ! Same function, same atom: S_mu_nu = 1.0
        do in = 1, 3
            dd(in + 10) = 1.0
        end do

    else

        if (ind + 2 > sktab%dimens(i, j)) then ! The index + 2 is outside bounds

            if (ind < sktab%dimens(i, j)) then

                ! Case 1:
                !        Index+2 is outside the bounds,
                !        But index is still inside bounds...
                !
                !        Do a 5th order spline extrapolation

                ! Get the 3 last distances
                x0 = (sktab%dimens(i, j) - 3)*sktab%sr(i, j)
                x1 = x0 + sktab%sr(i, j)
                x2 = x1 + sktab%sr(i, j)

                ! Set spline parameters
                xh = r - x1
                hl = x2 - x1

                ! Calculate the terms by extrapolating a 5th order spline.
                do in = inu, 10
                    f0 = sktab%skstab(in, sktab%dimens(i, j) - 2, i, j)
                    f1 = sktab%skstab(in, sktab%dimens(i, j) - 1, i, j)
                    f2 = sktab%skstab(in, sktab%dimens(i, j), i, j)
                    dd(in) = cubicspline(f0, f1, f2, x0, x1, xh, hl, sktab%sr(i, j))
                end do

            else ! if(ind < sktab%dimens(i,j))then

                if ((ind >= sktab%dimens(i, j)) .AND. (ind < mxind)) then

                    ! Case 2:
                    !        Index+2 is outside the bounds,
                    !        AND so is Index, but it is still
                    !        not too far from the last one.
                    !
                    !        Do a 5th order spline extrapolation.

                    ! Again, get the last 3 distances
                    x0 = (sktab%dimens(i, j) - 3)*sktab%sr(i, j)
                    x1 = x0 + sktab%sr(i, j)
                    x2 = x1 + sktab%sr(i, j)

                    ! set the spline parameters
                    xh = r - (mxind - 1)*sktab%sr(i, j)

                    do in = inu, 10
                        ! Calculate the 5th order spline extrapolation
                        f0 = sktab%skstab(in, sktab%dimens(i, j) - 2, i, j)
                        f1 = sktab%skstab(in, sktab%dimens(i, j) - 1, i, j)
                        f2 = sktab%skstab(in, sktab%dimens(i, j), i, j)
                        dd(in) = spline5th(f0, f1, f2, x0, x1, x2, xh, sktab%sr(i, j), mxind)
                    end do

                else  ! if( (ind >= sktab%dimens(i,j)) .AND. (ind < mxind) )then

                    ! Case 3:
                    !        Index+2 is outside the bounds,
                    !        AND so is Index ...
                    !
                    !        BUT, the distance is way too large.
                    !        All is zero.

                    do in = inu, 10
                        dd(in) = 0.0
                    end do

                end if ! if( (ind >= sktab%dimens(i,j)) .AND. (ind < mxind) )then

            end if ! if(ind < sktab%dimens(i,j))then

        else  ! if(ind+2 > sktab%dimens(i,j))then

            ! Case 4:
            !        Index+2 is inside the bounds
            !
            !        Do a simple interpolation.

            grdr = (r - (ind - 1.0d0)*sktab%sr(i, j))/sktab%sr(i, j)
            do in = inu, 10
                f0 = sktab%skstab(in, ind, i, j)
                f1 = sktab%skstab(in, ind + 1, i, j)
                f2 = sktab%skstab(in, ind + 2, i, j)
                dd(in) = f0 + (f1 - f0)*grdr + (f2 + f0 - 2.0*f1)*grdr*(grdr - 1.0)/2.0
            end do

        end if !  if(ind+2 > sktab%dimens(i,j))then

    end if ! if(r2 < 1e-8)then

end function skspar

! INTEGER FUNCTION SKHPAR
! =======================
!
! It is actually used to get the SK parameters,
! which are returned in 'dd' (_REAL_) after the call.
!
! That should be a subroutine, not a function!!
!
! See comments on SKSPAR above. (They are pretty similar)
!
integer function skhpar(i, j, r2, dd)

    use qm2_dftb_module, only : sktab, lmax

    implicit none
    integer:: in, ind, inu, maxmax, minmax
    _REAL_ :: cubicspline, grdr, r, spline5th

    integer :: i, j
    _REAL_ :: r2, dd(13)
    integer :: mxind
    _REAL_ :: x0, x1, x2, f0, f1, f2
    _REAL_ :: xh, hl
    skhpar = 0
    maxmax = max(lmax(i), lmax(j))
    minmax = min(lmax(i), lmax(j))
    if (maxmax <= 1) then
        inu = 10
    else
        if (maxmax <= 2) then
            if (minmax <= 1) then
                inu = 9
            else
                inu = 6
            end if
        else
            if (minmax <= 1) then
                inu = 8
            else
                if (minmax <= 2) then
                    inu = 4
                else
                    inu = 1
                end if
            end if
        end if
    end if
    mxind = sktab%dimens(i, j) + (0.3/sktab%sr(i, j) - 1.0)
    r = sqrt(r2)
    ind = r/sktab%sr(i, j) + 1.0
    if (r2 < 1e-8) then
        do in = 1, 3
            dd(in + 10) = sktab%skself(in, i)
        end do
    else
        if (ind + 2 > sktab%dimens(i, j)) then
            if (ind < sktab%dimens(i, j)) then
                x0 = (sktab%dimens(i, j) - 3)*sktab%sr(i, j)
                x1 = x0 + sktab%sr(i, j)
                x2 = x1 + sktab%sr(i, j)
                xh = r - x1
                hl = x2 - x1
                do in = inu, 10
                    f0 = sktab%skhtab(in, sktab%dimens(i, j) - 2, i, j)
                    f1 = sktab%skhtab(in, sktab%dimens(i, j) - 1, i, j)
                    f2 = sktab%skhtab(in, sktab%dimens(i, j), i, j)
                    dd(in) = cubicspline(f0, f1, f2, x0, x1, xh, hl, sktab%sr(i, j))
                end do
            else
                if ((ind >= sktab%dimens(i, j)) .AND. (ind < mxind)) then
                    x0 = (sktab%dimens(i, j) - 3)*sktab%sr(i, j)
                    x1 = x0 + sktab%sr(i, j)
                    x2 = x1 + sktab%sr(i, j)
                    xh = r - (mxind - 1)*sktab%sr(i, j)
                    do in = inu, 10
                        f0 = sktab%skhtab(in, sktab%dimens(i, j) - 2, i, j)
                        f1 = sktab%skhtab(in, sktab%dimens(i, j) - 1, i, j)
                        f2 = sktab%skhtab(in, sktab%dimens(i, j), i, j)
                        dd(in) = spline5th(f0, f1, f2, x0, x1, x2, xh, sktab%sr(i, j), mxind)
                    end do
                else
                    do in = inu, 10
                        dd(in) = 0.0
                    end do
                end if
            end if
        else
            grdr = (r - (ind - 1.0)*sktab%sr(i, j))/sktab%sr(i, j)
            do in = inu, 10
                f0 = sktab%skhtab(in, ind, i, j)
                f1 = sktab%skhtab(in, ind + 1, i, j)
                f2 = sktab%skhtab(in, ind + 2, i, j)
                dd(in) = f0 + (f1 - f0)*grdr + (f2 + f0 - 2.0*f1)*grdr*(grdr - 1.0)/2.0
            end do
        end if
    end if
end function skhpar

_REAL_ function cubicspline(f0, f1, f2, x0, x1, xh, hl, dr)
    implicit none
    _REAL_ :: f0, f1, f2, x0, x1, xh, hl, dr
    _REAL_ :: f1abl, f2abl, a, b, c, d
    f2abl = (f2 + f0 - 2.0*f1)/(dr*dr)
    f1abl = (f1 - f0)/dr + 0.5*f2abl*(x1 - x0)
    a = f1
    b = f1abl
    c = f2abl/2.0
    d = (f2 - a)/(hl*hl*hl) - b/(hl*hl) - c/hl
    cubicspline = a + b*xh + c*xh*xh + d*xh*xh*xh
end function cubicspline

_REAL_ function spline5th(f0, f1, f2, x0, x1, x2, xh, dr, mxind)
    implicit none
    _REAL_ :: f0, f1, f2, x0, x1, x2, xh, dr
    integer :: mxind
    _REAL_ :: hl, f1abl, f2abl, a, b, c, d, hsp, isp, jsp
    f2abl = (f2 + f0 - 2.0*f1)/(dr*dr)
    f1abl = (f1 - f0)/dr + 0.5*f2abl*(x1 - x0)
    a = f1
    b = f1abl
    c = f2abl/2.0
    hl = x2 - x1
    d = (f2 - a)/(hl*hl*hl) - b/(hl*hl) - c/hl
    f1abl = b + 2.0*c*hl + 3.0*d*hl*hl
    f2abl = 2.0*c + 6.0*d*hl
    hl = x2 - (mxind - 1)*dr
    hsp = 10.0*f2/(hl*hl*hl) - 4.0*f1abl/(hl*hl) + f2abl/(2.0*hl)
    isp = -15.0*f2/(hl*hl*hl*hl) + 7.0*f1abl/(hl*hl*hl) - f2abl/(hl*hl)
    jsp = 6.0*f2/(hl*hl*hl*hl*hl) - 3.0*f1abl/(hl*hl*hl*hl) + f2abl/(2.0*hl* &
        hl*hl)
    hl = xh*xh*xh
    spline5th = (hsp + isp*xh + jsp*xh*xh)*hl
end function spline5th
