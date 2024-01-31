! <compile=optimized>

! This file has the routines that deal with the repulsive
! part of the SCC-DFTB energy. This term is determined
! by doing "pure" SCC-DFTB calcualtions on model molecules,
! then doiong full B3LYP/6-31G* calculations with the same
! molecules. The difference between the two is fit to a
! set of pairwise interactions, described by splines over
! small distances.
!
! So, the repulsive potential is stored in the form of cubic splines
! over a range of intervals. The last spline (the last
! interval) is actually a 5th order spline.
!
! If we write the spline as:
!
!         Y(t) = SUM_{j=1,n+1} a_j . t^(j-1)
!
! a cubic spline (j=3) will be:
!
!         Y(t) = a_1 + a_2 . t + a_3 . t^2 + a_4 . t^3
!
! where 't' is a parameter in the [0,1] interval.

#include "dprec.fh"

!=============================================================
! Calculates the repulsive energy as a sum over paiwise
! contributions.
!=============================================================
subroutine dftb_repulsive(nquant_nlink, izp, qm_coords, erep)

    use constants, only : A_TO_BOHRS

    implicit none

!! Passed in:
    integer, intent(in)  :: nquant_nlink
    integer, intent(in)  :: izp(*)
    _REAL_, intent(in)  :: qm_coords(3, nquant_nlink)
    _REAL_, intent(out) :: erep

!! Locals:
    integer :: i, j, k, izpj, izpk
    _REAL_ :: dif(3), r, r2

!! Functions
    _REAL_ :: dftb_repen

    do j = 2, nquant_nlink
        izpj = izp(j)
        do k = 1, j - 1
            izpk = izp(k)
            r2 = 0.0d0
            do i = 1, 3
                dif(i) = qm_coords(i, k) - qm_coords(i, j)
            end do
            r2 = dif(1)**2 + dif(2)**2 + dif(3)**2
            r = sqrt(r2)*A_TO_BOHRS
            erep = erep + dftb_repen(r, izpj, izpk)
        end do ! k
    end do ! j
end subroutine dftb_repulsive

!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!

subroutine dftb_repulsivegrd(nquant_nlink, izp, qm_coords, grd)
!======================================================================
! Calculates the gradient of the repulsive energy,
! again over pairwise contributions.
!======================================================================

    use constants, only : A_TO_BOHRS
    implicit none

!! Passed in:
    integer, intent(in)  :: nquant_nlink
    integer, intent(in)  :: izp(*)
    _REAL_, intent(in)  :: qm_coords(3, nquant_nlink)
    _REAL_, intent(out) :: grd(3, nquant_nlink)

!! Locals
    _REAL_  :: dif(3), dgr, grdr, r, r2, oner
    _REAL_  :: dftb_grdrep
    integer :: i, j, k, izpj, izpk

    do j = 2, nquant_nlink
        izpj = izp(j)
        do k = 1, j - 1
            izpk = izp(k)
            do i = 1, 3
                dif(i) = (qm_coords(i, k) - qm_coords(i, j))*A_TO_BOHRS
            end do
            r2 = dif(1)**2 + dif(2)**2 + dif(3)**2
            oner = 1.0d0/sqrt(r2)
            r = oner*r2
            grdr = dftb_grdrep(r, izpj, izpk)*oner
            do i = 1, 3
                dgr = dif(i)*grdr
                grd(i, k) = grd(i, k) + dgr
                grd(i, j) = grd(i, j) - dgr
            end do
        end do
    end do
end subroutine dftb_repulsivegrd

! ===================================================
! Returs the repulsive energy between a pair of atoms
! ===================================================
_REAL_ function dftb_repen(r, izpj, izpk)

    use qm2_dftb_module, only : spltab

    implicit none

    _REAL_  :: r
    _REAL_  :: fhelp, xh, xv1
    integer :: izpj, izpk
    integer :: i, j

    ! The splines are valid only from a minimum distance.
    ! So, a 1.0e-2 minimum is set here.

    if (r < 1.0e-2) then
        fhelp = 0.0

    else

        ! If the distance is less than the smaller distance in the SK file, then use
        ! an exponential repulsion instead, of the kind E = A + exp(B .r + C)
        !
        if (r < spltab%xr(1, 1, izpj, izpk)) then
            fhelp = exp(-spltab%efkt(1, izpj, izpk)*r + spltab%efkt(2, izpj, izpk)) &
                + spltab%efkt(3, izpj, izpk)

        else

            ! is r > cutoff? No doughnut then.
            if (r > spltab%cutoff(izpj, izpk)) then
                fhelp = 0.0

            else ! OK, r is in the interval. Let's get to work.

                ! First, find THE spline that contains 'r' in its interval
                do i = 1, spltab%numint(izpj, izpk)
                    if (r >= spltab%xr(1, i, izpj, izpk) .AND. r <= spltab%xr(2, i, izpj, izpk)) exit
                end do

                ! Now, calculate the energy from that spline.
                xv1 = r - spltab%xr(1, i, izpj, izpk)      ! 't' must be between 0 and 1
                fhelp = spltab%coeff(1, i, izpj, izpk)   ! First coefficient doesn't depend on 't'
                xh = xv1
                !calculate the polynomial:
                if (i < spltab%numint(izpj, izpk)) then
                    ! Cubic spline
                    do j = 2, 4
                        fhelp = fhelp + spltab%coeff(j, i, izpj, izpk)*xh
                        xh = xh*xv1
                    end do
                else
                    ! 5th order (last) spline
                    do j = 2, 6
                        fhelp = fhelp + spltab%coeff(j, i, izpj, izpk)*xh
                        xh = xh*xv1
                    end do
                end if
            end if
        end if
    end if
    dftb_repen = fhelp
end function dftb_repen

!==================================================
! Calculates the gradient of the repulsive energy
! for a specific atom pair (types izpj and izpk)
! separated by a distance 'r'
!==================================================
_REAL_ function dftb_grdrep(r, izpj, izpk)

    use qm2_dftb_module, only : spltab

    implicit none

    _REAL_  :: r
    _REAL_  :: grdr, xv1, xh
    integer :: izpj, izpk
    integer :: i, j

    grdr = 0.0

    if (r < spltab%xr(1, 1, izpj, izpk)) then
        ! If r < the smallest interval, use an exponential repulsion
        grdr = -spltab%efkt(1, izpj, izpk)*exp(-spltab%efkt(1, izpj, izpk)*r + spltab%efkt(2, izpj, izpk))
    else
        ! if r > the maximum interval (cutoff) gr=0
        if (r > spltab%cutoff(izpj, izpk)) then
            grdr = 0.0
        else
            ! If r is in the interval described by the splines...
            !
            ! Search only THE ONE spline with min <= r <= max.
            do i = 1, spltab%numint(izpj, izpk)
                ! Found it. Brake the loop and calculate the derivative.
                if (r >= spltab%xr(1, i, izpj, izpk) .AND. r <= spltab%xr(2, i, izpj, izpk)) exit
            end do

            !
            ! Now that we have located THE spline that includes 'r' in its interval,
            ! we calculate the derivative of the spline in that interval.
            !
            ! The derivative is:
            !
            !  Y'(t) = SUM_{j=2,n} (j-1) a_j . t^(j-2)
            !
            xv1 = r - spltab%xr(1, i, izpj, izpk)
            xh = 1
            if (i < spltab%numint(izpj, izpk)) then
                ! 'i' is one of the cubic splines
                do j = 2, 4
                    grdr = grdr + (j - 1)*spltab%coeff(j, i, izpj, izpk)*xh
                    xh = xh*xv1
                end do
            else
                ! 'i' is the final (5th order) spline
                do j = 2, 6
                    grdr = grdr + (j - 1)*spltab%coeff(j, i, izpj, izpk)*xh
                    xh = xh*xv1
                end do
            end if
        end if
    end if
    dftb_grdrep = grdr
end function dftb_grdrep
