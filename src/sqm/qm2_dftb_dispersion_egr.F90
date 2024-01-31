!  -*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.fh"
! we caculate gradient gr, not force!
! dE/dx = -F_x

! ============================================
! This routine calculated BOTH the energy and
! gradient due to dispersion terms.
! ============================================
subroutine dispersion_energy(nquant_nlink, qm_coords)

    use qm2_dftb_module, only : disper, dispertmp
    use constants, only : A_TO_BOHRS, AU_TO_EV

    implicit none

!! Passed in:
    integer, intent(in) :: nquant_nlink
    _REAL_, intent(in) :: qm_coords(3, nquant_nlink)

!! Locals
    _REAL_  :: xm, rr
    _REAL_  :: C1, dgr, r, r2, dif(3)
    _REAL_  :: conv
    _REAL_  :: h1, h2, h3, fdamp, g612
    integer :: i, j, n, m

!! Pointers
    _REAL_, pointer  :: A, B, C, r0, rv, Rvdw(:, :), C6(:, :)

!! This is just to make the file more readable
    A => dispertmp%A
    B => dispertmp%B
    C => dispertmp%C
    r0 => dispertmp%r0
    rv => dispertmp%rv
    Rvdw => dispertmp%Rvdw
    C6 => dispertmp%C6

    ! calculate pairwise contribution i-j
    ! all values in eV and Angstr.
    ! convert r in Angstr
    ! Rvdw is already in Angstr
    ! then convert energy Edis from eV in H

    do i = 1, nquant_nlink
        do j = 1, i - 1
            rr = C/Rvdw(i, j)**A
            h3 = rv*0.5d0*Rvdw(i, j)**6
            r2 = 0.0d0
            do n = 1, 3
                dif(n) = qm_coords(n, i) - qm_coords(n, j)
                r2 = r2 + dif(n)**2
            end do
            r = sqrt(r2)

            ! Dispersion energy
            call dis_e(r, rr, h3, C1)
            disper%Edis = disper%Edis + C1*C6(i, j)/AU_TO_EV

        end do !j
    end do !i
    return
end subroutine dispersion_energy

subroutine dispersion_grad(qmmm_struct, nquant_nlink, gr)

    use qm2_dftb_module, only : disper, dispertmp
    use constants, only : BOHRS_TO_A, AU_TO_EV
    use qmmm_struct_module, only : qmmm_struct_type

    implicit none

!! Passed in:
    type(qmmm_struct_type), intent(in) :: qmmm_struct
    integer, intent(in) :: nquant_nlink
    _REAL_, intent(out) :: gr(3, *)

!! Locals
    _REAL_  :: xm, rr
    _REAL_  :: C1, dgr, r, r2, dif(3)
    _REAL_  :: conv
    _REAL_  :: h1, h2, h3, fdamp, g612
    integer :: i, j, n, m

!! Pointers
    _REAL_, pointer  :: A, B, C, r0, rv, Rvdw(:, :), C6(:, :)

!! This is just to make the file more readable
    A => dispertmp%A
    B => dispertmp%B
    C => dispertmp%C
    r0 => dispertmp%r0
    rv => dispertmp%rv
    Rvdw => dispertmp%Rvdw
    C6 => dispertmp%C6

    ! calculate pairwise contribution i-j
    ! all values in eV and Angstr.
    ! convert r in Angstr
    ! Rvdw is already in Angstr
    ! then convert energy Edis from eV in H

    do i = 1, nquant_nlink
        do j = 1, i - 1
            rr = C/Rvdw(i, j)**A
            h3 = rv*0.5d0*Rvdw(i, j)**6
            r2 = 0.0d0
            do n = 1, 3
                dif(n) = qmmm_struct%qm_coords(n, i) - qmmm_struct%qm_coords(n, j)
                r2 = r2 + dif(n)**2
            end do
            r = sqrt(r2)

            ! Dispersion energy - JUST TO GET C1 and H1 - NEEDS IMPROVING AT SOME POINT
            call dis_e(r, rr, h3, C1)

            ! Dispersion gradients
            do m = 1, 3
                xm = dif(m)
                call dis_gr(r, rr, xm, h3, C1, dgr)
                gr(m, i) = gr(m, i) + dgr*C6(i, j)
                gr(m, j) = gr(m, j) - dgr*C6(i, j)
            end do

        end do !j
    end do !i
    return
end subroutine dispersion_grad

! ====================================
! Calculate the pair dispersion energy
! ====================================
subroutine dis_e(r, rr, h3, C1)

    use qm2_dftb_module, only : dispertmp
    use constants, only : INVPI

    implicit none

!! Passed in:
    _REAL_, intent(in) :: r
    _REAL_, intent(in) :: rr
    _REAL_, intent(in) :: h3
    _REAL_, intent(out) :: C1

!! Locals
    _REAL_ :: h1, h2, fdamp, g612
    _REAL_ :: Catan, dgr

!! Pointers
    _REAL_, pointer :: A, B, C, r0, rv

!! This is just to make the file more readable
    A => dispertmp%A
    B => dispertmp%B
    C => dispertmp%C
    r0 => dispertmp%r0
    rv => dispertmp%rv

    if (C >= 0.0d0) then
        ! old scaling
        Catan = INVPI*atan(A*(r - r0 - rv))
        C1 = -B/ &
            ((r - rv)**6 &
            + (0.5d0 - Catan)*C*(r - r0 - rv)**2)
    else
        ! C is negative!! rv switches R**12, scaling function as in JCP!
        h1 = rr*r**(A)
        h2 = exp(h1)
        fdamp = 1 - h2
        g612 = (h3/r**12) - (1.0d0/r**6)
        C1 = g612*fdamp**(B)
    end if
    return
end subroutine dis_e

!==================================================
! Calculate a pairwise contribution to the gradient
! due to dispersion forces
!==================================================
subroutine dis_gr(r, rr, xm, h3, C1, dgr)

    use qm2_dftb_module, only : dispertmp
    use constants, only : BOHRS_TO_A, AU_TO_EV, INVPI

    implicit none
!! Passed in:
    _REAL_, intent(in) :: r
    _REAL_, intent(in) :: rr
    _REAL_, intent(in) :: xm
    _REAL_, intent(in) :: h3
    _REAL_, intent(in) :: C1
    _REAL_, intent(out) :: dgr

!! Locals
    _REAL_ :: Catan, CONV
    _REAL_ :: h1, h2, fdamp, g612

!! Pointers
    _REAL_, pointer :: A, B, C, r0, rv

!! This is just to make the file more readable
!! (Not that it really changes much ;-) )
    A => dispertmp%A
    B => dispertmp%B
    C => dispertmp%C
    r0 => dispertmp%r0
    rv => dispertmp%rv

    CONV = BOHRS_TO_A/AU_TO_EV

    if (C >= 0.0d0) then
        Catan = INVPI*atan(A*(r - r0 - rv))
        dgr = CONV*B*((6.0d0*xm*(-rv + r)**5)/ &
            r - (0.3184713375796178d0*A*C*xm*(-r0 - rv + r)**2)/ &
            (r*(1.0d0 + A**2*(-r0 - rv + r)**2)) &
            + (2.0d0*C*xm*(-r0 - rv + r)*(0.5d0 - Catan))/r)/ &
            (((-rv + r)**6 + C*(-r0 - rv + r)**2*(0.5d0 - Catan))**2)
    else
        h1 = rr*r**A
        h2 = exp(h1)
        fdamp = 1.0d0 - h2
        g612 = (h3/r**12) - (1.0d0/r**6)
        dgr = CONV*(xm/r)*(g612*B*(fdamp**(B - 1.0d0))*(-h2)*A*rr*r**(A - 1.0d0) &
            + (fdamp**B)*(6.0d0/r**7 - 12.0d0*h3/r**13))
    end if
    return
end subroutine dis_gr
