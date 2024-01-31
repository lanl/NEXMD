! <compile=optimized>
!-*- mode: f90; coding: iso-8859-15; -*-
#include "dprec.fh"

subroutine fermi(ndim, nel, ev, occ, efermi)

    use qm2_dftb_module, only : log_racc, fermi_str

    implicit none

    ! Passed in:
    integer, intent(in)  :: NDIM
    _REAL_, intent(in)  :: nel
    _REAL_, intent(in)  :: ev(NDIM)
    _REAL_, intent(out) :: occ(NDIM)
    _REAL_, intent(out) :: efermi

    ! Locals
    integer :: i, nef1, nef2, nup, ndown, nocc2, ndeg
    _REAL_, parameter :: degtol = 1.0e-4          ! Degeneracy tolerance
    _REAL_, parameter :: ckbol = 3.16683735319e-6 ! Boltzmann constant in H/K
    _REAL_  :: occdg, ef1, ef2, etol

    !Locals for use with telec
    integer :: istart, iend
    logical :: tzero, break_loop
    _REAL_  :: telec     ! Electronic temperature
    _REAL_  :: beta      ! beta = 1/ kT
    _REAL_  :: chleft
    _REAL_  :: ceps
    _REAL_  :: eeps
    _REAL_  :: charge
    _REAL_  :: fac

    !<--
    occ(1:ndim) = 0.0d0
    efermi = 0.0d0
    tzero = .true.

    ! telec: initially, will be equal to qmmm_nml%telec. BUT, later I want to make
    !        this 'self-adjustable' to achieve convergence in hard cases, so it
    !        will need to be stored in a separate variable.
    telec = fermi_str%telec

    ! etol defines the energy interval around the Fermi level where orbitals will
    ! be considered 'degenerate'. (Only if telec > 0)
    if (telec .gt. 5.0d0) then
        beta = 1.0d0/(ckbol*telec)
        etol = ckbol*telec*(log(beta) - log_racc)
        tzero = .false.
    else
        etol = degtol
        tzero = .true.
    end if

    ! Find the interval of orbitals to use in
    ! calculating the Fermi Energy.
    if (nel > int(nel)) then
        nef1 = int((nel + 2)/2)
        nef2 = int((nel + 2)/2)
    else
        nef1 = int((nel + 1)/2)
        nef2 = int((nel + 2)/2)
    end if

    ! Fermi Energy
    efermi = 0.5d0*(ev(nef1) + ev(nef2))

    ! Count the number of alpha (up) and beta(down) orbitals.
    ! start with n_alpha = n_beta = nef1
    nup = nef1
    ndown = nef1

    ! Number of alpha (up) electrons:
    ! Count the number of levels ABOVE 'nef1'
    ! within etol from the Fermi level.
    ! Recall that etol depends on telec.
    do while (nup < ndim)
        if (abs(ev(nup + 1) - efermi) > etol) exit
        nup = nup + 1
    end do

    ! Number of beta (down) electrons.
    ! This will actually eliminate from the
    ! count the beta electrons from levels
    ! within etol from Fermi level. In the end,
    ! ndown will be the number of doubly occupied
    ! orbitals.
    do while (ndown > 0)
        if (abs(ev(ndown) - efermi) > etol) exit
        ndown = ndown - 1
    end do

    ! ndeg: Total number of orbitals within 'etol' from Fermi level.
    !       If all orbitals are doubly occupied, ndeg = 0
    ndeg = nup - ndown
    nocc2 = ndown     ! number of doubly occupied orbitals.

    ! Fill the doubly occupied orbitals
    do i = 1, nocc2
        occ(i) = 2.0d0
    end do

    ! if ndeg =0, all orbitals are doubly occupied,
    ! and our work here is done.
    if (ndeg == 0) then
        return
    end if

    ! If there are electrons remaining to be distributed,
    ! we do it now.
    if (tzero) then
        ! telec = 0, occupy orbitals as usual.
        occdg = ndeg
        occdg = (nel - 2.0d0*nocc2)/occdg
        do i = nocc2 + 1, nocc2 + ndeg
            occ(i) = occdg
        end do
    else
        ! telec > 0, use Fermi distribution
        chleft = nel - 2.0d0*nocc2     ! Number of electrons still to be distributed.
        istart = nocc2 + 1               ! Start from the first unnocupied level
        iend = istart + ndeg - 1       ! End at the last degenerate level.

        ! The easy case
        if (ndeg .eq. 1) then
            occ(istart) = chleft
            return
        end if

        ! Now, iterate to fill the levels by bisection
        ef1 = efermi - etol - degtol
        ef2 = efermi + etol + degtol

!! This can sometimes cause an infinite loop.
!! The value of those variables should be adjustable.
!!      ceps = dacc * chleft
!!      eeps = dacc * max( abs(ef1), abs(ef2) )

        eeps = 1.0e-12
        ceps = 1.0e-12

        ! The contents of this loop must be done at least once.
        break_loop = .false.
        do while (.not. break_loop)

            efermi = 0.5d0*(ef1 + ef2)
            charge = 0.0d0
            do i = istart, iend
                occ(i) = 2.0d0/(1.0d0 + exp(beta*(ev(i) - efermi)))
                charge = charge + occ(i)
            end do
            if (charge .gt. chleft) then
                ef2 = efermi
            else
                ef1 = efermi
            end if

            break_loop = ((abs(charge - chleft) .lt. ceps) .or. (abs(ef1 - ef2) .lt. eeps))

        end do ! do while ( .not.((abs(charge-chleft) .lt. ceps) .or. (abs(ef1-ef2) .lt. eeps) ))

        ! If accuracy was not good enough, rescale occ(i)
        if (abs(charge - chleft) .gt. ceps) then
            fac = chleft/charge
            do i = istart, iend
                occ(i) = occ(i)*fac
            end do
        end if
    end if

end subroutine fermi
