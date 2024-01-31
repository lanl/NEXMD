! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"

!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine recip_reg here]
subroutine recip_reg(numatoms, charge, eer, vir, &
    mlimit, volume, recip, force, &
    ewaldcof, maxexp, iproc, nproc)

    use nblist, only : fraction
    use stack
    implicit none
    character(kind=1, len=9) :: routine = "recip_reg"
#  include "flocntrl.h"

    integer numatoms, mlimit(3), iproc, nproc
    _REAL_ charge(*), eer, vir(3, 3), &
        volume, recip(3, 3)
    _REAL_ force(3, *), ewaldcof, maxexp

    integer mmax, lcosf1, lcosf2, lcosf3, lsinf1, lsinf2, lsinf3, &
        lc, ls, lc12, ls12

    ! FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)
    if (do_rec == 0) return

    mmax = max(mlimit(1), mlimit(2), mlimit(3)) + 1
    call get_stack(lcosf1, mmax*numatoms, routine)
    call get_stack(lcosf2, mmax*numatoms, routine)
    call get_stack(lcosf3, mmax*numatoms, routine)
    call get_stack(lsinf1, mmax*numatoms, routine)
    call get_stack(lsinf2, mmax*numatoms, routine)
    call get_stack(lsinf3, mmax*numatoms, routine)
    call get_stack(lc, numatoms, routine)
    call get_stack(ls, numatoms, routine)
    call get_stack(lc12, numatoms, routine)
    call get_stack(ls12, numatoms, routine)
    if (.not. rstack_ok) then
        deallocate (r_stack)
        allocate (r_stack(1:lastrst), stat=alloc_ier)
        call reassign_rstack(routine)
    end if
    REQUIRE(rstack_ok)
    call do_recip_reg(numatoms, charge, eer, vir, &
        mlimit, volume, recip, force, fraction, ewaldcof, maxexp, &
        r_stack(lcosf1), r_stack(lcosf2), r_stack(lcosf3), &
        r_stack(lsinf1), r_stack(lsinf2), r_stack(lsinf3), &
        r_stack(lc12), r_stack(ls12), r_stack(lc), r_stack(ls), &
        mmax, nproc, iproc)
    call free_stack(ls12, routine)
    call free_stack(lc12, routine)
    call free_stack(ls, routine)
    call free_stack(lc, routine)
    call free_stack(lsinf3, routine)
    call free_stack(lsinf2, routine)
    call free_stack(lsinf1, routine)
    call free_stack(lcosf3, routine)
    call free_stack(lcosf2, routine)
    call free_stack(lcosf1, routine)
    return
end subroutine recip_reg
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine do_recip_reg here]
subroutine do_recip_reg(numatoms, charge, eer, vir, &
    mlimit, volume, recip, force, fraction, ewaldcof, maxexp, &
    cosf1, cosf2, cosf3, sinf1, sinf2, sinf3, c12, s12, c, s, &
    mmax, nproc, iproc)
    use constants, only : PI2, TWOPI, pi
    implicit none
#  include "def_time.h"

    integer numatoms, mlimit(3), mmax, iproc, nproc
    _REAL_ charge(*), eer, vir(3, 3), &
        volume, recip(3, 3)
    _REAL_ force(3, *), fraction(3, *), ewaldcof, maxexp
    _REAL_ cosf1(numatoms, mmax), cosf2(numatoms, mmax), &
        cosf3(numatoms, mmax), sinf1(numatoms, mmax), &
        sinf2(numatoms, mmax), sinf3(numatoms, mmax)
    _REAL_ c(numatoms), s(numatoms), &
        c12(numatoms), s12(numatoms)

    _REAL_ fac
    _REAL_ mhat1, mhat2, mhat3, denom, eterm, vterm, msq, maxexp2
    _REAL_ cstruct, sstruct, struc2, ene, mult, term
    integer n, m, m1, m2, m3, mp1, mp2, mp3, count

    call timer_start(TIME_EWFACTORS)
    fac = PI2/ewaldcof**2
    maxexp2 = maxexp**2
    ene = 0.d0
    do m2 = 1, 3
        do m1 = 1, 3
            vir(m1, m2) = 0.d0
        end do
    end do
    ! build the exponential factors for use in structure factors
    ! use sines,cosines since complex function is nonstandard
    do n = 1, numatoms
        cosf1(n, 1) = 1.d0
        cosf2(n, 1) = 1.d0
        cosf3(n, 1) = 1.d0
        sinf1(n, 1) = 0.d0
        sinf2(n, 1) = 0.d0
        sinf3(n, 1) = 0.d0
        cosf1(n, 2) = cos(twopi*fraction(1, n))
        cosf2(n, 2) = cos(twopi*fraction(2, n))
        cosf3(n, 2) = cos(twopi*fraction(3, n))
        sinf1(n, 2) = sin(twopi*fraction(1, n))
        sinf2(n, 2) = sin(twopi*fraction(2, n))
        sinf3(n, 2) = sin(twopi*fraction(3, n))
    end do
    ! get the higher factors by recursion, using trig addition rules
    ! negative values of m by complex conjugation, or even cosf, odd sinf
    do m = 3, mmax
        do n = 1, numatoms
            cosf1(n, m) = cosf1(n, m - 1)*cosf1(n, 2) - sinf1(n, m - 1)*sinf1(n, 2)
            cosf2(n, m) = cosf2(n, m - 1)*cosf2(n, 2) - sinf2(n, m - 1)*sinf2(n, 2)
            cosf3(n, m) = cosf3(n, m - 1)*cosf3(n, 2) - sinf3(n, m - 1)*sinf3(n, 2)
            sinf1(n, m) = sinf1(n, m - 1)*cosf1(n, 2) + cosf1(n, m - 1)*sinf1(n, 2)
            sinf2(n, m) = sinf2(n, m - 1)*cosf2(n, 2) + cosf2(n, m - 1)*sinf2(n, 2)
            sinf3(n, m) = sinf3(n, m - 1)*cosf3(n, 2) + cosf3(n, m - 1)*sinf3(n, 2)
        end do
    end do
    call timer_stop_start(TIME_EWFACTORS, TIME_EWRECSUM)
    count = -1
    do m1 = 0, mlimit(1)
        if (m1 == 0) mult = 1.d0
        if (m1 > 0) mult = 2.d0
        do m2 = -mlimit(2), mlimit(2)
            count = count + 1
            if (iproc /= mod(count, nproc)) goto 999
            mp1 = m1 + 1
            mp2 = abs(m2) + 1
            if (m2 < 0) then
                do n = 1, numatoms
                    c12(n) = cosf1(n, mp1)*cosf2(n, mp2) + sinf1(n, mp1)*sinf2(n, mp2)
                    s12(n) = sinf1(n, mp1)*cosf2(n, mp2) - cosf1(n, mp1)*sinf2(n, mp2)
                end do
            else
                do n = 1, numatoms
                    c12(n) = cosf1(n, mp1)*cosf2(n, mp2) - sinf1(n, mp1)*sinf2(n, mp2)
                    s12(n) = sinf1(n, mp1)*cosf2(n, mp2) + cosf1(n, mp1)*sinf2(n, mp2)
                end do
            end if
            do m3 = -mlimit(3), mlimit(3)
                ! columns of recip are reciprocal unit cell vecs
                ! thus mhat1,2,3 are cartesian components of reciprocal vector m
                mhat1 = recip(1, 1)*m1 + recip(1, 2)*m2 + recip(1, 3)*m3
                mhat2 = recip(2, 1)*m1 + recip(2, 2)*m2 + recip(2, 3)*m3
                mhat3 = recip(3, 1)*m1 + recip(3, 2)*m2 + recip(3, 3)*m3
                msq = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3
                denom = pi*volume*msq
                eterm = 0.d0
                vterm = 0.d0
                if (m1**2 + m2**2 + m3**2 > 0) then
                    eterm = exp(-fac*msq)/denom
                    vterm = 2.d0*(fac*msq + 1.d0)/msq
                end if
                ! mult takes care to double count for symmetry
                ! can take care of with eterm
                eterm = mult*eterm
                if (msq < maxexp2) then
                    mp3 = abs(m3) + 1
                    ! get the product of complex exponentials
                    if (m3 < 0) then
                        do n = 1, numatoms
                            c(n) = c12(n)*cosf3(n, mp3) + s12(n)*sinf3(n, mp3)
                            s(n) = s12(n)*cosf3(n, mp3) - c12(n)*sinf3(n, mp3)
                        end do
                    else
                        do n = 1, numatoms
                            c(n) = c12(n)*cosf3(n, mp3) - s12(n)*sinf3(n, mp3)
                            s(n) = s12(n)*cosf3(n, mp3) + c12(n)*sinf3(n, mp3)
                        end do
                    end if
                    ! get the structure factor
                    cstruct = 0.d0
                    sstruct = 0.d0
                    do n = 1, numatoms
                        cstruct = cstruct + charge(n)*c(n)
                        sstruct = sstruct + charge(n)*s(n)
                    end do
                    struc2 = cstruct**2 + sstruct**2
                    ene = ene + eterm*struc2
                    vir(1, 1) = vir(1, 1) + &
                        eterm*struc2*(vterm*mhat1*mhat1 - 1.d0)
                    vir(1, 2) = vir(1, 2) + eterm*struc2*(vterm*mhat1*mhat2)
                    vir(1, 3) = vir(1, 3) + eterm*struc2*(vterm*mhat1*mhat3)
                    vir(2, 1) = vir(2, 1) + eterm*struc2*(vterm*mhat2*mhat1)
                    vir(2, 2) = vir(2, 2) + &
                        eterm*struc2*(vterm*mhat2*mhat2 - 1.d0)
                    vir(2, 3) = vir(2, 3) + eterm*struc2*(vterm*mhat2*mhat3)
                    vir(3, 1) = vir(3, 1) + eterm*struc2*(vterm*mhat3*mhat1)
                    vir(3, 2) = vir(3, 2) + eterm*struc2*(vterm*mhat3*mhat2)
                    vir(3, 3) = vir(3, 3) + &
                        eterm*struc2*(vterm*mhat3*mhat3 - 1.d0)
                    ! get the forces
                    do n = 1, numatoms
                        term = eterm*twopi*charge(n)* &
                            (s(n)*cstruct - c(n)*sstruct)
                        force(1, n) = force(1, n) + term*mhat1
                        force(2, n) = force(2, n) + term*mhat2
                        force(3, n) = force(3, n) + term*mhat3
                    end do
                end if
            end do
999         continue
        end do
    end do
    call timer_stop(TIME_EWRECSUM)
    eer = 0.5d0*ene
    do m2 = 1, 3
        do m1 = 1, 3
            vir(m1, m2) = 0.5d0*vir(m1, m2)
        end do
    end do
    return
end subroutine do_recip_reg
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine recip_reg_dipole here]
subroutine recip_reg_dipole(numatoms, charge, eer, vir, &
    mlimit, volume, recip, force, fraction, ewaldcof, maxexp, &
    dipole, efield, &
    iproc, nproc, time, iscsum, stack_time)
    use stack
    implicit none
    character(kind=1, len=16) :: routine = "recip_reg_dipole"
#  include "flocntrl.h"

    integer numatoms, mlimit(3), iproc, nproc
    integer iscsum, stack_time
    _REAL_ charge(*), eer, vir(3, 3), &
        volume, recip(3, 3)
    _REAL_ force(3, *), fraction(3, *), &
        ewaldcof, maxexp, time(*), &
        efield(3, *), dipole(3, *)

    integer mmax, lcosf1, lcosf2, lcosf3, lsinf1, lsinf2, lsinf3, &
        lc, ls, lc12, ls12

    ! FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)
    if (do_rec == 0) return
    mmax = max(mlimit(1), mlimit(2), mlimit(3)) + 1
    call get_stack(lcosf1, mmax*numatoms, routine)
    call get_stack(lcosf2, mmax*numatoms, routine)
    call get_stack(lcosf3, mmax*numatoms, routine)
    call get_stack(lsinf1, mmax*numatoms, routine)
    call get_stack(lsinf2, mmax*numatoms, routine)
    call get_stack(lsinf3, mmax*numatoms, routine)
    call get_stack(lc, numatoms, routine)
    call get_stack(ls, numatoms, routine)
    call get_stack(lc12, numatoms, routine)
    call get_stack(ls12, numatoms, routine)
    if (.not. rstack_ok) then
        deallocate (r_stack)
        allocate (r_stack(1:lastrst), stat=alloc_ier)
        call reassign_rstack(routine)
    end if
    REQUIRE(rstack_ok)
    call do_recip_reg_dipole(numatoms, charge, eer, vir, &
        mlimit, volume, recip, force, fraction, ewaldcof, maxexp, &
        r_stack(lcosf1), r_stack(lcosf2), r_stack(lcosf3), &
        r_stack(lsinf1), r_stack(lsinf2), &
        r_stack(lsinf3), r_stack(lc12), r_stack(ls12), r_stack(lc), &
        r_stack(ls), mmax, iproc, nproc, &
        dipole, efield)
    call free_stack(ls12, routine)
    call free_stack(lc12, routine)
    call free_stack(ls, routine)
    call free_stack(lc, routine)
    call free_stack(lsinf3, routine)
    call free_stack(lsinf2, routine)
    call free_stack(lsinf1, routine)
    call free_stack(lcosf3, routine)
    call free_stack(lcosf2, routine)
    call free_stack(lcosf1, routine)
    return
end subroutine recip_reg_dipole
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine do_recip_reg_dipole here]
subroutine do_recip_reg_dipole(numatoms, charge, eer, vir, &
    mlimit, volume, recip, force, fraction, ewaldcof, maxexp, &
    cosf1, cosf2, cosf3, sinf1, sinf2, sinf3, c12, s12, c, s, mmax, &
    iproc, nproc, dipole, efield)
    use constants, only : twopi, pi2, pi
    implicit none
#  include "def_time.h"

    integer numatoms, mlimit(3), mmax, iproc, nproc
    _REAL_ charge(*), eer, vir(3, 3), &
        volume, recip(3, 3)
    _REAL_ force(3, *), fraction(3, *), ewaldcof, maxexp
    _REAL_ efield(3, *), dipole(3, *)
    _REAL_ cosf1(numatoms, mmax), cosf2(numatoms, mmax), &
        cosf3(numatoms, mmax), sinf1(numatoms, mmax), &
        sinf2(numatoms, mmax), sinf3(numatoms, mmax)
    _REAL_ c(numatoms), s(numatoms), &
        c12(numatoms), s12(numatoms)

    _REAL_ fac, vir1, vir2, vir3, vir4, vir5, vir6
    _REAL_ mhat1, mhat2, mhat3, denom, eterm, vterm, msq, maxexp2
    _REAL_ cstruct, sstruct, struc2, ene, mult
    _REAL_ m_dot_d, termr, termi, term2i
    integer n, m, m1, m2, m3, mp1, mp2, mp3, count

    call timer_start(TIME_EWFACTORS)
    fac = pi2/ewaldcof**2
    maxexp2 = maxexp**2
    ene = 0.d0
    do m2 = 1, 3
        do m1 = 1, 3
            vir(m1, m2) = 0.d0
        end do
    end do
    ! build the exponential factors for use in structure factors
    ! use sines,cosines since complex function is nonstandard
    do n = 1, numatoms
        cosf1(n, 1) = 1.d0
        cosf2(n, 1) = 1.d0
        cosf3(n, 1) = 1.d0
        sinf1(n, 1) = 0.d0
        sinf2(n, 1) = 0.d0
        sinf3(n, 1) = 0.d0
        cosf1(n, 2) = cos(twopi*fraction(1, n))
        cosf2(n, 2) = cos(twopi*fraction(2, n))
        cosf3(n, 2) = cos(twopi*fraction(3, n))
        sinf1(n, 2) = sin(twopi*fraction(1, n))
        sinf2(n, 2) = sin(twopi*fraction(2, n))
        sinf3(n, 2) = sin(twopi*fraction(3, n))
    end do
    ! get the higher factors by recursion, using trig addition rules
    ! negative values of m by complex conjugation, or even cosf, odd sinf
    do m = 3, mmax
        do n = 1, numatoms
            cosf1(n, m) = cosf1(n, m - 1)*cosf1(n, 2) - sinf1(n, m - 1)*sinf1(n, 2)
            cosf2(n, m) = cosf2(n, m - 1)*cosf2(n, 2) - sinf2(n, m - 1)*sinf2(n, 2)
            cosf3(n, m) = cosf3(n, m - 1)*cosf3(n, 2) - sinf3(n, m - 1)*sinf3(n, 2)
            sinf1(n, m) = sinf1(n, m - 1)*cosf1(n, 2) + cosf1(n, m - 1)*sinf1(n, 2)
            sinf2(n, m) = sinf2(n, m - 1)*cosf2(n, 2) + cosf2(n, m - 1)*sinf2(n, 2)
            sinf3(n, m) = sinf3(n, m - 1)*cosf3(n, 2) + cosf3(n, m - 1)*sinf3(n, 2)
        end do
    end do
    call timer_stop_start(TIME_EWFACTORS, TIME_EWRECSUM)
    count = -1
    do m1 = 0, mlimit(1)
        if (m1 == 0) mult = 1.d0
        if (m1 > 0) mult = 2.d0
        do m2 = -mlimit(2), mlimit(2)
            count = count + 1
            if (iproc /= mod(count, nproc)) goto 999
            mp1 = m1 + 1
            mp2 = abs(m2) + 1
            if (m2 < 0) then
                do n = 1, numatoms
                    c12(n) = cosf1(n, mp1)*cosf2(n, mp2) + sinf1(n, mp1)*sinf2(n, mp2)
                    s12(n) = sinf1(n, mp1)*cosf2(n, mp2) - cosf1(n, mp1)*sinf2(n, mp2)
                end do
            else
                do n = 1, numatoms
                    c12(n) = cosf1(n, mp1)*cosf2(n, mp2) - sinf1(n, mp1)*sinf2(n, mp2)
                    s12(n) = sinf1(n, mp1)*cosf2(n, mp2) + cosf1(n, mp1)*sinf2(n, mp2)
                end do
            end if
            do m3 = -mlimit(3), mlimit(3)
                ! columns of recip are reciprocal unit cell vecs
                ! thus mhat1,2,3 are cartesian components of reciprocal vector m
                mhat1 = recip(1, 1)*m1 + recip(1, 2)*m2 + recip(1, 3)*m3
                mhat2 = recip(2, 1)*m1 + recip(2, 2)*m2 + recip(2, 3)*m3
                mhat3 = recip(3, 1)*m1 + recip(3, 2)*m2 + recip(3, 3)*m3
                msq = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3
                denom = pi*volume*msq
                eterm = 0.d0
                vterm = 0.d0
                if (m1**2 + m2**2 + m3**2 > 0) then
                    eterm = exp(-fac*msq)/denom
                    vterm = 2.d0*(fac*msq + 1.d0)/msq
                end if
                ! mult takes care to double count for symmetry
                ! can take care of with eterm
                eterm = mult*eterm
                if (msq < maxexp2) then
                    mp3 = abs(m3) + 1
                    ! get the product of complex exponentials
                    if (m3 < 0) then
                        do n = 1, numatoms
                            c(n) = c12(n)*cosf3(n, mp3) + s12(n)*sinf3(n, mp3)
                            s(n) = s12(n)*cosf3(n, mp3) - c12(n)*sinf3(n, mp3)
                        end do
                    else
                        do n = 1, numatoms
                            c(n) = c12(n)*cosf3(n, mp3) - s12(n)*sinf3(n, mp3)
                            s(n) = s12(n)*cosf3(n, mp3) + c12(n)*sinf3(n, mp3)
                        end do
                    end if
                    ! get the structure factor
                    ! this time its the sum of (charge + 2pi*i*m dot dipole) times the complex exp
                    cstruct = 0.d0
                    sstruct = 0.d0
                    do n = 1, numatoms
                        m_dot_d = mhat1*dipole(1, n) + mhat2*dipole(2, n) + &
                            mhat3*dipole(3, n)
                        cstruct = cstruct + (charge(n)*c(n) - s(n)*twopi*m_dot_d)
                        sstruct = sstruct + (charge(n)*s(n) + c(n)*twopi*m_dot_d)
                    end do
                    struc2 = cstruct**2 + sstruct**2
                    ene = ene + eterm*struc2
                    ! virial not right yet. needs a contribution from field-dipole interaction
                    ! more conveniently done elsewhere
                    vir(1, 1) = vir(1, 1) + &
                        eterm*struc2*(vterm*mhat1*mhat1 - 1.d0)
                    vir(1, 2) = vir(1, 2) + eterm*struc2*(vterm*mhat1*mhat2)
                    vir(1, 3) = vir(1, 3) + eterm*struc2*(vterm*mhat1*mhat3)
                    vir(2, 1) = vir(2, 1) + eterm*struc2*(vterm*mhat2*mhat1)
                    vir(2, 2) = vir(2, 2) + &
                        eterm*struc2*(vterm*mhat2*mhat2 - 1.d0)
                    vir(2, 3) = vir(2, 3) + eterm*struc2*(vterm*mhat2*mhat3)
                    vir(3, 1) = vir(3, 1) + eterm*struc2*(vterm*mhat3*mhat1)
                    vir(3, 2) = vir(3, 2) + eterm*struc2*(vterm*mhat3*mhat2)
                    vir(3, 3) = vir(3, 3) + &
                        eterm*struc2*(vterm*mhat3*mhat3 - 1.d0)
                    ! get the fields and forces
                    do n = 1, numatoms
                        m_dot_d = mhat1*dipole(1, n) + mhat2*dipole(2, n) + &
                            mhat3*dipole(3, n)
                        termr = c(n)*cstruct + s(n)*sstruct
                        termi = s(n)*cstruct - c(n)*sstruct
                        term2i = twopi*m_dot_d*termr + charge(n)*termi
                        efield(1, n) = efield(1, n) + eterm*twopi*mhat1*termi
                        efield(2, n) = efield(2, n) + eterm*twopi*mhat2*termi
                        efield(3, n) = efield(3, n) + eterm*twopi*mhat3*termi
                        force(1, n) = force(1, n) + eterm*twopi*mhat1*term2i
                        force(2, n) = force(2, n) + eterm*twopi*mhat2*term2i
                        force(3, n) = force(3, n) + eterm*twopi*mhat3*term2i
                    end do
                end if
            end do
999         continue
        end do
    end do
    call timer_stop(TIME_EWRECSUM)
    eer = 0.5d0*ene
    do m2 = 1, 3
        do m1 = 1, 3
            vir(m1, m2) = 0.5d0*vir(m1, m2)
        end do
    end do
    return
end subroutine do_recip_reg_dipole
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine maxexp_from_mlim here]
subroutine maxexp_from_mlim(maxexp, mlimit, reclng, recip)
    implicit none
    _REAL_ maxexp
    integer mlimit(3)
    _REAL_ reclng(3), recip(3, 3)

    _REAL_ z1, z2, z3
    z1 = abs(mlimit(1)*recip(1, 1))
    z2 = abs(mlimit(2)*recip(2, 2))
    z3 = abs(mlimit(3)*recip(3, 3))
    maxexp = max(z1, z2, z3)
    return
end subroutine maxexp_from_mlim
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine find_maxexp here]
subroutine find_maxexp(ewaldcof, rtol, maxexp)
    use constants, only : PI, INVSQRTPI
    implicit none
    _REAL_ ewaldcof, rtol, maxexp

    integer i, n
    _REAL_ term, x, xlo, xhi, y, erfc

    write (6, *) 'ewaldcof,rtol = ', ewaldcof, rtol
    x = 0.5d0
    i = 0
30  x = 2.d0*x
    i = i + 1
    y = pi*x/ewaldcof
    call erfcfun(y, erfc)
    term = 2.d0*ewaldcof*erfc*INVSQRTPI
    if (term >= rtol) goto 30
    ! binary search tolerance is 2 to the -60th
    n = i + 60
    xlo = 0.d0
    xhi = x
    do i = 1, n
        x = (xlo + xhi)/2
        y = pi*x/ewaldcof
        call erfcfun(y, erfc)
        term = 2.d0*ewaldcof*erfc*INVSQRTPI
        if (term > rtol) then
            xlo = x
        else
            xhi = x
        end if
    end do
    maxexp = x

    return
end subroutine find_maxexp
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_mlim here]
subroutine get_mlim(maxexp, mlimit, eigmin, reclng, recip)
    implicit none
    _REAL_ maxexp
    _REAL_ eigmin
    integer mlimit(3)
    _REAL_ reclng(3), recip(3, 3)

    ! get coefficients for reciprocal space ewald sum

    integer mtop1, mtop2, mtop3, mlim1, mlim2, mlim3
    integer m1, m2, m3, nrecvecs
    _REAL_ z1, z2, z3, expo

    mtop1 = reclng(1)*maxexp/sqrt(eigmin)
    mtop2 = reclng(2)*maxexp/sqrt(eigmin)
    mtop3 = reclng(3)*maxexp/sqrt(eigmin)

    nrecvecs = 0
    mlim1 = 0
    mlim2 = 0
    mlim3 = 0
    do m1 = -mtop1, mtop1
        do m2 = -mtop2, mtop2
            do m3 = -mtop3, mtop3
                ! columns of recip are reciprocal unit cell vecs
                ! thus z1,z2,z3 are cartesian components of reciprocal vector m
                z1 = m1*recip(1, 1) + m2*recip(1, 2) + m3*recip(1, 3)
                z2 = m1*recip(2, 1) + m2*recip(2, 2) + m3*recip(2, 3)
                z3 = m1*recip(3, 1) + m2*recip(3, 2) + m3*recip(3, 3)
                expo = z1**2 + z2**2 + z3**2
                if (expo <= maxexp**2) then
                    nrecvecs = nrecvecs + 1
                    if (abs(m1) > mlim1) mlim1 = abs(m1)
                    if (abs(m2) > mlim2) mlim2 = abs(m2)
                    if (abs(m3) > mlim3) mlim3 = abs(m3)
                end if
            end do
        end do
    end do
    write (6, 66) nrecvecs
66  format(1x, 'number of reciprocal vecs = ', i10)
    write (6, 67) mlim1, mlim2, mlim3
67  format(1x, 'mlim1,2,3 = ', 3i6)
    mlimit(1) = mlim1
    mlimit(2) = mlim2
    mlimit(3) = mlim3
    return
end subroutine get_mlim
!---------------------------------------------------------------------
