! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

module relax_mat

    private caldis, calrate, corf, dinten, drates, indexn, kmat, remarc
    public noeread, noecalc

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Get distance-dependent terms for relaxation matrix analysis
#ifdef NMODE
    subroutine caldis(x, ddep, dddep, newf, amass)
#else
    subroutine caldis(x, ddep, dddep)
#endif

        !  Subroutine CALculate DIStances:

        !     Calculates distance information bewteen all hydrogen atoms:
        !           nath = number of H atoms
        !           x    = array of atom coordinates
        !           ddep  = distance dependent part of rate matrix element
        !           dddep = derivative of ddep with respect to proton m,i

        !      Treats methyl rotors as three equivalent protons, and
        !      aromatic D and E as equivalent.  The distances are
        !      calculated assuming fast motion for the methyl group
        !      and slow (compared to correlation time) motion for aromatic flip,
        !      which is the same as 1/6 averaging.

        !      Variables "i" and "j" are indexed in the nath scheme; "ii"
        !      and "jj" are the corresponding entries in the natmet scheme.
        use constants, only : zero
        implicit none
#  include "nmr.h"

!Passed in
        _REAL_ :: x(*), ddep(ma, ma), dddep(3, ma, ma)
#ifdef NMODE
        _REAL_ :: amass(*)
        logical newf
#endif

!Local
        _REAL_ :: dddepp(3)
        _REAL_ :: dmet(3, 3), rmet(3, 3, 3)
        _REAL_ :: sum2, sumeff, dprod, denom, term, diffij, work
        integer :: i, jj, ii, j, m2ij, kl, k, l, klmn, m, n, nn, lk, lkm

        ! --- constants:

        _REAL_, parameter :: x5o81 = (5.0d0/81.0d0)
        _REAL_, parameter :: x2o81 = (2.0d0/81.0d0)
        _REAL_, parameter :: x6o81 = (6.0d0/81.0d0)
        _REAL_, parameter :: x5o9 = (5.0d0/9.0d0)
        _REAL_, parameter :: x2o9 = (2.0d0/9.0d0)
        _REAL_, parameter :: x6o9 = (6.0d0/9.0d0)
        _REAL_, parameter :: x5o18 = (5.0d0/18.0d0)
        _REAL_, parameter :: x1o9 = (1.0d0/9.0d0)
        _REAL_, parameter :: x1o3 = (1.0d0/3.0d0)
        _REAL_, parameter :: x1o18 = (1.0d0/18.0d0)
        _REAL_, parameter :: x5o36 = (5.0d0/36.0d0)
        ! --- subscript arrays to collapse mutliple do-loops into one:
        integer, dimension(9), parameter :: ksub = (/1, 2, 3, 1, 2, 3, 1, 2, 3/)
        integer, dimension(9), parameter :: lsub = (/1, 1, 1, 2, 2, 2, 3, 3, 3/)
        integer, dimension(81), parameter ::  k1sub = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
            2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, &
            3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3/)
        integer, dimension(81), parameter ::  l1sub = (/1, 1, 1, 1, 1, 1, 1, 1, 1, &
            2, 2, 2, 2, 2, 2, 2, 2, 2, &
            3, 3, 3, 3, 3, 3, 3, 3, 3, &
            1, 1, 1, 1, 1, 1, 1, 1, 1, &
            2, 2, 2, 2, 2, 2, 2, 2, 2, &
            3, 3, 3, 3, 3, 3, 3, 3, 3, &
            1, 1, 1, 1, 1, 1, 1, 1, 1, &
            2, 2, 2, 2, 2, 2, 2, 2, 2, &
            3, 3, 3, 3, 3, 3, 3, 3, 3/)
        integer, dimension(81), parameter ::  m1sub = (/1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, &
            1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, &
            1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3/)
        integer, dimension(81), parameter ::  n1sub = (/1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, &
            1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, &
            1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, &
            1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, &
            1, 2, 3, 1, 2, 3, 1, 2, 3/)
        integer, dimension(6), parameter :: l2sub = (/1, 1, 1, 2, 2, 2/)
        integer, dimension(6), parameter :: k2sub = (/1, 2, 3, 1, 2, 3/)
        integer, dimension(4), parameter :: k3sub = (/1, 1, 2, 2/)
        integer, dimension(4), parameter :: l3sub = (/1, 2, 1, 2/)
        integer, dimension(18), parameter :: l4sub = (/1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2/)
        integer, dimension(18), parameter :: k4sub = (/1, 1, 1, 2, 2, 2, 3, 3, 3, 1, 1, 1, 2, 2, 2, 3, 3, 3/)
        integer, dimension(18), parameter ::  m4sub = (/1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3/)

        do i = 1, nath
            do jj = 1, natmet
                dddep(1, i, jj) = zero
                dddep(2, i, jj) = zero
                dddep(3, i, jj) = zero
            end do
        end do
        ii = 0
        do i = 1, nath

            !                **first 2 methyl protons of a methyl group
            !                   and first aromatic D or E skipped

            if (m2(i) == 1 .or. m2(i) == 4) cycle
            ii = ii + 1

            jj = 0
            do j = 1, nath

                !                  **first 2 methyl protons of a methyl group
                !                    and first aromatic D or E skipped

                if (m2(j) == 1 .or. m2(j) == 4) cycle
                jj = jj + 1

                m2ij = m2(i)*m2(j)
                if (m2ij == 9) then !             **methyl - methyl:

                    if (ii == jj) then !           **intramethyl:

                        ddep(ii, jj) = 8.13d-3  ! this is 1/(4*1.77**6)

                    else !                         **intermethyl:

                        do kl = 1, 9
                            k = ksub(kl)
                            l = lsub(kl)
                            rmet(k, l, 1) = x(3*(i - k) + 1) - x(3*(j - l) + 1)
                            sum2 = rmet(k, l, 1)**2
                            rmet(k, l, 2) = x(3*(i - k) + 2) - x(3*(j - l) + 2)
                            sum2 = sum2 + rmet(k, l, 2)**2
                            rmet(k, l, 3) = x(3*(i - k) + 3) - x(3*(j - l) + 3)
                            sum2 = sum2 + rmet(k, l, 3)**2
                            dmet(k, l) = sum2
                        end do
                        sumeff = 0.0d0
                        do klmn = 1, 81
                            k = k1sub(klmn)
                            l = l1sub(klmn)
                            m = m1sub(klmn)
                            n = n1sub(klmn)
                            dprod = (rmet(k, l, 1)*rmet(m, n, 1)) &
                                + (rmet(k, l, 2)*rmet(m, n, 2)) &
                                + (rmet(k, l, 3)*rmet(m, n, 3))
                            denom = 1.0d0/sqrt((dmet(k, l)*dmet(m, n))**5)
                            term = (3.0d0*dprod**2 - dmet(k, l)*dmet(m, n))*denom
                            sumeff = sumeff + term
                            do nn = 1, 3
                                dddep(nn, i - k + 1, jj) = dddep(nn, i - k + 1, jj) &
                                    - x5o81*term*rmet(k, l, nn)/dmet(k, l) &
                                    - (x2o81*dmet(m, n)*rmet(k, l, nn) - &
                                    x6o81*dprod*rmet(m, n, nn))*denom
                            end do
                        end do
                        ddep(ii, jj) = sumeff/162.0d0
                    end if

                else if (m2ij == 6) then !     **methyl to ordinary:

                    if (m2(i) == 3) then !       ** i is the methyl:
                        do k = 1, 3
                            rmet(k, 1, 1) = x(3*(i - k) + 1) - x(3*(j - 1) + 1)
                            sum2 = rmet(k, 1, 1)**2
                            rmet(k, 1, 2) = x(3*(i - k) + 2) - x(3*(j - 1) + 2)
                            sum2 = sum2 + rmet(k, 1, 2)**2
                            rmet(k, 1, 3) = x(3*(i - k) + 3) - x(3*(j - 1) + 3)
                            sum2 = sum2 + rmet(k, 1, 3)**2
                            dmet(k, 1) = sum2
                        end do
                        sumeff = 0.0d0
                        do kl = 1, 9
                            k = ksub(kl)
                            l = lsub(kl)
                            dprod = rmet(k, 1, 1)*rmet(l, 1, 1) &
                                + rmet(k, 1, 2)*rmet(l, 1, 2) &
                                + rmet(k, 1, 3)*rmet(l, 1, 3)
                            denom = 1.0d0/sqrt((dmet(k, 1)*dmet(l, 1))**5)
                            term = (3.0d0*dprod**2 - dmet(k, 1)*dmet(l, 1))*denom
                            sumeff = sumeff + term
                            do n = 1, 3
                                dddep(n, i - k + 1, jj) = dddep(n, i - k + 1, jj) &
                                    - x5o9*term*rmet(k, 1, n)/dmet(k, 1) &
                                    - (x2o9*dmet(l, 1)*rmet(k, 1, n) - &
                                    x6o9*dprod*rmet(l, 1, n))*denom
                            end do
                        end do
                        ddep(ii, jj) = sumeff/18.0d0

                    else !                     ** j is the methyl:
                        do k = 1, 3
                            rmet(k, 1, 1) = x(3*(i - 1) + 1) - x(3*(j - k) + 1)
                            sum2 = rmet(k, 1, 1)**2
                            rmet(k, 1, 2) = x(3*(i - 1) + 2) - x(3*(j - k) + 2)
                            sum2 = sum2 + rmet(k, 1, 2)**2
                            rmet(k, 1, 3) = x(3*(i - 1) + 3) - x(3*(j - k) + 3)
                            sum2 = sum2 + rmet(k, 1, 3)**2
                            dmet(k, 1) = sum2
                        end do
                        sumeff = 0.0d0
                        do kl = 1, 9
                            k = ksub(kl)
                            l = lsub(kl)
                            dprod = rmet(k, 1, 1)*rmet(l, 1, 1) &
                                + rmet(k, 1, 2)*rmet(l, 1, 2) &
                                + rmet(k, 1, 3)*rmet(l, 1, 3)
                            denom = 1.0d0/sqrt((dmet(k, 1)*dmet(l, 1))**5)
                            term = (3.0d0*dprod**2 - dmet(k, 1)*dmet(l, 1))*denom
                            sumeff = sumeff + term
                            do n = 1, 3
                                dddep(n, i, jj) = dddep(n, i, jj) &
                                    - x5o18*term*(rmet(k, 1, n)/dmet(k, 1) &
                                    + rmet(l, 1, n)/dmet(l, 1)) &
                                    - (dmet(l, 1)*rmet(k, 1, n) &
                                    + dmet(k, 1)*rmet(l, 1, n) &
                                    - 3.0d0*dprod*(rmet(k, 1, n) + rmet(l, 1, n)))*denom*x1o9
                            end do
                        end do
                        ddep(ii, jj) = sumeff/18.0d0

                    end if

                else if (m2ij == 15) then !        *** aromatic to methyl:

                    if (m2(i) == 3) then !           **i is methyl, j is aromatic:
                        do lk = 1, 6
                            l = l2sub(lk)
                            k = k2sub(lk)
                            rmet(k, l, 1) = x(3*(i - k) + 1) - x(3*(j - l) + 1)
                            sum2 = rmet(k, l, 1)**2
                            rmet(k, l, 2) = x(3*(i - k) + 2) - x(3*(j - l) + 2)
                            sum2 = sum2 + rmet(k, l, 2)**2
                            rmet(k, l, 3) = x(3*(i - k) + 3) - x(3*(j - l) + 3)
                            sum2 = sum2 + rmet(k, l, 3)**2
                            dmet(k, l) = sum2
                        end do
                        sumeff = 0.0d0
                        do lkm = 1, 18
                            l = l4sub(lkm)
                            k = k4sub(lkm)
                            m = m4sub(lkm)
                            dprod = (rmet(k, l, 1)*rmet(m, l, 1)) &
                                + (rmet(k, l, 2)*rmet(m, l, 2)) &
                                + (rmet(k, l, 3)*rmet(m, l, 3))
                            denom = 1.0d0/sqrt((dmet(k, l)*dmet(m, l))**5)
                            term = (3.0d0*dprod**2 - dmet(k, l)*dmet(m, l))*denom
                            sumeff = sumeff + term
                            do n = 1, 3
                                dddep(n, i - k + 1, jj) = dddep(n, i - k + 1, jj) &
                                    - x5o18*term*rmet(k, l, n)/dmet(k, l) &
                                    - (x1o9*dmet(m, l)*rmet(k, l, n) - &
                                    x1o3*dprod*rmet(m, l, n))*denom
                            end do
                        end do
                        ddep(ii, jj) = sumeff/36.0d0

                    else !                          **i is aromatic, j is methyl:

                        do lk = 1, 6
                            l = l2sub(lk)
                            k = k2sub(lk)
                            rmet(k, l, 1) = x(3*(i - l) + 1) - x(3*(j - k) + 1)
                            sum2 = rmet(k, l, 1)**2
                            rmet(k, l, 2) = x(3*(i - l) + 2) - x(3*(j - k) + 2)
                            sum2 = sum2 + rmet(k, l, 2)**2
                            rmet(k, l, 3) = x(3*(i - l) + 3) - x(3*(j - k) + 3)
                            sum2 = sum2 + rmet(k, l, 3)**2
                            dmet(k, l) = sum2
                        end do
                        sumeff = 0.0d0
                        do lkm = 1, 18
                            l = l4sub(lkm)
                            k = k4sub(lkm)
                            m = m4sub(lkm)
                            dprod = (rmet(k, l, 1)*rmet(m, l, 1)) &
                                + (rmet(k, l, 2)*rmet(m, l, 2)) &
                                + (rmet(k, l, 3)*rmet(m, l, 3))
                            denom = 1.d0/sqrt((dmet(k, l)*dmet(m, l))**5)
                            term = (3.0d0*dprod**2 - dmet(k, l)*dmet(m, l))*denom
                            sumeff = sumeff + term
                            do n = 1, 3
                                dddep(n, i - l + 1, jj) = dddep(n, i - l + 1, jj) &
                                    - x5o36*term*(rmet(k, l, n)/dmet(k, l) &
                                    + rmet(m, l, n)/dmet(m, l)) &
                                    - (dmet(m, l)*rmet(k, l, n) &
                                    + dmet(k, l)*rmet(m, l, n) &
                                    - 3.0d0*dprod*(rmet(k, l, n) + rmet(m, l, n)))*denom*x1o18
                            end do
                        end do
                        ddep(ii, jj) = sumeff/36.0d0

                    end if

                else if (m2ij == 25) then !           **aromatic to aromatic:

                    if (ii /= jj) then !                **inter
                        sumeff = 0.0d0
                        do kl = 1, 4
                            k = k3sub(kl)
                            l = l3sub(kl)
                            sum2 = 0.0d0
                            do n = 1, 3
                                diffij = x(3*(i - k) + n) - x(3*(j - l) + n)
                                sum2 = sum2 + diffij**2
                                dddepp(n) = -6.0d0*diffij
                            end do
                            sum2 = 1.0d0/sum2
                            term = sum2**3
                            sumeff = sumeff + term
                            dddep(1, i - k + 1, jj) = dddep(1, i - k + 1, jj) + &
                                0.25d0*term*dddepp(1)*sum2
                            dddep(2, i - k + 1, jj) = dddep(2, i - k + 1, jj) + &
                                0.25d0*term*dddepp(2)*sum2
                            dddep(3, i - k + 1, jj) = dddep(3, i - k + 1, jj) + &
                                0.25d0*term*dddepp(3)*sum2
                        end do
                        ddep(ii, jj) = 0.25d0*sumeff

                    else !                                **intra:

                        ddep(ii, jj) = 1.582d-4  !  this is 1/4.30**6

                    end if

                else if (m2ij == 10) then !              **aromatic to ordinary:

                    if (m2(i) /= 2) then !                **i is aromatic:
                        sumeff = 0.0d0
                        do k = 1, 2
                            sum2 = 0.0d0
                            do n = 1, 3
                                diffij = x(3*(i - k) + n) - x(3*(j - 1) + n)
                                sum2 = sum2 + diffij**2
                                dddep(n, i - k + 1, jj) = -6.0d0*diffij
                            end do
                            sum2 = 1.0d0/sum2
                            term = sum2**3
                            sumeff = sumeff + term
                            dddep(1, i - k + 1, jj) = 0.5d0*term*dddep(1, i - k + 1, jj)*sum2
                            dddep(2, i - k + 1, jj) = 0.5d0*term*dddep(2, i - k + 1, jj)*sum2
                            dddep(3, i - k + 1, jj) = 0.5d0*term*dddep(3, i - k + 1, jj)*sum2
                        end do
                        ddep(ii, jj) = 0.5d0*sumeff

                    else !                                **j is aromatic:

                        sumeff = 0.0d0
                        do k = 1, 2
                            sum2 = 0.0d0
                            do n = 1, 3
                                diffij = x(3*(i - 1) + n) - x(3*(j - k) + n)
                                sum2 = sum2 + diffij**2
                                dddepp(n) = -6.0d0*diffij
                            end do
                            sum2 = 1.0d0/sum2
                            term = sum2**3
                            sumeff = sumeff + term
                            dddep(1, i, jj) = dddep(1, i, jj) + 0.5d0*term*dddepp(1)*sum2
                            dddep(2, i, jj) = dddep(2, i, jj) + 0.5d0*term*dddepp(2)*sum2
                            dddep(3, i, jj) = dddep(3, i, jj) + 0.5d0*term*dddepp(3)*sum2
                        end do
                        ddep(ii, jj) = 0.5d0*sumeff
                    end if

                else if (m2(i) == 2 .and. m2(j) == 2) then !  *** both i and j are ordinary:

                    if (i == j) then !                   **intra:

                        ddep(ii, jj) = 0.0d0
                    else if (i < j) then !               **inter,first time:

                        sum2 = 0.0d0
                        do n = 1, 3
                            diffij = x(3*(i - 1) + n) - x(3*(j - 1) + n)
                            sum2 = sum2 + diffij**2
                            dddep(n, i, jj) = -6.0d0*diffij
                        end do
                        work = 1.0d0/sum2**4
                        dddep(1, i, jj) = dddep(1, i, jj)*work
                        dddep(2, i, jj) = dddep(2, i, jj)*work
                        dddep(3, i, jj) = dddep(3, i, jj)*work
                        ddep(ii, jj) = work*sum2
#ifdef NMODE

                        !  --- get correction factor for this spin pair:

                        call corf(x, i, j, ihyp(i), ihyp(j), newf, amass)
                        newf = .false.
                        ddep(ii, jj) = ddep(ii, jj)*gamma_nmr
                        dddep(1, i, jj) = gamma_nmr*dddep(1, i, jj) + &
                            ddep(ii, jj)*dgamma(1)
                        dddep(2, i, jj) = gamma_nmr*dddep(2, i, jj) + &
                            ddep(ii, jj)*dgamma(2)
                        dddep(3, i, jj) = gamma_nmr*dddep(3, i, jj) + &
                            ddep(ii, jj)*dgamma(3)
                        do n = 1, iscale
                            dratg(ii, jj, n) = ddep(ii, jj)*dgamma(n + 3)/gamma_nmr
                        end do
#endif
                    else !                               **inter, use previous:

                        dddep(1, i, jj) = -dddep(1, j, ii)
                        dddep(2, i, jj) = -dddep(2, j, ii)
                        dddep(3, i, jj) = -dddep(3, j, ii)
                        ddep(ii, jj) = ddep(jj, ii)
#ifdef NMODE
                        do n = 1, iscale
                            dratg(ii, jj, n) = dratg(jj, ii, n)
                        end do
#endif
                    end if

                else
                    write (6, *) 'caldis error:', i, j, m2(i), m2(j)
                    call mexit(6, 1)
                end if

            end do
        end do

        return
    end subroutine caldis

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine calrate here]
    subroutine calrate(ddep, rate, trp)

        !   Subroutine CALculate RATE matrix.

        !   Correlation times are calculated according to:
        !          tau(i,j) = tau(i)*tau(j)/(tau(i)+tau(j))
        !   That is they assume that the motion of the whole molecule, as well
        !   as internal motion is isotropic.

        !   The spectral densities default to a resonant frequency of 500
        !   MHz (see OMEGA in subroutine mdread).
        !   Units of time for correlation times are nsec;
        !   units for the rate matrix elements are sec**-1.

        !-----modified P F Yip 8/9/89, by dac 10/25/89, Ross Walker 08/02/05
        use constants, only : TWOPI
        implicit none
#  include "nmr.h"

        _REAL_ :: gsh, gyr, hbar, gsj, tenm16, con0, con1, con2
        parameter(gyr=2.6751965d4, hbar=1.0545919d0, gsh=gyr*hbar*gyr)
        parameter(tenm16=1.d-16, con0=gsh*gsh*tenm16)
        parameter(con1=3.d0*con0, con2=6*con0)

        _REAL_ :: rate(ma, ma), trp(ma, ma)
        _REAL_ :: spd(lst, 3)
        _REAL_ :: ddep(ma, ma)

!Local
        _REAL_ :: omega2, taui, tauj, tauw, tsws, ts2ws
        _REAL_ :: popisq, popjsq, popij, s0, s1, s2, ootc3, tc3
        integer :: ij, i, j, k, ii, jj

        !  --- square of Larmor frequency in (radians/nsec)**2 :

        omega2 = (TWOPI*omega*1.d-3)**2

        !---- the first block ( the 80 do loop ) calculates the spec densities

        ij = 0
        do i = 1, nath
            if (m2(i) == 1 .or. m2(i) == 4) cycle
            taui = tau(i)
            do j = 1, i
                if (m2(j) == 1 .or. m2(j) == 4) cycle
                ij = ij + 1
                tauj = tau(j)
                tauw = (taui*tauj)/(taui + tauj)
                tsws = tauw*tauw*omega2
                ts2ws = 4.d0*tsws
                spd(ij, 1) = con0*tauw
                spd(ij, 2) = con1*tauw/(1.d0 + tsws)
                spd(ij, 3) = con2*tauw/(1.d0 + ts2ws)
            end do
        end do

        ! Calculate transition probabilites from distances and spectral densities.

        ! Note: populations are actually sqrt(pop(i)) and the spectral densities
        ! are calculated to generate the symmetrized rate matrix:

        !                                -1
        !                     R_prime = P   R  P

        ! where P is a diagonal matrix containing the square roots of the
        ! populations of each group of equivalent nuclei.

        !      The fundamental equations are from Macura and Ernst, Mol. Phys.
        !       41, 95-117 (1980), see esp. section 4.  The symmetrization
        !       procedure has been derived by many people; an explicit
        !       presentation that corresponds to the way we do things (i.e.,
        !       with jump models for methyl groups) is given by Olejniczak,
        !       J. Magn. Res. 81, 392-394 (1989).

        !       The rate matrix R ( nonsymmetric )  is:

        !         rate(i,i)=selfrel(i,i)+{SUM(j not i )otherrel(i,j)}
        !         rate(i,j)=crossrel(i,j) = population(i)*(W2(i,j)-W0(i,j))
        !       [Note that population(i) is the same as popisq and also
        !       ( pop(i)**2 ).]

        !        After the transformation above,
        !              rate(i,j)=pop(i)*pop(j)*(W2(i,j)-W0(i,j))
        !        This is computed in the (else loop of if(ii.eq.jj))

        !        The transformed rate matrix is now symmetric and can be
        !        diagonalized.  Later, in subroutine remarc, when the
        !        intensities are calcualated, they are multiplied by
        !        population factors that effect the inverse transformation
        !        back to the original magnetization Bloch equations.

        k = 1
        ii = 0
        do i = 1, nath
            if (m2(i) == 1 .or. m2(i) == 4) cycle
            ii = ii + 1
            if (m2(i) == 3) then
                popisq = 3.0d0
            else if (m2(i) == 5) then
                popisq = 2.0d0
            else
                popisq = 1.0d0
            end if
            jj = 0
            do j = 1, i
                if (m2(j) == 1 .or. m2(j) == 4) cycle
                jj = jj + 1
                if (m2(j) == 3) then
                    popjsq = 3.0d0
                else if (m2(j) == 5) then
                    popjsq = 2.0d0
                else
                    popjsq = 1.0d0
                end if
                popij = pop(i)*pop(j)

                s0 = spd(k, 1)
                s1 = spd(k, 2)
                s2 = spd(k, 3)
                if (ii == jj) then

                    !       For the equivalent protons in the group self-relaxtion must
                    !       be accounted for: (W1 + W2)

                    if (popisq == 1.0d0) then
                        rate(ii, ii) = 0.0d0
                    else if (popisq == 2.0d0) then
                        if (iroesy > 0) then
                            rate(ii, ii) = 0.5d0*(9.d0*s0 + 3.d0*s1 + 2.*s2)*ddep(ii, ii)
                        else
                            rate(ii, ii) = 2.0d0*(0.50d0*s1 + s2)*ddep(ii, ii)
                        end if
                    else if (popisq == 3.0d0) then
                        if (iroesy > 0) then

                            ! --- just use the fast methyl motion limit here:

                            rate(ii, ii) = 0.5d0*(13.d0*s0 + 3.d0*s1 + &
                                3.d0*s2)*ddep(ii, ii)
                        else

                            !      For methyl relaxation,
                            !  --- use Kalk & Berendsen formula [JMR 24, 343 (1978)], Eq. 11:

                            taui = 0.5d0*tau(i)
                            ootc3 = 1./taui + 1./taumet
                            tc3 = 1./ootc3
                            rate(ii, ii) = 2.0d0*5.372d0*(0.25d0*( &
                                taui/(1.+omega2*taui*taui) &
                                + 4.d0*taui/(1.+omega2*4.d0*taui*taui)) &
                                + 0.75d0*(tc3/(1.d0 + omega2*tc3*tc3) &
                                + 4.d0*tc3/(1.d0 + omega2*4.d0*tc3*tc3)))
#if 0
                            write (6, *) 'for methyl: ', ii, rate(ii, ii), taui, tc3
                            term1 = 2.0d0*5.372d0*0.25d0* &
                                taui/(1.+omega2*taui*taui)
                            term2 = 2.0d0*5.372d0*0.25d0* &
                                4.d0*taui/(1.+omega2*4.d0*taui*taui)
                            term3 = 2.0d0*5.372d0* &
                                0.75d0*tc3/(1.d0 + omega2*tc3*tc3)
                            term4 = 2.0d0*5.372d0* &
                                0.75d0*4.d0*tc3/(1.d0 + omega2*4.d0*tc3*tc3)
                            write (6, *) '            ', term1, term2, term3, term4
#endif
                        end if
                    else
                        write (6, *) 'bad value for popisq:', popisq, i, ii
                        call mexit(6, 1)
                    end if  ! (popisq == 1.0d0)
                    trp(ii, ii) = 0.0d0
                else

                    !         Calculate the cross relaxation terms in rate: (W2 - W0)

                    if (iroesy > 0) then
                        rate(ii, jj) = ddep(ii, jj)*(s2 + 4.0d0*s0)*popij*0.5d0
                    else
                        rate(ii, jj) = ddep(ii, jj)*(s2 - s0)*popij
                    end if

                    !       Now the off-diagonal terms which contribute to the diagonal rate
                    !       elements: (W0 + 2W1 + W2). These are summed below.
                    !     Note that trp is not symmetric, therefore trp2 is needed;
                    !     trp is the lower diag and trp2 is the upper diag.

                    if (iroesy > 0) then
                        trp(ii, jj) = &
                            ddep(ii, jj)*(3.d0*s1 + s2 + 5.d0*s0)*popjsq*0.5d0
                        trp(jj, ii) = &
                            ddep(ii, jj)*(3.d0*s1 + s2 + 5.d0*s0)*popisq*0.5d0
                    else
                        trp(ii, jj) = ddep(ii, jj)*(s1 + s2 + s0)*popjsq
                        trp(jj, ii) = ddep(ii, jj)*(s1 + s2 + s0)*popisq
                    end if
#ifdef DEBUG_NMR
#endif

                end if  ! (ii == jj)
                k = k + 1
            end do
        end do

        !     The diagonal self-relaxation and all the cross-peak terms
        !     of the transition probabilities have already been assigned
        !     to the rate matrix.  Here the summation over i.ne.j is
        !     performed to generate the diagonal terms of the rate matrix.

        do ii = 1, natmet
            do jj = 1, natmet
                rate(ii, ii) = rate(ii, ii) + trp(ii, jj)
            end do
#ifdef DEBUG_NMR
            write (6, *) 'rate:', ii, ii, rate(ii, ii)
#endif
        end do

        ! --- filling up the whole rate matrix:

        do ii = 1, natmet
            do jj = ii + 1, natmet
                rate(ii, jj) = rate(jj, ii)
#ifdef DEBUG_NMR
                write (6, *) 'rate:', ii, jj, rate(ii, jj)
#endif
            end do
        end do

        return
    end subroutine calrate
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine corf here]
#ifdef NMODE
    subroutine corf(x, khyd, lhyd, kreal, lreal, newf, amass)
        implicit none
        integer :: khyd, lhyd, kreal, lreal
        logical :: newf
        _REAL_ :: x, amass

        !     ----- to calculate dipole-dipole correlation functions
        !           from normal modes, using Henry-Szabo Eq. (45)

        !     ----  note that the "x" array is really xhyd, so that "khyd"
        !           and "lhyd" are atom numbers in the nath scheme.
        !           The "vect" array is indexed by absolute atom numbers, and
        !           "kreal" and "lreal" are pointers to that.

        !           The first "iscale-1" frequencies will be the fitted values,
        !           held at the end of the x array; the remainder will be the
        !           "true" frequencies, held in "freq".
        !           These "true" frequencies will in fact be scaled by xkappa,
        !           which is held in x(3*natom + iscale).

        !           Input parameter "newf" is true if this is the first time
        !           that corf has been called since the effective frequencies
        !           have been changed.

        !     ----  returns:
        !            gamma_nmr: motional correction factor to cross-relaxation rate
        !            dgamma(1->3): dervative of gamma_nmr with respect to
        !                         (xk-xl), (yk-yl), (zk-zl)
        !            dgamma(4->iscale+2): derivative of gamma_nmr with respect to
        !                         frq(n), expressed in wavenumbers.
        !            dgamma(iscale+3): derivative of gamma_nmr with respect to
        !                         scaling factor for fixed frequencies

        !    ---- notes on other variables:
        !            dfac(n) is the frequency-dependent part of the formula to
        !            get cartesian fluctations for the n-th mode, along the
        !            khyd-lhyd direction.  The real cartesian fluctuation is
        !            given by the product of dfac and vfac, as shown below.

        !            ddfac(n) is the log derivative of dfac(n) with respect to the
        !            n-th scaling factor (which is the frequency for 1 -> iscale-1,
        !            and is the global scaling factor for n=iscale. Hence,
        !            to get the full derivative of the cartesian fluctuation
        !            with respect to the n-th variable parameter, multiply
        !            dfac * vfac * ddfac.  This product is stored in ddlijc.

        integer:: i, ia, ib, isnap, ivuse, j, k, ki, kj, l, li, lj, mu, n
        _REAL_ :: amp, argq, consq, ddfac, ddlijc, del, delijc, dfac, e, &
            fre, frq, hsfac1, hsfac2, hsfac3, one, qcorr, reff, req, suma, &
            two, twopi, weff
        _REAL_ kt
        logical first
#  include "nmr.h"
#  include "md.h"
#  include "memory.h"
#  include "def_time.h"
        _REAL_ amass(*)
        _REAL_ xcopy(3*natom)
        dimension x(*), e(3), del(3, 3)
        dimension ddlijc(3, 3, mxvect), dfac(mxvect), ddfac(mxvect)
        save dfac, ddfac, first
        data first/.true./

        !     functions
        k(i) = 3*(khyd - 1) + i
        l(i) = 3*(lhyd - 1) + i
        !     end functions

        one = 1.d0
        two = 2.d0
        if (newf) then
            kt = vtemp*0.002
            !                  ----CONSQ = hc/2kT in cm
            !                       (use for quantum, Bose statistics)
            consq = 0.71942/vtemp
            do n = 1, iscale - 1
                frq = x(3*natom + n)
                fre = frq*frq/11791.79
                if (bose) then
                    argq = frq*consq
                    qcorr = argq/tanh(argq)
                    dfac(n) = kt*(qcorr/fre)
                    ddfac(n) = -(one/frq + consq*tanh(argq)/sinh(argq)**2)
                else
                    dfac(n) = kt/fre
                    ddfac(n) = -two/frq
                end if
            end do

            !   --- set up fluctuations for the "fixed" frequencies:

            ddfac(iscale) = 0.0
            do n = iscale, nvect
                if (freq(n) < omegax) then
                    frq = x(3*natom + iscale)*abs(freq(n))
                else
                    frq = freq(n)
                end if
                fre = frq*frq/11791.79
                if (bose) then
                    argq = frq*consq
                    qcorr = argq/tanh(argq)
                    ddfac(n) = -(one/x(3*natom + iscale) &
                        + consq*freq(n)*tanh(argq)/sinh(argq)**2)
                else
                    qcorr = 1.0
                    ddfac(n) = -two/x(3*natom + iscale)
                end if
                dfac(n) = kt*(qcorr/fre)
                if (frq > 6000.d0) dfac(n) = 0.d0
                write (6, *) n, frq, dfac(n)
                if (freq(n) >= omegax) ddfac(n) = 0.0
            end do
            first = .false.
            if (nmsnap > 0) then
#ifdef DEBUG_NMR

                !       --- check orthogonality:

                dot77 = 0.0
                dot78 = 0.0
                dot89 = 0.0
                n = 0
                do i = 1, 3*natom
                    if (mod(i, 3) == 1) n = n + 1
                    dot77 = dot77 + vect(i, 7)*vect(i, 7)*amass(n)
                    dot78 = dot78 + vect(i, 7)*vect(i, 8)*amass(n)
                    dot89 = dot89 + vect(i, 8)*vect(i, 9)*amass(n)
                end do
                write (6, *) 'orthog check: ', dot77, dot78, dot89
                write (6, *) (amass(i), i=1, natom)
#endif

                call amrset(54185253)

                !       ---output "snapshots" of structures averaged over the modes:

                isnap = 0
                write (6, *) 'mode snapshot ', isnap
                write (6, '(i5)') natom
                write (6, '(6f12.7)') (x(i), i=1, 3*natom)
                twopi = 6.28318d0
                do isnap = 1, nmsnap
                    do i = 1, 3*natom
                        xcopy(i) = x(i)
                    end do
                    do n = 1, nvect

                        !    here is (possibly) quantum amplitude, but classical motion:

                        !        here is Gaussian with the quantum std. deviation:

                        call gauss(0.d0, sqrt(dfac(n)), amp)

                        do i = 1, 3*natom
                            xcopy(i) = xcopy(i) + amp*vect(i, n)
                        end do
                    end do
                    write (6, *) 'mode snapshot ', isnap
                    write (6, '(i5)') natom
                    write (6, '(6f12.7)') (xcopy(i), i=1, 3*natom)
                end do
                call mexit(6, 0)
            end if
        end if

        !    ---e(i) is the unit vector along the l-k bond:

        req = 0.0
        do i = 1, 3
            e(i) = x(k(i)) - x(l(i))
            req = req + e(i)*e(i)
        end do
        req = sqrt(req)
        do i = 1, 3
            e(i) = e(i)/req
        end do

        !     ----- calculate the correlation matrix for delta
        !        as in Eq. 7.16 of lamm and szabo j. chem. phys. 85, 7334 (1986)
        !        Note that the rhs of this eq. should be multiplied by kT

        ivuse = 1
43      continue
        do i = 1, 3
            ki = 3*(kreal - 1) + i
            li = 3*(lreal - 1) + i
            do j = 1, 3
                kj = 3*(kreal - 1) + j
                lj = 3*(lreal - 1) + j
                del(i, j) = 0.0
                ddlijc(i, j, iscale) = 0.0

                !  --- do fitted frequencies, for which we will need derivatives:

                do n = 1, iscale - 1
                    vfac = (vect(ki, n) - vect(li, n))*(vect(kj, n) - vect(lj, n))
                    delijc = dfac(n)*vfac
                    del(i, j) = del(i, j) + delijc
                    ddlijc(i, j, n) = delijc*ddfac(n)
                end do

                ! --- now add in remaining fixed frequency modes:

                if (per_mode) then
                    n = ivuse
                    vfac = (vect(ki, n) - vect(li, n))*(vect(kj, n) - vect(lj, n))
                    del(i, j) = del(i, j) + dfac(n)*vfac
                    ddlijc(i, j, iscale) = ddlijc(i, j, iscale) + &
                        dfac(n)*vfac*ddfac(n)
                else
                    do n = iscale, nvect
                        vfac = (vect(ki, n) - vect(li, n))*(vect(kj, n) - vect(lj, n))
                        del(i, j) = del(i, j) + dfac(n)*vfac
                        ddlijc(i, j, iscale) = ddlijc(i, j, iscale) + &
                            dfac(n)*vfac*ddfac(n)
                    end do
                end if
            end do
        end do

        !    ---- use eq. 45 of Henry & Szabo [JCP 82:4753(1985)]
        !         to get motional correction for dipole-dipole correlation:

        !         (The code with ihsful=1 uses the full Henry-Szabo
        !         formula; that without this defined leaves out the terms
        !         dependent upon A1 and A2, i.e. the bond-length dependent
        !         terms.)

        if (ihsful == 1) then
            hsfac1 = 15.d0
            hsfac2 = 5.d0
            hsfac3 = 10.d0
        else
            hsfac1 = 3.d0
            hsfac2 = 1.d0
            hsfac3 = 2.d0
        end if
        suma = 0.0
        reff = 0.0
        do i = 1, 3
            suma = suma - 3.0*del(i, i)
            reff = reff + del(i, i)
            do j = 1, 3
                suma = suma + hsfac1*del(i, j)*e(i)*e(j)
                reff = reff - del(i, j)*e(i)*e(j)
            end do
        end do
        weff = 1.0 + (0.5/req**2)*suma
        gamma_nmr = weff**2
        reff = 0.5*reff/req
        if (per_mode) then
            write (6, *) 'ivuse, gamma:', ivuse, gamma_nmr
            ivuse = ivuse + 1
            if (ivuse <= nvect) goto 43
            stop
        else
            write (6, *) 'gamma, reff:', gamma_nmr, reff
        end if

        !   ---  now get dgamma = derivative of gamma_nmr with respect to
        !          (xk-xl), (yk-yl), (zk-zl)

        do mu = 1, 3
            dgamma(mu) = 0.0
            do ib = 1, 3
                dgamma(mu) = dgamma(mu) + hsfac2*e(ib)*del(mu, ib) &
                    + e(mu)*del(ib, ib)
                do ia = 1, 3
                    dgamma(mu) = dgamma(mu) - hsfac3*e(ia)*e(ib)*e(mu)*del(ia, ib)
                end do
            end do
            dgamma(mu) = weff*6.*dgamma(mu)/req**3
        end do

        !  --- get more of dgamma: derivatives with respect to the weights

        do n = 1, iscale
            suma = 0.0
            do i = 1, 3
                suma = suma - 3.0*ddlijc(i, i, n)
                do j = 1, 3
                    suma = suma + hsfac1*ddlijc(i, j, n)*e(i)*e(j)
                end do
            end do
            dgamma(n + 3) = weff*suma/req**2
        end do

        return
    end subroutine corf
#endif

#ifdef ORIGDERIV

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dinten here]
    subroutine dinten(vecs, facj, ii, jj, c, dint)
#else

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dinten here]
    subroutine dinten(amat, ii, jj, ddrat, dorat, dint, taum)
#endif

        ! Subroutine Derivatives of INTENsities:

        !     This is a subroutine to calculate the derivatives of the
        !     intensites w.r.t the protons coordinates.  Dint contains
        !     these derivatives (indexed in the natmet scheme); all the
        !     remaining arguments are input variables.  The cross peak
        !     labels are ii and jj.

        implicit none
#  include "nmr.h"

        _REAL_ dint(3*ma + mxvar)
        integer ii, jj
#ifdef ORIGDERIV
        _REAL_ c(3*ma + mxvar, lst), facj(ma, ma), vecs(ma, ma)
        integer lm
        integer nath3
        integer r, u
        _REAL_ temp

        ! ---- construct product of K matrix, facj and one L, then
        !        multiply by other row of L(trans):

        nath3 = 3*nath + iscale

        ! --- (zeroing out of Dint is now handled in remarc!)

        lm = 0
        do r = 1, natmet

            ! ----  off diagonal terms:

            do u = 1, r - 1
                lm = lm + 1
                temp = facj(r, u) &
                    *(vecs(ii, r)*vecs(jj, u) + vecs(ii, u)*vecs(jj, r))
                call D_OR_S() axpy(nath3, temp, c(1, lm), 1, dint, 1)
            end do

            ! ---- diagonal terms:

            lm = lm + 1
            temp = vecs(ii, r)*vecs(jj, r)*facj(r, r)
            call D_OR_S() axpy(nath3, temp, c(1, lm), 1, dint, 1)
        end do

#else

        !  --- here use Ping Yips "fast" derivative routine,  J. Biomol. NMR
        !        3:361-365(1993):

        _REAL_ amat(ma, ma, 5), w(5)
        _REAL_ dorat(3*ma, ma), ddrat(3*ma, ma)
        _REAL_ cutoff
        parameter(cutoff=0.001d0)
        integer kk, ll, m
        integer muk, mm, mum
        integer nath3
        _REAL_ ail1, ail2, ail3, ail4, ail5
        _REAL_ ajl1, ajl2, ajl3, ajl4, ajl5
        _REAL_ f1
        _REAL_ taum

        !   --- five-point Gaussian quadrature:

        save w
        data w/0.1184634425d0, 0.2393143352d0, 0.2844444444d0, &
            0.2393143352d0, 0.1184634425d0/

        nath3 = 3*nath

        !   --- part arising from ddrat:

        do kk = 1, natmet
            if (amat(ii, kk, 5) < cutoff .and. amat(jj, kk, 5) < cutoff) &
                cycle
            f1 = -taum*(w(1)*(amat(ii, kk, 1)*amat(jj, kk, 5) &
                + amat(ii, kk, 5)*amat(jj, kk, 1)) &
                + w(2)*(amat(ii, kk, 2)*amat(jj, kk, 4) &
                + amat(ii, kk, 4)*amat(jj, kk, 2)) &
                + w(3)*(amat(ii, kk, 3)*amat(jj, kk, 3)))
            do muk = 1, nath3
                dint(muk) = dint(muk) + f1*ddrat(muk, kk)
            end do
        end do

        !   --- part arising from dorat:

        do ll = 1, natmet
            ajl5 = amat(jj, ll, 5)
            ail5 = amat(ii, ll, 5)
            if (ajl5 < cutoff .and. ail5 < cutoff) cycle
            ajl1 = amat(jj, ll, 1)
            ajl2 = amat(jj, ll, 2)
            ajl3 = amat(jj, ll, 3)
            ajl4 = amat(jj, ll, 4)
            ail1 = amat(ii, ll, 1)
            ail2 = amat(ii, ll, 2)
            ail3 = amat(ii, ll, 3)
            ail4 = amat(ii, ll, 4)
            mum = 0
            do m = 1, nath
                mm = inn(m)
                f1 = -taum*(w(1)*(amat(ii, mm, 5)*ajl1 &
                    + ail5*amat(jj, mm, 1) &
                    + amat(ii, mm, 1)*ajl5 &
                    + ail1*amat(jj, mm, 5)) &
                    + w(2)*(amat(ii, mm, 4)*ajl2 &
                    + ail4*amat(jj, mm, 2) &
                    + amat(ii, mm, 2)*ajl4 &
                    + ail2*amat(jj, mm, 4)) &
                    + w(3)*(amat(ii, mm, 3)*ajl3 &
                    + ail3*amat(jj, mm, 3)))

                dint(mum + 1) = dint(mum + 1) + f1*dorat(mum + 1, ll)
                dint(mum + 2) = dint(mum + 2) + f1*dorat(mum + 2, ll)
                dint(mum + 3) = dint(mum + 3) + f1*dorat(mum + 3, ll)
                mum = mum + 3
            end do
        end do
#endif
        return
    end subroutine dinten

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine drates here]
    subroutine drates(ddep, dddep, rate, trp, dorat, ddrat)

        !  Subroutine Derivatives of RATES:

        !     This subroutine calculates the derivatives
        !     of the transition rates with repect to the
        !     cordinates of the individual protons. On
        !     input, it needs the x, dsq2, rate, and trp
        !     arrays.  The output derivatives are placed
        !     in the ddrat and dorat arrays.

        implicit none
        integer:: i, ii, jj, n
        _REAL_ :: dddep, ddep, ddrat, dorat, rate, tii, trp
#  include "nmr.h"

        dimension ddep(ma, ma), dddep(3, ma, ma)
        dimension rate(ma, ma), trp(ma, ma)
        dimension dorat(3, ma, ma), ddrat(3, ma, ma)
        dimension tii(3)

        ! --- Now compute the derivatives: the diagonal elements
        !       of grad R go into ddrat, the off-diagonal elements
        !       (just a single column) into dorat.
        !       Here "i" and "j" are indexing in the nath scheme,
        !       "ii" and "jj" are the corresponding natmet indices.
        !       Note that the second index of ddrat and dorat are
        !       indexed by the nath scheme, while the third index
        !       goes by the natmet scheme; cf. the way these arrays
        !       are used in the kmat subroutine.

#ifdef DEBUG_NMR
        write (6, *) 'drates:'
#endif
        do i = 1, nath
            ii = inn(i)
            do n = 1, 3
                tii(n) = 0.0d0
            end do
            do jj = 1, natmet
                if (ii == jj) cycle
                do n = 1, 3
                    dorat(n, i, jj) = dddep(n, i, jj)*rate(ii, jj)/ddep(ii, jj)
                    ddrat(n, i, jj) = dddep(n, i, jj)*trp(jj, ii)/ddep(ii, jj)
                    tii(n) = tii(n) + dddep(n, i, jj)*trp(ii, jj)/ddep(ii, jj)
                end do
#ifdef DEBUG_NMR
                write (6, *) 'do:', i, jj, (dorat(n, i, jj), n=1, 3)
                write (6, *) 'dd:', i, jj, (ddrat(n, i, jj), n=1, 3)
#endif
            end do

            do n = 1, 3
                dorat(n, i, ii) = 0.0d0
                ddrat(n, i, ii) = tii(n)
            end do
#ifdef DEBUG_NMR
            write (6, *) i, ii, (dorat(n, i, ii), n=1, 3)
            write (6, *) i, ii, (ddrat(n, i, ii), n=1, 3)
#endif
        end do
#ifdef NMODE

        !   --- get derivatives of rate matrix elements with respect to
        !         normal mode frequencies:

        do ii = 1, natmet
            do jj = 1, natmet
                if (ii == jj) goto 40
                do n = 1, iscale
                    dratg(ii, ii, n) = dratg(ii, ii, n) + dratg(ii, jj, n)*trp(ii, jj)
                    dratg(ii, jj, n) = dratg(ii, jj, n)*rate(ii, jj)
                end do
            end do
        end do

#endif
        return
    end subroutine drates

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine indexn here]
    subroutine indexn(ix, ih, iin)

        !  Subroutine INDEX for Noe calculations:

        !   --- reads the tau and ihyd arrays that set up the sub-molecules,
        !       and sets up all the parameters in the /methyl/ common
        !       block that describes the molecule

        implicit none
        integer:: i, id1, idummy, ie1, ie2, iin, ix, iz, k, mtyp, n, &
            natomx, ngrp, nn
        _REAL_ :: dummy, one, sqrt2, sqrt3, three, two
#  include "nmr.h"
#  include "memory.h"
#  include "box.h"

        dimension ix(*)
        character(len=4) ih(*)

        character(len=3) rname
        character(len=4) aname
        parameter(three=3.0d0, one=1.0d0, two=2.0d0)
        if (natom > matom) then
            write (6, *) 'MATOM is not big enough!'
            call mexit(6, 1)
        end if
        sqrt2 = sqrt(two)
        sqrt3 = sqrt(three)
        natmet = 0
        nath = 0

        !   ---- call the usual AMBER group identifiers to determine the
        !        submolecule, then strip away all non-hydrogens

        natomx = natom
        do i = 1, natomx
            ihyd(i) = 0
            inatom(i) = 0
        end do
        call rgroup(natomx, natc, nres, ngrp, ix(i02), ih(m02), ih(m04), ih(m06), &
            ih(m08), ihyd, idummy, idummy, idummy, idummy, dummy, dummy, &
            .false., .false., .false., .false., iin, .true.)
        do i = 1, natomx
            if (ihyd(i) /= 1 .or. resat(i) (1:1) /= 'H') ihyd(i) = 0
            if ((id2o == 1) .and. &
                ((resat(i) (1:2) == 'H ') .or. &
                (resat(i) (1:2) == 'HO' .and. resat(i) (6:8) /= 'TYR') .or. &
                (resat(i) (1:3) == 'H1 ' .and. resat(i) (6:8) == 'GUA') .or. &
                (resat(i) (1:3) == 'H1 ' .and. resat(i) (6:7) == 'DG') .or. &
                (resat(i) (1:3) == 'H21' .and. resat(i) (6:8) == 'GUA') .or. &
                (resat(i) (1:3) == 'H21' .and. resat(i) (6:7) == 'DG') .or. &
                (resat(i) (1:3) == 'H22' .and. resat(i) (6:8) == 'GUA') .or. &
                (resat(i) (1:3) == 'H22' .and. resat(i) (6:7) == 'DG') .or. &
                (resat(i) (1:3) == 'H41' .and. resat(i) (6:8) == 'CYT') .or. &
                (resat(i) (1:3) == 'H41' .and. resat(i) (6:7) == 'DC') .or. &
                (resat(i) (1:3) == 'H42' .and. resat(i) (6:8) == 'CYT') .or. &
                (resat(i) (1:3) == 'H42' .and. resat(i) (6:7) == 'DC') .or. &
                (resat(i) (1:3) == 'H3 ' .and. resat(i) (6:8) == 'THY') .or. &
                (resat(i) (1:3) == 'H3 ' .and. resat(i) (6:7) == 'DT') .or. &
                (resat(i) (1:3) == 'H61' .and. resat(i) (6:8) == 'ADE') .or. &
                (resat(i) (1:3) == 'H61' .and. resat(i) (6:7) == 'DA') .or. &
                (resat(i) (1:3) == 'H62' .and. resat(i) (6:7) == 'DA') .or. &
                (resat(i) (1:3) == 'H62' .and. resat(i) (6:8) == 'ADE'))) &
                ihyd(i) = 0
#ifdef DEBUG_NMR
            if (ihyd(i) /= 0) write (6, '(a14)') resat(i)
#endif
        end do
        k = -3
        do i = 1, natomx
            k = k + 3
            if (ihyd(i) == 0) cycle
            nath = nath + 1
            ihyp(nath) = i

            ! ---    use resat array to determine type of this proton; this
            !        code is specific to the Amber all-atom IUPAC names

            mtyp = 2
            rname = resat(i) (6:8)
            aname = resat(i) (1:4)
            if (aname == 'HD11' .or. aname == 'HD12' .or. &
                aname == 'HG21' .or. aname == 'HG22' .or. &
                aname == 'H71 ' .or. aname == 'H72 ' .or. &
                aname == 'HNZ1' .or. aname == 'HNZ2' .or. &
                aname == 'HD21' .or. aname == 'HD22') then
                mtyp = 1
            else if (aname == 'HD13' .or. aname == 'HD23' .or. &
                aname == 'HNZ3' .or. &
                aname == 'HG23' .or. aname == 'H73 ') then
                mtyp = 3
            else if (rname == 'ACE') then
                if (aname == 'H1  ' .or. aname == 'H2  ') mtyp = 1
                if (aname == 'H3  ') mtyp = 3
            else if (rname == 'ALA') then
                if (aname == 'HB1 ' .or. aname == 'HB2 ') mtyp = 1
                if (aname == 'HB3 ') mtyp = 3
            else if (rname == 'ILE') then
                if (aname == 'HD1 ' .or. aname == 'HD2 ') mtyp = 1
                if (aname == 'HD3 ') mtyp = 3
            else if (rname == 'THR') then
                if (aname == 'HG1 ' .or. aname == 'HG2 ') mtyp = 1
                if (aname == 'HG3 ') mtyp = 3
            else if (rname == 'VAL') then
                if (aname == 'HG11' .or. aname == 'HG12') mtyp = 1
                if (aname == 'HG13') mtyp = 3
            else if (rname == 'MET') then
                if (aname == 'HE1 ' .or. aname == 'HE2 ') mtyp = 1
                if (aname == 'HE3 ') mtyp = 3
            else if (rname == 'PHE' .or. rname == 'TYR') then
                if (aname == 'HD1 ' .or. aname == 'HE1 ') mtyp = 4
                if (aname == 'HE2 ') mtyp = 5
                if (aname == 'HD2 ') mtyp = 6
            end if

            if (mtyp == 1) then
                m2(nath) = 1
                pop(nath) = 0.0d0
            else if (mtyp == 3) then
                m2(nath) = 3
                pop(nath) = sqrt3
                natmet = natmet + 1
                popn(natmet) = sqrt3
            else if (mtyp == 4) then
                m2(nath) = 4
                pop(nath) = 0.0d0
            else if (mtyp == 5 .or. mtyp == 6) then
                m2(nath) = 5
                pop(nath) = sqrt2
                natmet = natmet + 1
                popn(natmet) = sqrt2

                ! --- rearrange ring atoms to get HE1&HE2, HD1&HD2 next to each other:

                if (mtyp == 6) then
                    id1 = ihyp(nath - 4)
                    ie1 = ihyp(nath - 3)
                    iz = ihyp(nath - 2)
                    ie2 = ihyp(nath - 1)
                    ihyp(nath - 4) = iz
                    ihyp(nath - 3) = ie1
                    ihyp(nath - 2) = ie2
                    ihyp(nath - 1) = id1
                    m2(nath - 4) = 2
                    m2(nath - 3) = 4
                    m2(nath - 2) = 5
                    m2(nath - 1) = 4
                    pop(nath - 4) = 1.0d0
                    pop(nath - 3) = 0.0d0
                    pop(nath - 2) = sqrt2
                    pop(nath - 1) = 0.0d0
                end if

            else
                natmet = natmet + 1
                popn(natmet) = 1.0d0
                pop(nath) = one
                m2(nath) = 2
            end if  ! (mtyp == 1)
            inatom(i) = natmet
#ifdef DEBUG_NMR
            write (6, *) i, nath, natmet, resat(i), ihyp(nath), m2(nath), &
                pop(nath), popn(natmet)
#endif
        end do
#ifdef DEBUG_NMR
        write (6, *) 'nath, natmet are: ', nath, natmet
#endif
        if (nath > ma) then
            write (6, *) 'Maximum allowed value for nath is ', ma
            call mexit(6, 1)
        end if

        !  --- set up inn array that points from nath -> natmet scheme:

        nn = 1
        do n = 1, nath
            inn(n) = nn
            if (m2(n) == 2 .or. m2(n) == 3 .or. m2(n) == 5) nn = nn + 1
        end do
#ifdef DEBUG_NMR
        write (6, 30) (m2(i), i=1, nath)
30      format(' m2:'/(10i5))
        write (6, 31) (pop(i), i=1, nath)
31      format(' pop:'/(5f10.5))
        write (6, 32) (ihyp(i), i=1, nath)
32      format(' ihyp:'/(10i5))
        write (6, 33) (popn(i), i=1, natmet)
33      format(' popn:'/(5f10.5))
        write (6, 34) (inatom(i), i=1, natomx)
34      format(' inatom:'/(10i5))
#endif
        return
    end subroutine indexn
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine kmat here]
#ifdef ORIGDERIV
    subroutine kmat(vecs, dorat, ddrat, c)
        implicit none
        _REAL_ :: vecs, dorat, ddrat, c

        !  Subroutine K MATrix:

        !   Subroutine kmat calculates the K matrix from the derivatives of
        !     the rates given by drates. The output is the K matrix
        !     which will be used in the derivatives of the NOEs.
        !     N.B. K is called "C" here and elsewhere to conform
        !     to standard fortran typing conventions.

        integer:: ip, k, kf, l, lm, m, muk, n, nath3
        _REAL_ :: dol, rate, svecs1
#  include "nmr.h"

        dimension vecs(ma, ma), c(3*ma + mxvar, lst), rate(ma, ma)
        dimension dorat(3*ma, ma), ddrat(3*ma, ma)
        dimension dol(3*ma, ma)
        nath3 = 3*nath

        !-----the matrix K is constructed by the following do loops

        !   --- first, do some preliminary matrix multiplies:

        do m = 1, natmet
            do k = 1, 3*nath
                dol(k, m) = 0.0
            end do
            do ip = 1, natmet
                do k = 1, 3*nath
                    dol(k, m) = dol(k, m) + dorat(k, ip)*vecs(ip, m)
                end do
            end do
        end do

        !   --- now for the part that depends upon the diagonal rate elements:

        lm = 0
        kf = 3*nath + 1
        do l = 1, natmet
            do m = 1, l
                lm = lm + 1
                do muk = 1, 3*nath + 1
                    c(muk, lm) = 0.0
                end do
                do ip = 1, natmet
                    svecs1 = vecs(ip, l)*vecs(ip, m)
                    call D_OR_S() axpy(nath3, svecs1, ddrat(1, ip), 1, c(1, lm), 1)
                end do

                !  --- now do the part that depends on the off-diagonal rate elements:

                k = 0
                do n = 1, nath
                    c(k + 1, lm) = c(k + 1, lm) + vecs(inn(n), l)*dol(k + 1, m) &
                        + vecs(inn(n), m)*dol(k + 1, l)
                    c(k + 2, lm) = c(k + 2, lm) + vecs(inn(n), l)*dol(k + 2, m) &
                        + vecs(inn(n), m)*dol(k + 2, l)
                    c(k + 3, lm) = c(k + 3, lm) + vecs(inn(n), l)*dol(k + 3, m) &
                        + vecs(inn(n), m)*dol(k + 3, l)
                    k = k + 3
                end do
            end do
        end do
#  ifdef NMODE

        !   ---  add to the K matrix the derivatives with respect to
        !         the normal mode frequencies:

        do n = 1, iscale
            do jj = 1, natmet
                do ii = 1, natmet
                    dol(ii, jj) = 0.0
                end do
                do kk = 1, natmet
                    do ii = 1, natmet
                        dol(ii, jj) = dol(ii, jj) + dratg(ii, kk, n)*vecs(kk, jj)
                    end do
                end do
            end do

            lm = 0
            do ii = 1, natmet
                do jj = 1, ii
                    lm = lm + 1
                    sum = 0.0
                    do kk = 1, natmet
                        sum = sum + vecs(kk, ii)*dol(kk, jj)
                    end do
                    c(3*nath + n, lm) = sum
                end do
            end do
        end do
#  endif
        return
    end subroutine kmat
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine remarc here]
    subroutine remarc(ddep, dddep, f, xhyd, ksub, dorat, ddrat, &
        rate, trp, dint)

        !   Subroutine RElaxation MAtrix Refinement Code:

        !       Calculates theoretical intensities from diagonalized rate matrix.
        !       also calculates the derivatives of the intensities w.r.t the
        !       proton cordinates.  This routine is called once for each submole-
        !       cule, and mainly does the following:

        !       Call calrate(ddep,rate,trp)
        !       Call LAPACK for diagonalization
        !         [or: use perturbation estimates; see compilation flags, below]

        !       For the derivatives: calls drates, dinten, and kmat or amatg

        !       For mixing times less than tausw, the perturbation expansions are
        !           used; for mixing times greater than tausw, the exact expressions
        !           are used.

        !  Inputs:

        !     ddep     -- distance dependent portion of the rate matrix elements;
        !                     computed in caldis.
        !     dddep    -- derivative of ddep with respect to coordinates; also
        !                     computed in caldis.
        !     f        -- force vector, to be updated here
        !     xhyd     -- coordinates of the protons in the current submolecule
        !     enoe     -- energy penalty for NOE volumes; to be updated here
        !     ksub     -- number of the current submolecule; just for printing
        !     dorat
        !     ddrat    -- derivatives of the off-diagonal (dorat) and diagonal
        !                     (ddrat) elements of the rate matrix; computed
        !                     in drates (called below); included in the argument
        !                     list here just for the sake of memory management.
        !     rate     -- relaxation rate matrix; computed in calrate (called
        !                     below; included in argument list for memory management
        !     trp      -- work array
        !     Dint     -- derivatives of the computed intensities with respect
        !                     to variables (coordinates and other fitting variables)

        !  Outputs:

        !     enoe and f are updated.

        !  --- conditional compilation flags:

        !      THIRDO  to add 3rd order terms to perturbation expansion

        !      ORIGDERIV  for exact derivatives, if this is set, use the
        !                 original formulation given in Yip & Case, J.
        !                 Magn. Reson. 83:643 (1989);  if not set, use the
        !                 integral formulation given by Yip, J. Biomol. NMR
        !                 3:361 (1993).

        use constants, only : zero, one, two, half, third, fourth, sixth
        !implicit none
#  include "nmr.h"
#  include "memory.h"
#  include "md.h"
#  include "def_time.h"

        _REAL_ xhyd(*), f(*)
        integer ksub
        _REAL_ rate(ma, ma)
        _REAL_ dorat(3*ma, ma), ddrat(3*ma, ma), dint(3*ma + mxvar)
        _REAL_ trp(ma, ma), ddep(ma, ma), dddep(3, ma, ma)
        _REAL_ ddisc(3*ma), dwijdx(3*ma), etp(ma)
        character(len=1) uplo, jobz
#ifdef THIRDO
        _REAL_ delta(ma, ma)
#endif
#ifdef ORIGDERIV
        _REAL_ facj(ma, ma), c(3*ma + mxvar, lst)
#else
        _REAL_ amat(ma, ma, 5)
        integer iamat(ma)
#endif
        integer muki(ma), mukj(ma), mukl(ma)
#ifdef THIRDO
        integer mukm(ma)
#endif
        _REAL_ pr(lst), b(ma, 9), eig(ma), et(ma), vecs(ma, ma)

        _REAL_ acalc
        integer muk
        integer k
        integer j
        integer i
        integer ier
        integer imix
        _REAL_ taum
        integer n
        integer ipeak
        _REAL_ acalc1
        integer ii
        integer jj
        _REAL_ firsto
        _REAL_ delij
        _REAL_ disc
        _REAL_ disci
        integer kk
        _REAL_ secndo
        integer imuk
        integer jmuk
        integer imukmx
        integer jmukmx
        integer ll
        _REAL_ fac
        _REAL_ delil
        _REAL_ deljl
        _REAL_ secnda
        _REAL_ secndb
        _REAL_ secndc
        _REAL_ scont
        _REAL_ fac1
        _REAL_ fac2
        integer lmuk
        integer lmukmx
        _REAL_ facil
        _REAL_ facjl
        _REAL_ facij
        _REAL_ thirdo
        _REAL_ penalty
        _REAL_ aexpn
        _REAL_ arngen
        _REAL_ ainv
        _REAL_ weight
        _REAL_ dpen
        integer i3real
        integer i3
        _REAL_ fsx
        parameter(fsx=-5.d0/6.d0)
        data uplo, jobz/'U', 'V'/

        !=======================================================================

        !    Section 1:  Compute rate matrix and its derivatives; diagonalize
        !                this if "exact" volumes are needed.

        !=======================================================================

        !----initial zero-out, to allow for summing peaks together:

        acalc = zero
        do muk = 1, 3*nath + iscale
            dint(muk) = zero
        end do

        !----calculate the rate matrix

        call timer_start(TIME_CALRATE)
        call calrate(ddep, rate, trp)
        call timer_stop_start(TIME_CALRATE, TIME_DRATES)

        !-----calculate the derivatives of the rates w.r.t the H atoms coord.

        call drates(ddep, dddep, rate, trp, dorat, ddrat)
        call timer_stop_start(TIME_DRATES, TIME_DSPEV)
#ifdef THIRDO
        do i = 1, natmet
            do j = i + 1, natmet
                delta(i, j) = one/(-rate(i, i) + rate(j, j))
                delta(j, i) = -delta(i, j)
            end do
        end do
#endif
        if (tausw < 9.) then
            !                            (we will be calculating exact volumes)

            !-----pack the upper triangular rate matrix by column
            !      into an array pr

            k = 1
            do j = 1, natmet
                do i = 1, j
                    pr(k) = rate(i, j)
                    k = k + 1
                end do
            end do

            !-----calc eigenvalues and vectors

            call D_OR_S() spev(jobz, uplo, natmet, pr, eig, vecs, ma, b, ier)
            call timer_stop(TIME_DSPEV)
#ifdef ORIGDERIV
            call timer_start(TIME_KMAT)

            !------calc the K_prime matrix from the derivatives of the rates and the
            !      eigenvectors; this is essentially the scheme in the original
            !      1989 Yip & Case paper.

            call kmat(vecs, dorat, ddrat, c)
            call timer_stop(TIME_KMAT)
#endif
        end if

        !=======================================================================

        !   Section 2:  Grand loop over mixing times, then over peaks within
        !               each mixing time.

        !=======================================================================

        ! --- loop over mixing times

        do imix = 1, nummt
            taum = emix(imix)
            if (taum < tausw) then
                !                              (set up preliminary stuff for perturbation)

                call timer_start(TIME_REMARC)
                do n = 1, natmet
                    if (rate(n, n)*taum < 50.) then
                        etp(n) = exp(-rate(n, n)*taum)
                    else
                        etp(n) = zero
                    end if
                end do
                call timer_stop(TIME_REMARC)
            else
                call timer_start(TIME_DINTEN)
                !                              (set up preliminary stuff for exact calc.)

                do n = 1, natmet
                    if (eig(n)*taum < 50.) then
                        et(n) = exp(-eig(n)*taum)
                    else
                        et(n) = zero
                    end if
                end do
#ifdef ORIGDERIV

                !   --- compute function of eigenvalues to be used later:

                do i = 1, natmet
                    do j = i + 1, natmet
                        facj(i, j) = (et(i) - et(j))/(eig(i) - eig(j))
                        facj(j, i) = facj(i, j)
                    end do
                    facj(i, i) = -taum*et(i)
                end do
#else

                !   --- compute intensity matrix at Gaussian quadrature points:

                !   --- (get unique list of atoms referred to:)

                do i = 1, natmet
                    iamat(i) = 0
                end do
                do ipeak = 1, npeak(imix)
                    iamat(inatom(ihp(imix, ipeak))) = 1
                    iamat(inatom(jhp(imix, ipeak))) = 1
                end do

                call amatg(vecs, eig, taum, amat, iamat)
#endif
                call timer_stop(TIME_DINTEN)
            end if

            ! --- loop over all observed peaks with this mixing time

            do ipeak = 1, npeak(imix)

                !         acalc1 will be computed intensity for each peak, before summing
                !         for overlaps:

                acalc1 = zero

                ii = inatom(ihp(imix, ipeak))
                jj = inatom(jhp(imix, ipeak))
                if (ii <= 0 .or. jj <= 0) then
                    write (6, *) 'Bad submolecule: ', imix, ipeak, ii, jj, &
                        ihp(imix, ipeak), jhp(imix, ipeak)
                    call mexit(6, 1)
                end if
                if (taum < tausw) then
                    call timer_start(TIME_REMARC)

                    !=======================================================================

                    !   Section 2a:  For perturbation theory estimates, compute here
                    !                the peak intensities and derivatives, through
                    !                second or third order in the off-diagonal rate matrix
                    !                elements.

                    !=======================================================================

                    !  --- use P. Yips perturbation expansion to get an estimate
                    !       of the insensties:  see Chem. Phys. Lett. 161:50-54 (1989).
                    !       Diagonal peaks are expanded through second order, using
                    !       expressions derived by dac.

                    !  --- the first order term:

                    if (ii == jj) then

                        !  ---     diagonal:

                        firsto = etp(ii)
                        do muk = 1, 3*nath
                            dint(muk) = dint(muk) - taum*firsto*ddrat(muk, ii)
                        end do
                    else

                        !  ---      off-diagonal:

                        delij = one/(rate(jj, jj) - rate(ii, ii))
                        disc = sqrt(fourth/(delij*delij) + rate(ii, jj)**2)
                        disci = one/disc

                        !    ---- basic first order expression is:
                        !             firsto=-rate(ii,jj)
                        !    .          *exp(-Half*taum*(rate(ii,ii)+rate(jj,jj)))
                        !    .          *sinh(disc*taum)/disc

                        !         here we break it down to avoid sinh evaluation if we can,
                        !            and to reduce exponential overflow:

                        if (disc*taum < 0.2) then
                            firsto = -rate(ii, jj)*sqrt(etp(ii)*etp(jj))*taum
                        else if (disc*taum < 4.0) then
                            firsto = -rate(ii, jj)*sqrt(etp(ii)*etp(jj)) &
                                *sinh(disc*taum)*disci
                        else if (disc*taum > 50.0) then
                            write (0, *) 'remarc: ', ii, jj, rate(ii, ii), rate(jj, jj), &
                                rate(ii, jj), disc
                            write (0, *) ipeak, imix, ihp(imix, ipeak), jhp(imix, ipeak)
                            firsto = zero
                        else
                            firsto = -rate(ii, jj)*exp(taum*disc)*sqrt(etp(ii)*etp(jj)) &
                                *disci*half
                        end if

                        !      ---dwijdx(muk) will contain d(rate(ii,jj))/dx(muk):

                        do muk = 1, 3*nath
                            dwijdx(muk) = zero
                        end do
                        muk = -3
                        do k = 1, nath
                            muk = muk + 3
                            kk = inn(k)
                            if (kk == ii) then
                                dwijdx(muk + 1) = dorat(muk + 1, jj)
                                dwijdx(muk + 2) = dorat(muk + 2, jj)
                                dwijdx(muk + 3) = dorat(muk + 3, jj)
                            else if (kk == jj) then
                                dwijdx(muk + 1) = dorat(muk + 1, ii)
                                dwijdx(muk + 2) = dorat(muk + 2, ii)
                                dwijdx(muk + 3) = dorat(muk + 3, ii)
                            end if
                        end do

                        !      --- ddsic(muk) will contain d(disc)/dx(muk):

                        do muk = 1, 3*nath
                            ddisc(muk) = (fourth*(ddrat(muk, jj) - ddrat(muk, ii))/delij + &
                                rate(ii, jj)*dwijdx(muk))*disci
                        end do

                        !      --- now for the final expression:

                        do muk = 1, 3*nath
                            dint(muk) = dint(muk) + (dwijdx(muk)/rate(ii, jj) - &
                                half*taum*(ddrat(muk, ii) + ddrat(muk, jj)) + &
                                taum*ddisc(muk)/tanh(disc*taum) - &
                                ddisc(muk)*disci)*firsto
                        end do
                    end if

                    !  --- the second order term:

                    secndo = zero
                    muk = -3
                    imuk = 0
                    jmuk = 0
                    do k = 1, nath
                        muk = muk + 3
                        kk = inn(k)
                        if (kk == ii) then
                            imuk = imuk + 1
                            muki(imuk) = muk
                        else if (kk == jj) then
                            jmuk = jmuk + 1
                            mukj(jmuk) = muk
                        end if
                    end do
                    imukmx = imuk
                    jmukmx = jmuk

                    do ll = 1, natmet
                        if (ll == ii .or. ll == jj) cycle
                        fac = rate(ii, ll)*rate(ll, jj)
                        if (fac < 1.0e-3) cycle
                        delil = one/(rate(ll, ll) - rate(ii, ii))
                        deljl = one/(rate(ll, ll) - rate(jj, jj))
                        if (ii == jj) then
                            secnda = fac*delil*taum*etp(ii)
                            secndb = fac*delil*etp(ll)*delil
                            secndc = -fac*delil*etp(ii)*delil
                            scont = secnda + secndb + secndc
                            secndo = secndo + scont
                        else

                            secnda = fac*etp(ii)*(delil*delij)
                            secndb = fac*etp(ll)*(delil*deljl)
                            secndc = -fac*etp(jj)*(deljl*delij)
                            scont = secnda + secndb + secndc
                            secndo = secndo + scont
                        end if

                        !  --- the second order derivative:

                        fac1 = scont/rate(ii, ll)
                        !forcevector
                        do imuk = 1, imukmx
                            muk = muki(imuk)
                            dint(muk + 1) = dint(muk + 1) + fac1*dorat(muk + 1, ll)
                            dint(muk + 2) = dint(muk + 2) + fac1*dorat(muk + 2, ll)
                            dint(muk + 3) = dint(muk + 3) + fac1*dorat(muk + 3, ll)
                        end do
                        fac2 = scont/rate(jj, ll)
                        !forcevector
                        do jmuk = 1, jmukmx
                            muk = mukj(jmuk)
                            dint(muk + 1) = dint(muk + 1) + fac2*dorat(muk + 1, ll)
                            dint(muk + 2) = dint(muk + 2) + fac2*dorat(muk + 2, ll)
                            dint(muk + 3) = dint(muk + 3) + fac2*dorat(muk + 3, ll)
                        end do
                        muk = -3
                        lmuk = 0
                        do k = 1, nath
                            muk = muk + 3
                            kk = inn(k)
                            if (kk == ll) then
                                lmuk = lmuk + 1
                                mukl(lmuk) = muk
                            end if
                        end do
                        lmukmx = lmuk
                        !forcevector
                        do lmuk = 1, lmukmx
                            muk = mukl(lmuk)
                            dint(muk + 1) = dint(muk + 1) + fac1*dorat(muk + 1, ii) &
                                + fac2*dorat(muk + 1, jj)
                            dint(muk + 2) = dint(muk + 2) + fac1*dorat(muk + 2, ii) &
                                + fac2*dorat(muk + 2, jj)
                            dint(muk + 3) = dint(muk + 3) + fac1*dorat(muk + 3, ii) &
                                + fac2*dorat(muk + 3, jj)
                        end do

                        if (ii == jj) then
                            do muk = 1, 3*nath
                                facil = (ddrat(muk, ii) - ddrat(muk, ll))*delil
                                dint(muk) = dint(muk) &
                                    - secnda*(taum*ddrat(muk, ii) - facil) &
                                    - secndb*(taum*ddrat(muk, ll) - two*facil) &
                                    - secndc*(taum*ddrat(muk, ii) - two*facil)
                            end do
                        else
                            do muk = 1, 3*nath
                                facil = (ddrat(muk, ii) - ddrat(muk, ll))*delil
                                facjl = (ddrat(muk, jj) - ddrat(muk, ll))*deljl
                                facij = (ddrat(muk, ii) - ddrat(muk, jj))*delij
                                dint(muk) = dint(muk) &
                                    - secnda*(taum*ddrat(muk, ii) - facij - facil) &
                                    - secndb*(taum*ddrat(muk, ll) - facil - facjl) &
                                    - secndc*(taum*ddrat(muk, jj) - facjl - facij)
                            end do
                        end if
                    end do

                    !   --- now work on the third-order expression:

                    thirdo = zero
#ifdef THIRDO
                    if (ii == jj) goto 290
                    do ll = 1, natmet
                        if (ll == ii) cycle
                        muk = -3
                        lmuk = 0
                        do k = 1, nath
                            muk = muk + 3
                            kk = inn(k)
                            if (kk == ll) then
                                lmuk = lmuk + 1
                                mukl(lmuk) = muk
                            end if
                        end do
                        lmukmx = lmuk

                        do mm = 1, natmet
                            if (mm == jj .or. mm == ll) cycle
                            if (ll == jj) then
                                if (mm == ii) cycle
                                thir = -rate(jj, mm)*rate(mm, jj)*rate(jj, ii)* &
                                    (etp(ii)*delta(ii, jj)*delta(ii, mm)*delta(ii, jj) + &
                                    etp(mm)*delta(mm, jj)*delta(mm, jj)*delta(mm, ii) + &
                                    etp(jj)*delta(ii, jj)*delta(mm, jj)* &
                                    (delta(ii, jj) + delta(mm, jj)) + &
                                    taum*etp(jj)*delta(jj, ii)*delta(jj, mm))
                            else if (mm == ii) then
                                thir = -rate(ii, ll)*rate(ll, ii)*rate(ii, jj)* &
                                    (etp(jj)*delta(jj, ii)*delta(jj, ll)*delta(jj, ii) + &
                                    etp(ll)*delta(ll, ii)*delta(ll, ii)*delta(ll, jj) + &
                                    etp(ii)*delta(jj, ii)*delta(ll, ii)* &
                                    (delta(jj, ii) + delta(ll, ii)) + &
                                    taum*etp(ii)*delta(ii, jj)*delta(ii, ll))
                            else
                                cycle

                                !         --- to include the "four-body" terms, i.e. terms in the
                                !             third-order expression where ii,jj,ll and mm are all
                                !             different, remove the "go to 122" above, and uncomment
                                !             the following statement:

                            end if
                            thirdo = thirdo + thir

                            !  --- work on principal contribution of third order term to the
                            !          derivatives:

                            muk = -3
                            mmuk = 0
                            do k = 1, nath
                                muk = muk + 3
                                kk = inn(k)
                                if (kk == mm) then
                                    mmuk = mmuk + 1
                                    mukm(mmuk) = muk
                                end if
                            end do
                            mmukmx = mmuk

                            fac = thir/rate(ii, ll)
                            !forcevector
                            do imuk = 1, imukmx
                                muk = muki(imuk)
                                dint(muk + 1) = dint(muk + 1) + fac*dorat(muk + 1, ll)
                                dint(muk + 2) = dint(muk + 2) + fac*dorat(muk + 2, ll)
                                dint(muk + 3) = dint(muk + 3) + fac*dorat(muk + 3, ll)
                            end do

                            fac1 = thir/rate(ii, ll)
                            fac2 = thir/rate(ll, mm)
                            !forcevector
                            do lmuk = 1, lmukmx
                                muk = mukl(lmuk)
                                dint(muk + 1) = dint(muk + 1) + fac1*dorat(muk + 1, ii) &
                                    + fac2*dorat(muk + 1, mm)
                                dint(muk + 2) = dint(muk + 2) + fac1*dorat(muk + 2, ii) &
                                    + fac2*dorat(muk + 2, mm)
                                dint(muk + 3) = dint(muk + 3) + fac1*dorat(muk + 3, ii) &
                                    + fac2*dorat(muk + 3, mm)
                            end do

                            fac1 = thir/rate(ll, mm)
                            fac2 = thir/rate(mm, jj)
                            !forcevector
                            do mmuk = 1, mmukmx
                                muk = mukm(mmuk)
                                dint(muk + 1) = dint(muk + 1) + fac1*dorat(muk + 1, ll) &
                                    + fac2*dorat(muk + 1, jj)
                                dint(muk + 2) = dint(muk + 2) + fac1*dorat(muk + 2, ll) &
                                    + fac2*dorat(muk + 2, jj)
                                dint(muk + 3) = dint(muk + 3) + fac1*dorat(muk + 3, ll) &
                                    + fac2*dorat(muk + 3, jj)
                            end do

                            fac = thir/rate(mm, jj)
                            !forcevector
                            do jmuk = 1, jmukmx
                                muk = mukj(jmuk)
                                dint(muk + 1) = dint(muk + 1) + fac*dorat(muk + 1, mm)
                                dint(muk + 2) = dint(muk + 2) + fac*dorat(muk + 2, mm)
                                dint(muk + 3) = dint(muk + 3) + fac*dorat(muk + 3, mm)
                            end do

                        end do
                    end do
290                 continue
#endif
                    acalc1 = acalc1 + firsto + secndo + thirdo
                    call timer_stop(TIME_REMARC)
                else
                    call timer_start(TIME_DINTEN)

                    !=======================================================================

                    !   Section 2b: exact calculation for the derivatives of the NOEs
                    !       w.r.t. the proton coordinates

                    !=======================================================================

#ifdef ORIGDERIV
                    call dinten(vecs, facj, ii, jj, c, dint)
#else
                    call dinten(amat, ii, jj, ddrat, dorat, dint, taum)
#endif
                    call timer_stop(TIME_DINTEN)

                    ! --- get exact computed intensity

                    do k = 1, natmet
                        acalc1 = acalc1 + vecs(ii, k)*et(k)*vecs(jj, k)
                    end do
                end if
                acalc1 = acalc1*popn(ii)*popn(jj)
                acalc = acalc + acalc1

                !=======================================================================

                !  Section 3:  compute penalty function and update force vector

                !=======================================================================

                ! --- if the peak intensity value is negative, then just go on to the next
                !       peak, and keep summing until a positive weight is found:

                call timer_start(TIME_REMARC)
                if (aexp(imix, ipeak) < zero) then
#ifndef MPI
                    if (iprint /= 0 .and. penalty > pencut) &
                        write (81, 1001) ksub, imix, resat(ihp(imix, ipeak)), &
                        resat(jhp(imix, ipeak)), acalc1, peakid(ipeak)
#endif
                    call timer_stop(TIME_REMARC)
                    cycle
                end if

                ! --- we need to have the scaling factor in the same range as atomic
                !      coordinates.  Hence, an overall scaling factor "oscale"
                !      is applied

                aexpn = oscale*aexp(imix, ipeak)
                arngen = oscale*arange(imix, ipeak)

                !  --- If the invwt1,invwt2 values are not both 1.0, then
                !      set weights equal to inverse of normalized experimental
                !      intensity, with minimum of "invwt1" and maximum of "invwt2":

                if (invwt1 /= one .or. invwt2 /= one) then
                    if (aexpn /= zero) then
                        ainv = one/aexpn
                    else
                        ainv = invwt2
                    end if
                    awt(imix, ipeak) = max(min(ainv, invwt2), invwt1)
                else
                    awt(imix, ipeak) = one
                end if

                !  --- set up "exper" and "calc" arrays for later statistical analysis:
                !       here we ignore "diagonal" peaks greater than 0.5:

                if (acalc < half) then
                    ntot = ntot + 1
                    if (ntot > mtot) then
                        write (6, *) 'Too many peaks for correlation analysis!'
                        call mexit(6, 1)
                    end if
                    calc(ntot) = acalc
                    exper(ntot) = aexpn
                    ipmix(ntot) = imix

                    !  --- assign to intra or inter residue peak:

                    if (resat(ihp(imix, ipeak)) (11:13) == &
                        resat(jhp(imix, ipeak)) (11:13)) then
                        ntota = ntota + 1
                        calca(ntota) = acalc
                        expera(ntota) = aexpn
                    else
                        ntotb = ntotb + 1
                        calcb(ntotb) = acalc
                        experb(ntotb) = aexpn
                    end if
                end if

                weight = wnoesy*awt(imix, ipeak)
                if (ipnlty == 1) then

                    !   --- here the penalty is zero if the calculated volume is in
                    !       the range: aexpn - arngen < acalc < aexpn + arngen
                    !       Outside this range the penalty increases linearly.

                    if (acalc > aexpn + arngen) then
                        penalty = weight*(acalc - (aexpn + arngen))
                    else if (acalc < aexpn - arngen) then
                        penalty = weight*((aexpn - arngen) - acalc)
                    else
                        penalty = zero
                    end if
                else if (ipnlty == 2) then

                    !   --- here the penalty is zero if the calculated volume is in
                    !       the range: aexpn - arngen < acalc < aexpn + arngen
                    !       Outside this range the penalty increases quadratically.

                    if (acalc > aexpn + arngen) then
                        penalty = weight*(acalc - (aexpn + arngen))**2
                    else if (acalc < aexpn - arngen) then
                        penalty = weight*(acalc - (aexpn - arngen))**2
                    else
                        penalty = zero
                    end if
                else if (ipnlty == 3) then

                    !   --- here the penalty is zero if the calculated volume is in
                    !       the range: aexpn - arngen < acalc < aexpn + arngen
                    !       Outside this range the penalty increases quadratically,
                    !       but based on the Sixth-root of the observed intensities.

                    if (acalc > aexpn + arngen) then
                        penalty = weight*(acalc**sixth - (aexpn + arngen)**sixth)**2
                    else if (acalc < aexpn - arngen) then
                        penalty = weight*(acalc**sixth - (aexpn - arngen)**sixth)**2
                    else
                        penalty = zero
                    end if
                else
                    write (6, *) 'Bad value for ipnlty: ', ipnlty
                    call mexit(6, 1)
                end if  ! (ipnlty == 1)
                enoe = enoe + penalty
#ifndef MPI
                if (iprint /= 0 .and. penalty > pencut) then
                    if (acalc /= acalc1) then
                        write (81, 1001) ksub, imix, resat(ihp(imix, ipeak)), &
                            resat(jhp(imix, ipeak)), acalc1, peakid(ipeak)
                        write (81, 1002) ksub, imix, resat(ihp(imix, ipeak)), &
                            resat(jhp(imix, ipeak)), aexpn, acalc, penalty, peakid(ipeak)
                        write (81, '()')
                    else
                        write (81, 1000) ksub, imix, resat(ihp(imix, ipeak)), &
                            resat(jhp(imix, ipeak)), aexpn, acalc, penalty, peakid(ipeak)
                    end if
                end if
#endif

                ! --- set up overall factor for updating derivatives:

                if (ipnlty == 1) then
                    if (acalc > aexpn + arngen) then
                        dpen = weight*popn(ii)*popn(jj)
                    else if (acalc < aexpn - arngen) then
                        dpen = -weight*popn(ii)*popn(jj)
                    else
                        dpen = zero
                    end if
                else if (ipnlty == 2) then
                    if (acalc > aexpn + arngen) then
                        dpen = two*weight*(acalc - (aexpn + arngen))*popn(ii)*popn(jj)
                    else if (acalc < aexpn - arngen) then
                        dpen = two*weight*(acalc - (aexpn - arngen))*popn(ii)*popn(jj)
                    else
                        dpen = zero
                    end if
                else if (ipnlty == 3) then
                    if (acalc > aexpn + arngen) then
                        dpen = third*weight*(acalc**sixth - (aexpn + arngen)**sixth)* &
                            (acalc**fsx)*popn(ii)*popn(jj)
                    else if (acalc < aexpn - arngen) then
                        dpen = third*weight*(acalc**sixth - (aexpn - arngen)**sixth)* &
                            (acalc**fsx)*popn(ii)*popn(jj)
                    else
                        dpen = zero
                    end if
                end if
#ifdef NMODE

                !    ---- accumulate forces relating to normal mode frequencies:

                do n = 1, iscale
                    f(3*natom + n) = f(3*natom + n) - dpen*dint(3*nath + n)
                end do
#endif

                do i = 1, nath
                    i3real = 3*(ihyp(i) - 1)
                    i3 = 3*(i - 1)
                    f(i3real + 1) = f(i3real + 1) - dpen*dint(i3 + 1)
                    f(i3real + 2) = f(i3real + 2) - dpen*dint(i3 + 2)
                    f(i3real + 3) = f(i3real + 3) - dpen*dint(i3 + 3)
                end do

                !  --- zero out calculated intensity and Dint, for preparation
                !        for next peak; this code arises from the possibility of
                !        summing peaks together.

                acalc = zero
                do muk = 1, 3*nath + 1
                    dint(muk) = zero
                end do
                call timer_stop(TIME_REMARC)
            end do
        end do

        return
#ifndef MPI
1000    format('noe:', i4, i3, 2x, a13, 3x, a13, 3x, 3f10.5, i5)
1001    format('ovl:', i4, i3, 2x, a13, 3x, a13, 13x, f10.5, 10x, i5)
1002    format('tot:', i4, i3, 2x, a13, 3x, a13, 3x, 3f10.5, i5)
#endif
    end subroutine remarc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine remhet here]
    subroutine remhet(f, x, ksub, newf, amass)
        use constants, only : zero, third, sixth
        implicit none
        _REAL_ f(*), x(*)
        integer ksub
        logical newf
        _REAL_ amass(*)
#ifdef NMODE

        ! Subroutine RElaxation MAtrix for HETeronuclear relaxation (2-spin!)

        !     interprets input "aexp" from noecalc as a "gamma" correction
        !        value, and sets up forces to minimize the differences between
        !        calc. and exp. values.

#  include "nmr.h"
#  include "memory.h"

        integer imix, ipeak
        _REAL_ acalc, aexpn, arngen
        _REAL_ weight, penalty, dpen
        integer n
        _REAL_ fsx
        parameter(fsx=-5.d0/6.d0)

        imix = 1
        do ipeak = 1, npeak(1)
            call corf(x, ihp(1, ipeak), jhp(1, ipeak), ihp(1, ipeak), jhp(1, ipeak), &
                newf, amass)
            newf = .false.
            acalc = gamma_nmr

            ! --- compute penalty function and update force vector

            aexpn = aexp(imix, ipeak)
            arngen = arange(imix, ipeak)

            weight = wnoesy*awt(imix, ipeak)
            if (ipnlty == 1) then

                !   --- here the penalty is zero if the calculated volume is in
                !       the range: aexpn - arngen < acalc < aexpn + arngen
                !       Outside this range the penalty increases linearly.

                if (acalc > aexpn + arngen) then
                    penalty = weight*(acalc - (aexpn + arngen))
                else if (acalc < aexpn - arngen) then
                    penalty = weight*((aexpn - arngen) - acalc)
                else
                    penalty = zero
                end if
            else if (ipnlty == 2) then

                !   --- here the penalty is zero if the calculated volume is in
                !       the range: aexpn - arngen < acalc < aexpn + arngen
                !       Outside this range the penalty increases quadratically.

                if (acalc > aexpn + arngen) then
                    penalty = weight*(acalc - (aexpn + arngen))**2
                else if (acalc < aexpn - arngen) then
                    penalty = weight*(acalc - (aexpn - arngen))**2
                else
                    penalty = zero
                end if
            else if (ipnlty == 3) then

                !   --- here the penalty is zero if the calculated volume is in
                !       the range: aexpn - arngen < acalc < aexpn + arngen
                !       Outside this range the penalty increases quadratically,
                !       but based on the Sixth-root of the observed intensities.

                if (acalc > aexpn + arngen) then
                    penalty = weight*(acalc**sixth - (aexpn + arngen)**sixth)**2
                else if (acalc < aexpn - arngen) then
                    penalty = weight*(acalc**sixth - (aexpn - arngen)**sixth)**2
                else
                    penalty = zero
                end if
            else
                write (6, *) 'Bad value for ipnlty: ', ipnlty
                call mexit(6, 1)
            end if  ! (ipnlty == 1)
            enoe = enoe + penalty
#ifndef MPI
1000        format('het:', i3, '  1  ', a13, 3x, a13, 3x, 3f10.5)
            if (iprint /= 0 .and. penalty > pencut) then
                write (81, 1000) ksub, &
                    resat(ihp(1, ipeak)), &
                    resat(jhp(1, ipeak)), aexpn, acalc, penalty
            end if
#endif

            ! --- set up overall factor for updating derivatives:

            if (ipnlty == 1) then
                if (acalc > aexpn + arngen) then
                    dpen = weight
                else if (acalc < aexpn - arngen) then
                    dpen = -weight
                else
                    dpen = zero
                end if
            else if (ipnlty == 2) then
                if (acalc > aexpn + arngen) then
                    dpen = 2.*weight*(acalc - (aexpn + arngen))
                else if (acalc < aexpn - arngen) then
                    dpen = 2.*weight*(acalc - (aexpn - arngen))
                else
                    dpen = zero
                end if
            else if (ipnlty == 3) then
                if (acalc > aexpn + arngen) then
                    dpen = third*weight*(acalc**sixth - (aexpn + arngen)**sixth)* &
                        (acalc**fsx)
                else if (acalc < aexpn - arngen) then
                    dpen = third*weight*(acalc**sixth - (aexpn - arngen)**sixth)* &
                        (acalc**fsx)
                else
                    dpen = zero
                end if
            end if

            !    ---- accumulate forces relating to normal mode frequencies:

            do n = 1, iscale
                f(3*natom + n) = f(3*natom + n) - dpen*dgamma(n + 3)
            end do

        end do

#endif
        return
    end subroutine remhet
#ifndef ORIGDERIV

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine amatg here]
    subroutine amatg(vecs, eig, taum, amat, iamat)

        !   --- compute the intensity matrix at Gaussian points needed for
        !         the integral form of the exact derivative
        use constants, only : zero
        implicit none
#  include "nmr.h"
        _REAL_ vecs(ma, ma), eig(ma), amat(ma, ma, 5), s(5)
        _REAL_ et(ma, 5)
        integer iamat(ma)

        integer i, j, k
        integer ig, igmax
        _REAL_ taum, taus
        _REAL_ v
        save s, igmax
        data s/0.0469100770d0, 0.2307653449d0, 0.5d0, 0.7692346551d0, &
            0.9530899230d0/
        data igmax/5/

        do ig = 1, igmax
            taus = s(ig)*taum
            do k = 1, natmet
                et(k, ig) = exp(-taus*eig(k))
            end do
        end do
        do i = 1, natmet
            if (iamat(i) == 0) cycle
            do j = 1, natmet
                amat(i, j, 1) = zero
                amat(i, j, 2) = zero
                amat(i, j, 3) = zero
                amat(i, j, 4) = zero
                amat(i, j, 5) = zero
                do k = 1, natmet
                    v = vecs(i, k)*vecs(j, k)
                    amat(i, j, 1) = amat(i, j, 1) + v*et(k, 1)
                    amat(i, j, 2) = amat(i, j, 2) + v*et(k, 2)
                    amat(i, j, 3) = amat(i, j, 3) + v*et(k, 3)
                    amat(i, j, 4) = amat(i, j, 4) + v*et(k, 4)
                    amat(i, j, 5) = amat(i, j, 5) + v*et(k, 5)
                end do
            end do
        end do
        return
    end subroutine amatg
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine noecalc here]
    subroutine noecalc(x, f, xx, ix)

        !  Subroutine NOEsy intensity CALCulation:

        !  --- computes relaxation matrix approximation to noe
        !      intensities and their derivatives.  Places energy
        !      into enoe and the negative of the gradient is used
        !      to update f.

        !   Inputs:

        !         x        --  coordinates
        !         f        --  forces
        !         xx,ix    -- main storage arrays

        !   Outputs:

        !         enoe     -- contains the penaly energy for noe-related restraints
        !         f        -- is updated with new derivatives

        !   Overview:

        !      This is mostly a "driver" routine.  It reads input for the
        !      submolecules, in namelist "&noeexp", then calls indexn, caldis
        !      and remarc to do the brunt of the work.  The computation
        !      and printing of various statistics is also handled by this routine.

        !      There is a major additional section, indicated by #ifdef NMODE,
        !      that allows vibrational modes to be read in and used to compute
        !      motional correction factors for the NOEs.  This is still
        !      experimental code, and is *not* currently described in the
        !      users manual.

        !-----------------------------------------------------------------------------
        use file_io_dat
        implicit none
#ifdef NMODE
        integer:: line, n
        _REAL_ :: consq, dev, dpen, earg, edev, penalty
#endif
        integer:: i, im, imet, imix, indxs, ip, ix, k, ksub, l, nuse
        _REAL_ :: f, prob, r, rms, tauc, ten, x, xhyd, xmet, xrfac6, xrfact, xx, z
        logical newf

#  include "nmr.h"
#  include "memory.h"
#  include "md.h"
#  include "def_time.h"
#ifdef NMODE

        real freq4, xw(3*matom)
        _REAL_ kt
        parameter(kt=0.6d0)
        parameter(consq=2.39805d-3)
        !                  ----CONSQ = hc/2kT in cm, with T=300K
        !                       (use for quantum, Bose statistics)
#endif
#ifdef MPI
#  include "parallel.h"
#endif
#  include "extra.h"

        _REAL_ frespa
        dimension xx(*), ix(*)
        dimension xmet(isubr), imet(isubi)
        equivalence(nath, imet(1))
        equivalence(tau(1), xmet(1))
        dimension x(*), f(*)
        dimension xhyd(3*ma)

        ten = 10.d0
#ifndef MPI

        ! --- set up printer output:

1000    format(1x, 'Wnoesy = ', f8.3)
1010    format(1x, 79('-')/4x, 'sub mix    Proton 1', 8x, 'Proton 2', 6x, &
            '  Exp.      Calc.     Penalty   ID'/1x, 79('-'))
        if (iprint /= 0 .and. master) then
            rewind (81)
            write (81, 1000) wnoesy
            write (81, 1010)
            write (6, 1000) wnoesy
        end if
#endif

        ! --- zero out the forces arising from the "extra" degrees of freedom:

        do im = 1, iscale
            f(3*natom + im) = 0.0d0
        end do

        ! --- use respa impulses for this force every noeskp steps:

        enoe = 0.0d0
        if (mod(irespa, noeskp) /= 0) return
        frespa = noeskp

        ! --- read in relaxation matrix parameters, desired mixing times, and
        !      experimental intensities

        call timer_start(TIME_NOECALC1)
        ntot = 0
        ntota = 0
        ntotb = 0
#ifdef MPI
        ksub = mytaskid + 1 - numtasks
#else
        ksub = 0
#endif
        enoe = 0.0d0
        do i = 1, 3*natom + iscale
            xx(l110 - 1 + i) = 0.0d0
        end do

140     continue

        ! --- grab submolecule info from where it was stored in noeread:

#ifdef MPI
        ksub = ksub + numtasks
#else
        ksub = ksub + 1
#endif
        if (ksub > maxsub) then
            call timer_stop(TIME_NOECALC1)
            goto 280
        end if
        indxs = l105 + (ksub - 1)*isubr - 1
        do i = 1, isubr
            xmet(i) = xx(indxs + i)
        end do
        indxs = i65 + (ksub - 1)*isubi - 1
        do i = 1, isubi
            imet(i) = ix(indxs + i)
        end do

        ! --- initialize coordinates of the sub-molecule:

        if (ihet == 0) then
            l = 0
            do i = 1, nath
                k = 3*(ihyp(i) - 1)
                xhyd(l + 1) = x(k + 1)
                xhyd(l + 2) = x(k + 2)
                xhyd(l + 3) = x(k + 3)
                l = l + 3
            end do
        end if

        ! --- overall rotational correlation time in h2o is taurot; here
        !      we hard-wire in the assumption that this is increased by
        !      23% in d2o.

        tauc = taurot*2.d0
        if (id2o == 1) tauc = 1.23d0*tauc

        ! --- in calrate, the relative tumbling time for each proton pair
        !      is calculated as taui*tauj/(taui+tauj).  This gives a constant
        !      in the current code, but is included for (future?) models in
        !      which the effective correlation time is not the same for all
        !      atoms.  Here we set up the tau(i) array as twice the rotational
        !      correlation time.

        do im = 1, nath
            tau(im) = tauc
        end do

        !----calculate  penalty function and its gradient:  caldis computes the
        !       distance-dependent parts of the rate matrix, and their derivatives,
        !       and remarc completes the work (q.v.).  Subroutine remhet is
        !       called for vibrational calculations of Lipari-Szabo order
        !       parameters for heteronuclear relaxation.  Logical variable "newf"
        !       is true if the frequencies are changed since the last time
        !       order parameters were evaluated.

        call timer_stop(TIME_NOECALC1)
        if (ksub == 1) then
            newf = .true.
        else
            newf = .false.
        end if
        if (ihet == 1) then
            call remhet(xx(l110), xx(lcrd), ksub, newf, xx(lmass))
        else
            call timer_start(TIME_CALDIS)
#ifdef NMODE
            call caldis(xhyd, xx(l115), xx(l120), newf, xx(lmass))
#else
            call caldis(xhyd, xx(l115), xx(l120))
#endif
            call timer_stop(TIME_CALDIS)
            call remarc(xx(l115), xx(l120), xx(l110), xhyd, ksub, &
                xx(l125), xx(l130), xx(l135), xx(l140), xx(l145))
        end if

        !----redo calc. for another submolecule

        call timer_start(TIME_NOECALC1)
        goto 140

        !---if all sub mol. are done, compute and print some statistics,
        !     update the forces and return:

280     continue
#ifndef MPI

        ! --- get linear correlation coefficient between cacl. and obs.

1050    format(/50x, 'Total NOESY penalty: ', f7.2/)
1060    format(/21x, '   #   Pearson r  rms error   R1      Rr      ', &
            /1x, 75('-'))
        call timer_start(TIME_NOECALC2)
        if (iprint /= 0) then
            write (81, 1050) enoe
            if (ntot < 3) goto 350
            write (81, 1060)
            call pearsn(exper, calc, ntot, r, prob, z, rms, xrfact, xrfac6, nuse)
            write (81, '(a21,i5,4f10.5)') &
                'Full  Correlation: = ', nuse, r, rms, xrfact, xrfac6
            call pearsn(expera, calca, ntota, r, prob, z, rms, xrfact, xrfac6, nuse)
            write (81, '(a21,i5,4f10.5)') &
                'Intra Correlation: = ', nuse, r, rms, xrfact, xrfac6
            call pearsn(experb, calcb, ntotb, r, prob, z, rms, xrfact, xrfac6, nuse)
            write (81, '(a21,i5,4f10.5)') &
                'Inter Correlation: = ', nuse, r, rms, xrfact, xrfac6
            ntotb = 0
            do ip = 1, ntot
                if (exper(ip) < 0.1) cycle
                ntotb = ntotb + 1
                experb(ntotb) = exper(ip)
                calcb(ntotb) = calc(ip)
            end do
            call pearsn(experb, calcb, ntotb, r, prob, z, rms, xrfact, xrfac6, nuse)
            write (81, '(a21,i5,4f10.5)') &
                '>0.1  Correlation: = ', nuse, r, rms, xrfact, xrfac6

            ntotb = 0
            do ip = 1, ntot
                if (exper(ip) < 0.05 .or. exper(ip) > 0.1) cycle
                ntotb = ntotb + 1
                experb(ntotb) = exper(ip)
                calcb(ntotb) = calc(ip)
            end do
            call pearsn(experb, calcb, ntotb, r, prob, z, rms, xrfact, xrfac6, nuse)
            write (81, '(a21,i5,4f10.5)') &
                '0.05-0.1 Correl  : = ', nuse, r, rms, xrfact, xrfac6

            ntotb = 0
            do ip = 1, ntot
                if (exper(ip) > 0.05 .or. exper(ip) < 0.025) cycle
                ntotb = ntotb + 1
                experb(ntotb) = exper(ip)
                calcb(ntotb) = calc(ip)
            end do
            call pearsn(experb, calcb, ntotb, r, prob, z, rms, xrfact, xrfac6, nuse)
            write (81, '(a21,i5,4f10.5)') &
                '0.025-0.05 Correl: = ', nuse, r, rms, xrfact, xrfac6

            ntotb = 0
            do ip = 1, ntot
                if (exper(ip) > 0.025) cycle
                ntotb = ntotb + 1
                experb(ntotb) = exper(ip)
                calcb(ntotb) = calc(ip)
            end do
            call pearsn(experb, calcb, ntotb, r, prob, z, rms, xrfact, xrfac6, nuse)
            write (81, '(a21,i5,4f10.5)') &
                '<0.025 Correlation:= ', nuse, r, rms, xrfact, xrfac6

            !  --statistics for each mixing time:

            do imix = 1, nummt
                ntotb = 0
                do ip = 1, ntot
                    if (ipmix(ip) /= imix) cycle
                    ntotb = ntotb + 1
                    experb(ntotb) = exper(ip)
                    calcb(ntotb) = calc(ip)
                end do
                call pearsn(experb, calcb, ntotb, r, prob, z, rms, xrfact, xrfac6, nuse)
                write (81, '(i10,a11,i5,4f10.5)') &
                    imix, ' Correl: = ', nuse, r, rms, xrfact, xrfac6
            end do
350         continue
        end if
        call timer_stop(TIME_NOECALC2)
#endif

        ! --- update full derivatives

#ifdef NMODE
        do i = 1, 3*natom
            f(i) = 0.0
        end do
        do i = 3*natom + 1, 3*natom + iscale
            f(i) = f(i) + frespa*xx(l110 - 1 + i)
        end do

1080    format(' ', 78('-'))
        if (iprint /= 0) then
            write (6, 1080)
            write (6, *) 'Heteronuclear normal mode analysis:'
            write (6, *)
            if (iscale > 1) then
                write (6, *) 'frequencies:'
                write (6, '(6f12.6)') (x(3*natom + i), i=1, iscale - 1)
            end if
            write (6, 1080)
            rewind (81)
390         read (81, '(a80)', end=400) line
            write (6, '(a80)') line
            goto 390
400         rewind (81)
            write (6, 1080)
        end if

        !  --- now set up penalties to keep the effective frequencies close
        !         to their "true" values:  right now assumes first six
        !         frequencies are zero.

        !      set the penalty function to:

        !      penalty = (xdev/freq_true)(freq_ref - freq_true)**2
        !               + exp(-6*freq_ref)

        !      where freq_exp is the "true" frequency, and freq_ref is the
        !          current value of the "refined" frequency.  The first term
        !          is a usual quadratic penalty (with weight equal to the
        !          inverse of freq_true,) and the second term should serve to
        !          keep all frequencies positive.

        !     xdev = 0.05  is the default

        edev = 0.0
        do n = 1, iscale - 1
            dev = x(3*natom + n) - freq(n)
            earg = min(-6.0*x(3*natom + n), ten)
            earg = max(-ten - ten, earg)
            penalty = xdev*dev**2/freq(n) + exp(earg)
            enoe = enoe + penalty
            edev = edev + penalty
            dpen = 2.0*xdev*dev/freq(n) - 6.*exp(earg)
            f(3*natom + n) = f(3*natom + n) - dpen
        end do
        if (iprint /= 0) then
            write (6, '(a23,f12.5)') 'global scaling factor: ', &
                x(3*natom + iscale)
            write (6, '(a29,e14.6)') 'frequency deviation penalty: ', edev
        end if
#else
        do i = 1, 3*natom + iscale
            f(i) = f(i) + frespa*xx(l110 - 1 + i)
        end do
#endif

        return

#ifdef DEBUG_NMR
1030    format(1x, 79('=')/'Data for submolecule', i4, ': id2o = ', &
            i1, ', oscale = ', e12.5, ', taumet = ', e12.5, / &
            36x, 'omega = ', e12.5, 'taurot =', e12.5/)
#endif

    end subroutine noecalc

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine noeread here]
    subroutine noeread(xx, ix, ih)

        use file_io_dat
        implicit none
#ifdef NMODE
        integer:: im, ivec, ivform, k, l, lmax, lmin, n3, nf
        _REAL_ :: consq, omecut, omescl, xkappa
#endif
        integer:: i, ifind, iin, imet, imix, imol, indxs, ix, j, ksub, numpks
        _REAL_ :: xmet, xx
#  include "nmr.h"
#  include "memory.h"
#  include "md.h"
#  include "parallel.h"

        dimension xx(*), ix(*)
        character(len=4) ih(*)
        dimension xmet(isubr), imet(isubi)
        equivalence(nath, imet(1))
        equivalence(tau(1), xmet(1))
        logical hetero
        character(len=80) line
        namelist /noeexp/ npeak, emix, ihp, jhp, aexp, awt, arange, id2o, oscale, &
            taumet, omega, invwt1, invwt2, hetero, iroesy, taurot, peakid
#ifdef NMODE
        logical hsfull
        integer g98_vecs, jaguar_vecs
        character(len=4) star
        character(len=2) spacer
        real freq4
        _REAL_ kt
        namelist /nmode/ nvect, hsfull, bose, xdev, xkappa, ivform, omegax, &
            iusev, jaguar_vecs, g98_vecs, vtemp, per_mode, nmsnap
#endif

        maxsub = 0
        imol = 0
        ksub = 0
        id2o = 0
        iroesy = 0
        hetero = .false.
        taurot = 1.0d0
        oscale = 1.0d0
        invwt1 = 1.0d0
        invwt2 = 1.0d0
        if (iredir(4) /= 0) then
            call amopen(35, redir(4) (1:iredir(4)), 'O', 'F', 'R')
            iin = 35
            write (6, 1020) redir(4) (1:iredir(4))
        else
            return
        end if

        !  --- read and echo title from noeexp file:

        write (6, *) 'Here are comments from the NOEsy input file:'
42      read (iin, '(a)', end=140) line
        if (line(1:1) == '#') then
            write (6, *) line
            goto 42
        end if
        backspace (iin)
        write (6, *)

#ifdef NMODE

        !    --- first, read in namelist nmode, after setting defaults:

        nvect = 99999
        hsfull = .true.
        bose = .true.
        xkappa = 1.0
        omegax = 1000.
        xdev = 0.05
        ivform = 1
        jaguar_vecs = 0
        g98_vecs = 0
        per_mode = .false.
        nmsnap = 0
        vtemp = 300.

        call nmlsrc('nmode', iin, ifind)
        if (ifind == 0) then
            write (6, *) 'Unable to fine namelist "nmode"'
            call mexit(6, 1)
        end if
        read (iin, nml=nmode, err=460)
        write (6, *) 'read nmode: ', nvect, hsfull, bose, xdev, xkappa, ivform, &
            vtemp, omegax
        if (hsfull) then
            ihsful = 1
        else
            ihsful = 0
        end if
        kt = vtemp*0.002
        !                  ----CONSQ = hc/2kT in cm
        !                       (use for quantum, Bose statistics)
        consq = 0.71942/vtemp
        xx(lcrd - 1 + 3*natom + iscale) = xkappa

        !   --- read in the normal mode frequencies and vectors:

        nf = 60

        if (jaguar_vecs > 0) then
            call amopen(nf, vecs, 'O', 'F', 'R')
            do lmin = 1, nvect, 7
                lmax = min(nvect, lmin + 6)
                read (nf, '(13x,7f9.2)') (freq(l), l=lmin, lmax)
                write (6, '(13x,7f9.2)') (freq(l), l=lmin, lmax)
                do k = 1, 3*natom
                    read (nf, '(13x,7f9.5)') (vect(k, l), l=lmin, lmax)
                end do
                read (nf, '(a2)') spacer
            end do
        else if (g98_vecs > 0) then
            call amopen(nf, vecs, 'O', 'F', 'R')
            do lmin = 1, nvect, 5
                read (nf, '(a2)') spacer
                read (nf, '(a2)') spacer
                lmax = min(nvect, lmin + 4)
                read (nf, '(23x,5f10.4)') (freq(l), l=lmin, lmax)
                write (6, '(23x,5f10.4)') (freq(l), l=lmin, lmax)
                read (nf, '(a2)') spacer
                read (nf, '(a2)') spacer
                read (nf, '(a2)') spacer
                read (nf, '(a2)') spacer
                read (nf, '(a2)') spacer
                read (nf, '(a2)') spacer
                do k = 1, 3*natom
                    read (nf, '(23x,5f10.5)') (vect(k, l), l=lmin, lmax)
                end do
            end do
        else
            if (ivform == 0) then
                call amopen(nf, vecs, 'O', 'U', 'R')
                read (nf, err=470) n3
            else
                call amopen(nf, vecs, 'O', 'F', 'R')
                read (nf, '(a40)', err=470) title
                read (nf, '(i5)', err=470) n3
            end if
            if (n3 /= natom*3) then
                write (6, *) 'number of atoms wrong in vecs: ', n3, natom*3
                call mexit(6, 1)
            end if

            !  --- do not need coords, put temporarily into vect(,1):

            if (ivform == 0) then
                read (nf, err=470)
            else
                read (nf, '(7f11.5)', err=470) (vect(j, 1), j=1, n3)
            end if
            do l = 1, nvect
                if (ivform == 0) then
                    read (nf, end=70, err=470) ivec, freq4
                    freq(l) = freq4
                    read (nf, err=470) (vect(k, l), k=1, n3)
                else
                    read (nf, '(a4)', end=70, err=470) star
                    read (nf, '(i5,f12.6)', err=470) ivec, freq(l)
                    read (nf, '(7f11.5)', err=470) (vect(k, l), k=1, n3)
                end if
                if (freq(l) < 0.5) freq(l) = 10000.
            end do
            goto 80
70          nvect = l - 1
        end if
80      write (6, *) 'Found', nvect, ' eigenvectors in vecs file'
        close (nf)

        !  --- hard-wired: scale frequencies as rafael suggested:

        if (freqe /= 'dummy') then
            call amopen(nf, freqe, 'O', 'F', 'R')
            read (nf, *) omecut, omescl
            do i = 1, nvect
                if (freq(i) < omecut) freq(i) = omescl*freq(i)
            end do
        end if

        !  --- initially set scaling factors to "true" frequencies,

        do im = 1, iscale - 1
            xx(lcrd - 1 + 3*natom + im) = freq(im)
        end do

        !  --- kludge for now: zero out velocities of all real atoms:

        do i = 1, 3*natom
            xx(lvel - 1 + i) = 0.0
        end do
#endif

140     numpks = 0

        !   --zero out the arange and npeak arrays:

        do i = 1, mxtau
            npeak(i) = 0
            do j = 1, mxp
                arange(i, j) = 0.0d0
            end do
        end do
        call nmlsrc('noeexp', iin, ifind)
        if (ifind == 0) goto 280
        read (iin, nml=noeexp, end=280, err=450)
        imol = imol + 1
        if (hetero) then
            ihet = 1
        else
            ihet = 0
        end if
        do imix = 1, mxtau
            if (npeak(imix) < 0) goto 190
            if (npeak(imix) > mxp) then
                write (6, *) 'Npeak is too big.'
                call mexit(6, 1)
            end if

            do i = 1, npeak(imix)
                numpks = numpks + 1
            end do
        end do

        ! ------ if we fall out of the do-loop, set the number of mixing times
        !         to MXTAU:

        nummt = mxtau
        goto 200

        ! ------come here when negative peak number is encountered;
        !         total number of mixing times is nummt

190     nummt = imix - 1
200     write (6, 1040) nummt, numpks
        if (nummt <= 0) goto 280

        ! --- done reading experimental info; get various information about
        !       protons in this submolecule into the appropriate arrays:

        if (.not. hetero) call indexn(ix, ih, iin)

        ! --- Store info on this submolecule into appropriate locations in the
        !       xx and ix arrays:

        ksub = ksub + 1
        if (ksub > mxsub) then
            write (6, *) 'Too many sub-molecules!'
            call mexit(6, 1)
        end if
        indxs = l105 + (ksub - 1)*isubr - 1
        do i = 1, isubr
            xx(indxs + i) = xmet(i)
        end do
        indxs = i65 + (ksub - 1)*isubi - 1
        do i = 1, isubi
            ix(indxs + i) = imet(i)
        end do

        !   ---cycle back for another sub-molecule:

        goto 140

280     maxsub = ksub
        return

        ! namelist error on read:

450     write (6, *) 'Namelist reports error in reading noeexp'
        write (6, *) '-- Subscript out of range implies dimensioning problem'
        write (6, *) '-- (see nmr.h)'
        call mexit(6, 1)
#ifdef NMODE
460     write (6, *) 'Namelist reports error in reading nmode'
        call mexit(6, 1)
470     write (6, *) 'Error reported in reading frequency file'
        call mexit(6, 1)
#endif

1020    format(' Noesy volumes will be read from file: ', a)
1040    format(1x, 'Read ', i5, ' mixing times with ', i5, &
            ' total peaks.')

    end subroutine noeread

end module relax_mat
