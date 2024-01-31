! <compile=optimized>

module ew_bspline

#include "copyright.h"
#include "dprec.fh"

    integer, parameter :: bspl_maxorder = 20, BSPL_MAX_ATOMS = 50000, &
        BSPL_MAX_LAYERS = 200, maxnfft = 500

    integer nbsplist(BSPL_MAX_ATOMS), nbspstrt(BSPL_MAX_LAYERS)
    integer nstart, nend, nremain, ndelt, nst

    integer :: kbot, ktop, jbot, jtop
    integer :: igood
    integer :: my_ks(BSPL_MAX_ATOMS)
    integer :: nxtra, pe_xtra(BSPL_MAX_LAYERS), lget(BSPL_MAX_LAYERS)

    private
    public get_grid_weights, kbot, ktop, jbot, jtop, load_prefacs, &
        fill_bspline_0, fill_bspline_1, fill_bspline_2

!===========================================================================
contains
!===========================================================================

!--------------------------------------------------------
!       GET_GRID_WEIGHTS
!-----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_grid_weights here]
    subroutine get_grid_weights( &
        natom, crd, recip, nfft1, nfft2, nfft3, &
        fr1, fr2, fr3, order, theta1, theta2, theta3, &
        dtheta1, dtheta2, dtheta3, &
        d2theta1, d2theta2, d2theta3, &
        my_cg, nmine, nderiv, num_ks)
        implicit none
        integer natom, nfft1, nfft2, nfft3, order
        _REAL_ crd(3, natom), recip(3, 3)
        _REAL_ fr1(*), fr2(*), fr3(*)
        _REAL_ theta1(order, *), theta2(order, *), &
            theta3(order, *)
        _REAL_ dtheta1(order, *), dtheta2(order, *), &
            dtheta3(order, *)
        _REAL_ d2theta1(order, *), d2theta2(order, *), &
            d2theta3(order, *)
        integer my_cg(*), nmine, nderiv, num_ks
        _REAL_ fr3n, fr2n, fr1n, w, w1, w2, w3
        _REAL_ anint

#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"

        integer kbot0, ktop1, ido
#endif

        integer n, k00
        integer imine

        !     sanity check:

        if (order < 2 + nderiv) then
            write (6, *) 'too many B-spline derivs for order! '
            call mexit(6, 1)
        end if
        if (nderiv > 2) then
            write (6, *) 'More than 2 derivs of B-splines not implemented yet!'
            call mexit(6, 1)
        end if
#ifdef MPI
        kbot0 = mxystart(mytaskid)
        kbot = kbot0 + 1
        ktop = kbot0 + mxyslabs
        ktop1 = ktop + order - 2
        imine = 0
        !  ------------------------------------------------
        !     First filter the atoms and make a list (my_cg)
        !          of atoms needed for generating this part of
        !          the grid
        !     If there are too many atoms for the allotted space,
        !          return when nmine exceeds num_ks

        do n = 1, natom
            w = crd(1, n)*recip(1, 3) &
                + crd(2, n)*recip(2, 3) + crd(3, n)*recip(3, 3)
            fr3n = nfft3*(w - (anint(w) - 0.5d0))
            k00 = int(fr3n)
            ! code for filtering atoms. In single proc mode, do all atoms
            ido = 0
            if (ktop1 >= nfft3) then
                if (k00 >= kbot0 .or. k00 <= ktop1 - nfft3) ido = 1
            else
                if (k00 >= kbot0 .and. k00 <= ktop1) ido = 1
            end if
            if (ido == 1) then
                imine = imine + 1
                !           ----- do not fill my_cg past num_ks ------
                if (imine <= num_ks) my_cg(imine) = n
            end if
        end do
        nmine = imine
        !   ----- ERROR condition met --------------------
        if (imine > num_ks) return
        !   ----------------------------------------------
#else
        nmine = natom
#endif

        imine = 0
        do imine = 1, nmine
#ifdef MPI
            n = my_cg(imine)
#else
            n = imine
#endif
            w = crd(1, n)*recip(1, 3) &
                + crd(2, n)*recip(2, 3) + crd(3, n)*recip(3, 3)
            fr3n = nfft3*(w - (anint(w) - 0.5d0))
            w = crd(1, n)*recip(1, 1) &
                + crd(2, n)*recip(2, 1) + crd(3, n)*recip(3, 1)
            fr1n = nfft1*(w - (anint(w) - 0.5d0))
            w = crd(1, n)*recip(1, 2) &
                + crd(2, n)*recip(2, 2) + crd(3, n)*recip(3, 2)
            fr2n = nfft2*(w - (anint(w) - 0.5d0))
            fr1(imine) = fr1n
            fr2(imine) = fr2n
            fr3(imine) = fr3n
            w1 = fr1n - int(fr1n)
            w2 = fr2n - int(fr2n)
            w3 = fr3n - int(fr3n)
            if (nderiv == 0) then
                call fill_bspline_0(w1, order, theta1(1, imine))
                call fill_bspline_0(w2, order, theta2(1, imine))
                call fill_bspline_0(w3, order, theta3(1, imine))
            else if (nderiv == 1) then
                call fill_bspline_1(w1, order, theta1(1, imine), &
                    dtheta1(1, imine))
                call fill_bspline_1(w2, order, theta2(1, imine), &
                    dtheta2(1, imine))
                call fill_bspline_1(w3, order, theta3(1, imine), &
                    dtheta3(1, imine))
            else if (nderiv == 2) then
                call fill_bspline_2(w1, order, theta1(1, imine), &
                    dtheta1(1, imine), d2theta1(1, imine))
                call fill_bspline_2(w2, order, theta2(1, imine), &
                    dtheta2(1, imine), d2theta2(1, imine))
                call fill_bspline_2(w3, order, theta3(1, imine), &
                    dtheta3(1, imine), d2theta3(1, imine))
            end if
        end do  !  imine = 1,nmine
        return
    end subroutine get_grid_weights
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_bspline_0 here]
    subroutine fill_bspline_0(w, order, array)
        !---------- use standard B-spline recursions: see doc file
        !----- w is fraction between 0 and 1;
        !----- order is the order of interpolation
        !----- array is the array of shifted bsplines
        ! Using notation from Essmann et al; w = u-[u] and
        ! array(j) = M_n(w + order - j)  where n is order & w = u - [u];
        implicit none
        integer order
        _REAL_ w, array(order), div

        integer k
        ! init order 2
        array(2) = w
        array(1) = 1.d0 - w
        if (order == 2) return
        ! one pass to order 3
        array(3) = 0.5d0*w*array(2)
        array(2) = 0.5d0*((w + 1.d0)*array(1) + (2.d0 - w)*array(2))
        array(1) = 0.5d0*(1.d0 - w)*array(1)
        if (order == 3) return
        ! one pass to order 4
        div = 1.d0/3.d0
        array(4) = div*w*array(3)
        array(3) = div*((w + 1.d0)*array(2) + (3.d0 - w)*array(3))
        array(2) = div*((w + 2.d0)*array(1) + (2.d0 - w)*array(2))
        array(1) = div*(1.d0 - w)*array(1)
        if (order == 4) return
        do k = 5, order
            call one_pass_bspline(array, w, k)
        end do
        return
    end subroutine fill_bspline_0
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_bspline_1 here]
    subroutine fill_bspline_1(w, order, array, darray)
        !---------- use standard B-spline recursions: see doc file
        !----- w is fraction between 0 and 1;
        !----- order is the order of interpolation
        !----- array is the array of shifted bsplines
        !----- darray is array of derivs
        ! Using notation from Essmann et al; w = u-[u] and
        ! array(j) = M_n(w + order - j)  where n is order & w = u - [u];
        implicit none
        integer order
        _REAL_ w, array(order), &
            darray(order), div

        integer k
        ! init order 2
        array(2) = w
        array(1) = 1.d0 - w
        if (order == 3) then
            ! deriv
            darray(1) = -array(1)
            darray(2) = array(1) - array(2)
            darray(3) = array(2)
            ! one pass to order 3
            array(3) = 0.5d0*w*array(2)
            array(2) = 0.5d0*((w + 1.d0)*array(1) + (2.d0 - w)*array(2))
            array(1) = 0.5d0*(1.d0 - w)*array(1)
            return
        else if (order == 4) then
            ! one pass to order 3
            array(3) = 0.5d0*w*array(2)
            array(2) = 0.5d0*((w + 1.d0)*array(1) + (2.d0 - w)*array(2))
            array(1) = 0.5d0*(1.d0 - w)*array(1)
            ! diff to get darray
            darray(1) = -array(1)
            darray(2) = array(1) - array(2)
            darray(3) = array(2) - array(3)
            darray(4) = array(3)
            ! one final pass to order 4
            div = 1.d0/3.d0
            array(4) = div*w*array(3)
            array(3) = div*((w + 1.d0)*array(2) + (3.d0 - w)*array(3))
            array(2) = div*((w + 2.d0)*array(1) + (2.d0 - w)*array(2))
            array(1) = div*(1.d0 - w)*array(1)
            return
        else
            ! general order case
            ! one pass to order 3
            array(3) = 0.5d0*w*array(2)
            array(2) = 0.5d0*((w + 1.d0)*array(1) + (2.d0 - w)*array(2))
            array(1) = 0.5d0*(1.d0 - w)*array(1)
            ! another pass to order 4
            div = 1.d0/3.d0
            array(4) = div*w*array(3)
            array(3) = div*((w + 1.d0)*array(2) + (3.d0 - w)*array(3))
            array(2) = div*((w + 2.d0)*array(1) + (2.d0 - w)*array(2))
            array(1) = div*(1.d0 - w)*array(1)
            ! compute standard b-spline recursion
            do k = 5, order - 1
                call one_pass_bspline(array, w, k)
            end do
            call diff_bspline(array, darray, order)
            ! one more recursion
            call one_pass_bspline(array, w, order)
            return
        end if  ! ( order == 3 )
    end subroutine fill_bspline_1
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_bspline_2 here]
    subroutine fill_bspline_2(w, order, array, darray, d2array)
        !---------- use standard B-spline recursions: see doc file
        !----- w is fraction between 0 and 1;
        !----- order is the order of interpolation
        !----- array is the array of shifted bsplines
        !----- darray is array of derivs
        !----- d2array is array of 2nd derivs
        ! Using notation from Essmann et al; w = u-[u] and
        ! array(j) = M_n(w + order - j)  where n is order & w = u - [u];
        implicit none
        integer order
        _REAL_ w, array(order), &
            darray(order), d2array(order), div

        integer k
        if (order == 4) then
            ! init order 2
            array(2) = w
            array(1) = 1.d0 - w
            ! in short
            d2array(1) = array(1)
            d2array(2) = array(2) - 2.d0*array(1)
            d2array(3) = array(1) - 2.d0*array(2)
            d2array(4) = array(2)
            ! one pass to order 3
            array(3) = 0.5d0*w*array(2)
            array(2) = 0.5d0*((w + 1.d0)*array(1) + (2.d0 - w)*array(2))
            array(1) = 0.5d0*(1.d0 - w)*array(1)
            ! order four deriv
            darray(1) = -array(1)
            darray(2) = array(1) - array(2)
            darray(3) = array(2) - array(3)
            darray(4) = array(3)
            ! one pass to order 4
            div = 1.d0/3.d0
            array(4) = div*w*array(3)
            array(3) = div*((w + 1.d0)*array(2) + (3.d0 - w)*array(3))
            array(2) = div*((w + 2.d0)*array(1) + (2.d0 - w)*array(2))
            array(1) = div*(1.d0 - w)*array(1)
            return
        else
            ! init order 2
            array(2) = w
            array(1) = 1.d0 - w
            ! one pass to order 3
            array(3) = 0.5d0*w*array(2)
            array(2) = 0.5d0*((w + 1.d0)*array(1) + (2.d0 - w)*array(2))
            array(1) = 0.5d0*(1.d0 - w)*array(1)
            ! compute standard b-spline recursion
            do k = 4, order - 2
                call one_pass_bspline(array, w, k)
            end do
            ! perform standard b-spline differentiation
            call diff_bspline(array, darray, order - 1)
            !   deriv of deriv
            call diff_bspline(darray, d2array, order)
            ! perform standard b-spline differentiation
            call one_pass_bspline(array, w, order - 1)
            call diff_bspline(array, darray, order)
            ! one more recursion
            call one_pass_bspline(array, w, order)
            return
        end if  ! ( order == 4 )
    end subroutine fill_bspline_2
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine diff_bspline here]
    subroutine diff_bspline(c, d, n)
        implicit none
        _REAL_ c(*), d(*)
        integer n
        ! Using notation from Essmann et al; w = u-[u] and
        ! array(j) = M_n(w + order - j)  where n is order
        ! DERIVATIVE:    d/dw M_n(w) = M_n-1(w) - M_n-1(w-1)
        ! i.e.   d/dw M_n(w+n-j) = M_n-1(w+n-j) - M_n-1(w+n-j-1)
        ! i.e.   new(j) = old(j-1) - old(j)
        ! where old is array before one_pass (thus n->n-1) and new is array afterwards

        integer j
        d(1) = -c(1)
        do j = 2, n - 1
            d(j) = c(j - 1) - c(j)
        end do
        d(n) = c(n - 1)
        return
    end subroutine diff_bspline
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine one_pass_bspline here]
    subroutine one_pass_bspline(c, w, n)
        implicit none
        _REAL_ c(*), w
        integer n
        ! Using notation from Essmann et al; w = u-[u] and
        ! array(j) = M_n(w + order - j)  where n is order
        ! RECURSION:  M_n(w) = (w/(n-1))*M_n-1(w)+((n-w)/(n-1))*M_n-1(w-1)
        ! i.e.   M_n(w+n-j) = ((w+n-j)/(n-1))*M_n-1(w+n-j)+((j-w)/(n-1))*M_n-1(w+n-j-1)
        ! i.e.   new(j) = ((w+n-j)/(n-1))*old(j-1) + ((j-w)/(n-1))*old(j)
        ! where old is array before one_pass (thus n->n-1) and new is array afterwards
        ! write backwards to do it with one array

        _REAL_ div
        integer j

        div = 1.d0/(n - 1)
        c(n) = div*w*c(n - 1)
        do j = 1, n - 2
            c(n - j) = div*((w + j)*c(n - j - 1) + (n - j - w)*c(n - j))
        end do
        c(1) = div*(1 - w)*c(1)
        return
    end subroutine one_pass_bspline
!-------------------------------------------------------
!-------------------------------------------------------
!    LOAD PREFACS
!-------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine load_prefacs here]
    subroutine load_prefacs(prefac1, prefac2, prefac3, &
        nfft1, nfft2, nfft3, order, opt_infl)
        implicit none
        _REAL_ prefac1(*), prefac2(*), prefac3(*)
        integer nfft1, nfft2, nfft3, order, opt_infl

        _REAL_ array(bspl_maxorder), w
        _REAL_ bsp_arr(maxnfft), bsp_mod(maxnfft)
        integer i, maxn, num

        ! this routine loads the moduli of the inverse DFT of the B splines
        ! bsp_mod1-3 hold these values, nfft1-3 are the grid dimensions,
        ! Order is the order of the B spline approx.

        if (order > bspl_maxorder) then
            write (6, *) &
                'order too large! check on BSPL_MAXORDER in ew_bspline.f'
            call mexit(6, 1)
            stop
        end if
        maxn = max(nfft1, nfft2, nfft3)
        if (maxn > maxnfft) then
            write (6, *) 'nfft1-3 too large! check on MAXNFFT in ew_bspline.f'
            call mexit(6, 1)
            stop
        end if
        w = 0.d0
        num = 1
        call fill_bspline_0(w, order, array)
        do i = 1, maxn
            bsp_arr(i) = 0.d0
        end do
        do i = 1, order - 1
            bsp_arr(i) = array(order - i)
        end do
        call dftmod(bsp_mod, bsp_arr, nfft1)
        call factor_lambda(bsp_mod, nfft1, order, prefac1, opt_infl)
        call dftmod(bsp_mod, bsp_arr, nfft2)
        call factor_lambda(bsp_mod, nfft2, order, prefac2, opt_infl)
        call dftmod(bsp_mod, bsp_arr, nfft3)
        call factor_lambda(bsp_mod, nfft3, order, prefac3, opt_infl)
        return
    end subroutine load_prefacs
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dftmod here]
    subroutine dftmod(bsp_mod, bsp_arr, nfft)
        use constants, only : TWOPI
        implicit none
        integer nfft
        _REAL_ bsp_mod(nfft), bsp_arr(nfft)
        ! Computes the modulus of the discrete fourier transform of bsp_arr,
        !  storing it into bsp_mod

        integer j, k
        _REAL_ sum1, sum2, arg, tiny
        tiny = 1.d-7
        do k = 1, nfft
            sum1 = 0.d0
            sum2 = 0.d0
            do j = 1, nfft
                arg = TWOPI*(k - 1)*(j - 1)/nfft
                sum1 = sum1 + bsp_arr(j)*cos(arg)
                sum2 = sum2 + bsp_arr(j)*sin(arg)
            end do
            bsp_mod(k) = sum1**2 + sum2**2
        end do
        ! Fix the ONE case where exponential Euler spline interpolation fails
        ! this arbitrary assignment to avoid division by zero is okay
        ! since it happens with highest frequency reciprocal vectors out in tail
        ! of reciprocal sum
        do k = 1, nfft
            if (bsp_mod(k) < tiny) &
                bsp_mod(k) = 0.5d0*(bsp_mod(k - 1) + bsp_mod(k + 1))
        end do
        return
    end subroutine dftmod
!---------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine factor_lambda here]
    subroutine factor_lambda(bsp_mod, nfft, order, prefac, opt_infl)
        use constants, only : PI
        implicit none
        integer nfft, order, opt_infl
        _REAL_ bsp_mod(nfft), prefac(nfft)
        integer kcut
        parameter(kcut=50)

        ! factor in optimal lambda coefficient for bspline approximant
        ! of complex exponentials, thus modify influence function

        integer k, nf, m, order2
        _REAL_ lambda, x, gsum, gsum2
        !     _REAL_ tp(MAXNFFT),t2p(MAXNFFT)
        nf = nfft/2
        order2 = 2*order

        !     ---something wrong with get_tarray for large order.
        !     use clunky but reliable gamma_sum
        !     call get_tarray(tp,order,nfft)
        !     call get_tarray(t2p,order2,nfft)

        do k = 1, nfft
            if (opt_infl == 0) then
                lambda = 1.d0
            else
                m = k - 1
                if (k > nf) m = k - 1 - nfft
                x = (pi*m)/nfft
                order2 = 2*order
                if (m == 0) then
                    lambda = 1.d0
                else
                    call gamma_sum(gsum, m, nfft, order, kcut)
                    call gamma_sum(gsum2, m, nfft, order2, kcut)
                    lambda = gsum/gsum2
                end if
                !        lambda = tp(k)/t2p(k)
            end if
            prefac(k) = lambda**2/bsp_mod(k)
        end do
        return
    end subroutine factor_lambda
!------------------------------------------------------------------------
!------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_tarray here]
    subroutine get_tarray(tarray, order, nfft)
        use constants, only : PI
        implicit none
        _REAL_ tarray(*)
        integer order, nfft

        integer j, k, nf, m
        _REAL_ x, g0, gp, sum, prefac
        _REAL_ alpha(bspl_maxorder + 1)
        _REAL_ powcos(bspl_maxorder + 1)
        _REAL_ powsin(bspl_maxorder + 1), term

        nf = nfft/2

        call recur(alpha, order - 1)
        tarray(1) = 1.d0
        !  the above sum is obtained by differentiating cot(x) order-1 times
        do k = 2, nfft
            m = k - 1
            if (k > nf) m = k - 1 - nfft
            x = pi*dble(m)/nfft
            g0 = cos(x)/sin(x)
            gp = 1.d0
            powcos(1) = 1.d0
            powsin(1) = 1.d0
            do j = 1, order
                powcos(j + 1) = cos(x)*powcos(j)
                powsin(j + 1) = sin(x)*powsin(j)
            end do
            sum = 0.d0
            do j = 1, order + 1
                sum = sum + alpha(j)*powcos(j)*powsin(order + 2 - j)
            end do
            term = x/sin(x)
            prefac = term
            do j = 1, order - 1
                !         prefac = prefac * x / j
                prefac = prefac*(-term)
            end do
            tarray(k) = prefac*sum
        end do
        return
    end subroutine get_tarray
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine recur here]
    subroutine recur(alpha, n)
        !-----------------------------------------------------------
        !     obtains the coefficients alpha(k;n) for the expansion of g_n
        !     where g_0(u) is cos(u)/sin(u), and g_{n}(u) = (d^n/du^n g_0(u))/n!
        !     and the expansion of g_n in terms of g_0 is:

        !     g_n = (-1)^n \times sum_{k=0}^{n+1} alpha(k;n) g_0^k

        !     alpha satisfies the recursion
        !     alpha(0;n+1) = alpha(1;n)/(n+1)
        !     alpha(k;n+1) = ((k+1) alpha(k+1;n) + (k-1) alpha(k-1;n))/(n+1),
        !     for     1 \le k \le n
        !     alpha(n+1;n+1) = (n alpha(n;n))/(n+1) = 0
        !     alpha(n+2;n+1) = ((n+1) alpha(n+1;n))/(n+1) = alpha(n+1;n)
        !     for fortran purposes  alpha(k) = alpha(k-1;n)
        !-------------------------------------------------------------
        implicit none
        integer n
        _REAL_ alpha(*)

        integer m
        !     start things off at n = 1
        alpha(1) = 1.d0
        alpha(2) = 0.d0
        alpha(3) = 1.d0

        do m = 1, n - 1
            call recurstep(alpha, m)
        end do
        return
    end subroutine recur
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine recurstep here]
    subroutine recurstep(alpha, n)
        !-------------------------------------------------------------
        !     updates the coefficients alpha(k;n) for the expansion of g_n
        !     This routine performs one step of the recursion;
        !     input  alpha(k;n)    k = 0,1,...,n+1
        !     output alpha(k;n+1)  k = 0,1,...,n+2
        !     for fortran purposes  alpha(k) = alpha(k-1;n)
        !-------------------------------------------------------------
        implicit none
        integer n
        _REAL_ alpha(*)

        integer k
        _REAL_ beta(bspl_maxorder + 1)

        !     first do a straight copy; beta then holds old values
        do k = 1, n + 2
            beta(k) = alpha(k)
        end do
        !     now recur
        alpha(1) = beta(2)/(n + 1)
        do k = 1, n
            alpha(k + 1) = ((k + 1)*beta(k + 2) + (k - 1)*beta(k))/(n + 1)
        end do
        alpha(n + 2) = (n*beta(n + 1))/(n + 1)
        alpha(n + 3) = beta(n + 2)
        return
    end subroutine recurstep
!---------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine test_diff here]
    subroutine test_diff(w, order, delta)
        implicit none
        _REAL_ w, delta
        integer order, i
        _REAL_ array(25), arrayp(25), &
            arraym(25), darray(25), d2array(25), sav, &
            ndarray(25), nd2array(25), tmp(25), tmp2(25), &
            err

        sav = w
66      format(1x, 2f16.12, 1x, e12.4)
        call fill_bspline_2(w, order, array, darray, d2array)
        w = sav + delta
        call fill_bspline_2(w, order, arrayp, tmp, tmp2)
        w = sav - delta
        call fill_bspline_2(w, order, arraym, tmp, tmp2)
        do i = 1, order
            ndarray(i) = (arrayp(i) - arraym(i))/(2.d0*delta)
            nd2array(i) = &
                (arrayp(i) + arraym(i) - 2.d0*array(i))/delta**2
        end do
        write (6, *) 'analytical, numerical derivs: '
        do i = 1, order
            err = abs(darray(i) - ndarray(i))/abs(darray(i))
            write (6, 66) darray(i), ndarray(i), err
        end do
        write (6, *) 'analytical, numerical 2nd derivs: '
        do i = 1, order
            err = abs(d2array(i) - nd2array(i))/abs(d2array(i))
            write (6, 66) d2array(i), nd2array(i), err
        end do
        return
    end subroutine test_diff
!---------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine gamma_sum here]
    subroutine gamma_sum(gsum, m, nfft, order, kcut)
        use constants, only : PI
        implicit none
        _REAL_ gsum
        integer m, nfft, order, kcut

        _REAL_ frac, x
        integer k

        if (m == 0) then
            gsum = 1.d0
            return
        end if
        frac = dble(m)/nfft
        x = pi*frac
        gsum = 1.d0
        do k = 1, kcut
            gsum = gsum + (x/(x + pi*k))**order
        end do
        do k = 1, kcut
            gsum = gsum + (x/(x - pi*k))**order
        end do
        return
    end subroutine gamma_sum

end module ew_bspline
