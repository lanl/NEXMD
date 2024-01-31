! <compile=optimized>

module ew_recip

#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"

    _REAL_, allocatable :: m1_tbl(:)
    _REAL_, allocatable :: m2_tbl(:)
    _REAL_, allocatable :: m3_tbl(:)

    logical, save :: first_pme

    public deallocate_m1m2m3, do_pmesh_kspace

    _REAL_, dimension(3), save :: frcx

!private fill_charge_grid, grad_sumrc, scalar_sumrc, scalar_sumrc_orthog

!===========================================================================
contains
!===========================================================================

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Deallocate the module tables
    subroutine deallocate_m1m2m3()
        implicit none
        integer:: ier
        if (.not. allocated(m1_tbl)) return
        deallocate (m1_tbl, m2_tbl, m3_tbl, stat=ier)
        REQUIRE(ier == 0)
        return
    end subroutine deallocate_m1m2m3

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine do_pmesh_kspace here]
    subroutine do_pmesh_kspace(natom, crd, charge, &
        frc, prefac1, prefac2, prefac3, fftable, qm_pot_only)
        use nblist, only : recip, volume
        use stack
        use ew_bspline, only : get_grid_weights

        implicit none
        character(kind=1, len=15) :: routine = "do_pmesh_kspace"
#  include "flocntrl.h"
#  include "ew_pme_recip.h"
#  include "ew_frc.h"
#  include "def_time.h"
#  include "box.h"

#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
        include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif

        ! INPUT
        !       natom:  number of atoms
        !       crd   atomic coords
        !       charge  atomic charges

#ifdef MPI
        integer ierr
#endif

        integer natom, num_ks_trial
        _REAL_ crd(3, natom), charge(natom)
        integer nmine, nderiv

        logical, intent(in) :: qm_pot_only

        ! OUTPUT
        !       eer:  ewald reciprocal or k-space  energy
        !       frc forces incremented by k-space sum

        _REAL_ frc(3, natom)

        ! HEAP STORAGE:  These arrays need to be preserved throughout simulation

        _REAL_ prefac1(*), prefac2(*), prefac3(*), fftable(*)

        integer l_d2th1, l_d2th2, l_d2th3

        integer nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw
        integer l_fr1, l_fr2, l_fr3
        integer l_th1, l_th2, l_th3, l_dth1, l_dth2, l_dth3
        integer l_fftw, l_q
        integer imy_cg
        integer num_ks
        integer i, ii
        save num_ks
        data num_ks/0/
        integer l_tmpy, l_alpha, l_beta
        logical not_done

        nderiv = 1

        !     FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)

        if (do_rec == 0) return

        !     get some integer array dimensions
        call get_fftdims(nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw)
# ifdef MPI
#   ifdef LES
        num_ks = natom
#   else
        if (num_ks == 0) &
            num_ks = min((((nxyslab(0) + order - 1)*natom*4)/ &
            (3*nfft3)), natom)
#   endif
# else
        num_ks = natom
# endif
        call timer_start(TIME_BSPL)
        num_ks_trial = 0
        not_done = .true.
        do while (not_done)
            call get_stack(l_fftw, sizffwrk, routine)
            call get_stack(l_q, siz_q, routine)
            call get_stack(l_fr1, num_ks, routine)
            call get_stack(l_fr2, num_ks, routine)
            call get_stack(l_fr3, num_ks, routine)
            call get_stack(l_th1, num_ks*order, routine)
            call get_stack(l_th2, num_ks*order, routine)
            call get_stack(l_th3, num_ks*order, routine)
            if (nderiv == 1) then
                call get_stack(l_dth1, num_ks*order, routine)
                call get_stack(l_dth2, num_ks*order, routine)
                call get_stack(l_dth3, num_ks*order, routine)
                call get_stack(l_d2th1, order, routine)
                call get_stack(l_d2th2, order, routine)
                call get_stack(l_d2th3, order, routine)
            end if
            if (.not. rstack_ok) then
                deallocate (r_stack)
                allocate (r_stack(1:lastrst), stat=alloc_ier)
                call reassign_rstack(routine)
            end if
            REQUIRE(rstack_ok)

            call get_istack(imy_cg, num_ks, routine)
            if (.not. istack_ok) then
                deallocate (i_stack)
                allocate (i_stack(1:lastist), stat=alloc_ier)
                call reassign_istack(routine)
            end if
            REQUIRE(istack_ok)

            call get_grid_weights( &
                natom, crd, recip, nfft1, nfft2, nfft3, &
                r_stack(l_fr1), r_stack(l_fr2), r_stack(l_fr3), order, &
                r_stack(l_th1), r_stack(l_th2), r_stack(l_th3), &
                r_stack(l_dth1), r_stack(l_dth2), r_stack(l_dth3), &
                r_stack(l_d2th1), r_stack(l_d2th2), r_stack(l_d2th3), &
                i_stack(imy_cg), nmine, nderiv, num_ks)
            if (nmine > num_ks) then
#ifdef MPI
                write (6, '("*****  Processor ",i6)') mytaskid
#endif
                write (6, '("***** System must be very inhomogeneous.")')
                write (6, '("*****  Readjusting recip sizes.")')
                write (6, '(A,i9,A,i9/)') &
                    " In this slab, Atoms found: ", nmine, &
                    "  Allocated: ", num_ks
                if (num_ks_trial >= 2) call mexit(6, 1)
                if (nderiv == 1) then
                    call free_stack(l_d2th3, routine)
                    call free_stack(l_d2th2, routine)
                    call free_stack(l_d2th1, routine)
                    call free_stack(l_dth3, routine)
                    call free_stack(l_dth2, routine)
                    call free_stack(l_dth1, routine)
                end if
                call free_stack(l_th3, routine)
                call free_stack(l_th2, routine)
                call free_stack(l_th1, routine)
                call free_stack(l_fr3, routine)
                call free_stack(l_fr2, routine)
                call free_stack(l_fr1, routine)
                call free_stack(l_q, routine)
                call free_stack(l_fftw, routine)
                call free_istack(imy_cg, routine)
                num_ks = nmine*4/3
                num_ks_trial = num_ks_trial + 1
            else
                not_done = .false.
            end if
        end do
        call timer_stop_start(TIME_BSPL, TIME_FILLG)

        !........Fill Charge Grid
        !          charges are approximated on an even grid
        call fill_charge_grid(natom, charge, &
            r_stack(l_th1), r_stack(l_th2), r_stack(l_th3), &
            r_stack(l_fr1), r_stack(l_fr2), r_stack(l_fr3), &
            order, nfft1, nfft2, nfft3, nfftdim1, nfftdim2 &
# ifdef MPI
            , mxyslabs &
# else
            , nfftdim3 &
# endif
            , r_stack(l_q), i_stack(imy_cg), nmine)

        call timer_stop_start(TIME_FILLG, TIME_FFT)
        call get_stack(l_tmpy, 2*nfftdim1, routine)
        call get_stack(l_alpha, nfft1, routine)
        call get_stack(l_beta, nfft1, routine)
        if (.not. rstack_ok) then
            deallocate (r_stack)
            allocate (r_stack(1:lastrst), stat=alloc_ier)
            call reassign_rstack(routine)
        end if
        REQUIRE(rstack_ok)

        call fft_backrc( &
            r_stack(l_q), fftable, r_stack(l_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork &
            , r_stack(l_tmpy), &
            r_stack(l_alpha), r_stack(l_beta) &
            )

        call timer_stop_start(TIME_FFT, TIME_SCSUM)

        !           -------------SCALAR SUM------------------

        if (ifbox == 1) then
            call scalar_sumrc_orthog( &
                r_stack(l_q), &
                ew_coeff, volume, recip, prefac1, prefac2, prefac3, &
                nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, eer, rec_vir)
        else
            call scalar_sumrc( &
                r_stack(l_q), &
                ew_coeff, volume, recip, prefac1, prefac2, prefac3, &
                nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, eer, rec_vir)
        end if
!call dumpq(r_stack(l_q),nfft3,nfftdim1,nfft3,24)

        call timer_stop_start(TIME_SCSUM, TIME_FFT)

        !           -----------FFT FORWARD--------------------

        call fft_forwardrc( &
            r_stack(l_q), fftable, r_stack(l_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork &
            , r_stack(l_tmpy), &
            r_stack(l_alpha), r_stack(l_beta) &
            )
        call free_stack(l_beta, routine)
        call free_stack(l_alpha, routine)
        call free_stack(l_tmpy, routine)

        call timer_stop_start(TIME_FFT, TIME_GRADS)

        !           -----------GRAD SUM--------------------
        call grad_sumrc( &
            natom, charge, recip, &
            r_stack(l_th1), r_stack(l_th2), r_stack(l_th3), &
            r_stack(l_dth1), r_stack(l_dth2), r_stack(l_dth3), &
            frc, r_stack(l_fr1), r_stack(l_fr2), r_stack(l_fr3), &
            order, nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, &
            r_stack(l_q), &
            i_stack(imy_cg), nmine, &
            qm_pot_only)

        call timer_stop(TIME_GRADS)
        if (nderiv == 1) then
            call free_stack(l_d2th3, routine)
            call free_stack(l_d2th2, routine)
            call free_stack(l_d2th1, routine)
            call free_stack(l_dth3, routine)
            call free_stack(l_dth2, routine)
            call free_stack(l_dth1, routine)
        end if
        call free_stack(l_th3, routine)
        call free_stack(l_th2, routine)
        call free_stack(l_th1, routine)
        call free_stack(l_fr3, routine)
        call free_stack(l_fr2, routine)
        call free_stack(l_fr1, routine)
        call free_stack(l_q, routine)
        call free_stack(l_fftw, routine)

        call free_istack(imy_cg, routine)
#ifdef MPI
        call mpi_barrier(recip_comm, ierr)
#endif
        return
    end subroutine do_pmesh_kspace
!-------------------------------------------------------------------
!                                -----
!     --- FILL_CHARGE_GRID -- RC-------
!                                -----
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_charge_grid here]
    subroutine fill_charge_grid(natom, charge, &
        theta1, theta2, theta3, fr1, fr2, fr3, &
        order, nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, q, my_cg, nmine)

        !---------------------------------------------------------------------
        ! INPUT:
        !      natom:  number of atoms
        !      charge: the array of atomic charges
        !      theta1,theta2,theta3: the spline coeff arrays
        !      fr1,fr2,fr3 the scaled and shifted fractional coords
        !      nfft1,nfft2,nfft3: the charge grid dimensions
        !      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims
        !      order: the order of spline interpolation
        ! OUTPUT:
        !      Q the charge grid
        !---------------------------------------------------------------------

        use ew_bspline, only : kbot, ktop

        implicit none
        integer natom, order, nfft1, nfft2, nfft3
        integer nfftdim1, nfftdim2, nfftdim3
        _REAL_ fr1(natom), fr2(natom), fr3(natom)
        _REAL_ theta1(order, natom), theta2(order, natom), &
            theta3(order, natom), charge(natom)
        _REAL_ q(*)
        integer my_cg(*), nmine

#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
        integer kbot0
#endif

        integer n, ith1, ith2, ith3, i0, j0, k0, i, j, k, i00, j00, iqk, iqj
        _REAL_ prod
        integer im, kq, ntot

        !........Zero the Charge grids

        ntot = 2*nfftdim1*nfftdim2*nfftdim3
        call zero_array(q, ntot)
#ifdef MPI
        kbot0 = mxystart(mytaskid)
        kbot = kbot0 + 1
        ktop = kbot0 + mxyslabs
#endif

        do im = 1, nmine
#ifdef MPI
            n = my_cg(im)
#else
            n = im
#endif
            k0 = int(fr3(im)) - order
            j00 = int(fr2(im)) - order
            i00 = int(fr1(im)) - order
            do ith3 = 1, order
                k0 = k0 + 1
                if (k0 >= 0) then
                    k = k0 + 1
                else
                    k = k0 + 1 + nfft3
                end if

#ifdef MPI
                if (k >= kbot .and. k <= ktop) then
                    kq = k - kbot0
#else
                kq = k
#endif
                iqk = (kq - 1)*2*nfftdim1*nfftdim2
                j0 = j00
                do ith2 = 1, order
                    j0 = j0 + 1

                    iqj = iqk + j0*2*nfftdim1
                    if (j0 < 0) iqj = iqj + 2*nfft2*nfftdim1

                    prod = theta2(ith2, im)*theta3(ith3, im)*charge(n)
                    i0 = i00 + 1
                    do ith1 = 1, order
                        i0 = i0 + 1
                        if (i0 >= 1) then
                            q(i0 + iqj) = q(i0 + iqj) + theta1(ith1, im)*prod
                        else
                            q(i0 + nfft1 + iqj) = q(i0 + nfft1 + iqj) + theta1(ith1, im)*prod
                        end if
                    end do
                end do
#ifdef MPI
            end if
#endif
        end do  !  ith3 = 1,order
    end do  !  im = 1,nmine
    return
end subroutine fill_charge_grid

!-------------------------------------------------------------------
!  GRAD_SUM RC    **** REAL not complex ****
!-----------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine grad_sumrc here]
subroutine grad_sumrc( &
    natom, charge, recip, theta1, theta2, theta3, &
    dtheta1, dtheta2, dtheta3, frc, fr1, fr2, fr3, &
    order, nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, &
    q, my_cg, nmine, &
    qm_pot_only)
    use constants, only : INV_AMBER_ELECTROSTATIC, zero
    use qmmm_module, only : qmewald, qmmm_struct, qmmm_nml
    use ew_bspline, only : kbot, ktop
    use decomp, only : decpr, decpair
    use crg_reloc, only : ifcr, cr_dcdr_tbl, cr_add_dcdr_factor
    use file_io_dat

    implicit none
    integer natom, order, nfft1, nfft2, nfft3
    integer nfftdim1, nfftdim2, nfftdim3
    integer kq
    logical, intent(in) :: qm_pot_only
#include "box.h"
#include "md.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#  ifdef MPI_DOUBLE_PRECISION
#    undef MPI_DOUBLE_PRECISION
#  endif
    include 'mpif.h'
#  ifdef CRAY_PVP
#    define MPI_DOUBLE_PRECISION MPI_REAL8
#  endif
    _REAL_ q(*)
#else
    _REAL_ q(*)

    ! If MPI only some atoms are touched by the local process. Hence loop is over
    ! im=1,nmine, arrays are filled in order of being seen by process
    ! and n is obtained by a pointer.

#endif

    _REAL_ recip(3, 3), recip11, recip22, recip33
    _REAL_ fr1(natom), fr2(natom), fr3(natom)
    _REAL_ frc(3, natom)
    _REAL_ theta1(order, natom), theta2(order, natom), &
        theta3(order, natom), charge(natom)
    _REAL_ dtheta1(order, natom), dtheta2(order, natom), &
        dtheta3(order, natom)
    integer my_cg(*), nmine
    integer iqk, iqj, j00, i00

    integer n, nn, ith1, ith2, ith3, i0, j0, k0, i, j, k, im
    _REAL_ f1, f2, f3, dedc, term, chargen
    _REAL_ dfx, dfy, dfz, dnfft1, dnfft2, dnfft3
    _REAL_ f1fac, f2fac, f3fac
    _REAL_ dectmp
    logical do_atom, fluc_crg

    dnfft1 = nfft1
    dnfft2 = nfft2
    dnfft3 = nfft3

    recip11 = recip(1, 1)*dnfft1
    recip22 = recip(2, 2)*dnfft2
    recip33 = recip(3, 3)*dnfft3

    if (qm_pot_only) then
!  ===============  Below for qm grad sum yielding recip potential ======
        !1 to nquant_nlink
        qmewald%mmpot(1:qmmm_struct%nquant_nlink) = zero
        do im = 1, nmine
#ifdef MPI
            n = my_cg(im)
#else
            n = im
#endif
            do_atom = (qmmm_struct%atom_mask(n) .or. qmmm_struct%mm_link_mask(n))
            if (do_atom) then
                f1 = zero
                i00 = int(fr1(im)) - order
                j00 = int(fr2(im)) - order
                k0 = int(fr3(im)) - order
                do ith3 = 1, order
                    k0 = k0 + 1
                    if (k0 >= 0) then
                        k = k0 + 1
                    else
                        k = k0 + 1 + nfft3
                    end if

#ifdef MPI
                    if (k >= kbot .and. k <= ktop) then
                        kq = k - mxystart(mytaskid)
#else
                    kq = k
#endif
                    iqk = (kq - 1)*2*nfftdim1*nfftdim2
                    j0 = j00
                    do ith2 = 1, order
                        j0 = j0 + 1

                        iqj = iqk + j0*2*nfftdim1
                        if (j0 < 0) iqj = iqj + 2*nfft2*nfftdim1

                        i0 = i00 + 1
                        f1fac = theta2(ith2, im)*theta3(ith3, im)
                        do ith1 = 1, order
                            i0 = i0 + 1
                            if (i0 >= 1) then
                                term = q(i0 + iqj)
                            else
                                term = q(i0 + nfft1 + iqj)
                            end if
                            f1 = f1 + term*theta1(ith1, im)*f1fac
                        end do
                    end do
#ifdef MPI
                end if
#endif
            end do  !  ith3 = 1,order

            do nn = 1, qmmm_struct%nquant_nlink
                if (qmmm_struct%iqmatoms(nn) == n) &
                    qmewald%mmpot(nn) = f1*INV_AMBER_ELECTROSTATIC
            end do

        end if
    end do
!  ===============  Below for nornmal grad sum yielding forces ======

else
    do im = 1, nmine
#ifdef MPI
        n = my_cg(im)
#else
        n = im
#endif
        ! Does this charge fluctuate?
        if (ifcr /= 0) then
            if (cr_dcdr_tbl(n) /= 0) then
                fluc_crg = .true.
                dedc = zero
            else
                fluc_crg = .false.
            end if
        else
            fluc_crg = .false.
        end if

        f1 = zero
        f2 = zero
        f3 = zero
        i00 = int(fr1(im)) - order
        j00 = int(fr2(im)) - order
        k0 = int(fr3(im)) - order
        chargen = charge(n)
        do ith3 = 1, order
            k0 = k0 + 1
            if (k0 >= 0) then
                k = k0 + 1
            else
                k = k0 + 1 + nfft3
            end if

#ifdef MPI
            if (k >= kbot .and. k <= ktop) then
                kq = k - mxystart(mytaskid)
#else
            kq = k
#endif
            iqk = (kq - 1)*2*nfftdim1*nfftdim2
            j0 = j00
            do ith2 = 1, order
                j0 = j0 + 1

                iqj = iqk + j0*2*nfftdim1
                if (j0 < 0) iqj = iqj + 2*nfft2*nfftdim1

                i0 = i00 + 1

                f1fac = theta2(ith2, im)*theta3(ith3, im)*chargen
                f2fac = dtheta2(ith2, im)*theta3(ith3, im)*chargen
                f3fac = theta2(ith2, im)*dtheta3(ith3, im)*chargen
                do ith1 = 1, order
                    i0 = i0 + 1
                    if (i0 >= 1) then
                        term = q(i0 + iqj)
                    else
                        term = q(i0 + nfft1 + iqj)
                    end if

                    !               ---force is negative of grad
                    f1 = f1 - term*dtheta1(ith1, im)*f1fac
                    f2 = f2 - term*theta1(ith1, im)*f2fac
                    f3 = f3 - term*theta1(ith1, im)*f3fac
                    if (fluc_crg) then
                        dedc = dedc + term*theta1(ith1, im)*theta2(ith2, im) &
                            *theta3(ith3, im)
                    end if

                    ! -- ti decomp
                    if (decpr .and. idecomp > 0) then
                        dectmp = term*theta1(ith1, im)*theta2(ith2, im)*theta3(ith3, im)*chargen*0.5d0
                        call decpair(2, n, n, dectmp/(nstlim/ntpr))
                    end if
                end do
            end do
#ifdef MPI
        end if
#endif
    end do  !  ith3 = 1,order

    if (ifbox == 1) then  ! orthogonal unit cell
        dfx = recip11*f1
        dfy = recip22*f2
        dfz = recip33*f3
    else
        f1 = dnfft1*f1
        f2 = dnfft2*f2
        f3 = dnfft3*f3
        dfx = recip(1, 1)*f1 + recip(1, 2)*f2 + recip(1, 3)*f3
        dfy = recip(2, 1)*f1 + recip(2, 2)*f2 + recip(2, 3)*f3
        dfz = recip(3, 1)*f1 + recip(3, 2)*f2 + recip(3, 3)*f3
    end if
    frc(1, n) = frc(1, n) + dfx
    frc(2, n) = frc(2, n) + dfy
    frc(3, n) = frc(3, n) + dfz
    frcx(1) = frcx(1) + dfx
    frcx(2) = frcx(2) + dfy
    frcx(3) = frcx(3) + dfz

    if (fluc_crg) then
        call cr_add_dcdr_factor(n, dedc)
    end if
end do  !  im = 1,nmine
end if   !  qm_pot_only
return
end subroutine grad_sumrc
!=================================================================
!     --- SCALAR_SUM RC---
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine scalar_sumrc here]
subroutine scalar_sumrc( &
    q, &
    ewaldcof, volume, recip, prefac1, prefac2, prefac3, &
    nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, eer, rec_vir)
    use constants, only : pi, PI2
    implicit none
#  include "extra.h"
#  include "box.h"
#  include "md.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
    integer ier, i
#endif
    integer nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3
    _REAL_ prefac1(nfft1), prefac2(nfft2), &
        prefac3(nfft3), ewaldcof, volume
    _REAL_ eer, rec_vir(3, 3)
    _REAL_ recip(3, 3)
    integer k2q
    _REAL_ q(2, nfft3, nfftdim1, nfft2)
    _REAL_ fac, denom, eterm, vterm, energy
    integer k1, k2, k3, m1, m2, m3, nff, indtop
    integer nf1, nf2, nf3
    integer k10
    _REAL_ mhat1, mhat2, mhat3, msq, struc2, msq_inv, piv_inv
    integer k1s, k2s, k3s, m1s, m2s, m3s
    _REAL_ mhat1s, mhat2s, mhat3s, msqs
    _REAL_ struc2s, eterms, vterms, denoms, tmp1, tmp2

    indtop = nfft1*nfft2*nfft3
    piv_inv = 1.d0/(pi*volume)
    fac = PI2/(ewaldcof*ewaldcof)
    nff = nfft1*nfft2
    nf1 = nfft1/2
    if (2*nf1 < nfft1) nf1 = nf1 + 1
    nf2 = nfft2/2
    if (2*nf2 < nfft2) nf2 = nf2 + 1
    nf3 = nfft3/2
    if (2*nf3 < nfft3) nf3 = nf3 + 1
    energy = 0.d0
    k10 = 1
#ifndef noVIRIAL
    do m2 = 1, 3
        do m1 = 1, 3
            rec_vir(m1, m2) = 0.d0
        end do
    end do
#endif
    !........Insist that Q(1,1,1,1) is set to 0 (true already for neutral)

    if (master) then
        q(1, 1, 1, 1) = 0.d0
        q(2, 1, 1, 1) = 0.d0
    end if

    !======================================================================
    !        BIG LOOP
    !======================================================================

#ifdef MPI
    do k2q = 1, mxzslabs
        if (master) then
            k2 = k2q
            k2s = mod(nfft2 - k2 + 1, nfft2) + 1
        else
            k2 = k2q + mxzstart(mytaskid)
            k2s = mod(nfft2 - k2 + 1, nfft2) + 1
        end if
#else
    do k2q = 1, nfft2
        k2s = mod(nfft2 - k2q + 1, nfft2) + 1
        k2 = k2q
#endif
        m2 = k2 - 1
        if (k2 > nf2) m2 = k2 - 1 - nfft2

        do k3 = 1, nfft3
            k3s = mod(nfft3 - k3 + 1, nfft3) + 1
            m3 = k3 - 1
            if (k3 > nf3) m3 = k3 - 1 - nfft3
            k10 = 1
            if (master) then
                if (k3 + k2 == 2) k10 = 2
            end if
            do k1 = k10, nf1 + 1
                k1s = nfft1 - k1 + 2
                m1 = k1 - 1
                if (k1 > nf1) m1 = k1 - 1 - nfft1
                mhat1 = recip(1, 1)*m1 + recip(1, 2)*m2 + recip(1, 3)*m3
                mhat2 = recip(2, 1)*m1 + recip(2, 2)*m2 + recip(2, 3)*m3
                mhat3 = recip(3, 1)*m1 + recip(3, 2)*m2 + recip(3, 3)*m3
                msq = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3
                msq_inv = 1.d0/msq
                eterm = exp(-fac*msq)*prefac1(k1)*prefac2(k2)*prefac3(k3) &
                    *piv_inv*msq_inv
                struc2 = q(1, k3, k1, k2q)*q(1, k3, k1, k2q) + &
                    q(2, k3, k1, k2q)*q(2, k3, k1, k2q)
                tmp1 = eterm*struc2
                energy = energy + tmp1
#ifndef noVIRIAL
                vterm = 2.d0*(fac*msq + 1.d0)*msq_inv
                tmp2 = tmp1*vterm
                rec_vir(1, 1) = rec_vir(1, 1) + tmp1*(vterm*mhat1*mhat1 - 1.d0)
                rec_vir(1, 2) = rec_vir(1, 2) + tmp2*mhat1*mhat2
                rec_vir(1, 3) = rec_vir(1, 3) + tmp2*mhat1*mhat3
                rec_vir(2, 1) = rec_vir(2, 1) + tmp2*mhat2*mhat1
                rec_vir(2, 2) = rec_vir(2, 2) + tmp1*(vterm*mhat2*mhat2 - 1.d0)
                rec_vir(2, 3) = rec_vir(2, 3) + tmp2*mhat2*mhat3
                rec_vir(3, 1) = rec_vir(3, 1) + tmp2*mhat3*mhat1
                rec_vir(3, 2) = rec_vir(3, 2) + tmp2*mhat3*mhat2
                rec_vir(3, 3) = rec_vir(3, 3) + tmp1*(vterm*mhat3*mhat3 - 1.d0)
#endif

                if (k1 > 1 .and. k1 <= nfft1) then
                    m1s = k1s - 1
                    if (k1s > nf1) m1s = k1s - 1 - nfft1
                    m2s = k2s - 1
                    if (k2s > nf2) m2s = k2s - 1 - nfft2
                    m3s = k3s - 1
                    if (k3s > nf3) m3s = k3s - 1 - nfft3
                    mhat1s = recip(1, 1)*m1s + recip(1, 2)*m2s + recip(1, 3)*m3s
                    mhat2s = recip(2, 1)*m1s + recip(2, 2)*m2s + recip(2, 3)*m3s
                    mhat3s = recip(3, 1)*m1s + recip(3, 2)*m2s + recip(3, 3)*m3s
                    msqs = mhat1s*mhat1s + mhat2s*mhat2s + mhat3s*mhat3s
                    msq_inv = 1.d0/msqs
                    eterms = exp(-fac*msqs)*prefac1(k1s)*prefac2(k2s)* &
                        prefac3(k3s)*piv_inv*msq_inv
                    tmp1 = eterms*struc2
                    energy = energy + tmp1
#ifndef noVIRIAL
                    vterms = 2.d0*(fac*msqs + 1.d0)*msq_inv
                    tmp2 = tmp1*vterms
                    rec_vir(1, 1) = rec_vir(1, 1) + tmp1*(vterms*mhat1s*mhat1s - 1.d0)
                    rec_vir(1, 2) = rec_vir(1, 2) + tmp2*mhat1s*mhat2s
                    rec_vir(1, 3) = rec_vir(1, 3) + tmp2*mhat1s*mhat3s
                    rec_vir(2, 1) = rec_vir(2, 1) + tmp2*mhat2s*mhat1s
                    rec_vir(2, 2) = rec_vir(2, 2) + tmp1*(vterms*mhat2s*mhat2s - 1.d0)
                    rec_vir(2, 3) = rec_vir(2, 3) + tmp2*mhat2s*mhat3s
                    rec_vir(3, 1) = rec_vir(3, 1) + tmp2*mhat3s*mhat1s
                    rec_vir(3, 2) = rec_vir(3, 2) + tmp2*mhat3s*mhat2s
                    rec_vir(3, 3) = rec_vir(3, 3) + tmp1*(vterms*mhat3s*mhat3s - 1.d0)
#endif
                end if

                q(1, k3, k1, k2q) = eterm*q(1, k3, k1, k2q)
                q(2, k3, k1, k2q) = eterm*q(2, k3, k1, k2q)
            end do
        end do
    end do

    eer = 0.5d0*energy

#ifndef noVIRIAL
    do m2 = 1, 3
        do m1 = 1, 3
            rec_vir(m1, m2) = 0.5d0*rec_vir(m1, m2)
        end do
    end do
#endif

    return
end subroutine scalar_sumrc

!=================================================================
!     --- SCALAR_SUM RC---  for ORTHOGONAL CELL
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine scalar_sumrc_orthog here]
subroutine scalar_sumrc_orthog( &
    q, &
    ewaldcof, volume, recip, prefac1, prefac2, prefac3, &
    nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, eer, rec_vir)
    use constants, only : PI, PI2, HALF, ZERO
    implicit none
#  include "extra.h"
#  include "box.h"
#  include "md.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#endif
    integer ier, i, itot
    integer nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3
    _REAL_ prefac1(nfft1), prefac2(nfft2), &
        prefac3(nfft3), ewaldcof, volume
    _REAL_ eer, rec_vir(3, 3)
    _REAL_ recip(3, 3)
    integer k2q
    _REAL_ q(2, nfft3, nfftdim1, nfft2)
    _REAL_ fac, denom, eterm, vterm, energy
    integer k1, k2, k3, m1, m2, m3, nff, indtop
    integer nf1, nf2, nf3
    integer k10
    _REAL_ mhat1, mhat2, mhat3, msq, struc2, msq_inv, piv_inv
    integer k1s, k2s, k3s, m1s, m2s, m3s
    _REAL_ mhat1s, mhat2s, mhat3s, msqs
    _REAL_ struc2s, eterms, vterms, denoms
    _REAL_ tmp1, tmp2, m2_m3_tbl, m2_m3_tbls
    _REAL_ recip11, recip22, recip33
    _REAL_ recip11sq, recip22sq, recip33sq

    indtop = nfft1*nfft2*nfft3
    piv_inv = 1.d0/(pi*volume)
    fac = PI2/(ewaldcof*ewaldcof)
    nff = nfft1*nfft2
    nf1 = nfft1/2
    if (2*nf1 < nfft1) nf1 = nf1 + 1
    nf2 = nfft2/2
    if (2*nf2 < nfft2) nf2 = nf2 + 1
    nf3 = nfft3/2
    if (2*nf3 < nfft3) nf3 = nf3 + 1
    recip11 = recip(1, 1)
    recip22 = recip(2, 2)
    recip33 = recip(3, 3)

    if (first_pme) then
        allocate (m1_tbl(-(nfft1/2 + 1):nfft1/2 + 1), &
            m2_tbl(-(nfft2/2 + 1):nfft2/2 + 1), &
            m3_tbl(-(nfft3/2 + 1):nfft3/2 + 1), &
            stat=ier)
        REQUIRE(ier == 0)
    end if

    if (ntp > 0 .or. first_pme) then
        recip11sq = recip11*recip11
        recip22sq = recip22*recip22
        recip33sq = recip33*recip33
        itot = 0
        do i = -(nfft1/2 + 1), nfft1/2 + 1
            itot = itot + 1
            m1_tbl(i) = -fac*dble(i)*dble(i)*recip11sq
        end do
        call vdexp(itot, m1_tbl(-nfft1/2 - 1), m1_tbl(-nfft1/2 - 1))

        itot = 0
        do i = -(nfft2/2 + 1), nfft2/2 + 1
            itot = itot + 1
            m2_tbl(i) = -fac*dble(i)*dble(i)*recip22sq
        end do
        call vdexp(itot, m2_tbl(-nfft2/2 - 1), m2_tbl(-nfft2/2 - 1))

        itot = 0
        do i = -(nfft3/2 + 1), nfft3/2 + 1
            itot = itot + 1
            m3_tbl(i) = -fac*dble(i)*dble(i)*recip33sq
        end do
        call vdexp(itot, m3_tbl(-nfft3/2 - 1), m3_tbl(-nfft3/2 - 1))

        first_pme = .false.
    end if
    energy = 0.d0
    k10 = 1
#ifndef noVIRIAL
    rec_vir = zero
#endif

    !........Insist that Q(1,1,1,1) is set to 0 (true already for neutral)

    if (master) then
        q(1, 1, 1, 1) = zero
        q(2, 1, 1, 1) = zero
    end if

    !======================================================================
    !        BIG LOOP
    !======================================================================

#ifdef MPI
    do k2q = 1, mxzslabs
        if (master) then
            k2 = k2q
            k2s = mod(nfft2 - k2 + 1, nfft2) + 1
        else
            k2 = k2q + mxzstart(mytaskid)
            k2s = mod(nfft2 - k2 + 1, nfft2) + 1
        end if
#else
    do k2q = 1, nfft2
        k2s = mod(nfft2 - k2q + 1, nfft2) + 1
        k2 = k2q
#endif
        m2 = k2 - 1
        if (k2 > nf2) m2 = k2 - 1 - nfft2
        mhat2 = recip22*m2
        m2s = k2s - 1
        if (k2s > nf2) m2s = k2s - 1 - nfft2
        mhat2s = recip22*m2s

        do k3 = 1, nfft3
            k3s = mod(nfft3 - k3 + 1, nfft3) + 1
            m3 = k3 - 1
            if (k3 > nf3) m3 = k3 - 1 - nfft3
            mhat3 = recip33*m3
            m3s = k3s - 1
            if (k3s > nf3) m3s = k3s - 1 - nfft3
            mhat3s = recip33*m3s

            k10 = 1
            if (master) then
                if (k3 + k2 == 2) k10 = 2
            end if

            m2_m3_tbl = piv_inv*m2_tbl(m2)*m3_tbl(m3) &
                *prefac2(k2)*prefac3(k3)
            m2_m3_tbls = piv_inv*m2_tbl(m2s)*m3_tbl(m3s) &
                *prefac2(k2s)*prefac3(k3s)

            do k1 = k10, nf1 + 1
                k1s = nfft1 - k1 + 2
                m1 = k1 - 1
                if (k1 > nf1) m1 = k1 - 1 - nfft1
                mhat1 = recip11*m1
                msq = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3
                msq_inv = 1.d0/msq
                eterm = m1_tbl(m1)*m2_m3_tbl* &
                    prefac1(k1)*msq_inv
                struc2 = q(1, k3, k1, k2q)*q(1, k3, k1, k2q) + &
                    q(2, k3, k1, k2q)*q(2, k3, k1, k2q)
                tmp1 = eterm*struc2
                energy = energy + tmp1
#ifndef noVIRIAL
                vterm = 2.d0*(fac*msq + 1.d0)*msq_inv
                tmp2 = tmp1*vterm
                rec_vir(1, 1) = rec_vir(1, 1) + tmp1*(vterm*mhat1*mhat1 - 1.d0)
                rec_vir(1, 2) = rec_vir(1, 2) + tmp2*mhat1*mhat2
                rec_vir(1, 3) = rec_vir(1, 3) + tmp2*mhat1*mhat3
                rec_vir(2, 1) = rec_vir(2, 1) + tmp2*mhat2*mhat1
                rec_vir(2, 2) = rec_vir(2, 2) + tmp1*(vterm*mhat2*mhat2 - 1.d0)
                rec_vir(2, 3) = rec_vir(2, 3) + tmp2*mhat2*mhat3
                rec_vir(3, 1) = rec_vir(3, 1) + tmp2*mhat3*mhat1
                rec_vir(3, 2) = rec_vir(3, 2) + tmp2*mhat3*mhat2
                rec_vir(3, 3) = rec_vir(3, 3) + tmp1*(vterm*mhat3*mhat3 - 1.d0)
#endif

                if (k1 > 1 .and. k1 <= nfft1) then
                    m1s = k1s - 1
                    if (k1s > nf1) m1s = k1s - 1 - nfft1
                    mhat1s = recip11*m1s
                    msqs = mhat1s*mhat1s + mhat2s*mhat2s + mhat3s*mhat3s
                    msq_inv = 1.d0/msqs
                    eterms = m1_tbl(m1s)*m2_m3_tbls* &
                        prefac1(k1s)*msq_inv
                    tmp1 = eterms*struc2
                    energy = energy + tmp1
#ifndef noVIRIAL
                    vterms = 2.d0*(fac*msqs + 1.d0)*msq_inv
                    tmp2 = tmp1*vterms
                    rec_vir(1, 1) = rec_vir(1, 1) + tmp1*(vterms*mhat1s*mhat1s - 1.d0)
                    rec_vir(1, 2) = rec_vir(1, 2) + tmp2*mhat1s*mhat2s
                    rec_vir(1, 3) = rec_vir(1, 3) + tmp2*mhat1s*mhat3s
                    rec_vir(2, 1) = rec_vir(2, 1) + tmp2*mhat2s*mhat1s
                    rec_vir(2, 2) = rec_vir(2, 2) + tmp1*(vterms*mhat2s*mhat2s - 1.d0)
                    rec_vir(2, 3) = rec_vir(2, 3) + tmp2*mhat2s*mhat3s
                    rec_vir(3, 1) = rec_vir(3, 1) + tmp2*mhat3s*mhat1s
                    rec_vir(3, 2) = rec_vir(3, 2) + tmp2*mhat3s*mhat2s
                    rec_vir(3, 3) = rec_vir(3, 3) + tmp1*(vterms*mhat3s*mhat3s - 1.d0)
#endif
                end if

                q(1, k3, k1, k2q) = eterm*q(1, k3, k1, k2q)
                q(2, k3, k1, k2q) = eterm*q(2, k3, k1, k2q)
            end do
        end do
    end do

    eer = 0.5d0*energy

#ifndef noVIRIAL
    rec_vir = HALF*rec_vir
#endif

    return

end subroutine scalar_sumrc_orthog

end module ew_recip
