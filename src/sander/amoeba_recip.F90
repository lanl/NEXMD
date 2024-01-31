#include "dprec.fh"
#include "assert.fh"

!---------------------------------------------------------
module amoeba_recip

    implicit none
    private

    integer, save :: do_flag
    logical, save :: perm_field_done = .FALSE.
    _REAL_, allocatable, save :: G_func(:, :, :), Qperm(:, :, :), &
        Q1(:, :, :), Q2(:, :, :), &
        theta1(:, :, :), theta2(:, :, :), theta3(:, :, :), &
        fractional_multipole(:, :), perm_F_field(:, :)

    integer, allocatable, save :: init_grid_ind(:, :)
    integer, parameter :: Max_Bspline_order = 25
# include "amoeba_mpole_index.h"
    public AM_RECIP_perm_field, AM_RECIP_dipole_field, AM_RECIP_ene_frc, &
        AM_RECIP_set_user_bit, AM_RECIP_allocate
#ifdef MPI
    public AM_RECIP_bcast
#endif
contains
!---------------------------------------------------------
#ifdef MPI
    subroutine AM_RECIP_bcast
        implicit none
        integer ierr

        include 'mpif.h'
# include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
    end subroutine AM_RECIP_bcast
#endif

    subroutine AM_RECIP_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then ! do in all cases
            do_flag = ibset(do_flag, USER_INDUCE_BIT)
            do_flag = ibset(do_flag, USER_POSTINDUCE_BIT)
        elseif (do_this == 2) then ! do the induction, not the post-induction
            do_flag = ibset(do_flag, USER_INDUCE_BIT)
            do_flag = ibclr(do_flag, USER_POSTINDUCE_BIT)
        elseif (do_this == 3) then ! do the post-induction, not the induction
            do_flag = ibclr(do_flag, USER_INDUCE_BIT)
            do_flag = ibset(do_flag, USER_POSTINDUCE_BIT)
        elseif (do_this == 0) then
            do_flag = ibclr(do_flag, USER_INDUCE_BIT)
            do_flag = ibclr(do_flag, USER_POSTINDUCE_BIT)
        else
            write (6, *) 'AM_RECIP_set_user_bit: bad value of user do_this'
            call mexit(6, 1)
        end if
    end subroutine AM_RECIP_set_user_bit
!---------------------------------------------------------
    subroutine AM_RECIP_allocate(numatoms)
        integer, intent(in) :: numatoms

        integer :: nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw, &
            ier, dr_order, Bspline_order
#include "ew_pme_recip.h"
#include "do_flag.h"

        call get_fftdims(nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw)
        allocate (Qperm(2*nfftdim1, nfftdim2, nfftdim3), stat=ier)
        REQUIRE(ier == 0)
        allocate (Q1(2*nfftdim1, nfftdim2, nfftdim3), stat=ier)
        REQUIRE(ier == 0)
        allocate (Q2(2*nfftdim1, nfftdim2, nfftdim3), stat=ier)
        REQUIRE(ier == 0)
        allocate (G_func(nfft3, nfftdim1, nfft2), stat=ier)
        REQUIRE(ier == 0)
        dr_order = 3
        Bspline_order = order
        if (Bspline_order < dr_order + 2) then
            write (6, *) 'Spline order too small. Must be at least ', dr_order + 2
            call mexit(6, 1)
        end if
        if (Bspline_order > Max_Bspline_order) then
            write (6, *) 'Bspline_order too big! Max = ', Max_Bspline_order
            call mexit(6, 1)
        end if
        allocate (theta1(0:dr_order, Bspline_order, numatoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (theta2(0:dr_order, Bspline_order, numatoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (theta3(0:dr_order, Bspline_order, numatoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (init_grid_ind(3, numatoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (fractional_multipole(10, numatoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (perm_F_field(20, numatoms), stat=ier)
        REQUIRE(ier == 0)
        do_flag = ibset(do_flag, VALID_BIT)
    end subroutine AM_RECIP_allocate
!---------------------------------------------------------
    subroutine AM_RECIP_perm_field(numatoms, crd, cart_dipole_field, x)
        use amoeba_multipoles, only : global_multipole
        use amoeba_induced, only : sq_polinv, is_polarizable
        use nblist, only : recip, volume
        use stack
        character(kind=1, len=19) :: routine = "AM_RECIP_perm_field"
        integer, intent(in) :: numatoms
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: cart_dipole_field(3, *)
        _REAL_, intent(in) :: x(*)
#include "ew_pme_recip.h"
#include "do_flag.h"
#include "def_time.h"

        integer :: nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw
        integer :: Bspline_order, dr_order, p_tmpy, p_alpha, p_beta, p_fftw, p_FdipF
        _REAL_ :: scratch(Max_Bspline_order*Max_Bspline_order)

        if (iand(do_flag, PROCEED_INDUCE) /= PROCEED_INDUCE) return

        call timer_start(TIME_BSPL)
        call get_fftdims(nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw)
        ! next get the bspline coeffs--later for MPI fill the select array
        ! for slabs i.e. those atoms landing in this proc's slab
        Bspline_order = order
        dr_order = 3

        ! fill theta1-3. these are saved for use in induction scf & final energy,frc
        call AM_RECIP_Bspline_fill(numatoms, dr_order, Bspline_order, &
            nfft1, nfft2, nfft3, crd, recip, &
            scratch, init_grid_ind, theta1, theta2, theta3)

        call AM_RECIP_global_to_fractional(numatoms, recip, nfft1, nfft2, nfft3, &
            global_multipole, fractional_multipole)

        call timer_stop_start(TIME_BSPL, TIME_FILLG)

        call AM_RECIP_perm_fill_grid(numatoms, Bspline_order, dr_order, &
            fractional_multipole, theta1, theta2, theta3, &
            init_grid_ind, nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, Q1)

        call timer_stop_start(TIME_FILLG, TIME_FFT)

        call get_stack(p_tmpy, 2*nfftdim1, routine)
        call get_stack(p_alpha, nfft1, routine)
        call get_stack(p_beta, nfft1, routine)
        call get_stack(p_fftw, sffw, routine)
        if (.not. rstack_ok) then
            deallocate (r_stack)
            allocate (r_stack(1:lastrst), stat=alloc_ier)
            call reassign_rstack(routine)
        end if
        REQUIRE(rstack_ok)

        call fft_backrc( &
            Q1, x(lfftable), r_stack(p_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
            r_stack(p_tmpy), r_stack(p_alpha), r_stack(p_beta) &
            )

        ! copy Q1 to save for final virial,force calculation
        call array_copy(Q1, Qperm, siz_q)
        call timer_stop_start(TIME_FFT, TIME_SCSUM)
        call AM_RECIP_get_recip_Gfunc(x(lprefac1), x(lprefac2), x(lprefac3), &
            recip, volume, ew_coeff, &
            nfft1, nfftdim1, nfft2, nfft3, G_func)
        call AM_RECIP_G_times_Q(nfft1, nfftdim1, nfft2, nfft3, G_func, Q1)
        call timer_stop_start(TIME_SCSUM, TIME_FFT)
        call fft_forwardrc( &
            Q1, x(lfftable), r_stack(p_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
            r_stack(p_tmpy), r_stack(p_alpha), r_stack(p_beta) &
            )
        call free_stack(p_fftw, routine)
        call free_stack(p_beta, routine)
        call free_stack(p_alpha, routine)
        call free_stack(p_tmpy, routine)

        call timer_stop_start(TIME_FFT, TIME_GRADS)

        call get_stack(p_FdipF, 3*numatoms, routine)
        if (.not. rstack_ok) then
            deallocate (r_stack)
            allocate (r_stack(1:lastrst), stat=alloc_ier)
            call reassign_rstack(routine)
        end if
        REQUIRE(rstack_ok)
        call AM_RECIP_get_perm_F_field(numatoms, is_polarizable, &
            init_grid_ind, dr_order, Bspline_order, &
            theta1, theta2, theta3, &
            nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, &
            Q1, perm_F_field, r_stack(p_FdipF))
        call AM_RECIP_Fdip_to_Cdip_field(numatoms, is_polarizable, nfft1, nfft2, nfft3, &
            recip, r_stack(p_FdipF), cart_dipole_field)
        call free_stack(p_FdipF, routine)
        perm_field_done = .true.

        call timer_stop(TIME_GRADS)

    end subroutine AM_RECIP_perm_field
!---------------------------------------------------------
    subroutine AM_RECIP_dipole_field(numatoms, x, &
        ind_dip1, ind_dip2, dip_field1, dip_field2)
        use amoeba_induced, only : is_polarizable
        use nblist, only : recip, volume
        use stack
        character(kind=1, len=21) :: routine = "AM_RECIP_dipole_field"
        integer, intent(in) :: numatoms
        _REAL_, intent(in) :: x(*)
        _REAL_, intent(in) :: ind_dip1(3, *), ind_dip2(3, *)
        _REAL_, intent(out) :: dip_field1(3, *), dip_field2(3, *)
#include "ew_pme_recip.h"
#include "do_flag.h"
#include "def_time.h"

        integer :: nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw
        integer :: Bspline_order, dr_order, p_tmpy, p_alpha, p_beta, p_fftw, &
            p_fdip1, p_fdip2, p_frac_field1, p_frac_field2

        if (iand(do_flag, PROCEED_INDUCE) /= PROCEED_INDUCE) return

        call timer_start(TIME_BSPL)
        call get_fftdims(nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw)
        Bspline_order = order
        dr_order = 3

        call get_stack(p_fdip1, 3*numatoms, routine)
        call get_stack(p_fdip2, 3*numatoms, routine)
        call get_stack(p_frac_field1, 3*numatoms, routine)
        call get_stack(p_frac_field2, 3*numatoms, routine)
        if (.not. rstack_ok) then
            deallocate (r_stack)
            allocate (r_stack(1:lastrst), stat=alloc_ier)
            call reassign_rstack(routine)
        end if
        REQUIRE(rstack_ok)
        call AM_RECIP_dip_to_frac_dip(numatoms, is_polarizable, nfft1, nfft2, nfft3, &
            recip, ind_dip1, ind_dip2, r_stack(p_fdip1), r_stack(p_fdip2))

        call timer_stop_start(TIME_BSPL, TIME_FILLG)
        call AM_RECIP_dipole_fill_grids(numatoms, Bspline_order, dr_order, &
            is_polarizable, &
            r_stack(p_fdip1), r_stack(p_fdip2), &
            theta1, theta2, theta3, &
            init_grid_ind, nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, Q1, Q2)
        call timer_stop_start(TIME_FILLG, TIME_FFT)
        call get_stack(p_tmpy, 2*nfftdim1, routine)
        call get_stack(p_alpha, nfft1, routine)
        call get_stack(p_beta, nfft1, routine)
        call get_stack(p_fftw, sffw, routine)
        if (.not. rstack_ok) then
            deallocate (r_stack)
            allocate (r_stack(1:lastrst), stat=alloc_ier)
            call reassign_rstack(routine)
        end if
        REQUIRE(rstack_ok)
        call fft_backrc( &
            Q1, x(lfftable), r_stack(p_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
            r_stack(p_tmpy), r_stack(p_alpha), r_stack(p_beta) &
            )
        call fft_backrc( &
            Q2, x(lfftable), r_stack(p_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
            r_stack(p_tmpy), r_stack(p_alpha), r_stack(p_beta) &
            )
        call timer_stop_start(TIME_FFT, TIME_SCSUM)
        call AM_RECIP_G_times_Q(nfft1, nfftdim1, nfft2, nfft3, G_func, Q1)
        call AM_RECIP_G_times_Q(nfft1, nfftdim1, nfft2, nfft3, G_func, Q2)
        call timer_stop_start(TIME_SCSUM, TIME_FFT)
        call fft_forwardrc( &
            Q1, x(lfftable), r_stack(p_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
            r_stack(p_tmpy), r_stack(p_alpha), r_stack(p_beta) &
            )
        call fft_forwardrc( &
            Q2, x(lfftable), r_stack(p_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
            r_stack(p_tmpy), r_stack(p_alpha), r_stack(p_beta) &
            )
        call free_stack(p_fftw, routine)
        call free_stack(p_beta, routine)
        call free_stack(p_alpha, routine)
        call free_stack(p_tmpy, routine)

        call timer_stop_start(TIME_FFT, TIME_GRADS)
        call AM_RECIP_get_2dip_F_fields(numatoms, is_polarizable, &
            init_grid_ind, dr_order, Bspline_order, &
            theta1, theta2, theta3, &
            nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, &
            Q1, Q2, r_stack(p_frac_field1), r_stack(p_frac_field2))
        call AM_RECIP_Fdip_to_Cdip_field(numatoms, is_polarizable, nfft1, nfft2, nfft3, &
            recip, r_stack(p_frac_field1), dip_field1)
        call AM_RECIP_Fdip_to_Cdip_field(numatoms, is_polarizable, nfft1, nfft2, nfft3, &
            recip, r_stack(p_frac_field2), dip_field2)
        call free_stack(p_frac_field2, routine)
        call free_stack(p_frac_field1, routine)
        call free_stack(p_fdip2, routine)
        call free_stack(p_fdip1, routine)
        call timer_stop(TIME_GRADS)
    end subroutine AM_RECIP_dipole_field
!---------------------------------------------------------
    subroutine AM_RECIP_ene_frc(numatoms, crd, x, ind_dip1, ind_dip2, &
        ene_perm, ene_ind, frc, virial)
        use amoeba_multipoles, only : torque_field, global_multipole
        use amoeba_induced, only : is_polarizable
        use amoeba_multipoles, only : coulomb_const_kcal_per_mole
        use nblist, only : recip, volume
        use stack
        character(kind=1, len=16) :: routine = "AM_RECIP_ene_frc"
        integer, intent(in) :: numatoms
        _REAL_, intent(in) :: crd(3, *), x(*)
        _REAL_, intent(in) :: ind_dip1(3, *), ind_dip2(3, *)
        _REAL_, intent(out) :: ene_perm, ene_ind
        _REAL_, intent(inout) :: frc(3, *), virial(3, 3)
#include "ew_pme_recip.h"
#include "do_flag.h"
#include "def_time.h"

        integer :: nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw, j, k
        integer :: Bspline_order, dr_order, p_tmpy, p_alpha, p_beta, p_fftw, &
            p_fdip1, p_fdip2, p_frac_field1, p_frac_field2, p_cdf

        ene_perm = 0.d0
        ene_ind = 0.d0
        if (iand(do_flag, PROCEED_POSTINDUCE) /= PROCEED_POSTINDUCE) return

        if (.not. perm_field_done) then ! this occurs if AM_RECIP_perm_field
            ! was not called since last call to here
            ! i.e. no induced dipoles
            ! we need permanent field for ene_perm
            call get_stack(p_cdf, 3*numatoms, routine)
            if (.not. rstack_ok) then
                deallocate (r_stack)
                allocate (r_stack(1:lastrst), stat=alloc_ier)
                call reassign_rstack(routine)
            end if
            REQUIRE(rstack_ok)
            call AM_RECIP_perm_field(numatoms, crd, r_stack(p_cdf), x)
            call free_stack(p_cdf, routine)
        end if
        perm_field_done = .false.  ! reset for next force call

        call timer_start(TIME_BSPL)
        call get_fftdims(nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw)
        Bspline_order = order
        dr_order = 3

        call get_stack(p_fdip1, 3*numatoms, routine)
        call get_stack(p_fdip2, 3*numatoms, routine)
        call get_stack(p_frac_field1, 20*numatoms, routine)
        call get_stack(p_frac_field2, 20*numatoms, routine)
        if (.not. rstack_ok) then
            deallocate (r_stack)
            allocate (r_stack(1:lastrst), stat=alloc_ier)
            call reassign_rstack(routine)
        end if
        REQUIRE(rstack_ok)
        call AM_RECIP_dip_to_frac_dip(numatoms, is_polarizable, nfft1, nfft2, nfft3, &
            recip, ind_dip1, ind_dip2, r_stack(p_fdip1), r_stack(p_fdip2))

        call timer_stop_start(TIME_BSPL, TIME_FILLG)
        call AM_RECIP_dipole_fill_grids(numatoms, Bspline_order, dr_order, &
            is_polarizable, &
            r_stack(p_fdip1), r_stack(p_fdip2), &
            theta1, theta2, theta3, &
            init_grid_ind, nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, Q1, Q2)
        call timer_stop_start(TIME_FILLG, TIME_FFT)
        call get_stack(p_tmpy, 2*nfftdim1, routine)
        call get_stack(p_alpha, nfft1, routine)
        call get_stack(p_beta, nfft1, routine)
        call get_stack(p_fftw, sffw, routine)
        if (.not. rstack_ok) then
            deallocate (r_stack)
            allocate (r_stack(1:lastrst), stat=alloc_ier)
            call reassign_rstack(routine)
        end if
        REQUIRE(rstack_ok)
        call fft_backrc( &
            Q1, x(lfftable), r_stack(p_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
            r_stack(p_tmpy), r_stack(p_alpha), r_stack(p_beta) &
            )
        call fft_backrc( &
            Q2, x(lfftable), r_stack(p_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
            r_stack(p_tmpy), r_stack(p_alpha), r_stack(p_beta) &
            )
        call timer_stop_start(TIME_FFT, TIME_SCSUM)
        ! recall Qperm was saved during perm field calc
        call AM_RECIP_scalar_sum(recip, ew_coeff, nfft1, nfftdim1, nfft2, nfft3, &
            G_func, Qperm, Q1, Q2, virial)
        call timer_stop_start(TIME_SCSUM, TIME_FFT)
        call fft_forwardrc( &
            Q1, x(lfftable), r_stack(p_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
            r_stack(p_tmpy), r_stack(p_alpha), r_stack(p_beta) &
            )
        call fft_forwardrc( &
            Q2, x(lfftable), r_stack(p_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, &
            r_stack(p_tmpy), r_stack(p_alpha), r_stack(p_beta) &
            )
        call free_stack(p_fftw, routine)
        call free_stack(p_beta, routine)
        call free_stack(p_alpha, routine)
        call free_stack(p_tmpy, routine)

        call timer_stop_start(TIME_FFT, TIME_GRADS)
        ! get the field due to the two sets of dipoles
        ! (up to grad of quadrupole terms i.e. octupolar order)
        call AM_RECIP_get_ind_F_fields(numatoms, &
            init_grid_ind, dr_order, Bspline_order, &
            theta1, theta2, theta3, &
            nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, &
            Q1, Q2, r_stack(p_frac_field1), r_stack(p_frac_field2))
        call AM_RECIP_perm_ene_grad(numatoms, nfft1, nfft2, nfft3, recip, &
            perm_F_field, fractional_multipole, ene_perm, frc)
        call AM_RECIP_ind_ene_grad(numatoms, nfft1, nfft2, nfft3, &
            is_polarizable, recip, &
            perm_F_field, fractional_multipole, &
            r_stack(p_frac_field1), r_stack(p_frac_field2), &
            r_stack(p_fdip1), r_stack(p_fdip2), &
            ene_ind, frc)
        call AM_RECIP_field_torque_virial( &
            numatoms, nfft1, nfft2, nfft3, recip, &
            perm_F_field, &
            r_stack(p_frac_field1), r_stack(p_frac_field2), &
            torque_field, virial, global_multipole, &
            ind_dip1, ind_dip2)
        call free_stack(p_frac_field2, routine)
        call free_stack(p_frac_field1, routine)
        call free_stack(p_fdip2, routine)
        call free_stack(p_fdip1, routine)
        call timer_stop(TIME_GRADS)
        do j = 1, 3
            do k = 1, 3
                virial(j, k) = coulomb_const_kcal_per_mole*virial(j, k)
            end do
        end do
    end subroutine AM_RECIP_ene_frc
!---------------------------------------------------------
    subroutine AM_RECIP_perm_fill_grid(numatoms, Bspline_order, dr_order, &
        Fmpole, theta1, theta2, theta3, &
        init_grid_ind, nfft1, nfft2, nfft3, &
        nfftdim1, nfftdim2, nfftdim3, Q)
        ! note hardwired for quadrupoles--
        integer, intent(in) :: numatoms, Bspline_order, dr_order
        _REAL_, intent(in) :: Fmpole(10, *)
        _REAL_, intent(in) :: theta1(0:dr_order, Bspline_order, *), &
            theta2(0:dr_order, Bspline_order, *), &
            theta3(0:dr_order, Bspline_order, *)
        integer, intent(in) :: init_grid_ind(3, *)
        integer, intent(in) :: nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3
        _REAL_, intent(out) :: Q(2*nfftdim1, nfftdim2, nfftdim3)

        integer :: n, igrd0, jgrd0, kgrd0, ith1, ith2, ith3, i, i0, j, j0, k, k0, ntot
        _REAL_ :: t0, t1, t2, u0, u1, u2, v0, v1, v2, term0, term1, term2

        ntot = 2*nfftdim1*nfftdim2*nfftdim3
        call zero_array(Q, ntot)
        do n = 1, numatoms
            igrd0 = init_grid_ind(1, n) !begin index in 1st direction
            jgrd0 = init_grid_ind(2, n) !begin index in 2nd direction
            kgrd0 = init_grid_ind(3, n) !begin index in 3rd direction
            k0 = kgrd0
            do ith3 = 1, Bspline_order
                k0 = k0 + 1
                k = k0 + 1 + (nfft3 - isign(nfft3, k0))/2
                j0 = jgrd0
                v0 = theta3(0, ith3, n) !theta3
                v1 = theta3(1, ith3, n) !1st deriv of theta3
                v2 = theta3(2, ith3, n) !2nd deriv of theta3
                do ith2 = 1, Bspline_order
                    j0 = j0 + 1
                    j = j0 + 1 + (nfft2 - isign(nfft2, j0))/2
                    u0 = theta2(0, ith2, n) !theta2
                    u1 = theta2(1, ith2, n) !1st deriv of theta2
                    u2 = theta2(2, ith2, n) !2nd deriv of theta2
                    i0 = igrd0
! hardwire our knowledge of layout of theta1,2,3 to pre-assemble factors
                    term0 = Fmpole(Ind_000, n)*u0*v0 + &
                        Fmpole(Ind_010, n)*u1*v0 + &
                        Fmpole(Ind_001, n)*u0*v1 + &
                        Fmpole(Ind_020, n)*u2*v0 + &
                        Fmpole(Ind_002, n)*u0*v2 + &
                        Fmpole(Ind_011, n)*u1*v1
                    term1 = Fmpole(Ind_100, n)*u0*v0 + &
                        Fmpole(Ind_110, n)*u1*v0 + &
                        Fmpole(Ind_101, n)*u0*v1
                    term2 = Fmpole(Ind_200, n)*u0*v0
                    do ith1 = 1, Bspline_order
                        i0 = i0 + 1
                        i = i0 + 1 + (nfft1 - isign(nfft1, i0))/2
                        t0 = theta1(0, ith1, n) !theta1
                        t1 = theta1(1, ith1, n) !1st deriv of theta1
                        t2 = theta1(2, ith1, n) !2nd deriv of theta1
                        Q(i, j, k) = Q(i, j, k) + term0*t0 + term1*t1 + term2*t2
                    end do
                end do !ith2 = 1,Bspline_order
            end do !ith3 = 1,Bspline_order
        end do !n = 1,numatoms

    end subroutine AM_RECIP_perm_fill_grid
!---------------------------------------------------------
    subroutine AM_RECIP_dipole_fill_grids(numatoms, Bspline_order, dr_order, &
        is_polarizable, fdip1, fdip2, &
        theta1, theta2, theta3, &
        init_grid_ind, nfft1, nfft2, nfft3, &
        nfftdim1, nfftdim2, nfftdim3, Q1, Q2)
        integer, intent(in) :: numatoms, Bspline_order, dr_order
        logical, intent(in) :: is_polarizable(*)
        _REAL_, intent(in) :: fdip1(3, *), fdip2(3, *)
        _REAL_, intent(in) :: theta1(0:dr_order, Bspline_order, *), &
            theta2(0:dr_order, Bspline_order, *), &
            theta3(0:dr_order, Bspline_order, *)
        integer, intent(in) :: init_grid_ind(3, *)
        integer, intent(in) :: nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3
        _REAL_, intent(out) :: Q1(2*nfftdim1, nfftdim2, nfftdim3)
        _REAL_, intent(out) :: Q2(2*nfftdim1, nfftdim2, nfftdim3)

        integer :: n, igrd0, jgrd0, kgrd0, ith1, ith2, ith3, i, i0, j, j0, k, k0, ntot
        _REAL_ :: t0, t1, u0, u1, v0, v1, term1_0, term1_1, term2_0, term2_1

        ntot = 2*nfftdim1*nfftdim2*nfftdim3
        call zero_array(Q1, ntot)
        call zero_array(Q2, ntot)
        do n = 1, numatoms
            if (is_polarizable(n)) then
                igrd0 = init_grid_ind(1, n) !begin index in 1st direction
                jgrd0 = init_grid_ind(2, n) !begin index in 2nd direction
                kgrd0 = init_grid_ind(3, n) !begin index in 3rd direction
                k0 = kgrd0
                do ith3 = 1, Bspline_order
                    k0 = k0 + 1
                    k = k0 + 1 + (nfft3 - isign(nfft3, k0))/2
                    j0 = jgrd0
                    v0 = theta3(0, ith3, n) !theta3
                    v1 = theta3(1, ith3, n) !1st deriv of theta3
                    do ith2 = 1, Bspline_order
                        j0 = j0 + 1
                        j = j0 + 1 + (nfft2 - isign(nfft2, j0))/2
                        u0 = theta2(0, ith2, n) !theta2
                        u1 = theta2(1, ith2, n) !1st deriv of theta2
                        i0 = igrd0
                        term1_0 = fdip1(2, n)*u1*v0 + fdip1(3, n)*u0*v1
                        term2_0 = fdip2(2, n)*u1*v0 + fdip2(3, n)*u0*v1
                        term1_1 = fdip1(1, n)*u0*v0
                        term2_1 = fdip2(1, n)*u0*v0
                        do ith1 = 1, Bspline_order
                            i0 = i0 + 1
                            i = i0 + 1 + (nfft1 - isign(nfft1, i0))/2
                            t0 = theta1(0, ith1, n) !theta1
                            t1 = theta1(1, ith1, n) !1st deriv of theta1
                            Q1(i, j, k) = Q1(i, j, k) + term1_0*t0 + term1_1*t1
                            Q2(i, j, k) = Q2(i, j, k) + term2_0*t0 + term2_1*t1
                        end do
                    end do !ith2 = 1,Bspline_order
                end do !ith3 = 1,Bspline_order
            end if ! is_polarizable(n)
        end do !n = 1,numatoms
    end subroutine AM_RECIP_dipole_fill_grids
!---------------------------------------------------------
    subroutine AM_RECIP_Bspline_fill(numatoms, dr_order, Bspline_order, &
        nfft1, nfft2, nfft3, crd, recip, &
        scratch, init_grid_ind, theta1, theta2, theta3)
        integer, intent(in)  :: numatoms, dr_order, Bspline_order, nfft1, nfft2, nfft3
        _REAL_, intent(in)   :: crd(3, *), recip(3, 3)
        _REAL_, intent(out)  :: scratch(Bspline_order, Bspline_order)
        integer, intent(out) :: init_grid_ind(3, *)
        _REAL_, intent(out)  :: theta1(0:dr_order, Bspline_order, *), &
            theta2(0:dr_order, Bspline_order, *), &
            theta3(0:dr_order, Bspline_order, *)

        integer n, ifr
        _REAL_ w, fr

        do n = 1, numatoms
            w = crd(1, n)*recip(1, 1) + crd(2, n)*recip(2, 1) + crd(3, n)*recip(3, 1)
            fr = nfft1*(w - dnint(w) + 0.5d0)
            ifr = int(fr)
            w = fr - ifr
            init_grid_ind(1, n) = ifr - Bspline_order
            call AM_RECIP_bspline_fill_gen(w, Bspline_order, scratch, dr_order, &
                theta1(:, :, n))
            w = crd(1, n)*recip(1, 2) + crd(2, n)*recip(2, 2) + crd(3, n)*recip(3, 2)
            fr = nfft2*(w - dnint(w) + 0.5d0)
            ifr = int(fr)
            w = fr - ifr
            init_grid_ind(2, n) = ifr - Bspline_order
            call AM_RECIP_bspline_fill_gen(w, Bspline_order, scratch, dr_order, &
                theta2(:, :, n))
            w = crd(1, n)*recip(1, 3) + crd(2, n)*recip(2, 3) + crd(3, n)*recip(3, 3)
            fr = nfft3*(w - dnint(w) + 0.5d0)
            ifr = int(fr)
            w = fr - ifr
            init_grid_ind(3, n) = ifr - Bspline_order
            call AM_RECIP_bspline_fill_gen(w, Bspline_order, scratch, dr_order, &
                theta3(:, :, n))
        end do
    end subroutine AM_RECIP_Bspline_fill
!---------------------------------------------------------
    subroutine AM_RECIP_bspline_fill_gen(w, spline_order, array, dr_order, new_array)
        implicit none
        _REAL_ w
        integer spline_order
        _REAL_, intent(out)  :: array(spline_order, spline_order)
        integer, intent(in)  :: dr_order
        _REAL_, intent(out)  :: new_array(0:dr_order, spline_order)

        integer k, j

! init spline_order 2
        array(2, 2) = w
        array(1, 2) = 1.d0 - w
! one pass to spline_order 3
        array(3, 3) = 0.5d0*w*array(2, 2)
        array(2, 3) = 0.5d0*((w + 1.d0)*array(1, 2) + (2.d0 - w)*array(2, 2))
        array(1, 3) = 0.5d0*(1.d0 - w)*array(1, 2)
! compute standard b-spline recursion
        do k = 4, spline_order
            call AM_RECIP_bspline_one_pass_recur(w, k, array(1, k - 1), array(1, k))
        end do
! do derivatives
        if (dr_order > 0) then
            call AM_RECIP_bspline_diff(array(1, spline_order - 1), spline_order)
            if (dr_order > 1) then
                call AM_RECIP_bspline_diff(array(1, spline_order - 2), spline_order - 1)
                call AM_RECIP_bspline_diff(array(1, spline_order - 2), spline_order)
                if (dr_order > 2) then
                    call AM_RECIP_bspline_diff(array(1, spline_order - 3), spline_order - 2)
                    call AM_RECIP_bspline_diff(array(1, spline_order - 3), spline_order - 1)
                    call AM_RECIP_bspline_diff(array(1, spline_order - 3), spline_order)
                    if (dr_order > 3) then
                        call AM_RECIP_bspline_diff(array(1, spline_order - 4), spline_order - 3)
                        call AM_RECIP_bspline_diff(array(1, spline_order - 4), spline_order - 2)
                        call AM_RECIP_bspline_diff(array(1, spline_order - 4), spline_order - 1)
                        call AM_RECIP_bspline_diff(array(1, spline_order - 4), spline_order)
                        if (dr_order > 4) then
                            call AM_RECIP_bspline_diff(array(1, spline_order - 5), spline_order - 4)
                            call AM_RECIP_bspline_diff(array(1, spline_order - 5), spline_order - 3)
                            call AM_RECIP_bspline_diff(array(1, spline_order - 5), spline_order - 2)
                            call AM_RECIP_bspline_diff(array(1, spline_order - 5), spline_order - 1)
                            call AM_RECIP_bspline_diff(array(1, spline_order - 5), spline_order)
                            if (dr_order > 5) then
                                write (6, *) 'spline derivs of order > 5 not implemented!'
                                call mexit(6, 1)
                            end if !( dr_order > 5 )then
                        end if !( dr_order > 4 )then
                    end if !( dr_order > 3 )then
                end if !( dr_order > 2 )then
            end if !( dr_order > 1 )then
        end if !( dr_order > 0 )then
! re-arrange array
        do k = 1, spline_order
            do j = 0, dr_order
                new_array(j, k) = array(k, spline_order - j)
            end do
        end do
        return
    end subroutine AM_RECIP_bspline_fill_gen
!---------------------------------------------------
    subroutine AM_RECIP_bspline_one_pass_recur(w, n, old, new)
        implicit none
        _REAL_ old(*), new(*), w
        integer n
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order
! RECURSION:  M_n(w) = (w/(n-1))*M_n-1(w)+((n-w)/(n-1))*M_n-1(w-1)
! i.e.   M_n(w+n-j) = ((w+n-j)/(n-1))*M_n-1(w+n-j)+((j-w)/(n-1))*M_n-1(w+n-j-1)
! i.e.   new(j) = ((w+n-j)/(n-1))*old(j-1) + ((j-w)/(n-1))*old(j)
! where old is array before one_pass (thus n->n-1) and new is array afterwards

        _REAL_ div
        integer j

        div = 1.d0/(n - 1)
        new(n) = div*w*old(n - 1)
        do j = 1, n - 2
            new(n - j) = div*((w + j)*old(n - j - 1) + (n - j - w)*old(n - j))
        end do
        new(1) = div*(1 - w)*old(1)
        return
    end subroutine AM_RECIP_bspline_one_pass_recur
!-------------------------------------------------------------
    subroutine AM_RECIP_bspline_diff(c, n)
        implicit none
        _REAL_ c(*)
        integer n
! Using notation from Essmann et al; w = u-[u] and
! array(j) = M_n(w + order - j)  where n is order
! DERIVATIVE:    d/dw M_n(w) = M_n-1(w) - M_n-1(w-1)
! i.e.   d/dw M_n(w+n-j) = M_n-1(w+n-j) - M_n-1(w+n-j-1)
! i.e.   new(j) = old(j-1) - old(j)
! where old is array before one_pass (thus n->n-1) and new is array afterwards
! do backwards to do in place

        integer j
        c(n) = c(n - 1)
        do j = n - 1, 2, -1
            c(j) = c(j - 1) - c(j)
        end do
        c(1) = -c(1)
        return
    end subroutine AM_RECIP_bspline_diff
!---------------------------------------------------------
    subroutine AM_RECIP_get_recip_Gfunc(prefac1, prefac2, prefac3, &
        recip, volume, ewald_coeff, &
        nfft1, nfftdim1, nfft2, nfft3, G)

        use constants, only : PI

        integer, intent(in) :: nfft1, nfftdim1, nfft2, nfft3
        _REAL_, intent(in) :: prefac1(nfft1), prefac2(nfft2), prefac3(nfft3)
        _REAL_, intent(in) :: recip(3, 3), volume, ewald_coeff
        _REAL_, intent(out) :: G(nfft3, nfftdim1, nfft2)

        integer :: k1, k2, k3, m1, m2, m3, k10, nf1, nf2, nf3
        _REAL_  :: fac, mhat1, mhat2, mhat3, msq, denom

        fac = pi**2/ewald_coeff**2
        nf1 = nfft1/2
        if (2*nf1 .lt. nfft1) nf1 = nf1 + 1
        nf2 = nfft2/2
        if (2*nf2 .lt. nfft2) nf2 = nf2 + 1
        nf3 = nfft3/2
        if (2*nf3 .lt. nfft3) nf3 = nf3 + 1

        do k2 = 1, nfft2
            m2 = k2 - 1
            if (k2 > nf2) m2 = k2 - 1 - nfft2
            do k3 = 1, nfft3
                m3 = k3 - 1
                if (k3 > nf3) m3 = k3 - 1 - nfft3
                k10 = 1
                if (k3 + k2 .eq. 2) k10 = 2
                do k1 = k10, nf1 + 1
                    m1 = k1 - 1
                    if (k1 > nf1) m1 = k1 - 1 - nfft1
                    mhat1 = recip(1, 1)*m1 + recip(1, 2)*m2 + recip(1, 3)*m3
                    mhat2 = recip(2, 1)*m1 + recip(2, 2)*m2 + recip(2, 3)*m3
                    mhat3 = recip(3, 1)*m1 + recip(3, 2)*m2 + recip(3, 3)*m3
                    msq = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3
                    denom = pi*volume*msq
                    G(k3, k1, k2) = exp(-fac*msq)*prefac1(k1)* &
                        prefac2(k2)*prefac3(k3)/denom
                end do
            end do
        end do

    end subroutine AM_RECIP_get_recip_Gfunc
!--------------------------------------------------------------
    subroutine AM_RECIP_G_times_Q(nfft1, nfftdim1, nfft2, nfft3, G, Q)
        implicit none
        integer, intent(in) :: nfft1, nfftdim1, nfft2, nfft3
        _REAL_, intent(in) ::  G(nfft3, nfftdim1, nfft2)
        _REAL_, intent(out) :: Q(2, nfft3, nfftdim1, nfft2)

        integer :: k1, k2, k3, k10, nf1

        nf1 = nfft1/2
        if (2*nf1 < nfft1) nf1 = nf1 + 1

!........Insist that Q(1,1,1,1) is set to 0 (true already for neutral)
        Q(1, 1, 1, 1) = 0.d0
        Q(2, 1, 1, 1) = 0.d0
        do k2 = 1, nfft2
            do k3 = 1, nfft3
                k10 = 1
                if (k3 + k2 == 2) k10 = 2
                do k1 = k10, nf1 + 1
                    Q(1, k3, k1, k2) = G(k3, k1, k2)*Q(1, k3, k1, k2)
                    Q(2, k3, k1, k2) = G(k3, k1, k2)*Q(2, k3, k1, k2)
                end do
            end do
        end do

    end subroutine AM_RECIP_G_times_Q
!--------------------------------------------------------------
    subroutine AM_RECIP_scalar_sum(recip, ewald_coeff, &
        nfft1, nfftdim1, nfft2, nfft3, &
        G, Q, Q1, Q2, virial)

        use constants, only : pi

        _REAL_, intent(in) :: recip(3, 3), ewald_coeff
        integer, intent(in) :: nfft1, nfftdim1, nfft2, nfft3
        _REAL_, intent(in) :: G(nfft3, nfftdim1, nfft2)
        _REAL_, intent(in) :: Q(2, nfft3, nfftdim1, nfft2)
        _REAL_, intent(inout) :: Q1(2, nfft3, nfftdim1, nfft2)
        _REAL_, intent(inout) :: Q2(2, nfft3, nfftdim1, nfft2)
        _REAL_, intent(out) :: virial(3, 3)

        integer :: k1, k2, k3, m1, m2, m3, k10, nf1, nf2, nf3
        _REAL_  :: fac, mhat1, mhat2, mhat3, msq, struc2, eterm, vterm, &
            q11, q12, q21, q22, tmp1, tmp2, mult, vxx, vxy, vxz, vyy, vyz, vzz

        fac = pi**2/ewald_coeff**2
        nf1 = nfft1/2
        if (2*nf1 .lt. nfft1) nf1 = nf1 + 1
        nf2 = nfft2/2
        if (2*nf2 .lt. nfft2) nf2 = nf2 + 1
        nf3 = nfft3/2
        if (2*nf3 .lt. nfft3) nf3 = nf3 + 1

        vxx = 0.d0
        vxy = 0.d0
        vxz = 0.d0
        vyy = 0.d0
        vyz = 0.d0
        vzz = 0.d0
        do k2 = 1, nfft2
            m2 = k2 - 1
            if (k2 > nf2) m2 = k2 - 1 - nfft2
            do k3 = 1, nfft3
                m3 = k3 - 1
                if (k3 > nf3) m3 = k3 - 1 - nfft3
                k10 = 1
                if (k3 + k2 .eq. 2) k10 = 2
                do k1 = k10, nf1 + 1
                    if (k1 > 1) then
                        mult = 2.d0
                    else
                        mult = 1.d0
                    end if
                    m1 = k1 - 1
                    if (k1 > nf1) m1 = k1 - 1 - nfft1
                    mhat1 = recip(1, 1)*m1 + recip(1, 2)*m2 + recip(1, 3)*m3
                    mhat2 = recip(2, 1)*m1 + recip(2, 2)*m2 + recip(2, 3)*m3
                    mhat3 = recip(3, 1)*m1 + recip(3, 2)*m2 + recip(3, 3)*m3
                    msq = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3
                    eterm = mult*G(k3, k1, k2)
                    vterm = 2.d0*(fac*msq + 1.d0)/msq
                    q11 = Q(1, k3, k1, k2) + Q1(1, k3, k1, k2)
                    q12 = Q(2, k3, k1, k2) + Q1(2, k3, k1, k2)
                    q21 = Q(1, k3, k1, k2) + Q2(1, k3, k1, k2)
                    q22 = Q(2, k3, k1, k2) + Q2(2, k3, k1, k2)
                    struc2 = q11*q21 + q12*q22
                    tmp1 = eterm*struc2
                    tmp2 = tmp1*vterm
                    vxx = vxx + tmp2*mhat1*mhat1 - tmp1
                    vxy = vxy + tmp2*mhat1*mhat2
                    vxz = vxz + tmp2*mhat1*mhat3
                    vyy = vyy + tmp2*mhat2*mhat2 - tmp1
                    vyz = vyz + tmp2*mhat2*mhat3
                    vzz = vzz + tmp2*mhat3*mhat3 - tmp1
                    ! transform Q1,Q2. No need to do Q--done in permanent field routine
                    Q1(1, k3, k1, k2) = G(k3, k1, k2)*Q1(1, k3, k1, k2)
                    Q1(2, k3, k1, k2) = G(k3, k1, k2)*Q1(2, k3, k1, k2)
                    Q2(1, k3, k1, k2) = G(k3, k1, k2)*Q2(1, k3, k1, k2)
                    Q2(2, k3, k1, k2) = G(k3, k1, k2)*Q2(2, k3, k1, k2)
                end do
            end do
        end do
        virial(1, 1) = virial(1, 1) + 0.5d0*vxx
        virial(1, 2) = virial(1, 2) + 0.5d0*vxy
        virial(2, 1) = virial(2, 1) + 0.5d0*vxy
        virial(1, 3) = virial(1, 3) + 0.5d0*vxz
        virial(3, 1) = virial(3, 1) + 0.5d0*vxz
        virial(2, 2) = virial(2, 2) + 0.5d0*vyy
        virial(2, 3) = virial(2, 3) + 0.5d0*vyz
        virial(3, 2) = virial(3, 2) + 0.5d0*vyz
        virial(3, 3) = virial(3, 3) + 0.5d0*vzz
    end subroutine AM_RECIP_scalar_sum
!--------------------------------------------------------------
    subroutine AM_RECIP_get_perm_F_field(numatoms, is_polarizable, &
        init_grid_ind, dr_order, Bspline_order, &
        theta1, theta2, theta3, &
        nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, &
        Q_p, Fperm_field, Fperm_dipfield)
        integer, intent(in) :: numatoms, init_grid_ind(3, *), dr_order, Bspline_order
        logical, intent(in) :: is_polarizable(*)
        _REAL_, intent(in) :: theta1(0:dr_order, Bspline_order, *), &
            theta2(0:dr_order, Bspline_order, *), &
            theta3(0:dr_order, Bspline_order, *)
        integer, intent(in) :: nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3
        _REAL_, intent(in) :: Q_p(2*nfftdim1, nfftdim2, nfftdim3)
        _REAL_, intent(out) :: Fperm_field(20, *), Fperm_dipfield(3, *)

! Note order of field elements
! 1 Field   corresponds to 000 in spline indices
! 2-4 -> F_x, F_y, F_z respectively or to 100,010,001 in indices
! 5-10 -> F_xx,F_yy,F_zz,F_xy,F_xz, F_yz resp. or to
!            200, 020, 002, 110, 101, 011 in indices
! 11-20 -> F_xxx,F_yyy,F_zzz,F_xxy,F_xxz,F_xyy,F_yyz,F_xzz,F_yzz,F_xyz or to
!          300,030,003,210,201,120,021,102,012,111  in indices
        _REAL_ :: tq_p, t_p(0:3), u(0:3), v(0:3)
        _REAL_ :: tu_p(10), tuv_p(20)
        integer :: m, n, i, i0, j, j0, k, k0, igrd0, jgrd0, kgrd0, ith1, ith2, ith3

        do n = 1, numatoms
            igrd0 = init_grid_ind(1, n) !begin index in 1st direction
            jgrd0 = init_grid_ind(2, n) !begin index in 2nd direction
            kgrd0 = init_grid_ind(3, n) !begin index in 3rd direction
            do m = 1, 20
                tuv_p(m) = 0.d0
            end do
            k0 = kgrd0
            do ith3 = 1, Bspline_order
                do m = 1, 10
                    tu_p(m) = 0.d0
                end do
                do m = 0, 3
                    v(m) = theta3(m, ith3, n)
                end do
                k0 = k0 + 1
                k = k0 + 1 + (nfft3 - isign(nfft3, k0))/2
                j0 = jgrd0
                do ith2 = 1, Bspline_order
                    j0 = j0 + 1
                    j = j0 + 1 + (nfft2 - isign(nfft2, j0))/2
                    i0 = igrd0
                    do m = 0, 3
                        u(m) = theta2(m, ith2, n)
                        t_p(m) = 0.d0
                    end do
                    do ith1 = 1, Bspline_order
                        i0 = i0 + 1
                        i = i0 + 1 + (nfft1 - isign(nfft1, i0))/2
                        tq_p = Q_p(i, j, k)
                        t_p(0) = t_p(0) + tq_p*theta1(0, ith1, n)
                        t_p(1) = t_p(1) + tq_p*theta1(1, ith1, n)
                        t_p(2) = t_p(2) + tq_p*theta1(2, ith1, n)
                        t_p(3) = t_p(3) + tq_p*theta1(3, ith1, n)
                    end do !ith1 = 1,Bspline_order
                    tu_p(Ind_00) = tu_p(Ind_00) + t_p(0)*u(0)
                    tu_p(Ind_10) = tu_p(Ind_10) + t_p(1)*u(0)
                    tu_p(Ind_01) = tu_p(Ind_01) + t_p(0)*u(1)
                    tu_p(Ind_20) = tu_p(Ind_20) + t_p(2)*u(0)
                    tu_p(Ind_02) = tu_p(Ind_02) + t_p(0)*u(2)
                    tu_p(Ind_11) = tu_p(Ind_11) + t_p(1)*u(1)
                    tu_p(Ind_30) = tu_p(Ind_30) + t_p(3)*u(0)
                    tu_p(Ind_03) = tu_p(Ind_03) + t_p(0)*u(3)
                    tu_p(Ind_21) = tu_p(Ind_21) + t_p(2)*u(1)
                    tu_p(Ind_12) = tu_p(Ind_12) + t_p(1)*u(2)
                end do !ith2 = 1,Spline_order
                tuv_p(Ind_000) = tuv_p(Ind_000) + tu_p(Ind_00)*v(0)
                tuv_p(Ind_100) = tuv_p(Ind_100) + tu_p(Ind_10)*v(0)
                tuv_p(Ind_010) = tuv_p(Ind_010) + tu_p(Ind_01)*v(0)
                tuv_p(Ind_001) = tuv_p(Ind_001) + tu_p(Ind_00)*v(1)
                tuv_p(Ind_200) = tuv_p(Ind_200) + tu_p(Ind_20)*v(0)
                tuv_p(Ind_020) = tuv_p(Ind_020) + tu_p(Ind_02)*v(0)
                tuv_p(Ind_002) = tuv_p(Ind_002) + tu_p(Ind_00)*v(2)
                tuv_p(Ind_110) = tuv_p(Ind_110) + tu_p(Ind_11)*v(0)
                tuv_p(Ind_101) = tuv_p(Ind_101) + tu_p(Ind_10)*v(1)
                tuv_p(Ind_011) = tuv_p(Ind_011) + tu_p(Ind_01)*v(1)
                tuv_p(Ind_300) = tuv_p(Ind_300) + tu_p(Ind_30)*v(0)
                tuv_p(Ind_030) = tuv_p(Ind_030) + tu_p(Ind_03)*v(0)
                tuv_p(Ind_003) = tuv_p(Ind_003) + tu_p(Ind_00)*v(3)
                tuv_p(Ind_210) = tuv_p(Ind_210) + tu_p(Ind_21)*v(0)
                tuv_p(Ind_201) = tuv_p(Ind_201) + tu_p(Ind_20)*v(1)
                tuv_p(Ind_120) = tuv_p(Ind_120) + tu_p(Ind_12)*v(0)
                tuv_p(Ind_021) = tuv_p(Ind_021) + tu_p(Ind_02)*v(1)
                tuv_p(Ind_102) = tuv_p(Ind_102) + tu_p(Ind_10)*v(2)
                tuv_p(Ind_012) = tuv_p(Ind_012) + tu_p(Ind_01)*v(2)
                tuv_p(Ind_111) = tuv_p(Ind_111) + tu_p(Ind_11)*v(1)
            end do !ith3 = 1,Spline_order
            do m = 1, 20
                Fperm_field(m, n) = tuv_p(m)
            end do
            if (is_polarizable(n)) then
                Fperm_dipfield(1, n) = tuv_p(Ind_100)
                Fperm_dipfield(2, n) = tuv_p(Ind_010)
                Fperm_dipfield(3, n) = tuv_p(Ind_001)
            else
                Fperm_dipfield(1, n) = 0.d0
                Fperm_dipfield(2, n) = 0.d0
                Fperm_dipfield(3, n) = 0.d0
            end if
        end do !n = 1,numatoms
    end subroutine AM_RECIP_get_perm_F_field
!--------------------------------------------------------------
    subroutine AM_RECIP_get_2dip_F_fields(numatoms, is_polarizable, &
        init_grid_ind, dr_order, Bspline_order, &
        theta1, theta2, theta3, &
        nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, Q1, Q2, &
        Fdip_field1, Fdip_field2)
        integer, intent(in) :: numatoms, init_grid_ind(3, *), dr_order, Bspline_order
        logical, intent(in) :: is_polarizable(*)
        _REAL_, intent(in) :: theta1(0:dr_order, Bspline_order, *), &
            theta2(0:dr_order, Bspline_order, *), &
            theta3(0:dr_order, Bspline_order, *)
        integer, intent(in) :: nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3
        _REAL_, intent(in) :: Q1(2*nfftdim1, nfftdim2, nfftdim3)
        _REAL_, intent(in) :: Q2(2*nfftdim1, nfftdim2, nfftdim3)
        _REAL_, intent(out) :: Fdip_field1(3, *), Fdip_field2(3, *)

        _REAL_ :: tq1, tq2, t0_1, t0_2, t1_1, t1_2, u0, u1, v0, v1, &
            tu00_1, tu00_2, tu10_1, tu10_2, tu01_1, tu01_2, &
            tuv001_1, tuv001_2, tuv010_1, tuv010_2, tuv100_1, tuv100_2
        integer :: n, i, i0, j, j0, k, k0, igrd0, jgrd0, kgrd0, ith1, ith2, ith3

        do n = 1, numatoms
            if (is_polarizable(n)) then
                igrd0 = init_grid_ind(1, n) !begin index in 1st direction
                jgrd0 = init_grid_ind(2, n) !begin index in 2nd direction
                kgrd0 = init_grid_ind(3, n) !begin index in 3rd direction
                tuv001_1 = 0.d0
                tuv010_1 = 0.d0
                tuv100_1 = 0.d0
                tuv001_2 = 0.d0
                tuv010_2 = 0.d0
                tuv100_2 = 0.d0
                k0 = kgrd0
                do ith3 = 1, Bspline_order
                    v0 = theta3(0, ith3, n) !theta3
                    v1 = theta3(1, ith3, n) !1st deriv of theta3
                    tu00_1 = 0.d0
                    tu10_1 = 0.d0
                    tu01_1 = 0.d0
                    tu00_2 = 0.d0
                    tu10_2 = 0.d0
                    tu01_2 = 0.d0
                    k0 = k0 + 1
                    k = k0 + 1 + (nfft3 - isign(nfft3, k0))/2
                    j0 = jgrd0
                    do ith2 = 1, Bspline_order
                        j0 = j0 + 1
                        j = j0 + 1 + (nfft2 - isign(nfft2, j0))/2
                        i0 = igrd0
                        u0 = theta2(0, ith2, n) !theta2
                        u1 = theta2(1, ith2, n) !1st deriv of theta2
                        t0_1 = 0.d0
                        t1_1 = 0.d0
                        t0_2 = 0.d0
                        t1_2 = 0.d0
                        do ith1 = 1, Bspline_order
                            i0 = i0 + 1
                            i = i0 + 1 + (nfft1 - isign(nfft1, i0))/2
                            tq1 = Q1(i, j, k)
                            tq2 = Q2(i, j, k)
                            t0_1 = t0_1 + tq1*theta1(0, ith1, n)
                            t1_1 = t1_1 + tq1*theta1(1, ith1, n)
                            t0_2 = t0_2 + tq2*theta1(0, ith1, n)
                            t1_2 = t1_2 + tq2*theta1(1, ith1, n)
                        end do !ith1 = 1,Bspline_order
                        tu00_1 = tu00_1 + t0_1*u0
                        tu10_1 = tu10_1 + t1_1*u0
                        tu01_1 = tu01_1 + t0_1*u1
                        tu00_2 = tu00_2 + t0_2*u0
                        tu10_2 = tu10_2 + t1_2*u0
                        tu01_2 = tu01_2 + t0_2*u1
                    end do !ith2 = 1,Spline_order
                    tuv100_1 = tuv100_1 + tu10_1*v0
                    tuv010_1 = tuv010_1 + tu01_1*v0
                    tuv001_1 = tuv001_1 + tu00_1*v1
                    tuv100_2 = tuv100_2 + tu10_2*v0
                    tuv010_2 = tuv010_2 + tu01_2*v0
                    tuv001_2 = tuv001_2 + tu00_2*v1
                end do !ith3 = 1,Spline_order
                Fdip_field1(1, n) = tuv100_1
                Fdip_field1(2, n) = tuv010_1
                Fdip_field1(3, n) = tuv001_1
                Fdip_field2(1, n) = tuv100_2
                Fdip_field2(2, n) = tuv010_2
                Fdip_field2(3, n) = tuv001_2
            end if !is_polarizable(n)
        end do !n = 1,numatoms
    end subroutine AM_RECIP_get_2dip_F_fields
!--------------------------------------------------------------
    subroutine AM_RECIP_get_ind_F_fields(numatoms, &
        init_grid_ind, dr_order, Bspline_order, &
        theta1, theta2, theta3, &
        nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, &
        Q1, Q2, Fdip_field1, Fdip_field2)
        integer, intent(in) :: numatoms, init_grid_ind(3, *), dr_order, Bspline_order
        _REAL_, intent(in) :: theta1(0:dr_order, Bspline_order, *), &
            theta2(0:dr_order, Bspline_order, *), &
            theta3(0:dr_order, Bspline_order, *)
        integer, intent(in) :: nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3
        _REAL_, intent(in) :: Q1(2*nfftdim1, nfftdim2, nfftdim3)
        _REAL_, intent(in) :: Q2(2*nfftdim1, nfftdim2, nfftdim3)
        _REAL_, intent(out) :: Fdip_field1(20, *), Fdip_field2(20, *)

! Note order of field elements
! 1 Field   corresponds to 000 in spline indices
! 2-4 -> F_x, F_y, F_z respectively or to 100,010,001 in indices
! 5-10 -> F_xx,F_yy,F_zz,F_xy,F_xz, F_yz resp. or to
!            200, 020, 002, 110, 101, 011 in indices
! 11-20 -> F_xxx,F_yyy,F_zzz,F_xxy,F_xxz,F_xyy,F_yyz,F_xzz,F_yzz,F_xyz or to
!          300,030,003,210,201,120,021,102,012,111  in indices
        _REAL_ :: tq_1, tq_2, t_1(0:3), t_2(0:3), u(0:3), v(0:3)
        _REAL_ :: tu_1(10), tu_2(10), tuv_1(20), tuv_2(20)
        integer :: m, n, i, i0, j, j0, k, k0, igrd0, jgrd0, kgrd0, ith1, ith2, ith3

        do n = 1, numatoms
            igrd0 = init_grid_ind(1, n) !begin index in 1st direction
            jgrd0 = init_grid_ind(2, n) !begin index in 2nd direction
            kgrd0 = init_grid_ind(3, n) !begin index in 3rd direction
            do m = 1, 20
                tuv_1(m) = 0.d0
                tuv_2(m) = 0.d0
            end do
            k0 = kgrd0
            do ith3 = 1, Bspline_order
                do m = 1, 10
                    tu_1(m) = 0.d0
                    tu_2(m) = 0.d0
                end do
                do m = 0, 3
                    v(m) = theta3(m, ith3, n)
                end do
                k0 = k0 + 1
                k = k0 + 1 + (nfft3 - isign(nfft3, k0))/2
                j0 = jgrd0
                do ith2 = 1, Bspline_order
                    j0 = j0 + 1
                    j = j0 + 1 + (nfft2 - isign(nfft2, j0))/2
                    i0 = igrd0
                    do m = 0, 3
                        u(m) = theta2(m, ith2, n)
                        t_1(m) = 0.d0
                        t_2(m) = 0.d0
                    end do
                    do ith1 = 1, Bspline_order
                        i0 = i0 + 1
                        i = i0 + 1 + (nfft1 - isign(nfft1, i0))/2
                        tq_1 = Q1(i, j, k)
                        tq_2 = Q2(i, j, k)
                        t_1(0) = t_1(0) + tq_1*theta1(0, ith1, n)
                        t_1(1) = t_1(1) + tq_1*theta1(1, ith1, n)
                        t_1(2) = t_1(2) + tq_1*theta1(2, ith1, n)
                        t_1(3) = t_1(3) + tq_1*theta1(3, ith1, n)
                        t_2(0) = t_2(0) + tq_2*theta1(0, ith1, n)
                        t_2(1) = t_2(1) + tq_2*theta1(1, ith1, n)
                        t_2(2) = t_2(2) + tq_2*theta1(2, ith1, n)
                        t_2(3) = t_2(3) + tq_2*theta1(3, ith1, n)
                    end do !ith1 = 1,Bspline_order
                    tu_1(Ind_00) = tu_1(Ind_00) + t_1(0)*u(0)
                    tu_1(Ind_10) = tu_1(Ind_10) + t_1(1)*u(0)
                    tu_1(Ind_01) = tu_1(Ind_01) + t_1(0)*u(1)
                    tu_1(Ind_20) = tu_1(Ind_20) + t_1(2)*u(0)
                    tu_1(Ind_02) = tu_1(Ind_02) + t_1(0)*u(2)
                    tu_1(Ind_11) = tu_1(Ind_11) + t_1(1)*u(1)
                    tu_1(Ind_30) = tu_1(Ind_30) + t_1(3)*u(0)
                    tu_1(Ind_03) = tu_1(Ind_03) + t_1(0)*u(3)
                    tu_1(Ind_21) = tu_1(Ind_21) + t_1(2)*u(1)
                    tu_1(Ind_12) = tu_1(Ind_12) + t_1(1)*u(2)

                    tu_2(Ind_00) = tu_2(Ind_00) + t_2(0)*u(0)
                    tu_2(Ind_10) = tu_2(Ind_10) + t_2(1)*u(0)
                    tu_2(Ind_01) = tu_2(Ind_01) + t_2(0)*u(1)
                    tu_2(Ind_20) = tu_2(Ind_20) + t_2(2)*u(0)
                    tu_2(Ind_02) = tu_2(Ind_02) + t_2(0)*u(2)
                    tu_2(Ind_11) = tu_2(Ind_11) + t_2(1)*u(1)
                    tu_2(Ind_30) = tu_2(Ind_30) + t_2(3)*u(0)
                    tu_2(Ind_03) = tu_2(Ind_03) + t_2(0)*u(3)
                    tu_2(Ind_21) = tu_2(Ind_21) + t_2(2)*u(1)
                    tu_2(Ind_12) = tu_2(Ind_12) + t_2(1)*u(2)

                end do !ith2 = 1,Spline_order
                tuv_1(Ind_000) = tuv_1(Ind_000) + tu_1(Ind_00)*v(0)
                tuv_1(Ind_100) = tuv_1(Ind_100) + tu_1(Ind_10)*v(0)
                tuv_1(Ind_010) = tuv_1(Ind_010) + tu_1(Ind_01)*v(0)
                tuv_1(Ind_001) = tuv_1(Ind_001) + tu_1(Ind_00)*v(1)
                tuv_1(Ind_200) = tuv_1(Ind_200) + tu_1(Ind_20)*v(0)
                tuv_1(Ind_020) = tuv_1(Ind_020) + tu_1(Ind_02)*v(0)
                tuv_1(Ind_002) = tuv_1(Ind_002) + tu_1(Ind_00)*v(2)
                tuv_1(Ind_110) = tuv_1(Ind_110) + tu_1(Ind_11)*v(0)
                tuv_1(Ind_101) = tuv_1(Ind_101) + tu_1(Ind_10)*v(1)
                tuv_1(Ind_011) = tuv_1(Ind_011) + tu_1(Ind_01)*v(1)
                tuv_1(Ind_300) = tuv_1(Ind_300) + tu_1(Ind_30)*v(0)
                tuv_1(Ind_030) = tuv_1(Ind_030) + tu_1(Ind_03)*v(0)
                tuv_1(Ind_003) = tuv_1(Ind_003) + tu_1(Ind_00)*v(3)
                tuv_1(Ind_210) = tuv_1(Ind_210) + tu_1(Ind_21)*v(0)
                tuv_1(Ind_201) = tuv_1(Ind_201) + tu_1(Ind_20)*v(1)
                tuv_1(Ind_120) = tuv_1(Ind_120) + tu_1(Ind_12)*v(0)
                tuv_1(Ind_021) = tuv_1(Ind_021) + tu_1(Ind_02)*v(1)
                tuv_1(Ind_102) = tuv_1(Ind_102) + tu_1(Ind_10)*v(2)
                tuv_1(Ind_012) = tuv_1(Ind_012) + tu_1(Ind_01)*v(2)
                tuv_1(Ind_111) = tuv_1(Ind_111) + tu_1(Ind_11)*v(1)

                tuv_2(Ind_000) = tuv_2(Ind_000) + tu_2(Ind_00)*v(0)
                tuv_2(Ind_100) = tuv_2(Ind_100) + tu_2(Ind_10)*v(0)
                tuv_2(Ind_010) = tuv_2(Ind_010) + tu_2(Ind_01)*v(0)
                tuv_2(Ind_001) = tuv_2(Ind_001) + tu_2(Ind_00)*v(1)
                tuv_2(Ind_200) = tuv_2(Ind_200) + tu_2(Ind_20)*v(0)
                tuv_2(Ind_020) = tuv_2(Ind_020) + tu_2(Ind_02)*v(0)
                tuv_2(Ind_002) = tuv_2(Ind_002) + tu_2(Ind_00)*v(2)
                tuv_2(Ind_110) = tuv_2(Ind_110) + tu_2(Ind_11)*v(0)
                tuv_2(Ind_101) = tuv_2(Ind_101) + tu_2(Ind_10)*v(1)
                tuv_2(Ind_011) = tuv_2(Ind_011) + tu_2(Ind_01)*v(1)
                tuv_2(Ind_300) = tuv_2(Ind_300) + tu_2(Ind_30)*v(0)
                tuv_2(Ind_030) = tuv_2(Ind_030) + tu_2(Ind_03)*v(0)
                tuv_2(Ind_003) = tuv_2(Ind_003) + tu_2(Ind_00)*v(3)
                tuv_2(Ind_210) = tuv_2(Ind_210) + tu_2(Ind_21)*v(0)
                tuv_2(Ind_201) = tuv_2(Ind_201) + tu_2(Ind_20)*v(1)
                tuv_2(Ind_120) = tuv_2(Ind_120) + tu_2(Ind_12)*v(0)
                tuv_2(Ind_021) = tuv_2(Ind_021) + tu_2(Ind_02)*v(1)
                tuv_2(Ind_102) = tuv_2(Ind_102) + tu_2(Ind_10)*v(2)
                tuv_2(Ind_012) = tuv_2(Ind_012) + tu_2(Ind_01)*v(2)
                tuv_2(Ind_111) = tuv_2(Ind_111) + tu_2(Ind_11)*v(1)

            end do !ith3 = 1,Spline_order
            do m = 1, 20
                Fdip_field1(m, n) = tuv_1(m)
                Fdip_field2(m, n) = tuv_2(m)
            end do
        end do !n = 1,numatoms
    end subroutine AM_RECIP_get_ind_F_fields
!--------------------------------------------------------------
    subroutine AM_RECIP_perm_ene_grad(numatoms, nfft1, nfft2, nfft3, recip, &
        F_field, F_mpole, ene_perm, frc)
        use amoeba_multipoles, only : coulomb_const_kcal_per_mole
        integer, intent(in) :: numatoms, nfft1, nfft2, nfft3
        _REAL_, intent(in) :: recip(3, 3), F_field(20, *), F_mpole(10, *)
        _REAL_, intent(out) :: ene_perm
        _REAL_, intent(inout) :: frc(3, *)

        _REAL_ :: f1, f2, f3, dfx, dfy, dfz
        _REAL_ :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)
        integer :: k, n, j1, j2, j3

        call AM_RECIP_xform_matrices(nfft1, nfft2, nfft3, recip, &
            mpole_xform_3x3, field_xform_3x3)
        do n = 1, numatoms
            f1 = 0.d0
            f2 = 0.d0
            f3 = 0.d0
            do k = 1, 10
                ene_perm = ene_perm + F_mpole(k, n)*F_field(k, n)
                j1 = deriv1(k)
                j2 = deriv2(k)
                j3 = deriv3(k)
                f1 = f1 + F_mpole(k, n)*F_field(j1, n)
                f2 = f2 + F_mpole(k, n)*F_field(j2, n)
                f3 = f3 + F_mpole(k, n)*F_field(j3, n)
            end do
! force is negative of gradient
! transform from scaled fractional to cartesian--same as fields
            dfx = field_xform_3x3(1, 1)*f1 + field_xform_3x3(1, 2)*f2 + &
                field_xform_3x3(1, 3)*f3
            dfy = field_xform_3x3(2, 1)*f1 + field_xform_3x3(2, 2)*f2 + &
                field_xform_3x3(2, 3)*f3
            dfz = field_xform_3x3(3, 1)*f1 + field_xform_3x3(3, 2)*f2 + &
                field_xform_3x3(3, 3)*f3
            frc(1, n) = frc(1, n) - coulomb_const_kcal_per_mole*dfx
            frc(2, n) = frc(2, n) - coulomb_const_kcal_per_mole*dfy
            frc(3, n) = frc(3, n) - coulomb_const_kcal_per_mole*dfz
        end do !n = 1,numatoms
        ene_perm = 0.5d0*coulomb_const_kcal_per_mole*ene_perm
    end subroutine AM_RECIP_perm_ene_grad
!--------------------------------------------------------------
    subroutine AM_RECIP_ind_ene_grad(numatoms, nfft1, nfft2, nfft3, &
        is_polarizable, recip, &
        F_field, F_mpole, F_dip1_field, F_dip2_field, &
        F_dip1, F_dip2, ene_ind, frc)
        use amoeba_multipoles, only : coulomb_const_kcal_per_mole
        integer, intent(in) :: numatoms, nfft1, nfft2, nfft3
        logical, intent(in) :: is_polarizable(*)
        _REAL_, intent(in) :: recip(3, 3), F_field(20, *), F_mpole(10, *)
        _REAL_, intent(in) :: F_dip1_field(20, *), F_dip2_field(20, *)
        _REAL_, intent(in) :: F_dip1(3, *), F_dip2(3, *)
        _REAL_, intent(out) :: ene_ind
        _REAL_, intent(inout) :: frc(3, *)

        _REAL_ :: f1, f2, f3, dfx, dfy, dfz
        _REAL_ :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)
        integer :: k, n, j1, j2, j3

        call AM_RECIP_xform_matrices(nfft1, nfft2, nfft3, recip, &
            mpole_xform_3x3, field_xform_3x3)
        do n = 1, numatoms
            f1 = 0.d0
            f2 = 0.d0
            f3 = 0.d0
            ! add the force on mpoles due to ave of dip fields
            do k = 1, 10
                j1 = deriv1(k)
                j2 = deriv2(k)
                j3 = deriv3(k)
                f1 = f1 + 0.5d0*F_mpole(k, n)*(F_dip1_field(j1, n) + F_dip2_field(j1, n))
                f2 = f2 + 0.5d0*F_mpole(k, n)*(F_dip1_field(j2, n) + F_dip2_field(j2, n))
                f3 = f3 + 0.5d0*F_mpole(k, n)*(F_dip1_field(j3, n) + F_dip2_field(j3, n))
            end do
            if (is_polarizable(n)) then
                ! induced energy is interaction of direct dipoles with perm field
                ene_ind = ene_ind + F_dip1(1, n)*F_field(Ind_100, n) + &
                    F_dip1(2, n)*F_field(Ind_010, n) + &
                    F_dip1(3, n)*F_field(Ind_001, n)
                do k = 2, 4
                    j1 = deriv1(k)
                    j2 = deriv2(k)
                    j3 = deriv3(k)
                    ! add the force on ave of dipoles due to perm field
                    f1 = f1 + 0.5d0*(F_dip1(k - 1, n) + F_dip2(k - 1, n))*F_field(j1, n)
                    f2 = f2 + 0.5d0*(F_dip1(k - 1, n) + F_dip2(k - 1, n))*F_field(j2, n)
                    f3 = f3 + 0.5d0*(F_dip1(k - 1, n) + F_dip2(k - 1, n))*F_field(j3, n)
                    ! next the forces of dips on each other
                    f1 = f1 + 0.5d0*(F_dip1(k - 1, n)*F_dip2_field(j1, n) + &
                        F_dip2(k - 1, n)*F_dip1_field(j1, n))
                    f2 = f2 + 0.5d0*(F_dip1(k - 1, n)*F_dip2_field(j2, n) + &
                        F_dip2(k - 1, n)*F_dip1_field(j2, n))
                    f3 = f3 + 0.5d0*(F_dip1(k - 1, n)*F_dip2_field(j3, n) + &
                        F_dip2(k - 1, n)*F_dip1_field(j3, n))
                end do
            end if !is_polarizable(n)
! force is negative of gradient
! transform from scaled fractional to cartesian--same as fields
            dfx = field_xform_3x3(1, 1)*f1 + field_xform_3x3(1, 2)*f2 + &
                field_xform_3x3(1, 3)*f3
            dfy = field_xform_3x3(2, 1)*f1 + field_xform_3x3(2, 2)*f2 + &
                field_xform_3x3(2, 3)*f3
            dfz = field_xform_3x3(3, 1)*f1 + field_xform_3x3(3, 2)*f2 + &
                field_xform_3x3(3, 3)*f3
            frc(1, n) = frc(1, n) - coulomb_const_kcal_per_mole*dfx
            frc(2, n) = frc(2, n) - coulomb_const_kcal_per_mole*dfy
            frc(3, n) = frc(3, n) - coulomb_const_kcal_per_mole*dfz
        end do !n = 1,numatoms
        ene_ind = 0.5d0*coulomb_const_kcal_per_mole*ene_ind
    end subroutine AM_RECIP_ind_ene_grad
!--------------------------------------------------------------
    subroutine AM_RECIP_global_to_fractional(numatoms, recip, &
        nfft1, nfft2, nfft3, glob_mpole, frac_mpole)
        integer, intent(in) :: numatoms
        _REAL_, intent(in) :: recip(3, 3)
        integer, intent(in) :: nfft1, nfft2, nfft3
        _REAL_, intent(in) :: glob_mpole(10, *)
        _REAL_, intent(out) :: frac_mpole(10, *)

        _REAL_ :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)
        integer order, dimxy, n
        _REAL_ Mpole_xy(MAXMP*MAXMP)
        order = 10
        dimxy = 10
! first get mpole_xform_3x3
        call AM_RECIP_xform_matrices(nfft1, nfft2, nfft3, recip, &
            mpole_xform_3x3, field_xform_3x3)
        call XFORM_MPOLE_matrix(mpole_xform_3x3, Mpole_xy, order)
        do n = 1, numatoms
            call XFORM_MPOLE(Mpole_xy, dimxy, glob_mpole(:, n), frac_mpole(:, n), order)
        end do

    end subroutine AM_RECIP_global_to_fractional
!---------------------------------------------------------
    subroutine AM_RECIP_field_torque_virial( &
        numatoms, nfft1, nfft2, nfft3, recip, &
        F_field, F_dip1_field, F_dip2_field, torque_field, &
        virial, mpole_p, dipole1, dipole2)

        integer, intent(in) :: numatoms, nfft1, nfft2, nfft3
        _REAL_, intent(in)  :: recip(3, 3), F_field(20, *), &
            F_dip1_field(20, *), F_dip2_field(20, *)
        _REAL_, intent(in) :: mpole_p(10, *), dipole1(3, *), dipole2(3, *)
        _REAL_, intent(inout) :: torque_field(10, *)
        _REAL_, intent(inout) :: virial(3, 3)

        _REAL_ :: field_p(10), field_d1(10), field_d2(10)
        _REAL_ :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)
        integer order, dimxy
        _REAL_ Field_xy(MAXMP*MAXMP)
        integer :: i, j, k, n

        order = 10
        dimxy = 10
        ! transform fields from fractional to cartesian
        call AM_RECIP_xform_matrices(nfft1, nfft2, nfft3, recip, &
            mpole_xform_3x3, field_xform_3x3)
        ! get the higher order terms
        call XFORM_MPOLE_field_matrix(field_xform_3x3, Field_xy, order)
        do n = 1, numatoms
            call XFORM_FIELD(Field_xy, dimxy, F_field(:, n), field_p, order)
            call XFORM_FIELD(Field_xy, dimxy, F_dip1_field(:, n), field_d1, order)
            call XFORM_FIELD(Field_xy, dimxy, F_dip2_field(:, n), field_d2, order)
            ! torque is due to permanent field plus ave of dip fields
            do i = 1, 10
                torque_field(i, n) = torque_field(i, n) + &
                    field_p(i) + 0.5d0*(field_d1(i) + field_d2(i))
            end do
            ! extra virial terms due to dipolar contributions
            do j = 1, 3
                do k = 1, 3
                    virial(j, k) = virial(j, k) - &
                        (torque_field(j + 1, n)*mpole_p(k + 1, n) + &
                        0.5d0*field_p(j + 1)*(dipole1(k, n) + dipole2(k, n)) + &
                        0.5d0*(field_d1(j + 1)*dipole2(k, n) + field_d2(j + 1)*dipole1(k, n)))
                end do
            end do
            ! quadrupolar contribs
            virial(1, 1) = virial(1, 1) - (2.d0*mpole_p(5, n)*torque_field(5, n) + &
                mpole_p(8, n)*torque_field(8, n) + mpole_p(9, n)*torque_field(9, n))
            virial(2, 2) = virial(2, 2) - (2.d0*mpole_p(6, n)*torque_field(6, n) + &
                mpole_p(8, n)*torque_field(8, n) + mpole_p(10, n)*torque_field(10, n))
            virial(3, 3) = virial(3, 3) - (2.d0*mpole_p(7, n)*torque_field(7, n) + &
                mpole_p(9, n)*torque_field(9, n) + mpole_p(10, n)*torque_field(10, n))
            virial(1, 2) = virial(1, 2) - (2.d0*mpole_p(6, n)*torque_field(8, n) + &
                mpole_p(8, n)*torque_field(5, n) + mpole_p(10, n)*torque_field(9, n))
            virial(2, 1) = virial(2, 1) - (2.d0*mpole_p(5, n)*torque_field(8, n) + &
                mpole_p(8, n)*torque_field(6, n) + mpole_p(9, n)*torque_field(10, n))
            virial(1, 3) = virial(1, 3) - (2.d0*mpole_p(7, n)*torque_field(9, n) + &
                mpole_p(9, n)*torque_field(5, n) + mpole_p(10, n)*torque_field(8, n))
            virial(3, 1) = virial(3, 1) - (2.d0*mpole_p(5, n)*torque_field(9, n) + &
                mpole_p(8, n)*torque_field(10, n) + mpole_p(9, n)*torque_field(7, n))
            virial(2, 3) = virial(2, 3) - (2.d0*mpole_p(7, n)*torque_field(10, n) + &
                mpole_p(9, n)*torque_field(8, n) + mpole_p(10, n)*torque_field(6, n))
            virial(3, 2) = virial(3, 2) - (2.d0*mpole_p(6, n)*torque_field(10, n) + &
                mpole_p(8, n)*torque_field(9, n) + mpole_p(10, n)*torque_field(7, n))
        end do
        !symmetrize
        virial(1, 2) = 0.5d0*(virial(1, 2) + virial(2, 1))
        virial(1, 3) = 0.5d0*(virial(1, 3) + virial(3, 1))
        virial(2, 3) = 0.5d0*(virial(2, 3) + virial(3, 2))
        virial(2, 1) = virial(1, 2)
        virial(3, 1) = virial(1, 3)
        virial(3, 2) = virial(2, 3)

    end subroutine AM_RECIP_field_torque_virial
!---------------------------------------------------------
    subroutine AM_RECIP_Fdip_to_Cdip_field( &
        numatoms, is_polarizable, nfft1, nfft2, nfft3, &
        recip, frac_dipole_field, cart_dipole_field)
        integer, intent(in) :: numatoms, nfft1, nfft2, nfft3
        logical, intent(in) :: is_polarizable(*)
        _REAL_, intent(in) :: recip(3, 3), frac_dipole_field(3, *)
        _REAL_, intent(out) :: cart_dipole_field(3, *)

        integer :: n, i, j
        _REAL_ :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)

        call AM_RECIP_xform_matrices(nfft1, nfft2, nfft3, recip, &
            mpole_xform_3x3, field_xform_3x3)
        do n = 1, numatoms
            if (is_polarizable(n)) then
                do i = 1, 3
                    ! add recip contribution to cart_dipole_field
                    do j = 1, 3
                        cart_dipole_field(i, n) = cart_dipole_field(i, n) + &
                            field_xform_3x3(i, j)*frac_dipole_field(j, n)
                    end do
                end do
            end if
        end do
    end subroutine AM_RECIP_Fdip_to_Cdip_field
!--------------------------------------------------------------
    subroutine AM_RECIP_dip_to_frac_dip( &
        numatoms, is_polarizable, nfft1, nfft2, nfft3, &
        recip, ind_dip1, ind_dip2, f_dipole1, f_dipole2)
        integer, intent(in) :: numatoms, nfft1, nfft2, nfft3
        logical, intent(in) :: is_polarizable(*)
        _REAL_, intent(in)  :: recip(3, 3), ind_dip1(3, *), ind_dip2(3, *)
        _REAL_, intent(out) :: f_dipole1(3, *), f_dipole2(3, *)

        integer :: n, i, j
        _REAL_ :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)

        call AM_RECIP_xform_matrices(nfft1, nfft2, nfft3, recip, &
            mpole_xform_3x3, field_xform_3x3)
        do n = 1, numatoms
            if (is_polarizable(n)) then
                do i = 1, 3
                    f_dipole1(i, n) = 0.d0
                    f_dipole2(i, n) = 0.d0
                    do j = 1, 3
                        f_dipole1(i, n) = f_dipole1(i, n) + mpole_xform_3x3(i, j)*ind_dip1(j, n)
                        f_dipole2(i, n) = f_dipole2(i, n) + mpole_xform_3x3(i, j)*ind_dip2(j, n)
                    end do
                end do
            end if
        end do
    end subroutine AM_RECIP_dip_to_frac_dip
!-----------------------------------------------------------------------
    subroutine AM_RECIP_xform_matrices(nfft1, nfft2, nfft3, recip, &
        mpole_xform_3x3, field_xform_3x3)
        integer, intent(in) :: nfft1, nfft2, nfft3
        _REAL_, intent(in) :: recip(3, 3)
        _REAL_, intent(out) :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)

        _REAL_ du1_dx, du1_dy, du1_dz, du2_dx, du2_dy, du2_dz, du3_dx, du3_dy, du3_dz

        du1_dx = nfft1*recip(1, 1)
        du1_dy = nfft1*recip(2, 1)
        du1_dz = nfft1*recip(3, 1)
        du2_dx = nfft2*recip(1, 2)
        du2_dy = nfft2*recip(2, 2)
        du2_dz = nfft2*recip(3, 2)
        du3_dx = nfft3*recip(1, 3)
        du3_dy = nfft3*recip(2, 3)
        du3_dz = nfft3*recip(3, 3)

        field_xform_3x3(1, 1) = du1_dx
        field_xform_3x3(1, 2) = du2_dx
        field_xform_3x3(1, 3) = du3_dx
        field_xform_3x3(2, 1) = du1_dy
        field_xform_3x3(2, 2) = du2_dy
        field_xform_3x3(2, 3) = du3_dy
        field_xform_3x3(3, 1) = du1_dz
        field_xform_3x3(3, 2) = du2_dz
        field_xform_3x3(3, 3) = du3_dz

        mpole_xform_3x3(1, 1) = du1_dx
        mpole_xform_3x3(1, 2) = du1_dy
        mpole_xform_3x3(1, 3) = du1_dz
        mpole_xform_3x3(2, 1) = du2_dx
        mpole_xform_3x3(2, 2) = du2_dy
        mpole_xform_3x3(2, 3) = du2_dz
        mpole_xform_3x3(3, 1) = du3_dx
        mpole_xform_3x3(3, 2) = du3_dy
        mpole_xform_3x3(3, 3) = du3_dz
    end subroutine AM_RECIP_xform_matrices
!--------------------------------------------------------------
    subroutine AM_test_frac_cart_ene(numatoms, nfft1, nfft2, nfft3, recip, &
        F_mpole, F_field, C_mpole)
        use amoeba_multipoles, only : coulomb_const_kcal_per_mole
        integer, intent(in) :: numatoms, nfft1, nfft2, nfft3
        _REAL_, intent(in) :: recip(3, 3), F_field(20, *), F_mpole(10, *), C_mpole(10, *)

        _REAL_ ene_F, ene_C, C_field(10)
        _REAL_ :: mpole_xform_3x3(3, 3), field_xform_3x3(3, 3)
        integer k, n, dimxy, order
        _REAL_ Field_xy(MAXMP*MAXMP)

        ene_F = 0.d0
        ene_C = 0.d0
        call AM_RECIP_xform_matrices(nfft1, nfft2, nfft3, recip, &
            mpole_xform_3x3, field_xform_3x3)
        call XFORM_MPOLE_field_matrix(field_xform_3x3, Field_xy, order)
        do n = 1, numatoms
            call XFORM_FIELD(Field_xy, dimxy, F_field(:, n), C_field, order)
            do k = 1, 10
                ene_F = ene_F + F_mpole(k, n)*F_field(k, n)
                ene_C = ene_C + C_mpole(k, n)*C_field(k)
            end do
        end do
        ene_F = 0.5d0*coulomb_const_kcal_per_mole*ene_F
        ene_C = 0.5d0*coulomb_const_kcal_per_mole*ene_C
        write (6, '(a,2(1x,f14.4))') 'AM_test_frac_cart:  ene_F,ene_C, = ', ene_F, ene_C
    end subroutine AM_test_frac_cart_ene
!--------------------------------------------------------------
end module amoeba_recip
