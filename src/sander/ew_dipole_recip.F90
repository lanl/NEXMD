! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"

module ew_dipole_recip

    public do_pmesh_dipole_kspace

    private fill_dipole_grid, grad_sum_dipolerc, scalar_sum_dipolerc

contains

!     --- DO_PMESH_dipole_KSPACE ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine do_pmesh_dipole_kspace here]
    subroutine do_pmesh_dipole_kspace( &
        numatoms, crd, charge, &
        recip, volume, ewald_coeff, &
        eer, frc, virial, dipole, field, &
        prefac1, prefac2, prefac3, fftable, &
        frcx, frcy, frcz)
        use ew_bspline
        use stack
        implicit none
        character(kind=1, len=22):: routine = "do_pmesh_dipole_kspace"

#  include "flocntrl.h"
#  include "ew_pme_recip.h"
#  include "def_time.h"

#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
        include 'mpif.h'
        integer ierr
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif

        ! INPUT
        !       numatoms:  number of atoms
        !       x,y,z:   atomic coords
        !       crd   atomic coords
        !       charge  atomic charges
        !       recip: 3x3 array of reciprocal unit cell vectors (stored as columns)
        !       volume: the volume of the unit cell
        !       ewald_coeff:   ewald convergence parameter
        !       order: the order of Bspline interpolation. E.g. cubic is order 4
        !          fifth degree is order 6 etc. The order must be an even number
        !          and at least 4.
        !       nfft1,nfft2,nfft3: the dimensions of the charge grid array

        integer numatoms
        _REAL_ crd(3, numatoms), &
            charge(numatoms), recip(3, 3), volume, ewald_coeff
        _REAL_ dipole(3, numatoms), field(3, numatoms)
        integer nmine, nderiv
        ! OUTPUT
        !       eer:  ewald reciprocal or k-space  energy
        !       frc forces incremented by k-space sum
        !       virial:  virial due to k-space sum (valid for atomic scaling;
        !                rigid molecule virial needs a correction term not
        !                computed here

        _REAL_ eer, frc(3, numatoms), virial(3, 3)

        ! HEAP STORAGE:  These arrays need to be preserved throughout simulation

        _REAL_ prefac1(*), prefac2(*), prefac3(*), fftable(*)

        ! STACK STORAGE: These arrays can be tossed after leaving this routine

        _REAL_ frcx, frcy, frcz

        integer nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw
        integer l_fr1, l_fr2, l_fr3, num_ks
        integer l_th1, l_th2, l_th3, l_dth1, l_dth2, l_dth3
        integer l_d2th1, l_d2th2, l_d2th3
        integer l_fftw, l_q
        integer imy_cg
        integer l_tmpy, l_alpha, l_beta, l_rctable, l_qq

        !     FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)
        if (do_rec == 0) return

        !  get some integer array dimensions
        call get_fftdims(nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork, sfft, sffw)
        nderiv = 2

# ifdef MPI
        num_ks = min((((nxyslab(0) + order - 1)*numatoms*4)/ &
            (3*nfft3)), numatoms)
# else
        num_ks = numatoms
# endif
        call get_stack(l_fftw, sizffwrk, routine)
        call get_stack(l_q, siz_q, routine)
        call get_stack(l_fr1, num_ks, routine)
        call get_stack(l_fr2, num_ks, routine)
        call get_stack(l_fr3, num_ks, routine)
        call get_stack(l_th1, num_ks*order, routine)
        call get_stack(l_th2, num_ks*order, routine)
        call get_stack(l_th3, num_ks*order, routine)
        if (nderiv >= 1) then
            call get_stack(l_dth1, num_ks*order, routine)
            call get_stack(l_dth2, num_ks*order, routine)
            call get_stack(l_dth3, num_ks*order, routine)
        end if
        if (nderiv == 2) then
            call get_stack(l_d2th1, num_ks*order, routine)
            call get_stack(l_d2th2, num_ks*order, routine)
            call get_stack(l_d2th3, num_ks*order, routine)
        end if

        call get_istack(imy_cg, num_ks, routine)
        if (.not. rstack_ok) then
            deallocate (r_stack)
            allocate (r_stack(1:lastrst), stat=alloc_ier)
            call reassign_rstack(routine)
        end if
        if (.not. istack_ok) then
            deallocate (i_stack)
            allocate (i_stack(1:lastist), stat=alloc_ier)
            call reassign_istack(routine)
        end if
        REQUIRE(rstack_ok)
        REQUIRE(istack_ok)

        call timer_start(TIME_BSPL)
        call get_grid_weights( &
            numatoms, crd, recip, nfft1, nfft2, nfft3, &
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
            call mexit(6, 1)
        end if
        call timer_stop_start(TIME_BSPL, TIME_FILLG)

        !........Fill Charge Grid

        call fill_dipole_grid(numatoms, charge, dipole, recip, &
            r_stack(l_dth1), r_stack(l_dth2), r_stack(l_dth3), &
            r_stack(l_th1), r_stack(l_th2), r_stack(l_th3), &
            r_stack(l_fr1), r_stack(l_fr2), r_stack(l_fr3), &
            order, nfft1, nfft2, nfft3, nfftdim1, nfftdim2 &
# ifndef MPI
            , nfftdim3 &
# else
            , mxyslabs &
# endif
            , r_stack(l_q), i_stack(imy_cg), nmine)
        call timer_stop_start(TIME_FILLG, TIME_FFT)

#ifdef MPI
        call mpi_barrier(recip_comm, ierr)
#endif
        call get_stack(l_tmpy, 2*nfftdim1, routine)
        call get_stack(l_alpha, nfft1, routine)
        call get_stack(l_beta, nfft1, routine)
        if (.not. rstack_ok) then
            deallocate (r_stack)
            allocate (r_stack(1:lastrst), stat=alloc_ier)
            call reassign_rstack(routine)
        end if
        call fft_backrc( &
            r_stack(l_q), fftable, r_stack(l_fftw), &
            nfft1, nfft2, nfft3, &
            nfftdim1, nfftdim2, nfftdim3, nfftable, nffwork &
            , r_stack(l_tmpy), &
            r_stack(l_alpha), r_stack(l_beta) &
            )

        call timer_stop_start(TIME_FFT, TIME_SCSUM)
        call scalar_sum_dipolerc( &
            r_stack(l_q), &
            ewald_coeff, volume, recip, prefac1, prefac2, prefac3, &
            nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, eer, virial)
# ifdef MPI
        call mpi_barrier(recip_comm, ierr)
# endif /* MPI */

        call timer_stop_start(TIME_SCSUM, TIME_FFT)
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
        call grad_sum_dipolerc( &
            numatoms, charge, dipole, recip, &
            r_stack(l_th1), r_stack(l_th2), r_stack(l_th3), &
            r_stack(l_dth1), r_stack(l_dth2), r_stack(l_dth3), &
            r_stack(l_d2th1), r_stack(l_d2th2), r_stack(l_d2th3), &
            frc, field, r_stack(l_fr1), r_stack(l_fr2), r_stack(l_fr3), &
            frcx, frcy, frcz, &
            order, nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, &
            r_stack(l_q), i_stack(imy_cg), nmine)
        call timer_stop(TIME_GRADS)
        if (nderiv == 2) then
            call free_stack(l_d2th3, routine)
            call free_stack(l_d2th2, routine)
            call free_stack(l_d2th1, routine)
        end if
        if (nderiv >= 1) then
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
    end subroutine do_pmesh_dipole_kspace

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_dipole_grid here]
    subroutine fill_dipole_grid( &
        numatoms, charge, dipole, recip, &
        dtheta1, dtheta2, dtheta3, &
        theta1, theta2, theta3, fr1, fr2, fr3, &
        order, nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, &
        q, my_cg, nmine)

        !---------------------------------------------------------------------
        ! INPUT:
        !      numatoms:  number of atoms
        !      charge: the array of atomic charges
        !      dipole: the array of atomic dipoles
        !      theta1,theta2,theta3: the spline coeff arrays
        !      dtheta1,dtheta2,dtheta3: derivs of the spline coeff arrays
        !      fr1,fr2,fr3 the scaled and shifted fractional coords
        !      nfft1,nfft2,nfft3: the charge grid dimensions
        !      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims
        !      order: the order of spline interpolation
        ! OUTPUT:
        !      Q the charge grid
        !---------------------------------------------------------------------
        use ew_bspline, only : kbot, ktop
        implicit none
        integer numatoms, order, nfft1, nfft2, nfft3
        integer nfftdim1, nfftdim2, nfftdim3
        _REAL_ fr1(numatoms), fr2(numatoms), fr3(numatoms)
        _REAL_ theta1(order, numatoms), theta2(order, numatoms), &
            theta3(order, numatoms), charge(numatoms)
        _REAL_ dtheta1(order, numatoms), dtheta2(order, numatoms), &
            dtheta3(order, numatoms), dipole(3, numatoms)
        !      _REAL_ Q(2,nfftdim1,nfftdim2,nfftdim3)
        _REAL_ q(2*nfftdim1, nfftdim2, nfftdim3)
        _REAL_ recip(3, 3)
        integer my_cg(*), nmine

#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
#endif

        integer im, n, ith1, ith2, ith3, i0, j0, k0, i, j, k, kq
        integer ntot, kbot0
        _REAL_ prod, prod1, prod2, prod3, mu1, mu2, mu3

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

            !       ---note mu1 is grad of u1 dotted with dipole
            !       where grad u1 is given by nfft1 times 1st recip vector etc.

            mu1 = nfft1*(recip(1, 1)*dipole(1, n) + recip(2, 1)*dipole(2, n) + &
                recip(3, 1)*dipole(3, n))
            mu2 = nfft2*(recip(1, 2)*dipole(1, n) + recip(2, 2)*dipole(2, n) + &
                recip(3, 2)*dipole(3, n))
            mu3 = nfft3*(recip(1, 3)*dipole(1, n) + recip(2, 3)*dipole(2, n) + &
                recip(3, 3)*dipole(3, n))
            k0 = int(fr3(im)) - order
            do ith3 = 1, order
                k0 = k0 + 1
                k = k0 + 1 + (nfft3 - sign(nfft3, k0))/2
#ifdef MPI
                if (k >= kbot .and. k <= ktop) then
#endif
#ifdef MPI
                    kq = k - kbot0
#else
                    kq = k
#endif
                    j0 = int(fr2(im)) - order
                    do ith2 = 1, order
                        j0 = j0 + 1
                        j = j0 + 1 + (nfft2 - sign(nfft2, j0))/2

                        !           ---field is negative of grad of phi
                        !           grad phi dotted with dipole re-expressed using chain rule. See above
                        !           comments about mu1,mu2,mu3

                        prod1 = theta2(ith2, im)*theta3(ith3, im)*charge(n) + &
                            dtheta2(ith2, im)*theta3(ith3, im)*mu2 + &
                            theta2(ith2, im)*dtheta3(ith3, im)*mu3
                        prod2 = theta2(ith2, im)*theta3(ith3, im)*mu1
                        i0 = int(fr1(im)) - order
                        do ith1 = 1, order
                            i0 = i0 + 1
                            i = i0 + 1 + (nfft1 - sign(nfft1, i0))/2
                            q(i, j, kq) = q(i, j, kq) + theta1(ith1, im)*prod1 &
                                + dtheta1(ith1, im)*prod2
                        end do
                    end do
#ifdef MPI
                end if
#endif
            end do
        end do  !  im = 1,nmine
        return
    end subroutine fill_dipole_grid
!-------------------------------------------------------------------
!     --- GRAD_SUM    RC---
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine grad_sum_dipolerc here]
    subroutine grad_sum_dipolerc( &
        numatoms, charge, dipole, recip, &
        theta1, theta2, theta3, &
        dtheta1, dtheta2, dtheta3, &
        d2theta1, d2theta2, d2theta3, &
        frc, field, fr1, fr2, fr3, &
        frcx, frcy, frcz, &
        order, nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, &
        q, my_cg, nmine)
        use ew_bspline, only : kbot, ktop
        implicit none
        integer numatoms, order, nfft1, nfft2, nfft3
        integer nfftdim1, nfftdim2, nfftdim3
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
        include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
        _REAL_ q(2*nfftdim1, nfftdim2, mxyslabs)
#else
        _REAL_ q(2*nfftdim1, nfftdim2, nfftdim3)

        ! If MPI only some atoms are touched by the local process. Hence loop is over
        ! igoo=1,igood, arrays are filled in order of being seen by process
        ! and n is obtained by a pointer.
        ! Due to similarity,otherwise of code, MPI is merged by simply re-defining
        ! igoo, igood and the array dimension kq.

#endif

        _REAL_ recip(3, 3)
        _REAL_ fr1(numatoms), fr2(numatoms), fr3(numatoms)
        _REAL_ frc(3, numatoms), field(3, numatoms)
        _REAL_ theta1(order, numatoms), theta2(order, numatoms), &
            theta3(order, numatoms), charge(numatoms), &
            dipole(3, numatoms)
        _REAL_ dtheta1(order, numatoms), dtheta2(order, numatoms), &
            dtheta3(order, numatoms)
        _REAL_ d2theta1(order, numatoms), &
            d2theta2(order, numatoms), d2theta3(order, numatoms)
        integer my_cg(*), nmine

        integer im, kq
        integer n, ith1, ith2, ith3, i0, j0, k0, i, j, k
        _REAL_ frcx, frcy, frcz
        _REAL_ dfx, dfy, dfz
        _REAL_ du1_dx, du1_dy, du1_dz, du2_dx, du2_dy, &
            du2_dz, du3_dx, du3_dy, du3_dz
        _REAL_ mu1, mu2, mu3
        _REAL_ dphi_du1fi, dphi_du2fi, dphi_du3fi
        _REAL_ phi, dphi_du1, dphi_du2, dphi_du3, d2phi_du1du1, &
            d2phi_du1du2, d2phi_du1du3, d2phi_du2du2, &
            d2phi_du2du3, d2phi_du3du3, &
            d2phi_du2du1, d2phi_du3du1, d2phi_du3du2
        _REAL_ ene, de_du1, de_du2, de_du3, &
            ft1, ft2, ft3, s, s1, s2, s3, ft22, ft23, ft33, &
            s11, s21, s31, s12, s22, s32, s13, s23, s33

        !     ---columns of recip are reciprocal unit cell vectors

        ene = 0.d0

        !     ---columns of recip are reciprocal unit cell vecs

        du1_dx = nfft1*recip(1, 1)
        du1_dy = nfft1*recip(2, 1)
        du1_dz = nfft1*recip(3, 1)
        du2_dx = nfft2*recip(1, 2)
        du2_dy = nfft2*recip(2, 2)
        du2_dz = nfft2*recip(3, 2)
        du3_dx = nfft3*recip(1, 3)
        du3_dy = nfft3*recip(2, 3)
        du3_dz = nfft3*recip(3, 3)
        do im = 1, nmine
#ifdef MPI
            n = my_cg(im)
#else
            n = im
#endif
            phi = 0.d0
            dphi_du1 = 0.d0
            dphi_du2 = 0.d0
            dphi_du3 = 0.d0
            dphi_du1fi = 0.d0
            dphi_du2fi = 0.d0
            dphi_du3fi = 0.d0
            d2phi_du1du1 = 0.d0
            d2phi_du1du2 = 0.d0
            d2phi_du1du3 = 0.d0
            d2phi_du2du1 = 0.d0
            d2phi_du2du2 = 0.d0
            d2phi_du2du3 = 0.d0
            d2phi_du3du1 = 0.d0
            d2phi_du3du2 = 0.d0
            d2phi_du3du3 = 0.d0

            !      ---note mu1 is grad of u1 dotted with dipole

            mu1 = nfft1*(recip(1, 1)*dipole(1, n) + recip(2, 1)*dipole(2, n) + &
                recip(3, 1)*dipole(3, n))
            mu2 = nfft2*(recip(1, 2)*dipole(1, n) + recip(2, 2)*dipole(2, n) + &
                recip(3, 2)*dipole(3, n))
            mu3 = nfft3*(recip(1, 3)*dipole(1, n) + recip(2, 3)*dipole(2, n) + &
                recip(3, 3)*dipole(3, n))
            k0 = int(fr3(im)) - order
            do ith3 = 1, order
                k0 = k0 + 1
                k = k0 + 1 + (nfft3 - sign(nfft3, k0))*.5
#ifdef MPI
                if (k >= kbot .and. k <= ktop) then
                    kq = k - mxystart(mytaskid)
#else
                kq = k
#endif
                j0 = int(fr2(im)) - order
                do ith2 = 1, order
                    j0 = j0 + 1
                    j = j0 + 1 + (nfft2 - sign(nfft2, j0))*.5
                    i0 = int(fr1(im)) - order
                    ft1 = theta2(ith2, im)*theta3(ith3, im)
                    ft2 = dtheta2(ith2, im)*theta3(ith3, im)
                    ft3 = theta2(ith2, im)*dtheta3(ith3, im)
                    ft22 = d2theta2(ith2, im)*theta3(ith3, im)
                    ft23 = dtheta2(ith2, im)*dtheta3(ith3, im)
                    ft33 = theta2(ith2, im)*d2theta3(ith3, im)
                    s1 = 0.d0
                    s2 = 0.d0
                    s3 = 0.d0
                    do ith1 = 1, order
                        i0 = i0 + 1
                        i = i0 + 1 + (nfft1 - sign(nfft1, i0))*.5
                        s1 = s1 + theta1(ith1, im)*q(i, j, kq)
                        s2 = s2 + dtheta1(ith1, im)*q(i, j, kq)
                        s3 = s3 + d2theta1(ith1, im)*q(i, j, kq)
                    end do
                    phi = phi + s1*ft1
                    dphi_du1 = dphi_du1 + s2*ft1
                    dphi_du2 = dphi_du2 + s1*ft2
                    dphi_du3 = dphi_du3 + s1*ft3
                    d2phi_du1du1 = d2phi_du1du1 + s3*ft1
                    d2phi_du1du2 = d2phi_du1du2 + s2*ft2
                    d2phi_du1du3 = d2phi_du1du3 + s2*ft3
                    d2phi_du2du2 = d2phi_du2du2 + s1*ft22
                    d2phi_du2du3 = d2phi_du2du3 + s1*ft23
                    d2phi_du3du3 = d2phi_du3du3 + s1*ft33
                end do
#ifdef MPI
            end if
#endif
        end do  !  ith3 = 1,order

        !        ---field is negative of grad of phi
        !        ene = ene + charge(n)*phi - dipole .dot. Efield
        !        re-express dot product using chain rule on grad of phi

        ene = ene + charge(n)*phi + mu1*dphi_du1 + mu2*dphi_du2 + &
            mu3*dphi_du3
        de_du1 = charge(n)*dphi_du1 + mu1*d2phi_du1du1 + &
            mu2*d2phi_du1du2 + mu3*d2phi_du1du3
        de_du2 = charge(n)*dphi_du2 + mu1*d2phi_du1du2 + &
            mu2*d2phi_du2du2 + mu3*d2phi_du2du3
        de_du3 = charge(n)*dphi_du3 + mu1*d2phi_du1du3 + &
            mu2*d2phi_du2du3 + mu3*d2phi_du3du3

        !      ---field is negative of grad of phi

        field(1, n) = field(1, n) - (dphi_du1*du1_dx + dphi_du2*du2_dx + &
            dphi_du3*du3_dx)
        field(2, n) = field(2, n) - (dphi_du1*du1_dy + dphi_du2*du2_dy + &
            dphi_du3*du3_dy)
        field(3, n) = field(3, n) - (dphi_du1*du1_dz + dphi_du2*du2_dz + &
            dphi_du3*du3_dz)

        !      ---force is negative of grad of ene

        dfx = de_du1*du1_dx + de_du2*du2_dx + de_du3*du3_dx
        dfy = de_du1*du1_dy + de_du2*du2_dy + de_du3*du3_dy
        dfz = de_du1*du1_dz + de_du2*du2_dz + de_du3*du3_dz
        frc(1, n) = frc(1, n) - dfx
        frc(2, n) = frc(2, n) - dfy
        frc(3, n) = frc(3, n) - dfz
        frcx = frcx - dfx
        frcy = frcy - dfy
        frcz = frcz - dfz
    end do  !  im = 1,nmine
    ene = 0.5d0*ene
    return
end subroutine grad_sum_dipolerc
!-------------------------------------------------------------------
!     --- SCALAR_SUM ---
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine scalar_sum_dipolerc here]
subroutine scalar_sum_dipolerc( &
    q, ewaldcof, volume, recip, prefac1, prefac2, prefac3, &
    nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3, eer, vir)
    use constants, only : pi, pi2
    implicit none
#  include "extra.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
    integer ier
#endif
    integer nfft1, nfft2, nfft3, nfftdim1, nfftdim2, nfftdim3
    _REAL_ prefac1(nfft1), prefac2(nfft2), &
        prefac3(nfft3), ewaldcof, volume
    _REAL_ eer, vir(3, 3)
    _REAL_ recip(3, 3)
    _REAL_ q(2, nfft3, nfftdim1, nfft2)
    _REAL_ fac, denom, eterm, vterm, energy
    integer k, k1, k2, k3, m1, m2, m3, nff, ind, jnd, indtop, k2q
    integer nf1, nf2, nf3, k10
    _REAL_ mhat1, mhat2, mhat3, msq, struc2, t1, t2
    integer k1s, k2s, k3s, m1s, m2s, m3s
    _REAL_ mhat1s, mhat2s, mhat3s, msqs, qr, qi
    _REAL_ struc2s, eterms, vterms, denoms

    indtop = nfft1*nfft2*nfft3
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
    do m2 = 1, 3
        do m1 = 1, 3
            vir(m1, m2) = 0.d0
        end do
    end do

    !........Insist that Q(1,1,1,1) is set to 0 (true already for neutral)
    !      Note that for parallel, only the first processor has the
    !      Q(1,1,1,1) so only that pe should set it to zero.
    !      For all other pe's, Q(1,1,1,1) is really in a slab in the middle
    !      of the grid

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
        do k3 = 1, nfft3
            k3s = mod(nfft3 - k3 + 1, nfft3) + 1
            k10 = 1
            if (master) then
                if (k3 + k2 == 2) k10 = 2
            end if
            do k1 = k10, nf1 + 1
                k1s = nfft1 - k1 + 2
                m1 = k1 - 1
                if (k1 > nf1) m1 = k1 - 1 - nfft1
                m2 = k2 - 1
                if (k2 > nf2) m2 = k2 - 1 - nfft2
                m3 = k3 - 1
                if (k3 > nf3) m3 = k3 - 1 - nfft3
                mhat1 = recip(1, 1)*m1 + recip(1, 2)*m2 + recip(1, 3)*m3
                mhat2 = recip(2, 1)*m1 + recip(2, 2)*m2 + recip(2, 3)*m3
                mhat3 = recip(3, 1)*m1 + recip(3, 2)*m2 + recip(3, 3)*m3
                msq = mhat1*mhat1 + mhat2*mhat2 + mhat3*mhat3
                denom = pi*volume*msq
                eterm = exp(-fac*msq)*prefac1(k1)*prefac2(k2)*prefac3(k3)/ &
                    denom
                vterm = 2.d0*(fac*msq + 1.d0)/msq
                !#ifdef MPI
                struc2 = q(1, k3, k1, k2q)*q(1, k3, k1, k2q) + &
                    q(2, k3, k1, k2q)*q(2, k3, k1, k2q)
                energy = energy + eterm*struc2
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
                    denoms = pi*volume*msqs
                    eterms = exp(-fac*msqs)*prefac1(k1s)*prefac2(k2s)*prefac3(k3s)/ &
                        denoms
                    vterms = 2.d0*(fac*msqs + 1.d0)/msqs
                    struc2s = q(1, k3, k1, k2q)*q(1, k3, k1, k2q) + &
                        q(2, k3, k1, k2q)*q(2, k3, k1, k2q)
                    energy = energy + eterms*struc2s
                    vir(1, 1) = vir(1, 1) + &
                        eterms*struc2s*(vterms*mhat1s*mhat1s - 1.d0)
                    vir(1, 2) = vir(1, 2) + eterms*struc2s*(vterms*mhat1s*mhat2s)
                    vir(1, 3) = vir(1, 3) + eterms*struc2s*(vterms*mhat1s*mhat3s)
                    vir(2, 1) = vir(2, 1) + eterms*struc2s*(vterms*mhat2s*mhat1s)
                    vir(2, 2) = vir(2, 2) + &
                        eterms*struc2s*(vterms*mhat2s*mhat2s - 1.d0)
                    vir(2, 3) = vir(2, 3) + eterms*struc2s*(vterms*mhat2s*mhat3s)
                    vir(3, 1) = vir(3, 1) + eterms*struc2s*(vterms*mhat3s*mhat1s)
                    vir(3, 2) = vir(3, 2) + eterms*struc2s*(vterms*mhat3s*mhat2s)
                    vir(3, 3) = vir(3, 3) + &
                        eterms*struc2s*(vterms*mhat3s*mhat3s - 1.d0)
                end if

                qr = q(1, k3, k1, k2q)
                qi = q(2, k3, k1, k2q)
                q(1, k3, k1, k2q) = eterm*q(1, k3, k1, k2q)
                q(2, k3, k1, k2q) = eterm*q(2, k3, k1, k2q)
            end do
        end do
    end do

    eer = 0.5d0*energy
    do m2 = 1, 3
        do m1 = 1, 3
            vir(m1, m2) = 0.5d0*vir(m1, m2)
        end do
    end do
    return
end subroutine scalar_sum_dipolerc

end module ew_dipole_recip
