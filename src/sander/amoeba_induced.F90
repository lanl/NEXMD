#include "dprec.fh"
#include "assert.fh"

!------------------------------------------------------------
module amoeba_induced
    implicit none
    private
    integer, save :: do_flag

    _REAL_, save, allocatable :: polarizability(:)
    _REAL_, save, allocatable :: sq_polinv(:)
    logical, save, allocatable :: is_polarizable(:)
    _REAL_, save, allocatable :: ind_dip_d(:, :)
    _REAL_, save, allocatable :: ind_dip_p(:, :)
    _REAL_, save, allocatable :: dip_field_d(:, :)
    _REAL_, save, allocatable :: dip_field_p(:, :)
    _REAL_, save, allocatable :: cart_dipole_field(:, :)
    _REAL_, save, allocatable :: dip_d_perm(:, :)
    _REAL_, save, allocatable :: dip_p_perm(:, :)
    _REAL_, save, allocatable :: old_dip_d(:, :)
    _REAL_, save, allocatable :: old_dip_p(:, :)

#ifndef DISABLE_AMOEBA_CG
    private :: cg_init, cg_fini, cg_eval
    _REAL_, save, allocatable :: cg_r_d(:, :), cg_r_p(:, :), cg_d_d(:, :), cg_d_p(:, :)
#endif /* DISABLE_AMOEBA_CG */

    public AM_INDUCED_eval, AM_INDUCED_set_user_bit, AM_INDUCED_readparm, &
        polarizability, sq_polinv, is_polarizable, ind_dip_d, ind_dip_p
#ifdef MPI
    public AM_INDUCED_bcast
#endif
contains
!------------------------------------------------------------------------
    subroutine AM_INDUCED_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"

        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
            call AM_INDUCED_set_induced_flags()
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_INDUCED_set_user_bit
!------------------------------------------------------------------------
#ifdef MPI
    subroutine AM_INDUCED_bcast(num_atoms)
        implicit none
        integer num_atoms, ier, n, j
        _REAL_ polmin

        include 'mpif.h'
#  include "extra.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ier)

        polmin = 1.d-12

        if (.not. master) then
            allocate (polarizability(num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (sq_polinv(num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (is_polarizable(num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (ind_dip_d(3, num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (ind_dip_p(3, num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (dip_field_d(3, num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (dip_field_p(3, num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (cart_dipole_field(3, num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (dip_d_perm(3, num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (dip_p_perm(3, num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (old_dip_d(3, num_atoms), stat=ier)
            REQUIRE(ier == 0)
            allocate (old_dip_p(3, num_atoms), stat=ier)
            REQUIRE(ier == 0)
        end if

        if (.not. master) then
            do n = 1, num_atoms
                is_polarizable(n) = .false.
                polarizability(n) = 0.d0
                sq_polinv(n) = 1.d0/sqrt(polmin)
                do j = 1, 3
                    ind_dip_d(j, n) = 0.d0
                    ind_dip_p(j, n) = 0.d0
                end do
            end do
        end if

        call mpi_bcast(polarizability, num_atoms, MPI_DOUBLE_PRECISION, 0, commsander, ier)

        if (.not. master) then
            do n = 1, num_atoms
                if (polarizability(n) > polmin) then
                    is_polarizable(n) = .true.
                    sq_polinv(n) = 1.d0/sqrt(polarizability(n))
                end if
            end do
        end if

    end subroutine AM_INDUCED_bcast
#endif

    function AM_INDUCED_readparm(nf, num_atoms, mass)
        integer :: AM_INDUCED_readparm
        integer, intent(in) :: nf, num_atoms

#include "do_flag.h"
        integer :: j, n, ier, dim1, numpolar
        _REAL_ :: polmin
        _REAL_ :: mass(*)

        AM_INDUCED_readparm = 0
        polmin = 1.d-12

        ! initialize things
        allocate (polarizability(num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (sq_polinv(num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (is_polarizable(num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (ind_dip_d(3, num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (ind_dip_p(3, num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (dip_field_d(3, num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (dip_field_p(3, num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (cart_dipole_field(3, num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (dip_d_perm(3, num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (dip_p_perm(3, num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (old_dip_d(3, num_atoms), stat=ier)
        REQUIRE(ier == 0)
        allocate (old_dip_p(3, num_atoms), stat=ier)
        REQUIRE(ier == 0)
        do n = 1, num_atoms
            is_polarizable(n) = .false.
            polarizability(n) = 0.d0
            sq_polinv(n) = 1.d0/sqrt(polmin)
            do j = 1, 3
                ind_dip_d(j, n) = 0.d0
                ind_dip_p(j, n) = 0.d0
            end do
        end do

        call AMOEBA_get_numlist('AMOEBA_POLARIZABILITY_', nf, numpolar)
        if (numpolar <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if

        if (numpolar /= num_atoms) then
            write (6, *) 'mismatch between polarizability and num_atoms'
            call mexit(6, 1)
        end if

        dim1 = 1

        call AMOEBA_read_real_list_data('AMOEBA_POLARIZABILITY_', nf, &
            dim1, numpolar, polarizability)

        do n = 1, numpolar
            if (polarizability(n) > polmin) then
                is_polarizable(n) = .true.
                sq_polinv(n) = 1.d0/sqrt(polarizability(n))
                if (abs(mass(n) - 40.078) < 0.0001d0) &
                    write (6, *) "Ca damp", n, mass(n), polarizability(n)
                if (abs(mass(n) - 40.078) < 0.0001d0) &
                    sq_polinv(n) = 1/1.35d0**3*sq_polinv(n)
                if (abs(mass(n) - 24.305) < 0.0001d0) &
                    write (6, *) "Mg damp", n, mass(n), polarizability(n)
                if (abs(mass(n) - 24.305) < 0.0001d0) &
                    sq_polinv(n) = 1/1.60d0**3*sq_polinv(n)
                if (abs(mass(n) - 65.409) < 0.0001d0) &
                    write (6, *) "Zn damp", n, mass(n), polarizability(n)
                if (abs(mass(n) - 65.409) < 0.0001d0) &
                    sq_polinv(n) = 1/1.34d0**3*sq_polinv(n)
            end if
        end do
        AM_INDUCED_readparm = 1
        do_flag = ibset(do_flag, VALID_BIT)

#ifndef DISABLE_AMOEBA_CG
        call cg_init(num_atoms)
#endif /* DISABLE_AMOEBA_CG */

    end function AM_INDUCED_readparm
!------------------------------------------------------------------------
    subroutine AM_INDUCED_set_induced_flags()
#  include "ew_cntrl.h"
#  include "ew_mpole.h"
        induced = 1
        indmeth = 1
    end subroutine AM_INDUCED_set_induced_flags
!------------------------------------------------------------------------
    subroutine AM_INDUCED_eval(numatoms, crd, x, ipairs, diprms, dipiter)
        use amoeba_mdin, only : sor_coefficient, dipole_scf_tol, &
#ifndef DISABLE_AMOEBA_CG
            dipole_scf_use_cg, &
#endif /* DISABLE_AMOEBA_CG */
            dipole_scf_iter_max, amoeba_verbose
        integer, intent(in) :: numatoms
        _REAL_, intent(in) ::  crd(3, *)
        _REAL_, intent(in) :: x(*)
        integer, intent(in) :: ipairs(*)
        _REAL_, intent(out) :: diprms, dipiter

        integer :: iter
        _REAL_ rms1, rms2, rms, oldrms
#include "do_flag.h"

        if (do_flag /= PROCEED) return

#ifndef DISABLE_AMOEBA_CG
        if (dipole_scf_use_cg .ne. 0) then
            call cg_eval(numatoms, crd, x, ipairs, diprms, dipiter)
            return
        end if ! dipole_scf_use_cg.ne.0
#endif /* DISABLE_AMOEBA_CG */

        call AM_NonBond_perm_fields(numatoms, is_polarizable, crd, x, ipairs, &
            dip_field_d, dip_field_p, cart_dipole_field)

        call AM_INDUCED_fields_to_ind_dips(numatoms, is_polarizable, &
            dip_field_d, dip_field_p, &
            polarizability, ind_dip_d, ind_dip_p)
        ! save the dips due to permanent fields
        call array_copy(ind_dip_d, dip_d_perm, 3*numatoms)
        call array_copy(ind_dip_p, dip_p_perm, 3*numatoms)
        iter = 0
        rms = 1.d0

        do
            ! copy dipoles to old
            call array_copy(ind_dip_d, old_dip_d, 3*numatoms)
            call array_copy(ind_dip_p, old_dip_p, 3*numatoms)
            call AM_NonBond_dip_dip_fields(numatoms, x, ind_dip_d, ind_dip_p, &
                dip_field_d, dip_field_p)
            ! get dipoles due to dipole fields
            call AM_INDUCED_fields_to_ind_dips(numatoms, is_polarizable, &
                dip_field_d, dip_field_p, &
                polarizability, ind_dip_d, ind_dip_p)
            ! add dipoles due to permfields
            call array_add(ind_dip_d, dip_d_perm, 3*numatoms)
            call array_add(ind_dip_p, dip_p_perm, 3*numatoms)
            call AM_INDUCED_SOR_update(numatoms, is_polarizable, sor_coefficient, &
                ind_dip_d, ind_dip_p, &
                old_dip_d, old_dip_p)
            call AM_INDUCED_rms_diff(numatoms, is_polarizable, rms1, rms2, &
                ind_dip_d, ind_dip_p, &
                old_dip_d, old_dip_p)
            oldrms = rms
            rms = MAX(rms1, rms2)
            iter = iter + 1
            if (rms < dipole_scf_tol) exit
            if (rms > oldrms) then
                write (6, *) 'dipoles failed to converge: rms increasing!'
                call mexit(6, 1)
            end if
            if (iter > dipole_scf_iter_max) then
                write (6, *) 'dipoles failed to converge: iter max exceeded!'
                call mexit(6, 1)
            end if
        end do
        if (amoeba_verbose > 0) then
            write (6, '(a,i5,1x,e12.4)') 'out of loop, iter,rms = ', iter, rms
        end if
        dipiter = iter
        diprms = rms

    end subroutine AM_INDUCED_eval
!----------------------------------------------------------

#ifndef DISABLE_AMOEBA_CG

    subroutine cg_init(numatoms)

        use amoeba_mdin, only : dipole_scf_use_cg, dipole_scf_cg_niter

        implicit none

        integer, intent(in) :: numatoms

#if defined(MPI) && defined(AMOEBA_WORKS_WITH_IT)
        include 'mpif.h'
#  include "extra.h"
#  include "parallel.h"
#endif /* MPI */

        integer :: ierr

#ifdef MPI
#endif /* MPI */

        if (dipole_scf_use_cg .eq. 0) &
            return

        write (unit=6, fmt='(/a)') &
            ' ** Info ** : Using Conjugate Gradients to compute induced dipoles for AMOEBA'

        allocate (cg_r_d(3, numatoms), cg_r_p(3, numatoms), &
            cg_d_d(3, numatoms), cg_d_p(3, numatoms), &
            stat=ierr)

        if (ierr /= 0) then
            write (unit=6, fmt='(a)') &
                ' ** Error ** : out of memory in amoeba_induced.f%cg_init()'
            call mexit(6, 1)
        end if

        do ierr = 1, numatoms
            cg_r_d(1:3, ierr) = 0.0d0
            cg_r_p(1:3, ierr) = 0.0d0

            cg_d_d(1:3, ierr) = 0.0d0
            cg_d_p(1:3, ierr) = 0.0d0
        end do

        dipole_scf_cg_niter = max(2, dipole_scf_cg_niter)

    end subroutine cg_init

    subroutine cg_fini()

        use amoeba_mdin, only : dipole_scf_use_cg

        implicit none

        if (dipole_scf_use_cg .ne. 0) then
            REQUIRE(allocated(cg_r_d))
            deallocate (cg_r_d, cg_r_p, cg_d_d, cg_d_p)
        end if ! dipole_scf_use_cg

    end subroutine cg_fini

! this SUBROUTINE implements naive conjugate gradients
! with Jacobi preconditioner as described, e.g., here:
! http://www.math-linux.com/spip.php?article55
!
! there are actually two systems of equations (d, p); both
! are treated as one twice bigger system to distribute the
! error more equally among d and p (this seems to lead to
! better energy conservation [and slows down convergence a bit])

    subroutine cg_eval(numatoms, crd, x, ipairs, diprms, dipiter)

        use amoeba_mdin, only : dipole_scf_cg_niter, amoeba_verbose

        implicit none

        integer, intent(in) :: numatoms
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(in) :: x(*)
        integer, intent(in) :: ipairs(*)

        _REAL_, intent(inout) :: diprms, dipiter

        _REAL_ :: alpha, beta, z_dot_r, z_dot_r_new, r_norm_d, r_norm_p, tmp

        integer :: n, j, iter

        if (amoeba_verbose /= 0) &
            write (unit=6, fmt='(a)') '>> entering amoeba_induced.f%cg_eval() >>'

        call AM_NonBond_perm_fields(numatoms, is_polarizable, crd, x, ipairs, &
            dip_field_d, dip_field_p, cart_dipole_field)

        ! get the initial guess (ind_dip_[dp] from previous step are better
        ! for "absolute" convergence, but worse for energy conservation)
        call AM_INDUCED_fields_to_ind_dips &
            (numatoms, is_polarizable, &
            dip_field_d, dip_field_p, &
            polarizability, &
            ind_dip_d, ind_dip_p)

        ! get initial residues
        call AM_NonBond_dip_dip_fields &
            (numatoms, x, &
            ind_dip_d, ind_dip_p, &
            dip_field_d, dip_field_p)

        ! now dip_field_[dp] contain fields due to ind_dip_[dp]
        z_dot_r = 0.0d0

        do n = 1, numatoms
            tmp = polarizability(n)
            do j = 1, 3
                cg_r_d(j, n) = -dip_field_d(j, n)
                cg_r_p(j, n) = -dip_field_p(j, n)

                ! initial search direction
                cg_d_d(j, n) = tmp*cg_r_d(j, n)
                cg_d_p(j, n) = tmp*cg_r_p(j, n)

                z_dot_r = z_dot_r &
                    + cg_d_d(j, n)*cg_r_d(j, n) &
                    + cg_d_p(j, n)*cg_r_p(j, n)
            end do
        end do

        do iter = 1, dipole_scf_cg_niter

            ! compute A d
            call AM_NonBond_dip_dip_fields &
                (numatoms, x, &
                cg_d_d, cg_d_p, &
                dip_field_d, dip_field_p)

            alpha = 0.0d0

            ! store A d in dip_field_[dp]
            do n = 1, numatoms
                if (is_polarizable(n)) then
                    tmp = 1.0d0/polarizability(n)
                    do j = 1, 3
                        dip_field_d(j, n) = dip_field_d(j, n) + tmp*cg_d_d(j, n)
                        alpha = alpha + cg_d_d(j, n)*dip_field_d(j, n)

                        dip_field_p(j, n) = dip_field_p(j, n) + tmp*cg_d_p(j, n)
                        alpha = alpha + cg_d_p(j, n)*dip_field_p(j, n)
                    end do
                end if ! is_polarizable(n)
            end do

            alpha = z_dot_r/alpha

            r_norm_d = 0.0d0
            r_norm_p = 0.0d0

            do n = 1, numatoms
                do j = 1, 3
                    ! update solution/residues
                    ind_dip_d(j, n) = ind_dip_d(j, n) + alpha*cg_d_d(j, n)
                    ind_dip_p(j, n) = ind_dip_p(j, n) + alpha*cg_d_p(j, n)

                    cg_r_d(j, n) = cg_r_d(j, n) - alpha*dip_field_d(j, n)
                    cg_r_p(j, n) = cg_r_p(j, n) - alpha*dip_field_p(j, n)

                    ! get the residue norms
                    r_norm_d = max(r_norm_d, polarizability(n)*abs(cg_r_d(j, n)))
                    r_norm_p = max(r_norm_p, polarizability(n)*abs(cg_r_p(j, n)))
                end do
            end do

            r_norm_d = 4.8033324d0*r_norm_d
            r_norm_p = 4.8033324d0*r_norm_p

            if (amoeba_verbose /= 0) &
                write (unit=6, fmt='(a,i3,a,f14.10,a,f14.10)') &
                'CG : iter = ', iter, ', residue norms: _d = ', &
                r_norm_d, ', _p = ', r_norm_p

            z_dot_r_new = 0.0d0

            do n = 1, numatoms
                z_dot_r_new = z_dot_r_new &
                    + polarizability(n)*(cg_r_d(1, n)**2 &
                    + cg_r_d(2, n)**2 &
                    + cg_r_d(3, n)**2 &
                    + cg_r_p(1, n)**2 &
                    + cg_r_p(2, n)**2 &
                    + cg_r_p(3, n)**2)
            end do

            beta = z_dot_r_new/z_dot_r
            z_dot_r = z_dot_r_new

            ! update search direction
            do n = 1, numatoms
                tmp = polarizability(n)
                do j = 1, 3
                    cg_d_d(j, n) = tmp*cg_r_d(j, n) + beta*cg_d_d(j, n)
                    cg_d_p(j, n) = tmp*cg_r_p(j, n) + beta*cg_d_p(j, n)
                end do
            end do
        end do ! iter = 1, dipole_scf_cg_niter

        diprms = -max(r_norm_d, r_norm_p)
        dipiter = dipole_scf_cg_niter

        if (amoeba_verbose /= 0) &
            write (unit=6, fmt='(a)') '<< leaving amoeba_induced.f%cg_eval() <<'

    end subroutine cg_eval

#endif /* DISABLE_AMOEBA_CG */

end module amoeba_induced
!----------------------------------------------------------
subroutine AM_INDUCED_add_cart_to_dip(numatoms, is_polarizable, &
    cart_dipole_field, &
    dip_field_d, dip_field_p)
    integer, intent(in) :: numatoms
    logical, intent(in) :: is_polarizable(*)
    _REAL_, intent(in) :: cart_dipole_field(3, *)
    _REAL_, intent(inout) :: dip_field_d(3, *), dip_field_p(3, *)

    integer :: i, n
    do n = 1, numatoms
        if (is_polarizable(n)) then
            do i = 1, 3
                dip_field_d(i, n) = dip_field_d(i, n) + cart_dipole_field(i, n)
                dip_field_p(i, n) = dip_field_p(i, n) + cart_dipole_field(i, n)
            end do
        end if
    end do
end subroutine AM_INDUCED_add_cart_to_dip
!-------------------------------------------------------------------------------
subroutine AM_INDUCED_fields_to_ind_dips(numatoms, is_polarizable, &
    dip_field_d, dip_field_p, polarizability, &
    ind_dip_d, ind_dip_p)
    integer, intent(in) :: numatoms
    logical, intent(in) :: is_polarizable(*)
    _REAL_, intent(in)  :: dip_field_d(3, *), dip_field_p(3, *), polarizability(*)
    _REAL_, intent(out) :: ind_dip_d(3, *), ind_dip_p(3, *)

    integer :: i, n
    _REAL_ pol
    do n = 1, numatoms
        if (is_polarizable(n)) then
            pol = polarizability(n)
            do i = 1, 3
                !minus sign since our field is actually gradient of potential
                ind_dip_d(i, n) = -pol*dip_field_d(i, n)
                ind_dip_p(i, n) = -pol*dip_field_p(i, n)
            end do
        else
            do i = 1, 3
                ind_dip_d(i, n) = 0.d0
                ind_dip_p(i, n) = 0.d0
            end do
        end if
    end do
end subroutine AM_INDUCED_fields_to_ind_dips
!-------------------------------------------------------------------------------
subroutine AM_INDUCED_SOR_update(numatoms, is_polarizable, sor_coeff, &
    ind_dip_d, ind_dip_p, old_dip_d, old_dip_p)
    integer, intent(in) :: numatoms
    logical, intent(in) :: is_polarizable(*)
    _REAL_, intent(in) :: sor_coeff
    _REAL_, intent(inout) :: ind_dip_d(3, *), ind_dip_p(3, *)
    _REAL_, intent(in) :: old_dip_d(3, *), old_dip_p(3, *)

    _REAL_ c_sor
    integer i, n

    c_sor = 1.d0 - sor_coeff
    do n = 1, numatoms
        if (is_polarizable(n)) then
            do i = 1, 3
                ind_dip_d(i, n) = c_sor*old_dip_d(i, n) + sor_coeff*ind_dip_d(i, n)
                ind_dip_p(i, n) = c_sor*old_dip_p(i, n) + sor_coeff*ind_dip_p(i, n)
            end do
        end if
    end do
end subroutine AM_INDUCED_SOR_update
!-------------------------------------------------------------------------------
subroutine AM_INDUCED_rms_diff(numatoms, is_polarizable, rms1, rms2, &
    ind_dip_d, ind_dip_p, old_dip_d, old_dip_p)
    integer, intent(in) :: numatoms
    logical, intent(in) :: is_polarizable(*)
    _REAL_, intent(out) :: rms1, rms2
    _REAL_, intent(in) :: ind_dip_d(3, *), ind_dip_p(3, *)
    _REAL_, intent(in) :: old_dip_d(3, *), old_dip_p(3, *)

    _REAL_ :: debye
    integer :: n, num

    debye = 4.8033324d0
    num = 0
    rms1 = 0.d0
    rms2 = 0.d0
    do n = 1, numatoms
        if (is_polarizable(n)) then
            num = num + 1
            rms1 = rms1 + (ind_dip_d(1, n) - old_dip_d(1, n))**2 + &
                (ind_dip_d(2, n) - old_dip_d(2, n))**2 + &
                (ind_dip_d(3, n) - old_dip_d(3, n))**2
            rms2 = rms2 + (ind_dip_p(1, n) - old_dip_p(1, n))**2 + &
                (ind_dip_p(2, n) - old_dip_p(2, n))**2 + &
                (ind_dip_p(3, n) - old_dip_p(3, n))**2
        end if
    end do
    if (num == 0) then
        write (6, *) 'AM_INDUCED_rms_diff: num polarizable = 0!!'
        call mexit(6, 1)
    end if
    rms1 = debye*sqrt(rms1/num)
    rms2 = debye*sqrt(rms2/num)
end subroutine AM_INDUCED_rms_diff
!-------------------------------------------------------------
