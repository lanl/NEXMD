#include "dprec.fh"
#include "assert.fh"

!---------------------------------------------------------
module amoeba_bonds
    implicit none
    private

    integer, save :: num_list = 0, startlist, endlist
    _REAL_, save :: energy, virial(3, 3)

    integer, allocatable, save :: list(:, :)

    _REAL_, allocatable, save :: force_constant(:), equil_value(:)
    _REAL_, allocatable, save :: fn(:), dfn_darg(:), arg(:), darg_dcrd(:, :)
    integer, save :: num_params = 0, do_flag = 0
    integer, save :: ftable_degree
    _REAL_, allocatable, save :: ftable_coeff(:)
    public AM_BONDS_readparm, AM_BONDS_deallocate, AM_BONDS_set_user_bit, &
        AM_BONDS_eval
#ifdef MPI
    public AM_BONDS_bcast
#endif
contains
!-------------------------------------------------------------
#ifdef MPI
    subroutine AM_BONDS_bcast
        implicit none
        integer siztask, ierr

        include 'mpif.h'
#  include "extra.h"
#  include "do_flag.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_list, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_params, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(ftable_degree, 1, MPI_INTEGER, 0, commsander, ierr)

        if (num_list .eq. 0) return

        if (.not. master) then
            allocate (list(3, num_list), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (force_constant(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (equil_value(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (ftable_coeff(0:ftable_degree), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        call mpi_bcast(list, 3*num_list, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(force_constant, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(equil_value, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(ftable_coeff, 1 + ftable_degree, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

        if (.not. master) then
            call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
            allocate (fn(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (arg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (dfn_darg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (darg_dcrd(3, siztask), stat=ierr)
            REQUIRE(ierr == 0)
        end if
    end subroutine AM_BONDS_bcast
#endif

    function AM_BONDS_readparm(nf)
        integer :: AM_BONDS_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, siztask, ier1, ier2, ier3, ier4

        AM_BONDS_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_REGULAR_BOND_', nf, num_list)
        if (num_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        !write(6,*)'num_list = ',num_list
        !allocate
        dim1 = 3 ! 2 atoms plus param ptr
        allocate (list(dim1, num_list), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate list; num_list = ', num_list
            call mexit(6, 1)
        end if
        call AMOEBA_read_list_data('AMOEBA_REGULAR_BOND_', nf, dim1, num_list, list)
        call AM_VAL_get_num_params('AMOEBA_REGULAR_BOND_', nf, num_params)
        if (num_params <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        !allocate
        allocate (force_constant(num_params), stat=ier1)
        allocate (equil_value(num_params), stat=ier2)
        if ((ier1 /= 0) .or. (ier2 /= 0)) then
            write (6, *) 'failed to allocate parms; num_params = ', num_params
            call mexit(6, 1)
        end if
        call AM_VAL_read_force_constant('AMOEBA_REGULAR_BOND_', nf, &
            num_params, force_constant)
        call AM_VAL_read_equil_value('AMOEBA_REGULAR_BOND_', nf, &
            num_params, equil_value)
        call AM_VAL_get_ftab_degree('AMOEBA_REGULAR_BOND_', nf, ftable_degree)
        allocate (ftable_coeff(0:ftable_degree), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate ftable_coeff; ftable_degree = ', &
                ftable_degree
            call mexit(6, 1)
        end if
        call AM_VAL_read_ftable_coeffs('AMOEBA_REGULAR_BOND_', nf, &
            ftable_degree, ftable_coeff)
        call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
        !allocate scratch
        allocate (fn(siztask), stat=ier1)
        allocate (arg(siztask), stat=ier2)
        allocate (dfn_darg(siztask), stat=ier3)
        allocate (darg_dcrd(3, siztask), stat=ier4)
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0)) then
            write (6, *) 'problems allocating bond scratch space '
            call mexit(6, 1)
        end if

        do_flag = ibset(do_flag, VALID_BIT)
        AM_BONDS_readparm = 1
    end function AM_BONDS_readparm
!-------------------------------------------------------------
!-------------------------------------------------------------
    subroutine AM_BONDS_deallocate()
        if (allocated(list)) deallocate (list)
        if (allocated(force_constant)) deallocate (force_constant)
        if (allocated(equil_value)) deallocate (equil_value)
        if (allocated(fn)) deallocate (fn)
        if (allocated(arg)) deallocate (arg)
        if (allocated(dfn_darg)) deallocate (dfn_darg)
        if (allocated(darg_dcrd)) deallocate (darg_dcrd)
        if (allocated(ftable_coeff)) deallocate (ftable_coeff)
    end subroutine AM_BONDS_deallocate
!-------------------------------------------------------------
    subroutine AM_BONDS_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_BONDS_set_user_bit
!-------------------------------------------------------------
    subroutine AM_BONDS_eval(crd, frc, ene, vir)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), vir(3, 3)
        _REAL_, intent(out) :: ene
#include "do_flag.h"

        integer j, k, nlist
! initialize
        energy = 0.d0
        do k = 1, 3
            do j = 1, 3
                virial(j, k) = 0.d0
            end do
        end do
        if (do_flag /= PROCEED) then
            return
        end if
        call AM_BONDS_get_args(crd, list, equil_value, arg, darg_dcrd)
        nlist = endlist - startlist + 1
        call AM_VAL_FTAB_eval_f_df(nlist, ftable_degree, ftable_coeff, &
            arg, fn, dfn_darg)
        call AM_BONDS_get_ene_frc(list, force_constant, fn, dfn_darg, darg_dcrd, crd, frc)
        ene = energy
        do k = 1, 3
            do j = 1, 3
                vir(j, k) = vir(j, k) + virial(j, k)
            end do
        end do
    end subroutine AM_BONDS_eval
!-------------------------------------------------------------
    subroutine AM_BONDS_get_args(crd, blist, equil_value, arg, darg_dcrd)
!--------------------------------------------------------
! This routine calculates bond function argument and its derivatives with
! respect to atomic positions of atom i. This gets used
! in a bond energy routine as well as restraint energy and
! any cross terms involving bond stretching. Note that derivatives with respect
! to atom j are jusr the negative of those with respect to atom i
! NOTE that periodic boundary conditions are not used here
!      i.e. imaging is done on a per molecule basis
!--------------------------------------------------------
! INPUT variables:
!    startlist,endlist  the set of bond indices to process--now global variates
!    crd the atomic coord array
!    blist: 3 x nbond array giving for each bond the index of the first
!           atom, index of the second atom, and index into the bond
!           parameter table giving the force constant and equilibrium length
!    equil_value the list of ref bond lengths
! OUTPUT variables:
!    arg, array of bond function args
!    darg_dcrdi   derivs of arg wrt crds of atom i
!--------------------------------------------------------

        _REAL_, intent(in) :: crd(3, *)
        integer, intent(in) :: blist(3, *)
        _REAL_, intent(in) :: equil_value(*)
        _REAL_, intent(out) :: arg(*), darg_dcrd(3, *)

        integer i, j, n, ind, it
        _REAL_ bl, dx, dy, dz

        do n = startlist, endlist
            ind = n - startlist + 1
            i = blist(1, n)
            j = blist(2, n)
            it = blist(3, n)
            dx = crd(1, i) - crd(1, j)
            dy = crd(2, i) - crd(2, j)
            dz = crd(3, i) - crd(3, j)
            bl = sqrt(dx*dx + dy*dy + dz*dz)
! recall the reference bond length is 2nd parameter
            arg(ind) = bl - equil_value(it)
! differentiate bl to get darg_dcrdi
            darg_dcrd(1, ind) = dx/bl
            darg_dcrd(2, ind) = dy/bl
            darg_dcrd(3, ind) = dz/bl
        end do

        return
    end subroutine AM_BONDS_get_args
!--------------------------------------------------------
    subroutine AM_BONDS_get_ene_frc(blist, force_constant, &
        fn, dfn_darg, darg_dcrd, crd, frc)
        integer, intent(in) :: blist(3, *)
        _REAL_, intent(in) :: force_constant(*), fn(*), dfn_darg(*), darg_dcrd(3, *)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *)

        integer i, j, n, p, q, ind, it
        _REAL_ f(6), frc_K

        do n = startlist, endlist
            ind = n - startlist + 1
            i = blist(1, n)
            j = blist(2, n)
            it = blist(3, n)
            frc_K = force_constant(it)
            energy = energy + frc_K*fn(ind)
! apply chain rule to get deriv of energy with respect to crds of i
!  df holds the deriv of f with respect to its arg (i.e. b-b0)
!  while darg_dcrdi holds the derivs of b-b0 with respect to crds of i
            f(1) = frc_K*dfn_darg(ind)*darg_dcrd(1, ind)
            f(2) = frc_K*dfn_darg(ind)*darg_dcrd(2, ind)
            f(3) = frc_K*dfn_darg(ind)*darg_dcrd(3, ind)
! deriv wrt j is opposite to that wrt i
            f(4) = -f(1)
            f(5) = -f(2)
            f(6) = -f(3)
! recall force is negative of grad
            frc(1, i) = frc(1, i) - f(1)
            frc(2, i) = frc(2, i) - f(2)
            frc(3, i) = frc(3, i) - f(3)
            frc(1, j) = frc(1, j) - f(4)
            frc(2, j) = frc(2, j) - f(5)
            frc(3, j) = frc(3, j) - f(6)
! update the virial
            do q = 1, 3
                do p = 1, 3
                    virial(p, q) = virial(p, q) + f(p)*crd(q, i) + f(p + 3)*crd(q, j)
                end do
            end do
        end do

        return
    end subroutine AM_BONDS_get_ene_frc
!----------------------------------------------------
end module amoeba_bonds
!-----------------------------------------------------
module amoeba_ureyb
    implicit none
    private

    integer, save :: num_list = 0, startlist, endlist
    _REAL_, save :: energy, virial(3, 3)

    integer, allocatable, save :: list(:, :)

    _REAL_, allocatable, save :: force_constant(:), equil_value(:)
    _REAL_, allocatable, save :: fn(:), dfn_darg(:), arg(:), darg_dcrd(:, :)
    integer, save :: num_params = 0, do_flag = 0
    integer, save :: ftable_degree
    _REAL_, allocatable, save :: ftable_coeff(:)
    public AM_UREYB_readparm, AM_UREYB_deallocate, AM_UREYB_set_user_bit, &
        AM_UREYB_eval
#ifdef MPI
    public AM_UREYB_bcast
#endif
contains
!-------------------------------------------------------------
#ifdef MPI
    subroutine AM_UREYB_bcast
        implicit none
        integer siztask, ierr

        include 'mpif.h'
#  include "extra.h"
#  include "do_flag.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_list, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_params, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(ftable_degree, 1, MPI_INTEGER, 0, commsander, ierr)

        if (num_list .eq. 0) return

        if (.not. master) then
            allocate (list(3, num_list), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (force_constant(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (equil_value(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (ftable_coeff(0:ftable_degree), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        call mpi_bcast(list, 3*num_list, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(force_constant, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(equil_value, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(ftable_coeff, 1 + ftable_degree, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

        if (.not. master) then
            call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
            !allocate scratch
            allocate (fn(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (arg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (dfn_darg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (darg_dcrd(3, siztask), stat=ierr)
            REQUIRE(ierr == 0)
        end if
    end subroutine AM_UREYB_bcast
#endif

    function AM_UREYB_readparm(nf)
        integer :: AM_UREYB_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, siztask, ier1, ier2, ier3, ier4

        AM_UREYB_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_UREY_BRADLEY_BOND_', nf, num_list)
        if (num_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        !write(6,*)'num_list = ',num_list
        !allocate
        dim1 = 3 ! 2 atoms plus param ptr
        allocate (list(dim1, num_list), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate list; num_list = ', num_list
            call mexit(6, 1)
        end if
        call AMOEBA_read_list_data('AMOEBA_UREY_BRADLEY_BOND_', nf, dim1, num_list, list)
        call AM_VAL_get_num_params('AMOEBA_UREY_BRADLEY_BOND_', nf, num_params)
        if (num_params <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        !allocate
        allocate (force_constant(num_params), stat=ier1)
        allocate (equil_value(num_params), stat=ier2)
        if ((ier1 /= 0) .or. (ier2 /= 0)) then
            write (6, *) 'failed to allocate parms; num_params = ', num_params
            call mexit(6, 1)
        end if
        call AM_VAL_read_force_constant('AMOEBA_UREY_BRADLEY_BOND_', nf, &
            num_params, force_constant)
        call AM_VAL_read_equil_value('AMOEBA_UREY_BRADLEY_BOND_', nf, &
            num_params, equil_value)
        call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
        !allocate
        allocate (fn(siztask), stat=ier1)
        allocate (arg(siztask), stat=ier2)
        allocate (dfn_darg(siztask), stat=ier3)
        allocate (darg_dcrd(3, siztask), stat=ier4)
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0)) then
            write (6, *) 'problems allocating bond scratch space '
            call mexit(6, 1)
        end if
        call AM_VAL_get_ftab_degree('AMOEBA_UREY_BRADLEY_BOND_', nf, ftable_degree)
        allocate (ftable_coeff(0:ftable_degree), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate ftable_coeff; ftable_degree = ', &
                ftable_degree
            call mexit(6, 1)
        end if
        call AM_VAL_read_ftable_coeffs('AMOEBA_UREY_BRADLEY_BOND_', nf, &
            ftable_degree, ftable_coeff)
        do_flag = ibset(do_flag, VALID_BIT)
        AM_UREYB_readparm = 1
    end function AM_UREYB_readparm
!-------------------------------------------------------------
!-------------------------------------------------------------
    subroutine AM_UREYB_deallocate()
        if (allocated(list)) deallocate (list)
        if (allocated(force_constant)) deallocate (force_constant)
        if (allocated(equil_value)) deallocate (equil_value)
        if (allocated(fn)) deallocate (fn)
        if (allocated(arg)) deallocate (arg)
        if (allocated(dfn_darg)) deallocate (dfn_darg)
        if (allocated(darg_dcrd)) deallocate (darg_dcrd)
        if (allocated(ftable_coeff)) deallocate (ftable_coeff)
    end subroutine AM_UREYB_deallocate
!-------------------------------------------------------------
    subroutine AM_UREYB_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_UREYB_set_user_bit
!-------------------------------------------------------------
    subroutine AM_UREYB_eval(crd, frc, ene, vir)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), vir(3, 3)
        _REAL_, intent(out) :: ene
#include "do_flag.h"

        integer j, k, nlist
! initialize
        energy = 0.d0
        do k = 1, 3
            do j = 1, 3
                virial(j, k) = 0.d0
            end do
        end do
        if (do_flag /= PROCEED) then
            return
        end if
        call AM_UREYB_get_args(crd, list, equil_value, arg, darg_dcrd)
        nlist = endlist - startlist + 1
        call AM_VAL_FTAB_eval_f_df(nlist, ftable_degree, ftable_coeff, &
            arg, fn, dfn_darg)
        call AM_UREYB_get_ene_frc(list, force_constant, fn, dfn_darg, darg_dcrd, crd, frc)
        ene = energy
        do k = 1, 3
            do j = 1, 3
                vir(j, k) = vir(j, k) + virial(j, k)
            end do
        end do
    end subroutine AM_UREYB_eval
!---------------------------------------------------------------
    subroutine AM_UREYB_get_args(crd, blist, equil_value, arg, darg_dcrdi)
!--------------------------------------------------------
! These urey bradley routines are essentially identical to the regular bond
! routines---
!--------------------------------------------------------
        _REAL_, intent(in) :: crd(3, *)
        integer, intent(in) :: blist(3, *)
        _REAL_, intent(in) :: equil_value(*)
        _REAL_, intent(out) :: arg(*), darg_dcrdi(3, *)

        integer :: i, j, n, ind, it
        _REAL_ :: bl, dx, dy, dz

        do n = startlist, endlist
            ind = n - startlist + 1
            i = blist(1, n)
            j = blist(2, n)
            it = blist(3, n)
            dx = crd(1, i) - crd(1, j)
            dy = crd(2, i) - crd(2, j)
            dz = crd(3, i) - crd(3, j)
            bl = sqrt(dx*dx + dy*dy + dz*dz)
! recall the reference bond length is 2nd parameter
            arg(ind) = bl - equil_value(it)
! differentiate bl to get darg_dcrdi
            darg_dcrdi(1, ind) = dx/bl
            darg_dcrdi(2, ind) = dy/bl
            darg_dcrdi(3, ind) = dz/bl
        end do

        return
    end subroutine AM_UREYB_get_args
!--------------------------------------------------------
    subroutine AM_UREYB_get_ene_frc(blist, force_constant, fn, dfn, darg_dcrdi, crd, frc)
        integer, intent(in) :: blist(3, *)
        _REAL_, intent(in) :: force_constant(*), fn(*), dfn(*), &
            darg_dcrdi(3, *), crd(3, *)
        _REAL_, intent(inout) :: frc(3, *)

        integer :: i, j, n, p, q, it, ind
        _REAL_ f(6), frc_K

        do n = startlist, endlist
            ind = n - startlist + 1
            i = blist(1, n)
            j = blist(2, n)
            it = blist(3, n)
            frc_K = force_constant(it)
            energy = energy + frc_K*fn(ind)
! apply chain rule to get deriv of energy with respect to crds of i
!  df holds the deriv of f with respect to its arg (i.e. b-b0)
!  while darg_dcrdi holds the derivs of b-b0 with respect to crds of i
            f(1) = frc_K*dfn(ind)*darg_dcrdi(1, ind)
            f(2) = frc_K*dfn(ind)*darg_dcrdi(2, ind)
            f(3) = frc_K*dfn(ind)*darg_dcrdi(3, ind)
! deriv wrt j is opposite to that wrt i
            f(4) = -f(1)
            f(5) = -f(2)
            f(6) = -f(3)
! recall force is negative of grad
            frc(1, i) = frc(1, i) - f(1)
            frc(2, i) = frc(2, i) - f(2)
            frc(3, i) = frc(3, i) - f(3)
            frc(1, j) = frc(1, j) - f(4)
            frc(2, j) = frc(2, j) - f(5)
            frc(3, j) = frc(3, j) - f(6)
! update the virial
            do q = 1, 3
                do p = 1, 3
                    virial(p, q) = virial(p, q) + f(p)*crd(q, i) + f(p + 3)*crd(q, j)
                end do
            end do
        end do

        return
    end subroutine AM_UREYB_get_ene_frc
!-----------------------------------------------------
end module amoeba_ureyb
!-----------------------------------------------------
module amoeba_reg_angles

    use constants, only : PI

    implicit none
    private

    integer, save :: num_list = 0, startlist, endlist
    _REAL_, save :: energy, virial(3, 3)

    integer, allocatable, save :: list(:, :)

    _REAL_, allocatable, save :: force_constant(:), equil_value(:)
    _REAL_, allocatable, save :: fn(:), dfn_darg(:), arg(:), darg_dcrd(:, :)
    integer, save :: num_params = 0, do_flag = 0
    integer, save :: ftable_degree
    _REAL_, allocatable, save :: ftable_coeff(:)
    _REAL_, parameter ::  pt999999 = 0.999999d0
    _REAL_, parameter :: radians_to_degrees = 180.d0/pi
    _REAL_, parameter :: degrees_to_radians = pi/180.d0
    public AM_REG_ANGLES_readparm, AM_REG_ANGLES_deallocate, &
        AM_REG_ANGLES_set_user_bit, AM_REG_ANGLES_eval
#ifdef MPI
    public AM_REG_ANGLES_bcast
#endif
contains
!-----------------------------------------------------

#ifdef MPI
    subroutine AM_REG_ANGLES_bcast
        implicit none
        integer siztask, ierr

        include 'mpif.h'
#  include "extra.h"
#  include "do_flag.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_list, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_params, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(ftable_degree, 1, MPI_INTEGER, 0, commsander, ierr)

        if (num_list .eq. 0) return

        if (.not. master) then
            allocate (list(4, num_list), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (force_constant(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (equil_value(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (ftable_coeff(0:ftable_degree), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        call mpi_bcast(list, 4*num_list, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(force_constant, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(equil_value, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(ftable_coeff, 1 + ftable_degree, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

        if (.not. master) then
            call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
            allocate (fn(siztask), stat=ierr)
            allocate (arg(siztask), stat=ierr)
            allocate (dfn_darg(siztask), stat=ierr)
            allocate (darg_dcrd(9, siztask), stat=ierr)
        end if
    end subroutine AM_REG_ANGLES_bcast
#endif

    function AM_REG_ANGLES_readparm(nf)
        integer :: AM_REG_ANGLES_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, siztask, ier1, ier2, ier3, ier4

        AM_REG_ANGLES_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_REGULAR_ANGLE_', nf, num_list)
        if (num_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        dim1 = 4 ! 3 atoms plus param ptr
        allocate (list(dim1, num_list), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate list; num_list = ', num_list
            call mexit(6, 1)
        end if
        call AMOEBA_read_list_data('AMOEBA_REGULAR_ANGLE_', nf, dim1, num_list, list)
        call AM_VAL_get_num_params('AMOEBA_REGULAR_ANGLE_', nf, num_params)
        if (num_params <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        allocate (force_constant(num_params), stat=ier1)
        allocate (equil_value(num_params), stat=ier2)
        if ((ier1 /= 0) .or. (ier2 /= 0)) then
            write (6, *) 'failed to allocate parms; num_params = ', num_params
            call mexit(6, 1)
        end if
        call AM_VAL_read_force_constant('AMOEBA_REGULAR_ANGLE_', nf, &
            num_params, force_constant)
        call AM_VAL_read_equil_value('AMOEBA_REGULAR_ANGLE_', nf, &
            num_params, equil_value)
        call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
        allocate (fn(siztask), stat=ier1)
        allocate (arg(siztask), stat=ier2)
        allocate (dfn_darg(siztask), stat=ier3)
        allocate (darg_dcrd(9, siztask), stat=ier4) !3 atoms to get derivs wrt
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0)) then
            write (6, *) 'problems allocating angle scratch space '
            call mexit(6, 1)
        end if
        call AM_VAL_get_ftab_degree('AMOEBA_REGULAR_ANGLE_', nf, ftable_degree)
        allocate (ftable_coeff(0:ftable_degree), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate ftable_coeff; ftable_degree = ', &
                ftable_degree
            call mexit(6, 1)
        end if
        call AM_VAL_read_ftable_coeffs('AMOEBA_REGULAR_ANGLE_', nf, &
            ftable_degree, ftable_coeff)
        do_flag = ibset(do_flag, VALID_BIT)
        AM_REG_ANGLES_readparm = 1
    end function AM_REG_ANGLES_readparm
!-------------------------------------------------------------
!-------------------------------------------------------------
    subroutine AM_REG_ANGLES_deallocate()
        if (allocated(list)) deallocate (list)
        if (allocated(force_constant)) deallocate (force_constant)
        if (allocated(equil_value)) deallocate (equil_value)
        if (allocated(fn)) deallocate (fn)
        if (allocated(arg)) deallocate (arg)
        if (allocated(dfn_darg)) deallocate (dfn_darg)
        if (allocated(darg_dcrd)) deallocate (darg_dcrd)
        if (allocated(ftable_coeff)) deallocate (ftable_coeff)
    end subroutine AM_REG_ANGLES_deallocate
!-------------------------------------------------------------
    subroutine AM_REG_ANGLES_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_REG_ANGLES_set_user_bit
!-------------------------------------------------------------
    subroutine AM_REG_ANGLES_eval(crd, frc, ene, vir)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), vir(3, 3)
        _REAL_, intent(out) :: ene
#include "do_flag.h"

        integer j, k, nlist
! initialize
        energy = 0.d0
        do k = 1, 3
            do j = 1, 3
                virial(j, k) = 0.d0
            end do
        end do
        if (do_flag /= PROCEED) then
            return
        end if
        call AM_REG_ANGLES_get_args(crd, list, equil_value, arg, darg_dcrd)
        nlist = endlist - startlist + 1
        call AM_VAL_FTAB_eval_f_df(nlist, ftable_degree, ftable_coeff, &
            arg, fn, dfn_darg)
        call AM_REG_ANGLES_get_ene_frc(list, force_constant, &
            fn, dfn_darg, darg_dcrd, crd, frc)
        ene = energy
        do k = 1, 3
            do j = 1, 3
                vir(j, k) = vir(j, k) + virial(j, k)
            end do
        end do
    end subroutine AM_REG_ANGLES_eval
!-------------------------------------------------------------
    subroutine AM_REG_ANGLES_get_args(crd, alist, equil_value, arg, darg_dcrdijk)
!--------------------------------------------------------
! This routine calculates angle function argument and its derivatives with
! respect to atomic positions of atom i,j,k.
!--------------------------------------------------------
! INPUT variables:
!    startlist,endlist  the set of bond indices to process--now global variates
!    crd the atomic coord array
!    alist: 4 x nangles array giving for each angle the index of the first
!           atom, index of the second atom, index of the thirs,
!           and the param table pointer
!    equil_value the list of angle ref angle
! OUTPUT variables:
!    arg, array of angle function args
!    darg_dcrdi   derivs of arg wrt crds of atom i
!    darg_dcrdj   derivs of arg wrt crds of atom j
!    darg_dcrdk   derivs of arg wrt crds of atom k
!--------------------------------------------------------
        integer, intent(in) :: alist(4, *)
        _REAL_, intent(in) :: crd(3, *), equil_value(*)
        _REAL_, intent(out) :: arg(*), darg_dcrdijk(9, *)
! LOCAL variables
        integer :: i, j, k, n, ind, it
        _REAL_ :: xij, yij, zij, xkj, ykj, zkj, cosang, &
            dotp, lenij2, lenkj2, lenp, ang, dang_dcosang, &
            dcosang_dxij, dcosang_dyij, dcosang_dzij, &
            dcosang_dxkj, dcosang_dykj, dcosang_dzkj, ang0
! note units are degrees not radians...possibly change back later
        do n = startlist, endlist
            ind = n - startlist + 1
            i = alist(1, n)
            j = alist(2, n)
            k = alist(3, n)
            it = alist(4, n)
            ang0 = equil_value(it)
            xij = crd(1, i) - crd(1, j)
            yij = crd(2, i) - crd(2, j)
            zij = crd(3, i) - crd(3, j)
            xkj = crd(1, k) - crd(1, j)
            ykj = crd(2, k) - crd(2, j)
            zkj = crd(3, k) - crd(3, j)
! cosine of angle is given by dot product of rij and rkj
!   divided by the product of the lengths of rij and rkj
            dotp = xij*xkj + yij*ykj + zij*zkj
            lenij2 = xij**2 + yij**2 + zij**2
            lenkj2 = xkj**2 + ykj**2 + zkj**2
            lenp = sqrt(lenij2*lenkj2)
            cosang = dotp/lenp
! avoid angle of pi and 0; really you need to use arcsin formulation
! near those points. however due to severe strain you could get bad values
! even for reference angle not near 0 or pi.
            cosang = max(-pt999999, cosang)
            cosang = min(pt999999, cosang)
            ang = radians_to_degrees*acos(cosang)
            dang_dcosang = -radians_to_degrees/sqrt(1.d0 - cosang**2)
! again note units are degrees for now
!       ang = acos(cosang)
!       dang_dcosang = -1.d0 / sqrt(1.d0-cosang**2)
! angle function argument is angle ijk minus reference angle
! reference angle is aparm(2,at)
            arg(ind) = ang - ang0
! deriv of dotp wrt xij is xkj; deriv of lenp^-1 wrt xij is
!  lenkj^-1 * (-lenij^-2) * (xij/lenij) = -xij / (lenp*lenij2)
! similar for others
            dcosang_dxij = xkj/lenp - (dotp*xij)/(lenp*lenij2)
            dcosang_dyij = ykj/lenp - (dotp*yij)/(lenp*lenij2)
            dcosang_dzij = zkj/lenp - (dotp*zij)/(lenp*lenij2)
            dcosang_dxkj = xij/lenp - (dotp*xkj)/(lenp*lenkj2)
            dcosang_dykj = yij/lenp - (dotp*ykj)/(lenp*lenkj2)
            dcosang_dzkj = zij/lenp - (dotp*zkj)/(lenp*lenkj2)
! now use the chain rule
! first the i crds
            darg_dcrdijk(1, ind) = dang_dcosang*dcosang_dxij
            darg_dcrdijk(2, ind) = dang_dcosang*dcosang_dyij
            darg_dcrdijk(3, ind) = dang_dcosang*dcosang_dzij
! next the k crds
            darg_dcrdijk(7, ind) = dang_dcosang*dcosang_dxkj
            darg_dcrdijk(8, ind) = dang_dcosang*dcosang_dykj
            darg_dcrdijk(9, ind) = dang_dcosang*dcosang_dzkj
! finally the j crds
            darg_dcrdijk(4, ind) = -dang_dcosang*(dcosang_dxij + dcosang_dxkj)
            darg_dcrdijk(5, ind) = -dang_dcosang*(dcosang_dyij + dcosang_dykj)
            darg_dcrdijk(6, ind) = -dang_dcosang*(dcosang_dzij + dcosang_dzkj)
        end do

        return
    end subroutine AM_REG_ANGLES_get_args
!----------------------------------------------------------
    subroutine AM_REG_ANGLES_get_ene_frc(alist, force_constant, func, dfunc, &
        darg_dcrdijk, crd, frc)
        integer, intent(in) :: alist(4, *) ! 3 atoms plus param ptr
        _REAL_, intent(in) :: force_constant(*), func(*), dfunc(*), &
            darg_dcrdijk(9, *), crd(3, *)
        _REAL_, intent(inout) :: frc(3, *)

        integer i, j, k, it, n, p, q, ind
        _REAL_ term, Frc_K, f(9)
        _REAL_ factor
        ! correct bending constant for units in degrees instead of radians
        factor = degrees_to_radians**2

        do n = startlist, endlist
            ind = n - startlist + 1
            i = alist(1, n)
            j = alist(2, n)
            k = alist(3, n)
            it = alist(4, n)
            Frc_K = factor*force_constant(it)
            energy = energy + Frc_K*func(ind)
! apply chain rule to get deriv of energy with respect to crds of i,j,k
! df holds the deriv of f with respect to its arg (i.e. a-a0)
! while darg_dcrdijk holds the derivs of arg with respect to crds of i,j,k
! recall force is negative of grad
            term = Frc_K*dfunc(ind)
            f(1) = term*darg_dcrdijk(1, ind)
            f(2) = term*darg_dcrdijk(2, ind)
            f(3) = term*darg_dcrdijk(3, ind)
            f(4) = term*darg_dcrdijk(4, ind)
            f(5) = term*darg_dcrdijk(5, ind)
            f(6) = term*darg_dcrdijk(6, ind)
            f(7) = term*darg_dcrdijk(7, ind)
            f(8) = term*darg_dcrdijk(8, ind)
            f(9) = term*darg_dcrdijk(9, ind)
            frc(1, i) = frc(1, i) - f(1)
            frc(2, i) = frc(2, i) - f(2)
            frc(3, i) = frc(3, i) - f(3)
            frc(1, j) = frc(1, j) - f(4)
            frc(2, j) = frc(2, j) - f(5)
            frc(3, j) = frc(3, j) - f(6)
            frc(1, k) = frc(1, k) - f(7)
            frc(2, k) = frc(2, k) - f(8)
            frc(3, k) = frc(3, k) - f(9)
! now get virial
            do q = 1, 3
                do p = 1, 3
                    virial(p, q) = virial(p, q) + f(p)*crd(q, i) + f(p + 3)*crd(q, j) + &
                        f(p + 6)*crd(q, k)
                end do
            end do

        end do

        return
    end subroutine AM_REG_ANGLES_get_ene_frc
!-------------------------------------------------------------
end module amoeba_reg_angles
!-------------------------------------------------------------
module amoeba_trig_angles

    use constants, only : pi

    implicit none
    private

    integer, save :: num_list = 0, startlist, endlist
    _REAL_, save :: energy, virial(3, 3)

    integer, allocatable, save :: list(:, :)

    _REAL_, allocatable, save :: force_constant(:), equil_value(:)
    _REAL_, allocatable, save :: fn(:), dfn_darg(:), arg(:), darg_dcrd(:, :)
    integer, save :: num_params = 0, do_flag = 0
    integer, save :: ftable_degree
    _REAL_, allocatable, save :: ftable_coeff(:)
    _REAL_, parameter ::  pt999999 = 0.999999d0
    _REAL_, parameter :: radians_to_degrees = 180.d0/pi
    _REAL_, parameter :: degrees_to_radians = pi/180.d0
    _REAL_ ::     xic, yic, zic, xkc, ykc, zkc, &
        ang, dang_dxic, dang_dyic, dang_dzic, &
        dang_dxkc, dang_dykc, dang_dzkc
    public AM_TRIG_ANGLES_readparm, AM_TRIG_ANGLES_deallocate, &
        AM_TRIG_ANGLES_set_user_bit, AM_TRIG_ANGLES_eval
#ifdef MPI
    public AM_TRIG_ANGLES_bcast
#endif
contains

#ifdef MPI
    subroutine AM_TRIG_ANGLES_bcast
        implicit none
        integer siztask, ierr

        include 'mpif.h'
#  include "extra.h"
#  include "do_flag.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_list, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_params, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(ftable_degree, 1, MPI_INTEGER, 0, commsander, ierr)

        if (num_list .eq. 0) return

        if (.not. master) then
            allocate (list(5, num_list), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (force_constant(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (equil_value(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (ftable_coeff(0:ftable_degree), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        call mpi_bcast(list, 5*num_list, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(force_constant, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(equil_value, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(ftable_coeff, 1 + ftable_degree, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

        if (.not. master) then
            call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
            allocate (fn(siztask), stat=ierr)
            allocate (arg(siztask), stat=ierr)
            allocate (dfn_darg(siztask), stat=ierr)
            allocate (darg_dcrd(12, siztask), stat=ierr)
        end if
    end subroutine AM_TRIG_ANGLES_bcast
#endif

!-----------------------------------------------------
    function AM_TRIG_ANGLES_readparm(nf)
        integer :: AM_TRIG_ANGLES_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, siztask, ier1, ier2, ier3, ier4

        AM_TRIG_ANGLES_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_TRIGONAL_ANGLE_', nf, num_list)
        if (num_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        dim1 = 5 ! 4 atoms plus param ptr
        allocate (list(dim1, num_list), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate list; num_list = ', num_list
            call mexit(6, 1)
        end if
        call AMOEBA_read_list_data('AMOEBA_TRIGONAL_ANGLE_', nf, dim1, num_list, list)
        call AM_VAL_get_num_params('AMOEBA_TRIGONAL_ANGLE_', nf, num_params)
        if (num_params <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        allocate (force_constant(num_params), stat=ier1)
        allocate (equil_value(num_params), stat=ier2)
        if ((ier1 /= 0) .or. (ier2 /= 0)) then
            write (6, *) 'failed to allocate parms; num_params = ', num_params
            call mexit(6, 1)
        end if
        call AM_VAL_read_force_constant('AMOEBA_TRIGONAL_ANGLE_', nf, &
            num_params, force_constant)
        call AM_VAL_read_equil_value('AMOEBA_TRIGONAL_ANGLE_', nf, &
            num_params, equil_value)
        call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
        allocate (fn(siztask), stat=ier1)
        allocate (arg(siztask), stat=ier2)
        allocate (dfn_darg(siztask), stat=ier3)
        allocate (darg_dcrd(12, siztask), stat=ier4) !4 atoms to get derivs wrt
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0)) then
            write (6, *) 'problems allocating angle scratch space '
            call mexit(6, 1)
        end if
        call AM_VAL_get_ftab_degree('AMOEBA_TRIGONAL_ANGLE_', nf, ftable_degree)
        allocate (ftable_coeff(0:ftable_degree), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate ftable_coeff; ftable_degree = ', &
                ftable_degree
            call mexit(6, 1)
        end if
        call AM_VAL_read_ftable_coeffs('AMOEBA_TRIGONAL_ANGLE_', nf, &
            ftable_degree, ftable_coeff)
        do_flag = ibset(do_flag, VALID_BIT)
        AM_TRIG_ANGLES_readparm = 1
    end function AM_TRIG_ANGLES_readparm
!-------------------------------------------------------------
!-------------------------------------------------------------
    subroutine AM_TRIG_ANGLES_deallocate()
        if (allocated(list)) deallocate (list)
        if (allocated(force_constant)) deallocate (force_constant)
        if (allocated(equil_value)) deallocate (equil_value)
        if (allocated(fn)) deallocate (fn)
        if (allocated(arg)) deallocate (arg)
        if (allocated(dfn_darg)) deallocate (dfn_darg)
        if (allocated(darg_dcrd)) deallocate (darg_dcrd)
        if (allocated(ftable_coeff)) deallocate (ftable_coeff)
    end subroutine AM_TRIG_ANGLES_deallocate
!-------------------------------------------------------------
    subroutine AM_TRIG_ANGLES_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_TRIG_ANGLES_set_user_bit
!-------------------------------------------------------------
    subroutine AM_TRIG_ANGLES_eval(crd, frc, ene, vir)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), vir(3, 3)
        _REAL_, intent(out) :: ene
#include "do_flag.h"

        integer j, k, nlist
! initialize
        energy = 0.d0
        do k = 1, 3
            do j = 1, 3
                virial(j, k) = 0.d0
            end do
        end do
        if (do_flag /= PROCEED) then
            return
        end if
        call AM_TRIG_ANGLES_get_args(crd, list, equil_value, arg, darg_dcrd)
        nlist = endlist - startlist + 1
        call AM_VAL_FTAB_eval_f_df(nlist, ftable_degree, ftable_coeff, &
            arg, fn, dfn_darg)
        call AM_TRIG_ANGLES_get_ene_frc(list, force_constant, &
            fn, dfn_darg, darg_dcrd, crd, frc)
        ene = energy
        do k = 1, 3
            do j = 1, 3
                vir(j, k) = vir(j, k) + virial(j, k)
            end do
        end do
    end subroutine AM_TRIG_ANGLES_eval
!-------------------------------------------------------------
    subroutine AM_TRIG_ANGLES_get_args(crd, alist, equil_value, arg, darg_dcrdijkl)
!--------------------------------------------------------
! This routine calculates angle function argument and its derivatives with
! respect to atomic positions of atom i,j,k,l.
!--------------------------------------------------------
! INPUT variables:
!    nangles:  number of angles in list
!    crd the atomic coord array
!    alist: 5 x nangles array giving for each angle the index of the first
!           atom, index of the second atom, index of the third,
!           index of the 4th for calc of projection
!           and then the param pointer
!    equil_value the list of ref angles
! OUTPUT variables:
!    arg, array of angle function args
!    darg_dcrdijkl   derivs of arg wrt crds of atom i,j,k,l
!--------------------------------------------------------
        _REAL_, intent(in) :: crd(3, *)
        integer, intent(in) :: alist(5, *)
        _REAL_, intent(in) :: equil_value(*)
        _REAL_, intent(out) :: arg(*), darg_dcrdijkl(12, *)
! LOCAL variables
        integer i, j, k, l, n, m, it, ind
        _REAL_ dotp
        _REAL_ xcen, ycen, zcen

        _REAL_ ang0
        _REAL_ dang_dxcen, dang_dycen, dang_dzcen
        _REAL_ v1(3), v2(3), v3(3), p(3), dp_dv1_p(3), dp_dv2_p(3)
        _REAL_ dotp1, dotp2, dxcen_dv1_p, dycen_dv1_p, dzcen_dv1_p, &
            dxcen_dv2_p, dycen_dv2_p, dzcen_dv2_p, dang_dv1_p, dang_dv2_p
! note units are degrees not radians...possibly change back later
        do n = startlist, endlist
            ind = n - startlist + 1
            i = alist(1, n)
            j = alist(2, n)
            k = alist(3, n)
            l = alist(4, n)
            it = alist(5, n)
            ang0 = equil_value(it)
!------first get the projected center xcen,ycen,zcen
            do m = 1, 3
                v1(m) = crd(m, i) - crd(m, l)
                v2(m) = crd(m, k) - crd(m, l)
                v3(m) = crd(m, j) - crd(m, l)
            end do
            call AM_VAL_VEC3D_get_perp_to_vecs(v1, v2, p, dp_dv1_p, dp_dv2_p)
            dotp = v3(1)*p(1) + v3(2)*p(2) + v3(3)*p(3)
            xcen = crd(1, j) - dotp*p(1)
            ycen = crd(2, j) - dotp*p(2)
            zcen = crd(3, j) - dotp*p(3)
            dotp1 = v3(1)*dp_dv1_p(1) + v3(2)*dp_dv1_p(2) + v3(3)*dp_dv1_p(3)
            dotp2 = v3(1)*dp_dv2_p(1) + v3(2)*dp_dv2_p(2) + v3(3)*dp_dv2_p(3)
            dxcen_dv1_p = -(dotp1*p(1) + dotp*dp_dv1_p(1))
            dycen_dv1_p = -(dotp1*p(2) + dotp*dp_dv1_p(2))
            dzcen_dv1_p = -(dotp1*p(3) + dotp*dp_dv1_p(3))
            dxcen_dv2_p = -(dotp2*p(1) + dotp*dp_dv2_p(1))
            dycen_dv2_p = -(dotp2*p(2) + dotp*dp_dv2_p(2))
            dzcen_dv2_p = -(dotp2*p(3) + dotp*dp_dv2_p(3))
            xic = crd(1, i) - xcen
            yic = crd(2, i) - ycen
            zic = crd(3, i) - zcen
            xkc = crd(1, k) - xcen
            ykc = crd(2, k) - ycen
            zkc = crd(3, k) - zcen
!--------next get angle as in a normal angle function
            call AM_TRIG_ANGLES_getang_cenproj()
!  derivative of angle wrt position of projected center (base of ic & kc)
            dang_dxcen = -(dang_dxic + dang_dxkc)
            dang_dycen = -(dang_dyic + dang_dykc)
            dang_dzcen = -(dang_dzic + dang_dzkc)
! extra deriv component due to effect of v1 on ang through cen
            dang_dv1_p = dang_dxcen*dxcen_dv1_p + dang_dycen*dycen_dv1_p + &
                dang_dzcen*dzcen_dv1_p
            dang_dv2_p = dang_dxcen*dxcen_dv2_p + dang_dycen*dycen_dv2_p + &
                dang_dzcen*dzcen_dv2_p
            arg(ind) = ang - ang0
! coordiindate v1_p given by dot product v1 and p ==> dv1_p_dv1 = p
            darg_dcrdijkl(1, ind) = dang_dxic + dang_dv1_p*p(1)
            darg_dcrdijkl(2, ind) = dang_dyic + dang_dv1_p*p(2)
            darg_dcrdijkl(3, ind) = dang_dzic + dang_dv1_p*p(3)
! indext for atom j
            darg_dcrdijkl(4, ind) = dang_dxcen
            darg_dcrdijkl(5, ind) = dang_dycen
            darg_dcrdijkl(6, ind) = dang_dzcen
! indext for atom k
            darg_dcrdijkl(7, ind) = dang_dxkc + dang_dv2_p*p(1)
            darg_dcrdijkl(8, ind) = dang_dykc + dang_dv2_p*p(2)
            darg_dcrdijkl(9, ind) = dang_dzkc + dang_dv2_p*p(3)
! deriv wrt crds of l due eindtirely to effect of v1,v2 on ang through cen
            darg_dcrdijkl(10, ind) = -(dang_dv1_p + dang_dv2_p)*p(1)
            darg_dcrdijkl(11, ind) = -(dang_dv1_p + dang_dv2_p)*p(2)
            darg_dcrdijkl(12, ind) = -(dang_dv1_p + dang_dv2_p)*p(3)
        end do
    end subroutine AM_TRIG_ANGLES_get_args
!----------------------------------------------------------
    subroutine AM_TRIG_ANGLES_getang_cenproj()
        implicit none

        _REAL_ :: dotp, lenic2, lenkc2, lenp, cosang, dang_dcosang
        _REAL_ :: dcosang_dxic, dcosang_dyic, dcosang_dzic
        _REAL_ :: dcosang_dxkc, dcosang_dykc, dcosang_dzkc

! cosine of angle is given by dot product of rij and rkj
!   divided by the product of the lengths of rij and rkj
        dotp = xic*xkc + yic*ykc + zic*zkc
        lenic2 = xic**2 + yic**2 + zic**2
        lenkc2 = xkc**2 + ykc**2 + zkc**2
        lenp = sqrt(lenic2*lenkc2)
        cosang = dotp/lenp
! avoid angle of pi and 0; really you need to use cosangle formulation
! near those points. however due to severe strain you could get bad values
! even for reference angle not near 0 or pi.
        cosang = max(-pt999999, cosang)
        cosang = min(pt999999, cosang)
        ang = radians_to_degrees*acos(cosang)
        dang_dcosang = -radians_to_degrees/sqrt(1.d0 - cosang**2)
! deriv of dotp wrt xic is xkc; deriv of lenp^-1 wrt xic is
!  lenkc^-1 * (-lenic^-2) * (xic/lenic) = -xic / (lenp*lenic2)
! similar for others
        dcosang_dxic = xkc/lenp - (dotp*xic)/(lenp*lenic2)
        dcosang_dyic = ykc/lenp - (dotp*yic)/(lenp*lenic2)
        dcosang_dzic = zkc/lenp - (dotp*zic)/(lenp*lenic2)
        dcosang_dxkc = xic/lenp - (dotp*xkc)/(lenp*lenkc2)
        dcosang_dykc = yic/lenp - (dotp*ykc)/(lenp*lenkc2)
        dcosang_dzkc = zic/lenp - (dotp*zkc)/(lenp*lenkc2)
! now use the chain rule
        dang_dxic = dang_dcosang*dcosang_dxic
        dang_dyic = dang_dcosang*dcosang_dyic
        dang_dzic = dang_dcosang*dcosang_dzic
        dang_dxkc = dang_dcosang*dcosang_dxkc
        dang_dykc = dang_dcosang*dcosang_dykc
        dang_dzkc = dang_dcosang*dcosang_dzkc
    end subroutine AM_TRIG_ANGLES_getang_cenproj
!----------------------------------------------------------
    subroutine AM_TRIG_ANGLES_get_ene_frc(alist, force_constant, func, dfunc, &
        darg_dcrdijkl, crd, frc)
        integer, intent(in) :: alist(5, *) ! 4 atoms plus param ptr
        _REAL_, intent(in) :: force_constant(*), func(*), dfunc(*), &
            darg_dcrdijkl(12, *), crd(3, *)
        _REAL_, intent(inout) :: frc(3, *)

        integer i, j, k, l, n, m, p, q, it, ind
        _REAL_ term, Frc_K, f(12)
        _REAL_ factor
        ! correct bending constant for units in degrees instead of radians
        factor = degrees_to_radians**2
        do n = startlist, endlist
            ind = n - startlist + 1
            i = alist(1, n)
            j = alist(2, n)
            k = alist(3, n)
            l = alist(4, n)
            it = alist(5, n)
            Frc_K = factor*force_constant(it)
            energy = energy + Frc_K*func(ind)
! apply chain rule to get deriv of energy with respect to crds of i,j,k,l
!  df holds the deriv of f with respect to its arg
!  while darg_dcrdijkl holds the derivs of arg with respect to crds of i,j,k,l
! recall force is negative of grad
            term = Frc_K*dfunc(ind)
            do m = 1, 12
                f(m) = term*darg_dcrdijkl(m, ind)
            end do
            frc(1, i) = frc(1, i) - f(1)
            frc(2, i) = frc(2, i) - f(2)
            frc(3, i) = frc(3, i) - f(3)
            frc(1, j) = frc(1, j) - f(4)
            frc(2, j) = frc(2, j) - f(5)
            frc(3, j) = frc(3, j) - f(6)
            frc(1, k) = frc(1, k) - f(7)
            frc(2, k) = frc(2, k) - f(8)
            frc(3, k) = frc(3, k) - f(9)
            frc(1, l) = frc(1, l) - f(10)
            frc(2, l) = frc(2, l) - f(11)
            frc(3, l) = frc(3, l) - f(12)
! now get virial
            do q = 1, 3
                do p = 1, 3
                    virial(p, q) = virial(p, q) + f(p)*crd(q, i) + f(p + 3)*crd(q, j) + &
                        f(p + 6)*crd(q, k) + f(p + 9)*crd(q, l)
                end do
            end do
        end do
    end subroutine AM_TRIG_ANGLES_get_ene_frc
!--------------------------------------------------------------
end module amoeba_trig_angles
!-------------------------------------------------------------
module amoeba_opbend_angles

    use constants, only : pi

    implicit none
    private

    integer, save :: num_list = 0, startlist, endlist
    _REAL_, save :: energy, virial(3, 3)

    integer, allocatable, save :: list(:, :)

    _REAL_, allocatable, save :: force_constant(:), equil_value(:)
    _REAL_, allocatable, save :: fn(:), dfn_darg(:), arg(:), darg_dcrd(:, :)
    integer, save :: num_params = 0, do_flag = 0
    integer, save :: ftable_degree
    _REAL_, allocatable, save :: ftable_coeff(:)
    _REAL_, parameter ::  pt999999 = 0.999999d0
    _REAL_, parameter :: radians_to_degrees = 180.d0/pi
    _REAL_, parameter :: degrees_to_radians = pi/180.d0
    _REAL_, parameter :: opbend_unit = 0.02191418d0
    public AM_OPBEND_ANGLES_readparm, AM_OPBEND_ANGLES_deallocate, &
        AM_OPBEND_ANGLES_set_user_bit, AM_OPBEND_ANGLES_eval
#ifdef MPI
    public AM_OPBEND_ANGLES_bcast
#endif
contains
!-------------------------------------------------------------

#ifdef MPI
    subroutine AM_OPBEND_ANGLES_bcast
        implicit none
        integer siztask, ierr

        include 'mpif.h'
#  include "extra.h"
#  include "do_flag.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_list, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_params, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(ftable_degree, 1, MPI_INTEGER, 0, commsander, ierr)

        if (num_list .eq. 0) return

        if (.not. master) then
            allocate (list(5, num_list), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (force_constant(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (ftable_coeff(0:ftable_degree), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        call mpi_bcast(list, 5*num_list, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(force_constant, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(ftable_coeff, 1 + ftable_degree, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

        if (.not. master) then
            call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
            allocate (fn(siztask), stat=ierr)
            allocate (arg(siztask), stat=ierr)
            allocate (dfn_darg(siztask), stat=ierr)
            allocate (darg_dcrd(12, siztask), stat=ierr)
        end if
    end subroutine AM_OPBEND_ANGLES_bcast
#endif

    function AM_OPBEND_ANGLES_readparm(nf)
        integer :: AM_OPBEND_ANGLES_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, siztask, ier1, ier2, ier3, ier4

        AM_OPBEND_ANGLES_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_OPBEND_ANGLE_', nf, num_list)
        if (num_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        dim1 = 5 ! 4 atoms plus param ptr
        allocate (list(dim1, num_list), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate list; num_list = ', num_list
            call mexit(6, 1)
        end if
        call AMOEBA_read_list_data('AMOEBA_OPBEND_ANGLE_', nf, dim1, num_list, list)
        call AM_VAL_get_num_params('AMOEBA_OPBEND_ANGLE_', nf, num_params)
        if (num_params <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        allocate (force_constant(num_params), stat=ier)
        if ((ier /= 0)) then
            write (6, *) 'failed to allocate parms; num_params = ', num_params
            call mexit(6, 1)
        end if
        call AM_VAL_read_force_constant('AMOEBA_OPBEND_ANGLE_', nf, &
            num_params, force_constant)
        call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
        allocate (fn(siztask), stat=ier1)
        allocate (arg(siztask), stat=ier2)
        allocate (dfn_darg(siztask), stat=ier3)
        allocate (darg_dcrd(12, siztask), stat=ier4) !4 atoms to get derivs wrt
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0)) then
            write (6, *) 'problems allocating angle scratch space '
            call mexit(6, 1)
        end if
        call AM_VAL_get_ftab_degree('AMOEBA_OPBEND_ANGLE_', nf, ftable_degree)
        allocate (ftable_coeff(0:ftable_degree), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate ftable_coeff; ftable_degree = ', &
                ftable_degree
            call mexit(6, 1)
        end if
        call AM_VAL_read_ftable_coeffs('AMOEBA_OPBEND_ANGLE_', nf, &
            ftable_degree, ftable_coeff)
        do_flag = ibset(do_flag, VALID_BIT)
        AM_OPBEND_ANGLES_readparm = 1
    end function AM_OPBEND_ANGLES_readparm
!-------------------------------------------------------------
!-------------------------------------------------------------
    subroutine AM_OPBEND_ANGLES_deallocate()
        if (allocated(list)) deallocate (list)
        if (allocated(force_constant)) deallocate (force_constant)
        if (allocated(fn)) deallocate (fn)
        if (allocated(arg)) deallocate (arg)
        if (allocated(dfn_darg)) deallocate (dfn_darg)
        if (allocated(darg_dcrd)) deallocate (darg_dcrd)
        if (allocated(ftable_coeff)) deallocate (ftable_coeff)
    end subroutine AM_OPBEND_ANGLES_deallocate
!-------------------------------------------------------------
    subroutine AM_OPBEND_ANGLES_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_OPBEND_ANGLES_set_user_bit
!-------------------------------------------------------------
    subroutine AM_OPBEND_ANGLES_eval(crd, frc, ene, vir)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), vir(3, 3)
        _REAL_, intent(out) :: ene
#include "do_flag.h"

        integer j, k, nlist
! initialize
        energy = 0.d0
        do k = 1, 3
            do j = 1, 3
                virial(j, k) = 0.d0
            end do
        end do
        if (do_flag /= PROCEED) then
            return
        end if
        call AM_OPBEND_ANGLES_get_args(crd, list, arg, darg_dcrd)
        nlist = endlist - startlist + 1
        call AM_VAL_FTAB_eval_f_df(nlist, ftable_degree, ftable_coeff, &
            arg, fn, dfn_darg)
        call AM_OPBEND_ANGLES_get_ene_frc(list, force_constant, &
            fn, dfn_darg, darg_dcrd, crd, frc)
        ene = energy
        do k = 1, 3
            do j = 1, 3
                vir(j, k) = vir(j, k) + virial(j, k)
            end do
        end do
    end subroutine AM_OPBEND_ANGLES_eval
!-------------------------------------------------------------
    subroutine AM_OPBEND_ANGLES_get_args(crd, alist, arg, darg_dcrd_ijkl)
!--------------------------------------------------------
! This routine calculates angle function argument and its derivatives with
! respect to atomic positions of atom i,j,k,l.
!--------------------------------------------------------
! INPUT variables:
!    nangles:  number of angles in list
!    crd the atomic coord array
!    alist: 5 x nangles array giving for each angle the index of the first
!           atom, index of the second atom, index of the third,
!           index of the 4th for calc of projection
!           and parameter pointer
! OUTPUT variables:
!    arg, array of angle function args
!    darg_dcrdijkl   derivs of arg wrt crds of atom i,j,k,l
!--------------------------------------------------------
        integer, intent(in) :: alist(5, *)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(out) :: arg(*), darg_dcrd_ijkl(12, *)

! LOCAL variables
        integer :: i, j, k, l, n, m, ind
        _REAL_ :: v1(3), v2(3), v3(3), p(3), dp_dv1_p(3), dp_dv2_p(3)
        _REAL_ :: dotp, lenjl2, lenjl
        _REAL_ :: ang, darg_dang, dang_dsinang, dsinang_dv1(3), &
            dsinang_dv2(3), dsinang_dv3(3), dsinang_dp(3), term1, term2, sinang
! note units are degrees not radians...possibly change back later
        do n = startlist, endlist
            ind = n - startlist + 1
            i = alist(1, n)
            j = alist(2, n)
            k = alist(3, n)
            l = alist(4, n)
            do m = 1, 3
                v1(m) = crd(m, i) - crd(m, l)
                v2(m) = crd(m, k) - crd(m, l)
                v3(m) = crd(m, j) - crd(m, l)
            end do
            call AM_VAL_VEC3D_get_perp_to_vecs(v1, v2, p, dp_dv1_p, dp_dv2_p)
            ! use arcsin formulation since angle is near 0
            dotp = p(1)*v3(1) + p(2)*v3(2) + p(3)*v3(3)
            lenjl2 = v3(1)*v3(1) + v3(2)*v3(2) + v3(3)*v3(3)
            lenjl = sqrt(lenjl2)
            sinang = dotp/lenjl
            sinang = max(-pt999999, sinang)
            sinang = min(pt999999, sinang)
            ang = radians_to_degrees*asin(sinang)
            dang_dsinang = radians_to_degrees/sqrt(1.d0 - sinang**2)
            arg(ind) = dabs(ang)  ! reference angle is zero
            darg_dang = 1.d0
            if (ang < 0) darg_dang = -1.d0
            do m = 1, 3
                dsinang_dv3(m) = p(m)/lenjl - (dotp*v3(m))/(lenjl*lenjl2)
                dsinang_dp(m) = v3(m)/lenjl - (dotp*p(m))/lenjl
            end do
            term1 = dp_dv1_p(1)*dsinang_dp(1) + dp_dv1_p(2)*dsinang_dp(2) + &
                dp_dv1_p(3)*dsinang_dp(3)
            term2 = dp_dv2_p(1)*dsinang_dp(1) + dp_dv2_p(2)*dsinang_dp(2) + &
                dp_dv2_p(3)*dsinang_dp(3)
            do m = 1, 3
                dsinang_dv1(m) = term1*p(m)
                dsinang_dv2(m) = term2*p(m)
            end do
            do m = 1, 3
                darg_dcrd_ijkl(m, ind) = darg_dang*dang_dsinang*dsinang_dv1(m)
                darg_dcrd_ijkl(m + 3, ind) = darg_dang*dang_dsinang*dsinang_dv3(m)
                darg_dcrd_ijkl(m + 6, ind) = darg_dang*dang_dsinang*dsinang_dv2(m)
                darg_dcrd_ijkl(m + 9, ind) = -(darg_dcrd_ijkl(m, ind) + &
                    darg_dcrd_ijkl(m + 3, ind) + &
                    darg_dcrd_ijkl(m + 6, ind))
            end do
        end do
    end subroutine AM_OPBEND_ANGLES_get_args
!----------------------------------------------------------
    subroutine AM_OPBEND_ANGLES_get_ene_frc(alist, force_constant, func, dfunc, &
        darg_dcrdijkl, crd, frc)
        integer, intent(in) :: alist(5, *) ! 4 atoms plus param ptr
        _REAL_, intent(in) :: force_constant(*), func(*), dfunc(*), &
            darg_dcrdijkl(12, *), crd(3, *)
        _REAL_, intent(inout) :: frc(3, *)

        integer :: i, j, k, l, n, m, p, q, ind, it
        _REAL_ :: term, Frc_K, f(12)

        do n = startlist, endlist
            ind = n - startlist + 1
            i = alist(1, n)
            j = alist(2, n)
            k = alist(3, n)
            l = alist(4, n)
            it = alist(5, n)
            Frc_K = opbend_unit*force_constant(it)
            energy = energy + Frc_K*func(ind)
! apply chain rule to get deriv of energy with respect to crds of i,j,k,l
!  df holds the deriv of f with respect to its arg
!  while darg_dcrdijkl holds the derivs of arg with respect to crds of i,j,k,l
! recall force is negative of grad
            term = Frc_K*dfunc(ind)
            do m = 1, 12
                f(m) = term*darg_dcrdijkl(m, ind)
            end do
            frc(1, i) = frc(1, i) - f(1)
            frc(2, i) = frc(2, i) - f(2)
            frc(3, i) = frc(3, i) - f(3)
            frc(1, j) = frc(1, j) - f(4)
            frc(2, j) = frc(2, j) - f(5)
            frc(3, j) = frc(3, j) - f(6)
            frc(1, k) = frc(1, k) - f(7)
            frc(2, k) = frc(2, k) - f(8)
            frc(3, k) = frc(3, k) - f(9)
            frc(1, l) = frc(1, l) - f(10)
            frc(2, l) = frc(2, l) - f(11)
            frc(3, l) = frc(3, l) - f(12)
! now get virial
            do q = 1, 3
                do p = 1, 3
                    virial(p, q) = virial(p, q) + f(p)*crd(q, i) + f(p + 3)*crd(q, j) + &
                        f(p + 6)*crd(q, k) + f(p + 9)*crd(q, l)
                end do
            end do
        end do
    end subroutine AM_OPBEND_ANGLES_get_ene_frc
!--------------------------------------------------------------
end module amoeba_opbend_angles
!-------------------------------------------------------------
module amoeba_torsions
    implicit none
    private

    integer, save :: num_list = 0, startlist, endlist
    _REAL_, save :: energy, virial(3, 3)

    integer, allocatable, save :: list(:, :)

    _REAL_, allocatable, save :: force_constant(:), periodicity(:), &
        cosphase(:), sinphase(:)
    _REAL_, allocatable, save :: fn(:), dfn_darg(:), cosarg(:), sinarg(:), &
        darg_dcrd(:, :)
    integer, save :: num_params = 0, do_flag = 0
    _REAL_, parameter :: torsion_unit = 0.5d0
    public AM_TORSIONS_readparm, AM_TORSIONS_deallocate, &
        AM_TORSIONS_set_user_bit, AM_TORSIONS_eval
#ifdef MPI
    public AM_TORSIONS_bcast
#endif
contains

#ifdef MPI
    subroutine AM_TORSIONS_bcast
        implicit none
        integer siztask, ierr

        include 'mpif.h'
#  include "extra.h"
#  include "do_flag.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_list, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_params, 1, MPI_INTEGER, 0, commsander, ierr)

        if (num_list .eq. 0) return

        if (.not. master) then
            allocate (list(5, num_list), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (force_constant(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (periodicity(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (cosphase(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (sinphase(num_params), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        call mpi_bcast(list, 5*num_list, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(force_constant, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(periodicity, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(cosphase, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(sinphase, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

        if (.not. master) then
            call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
            allocate (fn(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (cosarg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (sinarg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (dfn_darg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (darg_dcrd(12, siztask), stat=ierr)
            REQUIRE(ierr == 0)
        end if
    end subroutine AM_TORSIONS_bcast
#endif

!-------------------------------------------------------------
    function AM_TORSIONS_readparm(nf)
        integer :: AM_TORSIONS_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, siztask, ier1, ier2, ier3, ier4, ier5, n
        _REAL_, allocatable :: phase(:)

        AM_TORSIONS_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_TORSION_', nf, num_list)
        if (num_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        dim1 = 5 ! 4 atoms plus param ptr
        allocate (list(dim1, num_list), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate list; num_list = ', num_list
            call mexit(6, 1)
        end if
        call AMOEBA_read_list_data('AMOEBA_TORSION_', nf, dim1, num_list, list)
        call AM_VAL_get_num_params('AMOEBA_TORSION_', nf, num_params)
        if (num_params <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        allocate (force_constant(num_params), stat=ier1)
        allocate (periodicity(num_params), stat=ier2)
        allocate (phase(num_params), stat=ier3)
        allocate (cosphase(num_params), stat=ier4)
        allocate (sinphase(num_params), stat=ier5)
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) .or. (ier5 /= 0)) then
            write (6, *) 'failed to allocate parms; num_params = ', num_params
            call mexit(6, 1)
        end if
        call AM_VAL_read_force_constant('AMOEBA_TORSION_', nf, &
            num_params, force_constant)
        call AM_VAL_read_periodicity('AMOEBA_TORSION_', nf, &
            num_params, periodicity)
        call AM_VAL_read_phase('AMOEBA_TORSION_', nf, &
            num_params, phase)
        do n = 1, num_params
            cosphase(n) = cos(phase(n))
            sinphase(n) = sin(phase(n))
        end do
        deallocate (phase)
        call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
        allocate (fn(siztask), stat=ier1)
        allocate (cosarg(siztask), stat=ier2)
        allocate (sinarg(siztask), stat=ier2)
        allocate (dfn_darg(siztask), stat=ier3)
        allocate (darg_dcrd(12, siztask), stat=ier4) !4 atoms to get derivs wrt
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0)) then
            write (6, *) 'problems allocating angle scratch space '
            call mexit(6, 1)
        end if
        do_flag = ibset(do_flag, VALID_BIT)
        AM_TORSIONS_readparm = 1
    end function AM_TORSIONS_readparm
!-----------------------------------------------
!-------------------------------------------------------------
    subroutine AM_TORSIONS_deallocate()
        if (allocated(list)) deallocate (list)
        if (allocated(force_constant)) deallocate (force_constant)
        if (allocated(periodicity)) deallocate (periodicity)
        if (allocated(cosphase)) deallocate (cosphase)
        if (allocated(sinphase)) deallocate (sinphase)
        if (allocated(fn)) deallocate (fn)
        if (allocated(cosarg)) deallocate (cosarg)
        if (allocated(sinarg)) deallocate (sinarg)
        if (allocated(dfn_darg)) deallocate (dfn_darg)
        if (allocated(darg_dcrd)) deallocate (darg_dcrd)
    end subroutine AM_TORSIONS_deallocate
!-------------------------------------------------------------
    subroutine AM_TORSIONS_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_TORSIONS_set_user_bit
!-------------------------------------------------------------
    subroutine AM_TORSIONS_eval(crd, frc, ene, vir)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), vir(3, 3)
        _REAL_, intent(out) :: ene
#include "do_flag.h"

        integer j, k
! initialize
        energy = 0.d0
        do k = 1, 3
            do j = 1, 3
                virial(j, k) = 0.d0
            end do
        end do
        if (do_flag /= PROCEED) then
            return
        end if
        call AM_TORSIONS_get_args(crd, list, cosarg, sinarg, darg_dcrd)
        call AM_TORSIONS_func(list, force_constant, periodicity, cosphase, sinphase, &
            cosarg, sinarg, fn, dfn_darg)
        call AM_TORSIONS_get_ene_frc(list, fn, dfn_darg, darg_dcrd, crd, frc)
        ene = energy
        do k = 1, 3
            do j = 1, 3
                vir(j, k) = vir(j, k) + virial(j, k)
            end do
        end do
    end subroutine AM_TORSIONS_eval
!-----------------------------------------------
    subroutine AM_TORSIONS_get_args(crd, list, cosarg, sinarg, darg_dcrdijkl)
!--------------------------------------------------------
! This routine calculates torsion function argument and its derivatives with
! respect to atomic positions of atom i,j,k,l.
!--------------------------------------------------------
! INPUT variables:
!    nangles:  number of angles in list
!    crd the atomic coord array
!    list: 4 x ntorsions array giving for each torsion the index of the first
!           atom, index of the second atom, index of the third, fourth,
!           and index into the torsion parameter table giving the
!           18 terms of 6th order fourier expansion
! OUTPUT variables:
!    cosarg, array of cosines of phi angles
!    sinarg, array of sines of phi angles
!    darg_dcrdijkl   derivs of arg wrt crds of atom i,j,k,l
!--------------------------------------------------------
        _REAL_, intent(in) :: crd(3, *)
        integer, intent(in) :: list(5, *)
        _REAL_, intent(out) :: cosarg(*), sinarg(*), darg_dcrdijkl(12, *)

! LOCAL variables
        integer i, j, k, l, n, m, ind
        _REAL_ crd_abcd(12), gradphi_abcd(12), cosphi, sinphi

        do n = startlist, endlist
            ind = n - startlist + 1
            i = list(1, n)
            j = list(2, n)
            k = list(3, n)
            l = list(4, n)
            do m = 1, 3
                crd_abcd(m) = crd(m, i)
                crd_abcd(m + 3) = crd(m, j)
                crd_abcd(m + 6) = crd(m, k)
                crd_abcd(m + 9) = crd(m, l)
            end do
            call AM_VAL_GEOM_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)
            cosarg(ind) = cosphi
            sinarg(ind) = sinphi
            do m = 1, 12
                darg_dcrdijkl(m, ind) = gradphi_abcd(m)
            end do
        end do
    end subroutine AM_TORSIONS_get_args
!----------------------------------------------------------
    subroutine AM_TORSIONS_func(list, force_constant, periodicity, &
        cosphase, sinphase, cosarg, sinarg, func, dfunc_darg)
        integer, intent(in) :: list(5, *)
        _REAL_, intent(in) :: force_constant(*), periodicity(*)
        _REAL_, intent(in) :: cosphase(*), sinphase(*)
        _REAL_, intent(in) :: cosarg(*), sinarg(*)
        _REAL_, intent(out) :: func(*), dfunc_darg(*)

        integer :: ind, n, m, it, degree
        _REAL_ :: cosine(10), sine(10) ! maximum of degree is 10
        _REAL_ :: amplitude, period, cos_phase, sin_phase

        do n = startlist, endlist
            ind = n - startlist + 1
            cosine(1) = cosarg(ind)
            sine(1) = sinarg(ind)
            it = list(5, n)
            amplitude = torsion_unit*force_constant(it)
            period = periodicity(it)
            cos_phase = cosphase(it)
            sin_phase = sinphase(it)
            degree = nint(period)
            if (degree == 2) then
                cosine(2) = cosine(1)*cosine(1) - sine(1)*sine(1)
                sine(2) = sine(1)*cosine(1) + cosine(1)*sine(1)
            elseif (degree == 3) then
                cosine(2) = cosine(1)*cosine(1) - sine(1)*sine(1)
                sine(2) = sine(1)*cosine(1) + cosine(1)*sine(1)
                cosine(3) = cosine(2)*cosine(1) - sine(2)*sine(1)
                sine(3) = sine(2)*cosine(1) + cosine(2)*sine(1)
            elseif (degree > 3) then
                if (degree > 10) then
                    write (6, *) 'AM_TORSIONS_func: degreee too big: ', degree
                    call mexit(6, 1)
                end if
                do m = 2, degree
                    cosine(m) = cosine(m - 1)*cosine(1) - sine(m - 1)*sine(1)
                    sine(m) = sine(m - 1)*cosine(1) + cosine(m - 1)*sine(1)
                end do
            end if
            func(ind) = amplitude*(1.d0 + cos_phase*cosine(degree) + &
                sin_phase*sine(degree))
            dfunc_darg(ind) = amplitude*period*(sin_phase*cosine(degree) - &
                cos_phase*sine(degree))
        end do
    end subroutine AM_TORSIONS_func
!----------------------------------------------------------
    subroutine AM_TORSIONS_get_ene_frc(list, func, dfunc, darg_dcrdijkl, crd, frc)
        integer, intent(in) :: list(5, *) ! 4 atoms plus param ptr
        _REAL_, intent(in) :: func(*), dfunc(*), &
            darg_dcrdijkl(12, *), crd(3, *)
        _REAL_, intent(inout) :: frc(3, *)

        integer :: i, j, k, l, n, m, p, q, ind
        _REAL_ :: term, f(12)

        do n = startlist, endlist
            ind = n - startlist + 1
            i = list(1, n)
            j = list(2, n)
            k = list(3, n)
            l = list(4, n)
            energy = energy + func(ind)
! apply chain rule to get deriv of energy with respect to crds of i,j,k,l
!  df holds the deriv of f with respect to its arg
!  while darg_dcrdijkl holds the derivs of arg with respect to crds of i,j,k,l
! recall force is negative of grad
            term = dfunc(ind)
            do m = 1, 12
                f(m) = term*darg_dcrdijkl(m, ind)
            end do
            frc(1, i) = frc(1, i) - f(1)
            frc(2, i) = frc(2, i) - f(2)
            frc(3, i) = frc(3, i) - f(3)
            frc(1, j) = frc(1, j) - f(4)
            frc(2, j) = frc(2, j) - f(5)
            frc(3, j) = frc(3, j) - f(6)
            frc(1, k) = frc(1, k) - f(7)
            frc(2, k) = frc(2, k) - f(8)
            frc(3, k) = frc(3, k) - f(9)
            frc(1, l) = frc(1, l) - f(10)
            frc(2, l) = frc(2, l) - f(11)
            frc(3, l) = frc(3, l) - f(12)
! now get virial
            do q = 1, 3
                do p = 1, 3
                    virial(p, q) = virial(p, q) + f(p)*crd(q, i) + f(p + 3)*crd(q, j) + &
                        f(p + 6)*crd(q, k) + f(p + 9)*crd(q, l)
                end do
            end do
        end do
    end subroutine AM_TORSIONS_get_ene_frc
!--------------------------------------------------------------
end module amoeba_torsions
!-------------------------------------------------------------
module amoeba_stretch_torsions
    implicit none
    private

    integer, save :: num_list = 0, startlist, endlist
    _REAL_, save :: energy, virial(3, 3)

    integer, allocatable, save :: list(:, :)

    _REAL_, allocatable, save :: force_constant(:), periodicity(:), &
        cosphase(:), sinphase(:), bond_equil_value(:)
    _REAL_, allocatable, save :: fn(:), dfn_darg(:), cosarg(:), sinarg(:), &
        darg1_dcrd(:, :), arg2(:), darg2_dcrd(:, :)
    integer, save :: num_params = 0, do_flag = 0
    _REAL_, parameter :: stretch_tor_unit = 1.d0
    public AM_STRETCH_TORSIONS_readparm, AM_STRETCH_TORSIONS_deallocate, &
        AM_STRETCH_TORSIONS_suser_bit, AM_STRETCH_TORSIONS_eval
#ifdef MPI
    public AM_STRETCH_TORSIONS_bcast
#endif
contains

#ifdef MPI
    subroutine AM_STRETCH_TORSIONS_bcast
        implicit none
        integer siztask, ierr

        include 'mpif.h'
#  include "extra.h"
#  include "do_flag.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_list, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_params, 1, MPI_INTEGER, 0, commsander, ierr)
        if (num_list .eq. 0) return

        if (.not. master) then
            allocate (list(5, num_list), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (force_constant(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (periodicity(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (cosphase(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (sinphase(num_params), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        call mpi_bcast(list, 5*num_list, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(force_constant, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(periodicity, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(cosphase, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(sinphase, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

        if (.not. master) then
            call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
            allocate (fn(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (cosarg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (sinarg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (dfn_darg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (darg1_dcrd(12, siztask), stat=ierr)
            REQUIRE(ierr == 0) !4 atoms to get derivs wrt
            allocate (arg2(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (darg2_dcrd(6, siztask), stat=ierr)
            REQUIRE(ierr == 0) !2 atoms to get derivs wrt
        end if
    end subroutine AM_STRETCH_TORSIONS_bcast
#endif

!-------------------------------------------------------------
    function AM_STRETCH_TORSIONS_readparm(nf)
        integer :: AM_STRETCH_TORSIONS_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, siztask, n
        _REAL_, allocatable :: phase(:)

        AM_STRETCH_TORSIONS_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_STRETCH_TORSION_', nf, num_list)
        if (num_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        dim1 = 5 ! 4 atoms plus param ptr
        allocate (list(dim1, num_list), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate list; num_list = ', num_list
            call mexit(6, 1)
        end if
        call AMOEBA_read_list_data('AMOEBA_STRETCH_TORSION_', nf, dim1, num_list, list)
        call AM_VAL_get_num_params('AMOEBA_STRETCH_TORSION_', nf, num_params)
        if (num_params <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        !allocate
        allocate (force_constant(num_params), stat=ier)
        REQUIRE(ier == 0)
        allocate (periodicity(num_params), stat=ier)
        REQUIRE(ier == 0)
        allocate (phase(num_params), stat=ier)
        REQUIRE(ier == 0)
        allocate (cosphase(num_params), stat=ier)
        REQUIRE(ier == 0)
        allocate (sinphase(num_params), stat=ier)
        REQUIRE(ier == 0)
        allocate (bond_equil_value(num_params), stat=ier)
        REQUIRE(ier == 0)
        call AM_VAL_read_force_constant('AMOEBA_STRETCH_TORSION_', nf, &
            num_params, force_constant)
        call AM_VAL_read_periodicity('AMOEBA_STRETCH_TORSION_', nf, &
            num_params, periodicity)
        call AM_VAL_read_phase('AMOEBA_STRETCH_TORSION_', nf, &
            num_params, phase)
        call AM_VAL_read_equil_value('AMOEBA_STRETCH_TORSION_BOND_', nf, &
            num_params, bond_equil_value)
        do n = 1, num_params
            cosphase(n) = cos(phase(n))
            sinphase(n) = sin(phase(n))
        end do
        deallocate (phase)
        call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
        allocate (fn(siztask), stat=ier)
        REQUIRE(ier == 0)
        allocate (cosarg(siztask), stat=ier)
        REQUIRE(ier == 0)
        allocate (sinarg(siztask), stat=ier)
        REQUIRE(ier == 0)
        allocate (dfn_darg(siztask), stat=ier)
        REQUIRE(ier == 0)
        allocate (darg1_dcrd(12, siztask), stat=ier)
        REQUIRE(ier == 0) !4 atoms to get derivs wrt
        allocate (arg2(siztask), stat=ier)
        REQUIRE(ier == 0)
        allocate (darg2_dcrd(6, siztask), stat=ier)
        REQUIRE(ier == 0) !2 atoms to get derivs wrt
        do_flag = ibset(do_flag, VALID_BIT)
        AM_STRETCH_TORSIONS_readparm = 1
    end function AM_STRETCH_TORSIONS_readparm
!-----------------------------------------------
    subroutine AM_STRETCH_TORSIONS_deallocate()
        if (allocated(list)) deallocate (list)
        if (allocated(force_constant)) deallocate (force_constant)
        if (allocated(periodicity)) deallocate (periodicity)
        if (allocated(cosphase)) deallocate (cosphase)
        if (allocated(sinphase)) deallocate (sinphase)
        if (allocated(bond_equil_value)) deallocate (bond_equil_value)
        if (allocated(fn)) deallocate (fn)
        if (allocated(cosarg)) deallocate (cosarg)
        if (allocated(sinarg)) deallocate (sinarg)
        if (allocated(dfn_darg)) deallocate (dfn_darg)
        if (allocated(darg1_dcrd)) deallocate (darg1_dcrd)
        if (allocated(arg2)) deallocate (arg2)
        if (allocated(darg2_dcrd)) deallocate (darg2_dcrd)
    end subroutine AM_STRETCH_TORSIONS_deallocate
!-------------------------------------------------------------
    subroutine AM_STRETCH_TORSIONS_suser_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_STRETCH_TORSIONS_suser_bit
!-------------------------------------------------------------
    subroutine AM_STRETCH_TORSIONS_eval(crd, frc, ene, vir)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), vir(3, 3)
        _REAL_, intent(out) :: ene
#include "do_flag.h"

        integer j, k
! initialize
        energy = 0.d0
        do k = 1, 3
            do j = 1, 3
                virial(j, k) = 0.d0
            end do
        end do
        if (do_flag /= PROCEED) then
            return
        end if
        call AM_STRETCH_TORSIONS_get_tor_arg(crd, list, cosarg, sinarg, darg1_dcrd)
        call AM_STRETCH_TORSIONS_bond_arg(crd, list, bond_equil_value, &
            arg2, darg2_dcrd)
        call AM_STRETCH_TORSIONS_func(list, force_constant, periodicity, &
            cosphase, sinphase, cosarg, sinarg, fn, dfn_darg)
        call AM_STRETCH_TORSIONS_get_ene_frc(list, fn, dfn_darg, &
            darg1_dcrd, arg2, darg2_dcrd, crd, frc)
        ene = energy
        do k = 1, 3
            do j = 1, 3
                vir(j, k) = vir(j, k) + virial(j, k)
            end do
        end do
    end subroutine AM_STRETCH_TORSIONS_eval
!-----------------------------------------------
    subroutine AM_STRETCH_TORSIONS_get_tor_arg(crd, list, cosarg, sinarg, darg_dcrdijkl)
!--------------------------------------------------------
! This routine calculates torsion function argument and its derivatives with
! respect to atomic positions of atom i,j,k,l.
!--------------------------------------------------------
! INPUT variables:
!    nangles:  number of angles in list
!    crd the atomic coord array
!    list: 4 x ntorsions array giving for each torsion the index of the first
!           atom, index of the second atom, index of the third, fourth,
!           and index into the torsion parameter table giving the
!           18 terms of 6th order fourier expansion
! OUTPUT variables:
!    cosarg, array of cosines of phi angles
!    sinarg, array of sines of phi angles
!    darg_dcrdijkl   derivs of arg wrt crds of atom i,j,k,l
!--------------------------------------------------------
        _REAL_, intent(in) :: crd(3, *)
        integer, intent(in) :: list(5, *)
        _REAL_, intent(out) :: cosarg(*), sinarg(*), darg_dcrdijkl(12, *)

! LOCAL variables
        integer i, j, k, l, n, m, ind
        _REAL_ crd_abcd(12), gradphi_abcd(12), cosphi, sinphi

        do n = startlist, endlist
            ind = n - startlist + 1
            i = list(1, n)
            j = list(2, n)
            k = list(3, n)
            l = list(4, n)
            do m = 1, 3
                crd_abcd(m) = crd(m, i)
                crd_abcd(m + 3) = crd(m, j)
                crd_abcd(m + 6) = crd(m, k)
                crd_abcd(m + 9) = crd(m, l)
            end do
            call AM_VAL_GEOM_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)
            cosarg(ind) = cosphi
            sinarg(ind) = sinphi
            do m = 1, 12
                darg_dcrdijkl(m, ind) = gradphi_abcd(m)
            end do
        end do
    end subroutine AM_STRETCH_TORSIONS_get_tor_arg
!----------------------------------------------------------
    subroutine AM_STRETCH_TORSIONS_bond_arg(crd, list, equil_value, arg, darg_dcrd)
!--------------------------------------------------------
! This routine calculates bond function argument and its derivatives with
! respect to atomic positions of atom j,k. (tor list atoms are i,j,k,l
! This gets used in a stretch-torsion energy routine
! Note that derivatives with respect
! to atom k are just the negative of those with respect to atom i
! NOTE that periodic boundary conditions are not used here
!      i.e. imaging is done on a per molecule basis
!--------------------------------------------------------
! INPUT variables:
!    startlist,endlist  the set of bond indices to process--now global variates
!    crd the atomic coord array
!    list: 3 x nbond array giving for each bond the index of the first
!           atom, index of the second atom, and index into the bond
!           parameter table giving the equilibrium length
!    equil_value the list of ref bond lengths
! OUTPUT variables:
!    arg, array of bond function args
!    darg_dcrd   derivs of arg wrt crds of atom j,k
!--------------------------------------------------------

        _REAL_, intent(in) :: crd(3, *)
        integer, intent(in) :: list(5, *)
        _REAL_, intent(in) :: equil_value(*)
        _REAL_, intent(out) :: arg(*), darg_dcrd(6, *)

        integer j, k, n, ind, it
        _REAL_ bl, dx, dy, dz

        do n = startlist, endlist
            ind = n - startlist + 1
            j = list(2, n)
            k = list(3, n)
            it = list(5, n)
            dx = crd(1, j) - crd(1, k)
            dy = crd(2, j) - crd(2, k)
            dz = crd(3, j) - crd(3, k)
            bl = sqrt(dx*dx + dy*dy + dz*dz)
! recall the reference bond length is 2nd parameter
            arg(ind) = bl - equil_value(it)
! differentiate bl to get darg_dcrd
            darg_dcrd(1, ind) = dx/bl
            darg_dcrd(2, ind) = dy/bl
            darg_dcrd(3, ind) = dz/bl
            darg_dcrd(4, ind) = -dx/bl
            darg_dcrd(5, ind) = -dy/bl
            darg_dcrd(6, ind) = -dz/bl
        end do

        return
    end subroutine AM_STRETCH_TORSIONS_bond_arg
!------------------------------------------------------------
    subroutine AM_STRETCH_TORSIONS_func(list, force_constant, periodicity, &
        cosphase, sinphase, cosarg, sinarg, func, dfunc_darg)
        integer, intent(in) :: list(5, *)
        _REAL_, intent(in) :: force_constant(*), periodicity(*)
        _REAL_, intent(in) :: cosphase(*), sinphase(*)
        _REAL_, intent(in) :: cosarg(*), sinarg(*)
        _REAL_, intent(out) :: func(*), dfunc_darg(*)

        integer :: ind, n, m, it, degree
        _REAL_ :: cosine(10), sine(10) ! maximum of degree is 10
        _REAL_ :: amplitude, period, cos_phase, sin_phase

        do n = startlist, endlist
            ind = n - startlist + 1
            cosine(1) = cosarg(ind)
            sine(1) = sinarg(ind)
            it = list(5, n)
            amplitude = force_constant(it)
            period = periodicity(it)
            cos_phase = cosphase(it)
            sin_phase = sinphase(it)
            degree = nint(period)
            if (degree == 2) then
                cosine(2) = cosine(1)*cosine(1) - sine(1)*sine(1)
                sine(2) = sine(1)*cosine(1) + cosine(1)*sine(1)
            elseif (degree == 3) then
                cosine(2) = cosine(1)*cosine(1) - sine(1)*sine(1)
                sine(2) = sine(1)*cosine(1) + cosine(1)*sine(1)
                cosine(3) = cosine(2)*cosine(1) - sine(2)*sine(1)
                sine(3) = sine(2)*cosine(1) + cosine(2)*sine(1)
            elseif (degree > 3) then
                if (degree > 10) then
                    write (6, *) 'AM_TORSIONS_func: degreee too big: ', degree
                    call mexit(6, 1)
                end if
                do m = 2, degree
                    cosine(m) = cosine(m - 1)*cosine(1) - sine(m - 1)*sine(1)
                    sine(m) = sine(m - 1)*cosine(1) + cosine(m - 1)*sine(1)
                end do
            end if
            func(ind) = amplitude*(1.d0 + cos_phase*cosine(degree) + &
                sin_phase*sine(degree))
            dfunc_darg(ind) = amplitude*period*(sin_phase*cosine(degree) - &
                cos_phase*sine(degree))
        end do
    end subroutine AM_STRETCH_TORSIONS_func
!----------------------------------------------------------
    subroutine AM_STRETCH_TORSIONS_get_ene_frc(list, fn, dfn_darg1, &
        darg1_dcrd, arg2, darg2_dcrd, crd, frc)
        integer, intent(in) :: list(5, *) ! 4 atoms plus param ptr
        _REAL_, intent(in) :: fn(*), dfn_darg1(*), darg1_dcrd(12, *), arg2(*), &
            darg2_dcrd(6, *), crd(3, *)
        _REAL_, intent(inout) :: frc(3, *)

        integer :: i, j, k, l, n, m, p, q, ind
        _REAL_ :: term, f(12)

        do n = startlist, endlist
            ind = n - startlist + 1
            i = list(1, n)
            j = list(2, n)
            k = list(3, n)
            l = list(4, n)
            energy = energy + stretch_tor_unit*arg2(ind)*fn(ind)
! apply product rule and then chain rule to get deriv of energy with respect
!  to crds of i,j,k,l
!  dfn_darg1 holds the deriv of fn with respect to its arg
!  while darg1_dcrd holds the derivs of arg with respect to crds of i,j,k,l
! recall force is negative of grad
            term = stretch_tor_unit*arg2(ind)*dfn_darg1(ind)
            do m = 1, 12
                f(m) = term*darg1_dcrd(m, ind)
            end do
            ! 2nd term from product rule
            term = stretch_tor_unit*fn(ind)
            do m = 1, 6
                f(m + 3) = f(m + 3) + term*darg2_dcrd(m, ind)
            end do
            frc(1, i) = frc(1, i) - f(1)
            frc(2, i) = frc(2, i) - f(2)
            frc(3, i) = frc(3, i) - f(3)
            frc(1, j) = frc(1, j) - f(4)
            frc(2, j) = frc(2, j) - f(5)
            frc(3, j) = frc(3, j) - f(6)
            frc(1, k) = frc(1, k) - f(7)
            frc(2, k) = frc(2, k) - f(8)
            frc(3, k) = frc(3, k) - f(9)
            frc(1, l) = frc(1, l) - f(10)
            frc(2, l) = frc(2, l) - f(11)
            frc(3, l) = frc(3, l) - f(12)
! now get virial
            do q = 1, 3
                do p = 1, 3
                    virial(p, q) = virial(p, q) + f(p)*crd(q, i) + f(p + 3)*crd(q, j) + &
                        f(p + 6)*crd(q, k) + f(p + 9)*crd(q, l)
                end do
            end do
        end do
    end subroutine AM_STRETCH_TORSIONS_get_ene_frc
!--------------------------------------------------------------
end module amoeba_stretch_torsions
!-----------------------------------------------
module amoeba_pitorsions
    implicit none
    private

    integer, save :: num_list = 0, startlist, endlist
    _REAL_, save :: energy, virial(3, 3)

    integer, allocatable, save :: list(:, :)

    _REAL_, allocatable, save :: force_constant(:), periodicity(:), &
        cosphase(:), sinphase(:)
    _REAL_, allocatable, save :: fn(:), dfn_darg(:), cosarg(:), sinarg(:), &
        darg_dcrd(:, :)
    integer, save :: num_params = 0, do_flag = 0
    public AM_PITORSIONS_readparm, AM_PITORSIONS_deallocate, &
        AM_PITORSIONS_set_user_bit, AM_PITORSIONS_eval
#ifdef MPI
    public AM_PITORSIONS_bcast
#endif
contains

#ifdef MPI
    subroutine AM_PITORSIONS_bcast
        implicit none
        integer siztask, ierr

        include 'mpif.h'
#  include "extra.h"
#  include "do_flag.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_list, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_params, 1, MPI_INTEGER, 0, commsander, ierr)
        if (num_list .eq. 0) return

        if (.not. master) then
            allocate (list(7, num_list), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (force_constant(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (periodicity(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (cosphase(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (sinphase(num_params), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        call mpi_bcast(list, 7*num_list, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(force_constant, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(periodicity, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(cosphase, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(sinphase, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

        if (.not. master) then
            call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
            allocate (fn(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (cosarg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (sinarg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (dfn_darg(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (darg_dcrd(18, siztask), stat=ierr)
            REQUIRE(ierr == 0)
        end if
    end subroutine AM_PITORSIONS_bcast
#endif

!-------------------------------------------------------------
    function AM_PITORSIONS_readparm(nf)
        integer :: AM_PITORSIONS_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, siztask, ier1, ier2, ier3, ier4, ier5, n
        _REAL_, allocatable :: phase(:)

        AM_PITORSIONS_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_PI_TORSION_', nf, num_list)
        if (num_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        dim1 = 7 ! 6 atoms plus param ptr
        allocate (list(dim1, num_list), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate list; num_list = ', num_list
            call mexit(6, 1)
        end if
        call AMOEBA_read_list_data('AMOEBA_PI_TORSION_', nf, dim1, num_list, list)
        call AM_VAL_get_num_params('AMOEBA_PI_TORSION_', nf, num_params)
        if (num_params <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        allocate (force_constant(num_params), stat=ier1)
        allocate (periodicity(num_params), stat=ier2)
        allocate (phase(num_params), stat=ier3)
        allocate (cosphase(num_params), stat=ier4)
        allocate (sinphase(num_params), stat=ier5)
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) .or. (ier5 /= 0)) then
            write (6, *) 'failed to allocate parms; num_params = ', num_params
            call mexit(6, 1)
        end if
        call AM_VAL_read_force_constant('AMOEBA_PI_TORSION_', nf, &
            num_params, force_constant)
        call AM_VAL_read_periodicity('AMOEBA_PI_TORSION_', nf, &
            num_params, periodicity)
        call AM_VAL_read_phase('AMOEBA_PI_TORSION_', nf, &
            num_params, phase)
        do n = 1, num_params
            cosphase(n) = cos(phase(n))
            sinphase(n) = sin(phase(n))
        end do
        deallocate (phase)
        call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
        allocate (fn(siztask), stat=ier1)
        allocate (cosarg(siztask), stat=ier2)
        allocate (sinarg(siztask), stat=ier2)
        allocate (dfn_darg(siztask), stat=ier3)
        allocate (darg_dcrd(18, siztask), stat=ier4) !6 atoms to get derivs wrt
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0)) then
            write (6, *) 'problems allocating angle scratch space '
            call mexit(6, 1)
        end if
        do_flag = ibset(do_flag, VALID_BIT)
        AM_PITORSIONS_readparm = 1
    end function AM_PITORSIONS_readparm
!-----------------------------------------------
!-------------------------------------------------------------
    subroutine AM_PITORSIONS_deallocate()
        if (allocated(list)) deallocate (list)
        if (allocated(force_constant)) deallocate (force_constant)
        if (allocated(periodicity)) deallocate (periodicity)
        if (allocated(cosphase)) deallocate (cosphase)
        if (allocated(sinphase)) deallocate (sinphase)
        if (allocated(fn)) deallocate (fn)
        if (allocated(cosarg)) deallocate (cosarg)
        if (allocated(sinarg)) deallocate (sinarg)
        if (allocated(dfn_darg)) deallocate (dfn_darg)
        if (allocated(darg_dcrd)) deallocate (darg_dcrd)
    end subroutine AM_PITORSIONS_deallocate
!-------------------------------------------------------------
    subroutine AM_PITORSIONS_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_PITORSIONS_set_user_bit
!-------------------------------------------------------------
    subroutine AM_PITORSIONS_eval(crd, frc, ene, vir)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), vir(3, 3)
        _REAL_, intent(out) :: ene
#include "do_flag.h"

        integer j, k
        energy = 0.d0
        do k = 1, 3
            do j = 1, 3
                virial(j, k) = 0.d0
            end do
        end do
        if (do_flag /= PROCEED) then
            return
        end if
        call AM_PITORSIONS_get_args(crd, list, cosarg, sinarg, darg_dcrd)
        call AM_PITORSIONS_func(list, force_constant, periodicity, &
            cosphase, sinphase, cosarg, sinarg, fn, dfn_darg)
        call AM_PITORSIONS_get_ene_frc(list, fn, dfn_darg, darg_dcrd, crd, frc)
        ene = energy
        do k = 1, 3
            do j = 1, 3
                vir(j, k) = vir(j, k) + virial(j, k)
            end do
        end do
    end subroutine AM_PITORSIONS_eval
!-----------------------------------------------
    subroutine AM_PITORSIONS_get_args(crd, list, cosarg, sinarg, darg_dcrd_ijklmn)
        _REAL_, intent(in) :: crd(3, *)
        integer, intent(in) :: list(7, *) !6 atoms plus parm ptr
        _REAL_, intent(out) :: cosarg(*), sinarg(*), darg_dcrd_ijklmn(18, *)

        integer :: h, i, j, k, l, m, n, p, ind
        _REAL_ :: ril(3), rjl(3), rlk(3), rmk(3), rnk(3), p_ijl(3), p_mnk(3), &
            dp_ijl_dril_p(3), dp_ijl_drjl_p(3), &
            dp_mnk_drmk_p(3), dp_mnk_drnk_p(3)
        _REAL_ :: crd_abcd(12), gradphi_abcd(12), cosphi, sinphi
        _REAL_ :: termi, termj, termm, termn

        do h = startlist, endlist
            ind = h - startlist + 1
            i = list(1, h)
            j = list(2, h)
            k = list(3, h)
            l = list(4, h)
            m = list(5, h)
            n = list(6, h)
            do p = 1, 3
                ril(p) = crd(p, i) - crd(p, l)
                rjl(p) = crd(p, j) - crd(p, l)
                rlk(p) = crd(p, l) - crd(p, k)
                rmk(p) = crd(p, m) - crd(p, k)
                rnk(p) = crd(p, n) - crd(p, k)
            end do
            call AM_VAL_VEC3D_get_perp_to_vecs(ril, rjl, p_ijl, dp_ijl_dril_p, dp_ijl_drjl_p)
            call AM_VAL_VEC3D_get_perp_to_vecs(rmk, rnk, p_mnk, dp_mnk_drmk_p, dp_mnk_drnk_p)
            ! use the regular torsion results for artificial sites a,b,c,d
            do p = 1, 3
                crd_abcd(p) = crd(p, k) + p_ijl(p)
                crd_abcd(p + 3) = crd(p, k)
                crd_abcd(p + 6) = crd(p, l)
                crd_abcd(p + 9) = crd(p, l) + p_mnk(p)
            end do
            call AM_VAL_GEOM_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)
            cosarg(ind) = cosphi
            sinarg(ind) = sinphi
            ! transfer grad wrt point a to i,j,l and grad wrt point d to m,n,k
            ! recall that movement of i parallel to p_ijk causes change in p_ijk.
            ! hence a in direction dp_ijl_dril_p
            termi = dp_ijl_dril_p(1)*gradphi_abcd(1) + &
                dp_ijl_dril_p(2)*gradphi_abcd(2) + &
                dp_ijl_dril_p(3)*gradphi_abcd(3)
            termj = dp_ijl_drjl_p(1)*gradphi_abcd(1) + &
                dp_ijl_drjl_p(2)*gradphi_abcd(2) + &
                dp_ijl_drjl_p(3)*gradphi_abcd(3)
            termm = dp_mnk_drmk_p(1)*gradphi_abcd(10) + &
                dp_mnk_drmk_p(2)*gradphi_abcd(11) + &
                dp_mnk_drmk_p(3)*gradphi_abcd(12)
            termn = dp_mnk_drnk_p(1)*gradphi_abcd(10) + &
                dp_mnk_drnk_p(2)*gradphi_abcd(11) + &
                dp_mnk_drnk_p(3)*gradphi_abcd(12)
            do p = 1, 3
                darg_dcrd_ijklmn(p, ind) = termi*p_ijl(p)  ! for i
                darg_dcrd_ijklmn(3 + p, ind) = termj*p_ijl(p) ! for j
                darg_dcrd_ijklmn(6 + p, ind) = gradphi_abcd(p) + gradphi_abcd(p + 3) - &
                    (termm + termn)*p_mnk(p) ! for k
                darg_dcrd_ijklmn(9 + p, ind) = gradphi_abcd(p + 6) + gradphi_abcd(p + 9) - &
                    (termi + termj)*p_ijl(p) ! for l
                darg_dcrd_ijklmn(12 + p, ind) = termm*p_mnk(p)  ! for m
                darg_dcrd_ijklmn(15 + p, ind) = termn*p_mnk(p)  ! for n
            end do
        end do
    end subroutine AM_PITORSIONS_get_args
!---------------------------------------------------------------
    subroutine AM_PITORSIONS_func(list, force_constant, periodicity, &
        cosphase, sinphase, cosarg, sinarg, func, dfunc_darg)
        integer, intent(in) :: list(7, *)! 6 atoms plus parm ptr
        _REAL_, intent(in) :: force_constant(*), periodicity(*)
        _REAL_, intent(in) :: cosphase(*), sinphase(*)
        _REAL_, intent(in) :: cosarg(*), sinarg(*)
        _REAL_, intent(out) :: func(*), dfunc_darg(*)

        integer :: ind, n, m, it, degree
        _REAL_ :: cosine(10), sine(10) ! maximum of degree is 10
        _REAL_ :: amplitude, period, cos_phase, sin_phase

        do n = startlist, endlist
            ind = n - startlist + 1
            cosine(1) = cosarg(ind)
            sine(1) = sinarg(ind)
            it = list(7, n)
            amplitude = force_constant(it)
            period = periodicity(it)
            cos_phase = cosphase(it)
            sin_phase = sinphase(it)
            degree = nint(period)
            if (degree == 2) then
                cosine(2) = cosine(1)*cosine(1) - sine(1)*sine(1)
                sine(2) = sine(1)*cosine(1) + cosine(1)*sine(1)
            elseif (degree == 3) then
                cosine(2) = cosine(1)*cosine(1) - sine(1)*sine(1)
                sine(2) = sine(1)*cosine(1) + cosine(1)*sine(1)
                cosine(3) = cosine(2)*cosine(1) - sine(2)*sine(1)
                sine(3) = sine(2)*cosine(1) + cosine(2)*sine(1)
            elseif (degree > 3) then
                if (degree > 10) then
                    write (6, *) 'AM_PITORSIONS_func: degreee too big: ', degree
                    call mexit(6, 1)
                end if
                do m = 2, degree
                    cosine(m) = cosine(m - 1)*cosine(1) - sine(m - 1)*sine(1)
                    sine(m) = sine(m - 1)*cosine(1) + cosine(m - 1)*sine(1)
                end do
            end if
            func(ind) = amplitude*(1.d0 + cos_phase*cosine(degree) + &
                sin_phase*sine(degree))
            dfunc_darg(ind) = amplitude*period*(sin_phase*cosine(degree) - &
                cos_phase*sine(degree))
        end do
    end subroutine AM_PITORSIONS_func
!----------------------------------------------------------
    subroutine AM_PITORSIONS_get_ene_frc(list, func, dfunc, darg_dcrd_ijklmn, crd, frc)
        integer, intent(in) :: list(7, *) ! 6 atoms plus param ptr
        _REAL_, intent(in) :: func(*), dfunc(*), &
            darg_dcrd_ijklmn(18, *), crd(3, *)
        _REAL_, intent(inout) :: frc(3, *)

        integer :: i, j, k, l, n, m, h, p, q, ind
        _REAL_ :: term, f(18)

        do h = startlist, endlist
            ind = h - startlist + 1
            i = list(1, h)
            j = list(2, h)
            k = list(3, h)
            l = list(4, h)
            m = list(5, h)
            n = list(6, h)
            energy = energy + func(ind)
! apply chain rule to get deriv of energy with respect to crds of i,j,k,l
!  df holds the deriv of f with respect to its arg
!  while darg_dcrdijkl holds the derivs of arg with respect to crds of i,j,k,l
! recall force is negative of grad
            term = dfunc(ind)
            do p = 1, 18
                f(p) = term*darg_dcrd_ijklmn(p, ind)
            end do
            frc(1, i) = frc(1, i) - f(1)
            frc(2, i) = frc(2, i) - f(2)
            frc(3, i) = frc(3, i) - f(3)
            frc(1, j) = frc(1, j) - f(4)
            frc(2, j) = frc(2, j) - f(5)
            frc(3, j) = frc(3, j) - f(6)
            frc(1, k) = frc(1, k) - f(7)
            frc(2, k) = frc(2, k) - f(8)
            frc(3, k) = frc(3, k) - f(9)
            frc(1, l) = frc(1, l) - f(10)
            frc(2, l) = frc(2, l) - f(11)
            frc(3, l) = frc(3, l) - f(12)
            frc(1, m) = frc(1, m) - f(13)
            frc(2, m) = frc(2, m) - f(14)
            frc(3, m) = frc(3, m) - f(15)
            frc(1, n) = frc(1, n) - f(16)
            frc(2, n) = frc(2, n) - f(17)
            frc(3, n) = frc(3, n) - f(18)
! now get virial
            do q = 1, 3
                do p = 1, 3
                    virial(p, q) = virial(p, q) + f(p)*crd(q, i) + f(p + 3)*crd(q, j) + &
                        f(p + 6)*crd(q, k) + f(p + 9)*crd(q, l) + &
                        f(p + 12)*crd(q, m) + f(p + 15)*crd(q, n)
                end do
            end do
        end do
    end subroutine AM_PITORSIONS_get_ene_frc
!---------------------------------------------------------------
end module amoeba_pitorsions
!-----------------------------------------------
module amoeba_stretch_bend

    use constants, only : pi
    implicit none
    private

    integer, save :: num_list = 0, startlist, endlist
    _REAL_, save :: energy, virial(3, 3)

    integer, allocatable, save :: list(:, :)

    _REAL_, allocatable, save :: force_constant(:), angle_equil_value(:)
    _REAL_, allocatable, save :: bond1_equil_value(:), bond2_equil_value(:)
    _REAL_, allocatable, save :: arg1(:), darg1_dcrd(:, :)
    _REAL_, allocatable, save :: arg2(:), darg2_dcrd(:, :)
    integer, save :: num_params = 0, do_flag = 0
    _REAL_, parameter ::  pt999999 = 0.999999d0
    _REAL_, parameter :: radians_to_degrees = 180.d0/pi
    _REAL_, parameter :: degrees_to_radians = pi/180.d0
    public AM_STRETCH_BEND_readparm, AM_STRETCH_BEND_deallocate, &
        AM_STRETCH_BEND_set_user_bit, AM_STRETCH_BEND_eval
#ifdef MPI
    public AM_STRETCH_BEND_bcast
#endif
contains

#ifdef MPI
    subroutine AM_STRETCH_BEND_bcast
        implicit none
        integer siztask, ierr

        include 'mpif.h'
#  include "extra.h"
#  include "do_flag.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_list, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_params, 1, MPI_INTEGER, 0, commsander, ierr)

        if (num_list .eq. 0) return

        if (.not. master) then
            allocate (list(4, num_list), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (force_constant(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (angle_equil_value(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (bond1_equil_value(num_params), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (bond2_equil_value(num_params), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        call mpi_bcast(list, 4*num_list, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(force_constant, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(angle_equil_value, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(bond1_equil_value, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(bond2_equil_value, num_params, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

        if (.not. master) then
            call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
            allocate (arg1(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (arg2(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (darg1_dcrd(9, siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (darg2_dcrd(9, siztask), stat=ierr)
            REQUIRE(ierr == 0)
        end if
    end subroutine AM_STRETCH_BEND_bcast
#endif

!-----------------------------------------------------
    function AM_STRETCH_BEND_readparm(nf)
        integer :: AM_STRETCH_BEND_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, siztask, ier1, ier2, ier3, ier4

        AM_STRETCH_BEND_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_STRETCH_BEND_', nf, num_list)
        if (num_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        dim1 = 4 ! 3 atoms plus param ptr
        allocate (list(dim1, num_list), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate list; num_list = ', num_list
            call mexit(6, 1)
        end if
        call AMOEBA_read_list_data('AMOEBA_STRETCH_BEND_', nf, dim1, num_list, list)
        call AM_VAL_get_num_params('AMOEBA_STRETCH_BEND_', nf, num_params)
        if (num_params <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        allocate (force_constant(num_params), stat=ier1)
        allocate (angle_equil_value(num_params), stat=ier2)
        allocate (bond1_equil_value(num_params), stat=ier3)
        allocate (bond2_equil_value(num_params), stat=ier4)
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0)) then
            write (6, *) 'failed to allocate parms; num_params = ', num_params
            call mexit(6, 1)
        end if
        call AM_VAL_read_force_constant('AMOEBA_STRETCH_BEND_', nf, &
            num_params, force_constant)
        call AM_VAL_read_equil_value('AMOEBA_STRETCH_BEND_ANGLE_', nf, &
            num_params, angle_equil_value)
        call AM_VAL_read_equil_value('AMOEBA_STRETCH_BEND_BOND1_', nf, &
            num_params, bond1_equil_value)
        call AM_VAL_read_equil_value('AMOEBA_STRETCH_BEND_BOND2_', nf, &
            num_params, bond2_equil_value)
        call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
        allocate (arg1(siztask), stat=ier1)
        allocate (arg2(siztask), stat=ier2)
        allocate (darg1_dcrd(9, siztask), stat=ier3) !3 atoms to get derivs wrt
        allocate (darg2_dcrd(9, siztask), stat=ier4) !3 atoms to get derivs wrt
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0)) then
            write (6, *) 'problems allocating stretch-bend scratch space '
            call mexit(6, 1)
        end if
        do_flag = ibset(do_flag, VALID_BIT)
        AM_STRETCH_BEND_readparm = 1
    end function AM_STRETCH_BEND_readparm
!-------------------------------------------------------------
!-------------------------------------------------------------
    subroutine AM_STRETCH_BEND_deallocate()
        if (allocated(list)) deallocate (list)
        if (allocated(force_constant)) deallocate (force_constant)
        if (allocated(angle_equil_value)) deallocate (angle_equil_value)
        if (allocated(bond1_equil_value)) deallocate (bond1_equil_value)
        if (allocated(bond2_equil_value)) deallocate (bond2_equil_value)
        if (allocated(arg1)) deallocate (arg1)
        if (allocated(darg1_dcrd)) deallocate (darg1_dcrd)
        if (allocated(arg2)) deallocate (arg2)
        if (allocated(darg2_dcrd)) deallocate (darg2_dcrd)
    end subroutine AM_STRETCH_BEND_deallocate
!-------------------------------------------------------------
    subroutine AM_STRETCH_BEND_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_STRETCH_BEND_set_user_bit
!-------------------------------------------------------------
    subroutine AM_STRETCH_BEND_eval(crd, frc, ene, vir)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), vir(3, 3)
        _REAL_, intent(out) :: ene
#include "do_flag.h"

        integer j, k
! initialize
        energy = 0.d0
        do k = 1, 3
            do j = 1, 3
                virial(j, k) = 0.d0
            end do
        end do
        if (do_flag /= PROCEED) then
            return
        end if
! first get the geometry terms
! first the stretch contribution
        call AM_STRETCH_BEND_get_stretchargs(crd, list, &
            bond1_equil_value, bond2_equil_value, arg1, darg1_dcrd)
! next the bend contribution.can use standard angle function
        call AM_STRETCH_BEND_get_bendargs(crd, list, angle_equil_value, arg2, darg2_dcrd)
! finally call the energy force evaluator. energy is given by a constant
!  times arg1 x arg2---so doesn't fit into the usual functable methodology
! which assumes functions of one argument
        call AM_STRETCH_BEND_get_ene_frc(list, force_constant, arg1, arg2, &
            darg1_dcrd, darg2_dcrd, crd, frc)
        ene = energy
        do k = 1, 3
            do j = 1, 3
                vir(j, k) = vir(j, k) + virial(j, k)
            end do
        end do
    end subroutine AM_STRETCH_BEND_eval
!---------------------------------------------------
    subroutine AM_STRETCH_BEND_get_stretchargs(crd, sblist, &
        bond1_equil_value, bond2_equil_value, &
        arg, darg_dcrdijk)
!--------------------------------------------------------
! This routine calculates ij & kj bond function arguments and
! their derivatives with
! respect to atomic positions of atom i,j,k. This gets used
! in a stretch bend energy routine
! NOTE that periodic boundary conditions are not used here
!      i.e. imaging is done on a per molecule basis
!--------------------------------------------------------
! INPUT variables:
!    mstrbend: the size of the stretchbend list
!    crd the atomic coord array
!    sblist: 4 x mstrbend array giving for each strbend the index of the first
!           atom, index of the second atom, index of the thirs,
!           and index into the stretch bend parameter table giving the
!           force constant and equilibrium angle (as well as ideal bond dists)
! OUTPUT variables:
!    arg, array of stretch bend function args
!    darg_dcrdi   derivs of arg wrt crds of atom i
!    darg_dcrdj   derivs of arg wrt crds of atom j
!    darg_dcrdk   derivs of arg wrt crds of atom k
!--------------------------------------------------------
        integer, intent(in) :: sblist(4, *)
        _REAL_, intent(in) :: crd(3, *), bond1_equil_value(*), bond2_equil_value(*)
        _REAL_, intent(out) :: arg(*), darg_dcrdijk(9, *)

        integer :: i, j, k, n, it, ind
        _REAL_ :: xij, yij, zij, xkj, ykj, zkj, bl

        do n = startlist, endlist
            ind = n - startlist + 1
            i = sblist(1, n)
            j = sblist(2, n)
            k = sblist(3, n)
            it = sblist(4, n)
! first the ij bond
            xij = crd(1, i) - crd(1, j)
            yij = crd(2, i) - crd(2, j)
            zij = crd(3, i) - crd(3, j)
            bl = sqrt(xij**2 + yij**2 + zij**2)
            arg(ind) = bl - bond1_equil_value(it)
! differentiate bl to get darg_dcrd of i
            darg_dcrdijk(1, ind) = xij/bl
            darg_dcrdijk(2, ind) = yij/bl
            darg_dcrdijk(3, ind) = zij/bl
! next add the kj bond contributions
            xkj = crd(1, k) - crd(1, j)
            ykj = crd(2, k) - crd(2, j)
            zkj = crd(3, k) - crd(3, j)
            bl = sqrt(xkj**2 + ykj**2 + zkj**2)
            arg(ind) = arg(ind) + bl - bond2_equil_value(it)
            darg_dcrdijk(7, ind) = xkj/bl
            darg_dcrdijk(8, ind) = ykj/bl
            darg_dcrdijk(9, ind) = zkj/bl
! j get the negative of the forces on i,k
            darg_dcrdijk(4, ind) = -(darg_dcrdijk(1, ind) + darg_dcrdijk(7, ind))
            darg_dcrdijk(5, ind) = -(darg_dcrdijk(2, ind) + darg_dcrdijk(8, ind))
            darg_dcrdijk(6, ind) = -(darg_dcrdijk(3, ind) + darg_dcrdijk(9, ind))
        end do
    end subroutine AM_STRETCH_BEND_get_stretchargs
!--------------------------------------------------------------------
    subroutine AM_STRETCH_BEND_get_bendargs(crd, sblist, angle_equil_value, &
        arg, darg_dcrdijk)
!--------------------------------------------------------
! This routine calculates angle function argument and its derivatives with
! respect to atomic positions of atom i,j,k.
!--------------------------------------------------------
! INPUT variables:
!    nangles:  number of angles in list
!    crd the atomic coord array
!    alist: 3 x nangles array giving for each angle the index of the first
!           atom, index of the second atom, index of the thirs,
! OUTPUT variables:
!    arg, array of angle function args
!    darg_dcrdijk   derivs of arg wrt crds of atom i,j,k
!--------------------------------------------------------
        integer, intent(in) :: sblist(4, *)
        _REAL_, intent(in) :: crd(3, *), angle_equil_value(*)
        _REAL_, intent(out) :: arg(*), darg_dcrdijk(9, *)

! LOCAL variables
        integer :: i, j, k, n, ind, it
        _REAL_ :: xij, yij, zij, xkj, ykj, zkj, cosang, &
            dotp, lenij2, lenkj2, lenp, ang, dang_dcosang, &
            dcosang_dxij, dcosang_dyij, dcosang_dzij, &
            dcosang_dxkj, dcosang_dykj, dcosang_dzkj

! note units are degrees not radians...possibly change back later
        do n = startlist, endlist
            ind = n - startlist + 1
            i = sblist(1, n)
            j = sblist(2, n)
            k = sblist(3, n)
            it = sblist(4, n)
            xij = crd(1, i) - crd(1, j)
            yij = crd(2, i) - crd(2, j)
            zij = crd(3, i) - crd(3, j)
            xkj = crd(1, k) - crd(1, j)
            ykj = crd(2, k) - crd(2, j)
            zkj = crd(3, k) - crd(3, j)
! cosine of angle is given by dot product of rij and rkj
!   divided by the product of the lengths of rij and rkj
            dotp = xij*xkj + yij*ykj + zij*zkj
            lenij2 = xij**2 + yij**2 + zij**2
            lenkj2 = xkj**2 + ykj**2 + zkj**2
            lenp = sqrt(lenij2*lenkj2)
            cosang = dotp/lenp
! avoid angle of pi and 0; really you need to use cosangle formulation
! near those points. however due to severe strain you could get bad values
! even for reference angle not near 0 or pi.
            cosang = max(-pt999999, cosang)
            cosang = min(pt999999, cosang)
            ang = radians_to_degrees*acos(cosang)
            dang_dcosang = -radians_to_degrees/sqrt(1.d0 - cosang**2)
! again note units are degrees for now
!       ang = acos(cosang)
!       dang_dcosang = -1.d0 / sqrt(1.d0-cosang**2)
! angle function argument is angle ijk minus reference angle
! reference angle is aparm(2,at)
            arg(ind) = ang - angle_equil_value(it)
! deriv of dotp wrt xij is xkj; deriv of lenp^-1 wrt xij is
!  lenkj^-1 * (-lenij^-2) * (xij/lenij) = -xij / (lenp*lenij2)
! similar for others
            dcosang_dxij = xkj/lenp - (dotp*xij)/(lenp*lenij2)
            dcosang_dyij = ykj/lenp - (dotp*yij)/(lenp*lenij2)
            dcosang_dzij = zkj/lenp - (dotp*zij)/(lenp*lenij2)
            dcosang_dxkj = xij/lenp - (dotp*xkj)/(lenp*lenkj2)
            dcosang_dykj = yij/lenp - (dotp*ykj)/(lenp*lenkj2)
            dcosang_dzkj = zij/lenp - (dotp*zkj)/(lenp*lenkj2)
! now use the chain rule
! first the i crds
            darg_dcrdijk(1, ind) = dang_dcosang*dcosang_dxij
            darg_dcrdijk(2, ind) = dang_dcosang*dcosang_dyij
            darg_dcrdijk(3, ind) = dang_dcosang*dcosang_dzij
! next the k crds
            darg_dcrdijk(7, ind) = dang_dcosang*dcosang_dxkj
            darg_dcrdijk(8, ind) = dang_dcosang*dcosang_dykj
            darg_dcrdijk(9, ind) = dang_dcosang*dcosang_dzkj
! finally the j crds
            darg_dcrdijk(4, ind) = -dang_dcosang*(dcosang_dxij + dcosang_dxkj)
            darg_dcrdijk(5, ind) = -dang_dcosang*(dcosang_dyij + dcosang_dykj)
            darg_dcrdijk(6, ind) = -dang_dcosang*(dcosang_dzij + dcosang_dzkj)
        end do
    end subroutine AM_STRETCH_BEND_get_bendargs
!--------------------------------------------------------------------
    subroutine AM_STRETCH_BEND_get_ene_frc(sblist, force_constant, arg1, arg2, &
        darg1_dcrdijk, darg2_dcrdijk, crd, frc)
        integer, intent(in) :: sblist(4, *)
        _REAL_, intent(in) :: force_constant(*)
        _REAL_, intent(in) :: arg1(*), arg2(*)
        _REAL_, intent(in) :: darg1_dcrdijk(9, *), darg2_dcrdijk(9, *)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *)

        integer  :: i, j, k, n, p, q, it, ind
        _REAL_ :: Frc_K, term1, term2, f(9)

        do n = startlist, endlist
            ind = n - startlist + 1
            i = sblist(1, n)
            j = sblist(2, n)
            k = sblist(3, n)
            it = sblist(4, n)
            Frc_K = degrees_to_radians*force_constant(it)
            energy = energy + Frc_K*arg1(ind)*arg2(ind) ! prod of stretch,bend terms
! apply product rule to get deriv of energy with respect to crds of i,j,k
!  while e.g. darg1_dcrdi holds the derivs of arg1 with respect to crds of i
! recall force is negative of grad
            term1 = Frc_K*arg1(ind)
            term2 = Frc_K*arg2(ind)
            f(1) = term2*darg1_dcrdijk(1, ind) + term1*darg2_dcrdijk(1, ind)
            f(2) = term2*darg1_dcrdijk(2, ind) + term1*darg2_dcrdijk(2, ind)
            f(3) = term2*darg1_dcrdijk(3, ind) + term1*darg2_dcrdijk(3, ind)
            f(4) = term2*darg1_dcrdijk(4, ind) + term1*darg2_dcrdijk(4, ind)
            f(5) = term2*darg1_dcrdijk(5, ind) + term1*darg2_dcrdijk(5, ind)
            f(6) = term2*darg1_dcrdijk(6, ind) + term1*darg2_dcrdijk(6, ind)
            f(7) = term2*darg1_dcrdijk(7, ind) + term1*darg2_dcrdijk(7, ind)
            f(8) = term2*darg1_dcrdijk(8, ind) + term1*darg2_dcrdijk(8, ind)
            f(9) = term2*darg1_dcrdijk(9, ind) + term1*darg2_dcrdijk(9, ind)
            frc(1, i) = frc(1, i) - f(1)
            frc(2, i) = frc(2, i) - f(2)
            frc(3, i) = frc(3, i) - f(3)
            frc(1, j) = frc(1, j) - f(4)
            frc(2, j) = frc(2, j) - f(5)
            frc(3, j) = frc(3, j) - f(6)
            frc(1, k) = frc(1, k) - f(7)
            frc(2, k) = frc(2, k) - f(8)
            frc(3, k) = frc(3, k) - f(9)
! now get virial
            do q = 1, 3
                do p = 1, 3
                    virial(p, q) = virial(p, q) + f(p)*crd(q, i) + f(p + 3)*crd(q, j) + &
                        f(p + 6)*crd(q, k)
                end do
            end do
        end do
    end subroutine AM_STRETCH_BEND_get_ene_frc
!---------------------------------------------------
end module amoeba_stretch_bend
!-----------------------------------------------
module amoeba_torsion_torsion

    use constants, only : pi
    implicit none
    private

    integer, save :: num_list = 0, startlist, endlist
    _REAL_, save :: energy, virial(3, 3)

    integer, allocatable, save :: list(:, :)

    _REAL_, allocatable, save :: arg1(:), arg2(:)
    _REAL_, allocatable, save :: darg1_dcrd(:, :), darg2_dcrd(:, :)
    _REAL_, allocatable, save :: func(:), dfunc_darg1(:), dfunc_darg2(:)
    integer, save :: num_params = 0, do_flag = 0
    type :: angle_angle_functable
        integer :: dim1 = 0, dim2 = 0
        _REAL_, pointer :: angle1(:) => null()
        _REAL_, pointer :: angle2(:) => null()
        _REAL_, pointer :: func(:, :) => null()
        _REAL_, pointer :: dfunc_dangle1(:, :) => null()
        _REAL_, pointer :: dfunc_dangle2(:, :) => null()
        _REAL_, pointer :: d2func_dangle1_dangle2(:, :) => null()
    end type angle_angle_functable
    type(angle_angle_functable), allocatable, save :: torsion_torsion_table(:)

    _REAL_, parameter :: radians_to_degrees = 180.d0/pi
    _REAL_, parameter :: degrees_to_radians = pi/180.d0

    public AM_TOR_TOR_readparm, AM_TOR_TOR_deallocate, &
        AM_TOR_TOR_set_user_bit, AM_TOR_TOR_eval
#ifdef MPI
    public AM_TOR_TOR_bcast
#endif
contains

#ifdef MPI
    subroutine AM_TOR_TOR_bcast
        implicit none
        integer siztask, ierr, dim1, dim2, i

        include 'mpif.h'
#  include "extra.h"
#  include "do_flag.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_list, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_params, 1, MPI_INTEGER, 0, commsander, ierr)

        if (num_list .eq. 0) return

        if (.not. master) then
            allocate (list(6, num_list), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (torsion_torsion_table(num_params), stat=ierr)
            REQUIRE(ierr == 0)

            do i = 1, num_params
                call mpi_bcast(torsion_torsion_table(i)%dim1, 1, MPI_INTEGER, 0, commsander, ierr)
                call mpi_bcast(torsion_torsion_table(i)%dim2, 1, MPI_INTEGER, 0, commsander, ierr)
                dim1 = torsion_torsion_table(i)%dim1
                dim2 = torsion_torsion_table(i)%dim2
                allocate (torsion_torsion_table(i)%func(dim1, dim2))
                allocate (torsion_torsion_table(i)%angle1(dim1))
                allocate (torsion_torsion_table(i)%angle2(dim2))
                allocate (torsion_torsion_table(i)%dfunc_dangle1(dim1, dim2))
                allocate (torsion_torsion_table(i)%dfunc_dangle2(dim1, dim2))
                allocate (torsion_torsion_table(i)%d2func_dangle1_dangle2(dim1, dim2))
            end do
        end if

        call mpi_bcast(list, 6*num_list, MPI_INTEGER, 0, commsander, ierr)

        do i = 1, num_params
            dim1 = torsion_torsion_table(i)%dim1
            dim2 = torsion_torsion_table(i)%dim2
            call mpi_bcast(torsion_torsion_table(i)%func, dim1*dim2, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
            call mpi_bcast(torsion_torsion_table(i)%angle1, dim1, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
            call mpi_bcast(torsion_torsion_table(i)%angle2, dim2, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
            call mpi_bcast(torsion_torsion_table(i)%dfunc_dangle1, dim1*dim2, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
            call mpi_bcast(torsion_torsion_table(i)%dfunc_dangle2, dim1*dim2, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
            call mpi_bcast(torsion_torsion_table(i)%d2func_dangle1_dangle2, dim1*dim2, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        end do

        if (.not. master) then
            call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
            allocate (arg1(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (arg2(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (darg1_dcrd(12, siztask), stat=ierr) !4 atoms to get derivs wrt
            REQUIRE(ierr == 0)
            allocate (darg2_dcrd(12, siztask), stat=ierr) !4 atoms to get derivs wrt
            REQUIRE(ierr == 0)
            allocate (func(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (dfunc_darg1(siztask), stat=ierr)
            REQUIRE(ierr == 0)
            allocate (dfunc_darg2(siztask), stat=ierr)
            REQUIRE(ierr == 0)
        end if
    end subroutine AM_TOR_TOR_bcast
#endif

!-------------------------------------------------------------
    function AM_TOR_TOR_readparm(nf)
        integer :: AM_TOR_TOR_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, n, siztask, ier1, ier2, ier3, ier4, ier5, ier6, ier7
        character(len=2) :: word

        AM_TOR_TOR_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_TORSION_TORSION_', nf, num_list)
        if (num_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        dim1 = 6 ! 5 atoms plus param ptr
        allocate (list(dim1, num_list), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate list; num_list = ', num_list
            call mexit(6, 1)
        end if
        call AMOEBA_read_list_data('AMOEBA_TORSION_TORSION_', nf, dim1, num_list, list)
        call AM_VAL_get_num_params('AMOEBA_TORSION_TORSION_', nf, num_params)
        if (num_params <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        allocate (torsion_torsion_table(num_params), stat=ier)
        if (ier /= 0) then
            write (6, *) 'failed to allocate torsion_torsion_table; num_params = ', &
                num_params
            call mexit(6, 1)
        end if
        do n = 1, num_params
            write (word, '(i2.2)') n
            call AM_TOR_TOR_read_tortor_ftable( &
                'AMOEBA_TORSION_TORSION_TORTOR_TABLE_'//word//'_', &
                nf, torsion_torsion_table(n))
        end do
        call AMOEBA_get_startlist_endlist(num_list, startlist, endlist, siztask)
        allocate (arg1(siztask), stat=ier1)
        allocate (arg2(siztask), stat=ier2)
        allocate (darg1_dcrd(12, siztask), stat=ier3) !4 atoms to get derivs wrt
        allocate (darg2_dcrd(12, siztask), stat=ier4) !4 atoms to get derivs wrt
        allocate (func(siztask), stat=ier5)
        allocate (dfunc_darg1(siztask), stat=ier6)
        allocate (dfunc_darg2(siztask), stat=ier7)
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) .or. &
            (ier5 /= 0) .or. (ier6 /= 0) .or. (ier7 /= 0)) then
            write (6, *) 'problems allocating torsion-torsion scratch space '
            call mexit(6, 1)
        end if
        do_flag = ibset(do_flag, VALID_BIT)
        AM_TOR_TOR_readparm = 1
    end function AM_TOR_TOR_readparm
!-------------------------------------------------------------
    subroutine AM_TOR_TOR_read_tortor_ftable(header, nf, tortor_table)
        character(len=*), intent(in) :: header
        integer, intent(in) :: nf
        type(angle_angle_functable) :: tortor_table

        integer :: iok, ionerr, ier1, ier2, ier3, ier4, ier5, ier6
        character(len=80) :: fmt
        character(len=80) :: fmtin, dtype

        ionerr = 0 !fatal if missing
        fmtin = '(10I8)'
        dtype = header//'DIMS'
        call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
        read (nf, fmt) tortor_table%dim1, tortor_table%dim2
        allocate (tortor_table%angle1(tortor_table%dim1), stat=ier1)
        allocate (tortor_table%angle2(tortor_table%dim2), stat=ier2)
        allocate (tortor_table%func(tortor_table%dim1, tortor_table%dim2), stat=ier3)
        allocate (tortor_table%dfunc_dangle1(tortor_table%dim1, tortor_table%dim2), &
            stat=ier4)
        allocate (tortor_table%dfunc_dangle2(tortor_table%dim1, tortor_table%dim2), &
            stat=ier5)
        allocate (tortor_table%d2func_dangle1_dangle2(tortor_table%dim1, &
            tortor_table%dim2), stat=ier6)
        if ((ier1 /= 0) .or. (ier2 /= 0) .or. (ier3 /= 0) .or. (ier4 /= 0) .or. &
            (ier5 /= 0) .or. (ier6 /= 0)) then
            write (6, *) 'AM_TOR_TOR_read_tortor_ftable:problems allocating table '
            call mexit(6, 1)
        end if
        fmtin = '(5E16.8)'
        dtype = header//'ANGLE1'
        call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
        read (nf, fmt) tortor_table%angle1
        dtype = header//'ANGLE2'
        call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
        read (nf, fmt) tortor_table%angle2
        dtype = header//'FUNC'
        call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
        read (nf, fmt) tortor_table%func
        dtype = header//'DFUNC_DANGLE1'
        call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
        read (nf, fmt) tortor_table%dfunc_dangle1
        dtype = header//'DFUNC_DANGLE2'
        call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
        read (nf, fmt) tortor_table%dfunc_dangle2
        dtype = header//'D2FUNC_DANGLE1_DANGLE2'
        call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
        read (nf, fmt) tortor_table%d2func_dangle1_dangle2

    end subroutine AM_TOR_TOR_read_tortor_ftable
!-----------------------------------------------
    subroutine AM_TOR_TOR_deallocate()
        integer :: n
        if (allocated(list)) deallocate (list)
        if (allocated(arg1)) deallocate (arg1)
        if (allocated(arg2)) deallocate (arg2)
        if (allocated(darg1_dcrd)) deallocate (darg1_dcrd)
        if (allocated(darg2_dcrd)) deallocate (darg2_dcrd)
        if (allocated(func)) deallocate (func)
        if (allocated(dfunc_darg1)) deallocate (dfunc_darg1)
        if (allocated(dfunc_darg2)) deallocate (dfunc_darg2)
        if (allocated(torsion_torsion_table)) then
            do n = 1, num_params
                if (associated(torsion_torsion_table(n)%angle1)) &
                    deallocate (torsion_torsion_table(n)%angle1)
                if (associated(torsion_torsion_table(n)%angle2)) &
                    deallocate (torsion_torsion_table(n)%angle2)
                if (associated(torsion_torsion_table(n)%func)) &
                    deallocate (torsion_torsion_table(n)%func)
                if (associated(torsion_torsion_table(n)%dfunc_dangle1)) &
                    deallocate (torsion_torsion_table(n)%dfunc_dangle1)
                if (associated(torsion_torsion_table(n)%dfunc_dangle2)) &
                    deallocate (torsion_torsion_table(n)%dfunc_dangle2)
                if (associated(torsion_torsion_table(n)%d2func_dangle1_dangle2)) &
                    deallocate (torsion_torsion_table(n)%d2func_dangle1_dangle2)
            end do
            deallocate (torsion_torsion_table)
        end if
    end subroutine AM_TOR_TOR_deallocate
!-------------------------------------------------------------
    subroutine AM_TOR_TOR_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"
        if (do_this == 1) then
            do_flag = ibset(do_flag, USER_BIT)
        else
            do_flag = ibclr(do_flag, USER_BIT)
        end if
    end subroutine AM_TOR_TOR_set_user_bit
!-------------------------------------------------------------
    subroutine AM_TOR_TOR_eval(crd, frc, ene, vir)
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(inout) :: frc(3, *), vir(3, 3)
        _REAL_, intent(out) :: ene
#include "do_flag.h"

        integer j, k
! initialize
        energy = 0.d0
        do k = 1, 3
            do j = 1, 3
                virial(j, k) = 0.d0
            end do
        end do
        if (do_flag /= PROCEED) then
            return
        end if
        call AM_TOR_TOR_get_args(list, crd, arg1, arg2, darg1_dcrd, darg2_dcrd)
        call AM_TOR_TOR_func(list, arg1, arg2, torsion_torsion_table, &
            func, dfunc_darg1, dfunc_darg2)
        call AM_TOR_TOR_get_ene_frc(list, crd, func, dfunc_darg1, dfunc_darg2, &
            darg1_dcrd, darg2_dcrd, frc)
        ene = energy
        do k = 1, 3
            do j = 1, 3
                vir(j, k) = vir(j, k) + virial(j, k)
            end do
        end do
    end subroutine AM_TOR_TOR_eval
!-----------------------------------------------
    subroutine AM_TOR_TOR_get_args(list, crd, arg1, arg2, &
        darg1_dcrdijkl, darg2_dcrdjklm)
!--------------------------------------------------------
! This routine calculates torsion-torsion function arguments and their
!  derivatives with  respect to atomic positions of atom i,j,k,l,m.
!--------------------------------------------------------
! INPUT variables:
!    crd the atomic coord array
!    list: 6 x ntortor array giving for each torsion-torsion
!           the index of the first,second,third,fourth and fifth atoms
!           and the param table pointer
! OUTPUT variables:
!    for each torsion-torsion in list
!    arg1--the torsion angle of i,j,k,l
!    arg2--the torsion angle of j,k,l,m
!    darg1_dcrdijkl   derivs of arg1 wrt crds of atom i,j,k,l
!    darg2_dcrdjklm   derivs of arg2 wrt crds of atom j,k,l,m
!--------------------------------------------------------
        integer, intent(in) :: list(6, *) !5 atoms plus param ptr
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(out) :: arg1(*), arg2(*)
        _REAL_, intent(out) :: darg1_dcrdijkl(12, *), darg2_dcrdjklm(12, *)

! LOCAL variables
        integer i, j, k, l, m, n, p, ind
        _REAL_ :: crd_abcd(12), gradphi_abcd(12), cosphi, sinphi
        _REAL_ :: phi

        do n = startlist, endlist
            ind = n - startlist + 1
            i = list(1, n)
            j = list(2, n)
            k = list(3, n)
            l = list(4, n)
            m = list(5, n)
            do p = 1, 3
                crd_abcd(p) = crd(p, i)
                crd_abcd(p + 3) = crd(p, j)
                crd_abcd(p + 6) = crd(p, k)
                crd_abcd(p + 9) = crd(p, l)
            end do
            call AM_VAL_GEOM_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)
            phi = radians_to_degrees*acos(cosphi)
            if (sinphi >= 0.d0) then
                arg1(ind) = phi
            else
                arg1(ind) = -phi
            end if
            do p = 1, 12
                darg1_dcrdijkl(p, ind) = radians_to_degrees*gradphi_abcd(p)
            end do
            do p = 1, 3
                crd_abcd(p) = crd(p, j)
                crd_abcd(p + 3) = crd(p, k)
                crd_abcd(p + 6) = crd(p, l)
                crd_abcd(p + 9) = crd(p, m)
            end do
            call AM_VAL_GEOM_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)
            phi = radians_to_degrees*acos(cosphi)
            if (sinphi >= 0.d0) then
                arg2(ind) = phi
            else
                arg2(ind) = -phi
            end if
            do p = 1, 12
                darg2_dcrdjklm(p, ind) = radians_to_degrees*gradphi_abcd(p)
            end do
        end do
    end subroutine AM_TOR_TOR_get_args
!-----------------------------------------------
    subroutine AM_TOR_TOR_func(list, arg1, arg2, tortor_table, &
        func, dfunc_darg1, dfunc_darg2)
        integer, intent(in) :: list(6, *) !5 atoms plus param ptr
        _REAL_, intent(in) :: arg1(*), arg2(*)
        type(angle_angle_functable), intent(in) :: tortor_table(*)
        _REAL_, intent(out) :: func(*), dfunc_darg1(*), dfunc_darg2(*)

        integer :: ind, n, it, ind1, ind2
        _REAL_ :: ang1, ang2, ang1_lo, ang1_hi, ang2_lo, ang2_hi
        _REAL_ :: f(4), df_da1(4), df_da2(4), d2f_da1_da2(4), &
            e, de_dang1, de_dang2
        integer :: AM_VAL_real_array_index

        do n = startlist, endlist
            ind = n - startlist + 1
            it = list(6, n)
            ang1 = arg1(ind)
            ang2 = arg2(ind)
            ind1 = AM_VAL_real_array_index(ang1, tortor_table(it)%angle1, &
                tortor_table(it)%dim1)
            ind2 = AM_VAL_real_array_index(ang2, tortor_table(it)%angle1, &
                tortor_table(it)%dim1)
            ang1_lo = tortor_table(it)%angle1(ind1)
            ang1_hi = tortor_table(it)%angle1(ind1 + 1)
            ang2_lo = tortor_table(it)%angle2(ind2)
            ang2_hi = tortor_table(it)%angle2(ind2 + 1)
            ! counter-clockwise order around surrounding table vertices
            f(1) = tortor_table(it)%func(ind1, ind2)
            f(2) = tortor_table(it)%func(ind1 + 1, ind2)
            f(3) = tortor_table(it)%func(ind1 + 1, ind2 + 1)
            f(4) = tortor_table(it)%func(ind1, ind2 + 1)
            df_da1(1) = tortor_table(it)%dfunc_dangle1(ind1, ind2)
            df_da1(2) = tortor_table(it)%dfunc_dangle1(ind1 + 1, ind2)
            df_da1(3) = tortor_table(it)%dfunc_dangle1(ind1 + 1, ind2 + 1)
            df_da1(4) = tortor_table(it)%dfunc_dangle1(ind1, ind2 + 1)
            df_da2(1) = tortor_table(it)%dfunc_dangle2(ind1, ind2)
            df_da2(2) = tortor_table(it)%dfunc_dangle2(ind1 + 1, ind2)
            df_da2(3) = tortor_table(it)%dfunc_dangle2(ind1 + 1, ind2 + 1)
            df_da2(4) = tortor_table(it)%dfunc_dangle2(ind1, ind2 + 1)
            d2f_da1_da2(1) = tortor_table(it)%d2func_dangle1_dangle2(ind1, ind2)
            d2f_da1_da2(2) = tortor_table(it)%d2func_dangle1_dangle2(ind1 + 1, ind2)
            d2f_da1_da2(3) = tortor_table(it)%d2func_dangle1_dangle2(ind1 + 1, ind2 + 1)
            d2f_da1_da2(4) = tortor_table(it)%d2func_dangle1_dangle2(ind1, ind2 + 1)
            call AM_VAL_bcuint1(f, df_da1, df_da2, d2f_da1_da2, ang1_lo, ang1_hi, &
                ang2_lo, ang2_hi, ang1, ang2, e, de_dang1, de_dang2)
            func(ind) = e
            dfunc_darg1(ind) = de_dang1
            dfunc_darg2(ind) = de_dang2
        end do
    end subroutine AM_TOR_TOR_func
!-----------------------------------------------
    subroutine AM_TOR_TOR_get_ene_frc(list, crd, fn, dfn_darg1, dfn_darg2, &
        darg1_dcrdijkl, darg2_dcrdjklm, frc)
        integer, intent(in) :: list(6, *) !5 atoms plus param ptr
        _REAL_, intent(in) :: crd(3, *)
        _REAL_, intent(in) :: fn(*), dfn_darg1(*), dfn_darg2(*)
        _REAL_, intent(in) :: darg1_dcrdijkl(12, *), darg2_dcrdjklm(12, *)
        _REAL_, intent(inout) :: frc(3, *)

        integer ::  i, j, k, l, m, n, p, q, ind
        _REAL_ :: f(12), g(12), term

        do n = startlist, endlist
            ind = n - startlist + 1
            i = list(1, n)
            j = list(2, n)
            k = list(3, n)
            l = list(4, n)
            m = list(5, n)
            energy = energy + fn(ind)
! apply chain rule to get deriv of energy with respect to crds of i,j,k,l,m
! dfn_darg1 holds the deriv of fn with respect to arg1
! dfn_darg2 holds the deriv of fn with respect to arg2
! while darg1_dcrdijkl holds the derivs of arg1 with respect to crds of i,j,k,l
! and darg2_dcrdjklm holds the derivs of arg2 with respect to crds of j,k,l,m
! recall force is negative of grad
            term = dfn_darg1(ind)
            do p = 1, 12
                f(p) = term*darg1_dcrdijkl(p, ind)
            end do
            term = dfn_darg2(ind)
            do p = 1, 12
                g(p) = term*darg2_dcrdjklm(p, ind)
            end do
            do p = 1, 3
                frc(p, i) = frc(p, i) - f(p)
                frc(p, j) = frc(p, j) - f(p + 3) - g(p)
                frc(p, k) = frc(p, k) - f(p + 6) - g(p + 3)
                frc(p, l) = frc(p, l) - f(p + 9) - g(p + 6)
                frc(p, m) = frc(p, m) - g(p + 9)
            end do
! now get virial
            do q = 1, 3
                do p = 1, 3
                    virial(p, q) = virial(p, q) + f(p)*crd(q, i) + f(p + 3)*crd(q, j) + &
                        f(p + 6)*crd(q, k) + f(p + 9)*crd(q, l) + &
                        g(p)*crd(q, j) + g(p + 3)*crd(q, k) + &
                        g(p + 6)*crd(q, l) + g(p + 9)*crd(q, m)
                end do
            end do
        end do
    end subroutine AM_TOR_TOR_get_ene_frc
!-----------------------------------------------
end module amoeba_torsion_torsion
!-----------------------------------------------
!-----------------------------------------------
subroutine AM_VAL_get_num_params(header, nf, num_params)
    implicit none
    character(len=*), intent(in) :: header
    integer, intent(in) :: nf
    integer, intent(out) :: num_params

    integer :: iok, ionerr
    character(len=80) :: fmt
    character(len=80) :: fmtin, ifmt, dtype

    ifmt = '(10I8)'
    dtype = header//'NUM_PARAMS'
    ionerr = 1 ! not fatal if missing
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    if (iok == 0) then !this data type found in prmtop
        read (nf, fmt) num_params
    else !either old style prmtop or data not found
        num_params = 0
    end if
end subroutine AM_VAL_get_num_params
!--------------------------------------------------------------------
subroutine AM_VAL_get_ftab_degree(header, nf, degree)
    implicit none
    character(len=*), intent(in) :: header
    integer, intent(in) :: nf
    integer, intent(out) :: degree

    integer :: iok, ionerr
    character(len=80) :: fmt
    character(len=80) :: fmtin, ifmt, dtype

    ifmt = '(10I8)'
    dtype = header//'FTAB_DEGREE'
    ionerr = 0 ! fatal if missing
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    if (iok == 0) then !this data type found in prmtop
        read (nf, fmt) degree
    else !either old style prmtop or data not found
        degree = 0
    end if
end subroutine AM_VAL_get_ftab_degree
!--------------------------------------------------------------------
subroutine AM_VAL_read_force_constant(header, nf, num_params, force_constant)
    implicit none
    character(len=*), intent(in) :: header
    integer, intent(in) :: nf
    integer, intent(in) :: num_params
    _REAL_, intent(out) :: force_constant(num_params)

    integer :: iok, ionerr, k
    character(len=80) :: fmt
    character(len=80) :: fmtin, dtype

    ionerr = 0 !fatal if missing
    fmtin = '(5E16.8)'
    dtype = header//'FORCE_CONSTANT'
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    read (nf, fmt) (force_constant(k), k=1, num_params)
end subroutine AM_VAL_read_force_constant
!----------------------------------------------------------
subroutine AM_VAL_read_equil_value(header, nf, num_params, equil_value)
    implicit none
    character(len=*), intent(in) :: header
    integer, intent(in) :: nf
    integer, intent(in) :: num_params
    _REAL_, intent(out) :: equil_value(num_params)

    integer :: iok, ionerr, k
    character(len=80) :: fmt
    character(len=80) :: fmtin, dtype

    ionerr = 0 !fatal if missing
    fmtin = '(5E16.8)'
    dtype = header//'EQUIL_VALUE'
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    read (nf, fmt) (equil_value(k), k=1, num_params)
end subroutine AM_VAL_read_equil_value
!----------------------------------------------------------
subroutine AM_VAL_read_periodicity(header, nf, num_params, periodicity)
    implicit none
    character(len=*), intent(in) :: header
    integer, intent(in) :: nf
    integer, intent(in) :: num_params
    _REAL_, intent(out) :: periodicity(num_params)

    integer :: iok, ionerr, k
    character(len=80) :: fmt
    character(len=80) :: fmtin, dtype

    ionerr = 0 !fatal if missing
    fmtin = '(5E16.8)'
    dtype = header//'PERIODICITY'
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    read (nf, fmt) (periodicity(k), k=1, num_params)
end subroutine AM_VAL_read_periodicity
!----------------------------------------------------------
subroutine AM_VAL_read_phase(header, nf, num_params, phase)
    implicit none
    character(len=*), intent(in) :: header
    integer, intent(in) :: nf
    integer, intent(in) :: num_params
    _REAL_, intent(out) :: phase(num_params)

    integer :: iok, ionerr, k
    character(len=80) :: fmt
    character(len=80) :: fmtin, dtype

    ionerr = 0 !fatal if missing
    fmtin = '(5E16.8)'
    dtype = header//'PHASE'
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    read (nf, fmt) (phase(k), k=1, num_params)
end subroutine AM_VAL_read_phase
!----------------------------------------------------------
subroutine AM_VAL_read_ftable_coeffs(header, nf, degree, coeff)
    implicit none
    character(len=*), intent(in) :: header
    integer, intent(in) :: nf
    integer, intent(in) :: degree
    _REAL_, intent(out) :: coeff(0:degree)

    integer :: iok, ionerr, j
    character(len=80) :: fmt
    character(len=80) :: fmtin, dtype

    ionerr = 0 !fatal if missing
    fmtin = '(5E16.8)'
    dtype = header//'FTAB_COEFFS'
    call nxtsec(nf, 6, ionerr, fmtin, dtype, fmt, iok)
    read (nf, fmt) (coeff(j), j=0, degree)
end subroutine AM_VAL_read_ftable_coeffs
!----------------------------------------------------------
function AM_VAL_real_array_index(val, array, num)
    implicit none
    integer :: AM_VAL_real_array_index
    integer, intent(in) :: num
    _REAL_, intent(in) :: val, array(num)

    integer indhi, indlo, ind
! bisection search for ind just before val (i.e. array(ind)<val<array(ind+1) )
! assume array is ordered in increasing order
    indlo = 1
    indhi = num
    do while (indhi - indlo > 1)
        ind = (indhi + indlo)/2
        if (array(ind) > val) then
            indhi = ind
        else
            indlo = ind
        end if
    end do
    AM_VAL_real_array_index = indlo
    return
end function AM_VAL_real_array_index
!----------------------------------------------------------
subroutine AM_VAL_FTAB_eval_f_df(nlist, degree, coeff, arg, func, dfunc_darg)
    implicit none
    integer, intent(in) :: nlist
    integer, intent(in) :: degree
    _REAL_, intent(in) :: coeff(0:degree)
    _REAL_, intent(in) :: arg(nlist)
    _REAL_, intent(out) :: func(nlist), dfunc_darg(nlist)

    integer :: n, deg
    _REAL_ :: dx

! setup in case of smaller degrees
    if (degree == 2) then
        do n = 1, nlist
            dx = arg(n)
            func(n) = coeff(0) + dx*(coeff(1) + dx*coeff(2))
            dfunc_darg(n) = coeff(1) + 2.d0*dx*coeff(2)
        end do
    elseif (degree == 3) then
        do n = 1, nlist
            dx = arg(n)
            func(n) = coeff(0) + dx*(coeff(1) + dx*(coeff(2) + dx*coeff(3)))
            dfunc_darg(n) = coeff(1) + dx*(2.d0*coeff(2) + dx*3.d0*coeff(3))
        end do
    elseif (degree == 4) then
        do n = 1, nlist
            dx = arg(n)
            func(n) = coeff(0) + dx*(coeff(1) + dx*(coeff(2) + dx*(coeff(3) + &
                dx*coeff(4))))
            dfunc_darg(n) = coeff(1) + dx*(2.d0*coeff(2) + dx*(3.d0*coeff(3) + &
                dx*4.d0*coeff(4)))
        end do
    else
        do n = 1, nlist
            dx = arg(n)
            func(n) = coeff(degree)
            dfunc_darg(n) = 0.d0
            do deg = degree - 1, 0, -1
                dfunc_darg(n) = func(n) + dx*dfunc_darg(n)
                func(n) = coeff(deg) + dx*func(n)
            end do
        end do
    end if
end subroutine AM_VAL_FTAB_eval_f_df
!-------------------------------------------------------------------------------
subroutine AM_VAL_GEOM_torsion(crd_abcd, gradphi_abcd, cosphi, sinphi)
    implicit none
    _REAL_, intent(in) :: crd_abcd(12)
    _REAL_, intent(out) :: gradphi_abcd(12), cosphi, sinphi
! given coords of points a,b,c,d this routine calculates
! cosine and sine of torsion phi
! as well as gradient of phi with respect to coords of a,b,c,d
! units are in radians

    _REAL_ :: rab(3), rcb(3), rdc(3), ucb(3), upab(3), updc(3), rcross(3), &
        upabc(3), upbcd(3), sizcb, dotp_ab_cb, sizpab, dotp_dc_cb, &
        sizpdc, S(3), dot
    integer :: m
    do m = 1, 3
        rab(m) = crd_abcd(m) - crd_abcd(m + 3)
        rcb(m) = crd_abcd(m + 6) - crd_abcd(m + 3)
        rdc(m) = crd_abcd(m + 9) - crd_abcd(m + 6)
    end do
    sizcb = sqrt(rcb(1)*rcb(1) + rcb(2)*rcb(2) + rcb(3)*rcb(3))
    ucb(1) = rcb(1)/sizcb
    ucb(2) = rcb(2)/sizcb
    ucb(3) = rcb(3)/sizcb
    dotp_ab_cb = rab(1)*ucb(1) + rab(2)*ucb(2) + rab(3)*ucb(3)
! upab is unit vector along component rab perp to ucb
    dot = rab(1)*ucb(1) + rab(2)*ucb(2) + rab(3)*ucb(3)
    upab(1) = rab(1) - dot*ucb(1)
    upab(2) = rab(2) - dot*ucb(2)
    upab(3) = rab(3) - dot*ucb(3)
    sizpab = sqrt(upab(1)*upab(1) + upab(2)*upab(2) + upab(3)*upab(3))
    upab(1) = upab(1)/sizpab
    upab(2) = upab(2)/sizpab
    upab(3) = upab(3)/sizpab
    dotp_dc_cb = rdc(1)*ucb(1) + rdc(2)*ucb(2) + rdc(3)*ucb(3)
! updc is unit vector along component rdc perp to ucb
    dot = rdc(1)*ucb(1) + rdc(2)*ucb(2) + rdc(3)*ucb(3)
    updc(1) = rdc(1) - dot*ucb(1)
    updc(2) = rdc(2) - dot*ucb(2)
    updc(3) = rdc(3) - dot*ucb(3)
    sizpdc = sqrt(updc(1)*updc(1) + updc(2)*updc(2) + updc(3)*updc(3))
    updc(1) = updc(1)/sizpdc
    updc(2) = updc(2)/sizpdc
    updc(3) = updc(3)/sizpdc
! cosine of phi is given by dot product of upab and updc
    cosphi = upab(1)*updc(1) + upab(2)*updc(2) + upab(3)*updc(3)
! sine of phi is given by dot product of ucb and upab x updc
    rcross(1) = upab(2)*updc(3) - upab(3)*updc(2)
    rcross(2) = upab(3)*updc(1) - upab(1)*updc(3)
    rcross(3) = upab(1)*updc(2) - upab(2)*updc(1)
    sinphi = rcross(1)*ucb(1) + rcross(2)*ucb(2) + rcross(3)*ucb(3)
! gradient of phi wrt ra is perp to abc plane---movement of ra by dr perp
! to abc plane results in dphi of dr/sizpab
! perp to abc given by upab x ucb  (these are orthogonal unit vectors)
    upabc(1) = upab(2)*ucb(3) - upab(3)*ucb(2)
    upabc(2) = upab(3)*ucb(1) - upab(1)*ucb(3)
    upabc(3) = upab(1)*ucb(2) - upab(2)*ucb(1)
! grad of phi wrt rd is perp to bcd plane--calc sim to grad phi wrt ra
! perp given by updc x ucb or ucb x updc
    upbcd(1) = ucb(2)*updc(3) - ucb(3)*updc(2)
    upbcd(2) = ucb(3)*updc(1) - ucb(1)*updc(3)
    upbcd(3) = ucb(1)*updc(2) - ucb(2)*updc(1)
! now have enough for gradphi for a and d
    do m = 1, 3
        gradphi_abcd(m) = upabc(m)/sizpab
        gradphi_abcd(9 + m) = upbcd(m)/sizpdc
    end do
! following chap 5 of thesis of Bekker we have grad phi wrt b = -grad phi wrt a
! plus some vec S and rad phi wrt c = -grad phi wrt d - S
! S is perp to rcb; using simple torque rule and identity for
! triple cross product he derives S (eqn 5.20)
    do m = 1, 3
        S(m) = (dotp_ab_cb/sizcb)*gradphi_abcd(m) + &
            (dotp_dc_cb/sizcb)*gradphi_abcd(m + 9)
        gradphi_abcd(m + 3) = S(m) - gradphi_abcd(m)
        gradphi_abcd(m + 6) = -S(m) - gradphi_abcd(m + 9)
    end do
end subroutine AM_VAL_GEOM_torsion
!-------------------------------------------------------
_REAL_ function VEC3D_unitperpto_unitvec(v, u, w)
    implicit none
! removes component of v along unit vector u --returns length of new w
! normalizes resulting w
    _REAL_ v(3), u(3), w(3), siz, dot
    dot = v(1)*u(1) + v(2)*u(2) + v(3)*u(3)
    w(1) = v(1) - dot*u(1)
    w(2) = v(2) - dot*u(2)
    w(3) = v(3) - dot*u(3)
    siz = sqrt(w(1)*w(1) + w(2)*w(2) + w(3)*w(3))
    w(1) = w(1)/siz
    w(2) = w(2)/siz
    w(3) = w(3)/siz
    VEC3D_unitperpto_unitvec = siz
    return
end function VEC3D_unitperpto_unitvec
!-------------------------------------------------------
subroutine AM_VAL_VEC3D_get_perp_to_vecs(v1, v2, p, dp_dv1_p, dp_dv2_p)
    implicit none
    _REAL_ :: v1(3), v2(3), p(3), dp_dv1_p(3), dp_dv2_p(3)
! given two vectors v1,v2, with associated unit vectors u1,u2, get the unit
! cross product p = v1 x v2 / |v1 x v2|
! a change dv1 in v1 within v1v2 plane has no effect on p
! only a change dv1_p along p can change p---to produce this change without
! changing v2 need to rotate about v2; this causes change in p along u2p
! where u2p is u2 x p; Change dv1_p is equiv to rotation of size
!  dv1_p / v1_dot_u2p about axis u2 = v2/|v2| , where v1_dot_u2p is dot
!  product of v1 with u2p; change in p is therefore
! dp = (-dv1_p/v1_dot_u2p)*u2p thus
! dp_dv1_p = -u2p / v1_dot_u2p
    _REAL_ :: siz, u1(3), u2(3), u1p(3), u2p(3), dotp

    siz = v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)
    u1(1) = v1(1)/siz
    u1(2) = v1(2)/siz
    u1(3) = v1(3)/siz
    siz = v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3)
    u2(1) = v2(1)/siz
    u2(2) = v2(2)/siz
    u2(3) = v2(3)/siz
! p = u1 x u2 normed since u1 not perp to u2 necessarily
    p(1) = u1(2)*u2(3) - u1(3)*u2(2)
    p(2) = u1(3)*u2(1) - u1(1)*u2(3)
    p(3) = u1(1)*u2(2) - u1(2)*u2(1)
    siz = sqrt(p(1)*p(1) + p(2)*p(2) + p(3)*p(3))
    p(1) = p(1)/siz
    p(2) = p(2)/siz
    p(3) = p(3)/siz
! u1p is unit vec perp to u1 and p : p x u1
    u1p(1) = p(2)*u1(3) - p(3)*u1(2)
    u1p(2) = p(3)*u1(1) - p(1)*u1(3)
    u1p(3) = p(1)*u1(2) - p(2)*u1(1)
! u2p is unit vec perp to u2 and p : u2 x p
    u2p(1) = u2(2)*p(3) - u2(3)*p(2)
    u2p(2) = u2(3)*p(1) - u2(1)*p(3)
    u2p(3) = u2(1)*p(2) - u2(2)*p(1)
! change in v1 along p gives dp along -u2p
    dotp = v1(1)*u2p(1) + v1(2)*u2p(2) + v1(3)*u2p(3)
    dp_dv1_p(1) = -u2p(1)/dotp
    dp_dv1_p(2) = -u2p(2)/dotp
    dp_dv1_p(3) = -u2p(3)/dotp
! change in v2 along p gives dp along -u1p
    dotp = v2(1)*u1p(1) + v2(2)*u1p(2) + v2(3)*u1p(3)
    dp_dv2_p(1) = -u1p(1)/dotp
    dp_dv2_p(2) = -u1p(2)/dotp
    dp_dv2_p(3) = -u1p(3)/dotp
    return
end subroutine AM_VAL_VEC3D_get_perp_to_vecs
!-------------------------------------------------------------------------------
!-------------------------------------------------------------
!     #################################################################
!     ##                                                             ##
!     ##  subroutine bcuint1  --  bicubic interpolation of gradient  ##
!     ##                                                             ##
!     #################################################################
!
!
!     "bcuint1" performs a bicubic interpolation of the function
!     value and gradient along the directions of a 2D spline grid
!
!     contributed by Jay William Ponder
!
!     literature reference:
!
!     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
!     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
!     University Press, 1992, Section 3.6
!
!-------------------------------------------------------------
subroutine AM_VAL_bcuint1(y, y1, y2, y12, x1l, x1u, x2l, x2u, x1, x2, ansy, ansy1, ansy2)
    _REAL_, intent(in) ::  y(4), y1(4), y2(4), y12(4)
    _REAL_, intent(in) :: x1, x1l, x1u
    _REAL_, intent(in) :: x2, x2l, x2u
    _REAL_, intent(out) :: ansy, ansy1, ansy2

    integer i
    _REAL_ t, u, c(4, 4)

! get coefficients, then perform bicubic interpolation
    call AM_VAL_bcucof(y, y1, y2, y12, x1u - x1l, x2u - x2l, c)
    t = (x1 - x1l)/(x1u - x1l)
    u = (x2 - x2l)/(x2u - x2l)
    ansy = 0.0d0
    ansy1 = 0.0d0
    ansy2 = 0.0d0
    do i = 4, 1, -1
        ansy = t*ansy + ((c(i, 4)*u + c(i, 3))*u + c(i, 2))*u + c(i, 1)
        ansy1 = u*ansy1 + (3.0d0*c(4, i)*t + 2.0d0*c(3, i))*t + c(2, i)
        ansy2 = t*ansy2 + (3.0d0*c(i, 4)*u + 2.0d0*c(i, 3))*u + c(i, 2)
    end do
    ansy1 = ansy1/(x1u - x1l)
    ansy2 = ansy2/(x2u - x2l)
end subroutine AM_VAL_bcuint1
!-------------------------------------------------------------
!     "bcucof" determines the coefficient matrix needed for bicubic
!     interpolation of a function, gradients and cross derivatives
!
subroutine AM_VAL_bcucof(y, y1, y2, y12, d1, d2, c)
    _REAL_, intent(in) ::  y(4), y1(4), y2(4), y12(4)
    _REAL_, intent(in) :: d1, d2
    _REAL_, intent(out) :: c(4, 4)

    _REAL_ xx, d1d2
    _REAL_ x(16), cl(16)
    _REAL_ wt(16, 16)
    integer i, j, k
    save wt
    data wt/1.0d0, 0.0d0, -3.0d0, 2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        -3.0d0, 0.0d0, 9.0d0, -6.0d0, 2.0d0, 0.0d0, -6.0d0, 4.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        3.0d0, 0.0d0, -9.0d0, 6.0d0, -2.0d0, 0.0d0, 6.0d0, -4.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 9.0d0, -6.0d0, 0.0d0, 0.0d0, -6.0d0, 4.0d0, &
        0.0d0, 0.0d0, 3.0d0, -2.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, -9.0d0, 6.0d0, 0.0d0, 0.0d0, 6.0d0, -4.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, -3.0d0, 2.0d0, &
        -2.0d0, 0.0d0, 6.0d0, -4.0d0, 1.0d0, 0.0d0, -3.0d0, 2.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        -1.0d0, 0.0d0, 3.0d0, -2.0d0, 1.0d0, 0.0d0, -3.0d0, 2.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, -3.0d0, 2.0d0, 0.0d0, 0.0d0, 3.0d0, -2.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 3.0d0, -2.0d0, &
        0.0d0, 0.0d0, -6.0d0, 4.0d0, 0.0d0, 0.0d0, 3.0d0, -2.0d0, &
        0.0d0, 1.0d0, -2.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, -3.0d0, 6.0d0, -3.0d0, 0.0d0, 2.0d0, -4.0d0, 2.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 3.0d0, -6.0d0, 3.0d0, 0.0d0, -2.0d0, 4.0d0, -2.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, -3.0d0, 3.0d0, 0.0d0, 0.0d0, 2.0d0, -2.0d0, &
        0.0d0, 0.0d0, -1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 3.0d0, -3.0d0, 0.0d0, 0.0d0, -2.0d0, 2.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, -2.0d0, 1.0d0, &
        0.0d0, -2.0d0, 4.0d0, -2.0d0, 0.0d0, 1.0d0, -2.0d0, 1.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, -1.0d0, 2.0d0, -1.0d0, 0.0d0, 1.0d0, -2.0d0, 1.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, 1.0d0, -1.0d0, 0.0d0, 0.0d0, -1.0d0, 1.0d0, &
        0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, -1.0d0, 1.0d0, &
        0.0d0, 0.0d0, 2.0d0, -2.0d0, 0.0d0, 0.0d0, -1.0d0, 1.0d0/
!
!     pack a temporary vector of corner values
!
    d1d2 = d1*d2
    do i = 1, 4
        x(i) = y(i)
        x(i + 4) = y1(i)*d1
        x(i + 8) = y2(i)*d2
        x(i + 12) = y12(i)*d1d2
    end do
!
!     matrix multiply by the stored weight table
!
    do i = 1, 16
        xx = 0.0d0
        do k = 1, 16
            xx = xx + wt(i, k)*x(k)
        end do
        cl(i) = xx
    end do
!
!     unpack the result into the coefficient table
!
    j = 0
    do i = 1, 4
        do k = 1, 4
            j = j + 1
            c(i, k) = cl(j)
        end do
    end do
end subroutine AM_VAL_bcucof
