#include "dprec.fh"
#include "assert.fh"

!------------------------------------------------------------
module amoeba_adjust
    implicit none
    private

#  include "amoeba_mpole_index.h"

    integer, save :: do_flag
    integer, save, allocatable :: adjust_list(:, :)
    integer, save :: num_adjust_list
    _REAL_, save, allocatable :: dipole_dipole_tensor(:, :)
! adjust list has 3rd component describing the relationship of the pair
! to each other: 9 possible values. Value is 1,2,3,4 if they are
! 1-2,1-3,1-4 or 1-5 pair; 5,6,7,8 if they are 1-2,1-3,1-4 or 1-5 pair
! within a polar group, or 9 if they are in a polar group but beyond 1-5
    _REAL_, save, dimension(9) :: vdw_weight, mpole_weight, direct_weight, &
        polar_weight, mutual_weight
    public AM_ADJUST_readparm, AM_ADJUST_permfield, AM_ADJUST_set_user_bit, &
        AM_ADJUST_dip_dip_fields, AM_ADJUST_ene_frc
#ifdef MPI
    public AM_ADJUST_bcast
#endif
contains
!-------------------------------------------------------
    subroutine AM_ADJUST_set_user_bit(do_this)
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
            write (6, *) 'AM_ADJUST_set_user_bit: bad value of user do_this'
            call mexit(6, 1)
        end if

    end subroutine AM_ADJUST_set_user_bit
!-------------------------------------------------------
#ifdef MPI
    subroutine AM_ADJUST_bcast
        implicit none
        integer ierr

        include 'mpif.h'
#  include "extra.h"
#  include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(num_adjust_list, 1, MPI_INTEGER, 0, commsander, ierr)

        if (.not. master) then
            allocate (adjust_list(3, num_adjust_list), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        call mpi_bcast(adjust_list, 3*num_adjust_list, MPI_INTEGER, 0, commsander, ierr)
        call mpi_bcast(vdw_weight, 9, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(mpole_weight, 9, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(direct_weight, 9, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(polar_weight, 9, MPI_DOUBLE_PRECISION, 0, commsander, ierr)
        call mpi_bcast(mutual_weight, 9, MPI_DOUBLE_PRECISION, 0, commsander, ierr)

        if (.not. master) then
            allocate (dipole_dipole_tensor(6, num_adjust_list), stat=ierr)
            REQUIRE(ierr == 0)
        end if
    end subroutine AM_ADJUST_bcast
#endif

    function AM_ADJUST_readparm(nf)
        integer :: AM_ADJUST_readparm
        integer, intent(in) :: nf
#include "do_flag.h"
        integer :: ier, dim1, num

        AM_ADJUST_readparm = 0
        call AMOEBA_get_numlist('AMOEBA_ADJUST_', nf, num_adjust_list)
        if (num_adjust_list <= 0) then
            do_flag = ibclr(do_flag, VALID_BIT)
            return
        end if
        !allocate
        allocate (adjust_list(3, num_adjust_list), stat=ier)
        REQUIRE(ier == 0)
        dim1 = 3
        call AMOEBA_read_list_data('AMOEBA_ADJUST_', nf, dim1, num_adjust_list, &
            adjust_list)
        ! read in weight lists
        dim1 = 1
        num = 9
        call AMOEBA_read_real_list_data('AMOEBA_ADJUST_VDW_WEIGHTS_', nf, &
            dim1, num, vdw_weight)
        call AMOEBA_read_real_list_data('AMOEBA_ADJUST_MPOLE_WEIGHTS_', nf, &
            dim1, num, mpole_weight)
        call AMOEBA_read_real_list_data('AMOEBA_ADJUST_DIRECT_WEIGHTS_', nf, &
            dim1, num, direct_weight)
        call AMOEBA_read_real_list_data('AMOEBA_ADJUST_POLAR_WEIGHTS_', nf, &
            dim1, num, polar_weight)
        call AMOEBA_read_real_list_data('AMOEBA_ADJUST_MUTUAL_WEIGHTS_', nf, &
            dim1, num, mutual_weight)
        ! allocate the dipole_dipole_tensors
        allocate (dipole_dipole_tensor(6, num_adjust_list), stat=ier)
        REQUIRE(ier == 0)
        AM_ADJUST_readparm = 1
        do_flag = ibset(do_flag, VALID_BIT)
    end function AM_ADJUST_readparm
!-------------------------------------------------------
    subroutine AM_ADJUST_permfield(crd, x, direct_gradphi, polar_gradphi)

        use amoeba_multipoles, only : global_multipole
        use amoeba_induced, only : sq_polinv, is_polarizable
        use amoeba_mdin, only : thole_expon_coeff

        _REAL_, intent(in) :: crd(3, *), x(*)
        _REAL_, intent(inout) :: direct_gradphi(3, *), polar_gradphi(3, *)

#  include "def_time.h"
#  include "ew_erfc_spline.h"
#  include "ew_pme_recip.h"

        call timer_start(TIME_ADJ)
        call AM_ADJUST_calc_permfield(crd, num_adjust_list, adjust_list, &
            ew_coeff, eedtbdns, x(leed_cub), x(leed_lin), &
            ee_type, eedmeth, dxdr, thole_expon_coeff, &
            sq_polinv, is_polarizable, global_multipole, &
            direct_weight, polar_weight, mutual_weight, &
            direct_gradphi, polar_gradphi, &
            dipole_dipole_tensor)
        call timer_stop(TIME_ADJ)
    end subroutine AM_ADJUST_permfield
!-------------------------------------------------------
    subroutine AM_ADJUST_calc_permfield(crd, num_adjust_list, adjust_list, &
        ewaldcof, eedtbdns, eed_cub, eed_lin, &
        ee_type, eedmeth, dxdr, thole_expon_coeff, &
        sq_polinv, is_polarizable, global_multipole, &
        direct_weight, polar_weight, mutual_weight, &
        direct_gradphi, polar_gradphi, &
        dipole_dipole_tensor)
        use constants, only : zero, third, half, one, two, three, five
        _REAL_, intent(in) :: crd(3, *)
        integer, intent(in) :: num_adjust_list, adjust_list(3, *)
        _REAL_, intent(in) :: ewaldcof, eedtbdns, eed_cub(4, *), eed_lin(2, *)
        integer, intent(in) :: ee_type, eedmeth
        _REAL_, intent(in) :: dxdr, thole_expon_coeff, sq_polinv(*)
        logical, intent(in) :: is_polarizable(*)
        _REAL_, intent(in) :: global_multipole(10, *)
        _REAL_, intent(in) :: direct_weight(*), polar_weight(*), mutual_weight(*)
        _REAL_, intent(inout) :: direct_gradphi(3, *), polar_gradphi(3, *)
        _REAL_, intent(out) :: dipole_dipole_tensor(6, *)

        _REAL_ :: delx, dely, delz, delr2, delr, delr2inv, x, dx, switch, d_switch_dx, xx
        _REAL_ :: D(3), A(0:3), B(0:3), BD(3), fac, fact, del, gphi_i(3), gphi_j(3)
        _REAL_  :: expon, expo, clam3, clam5, clam7, delr3inv, delr5inv, delr7inv
        _REAL_  :: Rn(1), Rn_1(4), Rn_2(10), Rn_3(20)

        integer :: i, j, k, n, n_adj, ind
#  include "do_flag.h"

        if (iand(do_flag, PROCEED_INDUCE) /= PROCEED_INDUCE) return

        fac = two*ewaldcof*ewaldcof
        del = one/eedtbdns
        do n_adj = 1, num_adjust_list
            i = adjust_list(1, n_adj)
            j = adjust_list(2, n_adj)
            if (is_polarizable(i) .or. is_polarizable(j)) then
                k = adjust_list(3, n_adj)
                delx = crd(1, j) - crd(1, i)
                dely = crd(2, j) - crd(2, i)
                delz = crd(3, j) - crd(3, i)
                delr2 = delx*delx + dely*dely + delz*delz
                delr = sqrt(delr2)
                delr2inv = one/delr2
                x = dxdr*delr
                if (eedmeth == 1) then
                    !           -- cubic spline on switch
                    ind = eedtbdns*x + 1
                    dx = x - (ind - one)*del
                    switch = eed_cub(1, ind) + dx*(eed_cub(2, ind) + &
                        dx*(eed_cub(3, ind) + dx*eed_cub(4, ind)*third)*half)
                    d_switch_dx = eed_cub(2, ind) + dx*(eed_cub(3, ind) + &
                        dx*eed_cub(4, ind)*half)
                else if (eedmeth == 2) then
                    !           ---linear lookup on switch, deriv
                    xx = eedtbdns*x + 1
                    ind = xx
                    dx = xx - ind
                    switch = (one - dx)*eed_lin(1, ind) + dx*eed_lin(1, ind + 1)
                    d_switch_dx = (one - dx)*eed_lin(2, ind) + dx*eed_lin(2, ind + 1)
                else if (eedmeth == 3) then
                    !           ---direct function call:
                    call get_ee_func(x, switch, d_switch_dx, ee_type)
                else if (eedmeth == 4) then
                    write (6, *) 'eedmeth is 4 inside AM_ADJUST_calc_permfield!!'
                    call mexit(6, 1)
                end if
                !------------------------------------------------------------
                ! McMurchie-Davidson recursion holds for any smooth function of r
                ! that is, to get the higher order derivs wrt x,y,z of g(r)
                ! define R(0,0,0,0) = g(r)
                ! next  R(0,0,0,n+1) = -(1/r)d/dr R(0,0,0,n)
                ! then denote  R(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
                ! quantities of interest obtained by setting n=0
                ! McMurchie-Davidson says that
                !  R(t+1,u,v,n) = t*R(t-1,u,v,n+1) + x*R(t,u,v,n)
                !  R(t,u+1,v,n) = u*R(t,u-1,v,n+1) + y*R(t,u,v,n)
                !  R(t,u,v+1,n) = v*R(t,u,v-1,n+1) + z*R(t,u,v,n)
                ! use similar data structures as in direct sum
                !-------------------------------------------------------------
                ! calc the contributions for cancelling reciprocal sum for ij
                ! -erf = erfc - 1.0 ---first get the erfc part as in direct sum
                B(0) = switch*delr*delr2inv
                fact = d_switch_dx*dxdr
                B(1) = (B(0) - fact)*delr2inv
                fact = fac*fact
                B(2) = (three*B(1) - fact)*delr2inv
                fact = fac*fact
                B(3) = (five*B(2) - fact)*delr2inv
                ! damping factors
                delr3inv = delr2inv/delr
                delr5inv = delr3inv*delr2inv
                delr7inv = delr5inv*delr2inv
                expon = thole_expon_coeff*delr2*delr*sq_polinv(i)*sq_polinv(j)
                expo = exp(-expon)
                ! clam3 = 1.d0-lam3, clam5 = 1.d0-lam5 etc. where
                ! lam is from ponder's paper
                clam3 = expo
                clam5 = (1.d0 + expon)*expo
                clam7 = (1.d0 + expon + 0.6d0*expon**2)*expo
                BD(1) = clam3*delr3inv
                BD(2) = 3.d0*clam5*delr5inv
                BD(3) = 15.d0*clam7*delr7inv
                ! get the 1/r part
                ! note that damped coulomb part given by A(j) - BD(j) while
                ! the erf part is given by A(j) - B(j)
                A(0) = one*delr*delr2inv
                A(1) = A(0)*delr2inv
                A(2) = three*A(1)*delr2inv
                A(3) = five*A(2)*delr2inv
                ! first do the mutual field to store dipole-dipole tensor
                ! smaller recursion than for permanent field
                D(1) = -(mutual_weight(k)*(A(1) - BD(1)) + (B(1) - A(1)))
                D(2) = mutual_weight(k)*(A(2) - BD(2)) + (B(2) - A(2))
                n = 2
                Rn(Ind_000) = D(n)
                Rn_1(Ind_000) = D(n - 1)
                Rn_1(Ind_100) = delx*Rn(Ind_000)
                Rn_1(Ind_010) = dely*Rn(Ind_000)
                Rn_1(Ind_001) = delz*Rn(Ind_000)
                Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
                Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
                Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
                Rn_2(Ind_110) = delx*Rn_1(Ind_010)
                Rn_2(Ind_101) = delx*Rn_1(Ind_001)
                Rn_2(Ind_011) = dely*Rn_1(Ind_001)
                ! store dipole_dipole tensor
                dipole_dipole_tensor(1, n_adj) = Rn_2(Ind_200)
                dipole_dipole_tensor(2, n_adj) = Rn_2(Ind_110)
                dipole_dipole_tensor(3, n_adj) = Rn_2(Ind_101)
                dipole_dipole_tensor(4, n_adj) = Rn_2(Ind_020)
                dipole_dipole_tensor(5, n_adj) = Rn_2(Ind_011)
                dipole_dipole_tensor(6, n_adj) = Rn_2(Ind_002)
                ! next do the direct field
                D(1) = -(direct_weight(k)*(A(1) - BD(1)) + (B(1) - A(1)))
                D(2) = direct_weight(k)*(A(2) - BD(2)) + (B(2) - A(2))
                D(3) = -(direct_weight(k)*(A(3) - BD(3)) + (B(3) - A(3)))
                n = 3
                Rn(Ind_000) = D(n)
                Rn_1(Ind_000) = D(n - 1)
                Rn_1(Ind_100) = delx*Rn(Ind_000)
                Rn_1(Ind_010) = dely*Rn(Ind_000)
                Rn_1(Ind_001) = delz*Rn(Ind_000)
                Rn_2(Ind_000) = D(n - 2)
                Rn_2(Ind_100) = delx*Rn_1(Ind_000)
                Rn_2(Ind_010) = dely*Rn_1(Ind_000)
                Rn_2(Ind_001) = delz*Rn_1(Ind_000)
                Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
                Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
                Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
                Rn_2(Ind_110) = delx*Rn_1(Ind_010)
                Rn_2(Ind_101) = delx*Rn_1(Ind_001)
                Rn_2(Ind_011) = dely*Rn_1(Ind_001)
                Rn_3(Ind_100) = delx*Rn_2(Ind_000)
                Rn_3(Ind_010) = dely*Rn_2(Ind_000)
                Rn_3(Ind_001) = delz*Rn_2(Ind_000)
                Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
                Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
                Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
                Rn_3(Ind_110) = delx*Rn_2(Ind_010)
                Rn_3(Ind_101) = delx*Rn_2(Ind_001)
                Rn_3(Ind_011) = dely*Rn_2(Ind_001)
                Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
                Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
                Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
                Rn_3(Ind_210) = dely*Rn_2(Ind_200)
                Rn_3(Ind_201) = delz*Rn_2(Ind_200)
                Rn_3(Ind_120) = delx*Rn_2(Ind_020)
                Rn_3(Ind_021) = delz*Rn_2(Ind_020)
                Rn_3(Ind_102) = delx*Rn_2(Ind_002)
                Rn_3(Ind_012) = dely*Rn_2(Ind_002)
                Rn_3(Ind_111) = delx*Rn_2(Ind_011)
                ! finally get field components
                if (is_polarizable(i)) then
                    gphi_i(1) = Rn_3(Ind_100)*global_multipole(Ind_000, j) + &
                        Rn_3(Ind_200)*global_multipole(Ind_100, j) + &
                        Rn_3(Ind_110)*global_multipole(Ind_010, j) + &
                        Rn_3(Ind_101)*global_multipole(Ind_001, j) + &
                        Rn_3(Ind_300)*global_multipole(Ind_200, j) + &
                        Rn_3(Ind_120)*global_multipole(Ind_020, j) + &
                        Rn_3(Ind_102)*global_multipole(Ind_002, j) + &
                        Rn_3(Ind_210)*global_multipole(Ind_110, j) + &
                        Rn_3(Ind_201)*global_multipole(Ind_101, j) + &
                        Rn_3(Ind_111)*global_multipole(Ind_011, j)
                    gphi_i(2) = Rn_3(Ind_010)*global_multipole(Ind_000, j) + &
                        Rn_3(Ind_110)*global_multipole(Ind_100, j) + &
                        Rn_3(Ind_020)*global_multipole(Ind_010, j) + &
                        Rn_3(Ind_011)*global_multipole(Ind_001, j) + &
                        Rn_3(Ind_210)*global_multipole(Ind_200, j) + &
                        Rn_3(Ind_030)*global_multipole(Ind_020, j) + &
                        Rn_3(Ind_012)*global_multipole(Ind_002, j) + &
                        Rn_3(Ind_120)*global_multipole(Ind_110, j) + &
                        Rn_3(Ind_111)*global_multipole(Ind_101, j) + &
                        Rn_3(Ind_021)*global_multipole(Ind_011, j)
                    gphi_i(3) = Rn_3(Ind_001)*global_multipole(Ind_000, j) + &
                        Rn_3(Ind_101)*global_multipole(Ind_100, j) + &
                        Rn_3(Ind_011)*global_multipole(Ind_010, j) + &
                        Rn_3(Ind_002)*global_multipole(Ind_001, j) + &
                        Rn_3(Ind_201)*global_multipole(Ind_200, j) + &
                        Rn_3(Ind_021)*global_multipole(Ind_020, j) + &
                        Rn_3(Ind_003)*global_multipole(Ind_002, j) + &
                        Rn_3(Ind_111)*global_multipole(Ind_110, j) + &
                        Rn_3(Ind_102)*global_multipole(Ind_101, j) + &
                        Rn_3(Ind_012)*global_multipole(Ind_011, j)
                    !negative for i since d/dx_i = -d/delx
                    direct_gradphi(1, i) = direct_gradphi(1, i) - gphi_i(1)
                    direct_gradphi(2, i) = direct_gradphi(2, i) - gphi_i(2)
                    direct_gradphi(3, i) = direct_gradphi(3, i) - gphi_i(3)
                end if ! is_polarizable(i)
                if (is_polarizable(j)) then
                    ! note negative contribs for dipoles of i since d/dx_i = -d/delx)
                    ! (look at contributions to electrostatic potential at j -these will
                    !  contain negative contributions due to dipolar contribs at i--the
                    !  extra deriv in tensor due to grad at j has no negative signs)
                    gphi_j(1) = Rn_3(Ind_100)*global_multipole(Ind_000, i) - &
                        Rn_3(Ind_200)*global_multipole(Ind_100, i) - &
                        Rn_3(Ind_110)*global_multipole(Ind_010, i) - &
                        Rn_3(Ind_101)*global_multipole(Ind_001, i) + &
                        Rn_3(Ind_300)*global_multipole(Ind_200, i) + &
                        Rn_3(Ind_120)*global_multipole(Ind_020, i) + &
                        Rn_3(Ind_102)*global_multipole(Ind_002, i) + &
                        Rn_3(Ind_210)*global_multipole(Ind_110, i) + &
                        Rn_3(Ind_201)*global_multipole(Ind_101, i) + &
                        Rn_3(Ind_111)*global_multipole(Ind_011, i)
                    gphi_j(2) = Rn_3(Ind_010)*global_multipole(Ind_000, i) - &
                        Rn_3(Ind_110)*global_multipole(Ind_100, i) - &
                        Rn_3(Ind_020)*global_multipole(Ind_010, i) - &
                        Rn_3(Ind_011)*global_multipole(Ind_001, i) + &
                        Rn_3(Ind_210)*global_multipole(Ind_200, i) + &
                        Rn_3(Ind_030)*global_multipole(Ind_020, i) + &
                        Rn_3(Ind_012)*global_multipole(Ind_002, i) + &
                        Rn_3(Ind_120)*global_multipole(Ind_110, i) + &
                        Rn_3(Ind_111)*global_multipole(Ind_101, i) + &
                        Rn_3(Ind_021)*global_multipole(Ind_011, i)
                    gphi_j(3) = Rn_3(Ind_001)*global_multipole(Ind_000, i) - &
                        Rn_3(Ind_101)*global_multipole(Ind_100, i) - &
                        Rn_3(Ind_011)*global_multipole(Ind_010, i) - &
                        Rn_3(Ind_002)*global_multipole(Ind_001, i) + &
                        Rn_3(Ind_201)*global_multipole(Ind_200, i) + &
                        Rn_3(Ind_021)*global_multipole(Ind_020, i) + &
                        Rn_3(Ind_003)*global_multipole(Ind_002, i) + &
                        Rn_3(Ind_111)*global_multipole(Ind_110, i) + &
                        Rn_3(Ind_102)*global_multipole(Ind_101, i) + &
                        Rn_3(Ind_012)*global_multipole(Ind_011, i)
                    direct_gradphi(1, j) = direct_gradphi(1, j) + gphi_j(1)
                    direct_gradphi(2, j) = direct_gradphi(2, j) + gphi_j(2)
                    direct_gradphi(3, j) = direct_gradphi(3, j) + gphi_j(3)
                end if !( is_polarizable(j) )then
                ! next do the polar field --negate the odd order terms
                D(1) = -(polar_weight(k)*(A(1) - BD(1)) + (B(1) - A(1)))
                D(2) = polar_weight(k)*(A(2) - BD(2)) + (B(2) - A(2))
                D(3) = -(polar_weight(k)*(A(3) - BD(3)) + (B(3) - A(3)))
                n = 3
                Rn(Ind_000) = D(n)
                Rn_1(Ind_000) = D(n - 1)
                Rn_1(Ind_100) = delx*Rn(Ind_000)
                Rn_1(Ind_010) = dely*Rn(Ind_000)
                Rn_1(Ind_001) = delz*Rn(Ind_000)
                Rn_2(Ind_000) = D(n - 2)
                Rn_2(Ind_100) = delx*Rn_1(Ind_000)
                Rn_2(Ind_010) = dely*Rn_1(Ind_000)
                Rn_2(Ind_001) = delz*Rn_1(Ind_000)
                Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
                Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
                Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
                Rn_2(Ind_110) = delx*Rn_1(Ind_010)
                Rn_2(Ind_101) = delx*Rn_1(Ind_001)
                Rn_2(Ind_011) = dely*Rn_1(Ind_001)
                Rn_3(Ind_100) = delx*Rn_2(Ind_000)
                Rn_3(Ind_010) = dely*Rn_2(Ind_000)
                Rn_3(Ind_001) = delz*Rn_2(Ind_000)
                Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
                Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
                Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
                Rn_3(Ind_110) = delx*Rn_2(Ind_010)
                Rn_3(Ind_101) = delx*Rn_2(Ind_001)
                Rn_3(Ind_011) = dely*Rn_2(Ind_001)
                Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
                Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
                Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
                Rn_3(Ind_210) = dely*Rn_2(Ind_200)
                Rn_3(Ind_201) = delz*Rn_2(Ind_200)
                Rn_3(Ind_120) = delx*Rn_2(Ind_020)
                Rn_3(Ind_021) = delz*Rn_2(Ind_020)
                Rn_3(Ind_102) = delx*Rn_2(Ind_002)
                Rn_3(Ind_012) = dely*Rn_2(Ind_002)
                Rn_3(Ind_111) = delx*Rn_2(Ind_011)
                ! finally get field components
                if (is_polarizable(i)) then
                    gphi_i(1) = Rn_3(Ind_100)*global_multipole(Ind_000, j) + &
                        Rn_3(Ind_200)*global_multipole(Ind_100, j) + &
                        Rn_3(Ind_110)*global_multipole(Ind_010, j) + &
                        Rn_3(Ind_101)*global_multipole(Ind_001, j) + &
                        Rn_3(Ind_300)*global_multipole(Ind_200, j) + &
                        Rn_3(Ind_120)*global_multipole(Ind_020, j) + &
                        Rn_3(Ind_102)*global_multipole(Ind_002, j) + &
                        Rn_3(Ind_210)*global_multipole(Ind_110, j) + &
                        Rn_3(Ind_201)*global_multipole(Ind_101, j) + &
                        Rn_3(Ind_111)*global_multipole(Ind_011, j)
                    gphi_i(2) = Rn_3(Ind_010)*global_multipole(Ind_000, j) + &
                        Rn_3(Ind_110)*global_multipole(Ind_100, j) + &
                        Rn_3(Ind_020)*global_multipole(Ind_010, j) + &
                        Rn_3(Ind_011)*global_multipole(Ind_001, j) + &
                        Rn_3(Ind_210)*global_multipole(Ind_200, j) + &
                        Rn_3(Ind_030)*global_multipole(Ind_020, j) + &
                        Rn_3(Ind_012)*global_multipole(Ind_002, j) + &
                        Rn_3(Ind_120)*global_multipole(Ind_110, j) + &
                        Rn_3(Ind_111)*global_multipole(Ind_101, j) + &
                        Rn_3(Ind_021)*global_multipole(Ind_011, j)
                    gphi_i(3) = Rn_3(Ind_001)*global_multipole(Ind_000, j) + &
                        Rn_3(Ind_101)*global_multipole(Ind_100, j) + &
                        Rn_3(Ind_011)*global_multipole(Ind_010, j) + &
                        Rn_3(Ind_002)*global_multipole(Ind_001, j) + &
                        Rn_3(Ind_201)*global_multipole(Ind_200, j) + &
                        Rn_3(Ind_021)*global_multipole(Ind_020, j) + &
                        Rn_3(Ind_003)*global_multipole(Ind_002, j) + &
                        Rn_3(Ind_111)*global_multipole(Ind_110, j) + &
                        Rn_3(Ind_102)*global_multipole(Ind_101, j) + &
                        Rn_3(Ind_012)*global_multipole(Ind_011, j)
                    !negative for i since d/dx_i = -d/delx
                    polar_gradphi(1, i) = polar_gradphi(1, i) - gphi_i(1)
                    polar_gradphi(2, i) = polar_gradphi(2, i) - gphi_i(2)
                    polar_gradphi(3, i) = polar_gradphi(3, i) - gphi_i(3)
                end if ! is_polarizable(i)
                if (is_polarizable(j)) then
                    ! note negative contribs for dipoles of i since d/dx_i = -d/delx)
                    ! (look at contributions to electrostatic potential at j -these will
                    !  contain negative contributions due to dipolar contribs at i--the
                    !  extra deriv in tensor due to grad at j has no negative signs)
                    gphi_j(1) = Rn_3(Ind_100)*global_multipole(Ind_000, i) - &
                        Rn_3(Ind_200)*global_multipole(Ind_100, i) - &
                        Rn_3(Ind_110)*global_multipole(Ind_010, i) - &
                        Rn_3(Ind_101)*global_multipole(Ind_001, i) + &
                        Rn_3(Ind_300)*global_multipole(Ind_200, i) + &
                        Rn_3(Ind_120)*global_multipole(Ind_020, i) + &
                        Rn_3(Ind_102)*global_multipole(Ind_002, i) + &
                        Rn_3(Ind_210)*global_multipole(Ind_110, i) + &
                        Rn_3(Ind_201)*global_multipole(Ind_101, i) + &
                        Rn_3(Ind_111)*global_multipole(Ind_011, i)
                    gphi_j(2) = Rn_3(Ind_010)*global_multipole(Ind_000, i) - &
                        Rn_3(Ind_110)*global_multipole(Ind_100, i) - &
                        Rn_3(Ind_020)*global_multipole(Ind_010, i) - &
                        Rn_3(Ind_011)*global_multipole(Ind_001, i) + &
                        Rn_3(Ind_210)*global_multipole(Ind_200, i) + &
                        Rn_3(Ind_030)*global_multipole(Ind_020, i) + &
                        Rn_3(Ind_012)*global_multipole(Ind_002, i) + &
                        Rn_3(Ind_120)*global_multipole(Ind_110, i) + &
                        Rn_3(Ind_111)*global_multipole(Ind_101, i) + &
                        Rn_3(Ind_021)*global_multipole(Ind_011, i)
                    gphi_j(3) = Rn_3(Ind_001)*global_multipole(Ind_000, i) - &
                        Rn_3(Ind_101)*global_multipole(Ind_100, i) - &
                        Rn_3(Ind_011)*global_multipole(Ind_010, i) - &
                        Rn_3(Ind_002)*global_multipole(Ind_001, i) + &
                        Rn_3(Ind_201)*global_multipole(Ind_200, i) + &
                        Rn_3(Ind_021)*global_multipole(Ind_020, i) + &
                        Rn_3(Ind_003)*global_multipole(Ind_002, i) + &
                        Rn_3(Ind_111)*global_multipole(Ind_110, i) + &
                        Rn_3(Ind_102)*global_multipole(Ind_101, i) + &
                        Rn_3(Ind_012)*global_multipole(Ind_011, i)
                    polar_gradphi(1, j) = polar_gradphi(1, j) + gphi_j(1)
                    polar_gradphi(2, j) = polar_gradphi(2, j) + gphi_j(2)
                    polar_gradphi(3, j) = polar_gradphi(3, j) + gphi_j(3)
                end if !( is_polarizable(j) )then
            end if ! (is_polarizable(i) .or. is_polarizable(j) )then
        end do
    end subroutine AM_ADJUST_calc_permfield
!-------------------------------------------------------
    subroutine AM_ADJUST_dip_dip_fields( &
        ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)
        _REAL_, intent(in) :: ind_dip_d(3, *), ind_dip_p(3, *)
        _REAL_, intent(inout) :: dip_field_d(3, *), dip_field_p(3, *)

#  include "do_flag.h"

        if (iand(do_flag, PROCEED_INDUCE) /= PROCEED_INDUCE) return

        call timer_start(TIME_ADJ)
        call AM_ADJUST_calc_dip_dip_fields( &
            num_adjust_list, adjust_list, &
            dipole_dipole_tensor, &
            ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)
        call timer_stop(TIME_ADJ)
    end subroutine AM_ADJUST_dip_dip_fields
!-------------------------------------------------------
    subroutine AM_ADJUST_calc_dip_dip_fields( &
        num_adjust_list, adjust_list, &
        dipole_dipole_tensor, &
        ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)
        integer, intent(in) :: num_adjust_list, adjust_list(3, *)
        _REAL_, intent(in) :: dipole_dipole_tensor(6, *)
        _REAL_, intent(in) :: ind_dip_d(3, *), ind_dip_p(3, *)
        _REAL_, intent(inout) :: dip_field_d(3, *), dip_field_p(3, *)

        integer :: i, j, n
        do n = 1, num_adjust_list
            i = adjust_list(1, n)
            j = adjust_list(2, n)
            ! minus signs due to deriv wrt position of i
            dip_field_d(1, i) = dip_field_d(1, i) - &
                (dipole_dipole_tensor(1, n)*ind_dip_d(1, j) + &
                dipole_dipole_tensor(2, n)*ind_dip_d(2, j) + &
                dipole_dipole_tensor(3, n)*ind_dip_d(3, j))
            dip_field_d(2, i) = dip_field_d(2, i) - &
                (dipole_dipole_tensor(2, n)*ind_dip_d(1, j) + &
                dipole_dipole_tensor(4, n)*ind_dip_d(2, j) + &
                dipole_dipole_tensor(5, n)*ind_dip_d(3, j))
            dip_field_d(3, i) = dip_field_d(3, i) - &
                (dipole_dipole_tensor(3, n)*ind_dip_d(1, j) + &
                dipole_dipole_tensor(5, n)*ind_dip_d(2, j) + &
                dipole_dipole_tensor(6, n)*ind_dip_d(3, j))
            dip_field_d(1, j) = dip_field_d(1, j) - &
                (dipole_dipole_tensor(1, n)*ind_dip_d(1, i) + &
                dipole_dipole_tensor(2, n)*ind_dip_d(2, i) + &
                dipole_dipole_tensor(3, n)*ind_dip_d(3, i))
            dip_field_d(2, j) = dip_field_d(2, j) - &
                (dipole_dipole_tensor(2, n)*ind_dip_d(1, i) + &
                dipole_dipole_tensor(4, n)*ind_dip_d(2, i) + &
                dipole_dipole_tensor(5, n)*ind_dip_d(3, i))
            dip_field_d(3, j) = dip_field_d(3, j) - &
                (dipole_dipole_tensor(3, n)*ind_dip_d(1, i) + &
                dipole_dipole_tensor(5, n)*ind_dip_d(2, i) + &
                dipole_dipole_tensor(6, n)*ind_dip_d(3, i))
            ! now do other set of dipoles
            dip_field_p(1, i) = dip_field_p(1, i) - &
                (dipole_dipole_tensor(1, n)*ind_dip_p(1, j) + &
                dipole_dipole_tensor(2, n)*ind_dip_p(2, j) + &
                dipole_dipole_tensor(3, n)*ind_dip_p(3, j))
            dip_field_p(2, i) = dip_field_p(2, i) - &
                (dipole_dipole_tensor(2, n)*ind_dip_p(1, j) + &
                dipole_dipole_tensor(4, n)*ind_dip_p(2, j) + &
                dipole_dipole_tensor(5, n)*ind_dip_p(3, j))
            dip_field_p(3, i) = dip_field_p(3, i) - &
                (dipole_dipole_tensor(3, n)*ind_dip_p(1, j) + &
                dipole_dipole_tensor(5, n)*ind_dip_p(2, j) + &
                dipole_dipole_tensor(6, n)*ind_dip_p(3, j))
            dip_field_p(1, j) = dip_field_p(1, j) - &
                (dipole_dipole_tensor(1, n)*ind_dip_p(1, i) + &
                dipole_dipole_tensor(2, n)*ind_dip_p(2, i) + &
                dipole_dipole_tensor(3, n)*ind_dip_p(3, i))
            dip_field_p(2, j) = dip_field_p(2, j) - &
                (dipole_dipole_tensor(2, n)*ind_dip_p(1, i) + &
                dipole_dipole_tensor(4, n)*ind_dip_p(2, i) + &
                dipole_dipole_tensor(5, n)*ind_dip_p(3, i))
            dip_field_p(3, j) = dip_field_p(3, j) - &
                (dipole_dipole_tensor(3, n)*ind_dip_p(1, i) + &
                dipole_dipole_tensor(5, n)*ind_dip_p(2, i) + &
                dipole_dipole_tensor(6, n)*ind_dip_p(3, i))
        end do
    end subroutine AM_ADJUST_calc_dip_dip_fields
!-------------------------------------------------------
    subroutine AM_ADJUST_ene_frc(crd, x, ind_dip_d, ind_dip_p, &
        ene_perm, ene_ind, ene_vdw, frc, virial)

        use amoeba_multipoles, only : global_multipole
        use amoeba_induced, only : sq_polinv, is_polarizable
        use amoeba_mdin, only : thole_expon_coeff
        use amoeba_vdw, only : AM_VDW_ADJUST_ene_frc
        use constants, only : zero

        _REAL_, intent(in) :: crd(3, *), x(*), ind_dip_d(3, *), ind_dip_p(3, *)
        _REAL_, intent(inout) :: ene_perm, ene_ind, ene_vdw, frc(3, *), virial(3, 3)

#  include "def_time.h"
#  include "ew_erfc_spline.h"
#  include "ew_pme_recip.h"

        ene_perm = zero
        ene_ind = zero
        ene_vdw = zero
        call timer_start(TIME_VDW)
        call AM_VDW_ADJUST_ene_frc(crd, num_adjust_list, adjust_list, vdw_weight, &
            ene_vdw, frc, virial)
        call timer_stop_start(TIME_VDW, TIME_ADJ)
        call AM_ADJUST_calc_ene_frc(crd, num_adjust_list, adjust_list, &
            ew_coeff, eedtbdns, x(leed_cub), x(leed_lin), &
            ee_type, eedmeth, dxdr, thole_expon_coeff, &
            sq_polinv, is_polarizable, global_multipole, &
            ind_dip_d, ind_dip_p, &
            mpole_weight, polar_weight, &
            direct_weight, mutual_weight, &
            ene_perm, ene_ind, frc, virial)
        call timer_stop(TIME_ADJ)
    end subroutine AM_ADJUST_ene_frc
!-------------------------------------------------------
    subroutine AM_ADJUST_calc_ene_frc(crd, num_adjust_list, adjust_list, &
        ewaldcof, eedtbdns, eed_cub, eed_lin, &
        ee_type, eedmeth, dxdr, thole_expon_coeff, &
        sq_polinv, is_polarizable, global_multipole, &
        ind_dip_d, ind_dip_p, &
        mpole_weight, polar_weight, &
        direct_weight, mutual_weight, &
        ene_perm, ene_ind, frc, virial)
        use amoeba_multipoles, only : coulomb_const_kcal_per_mole, torque_field
        use constants, only : zero, third, half, one, two, three, five, seven, nine
        _REAL_, intent(in) :: crd(3, *)
        integer, intent(in) :: num_adjust_list, adjust_list(3, *)
        _REAL_, intent(in) :: ewaldcof, eedtbdns, eed_cub(4, *), eed_lin(2, *)
        integer, intent(in) :: ee_type, eedmeth
        _REAL_, intent(in) :: dxdr, thole_expon_coeff, sq_polinv(*)
        logical, intent(in) :: is_polarizable(*)
        _REAL_, intent(in) :: global_multipole(10, *), ind_dip_d(3, *), ind_dip_p(3, *)
        _REAL_, intent(in) :: mpole_weight(*), polar_weight(*)
        _REAL_, intent(in) :: direct_weight(*), mutual_weight(*)
        _REAL_, intent(inout) :: ene_perm, ene_ind, frc(3, *), virial(3, 3)

        _REAL_ :: delx, dely, delz, delr2, delr, delr2inv, x, dx, switch, d_switch_dx, xx
        _REAL_ :: vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz
        _REAL_ :: A(0:5), B(0:5), C(0:5), BD(4), fac, fact, del
        _REAL_ :: expon, expo, clam3, clam5, clam7, clam9, &
            delr3inv, delr5inv, delr7inv, delr9inv
        _REAL_ :: Rn(1), Rn_1(4), Rn_2(10), Rn_3(20), Rn_4(35), Rn_5(56)
        _REAL_ :: gmi(10), gmj(10), phi(20), e_pp, g_pp(3), e_ind, g_ind(3), &
            i_di(3), i_dj(3), i_pi(3), i_pj(3)
        _REAL_, parameter :: const1 = 0.6d0, const2 = 18.d0/35.d0, const3 = 9.d0/35.d0

        integer :: i, j, k, n, n_adj, ind, jj
#  include "do_flag.h"

        if (iand(do_flag, PROCEED_POSTINDUCE) /= PROCEED_POSTINDUCE) return

        fac = two*ewaldcof*ewaldcof
        del = one/eedtbdns
        vxx = zero
        vxy = zero
        vxz = zero
        vyx = zero
        vyy = zero
        vyz = zero
        vzx = zero
        vzy = zero
        vzz = zero
        do n_adj = 1, num_adjust_list
            i = adjust_list(1, n_adj)
            j = adjust_list(2, n_adj)
            k = adjust_list(3, n_adj)
            delx = crd(1, j) - crd(1, i)
            dely = crd(2, j) - crd(2, i)
            delz = crd(3, j) - crd(3, i)
            delr2 = delx*delx + dely*dely + delz*delz
            delr = sqrt(delr2)
            delr2inv = one/delr2
            x = dxdr*delr
            if (eedmeth == 1) then
                !           -- cubic spline on switch
                ind = eedtbdns*x + 1
                dx = x - (ind - one)*del
                switch = eed_cub(1, ind) + dx*(eed_cub(2, ind) + &
                    dx*(eed_cub(3, ind) + dx*eed_cub(4, ind)*third)*half)
                d_switch_dx = eed_cub(2, ind) + dx*(eed_cub(3, ind) + &
                    dx*eed_cub(4, ind)*half)
            else if (eedmeth == 2) then
                !           ---linear lookup on switch, deriv
                xx = eedtbdns*x + 1
                ind = xx
                dx = xx - ind
                switch = (one - dx)*eed_lin(1, ind) + dx*eed_lin(1, ind + 1)
                d_switch_dx = (one - dx)*eed_lin(2, ind) + dx*eed_lin(2, ind + 1)
            else if (eedmeth == 3) then
                !           ---direct function call:
                call get_ee_func(x, switch, d_switch_dx, ee_type)
            else if (eedmeth == 4) then
                write (6, *) 'eedmeth is 4 inside AM_ADJUST_calc_ene_frc!!'
                call mexit(6, 1)
            end if
            !------------------------------------------------------------
            ! McMurchie-Davidson recursion holds for any smooth function of r
            ! that is, to get the higher order derivs wrt x,y,z of g(r)
            ! define R(0,0,0,0) = g(r)
            ! next  R(0,0,0,n+1) = -(1/r)d/dr R(0,0,0,n)
            ! then denote  R(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
            ! quantities of interest obtained by setting n=0
            ! McMurchie-Davidson says that
            !  R(t+1,u,v,n) = t*R(t-1,u,v,n+1) + x*R(t,u,v,n)
            !  R(t,u+1,v,n) = u*R(t,u-1,v,n+1) + y*R(t,u,v,n)
            !  R(t,u,v+1,n) = v*R(t,u,v-1,n+1) + z*R(t,u,v,n)
            ! use similar data structures as in direct sum
            !-------------------------------------------------------------
            ! calc the contributions for cancelling reciprocal sum for ij
            ! -erf = erfc - 1.0 ---first get the erfc part as in direct sum
            B(0) = switch*delr*delr2inv
            fact = d_switch_dx*dxdr
            B(1) = (B(0) - fact)*delr2inv
            fact = fac*fact
            B(2) = (three*B(1) - fact)*delr2inv
            fact = fac*fact
            B(3) = (five*B(2) - fact)*delr2inv
            fact = fac*fact
            B(4) = (seven*B(3) - fact)*delr2inv
            fact = fac*fact
            B(5) = (nine*B(4) - fact)*delr2inv
            ! the -erf part is given by B(j) - A(j); weighted direct by
            ! mpole_weight*A(j)
            A(0) = one*delr*delr2inv
            A(1) = A(0)*delr2inv
            A(2) = three*A(1)*delr2inv
            A(3) = five*A(2)*delr2inv
            A(4) = seven*A(3)*delr2inv
            A(5) = nine*A(4)*delr2inv
            ! damped part for permanent-induced or induced-induced interactions
            delr3inv = delr2inv/delr
            delr5inv = delr3inv*delr2inv
            delr7inv = delr5inv*delr2inv
            delr9inv = delr7inv*delr2inv
            expon = thole_expon_coeff*delr2*delr*sq_polinv(i)*sq_polinv(j)
            expo = exp(-expon)
            ! clam3 = 1.d0-lam3, clam5 = 1.d0-lam5 etc. where
            ! lam is from ponder's paper
            clam3 = expo
            clam5 = (one + expon)*expo
            clam7 = (one + expon + const1*expon**2)*expo
            clam9 = (one + expon + const2*expon**2 + const3*expon**3)*expo
            BD(1) = clam3*delr3inv
            BD(2) = three*clam5*delr5inv
            BD(3) = 15.d0*clam7*delr7inv
            BD(4) = 105.d0*clam9*delr9inv
            ! first handle the permanent-permanent interactions
            do jj = 0, 5
                C(jj) = (B(jj) - A(jj)) + mpole_weight(k)*A(jj)
            end do
            ! negate the odd order to get sign right
            C(1) = -C(1)
            C(3) = -C(3)
            C(5) = -C(5)
            ! get the interaction tensor
            n = 5
            Rn(Ind_000) = C(n)
            Rn_1(Ind_000) = C(n - 1)
            Rn_1(Ind_100) = delx*Rn(Ind_000)
            Rn_1(Ind_010) = dely*Rn(Ind_000)
            Rn_1(Ind_001) = delz*Rn(Ind_000)
            Rn_2(Ind_000) = C(n - 2)
            Rn_2(Ind_100) = delx*Rn_1(Ind_000)
            Rn_2(Ind_010) = dely*Rn_1(Ind_000)
            Rn_2(Ind_001) = delz*Rn_1(Ind_000)
            Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
            Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
            Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
            Rn_2(Ind_110) = delx*Rn_1(Ind_010)
            Rn_2(Ind_101) = delx*Rn_1(Ind_001)
            Rn_2(Ind_011) = dely*Rn_1(Ind_001)
            Rn_3(Ind_000) = C(n - 3)
            Rn_3(Ind_100) = delx*Rn_2(Ind_000)
            Rn_3(Ind_010) = dely*Rn_2(Ind_000)
            Rn_3(Ind_001) = delz*Rn_2(Ind_000)
            Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
            Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
            Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
            Rn_3(Ind_110) = delx*Rn_2(Ind_010)
            Rn_3(Ind_101) = delx*Rn_2(Ind_001)
            Rn_3(Ind_011) = dely*Rn_2(Ind_001)
            Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
            Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
            Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
            Rn_3(Ind_210) = dely*Rn_2(Ind_200)
            Rn_3(Ind_201) = delz*Rn_2(Ind_200)
            Rn_3(Ind_120) = delx*Rn_2(Ind_020)
            Rn_3(Ind_021) = delz*Rn_2(Ind_020)
            Rn_3(Ind_102) = delx*Rn_2(Ind_002)
            Rn_3(Ind_012) = dely*Rn_2(Ind_002)
            Rn_3(Ind_111) = delx*Rn_2(Ind_011)
            Rn_4(Ind_000) = C(n - 4)
            Rn_4(Ind_100) = delx*Rn_3(Ind_000)
            Rn_4(Ind_010) = dely*Rn_3(Ind_000)
            Rn_4(Ind_001) = delz*Rn_3(Ind_000)
            Rn_4(Ind_200) = Rn_3(Ind_000) + delx*Rn_3(Ind_100)
            Rn_4(Ind_020) = Rn_3(Ind_000) + dely*Rn_3(Ind_010)
            Rn_4(Ind_002) = Rn_3(Ind_000) + delz*Rn_3(Ind_001)
            Rn_4(Ind_110) = delx*Rn_3(Ind_010)
            Rn_4(Ind_101) = delx*Rn_3(Ind_001)
            Rn_4(Ind_011) = dely*Rn_3(Ind_001)
            Rn_4(Ind_300) = 2.d0*Rn_3(Ind_100) + delx*Rn_3(Ind_200)
            Rn_4(Ind_030) = 2.d0*Rn_3(Ind_010) + dely*Rn_3(Ind_020)
            Rn_4(Ind_003) = 2.d0*Rn_3(Ind_001) + delz*Rn_3(Ind_002)
            Rn_4(Ind_210) = dely*Rn_3(Ind_200)
            Rn_4(Ind_201) = delz*Rn_3(Ind_200)
            Rn_4(Ind_120) = delx*Rn_3(Ind_020)
            Rn_4(Ind_021) = delz*Rn_3(Ind_020)
            Rn_4(Ind_102) = delx*Rn_3(Ind_002)
            Rn_4(Ind_012) = dely*Rn_3(Ind_002)
            Rn_4(Ind_111) = delx*Rn_3(Ind_011)
            Rn_4(Ind_400) = 3.d0*Rn_3(Ind_200) + delx*Rn_3(Ind_300)
            Rn_4(Ind_040) = 3.d0*Rn_3(Ind_020) + dely*Rn_3(Ind_030)
            Rn_4(Ind_004) = 3.d0*Rn_3(Ind_002) + delz*Rn_3(Ind_003)
            Rn_4(Ind_310) = dely*Rn_3(Ind_300)
            Rn_4(Ind_301) = delz*Rn_3(Ind_300)
            Rn_4(Ind_130) = delx*Rn_3(Ind_030)
            Rn_4(Ind_031) = delz*Rn_3(Ind_030)
            Rn_4(Ind_103) = delx*Rn_3(Ind_003)
            Rn_4(Ind_013) = dely*Rn_3(Ind_003)
            Rn_4(Ind_220) = Rn_3(Ind_020) + delx*Rn_3(Ind_120)
            Rn_4(Ind_202) = Rn_3(Ind_002) + delx*Rn_3(Ind_102)
            Rn_4(Ind_022) = Rn_3(Ind_002) + dely*Rn_3(Ind_012)
            Rn_4(Ind_211) = dely*Rn_3(Ind_201)
            Rn_4(Ind_121) = delx*Rn_3(Ind_021)
            Rn_4(Ind_112) = delx*Rn_3(Ind_012)
            Rn_5(Ind_000) = C(n - 5)
            Rn_5(Ind_100) = delx*Rn_4(Ind_000)
            Rn_5(Ind_010) = dely*Rn_4(Ind_000)
            Rn_5(Ind_001) = delz*Rn_4(Ind_000)
            Rn_5(Ind_200) = Rn_4(Ind_000) + delx*Rn_4(Ind_100)
            Rn_5(Ind_020) = Rn_4(Ind_000) + dely*Rn_4(Ind_010)
            Rn_5(Ind_002) = Rn_4(Ind_000) + delz*Rn_4(Ind_001)
            Rn_5(Ind_110) = delx*Rn_4(Ind_010)
            Rn_5(Ind_101) = delx*Rn_4(Ind_001)
            Rn_5(Ind_011) = dely*Rn_4(Ind_001)
            Rn_5(Ind_300) = 2.d0*Rn_4(Ind_100) + delx*Rn_4(Ind_200)
            Rn_5(Ind_030) = 2.d0*Rn_4(Ind_010) + dely*Rn_4(Ind_020)
            Rn_5(Ind_003) = 2.d0*Rn_4(Ind_001) + delz*Rn_4(Ind_002)
            Rn_5(Ind_210) = dely*Rn_4(Ind_200)
            Rn_5(Ind_201) = delz*Rn_4(Ind_200)
            Rn_5(Ind_120) = delx*Rn_4(Ind_020)
            Rn_5(Ind_021) = delz*Rn_4(Ind_020)
            Rn_5(Ind_102) = delx*Rn_4(Ind_002)
            Rn_5(Ind_012) = dely*Rn_4(Ind_002)
            Rn_5(Ind_111) = delx*Rn_4(Ind_011)
            Rn_5(Ind_400) = 3.d0*Rn_4(Ind_200) + delx*Rn_4(Ind_300)
            Rn_5(Ind_040) = 3.d0*Rn_4(Ind_020) + dely*Rn_4(Ind_030)
            Rn_5(Ind_004) = 3.d0*Rn_4(Ind_002) + delz*Rn_4(Ind_003)
            Rn_5(Ind_310) = dely*Rn_4(Ind_300)
            Rn_5(Ind_301) = delz*Rn_4(Ind_300)
            Rn_5(Ind_130) = delx*Rn_4(Ind_030)
            Rn_5(Ind_031) = delz*Rn_4(Ind_030)
            Rn_5(Ind_103) = delx*Rn_4(Ind_003)
            Rn_5(Ind_013) = dely*Rn_4(Ind_003)
            Rn_5(Ind_220) = Rn_4(Ind_020) + delx*Rn_4(Ind_120)
            Rn_5(Ind_202) = Rn_4(Ind_002) + delx*Rn_4(Ind_102)
            Rn_5(Ind_022) = Rn_4(Ind_002) + dely*Rn_4(Ind_012)
            Rn_5(Ind_211) = dely*Rn_4(Ind_201)
            Rn_5(Ind_121) = delx*Rn_4(Ind_021)
            Rn_5(Ind_112) = delx*Rn_4(Ind_012)
            Rn_5(Ind_500) = 4.d0*Rn_4(Ind_300) + delx*Rn_4(Ind_400)
            Rn_5(Ind_050) = 4.d0*Rn_4(Ind_030) + dely*Rn_4(Ind_040)
            Rn_5(Ind_005) = 4.d0*Rn_4(Ind_003) + delz*Rn_4(Ind_004)
            Rn_5(Ind_410) = dely*Rn_4(Ind_400)
            Rn_5(Ind_401) = delz*Rn_4(Ind_400)
            Rn_5(Ind_140) = delx*Rn_4(Ind_040)
            Rn_5(Ind_041) = delz*Rn_4(Ind_040)
            Rn_5(Ind_104) = delx*Rn_4(Ind_004)
            Rn_5(Ind_014) = dely*Rn_4(Ind_004)
            Rn_5(Ind_320) = Rn_4(Ind_300) + dely*Rn_4(Ind_310)
            Rn_5(Ind_302) = Rn_4(Ind_300) + delz*Rn_4(Ind_301)
            Rn_5(Ind_230) = Rn_4(Ind_030) + delx*Rn_4(Ind_130)
            Rn_5(Ind_032) = Rn_4(Ind_030) + delz*Rn_4(Ind_031)
            Rn_5(Ind_203) = Rn_4(Ind_003) + delx*Rn_4(Ind_103)
            Rn_5(Ind_023) = Rn_4(Ind_003) + dely*Rn_4(Ind_013)
            Rn_5(Ind_311) = dely*Rn_4(Ind_301)
            Rn_5(Ind_131) = delx*Rn_4(Ind_031)
            Rn_5(Ind_113) = delx*Rn_4(Ind_013)
            Rn_5(Ind_221) = delz*Rn_4(Ind_220)
            Rn_5(Ind_212) = dely*Rn_4(Ind_202)
            Rn_5(Ind_122) = delx*Rn_4(Ind_022)
            ! phi array (electrostatic potential at i due to j permanent moments
            ! and derivs of that esp wrt r_i )
            ! minus signs arise due to derivs of r_j - r_i wrt r_i
            do jj = 1, 10
                gmj(jj) = global_multipole(jj, j)
                gmi(jj) = global_multipole(jj, i)
            end do
            phi(Ind_000) = Rn_5(Ind_000)*gmj(Ind_000) + Rn_5(Ind_100)*gmj(Ind_100) + &
                Rn_5(Ind_010)*gmj(Ind_010) + Rn_5(Ind_001)*gmj(Ind_001) + &
                Rn_5(Ind_200)*gmj(Ind_200) + Rn_5(Ind_020)*gmj(Ind_020) + &
                Rn_5(Ind_002)*gmj(Ind_002) + Rn_5(Ind_110)*gmj(Ind_110) + &
                Rn_5(Ind_101)*gmj(Ind_101) + Rn_5(Ind_011)*gmj(Ind_011)
            phi(Ind_100) = -(Rn_5(Ind_100)*gmj(Ind_000) + Rn_5(Ind_200)*gmj(Ind_100) + &
                Rn_5(Ind_110)*gmj(Ind_010) + Rn_5(Ind_101)*gmj(Ind_001) + &
                Rn_5(Ind_300)*gmj(Ind_200) + Rn_5(Ind_120)*gmj(Ind_020) + &
                Rn_5(Ind_102)*gmj(Ind_002) + Rn_5(Ind_210)*gmj(Ind_110) + &
                Rn_5(Ind_201)*gmj(Ind_101) + Rn_5(Ind_111)*gmj(Ind_011))
            phi(Ind_010) = -(Rn_5(Ind_010)*gmj(Ind_000) + Rn_5(Ind_110)*gmj(Ind_100) + &
                Rn_5(Ind_020)*gmj(Ind_010) + Rn_5(Ind_011)*gmj(Ind_001) + &
                Rn_5(Ind_210)*gmj(Ind_200) + Rn_5(Ind_030)*gmj(Ind_020) + &
                Rn_5(Ind_012)*gmj(Ind_002) + Rn_5(Ind_120)*gmj(Ind_110) + &
                Rn_5(Ind_111)*gmj(Ind_101) + Rn_5(Ind_021)*gmj(Ind_011))
            phi(Ind_001) = -(Rn_5(Ind_001)*gmj(Ind_000) + Rn_5(Ind_101)*gmj(Ind_100) + &
                Rn_5(Ind_011)*gmj(Ind_010) + Rn_5(Ind_002)*gmj(Ind_001) + &
                Rn_5(Ind_201)*gmj(Ind_200) + Rn_5(Ind_021)*gmj(Ind_020) + &
                Rn_5(Ind_003)*gmj(Ind_002) + Rn_5(Ind_111)*gmj(Ind_110) + &
                Rn_5(Ind_102)*gmj(Ind_101) + Rn_5(Ind_012)*gmj(Ind_011))
            phi(Ind_200) = Rn_5(Ind_200)*gmj(Ind_000) + Rn_5(Ind_300)*gmj(Ind_100) + &
                Rn_5(Ind_210)*gmj(Ind_010) + Rn_5(Ind_201)*gmj(Ind_001) + &
                Rn_5(Ind_400)*gmj(Ind_200) + Rn_5(Ind_220)*gmj(Ind_020) + &
                Rn_5(Ind_202)*gmj(Ind_002) + Rn_5(Ind_310)*gmj(Ind_110) + &
                Rn_5(Ind_301)*gmj(Ind_101) + Rn_5(Ind_211)*gmj(Ind_011)
            phi(Ind_020) = Rn_5(Ind_020)*gmj(Ind_000) + Rn_5(Ind_120)*gmj(Ind_100) + &
                Rn_5(Ind_030)*gmj(Ind_010) + Rn_5(Ind_021)*gmj(Ind_001) + &
                Rn_5(Ind_220)*gmj(Ind_200) + Rn_5(Ind_040)*gmj(Ind_020) + &
                Rn_5(Ind_022)*gmj(Ind_002) + Rn_5(Ind_130)*gmj(Ind_110) + &
                Rn_5(Ind_121)*gmj(Ind_101) + Rn_5(Ind_031)*gmj(Ind_011)
            phi(Ind_002) = Rn_5(Ind_002)*gmj(Ind_000) + Rn_5(Ind_102)*gmj(Ind_100) + &
                Rn_5(Ind_012)*gmj(Ind_010) + Rn_5(Ind_003)*gmj(Ind_001) + &
                Rn_5(Ind_202)*gmj(Ind_200) + Rn_5(Ind_022)*gmj(Ind_020) + &
                Rn_5(Ind_004)*gmj(Ind_002) + Rn_5(Ind_112)*gmj(Ind_110) + &
                Rn_5(Ind_103)*gmj(Ind_101) + Rn_5(Ind_013)*gmj(Ind_011)
            phi(Ind_110) = Rn_5(Ind_110)*gmj(Ind_000) + Rn_5(Ind_210)*gmj(Ind_100) + &
                Rn_5(Ind_120)*gmj(Ind_010) + Rn_5(Ind_111)*gmj(Ind_001) + &
                Rn_5(Ind_310)*gmj(Ind_200) + Rn_5(Ind_130)*gmj(Ind_020) + &
                Rn_5(Ind_112)*gmj(Ind_002) + Rn_5(Ind_220)*gmj(Ind_110) + &
                Rn_5(Ind_211)*gmj(Ind_101) + Rn_5(Ind_121)*gmj(Ind_011)
            phi(Ind_101) = Rn_5(Ind_101)*gmj(Ind_000) + Rn_5(Ind_201)*gmj(Ind_100) + &
                Rn_5(Ind_111)*gmj(Ind_010) + Rn_5(Ind_102)*gmj(Ind_001) + &
                Rn_5(Ind_301)*gmj(Ind_200) + Rn_5(Ind_121)*gmj(Ind_020) + &
                Rn_5(Ind_103)*gmj(Ind_002) + Rn_5(Ind_211)*gmj(Ind_110) + &
                Rn_5(Ind_202)*gmj(Ind_101) + Rn_5(Ind_112)*gmj(Ind_011)
            phi(Ind_011) = Rn_5(Ind_011)*gmj(Ind_000) + Rn_5(Ind_111)*gmj(Ind_100) + &
                Rn_5(Ind_021)*gmj(Ind_010) + Rn_5(Ind_012)*gmj(Ind_001) + &
                Rn_5(Ind_211)*gmj(Ind_200) + Rn_5(Ind_031)*gmj(Ind_020) + &
                Rn_5(Ind_013)*gmj(Ind_002) + Rn_5(Ind_121)*gmj(Ind_110) + &
                Rn_5(Ind_112)*gmj(Ind_101) + Rn_5(Ind_022)*gmj(Ind_011)
            phi(Ind_300) = -(Rn_5(Ind_300)*gmj(Ind_000) + Rn_5(Ind_400)*gmj(Ind_100) + &
                Rn_5(Ind_310)*gmj(Ind_010) + Rn_5(Ind_301)*gmj(Ind_001) + &
                Rn_5(Ind_500)*gmj(Ind_200) + Rn_5(Ind_320)*gmj(Ind_020) + &
                Rn_5(Ind_302)*gmj(Ind_002) + Rn_5(Ind_410)*gmj(Ind_110) + &
                Rn_5(Ind_401)*gmj(Ind_101) + Rn_5(Ind_311)*gmj(Ind_011))
            phi(Ind_030) = -(Rn_5(Ind_030)*gmj(Ind_000) + Rn_5(Ind_130)*gmj(Ind_100) + &
                Rn_5(Ind_040)*gmj(Ind_010) + Rn_5(Ind_031)*gmj(Ind_001) + &
                Rn_5(Ind_230)*gmj(Ind_200) + Rn_5(Ind_050)*gmj(Ind_020) + &
                Rn_5(Ind_032)*gmj(Ind_002) + Rn_5(Ind_140)*gmj(Ind_110) + &
                Rn_5(Ind_131)*gmj(Ind_101) + Rn_5(Ind_041)*gmj(Ind_011))
            phi(Ind_003) = -(Rn_5(Ind_003)*gmj(Ind_000) + Rn_5(Ind_103)*gmj(Ind_100) + &
                Rn_5(Ind_013)*gmj(Ind_010) + Rn_5(Ind_004)*gmj(Ind_001) + &
                Rn_5(Ind_203)*gmj(Ind_200) + Rn_5(Ind_023)*gmj(Ind_020) + &
                Rn_5(Ind_005)*gmj(Ind_002) + Rn_5(Ind_113)*gmj(Ind_110) + &
                Rn_5(Ind_104)*gmj(Ind_101) + Rn_5(Ind_014)*gmj(Ind_011))
            phi(Ind_210) = -(Rn_5(Ind_210)*gmj(Ind_000) + Rn_5(Ind_310)*gmj(Ind_100) + &
                Rn_5(Ind_220)*gmj(Ind_010) + Rn_5(Ind_211)*gmj(Ind_001) + &
                Rn_5(Ind_410)*gmj(Ind_200) + Rn_5(Ind_230)*gmj(Ind_020) + &
                Rn_5(Ind_212)*gmj(Ind_002) + Rn_5(Ind_320)*gmj(Ind_110) + &
                Rn_5(Ind_311)*gmj(Ind_101) + Rn_5(Ind_221)*gmj(Ind_011))
            phi(Ind_201) = -(Rn_5(Ind_201)*gmj(Ind_000) + Rn_5(Ind_301)*gmj(Ind_100) + &
                Rn_5(Ind_211)*gmj(Ind_010) + Rn_5(Ind_202)*gmj(Ind_001) + &
                Rn_5(Ind_401)*gmj(Ind_200) + Rn_5(Ind_221)*gmj(Ind_020) + &
                Rn_5(Ind_203)*gmj(Ind_002) + Rn_5(Ind_311)*gmj(Ind_110) + &
                Rn_5(Ind_302)*gmj(Ind_101) + Rn_5(Ind_212)*gmj(Ind_011))
            phi(Ind_120) = -(Rn_5(Ind_120)*gmj(Ind_000) + Rn_5(Ind_220)*gmj(Ind_100) + &
                Rn_5(Ind_130)*gmj(Ind_010) + Rn_5(Ind_121)*gmj(Ind_001) + &
                Rn_5(Ind_320)*gmj(Ind_200) + Rn_5(Ind_140)*gmj(Ind_020) + &
                Rn_5(Ind_122)*gmj(Ind_002) + Rn_5(Ind_230)*gmj(Ind_110) + &
                Rn_5(Ind_221)*gmj(Ind_101) + Rn_5(Ind_131)*gmj(Ind_011))
            phi(Ind_021) = -(Rn_5(Ind_021)*gmj(Ind_000) + Rn_5(Ind_121)*gmj(Ind_100) + &
                Rn_5(Ind_031)*gmj(Ind_010) + Rn_5(Ind_022)*gmj(Ind_001) + &
                Rn_5(Ind_221)*gmj(Ind_200) + Rn_5(Ind_041)*gmj(Ind_020) + &
                Rn_5(Ind_023)*gmj(Ind_002) + Rn_5(Ind_131)*gmj(Ind_110) + &
                Rn_5(Ind_122)*gmj(Ind_101) + Rn_5(Ind_032)*gmj(Ind_011))
            phi(Ind_102) = -(Rn_5(Ind_102)*gmj(Ind_000) + Rn_5(Ind_202)*gmj(Ind_100) + &
                Rn_5(Ind_112)*gmj(Ind_010) + Rn_5(Ind_103)*gmj(Ind_001) + &
                Rn_5(Ind_302)*gmj(Ind_200) + Rn_5(Ind_122)*gmj(Ind_020) + &
                Rn_5(Ind_104)*gmj(Ind_002) + Rn_5(Ind_212)*gmj(Ind_110) + &
                Rn_5(Ind_203)*gmj(Ind_101) + Rn_5(Ind_113)*gmj(Ind_011))
            phi(Ind_012) = -(Rn_5(Ind_012)*gmj(Ind_000) + Rn_5(Ind_112)*gmj(Ind_100) + &
                Rn_5(Ind_022)*gmj(Ind_010) + Rn_5(Ind_013)*gmj(Ind_001) + &
                Rn_5(Ind_212)*gmj(Ind_200) + Rn_5(Ind_032)*gmj(Ind_020) + &
                Rn_5(Ind_014)*gmj(Ind_002) + Rn_5(Ind_122)*gmj(Ind_110) + &
                Rn_5(Ind_113)*gmj(Ind_101) + Rn_5(Ind_023)*gmj(Ind_011))
            phi(Ind_111) = -(Rn_5(Ind_111)*gmj(Ind_000) + Rn_5(Ind_211)*gmj(Ind_100) + &
                Rn_5(Ind_121)*gmj(Ind_010) + Rn_5(Ind_112)*gmj(Ind_001) + &
                Rn_5(Ind_311)*gmj(Ind_200) + Rn_5(Ind_131)*gmj(Ind_020) + &
                Rn_5(Ind_113)*gmj(Ind_002) + Rn_5(Ind_221)*gmj(Ind_110) + &
                Rn_5(Ind_212)*gmj(Ind_101) + Rn_5(Ind_122)*gmj(Ind_011))
            e_pp = phi(Ind_000)*gmi(Ind_000) + phi(Ind_100)*gmi(Ind_100) + &
                phi(Ind_010)*gmi(Ind_010) + phi(Ind_001)*gmi(Ind_001) + &
                phi(Ind_200)*gmi(Ind_200) + phi(Ind_020)*gmi(Ind_020) + &
                phi(Ind_002)*gmi(Ind_002) + phi(Ind_110)*gmi(Ind_110) + &
                phi(Ind_101)*gmi(Ind_101) + phi(Ind_011)*gmi(Ind_011)
            ! gradient of e_pp wrt r_i
            g_pp(1) = phi(Ind_100)*gmi(Ind_000) + phi(Ind_200)*gmi(Ind_100) + &
                phi(Ind_110)*gmi(Ind_010) + phi(Ind_101)*gmi(Ind_001) + &
                phi(Ind_300)*gmi(Ind_200) + phi(Ind_120)*gmi(Ind_020) + &
                phi(Ind_102)*gmi(Ind_002) + phi(Ind_210)*gmi(Ind_110) + &
                phi(Ind_201)*gmi(Ind_101) + phi(Ind_111)*gmi(Ind_011)
            g_pp(2) = phi(Ind_010)*gmi(Ind_000) + phi(Ind_110)*gmi(Ind_100) + &
                phi(Ind_020)*gmi(Ind_010) + phi(Ind_011)*gmi(Ind_001) + &
                phi(Ind_210)*gmi(Ind_200) + phi(Ind_030)*gmi(Ind_020) + &
                phi(Ind_012)*gmi(Ind_002) + phi(Ind_120)*gmi(Ind_110) + &
                phi(Ind_111)*gmi(Ind_101) + phi(Ind_021)*gmi(Ind_011)
            g_pp(3) = phi(Ind_001)*gmi(Ind_000) + phi(Ind_101)*gmi(Ind_100) + &
                phi(Ind_011)*gmi(Ind_010) + phi(Ind_002)*gmi(Ind_001) + &
                phi(Ind_201)*gmi(Ind_200) + phi(Ind_021)*gmi(Ind_020) + &
                phi(Ind_003)*gmi(Ind_002) + phi(Ind_111)*gmi(Ind_110) + &
                phi(Ind_102)*gmi(Ind_101) + phi(Ind_012)*gmi(Ind_011)
            ! torque field at i due to permanent mpoles of j
            do jj = 1, 10
                torque_field(jj, i) = torque_field(jj, i) + phi(jj)
            end do
            ! next do field at j due to the permanent mpoles of i to get torque at j
            ! electrostatic potential at j due to permanent mpoles at i
            ! and derivatives of that with respect to r_j
            ! minus signs due to deriv of r_j-r_i wrt r_i entering into potential
            phi(Ind_000) = Rn_5(Ind_000)*gmi(Ind_000) - Rn_5(Ind_100)*gmi(Ind_100) - &
                Rn_5(Ind_010)*gmi(Ind_010) - Rn_5(Ind_001)*gmi(Ind_001) + &
                Rn_5(Ind_200)*gmi(Ind_200) + Rn_5(Ind_020)*gmi(Ind_020) + &
                Rn_5(Ind_002)*gmi(Ind_002) + Rn_5(Ind_110)*gmi(Ind_110) + &
                Rn_5(Ind_101)*gmi(Ind_101) + Rn_5(Ind_011)*gmi(Ind_011)
            phi(Ind_100) = Rn_5(Ind_100)*gmi(Ind_000) - Rn_5(Ind_200)*gmi(Ind_100) - &
                Rn_5(Ind_110)*gmi(Ind_010) - Rn_5(Ind_101)*gmi(Ind_001) + &
                Rn_5(Ind_300)*gmi(Ind_200) + Rn_5(Ind_120)*gmi(Ind_020) + &
                Rn_5(Ind_102)*gmi(Ind_002) + Rn_5(Ind_210)*gmi(Ind_110) + &
                Rn_5(Ind_201)*gmi(Ind_101) + Rn_5(Ind_111)*gmi(Ind_011)
            phi(Ind_010) = Rn_5(Ind_010)*gmi(Ind_000) - Rn_5(Ind_110)*gmi(Ind_100) - &
                Rn_5(Ind_020)*gmi(Ind_010) - Rn_5(Ind_011)*gmi(Ind_001) + &
                Rn_5(Ind_210)*gmi(Ind_200) + Rn_5(Ind_030)*gmi(Ind_020) + &
                Rn_5(Ind_012)*gmi(Ind_002) + Rn_5(Ind_120)*gmi(Ind_110) + &
                Rn_5(Ind_111)*gmi(Ind_101) + Rn_5(Ind_021)*gmi(Ind_011)
            phi(Ind_001) = Rn_5(Ind_001)*gmi(Ind_000) - Rn_5(Ind_101)*gmi(Ind_100) - &
                Rn_5(Ind_011)*gmi(Ind_010) - Rn_5(Ind_002)*gmi(Ind_001) + &
                Rn_5(Ind_201)*gmi(Ind_200) + Rn_5(Ind_021)*gmi(Ind_020) + &
                Rn_5(Ind_003)*gmi(Ind_002) + Rn_5(Ind_111)*gmi(Ind_110) + &
                Rn_5(Ind_102)*gmi(Ind_101) + Rn_5(Ind_012)*gmi(Ind_011)
            phi(Ind_200) = Rn_5(Ind_200)*gmi(Ind_000) - Rn_5(Ind_300)*gmi(Ind_100) - &
                Rn_5(Ind_210)*gmi(Ind_010) - Rn_5(Ind_201)*gmi(Ind_001) + &
                Rn_5(Ind_400)*gmi(Ind_200) + Rn_5(Ind_220)*gmi(Ind_020) + &
                Rn_5(Ind_202)*gmi(Ind_002) + Rn_5(Ind_310)*gmi(Ind_110) + &
                Rn_5(Ind_301)*gmi(Ind_101) + Rn_5(Ind_211)*gmi(Ind_011)
            phi(Ind_020) = Rn_5(Ind_020)*gmi(Ind_000) - Rn_5(Ind_120)*gmi(Ind_100) - &
                Rn_5(Ind_030)*gmi(Ind_010) - Rn_5(Ind_021)*gmi(Ind_001) + &
                Rn_5(Ind_220)*gmi(Ind_200) + Rn_5(Ind_040)*gmi(Ind_020) + &
                Rn_5(Ind_022)*gmi(Ind_002) + Rn_5(Ind_130)*gmi(Ind_110) + &
                Rn_5(Ind_121)*gmi(Ind_101) + Rn_5(Ind_031)*gmi(Ind_011)
            phi(Ind_002) = Rn_5(Ind_002)*gmi(Ind_000) - Rn_5(Ind_102)*gmi(Ind_100) - &
                Rn_5(Ind_012)*gmi(Ind_010) - Rn_5(Ind_003)*gmi(Ind_001) + &
                Rn_5(Ind_202)*gmi(Ind_200) + Rn_5(Ind_022)*gmi(Ind_020) + &
                Rn_5(Ind_004)*gmi(Ind_002) + Rn_5(Ind_112)*gmi(Ind_110) + &
                Rn_5(Ind_103)*gmi(Ind_101) + Rn_5(Ind_013)*gmi(Ind_011)
            phi(Ind_110) = Rn_5(Ind_110)*gmi(Ind_000) - Rn_5(Ind_210)*gmi(Ind_100) - &
                Rn_5(Ind_120)*gmi(Ind_010) - Rn_5(Ind_111)*gmi(Ind_001) + &
                Rn_5(Ind_310)*gmi(Ind_200) + Rn_5(Ind_130)*gmi(Ind_020) + &
                Rn_5(Ind_112)*gmi(Ind_002) + Rn_5(Ind_220)*gmi(Ind_110) + &
                Rn_5(Ind_211)*gmi(Ind_101) + Rn_5(Ind_121)*gmi(Ind_011)
            phi(Ind_101) = Rn_5(Ind_101)*gmi(Ind_000) - Rn_5(Ind_201)*gmi(Ind_100) - &
                Rn_5(Ind_111)*gmi(Ind_010) - Rn_5(Ind_102)*gmi(Ind_001) + &
                Rn_5(Ind_301)*gmi(Ind_200) + Rn_5(Ind_121)*gmi(Ind_020) + &
                Rn_5(Ind_103)*gmi(Ind_002) + Rn_5(Ind_211)*gmi(Ind_110) + &
                Rn_5(Ind_202)*gmi(Ind_101) + Rn_5(Ind_112)*gmi(Ind_011)
            phi(Ind_011) = Rn_5(Ind_011)*gmi(Ind_000) - Rn_5(Ind_111)*gmi(Ind_100) - &
                Rn_5(Ind_021)*gmi(Ind_010) - Rn_5(Ind_012)*gmi(Ind_001) + &
                Rn_5(Ind_211)*gmi(Ind_200) + Rn_5(Ind_031)*gmi(Ind_020) + &
                Rn_5(Ind_013)*gmi(Ind_002) + Rn_5(Ind_121)*gmi(Ind_110) + &
                Rn_5(Ind_112)*gmi(Ind_101) + Rn_5(Ind_022)*gmi(Ind_011)
            ! torque field at j due to permanent mpoles of i
            do jj = 1, 10
                torque_field(jj, j) = torque_field(jj, j) + phi(jj)
            end do
            e_ind = 0.d0
            g_ind(1) = 0.d0
            g_ind(2) = 0.d0
            g_ind(3) = 0.d0
            if (is_polarizable(i) .or. is_polarizable(j)) then
                do jj = 1, 3
                    i_di(jj) = ind_dip_d(jj, i)
                    i_dj(jj) = ind_dip_d(jj, j)
                    i_pi(jj) = ind_dip_p(jj, i)
                    i_pj(jj) = ind_dip_p(jj, j)
                end do
                ! first the direct dipoles interact with polar field
                ! recall A gives coulomb, BD damped and B the erfc*1/r contributions
                do jj = 1, 4
                    C(jj) = (B(jj) - A(jj)) + polar_weight(k)*(A(jj) - BD(jj))
                end do
                ! negate the odd order to get sign right
                C(1) = -C(1)
                C(3) = -C(3)
                ! get the interaction tensor
                n = 4
                Rn(Ind_000) = C(n)
                Rn_1(Ind_000) = C(n - 1)
                Rn_1(Ind_100) = delx*Rn(Ind_000)
                Rn_1(Ind_010) = dely*Rn(Ind_000)
                Rn_1(Ind_001) = delz*Rn(Ind_000)
                Rn_2(Ind_000) = C(n - 2)
                Rn_2(Ind_100) = delx*Rn_1(Ind_000)
                Rn_2(Ind_010) = dely*Rn_1(Ind_000)
                Rn_2(Ind_001) = delz*Rn_1(Ind_000)
                Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
                Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
                Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
                Rn_2(Ind_110) = delx*Rn_1(Ind_010)
                Rn_2(Ind_101) = delx*Rn_1(Ind_001)
                Rn_2(Ind_011) = dely*Rn_1(Ind_001)
                Rn_3(Ind_000) = C(n - 3)
                Rn_3(Ind_100) = delx*Rn_2(Ind_000)
                Rn_3(Ind_010) = dely*Rn_2(Ind_000)
                Rn_3(Ind_001) = delz*Rn_2(Ind_000)
                Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
                Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
                Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
                Rn_3(Ind_110) = delx*Rn_2(Ind_010)
                Rn_3(Ind_101) = delx*Rn_2(Ind_001)
                Rn_3(Ind_011) = dely*Rn_2(Ind_001)
                Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
                Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
                Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
                Rn_3(Ind_210) = dely*Rn_2(Ind_200)
                Rn_3(Ind_201) = delz*Rn_2(Ind_200)
                Rn_3(Ind_120) = delx*Rn_2(Ind_020)
                Rn_3(Ind_021) = delz*Rn_2(Ind_020)
                Rn_3(Ind_102) = delx*Rn_2(Ind_002)
                Rn_3(Ind_012) = dely*Rn_2(Ind_002)
                Rn_3(Ind_111) = delx*Rn_2(Ind_011)
                !Rn_4(Ind_000) = C(n-4) NOT NEEDED
                Rn_4(Ind_100) = delx*Rn_3(Ind_000)
                Rn_4(Ind_010) = dely*Rn_3(Ind_000)
                Rn_4(Ind_001) = delz*Rn_3(Ind_000)
                Rn_4(Ind_200) = Rn_3(Ind_000) + delx*Rn_3(Ind_100)
                Rn_4(Ind_020) = Rn_3(Ind_000) + dely*Rn_3(Ind_010)
                Rn_4(Ind_002) = Rn_3(Ind_000) + delz*Rn_3(Ind_001)
                Rn_4(Ind_110) = delx*Rn_3(Ind_010)
                Rn_4(Ind_101) = delx*Rn_3(Ind_001)
                Rn_4(Ind_011) = dely*Rn_3(Ind_001)
                Rn_4(Ind_300) = 2.d0*Rn_3(Ind_100) + delx*Rn_3(Ind_200)
                Rn_4(Ind_030) = 2.d0*Rn_3(Ind_010) + dely*Rn_3(Ind_020)
                Rn_4(Ind_003) = 2.d0*Rn_3(Ind_001) + delz*Rn_3(Ind_002)
                Rn_4(Ind_210) = dely*Rn_3(Ind_200)
                Rn_4(Ind_201) = delz*Rn_3(Ind_200)
                Rn_4(Ind_120) = delx*Rn_3(Ind_020)
                Rn_4(Ind_021) = delz*Rn_3(Ind_020)
                Rn_4(Ind_102) = delx*Rn_3(Ind_002)
                Rn_4(Ind_012) = dely*Rn_3(Ind_002)
                Rn_4(Ind_111) = delx*Rn_3(Ind_011)
                Rn_4(Ind_400) = 3.d0*Rn_3(Ind_200) + delx*Rn_3(Ind_300)
                Rn_4(Ind_040) = 3.d0*Rn_3(Ind_020) + dely*Rn_3(Ind_030)
                Rn_4(Ind_004) = 3.d0*Rn_3(Ind_002) + delz*Rn_3(Ind_003)
                Rn_4(Ind_310) = dely*Rn_3(Ind_300)
                Rn_4(Ind_301) = delz*Rn_3(Ind_300)
                Rn_4(Ind_130) = delx*Rn_3(Ind_030)
                Rn_4(Ind_031) = delz*Rn_3(Ind_030)
                Rn_4(Ind_103) = delx*Rn_3(Ind_003)
                Rn_4(Ind_013) = dely*Rn_3(Ind_003)
                Rn_4(Ind_220) = Rn_3(Ind_020) + delx*Rn_3(Ind_120)
                Rn_4(Ind_202) = Rn_3(Ind_002) + delx*Rn_3(Ind_102)
                Rn_4(Ind_022) = Rn_3(Ind_002) + dely*Rn_3(Ind_012)
                Rn_4(Ind_211) = dely*Rn_3(Ind_201)
                Rn_4(Ind_121) = delx*Rn_3(Ind_021)
                Rn_4(Ind_112) = delx*Rn_3(Ind_012)
                if (is_polarizable(i)) then
                    ! induced direct dipoles at i interact with permanent at j
                    phi(Ind_100) = &
                        -(Rn_4(Ind_100)*gmj(Ind_000) + Rn_4(Ind_200)*gmj(Ind_100) + &
                        Rn_4(Ind_110)*gmj(Ind_010) + Rn_4(Ind_101)*gmj(Ind_001) + &
                        Rn_4(Ind_300)*gmj(Ind_200) + Rn_4(Ind_120)*gmj(Ind_020) + &
                        Rn_4(Ind_102)*gmj(Ind_002) + Rn_4(Ind_210)*gmj(Ind_110) + &
                        Rn_4(Ind_201)*gmj(Ind_101) + Rn_4(Ind_111)*gmj(Ind_011))
                    phi(Ind_010) = &
                        -(Rn_4(Ind_010)*gmj(Ind_000) + Rn_4(Ind_110)*gmj(Ind_100) + &
                        Rn_4(Ind_020)*gmj(Ind_010) + Rn_4(Ind_011)*gmj(Ind_001) + &
                        Rn_4(Ind_210)*gmj(Ind_200) + Rn_4(Ind_030)*gmj(Ind_020) + &
                        Rn_4(Ind_012)*gmj(Ind_002) + Rn_4(Ind_120)*gmj(Ind_110) + &
                        Rn_4(Ind_111)*gmj(Ind_101) + Rn_4(Ind_021)*gmj(Ind_011))
                    phi(Ind_001) = &
                        -(Rn_4(Ind_001)*gmj(Ind_000) + Rn_4(Ind_101)*gmj(Ind_100) + &
                        Rn_4(Ind_011)*gmj(Ind_010) + Rn_4(Ind_002)*gmj(Ind_001) + &
                        Rn_4(Ind_201)*gmj(Ind_200) + Rn_4(Ind_021)*gmj(Ind_020) + &
                        Rn_4(Ind_003)*gmj(Ind_002) + Rn_4(Ind_111)*gmj(Ind_110) + &
                        Rn_4(Ind_102)*gmj(Ind_101) + Rn_4(Ind_012)*gmj(Ind_011))
                    phi(Ind_200) = &
                        Rn_4(Ind_200)*gmj(Ind_000) + Rn_4(Ind_300)*gmj(Ind_100) + &
                        Rn_4(Ind_210)*gmj(Ind_010) + Rn_4(Ind_201)*gmj(Ind_001) + &
                        Rn_4(Ind_400)*gmj(Ind_200) + Rn_4(Ind_220)*gmj(Ind_020) + &
                        Rn_4(Ind_202)*gmj(Ind_002) + Rn_4(Ind_310)*gmj(Ind_110) + &
                        Rn_4(Ind_301)*gmj(Ind_101) + Rn_4(Ind_211)*gmj(Ind_011)
                    phi(Ind_020) = &
                        Rn_4(Ind_020)*gmj(Ind_000) + Rn_4(Ind_120)*gmj(Ind_100) + &
                        Rn_4(Ind_030)*gmj(Ind_010) + Rn_4(Ind_021)*gmj(Ind_001) + &
                        Rn_4(Ind_220)*gmj(Ind_200) + Rn_4(Ind_040)*gmj(Ind_020) + &
                        Rn_4(Ind_022)*gmj(Ind_002) + Rn_4(Ind_130)*gmj(Ind_110) + &
                        Rn_4(Ind_121)*gmj(Ind_101) + Rn_4(Ind_031)*gmj(Ind_011)
                    phi(Ind_002) = &
                        Rn_4(Ind_002)*gmj(Ind_000) + Rn_4(Ind_102)*gmj(Ind_100) + &
                        Rn_4(Ind_012)*gmj(Ind_010) + Rn_4(Ind_003)*gmj(Ind_001) + &
                        Rn_4(Ind_202)*gmj(Ind_200) + Rn_4(Ind_022)*gmj(Ind_020) + &
                        Rn_4(Ind_004)*gmj(Ind_002) + Rn_4(Ind_112)*gmj(Ind_110) + &
                        Rn_4(Ind_103)*gmj(Ind_101) + Rn_4(Ind_013)*gmj(Ind_011)
                    phi(Ind_110) = &
                        Rn_4(Ind_110)*gmj(Ind_000) + Rn_4(Ind_210)*gmj(Ind_100) + &
                        Rn_4(Ind_120)*gmj(Ind_010) + Rn_4(Ind_111)*gmj(Ind_001) + &
                        Rn_4(Ind_310)*gmj(Ind_200) + Rn_4(Ind_130)*gmj(Ind_020) + &
                        Rn_4(Ind_112)*gmj(Ind_002) + Rn_4(Ind_220)*gmj(Ind_110) + &
                        Rn_4(Ind_211)*gmj(Ind_101) + Rn_4(Ind_121)*gmj(Ind_011)
                    phi(Ind_101) = &
                        Rn_4(Ind_101)*gmj(Ind_000) + Rn_4(Ind_201)*gmj(Ind_100) + &
                        Rn_4(Ind_111)*gmj(Ind_010) + Rn_4(Ind_102)*gmj(Ind_001) + &
                        Rn_4(Ind_301)*gmj(Ind_200) + Rn_4(Ind_121)*gmj(Ind_020) + &
                        Rn_4(Ind_103)*gmj(Ind_002) + Rn_4(Ind_211)*gmj(Ind_110) + &
                        Rn_4(Ind_202)*gmj(Ind_101) + Rn_4(Ind_112)*gmj(Ind_011)
                    phi(Ind_011) = &
                        Rn_4(Ind_011)*gmj(Ind_000) + Rn_4(Ind_111)*gmj(Ind_100) + &
                        Rn_4(Ind_021)*gmj(Ind_010) + Rn_4(Ind_012)*gmj(Ind_001) + &
                        Rn_4(Ind_211)*gmj(Ind_200) + Rn_4(Ind_031)*gmj(Ind_020) + &
                        Rn_4(Ind_013)*gmj(Ind_002) + Rn_4(Ind_121)*gmj(Ind_110) + &
                        Rn_4(Ind_112)*gmj(Ind_101) + Rn_4(Ind_022)*gmj(Ind_011)
                    e_ind = e_ind + 0.5d0* &
                        (phi(Ind_100)*i_di(1) + phi(Ind_010)*i_di(2) + &
                        phi(Ind_001)*i_di(3))
                    g_ind(1) = g_ind(1) + 0.5d0* &
                        (phi(Ind_200)*i_di(1) + phi(Ind_110)*i_di(2) + &
                        phi(Ind_101)*i_di(3))
                    g_ind(2) = g_ind(2) + 0.5d0* &
                        (phi(Ind_110)*i_di(1) + phi(Ind_020)*i_di(2) + &
                        phi(Ind_011)*i_di(3))
                    g_ind(3) = g_ind(3) + 0.5d0* &
                        (phi(Ind_101)*i_di(1) + phi(Ind_011)*i_di(2) + &
                        phi(Ind_002)*i_di(3))
                    ! calculate the potential at j due to direct induced at i
                    ! and derivatives of that with respect to r_j
                    ! minus signs due to deriv of r_j-r_i wrt r_i entering into potential
                    phi(Ind_000) = -(Rn_4(Ind_100)*i_di(1) + Rn_4(Ind_010)*i_di(2) + &
                        Rn_4(Ind_001)*i_di(3))
                    phi(Ind_100) = -(Rn_4(Ind_200)*i_di(1) + Rn_4(Ind_110)*i_di(2) + &
                        Rn_4(Ind_101)*i_di(3))
                    phi(Ind_010) = -(Rn_4(Ind_110)*i_di(1) + Rn_4(Ind_020)*i_di(2) + &
                        Rn_4(Ind_011)*i_di(3))
                    phi(Ind_001) = -(Rn_4(Ind_101)*i_di(1) + Rn_4(Ind_011)*i_di(2) + &
                        Rn_4(Ind_002)*i_di(3))
                    phi(Ind_200) = -(Rn_4(Ind_300)*i_di(1) + Rn_4(Ind_210)*i_di(2) + &
                        Rn_4(Ind_201)*i_di(3))
                    phi(Ind_020) = -(Rn_4(Ind_120)*i_di(1) + Rn_4(Ind_030)*i_di(2) + &
                        Rn_4(Ind_021)*i_di(3))
                    phi(Ind_002) = -(Rn_4(Ind_102)*i_di(1) + Rn_4(Ind_012)*i_di(2) + &
                        Rn_4(Ind_003)*i_di(3))
                    phi(Ind_110) = -(Rn_4(Ind_210)*i_di(1) + Rn_4(Ind_120)*i_di(2) + &
                        Rn_4(Ind_111)*i_di(3))
                    phi(Ind_101) = -(Rn_4(Ind_201)*i_di(1) + Rn_4(Ind_111)*i_di(2) + &
                        Rn_4(Ind_102)*i_di(3))
                    phi(Ind_011) = -(Rn_4(Ind_111)*i_di(1) + Rn_4(Ind_021)*i_di(2) + &
                        Rn_4(Ind_012)*i_di(3))
                    ! torque field at j due to direct induced mpoles of i
                    ! note factor of 1/2
                    do jj = 1, 10
                        torque_field(jj, j) = torque_field(jj, j) + 0.5d0*phi(jj)
                    end do
                end if !is_polarizable(i)
                if (is_polarizable(j)) then
                    ! induced direct dipoles at j interact with permanent at i
                    ! phi array (electrostatic potential at i due to j induced moments
                    ! and derivs of that esp wrt r_i )
                    ! first that due to ind_dip_d for energy contribution
                    ! minus signs arise due to derivs of r_j - r_i wrt r_i
                    phi(Ind_000) = Rn_4(Ind_100)*i_dj(1) + Rn_4(Ind_010)*i_dj(2) + &
                        Rn_4(Ind_001)*i_dj(3)
                    phi(Ind_100) = -(Rn_4(Ind_200)*i_dj(1) + Rn_4(Ind_110)*i_dj(2) + &
                        Rn_4(Ind_101)*i_dj(3))
                    phi(Ind_010) = -(Rn_4(Ind_110)*i_dj(1) + Rn_4(Ind_020)*i_dj(2) + &
                        Rn_4(Ind_011)*i_dj(3))
                    phi(Ind_001) = -(Rn_4(Ind_101)*i_dj(1) + Rn_4(Ind_011)*i_dj(2) + &
                        Rn_4(Ind_002)*i_dj(3))
                    phi(Ind_200) = Rn_4(Ind_300)*i_dj(1) + Rn_4(Ind_210)*i_dj(2) + &
                        Rn_4(Ind_201)*i_dj(3)
                    phi(Ind_020) = Rn_4(Ind_120)*i_dj(1) + Rn_4(Ind_030)*i_dj(2) + &
                        Rn_4(Ind_021)*i_dj(3)
                    phi(Ind_002) = Rn_4(Ind_102)*i_dj(1) + Rn_4(Ind_012)*i_dj(2) + &
                        Rn_4(Ind_003)*i_dj(3)
                    phi(Ind_110) = Rn_4(Ind_210)*i_dj(1) + Rn_4(Ind_120)*i_dj(2) + &
                        Rn_4(Ind_111)*i_dj(3)
                    phi(Ind_101) = Rn_4(Ind_201)*i_dj(1) + Rn_4(Ind_111)*i_dj(2) + &
                        Rn_4(Ind_102)*i_dj(3)
                    phi(Ind_011) = Rn_4(Ind_111)*i_dj(1) + Rn_4(Ind_021)*i_dj(2) + &
                        Rn_4(Ind_012)*i_dj(3)
                    phi(Ind_300) = -(Rn_4(Ind_400)*i_dj(1) + Rn_4(Ind_310)*i_dj(2) + &
                        Rn_4(Ind_301)*i_dj(3))
                    phi(Ind_030) = -(Rn_4(Ind_130)*i_dj(1) + Rn_4(Ind_040)*i_dj(2) + &
                        Rn_4(Ind_031)*i_dj(3))
                    phi(Ind_003) = -(Rn_4(Ind_103)*i_dj(1) + Rn_4(Ind_013)*i_dj(2) + &
                        Rn_4(Ind_004)*i_dj(3))
                    phi(Ind_210) = -(Rn_4(Ind_310)*i_dj(1) + Rn_4(Ind_220)*i_dj(2) + &
                        Rn_4(Ind_211)*i_dj(3))
                    phi(Ind_201) = -(Rn_4(Ind_301)*i_dj(1) + Rn_4(Ind_211)*i_dj(2) + &
                        Rn_4(Ind_202)*i_dj(3))
                    phi(Ind_120) = -(Rn_4(Ind_220)*i_dj(1) + Rn_4(Ind_130)*i_dj(2) + &
                        Rn_4(Ind_121)*i_dj(3))
                    phi(Ind_021) = -(Rn_4(Ind_121)*i_dj(1) + Rn_4(Ind_031)*i_dj(2) + &
                        Rn_4(Ind_022)*i_dj(3))
                    phi(Ind_102) = -(Rn_4(Ind_202)*i_dj(1) + Rn_4(Ind_112)*i_dj(2) + &
                        Rn_4(Ind_103)*i_dj(3))
                    phi(Ind_012) = -(Rn_4(Ind_112)*i_dj(1) + Rn_4(Ind_022)*i_dj(2) + &
                        Rn_4(Ind_013)*i_dj(3))
                    phi(Ind_111) = -(Rn_4(Ind_211)*i_dj(1) + Rn_4(Ind_121)*i_dj(2) + &
                        Rn_4(Ind_112)*i_dj(3))
                    e_ind = e_ind + 0.5d0* &
                        (phi(Ind_000)*gmi(Ind_000) + phi(Ind_100)*gmi(Ind_100) + &
                        phi(Ind_010)*gmi(Ind_010) + phi(Ind_001)*gmi(Ind_001) + &
                        phi(Ind_200)*gmi(Ind_200) + phi(Ind_020)*gmi(Ind_020) + &
                        phi(Ind_002)*gmi(Ind_002) + phi(Ind_110)*gmi(Ind_110) + &
                        phi(Ind_101)*gmi(Ind_101) + phi(Ind_011)*gmi(Ind_011))
                    g_ind(1) = g_ind(1) + 0.5d0* &
                        (phi(Ind_100)*gmi(Ind_000) + phi(Ind_200)*gmi(Ind_100) + &
                        phi(Ind_110)*gmi(Ind_010) + phi(Ind_101)*gmi(Ind_001) + &
                        phi(Ind_300)*gmi(Ind_200) + phi(Ind_120)*gmi(Ind_020) + &
                        phi(Ind_102)*gmi(Ind_002) + phi(Ind_210)*gmi(Ind_110) + &
                        phi(Ind_201)*gmi(Ind_101) + phi(Ind_111)*gmi(Ind_011))
                    g_ind(2) = g_ind(2) + 0.5d0* &
                        (phi(Ind_010)*gmi(Ind_000) + phi(Ind_110)*gmi(Ind_100) + &
                        phi(Ind_020)*gmi(Ind_010) + phi(Ind_011)*gmi(Ind_001) + &
                        phi(Ind_210)*gmi(Ind_200) + phi(Ind_030)*gmi(Ind_020) + &
                        phi(Ind_012)*gmi(Ind_002) + phi(Ind_120)*gmi(Ind_110) + &
                        phi(Ind_111)*gmi(Ind_101) + phi(Ind_021)*gmi(Ind_011))
                    g_ind(3) = g_ind(3) + 0.5d0* &
                        (phi(Ind_001)*gmi(Ind_000) + phi(Ind_101)*gmi(Ind_100) + &
                        phi(Ind_011)*gmi(Ind_010) + phi(Ind_002)*gmi(Ind_001) + &
                        phi(Ind_201)*gmi(Ind_200) + phi(Ind_021)*gmi(Ind_020) + &
                        phi(Ind_003)*gmi(Ind_002) + phi(Ind_111)*gmi(Ind_110) + &
                        phi(Ind_102)*gmi(Ind_101) + phi(Ind_012)*gmi(Ind_011))
                    ! torque field at i due to direct dipoles of j
                    ! note factor of 1/2
                    do jj = 1, 10
                        torque_field(jj, i) = torque_field(jj, i) + 0.5d0*phi(jj)
                    end do
                end if ! is_polarizable(j) )then
                ! next the polar dipoles interact with direct field in grad
                ! recall A gives coulomb, BD damped and B the erfc*1/r contributions
                do jj = 1, 4
                    C(jj) = (B(jj) - A(jj)) + direct_weight(k)*(A(jj) - BD(jj))
                end do
                ! negate the odd order to get sign right
                C(1) = -C(1)
                C(3) = -C(3)
                ! get the interaction tensor
                n = 4
                Rn(Ind_000) = C(n)
                Rn_1(Ind_000) = C(n - 1)
                Rn_1(Ind_100) = delx*Rn(Ind_000)
                Rn_1(Ind_010) = dely*Rn(Ind_000)
                Rn_1(Ind_001) = delz*Rn(Ind_000)
                Rn_2(Ind_000) = C(n - 2)
                Rn_2(Ind_100) = delx*Rn_1(Ind_000)
                Rn_2(Ind_010) = dely*Rn_1(Ind_000)
                Rn_2(Ind_001) = delz*Rn_1(Ind_000)
                Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
                Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
                Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
                Rn_2(Ind_110) = delx*Rn_1(Ind_010)
                Rn_2(Ind_101) = delx*Rn_1(Ind_001)
                Rn_2(Ind_011) = dely*Rn_1(Ind_001)
                Rn_3(Ind_000) = C(n - 3)
                Rn_3(Ind_100) = delx*Rn_2(Ind_000)
                Rn_3(Ind_010) = dely*Rn_2(Ind_000)
                Rn_3(Ind_001) = delz*Rn_2(Ind_000)
                Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
                Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
                Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
                Rn_3(Ind_110) = delx*Rn_2(Ind_010)
                Rn_3(Ind_101) = delx*Rn_2(Ind_001)
                Rn_3(Ind_011) = dely*Rn_2(Ind_001)
                Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
                Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
                Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
                Rn_3(Ind_210) = dely*Rn_2(Ind_200)
                Rn_3(Ind_201) = delz*Rn_2(Ind_200)
                Rn_3(Ind_120) = delx*Rn_2(Ind_020)
                Rn_3(Ind_021) = delz*Rn_2(Ind_020)
                Rn_3(Ind_102) = delx*Rn_2(Ind_002)
                Rn_3(Ind_012) = dely*Rn_2(Ind_002)
                Rn_3(Ind_111) = delx*Rn_2(Ind_011)
                Rn_4(Ind_100) = delx*Rn_3(Ind_000)
                Rn_4(Ind_010) = dely*Rn_3(Ind_000)
                Rn_4(Ind_001) = delz*Rn_3(Ind_000)
                Rn_4(Ind_200) = Rn_3(Ind_000) + delx*Rn_3(Ind_100)
                Rn_4(Ind_020) = Rn_3(Ind_000) + dely*Rn_3(Ind_010)
                Rn_4(Ind_002) = Rn_3(Ind_000) + delz*Rn_3(Ind_001)
                Rn_4(Ind_110) = delx*Rn_3(Ind_010)
                Rn_4(Ind_101) = delx*Rn_3(Ind_001)
                Rn_4(Ind_011) = dely*Rn_3(Ind_001)
                Rn_4(Ind_300) = 2.d0*Rn_3(Ind_100) + delx*Rn_3(Ind_200)
                Rn_4(Ind_030) = 2.d0*Rn_3(Ind_010) + dely*Rn_3(Ind_020)
                Rn_4(Ind_003) = 2.d0*Rn_3(Ind_001) + delz*Rn_3(Ind_002)
                Rn_4(Ind_210) = dely*Rn_3(Ind_200)
                Rn_4(Ind_201) = delz*Rn_3(Ind_200)
                Rn_4(Ind_120) = delx*Rn_3(Ind_020)
                Rn_4(Ind_021) = delz*Rn_3(Ind_020)
                Rn_4(Ind_102) = delx*Rn_3(Ind_002)
                Rn_4(Ind_012) = dely*Rn_3(Ind_002)
                Rn_4(Ind_111) = delx*Rn_3(Ind_011)
                Rn_4(Ind_400) = 3.d0*Rn_3(Ind_200) + delx*Rn_3(Ind_300)
                Rn_4(Ind_040) = 3.d0*Rn_3(Ind_020) + dely*Rn_3(Ind_030)
                Rn_4(Ind_004) = 3.d0*Rn_3(Ind_002) + delz*Rn_3(Ind_003)
                Rn_4(Ind_310) = dely*Rn_3(Ind_300)
                Rn_4(Ind_301) = delz*Rn_3(Ind_300)
                Rn_4(Ind_130) = delx*Rn_3(Ind_030)
                Rn_4(Ind_031) = delz*Rn_3(Ind_030)
                Rn_4(Ind_103) = delx*Rn_3(Ind_003)
                Rn_4(Ind_013) = dely*Rn_3(Ind_003)
                Rn_4(Ind_220) = Rn_3(Ind_020) + delx*Rn_3(Ind_120)
                Rn_4(Ind_202) = Rn_3(Ind_002) + delx*Rn_3(Ind_102)
                Rn_4(Ind_022) = Rn_3(Ind_002) + dely*Rn_3(Ind_012)
                Rn_4(Ind_211) = dely*Rn_3(Ind_201)
                Rn_4(Ind_121) = delx*Rn_3(Ind_021)
                Rn_4(Ind_112) = delx*Rn_3(Ind_012)
                if (is_polarizable(i)) then
                    ! induced direct dipoles at i interact with permanent at j
                    phi(Ind_100) = &
                        -(Rn_4(Ind_100)*gmj(Ind_000) + Rn_4(Ind_200)*gmj(Ind_100) + &
                        Rn_4(Ind_110)*gmj(Ind_010) + Rn_4(Ind_101)*gmj(Ind_001) + &
                        Rn_4(Ind_300)*gmj(Ind_200) + Rn_4(Ind_120)*gmj(Ind_020) + &
                        Rn_4(Ind_102)*gmj(Ind_002) + Rn_4(Ind_210)*gmj(Ind_110) + &
                        Rn_4(Ind_201)*gmj(Ind_101) + Rn_4(Ind_111)*gmj(Ind_011))
                    phi(Ind_010) = &
                        -(Rn_4(Ind_010)*gmj(Ind_000) + Rn_4(Ind_110)*gmj(Ind_100) + &
                        Rn_4(Ind_020)*gmj(Ind_010) + Rn_4(Ind_011)*gmj(Ind_001) + &
                        Rn_4(Ind_210)*gmj(Ind_200) + Rn_4(Ind_030)*gmj(Ind_020) + &
                        Rn_4(Ind_012)*gmj(Ind_002) + Rn_4(Ind_120)*gmj(Ind_110) + &
                        Rn_4(Ind_111)*gmj(Ind_101) + Rn_4(Ind_021)*gmj(Ind_011))
                    phi(Ind_001) = &
                        -(Rn_4(Ind_001)*gmj(Ind_000) + Rn_4(Ind_101)*gmj(Ind_100) + &
                        Rn_4(Ind_011)*gmj(Ind_010) + Rn_4(Ind_002)*gmj(Ind_001) + &
                        Rn_4(Ind_201)*gmj(Ind_200) + Rn_4(Ind_021)*gmj(Ind_020) + &
                        Rn_4(Ind_003)*gmj(Ind_002) + Rn_4(Ind_111)*gmj(Ind_110) + &
                        Rn_4(Ind_102)*gmj(Ind_101) + Rn_4(Ind_012)*gmj(Ind_011))
                    phi(Ind_200) = &
                        Rn_4(Ind_200)*gmj(Ind_000) + Rn_4(Ind_300)*gmj(Ind_100) + &
                        Rn_4(Ind_210)*gmj(Ind_010) + Rn_4(Ind_201)*gmj(Ind_001) + &
                        Rn_4(Ind_400)*gmj(Ind_200) + Rn_4(Ind_220)*gmj(Ind_020) + &
                        Rn_4(Ind_202)*gmj(Ind_002) + Rn_4(Ind_310)*gmj(Ind_110) + &
                        Rn_4(Ind_301)*gmj(Ind_101) + Rn_4(Ind_211)*gmj(Ind_011)
                    phi(Ind_020) = &
                        Rn_4(Ind_020)*gmj(Ind_000) + Rn_4(Ind_120)*gmj(Ind_100) + &
                        Rn_4(Ind_030)*gmj(Ind_010) + Rn_4(Ind_021)*gmj(Ind_001) + &
                        Rn_4(Ind_220)*gmj(Ind_200) + Rn_4(Ind_040)*gmj(Ind_020) + &
                        Rn_4(Ind_022)*gmj(Ind_002) + Rn_4(Ind_130)*gmj(Ind_110) + &
                        Rn_4(Ind_121)*gmj(Ind_101) + Rn_4(Ind_031)*gmj(Ind_011)
                    phi(Ind_002) = &
                        Rn_4(Ind_002)*gmj(Ind_000) + Rn_4(Ind_102)*gmj(Ind_100) + &
                        Rn_4(Ind_012)*gmj(Ind_010) + Rn_4(Ind_003)*gmj(Ind_001) + &
                        Rn_4(Ind_202)*gmj(Ind_200) + Rn_4(Ind_022)*gmj(Ind_020) + &
                        Rn_4(Ind_004)*gmj(Ind_002) + Rn_4(Ind_112)*gmj(Ind_110) + &
                        Rn_4(Ind_103)*gmj(Ind_101) + Rn_4(Ind_013)*gmj(Ind_011)
                    phi(Ind_110) = &
                        Rn_4(Ind_110)*gmj(Ind_000) + Rn_4(Ind_210)*gmj(Ind_100) + &
                        Rn_4(Ind_120)*gmj(Ind_010) + Rn_4(Ind_111)*gmj(Ind_001) + &
                        Rn_4(Ind_310)*gmj(Ind_200) + Rn_4(Ind_130)*gmj(Ind_020) + &
                        Rn_4(Ind_112)*gmj(Ind_002) + Rn_4(Ind_220)*gmj(Ind_110) + &
                        Rn_4(Ind_211)*gmj(Ind_101) + Rn_4(Ind_121)*gmj(Ind_011)
                    phi(Ind_101) = &
                        Rn_4(Ind_101)*gmj(Ind_000) + Rn_4(Ind_201)*gmj(Ind_100) + &
                        Rn_4(Ind_111)*gmj(Ind_010) + Rn_4(Ind_102)*gmj(Ind_001) + &
                        Rn_4(Ind_301)*gmj(Ind_200) + Rn_4(Ind_121)*gmj(Ind_020) + &
                        Rn_4(Ind_103)*gmj(Ind_002) + Rn_4(Ind_211)*gmj(Ind_110) + &
                        Rn_4(Ind_202)*gmj(Ind_101) + Rn_4(Ind_112)*gmj(Ind_011)
                    phi(Ind_011) = &
                        Rn_4(Ind_011)*gmj(Ind_000) + Rn_4(Ind_111)*gmj(Ind_100) + &
                        Rn_4(Ind_021)*gmj(Ind_010) + Rn_4(Ind_012)*gmj(Ind_001) + &
                        Rn_4(Ind_211)*gmj(Ind_200) + Rn_4(Ind_031)*gmj(Ind_020) + &
                        Rn_4(Ind_013)*gmj(Ind_002) + Rn_4(Ind_121)*gmj(Ind_110) + &
                        Rn_4(Ind_112)*gmj(Ind_101) + Rn_4(Ind_022)*gmj(Ind_011)
                    g_ind(1) = g_ind(1) + 0.5d0* &
                        (phi(Ind_200)*i_pi(1) + phi(Ind_110)*i_pi(2) + &
                        phi(Ind_101)*i_pi(3))
                    g_ind(2) = g_ind(2) + 0.5d0* &
                        (phi(Ind_110)*i_pi(1) + phi(Ind_020)*i_pi(2) + &
                        phi(Ind_011)*i_pi(3))
                    g_ind(3) = g_ind(3) + 0.5d0* &
                        (phi(Ind_101)*i_pi(1) + phi(Ind_011)*i_pi(2) + &
                        phi(Ind_002)*i_pi(3))
                    ! calculate the potential at j due to polar induced at i
                    ! and derivatives of that with respect to r_j
                    ! minus signs due to deriv of r_j-r_i wrt r_i entering into potential
                    phi(Ind_000) = -(Rn_4(Ind_100)*i_pi(1) + Rn_4(Ind_010)*i_pi(2) + &
                        Rn_4(Ind_001)*i_pi(3))
                    phi(Ind_100) = -(Rn_4(Ind_200)*i_pi(1) + Rn_4(Ind_110)*i_pi(2) + &
                        Rn_4(Ind_101)*i_pi(3))
                    phi(Ind_010) = -(Rn_4(Ind_110)*i_pi(1) + Rn_4(Ind_020)*i_pi(2) + &
                        Rn_4(Ind_011)*i_pi(3))
                    phi(Ind_001) = -(Rn_4(Ind_101)*i_pi(1) + Rn_4(Ind_011)*i_pi(2) + &
                        Rn_4(Ind_002)*i_pi(3))
                    phi(Ind_200) = -(Rn_4(Ind_300)*i_pi(1) + Rn_4(Ind_210)*i_pi(2) + &
                        Rn_4(Ind_201)*i_pi(3))
                    phi(Ind_020) = -(Rn_4(Ind_120)*i_pi(1) + Rn_4(Ind_030)*i_pi(2) + &
                        Rn_4(Ind_021)*i_pi(3))
                    phi(Ind_002) = -(Rn_4(Ind_102)*i_pi(1) + Rn_4(Ind_012)*i_pi(2) + &
                        Rn_4(Ind_003)*i_pi(3))
                    phi(Ind_110) = -(Rn_4(Ind_210)*i_pi(1) + Rn_4(Ind_120)*i_pi(2) + &
                        Rn_4(Ind_111)*i_pi(3))
                    phi(Ind_101) = -(Rn_4(Ind_201)*i_pi(1) + Rn_4(Ind_111)*i_pi(2) + &
                        Rn_4(Ind_102)*i_pi(3))
                    phi(Ind_011) = -(Rn_4(Ind_111)*i_pi(1) + Rn_4(Ind_021)*i_pi(2) + &
                        Rn_4(Ind_012)*i_pi(3))
                    ! torque field at j due to direct induced mpoles of i
                    ! note factor of 1/2
                    do jj = 1, 10
                        torque_field(jj, j) = torque_field(jj, j) + 0.5d0*phi(jj)
                    end do
                end if !is_polarizable(i)
                if (is_polarizable(j)) then
                    ! induced polar dipoles at j interact with permanent at i
                    ! phi array (electrostatic potential at i due to j induced moments
                    ! and derivs of that esp wrt r_i )
                    ! minus signs arise due to derivs of r_j - r_i wrt r_i
                    phi(Ind_000) = Rn_4(Ind_100)*i_pj(1) + Rn_4(Ind_010)*i_pj(2) + &
                        Rn_4(Ind_001)*i_pj(3)
                    phi(Ind_100) = -(Rn_4(Ind_200)*i_pj(1) + Rn_4(Ind_110)*i_pj(2) + &
                        Rn_4(Ind_101)*i_pj(3))
                    phi(Ind_010) = -(Rn_4(Ind_110)*i_pj(1) + Rn_4(Ind_020)*i_pj(2) + &
                        Rn_4(Ind_011)*i_pj(3))
                    phi(Ind_001) = -(Rn_4(Ind_101)*i_pj(1) + Rn_4(Ind_011)*i_pj(2) + &
                        Rn_4(Ind_002)*i_pj(3))
                    phi(Ind_200) = Rn_4(Ind_300)*i_pj(1) + Rn_4(Ind_210)*i_pj(2) + &
                        Rn_4(Ind_201)*i_pj(3)
                    phi(Ind_020) = Rn_4(Ind_120)*i_pj(1) + Rn_4(Ind_030)*i_pj(2) + &
                        Rn_4(Ind_021)*i_pj(3)
                    phi(Ind_002) = Rn_4(Ind_102)*i_pj(1) + Rn_4(Ind_012)*i_pj(2) + &
                        Rn_4(Ind_003)*i_pj(3)
                    phi(Ind_110) = Rn_4(Ind_210)*i_pj(1) + Rn_4(Ind_120)*i_pj(2) + &
                        Rn_4(Ind_111)*i_pj(3)
                    phi(Ind_101) = Rn_4(Ind_201)*i_pj(1) + Rn_4(Ind_111)*i_pj(2) + &
                        Rn_4(Ind_102)*i_pj(3)
                    phi(Ind_011) = Rn_4(Ind_111)*i_pj(1) + Rn_4(Ind_021)*i_pj(2) + &
                        Rn_4(Ind_012)*i_pj(3)
                    phi(Ind_300) = -(Rn_4(Ind_400)*i_pj(1) + Rn_4(Ind_310)*i_pj(2) + &
                        Rn_4(Ind_301)*i_pj(3))
                    phi(Ind_030) = -(Rn_4(Ind_130)*i_pj(1) + Rn_4(Ind_040)*i_pj(2) + &
                        Rn_4(Ind_031)*i_pj(3))
                    phi(Ind_003) = -(Rn_4(Ind_103)*i_pj(1) + Rn_4(Ind_013)*i_pj(2) + &
                        Rn_4(Ind_004)*i_pj(3))
                    phi(Ind_210) = -(Rn_4(Ind_310)*i_pj(1) + Rn_4(Ind_220)*i_pj(2) + &
                        Rn_4(Ind_211)*i_pj(3))
                    phi(Ind_201) = -(Rn_4(Ind_301)*i_pj(1) + Rn_4(Ind_211)*i_pj(2) + &
                        Rn_4(Ind_202)*i_pj(3))
                    phi(Ind_120) = -(Rn_4(Ind_220)*i_pj(1) + Rn_4(Ind_130)*i_pj(2) + &
                        Rn_4(Ind_121)*i_pj(3))
                    phi(Ind_021) = -(Rn_4(Ind_121)*i_pj(1) + Rn_4(Ind_031)*i_pj(2) + &
                        Rn_4(Ind_022)*i_pj(3))
                    phi(Ind_102) = -(Rn_4(Ind_202)*i_pj(1) + Rn_4(Ind_112)*i_pj(2) + &
                        Rn_4(Ind_103)*i_pj(3))
                    phi(Ind_012) = -(Rn_4(Ind_112)*i_pj(1) + Rn_4(Ind_022)*i_pj(2) + &
                        Rn_4(Ind_013)*i_pj(3))
                    phi(Ind_111) = -(Rn_4(Ind_211)*i_pj(1) + Rn_4(Ind_121)*i_pj(2) + &
                        Rn_4(Ind_112)*i_pj(3))
                    g_ind(1) = g_ind(1) + 0.5d0* &
                        (phi(Ind_100)*gmi(Ind_000) + phi(Ind_200)*gmi(Ind_100) + &
                        phi(Ind_110)*gmi(Ind_010) + phi(Ind_101)*gmi(Ind_001) + &
                        phi(Ind_300)*gmi(Ind_200) + phi(Ind_120)*gmi(Ind_020) + &
                        phi(Ind_102)*gmi(Ind_002) + phi(Ind_210)*gmi(Ind_110) + &
                        phi(Ind_201)*gmi(Ind_101) + phi(Ind_111)*gmi(Ind_011))
                    g_ind(2) = g_ind(2) + 0.5d0* &
                        (phi(Ind_010)*gmi(Ind_000) + phi(Ind_110)*gmi(Ind_100) + &
                        phi(Ind_020)*gmi(Ind_010) + phi(Ind_011)*gmi(Ind_001) + &
                        phi(Ind_210)*gmi(Ind_200) + phi(Ind_030)*gmi(Ind_020) + &
                        phi(Ind_012)*gmi(Ind_002) + phi(Ind_120)*gmi(Ind_110) + &
                        phi(Ind_111)*gmi(Ind_101) + phi(Ind_021)*gmi(Ind_011))
                    g_ind(3) = g_ind(3) + 0.5d0* &
                        (phi(Ind_001)*gmi(Ind_000) + phi(Ind_101)*gmi(Ind_100) + &
                        phi(Ind_011)*gmi(Ind_010) + phi(Ind_002)*gmi(Ind_001) + &
                        phi(Ind_201)*gmi(Ind_200) + phi(Ind_021)*gmi(Ind_020) + &
                        phi(Ind_003)*gmi(Ind_002) + phi(Ind_111)*gmi(Ind_110) + &
                        phi(Ind_102)*gmi(Ind_101) + phi(Ind_012)*gmi(Ind_011))
                    ! torque field at i due to polar dipoles of j
                    ! note factor of 1/2
                    do jj = 1, 10
                        torque_field(jj, i) = torque_field(jj, i) + 0.5d0*phi(jj)
                    end do
                end if ! is_polarizable(j) )then
                ! finally the induced-induced interactions
                if (is_polarizable(i) .and. is_polarizable(j)) then
                    ! recall A gives coulomb, BD damped and B the erfc*1/r contributions
                    do jj = 1, 3
                        C(jj) = (B(jj) - A(jj)) + mutual_weight(k)*(A(jj) - BD(jj))
                    end do
                    ! negate the odd order to get sign right
                    C(1) = -C(1)
                    C(3) = -C(3)
                    ! get the interaction tensor
                    n = 3
                    Rn(Ind_000) = C(n)
                    Rn_1(Ind_000) = C(n - 1)
                    Rn_1(Ind_100) = delx*Rn(Ind_000)
                    Rn_1(Ind_010) = dely*Rn(Ind_000)
                    Rn_1(Ind_001) = delz*Rn(Ind_000)
                    Rn_2(Ind_000) = C(n - 2)
                    Rn_2(Ind_100) = delx*Rn_1(Ind_000)
                    Rn_2(Ind_010) = dely*Rn_1(Ind_000)
                    Rn_2(Ind_001) = delz*Rn_1(Ind_000)
                    Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
                    Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
                    Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
                    Rn_2(Ind_110) = delx*Rn_1(Ind_010)
                    Rn_2(Ind_101) = delx*Rn_1(Ind_001)
                    Rn_2(Ind_011) = dely*Rn_1(Ind_001)
                    Rn_3(Ind_100) = delx*Rn_2(Ind_000)
                    Rn_3(Ind_010) = dely*Rn_2(Ind_000)
                    Rn_3(Ind_001) = delz*Rn_2(Ind_000)
                    Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
                    Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
                    Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
                    Rn_3(Ind_110) = delx*Rn_2(Ind_010)
                    Rn_3(Ind_101) = delx*Rn_2(Ind_001)
                    Rn_3(Ind_011) = dely*Rn_2(Ind_001)
                    Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
                    Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
                    Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
                    Rn_3(Ind_210) = dely*Rn_2(Ind_200)
                    Rn_3(Ind_201) = delz*Rn_2(Ind_200)
                    Rn_3(Ind_120) = delx*Rn_2(Ind_020)
                    Rn_3(Ind_021) = delz*Rn_2(Ind_020)
                    Rn_3(Ind_102) = delx*Rn_2(Ind_002)
                    Rn_3(Ind_012) = dely*Rn_2(Ind_002)
                    Rn_3(Ind_111) = delx*Rn_2(Ind_011)
                    ! phi array (electrostatic pot at i due to j induced direct
                    ! and derivs of that esp wrt r_i )
                    phi(Ind_200) = Rn_3(Ind_300)*i_dj(1) + Rn_3(Ind_210)*i_dj(2) + &
                        Rn_3(Ind_201)*i_dj(3)
                    phi(Ind_020) = Rn_3(Ind_120)*i_dj(1) + Rn_3(Ind_030)*i_dj(2) + &
                        Rn_3(Ind_021)*i_dj(3)
                    phi(Ind_002) = Rn_3(Ind_102)*i_dj(1) + Rn_3(Ind_012)*i_dj(2) + &
                        Rn_3(Ind_003)*i_dj(3)
                    phi(Ind_110) = Rn_3(Ind_210)*i_dj(1) + Rn_3(Ind_120)*i_dj(2) + &
                        Rn_3(Ind_111)*i_dj(3)
                    phi(Ind_101) = Rn_3(Ind_201)*i_dj(1) + Rn_3(Ind_111)*i_dj(2) + &
                        Rn_3(Ind_102)*i_dj(3)
                    phi(Ind_011) = Rn_3(Ind_111)*i_dj(1) + Rn_3(Ind_021)*i_dj(2) + &
                        Rn_3(Ind_012)*i_dj(3)
                    g_ind(1) = g_ind(1) + 0.5d0* &
                        (phi(Ind_200)*i_pi(1) + phi(Ind_110)*i_pi(2) + &
                        phi(Ind_101)*i_pi(3))
                    g_ind(2) = g_ind(2) + 0.5d0* &
                        (phi(Ind_110)*i_pi(1) + phi(Ind_020)*i_pi(2) + &
                        phi(Ind_011)*i_pi(3))
                    g_ind(3) = g_ind(3) + 0.5d0* &
                        (phi(Ind_101)*i_pi(1) + phi(Ind_011)*i_pi(2) + &
                        phi(Ind_002)*i_pi(3))
                    ! phi array (electrostatic pot at i due to j induced polar
                    ! and derivs of that esp wrt r_i )
                    phi(Ind_200) = Rn_3(Ind_300)*i_pj(1) + Rn_3(Ind_210)*i_pj(2) + &
                        Rn_3(Ind_201)*i_pj(3)
                    phi(Ind_020) = Rn_3(Ind_120)*i_pj(1) + Rn_3(Ind_030)*i_pj(2) + &
                        Rn_3(Ind_021)*i_pj(3)
                    phi(Ind_002) = Rn_3(Ind_102)*i_pj(1) + Rn_3(Ind_012)*i_pj(2) + &
                        Rn_3(Ind_003)*i_pj(3)
                    phi(Ind_110) = Rn_3(Ind_210)*i_pj(1) + Rn_3(Ind_120)*i_pj(2) + &
                        Rn_3(Ind_111)*i_pj(3)
                    phi(Ind_101) = Rn_3(Ind_201)*i_pj(1) + Rn_3(Ind_111)*i_pj(2) + &
                        Rn_3(Ind_102)*i_pj(3)
                    phi(Ind_011) = Rn_3(Ind_111)*i_pj(1) + Rn_3(Ind_021)*i_pj(2) + &
                        Rn_3(Ind_012)*i_pj(3)
                    g_ind(1) = g_ind(1) + 0.5d0* &
                        (phi(Ind_200)*i_di(1) + phi(Ind_110)*i_di(2) + &
                        phi(Ind_101)*i_di(3))
                    g_ind(2) = g_ind(2) + 0.5d0* &
                        (phi(Ind_110)*i_di(1) + phi(Ind_020)*i_di(2) + &
                        phi(Ind_011)*i_di(3))
                    g_ind(3) = g_ind(3) + 0.5d0* &
                        (phi(Ind_101)*i_di(1) + phi(Ind_011)*i_di(2) + &
                        phi(Ind_002)*i_di(3))
                end if !is_polarizable(i) .and. is_polarizable(j)
            end if !( is_polarizable(i) .or. is_polarizable(j) )then
            ! frc is negative gradient
            ene_perm = ene_perm + coulomb_const_kcal_per_mole*e_pp
            ene_ind = ene_ind + coulomb_const_kcal_per_mole*e_ind
            frc(1, j) = frc(1, j) + coulomb_const_kcal_per_mole*(g_pp(1) + g_ind(1))
            frc(2, j) = frc(2, j) + coulomb_const_kcal_per_mole*(g_pp(2) + g_ind(2))
            frc(3, j) = frc(3, j) + coulomb_const_kcal_per_mole*(g_pp(3) + g_ind(3))
            frc(1, i) = frc(1, i) - coulomb_const_kcal_per_mole*(g_pp(1) + g_ind(1))
            frc(2, i) = frc(2, i) - coulomb_const_kcal_per_mole*(g_pp(2) + g_ind(2))
            frc(3, i) = frc(3, i) - coulomb_const_kcal_per_mole*(g_pp(3) + g_ind(3))
            vxx = vxx - coulomb_const_kcal_per_mole*delx*(g_pp(1) + g_ind(1))
            vxy = vxy - coulomb_const_kcal_per_mole*delx*(g_pp(2) + g_ind(2))
            vxz = vxz - coulomb_const_kcal_per_mole*delx*(g_pp(3) + g_ind(3))
            vyx = vyx - coulomb_const_kcal_per_mole*dely*(g_pp(1) + g_ind(1))
            vyy = vyy - coulomb_const_kcal_per_mole*dely*(g_pp(2) + g_ind(2))
            vyz = vyz - coulomb_const_kcal_per_mole*dely*(g_pp(3) + g_ind(3))
            vzx = vzx - coulomb_const_kcal_per_mole*delz*(g_pp(1) + g_ind(1))
            vzy = vzy - coulomb_const_kcal_per_mole*delz*(g_pp(2) + g_ind(2))
            vzz = vzz - coulomb_const_kcal_per_mole*delz*(g_pp(3) + g_ind(3))
        end do !n_adj = 1,num_adjust_list
        virial(1, 1) = virial(1, 1) + vxx
        virial(1, 2) = virial(1, 2) + half*(vxy + vyx)
        virial(1, 3) = virial(1, 3) + half*(vxz + vzx)
        virial(2, 1) = virial(2, 1) + half*(vxy + vyx)
        virial(2, 2) = virial(2, 2) + vyy
        virial(2, 3) = virial(2, 3) + half*(vyz + vzy)
        virial(3, 1) = virial(3, 1) + half*(vxz + vzx)
        virial(3, 2) = virial(3, 2) + half*(vyz + vzy)
        virial(3, 3) = virial(3, 3) + vzz
    end subroutine AM_ADJUST_calc_ene_frc
!-------------------------------------------------------
end module amoeba_adjust
