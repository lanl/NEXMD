#include "dprec.fh"
#include "assert.fh"

!------------------------------------------------------------
module amoeba_direct
    implicit none
    private

#  include "amoeba_mpole_index.h"
    integer, save :: do_flag
    integer, save ::  num_pairs_in_ee_cut = 0, size_dipole_dipole_list = -1
    integer, save :: num_tensor
    _REAL_, save, allocatable :: dipole_dipole_tensor(:, :)
    integer, save, allocatable :: dipole_dipole_list(:, :)
    _REAL_, parameter :: safety = 1.25d0, checklist = 0.9d0

    public AM_DIRECT_permfield, AM_DIRECT_dip_dip_field, AM_DIRECT_ene_frc, &
        AM_DIRECT_set_user_bit
#ifdef MPI
    public AM_DIRECT_bcast
#endif

contains

#ifdef MPI
    subroutine AM_DIRECT_bcast
        implicit none
        integer ierr

        include 'mpif.h'
# include "parallel.h"
        call mpi_bcast(do_flag, 1, MPI_INTEGER, 0, commsander, ierr)
    end subroutine AM_DIRECT_bcast
#endif

!-----------------------------------------------------------
    subroutine AM_DIRECT_set_user_bit(do_this)
        integer, intent(in) :: do_this
#include "do_flag.h"

        ! set the valid bit---this part always since no parmread needed
        do_flag = ibset(do_flag, VALID_BIT)

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
            write (6, *) 'AM_DIRECT_set_user_bit: bad value of user do_this'
            call mexit(6, 1)
        end if
    end subroutine AM_DIRECT_set_user_bit
!-----------------------------------------------------------
    subroutine AM_DIRECT_permfield(ipairs, x, cart_dipole_field)

        use nblist, only : imagcrds, bckptr, nlogrid, nhigrid, numvdw, numhbnd, &
            myindexlo, myindexhi, numimg
        use amoeba_multipoles, only : global_multipole
        use amoeba_induced, only : sq_polinv, is_polarizable
        use amoeba_mdin, only : ee_dsum_cut, ee_damped_cut, thole_expon_coeff

        integer, intent(in) :: ipairs(*)
        _REAL_, intent(in) :: x(*)
        _REAL_, intent(inout) :: cart_dipole_field(3, *)

#  include "extra.h"
#  include "def_time.h"
#  include "ew_erfc_spline.h"
#  include "ew_pme_recip.h"

        integer :: i, k, numpack, index, ncell_lo, ncell_hi, ntot, ier
        _REAL_ xk, yk, zk

        call timer_start(TIME_SHORT_ENE)
        ! check if dipole_dipole_list big enough
        if (num_pairs_in_ee_cut > checklist*size_dipole_dipole_list) then
            if (allocated(dipole_dipole_tensor)) deallocate (dipole_dipole_tensor)
            if (allocated(dipole_dipole_list)) deallocate (dipole_dipole_list)
            call AM_DIRECT_count_num_ee_pairs( &
                ipairs, ee_dsum_cut, num_pairs_in_ee_cut)
            size_dipole_dipole_list = safety*num_pairs_in_ee_cut
            if (master) then
                write (6, '(a,i10,i10)') &
                    'num_pairs_in_ee_cut,size_dipole_dipole_list = ', &
                    num_pairs_in_ee_cut, size_dipole_dipole_list
            end if
            allocate (dipole_dipole_tensor(6, size_dipole_dipole_list), stat=ier)
            REQUIRE(ier == 0)
            allocate (dipole_dipole_list(2, size_dipole_dipole_list), stat=ier)
            REQUIRE(ier == 0)
        end if
        numpack = 1
        num_tensor = 0
        do index = myindexlo, myindexhi
            if (numimg(index) > 0) then
                ncell_lo = nlogrid(index)
                ncell_hi = nhigrid(index)
                do k = ncell_lo, ncell_hi
                    i = bckptr(k)
                    xk = imagcrds(1, k)
                    yk = imagcrds(2, k)
                    zk = imagcrds(3, k)
                    ntot = numvdw(i) + numhbnd(i)
                    if (ntot > 0) then
                        call AM_DIRECT_permfield_i(i, ipairs(numpack), ntot, &
                            xk, yk, zk, ew_coeff, eedtbdns, x(leed_cub), x(leed_lin), &
                            ee_type, eedmeth, dxdr, ee_dsum_cut, &
                            ee_damped_cut, thole_expon_coeff, &
                            sq_polinv, is_polarizable, &
                            dipole_dipole_tensor, dipole_dipole_list, &
                            size_dipole_dipole_list, num_tensor, &
                            global_multipole, cart_dipole_field)
                        numpack = numpack + ntot
                    end if  ! ( ntot > 0 )
                end do  !  k = ncell_lo,ncell_hi
            end if  ! ( numimg(k) > 0 )
        end do  !  index = myindexlo,myindexhi
        call timer_stop(TIME_SHORT_ENE)
        return
    end subroutine AM_DIRECT_permfield
!-------------------------------------------------------
    subroutine AM_DIRECT_ene_frc(ipairs, crd, x, ind_dip_d, ind_dip_p, &
        ene_perm, ene_ind, ene_vdw, frc, virial)

        use nblist, only : imagcrds, bckptr, nlogrid, nhigrid, numvdw, numhbnd, &
            myindexlo, myindexhi, numimg
        use amoeba_multipoles, only : global_multipole
        use amoeba_induced, only : sq_polinv, is_polarizable
        use amoeba_mdin, only : ee_dsum_cut, ee_damped_cut, thole_expon_coeff
        use amoeba_vdw, only : AM_VDW_DIRECT_ene_frc_i

        integer, intent(in) :: ipairs(*)
        _REAL_, intent(in) :: crd(3, *), x(*)
        _REAL_, intent(in) :: ind_dip_d(3, *), ind_dip_p(3, *)
        _REAL_, intent(inout) :: ene_perm, ene_ind, ene_vdw, frc(3, *), virial(3, 3)
#  include "parallel.h"
#  include "def_time.h"
#  include "ew_erfc_spline.h"
#  include "ew_pme_recip.h"

        integer :: i, k, numpack, index, ncell_lo, ncell_hi, ntot
        _REAL_ xk, yk, zk

        numpack = 1
        ene_perm = 0.d0
        ene_ind = 0.d0
        ene_vdw = 0.d0
        call timer_start(TIME_SHORT_ENE)
        do index = myindexlo, myindexhi
            if (numimg(index) > 0) then
                ncell_lo = nlogrid(index)
                ncell_hi = nhigrid(index)
                do k = ncell_lo, ncell_hi
                    i = bckptr(k)
                    xk = imagcrds(1, k)
                    yk = imagcrds(2, k)
                    zk = imagcrds(3, k)
                    ntot = numvdw(i) + numhbnd(i)
                    if (ntot > 0) then
                        call AM_DIRECT_ene_force_i(i, ipairs(numpack), ntot, xk, yk, zk, &
                            ew_coeff, eedtbdns, &
                            x(leed_cub), x(leed_lin), &
                            ee_type, eedmeth, dxdr, ee_dsum_cut, &
                            ee_damped_cut, thole_expon_coeff, &
                            sq_polinv, is_polarizable, &
                            ind_dip_d, ind_dip_p, &
                            global_multipole, ene_perm, &
                            ene_ind, frc, virial)
                        call AM_VDW_DIRECT_ene_frc_i(i, ipairs(numpack), ntot, xk, yk, zk, &
                            crd, ene_vdw, frc, virial)
                        numpack = numpack + ntot
                    end if
                end do
            end if
        end do
        call timer_stop(TIME_SHORT_ENE)
        return
    end subroutine AM_DIRECT_ene_frc
!-------------------------------------------------------
    subroutine AM_DIRECT_count_num_ee_pairs( &
        ipairs, ee_dsum_cut, num_pairs_in_ee_cut)
        use nblist, only : imagcrds, bckptr, nlogrid, nhigrid, numvdw, numhbnd, &
            myindexlo, myindexhi, numimg
        integer, intent(in) :: ipairs(*)
        _REAL_, intent(in) :: ee_dsum_cut
        integer, intent(out) :: num_pairs_in_ee_cut

        integer :: i, k, numpack, index, ncell_lo, ncell_hi, ntot
        _REAL_ xk, yk, zk

        numpack = 1
        num_pairs_in_ee_cut = 0
        do index = myindexlo, myindexhi
            if (numimg(index) > 0) then
                ncell_lo = nlogrid(index)
                ncell_hi = nhigrid(index)
                do k = ncell_lo, ncell_hi
                    i = bckptr(k)
                    xk = imagcrds(1, k)
                    yk = imagcrds(2, k)
                    zk = imagcrds(3, k)
                    ntot = numvdw(i) + numhbnd(i)
                    if (ntot > 0) then
                        call AM_DIRECT_increment_ee_pairs( &
                            ipairs(numpack), ntot, &
                            xk, yk, zk, ee_dsum_cut, &
                            num_pairs_in_ee_cut)
                        numpack = numpack + ntot
                    end if
                end do
            end if
        end do

    end subroutine AM_DIRECT_count_num_ee_pairs
!-------------------------------------------------------
    subroutine AM_DIRECT_increment_ee_pairs( &
        ipairs, numtot, xk, yk, zk, ee_dsum_cut, &
        num_pairs_in_ee_cut)
        use nblist, only : imagcrds, tranvec
        integer, intent(in) :: ipairs(*), numtot
        _REAL_, intent(in) :: xk, yk, zk, ee_dsum_cut
        integer, intent(inout) :: num_pairs_in_ee_cut

        _REAL_ :: ee_dsum_cut2
        _REAL_ :: xktran(3, 18)
        integer :: mask27
        integer :: m, np, itran
        _REAL_ :: delx, dely, delz, delr2

        mask27 = 2**27 - 1
        ee_dsum_cut2 = ee_dsum_cut*ee_dsum_cut
        do m = 1, 18
            xktran(1, m) = tranvec(1, m) - xk
            xktran(2, m) = tranvec(2, m) - yk
            xktran(3, m) = tranvec(3, m) - zk
        end do
        do m = 1, numtot
            np = ipairs(m)
            itran = ishft(np, -27)
            np = iand(np, mask27)
            delx = imagcrds(1, np) + xktran(1, itran)
            dely = imagcrds(2, np) + xktran(2, itran)
            delz = imagcrds(3, np) + xktran(3, itran)
            delr2 = delx*delx + dely*dely + delz*delz
            if (delr2 < ee_dsum_cut2) then
                num_pairs_in_ee_cut = num_pairs_in_ee_cut + 1
            end if
        end do
    end subroutine AM_DIRECT_increment_ee_pairs
!-------------------------------------------------------
    subroutine AM_DIRECT_permfield_i(i, ipairs, numtot, xk, yk, zk, &
        ewaldcof, eedtbdns, eed_cub, eed_lin, &
        ee_type, eedmeth, dxdr, ee_dsum_cut, &
        ee_damped_cut, thole_expon_coeff, &
        sq_polinv, is_polarizable, &
        dipole_dipole_tensor, dipole_dipole_list, &
        size_dipole_dipole_list, num_tensor, &
        global_multipole, gradphi)
        use nblist, only : bckptr, imagcrds, tranvec
        use constants, only : zero, third, half, one, two, three, five

        integer, intent(in) :: i, ipairs(*), numtot
        _REAL_, intent(in) :: xk, yk, zk, ewaldcof, eedtbdns, eed_cub(4, *), eed_lin(2, *)
        integer, intent(in) :: ee_type, eedmeth
        _REAL_, intent(in) :: dxdr, ee_dsum_cut, ee_damped_cut, &
            thole_expon_coeff, sq_polinv(*)
        logical, intent(in) :: is_polarizable(*)
        _REAL_, intent(out) :: dipole_dipole_tensor(6, *)
        integer, intent(out) :: dipole_dipole_list(2, *)
        integer, intent(in) :: size_dipole_dipole_list
        integer, intent(inout) :: num_tensor
        _REAL_, intent(in) :: global_multipole(10, *)
        _REAL_, intent(inout) :: gradphi(3, *)

        ! local variables
        _REAL_ :: ee_dsum_cut2, ee_damped_cut2
        _REAL_ :: xktran(3, 18)
        integer :: mask27
        integer :: m, n, np, itran, j, ind
        _REAL_ :: delx, dely, delz, delr2, delr, delr2inv, x, dx, switch, d_switch_dx, xx
        _REAL_ :: B(0:3), BD(3), fac, fact, del, gphi_i(3), gphi_j(3)
        _REAL_  :: asq, expon, expo, clam3, clam5, clam7, delr3inv, delr5inv, delr7inv
        _REAL_  :: Rn(1), Rn_1(4), Rn_2(10), Rn_3(20)
#  include "do_flag.h"

        if (iand(do_flag, PROCEED_INDUCE) /= PROCEED_INDUCE) return

        mask27 = 2**27 - 1
        fac = two*ewaldcof*ewaldcof
        del = one/eedtbdns
        ee_dsum_cut2 = ee_dsum_cut*ee_dsum_cut
        ee_damped_cut2 = ee_damped_cut*ee_damped_cut
        asq = thole_expon_coeff*sq_polinv(i)

        do m = 1, 18
            xktran(1, m) = tranvec(1, m) - xk
            xktran(2, m) = tranvec(2, m) - yk
            xktran(3, m) = tranvec(3, m) - zk
        end do
        do m = 1, numtot
            np = ipairs(m)
            itran = ishft(np, -27)
            np = iand(np, mask27)
            j = bckptr(np)
            delx = imagcrds(1, np) + xktran(1, itran)
            dely = imagcrds(2, np) + xktran(2, itran)
            delz = imagcrds(3, np) + xktran(3, itran)
            delr2 = delx*delx + dely*dely + delz*delz
            if ((delr2 < ee_dsum_cut2) .and. &
                (is_polarizable(i) .or. is_polarizable(j))) then
                !---------------------------------------------------------
                ! McMurchie-Davidson recursion for interaction tensor:
                ! interaction at point charge level given by complementary Boys
                ! B(0) = int_0^1 exp(-pr^2^t^2)dt
                ! complementary Boys is BC(0) = 1/r - B(0)
                ! (d/dr) B(0) = (-2p)*r int_0^1 t^2*exp(-pr^2t^2)dt = (-2p)*r B(1)
                ! and so if R(0,0,0,n) = (-2p)^n B(n) then we have
                ! (1/r)*(d/dr) R(0,0,0,n) = R(0,0,0,n+1)
                ! Now let R(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
                ! Then e.g. R(t+1,u,v,n) = t*R(t-1,u,v,n+1) + x*R(t,u,v,n)
                ! proof:
                ! R(t+1,u,v,n) = (d/dx)^t+1 (d/dy)^u (d/dz)^v R(0,0,0,n)
                !              = (d/dx)^t (d/dy)^u (d/dz)^v (d/dx) R(0,0,0,n)
                !              = (d/dx)^t (d/dy)^u (d/dz)^v x*(1/r)(d/dr)R(0,0,0,n)
                !              = (d/dx)^t (d/dy)^u (d/dz)^v [x*R(0,0,0,n+1)]
                !              = (d/dx)^t [x*R(0,u,v,n+1)]
                !              = t*R(t-1,u,v,n+1) + x*R(t,u,v,n+1) (Leibniz)
                ! similar recursions hold for R(t,u+1,v,n),R(t,u,v+1,n)
                ! R(t,u+1,v,n) = u*R(t,u-1,v,n+1) + y*R(t,u,v,n+1)
                ! R(t,u,v+1,n) = v*R(t,u,v-1,n+1) + z*R(t,u,v,n+1)
                ! below array is packed---hence use of Ind_tuv
                ! Rn(Ind_tuv) denotes R(t,u,v,n)
                ! Due to its form--we recur downwards in n
                !---------------------------------------------------------
                ! top n is 3 for dipole fields
                ! get boys and R(0,0,0,n), n=0,1,2,3
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
                    switch = (one - dx)*eed_lin(1, ind) + &
                        dx*eed_lin(1, ind + 1)
                    d_switch_dx = (one - dx)*eed_lin(2, ind) + &
                        dx*eed_lin(2, ind + 1)
                else if (eedmeth == 3) then
                    !           ---direct function call:
                    call get_ee_func(x, switch, d_switch_dx, ee_type)
                else if (eedmeth == 4) then
                    !            ---use un-modified Coulomb interaction, no switch
                    switch = one
                    d_switch_dx = zero
                else
                    write (6, *) 'bad eedmeth in ew_short_dip: ', eedmeth
                    call mexit(6, 1)
                end if  ! ( eedmeth == 1 )

                ! TD Got the idea for B_l from Walter Smith's CCP5 article 1982
                ! Ewald for point multipoles
                ! B_l satisfies grad_i B_l(|r_j - r_i|) = (r_j - r_i)B_{l+1}(|r_j-r_i|)
                ! grad_j B_l(|r_j - r_i|) = -grad_i B_l(|r_j - r_i|)

                B(0) = switch*delr*delr2inv
                fact = d_switch_dx*dxdr
                B(1) = (B(0) - fact)*delr2inv
                fact = fac*fact
                B(2) = (three*B(1) - fact)*delr2inv
                fact = fac*fact
                B(3) = (five*B(2) - fact)*delr2inv
                if (delr2 < ee_damped_cut2) then
                    !-------------------------------------------------------
                    ! McMurchie-Davidson holds for damped tensor as well---in fact,
                    ! BD below satisfies BD(n+1) = (1/r)(d/dr)BD(n)
                    ! RD(0,0,0,n) = BD(n), n=1,2,3
                    ! RD(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
                    !-------------------------------------------------------
                    delr3inv = delr2inv/delr
                    delr5inv = delr3inv*delr2inv
                    delr7inv = delr5inv*delr2inv
                    expon = asq*delr2*delr*sq_polinv(j)
                    expo = exp(-expon)
                    ! clam3 = 1.d0-lam3, clam5 = 1.d0-lam5 etc. where
                    ! lam is from ponder's paper
                    clam3 = expo
                    clam5 = (1.d0 + expon)*expo
                    clam7 = (1.d0 + expon + 0.6d0*expon**2)*expo
                    BD(1) = clam3*delr3inv
                    BD(2) = 3.d0*clam5*delr5inv
                    BD(3) = 15.d0*clam7*delr7inv
                    ! correct the boys factors by damped factors
                    ! ewald dsum tensor and damped tensor both satisfy
                    ! McMurchie-Davidson recursion-thus their sum does as well
                    ! use recur for sum--starting with summed factors
                    B(1) = B(1) - BD(1)
                    B(2) = B(2) - BD(2)
                    B(3) = B(3) - BD(3)
                end if !( delr2 < damped_cut2 )then
                ! negate the odd order boys factors
                B(1) = -B(1)
                B(3) = -B(3)
                n = 3
                Rn(Ind_000) = B(n)
                Rn_1(Ind_000) = B(n - 1)
                Rn_1(Ind_100) = delx*Rn(Ind_000)
                Rn_1(Ind_010) = dely*Rn(Ind_000)
                Rn_1(Ind_001) = delz*Rn(Ind_000)
                Rn_2(Ind_000) = B(n - 2)
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
                ! store the dipole-dipole component
                num_tensor = num_tensor + 1
                if (num_tensor > size_dipole_dipole_list) then
                    write (6, *) 'Too many dipole_dipole interactions for allocated'
                    call mexit(6, 1)
                end if
                dipole_dipole_list(1, num_tensor) = i
                dipole_dipole_list(2, num_tensor) = j
                dipole_dipole_tensor(1, num_tensor) = Rn_3(Ind_200)
                dipole_dipole_tensor(2, num_tensor) = Rn_3(Ind_110)
                dipole_dipole_tensor(3, num_tensor) = Rn_3(Ind_101)
                dipole_dipole_tensor(4, num_tensor) = Rn_3(Ind_020)
                dipole_dipole_tensor(5, num_tensor) = Rn_3(Ind_011)
                dipole_dipole_tensor(6, num_tensor) = Rn_3(Ind_002)
                ! next do the field components
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
                    ! minus sign due to derivs wrt crds of i
                    gradphi(1, i) = gradphi(1, i) - gphi_i(1)
                    gradphi(2, i) = gradphi(2, i) - gphi_i(2)
                    gradphi(3, i) = gradphi(3, i) - gphi_i(3)
                end if
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
                    gradphi(1, j) = gradphi(1, j) + gphi_j(1)
                    gradphi(2, j) = gradphi(2, j) + gphi_j(2)
                    gradphi(3, j) = gradphi(3, j) + gphi_j(3)
                end if !( is_polarizable(j) )then
            end if !( delr2 < ee_dsum_cut2 )then
        end do !m = 1,numlist
    end subroutine AM_DIRECT_permfield_i
!-------------------------------------------------------
    subroutine AM_DIRECT_dip_dip_field( &
        ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)
        _REAL_, intent(in) :: ind_dip_d(3, *), ind_dip_p(3, *)
        _REAL_, intent(inout) :: dip_field_d(3, *), dip_field_p(3, *)
#  include "do_flag.h"

        if (iand(do_flag, PROCEED_INDUCE) /= PROCEED_INDUCE) return

        call timer_start(TIME_SHORT_ENE)
        call AM_DIRECT_calc_dipdip_field( &
            num_tensor, &
            dipole_dipole_list, dipole_dipole_tensor, &
            ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)
        call timer_stop(TIME_SHORT_ENE)

    end subroutine AM_DIRECT_dip_dip_field
!-------------------------------------------------------
    subroutine AM_DIRECT_calc_dipdip_field( &
        num_tensor, &
        dipole_dipole_list, dipole_dipole_tensor, &
        ind_dip_d, ind_dip_p, dip_field_d, dip_field_p)
        integer, intent(in) :: num_tensor, dipole_dipole_list(2, *)
        _REAL_, intent(in) :: dipole_dipole_tensor(6, *)
        _REAL_, intent(in) :: ind_dip_d(3, *), ind_dip_p(3, *)
        _REAL_, intent(inout) :: dip_field_d(3, *), dip_field_p(3, *)

        integer :: i, j, n

        do n = 1, num_tensor
            i = dipole_dipole_list(1, n)
            j = dipole_dipole_list(2, n)
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
            ! other set of dipoles, fields
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
    end subroutine AM_DIRECT_calc_dipdip_field
!-------------------------------------------------------
    subroutine AM_DIRECT_ene_force_i(i, ipairs, numtot, xk, yk, zk, &
        ewaldcof, eedtbdns, eed_cub, eed_lin, &
        ee_type, eedmeth, dxdr, ee_dsum_cut, &
        ee_damped_cut, thole_expon_coeff, &
        sq_polinv, is_polarizable, &
        ind_dip_d, ind_dip_p, &
        global_multipole, ene_perm, &
        ene_ind, frc, virial)
        use nblist, only : bckptr, imagcrds, tranvec
        use amoeba_multipoles, only : coulomb_const_kcal_per_mole, torque_field
        use constants, only : zero, third, half, one, two, three, five, seven, nine

        integer, intent(in) :: i, ipairs(*), numtot
        _REAL_, intent(in) :: xk, yk, zk, ewaldcof, eedtbdns, eed_cub(4, *), eed_lin(2, *)
        integer, intent(in) :: ee_type, eedmeth
        _REAL_, intent(in) :: dxdr, ee_dsum_cut, ee_damped_cut, &
            thole_expon_coeff, sq_polinv(*)
        logical, intent(in) :: is_polarizable(*)
        _REAL_, intent(in) :: ind_dip_d(3, *), ind_dip_p(3, *)
        _REAL_, intent(in) :: global_multipole(10, *)
        _REAL_, intent(inout) :: ene_perm, ene_ind, frc(3, *), virial(3, 3)

        ! local variables
        _REAL_ :: ee_dsum_cut2, ee_damped_cut2
        _REAL_ :: xktran(3, 18)
        integer :: mask27
        integer :: m, n, np, itran, j, ind, jj
        _REAL_ :: delx, dely, delz, delr2, delr, delr2inv, x, dx, switch, d_switch_dx, xx
        _REAL_ :: vxx, vxy, vxz, vyx, vyy, vyz, vzx, vzy, vzz
        _REAL_ :: B(0:5), BD(5), fac, fact, del, phi(20), gmj(10), gmi(10), tmi(10)
        _REAL_ :: e_pp, e_ind, g_pp(3), g_ind(3), &
            i_di(3), i_pi(3), i_dj(3), i_pj(3), i_mi(3), i_mj(3)
        _REAL_  :: asq, expon, expo, clam3, clam5, clam7, clam9, &
            delr3inv, delr5inv, delr7inv, delr9inv
        _REAL_  :: Rn(1), Rn_1(4), Rn_2(10), Rn_3(20), Rn_4(35), Rn_5(56)
        _REAL_, parameter :: const1 = 0.6d0, const2 = 18.d0/35.d0, const3 = 9.d0/35.d0
#  include "do_flag.h"

        if (iand(do_flag, PROCEED_POSTINDUCE) /= PROCEED_POSTINDUCE) return

        mask27 = 2**27 - 1
        fac = two*ewaldcof*ewaldcof
        del = one/eedtbdns
        ee_dsum_cut2 = ee_dsum_cut*ee_dsum_cut
        ee_damped_cut2 = ee_damped_cut*ee_damped_cut
        asq = thole_expon_coeff*sq_polinv(i)

        vxx = zero
        vxy = zero
        vxz = zero
        vyx = zero
        vyy = zero
        vyz = zero
        vzx = zero
        vzy = zero
        vzz = zero
        do m = 1, 18
            xktran(1, m) = tranvec(1, m) - xk
            xktran(2, m) = tranvec(2, m) - yk
            xktran(3, m) = tranvec(3, m) - zk
        end do
        do m = 1, numtot
            np = ipairs(m)
            itran = ishft(np, -27)
            np = iand(np, mask27)
            j = bckptr(np)
            delx = imagcrds(1, np) + xktran(1, itran)
            dely = imagcrds(2, np) + xktran(2, itran)
            delz = imagcrds(3, np) + xktran(3, itran)
            delr2 = delx*delx + dely*dely + delz*delz
            if (delr2 < ee_dsum_cut2) then
                !---------------------------------------------------------
                ! McMurchie-Davidson recursion for interaction tensor:
                ! interaction at point charge level given by complementary Boys
                ! B(0) = int_0^1 exp(-pr^2^t^2)dt
                ! complementary Boys is BC(0) = 1/r - B(0)
                ! (d/dr) B(0) = (-2p)*r int_0^1 t^2*exp(-pr^2t^2)dt = (-2p)*r B(1)
                ! and so if R(0,0,0,n) = (-2p)^n B(n) then we have
                ! (1/r)*(d/dr) R(0,0,0,n) = R(0,0,0,n+1)
                ! Now let R(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
                ! Then e.g. R(t+1,u,v,n) = t*R(t-1,u,v,n+1) + x*R(t,u,v,n)
                ! proof:
                ! R(t+1,u,v,n) = (d/dx)^t+1 (d/dy)^u (d/dz)^v R(0,0,0,n)
                !              = (d/dx)^t (d/dy)^u (d/dz)^v (d/dx) R(0,0,0,n)
                !              = (d/dx)^t (d/dy)^u (d/dz)^v x*(1/r)(d/dr)R(0,0,0,n)
                !              = (d/dx)^t (d/dy)^u (d/dz)^v [x*R(0,0,0,n+1)]
                !              = (d/dx)^t [x*R(0,u,v,n+1)]
                !              = t*R(t-1,u,v,n+1) + x*R(t,u,v,n+1) (Leibniz)
                ! similar recursions hold for R(t,u+1,v,n),R(t,u,v+1,n)
                ! R(t,u+1,v,n) = u*R(t,u-1,v,n+1) + y*R(t,u,v,n+1)
                ! R(t,u,v+1,n) = v*R(t,u,v-1,n+1) + z*R(t,u,v,n+1)
                ! below array is packed---hence use of Ind_tuv
                ! Rn(Ind_tuv) denotes R(t,u,v,n)
                ! Due to its form--we recur downwards in n
                !---------------------------------------------------------
                ! top n is 5 for energy and forces
                ! get boys and R(0,0,0,n), n=0,...,5
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
                    switch = (one - dx)*eed_lin(1, ind) + &
                        dx*eed_lin(1, ind + 1)
                    d_switch_dx = (one - dx)*eed_lin(2, ind) + &
                        dx*eed_lin(2, ind + 1)
                else if (eedmeth == 3) then
                    !           ---direct function call:
                    call get_ee_func(x, switch, d_switch_dx, ee_type)
                else if (eedmeth == 4) then
                    !            ---use un-modified Coulomb interaction, no switch
                    switch = one
                    d_switch_dx = zero
                else
                    write (6, *) 'bad eedmeth in ew_short_dip: ', eedmeth
                    call mexit(6, 1)
                end if  ! ( eedmeth == 1 )

                ! TD Got the idea for B_l from Walter Smith's CCP5 article 1982
                ! Ewald for point multipoles
                ! B_l satisfies grad_i B_l(|r_j - r_i|) = (r_j - r_i)B_{l+1}(|r_j-r_i|)

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
                ! negate the odd order boys factors
                B(1) = -B(1)
                B(3) = -B(3)
                B(5) = -B(5)
                n = 5
                Rn(Ind_000) = B(n)
                Rn_1(Ind_000) = B(n - 1)
                Rn_1(Ind_100) = delx*Rn(Ind_000)
                Rn_1(Ind_010) = dely*Rn(Ind_000)
                Rn_1(Ind_001) = delz*Rn(Ind_000)
                Rn_2(Ind_000) = B(n - 2)
                Rn_2(Ind_100) = delx*Rn_1(Ind_000)
                Rn_2(Ind_010) = dely*Rn_1(Ind_000)
                Rn_2(Ind_001) = delz*Rn_1(Ind_000)
                Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
                Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
                Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
                Rn_2(Ind_110) = delx*Rn_1(Ind_010)
                Rn_2(Ind_101) = delx*Rn_1(Ind_001)
                Rn_2(Ind_011) = dely*Rn_1(Ind_001)
                Rn_3(Ind_000) = B(n - 3)
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
                Rn_4(Ind_000) = B(n - 4)
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
                Rn_5(Ind_000) = B(n - 5)
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
                ! phi array (electrostatic potential at i due to j permanent mpoles
                ! and derivs of that esp wrt r_i )
                ! minus signs arise due to derivs of r_j - r_i wrt r_i
                do jj = 1, 10
                    gmi(jj) = global_multipole(jj, i)
                    tmi(jj) = global_multipole(jj, i) !used for torque contrib
                    gmj(jj) = global_multipole(jj, j)
                end do
                do jj = 1, 3
                    g_ind(jj) = 0.d0
                    i_dj(jj) = ind_dip_d(jj, j)
                    i_di(jj) = ind_dip_d(jj, i)
                    i_pj(jj) = ind_dip_p(jj, j)
                    i_pi(jj) = ind_dip_p(jj, i)
                    i_mj(jj) = ind_dip_d(jj, j) + ind_dip_p(jj, j)
                    i_mi(jj) = ind_dip_d(jj, i) + ind_dip_p(jj, i)
                    tmi(jj + 1) = tmi(jj + 1) + 0.5d0*i_mi(jj) !used for torque contrib
                end do
                ! initialize induction contributions
                e_ind = 0.d0
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
                if (is_polarizable(i)) then
                    e_ind = e_ind + 0.5d0* &
                        (phi(Ind_100)*i_di(1) + phi(Ind_010)*i_di(2) + &
                        phi(Ind_001)*i_di(3))
                    g_ind(1) = g_ind(1) + 0.5d0* &
                        (phi(Ind_200)*i_mi(1) + phi(Ind_110)*i_mi(2) + &
                        phi(Ind_101)*i_mi(3))
                    g_ind(2) = g_ind(2) + 0.5d0* &
                        (phi(Ind_110)*i_mi(1) + phi(Ind_020)*i_mi(2) + &
                        phi(Ind_011)*i_mi(3))
                    g_ind(3) = g_ind(3) + 0.5d0* &
                        (phi(Ind_101)*i_mi(1) + phi(Ind_011)*i_mi(2) + &
                        phi(Ind_002)*i_mi(3))
                end if !is_polarizable(i)
                ! get the field at j due to permanent + induced mpoles of i
                ! for torque at j
                ! electrostatic potential at j due to permanent + induced mpoles at i
                ! and derivatives of that with respect to r_j
                ! minus signs due to deriv of r_j-r_i wrt r_i entering into potential
                phi(Ind_000) = Rn_5(Ind_000)*tmi(Ind_000) - Rn_5(Ind_100)*tmi(Ind_100) - &
                    Rn_5(Ind_010)*tmi(Ind_010) - Rn_5(Ind_001)*tmi(Ind_001) + &
                    Rn_5(Ind_200)*tmi(Ind_200) + Rn_5(Ind_020)*tmi(Ind_020) + &
                    Rn_5(Ind_002)*tmi(Ind_002) + Rn_5(Ind_110)*tmi(Ind_110) + &
                    Rn_5(Ind_101)*tmi(Ind_101) + Rn_5(Ind_011)*tmi(Ind_011)
                phi(Ind_100) = Rn_5(Ind_100)*tmi(Ind_000) - Rn_5(Ind_200)*tmi(Ind_100) - &
                    Rn_5(Ind_110)*tmi(Ind_010) - Rn_5(Ind_101)*tmi(Ind_001) + &
                    Rn_5(Ind_300)*tmi(Ind_200) + Rn_5(Ind_120)*tmi(Ind_020) + &
                    Rn_5(Ind_102)*tmi(Ind_002) + Rn_5(Ind_210)*tmi(Ind_110) + &
                    Rn_5(Ind_201)*tmi(Ind_101) + Rn_5(Ind_111)*tmi(Ind_011)
                phi(Ind_010) = Rn_5(Ind_010)*tmi(Ind_000) - Rn_5(Ind_110)*tmi(Ind_100) - &
                    Rn_5(Ind_020)*tmi(Ind_010) - Rn_5(Ind_011)*tmi(Ind_001) + &
                    Rn_5(Ind_210)*tmi(Ind_200) + Rn_5(Ind_030)*tmi(Ind_020) + &
                    Rn_5(Ind_012)*tmi(Ind_002) + Rn_5(Ind_120)*tmi(Ind_110) + &
                    Rn_5(Ind_111)*tmi(Ind_101) + Rn_5(Ind_021)*tmi(Ind_011)
                phi(Ind_001) = Rn_5(Ind_001)*tmi(Ind_000) - Rn_5(Ind_101)*tmi(Ind_100) - &
                    Rn_5(Ind_011)*tmi(Ind_010) - Rn_5(Ind_002)*tmi(Ind_001) + &
                    Rn_5(Ind_201)*tmi(Ind_200) + Rn_5(Ind_021)*tmi(Ind_020) + &
                    Rn_5(Ind_003)*tmi(Ind_002) + Rn_5(Ind_111)*tmi(Ind_110) + &
                    Rn_5(Ind_102)*tmi(Ind_101) + Rn_5(Ind_012)*tmi(Ind_011)
                phi(Ind_200) = Rn_5(Ind_200)*tmi(Ind_000) - Rn_5(Ind_300)*tmi(Ind_100) - &
                    Rn_5(Ind_210)*tmi(Ind_010) - Rn_5(Ind_201)*tmi(Ind_001) + &
                    Rn_5(Ind_400)*tmi(Ind_200) + Rn_5(Ind_220)*tmi(Ind_020) + &
                    Rn_5(Ind_202)*tmi(Ind_002) + Rn_5(Ind_310)*tmi(Ind_110) + &
                    Rn_5(Ind_301)*tmi(Ind_101) + Rn_5(Ind_211)*tmi(Ind_011)
                phi(Ind_020) = Rn_5(Ind_020)*tmi(Ind_000) - Rn_5(Ind_120)*tmi(Ind_100) - &
                    Rn_5(Ind_030)*tmi(Ind_010) - Rn_5(Ind_021)*tmi(Ind_001) + &
                    Rn_5(Ind_220)*tmi(Ind_200) + Rn_5(Ind_040)*tmi(Ind_020) + &
                    Rn_5(Ind_022)*tmi(Ind_002) + Rn_5(Ind_130)*tmi(Ind_110) + &
                    Rn_5(Ind_121)*tmi(Ind_101) + Rn_5(Ind_031)*tmi(Ind_011)
                phi(Ind_002) = Rn_5(Ind_002)*tmi(Ind_000) - Rn_5(Ind_102)*tmi(Ind_100) - &
                    Rn_5(Ind_012)*tmi(Ind_010) - Rn_5(Ind_003)*tmi(Ind_001) + &
                    Rn_5(Ind_202)*tmi(Ind_200) + Rn_5(Ind_022)*tmi(Ind_020) + &
                    Rn_5(Ind_004)*tmi(Ind_002) + Rn_5(Ind_112)*tmi(Ind_110) + &
                    Rn_5(Ind_103)*tmi(Ind_101) + Rn_5(Ind_013)*tmi(Ind_011)
                phi(Ind_110) = Rn_5(Ind_110)*tmi(Ind_000) - Rn_5(Ind_210)*tmi(Ind_100) - &
                    Rn_5(Ind_120)*tmi(Ind_010) - Rn_5(Ind_111)*tmi(Ind_001) + &
                    Rn_5(Ind_310)*tmi(Ind_200) + Rn_5(Ind_130)*tmi(Ind_020) + &
                    Rn_5(Ind_112)*tmi(Ind_002) + Rn_5(Ind_220)*tmi(Ind_110) + &
                    Rn_5(Ind_211)*tmi(Ind_101) + Rn_5(Ind_121)*tmi(Ind_011)
                phi(Ind_101) = Rn_5(Ind_101)*tmi(Ind_000) - Rn_5(Ind_201)*tmi(Ind_100) - &
                    Rn_5(Ind_111)*tmi(Ind_010) - Rn_5(Ind_102)*tmi(Ind_001) + &
                    Rn_5(Ind_301)*tmi(Ind_200) + Rn_5(Ind_121)*tmi(Ind_020) + &
                    Rn_5(Ind_103)*tmi(Ind_002) + Rn_5(Ind_211)*tmi(Ind_110) + &
                    Rn_5(Ind_202)*tmi(Ind_101) + Rn_5(Ind_112)*tmi(Ind_011)
                phi(Ind_011) = Rn_5(Ind_011)*tmi(Ind_000) - Rn_5(Ind_111)*tmi(Ind_100) - &
                    Rn_5(Ind_021)*tmi(Ind_010) - Rn_5(Ind_012)*tmi(Ind_001) + &
                    Rn_5(Ind_211)*tmi(Ind_200) + Rn_5(Ind_031)*tmi(Ind_020) + &
                    Rn_5(Ind_013)*tmi(Ind_002) + Rn_5(Ind_121)*tmi(Ind_110) + &
                    Rn_5(Ind_112)*tmi(Ind_101) + Rn_5(Ind_022)*tmi(Ind_011)
                ! torque field at j due to permanent + induced mpoles of i
                do jj = 1, 10
                    torque_field(jj, j) = torque_field(jj, j) + phi(jj)
                end do
                if (is_polarizable(j)) then
                    ! phi array (electrostatic potential at i due to j induced moments
                    ! and derivs of that esp wrt r_i )
                    ! first that due to ind_dip_d for energy contribution
                    ! minus signs arise due to derivs of r_j - r_i wrt r_i
                    phi(Ind_000) = Rn_5(Ind_100)*i_dj(1) + Rn_5(Ind_010)*i_dj(2) + &
                        Rn_5(Ind_001)*i_dj(3)
                    phi(Ind_100) = -(Rn_5(Ind_200)*i_dj(1) + Rn_5(Ind_110)*i_dj(2) + &
                        Rn_5(Ind_101)*i_dj(3))
                    phi(Ind_010) = -(Rn_5(Ind_110)*i_dj(1) + Rn_5(Ind_020)*i_dj(2) + &
                        Rn_5(Ind_011)*i_dj(3))
                    phi(Ind_001) = -(Rn_5(Ind_101)*i_dj(1) + Rn_5(Ind_011)*i_dj(2) + &
                        Rn_5(Ind_002)*i_dj(3))
                    phi(Ind_200) = Rn_5(Ind_300)*i_dj(1) + Rn_5(Ind_210)*i_dj(2) + &
                        Rn_5(Ind_201)*i_dj(3)
                    phi(Ind_020) = Rn_5(Ind_120)*i_dj(1) + Rn_5(Ind_030)*i_dj(2) + &
                        Rn_5(Ind_021)*i_dj(3)
                    phi(Ind_002) = Rn_5(Ind_102)*i_dj(1) + Rn_5(Ind_012)*i_dj(2) + &
                        Rn_5(Ind_003)*i_dj(3)
                    phi(Ind_110) = Rn_5(Ind_210)*i_dj(1) + Rn_5(Ind_120)*i_dj(2) + &
                        Rn_5(Ind_111)*i_dj(3)
                    phi(Ind_101) = Rn_5(Ind_201)*i_dj(1) + Rn_5(Ind_111)*i_dj(2) + &
                        Rn_5(Ind_102)*i_dj(3)
                    phi(Ind_011) = Rn_5(Ind_111)*i_dj(1) + Rn_5(Ind_021)*i_dj(2) + &
                        Rn_5(Ind_012)*i_dj(3)
                    e_ind = e_ind + 0.5d0* &
                        (phi(Ind_000)*gmi(Ind_000) + phi(Ind_100)*gmi(Ind_100) + &
                        phi(Ind_010)*gmi(Ind_010) + phi(Ind_001)*gmi(Ind_001) + &
                        phi(Ind_200)*gmi(Ind_200) + phi(Ind_020)*gmi(Ind_020) + &
                        phi(Ind_002)*gmi(Ind_002) + phi(Ind_110)*gmi(Ind_110) + &
                        phi(Ind_101)*gmi(Ind_101) + phi(Ind_011)*gmi(Ind_011))
                    ! next that due to ind_dip_d+ind_dip_p for force contribution
                    phi(Ind_000) = Rn_5(Ind_100)*i_mj(1) + Rn_5(Ind_010)*i_mj(2) + &
                        Rn_5(Ind_001)*i_mj(3)
                    phi(Ind_100) = -(Rn_5(Ind_200)*i_mj(1) + Rn_5(Ind_110)*i_mj(2) + &
                        Rn_5(Ind_101)*i_mj(3))
                    phi(Ind_010) = -(Rn_5(Ind_110)*i_mj(1) + Rn_5(Ind_020)*i_mj(2) + &
                        Rn_5(Ind_011)*i_mj(3))
                    phi(Ind_001) = -(Rn_5(Ind_101)*i_mj(1) + Rn_5(Ind_011)*i_mj(2) + &
                        Rn_5(Ind_002)*i_mj(3))
                    phi(Ind_200) = Rn_5(Ind_300)*i_mj(1) + Rn_5(Ind_210)*i_mj(2) + &
                        Rn_5(Ind_201)*i_mj(3)
                    phi(Ind_020) = Rn_5(Ind_120)*i_mj(1) + Rn_5(Ind_030)*i_mj(2) + &
                        Rn_5(Ind_021)*i_mj(3)
                    phi(Ind_002) = Rn_5(Ind_102)*i_mj(1) + Rn_5(Ind_012)*i_mj(2) + &
                        Rn_5(Ind_003)*i_mj(3)
                    phi(Ind_110) = Rn_5(Ind_210)*i_mj(1) + Rn_5(Ind_120)*i_mj(2) + &
                        Rn_5(Ind_111)*i_mj(3)
                    phi(Ind_101) = Rn_5(Ind_201)*i_mj(1) + Rn_5(Ind_111)*i_mj(2) + &
                        Rn_5(Ind_102)*i_mj(3)
                    phi(Ind_011) = Rn_5(Ind_111)*i_mj(1) + Rn_5(Ind_021)*i_mj(2) + &
                        Rn_5(Ind_012)*i_mj(3)
                    phi(Ind_300) = -(Rn_5(Ind_400)*i_mj(1) + Rn_5(Ind_310)*i_mj(2) + &
                        Rn_5(Ind_301)*i_mj(3))
                    phi(Ind_030) = -(Rn_5(Ind_130)*i_mj(1) + Rn_5(Ind_040)*i_mj(2) + &
                        Rn_5(Ind_031)*i_mj(3))
                    phi(Ind_003) = -(Rn_5(Ind_103)*i_mj(1) + Rn_5(Ind_013)*i_mj(2) + &
                        Rn_5(Ind_004)*i_mj(3))
                    phi(Ind_210) = -(Rn_5(Ind_310)*i_mj(1) + Rn_5(Ind_220)*i_mj(2) + &
                        Rn_5(Ind_211)*i_mj(3))
                    phi(Ind_201) = -(Rn_5(Ind_301)*i_mj(1) + Rn_5(Ind_211)*i_mj(2) + &
                        Rn_5(Ind_202)*i_mj(3))
                    phi(Ind_120) = -(Rn_5(Ind_220)*i_mj(1) + Rn_5(Ind_130)*i_mj(2) + &
                        Rn_5(Ind_121)*i_mj(3))
                    phi(Ind_021) = -(Rn_5(Ind_121)*i_mj(1) + Rn_5(Ind_031)*i_mj(2) + &
                        Rn_5(Ind_022)*i_mj(3))
                    phi(Ind_102) = -(Rn_5(Ind_202)*i_mj(1) + Rn_5(Ind_112)*i_mj(2) + &
                        Rn_5(Ind_103)*i_mj(3))
                    phi(Ind_012) = -(Rn_5(Ind_112)*i_mj(1) + Rn_5(Ind_022)*i_mj(2) + &
                        Rn_5(Ind_013)*i_mj(3))
                    phi(Ind_111) = -(Rn_5(Ind_211)*i_mj(1) + Rn_5(Ind_121)*i_mj(2) + &
                        Rn_5(Ind_112)*i_mj(3))
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
                    ! torque field at i due to induced moments of j
                    !  note the factor of 1/2
                    do jj = 1, 10
                        torque_field(jj, i) = torque_field(jj, i) + 0.5d0*phi(jj)
                    end do
                    if (is_polarizable(i)) then !i,j both polarizable
                        ! phi array (electrostatic potential at i due to j induced moments
                        ! and derivs of that esp wrt r_i )
                        phi(Ind_200) = Rn_5(Ind_300)*i_dj(1) + Rn_5(Ind_210)*i_dj(2) + &
                            Rn_5(Ind_201)*i_dj(3)
                        phi(Ind_020) = Rn_5(Ind_120)*i_dj(1) + Rn_5(Ind_030)*i_dj(2) + &
                            Rn_5(Ind_021)*i_dj(3)
                        phi(Ind_002) = Rn_5(Ind_102)*i_dj(1) + Rn_5(Ind_012)*i_dj(2) + &
                            Rn_5(Ind_003)*i_dj(3)
                        phi(Ind_110) = Rn_5(Ind_210)*i_dj(1) + Rn_5(Ind_120)*i_dj(2) + &
                            Rn_5(Ind_111)*i_dj(3)
                        phi(Ind_101) = Rn_5(Ind_201)*i_dj(1) + Rn_5(Ind_111)*i_dj(2) + &
                            Rn_5(Ind_102)*i_dj(3)
                        phi(Ind_011) = Rn_5(Ind_111)*i_dj(1) + Rn_5(Ind_021)*i_dj(2) + &
                            Rn_5(Ind_012)*i_dj(3)
                        g_ind(1) = g_ind(1) + 0.5d0* &
                            (phi(Ind_200)*i_pi(1) + phi(Ind_110)*i_pi(2) + &
                            phi(Ind_101)*i_pi(3))
                        g_ind(2) = g_ind(2) + 0.5d0* &
                            (phi(Ind_110)*i_pi(1) + phi(Ind_020)*i_pi(2) + &
                            phi(Ind_011)*i_pi(3))
                        g_ind(3) = g_ind(3) + 0.5d0* &
                            (phi(Ind_101)*i_pi(1) + phi(Ind_011)*i_pi(2) + &
                            phi(Ind_002)*i_pi(3))
                        phi(Ind_200) = Rn_5(Ind_300)*i_pj(1) + Rn_5(Ind_210)*i_pj(2) + &
                            Rn_5(Ind_201)*i_pj(3)
                        phi(Ind_020) = Rn_5(Ind_120)*i_pj(1) + Rn_5(Ind_030)*i_pj(2) + &
                            Rn_5(Ind_021)*i_pj(3)
                        phi(Ind_002) = Rn_5(Ind_102)*i_pj(1) + Rn_5(Ind_012)*i_pj(2) + &
                            Rn_5(Ind_003)*i_pj(3)
                        phi(Ind_110) = Rn_5(Ind_210)*i_pj(1) + Rn_5(Ind_120)*i_pj(2) + &
                            Rn_5(Ind_111)*i_pj(3)
                        phi(Ind_101) = Rn_5(Ind_201)*i_pj(1) + Rn_5(Ind_111)*i_pj(2) + &
                            Rn_5(Ind_102)*i_pj(3)
                        phi(Ind_011) = Rn_5(Ind_111)*i_pj(1) + Rn_5(Ind_021)*i_pj(2) + &
                            Rn_5(Ind_012)*i_pj(3)
                        g_ind(1) = g_ind(1) + 0.5d0* &
                            (phi(Ind_200)*i_di(1) + phi(Ind_110)*i_di(2) + &
                            phi(Ind_101)*i_di(3))
                        g_ind(2) = g_ind(2) + 0.5d0* &
                            (phi(Ind_110)*i_di(1) + phi(Ind_020)*i_di(2) + &
                            phi(Ind_011)*i_di(3))
                        g_ind(3) = g_ind(3) + 0.5d0* &
                            (phi(Ind_101)*i_di(1) + phi(Ind_011)*i_di(2) + &
                            phi(Ind_002)*i_di(3))
                    end if !is_polarizable(i)
                end if !is_polarizable(j)
                if ((delr2 < ee_damped_cut2) .and. &
                    (is_polarizable(i) .or. is_polarizable(j))) then
                    !-------------------------------------------------------
                    ! McMurchie-Davidson holds for damped tensor as well---in fact,
                    ! BD below satisfies BD(n+1) = (1/r)(d/dr)BD(n)
                    ! RD(0,0,0,n) = BD(n), n=1,2,3
                    ! RD(t,u,v,n) = (d/dx)^t (d/dy)^u (d/dz)^v R(0,0,0,n)
                    !-------------------------------------------------------
                    delr3inv = delr2inv/delr
                    delr5inv = delr3inv*delr2inv
                    delr7inv = delr5inv*delr2inv
                    delr9inv = delr7inv*delr2inv
                    expon = asq*delr2*delr*sq_polinv(j)
                    expo = exp(-expon)
                    ! clam3 = 1.d0-lam3, clam5 = 1.d0-lam5 etc. where
                    ! lam is from ponder's paper
                    clam3 = expo
                    clam5 = (1.d0 + expon)*expo
                    clam7 = (1.d0 + expon + const1*expon**2)*expo
                    clam9 = (1.d0 + expon + const2*expon**2 + const3*expon**3)*expo
                    BD(1) = -clam3*delr3inv
                    BD(2) = 3.d0*clam5*delr5inv
                    BD(3) = -15.d0*clam7*delr7inv
                    BD(4) = 105.d0*clam9*delr9inv
                    n = 4
                    Rn(Ind_000) = BD(n)
                    Rn_1(Ind_000) = BD(n - 1)
                    Rn_1(Ind_100) = delx*Rn(Ind_000)
                    Rn_1(Ind_010) = dely*Rn(Ind_000)
                    Rn_1(Ind_001) = delz*Rn(Ind_000)
                    Rn_2(Ind_000) = BD(n - 2)
                    Rn_2(Ind_100) = delx*Rn_1(Ind_000)
                    Rn_2(Ind_010) = dely*Rn_1(Ind_000)
                    Rn_2(Ind_001) = delz*Rn_1(Ind_000)
                    Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
                    Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
                    Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
                    Rn_2(Ind_110) = delx*Rn_1(Ind_010)
                    Rn_2(Ind_101) = delx*Rn_1(Ind_001)
                    Rn_2(Ind_011) = dely*Rn_1(Ind_001)
                    Rn_3(Ind_000) = BD(n - 3)
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
                        ! minus sign since we remove damped correction
                        e_ind = e_ind - 0.5d0* &
                            (phi(Ind_100)*i_di(1) + phi(Ind_010)*i_di(2) + &
                            phi(Ind_001)*i_di(3))
                        g_ind(1) = g_ind(1) - 0.5d0* &
                            (phi(Ind_200)*i_mi(1) + phi(Ind_110)*i_mi(2) + &
                            phi(Ind_101)*i_mi(3))
                        g_ind(2) = g_ind(2) - 0.5d0* &
                            (phi(Ind_110)*i_mi(1) + phi(Ind_020)*i_mi(2) + &
                            phi(Ind_011)*i_mi(3))
                        g_ind(3) = g_ind(3) - 0.5d0* &
                            (phi(Ind_101)*i_mi(1) + phi(Ind_011)*i_mi(2) + &
                            phi(Ind_002)*i_mi(3))
                        ! next do torque field at j due to induced at i
                        ! potential at j do to dipoles negative due to derivs wrt r_i
                        ! higher order are derivs of pot wrt r_j so no sign change
                        phi(Ind_000) = -(Rn_4(Ind_100)*i_mi(1) + Rn_4(Ind_010)*i_mi(2) + &
                            Rn_4(Ind_001)*i_mi(3))
                        phi(Ind_100) = -(Rn_4(Ind_200)*i_mi(1) + Rn_4(Ind_110)*i_mi(2) + &
                            Rn_4(Ind_101)*i_mi(3))
                        phi(Ind_010) = -(Rn_4(Ind_110)*i_mi(1) + Rn_4(Ind_020)*i_mi(2) + &
                            Rn_4(Ind_011)*i_mi(3))
                        phi(Ind_001) = -(Rn_4(Ind_101)*i_mi(1) + Rn_4(Ind_011)*i_mi(2) + &
                            Rn_4(Ind_002)*i_mi(3))
                        phi(Ind_200) = -(Rn_4(Ind_300)*i_mi(1) + Rn_4(Ind_210)*i_mi(2) + &
                            Rn_4(Ind_201)*i_mi(3))
                        phi(Ind_020) = -(Rn_4(Ind_120)*i_mi(1) + Rn_4(Ind_030)*i_mi(2) + &
                            Rn_4(Ind_021)*i_mi(3))
                        phi(Ind_002) = -(Rn_4(Ind_102)*i_mi(1) + Rn_4(Ind_012)*i_mi(2) + &
                            Rn_4(Ind_003)*i_mi(3))
                        phi(Ind_110) = -(Rn_4(Ind_210)*i_mi(1) + Rn_4(Ind_120)*i_mi(2) + &
                            Rn_4(Ind_111)*i_mi(3))
                        phi(Ind_101) = -(Rn_4(Ind_201)*i_mi(1) + Rn_4(Ind_111)*i_mi(2) + &
                            Rn_4(Ind_102)*i_mi(3))
                        phi(Ind_011) = -(Rn_4(Ind_111)*i_mi(1) + Rn_4(Ind_021)*i_mi(2) + &
                            Rn_4(Ind_012)*i_mi(3))
                        ! torque field at j due to induced moments of i
                        !  note the factor of 1/2 (term ~ 1/2 the total induced moment)
                        ! the minus sign is since we remove damped contributions
                        do jj = 1, 10
                            torque_field(jj, j) = torque_field(jj, j) - 0.5d0*phi(jj)
                        end do
                    end if !is_polarizable(i)
                    if (is_polarizable(j)) then
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
                        ! minus sign since we remove damped contributions
                        e_ind = e_ind - 0.5d0* &
                            (phi(Ind_000)*gmi(Ind_000) + phi(Ind_100)*gmi(Ind_100) + &
                            phi(Ind_010)*gmi(Ind_010) + phi(Ind_001)*gmi(Ind_001) + &
                            phi(Ind_200)*gmi(Ind_200) + phi(Ind_020)*gmi(Ind_020) + &
                            phi(Ind_002)*gmi(Ind_002) + phi(Ind_110)*gmi(Ind_110) + &
                            phi(Ind_101)*gmi(Ind_101) + phi(Ind_011)*gmi(Ind_011))
                        ! next that due to ind_dip_d+ind_dip_p for force contribution
                        phi(Ind_000) = Rn_4(Ind_100)*i_mj(1) + Rn_4(Ind_010)*i_mj(2) + &
                            Rn_4(Ind_001)*i_mj(3)
                        phi(Ind_100) = -(Rn_4(Ind_200)*i_mj(1) + Rn_4(Ind_110)*i_mj(2) + &
                            Rn_4(Ind_101)*i_mj(3))
                        phi(Ind_010) = -(Rn_4(Ind_110)*i_mj(1) + Rn_4(Ind_020)*i_mj(2) + &
                            Rn_4(Ind_011)*i_mj(3))
                        phi(Ind_001) = -(Rn_4(Ind_101)*i_mj(1) + Rn_4(Ind_011)*i_mj(2) + &
                            Rn_4(Ind_002)*i_mj(3))
                        phi(Ind_200) = Rn_4(Ind_300)*i_mj(1) + Rn_4(Ind_210)*i_mj(2) + &
                            Rn_4(Ind_201)*i_mj(3)
                        phi(Ind_020) = Rn_4(Ind_120)*i_mj(1) + Rn_4(Ind_030)*i_mj(2) + &
                            Rn_4(Ind_021)*i_mj(3)
                        phi(Ind_002) = Rn_4(Ind_102)*i_mj(1) + Rn_4(Ind_012)*i_mj(2) + &
                            Rn_4(Ind_003)*i_mj(3)
                        phi(Ind_110) = Rn_4(Ind_210)*i_mj(1) + Rn_4(Ind_120)*i_mj(2) + &
                            Rn_4(Ind_111)*i_mj(3)
                        phi(Ind_101) = Rn_4(Ind_201)*i_mj(1) + Rn_4(Ind_111)*i_mj(2) + &
                            Rn_4(Ind_102)*i_mj(3)
                        phi(Ind_011) = Rn_4(Ind_111)*i_mj(1) + Rn_4(Ind_021)*i_mj(2) + &
                            Rn_4(Ind_012)*i_mj(3)
                        phi(Ind_300) = -(Rn_4(Ind_400)*i_mj(1) + Rn_4(Ind_310)*i_mj(2) + &
                            Rn_4(Ind_301)*i_mj(3))
                        phi(Ind_030) = -(Rn_4(Ind_130)*i_mj(1) + Rn_4(Ind_040)*i_mj(2) + &
                            Rn_4(Ind_031)*i_mj(3))
                        phi(Ind_003) = -(Rn_4(Ind_103)*i_mj(1) + Rn_4(Ind_013)*i_mj(2) + &
                            Rn_4(Ind_004)*i_mj(3))
                        phi(Ind_210) = -(Rn_4(Ind_310)*i_mj(1) + Rn_4(Ind_220)*i_mj(2) + &
                            Rn_4(Ind_211)*i_mj(3))
                        phi(Ind_201) = -(Rn_4(Ind_301)*i_mj(1) + Rn_4(Ind_211)*i_mj(2) + &
                            Rn_4(Ind_202)*i_mj(3))
                        phi(Ind_120) = -(Rn_4(Ind_220)*i_mj(1) + Rn_4(Ind_130)*i_mj(2) + &
                            Rn_4(Ind_121)*i_mj(3))
                        phi(Ind_021) = -(Rn_4(Ind_121)*i_mj(1) + Rn_4(Ind_031)*i_mj(2) + &
                            Rn_4(Ind_022)*i_mj(3))
                        phi(Ind_102) = -(Rn_4(Ind_202)*i_mj(1) + Rn_4(Ind_112)*i_mj(2) + &
                            Rn_4(Ind_103)*i_mj(3))
                        phi(Ind_012) = -(Rn_4(Ind_112)*i_mj(1) + Rn_4(Ind_022)*i_mj(2) + &
                            Rn_4(Ind_013)*i_mj(3))
                        phi(Ind_111) = -(Rn_4(Ind_211)*i_mj(1) + Rn_4(Ind_121)*i_mj(2) + &
                            Rn_4(Ind_112)*i_mj(3))
                        ! minus sign since we remove damped contributions
                        g_ind(1) = g_ind(1) - 0.5d0* &
                            (phi(Ind_100)*gmi(Ind_000) + phi(Ind_200)*gmi(Ind_100) + &
                            phi(Ind_110)*gmi(Ind_010) + phi(Ind_101)*gmi(Ind_001) + &
                            phi(Ind_300)*gmi(Ind_200) + phi(Ind_120)*gmi(Ind_020) + &
                            phi(Ind_102)*gmi(Ind_002) + phi(Ind_210)*gmi(Ind_110) + &
                            phi(Ind_201)*gmi(Ind_101) + phi(Ind_111)*gmi(Ind_011))
                        g_ind(2) = g_ind(2) - 0.5d0* &
                            (phi(Ind_010)*gmi(Ind_000) + phi(Ind_110)*gmi(Ind_100) + &
                            phi(Ind_020)*gmi(Ind_010) + phi(Ind_011)*gmi(Ind_001) + &
                            phi(Ind_210)*gmi(Ind_200) + phi(Ind_030)*gmi(Ind_020) + &
                            phi(Ind_012)*gmi(Ind_002) + phi(Ind_120)*gmi(Ind_110) + &
                            phi(Ind_111)*gmi(Ind_101) + phi(Ind_021)*gmi(Ind_011))
                        g_ind(3) = g_ind(3) - 0.5d0* &
                            (phi(Ind_001)*gmi(Ind_000) + phi(Ind_101)*gmi(Ind_100) + &
                            phi(Ind_011)*gmi(Ind_010) + phi(Ind_002)*gmi(Ind_001) + &
                            phi(Ind_201)*gmi(Ind_200) + phi(Ind_021)*gmi(Ind_020) + &
                            phi(Ind_003)*gmi(Ind_002) + phi(Ind_111)*gmi(Ind_110) + &
                            phi(Ind_102)*gmi(Ind_101) + phi(Ind_012)*gmi(Ind_011))
                        ! torque field at i due to induced moments of j
                        !  note the factor of 1/2 (term ~ 1/2 the total induced moment)
                        ! the minus sign is since we remove damped contributions
                        do jj = 1, 10
                            torque_field(jj, i) = torque_field(jj, i) - 0.5d0*phi(jj)
                        end do
                        if (is_polarizable(i)) then ! i,j both polarizable
                            ! phi array (electrostatic pot at i due to j induced moments
                            ! and derivs of that esp wrt r_i )
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
                            ! minus sign since we remove damped contributions
                            g_ind(1) = g_ind(1) - 0.5d0* &
                                (phi(Ind_200)*i_pi(1) + phi(Ind_110)*i_pi(2) + &
                                phi(Ind_101)*i_pi(3))
                            g_ind(2) = g_ind(2) - 0.5d0* &
                                (phi(Ind_110)*i_pi(1) + phi(Ind_020)*i_pi(2) + &
                                phi(Ind_011)*i_pi(3))
                            g_ind(3) = g_ind(3) - 0.5d0* &
                                (phi(Ind_101)*i_pi(1) + phi(Ind_011)*i_pi(2) + &
                                phi(Ind_002)*i_pi(3))
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
                            ! minus sign since we remove damped contributions
                            g_ind(1) = g_ind(1) - 0.5d0* &
                                (phi(Ind_200)*i_di(1) + phi(Ind_110)*i_di(2) + &
                                phi(Ind_101)*i_di(3))
                            g_ind(2) = g_ind(2) - 0.5d0* &
                                (phi(Ind_110)*i_di(1) + phi(Ind_020)*i_di(2) + &
                                phi(Ind_011)*i_di(3))
                            g_ind(3) = g_ind(3) - 0.5d0* &
                                (phi(Ind_101)*i_di(1) + phi(Ind_011)*i_di(2) + &
                                phi(Ind_002)*i_di(3))
                        end if !is_polarizable(i)
                    end if !is_polarizable(j)
                end if !(delr2 < ee_damped_cut2) .and. &
                !(is_polarizable(i) .or. is_polarizable(j))
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
            end if !( delr2 < ee_dsum_cut2 )then
        end do !m = 1,numlist
        virial(1, 1) = virial(1, 1) + vxx
        virial(1, 2) = virial(1, 2) + half*(vxy + vyx)
        virial(1, 3) = virial(1, 3) + half*(vxz + vzx)
        virial(2, 1) = virial(2, 1) + half*(vxy + vyx)
        virial(2, 2) = virial(2, 2) + vyy
        virial(2, 3) = virial(2, 3) + half*(vyz + vzy)
        virial(3, 1) = virial(3, 1) + half*(vxz + vzx)
        virial(3, 2) = virial(3, 2) + half*(vyz + vzy)
        virial(3, 3) = virial(3, 3) + vzz
    end subroutine AM_DIRECT_ene_force_i
!-------------------------------------------------------
end module amoeba_direct
