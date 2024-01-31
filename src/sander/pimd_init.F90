#include "assert.fh"
#include "dprec.fh"

subroutine pimd_init(natom, mass, winv, v, pimdtype)
    use constants, only : kB
    use pimd_vars, only : nbead, natomCL, nbead_inv, tempbuf, &
        PRIPIMD, NMPIMD, CMD, RPMD, &
        equal_part, x_centroid, real_mass, nrg_all, &
        cartpos, cartvel
    use cmd_vars, only : file_nrg_cmd, file_pos_cmd, file_vel_cmd &
        , cart_to_nmode, nmode_to_cart, fcart_to_fnmode &
        , lambda_nmode, omega_nmode, adiab_param
#if defined(LES) && defined(MPI)
    use evb_pimd, only : bead_dcrypt
    use miller, only : ti_mass, nti_mass, do_ti_mass
#endif /* LES:MPI */
    use full_pimd_vars, only : mybeadid
    implicit none
# include "md.h"
# include "extra.h"
#ifdef LES
# include "les.h"
#endif

    integer iatom, natom, pimdtype, ierr
    _REAL_ mass(*), winv(*), kT, beta
    _REAL_ v(3, natom), sigma, dummy
    integer :: idim

#if defined(LES) && defined(MPI)
    integer :: m, n, mm, nn
    _REAL_  :: mass_old, mass_new
#endif /* LES:MPI */

    allocate (tempbuf(3, natom), stat=ierr)
    REQUIRE(ierr .eq. 0)

    allocate (x_centroid(3, natom), stat=ierr)
    REQUIRE(ierr .eq. 0)

    allocate (real_mass(1:natom), stat=ierr)
    REQUIRE(ierr == 0)

#ifdef LES
    call part_pimd_init()
#else
    call full_pimd_init()
#endif

#if defined(LES) && defined(MPI)

!  +---------------------------------------------------------------+
!  |  Modify mass for thermodynamic integration                    |
!  |                                                               |
!  |  M = (1-lambda) M_i + lambda M_f                              |
!  +---------------------------------------------------------------+

    if (do_ti_mass .and. ievb /= 0) then
        do n = 1, nti_mass
            mass_old = mass(ti_mass(n)%atm_ndx)
            if (mass_old /= ti_mass(n)%init) then
                write (6, '(A)') "ERROR: User requests TI by mass but the initial " &
                    //"mass is not equal to the mass in the topology file."
                write (6, '(5X,A,F10.5)') "Mass from topology file: ", mass_old
                write (6, '(5X,A,F10.5)') "Initial mass: ", ti_mass(n)%init
                call mexit(6, 1)
            end if

            mass_new = ti_mass(n)%init + ti_mass(n)%lambda &
                *(ti_mass(n)%final - ti_mass(n)%init)

            write (6, '(A)') "Performing thermodynamic integration by mass"
            write (6, '((5X,(A,I10),A,2X,3(A,F10.5)))') 'ti_mass(', ti_mass(n)%atm_ndx &
                , ') ::', 'lambda = ', ti_mass(n)%lambda &
                , ', M_i = ', ti_mass(n)%init, ', M_f = ', ti_mass(n)%final

            write (6, 1000) 'ATOM INDEX', 'ORIGINAL MASS', 'EFFECTIVE MASS'

            mm = ti_mass(n)%atm_ndx

            do nn = 1, nbead
                mass(bead_dcrypt(mm, nn)) = mass_new
                write (6, 2000) bead_dcrypt(mm, nn), mass_old, mass(bead_dcrypt(mm, nn))
            end do
        end do
    end if

1000 format(5X, 3(A, 4X))
2000 format(5X, I8, 4X, 2(F10.5, 4X))

#endif /* LES:MPI */

    real_mass(1:natom) = mass(1:natom)

    kT = kB*temp0
    beta = 1.0d0/kT
    equal_part = 1.5d0*dble(natomCL)/beta

    allocate (nrg_all(nbead), stat=ierr)
    REQUIRE(ierr .eq. 0)

    nbead_inv = 1.d0/nbead

    if (pimdtype /= CMD .and. adiab_param /= 1.d0) then
        write (6, *) "Stop. For ipimd = ", pimdtype, " adiab_param must be = 1.d0"
        stop
    end if

    if (pimdtype .eq. CMD) then
        if (adiab_param > 1.d0) then
            write (6, *) "Stop. adiab_param must be <= 1.d0"
            stop
        else
            open (file_pos_cmd, file='CMD_position.dat')
            open (file_vel_cmd, file='CMD_velocity.dat')
            write (file_pos_cmd, *) 'CMD position'
            write (file_vel_cmd, *) 'CMD velocity'
        end if
    end if

    if (pimdtype .eq. PRIPIMD) then
#ifdef LES
        tmass = 0.d0
        do iatom = 1, natom
            if (cnum(iatom) .ne. 0) then
                mass(iatom) = mass(iatom)/dble(nbead)
                winv(iatom) = 1.d0/mass(iatom)
                tmass = tmass + mass(iatom)
            else
                tmass = tmass + mass(iatom)
            end if
        end do
        tmassinv = 1.d0/tmass
#else
        tmass = 0.d0
        do iatom = 1, natom
            tmass = tmass + mass(iatom)
            mass(iatom) = mass(iatom)/dble(nbead)
            winv(iatom) = 1.d0/mass(iatom)
        end do
        tmassinv = 1.d0/tmass
#endif
    end if

    if (pimdtype .eq. NMPIMD .or. pimdtype .eq. CMD) then

        allocate (lambda_nmode(nbead), stat=ierr)
        REQUIRE(ierr == 0)
        allocate (omega_nmode(nbead), stat=ierr)
        REQUIRE(ierr == 0)
        allocate (nmode_to_cart(nbead, nbead), stat=ierr)
        REQUIRE(ierr == 0)
        allocate (cart_to_nmode(nbead, nbead), stat=ierr)
        REQUIRE(ierr == 0)
        allocate (fcart_to_fnmode(nbead, nbead), stat=ierr)
        REQUIRE(ierr == 0)

        allocate (cartpos(3, natom), stat=ierr)
        REQUIRE(ierr .eq. 0)

        allocate (cartvel(3, natom), stat=ierr)
        REQUIRE(ierr .eq. 0)

        call transformation_matrix(nbead)

#ifdef LES
        tmass = 0.d0
        do iatom = 1, natom
            if (cnum(iatom) .ne. 0) then
                tmass = tmass + mass(iatom)/dble(nbead)
            else
                tmass = tmass + mass(iatom)
            end if
            if (cnum(iatom) .ne. 0 .and. cnum(iatom) .ne. 1) then
                mass(iatom) = &
                    mass(iatom)*lambda_nmode(cnum(iatom))*adiab_param*adiab_param
            end if
        end do
        tmassinv = 1.d0/tmass
#else
        if (mybeadid > 1) then
            mass(1:natom) = &
                mass(1:natom)*lambda_nmode(mybeadid)*adiab_param*adiab_param
            winv(1:natom) = 1.d0/mass(1:natom)
        end if
#endif
    end if

#ifdef LES
    do iatom = 1, natom
        if (pimdtype > 0 .and. irest == 0) then
            sigma = sqrt(kT/mass(iatom))
            do idim = 1, 3
                call gauss(0.d0, sigma, v(idim, iatom))
            end do
        end if
    end do
#else
    if (pimdtype > 0 .and. irest == 0) then
        do iatom = 1, (mybeadid - 1)*natom*3
            call gauss(0.d0, 1.d0, dummy)
        end do

        do iatom = 1, natom
            sigma = sqrt(kT/mass(iatom))
            do idim = 1, 3
                call gauss(0.d0, sigma, v(idim, iatom))
            end do
        end do

        do iatom = mybeadid*natom*3 + 1, nbead*natom*3
            call gauss(0.d0, 1.d0, dummy)
        end do
    end if
#endif

    if (master) then
        if (pimdtype == PRIPIMD) then
            write (6, '("THIS IS A PRIMITIVE PIMD RUN")')
        elseif (pimdtype == NMPIMD) then
            write (6, '("THIS IS A NORMAL-MODE PIMD RUN")')
        elseif (pimdtype == CMD) then
            write (6, '("THIS IS A CENTROID MD RUN")')
        elseif (pimdtype == RPMD) then
            write (6, '("THIS IS A RING-POLYMER MD RUN")')
        end if
        write (6, '("Parameters:")')
        write (6, '("number of beads           = ",i6)') nbead
        write (6, '("number of classical atoms = ",i6)') natomCL
        write (6, '("temperature (Kelvin)      = ",f8.2)') temp0
    end if

end subroutine pimd_init

subroutine pimd_finalize(pimdtype)
    use pimd_vars, only : nbead, nbead_inv, tempbuf, NMPIMD, CMD, equal_part, &
        x_centroid, real_mass, nrg_all, cartpos, cartvel, dmdlm, itimass
    use cmd_vars, only : file_nrg_cmd, file_pos_cmd, file_vel_cmd, cart_to_nmode, nmode_to_cart, fcart_to_fnmode, &
        lambda_nmode, omega_nmode, adiab_param

    implicit none
    integer pimdtype, ierr

#ifdef LES
    call part_pimd_finalize()
#else
    call full_pimd_finalize()
#endif
    deallocate (tempbuf, stat=ierr)
    REQUIRE(ierr == 0)
    deallocate (x_centroid, stat=ierr)
    REQUIRE(ierr == 0)
    deallocate (real_mass, stat=ierr)
    REQUIRE(ierr == 0)
    if (itimass > 0) then
        deallocate (dmdlm, stat=ierr)
        REQUIRE(ierr == 0)
    end if
    deallocate (nrg_all, stat=ierr)
    REQUIRE(ierr == 0)

    if (pimdtype .eq. NMPIMD .or. pimdtype .eq. CMD) then
        deallocate (lambda_nmode, stat=ierr)
        REQUIRE(ierr == 0)
        deallocate (omega_nmode, stat=ierr)
        REQUIRE(ierr == 0)
        deallocate (nmode_to_cart, stat=ierr)
        REQUIRE(ierr == 0)
        deallocate (cart_to_nmode, stat=ierr)
        REQUIRE(ierr == 0)
        deallocate (fcart_to_fnmode, stat=ierr)
        REQUIRE(ierr == 0)
        deallocate (cartpos, stat=ierr)
        REQUIRE(ierr == 0)
        deallocate (cartvel, stat=ierr)
        REQUIRE(ierr == 0)
    end if

end subroutine pimd_finalize

subroutine neb_init()
    use full_pimd_vars, only : mybeadid
    use neb_vars, only : neb_nbead, tangents, springforce, xprev, xnext, &
        neb_force, neb_nrg_all
    use file_io_dat
    implicit none
#include "md.h"
#include "memory.h"
#include "parallel.h"
    integer i, ierr

    neb_nbead = numgroup
    mybeadid = worldrank/(worldsize/neb_nbead) + 1

    allocate (neb_nrg_all(neb_nbead), stat=ierr)
    REQUIRE(ierr .eq. 0)

    allocate (springforce(3*natom), stat=ierr)
    REQUIRE(ierr .eq. 0)
    allocate (tangents(3*natom), stat=ierr)
    REQUIRE(ierr .eq. 0)
    allocate (xprev(3*natom), stat=ierr)
    REQUIRE(ierr .eq. 0)
    allocate (xnext(3*natom), stat=ierr)
    REQUIRE(ierr .eq. 0)
    allocate (neb_force(3*natom), stat=ierr) !Added by Sally Pias
    REQUIRE(ierr .eq. 0)

end subroutine neb_init

subroutine neb_finalize()
    use neb_vars, only : tangents, springforce, xprev, xnext, neb_force, neb_nrg_all

    deallocate (neb_nrg_all)
    deallocate (tangents)
    deallocate (springforce)
    deallocate (xprev)
    deallocate (xnext)
    deallocate (neb_force)

end subroutine neb_finalize

subroutine full_pimd_init()
    use pimd_vars, only : nbead, natomCL
    use full_pimd_vars, only : xall, mybeadid
    use file_io_dat
    implicit none
#  include "md.h"
#  include "memory.h"
#  include "parallel.h"
    integer ierr

    nbead = numgroup
    natomCL = natom
    mybeadid = worldrank/(worldsize/nbead) + 1
    if (worldrank .eq. 0) open (unit=pimd_unit, file=pimdout)

    REQUIRE(.not. allocated(xall))

    allocate (xall(3, natom, nbead), stat=ierr)
    REQUIRE(ierr .eq. 0)

end subroutine full_pimd_init

subroutine full_pimd_finalize()
    use full_pimd_vars, only : xall

    implicit none
    integer pimdtype, ierr

    deallocate (xall, stat=ierr)
    REQUIRE(ierr .eq. 0)

end subroutine full_pimd_finalize

#ifdef  LES
!------------------------------------------------------------------------------
subroutine part_pimd_init()
! initialize parameters

    use constants, only : kB

    use pimd_vars, only : natomCL, nbead
    use part_pimd_vars, only : frcx_copy, frcx_temp, &
        xrep, vrep, frep, ftmp, pimd_mmchg

    use qmmm_module, only : qmmm_nml, qmmm_struct
    !..................................................

    implicit none
# include "memory.h"
# include "les.h"
    integer :: iatom, alloc_error

    nbead = ncopy
    natomCL = 0
    do iatom = 1, natom
        if (cnum(iatom) .eq. 0 .or. cnum(iatom) .eq. 1) natomCL = natomCL + 1
    end do

    allocate (frcx_copy(3, ncopy), stat=alloc_error)
    REQUIRE(alloc_error == 0)
    allocate (frcx_temp(3, ncopy), stat=alloc_error)
    REQUIRE(alloc_error == 0)

    allocate (xrep(3*natomCL), stat=alloc_error)
    REQUIRE(alloc_error == 0)
    allocate (frep(3*natomCL), stat=alloc_error)
    REQUIRE(alloc_error == 0)
    allocate (vrep(3*natomCL), stat=alloc_error)
    REQUIRE(alloc_error == 0)
    allocate (ftmp(3, natom), stat=alloc_error)
    REQUIRE(alloc_error == 0)
    allocate (pimd_mmchg(natom), stat=alloc_error)
    REQUIRE(alloc_error == 0)

end subroutine part_pimd_init

subroutine part_pimd_finalize()
    use part_pimd_vars, only : frcx_copy, frcx_temp, massCL, &
        xrep, vrep, frep, ftmp, pimd_mmchg

    use qmmm_module, only : qmmm_nml, qmmm_struct

    implicit none
    integer alloc_error

    deallocate (frcx_copy)
    deallocate (frcx_temp)

    deallocate (frep, stat=alloc_error)
    REQUIRE(alloc_error == 0)

    deallocate (vrep, stat=alloc_error)
    REQUIRE(alloc_error == 0)

    deallocate (xrep, stat=alloc_error)
    REQUIRE(alloc_error == 0)

    deallocate (ftmp, stat=alloc_error)
    REQUIRE(alloc_error == 0)

    deallocate (pimd_mmchg, stat=alloc_error)
    REQUIRE(alloc_error == 0)

end subroutine part_pimd_finalize

subroutine part_setup_cnst_press_pimd(Nkt, tau_vol)
    use constants, only : kB, hbar
    use pimd_vars, only : nbead, natomCL
    use nose_hoover_module, only : M, thermo_lnv, v_lnv, f_lnv_p, c2_lnv, mass_lnv, Thermostat_init
    implicit none
#include "md.h"
#include "memory.h"
    _REAL_ :: tau_vol, kT, Nkt, omegaP
    integer :: j

    kT = kB*temp0

    Nkt = natomCL*kT

    omegaP = kT*sqrt(dble(nbead))/hbar

    c2_lnv = 1 + 1.0/dble(natomCL)

    v_lnv = 0.d0

    f_lnv_p = 0.d0

    mass_lnv = 3*(natomCL + 1)*kT*tau_vol*tau_vol

    thermo_lnv%Q(1:M) = kT/omegaP/omegaP

    thermo_lnv%Q_inv(1:M) = 1.0/thermo_lnv%Q(1)

    thermo_lnv%kT = kT

    thermo_lnv%Ndof_kT = kT

    thermo_lnv%eta(:) = 0.0d0

    do j = 1, M
        call gauss(0.d0, sqrt(kT/thermo_lnv%Q(j)), thermo_lnv%v(j))
    end do

    do j = 2, M
        thermo_lnv%a(j) = thermo_lnv%Q_inv(j)*(thermo_lnv%Q(j - 1)*thermo_lnv%v(j - 1)**2 - thermo_lnv%kT)
    end do

end subroutine part_setup_cnst_press_pimd

#endif /* LES */

subroutine full_setup_cnst_press_pimd(Nkt, tau_vol)
    use constants, only : kB, hbar
    use pimd_vars, only : nbead
    use nose_hoover_module, only : M, thermo_lnv, v_lnv, f_lnv_p, c2_lnv, mass_lnv, Thermostat_init
    implicit none
#include "md.h"
#include "memory.h"
    _REAL_ :: tau_vol, kT, Nkt, omegaP
    integer :: j

    kT = kB*temp0

    Nkt = natom*kT

    omegaP = kT*sqrt(dble(nbead))/hbar

    c2_lnv = 1 + 1.0/dble(natom)

    v_lnv = 0.d0

    f_lnv_p = 0.d0

    mass_lnv = 3*(natom + 1)*kT*tau_vol*tau_vol

    thermo_lnv%Q(1:M) = kT/omegaP/omegaP

    thermo_lnv%Q_inv(1:M) = 1.0/thermo_lnv%Q(1)

    thermo_lnv%kT = kT

    thermo_lnv%Ndof_kT = kT

    thermo_lnv%eta(:) = 0.0d0

    do j = 1, M
        call gauss(0.d0, sqrt(kT/thermo_lnv%Q(j)), thermo_lnv%v(j))
    end do

    do j = 2, M
        thermo_lnv%a(j) = thermo_lnv%Q_inv(j)*(thermo_lnv%Q(j - 1)*thermo_lnv%v(j - 1)**2 - thermo_lnv%kT)
    end do

end subroutine full_setup_cnst_press_pimd
