! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++
!This module contains code for the linear
!response theory of electrostatic interactions
!+++++++++++++++++++++++++++++++++++++++++++++

module linear_response

    use findmask

    implicit none

    integer, save :: ilrt, calls, evals, lrt_interval, ier_alloc, numlrt
    integer, dimension(:), allocatable, save :: nlrt
    logical, save :: do_lrt
    character(len=256), save :: lrtmask
    _REAL_, save :: energy_w0, energy_m0, energy_vdw0, elrt, elrt2, sasa
    _REAL_, save :: elrtvdw, elrtvdw2, lrtsasa, lrtsasa2
    _REAL_, dimension(:), allocatable, save :: crg_w0, crg_m0, f_scratch
    _REAL_, dimension(:), allocatable, save :: cn1_lrt, cn2_lrt

contains

    !-------------------------------------------------------------------------------

    subroutine setup_linear_response(natom, nres, igraph, isymbl, ipres, lbres, crd, charge, ntypes, iac, ico, cn1, cn2, master)

        integer, intent(in) :: natom, nres, ipres(*), ntypes, iac(*), ico(ntypes, ntypes)
        character(len=4), intent(in) :: igraph(*), isymbl(*), lbres(*)
        _REAL_, intent(in) :: crd(*), charge(*), cn1(*), cn2(*)
        logical, intent(in) :: master

        integer i, j, index
        _REAL_ solvent_type(ntypes), solute_type(ntypes)

        do_lrt = .false.
        if (lrt_interval == 1) then
            do_lrt = .true.
        end if
        calls = 0
        evals = 0
        elrt = 0.0d0
        elrt2 = 0.0d0
        elrtvdw = 0.0d0
        elrtvdw2 = 0.0d0
        lrtsasa = 0.0d0
        lrtsasa2 = 0.0d0
        energy_vdw0 = 0

        allocate (nlrt(natom), stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('setup in linear_response.f', 'cant allocate nlrt', '')
        allocate (crg_m0(natom), stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('setup in linear_response.f', 'cant allocate crg_m0', '')
        allocate (crg_w0(natom), stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('setup in linear_response.f', 'cant allocate crg_w0', '')
        allocate (f_scratch(3*natom), stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('setup in linear_response.f', 'cant allocate f_scratch', '')

        allocate (cn1_lrt(ntypes*(ntypes + 1)/2), stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('setup in linear_response.f', 'cant allocate cn1_lrt', '')
        allocate (cn2_lrt(ntypes*(ntypes + 1)/2), stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('setup in linear_response.f', 'cant allocate cn2_lrt', '')

        call atommask(natom, nres, 0, igraph, isymbl, ipres, lbres, crd, lrtmask, nlrt)

        numlrt = sum(nlrt)

        if (master) then
            write (6, *)
            write (6, *) 'Setting up LRT calculation'
            write (6, *) '      '
            write (6, '(a,a,a,i5,a)') '     LRT Mask ', lrtmask(1:len_trim(lrtmask)), ' matches ', numlrt, ' atoms'
        end if

        if (master) then
            if (numlrt > 3500) then
                call sander_bomb('setup in linear_response.f', 'More than 3500 lrt atoms, exiting', '')
            end if
        end if

        do i = 1, natom
            if (nlrt(i) == 1) then
                crg_m0(i) = 0.0d0
                crg_w0(i) = charge(i)
            else
                crg_m0(i) = charge(i)
                crg_w0(i) = 0.0d0
            end if
        end do

        solute_type(1:ntypes) = 0
        solvent_type(1:ntypes) = 0

        do i = 1, natom
            if (nlrt(i) == 1) then
                solute_type(iac(i)) = 1
            else
                solvent_type(iac(i)) = 1
            end if
        end do

        do i = 1, ntypes
            if (solute_type(i) == solvent_type(i)) then
                call sander_bomb('setup in linear_response.f', 'an atom type occurs both in solvent and solute', '')
            end if
        end do

        cn1_lrt(1:ntypes*(ntypes + 1)/2) = 0.0d0
        cn2_lrt(1:ntypes*(ntypes + 1)/2) = 0.0d0

        do i = 1, ntypes
            do j = 1, ntypes
                index = ico(i, j)
                if (index > 0) then
                    if (solvent_type(i) == solvent_type(j)) then
                        cn1_lrt(index) = cn1(index)
                        cn2_lrt(index) = cn2(index)
                    else
                        cn1_lrt(index) = 0.0d0
                        cn2_lrt(index) = 0.0d0
                    end if
                else
                end if
            end do
        end do

        return

    end subroutine setup_linear_response

!-------------------------------------------------------------------------------

    subroutine ee_linear_response(energy, master)

        _REAL_ energy, elrt_tmp
        logical master

        if (do_lrt) then
            evals = evals + 1

            if (master) then
                elrt_tmp = 0.5d0*(energy - energy_m0 - energy_w0) ! yields EE interaction energy

                elrt = elrt + elrt_tmp
                elrt2 = elrt2 + elrt_tmp*elrt_tmp

                elrtvdw = elrtvdw + energy_vdw0
                elrtvdw2 = elrtvdw2 + (energy_vdw0*energy_vdw0)

                lrtsasa = lrtsasa + sasa
                lrtsasa2 = lrtsasa2 + (sasa*sasa)

                write (6, '(a,f12.4,a,f12.4,a,f12.4)') 'LIE: EE ', elrt_tmp, ' VDW ', energy_vdw0, ' SASA ', sasa
            end if
        end if

        calls = calls + 1

        if (mod(calls + 1, lrt_interval) == 0) then
            do_lrt = .true.
        else
            do_lrt = .false.
        end if

        return

    end subroutine ee_linear_response

!-------------------------------------------------------------------------------

    subroutine lrt_solute_sasa(x, natom, vdwrad)

        use nblist, only : a, b, c
        use constants, only : zero
        use icosasurf, only : icosa_init, icosa_sphere_approx

        _REAL_, intent(in) :: x(3, *), vdwrad(*)
        integer, intent(in) :: natom

        _REAL_ crd(3, numlrt), offset(3), radii(numlrt)
        _REAL_ dx, dy, dz, dist2, vdwrad2
        integer i, j, k, count, ineighbor(13000000), ineighborpt

        ! Take the solute atom coordinates out of x and tranform into fractionals
        k = 1
        do i = 1, natom
            if (nlrt(i) == 1) then
                crd(1, k) = x(1, i)/a
                crd(2, k) = x(2, i)/b
                crd(3, k) = x(3, i)/c
                radii(k) = vdwrad(i)
                k = k + 1
            end if
        end do

        ! image the solute coordinates around the first atom
        offset(1:3) = crd(1:3, 1)

        crd(1, 1:numlrt) = crd(1, 1:numlrt) - anint(crd(1, 1:natom) - offset(1))
        crd(2, 1:numlrt) = crd(2, 1:numlrt) - anint(crd(2, 1:natom) - offset(2))
        crd(3, 1:numlrt) = crd(3, 1:numlrt) - anint(crd(3, 1:natom) - offset(3))

        ! unfracture the coords

        crd(1, 1:numlrt) = a*crd(1, 1:numlrt)
        crd(2, 1:numlrt) = b*crd(2, 1:numlrt)
        crd(3, 1:numlrt) = c*crd(3, 1:numlrt)

        ! Prepare the icosasurf
        sasa = 0
        count = 0
        ! compute all intersolute distances and the SAS
        do i = 1, numlrt
            do j = 1, numlrt
                if (i /= j) then
                    dx = crd(1, i) - crd(1, j)
                    dy = crd(2, i) - crd(2, j)
                    dz = crd(3, i) - crd(3, j)
                    dist2 = dx*dx + dy*dy + dz*dz
                    vdwrad2 = (radii(i) + radii(j))**2
                    if (vdwrad2 > dist2) then
                        count = count + 1
                        ineighbor(count) = j
                    else
                    end if
                end if
            end do
            count = count + 1
            ineighbor(count) = 0
            if (i == 1) then
                ineighborpt = 1
                call icosa_init(2, 3, zero)
            end if
            sasa = sasa + icosa_sphere_approx(i, crd, radii, ineighborpt, ineighbor, 0)
        end do

        return

    end subroutine lrt_solute_sasa

!-------------------------------------------------------------------------------

    subroutine cleanup_linear_response(master)

        logical master

        if (master) then
            elrt = elrt/evals
            elrt2 = elrt2/evals - elrt*elrt

            elrtvdw = elrtvdw/evals
            elrtvdw2 = elrtvdw2/evals - elrtvdw*elrtvdw

            lrtsasa = lrtsasa/evals
            lrtsasa2 = lrtsasa2/evals - lrtsasa*lrtsasa

            ! Ignore numerical round off errors for single evaluations
            if (elrt2 < 0.0d0) elrt2 = 0.0d0
            if (elrtvdw2 < 0.0d0) elrtvdw2 = 0.0d0
            if (lrtsasa2 < 0.0d0) lrtsasa2 = 0.0d0

            elrt2 = sqrt(elrt2)
            elrtvdw2 = sqrt(elrtvdw2)
            lrtsasa2 = sqrt(lrtsasa2)

            write (6, '(a)') " Final LIE Evaluation:"
            write (6, '(a,i8,a,f8.2,a,f8.2,a,f8.2,a,f8.2,a,f8.2,a,f8.2)') &
                "Evals: ", evals, " E(LRT): ", elrt, " RMS ", elrt2, " E(Disp): ", &
                elrtvdw, " RMS ", elrtvdw2, " SASA: ", lrtsasa, " RMS ", lrtsasa2
        end if

        if (allocated(nlrt)) deallocate (nlrt, stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('cleanup in linear_response.f]', 'cant deallocate nlrt', '')
        if (allocated(crg_m0)) deallocate (crg_m0, stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('cleanup in linear_response.f]', 'cant deallocate crg_m0', '')
        if (allocated(crg_w0)) deallocate (crg_w0, stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('cleanup in linear_response.f]', 'cant deallocate crg_w0', '')
        if (allocated(f_scratch)) deallocate (f_scratch, stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('cleanup in linear_response.f]', 'cant deallocate f_scratch', '')
        if (allocated(cn1_lrt)) deallocate (cn1_lrt, stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('cleanup in linear_response.f]', 'cant deallocate cn1_lrt', '')
        if (allocated(cn2_lrt)) deallocate (cn2_lrt, stat=ier_alloc)
        if (ier_alloc /= 0) call sander_bomb('cleanup in linear_response.f]', 'cant deallocate cn2_lrt', '')

        return
    end subroutine cleanup_linear_response

end module linear_response
