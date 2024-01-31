! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Emit coordinates or velocities, r(istart:n), to unit nf.
subroutine corpac(r, istart, n, nf, loutfm)
    use bintraj

    implicit none
    integer, intent(in) :: istart, n, nf
    _REAL_, intent(in) ::  r(n)
    logical, intent(in) ::  loutfm  ! true for formatted output

    integer i, j
    logical three_digit_fractional_part
    _REAL_ rmax, rmin
    parameter(rmax=9999.99d0)
    parameter(rmin=-999.99d0)

    if (istart > n) return

    if (.not. loutfm) then

        ! Unformatted writes:

        call write_binary_traj(r, istart, n, nf)

    else

        ! Formatted writes:

        three_digit_fractional_part = .true.
        do i = istart, n
            if (r(i) > rmax .or. r(i) < rmin) then
                three_digit_fractional_part = .false.
                exit
            end if
        end do

        if (three_digit_fractional_part) then
            write (nf, '(10f8.3)') (r(j), j=istart, n)
        else
            write (nf, '(10f8.2)') (r(j), j=istart, n)
        end if

    end if

    return
end subroutine corpac
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calculate various center of mass related terms.
subroutine ekcmr(nspm, nsp, tma, ekcmt, xr, v, amass, istart, iend)
    use constants, only : zero, one
    implicit none

    !     ----- ROUTINE TO CALCULATE THE TOTAL KINETIC ENERGY OF THE
    !           CENTER OF MASS OF THE SUB-MOLECULES AND ALSO THE
    !           COORDINATES OF THE MOLECULES RELATIVE TO THE CENTER OF
    !           MASS -----

    integer nspm
    integer nsp(*)
    _REAL_ tma(*)
    _REAL_ ekcmt(*)
    _REAL_ xr(*)
    _REAL_ v(*)
    _REAL_ amass(*)
    integer istart
    integer iend

#  include "box.h"

    _REAL_ aamass
    integer i3
    integer iat
    integer j
    integer j3
    integer n
    integer nn
    _REAL_ tmn
    _REAL_ vcm(3)
    _REAL_ xcm(3)

    i3 = 0
    iat = 0
    ekcmt(1) = ZERO
    ekcmt(2) = ZERO
    ekcmt(3) = ZERO

    do n = 1, nspm
        nn = nsp(n)
        tmn = ONE/tma(n)

        !       ----- CALCULATE THE CENTER OF MASS AND THEN MOVE EACH
        !             SUB-MOLECULE TO ITS CENTER OF MASS -----

        j3 = i3
        xcm(1) = ZERO
        xcm(2) = ZERO
        xcm(3) = ZERO
        vcm(1) = ZERO
        vcm(2) = ZERO
        vcm(3) = ZERO
        do j = 1, nn
            aamass = amass(iat + j)
            xcm(1) = xcm(1) + xr(j3 + 1)*aamass
            xcm(2) = xcm(2) + xr(j3 + 2)*aamass
            xcm(3) = xcm(3) + xr(j3 + 3)*aamass

            ! each processor knows all coordinates, but only some velocities:
            if (iat + j >= istart .and. iat + j <= iend) then
                vcm(1) = vcm(1) + v(j3 + 1)*aamass
                vcm(2) = vcm(2) + v(j3 + 2)*aamass
                vcm(3) = vcm(3) + v(j3 + 3)*aamass
            end if

            j3 = j3 + 3
        end do

        xcm(1) = xcm(1)*tmn
        xcm(2) = xcm(2)*tmn
        xcm(3) = xcm(3)*tmn

        j3 = i3
        do j = 1, nn
            xr(j3 + 1) = xr(j3 + 1) - xcm(1)
            xr(j3 + 2) = xr(j3 + 2) - xcm(2)
            xr(j3 + 3) = xr(j3 + 3) - xcm(3)
            j3 = j3 + 3
        end do

        ekcmt(1) = ekcmt(1) + tmn*vcm(1)*vcm(1)
        ekcmt(2) = ekcmt(2) + tmn*vcm(2)*vcm(2)
        ekcmt(3) = ekcmt(3) + tmn*vcm(3)*vcm(3)

        !       ----- END OF CALCULATION FOR EACH SUB-MOLECULE -----

        i3 = j3
        iat = iat + nn
    end do

    return
end subroutine ekcmr
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Open the coordinate, velocity, energy and constant pH output files.
subroutine open_dump_files
    use bintraj, only : open_binary_files
    use file_io_dat

    implicit none
#  include "md.h"
#  include "extra.h"

    !     subr amopen(lun,fname,fstat,fform,facc)

    if (master) then
        if (ioutfm <= 0) then

            !     ----- FORMATTED DUMPING -----

            if (ntwx > 0) then
!              call amopen(MDCRD_UNIT,mdcrd,owrite,'F',facc)
                call amopen(MDCRD_UNIT, mdcrd, 'U', 'F', facc)
                if (facc /= 'A') write (MDCRD_UNIT, 1000) title
            end if
            if (ntwv > 0) then
                call amopen(MDVEL_UNIT, mdvel, owrite, 'F', 'W')
                write (MDVEL_UNIT, 1000) title
            end if
        else
            call open_binary_files
        end if  ! (ioutfm <= 0)
        if (icnstph /= 0) then
            call amopen(CPOUT_UNIT, cpout, owrite, 'F', 'W')
        end if
        if (ntwe > 0) then
            call amopen(MDEN_UNIT, mden, owrite, 'F', 'W')
        end if

    end if  ! ( master )
1000 format(a80)
    return
end subroutine open_dump_files
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Close the coordinate, velocity, energy and constant pH output files.
subroutine close_dump_files
    use bintraj, only : close_binary_files
    use file_io_dat

    implicit none
#  include "extra.h"
#  include "md.h"

    if (master) then
        if (ioutfm > 0) then
            call close_binary_files
        else
            if (ntwx > 0) close (MDCRD_UNIT)
            if (ntwv > 0) close (MDVEL_UNIT)
        end if
        if (ntwe > 0) close (MDEN_UNIT)
        if (imin == 5) close (INPTRAJ_UNIT)
        if (ntpr > 0) close (7)
        if (icnstph /= 0) close (CPOUT_UNIT)

    end if
    return
end subroutine close_dump_files
!----------------------------------------------------------------------

#ifndef PBSA
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Energy output for md, in human-readable form.
subroutine prntmd(nstep, nitp, nits, time, ener, onefac, iout7, rms)
    use pimd_vars, only : ipimd, nbead, natomCL, itimass
    use neb_vars, only : ineb
    use cmd_vars, only : adiab_param
    use file_io_dat

#if defined( MPI )
    use evb_parm, only : nevb, nbias, evb_dyn
    use evb_data, only : evb_frc, evb_vel0, evb_Hmat, evb_bias &
        , evb_nrg, evb_nrg_ave, evb_nrg_rms
    ! DAN ROE: Added for printing replica information
    use remd, only : rem, repnum, mdloop, stagid
#ifdef LES
    use evb_pimd, only : evb_vec0_bead, nbead_inv
#endif
#endif /* MPI */
#ifdef RISMSANDER
    use sander_rism_interface, only : rismprm, RISM_NONE, RISM_FULL, RISM_INTERP, &
        rism_calc_type, rism_thermo_print
#endif

    use qmmm_module, only : qmmm_nml, qmmm_struct
    use cns_xref
#ifdef _XRAY
    use xray_interface_module, only : xray_active, xray_write_md_state
#endif
    use state
    use charmm_mod, only : charmm_active
    use crg_reloc, only : ifcr
! SGLD
    use sgld, only : tsgld, trxsgld, sgft, tempsg
    use ff11_mod, only : cmap_active
    use nbips, only : ips
    use emap, only : temap

    implicit none

! Passed in
    integer, intent(in) :: nstep, nitp, nits, iout7
    _REAL_, intent(in) :: onefac(3), time
    type(state_rec), intent(in) :: ener
    logical, intent(in) :: rms ! true if this is a print of RMS fluctuations

!Local
    integer :: n, total_nstlim, total_nstep
    _REAL_ :: enb14, edihed, eangle, ebond, ehbond, eel, enonb, epot, egb, epb
    _REAL_ :: virx, viry, virz, virt, dipole_temp, dipiter, diprms, virvsene
    _REAL_ :: dvdl, esurf, avetot, aveind, aveper, epol, econst, eel14, ekcmt
    _REAL_ :: edisp, tempsgj, sgftj, escf, ekcmx, ekcmy, ekcmz, enemap
    _REAL_ :: press, presx, presy, presz, densit, volume, boxx, boxy, boxz
    _REAL_ :: rms_pbs, eksolv, eksolt, ektot, etot
#ifdef RISMSANDER
    _REAL_ :: erism
    _REAL_ :: pot_array(potential_energy_rec_len)
#endif /*RISMSANDER*/
    _REAL_ :: ect

#  include "md.h"
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "tgtmd.h"
#ifdef LES
#  include "les.h"
#endif

    !     ----- DEFINE VARIOUS TERMS -----

    etot = ener%tot
    ektot = ener%kin%tot
    if (ipimd == 3) then
        temp = ener%kin%tot*onefac(1)  !  for CMD, this is the centroid KE
    else
        temp = ener%kin%tot*onefac(1)  !  for PIMD, this is mean classical KE of beads
    end if
    eksolt = ener%kin%solt*onefac(2)
#ifdef LES

    ! modified for LES, ener(4) is LES KE now, not solvent
    ! so it should be reported along with temperature for LES region

    if (temp0les < 0.d0) then
        if (ntt == 5) then
            eksolv = ener%kin%solv*onefac(3)
        else
            eksolv = 0.0d0
        end if
    else
        eksolv = ener%kin%solv*onefac(3)
    end if

    if (ipimd > 0) then
        ektot = ener%kin%solv
    end if
#else
    if (ntt == 5) then
        eksolv = ener%kin%solv*onefac(3)
    else
        rms_pbs = ener%kin%solv
    end if
#endif

    boxx = ener%box(1)
    boxy = ener%box(2)
    boxz = ener%box(3)
    volume = ener%volume
    densit = ener%density

    presx = ener%pres(1)
    presy = ener%pres(2)
    presz = ener%pres(3)
    press = ener%pres(4)

    ekcmx = ener%cmt(1)
    ekcmy = ener%cmt(2)
    ekcmz = ener%cmt(3)
    ekcmt = ener%cmt(4)

    virx = ener%vir(1)
    viry = ener%vir(2)
    virz = ener%vir(3)
    virt = ener%vir(4)

    epot = ener%pot%tot
    enonb = ener%pot%vdw
    eel = ener%pot%elec
    ehbond = ener%pot%hbond
    egb = ener%pot%gb
    epb = ener%pot%pb
    ebond = ener%pot%bond
    eangle = ener%pot%angle
    edihed = ener%pot%dihedral
    enb14 = ener%pot%vdw_14
    eel14 = ener%pot%elec_14
    econst = ener%pot%constraint
    epol = ener%pot%polar
    virvsene = ener%virvsene
    aveper = ener%aveper
    aveind = ener%aveind
    avetot = ener%avetot
    esurf = ener%pot%surf
    dvdl = ener%pot%dvdl
    diprms = ener%diprms
    dipiter = ener%dipiter
    dipole_temp = ener%dipole_temp
    escf = ener%pot%scf
    edisp = ener%pot%disp
    enemap = ener%pot%emap
#ifdef RISMSANDER
    erism = ener%pot%rism
#endif /*RISMSANDER*/
    ect = ener%pot%ct

    write (6, 9018) nstep, time, temp, press
    write (6, 9028) etot, ektot, epot
    write (6, 9038) ebond, eangle, edihed
! CHARMM SPECIFIC ENERGY TERMS
    if (charmm_active) write (6, 9039) ener%pot%angle_ub, &
        ener%pot%imp, &
        ener%pot%cmap

    write (6, 9048) enb14, eel14, enonb

#ifdef RISMSANDER
    if (igb == 0 .and. ipb == 0 .and. rismprm%irism == 0) then
#else
    if (igb == 0 .and. ipb == 0) then
#endif
        write (6, 9058) eel, ehbond, econst
    else if (igb == 10 .or. ipb /= 0) then
        write (6, 9060) eel, epb, econst
#ifdef APBS
    else if (igb == 6 .and. mdin_apbs) then
        write (6, 9060) eel, epb, econst
#endif /* APBS */
#ifdef RISMSANDER
    else if (rismprm%irism == 1) then
        write (6, 9061) eel, erism, econst
#endif

    else
        write (6, 9059) eel, egb, econst
    end if

    if (ifcr /= 0) write (6, 9099) ect

#ifdef MPI
    if (ievb /= 0) then
        write (6, '(A)')
        write (6, '(A)') ' EVB:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        if (evb_frc%evb_ave) then
            write (6, 888) evb_nrg_ave(1), evb_nrg_ave(2), evb_nrg_ave(3)
            evb_frc%evb_ave = .false.
        else if (evb_frc%evb_rms) then
            write (6, 888) evb_nrg_rms(1), evb_nrg_rms(2), evb_nrg_rms(3)
            evb_frc%evb_rms = .false.
        else
            if (trim(adjustl(evb_dyn)) == "evb_map") then
                write (6, 777) evb_nrg(1), evb_nrg(2), evb_nrg(3)
            else
                write (6, 888) evb_nrg(1), evb_nrg(2), evb_nrg(3)
            end if
#ifdef LES
            write (6, 999) 'C_0^2  = ', (sum(evb_vec0_bead(n, :)**2)*nbead_inv, n=1, nevb)
#else
            write (6, 999) 'C_0^2  = ', (evb_Hmat%evb_vec0(n)**2, n=1, nevb)
#endif
            if (nbias > 0) write (6, 999) 'EVB RC = ', (evb_bias%RC(n), n=1, nbias)
            if (nbias > 1) write (6, 999) 'Vumb_i = ', (evb_bias%nrg_bias(n), n=1, nbias)
        end if
    end if
#endif /* MPI */

    if (qmmm_nml%ifqnt) then
        !write the SCF energy
        if (qmmm_nml%qmtheory%PM3) then
            if (qmmm_nml%qmmm_int == 3) then
                write (6, 9090) escf ! PM3-MM*
            else if (qmmm_nml%qmmm_int == 4) then
                write (6, 9096) escf ! PM3/MMX2
            else
                write (6, 9080) escf
            end if
        else if (qmmm_nml%qmtheory%AM1) then
            write (6, 9081) escf
        else if (qmmm_nml%qmtheory%AM1D) then
            write (6, 9981) escf
        else if (qmmm_nml%qmtheory%MNDO) then
            write (6, 9082) escf
        else if (qmmm_nml%qmtheory%MNDOD) then
            write (6, 9982) escf
        else if (qmmm_nml%qmtheory%PDDGPM3) then
            write (6, 9083) escf
        else if (qmmm_nml%qmtheory%PDDGMNDO) then
            write (6, 9084) escf
        else if (qmmm_nml%qmtheory%PM3CARB1) then
            write (6, 9085) escf
        else if (qmmm_nml%qmtheory%DFTB) then
            write (6, 9086) escf
        else if (qmmm_nml%qmtheory%RM1) then
            write (6, 9087) escf
        else if (qmmm_nml%qmtheory%PDDGPM3_08) then
            write (6, 9088) escf
        else if (qmmm_nml%qmtheory%PM6) then
            write (6, 9089) escf
        else if (qmmm_nml%qmtheory%PM3ZNB) then
            write (6, 9091) escf
        else if (qmmm_nml%qmtheory%EXTERN) then
            write (6, 9092) escf
        else if (qmmm_nml%qmtheory%PM3MAIS) then
            write (6, 9093) escf
        else
            write (6, '(" ERROR - UNKNOWN QM THEORY")')
        end if
    end if
#ifdef PUPIL_SUPPORT
!jtc ========================= PUPIL INTERFACE =========================
    write (6, 9900) escf
!jtc ========================= PUPIL INTERFACE =========================
#endif /*  !PUPIL_SUPPORT */
    if (gbsa > 0) write (6, 9077) esurf
    if (igb == 10 .or. ipb /= 0) write (6, 9074) esurf, edisp
#ifdef APBS
    if (igb == 6 .and. mdin_apbs) write (6, 9097) esurf
#endif /* APBS */
    if (econst /= 0.0) write (6, 9076) epot - econst
    if ((icfe > 0) .OR. (itimass > 0)) write (6, 9100) dvdl
    if (ineb > 0) call pimd_neb_energy_report(6)

#ifndef LES
    if (rms .and. ntt == 0) write (6, 9075) rms_pbs
#endif
    if (volume /= 0.0) write (6, 9078) ekcmt, virt, volume
#ifdef LES

    ! if LES T bath is used, write TEMPERATURE for LES and non-LES

    if (temp0les >= 0.d0) then
        write (6, 9067) eksolt, eksolv
    end if
#endif

    if (csurften > 0) &
        write (6, 9072) ener%surface_ten

    if (cmap_active .and. epol /= 0.0) then
        write (6, 9068) epol, ener%pot%cmap
    else
        if (cmap_active) write (6, 9069) ener%pot%cmap
        if (epol /= 0.0) write (6, 9070) epol
    end if
    !     if (induced.gt.0) WRITE(6,9288) aveper,aveIND,avetot
    if (induced > 0 .and. indmeth < 3) write (6, 9190) diprms, dipiter
    if (induced > 0 .and. indmeth == 3) write (6, 9191) diprms, &
        dipole_temp
    if (volume /= 0.0) write (6, 9079) densit
#ifndef noVIRIAL
    if (igb == 0 .and. ipb == 0 .and. iyammp == 0 .and. induced == 0 &
        .and. use_pme == 1) write (6, 9188) virvsene
#endif

    if (itgtmd == 1) then
        write (6, '(a,f8.3)') "Current RMSD from reference: ", &
            rmsdvalue
        write (6, '(a,f8.3)') "Current target RMSD:         ", &
            tgtrmsd
    end if
    !  Printout SGLD guiding information
    if (tsgld) then
!  SGLD properties
        write (6, 1005) ener%sgld%sgft, ener%sgld%tempsg, ener%sgld%templf, ener%sgld%treflf, &
            ener%sgld%frclf, ener%sgld%epotlf, ener%sgld%sgwt
        write (6, 1006) ener%sgld%sgff, ener%sgld%sgscal, ener%sgld%temphf, ener%sgld%trefhf, &
            ener%sgld%frchf, ener%sgld%epothf, ener%sgld%virsg
    end if

    call write_cns_xref_md_energies(ener)

#ifdef _XRAY
    if (xray_active) call xray_write_md_state(6)
#endif

#ifdef MPI
! DAN ROE: Print current REMD info (replica#, temp0, excgh#) only for iout7>0
!  (Not for average/rms)
    if (rem > 0 .and. rem /= 4 .and. iout7 > 0) then
        if (trxsgld) then
            write (6, 9064) temp0, sgft, tempsg, stagid, repnum, mdloop
        else
            write (6, 9065) temp0, repnum, mdloop
        end if
    else if (rem == 4) then
        write (6, 9066) solvph, repnum, mdloop
    end if
#endif
!  wxw: EMAP energy
    if (temap) then
        write (6, 9062) enemap
    end if

    write (6, 8088)

    !     --- flush i/o buffer ---

    call amflsh(6)
    if (iout7 == 0) return

    !       ----- OUTPUT THE INFO FILE if requested -----

    write (7, 9018) nstep, time, temp, press
    write (7, 9028) etot, ektot, epot
    write (7, 9038) ebond, eangle, edihed
! CHARMM SPECIFIC ENERGY TERMS
    if (charmm_active) write (7, 9039) ener%pot%angle_ub, &
        ener%pot%imp, &
        ener%pot%cmap

    write (7, 9048) enb14, eel14, enonb
#ifdef RISMSANDER
    if (igb == 0 .and. ipb == 0 .and. rismprm%irism == 0) then
#else
    if (igb == 0 .and. ipb == 0) then
#endif
        write (7, 9058) eel, ehbond, econst
    else if (igb == 10 .or. ipb /= 0) then
        write (7, 9060) eel, epb, econst
#ifdef APBS
    else if (igb == 6 .and. mdin_apbs) then
        write (7, 9060) eel, epb, econst
#endif /* APBS */
#ifdef RISMSANDER
    else if (rismprm%irism == 1) then
        write (7, 9061) eel, erism, econst
#endif
    else
        write (7, 9059) eel, egb, econst
    end if

#ifdef MPI
    if (ievb /= 0) then
        write (7, '(A)')
        write (7, '(A)') ' EVB:~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        if (evb_frc%evb_ave) then
            write (7, 888) evb_nrg_ave(1), evb_nrg_ave(2), evb_nrg_ave(3)
            evb_frc%evb_ave = .false.
        else if (evb_frc%evb_rms) then
            write (7, 888) evb_nrg_rms(1), evb_nrg_rms(2), evb_nrg_rms(3)
            evb_frc%evb_rms = .false.
        else
            if (trim(adjustl(evb_dyn)) == "evb_map") then
                write (7, 777) evb_nrg(1), evb_nrg(2), evb_nrg(3)
            else
                write (7, 888) evb_nrg(1), evb_nrg(2), evb_nrg(3)
            end if
#ifdef LES
            write (7, 999) 'C_0^2  = ', (sum(evb_vec0_bead(n, :)**2)*nbead_inv, n=1, nevb)
#else
            write (7, 999) 'C_0^2  = ', (evb_Hmat%evb_vec0(n)**2, n=1, nevb)
#endif
            if (nbias > 0) write (7, 999) 'EVB RC = ', (evb_bias%RC(n), n=1, nbias)
            if (nbias > 1) write (7, 999) 'Vumb_i = ', (evb_bias%nrg_bias(n), n=1, nbias)
        end if
    end if

! DAN ROE: Print current REMD info (replica#, temp0, excgh#) only for iout7>0
!  (Not for average/rms)
    if (rem > 0 .and. iout7 > 0) then
        if (trxsgld) then
            write (7, 9064) temp0, sgft, tempsg, stagid, repnum, mdloop
        else
            write (7, 9065) temp0, repnum, mdloop
        end if
    end if
#endif /* MPI */
!  wxw: EMAP energy
    if (temap .and. iout7 > 0) then
        write (7, 9062) enemap
    end if

    if (qmmm_nml%ifqnt) then
        !write the SCF energy
        if (qmmm_nml%qmtheory%PM3) then
            if (qmmm_nml%qmmm_int == 3) then
                write (7, 9090) escf ! PM3-MM*
            else if (qmmm_nml%qmmm_int == 4) then
                write (7, 9096) escf ! PM3/MMX2
            else
                write (7, 9080) escf
            end if
        else if (qmmm_nml%qmtheory%AM1) then
            write (7, 9081) escf
        else if (qmmm_nml%qmtheory%AM1D) then
            write (7, 9981) escf
        else if (qmmm_nml%qmtheory%MNDO) then
            write (7, 9082) escf
        else if (qmmm_nml%qmtheory%MNDOD) then
            write (7, 9982) escf
        else if (qmmm_nml%qmtheory%PDDGPM3) then
            write (7, 9083) escf
        else if (qmmm_nml%qmtheory%PDDGMNDO) then
            write (7, 9084) escf
        else if (qmmm_nml%qmtheory%PM3CARB1) then
            write (7, 9085) escf
        else if (qmmm_nml%qmtheory%DFTB) then
            write (7, 9086) escf
        else if (qmmm_nml%qmtheory%RM1) then
            write (7, 9087) escf
        else if (qmmm_nml%qmtheory%PDDGPM3_08) then
            write (7, 9088) escf
        else if (qmmm_nml%qmtheory%PM6) then
            write (7, 9089) escf
        else if (qmmm_nml%qmtheory%PM3ZNB) then
            write (7, 9091) escf
        else if (qmmm_nml%qmtheory%EXTERN) then
            write (7, 9092) escf
        else if (qmmm_nml%qmtheory%PM3MAIS) then
            write (7, 9093) escf
        else
            write (7, '(" ERROR - UNKNOWN QM THEORY")')
        end if
    end if
#ifdef PUPIL_SUPPORT
!jtc ========================= PUPIL INTERFACE =========================
    write (7, 9900) escf
!jtc ========================= PUPIL INTERFACE =========================
#endif /*  !PUPIL_SUPPORT */
    if (gbsa > 0) write (7, 9077) esurf
    if (igb == 10 .or. ipb /= 0) write (7, 9074) esurf, edisp
#ifdef APBS
    if (igb == 6 .and. mdin_apbs) write (7, 9097) esurf
#endif /* APBS */
    if (econst /= 0.0) write (7, 9076) epot - econst
    if ((icfe > 0) .OR. (itimass > 0)) write (7, 9100) dvdl
#ifndef LES
    if (rms .and. ntt == 0) write (6, 9075) rms_pbs
#endif
    if (volume /= 0.0) write (7, 9078) ekcmt, virt, volume
#ifdef LES

    ! if LES T bath is used, write TEMPERATURE for LES and non-LES

    if (temp0les >= 0.d0) then
        write (7, 9067) eksolt, eksolv
    end if
#endif

    if (csurften .gt. 0) &
        write (7, 9072) ener%surface_ten

    if (cmap_active .and. epol /= 0.0) then
        write (7, 9068) epol, ener%pot%cmap
    else
        if (cmap_active) write (7, 9069) ener%pot%cmap
        if (epol /= 0.0) write (7, 9070) epol
    end if
    !     if (induced.gt.0) WRITE(7,9088) aveper,aveIND,avetot
    if (induced > 0 .and. indmeth < 3) write (7, 9190) diprms, dipiter
    if (induced > 0 .and. indmeth == 3) write (7, 9191) diprms, &
        dipole_temp
    if (volume /= 0.0) write (7, 9079) densit
#ifndef noVIRIAL
    if (igb == 0 .and. ipb == 0 .and. iyammp == 0 .and. induced == 0 &
        .and. ips == 0 .and. use_pme == 1) write (7, 9188) virvsene
#endif
    !  Printout SGLD guiding information
    if (tsgld) then
!  SGLD properties
        write (7, 1005) ener%sgld%sgft, ener%sgld%tempsg, ener%sgld%templf, ener%sgld%treflf, &
            ener%sgld%frclf, ener%sgld%epotlf, ener%sgld%sgwt
        write (7, 1006) ener%sgld%sgff, ener%sgld%sgscal, ener%sgld%temphf, ener%sgld%trefhf, &
            ener%sgld%frchf, ener%sgld%epothf, ener%sgld%virsg
    end if
    call nmrptx(7)

#ifdef RISMSANDER
    if (rismprm%irism == 1 .and. rismprm%write_thermo == 1) then
        if (rism_calc_type(nstep) == RISM_FULL) &
            call rism_thermo_print(.false., transfer(ener%pot, pot_array))
    end if
#endif /*RISM*/

#ifndef NO_DETAILED_TIMINGS
    !Print Timing estimates to mdinfo.
#  ifdef MPI
    if (rem .ne. 0) then
        total_nstlim = nstlim*numexchg
        total_nstep = nstlim*(mdloop - 1) + nstep
    else
        total_nstlim = nstlim
        total_nstep = nstep
    end if
#  else
    total_nstlim = nstlim
    total_nstep = nstep
#  endif
    if (total_nstep /= 0) &
        call print_ongoing_time_summary(total_nstlim, nstep, dt, 7)
#endif

777 format(1x, 'V_TOT  = ', f14.4, 2x, 'V_EVB   = ', f14.4, 2x &
        , 'V_MAP      = ', f14.4)
888 format(1x, 'V_TOT  = ', f14.4, 2x, 'V_EVB   = ', f14.4, 2x &
        , 'V_UMB      = ', f14.4)
999 format(1x, A, (4(2x, f14.4)))
8088 format(t2, 78('-'),/)
9018 format(/1x, 'NSTEP =', i9, 3x, 'TIME(PS) =', f12.3, 2x, &
        'TEMP(K) =', f9.2, 2x, 'PRESS =', f8.1)
9028 format(1x, 'Etot   = ', f14.4, 2x, 'EKtot   = ', f14.4, 2x, &
        'EPtot      = ', f14.4)
9038 format(1x, 'BOND   = ', f14.4, 2x, 'ANGLE   = ', f14.4, 2x, &
        'DIHED      = ', f14.4)
9039 format(1x, 'UB     = ', f14.4, 2x, 'IMP     = ', f14.4, 2x, &
        'CMAP       = ', f14.4)
9048 format(1x, '1-4 NB = ', f14.4, 2x, '1-4 EEL = ', f14.4, 2x, &
        'VDWAALS    = ', f14.4)
9058 format(1x, 'EELEC  = ', f14.4, 2x, 'EHBOND  = ', f14.4, 2x, &
        'RESTRAINT  = ', f14.4)
9059 format(1x, 'EELEC  = ', f14.4, 2x, 'EGB     = ', f14.4, 2x, &
        'RESTRAINT  = ', f14.4)
9060 format(1x, 'EELEC  = ', f14.4, 2x, 'EPB     = ', f14.4, 2x, &
        'RESTRAINT  = ', f14.4)
#ifdef RISMSANDER
9061 format(1x, 'EELEC  = ', f14.4, 2x, 'ERISM   = ', f14.4, 2x, &
        'RESTRAINT  = ', f14.4)
#endif
! Xiongwu: add sgft, tempsg for RXSGLD
9064 format(1x, 'TEMP0= ', f6.1, 1x, 'SGFT= ', f4.2, 1x, 'TEMPSG= ', f6.1, 1x, &
        'STAGE= ', i3, 1x, 'REPNUM= ', i3, 1x, 'EXCHANGE=', i6)
9062 format(1x, 'EMAP   = ', f14.4)
! DAN ROE: Added Temp, Rep, Exchange
9065 format(1x, 'TEMP0  = ', f14.4, 2x, 'REPNUM  = ', i14, 2x, &
        'EXCHANGE#  = ', i14)
9066 format(1x, 'SOLVPH = ', f14.4, 2x, 'REPNUM  = ', i14, 2x, &
        'EXCHANGE#  = ', i14)
9072 format(52x, 'SURFTEN    = ', f14.4)
9074 format(1x, 'ECAVITY= ', f14.4, 2x, 'EDISPER = ', f14.4)
9075 format('|E(PBS) = ', f14.4)
9076 format(1x, 'EAMBER (non-restraint)  = ', f14.4)
9077 format(1x, 'ESURF= ', f14.4)
9078 format(1x, 'EKCMT  = ', f14.4, 2x, 'VIRIAL  = ', f14.4, 2x, &
        'VOLUME     = ', f14.4)
9079 format(52x, 'Density    = ', f14.4)
9080 format(1x, 'PM3ESCF= ', f14.4)
9081 format(1x, 'AM1ESCF= ', f14.4)
9981 format(1x, 'AM1DESCF= ', f14.4)
9082 format(1x, 'MNDOESCF= ', f13.4)
9982 format(1x, 'MNDODESCF= ', f13.4)
9083 format(1x, 'PDDGPM3-ESCF= ', f12.4)
9084 format(1x, 'PDDGMNDO-ESCF= ', f11.4)
9085 format(1x, 'PM3CARB1-ESCF= ', f11.4)
9086 format(1x, 'DFTBESCF= ', f13.4)
9087 format(1x, 'RM1ESCF= ', f14.4)
9088 format(1x, 'PDDGPM3_08-ESCF= ', f12.4)
9089 format(1x, 'PM6ESCF= ', f14.4)
9090 format(1x, 'PM3MMXESCF= ', f14.4)
9091 format(1x, 'PM3ZNBESCF= ', f14.4)
9092 format(1x, 'EXTERNESCF= ', f14.4)
9093 format(1x, 'PM3MAISESCF= ', f14.4)
9096 format(1x, 'PM3MMX2ESCF= ', f14.4)
9097 format(1x, 'ENPOLAR= ', f14.4)
9099 format(1x, 'ECRG   = ', f14.4)

#ifdef LES
    ! LES and non-LES temperatures (no solvent/solute)
9067 format(1x, 'T_non-LES =', f10.4, 2x, 'T_LES     = ', f10.4)
#endif
9068 format(1x, 'EPOLZ  = ', f14.4, 2x, 'CMAP    = ', f14.4)
9069 format(1x, 'CMAP   = ', f14.4)
9070 format(1x, 'EPOLZ  = ', f14.4)

9288 format(1x, 'DIPOLE MOMENTS/RESIDUE:', / &
        , 1x, 'PERMEN = ', f8.3, 2x, 'INDUCED = ', f8.3, &
        2x, 'VECTOR SUM =', f8.3)
9190 format(1x, 'Dipole convergence: rms = ', e10.3, ' iters = ', f6.2)
9191 format(1x, 'Dipole convergence: rms = ', e10.3, &
        ' temperature = ', f6.2)
9100 format(1x, 'DV/DL  = ', f14.4)
9188 format(1x, 'Ewald error estimate: ', e12.4)
1005 format(" SGLF = ", F8.4, X, F8.2, X, F9.4, X, F9.4, X, F7.4, X, F14.4, X, F10.4)
1006 format(" SGHF = ", F8.4, X, F8.6, X, F9.4, X, F9.4, X, F7.4, X, F14.4, X, F10.4)
#ifdef PUPIL_SUPPORT
!jtc ========================= PUPIL INTERFACE =========================
9900 format(1x, 'PUPESCF= ', f14.4)
!jtc ========================= PUPIL INTERFACE =========================
#endif /* !PUPIL_SUPPORT  */
    return
end subroutine prntmd
#endif /*ifndef PBSA*/
!-----------------------------------------------------------------------

subroutine pimd_report(nstep, time, ounit, ener, onefac)
    use state
    use pimd_vars, only : nbead, ipimd, itimass
    use neb_vars, only : ineb
    use cmd_vars, only : adiab_param
    use qmmm_module, only : qmmm_nml, qmmm_struct
    use charmm_mod, only : charmm_active
    use file_io_dat

    implicit none

#  include "md.h"
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#  include "tgtmd.h"
#ifdef LES
#  include "les.h"
#endif
    integer ounit, nstep
    type(state_rec) :: ener
    _REAL_ time, onefac(*)
    _REAL_ etot, ektot, eksolt, eksolv, rms_pbs
    _REAL_ boxx, boxy, boxz, volume, densit, presx, presy, presz, press
    _REAL_ ekcmx, ekcmy, ekcmz, ekcmt, virx, viry, virz, virt
    _REAL_ epot, enonb, eel, ehbond, ebond, eangle, edihed, enb14, eel14, econst
    _REAL_ epol, aveper, aveind, avetot, esurf, dvdl, virvsene, diprms, dipiter, egb, epb
    _REAL_ dipole_temp, escf, sgftj, tempsgj, edisp

    etot = ener%tot
    ektot = ener%kin%solv ! for PIMD, this is "virial" expression
    if (ipimd == 3 .or. neb > 0) then
        temp = ener%kin%tot*onefac(1) ! this is for CMD
    else
        temp = ener%kin%tot*onefac(1)/nbead  !  for PIMD, this is mean classical KE of beads
    end if

    eksolt = ener%kin%solt*onefac(2)
    boxx = ener%box(1)
    boxy = ener%box(2)
    boxz = ener%box(3)
    volume = ener%volume
    densit = ener%density

    presx = ener%pres(1)
    presy = ener%pres(2)
    presz = ener%pres(3)
    press = ener%pres(4)

    ekcmx = ener%cmt(1)
    ekcmy = ener%cmt(2)
    ekcmz = ener%cmt(3)
    ekcmt = ener%cmt(4)

    virx = ener%vir(1)
    viry = ener%vir(2)
    virz = ener%vir(3)
    virt = ener%vir(4)

    epot = ener%pot%tot
    enonb = ener%pot%vdw
    eel = ener%pot%elec
    ehbond = ener%pot%hbond
    egb = ener%pot%gb
    epb = ener%pot%pb
    ebond = ener%pot%bond
    eangle = ener%pot%angle
    edihed = ener%pot%dihedral
    enb14 = ener%pot%vdw_14
    eel14 = ener%pot%elec_14
    econst = ener%pot%constraint
    epol = ener%pot%polar
    virvsene = ener%virvsene
    aveper = ener%aveper
    aveind = ener%aveind
    avetot = ener%avetot
    esurf = ener%pot%surf
    dvdl = ener%pot%dvdl
    diprms = ener%diprms
    dipiter = ener%dipiter
    dipole_temp = ener%dipole_temp
    escf = ener%pot%scf
    edisp = ener%pot%disp

    write (ounit, 8088)
    write (ounit, 9018) nstep, time, temp, press
    write (ounit, 9028) etot, ektot, epot
    write (ounit, 9038) ebond, eangle, edihed
! CHARMM SPECIFIC ENERGY TERMS
    if (charmm_active) write (6, 9039) ener%pot%angle_ub, &
        ener%pot%imp, &
        ener%pot%cmap

    write (ounit, 9048) enb14, eel14, enonb
    if (igb == 0 .and. ipb == 0) then
        write (ounit, 9058) eel, ehbond, econst
    else if (igb == 10 .or. ipb /= 0) then
        write (ounit, 9060) eel, epb, econst
    else
        write (ounit, 9059) eel, egb, econst
    end if
    ! For PIMD/NMPIMD/CMD/RPMD output
    if (volume /= 0.0) write (ounit, 9078) ekcmt, virt, volume

    if (volume /= 0.0) write (ounit, 9079) densit
    ! Report of dV/dl for TI w.r.t. mass
    if (itimass > 0) write (ounit, 9100) dvdl

    if (qmmm_nml%ifqnt) then
        !write the SCF energy
        if (qmmm_nml%qmtheory%PM3) then
            write (ounit, 9080) escf
        else if (qmmm_nml%qmtheory%AM1) then
            write (ounit, 9081) escf
        else if (qmmm_nml%qmtheory%AM1D) then
            write (ounit, 9981) escf
        else if (qmmm_nml%qmtheory%MNDO) then
            write (ounit, 9082) escf
        else if (qmmm_nml%qmtheory%MNDOD) then
            write (ounit, 9982) escf
        else if (qmmm_nml%qmtheory%PDDGPM3) then
            write (ounit, 9083) escf
        else if (qmmm_nml%qmtheory%PDDGMNDO) then
            write (ounit, 9084) escf
        else if (qmmm_nml%qmtheory%PM3CARB1) then
            write (ounit, 9085) escf
        else if (qmmm_nml%qmtheory%DFTB) then
            write (ounit, 9086) escf
        else if (qmmm_nml%qmtheory%RM1) then
            write (ounit, 9087) escf
        else if (qmmm_nml%qmtheory%PDDGPM3_08) then
            write (ounit, 9088) escf
        else if (qmmm_nml%qmtheory%PM6) then
            write (ounit, 9089) escf
        else if (qmmm_nml%qmtheory%PM3ZNB) then
            write (ounit, 9091) escf
        else if (qmmm_nml%qmtheory%EXTERN) then
            write (ounit, 9092) escf
        else if (qmmm_nml%qmtheory%PM3MAIS) then
            write (ounit, 9093) escf
        else
            write (7, '(" ERROR - UNKNOWN QM THEORY")')
        end if
    end if
    if (gbsa > 0) write (7, 9077) esurf
    if (igb == 10 .or. ipb /= 0) write (7, 9074) esurf, edisp
    if (econst /= 0.0) write (7, 9076) epot - econst
    if (ineb > 0) call pimd_neb_energy_report(ounit)

    ! For PIMD/NMPIMD/CMD/RPMD output
8088 format(t2, 78('-'),/)
9018 format(/1x, 'NSTEP =', i9, 3x, 'TIME(PS) =', f12.5, 2x, &
        'TEMP(K) =', f9.2, 2x, 'PRESS =', f8.1)
9028 format(1x, 'Etot   = ', f14.4, 2x, 'EKtot   = ', f14.4, 2x, &
        'EPtot      = ', f14.4)
9038 format(1x, 'BOND   = ', f14.4, 2x, 'ANGLE   = ', f14.4, 2x, &
        'DIHED      = ', f14.4)
9039 format(1x, 'UB     = ', f14.4, 2x, 'IMP     = ', f14.4, 2x, &
        'CMAP       = ', f14.4)

9048 format(1x, '1-4 NB = ', f14.4, 2x, '1-4 EEL = ', f14.4, 2x, &
        'VDWAALS    = ', f14.4)
9058 format(1x, 'EELEC  = ', f14.4, 2x, 'EHBOND  = ', f14.4, 2x, &
        'RESTRAINT  = ', f14.4)
9059 format(1x, 'EELEC  = ', f14.4, 2x, 'EGB     = ', f14.4, 2x, &
        'RESTRAINT  = ', f14.4)
9060 format(1x, 'EELEC  = ', f14.4, 2x, 'EPB     = ', f14.4, 2x, &
        'RESTRAINT  = ', f14.4)
9074 format(1x, 'ECAVITY= ', f14.4, 2x, 'EDISPER = ', f14.4)
9075 format('|E(PBS) = ', f14.4)
9076 format(1x, 'EAMBER (non-restraint)  = ', f14.4)
9077 format(1x, 'ESURF= ', f14.4)
9078 format(1x, 'EKCMT  = ', f14.4, 2x, 'VIRIAL  = ', f14.4, 2x, &
        'VOLUME     = ', f14.4)
9079 format(52x, 'Density    = ', f14.4)
9080 format(1x, 'PM3ESCF= ', f14.4)
9081 format(1x, 'AM1ESCF= ', f14.4)
9981 format(1x, 'AM1DESCF= ', f14.4)
9082 format(1x, 'MNDOESCF= ', f13.4)
9982 format(1x, 'MNDODESCF= ', f13.4)
9083 format(1x, 'PDDGPM3-ESCF= ', f12.4)
9084 format(1x, 'PDDGMNDO-ESCF= ', f11.4)
9085 format(1x, 'PM3CARB1-ESCF= ', f11.4)
9086 format(1x, 'DFTBESCF= ', f13.4)
9087 format(1x, 'RM1ESCF= ', f14.4)
9088 format(1x, 'PDDGPM3_08-ESCF= ', f12.4)
9089 format(1x, 'PM6ESCF= ', f14.4)
9091 format(1x, 'PM3ZNBESCF= ', f14.4)
9092 format(1x, 'EXTERNESCF= ', f14.4)
9093 format(1x, 'PM3MAISESCF= ', f14.4)
9100 format(1x, 'DV/DL  = ', f14.4)
end subroutine pimd_report

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Assign velocities from a Maxwellian distribution.
subroutine setvel(nr, v, winv, tempi, init, iscale, scalm)
    use amoeba_mdin, only : iamoeba, beeman_integrator
    implicit none
#  include "memory.h"

! Variable descriptions
!
! Passed variables
!  nr    : number of atoms
!  v     : velocity array
!  winv  : inverse mass array
!  tempi : temperature to randomize velocities to
!  init  : defines what stage of simulation we're in
!  iscale: nmr restraint-related variable
!  scalm : nmr restraint-related variable
!
! Local variables
!  i, j  : loop indices
!  nr3   : nr * 3
!  boltz : boltzmann's constant in appropriate units
!  sd    : maxwell-boltzmann factor from which to derive temperature

    ! Passed variables
    _REAL_, intent(out) :: v(*)
    _REAL_, intent(in)  :: winv(*), tempi, scalm
    integer, intent(in) :: nr, init, iscale

    ! Local variables
    integer             :: i, j, nr3
    _REAL_              :: boltz, sd

    nr3 = 3*nr

    !     ----- Assign velocities from a Maxwellian distribution

    if (tempi < 1.d-6) then
        v(1:nr3 + iscale) = 0.d0
        return
    end if

    boltz = 8.31441d-3*tempi/4.184d0
    i = 0
    do j = 1, nr
        sd = sqrt(boltz*winv(j))
        call gauss(0.d0, sd, v(i + 1))
        call gauss(0.d0, sd, v(i + 2))
        call gauss(0.d0, sd, v(i + 3))
        i = i + 3
    end do
    if (iscale > 0) then
        sd = sqrt(boltz/scalm)
        do j = 1, iscale
            call gauss(0.d0, sd, v(i + j))
        end do
    end if
    if (iamoeba == 1 .and. beeman_integrator == 1) then
        do j = 1, nr3
            v(j) = 20.455d0*v(j) ! beeman uses time in ps units
        end do
    end if
    return
end subroutine setvel

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mdeng here]
subroutine mdeng(nf, nstep, time, ener, onefac, ntp, csurften)
    use state
    use constants, only : zero
    implicit none
#  include "box.h"
    integer nf, nstep, ntp, i
    integer, intent(in) :: csurften
    _REAL_ onefac(3), time
    type(state_rec), intent(in) :: ener
    logical first
    character(len=16) labs(42)
    save first
    save labs
    data first/.true./
    data labs/'Nsteps  ', 'time(ps)  ', 'Etot  ', 'EKinetic  ', &
        'Temp  ', 'T_solute ', 'T_solv  ', 'Pres_scal_solu ', &
        'Pres_scal_solv ', 'BoxX  ', 'BoxY  ', 'BoxZ  ', &
        'volume  ', 'pres_X  ', 'pres_Y  ', 'pres_Z  ', &
        'Pressure ', 'EKCoM_x ', 'EKCoM_y ', 'EKCoM_z', &
        'EKComTot ', 'VIRIAL_x ', 'VIRIAL_y ', 'VIRIAL_z ', &
        'VIRIAL_tot ', 'E_pot  ', 'E_vdw  ', 'E_el  ', &
        'E_hbon  ', 'E_bon  ', 'E_angle  ', 'E_dih  ', &
        'E_14vdw  ', 'E_14el  ', 'E_const  ', 'E_pol  ', &
        'AV_permMoment ', 'AV_indMoment ', 'AV_totMoment ', &
        'Density', 'dV/dlambda', 'surften'/

    !     ----- DEFINE VARIOUS TERMS -----

    if (first) then
        !       -- up to Ekinetic:
        write (nf, 1) 'L0 ', (labs(i), i=1, 4)
        !       -- up to Pres_scal_solu:
        write (nf, 1) 'L1 ', (labs(i), i=5, 8)
        !       -- up to boxZ:
        write (nf, 1) 'L2 ', (labs(i), i=9, 12)
        !       -- up to pres_Z:
        write (nf, 1) 'L3 ', (labs(i), i=13, 16)
        !       -- up to EKCoM_z:
        write (nf, 1) 'L4 ', (labs(i), i=17, 20)
        !       -- up to VIRIAL_z:
        write (nf, 1) 'L5 ', (labs(i), i=21, 24)
        !       -- up to E_el:
        write (nf, 1) 'L6 ', (labs(i), i=25, 28)
        !       -- up to E_dih:
        write (nf, 1) 'L7 ', (labs(i), i=29, 32)
        !       -- up to E_pol:
        write (nf, 1) 'L8 ', (labs(i), i=33, 36)
        !       -- up to Density:
        write (nf, 1) 'L9 ', (labs(i), i=37, 41)
        ! surface tension info if constant surface tension in use.
        if (csurften > 0) &
            write (nf, 1) 'L10 ', labs(42)
1       format(a, 10(1x, a))
        first = .false.
    end if

    !     ----- write values for this step -----

    !     -- up to Ekinetic:
    write (nf, 2) 'L0 ', nstep, time, ener%tot, ener%kin%tot
    !     -- up to Pres_scal_solu:
    write (nf, 3) 'L1 ', ener%kin%tot*onefac(1), ener%kin%solt*onefac(2), &
        ener%kin%solv*onefac(3), ener%kin%pres_scale_solt
    !     -- up to boxZ:
    write (nf, 3) 'L2 ', ener%kin%pres_scale_solv, ener%box(1), ener%box(2), ener%box(3)
    !     -- up to pres_Z:
    write (nf, 3) 'L3 ', ener%volume, ener%pres(1), ener%pres(2), ener%pres(3)
    !     -- up to EKCoM_z:
    write (nf, 3) 'L4 ', ener%pres(4), ener%cmt(1), ener%cmt(2), ener%cmt(3)
    !     -- up to VIRIAL_z:
    if (ntp > 0) then
        write (nf, 3) 'L5 ', ener%cmt(4), ener%vir(1), ener%vir(2), ener%vir(3)
    else
        write (nf, 3) 'L5 ', ener%cmt(4), ZERO, ZERO, ZERO
    end if
    !     -- up to E_el:
    write (nf, 3) 'L6 ', ener%vir(4), ener%pot%tot, ener%pot%vdw, ener%pot%elec
    !     -- up to E_dih:
    write (nf, 3) 'L7 ', ener%pot%hbond, ener%pot%bond, ener%pot%angle, ener%pot%dihedral
    !      -- up to E_pol:
    write (nf, 3) 'L8 ', ener%pot%vdw_14, ener%pot%elec_14, ener%pot%constraint, &
        ener%pot%polar
    !     -- up to dV/dlambda:
    write (nf, 3) 'L9 ', ener%aveper, ener%aveind, ener%avetot, &
        ener%density, ener%pot%dvdl

    ! Constant surface tension info if running with constant surface tension
    if (csurften > 0) &
        write (nf, 3) 'L10 ', ener%surface_ten

2   format(a, i8, 20(2x, e16.10))
3   format(a, 20(e16.10, 2x))
    return
end subroutine mdeng

!+ [Enter a one-line description of subroutine mden2 here]
subroutine mden2(nf, nstep, time, ener, onefac)

    use state

    implicit none
    integer                     :: nf, nstep
    _REAL_                      :: onefac, time
    dimension                   :: onefac(3)
    type(state_rec), intent(in) :: ener

    !     ----- DEFINE VARIOUS TERMS -----

    write (nf, '(''Nsteps = '',i10)') nstep
    write (nf, '(''Time               = '',f20.10,'' p.s. '')') time
    write (nf, '(''Etotal             = '',f20.10)') ener%tot
    write (nf, '(''EKinetic           = '',f20.10)') ener%kin%tot
    write (nf, '(''Temperature        = '',f20.10)') ener%kin%tot*onefac(1)
    write (nf, '(''Temp_solute        = '',f20.10)') ener%kin%solt*onefac(2)
    write (nf, '(''Temp_solvent       = '',f20.10)') ener%kin%solv*onefac(3)
    write (nf, '(''Press_SCALE_solute = '',f20.10)') ener%kin%pres_scale_solt
    write (nf, '(''Press_SCALE_solvent= '',f20.10)') ener%kin%pres_scale_solv
    write (nf, '(''BOX_x              = '',f20.10)') ener%box(1)
    write (nf, '(''BOX_y              = '',f20.10)') ener%box(2)
    write (nf, '(''BOX_z              = '',f20.10)') ener%box(3)
    write (nf, '(''VOLUME             = '',f20.10)') ener%volume
    write (nf, '(''PRES_x             = '',f20.10)') ener%pres(1)
    write (nf, '(''PRES_y             = '',f20.10)') ener%pres(2)
    write (nf, '(''PRES_z             = '',f20.10)') ener%pres(3)
    write (nf, '(''PRESSURE           = '',f20.10)') ener%pres(4)
    write (nf, '(''EKCoM_x            = '',f20.10)') ener%cmt(1)
    write (nf, '(''EKCoM_y            = '',f20.10)') ener%cmt(2)
    write (nf, '(''EKCoM_z            = '',f20.10)') ener%cmt(3)
    write (nf, '(''EKCoMTotal         = '',f20.10)') ener%cmt(4)
    write (nf, '(''VIRIAL_x           = '',f20.10)') ener%vir(1)
    write (nf, '(''VIRIAL_y           = '',f20.10)') ener%vir(2)
    write (nf, '(''VIRIAL_z           = '',f20.10)') ener%vir(3)
    write (nf, '(''VIRIAL_Total       = '',f20.10)') ener%vir(4)
    write (nf, '(''Epotential         = '',f20.10)') ener%pot%tot
    write (nf, '(''vanderwaals        = '',f20.10)') ener%pot%vdw
    write (nf, '(''electrostatic      = '',f20.10)') ener%pot%elec
    write (nf, '(''h-bond             = '',f20.10)') ener%pot%hbond
    write (nf, '(''bond               = '',f20.10)') ener%pot%bond
    write (nf, '(''angle              = '',f20.10)') ener%pot%angle
    write (nf, '(''dihedral           = '',f20.10)') ener%pot%dihedral
    write (nf, '(''1-4v.d.w.          = '',f20.10)') ener%pot%vdw_14
    write (nf, '(''1-4electrostatic   = '',f20.10)') ener%pot%elec_14
    write (nf, '(''restraint          = '',f20.10)') ener%pot%constraint
    write (nf, '(''polarization       = '',f20.10)') ener%pot%polar
    write (nf, '(''ave perm moment    = '',f20.10)') ener%aveper
    write (nf, '(''ave ind  moment    = '',f20.10)') ener%aveind
    write (nf, '(''ave total moment   = '',f20.10)') ener%avetot
    write (nf, '(''density            = '',f20.10)') ener%density
    write (nf, '(''SGFT               = '',f20.10)') ener%kin%sgft
    write (nf, '(''TEMPSGI            = '',f20.10)') ener%kin%tempsg
    write (nf, '(''-----------------------------------------'')')

    return
end subroutine mden2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cenmas here]
subroutine cenmas(natom, x, v, amass, ekcm, xcm, vcm, acm, ekrot, ocm, icm)
    implicit none
    integer:: i, i0, icm, j, lh, m, mh, n, natom
    _REAL_ :: aamass, acm, amass, comvel, crit, d, ekcm, ekrot, ocm, &
        tcm, tmass, tmassinv, v, vcm, x, x1, x2, x3, xcm, &
        xx, xy, xz, yy, yz, zz
#  include "extra.h"

    !     ----- ROUTINE TO CALCULATE THE TRANSLATIONAL AND ROTATIONAL
    !           KINETIC ENERGIES AND VELOCITIES -----

    !     icm=0   return, do nothing
    !     icm=1   just return com coords in XCM
    !     icm=2   also calculate energy of com motion
    !     icm=3   also calcualte the angular momemtum
    !     icm=4   also calculate the inertia tensor
    !     icm=5   same as icm=4, but does not write anything to file - used by neb

    dimension x(*), v(*), amass(*), xcm(*), vcm(*), acm(*), ocm(*)
    dimension tcm(3, 3), lh(3), mh(3)

    data crit/1.0d-06/

    if (icm == 0) return

    i0 = 3*natom

    !     ----- CALCULATE THE CENTER OF MASS COORDINATES -----

    xcm(1) = 0.d0
    xcm(2) = 0.d0
    xcm(3) = 0.d0

    i = 0
    tmass = 0.d0
    do j = 1, natom
        aamass = amass(j)
        tmass = tmass + aamass
        xcm(1) = xcm(1) + x(i + 1)*aamass
        xcm(2) = xcm(2) + x(i + 2)*aamass
        xcm(3) = xcm(3) + x(i + 3)*aamass
        i = i + 3
    end do
    tmassinv = 1.d0/tmass
    xcm(1) = xcm(1)*tmassinv
    xcm(2) = xcm(2)*tmassinv
    xcm(3) = xcm(3)*tmassinv

    if (iabs(icm) < 2) return

    !     ----- CALCULATE THE VELOCITY AND THE TRANSLATIONAL
    !           KINETIC ENERGY OF THE CENTRE OF MASS -----

    ekcm = 0.d0
    vcm(1) = 0.0d0
    vcm(2) = 0.0d0
    vcm(3) = 0.0d0

    i = 0
    do j = 1, natom
        aamass = amass(j)
        vcm(1) = vcm(1) + v(i + 1)*aamass
        vcm(2) = vcm(2) + v(i + 2)*aamass
        vcm(3) = vcm(3) + v(i + 3)*aamass
        i = i + 3
    end do
    do i = 1, 3
        vcm(i) = vcm(i)*tmassinv
        ekcm = ekcm + vcm(i)*vcm(i)
    end do
    ekcm = ekcm*tmass*0.5d0
    comvel = sqrt(vcm(1)*vcm(1) + vcm(2)*vcm(2) + vcm(3)*vcm(3))

    if (iabs(icm) < 3) return

    !     ----- CALCULATE THE ANGULAR MOMENTUM ABOUT THE
    !           CENTER OF MASS ----

    acm(1) = 0.0d0
    acm(2) = 0.0d0
    acm(3) = 0.0d0

    i = 0
    do j = 1, natom
        aamass = amass(j)
        acm(1) = acm(1) + (x(i + 2)*v(i + 3) - x(i + 3)*v(i + 2))*aamass
        acm(2) = acm(2) + (x(i + 3)*v(i + 1) - x(i + 1)*v(i + 3))*aamass
        acm(3) = acm(3) + (x(i + 1)*v(i + 2) - x(i + 2)*v(i + 1))*aamass
        i = i + 3
    end do

    acm(1) = acm(1) - (xcm(2)*vcm(3) - xcm(3)*vcm(2))*tmass
    acm(2) = acm(2) - (xcm(3)*vcm(1) - xcm(1)*vcm(3))*tmass
    acm(3) = acm(3) - (xcm(1)*vcm(2) - xcm(2)*vcm(1))*tmass

    if (iabs(icm) < 4) return

    !     ----- CALCULATE THE INERTIA TENSOR -----

    xx = 0.d0
    xy = 0.d0
    xz = 0.d0
    yy = 0.d0
    yz = 0.d0
    zz = 0.d0

    i = 0
    do j = 1, natom
        x1 = x(i + 1) - xcm(1)
        x2 = x(i + 2) - xcm(2)
        x3 = x(i + 3) - xcm(3)
        aamass = amass(j)
        xx = xx + x1*x1*aamass
        xy = xy + x1*x2*aamass
        xz = xz + x1*x3*aamass
        yy = yy + x2*x2*aamass
        yz = yz + x2*x3*aamass
        zz = zz + x3*x3*aamass
        i = i + 3
    end do
    tcm(1, 1) = yy + zz
    tcm(2, 1) = -xy
    tcm(3, 1) = -xz
    tcm(1, 2) = -xy
    tcm(2, 2) = xx + zz
    tcm(3, 2) = -yz
    tcm(1, 3) = -xz
    tcm(2, 3) = -yz
    tcm(3, 3) = xx + yy

    !     ----- INVERT THE INERTIA TENSOR -----

    call matinv(tcm, 3, d, lh, mh)
    if (abs(d) <= crit) then
        write (6, 307)
307     format(/5x, '%CENMAS-F-INERTIA_TENSOR,  determinant', &
            ' is zero ... stop')
        call mexit(6, 1)
    end if

    !     ----- CALCULATE THE ANGULAR VELOCITY ABOUT THE CENTER OF
    !           MASS AND THE ROTATIONAL KINETIC ENERGY -----

    ekrot = 0.d0
    do m = 1, 3
        ocm(m) = 0.d0
        do n = 1, 3
            ocm(m) = ocm(m) + tcm(m, n)*acm(n)
        end do
        ekrot = ekrot + ocm(m)*acm(m)
    end do
    ekrot = ekrot*0.5d0

    if (icm < 0) return
    if (master .AND. icm .NE. 5) then
        write (6, '(/3x,a,f11.4,3x,a,f11.4,3x,a,f12.6)') 'KE Trans =', &
            ekcm, 'KE Rot =', ekrot, 'C.O.M. Vel =', comvel
        call amflsh(6)
    end if
    return
end subroutine cenmas

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Remove Center of Mass transrotational motion.
subroutine stopcm(nr, x, v, xcm, vcm, ocm, printmsg)

    implicit none
    integer nr
    _REAL_ x(*), v(*), xcm(*), vcm(*), ocm(*)
    logical printmsg !if true then a message is printed stating translational
    !and rotational motion removed. Else no message is printed
    !Useful for neb to stop too many messages being printed to screen.

#ifdef MPI
#  include "parallel.h"
#endif
#  include "extra.h"

    integer i, j, m
    _REAL_ x1, x2, x3

    !     ----- STOP THE CENTER OF MASS TRANSLATION -----

    i = 0
    do j = 1, nr
        do m = 1, 3
            i = i + 1
            v(i) = v(i) - vcm(m)
        end do
    end do

    !     ----- STOP THE ROTATION ABOUT THE CENTER OF MASS -----

    i = 0
    do j = 1, nr
        x1 = x(i + 1) - xcm(1)
        x2 = x(i + 2) - xcm(2)
        x3 = x(i + 3) - xcm(3)
        v(i + 1) = v(i + 1) - ocm(2)*x3 + ocm(3)*x2
        v(i + 2) = v(i + 2) - ocm(3)*x1 + ocm(1)*x3
        v(i + 3) = v(i + 3) - ocm(1)*x2 + ocm(2)*x1
        i = i + 3
    end do
    if (master .AND. printmsg) write (6, 9008)
9008 format(/3x, 'Translational and rotational motion removed')
    return
end subroutine stopcm
