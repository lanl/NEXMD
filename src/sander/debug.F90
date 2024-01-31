! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine debug_mem here]
subroutine debug_mem(natom, ntypes, &
    istartr, iendr, istarti, iendi)

    ! code for debugging the forces
    ! also for flow control (just for force debug so far)
    !------------------------------------------------------------------
    use parms, only : nttyp
    implicit none
#  include "debug.h"
    ! sets up some space for saved arrays
    integer natom, ntypes, istartr, iendr, istarti, iendi
    integer mem_ptr
    mem_ptr = istartr
    call adj_mem_ptr(mem_ptr, lscg, natom)
    call adj_mem_ptr(mem_ptr, lsrms, natom)
    call adj_mem_ptr(mem_ptr, lsdip, 3*natom)
    call adj_mem_ptr(mem_ptr, lsind, 3*natom)
    call adj_mem_ptr(mem_ptr, lscn1, nttyp)
    call adj_mem_ptr(mem_ptr, lscn2, nttyp)
    call adj_mem_ptr(mem_ptr, lf1, 3*natom)
    call adj_mem_ptr(mem_ptr, lf2, 3*natom)
    call adj_mem_ptr(mem_ptr, lf3, 3*natom)
    call adj_mem_ptr(mem_ptr, lf4, 3*natom)
    call adj_mem_ptr(mem_ptr, lf5, 3*natom)
    call adj_mem_ptr(mem_ptr, lf6, 3*natom)
    iendr = mem_ptr
    return
end subroutine debug_mem
!------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine load_debug here]
subroutine load_debug(nf)
    use nblist, only : nbflag, nbfilter, cutlist
    use charmm_mod, only : do_charmm_dump_gold
    implicit none
    integer nf
#  include "flocntrl.h"
#  include "debug.h"
#  include "md.h"

    integer ifind, j
    namelist /debugf/ do_dir, do_rec, do_adj, do_self, do_bond, &
        do_angle, do_ephi, doxconst, do_cap, do_14, &
        do_debugf, neglgdel, zerochg, zerovdw, zerodip, &
        atomn, nranatm, &
        ranseed, chkvir, dumpfrc, rmsfrc, do_tgt, &
        do_charmm_dump_gold, &
        do_pbdir, do_pbnp, do_pbfd

    ! default flow control all force routines turned on
    do_dir = 1
    do_rec = 1
    do_adj = 1
    do_self = 1
    do_bond = 1
    do_angle = 1
    do_ephi = 1
    doxconst = 1
    do_cap = 1
    do_14 = 1
    do_tgt = 1
    do_pbdir = 1
    do_pbnp = 1
    do_pbfd = 1

    ! debug off by default
    do_debugf = 0
    neglgdel = 5
    zerochg = 0
    zerovdw = 0
    zerodip = 0
    nranatm = 0
    ranseed = 71277
    chkvir = 0
    dumpfrc = 0
    rmsfrc = 0
    do j = 1, natomn
        atomn(j) = 0
    end do
    rewind (5)
    call nmlsrc('debugf', nf, ifind)
    if (ifind /= 0) then
        read (5, nml=debugf, err=190)
    end if

    if (do_debugf == 0) then
        do_dir = 1
        do_rec = 1
        do_adj = 1
        do_self = 1
        do_bond = 1
        do_angle = 1
        do_ephi = 1
        doxconst = 1
        do_tgt = 1
        do_cap = 1
        do_pbdir = 1
        do_pbnp = 1
        do_pbfd = 1

        neglgdel = 5
        zerochg = 0
        zerovdw = 0
        zerodip = 0
        nranatm = 0
        ranseed = 71277
        chkvir = 0
        dumpfrc = 0
        rmsfrc = 0
    end if
    if (do_debugf == 1) then
        ! disable list calls, use old list logic
        ! this prevents strange scaling behavior when
        nbflag = 0
        ntnb = 0
        nbfilter = 1.5d0*cutlist
    end if
    return
190 write (6, *) 'Error in &debugf namelist'
    call mexit(6, 1)
end subroutine load_debug
!----------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine debug_frc here]
subroutine debug_frc(xx, ix, ih, ipairs, x, f, &
    cn1, cn2, qsetup)
    use nblist, only : volume
    use parms, only : nttyp
    use state
    implicit none
#  include "flocntrl.h"
#  include "debug.h"
#  include "memory.h"
#  include "extra.h"
#  include "md.h"
#  include "ew_frc.h"
#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "ew_mpole.h"
    logical qsetup
    _REAL_ xx(*), x(*), f(3, *)
    _REAL_ cn1(*), cn2(*)
    integer ipairs(*), ix(*)
    character(len=4) ih(*)
    _REAL_ vir(4), time, onefac(3), rms
    type(state_rec) :: ener
    integer iout7, nstep, nitp, nits, j, k, atomnum
    _REAL_ apfrc(3), del, y
    _REAL_ apavir, apmvir, exavir, exmvir, exvir
    _REAL_ cpfrc1(3, natomn), &
        cpfrc2(3, natomn)
    _REAL_ dudv
    type(state_rec) :: ene
    _REAL_ intvir(3, 3), duda(3, 3), apduda(3, 3)
    _REAL_ mduda(3, 3), mapduda(3, 3)
    _REAL_ mvten(3, 3), avten(3, 3), diff
    integer ranatm(natomn), type

    ene = null_state_rec
    ener = null_state_rec

    if (do_debugf == 0) return
    ! save charge and cn1,cn2 arrays for later use
    call array_copy(xx(l15), xx(lscg), natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(xx(linddip), xx(lsind), 3*natom)
    end if
    call array_copy(cn1, xx(lscn1), nttyp)
    call array_copy(cn2, xx(lscn2), nttyp)
    if (zerochg == 1) then
        call zero_array(xx(l15), natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. zerodip == 1) then
        call zero_array(xx(linddip), 3*natom)
    end if
    if (zerovdw == 1) then
        call zero_array(cn1, nttyp)
        call zero_array(cn2, nttyp)
    end if

    ! first call force for analytic result
    if (master) then
        write (6, *) 'DEBUG FORCE!; calling force routine'
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        vir, ener, qsetup)
    call get_intvir(natom, x, f, subvir, intvir)
    do k = 1, 3
        do j = 1, 3
            atvir(j, k) = atvir(j, k) + intvir(j, k)
            ! now save molvir and atvir for virial checking, otherwise repeated
            ! calls to force overwrites correct values
            avten(j, k) = atvir(j, k)
            mvten(j, k) = molvir(j, k)
        end do
    end do
    exmvir = molvir(1, 1) + molvir(2, 2) + molvir(3, 3)
    exavir = atvir(1, 1) + atvir(2, 2) + atvir(3, 3)
    if (master) then
        write (6, *) 'DEBUG FORCE!; back from force routine'
    end if

    ene = ener ! ener has been populated by get_analfrc
    ! now copy it to ene

    iout7 = 0
    onefac(1) = 1.d0
    onefac(2) = 1.d0
    onefac(3) = 1.d0
    nstep = 0
    nitp = 0
    nits = 0
    ene%kin%tot = 0.d0
    ene%kin%solt = 0.d0
    ene%tot = ener%pot%tot
    ene%density = tmass/(0.602204d0*volume)
    ene%volume = volume
    ene%vir(4) = vir(1) + vir(2) + vir(3)
    if (master) then
        call prntmd(nstep, nitp, nits, t, ene, onefac, iout7, .false.)
    end if
    if (dumpfrc == 1) then
        call force_dump(natom, xx, ix, ih, ipairs, x, f, &
            ener, xx(lscg), xx(lscn1), xx(lscn2), xx(lsdip), xx(l15), &
            cn1, cn2, xx(lfixdip), xx(linddip), xx(lsind), qsetup)
    end if
    if (rmsfrc == 1) then
        call rms_check(natom, xx, ix, ih, ipairs, x, f, &
            ener, xx(lsrms), xx(lf2), xx(lf3), xx(lf4), xx(lf5), xx(lf6), &
            xx(lscg), xx(lscn1), xx(lscn2), xx(lsdip), &
            xx(l15), cn1, cn2, xx(lfixdip), &
            xx(linddip), xx(lsind), qsetup)
    end if

    ! get the delta for numerical force, virial calcs
    del = 1.d0
    do j = 1, neglgdel
        del = del/1.d1
    end do
    ! save copies of analytic forces for user defined and random atoms
    do j = 1, natomn
        atomnum = atomn(j)
        if (atomnum == 0) goto 50
        cpfrc1(1, j) = f(1, atomnum)
        cpfrc1(2, j) = f(2, atomnum)
        cpfrc1(3, j) = f(3, atomnum)
    end do
50  continue
    if (nranatm > natomn) then
        if (master) write (6, 166) natomn
166     format(1x, 'MAX NUM Random atoms: ', i6)
        call mexit(6, 1)
    end if
    if (nranatm > 0) then
        call amrset(ranseed)
        do j = 1, nranatm
            call amrand(y)
            atomnum = y*natom
            ranatm(j) = atomnum
            cpfrc2(1, j) = f(1, atomnum)
            cpfrc2(2, j) = f(2, atomnum)
            cpfrc2(3, j) = f(3, atomnum)
        end do
    end if
    ! now do user defined atoms
    if (atomn(1) > 0) then
        if (master) write (6, 168)
        if (master) &
            write (6, *) '----------------------------------------------'
    end if
    do j = 1, natomn
        atomnum = atomn(j)
        if (atomnum == 0) goto 100
        call get_numfrc(xx, ix, ih, ipairs, x, f, &
            vir, ener, atomnum, del, apfrc, qsetup)
        call rmsdiff(apfrc, cpfrc1(1, j), rms)
        if (master) then
            write (6, 60) atomnum
            do k = 1, 3
                diff = apfrc(k) - cpfrc1(k, j)
                write (6, 66) k, apfrc(k), cpfrc1(k, j), diff
            end do
            write (6, 70) rms
        end if
    end do
100 continue
    if (master) &
        write (6, *) '--------------------------------------------'
    if (nranatm > 0) then
        if (master) write (6, 167)
        if (master) &
            write (6, *) '--------------------------------------------'
        do j = 1, nranatm
            atomnum = ranatm(j)
            call get_numfrc(xx, ix, ih, ipairs, x, f, &
                vir, ener, atomnum, del, apfrc, qsetup)
            call rmsdiff(apfrc, cpfrc2(1, j), rms)
            if (master) then
                write (6, 60) atomnum
                do k = 1, 3
                    diff = apfrc(k) - cpfrc2(k, j)
                    write (6, 66) k, apfrc(k), cpfrc2(k, j), diff
                end do
                write (6, 70) rms
            end if
        end do
        if (master) &
            write (6, *) '--------------------------------------------'
    end if
    ! now check the virials. First the molvir, type = 1
    if (chkvir /= 0) then
        type = 1
        call check_virial(xx, ix, ih, ipairs, x, f, &
            ene, del, dudv, type, qsetup)
        apmvir = 3*volume*dudv
        call check_vtens(xx, ix, ih, ipairs, x, f, &
            ene, del, mduda, molvir, mapduda, type, qsetup)
        ! Next the atvir, type = 2
        type = 2
        call check_virial(xx, ix, ih, ipairs, x, f, &
            ene, del, dudv, type, qsetup)
        apavir = 3*volume*dudv
        call check_vtens(xx, ix, ih, ipairs, x, f, &
            ene, del, duda, atvir, apduda, type, qsetup)
        if (master) then
            write (6, *) 'Checking analytic virial trace versus'
            write (6, *) 'Numerical calculation of 3V dU/dV'
            write (6, *) '--------------------------------------------'
            call compare(exmvir, apmvir, 'Molecular virial:      ')
            call compare(exavir, apavir, 'Atomic virial:         ')
            write (6, *) '--------------------------------------------'
            write (6, *) 'Checking numerical calculation of DU/da_ij'
            write (6, *) 'where a is the unit cell matrix, against'
            write (6, *) 'VTa^-1, where T is molecular or atomic virial tensor'
            write (6, *) 'See eqn. 2.6-2.8 in Essmann et al. JCP 103,8577'
            write (6, *) '--------------------------------------------'
            call compare(mduda(1, 1), mapduda(1, 1), 'Molec.   dUda_(1,1)    ')
            call compare(mduda(1, 2), mapduda(1, 2), 'Molec.   dUda_(1,2)    ')
            call compare(mduda(1, 3), mapduda(1, 3), 'Molec.   dUda_(1,3)    ')
            call compare(mduda(2, 1), mapduda(2, 1), 'Molec.   dUda_(2,1)    ')
            call compare(mduda(2, 2), mapduda(2, 2), 'Molec.   dUda_(2,2)    ')
            call compare(mduda(2, 3), mapduda(2, 3), 'Molec.   dUda_(2,3)    ')
            call compare(mduda(3, 1), mapduda(3, 1), 'Molec.   dUda_(3,1)    ')
            call compare(mduda(3, 2), mapduda(3, 2), 'Molec.   dUda_(3,2)    ')
            call compare(mduda(3, 3), mapduda(3, 3), 'Molec.   dUda_(3,3)    ')
            call compare(duda(1, 1), apduda(1, 1), 'Atomic   dUda_(1,1)    ')
            call compare(duda(1, 2), apduda(1, 2), 'Atomic   dUda_(1,2)    ')
            call compare(duda(1, 3), apduda(1, 3), 'Atomic   dUda_(1,3)    ')
            call compare(duda(2, 1), apduda(2, 1), 'Atomic   dUda_(2,1)    ')
            call compare(duda(2, 2), apduda(2, 2), 'Atomic   dUda_(2,2)    ')
            call compare(duda(2, 3), apduda(2, 3), 'Atomic   dUda_(2,3)    ')
            call compare(duda(3, 1), apduda(3, 1), 'Atomic   dUda_(3,1)    ')
            call compare(duda(3, 2), apduda(3, 2), 'Atomic   dUda_(3,2)    ')
            call compare(duda(3, 3), apduda(3, 3), 'Atomic   dUda_(3,3)    ')
            write (6, *) '--------------------------------------------'
        end if  ! ( master )
    end if  ! ( chkvir /= 0)
    if (master) then
        call mexit(6, 0)
    else
        call mexit(0, 0)
    end if
60  format(9x, 'NUMERICAL, ANALYTICAL FORCES (diff) from atom ', i9)
66  format(1x, i6, f14.8, 1x, f14.8, 1x, f14.8)
67  format(1x, 3(1x, f14.8))
70  format(1x, 'RMS force error = ', e10.3)
167 format(1x, 'Checking numerical force for random atoms')
168 format(1x, 'Checking numerical force for user chosen atoms')
end subroutine debug_frc
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine crd_mod here]
subroutine crd_mod(crd, atomn, j, del, save)
    implicit none
    _REAL_ crd(3, *), del, save
    integer atomn, j

    save = crd(j, atomn)
    crd(j, atomn) = save + del
    return
end subroutine crd_mod
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine crd_rst here]
subroutine crd_rst(crd, atomn, j, del, save)
    implicit none
    _REAL_ crd(3, *), del, save
    integer atomn, j

    crd(j, atomn) = save
    return
end subroutine crd_rst
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_analfrc here]
subroutine get_analfrc(xx, ix, ih, ipairs, x, f, &
    vir, ene, qsetup)
    use stack
    use state
    implicit none
    character(kind=1, len=11) :: routine = "get_analfrc"
#  include "memory.h"
#  include "md.h"
#  include "ew_mpole.h"
#  include "ew_cntrl.h"
#ifdef MPI
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
    include 'mpif.h'
    integer ierr
#  include "parallel.h"
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif
    logical qsetup
    _REAL_ xx(*), x(*), f(3, *), vir(*)

    integer ipairs(*), ix(*)
    type(state_rec) :: ene
    character(len=4) ih(*)
    _REAL_ ekcmt(4)
    integer ltmp
    logical, parameter :: do_list_update = .false.

#ifdef MPI
    call mpi_barrier(commsander, ierr)
#endif
    call fix_xr(x, natom, nspm, ix(i70), xx(l75), &
        ekcmt, xx(l45), xx(lvel), xx(lmass))
    call force(xx, ix, ih, ipairs, x, f, ene, vir, &
        xx(l96), xx(l97), xx(l98), xx(l99), qsetup, do_list_update, 0)
#ifdef MPI
    call get_stack(ltmp, 3*natom, routine)
    if (.not. rstack_ok) then
        deallocate (r_stack)
        allocate (r_stack(1:lastrst), stat=alloc_ier)
        call reassign_rstack(routine)
    end if
    REQUIRE(rstack_ok)
    call merge_forces(f, r_stack(ltmp), natom)
    ! merge fields if dipoles
    if (mpoltype >= 1) then
        call merge_forces(xx(lfield), r_stack(ltmp), natom)
    end if
    call free_stack(ltmp, routine)
#endif
    return
end subroutine get_analfrc
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fix_xr here]
subroutine fix_xr(x, natom, nspm, nsp, tma, ekcmt, xr, v, amass)
    implicit none
    _REAL_ x(3, *), tma(*), ekcmt(*), &
        xr(3, *), v(*), amass(*)
    integer natom, nspm, nsp(*), i

    do i = 1, natom
        xr(1, i) = x(1, i)
        xr(2, i) = x(2, i)
        xr(3, i) = x(3, i)
    end do
!#ifdef MPI
    call ekcmr(nspm, nsp, tma, ekcmt, xr, v, amass, 1, natom)
!#else
!   call ekcmr(nspm,nsp,tma,ekcmt,xr,v,amass)
!#endif
    return
end subroutine fix_xr
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_numfrc here]
subroutine get_numfrc(xx, ix, ih, ipairs, x, f, &
    vir, ene, atomn, del, apfrc, qsetup)
    use state
    implicit none
#  include "memory.h"
#  include "md.h"

    logical qsetup
    _REAL_ xx(*), x(*), f(*), &
        vir(*), del, apfrc(3)

    type(state_rec) :: ene
    integer ipairs(*), ix(*), atomn
    character(len=4) ih(*)
    _REAL_ enep, enem, save, delta
    integer j
    logical, parameter :: do_list_update = .false.

    do j = 1, 3
        delta = del
        call crd_mod(x, atomn, j, delta, save)
        call force(xx, ix, ih, ipairs, x, f, ene, vir, &
            xx(l96), xx(l97), xx(l98), xx(l99), qsetup, do_list_update, 0)
        enep = ene%pot%tot
        call crd_rst(x, atomn, j, delta, save)
        delta = -del
        call crd_mod(x, atomn, j, delta, save)
        call force(xx, ix, ih, ipairs, x, f, ene, vir, &
            xx(l96), xx(l97), xx(l98), xx(l99), qsetup, do_list_update, 0)
        enem = ene%pot%tot
        call crd_rst(x, atomn, j, delta, save)
        ! force is negative of gradient
        apfrc(j) = -(enep - enem)/(2.d0*del)
    end do
    return
end subroutine get_numfrc
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rmsdiff here]
subroutine rmsdiff(apfrc, frc, rms)
    implicit none
    _REAL_ apfrc(3), frc(3)
    _REAL_ num, den, small, rms
    small = 1.d-6
    den = frc(1)**2 + frc(2)**2 + frc(3)**2 + small
    num = (apfrc(1) - frc(1))**2 + (apfrc(2) - frc(2))**2 + &
        (apfrc(3) - frc(3))**2
    rms = sqrt(num/den)
    return
end subroutine rmsdiff
!----------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rms_check here]
subroutine rms_check(natom, xx, ix, ih, ipairs, x, f, &
    ener, srms, f2, f3, f4, f5, f6, scg, scn1, scn2, sdip, cg, cn1, cn2, &
    fixdip, inddip, sind, qsetup)
    use ew_recip, only : frcx
    use parms, only : nttyp
    use state
    implicit none
#  include "ew_frc.h"
#  include "ew_cntrl.h"
#  include "ew_mpole.h"
#  include "extra.h"
#  include "flocntrl.h"
    _REAL_ xx(*), x(3, *), f(3, *)
    _REAL_ srms(*), f2(3, *), f3(3, *), f4(3, *), &
        f5(3, *), f6(3, *)
    _REAL_ scg(*), scn1(*), scn2(*), sdip(3, *), &
        cg(*), cn1(*), cn2(*), fixdip(3, *)
    type(state_rec) :: ener, ene
    integer ipairs(*), ix(*)
    character(len=4) ih(*)
    integer j, k, natom, jmax, matom, jmin
    _REAL_ tvec(16), teer, teed, teea, tees, &
        tevd, tevr, tfrcx, tfrcy, tfrcz
    _REAL_ vir, tvir
    _REAL_ trec_vir(3, 3), tdir_vir(3, 3), tadj_vir(3, 3)
    _REAL_ inddip(3, *), sind(3, *)
    _REAL_ trec_vird(3, 3), tatvir(3, 3), tmolvir(3, 3)
    _REAL_ tsubvir(3, 3), tf(3), tx(3)
    _REAL_ maxrms, rms, virial(4)
    integer sdo_dir, sdo_rec, sdo_adj, sdo_self, &
        sdo_bond, sdo_angle, sdo_ephi
    integer szerochg, szerovdw, szerodip, count, mycount
    logical qsetup

    mycount = 0
    if (master) then
        call amopen(30, 'forcedump.dat', 'O', 'F', 'R')
        write (6, *) 'CHECKING energies,forces against forcedump.dat'
        read (30, 299) matom
        if (matom /= natom) then
            write (6, *) 'RMS CHECK: atom numbers do not match with file'
            return
        end if
        maxrms = -1.d0
        jmax = -1
        do j = 1, natom
            read (30, 301) tx(1), tx(2), tx(3)
            rms = (tx(1) - x(1, j))**2 + (tx(2) - x(2, j))**2 + &
                (tx(3) - x(3, j))**2
            if (rms > maxrms) then
                maxrms = rms
                jmax = j
            end if
        end do
        write (6, *) 'First check new, old coordinates'
        write (6, *) '-----------------------------------------'
        write (6, 51) sqrt(maxrms), jmax
51      format(1x, 'Maximum atomic RMS coordinate error: ', e9.3, &
            ' at atom number ', i9)

        read (30, 200) count
        if (mycount /= count) goto 666
        do j = 1, 16
            read (30, 300) tvec(j)
        end do
        read (30, 301) teer, teed, teea
        read (30, 301) tees, tevd, tevr
        read (30, 301) tfrcx, tfrcy, tfrcz
        write (6, *) '-----------------------------------------'
        write (6, *) 'ERRORS in Energies, virials and forces'
        write (6, *) 'Compared to stored results in forcedump.dat'
        write (6, *) '-----------------------------------------'
        write (6, *) 'Relative errors in energy components:'
        write (6, *) '-----------------------------------------'
        call compare(ener%pot%tot, tvec(1), 'Potential energy:      ')
        call compare(ener%pot%vdw, tvec(2), 'Nonbond energy:        ')
        call compare(ener%pot%elec, tvec(3), 'Electrostatic energy:  ')
        call compare(ener%pot%bond, tvec(5), 'Bond energy:           ')
        call compare(ener%pot%angle, tvec(6), 'Angle energy:          ')
        call compare(ener%pot%dihedral, tvec(7), 'Dihedral energy:       ')
        call compare(ener%pot%vdw_14, tvec(8), '1-4VDW energy:         ')
        call compare(ener%pot%elec_14, tvec(9), '1-4EE energy:          ')
        call compare(ener%pot%constraint, tvec(10), 'Constraint energy:     ')
        call compare(ener%pot%polar, tvec(11), 'Polarization energy:   ')
        write (6, *) '-----------------------------------------'
        write (6, *) 'Relative errors in ewald energy components:'
        write (6, *) '-----------------------------------------'
        call compare(eed, teed, 'Direct ee:             ')
        call compare(eer, teer, 'Reciprocal ee:         ')
        call compare(eea, teea, 'Adjust ee:             ')
        call compare(ees, tees, 'Self ee:               ')
        call compare(evdw - evdwr, tevd, 'Direct vdw:            ')
        call compare(evdwr, tevr, 'Long-range vdw:        ')
        write (6, *) '-----------------------------------------'
        write (6, *) 'Relative errors in total force components:'
        write (6, *) '-----------------------------------------'
        call compare(frcx(1), tfrcx, 'FRCX:                  ')
        call compare(frcx(2), tfrcy, 'FRCY:                  ')
        call compare(frcx(3), tfrcz, 'FRCZ:                  ')

        read (30, 200) count
        mycount = mycount + 1
        if (mycount /= count) goto 666
        read (30, 301) trec_vir(1, 1), trec_vir(1, 2), trec_vir(1, 3)
        read (30, 301) trec_vir(2, 1), trec_vir(2, 2), trec_vir(2, 3)
        read (30, 301) trec_vir(3, 1), trec_vir(3, 2), trec_vir(3, 3)
        read (30, 301) trec_vird(1, 1), trec_vird(1, 2), trec_vird(1, 3)
        read (30, 301) trec_vird(2, 1), trec_vird(2, 2), trec_vird(2, 3)
        read (30, 301) trec_vird(3, 1), trec_vird(3, 2), trec_vird(3, 3)
        read (30, 301) tdir_vir(1, 1), tdir_vir(1, 2), tdir_vir(1, 3)
        read (30, 301) tdir_vir(2, 1), tdir_vir(2, 2), tdir_vir(2, 3)
        read (30, 301) tdir_vir(3, 1), tdir_vir(3, 2), tdir_vir(3, 3)
        read (30, 301) tadj_vir(1, 1), tadj_vir(1, 2), tadj_vir(1, 3)
        read (30, 301) tadj_vir(2, 1), tadj_vir(2, 2), tadj_vir(2, 3)
        read (30, 301) tadj_vir(3, 1), tadj_vir(3, 2), tadj_vir(3, 3)
        read (30, 301) tatvir(1, 1), tatvir(1, 2), tatvir(1, 3)
        read (30, 301) tatvir(2, 1), tatvir(2, 2), tatvir(2, 3)
        read (30, 301) tatvir(3, 1), tatvir(3, 2), tatvir(3, 3)
        read (30, 301) tmolvir(1, 1), tmolvir(1, 2), tmolvir(1, 3)
        read (30, 301) tmolvir(2, 1), tmolvir(2, 2), tmolvir(2, 3)
        read (30, 301) tmolvir(3, 1), tmolvir(3, 2), tmolvir(3, 3)
        read (30, 301) tsubvir(1, 1), tsubvir(1, 2), tsubvir(1, 3)
        read (30, 301) tsubvir(2, 1), tsubvir(2, 2), tsubvir(2, 3)
        read (30, 301) tsubvir(3, 1), tsubvir(3, 2), tsubvir(3, 3)
        write (6, *) '-----------------------------------------'
        write (6, *) 'Relative errors in Virial Tensors:'
        write (6, *) '-----------------------------------------'
        call compare_vir(rec_vir, trec_vir, &
            'Reciprocal ee virial:  ')
        call compare_vir(dir_vir, tdir_vir, &
            'Direct virial:         ')
        call compare_vir(adj_vir, tadj_vir, &
            'Adjust ee virial:      ')
        call compare_vir(atvir, tatvir, &
            'Nonbond atomic virial: ')
        call compare_vir(molvir, tmolvir, &
            'Nonbond Molec. virial: ')
        write (6, *) '-----------------------------------------'
    end if  ! ( master )
    ! end master only
    ! restore orig cg,cn1,cn2
    call save_flow(sdo_dir, sdo_rec, sdo_adj, sdo_self, &
        sdo_bond, sdo_angle, sdo_ephi, szerochg, szerodip, szerovdw)
    call array_copy(scg, cg, natom)
    call array_copy(scn1, cn1, nttyp)
    call array_copy(scn2, cn2, nttyp)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(sind, inddip, 3*natom)
    end if
    if (szerochg == 1) then
        call zero_array(cg, natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. szerodip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    if (szerovdw == 1) then
        call zero_array(cn1, nttyp)
        call zero_array(cn2, nttyp)
    end if
    ! redo forces
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    if (master) then
        call array_copy(f, f3, 3*natom)
        read (30, 200) count
        mycount = mycount + 1
        if (mycount /= count) goto 666
        do j = 1, natom
            read (30, 301) f4(1, j), f4(2, j), f4(3, j)
        end do
        write (6, *) 'TOTAL FORCES:     '
        write (6, *) '-----------------------------------------'
        call frc_compare(f3, f4, f3, f4, srms, natom)
        ! field if dipoles
        if (mpoltype >= 1) then
            call array_copy(xx(lfield), f5, 3*natom)
            read (30, 200) count
            mycount = mycount + 1
            if (mycount /= count) goto 666
            do j = 1, natom
                read (30, 301) f6(1, j), f6(2, j), f6(3, j)
            end do
            write (6, *) 'TOTAL FIELDS:     '
            write (6, *) '-----------------------------------------'
            call frc_compare(f5, f6, f5, f6, srms, natom)
        end if
    end if
    ! .. DO FORCE COMPONENTS
    ! bonds
    call zero_flow()
    do_bond = sdo_bond
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    if (master) then
        read (30, 200) count
        mycount = mycount + 1
        if (mycount /= count) goto 666
        do j = 1, natom
            read (30, 301) f2(1, j), f2(2, j), f2(3, j)
        end do
        write (6, *) 'FORCES from BOND ROUTINE:     '
        write (6, *) '-----------------------------------------'
        call frc_compare(f, f2, f3, f4, srms, natom)
    end if
    ! angles
    call zero_flow()
    do_angle = sdo_angle
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    if (master) then
        read (30, 200) count
        mycount = mycount + 1
        if (mycount /= count) goto 666
        do j = 1, natom
            read (30, 301) f2(1, j), f2(2, j), f2(3, j)
        end do
        write (6, *) 'FORCES from ANGLE ROUTINE:     '
        write (6, *) '-----------------------------------------'
        call frc_compare(f, f2, f3, f4, srms, natom)
    end if
    ! dihedrals
    call zero_flow()
    do_ephi = sdo_ephi
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    if (master) then
        read (30, 200) count
        mycount = mycount + 1
        if (mycount /= count) goto 666
        do j = 1, natom
            read (30, 301) f2(1, j), f2(2, j), f2(3, j)
        end do
        write (6, *) 'FORCES from DIHEDRAL ROUTINE:     '
        write (6, *) '-----------------------------------------'
        call frc_compare(f, f2, f3, f4, srms, natom)
    end if
    ! direct vdw sum
    call zero_flow()
    do_dir = sdo_dir
    ! zero the charges
    call zero_array(cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    ! restore the charges
    call array_copy(scg, cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(sind, inddip, 3*natom)
    end if
    if (master) then
        read (30, 200) count
        mycount = mycount + 1
        if (mycount /= count) goto 666
        do j = 1, natom
            read (30, 301) f2(1, j), f2(2, j), f2(3, j)
        end do
        write (6, *) 'FORCES from Direct VDW ROUTINE:     '
        write (6, *) '-----------------------------------------'
        call frc_compare(f, f2, f3, f4, srms, natom)
    end if
    ! direct ee sum
    call zero_flow()
    do_dir = sdo_dir
    ! zero the vdw
    call zero_array(cn1, nttyp)
    call zero_array(cn2, nttyp)
    if (szerochg == 1) then
        call zero_array(cg, natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. szerodip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    ! restore the vdw
    call array_copy(scn1, cn1, nttyp)
    call array_copy(scn2, cn2, nttyp)
    ! restore the charges
    call array_copy(scg, cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(sind, inddip, 3*natom)
    end if
    if (master) then
        read (30, 200) count
        mycount = mycount + 1
        if (mycount /= count) goto 666
        do j = 1, natom
            read (30, 301) f2(1, j), f2(2, j), f2(3, j)
        end do
        write (6, *) 'FORCES from Direct EE ROUTINE:     '
        write (6, *) '-----------------------------------------'
        call frc_compare(f, f2, f3, f4, srms, natom)
        ! field if dipoles
        if (mpoltype >= 1) then
            read (30, 200) count
            mycount = mycount + 1
            if (mycount /= count) goto 666
            do j = 1, natom
                read (30, 301) f2(1, j), f2(2, j), f2(3, j)
            end do
            write (6, *) 'FIELDS from Direct EE ROUTINE:     '
            write (6, *) '-----------------------------------------'
            call frc_compare(xx(lfield), f2, f5, f6, srms, natom)
        end if
    end if
    ! reciprocal ee sum
    call zero_flow()
    do_rec = sdo_rec
    if (szerochg == 1) then
        call zero_array(cg, natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. szerodip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    ! restore the charges
    call array_copy(scg, cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(sind, inddip, 3*natom)
    end if
    if (master) then
        read (30, 200) count
        mycount = mycount + 1
        if (mycount /= count) goto 666
        do j = 1, natom
            read (30, 301) f2(1, j), f2(2, j), f2(3, j)
        end do
        write (6, *) 'FORCES from Reciprocal EE ROUTINE:     '
        write (6, *) '-----------------------------------------'
        call frc_compare(f, f2, f3, f4, srms, natom)
        ! field if dipoles
        if (mpoltype >= 1) then
            read (30, 200) count
            mycount = mycount + 1
            if (mycount /= count) goto 666
            do j = 1, natom
                read (30, 301) f2(1, j), f2(2, j), f2(3, j)
            end do
            write (6, *) 'FIELDS from Reciprocal EE ROUTINE:     '
            write (6, *) '-----------------------------------------'
            call frc_compare(xx(lfield), f2, f5, f6, srms, natom)
        end if
    end if
    ! adj ee sum
    call zero_flow()
    do_adj = sdo_adj
    if (szerochg == 1) then
        call zero_array(cg, natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. szerodip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    if (master) then
        read (30, 200) count
        mycount = mycount + 1
        if (mycount /= count) goto 666
        do j = 1, natom
            read (30, 301) f2(1, j), f2(2, j), f2(3, j)
        end do
        write (6, *) 'FORCES from Adjust EE ROUTINE:     '
        write (6, *) '-----------------------------------------'
        call frc_compare(f, f2, f3, f4, srms, natom)
        ! field if dipoles
        if (mpoltype >= 1) then
            read (30, 200) count
            mycount = mycount + 1
            if (mycount /= count) goto 666
            do j = 1, natom
                read (30, 301) f2(1, j), f2(2, j), f2(3, j)
            end do
            write (6, *) 'FIELDS from Adjust EE ROUTINE:     '
            write (6, *) '-----------------------------------------'
            call frc_compare(xx(lfield), f2, f5, f6, srms, natom)
        end if
    end if
    ! fix back again (side effects in force routine)
    call restore_flow(sdo_dir, sdo_rec, sdo_adj, sdo_self, &
        sdo_bond, sdo_angle, sdo_ephi, szerochg, szerodip, szerovdw)
    call array_copy(scg, cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(sind, inddip, 3*natom)
    end if
    call array_copy(scn1, cn1, nttyp)
    call array_copy(scn2, cn2, nttyp)
    if (szerochg == 1) then
        call zero_array(cg, natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. szerodip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    if (szerovdw == 1) then
        call zero_array(cn1, nttyp)
        call zero_array(cn2, nttyp)
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)

200 format(1x, i2)
299 format(1x, i8)
300 format(1x, e23.16)
301 format(3(1x, e23.16))
305 format(1x, 7i5)
306 format(1x, 2i5)
    close (30)
    return
666 continue
    write (6, *) 'Misalignment of forcedump (diff versions?)'
    close (30)
    return
end subroutine rms_check
!--------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine compare here]
subroutine compare(a, b, string)
    implicit none
    _REAL_ a, b, small, err, den
    character(len=*) string
    small = 1.d-11
    den = max(abs(a), abs(b)) + small
    err = abs(a - b)/den
    write (6, 60) string, a, b, err
60  format(1x, a, e12.5, ' vs ', e12.5, ' relative error: ', e9.3)
    return
end subroutine compare
!--------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine compare_vir here]
subroutine compare_vir(u, v, string)
    implicit none
    _REAL_ u(3, 3), v(3, 3), small, err, den
    integer i, j
    character(len=*) string
    small = 1.d-11
    do j = 1, 3
        do i = 1, 3
            den = max(abs(u(i, j)), abs(v(i, j))) + small
            err = abs(u(i, j) - v(i, j))/den
            if (i == 1 .and. j == 1) then
                write (6, 60) string, u(i, j), v(i, j), err
            else
                write (6, 61) u(i, j), v(i, j), err
            end if
        end do
    end do
60  format(1x, a, e12.5, ' vs ', e12.5, ' relative error: ', e9.3)
61  format(24x, e12.5, ' vs ', e12.5, ' relative error: ', e9.3)
    return
end subroutine compare_vir
!--------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine force_dump here]
subroutine force_dump(natom, xx, ix, ih, ipairs, x, f, &
    ener, scg, scn1, scn2, sdip, cg, cn1, cn2, fixdip, &
    inddip, sind, qsetup)
    use ew_recip, only : frcx
    use parms, only : nttyp
    use state
    implicit none
#  include "ew_frc.h"
#  include "ew_cntrl.h"
#  include "ew_mpole.h"
#  include "extra.h"
#  include "debug.h"
#  include "flocntrl.h"
    _REAL_ xx(*), x(3, *), f(3, *)
    _REAL_ scg(*), scn1(*), scn2(*), sdip(3, *), &
        cg(*), cn1(*), cn2(*), fixdip(3, *)
    integer ipairs(*), ix(*)
    character(len=4) ih(*)
    integer j, natom, kk
    _REAL_ inddip(3, *), sind(3, *)
    integer sdo_dir, sdo_rec, sdo_adj, sdo_self, &
        sdo_bond, sdo_angle, sdo_ephi
    integer szerochg, szerovdw, szerodip, count
    _REAL_ virial(4)
    type(state_rec) :: ene
    type(state_rec) :: ener
    logical qsetup

    if (master) then
        call amopen(30, 'forcedump.dat', 'N', 'F', 'W')
        write (30, 299) natom
        do j = 1, natom
            write (30, 301) x(1, j), x(2, j), x(3, j)
        end do
        count = 0
        write (30, 199) count
199     format(1x, i2, ' START of Energies')

        write (30, 300) ener%pot%tot
        write (30, 300) ener%pot%vdw
        write (30, 300) ener%pot%elec
        write (30, 300) ener%pot%gb
        write (30, 300) ener%pot%bond
        write (30, 300) ener%pot%angle
        write (30, 300) ener%pot%dihedral
        write (30, 300) ener%pot%vdw_14
        write (30, 300) ener%pot%elec_14
        write (30, 300) ener%pot%constraint !10
        write (30, 300) ener%pot%polar
        write (30, 300) ener%aveper
        write (30, 300) ener%aveind
        write (30, 300) ener%avetot
        write (30, 300) ener%pot%surf
        write (30, 300) ener%pot%dvdl       !16

        write (30, 301) eer, eed, eea
        write (30, 301) ees, evdw - evdwr, evdwr
        write (30, 301) frcx(1), frcx(2), frcx(3)
        count = 1
        write (30, 200) count
200     format(1x, i2, ' START of VIRIALS')
        write (30, 301) rec_vir(1, 1), rec_vir(1, 2), rec_vir(1, 3)
        write (30, 301) rec_vir(2, 1), rec_vir(2, 2), rec_vir(2, 3)
        write (30, 301) rec_vir(3, 1), rec_vir(3, 2), rec_vir(3, 3)
        write (30, 301) rec_vird(1, 1), rec_vird(1, 2), rec_vird(1, 3)
        write (30, 301) rec_vird(2, 1), rec_vird(2, 2), rec_vird(2, 3)
        write (30, 301) rec_vird(3, 1), rec_vird(3, 2), rec_vird(3, 3)
        write (30, 301) dir_vir(1, 1), dir_vir(1, 2), dir_vir(1, 3)
        write (30, 301) dir_vir(2, 1), dir_vir(2, 2), dir_vir(2, 3)
        write (30, 301) dir_vir(3, 1), dir_vir(3, 2), dir_vir(3, 3)
        write (30, 301) adj_vir(1, 1), adj_vir(1, 2), adj_vir(1, 3)
        write (30, 301) adj_vir(2, 1), adj_vir(2, 2), adj_vir(2, 3)
        write (30, 301) adj_vir(3, 1), adj_vir(3, 2), adj_vir(3, 3)
        write (30, 301) atvir(1, 1), atvir(1, 2), atvir(1, 3)
        write (30, 301) atvir(2, 1), atvir(2, 2), atvir(2, 3)
        write (30, 301) atvir(3, 1), atvir(3, 2), atvir(3, 3)
        write (30, 301) molvir(1, 1), molvir(1, 2), molvir(1, 3)
        write (30, 301) molvir(2, 1), molvir(2, 2), molvir(2, 3)
        write (30, 301) molvir(3, 1), molvir(3, 2), molvir(3, 3)
        write (30, 301) subvir(1, 1), subvir(1, 2), subvir(1, 3)
        write (30, 301) subvir(2, 1), subvir(2, 2), subvir(2, 3)
        write (30, 301) subvir(3, 1), subvir(3, 2), subvir(3, 3)
    end if  ! ( master )
    ! save current values... do force components
    call save_flow(sdo_dir, sdo_rec, sdo_adj, sdo_self, &
        sdo_bond, sdo_angle, sdo_ephi, szerochg, szerodip, szerovdw)
    ! restore orig cg,cn1,cn2
    call array_copy(scg, cg, natom)
    call array_copy(scn1, cn1, nttyp)
    call array_copy(scn2, cn2, nttyp)
    if (szerochg == 1) then
        call zero_array(cg, natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. szerodip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    if (szerovdw == 1) then
        call zero_array(cn1, nttyp)
        call zero_array(cn2, nttyp)
    end if
    ! redo forces
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    if (master) then
        count = count + 1
        write (30, 201) count
201     format(1x, i2, ' Total Force')
        do j = 1, natom
            write (30, 301) f(1, j), f(2, j), f(3, j)
        end do
        ! field if dipoles
        if (mpoltype >= 1) then
            count = count + 1
            write (30, 2010) count
2010        format(1x, i2, ' Total Field')
            kk = 0
            do j = 1, natom
                write (30, 301) xx(lfield + kk), xx(lfield + kk + 1), xx(lfield + kk + 2)
                kk = kk + 3
            end do
        end if
    end if
    ! bonds
    call zero_flow()
    do_bond = sdo_bond
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    if (master) then
        count = count + 1
        write (30, 202) count
202     format(1x, i2, ' Bond Force')
        do j = 1, natom
            write (30, 301) f(1, j), f(2, j), f(3, j)
        end do
    end if
    ! angles
    call zero_flow()
    do_angle = sdo_angle
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    if (master) then
        count = count + 1
        write (30, 203) count
203     format(1x, i2, ' Angle Force')
        do j = 1, natom
            write (30, 301) f(1, j), f(2, j), f(3, j)
        end do
    end if
    ! dihedrals
    call zero_flow()
    do_ephi = sdo_ephi
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    if (master) then
        count = count + 1
        write (30, 204) count
204     format(1x, i2, ' Dihedral Force')
        do j = 1, natom
            write (30, 301) f(1, j), f(2, j), f(3, j)
        end do
    end if
    ! direct vdw sum
    call zero_flow()
    do_dir = sdo_dir
    ! zero the charges
    call zero_array(cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    ! restore the charges
    call array_copy(scg, cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(sind, inddip, 3*natom)
    end if
    if (master) then
        count = count + 1
        write (30, 205) count
205     format(1x, i2, ' Van der Waals Force')
        do j = 1, natom
            write (30, 301) f(1, j), f(2, j), f(3, j)
        end do
    end if
    ! direct ee sum
    call zero_flow()
    do_dir = sdo_dir
    ! zero the vdw
    call zero_array(cn1, nttyp)
    call zero_array(cn2, nttyp)
    ! zero the charges
    if (szerochg == 1) then
        call zero_array(cg, natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. szerodip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    ! restore the vdw
    call array_copy(scn1, cn1, nttyp)
    call array_copy(scn2, cn2, nttyp)
    ! restore the charges
    call array_copy(scg, cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(sind, inddip, 3*natom)
    end if
    if (master) then
        count = count + 1
        write (30, 206) count
206     format(1x, i2, ' Direct sum electrostatic Force')
        do j = 1, natom
            write (30, 301) f(1, j), f(2, j), f(3, j)
        end do
        ! field if dipoles
        if (mpoltype >= 1) then
            count = count + 1
            write (30, 2060) count
2060        format(1x, i2, ' Direct sum electrostatic Field')
            kk = 0
            do j = 1, natom
                write (30, 301) xx(lfield + kk), xx(lfield + kk + 1), xx(lfield + kk + 2)
                kk = kk + 3
            end do
        end if
    end if
    ! reciprocal ee sum
    call zero_flow()
    do_rec = sdo_rec
    ! zero the charges
    if (szerochg == 1) then
        call zero_array(cg, natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. szerodip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    ! restore the charges
    call array_copy(scg, cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(sind, inddip, 3*natom)
    end if
    if (master) then
        count = count + 1
        write (30, 207) count
207     format(1x, i2, ' Reciprocal sum electrostatic Force')
        do j = 1, natom
            write (30, 301) f(1, j), f(2, j), f(3, j)
        end do
        ! field if dipoles
        if (mpoltype >= 1) then
            count = count + 1
            write (30, 2070) count
2070        format(1x, i2, ' Reciprocal sum electrostatic Field')
            kk = 0
            do j = 1, natom
                write (30, 301) xx(lfield + kk), xx(lfield + kk + 1), xx(lfield + kk + 2)
                kk = kk + 3
            end do
        end if
    end if
    ! adj ee sum
    call zero_flow()
    do_adj = sdo_adj
    ! zero the charges
    if (szerochg == 1) then
        call zero_array(cg, natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. szerodip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
    ! restore the charges
    call array_copy(scg, cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(sind, inddip, 3*natom)
    end if
    if (master) then
        count = count + 1
        write (30, 208) count
208     format(1x, i2, ' Adjust sum electrostatic Force')
        do j = 1, natom
            write (30, 301) f(1, j), f(2, j), f(3, j)
        end do
        ! field if dipoles
        if (mpoltype >= 1) then
            count = count + 1
            write (30, 2080) count
2080        format(1x, i2, ' Adjust sum electrostatic Field')
            kk = 0
            do j = 1, natom
                write (30, 301) xx(lfield + kk), xx(lfield + kk + 1), xx(lfield + kk + 2)
                kk = kk + 3
            end do
        end if
    end if
    ! end components
    if (master) then
        close (30)
    end if
    ! fix back again (side effects in force routine)
    call restore_flow(sdo_dir, sdo_rec, sdo_adj, sdo_self, &
        sdo_bond, sdo_angle, sdo_ephi, szerochg, szerodip, szerovdw)
    call array_copy(scg, cg, natom)
    if (indmeth == 3 .and. irstdip == 1) then
        call array_copy(sind, inddip, 3*natom)
    end if
    call array_copy(scn1, cn1, nttyp)
    call array_copy(scn2, cn2, nttyp)
    if (szerochg == 1) then
        call zero_array(cg, natom)
    end if
    if (indmeth == 3 .and. irstdip == 1 .and. szerodip == 1) then
        call zero_array(inddip, 3*natom)
    end if
    if (szerovdw == 1) then
        call zero_array(cn1, nttyp)
        call zero_array(cn2, nttyp)
    end if
    call get_analfrc(xx, ix, ih, ipairs, x, f, &
        virial, ene, qsetup)
299 format(1x, i8)
300 format(1x, e23.16)
301 format(3(1x, e23.16))
305 format(1x, 7i5)
306 format(1x, 2i5)
    return
end subroutine force_dump
!--------------------------------------------------
#ifdef MPI

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine merge_forces here]
subroutine merge_forces(f, ftmp, natom)
    implicit none
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
    include 'mpif.h'
    integer ierr
#  include "parallel.h"
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
    ! needed in mpi case to put forces together
    integer natom, init
    _REAL_ f(*), ftmp(*)
    integer i
    call MPI_Gatherv(f(iparpt3(sanderrank) + 1), &
        rcvcnt3(sanderrank), &
        MPI_DOUBLE_PRECISION, &
        f, rcvcnt3(0:sandersize), iparpt3(0:sandersize), &
        MPI_DOUBLE_PRECISION, &
        0, commsander, ierr)
    return
end subroutine merge_forces
#endif
!--------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine check_virial here]
subroutine check_virial(xx, ix, ih, ipairs, x, f, &
    ene, del, dudv, type, qsetup)
    use nblist, only : volume
    implicit none
#  include "memory.h"
#  include "md.h"

    _REAL_ xx(*), x(*), f(*), &
        ene(*), del, dudv
    integer ipairs(*), ix(*), type
    character(len=4) ih(*)
    _REAL_ factor(3), term, vir(4)
    _REAL_ enep, enem, volp, volm
    logical, intent(in) :: qsetup
    logical, parameter :: do_list_update = .false.

    ! scale cell up
    term = (1.d0 + del)**(1.d0/3.d0)
    factor(1) = term
    factor(2) = term
    factor(3) = term
    call redo_ucell(factor)
    call ew_pscale(natom, x, xx(lmass), nspm, ix(i70), type)
    call force(xx, ix, ih, ipairs, x, f, ene(23), vir, &
        xx(l96), xx(l97), xx(l98), xx(l99), qsetup, do_list_update, 0)
    enep = ene(23)
    volp = volume
    ! put cell back
    term = 1.d0/term
    factor(1) = term
    factor(2) = term
    factor(3) = term
    call redo_ucell(factor)
    call ew_pscale(natom, x, xx(lmass), nspm, ix(i70), type)
    ! now scale cell down
    term = (1.d0 - del)**(1.d0/3.d0)
    factor(1) = term
    factor(2) = term
    factor(3) = term
    call redo_ucell(factor)
    call ew_pscale(natom, x, xx(lmass), nspm, ix(i70), type)
    call force(xx, ix, ih, ipairs, x, f, ene(23), vir, &
        xx(l96), xx(l97), xx(l98), xx(l99), qsetup, do_list_update, 0)
    enem = ene(23)
    volm = volume
    dudv = (enep - enem)/(volp - volm)
    ! put cell back
    term = 1.d0/term
    factor(1) = term
    factor(2) = term
    factor(3) = term
    call redo_ucell(factor)
    call ew_pscale(natom, x, xx(lmass), nspm, ix(i70), type)
    call force(xx, ix, ih, ipairs, x, f, ene(23), vir, &
        xx(l96), xx(l97), xx(l98), xx(l99), qsetup, do_list_update, 0)
    return
end subroutine check_virial
!--------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine check_vtens here]
subroutine check_vtens(xx, ix, ih, ipairs, x, f, &
    ene, del, duda, vten, apduda, type, qsetup)
    use nblist, only : ucell, recip
    implicit none
#  include "memory.h"
#  include "md.h"

    _REAL_ xx(*), x(*), f(*), &
        ene(*), del, duda(3, 3), vir(4), vten(3, 3), apduda(3, 3)
    integer ipairs(*), ix(*), type
    character(len=4) ih(*)
    _REAL_ enep, enem, save(3, 3), new(3, 3)
    integer i, j, k, l
    logical, intent(in) :: qsetup

    logical, parameter :: do_list_update = .false.
    ! save existing cell
    do i = 1, 3
        do j = 1, 3
            save(i, j) = ucell(i, j)
        end do
    end do
    do i = 1, 3
        do j = 1, 3
            ! get apduda(i,j)
            do k = 1, 3
                do l = 1, 3
                    new(k, l) = save(k, l)
                end do
            end do
            new(i, j) = save(i, j) + del
            call new_ucell_gen(new)
            call ew_pscale(natom, x, xx(lmass), nspm, ix(i70), type)
            call force(xx, ix, ih, ipairs, x, f, ene(23), vir, &
                xx(l96), xx(l97), xx(l98), xx(l99), qsetup, do_list_update, 0)
            enep = ene(23)
            new(i, j) = save(i, j) - del
            call new_ucell_gen(new)
            call ew_pscale(natom, x, xx(lmass), nspm, ix(i70), type)
            call force(xx, ix, ih, ipairs, x, f, ene(23), vir, &
                xx(l96), xx(l97), xx(l98), xx(l99), qsetup, do_list_update, 0)
            enem = ene(23)
            apduda(i, j) = (enep - enem)/(2.d0*del)
        end do
    end do
    ! put cell back
    do k = 1, 3
        do l = 1, 3
            new(k, l) = save(k, l)
        end do
    end do
    call new_ucell_gen(new)
    call ew_pscale(natom, x, xx(lmass), nspm, ix(i70), type)
    call force(xx, ix, ih, ipairs, x, f, ene(23), vir, &
        xx(l96), xx(l97), xx(l98), xx(l99), qsetup, do_list_update, 0)
    ! calc approx duda
    do i = 1, 3
        do j = 1, 3
            duda(i, j) = 0.d0
            do k = 1, 3
                ! recip is transpose of inverse of ucell
                duda(i, j) = duda(i, j) + vten(i, k)*recip(k, j)
            end do
        end do
    end do
    return
end subroutine check_vtens
!--------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_intvir here]
subroutine get_intvir(natom, x, f, subvir, intvir)
    implicit none
    _REAL_ x(3, *), f(3, *), subvir(3, 3), &
        intvir(3, 3)
    integer natom
    integer i, j, n
    _REAL_ mecvir(3, 3)
    do j = 1, 3
        do i = 1, 3
            mecvir(i, j) = 0.d0
        end do
    end do
    do n = 1, natom
        do j = 1, 3
            do i = 1, 3
                mecvir(i, j) = mecvir(i, j) - f(i, n)*x(j, n)
            end do
        end do
    end do
    do j = 1, 3
        do i = 1, 3
            intvir(i, j) = mecvir(i, j) - subvir(i, j)
        end do
    end do
    return
end subroutine get_intvir
!--------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine frc_compare here]
subroutine frc_compare(f1, f2, f3, f4, srms, natom)
    implicit none
    _REAL_ f1(3, *), f2(3, *), f3(3, *), &
        f4(3, *), srms(*)
    integer natom

    integer i, j, k, jmax, jmin
    _REAL_ num, den, totnum, totden, &
        totden1, totden2, den1, den2, rms, maxrms, minrms, small

    small = 1.d-10
    maxrms = -1.d0
    minrms = 10.d0
    jmax = -1
    totnum = 0.d0
    write (6, *) 'ABSOLUTE RMS deviations: (in kcal/(mol Ang))'
    do j = 1, natom
        num = 0.d0
        do i = 1, 3
            num = num + (f1(i, j) - f2(i, j))**2
        end do
        rms = sqrt(num/3)
        srms(j) = rms
        totnum = totnum + num
        if (rms > maxrms) then
            jmax = j
            maxrms = rms
        end if
        if (rms < minrms) then
            jmin = j
            minrms = rms
        end if
    end do
    write (6, 61) maxrms, jmax
    write (6, 62) minrms, jmin
    call do_percentiles(srms, maxrms, minrms, natom)
    write (6, 60) sqrt(totnum/(3*natom))
61  format(1x, 'Maximum atomic RMS force error: ', e9.3, &
        ' at atom number ', i9)
62  format(1x, 'Minimum atomic RMS force error: ', e9.3, &
        ' at atom number ', i9)
60  format(1x, 'Overall RMS FORCE ERROR:  ', e9.3)
    write (6, *) '-----------------------------------------'
    maxrms = -1.d0
    minrms = 10.d0
    jmax = -1
    totnum = 0.d0
    totden1 = 0.d0
    totden2 = 0.d0
    write (6, *) 'RELATIVE RMS deviations: '
    do j = 1, natom
        num = 0.d0
        den1 = 0.d0
        den2 = 0.d0
        do i = 1, 3
            num = num + (f1(i, j) - f2(i, j))**2
            den1 = den1 + f3(i, j)**2
            den2 = den2 + f4(i, j)**2
        end do
        den = small
        if (den1 > den) den = den1
        if (den2 > den) den = den2
        rms = sqrt(num/den)
        srms(j) = rms
        totnum = totnum + num
        totden1 = totden1 + den1
        totden2 = totden2 + den2
        if (rms > maxrms) then
            jmax = j
            maxrms = rms
        end if
        if (rms < minrms) then
            jmin = j
            minrms = rms
        end if
    end do
    totden = small
    if (totden1 > totden) totden = totden1
    if (totden2 > totden) totden = totden2
    write (6, 61) maxrms, jmax
    write (6, 62) minrms, jmin
    call do_percentiles(srms, maxrms, minrms, natom)
    write (6, 60) sqrt(totnum/totden)
    write (6, *) '-----------------------------------------'

    return
end subroutine frc_compare
!------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine do_percentiles here]
subroutine do_percentiles(srms, max, min, natom)
    implicit none
    _REAL_ srms(*), max, min
    integer natom
    integer count(500), totcount, med, sev5, ninety, ninety5, ninety9
    integer imed, isev5, ininety, ininety5, ininety9, i, j, k
    _REAL_ apmed, apsev5, apninety, &
        apninety5, apninety9, range, up, lo, frac, tiny

    tiny = 1.d-16
    if (max < tiny) then
        write (6, 71)
        write (6, *) 'RMS error range too small to get percentiles!'
        return
    end if
    med = natom/2
    sev5 = (3*natom)/4
    ninety = (9*natom)/10
    ninety5 = (95*natom)/100
    ninety9 = (99*natom)/100
    ! bin rms to quick estimate median,75th,90th,95th percentiles
    range = max - min
    do j = 1, 500
        count(j) = 0
    end do
    do j = 1, natom
        k = 500*(srms(j) - min)/range + 1
        count(k) = count(k) + 1
    end do
    totcount = 0
    imed = 0
    isev5 = 0
    ininety = 0
    ininety5 = 0
    ininety9 = 0
    do j = 1, 500
        up = min + j*range/500
        lo = up - range/500
        if (imed == 0 .and. totcount + count(j) > med) then
            imed = 1
            frac = (dble(med) - dble(totcount))/count(j)
            apmed = lo + frac*range/500
        end if
        if (isev5 == 0 .and. totcount + count(j) > sev5) then
            isev5 = 1
            frac = (dble(sev5) - dble(totcount))/count(j)
            apsev5 = lo + frac*range/500
        end if
        if (ininety == 0 .and. totcount + count(j) > ninety) then
            ininety = 1
            frac = (dble(ninety) - dble(totcount))/count(j)
            apninety = lo + frac*range/500
        end if
        if (ininety5 == 0 .and. totcount + count(j) > ninety5) then
            ininety5 = 1
            frac = (dble(ninety5) - dble(totcount))/count(j)
            apninety5 = lo + frac*range/500
        end if
        if (ininety9 == 0 .and. totcount + count(j) > ninety9) then
            ininety9 = 1
            frac = (dble(ninety9) - dble(totcount))/count(j)
            apninety9 = lo + frac*range/500
        end if
        totcount = totcount + count(j)
    end do
    write (6, 71)
71  format(1x, 'The 50th,75th,90th,95th and 99th percentile', &
        '  Atomic RMS force errors:')
    write (6, 70) apmed, apsev5, apninety, apninety5, apninety9
70  format(5(1x, e9.3))
    return
end subroutine do_percentiles
!------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine save_flow here]
subroutine save_flow(sdo_dir, sdo_rec, sdo_adj, sdo_self, &
    sdo_bond, sdo_angle, sdo_ephi, szerochg, szerodip, szerovdw)
    implicit none
#  include "flocntrl.h"
#  include "debug.h"
    integer sdo_dir, sdo_rec, sdo_adj, sdo_self, &
        sdo_bond, sdo_angle, sdo_ephi
    integer szerochg, szerovdw, szerodip
    sdo_dir = do_dir
    sdo_rec = do_rec
    sdo_adj = do_adj
    sdo_self = do_self
    sdo_bond = do_bond
    sdo_angle = do_angle
    sdo_ephi = do_ephi
    szerovdw = zerovdw
    szerochg = zerochg
    szerodip = zerodip
    return
end subroutine save_flow
!------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine restore_flow here]
subroutine restore_flow(sdo_dir, sdo_rec, sdo_adj, sdo_self, &
    sdo_bond, sdo_angle, sdo_ephi, szerochg, szerodip, szerovdw)
    implicit none
#  include "flocntrl.h"
#  include "debug.h"
    integer sdo_dir, sdo_rec, sdo_adj, sdo_self, &
        sdo_bond, sdo_angle, sdo_ephi
    integer szerochg, szerovdw, szerodip
    do_dir = sdo_dir
    do_rec = sdo_rec
    do_adj = sdo_adj
    do_self = sdo_self
    do_bond = sdo_bond
    do_angle = sdo_angle
    do_ephi = sdo_ephi
    zerovdw = szerovdw
    zerochg = szerochg
    zerodip = szerodip
    return
end subroutine restore_flow
!------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine zero_flow here]
subroutine zero_flow()
    implicit none
#  include "flocntrl.h"
#  include "debug.h"
    do_dir = 0
    do_rec = 0
    do_adj = 0
    do_self = 0
    do_bond = 0
    do_angle = 0
    do_ephi = 0
    return
end subroutine zero_flow
!------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dump_flow here]
subroutine dump_flow()
    implicit none
#  include "flocntrl.h"
#  include "debug.h"
    write (6, 60) do_bond, do_angle, do_ephi, do_dir, do_rec, &
        do_adj, do_self, zerovdw, zerochg
60  format(1x, 9i5)
    return
end subroutine dump_flow
!------------------------------------------------------------
