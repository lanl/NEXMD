! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_timer_array here]
subroutine fill_timer_array()
    implicit none
#  include "def_time.h"

    ! after adding the defined integer to "def_time.h"
    ! add a call to add_timer
    ! the arguments are the index of the timer, the
    ! index of its parent timer (must have been added before
    ! unless its -1) and an output string <= 20 characters

    call add_timer(TIME_TOTAL, -1, 'Total time')
    call add_timer(TIME_RDPARM, TIME_TOTAL, 'Read topology time')
    call add_timer(TIME_RDCRD, TIME_TOTAL, 'Read coords time')
    call add_timer(TIME_FASTWT, TIME_TOTAL, 'Fast Water setup')
    call add_timer(TIME_RUNMD, TIME_TOTAL, 'Runmd Time')
    call add_timer(TIME_FORCE, TIME_RUNMD, 'Force time')
    call add_timer(TIME_NONBON, TIME_FORCE, 'Nonbond force')
    !-- list
    call add_timer(TIME_LIST, TIME_NONBON, 'List time')
    call add_timer(TIME_BLDLST, TIME_LIST, 'Build the list')
    !--ewald
    call add_timer(TIME_EWALD, TIME_NONBON, 'Ewald time')
    call add_timer(TIME_DIR, TIME_EWALD, 'Direct Ewald time')
    call add_timer(TIME_SHORT_ENE, TIME_DIR, 'Short_ene time')
    ! vdw
    call add_timer(TIME_VDW, TIME_DIR, 'VDW time')
    call add_timer(TIME_ADJ, TIME_EWALD, 'Adjust Ewald time')
    call add_timer(TIME_SELF, TIME_EWALD, 'Self Ewald time')
    call add_timer(TIME_REC, TIME_EWALD, 'Recip Ewald time')
    call add_timer(TIME_EWFRCADJ, TIME_EWALD, "Force Adjust")
    call add_timer(TIME_EWVIRIAL, TIME_EWALD, "Virial junk")
    call add_timer(TIME_EWFSTRT, TIME_EWALD, "Start synchronization")
    call add_timer(TIME_BSPL, TIME_REC, 'Fill Bspline coeffs')
    call add_timer(TIME_FILLG, TIME_REC, 'Fill charge grid')
    call add_timer(TIME_SCSUM, TIME_REC, 'Scalar sum')
    call add_timer(TIME_GRADS, TIME_REC, 'Grad sum')
    call add_timer(TIME_FFT, TIME_REC, 'FFT time')

    call add_timer(TIME_FFTCOMM, &
        TIME_FFT, 'FFT back comm time')
    call add_timer(TIME_FFTTRANS1, &
        TIME_FFTCOMM, 'FFT back trnsp time')
    call add_timer(TIME_FFTXTRA1, &
        TIME_FFTCOMM, 'FFT back waitall time')
    call add_timer(TIME_FFTXTRA2, &
        TIME_FFTCOMM, 'FFT call transp')

    call add_timer(TIME_FFTCOMM3, &
        TIME_FFT, 'FFT fwd comm')
    call add_timer(TIME_FFTTRANS3, &
        TIME_FFTCOMM3, 'FFT fwd trnsp')
    call add_timer(TIME_FFTXTRA3, &
        TIME_FFTCOMM3, 'FFT fwd waitall')
    call add_timer(TIME_FFTXTRA4, &
        TIME_FFTCOMM3, 'FFT fwd waitany')
    call add_timer(TIME_FFTXTRA5, &
        TIME_FFTCOMM3, 'FFT fwd postrcv')
    call add_timer(TIME_FFTXTRA6, &
        TIME_FFTCOMM3, 'FFT fwd send')
    call add_timer(TIME_FFTCOMM4, &
        TIME_FFTCOMM3, 'FFT packsend')

    call add_timer(TIME_EWFACTORS, TIME_REC, 'Calc trig tables')
    call add_timer(TIME_EWRECSUM, &
        TIME_REC, 'Regular Ewald sum')
    call add_timer(TIME_DISTDIP, &
        TIME_EWALD, 'dipole distribute time')
    call add_timer(TIME_COLLFIELD, &
        TIME_EWALD, 'Field Collect time')
    call add_timer(TIME_LESADJ, TIME_EWALD, 'LES adjust time')

#ifdef RISMSANDER
    !--3d-rism
    call add_timer(TIME_RISM, TIME_NONBON, '3D-RISM time')
    call add_timer(TIME_ULJUV, TIME_RISM, 'LJ Grid time')
    call add_timer(TIME_UCOULU, TIME_RISM, 'Ewald Grid time')
    call add_timer(TIME_ASYMP, TIME_RISM, 'Asymptotics time')
    call add_timer(TIME_RXRISM, TIME_RISM, 'RXRISM time')
    call add_timer(TIME_R1RISM, TIME_RXRISM, 'R1RISM time')
    call add_timer(TIME_RISMFFT, TIME_R1RISM, 'FFT time')
    call add_timer(TIME_MDIIS, TIME_R1RISM, 'MDIIS time')
    call add_timer(TIME_MDIIS_LAPACK, TIME_MDIIS, 'LAPACK time')
    call add_timer(TIME_MDIIS_DATA, TIME_MDIIS, 'DATA time')
    call add_timer(TIME_EXCHEM, TIME_RISM, 'EXCHEM time')
    call add_timer(TIME_FF, TIME_RISM, 'FF time')
    call add_timer(TIME_SAVECRDINTERP, TIME_RISM, 'Save Interpolation Data time')
    call add_timer(TIME_CRDINTERP, TIME_RISM, 'CRD Interpolation time')
    call add_timer(TIME_CRDCUTLIST, TIME_CRDINTERP, 'Interpolation Neighbour List')
    call add_timer(TIME_REORIENT, TIME_RISM, 'Reorient Solute time')
    call add_timer(TIME_RESIZE, TIME_RISM, 'Resize Solvent Box time')
    call add_timer(TIME_PMV, TIME_RISM, 'Partial Molar Volume time')
    call add_timer(TIME_CUVPROP, TIME_RISM, 'Solution Propagation time')
#endif

    !-- egb
    call add_timer(TIME_EGB, TIME_NONBON, 'Gen Born time')
    call add_timer(TIME_GBRAD1, TIME_EGB, 'Calc gb radii')
    call add_timer(TIME_GBRADDIST, &
        TIME_EGB, 'Communicate gb radii')
    call add_timer(TIME_GBRAD2, TIME_EGB, 'Calc gb diag')
    call add_timer(TIME_GBFRC, TIME_EGB, 'Calc gb off-diag')
    call add_timer(TIME_GBSA, TIME_EGB, 'Surface area energy')
    !-- pb
    call add_timer(TIME_PBFORCE, TIME_NONBON, 'PB Nonbond')
    call add_timer(TIME_PBLIST, TIME_PBFORCE, 'PB NB list')
    call add_timer(TIME_PBSETUP, TIME_PBFORCE, 'PB FD Setup')
    call add_timer(TIME_PBFDFRC, TIME_PBFORCE, 'PB FD Force')
    call add_timer(TIME_PBSOLV, TIME_PBFDFRC, 'PB Solver')
    call add_timer(TIME_PBEPS, TIME_PBFDFRC, 'PB Sas/Eps')
    call add_timer(TIME_PBDIRECT, TIME_PBFORCE, 'PB Direct')
    call add_timer(TIME_PBMP, TIME_PBFORCE, 'PB Multiple')
    !-- np
    call add_timer(TIME_NPFORCE, TIME_NONBON, 'NP Nonbond')
    call add_timer(TIME_NPSAS, TIME_NPFORCE, 'NP Sas')
    call add_timer(TIME_NPCAV, TIME_NPFORCE, 'NP Cavity')
    call add_timer(TIME_NPDIS, TIME_NPFORCE, 'NP Dispersion')
    !-- misc nonbond
    call add_timer(TIME_TRIPDIP, TIME_NONBON, 'Triple dipole time')
    call add_timer(TIME_EEXIPS, TIME_NONBON, 'IPS excludes')
    call add_timer(TIME_YAMMPNB, TIME_NONBON, 'Yammp nonbond time')
    !-- bonded
    call add_timer(TIME_QMMM, TIME_FORCE, 'QMMM')
    call add_timer(TIME_QMMMSETUP, TIME_QMMM, 'QMMM setup')
    call add_timer(TIME_QMMMEWALDSETUP, TIME_QMMMSETUP, 'QMMM ewald setup')
    call add_timer(TIME_QMMMVARIABLESOLVCALC, TIME_QMMM, 'QMMM Var Solv Calc')
    call add_timer(TIME_QMMMEWALDKTABLE, TIME_QMMM, 'QMMM Ewald KTable')
    call add_timer(TIME_QMMMLISTBUILD, TIME_QMMM, 'QMMM list build')
    call add_timer(TIME_QMMMCOORDSX, TIME_QMMM, 'QMMM prep coords')
    call add_timer(TIME_QMMMRIJEQNS, TIME_QMMM, 'QMMM RIJ Eqns Calc')
    call add_timer(TIME_QMMMENERGY, TIME_QMMM, 'QMMM energy')
    call add_timer(TIME_QMMMENERGYHCORE, TIME_QMMMENERGY, 'QMMM hcore calc')
    call add_timer(TIME_QMMMENERGYHCOREQM, TIME_QMMMENERGYHCORE, 'QMMM hcore QM-QM')
    call add_timer(TIME_QMMMENERGYHCOREQMMM, TIME_QMMMENERGYHCORE, 'QMMM hcore QM-MM')
    call add_timer(TIME_QMMMENERGYSCF, TIME_QMMMENERGY, 'QMMM scf')
    call add_timer(TIME_QMMMENERGYSCFDENPRED, TIME_QMMMENERGYSCF, 'QMMM Density Predict')
    call add_timer(TIME_QMMMENERGYSCFFOCKPRED, TIME_QMMMENERGYSCF, 'QMMM Fock Predict')
    call add_timer(TIME_QMMMENERGYSCFFOCK, TIME_QMMMENERGYSCF, 'QMMM fock build')
    call add_timer(TIME_QMMMENERGYSCFFOCKRED, TIME_QMMMENERGYSCF, 'QMMM fock dist')
    call add_timer(TIME_QMMMENERGYSCFFOCKEWALD, TIME_QMMMENERGYSCFFOCK, 'QMMM Ewald Contrib')
    call add_timer(TIME_QMMMENERGYSCFFOCKGB, TIME_QMMMENERGYSCFFOCK, 'QMMM GB Fock Terms')
    call add_timer(TIME_QMMMENERGYSCFELEC, TIME_QMMMENERGYSCF, 'QMMM elec-energy calc')
    call add_timer(TIME_QMMMENERGYSCFDIAG, TIME_QMMMENERGYSCF, 'QMMM full matrix diag')
    call add_timer(TIME_QMMMENERGYSCFPSEUDO, TIME_QMMMENERGYSCF, 'QMMM pseudo matrix diag')
    call add_timer(TIME_QMMMENERGYSCFDEN, TIME_QMMMENERGYSCF, 'QMMM density build')
    call add_timer(TIME_QMMMENERGYSCFDENBCAST, TIME_QMMMENERGYSCF, 'QMMM density dist')
    call add_timer(TIME_QMMMEWALDENERGY, TIME_QMMMENERGY, 'QMMM ewald energy')
    call add_timer(TIME_QMMMGBENERGY, TIME_QMMMENERGY, 'QMMM GB energy')
    call add_timer(TIME_QMMMDFTBDISPE, TIME_QMMMENERGY, 'QMMM DFTB Disp E')
    call add_timer(TIME_QMMMFQM, TIME_QMMM, 'QMMM QM-QM force')
    call add_timer(TIME_QMMMDFTBDISPF, TIME_QMMMFQM, 'QMMM DFTB Disp Grad')
    call add_timer(TIME_QMMMDFTBREPULF, TIME_QMMMFQM, 'QMMM DFTB Repul Grad')
    call add_timer(TIME_QMMMDFTBHZEROF, TIME_QMMMFQM, 'QMMM DFTB Hzero Grad')
    call add_timer(TIME_QMMMDFTBGAMMAF, TIME_QMMMFQM, 'QMMM DFTB Gamma Grad')
    call add_timer(TIME_QMMMFQMMM, TIME_QMMM, 'QMMM QM-MM force')
    call add_timer(TIME_QMMMFQMEWALD, TIME_QMMM, 'QMMM Ewald force')
    call add_timer(TIME_QMMMMULLIK, TIME_QMMM, 'QMMM Mulliken Chgs')
    call add_timer(TIME_QMMMCOLLATEF, TIME_QMMM, 'QMMM Collate Forces')
    call add_timer(TIME_BOND, TIME_FORCE, 'Bond/Angle/Dihedral')
    !-- force communications
    call add_timer(TIME_COLLFRC, TIME_FORCE, 'FRC Collect time')
    !-- more runmd stuff
    call add_timer(TIME_SHAKE, TIME_RUNMD, 'Shake time')
    call add_timer(TIME_VERLET, TIME_RUNMD, 'Verlet update time')
    call add_timer(TIME_DIPUP, TIME_RUNMD, 'Dipole update time')
    call add_timer(TIME_EKCMR, TIME_RUNMD, 'Ekcmr time')
    !-- Coord communications
    call add_timer(TIME_DISTCRD, TIME_RUNMD, 'CRD distribute time')

    !-- NOE
    call add_timer(TIME_NOE, TIME_FORCE, 'Noe calc time ')
    call add_timer(TIME_CALRATE, TIME_NOE, 'Calrate time')
    call add_timer(TIME_DSPEV, TIME_NOE, 'Dspev time')
    call add_timer(TIME_SHFDER, TIME_NOE, 'shift der time')
    call add_timer(TIME_KMAT, TIME_NOE, 'Kmat time')
    call add_timer(TIME_NOECALC1, TIME_NOE, 'Noecalc1 time')
    call add_timer(TIME_NOECALC2, TIME_NOE, 'Noecalc2 time')
    call add_timer(TIME_CALDIS, TIME_NOE, 'Calc dis time')
    call add_timer(TIME_REMARC, TIME_NOE, 'Remarc time')
    call add_timer(TIME_DINTEN, TIME_NOE, 'Dinten time')
    call add_timer(TIME_RINGCURR, TIME_NOE, 'Ringcurr time')
    call add_timer(TIME_ELECNOE, TIME_NOE, 'Electro. noe time')
    call add_timer(TIME_ANISO, TIME_NOE, 'Anisotr. noe time')
    call add_timer(TIME_DRATES, TIME_NOE, 'Anisotr. noe time')
    !--IPS
    call add_timer(TIME_AIPS, TIME_NONBON, 'AIPS time')
    call add_timer(TIME_AIPS_FUNC, TIME_AIPS, 'AIPS function')
    call add_timer(TIME_AIPS_GRID, TIME_AIPS, 'AIPS grid')
    call add_timer(TIME_AIPS_FFT, TIME_AIPS, 'AIPS FFT')
    call add_timer(TIME_AIPS_SUM, TIME_AIPS, 'AIPS sum')
    call add_timer(TIME_AIPS_FRC, TIME_AIPS, 'AIPS force')

    return
end subroutine fill_timer_array
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine init_timers here]
subroutine init_timers()
    implicit none
#  include "new_time.h"
#  include "def_time.h"
    integer i, success
    do i = 1, maxtime
        tpar_p(i) = 0
        tchild_p(i) = 0
        tsib_p(i) = 0
        t_level(i) = 0
        t_state(i) = 0
        t_added(i) = 0
        tnext_p(i) = 0
        tprev_p(i) = 0
    end do
    do i = 1, maxtime
        t_accum(i) = 0.d0
        tch_acc(i) = 0.d0
        t_curr(i) = 0.d0
    end do
    do i = 1, maxtime
        t_string(i) = ' '
    end do
    call fill_timer_array()

    ! establish what level each timer is at

    do i = 1, 10000
        call get_level(success)
        if (success == 1) goto 100
    end do
100 continue

    ! build the lists(s) for output
    ! should imitate a depth-first tree walk

    call build_lists()
    return
end subroutine init_timers
!----------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine add_timer here]
subroutine add_timer(index, par_index, string)
    implicit none
#  include "new_time.h"
    integer index, par_index
    character(len=*) string
    integer i

    if (index > maxtime) then
        write (6, *) 'index ', index, ' bigger than MAXTIME ', maxtime
        write (6, *) 'attempt to add timer ', string
        call mexit(6, 1)
    end if
    if (t_added(index) /= 0) then
        write (6, 8) index, string
8       format(1x, 'add_timer: timer (', i2, 1x, a, ') already added! ')
        call mexit(6, 1)
    end if
    if (par_index > 0) then
        if (t_added(par_index) == 0) then
            write (6, *) 'add_timer: parent timer not yet added! '
            write (6, 9) index, par_index, string
9           format(1x, 'index,par_index,string = ', i2, 1x, i2, 1x, a)
            call mexit(6, 1)
        end if
    end if
    t_added(index) = 1
    t_string(index) = string
    tpar_p(index) = par_index
    if (par_index > 0) then

        ! fortran linked list of sibs

        if (tchild_p(par_index) == 0) then
            tchild_p(par_index) = index
        else
            i = tchild_p(par_index)
10          continue
            if (tsib_p(i) > 0) then
                i = tsib_p(i)
                goto 10
            else
                tsib_p(i) = index
            end if
        end if
    end if
    return
end subroutine add_timer
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Synchronize processors to ensure timer accuracy.
subroutine timer_barrier(communicator)

    implicit none
    integer communicator
    logical barrier_active

#ifdef MPI
    include 'mpif.h'
    integer ierr
#endif

    barrier_active = .true.
    if (barrier_active) then
#ifdef MPI
        call mpi_barrier(communicator, ierr)
#endif
    end if
    return
end subroutine timer_barrier
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine timer_start here]
subroutine timer_start(istart)
    implicit none
#  include "new_time.h"
    _REAL_ tim
    integer istart
    call wallclock(tim)
    if (istart > 0 .and. istart < maxtime) then
        if (t_state(istart) /= 0) then
            write (6, *) 'Invalid start to time: ', &
                istart, ' ', t_string(istart)
            call mexit(6, 1)
        end if
        t_curr(istart) = tim
        t_state(istart) = 1
    end if
    return
end subroutine timer_start
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine timer_stop here]
subroutine timer_stop(istop)
    implicit none
#  include "new_time.h"
    _REAL_ tim
    integer istop
    call wallclock(tim)
    if (istop > 0 .and. istop < maxtime) then
        if (t_state(istop) == 0) then
            write (6, *) 'Invalid stop to time: ', &
                istop, ' ', t_string(istop)
            call mexit(6, 1)
        end if
        t_accum(istop) = t_accum(istop) + tim - t_curr(istop)
        ! reset, so can't stop again without start
        t_state(istop) = 0
    end if
    return
end subroutine timer_stop
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine timer_stop_start here]
subroutine timer_stop_start(istop, istart)
    implicit none
#  include "new_time.h"
    _REAL_ tim
    integer istop, istart

    ! TEC suggestion... combine above 2 to save calls

    call wallclock(tim)
    if (istop > 0 .and. istop < maxtime) then
        if (t_state(istop) == 0) then
            write (6, *) 'Invalid stop to time: ', t_string(istop)
            call mexit(6, 1)
        end if
        t_accum(istop) = t_accum(istop) + tim - t_curr(istop)

        ! reset, so can't stop again without start

        t_state(istop) = 0
    end if
    if (istart > 0 .and. istart < maxtime) then
        if (t_state(istart) /= 0) then
            write (6, *) 'Invalid start to time: ', t_string(istart)
            call mexit(6, 1)
        end if
        t_curr(istart) = tim
        t_state(istart) = 1
    end if
    return
end subroutine timer_stop_start
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_level here]
subroutine get_level(success)
    implicit none
    integer success
#  include "new_time.h"
    integer i
    do i = 1, maxtime

        ! only look at undetermined ones

        if (t_added(i) == 1) then
            if (t_level(i) == 0) then
                if (tpar_p(i) == -1) then
                    t_level(i) = 1
                else if (t_level(tpar_p(i)) /= 0) then
                    t_level(i) = t_level(tpar_p(i)) + 1
                end if
            end if
        end if
    end do

    ! check if all defined

    success = 1
    do i = 1, maxtime
        if (t_added(i) == 1 .and. t_level(i) == 0) success = 0
    end do
    return
end subroutine get_level
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine build_lists here]
subroutine build_lists()
    implicit none
#  include "new_time.h"
    integer i, j, numfirst, numlast
    do i = 1, maxtime
        if (t_added(i) == 1) then
            if (tsib_p(i) > 0) then
                if (tchild_p(tsib_p(i)) > 0) then

                    ! follow the children to get Tnext_P(i)

                    j = tchild_p(tsib_p(i))
5                   continue
                    if (tchild_p(j) /= 0) then
                        j = tchild_p(j)
                        goto 5
                    end if
                    tnext_p(i) = j
                else
                    tnext_p(i) = tsib_p(i)
                end if
            else if (tpar_p(i) > 0) then
                tnext_p(i) = tpar_p(i)
            end if
            if (tnext_p(i) /= 0) then
                tprev_p(tnext_p(i)) = i
            end if
        end if
    end do
    t_numtree = 0
    do i = 1, maxtime
        if (t_added(i) == 1) then
            if (tnext_p(i) == 0) then
                t_numtree = t_numtree + 1
                if (t_numtree > t_maxtree) then
                    write (6, *) 'Too many timer trees. Check T_MAXTREE'
                    call mexit(6, 1)
                end if
                t_last(t_numtree) = i
            end if
        end if
    end do

    ! follow ends back to get start of lists

    do j = 1, t_numtree
        i = t_last(j)
10      continue
        if (tprev_p(i) /= 0) then
            i = tprev_p(i)
            goto 10
        end if
        t_first(j) = i
    end do
    return
end subroutine build_lists
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine finish_timers here]
subroutine finish_timers(all)
    implicit none
    _REAL_ all
#  include "new_time.h"
    _REAL_ t1, p1, p2, p3, minp
    integer i, j, k, max

    ! get the max level

    max = 0
    do i = 1, maxtime
        if (t_level(i) > max) max = t_level(i)
    end do

    ! accumulate child times
    ! do from maxlevel up to level 2 (level 1 has no parent)
    ! first zero them

    do i = 1, maxtime
        tch_acc(i) = 0.d0
    end do
    do j = 1, max - 1
        k = max - j + 1
        do i = 1, maxtime
            if (t_level(i) == k) then

                ! make sure accumulated time is at least child times (roundoff)

                if (tch_acc(i) > t_accum(i)) t_accum(i) = tch_acc(i)
                if (tpar_p(i) > 0) then
                    tch_acc(tpar_p(i)) = tch_acc(tpar_p(i)) + t_accum(i)
                end if
            end if
        end do
    end do

    ! get other times

    do i = 1, maxtime
        t_other(i) = t_accum(i) - tch_acc(i)
    end do

    ! minimal percent of parent to be printed

    minp = 0.005d0

    ! T_print = 0 means don't print
    ! T_print = 1 means print
    ! T_print = 2 means there is time unaccounted for by child times
    !             that needs to be printed as well
    ! find who gets printed. Need to go down by levels
    ! you don't get printed unless your parent is

    do i = 1, maxtime
        t_print(i) = 0
    end do
    do j = 1, max
        do i = 1, maxtime
            if (t_level(i) == j) then
                if (tpar_p(i) > 0) then
                    t1 = t_accum(tpar_p(i))
                else
                    t1 = all
                end if
                if (t1 > 1.d-4) then
                    p1 = 100.d0*t_accum(i)/t1
                else
                    p1 = 0.d0
                end if
                if (t_accum(i) > 1.d-4) then
                    p2 = 100.d0*t_other(i)/t_accum(i)
                    p3 = 100.d0*tch_acc(i)/t_accum(i)
                else
                    p2 = 0.d0
                    p3 = 0.d0
                end if
                if (tpar_p(i) > 0) then

                    ! is your parent printable and are you a big enough percent of it??

                    if (t_print(tpar_p(i)) > 0 .and. p1 > minp) then

                        ! do you have children and is there some time not in child times??

                        if (tchild_p(i) > 0 .and. p2 > minp &
                            .and. p3 > minp) then
                            t_print(i) = 2
                        else
                            t_print(i) = 1
                        end if
                    end if
                else
                    if (p1 > minp) then
                        if (tchild_p(i) > 0 .and. p2 > minp) then
                            t_print(i) = 2
                        else
                            t_print(i) = 1
                        end if
                    end if
                end if
            end if
        end do
    end do

    return
end subroutine finish_timers
!---------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine write_timers here]
subroutine write_timers(all, nf)
    implicit none
    _REAL_ all
    integer nf
#  include "new_time.h"
#  include "def_time.h"
    integer i, j, k
    integer length(10)
    _REAL_ t, p
    _REAL_ total
    character(len=35) space
    character(len=5) short
    character(len=20) other
    other = 'Other'
    space = '                         '
    length(1) = 1
    do j = 2, 10
        length(j) = length(j - 1) + 3
    end do
    do k = 1, t_numtree
        i = t_first(k)
10      continue
        if (t_print(i) > 0) then
            if (t_print(i) == 2) then

                ! need to first print leftover time from sum of children

                short = t_string(i) (1:5)
                t = t_other(i)
                p = 100.d0*t_other(i)/t_accum(i)
                j = t_level(i) + 1
                if (p > 99.9) then
                    write (nf, 29) space(1:length(j)), other, t, short
                else
                    write (nf, 30) space(1:length(j)), other, t, p, short
                end if
            end if

            ! now print out entry for this time

            if (tpar_p(i) > 0) then
                short = t_string(tpar_p(i)) (1:5)
                j = tpar_p(i)
                total = t_accum(j)
            else
                short = 'ALL '
                total = all
            end if
            t = t_accum(i)
            p = 100.d0*t/total
            j = t_level(i)
            if (j > 10) then
                write (nf, *) 'level too big:', t_string(i)
            end if
            if (p > 99.9) then
                write (nf, 29) space(1:length(j)), t_string(i), t, short
            else
                write (nf, 30) space(1:length(j)), t_string(i), t, p, short
            end if
        end if
        if (tnext_p(i) /= 0) then
            i = tnext_p(i)
            goto 10
        end if
    end do
29  format('|', a, a, 1x, f10.2, ' (100.0% of ', a, ')')
30  format('|', a, a, 1x, f10.2, ' (', f5.2, '% of ', a, ')')
    return
end subroutine write_timers
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine write_timerstats here]
subroutine write_timerstats(tim, tim2, mint, maxt, &
    to, to2, minto, maxto, nf)
    implicit none
    _REAL_ all
    _REAL_ tim(*), tim2(*), mint(*), maxt(*), &
        to(*), to2(*), minto(*), maxto(*)
    integer nf
#  include "new_time.h"
#  include "def_time.h"
    integer i, j, k
    integer length(10)
    _REAL_ total
    character(len=35) space
    character(len=20) other
    other = 'Other'
    space = '                         '
    length(1) = 1
    do j = 2, 10
        length(j) = length(j - 1) + 3
    end do
    do k = 1, t_numtree
        i = t_first(k)
10      continue
        if (t_print(i) > 0) then
            if (t_print(i) == 2) then
                ! need to first print leftover time from sum of children
                j = t_level(i) + 1
                write (nf, 30) space(1:length(j)), other, to(i), &
                    minto(i), maxto(i), to2(i)
            end if
            ! now print out entry for this time
            j = t_level(i)
            if (j > 10) then
                write (nf, *) 'level too big:', t_string(i)
            end if
            write (nf, 30) space(1:length(j)), t_string(i), &
                tim(i), mint(i), maxt(i), tim2(i)
        end if
        if (tnext_p(i) /= 0) then
            i = tnext_p(i)
            goto 10
        end if
    end do
30  format('|', a, a, 1x, f10.2, ' (', 3(f10.2), ')')
    return
end subroutine write_timerstats
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine profile_time here]
subroutine profile_time(all, num_calls_nblist, profile_mpi)
    use stack, only : ihighest_stk, highest_stk
    implicit none
    _REAL_ all
    integer, intent(in) :: num_calls_nblist, profile_mpi
#  include "new_time.h"
#  include "extra.h"

#ifdef MPI
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
    include 'mpif.h'
    integer ierr
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
    integer ist(mpi_status_size), ier
#endif
    _REAL_ buf(maxtime + 1), tim(maxtime), &
        to(maxtime), to2(maxtime), minto(maxtime), &
        maxto(maxtime), &
        tim2(maxtime), mint(maxtime), maxt(maxtime), &
        allave
    integer i, j
    if (master) write (6, '(/80("-")/,"   5.  TIMINGS",/80("-")/)')
#ifndef NO_DETAILED_TIMINGS

# ifndef MPI
    ! single processor code
    call finish_timers(all)
    call write_timers(all, 6)
    write (6, '(/,a,i10)') '| Highest rstack allocated: ', highest_stk
    write (6, '(a,i10)') '| Highest istack allocated: ', ihighest_stk
    return
# else
    if (numtasks == 1 .or. mpi_orig) then
        call finish_timers(all)
        call write_timers(all, 6)
        write (6, '(/,a,i10)') '| Highest rstack allocated: ', highest_stk
        write (6, '(a,i10)') '| Highest istack allocated: ', ihighest_stk
        return
    end if

    ! numtasks > 1...get processor by processor info, do statistics

    if (.not. master) then
        do i = 1, maxtime
            buf(i) = t_accum(i)
        end do
        buf(maxtime + 1) = all

        ! send to master

        call mpi_send(buf, maxtime + 1, MPI_DOUBLE_PRECISION, &
            0, mytaskid, commsander, ier)
        return
    else

        ! go through processors. first yourself
        ! only write to profile_mpi if profile_mpi=1 in cntrl namelist
        ! note this can really slow down multisander jobs.
        if (profile_mpi == 1) then
            open (unit=8, file='profile_mpi')
            j = 0
            write (8, 16) j
            call finish_timers(all)
            call write_timers(all, 8)
        else
            call finish_timers(all) !RCW not sure if we need this.
        end if

        ! update for average,stats

        do i = 1, maxtime
            tim(i) = t_accum(i)
            tim2(i) = t_accum(i)**2
            mint(i) = t_accum(i)
            maxt(i) = t_accum(i)
            to(i) = t_other(i)
            to2(i) = t_other(i)**2
            minto(i) = t_other(i)
            maxto(i) = t_other(i)
        end do
        allave = all
        do j = 1, numtasks - 1
            call mpi_recv(buf, maxtime + 1, MPI_DOUBLE_PRECISION, &
                j, j, commsander, ist, ierr)
            do i = 1, maxtime
                t_accum(i) = buf(i)
            end do
            all = buf(maxtime + 1)
            if (profile_mpi == 1) then
                write (8, 16) j
                call finish_timers(all)
                call write_timers(all, 8)
            end if

            ! update for average,stats

            do i = 1, maxtime
                tim(i) = tim(i) + t_accum(i)
                tim2(i) = tim2(i) + t_accum(i)**2
                if (t_accum(i) < mint(i)) mint(i) = t_accum(i)
                if (t_accum(i) > maxt(i)) maxt(i) = t_accum(i)
                to(i) = to(i) + t_other(i)
                to2(i) = to2(i) + t_other(i)**2
                if (t_other(i) < minto(i)) minto(i) = t_other(i)
                if (t_other(i) > maxto(i)) maxto(i) = t_other(i)
            end do
            allave = allave + all
        end do
        do i = 1, maxtime
            t_accum(i) = tim(i)/numtasks
        end do
        all = allave/numtasks

        ! the averages are run through finish_timers.. this sets print options

        write (6, 17)
        call finish_timers(all)
        call write_timers(all, 6)

        if (profile_mpi == 1) then
            ! do the statistics

            do i = 1, maxtime
                tim(i) = tim(i)/numtasks
                tim2(i) = sqrt(tim2(i)/numtasks - tim(i)**2)
                to(i) = to(i)/numtasks
                to2(i) = sqrt(to2(i)/numtasks - to(i)**2)
            end do
            write (8, 18)
            write (8, 19)
            call write_timerstats(tim, tim2, mint, maxt, to, to2, minto, maxto, 8)
            close (unit=8)
        end if
        write (6, '(/,a,i10)') '| Number of list builds   : ', num_calls_nblist
        write (6, '(/,a,i10)') '| Highest rstack allocated: ', highest_stk
        write (6, '(a,i10)') '| Highest istack allocated: ', ihighest_stk
    end if  ! ( .not. master)
    return
# endif

#endif

16  format('|>>>>>>>>PROFILE of TIMES  for process ', i4)
17  format('|>>>>>>>>PROFILE of Average TIMES>>>>>>>>> ')
18  format('|>>>>>>>>Statistics of TIMES>>>>>>>>> ')
19  format('|>>>>>>>>Printed as average time ', &
        '(min,max,sd) >>>>>>>>> ')
end subroutine profile_time
!------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Prints a summary of the estimated runtime remaining etc to mdinfo
!+ and optionally mdout.
!+ Written by Ross Walker (SDSC, 2009)
subroutine print_ongoing_time_summary(total_steps, current_step, time_step, write_unit)
    use qmmm_module, only : qmmm_nml

! Make the ongoing time summary work with explicit Constant pH simulations JMS 2/2011
    use constantph, only : relaxations, tot_relax

    implicit none

! Passed in
    integer :: total_steps  !Total steps in simulation (nstlim)
    integer :: current_step  !Current step number.
    integer :: write_unit  !Unit number to write the output to.
    _REAL_ :: time_step  !Time step in ps.

! Local
    integer, save :: last_step_count = 0 !Step count the last time we were called.
    integer :: steps_remaining, step_interval, est_steps_remaining ! estimated steps remaining JMS 2/2011
    integer :: full_current_step, full_step_interval ! accounts for solvent relaxation JMS 2/2011
    _REAL_, save :: previous_time
    _REAL_, save :: start_time
    _REAL_ :: current_time
    _REAL_ :: elapsed_time, total_elapsed_time, time_per_step, time_remaining
    _REAL_ :: avg_time_per_step
    _REAL_ :: ns_per_day, sec_per_ns, total_sec_per_ns, total_ns_per_day

! Needed for ntrelax. ntrelax = 0 if icnstph /= 2 in mdread.f
#include "md.h"

    call wallclock(current_time)

    !The following is for initialization only.
    if (current_step == 0) then
        previous_time = current_time
        start_time = current_time
        return
    end if

    steps_remaining = total_steps - current_step

    ! estimate how many steps we have left by estimating how many relaxations we
    ! will do from here on out based on how many have been done already, then
    ! multiplying it by the number of relaxation steps we perform: JMS 2/2011
    est_steps_remaining = tot_relax/current_step*steps_remaining*ntrelax + steps_remaining
    step_interval = current_step - last_step_count

    ! Adjust current_step and step_interval for relaxation steps taken (JMS 2/2011)
    full_current_step = current_step + relaxations*ntrelax
    full_step_interval = step_interval + relaxations*ntrelax

    ! Elapsed time in seconds for last interval and over whole sim
    elapsed_time = current_time - previous_time
    total_elapsed_time = current_time - start_time

    ! Time per step in ms for recent chunk and avg over all completed steps
    time_per_step = 1000.0d0*elapsed_time/dble(full_step_interval)
    avg_time_per_step = 1000.0d0*total_elapsed_time/dble(full_current_step)

    time_remaining = dble(steps_remaining)*total_elapsed_time/dble(current_step)

    !Current elapsed time is time for time_step * step_interval picoseconds
    !so sec_per_ps = elapsed_time / (step_interval * time_step)
    !   sec_per_ns = sec_per_ps / 1000
    sec_per_ns = 1000.0d0*elapsed_time/(dble(step_interval)*time_step)

    !ns_per_day = (60*60*24) / sec_per_ns
    !           = 86400 / sec_per_ns
    ns_per_day = 86400.0d0/sec_per_ns
    total_sec_per_ns = 1000.0d0*total_elapsed_time/(dble(current_step)*time_step)
    total_ns_per_day = 86400.0d0/total_sec_per_ns

    if (current_step /= total_steps) then
        !Print ongoing timing info.
        write (write_unit, 8088)
        write (write_unit, '(a)') &
            '| Current Timing Info'
        write (write_unit, '(a)') &
            '| -------------------'
        write (write_unit, '(a,i9,a,i9,a,i9)') '| Total steps : ', total_steps, ' | Completed : ', current_step, &
            ' | Remaining : ', steps_remaining
        write (write_unit, '(a)') '|'
        write (write_unit, '(a,i7,a)') '| Average timings for last ', step_interval, ' steps:'
        if (icnstph .gt. 1) then
            write (write_unit, '(a)') '| Solvent relaxation steps only included in per-step timings'
            write (write_unit, '(a,i7,a)') '| Solvent relaxation timings for last ', &
                step_interval, ' steps:'
            write (write_unit, '(a,i9)') '|     Steps per relaxation  = ', ntrelax
            write (write_unit, '(a,i9)') '|     Relaxation cycles     = ', relaxations
            write (write_unit, '(a,i9)') '|     Relaxation time steps = ', ntrelax*relaxations
        end if
        write (write_unit, '(a,f10.2,a,f10.2)') &
            '|     Elapsed(s) = ', elapsed_time, ' Per Step(ms) = ', time_per_step
        if (qmmm_nml%ifqnt .and. qmmm_nml%qmtheory%EXTERN) then
            write (write_unit, '(a,f10.3,a,f10.2)') &
                '|         ps/day = ', ns_per_day*1000.0d0, '   seconds/ps = ', sec_per_ns/1000.0d0
        else
            write (write_unit, '(a,f10.2,a,f10.2)') &
                '|         ns/day = ', ns_per_day, '   seconds/ns = ', sec_per_ns
        end if
        write (write_unit, '(a)') '|'
        write (write_unit, '(a)') '| Average timings for all steps:'
        if (icnstph .gt. 1) then
            write (write_unit, '(a,i7,a)') '| Solvent relaxation timings for last ', &
                step_interval, ' steps:'
            write (write_unit, '(a,i9)') '|     Steps per relaxation  = ', ntrelax
            write (write_unit, '(a,i9)') '|     Relaxation cycles     = ', tot_relax
            write (write_unit, '(a,i9)') '|     Relaxation time steps = ', ntrelax*tot_relax
        end if
        write (write_unit, '(a,f10.2,a,f10.2)') &
            '|     Elapsed(s) = ', total_elapsed_time, ' Per Step(ms) = ', avg_time_per_step
        if (qmmm_nml%ifqnt .and. qmmm_nml%qmtheory%EXTERN) then
            write (write_unit, '(a,f10.3,a,f10.2)') &
                '|         ps/day = ', total_ns_per_day*1000.0d0, '   seconds/ps = ', total_sec_per_ns/1000.0d0
        else
            write (write_unit, '(a,f10.2,a,f10.2)') &
                '|         ns/day = ', total_ns_per_day, '   seconds/ns = ', total_sec_per_ns
        end if
        write (write_unit, '(a)') '|'
        if (time_remaining < 60) then
            write (write_unit, '(a,f9.1,a)') '| Estimated time remaining: ', time_remaining, ' seconds.'
        else if (time_remaining < 3600) then
            write (write_unit, '(a,f9.1,a)') '| Estimated time remaining: ', time_remaining/60.0d0, ' minutes.'
        else
            write (write_unit, '(a,f9.1,a)') '| Estimated time remaining: ', time_remaining/3600.0d0, ' hours.'
        end if
        write (write_unit, 8088)
    else
        !Write final timing info
        write (write_unit, '(/, a)') &
            '| Final Performance Info:'
        write (write_unit, '(a)') '| -----------------------------------------------------'
        if (write_unit /= 6) then
            write (write_unit, '(a,i7,a)') '| Average timings for last ', step_interval, ' steps:'
            if (icnstph .gt. 1) then
                write (write_unit, '(a)') '| Solvent relaxation steps only included in per-step timings'
                write (write_unit, '(a,i7,a)') '| Solvent relaxation timings for last ', &
                    step_interval, ' steps:'
                write (write_unit, '(a,i9)') '|     Steps per relaxation  = ', ntrelax
                write (write_unit, '(a,i9)') '|     Relaxation cycles     = ', relaxations
                write (write_unit, '(a,i9)') '|     Relaxation time steps = ', ntrelax*relaxations
            end if
            write (write_unit, '(a,f10.2,a,f10.2)') &
                '|     Elapsed(s) = ', elapsed_time, ' Per Step(ms) = ', time_per_step
            if (qmmm_nml%ifqnt .and. qmmm_nml%qmtheory%EXTERN) then
                write (write_unit, '(a,f10.3,a,f10.2)') &
                    '|         ps/day = ', ns_per_day*1000.0d0, '   seconds/ps = ', sec_per_ns/1000.0d0
            else
                write (write_unit, '(a,f10.2,a,f10.2)') &
                    '|         ns/day = ', ns_per_day, '   seconds/ns = ', sec_per_ns
            end if
            write (write_unit, '(a)') '|'
        end if
        write (write_unit, '(a)') '| Average timings for all steps:'
        write (write_unit, '(a,f10.2,a,f10.2)') &
            '|     Elapsed(s) = ', total_elapsed_time, ' Per Step(ms) = ', avg_time_per_step
        if (qmmm_nml%ifqnt .and. qmmm_nml%qmtheory%EXTERN) then
            write (write_unit, '(a,f10.3,a,f10.2)') &
                '|         ps/day = ', total_ns_per_day*1000.0d0, '   seconds/ns = ', total_sec_per_ns/1000.0d0
        else
            write (write_unit, '(a,f10.2,a,f10.2)') &
                '|         ns/day = ', total_ns_per_day, '   seconds/ns = ', total_sec_per_ns
        end if
        write (write_unit, '(a, /)') '| -----------------------------------------------------'
    end if

    last_step_count = current_step
    previous_time = current_time
    relaxations = 0 ! JMS 2/2011 -- make timing facility work with explicit constant pH simulations

    return

8088 format(t2, 78('-'))

end subroutine print_ongoing_time_summary
