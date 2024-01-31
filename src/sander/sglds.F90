! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

module sgld

! MODULE: sgld
! ================== SELF-GUIDED MOLECULAR/LANGEVIN DYNAMICS ====================
! Xiongwu Wu, 2011

    implicit none

!     Head file for the self-guided Langevin Dynamics simulation
!
!      variables for SGLD simulation
!

! ... integers:
!
!
!    SGMD/SGLD applying range
!      ISGSTA      Begining atom index applying SGLD
!      ISGEND      Ending atom index applying SGLD
!
!
    integer, save:: ISGSTA, ISGEND, rndfsg, FIXCOM
!
! ... floats:
!
!
!    SGMD/SGLD VARIABLES
!     SGFT    !  Guiding factor
!     TSGAVG  !  Local average time, ps
!     TEMPSG  !  Guiding temperature, K
!     TREFLF  !  Reference low frequency temperature, K
!     TSGAVP  !  Convergence time, ps
!
!
!     SGAVG0  !  Local average remains
!     SGAVG1  !  Local average factor, SGAVG1=1-SGAVG0
!     SGAVP0  !  Convergence average remains
!     SGAVP1  !  Convergency average factor, SGAVP1=1-SGAVP0
!     GAMMAS  !  friction coefficient
!     TEMPSG  !  Target guiding temperature
!     TREFLF  !  Reference low frequency temperature
!

    _REAL_ SGFT, SGFF, SGFD, TSGAVG, TEMPSG, TSGSET, TREFLF, TSGAVP, GAMMAS, SGMASS, &
        SGAVG0, SGAVG1, SGAVP0, SGAVP1, SGFTI, SGFFI, SGFDI, TEMPSGI, EPOTLF, EPOTHF, &
        SGEFLF, SGEFHF, SGCFLF, AVGDD, AVGFF, AVGGD, AVGGF, AVGPP, AVGGP, &
        TEMPLFI, TEMPHFI, TREFLFI, TREFHFI, FRCLF, FRCHF, VIRSG, &
        AVGELF, AVGTLF, AVGEFLF, AVGEFHF, AVGCFLF, AVGCFHF, SGWT, &
        TEMPRXLF, MYSCALSG, SGSCALE

!  common block for parallel broadcast
    common/sgldr/SGFT, SGFF, SGFD, TSGAVG, TEMPSG, TSGSET, TREFLF, &
        SGFTI, SGFFI, SGFDI, TEMPSGI, SGMASS, &
        AVGDD, AVGFF, AVGGD, AVGGF, AVGPP, AVGGP, &
        TEMPLFI, TEMPHFI, TREFLFI, TREFHFI, &
        AVGTLF, AVGEFLF, AVGEFHF, AVGCFLF, AVGCFHF, &
        TEMPRXLF, MYSCALSG
!  Number of broadcasting variables
    integer, parameter :: nsgld_real = 29

!
! ... flags:
!
!
!     TSGLD   ! Perform SGLD
!     TSGLDFP ! Perform SGLDfp to mantain canonical ensemble distribution
!     TLANGV  ! The simulation is a Langevin dynamics simulation
!
!
    logical, save:: TSGLD, TSGLDFP, TLANGV, TRXSGLD

!
! ...allocatable arrays:
!
!     VSG     ! local averages of velocity
!     DSG     ! interaction forces
!     FSG     ! local averages of force
!     GSG     ! guiding force
!     HSG     ! local averages of guiding force
!     SHKSG   ! constraint forces
!
    _REAL_, dimension(:), allocatable, private, save :: VSG, DSG, FSG, GSG, HSG, SHKSG

contains

    subroutine PSGLD(NATOM, AMASS, V)
!-----------------------------------------------------------------------
!     This routine perform initiation for the self-guided
!       Langevin Dynamcs (SGLD) simulaiton
!
        implicit none
#include "md.h"
#ifdef MPI
#  include "parallel.h"
#endif
        integer NATOM, I, I3, M, ierror
        _REAL_ AMASS(*), V(*)
        _REAL_ AMASSI, VI3, GAMM
!
!  Check for invalid sgld setting
        if (ISGSTA < 1) ISGSTA = 1
        if (ISGEND > NATOM .OR. ISGEND < 1) ISGEND = NATOM
        if (TSGAVG .LT. DT) TSGAVG = DT
        if (TSGAVP .LT. DT) TSGAVP = 10.0D0*TSGAVG
        SGAVG1 = DT/TSGAVG
        SGAVG0 = 1.0D0 - SGAVG1
        SGAVP1 = DT/TSGAVP
        SGAVP0 = 1.0D0 - SGAVP1
        if (TLANGV) then
            if (SGFT < -1.0D1) SGFT = 1.0D0
        else
            if (SGFT < -1.0D1) SGFT = 0.2D0
        end if
#ifdef MPI
        if (mytaskid .eq. 0) then
#endif
            write (6, 910) ISGSTA, ISGEND, TSGAVG
            if (TEMPSG > 0.0D0) then
                write (6, 920) TEMPSG
                if (TEMPSG < 100.0D0) write (6, 921)
                if (TSGLDFP) then
                    write (6, 922)
                end if
            else
                write (6, 925) SGFT
                if (SGFF > -1.0D2) write (6, 926) SGFF
            end if
            if (TSGLDFP) then
                write (6, 928)
            end if
            if (TLANGV) then
                write (6, 929) "  SGLD"
                write (6, 930) GAMMA_LN
                if (TSGLDFP) then
                    write (6, 928) "  SGLDfp"
                    if (TREFLF < 1.0D-6) write (6, 941)
                end if
            else
                write (6, 929) "  SGMD"
                if (TSGLDFP) then
                    write (6, 928) "  SGMDfp"
                    if (TREFLF < 1.0D-6) write (6, 951)
                end if
            end if
            if (FIXCOM > 0) write (6, 931)
            write (6, 935)
#ifdef MPI
        end if
#endif
        if (GAMMA_LN > 0.0D0) then
            GAMMAS = GAMMA_LN/20.455d0
        else
            GAMMAS = 1.0D0/20.455d0
        end if
        TSGSET = TEMP0
        if (TSGSET < 1.0D-6) TSGSET = MAX(TEMPSG, 300.0D0)
        GAMM = SQRT(DT/TSGAVG)
        if (TREFLF > 1.0D-6) GAMM = SQRT(TREFLF/TSGSET)
!     Allocate working arrays
        allocate (VSG(3*natom), DSG(3*natom), FSG(3*natom), &
            GSG(3*natom), HSG(3*natom), SHKSG(3*natom), &
            stat=ierror)
        REQUIRE(ierror == 0)
!    Initialize arrays
        SGMASS = 0.0D0
        do I = 1, NATOM
            AMASSI = AMASS(I)
            SGMASS = SGMASS + AMASSI
            if ((I >= ISGSTA) .AND. (I <= ISGEND)) then
            end if
            do M = 1, 3
                I3 = 3*I - 3 + M
                VI3 = V(I3)
                VSG(I3) = GAMM*VI3
                DSG(I3) = VI3
                FSG(I3) = 0.0D0
                GSG(I3) = 0.0D0
                HSG(I3) = 0.0D0
                SHKSG(I3) = 0.0D0
            end do

        end do
        SGFTI = SGFT
        if (SGFF <= -1.0D2) then
            SGFFI = 0.0D0
        else
            SGFFI = SGFF
            write (6, 961) SGFF
        end if
        if (SGFD <= -1.0D2) then
            SGFDI = 0.0D0
        else
            SGFDI = SGFD
            write (6, 962) SGFD
        end if
        TEMPSGI = TEMPSG
        if (TEMPSG < 1.0D-6) TEMPSGI = TSGSET
        TREFLFI = TREFLF
        if (TREFLF < 1.0D-6) TREFLFI = GAMM*TSGSET
        TEMPLFI = TREFLFI
        TREFHFI = TSGSET - TREFLFI
        TEMPHFI = TSGSET - TEMPLFI
        EPOTLF = 2.0D10
        SGWT = 0.0D0
        VIRSG = 0.0D0
        FRCLF = 1.0D0
        FRCHF = 1.0D0
        AVGELF = 0.0D0
        AVGTLF = TREFLFI
        AVGEFLF = 1.0D0
        AVGEFHF = 1.0D0
        AVGCFLF = 1.0D0
        AVGCFHF = 1.0D0
        AVGFF = 1.0D-6
        AVGDD = 1.0D-6
        AVGGF = 0.0D0
        AVGGD = 0.0D0
        AVGPP = 1.0D-6
        AVGGP = 0.0D0
        SGEFLF = 1.0D0
        SGEFHF = 1.0D0
        SGCFLF = 1.0D0
910     format("  _________________ SGLD parameters _________________"/ &
            "  Parameters for self-guided Langevin dynamics (SGLD) simulation"// &
            "  Guiding range from ", I5, "  TO ", I5/ &
            "  Local averaging time: ", F10.4, " ps ")
920     format("  Guiding temperature:", F8.2, " K")
921     format("  *** WARNING: Guiding temperature is redefined since AMBER12 ***"/ &
            "  Set TEMPSG > TEMP0 to enhance conformational search!"/ &
            "  TEMPSG defines a seaching ability comparable to a simulation at T=TEMPSG"/)
922     format("  *** WARNING: SGFT, instead of TEMPSG, should be set for SGMDfp/SGLDfp! "/ &
            "  TEMPSG defines a seaching ability comparable to a simulation at temperature TEMPSG"/)
925     format("  Momentum guiding factor: ", F8.4)
926     format("  Force guiding factor: ", F8.4)
928     format(A8, "  method is used to mantain a canonical distribution. ")
929     format(A8, "  method is used to enhance conformational search. ")
930     format("  Collision frequency:", F8.2, " 1/ps")
931     format("  Translation of COM is freezed!")
935     format("  SGMD/SGLD output properties:"/ &
            "  SGLF=  SGFT   TEMPSG   TEMPLF   TREFLF   FRCLF   EPOTLF    SGWT"/ &
            "  SGHF=  SGFF   SGFD     TEMPHF   TREFHF   FRCHF   EPOTHF   VIRSG"/ &
            "         SGMD/SGLD weighting factor=exp(SGWT)"/ &
            " _______________________________________________________"/)
941     format("  WARNING: TREFLF is not defined and will be estimated from simulation. ")
951     format("  WARNING: TREFLF must be defined for SGMDfp. ")
961     format("  SGFF is fixed at: ", F8.4)
962     format("  SGFD is fixed at: ", F8.4)
        return
    end subroutine PSGLD

    subroutine SGFSHAKE(ISTART, IEND, DT, AMASS, X, QCALC)
!-----------------------------------------------------------------------
!     This routine perform initiation for the self-guided
!       Langevin Dynamcs (SGLD) simulaiton
!
        implicit none
        integer ISTART, IEND, I, I3, M
        logical QCALC
        _REAL_ DT, AMASS(*), X(*)
        _REAL_ DT2, FACT
        if (QCALC) then
!  Calculate shake forces
            DT2 = 1.0D0/(DT*DT)
            do I = ISTART, IEND
                FACT = DT2*AMASS(I)
                I3 = (I - 1)*3
                do M = 1, 3
                    I3 = I3 + 1
                    SHKSG(I3) = FACT*(X(I3) - DSG(I3))
                end do
            end do
        else
!  Save old positions
            do I = ISTART, IEND
                I3 = (I - 1)*3
                do M = 1, 3
                    I3 = I3 + 1
                    DSG(I3) = X(I3)
                end do
            end do
        end if
        return
    end subroutine SGFSHAKE

    subroutine SGLDW(NATOM, ISTART, IEND, NTP, &
        DTX, TEMP0, ENER, AMASS, WINV, X, F, V)
!-----------------------------------------------------------------------
!     This routine perform SGLD integration
!
        use state
        implicit none
#ifdef MPI
        include 'mpif.h'
        integer ierr
# include "parallel.h"
        _REAL_ temp1(20)
# ifndef USE_MPI_IN_PLACE
        _REAL_ :: temp2(20)
# endif
#endif
        integer NATOM, ISTART, IEND, NTP
        _REAL_ DTX, TEMP0
        type(state_rec) :: ener
        _REAL_ AMASS(*), WINV(*), X(*), F(*), V(*)
!
        integer I, M, I3, JSTA, JEND
        _REAL_ BOLTZ, AMASSI, TEMPI
        _REAL_ FACT, WFAC, DRAGI, PV1, PV2, RSD, FLN
        _REAL_ PTOT(3), FTOT(3), EKIN, EKINSG
        _REAL_ SUMDD, SUMFF, SUMGD, SUMGF, SUMPP, SUMGP
        _REAL_ VIT, VI3, FI3, GI3, HI3, VSGI, DSGI, FSGI, PSGI, FRICI
        parameter(BOLTZ=1.987192d-3)
!
        JSTA = ISTART
        JEND = IEND
        if (JSTA < ISGSTA) JSTA = ISGSTA
        if (JEND > ISGEND) JEND = ISGEND
!   Estimate guiding factor or guiding temperature
        do M = 1, 3
            PTOT(M) = 0.0D0
            FTOT(M) = 0.0D0
        end do
        EKIN = 0.0D0
        EKINSG = 0.0D0
        PV1 = 0.0d0
        PV2 = 0.0d0
        do I = 1, NATOM
            AMASSI = AMASS(I)
            WFAC = DTX*0.5D0*WINV(I)
            RSD = SQRT(2.D0*GAMMAS*BOLTZ*TEMP0*AMASSI/DTX)
            I3 = 3*(I - 1)
            do M = 1, 3
!   Keep random number series the same as that in a single cpu simulation
                call GAUSS(0.D0, RSD, FLN)
                I3 = 3*(I - 1) + M
                VI3 = V(I3)
                if (I >= JSTA .AND. I <= JEND) then
                    VSGI = SGAVG0*VSG(I3) + SGAVG1*VI3
                    VSG(I3) = VSGI
                    FI3 = F(I3)
                    DSGI = FI3 + SHKSG(I3)
                    FSGI = SGAVG0*FSG(I3) + SGAVG1*DSGI
                    FSG(I3) = FSGI
                    DSGI = DSGI - FSGI
                    SHKSG(I3) = DSGI
                    PSGI = GAMMAS*AMASSI*VSGI
                    DRAGI = SGFTI*PSGI + SGFFI*FSGI + SGFDI*DSGI
                    DSG(I3) = DRAGI
                    FI3 = FI3 + fln + DRAGI
                    F(I3) = FI3
                    VIT = VI3 + FI3*wfac
                    PV1 = PV1 + VIT*DRAGI
                    PV2 = PV2 + AMASSI*VIT*VIT
                    EKIN = EKIN + AMASSI*VI3*VI3
                    EKINSG = EKINSG + AMASSI*VSGI*VSGI
                    PTOT(M) = PTOT(M) + AMASSI*VI3
                    FTOT(M) = FTOT(M) + FI3
                else if (I >= ISTART .AND. I <= IEND) then
                    FI3 = F(I3) + fln
                    F(I3) = FI3
                    PTOT(M) = PTOT(M) + AMASSI*VI3
                    FTOT(M) = FTOT(M) + FI3
                end if
            end do
        end do
#ifdef MPI
        if (SANDERSIZE > 1) then
!  Combining all node results
!
            TEMP1(1) = PV1
            TEMP1(2) = PV2
            TEMP1(3) = PTOT(1)
            TEMP1(4) = PTOT(2)
            TEMP1(5) = PTOT(3)
            TEMP1(6) = FTOT(1)
            TEMP1(7) = FTOT(2)
            TEMP1(8) = FTOT(3)
            TEMP1(9) = EKIN
            TEMP1(10) = EKINSG
# ifdef USE_MPI_IN_PLACE
            call mpi_allreduce(MPI_IN_PLACE, temp1, 10, &
                MPI_DOUBLE_PRECISION, MPI_SUM, commsander, ierr)
            PV1 = TEMP1(1)
            PV2 = TEMP1(2)
            PTOT(1) = TEMP1(3)
            PTOT(2) = TEMP1(4)
            PTOT(3) = TEMP1(5)
            FTOT(1) = TEMP1(6)
            FTOT(2) = TEMP1(7)
            FTOT(3) = TEMP1(8)
            EKIN = TEMP1(9)
            EKINSG = TEMP1(10)
#else
            call MPI_ALLREDUCE(TEMP1, TEMP2, 10, &
                MPI_DOUBLE_PRECISION, MPI_SUM, COMMSANDER, IERR)
            PV1 = TEMP2(1)
            PV2 = TEMP2(2)
            PTOT(1) = TEMP2(3)
            PTOT(2) = TEMP2(4)
            PTOT(3) = TEMP2(5)
            FTOT(1) = TEMP2(6)
            FTOT(2) = TEMP2(7)
            FTOT(3) = TEMP2(8)
            EKIN = TEMP2(9)
            EKINSG = TEMP2(10)
# endif
        end if
#endif
        if (FIXCOM > 0) then
            do M = 1, 3
                PTOT(M) = PTOT(M)/SGMASS
                FTOT(M) = FTOT(M)/SGMASS
            end do
        else
            do M = 1, 3
                PTOT(M) = 0.0D0
                FTOT(M) = 0.0D0
            end do
        end if
        SGSCALE = (1.0D0 + 0.5D0*GAMMAS*DTX)*PV1/(PV2 - 0.5D0*DTX*PV1)
        SUMDD = 0.0d0
        SUMFF = 0.0d0
        SUMGD = 0.0d0
        SUMGF = 0.0d0
        SUMPP = 0.0d0
        SUMGP = 0.0d0
        do I = ISTART, IEND
            WFAC = WINV(I)*DTX
            AMASSI = AMASS(I)
            if (I < JSTA .OR. I > JEND) then
                FACT = 0.5D0*DTX*GAMMAS
            else
                FACT = 0.5D0*DTX*(GAMMAS + SGSCALE)
            end if
            do M = 1, 3
                I3 = 3*(I - 1) + M
                VI3 = V(I3) - PTOT(M)
                FI3 = F(I3) - AMASSI*FTOT(M)
                F(I3) = FI3
                VIT = ((1.0D0 - FACT)*VI3 + FI3*WFAC)/(1.0D0 + FACT)
                V(I3) = VIT
                VIT = 0.5D0*(VI3 + VIT)
                VSGI = VSG(I3)
                DSGI = SHKSG(I3)
                FSGI = FSG(I3)
                FRICI = AMASSI*SGSCALE*VIT
                PSGI = AMASSI*GAMMAS*VSGI
                FI3 = SGFTI*PSGI - FRICI
                DRAGI = DSG(I3) - FRICI
                DSG(I3) = DRAGI
                GI3 = SGAVG0*GSG(I3) + SGAVG1*FI3
                HI3 = SGAVG0*HSG(I3) + SGAVG1*DRAGI
                GSG(I3) = GI3
                HSG(I3) = HI3
                SUMDD = SUMDD + DSGI*DSGI
                SUMFF = SUMFF + FSGI*FSGI
                SUMGD = SUMGD + (FI3 - GI3)*DSGI
                SUMGF = SUMGF + (GI3)*FSGI
                SUMPP = SUMPP + PSGI*PSGI
                SUMGP = SUMGP + (HI3)*PSGI
            end do
        end do
        ! Estimate the viral due to the guiding forces
        if (NTP > 0) then
            VIRSG = 0.0D0
            do I = ISTART, IEND
                do M = 1, 3
                    I3 = 3*(I - 1) + M
                    VIRSG = VIRSG + X(I3)*DSG(I3)
                end do
            end do
            VIRSG = VIRSG/3.0D0
        end if
#ifdef MPI
        if (SANDERSIZE > 1) then
!  Combining all node results
!
            TEMP1(1) = SUMDD
            TEMP1(2) = SUMFF
            TEMP1(3) = SUMGD
            TEMP1(4) = SUMGF
            TEMP1(5) = SUMPP
            TEMP1(6) = SUMGP
            TEMP1(7) = VIRSG
# ifdef USE_MPI_IN_PLACE
            call mpi_allreduce(MPI_IN_PLACE, temp1, 7, &
                MPI_DOUBLE_PRECISION, MPI_SUM, commsander, ierr)
            SUMDD = TEMP1(1)
            SUMFF = TEMP1(2)
            SUMGD = TEMP1(3)
            SUMGF = TEMP1(4)
            SUMPP = TEMP1(5)
            SUMGP = TEMP1(6)
            VIRSG = TEMP1(7)
#else
            call MPI_ALLREDUCE(TEMP1, TEMP2, 7, &
                MPI_DOUBLE_PRECISION, MPI_SUM, COMMSANDER, IERR)
            SUMDD = TEMP2(1)
            SUMFF = TEMP2(2)
            SUMGD = TEMP2(3)
            SUMGF = TEMP2(4)
            SUMPP = TEMP2(5)
            SUMGP = TEMP2(6)
            VIRSG = TEMP2(7)
# endif
        end if
#endif
        ! Estimate momentum guiding factor and guiding temperature
        AVGFF = SGAVP0*AVGFF + SGAVP1*SUMFF
        AVGDD = SGAVP0*AVGDD + SGAVP1*SUMDD
        AVGGF = SGAVP0*AVGGF + SGAVP1*SUMGF
        AVGGD = SGAVP0*AVGGD + SGAVP1*SUMGD
        AVGPP = SGAVP0*AVGPP + SGAVP1*SUMPP
        AVGGP = SGAVP0*AVGGP + SGAVP1*SUMGP
        ! Estimate SGLD factors
        SGEFLF = 1.0D0 + AVGGF/AVGFF
        SGEFHF = 1.0D0 + AVGGD/AVGDD
        SGCFLF = 1.0D0 - AVGGP/AVGPP
        ! Estimate low frequency temperatures
        TEMPI = EKIN/RNDFSG/BOLTZ
        TEMPLFI = SGAVP0*TEMPLFI + SGAVP1*EKINSG/RNDFSG/BOLTZ
        TEMPHFI = TSGSET - TEMPLFI
        ! Estimate reference temperatures
        TREFLFI = TREFLF
        if (TREFLF < 1.0D-5) then
            TREFLFI = MAX(TEMPLFI*SGCFLF, TEMPLFI/1.0D1)
            if (TRXSGLD .AND. TEMPRXLF > 1.0D-5) TREFLFI = TEMPRXLF
        end if
        TREFHFI = TSGSET - TREFLFI
        ! Estimate guiding temperatures
        TEMPSGI = SGAVP0*TEMPSGI + SGAVP1*TEMPI*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
        ! Estimate weighting factor
        FRCLF = SGEFLF + SGFFI
        FRCHF = SGEFHF + SGFDI
        ! Adjust guiding factors if wanted
        if (TEMPSG > 1.0D-5) then
            SGFTI = SGFTI + SGAVP1*(TEMPSG - TEMPSGI)*(TEMPSG - TSGSET) &
                /TSGSET/(TEMPSG + TEMPSGI)
            SGFTI = MIN(SGFTI, 1.0D1)
            SGFTI = MAX(SGFTI, -1.0D1)
        end if
        if (TSGLDFP) then
            ! Estimate force guiding factor
            if (SGFD < -1.0D2 .AND. SGFF < -1.0D2) then
                SGFDI = SGAVP0*SGFDI + SGAVP1*(TEMPHFI/TREFHFI - SGEFHF)
                SGFFI = SGAVP0*SGFFI + SGAVP1*(TEMPLFI/TREFLFI - SGEFLF)
            else if (SGFD < -1.0D2) then
                SGFDI = SGAVG0*SGFDI + SGAVP1*(TEMPLFI/TREFLFI - SGEFLF - SGFF)
            else if (SGFF < -1.0D2) then
                SGFFI = SGAVP0*SGFFI + SGAVP1*(TEMPLFI/TREFLFI - SGEFLF)
            end if
            SGSCALE = SGFDI
        else
            ! Convert SGSCALE to unitless
            SGSCALE = SGSCALE/GAMMAS
        end if
        ! update accumulators
        call SGENERGY(ENER)
        return
    end subroutine SGLDW

    subroutine SGMDW(NATOM, ISTART, IEND, NTP, &
        DTX, ENER, AMASS, WINV, X, F, V)
!-----------------------------------------------------------------------
!     This routine calculate guiding force using SGLD method
!     for MD simulation
!
        use state
        implicit none
#ifdef MPI
        include 'mpif.h'
        integer ierr
# include "parallel.h"
        _REAL_ temp1(20)
# ifndef USE_MPI_IN_PLACE
        _REAL_ :: temp2(20)
# endif
#endif
        integer NATOM, ISTART, IEND, NTP
        _REAL_ DTX
        type(state_rec) :: ener
        _REAL_ AMASS(*), WINV(*), X(*), F(*), V(*)
!
        integer JSTA, JEND, I, M, I3
        _REAL_ BOLTZ, AMASSI, TEMPI
        _REAL_ FACT, WFAC, FRICI, DRAGI, PV1, PV2
        _REAL_ PTOT(3), FTOT(3), EKIN, EKINSG
        _REAL_ SUMDD, SUMFF, SUMGD, SUMGF
        _REAL_ VI3, VIT, FI3, GI3, HI3, VSGI, DSGI, FSGI, PSGI
        parameter(BOLTZ=1.987192d-3)
!
        JSTA = ISTART
        JEND = IEND
        if (JSTA < ISGSTA) JSTA = ISGSTA
        if (JEND > ISGEND) JEND = ISGEND
!    Estimate guiding factor or guiding temperature
        do M = 1, 3
            PTOT(M) = 0.0D0
            FTOT(M) = 0.0D0
        end do
        EKIN = 0.0D0
        EKINSG = 0.0D0
        PV1 = 0.0d0
        PV2 = 0.0d0
        do I = ISTART, IEND
            AMASSI = AMASS(I)
            wfac = DTX*0.5d0*WINV(I)
            do M = 1, 3
                I3 = 3*(I - 1) + M
                VI3 = V(I3)
                FI3 = F(I3)
                if (I >= JSTA .AND. I <= JEND) then
                    VSGI = SGAVG0*VSG(I3) + SGAVG1*VI3
                    VSG(I3) = VSGI
                    DSGI = FI3 + SHKSG(I3)
                    FSGI = SGAVG0*FSG(I3) + SGAVG1*DSGI
                    FSG(I3) = FSGI
                    DSGI = DSGI - FSGI
                    SHKSG(I3) = DSGI
                    PSGI = GAMMAS*AMASSI*VSGI
                    DRAGI = SGFTI*PSGI + SGFFI*FSGI + SGFDI*DSGI
                    DSG(I3) = DRAGI
                    FI3 = FI3 + DRAGI
                    F(I3) = FI3
                    VIT = VI3 + FI3*wfac
                    PV1 = PV1 + VIT*DRAGI
                    PV2 = PV2 + AMASSI*VIT*VIT
                    EKIN = EKIN + AMASSI*VI3*VI3
                    EKINSG = EKINSG + AMASSI*VSGI*VSGI
                end if
                PTOT(M) = PTOT(M) + AMASSI*VI3
                FTOT(M) = FTOT(M) + FI3
            end do
        end do
#ifdef MPI
        if (SANDERSIZE > 1) then
!  Combining all node results
!
            TEMP1(1) = PV1
            TEMP1(2) = PV2
            TEMP1(3) = PTOT(1)
            TEMP1(4) = PTOT(2)
            TEMP1(5) = PTOT(3)
            TEMP1(6) = FTOT(1)
            TEMP1(7) = FTOT(2)
            TEMP1(8) = FTOT(3)
            TEMP1(9) = EKIN
            TEMP1(10) = EKINSG
# ifdef USE_MPI_IN_PLACE
            call mpi_allreduce(MPI_IN_PLACE, temp1, 10, &
                MPI_DOUBLE_PRECISION, MPI_SUM, commsander, ierr)
            PV1 = TEMP1(1)
            PV2 = TEMP1(2)
            PTOT(1) = TEMP1(3)
            PTOT(2) = TEMP1(4)
            PTOT(3) = TEMP1(5)
            FTOT(1) = TEMP1(6)
            FTOT(2) = TEMP1(7)
            FTOT(3) = TEMP1(8)
            EKIN = TEMP1(9)
            EKINSG = TEMP1(10)
#else
            call MPI_ALLREDUCE(TEMP1, TEMP2, 10, &
                MPI_DOUBLE_PRECISION, MPI_SUM, COMMSANDER, IERR)
            PV1 = TEMP2(1)
            PV2 = TEMP2(2)
            PTOT(1) = TEMP2(3)
            PTOT(2) = TEMP2(4)
            PTOT(3) = TEMP2(5)
            FTOT(1) = TEMP2(6)
            FTOT(2) = TEMP2(7)
            FTOT(3) = TEMP2(8)
            EKIN = TEMP2(9)
            EKINSG = TEMP2(10)
# endif
        end if
#endif
        if (FIXCOM > 0) then
            do M = 1, 3
                PTOT(M) = PTOT(M)/SGMASS
                FTOT(M) = FTOT(M)/SGMASS
            end do
        else
            do M = 1, 3
                PTOT(M) = 0.0D0
                FTOT(M) = 0.0D0
            end do
        end if
        SGSCALE = PV1/(PV2 - 0.5D0*DTX*PV1)
        FACT = SGSCALE/(1.0D0 + 0.5D0*SGSCALE*DTX)
        SUMDD = 0.0d0
        SUMFF = 0.0d0
        SUMGD = 0.0d0
        SUMGF = 0.0d0
        do I = ISTART, IEND
            WFAC = 0.5D0*WINV(I)*DTX
            AMASSI = AMASS(I)
            do M = 1, 3
                I3 = 3*(I - 1) + M
                VI3 = V(I3) - PTOT(M)
                V(I3) = VI3
                FI3 = F(I3) - AMASSI*FTOT(M)
                if (I < JSTA .OR. I > JEND) then
                    F(I3) = FI3
                    cycle
                end if
                VIT = VI3 + FI3*WFAC
                FRICI = FACT*AMASSI*VIT
                F(I3) = FI3 - FRICI
                DRAGI = DSG(I3) - FRICI
                DSG(I3) = DRAGI
                VSGI = VSG(I3)
                DSGI = SHKSG(I3)
                FSGI = FSG(I3)
                PSGI = AMASSI*GAMMAS*VSGI
                FI3 = SGFTI*PSGI - FRICI
                GI3 = SGAVG0*GSG(I3) + SGAVG1*FI3
                HI3 = SGAVG0*HSG(I3) + SGAVG1*DRAGI
                GSG(I3) = GI3
                HSG(I3) = HI3
                SUMDD = SUMDD + DSGI*DSGI
                SUMFF = SUMFF + FSGI*FSGI
                SUMGD = SUMGD + (FI3 - GI3)*DSGI
                SUMGF = SUMGF + (GI3)*FSGI
            end do
        end do
        ! Estimate the viral due to the guiding forces
        if (NTP > 0) then
            VIRSG = 0.0D0
            do I = JSTA, JEND
                do M = 1, 3
                    I3 = 3*(I - 1) + M
                    VIRSG = VIRSG + X(I3)*DSG(I3)
                end do
            end do
            VIRSG = VIRSG/3.0D0
        end if
#ifdef MPI
        if (SANDERSIZE > 1) then
!  Combining all node results
!
            TEMP1(1) = SUMDD
            TEMP1(2) = SUMFF
            TEMP1(3) = SUMGD
            TEMP1(4) = SUMGF
            TEMP1(5) = VIRSG
# ifdef USE_MPI_IN_PLACE
            call mpi_allreduce(MPI_IN_PLACE, temp1, 5, &
                MPI_DOUBLE_PRECISION, MPI_SUM, commsander, ierr)
            SUMDD = TEMP1(1)
            SUMFF = TEMP1(2)
            SUMGD = TEMP1(3)
            SUMGF = TEMP1(4)
            VIRSG = TEMP1(5)
#else
            call MPI_ALLREDUCE(TEMP1, TEMP2, 5, &
                MPI_DOUBLE_PRECISION, MPI_SUM, COMMSANDER, IERR)
            SUMDD = TEMP2(1)
            SUMFF = TEMP2(2)
            SUMGD = TEMP2(3)
            SUMGF = TEMP2(4)
            VIRSG = TEMP2(5)
# endif
        end if
#endif
        ! Estimate momentum guiding factor and guiding temperature
        AVGFF = SGAVP0*AVGFF + SGAVP1*SUMFF
        AVGDD = SGAVP0*AVGDD + SGAVP1*SUMDD
        AVGGF = SGAVP0*AVGGF + SGAVP1*SUMGF
        AVGGD = SGAVP0*AVGGD + SGAVP1*SUMGD
        ! Estimate SGLD factors
        SGEFLF = SGAVP0*SGEFLF + SGAVP1*(AVGGF/AVGFF + 1.0D0)
        SGEFHF = SGAVP0*SGEFHF + SGAVP1*(AVGGD/AVGDD + 1.0D0)
        ! Estimate low frequency temperatures
        TEMPI = EKIN/RNDFSG/BOLTZ
        TEMPLFI = SGAVP0*TEMPLFI + SGAVP1*EKINSG/RNDFSG/BOLTZ
        TEMPHFI = TSGSET - TEMPLFI
        ! Estimate reference temperatures
        TREFLFI = TREFLF
        if (TREFLF < 1.0D-5) then
            TREFLFI = TEMPLFI
            if (TRXSGLD .AND. TEMPRXLF > 1.0D-5) TREFLFI = TEMPRXLF
        end if
        TREFHFI = TSGSET - TREFLFI
        ! Estimate guiding temperatures
        TEMPSGI = SGAVP0*TEMPSGI + SGAVP1*TEMPI*TREFHFI*TEMPLFI/TREFLFI/TEMPHFI
        ! Estimate weighting factor
        FRCLF = SGEFLF + SGFFI
        FRCHF = SGEFHF + SGFDI
        ! Adjust guiding factors if wanted
        if (TEMPSG > 1.0D-5) then
            SGFTI = SGFTI + SGAVP1*(TEMPSG - TEMPSGI)*(TEMPSG - TSGSET) &
                /TSGSET/(TEMPSG + TEMPSGI)
            SGFTI = MIN(SGFTI, 1.0D1)
            SGFTI = MAX(SGFTI, -1.0D1)
        end if
        if (TSGLDFP) then
            ! Estimate force guiding factor
            if (SGFD < -1.0D2 .AND. SGFF < -1.0D2) then
                SGFDI = SGAVP0*SGFDI + SGAVP1*(TEMPHFI/TREFHFI - SGEFHF)
                SGFFI = SGAVP0*SGFFI + SGAVP1*(TEMPLFI/TREFLFI - SGEFLF)
            else if (SGFD < -1.0D2) then
                SGFDI = SGAVG0*SGFDI + SGAVP1*(TEMPLFI/TREFLFI - SGEFLF - SGFF)
            else if (SGFF < -1.0D2) then
                SGFFI = SGAVP0*SGFFI + SGAVP1*(TEMPLFI/TREFLFI - SGEFLF)
            end if
            SGSCALE = SGFDI
        else
            ! Convert SGSCALE to unitless
            SGSCALE = SGSCALE/GAMMAS
        end if
        ! update accumulators
        call SGENERGY(ENER)
        return
    end subroutine SGMDW

    subroutine SGENERGY(ENER)
!-----------------------------------------------------------------------
!     This routine set the ener fields of SGLD variables
!
        use state
        implicit none
        type(state_rec) :: ener
        _REAL_ BOLTZ
        parameter(BOLTZ=1.987192d-3)
        _REAL_ EPOTI
        ! Weighting accumulators
        EPOTI = ENER%POT%TOT
        if (EPOTLF > 1.0D10) then
            EPOTLF = EPOTI
            AVGEFLF = FRCLF
            AVGEFHF = FRCHF
            AVGCFLF = TREFLFI/TEMPLFI
            AVGCFHF = TREFHFI/(TSGSET - TEMPLFI)
            AVGELF = EPOTLF
            AVGTLF = TEMPLFI
        else
            EPOTLF = SGAVG0*EPOTLF + SGAVG1*EPOTI
            AVGEFLF = SGAVP0*AVGEFLF + SGAVP1*FRCLF
            AVGEFHF = SGAVP0*AVGEFHF + SGAVP1*FRCHF
            AVGCFLF = SGAVP0*AVGCFLF + SGAVP1*TREFLFI/TEMPLFI
            AVGCFHF = SGAVP0*AVGCFHF + SGAVP1*TREFHFI/(TSGSET - TEMPLFI)
            AVGELF = SGAVP0*AVGELF + SGAVP1*EPOTLF
            AVGTLF = SGAVP0*AVGTLF + SGAVP1*TEMPLFI
        end if
        EPOTHF = EPOTI - EPOTLF
        SGWT = ((AVGEFLF*AVGCFLF - 1.0D0)*(EPOTLF - AVGELF) + &
            (AVGEFHF*AVGCFHF - 1.0D0)*EPOTHF + VIRSG)/(BOLTZ*TSGSET)
        ! Update ENER structure
        ENER%SGLD%SGFT = SGFTI
        ENER%SGLD%SGFF = SGFFI
        ENER%SGLD%SGSCAL = SGSCALE
        ENER%SGLD%TEMPSG = TEMPSGI
        ENER%SGLD%TEMPLF = TEMPLFI
        ENER%SGLD%TEMPHF = TEMPHFI
        ENER%SGLD%TREFLF = TREFLFI
        ENER%SGLD%TREFHF = TREFHFI
        ENER%SGLD%FRCLF = FRCLF
        ENER%SGLD%FRCHF = FRCHF
        ENER%SGLD%EPOTLF = EPOTLF
        ENER%SGLD%EPOTHF = EPOTHF
        ENER%SGLD%VIRSG = VIRSG
        ENER%SGLD%SGWT = SGWT
        return
    end subroutine SGENERGY

#ifdef MPI

!*********************************************************************
!               SUBROUTINE SGLD_EXCHG
!*********************************************************************
!  exchange all data in the SGLDR common block
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine sgld_exchg(irep)

        implicit none
        include 'mpif.h'
#  include "parallel.h"

        integer, intent(in) :: irep
        integer ierror, istatus(mpi_status_size)
        call mpi_sendrecv_replace(sgft, nsgld_real, mpi_double_precision, &
            irep, 511, irep, 511, commmaster, istatus, ierror)
        return
    end subroutine sgld_exchg

!*********************************************************************
!               SUBROUTINE REMD_SCALE_VELO
!*********************************************************************
! Scale velocities based on new temps after exchange
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine rxsgld_scale(stagid, nr, myscaling, amass, v)

        implicit none
        include 'mpif.h'
#  include "parallel.h"

        integer, intent(in) :: stagid, nr
        _REAL_, intent(in) :: myscaling
        _REAL_, dimension(*), intent(inout) :: v, amass
        integer ierror
        integer i, j, jsta, jend
        _REAL_ ek, elf, ehf, elh, amassi, vi, vlf, vhf, chk, hfscale
        _REAL_ temp1(10), temp2(10)
!--------------------
        call mpi_bcast(stagid, 1, mpi_integer, 0, commsander, ierror)

        ! All processes scale velocities.
        ! DAN ROE: This could potentially be divided up as in runmd
        !  since when there are mutiple threads per group each thread
        !  only ever knows about its own subset of velocities anyway.
#ifdef VERBOSE_REMD
        if (sanderrank == 0) then
            write (6, '(a,f8.3,a,f8.3)') &
                "RXSGLD: scaling guiding properties by ", myscalsg, &
                " to match a new guiding factor ", sgfti
        end if
#endif
! ---=== Broadcast RXSGLD guiding effect ===---
        if (SANDERSIZE > 1) call mpi_bcast(sgft, nsgld_real, mpi_double_precision, &
            0, commsander, ierror)
        if (myscalsg > 0.0d0) then
            jsta = iparpt(mytaskid) + 1
            jend = iparpt(mytaskid + 1)
            if (JSTA < ISGSTA) JSTA = ISGSTA
            if (JEND > ISGEND) JEND = ISGEND
            ek = 0.0d0
            ehf = 0.0d0
            elf = 0.0d0
            elh = 0.0d0
            do i = jsta, jend
                amassi = amass(i)
                do j = 3*i - 2, 3*i
                    vi = v(j)
                    vlf = vsg(j)
                    vhf = vi - vlf
                    ek = ek + amassi*vi*vi
                    ehf = ehf + amassi*vhf*vhf
                    elf = elf + amassi*vlf*vlf
                    elh = elh + amassi*vhf*vlf
                end do
            end do
            if (SANDERSIZE > 1) then
!  Combining all node results
!
                TEMP1(1) = ek
                TEMP1(2) = ehf
                TEMP1(3) = elf
                TEMP1(4) = elh
# ifdef USE_MPI_IN_PLACE
                call mpi_allreduce(MPI_IN_PLACE, temp1, 4, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, commsander, ierror)
                ek = TEMP1(1)
                ehf = TEMP1(2)
                elf = TEMP1(3)
                elh = TEMP1(4)
#else
                call MPI_ALLREDUCE(TEMP1, TEMP2, 4, &
                    MPI_DOUBLE_PRECISION, MPI_SUM, COMMSANDER, IERROR)
                ek = TEMP2(1)
                ehf = TEMP2(2)
                elf = TEMP2(3)
                elh = TEMP2(4)
# endif
            end if
            ! solve high frequency scaling factor
            chk = elh*elh + ehf*(myscaling*myscaling*ek - myscalsg*myscalsg*elf)
            hfscale = 1.0d0
            if (chk > 0.0d0) hfscale = (sqrt(chk) - elh)/ehf
            hfscale = max(hfscale, 0.5d0)
            hfscale = min(hfscale, 2.0d0)
            do i = 1, nr
                do j = 3*i - 2, 3*i
                    if (i < jsta .or. i > jend) then
                        v(j) = myscaling*v(j)
                    else
                        vsg(j) = vsg(j)*myscalsg
                        gsg(j) = gsg(j)*myscalsg
                        hsg(j) = hsg(j)*myscalsg
                        v(j) = hfscale*v(j) + (myscalsg - hfscale)*vsg(j)
                    end if
                end do
            end do
        end if
        return

    end subroutine rxsgld_scale

!*********************************************************************
!               FUNCTION TEMPSGLOOKUP
!*********************************************************************
! lookup temp in templist and return its index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer function tempsglookup(numreps, temp, tsg, sgft, temps, tsgs, sgfts)

        implicit none
#  include "parallel.h"

        integer numreps
        _REAL_, intent(in) :: temp, tsg, sgft
        _REAL_, dimension(numreps), intent(in) :: temps, tsgs, sgfts

        integer i

        tempsglookup = 0
        do i = 1, numreps
            if (abs(temp - temps(i)) < 1.0d-6 &
                .and. abs(tsg - tsgs(i)) < 1.0d-6 &
                .and. abs(sgft - sgfts(i)) < 1.0d-6) then
                if (tempsglookup > 0) then
                    write (6, *) "================================"
                    write (6, *) "Two replicas are the same: ", tempsglookup, i
                    write (6, *) "================================"
                    call mexit(6, 1)
                end if
                tempsglookup = i
            end if
        end do
        return
    end function tempsglookup

!*********************************************************************
!               FUNCTION STAGIDLOOKUP
!*********************************************************************
! lookup stagid in stagidlist and return its neighboring index
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    integer function stagidlookup(numreps, id, idtable)

        implicit none
#  include "parallel.h"

        integer numreps
        integer, dimension(numreps), intent(in) :: idtable
        integer, intent(in) :: id

        integer i

        stagidlookup = -1
        do i = 1, numreps
            if (id == idtable(i)) stagidlookup = i
        end do
        return
    end function stagidlookup

!*********************************************************************
!               SUBROUTINE SORTTEMPSG
!*********************************************************************
! sort temp ascendingly
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine sorttempsg(numreps, temp, tsgs, sgfts)

        implicit none

#  include "parallel.h"

        integer numreps
        _REAL_, dimension(numreps), intent(inout) :: temp, tsgs, sgfts
        _REAL_, dimension(numreps) :: tmp
        integer, dimension(numreps) :: tmpid

        _REAL_ tempt
        integer i, j, ii

        do i = 1, numreps
            tmp(i) = 1000000*temp(i) + 100*tsgs(i) + sgfts(i)
            tmpid(i) = i
        end do
        do i = 1, numreps
            do j = i + 1, numreps
                if (tmp(j) < tmp(i)) then
                    tempt = tmp(i)
                    tmp(i) = tmp(j)
                    tmp(j) = tempt
                    ii = tmpid(i)
                    tmpid(i) = tmpid(j)
                    tmpid(j) = ii
                end if
            end do
        end do
        do i = 1, numreps
            ii = tmpid(i)
            tmp(i) = temp(ii)
        end do
        do i = 1, numreps
            ii = tmpid(i)
            temp(i) = tmp(i)
            tmp(i) = tsgs(ii)
        end do
        do i = 1, numreps
            ii = tmpid(i)
            tsgs(i) = tmp(i)
            tmp(i) = sgfts(ii)
        end do
        do i = 1, numreps
            sgfts(i) = tmp(i)
        end do
        return
    end subroutine sorttempsg

#endif /* MPI */

end module sgld
