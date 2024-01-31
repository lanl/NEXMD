#include "copyright.h"
#include "dprec.fh"

module EDISPMOD

    implicit none

contains

    subroutine VDW_DISP_IJ( &
    &     QI, Z0I, ZQI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I, &
    &     QJ, Z0J, ZQJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J, &
    &     RIJ,  &
    &     EDISP, dEdQI, dEdQJ)

        implicit none

        _REAL_, intent(IN) :: QI, Z0I, ZQI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I
        _REAL_, intent(IN) :: QJ, Z0J, ZQJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J
        _REAL_, intent(IN) :: RIJ
        _REAL_, intent(OUT) :: EDISP, dEdQI, dEdQJ

        _REAL_ :: TT(6), DTDB(6)
        _REAL_ :: ZERO
        _REAL_ :: DPI, QPI, ETA1I, ETA2I, DDPI, DQPI, DETA1I, DETA2I
        _REAL_ :: DPJ, QPJ, ETA1J, ETA2J, DDPJ, DQPJ, DETA1J, DETA2J

        _REAL_ :: R_11, dR_11I, dR_11J
        _REAL_ :: R_12, dR_12I, dR_12J
        _REAL_ :: R_21, dR_21I, dR_21J

        _REAL_ :: C6IJ, C8IJ, DC6IJ_I, DC6IJ_J, DC8IJ_I, DC8IJ_J
        _REAL_ :: ZI, ZJ, dZIdQI, dZJdQJ
        _REAL_ :: B, DBDZI, DBDZJ
        _REAL_ :: DBDQI, DBDQJ
        _REAL_ :: RIJ6, RIJ8

        ZERO = 0.D0
!C INTENT OUT VARIABLES INITIALIZED
        EDISP = ZERO
        dEdQI = ZERO
        dEdQJ = ZERO

        call Get_EtaAndPol( &
        & QI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I, &
        & DPI, QPI, ETA1I, ETA2I, &
        & DDPI, DQPI, DETA1I, DETA2I)

        call Get_EtaAndPol( &
        & QJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J, &
        & DPJ, QPJ, ETA1J, ETA2J, &
        & DDPJ, DQPJ, DETA1J, DETA2J)

        R_11 = ETA1I*ETA1J/(ETA1I + ETA1J)
        R_12 = ETA1I*ETA2J/(ETA1I + ETA2J)
        R_21 = ETA2I*ETA1J/(ETA2I + ETA1J)

        dR_11I = DETA1I*ETA1J/(ETA1I + ETA1J)  &
        &   - (ETA1I*ETA1J/(ETA1I + ETA1J)**2)*DETA1I
        dR_11J = ETA1I*DETA1J/(ETA1I + ETA1J) &
        &     - (ETA1I*ETA1J/(ETA1I + ETA1J)**2)*DETA1J

        DR_12I = DETA1I*ETA2J/(ETA1I + ETA2J) &
        &     - (ETA1I*ETA2J/(ETA1I + ETA2J)**2)*DETA1I
        DR_12J = ETA1I*DETA2J/(ETA1I + ETA2J) &
        &     - (ETA1I*ETA2J/(ETA1I + ETA2J)**2)*DETA2J

        DR_21I = DETA2I*ETA1J/(ETA2I + ETA1J) &
        &     - (ETA2I*ETA1J/(ETA2I + ETA1J)**2)*DETA2I
        DR_21J = ETA2I*DETA1J/(ETA2I + ETA1J) &
        &     - (ETA2I*ETA1J/(ETA2I + ETA1J)**2)*DETA1J

        C6IJ = 1.5D0*DPI*DPJ*R_11
        C8IJ = (15.D0/4.D0)*(DPI*QPJ*R_12 + QPI*DPJ*R_21)

        DC6IJ_I = 1.5D0*(DDPI*DPJ*R_11 + DPI*DPJ*DR_11I)
        DC6IJ_J = 1.5D0*(DPI*DDPJ*R_11 + DPI*DPJ*DR_11J)

        DC8IJ_I = (15.D0/4.D0)*( &
        &     DDPI*QPJ*R_12 + DPI*QPJ*DR_12I &
        &     + DQPI*DPJ*R_21 + QPI*DPJ*DR_21I)

        DC8IJ_J = (15.D0/4.D0)*( &
        &     DPI*DQPJ*R_12 + DPI*QPJ*DR_12J &
        &     + QPI*DDPJ*R_21 + QPI*DPJ*DR_21J)

        ZI = Z0I*EXP(-ZQI*QI)
        ZJ = Z0J*EXP(-ZQJ*QJ)

        dZIdQI = ZI*(-ZQI)
        dZJdQJ = ZJ*(-ZQJ)

        call BornMayerExp(ZI, ZJ, RIJ, B)
        call BornMayerExp_DZA(ZI, ZJ, RIJ, DBDZI)
        call BornMayerExp_DZB(ZI, ZJ, RIJ, DBDZJ)

        DBDQI = DBDZI*DZIDQI
        DBDQJ = DBDZJ*DZJDQJ

        call TangToennies(RIJ, B, TT)
        call TangToennies_DB(RIJ, B, DTDB)

        RIJ6 = RIJ**6
        RIJ8 = RIJ**8

        EDISP = -TT(1)*C6IJ/RIJ6 - TT(2)*C8IJ/RIJ8

        dEdQI = -DTDB(1)*DBDQI*C6IJ/RIJ6 &
        &   - TT(1)*DC6IJ_I/RIJ6 &
        &     - DTDB(2)*DBDQI*C8IJ/RIJ8 &
        &     - TT(2)*DC8IJ_I/RIJ8

        dEdQJ = -DTDB(1)*DBDQJ*C6IJ/RIJ6 &
        &     - TT(1)*DC6IJ_J/RIJ6 &
        &     - DTDB(2)*DBDQJ*C8IJ/RIJ8 &
        &     - TT(2)*DC8IJ_J/RIJ8

    end subroutine VDW_DISP_IJ

    subroutine VDW_DISP_IJ_DRI( &
    &     QI, Z0I, ZQI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I, &
    &     QJ, Z0J, ZQJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J, &
    &     CXI, CYI, CZI, CXJ, CYJ, CZJ, &
    &     GXI, GYI, GZI)
        implicit none

        _REAL_, intent(IN) :: QI, Z0I, ZQI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I
        _REAL_, intent(IN) :: QJ, Z0J, ZQJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J
        _REAL_, intent(IN) :: CXI, CYI, CZI, CXJ, CYJ, CZJ
        _REAL_, intent(OUT) :: GXI, GYI, GZI

        _REAL_ ::  DPI, QPI, ETA1I, ETA2I
        _REAL_ ::  DDPI, DQPI, DETA1I, DETA2I

        _REAL_ ::  DPJ, QPJ, ETA1J, ETA2J
        _REAL_ ::  DDPJ, DQPJ, DETA1J, DETA2J

        _REAL_ ::  DX, DY, DZ, RIJ2, RIJ, UX, UY, UZ

        _REAL_ ::  R_11, R_12, R_21
        _REAL_ ::  C6IJ, C8IJ, ZI, ZJ

        _REAL_ ::  B, DBDR, DEDR
        _REAL_ ::  RIJ6, RIJ7, RIJ8, RIJ9

        _REAL_ ::  ZERO

        _REAL_ ::  TT(6), DTDR(6)

        ZERO = 0.D0
!C INTENT OUT VARIABLES INITIALIZED
        GXI = ZERO
        GYI = ZERO
        GZI = ZERO

        DX = CXI - CXJ
        DY = CYI - CYJ
        DZ = CZI - CZJ
        RIJ2 = DX*DX + DY*DY + DZ*DZ
        RIJ = SQRT(RIJ2)

        UX = DX/RIJ
        UY = DY/RIJ
        UZ = DZ/RIJ

        call Get_EtaAndPol( &
        &     QI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I, &
        &     DPI, QPI, ETA1I, ETA2I, &
        &     DDPI, DQPI, DETA1I, DETA2I)

        call Get_EtaAndPol( &
        &     QJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J, &
        &     DPJ, QPJ, ETA1J, ETA2J, &
        &     DDPJ, DQPJ, DETA1J, DETA2J)

        R_11 = ETA1I*ETA1J/(ETA1I + ETA1J)
        R_12 = ETA1I*ETA2J/(ETA1I + ETA2J)
        R_21 = ETA2I*ETA1J/(ETA2I + ETA1J)

        C6IJ = 1.5D0*DPI*DPJ*R_11
        C8IJ = (15.D0/4.D0)*(DPI*QPJ*R_12 + QPI*DPJ*R_21)

        ZI = Z0I*EXP(-ZQI*QI)
        ZJ = Z0J*EXP(-ZQJ*QJ)

        call BornMayerExp(ZI, ZJ, RIJ, B)
        call BornMayerExp_DR(ZI, ZJ, RIJ, DBDR)
        call TangToennies(RIJ, B, TT)
        call TangToennies_DR(RIJ, B, DBDR, DTDR)

        RIJ6 = RIJ**6
        RIJ7 = RIJ6*RIJ
        RIJ8 = RIJ7*RIJ
        RIJ9 = RIJ8*RIJ

        DEDR = -DTDR(1)*C6IJ/RIJ6 &
        &     - (-6.D0)*TT(1)*C6IJ/RIJ7 &
        &     - DTDR(2)*C8IJ/RIJ8 &
        &     - (-8.D0)*TT(2)*C8IJ/RIJ9

        GXI = DEDR*UX
        GYI = DEDR*UY
        GZI = DEDR*UZ

    end subroutine VDW_DISP_IJ_DRI

!C====================================================================
!C====================================================================
!C====================================================================

    subroutine Get_EtaAndPol( &
    &     QI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I, &
    &     DPI, QPI, ETA1I, ETA2I, &
    &     DDPI, DQPI, DETA1I, DETA2I)
        implicit none

        _REAL_, intent(IN) :: QI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I
        _REAL_, intent(OUT) :: DPI, QPI, ETA1I, ETA2I
        _REAL_, intent(OUT) :: DDPI, DQPI, DETA1I, DETA2I

        _REAL_ :: NVALQI, NEFI, DNEFI, TMP
        _REAL_ :: ZERO, DTMP

        ZERO = 0.D0

!C     NUMBER OF EFFECTIVE ELECTRONS IN ATOM I
        NVALQI = NVAL0I - QI
        if (NVALQI .LE. ZERO) then
            NEFI = ZERO
            DNEFI = ZERO
        else
            NEFI = NVALQI*NEFF0I/NVAL0I
            DNEFI = -NEFF0I/NVAL0I
        end if

!C     DIPOLE POLARIZABILITY OF ATOM I
        DPI = DP0I*EXP(-DPQI*QI)
        DDPI = DPI*(-DPQI)

!C     QUADRUPOLAR POLARIZABILITY OF ATOM I
        QPI = QP0I*EXP(-QPQI*QI)
        DQPI = QPI*(-QPQI)

        if (DPI .GT. ZERO) then
            ETA1I = SQRT(NEFI/DPI)
            if (ETA1I .GT. 1.D-11) then
                DETA1I = (0.5D0/ETA1I)*(DNEFI/DPI - NEFI/DPI**2*DDPI)
            else
                DETA1I = ZERO
            end if
        else
            ETA1I = ZERO
            DETA1I = ZERO
        end if

        if (QPI .GT. ZERO) then
            TMP = SQRT(9.D0*NEFI*DPI)
            if (TMP .GT. 1.D-11) then
                DTMP = (0.5D0/TMP)*(9.D0*DNEFI*DPI + 9.D0*NEFI*DDPI)
            else
                DTMP = ZERO
            end if
            ETA2I = SQRT(TMP/QPI)
            DETA2I = (0.5D0/ETA2I)*(DTMP/QPI - TMP/QPI**2*DQPI)
        else
            ETA2I = ZERO
            DETA2I = ZERO
        end if

    end subroutine Get_EtaAndPol

    subroutine BornMayerExp(ZI, ZJ, RIJ, BIJ)

        use ErepMod, only : SS_Overlap, SS_RadDer
        implicit none

        _REAL_, intent(IN) :: ZI, ZJ, RIJ
        _REAL_, intent(OUT) :: BIJ

        _REAL_ :: ZERO, SIJ, DER

        ZERO = 0.D0

        call SS_Overlap(ZI, ZJ, ZERO, ZERO, ZERO, ZERO, ZERO, RIJ, SIJ)
        !C----- BORN MAYER EXPONENT
        call SS_RadDer(ZI, ZJ, RIJ, DER)
        BIJ = 1.0D0
        if (SIJ .GT. 1.D-30) then
            BIJ = DER/SIJ
        end if

    end subroutine BornMayerExp

    subroutine BornMayerExp_DR(ZI, ZJ, RIJ, DBDR)

        use ErepMod, only : SS_Overlap, SS_Overlap_DRA, SS_RadDer, SS_RadDer_DR
        implicit none

        _REAL_, intent(IN) :: ZI, ZJ, RIJ
        _REAL_, intent(OUT) :: DBDR

        _REAL_ :: ZERO
        _REAL_ :: DS(3), DSDR, S, DER, DDER

        ZERO = 0.D0

        call SS_Overlap(ZI, ZJ, ZERO, ZERO, RIJ, ZERO, ZERO, ZERO, S)
        call SS_Overlap_DRA(ZI, ZJ, ZERO, ZERO, RIJ, ZERO, ZERO, ZERO, DS)

        DSDR = DS(3)

        !C----- BORN MAYER EXPONENT
        call SS_RadDer(ZI, ZJ, RIJ, DER)
        call SS_RadDer_DR(ZI, ZJ, RIJ, DDER)

        DBDR = 0.0D0
        if (S .GT. 1.D-15) then
            DBDR = DDER/S - (DER/S**2)*DSDR
        end if

    end subroutine BornMayerExp_DR

    subroutine BornMayerExp_DZA(ZI, ZJ, RIJ, BIJ_DZA)

        use ErepMod, only : SS_Overlap, SS_RadDer, SS_RadDer_DZA, SS_Overlap_DZA
        implicit none

        _REAL_, intent(IN) :: ZI, ZJ, RIJ
        _REAL_, intent(OUT) :: BIJ_DZA

        _REAL_ :: ZERO, SIJ, DSIJ, DER, DDER

        ZERO = 0.D0

        call SS_Overlap(ZI, ZJ, ZERO, ZERO, ZERO, ZERO, ZERO, RIJ, SIJ)
        call SS_Overlap_DZA(ZI, ZJ, ZERO, ZERO, ZERO, ZERO, ZERO, RIJ, DSIJ)
        call SS_RadDer(ZI, ZJ, RIJ, DER)
        call SS_RadDer_DZA(ZI, ZJ, RIJ, DDER)

        BIJ_DZA = 0.0D0
        if (SIJ .GT. 1.D-15) then
            BIJ_DZA = DDER/SIJ - (DER/SIJ**2)*DSIJ
        end if

    end subroutine BornMayerExp_DZA

    subroutine BornMayerExp_DZB(ZI, ZJ, RIJ, BIJ_DZB)

        implicit none

        _REAL_, intent(IN) :: ZI, ZJ, RIJ
        _REAL_, intent(OUT) :: BIJ_DZB

        call BornMayerExp_DZA(ZJ, ZI, RIJ, BIJ_DZB)

    end subroutine BornMayerExp_DZB

    subroutine TangToennies(R, B, TT)

        implicit none

        _REAL_, intent(IN) :: R, B
        _REAL_, intent(OUT) :: TT(6)

        _REAL_ :: ZERO, ONE, TWO, THREE, FOUR, PT5, PT25
        integer :: I
        _REAL_ :: BR, EXPBR, KFACT, DFACT

        ZERO = 0.D0
        ONE = 1.D0
        TWO = 2.D0
        THREE = 3.D0
        FOUR = 4.D0
        PT5 = 0.5D0
        PT25 = 0.25D0

        do I = 1, 6
            TT(I) = ZERO
        end do
        BR = B*R
        EXPBR = EXP(-BR)
        KFACT = ONE
        DFACT = ONE
        do I = 1, 6, 1
            KFACT = KFACT*I
            DFACT = DFACT*BR
            TT(1) = TT(1) + DFACT/KFACT
        end do
        !C     THIS IS THE K=0 PART
        TT(1) = TT(1) + ONE

        TT(2) = TT(1)
        do I = 7, 8, 1
            KFACT = KFACT*I
            DFACT = DFACT*BR
            TT(2) = TT(2) + DFACT/KFACT
        end do
        TT(3) = TT(2)
        do I = 9, 10, 1
            KFACT = KFACT*I
            DFACT = DFACT*BR
            TT(3) = TT(3) + DFACT/KFACT
        end do
        TT(4) = TT(3)
        do I = 11, 12, 1
            KFACT = KFACT*I
            DFACT = DFACT*BR
            TT(4) = TT(4) + DFACT/KFACT
        end do
        TT(5) = TT(4)
        do I = 13, 14, 1
            KFACT = KFACT*I
            DFACT = DFACT*BR
            TT(5) = TT(5) + DFACT/KFACT
        end do
        TT(6) = TT(5)
        do I = 15, 16, 1
            KFACT = KFACT*I
            DFACT = DFACT*BR
            TT(6) = TT(6) + DFACT/KFACT
        end do
        do I = 1, 6
            TT(I) = ONE - TT(I)*EXPBR
            if (TT(I) .LT. ZERO) TT(I) = ZERO
        end do
    end subroutine TangToennies

    subroutine TangToennies_DB(R, B, DTT)
        implicit none

        _REAL_, intent(IN) :: R, B
        _REAL_, intent(OUT) :: DTT(6)

        _REAL_ :: TT(6)
        _REAL_ :: ZERO, ONE, TWO, THREE, FOUR, PT5, PT25
        integer :: I
        _REAL_ :: BR, EXPBR, KFACT, DFACT

        _REAL_ :: DEXPBR, DDFACT

        ZERO = 0.D0
        ONE = 1.D0
        TWO = 2.D0
        THREE = 3.D0
        FOUR = 4.D0
        PT5 = 0.5D0
        PT25 = 0.25D0

        do I = 1, 6
            TT(I) = ZERO
            DTT(I) = ZERO
        end do
        BR = B*R
        EXPBR = EXP(-BR)
        DEXPBR = -R*EXPBR

        KFACT = ONE
        DFACT = ONE
        DDFACT = ZERO
        do I = 1, 6, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*R
            DFACT = DFACT*BR
            TT(1) = TT(1) + DFACT/KFACT
            DTT(1) = DTT(1) + DDFACT/KFACT
        end do
!C     THIS IS THE K=0 PART
        TT(1) = TT(1) + ONE

        TT(2) = TT(1)
        DTT(2) = DTT(1)
        do I = 7, 8, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*R
            DFACT = DFACT*BR
            TT(2) = TT(2) + DFACT/KFACT
            DTT(2) = DTT(2) + DDFACT/KFACT
        end do
        TT(3) = TT(2)
        DTT(3) = DTT(2)
        do I = 9, 10, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*R
            DFACT = DFACT*BR
            TT(3) = TT(3) + DFACT/KFACT
            DTT(3) = DTT(3) + DDFACT/KFACT
        end do
        TT(4) = TT(3)
        DTT(4) = DTT(3)
        do I = 11, 12, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*R
            DFACT = DFACT*BR
            TT(4) = TT(4) + DFACT/KFACT
            DTT(4) = DTT(4) + DDFACT/KFACT
        end do
        TT(5) = TT(4)
        DTT(5) = DTT(4)
        do I = 13, 14, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*R
            DFACT = DFACT*BR
            TT(5) = TT(5) + DFACT/KFACT
            DTT(5) = DTT(5) + DDFACT/KFACT
        end do
        TT(6) = TT(5)
        DTT(6) = DTT(5)
        do I = 15, 16, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*R
            DFACT = DFACT*BR
            TT(6) = TT(6) + DFACT/KFACT
            DTT(6) = DTT(6) + DDFACT/KFACT
        end do
        do I = 1, 6
            DTT(I) = -TT(I)*(-R*EXPBR) - DTT(I)*EXPBR
            TT(I) = ONE - TT(I)*EXPBR
            if (TT(I) .LT. ZERO) then
                TT(I) = ZERO
                DTT(I) = ZERO
            end if
        end do
    end subroutine TangToennies_DB

    subroutine TangToennies_DR(R, B, DBDR, DTT)

        implicit none

        _REAL_, intent(IN) :: R, B, DBDR
        _REAL_, intent(OUT) :: DTT(6)

        _REAL_ :: TT(6)
        _REAL_ :: ZERO, ONE, TWO, THREE, FOUR, PT5, PT25
        integer:: I
        _REAL_ :: BR, EXPBR, KFACT, DFACT

        _REAL_ :: DDFACT, DBR, DEXP

        ZERO = 0.D0
        ONE = 1.D0
        TWO = 2.D0
        THREE = 3.D0
        FOUR = 4.D0
        PT5 = 0.5D0
        PT25 = 0.25D0

        do I = 1, 6
            TT(I) = ZERO
            DTT(I) = ZERO
        end do
        BR = B*R
        DBR = DBDR*R + B
        EXPBR = EXP(-BR)
        DEXP = (-DBR)*EXP(-BR)

        KFACT = ONE
        DFACT = ONE
        DDFACT = ZERO
        do I = 1, 6, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*DBR
            DFACT = DFACT*BR
            DTT(1) = DTT(1) + DDFACT/KFACT
            TT(1) = TT(1) + DFACT/KFACT
        end do
!C     THIS IS THE K=0 PART
        TT(1) = TT(1) + ONE

        DTT(2) = DTT(1)
        TT(2) = TT(1)
        do I = 7, 8, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*DBR
            DFACT = DFACT*BR
            DTT(2) = DTT(2) + DDFACT/KFACT
            TT(2) = TT(2) + DFACT/KFACT
        end do
        DTT(3) = DTT(2)
        TT(3) = TT(2)
        do I = 9, 10, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*DBR
            DFACT = DFACT*BR
            DTT(3) = DTT(3) + DDFACT/KFACT
            TT(3) = TT(3) + DFACT/KFACT
        end do
        DTT(4) = DTT(3)
        TT(4) = TT(3)
        do I = 11, 12, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*DBR
            DFACT = DFACT*BR
            DTT(4) = DTT(4) + DDFACT/KFACT
            TT(4) = TT(4) + DFACT/KFACT
        end do
        DTT(5) = DTT(4)
        TT(5) = TT(4)
        do I = 13, 14, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*DBR
            DFACT = DFACT*BR
            DTT(5) = DTT(5) + DDFACT/KFACT
            TT(5) = TT(5) + DFACT/KFACT
        end do
        DTT(6) = DTT(5)
        TT(6) = TT(5)
        do I = 15, 16, 1
            KFACT = KFACT*I
            DDFACT = DDFACT*BR + DFACT*DBR
            DFACT = DFACT*BR
            DTT(6) = DTT(6) + DDFACT/KFACT
            TT(6) = TT(6) + DFACT/KFACT
        end do
        do I = 1, 6
            DTT(I) = -DTT(I)*EXPBR - TT(I)*DEXP
            TT(I) = ONE - TT(I)*EXPBR
            if (TT(I) .LE. ZERO) then
                TT(I) = ZERO
                DTT(I) = ZERO
            end if
        end do
    end subroutine TangToennies_DR

end module EDISPMOD
