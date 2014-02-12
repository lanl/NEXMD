#include "copyright.h"
#include "dprec.fh"

MODULE EDISPMOD

  IMPLICIT NONE

CONTAINS


  SUBROUTINE VDW_DISP_IJ( &
     &     QI, Z0I,ZQI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I, &
     &     QJ, Z0J,ZQJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J, &
     &     RIJ,  &
     &     EDISP,dEdQI,dEdQJ)

    IMPLICIT NONE

    _REAL_,INTENT(IN) :: QI, Z0I,ZQI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I
    _REAL_,INTENT(IN) :: QJ, Z0J,ZQJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J
    _REAL_,INTENT(IN) :: RIJ
    _REAL_,INTENT(OUT) :: EDISP,dEdQI,dEdQJ
      


    _REAL_ :: TT(6),DTDB(6)
    _REAL_ :: ZERO
    _REAL_ :: DPI,QPI,ETA1I,ETA2I,DDPI,DQPI,DETA1I,DETA2I
    _REAL_ :: DPJ,QPJ,ETA1J,ETA2J,DDPJ,DQPJ,DETA1J,DETA2J

    _REAL_ :: R_11,dR_11I,dR_11J
    _REAL_ :: R_12,dR_12I,dR_12J
    _REAL_ :: R_21,dR_21I,dR_21J

    _REAL_ :: C6IJ,C8IJ, DC6IJ_I,DC6IJ_J, DC8IJ_I,DC8IJ_J
    _REAL_ :: ZI,ZJ, dZIdQI, dZJdQJ
    _REAL_ :: B,DBDZI,DBDZJ
    _REAL_ :: DBDQI,DBDQJ
    _REAL_ :: RIJ6,RIJ8

    ZERO  = 0.D0
!C INTENT OUT VARIABLES INITIALIZED
    EDISP = ZERO
    dEdQI = ZERO
    dEdQJ = ZERO
    
    CALL Get_EtaAndPol( &
         & QI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I, &
         & DPI,QPI,ETA1I,ETA2I, &
         & DDPI,DQPI,DETA1I,DETA2I)
    
    CALL Get_EtaAndPol( &
         & QJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J, &
         & DPJ,QPJ,ETA1J,ETA2J, &
         & DDPJ,DQPJ,DETA1J,DETA2J)
    

    R_11 = ETA1I*ETA1J/(ETA1I+ETA1J)
    R_12 = ETA1I*ETA2J/(ETA1I+ETA2J)
    R_21 = ETA2I*ETA1J/(ETA2I+ETA1J)
    
    dR_11I = DETA1I*ETA1J/(ETA1I+ETA1J)  &
         &   - (ETA1I*ETA1J/(ETA1I+ETA1J)**2) * DETA1I
    dR_11J = ETA1I*DETA1J/(ETA1I+ETA1J) &
         &     - (ETA1I*ETA1J/(ETA1I+ETA1J)**2) * DETA1J

    DR_12I = DETA1I*ETA2J/(ETA1I+ETA2J) &
         &     - (ETA1I*ETA2J/(ETA1I+ETA2J)**2) * DETA1I
    DR_12J = ETA1I*DETA2J/(ETA1I+ETA2J) &
         &     - (ETA1I*ETA2J/(ETA1I+ETA2J)**2) * DETA2J

    DR_21I = DETA2I*ETA1J/(ETA2I+ETA1J) &
         &     - (ETA2I*ETA1J/(ETA2I+ETA1J)**2) * DETA2I
    DR_21J = ETA2I*DETA1J/(ETA2I+ETA1J) &
         &     - (ETA2I*ETA1J/(ETA2I+ETA1J)**2) * DETA1J

    C6IJ = 1.5D0 * DPI*DPJ*R_11
    C8IJ = (15.D0/4.D0)*(DPI*QPJ*R_12 + QPI*DPJ*R_21)
      

    DC6IJ_I = 1.5D0 * ( DDPI*DPJ*R_11 + DPI*DPJ*DR_11I )
    DC6IJ_J = 1.5D0 * ( DPI*DDPJ*R_11 + DPI*DPJ*DR_11J )
    
    DC8IJ_I = (15.D0/4.D0)*( &
         &     DDPI*QPJ*R_12 + DPI*QPJ*DR_12I &
         &     + DQPI*DPJ*R_21 + QPI*DPJ*DR_21I )

    DC8IJ_J = (15.D0/4.D0)*( &
         &     DPI*DQPJ*R_12 + DPI*QPJ*DR_12J &
         &     + QPI*DDPJ*R_21 + QPI*DPJ*DR_21J )


    ZI = Z0I * EXP(-ZQI*QI)
    ZJ = Z0J * EXP(-ZQJ*QJ)
    
    dZIdQI = ZI * (-ZQI)
    dZJdQJ = ZJ * (-ZQJ)

    CALL BornMayerExp(ZI,ZJ,RIJ,B)
    CALL BornMayerExp_DZA(ZI,ZJ,RIJ,DBDZI)
    CALL BornMayerExp_DZB(ZI,ZJ,RIJ,DBDZJ)
    
    DBDQI = DBDZI * DZIDQI
    DBDQJ = DBDZJ * DZJDQJ
    
    CALL TangToennies(RIJ,B,TT)
    CALL TangToennies_DB(RIJ,B,DTDB)
    
    RIJ6 = RIJ**6
    RIJ8 = RIJ**8
    
    EDISP = - TT(1)*C6IJ/RIJ6 - TT(2)*C8IJ/RIJ8
    
    dEdQI = - DTDB(1)*DBDQI * C6IJ/RIJ6 &
         &   - TT(1)*DC6IJ_I/RIJ6 &
         &     - DTDB(2)*DBDQI * C8IJ/RIJ8 &
         &     - TT(2)*DC8IJ_I/RIJ8
    
    dEdQJ = - DTDB(1)*DBDQJ * C6IJ/RIJ6 &
         &     - TT(1)*DC6IJ_J/RIJ6 &
         &     - DTDB(2)*DBDQJ * C8IJ/RIJ8 &
         &     - TT(2)*DC8IJ_J/RIJ8

  END SUBROUTINE VDW_DISP_IJ




  SUBROUTINE VDW_DISP_IJ_DRI( &
       &     QI, Z0I,ZQI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I, &
       &     QJ, Z0J,ZQJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J, &
       &     CXI,CYI,CZI,  CXJ,CYJ,CZJ, &
       &     GXI,GYI,GZI)
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: QI, Z0I,ZQI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I
    _REAL_,INTENT(IN) :: QJ, Z0J,ZQJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J
    _REAL_,INTENT(IN) :: CXI,CYI,CZI,  CXJ,CYJ,CZJ
    _REAL_,INTENT(OUT) :: GXI,GYI,GZI


    _REAL_ ::  DPI,QPI,ETA1I,ETA2I
    _REAL_ ::  DDPI,DQPI,DETA1I,DETA2I

    _REAL_ ::  DPJ,QPJ,ETA1J,ETA2J
    _REAL_ ::  DDPJ,DQPJ,DETA1J,DETA2J

    _REAL_ ::  DX,DY,DZ,RIJ2,RIJ,UX,UY,UZ

    _REAL_ ::  R_11,R_12,R_21
    _REAL_ ::  C6IJ,C8IJ,ZI,ZJ

    _REAL_ ::  B,DBDR,DEDR
    _REAL_ ::  RIJ6,RIJ7,RIJ8,RIJ9

    _REAL_ ::  ZERO

    _REAL_ ::  TT(6),DTDR(6)


    ZERO  = 0.D0
!C INTENT OUT VARIABLES INITIALIZED
    GXI = ZERO
    GYI = ZERO
    GZI = ZERO

    DX   = CXI-CXJ
    DY   = CYI-CYJ
    DZ   = CZI-CZJ
    RIJ2 = DX*DX + DY*DY + DZ*DZ
    RIJ  = SQRT(RIJ2)

    UX = DX/RIJ
    UY = DY/RIJ
    UZ = DZ/RIJ



    CALL Get_EtaAndPol( &
         &     QI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I, &
         &     DPI,QPI,ETA1I,ETA2I, &
         &     DDPI,DQPI,DETA1I,DETA2I)


    CALL Get_EtaAndPol( &
         &     QJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J, &
         &     DPJ,QPJ,ETA1J,ETA2J, &
         &     DDPJ,DQPJ,DETA1J,DETA2J)


    R_11 = ETA1I*ETA1J/(ETA1I+ETA1J)
    R_12 = ETA1I*ETA2J/(ETA1I+ETA2J)
    R_21 = ETA2I*ETA1J/(ETA2I+ETA1J)

    C6IJ = 1.5D0 * DPI*DPJ*R_11
    C8IJ = (15.D0/4.D0)*(DPI*QPJ*R_12 + QPI*DPJ*R_21)

    ZI = Z0I * EXP(-ZQI*QI)
    ZJ = Z0J * EXP(-ZQJ*QJ)

    CALL BornMayerExp(ZI,ZJ,RIJ,B)
    CALL BornMayerExp_DR(ZI,ZJ,RIJ,DBDR)
    CALL TangToennies(RIJ,B,TT)
    CALL TangToennies_DR(RIJ,B,DBDR,DTDR)

    RIJ6 = RIJ**6
    RIJ7 = RIJ6 * RIJ
    RIJ8 = RIJ7 * RIJ
    RIJ9 = RIJ8 * RIJ

    DEDR = -DTDR(1)*C6IJ/RIJ6 &
         &     - (-6.D0) * TT(1)*C6IJ/RIJ7 &
         &     -DTDR(2)*C8IJ/RIJ8 &
         &     - (-8.D0) * TT(2)*C8IJ/RIJ9
         
         
    GXI = DEDR * UX
    GYI = DEDR * UY
    GZI = DEDR * UZ

  END SUBROUTINE VDW_DISP_IJ_DRI


!C====================================================================
!C====================================================================
!C====================================================================




  SUBROUTINE Get_EtaAndPol( &
       &     QI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I, &
       &     DPI,QPI,ETA1I,ETA2I, &
       &     DDPI,DQPI,DETA1I,DETA2I)
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: QI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I
    _REAL_,INTENT(OUT) :: DPI,QPI,ETA1I,ETA2I
    _REAL_,INTENT(OUT) :: DDPI,DQPI,DETA1I,DETA2I

    _REAL_ :: NVALQI,NEFI,DNEFI,TMP
    _REAL_ :: ZERO,DTMP

    ZERO = 0.D0

!C     NUMBER OF EFFECTIVE ELECTRONS IN ATOM I
    NVALQI = NVAL0I - QI
    IF ( NVALQI .LE. ZERO ) THEN
       NEFI = ZERO
       DNEFI = ZERO
    ELSE
       NEFI = NVALQI * NEFF0I / NVAL0I
       DNEFI = - NEFF0I / NVAL0I
    END IF
    
!C     DIPOLE POLARIZABILITY OF ATOM I
    DPI  = DP0I * EXP(-DPQI*QI)
    DDPI = DPI * (-DPQI)
      
!C     QUADRUPOLAR POLARIZABILITY OF ATOM I
    QPI  = QP0I * EXP(-QPQI*QI)
    DQPI = QPI * (-QPQI)
    
    IF ( DPI .GT. ZERO ) THEN
       ETA1I = SQRT( NEFI / DPI )
       IF ( ETA1I .GT. 1.D-11 ) THEN
          DETA1I = (0.5D0/ETA1I) * ( DNEFI/DPI - NEFI/DPI**2 * DDPI )
       ELSE
          DETA1I = ZERO
       END IF
    ELSE
       ETA1I = ZERO
       DETA1I = ZERO
    END IF
    
    IF ( QPI .GT. ZERO ) THEN
       TMP = SQRT( 9.D0 * NEFI * DPI )
       IF ( TMP .GT. 1.D-11 ) THEN
          DTMP = (0.5D0/TMP) * (9.D0*DNEFI*DPI + 9.D0*NEFI*DDPI)
       ELSE
          DTMP = ZERO
       END IF
       ETA2I = SQRT( TMP / QPI )
       DETA2I = (0.5D0/ETA2I) * ( DTMP/QPI - TMP/QPI**2 * DQPI )
    ELSE
       ETA2I  = ZERO
       DETA2I = ZERO
    END IF
    
  END SUBROUTINE Get_EtaAndPol






  SUBROUTINE BornMayerExp(ZI,ZJ,RIJ,BIJ)

    USE ErepMod, ONLY : SS_Overlap,SS_RadDer
    IMPLICIT NONE
    
    _REAL_,INTENT(IN) :: ZI,ZJ,RIJ
    _REAL_,INTENT(OUT) :: BIJ

    _REAL_ :: ZERO,SIJ,DER

    ZERO=0.D0

    CALL SS_Overlap(ZI,ZJ,ZERO,ZERO,ZERO,ZERO,ZERO,RIJ,SIJ)
    !C----- BORN MAYER EXPONENT
    CALL SS_RadDer(ZI,ZJ,RIJ,DER)
    BIJ=1.0D0
    IF ( SIJ .GT. 1.D-30 ) THEN
       BIJ = DER/SIJ
    END IF
    
  END SUBROUTINE BornMayerExp


  SUBROUTINE BornMayerExp_DR(ZI,ZJ,RIJ,DBDR)

    USE ErepMod, ONLY : SS_Overlap,SS_Overlap_DRA,SS_RadDer,SS_RadDer_DR
    IMPLICIT NONE
    
    _REAL_,INTENT(IN) :: ZI,ZJ,RIJ
    _REAL_,INTENT(OUT) :: DBDR

    _REAL_ :: ZERO
    _REAL_ :: DS(3),DSDR,S,DER,DDER

    ZERO=0.D0
    
    CALL SS_Overlap(ZI,ZJ,ZERO,ZERO,RIJ,ZERO,ZERO,ZERO,S)
    CALL SS_Overlap_DRA(ZI,ZJ,ZERO,ZERO,RIJ,ZERO,ZERO,ZERO,DS)
    
    DSDR = DS(3)
    
    !C----- BORN MAYER EXPONENT
    CALL SS_RadDer(ZI,ZJ,RIJ,DER)
    CALL SS_RadDer_DR(ZI,ZJ,RIJ,DDER)
    
    DBDR=0.0D0
    IF ( S .GT. 1.D-15 ) THEN
       DBDR = DDER/S - (DER/S**2)*DSDR
    END IF
    
  END SUBROUTINE BornMayerExp_DR




  SUBROUTINE BornMayerExp_DZA(ZI,ZJ,RIJ,BIJ_DZA)

    USE ErepMod, ONLY : SS_Overlap,SS_RadDer,SS_RadDer_DZA,SS_Overlap_DZA
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: ZI,ZJ,RIJ
    _REAL_,INTENT(OUT) :: BIJ_DZA

    _REAL_ :: ZERO,SIJ,DSIJ,DER,DDER


    ZERO=0.D0
    
    CALL SS_Overlap(ZI,ZJ,ZERO,ZERO,ZERO,ZERO,ZERO,RIJ,SIJ)
    CALL SS_Overlap_DZA(ZI,ZJ,ZERO,ZERO,ZERO,ZERO,ZERO,RIJ,DSIJ)
    CALL SS_RadDer(ZI,ZJ,RIJ,DER)
    CALL SS_RadDer_DZA(ZI,ZJ,RIJ,DDER)
    
    BIJ_DZA=0.0D0
    IF ( SIJ .GT. 1.D-15 ) THEN
       BIJ_DZA = DDER/SIJ - (DER/SIJ**2) * DSIJ
    END IF
    
  END SUBROUTINE BornMayerExp_DZA


  SUBROUTINE BornMayerExp_DZB(ZI,ZJ,RIJ,BIJ_DZB)

    IMPLICIT NONE
    
    _REAL_,INTENT(IN) :: ZI,ZJ,RIJ
    _REAL_,INTENT(OUT) :: BIJ_DZB

    CALL BornMayerExp_DZA(ZJ,ZI,RIJ,BIJ_DZB)
    
  END SUBROUTINE BornMayerExp_DZB





      
  SUBROUTINE TangToennies(R,B,TT)

    IMPLICIT NONE

    _REAL_,INTENT(IN) :: R,B
    _REAL_,INTENT(OUT) :: TT(6)

    _REAL_ :: ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
    INTEGER :: I
    _REAL_ :: BR,EXPBR,KFACT,DFACT

    ZERO  =0.D0
    ONE   =1.D0
    TWO   =2.D0
    THREE =3.D0
    FOUR  =4.D0
    PT5   =0.5D0
    PT25  =0.25D0

    
    DO I=1,6
       TT(I) = ZERO
    END DO
    BR = B*R
    EXPBR = EXP( -BR )
    KFACT = ONE
    DFACT = ONE
    DO I=1,6,1
       KFACT = KFACT*I
       DFACT = DFACT*BR
       TT(1) = TT(1) + DFACT/KFACT
    END DO
    !C     THIS IS THE K=0 PART
    TT(1) = TT(1) + ONE

    TT(2)=TT(1)
    DO I=7,8,1
       KFACT = KFACT*I
       DFACT = DFACT*BR
       TT(2) = TT(2) + DFACT/KFACT
    END DO
    TT(3) = TT(2)
    DO I=9,10,1
       KFACT = KFACT*I
       DFACT = DFACT*BR
       TT(3) = TT(3) + DFACT/KFACT
    END DO
    TT(4) = TT(3)
    DO I=11,12,1
       KFACT = KFACT*I
       DFACT = DFACT*BR
       TT(4) = TT(4) + DFACT/KFACT
    END DO
    TT(5) = TT(4)
    DO I=13,14,1
       KFACT = KFACT*I
       DFACT = DFACT*BR
       TT(5) = TT(5) + DFACT/KFACT
    END DO
    TT(6) = TT(5)
    DO I=15,16,1
       KFACT = KFACT*I
       DFACT = DFACT*BR
       TT(6) = TT(6) + DFACT/KFACT
    END DO
    DO I=1,6
       TT(I) = ONE - TT(I)*EXPBR
       IF ( TT(I) .LT. ZERO ) TT(I) = ZERO
    END DO
  END SUBROUTINE TangToennies



  SUBROUTINE TangToennies_DB(R,B,DTT)
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: R,B
    _REAL_,INTENT(OUT) :: DTT(6)

    _REAL_ :: TT(6)
    _REAL_ :: ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
    INTEGER :: I
    _REAL_ :: BR,EXPBR,KFACT,DFACT

    _REAL_ :: DEXPBR,DDFACT

    ZERO  =0.D0
    ONE   =1.D0
    TWO   =2.D0
    THREE =3.D0
    FOUR  =4.D0
    PT5   =0.5D0
    PT25  =0.25D0


    DO I=1,6
       TT(I) = ZERO
       DTT(I)= ZERO
    END DO
    BR = B*R
    EXPBR = EXP( -BR )
    DEXPBR = -R*EXPBR

    KFACT = ONE
    DFACT = ONE
    DDFACT= ZERO
    DO I=1,6,1
       KFACT = KFACT*I
       DDFACT = DDFACT*BR + DFACT*R
       DFACT = DFACT*BR
       TT(1) = TT(1) + DFACT/KFACT
       DTT(1) = DTT(1) + DDFACT/KFACT
    END DO
!C     THIS IS THE K=0 PART
    TT(1) = TT(1) + ONE

    TT(2)=TT(1)
    DTT(2) = DTT(1)
    DO I=7,8,1
       KFACT = KFACT*I
       DDFACT = DDFACT*BR + DFACT*R
       DFACT = DFACT*BR
       TT(2) = TT(2) + DFACT/KFACT
       DTT(2) = DTT(2) + DDFACT/KFACT
    END DO
    TT(3) = TT(2)
    DTT(3) = DTT(2)
    DO I=9,10,1
       KFACT = KFACT*I
       DDFACT = DDFACT*BR + DFACT*R
       DFACT = DFACT*BR
       TT(3) = TT(3) + DFACT/KFACT
       DTT(3) = DTT(3) + DDFACT/KFACT
    END DO
    TT(4) = TT(3)
    DTT(4) = DTT(3)
    DO I=11,12,1
       KFACT = KFACT*I
       DDFACT = DDFACT*BR + DFACT*R
       DFACT = DFACT*BR
       TT(4) = TT(4) + DFACT/KFACT
       DTT(4) = DTT(4) + DDFACT/KFACT
    END DO
    TT(5) = TT(4)
    DTT(5) = DTT(4)
    DO I=13,14,1
       KFACT = KFACT*I
       DDFACT = DDFACT*BR + DFACT*R
       DFACT = DFACT*BR
       TT(5) = TT(5) + DFACT/KFACT
       DTT(5) = DTT(5) + DDFACT/KFACT
    END DO
    TT(6) = TT(5)
    DTT(6) = DTT(5)
    DO I=15,16,1
       KFACT = KFACT*I
       DDFACT = DDFACT*BR + DFACT*R
       DFACT = DFACT*BR
       TT(6) = TT(6) + DFACT/KFACT
       DTT(6) = DTT(6) + DDFACT/KFACT
    END DO
    DO I=1,6
       DTT(I) = - TT(I)*(-R*EXPBR) - DTT(I)*EXPBR
       TT(I) = ONE - TT(I)*EXPBR
       IF ( TT(I) .LT. ZERO ) THEN
          TT(I) = ZERO
          DTT(I) = ZERO
       END IF
    END DO
  END SUBROUTINE TangToennies_DB




  SUBROUTINE TangToennies_DR(R,B,DBDR,DTT)

    IMPLICIT NONE

    _REAL_,INTENT(IN) :: R,B,DBDR
    _REAL_,INTENT(OUT) :: DTT(6)
    
    _REAL_ :: TT(6)
    _REAL_ :: ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
    INTEGER:: I
    _REAL_ :: BR,EXPBR,KFACT,DFACT
    
    _REAL_ :: DDFACT,DBR,DEXP
    
    ZERO  =0.D0
    ONE   =1.D0
    TWO   =2.D0
    THREE =3.D0
    FOUR  =4.D0
    PT5   =0.5D0
    PT25  =0.25D0
    

    DO I=1,6
       TT(I) = ZERO
       DTT(I)= ZERO
    END DO
    BR = B*R
    DBR   = DBDR*R + B
    EXPBR = EXP( -BR )
    DEXP  = (-DBR)*EXP(-BR)

    KFACT = ONE
    DFACT = ONE
    DDFACT= ZERO
    DO I=1,6,1
       KFACT  = KFACT*I
       DDFACT = DDFACT*BR + DFACT*DBR
       DFACT  = DFACT*BR
       DTT(1)  = DTT(1) + DDFACT/KFACT
       TT(1)  = TT(1) + DFACT/KFACT
    END DO
!C     THIS IS THE K=0 PART
    TT(1) = TT(1) + ONE

    DTT(2)=DTT(1)
    TT(2)=TT(1)
    DO I=7,8,1
       KFACT = KFACT*I
       DDFACT= DDFACT*BR + DFACT*DBR
       DFACT = DFACT*BR
       DTT(2) = DTT(2) + DDFACT/KFACT
       TT(2)  = TT(2) + DFACT/KFACT
    END DO
    DTT(3) = DTT(2)
    TT(3) = TT(2)
    DO I=9,10,1
       KFACT = KFACT*I
       DDFACT= DDFACT*BR + DFACT*DBR
       DFACT = DFACT*BR
       DTT(3) = DTT(3) + DDFACT/KFACT
       TT(3) = TT(3) + DFACT/KFACT
    END DO
    DTT(4) = DTT(3)
    TT(4) = TT(3)
    DO I=11,12,1
       KFACT = KFACT*I
       DDFACT= DDFACT*BR + DFACT*DBR
       DFACT = DFACT*BR
       DTT(4) = DTT(4) + DDFACT/KFACT
       TT(4) = TT(4) + DFACT/KFACT
    END DO
    DTT(5) = DTT(4)
    TT(5) = TT(4)
    DO I=13,14,1
       KFACT = KFACT*I
       DDFACT= DDFACT*BR + DFACT*DBR
       DFACT = DFACT*BR
       DTT(5) = DTT(5) + DDFACT/KFACT
       TT(5) = TT(5) + DFACT/KFACT
    END DO
    DTT(6) = DTT(5)
    TT(6) = TT(5)
    DO I=15,16,1
       KFACT = KFACT*I
       DDFACT= DDFACT*BR + DFACT*DBR
       DFACT = DFACT*BR
       DTT(6) = DTT(6) + DDFACT/KFACT
       TT(6) = TT(6) + DFACT/KFACT
    END DO
    DO I=1,6
       DTT(I) = - DTT(I)*EXPBR - TT(I)*DEXP
       TT(I) = ONE - TT(I)*EXPBR
       IF ( TT(I) .LE. ZERO ) THEN
          TT(I) = ZERO
          DTT(I)= ZERO
       END IF
    END DO
  END SUBROUTINE TangToennies_DR



END MODULE EDISPMOD
