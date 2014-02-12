#include "copyright.h"
#include "dprec.fh"

MODULE EVDWMOD

  IMPLICIT NONE

CONTAINS

  SUBROUTINE VDW_IJ( &
       &     QI, SI,Z0I,ZQI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I, &
       &     QJ, SJ,Z0J,ZQJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J, &
       &     RIJ,  &
       &     EVDW,dEdQI,dEdQJ)

    
    USE EREPMOD, ONLY : VDW_REP_IJ
    USE EDISPMOD, ONLY : VDW_DISP_IJ
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: QI, SI,Z0I,ZQI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I
    _REAL_,INTENT(IN) :: QJ, SJ,Z0J,ZQJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J
    _REAL_,INTENT(IN) :: RIJ
    _REAL_,INTENT(OUT) :: EVDW,dEdQI,dEdQJ

    _REAL_ :: Erep,Edisp,dErepdQi,dErepdQj,dEdispdQi,dEdispdQj

    Erep      = 0.D0
    Edisp     = 0.D0
    dErepdQi  = 0.D0
    dErepdQj  = 0.D0
    dEdispdQi = 0.D0
    dEdispdQj = 0.D0

    CALL VDW_REP_IJ( &
         &     QI,SI,Z0I,ZQI, &
         &     QJ,SJ,Z0J,ZQJ, &
         &     RIJ,  &
         &     Erep,dErepdQI,dErepdQJ)

    CALL VDW_DISP_IJ( &
         &     QI, Z0I,ZQI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I, &
         &     QJ, Z0J,ZQJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J, &
         &     RIJ, &
         &     Edisp,dEdispdQI,dEdispdQJ)
    
    Evdw  = Erep + Edisp
    dEdQI = dErepdQI + dEdispdQI
    dEdQJ = dErepdQJ + dEdispdQJ
    
  END SUBROUTINE VDW_IJ


  SUBROUTINE VDW_IJ_DRI( &
       &     QI, SI,Z0I,ZQI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I, &
       &     QJ, SJ,Z0J,ZQJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J, &
       &     CXI,CYI,CZI,  CXJ,CYJ,CZJ, &
       &     GXI,GYI,GZI)
       
    
    USE EREPMOD, ONLY : VDW_REP_IJ_DRI
    USE EDISPMOD, ONLY : VDW_DISP_IJ_DRI
    IMPLICIT NONE

    _REAL_,INTENT(IN) :: QI, SI,Z0I,ZQI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I
    _REAL_,INTENT(IN) :: QJ, SJ,Z0J,ZQJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J
    _REAL_,INTENT(IN) :: CXI,CYI,CZI,  CXJ,CYJ,CZJ
    _REAL_,INTENT(OUT) :: GXI,GYI,GZI

    _REAL_ :: GXID,GYID,GZID,GXIR,GYIR,GZIR


    GXI=0.D0
    GYI=0.D0
    GZI=0.D0
    GXIR=0.D0
    GYIR=0.D0
    GZIR=0.D0
    GXID=0.D0
    GYID=0.D0
    GZID=0.D0
    
    CALL VDW_REP_IJ_DRI( &
         &    QI,SI,Z0I,ZQI, &
         &    QJ,SJ,Z0J,ZQJ, &
         &    CXI,CYI,CZI,  CXJ,CYJ,CZJ, &
         &    GXIR,GYIR,GZIR)
    
    
    CALL VDW_DISP_IJ_DRI( &
         &    QI, Z0I,ZQI, DP0I,DPQI, QP0I,QPQI, NEFF0I,NVAL0I, &
         &    QJ, Z0J,ZQJ, DP0J,DPQJ, QP0J,QPQJ, NEFF0J,NVAL0J, &
         &    CXI,CYI,CZI,  CXJ,CYJ,CZJ, &
         &    GXID,GYID,GZID)
    

    GXI = GXIR + GXID
    GYI = GYIR + GYID
    GZI = GZIR + GZID
    
  END SUBROUTINE VDW_IJ_DRI
      


END MODULE EVDWMOD
