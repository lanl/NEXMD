#include "copyright.h"
#include "dprec.fh"

module EVDWMOD

    implicit none

contains

    subroutine VDW_IJ( &
    &     QI, SI, Z0I, ZQI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I, &
    &     QJ, SJ, Z0J, ZQJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J, &
    &     RIJ,  &
    &     EVDW, dEdQI, dEdQJ)

        use EREPMOD, only : VDW_REP_IJ
        use EDISPMOD, only : VDW_DISP_IJ
        implicit none

        _REAL_, intent(IN) :: QI, SI, Z0I, ZQI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I
        _REAL_, intent(IN) :: QJ, SJ, Z0J, ZQJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J
        _REAL_, intent(IN) :: RIJ
        _REAL_, intent(OUT) :: EVDW, dEdQI, dEdQJ

        _REAL_ :: Erep, Edisp, dErepdQi, dErepdQj, dEdispdQi, dEdispdQj

        Erep = 0.D0
        Edisp = 0.D0
        dErepdQi = 0.D0
        dErepdQj = 0.D0
        dEdispdQi = 0.D0
        dEdispdQj = 0.D0

        call VDW_REP_IJ( &
        &     QI, SI, Z0I, ZQI, &
        &     QJ, SJ, Z0J, ZQJ, &
        &     RIJ,  &
        &     Erep, dErepdQI, dErepdQJ)

        call VDW_DISP_IJ( &
        &     QI, Z0I, ZQI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I, &
        &     QJ, Z0J, ZQJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J, &
        &     RIJ, &
        &     Edisp, dEdispdQI, dEdispdQJ)

        Evdw = Erep + Edisp
        dEdQI = dErepdQI + dEdispdQI
        dEdQJ = dErepdQJ + dEdispdQJ

    end subroutine VDW_IJ

    subroutine VDW_IJ_DRI( &
    &     QI, SI, Z0I, ZQI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I, &
    &     QJ, SJ, Z0J, ZQJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J, &
    &     CXI, CYI, CZI, CXJ, CYJ, CZJ, &
    &     GXI, GYI, GZI)

        use EREPMOD, only : VDW_REP_IJ_DRI
        use EDISPMOD, only : VDW_DISP_IJ_DRI
        implicit none

        _REAL_, intent(IN) :: QI, SI, Z0I, ZQI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I
        _REAL_, intent(IN) :: QJ, SJ, Z0J, ZQJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J
        _REAL_, intent(IN) :: CXI, CYI, CZI, CXJ, CYJ, CZJ
        _REAL_, intent(OUT) :: GXI, GYI, GZI

        _REAL_ :: GXID, GYID, GZID, GXIR, GYIR, GZIR

        GXI = 0.D0
        GYI = 0.D0
        GZI = 0.D0
        GXIR = 0.D0
        GYIR = 0.D0
        GZIR = 0.D0
        GXID = 0.D0
        GYID = 0.D0
        GZID = 0.D0

        call VDW_REP_IJ_DRI( &
        &    QI, SI, Z0I, ZQI, &
        &    QJ, SJ, Z0J, ZQJ, &
        &    CXI, CYI, CZI, CXJ, CYJ, CZJ, &
        &    GXIR, GYIR, GZIR)

        call VDW_DISP_IJ_DRI( &
        &    QI, Z0I, ZQI, DP0I, DPQI, QP0I, QPQI, NEFF0I, NVAL0I, &
        &    QJ, Z0J, ZQJ, DP0J, DPQJ, QP0J, QPQJ, NEFF0J, NVAL0J, &
        &    CXI, CYI, CZI, CXJ, CYJ, CZJ, &
        &    GXID, GYID, GZID)

        GXI = GXIR + GXID
        GYI = GYIR + GYID
        GZI = GZIR + GZID

    end subroutine VDW_IJ_DRI

end module EVDWMOD
