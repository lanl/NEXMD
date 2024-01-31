! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
subroutine qm2_dihed(XYZ, I, J, K, L, ANGLE)
!********************************************************************
!
!      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
!            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
!            ARE IN ARRAY XYZ.
!
!     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
!           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
!
!********************************************************************
    implicit none

!Passed in
    _REAL_, intent(in) :: XYZ(3, *)
    integer, intent(in) :: i, j, k, l
    _REAL_, intent(out) :: angle

!Local
    _REAL_ xi1, xj1, xl1, yi1, yj1, yl1, zi1, zj1, zl1
    _REAL_ xi2, xl2, yi2, yl2, costh, sinth, cosph, sinph
    _REAL_ yj2, yi3, yl3
    _REAL_ dist, cosa, ddd, YXDIST

    XI1 = XYZ(1, I) - XYZ(1, K)
    XJ1 = XYZ(1, J) - XYZ(1, K)
    XL1 = XYZ(1, L) - XYZ(1, K)
    YI1 = XYZ(2, I) - XYZ(2, K)
    YJ1 = XYZ(2, J) - XYZ(2, K)
    YL1 = XYZ(2, L) - XYZ(2, K)
    ZI1 = XYZ(3, I) - XYZ(3, K)
    ZJ1 = XYZ(3, J) - XYZ(3, K)
    ZL1 = XYZ(3, L) - XYZ(3, K)
!      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
    DIST = 1.0d0/SQRT(XJ1**2 + YJ1**2 + ZJ1**2)
    COSA = ZJ1*DIST
    if (COSA .GT. 1.0D0) COSA = 1.0D0
    if (COSA .LT. -1.0D0) COSA = -1.0D0
    DDD = 1.0D0 - COSA**2
    if (DDD .LE. 0.0) GO TO 10
    YXDIST = DIST*(1.0d0/SQRT(DDD))
    if (YXDIST .LT. 1.0D6) GO TO 20
10  continue
    XI2 = XI1
    XL2 = XL1
    YI2 = YI1
    YL2 = YL1
    COSTH = COSA
    SINTH = 0.D0
    GO TO 30
20  COSPH = YJ1*YXDIST
    SINPH = XJ1*YXDIST
    XI2 = XI1*COSPH - YI1*SINPH
    XL2 = XL1*COSPH - YL1*SINPH
    YI2 = XI1*SINPH + YI1*COSPH
    YJ2 = XJ1*SINPH + YJ1*COSPH
    YL2 = XL1*SINPH + YL1*COSPH
!      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
    COSTH = COSA
    SINTH = YJ2*DIST
30  continue
    YI3 = YI2*COSTH - ZI1*SINTH
    YL3 = YL2*COSTH - ZL1*SINTH
    call qm2_dang(XL2, YL3, XI2, YI3, ANGLE)
    if (ANGLE .LT. 0.) ANGLE = 4.0D0*ASIN(1.0D00) + ANGLE
    if (ANGLE .GE. 6.2831853D0) ANGLE = 0.D0
    return
end subroutine qm2_dihed

!------------------------------------------------------------
subroutine qm2_dang(A1, A2, B1, B2, RCOS)
!*********************************************************************
!
!    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
!          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
!
!*********************************************************************

    implicit none

!Passed in
    _REAL_, intent(inout) :: a1, a2, b1, b2
    _REAL_, intent(out) :: rcos

!Local
    _REAL_ zero, anorm, bnorm, sinth, costh

    ZERO = 1.0D-6
    if (ABS(A1) .LT. ZERO .AND. ABS(A2) .LT. ZERO) GO TO 10
    if (ABS(B1) .LT. ZERO .AND. ABS(B2) .LT. ZERO) GO TO 10
    ANORM = 1.0D0/SQRT(A1**2 + A2**2)
    BNORM = 1.0D0/SQRT(B1**2 + B2**2)
    A1 = A1*ANORM
    A2 = A2*ANORM
    B1 = B1*BNORM
    B2 = B2*BNORM
    SINTH = (A1*B2) - (A2*B1)
    COSTH = A1*B1 + A2*B2
    if (COSTH .GT. 1.0D0) COSTH = 1.0D0
    if (COSTH .LT. -1.0D0) COSTH = -1.0D0
    RCOS = ACOS(COSTH)
    if (ABS(RCOS) .LT. 4.0D-4) GO TO 10
    if (SINTH .GT. 0.D0) RCOS = 4.0D0*ASIN(1.0D00) - RCOS
    RCOS = -RCOS
    return
10  RCOS = 0.0D0
    return
end subroutine qm2_dang
