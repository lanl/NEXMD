#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++
!This module contains rotation utilities
!used in MNDO-d implementation
!Created by Taisung Lee, 02/15/2011 
!+++++++++++++++++++++++++++++++++++++++

module Rotation

public GetRotationMatrix
public Rotate2Center2Electron, Rotate1Elec, RotateCore



private GenerateRotationMatrix, rotationMatrix, xij_saved, tolerance

    _REAL_, save::matrix_saved(15,45)
    _REAL_, save::xij_saved(3)=(/ 1.0D9, 1.0D9, 1.0D9 /) 
    _REAL_, parameter::tolerance=1.0D-5
    logical, parameter::keepOld=.false.
    
contains

subroutine GetRotationMatrix(xij, rotationMatrix, hasDOrbital)

    _REAL_, intent(in)::xij(3)
    _REAL_, intent(inout)::rotationMatrix(15,45)
    logical, intent(in)::hasDOrbital
        
    if ((sum( (xij-xij_saved)**2) >= tolerance).or. .not. keepOld) then
       call GenerateRotationMatrix(xij, rotationMatrix, hasDOrbital)
    else 
       rotationMatrix=matrix_saved 
    end if

    return
    
end subroutine GetRotationMatrix
    
    
subroutine GenerateRotationMatrix(xij, matrix, hasDOrbital)
      
      use constants, only : zero, one, two, half, fourth, &
                            A_TO_BOHRS
      
!     *                                                                 
!     ROTATION MATRIX FOR A GIVEN ATOM PAIR I-J (I.GT.J).               
!     R         INTERATOMIC DISTANCE, IN ATOMIC UNITS (O).              
!     matrix()      PRECOMBINED ELEMENTS OF THE ROTATION MATRIX (O).        
!     *                                                                 
      implicit none
      _REAL_,intent(in)::xij(3)
      _REAL_,intent(inout)::matrix(15,45)
      logical, intent(in)::hasDOrbital
      
      
      integer, parameter::INDX(9) = (/ 0,1,3,6,10,15,21,28,36 /) 
      _REAL_,PARAMETER::SMALL=1.0D-07
      _REAL_,PARAMETER::PT5SQ3=0.8660254037841D0, PT5=half, PT25=fourth

      integer::i,j,k,l, kl,ij
      _REAL_::P(3,3),D(5,5) 
      _REAL_::x11, x22, x33, B, R
      _REAL_:: CA, CB, C2A, C2B, SA, SB, S2A, S2B
      _REAL_::SQB

               
      X11    = xij(1) 
      X22    = xij(2)
      X33    = xij(3) 
      B      = X11*X11+X22*X22 
      R      = SQRT(B+X33*X33) 
      SQB    = SQRT(B) 
      SB     = SQB/R 
!     CHECK FOR SPECIAL CASE (BOTH ATOMS ON Z AXIS).                    
      IF(SB.GT.SMALL) THEN 
         CA  = X11/SQB 
         SA  = X22/SQB 
         CB  = X33/R 
      ELSE 
         SA  = ZERO 
         SB  = ZERO 
         IF(X33.LT.ZERO) THEN 
            CA  =-ONE 
            CB  =-ONE 
         ELSE IF(X33.GT.ZERO) THEN 
            CA  = ONE 
            CB  = ONE 
         ELSE 
            CA  = ZERO 
            CB  = ZERO 
         ENDIF 
      ENDIF 
!     CONVERT DISTANCE TO ATOMIC UNITS.                                 
      R      = R*A_TO_BOHRS 
! *** CALCULATE ROTATION MATRIX ELEMENTS.                               
      P(1,1) = CA*SB 
      P(2,1) = CA*CB 
      P(3,1) =-SA 
      P(1,2) = SA*SB 
      P(2,2) = SA*CB 
      P(3,2) = CA 
      P(1,3) = CB 
      P(2,3) =-SB 
      P(3,3) = ZERO 
      IF(hasDOrbital) THEN 
         C2A    = two*CA*CA-ONE 
         C2B    = two*CB*CB-ONE 
         S2A    = two*SA*CA 
         S2B    = two*SB*CB 
         D(1,1) = PT5SQ3*C2A*SB*SB 
         D(2,1) = PT5*C2A*S2B 
         D(3,1) =-S2A*SB 
         D(4,1) = C2A*(CB*CB+PT5*SB*SB) 
         D(5,1) =-S2A*CB 
         D(1,2) = PT5SQ3*CA*S2B 
         D(2,2) = CA*C2B 
         D(3,2) =-SA*CB 
         D(4,2) =-PT5*CA*S2B 
         D(5,2) = SA*SB 
         D(1,3) = CB*CB-PT5*SB*SB 
         D(2,3) =-PT5SQ3*S2B 
         D(3,3) = ZERO 
         D(4,3) = PT5SQ3*SB*SB 
         D(5,3) = ZERO 
         D(1,4) = PT5SQ3*SA*S2B 
         D(2,4) = SA*C2B 
         D(3,4) = CA*CB 
         D(4,4) =-PT5*SA*S2B 
         D(5,4) =-CA*SB 
         D(1,5) = PT5SQ3*S2A*SB*SB 
         D(2,5) = PT5*S2A*S2B 
         D(3,5) = C2A*SB 
         D(4,5) = S2A*(CB*CB+PT5*SB*SB) 
         D(5,5) = C2A*CB 
      ENDIF 
!     *                                                                 
!     PRECOMBINE ROTATION MATRIX ELEMENTS.                              
!     *                                                                 
!     THE FIRST INDEX OF matrix(IJ,KL) IS CONSECUTIVE. AS MANY ELEMENTS     
!     AS NEEDED ARE DEFINED (1 FOR KL=SS, 3 FOR KL=PS, 6 FOR KL=PP,     
!     5 FOR KL=DS, 15 FOR KL=DP, AND 15 FOR KL=DD).                     
!     THE SECOND INDEX OF matrix(IJ,KL) IS A STANDARD PAIR INDEX.           
!     KL=(K*(K-1)/2+L , ORDER OF K AND L AS IN INTEGRAL EVALUATION.     
!     *                                                                 
!     S-S                                                               
      matrix(1,1)   = ONE 
!     P-S                                                               
      DO 10 K=1,3 
      KL        = INDX(K+1)+1 
      matrix(1,KL)  = P(K,1) 
      matrix(2,KL)  = P(K,2) 
   10 matrix(3,KL)  = P(K,3) 
!     P-P                                                               
      DO 20 K=1,3 
      KL        = INDX(K+1)+K+1 
      matrix(1,KL)  = P(K,1)*P(K,1) 
      matrix(2,KL)  = P(K,1)*P(K,2)  
      matrix(3,KL)  = P(K,2)*P(K,2) 
      matrix(4,KL)  = P(K,1)*P(K,3) 
      matrix(5,KL)  = P(K,2)*P(K,3) 
   20 matrix(6,KL)  = P(K,3)*P(K,3) 
      DO 30 K=2,3 
      DO 30 L=1,K-1 
      KL        = INDX(K+1)+L+1 
      matrix(1,KL)  = P(K,1)*P(L,1)*two 
      matrix(2,KL)  = P(K,1)*P(L,2)+P(K,2)*P(L,1) 
      matrix(3,KL)  = P(K,2)*P(L,2)*two 
      matrix(4,KL)  = P(K,1)*P(L,3)+P(K,3)*P(L,1) 
      matrix(5,KL)  = P(K,2)*P(L,3)+P(K,3)*P(L,2) 
   30 matrix(6,KL)  = P(K,3)*P(L,3)*two 
   
      IF(hasDOrbital) then 
!     D-S                                                               
      DO 40 K=1,5 
      KL        = INDX(K+4)+1 
      matrix(1,KL)  = D(K,1) 
      matrix(2,KL)  = D(K,2) 
      matrix(3,KL)  = D(K,3) 
      matrix(4,KL)  = D(K,4) 
   40 matrix(5,KL)  = D(K,5) 
!     D-P                                                               
      DO 50 K=1,5 
      DO 50 L=1,3 
      KL        = INDX(K+4)+L+1 
      matrix(1,KL)  = D(K,1)*P(L,1) 
      matrix(2,KL)  = D(K,1)*P(L,2) 
      matrix(3,KL)  = D(K,1)*P(L,3) 
      matrix(4,KL)  = D(K,2)*P(L,1) 
      matrix(5,KL)  = D(K,2)*P(L,2) 
      matrix(6,KL)  = D(K,2)*P(L,3) 
      matrix(7,KL)  = D(K,3)*P(L,1) 
      matrix(8,KL)  = D(K,3)*P(L,2) 
      matrix(9,KL)  = D(K,3)*P(L,3) 
      matrix(10,KL) = D(K,4)*P(L,1) 
      matrix(11,KL) = D(K,4)*P(L,2) 
      matrix(12,KL) = D(K,4)*P(L,3) 
      matrix(13,KL) = D(K,5)*P(L,1) 
      matrix(14,KL) = D(K,5)*P(L,2) 
   50 matrix(15,KL) = D(K,5)*P(L,3) 
!     D-D                                                               
      DO 60 K=1,5 
      KL        = INDX(K+4)+K+4 
      matrix(1,KL)  = D(K,1)*D(K,1) 
      matrix(2,KL)  = D(K,1)*D(K,2) 
      matrix(3,KL)  = D(K,2)*D(K,2) 
      matrix(4,KL)  = D(K,1)*D(K,3) 
      matrix(5,KL)  = D(K,2)*D(K,3) 
      matrix(6,KL)  = D(K,3)*D(K,3) 
      matrix(7,KL)  = D(K,1)*D(K,4) 
      matrix(8,KL)  = D(K,2)*D(K,4) 
      matrix(9,KL)  = D(K,3)*D(K,4) 
      matrix(10,KL) = D(K,4)*D(K,4) 
      matrix(11,KL) = D(K,1)*D(K,5) 
      matrix(12,KL) = D(K,2)*D(K,5) 
      matrix(13,KL) = D(K,3)*D(K,5) 
      matrix(14,KL) = D(K,4)*D(K,5) 
   60 matrix(15,KL) = D(K,5)*D(K,5) 
      DO 70 K=2,5 
      DO 70 L=1,K-1 
      KL        = INDX(K+4)+L+4 
      matrix(1,KL)  = D(K,1)*D(L,1)*two 
      matrix(2,KL)  = D(K,1)*D(L,2)+D(K,2)*D(L,1) 
      matrix(3,KL)  = D(K,2)*D(L,2)*two
      matrix(4,KL)  = D(K,1)*D(L,3)+D(K,3)*D(L,1) 
      matrix(5,KL)  = D(K,2)*D(L,3)+D(K,3)*D(L,2) 
      matrix(6,KL)  = D(K,3)*D(L,3)*two 
      matrix(7,KL)  = D(K,1)*D(L,4)+D(K,4)*D(L,1) 
      matrix(8,KL)  = D(K,2)*D(L,4)+D(K,4)*D(L,2) 
      matrix(9,KL)  = D(K,3)*D(L,4)+D(K,4)*D(L,3) 
      matrix(10,KL) = D(K,4)*D(L,4)*two 
      matrix(11,KL) = D(K,1)*D(L,5)+D(K,5)*D(L,1) 
      matrix(12,KL) = D(K,2)*D(L,5)+D(K,5)*D(L,2) 
      matrix(13,KL) = D(K,3)*D(L,5)+D(K,5)*D(L,3) 
      matrix(14,KL) = D(K,4)*D(L,5)+D(K,5)*D(L,4) 
   70 matrix(15,KL) = D(K,5)*D(L,5)*two 
   
      end if ! hasDOribtal
   
      if (keepOld) matrix_saved=matrix
      RETURN 
                
end subroutine GenerateRotationMatrix


SUBROUTINE Rotate2Center2Electron(W,LIMIJ,LImkl,rotationMatrix) 
!     *                                                                 
!     TRANSFORMATION OF TWO-ELECTRON TWO-CENTER INTEGRALS (SPD BASIS)   
!     FROM THE LOCAL TO THE MOLECULAR COORDINATE SYSTEM.                
!     A TWO-STEP PROCEDURE IS USED.                                     
!     *             

      use constants          , only : zero
      use ElementOrbitalIndex, only : Is2CenterIntegralZero
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
!SET  IROTD : CHOICE OF METHOD FOR ROTATING TWO-ELECTRON SPD INTEGRALS. 
!SET  IROTV : DEFAULT VALUE FOR OPTION IROTD.                           
!SET  IROTV : MACHINE-SPECIFIC VALUE WHICH MAY BE REDEFINED.            
!SET  IROTV = 1  UNROLLED SCALAR CODE.                                  
!SET  IROTV = 2  VECTORIZED CODE WITH SMALL LOOPS.                      
!SET  IROTV = 3  STRAIGHT MATRIX MULTIPLICATIONS WITH MANY ZERO VALUES. 
      PARAMETER (IROTV=1) 
!RAY  PARAMETER (IROTV=3)                                               
      PARAMETER (LIMY=45) 
      PARAMETER (LIMV=45) 
      PARAMETER (VSMALL=1.0D-10) 
      LOGICAL R2CENT(45,45) 
      DIMENSION W(LImkl,LIMIJ) 
      DIMENSION V(LIMV,LIMV) 
      DIMENSION YM(LIMY,LIMY) 
      DIMENSION rotationMatrix(15,45) 
      DIMENSION WKL(45) 
      DIMENSION MET(45) 
      DIMENSION META(6),METB(6) 
      DIMENSION METI(15,6) 
      DIMENSION IKLV(45) 
!     MET ASSIGNS A PAIR TYPE (1 SS, 2 PS, 3 PP, 4 DS, 5 DP, 6 DD)      
!     MET(KL) TO EACH STANDARD PAIR INDEX KL=INDX(K)+L.                 
      DATA MET /  1, 2, 3, 2, 3, 3, 2, 3, 3, 3, 4, 5, 5, 5, 6,          &
     &            4, 5, 5, 5, 6, 6, 4, 5, 5, 5, 6, 6, 6, 4, 5,          &
     &            5, 5, 6, 6, 6, 6, 4, 5, 5, 5, 6, 6, 6, 6, 6/          
!     META AND METB DEFINE THE FIRST AND LAST ADDRESS OF A SPECIAL      
!     PAIR INDEX (KLV) FOR EACH PAIR TYPE (1-1 SS, 2-4 PS, 5-10 PP,     
!     11-15 DS, 16-30 DP, 31-45 DD). WITHIN EACH SUCH RANGE, THE PAIR   
!     INDICES ARE ORDERED IN THE STANDARD MANNER.                       
      DATA META/  1, 2, 5,11,16,31 / 
      DATA METB/  1, 4,10,15,30,45 / 
!     METI(I,mkl) DEFINES THE STANDARD PAIR INDICES IJ BELONGING        
!     TO A GIVEN VALUE OF mkl=MET(KL).                                  
      DATA METI/  1, 14*0,                                              &
     &            2, 4, 7, 12*0,                                        &
     &            3, 5, 6, 8, 9,10,  9*0,                               &
     &           11,16,22,29,37,    10*0,                               &
     &           12,13,14,17,18,19,23,24,25,30,31,32,38,39,40,          &
     &           15,20,21,26,27,28,33,34,35,36,41,42,43,44,45/          
!     IKLV REFERS A GIVEN SPECIAL PAIR INDEX (KLV) - USED FOR STORING   
!     INTERMEDIATE INTEGRALS - TO THE STANDARD PAIR INDEX KL=INDX(K)+L  
!     SUCH THAT KL=IKLV(KLV).                                           
      DATA IKLV/  1, 2, 4, 7, 3, 5, 6, 8, 9,10,11,16,22,29,37,          &
     &           12,13,14,17,18,19,23,24,25,30,31,32,38,39,40,          &
     &           15,20,21,26,27,28,33,34,35,36,41,42,43,44,45/          
!     *                                                                 
! *** INITIALIZATION.                                                   
!     *                 

      call Is2CenterIntegralZero(r2cent)
      DO 10 IJ=1,LIMIJ 
      DO 10 KL=1,LImkl 
   10 V(KL,IJ) = ZERO 
      IROTD    = IROTV 
      IF(IROTD.EQ.2) GO TO 300 
      IF(IROTD.GT.2) GO TO 500 
!     *                                                                 
! *** FIRST SECTION. UNROLLED SCALAR CODE.                              
!     *                                                                 
! *** DEFINITION OF INDICES FOR INTERMEDIATE INTEGRALS V(KL,IJ),        
!     FINAL INTEGRALS W(KL,IJ), AND TRANSFORMATION MATRIX rotationMatrix(M,KL).     
!     IJ    STANDARD PAIR INDEX.                                        
!     KL    STANDARD PAIR INDEX.                                        
!     M     CONSECUTIVE INDEX (UP TO 1 FOR KL=SS, 3 FOR KL=PS,          
!           6 FOR KL=PP, 5 FOR KL=DS, 15 FOR KL=PD AND KL=DD).          
! *** FIRST STEP - TRANSFORMATION FOR INDICES KL IN (KL,IJ).            
!     LOOP OVER OUTER TWO INDICES IJ (NOT TRANSFORMED).                 
!$DIR SCALAR                                                            
!VDIR NOVECTOR                                                          
!VOCL LOOP,SCALAR                                                       
      DO 180 IJ=1,LIMIJ 
!     LOOP OVER INNER TWO INDICES KL (TRANSFORMED).                     
!DIR$ NEXTSCALAR                                                        
!$DIR SCALAR                                                            
!VDIR NOVECTOR                                                          
!VOCL LOOP,SCALAR                                                       
      DO 170 KL=1,LImkl 
      IF(R2CENT(KL,IJ)) THEN 
         WREPP    = W(KL,IJ) 
         W(KL,IJ) = ZERO 
         GO TO (110,120,130,140,150,160), MET(KL) 
  110    V(1,IJ)  =  WREPP 
         GO TO 170 
  120    CONTINUE 
         V(2,IJ)  = V(2,IJ)  + rotationMatrix(1,KL) * WREPP 
         V(4,IJ)  = V(4,IJ)  + rotationMatrix(2,KL) * WREPP 
         V(7,IJ)  = V(7,IJ)  + rotationMatrix(3,KL) * WREPP 
         GO TO 170 
  130    CONTINUE 
         V(3,IJ)  = V(3,IJ)  + rotationMatrix(1,KL) * WREPP 
         V(5,IJ)  = V(5,IJ)  + rotationMatrix(2,KL) * WREPP 
         V(6,IJ)  = V(6,IJ)  + rotationMatrix(3,KL) * WREPP 
         V(8,IJ)  = V(8,IJ)  + rotationMatrix(4,KL) * WREPP 
         V(9,IJ)  = V(9,IJ)  + rotationMatrix(5,KL) * WREPP 
         V(10,IJ) = V(10,IJ) + rotationMatrix(6,KL) * WREPP 
         GO TO 170 
  140    CONTINUE 
         V(11,IJ) = V(11,IJ) + rotationMatrix(1,KL) * WREPP 
         V(16,IJ) = V(16,IJ) + rotationMatrix(2,KL) * WREPP 
         V(22,IJ) = V(22,IJ) + rotationMatrix(3,KL) * WREPP 
         V(29,IJ) = V(29,IJ) + rotationMatrix(4,KL) * WREPP 
         V(37,IJ) = V(37,IJ) + rotationMatrix(5,KL) * WREPP 
         GO TO 170 
  150    CONTINUE 
         V(12,IJ) = V(12,IJ) + rotationMatrix( 1,KL) * WREPP 
         V(13,IJ) = V(13,IJ) + rotationMatrix( 2,KL) * WREPP 
         V(14,IJ) = V(14,IJ) + rotationMatrix( 3,KL) * WREPP 
         V(17,IJ) = V(17,IJ) + rotationMatrix( 4,KL) * WREPP 
         V(18,IJ) = V(18,IJ) + rotationMatrix( 5,KL) * WREPP 
         V(19,IJ) = V(19,IJ) + rotationMatrix( 6,KL) * WREPP 
         V(23,IJ) = V(23,IJ) + rotationMatrix( 7,KL) * WREPP 
         V(24,IJ) = V(24,IJ) + rotationMatrix( 8,KL) * WREPP 
         V(25,IJ) = V(25,IJ) + rotationMatrix( 9,KL) * WREPP 
         V(30,IJ) = V(30,IJ) + rotationMatrix(10,KL) * WREPP 
         V(31,IJ) = V(31,IJ) + rotationMatrix(11,KL) * WREPP 
         V(32,IJ) = V(32,IJ) + rotationMatrix(12,KL) * WREPP 
         V(38,IJ) = V(38,IJ) + rotationMatrix(13,KL) * WREPP 
         V(39,IJ) = V(39,IJ) + rotationMatrix(14,KL) * WREPP 
         V(40,IJ) = V(40,IJ) + rotationMatrix(15,KL) * WREPP 
         GO TO 170 
  160    CONTINUE 
         V(15,IJ) = V(15,IJ) + rotationMatrix( 1,KL) * WREPP 
         V(20,IJ) = V(20,IJ) + rotationMatrix( 2,KL) * WREPP 
         V(21,IJ) = V(21,IJ) + rotationMatrix( 3,KL) * WREPP 
         V(26,IJ) = V(26,IJ) + rotationMatrix( 4,KL) * WREPP 
         V(27,IJ) = V(27,IJ) + rotationMatrix( 5,KL) * WREPP 
         V(28,IJ) = V(28,IJ) + rotationMatrix( 6,KL) * WREPP 
         V(33,IJ) = V(33,IJ) + rotationMatrix( 7,KL) * WREPP 
         V(34,IJ) = V(34,IJ) + rotationMatrix( 8,KL) * WREPP 
         V(35,IJ) = V(35,IJ) + rotationMatrix( 9,KL) * WREPP 
         V(36,IJ) = V(36,IJ) + rotationMatrix(10,KL) * WREPP 
         V(41,IJ) = V(41,IJ) + rotationMatrix(11,KL) * WREPP 
         V(42,IJ) = V(42,IJ) + rotationMatrix(12,KL) * WREPP 
         V(43,IJ) = V(43,IJ) + rotationMatrix(13,KL) * WREPP 
         V(44,IJ) = V(44,IJ) + rotationMatrix(14,KL) * WREPP 
         V(45,IJ) = V(45,IJ) + rotationMatrix(15,KL) * WREPP 
      ENDIF 
  170 END DO 
  180 END DO 
! *** SECOND STEP - TRANSFORMATION FOR OUTER TWO INDICES.               
!     LOOP OVER OUTER TWO INDICES IJ (TRANSFORMED).                     
!$DIR SCALAR                                                            
!VDIR NOVECTOR                                                          
!VOCL LOOP,SCALAR                                                       
      DO 280 IJ=1,LIMIJ 
      MIJ    = MET(IJ) 
!     LOOP OVER INNER TWO INDICES KL (NOT TRANSFORMED).                 
!DIR$ NEXTSCALAR                                                        
!$DIR SCALAR                                                            
!VDIR NOVECTOR                                                          
!VOCL LOOP,SCALAR                                                       
      DO 270 KL=1,LImkl 
      IF(ABS(V(KL,IJ)).GT.VSMALL) THEN 
         WREPP  = V(KL,IJ) 
         GO TO (210,220,230,240,250,260), MIJ 
  210    W(KL,1)  =  WREPP 
         GO TO 270 
  220    CONTINUE 
         W(KL,2)  = W(KL,2)  + rotationMatrix(1,IJ) * WREPP 
         W(KL,4)  = W(KL,4)  + rotationMatrix(2,IJ) * WREPP 
         W(KL,7)  = W(KL,7)  + rotationMatrix(3,IJ) * WREPP 
         GO TO 270 
  230    CONTINUE 
         W(KL,3)  = W(KL,3)  + rotationMatrix(1,IJ) * WREPP 
         W(KL,5)  = W(KL,5)  + rotationMatrix(2,IJ) * WREPP 
         W(KL,6)  = W(KL,6)  + rotationMatrix(3,IJ) * WREPP 
         W(KL,8)  = W(KL,8)  + rotationMatrix(4,IJ) * WREPP 
         W(KL,9)  = W(KL,9)  + rotationMatrix(5,IJ) * WREPP 
         W(KL,10) = W(KL,10) + rotationMatrix(6,IJ) * WREPP 
         GO TO 270 
  240    CONTINUE 
         W(KL,11) = W(KL,11) + rotationMatrix(1,IJ) * WREPP 
         W(KL,16) = W(KL,16) + rotationMatrix(2,IJ) * WREPP 
         W(KL,22) = W(KL,22) + rotationMatrix(3,IJ) * WREPP 
         W(KL,29) = W(KL,29) + rotationMatrix(4,IJ) * WREPP 
         W(KL,37) = W(KL,37) + rotationMatrix(5,IJ) * WREPP 
         GO TO 270 
  250    CONTINUE 
         W(KL,12) = W(KL,12) + rotationMatrix( 1,IJ) * WREPP 
         W(KL,13) = W(KL,13) + rotationMatrix( 2,IJ) * WREPP 
         W(KL,14) = W(KL,14) + rotationMatrix( 3,IJ) * WREPP 
         W(KL,17) = W(KL,17) + rotationMatrix( 4,IJ) * WREPP 
         W(KL,18) = W(KL,18) + rotationMatrix( 5,IJ) * WREPP 
         W(KL,19) = W(KL,19) + rotationMatrix( 6,IJ) * WREPP 
         W(KL,23) = W(KL,23) + rotationMatrix( 7,IJ) * WREPP 
         W(KL,24) = W(KL,24) + rotationMatrix( 8,IJ) * WREPP 
         W(KL,25) = W(KL,25) + rotationMatrix( 9,IJ) * WREPP 
         W(KL,30) = W(KL,30) + rotationMatrix(10,IJ) * WREPP 
         W(KL,31) = W(KL,31) + rotationMatrix(11,IJ) * WREPP 
         W(KL,32) = W(KL,32) + rotationMatrix(12,IJ) * WREPP 
         W(KL,38) = W(KL,38) + rotationMatrix(13,IJ) * WREPP 
         W(KL,39) = W(KL,39) + rotationMatrix(14,IJ) * WREPP 
         W(KL,40) = W(KL,40) + rotationMatrix(15,IJ) * WREPP 
         GO TO 270 
  260    CONTINUE 
         W(KL,15) = W(KL,15) + rotationMatrix( 1,IJ) * WREPP 
         W(KL,20) = W(KL,20) + rotationMatrix( 2,IJ) * WREPP 
         W(KL,21) = W(KL,21) + rotationMatrix( 3,IJ) * WREPP 
         W(KL,26) = W(KL,26) + rotationMatrix( 4,IJ) * WREPP 
         W(KL,27) = W(KL,27) + rotationMatrix( 5,IJ) * WREPP 
         W(KL,28) = W(KL,28) + rotationMatrix( 6,IJ) * WREPP 
         W(KL,33) = W(KL,33) + rotationMatrix( 7,IJ) * WREPP 
         W(KL,34) = W(KL,34) + rotationMatrix( 8,IJ) * WREPP 
         W(KL,35) = W(KL,35) + rotationMatrix( 9,IJ) * WREPP 
         W(KL,36) = W(KL,36) + rotationMatrix(10,IJ) * WREPP 
         W(KL,41) = W(KL,41) + rotationMatrix(11,IJ) * WREPP 
         W(KL,42) = W(KL,42) + rotationMatrix(12,IJ) * WREPP 
         W(KL,43) = W(KL,43) + rotationMatrix(13,IJ) * WREPP 
         W(KL,44) = W(KL,44) + rotationMatrix(14,IJ) * WREPP 
         W(KL,45) = W(KL,45) + rotationMatrix(15,IJ) * WREPP 
      ENDIF 
  270 END DO 
  280 END DO 
      RETURN 
!     *                                                                 
! *** SECOND SECTION. ALTERNATIVE VECTORIZED CODE.                      
!     *                                                                 
! *** DEFINITION OF INDICES FOR INTERMEDIATE INTEGRALS V(KLV,IJ),       
!     FINAL INTEGRALS W(KL,IJ), AND TRANSFORMATION MATRIX rotationMatrix(M,KL).     
!     IJ    STANDARD PAIR INDEX.                                        
!     KL    STANDARD PAIR INDEX.                                        
!     KLV   SPECIAL INDEX SUCH THAT THE INTERMEDIATE ARRAY IS           
!           ACCESSED SEQUENTIALLY WITH STRIDE 1 (SEE BELOW).            
!           THE CORRESPONDENCE OF KLV WITH THE STANDARD PAIR            
!           INDEX IS DEFINED IN IKLV (SEE DATA STATEMENT).              
!     M     CONSECUTIVE INDEX (UP TO 1 FOR KL=SS, 3 FOR KL=PS,          
!           6 FOR KL=PP, 5 FOR KL=DS, 15 FOR KL=PD AND KL=DD).          
!                                                                       
!     KLV   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15                 
!     K     1  2  3  4  2  3  3  4  4  4  5  6  7  8  9                 
!     L     1  1  1  1  2  2  3  2  3  4  1  1  1  1  1                 
!     TYPE SS PS PS PS PP PP PP PP PP PP DS DS DS DS DS                 
!                                                                       
!     KLV  16 17 18 19 20 21 22 23 24 25 26 27 28 29 30                 
!     K     5  5  5  6  6  6  7  7  7  8  8  8  9  9  9                 
!     L     2  3  4  2  3  4  2  3  4  2  3  4  2  3  4                 
!     TYPE DP DP DP DP DP DP DP DP DP DP DP DP DP DP DP                 
!                                                                       
!     KLV  31 32 33 34 35 36 37 38 39 40 41 42 43 44 45                 
!     K     5  6  6  7  7  7  8  8  8  8  9  9  9  9  9                 
!     L     5  5  6  5  6  7  5  6  7  8  5  6  7  8  9                 
!     TYPE DD DD DD DD DD DD DD DD DD DD DD DD DD DD DD                 
!                                                                       
!     TYPE SS PS PP DS DP DD                                            
!     MET   1  2  3  4  5  6    LABEL FOR GIVEN TYPE                    
!     META  1  2  5 11 16 31    FIRST KLV VALUE FOR GIVEN TYPE          
!     METB  1  4 10 15 30 45    LAST  KLV VALUE FOR GIVEN TYPE          
!                                                                       
!     *                                                                 
! *** FIRST STEP - TRANSFORMATION FOR INDICES KL IN (KL,IJ).            
  300 CONTINUE 
!     LOOP OVER OUTER TWO INDICES IJ (NOT TRANSFORMED).                 
      DO 330 IJ=1,LIMIJ 
!     LOOP OVER INNER TWO INDICES KL (TRANSFORMED).                     
      DO 320 KL=1,LImkl 
      IF(R2CENT(KL,IJ)) THEN 
         WREPP  = W(KL,IJ) 
         mkl    = MET(KL) 
         KLA    = META(mkl) 
         KLB    = METB(mkl) 
         DO 310 KLV=KLA,KLB 
  310    V(KLV,IJ) = V(KLV,IJ) + rotationMatrix(KLV-KLA+1,KL)*WREPP 
      ENDIF 
  320 END DO 
  330 END DO 
! *** SECOND STEP - TRANSFORMATION FOR OUTER TWO INDICES.               
!     LOOP OVER INNER TWO INDICES KL (NOT TRANSFORMED).                 
      DO 450 KLV=1,LImkl 
      KL     = IKLV(KLV) 
!     INITIALIZE BUFFER.                                                
      DO 410 IJ=1,LIMIJ 
  410 WKL(IJ)= ZERO 
!     LOOP OVER OUTER TWO INDICES IJ (TRANSFORMED).                     
      DO 430 IJ=1,LIMIJ 
      IF(ABS(V(KLV,IJ)).GT.VSMALL) THEN 
         MIJ    = MET(IJ) 
         MIJA   = META(MIJ) 
         MIJB   = METB(MIJ) 
         WREPP  = V(KLV,IJ) 
         DO 420 M=MIJA,MIJB 
  420    WKL(M) = WKL(M) + rotationMatrix(M-MIJA+1,IJ)*WREPP 
      ENDIF 
  430 END DO 
!     TRANSFER RESULTS FROM BUFFER TO W(KL,IJ).                         
      DO 440 IJ=1,LIMIJ 
  440 W(KL,IKLV(IJ)) = WKL(IJ) 
  450 END DO 
      RETURN 
!     *                                                                 
! *** THIRD SECTION. STRAIGHT MATRIX MANIPULATIONS.                     
!     W(MOL,TRANSPOSE) = Y(TRANSPOSE) * W(LOC,TRANSPOSE)  * Y           
!     THE INTEGRALS IN THE MOLECULAR COORDINATE SYSTEM ARE NEEDED       
!     AS A TRANSPOSED MATRIX. THEREFORE THE INTEGRALS IN THE LOCAL      
!     COORDINATE SYSTEM ARE ALSO PROVIDED AS A TRANSPOSED MATRIX.       
!     THE TRANSFORMATION MATRIX YM IS GENERATED IN TRANSPOSED FORM,     
!     Y = YM(TRANSPOSE), Y(TRANSPOSE) = YM,                             
!     FROM THE ELEMENTS OF ARRAY rotationMatrix.                                    
!     *                                                                 
!     STANDARD PAIR INDICES ARE USED THROUGHOUT. THEY REFER TO THE      
!     ORDER OF ATOMIC ORBITALS IN INTEGRAL EVALUATION FOR W(LOC),       
!     AND TO THAT IN THE SCF CALCULATION FOR W(MOL). THE MATRIX Y       
!     EFFECTS THE SWITCH FROM ONE ORDERING SCHEME TO THE OTHER.         
!     *                                                                 
  500 CONTINUE 
! *** INITIALIZATION.                                                   
      LIMYM  = MAX(LIMIJ,LImkl) 
      DO 510 IJ=1,LIMYM 
      DO 510 KL=1,LIMYM 
  510 YM(KL,IJ) = ZERO 
      YM( 1, 1) = ONE 
      DO 530 KL=2,LIMYM 
      mkl    = MET(KL) 
      NKL    = METB(mkl)-META(mkl)+1 
      DO 520 I=1,NKL 
  520 YM(METI(I,mkl),KL) = rotationMatrix(I,KL) 
  530 END DO 
! *** FIRST IMPLEMENTATION OF THE TWO-STEP TRANSFORMATION.              
!     V(INTERMEDIATE)  = YM * W(LOC,TRANSPOSE)                          
!     W(MOL,TRANSPOSE) = V(INTERMEDIATE) * YM(TRANSPOSE).               
!     FORTRAN VERSION OF DGEMM CHECKS ON SPARSITY OF SECOND MATRIX.     
      IF(IROTD.EQ.3) THEN 
         CALL DGEMM ('N','N',LImkl,LIMIJ,LImkl,ONE,YM,LIMY,W,LImkl,     &
     &               ZERO,V,LIMV)                                       
         CALL DGEMM ('N','T',LImkl,LIMIJ,LIMIJ,ONE,V,LIMV,YM,LIMY,      &
     &               ZERO,W,LImkl)                                      
      ENDIF 
! *** SECOND IMPLEMENTATION OF THE TWO-STEP TRANSFORMATION.             
!     V(INTERMEDIATE)  = W(LOC,TRANSPOSE) * YM(TRANSPOSE).              
!     W(MOL,TRANSPOSE) = YM * V(INTERMEDIATE).                          
!     FORTRAN VERSION OF DGEMM CHECKS ON SPARSITY OF SECOND MATRIX.     
      IF(IROTD.EQ.4) THEN 
         CALL DGEMM ('N','T',LImkl,LIMIJ,LIMIJ,ONE,W,LImkl,YM,LIMY,     &
     &               ZERO,V,LIMV)                                       
         CALL DGEMM ('N','N',LImkl,LIMIJ,LImkl,ONE,YM,LIMY,V,LIMV,      &
     &               ZERO,W,LImkl)                                      
      ENDIF 
      RETURN 
end subroutine Rotate2Center2Electron    


SUBROUTINE Rotate1Elec(IA,JA,IORBS,JORBS,T,matrix,norbs, H, matsize) 
!     *                                                                 
!     THIS ROUTINE TRANSFORMS two-CENTER ONE-ELECTRON INTEGRALS FROM    
!     LOCAL TO MOLECULAR COORDINATES, AND INCLUDES THEM IN H(LM4).      
!     USEFUL FOR RESONANCE INTEGRALS AND FOR OVERLAP INTEGRALS.         
!     *                                                                 
!     NOTATION. I=INPUT, O=OUTPUT.                                      
!     IA        INDEX OF FIRST ORBITAL AT ATOM I (I).                   
!     JA        INDEX OF FIRST ORBITAL AT ATOM J (I).                   
!     IORBS     NUMBER OF ORBITALS AT ATOM I (I).                       
!     JORBS     NUMBER OF ORBITALS AT ATOM J (I).                       
!     T(14)     LOCAL two-CENTER ONE-ELECTRON INTEGRALS (I).            
!     matrix()      PRECOMPUTED COMBINATION OF ROTATION MATRIX ELEMENTS (I).
!     H(LM4)    ONE-ELECTRON MATRIX IN MOLECULAR COORDINATES (O).       
!     *                                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      integer, intent(in) :: norbs, matsize
      _REAL_, intent(inout):: T(*),H(matsize),matrix(15,45)

! local
      integer ::indx(1:norbs)
      _REAL_:: HDD(15)                                     
 
!  
!     
      do i=1, norbs
         indx(i)=i*(i-1)/2
      end do


! *** SECTION FOR AN SP-BASIS.                                          
!     S(I)-S(J)                                                         
      IS     = INDX(IA)+JA 
      H(IS)  = T(1) 
      IF(IORBS.EQ.1 .AND. JORBS.EQ.1) RETURN 
!     S(I)-P(J)                                                         
      IF(JORBS.GE.4) THEN 
         H(IS+1) = T(2)*matrix(1,2) 
         H(IS+2) = T(2)*matrix(2,2) 
         H(IS+3) = T(2)*matrix(3,2) 
      ENDIF 
!     P(I)-S(J)                                                         
      IF(IORBS.GE.4) THEN 
         IX      = INDX(IA+1)+JA 
         IY      = INDX(IA+2)+JA 
         IZ      = INDX(IA+3)+JA 
         H(IX)   = T(3)*matrix(1,2) 
         H(IY)   = T(3)*matrix(2,2) 
         H(IZ)   = T(3)*matrix(3,2) 
!     P(I)-P(J).                                                        
         IF(JORBS.GE.4) THEN 
            T45     = T(4)-T(5) 
            H(IX+1) = matrix(1,3)*T45+T(5) 
            H(IX+2) = matrix(2,3)*T45 
            H(IX+3) = matrix(4,3)*T45 
            H(IY+1) = H(IX+2) 
            H(IY+2) = matrix(3,3)*T45+T(5) 
            H(IY+3) = matrix(5,3)*T45 
            H(IZ+1) = H(IX+3) 
            H(IZ+2) = H(IY+3) 
            H(IZ+3) = matrix(6,3)*T45+T(5) 
         ENDIF 
      ENDIF 
! *** SECTION INVOLVING D ORBITALS.                                     
!     D(I)-S(J)                                                         
      IF(IORBS.GE.9) THEN 
         DO 10 I=1,5 
         M      = INDX(IA+3+I)+JA 
         H(M)   = T(6)*matrix(I,11) 
   10    CONTINUE 
!     D(I)-P(J)                                                         
         IF(JORBS.GE.4) THEN 
            IJ     = 0 
            DO 30 I=1,5 
            M      = INDX(IA+3+I)+JA 
            DO 20 J=1,3 
            IJ     = IJ+1 
            H(M+J) = T(8) * matrix(IJ,12)                                   &
     &              +T(10)*(matrix(IJ,18)+matrix(IJ,25))                        
   20       CONTINUE 
   30       CONTINUE 
         ENDIF 
      ENDIF 
!     S(I)-D(J)                                                         
      IF(JORBS.GE.9) THEN 
         M      = INDX(IA)+JA+3 
         DO 40 I=1,5 
         H(M+I) = T(7)*matrix(I,11) 
   40    CONTINUE 
         IF(IORBS.GE.4) THEN 
!     P(I)-D(J)                                                         
            DO 60 I=1,3 
            M      = INDX(IA+I)+JA+3 
            DO 50 J=1,5 
            IJ     = 3*(J-1)+I 
            H(M+J) = T(9) * matrix(IJ,12)                                   &
     &              +T(11)*(matrix(IJ,18)+matrix(IJ,25))                        
   50       CONTINUE 
   60       CONTINUE 
!     D(I)-D(J)                                                         
            IF(IORBS.GE.9) THEN 
               DO 70 I=1,15 
               HDD(I) = T(12)* matrix(I,15)                                 &
     &                 +T(13)*(matrix(I,21)+matrix(I,28))                       &
     &                 +T(14)*(matrix(I,36)+matrix(I,45))                       
   70          CONTINUE 
               DO 90 I=1,5 
               M      = INDX(IA+3+I)+JA+3 
               DO 80 J=1,5 
               IF(I.GE.J) THEN 
                  IJ  = INDX(I)+J 
               ELSE 
                  IJ  = INDX(J)+I 
               ENDIF 
               H(M+J) = HDD(IJ) 
   80          CONTINUE 
   90          CONTINUE 
            ENDIF 
         ENDIF 
      ENDIF 
      RETURN 
END subroutine Rotate1Elec                                          

SUBROUTINE RotateCore(IA,JA,IORBS,JORBS,IP,JP,CORE,matrix,H)
!C     *
!C     THIS ROUTINE TRANSFORMS THE CORE ELECTRON ATTRACTION INTEGRALS
!C     FROM LOCAL TO MOLECULAR COORDINATES, AND INCLUDES THEM IN THE
!C     CORE HAMILTONIAN.
!C     *
!C     NOTATION. I=INPUT, O=OUTPUT, S=SCRATCH.
!C     IA        INDEX OF FIRST BASIS ORBITAL AT ATOM I (I).
!C     JA        INDEX OF FIRST BASIS ORBITAL AT ATOM J (I).
!C     IORBS     NUMBER OF BASIS ORBITALS AT ATOM I (I).
!C     JORBS     NUMBER OF BASIS ORBITALS AT ATOM J (I).
!C     IP        INDEX FOR (S,S) OF ATOM I IN LINEAR ARRAY H (I,S).
!C     JP        INDEX FOR (S,S) OF ATOM J IN LINEAR ARRAY H (I,S).
!C     CORE()    LOCAL CORE ELECTRON ATTRACTION INTEGRALS (I).
!C     YY()      PRECOMBINED ELEMENTS OF ROTATION MATRIX (I).
!C     H(LMH)    CORE HAMILTONIAN MATRIX (O).
!C     *
!C     DEPENDING ON THE INPUT DATA IN THE ARGUMENT LIST, THE INTEGRALS
!C     ARE INCLUDED EITHER IN THE FULL CORE HAMILTONIAN H(LM4) OR IN
!C     THE ONE-CENTER PART H(LM6) OF THE CORE HAMILTONIAN.
!C
!C     ARGUMENT       FIRST CASE              SECOND CASE
!C     IA             NFIRST(I)               1
!C     JA             NFIRST(J)               1
!C     IP             INDX(IA)+IA             NW(I)
!C     JP             INDX(JA)+JA             NW(J)
!C     LMH            LM4                      LM6
!C     *
!C     FOR IORBS.LE.0 OR JORBS.LE.0, THERE ARE NO BASIS FUNCTIONS AT
!C     ATOMS I OR J, RESPECTIVELY, AND HENCE NO CONTRIBUTIONS TO THE
!C     CORE HAMILTONIAN. THIS CASE MAY ARISE IN CALCULATIONS WITH
!C     EXTERNAL POINT CHARGES.
!C     *
      IMPLICIT none
      _REAL_, intent(in)::CORE(10,2),matrix(15,45)
      _REAL_, intent(inout)::H(*)
      integer, intent(in)::IA,JA,IORBS,JORBS,IP,JP

!local

      integer::i,j,k,l,is, ii,jj,kk
      integer::ix,iy,iz,idp,idd,id
      _REAL_::HPP(6),HDP(15),HDD(15)
      
!C$DIR SCALAR
!*VDIR NOVECTOR
!*VOCL LOOP,SCALAR
      DO 60 KK=1,2
      IF(KK.EQ.1) THEN
         IS  = IP
         K   = IA-1
         L   = IORBS
      ELSE
         IS  = JP
         K   = JA-1
         L   = JORBS
      ENDIF
      IF(L.LE.0) GO TO 60
!C     S-S
      H(IS)  = H(IS)+CORE(1,KK)
      IF(L.EQ.1) GO TO 60
!C     INTERMEDIATE RESULTS FOR P-P
      DO 10 I=1,6
      HPP(I) = CORE(3,KK)*matrix(I,3)+CORE(4,KK)*(matrix(I,6)+matrix(I,10))
   10 CONTINUE
!C     P-S
      IX     = IS+1+K
      IY     = IX+2+K
      IZ     = IY+3+K
      H(IX)  = H(IX)+CORE(2,KK)*matrix(1,2)
      H(IY)  = H(IY)+CORE(2,KK)*matrix(2,2)
      H(IZ)  = H(IZ)+CORE(2,KK)*matrix(3,2)
!C     P-P
      H(IX+1)= H(IX+1)+HPP(1)
      H(IY+1)= H(IY+1)+HPP(2)
      H(IY+2)= H(IY+2)+HPP(3)
      H(IZ+1)= H(IZ+1)+HPP(4)
      H(IZ+2)= H(IZ+2)+HPP(5)
      H(IZ+3)= H(IZ+3)+HPP(6)
!C     INTERMEDIATE RESULTS FOR D-P AND D-D
      IF(L.EQ.4) GO TO 60
      DO 20 I=1,15
      HDP(I) = CORE(6,KK)*matrix(I,12)+CORE( 8,KK)*(matrix(I,18)+matrix(I,25))
      HDD(I) = CORE(7,KK)*matrix(I,15)+CORE( 9,KK)*(matrix(I,21)+matrix(I,28)) &
                                 +CORE(10,KK)*(matrix(I,36)+matrix(I,45))
   20 CONTINUE
!C     D-S
      IDP    = 0
      IDD    = 0
      ID     = IZ+3+K
      DO 50 I=5,9
      H(ID+1)= H(ID+1)+CORE(5,KK)*matrix(I-4,11)
!C     D-P
      DO 30 J=2,4
      IDP    = IDP+1
      H(ID+J)= H(ID+J)+HDP(IDP)
   30 CONTINUE
!C     D-D
      DO 40 J=5,I
      IDD    = IDD+1
      H(ID+J)= H(ID+J)+HDD(IDD)
   40 CONTINUE
      ID     = ID+I+K
   50 CONTINUE
   60 CONTINUE
      RETURN
END subroutine RotateCore
      

end module Rotation
