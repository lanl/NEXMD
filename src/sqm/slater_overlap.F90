#include "dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! module contains
! 1. wrapper for Thiel's slater overlap routines
! 2. Slater-Condon parameters
!
! by: Taisung Lee (Rutgers, 2011)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module SlaterOverlap

use constants, only : zero, half, one, two, &
                      FC, logFC, CC, binomialCoefficient, &
                      AU_TO_EV

implicit none

public::GetSlaterOverlap, GetSlaterCondonParameter
private::SetupSlaterAuxiliary,CalculateOverlap, A, B, B0 

! Constants
!     THE ARRAY B0(I) CONTAINS THE B INTEGRALS FOR ZERO ARGUMENT. 
_REAL_, parameter :: B0(1:15)= &
     &      (/ 2.0D0,0.0D0,0.666666666666667D0,0.0D0,0.4D0,0.0D0,       &
     &         0.285714285714286D0,0.0D0,0.222222222222222D0,0.0D0,     &
     &         0.181818181818182D0,0.0D0,0.153846153846154D0,0.0D0,     &
     &         0.133333333333333D0 /)

                                                      

!local variables
_REAL_, save::A(15),B(15) 

contains

function GetSlaterOverlap(na, la, nb, lb, mm, zeta_a, zeta_b, rab) result (overlap)
    implicit none
    integer, intent(in)::na, la, nb, lb, mm
    _REAL_, intent(in)::zeta_a,zeta_b, rab
    _REAL_::overlap
    
    ! local variables
    integer, save::ntotal=-1
    _REAL_, save:: zeta_a_old=-1.0D99, zeta_b_old=-1.0D99
    _REAL_, save:: rab_old=-1.0D99
    _REAL_, parameter:: tolerance=1.D0-16
    
    logical::resetup
    
    if ((la.ge.na) .or. (lb.ge.nb)) then
       overlap=0.D0
       return
    end if
    resetup=.true.
    if (ntotal.eq. (na+nb)) then
       if (abs(rab_old/rab-1.D0)< tolerance) then
          if (abs(zeta_a_old/zeta_a-1.D0) < tolerance) then
            if (abs(zeta_b_old/zeta_b-1.D0) < tolerance) then
               resetup=.false.
            end if
          end if
       end if
    end if
    
    if (resetup) then
       ntotal=(na+nb)
       rab_old=rab
       zeta_a_old=zeta_a
       zeta_b_old=zeta_b
       call SetupSlaterAuxiliary(ntotal, zeta_a, zeta_b, rab)
    end if 
       
    overlap=CalculateOverlap (na,la,mm,nb,lb,zeta_a*rab,zeta_b*rab)
    
end function GetSlaterOverlap


SUBROUTINE SetupSlaterAuxiliary(N,SA,SB,RAB) 
!     *                                                                 
!     CALCULATION OF AUXILIARY INTEGRALS FOR STO OVERLAPS.              
!     *     
                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      implicit integer (I-N)
      !explicit type to satisfy PGI v8 compiler
      DOUBLE PRECISION BETPOW(17) 
      
                       
! *** INITIALIZATION.                                                   
      ALPHA  = 0.5D0*RAB*(SA+SB) 
      BETA   = 0.5D0*RAB*(SA-SB) 
! *** AUXILIARY A INTEGRALS FOR CALCULATION OF OVERLAPS.                
      C      = EXP(-ALPHA) 
      RALPHA = 1.0D0/ALPHA 
      A(1)   = C*RALPHA 
      DO 10 I=1,N 
      A(I+1) = (A(I)*I+C)*RALPHA 
   10 END DO 
! *** AUXILIARY B INTEGRALS FOR CALCULATION OF OVERLAPS.                
!     THE CODE IS VALID ONLY FOR N.LE.14, I.E. FOR OVERLAPS             
!     INVOLVING ORBITALS WITH MAIN QUANTUM NUMBERS UP TO 7.             
!     BRANCHING DEPENDING ON ABSOLUTE VALUE OF THE ARGUMENT.            
      ABSX   = ABS(BETA) 
!     ZERO ARGUMENT.                                                    
      IF(ABSX.LT.1.0D-06) THEN 
!DIR$ SHORTLOOP                                                         
!$DIR MAX_TRIPS(64)                                                     
!VDIR LOOPCNT=64                                                        
!VOCL LOOP,REPEAT(63)                                                   
         DO 20 I=1,N+1 
   20    B(I) = B0(I) 
         RETURN 
      ENDIF 
!     LARGE ARGUMENT.                                                   
      IF((ABSX.GT.0.5D0 .AND. N.LE.5) .OR.                              &
     &   (ABSX.GT.1.0D0 .AND. N.LE.7) .OR.                              &
     &   (ABSX.GT.2.0D0 .AND. N.LE.10).OR.                              &
     &    ABSX.GT.3.0D0) THEN                                           
         EXPX   = EXP(BETA) 
         EXPMX  = 1.0D0/EXPX 
         RX     = 1.0D0/BETA 
         B(1)   = (EXPX-EXPMX)*RX 
         DO 30 I=1,N 
         EXPX   = -EXPX 
   30    B(I+1)= (I*B(I)+EXPX-EXPMX)*RX 
         RETURN 
      ENDIF 
!     SMALL ARGUMENT.                                                   
      IF(ABSX.LE.0.5D0) THEN 
         LAST = 6 
      ELSE IF(ABSX.LE.1.0D0) THEN 
         LAST = 7 
      ELSE IF(ABSX.LE.2.0D0) THEN 
         LAST = 12 
      ELSE 
         LAST = 15 
      ENDIF 
      BETPOW(1) = 1.0D0 
      DO 40 M=1,LAST 
   40 BETPOW(M+1) = -BETA*BETPOW(M) 
      DO 60 I=1,N+1 
      Y      = 0.0D0 
      MA     = 1-MOD(I,2) 
      DO 50 M=MA,LAST,2 
   50 Y      = Y+BETPOW(M+1)/(FC(M+1)*(M+I)) 
   60 B(I)   = Y*2.0D0 
      RETURN 
END SUBROUTINE SetupSlaterAuxiliary                                          

FUNCTION CalculateOverlap (NA,LA,MM,NB,LB,ALPHA,BETA) 
!     *                                                                 
!     OVERLAP INTEGRALS BETWEEN SLATER TYPE ORBITALS.                   
!     *                                                                 
!     QUANTUM NUMBERS (NA,LA,MM) AND (NB,LB,MM).                        
!     NA AND NB MUST BE POSITIVE AND LESS THAN OR EQUAL TO 7.           
!     LA, LB AND ABS(MM) MUST BE LESS THAN OR EQUAL TO 5.               
!     FURTHER RESTRICTIONS ARE LA.LE.NA, LB.LE.NB,                      
!     MM.LE.LA, AND MM.LE.LB.                                           
!     *                                                                 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
      implicit integer (I-N)
      
!     DEFINE ADDRESSES FOR INDEX PAIRS (00,10,20,30,40,50,60,70).       
      integer,parameter:: IAD(1:8)= (/ 1,2,4,7,11,16,22,29 /) 
!     DEFINE BINOMIAL COEFFICIENTS (00,10,11,20,...,77).                
      integer,parameter:: IBINOM(1:36)= (/&
     &            1,1,1,1,2,1,1,3,3,1,1,4,6,4,1,1,5,10,10,5,1,          &
     &            1,6,15,20,15,6,1,1,7,21,35,35,21,7,1 /)                 

              
! *** INITIALIZATION.                                                   
      M      = ABS(MM) 
      NAB    = NA+NB+1 
      X      = 0.0D0 
! *** FIND A AND B INTEGRALS.                                           
!     P      = (ALPHA + BETA)*0.5                                       
!     PT     = (ALPHA - BETA)*0.5                                       
!     CALL AINTGS(A,P,NA+NB)                                            
!     CALL BINTGS(B,PT,NA+NB)                                           
! *** SECTION USED FOR OVERLAP INTEGRALS INVOLVING S FUNCTIONS.         
      IF((LA.GT.0).OR.(LB.GT.0)) GO TO 20 
      IADA   = IAD(NA+1) 
      IADB   = IAD(NB+1) 
      DO 10 I=0,NA 
      IBA    = IBINOM(IADA+I) 
      DO 10 J=0,NB 
      IBB    = IBA*IBINOM(IADB+J) 
      IF(MOD(J,2).EQ.1) IBB=-IBB 
      IJ     = I+J 
      X      = X+IBB*A(NAB-IJ)*B(IJ+1) 
   10 CONTINUE 
      SS     = X  * 0.5D0 
      SS     = SS * SQRT( ALPHA**(2*NA+1)*BETA**(2*NB+1)/               &
     &                    (FC(2*NA+1)*FC(2*NB+1)) )                     
      CalculateOverlap = SS 
      RETURN 
! *** SECTION USED FOR OVERLAP INTEGRALS INVOLVING P FUNCTIONS.         
! *** SPECIAL CASE M=0, S-P(SIGMA), P(SIGMA)-S, P(SIGMA)-P(SIGMA).      
   20 IF(LA.GT.1 .OR. LB.GT.1) GO TO 320 
      IF(M.GT.0) GO TO 220 
      IU     = MOD(LA,2) 
      IV     = MOD(LB,2) 
      NAMU   = NA-IU 
      NBMV   = NB-IV 
      IADNA  = IAD(NAMU+1) 
      IADNB  = IAD(NBMV+1) 
      DO 130 KC=0,IU 
      IC     = NAB-IU-IV+KC 
      JC     = 1+KC 
      DO 130 KD=0,IV 
      ID     = IC+KD 
      JD     = JC+KD 
      DO 130 KE=0,NAMU 
      IBE    = IBINOM(IADNA+KE) 
      IE     = ID-KE 
      JE     = JD+KE 
      DO 130 KF=0,NBMV 
      IBF    = IBE*IBINOM(IADNB+KF) 
      IF(MOD(KD+KF,2).EQ.1) IBF=-IBF 
      X      = X+IBF*A(IE-KF)*B(JE+KF) 
  130 CONTINUE 
      SS     = X  * SQRT( (2*LA+1)*(2*LB+1)*0.25D0 ) 
!     COMPUTE OVERLAP INTEGRAL FROM REDUCED OVERLAP INTEGRAL.           
      SS     = SS * SQRT( ALPHA**(2*NA+1)*BETA**(2*NB+1)/               &
     &                    (FC(2*NA+1)*FC(2*NB+1)) )                     
      IF(MOD(LB,2).EQ.1) SS=-SS 
      CalculateOverlap = SS 
      RETURN 
! *** SECTION USED FOR OVERLAP INTEGRALS INVOLVING P FUNCTIONS.         
! *** SPECIAL CASE LA=LB=M=1, P(PI)-P(PI).                              
  220 IADNA  = IAD(NA) 
      IADNB  = IAD(NB) 
      DO 230 KE=0,NA-1 
      IBE    = IBINOM(IADNA+KE) 
      IE     = NAB-KE 
      JE     = KE+1 
      DO 230 KF=0,NB-1 
      IBF    = IBE*IBINOM(IADNB+KF) 
      IF(MOD(KF,2).EQ.1) IBF=-IBF 
      I      = IE-KF 
      J      = JE+KF 
      X      = X+IBF*(A(I)*B(J)-A(I)*B(J+2)-A(I-2)*B(J)+A(I-2)*B(J+2)) 
  230 CONTINUE 
      SS     = X  * 0.75D0 
!     COMPUTE OVERLAP INTEGRAL FROM REDUCED OVERLAP INTEGRAL.           
      SS     = SS * SQRT( ALPHA**(2*NA+1)*BETA**(2*NB+1)/               &
     &                    (FC(2*NA+1)*FC(2*NB+1)) )                     
      IF(MOD(LB+MM,2).EQ.1) SS=-SS 
      CalculateOverlap = SS 
      RETURN 
! *** SECTION USED FOR OVERLAP INTEGRALS INVOLVING NON-S FUNCTIONS.     
! *** GENERAL CASE LA.GT.1 OR LB.GT.1, M.GE.0.                          
  320 LAM    = LA-M 
      LBM    = LB-M 
      IADA   = IAD(LA+1)+M 
      IADB   = IAD(LB+1)+M 
      IADM   = IAD(M+1) 
      IU1    = MOD(LAM,2) 
      IV1    = MOD(LBM,2) 
      IUC    = 0 
      DO 340 IU=IU1,LAM,2 
      IUC    = IUC+1 
      CU     = CC(IADA,IUC) 
      NAMU   = NA-M-IU 
      IADNA  = IAD(NAMU+1) 
      IADU   = IAD(IU+1) 
      IVC    = 0 
      DO 340 IV=IV1,LBM,2 
      IVC    = IVC+1 
      NBMV   = NB-M-IV 
      IADNB  = IAD(NBMV+1) 
      IADV   = IAD(IV+1) 
      SUM    = 0.0D0 
      DO 330 KC=0,IU 
      IBC    = IBINOM(IADU+KC) 
      IC     = NAB-IU-IV+KC 
      JC     = 1+KC 
      DO 330 KD=0,IV 
      IBD    = IBC*IBINOM(IADV+KD) 
      ID     = IC+KD 
      JD     = JC+KD 
      DO 330 KE=0,NAMU 
      IBE    = IBD*IBINOM(IADNA+KE) 
      IE     = ID-KE 
      JE     = JD+KE 
      DO 330 KF=0,NBMV 
      IBF    = IBE*IBINOM(IADNB+KF) 
      IFF    = IE-KF 
      JFF    = JE+KF 
      DO 330 KA=0,M 
      IBA    = IBF*IBINOM(IADM+KA) 
      I      = IFF-2*KA 
      DO 330 KB=0,M 
      IBB    = IBA*IBINOM(IADM+KB) 
      IF(MOD(KA+KB+KD+KF,2).EQ.1) IBB=-IBB 
      J      = JFF+2*KB 
      SUM    = SUM+IBB*A(I)*B(J) 
  330 CONTINUE 
      X      = X+SUM*CU*CC(IADB,IVC) 
  340 CONTINUE 
      SS     = X*(FC(M+2)/8.0D0)**2* SQRT( (2*LA+1)*FC(LA-M+1)*         &
     &         (2*LB+1)*FC(LB-M+1)/(4.0D0*FC(LA+M+1)*FC(LB+M+1)))       
!     COMPUTE OVERLAP INTEGRAL FROM REDUCED OVERLAP INTEGRAL.           
      SS     = SS * SQRT( ALPHA**(2*NA+1)*BETA**(2*NB+1)/               &
     &                    (FC(2*NA+1)*FC(2*NB+1)) )                     

       IF(MOD(LB+MM,2).EQ.1) SS=-SS 
      CalculateOverlap = SS 
      RETURN 
END function CalculateOverlap                                           

function GetSlaterCondonParameter(K,NA,EA,NB,EB,NC,EC,ND,ED) result(slaterCondon)
      
!C     *
!C     CALCULATE THE RADIAL PART OF ONE-CENTER TWO-ELECTRON INTEGRALS
!C     (SLATER-CONDON PARAMETER).
!c     K     - TYPE OF INTEGRAL, CAN BE EQUAL TO 0,1,2,3,4 IN SPD-BASIS
!C     NA,NB - PRINCIPLE QUANTUM NUMBER OF AO, ELECTRON 1
!C     EA,EB - EXPONENTS OF AO, ELECTRON 1
!C     NC,ND - PRINCIPLE QUANTUM NUMBER OF AO, ELECTRON 2
!C     EC,ED - EXPONENTS OF AO, ELECTRON 2
!C     *
     integer, intent(in)::K,NA,NB,NC,ND
     _REAL_, intent(in)::EA,EB,EC,ED
     _REAL_::slaterCondon
     
     !local
     integer::NAB, NCD, N, I, M, M1, M2
     _REAL_::EAB, ECD, E, C, S0, S1, S2, S3
     _REAL_::AEA, AEB, AEC, AED, AE, A2, ACD, AAB
!      COMMON
!     ./CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP
!     ./CONSTN/ ZERO,ONE,TWO,THREE,FOUR,PT5,PT25
!     ./PSC   / F(30),B(30,30)
      AEA    = LOG(EA)
      AEB    = LOG(EB)
      AEC    = LOG(EC)
      AED    = LOG(ED)
      NAB    = NA+NB
      NCD    = NC+ND
      ECD    = EC+ED
      EAB    = EA+EB
      E      = ECD+EAB
      N      = NAB+NCD
      AE     = LOG(E)
      A2     = LOG(TWO)
      ACD    = LOG(ECD)
      AAB    = LOG(EAB)
      C      = EXP(logFC(N)+NA*AEA+NB*AEB+NC*AEC+ND*AED &
                   +half*(AEA+AEB+AEC+AED)+A2*(N+2) &
                   -half*(logFC(2*NA+1)+logFC(2*NB+1)   &
                   +logFC(2*NC+1)+logFC(2*ND+1))-AE*N)
      C      = C*AU_TO_EV
      S0     = ONE/E
      S1     = ZERO
      S2     = ZERO
      M      = NCD-K
      DO 10 I=1,M
      S0     = S0*E/ECD
   10 S1     = S1+S0*(BinomialCoefficient(NCD-K,I)-BinomialCoefficient(NCD+K+1,I))/BinomialCoefficient(N,I)
      M1     = M+1
      M2     = NCD+K+1
      DO 20 I=M1,M2
      S0     = S0*E/ECD
   20 S2     = S2+S0*BinomialCoefficient(M2,I)/BinomialCoefficient(N,I)
      S3     = EXP(AE*N-ACD*M2-AAB*(NAB-K))/BinomialCoefficient(N,M2)
      slaterCondon = C*(S1-S2+S3)
      RETURN
end function GetSlaterCondonParameter

end module SlaterOverlap
