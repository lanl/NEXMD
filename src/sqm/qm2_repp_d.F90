! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
!
!  this subroutine is modified from Thiel's program
!
! (Taisung Lee, 2011)  

      SUBROUTINE qm2_repp_d(qmtypei,qmtypej,R,RI,CORE,W,LIMIJ,LImkl,MODESP) 
!     *                                                                 
!     TWO-CENTER TWO-ELECTRON REPULSION INTEGRALS AND                   
!     TWO-CENTER ONE-ELECTRON ATTRACTIONS IN LOCAL COORDINATES.         
!     GENERAL VERSION FOR MNDO/d.                                       
!     *                                                                 
!     THE TWO-ELECTRON INTEGRALS ARE RETURNED AS W(LImkl,LIMIJ)         
!     USING STANDARD PAIR INDICES IJ=(I*(I-1))/2+J FOR ADDRESSING.      
!     NOTE THAT THE ORBITAL ORDER FOR INTEGRAL EVALUATION APPLIES.      
!     AO     1     2     3     4     5     6     7     8     9          
!     L      0     1     1     1     2     2     2     2     2          
!     M      0     0     1    -1     0     1    -1     2    -2          
!     TYPE   S    PZ    PX    PY   DZ2   DXZ   DYZ DX2Y2   DXY          
!     *                                                                 
!     INPUT DATA VIA ARGUMENT LIST.                                     
!     NI       ATOMIC NUMBER OF ATOM A (ORBITALS I AND J, PAIRS IJ).    
!     NJ       ATOMIC NUMBER OF ATOM B (ORBITALS K AND L, PAIRS KL). 
!
!
!     R        INTERATOMIC DISTANCE (IN ATOMIC UNITS).                  
!     RI       FOR MODESP=1: LOCAL TWO-ELECTRON INTEGRALS FOR SP BASIS. 
!     CORE     FOR MODESP=1: LOCAL ONE-ELECTRON INTEGRALS FOR SP BASIS. 
!     LIMIJ    COLUMN DIMENSION OF W. 1 FOR S, 10 FOR SP, 45 FOR SPD.   
!     LImkl    ROW    DIMENSION OF W. 1 FOR S, 10 FOR SP, 45 FOR SPD.   
!     MODESP   OPTION FOR INTEGRAL EVALUATION.                          
!              =0   EVALUATE ALL INTEGRALS (UP TO SPD BASIS).           
!              =1   EVALUATE ONLY INTEGRALS INVOLVING D ORBITALS.       
!     *                                                                 
!     SPECIAL CONVENTION: NI=0 OR NJ=0 DENOTES AN EXTERNAL POINT        
!     CHARGE WITHOUT BASIS ORBITALS. THE CHARGE IS 1 ATOMIC UNIT.       
!     THE VALUES OF DD(I,0) AND PO(I,0) ARE DEFINED TO BE ZERO.         
!     *                                                                 
!     OUTPUT DATA VIA ARGUMENT LIST.                                    
!     CORE     COMPLETE SET OF LOCAL ONE-ELECTRON INTEGRALS.            
!     W        COMPLETE SET OF LOCAL TWO-ELECTRON INTEGRALS.            
!     *                                                                 
!     CHARGE SEPARATIONS DD(6,LMZ) AND ADDITIVE TERMS PO(9,LMZ)         
!     ARE TAKEN FROM COMMON BLOCK MULTIP.                               
!     INDEX OF DD AND PO : SS 1, SP 2, PP 3, SD 4, PD 5, DD 6.          
!     MULTIPOLE          :  L=0,  L=1,  L=2,  L=2,  L=1,  L=2.          
!     SPECIAL INDEX OF PO: PP 7, DD 8.                                  
!     MULTIPOLE          :  L=0,  L=0.                                  
!     FOR ATOMIC CORE    : ADDITIVE TERM PO(9,NI)                       
!                                                                       
!     NOTE THAT DD AND PO ARE DEFINED IN ACCORDANCE WITH THE TCA92      
!     PAPER, WITH THE FOLLOWING EXCEPTION: DD(4,N) AND DD(6,N)          
!     CONTAIN AN EXTRA FACTOR OF SQRT2 TO SIMPLIFY THE CALCULATIONS     
!     INVOLVING THE SQUARE QUADRUPOLES (L=2). SUCH A FACTOR WOULD       
!     ALSO BE HELPFUL IN DD(3,N), BUT WOULD CONFLICT WITH THE CODE      
!     FOR THE TCA77 MULTIPOLES IN SUBROUTINE REPP. THEREFORE, THE       
!     EXTRA FACTOR OF SQRT2 IS INCLUDED EXPLICITLY FOR DD(3,N) IN       
!     THIS ROUTINE (SEE BELOW).                                         
!     *                                                                 
!     THE COEFFICIENTS RELATING ANALYTICAL AND POINT-CHARGE MULTIPOLE   
!     MOMENTS, SEE EQUATION (17) AND TABLE 2 OF TCA PAPER, ARE GIVEN    
!     BELOW AS PARAMETER STATEMENTS (VARIABLES CLM..) WHERE THE DIGITS  
!     IN CLM.. REFER TO THE STANDARD PAIR INDEX (SEE ABOVE).            
!     THE CURRENT DEFINITIONS INCLUDE THE SIGNS OF THE COEFFICIENTS     
!     THAT ARE MISSING IN TABLE 2 OF THE TCA PAPER.                     
!     ONLY THOSE NONZERO COEFFICIENTS ARE DEFINED WHICH ARE NEEDED      
!     AND WHOSE ABSOLUTE VALUE IS NOT EQUAL TO 1.         
!     *                                                                 
      use constants  , only : zero, fourth, half, one, two, sqrt2, AU_TO_EV
      use qmmm_module, only : qm2_params, qmmm_struct
      use QM2_parameters
      
      implicit none
  
      integer, intent(in)::qmtypei, qmtypej, LIMIJ,LImkl,MODESP
      _REAL_, intent(in)::R      
      _REAL_, intent(out)::RI(22), CORE(10,2)
      _REAL_, intent(out)::W(LImkl,LIMIJ)
 
      _REAL_, parameter :: CLM3  = 0.13333333333333D+01 
      _REAL_, parameter :: CLM6  =-0.66666666666667D+00 
      _REAL_, parameter :: CLM10 =-0.66666666666667D+00 
      _REAL_, parameter :: CLM11 = 0.11547005383793D+01 
      _REAL_, parameter :: CLM12 = 0.11547005383793D+01 
      _REAL_, parameter :: CLM13 =-0.57735026918963D+00 
      _REAL_, parameter :: CLM15 = 0.13333333333333D+01 
      _REAL_, parameter :: CLM20 = 0.57735026918963D+00 
      _REAL_, parameter :: CLM21 = 0.66666666666667D+00 
      _REAL_, parameter :: CLM28 = 0.66666666666667D+00 
      _REAL_, parameter :: CLM33 =-0.11547005383793D+01 
      _REAL_, parameter :: CLM36 =-0.13333333333333D+01 
      _REAL_, parameter :: CLM45 =-0.13333333333333D+01 
      _REAL_, parameter :: SMALL = 1.0D-06 
      integer,parameter :: LCW(10)=(/  1,2,3,6,11,12,15,18,21,36 /)
      
      !      COMMON                                                            &
!     &/CONSTF/ A0,AFACT,EV,EVCAL,PI,W1,W2,BIGEXP                        &
!     &/CONSTN/ ZERO,ONE,TWO,THREE,FOUR,half,fourth                         &
!     &/MULTIP/ DD(6,0:LMZ),PO(9,0:LMZ)                                  &
!     &/PARDER/ TORE(LMZ),EHEAT(LMZ),EISOL(LMZ)                          


 !local    
      _REAL_:: LX(6,6) 
      _REAL_:: X(405) 
      _REAL_:: XM(14) 
      _REAL_::dd(6,2), po(9,2), tore(2)
      _REAL_::R2, DA, QA, QAS, DA2, RMDA2, RPDA2, ADD
      _REAL_::QB, DB, QBS, AB, ADDAS, ADDBS, AA, ADDSS
      _REAL_::ZZZZ, XYXY, SIG, TORENC, FACTOR
      integer:: N, NI, NJ, I, LIJ, IJ, KL, LIJMAX, LKLMAX
      integer::NC, NE, LE, ICMAX, IFIRST, ILAST 
    
      NI=1
      NJ=2     
      if (qmtypei .eq. 0) then
        dd(1:6, 1)=0.D0
        po(1:9, 1)=0.D0
        tore(1)=one
      else 
        dd(1:6, 1)= qm2_params%dd(1:6,qmtypei)
        po(1:9, 1)= qm2_params%po(1:9,qmtypei)
        tore(1)=core_chg(qmmm_struct%qm_type_id(qmtypei))
      end if
      
      if (qmtypej .eq. 0) then
        dd(1:6, 2)=0.D0
        po(1:9, 2)=0.D0
        tore(2)=one
      else      
        dd(1:6, 2)= qm2_params%dd(1:6,qmtypej)
        po(1:9, 2)= qm2_params%po(1:9,qmtypej)
        tore(2)=core_chg(qmmm_struct%qm_type_id(qmtypej))
      end if
      
      
      

!     *                                                                 
! *** INITIALIZATION.                                                   
!     *  

        RI(2)=-RI(2)
        RI(5)=-RI(5)
        RI(8:10)=-RI(8:10)
        RI(13:15)=-RI(13:15)  
                                                                               
      DO 10 IJ=1,LIMIJ 
      DO 10 KL=1,LImkl 
   10 W(KL,IJ) = ZERO 
      IF(LIMIJ.EQ.1) THEN 
         LIJMAX = 1 
      ELSE IF(LIMIJ.EQ.10) THEN 
         LIJMAX = 3 
      ELSE 
         LIJMAX = 6 
      ENDIF 
      IF(LImkl.EQ.1) THEN 
         LKLMAX = 1 
      ELSE IF(LImkl.EQ.10) THEN 
         LKLMAX = 3 
      ELSE 
         LKLMAX = 6 
      ENDIF 
!     *                                                                 
! *** DISTANCES BETWEEN POINT CHARGES.                                  
!     *                                                                 
      R2     = R*R 
      I      = 0 
      DO 180 LIJ=1,LIJMAX 
      DA     = DD(LIJ,NI) 
      QA     = PO(LIJ,NI) 
      IF(LIJ.EQ.3) THEN 
         QAS = PO(7,NI) 
         DA  = DA*SQRT2 
      ELSE IF(LIJ.EQ.6) THEN 
         QAS = PO(8,NI) 
      ELSE 
         QAS = PO(1,NI) 
      ENDIF 
      DA2    = DA*DA 
      RMDA2  = (R-DA)**2 
      RPDA2  = (R+DA)**2 
      IF(MODESP.EQ.0) THEN 
         GO TO (110,120,130,140,120,130), LIJ 
      ELSE 
         GO TO (115,125,135,140,120,130), LIJ 
      ENDIF 
! *** (SS,KL), LIJ=1.                                                   
  110 CONTINUE 
!     KL=SS.                                                            
      LX(1,LIJ) = I 
      QB     = PO(1,NJ) 
      ADD    = (QA+QB)**2 
!     L1=0, L2=0, M=0, LABEL=1.                                         
      X(I+1) = R2+ADD 
      I      = I+1 
      IF(LKLMAX.GT.1) THEN 
!        KL=PS.                                                         
         LX(2,LIJ) = I 
         DB     = DD(2,NJ) 
         QB     = PO(2,NJ) 
         ADD    = (QA+QB)**2 
!        L1=0, L2=1, M=0, LABEL=3.                                      
         X(I+1) = (R+DB)**2+ADD 
         X(I+2) = (R-DB)**2+ADD 
         I      = I+2 
!        KL=PP.                                                         
         LX(3,LIJ) = I 
         DB     = DD(3,NJ)*SQRT2 
         QB     = PO(3,NJ) 
         QBS    = PO(7,NJ) 
         ADD    = (QA+QB)**2 
!        L1=0, L2=0, M=0, LABEL=1.                                      
         X(I+1) = R2+(QA+QBS)**2 
!        L1=0, L2=2, M=0, LABEL=6.                                      
         X(I+2) = (R-DB)**2+ADD 
         X(I+3) = R2+DB**2+ADD 
         X(I+4) = (R+DB)**2+ADD 
         I      = I+4 
      ENDIF 
  115 CONTINUE 
      IF(LKLMAX.GT.3) THEN 
!        KL=DS.                                                         
         LX(4,LIJ) = I 
         DB     = DD(4,NJ) 
         QB     = PO(4,NJ) 
         ADD    = (QA+QB)**2 
!        L1=0, L2=2, M=0, LABEL=6.                                      
         X(I+1) = (R-DB)**2+ADD 
         X(I+2) = R2+DB**2+ADD 
         X(I+3) = (R+DB)**2+ADD 
         I      = I+3 
!        KL=DP.                                                         
         LX(5,LIJ) = I 
         DB     = DD(5,NJ) 
         QB     = PO(5,NJ) 
         ADD    = (QA+QB)**2 
!        L1=0, L2=1, M=0, LABEL=3.                                      
         X(I+1) = (R+DB)**2+ADD 
         X(I+2) = (R-DB)**2+ADD 
         I      = I+2 
!        KL=DD.                                                         
         LX(6,LIJ) = I 
         DB     = DD(6,NJ) 
         QB     = PO(6,NJ) 
         QBS    = PO(8,NJ) 
         ADD    = (QA+QB)**2 
!        L1=0, L2=0, M=0, LABEL=1.                                      
         X(I+1) = R2+(QA+QBS)**2 
!        L1=0, L2=2, M=0, LABEL=6.                                      
         X(I+2) = (R-DB)**2+ADD 
         X(I+3) = R2+DB**2+ADD 
         X(I+4) = (R+DB)**2+ADD 
         I      = I+4 
      ENDIF 
      GO TO 180 
! *** (PS,KL), LIJ=2, OR (DP,KL), LIJ=5.                                
  120 CONTINUE 
!     KL=SS.                                                            
      LX(1,LIJ) = I 
      QB     = PO(1,NJ) 
      ADD    = (QA+QB)**2 
!     L1=1, L2=0, M=0, LABEL=2.                                         
      X(I+1) = RPDA2+ADD 
      X(I+2) = RMDA2+ADD 
      I      = I+2 
      IF(LKLMAX.GT.1) THEN 
!        KL=PS.                                                         
         LX(2,LIJ) = I 
         DB     = DD(2,NJ) 
         QB     = PO(2,NJ) 
         ADD    = (QA+QB)**2 
!        L1=1, L2=1, M=0, LABEL=4.                                      
         X(I+1) = (R+DA-DB)**2+ADD 
         X(I+2) = (R-DA+DB)**2+ADD 
         X(I+3) = (R-DA-DB)**2+ADD 
         X(I+4) = (R+DA+DB)**2+ADD 
!        L1=1, L2=1, M=1, LABEL=5.                                      
         X(I+5) = R2+(DA-DB)**2+ADD 
         X(I+6) = R2+(DA+DB)**2+ADD 
         I      = I+6 
!        KL=PP.                                                         
         LX(3,LIJ) = I 
         DB     = DD(3,NJ)*SQRT2 
         QB     = PO(3,NJ) 
         QBS    = PO(7,NJ) 
         ADD    = (QA+QB)**2 
         ADDBS  = (QA+QBS)**2 
!        L1=1, L2=0, M=0, LABEL=2.                                      
         X(I+1) = RPDA2+ADDBS 
         X(I+2) = RMDA2+ADDBS 
!        L1=1, L2=2, M=0, LABEL=8.                                      
         X(I+3) = (R-DA-DB)**2+ADD 
         X(I+4) = RMDA2+DB**2+ADD 
         X(I+5) = (R-DA+DB)**2+ADD 
         X(I+6) = (R+DA-DB)**2+ADD 
         X(I+7) = RPDA2+DB**2+ADD 
         X(I+8) = (R+DA+DB)**2+ADD 
!        L1=1, L2=2, M=1, LABEL=11.                                     
         AB     = DB/SQRT2 
         X(I+9) = (R-AB)**2+(DA-AB)**2+ADD 
         X(I+10)= (R+AB)**2+(DA-AB)**2+ADD 
         X(I+11)= (R-AB)**2+(DA+AB)**2+ADD 
         X(I+12)= (R+AB)**2+(DA+AB)**2+ADD 
         I      = I+12 
      ENDIF 
  125 CONTINUE 
      IF(LKLMAX.GT.3) THEN 
!        KL=DS.                                                         
         LX(4,LIJ) = I 
         DB     = DD(4,NJ) 
         QB     = PO(4,NJ) 
         ADD    = (QA+QB)**2 
!        L1=1, L2=2, M=0, LABEL=8.                                      
         X(I+1) = (R-DA-DB)**2+ADD 
         X(I+2) = RMDA2+DB**2+ADD 
         X(I+3) = (R-DA+DB)**2+ADD 
         X(I+4) = (R+DA-DB)**2+ADD 
         X(I+5) = RPDA2+DB**2+ADD 
         X(I+6) = (R+DA+DB)**2+ADD 
!        L1=1, L2=2, M=1, LABEL=11.                                     
         AB     = DB/SQRT2 
         X(I+7) = (R-AB)**2+(DA-AB)**2+ADD 
         X(I+8) = (R+AB)**2+(DA-AB)**2+ADD 
         X(I+9) = (R-AB)**2+(DA+AB)**2+ADD 
         X(I+10)= (R+AB)**2+(DA+AB)**2+ADD 
         I      = I+10 
!        KL=DP.                                                         
         LX(5,LIJ) = I 
         DB     = DD(5,NJ) 
         QB     = PO(5,NJ) 
         ADD    = (QA+QB)**2 
!        L1=1, L2=1, M=0, LABEL=4.                                      
         X(I+1) = (R+DA-DB)**2+ADD 
         X(I+2) = (R-DA+DB)**2+ADD 
         X(I+3) = (R-DA-DB)**2+ADD 
         X(I+4) = (R+DA+DB)**2+ADD 
!        L1=1, L2=1, M=1, LABEL=5.                                      
         X(I+5) = R2+(DA-DB)**2+ADD 
         X(I+6) = R2+(DA+DB)**2+ADD 
         I      = I+6 
!        KL=DD.                                                         
         LX(6,LIJ) = I 
         DB     = DD(6,NJ) 
         QB     = PO(6,NJ) 
         QBS    = PO(8,NJ) 
         ADD    = (QA+QB)**2 
         ADDBS  = (QA+QBS)**2 
!        L1=1, L2=0, M=0, LABEL=2.                                      
         X(I+1) = RPDA2+ADDBS 
         X(I+2) = RMDA2+ADDBS 
!        L1=1, L2=2, M=0, LABEL=8.                                      
         X(I+3) = (R-DA-DB)**2+ADD 
         X(I+4) = RMDA2+DB**2+ADD 
         X(I+5) = (R-DA+DB)**2+ADD 
         X(I+6) = (R+DA-DB)**2+ADD 
         X(I+7) = RPDA2+DB**2+ADD 
         X(I+8) = (R+DA+DB)**2+ADD 
!        L1=1, L2=2, M=1, LABEL=11.                                     
         AB     = DB/SQRT2 
         X(I+9) = (R-AB)**2+(DA-AB)**2+ADD 
         X(I+10)= (R+AB)**2+(DA-AB)**2+ADD 
         X(I+11)= (R-AB)**2+(DA+AB)**2+ADD 
         X(I+12)= (R+AB)**2+(DA+AB)**2+ADD 
         I      = I+12 
      ENDIF 
      GO TO 180 
! *** (PP,KL), LIJ=3, OR (DD,KL), LIJ=6.                                
  130 CONTINUE 
!     KL=SS.                                                            
      LX(1,LIJ) = I 
      QB     = PO(1,NJ) 
      ADD    = (QA+QB)**2 
      ADDAS  = (QAS+QB)**2 
!     L1=0, L2=0, M=0, LABEL=1.                                         
      X(I+1) = R2+ADDAS 
!     L1=2, L2=0, M=0, LABEL=7.                                         
      X(I+2) = RMDA2+ADD 
      X(I+3) = R2+DA2+ADD 
      X(I+4) = RPDA2+ADD 
      I      = I+4 
      IF(LKLMAX.GT.1) THEN 
!        KL=PS.                                                         
         LX(2,LIJ) = I 
         DB     = DD(2,NJ) 
         QB     = PO(2,NJ) 
         ADD    = (QA+QB)**2 
         ADDAS  = (QAS+QB)**2 
!        L1=0, L2=1, M=0, LABEL=3.                                      
         X(I+1) = (R+DB)**2+ADDAS 
         X(I+2) = (R-DB)**2+ADDAS 
!        L1=2, L2=1, M=0, LABEL=9.                                      
         X(I+3) = (R-DA-DB)**2+ADD 
         X(I+4) = (R-DB)**2+DA2+ADD 
         X(I+5) = (R+DA-DB)**2+ADD 
         X(I+6) = (R-DA+DB)**2+ADD 
         X(I+7) = (R+DB)**2+DA2+ADD 
         X(I+8) = (R+DA+DB)**2+ADD 
!        L1=2, L2=1, M=1, LABEL=12.                                     
         AA     = DA/SQRT2 
         X(I+9) = (R+AA)**2+(AA-DB)**2+ADD 
         X(I+10)= (R-AA)**2+(AA-DB)**2+ADD 
         X(I+11)= (R+AA)**2+(AA+DB)**2+ADD 
         X(I+12)= (R-AA)**2+(AA+DB)**2+ADD 
         I      = I+12 
!        KL=PP.                                                         
         LX(3,LIJ) = I 
         DB     = DD(3,NJ)*SQRT2 
         QB     = PO(3,NJ) 
         QBS    = PO(7,NJ) 
         ADD    = (QA+QB)**2 
         ADDAS  = (QAS+QB)**2 
         ADDBS  = (QA+QBS)**2 
         ADDSS  = (QAS+QBS)**2 
!        L1=0, L2=0, M=0, LABEL=1.                                      
         X(I+1) = R2+ADDSS 
!        L1=0, L2=2, M=0, LABEL=6.                                      
         X(I+2) = (R-DB)**2+ADDAS 
         X(I+3) = R2+DB**2+ADDAS 
         X(I+4) = (R+DB)**2+ADDAS 
!        L1=2, L2=0, M=0, LABEL=7.                                      
         X(I+5) = RMDA2+ADDBS 
         X(I+6) = R2+DA2+ADDBS 
         X(I+7) = RPDA2+ADDBS 
!        L1=2, L2=2, M=0, LABEL=10.                                     
         X(I+8) = (R-DA-DB)**2+ADD 
         X(I+9) = (R-DA+DB)**2+ADD 
         X(I+10)= (R+DA-DB)**2+ADD 
         X(I+11)= (R+DA+DB)**2+ADD 
         X(I+12)= RMDA2+DB**2+ADD 
         X(I+13)= (R-DB)**2+DA2+ADD 
         X(I+14)= RPDA2+DB**2+ADD 
         X(I+15)= (R+DB)**2+DA2+ADD 
         X(I+16)= R2+(DA-DB)**2+ADD 
         X(I+17)= R2+(DA+DB)**2+ADD 
         X(I+18)= R2+DA2+DB**2+ADD 
!        L1=2, L2=2, M=1, LABEL=13.                                     
         AA     = DA/SQRT2 
         AB     = DB/SQRT2 
         X(I+19)= (R+AA-AB)**2+(AA-AB)**2+ADD 
         X(I+20)= (R+AA+AB)**2+(AA-AB)**2+ADD 
         X(I+21)= (R-AA-AB)**2+(AA-AB)**2+ADD 
         X(I+22)= (R-AA+AB)**2+(AA-AB)**2+ADD 
         X(I+23)= (R+AA-AB)**2+(AA+AB)**2+ADD 
         X(I+24)= (R+AA+AB)**2+(AA+AB)**2+ADD 
         X(I+25)= (R-AA-AB)**2+(AA+AB)**2+ADD 
         X(I+26)= (R-AA+AB)**2+(AA+AB)**2+ADD 
         I      = I+26 
      ENDIF 
  135 CONTINUE 
      IF(LKLMAX.GT.3) THEN 
!        KL=DS.                                                         
         LX(4,LIJ) = I 
         DB     = DD(4,NJ) 
         QB     = PO(4,NJ) 
         ADD    = (QA+QB)**2 
         ADDAS  = (QAS+QB)**2 
!        L1=0, L2=2, M=0, LABEL=6.                                      
         X(I+1) = (R-DB)**2+ADDAS 
         X(I+2) = R2+DB**2+ADDAS 
         X(I+3) = (R+DB)**2+ADDAS 
!        L1=2, L2=2, M=0, LABEL=10.                                     
         X(I+4) = (R-DA-DB)**2+ADD 
         X(I+5) = (R-DA+DB)**2+ADD 
         X(I+6) = (R+DA-DB)**2+ADD 
         X(I+7) = (R+DA+DB)**2+ADD 
         X(I+8) = RMDA2+DB**2+ADD 
         X(I+9) = (R-DB)**2+DA2+ADD 
         X(I+10)= RPDA2+DB**2+ADD 
         X(I+11)= (R+DB)**2+DA2+ADD 
         X(I+12)= R2+(DA-DB)**2+ADD 
         X(I+13)= R2+(DA+DB)**2+ADD 
         X(I+14)= R2+DA2+DB**2+ADD 
!        L1=2, L2=2, M=1, LABEL=13.                                     
         AA     = DA/SQRT2 
         AB     = DB/SQRT2 
         X(I+15)= (R+AA-AB)**2+(AA-AB)**2+ADD 
         X(I+16)= (R+AA+AB)**2+(AA-AB)**2+ADD 
         X(I+17)= (R-AA-AB)**2+(AA-AB)**2+ADD 
         X(I+18)= (R-AA+AB)**2+(AA-AB)**2+ADD 
         X(I+19)= (R+AA-AB)**2+(AA+AB)**2+ADD 
         X(I+20)= (R+AA+AB)**2+(AA+AB)**2+ADD 
         X(I+21)= (R-AA-AB)**2+(AA+AB)**2+ADD 
         X(I+22)= (R-AA+AB)**2+(AA+AB)**2+ADD 
         I      = I+22 
!        KL=DP.                                                         
         LX(5,LIJ) = I 
         DB     = DD(5,NJ) 
         QB     = PO(5,NJ) 
         ADD    = (QA+QB)**2 
         ADDAS  = (QAS+QB)**2 
!        L1=0, L2=1, M=0, LABEL=3.                                      
         X(I+1) = (R+DB)**2+ADDAS 
         X(I+2) = (R-DB)**2+ADDAS 
!        L1=2, L2=1, M=0, LABEL=9.                                      
         X(I+3) = (R-DA-DB)**2+ADD 
         X(I+4) = (R-DB)**2+DA2+ADD 
         X(I+5) = (R+DA-DB)**2+ADD 
         X(I+6) = (R-DA+DB)**2+ADD 
         X(I+7) = (R+DB)**2+DA2+ADD 
         X(I+8) = (R+DA+DB)**2+ADD 
!        L1=2, L2=1, M=1, LABEL=12.                                     
         AA     = DA/SQRT2 
         X(I+9) = (R+AA)**2+(AA-DB)**2+ADD 
         X(I+10)= (R-AA)**2+(AA-DB)**2+ADD 
         X(I+11)= (R+AA)**2+(AA+DB)**2+ADD 
         X(I+12)= (R-AA)**2+(AA+DB)**2+ADD 
         I      = I+12 
!        KL=DD.                                                         
         LX(6,LIJ) = I 
         DB     = DD(6,NJ) 
         QB     = PO(6,NJ) 
         QBS    = PO(8,NJ) 
         ADD    = (QA+QB)**2 
         ADDAS  = (QAS+QB)**2 
         ADDBS  = (QA+QBS)**2 
         ADDSS  = (QAS+QBS)**2 
!        L1=0, L2=0, M=0, LABEL=1.                                      
         X(I+1) = R2+ADDSS 
!        L1=0, L2=2, M=0, LABEL=6.                                      
         X(I+2) = (R-DB)**2+ADDAS 
         X(I+3) = R2+DB**2+ADDAS 
         X(I+4) = (R+DB)**2+ADDAS 
!        L1=2, L2=0, M=0, LABEL=7.                                      
         X(I+5) = RMDA2+ADDBS 
         X(I+6) = R2+DA2+ADDBS 
         X(I+7) = RPDA2+ADDBS 
!        L1=2, L2=2, M=0, LABEL=10.                                     
         X(I+8) = (R-DA-DB)**2+ADD 
         X(I+9) = (R-DA+DB)**2+ADD 
         X(I+10)= (R+DA-DB)**2+ADD 
         X(I+11)= (R+DA+DB)**2+ADD 
         X(I+12)= RMDA2+DB**2+ADD 
         X(I+13)= (R-DB)**2+DA2+ADD 
         X(I+14)= RPDA2+DB**2+ADD 
         X(I+15)= (R+DB)**2+DA2+ADD 
         X(I+16)= R2+(DA-DB)**2+ADD 
         X(I+17)= R2+(DA+DB)**2+ADD 
         X(I+18)= R2+DA2+DB**2+ADD 
!        L1=2, L2=2, M=1, LABEL=13.                                     
         AA     = DA/SQRT2 
         AB     = DB/SQRT2 
         X(I+19)= (R+AA-AB)**2+(AA-AB)**2+ADD 
         X(I+20)= (R+AA+AB)**2+(AA-AB)**2+ADD 
         X(I+21)= (R-AA-AB)**2+(AA-AB)**2+ADD 
         X(I+22)= (R-AA+AB)**2+(AA-AB)**2+ADD 
         X(I+23)= (R+AA-AB)**2+(AA+AB)**2+ADD 
         X(I+24)= (R+AA+AB)**2+(AA+AB)**2+ADD 
         X(I+25)= (R-AA-AB)**2+(AA+AB)**2+ADD 
         X(I+26)= (R-AA+AB)**2+(AA+AB)**2+ADD 
         I      = I+26 
      ENDIF 
      GO TO 180 
! *** (DS,KL), LIJ=4.                                                   
  140 CONTINUE 
!     KL=SS.                                                            
      LX(1,LIJ) = I 
      QB     = PO(1,NJ) 
      ADD    = (QA+QB)**2 
!     L1=2, L2=0, M=0, LABEL=7.                                         
      X(I+1) = RMDA2+ADD 
      X(I+2) = R2+DA2+ADD 
      X(I+3) = RPDA2+ADD 
      I      = I+3 
      IF(LKLMAX.GT.1) THEN 
!        KL=PS.                                                         
         LX(2,LIJ) = I 
         DB     = DD(2,NJ) 
         QB     = PO(2,NJ) 
         ADD    = (QA+QB)**2 
!        L1=2, L2=1, M=0, LABEL=9.                                      
         X(I+1) = (R-DA-DB)**2+ADD 
         X(I+2) = (R-DB)**2+DA2+ADD 
         X(I+3) = (R+DA-DB)**2+ADD 
         X(I+4) = (R-DA+DB)**2+ADD 
         X(I+5) = (R+DB)**2+DA2+ADD 
         X(I+6) = (R+DA+DB)**2+ADD 
!        L1=2, L2=1, M=1, LABEL=12.                                     
         AA     = DA/SQRT2 
         X(I+7) = (R+AA)**2+(AA-DB)**2+ADD 
         X(I+8) = (R-AA)**2+(AA-DB)**2+ADD 
         X(I+9) = (R+AA)**2+(AA+DB)**2+ADD 
         X(I+10)= (R-AA)**2+(AA+DB)**2+ADD 
         I      = I+10 
!        KL=PP.                                                         
         LX(3,LIJ) = I 
         DB     = DD(3,NJ)*SQRT2 
         QB     = PO(3,NJ) 
         QBS    = PO(7,NJ) 
         ADD    = (QA+QB)**2 
         ADDBS  = (QA+QBS)**2 
!        L1=2, L2=0, M=0, LABEL=7.                                      
         X(I+1) = RMDA2+ADDBS 
         X(I+2) = R2+DA2+ADDBS 
         X(I+3) = RPDA2+ADDBS 
!        L1=2, L2=2, M=0, LABEL=10.                                     
         X(I+4) = (R-DA-DB)**2+ADD 
         X(I+5) = (R-DA+DB)**2+ADD 
         X(I+6) = (R+DA-DB)**2+ADD 
         X(I+7) = (R+DA+DB)**2+ADD 
         X(I+8) = RMDA2+DB**2+ADD 
         X(I+9) = (R-DB)**2+DA2+ADD 
         X(I+10)= RPDA2+DB**2+ADD 
         X(I+11)= (R+DB)**2+DA2+ADD 
         X(I+12)= R2+(DA-DB)**2+ADD 
         X(I+13)= R2+(DA+DB)**2+ADD 
         X(I+14)= R2+DA2+DB**2+ADD 
!        L1=2, L2=2, M=1, LABEL=13.                                     
         AA     = DA/SQRT2 
         AB     = DB/SQRT2 
         X(I+15)= (R+AA-AB)**2+(AA-AB)**2+ADD 
         X(I+16)= (R+AA+AB)**2+(AA-AB)**2+ADD 
         X(I+17)= (R-AA-AB)**2+(AA-AB)**2+ADD 
         X(I+18)= (R-AA+AB)**2+(AA-AB)**2+ADD 
         X(I+19)= (R+AA-AB)**2+(AA+AB)**2+ADD 
         X(I+20)= (R+AA+AB)**2+(AA+AB)**2+ADD 
         X(I+21)= (R-AA-AB)**2+(AA+AB)**2+ADD 
         X(I+22)= (R-AA+AB)**2+(AA+AB)**2+ADD 
         I      = I+22 
      ENDIF 
      IF(LKLMAX.GT.3) THEN 
!        KL=DS.                                                         
         LX(4,LIJ) = I 
         DB     = DD(4,NJ) 
         QB     = PO(4,NJ) 
         ADD    = (QA+QB)**2 
!        L1=2, L2=2, M=0, LABEL=10.                                     
         X(I+1) = (R-DA-DB)**2+ADD 
         X(I+2) = (R-DA+DB)**2+ADD 
         X(I+3) = (R+DA-DB)**2+ADD 
         X(I+4) = (R+DA+DB)**2+ADD 
         X(I+5) = RMDA2+DB**2+ADD 
         X(I+6) = (R-DB)**2+DA2+ADD 
         X(I+7) = RPDA2+DB**2+ADD 
         X(I+8) = (R+DB)**2+DA2+ADD 
         X(I+9) = R2+(DA-DB)**2+ADD 
         X(I+10)= R2+(DA+DB)**2+ADD 
         X(I+11)= R2+DA2+DB**2+ADD 
!        L1=2, L2=2, M=1, LABEL=13.                                     
         AA     = DA/SQRT2 
         AB     = DB/SQRT2 
         X(I+12)= (R+AA-AB)**2+(AA-AB)**2+ADD 
         X(I+13)= (R+AA+AB)**2+(AA-AB)**2+ADD 
         X(I+14)= (R-AA-AB)**2+(AA-AB)**2+ADD 
         X(I+15)= (R-AA+AB)**2+(AA-AB)**2+ADD 
         X(I+16)= (R+AA-AB)**2+(AA+AB)**2+ADD 
         X(I+17)= (R+AA+AB)**2+(AA+AB)**2+ADD 
         X(I+18)= (R-AA-AB)**2+(AA+AB)**2+ADD 
         X(I+19)= (R-AA+AB)**2+(AA+AB)**2+ADD 
         I      = I+19 
!        KL=DP.                                                         
         LX(5,LIJ) = I 
         DB     = DD(5,NJ) 
         QB     = PO(5,NJ) 
         ADD    = (QA+QB)**2 
!        L1=2, L2=1, M=0, LABEL=9.                                      
         X(I+1) = (R-DA-DB)**2+ADD 
         X(I+2) = (R-DB)**2+DA2+ADD 
         X(I+3) = (R+DA-DB)**2+ADD 
         X(I+4) = (R-DA+DB)**2+ADD 
         X(I+5) = (R+DB)**2+DA2+ADD 
         X(I+6) = (R+DA+DB)**2+ADD 
!        L1=2, L2=1, M=1, LABEL=12.                                     
         AA     = DA/SQRT2 
         X(I+7) = (R+AA)**2+(AA-DB)**2+ADD 
         X(I+8) = (R-AA)**2+(AA-DB)**2+ADD 
         X(I+9) = (R+AA)**2+(AA+DB)**2+ADD 
         X(I+10)= (R-AA)**2+(AA+DB)**2+ADD 
         I      = I+10 
!        KL=DD.                                                         
         LX(6,LIJ) = I 
         DB     = DD(6,NJ) 
         QB     = PO(6,NJ) 
         QBS    = PO(8,NJ) 
         ADD    = (QA+QB)**2 
         ADDBS  = (QA+QBS)**2 
!        L1=2, L2=0, M=0, LABEL=7.                                      
         X(I+1) = RMDA2+ADDBS 
         X(I+2) = R2+DA2+ADDBS 
         X(I+3) = RPDA2+ADDBS 
!        L1=2, L2=2, M=0, LABEL=10.                                     
         X(I+4) = (R-DA-DB)**2+ADD 
         X(I+5) = (R-DA+DB)**2+ADD 
         X(I+6) = (R+DA-DB)**2+ADD 
         X(I+7) = (R+DA+DB)**2+ADD 
         X(I+8) = RMDA2+DB**2+ADD 
         X(I+9) = (R-DB)**2+DA2+ADD 
         X(I+10)= RPDA2+DB**2+ADD 
         X(I+11)= (R+DB)**2+DA2+ADD 
         X(I+12)= R2+(DA-DB)**2+ADD 
         X(I+13)= R2+(DA+DB)**2+ADD 
         X(I+14)= R2+DA2+DB**2+ADD 
!        L1=2, L2=2, M=1, LABEL=13.                                     
         AA     = DA/SQRT2 
         AB     = DB/SQRT2 
         X(I+15)= (R+AA-AB)**2+(AA-AB)**2+ADD 
         X(I+16)= (R+AA+AB)**2+(AA-AB)**2+ADD 
         X(I+17)= (R-AA-AB)**2+(AA-AB)**2+ADD 
         X(I+18)= (R-AA+AB)**2+(AA-AB)**2+ADD 
         X(I+19)= (R+AA-AB)**2+(AA+AB)**2+ADD 
         X(I+20)= (R+AA+AB)**2+(AA+AB)**2+ADD 
         X(I+21)= (R-AA-AB)**2+(AA+AB)**2+ADD 
         X(I+22)= (R-AA+AB)**2+(AA+AB)**2+ADD 
         I      = I+22 
      ENDIF 
  180 END DO 
      ILAST  = I 
! *** CALCULATE INVERSE DISTANCES.                                      
!     NUMERATOR ONE IMPLIES ATOMIC UNITS.                               
!     NUMERATOR EV  IMPLIES CONVERSION TO EV AT THIS POINT ALREADY.     
!     THE CONVERSION IS DONE HERE FOR COMPUTATIONAL EFFICIENCY.         
!     ALTERNATIVELY IT COULD BE INCLUDED AFTER THE DO 280 LOOP          
!     AT THE LEVEL OF THE TWO-ELECTRON INTEGRALS W(LImkl,LIMIJ).        
      DO 190 I=1,ILAST 
  190 X(I)   = AU_TO_EV/SQRT(X(I)) 
!     *                                                                 
!     EVALUATE SEMIEMPIRICAL MULTIPOLE-MULTIPOLE INTERACTIONS XM(LABEL).
!     EVALUATE THE UNIQUE INTEGRALS W(KL,IJ) FROM XM(LABEL).            
!     DEFINE THE REMAINING SYMMETRY-RELATED INTEGRALS W(KL,IJ).         
!     *                                                                 
      DO 280 LIJ=1,LIJMAX 
      IF(MODESP.EQ.0) THEN 
         GO TO (210,220,230,240,250,260), LIJ 
      ELSE 
         GO TO (215,225,235,240,250,260), LIJ 
      ENDIF 
! *** (SS,KL), LIJ=1.                                                   
  210 CONTINUE 
!     KL=SS.                                                            
      I      = LX(1,LIJ) 
      XM(1) = X(I+1) 
      W( 1, 1) = XM( 1) 
      IF(LKLMAX.GT.1) THEN 
!        KL=PS.                                                         
         I   = LX(2,LIJ) 
         XM(3) = (X(I+1)-X(I+2))*half 
         W( 2, 1) = XM( 3) 
!        KL=PP.                                                         
         I   = LX(3,LIJ) 
         XM(1) =  X(I+1) 
         XM(6) = (X(I+2)-X(I+3)*TWO+X(I+4))*fourth 
         W( 3, 1) = XM( 1)                                              &
     &             +XM( 6) * CLM3                                       
         W( 6, 1) = XM( 1)                                              &
     &             +XM( 6) * CLM6                                       
         W(10, 1) = W( 6, 1) 
      ENDIF 
  215 CONTINUE 
      IF(LKLMAX.GT.3) THEN 
!        KL=DS.                                                         
         I   = LX(4,LIJ) 
         XM(6) = (X(I+1)-X(I+2)*TWO+X(I+3))*fourth 
         W(11, 1) = XM( 6) * CLM11 
!        KL=DP.                                                         
         I   = LX(5,LIJ) 
         XM(3) = (X(I+1)-X(I+2))*half 
         W(12, 1) = XM( 3) * CLM12 
         W(18, 1) = XM( 3) 
         W(25, 1) = W(18, 1) 
!        KL=DD.                                                         
         I   = LX(6,LIJ) 
         XM(1) =  X(I+1) 
         XM(6) = (X(I+2)-X(I+3)*TWO+X(I+4))*fourth 
         W(15, 1) = XM( 1)                                              &
     &             +XM( 6) * CLM15                                      
         W(21, 1) = XM( 1)                                              &
     &             +XM( 6) * CLM21                                      
         W(36, 1) = XM( 1)                                              &
     &             +XM( 6) * CLM36                                      
         W(28, 1) = W(21, 1) 
         W(45, 1) = W(36, 1) 
      ENDIF 
      GO TO 280 
! *** (PS,KL), LIJ=2.                                                   
  220 CONTINUE 
!     KL=SS.                                                            
      I      = LX(1,LIJ) 
      XM(2) =-(X(I+1)-X(I+2))*half 
      W( 1, 2) = XM( 2) 
      IF(LKLMAX.GT.1) THEN 
!        KL=PS.                                                         
         I  = LX(2,LIJ) 
         XM(4) = (X(I+1)+X(I+2)-X(I+3)-X(I+4))*fourth 
         XM(5) = (X(I+5)-X(I+6))*half 
         W( 2, 2) = XM( 4) 
         W( 4, 4) = XM( 5) 
         W( 7, 7) = W( 4, 4) 
!        KL=PP.                                                         
         I   = LX(3,LIJ) 
         XM(2) =-(X(I+1)-X(I+2))*half 
         XM(8) = (X(I+3)-X(I+4)*TWO+X(I+5)                              &
     &           -X(I+6)+X(I+7)*TWO-X(I+8))/8.0D0                       
         XM(11)=-(X(I+9)-X(I+10)-X(I+11)+X(I+12))*fourth 
         W( 3, 2) = XM( 2)                                              &
     &             +XM( 8) * CLM3                                       
         W( 6, 2) = XM( 2)                                              &
     &             +XM( 8) * CLM6                                       
         W( 5, 4) = XM(11) 
         W(10, 2) = W( 6, 2) 
         W( 8, 7) = W( 5, 4) 
      ENDIF 
  225 CONTINUE 
      IF(LKLMAX.GT.3) THEN 
!        KL=DS.                                                         
         I   = LX(4,LIJ) 
         XM(8) = (X(I+1)-X(I+2)*TWO+X(I+3)                              &
     &           -X(I+4)+X(I+5)*TWO-X(I+6))/8.0D0                       
         XM(11)=-(X(I+7)-X(I+8)-X(I+9)+X(I+10))*fourth 
         W(11, 2) = XM( 8) * CLM11 
         W(16, 4) = XM(11) 
         W(22, 7) = W(16, 4) 
!        KL=DP.                                                         
         I   = LX(5,LIJ) 
         XM(4) = (X(I+1)+X(I+2)-X(I+3)-X(I+4))*fourth 
         XM(5) = (X(I+5)-X(I+6))*half 
         W(12, 2) = XM( 4) * CLM12 
         W(18, 2) = XM( 4) 
         W(13, 4) = XM( 5) * CLM13 
         W(17, 4) = XM( 5) 
         W(31, 4) = XM( 5) 
         W(25, 2) = W(18, 2) 
         W(40, 4) = W(31, 4) 
         W(14, 7) = W(13, 4) 
         W(23, 7) = W(17, 4) 
         W(32, 7) =-W(31, 4) 
         W(39, 7) = W(31, 4) 
!        KL=DD.                                                         
         I   = LX(6,LIJ) 
         XM(2) =-(X(I+1)-X(I+2))*half 
         XM(8) = (X(I+3)-X(I+4)*TWO+X(I+5)                              &
     &           -X(I+6)+X(I+7)*TWO-X(I+8))/8.0D0                       
         XM(11)=-(X(I+9)-X(I+10)-X(I+11)+X(I+12))*fourth 
         W(15, 2) = XM( 2)                                              &
     &             +XM( 8) * CLM15                                      
         W(21, 2) = XM( 2)                                              &
     &             +XM( 8) * CLM21                                      
         W(36, 2) = XM( 2)                                              &
     &             +XM( 8) * CLM36                                      
         W(20, 4) = XM(11) * CLM20 
         W(34, 4) = XM(11) 
         W(28, 2) = W(21, 2) 
         W(45, 2) = W(36, 2) 
         W(43, 4) = W(34, 4) 
         W(26, 7) = W(20, 4) 
         W(35, 7) =-W(34, 4) 
         W(42, 7) = W(34, 4) 
      ENDIF 
      GO TO 280 
! *** (PP,KL), LIJ=3.                                                   
  230 CONTINUE 
!     KL=SS.                                                            
      I      = LX(1,LIJ) 
      XM(1) =  X(I+1) 
      XM(7) = (X(I+2)-X(I+3)*TWO+X(I+4))*fourth 
      W( 1, 3) = XM( 1)                                                 &
     &          +XM( 7) * CLM3                                          
      W( 1, 6) = XM( 1)                                                 &
     &          +XM( 7) * CLM6                                          
      W( 1,10) = W( 1, 6) 
      IF(LKLMAX.GT.1) THEN 
!        KL=PS.                                                         
         I   = LX(2,LIJ) 
         XM(3) = (X(I+1)-X(I+2))*half 
         XM(9) =-(X(I+3)-X(I+4)*TWO+X(I+5)                              &
     &           -X(I+6)+X(I+7)*TWO-X(I+8))/8.0D0                       
         XM(12)=-(X(I+9)-X(I+10)-X(I+11)+X(I+12))*fourth 
         W( 2, 3) = XM( 3)                                              &
     &             +XM( 9) * CLM3                                       
         W( 4, 5) = XM(12) 
         W( 2, 6) = XM( 3)                                              &
     &             +XM( 9) * CLM6                                       
         W( 7, 8) = W( 4, 5) 
         W( 2,10) = W( 2, 6) 
!        KL=PP.                                                         
         I   = LX(3,LIJ) 
         XM(1) =  X(I+1) 
         XM(6) = (X(I+2)-X(I+3)*TWO+X(I+4))*fourth 
         XM(7) = (X(I+5)-X(I+6)*TWO+X(I+7))*fourth 
         ZZZZ  =  X(I+8)+X(I+9)+X(I+10)+X(I+11) - TWO*                  &
     &           (X(I+12)+X(I+13)+X(I+14)+X(I+15)-X(I+16)-X(I+17))      
         XYXY  =  X(I+16)+X(I+17)-X(I+18)*TWO 
         XM(10)= (ZZZZ-XYXY)/16.0D0 
         XM(13)= (X(I+19)-X(I+20)-X(I+21)+X(I+22)                       &
     &           -X(I+23)+X(I+24)+X(I+25)-X(I+26))/8.0D0                
         XM(14)= XYXY*fourth 
         W( 3, 3) = XM( 1)                                              &
     &             +XM( 6) * CLM3                                       &
     &             +XM( 7) * CLM3                                       &
     &             +XM(10) * CLM3*CLM3                                  
         W( 6, 3) = XM( 1)                                              &
     &             +XM( 6) * CLM6                                       &
     &             +XM( 7) * CLM3                                       &
     &             +XM(10) * CLM3*CLM6                                  
         W( 5, 5) = XM(13) 
         W( 3, 6) = XM( 1)                                              &
     &             +XM( 6) * CLM3                                       &
     &             +XM( 7) * CLM6                                       &
     &             +XM(10) * CLM6*CLM3                                  
         W( 6, 6) = XM( 1)                                              &
     &             +XM( 6) * CLM6                                       &
     &             +XM( 7) * CLM6                                       &
     &             +XM(10) * CLM6*CLM6                                  &
     &             +XM(14)                                              
         W(10, 6) = XM( 1)                                              &
     &             +XM( 6) * CLM10                                      &
     &             +XM( 7) * CLM6                                       &
     &             +XM(10) * CLM6*CLM10                                 &
     &             -XM(14)                                              
         W( 9, 9) = XM(14) 
         W(10, 3) = W( 6, 3) 
         W( 8, 8) = W( 5, 5) 
         W( 3,10) = W( 3, 6) 
         W( 6,10) = W(10, 6) 
         W(10,10) = W( 6, 6) 
      ENDIF 
  235 CONTINUE 
      IF(LKLMAX.GT.3) THEN 
!        KL=DS.                                                         
         I   = LX(4,LIJ) 
         XM(6) = (X(I+1)-X(I+2)*TWO+X(I+3))*fourth 
         ZZZZ  =  X(I+4)+X(I+5)+X(I+6)+X(I+7) - TWO*                    &
     &           (X(I+8)+X(I+9)+X(I+10)+X(I+11)-X(I+12)-X(I+13))        
         XYXY  =  X(I+12)+X(I+13)-X(I+14)*TWO 
         XM(10)= (ZZZZ-XYXY)/16.0D0 
         XM(13)= (X(I+15)-X(I+16)-X(I+17)+X(I+18)                       &
     &           -X(I+19)+X(I+20)+X(I+21)-X(I+22))/8.0D0                
         XM(14)= XYXY*fourth 
         W(11, 3) = XM( 6) * CLM11                                      &
     &             +XM(10) * CLM3*CLM11                                 
         W(16, 5) = XM(13) 
         W(11, 6) = XM( 6) * CLM11                                      &
     &             +XM(10) * CLM6*CLM11                                 
         W(29, 6) = XM(14) 
         W(22, 8) = W(16, 5) 
         W(37, 9) = W(29, 6) 
         W(11,10) = W(11, 6) 
         W(29,10) =-W(29, 6) 
!        KL=DP.                                                         
         I   = LX(5,LIJ) 
         XM(3) = (X(I+1)-X(I+2))*half 
         XM(9) =-(X(I+3)-X(I+4)*TWO+X(I+5)                              &
     &           -X(I+6)+X(I+7)*TWO-X(I+8))/8.0D0                       
         XM(12)=-(X(I+9)-X(I+10)-X(I+11)+X(I+12))*fourth 
         W(12, 3) = XM( 3) * CLM12                                      &
     &             +XM( 9) * CLM3*CLM12                                 
         W(18, 3) = XM( 3)                                              &
     &             +XM( 9) * CLM3                                       
         W(13, 5) = XM(12) * CLM13 
         W(17, 5) = XM(12) 
         W(31, 5) = XM(12) 
         W(12, 6) = XM( 3) * CLM12                                      &
     &             +XM( 9) * CLM6*CLM12                                 
         W(18, 6) = XM( 3)                                              &
     &             +XM( 9) * CLM6                                       
         W(25, 6) = XM( 3)                                              &
     &             +XM( 9) * CLM6                                       
         W(25, 3) = W(18, 3) 
         W(40, 5) = W(31, 5) 
         W(14, 8) = W(13, 5) 
         W(23, 8) = W(17, 5) 
         W(32, 8) =-W(31, 5) 
         W(39, 8) = W(31, 5) 
         W(19, 9) = W(30, 6) 
         W(24, 9) = W(30, 6) 
         W(38, 9) = W(30, 6) 
         W(12,10) = W(12, 6) 
         W(18,10) = W(25, 6) 
         W(25,10) = W(18, 6) 
         W(30,10) =-W(30, 6) 
!        KL=DD.                                                         
         I   = LX(6,LIJ) 
         XM(1) =  X(I+1) 
         XM(6) = (X(I+2)-X(I+3)*TWO+X(I+4))*fourth 
         XM(7) = (X(I+5)-X(I+6)*TWO+X(I+7))*fourth 
         ZZZZ  =  X(I+8)+X(I+9)+X(I+10)+X(I+11) - TWO*                  &
     &           (X(I+12)+X(I+13)+X(I+14)+X(I+15)-X(I+16)-X(I+17))      
         XYXY  =  X(I+16)+X(I+17)-X(I+18)*TWO 
         XM(10)= (ZZZZ-XYXY)/16.0D0 
         XM(13)= (X(I+19)-X(I+20)-X(I+21)+X(I+22)                       &
     &           -X(I+23)+X(I+24)+X(I+25)-X(I+26))/8.0D0                
         XM(14)= XYXY*fourth 
         W(15, 3) = XM( 1)                                              &
     &             +XM( 6) * CLM15                                      &
     &             +XM( 7) * CLM3                                       &
     &             +XM(10) * CLM3*CLM15                                 
         W(21, 3) = XM( 1)                                              &
     &             +XM( 6) * CLM21                                      &
     &             +XM( 7) * CLM3                                       &
     &             +XM(10) * CLM3*CLM21                                 
         W(36, 3) = XM( 1)                                              &
     &             +XM( 6) * CLM36                                      &
     &             +XM( 7) * CLM3                                       &
     &             +XM(10) * CLM3*CLM36                                 
         W(20, 5) = XM(13) * CLM20 
         W(34, 5) = XM(13) 
         W(15, 6) = XM( 1)                                              &
     &             +XM( 6) * CLM15                                      &
     &             +XM( 7) * CLM6                                       &
     &             +XM(10) * CLM6*CLM15                                 
         W(21, 6) = XM( 1)                                              &
     &             +XM( 6) * CLM21                                      &
     &             +XM( 7) * CLM6                                       &
     &             +XM(10) * CLM6*CLM21                                 &
     &             +XM(14)                                              
         W(28, 6) = XM( 1)                                              &
     &             +XM( 6) * CLM28                                      &
     &             +XM( 7) * CLM6                                       &
     &             +XM(10) * CLM6*CLM28                                 &
     &             -XM(14)                                              
         W(33, 6) = XM(14) * CLM33 
         W(36, 6) = XM( 1)                                              &
     &             +XM( 6) * CLM36                                      &
     &             +XM( 7) * CLM6                                       &
     &             +XM(10) * CLM6*CLM36                                 
         W(27, 9) = XM(14) 
         W(28, 3) = W(21, 3) 
         W(45, 3) = W(36, 3) 
         W(43, 5) = W(34, 5) 
         W(45, 6) = W(36, 6) 
         W(26, 8) = W(20, 5) 
         W(35, 8) =-W(34, 5) 
         W(42, 8) = W(34, 5) 
         W(41, 9) = W(33, 6) 
         W(15,10) = W(15, 6) 
         W(21,10) = W(28, 6) 
         W(28,10) = W(21, 6) 
         W(33,10) =-W(33, 6) 
         W(36,10) = W(36, 6) 
         W(45,10) = W(36, 6) 
      ENDIF 
      GO TO 280 
! *** (DS,KL), LIJ=4.                                                   
  240 CONTINUE 
!     KL=SS.                                                            
      I      = LX(1,LIJ) 
      XM(7) = (X(I+1)-X(I+2)*TWO+X(I+3))*fourth 
      W( 1,11) = XM( 7) * CLM11 
      IF(LKLMAX.GT.1) THEN 
!        KL=PS.                                                         
         I   = LX(2,LIJ) 
         XM(9) =-(X(I+1)-X(I+2)*TWO+X(I+3)                              &
     &                 -X(I+4)+X(I+5)*TWO-X(I+6))/8.0D0                 
         XM(12)=-(X(I+7)-X(I+8)-X(I+9)+X(I+10))*fourth 
         W( 2,11) = XM( 9) * CLM11 
         W( 4,16) = XM(12) 
         W( 7,22) = W( 4,16) 
!        KL=PP.                                                         
         I   = LX(3,LIJ) 
         XM(7) = (X(I+1)-X(I+2)*TWO+X(I+3))*fourth 
         ZZZZ  =  X(I+4)+X(I+5)+X(I+6)+X(I+7) - TWO*                    &
     &           (X(I+8)+X(I+9)+X(I+10)+X(I+11)-X(I+12)-X(I+13))        
         XYXY  =  X(I+12)+X(I+13)-X(I+14)*TWO 
         XM(10)= (ZZZZ-XYXY)/16.0D0 
         XM(13)= (X(I+15)-X(I+16)-X(I+17)+X(I+18)                       &
     &           -X(I+19)+X(I+20)+X(I+21)-X(I+22))/8.0D0                
         XM(14)= XYXY*fourth 
         W( 3,11) = XM( 7) * CLM11                                      &
     &             +XM(10) * CLM11*CLM3                                 
         W( 6,11) = XM( 7) * CLM11                                      &
     &             +XM(10) * CLM11*CLM6                                 
         W( 5,16) = XM(13) 
         W( 6,29) = XM(14) 
         W(10,11) = W( 6,11) 
         W( 8,22) = W( 5,16) 
         W(10,29) =-W( 6,29) 
         W( 9,37) = W( 6,29) 
      ENDIF 
      IF(LKLMAX.GT.3) THEN 
!        KL=DS.                                                         
         I   = LX(4,LIJ) 
         ZZZZ  =  X(I+1)+X(I+2)+X(I+3)+X(I+4) - TWO*                    &
     &           (X(I+5)+X(I+6)+X(I+7)+X(I+8)-X(I+9)-X(I+10))           
         XYXY  =  X(I+9)+X(I+10)-X(I+11)*TWO 
         XM(10)= (ZZZZ-XYXY)/16.0D0 
         XM(13)= (X(I+12)-X(I+13)-X(I+14)+X(I+15)                       &
     &           -X(I+16)+X(I+17)+X(I+18)-X(I+19))/8.0D0                
         XM(14)= XYXY*fourth 
         W(11,11) = XM(10) * CLM11*CLM11 
         W(16,16) = XM(13) 
         W(29,29) = XM(14) 
         W(22,22) = W(16,16) 
         W(37,37) = W(29,29) 
!        KL=DP.                                                         
         I   = LX(5,LIJ) 
         XM(9) =-(X(I+1)-X(I+2)*TWO+X(I+3)                              &
     &           -X(I+4)+X(I+5)*TWO-X(I+6))/8.0D0                       
         XM(12)=-(X(I+7)-X(I+8)-X(I+9)+X(I+10))*fourth 
         W(12,11) = XM( 9) * CLM11*CLM12 
         W(18,11) = XM( 9) * CLM11 
         W(13,16) = XM(12) * CLM13 
         W(17,16) = XM(12) 
         W(31,16) = XM(12) 
         W(25,11) = W(18,11) 
         W(40,16) = W(31,16) 
         W(14,22) = W(13,16) 
         W(23,22) = W(17,16) 
         W(32,22) =-W(31,16) 
         W(39,22) = W(31,16) 
         W(25,29) =-W(18,29) 
         W(30,29) = W(18,29) 
         W(19,37) = W(18,29) 
         W(24,37) = W(18,29) 
         W(38,37) = W(18,29) 
!        KL=DD.                                                         
         I   = LX(6,LIJ) 
         XM(7) = (X(I+1)-X(I+2)*TWO+X(I+3))*fourth 
         ZZZZ  =  X(I+4)+X(I+5)+X(I+6)+X(I+7) - TWO*                    &
     &           (X(I+8)+X(I+9)+X(I+10)+X(I+11)-X(I+12)-X(I+13))        
         XYXY  =  X(I+12)+X(I+13)-X(I+14)*TWO 
         XM(10)= (ZZZZ-XYXY)/16.0D0 
         XM(13)= (X(I+15)-X(I+16)-X(I+17)+X(I+18)                       &
     &           -X(I+19)+X(I+20)+X(I+21)-X(I+22))/8.0D0                
         XM(14)= XYXY*fourth 
         W(15,11) = XM( 7) * CLM11                                      &
     &             +XM(10) * CLM11*CLM15                                
         W(21,11) = XM( 7) * CLM11                                      &
     &             +XM(10) * CLM11*CLM21                                
         W(36,11) = XM( 7) * CLM11                                      &
     &             +XM(10) * CLM11*CLM36                                
         W(20,16) = XM(13) * CLM20 
         W(34,16) = XM(13) 
         W(21,29) = XM(14) 
         W(33,29) = XM(14) * CLM33 
         W(28,11) = W(21,11) 
         W(45,11) = W(36,11) 
         W(43,16) = W(34,16) 
         W(26,22) = W(20,16) 
         W(35,22) =-W(34,16) 
         W(42,22) = W(34,16) 
         W(28,29) =-W(21,29) 
         W(27,37) = W(21,29) 
         W(41,37) = W(33,29) 
      ENDIF 
      GO TO 280 
! *** (DP,KL), LIJ=5.                                                   
  250 CONTINUE 
!     KL=SS.                                                            
      I      = LX(1,LIJ) 
      XM(2) =-(X(I+1)-X(I+2))*half 
      W( 1,12) = XM( 2) * CLM12 
      W( 1,18) = XM( 2) 
      W( 1,25) = W( 1,18) 
      IF(LKLMAX.GT.1) THEN 
!        KL=PS.                                                         
         I   = LX(2,LIJ) 
         XM(4) = (X(I+1)+X(I+2)-X(I+3)-X(I+4))*fourth 
         XM(5) = (X(I+5)-X(I+6))*half 
         W( 2,12) = XM( 4) * CLM12 
         W( 4,13) = XM( 5) * CLM13 
         W( 4,17) = XM( 5) 
         W( 2,18) = XM( 4) 
         W( 4,31) = XM( 5) 
         W( 7,14) = W( 4,13) 
         W( 7,23) = W( 4,17) 
         W( 2,25) = W( 2,18) 
         W( 7,32) =-W( 4,31) 
         W( 7,39) = W( 4,31) 
         W( 4,40) = W( 4,31) 
!        KL=PP.                                                         
         I   = LX(3,LIJ) 
         XM(2) =-(X(I+1)-X(I+2))*half 
         XM(8) = (X(I+3)-X(I+4)*TWO+X(I+5)                              &
     &           -X(I+6)+X(I+7)*TWO-X(I+8))/8.0D0                       
         XM(11)=-(X(I+9)-X(I+10)-X(I+11)+X(I+12))*fourth 
         W( 3,12) = XM( 2) * CLM12                                      &
     &             +XM( 8) * CLM12*CLM3                                 
         W( 6,12) = XM( 2) * CLM12                                      &
     &             +XM( 8) * CLM12*CLM6                                 
         W( 5,13) = XM(11) * CLM13 
         W( 5,17) = XM(11) 
         W( 3,18) = XM( 2)                                              &
     &             +XM( 8) * CLM3                                       
         W( 6,18) = XM( 2)                                              &
     &             +XM( 8) * CLM6                                       
         W(10,18) = XM( 2)                                              &
     &             +XM( 8) * CLM10                                      
         W( 5,31) = XM(11) 
         W(10,12) = W( 6,12) 
         W( 8,14) = W( 5,13) 
         W( 8,23) = W( 5,17) 
         W( 9,24) = W( 9,19) 
         W( 3,25) = W( 3,18) 
         W( 6,25) = W(10,18) 
         W(10,25) = W( 6,18) 
         W( 6,30) = W( 9,19) 
         W(10,30) =-W( 9,19) 
         W( 8,32) =-W( 5,31) 
         W( 9,38) = W( 9,19) 
         W( 8,39) = W( 5,31) 
         W( 5,40) = W( 5,31) 
      ENDIF 
      IF(LKLMAX.GT.3) THEN 
!        KL=DS.                                                         
         I   = LX(4,LIJ) 
         XM(8) = (X(I+1)-X(I+2)*TWO+X(I+3)                              &
     &           -X(I+4)+X(I+5)*TWO-X(I+6))/8.0D0                       
         XM(11)=-(X(I+7)-X(I+8)-X(I+9)+X(I+10))*fourth 
         W(11,12) = XM( 8) * CLM12*CLM11 
         W(16,13) = XM(11) * CLM13 
         W(16,17) = XM(11)  
         W(11,18) = XM( 8) * CLM11 
         W(16,31) = XM(11) 
         W(22,14) = W(16,13) 
         W(37,19) = W(29,18) 
         W(22,23) = W(16,17) 
         W(37,24) = W(29,18) 
         W(11,25) = W(11,18) 
         W(29,25) =-W(29,18) 
         W(29,30) = W(29,18) 
         W(22,32) =-W(16,31) 
         W(37,38) = W(29,18) 
         W(22,39) = W(16,31) 
         W(16,40) = W(16,31) 
!        KL=DP.                                                         
         I   = LX(5,LIJ) 
         XM(4) = (X(I+1)+X(I+2)-X(I+3)-X(I+4))*fourth 
         XM(5) = (X(I+5)-X(I+6))*half 
         W(12,12) = XM( 4) * CLM12*CLM12 
         W(18,12) = XM( 4) * CLM12 
         W(13,13) = XM( 5) * CLM13*CLM13 
         W(17,13) = XM( 5) * CLM13 
         W(31,13) = XM( 5) * CLM13 
         W(13,17) = XM( 5) * CLM13 
         W(17,17) = XM( 5) 
         W(31,17) = XM( 5) 
         W(12,18) = XM( 4) * CLM12 
         W(18,18) = XM( 4) 
         W(25,18) = XM( 4) 
         W(13,31) = XM( 5) * CLM13 
         W(17,31) = XM( 5) 
         W(31,31) = XM( 5) 
         W(40,31) = XM( 5) 
         W(25,12) = W(18,12) 
         W(40,13) = W(31,13) 
         W(14,14) = W(13,13) 
         W(23,14) = W(17,13) 
         W(32,14) =-W(31,13) 
         W(39,14) = W(31,13) 
         W(40,17) = W(31,17) 
         W(19,19) = W(30,18) 
         W(24,19) = W(30,18) 
         W(38,19) = W(30,18) 
         W(14,23) = W(13,17) 
         W(23,23) = W(17,17) 
         W(32,23) =-W(31,17) 
         W(39,23) = W(31,17) 
         W(19,24) = W(30,18) 
         W(24,24) = W(30,18) 
         W(38,24) = W(30,18) 
         W(12,25) = W(12,18) 
         W(18,25) = W(25,18) 
         W(25,25) = W(18,18) 
         W(30,25) =-W(30,18) 
         W(18,30) = W(30,18) 
         W(25,30) =-W(30,18) 
         W(30,30) = W(30,18) 
         W(14,32) =-W(13,31) 
         W(23,32) =-W(17,31) 
         W(32,32) = W(31,31) 
         W(39,32) =-W(40,31) 
         W(19,38) = W(30,18) 
         W(24,38) = W(30,18) 
         W(38,38) = W(30,18) 
         W(14,39) = W(13,31) 
         W(23,39) = W(17,31) 
         W(32,39) =-W(40,31) 
         W(39,39) = W(31,31) 
         W(13,40) = W(13,31) 
         W(17,40) = W(17,31) 
         W(31,40) = W(40,31) 
         W(40,40) = W(31,31) 
!        KL=DD.                                                         
         I   = LX(6,LIJ) 
         XM(2) =-(X(I+1)-X(I+2))*half 
         XM(8) = (X(I+3)-X(I+4)*TWO+X(I+5)                              &
     &           -X(I+6)+X(I+7)*TWO-X(I+8))/8.0D0                       
         XM(11)=-(X(I+9)-X(I+10)-X(I+11)+X(I+12))*fourth 
         W(15,12) = XM( 2) * CLM12                                      &
     &             +XM( 8) * CLM12*CLM15                                
         W(21,12) = XM( 2) * CLM12                                      &
     &             +XM( 8) * CLM12*CLM21                                
         W(36,12) = XM( 2) * CLM12                                      &
     &             +XM( 8) * CLM12*CLM36                                
         W(20,13) = XM(11) * CLM13*CLM20 
         W(34,13) = XM(11) * CLM13 
         W(20,17) = XM(11) * CLM20 
         W(34,17) = XM(11) 
         W(15,18) = XM( 2)                                              &
     &             +XM( 8) * CLM15                                      
         W(21,18) = XM( 2)                                              &
     &             +XM( 8) * CLM21                                      
         W(28,18) = XM( 2)                                              &
     &             +XM( 8) * CLM28                                      
         W(36,18) = XM( 2)                                              &
     &             +XM( 8) * CLM36                                      
         W(20,31) = XM(11) * CLM20 
         W(34,31) = XM(11) 
         W(43,31) = XM(11) 
         W(28,12) = W(21,12) 
         W(45,12) = W(36,12) 
         W(43,13) = W(34,13) 
         W(26,14) = W(20,13) 
         W(35,14) =-W(34,13) 
         W(42,14) = W(34,13) 
         W(43,17) = W(34,17) 
         W(45,18) = W(36,18) 
         W(41,19) = W(33,18) 
         W(26,23) = W(20,17) 
         W(35,23) =-W(34,17) 
         W(42,23) = W(34,17) 
         W(27,24) = W(27,19) 
         W(41,24) = W(33,18) 
         W(15,25) = W(15,18) 
         W(21,25) = W(28,18) 
         W(28,25) = W(21,18) 
         W(33,25) =-W(33,18) 
         W(36,25) = W(36,18) 
         W(45,25) = W(36,18) 
         W(21,30) = W(27,19) 
         W(28,30) =-W(27,19) 
         W(33,30) = W(33,18) 
         W(26,32) =-W(20,31) 
         W(35,32) = W(34,31) 
         W(42,32) =-W(43,31) 
         W(27,38) = W(27,19) 
         W(41,38) = W(33,18) 
         W(26,39) = W(20,31) 
         W(35,39) =-W(43,31) 
         W(42,39) = W(34,31) 
         W(20,40) = W(20,31) 
         W(34,40) = W(43,31) 
         W(43,40) = W(34,31) 
      ENDIF 
      GO TO 280 
! *** (DD,KL), LIJ=6.                                                   
  260 CONTINUE 
!     KL=SS.                                                            
      I      = LX(1,LIJ) 
      XM(1) =  X(I+1) 
      XM(7) = (X(I+2)-X(I+3)*TWO+X(I+4))*fourth 
      W( 1,15) = XM( 1)                                                 &
     &          +XM( 7) * CLM15                                         
      W( 1,21) = XM( 1)                                                 &
     &          +XM( 7) * CLM21                                         
      W( 1,36) = XM( 1)                                                 &
     &          +XM( 7) * CLM36                                         
      W( 1,28) = W( 1,21) 
      W( 1,45) = W( 1,36) 
      IF(LKLMAX.GT.1) THEN 
!        KL=PS.                                                         
         I   = LX(2,LIJ) 
         XM(3) = (X(I+1)-X(I+2))*half 
         XM(9) =-(X(I+3)-X(I+4)*TWO+X(I+5)                              &
     &           -X(I+6)+X(I+7)*TWO-X(I+8))/8.0D0                       
         XM(12)=-(X(I+9)-X(I+10)-X(I+11)+X(I+12))*fourth 
         W( 2,15) = XM( 3)                                              &
     &             +XM( 9) * CLM15                                      
         W( 4,20) = XM(12) * CLM20 
         W( 2,21) = XM( 3)                                              &
     &             +XM( 9) * CLM21                                      
         W( 4,34) = XM(12) 
         W( 2,36) = XM( 3)                                              &
     &             +XM( 9) * CLM36                                      
         W( 7,26) = W( 4,20) 
         W( 2,28) = W( 2,21) 
         W( 7,35) =-W( 4,34) 
         W( 7,42) = W( 4,34) 
         W( 4,43) = W( 4,34) 
         W( 2,45) = W( 2,36) 
!        KL=PP.                                                         
         I   = LX(3,LIJ) 
         XM(1) =  X(I+1) 
         XM(6) = (X(I+2)-X(I+3)*TWO+X(I+4))*fourth 
         XM(7) = (X(I+5)-X(I+6)*TWO+X(I+7))*fourth 
         ZZZZ  =  X(I+8)+X(I+9)+X(I+10)+X(I+11) - TWO*                  &
     &           (X(I+12)+X(I+13)+X(I+14)+X(I+15)-X(I+16)-X(I+17))      
         XYXY  =  X(I+16)+X(I+17)-X(I+18)*TWO 
         XM(10)= (ZZZZ-XYXY)/16.0D0 
         XM(13)= (X(I+19)-X(I+20)-X(I+21)+X(I+22)                       &
     &           -X(I+23)+X(I+24)+X(I+25)-X(I+26))/8.0D0                
         XM(14)= XYXY*fourth 
         W( 3,15) = XM( 1)                                              &
     &             +XM( 6) * CLM3                                       &
     &             +XM( 7) * CLM15                                      &
     &             +XM(10) * CLM15*CLM3                                 
         W( 6,15) = XM( 1)                                              &
     &             +XM( 6) * CLM6                                       &
     &             +XM( 7) * CLM15                                      &
     &             +XM(10) * CLM15*CLM6                                 
         W( 5,20) = XM(13) * CLM20 
         W( 3,21) = XM( 1)                                              &
     &             +XM( 6) * CLM3                                       &
     &             +XM( 7) * CLM21                                      &
     &             +XM(10) * CLM21*CLM3                                 
         W( 6,21) = XM( 1)                                              &
     &             +XM( 6) * CLM6                                       &
     &             +XM( 7) * CLM21                                      &
     &             +XM(10) * CLM21*CLM6                                 &
     &             +XM(14)                                              
         W(10,21) = XM( 1)                                              &
     &             +XM( 6) * CLM10                                      &
     &             +XM( 7) * CLM21                                      &
     &             +XM(10) * CLM21*CLM10                                &
     &             -XM(14)                                              
         W( 9,27) = XM(14) 
         W( 6,33) = XM(14) * CLM33 
         W( 5,34) = XM(13) 
         W( 3,36) = XM( 1)                                              &
     &             +XM( 6) * CLM3                                       &
     &             +XM( 7) * CLM36                                      &
     &             +XM(10) * CLM36*CLM3                                 
         W( 6,36) = XM( 1)                                              &
     &             +XM( 6) * CLM6                                       &
     &             +XM( 7) * CLM36                                      &
     &             +XM(10) * CLM36*CLM6                                 
         W(10,15) = W( 6,15) 
         W( 8,26) = W( 5,20) 
         W( 3,28) = W( 3,21) 
         W( 6,28) = W(10,21) 
         W(10,28) = W( 6,21) 
         W(10,33) =-W( 6,33) 
         W( 8,35) =-W( 5,34) 
         W(10,36) = W( 6,36) 
         W( 9,41) = W( 6,33) 
         W( 8,42) = W( 5,34) 
         W( 5,43) = W( 5,34) 
         W( 3,45) = W( 3,36) 
         W( 6,45) = W( 6,36) 
         W(10,45) = W( 6,36) 
      ENDIF 
      IF(LKLMAX.GT.3) THEN 
!        KL=DS.                                                         
         I   = LX(4,LIJ) 
         XM(6) = (X(I+1)-X(I+2)*TWO+X(I+3))*fourth 
         ZZZZ  =  X(I+4)+X(I+5)+X(I+6)+X(I+7) - TWO*                    &
     &           (X(I+8)+X(I+9)+X(I+10)+X(I+11)-X(I+12)-X(I+13))        
         XYXY  =  X(I+12)+X(I+13)-X(I+14)*TWO 
         XM(10)= (ZZZZ-XYXY)/16.0D0 
         XM(13)= (X(I+15)-X(I+16)-X(I+17)+X(I+18)                       &
     &           -X(I+19)+X(I+20)+X(I+21)-X(I+22))/8.0D0                
         XM(14)= XYXY*fourth 
         W(11,15) = XM( 6) * CLM11                                      &
     &             +XM(10) * CLM15*CLM11                                
         W(16,20) = XM(13) * CLM20 
         W(11,21) = XM( 6) * CLM11                                      &
     &             +XM(10) * CLM21*CLM11                                
         W(29,21) = XM(14) 
         W(29,33) = XM(14) * CLM33 
         W(16,34) = XM(13) 
         W(11,36) = XM( 6) * CLM11                                      &
     &             +XM(10) * CLM36*CLM11                                
         W(22,26) = W(16,20) 
         W(37,27) = W(29,21) 
         W(11,28) = W(11,21) 
         W(29,28) =-W(29,21) 
         W(22,35) =-W(16,34) 
         W(37,41) = W(29,33) 
         W(22,42) = W(16,34) 
         W(16,43) = W(16,34) 
         W(11,45) = W(11,36) 
!        KL=DP.                                                         
         I   = LX(5,LIJ) 
         XM(3) = (X(I+1)-X(I+2))*half 
         XM(9) =-(X(I+3)-X(I+4)*TWO+X(I+5)                              &
     &           -X(I+6)+X(I+7)*TWO-X(I+8))/8.0D0                       
         XM(12)=-(X(I+9)-X(I+10)-X(I+11)+X(I+12))*fourth 
         W(12,15) = XM( 3) * CLM12                                      &
     &             +XM( 9) * CLM15*CLM12                                
         W(18,15) = XM( 3)                                              &
     &             +XM( 9) * CLM15                                      
         W(13,20) = XM(12) * CLM20*CLM13 
         W(17,20) = XM(12) * CLM20 
         W(31,20) = XM(12) * CLM20 
         W(12,21) = XM( 3) * CLM12                                      &
     &             +XM( 9) * CLM21*CLM12                                
         W(18,21) = XM( 3)                                              &
     &             +XM( 9) * CLM21                                      
         W(25,21) = XM( 3)                                              &
     &             +XM( 9) * CLM21                                      
         W(13,34) = XM(12) * CLM13 
         W(17,34) = XM(12) 
         W(31,34) = XM(12) 
         W(40,34) = XM(12) 
         W(12,36) = XM( 3) * CLM12                                      &
     &             +XM( 9) * CLM36*CLM12                                
         W(18,36) = XM( 3)                                              &
     &             +XM( 9) * CLM36                                      
         W(25,15) = W(18,15) 
         W(40,20) = W(31,20) 
         W(14,26) = W(13,20) 
         W(23,26) = W(17,20) 
         W(32,26) =-W(31,20) 
         W(39,26) = W(31,20) 
         W(19,27) = W(30,21) 
         W(24,27) = W(30,21) 
         W(38,27) = W(30,21) 
         W(12,28) = W(12,21) 
         W(18,28) = W(25,21) 
         W(25,28) = W(18,21) 
         W(30,28) =-W(30,21) 
         W(25,33) =-W(18,33) 
         W(30,33) = W(18,33) 
         W(14,35) =-W(13,34) 
         W(23,35) =-W(17,34) 
         W(32,35) = W(31,34) 
         W(39,35) =-W(40,34) 
         W(25,36) = W(18,36) 
         W(19,41) = W(18,33) 
         W(24,41) = W(18,33) 
         W(38,41) = W(18,33) 
         W(14,42) = W(13,34) 
         W(23,42) = W(17,34) 
         W(32,42) =-W(40,34) 
         W(39,42) = W(31,34) 
         W(13,43) = W(13,34) 
         W(17,43) = W(17,34) 
         W(31,43) = W(40,34) 
         W(40,43) = W(31,34) 
         W(12,45) = W(12,36) 
         W(18,45) = W(18,36) 
         W(25,45) = W(18,36) 
!        KL=DD.                                                         
         I   = LX(6,LIJ) 
         XM(1) =  X(I+1) 
         XM(6) = (X(I+2)-X(I+3)*TWO+X(I+4))*fourth 
         XM(7) = (X(I+5)-X(I+6)*TWO+X(I+7))*fourth 
         ZZZZ  =  X(I+8)+X(I+9)+X(I+10)+X(I+11) - TWO*                  &
     &           (X(I+12)+X(I+13)+X(I+14)+X(I+15)-X(I+16)-X(I+17))      
         XYXY  =  X(I+16)+X(I+17)-X(I+18)*TWO 
         XM(10)= (ZZZZ-XYXY)/16.0D0 
         XM(13)= (X(I+19)-X(I+20)-X(I+21)+X(I+22)                       &
     &           -X(I+23)+X(I+24)+X(I+25)-X(I+26))/8.0D0                
         XM(14)= XYXY*fourth 
         W(15,15) = XM( 1)                                              &
     &             +XM( 6) * CLM15                                      &
     &             +XM( 7) * CLM15                                      &
     &             +XM(10) * CLM15*CLM15                                
         W(21,15) = XM( 1)                                              &
     &             +XM( 6) * CLM21                                      &
     &             +XM( 7) * CLM15                                      &
     &             +XM(10) * CLM15*CLM21                                
         W(36,15) = XM( 1)                                              &
     &             +XM( 6) * CLM36                                      &
     &             +XM( 7) * CLM15                                      &
     &             +XM(10) * CLM15*CLM36                                
         W(20,20) = XM(13) * CLM20*CLM20 
         W(34,20) = XM(13) * CLM20 
         W(15,21) = XM( 1)                                              &
     &             +XM( 6) * CLM15                                      &
     &             +XM( 7) * CLM21                                      &
     &             +XM(10) * CLM21*CLM15                                
         W(21,21) = XM( 1)                                              &
     &             +XM( 6) * CLM21                                      &
     &             +XM( 7) * CLM21                                      &
     &             +XM(10) * CLM21*CLM21                                &
     &             +XM(14)                                              
         W(28,21) = XM( 1)                                              &
     &             +XM( 6) * CLM28                                      &
     &             +XM( 7) * CLM21                                      &
     &             +XM(10) * CLM21*CLM28                                &
     &             -XM(14)                                              
         W(33,21) = XM(14) * CLM33 
         W(36,21) = XM( 1)                                              &
     &             +XM( 6) * CLM36                                      &
     &             +XM( 7) * CLM21                                      &
     &             +XM(10) * CLM21*CLM36                                
         W(27,27) = XM(14) 
         W(21,33) = XM(14) * CLM33 
         W(33,33) = XM(14) * CLM33*CLM33 
         W(20,34) = XM(13) * CLM20 
         W(34,34) = XM(13) 
         W(43,34) = XM(13) 
         W(15,36) = XM( 1)                                              &
     &             +XM( 6) * CLM15                                      &
     &             +XM( 7) * CLM36                                      &
     &             +XM(10) * CLM36*CLM15                                
         W(21,36) = XM( 1)                                              &
     &             +XM( 6) * CLM21                                      &
     &             +XM( 7) * CLM36                                      &
     &             +XM(10) * CLM36*CLM21                                
         W(36,36) = XM( 1)                                              &
     &             +XM( 6) * CLM36                                      &
     &             +XM( 7) * CLM36                                      &
     &             +XM(10) * CLM36*CLM36                                
         W(45,36) = XM( 1)                                              &
     &             +XM( 6) * CLM45                                      &
     &             +XM( 7) * CLM36                                      &
     &             +XM(10) * CLM36*CLM45                                
         W(28,15) = W(21,15) 
         W(45,15) = W(36,15) 
         W(43,20) = W(34,20) 
         W(45,21) = W(36,21) 
         W(26,26) = W(20,20) 
         W(35,26) =-W(34,20) 
         W(42,26) = W(34,20) 
         W(41,27) = W(33,21) 
         W(15,28) = W(15,21) 
         W(21,28) = W(28,21) 
         W(28,28) = W(21,21) 
         W(33,28) =-W(33,21) 
         W(36,28) = W(36,21) 
         W(45,28) = W(36,21) 
         W(28,33) =-W(21,33) 
         W(26,35) =-W(20,34) 
         W(35,35) = W(34,34) 
         W(42,35) =-W(43,34) 
         W(28,36) = W(21,36) 
         W(27,41) = W(21,33) 
         W(41,41) = W(33,33) 
         W(26,42) = W(20,34) 
         W(35,42) =-W(43,34) 
         W(42,42) = W(34,34) 
         W(20,43) = W(20,34) 
         W(34,43) = W(43,34) 
         W(43,43) = W(34,34) 
         W(15,45) = W(15,36) 
         W(21,45) = W(21,36) 
         W(28,45) = W(21,36) 
         W(36,45) = W(45,36) 
         W(45,45) = W(36,36) 
      ENDIF 
  280 END DO 
!     *                                                                 
! *** INCLUDE TWO-ELECTRON INTEGRALS INVOLVING ONLY S AND P ORBITALS.   
!     OPTION MODESP=1. LOCAL INTEGRALS RI(22) FROM SUBROUTINE REPP.     
!     *                                                                 
      IF(MODESP.EQ.1) THEN 
         W( 1, 1) = RI(1) 
         IF(LKLMAX.GT.1) THEN 
            W( 2, 1) = RI( 5) 
            W( 3, 1) = RI(11) 
            W( 6, 1) = RI(12) 
            W(10, 1) = RI(12) 
         ENDIF 
         IF(LIJMAX.GT.1) THEN 
            W( 1, 2) = RI( 2) 
            W( 1, 3) = RI( 3) 
            W( 1, 6) = RI( 4) 
            W( 1,10) = RI( 4) 
         ENDIF 
         IF(LIJMAX.GT.1 .AND. LKLMAX.GT.1) THEN 
            W( 2, 2) = RI( 6) 
            W( 3, 2) = RI(13) 
            W( 6, 2) = RI(14) 
            W(10, 2) = RI(14) 
            W( 2, 3) = RI( 8) 
            W( 3, 3) = RI(16) 
            W( 6, 3) = RI(18) 
            W(10, 3) = RI(18) 
            W( 4, 4) = RI( 7) 
            W( 5, 4) = RI(15) 
            W( 4, 5) = RI(10) 
            W( 5, 5) = RI(20) 
            W( 2, 6) = RI( 9) 
            W( 3, 6) = RI(17) 
            W( 6, 6) = RI(19) 
            W(10, 6) = RI(21) 
            W( 7, 7) = RI( 7) 
            W( 8, 7) = RI(15) 
            W( 7, 8) = RI(10) 
            W( 8, 8) = RI(20) 
            W( 9, 9) = RI(22) 
            W( 2,10) = RI( 9) 
            W( 3,10) = RI(17) 
            W( 6,10) = RI(21) 
            W(10,10) = RI(19) 
         ENDIF 
      ENDIF 
!     *                                                                 
! *** CORE-ELECTRON ATTRACTION INTEGRALS.                               
!     CORE(I,1) - ELECTRONS (IJ) AT ATOM A AND CORE OF ATOM B (NC=NJ).  
!     CORE(I,2) - ELECTRONS (KL) AT ATOM B AND CORE OF ATOM A (NC=NI).  
!     *                                                                 
!     INDEX I IN CORE(I,1) AND CORE(I,2).                               
!     (S S )=1, (S P )=2, (P P )=3, (P+P+)=4                            
!     (D S )=5, (D P )=6, (D D )=7, (D+P+)=8, (D+D+)=9, (D#D#)=10.      
!     NOTATION: P-SIGMA=P , P-PI=P+, D-SIGMA=D , D-PI=D+, D-DELTA=D#.   
!     *                                                                 
      IF(MODESP.EQ.1) THEN 
         IFIRST = 8 
      ELSE 
         IFIRST = 1 
      ENDIF 
!     LOOP OVER THE TWO POSSIBLE CASES.                                 
!     NC   ATOMIC NUMBER OF ATOM CARRYING THE CORE.                     
!     NE   ATOMIC NUMBER OF ATOM CARRYING THE ELECTRONS.                
!     LE   LIMIT FOR L PAIR INDEX, 1 FOR S, 3 FOR SP, 6 FOR SPD.        
!     SIG  SIGN FACTOR FOR INTEGRALS, ALWAYS POSITIVE EXCEPT FOR        
!          INTERACTIONS BETWEEN (IJ) AT ATOM A INVOLVING ONE P-SIGMA    
!          OR D-PI ORBITAL AND THE CORE OF ATOM B (CASE N=1).           
      DO 380 N=1,2 
      IF(N.EQ.1) THEN 
         NC  = NJ 
         NE  = NI 
         LE  = LIJMAX 
         SIG =-ONE 
      ELSE 
         NC  = NI 
         NE  = NJ 
         LE  = LKLMAX 
         SIG = ONE 
      ENDIF 
!     CHECK FOR EXTERNAL POINT CHARGE.                                  
      IF(NE.EQ.0) GO TO 380 
      IF(NC.EQ.0) THEN 
         TORENC = ONE 
      ELSE 
         TORENC = TORE(NC) 
      ENDIF 
!     ADDITIVE TERM FOR CORE.                                           
      QA     = PO(9,NC) 
!     CHECK FOR NEGLECT OF PENETRATION INTEGRALS.                       
!     IF THE ADDITIVE TERMS FOR SS (PO(1,NC)) AND THE CORE (PO(9,NC))   
!     ARE EQUAL, THE CORE-ELECTRON ATTRACTION INTEGRALS ARE EXPRESSED   
!     IN TERMS OF THE TWO-ELECTRON INTEGRALS (SS,KL) AND (IJ,SS).       
!     LCW(I) IS THE STANDARD PAIR INDEX CORRESPONDING TO I IN CORE(I,N).
      IF(ABS(QA-PO(1,NC)).LT.SMALL) THEN 
         IF(LE.EQ.1) THEN 
            ICMAX = 1 
         ELSE IF(LE.EQ.3) THEN 
            ICMAX = 4 
         ELSE 
            ICMAX = 10 
         ENDIF 
         IF(N.EQ.1) THEN 
            DO 310 I=1,ICMAX 
            CORE(I,1) = -TORENC*W(1,LCW(I)) 
  310       CONTINUE 
         ELSE 
            DO 320 I=1,ICMAX 
            CORE(I,2) = -TORENC*W(LCW(I),1) 
  320       CONTINUE 
         ENDIF 
         GO TO 380 
      ENDIF 
!     EVALUATE RELEVANT DISTANCES BETWEEN POINT CHARGES.                
      IF(MODESP.EQ.1) GO TO 330 
!     SS.                                                               
      QB     = PO(1,NE) 
      ADD    = (QA+QB)**2 
!     L1=0, L2=0, M=0, LABEL=1.                                         
      X(1)   = R2+ADD 
      ILAST  = 1 
      IF(LE.GT.1) THEN 
!        PS.                                                            
         DB     = DD(2,NE) 
         QB     = PO(2,NE) 
         ADD    = (QA+QB)**2 
!        L1=0, L2=1, M=0, LABEL=3.                                      
         X(2)   = (R+DB)**2+ADD 
         X(3)   = (R-DB)**2+ADD 
!        PP.                                                            
         DB     = DD(3,NE)*SQRT2 
         QB     = PO(3,NE) 
         QBS    = PO(7,NE) 
         ADD    = (QA+QB)**2 
!        L1=0, L2=0, M=0, LABEL=1.                                      
         X(4)   = R2+(QA+QBS)**2 
!        L1=0, L2=2, M=0, LABEL=6.                                      
         X(5)   = (R-DB)**2+ADD 
         X(6)   = R2+DB**2+ADD 
         X(7)   = (R+DB)**2+ADD 
         ILAST  = 7 
      ENDIF 
  330 CONTINUE 
      IF(LE.GT.3) THEN 
!        DS.                                                            
         DB     = DD(4,NE) 
         QB     = PO(4,NE) 
         ADD    = (QA+QB)**2 
!        L1=0, L2=2, M=0, LABEL=6.                                      
         X(8)   = (R-DB)**2+ADD 
         X(9)   = R2+DB**2+ADD 
         X(10)  = (R+DB)**2+ADD 
!        DP.                                                            
         DB     = DD(5,NE) 
         QB     = PO(5,NE) 
         ADD    = (QA+QB)**2 
!        L1=0, L2=1, M=0, LABEL=3.                                      
         X(11)  = (R+DB)**2+ADD 
         X(12)  = (R-DB)**2+ADD 
!        DD.                                                            
         DB     = DD(6,NE) 
         QB     = PO(6,NE) 
         QBS    = PO(8,NE) 
         ADD    = (QA+QB)**2 
!        L1=0, L2=0, M=0, LABEL=1.                                      
         X(13)  = R2+(QA+QBS)**2 
!        L1=0, L2=2, M=0, LABEL=6.                                      
         X(14)  = (R-DB)**2+ADD 
         X(15)  = R2+DB**2+ADD 
         X(16)  = (R+DB)**2+ADD 
         ILAST  = 16 
      ENDIF 
!     CALCULATE INVERSE DISTANCES.                                      
!     THE NUMERATOR IS CHOSEN SUCH THAT THE NUCLEAR CHARGES ARE         
!     INCLUDED ALREADY AT THIS POINT, AS WELL AS THE CONVERSION TO EV.  
      FACTOR = -TORENC*AU_TO_EV 
      DO 340 I=IFIRST,ILAST 
  340 X(I)   = FACTOR/SQRT(X(I)) 
!     EVALUATE SEMIEMPIRICAL MULTIPOLE-MULTIPOLE INTERACTIONS XM(LABEL).
!     EVALUATE THE UNIQUE CORE-ELECTRON ATTRACTION INTEGRALS CORE(I,N)  
!     FROM XM(LABEL).                                                   
      IF(MODESP.EQ.1) GO TO 350 
!     SS.                                                               
      XM(1)  = X(1) 
      CORE(1,N) = XM(1) 
      IF(LE.GT.1) THEN 
!        PS.                                                            
         XM(3) = (X(2)-X(3))*half 
         CORE(2,N) = XM(3) * SIG 
!        PP.                                                            
         XM(1) =  X(4) 
         XM(6) = (X(5)-X(6)*TWO+X(7))*fourth 
         CORE(3,N) = XM(1)                                              &
     &              +XM(6) * CLM3                                       
         CORE(4,N) = XM(1)                                              &
     &              +XM(6) * CLM6                                       
      ENDIF 
  350 CONTINUE 
      IF(LE.GT.3) THEN 
!        DS.                                                            
         XM(6) = (X(8)-X(9)*TWO+X(10))*fourth 
         CORE(5,N) = XM(6) * CLM11 
!        DP.                                                            
         XM(3) = (X(11)-X(12))*half 
         CORE(6,N) = XM(3) * CLM12 * SIG 
         CORE(8,N) = XM(3) * SIG 
!        DD.                                                            
         XM(1) =  X(13) 
         XM(6) = (X(14)-X(15)*TWO+X(16))*fourth 
         CORE(7,N) = XM(1)                                              &
     &              +XM(6) * CLM15                                      
         CORE(9,N) = XM(1)                                              &
     &              +XM(6) * CLM21                                      
         CORE(10,N)= XM(1)                                              &
     &              +XM(6) * CLM36                                      
      ENDIF 
  380 END DO 
!     WRITE(6,'('' CORE:'',/(2X,10F12.5))') CORE                        
      RETURN 
      END subroutine qm2_repp_d                                          
