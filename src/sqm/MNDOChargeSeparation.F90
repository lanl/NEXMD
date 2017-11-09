#include "dprec.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! one-center two electron integrals calculated from Slater-Condon parameters
! This modules is a wrapper for Thiel's routines
!
! by: Taisung Lee (Rutgers, 2011)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module MNDOChargeSeparation

public GetDDAndPho, GetOneCenter2Electron

private AIJL, POIJ


contains

subroutine GetDDAndPho(qmmm_struct, qmtype, DD, PO)

!C     *
!C     CALCULATION OF CHARGE SEPARATIONS AND ADDITIVE TERMS USED
!C     TO COMPUTE THE TWO-CENTER TWO-ELECTRON INTEGRALS IN MNDO/D.
!C     NI IS THE ATOMIC NUMBER OF THE CURRENT ELEMENT.
!C     *
!C     IN THE CASE OF AN SP BASIS, THE CHARGE SEPARATIONS AND ADDITIVE
!C     TERMS ARE EQUIVALENT IN MNDO/d AND MNDO.
!C     *
!C     CHARGE SEPARATIONS DD(6,LMZ) FROM FORMULA FUNCTION AIJL.
!C     ADDITIVE TERMS     PO(9,LMZ) FROM FUNCTION POIJ.
!C     SEE EQUATIONS (12)-(16) OF TCA PAPER FOR DD.
!C     SEE EQUATIONS (19)-(26) OF TCA PAPER FOR PO.
!C     *
!C     INDEX OF DD AND PO : SS 1, SP 2, PP 3, SD 4, PD 5, DD 6.
!C     MULTIPOLE          :  L=0,  L=1,  L=2,  L=2,  L=1,  L=2.
!C     SPECIAL INDEX OF PO: PP 7, DD 8.
!C     MULTIPOLE          :  L=0,  L=0.
!C     FOR ATOMIC CORE    : ADDITIVE TERM PO(9,NI) 
!C     *
!C     THE HYBRIDIZATION CONTRIBUTIONS TO THE DIPOLE MOMENT CONTAIN
!C     FACTORS WHICH DEPEND ON THE CHARGE SEPARATIONS AS FOLLOWS.
!C     HYFSP(NI) = DD(2,NI)*5.0832D0
!C     HYFPD(NI) = DD(5,NI)*5.0832D0
!C     THESE FACTORS ARE COMPUTED AND SAVED IN COMMON BLOCK DIPOL1.
!C     *

  use constants  , only : zero, half, one, two, three, AU_TO_EV
  use qmmm_module, only : qm2_params
  use qmmm_struct_module, only : qmmm_struct_type 
  implicit none

  type(qmmm_struct_type), intent(in) :: qmmm_struct
  integer, intent(in)::qmtype
  _REAL_, intent(out)::DD(6), PO(9)
  
   !local 
  integer, save:: qmtype_saved=-1
  _REAL_, save::DD_saved(6), PO_saved(9)
  logical, save::initialized=.false. 
  
  integer::NS, NP, ND, atomic_number
  _REAL_::ZS,ZP,ZD
  _REAL_::D, FG, FG1, FG2, AIJ22, AIJ43, AIJ52, AIJ63
  _REAL_::hyf, dipfac, hyfpd
  
      atomic_number=qmmm_struct%qm_type_id(qmtype)

      NS     = qm2_params%sp_quantum_number(qmtype)
      NP     = max(2,qm2_params%sp_quantum_number(qmtype)      )
      ND     = qm2_params%d_quantum_number(qmtype)
      ZS     = qm2_params%s_orb_exp_by_type(qmType)
      ZP     = qm2_params%p_orb_exp_by_type(qmType)
      ZD     = qm2_params%d_orb_exp_by_type(qmType)

!C
!
!C     AIJ-VALUES ARE COMPUTED AS DEFINED IN EQUATION (7) OF TCA PAPER
!C     AND THEN STORED AS VARIABLES AIJ.. WITH TWO DIGITS AT THE END.
!C     FIRST  DIGIT: 1 SS, 2 SP, 3 PP, 4 SD, 5 PD, 6 DD.
!C     SECOND DIGIT: L+1 FROM DEFINITION OF MULTIPOLE.
!C



    if ((qmtype_saved.ne.qmtype) .or. (.not.initialized)) then
       
!C *** SECTION FOR AN S BASIS.
!C     THERE IS ONLY ONE ADDITIVE TERM.   

      PO(1) = half*AU_TO_EV/qm2_params%gss(qmtype)
!      IF(atomic_number.EQ.1) RETURN
!C
!C *** SECTION FOR AN SP BASIS SET.
!C     CHARGE SEPARATIONS AND ADDITIVE TERMS MUST BE COMPUTED.


      AIJ22    = AIJL(ZS,ZP,NS,NP,1)
      DD(2) = AIJ22/SQRT(THREE)
      PO(2) = POIJ(1,DD(2), qm2_params%hsp(qmtype))

      DD(3) = SQRT((2*NP+1)*(2*NP+2)/20.0D0) / ZP
      PO(3) = POIJ(2,DD(3)*SQRT(TWO),qm2_params%hpp(qmtype))

      PO(7) = PO(1)
      HYF  = DD(2)*DIPFAC
      HYFPD= ZERO
      
!C
!C *** SECTION FOR AN SPD BASIS SET.
!C     ADDITIONAL TERMS INVOLVING D ORBITALS.
!C     NOTE THE EXTRA FACTOR OF SQRT2 FOR THE CHARGE SEPARATIONS
!C     DD(I) WITH I=4,6, WHICH REFER TO SQUARE QUADRUPOLES,
!C     FOR SIMPLIFICATION OF THE CODE IN REPPD.
      AIJ52    = AIJL(ZP,ZD,NP,ND,1)
      AIJ43    = AIJL(ZS,ZD,NS,ND,2)
      AIJ63    = AIJL(ZD,ZD,ND,ND,2)
!     SD

      D        = SQRT(AIJ43*SQRT(ONE/15.0D0))*SQRT(TWO)
      DD(4) = D
      FG=GetOneCenter2Electron(qmtype, 19)
      PO(4) = POIJ(2,D,FG)

!C     PD
      D        = AIJ52/SQRT(5.0D0)
!C     PREVIOUS STATEMENT AS IN THE TCA PAPER,
!C     NEXT STATEMENT AS A POSSIBLE ALTERNATIVE.
!C     FG       = REPD(33)-1.6D0*REPD(35) 
      DD(5) = D
      FG=GetOneCenter2Electron(qmtype, 23)
      FG1=GetOneCenter2Electron(qmtype, 35)
      PO(5) = POIJ(1,D,FG-1.8D0*FG1)
      HYFPD= D*DIPFAC
      
!C     DD
      FG=GetOneCenter2Electron(qmtype, 29)
      FG1=GetOneCenter2Electron(qmtype, 30)
      FG2=GetOneCenter2Electron(qmtype, 31)     
      PO(8) = POIJ(0,ONE,0.2D0*(FG+TWO*FG1+TWO*FG2))
!CDEC  NEXT TWO STATEMENTS RUN INTO COMPILER BUG ON DEC ALPHA (OSF1).
!C     THEY ARE REPLACED BY TWO MATHEMATICALLY EQUIVALENT STATEMENTS.
!C     D        = SQRT(AIJ63/7.0D0)
!C     D        = D*SQRT(TWO)
      D        = AIJ63/7.0D0
      D        = SQRT(TWO*D)
      DD(6) = D
      FG=GetOneCenter2Electron(qmtype, 44)      
      FG1=GetOneCenter2Electron(qmtype, 52)
      PO(6) = POIJ(2,D,FG-(20.0D0/35.0D0)*FG1)
      
      !  this one was added in param.f in Thiel's mndo implementation
      PO(9)=PO(1)
      
      qmtype_saved=qmtype
      DD_saved=DD
      PO_saved=PO
      initialized=.true.
      
    end if

    DD=DD_saved
    PO=PO_saved    
    

end subroutine

function AIJL(Z1,Z2,N1,N2,L)
!C *** GENERAL FORMULA FOR AIJ-VALUES (SEE BELOW).
!C     SEE EQUATION (7) OF TCA PAPER.

  use constants, only : FC

  implicit none

  _REAL_::AIJL
  _REAL_, intent(in)::Z1,Z2
  integer, intent(in)::N1, N2, L
  
  AIJL=FC(N1+N2+L+1) / SQRT( FC(2*N1+1)*FC(2*N2+1) )  &
        *  ( 2*Z1/(Z1+Z2) ) ** N1 * SQRT( 2*Z1/(Z1+Z2) )  &
        *  ( 2*Z2/(Z1+Z2) ) ** N2 * SQRT( 2*Z2/(Z1+Z2) )  &
        /  (Z1+Z2) ** L
        
end function AIJL

function GetOneCenter2Electron(qmType, index) result(integral)

! The one-center two-electron integrals should be calculated only
! once for each atom type--Need to be changed later..
! Taisung Lee, Rutgers, 2011
!
   use constants, only : half, one, two, three, four
   use SlaterOverlap
   use qmmm_module, only :  qm2_params
   
   implicit none
   
   integer, intent(in)::qmType, index
   _REAL_::integral
   
   !local
   
   _REAL_, parameter::S3 =0.17320508075689D+01
   _REAL_, parameter::S5 =0.22360679774998D+01   
   _REAL_, parameter::S15=0.38729833462074D+01
    
   _REAL_, save::REPD(52)=0.0D0
   logical, save::initialized=.false.
   integer, save::qmType_saved=-1
 
   integer::NS, ND
   _REAL_::ES, EP, ED
   _REAL_::R016, R036, R066, R155, R125, R244, R236, R266, R234, R246
   _REAL_::R355, R466
   _REAL_::F0DD, F2DD, F0PD, F2FD, F4DD, F2PD, G1PD, G3PD
  
   if((qmType_saved.ne.qmType) .or. (.not.initialized)) then

      NS     = qm2_params%sp_quantum_number(qmtype)
      ND     = qm2_params%d_quantum_number(qmtype)
      ES     = qm2_params%s_orb_exp_tail_by_type(qmType)
      EP     = qm2_params%p_orb_exp_tail_by_type(qmType)
      ED     = qm2_params%d_orb_exp_tail_by_type(qmType)
! *** SLATER-CONDON PARAMETERS (Rlij).
!     FIRST  DIGIT (l)  L QUANTUM NUMBER OF SLATER-CONDON PARAMETER.
!     SECOND DIGIT (i)  SS 1, SP 2, PP 3, SD 4, PD 5, DD 6 - ELECTRON 1.
!     SECOND DIGIT (j)  SS 1, SP 2, PP 3, SD 4, PD 5, DD 6 - ELECTRON 2.
      R016   = GetSlaterCondonParameter(0,NS,ES,NS,ES,ND,ED,ND,ED)
      R036   = GetSlaterCondonParameter(0,NS,EP,NS,EP,ND,ED,ND,ED)
      R066   = GetSlaterCondonParameter(0,ND,ED,ND,ED,ND,ED,ND,ED)
      R155   = GetSlaterCondonParameter(1,NS,EP,ND,ED,NS,EP,ND,ED)        
      R125   = GetSlaterCondonParameter(1,NS,ES,NS,EP,NS,EP,ND,ED)        
      R244   = GetSlaterCondonParameter(2,NS,ES,ND,ED,NS,ES,ND,ED)        
      R236   = GetSlaterCondonParameter(2,NS,EP,NS,EP,ND,ED,ND,ED)        
      R266   = GetSlaterCondonParameter(2,ND,ED,ND,ED,ND,ED,ND,ED)        
      R234   = GetSlaterCondonParameter(2,NS,EP,NS,EP,NS,ES,ND,ED)        
      R246   = GetSlaterCondonParameter(2,NS,ES,ND,ED,ND,ED,ND,ED)        
      R355   = GetSlaterCondonParameter(3,NS,EP,ND,ED,NS,EP,ND,ED)        
      R466   = GetSlaterCondonParameter(4,ND,ED,ND,ED,ND,ED,ND,ED)        
!     SAVE SLATER-CONDON PARAMETERS.
      F0DD = R066
      F2DD = R266
      F4DD = R466
      F0PD = R036
      F2PD = R236
      G1PD = R155
      G3PD = R355
!AWG  For PM6 we need to modify F0SD and G2SD for some elements
      if ( abs(qm2_params%F0SD(qmType)) > 1.0D-09) then
         R016 = qm2_params%F0SD(qmType)
      end if
      if ( abs(qm2_params%G2SD(qmType)) > 1.0D-09) then
         R244 = qm2_params%G2SD(qmType)
      end if
!!     KEEP PREDEFINED SLATER-CONDON PARAMETERS (IF REQUESTED).
!      IF(IF0SD.GT.0) THEN
!         R016  = F0SD
!      ELSE
!         F0SD = R016
!      ENDIF
!      IF(IG2SD.GT.0) THEN
!         R244  = G2SD
!      ELSE
!         G2SD = R244
!      ENDIF
! *** COMPUTE ONE-CENTER TWO-ELECTRON INTEGRALS
!     FROM THE SLATER-CONDON PARAMETERS.
      REPD( 1) = R016
      REPD( 2) = (TWO/(THREE*S5))*R125
      REPD( 3) = (ONE/S15)*R125
      REPD( 4) = (TWO/(5.D0*S5))*R234
      REPD( 5) = R036 + (FOUR/35.D0)*R236
      REPD( 6) = R036 + (TWO /35.D0)*R236
      REPD( 7) = R036 - (FOUR/35.D0)*R236
      REPD( 8) = -(ONE/(THREE*S5))*R125
      REPD( 9) = SQRT(THREE/125.D0)*R234
      REPD(10) = (S3   /35.D0)*R236
      REPD(11) = (THREE/35.D0)*R236
      REPD(12) = -(0.2D0/S5)*R234
      REPD(13) = R036 - (TWO/35.D0)*R236
      REPD(14) = -(TWO*S3/35.D0)*R236
      REPD(15) = -REPD( 3)
      REPD(16) = -REPD(11)
      REPD(17) = -REPD( 9)
      REPD(18) = -REPD(14)
      REPD(19) = 0.2D0*R244
      REPD(20) = (TWO/(7.D0*S5))*R246
      REPD(21) =  REPD(20)*half
      REPD(22) = -REPD(20)
      REPD(23) = (FOUR  /15.D0)*R155 + (27.D0  /245.D0)*R355
      REPD(24) = (TWO*S3/15.D0)*R155 - (9.D0*S3/245.D0)*R355
      REPD(25) = (ONE/15.D0)*R155 + (18.D0   /245.D0)*R355
      REPD(26) = -(S3/15.D0)*R155 + (12.D0*S3/245.D0)*R355
      REPD(27) = -(S3/15.D0)*R155 - (THREE*S3/245.D0)*R355
      REPD(28) = -REPD(27)
      REPD(29) = R066 + (FOUR/49.D0)*R266 + (FOUR / 49.D0)*R466
      REPD(30) = R066 + (TWO /49.D0)*R266 - (24.D0/441.D0)*R466
      REPD(31) = R066 - (FOUR/49.D0)*R266 + ( 6.D0/441.D0)*R466
      REPD(32) = SQRT(THREE/245.D0)*R246
      REPD(33) = 0.2D0*R155 + (24.D0/245.D0)*R355
      REPD(34) = 0.2D0*R155 - ( 6.D0/245.D0)*R355
      REPD(35) = (THREE/49.D0)*R355
      REPD(36) = (ONE/49.D0)*R266 + (30.D0  /441.D0)*R466
      REPD(37) = (S3 /49.D0)*R266 - (5.D0*S3/441.D0)*R466
      REPD(38) = R066 - (TWO/49.D0)*R266 -  (FOUR/441.D0)*R466
      REPD(39) = -(TWO*S3/49.D0)*R266 + (10.D0*S3/441.D0)*R466
      REPD(40) = -REPD(32)
      REPD(41) = -REPD(34)
      REPD(42) = -REPD(35)
      REPD(43) = -REPD(37)
      REPD(44) = (THREE/49.D0)*R266 + (20.D0/441.D0)*R466
      REPD(45) = -REPD(39)
      REPD(46) = 0.2D0*R155-(THREE/35.D0)*R355
      REPD(47) = -REPD(46)
      REPD(48) = (FOUR /49.D0)*R266 + (15.D0/441.D0)*R466
      REPD(49) = (THREE/49.D0)*R266 - ( 5.D0/147.D0)*R466
      REPD(50) = -REPD(49)
      REPD(51) = R066 + (FOUR/49.D0)*R266 - (34.D0/441.D0)*R466
      REPD(52) = (35.D0/441.D0)*R466
      
      qmType_saved=qmType
      initialized=.true.
      
   end if   
      
  integral=REPD(index)
  return
   
end function GetOneCenter2Electron 


FUNCTION POIJ (L,D,FG)

!C     *
!C     DETERMINE ADDITIVE TERMS RHO=POIJ FOR TWO-CENTER TWO-ELECTRON
!C     INTEGRALS FROM THE REQUIREMENT THAT THE APPROPRIATE ONE-CENTER
!C     TWO-ELECTRON INTEGRALS ARE REPRODUCED.
!C     *
!C     INPUT DATA (SEE EQUATIONS (19)-(26) IN TCA PAPER).
!C     L     L QUANTUM NUMBER OF ADDITIVE TERM.
!C     D     CHARGE SEPARATION.
!C     FG    REFERENCE ONE-CENTER INTEGRAL (OR SLATER-CONDON PARAMETER).
!C     *
!C     SPECIAL CONVENTION IN THE CASE L=2.
!C     THE INPUT VALUE OF D (AS DEFINED IN THE CALLING ROUTINE) EQUALS
!C     D2*SQRT(TWO), WITH D2 DEFINED IN EQUATIONS (14)-(16) OF THE TCA
!C     PAPER. THIS SPECIAL CONVENTION FOR L=2 SHOULD BE KEPT IN MIND
!C     WHEN COMPARING THE CODE BELOW WITH EQUATIONS (24)-(26) OF THE
!C     TCA PAPER. THE CONVENTION IS MOTIVATED BY CODE SIMPLIFICATIONS
!C     IN SUBROUTINE REPPD WHICH ARISE WHEN VALUES OF D2*SQRT(TWO) ARE
!C     STORED IN DD(4,NI) AND DD(6,NI).
!C     *
      use constants, only : fourth, half, one, two, AU_TO_EV
      
      implicit none
      
      _REAL_::POIJ
      _REAL_,intent(in)::D, FG
      integer,intent(in)::L
      
      _REAL_,PARAMETER::EPSIL=1.0D-08, G1=0.382D0, G2=0.618D0
      integer ,PARAMETER::NITER=100   

! local

      integer::i
      _REAL_::DSQ, EV4, EV8, A1, A2, Y1, Y2, DELTA, F1, F2

      IF(L.EQ.0) THEN
         POIJ = half*AU_TO_EV/FG
         RETURN
      ENDIF
! *** HIGHER TERMS.
      DSQ    = D*D
      EV4    = AU_TO_EV*fourth
      EV8    = AU_TO_EV/8.0D0
      A1     = 0.1D0
      A2     = 5.0D0
      DO 10 I=1,NITER
      DELTA  = A2-A1
      IF(DELTA.LT.EPSIL) GO TO 20
      Y1     = A1 + DELTA*G1
      Y2     = A1 + DELTA*G2
      IF(L.EQ.1) THEN
         F1  = (EV4*(ONE/Y1-ONE/SQRT(Y1**2+DSQ)) - FG)**2
         F2  = (EV4*(ONE/Y2-ONE/SQRT(Y2**2+DSQ)) - FG)**2
      ELSE IF(L.EQ.2) THEN
         F1  = (EV8*(ONE/Y1-TWO/SQRT(Y1**2+DSQ*half)   &
                           +ONE/SQRT(Y1**2+DSQ)) - FG)**2
         F2  = (EV8*(ONE/Y2-TWO/SQRT(Y2**2+DSQ*half)   &
                           +ONE/SQRT(Y2**2+DSQ)) - FG)**2
      ENDIF
      IF(F1.LT.F2) THEN
         A2  = Y2
      ELSE
         A1  = Y1
      ENDIF
  10  CONTINUE
!C     DEFINE ADDITIVE TERM AFTER CONVERGENCE OF ITERATIONS.
  20  IF(F1.GE.F2) THEN
         POIJ = A2
      ELSE
         POIJ = A1
      ENDIF

END function POIJ

end module MNDOChargeSeparation
