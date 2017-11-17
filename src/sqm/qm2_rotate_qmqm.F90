! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
subroutine qm2_rotate_qmqm(qm2_struct, qmmm_struct, loop_count,IQM, JQM,NI,NJ,XI,XJ,W,KI,&
                       RI, core)

!********************************************************
! Current routine maintained by: Ross Walker (TSRI, 2005)
!                           Andreas W. Goetz (SDSC 2010)
! Inlining and Optimising by:    Ross Walker (TSRI, 2005)
!********************************************************

!********************************************************
!   ROTATE CALCULATES THE TWO-PARTICLE INTERACTIONS.
!
!   ON INPUT
!             IQM    = qm atom I number (in 1 to nquant loop)
!             JQM    = qm atom J number in inner loop
!             NI     = ATOMIC NUMBER OF FIRST ATOM.
!             NJ     = ATOMIC NUMBER OF SECOND ATOM.
!             XI     = COORDINATE OF FIRST ATOM.
!             XJ     = COORDINATE OF SECOND ATOM.
!
! ON OUTPUT W      = ARRAY OF TWO-ELECTRON REPULSION INTEGRALS.
!                    (rotated into molecular frame)
!           core   = nuclear attraction integrals
!                    (in local coordinates, not rotated to molecular frame)
!
! *** THIS ROUTINE COMPUTES THE REPULSION AND NUCLEAR ATTRACTION
!     INTEGRALS OVER MOLECULAR-FRAME COORDINATES.  THE INTEGRALS OVER
!     LOCAL FRAME COORDINATES ARE EVALUATED BY SUBROUTINE REPP AND
!     STORED AS FOLLOWS (WHERE P-SIGMA = O,   AND P-PI = P AND P* )
!     IN RI
!     (SS/SS)=1,   (SO/SS)=2,   (OO/SS)=3,   (PP/SS)=4,   (SS/OS)=5,
!     (SO/SO)=6,   (SP/SP)=7,   (OO/SO)=8,   (PP/SO)=9,   (PO/SP)=10,
!     (SS/OO)=11,  (SS/PP)=12,  (SO/OO)=13,  (SO/PP)=14,  (SP/OP)=15,
!     (OO/OO)=16,  (PP/OO)=17,  (OO/PP)=18,  (PP/PP)=19,  (PO/PO)=20,
!     (PP/P*P*)=21,   (P*P/P*P)=22.
!
!***********************************************************************
      use constants  , only : one, A_TO_BOHRS, A2_TO_BOHRS2 
      use qmmm_module, only : qmmm_nml,  qm2_structure, qm2_params, &
                              qm2_rij_eqns, EXPONENTIAL_CUTOFF
      use qmmm_struct_module, only : qmmm_struct_type

      implicit none
!Passed in
      type(qm2_structure),intent(inout) :: qm2_struct
      type(qmmm_struct_type), intent(in) :: qmmm_struct
      integer, intent(in) :: loop_count, iqm, jqm, ni, nj
      integer, intent(inout) :: ki
      _REAL_, intent(in) :: xi(3), xj(3)
      _REAL_, intent(out) :: W(10*10)
      _REAL_, intent(out) ::RI(22), core(10,2)

!Local
      _REAL_ :: X(3),Y(3),Z(3)
      _REAL_ :: temp_real, temp_real2, anam1, C1, oneRIJ, RIJ
      _REAL_ :: rr,rr2, exp1i, exp1j, sqrtaee, BDD1i, BDD1j, bdd1ij
      _REAL_ :: a, xx11, xx21, xx22, xx31, xx32, xx33, yy11, yy21, yy22
      _REAL_ :: zz11, zz21, zz22, zz31, zz32, zz33, yyzz11, yyzz21, yyzz22
      _REAL_ :: xy11, xy21, xy22, xy31, xy32, xz11, xz21, xz22, xz31, xz32, xz33
      _REAL_ :: YZ11, yz21, yz22, yz31, yz32, css1, css2, csp1, cpps1, cppp1, csp2, cpps2, cppp2, scale
      _REAL_ :: PDDG_EXP1, PDDG_EXP2, PDDG_EXP3, PDDG_EXP4, PDDG_CORR
      integer :: i, j,k
      integer :: i_dimension, j_dimension
      logical :: I_IsSAtom, J_IsSAtom, I_IsSPAtom, J_IsSPAtom, I_IsSPDAtom, J_IsSPDAtom 

      X(1:3) = XI(1:3)-XJ(1:3)
      RIJ = X(1)*X(1) + X(2)*X(2) + X(3)*X(3) 
      rr2 = RIJ*A2_TO_BOHRS2
      oneRIJ = one/SQRT(RIJ)
      RIJ = RIJ*oneRIJ !one/oneRIJ
      rr = RIJ*A_TO_BOHRS
      BDD1i = qm2_params%multip_2c_elec_params(3,iqm) 
      BDD1j = qm2_params%multip_2c_elec_params(3,jqm) 
      BDD1ij = bdd1i+bdd1j
      BDD1ij = bdd1ij*bdd1ij
      SQRTAEE = 1.0d0/sqrt(RR2+bdd1ij)
      X(1:3) = X(1:3)*oneRIJ

      call qm2_repp(qmmm_struct, iqm,jqm,rr,rr2,RI,core,SQRTAEE)
      if(qmmm_nml%qmqm_erep_incore .and. loop_count>0 ) then
! Ross Walker - We will need these repulsion integrals later to
!               calculate qm-qm analytical derivatives so store them
!               in our array.
        
        qm2_struct%qm_qm_e_repul(1:22,loop_count) = RI(1:22)
      end if

      if (ABS(X(3)).GT.0.99999999D0) then
         X(3) = SIGN(1.D0,X(3))
         Y(1) = 0.D0
         Y(2) = 1.D0
         Y(3) = 0.D0
         Z(1) = 1.D0
         Z(2) = 0.D0
         Z(3) = 0.D0
      else
         Z(3)=SQRT(1.D0-X(3)*X(3))
         A=1.D0/Z(3)
         Y(1)=-A*X(2)*SIGN(1.D0,X(1))
         Y(2)=ABS(A*X(1))
         Y(3)=0.D0
         Z(1)=-A*X(1)*X(3)
         Z(2)=-A*X(2)*X(3)
      endif

      I_IsSAtom = qm2_params%natomic_orbs(iqm) == 1
      I_IsSPAtom = qm2_params%natomic_orbs(iqm) == 4
      I_IsSPDAtom = qm2_params%natomic_orbs(iqm) == 9
 
      J_IsSAtom = qm2_params%natomic_orbs(jqm) == 1
      J_IsSPAtom = qm2_params%natomic_orbs(jqm) == 4
      J_IsSPDAtom = qm2_params%natomic_orbs(jqm) == 9

      if ( (.not.I_IsSAtom) .OR. (.not.J_IsSAtom)) then
         XX11 = X(1)*X(1)
         XX21 = X(2)*X(1)
         XX22 = X(2)*X(2)
         XX31 = X(3)*X(1)
         XX32 = X(3)*X(2)
         XX33 = X(3)*X(3)
         YY11 = Y(1)*Y(1)
         YY21 = Y(2)*Y(1)
         YY22 = Y(2)*Y(2)
         ZZ11 = Z(1)*Z(1)
         ZZ21 = Z(2)*Z(1)
         ZZ22 = Z(2)*Z(2)
         ZZ31 = Z(3)*Z(1)
         ZZ32 = Z(3)*Z(2)
         ZZ33 = Z(3)*Z(3)
         YYZZ11 = YY11+ZZ11
         YYZZ21 = YY21+ZZ21
         YYZZ22 = YY22+ZZ22
         XY11 = 2.D0*X(1)*Y(1)
         XY21 =      X(1)*Y(2)+X(2)*Y(1)
         XY22 = 2.D0*X(2)*Y(2)
         XY31 =      X(3)*Y(1)
         XY32 =      X(3)*Y(2)
         XZ11 = 2.D0*X(1)*Z(1)
         XZ21 =      X(1)*Z(2)+X(2)*Z(1)
         XZ22 = 2.D0*X(2)*Z(2)
         XZ31 =      X(1)*Z(3)+X(3)*Z(1)
         XZ32 =      X(2)*Z(3)+X(3)*Z(2)
         XZ33 = 2.D0*X(3)*Z(3)
         YZ11 = 2.D0*Y(1)*Z(1)
         YZ21 =      Y(1)*Z(2)+Y(2)*Z(1)
         YZ22 = 2.D0*Y(2)*Z(2)
         YZ31 =      Y(1)*Z(3)
         YZ32 =      Y(2)*Z(3)
      endif
!     (S S/S S)
      W(1)=RI(1)
      KI = 1
      if (.not.J_IsSAtom) then
!     (S S/PX S)
         W(2)=RI(5)*X(1)
!     (S S/PX PX)
         W(3)=RI(11)*XX11+RI(12)*YYZZ11
!     (S S/PY S)
         W(4)=RI(5)*X(2)
!     (S S/PY PX)
         W(5)=RI(11)*XX21+RI(12)*YYZZ21
!     (S S/PY PY)
         W(6)=RI(11)*XX22+RI(12)*YYZZ22
!     (S S/PZ S)
         W(7)=RI(5)*X(3)
!     (S S/PZ PX)
         W(8)=RI(11)*XX31+RI(12)*ZZ31
!     (S S/PZ PY)
         W(9)=RI(11)*XX32+RI(12)*ZZ32
!     (S S/PZ PZ)
         W(10)=RI(11)*XX33+RI(12)*ZZ33
         KI = 10
      endif
      if (.not.I_IsSAtom) then
         if (.not.J_IsSAtom) then
!     (PX S/S S)
            W(11)=RI(2)*X(1)
!     (PX S/PX S)
            W(12)=RI(6)*XX11+RI(7)*YYZZ11
!     (PX S/PX PX)
            W(13)=X(1)*(RI(13)*XX11+RI(14)*YYZZ11)+RI(15)*(Y(1)*XY11+Z(1)*XZ11)
!     (PX S/PY S)
            W(14)=RI(6)*XX21+RI(7)*YYZZ21
!     (PX S/PY PX)
            W(15)=X(1)*(RI(13)*XX21+RI(14)*YYZZ21)+RI(15)*(Y(1)*XY21+Z(1)*XZ21)
!     (PX S/PY PY)
            W(16)=X(1)*(RI(13)*XX22+RI(14)*YYZZ22)+RI(15)*(Y(1)*XY22+Z(1)*XZ22)
!     (PX S/PZ S)
            W(17)=RI(6)*XX31+RI(7)*ZZ31
!     (PX S/PZ PX)
            W(18)=X(1)*(RI(13)*XX31+RI(14)*ZZ31)+RI(15)*(Y(1)*XY31+Z(1)*XZ31)
!     (PX S/PZ PY)
            W(19)=X(1)*(RI(13)*XX32+RI(14)*ZZ32)+RI(15)*(Y(1)*XY32+Z(1)*XZ32)
!     (PX S/PZ PZ)
            W(20)=X(1)*(RI(13)*XX33+RI(14)*ZZ33)+RI(15)*(          Z(1)*XZ33)
!     (PX PX/S S)
            W(21)=RI(3)*XX11+RI(4)*YYZZ11
!     (PX PX/PX S)
            W(22)=X(1)*(RI(8)*XX11+RI(9)*YYZZ11)+RI(10)*(Y(1)*XY11+Z(1)*XZ11)
!     (PX PX/PX PX)
            W(23) =                                             &
           (RI(16)*XX11+RI(17)*YYZZ11)*XX11+RI(18)*XX11*YYZZ11  &
           +RI(19)*(YY11*YY11+ZZ11*ZZ11)                        &
           +RI(20)*(XY11*XY11+XZ11*XZ11)+RI(21)*(YY11*ZZ11+ZZ11*YY11) &
           +RI(22)*YZ11*YZ11
!     (PX PX/PY S)
            W(24)=X(2)*(RI(8)*XX11+RI(9)*YYZZ11)+RI(10)*(Y(2)*XY11+Z(2)*XZ11)
!     (PX PX/PY PX)
            W(25) =                                                   &
           (RI(16)*XX11+RI(17)*YYZZ11)*XX21+RI(18)*XX11*YYZZ21        &
           +RI(19)*(YY11*YY21+ZZ11*ZZ21)+RI(20)*(XY11*XY21+XZ11*XZ21) &
           +RI(21)*(YY11*ZZ21+ZZ11*YY21)+RI(22)*YZ11*YZ21
!     (PX PX/PY PY)
            W(26) =                                                   &
           (RI(16)*XX11+RI(17)*YYZZ11)*XX22+RI(18)*XX11*YYZZ22        &
           +RI(19)*(YY11*YY22+ZZ11*ZZ22)+RI(20)*(XY11*XY22+XZ11*XZ22) &
           +RI(21)*(YY11*ZZ22+ZZ11*YY22)+RI(22)*YZ11*YZ22
!     (PX PX/PZ S)
            W(27)=X(3)*(RI(8)*XX11+RI(9)*YYZZ11)+RI(10)*(         +Z(3)*XZ11)
!     (PX PX/PZ PX)
            W(28) =                                         &
            (RI(16)*XX11+RI(17)*YYZZ11)*XX31                &
           +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ31      &
           +RI(20)*(XY11*XY31+XZ11*XZ31)+RI(22)*YZ11*YZ31
!     (PX PX/PZ PY)
            W(29) =                                         &
            (RI(16)*XX11+RI(17)*YYZZ11)*XX32                &
           +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ32      &
           +RI(20)*(XY11*XY32+XZ11*XZ32)+RI(22)*YZ11*YZ32
!     (PX PX/PZ PZ)
            W(30) =                                         &
            (RI(16)*XX11+RI(17)*YYZZ11)*XX33                &
           +(RI(18)*XX11+RI(19)*ZZ11+RI(21)*YY11)*ZZ33      &
           +RI(20)*XZ11*XZ33
!     (PY S/S S)
            W(31)=RI(2)*X(2)
!     (PY S/PX S)
            W(32)=RI(6)*XX21+RI(7)*YYZZ21
!     (PY S/PX PX)
            W(33)=X(2)*(RI(13)*XX11+RI(14)*YYZZ11)+RI(15)*(Y(2)*XY11+Z(2)*XZ11)
!     (PY S/PY S)
            W(34)=RI(6)*XX22+RI(7)*YYZZ22
!     (PY S/PY PX)
            W(35)=X(2)*(RI(13)*XX21+RI(14)*YYZZ21)+RI(15)*(Y(2)*XY21+Z(2)*XZ21)
!     (PY S/PY PY)
            W(36)=X(2)*(RI(13)*XX22+RI(14)*YYZZ22)+RI(15)*(Y(2)*XY22+Z(2)*XZ22)
!     (PY S/PZ S)
            W(37)=RI(6)*XX32+RI(7)*ZZ32
!     (PY S/PZ PX)
            W(38)=X(2)*(RI(13)*XX31+RI(14)*ZZ31)+RI(15)*(Y(2)*XY31+Z(2)*XZ31)
!     (PY S/PZ PY)
            W(39)=X(2)*(RI(13)*XX32+RI(14)*ZZ32)+RI(15)*(Y(2)*XY32+Z(2)*XZ32)
!     (PY S/PZ PZ)
            W(40)=X(2)*(RI(13)*XX33+RI(14)*ZZ33)+RI(15)*(         +Z(2)*XZ33)
!     (PY PX/S S)
            W(41)=RI(3)*XX21+RI(4)*YYZZ21
!     (PY PX/PX S)
            W(42)=X(1)*(RI(8)*XX21+RI(9)*YYZZ21)+RI(10)*(Y(1)*XY21+Z(1)*XZ21)
!     (PY PX/PX PX)
            W(43) =                                             &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX11+RI(18)*XX21*YYZZ11  &
           +RI(19)*(YY21*YY11+ZZ21*ZZ11)+RI(20)*(XY21*XY11+XZ21*XZ11) &
           +RI(21)*(YY21*ZZ11+ZZ21*YY11)+RI(22)*YZ21*YZ11
!     (PY PX/PY S)
            W(44)=X(2)*(RI(8)*XX21+RI(9)*YYZZ21)+RI(10)*(Y(2)*XY21+Z(2)*XZ21)
!     (PY PX/PY PX)
            W(45) =                                             &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX21+RI(18)*XX21*YYZZ21  &
           +RI(19)*(YY21*YY21+ZZ21*ZZ21)+RI(20)*(XY21*XY21+XZ21*XZ21) &
           +RI(21)*(YY21*ZZ21+ZZ21*YY21)+RI(22)*YZ21*YZ21
!     (PY PX/PY PY)
            W(46) =                                             &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX22+RI(18)*XX21*YYZZ22  &
           +RI(19)*(YY21*YY22+ZZ21*ZZ22)+RI(20)*(XY21*XY22+XZ21*XZ22) &
           +RI(21)*(YY21*ZZ22+ZZ21*YY22)+RI(22)*YZ21*YZ22
!     (PY PX/PZ S)
            W(47)=X(3)*(RI(8)*XX21+RI(9)*YYZZ21)+RI(10)*(         +Z(3)*XZ21)
!      (PY PX/PZ PX)
            W(48) =                                    &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX31            & 
           +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ31 &
           +RI(20)*(XY21*XY31+XZ21*XZ31)+RI(22)*YZ21*YZ31
!      (PY PX/PZ PY)
            W(49) =                                    &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX32            & 
           +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ32 &
           +RI(20)*(XY21*XY32+XZ21*XZ32)+RI(22)*YZ21*YZ32
!      (PY PX/PZ PZ)
            W(50) =                                    &
           (RI(16)*XX21+RI(17)*YYZZ21)*XX33            &
           +(RI(18)*XX21+RI(19)*ZZ21+RI(21)*YY21)*ZZ33 &
           +RI(20)*XZ21*XZ33 
!     (PY PY/S S)
            W(51)=RI(3)*XX22+RI(4)*YYZZ22
!     (PY PY/PX S)
            W(52)=X(1)*(RI(8)*XX22+RI(9)*YYZZ22)+RI(10)*(Y(1)*XY22+Z(1)*XZ22)
!      (PY PY/PX PX)
            W(53) =                                             &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX11+RI(18)*XX22*YYZZ11  &
           +RI(19)*(YY22*YY11+ZZ22*ZZ11)+RI(20)*(XY22*XY11+XZ22*XZ11) &
           +RI(21)*(YY22*ZZ11+ZZ22*YY11)+RI(22)*YZ22*YZ11
!     (PY PY/PY S)
            W(54)=X(2)*(RI(8)*XX22+RI(9)*YYZZ22)+RI(10)*(Y(2)*XY22+Z(2)*XZ22)
!      (PY PY/PY PX)
            W(55) =                                                   &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX21+RI(18)*XX22*YYZZ21        &
           +RI(19)*(YY22*YY21+ZZ22*ZZ21)+RI(20)*(XY22*XY21+XZ22*XZ21) &
           +RI(21)*(YY22*ZZ21+ZZ22*YY21)+RI(22)*YZ22*YZ21
!      (PY PY/PY PY)
            W(56) =                                                   &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX22+RI(18)*XX22*YYZZ22        &
           +RI(19)*(YY22*YY22+ZZ22*ZZ22)+RI(20)*(XY22*XY22+XZ22*XZ22) &
           +RI(21)*(YY22*ZZ22+ZZ22*YY22)+RI(22)*YZ22*YZ22
!     (PY PY/PZ S)
            W(57)=X(3)*(RI(8)*XX22+RI(9)*YYZZ22)+RI(10)*(         +Z(3)*XZ22)
!      (PY PY/PZ PX)
            W(58) =                                                   &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX31                           &
           +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ31                &
           +RI(20)*(XY22*XY31+XZ22*XZ31)+RI(22)*YZ22*YZ31
!      (PY PY/PZ PY)
            W(59) =                                                   &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX32                           &
           +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ32                &
           +RI(20)*(XY22*XY32+XZ22*XZ32)+RI(22)*YZ22*YZ32
!      (PY PY/PZ PZ)
            W(60) =                                                   &
           (RI(16)*XX22+RI(17)*YYZZ22)*XX33                           &
           +(RI(18)*XX22+RI(19)*ZZ22+RI(21)*YY22)*ZZ33                &
           +RI(20)*XZ22*XZ33
!     (PZ S/SS)
            W(61)=RI(2)*X(3)
!     (PZ S/PX S)
            W(62)=RI(6)*XX31+RI(7)*ZZ31
!     (PZ S/PX PX)
            W(63)=X(3)*(RI(13)*XX11+RI(14)*YYZZ11)+RI(15)*(         +Z(3)*XZ11) 
!     (PZ S/PY S)
            W(64)=RI(6)*XX32+RI(7)*ZZ32
!     (PZ S/PY PX)
            W(65)=X(3)*(RI(13)*XX21+RI(14)*YYZZ21)+RI(15)*(         +Z(3)*XZ21)
!     (PZ S/PY PY)
            W(66)=X(3)*(RI(13)*XX22+RI(14)*YYZZ22)+RI(15)*(         +Z(3)*XZ22) 
!     (PZ S/PZ S)
            W(67)=RI(6)*XX33+RI(7)*ZZ33
!     (PZ S/PZ PX)
            W(68)=X(3)*(RI(13)*XX31+RI(14)*ZZ31)+RI(15)*(         +Z(3)*XZ31)
!     (PZ S/PZ PY)
            W(69)=X(3)*(RI(13)*XX32+RI(14)*ZZ32)+RI(15)*(         +Z(3)*XZ32)
!     (PZ S/PZ PZ)
            W(70)=X(3)*(RI(13)*XX33+RI(14)*ZZ33)+RI(15)*(         +Z(3)*XZ33)
!     (PZ PX/S S)
            W(71)=RI(3)*XX31+RI(4)*ZZ31
!     (PZ PX/PX S)
            W(72)=X(1)*(RI(8)*XX31+RI(9)*ZZ31)+RI(10)*(Y(1)*XY31+Z(1)*XZ31)
!      (PZ PX/PX PX)
            W(73) =                                                     &
           (RI(16)*XX31+RI(17)*ZZ31)*XX11+RI(18)*XX31*YYZZ11            &
           +RI(19)*ZZ31*ZZ11+RI(20)*(XY31*XY11+XZ31*XZ11)               &
           +RI(21)*ZZ31*YY11+RI(22)*YZ31*YZ11
!     (PZ PX/PY S)                                                              
            W(74)=X(2)*(RI(8)*XX31+RI(9)*ZZ31)+RI(10)*(Y(2)*XY31+Z(2)*XZ31)                                  
!      (PZ PX/PY PX)                                                            
            W(75) =                                                     &
           (RI(16)*XX31+RI(17)*ZZ31)*XX21+RI(18)*XX31*YYZZ21            &
           +RI(19)*ZZ31*ZZ21+RI(20)*(XY31*XY21+XZ31*XZ21)               &
           +RI(21)*ZZ31*YY21+RI(22)*YZ31*YZ21
!      (PZ PX/PY PY)                                                            
            W(76) =                                                     &
           (RI(16)*XX31+RI(17)*ZZ31)*XX22+RI(18)*XX31*YYZZ22            &
           +RI(19)*ZZ31*ZZ22+RI(20)*(XY31*XY22+XZ31*XZ22)               &
           +RI(21)*ZZ31*YY22+RI(22)*YZ31*YZ22
!     (PZ PX/PZ S)                                                              
            W(77)=X(3)*(RI(8)*XX31+RI(9)*ZZ31)+RI(10)*(         +Z(3)*XZ31)                                  
!     (PZ PX/PZ PX)                                                     
            W(78) =                                                     &
            (RI(16)*XX31+RI(17)*ZZ31)*XX31                              &        
           +(RI(18)*XX31+RI(19)*ZZ31)*ZZ31                              &        
           +RI(20)*(XY31*XY31+XZ31*XZ31)                                &        
           +RI(22)*YZ31*YZ31                                                    
!      (PZ PX/PZ PY)                                                            
            W(79) =                                                     &
           (RI(16)*XX31+RI(17)*ZZ31)*XX32                               &
           +(RI(18)*XX31+RI(19)*ZZ31)*ZZ32                              &
           +RI(20)*(XY31*XY32+XZ31*XZ32)                                &
           +RI(22)*YZ31*YZ32
!      (PZ PX/PZ PZ)
            W(80) =                                                     &
            (RI(16)*XX31+RI(17)*ZZ31)*XX33                              &
           +(RI(18)*XX31+RI(19)*ZZ31)*ZZ33                              &
           +RI(20)*XZ31*XZ33
!     (PZ PY/S S)
            W(81)=RI(3)*XX32+RI(4)*ZZ32
!     (PZ PY/PX S)                                                              
            W(82)=X(1)*(RI(8)*XX32+RI(9)*ZZ32)+RI(10)*(Y(1)*XY32+Z(1)*XZ32)                                 
!      (PZ PY/PX PX)                                                            
            W(83) =                                                      &
           (RI(16)*XX32+RI(17)*ZZ32)*XX11+RI(18)*XX32*YYZZ11             &
           +RI(19)*ZZ32*ZZ11+RI(20)*(XY32*XY11+XZ32*XZ11)                &
           +RI(21)*ZZ32*YY11+RI(22)*YZ32*YZ11                                                    
!     (PZ PY/PY S)                                                              
            W(84)=X(2)*(RI(8)*XX32+RI(9)*ZZ32)+RI(10)*(Y(2)*XY32+Z(2)*XZ32)                                  
!      (PZ PY/PY PX)
            W(85) =                                                      &
           (RI(16)*XX32+RI(17)*ZZ32)*XX21+RI(18)*XX32*YYZZ21             &
           +RI(19)*ZZ32*ZZ21+RI(20)*(XY32*XY21+XZ32*XZ21)                &
           +RI(21)*ZZ32*YY21+RI(22)*YZ32*YZ21
!      (PZ PY/PY PY)
            W(86) =                                                      &
           (RI(16)*XX32+RI(17)*ZZ32)*XX22+RI(18)*XX32*YYZZ22             &
           +RI(19)*ZZ32*ZZ22+RI(20)*(XY32*XY22+XZ32*XZ22)                &
           +RI(21)*ZZ32*YY22+RI(22)*YZ32*YZ22                                                    
!     (PZ PY/PZ S)                                                              
            W(87)=X(3)*(RI(8)*XX32+RI(9)*ZZ32)+RI(10)*(         +Z(3)*XZ32)                                  
!      (PZ PY/PZ PX)                                                            
            W(88) =                                                       &
            (RI(16)*XX32+RI(17)*ZZ32)*XX31+(RI(18)*XX32+RI(19)*ZZ32)*ZZ31 &
           +RI(20)*(XY32*XY31+XZ32*XZ31)+RI(22)*YZ32*YZ31                                                    
!      (PZ PY/PZ PY)
            W(89) =                                                       &
            (RI(16)*XX32+RI(17)*ZZ32)*XX32+(RI(18)*XX32+RI(19)*ZZ32)*ZZ32 &
           +RI(20)*(XY32*XY32+XZ32*XZ32)+RI(22)*YZ32*YZ32                                                    
!       (PZ PY/PZ PZ)
            W(90) =                                                       &
            (RI(16)*XX32+RI(17)*ZZ32)*XX33+(RI(18)*XX32+RI(19)*ZZ32)*ZZ33 &
           +RI(20)*XZ32*XZ33                                                    
!     (PZ PZ/S S)
            W(91)=RI(3)*XX33+RI(4)*ZZ33
!     (PZ PZ/PX S)
            W(92)=X(1)*(RI(8)*XX33+RI(9)*ZZ33)+RI(10)*(          Z(1)*XZ33)
!       (PZ PZ/PX PX)
            W(93) =                                                       &
           (RI(16)*XX33+RI(17)*ZZ33)*XX11+RI(18)*XX33*YYZZ11              &
           +RI(19)*ZZ33*ZZ11+RI(20)*XZ33*XZ11                             &
           +RI(21)*ZZ33*YY11                                                    
!     (PZ PZ/PY S)                                                              
            W(94)=X(2)*(RI(8)*XX33+RI(9)*ZZ33)+RI(10)*(         +Z(2)*XZ33)                                  
!       (PZ PZ/PY PX)                                                           
            W(95) =                                                       &
           (RI(16)*XX33+RI(17)*ZZ33)*XX21+RI(18)*XX33*YYZZ21              &
           +RI(19)*ZZ33*ZZ21+RI(20)*XZ33*XZ21                             &
           +RI(21)*ZZ33*YY21                                                    
!       (PZ PZ/PY PY)
            W(96) =                                                       &
           (RI(16)*XX33+RI(17)*ZZ33)*XX22+RI(18)*XX33*YYZZ22              &
           +RI(19)*ZZ33*ZZ22+RI(20)*XZ33*XZ22+RI(21)*ZZ33*YY22                                                    
!     (PZ PZ/PZ S)
            W(97)=X(3)*(RI(8)*XX33+RI(9)*ZZ33)+RI(10)*(         +Z(3)*XZ33)                                  
!       (PZ PZ/PZ PX)
            W(98) =                                                       &
            (RI(16)*XX33+RI(17)*ZZ33)*XX31+(RI(18)*XX33+RI(19)*ZZ33)*ZZ31 &
            +RI(20)*XZ33*XZ31                                                    
!       (PZ PZ/PZ PY)
            W(99) =                                                       &
            (RI(16)*XX33+RI(17)*ZZ33)*XX32+(RI(18)*XX33+RI(19)*ZZ33)*ZZ32 &
           +RI(20)*XZ33*XZ32                                                    
!       (PZ PZ/PZ PZ)
            W(100) =                                                      &
            (RI(16)*XX33+RI(17)*ZZ33)*XX33+(RI(18)*XX33+RI(19)*ZZ33)*ZZ33 &
           +RI(20)*XZ33*XZ33
            KI = 100                                                         
         else
!     (PX S/S S)
            W(2)=RI(2)*X(1)
!     (PX PX/S S)
            W(3)=RI(3)*XX11+RI(4)*YYZZ11
!     (PY S/S S)
            W(4)=RI(2)*X(2)
!     (PY PX/S S)
            W(5)=RI(3)*XX21+RI(4)*YYZZ21
!     (PY PY/S S)
            W(6)=RI(3)*XX22+RI(4)*YYZZ22
!     (PZ S/SS)
            W(7)=RI(2)*X(3)
!     (PZ PX/S S)
            W(8)=RI(3)*XX31+RI(4)*ZZ31
!     (PZ PY/S S)
            W(9)=RI(3)*XX32+RI(4)*ZZ32
!     (PZ PZ/S S)
            W(10)=RI(3)*XX33+RI(4)*ZZ33
            KI = 10
         end if
      end if
  
end subroutine qm2_rotate_qmqm

