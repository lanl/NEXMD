      subroutine Vxi (xi,eta)
      
      implicit none
      include 'parH.par'
      include 'int.var'
      include 'parvar.var'
      include 'parvar.cmn'

      real*8 xi(Nb,Nb),eta(Nb,Nb)
      real*8 xis(Lt_M),etas(Lt_M)


! symmetric part:
      l = 0
      do i = 1,Nb
         do j = 1,i
            l = l + 1
            xis(l) = .5d0*(xi(i,j)+xi(j,i))
         enddo
      enddo
      call Vxi_pack (xis,etas)
      l = 0
      do i = 1,Nb
         do j = 1,i-1
            l = l + 1
            eta(i,j) = etas(l)
            eta(j,i) = etas(l)
         enddo
         l = l + 1
         eta(i,i) = etas(l)
      enddo

! antisymmetric part:
      l = 0
      do i = 1,Nb
         do j = 1,i
            l = l + 1
            xis(l) = .5d0*(xi(i,j)-xi(j,i))
         enddo
      enddo
      call Vxi_packA (xis,etas)
      l = 0
      do i = 1,Nb
         do j = 1,i-1
            l = l + 1
            eta(i,j) = eta(i,j) + etas(l)
            eta(j,i) = eta(j,i) - etas(l)
         enddo
         l = l + 1
!         if (abs(etas(l)).gt.1d-10) stop 'Vxi: design'         
      enddo
      return


      entry Vxi_symm (xi,eta)   ! eta = Vxi

! pack:
      l = 0
      do i = 1,Nb
         do j = 1,i
            l = l + 1
            xis(l) = xi(j,i)
         enddo
      enddo

! multiply:
      call Vxi_pack (xis,etas)

! unpack:
      l = 0
      do i = 1,Nb
         do j = 1,i
            l = l + 1
            eta(i,j) = etas(l)
            eta(j,i) = etas(l)
         enddo
      enddo
      end


      subroutine Vxi_pack (xi,eta)
      
      implicit none

      include 'parH.par'
      include 'int.var'
      include 'coulombM.var'
      include 'parvar.var'
      include 'cart.var'
      include 'coulombM.cmn'
      include 'parvar.cmn'
      include 'cart.cmn'
      include 'derivfl.cmn'
	
      real*8 xi(Lt),eta(Lt),fv
      integer im,in,kn
	
	character keywr*6
	common /keywr/ keywr
	logical first
      data first /.true./
      save first

! --- if this is the first time in this routine, load coulomb matrix
      if (first) then
	if (index(keywr,'INDO').NE.0.AND.index(keywr,'MINDO').EQ.0) then
	  print *, keywr,' hamiltonian requested'
	  print *, 'Use *Z program for ', keywr
	  stop
	endif	 
!      call wrb_arr2(Nc,WJ,WK,GSS,GSP,GPP,GP2,
!     +        HSP,GSD,GPD,GDD,'vvv.b','r','s')
         first=.false.
      endif
 

      do i=1,Lt
      eta(i)=0.0
      enddo

      CALL FOCK2(Nb_M,Nat,eta,xi,W,WJ,WK,NFIRST, &
           NMIDLE,NLAST,ID,ITYPE )
	
	if (iderivfl.eq.0) then ! We are not in analytic derivatives     
      CALL FOCK1(Nat,atm,eta,xi,NFIRST,NMIDLE,NLAST, &
       GSS,GSP,GPP,GP2,HSP,GSD,GPD,GDD)
      endif
	
      do i=1,Lt
      eta(i)=eta(i)*2.0
      enddo

	return
	end

      subroutine Vxi_packA (xi,eta)
      
      implicit none

      include 'parH.par'
      include 'int.var'
      include 'coulombM.var'
      include 'parvar.var'
      include 'cart.var'
      include 'coulombM.cmn'
      include 'parvar.cmn'
      include 'cart.cmn'
      include 'derivfl.cmn'
	
      real*8 xi(Lt),eta(Lt),fv
      integer im,in,kn
	
	character keywr*6
	common /keywr/ keywr
	logical first
      data first /.true./
      save first

! --- if this is the first time in this routine, load coulomb matrix
      if (first) then
	if (index(keywr,'INDO').NE.0.AND.index(keywr,'MINDO').EQ.0) then
	  print *, keywr,' hamiltonian requested'
	  print *, 'Use *Z program for ', keywr
	  stop
	endif	 
!      call wrb_arr2(Nc,WJ,WK,GSS,GSP,GPP,GP2,
!     +        HSP,GSD,GPD,GDD,'vvv.b','r','s')
         first=.false.
      endif
 

      do i=1,Lt
      eta(i)=0.0
      enddo

      CALL FOCK2(Nb_M,Nat,eta,xi,W,WJ,WK,NFIRST, &
           NMIDLE,NLAST,ID,ITYPE )     

	if (iderivfl.eq.0) then ! We are not in analytic derivatives
      CALL FOCK11(Nat,atm,eta,xi,NFIRST,NMIDLE,NLAST, &
       GSS,GSP,GPP,GP2,HSP,GSD,GPD,GDD)
      endif

      do i=1,Lt
      eta(i)=eta(i)*2.0
      enddo

	return
	end

      SUBROUTINE FOCK1(NUMAT,atm,F,PTOT,NFIRST,NMIDLE,NLAST, &
       GSS,GSP,GPP,GP2,HSP,GSD,GPD,GDD)

! *********************************************************************
!
! *** COMPUTE THE REMAINING CONTRIBUTIONS TO THE ONE-CENTRE ELEMENTS.
!
! *********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(*),PTOT(*),NFIRST(*),NMIDLE(*),NLAST(*), &
       GSS(*),GSP(*),GPP(*),GP2(*),HSP(*),GSD(*),GPD(*),GDD(*)
       integer atm(*)
       
      DO 100 II=1,NUMAT
         IA=NFIRST(II)
         IB=NMIDLE(II)
         IC=NLAST(II)
	 NI=atm(II)
         DTPOP=0.D0
         DAPOP=0.D0
         PTPOP=0.D0
         PAPOP=0.D0
         GOTO (100,40,30,30,30,20,20,20,20,20)IC-IA+2
   20    DTPOP=PTOT((IC*(IC+1))/2)+PTOT(((IC-1)*(IC))/2) &
              +PTOT(((IC-2)*(IC-1))/2)+PTOT(((IC-3)*(IC-2))/2) &
              +PTOT(((IC-4)*(IC-3))/2)
         DAPOP=0.5*PTOT((IC*(IC+1))/2)+0.5*PTOT(((IC-1)*(IC))/2) &
              +0.5*PTOT(((IC-2)*(IC-1))/2)+0.5*PTOT(((IC-3)*(IC-2))/2) &
              +0.5*PTOT(((IC-4)*(IC-3))/2)
   30    PTPOP=PTOT((IB*(IB+1))/2)+PTOT(((IB-1)*(IB))/2) &
              +PTOT(((IB-2)*(IB-1))/2)
         PAPOP=0.5*PTOT((IB*(IB+1))/2)+0.5*PTOT(((IB-1)*(IB))/2) &
              +0.5*PTOT(((IB-2)*(IB-1))/2)
   40    CONTINUE

!
!     F(S,S)
!
         KA=(IA*(IA+1))/2
         F(KA)=F(KA)+0.5*PTOT(KA)*GSS(NI)+PTPOP*GSP(NI) &
               -PAPOP*HSP(NI) + DTPOP*GSD(NI)
         IF (NI.LT.3) GO TO 100
         IPLUS=IA+1
         L=KA
         DO 70 J=IPLUS,IB
            M=L+IA
            L=L+J
!
!     F(P,P)
!
            F(L)=F(L)+PTOT(KA)*GSP(NI)-0.5*PTOT(KA)*HSP(NI)+ & 
             .5*PTOT(L)*GPP(NI)+(PTPOP-PTOT(L))*GP2(NI) &
            -0.5D0*(PAPOP-0.5*PTOT(L))*(GPP(NI)-GP2(NI)) &
            +DTPOP*GPD(NI)
!
!     F(S,P)
!
   70    F(M)=F(M)+2.D0*PTOT(M)*HSP(NI)-0.5*PTOT(M)*(HSP(NI)+GSP(NI))
!
!     F(P,P*)
!
         IMINUS=IB-1
         DO 80 J=IPLUS,IMINUS
            ICC=J+1
            DO 80 L=ICC,IB
               M=(L*(L-1))/2+J
   80    F(M)=F(M)+PTOT(M)*(GPP(NI)-GP2(NI)) &
            -0.5D0*0.5*PTOT(M)*(GPP(NI)+GP2(NI))
         DO 90 J=IB+1,IC
            M=(J*(J+1))/2
   90    F(M)=F(M)+PTOT(KA)*GSD(NI) &
               +PTPOP*GPD(NI) &
               +(DTPOP-0.5*PTOT(M))*GDD(NI)
  100 CONTINUE
      RETURN
      END

      SUBROUTINE FOCK11(NUMAT,atm,F,PTOT,NFIRST,NMIDLE,NLAST, &
       GSS,GSP,GPP,GP2,HSP,GSD,GPD,GDD)

! *********************************************************************
!
! *** COMPUTE THE REMAINING CONTRIBUTIONS TO THE ONE-CENTRE ELEMENTS.
!
! *********************************************************************
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION F(*),PTOT(*),NFIRST(*),NMIDLE(*),NLAST(*), &
       GSS(*),GSP(*),GPP(*),GP2(*),HSP(*),GSD(*),GPD(*),GDD(*)
       integer atm(*)
       
      DO 100 II=1,NUMAT
         IA=NFIRST(II)
         IB=NMIDLE(II)
         IC=NLAST(II)
         NI=atm(II)
         DTPOP=0.D0
         DAPOP=0.D0
         PTPOP=0.D0
         PAPOP=0.D0
         GOTO (100,40,30,30,30,20,20,20,20,20)IC-IA+2
   20    DTPOP=PTOT((IC*(IC+1))/2)+PTOT(((IC-1)*(IC))/2) &
              +PTOT(((IC-2)*(IC-1))/2)+PTOT(((IC-3)*(IC-2))/2) &
              +PTOT(((IC-4)*(IC-3))/2)
         DAPOP=0.5*PTOT((IC*(IC+1))/2)+0.5*PTOT(((IC-1)*(IC))/2) &
              +0.5*PTOT(((IC-2)*(IC-1))/2)+0.5*PTOT(((IC-3)*(IC-2))/2) &
              +0.5*PTOT(((IC-4)*(IC-3))/2)
   30    PTPOP=PTOT((IB*(IB+1))/2)+PTOT(((IB-1)*(IB))/2) &
              +PTOT(((IB-2)*(IB-1))/2)
         PAPOP=0.5*PTOT((IB*(IB+1))/2)+0.5*PTOT(((IB-1)*(IB))/2) &
              +0.5*PTOT(((IB-2)*(IB-1))/2)
   40    continue
!
!     F(S,S)
!
         KA=(IA*(IA+1))/2
! 0 -diagonal for assymetric
         IF (NI.LT.3) GO TO 100
         IPLUS=IA+1
         L=KA
         DO 70 J=IPLUS,IB
            M=L+IA
            L=L+J
!
!     F(P,P)
!
! 0 -diagonal for assymetric
!
!     F(S,P)
!
   70    F(M)=F(M)+0.5*PTOT(M)*(HSP(NI)-GSP(NI))
!
!     F(P,P*)
!
         IMINUS=IB-1
         DO 80 J=IPLUS,IMINUS
            ICC=J+1
            DO 80 L=ICC,IB
               M=(L*(L-1))/2+J
   80    F(M)=F(M)+0.5D0*0.5*PTOT(M)*(GPP(NI)-3D0*GP2(NI))
         DO 90 J=IB+1,IC
            M=(J*(J+1))/2
   90    F(M)=F(M)+PTOT(KA)*GSD(NI) &
               +PTPOP*GPD(NI) &
               +(DTPOP-0.5*PTOT(M))*GDD(NI)
  100 CONTINUE
      RETURN
      END
      
      SUBROUTINE FOCK2(MAXORB, NUMAT, F, PTOT, W, WJ, WK,  &
        NFIRST, NMIDLE, NLAST, ID, ITYPE )
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'parH.par'
      DIMENSION F(*), PTOT(*), WJ(*), WK(*), NFIRST(*), NMIDLE(*), &
                NLAST(*), W(*)

!***********************************************************************
!
! FOCK2 FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
! MATRIX
! ON INPUT  PTOT = TOTAL DENSITY MATRIX.
!           P    = ALPHA OR BETA DENSITY MATRIX.
!           W    = TWO-ELECTRON INTEGRAL MATRIX.
!
!  ON OUTPUT F   = PARTIAL FOCK MATRIX
!***********************************************************************

      integer ID,ITYPE, NUMCAL
      SAVE IFACT,I1FACT, IONE, LID, ICALCN, JINDEX, JJNDEX, KINDEX
      DIMENSION IFACT(Nb_M), &
      I1FACT(Nb_M), JINDEX(256), KINDEX(256), IJPERM(10), LLPERM(10), &
      PK(16), PJA(16), PJB(16), MMPERM(10), &
      PTOT2(Nb_M,16), JJNDEX(256)
      LOGICAL LID
      
      IF(ICALCN.EQ.0) THEN
         ICALCN=ICALCN+1
!
!   SET UP ARRAY OF LOWER HALF TRIANGLE INDICES (PASCAL'S TRIANGLE)
!
         DO 10 I=1,MAXORB
            IFACT(I)=(I*(I-1))/2
   10    I1FACT(I)=IFACT(I)+I
!
!   SET UP GATHER-SCATTER TYPE ARRAYS FOR USE WITH TWO-ELECTRON
!   INTEGRALS.  JINDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM I
!   INTEGRALS.  JJNDEX ARE THE INDICES OF THE J-INTEGRALS FOR ATOM J
!               KINDEX ARE THE INDICES OF THE K-INTEGRALS
!
         M=0
         DO 20 I=1,4
            DO 20 J=1,4
               IJ=MIN(I,J)
               JI=I+J-IJ
               DO 20 K=1,4
                  IK=MIN(I,K)
                  KI=I+K-IK
                  DO 20 L=1,4
                     M=M+1
                     KL=MIN(K,L)
                     LK=K+L-KL
                     JL=MIN(J,L)
                     LJ=J+L-JL
                     KINDEX(M)= IFACT(LJ) +JL + 10*( IFACT(KI) +IK) -10
   20    JINDEX(M)=(IFACT(JI) + IJ)*10 + IFACT(LK) + KL - 10
         L=0
         DO 30 I=1,4
            I1=(I-1)*4
            DO 30 J=1,I
               I1=I1+1
               L=L+1
               IJPERM(L)=I1
               MMPERM(L)=IJPERM(L)-16
               LLPERM(L)=(I1-1)*16
   30    CONTINUE
         L=0
         DO 40 I=1,10
            M=MMPERM(I)
            L=LLPERM(I)
            DO 40 K=1,16
               L=L+1
               M=M+16
   40    JJNDEX(L)=JINDEX(M)
         LID=(ID.EQ.0)
         IONE=1
         IF(ID.NE.0)IONE=0
!
!      END OF INITIALIZATION
!
      ENDIF
      IF(ITYPE.EQ.4) GOTO 200
!
!     START OF MNDO, AM1, OR PM3 OPTION
!
      KK=0
      L=0
      DO 60 I=1,NUMAT
         IA=NFIRST(I)
         IB=NLAST(I)
         M=0
         DO 50 J=IA,IB
            DO 50 K=IA,IB
               M=M+1
               JK=MIN(J,K)
               KJ=K+J-JK
               JK=JK+(KJ*(KJ-1))/2
               PTOT2(I,M)=PTOT(JK)
   50    CONTINUE
   60 CONTINUE
      DO 190 II=1,NUMAT
         IA=NFIRST(II)
         IB=NLAST(II)
!
!  IF NUMAT=2 THEN WE ARE IN A DERIVATIVE OR IN A MOLECULE CALCULATION
!
         IF(NUMAT.NE.2)THEN
            IMINUS=II-IONE
         ELSE
            IMINUS=II-1
         ENDIF
         DO 180 JJ=1,IMINUS
            JA=NFIRST(JJ)
            JB=NLAST(JJ)
            JC=NMIDLE(JJ)
            IF(LID) THEN
               IF(IB-IA.GE.3.AND.JB-JA.GE.3)THEN
!
!                         HEAVY-ATOM  - HEAVY-ATOM
!
!   EXTRACT COULOMB TERMS
!
                  DO 70 I=1,16
                     PJA(I)=PTOT2(II,I)
   70             PJB(I)=PTOT2(JJ,I)
!
!  COULOMB TERMS
!
                  CALL JAB(IA,JA,LLPERM,JINDEX, JJNDEX, PJA,PJB,W(KK+1), &
      F)
!
!  EXCHANGE TERMS
!
!
!  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX
!
                  L=0
                  DO 80 I=IA,IB
                     I1=IFACT(I)+JA
                     DO 80 J=I1,I1+3
                        L=L+1
   80             PK(L)=0.5*PTOT(J)
                  CALL KAB(IA,JA, PK, W(KK+1), KINDEX, F)
                  KK=KK+100
               ELSEIF(IB-IA.GE.3.AND.JA.EQ.JB)THEN
!
!                         LIGHT-ATOM  - HEAVY-ATOM
!
!
!   COULOMB TERMS
!
                  SUMDIA=0.D0
                  SUMOFF=0.D0
                  LL=I1FACT(JA)
                  K=0
                  DO 100 I=0,3
                     J1=IFACT(IA+I)+IA-1
                     DO 90 J=0,I-1
                        K=K+1
                        J1=J1+1
                        F(J1)=F(J1)+PTOT(LL)*W(KK+K)
   90                SUMOFF=SUMOFF+PTOT(J1)*W(KK+K)
                     J1=J1+1
                     K=K+1
                     F(J1)=F(J1)+PTOT(LL)*W(KK+K)
  100             SUMDIA=SUMDIA+PTOT(J1)*W(KK+K)
                  F(LL)=F(LL)+SUMOFF*2.D0+SUMDIA
!
!  EXCHANGE TERMS
!
!
!  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX
!
                  K=0
                  DO 120 I=IA,IB
                     I1=IFACT(I)+JA
                     SUM=0.D0
                     DO 110 J=IA,IB
                        K=K+1
                        J1=IFACT(J)+JA
  110                SUM=SUM+0.5*PTOT(J1)*W(KK+JINDEX(K))
  120             F(I1)=F(I1)-SUM
                  KK=KK+10
               ELSEIF(JB-JA.GE.3.AND.IA.EQ.IB)THEN
!
!                         HEAVY-ATOM - LIGHT-ATOM
!
!
!   COULOMB TERMS
!
                  SUMDIA=0.D0
                  SUMOFF=0.D0
                  LL=I1FACT(IA)
                  K=0
                  DO 140 I=0,3
                     J1=IFACT(JA+I)+JA-1
                     DO 130 J=0,I-1
                        K=K+1
                        J1=J1+1
                        F(J1)=F(J1)+PTOT(LL)*W(KK+K)
  130                SUMOFF=SUMOFF+PTOT(J1)*W(KK+K)
                     J1=J1+1
                     K=K+1
                     F(J1)=F(J1)+PTOT(LL)*W(KK+K)
  140             SUMDIA=SUMDIA+PTOT(J1)*W(KK+K)
                  F(LL)=F(LL)+SUMOFF*2.D0+SUMDIA
!
!  EXCHANGE TERMS
!
!
!  EXTRACT INTERSECTION OF ATOMS II AND JJ IN THE SPIN DENSITY MATRIX
!
                  K=IFACT(IA)+JA
                  J=0
                  DO 160 I=K,K+3
                     SUM=0.D0
                     DO 150 L=K,K+3
                        J=J+1
  150                SUM=SUM+0.5*PTOT(L)*W(KK+JINDEX(J))
  160             F(I)=F(I)-SUM
                  KK=KK+10
               ELSEIF(JB.EQ.JA.AND.IA.EQ.IB)THEN
!
!                         LIGHT-ATOM - LIGHT-ATOM
!
                  I1=I1FACT(IA)
                  J1=I1FACT(JA)
                  IJ=I1+JA-IA
                  F(I1)=F(I1)+PTOT(J1)*W(KK+1)
                  F(J1)=F(J1)+PTOT(I1)*W(KK+1)
                  F(IJ)=F(IJ)-0.5*PTOT(IJ)*W(KK+1)
                  KK=KK+1
               ENDIF
            ELSE
               DO 170 I=IA,IB
                  KA=IFACT(I)
                  DO 170 J=IA,I
                     KB=IFACT(J)
                     IJ=KA+J
                     AA=2.0D00
                     IF (I.EQ.J) AA=1.0D00
                     DO 170 K=JA,JC
                        KC=IFACT(K)
                        IF(I.GE.K) THEN
                           IK=KA+K
                        ELSE
                           IK=0
                        ENDIF
                        IF(J.GE.K) THEN
                           JK=KB+K
                        ELSE
                           JK=0
                        ENDIF
                        DO 170 L=JA,K
                           IF(I.GE.L) THEN
                              IL=KA+L
                           ELSE
                              IL=0
                           ENDIF
                           IF(J.GE.L) THEN
                              JL=KB+L
                           ELSE
                              JL=0
                           ENDIF
                           KL=KC+L
                           BB=2.0D00
                           IF (K.EQ.L) BB=1.0D00
                           KK=KK+1
                           AJ=WJ(KK)
                           AK=WK(KK)
!
!     A  IS THE REPULSION INTEGRAL (I,J/K,L) WHERE ORBITALS I AND J ARE
!     ON ATOM II, AND ORBITALS K AND L ARE ON ATOM JJ.
!     AA AND BB ARE CORRECTION FACTORS SINCE
!     (I,J/K,L)=(J,I/K,L)=(I,J/L,K)=(J,I/L,K)
!     IJ IS THE LOCATION OF THE MATRIX ELEMENTS BETWEEN ATOMIC ORBITALS
!     I AND J.  SIMILARLY FOR IK ETC.
!
! THIS FORMS THE TWO-ELECTRON TWO-CENTER REPULSION PART OF THE FOCK
! MATRIX.  THE CODE HERE IS HARD TO FOLLOW, AND IMPOSSIBLE TO MODIFY!,
! BUT IT WORKS,
                           IF(KL.LE.IJ)THEN
                              IF(I.EQ.K.AND.AA+BB.LT.2.1D0)THEN
                                 BB=BB*0.5D0
                                 AA=AA*0.5D0
                                 F(IJ)=F(IJ)+BB*AJ*PTOT(KL)
                                 F(KL)=F(KL)+AA*AJ*PTOT(IJ)
                              ELSE
                                 F(IJ)=F(IJ)+BB*AJ*PTOT(KL)
                                 F(KL)=F(KL)+AA*AJ*PTOT(IJ)
                                 A=AK*AA*BB*0.25D0
                                 F(IK)=F(IK)-A*0.5*PTOT(JL)
                                 F(IL)=F(IL)-A*0.5*PTOT(JK)
                                 F(JK)=F(JK)-A*0.5*PTOT(IL)
                                 F(JL)=F(JL)-A*0.5*PTOT(IK)
                              ENDIF
                           ENDIF
  170          CONTINUE
            ENDIF
  180    CONTINUE
  190 CONTINUE
      RETURN
!
!                    START OF MINDO/3 OPTION
!
  200 KR=0
      DO 230 II=1,NUMAT
         IA=NFIRST(II)
         IB=NLAST(II)
         IM1=II-IONE
         DO 220 JJ=1,IM1
            KR=KR+1
            IF(LID)THEN
               ELREP=W(KR)
               ELEXC=ELREP
            ELSE
               ELREP=WJ(KR)
               ELEXC=WK(KR)
            ENDIF
            JA=NFIRST(JJ)
            JB=NLAST(JJ)
            DO 210 I=IA,IB
               KA=IFACT(I)
               KK=KA+I
               DO 210 K=JA,JB
                  LL=I1FACT(K)
                  IK=KA+K
                  F(KK)=F(KK)+PTOT(LL)*ELREP
                  F(LL)=F(LL)+PTOT(KK)*ELREP
  210       F(IK)=F(IK)-0.5*PTOT(IK)*ELEXC
  220    CONTINUE
  230 CONTINUE
      RETURN
      END
      SUBROUTINE JAB(IA,JA,LLPERM,JINDEX, JJNDEX,PJA,PJB,W, F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION LLPERM(10), PJA(16), PJB(16), W(*), F(*), &
      JINDEX(256), JJNDEX(256), SUMA(10), SUMB(10)
!
!  FOR VECTOR MACHINES, REMOVE THE ARRAYS  SUMA AND SUMB, UNCOMMENT
!  THE LINES MARKED CVECTOR, AND COMMENT OUT THE SECOND WHOLE PART
!  OF THE SUBROUTINE
!VECTOR                  I=0
!VECTOR                  DO 100 I5=1,4
!VECTOR                  IIA=IA+I5-1
!VECTOR                  IJA=JA+I5-1
!VECTOR                  IOFF=(IIA*(IIA-1))/2+IA-1
!VECTOR                  JOFF=(IJA*(IJA-1))/2+JA-1
!VECTOR                  DO 100 I6=1,I5
!VECTOR                  IOFF=IOFF+1
!VECTOR                  JOFF=JOFF+1
!VECTOR                        I=I+1
!VECTOR                        L=LLPERM(I)
!VECTOR                        SUMA=0
!VECTOR                        SUMB=0
!VECTOR                        DO 90 K=1,16
!VECTOR                           L=L+1
!VECTOR                           SUMB=SUMB+PJA(K)*W(JJNDEX(L))
!VECTOR   90                   SUMA=SUMA+PJB(K)*W(JINDEX(L))
!VECTOR                        F(IOFF)=F(IOFF)+SUMA
!VECTOR  100             F(JOFF)=F(JOFF)+SUMB
      SUMA( 1)= &
      +PJA( 1)*W(  1)+PJA( 2)*W( 11)+PJA( 3)*W( 31)+PJA( 4)*W( 61) &
      +PJA( 5)*W( 11)+PJA( 6)*W( 21)+PJA( 7)*W( 41)+PJA( 8)*W( 71) &
      +PJA( 9)*W( 31)+PJA(10)*W( 41)+PJA(11)*W( 51)+PJA(12)*W( 81) &
      +PJA(13)*W( 61)+PJA(14)*W( 71)+PJA(15)*W( 81)+PJA(16)*W( 91)
      SUMA( 2)= &
      +PJA( 1)*W(  2)+PJA( 2)*W( 12)+PJA( 3)*W( 32)+PJA( 4)*W( 62) &
      +PJA( 5)*W( 12)+PJA( 6)*W( 22)+PJA( 7)*W( 42)+PJA( 8)*W( 72) &
      +PJA( 9)*W( 32)+PJA(10)*W( 42)+PJA(11)*W( 52)+PJA(12)*W( 82) &
      +PJA(13)*W( 62)+PJA(14)*W( 72)+PJA(15)*W( 82)+PJA(16)*W( 92)
      SUMA( 3)= &
      +PJA( 1)*W(  3)+PJA( 2)*W( 13)+PJA( 3)*W( 33)+PJA( 4)*W( 63) &
      +PJA( 5)*W( 13)+PJA( 6)*W( 23)+PJA( 7)*W( 43)+PJA( 8)*W( 73) &
      +PJA( 9)*W( 33)+PJA(10)*W( 43)+PJA(11)*W( 53)+PJA(12)*W( 83) &
      +PJA(13)*W( 63)+PJA(14)*W( 73)+PJA(15)*W( 83)+PJA(16)*W( 93)
      SUMA( 4)= &
      +PJA( 1)*W(  4)+PJA( 2)*W( 14)+PJA( 3)*W( 34)+PJA( 4)*W( 64) &
      +PJA( 5)*W( 14)+PJA( 6)*W( 24)+PJA( 7)*W( 44)+PJA( 8)*W( 74) &
      +PJA( 9)*W( 34)+PJA(10)*W( 44)+PJA(11)*W( 54)+PJA(12)*W( 84) &
      +PJA(13)*W( 64)+PJA(14)*W( 74)+PJA(15)*W( 84)+PJA(16)*W( 94)
      SUMA( 5)= &
      +PJA( 1)*W(  5)+PJA( 2)*W( 15)+PJA( 3)*W( 35)+PJA( 4)*W( 65) &
      +PJA( 5)*W( 15)+PJA( 6)*W( 25)+PJA( 7)*W( 45)+PJA( 8)*W( 75) &
      +PJA( 9)*W( 35)+PJA(10)*W( 45)+PJA(11)*W( 55)+PJA(12)*W( 85) &
      +PJA(13)*W( 65)+PJA(14)*W( 75)+PJA(15)*W( 85)+PJA(16)*W( 95)
      SUMA( 6)= &
      +PJA( 1)*W(  6)+PJA( 2)*W( 16)+PJA( 3)*W( 36)+PJA( 4)*W( 66) &
      +PJA( 5)*W( 16)+PJA( 6)*W( 26)+PJA( 7)*W( 46)+PJA( 8)*W( 76) &
      +PJA( 9)*W( 36)+PJA(10)*W( 46)+PJA(11)*W( 56)+PJA(12)*W( 86) &
      +PJA(13)*W( 66)+PJA(14)*W( 76)+PJA(15)*W( 86)+PJA(16)*W( 96)
      SUMA( 7)= &
      +PJA( 1)*W(  7)+PJA( 2)*W( 17)+PJA( 3)*W( 37)+PJA( 4)*W( 67) &
      +PJA( 5)*W( 17)+PJA( 6)*W( 27)+PJA( 7)*W( 47)+PJA( 8)*W( 77) &
      +PJA( 9)*W( 37)+PJA(10)*W( 47)+PJA(11)*W( 57)+PJA(12)*W( 87) &
      +PJA(13)*W( 67)+PJA(14)*W( 77)+PJA(15)*W( 87)+PJA(16)*W( 97)
      SUMA( 8)= &
      +PJA( 1)*W(  8)+PJA( 2)*W( 18)+PJA( 3)*W( 38)+PJA( 4)*W( 68) &
      +PJA( 5)*W( 18)+PJA( 6)*W( 28)+PJA( 7)*W( 48)+PJA( 8)*W( 78) &
      +PJA( 9)*W( 38)+PJA(10)*W( 48)+PJA(11)*W( 58)+PJA(12)*W( 88) &
      +PJA(13)*W( 68)+PJA(14)*W( 78)+PJA(15)*W( 88)+PJA(16)*W( 98)
      SUMA( 9)= &
      +PJA( 1)*W(  9)+PJA( 2)*W( 19)+PJA( 3)*W( 39)+PJA( 4)*W( 69) &
      +PJA( 5)*W( 19)+PJA( 6)*W( 29)+PJA( 7)*W( 49)+PJA( 8)*W( 79) &
      +PJA( 9)*W( 39)+PJA(10)*W( 49)+PJA(11)*W( 59)+PJA(12)*W( 89) &
      +PJA(13)*W( 69)+PJA(14)*W( 79)+PJA(15)*W( 89)+PJA(16)*W( 99)
      SUMA(10)= &
      +PJA( 1)*W( 10)+PJA( 2)*W( 20)+PJA( 3)*W( 40)+PJA( 4)*W( 70) &
      +PJA( 5)*W( 20)+PJA( 6)*W( 30)+PJA( 7)*W( 50)+PJA( 8)*W( 80) &
      +PJA( 9)*W( 40)+PJA(10)*W( 50)+PJA(11)*W( 60)+PJA(12)*W( 90) &
      +PJA(13)*W( 70)+PJA(14)*W( 80)+PJA(15)*W( 90)+PJA(16)*W(100)
      SUMB( 1)= &
      +PJB( 1)*W(  1)+PJB( 2)*W(  2)+PJB( 3)*W(  4)+PJB( 4)*W(  7) &
      +PJB( 5)*W(  2)+PJB( 6)*W(  3)+PJB( 7)*W(  5)+PJB( 8)*W(  8) &
      +PJB( 9)*W(  4)+PJB(10)*W(  5)+PJB(11)*W(  6)+PJB(12)*W(  9) &
      +PJB(13)*W(  7)+PJB(14)*W(  8)+PJB(15)*W(  9)+PJB(16)*W( 10)
      SUMB( 2)= &
      +PJB( 1)*W( 11)+PJB( 2)*W( 12)+PJB( 3)*W( 14)+PJB( 4)*W( 17) &
      +PJB( 5)*W( 12)+PJB( 6)*W( 13)+PJB( 7)*W( 15)+PJB( 8)*W( 18) &
      +PJB( 9)*W( 14)+PJB(10)*W( 15)+PJB(11)*W( 16)+PJB(12)*W( 19) &
      +PJB(13)*W( 17)+PJB(14)*W( 18)+PJB(15)*W( 19)+PJB(16)*W( 20)
      SUMB( 3)= &
      +PJB( 1)*W( 21)+PJB( 2)*W( 22)+PJB( 3)*W( 24)+PJB( 4)*W( 27) &
      +PJB( 5)*W( 22)+PJB( 6)*W( 23)+PJB( 7)*W( 25)+PJB( 8)*W( 28) &
      +PJB( 9)*W( 24)+PJB(10)*W( 25)+PJB(11)*W( 26)+PJB(12)*W( 29) &
      +PJB(13)*W( 27)+PJB(14)*W( 28)+PJB(15)*W( 29)+PJB(16)*W( 30)
      SUMB( 4)= &
      +PJB( 1)*W( 31)+PJB( 2)*W( 32)+PJB( 3)*W( 34)+PJB( 4)*W( 37) &
      +PJB( 5)*W( 32)+PJB( 6)*W( 33)+PJB( 7)*W( 35)+PJB( 8)*W( 38) &
      +PJB( 9)*W( 34)+PJB(10)*W( 35)+PJB(11)*W( 36)+PJB(12)*W( 39) &
      +PJB(13)*W( 37)+PJB(14)*W( 38)+PJB(15)*W( 39)+PJB(16)*W( 40)
      SUMB( 5)= &
      +PJB( 1)*W( 41)+PJB( 2)*W( 42)+PJB( 3)*W( 44)+PJB( 4)*W( 47) &
      +PJB( 5)*W( 42)+PJB( 6)*W( 43)+PJB( 7)*W( 45)+PJB( 8)*W( 48) &
      +PJB( 9)*W( 44)+PJB(10)*W( 45)+PJB(11)*W( 46)+PJB(12)*W( 49) &
      +PJB(13)*W( 47)+PJB(14)*W( 48)+PJB(15)*W( 49)+PJB(16)*W( 50)
      SUMB( 6)= &
      +PJB( 1)*W( 51)+PJB( 2)*W( 52)+PJB( 3)*W( 54)+PJB( 4)*W( 57) &
      +PJB( 5)*W( 52)+PJB( 6)*W( 53)+PJB( 7)*W( 55)+PJB( 8)*W( 58) &
      +PJB( 9)*W( 54)+PJB(10)*W( 55)+PJB(11)*W( 56)+PJB(12)*W( 59) &
      +PJB(13)*W( 57)+PJB(14)*W( 58)+PJB(15)*W( 59)+PJB(16)*W( 60)
      SUMB( 7)= &
      +PJB( 1)*W( 61)+PJB( 2)*W( 62)+PJB( 3)*W( 64)+PJB( 4)*W( 67) &
      +PJB( 5)*W( 62)+PJB( 6)*W( 63)+PJB( 7)*W( 65)+PJB( 8)*W( 68) &
      +PJB( 9)*W( 64)+PJB(10)*W( 65)+PJB(11)*W( 66)+PJB(12)*W( 69) &
      +PJB(13)*W( 67)+PJB(14)*W( 68)+PJB(15)*W( 69)+PJB(16)*W( 70)
      SUMB( 8)= &
      +PJB( 1)*W( 71)+PJB( 2)*W( 72)+PJB( 3)*W( 74)+PJB( 4)*W( 77) &
      +PJB( 5)*W( 72)+PJB( 6)*W( 73)+PJB( 7)*W( 75)+PJB( 8)*W( 78) &
      +PJB( 9)*W( 74)+PJB(10)*W( 75)+PJB(11)*W( 76)+PJB(12)*W( 79) &
      +PJB(13)*W( 77)+PJB(14)*W( 78)+PJB(15)*W( 79)+PJB(16)*W( 80)
      SUMB( 9)= &
      +PJB( 1)*W( 81)+PJB( 2)*W( 82)+PJB( 3)*W( 84)+PJB( 4)*W( 87) &
      +PJB( 5)*W( 82)+PJB( 6)*W( 83)+PJB( 7)*W( 85)+PJB( 8)*W( 88) &
      +PJB( 9)*W( 84)+PJB(10)*W( 85)+PJB(11)*W( 86)+PJB(12)*W( 89) &
      +PJB(13)*W( 87)+PJB(14)*W( 88)+PJB(15)*W( 89)+PJB(16)*W( 90)
      SUMB(10)= &
      +PJB( 1)*W( 91)+PJB( 2)*W( 92)+PJB( 3)*W( 94)+PJB( 4)*W( 97) &
      +PJB( 5)*W( 92)+PJB( 6)*W( 93)+PJB( 7)*W( 95)+PJB( 8)*W( 98) &
      +PJB( 9)*W( 94)+PJB(10)*W( 95)+PJB(11)*W( 96)+PJB(12)*W( 99) &
      +PJB(13)*W( 97)+PJB(14)*W( 98)+PJB(15)*W( 99)+PJB(16)*W(100)
      I=0
      DO 10 I5=1,4
         IIA=IA+I5-1
         IJA=JA+I5-1
         IOFF=(IIA*(IIA-1))/2+IA-1
         JOFF=(IJA*(IJA-1))/2+JA-1
         DO 10 I6=1,I5
            IOFF=IOFF+1
            JOFF=JOFF+1
            I=I+1
            F(IOFF)=F(IOFF)+SUMB(I)
   10 F(JOFF)=F(JOFF)+SUMA(I)
      RETURN
      END
      SUBROUTINE KAB(IA,JA, PK, W, KINDEX, F)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PK(*), W(*), F(*),KINDEX(256), SUM(16)
!
!  FOR VECTOR MACHINES, REMOVE THE ARRAY SUM, UNCOMMENT THE LINES
!  MARKED CVECTOR, AND COMMENT OUT THE SECOND WHOLE PART OF THE
!  SUBROUTINE
!
!VECTOR                  L=0
!VECTOR                  M=0
!VECTOR                  DO 130 J1=IA,IA+3
!VECTOR                  J=(J1*(J1-1))/2
!VECTOR                  DO 130 J2=JA,JA+3
!VECTOR                  M=M+1
!VECTOR                  J3=J+J2
!VECTOR                     SUM=0
!VECTOR                     DO 120 I=1,16
!VECTOR                        L=L+1
!VECTOR  120                SUM=SUM+PK(I)*W(KINDEX(L))
!VECTOR  130             F(J3)=F(J3)-SUM
      SUM( 1)= &
      +PK( 1)*W(  1)+PK( 2)*W(  2)+PK( 3)*W(  4)+PK( 4)*W(  7) &
      +PK( 5)*W( 11)+PK( 6)*W( 12)+PK( 7)*W( 14)+PK( 8)*W( 17) &
      +PK( 9)*W( 31)+PK(10)*W( 32)+PK(11)*W( 34)+PK(12)*W( 37) &
      +PK(13)*W( 61)+PK(14)*W( 62)+PK(15)*W( 64)+PK(16)*W( 67)
      SUM( 2)= &
      +PK( 1)*W(  2)+PK( 2)*W(  3)+PK( 3)*W(  5)+PK( 4)*W(  8) &
      +PK( 5)*W( 12)+PK( 6)*W( 13)+PK( 7)*W( 15)+PK( 8)*W( 18) &
      +PK( 9)*W( 32)+PK(10)*W( 33)+PK(11)*W( 35)+PK(12)*W( 38) &
      +PK(13)*W( 62)+PK(14)*W( 63)+PK(15)*W( 65)+PK(16)*W( 68)
      SUM( 3)= &
      +PK( 1)*W(  4)+PK( 2)*W(  5)+PK( 3)*W(  6)+PK( 4)*W(  9) &
      +PK( 5)*W( 14)+PK( 6)*W( 15)+PK( 7)*W( 16)+PK( 8)*W( 19) &
      +PK( 9)*W( 34)+PK(10)*W( 35)+PK(11)*W( 36)+PK(12)*W( 39) &
      +PK(13)*W( 64)+PK(14)*W( 65)+PK(15)*W( 66)+PK(16)*W( 69)
      SUM( 4)= &
      +PK( 1)*W(  7)+PK( 2)*W(  8)+PK( 3)*W(  9)+PK( 4)*W( 10) &
      +PK( 5)*W( 17)+PK( 6)*W( 18)+PK( 7)*W( 19)+PK( 8)*W( 20) &
      +PK( 9)*W( 37)+PK(10)*W( 38)+PK(11)*W( 39)+PK(12)*W( 40) &
      +PK(13)*W( 67)+PK(14)*W( 68)+PK(15)*W( 69)+PK(16)*W( 70)
      SUM( 5)= &
      +PK( 1)*W( 11)+PK( 2)*W( 12)+PK( 3)*W( 14)+PK( 4)*W( 17) &
      +PK( 5)*W( 21)+PK( 6)*W( 22)+PK( 7)*W( 24)+PK( 8)*W( 27) &
      +PK( 9)*W( 41)+PK(10)*W( 42)+PK(11)*W( 44)+PK(12)*W( 47) &
      +PK(13)*W( 71)+PK(14)*W( 72)+PK(15)*W( 74)+PK(16)*W( 77)
      SUM( 6)= &
      +PK( 1)*W( 12)+PK( 2)*W( 13)+PK( 3)*W( 15)+PK( 4)*W( 18) &
      +PK( 5)*W( 22)+PK( 6)*W( 23)+PK( 7)*W( 25)+PK( 8)*W( 28) &
      +PK( 9)*W( 42)+PK(10)*W( 43)+PK(11)*W( 45)+PK(12)*W( 48) &
      +PK(13)*W( 72)+PK(14)*W( 73)+PK(15)*W( 75)+PK(16)*W( 78)
      SUM( 7)= &
      +PK( 1)*W( 14)+PK( 2)*W( 15)+PK( 3)*W( 16)+PK( 4)*W( 19) &
      +PK( 5)*W( 24)+PK( 6)*W( 25)+PK( 7)*W( 26)+PK( 8)*W( 29) &
      +PK( 9)*W( 44)+PK(10)*W( 45)+PK(11)*W( 46)+PK(12)*W( 49) &
      +PK(13)*W( 74)+PK(14)*W( 75)+PK(15)*W( 76)+PK(16)*W( 79)
      SUM( 8)= &
      +PK( 1)*W( 17)+PK( 2)*W( 18)+PK( 3)*W( 19)+PK( 4)*W( 20) &
      +PK( 5)*W( 27)+PK( 6)*W( 28)+PK( 7)*W( 29)+PK( 8)*W( 30) &
      +PK( 9)*W( 47)+PK(10)*W( 48)+PK(11)*W( 49)+PK(12)*W( 50) &
      +PK(13)*W( 77)+PK(14)*W( 78)+PK(15)*W( 79)+PK(16)*W( 80)
      SUM( 9)= &
      +PK( 1)*W( 31)+PK( 2)*W( 32)+PK( 3)*W( 34)+PK( 4)*W( 37) &
      +PK( 5)*W( 41)+PK( 6)*W( 42)+PK( 7)*W( 44)+PK( 8)*W( 47) &
      +PK( 9)*W( 51)+PK(10)*W( 52)+PK(11)*W( 54)+PK(12)*W( 57) &
      +PK(13)*W( 81)+PK(14)*W( 82)+PK(15)*W( 84)+PK(16)*W( 87)
      SUM(10)= &
      +PK( 1)*W( 32)+PK( 2)*W( 33)+PK( 3)*W( 35)+PK( 4)*W( 38) &
      +PK( 5)*W( 42)+PK( 6)*W( 43)+PK( 7)*W( 45)+PK( 8)*W( 48) &
      +PK( 9)*W( 52)+PK(10)*W( 53)+PK(11)*W( 55)+PK(12)*W( 58) &
      +PK(13)*W( 82)+PK(14)*W( 83)+PK(15)*W( 85)+PK(16)*W( 88)
      SUM(11)= &
      +PK( 1)*W( 34)+PK( 2)*W( 35)+PK( 3)*W( 36)+PK( 4)*W( 39) &
      +PK( 5)*W( 44)+PK( 6)*W( 45)+PK( 7)*W( 46)+PK( 8)*W( 49) &
      +PK( 9)*W( 54)+PK(10)*W( 55)+PK(11)*W( 56)+PK(12)*W( 59) &
      +PK(13)*W( 84)+PK(14)*W( 85)+PK(15)*W( 86)+PK(16)*W( 89)
      SUM(12)= &
      +PK( 1)*W( 37)+PK( 2)*W( 38)+PK( 3)*W( 39)+PK( 4)*W( 40) &
      +PK( 5)*W( 47)+PK( 6)*W( 48)+PK( 7)*W( 49)+PK( 8)*W( 50) &
      +PK( 9)*W( 57)+PK(10)*W( 58)+PK(11)*W( 59)+PK(12)*W( 60) &
      +PK(13)*W( 87)+PK(14)*W( 88)+PK(15)*W( 89)+PK(16)*W( 90)
      SUM(13)= &
      +PK( 1)*W( 61)+PK( 2)*W( 62)+PK( 3)*W( 64)+PK( 4)*W( 67) &
      +PK( 5)*W( 71)+PK( 6)*W( 72)+PK( 7)*W( 74)+PK( 8)*W( 77) &
      +PK( 9)*W( 81)+PK(10)*W( 82)+PK(11)*W( 84)+PK(12)*W( 87) &
      +PK(13)*W( 91)+PK(14)*W( 92)+PK(15)*W( 94)+PK(16)*W( 97)
      SUM(14)= &
      +PK( 1)*W( 62)+PK( 2)*W( 63)+PK( 3)*W( 65)+PK( 4)*W( 68) &
      +PK( 5)*W( 72)+PK( 6)*W( 73)+PK( 7)*W( 75)+PK( 8)*W( 78) &
      +PK( 9)*W( 82)+PK(10)*W( 83)+PK(11)*W( 85)+PK(12)*W( 88) &
      +PK(13)*W( 92)+PK(14)*W( 93)+PK(15)*W( 95)+PK(16)*W( 98)
      SUM(15)= &
      +PK( 1)*W( 64)+PK( 2)*W( 65)+PK( 3)*W( 66)+PK( 4)*W( 69) &
      +PK( 5)*W( 74)+PK( 6)*W( 75)+PK( 7)*W( 76)+PK( 8)*W( 79) &
      +PK( 9)*W( 84)+PK(10)*W( 85)+PK(11)*W( 86)+PK(12)*W( 89) &
      +PK(13)*W( 94)+PK(14)*W( 95)+PK(15)*W( 96)+PK(16)*W( 99)
      SUM(16)= &
      +PK( 1)*W( 67)+PK( 2)*W( 68)+PK( 3)*W( 69)+PK( 4)*W( 70) &
      +PK( 5)*W( 77)+PK( 6)*W( 78)+PK( 7)*W( 79)+PK( 8)*W( 80) &
      +PK( 9)*W( 87)+PK(10)*W( 88)+PK(11)*W( 89)+PK(12)*W( 90) &
      +PK(13)*W( 97)+PK(14)*W( 98)+PK(15)*W( 99)+PK(16)*W(100)
      IF(IA.GT.JA)THEN
         M=0
         DO 10 J1=IA,IA+3
            J=(J1*(J1-1))/2
            DO 10 J2=JA,JA+3
               M=M+1
               J3=J+J2
   10    F(J3)=F(J3)-SUM(M)
      ELSE
!
!   IA IS LESS THAN JA, THEREFORE USE OTHER HALF OF TRIANGLE
!
         M=0
         DO 20 J1=IA,IA+3
            DO 20 J2=JA,JA+3
               M=M+1
               J3=(J2*(J2-1))/2+J1
   20    F(J3)=F(J3)-SUM(M)
      ENDIF
      RETURN
      END
