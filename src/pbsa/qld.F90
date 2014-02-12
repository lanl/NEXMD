#define _REAL_ double precision

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ql0001 here]
subroutine ql0001(m,me,mmax,n,nmax,mnn,c,d,a,b,xl,xu, &
      x,u,iout,ifail,iprint,war,lwar,iwar,liwar,eps)
   implicit none
   
   !**************************************************************************
   
   
   !             SOLUTION OF QUADRATIC PROGRAMMING PROBLEMS
   
   
   
   !   QL0001 SOLVES THE QUADRATIC PROGRAMMING PROBLEM
   
   !   MINIMIZE        .5*X'*C*X + D'*X
   !   SUBJECT TO      A(J)*X  +  B(J)   =  0  ,  J=1,...,ME
   !                   A(J)*X  +  B(J)  >=  0  ,  J=ME+1,...,M
   !                   XL  <=  X  <=  XU
   
   !   HERE C MUST BE AN N BY N SYMMETRIC AND POSITIVE MATRIX, D AN N-DIMENSIONAL
   !   VECTOR, A AN M BY N MATRIX AND B AN M-DIMENSIONAL VECTOR. THE ABOVE
   !   SITUATION IS INDICATED BY IWAR(1)=1. ALTERNATIVELY, I.E. IF IWAR(1)=0,
   !   THE OBJECTIVE FUNCTION MATRIX CAN ALSO BE PROVIDED IN FACTORIZED FORM.
   !   IN THIS CASE, C IS AN UPPER TRIANGULAR MATRIX.
   
   !   THE SUBROUTINE REORGANIZES SOME DATA SO THAT THE PROBLEM CAN BE SOLVED
   !   BY A MODIFICATION OF AN ALGORITHM PROPOSED BY POWELL (1983).
   
   
   !   USAGE:
   
   !      QL0001(M,ME,MMAX,N,NMAX,MNN,C,D,A,B,XL,XU,X,U,IOUT,IFAIL,IPRINT,
   !             WAR,LWAR,IWAR,LIWAR)
   
   
   !   DEFINITION OF THE PARAMETERS:
   
   !   M :        TOTAL NUMBER OF CONSTRAINTS.
   !   ME :       NUMBER OF EQUALITY CONSTRAINTS.
   !   MMAX :     ROW DIMENSION OF A. MMAX MUST BE AT LEAST ONE AND GREATER
   !              THAN M.
   !   N :        NUMBER OF VARIABLES.
   !   NMAX :     ROW DIMENSION OF C. NMAX MUST BE GREATER OR EQUAL TO N.
   !   MNN :      MUST BE EQUAL TO M + N + N.
   !   C(NMAX,NMAX): OBJECTIVE FUNCTION MATRIX WHICH SHOULD BE SYMMETRIC AND
   !              POSITIVE DEFINITE. IF IWAR(1) = 0, C IS SUPPOSED TO BE THE
   !              CHOLESKEY-FACTOR OF ANOTHER MATRIX, I.E. C IS UPPER
   !              TRIANGULAR.
   !   D(NMAX) :  CONTAINS THE CONSTANT VECTOR OF THE OBJECTIVE FUNCTION.
   !   A(MMAX,NMAX): CONTAINS THE DATA MATRIX OF THE LINEAR CONSTRAINTS.
   !   B(MMAX) :  CONTAINS THE CONSTANT DATA OF THE LINEAR CONSTRAINTS.
   !   XL(N),XU(N): CONTAIN THE LOWER AND UPPER BOUNDS FOR THE VARIABLES.
   !   X(N) :     ON RETURN, X CONTAINS THE OPTIMAL SOLUTION VECTOR.
   !   U(MNN) :   ON RETURN, U CONTAINS THE LAGRANGE MULTIPLIERS. THE FIRST
   !              M POSITIONS ARE RESERVED FOR THE MULTIPLIERS OF THE M
   !              LINEAR CONSTRAINTS AND THE SUBSEQUENT ONES FOR THE
   !              MULTIPLIERS OF THE LOWER AND UPPER BOUNDS. ON SUCCESSFUL
   !              TERMINATION, ALL VALUES OF U WITH RESPECT TO INEQUALITIES
   !              AND BOUNDS SHOULD BE GREATER OR EQUAL TO ZERO.
   !   IOUT :     INTEGER INDICATING THE DESIRED OUTPUT UNIT NUMBER, I.E.
   !              ALL WRITE-STATEMENTS START WITH 'WRITE(IOUT,... '.
   !   IFAIL :    SHOWS THE TERMINATION REASON.
   !      IFAIL = 0 :   SUCCESSFUL RETURN.
   !      IFAIL = 1 :   TOO MANY ITERATIONS (MORE THAN 40*(N+M)).
   !      IFAIL = 2 :   ACCURACY INSUFFICIENT TO SATISFY CONVERGENCE
   !                    CRITERION.
   !      IFAIL = 5 :   LENGTH OF A WORKING ARRAY IS TOO SHORT.
   !      IFAIL > 10 :  THE CONSTRAINTS ARE INCONSISTENT.
   !   IPRINT :   OUTPUT CONTROL.
   !      IPRINT = 0 :  NO OUTPUT OF QL0001.
   !      IPRINT > 0 :  BRIEF OUTPUT IN ERROR CASES.
   !   WAR(LWAR) : _REAL_ WORKING ARRAY. THE LENGTH LWAR SHOULD BE GRATER THAN
   !               3*NMAX*NMAX/2 + 10*NMAX + 2*MMAX.
   !   IWAR(LIWAR): INTEGER WORKING ARRAY. THE LENGTH LIWAR SHOULD BE AT
   !              LEAST N.
   !              IF IWAR(1)=1 INITIALLY, THEN THE CHOLESKY DECOMPOSITION
   !              WHICH IS REQUIRED BY THE DUAL ALGORITHM TO GET THE FIRST
   !              UNCONSTRAINED MINIMUM OF THE OBJECTIVE FUNCTION, IS
   !              PERFORMED INTERNALLY. OTHERWISE, I.E. IF IWAR(1)=0, THEN
   !              IT IS ASSUMED THAT THE USER PROVIDES THE INITIAL FAC-
   !              TORIZATION BY HIMSELF AND STORES IT IN THE UPPER TRIAN-
   !              GULAR PART OF THE ARRAY C.
   
   !   A NAMED COMMON-BLOCK  /CMACHE/EPS   MUST BE PROVIDED BY THE USER,
   !   WHERE EPS DEFINES A GUESS FOR THE UNDERLYING MACHINE PRECISION.
   
   
   !   AUTHOR (C): K. SCHITTKOWSKI,
   !               MATHEMATISCHES INSTITUT,
   !               UNIVERSITAET BAYREUTH,
   !               95440 BAYREUTH,
   !               GERMANY, F.R.
   
   
   !   VERSION:    1.5  (JUNE, 1991)
   
   
   !*********************************************************************
   
   
   !Passed variables  lw ??
   integer m,me,nmax,mmax,n,mnn,lwar,liwar,iout,ifail,iprint,iwar(liwar)
   _REAL_  c(nmax,n),d(n),a(mmax,n),b(mmax), &
           xl(n),xu(n),x(n),u(mnn),war(lwar)
   _REAL_  diag,zero,eps,qpeps,ten
   integer inw1,inw2,j,lw,mn,i,idiag,info,nact,maxit,in
   logical lql
   


   !     INTRINSIC FUNCTIONS:  sqrt
   
   !common /cmache/eps
   
   !     CONSTANT DATA
   


   lql=.false.
   if (iwar(1) == 1) lql=.true.
   zero=0.0d+0
   ten=1.d+1
   maxit=40*(m+n)
   qpeps=eps
   inw1=1
   inw2=inw1+mmax
   
   !     PREPARE PROBLEM DATA FOR EXECUTION
   
   if (m <= 0) goto 20
   in=inw1
   do j=1,m
      war(in)=-b(j)
      in=in+1
   end do
   20 lw=3*nmax*nmax/2 + 10*nmax + m
   if ((inw2+lw) > lwar) goto 80
   if (liwar < n) goto 81
   if (mnn < m+n+n) goto 82
   mn=m+n
   
   !     CALL OF QL0002
   
   call ql0002(n,m,me,mmax,mn,mnn,nmax,lql,a,war(inw1), &
         d,c,xl,xu,x,nact,iwar,maxit,qpeps,info,diag, &
         war(inw2),lw)
   
   !     TEST OF MATRIX CORRECTIONS
   
   ifail=0
   if (info == 1) goto 40
   if (info == 2) goto 90
   idiag=0
   if ((diag > zero).and.(diag < 1000.0)) idiag=diag
   if ((iprint > 0).and.(idiag > 0)) &
         write(iout,1000) idiag
   if (info < 0) goto  70
   
   !     REORDER MULTIPLIER
   
   do j=1,mnn
      u(j)=zero
   end do
   in=inw2-1
   if (nact == 0) goto 30
   do i=1,nact
      j=iwar(i)
      u(j)=war(in+i)
   end do
   30 continue
   return
   
   !     ERROR MESSAGES
   
   70 ifail=-info+10
   if ((iprint > 0).and.(nact > 0)) &
         write(iout,1100) -info,(iwar(i),i=1,nact)
   return
   80 ifail=5
   if (iprint > 0) write(iout,1200)
   return
   81 ifail=5
   if (iprint > 0) write(iout,1210)
   return
   82 ifail=5
   if (iprint > 0) write(iout,1220)
   return
   40 ifail=1
   if (iprint > 0) write(iout,1300) maxit
   return
   90 ifail=2
   if (iprint > 0) write(iout,1400)
   return
   
   !     FORMAT-INSTRUCTIONS
   
   1000 format(/8x,28h***ql: matrix g was enlarged,i3, &
         20h-times by unitmatrix)
   1100 format(/8x,18h***ql: constraint ,i5, &
         19h not consistent to ,/,(10x,10i5))
   1200 format(/8x,21h***ql: lwar too small)
   1210 format(/8x,22h***ql: liwar too small)
   1220 format(/8x,20h***ql: mnn too small)
   1300 format(/8x,37h***ql: too many iterations (more than,i6,1h))
   1400 format(/8x,50h***ql: accuracy insufficient to attain convergence)
end subroutine ql0001 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ql0002 here]
subroutine ql0002(n,m,meq,mmax,mn,mnn,nmax,lql,a,b,grad,g, &
      xl,xu,x,nact,iact,maxit,vsmall,info,diag,w,lw)
   implicit none
   
   !**************************************************************************
   
   
   !   THIS SUBROUTINE SOLVES THE QUADRATIC PROGRAMMING PROBLEM
   
   !       MINIMIZE      GRAD'*X  +  0.5 * X*G*X
   !       SUBJECT TO    A(K)*X  =  B(K)   K=1,2,...,MEQ,
   !                     A(K)*X >=  B(K)   K=MEQ+1,...,M,
   !                     XL  <=  X  <=  XU
   
   !   THE QUADRATIC PROGRAMMING METHOD PROCEEDS FROM AN INITIAL CHOLESKY-
   !   DECOMPOSITION OF THE OBJECTIVE FUNCTION MATRIX, TO CALCULATE THE
   !   UNIQUELY DETERMINED MINIMIZER OF THE UNCONSTRAINED PROBLEM.
   !   SUCCESSIVELY ALL VIOLATED CONSTRAINTS ARE ADDED TO A WORKING SET
   !   AND A MINIMIZER OF THE OBJECTIVE FUNCTION SUBJECT TO ALL CONSTRAINTS
   !   IN THIS WORKING SET IS COMPUTED. IT IS POSSIBLE THAT CONSTRAINTS
   !   HAVE TO LEAVE THE WORKING SET.
   
   
   !   DESCRIPTION OF PARAMETERS:
   
   !     N        : IS THE NUMBER OF VARIABLES.
   !     M        : TOTAL NUMBER OF CONSTRAINTS.
   !     MEQ      : NUMBER OF EQUALITY CONTRAINTS.
   !     MMAX     : ROW DIMENSION OF A, DIMENSION OF B. MMAX MUST BE AT
   !                LEAST ONE AND GREATER OR EQUAL TO M.
   !     MN       : MUST BE EQUAL M + N.
   !     MNN      : MUST BE EQUAL M + N + N.
   !     NMAX     : ROW DIEMSION OF G. MUST BE AT LEAST N.
   !     LQL      : DETERMINES INITIAL DECOMPOSITION.
   !        LQL = .FALSE.  : THE UPPER TRIANGULAR PART OF THE MATRIX G
   !                         CONTAINS INITIALLY THE CHOLESKY-FACTOR OF A SUITABLE
   !                         DECOMPOSITION.
   !        LQL = .TRUE.   : THE INITIAL CHOLESKY-FACTORISATION OF G IS TO BE
   !                         PERFORMED BY THE ALGORITHM.
   !     A(MMAX,NMAX) : A IS A MATRIX WHOSE COLUMNS ARE THE CONSTRAINTS NORMALS.
   !     B(MMAX)  : CONTAINS THE RIGHT HAND SIDES OF THE CONSTRAINTS.
   !     GRAD(N)  : CONTAINS THE OBJECTIVE FUNCTION VECTOR GRAD.
   !     G(NMAX,N): CONTAINS THE SYMMETRIC OBJECTIVE FUNCTION MATRIX.
   !     XL(N), XU(N): CONTAIN THE LOWER AND UPPER BOUNDS FOR X.
   !     X(N)     : VECTOR OF VARIABLES.
   !     NACT     : FINAL NUMBER OF ACTIVE CONSTRAINTS.
   !     IACT(K) (K=1,2,...,NACT): INDICES OF THE FINAL ACTIVE CONSTRAINTS.
   !     INFO     : REASON FOR THE RETURN FROM THE SUBROUTINE.
   !         INFO = 0 : CALCULATION WAS TERMINATED SUCCESSFULLY.
   !         INFO = 1 : MAXIMUM NUMBER OF ITERATIONS ATTAINED.
   !         INFO = 2 : ACCURACY IS INSUFFICIENT TO MAINTAIN INCREASING
   !                    FUNCTION VALUES.
   !         INFO < 0 : THE CONSTRAINT WITH INDEX ABS(INFO) AND THE CON-
   !                    STRAINTS WHOSE INDICES ARE IACT(K), K=1,2,...,NACT,
   !                    ARE INCONSISTENT.
   !     MAXIT    : MAXIMUM NUMBER OF ITERATIONS.
   !     VSMALL   : REQUIRED ACCURACY TO BE ACHIEVED (E.G. IN THE ORDER OF THE
   !                MACHINE PRECISION FOR SMALL AND WELL-CONDITIONED PROBLEMS).
   !     DIAG     : ON RETURN DIAG IS EQUAL TO THE MULTIPLE OF THE UNIT MATRIX
   !                THAT WAS ADDED TO G TO ACHIEVE POSITIVE DEFINITENESS.
   !     W(LW)    : THE ELEMENTS OF W(.) ARE USED FOR WORKING SPACE. THE LENGTH
   !                OF W MUST NOT BE LESS THAN (1.5*NMAX*NMAX + 10*NMAX + M).
   !                WHEN INFO = 0 ON RETURN, THE LAGRANGE MULTIPLIERS OF THE
   !                FINAL ACTIVE CONSTRAINTS ARE HELD IN W(K), K=1,2,...,NACT.
   !   THE VALUES OF N, M, MEQ, MMAX, MN, MNN AND NMAX AND THE ELEMENTS OF
   !   A, B, GRAD AND G ARE NOT ALTERED.
   
   !   THE FOLLOWING INTEGERS ARE USED TO PARTITION W:
   !     THE FIRST N ELEMENTS OF W HOLD LAGRANGE MULTIPLIER ESTIMATES.
   !     W(IWZ+I+(N-1)*J) HOLDS THE MATRIX ELEMENT Z(I,J).
   !     W(IWR+I+0.5*J*(J-1)) HOLDS THE UPPER TRIANGULAR MATRIX
   !       ELEMENT R(I,J). THE SUBSEQUENT N COMPONENTS OF W MAY BE
   !       TREATED AS AN EXTRA COLUMN OF R(.,.).
   !     W(IWW-N+I) (I=1,2,...,N) ARE USED FOR TEMPORARY STORAGE.
   !     W(IWW+I) (I=1,2,...,N) ARE USED FOR TEMPORARY STORAGE.
   !     W(IWD+I) (I=1,2,...,N) HOLDS G(I,I) DURING THE CALCULATION.
   !     W(IWX+I) (I=1,2,...,N) HOLDS VARIABLES THAT WILL BE USED TO
   !       TEST THAT THE ITERATIONS INCREASE THE OBJECTIVE FUNCTION.
   !     W(IWA+K) (K=1,2,...,M) USUALLY HOLDS THE RECIPROCAL OF THE
   !       LENGTH OF THE K-TH CONSTRAINT, BUT ITS SIGN INDICATES
   !       WHETHER THE CONSTRAINT IS ACTIVE.
   
   
   !   AUTHOR:    K. SCHITTKOWSKI,
   !              MATHEMATISCHES INSTITUT,
   !              UNIVERSITAET BAYREUTH,
   !              8580 BAYREUTH,
   !              GERMANY, F.R.
   
   !   AUTHOR OF ORIGINAL VERSION:
   !              M.J.D. POWELL, DAMTP,
   !              UNIVERSITY OF CAMBRIDGE, SILVER STREET
   !              CAMBRIDGE,
   !              ENGLAND
   
   
   !   REFERENCE: M.J.D. POWELL: ZQPCVX, A FORTRAN SUBROUTINE FOR CONVEX
   !              PROGRAMMING, REPORT DAMTP/1983/NA17, UNIVERSITY OF
   !              CAMBRIDGE, ENGLAND, 1983.
   
   
   !   VERSION :  2.0 (MARCH, 1987)
   
   
   !*************************************************************************
   
!subroutine ql0002(n,m,meq,mmax,mn,mnn,nmax,lql,a,b,grad,g, &
!      xl,xu,x,nact,iact,maxit,vsmall,info,diag,w,lw)

   integer  mmax,nmax,n,lw,nflag,iwwn,iact(n)
   _REAL_   a(mmax,n),b(mmax),grad(n),g(nmax,n),x(n), &
            w(lw),xl(n),xu(n),fmax,fmin
   integer m,meq,mn,mnn,nact,info,maxit
   _REAL_ cvmax,diag,diagr,fdiff,fdiffa,ga,gb,parinc,parnew &
         ,ratio,res,step,sum,sumx,sumy,suma,sumb,sumc,temp,tempa, &
         vsmall,xmag,xmagr,zero,one,two,onha,vfact
   
   !   INTRINSIC FUNCTIONS:   DMAX1,sqrt,abs,fmin
   
   integer iwz,iwr,iww,iwd,iwa,ifinc,kfinc,k,i,ia,id,ii,ir,ira, &
         irb,j,nm,iz,iza,iterc,itref,jfinc,iflag,iws,is,k1,iw,kk,ip, &
         ipp,il,iu,ju,kflag,lflag,jflag,kdrop,nu,mflag,knext,ix,iwx, &
         iwy,iy,jl
   logical lql,lower
   
   !   INITIAL ADDRESSES
   
   iwz=nmax
   iwr=iwz+nmax*nmax
   iww=iwr+(nmax*(nmax+3))/2
   iwd=iww+nmax
   iwx=iwd+nmax
   iwa=iwx+nmax
   
   !     SET SOME CONSTANTS.
   
   zero=0.d+0
   one=1.d+0
   two=2.d+0
   onha=1.5d+0
   vfact=1.d+0
   
   !     SET SOME PARAMETERS.
   !     NUMBER LESS THAN VSMALL ARE ASSUMED TO BE NEGLIGIBLE.
   !     THE MULTIPLE OF I THAT IS ADDED TO G IS AT MOST DIAGR TIMES
   !       THE LEAST MULTIPLE OF I THAT GIVES POSITIVE DEFINITENESS.
   !     X IS RE-INITIALISED IF ITS MAGNITUDE IS REDUCED BY THE
   !       FACTOR XMAGR.
   !     A CHECK IS MADE FOR AN INCREASE IN F EVERY IFINC ITERATIONS,
   !       AFTER KFINC ITERATIONS ARE COMPLETED.
   
   diagr=two
   diag=zero
   xmagr=1.0d-2
   ifinc=3
   kfinc=max0(10,n)
   
   !     FIND THE RECIPROCALS OF THE LENGTHS OF THE CONSTRAINT NORMALS.
   !     RETURN IF A CONSTRAINT IS INFEASIBLE DUE TO A ZERO NORMAL.
   
   nact=0
   if (m <= 0) goto 45
   do k=1,m
      sum=zero
      do i=1,n
         sum=sum+a(k,i)**2
      end do
      if (sum > zero) goto 20
      if (b(k) == zero) goto 30
      info=-k
      if (k <= meq) goto 730
      if (b(k)) 30,30,730
      20 sum=one/sqrt(sum)
      30 ia=iwa+k
      w(ia)=sum
   end do
   45 do k=1,n
      ia=iwa+m+k
      w(ia)=one
   end do
   
   !     IF NECESSARY INCREASE THE DIAGONAL ELEMENTS OF G.
   
   if (.not. lql) goto 165
   do i=1,n
      id=iwd+i
      w(id)=g(i,i)
      diag=fmax(diag,vsmall-w(id))
      if (i == n) cycle
      ii=i+1
      do j=ii,n
         ga=-fmin(w(id),g(j,j))
         gb=abs(w(id)-g(j,j))+abs(g(i,j))
         if (gb > zero) ga=ga+g(i,j)**2/gb
         diag=fmax(diag,ga)
      end do
   end do
   if (diag <= zero) goto 90
   70 diag=diagr*diag
   do i=1,n
      id=iwd+i
      g(i,i)=diag+w(id)
   end do
   
   !     FORM THE CHOLESKY FACTORISATION OF G. THE TRANSPOSE
   !     OF THE FACTOR WILL BE PLACED IN THE R-PARTITION OF W.
   
   90 ir=iwr
   do j=1,n
      ira=iwr
      irb=ir+1
      do i=1,j
         temp=g(i,j)
         if (i == 1) goto 110
         do k=irb,ir
            ira=ira+1
            temp=temp-w(k)*w(ira)
         end do
         110 ir=ir+1
         ira=ira+1
         if (i < j) w(ir)=temp/w(ira)
      end do
      if (temp < vsmall) goto 140
      w(ir)=sqrt(temp)
   end do
   goto 170
   
   !     INCREASE FURTHER THE DIAGONAL ELEMENT OF G.
   
   140 w(j)=one
   sumx=one
   k=j
   150 sum=zero
   ira=ir-1
   do i=k,j
      sum=sum-w(ira)*w(i)
      ira=ira+i
   end do
   ir=ir-k
   k=k-1
   w(k)=sum/w(ir)
   sumx=sumx+w(k)**2
   if (k >= 2) goto 150
   diag=diag+vsmall-temp/sumx
   goto 70
   
   !     STORE THE CHOLESKY FACTORISATION IN THE R-PARTITION
   !     OF W.
   
   165 ir=iwr
   do i=1,n
      do j=1,i
         ir=ir+1
         w(ir)=g(j,i)
      end do
   end do
      
   !     SET Z THE INVERSE OF THE MATRIX IN R.
      
   170 nm=n-1
   do i=1,n
      iz=iwz+i
      if (i == 1) goto 190
      do j=2,i
         w(iz)=zero
         iz=iz+n
      end do
      190 ir=iwr+(i+i*i)/2
      w(iz)=one/w(ir)
      if (i == n) cycle
      iza=iz
      do j=i,nm
         ir=ir+i
         sum=zero
         do k=iza,iz,n
            sum=sum+w(k)*w(ir)
            ir=ir+1
         end do
         iz=iz+n
         w(iz)=-sum/w(ir)
      end do
   end do
   
   !     SET THE INITIAL VALUES OF SOME VARIABLES.
   !     ITERC COUNTS THE NUMBER OF ITERATIONS.
   !     ITREF IS SET TO ONE WHEN ITERATIVE REFINEMENT IS REQUIRED.
   !     JFINC INDICATES WHEN TO TEST FOR AN INCREASE IN F.
   
   iterc=1
   itref=0
   jfinc=-kfinc
   
   !     SET X TO ZERO AND SET THE CORRESPONDING RESIDUALS OF THE
   !     KUHN-TUCKER CONDITIONS.
   
   230 iflag=1
   iws=iww-n
   do i=1,n
      x(i)=zero
      iw=iww+i
      w(iw)=grad(i)
      if (i > nact) cycle
      w(i)=zero
      is=iws+i
      k=iact(i)
      if (k <= m) goto 235
      if (k > mn) goto 234
      k1=k-m
      w(is)=xl(k1)
      cycle
      234 k1=k-mn
      w(is)=-xu(k1)
      cycle
      235 w(is)=b(k)
   end do
   xmag=zero
   vfact=1.d+0
   if (nact) 340,340,280
   
   !     SET THE RESIDUALS OF THE KUHN-TUCKER CONDITIONS FOR GENERAL X.
   
   250 iflag=2
   iws=iww-n
   do i=1,n
      iw=iww+i
      w(iw)=grad(i)
      if (lql) goto 259
      id=iwd+i
      w(id)=zero
      do j=i,n
         w(id)=w(id)+g(i,j)*x(j)
      end do
      do j=1,i
         id=iwd+j
         w(iw)=w(iw)+g(j,i)*w(id)
      end do
      cycle
      259 do j=1,n
         w(iw)=w(iw)+g(i,j)*x(j)
      end do
   end do
   if (nact == 0) goto 340
   do k=1,nact
      kk=iact(k)
      is=iws+k
      if (kk > m) goto 265
      w(is)=b(kk)
      do i=1,n
         iw=iww+i
         w(iw)=w(iw)-w(k)*a(kk,i)
         w(is)=w(is)-x(i)*a(kk,i)
      end do
      cycle
      265 if (kk > mn) goto 266
      k1=kk-m
      iw=iww+k1
      w(iw)=w(iw)-w(k)
      w(is)=xl(k1)-x(k1)
      cycle
      266 k1=kk-mn
      iw=iww+k1
      w(iw)=w(iw)+w(k)
      w(is)=-xu(k1)+x(k1)
   end do
   
   !     PRE-MULTIPLY THE VECTOR IN THE S-PARTITION OF W BY THE
   !     INVERS OF R TRANSPOSE.
   
   280 ir=iwr
   ip=iww+1
   ipp=iww+n
   il=iws+1
   iu=iws+nact
   do i=il,iu
      sum=zero
      if (i == il) goto 300
      ju=i-1
      do j=il,ju
         ir=ir+1
         sum=sum+w(ir)*w(j)
      end do
      300 ir=ir+1
      w(i)=(w(i)-sum)/w(ir)
   end do
   
   !     SHIFT X TO SATISFY THE ACTIVE CONSTRAINTS AND MAKE THE
   !     CORRESPONDING CHANGE TO THE GRADIENT RESIDUALS.
   
   do i=1,n
      iz=iwz+i
      sum=zero
      do j=il,iu
         sum=sum+w(j)*w(iz)
         iz=iz+n
      end do
      x(i)=x(i)+sum
      if (lql) goto 329
      id=iwd+i
      w(id)=zero
      do j=i,n
         w(id)=w(id)+g(i,j)*sum
      end do
      iw=iww+i
      do j=1,i
         id=iwd+j
         w(iw)=w(iw)+g(j,i)*w(id)
      end do
      cycle
      329 do j=1,n
         iw=iww+j
         w(iw)=w(iw)+sum*g(i,j)
      end do
   end do
   
   !     FORM THE SCALAR PRODUCT OF THE CURRENT GRADIENT RESIDUALS
   !     WITH EACH COLUMN OF Z.
   
   340 kflag=1
   goto 930
   350 if (nact == n) goto 380
   
   !     SHIFT X SO THAT IT SATISFIES THE REMAINING KUHN-TUCKER
   !     CONDITIONS.
   
   il=iws+nact+1
   iza=iwz+nact*n
   do i=1,n
      sum=zero
      iz=iza+i
      do j=il,iww
         sum=sum+w(iz)*w(j)
         iz=iz+n
      end do
      x(i)=x(i)-sum
   end do
   info=0
   if (nact == 0) goto 410
   
   !     UPDATE THE LAGRANGE MULTIPLIERS.
   
   380 lflag=3
   goto 740
   390 do k=1,nact
      iw=iww+k
      w(k)=w(k)+w(iw)
   end do
   
   !     REVISE THE VALUES OF XMAG.
   !     BRANCH IF ITERATIVE REFINEMENT IS REQUIRED.
   
   410 jflag=1
   goto 910
   420 if (iflag == itref) goto 250
   
   !     DELETE A CONSTRAINT IF A LAGRANGE MULTIPLIER OF AN
   !     INEQUALITY CONSTRAINT IS NEGATIVE.
   
   kdrop=0
   goto 440
   430 kdrop=kdrop+1
   if (w(kdrop) >= zero) goto 440
   if (iact(kdrop) <= meq) goto 440
   nu=nact
   mflag=1
   goto 800
   440 if (kdrop < nact) goto 430
   
   !     SEEK THE GREATEAST NORMALISED CONSTRAINT VIOLATION, DISREGARDING
   !     ANY THAT MAY BE DUE TO COMPUTER ROUNDING ERRORS.
   
   450 cvmax=zero
   if (m <= 0) goto 481
   do k=1,m
      ia=iwa+k
      if (w(ia) <= zero) cycle
      sum=-b(k)
      do i=1,n
         sum=sum+x(i)*a(k,i)
      end do
      sumx=-sum*w(ia)
      if (k <= meq) sumx=abs(sumx)
      if (sumx <= cvmax) cycle
      temp=abs(b(k))
      do i=1,n
         temp=temp+abs(x(i)*a(k,i))
      end do
      tempa=temp+abs(sum)
      if (tempa <= temp) cycle
      temp=temp+onha*abs(sum)
      if (temp <= tempa) cycle
      cvmax=sumx
      res=sum
      knext=k
   end do
   481 do k=1,n
      lower=.true.
      ia=iwa+m+k
      if (w(ia) <= zero) goto 485
      sum=xl(k)-x(k)
      if (sum) 482,485,483
      482 sum=x(k)-xu(k)
      lower=.false.
      483 if (sum <= cvmax) goto 485
      cvmax=sum
      res=-sum
      knext=k+m
      if (lower) goto 485
      knext=k+mn
   end do
   485 continue
   
   !     TEST FOR CONVERGENCE
   
   info=0
   if (cvmax <= vsmall) goto 700
   
   !     RETURN IF, DUE TO ROUNDING ERRORS, THE ACTUAL CHANGE IN
   !     X MAY NOT INCREASE THE OBJECTIVE FUNCTION
   
   jfinc=jfinc+1
   if (jfinc == 0) goto 510
   if (jfinc /= ifinc) goto 530
   fdiff=zero
   fdiffa=zero
   do i=1,n
      sum=two*grad(i)
      sumx=abs(sum)
      if (lql) goto 489
      id=iwd+i
      w(id)=zero
      do j=i,n
         ix=iwx+j
         w(id)=w(id)+g(i,j)*(w(ix)+x(j))
      end do
      do j=1,i
         id=iwd+j
         temp=g(j,i)*w(id)
         sum=sum+temp
         sumx=sumx+abs(temp)
      end do
      goto 495
      489 do j=1,n
         ix=iwx+j
         temp=g(i,j)*(w(ix)+x(j))
         sum=sum+temp
         sumx=sumx+abs(temp)
      end do
      495 ix=iwx+i
      fdiff=fdiff+sum*(x(i)-w(ix))
      fdiffa=fdiffa+sumx*abs(x(i)-w(ix))
   end do
   info=2
   sum=fdiffa+fdiff
   if (sum <= fdiffa) goto 700
   temp=fdiffa+onha*fdiff
   if (temp <= sum) goto 700
   jfinc=0
   info=0
   510 do i=1,n
      ix=iwx+i
      w(ix)=x(i)
   end do
   
   !     FORM THE SCALAR PRODUCT OF THE NEW CONSTRAINT NORMAL WITH EACH
   !     COLUMN OF Z. PARNEW WILL BECOME THE LAGRANGE MULTIPLIER OF
   !     THE NEW CONSTRAINT.
   
   530 iterc=iterc+1
   if (iterc <= maxit) goto 531
   info=1
   goto 710
   531 continue
   iws=iwr+(nact+nact*nact)/2
   if (knext > m) goto 541
   do i=1,n
      iw=iww+i
      w(iw)=a(knext,i)
   end do
   goto 549
   541 do i=1,n
      iw=iww+i
      w(iw)=zero
   end do
   k1=knext-m
   if (k1 > n) goto 545
   iw=iww+k1
   w(iw)=one
   iz=iwz+k1
   do i=1,n
      is=iws+i
      w(is)=w(iz)
      iz=iz+n
   end do
   goto 550
   545 k1=knext-mn
   iw=iww+k1
   w(iw)=-one
   iz=iwz+k1
   do i=1,n
      is=iws+i
      w(is)=-w(iz)
      iz=iz+n
   end do
   goto 550
   549 kflag=2
   goto 930
   550 parnew=zero
   
   !     APPLY GIVENS ROTATIONS TO MAKE THE LAST (N-NACT-2) SCALAR
   !     PRODUCTS EQUAL TO ZERO.
   
   if (nact == n) goto 570
   nu=n
   nflag=1
   goto 860
   
   !     BRANCH IF THERE IS NO NEED TO DELETE A CONSTRAINT.
   
   560 is=iws+nact
   if (nact == 0) goto 640
   suma=zero
   sumb=zero
   sumc=zero
   iz=iwz+nact*n
   do i=1,n
      iz=iz+1
      iw=iww+i
      suma=suma+w(iw)*w(iz)
      sumb=sumb+abs(w(iw)*w(iz))
      sumc=sumc+w(iz)**2
   end do
   temp=sumb+.1d+0*abs(suma)
   tempa=sumb+.2d+0*abs(suma)
   if (temp <= sumb) goto 570
   if (tempa <= temp) goto 570
   if (sumb > vsmall) goto 5
   goto 570
   5 sumc=sqrt(sumc)
   ia=iwa+knext
   if (knext <= m) sumc=sumc/w(ia)
   temp=sumc+.1d+0*abs(suma)
   tempa=sumc+.2d+0*abs(suma)
   if (temp <= sumc) goto 567
   if (tempa <= temp) goto 567
   goto 640
   
   !     CALCULATE THE MULTIPLIERS FOR THE NEW CONSTRAINT NORMAL
   !     EXPRESSED IN TERMS OF THE ACTIVE CONSTRAINT NORMALS.
   !     THEN WORK OUT WHICH CONTRAINT TO DROP.
   
   567 lflag=4
   goto 740
   570 lflag=1
   goto 740
   
   !     COMPLETE THE TEST FOR LINEARLY DEPENDENT CONSTRAINTS.
   
   571 if (knext > m) goto 574
   do i=1,n
      suma=a(knext,i)
      sumb=abs(suma)
      if (nact == 0) goto 581
      do k=1,nact
         kk=iact(k)
         if (kk <= m) goto 568
         kk=kk-m
         temp=zero
         if (kk == i) temp=w(iww+kk)
         kk=kk-n
         if (kk == i) temp=-w(iww+kk)
         goto 569
         568 continue
         iw=iww+k
         temp=w(iw)*a(kk,i)
         569 continue
         suma=suma-temp
         sumb=sumb+abs(temp)
      end do
      581 if (suma <= vsmall) cycle
      temp=sumb+.1d+0*abs(suma)
      tempa=sumb+.2d+0*abs(suma)
      if (temp <= sumb) cycle
      if (tempa <= temp) cycle
      goto 630
   end do
   lflag=1
   goto 775
   574 k1=knext-m
   if (k1 > n) k1=k1-n
   do i=1,n
      suma=zero
      if (i /= k1) goto 575
      suma=one
      if (knext > mn) suma=-one
      575 sumb=abs(suma)
      if (nact == 0) goto 582
      do k=1,nact
         kk=iact(k)
         if (kk <= m) goto 579
         kk=kk-m
         temp=zero
         if (kk == i) temp=w(iww+kk)
         kk=kk-n
         if (kk == i) temp=-w(iww+kk)
         goto 576
         579 iw=iww+k
         temp=w(iw)*a(kk,i)
         576 suma=suma-temp
      end do
      577 sumb=sumb+abs(temp)
      582 temp=sumb+.1d+0*abs(suma)
      tempa=sumb+.2d+0*abs(suma)
      if (temp <= sumb) cycle
      if (tempa <= temp) cycle
      goto 630
   end do
   lflag=1
   goto 775
   
   !     BRANCH IF THE CONTRAINTS ARE INCONSISTENT.
   
   580 info=-knext
   if (kdrop == 0) goto 700
   parinc=ratio
   parnew=parinc
   
   !     REVISE THE LAGRANGE MULTIPLIERS OF THE ACTIVE CONSTRAINTS.
   
   590 if (nact == 0) goto 601
   do k=1,nact
      iw=iww+k
      w(k)=w(k)-parinc*w(iw)
      if (iact(k) > meq) w(k)=fmax(zero,w(k))
   end do
   601 if (kdrop == 0) goto 680
   
   !     DELETE THE CONSTRAINT TO BE DROPPED.
   !     SHIFT THE VECTOR OF SCALAR PRODUCTS.
   !     THEN, IF APPROPRIATE, MAKE ONE MORE SCALAR PRODUCT ZERO.
   
   nu=nact+1
   mflag=2
   goto 800
   610 iws=iws-nact-1
   nu=min0(n,nu)
   do i=1,nu
      is=iws+i
      j=is+nact
      w(is)=w(j+1)
   end do
   nflag=2
   goto 860
   
   !     CALCULATE THE STEP TO THE VIOLATED CONSTRAINT.
   
   630 is=iws+nact
   640 sumy=w(is+1)
   step=-res/sumy
   parinc=step/sumy
   if (nact == 0) goto 660
   
   !     CALCULATE THE CHANGES TO THE LAGRANGE MULTIPLIERS, AND REDUCE
   !     THE STEP ALONG THE NEW SEARCH DIRECTION IF NECESSARY.
   
   lflag=2
   goto 740
   650 if (kdrop == 0) goto 660
   temp=one-ratio/parinc
   if (temp <= zero) kdrop=0
   if (kdrop == 0) goto 660
   step=ratio*sumy
   parinc=ratio
   res=temp*res
   
   !     UPDATE X AND THE LAGRANGE MULTIPIERS.
   !     DROP A CONSTRAINT IF THE FULL STEP IS NOT TAKEN.
   
   660 iwy=iwz+nact*n
   do i=1,n
      iy=iwy+i
      x(i)=x(i)+step*w(iy)
   end do
   parnew=parnew+parinc
   if (nact >= 1) goto 590
   
   !     ADD THE NEW CONSTRAINT TO THE ACTIVE SET.
   
   680 nact=nact+1
   w(nact)=parnew
   iact(nact)=knext
   ia=iwa+knext
   if (knext > mn) ia=ia-n
   w(ia)=-w(ia)
   
   !     ESTIMATE THE MAGNITUDE OF X. THEN BEGIN A NEW ITERATION,
   !     RE-INITILISING X IF THIS MAGNITUDE IS SMALL.
   
   jflag=2
   goto 910
   690 if (sum < (xmagr*xmag)) goto 230
   if (itref) 450,450,250
   
   !     INITIATE ITERATIVE REFINEMENT IF IT HAS NOT YET BEEN USED,
   !     OR RETURN AFTER RESTORING THE DIAGONAL ELEMENTS OF G.
   
   700 if (iterc == 0) goto 710
   itref=itref+1
   jfinc=-1
   if (itref == 1) goto 250
   710 if (.not. lql) return
   do i=1,n
      id=iwd+i
      g(i,i)=w(id)
   end do
   730 return
   
   
   !     THE REMAINING INSTRUCTIONS ARE USED AS SUBROUTINES.
   
   
   !********************************************************************
   
   
   !     CALCULATE THE LAGRANGE MULTIPLIERS BY PRE-MULTIPLYING THE
   !     VECTOR IN THE S-PARTITION OF W BY THE INVERSE OF R.
   
   740 ir=iwr+(nact+nact*nact)/2
   i=nact
   sum=zero
   goto 770
   750 ira=ir-1
   sum=zero
   if (nact == 0) goto 761
   do j=i,nact
      iw=iww+j
      sum=sum+w(ira)*w(iw)
      ira=ira+j
   end do
   761 ir=ir-i
   i=i-1
   770 iw=iww+i
   is=iws+i
   w(iw)=(w(is)-sum)/w(ir)
   if (i > 1) goto 750
   if (lflag == 3) goto 390
   if (lflag == 4) goto 571
   
   !     CALCULATE THE NEXT CONSTRAINT TO DROP.
   
   775 ip=iww+1
   ipp=iww+nact
   kdrop=0
   if (nact == 0) goto 791
   do k=1,nact
      if (iact(k) <= meq) cycle
      iw=iww+k
      if ((res*w(iw)) >= zero) cycle
      temp=w(k)/w(iw)
      if (kdrop == 0) goto 780
      if (abs(temp) >= abs(ratio)) cycle
      780 kdrop=k
      ratio=temp
   end do
   791 goto (580,650), lflag
   
   
   !********************************************************************
   
   
   !     DROP THE CONSTRAINT IN POSITION KDROP IN THE ACTIVE SET.
   
   800 ia=iwa+iact(kdrop)
   if (iact(kdrop) > mn) ia=ia-n
   w(ia)=-w(ia)
   if (kdrop == nact) goto 850
   
   !     SET SOME INDICES AND CALCULATE THE ELEMENTS OF THE NEXT
   !     GIVENS ROTATION.
   
   iz=iwz+kdrop*n
   ir=iwr+(kdrop+kdrop*kdrop)/2
   810 ira=ir
   ir=ir+kdrop+1
   temp=fmax(abs(w(ir-1)),abs(w(ir)))
   sum=temp*sqrt((w(ir-1)/temp)**2+(w(ir)/temp)**2)
   ga=w(ir-1)/sum
   gb=w(ir)/sum
   
   !     EXCHANGE THE COLUMNS OF R.
   
   do i=1,kdrop
      ira=ira+1
      j=ira-kdrop
      temp=w(ira)
      w(ira)=w(j)
      w(j)=temp
   end do
   w(ir)=zero
   
   !     APPLY THE ROTATION TO THE ROWS OF R.
   
   w(j)=sum
   kdrop=kdrop+1
   do i=kdrop,nu
      temp=ga*w(ira)+gb*w(ira+1)
      w(ira+1)=ga*w(ira+1)-gb*w(ira)
      w(ira)=temp
      ira=ira+i
   end do
   
   !     APPLY THE ROTATION TO THE COLUMNS OF Z.
   
   do i=1,n
      iz=iz+1
      j=iz-n
      temp=ga*w(j)+gb*w(iz)
      w(iz)=ga*w(iz)-gb*w(j)
      w(j)=temp
   end do
   
   !     REVISE IACT AND THE LAGRANGE MULTIPLIERS.
   
   iact(kdrop-1)=iact(kdrop)
   w(kdrop-1)=w(kdrop)
   if (kdrop < nact) goto 810
   850 nact=nact-1
   goto (250,610), mflag
   
   
   !********************************************************************
   
   
   !     APPLY GIVENS ROTATION TO REDUCE SOME OF THE SCALAR
   !     PRODUCTS IN THE S-PARTITION OF W TO ZERO.
   
   860 iz=iwz+nu*n
   870 iz=iz-n
   880 is=iws+nu
   nu=nu-1
   if (nu == nact) goto 900
   if (w(is) == zero) goto 870
   temp=fmax(abs(w(is-1)),abs(w(is)))
   sum=temp*sqrt((w(is-1)/temp)**2+(w(is)/temp)**2)
   ga=w(is-1)/sum
   gb=w(is)/sum
   w(is-1)=sum
   do i=1,n
      k=iz+n
      temp=ga*w(iz)+gb*w(k)
      w(k)=ga*w(k)-gb*w(iz)
      w(iz)=temp
      iz=iz-1
   end do
   goto 880
   900 goto (560,630), nflag
   
   
   !********************************************************************
   
   
   !     CALCULATE THE MAGNITUDE OF X AN REVISE XMAG.
   
   910 sum=zero
   do i=1,n
      sum=sum+abs(x(i))*vfact*(abs(grad(i))+abs(g(i,i)*x(i)))
      if (lql) cycle
      if (sum < 1.d-30) cycle
      vfact=1.d-10*vfact
      sum=1.d-10*sum
      xmag=1.d-10*xmag
   end do
   925 xmag=fmax(xmag,sum)
   goto (420,690), jflag
   
   
   !********************************************************************
   
   
   !     PRE-MULTIPLY THE VECTOR IN THE W-PARTITION OF W BY Z TRANSPOSE.
   
   930 jl=iww+1
   iz=iwz
   do i=1,n
      is=iws+i
      w(is)=zero
      iwwn=iww+n
      do j=jl,iwwn
         iz=iz+1
         w(is)=w(is)+w(iz)*w(j)
      end do
   end do
   goto (350,550), kflag
   return
   end subroutine ql0002 
