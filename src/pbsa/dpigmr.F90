!DECK DPIGMR
      SUBROUTINE DPIGMR (N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP, NRSTS,&
         JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR, RPAR, IPAR, WK,&
         DL, RHOL, NRMAX, B, BNRM, X, XL, ITOL, TOL, NELT, IA, JA, A, &
         ISYM, IUNIT, IFLAG, ERR)
!***BEGIN PROLOGUE  DPIGMR
!***SUBSIDIARY
!***PURPOSE  Internal routine for DGMRES.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      DOUBLE PRECISION (SPIGMR-S, DPIGMR-D)
!***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
!***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
!           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
!           Seager, Mark K., (LLNL), seager@llnl.gov
!             Lawrence Livermore National Laboratory
!             PO Box 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!***DESCRIPTION
!         This routine solves the linear system A * Z = R0 using a
!         scaled preconditioned version of the generalized minimum
!         residual method.  An initial guess of Z = 0 is assumed.
!
! *Usage:
!      INTEGER N, JSCAL, MAXL, MAXLP1, KMP, NRSTS, JPRE, NMSL, LGMR
!      INTEGER IPAR(USER DEFINED), NRMAX, ITOL, NELT, IA(NELT), JA(NELT)
!      INTEGER ISYM, IUNIT, IFLAG
!      DOUBLE PRECISION R0(N), SR(N), SZ(N), Z(N), V(N,MAXLP1),
!     $                 HES(MAXLP1,MAXL), Q(2*MAXL), RPAR(USER DEFINED),
!     $                 WK(N), DL(N), RHOL, B(N), BNRM, X(N), XL(N),
!     $                 TOL, A(NELT), ERR
!      EXTERNAL MATVEC, MSOLVE
!
!      CALL DPIGMR(N, R0, SR, SZ, JSCAL, MAXL, MAXLP1, KMP,
!     $     NRSTS, JPRE, MATVEC, MSOLVE, NMSL, Z, V, HES, Q, LGMR,
!     $     RPAR, IPAR, WK, DL, RHOL, NRMAX, B, BNRM, X, XL,
!     $     ITOL, TOL, NELT, IA, JA, A, ISYM, IUNIT, IFLAG, ERR)
!
! *Arguments:
! N      :IN       Integer
!         The order of the matrix A, and the lengths
!         of the vectors SR, SZ, R0 and Z.
! R0     :IN       Double Precision R0(N)
!         R0 = the right hand side of the system A*Z = R0.
!         R0 is also used as workspace when computing
!         the final approximation.
!         (R0 is the same as V(*,MAXL+1) in the call to DPIGMR.)
! SR     :IN       Double Precision SR(N)
!         SR is a vector of length N containing the non-zero
!         elements of the diagonal scaling matrix for R0.
! SZ     :IN       Double Precision SZ(N)
!         SZ is a vector of length N containing the non-zero
!         elements of the diagonal scaling matrix for Z.
! JSCAL  :IN       Integer
!         A flag indicating whether arrays SR and SZ are used.
!         JSCAL=0 means SR and SZ are not used and the
!                 algorithm will perform as if all
!                 SR(i) = 1 and SZ(i) = 1.
!         JSCAL=1 means only SZ is used, and the algorithm
!                 performs as if all SR(i) = 1.
!         JSCAL=2 means only SR is used, and the algorithm
!                 performs as if all SZ(i) = 1.
!         JSCAL=3 means both SR and SZ are used.
! MAXL   :IN       Integer
!         The maximum allowable order of the matrix H.
! MAXLP1 :IN       Integer
!         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
! KMP    :IN       Integer
!         The number of previous vectors the new vector VNEW
!         must be made orthogonal to.  (KMP .le. MAXL)
! NRSTS  :IN       Integer
!         Counter for the number of restarts on the current
!         call to DGMRES.  If NRSTS .gt. 0, then the residual
!         R0 is already scaled, and so scaling of it is
!         not necessary.
! JPRE   :IN       Integer
!         Preconditioner type flag.
! MATVEC :EXT      External.
!         Name of a routine which performs the matrix vector multiply
!         Y = A*X given A and X.  The name of the MATVEC routine must
!         be declared external in the calling program.  The calling
!         sequence to MATVEC is:
!             CALL MATVEC(N, X, Y, NELT, IA, JA, A, ISYM)
!         where N is the number of unknowns, Y is the product A*X
!         upon return, X is an input vector, and NELT is the number of
!         non-zeros in the SLAP IA, JA, A storage for the matrix A.
!         ISYM is a flag which, if non-zero, denotes that A is
!         symmetric and only the lower or upper triangle is stored.
! MSOLVE :EXT      External.
!         Name of the routine which solves a linear system Mz = r for
!         z given r with the preconditioning matrix M (M is supplied via
!         RPAR and IPAR arrays.  The name of the MSOLVE routine must
!         be declared external in the calling program.  The calling
!         sequence to MSOLVE is:
!             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
!         Where N is the number of unknowns, R is the right-hand side
!         vector and Z is the solution upon return.  NELT, IA, JA, A and
!         ISYM are defined as below.  RPAR is a double precision array
!         that can be used to pass necessary preconditioning information
!         and/or workspace to MSOLVE.  IPAR is an integer work array
!         for the same purpose as RPAR.
! NMSL   :OUT      Integer
!         The number of calls to MSOLVE.
! Z      :OUT      Double Precision Z(N)
!         The final computed approximation to the solution
!         of the system A*Z = R0.
! V      :OUT      Double Precision V(N,MAXLP1)
!         The N by (LGMR+1) array containing the LGMR
!         orthogonal vectors V(*,1) to V(*,LGMR).
! HES    :OUT      Double Precision HES(MAXLP1,MAXL)
!         The upper triangular factor of the QR decomposition
!         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
!         entries are the scaled inner-products of A*V(*,I)
!         and V(*,K).
! Q      :OUT      Double Precision Q(2*MAXL)
!         A double precision array of length 2*MAXL containing the
!         components of the Givens rotations used in the QR
!         decomposition of HES.  It is loaded in DHEQR and used in
!         DHELS.
! LGMR   :OUT      Integer
!         The number of iterations performed and
!         the current order of the upper Hessenberg
!         matrix HES.
! RPAR   :IN       Double Precision RPAR(USER DEFINED)
!         Double Precision workspace passed directly to the MSOLVE
!         routine.
! IPAR   :IN       Integer IPAR(USER DEFINED)
!         Integer workspace passed directly to the MSOLVE routine.
! WK     :IN       Double Precision WK(N)
!         A double precision work array of length N used by routines
!         MATVEC and MSOLVE.
! DL     :INOUT    Double Precision DL(N)
!         On input, a double precision work array of length N used for
!         calculation of the residual norm RHO when the method is
!         incomplete (KMP.lt.MAXL), and/or when using restarting.
!         On output, the scaled residual vector RL.  It is only loaded
!         when performing restarts of the Krylov iteration.
! RHOL   :OUT      Double Precision
!         A double precision scalar containing the norm of the final
!         residual.
! NRMAX  :IN       Integer
!         The maximum number of restarts of the Krylov iteration.
!         NRMAX .gt. 0 means restarting is active, while
!         NRMAX = 0 means restarting is not being used.
! B      :IN       Double Precision B(N)
!         The right hand side of the linear system A*X = b.
! BNRM   :IN       Double Precision
!         The scaled norm of b.
! X      :IN       Double Precision X(N)
!         The current approximate solution as of the last
!         restart.
! XL     :IN       Double Precision XL(N)
!         An array of length N used to hold the approximate
!         solution X(L) when ITOL=11.
! ITOL   :IN       Integer
!         A flag to indicate the type of convergence criterion
!         used.  See the driver for its description.
! TOL    :IN       Double Precision
!         The tolerance on residuals R0-A*Z in scaled norm.
! NELT   :IN       Integer
!         The length of arrays IA, JA and A.
! IA     :IN       Integer IA(NELT)
!         An integer array of length NELT containing matrix data.
!         It is passed directly to the MATVEC and MSOLVE routines.
! JA     :IN       Integer JA(NELT)
!         An integer array of length NELT containing matrix data.
!         It is passed directly to the MATVEC and MSOLVE routines.
! A      :IN       Double Precision A(NELT)
!         A double precision array of length NELT containing matrix
!         data. It is passed directly to the MATVEC and MSOLVE routines.
! ISYM   :IN       Integer
!         A flag to indicate symmetric matrix storage.
!         If ISYM=0, all non-zero entries of the matrix are
!         stored.  If ISYM=1, the matrix is symmetric and
!         only the upper or lower triangular part is stored.
! IUNIT  :IN       Integer
!         The i/o unit number for writing intermediate residual
!         norm values.
! IFLAG  :OUT      Integer
!         An integer error flag..
!         0 means convergence in LGMR iterations, LGMR.le.MAXL.
!         1 means the convergence test did not pass in MAXL
!           iterations, but the residual norm is .lt. norm(R0),
!           and so Z is computed.
!         2 means the convergence test did not pass in MAXL
!           iterations, residual .ge. norm(R0), and Z = 0.
! ERR    :OUT      Double Precision.
!         Error estimate of error in final approximate solution, as
!         defined by ITOL.
!
! *Cautions:
!     This routine will attempt to write to the Fortran logical output
!     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
!     this logical unit is attached to a file or terminal before calling
!     this routine with a non-zero value for IUNIT.  This routine does
!     not check for the validity of a non-zero IUNIT unit number.
!
!***SEE ALSO  DGMRES
!***ROUTINES CALLED  DAXPY, DCOPY, DHELS, DHEQR, DNRM2, DORTH, DRLCAL,
!                    DSCAL, ISDGMR
!***REVISION HISTORY  (YYMMDD)
!   890404  DATE WRITTEN
!   890404  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Removed MATVEC and MSOLVE from ROUTINES CALLED list.  (FNF)
!   910506  Made subsidiary to DGMRES.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!***END PROLOGUE  DPIGMR
!         The following is for optimized compilation on LLNL/LTSS Crays.
!LLL. OPTIMIZE
!     .. Scalar Arguments ..
      DOUBLE PRECISION BNRM, ERR, RHOL, TOL
      INTEGER IFLAG, ISYM, ITOL, IUNIT, JPRE, JSCAL, KMP, LGMR, MAXL, &
              MAXLP1, N, NELT, NMSL, NRMAX, NRSTS
!     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), B(*), DL(*), HES(MAXLP1,*), Q(*), R0(*), &
                       RPAR(*), SR(*), SZ(*), V(N,*), WK(*), X(*), &
                       XL(*), Z(*)
      INTEGER IA(NELT), IPAR(*), JA(NELT)
!     .. Subroutine Arguments ..
      EXTERNAL MATVEC, MSOLVE
!     .. Local Scalars ..
      DOUBLE PRECISION C, DLNRM, PROD, R0NRM, RHO, S, SNORMW, TEM
      INTEGER I, I2, INFO, IP1, ITER, ITMAX, J, K, LL, LLP1
!     .. External Functions ..
      DOUBLE PRECISION DNRM2
      INTEGER ISDGMR
      EXTERNAL DNRM2, ISDGMR
!     .. External Subroutines ..
      EXTERNAL DAXPY, DCOPY, DHELS, DHEQR, DORTH, DRLCAL, DSCAL
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!***FIRST EXECUTABLE STATEMENT  DPIGMR
!
!         Zero out the Z array.
!
      DO 5 I = 1,N
         Z(I) = 0
 5    CONTINUE
!
      IFLAG = 0
      LGMR = 0
      NMSL = 0
!         Load ITMAX, the maximum number of iterations.
      ITMAX =(NRMAX+1)*MAXL
!   -------------------------------------------------------------------
!         The initial residual is the vector R0.
!         Apply left precon. if JPRE < 0 and this is not a restart.
!         Apply scaling to R0 if JSCAL = 2 or 3.
!   -------------------------------------------------------------------
      IF ((JPRE .LT. 0) .AND.(NRSTS .EQ. 0)) THEN
         CALL DCOPY(N, R0, 1, WK, 1)
         CALL MSOLVE(N, WK, R0, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      IF (((JSCAL.EQ.2) .OR.(JSCAL.EQ.3)) .AND.(NRSTS.EQ.0)) THEN
         DO 10 I = 1,N
            V(I,1) = R0(I)*SR(I)
 10      CONTINUE
      ELSE
         DO 20 I = 1,N
            V(I,1) = R0(I)
 20      CONTINUE
      ENDIF
      R0NRM = DNRM2(N, V, 1)
      ITER = NRSTS*MAXL
!
!         Call stopping routine ISDGMR.
!
      IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE, &
          NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, V(1,1), Z, WK, &
          RPAR, IPAR, R0NRM, BNRM, SR, SZ, JSCAL, &
          KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM, &
          HES, JPRE) .NE. 0) RETURN
      TEM = 1.0D0/R0NRM
      CALL DSCAL(N, TEM, V(1,1), 1)
!
!         Zero out the HES array.
!
      DO 50 J = 1,MAXL
         DO 40 I = 1,MAXLP1
            HES(I,J) = 0
 40      CONTINUE
 50   CONTINUE
!   -------------------------------------------------------------------
!         Main loop to compute the vectors V(*,2) to V(*,MAXL).
!         The running product PROD is needed for the convergence test.
!   -------------------------------------------------------------------
      PROD = 1
      DO 90 LL = 1,MAXL
         LGMR = LL
!   -------------------------------------------------------------------
!        Unscale  the  current V(LL)  and store  in WK.  Call routine
!        MSOLVE    to   compute(M-inverse)*WK,   where    M   is  the
!        preconditioner matrix.  Save the answer in Z.   Call routine
!        MATVEC to compute  VNEW  = A*Z,  where  A is  the the system
!        matrix.  save the answer in  V(LL+1).  Scale V(LL+1).   Call
!        routine DORTH  to  orthogonalize the    new vector VNEW   =
!        V(*,LL+1).  Call routine DHEQR to update the factors of HES.
!   -------------------------------------------------------------------
        IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
           DO 60 I = 1,N
              WK(I) = V(I,LL)/SZ(I)
 60        CONTINUE
        ELSE
           CALL DCOPY(N, V(1,LL), 1, WK, 1)
        ENDIF
        IF (JPRE .GT. 0) THEN
           CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
           NMSL = NMSL + 1
           CALL MATVEC(N, Z, V(1,LL+1), NELT, IA, JA, A, ISYM)
        ELSE
           CALL MATVEC(N, WK, V(1,LL+1), NELT, IA, JA, A, ISYM)
        ENDIF
        IF (JPRE .LT. 0) THEN
           CALL DCOPY(N, V(1,LL+1), 1, WK, 1)
           CALL MSOLVE(N,WK,V(1,LL+1),NELT,IA,JA,A,ISYM,RPAR,IPAR)
           NMSL = NMSL + 1
        ENDIF
        IF ((JSCAL .EQ. 2) .OR.(JSCAL .EQ. 3)) THEN
           DO 65 I = 1,N
              V(I,LL+1) = V(I,LL+1)*SR(I)
 65        CONTINUE
        ENDIF
        CALL DORTH(V(1,LL+1), V, HES, N, LL, MAXLP1, KMP, SNORMW)
        HES(LL+1,LL) = SNORMW
        CALL DHEQR(HES, MAXLP1, LL, Q, INFO, LL)
        IF (INFO .EQ. LL) GO TO 120
!   -------------------------------------------------------------------
!         Update RHO, the estimate of the norm of the residual R0-A*ZL.
!         If KMP <  MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
!         necessarily orthogonal for LL > KMP.  The vector DL must then
!         be computed, and its norm used in the calculation of RHO.
!   -------------------------------------------------------------------
        PROD = PROD*Q(2*LL)
        RHO = ABS(PROD*R0NRM)
        IF ((LL.GT.KMP) .AND.(KMP.LT.MAXL)) THEN
           IF (LL .EQ. KMP+1) THEN
              CALL DCOPY(N, V(1,1), 1, DL, 1)
              DO 75 I = 1,KMP
                 IP1 = I + 1
                 I2 = I*2
                 S = Q(I2)
                 C = Q(I2-1)
                 DO 70 K = 1,N
                    DL(K) = S*DL(K) + C*V(K,IP1)
 70              CONTINUE
 75           CONTINUE
           ENDIF
           S = Q(2*LL)
           C = Q(2*LL-1)/SNORMW
           LLP1 = LL + 1
           DO 80 K = 1,N
              DL(K) = S*DL(K) + C*V(K,LLP1)
 80        CONTINUE
           DLNRM = DNRM2(N, DL, 1)
           RHO = RHO*DLNRM
        ENDIF
        RHOL = RHO
!   -------------------------------------------------------------------
!         Test for convergence.  If passed, compute approximation ZL.
!         If failed and LL < MAXL, then continue iterating.
!   -------------------------------------------------------------------
        ITER = NRSTS*MAXL + LGMR
        IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE, &
            NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, DL, Z, WK, &
            RPAR, IPAR, RHOL, BNRM, SR, SZ, JSCAL, &
            KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM, &
            HES, JPRE) .NE. 0) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
!   -------------------------------------------------------------------
!         Rescale so that the norm of V(1,LL+1) is one.
!   -------------------------------------------------------------------
        TEM = 1.0D0/SNORMW
        CALL DSCAL(N, TEM, V(1,LL+1), 1)
 90   CONTINUE
 100  CONTINUE
      IF (RHO .LT. R0NRM) GO TO 150
 120  CONTINUE
      IFLAG = 2
!
!         Load approximate solution with zero.
!
      DO 130 I = 1,N
         Z(I) = 0
 130  CONTINUE
      RETURN
 150  IFLAG = 1
!
!         Tolerance not met, but residual norm reduced.
!
      IF (NRMAX .GT. 0) THEN
!
!        If performing restarting (NRMAX > 0)  calculate the residual
!        vector RL and  store it in the DL  array.  If the incomplete
!        version is being used (KMP < MAXL) then DL has  already been
!        calculated up to a scaling factor.   Use DRLCAL to calculate
!        the scaled residual vector.
!
         CALL DRLCAL(N, KMP, MAXL, MAXL, V, Q, DL, SNORMW, PROD, &
              R0NRM)
      ENDIF
!   -------------------------------------------------------------------
!         Compute the approximation ZL to the solution.  Since the
!         vector Z was used as workspace, and the initial guess
!         of the linear iteration is zero, Z must be reset to zero.
!   -------------------------------------------------------------------
 200  CONTINUE
      LL = LGMR
      LLP1 = LL + 1
      DO 210 K = 1,LLP1
         R0(K) = 0
 210  CONTINUE
      R0(1) = R0NRM
      CALL DHELS(HES, MAXLP1, LL, Q, R0)
      DO 220 K = 1,N
         Z(K) = 0
 220  CONTINUE
      DO 230 I = 1,LL
         CALL DAXPY(N, R0(I), V(1,I), 1, Z, 1)
 230  CONTINUE
      IF ((JSCAL .EQ. 1) .OR.(JSCAL .EQ. 3)) THEN
         DO 240 I = 1,N
            Z(I) = Z(I)/SZ(I)
 240     CONTINUE
      ENDIF
      IF (JPRE .GT. 0) THEN
         CALL DCOPY(N, Z, 1, WK, 1)
         CALL MSOLVE(N, WK, Z, NELT, IA, JA, A, ISYM, RPAR, IPAR)
         NMSL = NMSL + 1
      ENDIF
      RETURN
!------------- LAST LINE OF DPIGMR FOLLOWS ----------------------------
      END
