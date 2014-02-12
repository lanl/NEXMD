!DECK ISDGMR
      INTEGER FUNCTION ISDGMR (N, B, X, XL, NELT, IA, JA, A, ISYM, &
         MSOLVE, NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ, &
         RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL, KMP, LGMR, MAXL, &
         MAXLP1, V, Q, SNORMW, PROD, R0NRM, HES, JPRE)
!***BEGIN PROLOGUE  ISDGMR
!***SUBSIDIARY
!***PURPOSE  Generalized Minimum Residual Stop Test.
!            This routine calculates the stop test for the Generalized
!            Minimum RESidual (GMRES) iteration scheme.  It returns a
!            non-zero if the error estimate (the type of which is
!            determined by ITOL) is less than the user specified
!            tolerance TOL.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      DOUBLE PRECISION (ISSGMR-S, ISDGMR-D)
!***KEYWORDS  GMRES, LINEAR SYSTEM, SLAP, SPARSE, STOP TEST
!***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
!           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
!           Seager, Mark K., (LLNL), seager@llnl.gov
!             Lawrence Livermore National Laboratory
!             PO Box 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!***DESCRIPTION
!
! *Usage:
!      INTEGER N, NELT, IA(NELT), JA(NELT), ISYM, NMSL, ITOL
!      INTEGER ITMAX, ITER, IUNIT, IWORK(USER DEFINED), JSCAL
!      INTEGER KMP, LGMR, MAXL, MAXLP1, JPRE
!      DOUBLE PRECISION B(N), X(N), XL(MAXL), A(NELT), TOL, ERR,
!     $                 R(N), Z(N), DZ(N), RWORK(USER DEFINED),
!     $                 RNRM, BNRM, SB(N), SX(N), V(N,MAXLP1),
!     $                 Q(2*MAXL), SNORMW, PROD, R0NRM,
!     $                 HES(MAXLP1,MAXL)
!      EXTERNAL MSOLVE
!
!      IF (ISDGMR(N, B, X, XL, NELT, IA, JA, A, ISYM, MSOLVE,
!     $     NMSL, ITOL, TOL, ITMAX, ITER, ERR, IUNIT, R, Z, DZ,
!     $     RWORK, IWORK, RNRM, BNRM, SB, SX, JSCAL,
!     $     KMP, LGMR, MAXL, MAXLP1, V, Q, SNORMW, PROD, R0NRM,
!     $     HES, JPRE) .NE. 0) THEN ITERATION DONE
!
! *Arguments:
! N      :IN       Integer.
!         Order of the Matrix.
! B      :IN       Double Precision B(N).
!         Right-hand-side vector.
! X      :IN       Double Precision X(N).
!         Approximate solution vector as of the last restart.
! XL     :OUT      Double Precision XL(N)
!         An array of length N used to hold the approximate
!         solution as of the current iteration.  Only computed by
!         this routine when ITOL=11.
! NELT   :IN       Integer.
!         Number of Non-Zeros stored in A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Double Precision A(NELT).
!         These arrays contain the matrix data structure for A.
!         It could take any form.  See "Description", in the DGMRES,
!         DSLUGM and DSDGMR routines for more details.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the upper
!         or lower triangle of the matrix is stored.
! MSOLVE :EXT      External.
!         Name of a routine which solves a linear system Mz = r for  z
!         given r with the preconditioning matrix M (M is supplied via
!         RWORK and IWORK arrays.  The name of the MSOLVE routine must
!         be declared external in the calling program.  The calling
!         sequence to MSOLVE is:
!             CALL MSOLVE(N, R, Z, NELT, IA, JA, A, ISYM, RWORK, IWORK)
!         Where N is the number of unknowns, R is the right-hand side
!         vector and Z is the solution upon return.  NELT, IA, JA, A and
!         ISYM are defined as above.  RWORK is a double precision array
!         that can be used to pass necessary preconditioning information
!         and/or workspace to MSOLVE.  IWORK is an integer work array
!         for the same purpose as RWORK.
! NMSL   :INOUT    Integer.
!         A counter for the number of calls to MSOLVE.
! ITOL   :IN       Integer.
!         Flag to indicate the type of convergence criterion used.
!         ITOL=0  Means the  iteration stops when the test described
!                 below on  the  residual RL  is satisfied.  This is
!                 the  "Natural Stopping Criteria" for this routine.
!                 Other values  of   ITOL  cause  extra,   otherwise
!                 unnecessary, computation per iteration and     are
!                 therefore much less efficient.
!         ITOL=1  Means   the  iteration stops   when the first test
!                 described below on  the residual RL  is satisfied,
!                 and there  is either right  or  no preconditioning
!                 being used.
!         ITOL=2  Implies     that   the  user    is   using    left
!                 preconditioning, and the second stopping criterion
!                 below is used.
!         ITOL=3  Means the  iteration stops   when  the  third test
!                 described below on Minv*Residual is satisfied, and
!                 there is either left  or no  preconditioning begin
!                 used.
!         ITOL=11 is    often  useful  for   checking  and comparing
!                 different routines.  For this case, the  user must
!                 supply  the  "exact" solution or  a  very accurate
!                 approximation (one with  an  error much less  than
!                 TOL) through a common block,
!                     COMMON /DSLBLK/ SOLN( )
!                 If ITOL=11, iteration stops when the 2-norm of the
!                 difference between the iterative approximation and
!                 the user-supplied solution  divided by the  2-norm
!                 of the  user-supplied solution  is  less than TOL.
!                 Note that this requires  the  user to  set up  the
!                 "COMMON     /DSLBLK/ SOLN(LENGTH)"  in the calling
!                 routine.  The routine with this declaration should
!                 be loaded before the stop test so that the correct
!                 length is used by  the loader.  This procedure  is
!                 not standard Fortran and may not work correctly on
!                 your   system (although  it  has  worked  on every
!                 system the authors have tried).  If ITOL is not 11
!                 then this common block is indeed standard Fortran.
! TOL    :IN       Double Precision.
!         Convergence criterion, as described above.
! ITMAX  :IN       Integer.
!         Maximum number of iterations.
! ITER   :IN       Integer.
!         The iteration for which to check for convergence.
! ERR    :OUT      Double Precision.
!         Error estimate of error in final approximate solution, as
!         defined by ITOL.  Letting norm() denote the Euclidean
!         norm, ERR is defined as follows..
!
!         If ITOL=0, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
!                               for right or no preconditioning, and
!                         ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
!                                norm(SB*(M-inverse)*B),
!                               for left preconditioning.
!         If ITOL=1, then ERR = norm(SB*(B-A*X(L)))/norm(SB*B),
!                               since right or no preconditioning
!                               being used.
!         If ITOL=2, then ERR = norm(SB*(M-inverse)*(B-A*X(L)))/
!                                norm(SB*(M-inverse)*B),
!                               since left preconditioning is being
!                               used.
!         If ITOL=3, then ERR =  Max  |(Minv*(B-A*X(L)))(i)/x(i)|
!                               i=1,n
!         If ITOL=11, then ERR = norm(SB*(X(L)-SOLN))/norm(SB*SOLN).
! IUNIT  :IN       Integer.
!         Unit number on which to write the error at each iteration,
!         if this is desired for monitoring convergence.  If unit
!         number is 0, no writing will occur.
! R      :INOUT    Double Precision R(N).
!         Work array used in calling routine.  It contains
!         information necessary to compute the residual RL = B-A*XL.
! Z      :WORK     Double Precision Z(N).
!         Workspace used to hold the pseudo-residual M z = r.
! DZ     :WORK     Double Precision DZ(N).
!         Workspace used to hold temporary vector(s).
! RWORK  :WORK     Double Precision RWORK(USER DEFINED).
!         Double Precision array that can be used by MSOLVE.
! IWORK  :WORK     Integer IWORK(USER DEFINED).
!         Integer array that can be used by MSOLVE.
! RNRM   :IN       Double Precision.
!         Norm of the current residual.  Type of norm depends on ITOL.
! BNRM   :IN       Double Precision.
!         Norm of the right hand side.  Type of norm depends on ITOL.
! SB     :IN       Double Precision SB(N).
!         Scaling vector for B.
! SX     :IN       Double Precision SX(N).
!         Scaling vector for X.
! JSCAL  :IN       Integer.
!         Flag indicating if scaling arrays SB and SX are being
!         used in the calling routine DPIGMR.
!         JSCAL=0 means SB and SX are not used and the
!                 algorithm will perform as if all
!                 SB(i) = 1 and SX(i) = 1.
!         JSCAL=1 means only SX is used, and the algorithm
!                 performs as if all SB(i) = 1.
!         JSCAL=2 means only SB is used, and the algorithm
!                 performs as if all SX(i) = 1.
!         JSCAL=3 means both SB and SX are used.
! KMP    :IN       Integer
!         The number of previous vectors the new vector VNEW
!         must be made orthogonal to.  (KMP .le. MAXL)
! LGMR   :IN       Integer
!         The number of GMRES iterations performed on the current call
!         to DPIGMR (i.e., # iterations since the last restart) and
!         the current order of the upper Hessenberg
!         matrix HES.
! MAXL   :IN       Integer
!         The maximum allowable order of the matrix H.
! MAXLP1 :IN       Integer
!         MAXPL1 = MAXL + 1, used for dynamic dimensioning of HES.
! V      :IN       Double Precision V(N,MAXLP1)
!         The N by (LGMR+1) array containing the LGMR
!         orthogonal vectors V(*,1) to V(*,LGMR).
! Q      :IN       Double Precision Q(2*MAXL)
!         A double precision array of length 2*MAXL containing the
!         components of the Givens rotations used in the QR
!         decomposition of HES.
! SNORMW :IN       Double Precision
!         A scalar containing the scaled norm of VNEW before it
!         is renormalized in DPIGMR.
! PROD   :IN       Double Precision
!         The product s1*s2*...*sl = the product of the sines of the
!         Givens rotations used in the QR factorization of the
!         Hessenberg matrix HES.
! R0NRM  :IN       Double Precision
!         The scaled norm of initial residual R0.
! HES    :IN       Double Precision HES(MAXLP1,MAXL)
!         The upper triangular factor of the QR decomposition
!         of the (LGMR+1) by LGMR upper Hessenberg matrix whose
!         entries are the scaled inner-products of A*V(*,I)
!         and V(*,K).
! JPRE   :IN       Integer
!         Preconditioner type flag.
!         (See description of IGWK(4) in DGMRES.)
!
! *Description
!       When using the GMRES solver,  the preferred value  for ITOL
!       is 0.  This is due to the fact that when ITOL=0 the norm of
!       the residual required in the stopping test is  obtained for
!       free, since this value is already  calculated  in the GMRES
!       algorithm.   The  variable  RNRM contains the   appropriate
!       norm, which is equal to norm(SB*(RL - A*XL))  when right or
!       no   preconditioning is  being  performed,   and equal   to
!       norm(SB*Minv*(RL - A*XL))  when using left preconditioning.
!       Here, norm() is the Euclidean norm.  Nonzero values of ITOL
!       require  additional work  to  calculate the  actual  scaled
!       residual  or its scaled/preconditioned  form,  and/or   the
!       approximate solution XL.  Hence, these values of  ITOL will
!       not be as efficient as ITOL=0.
!
! *Cautions:
!     This routine will attempt to write to the Fortran logical output
!     unit IUNIT, if IUNIT .ne. 0.  Thus, the user must make sure that
!     this logical unit is attached to a file or terminal before calling
!     this routine with a non-zero value for IUNIT.  This routine does
!     not check for the validity of a non-zero IUNIT unit number.
!
!     This routine does not verify that ITOL has a valid value.
!     The calling routine should make such a test before calling
!     ISDGMR, as is done in DGMRES.
!
!***SEE ALSO  DGMRES
!***ROUTINES CALLED  D1MACH, DCOPY, DNRM2, DRLCAL, DSCAL, DXLCAL
!***COMMON BLOCKS    DSLBLK
!***REVISION HISTORY  (YYMMDD)
!   890404  DATE WRITTEN
!   890404  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Corrected conversion errors, etc.  (FNF)
!   910502  Removed MSOLVE from ROUTINES CALLED list.  (FNF)
!   910506  Made subsidiary to DGMRES.  (FNF)
!   920407  COMMON BLOCK renamed DSLBLK.  (WRB)
!   920511  Added complete declaration section.  (WRB)
!   921026  Corrected D to E in output format.  (FNF)
!   921113  Corrected C***CATEGORY line.  (FNF)
!***END PROLOGUE  ISDGMR
!     .. Scalar Arguments ..
      DOUBLE PRECISION BNRM, ERR, PROD, R0NRM, RNRM, SNORMW, TOL
      INTEGER ISYM, ITER, ITMAX, ITOL, IUNIT, JPRE, JSCAL, KMP, LGMR, &
              MAXL, MAXLP1, N, NELT, NMSL
!     .. Array Arguments ..
      DOUBLE PRECISION A(*), B(*), DZ(*), HES(MAXLP1, MAXL), Q(*), R(*), &
                       RWORK(*), SB(*), SX(*), V(N,*), X(*), XL(*), Z(*)
      INTEGER IA(*), IWORK(*), JA(*)
!     .. Subroutine Arguments ..
      EXTERNAL MSOLVE
!     .. Arrays in Common ..
      DOUBLE PRECISION SOLN(1)
!     .. Local Scalars ..
      DOUBLE PRECISION DXNRM, FUZZ, RAT, RATMAX, SOLNRM, TEM
      INTEGER I, IELMAX
!     .. External Functions ..
      DOUBLE PRECISION D1MACH, DNRM2
      EXTERNAL D1MACH, DNRM2
!     .. External Subroutines ..
      EXTERNAL DCOPY, DRLCAL, DSCAL, DXLCAL
!     .. Intrinsic Functions ..
      INTRINSIC ABS, MAX, SQRT
!     .. Common blocks ..
      COMMON /DSLBLK/ SOLN
!     .. Save statement ..
      SAVE SOLNRM
!***FIRST EXECUTABLE STATEMENT  ISDGMR
      ISDGMR = 0
      IF ( ITOL.EQ.0 ) THEN
!
!       Use input from DPIGMR to determine if stop conditions are met.
!
         ERR = RNRM/BNRM
      ENDIF
      IF ( (ITOL.GT.0) .AND. (ITOL.LE.3) ) THEN
!
!       Use DRLCAL to calculate the scaled residual vector.
!       Store answer in R.
!
         IF ( LGMR.NE.0 ) CALL DRLCAL(N, KMP, LGMR, MAXL, V, Q, R, &
                                      SNORMW, PROD, R0NRM)
         IF ( ITOL.LE.2 ) THEN
!         err = ||Residual||/||RightHandSide||(2-Norms).
            ERR = DNRM2(N, R, 1)/BNRM
!
!         Unscale R by R0NRM*PROD when KMP < MAXL.
!
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN
               TEM = 1.0D0/(R0NRM*PROD)
               CALL DSCAL(N, TEM, R, 1)
            ENDIF
         ELSEIF ( ITOL.EQ.3 ) THEN
!         err = Max |(Minv*Residual)(i)/x(i)|
!         When JPRE .lt. 0, R already contains Minv*Residual.
            IF ( JPRE.GT.0 ) THEN
               CALL MSOLVE(N, R, DZ, NELT, IA, JA, A, ISYM, RWORK, &
                    IWORK)
               NMSL = NMSL + 1
            ENDIF
!
!         Unscale R by R0NRM*PROD when KMP < MAXL.
!
            IF ( (KMP.LT.MAXL) .AND. (LGMR.NE.0) ) THEN
               TEM = 1.0D0/(R0NRM*PROD)
               CALL DSCAL(N, TEM, R, 1)
            ENDIF
!
            FUZZ = D1MACH(1)
            IELMAX = 1
            RATMAX = ABS(DZ(1))/MAX(ABS(X(1)),FUZZ)
            DO 25 I = 2, N
               RAT = ABS(DZ(I))/MAX(ABS(X(I)),FUZZ)
               IF( RAT.GT.RATMAX ) THEN
                  IELMAX = I
                  RATMAX = RAT
               ENDIF
 25         CONTINUE
            ERR = RATMAX
            IF( RATMAX.LE.TOL ) ISDGMR = 1
            IF( IUNIT.GT.0 ) WRITE(IUNIT,1020) ITER, IELMAX, RATMAX
            RETURN
         ENDIF
      ENDIF
      IF ( ITOL.EQ.11 ) THEN
!
!       Use DXLCAL to calculate the approximate solution XL.
!
         IF ( (LGMR.NE.0) .AND. (ITER.GT.0) ) THEN
            CALL DXLCAL(N, LGMR, X, XL, XL, HES, MAXLP1, Q, V, R0NRM, &
                 DZ, SX, JSCAL, JPRE, MSOLVE, NMSL, RWORK, IWORK, &
                 NELT, IA, JA, A, ISYM)
         ELSEIF ( ITER.EQ.0 ) THEN
!         Copy X to XL to check if initial guess is good enough.
            CALL DCOPY(N, X, 1, XL, 1)
         ELSE
!         Return since this is the first call to DPIGMR on a restart.
            RETURN
         ENDIF
!
         IF ((JSCAL .EQ. 0) .OR.(JSCAL .EQ. 2)) THEN
!         err = ||x-TrueSolution||/||TrueSolution||(2-Norms).
            IF ( ITER.EQ.0 ) SOLNRM = DNRM2(N, SOLN, 1)
            DO 30 I = 1, N
               DZ(I) = XL(I) - SOLN(I)
 30         CONTINUE
            ERR = DNRM2(N, DZ, 1)/SOLNRM
         ELSE
            IF (ITER .EQ. 0) THEN
               SOLNRM = 0
               DO 40 I = 1,N
                  SOLNRM = SOLNRM + (SX(I)*SOLN(I))**2
 40            CONTINUE
               SOLNRM = SQRT(SOLNRM)
            ENDIF
            DXNRM = 0
            DO 50 I = 1,N
               DXNRM = DXNRM + (SX(I)*(XL(I)-SOLN(I)))**2
 50         CONTINUE
            DXNRM = SQRT(DXNRM)
!         err = ||SX*(x-TrueSolution)||/||SX*TrueSolution|| (2-Norms).
            ERR = DXNRM/SOLNRM
         ENDIF
      ENDIF
!
      IF( IUNIT.NE.0 ) THEN
         IF( ITER.EQ.0 ) THEN
            WRITE(IUNIT,1000) N, ITOL, MAXL, KMP
         ENDIF
         WRITE(IUNIT,1010) ITER, RNRM/BNRM, ERR
      ENDIF
      IF ( ERR.LE.TOL ) ISDGMR = 1
!
      RETURN
 1000 FORMAT(' Generalized Minimum Residual(',I3,I3,') for ', &
           'N, ITOL = ',I5, I5, &
           /' ITER','   Natural Err Est','   Error Estimate')
 1010 FORMAT(1X,I4,1X,D16.7,1X,D16.7)
 1020 FORMAT(1X,' ITER = ',I5, ' IELMAX = ',I5, &
           ' |R(IELMAX)/X(IELMAX)| = ',D12.5)
!------------- LAST LINE OF ISDGMR FOLLOWS ----------------------------
      END
