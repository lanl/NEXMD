!DECK DHELS
      SUBROUTINE DHELS (A, LDA, N, Q, B)
!***BEGIN PROLOGUE  DHELS
!***SUBSIDIARY
!***PURPOSE  Internal routine for DGMRES.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      DOUBLE PRECISION (SHELS-S, DHELS-D)
!***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
!***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
!           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
!           Seager, Mark K., (LLNL), seager@llnl.gov
!             Lawrence Livermore National Laboratory
!             PO Box 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!***DESCRIPTION
!        This routine is extracted from the LINPACK routine SGESL with
!        changes due to the fact that A is an upper Hessenberg matrix.
!
!        DHELS solves the least squares problem:
!
!                   MIN(B-A*X,B-A*X)
!
!        using the factors computed by DHEQR.
!
! *Usage:
!      INTEGER LDA, N
!      DOUBLE PRECISION A(LDA,N), Q(2*N), B(N+1)
!
!      CALL DHELS(A, LDA, N, Q, B)
!
! *Arguments:
! A       :IN       Double Precision A(LDA,N)
!          The output from DHEQR which contains the upper
!          triangular factor R in the QR decomposition of A.
! LDA     :IN       Integer
!          The leading dimension of the array A.
! N       :IN       Integer
!          A is originally an (N+1) by N matrix.
! Q       :IN       Double Precision Q(2*N)
!          The coefficients of the N Givens rotations
!          used in the QR factorization of A.
! B       :INOUT    Double Precision B(N+1)
!          On input, B is the right hand side vector.
!          On output, B is the solution vector X.
!
!***SEE ALSO  DGMRES
!***ROUTINES CALLED  DAXPY
!***REVISION HISTORY  (YYMMDD)
!   890404  DATE WRITTEN
!   890404  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Added C***FIRST EXECUTABLE STATEMENT line.  (FNF)
!   910506  Made subsidiary to DGMRES.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!***END PROLOGUE  DHELS
!         The following is for optimized compilation on LLNL/LTSS Crays.
!LLL. OPTIMIZE
!     .. Scalar Arguments ..
      INTEGER LDA, N
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), B(*), Q(*)
!     .. Local Scalars ..
      DOUBLE PRECISION C, S, T, T1, T2
      INTEGER IQ, K, KB, KP1
!     .. External Subroutines ..
      EXTERNAL DAXPY
!***FIRST EXECUTABLE STATEMENT  DHELS
!
!         Minimize(B-A*X,B-A*X).  First form Q*B.
!
      DO 20 K = 1, N
         KP1 = K + 1
         IQ = 2*(K-1) + 1
         C = Q(IQ)
         S = Q(IQ+1)
         T1 = B(K)
         T2 = B(KP1)
         B(K) = C*T1 - S*T2
         B(KP1) = S*T1 + C*T2
 20   CONTINUE
!
!         Now solve  R*X = Q*B.
!
      DO 40 KB = 1, N
         K = N + 1 - KB
         B(K) = B(K)/A(K,K)
         T = -B(K)
         CALL DAXPY(K-1, T, A(1,K), 1, B(1), 1)
 40   CONTINUE
      RETURN
!------------- LAST LINE OF DHELS FOLLOWS ----------------------------
      END
