!DECK DHEQR
      SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)
!***BEGIN PROLOGUE  DHEQR
!***SUBSIDIARY
!***PURPOSE  Internal routine for DGMRES.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      DOUBLE PRECISION (SHEQR-S, DHEQR-D)
!***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
!***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
!           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
!           Seager, Mark K., (LLNL), seager@llnl.gov
!             Lawrence Livermore National Laboratory
!             PO Box 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!***DESCRIPTION
!        This   routine  performs  a QR   decomposition  of an  upper
!        Hessenberg matrix A using Givens  rotations.  There  are two
!        options  available: 1)  Performing  a fresh decomposition 2)
!        updating the QR factors by adding a row and  a column to the
!        matrix A.
!
! *Usage:
!      INTEGER LDA, N, INFO, IJOB
!      DOUBLE PRECISION A(LDA,N), Q(2*N)
!
!      CALL DHEQR(A, LDA, N, Q, INFO, IJOB)
!
! *Arguments:
! A      :INOUT    Double Precision A(LDA,N)
!         On input, the matrix to be decomposed.
!         On output, the upper triangular matrix R.
!         The factorization can be written Q*A = R, where
!         Q is a product of Givens rotations and R is upper
!         triangular.
! LDA    :IN       Integer
!         The leading dimension of the array A.
! N      :IN       Integer
!         A is an (N+1) by N Hessenberg matrix.
! Q      :OUT      Double Precision Q(2*N)
!         The factors c and s of each Givens rotation used
!         in decomposing A.
! INFO   :OUT      Integer
!         = 0  normal value.
!         = K  if  A(K,K) .eq. 0.0 .  This is not an error
!           condition for this subroutine, but it does
!           indicate that DHELS will divide by zero
!           if called.
! IJOB   :IN       Integer
!         = 1     means that a fresh decomposition of the
!                 matrix A is desired.
!         .ge. 2  means that the current decomposition of A
!                 will be updated by the addition of a row
!                 and a column.
!
!***SEE ALSO  DGMRES
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   890404  DATE WRITTEN
!   890404  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910506  Made subsidiary to DGMRES.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!***END PROLOGUE  DHEQR
!         The following is for optimized compilation on LLNL/LTSS Crays.
!LLL. OPTIMIZE
!     .. Scalar Arguments ..
      INTEGER IJOB, INFO, LDA, N
!     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*), Q(*)
!     .. Local Scalars ..
      DOUBLE PRECISION C, S, T, T1, T2
      INTEGER I, IQ, J, K, KM1, KP1, NM1
!     .. Intrinsic Functions ..
      INTRINSIC ABS, SQRT
!***FIRST EXECUTABLE STATEMENT  DHEQR
      IF (IJOB .GT. 1) GO TO 70
!   -------------------------------------------------------------------
!         A new factorization is desired.
!   -------------------------------------------------------------------
!         QR decomposition without pivoting.
!
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
!
!           Compute K-th column of R.
!           First, multiply the K-th column of A by the previous
!           K-1 Givens rotations.
!
         IF (KM1 .LT. 1) GO TO 20
         DO 10 J = 1, KM1
            I = 2*(J-1) + 1
            T1 = A(J,K)
            T2 = A(J+1,K)
            C = Q(I)
            S = Q(I+1)
            A(J,K) = C*T1 - S*T2
            A(J+1,K) = S*T1 + C*T2
 10      CONTINUE
!
!         Compute Givens components C and S.
!
 20      CONTINUE
         IQ = 2*KM1 + 1
         T1 = A(K,K)
         T2 = A(KP1,K)
         IF( T2.EQ.0.0D0 ) THEN
            C = 1
            S = 0
         ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
            T = T1/T2
            S = -1.0D0/SQRT(1.0D0+T*T)
            C = -S*T
         ELSE
            T = T2/T1
            C = 1.0D0/SQRT(1.0D0+T*T)
            S = -C*T
         ENDIF
         Q(IQ) = C
         Q(IQ+1) = S
         A(K,K) = C*T1 - S*T2
         IF( A(K,K).EQ.0.0D0 ) INFO = K
 60   CONTINUE
      RETURN
!   -------------------------------------------------------------------
!         The old factorization of a will be updated.  A row and a
!         column has been added to the matrix A.  N by N-1 is now
!         the old size of the matrix.
!   -------------------------------------------------------------------
 70   CONTINUE
      NM1 = N - 1
!   -------------------------------------------------------------------
!         Multiply the new column by the N previous Givens rotations.
!   -------------------------------------------------------------------
      DO 100 K = 1,NM1
         I = 2*(K-1) + 1
         T1 = A(K,N)
         T2 = A(K+1,N)
         C = Q(I)
         S = Q(I+1)
         A(K,N) = C*T1 - S*T2
         A(K+1,N) = S*T1 + C*T2
 100  CONTINUE
!   -------------------------------------------------------------------
!         Complete update of decomposition by forming last Givens
!         rotation, and multiplying it times the column
!         vector(A(N,N),A(NP1,N)).
!   -------------------------------------------------------------------
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF ( T2.EQ.0.0D0 ) THEN
         C = 1
         S = 0
      ELSEIF( ABS(T2).GE.ABS(T1) ) THEN
         T = T1/T2
         S = -1.0D0/SQRT(1.0D0+T*T)
         C = -S*T
      ELSE
         T = T2/T1
         C = 1.0D0/SQRT(1.0D0+T*T)
         S = -C*T
      ENDIF
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
!------------- LAST LINE OF DHEQR FOLLOWS ----------------------------
      END
