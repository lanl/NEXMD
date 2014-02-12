!DECK DORTH
      SUBROUTINE DORTH (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
!***BEGIN PROLOGUE  DORTH
!***SUBSIDIARY
!***PURPOSE  Internal routine for DGMRES.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2A4, D2B4
!***TYPE      DOUBLE PRECISION (SORTH-S, DORTH-D)
!***KEYWORDS  GENERALIZED MINIMUM RESIDUAL, ITERATIVE PRECONDITION,
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
!***AUTHOR  Brown, Peter, (LLNL), pnbrown@llnl.gov
!           Hindmarsh, Alan, (LLNL), alanh@llnl.gov
!           Seager, Mark K., (LLNL), seager@llnl.gov
!             Lawrence Livermore National Laboratory
!             PO Box 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!***DESCRIPTION
!        This routine  orthogonalizes  the  vector  VNEW  against the
!        previous KMP  vectors in the   V array.  It uses  a modified
!        Gram-Schmidt   orthogonalization procedure with  conditional
!        reorthogonalization.
!
! *Usage:
!      INTEGER N, LL, LDHES, KMP
!      DOUBLE PRECISION VNEW(N), V(N,LL), HES(LDHES,LL), SNORMW
!
!      CALL DORTH(VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
!
! *Arguments:
! VNEW   :INOUT    Double Precision VNEW(N)
!         On input, the vector of length N containing a scaled
!         product of the Jacobian and the vector V(*,LL).
!         On output, the new vector orthogonal to V(*,i0) to V(*,LL),
!         where i0 = max(1, LL-KMP+1).
! V      :IN       Double Precision V(N,LL)
!         The N x LL array containing the previous LL
!         orthogonal vectors V(*,1) to V(*,LL).
! HES    :INOUT    Double Precision HES(LDHES,LL)
!         On input, an LL x LL upper Hessenberg matrix containing,
!         in HES(I,K), K.lt.LL, the scaled inner products of
!         A*V(*,K) and V(*,i).
!         On return, column LL of HES is filled in with
!         the scaled inner products of A*V(*,LL) and V(*,i).
! N      :IN       Integer
!         The order of the matrix A, and the length of VNEW.
! LL     :IN       Integer
!         The current order of the matrix HES.
! LDHES  :IN       Integer
!         The leading dimension of the HES array.
! KMP    :IN       Integer
!         The number of previous vectors the new vector VNEW
!         must be made orthogonal to (KMP .le. MAXL).
! SNORMW :OUT      DOUBLE PRECISION
!         Scalar containing the l-2 norm of VNEW.
!
!***SEE ALSO  DGMRES
!***ROUTINES CALLED  DAXPY, DDOT, DNRM2
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
!***END PROLOGUE  DORTH
!         The following is for optimized compilation on LLNL/LTSS Crays.
!LLL. OPTIMIZE
!     .. Scalar Arguments ..
      DOUBLE PRECISION SNORMW
      INTEGER KMP, LDHES, LL, N
!     .. Array Arguments ..
      DOUBLE PRECISION HES(LDHES,*), V(N,*), VNEW(*)
!     .. Local Scalars ..
      DOUBLE PRECISION ARG, SUMDSQ, TEM, VNRM
      INTEGER I, I0
!     .. External Functions ..
      DOUBLE PRECISION DDOT, DNRM2
      EXTERNAL DDOT, DNRM2
!     .. External Subroutines ..
      EXTERNAL DAXPY
!     .. Intrinsic Functions ..
      INTRINSIC MAX, SQRT
!***FIRST EXECUTABLE STATEMENT  DORTH
!
!         Get norm of unaltered VNEW for later use.
!
      VNRM = DNRM2(N, VNEW, 1)
!   -------------------------------------------------------------------
!         Perform the modified Gram-Schmidt procedure on VNEW =A*V(LL).
!         Scaled inner products give new column of HES.
!         Projections of earlier vectors are subtracted from VNEW.
!   -------------------------------------------------------------------
      I0 = MAX(1,LL-KMP+1)
      DO 10 I = I0,LL
         HES(I,LL) = DDOT(N, V(1,I), 1, VNEW, 1)
         TEM = -HES(I,LL)
         CALL DAXPY(N, TEM, V(1,I), 1, VNEW, 1)
 10   CONTINUE
!   -------------------------------------------------------------------
!         Compute SNORMW = norm of VNEW.  If VNEW is small compared
!         to its input value (in norm), then reorthogonalize VNEW to
!         V(*,1) through V(*,LL).  Correct if relative correction
!         exceeds 1000*(unit roundoff).  Finally, correct SNORMW using
!         the dot products involved.
!   -------------------------------------------------------------------
      SNORMW = DNRM2(N, VNEW, 1)
      IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0
      DO 30 I = I0,LL
         TEM = -DDOT(N, V(1,I), 1, VNEW, 1)
         IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GO TO 30
         HES(I,LL) = HES(I,LL) - TEM
         CALL DAXPY(N, TEM, V(1,I), 1, VNEW, 1)
         SUMDSQ = SUMDSQ + TEM**2
 30   CONTINUE
      IF (SUMDSQ .EQ. 0.0D0) RETURN
      ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
!
      RETURN
!------------- LAST LINE OF DORTH FOLLOWS ----------------------------
      END
