!DECK QS2I1D
      SUBROUTINE QS2I1D (IA, JA, A, N, KFLAG)
!***BEGIN PROLOGUE  QS2I1D
!***SUBSIDIARY
!***PURPOSE  Sort an integer array, moving an integer and DP array.
!            This routine sorts the integer array IA and makes the same
!            interchanges in the integer array JA and the double pre-
!            cision array A.  The array IA may be sorted in increasing
!            order or decreasing order.  A slightly modified QUICKSORT
!            algorithm is used.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  N6A2A
!***TYPE      DOUBLE PRECISION (QS2I1R-S, QS2I1D-D)
!***KEYWORDS  SINGLETON QUICKSORT, SLAP, SORT, SORTING
!***AUTHOR  Jones, R. E., (SNLA)
!           Kahaner, D. K., (NBS)
!           Seager, M. K., (LLNL) seager@llnl.gov
!           Wisniewski, J. A., (SNLA)
!***DESCRIPTION
!     Written by Rondall E Jones
!     Modified by John A. Wisniewski to use the Singleton QUICKSORT
!     algorithm. date 18 November 1976.
!
!     Further modified by David K. Kahaner
!     National Bureau of Standards
!     August, 1981
!
!     Even further modification made to bring the code up to the
!     Fortran 77 level and make it more readable and to carry
!     along one integer array and one double precision array during
!     the sort by
!     Mark K. Seager
!     Lawrence Livermore National Laboratory
!     November, 1987
!     This routine was adapted from the ISORT routine.
!
!     ABSTRACT
!         This routine sorts an integer array IA and makes the same
!         interchanges in the integer array JA and the double precision
!         array A.
!         The array IA may be sorted in increasing order or decreasing
!         order.  A slightly modified quicksort algorithm is used.
!
!     DESCRIPTION OF PARAMETERS
!        IA - Integer array of values to be sorted.
!        JA - Integer array to be carried along.
!         A - Double Precision array to be carried along.
!         N - Number of values in integer array IA to be sorted.
!     KFLAG - Control parameter
!           = 1 means sort IA in INCREASING order.
!           =-1 means sort IA in DECREASING order.
!
!***SEE ALSO  DS2Y
!***REFERENCES  R. C. Singleton, Algorithm 347, An Efficient Algorithm
!                 for Sorting With Minimal Storage, Communications ACM
!                 12:3 (1969), pp.185-7.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   761118  DATE WRITTEN
!   890125  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   900805  Changed XERROR calls to calls to XERMSG.  (RWC)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910506  Made subsidiary to DS2Y and corrected reference.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!   920929  Corrected format of reference.  (FNF)
!   921012  Corrected all f.p. constants to double precision.  (FNF)
!***END PROLOGUE  QS2I1D
!VD$R NOVECTOR
!VD$R NOCONCUR
!     .. Scalar Arguments ..
      INTEGER KFLAG, N
!     .. Array Arguments ..
      DOUBLE PRECISION A(N)
      INTEGER IA(N), JA(N)
!     .. Local Scalars ..
      DOUBLE PRECISION R, TA, TTA
      INTEGER I, IIT, IJ, IT, J, JJT, JT, K, KK, L, M, NN
!     .. Local Arrays ..
      INTEGER IL(21), IU(21)
!     .. External Subroutines ..
      EXTERNAL XERMSG
!     .. Intrinsic Functions ..
      INTRINSIC ABS, INT
!***FIRST EXECUTABLE STATEMENT  QS2I1D
      NN = N
      IF (NN.LT.1) THEN
         CALL XERMSG ('SLATEC', 'QS2I1D', &
            'The number of values to be sorted was not positive.', 1, 1)
         RETURN
      ENDIF
      IF( N.EQ.1 ) RETURN
      KK = ABS(KFLAG)
      IF ( KK.NE.1 ) THEN
         CALL XERMSG ('SLATEC', 'QS2I1D', &
            'The sort control parameter, K, was not 1 or -1.', 2, 1)
         RETURN
      ENDIF
!
!     Alter array IA to get decreasing order if needed.
!
      IF( KFLAG.LT.1 ) THEN
         DO 20 I=1,NN
            IA(I) = -IA(I)
 20      CONTINUE
      ENDIF
!
!     Sort IA and carry JA and A along.
!     And now...Just a little black magic...
      M = 1
      I = 1
      J = NN
      R = .375D0
 210  IF( R.LE.0.5898437D0 ) THEN
         R = R + 3.90625D-2
      ELSE
         R = R-.21875D0
      ENDIF
 225  K = I
!
!     Select a central element of the array and save it in location
!     it, jt, at.
!
      IJ = I + INT ((J-I)*R)
      IT = IA(IJ)
      JT = JA(IJ)
      TA = A(IJ)
!
!     If first element of array is greater than it, interchange with it.
!
      IF( IA(I).GT.IT ) THEN
         IA(IJ) = IA(I)
         IA(I)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(I)
         JA(I)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(I)
         A(I)   = TA
         TA     = A(IJ)
      ENDIF
      L=J
!
!     If last element of array is less than it, swap with it.
!
      IF( IA(J).LT.IT ) THEN
         IA(IJ) = IA(J)
         IA(J)  = IT
         IT     = IA(IJ)
         JA(IJ) = JA(J)
         JA(J)  = JT
         JT     = JA(IJ)
         A(IJ)  = A(J)
         A(J)   = TA
         TA     = A(IJ)
!
!     If first element of array is greater than it, swap with it.
!
         IF ( IA(I).GT.IT ) THEN
            IA(IJ) = IA(I)
            IA(I)  = IT
            IT     = IA(IJ)
            JA(IJ) = JA(I)
            JA(I)  = JT
            JT     = JA(IJ)
            A(IJ)  = A(I)
            A(I)   = TA
            TA     = A(IJ)
         ENDIF
      ENDIF
!
!     Find an element in the second half of the array which is
!     smaller than it.
!
  240 L=L-1
      IF( IA(L).GT.IT ) GO TO 240
!
!     Find an element in the first half of the array which is
!     greater than it.
!
  245 K=K+1
      IF( IA(K).LT.IT ) GO TO 245
!
!     Interchange these elements.
!
      IF( K.LE.L ) THEN
         IIT   = IA(L)
         IA(L) = IA(K)
         IA(K) = IIT
         JJT   = JA(L)
         JA(L) = JA(K)
         JA(K) = JJT
         TTA   = A(L)
         A(L)  = A(K)
         A(K)  = TTA
         GOTO 240
      ENDIF
!
!     Save upper and lower subscripts of the array yet to be sorted.
!
      IF( L-I.GT.J-K ) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 260
!
!     Begin again on another portion of the unsorted array.
!
  255 M = M-1
      IF( M.EQ.0 ) GO TO 300
      I = IL(M)
      J = IU(M)
  260 IF( J-I.GE.1 ) GO TO 225
      IF( I.EQ.J ) GO TO 255
      IF( I.EQ.1 ) GO TO 210
      I = I-1
  265 I = I+1
      IF( I.EQ.J ) GO TO 255
      IT = IA(I+1)
      JT = JA(I+1)
      TA =  A(I+1)
      IF( IA(I).LE.IT ) GO TO 265
      K=I
  270 IA(K+1) = IA(K)
      JA(K+1) = JA(K)
      A(K+1)  =  A(K)
      K = K-1
      IF( IT.LT.IA(K) ) GO TO 270
      IA(K+1) = IT
      JA(K+1) = JT
      A(K+1)  = TA
      GO TO 265
!
!     Clean up, if necessary.
!
  300 IF( KFLAG.LT.1 ) THEN
         DO 310 I=1,NN
            IA(I) = -IA(I)
 310     CONTINUE
      ENDIF
      RETURN
!------------- LAST LINE OF QS2I1D FOLLOWS ----------------------------
      END
