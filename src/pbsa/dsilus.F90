!DECK DSILUS
      SUBROUTINE DSILUS (N, NELT, IA, JA, A, ISYM, NL, IL, JL, L, DINV, &
         NU, IU, JU, U, NROW, NCOL)
!***BEGIN PROLOGUE  DSILUS
!***PURPOSE  Incomplete LU Decomposition Preconditioner SLAP Set Up.
!            Routine to generate the incomplete LDU decomposition of a
!            matrix.  The unit lower triangular factor L is stored by
!            rows and the unit upper triangular factor U is stored by
!            columns.  The inverse of the diagonal matrix D is stored.
!            No fill in is allowed.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  D2E
!***TYPE      DOUBLE PRECISION (SSILUS-S, DSILUS-D)
!***KEYWORDS  INCOMPLETE LU FACTORIZATION, ITERATIVE PRECONDITION,
!             NON-SYMMETRIC LINEAR SYSTEM, SLAP, SPARSE
!***AUTHOR  Greenbaum, Anne, (Courant Institute)
!           Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
! *Usage:
!     INTEGER N, NELT, IA(NELT), JA(NELT), ISYM
!     INTEGER NL, IL(NL), JL(NL), NU, IU(NU), JU(NU)
!     INTEGER NROW(N), NCOL(N)
!     DOUBLE PRECISION A(NELT), L(NL), DINV(N), U(NU)
!
!     CALL DSILUS( N, NELT, IA, JA, A, ISYM, NL, IL, JL, L,
!    $    DINV, NU, IU, JU, U, NROW, NCOL )
!
! *Arguments:
! N      :IN       Integer
!         Order of the Matrix.
! NELT   :IN       Integer.
!         Number of elements in arrays IA, JA, and A.
! IA     :IN       Integer IA(NELT).
! JA     :IN       Integer JA(NELT).
! A      :IN       Double Precision A(NELT).
!         These arrays should hold the matrix A in the SLAP Column
!         format.  See "Description", below.
! ISYM   :IN       Integer.
!         Flag to indicate symmetric storage format.
!         If ISYM=0, all non-zero entries of the matrix are stored.
!         If ISYM=1, the matrix is symmetric, and only the lower
!         triangle of the matrix is stored.
! NL     :OUT      Integer.
!         Number of non-zeros in the L array.
! IL     :OUT      Integer IL(NL).
! JL     :OUT      Integer JL(NL).
! L      :OUT      Double Precision L(NL).
!         IL, JL, L  contain the unit lower triangular factor of  the
!         incomplete decomposition  of some  matrix stored  in   SLAP
!         Row format.     The   Diagonal  of ones  *IS*  stored.  See
!         "DESCRIPTION", below for more details about the SLAP format.
! NU     :OUT      Integer.
!         Number of non-zeros in the U array.
! IU     :OUT      Integer IU(NU).
! JU     :OUT      Integer JU(NU).
! U      :OUT      Double Precision     U(NU).
!         IU, JU, U contain   the unit upper triangular factor of the
!         incomplete  decomposition    of some matrix  stored in SLAP
!         Column  format.   The Diagonal of ones   *IS*  stored.  See
!         "Description", below  for  more  details  about  the   SLAP
!         format.
! NROW   :WORK     Integer NROW(N).
!         NROW(I) is the number of non-zero elements in the I-th row
!         of L.
! NCOL   :WORK     Integer NCOL(N).
!         NCOL(I) is the number of non-zero elements in the I-th
!         column of U.
!
! *Description
!       IL, JL, L should contain the unit  lower triangular factor of
!       the incomplete decomposition of the A matrix  stored in SLAP
!       Row format.  IU, JU, U should contain  the unit upper factor
!       of the  incomplete decomposition of  the A matrix  stored in
!       SLAP Column format This ILU factorization can be computed by
!       the DSILUS routine. The diagonals (which are all one's) are
!       stored.
!
!       =================== S L A P Column format ==================
!
!       This routine  requires that  the matrix A  be stored in  the
!       SLAP Column format.  In this format the non-zeros are stored
!       counting down columns (except for  the diagonal entry, which
!       must appear first in each  "column")  and are stored  in the
!       double precision array A.   In other words,  for each column
!       in the matrix put the diagonal entry in  A.  Then put in the
!       other non-zero  elements going down  the column (except  the
!       diagonal) in order.   The  IA array holds the  row index for
!       each non-zero.  The JA array holds the offsets  into the IA,
!       A arrays  for  the  beginning  of each   column.   That  is,
!       IA(JA(ICOL)),  A(JA(ICOL)) points   to the beginning  of the
!       ICOL-th   column    in    IA and   A.      IA(JA(ICOL+1)-1),
!       A(JA(ICOL+1)-1) points to  the  end of the   ICOL-th column.
!       Note that we always have  JA(N+1) = NELT+1,  where N is  the
!       number of columns in  the matrix and NELT  is the number  of
!       non-zeros in the matrix.
!
!       Here is an example of the  SLAP Column  storage format for a
!       5x5 Matrix (in the A and IA arrays '|'  denotes the end of a
!       column):
!
!           5x5 Matrix      SLAP Column format for 5x5 matrix on left.
!                              1  2  3    4  5    6  7    8    9 10 11
!       |11 12  0  0 15|   A: 11 21 51 | 22 12 | 33 53 | 44 | 55 15 35
!       |21 22  0  0  0|  IA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
!       | 0  0 33  0 35|  JA:  1  4  6    8  9   12
!       | 0  0  0 44  0|
!       |51  0 53  0 55|
!
!       ==================== S L A P Row format ====================
!
!       This routine requires  that the matrix A  be  stored  in the
!       SLAP  Row format.   In this format  the non-zeros are stored
!       counting across  rows (except for the diagonal  entry, which
!       must  appear first  in each  "row")  and  are stored  in the
!       double precision  array A.  In other words, for each row  in
!       the matrix  put the diagonal  entry in A.   Then put in  the
!       other  non-zero elements  going across  the row  (except the
!       diagonal) in order.  The JA array holds the column index for
!       each non-zero.  The IA array holds the offsets  into the JA,
!       A  arrays  for  the   beginning  of  each  row.    That  is,
!       JA(IA(IROW)),A(IA(IROW)) are the first elements of the IROW-
!       th row in  JA and A,  and  JA(IA(IROW+1)-1), A(IA(IROW+1)-1)
!       are  the last elements  of the  IROW-th row.   Note  that we
!       always have  IA(N+1) = NELT+1, where N is the number of rows
!       in the matrix  and  NELT is the  number of non-zeros  in the
!       matrix.
!
!       Here is an example of the SLAP Row storage format for a  5x5
!       Matrix (in the A and JA arrays '|' denotes the end of a row):
!
!           5x5 Matrix         SLAP Row format for 5x5 matrix on left.
!                              1  2  3    4  5    6  7    8    9 10 11
!       |11 12  0  0 15|   A: 11 12 15 | 22 21 | 33 35 | 44 | 55 51 53
!       |21 22  0  0  0|  JA:  1  2  5 |  2  1 |  3  5 |  4 |  5  1  3
!       | 0  0 33  0 35|  IA:  1  4  6    8  9   12
!       | 0  0  0 44  0|
!       |51  0 53  0 55|
!
!***SEE ALSO  SILUR
!***REFERENCES  1. Gene Golub and Charles Van Loan, Matrix Computations,
!                  Johns Hopkins University Press, Baltimore, Maryland,
!                  1983.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   890404  DATE WRITTEN
!   890404  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   920511  Added complete declaration section.  (WRB)
!   920929  Corrected format of reference.  (FNF)
!   930701  Updated CATEGORY section.  (FNF, WRB)
!***END PROLOGUE  DSILUS
!     .. Scalar Arguments ..
      INTEGER ISYM, N, NELT, NL, NU
!     .. Array Arguments ..
      DOUBLE PRECISION A(NELT), DINV(N), L(NL), U(NU)
      INTEGER IA(NELT), IL(NL), IU(NU), JA(NELT), JL(NL), JU(NU), &
              NCOL(N), NROW(N)
!     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I, IBGN, ICOL, IEND, INDX, INDX1, INDX2, INDXC1, INDXC2, &
              INDXR1, INDXR2, IROW, ITEMP, J, JBGN, JEND, JTEMP, K, KC, &
              KR
!***FIRST EXECUTABLE STATEMENT  DSILUS
!
!         Count number of elements in each row of the lower triangle.
!
      DO 10 I=1,N
         NROW(I) = 0
         NCOL(I) = 0
 10   CONTINUE
!VD$R NOCONCUR
!VD$R NOVECTOR
      DO 30 ICOL = 1, N
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 20 J = JBGN, JEND
               IF( IA(J).LT.ICOL ) THEN
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
                  NROW(IA(J)) = NROW(IA(J)) + 1
                  IF( ISYM.NE.0 ) NCOL(IA(J)) = NCOL(IA(J)) + 1
               ENDIF
 20         CONTINUE
         ENDIF
 30   CONTINUE
      JU(1) = 1
      IL(1) = 1
      DO 40 ICOL = 1, N
         IL(ICOL+1) = IL(ICOL) + NROW(ICOL)
         JU(ICOL+1) = JU(ICOL) + NCOL(ICOL)
         NROW(ICOL) = IL(ICOL)
         NCOL(ICOL) = JU(ICOL)
 40   CONTINUE
!
!         Copy the matrix A into the L and U structures.
      DO 60 ICOL = 1, N
         DINV(ICOL) = A(JA(ICOL))
         JBGN = JA(ICOL)+1
         JEND = JA(ICOL+1)-1
         IF( JBGN.LE.JEND ) THEN
            DO 50 J = JBGN, JEND
               IROW = IA(J)
               IF( IROW.LT.ICOL ) THEN
!         Part of the upper triangle.
                  IU(NCOL(ICOL)) = IROW
                  U(NCOL(ICOL)) = A(J)
                  NCOL(ICOL) = NCOL(ICOL) + 1
               ELSE
!         Part of the lower triangle (stored by row).
                  JL(NROW(IROW)) = ICOL
                  L(NROW(IROW)) = A(J)
                  NROW(IROW) = NROW(IROW) + 1
                  IF( ISYM.NE.0 ) THEN
!         Symmetric...Copy lower triangle into upper triangle as well.
                     IU(NCOL(IROW)) = ICOL
                     U(NCOL(IROW)) = A(J)
                     NCOL(IROW) = NCOL(IROW) + 1
                  ENDIF
               ENDIF
 50         CONTINUE
         ENDIF
 60   CONTINUE
!
!         Sort the rows of L and the columns of U.
      DO 110 K = 2, N
         JBGN = JU(K)
         JEND = JU(K+1)-1
         IF( JBGN.LT.JEND ) THEN
            DO 80 J = JBGN, JEND-1
               DO 70 I = J+1, JEND
                  IF( IU(J).GT.IU(I) ) THEN
                     ITEMP = IU(J)
                     IU(J) = IU(I)
                     IU(I) = ITEMP
                     TEMP = U(J)
                     U(J) = U(I)
                     U(I) = TEMP
                  ENDIF
 70            CONTINUE
 80         CONTINUE
         ENDIF
         IBGN = IL(K)
         IEND = IL(K+1)-1
         IF( IBGN.LT.IEND ) THEN
            DO 100 I = IBGN, IEND-1
               DO 90 J = I+1, IEND
                  IF( JL(I).GT.JL(J) ) THEN
                     JTEMP = JU(I)
                     JU(I) = JU(J)
                     JU(J) = JTEMP
                     TEMP = L(I)
                     L(I) = L(J)
                     L(J) = TEMP
                  ENDIF
 90            CONTINUE
 100        CONTINUE
         ENDIF
 110  CONTINUE
!
!         Perform the incomplete LDU decomposition.
      DO 300 I=2,N
!
!           I-th row of L
         INDX1 = IL(I)
         INDX2 = IL(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 200
         DO 190 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 180
            INDXR1 = INDX1
            INDXR2 = INDX - 1
            INDXC1 = JU(JL(INDX))
            INDXC2 = JU(JL(INDX)+1) - 1
            IF(INDXC1 .GT. INDXC2) GO TO 180
 160        KR = JL(INDXR1)
 170        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 170
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 160
            ELSEIF(KR .EQ. KC) THEN
               L(INDX) = L(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 160
            ENDIF
 180        L(INDX) = L(INDX)/DINV(JL(INDX))
 190     CONTINUE
!
!         I-th column of U
 200     INDX1 = JU(I)
         INDX2 = JU(I+1) - 1
         IF(INDX1 .GT. INDX2) GO TO 260
         DO 250 INDX=INDX1,INDX2
            IF(INDX .EQ. INDX1) GO TO 240
            INDXC1 = INDX1
            INDXC2 = INDX - 1
            INDXR1 = IL(IU(INDX))
            INDXR2 = IL(IU(INDX)+1) - 1
            IF(INDXR1 .GT. INDXR2) GO TO 240
 210        KR = JL(INDXR1)
 220        KC = IU(INDXC1)
            IF(KR .GT. KC) THEN
               INDXC1 = INDXC1 + 1
               IF(INDXC1 .LE. INDXC2) GO TO 220
            ELSEIF(KR .LT. KC) THEN
               INDXR1 = INDXR1 + 1
               IF(INDXR1 .LE. INDXR2) GO TO 210
            ELSEIF(KR .EQ. KC) THEN
               U(INDX) = U(INDX) - L(INDXR1)*DINV(KC)*U(INDXC1)
               INDXR1 = INDXR1 + 1
               INDXC1 = INDXC1 + 1
               IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 210
            ENDIF
 240        U(INDX) = U(INDX)/DINV(IU(INDX))
 250     CONTINUE
!
!         I-th diagonal element
 260     INDXR1 = IL(I)
         INDXR2 = IL(I+1) - 1
         IF(INDXR1 .GT. INDXR2) GO TO 300
         INDXC1 = JU(I)
         INDXC2 = JU(I+1) - 1
         IF(INDXC1 .GT. INDXC2) GO TO 300
 270     KR = JL(INDXR1)
 280     KC = IU(INDXC1)
         IF(KR .GT. KC) THEN
            INDXC1 = INDXC1 + 1
            IF(INDXC1 .LE. INDXC2) GO TO 280
         ELSEIF(KR .LT. KC) THEN
            INDXR1 = INDXR1 + 1
            IF(INDXR1 .LE. INDXR2) GO TO 270
         ELSEIF(KR .EQ. KC) THEN
            DINV(I) = DINV(I) - L(INDXR1)*DINV(KC)*U(INDXC1)
            INDXR1 = INDXR1 + 1
            INDXC1 = INDXC1 + 1
            IF(INDXR1 .LE. INDXR2 .AND. INDXC1 .LE. INDXC2) GO TO 270
         ENDIF
!
 300  CONTINUE
!
!         Replace diagonal elements by their inverses.
!VD$ VECTOR
      DO 430 I=1,N
         DINV(I) = 1.0D0/DINV(I)
 430  CONTINUE
!
      RETURN
!------------- LAST LINE OF DSILUS FOLLOWS ----------------------------
      END
