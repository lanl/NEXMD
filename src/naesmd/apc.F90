!     ******* SAMPLE CALLING PROGRAM FOR SUBROUTINE APC           *******
!     ***     (MIN-COST        ASSIGNMENT PROBLEM)                       ***
!     ***                                                       ***
!     ***     THE PROGRAM IS BASED ON THE PAPER                       ***
!     ***     G. CARPANETO, S. MARTELLO, P. TOTH "ALGORITHMS   ***
!     ***       AND CODES FOR THE ASSIGNMENT PROBLEM",               ***
!     ***        ANNALS OF OPERATIONS RESEARCH 7, 1988.               ***
!     ***                                                       ***
!     ***     ALL THE SUBROUTINES ARE WRITTEN IN AMERICAN      ***
!     ***        STANDARD FORTRAN AND ARE ACCEPTED BY THE       ***
!     ***        PFORT VERIFIER.                                       ***
!     ***                                                       ***
!     ***     QUESTIONS        AND COMMENTS SHOULD BE DIRECTED        TO     ***
!     ***     SILVANO MARTELLO AND PAOLO TOTH                       ***
!     ***     D.E.I.S.,        UNIVERSITA DI BOLOGNA,                ***
!     ***     VIALE RISORGIMENTO 2,                            ***
!     ***     40136, BOLOGNA, ITALY.                           ***
!     ************************************************************
!
!      INTEGER A(260,260),F(260),Z
!
!      double precision b(260,260)
!
!      READ(5,10) N
!   10 FORMAT(20I4)
!
!!      DO 20 I=1,N
!        READ(5,10) (A(I,J),J=1,N)
!!   20 CONTINUE
!
!      DO 20 I=1,N
!        READ(5,*) (b(I,J),J=1,N)
!   20 CONTINUE
!
!      do i=1,n
!         do j=1,n
!!            a(i,j)=ifix(sngl(b(i,j)**2*1.d4))
!            a(i,j)=int((b(i,j)**2*1.d5))
!         enddo
!      enddo
!
!      DO I=1,N
!!        print*, (a(I,J),J=1,N)
!      enddo
!
!!     ----------------------
!
!      do i=1,n
!         do j=1,n
!            a(i,j)=-1*a(i,j)
!         enddo
!      enddo

!      CALL APC(N,A,F,Z)
!      WRITE(6,30) Z
!   30 FORMAT(26H  COST OF THE ASSIGNMENT =,I10/)
!      WRITE(6,40) (F(I),I=1,N)
!   40 FORMAT(11H ASSIGNMENT,20I4)
!
!      do i=1,n
!         print*,i,b(i,f(i))
!      enddo
!

!
!      STOP
!      END
subroutine APC(apc_data, N, A, F, Z)
    !
    ! SOLUTION OF THE LINEAR MIN-SUM ASSIGNMENT PROBLEM.
    !
    ! HUNGARIAN METHOD. COMPLEXITY O(N**3).
    !
    !
    ! MEANING OF THE INPUT PARAMETERS:
    ! N      = NUMBER OF ROWS AND COLUMNS OF THE COST MATRIX.
    ! A(I,J) = COST OF THE ASSIGNMENT OF ROW  I  TO COLUMN  J .
    ! ON RETURN, THE INPUT PARAMETERS ARE UNCHANGED.
    !
    ! MEANING OF THE OUTPUT PARAMETERS:
    ! F(I) = COLUMN ASSIGNED TO ROW  I .
    ! Z    = COST OF THE OPTIMAL ASSIGNMENT =
    !      = A(1,F(1)) + A(2,F(2)) + ... + A(N,F(N)) .
    !
    ! ALL THE PARAMETERS ARE INTEGERS.
    ! VECTOR  F  MUST BE DIMENSIONED AT LEAST AT  N , MATRIX  A
    ! AT LEAST AT  (N,N) . AS CURRENTLY DIMENSIONED, THE SIZE
    ! LIMITATION IS  N .LE. 260 . IN ALL THE SUBROUTINES, THE
    ! INTERNAL VARIABLES WHICH ARE PRESENTLY DIMENSIONED AT
    ! 260 MUST BE DIMENSIONED AT LEAST AT  N .
    !
    ! THE ONLY MACHINE-DEPENDENT CONSTANT USED IS  INF (DEFINED
    ! BY THE FIRST EXECUTABLE STATEMENT OF THIS SUBROUTINE). INF
    ! MUST BE SET TO A VERY LARGE INTEGER VALUE.
    !
    ! THE CODE IS BASED ON THE HUNGARIAN METHOD AS DESCRIBED BY
    ! LAWLER (COMBINATORIAL OPTIMIZATION : NETWORKS AND
    ! MATROIDS, HOLT, RINEHART AND WINSTON, NEW YORK, 1976).
    ! THE ALGORITHMIC PASCAL-LIKE DESCRIPTION OF THE CODE IS
    ! GIVEN IN G.CARPANETO, S.MARTELLO AND P.TOTH, ALGORITHMS AND
    ! CODES FOR THE ASSIGNMENT PROBLEM, ANNALS OF OPERATIONS
    ! RESEARCH 7, 1988.
    !
    ! SUBROUTINE APC DETERMINES THE INITIAL DUAL AND PARTIAL
    ! PRIMAL SOLUTIONS AND THEN SEARCHES FOR AUGMENTING PATHS
    ! UNTIL ALL ROWS AND COLUMNS ARE ASSIGNED.
    !
    ! MEANING OF THE MAIN INTERNAL VARIABLES:
    ! FB(J) = ROW ASSIGNED TO COLUMN  J .
    ! M     = NUMBER OF INITIAL ASSIGNMENTS.
    ! U(I)  = DUAL VARIABLE ASSOCIATED WITH ROW  I .
    ! V(J)  = DUAL VARIABLE ASSOCIATED WITH COLUMN  J .
    !
    ! APC NEEDS THE FOLLOWING SUBROUTINES: INCR
    !                                      INIT
    !                                      PATH
    !
    ! ALL THE SUBROUTINES ARE WRITTEN IN AMERICAN NATIONAL
    ! STANDARD FORTRAN AND ARE ACCEPTED BY THE PFORT VERIFIER.
    !
    !
    ! THIS WORK WAS SUPPORTED BY  C.N.R. , ITALY.
    !
    use naesmd_module, only : apc_common_struct
    type(apc_common_struct) :: apc_data

    integer A(260, 260), F(260), Z
    INF = 10**9
    ! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
    call INIT(apc_data, N, A, F, M, INF)
    if (M .EQ. N) GO TO 20
    ! SOLUTION OF THE REDUCED PROBLEM.
    do 10 I = 1, N
        if (F(I) .GT. 0) GO TO 10
        ! DETERMINATION OF AN AUGMENTING PATH STARTING FROM ROW  I .
        call PATH_AUG(apc_data, N, A, I, F, INF, J)
        ! ASSIGNMENT OF ROW  I  AND COLUMN  J .
        call INCR(apc_data, F, J)
10  continue
    ! COMPUTATION OF THE SOLUTION COST  Z .
20  Z = 0
    do 30 K = 1, N
        Z = Z + apc_data%U(K) + apc_data%V(K)
30  continue
    return
end subroutine APC
subroutine INCR(apc_data, F, J)
    !
    ! ASSIGNMENT OF COLUMN  J .
    !
    use naesmd_module, only : apc_common_struct
    type(apc_common_struct) :: apc_data
    integer F(260)
10  I = apc_data%RC(J)
    apc_data%FB(J) = I
    JJ = F(I)
    F(I) = J
    J = JJ
    if (J .GT. 0) GO TO 10
    return
end subroutine INCR
subroutine INIT(apc_data, N, A, F, M, INF)
    !
    ! SEARCH FOR THE INITIAL DUAL AND PARTIAL PRIMAL SOLUTIONS.
    !
    ! P(I) = FIRST UNSCANNED COLUMN OF ROW  I .
    !
    use naesmd_module, only : apc_common_struct
    type(apc_common_struct) :: apc_data
    integer A(260, 260), F(260)
    integer R
    ! PHASE 1 .
    M = 0
    do 10 K = 1, N
        F(K) = 0
        apc_data%FB(K) = 0
10  continue
    ! SCANNING OF THE COLUMNS ( INITIALIZATION OF  V(J) ).
    do 40 J = 1, N
        MIN = INF
        do 30 I = 1, N
            IA = A(I, J)
            if (IA .GT. MIN) GO TO 30
            if (IA .LT. MIN) GO TO 20
            if (F(I) .NE. 0) GO TO 30
20          MIN = IA
            R = I
30      continue
        apc_data%V(J) = MIN
        if (F(R) .NE. 0) GO TO 40
        ! ASSIGNMENT OF COLUMN  J  TO ROW  R .
        M = M + 1
        apc_data%FB(J) = R
        F(R) = J
        apc_data%U(R) = 0
        apc_data%P(R) = J + 1
40  continue
    ! PHASE 2 .
    ! SCANNING OF THE UNASSIGNED ROWS ( UPDATING OF  U(I) ).
    do 110 I = 1, N
        if (F(I) .NE. 0) GO TO 110
        MIN = INF
        do 60 K = 1, N
            IA = A(I, K) - apc_data%V(K)
            if (IA .GT. MIN) GO TO 60
            if (IA .LT. MIN) GO TO 50
            if (apc_data%FB(K) .NE. 0) GO TO 60
            if (apc_data%FB(J) .EQ. 0) GO TO 60
50          MIN = IA
            J = K
60      continue
        apc_data%U(I) = MIN
        JMIN = J
        if (apc_data%FB(J) .EQ. 0) GO TO 100
        do 80 J = JMIN, N
            if (A(I, J) - apc_data%V(J) .GT. MIN) GO TO 80
            R = apc_data%FB(J)
            KK = apc_data%P(R)
            if (KK .GT. N) GO TO 80
            do 70 K = KK, N
                if (apc_data%FB(K) .GT. 0) GO TO 70
                if (A(R, K) - apc_data%U(R) - apc_data%V(K) .EQ. 0) GO TO 90
70          continue
            apc_data%P(R) = N + 1
80      continue
        GO TO 110
        ! REASSIGNMENT OF ROW  R  AND COLUMN  K .
90      F(R) = K
        apc_data%FB(K) = R
        apc_data%P(R) = K + 1
        ! ASSIGNMENT OF COLUMN  J  TO ROW  I .
100     M = M + 1
        F(I) = J
        apc_data%FB(J) = I
        apc_data%P(I) = J + 1
110 continue
    return
end subroutine INIT
subroutine PATH_AUG(apc_data, N, A, II, F, INF, JJ)
    !
    ! DETERMINATION OF AN AUGMENTING PATH STARTING FROM
    ! UNASSIGNED ROW  II  AND TERMINATING AT UNASSIGNED COLUMN
    ! JJ , WITH UPDATING OF DUAL VARIABLES  U(I)  AND  V(J) .
    !
    ! MEANING OF THE MAIN INTERNAL VARIABLES:
    ! LR(L) = L-TH LABELLED ROW ( L=1,NLR ).
    ! PI(J) = MIN ( A(I,J) - U(I) - V(J) , SUCH THAT ROW  I  IS
    !         LABELLED AND NOT EQUAL TO  FB(J) ).
    ! RC(J) = ROW PRECEDING COLUMN  J  IN THE CURRENT
    !         ALTERNATING PATH.
    ! UC(L) = L-TH UNLABELLED COLUMN ( L=1,NUC ).
    !
    use naesmd_module, only : apc_common_struct
    type(apc_common_struct) :: apc_data
    integer A(260, 260), F(260), Z
    integer PI(260), LR(260), UC(260)
    integer R
    ! INITIALIZATION.
    LR(1) = II
    do 10 K = 1, N
        PI(K) = A(II, K) - apc_data%U(II) - apc_data%V(K)
        apc_data%RC(K) = II
        UC(K) = K
10  continue
    NUC = N
    NLR = 1
    GO TO 40
    ! SCANNING OF THE LABELLED ROWS.
20  R = LR(NLR)
    do 30 L = 1, NUC
        J = UC(L)
        IA = A(R, J) - apc_data%U(R) - apc_data%V(J)
        if (IA .GE. PI(J)) GO TO 30
        PI(J) = IA
        apc_data%RC(J) = R
30  continue
    ! SEARCH FOR A ZERO ELEMENT IN AN UNLABELLED COLUMN.
40  do 50 L = 1, NUC
        J = UC(L)
        if (PI(J) .EQ. 0) GO TO 100
50  continue
    ! UPDATING OF THE DUAL VARIABLES  U(I)  AND  V(J) .
    MIN = INF
    do 60 L = 1, NUC
        J = UC(L)
        if (MIN .GT. PI(J)) MIN = PI(J)
60  continue
    do 70 L = 1, NLR
        R = LR(L)
        apc_data%U(R) = apc_data%U(R) + MIN
70  continue
    do 90 J = 1, N
        if (PI(J) .EQ. 0) GO TO 80
        PI(J) = PI(J) - MIN
        GO TO 90
80      apc_data%V(J) = apc_data%V(J) - MIN
90  continue
    GO TO 40
100 if (apc_data%FB(J) .EQ. 0) GO TO 110
    ! LABELLING OF ROW  FB(J)  AND REMOVAL OF THE LABEL  OF
    ! COLUMN  J .
    NLR = NLR + 1
    LR(NLR) = apc_data%FB(J)
    UC(L) = UC(NUC)
    NUC = NUC - 1
    GO TO 20
    ! DETERMINATION OF THE UNASSIGNED COLUMN  J .
110 JJ = J
    return
end subroutine PATH_AUG
