!-----------------------------------------------------------------------
!  IMSL Name:  IVPRK/DIVPRK (Single/Double precision version)
!
!  Computer:   pcdsms/DOUBLE
!
!  Revised:    October 1, 1985
!
!  Purpose:    Solve an initial-value problem for ordinary differential
!              equations using the Runge-Kutta-Verner fifth- and
!              sixth-order method.
!
!  Usage:      CALL IVPRK (IDO, NEQ, FCN, X, XEND, TOL, PARAM, Y)
!
!  Arguments:
!     IDO    - Flag indicating the state of the computation.
!              (Input/Output)
!              1        Initial entry
!              2        Normal reentry
!              3        Final call to release workspace
!              4        Return because of interrupt 1
!              5        Return because of interrupt 2 with step accepted
!              6        Return because of interrupt 2 with step rejected
!              Normally, the initial call is made with IDO=1.  The
!              routine then sets IDO=2 and this value is then used for
!              all but the last call which is made with IDO=3.  This
!              final call is only used to release workspace, which was
!              automatically allocated by the initial call with IDO=1.
!     NEQ    - Number of differential equations.  (Input)
!     FCN    - User-supplied SUBROUTINE to evaluate functions.
!              The usage is
!              CALL FCN (NEQ, X, Y, YPRIME), where
!              NEQ    - Number of equations.  (Input)
!              X      - Independent variable.  (Input)
!              Y      - Array of length NEQ containing the dependent
!                       variable values.  (Input)
!              YPRIME - Array of length NEQ containing the values of
!                       dY/dX at (X,Y).  (Output)
!              FCN must be declared EXTERNAL in the calling program.
!     X      - Independent variable.  (Input/Output)
!              On input, X supplies the initial value.
!              On output, X is replaced by XEND unless error conditions
!              arise.  See IDO for details.
!     XEND   - Value of X at which the solution is desired.  (Input)
!              XEND may be less than the initial value of X.
!     TOL    - Tolerance for error control.  (Input)
!              An attempt is made to control the norm of the local error
!              such that the global error is proportional to TOL.
!              More than one run, with different values of TOL, can be
!              used to estimate the global error.
!              Generally, it should not be greater than 0.001.
!     PARAM  - Real vector of length 50 containing optional parameters.
!              (Input/Output)
!              If a parameter is zero then a default value is used.
!              The following parameters must be set by the user.
!                 PARAM              Meaning
!                   1   HINIT  - Initial value of the step-size H.
!                                Default: See Algorithm section.
!                   2   HMIN   - Minimum value of the step-size H.
!                                Default: 0.0
!                   3   HMAX   - Maximum value of the step-size H.
!                                Default: No limit
!                   4   MXSTEP - Maximum number of steps allowed
!                                Default: 500
!                   5   MXFCN  - Maximum number of function evaluations
!                                allowed.
!                                Default: No limit
!                   6          - Not used.
!                   7   INTRP1 - If nonzero then return with IDO=4,
!                                before each step.
!                                See Remark 3.
!                                Default: 0.
!                   8   INTRP2 - If nonzero then return with IDO=5,
!                                after every successful step and with
!                                IDO=6 after every unsuccessful step.
!                                See Remark 3.
!                                Default: 0.
!                   9   SCALE  - A measure of the scale of the problem,
!                                such as an approximation to the average
!                                value of a norm of the Jacobian along
!                                the trajectory.
!                                Default: 1.0
!                  10   INORM  - Switch determining error norm.
!                                In the following Ei is the absolute
!                                value of an estimate of the error in
!                                Y(i), called Yi here.
!                                0 - min(absolute error, relative error)
!                                    = max(Ei/Wi), i=1,2,...,NEQ, where
!                                    Wi = max(abs(Yi,1.0),
!                                1 - absolute error = max(Ei), i=1,2,...
!                                2 - max(Ei/Wi), i=1,2,..., where
!                                    Wi = max(abs(Yi),FLOOR),
!                                    and FLOOR is PARAM(11).
!                                3 - Euclidian norm scaled by YMAX
!                                    = sqrt(sum(Ei**2/Wi**2)), where
!                                    Wi = max(abs(Yi),1.0); for YMAX,
!                                    see Remark 1.
!                  11   FLOOR  - Used in the norm computation.
!                                Default: 1.0
!                  12-30       - Not used.
!              The following entries in PARAM are set by the program.
!                  31   HTRIAL - Current trial step size.
!                  32   HMIN   - Minimum step size allowed.
!                  33   HMAX   - Maximum step size allowed.
!                  34   NSTEP  - Number of steps taken.
!                  35   NFCN   - Number of function evaluations used.
!                  36-50       - Not used.
!     Y      - Vector of length NEQ of dependent variables.
!              (Input/Output)
!              On input, Y contains the initial values.  On output,
!              Y contains the approximate solution.
!
!  Remarks:
!  1. Automatic workspace usage is
!              IVPRK    10*NEQ units, or
!              DIVPRK   20*NEQ units.
!     Workspace may be explicitly provided, if desired, by use of
!     I2PRK/DI2PRK.  The reference is
!              CALL I2PRK (IDO, NEQ, FCN, X, XEND, TOL, PARAM, Y,
!                          VNORM, WK)
!     The additional arguments are as follows:
!     VNORM  - User-supplied SUBROUTINE to compute the norm of the
!              error.  (Input)
!              The routine may be provided by the user, or the IMSL
!              routine I3PRK/DI3PRK may be used.
!              The usage is
!              CALL VNORM (NEQ, V, Y, YMAX, ENORM), where
!              NEQ    - Number of equations.  (Input)
!              V      - Vector of length NEQ containing the vector whose
!                       norm is to be computed.  (Input)
!              Y      - Vector of length NEQ containing the values of
!                       the dependent variable.  (Input)
!              YMAX   - Vector of length NEQ containing the maximum Y
!                       values computed so far.  (Input)
!              ENORM  - Norm of the vector V.  (Output)
!              VNORM must be declared EXTERNAL in the calling program.
!     WK     - Real work array of length 10*NEQ.  WK must not be changed
!              from the first call with IDO=1 until after the final call
!              with IDO=3.
!
!  2. Informational errors
!     Type Code
!       4   1  Cannot satisfy error condition.  TOL may be too small.
!       4   2  Too many function evaluations needed.
!       4   3  Too many steps needed.  The problem may be stiff.
!
!  3. If PARAM(7) is nonzero, the subroutine returns with
!     IDO = 4, and will resume calculation at the point of interruption
!     if reentered with IDO = 4.  If PARAM(8) is nonzero, the
!     subroutine will interrupt the calculations immediately after it
!     decides whether or not to accept the result of the most
!     recent trial step.  IDO = 5 if the routine plans to accept,
!     or IDO = 6 if it plans to reject.  IDO may be changed by the user
!     in order to force acceptance of a step (by changing IDO from 6
!     to 5) that would otherwise be rejected, or vice versa.
!     Relevant parameters to observe after return from an interrupt
!     are IDO, HTRIAL, NSTEP, NFCN, and Y.  Y is the newly computed
!     trial value, accepted or not.
!
!  GAMS:       I1a1a
!
!  Chapter:    MATH/LIBRARY Differential Equations
!
!  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
!
!  Warranty:   IMSL warrants only that IMSL testing has been applied
!              to this code.  No other warranty, expressed or implied,
!              is applicable.
!
!-----------------------------------------------------------------------
!
SUBROUTINE DIVPRK (IDO, NEQ, FCN, X, XEND, TOL, PARAM, Y)
    !                                  SPECIFICATIONS FOR ARGUMENTS
    INTEGER    IDO, NEQ
    DOUBLE PRECISION X, XEND, TOL, PARAM(*), Y(*)
    EXTERNAL   FCN
    !                                  SPECIFICATIONS FOR SAVE VARIABLES
    INTEGER    IRW
    SAVE       IRW
    !                                  SPECIFICATIONS FOR SPECIAL CASES
    !                                  SPECIFICATIONS FOR COMMON /WORKSP/
    REAL       RWKSP(5000)
    DOUBLE PRECISION RDWKSP(2500)
    DOUBLE PRECISION DWKSP(2500)
    COMPLEX    CWKSP(2500)
    COMPLEX    *16 CZWKSP(1250)
    COMPLEX    *16 ZWKSP(1250)
    INTEGER    IWKSP(5000)
    LOGICAL    LWKSP(5000)
    EQUIVALENCE (DWKSP(1), RWKSP(1))
    EQUIVALENCE (CWKSP(1), RWKSP(1)), (ZWKSP(1), RWKSP(1))
    EQUIVALENCE (IWKSP(1), RWKSP(1)), (LWKSP(1), RWKSP(1))
    EQUIVALENCE (RDWKSP(1), RWKSP(1)), (CZWKSP(1), RWKSP(1))
    COMMON     /WORKSP/ RWKSP
    !                                  SPECIFICATIONS FOR SUBROUTINES
    EXTERNAL   E1MES, E1POP, E1PSH, E1STI, I1KNR, I1KRL, DI2PRK
    !                                  SPECIFICATIONS FOR FUNCTIONS
    EXTERNAL   I1KGT, N1RTY, DI3PRK
    INTEGER    I1KGT, N1RTY
    !
    WRITE(6,*)'DIVPRK CALLED IDO=',IDO
    CALL E1PSH ('DIVPRK ')
    !                                  Check IDO
    IF (IDO.LT.1 .OR. IDO.GT.6) THEN
        CALL E1STI (1, IDO)
        CALL E1MES (5, 2, 'The status indicator IDO = %(I1).  It '// &
            'must be in the range 1 to 6.')
        GO TO 9000
    END IF
    !
    IF (IDO .EQ. 1) THEN
        !                                  Check NEQ
        IF (NEQ .LT. 1) THEN
            CALL E1STI (1, NEQ)
            CALL E1MES (5, 1, 'The number of equations NEQ = %(I1).  '// &
                'It must be at least 1.')
            GO TO 9000
        END IF
        !                                  Allocate workspace
        IRW = I1KGT(10*NEQ,4)
        IF (N1RTY(0) .NE. 0) THEN
            CALL E1MES (5, 8, ' ')
            CALL E1STI (1, NEQ)
            CALL E1MES (5, 8, 'Workspace allocation based on '// &
                'NEQ = %(I1).')
            GO TO 9000
        END IF
        CALL I1KNR
    END IF
    !
    CALL DI2PRK (IDO, NEQ, FCN, X, XEND, TOL, PARAM, Y, DI3PRK, &
        RDWKSP(IRW))
!
9000 CONTINUE
     CALL E1POP ('DIVPRK ')
     IF (IDO .EQ. 3) CALL I1KRL (1)
     RETURN
 END

 !-----------------------------------------------------------------------
 !  IMSL Name:  I2PRK/DI2PRK (Single/Double precision version)
 !
 !  Computer:   pcdsms/DOUBLE
 !
 !  Revised:    January 29, 1985
 !
 !  Purpose:    Solve an initial value problem for ordinary differential
 !              equations using the Runge-Kutta-Verner fifth- and
 !              sixth-order method.
 !
 !  Usage:      CALL I2PRK (IDO, NEQ, FCN, X, XEND, TOL, PARAM, Y,
 !                          VNORM, WK)
 !
 !  Arguments:  (See IVPRK)
 !
 !  Chapter:    MATH/LIBRARY Integration and Differentiation
 !
 !  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
 !
 !  Warranty:   IMSL warrants only that IMSL testing has been applied
 !              to this code.  No other warranty, expressed or implied,
 !              is applicable.
 !
 !-----------------------------------------------------------------------
 !
 SUBROUTINE DI2PRK (IDO, NEQ, FCN, X, XEND, TOL, PARAM, Y, VNORM, &
     WK)
     !                                  SPECIFICATIONS FOR ARGUMENTS
     INTEGER    IDO, NEQ
     DOUBLE PRECISION X, XEND, TOL, PARAM(*), Y(*), WK(NEQ,*)
     EXTERNAL   FCN, VNORM
     !                                  SPECIFICATIONS FOR LOCAL VARIABLES
     INTEGER    I, K, NEWFCN
     DOUBLE PRECISION TEMP
     !                                  SPECIFICATIONS FOR SAVE VARIABLES
     INTEGER    INTRP1, INTRP2, IXEND, MXFCN, MXSTEP, NFCN, NSTEP, &
         NSUCFL, NSUCST
     DOUBLE PRECISION EPS, EST, HINIT, HMAG, HMAX, HMIN, HTRIAL, &
         RK(43), SCALE, TINY, WNORMY, XENDPV, XTRIAL
     LOGICAL    INIT
     SAVE
     !                                  SPECIFICATIONS FOR COMMON /DI5PRK/
     COMMON     /DI5PRK/ FLOOR, INORM
     INTEGER    INORM
     DOUBLE PRECISION FLOOR
     !                                  SPECIFICATIONS FOR INTRINSICS
     !     INTRINSIC  DABS,DMAX1,DMIN1,DSIGN
     INTRINSIC  DABS, DMAX1, DMIN1, DSIGN
     DOUBLE PRECISION DABS, DMAX1, DMIN1, DSIGN
     !                                  SPECIFICATIONS FOR SUBROUTINES
     EXTERNAL   E1MES, E1STI, E1STD, E1USR, DAXPY, DCOPY, DI4PRK, &
         DSWAP
     !                                  SPECIFICATIONS FOR FUNCTIONS
     EXTERNAL   DMACH, N1RTY, DDOT
     INTEGER    N1RTY
     DOUBLE PRECISION DMACH, DDOT
     !
     DATA INIT/.FALSE./, ISTAT/0/
     !
     CALL E1PSH ('DI2PRK ')
     !
     IF (IDO.LT.1 .OR. IDO.GT.6) THEN
         CALL E1STI (1, IDO)
         CALL E1MES (5, 2, 'The status indicator IDO = %(I1).  It '// &
             'must be in the range 1 to 6.')
         GO TO 9000
     END IF
     !
     IF (IDO .EQ. 1) THEN
         IF (ISTAT .EQ. 0) THEN
             ISTAT = 1
         ELSE
             CALL E1MES (5, 8, 'IDO = 1.  IDO can only be set to 1 '// &
                 'in the initial call to the routine, or '// &
                 'if the previous call was made with IDO = 3.')
             GO TO 9000
         END IF
     ELSE IF (IDO.EQ.2 .OR. IDO.EQ.3) THEN
         IF (ISTAT .EQ. 1) THEN
             IF (IDO .EQ. 3) THEN
                 ISTAT = 0
                 GO TO 9000
             END IF
         ELSE
             CALL E1STI (1, IDO)
             CALL E1MES (5, 9, 'IDO has been set to %(I1), but the '// &
                 'routine has not been initialized in a call '// &
                 'with IDO = 1.')
             GO TO 9000
         END IF
     END IF
     !                                  Cases - Initial entry, normal
     !                                  re-entry, interrupt re-entries
     GO TO (10, 40, 9000, 70, 160, 160), IDO
 !                                  Case 1 - initial entry
10 CONTINUE
   !                                  Check NEQ
   IF (NEQ .LT. 1) THEN
       CALL E1STI (1, NEQ)
       CALL E1MES (5, 1, 'The number of equations NEQ = %(I1).  '// &
           'It must be at least 1.')
   END IF
   !                                  Check TOL
   IF (TOL .LE. 0.0D0) THEN
       CALL E1STD (1, TOL)
       CALL E1MES (5, 3, 'The tolerance TOL = %(D1).  It must '// &
           'be greater than 0.0.')
   END IF
   IF (N1RTY(0) .GT. 0) GO TO 9000
   !                                  Check INORM
   IF (NINT(PARAM(10)) .GT. 3) THEN
       CALL E1STD (1, PARAM(10))
       CALL E1MES (5, 10, 'The argument PARAM(10) = INORM = %(D1).'// &
           '  It must be less than four.')
   END IF
   !                                  Get constants
   IF (.NOT.INIT) THEN
       CALL DI4PRK (RK)
       TINY = DMACH(1)
       EPS = DMACH(4)
       INIT = .TRUE.
   END IF
   !                                  Initialize
   INORM = 0
   HMIN = 0.0D0
   HINIT = 0.0D0
   SCALE = 1.0D0
   HMAG = 0.0D0
   HMAX = 2.0D0
   MXSTEP = 500
   NSTEP = 0
   MXFCN = 0
   NFCN = 0
   INTRP1 = 0
   INTRP2 = 0
   INORM = 0
   FLOOR = 1.0D0
   DO 20  K=1, 11
       IF (PARAM(K) .LT. 0.0D0) THEN
           CALL E1STI (1, K)
           CALL E1STD (1, PARAM(K))
           CALL E1MES (5, 4, 'Argument PARAM(%(I1)) = %(D1).  It '// &
               'cannot be negative.')
       END IF
20 CONTINUE
   IF (N1RTY(0) .NE. 0) GO TO 9000
   !                                  Initial WK(*,10) = YMAX(*)
   DO 30  I=1, NEQ
       WK(I,10) = DABS(Y(I))
30 CONTINUE
   !
   IF (PARAM(1) .GT. 0.0D0) HINIT = PARAM(1)
   IF (PARAM(2) .GT. 0.0D0) HMIN = PARAM(2)
   IF (PARAM(3) .GT. 0.0D0) HMAX = PARAM(3)
   IF (PARAM(4) .GT. 0.0D0) MXSTEP = PARAM(4)
   IF (PARAM(5) .GT. 0.0D0) MXFCN = PARAM(5)
   IF (PARAM(7) .GT. 0.0D0) INTRP1 = PARAM(7)
   IF (PARAM(8) .GT. 0.0D0) INTRP2 = PARAM(8)
   IF (PARAM(9) .GT. 0.0D0) THEN
       SCALE = PARAM(9)
       HMAX = DMIN1(HMAX,2.0D0/SCALE)
   END IF
   IF (PARAM(10) .GT. 0.0D0) INORM = PARAM(10)
   IF (PARAM(11) .GT. 0.0D0) FLOOR = PARAM(11)
   !                                  Initialize EPS and HMAG
   EST = 0.0D0
   IF (HINIT .EQ. 0.0D0) THEN
       HMAG = 0.5D0*HMAX*TOL**(1.0D0/6.0D0)
   ELSE
       HMAG = 0.5D0*HINIT
   END IF
   !                                  Set previous XEND initially to
   !                                  initial value of X
   XENDPV = X
   IXEND = 0
   NSUCST = 0
   NSUCFL = 0
   GO TO 50
   !                                  Case 2 - Normal re-entry (IDO.EQ.2)
   !                                  abort if XEND reached, and either
   !                                  X changed or XEND not changed
40 CONTINUE
   IF (IXEND.NE.0 .AND. X.NE.XENDPV) THEN
       CALL E1STD (1, X)
       CALL E1MES (5, 5, 'Independent variable X = %(D1) is '// &
           'changed from the previous call.')
       GO TO 9000
   ELSE IF (IXEND.NE.0 .AND. XEND.EQ.XENDPV) THEN
       CALL E1STD (1, XEND)
       CALL E1MES (5, 6, 'Final point XEND = %(D1) is not '// &
           'changed from the previous call.')
       GO TO 9000
   END IF
   !                                  Re-initialize flag IXEND
   IXEND = 0
   !                                  Case 3 - Re-entry following an
   !                                  interrupt (IDO .EQ. 4 to 6)
50 CONTINUE
   !                                  Loop through the following four
   !                                  stages, once for each trial step
   !                                  until the occurrence of one of the
   !                                  following (A) Normal return (with
   !                                  IDO .EQ. 2) on reaching XEND in
   !                                  stage 4, or (B) An error return in
   !                                  stage 1 or 4, or (C) An interrupt
   !                                  return (with IDO .EQ. 4, 5 or 6)
   !                                  in stage 1 or 4.
60 CONTINUE
   !                                  Stage 1 -- Prepare Do calculations
   !                                  of HMIN, HMAX, etc., and some
   !                                  parameter checking. End up with
   !                                  suitable values of HMAG, XTRIAL
   !                                  and HTRIAL for use in the
   !                                  integration step. Check MXFCN
   NEWFCN = 7
   IF (IDO .NE. 6) NEWFCN = NEWFCN + 1
   IF (MXFCN.GT.0 .AND. NFCN+NEWFCN.GT.MXFCN) THEN
       CALL E1STI (1, NFCN+NEWFCN)
       CALL E1STI (2, MXFCN)
       CALL E1MES (4, 2, 'Completion of the next step would make '// &
           'the number of function evaluations '// &
           '%(I1), but only %(I2) evaluations are allowed.')
       GO TO 9000
   END IF
   !                                  Check MXSTEP
   IF (MXSTEP.GT.0 .AND. NSTEP.GE.MXSTEP) THEN
       CALL E1STI (1, MXSTEP)
       CALL E1MES (4, 3, 'Maximum number of steps allowed, '// &
           '%(I1), used.  The problem may be stiff.')
       GO TO 9000
   END IF
   !                                  Calculate slope
   IF (IDO .NE. 6) THEN
       CALL E1USR ('ON')
       CALL FCN (NEQ, X, Y, WK(1,1))
       CALL E1USR ('OFF')
       NFCN = NFCN + 1
       IF (N1RTY(0) .GT. 3) GO TO 9000
   END IF
   !                                  Calculate HMIN - Use default unless
   !                                    value prescribed
   IF (HMIN .EQ. 0.0D0) THEN
       !                                  Calculate default value of HMIN
       HMIN = 10.0D0*DMAX1(TINY,EPS*DMAX1(DABS(XEND),DABS(X)))
   END IF
   !                                  Error return if HMIN .GT. HMAX
   IF (HMIN .GT. HMAX) THEN
       CALL E1STD (1, HMIN)
       CALL E1STD (2, HMAX)
       CALL E1MES (5, 7, 'HMIN = %(D1) is greater than HMAX = '// &
           '%(D2).')
       GO TO 9000
   END IF
   !                                  Calculate preliminary HMAG -
   IF (NSUCFL .LE. 1) THEN
       !                                  After a successful step, or at most
       !                                  one failure
       TEMP = 2.0D0*HMAG
       !                                  Avoid possible overflow
       !
       IF (TOL .LT. (2.0D0/.9D0)**6*EST) TEMP = 0.9D0* &
           (TOL/EST)**(1.0D0/6.0D0)*HMAG
       !                                  Avoid reduction by more than half
       HMAG = DMAX1(TEMP,0.5D0*HMAG)
   ELSE
       !                                  After two or more successive
       !                                  failures
       HMAG = 0.5D0*HMAG
   END IF
   !                                  Check against HMAX and HMIN
   HMAG = DMAX1(DMIN1(HMAG,HMAX),HMIN)
   !                                  Interrupt 1 (with IDO=4) if
   !                                  requested
   IF (INTRP1 .NE. 0) THEN
       IDO = 4
       GO TO 9000
   END IF
   !                                  Resume here on re-entry
70 CONTINUE
   !                                  Calculate HMAG, XTRIAL - depending
   !                                  on preliminary HMAG, XEND
   IF (HMAG .LT. DABS(XEND-X)) THEN
       !                                  Do not step more than half way to
       !                                  XEND
       HMAG = DMIN1(HMAG,0.5D0*DABS(XEND-X))
       XTRIAL = X + DSIGN(HMAG,XEND-X)
   ELSE
       !                                  Hit XEND exactly
       HMAG = DABS(XEND-X)
       XTRIAL = XEND
   END IF
   !                                  Calculate HTRIAL
   HTRIAL = XTRIAL - X
   !                                  Stage 2 -- Calculate YTRIAL WK(*,2),
   !                                  ..., WK(*,8) hold intermediate
   !                                  results needed in stage 3. WK(*,9)
   !                                  is temporary storage until finally
   !                                  it holds YTRIAL.
   CALL DCOPY (NEQ, Y, 1, WK(1,9), 1)
   CALL DAXPY (NEQ, HTRIAL*RK(1), WK(1,1), 1, WK(1,9), 1)
   CALL E1USR ('ON')
   CALL FCN (NEQ, X+HTRIAL/6.0D0, WK(1,9), WK(1,2))
   CALL E1USR ('OFF')
   NFCN = NFCN + 1
   IF (N1RTY(0) .GT. 3) GO TO 9000
   !
   DO 80  K=1, NEQ
       WK(K,9) = Y(K) + HTRIAL*DDOT(2,WK(K,1),NEQ,RK(2),1)
80 CONTINUE
   CALL E1USR ('ON')
   CALL FCN (NEQ, X+4.0D0*HTRIAL/15.0D0, WK(1,9), WK(1,3))
   CALL E1USR ('OFF')
   NFCN = NFCN + 1
   IF (N1RTY(0) .GT. 3) GO TO 9000
   !
   DO 90  K=1, NEQ
       WK(K,9) = Y(K) + HTRIAL*DDOT(3,WK(K,1),NEQ,RK(4),1)
90 CONTINUE
   CALL E1USR ('ON')
   CALL FCN (NEQ, X+2.0D0*HTRIAL/3.0D0, WK(1,9), WK(1,4))
   CALL E1USR ('OFF')
   NFCN = NFCN + 1
   IF (N1RTY(0) .GT. 3) GO TO 9000
   !
   DO 100  K=1, NEQ
       WK(K,9) = Y(K) + HTRIAL*DDOT(4,WK(K,1),NEQ,RK(7),1)
100 CONTINUE
    CALL E1USR ('ON')
    CALL FCN (NEQ, X+5.0D0*HTRIAL/6.0D0, WK(1,9), WK(1,5))
    CALL E1USR ('OFF')
    NFCN = NFCN + 1
    IF (N1RTY(0) .GT. 3) GO TO 9000
    !
    DO 110  K=1, NEQ
        WK(K,9) = Y(K) + HTRIAL*DDOT(5,WK(K,1),NEQ,RK(11),1)
110 CONTINUE
    CALL E1USR ('ON')
    CALL FCN (NEQ, X+HTRIAL, WK(1,9), WK(1,6))
    CALL E1USR ('OFF')
    NFCN = NFCN + 1
    IF (N1RTY(0) .GT. 3) GO TO 9000
    !
    DO 120  K=1, NEQ
        WK(K,9) = Y(K) + HTRIAL*DDOT(5,WK(K,1),NEQ,RK(16),1)
120 CONTINUE
    CALL E1USR ('ON')
    CALL FCN (NEQ, X+HTRIAL/15.0D0, WK(1,9), WK(1,7))
    CALL E1USR ('OFF')
    NFCN = NFCN + 1
    IF (N1RTY(0) .GT. 3) GO TO 9000
    !
    DO 130  K=1, NEQ
        WK(K,9) = Y(K) + HTRIAL*DDOT(7,WK(K,1),NEQ,RK(21),1)
130 CONTINUE
    CALL E1USR ('ON')
    CALL FCN (NEQ, X+HTRIAL, WK(1,9), WK(1,8))
    CALL E1USR ('OFF')
    NFCN = NFCN + 1
    IF (N1RTY(0) .GT. 3) GO TO 9000
    !                                  Calculate YTRIAL, the extrapolated
    !                                  approximation and store in WK(*,9)
    DO 140  K=1, NEQ
        WK(K,9) = Y(K) + HTRIAL*DDOT(8,WK(K,1),NEQ,RK(28),1)
140 CONTINUE
    !                                  Stage 3 -- Calculate the error
    !                                  estimate EST Calculate the
    !                                  unweighted absolute error estimate
    !                                  vector
    DO 150  K=1, NEQ
        WK(K,2) = DDOT(8,WK(K,1),NEQ,RK(36),1)
150 CONTINUE
    !                                  Calculate the weighted max norm of
    !                                  WK(*,2) as specified by the error
    !                                  control indicator INORM
    CALL E1USR ('ON')
    CALL VNORM (NEQ, WK(1,2), Y, WK(1,10), TEMP)
    CALL E1USR ('OFF')
    IF (N1RTY(0) .GT. 3) GO TO 9000
    !                                  Calculate EST - (The weighted max
    !                                  norm of WK(*,2))*HMAG*SCALE - EST
    !                                  is intended to be a measure of the
    !                                  error per unit step in YTRIAL
    EST = TEMP*HMAG*SCALE
    !                                  Stage 4 -- Make decisions Set IDO=5
    !                                  if step acceptable, else set IDO=6
    IF (EST .LE. TOL) THEN
        IDO = 5
    ELSE
        IDO = 6
    END IF
    !                                  Interrupt 2 if requested
    IF (INTRP2 .NE. 0) THEN
        CALL DSWAP (NEQ, Y, 1, WK(1,9), 1)
        GO TO 9000
    END IF
!                                  Resume here on re-entry
160 CONTINUE
    IF (INTRP2 .NE. 0) CALL DSWAP (NEQ, Y, 1, WK(1,9), 1)
    !
    IF (IDO .EQ. 5) THEN
        !                                  Step accepted, so update X, Y from
        !                                  XTRIAL, YTRIAL
        X = XTRIAL
        CALL DCOPY (NEQ, WK(1,9), 1, Y, 1)
        !                                  Update YMAX values
        DO 170  I=1, NEQ
            WK(I,10) = DMAX1(WK(I,10),DABS(Y(I)))
170     CONTINUE
        NSUCST = NSUCST + 1
        NSUCFL = 0
        NSTEP = NSTEP + 1
        !                                  Return with IDO=2, XEND saved, flag
        !                                  set
        IF (X .EQ. XEND) THEN
            IDO = 2
            XENDPV = XEND
            IXEND = 1
            GO TO 9000
        END IF
    ELSE IF (IDO .EQ. 6) THEN
        !                                  Step not accepted - Add 1 to number
        !                                  successive failures
        NSUCFL = NSUCFL + 1
        !                                  Error return if HMAG .LE. HMIN
        IF (HMAG .LE. HMIN) THEN
            CALL E1STD (1, TOL)
            CALL E1MES (4, 1, 'Unable to satisfy the error '// &
                'requirement.  TOL = %(D1) may be too small.')
            GO TO 9000
        END IF
    END IF
    !                                  End stage 4
    GO TO 60
!                                  End loop
9000 CONTINUE
     !                                  Set PARAMs for user
     PARAM(31) = HTRIAL
     PARAM(32) = HMIN
     PARAM(33) = HMAX
     PARAM(34) = NSTEP
     PARAM(35) = NFCN
     !
     CALL E1POP ('DI2PRK ')
     RETURN
 END
 !-----------------------------------------------------------------------
 !  IMSL Name:  I3PRK/DI3PRK (Single/Double precision version)
 !
 !  Computer:   pcdsms/DOUBLE
 !
 !  Revised:    January 29, 1985
 !
 !  Purpose:    Solve an initial value problem for ordinary differential
 !              equations using the Runge-Kutta-Verner fifth and sixth
 !              order method.
 !
 !  Usage:      CALL I3PRK (NEQ, V, Y, YMAX, ENORM)
 !
 !  Arguments:
 !     NEQ    - Number of equations.  (Input)
 !     V      - Vector of length NEQ containing the vector whose
 !              norm is to be computed.  (Input)
 !     Y      - Vector of length NEQ containing the values of the
 !              dependent variable.  (Input)
 !     YMAX   - Vector of length NEQ containing the maximum Y values
 !              computed so far.  (Input)
 !     ENORM  - Norm of the vector V.  (Output)
 !
 !  Chapter:    MATH/LIBRARY Integration and Differentiation
 !
 !  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
 !
 !  Warranty:   IMSL warrants only that IMSL testing has been applied
 !              to this code.  No other warranty, expressed or implied,
 !              is applicable.
 !
 !-----------------------------------------------------------------------
 !
 SUBROUTINE DI3PRK (NEQ, V, Y, YMAX, ENORM)
     !                                  SPECIFICATIONS FOR ARGUMENTS
     INTEGER    NEQ
     DOUBLE PRECISION ENORM, V(*), Y(*), YMAX(*)
     !                                  SPECIFICATIONS FOR LOCAL VARIABLES
     INTEGER    K
     DOUBLE PRECISION WEIGHT
     !                                  SPECIFICATIONS FOR COMMON /DI5PRK/
     COMMON     /DI5PRK/ FLOOR, INORM
     INTEGER    INORM
     DOUBLE PRECISION FLOOR
     !                                  SPECIFICATIONS FOR INTRINSICS
     !     INTRINSIC  DABS,DMAX1
     INTRINSIC  DABS, DMAX1
     DOUBLE PRECISION DABS, DMAX1
     !                                  SPECIFICATIONS FOR FUNCTIONS
     EXTERNAL   IDAMAX
     INTEGER    IDAMAX
     !
     IF (INORM .EQ. 0) THEN
         !                                  max (absolute, relative)
         ENORM = 0.0D0
         DO 10  K=1, NEQ
             WEIGHT = DMAX1(1.0D0,DABS(Y(K)))
             ENORM = DMAX1(ENORM,DABS(V(K))/WEIGHT)
10       CONTINUE
     !                                  Absolute error control
     ELSE IF (INORM .EQ. 1) THEN
         ENORM = DABS(V(IDAMAX(NEQ,V,1)))
     !                                  Relative error control
     ELSE IF (INORM .EQ. 2) THEN
         ENORM = 0.0D0
         DO 20  K=1, NEQ
             WEIGHT = DMAX1(DABS(Y(K)),FLOOR)
             ENORM = DMAX1(ENORM,DABS(V(K))/WEIGHT)
20       CONTINUE
     ELSE IF (INORM .EQ. 3) THEN
         !                                  Same as DGEAR's error control
         ENORM = 0.0D0
         DO 30  K=1, NEQ
             WEIGHT = DMAX1(1.0D0,DABS(YMAX(K)))
             ENORM = ENORM + (V(K)/WEIGHT)**2
30       CONTINUE
         ENORM = DSQRT(ENORM)
     END IF
 !
9000 CONTINUE
     RETURN
 END
 !-----------------------------------------------------------------------
 !  IMSL Name:  I4PRK/DI4PRK (Single/Double precision version)
 !
 !  Computer:   pcdsms/DOUBLE
 !
 !  Revised:    January 29, 1985
 !
 !  Purpose:    Solve an initial value problem for ordinary differential
 !              equations using the Runge-Kutta-Verner fifth and sixth
 !              order method.
 !
 !  Usage:      CALL I4PRK (RK)
 !
 !  Argument:
 !     RK     - Vector of length 43 to be initialized.  (Output)
 !
 !  Chapter:    MATH/LIBRARY Integration and Differentiation
 !
 !  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
 !
 !  Warranty:   IMSL warrants only that IMSL testing has been applied
 !              to this code.  No other warranty, expressed or implied,
 !              is applicable.
 !
 !-----------------------------------------------------------------------
 !
 SUBROUTINE DI4PRK (RK)
     !                                  SPECIFICATIONS FOR ARGUMENTS
     DOUBLE PRECISION RK(*)
     !
     RK(1) = 1.0D0/6.0D0
     RK(2) = 4.0D0/75.0D0
     RK(3) = 16.0D0/75.0D0
     RK(4) = 5.0D0/6.0D0
     RK(5) = -8.0D0/3.0D0
     RK(6) = 5.0D0/2.0D0
     RK(7) = -165.0D0/64.0D0
     RK(8) = 55.0D0/6.0D0
     RK(9) = -425.0D0/64.0D0
     RK(10) = 85.0D0/96.0D0
     RK(11) = 12.0D0/5.0D0
     RK(12) = -8.0D0
     RK(13) = 4015.0D0/612.0D0
     RK(14) = -11.0D0/36.0D0
     RK(15) = 88.0D0/255.0D0
     RK(16) = -8263.0D0/15000.0D0
     RK(17) = 124.0D0/75.0D0
     RK(18) = -4501.0D0/4760.0D0
     RK(19) = -81.0D0/250.0D0
     RK(20) = 2484.0D0/10625.0D0
     RK(21) = 3501.0D0/1720.0D0
     RK(22) = -300.0D0/43.0D0
     RK(23) = 297275.0D0/52632.0D0
     RK(24) = -319.0D0/2322.0D0
     RK(25) = 24068.0D0/84065.0D0
     RK(26) = 0.0D0
     RK(27) = 3850.0D0/26703.0D0
     RK(28) = 3.0D0/40.0D0
     RK(29) = 0.0D0
     RK(30) = 875.0D0/2244.0D0
     RK(31) = 23.0D0/72.0D0
     RK(32) = 264.0D0/1955.0D0
     RK(33) = 0.0D0
     RK(34) = 125.0D0/11592.0D0
     RK(35) = 43.0D0/616.0D0
     RK(36) = 1.0D0/160.0D0
     RK(37) = 0.0D0
     RK(38) = 125.0D0/17952.0D0
     RK(39) = -1.0D0/144.0D0
     RK(40) = 84.0D0/13685.0D0
     RK(41) = 3.0D0/44.0D0
     RK(42) = -125.0D0/11592.0D0
     RK(43) = -43.0D0/616.0D0
     !
     RETURN
 END
 !-----------------------------------------------------------------------
 !  IMSL Name:  DMACH (Double precision version)
 !
 !  Computer:   pcdsms/DOUBLE
 !
 !  Revised:    March 15, 1984
 !
 !  Purpose:    Generate double precision machine constants.
 !
 !  Usage:      DMACH(N)
 !
 !  Arguments:
 !     N      - Index of desired constant.  (Input)
 !     DMACH  - Machine constant.  (Output)
 !              DMACH(1) = B**(EMIN-1), the smallest positive magnitude.
 !              DMACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
 !              DMACH(3) = B**(-T), the smallest relative spacing.
 !              DMACH(4) = B**(1-T), the largest relative spacing.
 !              DMACH(5) = LOG10(B), the log, base 10, of the radix.
 !              DMACH(6) = not-a-number.
 !              DMACH(7) = positive machine infinity.
 !              DMACH(8) = negative machine infinity.
 !
 !  GAMS:       R1
 !
 !  Chapters:   MATH/LIBRARY Reference Material
 !              STAT/LIBRARY Reference Material
 !              SFUN/LIBRARY Reference Material
 !
 !  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
 !
 !  Warranty:   IMSL warrants only that IMSL testing has been applied
 !              to this code.  No other warranty, expressed or implied,
 !              is applicable.
 !
 !-----------------------------------------------------------------------
 !
 DOUBLE PRECISION FUNCTION DMACH (N)
     !                                  SPECIFICATIONS FOR ARGUMENTS
     INTEGER    N
     !                                  SPECIFICATIONS FOR SAVE VARIABLES
     DOUBLE PRECISION RMACH(8)
     SAVE       RMACH
     !                                  SPECIFICATIONS FOR SUBROUTINES
     EXTERNAL   E1MES, E1POP, E1PSH, E1STI
     !                                  SPECIFICATIONS FOR LOCAL VARIABLES
     INTEGER    IRMACH(16)
     !
     !!      EQUIVALENCE (RMACH, IRMACH)
     !                                  DEFINE CONSTANTS
     DATA RMACH(1)/2.22559D-308/
     DATA RMACH(2)/1.79728D308/
     DATA RMACH(3)/1.11048D-16/
     DATA RMACH(4)/2.22096D-16/
     DATA RMACH(5)/.3010299956639811952137388947245D0/
     DATA IRMACH(11)/0/
     DATA IRMACH(12)/1206910591/
     DATA RMACH(7)/1.79728D308/
     DATA RMACH(8)/-1.79728D308/
     !
     IF (N.LT.1 .OR. N.GT.8) THEN
         CALL E1PSH ('DMACH ')
         DMACH = RMACH(6)
         CALL E1STI (1, N)
         CALL E1MES (5, 5, 'The argument must be between 1 '// &
             'and 8 inclusive. N = %(I1)')
         CALL E1POP ('DMACH ')
     ELSE
         DMACH = RMACH(N)
     END IF
     !
     RETURN
 END
 !-----------------------------------------------------------------------
 !  IMSL Name:  E1STD
 !
 !  Computer:   pcdsms/SINGLE
 !
 !  Revised:    March 6, 1984
 !
 !  Purpose:    To store a real number for subsequent use within an error
 !              message.
 !
 !  Usage:      CALL E1STD(ID, DVALUE)
 !
 !  Arguments:
 !     ID     - Integer specifying the substitution index.  ID must be
 !              between 1 and 9.  (Input)
 !     DVALUE - The double precision number to be stored.  (Input)
 !
 !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
 !
 !  Warranty:   IMSL warrants only that IMSL testing has been applied
 !              to this code.  No other warranty, expressed or implied,
 !              is applicable.
 !
 !-----------------------------------------------------------------------
 !
 SUBROUTINE E1STD (ID, DVALUE)
     !                                  SPECIFICATIONS FOR ARGUMENTS
     INTEGER    ID
     DOUBLE PRECISION DVALUE
     !                                  SPECIFICATIONS FOR LOCAL VARIABLES
     INTEGER    I, IBEG, ILEN
     CHARACTER  ARRAY(24), SAVE*24
     !                                  SPECIFICATIONS FOR SAVE VARIABLES
     INTEGER    IFINIT
     CHARACTER  BLANK(1)
     SAVE       BLANK, IFINIT
     !                                  SPECIFICATIONS FOR SPECIAL CASES
     !                              SPECIFICATIONS FOR COMMON /ERCOM1/
     INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
         PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
         IALLOC(51), HDRFMT(7), TRACON(7)
     COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
         PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
         TRACON
     SAVE       /ERCOM1/
     !                              SPECIFICATIONS FOR COMMON /ERCOM2/
     CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
     COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
     SAVE       /ERCOM2/
     !                              SPECIFICATIONS FOR COMMON /ERCOM3/
     DOUBLE PRECISION ERCKSM
     COMMON     /ERCOM3/ ERCKSM
     SAVE       /ERCOM3/
     !                              SPECIFICATIONS FOR COMMON /ERCOM4/
     LOGICAL    ISUSER(51)
     COMMON     /ERCOM4/ ISUSER
     SAVE       /ERCOM4/
     !                                  SPECIFICATIONS FOR SUBROUTINES
     EXTERNAL   E1INIT, E1INPL
     !                                  SPECIFICATIONS FOR FUNCTIONS
     EXTERNAL   I1ERIF
     INTEGER    I1ERIF
     !
     DATA BLANK/' '/, IFINIT/0/
     !                                  INITIALIZE IF NECESSARY
     IF (IFINIT .EQ. 0) THEN
         CALL E1INIT
         IFINIT = 1
     END IF
     IF (DVALUE .EQ. 0.0D0) THEN
         WRITE (SAVE,'(D24.15)') DVALUE
     ELSE
         WRITE (SAVE,'(1PE24.15E4)') DVALUE
     END IF
     DO 40  I=1, 24
40       ARRAY(I) = SAVE(I:I)
         IBEG = I1ERIF(ARRAY,24,BLANK,1)
         IF (ID.GE.1 .AND. ID.LE.9) THEN
             ILEN = 25 - IBEG
             CALL E1INPL ('D', ID, ILEN, ARRAY(IBEG))
         END IF
         !
         RETURN
     END
     !-----------------------------------------------------------------------
     !  IMSL Name:  E1USR
     !
     !  Computer:   pcdsms/SINGLE
     !
     !  Revised:    November 2, 1984
     !
     !  Purpose:    Set USER CODE switch.
     !
     !  Usage:      CALL E1USR(SWITCH)
     !
     !  Arguments:
     !     SWITCH - Character string.  (Input)
     !                'ON'  Indicates that USER CODE mode is being entered.
     !                'OFF' Indicates that USER CODE mode is being exited.
     !  Remarks:
     !     When E1POP is called from a routine while in USER CODE mode,
     !     then an error message of type 1-4 will be printed (if an error
     !     condition is in effect and the print table allows it).
     !     However, an error message of type 1-4 will never be printed
     !     if USER CODE mode is not in effect.
     !
     !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
     !
     !  Warranty:   IMSL warrants only that IMSL testing has been applied
     !              to this code.  No other warranty, expressed or implied,
     !              is applicable.
     !
     !-----------------------------------------------------------------------
     !
     SUBROUTINE E1USR (SWITCH)
         !                                  SPECIFICATIONS FOR ARGUMENTS
         CHARACTER  SWITCH*(*)
         !                                  SPECIFICATIONS FOR SAVE VARIABLES
         INTEGER    IFINIT
         SAVE       IFINIT
         !                                  SPECIFICATIONS FOR SPECIAL CASES
         !                              SPECIFICATIONS FOR COMMON /ERCOM1/
         INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
             PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
             IALLOC(51), HDRFMT(7), TRACON(7)
         COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
             PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
             TRACON
         SAVE       /ERCOM1/
         !                              SPECIFICATIONS FOR COMMON /ERCOM2/
         CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
         COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
         SAVE       /ERCOM2/
         !                              SPECIFICATIONS FOR COMMON /ERCOM3/
         DOUBLE PRECISION ERCKSM
         COMMON     /ERCOM3/ ERCKSM
         SAVE       /ERCOM3/
         !                              SPECIFICATIONS FOR COMMON /ERCOM4/
         LOGICAL    ISUSER(51)
         COMMON     /ERCOM4/ ISUSER
         SAVE       /ERCOM4/
         !                                  SPECIFICATIONS FOR SUBROUTINES
         EXTERNAL   E1INIT, E1MES, E1STL
         !
         DATA IFINIT/0/
         !                                  INITIALIZE ERROR TABLE IF NECESSARY
         IF (IFINIT .EQ. 0) THEN
             CALL E1INIT
             IFINIT = 1
         END IF
         IF (SWITCH.EQ.'ON' .OR. SWITCH.EQ.'on') THEN
             ISUSER(CALLVL) = .TRUE.
         ELSE IF (SWITCH.EQ.'OFF' .OR. SWITCH.EQ.'off') THEN
             ISUSER(CALLVL) = .FALSE.
         ELSE
             CALL E1STL (1, SWITCH)
             CALL E1MES (5, 1, 'Invalid value for SWITCH in call to'// &
                 ' E1USR.  SWITCH must be set to ''ON'' or '// &
                 '''OFF''.  SWITCH = ''%(L1)'' ')
         END IF
         !
         RETURN
     END
     !-----------------------------------------------------------------------
     !  IMSL Name:  I1KNR
     !
     !  Computer:   pcdsms/SINGLE
     !
     !  Revised:    August 29, 1983
     !
     !  Purpose:    To cause the workspace which was allocated at the
     !              current level to NOT be released at that level.
     !
     !  Usage:      CALL I1KNR
     !
     !  Arguments:  None
     !
     !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
     !
     !  Warranty:   IMSL warrants only that IMSL testing has been applied
     !              to this code.  No other warranty, expressed or implied,
     !              is applicable.
     !
     !-----------------------------------------------------------------------
     !
     SUBROUTINE I1KNR
         !                                  SPECIFICATIONS FOR SPECIAL CASES
         !                              SPECIFICATIONS FOR COMMON /ERCOM1/
         INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
             PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
             IALLOC(51), HDRFMT(7), TRACON(7)
         COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
             PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
             TRACON
         SAVE       /ERCOM1/
         !                              SPECIFICATIONS FOR COMMON /ERCOM2/
         CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
         COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
         SAVE       /ERCOM2/
         !                              SPECIFICATIONS FOR COMMON /ERCOM3/
         DOUBLE PRECISION ERCKSM
         COMMON     /ERCOM3/ ERCKSM
         SAVE       /ERCOM3/
         !                              SPECIFICATIONS FOR COMMON /ERCOM4/
         LOGICAL    ISUSER(51)
         COMMON     /ERCOM4/ ISUSER
         SAVE       /ERCOM4/
         !                                  SPECIFICATIONS FOR FUNCTIONS
         EXTERNAL   I1KST
         INTEGER    I1KST
         !
         IALLOC(CALLVL-1) = I1KST(1)
         !
         RETURN
     END
     !-----------------------------------------------------------------------
     !  IMSL Name:  IDAMAX (Single precision version)
     !
     !  Computer:   pcdsms/SINGLE
     !
     !  Revised:    August 9, 1986
     !
     !  Purpose:    Find the smallest index of the component of a
     !              double-precision vector having maximum absolute value.
     !
     !  Usage:      IDAMAX(N, DX, INCX)
     !
     !  Arguments:
     !     N      - Length of vector X.  (Input)
     !     DX     - Double precision vector of length N*INCX.  (Input)
     !     INCX   - Displacement between elements of DX.  (Input)
     !              X(I) is defined to be DX(1+(I-1)*INCX). INCX must be
     !              greater than zero.
     !     IDAMAX - The smallest index I such that DABS(X(I)) is the maximum
     !              of DABS(X(J)) for J=1 to N.  (Output)
     !              X(I) refers to a specific element of DX. See INCX
     !              argument description.
     !
     !  GAMS:       D1a2
     !
     !  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
     !              STAT/LIBRARY Mathematical Support
     !
     !  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
     !
     !  Warranty:   IMSL warrants only that IMSL testing has been applied
     !              to this code.  No other warranty, expressed or implied,
     !              is applicable.
     !
     !-----------------------------------------------------------------------
     !
     INTEGER FUNCTION IDAMAX (N, DX, INCX)
         !                                  SPECIFICATIONS FOR ARGUMENTS
         INTEGER    N, INCX
         DOUBLE PRECISION DX(*)
         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
         INTEGER    I, II, NS
         DOUBLE PRECISION DMAX, XMAG
         !                                  SPECIFICATIONS FOR INTRINSICS
         !     INTRINSIC  DABS
         INTRINSIC  DABS
         DOUBLE PRECISION DABS
         !
         IDAMAX = 0
         IF (N .GE. 1) THEN
             IDAMAX = 1
             IF (N .GT. 1) THEN
                 IF (INCX .NE. 1) THEN
                     !                                  CODE FOR INCREMENTS NOT EQUAL TO 1.
                     DMAX = DABS(DX(1))
                     NS = N*INCX
                     II = 1
                     DO 10  I=1, NS, INCX
                         XMAG = DABS(DX(I))
                         IF (XMAG .GT. DMAX) THEN
                             IDAMAX = II
                             DMAX = XMAG
                         END IF
                         II = II + 1
10                   CONTINUE
                 ELSE
                     !                                  CODE FOR INCREMENTS EQUAL TO 1.
                     DMAX = DABS(DX(1))
                     DO 20  I=2, N
                         XMAG = DABS(DX(I))
                         IF (XMAG .GT. DMAX) THEN
                             IDAMAX = I
                             DMAX = XMAG
                         END IF
20                   CONTINUE
                 END IF
             END IF
         END IF
         RETURN
     END
     !-----------------------------------------------------------------------
     !  IMSL Name:  SSET (Single precision version)
     !
     !  Computer:   pcdsms/SINGLE
     !
     !  Revised:    August 9, 1986
     !
     !  Purpose:    Set the components of a vector to a scalar, all
     !              single precision.
     !
     !  Usage:      CALL SSET (N, SA, SX, INCX)
     !
     !  Arguments:
     !     N      - Length of vector X.  (Input)
     !     SA     - Real scalar.  (Input)
     !     SX     - Real vector of length N*INCX.  (Input/Output)
     !              SSET replaces X(I) with SA for I=1,...,N. X(I) refers to
     !              a specific element of SX. See INCX argument description.
     !     INCX   - Displacement between elements of SX.  (Input)
     !              X(I) is defined to be SX(1+(I-1)*INCX). INCX must be
     !              greater than zero.
     !
     !  GAMS:       D1a1
     !
     !  Chapters:   MATH/LIBRARY Basic Matrix/Vector Operations
     !              STAT/LIBRARY Mathematical Support
     !
     !  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
     !
     !  Warranty:   IMSL warrants only that IMSL testing has been applied
     !              to this code.  No other warranty, expressed or implied,
     !              is applicable.
     !
     !-----------------------------------------------------------------------
     !
     SUBROUTINE SSET (N, SA, SX, INCX)
         !                                  SPECIFICATIONS FOR ARGUMENTS
         INTEGER    N, INCX
         REAL       SA, SX(*)
         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
         INTEGER    I, M, MP1, NINCX
         !                                  SPECIFICATIONS FOR SPECIAL CASES
         !     INTRINSIC  MOD
         INTRINSIC  MOD
         INTEGER    MOD
         !
         IF (N .GT. 0) THEN
             IF (INCX .NE. 1) THEN
                 !                                  CODE FOR INCREMENT NOT EQUAL TO 1
                 NINCX = N*INCX
                 DO 10  I=1, NINCX, INCX
                     SX(I) = SA
10               CONTINUE
             ELSE
                 !                                  CODE FOR INCREMENT EQUAL TO 1
                 M = MOD(N,8)
                 !                                  CLEAN-UP LOOP
                 DO 30  I=1, M
                     SX(I) = SA
30               CONTINUE
                 MP1 = M + 1
                 DO 40  I=MP1, N, 8
                     SX(I) = SA
                     SX(I+1) = SA
                     SX(I+2) = SA
                     SX(I+3) = SA
                     SX(I+4) = SA
                     SX(I+5) = SA
                     SX(I+6) = SA
                     SX(I+7) = SA
40               CONTINUE
             END IF
         END IF
         RETURN
     END

     !-----------------------------------------------------------------------
     !  IMSL Name:  E1INIT
     !
     !  Computer:   pcdsms/SINGLE
     !
     !  Revised:    March 13, 1984
     !
     !  Purpose:    Initialization.
     !
     !  Usage:      CALL E1INIT
     !
     !  Arguments:  None
     !
     !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
     !
     !  Warranty:   IMSL warrants only that IMSL testing has been applied
     !              to this code.  No other warranty, expressed or implied,
     !              is applicable.
     !
     !-----------------------------------------------------------------------
     !
     SUBROUTINE E1INIT
         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
         INTEGER    L
         !                                  SPECIFICATIONS FOR SAVE VARIABLES
         INTEGER    ISINIT
         SAVE       ISINIT
         !                                  SPECIFICATIONS FOR SPECIAL CASES
         !                              SPECIFICATIONS FOR COMMON /ERCOM1/
         INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
             PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
             IALLOC(51), HDRFMT(7), TRACON(7)
         COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
             PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
             TRACON
         SAVE       /ERCOM1/
         !                              SPECIFICATIONS FOR COMMON /ERCOM2/
         CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
         COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
         SAVE       /ERCOM2/
         !                              SPECIFICATIONS FOR COMMON /ERCOM3/
         DOUBLE PRECISION ERCKSM
         COMMON     /ERCOM3/ ERCKSM
         SAVE       /ERCOM3/
         !                              SPECIFICATIONS FOR COMMON /ERCOM4/
         LOGICAL    ISUSER(51)
         COMMON     /ERCOM4/ ISUSER
         SAVE       /ERCOM4/
         !                              SPECIFICATIONS FOR COMMON /ERCOM8/
         INTEGER    PROLVL, XXLINE(10), XXPLEN(10), ICALOC(10), INALOC(10)
         COMMON     /ERCOM8/ PROLVL, XXLINE, XXPLEN, ICALOC, INALOC
         SAVE       /ERCOM8/
         !                              SPECIFICATIONS FOR COMMON /ERCOM9/
         CHARACTER  XXPROC(10)*31
         COMMON     /ERCOM9/ XXPROC
         SAVE       /ERCOM9/
         !
         DATA ISINIT/0/
         !
         IF (ISINIT .EQ. 0) THEN
             !                                  INITIALIZE
             CALLVL = 1
             ERCODE(1) = 0
             ERTYPE(1) = 0
             IALLOC(1) = 0
             ISUSER(1) = .TRUE.
             IFERR6 = 0
             IFERR7 = 0
             PLEN = 1
             MAXLEV = 50
             DO 10  L=2, 51
                 ERTYPE(L) = -1
                 ERCODE(L) = -1
                 IALLOC(L) = 0
                 ISUSER(L) = .FALSE.
10           CONTINUE
             DO 20  L=1, 7
                 HDRFMT(L) = 1
                 TRACON(L) = 1
20           CONTINUE
             PROLVL = 1
             DO 30  L=1, 10
30               ICALOC(L) = 0
                 XXLINE(1) = 0
                 XXPLEN(1) = 1
                 XXPROC(1) = '?'
                 RNAME(1) = 'USER'
                 PRINTB(1) = 0
                 PRINTB(2) = 0
                 DO 40  L=3, 7
40                   PRINTB(L) = 1
                     STOPTB(1) = 0
                     STOPTB(2) = 0
                     STOPTB(3) = 0
                     STOPTB(4) = 1
                     STOPTB(5) = 1
                     STOPTB(6) = 0
                     STOPTB(7) = 1
                     ERCKSM = 0.0D0
                     !                                  SET FLAG TO INDICATE THAT
                     !                                    INITIALIZATION HAS OCCURRED
                     ISINIT = 1
                 END IF
                 !
                 RETURN
             END
             !-----------------------------------------------------------------------
             !  IMSL Name:  E1INPL
             !
             !  Computer:   pcdsms/SINGLE
             !
             !  Revised:    March 2, 1984
             !
             !  Purpose:    To store a character string in the parameter list PLIST
             !              for use by the error message handler.
             !
             !  Usage:      CALL E1INPL(FORM,NUM,SLEN,STRUP)
             !
             !  Arguments:
             !     FORM   - A character string of length one to be inserted into
             !              PLIST which specifies the form of the string.  (Input)
             !              For example, 'L' for string, 'A' for character array,
             !              'I' for integer, 'K' for keyword (PROTRAN only).  An
             !              asterisk is inserted into PLIST preceding FORM.
             !     NUM    - Integer to be inserted as a character into PLIST
             !              immediately following FORM.  (Input)  NUM must be between
             !              1 and 9.
             !     SLEN   - The number of characters in STRUP.  (Input)  LEN must be
             !              less than or equal to 255.  The character representation
             !              of SLEN is inserted into PLIST after NUM and an asterisk.
             !     STRUP  - A character string of length LEN which is to be inserted
             !              into PLIST.  (Input)  Trailing blanks are ignored.
             !
             !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
             !
             !  Warranty:   IMSL warrants only that IMSL testing has been applied
             !              to this code.  No other warranty, expressed or implied,
             !              is applicable.
             !
             !-----------------------------------------------------------------------
             !
             SUBROUTINE E1INPL (FORM, NUM, SLEN, STRUP)
                 !                                  SPECIFICATIONS FOR ARGUMENTS
                 INTEGER    NUM, SLEN
                 CHARACTER  FORM, STRUP(*)
                 !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                 INTEGER    IER, L, LEN2, LENCK, LOC, NLEN, NNUM
                 CHARACTER  STRNCH(3)
                 !                                  SPECIFICATIONS FOR SAVE VARIABLES
                 CHARACTER  BLANK, PRCNT(1), TEMP(4)
                 SAVE       BLANK, PRCNT, TEMP
                 !                                  SPECIFICATIONS FOR SPECIAL CASES
                 !                              SPECIFICATIONS FOR COMMON /ERCOM1/
                 INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                     PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                     IALLOC(51), HDRFMT(7), TRACON(7)
                 COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                     PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                     TRACON
                 SAVE       /ERCOM1/
                 !                              SPECIFICATIONS FOR COMMON /ERCOM2/
                 CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
                 COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
                 SAVE       /ERCOM2/
                 !                              SPECIFICATIONS FOR COMMON /ERCOM3/
                 DOUBLE PRECISION ERCKSM
                 COMMON     /ERCOM3/ ERCKSM
                 SAVE       /ERCOM3/
                 !                              SPECIFICATIONS FOR COMMON /ERCOM4/
                 LOGICAL    ISUSER(51)
                 COMMON     /ERCOM4/ ISUSER
                 SAVE       /ERCOM4/
                 !                                  SPECIFICATIONS FOR INTRINSICS
                 !     INTRINSIC  IABS
                 INTRINSIC  IABS
                 INTEGER    IABS
                 !                                  SPECIFICATIONS FOR SUBROUTINES
                 EXTERNAL   C1TIC, M1VE
                 !
                 DATA TEMP/'*', ' ', ' ', '*'/, PRCNT/'%'/, BLANK/' '/
                 !
                 NNUM = IABS(NUM)
                 LENCK = PLEN + SLEN + 8
                 IF (NNUM.GE.1 .AND. NNUM.LE.9 .AND. LENCK.LE.300) THEN
                     TEMP(2) = FORM
                     CALL C1TIC (NNUM, TEMP(3), 1, IER)
                     LOC = PLEN + 1
                     IF (LOC .EQ. 2) LOC = 1
                     CALL M1VE (TEMP, 1, 4, 4, PLIST(LOC), 1, 4, 262, IER)
                     LOC = LOC + 4
                     IF (NUM .LT. 0) THEN
                         LEN2 = SLEN
                     ELSE
                         DO 10  L=1, SLEN
                             LEN2 = SLEN - L + 1
                             IF (STRUP(LEN2) .NE. BLANK) GO TO 20
10                       CONTINUE
                         LEN2 = 1
20                   CONTINUE
                 END IF
                 NLEN = 1
                 IF (LEN2 .GE. 10) NLEN = 2
                 IF (LEN2 .GE. 100) NLEN = 3
                 CALL C1TIC (LEN2, STRNCH, NLEN, IER)
                 CALL M1VE (STRNCH, 1, NLEN, 3, PLIST(LOC), 1, NLEN, 262, IER)
                 LOC = LOC + NLEN
                 CALL M1VE (PRCNT, 1, 1, 1, PLIST(LOC), 1, 1, 262, IER)
                 LOC = LOC + 1
                 CALL M1VE (STRUP, 1, LEN2, LEN2, PLIST(LOC), 1, LEN2, 262, &
                     IER)
                 PLEN = LOC + LEN2 - 1
             END IF
             !
             RETURN
         END
         !-----------------------------------------------------------------------
         !  IMSL Name:  E1MES
         !
         !  Computer:   pcdsms/SINGLE
         !
         !  Revised:    March 2, 1984
         !
         !  Purpose:    Set an error state for the current level in the stack.
         !              The message is printed immediately if the error type is
         !              5, 6, or 7 and the print attribute for that type is YES.
         !
         !  Usage:      CALL E1MES(IERTYP,IERCOD,MSGPKD)
         !
         !  Arguments:
         !     IERTYP - Integer specifying the error type.  (Input)
         !                IERTYP=1,  informational/note
         !                IERTYP=2,  informational/alert
         !                IERTYP=3,  informational/warning
         !                IERTYP=4,  informational/fatal
         !                IERTYP=5,  terminal
         !                IERTYP=6,  PROTRAN/warning
         !                IERTYP=7,  PROTRAN/fatal
         !     IERCOD - Integer specifying the error code.  (Input)
         !     MSGPKD - A character string containing the message.
         !              (Input)  Within the message, any of following may appear
         !                %(A1),%(A2),...,%(A9) for character arrays
         !                %(C1),%(C2),...,%(C9) for complex numbers
         !                %(D1),%(D2),...,%(D9) for double precision numbers
         !                %(I1),%(I2),...,%(I9) for integer numbers
         !                %(K1),%(K2),...,%(K9) for keywords
         !                %(L1),%(L2),...,%(L9) for literals (strings)
         !                %(R1),%(R2),...,%(R9) for real numbers
         !                %(Z1),%(Z2),...,%(Z9) for double complex numbers
         !              This provides a way to insert character arrays, strings,
         !              numbers, and keywords into the message.  See remarks
         !              below.
         !
         !  Remarks:
         !     The number of characters in the message after the insertion of
         !     the corresponding strings, etc. should not exceed 255.  If the
         !     limit is exceeded, only the first 255 characters will be used.
         !     The appropriate strings, etc. need to have been previously stored
         !     in common via calls to E1STA, E1STD, etc.  Line breaks may be
         !     specified by inserting the two characters '%/' into the message
         !     at the desired locations.
         !
         !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
         !
         !  Warranty:   IMSL warrants only that IMSL testing has been applied
         !              to this code.  No other warranty, expressed or implied,
         !              is applicable.
         !
         !-----------------------------------------------------------------------
         !
         SUBROUTINE E1MES (IERTYP, IERCOD, MSGPKD)
             !                                  SPECIFICATIONS FOR ARGUMENTS
             INTEGER    IERTYP, IERCOD
             CHARACTER  MSGPKD*(*)
             !                                  SPECIFICATIONS FOR LOCAL VARIABLES
             INTEGER    ERTYP2, I, IER, IPLEN, ISUB, LAST, LEN2, LOC, M, MS, &
                 NLOC, NUM, PBEG
             CHARACTER  MSGTMP(255)
             !                                  SPECIFICATIONS FOR SAVE VARIABLES
             INTEGER    IFINIT, NFORMS
             CHARACTER  BLNK, DBB(3), FIND(4), FORMS(9), INREF(25), LPAR, &
                 NCHECK(3), PERCNT, RPAR
             SAVE       BLNK, DBB, FIND, FORMS, IFINIT, INREF, LPAR, NCHECK, &
                 NFORMS, PERCNT, RPAR
             !                                  SPECIFICATIONS FOR SPECIAL CASES
             !                              SPECIFICATIONS FOR COMMON /ERCOM1/
             INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                 PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                 IALLOC(51), HDRFMT(7), TRACON(7)
             COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                 PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                 TRACON
             SAVE       /ERCOM1/
             !                              SPECIFICATIONS FOR COMMON /ERCOM2/
             CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
             COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
             SAVE       /ERCOM2/
             !                              SPECIFICATIONS FOR COMMON /ERCOM3/
             DOUBLE PRECISION ERCKSM
             COMMON     /ERCOM3/ ERCKSM
             SAVE       /ERCOM3/
             !                              SPECIFICATIONS FOR COMMON /ERCOM4/
             LOGICAL    ISUSER(51)
             COMMON     /ERCOM4/ ISUSER
             SAVE       /ERCOM4/
             !                                  SPECIFICATIONS FOR INTRINSICS
             !     INTRINSIC  LEN,MIN0
             INTRINSIC  LEN, MIN0
             INTEGER    LEN, MIN0
             !                                  SPECIFICATIONS FOR SUBROUTINES
             EXTERNAL   C1TCI, E1INIT, E1PRT, E1UCS, M1VE, M1VECH
             !                                  SPECIFICATIONS FOR FUNCTIONS
             EXTERNAL   I1DX
             INTEGER    I1DX
             !
             DATA FORMS/'A', 'C', 'D', 'I', 'K', 'L', 'R', 'S', 'Z'/, &
                 NFORMS/9/
             DATA PERCNT/'%'/, LPAR/'('/, RPAR/')'/, BLNK/' '/
             DATA INREF/' ', 'i', 'n', ' ', 'r', 'e', 'f', 'e', 'r', &
                 'e', 'n', 'c', 'e', ' ', 't', 'o', ' ', 'k', 'e', &
                 'y', 'w', 'o', 'r', 'd', ' '/
             DATA NCHECK/'N', '1', '*'/, DBB/'.', ' ', ' '/
             DATA FIND/'*', ' ', ' ', '*'/
             DATA IFINIT/0/
             !                                  INITIALIZE ERROR TABLE IF NECESSARY
             IF (IFINIT .EQ. 0) THEN
                 CALL E1INIT
                 IFINIT = 1
             END IF
             !                                  CHECK AND SET ERROR TYPE IF NECESSARY
             IF (IERTYP .NE. -1) THEN
                 ERTYPE(CALLVL) = IERTYP
             ELSE IF (IERTYP.LT.-1 .OR. IERTYP.GT.7) THEN
                 MSGLEN = 51
                 CALL M1VECH ('.  Error from E1MES.  Illegal error type'// &
                     ' specified. ', MSGLEN, MSGSAV, MSGLEN)
                 CALL E1PRT
                 STOP
             END IF
             !
             ERTYP2 = ERTYPE(CALLVL)
             !                                  SET ERROR CODE IF NECESSARY
             IF (IERCOD .GT. -1) ERCODE(CALLVL) = IERCOD
             LEN2 = LEN(MSGPKD)
             !
             IF (IERTYP.EQ.0 .OR. IERCOD.EQ.0) THEN
                 !                                  REMOVE THE ERROR STATE
                 MSGLEN = 0
             ELSE IF (LEN2.EQ.0 .OR. (LEN2.EQ.1.AND.MSGPKD(1:1).EQ.BLNK)) THEN
                 IF (ERTYP2 .EQ. 6) IFERR6 = 1
                 IF (ERTYP2 .EQ. 7) IFERR7 = 1
                 !                                  UPDATE CHECKSUM PARAMETER ERCKSM
                 CALL E1UCS
                 !                                  PRINT MESSAGE IF NECESSARY
                 IF (ERTYP2.GE.5 .AND. PRINTB(ERTYP2).EQ.1) CALL E1PRT
             ELSE
                 !                                  FILL UP MSGSAV WITH EXPANDED MESSAGE
                 LEN2 = MIN0(LEN2,255)
                 DO 10  I=1, LEN2
                     MSGTMP(I) = MSGPKD(I:I)
10               CONTINUE
                 MS = 0
                 M = 0
                 !                                  CHECK PLIST FOR KEYWORD NAME
                 NLOC = I1DX(PLIST,PLEN,NCHECK,3)
                 IF (NLOC.GT.0 .AND. HDRFMT(ERTYP2).EQ.3) THEN
                     !                                  M1VE INREF INTO MSGSAV
                     CALL M1VE (INREF, 1, 25, 25, MSGSAV, 1, 25, 25, IER)
                     !                                  GET LENGTH OF KEYWORD NAME
                     CALL C1TCI (PLIST(NLOC+3), 3, IPLEN, IER)
                     PBEG = NLOC + 3 + IER
                     !                                  M1VE KEYWORD NAME INTO MSGSAV
                     CALL M1VE (PLIST, PBEG, PBEG+IPLEN-1, PLEN, MSGSAV, 26, &
                         IPLEN+25, 255, IER)
                     !                                  UPDATE POINTER
                     MS = IPLEN + 25
                 END IF
                 !                                  INSERT DOT, BLANK, BLANK
                 CALL M1VE (DBB, 1, 3, 3, MSGSAV, MS+1, MS+3, 255, IER)
                 MS = MS + 3
                 !                                  LOOK AT NEXT CHARACTER
20               M = M + 1
                 ISUB = 0
                 IF (M .GT. LEN2-4) THEN
                     LAST = LEN2 - M + 1
                     DO 30  I=1, LAST
30                       MSGSAV(MS+I) = MSGTMP(M+I-1)
                         MSGLEN = MS + LAST
                         GO TO 40
                     ELSE IF (MSGTMP(M).EQ.PERCNT .AND. MSGTMP(M+1).EQ.LPAR .AND. &
                         MSGTMP(M+4).EQ.RPAR) THEN
                         CALL C1TCI (MSGTMP(M+3), 1, NUM, IER)
                         IF (IER.EQ.0 .AND. NUM.NE.0 .AND. I1DX(FORMS,NFORMS, &
                             MSGTMP(M+2),1).NE.0) THEN
                             !                                  LOCATE THE ITEM IN THE PARAMETER LIST
                             CALL M1VE (MSGTMP(M+2), 1, 2, 2, FIND, 2, 3, 4, IER)
                             LOC = I1DX(PLIST,PLEN,FIND,4)
                             IF (LOC .GT. 0) THEN
                                 !                                  SET IPLEN = LENGTH OF STRING
                                 CALL C1TCI (PLIST(LOC+4), 4, IPLEN, IER)
                                 PBEG = LOC + 4 + IER
                                 !                                  ADJUST IPLEN IF IT IS TOO BIG
                                 IPLEN = MIN0(IPLEN,255-MS)
                                 !                                  M1VE STRING FROM PLIST INTO MSGSAV
                                 CALL M1VE (PLIST, PBEG, PBEG+IPLEN-1, PLEN, MSGSAV, &
                                     MS+1, MS+IPLEN, 255, IER)
                                 IF (IER.GE.0 .AND. IER.LT.IPLEN) THEN
                                     !                                  UPDATE POINTERS
                                     M = M + 4
                                     MS = MS + IPLEN - IER
                                     !                                  BAIL OUT IF NO MORE ROOM
                                     IF (MS .GE. 255) THEN
                                         MSGLEN = 255
                                         GO TO 40
                                     END IF
                                     !                                  SET FLAG TO SHOW SUBSTITION WAS MADE
                                     ISUB = 1
                                 END IF
                             END IF
                         END IF
                     END IF
                     IF (ISUB .EQ. 0) THEN
                         MS = MS + 1
                         MSGSAV(MS) = MSGTMP(M)
                     END IF
                     GO TO 20
40                   ERTYP2 = ERTYPE(CALLVL)
                     IF (ERTYP2 .EQ. 6) IFERR6 = 1
                     IF (ERTYP2 .EQ. 7) IFERR7 = 1
                     !                                  UPDATE CHECKSUM PARAMETER ERCKSM
                     CALL E1UCS
                     !                                  PRINT MESSAGE IF NECESSARY
                     IF (ERTYP2.GE.5 .AND. PRINTB(ERTYP2).EQ.1) CALL E1PRT
                 END IF
                 !                                  CLEAR PARAMETER LIST
                 PLEN = 1
                 !
                 RETURN
             END
             !-----------------------------------------------------------------------
             !  IMSL Name:  E1POP
             !
             !  Computer:   pcdsms/SINGLE
             !
             !  Revised:    March 13, 1984
             !
             !  Purpose:    To pop a subroutine name from the error control stack.
             !
             !  Usage:      CALL E1POP(NAME)
             !
             !  Arguments:
             !     NAME   - A character string of length six specifying the name
             !              of the subroutine.  (Input)
             !
             !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
             !
             !  Warranty:   IMSL warrants only that IMSL testing has been applied
             !              to this code.  No other warranty, expressed or implied,
             !              is applicable.
             !
             !-----------------------------------------------------------------------
             !
             SUBROUTINE E1POP (NAME)
                 !                                  SPECIFICATIONS FOR ARGUMENTS
                 CHARACTER  NAME*(*)
                 !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                 INTEGER    IERTYP, IR
                 !                                  SPECIFICATIONS FOR SPECIAL CASES
                 !                              SPECIFICATIONS FOR COMMON /ERCOM1/
                 INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                     PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                     IALLOC(51), HDRFMT(7), TRACON(7)
                 COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                     PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                     TRACON
                 SAVE       /ERCOM1/
                 !                              SPECIFICATIONS FOR COMMON /ERCOM2/
                 CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
                 COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
                 SAVE       /ERCOM2/
                 !                              SPECIFICATIONS FOR COMMON /ERCOM3/
                 DOUBLE PRECISION ERCKSM
                 COMMON     /ERCOM3/ ERCKSM
                 SAVE       /ERCOM3/
                 !                              SPECIFICATIONS FOR COMMON /ERCOM4/
                 LOGICAL    ISUSER(51)
                 COMMON     /ERCOM4/ ISUSER
                 SAVE       /ERCOM4/
                 !                                  SPECIFICATIONS FOR SUBROUTINES
                 EXTERNAL   E1MES, E1PRT, E1PSH, E1STI, E1STL, I1KRL
                 !                                  SPECIFICATIONS FOR FUNCTIONS
                 EXTERNAL   I1KST
                 INTEGER    I1KST
                 !
                 IF (CALLVL .LE. 1) THEN
                     CALL E1PSH ('E1POP ')
                     CALL E1STL (1, NAME)
                     CALL E1MES (5, 1, 'Error condition in E1POP.  Cannot pop '// &
                         'from %(L1) because stack is empty.')
                     STOP
                 ELSE IF (NAME .NE. RNAME(CALLVL)) THEN
                     CALL E1STL (1, NAME)
                     CALL E1STL (2, RNAME(CALLVL))
                     CALL E1MES (5, 2, 'Error condition in E1POP.  %(L1) does '// &
                         'not match the name %(L2) in the stack.')
                     STOP
                 ELSE
                     IERTYP = ERTYPE(CALLVL)
                     IF (IERTYP .NE. 0) THEN
                         !                                  M1VE ERROR TYPE AND ERROR CODE TO
                         !                                    PREVIOUS LEVEL FOR ERROR TYPES 2-7
                         IF (IERTYP.GE.2 .AND. IERTYP.LE.7) THEN
                             ERTYPE(CALLVL-1) = ERTYPE(CALLVL)
                             ERCODE(CALLVL-1) = ERCODE(CALLVL)
                         END IF
                         !                                  CHECK PRINT TABLE TO DETERMINE
                         !                                    WHETHER TO PRINT STORED MESSAGE
                         IF (IERTYP .LE. 4) THEN
                             IF (ISUSER(CALLVL-1) .AND. PRINTB(IERTYP).EQ.1) &
                                 CALL E1PRT
                         ELSE
                             IF (PRINTB(IERTYP) .EQ. 1) CALL E1PRT
                         END IF
                         !                                  CHECK STOP TABLE AND ERROR TYPE TO
                         !                                    DETERMINE WHETHER TO STOP
                         IF (IERTYP .LE. 4) THEN
                             IF (ISUSER(CALLVL-1) .AND. STOPTB(IERTYP).EQ.1) THEN
                                 STOP
                             END IF
                         ELSE IF (IERTYP .EQ. 5) THEN
                             IF (STOPTB(IERTYP) .EQ. 1) THEN
                                 STOP
                             END IF
                         ELSE IF (HDRFMT(IERTYP) .EQ. 1) THEN
                             IF (ISUSER(CALLVL-1)) THEN
                                 IF (N1RGB(0) .NE. 0) THEN
                                     STOP
                                 END IF
                             END IF
                         END IF
                     END IF
                     !                                  SET ERROR TYPE AND CODE
                     IF (CALLVL .LT. MAXLEV) THEN
                         ERTYPE(CALLVL+1) = -1
                         ERCODE(CALLVL+1) = -1
                     END IF
                     !                                  SET IR = AMOUNT OF WORKSPACE
                     !                                  ALLOCATED AT THIS LEVEL
                     IR = I1KST(1) - IALLOC(CALLVL-1)
                     IF (IR .GT. 0) THEN
                         !                                  RELEASE WORKSPACE
                         CALL I1KRL (IR)
                         IALLOC(CALLVL) = 0
                     ELSE IF (IR .LT. 0) THEN
                         CALL E1STI (1, CALLVL)
                         CALL E1STI (2, IALLOC(CALLVL-1))
                         CALL E1STI (3, I1KST(1))
                         CALL E1MES (5, 3, 'Error condition in E1POP. '// &
                             ' The number of workspace allocations at '// &
                             'level %(I1) is %(I2).  However, the total '// &
                             'number of workspace allocations is %(I3).')
                         STOP
                     END IF
                     !                                  DECREASE THE STACK POINTER BY ONE
                     CALLVL = CALLVL - 1
                 END IF
                 !
                 RETURN
             END
             !-----------------------------------------------------------------------
             !  IMSL Name:  E1POS
             !
             !  Computer:   pcdsms/SINGLE
             !
             !  Revised:    March 2, 1984
             !
             !  Purpose:    Set or retrieve print and stop attributes.
             !
             !  Usage:      CALL E1POS(IERTYP,IPATT,ISATT)
             !
             !  Arguments:
             !     IERTYP - Integer specifying the error type for which print and
             !              stop attributes are to be set or retrieved.  (Input)  If
             !              IERTYP is 0 then the settings apply to all error types.
             !              If IERTYP is between 1 and 7, then the settings only
             !              apply to that specified error type.  If IERTYP is
             !              negative then the current print and stop attributes will
             !              be returned in IPATT and ISATT.
             !     IPATT  - If IERTYP is positive, IPATT is an integer specifying the
             !              desired print attribute as follows: -1 means no change,
             !              0 means NO, 1 means YES, and 2 means assign the default
             !              setting.  (Input)  If IERTYP is negative, IPATT is
             !              returned as 1 if print is YES or 0 if print is NO for
             !              error type IABS(IERTYP).  (Output)
             !     ISATT  - If IERTYP is positive, ISATT is an integer specifying the
             !              desired stop attribute as follows: -1 means no change,
             !              0 means NO, 1 means YES, and 2 means assign the default
             !              setting.  (Input)  If IERTYP is negative, ISATT is
             !              returned as 1 if print is YES or 0 if print is NO for
             !              error type IABS(IERTYP).  (Output)
             !
             !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
             !
             !  Warranty:   IMSL warrants only that IMSL testing has been applied
             !              to this code.  No other warranty, expressed or implied,
             !              is applicable.
             !
             !-----------------------------------------------------------------------
             !
             SUBROUTINE E1POS (IERTYP, IPATT, ISATT)
                 !                                  SPECIFICATIONS FOR ARGUMENTS
                 INTEGER    IERTYP, IPATT, ISATT
                 !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                 INTEGER    I, IER
                 !                                  SPECIFICATIONS FOR SAVE VARIABLES
                 INTEGER    DEFLTP(7), DEFLTS(7), IFINIT
                 SAVE       DEFLTP, DEFLTS, IFINIT
                 !                                  SPECIFICATIONS FOR SPECIAL CASES
                 !                              SPECIFICATIONS FOR COMMON /ERCOM1/
                 INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                     PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                     IALLOC(51), HDRFMT(7), TRACON(7)
                 COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                     PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                     TRACON
                 SAVE       /ERCOM1/
                 !                              SPECIFICATIONS FOR COMMON /ERCOM2/
                 CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
                 COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
                 SAVE       /ERCOM2/
                 !                              SPECIFICATIONS FOR COMMON /ERCOM3/
                 DOUBLE PRECISION ERCKSM
                 COMMON     /ERCOM3/ ERCKSM
                 SAVE       /ERCOM3/
                 !                              SPECIFICATIONS FOR COMMON /ERCOM4/
                 LOGICAL    ISUSER(51)
                 COMMON     /ERCOM4/ ISUSER
                 SAVE       /ERCOM4/
                 !                                  SPECIFICATIONS FOR INTRINSICS
                 !     INTRINSIC  IABS
                 INTRINSIC  IABS
                 INTEGER    IABS
                 !                                  SPECIFICATIONS FOR SUBROUTINES
                 EXTERNAL   E1INIT, E1MES, E1STI
                 !
                 DATA IFINIT/0/
                 DATA DEFLTP/0, 0, 1, 1, 1, 1, 1/, DEFLTS/0, 0, 0, 1, 1, 0, 1/
                 !                                  INITIALIZE ERROR TABLE IF NECESSARY
                 IF (IFINIT .EQ. 0) THEN
                     CALL E1INIT
                     IFINIT = 1
                 END IF
                 IER = 0
                 IF (IERTYP .GE. 0) THEN
                     IF (IPATT.LT.-1 .OR. IPATT.GT.2) THEN
                         CALL E1STI (1, IPATT)
                         CALL E1MES (5, 1, 'Invalid value specified for print '// &
                             'table attribute.  IPATT must be -1, 0, 1, '// &
                             'or 2.  IPATT = %(I1)')
                         IER = 1
                     END IF
                     IF (ISATT.LT.-1 .OR. ISATT.GT.2) THEN
                         CALL E1STI (1, ISATT)
                         CALL E1MES (5, 1, 'Invalid value specified for stop '// &
                             'table attribute.  ISATT must be -1, 0, 1, '// &
                             'or 2.  ISATT = %(I1)')
                         IER = 1
                     END IF
                 END IF
                 IF (IER .EQ. 0) THEN
                     IF (IERTYP .EQ. 0) THEN
                         IF (IPATT.EQ.0 .OR. IPATT.EQ.1) THEN
                             DO 10  I=1, 7
10                               PRINTB(I) = IPATT
                             ELSE IF (IPATT .EQ. 2) THEN
                                 !                                  ASSIGN DEFAULT SETTINGS
                                 DO 20  I=1, 7
20                                   PRINTB(I) = DEFLTP(I)
                                 END IF
                                 IF (ISATT.EQ.0 .OR. ISATT.EQ.1) THEN
                                     DO 30  I=1, 7
30                                       STOPTB(I) = ISATT
                                     ELSE IF (ISATT .EQ. 2) THEN
                                         !                                  ASSIGN DEFAULT SETTINGS
                                         DO 40  I=1, 7
40                                           STOPTB(I) = DEFLTS(I)
                                         END IF
                                     ELSE IF (IERTYP.GE.1 .AND. IERTYP.LE.7) THEN
                                         IF (IPATT.EQ.0 .OR. IPATT.EQ.1) THEN
                                             PRINTB(IERTYP) = IPATT
                                         ELSE IF (IPATT .EQ. 2) THEN
                                             !                                  ASSIGN DEFAULT SETTING
                                             PRINTB(IERTYP) = DEFLTP(IERTYP)
                                         END IF
                                         IF (ISATT.EQ.0 .OR. ISATT.EQ.1) THEN
                                             STOPTB(IERTYP) = ISATT
                                         ELSE IF (ISATT .EQ. 2) THEN
                                             !                                  ASSIGN DEFAULT SETTING
                                             STOPTB(IERTYP) = DEFLTS(IERTYP)
                                         END IF
                                     ELSE IF (IERTYP.LE.-1 .AND. IERTYP.GE.-7) THEN
                                         I = IABS(IERTYP)
                                         IPATT = PRINTB(I)
                                         ISATT = STOPTB(I)
                                     END IF
                                 END IF
                                 !
                                 RETURN
                             END
                             !-----------------------------------------------------------------------
                             !  IMSL Name:  E1PRT
                             !
                             !  Computer:   pcdsms/SINGLE
                             !
                             !  Revised:    March 14, 1984
                             !
                             !  Purpose:    To print an error message.
                             !
                             !  Usage:      CALL E1PRT
                             !
                             !  Arguments:  None
                             !
                             !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
                             !
                             !  Warranty:   IMSL warrants only that IMSL testing has been applied
                             !              to this code.  No other warranty, expressed or implied,
                             !              is applicable.
                             !
                             !-----------------------------------------------------------------------
                             !
                             SUBROUTINE E1PRT
                                 !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                                 INTEGER    ALL, I, IBEG, IBLOC, IBLOC2, IEND, IER, IHDR, J, &
                                     LERTYP, LOC, LOCM1, LOCX, MAXLOC, MAXTMP, MLOC, MOD, &
                                     NCBEG, NLOC, NOUT
                                 CHARACTER  MSGTMP(70), STRING(10)
                                 !                                  SPECIFICATIONS FOR SAVE VARIABLES
                                 CHARACTER  ATLINE(9), BLANK(1), DBB(3), FROM(6), MSGTYP(8,7), &
                                     PERSLA(2), QMARK, UNKNOW(8)
                                 !                                  SPECIFICATIONS FOR SPECIAL CASES
                                 !                              SPECIFICATIONS FOR COMMON /ERCOM1/
                                 INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                                     PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                                     IALLOC(51), HDRFMT(7), TRACON(7)
                                 COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                                     PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                                     TRACON
                                 SAVE       /ERCOM1/
                                 !                              SPECIFICATIONS FOR COMMON /ERCOM2/
                                 CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
                                 COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
                                 SAVE       /ERCOM2/
                                 !                              SPECIFICATIONS FOR COMMON /ERCOM3/
                                 DOUBLE PRECISION ERCKSM
                                 COMMON     /ERCOM3/ ERCKSM
                                 SAVE       /ERCOM3/
                                 !                              SPECIFICATIONS FOR COMMON /ERCOM4/
                                 LOGICAL    ISUSER(51)
                                 COMMON     /ERCOM4/ ISUSER
                                 SAVE       /ERCOM4/
                                 !                              SPECIFICATIONS FOR COMMON /ERCOM8/
                                 INTEGER    PROLVL, XXLINE(10), XXPLEN(10), ICALOC(10), INALOC(10)
                                 COMMON     /ERCOM8/ PROLVL, XXLINE, XXPLEN, ICALOC, INALOC
                                 SAVE       /ERCOM8/
                                 !                              SPECIFICATIONS FOR COMMON /ERCOM9/
                                 CHARACTER  XXPROC(10)*31
                                 COMMON     /ERCOM9/ XXPROC
                                 SAVE       /ERCOM9/
                                 SAVE       ATLINE, BLANK, DBB, FROM, MSGTYP, PERSLA, QMARK, &
                                     UNKNOW
                                 !                                  SPECIFICATIONS FOR INTRINSICS
                                 !     INTRINSIC  MIN0
                                 INTRINSIC  MIN0
                                 INTEGER    MIN0
                                 !                                  SPECIFICATIONS FOR SUBROUTINES
                                 EXTERNAL   C1TIC, M1VE, UMACH
                                 !                                  SPECIFICATIONS FOR FUNCTIONS
                                 EXTERNAL   I1DX, I1ERIF
                                 INTEGER    I1DX, I1ERIF
                                 !
                                 DATA MSGTYP/'N', 'O', 'T', 'E', ' ', ' ', ' ', ' ', 'A', &
                                     'L', 'E', 'R', 'T', ' ', ' ', ' ', 'W', 'A', 'R', &
                                     'N', 'I', 'N', 'G', ' ', 'F', 'A', 'T', 'A', 'L', &
                                     ' ', ' ', ' ', 'T', 'E', 'R', 'M', 'I', 'N', 'A', &
                                     'L', 'W', 'A', 'R', 'N', 'I', 'N', 'G', ' ', 'F', &
                                     'A', 'T', 'A', 'L', ' ', ' ', ' '/
                                 DATA UNKNOW/'U', 'N', 'K', 'N', 'O', 'W', 'N', ' '/
                                 DATA ATLINE/' ', 'a', 't', ' ', 'l', 'i', 'n', 'e', ' '/
                                 DATA BLANK/' '/, FROM/' ', 'f', 'r', 'o', 'm', ' '/
                                 DATA DBB/'.', ' ', ' '/, PERSLA/'%', '/'/
                                 DATA QMARK/'?'/
                                 !
                                 IF (MSGLEN .LE. 0) RETURN
                                 CALL UMACH (2, NOUT)
                                 MAXTMP = 70
                                 MOD = 0
                                 LERTYP = ERTYPE(CALLVL)
                                 IHDR = HDRFMT(LERTYP)
                                 IF (IHDR .EQ. 3) THEN
                                     IF (XXPROC(PROLVL)(1:1).EQ.QMARK .AND. XXLINE(PROLVL).EQ.0) &
                                         THEN
                                         IHDR = 1
                                     END IF
                                 END IF
                                 IEND = 0
                                 IF (IHDR.EQ.1 .AND. ERTYPE(CALLVL).LE.4) THEN
                                     MSGTMP(1) = BLANK(1)
                                     IEND = 1
                                     !                                  CONVERT ERROR CODE INTO CHAR STRING
                                     CALL C1TIC (ERCODE(CALLVL), STRING, 10, IER)
                                     !                                  LOCATE START OF NON-BLANK CHARACTERS
                                     IBEG = I1ERIF(STRING,10,BLANK,1)
                                     !                                  M1VE IT TO MSGTMP
                                     CALL M1VE (STRING, IBEG, 10, 10, MSGTMP, IEND+1, &
                                         IEND+11-IBEG, MAXTMP, IER)
                                     IEND = IEND + 11 - IBEG
                                 END IF
                                 IF (IHDR .NE. 2) THEN
                                     CALL M1VE (FROM, 1, 6, 6, MSGTMP, IEND+1, IEND+6, MAXTMP, IER)
                                     IEND = IEND + 6
                                 END IF
                                 IF (IHDR .EQ. 3) THEN
                                     !                                  THIS IS A PROTRAN RUN TIME ERROR MSG.
                                     !                                  RETRIEVE THE PROCEDURE NAME
                                     CALL M1VE (XXPROC(PROLVL), 1, XXPLEN(PROLVL), 31, MSGTMP, &
                                         IEND+1, IEND+XXPLEN(PROLVL), MAXTMP, IER)
                                     MLOC = IEND + XXPLEN(PROLVL) + 1
                                     MSGTMP(MLOC) = BLANK(1)
                                     IEND = IEND + I1DX(MSGTMP(IEND+1),XXPLEN(PROLVL)+1,BLANK,1) - &
                                         1
                                     IF (XXLINE(PROLVL) .GT. 0) THEN
                                         !                                  INSERT ATLINE
                                         CALL M1VE (ATLINE, 1, 9, 9, MSGTMP, IEND+1, IEND+9, &
                                             MAXTMP, IER)
                                         IEND = IEND + 9
                                         !                                  CONVERT PROTRAN GLOBAL LINE NUMBER
                                         CALL C1TIC (XXLINE(PROLVL), STRING, 10, IER)
                                         !                                  LOCATE START OF NON-BLANK CHARACTERS
                                         IBEG = I1ERIF(STRING,10,BLANK,1)
                                         !                                  M1VE GLOBAL LINE NUMBER TO MSGTMP
                                         CALL M1VE (STRING, IBEG, 10, 10, MSGTMP, IEND+1, &
                                             IEND+11-IBEG, MAXTMP, IER)
                                         IEND = IEND + 11 - IBEG
                                     END IF
                                 ELSE
                                     !                                  THIS IS EITHER A LIBRARY ERROR MSG
                                     !                                  OR A PROTRAN PREPROCESSOR ERROR MSG
                                     IF (IHDR .EQ. 1) THEN
                                         !                                  THIS IS A LIBRARY ERROR MESSAGE.
                                         !                                  RETRIEVE ROUTINE NAME
                                         CALL M1VE (RNAME(CALLVL), 1, 6, 6, MSGTMP, IEND+1, IEND+6, &
                                             MAXTMP, IER)
                                         MSGTMP(IEND+7) = BLANK(1)
                                         IEND = IEND + I1DX(MSGTMP(IEND+1),7,BLANK,1) - 1
                                     END IF
                                     !                                  ADD DOT, BLANK, BLANK IF NEEDED
                                     IF (I1DX(MSGSAV,3,DBB,3) .NE. 1) THEN
                                         CALL M1VE (DBB, 1, 3, 3, MSGTMP, IEND+1, IEND+3, MAXTMP, &
                                             IER)
                                         IEND = IEND + 3
                                         MOD = 3
                                     END IF
                                 END IF
                                 !                                  MSGTMP AND MSGSAV NOW CONTAIN THE
                                 !                                   ERROR MESSAGE IN FINAL FORM.
                                 NCBEG = 59 - IEND - MOD
                                 ALL = 0
                                 IBLOC = I1DX(MSGSAV,MSGLEN,PERSLA,2)
                                 IF (IBLOC.NE.0 .AND. IBLOC.LT.NCBEG) THEN
                                     LOCM1 = IBLOC - 1
                                     LOC = IBLOC + 1
                                 ELSE IF (MSGLEN .LE. NCBEG) THEN
                                     LOCM1 = MSGLEN
                                     ALL = 1
                                 ELSE
                                     LOC = NCBEG
                                 !                                  CHECK FOR APPROPRIATE PLACE TO SPLIT
10                               CONTINUE
                                 IF (MSGSAV(LOC) .NE. BLANK(1)) THEN
                                     LOC = LOC - 1
                                     IF (LOC .GT. 1) GO TO 10
                                     LOC = NCBEG + 1
                                 END IF
                                 LOCM1 = LOC - 1
                             END IF
                             !                                  NO BLANKS FOUND IN FIRST NCBEG CHARS
                             IF (LERTYP.GE.1 .AND. LERTYP.LE.7) THEN
                                 WRITE (NOUT,99995) (MSGTYP(I,LERTYP),I=1,8), &
                                     (MSGTMP(I),I=1,IEND), (MSGSAV(I),I=1,LOCM1)
                             ELSE
                                 WRITE (NOUT,99995) (UNKNOW(I),I=1,8), (MSGTMP(I),I=1,IEND), &
                                     (MSGSAV(I),I=1,LOCM1)
                             END IF
                             IF (ALL .EQ. 0) THEN
                                 !                                  PREPARE TO WRITE CONTINUATION OF
                                 !                                    MESSAGE
                                 !
                                 !                                  FIND WHERE TO BREAK MESSAGE
                                 !                                    LOC = NUMBER OF CHARACTERS OF
                                 !                                          MESSAGE WRITTEN SO FAR
20                               LOCX = LOC + 64
                                 NLOC = LOC + 1
                                 IBLOC2 = IBLOC
                                 MAXLOC = MIN0(MSGLEN-LOC,64)
                                 IBLOC = I1DX(MSGSAV(NLOC),MAXLOC,PERSLA,2)
                                 IF (MSGSAV(NLOC).EQ.BLANK(1) .AND. IBLOC2.EQ.0) NLOC = NLOC + &
                                     1
                                 IF (IBLOC .GT. 0) THEN
                                     !                                  PAGE BREAK FOUND AT IBLOC
                                     LOCX = NLOC + IBLOC - 2
                                     WRITE (NOUT,99996) (MSGSAV(I),I=NLOC,LOCX)
                                     LOC = NLOC + IBLOC
                                     GO TO 20
                                 !                                  DON'T BOTHER LOOKING FOR BLANK TO
                                 !                                    BREAK AT IF LOCX .GE. MSGLEN
                                 ELSE IF (LOCX .LT. MSGLEN) THEN
                                 !                                  CHECK FOR BLANK TO BREAK THE LINE
30                               CONTINUE
                                 IF (MSGSAV(LOCX) .EQ. BLANK(1)) THEN
                                     !                                  BLANK FOUND AT LOCX
                                     WRITE (NOUT,99996) (MSGSAV(I),I=NLOC,LOCX)
                                     LOC = LOCX
                                     GO TO 20
                                 END IF
                                 LOCX = LOCX - 1
                                 IF (LOCX .GT. NLOC) GO TO 30
                                 LOCX = LOC + 64
                                 !                                  NO BLANKS FOUND IN NEXT 64 CHARS
                                 WRITE (NOUT,99996) (MSGSAV(I),I=NLOC,LOCX)
                                 LOC = LOCX
                                 GO TO 20
                             ELSE
                                 !                                  ALL THE REST WILL FIT ON 1 LINE
                                 LOCX = MSGLEN
                                 WRITE (NOUT,99996) (MSGSAV(I),I=NLOC,LOCX)
                             END IF
                         END IF
                         !                                  SET LENGTH OF MSGSAV AND PLEN
                         !                                    TO SHOW THAT MESSAGE HAS
                         !                                    ALREADY BEEN PRINTED
9000                     MSGLEN = 0
                         PLEN = 1
                         IF (TRACON(LERTYP).EQ.1 .AND. CALLVL.GT.2) THEN
                             !                                  INITIATE TRACEBACK
                             WRITE (NOUT,99997)
                             DO 9005  J=CALLVL, 1, -1
                                 IF (J .GT. 1) THEN
                                     IF (ISUSER(J-1)) THEN
                                         WRITE (NOUT,99998) RNAME(J), ERTYPE(J), ERCODE(J)
                                     ELSE
                                         WRITE (NOUT,99999) RNAME(J), ERTYPE(J), ERCODE(J)
                                     END IF
                                 ELSE
                                     WRITE (NOUT,99998) RNAME(J), ERTYPE(J), ERCODE(J)
                                 END IF
9005                         CONTINUE
                         END IF
                         !
                         RETURN
99995                    FORMAT (/, ' *** ', 8A1, ' ERROR', 59A1)
99996                    FORMAT (' *** ', 9X, 64A1)
99997                    FORMAT (14X, 'Here is a traceback of subprogram calls', &
                             ' in reverse order:', /, 14X, '      Routine    Error ', &
                             'type    Error code', /, 14X, '      -------    ', &
                             '----------    ----------')
99998                    FORMAT (20X, A6, 5X, I6, 8X, I6)
99999                    FORMAT (20X, A6, 5X, I6, 8X, I6, 4X, '(Called internally)')
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  E1PSH
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    March 2, 1984
                     !
                     !  Purpose:    To push a subroutine name onto the error control stack.
                     !
                     !  Usage:      CALL E1PSH(NAME)
                     !
                     !  Arguments:
                     !     NAME   - A character string of length six specifing the name of
                     !              the subroutine.  (Input)
                     !
                     !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     SUBROUTINE E1PSH (NAME)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         CHARACTER  NAME*(*)
                         !                                  SPECIFICATIONS FOR SAVE VARIABLES
                         INTEGER    IFINIT
                         SAVE       IFINIT
                         !                                  SPECIFICATIONS FOR SPECIAL CASES
                         !                              SPECIFICATIONS FOR COMMON /ERCOM1/
                         INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                             PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                             IALLOC(51), HDRFMT(7), TRACON(7)
                         COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                             PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                             TRACON
                         SAVE       /ERCOM1/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM2/
                         CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
                         COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
                         SAVE       /ERCOM2/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM3/
                         DOUBLE PRECISION ERCKSM
                         COMMON     /ERCOM3/ ERCKSM
                         SAVE       /ERCOM3/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM4/
                         LOGICAL    ISUSER(51)
                         COMMON     /ERCOM4/ ISUSER
                         SAVE       /ERCOM4/
                         !                                  SPECIFICATIONS FOR SUBROUTINES
                         EXTERNAL   E1INIT, E1MES, E1STI
                         !                                  SPECIFICATIONS FOR FUNCTIONS
                         EXTERNAL   I1KST
                         INTEGER    I1KST
                         !
                         DATA IFINIT/0/
                         !                                  INITIALIZE ERROR TABLE IF NECESSARY
                         IF (IFINIT .EQ. 0) THEN
                             CALL E1INIT
                             IFINIT = 1
                         END IF
                         IF (CALLVL .GE. MAXLEV) THEN
                             CALL E1STI (1, MAXLEV)
                             CALL E1MES (5, 1, 'Error condition in E1PSH.  Push would '// &
                                 'cause stack level to exceed %(I1). ')
                             STOP
                         ELSE
                             !                                  STORE ALLOCATION LEVEL
                             IALLOC(CALLVL) = I1KST(1)
                             !                                  INCREMENT THE STACK POINTER BY ONE
                             CALLVL = CALLVL + 1
                             !                                  PUT SUBROUTINE NAME INTO STACK
                             RNAME(CALLVL) = NAME
                             !                                  SET ERROR TYPE AND ERROR CODE
                             ERTYPE(CALLVL) = 0
                             ERCODE(CALLVL) = 0
                         END IF
                         !
                         RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  E1STI
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    March 6, 1984
                     !
                     !  Purpose:    To store an integer for subsequent use within an error
                     !              message.
                     !
                     !  Usage:      CALL E1STI(II, IVALUE)
                     !
                     !  Arguments:
                     !     II     - Integer specifying the substitution index.  II must be
                     !              between 1 and 9.  (Input)
                     !     IVALUE - The integer to be stored.  (Input)
                     !
                     !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     SUBROUTINE E1STI (II, IVALUE)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         INTEGER    II, IVALUE
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    IBEG, IER, ILEN
                         CHARACTER  ARRAY(14)
                         !                                  SPECIFICATIONS FOR SAVE VARIABLES
                         INTEGER    IFINIT
                         CHARACTER  BLANK(1)
                         SAVE       BLANK, IFINIT
                         !                                  SPECIFICATIONS FOR SPECIAL CASES
                         !                              SPECIFICATIONS FOR COMMON /ERCOM1/
                         INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                             PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                             IALLOC(51), HDRFMT(7), TRACON(7)
                         COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                             PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                             TRACON
                         SAVE       /ERCOM1/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM2/
                         CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
                         COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
                         SAVE       /ERCOM2/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM3/
                         DOUBLE PRECISION ERCKSM
                         COMMON     /ERCOM3/ ERCKSM
                         SAVE       /ERCOM3/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM4/
                         LOGICAL    ISUSER(51)
                         COMMON     /ERCOM4/ ISUSER
                         SAVE       /ERCOM4/
                         !                                  SPECIFICATIONS FOR SUBROUTINES
                         EXTERNAL   C1TIC, E1INIT, E1INPL
                         !                                  SPECIFICATIONS FOR FUNCTIONS
                         EXTERNAL   I1ERIF
                         INTEGER    I1ERIF
                         !
                         DATA BLANK/' '/, IFINIT/0/
                         !                                  INITIALIZE IF NECESSARY
                         IF (IFINIT .EQ. 0) THEN
                             CALL E1INIT
                             IFINIT = 1
                         END IF
                         CALL C1TIC (IVALUE, ARRAY, 14, IER)
                         IBEG = I1ERIF(ARRAY,14,BLANK,1)
                         IF (II.GE.1 .AND. II.LE.9 .AND. IER.EQ.0) THEN
                             ILEN = 15 - IBEG
                             CALL E1INPL ('I', II, ILEN, ARRAY(IBEG))
                         END IF
                         !
                         RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  E1STL
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    November 8, 1985
                     !
                     !  Purpose:    To store a string for subsequent use within an error
                     !              message.
                     !
                     !  Usage:      CALL E1STL(IL,STRING)
                     !
                     !  Arguments:
                     !     IL     - Integer specifying the substitution index.  IL must be
                     !              between 1 and 9.  (Input)
                     !     STRING - A character string.  (Input)
                     !
                     !  Copyright:  1985 by IMSL, Inc.  All rights reserved.
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     SUBROUTINE E1STL (IL, STRING)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         INTEGER    IL
                         CHARACTER  STRING*(*)
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    I, LEN2
                         CHARACTER  STRGUP(255)
                         !                                  SPECIFICATIONS FOR SAVE VARIABLES
                         INTEGER    IFINIT
                         SAVE       IFINIT
                         !                                  SPECIFICATIONS FOR SPECIAL CASES
                         !                              SPECIFICATIONS FOR COMMON /ERCOM1/
                         INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                             PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                             IALLOC(51), HDRFMT(7), TRACON(7)
                         COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                             PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                             TRACON
                         SAVE       /ERCOM1/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM2/
                         CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
                         COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
                         SAVE       /ERCOM2/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM3/
                         DOUBLE PRECISION ERCKSM
                         COMMON     /ERCOM3/ ERCKSM
                         SAVE       /ERCOM3/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM4/
                         LOGICAL    ISUSER(51)
                         COMMON     /ERCOM4/ ISUSER
                         SAVE       /ERCOM4/
                         !                                  SPECIFICATIONS FOR INTRINSICS
                         !     INTRINSIC  IABS,LEN,MIN0
                         INTRINSIC  IABS, LEN, MIN0
                         INTEGER    IABS, LEN, MIN0
                         !                                  SPECIFICATIONS FOR SUBROUTINES
                         EXTERNAL   E1INIT, E1INPL
                         !
                         DATA IFINIT/0/
                         !                                  INITIALIZE IF NECESSARY
                         IF (IFINIT .EQ. 0) THEN
                             CALL E1INIT
                             IFINIT = 1
                         END IF
                         LEN2 = LEN(STRING)
                         LEN2 = MIN0(LEN2,255)
                         DO 10  I=1, LEN2
                             STRGUP(I) = STRING(I:I)
10                       CONTINUE
                         IF (IABS(IL).GE.1 .AND. IABS(IL).LE.9) THEN
                             CALL E1INPL ('L', IL, LEN2, STRGUP)
                         END IF
                         !
                         RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  E1UCS
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    March 8, 1984
                     !
                     !  Purpose:    To update the checksum number for error messages.
                     !
                     !  Usage:      CALL E1UCS
                     !
                     !  Arguments:  None
                     !
                     !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     SUBROUTINE E1UCS
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    I, IBEG, IBEG2, IEND, ILOC, IPOS, JLOC, NCODE, NLEN
                         DOUBLE PRECISION DNUM
                         !                                  SPECIFICATIONS FOR SAVE VARIABLES
                         DOUBLE PRECISION DMAX
                         CHARACTER  BLANK(1), COMMA(1), EQUAL(1), LPAR(1)
                         SAVE       BLANK, COMMA, DMAX, EQUAL, LPAR
                         !                                  SPECIFICATIONS FOR SPECIAL CASES
                         !                              SPECIFICATIONS FOR COMMON /ERCOM1/
                         INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                             PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                             IALLOC(51), HDRFMT(7), TRACON(7)
                         COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                             PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                             TRACON
                         SAVE       /ERCOM1/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM2/
                         CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
                         COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
                         SAVE       /ERCOM2/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM3/
                         DOUBLE PRECISION ERCKSM
                         COMMON     /ERCOM3/ ERCKSM
                         SAVE       /ERCOM3/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM4/
                         LOGICAL    ISUSER(51)
                         COMMON     /ERCOM4/ ISUSER
                         SAVE       /ERCOM4/
                         !                                  SPECIFICATIONS FOR INTRINSICS
                         !     INTRINSIC  DMOD
                         INTRINSIC  DMOD
                         DOUBLE PRECISION DMOD
                         !                                  SPECIFICATIONS FOR SUBROUTINES
                         EXTERNAL   S1ANUM
                         !                                  SPECIFICATIONS FOR FUNCTIONS
                         EXTERNAL   ICASE, I1X
                         INTEGER    ICASE, I1X
                         !
                         DATA BLANK(1)/' '/, COMMA(1)/','/, LPAR(1)/'('/
                         DATA EQUAL(1)/'='/, DMAX/1.0D+9/
                         !
                         IF (MSGLEN .GT. 1) THEN
                             IPOS = 0
                             IBEG2 = 1
10                           IBEG = IBEG2
                             IEND = MSGLEN
                             !                                  LOOK FOR BLANK, COMMA, LEFT PAREN.,
                             !                                  OR EQUAL SIGN
                             ILOC = I1X(MSGSAV(IBEG),IEND-IBEG+1,BLANK,1)
                             JLOC = I1X(MSGSAV(IBEG),IEND-IBEG+1,COMMA,1)
                             IF (ILOC.EQ.0 .OR. (JLOC.GT.0.AND.JLOC.LT.ILOC)) ILOC = JLOC
                             JLOC = I1X(MSGSAV(IBEG),IEND-IBEG+1,LPAR,1)
                             IF (ILOC.EQ.0 .OR. (JLOC.GT.0.AND.JLOC.LT.ILOC)) ILOC = JLOC
                             JLOC = I1X(MSGSAV(IBEG),IEND-IBEG+1,EQUAL,1)
                             IF (ILOC.EQ.0 .OR. (JLOC.GT.0.AND.JLOC.LT.ILOC)) ILOC = JLOC
                             IF (ILOC .GE. 1) THEN
                                 CALL S1ANUM (MSGSAV(IBEG+ILOC), IEND-IBEG-ILOC+1, NCODE, &
                                     NLEN)
                                 IF (NCODE.EQ.2 .OR. NCODE.EQ.3) THEN
                                     !                                  FLOATING POINT NUMBER FOUND.
                                     !                                  SET POINTERS TO SKIP OVER IT
                                     IBEG2 = IBEG + ILOC + NLEN
                                     IF (IBEG2 .LE. MSGLEN) THEN
                                         CALL S1ANUM (MSGSAV(IBEG2), IEND-IBEG2+1, NCODE, &
                                             NLEN)
                                         IF ((MSGSAV(IBEG2).EQ.'+'.OR.MSGSAV(IBEG2).EQ. &
                                             '-') .AND. NCODE.EQ.1) THEN
                                             !                                  INTEGER IMMEDIATELY FOLLOWS A REAL AS
                                             !                                  WITH SOME CDC NOS. LIKE 1.2345678+123
                                             !                                  SET POINTERS TO SKIP OVER IT
                                             IBEG2 = IBEG2 + NLEN
                                         END IF
                                     END IF
                                 ELSE
                                     IBEG2 = IBEG + ILOC
                                 END IF
                                 IEND = IBEG + ILOC - 1
                             END IF
                             !                                  UPDATE CKSUM USING PART OF MESSAGE
                             DO 20  I=IBEG, IEND
                                 IPOS = IPOS + 1
                                 DNUM = ICASE(MSGSAV(I))
                                 ERCKSM = DMOD(ERCKSM+DNUM*IPOS,DMAX)
20                           CONTINUE
                             !                                  GO BACK FOR MORE IF NEEDED
                             IF (IEND.LT.MSGLEN .AND. IBEG2.LT.MSGLEN) GO TO 10
                             !                                  UPDATE CKSUM USING ERROR TYPE
                             DNUM = ERTYPE(CALLVL)
                             ERCKSM = DMOD(ERCKSM+DNUM*(IPOS+1),DMAX)
                             !                                  UPDATE CKSUM USING ERROR CODE
                             DNUM = ERCODE(CALLVL)
                             ERCKSM = DMOD(ERCKSM+DNUM*(IPOS+2),DMAX)
                         END IF
                         !
                         RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  I1CSTR (Single precision version)
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    September 10, 1985
                     !
                     !  Purpose:    Case insensitive comparison of two character arrays.
                     !
                     !  Usage:      I1CSTR(STR1, LEN1, STR2, LEN2)
                     !
                     !  Arguments:
                     !     STR1   - First character array.  (Input)
                     !     LEN1   - Length of STR1.  (Input)
                     !     STR2   - Second character array.  (Input)
                     !     LEN2   - Length of STR2.  (Input)
                     !     I1CSTR - Integer function.  (Output) Where
                     !              I1CSTR = -1  if STR1 .LT. STR2,
                     !              I1CSTR =  0  if STR1 .EQ. STR2,
                     !              I1CSTR =  1  if STR1 .GT. STR2.
                     !
                     !  Remarks:
                     !  1. If the two arrays, STR1 and STR2,  are of unequal length, the
                     !     shorter array is considered as if it were extended with blanks
                     !     to the length of the longer array.
                     !
                     !  2. If one or both lengths are zero or negative the I1CSTR output is
                     !     based on comparison of the lengths.
                     !
                     !  GAMS:       N5c
                     !
                     !  Chapter:    MATH/LIBRARY Utilities
                     !
                     !  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     INTEGER FUNCTION I1CSTR (STR1, LEN1, STR2, LEN2)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         INTEGER    LEN1, LEN2
                         CHARACTER  STR1(LEN1), STR2(LEN2)
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    IC1, IC2, ICB, IS, L, LENM
                         !                                  SPECIFICATIONS FOR INTRINSICS
                         !     INTRINSIC  ISIGN,MIN0
                         INTRINSIC  ISIGN, MIN0
                         INTEGER    ISIGN, MIN0
                         !                                  SPECIFICATIONS FOR FUNCTIONS
                         EXTERNAL   ICASE
                         INTEGER    ICASE
                         !
                         IF (LEN1.GT.0 .AND. LEN2.GT.0) THEN
                             !                                  COMPARE FIRST LENM CHARACTERS
                             LENM = MIN0(LEN1,LEN2)
                             DO 10  L=1, LENM
                                 IC1 = ICASE(STR1(L))
                                 IC2 = ICASE(STR2(L))
                                 IF (IC1 .NE. IC2) THEN
                                     I1CSTR = ISIGN(1,IC1-IC2)
                                     RETURN
                                 END IF
10                           CONTINUE
                         END IF
                         !                                  COMPARISON BASED ON LENGTH OR
                         !                                  TRAILING BLANKS
                         IS = LEN1 - LEN2
                         IF (IS .EQ. 0) THEN
                             I1CSTR = 0
                         ELSE
                             IF (LEN1.LE.0 .OR. LEN2.LE.0) THEN
                                 !                                  COMPARISON BASED ON LENGTH
                                 I1CSTR = ISIGN(1,IS)
                             ELSE
                                 !                                  COMPARISON BASED ON TRAILING BLANKS
                                 !                                  TO EXTEND SHORTER ARRAY
                                 LENM = LENM + 1
                                 ICB = ICASE(' ')
                                 IF (IS .GT. 0) THEN
                                     !                                  EXTEND STR2 WITH BLANKS
                                     DO 20  L=LENM, LEN1
                                         IC1 = ICASE(STR1(L))
                                         IF (IC1 .NE. ICB) THEN
                                             I1CSTR = ISIGN(1,IC1-ICB)
                                             RETURN
                                         END IF
20                                   CONTINUE
                                 ELSE
                                     !                                  EXTEND STR1 WITH BLANKS
                                     DO 30  L=LENM, LEN2
                                         IC2 = ICASE(STR2(L))
                                         IF (ICB .NE. IC2) THEN
                                             I1CSTR = ISIGN(1,ICB-IC2)
                                             RETURN
                                         END IF
30                                   CONTINUE
                                 END IF
                                 !
                                 I1CSTR = 0
                             END IF
                         END IF
                         !
                         RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  I1DX (Single precision version)
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    September 9, 1985
                     !
                     !  Purpose:    Determine the array subscript indicating the starting
                     !              element at which a key character sequence begins.
                     !              (Case-insensitive version)
                     !
                     !  Usage:      I1DX(CHRSTR, I1LEN, KEY, KLEN)
                     !
                     !  Arguments:
                     !     CHRSTR - Character array to be searched.  (Input)
                     !     I1LEN  - Length of CHRSTR.  (Input)
                     !     KEY    - Character array that contains the key sequence.  (Input)
                     !     KLEN   - Length of KEY.  (Input)
                     !     I1DX   - Integer function.  (Output)
                     !
                     !  Remarks:
                     !  1. Returns zero when there is no match.
                     !
                     !  2. Returns zero if KLEN is longer than ISLEN.
                     !
                     !  3. Returns zero when any of the character arrays has a negative or
                     !     zero length.
                     !
                     !  GAMS:       N5c
                     !
                     !  Chapter:    MATH/LIBRARY Utilities
                     !
                     !  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     INTEGER FUNCTION I1DX (CHRSTR, I1LEN, KEY, KLEN)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         INTEGER    I1LEN, KLEN
                         CHARACTER  CHRSTR(*), KEY(*)
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    I, II, J
                         !                                  SPECIFICATIONS FOR FUNCTIONS
                         EXTERNAL   ICASE, I1CSTR
                         INTEGER    ICASE, I1CSTR
                         !
                         I1DX = 0
                         IF (KLEN.LE.0 .OR. I1LEN.LE.0) GO TO 9000
                         IF (KLEN .GT. I1LEN) GO TO 9000
                         !
                         I = 1
                         II = I1LEN - KLEN + 1
10                       IF (I .LE. II) THEN
                             IF (ICASE(CHRSTR(I)) .EQ. ICASE(KEY(1))) THEN
                                 IF (KLEN .NE. 1) THEN
                                     J = KLEN - 1
                                     IF (I1CSTR(CHRSTR(I+1),J,KEY(2),J) .EQ. 0) THEN
                                         I1DX = I
                                         GO TO 9000
                                     END IF
                                 ELSE
                                     I1DX = I
                                     GO TO 9000
                                 END IF
                             END IF
                             I = I + 1
                             GO TO 10
                         END IF
                         !
9000                     RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  I1ERIF
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    March 13, 1984
                     !
                     !  Purpose:    Return the position of the first element of a given
                     !              character array which is not an element of another
                     !              character array.
                     !
                     !  Usage:      I1ERIF(STR1, LEN1, STR2, LEN2)
                     !
                     !  Arguments:
                     !     STR1   - Character array to be searched.  (Input)
                     !     LEN1   - Length of STR1.  (Input)
                     !     STR2   - Character array to be searched for.  (Input)
                     !     LEN2   - Length of STR2.  (Input)
                     !     I1ERIF - Integer function.  (Output)
                     !
                     !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     INTEGER FUNCTION I1ERIF (STR1, LEN1, STR2, LEN2)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         INTEGER    LEN1, LEN2
                         CHARACTER  STR1(*), STR2(*)
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    I
                         !                                  SPECIFICATIONS FOR FUNCTIONS
                         EXTERNAL   I1X
                         INTEGER    I1X
                         !                              FIRST EXECUTABLE STATEMENT
                         IF (LEN1.LE.0 .OR. LEN2.LE.0) THEN
                             I1ERIF = 1
                         ELSE
                             DO 10  I=1, LEN1
                                 IF (I1X(STR2,LEN2,STR1(I),1) .EQ. 0) THEN
                                     I1ERIF = I
                                     RETURN
                                 END IF
10                           CONTINUE
                             I1ERIF = 0
                         END IF
                         !
                         RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  I1KGT
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    January 17, 1984
                     !
                     !  Purpose:    Allocate numerical workspace.
                     !
                     !  Usage:      I1KGT(NELMTS,ITYPE)
                     !
                     !  Arguments:
                     !     NELMTS - Number of elements of data type ITYPE to be
                     !              allocated.  (Input)
                     !     ITYPE  - Data type of array to be allocated.  (Input)
                     !                 1 - logical
                     !                 2 - integer
                     !                 3 - real
                     !                 4 - double precision
                     !                 5 - complex
                     !                 6 - double complex
                     !     I1KGT  - Integer function.  (Output)  Returns the index of the
                     !              first element in the current allocation.
                     !
                     !  Remarks:
                     !  1. On return, the array will occupy
                     !     WKSP(I1KGT), WKSP(I1KGT+1), ..., WKSP(I1KGT+NELMTS-1) where
                     !     WKSP is an array of data type ITYPE equivalenced to RWKSP.
                     !
                     !  2. If I1KGT is negative, the absolute value of I1KGT is the
                     !     additional workspace needed for the current allocation.
                     !
                     !  3. The allocator reserves the first sixteen integer locations of
                     !     the stack for its own internal bookkeeping.  These are initialized
                     !     by the function IWKIN upon the first call to the allocation
                     !     package.
                     !
                     !  4. The use of the first ten integer locations is as follows:
                     !      WKSP( 1) - LOUT    The number of current allocations
                     !      WKSP( 2) - LNOW    The current active length of the stack
                     !      WKSP( 3) - LUSED   The maximum value of WKSP(2) achieved
                     !                         thus far
                     !      WKSP( 4) - LBND    The lower bound of permanent storage which
                     !                         is one numeric storage unit more than the
                     !                         maximum allowed length of the stack.
                     !      WKSP( 5) - LMAX    The maximum length of the storage array
                     !      WKSP( 6) - LALC    The total number of allocations handled by
                     !                         I1KGT
                     !      WKSP( 7) - LNEED   The number of numeric storage units by which
                     !                         the array size must be increased for all past
                     !                         allocations to succeed
                     !      WKSP( 8) - LBOOK   The number of numeric storage units used for
                     !                         bookkeeping
                     !      WKSP( 9) - LCHAR   The pointer to the portion of the permanent
                     !                         stack which contains the bookkeeping and
                     !                         pointers for the character workspace
                     !                         allocation.
                     !      WKSP(10) - LLCHAR  The length of the array beginning at LCHAR
                     !                         set aside for character workspace bookkeeping
                     !                         and pointers.
                     !                 NOTE -  If character workspace is not being used,
                     !                         LCHAR and LLCHAR can be ignored.
                     !  5. The next six integer locations contain values describing the
                     !     amount of storage allocated by the allocation system to the
                     !     various data types.
                     !      WKSP(11) - Numeric storage units allocated to LOGICAL
                     !      WKSP(12) - Numeric storage units allocated to INTEGER
                     !      WKSP(13) - Numeric storage units allocated to REAL
                     !      WKSP(14) - Numeric storage units allocated to DOUBLE PRECISION
                     !      WKSP(15) - Numeric storage units allocated to COMPLEX
                     !      WKSP(16) - Numeric storage units allocated to DOUBLE COMPLEX
                     !
                     !  Copyright:  1984 by IMSL, Inc. All Rights Reserved.
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     INTEGER FUNCTION I1KGT (NELMTS, ITYPE)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         INTEGER    NELMTS, ITYPE
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    I, IDUMAL, IGAP, ILEFT, IPA, IPA7, ISA, ISA7, &
                             ISIZE(6), JTYPE, LALC, LBND, LBOOK, LMAX, LNEED, &
                             LNEED1, LNOW, LOUT, LUSED
                         !                                  SPECIFICATIONS FOR SAVE VARIABLES
                         LOGICAL    FIRST
                         SAVE       FIRST
                         !                                  SPECIFICATIONS FOR SPECIAL CASES
                         !                              SPECIFICATIONS FOR COMMON /ERCOM8/
                         INTEGER    PROLVL, XXLINE(10), XXPLEN(10), ICALOC(10), INALOC(10)
                         COMMON     /ERCOM8/ PROLVL, XXLINE, XXPLEN, ICALOC, INALOC
                         SAVE       /ERCOM8/
                         !                              SPECIFICATIONS FOR COMMON /ERCOM9/
                         CHARACTER  XXPROC(10)*31
                         COMMON     /ERCOM9/ XXPROC
                         SAVE       /ERCOM9/
                         !                                  SPECIFICATIONS FOR COMMON /WORKSP/
                         REAL       RWKSP(5000)
                         REAL       RDWKSP(5000)
                         DOUBLE PRECISION DWKSP(2500)
                         COMPLEX    CWKSP(2500)
                         COMPLEX    CZWKSP(2500)
                         COMPLEX    *16 ZWKSP(1250)
                         INTEGER    IWKSP(5000)
                         LOGICAL    LWKSP(5000)
                         EQUIVALENCE (DWKSP(1), RWKSP(1))
                         EQUIVALENCE (CWKSP(1), RWKSP(1)), (ZWKSP(1), RWKSP(1))
                         EQUIVALENCE (IWKSP(1), RWKSP(1)), (LWKSP(1), RWKSP(1))
                         EQUIVALENCE (RDWKSP(1), RWKSP(1)), (CZWKSP(1), RWKSP(1))
                         COMMON     /WORKSP/ RWKSP
                         !                                  SPECIFICATIONS FOR EQUIVALENCE
                         EQUIVALENCE (LOUT, IWKSP(1))
                         EQUIVALENCE (LNOW, IWKSP(2))
                         EQUIVALENCE (LUSED, IWKSP(3))
                         EQUIVALENCE (LBND, IWKSP(4))
                         EQUIVALENCE (LMAX, IWKSP(5))
                         EQUIVALENCE (LALC, IWKSP(6))
                         EQUIVALENCE (LNEED, IWKSP(7))
                         EQUIVALENCE (LBOOK, IWKSP(8))
                         EQUIVALENCE (ISIZE(1), IWKSP(11))
                         !                                  SPECIFICATIONS FOR INTRINSICS
                         !     INTRINSIC  IABS,MAX0,MOD
                         INTRINSIC  IABS, MAX0, MOD
                         INTEGER    IABS, MAX0, MOD
                         !                                  SPECIFICATIONS FOR SUBROUTINES
                         EXTERNAL   E1MES, E1POP, E1POS, E1PSH, E1STI, IWKIN
                         !                                  SPECIFICATIONS FOR FUNCTIONS
                         EXTERNAL   I1KQU
                         INTEGER    I1KQU
                         !
                         DATA FIRST/.TRUE./
                         !
                         CALL E1PSH ('I1KGT ')
                         !
                         IF (FIRST) THEN
                             !                                  INITIALIZE WORKSPACE IF NEEDED
                             FIRST = .FALSE.
                             CALL IWKIN (0)
                         END IF
                         !                                  NUMBER OF ELEMENTS LESS THAN 0
                         IF (NELMTS .LT. 0) THEN
                             CALL E1STI (1, NELMTS)
                             CALL E1MES (5, 2, 'Number of elements is not positive.%/'// &
                                 'NELMTS = %(I1).')
                             CALL E1POP ('I1KGT ')
                             GO TO 9000
                         END IF
                         !                                  ILLEGAL DATA TYPE REQUESTED
                         IF (ITYPE.EQ.0 .OR. IABS(ITYPE).GE.7) THEN
                             CALL E1MES (5, 3, 'Illegal data type requested.')
                             CALL E1POP ('I1KGT ')
                             GO TO 9000
                         END IF
                         !                                  BOOKKEEPING OVERWRITTEN
                         IF (LNOW.LT.LBOOK .OR. LNOW.GT.LUSED .OR. LUSED.GT.LMAX .OR. &
                             LNOW.GE.LBND .OR. LOUT.GT.LALC) THEN
                             CALL E1MES (5, 4, 'One or more of the first eight '// &
                                 'bookkeeping locations in IWKSP have been '// &
                                 'overwritten.')
                             CALL E1POP ('I1KGT ')
                             GO TO 9000
                         END IF
                         !
                         CALL E1POP ('I1KGT ')
                         !                                  DETERMINE NUMBER OF LOCATIONS STILL
                         !                                  AVAILABLE FOR DATA TYPE ITYPE
                         !                                  NOTE: I1KQU ALLOWS FOR 2 INTEGER
                         !                                        POINTERS WHICH MUST BE HANDLED
                         !                                        ARTIFICIALLY IF ILEFT = 0.
                         ILEFT = I1KQU(IABS(ITYPE))
                         !
                         IF (ITYPE .GT. 0) THEN
                             !                                  RELEASABLE STORAGE
                             IF (ILEFT .GE. NELMTS) THEN
                                 I1KGT = (LNOW*ISIZE(2)-1)/ISIZE(ITYPE) + 2
                                 I = ((I1KGT-1+NELMTS)*ISIZE(ITYPE)-1)/ISIZE(2) + 3
                                 !                                  IWKSP(I-1) CONTAINS THE DATA TYPE FOR
                                 !                                  THIS ALLOCATION. IWKSP(I) CONTAINS
                                 !                                  LNOW FOR THE PREVIOUS ALLOCATION.
                                 IWKSP(I-1) = ITYPE
                                 IWKSP(I) = LNOW
                                 LOUT = LOUT + 1
                                 LALC = LALC + 1
                                 LNOW = I
                                 LUSED = MAX0(LUSED,LNOW)
                                 LNEED = 0
                             ELSE
                                 !                                  RELEASABLE STORAGE WAS REQUESTED
                                 !                                  BUT THE STACK WOULD OVERFLOW.
                                 !                                  THEREFORE, ALLOCATE RELEASABLE
                                 !                                  SPACE THROUGH THE END OF THE STACK
                                 IF (LNEED .EQ. 0) THEN
                                     IDUMAL = (LNOW*ISIZE(2)-1)/ISIZE(ITYPE) + 2
                                     I = ((IDUMAL-1+ILEFT)*ISIZE(ITYPE)-1)/ISIZE(2) + 3
                                     !                                  ADVANCE COUNTERS AND STORE POINTERS
                                     !                                  IF THERE IS ROOM TO DO SO
                                     IF (I .LT. LBND) THEN
                                         !                                  IWKSP(I-1) CONTAINS THE DATA TYPE FOR
                                         !                                  THIS ALLOCATION. IWKSP(I) CONTAINS
                                         !                                  LNOW FOR THE PREVIOUS ALLOCATION.
                                         IWKSP(I-1) = ITYPE
                                         IWKSP(I) = LNOW
                                         LOUT = LOUT + 1
                                         LALC = LALC + 1
                                         LNOW = I
                                         LUSED = MAX0(LUSED,LNOW)
                                     END IF
                                 END IF
                                 !                                  CALCULATE AMOUNT NEEDED TO ACCOMODATE
                                 !                                  THIS ALLOCATION REQUEST
                                 LNEED1 = (NELMTS-ILEFT)*ISIZE(ITYPE)
                                 IF (ILEFT .EQ. 0) THEN
                                     IGAP = ISIZE(ITYPE) - MOD(LNOW+LNEED,ISIZE(ITYPE))
                                     IF (IGAP .EQ. ISIZE(ITYPE)) IGAP = 0
                                     LNEED1 = LNEED1 + 2*ISIZE(2) + IGAP
                                 END IF
                                 !                                  MODIFY LNEED ACCORDING TO THE SIZE
                                 !                                  OF THE BASE BEING USED (D.P. HERE)
                                 LNEED = LNEED + ((LNEED1+ISIZE(3)-1)/ISIZE(3))
                                 !                                  SINCE CURRENT ALLOCATION IS ILLEGAL,
                                 !                                  RETURN THE NEGATIVE OF THE ADDITIONAL
                                 !                                  AMOUNT NEEDED TO MAKE IT LEGAL
                                 I1KGT = -LNEED
                             END IF
                         ELSE
                             !                                  PERMANENT STORAGE
                             IF (ILEFT .GE. NELMTS) THEN
                                 JTYPE = -ITYPE
                                 I1KGT = (LBND*ISIZE(2)-1)/ISIZE(JTYPE) + 1 - NELMTS
                                 I = ((I1KGT-1)*ISIZE(JTYPE))/ISIZE(2) - 1
                                 !                                  IWKSP(I) CONTAINS LBND FOR PREVIOUS
                                 !                                  PERMANENT STORAGE ALLOCATION.
                                 !                                  IWKSP(I+1) CONTAINS THE DATA TYPE FOR
                                 !                                  THIS ALLOCATION.
                                 IWKSP(I) = LBND
                                 IWKSP(I+1) = JTYPE
                                 LALC = LALC + 1
                                 LBND = I
                                 LNEED = 0
                             ELSE
                                 !                                  PERMANENT STORAGE WAS REQUESTED
                                 !                                  BUT THE STACK WOULD OVERFLOW,
                                 !                                  THEREFORE, ALLOCATE RELEASABLE
                                 !                                  SPACE THROUGH THE END OF THE STACK
                                 IF (LNEED .EQ. 0) THEN
                                     JTYPE = -ITYPE
                                     IDUMAL = (LNOW*ISIZE(2)-1)/ISIZE(JTYPE) + 2
                                     I = ((IDUMAL-1+ILEFT)*ISIZE(JTYPE)-1)/ISIZE(2) + 3
                                     !                                  ADVANCE COUNTERS AND STORE POINTERS
                                     !                                  IF THERE IS ROOM TO DO SO
                                     IF (I .LT. LBND) THEN
                                         !                                  IWKSP(I-1) CONTAINS THE DATA TYPE FOR
                                         !                                  THIS ALLOCATION. IWKSP(I) CONTAINS
                                         !                                  LNOW FOR THE PREVIOUS ALLOCATION.
                                         IWKSP(I-1) = JTYPE
                                         IWKSP(I) = LNOW
                                         LOUT = LOUT + 1
                                         LALC = LALC + 1
                                         LNOW = I
                                         LUSED = MAX0(LUSED,LNOW)
                                     END IF
                                 END IF
                                 !                                  CALCULATE AMOUNT NEEDED TO ACCOMODATE
                                 !                                  THIS ALLOCATION REQUEST
                                 LNEED1 = (NELMTS-ILEFT)*ISIZE(-ITYPE)
                                 IF (ILEFT .EQ. 0) THEN
                                     IGAP = ISIZE(-ITYPE) - MOD(LNOW+LNEED,ISIZE(-ITYPE))
                                     IF (IGAP .EQ. ISIZE(-ITYPE)) IGAP = 0
                                     LNEED1 = LNEED1 + 2*ISIZE(2) + IGAP
                                 END IF
                                 !                                  MODIFY LNEED ACCORDING TO THE SIZE
                                 !                                  OF THE BASE BEING USED (D.P. HERE)
                                 LNEED = LNEED + ((LNEED1+ISIZE(3)-1)/ISIZE(3))
                                 !                                  SINCE CURRENT ALLOCATION IS ILLEGAL,
                                 !                                  RETURN THE NEGATIVE OF THE ADDITIONAL
                                 !                                  AMOUNT NEEDED TO MAKE IT LEGAL
                                 I1KGT = -LNEED
                             END IF
                         END IF
                         !                                  STACK OVERFLOW - UNRECOVERABLE ERROR
9000                     IF (LNEED .GT. 0) THEN
                             CALL E1POS (-5, IPA, ISA)
                             CALL E1POS (5, 0, 0)
                             CALL E1POS (-7, IPA7, ISA7)
                             CALL E1POS (7, 0, 0)
                             CALL E1PSH ('I1KGT ')
                             CALL E1STI (1, LNEED+(LMAX/ISIZE(3)))
                             IF (XXLINE(PROLVL).GE.1 .AND. XXLINE(PROLVL).LE.999) THEN
                                 CALL E1MES (7, 1, 'Insufficient workspace for current '// &
                                     'allocation(s).  Correct by inserting the '// &
                                     'following PROTRAN line: $OPTIONS;WORKSPACE=%'// &
                                     '(I1)')
                             ELSE
                                 CALL E1MES (5, 5, 'Insufficient workspace for current '// &
                                     'allocation(s). Correct by calling IWKIN '// &
                                     'from main program with the three following '// &
                                     'statements:  (REGARDLESS OF PRECISION)%/'// &
                                     '      COMMON /WORKSP/  RWKSP%/      REAL '// &
                                     'RWKSP(%(I1))%/      CALL IWKIN(%(I1))')
                             END IF
                             CALL E1POP ('I1KGT ')
                             CALL E1POS (5, IPA, ISA)
                             CALL E1POS (7, IPA7, ISA7)
                         END IF
                         !
                         RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  I1KQU
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    January 17, 1984
                     !
                     !  Purpose:    Return number of elements of data type ITYPE that
                     !              remain to be allocated in one request.
                     !
                     !  Usage:      I1KQU(ITYPE)
                     !
                     !  Arguments:
                     !     ITYPE  - Type of storage to be checked (Input)
                     !                 1 - logical
                     !                 2 - integer
                     !                 3 - real
                     !                 4 - double precision
                     !                 5 - complex
                     !                 6 - double complex
                     !     I1KQU  - Integer function. (Output) Returns number of elements
                     !              of data type ITYPE remaining in the stack.
                     !
                     !  Copyright:  1983 by IMSL, Inc.  All Rights Reserved
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     INTEGER FUNCTION I1KQU (ITYPE)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         INTEGER    ITYPE
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    ISIZE(6), LALC, LBND, LBOOK, LMAX, LNEED, LNOW, LOUT, &
                             LUSED
                         !                                  SPECIFICATIONS FOR SAVE VARIABLES
                         LOGICAL    FIRST
                         SAVE       FIRST
                         !                                  SPECIFICATIONS FOR SPECIAL CASES
                         !                                  SPECIFICATIONS FOR COMMON /WORKSP/
                         REAL       RWKSP(5000)
                         REAL       RDWKSP(5000)
                         DOUBLE PRECISION DWKSP(2500)
                         COMPLEX    CWKSP(2500)
                         COMPLEX    CZWKSP(2500)
                         COMPLEX    *16 ZWKSP(1250)
                         INTEGER    IWKSP(5000)
                         LOGICAL    LWKSP(5000)
                         EQUIVALENCE (DWKSP(1), RWKSP(1))
                         EQUIVALENCE (CWKSP(1), RWKSP(1)), (ZWKSP(1), RWKSP(1))
                         EQUIVALENCE (IWKSP(1), RWKSP(1)), (LWKSP(1), RWKSP(1))
                         EQUIVALENCE (RDWKSP(1), RWKSP(1)), (CZWKSP(1), RWKSP(1))
                         COMMON     /WORKSP/ RWKSP
                         !                                  SPECIFICATIONS FOR EQUIVALENCE
                         EQUIVALENCE (LOUT, IWKSP(1))
                         EQUIVALENCE (LNOW, IWKSP(2))
                         EQUIVALENCE (LUSED, IWKSP(3))
                         EQUIVALENCE (LBND, IWKSP(4))
                         EQUIVALENCE (LMAX, IWKSP(5))
                         EQUIVALENCE (LALC, IWKSP(6))
                         EQUIVALENCE (LNEED, IWKSP(7))
                         EQUIVALENCE (LBOOK, IWKSP(8))
                         EQUIVALENCE (ISIZE(1), IWKSP(11))
                         !                                  SPECIFICATIONS FOR INTRINSICS
                         !     INTRINSIC  MAX0
                         INTRINSIC  MAX0
                         INTEGER    MAX0
                         !                                  SPECIFICATIONS FOR SUBROUTINES
                         EXTERNAL   E1MES, E1POP, E1PSH, IWKIN
                         !
                         DATA FIRST/.TRUE./
                         !
                         CALL E1PSH ('I1KQU ')
                         !
                         IF (FIRST) THEN
                             !                                  INITIALIZE WORKSPACE IF NEEDED
                             FIRST = .FALSE.
                             CALL IWKIN (0)
                         END IF
                         !                                  BOOKKEEPING OVERWRITTEN
                         IF (LNOW.LT.LBOOK .OR. LNOW.GT.LUSED .OR. LUSED.GT.LMAX .OR. &
                             LNOW.GE.LBND .OR. LOUT.GT.LALC) THEN
                             CALL E1MES (5, 7, 'One or more of the first eight '// &
                                 'bookkeeping locations in IWKSP have been '// &
                                 'overwritten.')
                         ELSE IF (ITYPE.LE.0 .OR. ITYPE.GE.7) THEN
                             !                                  ILLEGAL DATA TYPE REQUESTED
                             CALL E1MES (5, 8, 'Illegal data type requested.')
                         ELSE
                             !                                  THIS CALCULATION ALLOWS FOR THE
                             !                                  TWO POINTER LOCATIONS IN THE STACK
                             !                                  WHICH ARE ASSIGNED TO EACH ALLOCATION
                             I1KQU = MAX0(((LBND-3)*ISIZE(2))/ISIZE(ITYPE)-(LNOW*ISIZE(2)- &
                                 1)/ISIZE(ITYPE)-1,0)
                         END IF
                         !
                         CALL E1POP ('I1KQU ')
                         RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  I1KRL
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    August 9, 1983
                     !
                     !  Purpose:    Deallocate the last N allocations made in the workspace.
                     !              stack by I1KGT
                     !
                     !  Usage:      CALL I1KRL(N)
                     !
                     !  Arguments:
                     !     N      - Number of allocations to be released top down (Input)
                     !
                     !  Copyright:  1983 by IMSL, Inc.  All Rights Reserved
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     SUBROUTINE I1KRL (N)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         INTEGER    N
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    I, IN, LALC, LBND, LBOOK, LMAX, LNEED, LNOW, LOUT, &
                             LUSED, NDX, NEXT
                         !                                  SPECIFICATIONS FOR SAVE VARIABLES
                         LOGICAL    FIRST
                         SAVE       FIRST
                         !                                  SPECIFICATIONS FOR SPECIAL CASES
                         !                                  SPECIFICATIONS FOR COMMON /WORKSP/
                         REAL       RWKSP(5000)
                         REAL       RDWKSP(5000)
                         DOUBLE PRECISION DWKSP(2500)
                         COMPLEX    CWKSP(2500)
                         COMPLEX    CZWKSP(2500)
                         COMPLEX    *16 ZWKSP(1250)
                         INTEGER    IWKSP(5000)
                         LOGICAL    LWKSP(5000)
                         EQUIVALENCE (DWKSP(1), RWKSP(1))
                         EQUIVALENCE (CWKSP(1), RWKSP(1)), (ZWKSP(1), RWKSP(1))
                         EQUIVALENCE (IWKSP(1), RWKSP(1)), (LWKSP(1), RWKSP(1))
                         EQUIVALENCE (RDWKSP(1), RWKSP(1)), (CZWKSP(1), RWKSP(1))
                         COMMON     /WORKSP/ RWKSP
                         !                                  SPECIFICATIONS FOR EQUIVALENCE
                         EQUIVALENCE (LOUT, IWKSP(1))
                         EQUIVALENCE (LNOW, IWKSP(2))
                         EQUIVALENCE (LUSED, IWKSP(3))
                         EQUIVALENCE (LBND, IWKSP(4))
                         EQUIVALENCE (LMAX, IWKSP(5))
                         EQUIVALENCE (LALC, IWKSP(6))
                         EQUIVALENCE (LNEED, IWKSP(7))
                         EQUIVALENCE (LBOOK, IWKSP(8))
                         !                                  SPECIFICATIONS FOR SUBROUTINES
                         EXTERNAL   E1MES, E1STI, IWKIN
                         !
                         DATA FIRST/.TRUE./
                         !
                         IF (FIRST) THEN
                             !                                  INITIALIZE WORKSPACE IF NEEDED
                             FIRST = .FALSE.
                             CALL IWKIN (0)
                         END IF
                         !                                  CALLING I1KRL(0) WILL CONFIRM
                         !                                  INTEGRITY OF SYSTEM AND RETURN
                         IF (N .LT. 0) THEN
                             CALL E1MES (5, 10, 'Error from subroutine I1KRL:  Attempt'// &
                                 ' to release a negative number of workspace'// &
                                 ' allocations. ')
                             GO TO 9000
                         END IF
                         !                                  BOOKKEEPING OVERWRITTEN
                         IF (LNOW.LT.LBOOK .OR. LNOW.GT.LUSED .OR. LUSED.GT.LMAX .OR. &
                             LNOW.GE.LBND .OR. LOUT.GT.LALC) THEN
                             CALL E1MES (5, 11, 'Error from subroutine I1KRL:  One or '// &
                                 'more of the first eight bookkeeping locations '// &
                                 'in IWKSP have been overwritten.  ')
                             GO TO 9000
                         END IF
                         !                                  CHECK ALL THE POINTERS IN THE
                         !                                  PERMANENT STORAGE AREA.  THEY MUST
                         !                                  BE MONOTONE INCREASING AND LESS THAN
                         !                                  OR EQUAL TO LMAX, AND THE INDEX OF
                         !                                  THE LAST POINTER MUST BE LMAX+1.
                         NDX = LBND
                         IF (NDX .NE. LMAX+1) THEN
                             DO 10  I=1, LALC
                                 NEXT = IWKSP(NDX)
                                 IF (NEXT .EQ. LMAX+1) GO TO 20
                                 !
                                 IF (NEXT.LE.NDX .OR. NEXT.GT.LMAX) THEN
                                     CALL E1MES (5, 12, 'Error from subroutine I1KRL:  '// &
                                         'A pointer in permanent storage has been '// &
                                         ' overwritten. ')
                                     GO TO 9000
                                 END IF
                                 NDX = NEXT
10                           CONTINUE
                             CALL E1MES (5, 13, 'Error from subroutine I1KRL:  A '// &
                                 'pointer in permanent storage has been '// &
                                 'overwritten. ')
                             GO TO 9000
                         END IF
20                       IF (N .GT. 0) THEN
                             DO 30  IN=1, N
                                 IF (LNOW .LE. LBOOK) THEN
                                     CALL E1MES (5, 14, 'Error from subroutine I1KRL:  '// &
                                         'Attempt to release a nonexistant '// &
                                         'workspace  allocation. ')
                                     GO TO 9000
                                 ELSE IF (IWKSP(LNOW).LT.LBOOK .OR. IWKSP(LNOW).GE.LNOW-1) &
                                     THEN
                                     !                                  CHECK TO MAKE SURE THE BACK POINTERS
                                     !                                  ARE MONOTONE.
                                     CALL E1STI (1, LNOW)
                                     CALL E1MES (5, 15, 'Error from subroutine I1KRL:  '// &
                                         'The pointer at IWKSP(%(I1)) has been '// &
                                         'overwritten.  ')
                                     GO TO 9000
                                 ELSE
                                     LOUT = LOUT - 1
                                     LNOW = IWKSP(LNOW)
                                 END IF
30                           CONTINUE
                         END IF
                         !
9000                     RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  I1KST
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    August 9, 1983
                     !
                     !  Purpose:    Return control information about the workspace stack.
                     !
                     !  Usage:      I1KST(NFACT)
                     !
                     !  Arguments:
                     !     NFACT  - Integer value between 1 and 6 inclusive returns the
                     !                 following information: (Input)
                     !                   NFACT = 1 - LOUT: number of current allocations
                     !                               excluding permanent storage. At the
                     !                               end of a run, there should be no
                     !                               active allocations.
                     !                   NFACT = 2 - LNOW: current active length
                     !                   NFACT = 3 - LTOTAL: total storage used thus far
                     !                   NFACT = 4 - LMAX: maximum storage allowed
                     !                   NFACT = 5 - LALC: total number of allocations made
                     !                               by I1KGT thus far
                     !                   NFACT = 6 - LNEED: number of numeric storage units
                     !                               by which the stack size must be
                     !                               increased for all past allocations
                     !                               to succeed
                     !     I1KST  - Integer function. (Output) Returns a workspace stack
                     !              statistic according to value of NFACT.
                     !
                     !  Copyright:  1983 by IMSL, Inc.  All Rights Reserved
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     INTEGER FUNCTION I1KST (NFACT)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         INTEGER    NFACT
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    ISTATS(7)
                         !                                  SPECIFICATIONS FOR SAVE VARIABLES
                         LOGICAL    FIRST
                         SAVE       FIRST
                         !                                  SPECIFICATIONS FOR SPECIAL CASES
                         !                                  SPECIFICATIONS FOR COMMON /WORKSP/
                         REAL       RWKSP(5000)
                         REAL       RDWKSP(5000)
                         DOUBLE PRECISION DWKSP(2500)
                         COMPLEX    CWKSP(2500)
                         COMPLEX    CZWKSP(2500)
                         COMPLEX    *16 ZWKSP(1250)
                         INTEGER    IWKSP(5000)
                         LOGICAL    LWKSP(5000)
                         EQUIVALENCE (DWKSP(1), RWKSP(1))
                         EQUIVALENCE (CWKSP(1), RWKSP(1)), (ZWKSP(1), RWKSP(1))
                         EQUIVALENCE (IWKSP(1), RWKSP(1)), (LWKSP(1), RWKSP(1))
                         EQUIVALENCE (RDWKSP(1), RWKSP(1)), (CZWKSP(1), RWKSP(1))
                         COMMON     /WORKSP/ RWKSP
                         !                                  SPECIFICATIONS FOR EQUIVALENCE
                         EQUIVALENCE (ISTATS(1), IWKSP(1))
                         !                                  SPECIFICATIONS FOR SUBROUTINES
                         EXTERNAL   E1MES, IWKIN
                         !
                         DATA FIRST/.TRUE./
                         !
                         IF (FIRST) THEN
                             !                                  INITIALIZE WORKSPACE IF NEEDED
                             FIRST = .FALSE.
                             CALL IWKIN (0)
                         END IF
                         !
                         IF (NFACT.LE.0 .OR. NFACT.GE.7) THEN
                             CALL E1MES (5, 9, 'Error from subroutine I1KST:  Argument'// &
                                 ' for I1KST must be between 1 and 6 inclusive.')
                         ELSE IF (NFACT .EQ. 1) THEN
                             !                                  LOUT
                             I1KST = ISTATS(1)
                         ELSE IF (NFACT .EQ. 2) THEN
                             !                                  LNOW + PERMANENT
                             I1KST = ISTATS(2) + (ISTATS(5)-ISTATS(4)+1)
                         ELSE IF (NFACT .EQ. 3) THEN
                             !                                  LUSED + PERMANENT
                             I1KST = ISTATS(3) + (ISTATS(5)-ISTATS(4)+1)
                         ELSE IF (NFACT .EQ. 4) THEN
                             !                                  LMAX
                             I1KST = ISTATS(5)
                         ELSE IF (NFACT .EQ. 5) THEN
                             !                                  LALC
                             I1KST = ISTATS(6)
                         ELSE IF (NFACT .EQ. 6) THEN
                             !                                  LNEED
                             I1KST = ISTATS(7)
                         END IF
                         !
                         RETURN
                     END
                     !-----------------------------------------------------------------------
                     !  IMSL Name:  I1X (Single precision version)
                     !
                     !  Computer:   pcdsms/SINGLE
                     !
                     !  Revised:    August 30, 1985
                     !
                     !  Purpose:    Determine the array subscript indicating the starting
                     !              element at which a key character sequence begins.
                     !              (Case-sensitive version)
                     !
                     !  Usage:      I1X(CHRSTR, I1LEN, KEY, KLEN)
                     !
                     !  Arguments:
                     !     CHRSTR - Character array to be searched.  (Input)
                     !     I1LEN  - Length of CHRSTR.  (Input)
                     !     KEY    - Character array that contains the key sequence.  (Input)
                     !     KLEN   - Length of KEY.  (Input)
                     !     I1X    - Integer function.  (Output)
                     !
                     !  Remarks:
                     !  1. Returns zero when there is no match.
                     !
                     !  2. Returns zero if KLEN is longer than ISLEN.
                     !
                     !  3. Returns zero when any of the character arrays has a negative or
                     !     zero length.
                     !
                     !  GAMS:       N5c
                     !
                     !  Chapter:    MATH/LIBRARY Utilities
                     !
                     !  Copyright:  1985 by IMSL, Inc.  All Rights Reserved.
                     !
                     !  Warranty:   IMSL warrants only that IMSL testing has been applied
                     !              to this code.  No other warranty, expressed or implied,
                     !              is applicable.
                     !
                     !-----------------------------------------------------------------------
                     !
                     INTEGER FUNCTION I1X (CHRSTR, I1LEN, KEY, KLEN)
                         !                                  SPECIFICATIONS FOR ARGUMENTS
                         INTEGER    I1LEN, KLEN
                         CHARACTER  CHRSTR(*), KEY(*)
                         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                         INTEGER    I, II, J
                         !
                         I1X = 0
                         IF (KLEN.LE.0 .OR. I1LEN.LE.0) GO TO 9000
                         IF (KLEN .GT. I1LEN) GO TO 9000
                         !
                         I = 1
                         II = I1LEN - KLEN + 1
10                       IF (I .LE. II) THEN
                             IF (CHRSTR(I) .EQ. KEY(1)) THEN
                                 DO 20  J=2, KLEN
                                     IF (CHRSTR(I+J-1) .NE. KEY(J)) GO TO 30
20                               CONTINUE
                                 I1X = I
                                 GO TO 9000
30                           CONTINUE
                         END IF
                         I = I + 1
                         GO TO 10
                     END IF
                     !
9000                 RETURN
                 END
                 !-----------------------------------------------------------------------
                 !  IMSL Name:  IACHAR (Single precision version)
                 !
                 !  Computer:   pcdsms/SINGLE
                 !
                 !  Revised:    September 9, 1985
                 !
                 !  Purpose:    Return the integer ASCII value of a character argument.
                 !
                 !  Usage:      IACHAR(CH)
                 !
                 !  Arguments:
                 !     CH     - Character argument for which the integer ASCII value
                 !              is desired.  (Input)
                 !     IACHAR - Integer ASCII value for CH.  (Output)
                 !              The character CH is in the IACHAR-th position of the
                 !              ASCII collating sequence.
                 !
                 !  GAMS:       N3
                 !
                 !  Chapter:    MATH/LIBRARY Utilities
                 !              STAT/LIBRARY Utilities
                 !
                 !  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
                 !
                 !  Warranty:   IMSL warrants only that IMSL testing has been applied
                 !              to this code.  No other warranty, expressed or implied,
                 !              is applicable.
                 !
                 !-----------------------------------------------------------------------
                 !
                 INTEGER FUNCTION IACHAR (CH)
                     !                                  SPECIFICATIONS FOR ARGUMENTS
                     CHARACTER  CH
                     !                                  SPECIFICATIONS FOR SAVE VARIABLES
                     IACHAR = ICHAR(CH)
                     RETURN
                 END
                 !-----------------------------------------------------------------------
                 !  IMSL Name:  ICASE (Single precision version)
                 !
                 !  Computer:   pcdsms/SINGLE
                 !
                 !  Revised:    September 9, 1985
                 !
                 !  Purpose:    Convert from character to the integer ASCII value without
                 !              regard to case.
                 !
                 !  Usage:      ICASE(CH)
                 !
                 !  Arguments:
                 !     CH     - Character to be converted.  (Input)
                 !     ICASE  - Integer ASCII value for CH without regard to the case
                 !              of CH.  (Output)
                 !              ICASE returns the same value as IMSL routine IACHAR for
                 !              all but lowercase letters.  For these, it returns the
                 !              IACHAR value for the corresponding uppercase letter.
                 !
                 !  GAMS:       N3
                 !
                 !  Chapter:    MATH/LIBRARY Utilities
                 !              STAT/LIBRARY Utilities
                 !
                 !  Copyright:  1986 by IMSL, Inc.  All Rights Reserved.
                 !
                 !  Warranty:   IMSL warrants only that IMSL testing has been applied
                 !              to this code.  No other warranty, expressed or implied,
                 !              is applicable.
                 !
                 !-----------------------------------------------------------------------
                 !
                 INTEGER FUNCTION ICASE (CH)
                     !                                  SPECIFICATIONS FOR ARGUMENTS
                     CHARACTER  CH
                     !                                  SPECIFICATIONS FOR FUNCTIONS
                     EXTERNAL   IACHAR
                     INTEGER    IACHAR
                     !
                     ICASE = IACHAR(CH)
                     IF (ICASE.GE.97 .AND. ICASE.LE.122) ICASE = ICASE - 32
                     !
                     RETURN
                 END
                 !-----------------------------------------------------------------------
                 !  IMSL Name:  IMACH (Single precision version)
                 !
                 !  Computer:   pcdsms/SINGLE
                 !
                 !  Revised:    March 26, 1984
                 !
                 !  Purpose:    Retrieve integer machine constants.
                 !
                 !  Usage:      IMACH(N)
                 !
                 !  Arguments:
                 !     N      - Index of desired constant.  (Input)
                 !     IMACH  - Machine constant.  (Output)
                 !
                 !  Remark:
                 !     Following is a description of the assorted integer machine
                 !     constants.
                 !
                 !     Words
                 !
                 !        IMACH( 1) = Number of bits per integer storage unit.
                 !        IMACH( 2) = Number of characters per integer storage unit.
                 !
                 !     Integers
                 !
                 !        Assume integers are represented in the S-DIGIT, BASE-A form
                 !        SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
                 !        where 0 .LE. X(I) .LT. A for I=0,...,S-1.  Then
                 !
                 !        IMACH( 3) = A, the base.
                 !        IMACH( 4) = S, number of BASE-A digits.
                 !        IMACH( 5) = A**S - 1, largest magnitude.
                 !
                 !     Floating-point numbers
                 !
                 !        Assume floating-point numbers are represented in the T-DIGIT,
                 !        BASE-B form SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
                 !        where 0 .LE. X(I) .LT. B for I=1,...,T,
                 !        0 .LT. X(1), and EMIN .LE. E .LE. EMAX.  Then
                 !
                 !        IMACH( 6) = B, the base.
                 !
                 !        Single precision
                 !
                 !           IMACH( 7) = T, number of BASE-B digits.
                 !           IMACH( 8) = EMIN, smallest exponent E.
                 !           IMACH( 9) = EMAX, largest exponent E.
                 !
                 !        Double precision
                 !
                 !           IMACH(10) = T, number of BASE-B digits.
                 !           IMACH(11) = EMIN, smallest exponent E.
                 !           IMACH(12) = EMAX, largest exponent E.
                 !
                 !  GAMS:       R1
                 !
                 !  Chapters:   MATH/LIBRARY Reference Material
                 !              STAT/LIBRARY Reference Material
                 !              SFUN/LIBRARY Reference Material
                 !
                 !  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
                 !
                 !  Warranty:   IMSL warrants only that IMSL testing has been applied
                 !              to this code.  No other warranty, expressed or implied,
                 !              is applicable.
                 !
                 !-----------------------------------------------------------------------
                 !
                 INTEGER FUNCTION IMACH (N)
                     !                                  SPECIFICATIONS FOR ARGUMENTS
                     INTEGER    N
                     !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                     INTEGER    NOUT
                     !                                  SPECIFICATIONS FOR SAVE VARIABLES
                     INTEGER    IMACHV(12)
                     SAVE       IMACHV
                     !                                  SPECIFICATIONS FOR SUBROUTINES
                     EXTERNAL   UMACH
                     !                                  DEFINE CONSTANTS
                     DATA IMACHV(1)/32/
                     DATA IMACHV(2)/4/
                     DATA IMACHV(3)/2/
                     DATA IMACHV(4)/31/
                     DATA IMACHV(5)/2147483647/
                     DATA IMACHV(6)/2/
                     DATA IMACHV(7)/24/
                     DATA IMACHV(8)/-125/
                     DATA IMACHV(9)/128/
                     DATA IMACHV(10)/53/
                     DATA IMACHV(11)/-1021/
                     DATA IMACHV(12)/1024/
                     !
                     IF (N.LT.1 .OR. N.GT.12) THEN
                         !                                  ERROR.  INVALID RANGE FOR N.
                         CALL UMACH (2, NOUT)
                         WRITE (NOUT,99999) N
99999                    FORMAT (/, ' *** TERMINAL ERROR 5 from IMACH.  The argument', &
                             /, ' ***          must be between 1 and 12 inclusive.' &
                             , /, ' ***          N = ', I6, '.', /)
                         IMACH = 0
                         STOP
                     !
                     ELSE
                         IMACH = IMACHV(N)
                     END IF
                     !
                     RETURN
                 END
                 !-----------------------------------------------------------------------
                 !  IMSL Name:  IWKIN (Single precision version)
                 !
                 !  Computer:   pcdsms/SINGLE
                 !
                 !  Revised:    January 17, 1984
                 !
                 !  Purpose:    Initialize bookkeeping locations describing the
                 !              workspace stack.
                 !
                 !  Usage:      CALL IWKIN (NSU)
                 !
                 !  Argument:
                 !     NSU    - Number of numeric storage units to which the workspace
                 !              stack is to be initialized
                 !
                 !  GAMS:       N4
                 !
                 !  Chapters:   MATH/LIBRARY Reference Material
                 !              STAT/LIBRARY Reference Material
                 !
                 !  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
                 !
                 !  Warranty:   IMSL warrants only that IMSL testing has been applied
                 !              to this code.  No other warranty, expressed or implied,
                 !              is applicable.
                 !
                 !-----------------------------------------------------------------------
                 !
                 SUBROUTINE IWKIN (NSU)
                     !                                  SPECIFICATIONS FOR ARGUMENTS
                     INTEGER    NSU
                     !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                     INTEGER    ISIZE(6), LALC, LBND, LBOOK, LMAX, LNEED, LNOW, LOUT, &
                         LUSED, MELMTS, MTYPE
                     !                                  SPECIFICATIONS FOR SAVE VARIABLES
                     LOGICAL    FIRST
                     SAVE       FIRST
                     !                                  SPECIFICATIONS FOR SPECIAL CASES
                     !                                  SPECIFICATIONS FOR COMMON /WORKSP/
                     REAL       RWKSP(5000)
                     REAL       RDWKSP(5000)
                     DOUBLE PRECISION DWKSP(2500)
                     COMPLEX    CWKSP(2500)
                     COMPLEX    CZWKSP(2500)
                     COMPLEX    *16 ZWKSP(1250)
                     INTEGER    IWKSP(5000)
                     LOGICAL    LWKSP(5000)
                     EQUIVALENCE (DWKSP(1), RWKSP(1))
                     EQUIVALENCE (CWKSP(1), RWKSP(1)), (ZWKSP(1), RWKSP(1))
                     EQUIVALENCE (IWKSP(1), RWKSP(1)), (LWKSP(1), RWKSP(1))
                     EQUIVALENCE (RDWKSP(1), RWKSP(1)), (CZWKSP(1), RWKSP(1))
                     COMMON     /WORKSP/ RWKSP
                     !                                  SPECIFICATIONS FOR EQUIVALENCE
                     EQUIVALENCE (LOUT, IWKSP(1))
                     EQUIVALENCE (LNOW, IWKSP(2))
                     EQUIVALENCE (LUSED, IWKSP(3))
                     EQUIVALENCE (LBND, IWKSP(4))
                     EQUIVALENCE (LMAX, IWKSP(5))
                     EQUIVALENCE (LALC, IWKSP(6))
                     EQUIVALENCE (LNEED, IWKSP(7))
                     EQUIVALENCE (LBOOK, IWKSP(8))
                     EQUIVALENCE (ISIZE(1), IWKSP(11))
                     !                                  SPECIFICATIONS FOR INTRINSICS
                     !     INTRINSIC  MAX0
                     INTRINSIC  MAX0
                     INTEGER    MAX0
                     !                                  SPECIFICATIONS FOR SUBROUTINES
                     EXTERNAL   E1MES, E1STI
                     !
                     DATA FIRST/.TRUE./
                     !
                     IF (.NOT.FIRST) THEN
                         IF (NSU .NE. 0) THEN
                             CALL E1STI (1, LMAX)
                             CALL E1MES (5, 100, 'Error from subroutine IWKIN:  '// &
                                 'Workspace stack has previously been '// &
                                 'initialized to %(I1). Correct by making the '// &
                                 'call to IWKIN the first executable '// &
                                 'statement in the main program.  ')
                             !
                             STOP
                         !
                         ELSE
                             RETURN
                         END IF
                     END IF
                     !
                     IF (NSU .EQ. 0) THEN
                         !                                  IF NSU=0 USE DEFAULT SIZE 5000
                         MELMTS = 5000
                     ELSE
                         MELMTS = NSU
                     END IF
                     !                                  NUMBER OF ITEMS .LT. 0
                     IF (MELMTS .LE. 0) THEN
                         CALL E1STI (1, MELMTS)
                         CALL E1MES (5, 1, 'Error from subroutine IWKIN:  Number '// &
                             'of numeric storage units is not positive. NSU '// &
                             '= %(I1) ')
                     ELSE
                         !
                         FIRST = .FALSE.
                         !                                  HERE TO INITIALIZE
                         !
                         !                                  SET DATA SIZES APPROPRIATE FOR A
                         !                                  STANDARD CONFORMING FORTRAN SYSTEM
                         !                                  USING THE FORTRAN
                         !                                  *NUMERIC STORAGE UNIT* AS THE
                         !                                  MEASURE OF SIZE.
                         !
                         !                                  TYPE IS REAL
                         MTYPE = 3
                         !                                  LOGICAL
                         ISIZE(1) = 1
                         !                                  INTEGER
                         ISIZE(2) = 1
                         !                                  REAL
                         ISIZE(3) = 1
                         !                                  DOUBLE PRECISION
                         ISIZE(4) = 2
                         !                                  COMPLEX
                         ISIZE(5) = 2
                         !                                  DOUBLE COMPLEX
                         ISIZE(6) = 4
                         !                                  NUMBER OF WORDS USED FOR BOOKKEEPING
                         LBOOK = 16
                         !                                  CURRENT ACTIVE LENGTH OF THE STACK
                         LNOW = LBOOK
                         !                                  MAXIMUM VALUE OF LNOW ACHIEVED THUS
                         !                                  FAR
                         LUSED = LBOOK
                         !                                  MAXIMUM LENGTH OF THE STORAGE ARRAY
                         LMAX = MAX0(MELMTS,((LBOOK+2)*ISIZE(2)+ISIZE(3)-1)/ISIZE(3))
                         !                                  LOWER BOUND OF THE PERMANENT STORAGE
                         !                                  WHICH IS ONE WORD MORE THAN THE
                         !                                  MAXIMUM ALLOWED LENGTH OF THE STACK
                         LBND = LMAX + 1
                         !                                  NUMBER OF CURRENT ALLOCATIONS
                         LOUT = 0
                         !                                  TOTAL NUMBER OF ALLOCATIONS MADE
                         LALC = 0
                         !                                  NUMBER OF WORDS BY WHICH THE ARRAY
                         !                                  SIZE MUST BE INCREASED FOR ALL PAST
                         !                                  ALLOCATIONS TO SUCCEED
                         LNEED = 0
                     END IF
                     !
                     RETURN
                 END
                 !-----------------------------------------------------------------------
                 !  IMSL Name:  M1VE
                 !
                 !  Computer:   pcdsms/SINGLE
                 !
                 !  Revised:    March 5, 1984
                 !
                 !  Purpose:    Move a subset of one character array to another.
                 !
                 !  Usage:      CALL M1VE(INSTR, INBEG, INEND, INLEN, OUTSTR, OUTBEG,
                 !                         OUTEND, OUTLEN, IER)
                 !
                 !  Arguments:
                 !     INSTR  - Source character array.  (Input)
                 !     INBEG  - First element of INSTR to be moved.  (Input)
                 !     INEND  - Last element of INSTR to be moved.  (Input)
                 !              The source subset is INSTR(INBEG),...,INSTR(INEND).
                 !     INLEN  - Length of INSTR.  (Input)
                 !     OUTSTR - Destination character array.  (Output)
                 !     IUTBEG - First element of OUTSTR destination.  (Input)
                 !     IUTEND - Last element of OUTSTR  destination.  (Input)
                 !              The destination subset is OUTSRT(IUTBEG),...,
                 !              OUTSTR(IUTEND).
                 !     IUTLEN - Length of OUTSTR.  (Input)
                 !     IER    - Completion code.  (Output)
                 !              IER = -2  indicates that the input parameters, INBEG,
                 !                        INEND, INLEN, IUTBEG, IUTEND are not
                 !                        consistent.  One of the conditions
                 !                        INBEG.GT.0, INEND.GE.INBEG, INLEN.GE.INEND,
                 !                        IUTBEG.GT.0, or IUTEND.GE.IUTBEG is not
                 !                        satisfied.
                 !              IER = -1  indicates that the length of OUTSTR is
                 !                        insufficient to hold the subset of INSTR.
                 !                        That is, IUTLEN is less than IUTEND.
                 !              IER =  0  indicates normal completion
                 !              IER >  0  indicates that the specified subset of OUTSTR,
                 !                        OUTSTR(IUTBEG),...,OUTSTR(IUTEND) is not long
                 !                        enough to hold the subset INSTR(INBEG),...,
                 !                        INSTR(INEND) of INSTR.  IER is set to the
                 !                        number of characters that were not moved.
                 !
                 !  Remarks:
                 !  1. If the subset of OUTSTR is longer than the subset of INSTR,
                 !     trailing blanks are moved to OUTSTR.
                 !  2. If the subset of INSTR is longer than the subset of OUTSTR,
                 !     the shorter subset is moved to OUTSTR and IER is set to the number
                 !     of characters that were not moved to OUTSTR.
                 !  3. If the length of OUTSTR is insufficient to hold the subset,
                 !     IER is set to -2 and nothing is moved.
                 !
                 !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
                 !
                 !  Warranty:   IMSL warrants only that IMSL testing has been applied
                 !              to this code.  No other warranty, expressed or implied,
                 !              is applicable.
                 !
                 !-----------------------------------------------------------------------
                 !
                 SUBROUTINE M1VE (INSTR, INBEG, INEND, INLEN, OUTSTR, IUTBEG, &
                     IUTEND, IUTLEN, IER)
                     !                                  SPECIFICATIONS FOR ARGUMENTS
                     INTEGER    INBEG, INEND, INLEN, IUTBEG, IUTEND, IUTLEN, IER
                     CHARACTER  INSTR(*), OUTSTR(*)
                     !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                     INTEGER    IUTLAS, KI, KO
                     !                                  SPECIFICATIONS FOR SAVE VARIABLES
                     CHARACTER  BLANK
                     SAVE       BLANK
                     !                                  SPECIFICATIONS FOR INTRINSICS
                     !     INTRINSIC  MIN0
                     INTRINSIC  MIN0
                     INTEGER    MIN0
                     !
                     DATA BLANK/' '/
                     !                                  CHECK INBEG, INEND, INLEN, IUTBEG,
                     !                                  AND IUTEND
                     !
                     IF (INBEG.LE.0 .OR. INEND.LT.INBEG .OR. INLEN.LT.INEND .OR. &
                         IUTBEG.LE.0 .OR. IUTEND.LT.IUTBEG) THEN
                         IER = -2
                         RETURN
                     ELSE IF (IUTLEN .LT. IUTEND) THEN
                         IER = -1
                         RETURN
                     END IF
                     !                                  DETERMINE LAST CHARACTER TO M1VE
                     IUTLAS = IUTBEG + MIN0(INEND-INBEG,IUTEND-IUTBEG)
                     !                                  M1VE CHARACTERS
                     KI = INBEG
                     DO 10  KO=IUTBEG, IUTLAS
                         OUTSTR(KO) = INSTR(KI)
                         KI = KI + 1
10                   CONTINUE
                     !                                   SET IER TO NUMBER OF CHARACTERS THAT
                     !                                   WHERE NOT MOVED
                     IER = KI - INEND - 1
                     !                                   APPEND BLANKS IF NECESSARY
                     DO 20  KO=IUTLAS + 1, IUTEND
                         OUTSTR(KO) = BLANK
20                   CONTINUE
                     !
                     RETURN
                 END
                 !-----------------------------------------------------------------------
                 !  IMSL Name:  M1VECH
                 !
                 !  Computer:   pcdsms/SINGLE
                 !
                 !  Revised:    December 31, 1984
                 !
                 !  Purpose:    Character substring assignment.
                 !
                 !  Usage:      CALL M1VECH (STR1, LEN1, STR2, LEN2)
                 !
                 !  Arguments:
                 !     STR1   - Source substring.  (Input)
                 !              The source substring is STR1(1:LEN1).
                 !     LEN1   - Length of STR1.  (Input)
                 !     STR2   - Destination substring.  (Output)
                 !              The destination substring is STR2(1:LEN2).
                 !     LEN2   - Length of STR2.  (Input)
                 !
                 !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
                 !
                 !  Warranty:   IMSL warrants only that IMSL testing has been applied
                 !              to this code.  No other warranty, expressed or implied,
                 !              is applicable.
                 !
                 !-----------------------------------------------------------------------
                 !
                 SUBROUTINE M1VECH (STR1, LEN1, STR2, LEN2)
                     !                                  SPECIFICATIONS FOR ARGUMENTS
                     INTEGER    LEN1, LEN2
                     CHARACTER  STR1*(*), STR2*(*)
                     !
                     STR2(1:LEN2) = STR1(1:LEN1)
                     !
                     RETURN
                 END
                 !-----------------------------------------------------------------------
                 !  IMSL Name:  N1RGB
                 !
                 !  Computer:   pcdsms/SINGLE
                 !
                 !  Revised:    March 2, 1984
                 !
                 !  Purpose:    Return a positive number as a flag to indicated that a
                 !              stop should occur due to one or more global errors.
                 !
                 !  Usage:      N1RGB(IDUMMY)
                 !
                 !  Arguments:
                 !     IDUMMY - Integer scalar dummy argument.
                 !
                 !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
                 !
                 !  Warranty:   IMSL warrants only that IMSL testing has been applied
                 !              to this code.  No other warranty, expressed or implied,
                 !              is applicable.
                 !
                 !-----------------------------------------------------------------------
                 !
                 INTEGER FUNCTION N1RGB (IDUMMY)
                     !                                  SPECIFICATIONS FOR ARGUMENTS
                     INTEGER    IDUMMY
                     !                                  SPECIFICATIONS FOR SPECIAL CASES
                     !                              SPECIFICATIONS FOR COMMON /ERCOM1/
                     INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                         PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                         IALLOC(51), HDRFMT(7), TRACON(7)
                     COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                         PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                         TRACON
                     SAVE       /ERCOM1/
                     !                              SPECIFICATIONS FOR COMMON /ERCOM2/
                     CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
                     COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
                     SAVE       /ERCOM2/
                     !                              SPECIFICATIONS FOR COMMON /ERCOM3/
                     DOUBLE PRECISION ERCKSM
                     COMMON     /ERCOM3/ ERCKSM
                     SAVE       /ERCOM3/
                     !                              SPECIFICATIONS FOR COMMON /ERCOM4/
                     LOGICAL    ISUSER(51)
                     COMMON     /ERCOM4/ ISUSER
                     SAVE       /ERCOM4/
                     !                                  INITIALIZE FUNCTION
                     N1RGB = 0
                     !                                  CHECK FOR GLOBAL ERROR TYPE 6
                     IF (IFERR6 .GT. 0) THEN
                         N1RGB = STOPTB(6)
                         IFERR6 = 0
                     END IF
                     !                                  CHECK FOR GLOBAL ERROR TYPE 7
                     IF (IFERR7 .GT. 0) THEN
                         N1RGB = STOPTB(7)
                         IFERR7 = 0
                     END IF
                     !
                     RETURN
                 END
                 !-----------------------------------------------------------------------
                 !  IMSL Name:  N1RTY
                 !
                 !  Computer:   pcdsms/SINGLE
                 !
                 !  Revised:    March 6, 1984
                 !
                 !  Purpose:    Retrieve an error type.
                 !
                 !  Usage:      N1RTY(IOPT)
                 !
                 !  Arguments:
                 !     IOPT   - Integer specifying the level.  (Input)
                 !              If IOPT=0 the error type for the current level is
                 !              returned.  If IOPT=1 the error type for the most
                 !              recently called routine (last pop) is returned.
                 !
                 !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
                 !
                 !  Warranty:   IMSL warrants only that IMSL testing has been applied
                 !              to this code.  No other warranty, expressed or implied,
                 !              is applicable.
                 !
                 !-----------------------------------------------------------------------
                 !
                 INTEGER FUNCTION N1RTY (IOPT)
                     !                                  SPECIFICATIONS FOR ARGUMENTS
                     INTEGER    IOPT
                     !                                  SPECIFICATIONS FOR SPECIAL CASES
                     !                              SPECIFICATIONS FOR COMMON /ERCOM1/
                     INTEGER    CALLVL, MAXLEV, MSGLEN, ERTYPE(51), ERCODE(51), &
                         PRINTB(7), STOPTB(7), PLEN, IFERR6, IFERR7, &
                         IALLOC(51), HDRFMT(7), TRACON(7)
                     COMMON     /ERCOM1/ CALLVL, MAXLEV, MSGLEN, ERTYPE, ERCODE, &
                         PRINTB, STOPTB, PLEN, IFERR6, IFERR7, IALLOC, HDRFMT, &
                         TRACON
                     SAVE       /ERCOM1/
                     !                              SPECIFICATIONS FOR COMMON /ERCOM2/
                     CHARACTER  MSGSAV(255), PLIST(300), RNAME(51)*6
                     COMMON     /ERCOM2/ MSGSAV, PLIST, RNAME
                     SAVE       /ERCOM2/
                     !                              SPECIFICATIONS FOR COMMON /ERCOM3/
                     DOUBLE PRECISION ERCKSM
                     COMMON     /ERCOM3/ ERCKSM
                     SAVE       /ERCOM3/
                     !                              SPECIFICATIONS FOR COMMON /ERCOM4/
                     LOGICAL    ISUSER(51)
                     COMMON     /ERCOM4/ ISUSER
                     SAVE       /ERCOM4/
                     !                                  SPECIFICATIONS FOR SUBROUTINES
                     EXTERNAL   E1PRT, M1VECH
                     !
                     IF (IOPT.NE.0 .AND. IOPT.NE.1) THEN
                         ERTYPE(CALLVL) = 5
                         ERCODE(CALLVL) = 1
                         MSGLEN = 47
                         CALL M1VECH ('.  The argument passed to N1RTY must be 0 or '// &
                             '1. ', MSGLEN, MSGSAV, MSGLEN)
                         CALL E1PRT
                         STOP
                     ELSE
                         N1RTY = ERTYPE(CALLVL+IOPT)
                     END IF
                     !
                     RETURN
                 END
                 !-----------------------------------------------------------------------
                 !  IMSL Name:  S1ANUM
                 !
                 !  Computer:   pcdsms/SINGLE
                 !
                 !  Revised:    March 28, 1984
                 !
                 !  Purpose:    Scan a token and identify it as follows: integer, real
                 !              number (single/double), FORTRAN relational operator,
                 !              FORTRAN logical operator, or FORTRAN logical constant.
                 !
                 !  Usage:      CALL S1ANUM(INSTR, SLEN, CODE, OLEN)
                 !
                 !  Arguments:
                 !     INSTR  - Character string to be scanned.  (Input)
                 !     SLEN   - Length of INSTR.  (Input)
                 !     CODE   - Token code.  (Output)  Where
                 !                 CODE =  0  indicates an unknown token,
                 !                 CODE =  1  indicates an integer number,
                 !                 CODE =  2  indicates a (single precision) real number,
                 !                 CODE =  3  indicates a (double precision) real number,
                 !                 CODE =  4  indicates a logical constant (.TRUE. or
                 !                               .FALSE.),
                 !                 CODE =  5  indicates the relational operator .EQ.,
                 !                 CODE =  6  indicates the relational operator .NE.,
                 !                 CODE =  7  indicates the relational operator .LT.,
                 !                 CODE =  8  indicates the relational operator .LE.,
                 !                 CODE =  9  indicates the relational operator .GT.,
                 !                 CODE = 10  indicates the relational operator .GE.,
                 !                 CODE = 11  indicates the logical operator .AND.,
                 !                 CODE = 12  indicates the logical operator .OR.,
                 !                 CODE = 13  indicates the logical operator .EQV.,
                 !                 CODE = 14  indicates the logical operator .NEQV.,
                 !                 CODE = 15  indicates the logical operator .NOT..
                 !     OLEN   - Length of the token as counted from the first character
                 !              in INSTR.  (Output)  OLEN returns a zero for an unknown
                 !              token (CODE = 0).
                 !
                 !  Remarks:
                 !  1. Blanks are considered significant.
                 !  2. Lower and upper case letters are not significant.
                 !
                 !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
                 !
                 !  Warranty:   IMSL warrants only that IMSL testing has been applied
                 !              to this code.  No other warranty, expressed or implied,
                 !              is applicable.
                 !
                 !-----------------------------------------------------------------------
                 !
                 SUBROUTINE S1ANUM (INSTR, SLEN, CODE, OLEN)
                     !                                  SPECIFICATIONS FOR ARGUMENTS
                     INTEGER    SLEN, CODE, OLEN
                     CHARACTER  INSTR(*)
                     !                                  SPECIFICATIONS FOR LOCAL VARIABLES
                     INTEGER    I, IBEG, IIBEG, J
                     LOGICAL    FLAG
                     CHARACTER  CHRSTR(6)
                     !                                  SPECIFICATIONS FOR SAVE VARIABLES
                     INTEGER    TABPTR(16), TDCNST, TICNST, TOKEN(13), TRCNST, TZERR
                     CHARACTER  DIGIT(10), LETTER(52), MINUS, PERIOD, PLUS, TABLE(38)
                     SAVE       DIGIT, LETTER, MINUS, PERIOD, PLUS, TABLE, TABPTR, &
                         TDCNST, TICNST, TOKEN, TRCNST, TZERR
                     !                                  SPECIFICATIONS FOR FUNCTIONS
                     EXTERNAL   I1X, I1CSTR
                     INTEGER    I1X, I1CSTR
                     !
                     DATA TOKEN/5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 4, 4/
                     DATA TABLE/'D', 'E', 'E', 'Q', 'N', 'E', 'L', 'T', 'L', &
                         'E', 'G', 'T', 'G', 'E', 'A', 'N', 'D', 'O', 'R', &
                         'E', 'Q', 'V', 'N', 'E', 'Q', 'V', 'N', 'O', 'T', &
                         'T', 'R', 'U', 'E', 'F', 'A', 'L', 'S', 'E'/
                     DATA TABPTR/1, 2, 3, 5, 7, 9, 11, 13, 15, 18, 20, 23, 27, 30, &
                         34, 39/
                     DATA DIGIT/'0', '1', '2', '3', '4', '5', '6', '7', '8', &
                         '9'/
                     DATA LETTER/'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', &
                         'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', &
                         'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'a', 'b', 'c', &
                         'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', &
                         'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', &
                         'x', 'y', 'z'/
                     DATA PERIOD/'.'/, PLUS/'+'/, MINUS/'-'/
                     DATA TZERR/0/, TICNST/1/
                     DATA TRCNST/2/, TDCNST/3/
                     !
                     IF (SLEN .LE. 0) THEN
                         CODE = 0
                         OLEN = 0
                         RETURN
                     END IF
                     !                                  STATE 0 - ASSUME ERROR TOKEN
                     IBEG = 1
                     CODE = TZERR
                     !                                  CHECK SIGN
                     IF (INSTR(IBEG).EQ.MINUS .OR. INSTR(IBEG).EQ.PLUS) THEN
                         FLAG = .TRUE.
                         IIBEG = IBEG
                         IBEG = IBEG + 1
                     ELSE
                         FLAG = .FALSE.
                     END IF
                     !                                  STATE 1 - ASSUME INTEGER CONSTANT
                     IF (I1X(DIGIT,10,INSTR(IBEG),1) .NE. 0) THEN
                         CODE = TICNST
                         IIBEG = IBEG
                         IBEG = IBEG + 1
                         !
10                       IF (IBEG .LE. SLEN) THEN
                             !
                             IF (I1X(DIGIT,10,INSTR(IBEG),1) .NE. 0) THEN
                                 IIBEG = IBEG
                                 IBEG = IBEG + 1
                                 GO TO 10
                             !
                             END IF
                         !
                         ELSE
                             GO TO 80
                         !
                         END IF
                         !
                         IF (INSTR(IBEG) .NE. PERIOD) GO TO 80
                     END IF
                     !                                  STATE 2 - ASSUME REAL CONSTANT
                     IF (CODE .EQ. TICNST) THEN
                         CODE = TRCNST
                         IIBEG = IBEG
                         IBEG = IBEG + 1
                         IF (IBEG .GT. SLEN) GO TO 80
                     ELSE IF (INSTR(IBEG).EQ.PERIOD .AND. SLEN.GE.2) THEN
                         IF (I1X(DIGIT,10,INSTR(IBEG+1),1) .NE. 0) THEN
                             CODE = TRCNST
                             IIBEG = IBEG + 1
                             IBEG = IBEG + 2
                             IF (IBEG .GT. SLEN) GO TO 80
                         END IF
                     END IF
                     !
                     IF (I1X(DIGIT,10,INSTR(IBEG),1) .NE. 0) THEN
                         CODE = TRCNST
                         IIBEG = IBEG
                         IBEG = IBEG + 1
                         !
20                       IF (IBEG .LE. SLEN) THEN
                             !
                             IF (I1X(DIGIT,10,INSTR(IBEG),1) .NE. 0) THEN
                                 IIBEG = IBEG
                                 IBEG = IBEG + 1
                                 GO TO 20
                             !
                             END IF
                         !
                         ELSE
                             GO TO 80
                         !
                         END IF
                     !
                     END IF
                     !
                     IF (CODE .EQ. TZERR) THEN
                         IF (INSTR(IBEG) .NE. PERIOD) GO TO 80
                         IBEG = IBEG + 1
                         IF (IBEG .GT. SLEN) GO TO 80
                     END IF
                     !
                     IF (I1X(LETTER,52,INSTR(IBEG),1) .EQ. 0) GO TO 80
                     CHRSTR(1) = INSTR(IBEG)
                     !
                     DO 30  I=2, 6
                         IBEG = IBEG + 1
                         IF (IBEG .GT. SLEN) GO TO 80
                         IF (I1X(LETTER,52,INSTR(IBEG),1) .EQ. 0) GO TO 40
                         CHRSTR(I) = INSTR(IBEG)
30                   CONTINUE
                     !
                     GO TO 80
                 !
40               CONTINUE
                 !
                 DO 50  J=1, 15
                     IF (I1CSTR(CHRSTR,I-1,TABLE(TABPTR(J)),TABPTR(J+1)-TABPTR(J)) &
                         .EQ. 0) GO TO 60
50               CONTINUE
                 !
                 GO TO 80
                 !                                  STATE 4 - LOGICAL OPERATOR
60               IF (J .GT. 2) THEN
                     !
                     IF (CODE .EQ. TRCNST) THEN
                         !
                         IF (INSTR(IBEG) .EQ. PERIOD) THEN
                             CODE = TICNST
                             IIBEG = IIBEG - 1
                         END IF
                         !
                         GO TO 80
                     !
                     ELSE IF (INSTR(IBEG) .NE. PERIOD) THEN
                         GO TO 80
                     !
                     ELSE IF (FLAG) THEN
                         GO TO 80
                     !
                     ELSE
                         CODE = TOKEN(J-2)
                         IIBEG = IBEG
                         GO TO 80
                     !
                     END IF
                 !
                 END IF
                 !                                  STATE 5 - DOUBLE PRECISION CONSTANT
                 IF (CODE .NE. TRCNST) GO TO 80
                 IF (INSTR(IBEG).EQ.MINUS .OR. INSTR(IBEG).EQ.PLUS) IBEG = IBEG + &
                     1
                 IF (IBEG .GT. SLEN) GO TO 80
                 !
                 IF (I1X(DIGIT,10,INSTR(IBEG),1) .EQ. 0) THEN
                     GO TO 80
                 !
                 ELSE
                     IIBEG = IBEG
                     IBEG = IBEG + 1
                     !
70                   IF (IBEG .LE. SLEN) THEN
                         !
                         IF (I1X(DIGIT,10,INSTR(IBEG),1) .NE. 0) THEN
                             IIBEG = IBEG
                             IBEG = IBEG + 1
                             GO TO 70
                         !
                         END IF
                     !
                     END IF
                 !
                 END IF
                 !
                 IF (J .EQ. 1) CODE = TDCNST
             !
80           CONTINUE
             !
             IF (CODE .EQ. TZERR) THEN
                 OLEN = 0
             !
             ELSE
                 OLEN = IIBEG
             END IF
             !
             RETURN
         END
         !-----------------------------------------------------------------------
         !  IMSL Name:  UMACH (Single precision version)
         !
         !  Computer:   pcdsms/SINGLE
         !
         !  Revised:    March 21, 1984
         !
         !  Purpose:    Set or retrieve input or output device unit numbers.
         !
         !  Usage:      CALL UMACH (N, NUNIT)
         !
         !  Arguments:
         !     N      - Index of desired unit.  (Input)
         !              The values of N are defined as follows:
         !              N = 1, corresponds to the standard input unit.
         !              N = 2, corresponds to the standard output unit.
         !     NUNIT  - I/O unit.  (Input or Output)
         !              If the value of N is negative, the unit corresponding
         !              to the index is reset to the value given in NUNIT.
         !              Otherwise, the value corresponding to the index is
         !              returned in NUNIT.
         !
         !  GAMS:       R1
         !
         !  Chapters:   MATH/LIBRARY Reference Material
         !              STAT/LIBRARY Reference Material
         !              SFUN/LIBRARY Reference Material
         !
         !  Copyright:  1984 by IMSL, Inc.  All Rights Reserved.
         !
         !  Warranty:   IMSL warrants only that IMSL testing has been applied
         !              to this code.  No other warranty, expressed or implied,
         !              is applicable.
         !
         !-----------------------------------------------------------------------
         !
         SUBROUTINE UMACH (N, NUNIT)
             !                                  SPECIFICATIONS FOR ARGUMENTS
             INTEGER    N, NUNIT
             !                                  SPECIFICATIONS FOR LOCAL VARIABLES
             INTEGER    NN, NOUT
             !                                  SPECIFICATIONS FOR SAVE VARIABLES
             INTEGER    UNIT(2)
             SAVE       UNIT
             !                                  SPECIFICATIONS FOR INTRINSICS
             !     INTRINSIC  IABS
             INTRINSIC  IABS
             INTEGER    IABS
             !
             DATA UNIT(1)/5/
             DATA UNIT(2)/6/
             !
             NN = IABS(N)
             IF (NN.NE.1 .AND. NN.NE.2) THEN
                 !                                  ERROR.  INVALID RANGE FOR N.
                 NOUT = UNIT(2)
                 WRITE (NOUT,99999) NN
99999            FORMAT (/, ' *** TERMINAL ERROR 5 from UMACH.  The absolute', &
                     /, ' ***          value of the index variable must be' &
                     , /, ' ***          1 or 2.  IABS(N) = ', I6, &
                     '.', /)
                 STOP
             !                                  CHECK FOR RESET OR RETRIEVAL
             ELSE IF (N .LT. 0) THEN
                 !                                  RESET
                 UNIT(NN) = NUNIT
             ELSE
                 !                                  RETRIEVE
                 NUNIT = UNIT(N)
             END IF
             !
             RETURN
         END

         !-----------------------------------------------------------------------
         !  IMSL Name:  C1TCI
         !
         !  Computer:   pcdsms/SINGLE
         !
         !  Revised:    August 13, 1984
         !
         !  Purpose:    Convert character string into corresponding integer
         !              form.
         !
         !  Usage:      CALL C1TCI (CHRSTR, SLEN, NUM, IER)
         !
         !  Arguments:
         !   CHRSTR  - Character array that contains the number description.
         !             (Input)
         !   SLEN    - Length of the character array.  (Input)
         !   NUM     - The answer.  (Output)
         !   IER     - Completion code.  (Output)  Where
         !                IER =-2  indicates that the number is too large to
         !                         be converted;
         !                IER =-1  indicates that SLEN <= 0;
         !                IER = 0  indicates normal completion;
         !                IER > 0  indicates that the input string contains a
         !                         nonnumeric character.  IER is the index of
         !                         the first nonnumeric character in CHRSTR.
         !
         !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
         !
         !  Warranty:   IMSL warrants only that IMSL testing has been applied
         !              to this code.  No other warranty, expressed or implied,
         !              is applicable.
         !
         !-----------------------------------------------------------------------
         !
         SUBROUTINE C1TCI (CHRSTR, SLEN, NUM, IER)
             !                                  SPECIFICATIONS FOR ARGUMENTS
             INTEGER    SLEN, NUM, IER
             CHARACTER  CHRSTR(*)
             !                                  SPECIFICATIONS FOR LOCAL VARIABLES
             INTEGER    COUNT, I, IMACH5, J, N, S, SIGN
             CHARACTER  ZERO
             !                                  SPECIFICATIONS FOR SAVE VARIABLES
             CHARACTER  BLANK, DIGIT*10, MINUS, PLUS
             SAVE       BLANK, DIGIT, MINUS, PLUS
             !                                  SPECIFICATIONS FOR EQUIVALENCE
             EQUIVALENCE (DIGIT, ZERO)
             !                                  SPECIFICATIONS FOR INTRINSICS
             !     INTRINSIC  INDEX
             INTRINSIC  INDEX
             INTEGER    INDEX
             !                                  SPECIFICATIONS FOR FUNCTIONS
             EXTERNAL   IMACH
             INTEGER    IMACH
             !
             DATA DIGIT/'0123456789'/
             DATA BLANK/' '/, MINUS/'-'/, PLUS/'+'/
             !
             !                                  CHECK SLEN
             NUM = 0
             IF (SLEN .LE. 0) THEN
                 IER = -1
                 GO TO 50
             END IF
             !                                  HANDLE LEADING BLANKS
             SIGN = 1
             I = 1
10           IF (I .LE. SLEN) THEN
                 IF (CHRSTR(I) .EQ. BLANK) THEN
                     I = I + 1
                     GO TO 10
                 END IF
             ELSE
                 IER = 1
                 GO TO 50
             END IF
             !                                  CHECK FOR SIGN, IF ANY
             S = I
             IF (CHRSTR(I) .EQ. MINUS) THEN
                 SIGN = -1
                 I = I + 1
             ELSE IF (CHRSTR(I) .EQ. PLUS) THEN
                 I = I + 1
             END IF
20           IF (I .LE. SLEN) THEN
                 IF (CHRSTR(I) .EQ. BLANK) THEN
                     I = I + 1
                     GO TO 20
                 END IF
             ELSE
                 IER = S
                 GO TO 50
             END IF
             !                                  SKIP LEADING ZERO
             J = I
30           IF (I .LE. SLEN) THEN
                 IF (CHRSTR(I) .EQ. ZERO) THEN
                     I = I + 1
                     GO TO 30
                 END IF
             ELSE
                 IER = 0
                 GO TO 50
             END IF
             !                                  CHECK FIRST NONBLANK CHARACTER
             COUNT = 0
             !                                  CHECK NUMERIC CHARACTERS
             IMACH5 = IMACH(5)
40           N = INDEX(DIGIT,CHRSTR(I))
             IF (N .NE. 0) THEN
                 COUNT = COUNT + 1
                 IF (NUM .GT. ((IMACH5-N)+1)/10) THEN
                     IER = -2
                     GO TO 50
                 ELSE
                     NUM = NUM*10 - 1 + N
                     I = I + 1
                     IF (I .LE. SLEN) GO TO 40
                 END IF
             END IF
             !
             IF (COUNT .EQ. 0) THEN
                 IF (I .GT. J) THEN
                     IER = I
                 ELSE
                     IER = S
                 END IF
             ELSE IF (I .GT. SLEN) THEN
                 NUM = SIGN*NUM
                 IER = 0
             ELSE
                 NUM = SIGN*NUM
                 IER = I
             END IF
         !
50       CONTINUE
         RETURN
     END
     !-----------------------------------------------------------------------
     !  IMSL Name:  C1TIC
     !
     !  Computer:   pcdsms/SINGLE
     !
     !  Revised:    March 9, 1984
     !
     !  Purpose:    Convert an integer to its corresponding character form.
     !              (Right justified)
     !
     !  Usage:      CALL C1TIC(NUM, CHRSTR, SLEN, IER)
     !
     !  Arguments:
     !     NUM    - Integer number.  (Input)
     !     CHRSTR - Character array that receives the result.  (Output)
     !     SLEN   - Length of the character array.  (Input)
     !     IER    - Completion code.  (Output) Where
     !                 IER < 0  indicates that SLEN <= 0,
     !                 IER = 0  indicates normal completion,
     !                 IER > 0  indicates that the character array is too
     !                       small to hold the complete number.  IER
     !                       indicates how many significant digits are
     !                       being truncated.
     !
     !  Remarks:
     !  1. The character array is filled in a right justified manner.
     !  2. Leading zeros are replaced by blanks.
     !  3. Sign is inserted only for negative number.
     !
     !  Copyright:  1984 by IMSL, Inc.  All rights reserved.
     !
     !  Warranty:   IMSL warrants only that IMSL testing has been applied
     !              to this code.  No other warranty, expressed or implied,
     !              is applicable.
     !
     !-----------------------------------------------------------------------
     !
     SUBROUTINE C1TIC (NUM, CHRSTR, SLEN, IER)
         !                                  SPECIFICATIONS FOR ARGUMENTS
         INTEGER    NUM, SLEN, IER
         CHARACTER  CHRSTR(*)
         !                                  SPECIFICATIONS FOR LOCAL VARIABLES
         INTEGER    I, J, K, L
         !                                  SPECIFICATIONS FOR SAVE VARIABLES
         CHARACTER  BLANK(1), DIGIT(10), MINUS(1)
         SAVE       BLANK, DIGIT, MINUS
         !                                  SPECIFICATIONS FOR INTRINSICS
         !     INTRINSIC  IABS
         INTRINSIC  IABS
         INTEGER    IABS
         !                                  SPECIFICATIONS FOR SUBROUTINES
         EXTERNAL   M1VE
         !
         DATA DIGIT/'0', '1', '2', '3', '4', '5', '6', '7', '8', &
             '9'/
         DATA BLANK/' '/, MINUS/'-'/
         !                                  CHECK SLEN
         IF (SLEN .LE. 0) THEN
             IER = -1
             RETURN
         END IF
         !                                  THE NUMBER IS ZERO
         IF (NUM .EQ. 0) THEN
             CALL M1VE (BLANK, 1, 1, 1, CHRSTR, 1, SLEN-1, SLEN, I)
             CHRSTR(SLEN) = DIGIT(1)
             IER = 0
             RETURN
         END IF
         !                                  CONVERT NUMBER DIGIT BY DIGIT TO
         !                                  CHARACTER FORM
         J = SLEN
         K = IABS(NUM)
10       IF (K.GT.0 .AND. J.GE.1) THEN
             L = K
             K = K/10
             L = L - K*10
             CHRSTR(J) = DIGIT(L+1)
             J = J - 1
             GO TO 10
         END IF
         !
20       IF (K .EQ. 0) THEN
             IF (NUM .LT. 0) THEN
                 CALL M1VE (MINUS, 1, 1, 1, CHRSTR, J, J, SLEN, I)
                 IF (I .NE. 0) THEN
                     IER = 1
                     RETURN
                 END IF
                 J = J - 1
             END IF
             IER = 0
             CALL M1VE (BLANK, 1, 1, 1, CHRSTR, 1, J, SLEN, I)
             RETURN
         END IF
         !                                  DETERMINE THE NUMBER OF SIGNIFICANT
         !                                  DIGITS BEING TRUNCATED
         I = 0
30       IF (K .GT. 0) THEN
             K = K/10
             I = I + 1
             GO TO 30
         END IF
         !
         IF (NUM .LT. 0) I = I + 1
         IER = I
         !
         RETURN
     END
