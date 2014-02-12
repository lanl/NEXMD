!DECK DCHKW
      SUBROUTINE DCHKW (NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR)
!***BEGIN PROLOGUE  DCHKW
!***SUBSIDIARY
!***PURPOSE  SLAP WORK/IWORK Array Bounds Checker.
!            This routine checks the work array lengths and interfaces
!            to the SLATEC error handler if a problem is found.
!***LIBRARY   SLATEC (SLAP)
!***CATEGORY  R2
!***TYPE      DOUBLE PRECISION (SCHKW-S, DCHKW-D)
!***KEYWORDS  ERROR CHECKING, SLAP, WORKSPACE CHECKING
!***AUTHOR  Seager, Mark K., (LLNL)
!             Lawrence Livermore National Laboratory
!             PO BOX 808, L-60
!             Livermore, CA 94550 (510) 423-3141
!             seager@llnl.gov
!***DESCRIPTION
!
! *Usage:
!     CHARACTER*(*) NAME
!     INTEGER LOCIW, LENIW, LOCW, LENW, IERR, ITER
!     DOUBLE PRECISION ERR
!
!     CALL DCHKW( NAME, LOCIW, LENIW, LOCW, LENW, IERR, ITER, ERR )
!
! *Arguments:
! NAME   :IN       Character*(*).
!         Name of the calling routine.  This is used in the output
!         message, if an error is detected.
! LOCIW  :IN       Integer.
!         Location of the first free element in the integer workspace
!         array.
! LENIW  :IN       Integer.
!         Length of the integer workspace array.
! LOCW   :IN       Integer.
!         Location of the first free element in the double precision
!         workspace array.
! LENRW  :IN       Integer.
!         Length of the double precision workspace array.
! IERR   :OUT      Integer.
!         Return error flag.
!               IERR = 0 => All went well.
!               IERR = 1 => Insufficient storage allocated for
!                           WORK or IWORK.
! ITER   :OUT      Integer.
!         Set to zero on return.
! ERR    :OUT      Double Precision.
!         Set to the smallest positive magnitude if all went well.
!         Set to a very large number if an error is detected.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   880225  DATE WRITTEN
!   881213  Previous REVISION DATE
!   890915  Made changes requested at July 1989 CML Meeting.  (MKS)
!   890922  Numerous changes to prologue to make closer to SLATEC
!           standard.  (FNF)
!   890929  Numerous changes to reduce SP/DP differences.  (FNF)
!   900805  Changed XERRWV calls to calls to XERMSG.  (RWC)
!   910411  Prologue converted to Version 4.0 format.  (BAB)
!   910502  Corrected XERMSG calls to satisfy Section 6.2.2 of ANSI
!           X3.9-1978.  (FNF)
!   910506  Made subsidiary.  (FNF)
!   920511  Added complete declaration section.  (WRB)
!   921015  Added code to initialize ITER and ERR when IERR=0.  (FNF)
!***END PROLOGUE  DCHKW
!     .. Scalar Arguments ..
      DOUBLE PRECISION ERR
      INTEGER IERR, ITER, LENIW, LENW, LOCIW, LOCW
      CHARACTER NAME*(*)
!     .. Local Scalars ..
      CHARACTER XERN1*8, XERN2*8, XERNAM*8
!     .. External Functions ..
      DOUBLE PRECISION D1MACH
      EXTERNAL D1MACH
!     .. External Subroutines ..
      EXTERNAL XERMSG
!***FIRST EXECUTABLE STATEMENT  DCHKW
!
!         Check the Integer workspace situation.
!
      IERR = 0
      ITER = 0
      ERR = D1MACH(1)
      IF( LOCIW.GT.LENIW ) THEN
         IERR = 1
         ERR = D1MACH(2)
         XERNAM = NAME
         WRITE (XERN1, '(I8)') LOCIW
         WRITE (XERN2, '(I8)') LENIW
         CALL XERMSG ('SLATEC', 'DCHKW', &
            'In ' // XERNAM // ', INTEGER work array too short.  ' // &
            'IWORK needs ' // XERN1 // '; have allocated ' // XERN2, &
            1, 1)
      ENDIF
!
!         Check the Double Precision workspace situation.
      IF( LOCW.GT.LENW ) THEN
         IERR = 1
         ERR = D1MACH(2)
         XERNAM = NAME
         WRITE (XERN1, '(I8)') LOCW
         WRITE (XERN2, '(I8)') LENW
         CALL XERMSG ('SLATEC', 'DCHKW', &
            'In ' // XERNAM // ', DOUBLE PRECISION work array too ' // &
            'short.  RWORK needs ' // XERN1 // '; have allocated ' // &
            XERN2, 1, 1)
      ENDIF
      RETURN
!------------- LAST LINE OF DCHKW FOLLOWS ----------------------------
      END
