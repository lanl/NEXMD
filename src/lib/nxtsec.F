      SUBROUTINE NXTSEC(IUNIT,IOUT,IONERR,FMTOLD,FLAG,FMT,IOK)
c
c   Subroutine NeXT SECtion
c
c   This routine reads data from a new-format PARM file. It
c   searches for the section with a %FLAG header of FLAG. It returns
c   the format for the section of data and places the file pointer on
c   the first line of the data block. The actual data read is performed
c   by the calling routine.
c
c   Data are read from the file on unit IUNIT, which is assumed
c   to already be open.
c
c   IOK: 0, flag found and data read
c       -1, then no %VERSION line found. This is an old-format PARM file.
c           In this case, any call to NXTSEC will merely replace FMT with
c           FMTOLD. This simplifies the calling procedure in the main
c           routine, since FMT will contain the approprate FMT regardless
c           of whether a new or old format PARM file is used (as long
c           as FMTOLD was appropriately set in the call list).
c       -2, then this is a new-format PARM file, but the requested
c           FLAG was not found. (Only if IONERR = 1).
c           
c    Program stops if a specified flag is not found and this is a new-format
c    PARM file.
c
c   IUNIT: Unit for reads, assumed to already be open.
c   IOUT: Unit for info/error writes
c   IONERR: 0, then if a requested flag is not found, the program
c              stops with an appropriate error
c           1, then if a requested flag is not found, the routine
c              returns with IOK set to -2.
c   FMTOLD: Format to use if read takes place from an old-style PARM file
c   FLAG: Flag for data section to read. Must be large enough to hold
c         any format string. Suggested length = char*255.
c   FMT: Returned with format to use for data. File pointer will be
c        at first line of data to be read upon return
c   IOK: see above.
c
c   IOUT: Unit for error prints 
c
c   Author: David Pearlman
c   Date: 09/00
c
c   Scott Brozell June 2004
c   Converted loop control to Fortran 90; these changes are g77 compatible.
c
c   The PARM file has the following format. 
c
c   %VERSION  VERSION_STAMP = Vxxxx.yyy  DATE = mm:dd:yy hh:mm:ss 
c
c      This line should appear as the first line in the file, but this
c      is not absolutely required. A search will be made for a line starting
c      with %VERSION and followed by the VERSION_STAMP field.
c      The version stamp is expected to be an F8.3 format field with
c      leading 0's in place. E.g. V0003.22. Leading 0s should also
c      be used for mm, dd, yy, hh, mm or ss fields that are < 10.
c
c   %FLAG flag
c      This line specifies the name for the block of data to follow
c      FLAGS MUST NOT HAVE ANY EMBEDDED BLANKS. Use underscore characters
c      in place of blanks, e.g. "BOND_PARMS" not "BOND PARMS".
c   %FORMAT format
c      This line provides the FORTRAN format for the data to follow.
c      This should be specified using standard FORTRAN rules, and with the
c      surrounding brackets intact. E.g. 
c         %FORMAT (8F10.3)
c      **> Data starts with the line immediately following the %FORMAT line.
c      The %FORMAT line and the data that follow will be associated with the
c      flag on the most recent %FLAG line read. 
c      The actual data read is performed by the calling routine.
c      All text following the %FORMAT flag is considered the format string
c      and the string CAN have embedded blanks.
c   %COMMENT comment
c      Comment line. Will be ignored in parsing the file. A %COMMENT line
c      can appear anywhere in the file EXCEPT A) between the %FORMAT
c      line and the data; or B) interspersed with the data lines.
c      While it recommended you use the %COMMENT line for clarity, it is
c      not technically required for comment lines. Any line without
c      a type specifier at the beginning of the line and which does not
c      appear within the data block is assumed to be a comment.
c
c   Note that in order to avoid confusion/mistakes, the above flags must
c   be left justified (start in column one) on a line to be recognized.
c
c   On the first call to this routine, it will search the file for
c   %FLAG cards and store the lines they appear on. That way, on
c   subsequent calls we'll know immediately if we should read further
c   down the file, rewind, or exit with an error (flag not found).

      implicit none
      integer  IUNIT
      integer  IOUT
      integer  IONERR
      character*(*) FMTOLD,FMT,FLAG
      integer  IOK

      integer  NNBCHR

      logical  FIRST
      SAVE FIRST
      DATA FIRST/.TRUE./
c
c   MXNXFL is maximum number of %FLAG cards that can be specified
c
      integer  MXNXFL
      PARAMETER (MXNXFL = 500)

      CHARACTER*80 NXTFLG
      CHARACTER*8 PRDAT,PRTIM
      CHARACTER*255 AA
      integer  IBLOCK
      integer  INXTFL
      integer  IPRVRR
      integer  NUMFLG
      real     RPVER
      COMMON /NXTLC1/INXTFL(2,MXNXFL),IPRVRR,NUMFLG,IBLOCK
      COMMON /NXTLC2/RPVER
      COMMON /NXTLC3/NXTFLG(MXNXFL),PRDAT,PRTIM

! carlos added this to avoid problem with data statement initilization
! need initprmtop here.  Now located in a header file. (TL 2010/01/29)
#include "nxtsec.h"      

      integer  I
      integer  IPT
      integer  IPT2
      integer  IPT3
      integer  IPT4
      integer  IPT5
      integer  IPT6
      integer  IPT7
      integer  IPT8
      integer  IPT9
      integer  IPT10
      integer  LFLAG
      integer  IL2US
      integer  IFIND
      integer  MBLOCK
      integer  ILFO

      IOK = 0
      IF (FIRST.or.initprmtop) THEN
c
      REWIND(IUNIT)
c
c   First, see if this is a new format PARM file. That is, if the %VERSION
c   line exists. If not, then we assume it's an old format PARM file. In
c   this case, every call to NXTSEC will simply result in an immediate
c   return. This means all reads from the calling routine will be done
c   sequentially from the PARM file. Store the version number as a real
c   in RPVER. Store the date and time strings as character strings in
c   PRDAT and PRTIM.
c
      do
         READ(IUNIT,11,END=20) AA
   11    FORMAT(A)
         IF (AA(1:8).NE.'%VERSION') cycle
c
         IPT = INDEX(AA,'VERSION_STAMP')
         IF (IPT.LE.0) cycle
c
         IPT2 = NNBCHR(AA,IPT+13,0,0)
         IF (AA(IPT2:IPT2).NE.'=') GO TO 9000
c
         IPT3 = NNBCHR(AA,IPT2+1,0,0)
         IF (AA(IPT3:IPT3).NE.'V') GO TO 9001
c
         IPT4 = NNBCHR(AA,IPT3+1,0,1)
         IF (IPT4-1 - (IPT3+1) + 1 .NE.8) GO TO 9002
         READ(AA(IPT3+1:IPT4-1),'(F8.3)') RPVER
c
         IPT5 = INDEX(AA,'DATE')
         IF (IPT5.LE.0) THEN
           PRDAT = 'xx/xx/xx'
           PRTIM = 'xx:xx:xx'
           GO TO 50
         END IF
         IPT6 = NNBCHR(AA,IPT5+4,0,0)
         IF (AA(IPT6:IPT6).NE.'=') GO TO 9003
         IPT7 = NNBCHR(AA,IPT6+1,0,0)
         IPT8 = NNBCHR(AA,IPT7+1,0,1)
         IF (IPT8-1 - IPT7 + 1 .NE. 8) GO TO 9004
         PRDAT = AA(IPT7:IPT8-1)

         IPT9 = NNBCHR(AA,IPT8+1,0,0)
         IPT10 = NNBCHR(AA,IPT9+1,0,1)
         IF (IPT10-1 - IPT9 + 1 .NE. 8) GO TO 9005
         PRTIM = AA(IPT9:IPT10-1)
         WRITE(IOUT,15) RPVER,PRDAT,PRTIM
   15    FORMAT('| New format PARM file being parsed.',/,
     *          '| Version = ',F8.3,' Date = ',A,' Time = ',A)
         IPRVRR = 0
         GO TO 50
      end do
c
c   Get here if no VERSION flag read. Set IPRVRR = 1 and return.
c   On subsequent calls, if IPRVRR = 1, we return immediately.
c
   20 IPRVRR = 1
      IOK = -1
      WRITE(IOUT,21)
   21 FORMAT('|  INFO: Old style PARM file read',/)
      fmt = fmtold
      rewind(iunit)
      first = .false.
      ! write (6,*) "setting initprmtop to F"
      initprmtop=.false.

      RETURN
c
c   %VERSION line successfully read. Now load the flags into NXTFLG(I)
c   and the line pointer and lengths of the flags into 
c   INXTFL(1,I) and INXTFL(2,I), respectively. NUMFLG will be the 
c   total number of flags read.
c
   50 REWIND(IUNIT)

      NUMFLG = 0
      I = 1
      do
         READ(IUNIT,11,END=99) AA
         IF (AA(1:5).EQ.'%FLAG') THEN
           NUMFLG = NUMFLG + 1
           IPT2 = NNBCHR(AA,6,0,0)
           IF (IPT2.EQ.-1) GO TO 9006
           IPT3 = NNBCHR(AA,IPT2,0,1)-1

           INXTFL(1,NUMFLG) = I
           INXTFL(2,NUMFLG) = IPT3-IPT2+1
           NXTFLG(NUMFLG) = AA(IPT2:IPT3)
         END IF
         I = I + 1
      end do
   99 REWIND(IUNIT)
      IBLOCK = 0
      FIRST = .FALSE.
      initprmtop=.false.
      END IF
c
c   Start search for passed flag name
c
c   If this is an old-style PARM file, we can't do the search. Simply
c   set IOK = -1, FMT to FMTOLD, and return
c
      IF (IPRVRR.EQ.1) THEN
         IOK = -1
         FMT = FMTOLD
         RETURN
      END IF
c
      LFLAG = NNBCHR(FLAG,1,0,1)-1
      IF (LFLAG.EQ.-2) LFLAG = LEN(FLAG)
      DO I = 1,NUMFLG
         IF (LFLAG.EQ.INXTFL(2,I)) THEN
            IF (FLAG(1:LFLAG).EQ.NXTFLG(I)(1:LFLAG)) THEN
               IL2US = INXTFL(1,I)
               GO TO 120
            END IF
         END IF
      END DO
c
c   Get here if flag does not correspond to any stored. Either stop
c   or return depending on IONERR flag.
c
      IF (IONERR.EQ.0) THEN
         GO TO 9007
      ELSE IF (IONERR.EQ.1) THEN
         IOK = -2
         RETURN
      END IF
c
c   Flag found. Set file pointer to the first line following the appropriate
c   %FLAG line and then search for %FORMAT field.
c
c   IBLOCK keeps track of the last %FLAG found. If this preceeded the
c   one being read now, we read forward to find the current requested FLAG.
c   If this followed the current request, rewind and read forward the
c   necessary number of lines. This should speed things up a bit.
c   
  120 IFIND = I
      MBLOCK = IBLOCK
      IF (IFIND.GT.IBLOCK) THEN
         do
            READ(IUNIT,11,END=9008) AA
            IF (AA(1:5).EQ.'%FLAG') THEN
               MBLOCK = MBLOCK + 1
               IF (MBLOCK.EQ.IFIND) exit
            END IF
         end do
      ELSE
         REWIND(IUNIT)
         DO I = 1,IL2US
            READ(IUNIT,11,END=9008)
         END DO
      END IF

      DO
         READ(IUNIT,11,END=9009) AA
         IF (AA(1:7).EQ.'%FORMAT') exit
      END DO
c
c   First %FORMAT found following appropriate %FLAG. Extract the
c   format and return. All non-blank characters following %FORMAT
c   comprise the format string (embedded blanks allowed).
c
      IPT2 = NNBCHR(AA,8,0,0)
      IF (IPT2.EQ.-1) GO TO 9010
      DO I = LEN(AA),IPT2,-1
         IF (AA(I:I).NE.' ') exit
      END DO
      IPT3 = I
c
c   Format string is in IPT2:IPT3. Make sure passed FMT string is large
c   enought to hold this and then return.
c
      ILFO = IPT3-IPT2+1
      IF (ILFO.GT.LEN(FMT)) GO TO 9011
      FMT = ' '
      FMT(1:ILFO) = AA(IPT2:IPT3)
c
c   Update IBLOCK pointer and return
c
      IBLOCK = IFIND
      RETURN
c
c   Errors:
c
 9000 WRITE(IOUT,9500)
 9500 FORMAT('ERROR: No = sign after VERSION_STAMP field in PARM')
      STOP
 9001 WRITE(IOUT,9501)
 9501 FORMAT('ERROR: Version number in PARM does not start with V')
      STOP
 9002 WRITE(IOUT,9502)
 9502 FORMAT('ERROR: Mal-formed version number in PARM. ',
     *       'Should be 8 chars')    
      STOP
 9003 WRITE(IOUT,9503)
 9503 FORMAT('ERROR: No = sign after DATE field in PARM')
      STOP
 9004 WRITE(IOUT,9504)
 9504 FORMAT('ERROR: Mal-formed date string in PARM. ',
     *       'Should be 8 characters & no embedded spaces.')
      STOP
 9005 WRITE(IOUT,9505)
 9505 FORMAT('ERROR: Mal-formed time string in PARM. ',
     *       'Should be 8 characters & no embedded spaces.')
      STOP
 9006 WRITE(IOUT,9506)
 9506 FORMAT('ERROR: No flag found following a %FLAG line in PARM')
      STOP
 9007 WRITE(IOUT,9507) FLAG(1:LFLAG)
 9507 FORMAT('ERROR: Flag "',A,'" not found in PARM file')
      STOP
 9008 WRITE(IOUT,9508) FLAG(1:LFLAG)
 9508 FORMAT('ERROR: Programming error in routine NXTSEC at "',A,'"')
      STOP
 9009 WRITE(IOUT,9509) FLAG(1:LFLAG)
 9509 FORMAT('ERROR: No %FORMAT field found following flag "',A,'"')
      STOP
 9010 WRITE(IOUT,9510) FLAG(1:LFLAG)
 9510 FORMAT('ERROR: No format string found following a %FORMAT ',
     *       'line in PARM',/,
     *       'Corresponding %FLAG is "',A,'"')
      STOP
 9011 WRITE(IOUT,9511) FLAG(1:LFLAG)
 9511 FORMAT('ERROR: Format string for flag "',A,'" too large',/,
     *       '       for FMT call-list parameter')
      STOP
c
      END
c
      FUNCTION NNBCHR(AA,IBEG,IEND,IOPER)
c
c   IOPER = 0: Find next non-blank character
c   IOPER = 1: Find next blank character
c
c   On return, NNBCHR is set to the appropriate pointer, or to -1
c      if no non-blank character found (IOPER = 0) or no blank
c      character found (IOPER = 1).
c
      implicit none
      integer  NNBCHR
      character*(*) AA
      integer  IBEG
      integer  IEND
      integer  IOPER

      integer  I
      integer  IBG
      integer  IEN

      IBG = IBEG
      IEN = IEND
      IF (IBEG.LE.0) IBG = 1
      IF (IEND.LE.0) IEN = LEN(AA)
c
      IF (IOPER.EQ.0) THEN
         DO I = IBG,IEN
           IF (AA(I:I).NE.' ') THEN
             NNBCHR = I
             RETURN
           END IF
         end do
         NNBCHR = -1
      ELSE IF (IOPER.EQ.1) THEN
         do I = IBG,IEN
           IF (AA(I:I).EQ.' ') THEN
             NNBCHR = I
             RETURN
           END IF
         end do
         NNBCHR = -1
      END IF
c
      RETURN
      END
!--------------------------------------------------------------
      SUBROUTINE NXTSEC_crd(IUNIT,IOUT,IONERR,FMTOLD,FLAG,FMT,IOK)
c
c   Subroutine NeXT SECtion
c
c   This routine reads data from a new-format COORD file. It
c   searches for the section with a %FLAG header of FLAG. It returns
c   the format for the section of data and places the file pointer on
c   the first line of the data block. The actual data read is performed
c   by the calling routine.
c
c   Data are read from the file on unit IUNIT, which is assumed
c   to already be open.
c
c   IOK: 0, flag found and data read
c       -1, then no %VERSION line found. This is an old-format COORD file.
c           In this case, any call to NXTSEC will merely replace FMT with
c           FMTOLD. This simplifies the calling procedure in the main
c           routine, since FMT will contain the approprate FMT regardless
c           of whether a new or old format COORD file is used (as long
c           as FMTOLD was appropriately set in the call list).
c       -2, then this is a new-format COORD file, but the requested
c           FLAG was not found. (Only if IONERR = 1).
c           
c    Program stops if a specified flag is not found and this is a new-format
c    COORD file.
c
c   IUNIT: Unit for reads, assumed to already be open.
c   IOUT: Unit for info/error writes
c   IONERR: 0, then if a requested flag is not found, the program
c              stops with an appropriate error
c           1, then if a requested flag is not found, the routine
c              returns with IOK set to -2.
c   FMTOLD: Format to use if read takes place from an old-style COORD file
c   FLAG: Flag for data section to read. Must be large enough to hold
c         any format string. Suggested length = char*255.
c   FMT: Returned with format to use for data. File pointer will be
c        at first line of data to be read upon return
c   IOK: see above.
c
c   IOUT: Unit for error prints 
c
c   Author: David Pearlman
c   Date: 09/00
c
c   Scott Brozell June 2004
c   Converted loop control to Fortran 90; these changes are g77 compatible.
c
c   The COORD file has the following format. 
c
c   %VERSION  VERSION_STAMP = Vxxxx.yyy  DATE = mm:dd:yy hh:mm:ss 
c
c      This line should appear as the first line in the file, but this
c      is not absolutely required. A search will be made for a line starting
c      with %VERSION and followed by the VERSION_STAMP field.
c      The version stamp is expected to be an F8.3 format field with
c      leading 0's in place. E.g. V0003.22. Leading 0s should also
c      be used for mm, dd, yy, hh, mm or ss fields that are < 10.
c
c   %FLAG flag
c      This line specifies the name for the block of data to follow
c      FLAGS MUST NOT HAVE ANY EMBEDDED BLANKS. Use underscore characters
c      in place of blanks, e.g. "BOND_COORDS" not "BOND COORDS".
c   %FORMAT format
c      This line provides the FORTRAN format for the data to follow.
c      This should be specified using standard FORTRAN rules, and with the
c      surrounding brackets intact. E.g. 
c         %FORMAT (8F10.3)
c      **> Data starts with the line immediately following the %FORMAT line.
c      The %FORMAT line and the data that follow will be associated with the
c      flag on the most recent %FLAG line read. 
c      The actual data read is performed by the calling routine.
c      All text following the %FORMAT flag is considered the format string
c      and the string CAN have embedded blanks.
c   %COMMENT comment
c      Comment line. Will be ignored in parsing the file. A %COMMENT line
c      can appear anywhere in the file EXCEPT A) between the %FORMAT
c      line and the data; or B) interspersed with the data lines.
c      While it recommended you use the %COMMENT line for clarity, it is
c      not technically required for comment lines. Any line without
c      a type specifier at the beginning of the line and which does not
c      appear within the data block is assumed to be a comment.
c
c   Note that in order to avoid confusion/mistakes, the above flags must
c   be left justified (start in column one) on a line to be recognized.
c
c   On the first call to this routine, it will search the file for
c   %FLAG cards and store the lines they appear on. That way, on
c   subsequent calls we'll know immediately if we should read further
c   down the file, rewind, or exit with an error (flag not found).

      implicit none
      integer  IUNIT
      integer  IOUT
      integer  IONERR
      character*(*) FMTOLD,FMT,FLAG
      integer  IOK

      integer  NNBCHR

      logical  FIRST
      save FIRST
      DATA FIRST/.TRUE./
c
c   MXNXFL is maximum number of %FLAG cards that can be specified
c
      integer  MXNXFL
      PARAMETER (MXNXFL = 500)

      CHARACTER*80 NXTFLG
      CHARACTER*8 PRDAT,PRTIM
      CHARACTER*255 AA
      integer  IBLOCK
      integer  INXTFL
      integer  IPRVRR
      integer  NUMFLG
      real     RPVER
      COMMON /NXTLC1_crd/INXTFL(2,MXNXFL),IPRVRR,NUMFLG,IBLOCK
      COMMON /NXTLC2_crd/RPVER
      COMMON /NXTLC3_crd/NXTFLG(MXNXFL),PRDAT,PRTIM

      integer  I
      integer  IPT
      integer  IPT2
      integer  IPT3
      integer  IPT4
      integer  IPT5
      integer  IPT6
      integer  IPT7
      integer  IPT8
      integer  IPT9
      integer  IPT10
      integer  LFLAG
      integer  IL2US
      integer  IFIND
      integer  MBLOCK
      integer  ILFO

      IOK = 0
      IF (FIRST) THEN
c
      REWIND(IUNIT)
c
c   First, see if this is a new format COORD file. That is, if the %VERSION
c   line exists. If not, then we assume it's an old format COORD file. In
c   this case, every call to NXTSEC will simply result in an immediate
c   return. This means all reads from the calling routine will be done
c   sequentially from the COORD file. Store the version number as a real
c   in RPVER. Store the date and time strings as character strings in
c   PRDAT and PRTIM.
c
      do
         READ(IUNIT,11,END=20) AA
   11    FORMAT(A)
         IF (AA(1:8).NE.'%VERSION') cycle
c
         IPT = INDEX(AA,'VERSION_STAMP')
         IF (IPT.LE.0) cycle
c
         IPT2 = NNBCHR(AA,IPT+13,0,0)
         IF (AA(IPT2:IPT2).NE.'=') GO TO 9000
c
         IPT3 = NNBCHR(AA,IPT2+1,0,0)
         IF (AA(IPT3:IPT3).NE.'V') GO TO 9001
c
         IPT4 = NNBCHR(AA,IPT3+1,0,1)
         IF (IPT4-1 - (IPT3+1) + 1 .NE.8) GO TO 9002
         READ(AA(IPT3+1:IPT4-1),'(F8.3)') RPVER
c
         IPT5 = INDEX(AA,'DATE')
         IF (IPT5.LE.0) THEN
           PRDAT = 'xx/xx/xx'
           PRTIM = 'xx:xx:xx'
           GO TO 50
         END IF
         IPT6 = NNBCHR(AA,IPT5+4,0,0)
         IF (AA(IPT6:IPT6).NE.'=') GO TO 9003
         IPT7 = NNBCHR(AA,IPT6+1,0,0)
         IPT8 = NNBCHR(AA,IPT7+1,0,1)
         IF (IPT8-1 - IPT7 + 1 .NE. 8) GO TO 9004
         PRDAT = AA(IPT7:IPT8-1)

         IPT9 = NNBCHR(AA,IPT8+1,0,0)
         IPT10 = NNBCHR(AA,IPT9+1,0,1)
         IF (IPT10-1 - IPT9 + 1 .NE. 8) GO TO 9005
         PRTIM = AA(IPT9:IPT10-1)
         WRITE(IOUT,15) RPVER,PRDAT,PRTIM
   15    FORMAT('| New format inpcrd file being parsed.',/,
     *          '| Version = ',F8.3,' Date = ',A,' Time = ',A)
         IPRVRR = 0
         GO TO 50
      end do
c
c   Get here if no VERSION flag read. Set IPRVRR = 1 and return.
c   On subsequent calls, if IPRVRR = 1, we return immediately.
c
   20 IPRVRR = 1
      IOK = -1
      WRITE(IOUT,21)
   21 FORMAT('|  INFO: Old style inpcrd file read',/)
      fmt = fmtold
      rewind(iunit)
      first = .false.
      RETURN
c
c   %VERSION line successfully read. Now load the flags into NXTFLG(I)
c   and the line pointer and lengths of the flags into 
c   INXTFL(1,I) and INXTFL(2,I), respectively. NUMFLG will be the 
c   total number of flags read.
c
   50 REWIND(IUNIT)
      NUMFLG = 0
      I = 1
      do
         READ(IUNIT,11,END=99) AA
         IF (AA(1:5).EQ.'%FLAG') THEN
           NUMFLG = NUMFLG + 1
           IPT2 = NNBCHR(AA,6,0,0)
           IF (IPT2.EQ.-1) GO TO 9006
           IPT3 = NNBCHR(AA,IPT2,0,1)-1

           INXTFL(1,NUMFLG) = I
           INXTFL(2,NUMFLG) = IPT3-IPT2+1
           NXTFLG(NUMFLG) = AA(IPT2:IPT3)
         END IF
         I = I + 1
      end do
   99 REWIND(IUNIT)
      IBLOCK = 0
      FIRST = .FALSE.
      END IF
c
c   Start search for passed flag name
c
c   If this is an old-style COORD file, we can't do the search. Simply
c   set IOK = -1, FMT to FMTOLD, and return
c
      IF (IPRVRR.EQ.1) THEN
         IOK = -1
         FMT = FMTOLD
         RETURN
      END IF
c
      LFLAG = NNBCHR(FLAG,1,0,1)-1
      IF (LFLAG.EQ.-2) LFLAG = LEN(FLAG)
      DO I = 1,NUMFLG
         IF (LFLAG.EQ.INXTFL(2,I)) THEN
            IF (FLAG(1:LFLAG).EQ.NXTFLG(I)(1:LFLAG)) THEN
               IL2US = INXTFL(1,I)
               GO TO 120
            END IF
         END IF
      END DO
c
c   Get here if flag does not correspond to any stored. Either stop
c   or return depending on IONERR flag.
c
      IF (IONERR.EQ.0) THEN
         GO TO 9007
      ELSE IF (IONERR.EQ.1) THEN
         IOK = -2
         RETURN
      END IF
c
c   Flag found. Set file pointer to the first line following the appropriate
c   %FLAG line and then search for %FORMAT field.
c
c   IBLOCK keeps track of the last %FLAG found. If this preceeded the
c   one being read now, we read forward to find the current requested FLAG.
c   If this followed the current request, rewind and read forward the
c   necessary number of lines. This should speed things up a bit.
c   
  120 IFIND = I
      MBLOCK = IBLOCK
      IF (IFIND.GT.IBLOCK) THEN
         do
            READ(IUNIT,11,END=9008) AA
            IF (AA(1:5).EQ.'%FLAG') THEN
               MBLOCK = MBLOCK + 1
               IF (MBLOCK.EQ.IFIND) exit
            END IF
         end do
      ELSE
         REWIND(IUNIT)
         DO I = 1,IL2US
            READ(IUNIT,11,END=9008)
         END DO
      END IF

      DO
         READ(IUNIT,11,END=9009) AA
         IF (AA(1:7).EQ.'%FORMAT') exit
      END DO
c
c   First %FORMAT found following appropriate %FLAG. Extract the
c   format and return. All non-blank characters following %FORMAT
c   comprise the format string (embedded blanks allowed).
c
      IPT2 = NNBCHR(AA,8,0,0)
      IF (IPT2.EQ.-1) GO TO 9010
      DO I = LEN(AA),IPT2,-1
         IF (AA(I:I).NE.' ') exit
      END DO
      IPT3 = I
c
c   Format string is in IPT2:IPT3. Make sure passed FMT string is large
c   enought to hold this and then return.
c
      ILFO = IPT3-IPT2+1
      IF (ILFO.GT.LEN(FMT)) GO TO 9011
      FMT = ' '
      FMT(1:ILFO) = AA(IPT2:IPT3)
c
c   Update IBLOCK pointer and return
c
      IBLOCK = IFIND
      RETURN
c
c   Errors:
c
 9000 WRITE(IOUT,9500)
 9500 FORMAT('ERROR: No = sign after VERSION_STAMP field in COORD')
      STOP
 9001 WRITE(IOUT,9501)
 9501 FORMAT('ERROR: Version number in COORD does not start with V')
      STOP
 9002 WRITE(IOUT,9502)
 9502 FORMAT('ERROR: Mal-formed version number in COORD. ',
     *       'Should be 8 chars')    
      STOP
 9003 WRITE(IOUT,9503)
 9503 FORMAT('ERROR: No = sign after DATE field in COORD')
      STOP
 9004 WRITE(IOUT,9504)
 9504 FORMAT('ERROR: Mal-formed date string in COORD. ',
     *       'Should be 8 characters & no embedded spaces.')
      STOP
 9005 WRITE(IOUT,9505)
 9505 FORMAT('ERROR: Mal-formed time string in COORD. ',
     *       'Should be 8 characters & no embedded spaces.')
      STOP
 9006 WRITE(IOUT,9506)
 9506 FORMAT('ERROR: No flag found following a %FLAG line in COORD')
      STOP
 9007 WRITE(IOUT,9507) FLAG(1:LFLAG)
 9507 FORMAT('ERROR: Flag "',A,'" not found in COORD file')
      STOP
 9008 WRITE(IOUT,9508) FLAG(1:LFLAG)
 9508 FORMAT('ERROR: Programming error in routine NXTSEC at "',A,'"')
      STOP
 9009 WRITE(IOUT,9509) FLAG(1:LFLAG)
 9509 FORMAT('ERROR: No %FORMAT field found following flag "',A,'"')
      STOP
 9010 WRITE(IOUT,9510) FLAG(1:LFLAG)
 9510 FORMAT('ERROR: No format string found following a %FORMAT ',
     *       'line in COORD',/,
     *       'Corresponding %FLAG is "',A,'"')
      STOP
 9011 WRITE(IOUT,9511) FLAG(1:LFLAG)
 9511 FORMAT('ERROR: Format string for flag "',A,'" too large',/,
     *       '       for FMT call-list parameter')
      STOP
c
      END

      subroutine nxtsec_crd_reset()
      implicit none
      integer  INXTFL
      integer  MXNXFL
      PARAMETER (MXNXFL = 500)
      integer  IPRVRR
      integer  NUMFLG
      integer  IBLOCK
      COMMON /NXTLC1_crd/INXTFL(2,MXNXFL),IPRVRR,NUMFLG,IBLOCK
      
      iblock = 0
      return
      end
