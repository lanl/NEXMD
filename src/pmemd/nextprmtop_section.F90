!*******************************************************************************!
! Module: nextprmtop_section_mod
!
! Description: <TBS>
!
!*******************************************************************************                                                                                
module nextprmtop_section_mod

  integer, save         :: iblock

  private nnbchr, iblock

contains

!*******************************************************************************!
! Subroutine:  nxtsec (NeXT SECtion)
!
! Description:
!
!   This routine reads data from a new-format PARM file. It
!   searches for the section with a %FLAG header of FLAG. It returns
!   the format for the section of data and places the file pointer on
!   the first line of the data block. The actual data read is performed
!   by the calling routine.
!
!   Data are read from the file on unit IUNIT, which is assumed
!   to already be open.
!
!   IOK: 0, flag found and data read
!       -1, then no %VERSION line found. This is an old-format PARM file.
!           In this case, any call to NXTSEC will merely replace FMT with
!           FMTOLD. This simplifies the calling procedure in the main
!           routine, since FMT will contain the approprate FMT regardless
!           of whether a new or old format PARM file is used (as long
!           as FMTOLD was appropriately set in the call list).
!       -2, then this is a new-format PARM file, but the requested
!           FLAG was not found. (Only if IONERR = 1).
!           
!    Program stops if a specified flag is not found and this is a new-format
!    PARM file.
!
!   IUNIT: Unit for reads, assumed to already be open.
!   IOUT: Unit for info/error writes
!   IONERR: 0, then if a requested flag is not found, the program
!              stops with an appropriate error
!           1, then if a requested flag is not found, the routine
!              returns with IOK set to -2.
!   FMTOLD: Format to use if read takes place from an old-style PARM file
!   FLAG: Flag for data section to read. Must be large enough to hold
!         any format string. Suggested length = char*255.
!   FMT: Returned with format to use for data. File pointer will be
!        at first line of data to be read upon return
!   IOK: see above.
!
!   IOUT: Unit for error prints 
!
!   Author: David Pearlman
!   Date: 09/00
!
!   Scott Brozell June 2004
!   Converted loop control to Fortran 90; these changes are g77 compatible.
!
!   The PARM file has the following format. 
!
!   %VERSION  VERSION_STAMP = Vxxxx.yyy  DATE = mm:dd:yy hh:mm:ss 
!
!      This line should appear as the first line in the file, but this
!      is not absolutely required. A search will be made for a line starting
!      with %VERSION and followed by the VERSION_STAMP field.
!      The version stamp is expected to be an F8.3 format field with
!      leading 0's in place. E.g. V0003.22. Leading 0s should also
!      be used for mm, dd, yy, hh, mm or ss fields that are < 10.
!
!   %FLAG flag
!      This line specifies the name for the block of data to follow
!      FLAGS MUST NOT HAVE ANY EMBEDDED BLANKS. Use underscore characters
!      in place of blanks, e.g. "BOND_PARMS" not "BOND PARMS".
!   %FORMAT format
!      This line provides the FORTRAN format for the data to follow.
!      This should be specified using standard FORTRAN rules, and with the
!      surrounding brackets intact. E.g. 
!         %FORMAT (8F10.3)
!      **> Data starts with the line immediately following the %FORMAT line.
!      The %FORMAT line and the data that follow will be associated with the
!      flag on the most recent %FLAG line read. 
!      The actual data read is performed by the calling routine.
!      All text following the %FORMAT flag is considered the format string
!      and the string CAN have embedded blanks.
!   %COMMENT comment
!      Comment line. Will be ignored in parsing the file. A %COMMENT line
!      can appear anywhere in the file EXCEPT A) between the %FORMAT
!      line and the data; or B) interspersed with the data lines.
!      While it recommended you use the %COMMENT line for clarity, it is
!      not technically required for comment lines. Any line without
!      a type specifier at the beginning of the line and which does not
!      appear within the data block is assumed to be a comment.
!
!   Note that in order to avoid confusion/mistakes, the above flags must
!   be left justified (start in column one) on a line to be recognized.
!
!   On the first call to this routine, it will search the file for
!   %FLAG cards and store the lines they appear on. That way, on
!   subsequent calls we'll know immediately if we should read further
!   down the file, rewind, or exit with an error (flag not found).
!
!*******************************************************************************

subroutine nxtsec(iunit, iout, ionerr, fmtold, flag, fmt, iok)

  implicit none

  ! Formal arguments:

  integer       :: iunit
  integer       :: iout
  integer       :: ionerr
  character*(*) :: fmtold,fmt,flag
  integer       :: iok

  ! Local variables:

! mxnxfl is maximum number of %flag cards that can be specified

  integer, parameter    :: mxnxfl = 500

  integer               :: i
  integer               :: ipt, ipt2, ipt3, ipt4, ipt5, ipt6, ipt7, ipt8, &
                           ipt9, ipt10
  integer               :: lflag
  integer               :: il2us
  integer               :: ifind
  integer               :: mblock
  integer               :: ilfo

  character*255         :: aa
  character*80, save    :: nxtflg(mxnxfl)
  character*8, save     :: prdat, prtim

  integer, save         :: inxtfl(2, mxnxfl), iprvrr, numflg
  logical, save         :: first = .true.
  real, save            :: rpver

  iok = 0

  if (first) then

    rewind(iunit)

    ! First, see if this is a new format PARM file. That is, if the %VERSION
    ! line exists. If not, then we assume it's an old format PARM file. In
    ! this case, every call to NXTSEC will simply result in an immediate
    ! return. This means all reads from the calling routine will be done
    ! sequentially from the PARM file. Store the version number as a real
    ! in RPVER. Store the date and time strings as character strings in
    ! PRDAT and PRTIM.

    do

      read(iunit,11,end=20) aa
11    format(a)
      if (aa(1:8).ne.'%VERSION') cycle

      ipt = index(aa,'VERSION_STAMP')
      if (ipt.le.0) cycle

      ipt2 = nnbchr(aa,ipt+13,0,0)
      if (aa(ipt2:ipt2).ne.'=') go to 9000

      ipt3 = nnbchr(aa,ipt2+1,0,0)
      if (aa(ipt3:ipt3).ne.'V') go to 9001

      ipt4 = nnbchr(aa,ipt3+1,0,1)
      if (ipt4-1 - (ipt3+1) + 1 .ne. 8) go to 9002
      read(aa(ipt3+1:ipt4-1),'(f8.3)') rpver

      ipt5 = index(aa,'DATE')
      if (ipt5.le.0) then
        prdat = 'xx/xx/xx'
        prtim = 'xx:xx:xx'
      go to 50
      end if
      ipt6 = nnbchr(aa,ipt5+4,0,0)
      if (aa(ipt6:ipt6).ne.'=') go to 9003
      ipt7 = nnbchr(aa,ipt6+1,0,0)
      ipt8 = nnbchr(aa,ipt7+1,0,1)
      if (ipt8-1 - ipt7 + 1 .ne. 8) go to 9004
      prdat = aa(ipt7:ipt8-1)

      ipt9 = nnbchr(aa,ipt8+1,0,0)
      ipt10 = nnbchr(aa,ipt9+1,0,1)
      if (ipt10-1 - ipt9 + 1 .ne. 8) go to 9005
      prtim = aa(ipt9:ipt10-1)
      write(iout,15) rpver,prdat,prtim
15    format('| New format PARM file being parsed.',/, &
             '| Version = ',F8.3,' Date = ',A,' Time = ',A)
      iprvrr = 0
      go to 50

    end do

! Get here if no VERSION flag read. Set IPRVRR = 1 and return.
! On subsequent calls, if IPRVRR = 1, we return immediately.

20  iprvrr = 1
    iok = -1
    write(iout,21)
21  format('|  INFO: Old style PARM file read',/)
    fmt = fmtold
    rewind(iunit)
    first = .false.
    return

! %VERSION line successfully read. Now load the flags into NXTFLG(I)
! and the line pointer and lengths of the flags into 
! INXTFL(1,I) and INXTFL(2,I), respectively. NUMFLG will be the 
! total number of flags read.

50  rewind(iunit)
    numflg = 0
    i = 1
    do
      read(iunit,11,end=99) aa
      if (aa(1:5).eq.'%FLAG') then
        numflg = numflg + 1
        ipt2 = nnbchr(aa,6,0,0)
        if (ipt2.eq.-1) go to 9006
        ipt3 = nnbchr(aa,ipt2,0,1)-1

        inxtfl(1,numflg) = i
        inxtfl(2,numflg) = ipt3-ipt2+1
        nxtflg(numflg) = aa(ipt2:ipt3)
      end if
      i = i + 1
    end do
99  rewind(iunit)
    iblock = 0
    first = .false.
  end if

! Start search for passed flag name
!
! If this is an old-style PARM file, we can't do the search. Simply
! set IOK = -1, FMT to FMTOLD, and return

  if (iprvrr .eq. 1) then
    iok = -1
    fmt = fmtold
    return
  end if

  lflag = nnbchr(flag,1,0,1)-1
  if (lflag.eq.-2) lflag = len(flag)
  do i = 1,numflg
    if (lflag.eq.inxtfl(2,i)) then
      if (flag(1:lflag).eq.nxtflg(i)(1:lflag)) then
        il2us = inxtfl(1,i)
        go to 120
      end if
    end if
  end do

! Get here if flag does not correspond to any stored. Either stop
! or return depending on IONERR flag.

  if (ionerr.eq.0) then
    go to 9007
  else if (ionerr.eq.1) then
    iok = -2
    return
  end if

! Flag found. Set file pointer to the first line following the appropriate
! %FLAG line and then search for %FORMAT field.
!
! IBLOCK keeps track of the last %FLAG found. If this preceeded the
! one being read now, we read forward to find the current requested FLAG.
! If this followed the current request, rewind and read forward the
! necessary number of lines. This should speed things up a bit.

120 &
  ifind = i
  mblock = iblock
  if (ifind .gt. iblock) then
    do
      read(iunit,11,end=9008) aa
      if (aa(1:5) .eq. '%FLAG') then
        mblock = mblock + 1
        if (mblock .eq. ifind) exit
      end if
    end do
  else
    rewind(iunit)
    do i = 1,il2us
      read(iunit,11,end=9008)
    end do
  end if

  do
    read(iunit,11,end=9009) aa
    if (aa(1:7).eq.'%FORMAT') exit
  end do

! First %FORMAT found following appropriate %FLAG. Extract the
! format and return. All non-blank characters following %FORMAT
! comprise the format string (embedded blanks allowed).

  ipt2 = nnbchr(aa,8,0,0)
  if (ipt2.eq.-1) go to 9010
  do i = len(aa),ipt2,-1
    if (aa(i:i).ne.' ') exit
  end do
  ipt3 = i

! Format string is in IPT2:IPT3. Make sure passed FMT string is large
! enought to hold this and then return.

  ilfo = ipt3-ipt2+1
  if (ilfo.gt.len(fmt)) go to 9011
  fmt = ' '
  fmt(1:ilfo) = aa(ipt2:ipt3)

! Update IBLOCK pointer and return

  iblock = ifind
  return

! Errors:

9000 write(iout,9500)
9500 format('ERROR: No = sign after VERSION_STAMP field in PARM')
  stop

9001 write(iout,9501)
9501 format('ERROR: Version number in PARM does not start with V')
  stop

9002 write(iout,9502)
9502 format('ERROR: Mal-formed version number in PARM. ', &
             'Should be 8 chars')    
  stop

9003 write(iout,9503)
9503 format('ERROR: No = sign after DATE field in PARM')
  stop

9004 write(iout,9504)
9504 format('ERROR: Mal-formed date string in PARM. ', &
             'Should be 8 characters & no embedded spaces.')
  stop

9005 write(iout,9505)
9505 format('ERROR: Mal-formed time string in PARM. ', &
             'Should be 8 characters & no embedded spaces.')
  stop

9006 write(iout,9506)
9506 format('ERROR: No flag found following a %FLAG line in PARM')
  stop

9007 write(iout,9507) flag(1:lflag)
9507 format('ERROR: Flag "',A,'" not found in PARM file')
  stop

9008 write(iout,9508) flag(1:lflag)
9508 format('ERROR: Programming error in routine NXTSEC')
  stop

9009 write(iout,9509) flag(1:lflag)
9509 format('ERROR: No %FORMAT field found following flag "',a,'"')
  stop

9010 write(iout,9510) flag(1:lflag)
9510 format('ERROR: No format string found following a %FORMAT ', &
             'line in PARM',/, &
             'Corresponding %FLAG is "',a,'"')
  stop

9011 write(iout,9511) flag(1:lflag)
9511 format('ERROR: Format string for flag "',a,'" too large',/, &
             '       for FMT call-list parameter')
  stop

end subroutine nxtsec


!*******************************************************************************!
! Function:  nnbchr
!
! Description:
!
! IOPER = 0: Find next non-blank character
! IOPER = 1: Find next blank character
!
! On return, NNBCHR is set to the appropriate pointer, or to -1
!   if no non-blank character found (IOPER = 0) or no blank
!   character found (IOPER = 1).
!*******************************************************************************

function nnbchr(aa, ibeg, iend, ioper)

  implicit none

  ! Formal arguments:

  character*(*) :: aa
  integer       :: ibeg
  integer       :: iend
  integer       :: ioper

  ! Local variables:

  integer       :: nnbchr
  integer       :: i
  integer       :: ibg
  integer       :: ien

  ibg = ibeg
  ien = iend
  if (ibeg .le. 0) ibg = 1
  if (iend .le. 0) ien = len(aa)

  if (ioper .eq. 0) then
    do i = ibg,ien
      if (aa(i:i).ne.' ') then
        nnbchr = i
        return
      end if
    end do
    nnbchr = -1
  else if (ioper .eq. 1) then
    do i = ibg,ien
      if (aa(i:i) .eq. ' ') then
        nnbchr = i
        return
      end if
    end do
    nnbchr = -1
  end if

  return

end function nnbchr

!*******************************************************************************
!
! Subroutine:  nxtsec_reset
!
! Description: <TBS>
!
!*******************************************************************************

subroutine nxtsec_reset()

  implicit none

  iblock = 0

  return

end subroutine nxtsec_reset

end module nextprmtop_section_mod
