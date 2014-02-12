#include "copyright.i"
!*******************************************************************************
!
! Module: file_io_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module file_io_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  amopen
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine amopen(lun, fname, fstat, ffmt, facc)

  use pmemd_lib_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: lun     ! logical unit number
  character(*)          :: fname   ! file name 
  character(1)          :: fstat   ! status code:
                                   ! "N", "O", or "U" = new, old, unk.
  character(1)          :: ffmt    ! format code: "U", "F" = unform., form.
  character(1)          :: facc    ! access code:
                                   ! "R", "W", "A" = read, read/write, append
! Local variables:

  integer               :: ios     ! i/o status variable
  character(11)         :: kform   ! form keyword
  character(7)          :: stat    ! status keyword
 
  if (fstat .eq. 'N') then
    stat = 'NEW'
  else if (fstat .eq. 'O') then
    stat = 'OLD'
  else if (fstat .eq. 'U') then
    stat = 'UNKNOWN'
  else
    write(mdout, '(/,2x,a,i4)') 'amopen: bogus fstat, unit ', lun
    call mexit(6, 1)
  end if

  if (ffmt .eq. 'U') then
    kform = 'UNFORMATTED'
  else if (ffmt .eq. 'F') then
    kform = 'FORMATTED'
  else
    write(mdout, '(/,2x,a,i4)') 'amopen: bogus ffmt, unit', lun
    call mexit(6, 1)
  end if

  open(lun, file=fname, status=stat, form=kform, iostat=ios)

  if (ios .ne. 0) then

    if (lun .eq. 6) then
#ifndef DUMB
      write(0, '(/,2x,a,i4,a,a)') 'Unit ', lun, ' Error on OPEN: ', fname
#endif
      close(0)
    else
      write(mdout, '(/,2x,a,i4,a,a)') 'Unit ', lun, ' Error on OPEN: ', fname
      close(unit = 6)
#ifndef DUMB
      write(0, '(/,2x,a,i4,a,a)') 'Unit ', lun, ' Error on OPEN: ', fname
#endif
      close(0)
    end if
    
    call mexit(6, 1)

  end if

  rewind(lun)

  return

end subroutine amopen

!*******************************************************************************
!
! Subroutine:  echoin (ECHO INput)
!
! Description:
!              
! This routine is called by the mdread routine for an Amber suite program.
! It reads and echos the information in the input file to the user.
!
! Special provision is made for the lines written to the file by the
! Amber/Interface to log the Interface lines used to generate the input
! file.
!
! Before this routine returns, the file on unit IN will be rewound.
!
! Author: David A. Pearlman
! Date: 8/92
!
! in:           The unit the input file is read from.
! iout:         The unit to which the echoed information is to be written.
!
!*******************************************************************************

subroutine echoin(in,iout)
   
   implicit none

! Formal arguments:

  integer               :: in
  integer               :: iout

! Local variables:

  character(80)         :: aa
  character(14)         :: bgflg
  character(14)         :: enflg
  integer               :: i
  integer               :: inter
  integer               :: intlin

  data bgflg /'!! BEGIN Amber'/
  data enflg /'!! END   Amber'/

! First echo the Interface script, if any:

      inter = 0
      intlin = 0
      do 10 i = 1, 999999
        read(in, 500, end = 20) aa
        if (aa(1:14) .eq. enflg) inter = 0
        if (inter .eq. 1) write(iout, 500) aa(1:79)
        if (aa(1:14) .eq. bgflg) then
          inter = 1
          if (intlin .eq. 0) then
            write(iout, *)
            write(iout, 1000)
            write(iout, *)
          end if
          intlin = 1
        end if
   10 continue

! Now echo the standard Amber input lines:

   20 if (intlin .gt. 0) write(iout, 501)
      rewind(in)
      write(iout, *)
      write(iout, 1010)
      write(iout, *)
      inter = 0
      do 30 i = 1, 999999
        read(in, 500, end = 40) aa
        if (aa(1:14) .eq. bgflg) inter = 1
        if (inter .eq. 0) write(iout, 500) aa(1:79)
        if (aa(1:14) .eq. enflg) inter = 0
   30 continue

   40 rewind(in)

      return

   ! Format statements
   
   500 format(a)
   501 format(79('-'))
   1000 format(' The Interface script used to generate the input file:')
   1010 format(' Here is the input file:')

end subroutine echoin 

!*******************************************************************************
!
! Subroutine:   nmlsrc (NaMeList SeaRCh)
!
! Description:
!
! This routine searches the file on unit iun for the namelist with name
! name. The namelist is denoted by a line where the first non-blank character
! string is either &NAME or $NAME.
!
! ifind: 0 Namelist was not found. File will be rewound before return.
!        1 Namelist was found. File will be backspaced to the line that
!          contained the namelist flag before returning.
!
! Limitations: the namelist flat must fit entirely in the first 80 columns
!              of the line on which it appears.
!
! Author: David A. Pearlman
! Date: 1/91
!
!*******************************************************************************

subroutine nmlsrc(name, iun, ifind)

  implicit none

! Formal arguments:

      character(*)      name
      integer           iun
      integer           ifind

! Local variables:

      character(80)     aline
      integer           i
      integer           ilen
      integer           j

      ilen = len(name)

      do 10 i = 1, 9999999
          read(iun, 1000, end = 500, err = 500) aline
          do 20 j = 1, 80
             if (aline(j:j) .ne. ' ') then
                if (aline(j:j) .eq. '&' .or. aline(j:j) .eq. '$') then
                   if (80 - j + 1 .ge. ilen) then
                      if (aline(j + 1:j + ilen) .eq. name) goto 200
                   else
                      goto 10
                   end if
                else
                   goto 10
                end if
             end if
   20     continue
   10 continue

! Namelist flag found:

  200 ifind = 1
      backspace(iun)
      return

! Namelist flag not found:

  500 ifind = 0
      rewind(iun)
      return

! Format statements:

 1000 format(a)

end subroutine nmlsrc

!*******************************************************************************
!
! Subroutine:  rfree
!
! Description:
!              
! Author:  George Seibel
!
! This is a free format reader for mixed Hollerith and numeric data.
! This code is a complete re-write of the old Amber rfree, and is
! now machine independent ANSI fortran 77.
! Rev 02-May-89: changed return on EOF back to stop. (Edit no longer
!                needs this bogus feature.)
! Rev 14-Mar-89: add else to else if check on ifld()
! Rev 01-Mar-89: initialize ierr to 0
! Rev 22-Feb-89: fixed bug in ifld() interpretation
! Rev 20-Feb-89: changed stop on EOF to return
! Rev 20-Jan-92: made ch2int() ebcdic-capable (BR)
!
!*******************************************************************************

subroutine rfree(ifld, ihol, ivar, fvar, in, iout)

  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

! (INPUT)

! code for field types to be read:  1 = Hollerith (4 byte int), 2 = integer,
!                                   3 = float,  0 = no more fields.

      integer           ifld(*)

! input and output logical unit numbers:

      integer           in, iout

! (OUTPUT)

! extracted Hollerith data (4 byte max for Amber):

      character(4)      ihol(*)

! extracted integer data:

      integer           ivar(*)

! extracted floating pt data:

      double precision  fvar(*)

! Local variables:

! length of character buffer for input line:

      integer           lenbuf

      parameter (lenbuf = 132)

! input line buffer, temp for undecoded tokens:

      character(lenbuf) buf, token

! just 4 bytes of space:

      character(4)      blank

! number of variables to read:

      integer           nvar

! counters for tokens, int, hol, and real variables:

      integer           ntoken, nint, nhol, nflt

! true if tokenizer is in a word:

      logical           inword

! pointer into char buffer, loop index:

      integer           ipt, i

! buf indices = beginning and end of current token:

      integer           ibeg, iend

! temp for integer values, error return for ch2int():

      integer           ival, err_code

      integer           lenstr

      err_code = 0
      nvar = 0
      ntoken = 0
      nint = 0
      nhol = 0
      nflt = 0
      blank = '    '
      ibeg = 1
      iend = lenbuf
      inword = .false.

! Initialize the output arrays:

      do 100 i = 1, 80
           if (ifld(i) .le. 0) goto 110
           nvar = nvar + 1
           if (ifld(i) .eq. 1) then
                nhol = nhol + 1
                read(blank, '(a4)') ihol(nhol)
           else if (ifld(i) .eq. 2) then
                nint = nint + 1
                ivar(nint) = 0
           else if (ifld(i) .eq. 3) then
                nflt = nflt + 1
                fvar(nflt) = 0.0d0
           else
                write(iout, '(5x,a)') 'rfree: bogus ifld()'
                call mexit(iout, 1)
           end if
  100 continue

  110 continue

! Read entire line into character buffer: 

      read(in, '(a)', end = 1000) buf

! Tokenize buf using any whitespace as delimiter:

      nint = 0
      nhol = 0
      nflt = 0

      do 200 ipt = 1, lenbuf
           if (ntoken .ge. nvar) return
           if (.not. inword) then

! Look for start of word = non-whitespace:

                if (buf(ipt:ipt) .gt. ' ') then
                     inword = .true.
                     ibeg = ipt
                end if
           else

! Look for end of word = whitespace or end of buf:

                if (buf(ipt:ipt) .le. ' ' .or. ipt .ge. lenbuf) then
                     inword = .false.
                     ntoken = ntoken + 1
                     iend = ipt
                     token = buf(ibeg:iend)
                     lenstr = iend - ibeg

! Decode according to ifld():

                     if (ifld(ntoken) .eq. 1) then

! Hollerith:

                          nhol = nhol + 1
                          read(token, '(a4)') ihol(nhol)
                     else if (ifld(ntoken) .eq. 2) then

! Integer:

                          nint = nint + 1
                          call ch2int(lenstr, token, ival, err_code)
                          if (err_code .ne. 0) goto 900
                          ivar(nint) = ival
                     else if (ifld(ntoken) .eq. 3) then

! Floating Point:

                          nflt = nflt + 1
                          if (index(token, '.') .gt. 0) then

! If decimal pt, use internal read:

                               read(token, '(f20.0)', err = 900) fvar(nflt)
                          else

! No decimal, use char to int routine:

                               call ch2int(lenstr, token, ival, err_code)
                               if (err_code .ne. 0) goto 900
                               fvar(nflt) = dble(ival)
                          end if
                     end if
                end if
           end if
  200 continue
      return
  900 continue

! Token could not be decoded:

      write(iout, '(/5x,a,i3,i3,a,/,a)')'rfree: Error decoding variable', &
            ntoken, ifld(ntoken), ' from:', buf(1:iend)
      call mexit(iout, 1)
 1000 continue

! Hit EOF:

      write(iout, '(/5x,a,i3)') 'rfree: End of file on unit ', in
      call mexit(iout, 1)

contains

!*******************************************************************************
!
! Subroutine:  ch2int
!
! Description:
!              
! Converts character representations of cardinal numbers to integer.
! Author:  George Seibel
! Rev 01-Mar-89: initialize err_code to 0
! Initial Rev: 19-Dec-87
!
!*******************************************************************************

subroutine ch2int(lenstr, string, ival, err_code)

  implicit none

! Formal arguments:

! (INPUT)

! length of string:

      integer           lenstr

! character representation of legitimate integer:

      character         string*(*)

! (OUTPUT)

! the integer result:

      integer           ival

! returned as one if any decoding error, zero if ok:

      integer           err_code

! Local variables:

      integer           i, j, num, ifirst, last, idec
      logical           isneg

      err_code = 0

! Look for minus sign:

      isneg = (index(string, '-') .gt. 0)

! Find first and last numeric character:

      last = lenstr
      do 200 i = 1, lenstr
           if (string(i:i) .ge. '0'.and. string(i:i) .le. '9') then
                ifirst = i
                do 100 j = ifirst, lenstr
                     if (string(i:i) .lt. '0' .or. string(i:i) .gt. '9') then
                          last = j - 1
                          goto 300
                     end if
  100           continue
                goto 300
           end if
  200 continue

! No numerics found - error return:

      err_code = 1
      return

! Crunch the number:

  300 continue
      num = 0
      idec = 0
      do 400 i = last, ifirst, - 1
           num = num + (ichar(string(i:i)) - ichar('0')) * 10**idec
           idec = idec + 1
  400 continue
      if (isneg) num = - num
      ival = num

      return

end subroutine ch2int

end subroutine rfree

!*******************************************************************************
!
! Subroutine:  rgroup
!
! Description:
!             
! Read in groups each atom is in.
!
! Will keep reading in groups of cards until a blank card is read.  All atoms
! in this card group will be put in the same group.  the only exception is a
! res card which starts with a - res no. The residues in this group will all be
! put in indiv.  groups.  Any other groups in this section must also start
! with a  i.e. they must also be indiv. groups.  Any time a "find" card is read,
! only those atoms that meet the specifications on at least 1 of the find
! cards will be included in the group.  This option may be ended by reading in
! "NFIND", terminating the input for the groups, or by reading in another
! "FIND".
!
! Mods for Rev A by gls:
!   - cpp selectable precision
!   - wildcard functionality fixed. (see findp)
!   - format change to print less garbage
!   - unreferenced statement labels removed
!   - changed title write from 12 to 19 A4. (not 20 to avoid wraps)
!   - changed to read from unit nf -- dac 12/90
!
!*******************************************************************************

subroutine rgroup(natom, natc, nres, ipres, lbres, igraph, &
                  isymbl, itree, igroup, weit, konst, belly, nfu)

  use pmemd_lib_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

      integer           natom
      integer           natc
      integer           nres
      integer           ipres(*)
      character(4)      lbres(*)
      character(4)      igraph(*)
      character(4)      isymbl(*)
      character(4)      itree(*)
      integer           igroup(*)
      double precision  weit(*)
      logical           konst
      logical           belly
      integer           nfu

! Local variables:

      double precision  fvar(20)
      integer           i
      integer           i1
      integer           i2
      integer           iend
      integer           iii
      character(4)      ifind
      integer           ifld(20)
      integer           igrp(8)
      character(4)      ihol(20)
      character(4)      iiend
      integer           isrch
      integer           istart
      integer           itime
      integer           ivar(20)
      integer           izero
      integer           j
      integer           jfld(20)
      character(4)      jgraph(11)
      integer           jgrp(8)
      character(4)      jresnm(11)
      character(4)      jsymbl(11)
      character(4)      jtree(11)
      integer           k
      character(4)      katn
      character(4)      ksear
      character(4)      ktypg
      integer           last_namelist_line_ctr
      integer           lfind
      integer           line_ctr
      integer           lsign
      logical           misc
      integer           natmg
      integer           nf
      character(4)      nfind
      integer           ngrp
      character(1)      readbuf(80)     ! Buffer for finding start of group.
      character(4)      title(20)
      double precision  wt

      common /propf/ jgraph, jresnm, jsymbl, jtree, isrch

      data ifld / 2*0, 13*2 , 5*0 /
      data jfld / 4*1, 16*0 /
      data katn  / 'ATOM' /
      data ksear / 'SEAR' /
      data ifind / 'FIND' /
      data nfind / 'NFIN' /
      data iiend / 'END ' /

! Res card lists residue groups igrp(i) to jgrp(i):
!
! If 1st value is negative, then every residue from igrp(i) to jgrp(i) is
! considered a separate group.
! If 2nd value = 0, then subgroup consists of 1 residue.
!
! Atom card reads in groups of atoms:
!
! If 2nd atom in each pair = 0 , then subgroup consists of 1 atom res and
! atom cards may be read in any order.
! End input with an "END " card.

! UPDATE NOTE:  This routine has been essentially broke in that it did not work
! unless you were positioned in the file at the beginning of the title for the
! group.  With the introduction of ewald and debugf namelists, as well as with 
! the use of wt namelists, there are lots of ways for this to not work out. We
! therefore introduce code that specifically places us after the last line that
! begins with 0 or more blanks followed by either a '&', '$', or '/'.  These
! characters should cover all conventions for namelist termination, and should
! allow us to start reading at the line following the last namelist.  If the
! title line for the group begins with one of these characters, things will not
! work out! - RED

  rewind nfu

  line_ctr = 0
  last_namelist_line_ctr = 0

  do
    read(nfu, end = 1000, err = 1010, fmt = '(80a1)') readbuf
    line_ctr = line_ctr + 1

    do i = 1, 80
      if (readbuf(i) .ne. ' ') then
        if (readbuf(i) .eq. '&' .or. &
            readbuf(i) .eq. '/' .or. &
            readbuf(i) .eq. '$') last_namelist_line_ctr = line_ctr
        exit
      end if
    end do
  end do

1000 continue   ! Completed the first read.  Reposition to correct spot:

  rewind nfu

  ! And read forward to the group title:

  do i = 1, last_namelist_line_ctr
    read(nfu, end = 1010, err = 1010, fmt = '(80a1)') readbuf
  end do

! Zero ngrp before calling this routine.  This is the number of groups read.

      ngrp = 0
      nf = nfu
      if (nf .le. 0) nf = 5

! Initialise the group array: 

      do 100 i = 1, natom
        igroup(i) = 0
        if (konst) weit(i) = 0.0d0
  100 continue

   22 continue

        itime = 0
        lsign = 0
        izero = 0
        isrch = 0
        natmg = 0

! Read different groups:

        read(nf, 9208) (title(k), k = 1, 19)

        if (title(1) .eq. iiend) then
          write(mdout, '(4x,a)') '----- END OF GROUP READ -----'
          goto 900
        end if

        ngrp = ngrp + 1
        write(mdout, '(4x,a,i5,a)') '----- READING GROUP ', &
                             ngrp, '; TITLE:'
        write(mdout, 9218) (title(k), k = 1, 19)

! If constrained groups read the weight for each group: 

        if (konst) then
          ifld(1) = 3
          ifld(2) = 0
          call rfree(ifld, ihol, ivar, fvar, nf, 6)
          wt = fvar(1)
          write(mdout, 9018) ngrp, wt
        end if
   10   continue

! Read the group cards: 

        ifld(1) = 1
        ifld(2) = 2
        call rfree(ifld, ihol, ivar, fvar, nf, 6)
        ktypg = ihol(1)
        k = 1
        do 120 i = 1, 7
          igrp(i) = ivar(k)
          jgrp(i) = ivar(k + 1)
          k = k + 2
  120   continue
        if (ktypg .eq. iiend) goto 16
        if (ktypg .eq. nfind) then
          write(mdout, 199)
          isrch = 0
          goto 10
        end if
        if (ktypg .eq. ifind) then

! Find option ... Read the atom specifics: 

          write(mdout, 200)
          do 64 iii = 1, 12
            call rfree(jfld, ihol, ivar, fvar, nf, 6)
            jgraph(iii) = ihol(1)
            jsymbl(iii) = ihol(2)
            jtree(iii) = ihol(3)
            jresnm(iii) = ihol(4)
            if (jgraph(iii) .eq. ksear) goto 65
            write(mdout, 202) jgraph(iii), jsymbl(iii), jtree(iii), jresnm(iii)
   64     continue
   65     continue

          isrch = iii - 1
          if (isrch .gt. 10) write(mdout, 66) isrch

! Now read in res and atoms to be searched: 

          goto 10
        end if
        itime = itime+1
        if (ktypg .ne. katn) then

! 1st group of 1st card must be - if any - numbers are found:

          if (itime .eq. 1 .and. igrp(1) .lt. 0) lsign = 1
          do 12 i = 1, 7
            i1 = igrp(i)
            if (i1 .eq. 0) goto 10
            i2 = jgrp(i)
            if (i2 .gt. nres) i2 = nres
            if (i1 .gt. 0 .and. lsign .eq. 1) goto 36
            if (i1 .lt. 0 .and. lsign .eq. 0) goto 36
            if (i1 .lt. 0) i1 = - i1
            if (i2 .le. 0) i2 = i1
            if (lsign .eq. 0) write(mdout, 14) ngrp, i1, i2
            do 13 j = i1, i2
              istart = ipres(j)
              iend = ipres(j + 1) - 1
              do 45 k = istart, iend
                if (isrch .gt. 0)  &
                  call findp(k, j, lfind, nres, ipres, lbres, isymbl, &
                             itree, igraph)
                if (isrch .gt. 0 .and. lfind .eq. 0) goto 45
                igroup(k) = ngrp
                if (konst) weit(k) = wt
                natmg = natmg + 1
   45         continue
              if (lsign .eq. 1) write(mdout, 46) ngrp, j
              if (lsign .eq. 1) ngrp = ngrp + 1
   13       continue
   12     continue
          goto 10
        end if

! Atom type constraints:

        if (lsign .eq. 1) goto 36
        write(mdout, 51) ngrp
        do 17 i = 1, 7
          i1 = igrp(i)
          if (i1 .lt. 0) goto 36
          if (i1 .eq. 0) goto 10
          i2 = jgrp(i)
          if (i2 .gt. natom) i2 = natom
          if (i2 .le. 0) i2 = i1
          write(mdout, 52) i1, i2
          do 18 j = i1, i2
            if (isrch .gt. 0)  &
              call findp(j, izero, lfind, nres, ipres, lbres, isymbl, &
                         itree, igraph)
            if (isrch .gt. 0 .and. lfind .eq. 0) goto 18
            natmg = natmg + 1
            igroup(j) = ngrp
            if (konst) weit(j) = wt
   18     continue
   17   continue
        goto 10

   16   if (lsign .eq. 1) ngrp = ngrp - 1
        if (itime .eq. 0) ngrp = ngrp - 1
        write(mdout, 222) natmg
      goto 22

   36 continue
      write(mdout, 127) ktypg, (igrp(i), jgrp(i), i = 1, 7)
      goto 10

! All groups are read;  return.

  900 continue

      if (konst) then

! Gather all the constrained atoms together:

        natc = 0
        do 920 i = 1, natom
          if (igroup(i) .le. 0) goto 920
          natc = natc + 1
          igroup(natc) = i
          weit(natc) = weit(i)
  920   continue

! Do not print the history of constrained atoms:

      else if (.not.belly) then

! Put the atoms which are not in the defined groups as the last group:

        misc = .false.
        do 820 i = 1, natom
          if (igroup(i) .ne. 0) goto 820
          misc = .true.
          igroup(i) = ngrp + 1
  820   continue
        if (misc) ngrp = ngrp + 1
      end if

      return

1010 continue

  write(mdout, *) '     ERROR READING GROUP INPUT!'
  call mexit(6, 1)

  199 format(' ', 5x, 'END OF ATOM SPECIFICATION', /)
  200 format(' ', 5x, 'ALL ATOMS THAT MEET 1 OF THE FOLLOWING', &
          ' SPECIFICATIONS WILL BE INCLUDED IN GROUP BELOW', /)
  202 format(' ', 5x, 'GRAPH NAME  = ', a4, 2x, 'SYMBOL  = ', a2, 4x, &
              'TREE SYMBOL  = ', a1, 5x, 'RESIDUE TYPE  = ', a4, /)
   66 format(' ', 5x, '**** NUMBER OF FIND CARDS  = ', i5, 2x, &
              'IS TOO BIG ******', /)
   14 format(' GRP', i5, ' RES', i5, ' TO ', i5)
   46 format(' ', 5x, 'GROUP', i5, 3x, 'CONSISTS OF RESIDUE', i5, /)
   51 format(' ', 5x, 'GROUP', i5, 3x, 'CONSISTS OF ATOMS -', /)
   52 format(' ', 34x, i5, 2x, 'TO', i5)
  222 format(' ', 5x, 'Number of atoms in this group  = ', i5)
  127 format(' ', 5x, '***PROBLEMS WITH GROUP', a4, 14i5, '*******', /)
 9018 format(/5x, 'GROUP', i5, ' HAS HARMONIC CONSTRAINTS', f12.5)
 9208 format(20a4)
 9218 format(1x, 20a4)

contains

!*******************************************************************************
!
! Internal Subroutine:  findp
!
! Description:
!             
! Checks if a given atom has certain characteristics.
!
! Rev A mod (G. Seibel, Apr 89)
! Changed iblank, jblank to iwild, jwild, to give the wildcard
! functionality promised by the doc.
! isymbl() dimensioned. (this was long-standing bug)
!
!*******************************************************************************

subroutine findp(iatom, ires, iloc, nres, ipres, lbres, isymbl, itree, igraph)

  implicit none

! Formal arguments:

      integer           iatom
      integer           ires
      integer           iloc
      integer           nres
      integer           ipres(*)
      character(4)      lbres(*)
      character(4)      isymbl(*)
      character(4)      itree(*)
      character(4)      igraph(*)

! Local variables:

      integer           i
      integer           isrch
      character(4)      iwild
      character(4)      jgraph(11)
      character(4)      jresnm(11)
      character(4)      jsymbl(11)
      character(4)      jtree(11)
      character(4)      jwild

      common /propf/ jgraph, jresnm, jsymbl, jtree, isrch

      data iwild /'*   '/
      data jwild /'* '/

! BUGBUG - The 16 bit jwild above is probably at the least unnecessary.

      iloc = 0
      if (ires .eq. 0) call findrs(iatom, ires, nres, ipres)

      do 10 i = 1, isrch
        if ((jresnm(i) .ne. iwild) .and. (jresnm(i) .ne. lbres(ires))) &
                      goto 10
        if ((jgraph(i) .ne. iwild) .and. (jgraph(i) .ne. igraph(iatom))) &
                      goto 10
        if ((jtree(i) .ne. jwild) .and. (jtree(i) .ne. itree(iatom))) &
                      goto 10
        if ((jsymbl(i) .ne. jwild) .and. (jsymbl(i) .ne. isymbl(iatom))) &
                      goto 10
        iloc = 1

        goto 20

   10 continue
   20 continue

      return

end subroutine findp

!*******************************************************************************
!
! Internal Subroutine:  findrs
!
! Description: <TBS>
!             
!*******************************************************************************

subroutine findrs(numa, ires, nres, ipres)

  use pmemd_lib_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

      integer           numa
      integer           ires
      integer           nres
      integer           ipres(*)

! Local variables:

      integer           i
      integer           im

      if (numa .le. 0) goto 11
      if (numa .ge. ipres(nres)) ires = nres
      if (numa .ge. ipres(nres)) return
      im = nres - 1

      do 10 i = 1, im
        if (numa .ge. ipres(i) .and. numa .lt. ipres(i + 1)) goto 12
   10 continue

   11 write(mdout, 100) numa
      goto 200
   12 ires = i
      return
  200 continue
  100 format(/2x, 'PROBLEMS FINDING RESIDUE OF ATOM ', i5)
      call mexit(6, 1)

end subroutine findrs

end subroutine rgroup

#ifdef AMOEBA
!*******************************************************************************!
! Subroutine:  get_crdfile_type
!
! Description: <TBS>
!
!*******************************************************************************

subroutine get_crdfile_type(crdfile_name, crdfile, crdfile_type)

  use file_io_dat_mod
  use nextinpcrd_section_mod

  implicit none

! Formal arguments:

  character(max_fn_len), intent(in)     :: crdfile_name
  integer, intent (in)                  :: crdfile
  integer, intent (out)                 :: crdfile_type

! Local variables:

  integer                       :: iok, n, ionerr
  character(len = 80)           :: fmt, fmtin, dtype

! BUGBUG - We should probably be basing our decision as to whether this is an
!          old or new style inpcrd file on looking for a version tag instead of
!          looking for a tagged title...

  call amopen(crdfile, crdfile_name, 'O', 'F', 'R')
  call nxtsec_crd_reset()

  fmtin = '(a)'
  dtype = 'TITLE'
  ionerr = 1 ! not fatal if missing
  call nxtsec_crd(crdfile, mdout, ionerr, fmtin, dtype, fmt, iok)

  if (iok .eq. 0) then 
    crdfile_type = crdfile_type_new
  else
    crdfile_type = crdfile_type_old  ! old style crdfile, no section hdrs.
  end if

  call nxtsec_crd_reset()
  close(crdfile)

  return 

end subroutine get_crdfile_type
#endif /* AMOEBA */

end module file_io_mod
