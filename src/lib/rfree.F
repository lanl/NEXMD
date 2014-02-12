c-----------------------------------------------------------------------
      subroutine rfree(ifld,ihol,ivar,fvar,in,iout)
C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
C
      implicit double precision (a-h,o-z)
c
c     Author:  George Seibel
c
c     This is a free format reader for mixed Hollerith and numeric data.
c     This code is a complete re-write of the old Amber rfree, and is
c     now machine independent ANSI fortran 77.
c     Rev 02-May-89: changed return on EOF back to stop. (Edit no longer
c                    needs this bogus feature.)
c     Rev 14-Mar-89: add else to elseif check on ifld()
c     Rev 01-Mar-89: initialize ierr to 0
c     Rev 22-Feb-89: fixed bug in ifld() interpretation
c     Rev 20-Feb-89: changed stop on EOF to return
c     Rev 20-Jan-92: made ch2int() ebcdic-capable (BR)
c
c     PARAMETERS:
c
      integer LENBUF
c        ... length of character buffer for input line
#ifdef IBM_VM_CMS
c     vm/cms/mvs want fixed records, so we give the obvious length
      parameter (LENBUF=80)
#else
      parameter (LENBUF=132)
#endif
c
c     INPUT:
c
      integer ifld(*)
c        ... code for field types to be read:  1 = Hollerith (4 byte int)
c            2 = integer   3 = float  0 = no more fields
      integer in, iout
c        ... input and output logical unit numbers
c
c     OUTPUT:
c
      character(len=4) ihol(*)
c        ... extracted Hollerith data (4 byte max for Amber)
      integer ivar(*)
c        ... extracted integer data
      dimension fvar(*)
c        ... extracted floating pt data
c
c     INTERNAL:
c
      character*(LENBUF) buf, token
c        ... input line buffer,  temp for undecoded tokens
      character*4 blank
c        ... just 4 bytes of space
      integer nvar
c        ... number of variables to read
      integer ntoken, nint, nhol, nflt
c        ... counters for tokens, int, hol, and real variables
      logical inword
c        ... true if tokenizer is in a word
      integer ipt, i
c        ... pointer into char buffer, loop index
      integer ibeg, iend
c        ... buf indices = beginning and end of current token
      integer ival, ierr
c        ... temp for integer values, error return for ch2int()
c
c
      ierr = 0
      nvar = 0
      ntoken = 0
      nint = 0
      nhol = 0
      nflt = 0
      blank = ' '
      ibeg = 1
      iend = LENBUF
      inword = .false.
c
c     --- initialize the output arrays ---
c
      do 100 i = 1, 80
           if (ifld(i) .le. 0) go to 110
           nvar = nvar+1
           if (ifld(i) .eq. 1) then
                nhol = nhol + 1
                read(blank,'(a4)') ihol(nhol)
           elseif (ifld(i) .eq. 2) then
                nint = nint + 1
                ivar(nint) = 0
           elseif (ifld(i) .eq. 3) then
                nflt = nflt + 1
                fvar(nflt) = 0.0d0
           else
                write(iout,'(5x,a)') 'rfree: bogus ifld()'
                call mexit(iout, 1)
           endif
  100 continue
c
  110 continue
c
c     --- read entire line into character buffer ---
c
      read(in,'(a)',end=1000) buf
c
c     --- tokenize buf using any whitespace as delimitter ---
c
      nint = 0
      nhol = 0
      nflt = 0
      do 200 ipt = 1, LENBUF
           if (ntoken .ge. nvar) return
           if (.not. inword) then
c               --- look for start of word = non-whitespace --
                if (buf(ipt:ipt) .gt. ' ') then
                     inword = .true.
                     ibeg = ipt
                endif
           else
c               --- look for end of word = whitespace or end of buf ---
                if (buf(ipt:ipt) .le. ' ' .or. ipt .ge. LENBUF) then
                     inword = .false.
                     ntoken = ntoken + 1
                     iend = ipt
                     token = buf(ibeg:iend)
                     lenstr = iend - ibeg
c
c                    --- decode according to ifld() ---
c
                     if (ifld(ntoken) .eq. 1) then
c                         --- Hollerith was a Great Man ---
                          nhol = nhol + 1
                          read(token,'(a4)') ihol(nhol)
                     elseif (ifld(ntoken) .eq. 2) then
c                         --- Integer Field ---
                          nint = nint + 1
                          call ch2int(lenstr,token,ival,ierr)
                          if (ierr .ne. 0) go to 900
                          ivar(nint) = ival
                     elseif (ifld(ntoken) .eq. 3) then
c                         --- Floating Point field ---
                          nflt = nflt + 1
                          if (index(token,'.') .gt. 0) then
c                              --- if decimal pt, use internal read ---
                               read(token,'(f20.0)',err=900) fvar(nflt)
                          else
c                              --- no decimal, use char to int routine ---
                               call ch2int(lenstr,token,ival,ierr)
                               if (ierr .ne. 0) go to 900
                               fvar(nflt) = float(ival)
                          endif
                     endif
                endif
           endif
  200 continue
      return
  900 continue
c
c     --- token could not be decoded ---
c
      write(iout,'(/5x,a,i3,i3,a,/,a)')'rfree: Error decoding variable',
     +      ntoken, ifld(ntoken), ' from:', buf(1:iend)
      write(iout,'(/5x,a)')'this indicates that your input contains',
     + ' incorrect information'
      write(iout,'(/5x,a,i3,a)') 'field ',ntoken,
     + ' was supposed to have',
     + ' a (1=character, 2=integer, 3=decimal) value'

      call mexit(iout, 1)
 1000 continue
c
c     --- hit EOF ---
c
      write(iout,'(/5x,a,i3)') 'rfree: End of file on unit ', in
      call mexit(iout, 1)
      end
c-----------------------------------------------------------------------
      subroutine ch2int(lenstr, string,ival,ierr)
C
C************************************************************************
C                              AMBER                                   **
C                                                                      **
C                  Copyright (c) 1986, 1991, 1995                      **
C             Regents of the University of California                  **
C                       All Rights Reserved.                           ** 
C                                                                      **
C  This software provided pursuant to a license agreement containing   **
C  restrictions on its disclosure, duplication, and use. This software **
C  contains confidential and proprietary information, and may not be   **
C  extracted or distributed, in whole or in part, for any purpose      **
C  whatsoever, without the express written permission of the authors.  **
C  This notice, and the associated author list, must be attached to    **
C  all copies, or extracts, of this software. Any additional           **
C  restrictions set forth in the license agreement also apply to this  **
C  software.                                                           **
C************************************************************************
C
c     converts character representations of cardinal numbers to integer
c     Author:  George Seibel
c     Rev 01-Mar-89: initialize ierr to 0
c     Initial Rev: 19-Dec-87
*     implicit none
c
c     INPUT:
c
      integer lenstr
c        ... length of string 
      character string*(*)
c        ... character representation of legitimate integer
c
c     OUTPUT:
c
      integer ival
c        ... the integer result
      integer ierr
c        ... returned as one if any decoding error, zero if ok
c
c     INTERNAL:
c
      integer i, j, num, ifirst, last, idec
      logical isneg
c
      ierr = 0
c
c     --- look for minus sign ---
c
      isneg = (index(string,'-') .gt. 0)
c
c     --- find first and last numeric character ---
c
      last = lenstr
      do 200 i = 1, lenstr
           if (string(i:i).ge.'0'.and. string(i:i).le.'9') then
                ifirst = i
                do 100 j = ifirst, lenstr
                     if (string(i:i).lt.'0'.or.string(i:i).gt.'9') then
                          last = j - 1
                          go to 300
                     endif
  100           continue
                go to 300
           endif
  200 continue
c
c     --- no numerics found - error return ---
c
      ierr = 1
      return
c
c     --- crunch the number ---
c
  300 continue
      num = 0
      idec = 0
      do 400 i = last, ifirst, -1
           num = num + (ichar(string(i:i)) - ichar('0')) * 10**idec
           idec = idec + 1
  400 continue
      if (isneg) num = -num
      ival = num
      return
      end
c-----------------------------------------------------------------------
