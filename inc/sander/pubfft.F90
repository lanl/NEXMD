#include "copyright.h"
#include "dprec.fh"


!     -----------------------------------------------------------------
!     All of the ewald code and supporting routines were written and
!     contributed by Tom Darden from the National Institute of
!     Environmental Health Sciences division of the NIH.
!     Originally written with a modified version of AMBER 3A, the code
!     was updated during the summer of 1994 to be compatible with
!     AMBER 4.1.

!     The 1D FFT code is a double precision version of fftpack from
!     netlib, written by Paul N. Swartztrauber at NCAR Boulder Colorado.
!     This file once contained code necessary to perform 3D FFTs where
!     libraries were not available.  It was based on piecing together a
!     series of 1D FFTs and was probably not super efficient.
!     In July 2007 it was unused and was removed by Scott Brozell.

!     See www.netlib.org/fftpack

!     -----------------------------------------------------------------

!     Revisions: tec3 (added comments, modified source to fit it with
!     default Makefile, Compile, MACHINE source handling scheme, added
!     CPP selectable precision, added CPP control of SGI multiprocessing,
!     changed output formats to rid write(6,*), streamlined source,
!     merged routines, etc...)

!     NOTE: Some of the routine names within here are rather long so
!     may choke lame compilers.  This code was developed on SGI
!     computers so should surely work on these; it has also been
!     limitedly tested on Cray and HP machines.

!     The following C preprocessor code (only visible in the
!     untransformed source) is a hack to allow single precision
!     versions of an originally all double precision code.
!     It is recommended however that users run with double
!     precision... (by default sander is double precision)


!     The following routines are defined, in alphabetical order:
!     CFFTB
!     CFFTB1
!     CFFTF
!     CFFTF1
!     CFFTI
!     CFFTI1
!     PASSB
!     PASSB2
!     PASSB3
!     PASSB4
!     PASSB5
!     PASSF
!     PASSF2
!     PASSF3
!     PASSF4
!     PASSF5


!     --- CFFTB ---

! Backward complex Transform


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cfftb here]
subroutine cfftb (n,c,wsave)
   implicit none
   integer:: iw1, iw2, n
   _REAL_ :: c, wsave
   dimension       c(*)       ,wsave(*)
   if (n == 1) return
   iw1 = n+n+1
   iw2 = iw1+n+n
   call cfftb1 (n,c,wsave,wsave(iw1),wsave(iw2))
   return
end subroutine cfftb 

!     --- CFFTB1 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cfftb1 here]
subroutine cfftb1 (n,c,ch,wa,ifac)
   implicit none
   integer:: i, idl1, ido, idot, ifac, ip, iw, ix2, ix3, ix4, k1, &
        l1, l2, n, n2, na, nac, nf
   _REAL_ :: c, ch, wa
   dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
   nf = ifac(2)
   na = 0
   l1 = 1
   iw = 1
   do 116 k1=1,nf
      ip = ifac(k1+2)
      l2 = ip*l1
      ido = n/l2
      idot = ido+ido
      idl1 = idot*l1
      if (ip /= 4) goto 103
      ix2 = iw+idot
      ix3 = ix2+idot
      if (na /= 0) goto 101
      call passb4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
      goto 102
      101 call passb4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
      102 na = 1-na
      goto 115
      103 if (ip /= 2) goto 106
      if (na /= 0) goto 104
      call passb2 (idot,l1,c,ch,wa(iw))
      goto 105
      104 call passb2 (idot,l1,ch,c,wa(iw))
      105 na = 1-na
      goto 115
      106 if (ip /= 3) goto 109
      ix2 = iw+idot
      if (na /= 0) goto 107
      call passb3 (idot,l1,c,ch,wa(iw),wa(ix2))
      goto 108
      107 call passb3 (idot,l1,ch,c,wa(iw),wa(ix2))
      108 na = 1-na
      goto 115
      109 if (ip /= 5) goto 112
      ix2 = iw+idot
      ix3 = ix2+idot
      ix4 = ix3+idot
      if (na /= 0) goto 110
      call passb5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      goto 111
      110 call passb5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      111 na = 1-na
      goto 115
      112 if (na /= 0) goto 113
      call passb (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
      goto 114
      113 call passb (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
      114 if (nac /= 0) na = 1-na
      115 l1 = l2
      iw = iw+(ip-1)*idot
   116 continue
   if (na == 0) return
   n2 = n+n
   do 117 i=1,n2
      c(i) = ch(i)
   117 continue
   return
end subroutine cfftb1 

!     --- CFFTF ---

! Forward complex Transform


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cfftf here]
subroutine cfftf (n,c,wsave)
   implicit none
   integer:: iw1, iw2, n
   _REAL_ :: c, wsave
   dimension       c(*)       ,wsave(*)
   if (n == 1) return
   iw1 = n+n+1
   iw2 = iw1+n+n
   call cfftf1 (n,c,wsave,wsave(iw1),wsave(iw2))
   return
end subroutine cfftf 

!     --- CFFTF1 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cfftf1 here]
subroutine cfftf1 (n,c,ch,wa,ifac)
   implicit none
   integer:: i, idl1, ido, idot, ifac, ip, iw, ix2, ix3, ix4, k1, &
        l1, l2, n, n2, na, nac, nf
   _REAL_ :: c, ch, wa
   dimension       ch(*)      ,c(*)       ,wa(*)      ,ifac(*)
   nf = ifac(2)
   na = 0
   l1 = 1
   iw = 1
   do 116 k1=1,nf
      ip = ifac(k1+2)
      l2 = ip*l1
      ido = n/l2
      idot = ido+ido
      idl1 = idot*l1
      if (ip /= 4) goto 103
      ix2 = iw+idot
      ix3 = ix2+idot
      if (na /= 0) goto 101
      call passf4 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3))
      goto 102
      101 call passf4 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3))
      102 na = 1-na
      goto 115
      103 if (ip /= 2) goto 106
      if (na /= 0) goto 104
      call passf2 (idot,l1,c,ch,wa(iw))
      goto 105
      104 call passf2 (idot,l1,ch,c,wa(iw))
      105 na = 1-na
      goto 115
      106 if (ip /= 3) goto 109
      ix2 = iw+idot
      if (na /= 0) goto 107
      call passf3 (idot,l1,c,ch,wa(iw),wa(ix2))
      goto 108
      107 call passf3 (idot,l1,ch,c,wa(iw),wa(ix2))
      108 na = 1-na
      goto 115
      109 if (ip /= 5) goto 112
      ix2 = iw+idot
      ix3 = ix2+idot
      ix4 = ix3+idot
      if (na /= 0) goto 110
      call passf5 (idot,l1,c,ch,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      goto 111
      110 call passf5 (idot,l1,ch,c,wa(iw),wa(ix2),wa(ix3),wa(ix4))
      111 na = 1-na
      goto 115
      112 if (na /= 0) goto 113
      call passf (nac,idot,ip,l1,idl1,c,c,c,ch,ch,wa(iw))
      goto 114
      113 call passf (nac,idot,ip,l1,idl1,ch,ch,ch,c,c,wa(iw))
      114 if (nac /= 0) na = 1-na
      115 l1 = l2
      iw = iw+(ip-1)*idot
   116 continue
   if (na == 0) return
   n2 = n+n
   do 117 i=1,n2
      c(i) = ch(i)
   117 continue
   return
end subroutine cfftf1 

!     --- CFFTI ---

! Initialization for CFFTF and CFFTB

!******************************************************************


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cffti here]
subroutine cffti(n,wsave)
   
   !******************************************************************
   
   ! subroutine cffti initializes the array wsave which is used in
   ! both cfftf and cfftb. the prime factorization of n together with
   ! a tabulation of the trigonometric functions are computed and
   ! stored in wsave.
   
   ! input parameter
   
   ! n       the length of the sequence to be transformed
   
   ! output parameter
   
   ! wsave   a work array which must be dimensioned at least 4*n+15
   !         the same work array can be used for both cfftf and cfftb
   !         as long as n remains unchanged. different wsave arrays
   !         are required for different values of n. the contents of
   !         wsave must not be changed between calls of cfftf or cfftb.

   implicit none
   integer      n
   _REAL_       wsave(*)

   integer      iw1
   integer      iw2

   if (n == 1) return
   iw1 = n+n+1
   iw2 = iw1+n+n
   call cffti1 (n,wsave(iw1),wsave(iw2))
   return
end subroutine cffti 

!     --- CFFTI1 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cffti1 here]
subroutine cffti1 (n,wa,ifac)
   implicit none
   integer:: i, i1, ib, ido, idot, ifac, ii, ip, ipm, j, k1, l1, l2, &
        ld, n, nf, nl, nq, nr, ntry, ntryh
   _REAL_ :: arg, argh, argld, fi, tpi, wa
   dimension       wa(*)      ,ifac(*)    ,ntryh(4)
   data ntryh(1),ntryh(2),ntryh(3),ntryh(4)/3,4,2,5/
   nl = n
   nf = 0
   j = 0
   101 j = j+1
   if (j-4) 102,102,103
   102 ntry = ntryh(j)
   goto 104
   103 ntry = ntry+2
   104 nq = nl/ntry
   nr = nl-ntry*nq
   if (nr) 101,105,101
   105 nf = nf+1
   ifac(nf+2) = ntry
   nl = nq
   if (ntry /= 2) goto 107
   if (nf == 1) goto 107
   do 106 i=2,nf
      ib = nf-i+2
      ifac(ib+2) = ifac(ib+1)
   106 continue
   ifac(3) = 2
   107 if (nl /= 1) goto 104
   ifac(1) = n
   ifac(2) = nf
   tpi = 6.28318530717959d0
   argh = tpi/dble(n)
   i = 2
   l1 = 1
   do 110 k1=1,nf
      ip = ifac(k1+2)
      ld = 0
      l2 = l1*ip
      ido = n/l2
      idot = ido+ido+2
      ipm = ip-1
      do 109 j=1,ipm
         i1 = i
         wa(i-1) = 1.d0
         wa(i) = 0.d0
         ld = ld+l1
         fi = 0.d0
         argld = dble(ld)*argh
         do 108 ii=4,idot,2
            i = i+2
            fi = fi+1.d0
            arg = fi*argld
            wa(i-1) = cos(arg)
            wa(i) = sin(arg)
         108 continue
         if (ip <= 5) goto 109
         wa(i1-1) = wa(i-1)
         wa(i1) = wa(i)
      109 continue
      l1 = l2
   110 continue
   return
end subroutine cffti1 

!     --- PASSB ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passb here]
subroutine passb (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
   implicit none
   integer:: i, idij, idj, idl, idl1, idlj, ido, idot, idp, ik, inc, &
        ip, ipp2, ipph, j, jc, k, l, l1, lc, nac, nt
   _REAL_ :: c1, c2, cc, ch, ch2, wa, wai, war
   dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          , &
         c1(ido,l1,ip)          ,wa(*)      ,c2(idl1,ip), &
         ch2(idl1,ip)
   idot = ido/2
   nt = ip*idl1
   ipp2 = ip+2
   ipph = (ip+1)/2
   idp = ip*ido
   
   if (ido < l1) goto 106
   do 103 j=2,ipph
      jc = ipp2-j
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
            ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
         101 continue
      102 continue
   103 continue
   do 105 k=1,l1
      do 104 i=1,ido
         ch(i,k,1) = cc(i,1,k)
      104 continue
   105 continue
   goto 112
   106 do 109 j=2,ipph
      jc = ipp2-j
      do 108 i=1,ido
         do 107 k=1,l1
            ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
            ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
         107 continue
      108 continue
   109 continue
   do 111 i=1,ido
      do 110 k=1,l1
         ch(i,k,1) = cc(i,1,k)
      110 continue
   111 continue
   112 idl = 2-ido
   inc = 0
   do 116 l=2,ipph
      lc = ipp2-l
      idl = idl+ido
      do 113 ik=1,idl1
         c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
         c2(ik,lc) = wa(idl)*ch2(ik,ip)
      113 continue
      idlj = idl
      inc = inc+ido
      do 115 j=3,ipph
         jc = ipp2-j
         idlj = idlj+inc
         if (idlj > idp) idlj = idlj-idp
         war = wa(idlj-1)
         wai = wa(idlj)
         do 114 ik=1,idl1
            c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
            c2(ik,lc) = c2(ik,lc)+wai*ch2(ik,jc)
         114 continue
      115 continue
   116 continue
   do 118 j=2,ipph
      do 117 ik=1,idl1
         ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
      117 continue
   118 continue
   do 120 j=2,ipph
      jc = ipp2-j
      do 119 ik=2,idl1,2
         ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
         ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
         ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
         ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
      119 continue
   120 continue
   nac = 1
   if (ido == 2) return
   nac = 0
   do 121 ik=1,idl1
      c2(ik,1) = ch2(ik,1)
   121 continue
   do 123 j=2,ip
      do 122 k=1,l1
         c1(1,k,j) = ch(1,k,j)
         c1(2,k,j) = ch(2,k,j)
      122 continue
   123 continue
   if (idot > l1) goto 127
   idij = 0
   do 126 j=2,ip
      idij = idij+2
      do 125 i=4,ido,2
         idij = idij+2
         do 124 k=1,l1
            c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
            c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
         124 continue
      125 continue
   126 continue
   return
   127 idj = 2-ido
   do 130 j=2,ip
      idj = idj+ido
      do 129 k=1,l1
         idij = idj
         do 128 i=4,ido,2
            idij = idij+2
            c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)-wa(idij)*ch(i,k,j)
            c1(i,k,j) = wa(idij-1)*ch(i,k,j)+wa(idij)*ch(i-1,k,j)
         128 continue
      129 continue
   130 continue
   return
end subroutine passb 

!     ---- PASSB2 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passb2 here]
subroutine passb2 (ido,l1,cc,ch,wa1)
   implicit none
   integer:: i, ido, k, l1
   _REAL_ :: cc, ch, ti2, tr2, wa1
   dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           , &
         wa1(*)
   if (ido > 2) goto 102
   do 101 k=1,l1
      ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
      ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
      ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
      ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
         tr2 = cc(i-1,1,k)-cc(i-1,2,k)
         ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
         ti2 = cc(i,1,k)-cc(i,2,k)
         ch(i,k,2) = wa1(i-1)*ti2+wa1(i)*tr2
         ch(i-1,k,2) = wa1(i-1)*tr2-wa1(i)*ti2
      103 continue
   104 continue
   return
end subroutine passb2 

!     --- PASSB3 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passb3 here]
subroutine passb3 (ido,l1,cc,ch,wa1,wa2)
   implicit none
   integer:: i, ido, k, l1
   _REAL_ :: cc, ch, ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, &
        taur, ti2, tr2, wa1, wa2
   dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           , &
         wa1(*)     ,wa2(*)
   data taur,taui /-.5d0,.866025403784439d0/
   if (ido /= 2) goto 102
   do 101 k=1,l1
      tr2 = cc(1,2,k)+cc(1,3,k)
      cr2 = cc(1,1,k)+taur*tr2
      ch(1,k,1) = cc(1,1,k)+tr2
      ti2 = cc(2,2,k)+cc(2,3,k)
      ci2 = cc(2,1,k)+taur*ti2
      ch(2,k,1) = cc(2,1,k)+ti2
      cr3 = taui*(cc(1,2,k)-cc(1,3,k))
      ci3 = taui*(cc(2,2,k)-cc(2,3,k))
      ch(1,k,2) = cr2-ci3
      ch(1,k,3) = cr2+ci3
      ch(2,k,2) = ci2+cr3
      ch(2,k,3) = ci2-cr3
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         tr2 = cc(i-1,2,k)+cc(i-1,3,k)
         cr2 = cc(i-1,1,k)+taur*tr2
         ch(i-1,k,1) = cc(i-1,1,k)+tr2
         ti2 = cc(i,2,k)+cc(i,3,k)
         ci2 = cc(i,1,k)+taur*ti2
         ch(i,k,1) = cc(i,1,k)+ti2
         cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
         ci3 = taui*(cc(i,2,k)-cc(i,3,k))
         dr2 = cr2-ci3
         dr3 = cr2+ci3
         di2 = ci2+cr3
         di3 = ci2-cr3
         ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
         ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
         ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
         ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
      103 continue
   104 continue
   return
end subroutine passb3 

!     --- PASSB4 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passb4 here]
subroutine passb4 (ido,l1,cc,ch,wa1,wa2,wa3)
   implicit none
   integer:: i, ido, k, l1
   _REAL_ :: cc, ch, ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, &
        ti4, tr1, tr2, tr3, tr4, wa1, wa2, wa3
   dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           , &
         wa1(*)     ,wa2(*)     ,wa3(*)
   if (ido /= 2) goto 102
   do 101 k=1,l1
      ti1 = cc(2,1,k)-cc(2,3,k)
      ti2 = cc(2,1,k)+cc(2,3,k)
      tr4 = cc(2,4,k)-cc(2,2,k)
      ti3 = cc(2,2,k)+cc(2,4,k)
      tr1 = cc(1,1,k)-cc(1,3,k)
      tr2 = cc(1,1,k)+cc(1,3,k)
      ti4 = cc(1,2,k)-cc(1,4,k)
      tr3 = cc(1,2,k)+cc(1,4,k)
      ch(1,k,1) = tr2+tr3
      ch(1,k,3) = tr2-tr3
      ch(2,k,1) = ti2+ti3
      ch(2,k,3) = ti2-ti3
      ch(1,k,2) = tr1+tr4
      ch(1,k,4) = tr1-tr4
      ch(2,k,2) = ti1+ti4
      ch(2,k,4) = ti1-ti4
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ti1 = cc(i,1,k)-cc(i,3,k)
         ti2 = cc(i,1,k)+cc(i,3,k)
         ti3 = cc(i,2,k)+cc(i,4,k)
         tr4 = cc(i,4,k)-cc(i,2,k)
         tr1 = cc(i-1,1,k)-cc(i-1,3,k)
         tr2 = cc(i-1,1,k)+cc(i-1,3,k)
         ti4 = cc(i-1,2,k)-cc(i-1,4,k)
         tr3 = cc(i-1,2,k)+cc(i-1,4,k)
         ch(i-1,k,1) = tr2+tr3
         cr3 = tr2-tr3
         ch(i,k,1) = ti2+ti3
         ci3 = ti2-ti3
         cr2 = tr1+tr4
         cr4 = tr1-tr4
         ci2 = ti1+ti4
         ci4 = ti1-ti4
         ch(i-1,k,2) = wa1(i-1)*cr2-wa1(i)*ci2
         ch(i,k,2) = wa1(i-1)*ci2+wa1(i)*cr2
         ch(i-1,k,3) = wa2(i-1)*cr3-wa2(i)*ci3
         ch(i,k,3) = wa2(i-1)*ci3+wa2(i)*cr3
         ch(i-1,k,4) = wa3(i-1)*cr4-wa3(i)*ci4
         ch(i,k,4) = wa3(i-1)*ci4+wa3(i)*cr4
      103 continue
   104 continue
   return
end subroutine passb4 

!     --- PASSB5 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passb5 here]
subroutine passb5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
   implicit none
   integer:: i, ido, k, l1
   _REAL_ :: cc, ch, ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, &
        di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3, ti4, &
        ti5, tr11, tr12, tr2, tr3, tr4, tr5, wa1, wa2, wa3, wa4
   dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           , &
         wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
   data tr11,ti11,tr12,ti12 /.309016994374947d0, &
         .951056516295154d0, &
         -.809016994374947d0,.587785252292473d0/
   if (ido /= 2) goto 102
   do 101 k=1,l1
      ti5 = cc(2,2,k)-cc(2,5,k)
      ti2 = cc(2,2,k)+cc(2,5,k)
      ti4 = cc(2,3,k)-cc(2,4,k)
      ti3 = cc(2,3,k)+cc(2,4,k)
      tr5 = cc(1,2,k)-cc(1,5,k)
      tr2 = cc(1,2,k)+cc(1,5,k)
      tr4 = cc(1,3,k)-cc(1,4,k)
      tr3 = cc(1,3,k)+cc(1,4,k)
      ch(1,k,1) = cc(1,1,k)+tr2+tr3
      ch(2,k,1) = cc(2,1,k)+ti2+ti3
      cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
      ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
      cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
      ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
      cr5 = ti11*tr5+ti12*tr4
      ci5 = ti11*ti5+ti12*ti4
      cr4 = ti12*tr5-ti11*tr4
      ci4 = ti12*ti5-ti11*ti4
      ch(1,k,2) = cr2-ci5
      ch(1,k,5) = cr2+ci5
      ch(2,k,2) = ci2+cr5
      ch(2,k,3) = ci3+cr4
      ch(1,k,3) = cr3-ci4
      ch(1,k,4) = cr3+ci4
      ch(2,k,4) = ci3-cr4
      ch(2,k,5) = ci2-cr5
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ti5 = cc(i,2,k)-cc(i,5,k)
         ti2 = cc(i,2,k)+cc(i,5,k)
         ti4 = cc(i,3,k)-cc(i,4,k)
         ti3 = cc(i,3,k)+cc(i,4,k)
         tr5 = cc(i-1,2,k)-cc(i-1,5,k)
         tr2 = cc(i-1,2,k)+cc(i-1,5,k)
         tr4 = cc(i-1,3,k)-cc(i-1,4,k)
         tr3 = cc(i-1,3,k)+cc(i-1,4,k)
         ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
         ch(i,k,1) = cc(i,1,k)+ti2+ti3
         cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
         ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
         cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
         ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         dr3 = cr3-ci4
         dr4 = cr3+ci4
         di3 = ci3+cr4
         di4 = ci3-cr4
         dr5 = cr2+ci5
         dr2 = cr2-ci5
         di5 = ci2-cr5
         di2 = ci2+cr5
         ch(i-1,k,2) = wa1(i-1)*dr2-wa1(i)*di2
         ch(i,k,2) = wa1(i-1)*di2+wa1(i)*dr2
         ch(i-1,k,3) = wa2(i-1)*dr3-wa2(i)*di3
         ch(i,k,3) = wa2(i-1)*di3+wa2(i)*dr3
         ch(i-1,k,4) = wa3(i-1)*dr4-wa3(i)*di4
         ch(i,k,4) = wa3(i-1)*di4+wa3(i)*dr4
         ch(i-1,k,5) = wa4(i-1)*dr5-wa4(i)*di5
         ch(i,k,5) = wa4(i-1)*di5+wa4(i)*dr5
      103 continue
   104 continue
   return
end subroutine passb5 

!     --- PASSF ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passf here]
subroutine passf (nac,ido,ip,l1,idl1,cc,c1,c2,ch,ch2,wa)
   implicit none
   integer:: i, idij, idj, idl, idl1, idlj, ido, idot, idp, ik, inc, &
        ip, ipp2, ipph, j, jc, k, l, l1, lc, nac, nt
   _REAL_ :: c1, c2, cc, ch, ch2, wa, wai, war
   dimension       ch(ido,l1,ip)          ,cc(ido,ip,l1)          , &
         c1(ido,l1,ip)          ,wa(*)      ,c2(idl1,ip), &
         ch2(idl1,ip)
   idot = ido/2
   nt = ip*idl1
   ipp2 = ip+2
   ipph = (ip+1)/2
   idp = ip*ido
   
   if (ido < l1) goto 106
   do 103 j=2,ipph
      jc = ipp2-j
      do 102 k=1,l1
         do 101 i=1,ido
            ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
            ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
         101 continue
      102 continue
   103 continue
   do 105 k=1,l1
      do 104 i=1,ido
         ch(i,k,1) = cc(i,1,k)
      104 continue
   105 continue
   goto 112
   106 do 109 j=2,ipph
      jc = ipp2-j
      do 108 i=1,ido
         do 107 k=1,l1
            ch(i,k,j) = cc(i,j,k)+cc(i,jc,k)
            ch(i,k,jc) = cc(i,j,k)-cc(i,jc,k)
         107 continue
      108 continue
   109 continue
   do 111 i=1,ido
      do 110 k=1,l1
         ch(i,k,1) = cc(i,1,k)
      110 continue
   111 continue
   112 idl = 2-ido
   inc = 0
   do 116 l=2,ipph
      lc = ipp2-l
      idl = idl+ido
      do 113 ik=1,idl1
         c2(ik,l) = ch2(ik,1)+wa(idl-1)*ch2(ik,2)
         c2(ik,lc) = -wa(idl)*ch2(ik,ip)
      113 continue
      idlj = idl
      inc = inc+ido
      do 115 j=3,ipph
         jc = ipp2-j
         idlj = idlj+inc
         if (idlj > idp) idlj = idlj-idp
         war = wa(idlj-1)
         wai = wa(idlj)
         do 114 ik=1,idl1
            c2(ik,l) = c2(ik,l)+war*ch2(ik,j)
            c2(ik,lc) = c2(ik,lc)-wai*ch2(ik,jc)
         114 continue
      115 continue
   116 continue
   do 118 j=2,ipph
      do 117 ik=1,idl1
         ch2(ik,1) = ch2(ik,1)+ch2(ik,j)
      117 continue
   118 continue
   do 120 j=2,ipph
      jc = ipp2-j
      do 119 ik=2,idl1,2
         ch2(ik-1,j) = c2(ik-1,j)-c2(ik,jc)
         ch2(ik-1,jc) = c2(ik-1,j)+c2(ik,jc)
         ch2(ik,j) = c2(ik,j)+c2(ik-1,jc)
         ch2(ik,jc) = c2(ik,j)-c2(ik-1,jc)
      119 continue
   120 continue
   nac = 1
   if (ido == 2) return
   nac = 0
   do 121 ik=1,idl1
      c2(ik,1) = ch2(ik,1)
   121 continue
   do 123 j=2,ip
      do 122 k=1,l1
         c1(1,k,j) = ch(1,k,j)
         c1(2,k,j) = ch(2,k,j)
      122 continue
   123 continue
   if (idot > l1) goto 127
   idij = 0
   do 126 j=2,ip
      idij = idij+2
      do 125 i=4,ido,2
         idij = idij+2
         do 124 k=1,l1
            c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
            c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
         124 continue
      125 continue
   126 continue
   return
   127 idj = 2-ido
   do 130 j=2,ip
      idj = idj+ido
      do 129 k=1,l1
         idij = idj
         do 128 i=4,ido,2
            idij = idij+2
            c1(i-1,k,j) = wa(idij-1)*ch(i-1,k,j)+wa(idij)*ch(i,k,j)
            c1(i,k,j) = wa(idij-1)*ch(i,k,j)-wa(idij)*ch(i-1,k,j)
         128 continue
      129 continue
   130 continue
   return
end subroutine passf 

!     --- PASSF2 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passf2 here]
subroutine passf2 (ido,l1,cc,ch,wa1)
   implicit none
   integer:: i, ido, k, l1
   _REAL_ :: cc, ch, ti2, tr2, wa1
   dimension       cc(ido,2,l1)           ,ch(ido,l1,2)           , &
         wa1(*)
   if (ido > 2) goto 102
   do 101 k=1,l1
      ch(1,k,1) = cc(1,1,k)+cc(1,2,k)
      ch(1,k,2) = cc(1,1,k)-cc(1,2,k)
      ch(2,k,1) = cc(2,1,k)+cc(2,2,k)
      ch(2,k,2) = cc(2,1,k)-cc(2,2,k)
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ch(i-1,k,1) = cc(i-1,1,k)+cc(i-1,2,k)
         tr2 = cc(i-1,1,k)-cc(i-1,2,k)
         ch(i,k,1) = cc(i,1,k)+cc(i,2,k)
         ti2 = cc(i,1,k)-cc(i,2,k)
         ch(i,k,2) = wa1(i-1)*ti2-wa1(i)*tr2
         ch(i-1,k,2) = wa1(i-1)*tr2+wa1(i)*ti2
      103 continue
   104 continue
   return
end subroutine passf2 

!     --- PASSF3 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passf3 here]
subroutine passf3 (ido,l1,cc,ch,wa1,wa2)
   implicit none
   integer:: i, ido, k, l1
   _REAL_ :: cc, ch, ci2, ci3, cr2, cr3, di2, di3, dr2, dr3, taui, &
        taur, ti2, tr2, wa1, wa2
   dimension       cc(ido,3,l1)           ,ch(ido,l1,3)           , &
         wa1(*)     ,wa2(*)
   data taur,taui /-.5d0,-.866025403784439d0/
   if (ido /= 2) goto 102
   do 101 k=1,l1
      tr2 = cc(1,2,k)+cc(1,3,k)
      cr2 = cc(1,1,k)+taur*tr2
      ch(1,k,1) = cc(1,1,k)+tr2
      ti2 = cc(2,2,k)+cc(2,3,k)
      ci2 = cc(2,1,k)+taur*ti2
      ch(2,k,1) = cc(2,1,k)+ti2
      cr3 = taui*(cc(1,2,k)-cc(1,3,k))
      ci3 = taui*(cc(2,2,k)-cc(2,3,k))
      ch(1,k,2) = cr2-ci3
      ch(1,k,3) = cr2+ci3
      ch(2,k,2) = ci2+cr3
      ch(2,k,3) = ci2-cr3
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         tr2 = cc(i-1,2,k)+cc(i-1,3,k)
         cr2 = cc(i-1,1,k)+taur*tr2
         ch(i-1,k,1) = cc(i-1,1,k)+tr2
         ti2 = cc(i,2,k)+cc(i,3,k)
         ci2 = cc(i,1,k)+taur*ti2
         ch(i,k,1) = cc(i,1,k)+ti2
         cr3 = taui*(cc(i-1,2,k)-cc(i-1,3,k))
         ci3 = taui*(cc(i,2,k)-cc(i,3,k))
         dr2 = cr2-ci3
         dr3 = cr2+ci3
         di2 = ci2+cr3
         di3 = ci2-cr3
         ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
         ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
         ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
         ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
      103 continue
   104 continue
   return
end subroutine passf3 

!     --- PASSF4 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passf4 here]
subroutine passf4 (ido,l1,cc,ch,wa1,wa2,wa3)
   implicit none
   integer:: i, ido, k, l1
   _REAL_ :: cc, ch, ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, &
        ti4, tr1, tr2, tr3, tr4, wa1, wa2, wa3
   dimension       cc(ido,4,l1)           ,ch(ido,l1,4)           , &
         wa1(*)     ,wa2(*)     ,wa3(*)
   if (ido /= 2) goto 102
   do 101 k=1,l1
      ti1 = cc(2,1,k)-cc(2,3,k)
      ti2 = cc(2,1,k)+cc(2,3,k)
      tr4 = cc(2,2,k)-cc(2,4,k)
      ti3 = cc(2,2,k)+cc(2,4,k)
      tr1 = cc(1,1,k)-cc(1,3,k)
      tr2 = cc(1,1,k)+cc(1,3,k)
      ti4 = cc(1,4,k)-cc(1,2,k)
      tr3 = cc(1,2,k)+cc(1,4,k)
      ch(1,k,1) = tr2+tr3
      ch(1,k,3) = tr2-tr3
      ch(2,k,1) = ti2+ti3
      ch(2,k,3) = ti2-ti3
      ch(1,k,2) = tr1+tr4
      ch(1,k,4) = tr1-tr4
      ch(2,k,2) = ti1+ti4
      ch(2,k,4) = ti1-ti4
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ti1 = cc(i,1,k)-cc(i,3,k)
         ti2 = cc(i,1,k)+cc(i,3,k)
         ti3 = cc(i,2,k)+cc(i,4,k)
         tr4 = cc(i,2,k)-cc(i,4,k)
         tr1 = cc(i-1,1,k)-cc(i-1,3,k)
         tr2 = cc(i-1,1,k)+cc(i-1,3,k)
         ti4 = cc(i-1,4,k)-cc(i-1,2,k)
         tr3 = cc(i-1,2,k)+cc(i-1,4,k)
         ch(i-1,k,1) = tr2+tr3
         cr3 = tr2-tr3
         ch(i,k,1) = ti2+ti3
         ci3 = ti2-ti3
         cr2 = tr1+tr4
         cr4 = tr1-tr4
         ci2 = ti1+ti4
         ci4 = ti1-ti4
         ch(i-1,k,2) = wa1(i-1)*cr2+wa1(i)*ci2
         ch(i,k,2) = wa1(i-1)*ci2-wa1(i)*cr2
         ch(i-1,k,3) = wa2(i-1)*cr3+wa2(i)*ci3
         ch(i,k,3) = wa2(i-1)*ci3-wa2(i)*cr3
         ch(i-1,k,4) = wa3(i-1)*cr4+wa3(i)*ci4
         ch(i,k,4) = wa3(i-1)*ci4-wa3(i)*cr4
      103 continue
   104 continue
   return
end subroutine passf4 

!     --- PASSF5 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine passf5 here]
subroutine passf5 (ido,l1,cc,ch,wa1,wa2,wa3,wa4)
   implicit none
   integer:: i, ido, k, l1
   _REAL_ :: cc, ch, ci2, ci3, ci4, ci5, cr2, cr3, cr4, cr5, di2, &
        di3, di4, di5, dr2, dr3, dr4, dr5, ti11, ti12, ti2, ti3, ti4, &
        ti5, tr11, tr12, tr2, tr3, tr4, tr5, wa1, wa2, wa3, wa4
   dimension       cc(ido,5,l1)           ,ch(ido,l1,5)           , &
         wa1(*)     ,wa2(*)     ,wa3(*)     ,wa4(*)
   data tr11,ti11,tr12,ti12 /.309016994374947d0, &
         -.951056516295154d0, &
         -.809016994374947d0,-.587785252292473d0/
   if (ido /= 2) goto 102
   do 101 k=1,l1
      ti5 = cc(2,2,k)-cc(2,5,k)
      ti2 = cc(2,2,k)+cc(2,5,k)
      ti4 = cc(2,3,k)-cc(2,4,k)
      ti3 = cc(2,3,k)+cc(2,4,k)
      tr5 = cc(1,2,k)-cc(1,5,k)
      tr2 = cc(1,2,k)+cc(1,5,k)
      tr4 = cc(1,3,k)-cc(1,4,k)
      tr3 = cc(1,3,k)+cc(1,4,k)
      ch(1,k,1) = cc(1,1,k)+tr2+tr3
      ch(2,k,1) = cc(2,1,k)+ti2+ti3
      cr2 = cc(1,1,k)+tr11*tr2+tr12*tr3
      ci2 = cc(2,1,k)+tr11*ti2+tr12*ti3
      cr3 = cc(1,1,k)+tr12*tr2+tr11*tr3
      ci3 = cc(2,1,k)+tr12*ti2+tr11*ti3
      cr5 = ti11*tr5+ti12*tr4
      ci5 = ti11*ti5+ti12*ti4
      cr4 = ti12*tr5-ti11*tr4
      ci4 = ti12*ti5-ti11*ti4
      ch(1,k,2) = cr2-ci5
      ch(1,k,5) = cr2+ci5
      ch(2,k,2) = ci2+cr5
      ch(2,k,3) = ci3+cr4
      ch(1,k,3) = cr3-ci4
      ch(1,k,4) = cr3+ci4
      ch(2,k,4) = ci3-cr4
      ch(2,k,5) = ci2-cr5
   101 continue
   return
   102 do 104 k=1,l1
      do 103 i=2,ido,2
         ti5 = cc(i,2,k)-cc(i,5,k)
         ti2 = cc(i,2,k)+cc(i,5,k)
         ti4 = cc(i,3,k)-cc(i,4,k)
         ti3 = cc(i,3,k)+cc(i,4,k)
         tr5 = cc(i-1,2,k)-cc(i-1,5,k)
         tr2 = cc(i-1,2,k)+cc(i-1,5,k)
         tr4 = cc(i-1,3,k)-cc(i-1,4,k)
         tr3 = cc(i-1,3,k)+cc(i-1,4,k)
         ch(i-1,k,1) = cc(i-1,1,k)+tr2+tr3
         ch(i,k,1) = cc(i,1,k)+ti2+ti3
         cr2 = cc(i-1,1,k)+tr11*tr2+tr12*tr3
         ci2 = cc(i,1,k)+tr11*ti2+tr12*ti3
         cr3 = cc(i-1,1,k)+tr12*tr2+tr11*tr3
         ci3 = cc(i,1,k)+tr12*ti2+tr11*ti3
         cr5 = ti11*tr5+ti12*tr4
         ci5 = ti11*ti5+ti12*ti4
         cr4 = ti12*tr5-ti11*tr4
         ci4 = ti12*ti5-ti11*ti4
         dr3 = cr3-ci4
         dr4 = cr3+ci4
         di3 = ci3+cr4
         di4 = ci3-cr4
         dr5 = cr2+ci5
         dr2 = cr2-ci5
         di5 = ci2-cr5
         di2 = ci2+cr5
         ch(i-1,k,2) = wa1(i-1)*dr2+wa1(i)*di2
         ch(i,k,2) = wa1(i-1)*di2-wa1(i)*dr2
         ch(i-1,k,3) = wa2(i-1)*dr3+wa2(i)*di3
         ch(i,k,3) = wa2(i-1)*di3-wa2(i)*dr3
         ch(i-1,k,4) = wa3(i-1)*dr4+wa3(i)*di4
         ch(i,k,4) = wa3(i-1)*di4-wa3(i)*dr4
         ch(i-1,k,5) = wa4(i-1)*dr5+wa4(i)*di5
         ch(i,k,5) = wa4(i-1)*di5-wa4(i)*dr5
      103 continue
   104 continue
   return
end subroutine passf5 

#ifndef MPI
!     --- PUBZ3D ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pubz3d here]
subroutine pubz3d(isign,n1,n2,n3,w,ld1,ld2,table,ntable, &
      work,nwork)
   implicit none

   integer n1,n2,n3,ld1,ld2,isign,ntable,nwork
#ifdef DPREC
   double complex w(ld1,ld2,n3)
   double complex work( nwork)
#else
   complex w(ld1,ld2,n3)
   complex work( nwork)
#endif
   _REAL_ table(ntable,3)

   integer i,j,k
   !     ...note: ntable should be 4*max(n1,n2,n3) +15
   !              nwork should be max(n1,n2,n3)

   !     ...transform along X  first
   
   do 100 k = 1, n3
      do 90 j = 1, n2
         do 70 i = 1,n1
            work(i) = w(i,j,k)
         70 continue
         if ( isign == -1) call cfftf(n1,work,table(1,1))
         if ( isign == 1) call cfftb(n1,work,table(1,1))
         do 80 i = 1,n1
            w(i,j,k) = work(i)
         80 continue
      90 continue
   100 continue
   
   !     ...transform along Y then
   
   do 200 k = 1,n3
      do 190 i = 1,n1
         do 170 j = 1,n2
            work(j) = w(i,j,k)
         170 continue
         if ( isign == -1) call cfftf(n2,work,table(1,2))
         if ( isign == 1) call cfftb(n2,work,table(1,2))
         do 180 j = 1,n2
            w(i,j,k) = work(j)
         180 continue
      190 continue
   200 continue
   
   !     ...transform along Z finally
   
   do 300 i = 1, n1
      do 290 j = 1, n2
         do 270 k = 1,n3
            work(k) = w(i,j,k)
         270 continue
         if ( isign == -1) call cfftf(n3,work,table(1,3))
         if ( isign == 1) call cfftb(n3,work,table(1,3))
         do 280 k = 1,n3
            w(i,j,k) = work(k)
         280 continue
      290 continue
   300 continue

   return
end subroutine pubz3d 

!     --- PUBZ3DI ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pubz3di here]
subroutine pubz3di(n1,n2,n3,table,ntable)
   implicit none
   integer n1,n2,n3,ntable
   _REAL_ table(ntable,3)
   !     NOTE: ntable should be 4*max(n1,n2,n3) +15

   call cffti(n1,table(1,1))
   call cffti(n2,table(1,2))
   call cffti(n3,table(1,3))

   return
end subroutine pubz3di 
#endif /* MPI */

!     --- RADF2 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine radf2 here]
SUBROUTINE RADF2 (IDO,L1,CC,CH,WA1)
   implicit _REAL_ (a-h,o-z)
   DIMENSION       CH(IDO,2,L1)           ,CC(IDO,L1,2)           , &
      WA1(*)
   DO 101 K=1,L1
      CH(1,1,K) = CC(1,K,1)+CC(1,K,2)
         CH(IDO,2,K) = CC(1,K,1)-CC(1,K,2)
   101 end do
   IF (IDO-2) 107,105,102
   102 IDP2 = IDO+2
   DO 104 K=1,L1
      DO 103 I=3,IDO,2
         IC = IDP2-I
         TR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
         TI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
         CH(I,1,K) = CC(I,K,1)+TI2
         CH(IC,2,K) = TI2-CC(I,K,1)
         CH(I-1,1,K) = CC(I-1,K,1)+TR2
         CH(IC-1,2,K) = CC(I-1,K,1)-TR2
      103 end do
   104 end do
   IF (MOD(IDO,2) == 1) RETURN
   105 DO 106 K=1,L1
       CH(1,2,K) = -CC(IDO,K,2)
       CH(IDO,1,K) = CC(IDO,K,1)
   106 end do
   107 RETURN
end SUBROUTINE RADF2

SUBROUTINE RADF3 (IDO,L1,CC,CH,WA1,WA2)
   implicit _REAL_ (a-h,o-z)
   DIMENSION       CH(IDO,3,L1)           ,CC(IDO,L1,3)           , &
   WA1(*)     ,WA2(*)
   DATA TAUR,TAUI /-.5,.866025403784439/
   DO 101 K=1,L1
      CR2 = CC(1,K,2)+CC(1,K,3)
      CH(1,1,K) = CC(1,K,1)+CR2
      CH(1,3,K) = TAUI*(CC(1,K,3)-CC(1,K,2))
      CH(IDO,2,K) = CC(1,K,1)+TAUR*CR2
   101 end do
   IF (IDO == 1) RETURN
   IDP2 = IDO+2
   DO 103 K=1,L1
      DO 102 I=3,IDO,2
         IC = IDP2-I
         DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
         DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
         DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
         DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
         CR2 = DR2+DR3
         CI2 = DI2+DI3
         CH(I-1,1,K) = CC(I-1,K,1)+CR2
         CH(I,1,K) = CC(I,K,1)+CI2
         TR2 = CC(I-1,K,1)+TAUR*CR2
         TI2 = CC(I,K,1)+TAUR*CI2
         TR3 = TAUI*(DI2-DI3)
         TI3 = TAUI*(DR3-DR2)
         CH(I-1,3,K) = TR2+TR3
         CH(IC-1,2,K) = TR2-TR3
         CH(I,3,K) = TI2+TI3
         CH(IC,2,K) = TI3-TI2
      102 end do
   103 end do
   RETURN
end SUBROUTINE RADF3

!     --- RADF3 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine radf3 here]
SUBROUTINE RADF4 (IDO,L1,CC,CH,WA1,WA2,WA3)
   implicit _REAL_ (a-h,o-z)
   DIMENSION       CC(IDO,L1,4)           ,CH(IDO,4,L1)           , &
   WA1(*)     ,WA2(*)     ,WA3(*)
   DATA HSQT2 /.7071067811865475/
   DO 101 K=1,L1
      TR1 = CC(1,K,2)+CC(1,K,4)
      TR2 = CC(1,K,1)+CC(1,K,3)
      CH(1,1,K) = TR1+TR2
      CH(IDO,4,K) = TR2-TR1
      CH(IDO,2,K) = CC(1,K,1)-CC(1,K,3)
      CH(1,3,K) = CC(1,K,4)-CC(1,K,2)
   101 end do
   IF (IDO-2) 107,105,102
   102 IDP2 = IDO+2
   DO 104 K=1,L1
      DO 103 I=3,IDO,2
         IC = IDP2-I
         CR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
         CI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
         CR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
         CI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
         CR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
         CI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
         TR1 = CR2+CR4
         TR4 = CR4-CR2
         TI1 = CI2+CI4
         TI4 = CI2-CI4
         TI2 = CC(I,K,1)+CI3
         TI3 = CC(I,K,1)-CI3
         TR2 = CC(I-1,K,1)+CR3
         TR3 = CC(I-1,K,1)-CR3
         CH(I-1,1,K) = TR1+TR2
         CH(IC-1,4,K) = TR2-TR1
         CH(I,1,K) = TI1+TI2
         CH(IC,4,K) = TI1-TI2
         CH(I-1,3,K) = TI4+TR3
         CH(IC-1,2,K) = TR3-TI4
         CH(I,3,K) = TR4+TI3
         CH(IC,2,K) = TR4-TI3
      103 end do
   104 end do
   IF (MOD(IDO,2) == 1) RETURN
   105 CONTINUE
   DO 106 K=1,L1
      TI1 = -HSQT2*(CC(IDO,K,2)+CC(IDO,K,4))
      TR1 = HSQT2*(CC(IDO,K,2)-CC(IDO,K,4))
      CH(IDO,1,K) = TR1+CC(IDO,K,1)
      CH(IDO,3,K) = CC(IDO,K,1)-TR1
      CH(1,2,K) = TI1-CC(IDO,K,3)
      CH(1,4,K) = TI1+CC(IDO,K,3)
   106 end do
   107 RETURN
end SUBROUTINE RADF4

!     --- RADF4 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine radf4 here]
SUBROUTINE RADF5 (IDO,L1,CC,CH,WA1,WA2,WA3,WA4)
   implicit _REAL_ (a-h,o-z)
   DIMENSION       CC(IDO,L1,5)           ,CH(IDO,5,L1)           , &
   WA1(*)     ,WA2(*)     ,WA3(*)     ,WA4(*)
   DATA TR11,TI11,TR12,TI12 /.309016994374947,.951056516295154, &
   -.809016994374947,.587785252292473/
   DO 101 K=1,L1
      CR2 = CC(1,K,5)+CC(1,K,2)
      CI5 = CC(1,K,5)-CC(1,K,2)
      CR3 = CC(1,K,4)+CC(1,K,3)
      CI4 = CC(1,K,4)-CC(1,K,3)
      CH(1,1,K) = CC(1,K,1)+CR2+CR3
      CH(IDO,2,K) = CC(1,K,1)+TR11*CR2+TR12*CR3
      CH(1,3,K) = TI11*CI5+TI12*CI4
      CH(IDO,4,K) = CC(1,K,1)+TR12*CR2+TR11*CR3
      CH(1,5,K) = TI12*CI5-TI11*CI4
   101 end do
   IF (IDO == 1) RETURN
   IDP2 = IDO+2
   DO 103 K=1,L1
      DO 102 I=3,IDO,2
         IC = IDP2-I
         DR2 = WA1(I-2)*CC(I-1,K,2)+WA1(I-1)*CC(I,K,2)
         DI2 = WA1(I-2)*CC(I,K,2)-WA1(I-1)*CC(I-1,K,2)
         DR3 = WA2(I-2)*CC(I-1,K,3)+WA2(I-1)*CC(I,K,3)
         DI3 = WA2(I-2)*CC(I,K,3)-WA2(I-1)*CC(I-1,K,3)
         DR4 = WA3(I-2)*CC(I-1,K,4)+WA3(I-1)*CC(I,K,4)
         DI4 = WA3(I-2)*CC(I,K,4)-WA3(I-1)*CC(I-1,K,4)
         DR5 = WA4(I-2)*CC(I-1,K,5)+WA4(I-1)*CC(I,K,5)
         DI5 = WA4(I-2)*CC(I,K,5)-WA4(I-1)*CC(I-1,K,5)
         CR2 = DR2+DR5
         CI5 = DR5-DR2
         CR5 = DI2-DI5
         CI2 = DI2+DI5
         CR3 = DR3+DR4
         CI4 = DR4-DR3
         CR4 = DI3-DI4
         CI3 = DI3+DI4
         CH(I-1,1,K) = CC(I-1,K,1)+CR2+CR3
         CH(I,1,K) = CC(I,K,1)+CI2+CI3
         TR2 = CC(I-1,K,1)+TR11*CR2+TR12*CR3
         TI2 = CC(I,K,1)+TR11*CI2+TR12*CI3
         TR3 = CC(I-1,K,1)+TR12*CR2+TR11*CR3
         TI3 = CC(I,K,1)+TR12*CI2+TR11*CI3
         TR5 = TI11*CR5+TI12*CR4
         TI5 = TI11*CI5+TI12*CI4
         TR4 = TI12*CR5-TI11*CR4
         TI4 = TI12*CI5-TI11*CI4
         CH(I-1,3,K) = TR2+TR5
         CH(IC-1,2,K) = TR2-TR5
         CH(I,3,K) = TI2+TI5
         CH(IC,2,K) = TI5-TI2
         CH(I-1,5,K) = TR3+TR4
         CH(IC-1,4,K) = TR3-TR4
         CH(I,5,K) = TI3+TI4
         CH(IC,4,K) = TI4-TI3
      102 end do
   103 end do
   RETURN
 end SUBROUTINE RADF5
 
!     --- RADF5 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine radf5 here]
SUBROUTINE RADFG (IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,WA)
   implicit _REAL_ (a-h,o-z)
   DIMENSION       CH(IDO,L1,IP)          ,CC(IDO,IP,L1)          , &
   C1(IDO,L1,IP)          ,C2(IDL1,IP), &
   CH2(IDL1,IP)           ,WA(*)
   DATA TPI/6.28318530717959/
   ARG = TPI/FLOAT(IP)
   DCP = COS(ARG)
   DSP = SIN(ARG)
   IPPH = (IP+1)/2
   IPP2 = IP+2
   IDP2 = IDO+2
   NBD = (IDO-1)/2
   IF (IDO == 1) GO TO 119
   DO 101 IK=1,IDL1
       CH2(IK,1) = C2(IK,1)
   101 end do
   DO 103 J=2,IP
      DO 102 K=1,L1
         CH(1,K,J) = C1(1,K,J)
      102 end do
   103 end do
   IF (NBD > L1) GO TO 107
   IS = -IDO
   DO 106 J=2,IP
      IS = IS+IDO
      IDIJ = IS
      DO 105 I=3,IDO,2
         IDIJ = IDIJ+2
         DO 104 K=1,L1
            CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
            CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
         104 end do
      105 end do
   106 end do
   GO TO 111
   107 IS = -IDO
   DO 110 J=2,IP
      IS = IS+IDO
      DO 109 K=1,L1
         IDIJ = IS
         DO 108 I=3,IDO,2
            IDIJ = IDIJ+2
            CH(I-1,K,J) = WA(IDIJ-1)*C1(I-1,K,J)+WA(IDIJ)*C1(I,K,J)
            CH(I,K,J) = WA(IDIJ-1)*C1(I,K,J)-WA(IDIJ)*C1(I-1,K,J)
         108 end do
      109 end do
   110 end do
   111 IF (NBD < L1) GO TO 115
   DO 114 J=2,IPPH
      JC = IPP2-J
      DO 113 K=1,L1
         DO 112 I=3,IDO,2
            C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
            C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
            C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
            C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
         112 end do
      113 end do
   114 end do
   GO TO 121
   115 DO 118 J=2,IPPH
      JC = IPP2-J
      DO 117 I=3,IDO,2
         DO 116 K=1,L1
            C1(I-1,K,J) = CH(I-1,K,J)+CH(I-1,K,JC)
            C1(I-1,K,JC) = CH(I,K,J)-CH(I,K,JC)
            C1(I,K,J) = CH(I,K,J)+CH(I,K,JC)
            C1(I,K,JC) = CH(I-1,K,JC)-CH(I-1,K,J)
         116 end do
     117 end do
   118 end do
   GO TO 121
   119 DO 120 IK=1,IDL1
       C2(IK,1) = CH2(IK,1)
   120 end do
   121 DO 123 J=2,IPPH
      JC = IPP2-J
      DO 122 K=1,L1
         C1(1,K,J) = CH(1,K,J)+CH(1,K,JC)
         C1(1,K,JC) = CH(1,K,JC)-CH(1,K,J)
      122 end do
   123 end do

   AR1 = 1.
   AI1 = 0.
   DO 127 L=2,IPPH
      LC = IPP2-L
      AR1H = DCP*AR1-DSP*AI1
      AI1 = DCP*AI1+DSP*AR1
      AR1 = AR1H
      DO 124 IK=1,IDL1
         CH2(IK,L) = C2(IK,1)+AR1*C2(IK,2)
         CH2(IK,LC) = AI1*C2(IK,IP)
      124 end do
      DC2 = AR1
      DS2 = AI1
      AR2 = AR1
      AI2 = AI1
      DO 126 J=3,IPPH
         JC = IPP2-J
         AR2H = DC2*AR2-DS2*AI2
         AI2 = DC2*AI2+DS2*AR2
         AR2 = AR2H
         DO 125 IK=1,IDL1
            CH2(IK,L) = CH2(IK,L)+AR2*C2(IK,J)
            CH2(IK,LC) = CH2(IK,LC)+AI2*C2(IK,JC)
         125 end do
      126 end do
   127 end do
   DO 129 J=2,IPPH
      DO 128 IK=1,IDL1
         CH2(IK,1) = CH2(IK,1)+C2(IK,J)
      128 end do
   129 end do

   IF (IDO < L1) GO TO 132
   DO 131 K=1,L1
      DO 130 I=1,IDO
          CC(I,1,K) = CH(I,K,1)
      130 end do
   131 end do
   GO TO 135
   132 DO 134 I=1,IDO
      DO 133 K=1,L1
         CC(I,1,K) = CH(I,K,1)
      133 end do
   134 end do
   135 DO 137 J=2,IPPH
      JC = IPP2-J
      J2 = J+J
      DO 136 K=1,L1
         CC(IDO,J2-2,K) = CH(1,K,J)
         CC(1,J2-1,K) = CH(1,K,JC)
      136 end do
   137 end do
   IF (IDO == 1) RETURN
   IF (NBD < L1) GO TO 141
   DO 140 J=2,IPPH
      JC = IPP2-J
      J2 = J+J
      DO 139 K=1,L1
         DO 138 I=3,IDO,2
            IC = IDP2-I
            CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
            CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
            CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
            CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
         138 end do
      139 end do
   140 end do
   RETURN
   141 DO 144 J=2,IPPH
      JC = IPP2-J
      J2 = J+J
      DO 143 I=3,IDO,2
         IC = IDP2-I
         DO 142 K=1,L1
            CC(I-1,J2-1,K) = CH(I-1,K,J)+CH(I-1,K,JC)
            CC(IC-1,J2-2,K) = CH(I-1,K,J)-CH(I-1,K,JC)
            CC(I,J2-1,K) = CH(I,K,J)+CH(I,K,JC)
            CC(IC,J2-2,K) = CH(I,K,JC)-CH(I,K,J)
         142 end do
      143 end do
   144 end do
   RETURN
end SUBROUTINE RADFG

!     --- RFFTF1 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rfftf1 here]
SUBROUTINE RFFTF1 (N,C,CH,WA,IFAC)
   implicit _REAL_ (a-h,o-z)
   DIMENSION       CH(*)      ,C(*)       ,WA(*)      ,IFAC(*)
   NF = IFAC(2)
   NA = 1
   L2 = N
   IW = N
   DO 111 K1=1,NF
      KH = NF-K1
      IP = IFAC(KH+3)
      L1 = L2/IP
      IDO = N/L2
      IDL1 = IDO*L1
      IW = IW-(IP-1)*IDO
      NA = 1-NA
      IF (IP /= 4) GO TO 102
      IX2 = IW+IDO
      IX3 = IX2+IDO
      IF (NA /= 0) GO TO 101
      CALL RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
      GO TO 110
      101 CALL RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
      GO TO 110
      102 IF (IP /= 2) GO TO 104
      IF (NA /= 0) GO TO 103
      CALL RADF2 (IDO,L1,C,CH,WA(IW))
      GO TO 110
      103 CALL RADF2 (IDO,L1,CH,C,WA(IW))
      GO TO 110
      104 IF (IP /= 3) GO TO 106
      IX2 = IW+IDO
      IF (NA /= 0) GO TO 105
      CALL RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
      GO TO 110
      105 CALL RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
      GO TO 110
      106 IF (IP /= 5) GO TO 108
      IX2 = IW+IDO
      IX3 = IX2+IDO
      IX4 = IX3+IDO
      IF (NA /= 0) GO TO 107
      CALL RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
      GO TO 110
      107 CALL RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
      GO TO 110
      108 IF (IDO == 1) NA = 1-NA
      IF (NA /= 0) GO TO 109
      CALL RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
      NA = 1
      GO TO 110
      109 CALL RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
      NA = 0
      110 L2 = L1
   111 end do
   IF (NA == 1) RETURN
   DO 112 I=1,N
      C(I) = CH(I)
   112 end do
   RETURN
end SUBROUTINE RFFTF1

!     --- RFFTI ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rffti here]
SUBROUTINE RFFTI (N,WSAVE)
   implicit none
   integer N
   _REAL_       WSAVE(*)
   IF (N == 1) RETURN
   CALL RFFTI1 (N,WSAVE(N+1),WSAVE(2*N+1))
   RETURN
end SUBROUTINE RFFTI

!     --- RFFTI1 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rffti1 here]
SUBROUTINE RFFTI1 (N,WA,IFAC)
   implicit _REAL_ (a-h,o-z)
   DIMENSION       WA(*)      ,IFAC(*)    ,NTRYH(4)
   DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
   NL = N
   NF = 0
   J = 0
   101 J = J+1
   IF (J-4) 102,102,103
   102 NTRY = NTRYH(J)
   GO TO 104
   103 NTRY = NTRY+2
   104 NQ = NL/NTRY
   NR = NL-NTRY*NQ
   IF (NR) 101,105,101
   105 NF = NF+1
   IFAC(NF+2) = NTRY
   NL = NQ
   IF (NTRY /= 2) GO TO 107
   IF (NF == 1) GO TO 107
   DO 106 I=2,NF
      IB = NF-I+2
      IFAC(IB+2) = IFAC(IB+1)
   106 end do
   IFAC(3) = 2
   107 IF (NL /= 1) GO TO 104
   IFAC(1) = N
   IFAC(2) = NF
   TPI = 6.28318530717959
   ARGH = TPI/FLOAT(N)
   IS = 0
   NFM1 = NF-1
   L1 = 1
   IF (NFM1 == 0) RETURN
   DO 110 K1=1,NFM1
      IP = IFAC(K1+2)
      LD = 0
      L2 = L1*IP
      IDO = N/L2
      IPM = IP-1
      DO 109 J=1,IPM
         LD = LD+L1
         I = IS
         ARGLD = FLOAT(LD)*ARGH
         FI = 0.
         DO 108 II=3,IDO,2
            I = I+2
            FI = FI+1.
            ARG = FI*ARGLD
            WA(I-1) = COS(ARG)
            WA(I) = SIN(ARG)
         108 end do
         IS = IS+IDO
      109 end do
      L1 = L2
   110 end do
   RETURN
end SUBROUTINE RFFTI1

!     --- SINT ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine sint here]
SUBROUTINE SINT (N,X,WSAVE)
   implicit _REAL_ (a-h,o-z)
   DIMENSION       X(*)       ,WSAVE(*)
   NP1 = N+1
   IW1 = N/2+1
   IW2 = IW1+NP1
   IW3 = IW2+NP1
   CALL SINT1(N,X,WSAVE,WSAVE(IW1),WSAVE(IW2),WSAVE(IW3))
   RETURN
end SUBROUTINE SINT

!     --- SINT1 ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine sint1 here]
SUBROUTINE SINT1(N,WAR,WAS,XH,X,IFAC)
   implicit _REAL_ (a-h,o-z)
   DIMENSION WAR(*),WAS(*),X(*),XH(*),IFAC(*)
   DATA SQRT3 /1.73205080756888/
   DO 100 I=1,N
      XH(I) = WAR(I)
      WAR(I) = X(I)
   100 end do
   IF (N-2) 101,102,103
   101 XH(1) = XH(1)+XH(1)
   GO TO 106
   102 XHOLD = SQRT3*(XH(1)+XH(2))
   XH(2) = SQRT3*(XH(1)-XH(2))
   XH(1) = XHOLD
   GO TO 106
   103 NP1 = N+1
   NS2 = N/2
   X(1) = 0.
   DO 104 K=1,NS2
      KC = NP1-K
      T1 = XH(K)-XH(KC)
      T2 = WAS(K)*(XH(K)+XH(KC))
      X(K+1) = T1+T2
      X(KC+1) = T2-T1
   104 end do
   MODN = MOD(N,2)
   IF (MODN /= 0) X(NS2+2) = 4.*XH(NS2+1)
   CALL RFFTF1 (NP1,X,XH,WAR,IFAC)
   XH(1) = .5*X(1)
   DO 105 I=3,N,2
      XH(I-1) = -X(I)
      XH(I) = XH(I-2)+X(I-1)
   105 end do
   IF (MODN /= 0) GO TO 106
   XH(N) = -X(N+1)
   106 DO 107 I=1,N
      X(I) = WAR(I)
      WAR(I) = XH(I)
   107 end do
   RETURN
end SUBROUTINE SINT1

!     --- SINTI ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine sinti here]
SUBROUTINE SINTI (N,WSAVE)

!******************************************************************

!subroutine sinti initializes the array wsave which is used in
!subroutine sint. the prime factorization of n together with
!a tabulation of the trigonometric functions are computed and
!stored in wsave.

!input parameter
!n       the length of the sequence to be transformed.  the method
!        is most efficient when n+1 is a product of small primes.
!output parameter
!wsave   a work array with at least int(2.5*n+15) locations.
!        different wsave arrays are required for different values
!        of n. the contents of wsave must not be changed between
!        calls of sint

   implicit none
   integer n
   _REAL_  WSAVE(*)
   _REAL_ PI
   DATA PI /3.14159265358979/

   integer NS2, NP1, K
   _REAL_ DT

   IF (N <= 1) RETURN
   NS2 = N/2
   NP1 = N+1
   DT = PI/FLOAT(NP1)
   DO 101 K=1,NS2
      WSAVE(K) = 2.*SIN(K*DT)
   101 end do
   CALL RFFTI (NP1,WSAVE(NS2+1))
   RETURN
end SUBROUTINE SINTI
