#define _REAL_ double precision

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine prodis here]
subroutine prodis(l,m,n,nirreg,bcopt, bi,bo,atmfirst, atmlast, maxirr, cirreg, wp, qp)

   implicit none

   ! passed variables

   integer bcopt, atmfirst, atmlast, maxirr
   _REAL_ cirreg(maxirr, 15)
   _REAL_ wp(maxirr), qp(maxirr)
   _REAL_ t(3,3)
   _REAL_ bi,bo

   ! common variables

   integer l, m, n, nirreg
   !common /lmn/l, m, n, nirreg

   ! local variables

   integer ir
   _REAL_ xx, yy, zz

   ! external functions

   _REAL_ fw, fq

   do ir = 1, nirreg

      xx = cirreg(ir, 1)
      yy = cirreg(ir, 2)
      zz = cirreg(ir, 3)
   !if(ir == 1) write(10,*) 'xyz',xx,yy,zz
      t(1,1) = cirreg(ir, 7)
      t(1,2) = cirreg(ir, 8)
      t(1,3) = cirreg(ir, 9)
      t(2,1) = cirreg(ir, 10)
      t(2,2) = cirreg(ir, 11)
      t(2,3) = cirreg(ir, 12)
      t(3,1) = cirreg(ir, 13)
      t(3,2) = cirreg(ir, 14)
      t(3,3) = cirreg(ir, 15)
      wp(ir) = fw(bcopt,bi,bo,atmfirst,atmlast,xx,yy,zz)
      qp(ir) = fq(bcopt,bi,bo,atmfirst,atmlast,xx,yy,zz,t)

   end do


end subroutine prodis 
 
