#define _REAL_ double precision
      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      implicit none
#     include "pb_constants.h"
      INTEGER m,mp,n,np,NMAX
      _REAL_ b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=500)
      INTEGER i,j,jj
      _REAL_ s,tmp(NMAX)
      do 12 j=1,n
        s=0.
        if(w(j).ne.ZERO)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
