#define _REAL_ double precision
      FUNCTION pythag(a,b)
      implicit none
#     include "pb_constants.h"
      _REAL_ a,b,pythag
      _REAL_ absa,absb
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
        pythag=absa*sqrt(ONE+(absb/absa)**2)
      else
        if(absb.eq.ZERO)then
          pythag=ZERO
        else
          pythag=absb*sqrt(ONE+(absa/absb)**2)
        endif
      endif
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.
