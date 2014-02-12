#include "copyright.h"
#include "dprec.fh"
subroutine qm2_smallest_number(SMALL,SMALLSUM)

! Calculates the smallest representable number SMALL and the
! smallest number SMALLSUM for which 1+SMALLSUM != 1

! Written by Ross Walker (TSRI, 2005)

      implicit none

      _REAL_, intent(out) :: SMALL, SMALLSUM

      SMALL = 1.0D0
      do while ((SMALL*0.5D0) /= 0.0D0)
        SMALL = SMALL * 0.5D0
      end do
      SMALLSUM=1.0D0
      do while ((1.0D0+(SMALLSUM*0.5D0)) /= 1.0D0)
        SMALLSUM=SMALLSUM*0.5D0
      end do
      return
end subroutine qm2_smallest_number

