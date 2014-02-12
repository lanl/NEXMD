#include "dprec.fh"
c-----------------------------------------------------------------------
      subroutine amflsh(lun)
c     wrapper for i/o buffer flushing routine
c     Author: George Seibel
c
      implicit none
c     INPUT:
c
      integer lun
c        ... logical unit number to flush
#ifdef SGI
      integer istat
#endif
c        ... return status from flush
c
c     --- most Unix (BSD, Convex, Sun, Stellar, SGI Iris...) ---
c
#if defined(AIX) || defined( XLF90)
c     --new for Version 2.3 of XLF; page 222 in the V2.3 Language Reference:
      CALL FLUSH_(LUN)
#else
#  ifdef NO_FLUSH_STDOUT
      if( lun /= 6 ) call flush(lun)
#  else
#     ifdef SGI
         call flush(lun,istat)
#     else
         call flush(lun)
#     endif
#  endif
#endif
      return
      end
c-----------------------------------------------------------------------

#ifdef F90_TIMER
      subroutine wallclock( wallc )
      implicit none
      _REAL_ wallc

      integer ncalls,n
      data ncalls /0/

#  ifdef NO_DETAILED_TIMINGS
      wallc = 0.d0
#  else
      integer count, rate
      call system_clock( COUNT=count, COUNT_RATE=rate)
      wallc = dble(count)/dble(rate)
      ncalls = ncalls + 1
#  endif
      return
c
      entry nwallclock ( n )
      n = ncalls
      return
c
      end
#endif
