#include "copyright.h"
#  define _REAL_ double precision
#include "timer.h"

#ifndef SANDER
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ SYSTEM CLOCK
subroutine wallclock( wallc )
   implicit none
   _REAL_ wallc
   integer ncalls,n
   integer mycount, myrate

   data ncalls /0/

   call system_clock( COUNT=mycount, COUNT_RATE=myrate)
   wallc = dble(mycount)/dble(myrate)
   ncalls = ncalls + 1
   return
entry nwallclock ( n )
   n = ncalls
   return
end subroutine wallclock
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PBSA TIMER
module pbtimer_module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Meng-Juei Hsieh
!  The Luo Research Group
!  University of California, Irvine
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none
   logical pbtimerswitch(pbsamaxtime)
   _REAL_ pbtime_bulletin (pbsamaxtime)
   character(len=25) PBTIME_desc (pbsamaxtime)

contains

subroutine pbtimer_init
   implicit none
   pbtimerswitch=.false.
   pbtime_bulletin=0d0

   PBTIME_desc(PBTIME_TOTAL)   ='Total time'
   PBTIME_desc(PBTIME_READ)    ='Read prm/crd time'
   PBTIME_desc(PBTIME_RUNMD)   ='Runmd Time'
   PBTIME_desc(PBTIME_FORCE)   ='Force time'
   PBTIME_desc(PBTIME_BOND)    ='Bond/Angle/Dihedral'
   PBTIME_desc(PBTIME_NONBON)  ='Nonbond force'
   !-- pb
   PBTIME_desc(PBTIME_PBFORCE)    = 'PB Nonbond'
   PBTIME_desc(PBTIME_PBLIST)     = 'PB NB list'
   PBTIME_desc(PBTIME_PBSETUP)    = 'PB FD grid'
   PBTIME_desc(PBTIME_PBSAS)      = 'PB Sasa'
   PBTIME_desc(PBTIME_PBFDFRC)    = 'PB FD force'
   PBTIME_desc(PBTIME_PBBUILDSYS) = 'PB Set linear sys'
   PBTIME_desc(PBTIME_PBSOLV)     = 'PB Solver'
   PBTIME_desc(PBTIME_PBITR)      = 'PB Iteration'
   PBTIME_desc(PBTIME_SINH)       = 'PB Sinh eva'
   PBTIME_desc(PBTIME_PBDBE)      = 'PB DB force'
   PBTIME_desc(PBTIME_PBMP)       = 'PB Multiple'
   PBTIME_desc(PBTIME_PBDIRECT)   = 'PB Direct'

   PBTIME_desc(PBTIME_PBSASRF)         = 'PB SA srf'

   PBTIME_desc(PBTIME_PBSAARC)         = 'PB SA arc'
   PBTIME_desc(PBTIME_PBSAARC_SETUP)   = 'PB SA arc setup'
   PBTIME_desc(PBTIME_PBCIRCLE)        = 'PB circle'
   PBTIME_desc(PBTIME_PBEXCLUDE)       = 'PB exclude'

   PBTIME_desc(PBTIME_PBEXMOL)         = 'PB exmol'
   PBTIME_desc(PBTIME_PBEXMOL_SETUP)   = 'PB exmol setup'
   PBTIME_desc(PBTIME_PBEXMOL_PARTA)   = 'PB exmol part a'
   PBTIME_desc(PBTIME_PBEXMOL_PARTB)   = 'PB exmol part b'
   PBTIME_desc(PBTIME_PBEXMOL_PARTC)   = 'PB exmol part c'
   PBTIME_desc(PBTIME_PBEXMOL_PARTD)   = 'PB exmol part d'
   PBTIME_desc(PBTIME_PBEXMOL_PARTE)   = 'PB exmol part e'
   PBTIME_desc(PBTIME_PBEXMOL_PARTF)   = 'PB exmol part f'
   PBTIME_desc(PBTIME_PBEPSBND)        = 'PB epsbnd'
   PBTIME_desc(PBTIME_PBEPSMAP)        = 'PB epsmap'

   PBTIME_desc(PBTIME_PBCALSA)         = 'PB calsa'
   !-- np
   PBTIME_desc(PBTIME_NPFORCE)    = 'NP Nonbond'
   PBTIME_desc(PBTIME_NPSAS)      = 'NP Sasa'
   PBTIME_desc(PBTIME_NPCAV)      = 'NP Cavity'
   PBTIME_desc(PBTIME_NPDIS)      = 'NP Dispersion'
   return
end subroutine pbtimer_init

subroutine pbtimer_start( num_event )
   implicit none

   integer num_event
   _REAL_ mytime

   call wallclock(mytime)
   if ( num_event > pbsamaxtime ) then
      write(6,*)'index ',num_event,' bigger than pbsamaxtime ',pbsamaxtime
      write(6,*)'attempt to add timer '
      call mexit(6,1)
   else if ( pbtimerswitch(num_event) ) then
      write(6,*)'timer: need to reset the time before start'
      call mexit(6,1)
   else if (pbtime_bulletin(num_event) <= 0) then
      pbtime_bulletin(num_event) = mytime
   else
      pbtime_bulletin(num_event) = mytime-pbtime_bulletin(num_event)
   endif
   pbtimerswitch(num_event) = .true.
   return
end subroutine pbtimer_start

subroutine pbtimer_stop( num_event )
   implicit none
   integer num_event
   _REAL_ mytime

   call wallclock(mytime)
   if ( num_event > pbsamaxtime )then
      write(6,*)'index ',num_event,' bigger than pbsamaxtime ',pbsamaxtime
      write(6,*)'attempt to close timer '
      call mexit(6,1)
   else if ( .not. pbtimerswitch(num_event) ) then
      write(6,*)'timer: need to start the time before stop'
      call mexit(6,1)
   end if
   pbtime_bulletin(num_event) = mytime-pbtime_bulletin(num_event)
   pbtimerswitch(num_event) = .false.
   return
end subroutine pbtimer_stop

subroutine pbtimer_summary
   implicit none
   integer i

#ifndef SANDER
   write(6,'(/80(1H-)/,''   5.  TIMINGS'',/80(1H-)/)')
#endif
   do i = 1, pbsamaxtime
      if ( .not. pbtimerswitch(i) ) then
         if (pbtime_bulletin(i) > 0) &
            write(6,"('|',1x,a,f10.2)") PBTIME_desc(i),pbtime_bulletin(i)
      else
         write(6,*) "Warning, this timer is not stopped: ",PBTIME_desc(i)
      endif
   enddo
   return
end subroutine pbtimer_summary


end module pbtimer_module
