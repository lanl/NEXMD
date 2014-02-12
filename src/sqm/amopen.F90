#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine amopen here]
subroutine amopen(lun,fname,fstat,fform,facc)
   
   !  When this is converted to Fortran 90
   !  these codes can be replaced with module types, eg,
   !       Character(*), public, parameter :: unknown = 'unknown'
   !  facc is used for file appending

   implicit none
   
   !     INPUT:
   
   integer lun
   !        ... logical unit number
   character(len=*) fname
   !        ... file name 
   character(len=1) fstat
   !        ... status code: 'N', 'O', 'R', 'U' = new, old, replace, unknown
   character(len=1) fform
   !        ... format code: 'F', 'U' = formatted, unformatted
   character(len=1) facc
   !        ... access code: 'R', 'W', 'A' = read, read/write, append
   
   !     INTERNAL:

   character(len=7) stat
   !        ... status keyword
   character(len=11) kform
   !        ... form keyword
   character(len=11) pos
   !        ... position keyword
   integer ios 
   !        ... i/o status variable

   if (fstat == 'N') then
      stat = 'NEW'
   else if (fstat == 'O') then
      stat = 'OLD'
   else if (fstat == 'R') then
      stat = 'REPLACE'
   else if (fstat == 'U') then
      stat = 'UNKNOWN'
   else
      write(6,'(/,2x,a,i4)') &
            'amopen: bogus fstat, unit ', lun
      call mexit(6, 1)
   end if
   
   if (fform == 'U') then
      kform = 'UNFORMATTED'
   else if (fform == 'F') then
      kform = 'FORMATTED'
   else
      write(6,'(/,2x,a,i4)') &
            'amopen: bogus fform, unit', lun
      call mexit(6, 1)
   end if

   if (facc == 'A') then
      pos = "APPEND"
   else
      pos = "ASIS" ! default f90
   end if
   
   open(unit=lun,file=fname,status=stat,form=kform,iostat=ios,position=pos)
   
   if (ios /= 0) then
      if (lun == 6) then
#ifndef DUMB
         write(0,'(/,2x,a,i4,a,a)') 'Unit ', lun, &
               ' Error on OPEN: ',fname
#endif
      else
         write(6,'(/,2x,a,i4,a,a)') 'Unit ', lun, &
               ' Error on OPEN: ',fname
         close(unit=6)
#ifndef DUMB
         write(0,'(/,2x,a,i4,a,a)') 'Unit ', lun, &
               ' Error on OPEN: ',fname
#endif
      end if
      call mexit(6, 1)
   end if
   if (pos /= "APPEND") rewind(lun)
   return
end subroutine amopen 
