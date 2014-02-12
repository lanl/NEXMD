#include "copyright.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmlsrc here]
subroutine nmlsrc(name,iun,ifind)

   implicit none
   integer:: i, ifind, ilen, iun, j
   ! Subroutine NaMeList SeaRCh

   ! This routine searches the file on unit IUN for the namelist with name
   ! NAME. The namelist is denoted by a line where the first non-blank character
   ! string is either &NAME or $NAME.

   ! IFIND: 0 Namelist was not found. File will be rewound before return.
   !        1 Namelist was found. File will be backspaced to the line that
   !          contained the namelist flag before returning.

   ! Limitations: the namelist flat must fit entirely in the first 80 columns
   !              of the line on which it appears.

   ! Author: David A. Pearlman
   ! Date: 1/91

   character(len=*) name
   character(len=80) aline
   ilen = len(name)

   do i = 1,999
      read(iun,1000,end=500,err=500) aline
      do j = 1,80
         if (aline(j:j) /= ' ') then
            if (aline(j:j) == '&' .or. aline(j:j) == '$') then
               if (80-j+1 >= ilen) then
                  ! This includes a trailing space to match the word length.
                  if (aline(j+1:j+ilen+1) == name) goto 200
               else
                  exit
               end if
            else
               exit
            end if
         end if
      end do
   end do

   ! Namelist flag not found:

   500 ifind = 0
   rewind(iun)
   return

   ! Namelist flag found:

   200 ifind = 1
   backspace(iun)
   return


   ! Format statements:

   1000 format(a)
end subroutine nmlsrc 
