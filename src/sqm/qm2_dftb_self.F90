! <compile=optimized>
!  -*- mode: f90; coding: iso-8859-15; -*-
#include "dprec.fh"

subroutine selfs(i,j,r2,iovpar,em,ne)
   implicit none

   integer, intent(in ) :: i,j
   _REAL_ , intent(in ) :: r2
   integer, intent(in ) :: ne
   _REAL_ , intent(out) :: em(ne,1)
   _REAL_ :: parm(13)
   integer :: id
   integer, external :: iovpar  

   id=iovpar(i,j,r2,parm)
   em(1,1)=parm(13)
   return
end subroutine selfs


subroutine selfp(i,j,r2,iovpar,em,ne)
   implicit none
   integer:: i, id, iovpar, j, l, m, ne
   _REAL_ :: em, parm, r2
   dimension parm(13),em(ne,3)
   external iovpar

   id=iovpar(i,j,r2,parm)
   do l=1,3
      do m=1,3
         em(l,m)=0.0
      enddo
      em(l,l)=parm(12)
   enddo
   return
end subroutine selfp


subroutine selfd(i,j,r2,iovpar,em,ne)
   implicit none
   integer:: i, id, iovpar, j, l, m, ne
   _REAL_ :: em, parm, r2
   dimension parm(13),em(ne,5)
   external iovpar

   id=iovpar(i,j,r2,parm)
   do l=1,5
      do m=1,5
         em(l,m)=0.0
      enddo
      em(l,l)=parm(11)
   enddo
   return
end subroutine selfd
