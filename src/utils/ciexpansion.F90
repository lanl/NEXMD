!
!********************************************************************
!
!  Subroutine to print 2 electron repulsion integrals
!
!********************************************************************
!
!--------------------------------------------------------------------
!
!  Josiah A. Bjorgaard NDSU 2013
!
!--------------------------------------------------------------------
!
   subroutine ciexpansion()
   use qmmm_module
   use davidson_module
   implicit none
   
   do n=1:nocc
       do m=1:nvirt
       coeff0=qm2ds%v0
       enddo
   enddo
