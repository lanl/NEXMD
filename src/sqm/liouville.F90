#include "dprec.fh"
!
!********************************************************************
!
!  Subroutines to work in Liouville space
!
!  Converted to fortran 90 and interfaced with sqm(AmberTools) by
!  Kirill Velizhanin (kirill@lanl.gov)
!
!********************************************************************
!
!--------------------------------------------------------------------
!
!  Acting by TDHF tetradic operator L on an arbitrary single-particle
!  matrix xi.
!       
!  Notation closely resembles one in S. Tretiak and S. Mukamel, Chem. Rev.,
!  102, 3171-3212 (2002)
!
!--------------------------------------------------------------------
!
   subroutine site2mo(zz,xi,v1)
   use qm2_davidson_module
   implicit none
   _REAL_ f0,f1
   parameter (f0 = 0.0)
   parameter (f1 = 1.0)
   _REAL_ xi(qm2ds%Nb,qm2ds%Nb),v1(qm2ds%Nrpa)
   _REAL_ zz(qm2ds%Nb,qm2ds%Nb)
   call dgemm ('T','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,f1,xi,qm2ds%Nb, &
      qm2ds%vhf,qm2ds%Nb,f0,zz,qm2ds%Nb)
   call dgemm ('T','N',qm2ds%Nh,qm2ds%Np,qm2ds%Nb,f1, &
      qm2ds%vhf(:,qm2ds%Np+1:qm2ds%Nb),qm2ds%Nb,zz,qm2ds%Nb,f0,v1,qm2ds%Nh)
   call dgemm ('T','N',qm2ds%Nh,qm2ds%Np,qm2ds%Nb,f1, &
      zz(1,qm2ds%Np+1),qm2ds%Nb,qm2ds%vhf,qm2ds%Nb,f0, &
      v1(qm2ds%Ncis+1),qm2ds%Nh)      
   return
   end subroutine
!
   subroutine mo2site (v1,xi,zz)
   use qm2_davidson_module
   implicit none
   _REAL_ f0,f1
   parameter (f0 = 0.0)
   parameter (f1 = 1.0)
   _REAL_ xi(qm2ds%Nb,qm2ds%Nb),v1(qm2ds%Nrpa)
   _REAL_ zz(qm2ds%Nb,qm2ds%Nb)
   call dgemm ('T','T',qm2ds%Np,qm2ds%Nb,qm2ds%Nh,f1,v1(1),qm2ds%Nh, &
      qm2ds%vhf(:,qm2ds%Np+1:qm2ds%Nb),qm2ds%Nb,f0,zz,qm2ds%Nb)
   call dgemm ('N','T',qm2ds%Nh,qm2ds%Nb,qm2ds%Np,f1,v1(qm2ds%Ncis+1), &
      qm2ds%Nh,qm2ds%vhf,qm2ds%Nb,f0,zz(qm2ds%Np+1,1),qm2ds%Nb)
   call dgemm ('N','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,f1,qm2ds%vhf, &
      qm2ds%Nb,zz,qm2ds%Nb,f0,xi,qm2ds%Nb)
   return
   end subroutine
!
!********************************************************************
!
   subroutine site2moph (Nb,Np,Nh,M4,uu,xi,v1,zz)
   implicit none
   integer Nb,Np,Nh,M4
   _REAL_ uu(Nb*Nb),xi(Nb,Nb),zz(Nb,Nb),v1(2*M4)
   _REAL_ f0,f1
   parameter (f0 = 0.0)
   parameter (f1 = 1.0)
   call dgemm ('T','N',Nb,Nb,Nb,f1,xi,Nb,uu,Nb,f0,zz,Nb)
   call dgemm ('T','N',Nh,Np,Nb,f1,uu(Nb*Np+1),Nb,zz,Nb,f0,v1,Nh)
   call dgemm ('T','N',Nh,Np,Nb,f1,zz(1,Np+1),Nb,uu,Nb,f0,v1(M4+1),Nh)      
   return
   end subroutine
        
   subroutine mo2siteph (Nb,Np,Nh,M4,uu,v1,xi,zz)
   implicit none
   integer Nb,Np,Nh,M4
   _REAL_ uu(Nb*Nb),xi(Nb,Nb),zz(Nb,Nb),v1(2*M4)
   _REAL_ f0,f1
   call dgemm ('T','T',Np,Nb,Nh,f1,v1(1),Nh,uu(Nb*Np+1),Nb,f0,zz,Nb)
   call dgemm ('N','T',Nh,Nb,Np,f1,v1(M4+1),Nh,uu,Nb,f0,zz(Np+1,1),Nb)
   call dgemm ('N','N',Nb,Nb,Nb,f1,uu,Nb,zz,Nb,f0,xi,Nb)
   return
   end subroutine
!
!********************************************************************
!
   subroutine site2mof (Nb,uu,xi,eta,zz)
   implicit none
   integer Nb
   _REAL_ uu(Nb,Nb),xi(Nb,Nb),zz(Nb,Nb),eta(Nb,Nb)
   _REAL_ f0,f1
   parameter (f0 = 0.0)
   parameter (f1 = 1.0)
   call dgemm ('N','N',Nb,Nb,Nb,f1,xi,Nb,uu,Nb,f0,zz,Nb)
   call dgemm ('T','N',Nb,Nb,Nb,f1,uu,Nb,zz,Nb,f0,eta,Nb)
   return
   entry mo2sitef (Nb,uu,xi,eta,zz)
   call dgemm ('N','T',Nb,Nb,Nb,f1,xi,Nb,uu,Nb,f0,zz,Nb)
   call dgemm ('N','N',Nb,Nb,Nb,f1,uu,Nb,zz,Nb,f0,eta,Nb)
   return
   end
!
!********************************************************************
!
!********************************************************************
!
   subroutine Vxi(xi,eta)
   use qm2_davidson_module
   implicit none
   _REAL_ xi(qm2ds%Nb,qm2ds%Nb),eta(qm2ds%Nb,qm2ds%Nb)
   _REAL_ coef
   integer i,j,l
   qm2ds%xis(:)=0.d0
   qm2ds%etas(:)=0.d0
   coef = 1.000000;
! symmetric part:
   l=0
   do i=1,qm2ds%Nb
    do j=1,i
       l=l+1
       qm2ds%xis(l)=coef*0.5*(xi(i,j)+xi(j,i))
    end do
   end do

   call Vxi_pack(qm2ds%xis,qm2ds%etas)
    l=0
    do i=1,qm2ds%Nb
      do j=1,i-1
         l = l + 1
         eta(i,j)=eta(i,j)+qm2ds%etas(l)
         eta(j,i)=eta(j,i)+qm2ds%etas(l)
      end do
      l=l+1
      eta(i,i)=qm2ds%etas(l)
   end do

!  antisymmetric part:
   l=0
   do i=1,qm2ds%Nb
    do j = 1,i
      l = l + 1
      qm2ds%xis(l)=coef*0.5*(xi(i,j)-xi(j,i))
    end do
   end do

   call Vxi_packA(qm2ds%xis,qm2ds%etas)
   l=0
   do i=1,qm2ds%Nb
      do j = 1,i-1
         l = l + 1
         eta(i,j) = eta(i,j) + qm2ds%etas(l)
         eta(j,i) = eta(j,i) - qm2ds%etas(l)
      end do

      l = l + 1
!     if (abs(etas(l)).gt.1d-10) stop 'Vxi: design'
   end do
   return

   entry Vxi_symm (xi,eta)   ! eta = Vxi
!  pack:
   l = 0
   do i = 1,qm2ds%Nb
    do j = 1,i
       l = l + 1
       qm2ds%xis(l) = xi(j,i)
    end do
   end do

!  multiply:
   call Vxi_pack(qm2ds%xis,qm2ds%etas)
!  unpack:
   l = 0
   do i = 1,qm2ds%Nb
    do j = 1,i
       l = l + 1
       eta(i,j) = qm2ds%etas(l)
       eta(j,i) = qm2ds%etas(l)
    end do
   end do
   end
!
!********************************************************************
!
!********************************************************************
!
   subroutine Vxi_pack(xi,eta)
   use qmmm_module,only: qm2_params
   use qm2_davidson_module

   implicit none
   _REAL_ xi(qm2ds%Lt),eta(qm2ds%Lt)
   character keywr*6
   common /keywr/ keywr
   logical first
      data first /.true./
      save first

! --- if this is the first time in this routine, load coulomb matrix
      if (first) then
   if (index(keywr,'INDO').NE.0.AND.index(keywr,'MINDO').EQ.0) then
     write(6,*)  keywr,' hamiltonian requested'
     write(6,*)  'Use *Z program for ', keywr
     stop
   endif 
         first=.false.
      endif
      eta(:)=0.0
      call qm2_fock2(eta,xi,qm2ds%W,qm2_params%orb_loc)
   if (qm2ds%iderivfl.eq.0) then ! We are not in analytic derivatives     
      call qm2_fock1(eta,xi) 
      endif
      eta(:)=eta(:)*2.0 !Why *2.0? Is it for the commutator in L(xi)?
   return
   end
!
!********************************************************************
!
!********************************************************************
!
   subroutine Vxi_packA (xi,eta)
   use qmmm_module,only: qm2_params
   use qm2_davidson_module
   implicit none
   _REAL_ xi(qm2ds%Lt),eta(qm2ds%Lt)
   character keywr*6
   common /keywr/ keywr
   logical first
      data first /.true./
      save first

! --- if this is the first time in this routine, load coulomb matrix
      if (first) then
   if (index(keywr,'INDO').NE.0.AND.index(keywr,'MINDO').EQ.0) then
     write(6,*)  keywr,' hamiltonian requested'
     write(6,*)  'Use *Z program for ', keywr
     stop
   endif 
         first=.false.
      endif
      eta(:)=0.0
      call qm2_fock2(eta,xi,qm2ds%W,qm2_params%orb_loc)
   if (qm2ds%iderivfl.eq.0) then ! We are not in analytic derivatives
      call qm2_fock1_skew(eta,xi)
      endif
      eta(:)=eta(:)*2.0
   return
   end
!

