#include "dprec.fh"
#include "assert.fh"

module quantum_prop_add
use naesmd_constants
use communism

contains

   subroutine checkcrossing(sim,Na,Nm,cross,lprint)
   implicit none

   type(simulation_t),pointer::sim
   integer k,j,i,ii,iii,iimdqt,lprint
! modified by Seba################
!      integer cross
! end modified by Seba################
   integer mdflag,Na,Nm,imax
   include 'md.par'
   include 'parH.par'
   include 'sizes'
   include 'md.cmn'
   include 'common'
! added by Seba################
   integer cross(sim%excN)
! end added by Seba################
   integer ascpr(260,260),z
   _REAL_ scprold,xmax
   _REAL_ normcheck(sim%excN),normcheckhop 
!      double precision normcheckhop 

!--------------------------------------------------------------------
!
!  Following the nonavoiding crossing of states
!
!  scpr is the overlap between new and old transition densities
!  if there is no crossing, the matrix is supposed to be close
!  to identity matrix
!
!--------------------------------------------------------------------
!
   scpr(1:sim%excN,1:sim%excN)=0.d0

   ! kav: FIXME
   ! Here and below, in the old code the size Ncis is used as a size of
   ! cmdqt vectors. This is not entirely correct, since the size is 
   ! Nrpa=2*Ncis (check it!!!). However, it can be very accurate becuase
   ! the contribution of negative energies to a positive excitation
   ! energy can be very small (would be great to check it)
   do i=1,sim%excN
      do j=1,sim%excN
         do k=1,sim%dav%Ncis
            scpr(i,j)=scpr(i,j)+cmdqtold(k,i)*cmdqtnew(k,j)
         end do
      end do
   end do

! modified by Seba######################################
! this part of the code that was previously commented, I have discommented it.
   do i=1,sim%excN
      do j=1,sim%excN
         ascpr(i,j)=int(scpr(i,j)**2*1.d5)
      end do
   end do
! window**************************************
   do i=1,sim%excN
      do j=1,sim%excN
         if((j.lt.(i-2)).or.(j.gt.(i+2))) then
            ascpr(i,j)=-1*ifix(sngl(1.d5))
         end if
      end do
   end do

!************************************

   do i=1,sim%excN
      do j=1,sim%excN
         ascpr(i,j)=-1*ascpr(i,j)
      enddo
   end do

   call apc(sim%excN,ascpr,iorden,z)

! end modified by Seba##############################

! modified by Seba#################################
!
!      xmax=0.0d0
!      xmax=0.0d0
!      imax=ihop
!      do i=1,sim%excN
!         if(dabs(scpr(ihop,i)).gt.xmax) then
!           xmax=dabs(scpr(ihop,i))
!           imax=i
!         endif
!      enddo
!      iorden(ihop)=imax
!
!      if(iorden(ihop).ne.ihop) then 
!!         normcheckhop=0.0d0
!!         do j=1,sim%excN
!!            normcheckhop=normcheckhop+scpr(j,iorden(ihop))**2
!!         enddo
!         if(dabs(scpr(ihop,iorden(ihop))).lt.0.9d0) then
!           cross=1
!         else
!           cross=2
!         endif
!         iordenhop=iorden(ihop)
!         if(lprint.ge.2) then
!            write(101,687) tfemto,cross,ihop,iordenhop, &
!      scpr(ihop,iorden(ihop))
!            call flush(101)
!         endif
!      else
!         cross=0
!      endif

   do i=1,sim%excN
      if(iorden(i).ne.i) then
         if(i.lt.iorden(i).or.i.eq.ihop) then
            if(dabs(scpr(i,iorden(i))).lt.0.9d0) then
               cross(i)=1
            else
               cross(i)=2
            end if

            iordenhop(i)=iorden(i)
            if(lprint.ge.2) then
!                 write(101,687) tfemto,cross(i),i,iordenhop(i),
!     $scpr(i,iorden(i))
!                 call flush(101)
            endif
         else
            cross(i)=0
         end if
      else
         cross(i)=0
      end if
   end do

! end modified by Seba#########################################

! Store the overlap matrix to print it if hops or cross

   do i=1,sim%excN
      do j=1,sim%excN
         scprreal(i,j)=scpr(i,j)
      end do
   end do

!***********************************

687   format(F18.10,1X,I4,1X,I4,1X,I4,1000(1x,f18.10))
688   FORMAT(f18.10,10000(1X,I4))
889   format(1000(1x,f18.10))
886   FORMAT(10000(1X,F7.4))

   return
   end subroutine

   subroutine checkcrossingmiddle(sim,Na,Nm,cross)

   implicit none

   type(simulation_t),pointer::sim
! modified by Seba
!      integer k,j,i,ii,iii,iimdqt,cross
   integer k,j,i,ii,iii,iimdqt
! end modified by Seba
   integer mdflag,Na,Nm
   include 'md.par'
   include 'parH.par'
   include 'sizes'
   include 'md.cmn'
   include 'common'
! modified by Seba
   integer cross(sim%excN)
! end modified by Seba
   integer ascpr(260,260),z
   double precision scprold
   double precision normcheck(sim%excN),normcheckhop 
!      double precision normcheckhop 


!***************************************************
! following the nonavoiding crossing of states
!***************************************************
   do i=1,sim%excN
      do j=1,sim%excN
         scpr(i,j)=0.0d0
         do k=1,sim%dav%Ncis
            scpr(i,j)=scpr(i,j)+cmdqtmiddleold(k,i)*cmdqtmiddle(k,j)
         end do
      end do
   end do
! modified by Seba
!      if(dabs(scpr(ihop,iordenhop)).ge.0.9d0) then
!         cross=2
!      endif
   do i=1,sim%excN
      if(i.lt.iorden(i).or.i.eq.ihop) then
         if(i.ne.iorden(ihop)) then
            if(dabs(scpr(i,iorden(i))).ge.0.9d0) then
               cross(i)=2
            end if
         end if
      end if
   end do
! end modified by Seba




!***********************************

687   format(F18.10,1X,I4,1X,I4,1X,I4,1000(1x,f18.10))
688   FORMAT(f18.10,10000(1X,I4))
889   format(1000(1x,f18.10))
886   FORMAT(10000(1X,F7.4))
686   FORMAT(10000(F18.10,1X,F7.4))

   return
   end subroutine 

! modified by Seba
!      SUBROUTINE vmdqtmiddlecalc(sim, iimdqt,Na,Nm,cross)
   subroutine vmdqtmiddlecalc(sim, iimdqt,Na,Nm)
! end modified by Seba

   implicit none

   type(simulation_t), pointer :: sim
! modified by Seba
!      integer k,j,i,ii,iii,iimdqt,cross
   integer k,j,i,ii,iii,iimdqt
! end modified by Seba
   integer mdflag,Na,Nm
   double precision E0 
   include 'md.par'
   include 'parH.par'
   real*8 Omega(Mx_M),fosc(Mx_M)
   real*8 xx(Na_M),yy(Na_M),zz(Na_M)
   real*8 xxp(Na_M),yyp(Na_M),zzp(Na_M)
   real*8 xxm(Na_M),yym(Na_M),zzm(Na_M)
   include 'sizes'
   include 'md.cmn'
   include 'common'

   if(iimdqt.eq.1) then
      do i=1,sim%excN
         vmdqtmiddleold(i)=vmdqtold(i) 
      end do

      do i=1,sim%dav%Ncis
         do j=1,sim%excN
            cmdqtmiddleold(i,j)=cmdqtold(i,j)
         end do
      end do

   else
      do i=1,sim%excN
         vmdqtmiddleold(i)=vmdqtmiddle(i) 
      end do

      do i=1,sim%dav%Ncis
         do j=1,sim%excN
            cmdqtmiddleold(i,j)=cmdqtmiddle(i,j)
         end do
      end do
   end if

   if(iimdqt.eq.nquantumstep) then

      do i=1,sim%excN
         vmdqtmiddle(i)=vmdqtnew(i)
      end do

      do i=1,sim%dav%Ncis
         do j=1,sim%excN
            cmdqtmiddle(i,j)=cmdqtnew(i,j)
         end do
      end do

   else

      if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
         do j=1,natom
            xx(j)=rxold(j) &
               +vxold(j)*dtquantum*dfloat(iimdqt) &
               +axold(j)*0.5d0*(dtquantum*dfloat(iimdqt))**2

            yy(j)=ryold(j) &
               +vyold(j)*dtquantum*dfloat(iimdqt) &
               +ayold(j)*0.5d0*(dtquantum*dfloat(iimdqt))**2

            zz(j)=rzold(j) &
               +vzold(j)*dtquantum*dfloat(iimdqt) &
               +azold(j)*0.5d0*(dtquantum*dfloat(iimdqt))**2
         end do
      end if

      if(ensemble.eq.'langev') then
         do j=1,natom

            xx(j)=rxold(j) &
               +vxold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt) &
               +axold(j)*afric(j)/(dtmdqt*dtmdqt) &
               *(dtquantum*dfloat(iimdqt))**2 &
               +prand(1,j)/dtmdqt*dtquantum*dfloat(iimdqt)

            yy(j)=ryold(j) &
               +vyold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt) &
               +ayold(j)*afric(j)/(dtmdqt*dtmdqt) &
               *(dtquantum*dfloat(iimdqt))**2 &
               +prand(2,j)/dtmdqt*dtquantum*dfloat(iimdqt)

            zz(j)=rzold(j) &
               +vzold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt) &
               +azold(j)*afric(j)/(dtmdqt*dtmdqt) &
               *(dtquantum*dfloat(iimdqt))**2 &
               +prand(3,j)/dtmdqt*dtquantum*dfloat(iimdqt)
         end do
      end if


      call do_sqm_davidson_update(sim,cmdqt=cmdqtmiddle, &
         vmdqt=vmdqtmiddle,vgs=vgs,rx=xx,ry=yy,rz=zz)        
!          !KGB
!          ! call ceo(Na,xx,yy,zz,atoms,sim%excN,E0,Omega,fosc,mdflag)
!          sim%dav%mdflag=2
!          call do_sqm_and_davidson(sim, xx, yy, zz)
!          call dav2naesmd_Omega(sim)
! 
!          vgs=E0
!          do i=1,sim%excN
!             vmdqtmiddle(i)=(Omega(i)+E0) 
!          enddo
! 
!          call dav2cmdqt(sim, cmdqtmiddle)
!          do i=1,nbasis
!             do j=1,sim%excN
!                cmdqtmiddle(i,j)=cmdqt(i,j)
!             enddo
!          enddo


   end if

! modified by Seba
!      if(lowvalue.gt.dabs(vmdqtmiddle(ihop)- &
!      vmdqtmiddle(iordenhop))) then
!         lowvalue=dabs(vmdqtmiddle(ihop)- &
!      vmdqtmiddle(iordenhop))
!         lowvaluestep=iimdqt
!      endif

   do j=1,sim%excN
      if(lowvalue(j).gt.dabs(vmdqtmiddle(j)- &
         vmdqtmiddle(iordenhop(j)))) then

         lowvalue(j)=dabs(vmdqtmiddle(j)- &
            vmdqtmiddle(iordenhop(j)))
         lowvaluestep(j)=iimdqt
      end if
   end do

! end modified by Seba

!      write(103,687) tfemtoquantum, vmdqtmiddle(ihop)
!     $-vmdqtmiddle(iordenhop)
!      call flush(103)

889   FORMAT(10000(1X,F18.10))
687   format(1000(1x,f18.10))

   return
   end subroutine

! modified by Seba
!      SUBROUTINE vmdqtlowvalue(sim,win,iimdqt,Na,Nm,cross)
   subroutine vmdqtlowvalue(sim,win,iimdqt,Na,Nm)
! end modified by Seba

   implicit none

   type(simulation_t), pointer :: sim
! modified by Seba
!      integer k,j,i,ii,iii,iimdqt,cross,win
   integer k,j,i,ii,iii,iimdqt,win
! end modified by Seba
   integer mdflag,Na,Nm
   double precision E0 
   include 'md.par'
   include 'parH.par'
   real*8 Omega(Mx_M),fosc(Mx_M)
   real*8 xx(Na_M),yy(Na_M),zz(Na_M)
   real*8 xxp(Na_M),yyp(Na_M),zzp(Na_M)
   real*8 xxm(Na_M),yym(Na_M),zzm(Na_M)
   include 'sizes'
   include 'md.cmn'
   include 'common'

   if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
      do j=1,natom
         xx(j)=rxold(j) &
            +vxold(j)*dtquantum*dfloat(iimdqt+win) &
            +axold(j)*0.5d0*(dtquantum*dfloat(iimdqt+win))**2

         yy(j)=ryold(j) &
            +vyold(j)*dtquantum*dfloat(iimdqt+win) &
            +ayold(j)*0.5d0*(dtquantum*dfloat(iimdqt+win))**2

         zz(j)=rzold(j) &
            +vzold(j)*dtquantum*dfloat(iimdqt+win) &
            +azold(j)*0.5d0*(dtquantum*dfloat(iimdqt+win))**2
      end do
   end if

   if(ensemble.eq.'langev') then
      do j=1,natom
         xx(j)=rxold(j) &
            +vxold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt+win) &
            +axold(j)*afric(j)/(dtmdqt*dtmdqt) &
            *(dtquantum*dfloat(iimdqt+win))**2 &
            +prand(1,j)/dtmdqt*dtquantum*dfloat(iimdqt+win)

         yy(j)=ryold(j) &
            +vyold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt+win) &
            +ayold(j)*afric(j)/(dtmdqt*dtmdqt) &
            *(dtquantum*dfloat(iimdqt+win))**2 &
            +prand(2,j)/dtmdqt*dtquantum*dfloat(iimdqt+win)

         zz(j)=rzold(j) &
            +vzold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt+win) &
            +azold(j)*afric(j)/(dtmdqt*dtmdqt) &
            *(dtquantum*dfloat(iimdqt+win))**2 &
            +prand(3,j)/dtmdqt*dtquantum*dfloat(iimdqt+win)
      end do
   end if
         


   call do_sqm_davidson_update(sim,cmdqt=cmdqtmiddle, &
      vmdqt=vmdqtmiddle,vgs=vgs,rx=xx,ry=yy,rz=zz)        

         ! KGB
         ! call ceo(Na,xx,yy,zz,atoms,sim%excN,E0,Omega,fosc,mdflag)
!          sim%dav%mdflag=2
!          call do_sqm_and_davidson(sim, xx, yy, zz)
!          call dav2naesmd_Omega(sim)
!          call dav2cmdqt(sim, cmdqtmiddle)

!          do i=1,nbasis
!             do j=1,sim%excN
!                cmdqtmiddle(i,j)=cmdqt(i,j) !FIXME
!             enddo
!          enddo


   if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
      do j=1,natom
         xx(j)=rxold(j) &
            +vxold(j)*dtquantum*dfloat(iimdqt-win) &
            +axold(j)*0.5d0*(dtquantum*dfloat(iimdqt-win))**2

         yy(j)=ryold(j) &
            +vyold(j)*dtquantum*dfloat(iimdqt-win) &
            +ayold(j)*0.5d0*(dtquantum*dfloat(iimdqt-win))**2

         zz(j)=rzold(j) &
            +vzold(j)*dtquantum*dfloat(iimdqt-win) &
            +azold(j)*0.5d0*(dtquantum*dfloat(iimdqt-win))**2
      end do
   end if

   if(ensemble.eq.'langev') then
      do j=1,natom
         xx(j)=rxold(j) &
            +vxold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt-win) &
            +axold(j)*afric(j)/(dtmdqt*dtmdqt) &
            *(dtquantum*dfloat(iimdqt-win))**2 &
            +prand(1,j)/dtmdqt*dtquantum*dfloat(iimdqt-win)

         yy(j)=ryold(j) &
            +vyold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt-win) &
            +ayold(j)*afric(j)/(dtmdqt*dtmdqt) &
            *(dtquantum*dfloat(iimdqt-win))**2 &
            +prand(2,j)/dtmdqt*dtquantum*dfloat(iimdqt-win)

         zz(j)=rzold(j) &
            +vzold(j)*vfric(j)/dtmdqt*dtquantum*dfloat(iimdqt-win) &
            +azold(j)*afric(j)/(dtmdqt*dtmdqt) &
            *(dtquantum*dfloat(iimdqt-win))**2 &
            +prand(3,j)/dtmdqt*dtquantum*dfloat(iimdqt-win)
      end do
   end if

!          do j=1,natom
!             xx(j)=xx(j)*convl
!             yy(j)=yy(j)*convl
!             zz(j)=zz(j)*convl
!          enddo

         ! call ceo(Na,xx,yy,zz,atoms,sim%excN,E0,Omega,fosc,mdflag)
   
   call do_sqm_davidson_update(sim,cmdqtmiddleold,rx=xx,ry=yy,rz=zz)

         ! FIXME
!          do i=1,nbasis
!             do j=1,sim%excN
!                cmdqtmiddleold(i,j)=cmdqt(i,j)
!             enddo
!          enddo
! 
!          do j=1,natom
!             xx(j)=xx(j)/convl
!             yy(j)=yy(j)/convl
!             zz(j)=zz(j)/convl
!          enddo

   return
   end subroutine
!
end module
!
