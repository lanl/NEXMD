#include "dprec.fh"
#include "assert.fh"

module nacT_analytic_module

   use naesmd_constants
   use langevin_temperature
   use communism
   use naesmd_space_module
   use dcart_xpm_module

   implicit none

   type xstep_t
      _REAL_,pointer::R(:,:) ! coordinates before any increment 
      _REAL_,pointer:: Rp(:,:) ! coordinates at + step in (3,:) shape
      _REAL_,pointer::Rm(:,:) ! coordinates at - step 
   end type xstep_t

   contains
!
   subroutine xstep_A2au(xs) ! conversion from Angstrom to atomic units
   implicit none

   type(xstep_t), intent(inout) :: xs

   xs%Rp(:,:)=xs%Rp(:,:)/convl
   xs%Rm(:,:)=xs%Rm(:,:)/convl
   xs%R(:,:)=xs%R(:,:)/convl

   return
   end subroutine
!
   subroutine xstep_au2A(xs) ! conversion from atomic units to Angstrom
   implicit none

   type(xstep_t),intent(inout)::xs

   xs%Rp(:,:)=xs%Rp(:,:)*convl
   xs%Rm(:,:)=xs%Rm(:,:)*convl
   xs%R(:,:)=xs%R(:,:)*convl

   return
   end subroutine
!
!********************************************************************
!
!********************************************************************
!
   function new_xstep(sim, xx, yy, zz, xxp, yyp, zzp, xxm, yym, zzm) result(xs)
   implicit none

   type(xstep_t)::xs

   type(simulation_t),pointer::sim
   _REAL_,intent(in)::xx(sim%Na), yy(sim%Na), zz(sim%Na)
   _REAL_,intent(in)::xxp(sim%Na),yyp(sim%Na),zzp(sim%Na)
   _REAL_,intent(in)::xxm(sim%Na),yym(sim%Na),zzm(sim%Na)

   allocate(xs%Rp(3,sim%Na))
   allocate(xs%Rm(3,sim%Na))
   allocate(xs%R(3,sim%Na))

   xs%R(1,:)=xx(:)
   xs%R(2,:)=yy(:)
   xs%R(3,:)=zz(:)

   xs%Rp(1,:)=xxp(:)
   xs%Rp(2,:)=yyp(:)
   xs%Rp(3,:)=zzp(:)

   xs%Rm(1,:)=xxm(:)
   xs%Rm(2,:)=yym(:)
   xs%Rm(3,:)=zzm(:)

   return
   end function
!
!********************************************************************
!
!  This is a new version of subrotine xxpxxm.
!  It calculates coordinates at time +/-dtnact
!
!********************************************************************
!
   function new_xstep_dtnact(sim,xx,yy,zz) result(xs)
   implicit none

   type(simulation_t),pointer::sim
   type(xstep_t)::xs
   _REAL_,target,intent(in)::xx(sim%Na),yy(sim%Na),zz(sim%Na)
   integer k,j,i,ii,iii,iimdqt
   _REAL_ x 

   include 'md.par'
   include 'sizes'
   include 'common'

   !_REAL_ :: dtnact, dtmdqt

   type(realp_t),pointer::v(:),a(:)
   type(realp_t)::r(3)


   !dtnact = sim%naemsd%dtnact
   !dtmdqt = sim%naemsd%dtmdqt

   v=>sim%naesmd%v%vold
   a=>sim%naesmd%a%vold
   r(1)%p=>xx
   r(2)%p=>yy
   r(3)%p=>zz

   allocate(xs%Rp(3,sim%Na))
   allocate(xs%Rm(3,sim%Na))
   allocate(xs%R(3,sim%Na))


   do j=1,natom
      if(ensemble.eq.'energy'.or.ensemble.eq.'temper') then
         do k=1,3
            xs%Rp(k,j)=r(k)%p(j)+v(k)%p(j)*dtnact &
               +a(k)%p(j)*0.5d0*dtnact*dtnact

            xs%Rm(k,j)=r(k)%p(j)-v(k)%p(j)*dtnact &
               -a(k)%p(j)*0.5d0*dtnact*dtnact
         end do

      else if(ensemble.eq.'langev') then
         do k=1,3
            xs%Rp(k,j)=r(k)%p(j)+v(k)%p(j)*vfric(j)/dtmdqt*dtnact  &
               +a(k)%p(j)*afric(j)/(dtmdqt*dtmdqt)*dtnact*dtnact     &
               +prand(k,j)/dtmdqt*dtnact

            xs%Rm(k,j)=r(k)%p(j)-v(k)%p(j)*vfric(j)/dtmdqt*dtnact  &
               -a(k)%p(j)*afric(j)/(dtmdqt*dtmdqt)*dtnact*dtnact     &
               -prand(k,j)/dtmdqt*dtnact
         end do
      end if
   end do

   xs%R(1,:)=xx(:)
   xs%R(2,:)=yy(:)
   xs%R(3,:)=zz(:)

   return
   end function
!
!********************************************************************
!
!********************************************************************
!
   function new_xstep_dtnact_r3(sim, r) result(xs)
   implicit none

   type(xstep_t)::xs

   type(simulation_t), pointer::sim
   type(realp_t), intent(in)::r(3)

   xs=new_xstep_dtnact(sim,r(1)%p,r(2)%p,r(3)%p)

   return
   end function
!
!********************************************************************
!
!  KGB this version calculates nacR and multiplies it by nuclear velocities 
!  hovewer the result is diffirent from 'direct' version. They should
!  coincide when dt -> 0. For the direct version the +/- coordinates are
!  calculated with accelerations and thermostat.
! 
!  JAKB This is not currently used or tested
!********************************************************************
!
   subroutine nacT_v(sim,nact)
   use qmmm_module,only:qmmm_nml

   implicit none

   type(simulation_t),pointer :: sim
   _REAL_,intent(inout)::nact(:,:)
   integer::i,j,k
   _REAL_::s

   call qmmm2coords_r(sim)
   nact = 0.d0
   write(6,*)'Called nacT_v'

   do i=2,sim%naesmd%npot
      do j=1,i-1
         call nacR_analytic(sim%coords,j,i) ! notice j corresponds to 
                                                   ! ihop, i - icheck
         do k=1,sim%Na
            s=sim%naesmd%v%xold(k)*sim%dav%dij(k*3-2) &
               +sim%naesmd%v%yold(k)*sim%dav%dij(k*3-1) &
               +sim%naesmd%v%zold(k)*sim%dav%dij(k*3)

               s=s*convl !FIXME unit conversion for dij?
               nact(i,j) = nact(i,j) + s
         end do

         nact(j,i)=-nact(i,j)
      end do
   end do

   k=1

   return
   end subroutine
!
!********************************************************************
!
!  KGB moved this subroutine from cadiab.F90
!  where it was formerly used.
!  Calculates coordinates at time +/-dtnact.
!
!********************************************************************
!
   subroutine xxpxxm(xx,yy,zz,xxp,yyp,zzp,xxm,yym,zzm)
   implicit none

   include 'md.par'
   include 'sizes'
   include 'common'

   integer k,j,i,ii,iii,iimdqt
   _REAL_ x 

   _REAL_,intent(in)::xx(Na_M),yy(Na_M),zz(Na_M)
   _REAL_ xxp(Na_M),yyp(Na_M),zzp(Na_M)
   _REAL_ xxm(Na_M),yym(Na_M),zzm(Na_M)

   do j=1,natom
      if(ensemble.eq.'energy'.or. ensemble.eq.'temper') then
         xxp(j)=xx(j)+vxold(j)*dtnact+axold(j)*0.5d0*dtnact*dtnact
         yyp(j)=yy(j)+vyold(j)*dtnact+ayold(j)*0.5d0*dtnact*dtnact
         zzp(j)=zz(j)+vzold(j)*dtnact + azold(j)*0.5d0*dtnact*dtnact

! xxm,m,yym, and zzm are xyz at t - dtnact
         xxm(j)=xx(j)-vxold(j)*dtnact-axold(j)*0.5d0*dtnact*dtnact
         yym(j)=yy(j)-vyold(j)*dtnact-ayold(j)*0.5d0*dtnact*dtnact
         zzm(j)=zz(j)-vzold(j)*dtnact-azold(j)*0.5d0*dtnact*dtnact

      else if(ensemble.eq.'langev') then
         xxp(j)=xx(j)+vxold(j)*vfric(j)/dtmdqt*dtnact      &
            +axold(j)*afric(j)/(dtmdqt*dtmdqt)*dtnact*dtnact &
            +prand(1,j)/dtmdqt*dtnact

         yyp(j)=yy(j)+vyold(j)*vfric(j)/dtmdqt*dtnact &
            +ayold(j)*afric(j)/(dtmdqt*dtmdqt)*dtnact*dtnact &
            +prand(2,j)/dtmdqt*dtnact

         zzp(j)=zz(j)+vzold(j)*vfric(j)/dtmdqt*dtnact &
            +azold(j)*afric(j)/(dtmdqt*dtmdqt)*dtnact*dtnact &
            +prand(3,j)/dtmdqt*dtnact

         xxm(j)=xx(j)-vxold(j)*vfric(j)/dtmdqt*dtnact &
            -axold(j)*afric(j)/(dtmdqt*dtmdqt)*dtnact*dtnact &
            -prand(1,j)/dtmdqt*dtnact

         yym(j)=yy(j)-vyold(j)*vfric(j)/dtmdqt*dtnact &
            -ayold(j)*afric(j)/(dtmdqt*dtmdqt)*dtnact*dtnact &
            -prand(2,j)/dtmdqt*dtnact

         zzm(j)=zz(j)-vzold(j)*vfric(j)/dtmdqt*dtnact &
            -azold(j)*afric(j)/(dtmdqt*dtmdqt)*dtnact*dtnact &
            -prand(3,j)/dtmdqt*dtnact

      end if

      xxp(j)=xxp(j)*convl
      yyp(j)=yyp(j)*convl
      zzp(j)=zzp(j)*convl
      xxm(j)=xxm(j)*convl
      yym(j)=yym(j)*convl
      zzm(j)=zzm(j)*convl
   end do

   return
   end subroutine
!
!********************************************************************
!
!  This version is similar to the nacT_analytic in the old code
!  main body is copied from sqm/nacr.F90
!
!********************************************************************
!
   function nacT_direct_ihc(sim,ihop,icheck,xstep)

   use qmmm_module,only:qmmm_struct,qm2_struct
   use qm2_davidson_module
   use naesmd_constants, only : kcalev

   implicit none

   _REAL_ nacT_direct_ihc ! function result
   type(simulation_t),pointer::sim
   integer,intent(in)::ihop,icheck
   type(xstep_t),intent(inout)::xstep

   ! NEW OR MODIFIED VARIABLES
   integer M4_M,M2_M,Np,Nh,Nb,Mx_M,Mx
   integer ideriv ! This will come from a module later
   !_REAL_, intent(in) :: xyz_in(3*qmmm_struct%nquant_nlink)

   ! OLD VARIABLES (MAY BE DEPRECATED)
   _REAL_ fPi,fbar,f,t,d,d1,ff,ff0,ff1,ff11,ddot
   integer i,ii,j,k,im,one,istate,mdflag,ip,ih,Nm,N3

   parameter(fPi=3.1415926535898d0)
   parameter(fbar=0.05d0)  ! maximum dE in numerical derivative, eV.
   integer Na,npot
        
   !include 'md.par'

   Na=qmmm_struct%nquant_nlink ! number of atoms
   N3=3*qmmm_struct%nquant_nlink ! number of degrees of freedom

   one=1
   ff0=0.0
   ff1=1.0
   ff11=-1.0
        
   Na=qmmm_struct%nquant_nlink
   Np=qm2ds%Np
   Nh=qm2ds%Nh
   Nb=qm2ds%Nb
   M4_M=qm2ds%Np*qm2ds%Nh
   M2_M=M4_M*2
   Mx=qm2ds%Mx 
   Mx_M=Mx
   npot=Mx

   ! Note, Davidson must be run prior to this for these assignments
   if(qmmm_struct%qm_mm_first_call) then
      write(6,*)  'sqm_energy() must be run once before executing this procedure!'

      if(qm2ds%Mx==0) write(6,*)  'excN must be > 0 to run this procedure!'
      call mexit(6,1)
   end if

   call getmodef(M2_M,Mx_M,Np,Nh,ihop,qm2ds%cmdqt,qm2ds%nacr_scratch)
   call getmodef(M2_M,Mx_M,Np,Nh,icheck,qm2ds%cmdqt,qm2ds%eta_scratch)
        
   call dgemm('N','T',Nb,Nb,Nb,ff1,qm2ds%nacr_scratch,Nb, &
      qm2ds%eta_scratch,Nb,ff0,qm2ds%eta,Nb)
   call dgemm('T','N',Nb,Nb,Nb,ff1,qm2ds%eta_scratch,Nb, &
      qm2ds%nacr_scratch,Nb,ff1,qm2ds%eta,Nb)
   call Iminus2rho(Nb,Np,qm2ds%eta,qm2ds%xi)
   call mo2sitef (Nb,qm2ds%vhf,qm2ds%xi,qm2ds%eta,qm2ds%xi_scratch)
        
   ! Above eta contains transition density martix between state 
   ! ihop and icheck in AO
   write(6,*)ihop,icheck,'nacT AO density matrix', qm2ds%eta

   call getmodef(M2_M,Mx_M,Np,Nh,ihop,qm2ds%cmdqt,qm2ds%nacr_scratch)
   call mo2sitef (Nb,qm2ds%vhf,qm2ds%nacr_scratch,qm2ds%xi,qm2ds%eta_scratch)
   call getmodef(M2_M,Mx_M,Np,Nh,icheck,qm2ds%cmdqt,qm2ds%nacr_scratch)
   call mo2sitef (Nb,qm2ds%vhf,qm2ds%nacr_scratch, &
      qm2ds%xi_scratch_2,qm2ds%eta_scratch)

   ! Above xi and xi_scratch_2 contain transition density martices 
   ! for states ihop and icheck in AO

   ! Term Tr(F^x rho_ij) (only symmetric part contributes)
   call packing(Nb,qm2ds%eta,qm2ds%nacr_scratch,'s')


   nacT_direct_ihc=dcart1_xpm(qm2_struct%den_matrix, &
      qm2ds%nacr_scratch,xstep%Rp,xstep%Rm) 

   write(6,*)'nacT_deriv+denmat',nacT_direct_ihc

   nacT_direct_ihc=nacT_direct_ihc*kcalev &
      /(qm2ds%e0(qm2ds%kx(icheck))-qm2ds%e0(qm2ds%kx(ihop))) &
      / sim%naesmd%dtnact/2.d0 
   ! factor of 2. is because dt = 2.0 * dtnact

   write(6,*)'nacT',nacT_direct_ihc


   flush(6)

   return
   end function nacT_direct_ihc
!
!********************************************************************
!  Wrapper for NACT called from the program 
!********************************************************************
!
   subroutine nacT_direct(sim,nact,xs)
   implicit none

   type(simulation_t),pointer::sim
   _REAL_,intent(inout)::nact(:,:)

   integer i,j
   _REAL_ s
   type(xstep_t),intent(inout)::xs

   !xs = new_xstep_dtnact_r3(sim, sim%naesmd%r%vold)

   !call qmmm2coords_r(sim)
   nact=0.d0

   call xstep_au2A(xs)
   !FIXME nacT_direct_ihc calls derivatives each time.. not necessary since
   !always the same
   do i=2,sim%naesmd%npot
      do j=1,i-1
         nact(i,j)=nacT_direct_ihc(sim,i,j,xs) ! notice i corresponds to 
                                                        ! ihop, i - icheck
         nact(j,i)=-nact(i,j)
      end do
   end do

   call xstep_A2au(xs)
   
   return
   end subroutine nacT_direct
!
!********************************************************************
!
!********************************************************************
!
   subroutine nacT_analytic(sim,nact,xstep,v_version)
   implicit none

   type(simulation_t),pointer::sim
   type(xstep_t),intent(inout)::xstep
   _REAL_,intent(inout)::nact(:,:)

   integer,optional::v_version
   _REAL_ t_start,t_finish   

   !FIXME as of 6/2016 this is not used in the program JAKB     
   call cpu_time(t_start)
   write(6,*)'nacT_analytic t_start=',t_start
   if(present(v_version)) then
      if(v_version.ne.0) then
         call nacT_v(sim,nact)
         write(6,*)'nacT_v called:',nact
         return
      end if
   end if

   call nacT_direct(sim,nact,xstep)
   write(6,*)'nacT_direct called:',nact

   call cpu_time(t_finish)
   write(6,*)'nacT_analytic t_finish=',t_finish
   sim%time_nact_took=sim%time_nact_took+t_finish-t_start 

   return
   end subroutine
!
end module
!
