#include "dprec.fh"
#include "assert.fh"

module nacT_analytic_module

   use naesmd_constants
   use langevin_temperature
   use communism
   use dcart_xpm_module
   implicit none

   type xstep_t
      _REAL_,allocatable::R(:,:) ! coordinates before any increment 
      _REAL_,allocatable:: Rp(:,:) ! coordinates at + step in (3,:) shape
      _REAL_,allocatable::Rm(:,:) ! coordinates at - step 
   end type xstep_t
   

   contains
!
   subroutine xstep_A2au(xs) ! conversion from Angstrom to atomic units
   implicit none

   type(xstep_t),pointer :: xs

   xs%Rp(:,:)=xs%Rp(:,:)/convl
   xs%Rm(:,:)=xs%Rm(:,:)/convl
   xs%R(:,:)=xs%R(:,:)/convl

   return
   end subroutine
!
   subroutine xstep_au2A(xs) ! conversion from atomic units to Angstrom
   implicit none

   type(xstep_t),pointer::xs

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
   subroutine new_xstep(sim, xx, yy, zz, xxp, yyp, zzp, xxm, yym, zzm, xs)
   implicit none

   type(xstep_t),pointer::xs

   type(simulation_t),pointer::sim
   _REAL_,intent(in)::xx(sim%Na), yy(sim%Na), zz(sim%Na)
   _REAL_,intent(in)::xxp(sim%Na),yyp(sim%Na),zzp(sim%Na)
   _REAL_,intent(in)::xxm(sim%Na),yym(sim%Na),zzm(sim%Na)


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
   end subroutine
!
!********************************************************************
!
!  This is a new version of subrotine xxpxxm.
!  It calculates coordinates at time +/-dtnact
!  JAKB, there are many pointers in here making it confusing
!********************************************************************
!
   subroutine new_xstep_dtnact(sim,xxx,yyy,zzz,xs)
   use naesmd_module
   implicit none

   type(simulation_t),pointer::sim
   type(xstep_t),pointer::xs
   _REAL_,target,intent(in)::xxx(sim%Na),yyy(sim%Na),zzz(sim%Na)
   integer k,j
   type(realp_t) :: vv(3), aa(3)
   type(realp_t) :: r(3)

   vv(1)%p=>sim%naesmd%vxold
   vv(2)%p=>sim%naesmd%vyold
   vv(3)%p=>sim%naesmd%vzold
   aa(1)%p=>sim%naesmd%axold
   aa(2)%p=>sim%naesmd%ayold
   aa(3)%p=>sim%naesmd%azold
   r(1)%p=>xxx
   r(2)%p=>yyy
   r(3)%p=>zzz


   do j=1,sim%naesmd%natom
      if(sim%naesmd%ensemble.eq.'energy'.or.sim%naesmd%ensemble.eq.'temper') then
         do k=1,3
            xs%Rp(k,j)=r(k)%p(j)+vv(k)%p(j)*sim%naesmd%dtnact &
               +aa(k)%p(j)*0.5d0*sim%naesmd%dtnact*sim%naesmd%dtnact

            xs%Rm(k,j)=r(k)%p(j)-vv(k)%p(j)*sim%naesmd%dtnact &
               -aa(k)%p(j)*0.5d0*sim%naesmd%dtnact*sim%naesmd%dtnact
         end do

      else if(sim%naesmd%ensemble.eq.'langev') then
         do k=1,3
            xs%Rp(k,j)=r(k)%p(j)+vv(k)%p(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtnact  &
               +aa(k)%p(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt)*sim%naesmd%dtnact*sim%naesmd%dtnact     &
               +sim%naesmd%prand(k,j)/sim%naesmd%dtmdqt*sim%naesmd%dtnact

            xs%Rm(k,j)=r(k)%p(j)-vv(k)%p(j)*sim%naesmd%vfric(j)/sim%naesmd%dtmdqt*sim%naesmd%dtnact  &
               -aa(k)%p(j)*sim%naesmd%afric(j)/(sim%naesmd%dtmdqt*sim%naesmd%dtmdqt)*sim%naesmd%dtnact*sim%naesmd%dtnact     &
               -sim%naesmd%prand(k,j)/sim%naesmd%dtmdqt*sim%naesmd%dtnact
         end do
      end if
   end do

   xs%R(1,:)=xxx(:)
   xs%R(2,:)=yyy(:)
   xs%R(3,:)=zzz(:)


   return
   end subroutine
!
!********************************************************************
!
!********************************************************************
!
   subroutine new_xstep_dtnact_r3(sim, r, xs)
   use naesmd_module, only : realp_t
   implicit none

   type(xstep_t),pointer::xs

   type(simulation_t), pointer::sim
   type(realp_t), intent(in)::r(3)

   call new_xstep_dtnact(sim,r(1)%p,r(2)%p,r(3)%p, xs)

   return
   end subroutine

!********************************************************************
!
!  This version is similar to the nacT_analytic in the old code
!  main body is copied from sqm/nacr.F90
!
!********************************************************************
!
   function nacT_direct_ihc(sim,ihop,icheck,xstep)

   use qm2_davidson_module
   use naesmd_constants, only : kcalev
  
   implicit none

   _REAL_ nacT_direct_ihc ! function result
   type(simulation_t),pointer::sim
   integer,intent(in)::ihop,icheck
   type(xstep_t),pointer::xstep

   ! NEW OR MODIFIED VARIABLES
   integer M4_M,M2_M,Np,Nh,Nb,Mx_M,Mx

   ! OLD VARIABLES (MAY BE DEPRECATED)
   _REAL_ fPi,fbar,ff0,ff1,ff11
   integer one,N3

   parameter(fPi=3.1415926535898d0)
   parameter(fbar=0.05d0)  ! maximum dE in numerical derivative, eV.
   integer Na
        
   Na=sim%qmmm%nquant_nlink ! number of atoms
   N3=3*sim%qmmm%nquant_nlink ! number of degrees of freedom

   one=1
   ff0=0.0
   ff1=1.0
   ff11=-1.0
        
   Na=sim%qmmm%nquant_nlink
   Np=sim%dav%Np
   Nh=sim%dav%Nh
   Nb=sim%dav%Nb
   M4_M=sim%dav%Np*sim%dav%Nh
   M2_M=M4_M*2
   Mx=sim%dav%Mx 
   Mx_M=Mx

   ! Note, Davidson must be run prior to this for these assignments
   if(sim%qmmm%qm_mm_first_call) then
      write(6,*)  'sqm_energy() must be run once before executing this procedure!'

      if(sim%dav%Mx==0) write(6,*)  'excN must be > 0 to run this procedure!'
      call mexit(6,1)
   end if

   call getmodef(M2_M,Mx_M,Np,Nh,ihop,sim%dav%cmdqt,sim%dav%nacr_scratch)
   call getmodef(M2_M,Mx_M,Np,Nh,icheck,sim%dav%cmdqt,sim%dav%eta_scratch)
        
   call dgemm('N','T',Nb,Nb,Nb,ff1,sim%dav%nacr_scratch,Nb, &
      sim%dav%eta_scratch,Nb,ff0,sim%dav%eta,Nb)
   call dgemm('T','N',Nb,Nb,Nb,ff1,sim%dav%eta_scratch,Nb, &
      sim%dav%nacr_scratch,Nb,ff1,sim%dav%eta,Nb)
   call Iminus2rho(Nb,Np,sim%dav%eta,sim%dav%xi)
   call mo2sitef (Nb,sim%dav%vhf,sim%dav%xi,sim%dav%eta,sim%dav%xi_scratch)
        
   ! Above eta contains transition density martix between sim%naesmd%state 
   ! ihop and icheck in AO

   ! Above xi and xi_scratch_2 contain transition density martices 
   ! for states ihop and icheck in AO

   ! Term Tr(F^x rho_ij) (only symmetric part contributes)
   call packing(Nb,sim%dav%eta,sim%dav%nacr_scratch,'s')


   nacT_direct_ihc=dcart1_xpm(sim%qnml, sim%qparams, sim%qmpi, sim%qm2, sim%dav, sim%qmmm, sim%dav%nacr_scratch,xstep%Rp,xstep%Rm) 


   nacT_direct_ihc=nacT_direct_ihc*kcalev &
      /(sim%dav%e0(sim%dav%kx(icheck))-sim%dav%e0(sim%dav%kx(ihop))) &
      / sim%naesmd%dtnact/2.d0 
   ! factor of 2. is because dt = 2.0 * sim%naesmd%dtnact



   flush(6)

   return
   end function nacT_direct_ihc
!
!********************************************************************
!  Wrapper for NACT called from the program 
!********************************************************************
!
   subroutine nacT_direct(sim,nact,xs)
   
   use dcart_xpm_module
   use constants          , only : EV_TO_KCAL
   use ElementOrbitalIndex, only: MaxValenceOrbitals
   use qm2_pm6_hof_module
   use dh_correction_module, only : dh_correction_grad
   use qm2_davidson_module ! CML 7/13/12
   
   implicit none

   type(simulation_t),pointer::sim
   _REAL_,intent(inout)::nact(:,:)
  
   _REAL_ psum(MaxValenceOrbitals**2*3)
   _REAL_ xyz_qmi(3), xyz_qmj(3)
   _REAL_ :: Fm(MaxValenceOrbitals*(MaxValenceOrbitals*2+1)) !BTN 10/08/2017 place to store fock matrix

   integer natqmi, natqmj, qmitype, qmjtype
   integer ii, iif, iil, jj, jjf, jjl, ij
   integer i,j,k,l
   integer n_atomic_orbi, n_atomic_orbj
   integer jstart, jend
   type(xstep_t),pointer::xs

   nact=0.d0
   
   call xstep_au2A(xs)
   
   !BTN Start here
   
#ifdef MPI
   do ii = sim%qmpi%nquant_nlink_istart, sim%qmpi%nquant_nlink_iend
      jstart =  sim%qmpi%nquant_nlink_jrange(1,ii)
      jend = sim%qmpi%nquant_nlink_jrange(2,ii)
#else
   do II=2,sim%qmmm%nquant_nlink
       jstart = 1
       jend = ii-1
#endif
       !Loop over all pairs of quantum atoms
       iif=sim%qparams%orb_loc(1,II)
       iil=sim%qparams%orb_loc(2,II)
       qmitype = sim%qmmm%qm_atom_type(ii)
       natqmi=sim%qmmm%iqm_atomic_numbers(II)
       do JJ=jstart,jend !jj=1,ii-1
           !  FORM DIATOMIC MATRICES
           jjf=sim%qparams%orb_loc(1,JJ)
           jjl=sim%qparams%orb_loc(2,JJ)
           !   GET FIRST ATOM
           qmjtype = sim%qmmm%qm_atom_type(jj)
           natqmj=sim%qmmm%iqm_atomic_numbers(JJ)
           IJ=0
           do I=jjf,jjl
               K=sim%qparams%pascal_tri1(i)+jjf-1
               do J=jjf,I
                   IJ=IJ+1
                   K=K+1
                   psum(IJ)=sim%qm2%den_matrix(K)               
               end do
           end do
           ! GET SECOND ATOM FIRST ATOM INTERSECTION
           do I=iif,iil
               L=sim%qparams%pascal_tri1(i)
               K=L+jjf-1
               do J=jjf,jjl
                   IJ=IJ+1
                   K=K+1
                   psum(IJ)=sim%qm2%den_matrix(K)
               end do
               K=L+iif-1
               do L=iif,I
                   K=K+1
                   IJ=IJ+1
                   psum(IJ)=sim%qm2%den_matrix(K)
               end do
           end do
           
           !BTN Changed to remove fock calculation from multiplication with transition density
           
           !BTN zero temp array
           Fm=0.0d0
           
           xyz_qmi(:)=xs%rm(:,ii)
           xyz_qmj(:)=xs%rm(:,jj)
           call qm2_dhc1(sim%rij,sim%qparams, sim%qnml, sim%qm2, sim%qmmm, psum,ii,jj, &
              qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi,natqmj,iif,iil,jjf, &
              jjl,Fm)
                    
           xyz_qmi(:)=xs%rp(:,ii)
           xyz_qmj(:)=xs%rp(:,jj)
           call qm2_dhc1(sim%rij,sim%qparams, sim%qnml, sim%qm2, sim%qmmm, psum,ii,jj,&
              qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi,natqmj,iif,iil,jjf, &
              jjl,sim%qm2%fock_matrix_dp(sim%qparams%pascal_tri1(ii-1)+jj,:))
              
           sim%qm2%fock_matrix_dp(sim%qparams%pascal_tri1(ii-1)+jj,:)= &
              sim%qm2%fock_matrix_dp(sim%qparams%pascal_tri1(ii-1)+jj,:)-Fm(:)
           
       end do
   end do
   
   
   !BTN End here

   !FIXME nacT_direct_ihc calls part of derivative each time.. is this necessary
   !or calculated previously?
   do i=2,sim%excN
       do j=1,i-1
           nact(i,j)=nacT_direct_ihc(sim,i,j,xs) ! notice i corresponds to 
                                                       ! ihop, i - icheck
           nact(j,i)=-nact(i,j)
       end do
   end do

   call xstep_A2au(xs)
   
   return
   end subroutine nacT_direct


  ! Kind of a useless wrapper for nact_direct isn't it?
   subroutine nacT_analytic(sim,nact,xstep)
        implicit none
        type(simulation_t),pointer::sim
        type(xstep_t),pointer::xstep
        _REAL_,intent(inout)::nact(:,:)
        _REAL_ t_start,t_finish
        call cpu_time(t_start)
        call nacT_direct(sim,nact,xstep)
        call cpu_time(t_finish)
        sim%time_nact_took=sim%time_nact_took+t_finish-t_start
        return
   end subroutine

end module!
