! <compile=optimized>

module ew_recip_spatial

#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"

_REAL_, allocatable :: m1_tbl(:)
_REAL_, allocatable :: m2_tbl(:)
_REAL_, allocatable :: m3_tbl(:)
integer,save :: sp_len1,sp_len2,sp_len3

interface 
   subroutine dumpq(a,n1,n2,n3,lun,p_type,n10,n20,n30)
     implicit none
     integer,intent(in) :: n1,n2,n3,lun
     integer,optional :: p_type,n10,n20,n30
     double complex,intent(in) :: a(n1,n2,n3)
   end subroutine dumpq
end interface
character(len=40) :: fmt0='(i4,1x,a,10i5)', &
                     fmt1='(2i5,2e15.5,10i5)'


#ifdef MPI


public deallocate_m1m2m3, spatial_do_pmesh_kspace

!===========================================================================
contains
!===========================================================================


!------------------------------------------------------------------
!           SPATIAL  DO_PMESH_KSPACE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spatial_do_pmesh_kspace( natom,crd,charge, &
      frc,prefac1,prefac2,prefac3,fftable,qm_pot_only)
   use nblist, only: recip,volume
   use qmmm_module, only:qmmm_nml,qmmm_struct, qmewald
   use stack
   use fft, only:backward_rc_fft,forward_rc_fft,get_fft_limits, &
         YZ_X_PARTITION,XY_Z_PARTITION
   use ew_recip
   
   implicit none
   character(kind=1,len=*), parameter :: routine="spatial_do_pmesh_kspace"
#  include "flocntrl.h"
#  include "ew_pme_recip.h"
#  include "ew_frc.h"
#  include "def_time.h"
#  include "box.h"

#  include "parallel.h"
#  include "ew_parallel.h"
   include 'mpif.h'


   ! INPUT
   !       natom:  number of atoms
   !       crd   atomic coords
   !       charge  atomic charges

   integer ierr
   
   integer natom,num_ks_trial
   _REAL_ crd(3,natom),charge(natom)
   integer nmine,nderiv
   
   logical, intent(in) :: qm_pot_only
   
   ! OUTPUT
   !       eer:  ewald reciprocal or k-space  energy
   !       frc forces incremented by k-space sum
   
   _REAL_ frc(3,natom)
   
   ! HEAP STORAGE:  These arrays need to be preserved throughout simulation
   
   _REAL_ prefac1(*),prefac2(*),prefac3(*),fftable(*)
   
   integer l_d2th1,l_d2th2,l_d2th3
   
   integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw
   integer l_fr1,l_fr2,l_fr3
   integer l_th1,l_th2,l_th3,l_dth1,l_dth2,l_dth3
   integer l_fftw,l_q
   integer imy_cg
   integer num_ks
   integer i,ii
   save num_ks
   data num_ks/0/
   integer l_tmpy,l_alpha,l_beta
   logical not_done
   logical not_done_sp
   integer num_ks_sp
   double complex,dimension(:),allocatable :: q_spat
   integer :: xgridmin,xgridmax,ygridmin,ygridmax,zgridmin,zgridmax

   not_done_sp = .true.
   nderiv = 1
   
   !     FLOW CONTROL FLAG (debug, possible vacuum sims and future respa)
   
   if ( do_rec == 0 )return
   
   !     get some integer array dimensions
   call get_fftdims(nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
   if(num_ks == 0) &
         num_ks = min((((nxyslab(0)+order-1)*natom*4)/ &
         (3*nfft3)),natom)
   call timer_start(TIME_BSPL)
   num_ks_trial = 0
   not_done = .true.
   do while(not_done)
      call get_stack(l_fftw,sizffwrk,routine)
      call get_stack(l_q,siz_q,routine)
      call get_stack(l_fr1,num_ks,routine)
      call get_stack(l_fr2,num_ks,routine)
      call get_stack(l_fr3,num_ks,routine)
      call get_stack(l_th1,num_ks*order,routine)
      call get_stack(l_th2,num_ks*order,routine)
      call get_stack(l_th3,num_ks*order,routine)
      if ( nderiv == 1 )then
         call get_stack(l_dth1,num_ks*order,routine)
         call get_stack(l_dth2,num_ks*order,routine)
         call get_stack(l_dth3,num_ks*order,routine)
         call get_stack(l_d2th1,order,routine)
         call get_stack(l_d2th2,order,routine)
         call get_stack(l_d2th3,order,routine)
      end if
      if(.not. rstack_ok)then
         deallocate(r_stack)
         allocate(r_stack(1:lastrst),stat=alloc_ier)
         call reassign_rstack(routine)
      endif
      REQUIRE(rstack_ok)

      call get_istack(imy_cg,num_ks,routine)
      if(.not. istack_ok)then
         deallocate(i_stack)
         allocate(i_stack(1:lastist),stat=alloc_ier)
         call reassign_istack(routine)
      endif
      REQUIRE(istack_ok)

      call spatial_get_grid_weights( &
            natom,crd,recip,nfft1,nfft2,nfft3, &
            r_stack(l_fr1),r_stack(l_fr2),r_stack(l_fr3),order, &
            r_stack(l_th1),r_stack(l_th2),r_stack(l_th3), &
            r_stack(l_dth1),r_stack(l_dth2),r_stack(l_dth3), &
            r_stack(l_d2th1),r_stack(l_d2th2),r_stack(l_d2th3), &
            i_stack(imy_cg),nmine,nderiv,num_ks)
      if(nmine > num_ks)then
         write(6,'("*****  Processor ",i6)') mytaskid
         write(6,'("***** System must be very inhomogeneous.")')
         write(6,'("*****  Readjusting recip sizes.")')
         write(6,'(A,i9,A,i9/)') &
               " In this slab, Atoms found: ",nmine, &
               "  Allocated: ",num_ks
         if( num_ks_trial >= 2 ) call mexit( 6,1 )
         if ( nderiv == 1 )then
            call free_stack(l_d2th3,routine)
            call free_stack(l_d2th2,routine)
            call free_stack(l_d2th1,routine)
            call free_stack(l_dth3,routine)
            call free_stack(l_dth2,routine)
            call free_stack(l_dth1,routine)
         end if
         call free_stack(l_th3,routine)
         call free_stack(l_th2,routine)
         call free_stack(l_th1,routine)
         call free_stack(l_fr3,routine)
         call free_stack(l_fr2,routine)
         call free_stack(l_fr1,routine)
         call free_stack(l_q,routine)
         call free_stack(l_fftw,routine)
         call free_istack(imy_cg,routine)
         num_ks= nmine*4/3
      else
         not_done = .false.
      end if
   enddo
   call timer_stop_start(TIME_BSPL,TIME_FILLG)
   
   !........Fill Charge Grid
   !          charges are approximated on an even grid
   ! Using nfftdim1-1 for the stride through x data lines to account
   !  for the nfft1/2+1 size after an r-c fft. The other stides are
   !  just nfft2 and nfft3
   if(.not. allocated(q_spat))then
      allocate(q_spat(siz_q),stat=ierr)
      REQUIRE(ierr == 0)
   endif
   call spatial_fill_charge_grid(natom,charge, &
         r_stack(l_th1),r_stack(l_th2),r_stack(l_th3), &
         r_stack(l_fr1),r_stack(l_fr2),r_stack(l_fr3), &
         order, &
         nfft1,nfft2,nfft3, &
         nfftdim1-1,nfft2,nfft3, &
         q_spat,i_stack(imy_cg),nmine,siz_q/2)
   call timer_stop_start(TIME_FILLG,TIME_FFT)
   call get_stack(l_tmpy,2*nfftdim1,routine)
   call get_stack(l_alpha,nfft1,routine)
   call get_stack(l_beta,nfft1,routine)
   if(.not. rstack_ok)then
      deallocate(r_stack)
      allocate(r_stack(1:lastrst),stat=alloc_ier)
      call reassign_rstack(routine)
   endif
   REQUIRE(rstack_ok)
   
   call get_fft_limits(YZ_X_PARTITION,xgridmin,xgridmax, &
         ygridmin,ygridmax,zgridmin,zgridmax,mytaskid)
   xgridmin=0
   xgridmax=nfftdim1
   call backward_rc_fft(q_spat)
   call get_fft_limits(XY_Z_PARTITION,zgridmin,zgridmax,xgridmin,xgridmax, &
         ygridmin,ygridmax,mytaskid)
   
   call timer_stop_start(TIME_FFT,TIME_SCSUM)
   
   !           -------------SCALAR SUM------------------
   
   if( ifbox == 1 ) then
      call spatial_scalar_sumrc_orthog( &
            q_spat,ew_coeff,volume,recip, &
            prefac1,prefac2,prefac3, &
            nfft1,nfft2,nfft3,xgridmax-xgridmin+1,ygridmax-ygridmin+1, & 
            eer,rec_vir)
   else
      call sander_bomb("Broken","broken","broken")
            call spatial_scalar_sumrc( &
            q_spat,ew_coeff,volume,recip, &
            prefac1,prefac2,prefac3, &
            nfft1,nfft2,nfft3,xgridmax-xgridmin+1,ygridmax-ygridmin+1, &
            eer,rec_vir)
   end if
   call timer_stop_start(TIME_SCSUM,TIME_FFT)
   
   !-----------FFT FORWARD--------------------
   
   call forward_rc_fft(q_spat,nfftdim1-1,nfft2,nfft3)
   
   call free_stack(l_beta,routine)
   call free_stack(l_alpha,routine)
   call free_stack(l_tmpy,routine)

   call timer_stop_start(TIME_FFT,TIME_GRADS)
   
   !-----------SPATIAL GRAD SUM--------------------
   call get_fft_limits(YZ_X_PARTITION,xgridmin,xgridmax,ygridmin,ygridmax, &
         zgridmin,zgridmax,mytaskid)
   call spatial_grad_sumrc( &
         natom,charge,recip, &
         r_stack(l_th1),r_stack(l_th2),r_stack(l_th3), &
         r_stack(l_dth1),r_stack(l_dth2),r_stack(l_dth3), &
         frc,r_stack(l_fr1),r_stack(l_fr2),r_stack(l_fr3), &
         order,nfft1,nfft2,nfft3, &
         xgridmax-xgridmin+1, ygridmax-ygridmin+1, zgridmax-zgridmin+1, &
         q_spat, &
         i_stack(imy_cg),nmine, &
         qm_pot_only )
         
   call timer_stop(TIME_GRADS)
   if ( nderiv == 1 )then
      call free_stack(l_d2th3,routine)
      call free_stack(l_d2th2,routine)
      call free_stack(l_d2th1,routine)
      call free_stack(l_dth3,routine)
      call free_stack(l_dth2,routine)
      call free_stack(l_dth1,routine)
   end if
   call free_stack(l_th3,routine)
   call free_stack(l_th2,routine)
   call free_stack(l_th1,routine)
   call free_stack(l_fr3,routine)
   call free_stack(l_fr2,routine)
   call free_stack(l_fr1,routine)
   call free_stack(l_q,routine)
   call free_stack(l_fftw,routine)

   call free_istack(imy_cg,routine)
   call mpi_barrier(recip_comm,ierr)

   deallocate(q_spat)
   return
end subroutine spatial_do_pmesh_kspace 





!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!                                -----
!     --- FILL_CHARGE_GRID -- SPATIAL-------
!                                -----
!-------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_charge_grid here]
subroutine spatial_fill_charge_grid(natom,charge, &
      theta1,theta2,theta3,fr1,fr2,fr3, &
      order,nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3,q,my_cg,nmine,siz_q)

   !---------------------------------------------------------------------
   ! INPUT:
   !      natom:  number of atoms
   !      charge: the array of atomic charges
   !      theta1,theta2,theta3: the spline coeff arrays
   !      fr1,fr2,fr3 the scaled and shifted fractional coords
   !      nfft1,nfft2,nfft3: the charge grid dimensions
   !      nfftdim1,nfftdim2,nfftdim3: physical charge grid dims
   !      order: the order of spline interpolation
   ! OUTPUT:
   !      Q the charge grid
   !---------------------------------------------------------------------
  use fft, only:get_fft_limits,YZ_X_PARTITION
  use constants,only:ZERO,ONE
  use ew_bspline,only:kbot,ktop
  
   implicit none
   integer,intent(in) :: natom,order,nfft1,nfft2,nfft3,siz_q
   integer,intent(in) :: nfftdim1,nfftdim2,nfftdim3
   _REAL_ fr1(natom),fr2(natom),fr3(natom)
   _REAL_ theta1(order,natom),theta2(order,natom), &
         theta3(order,natom),charge(natom)
   double complex q(siz_q)
   integer my_cg(*),nmine

#  include "parallel.h"
#  include "ew_parallel.h"
   integer kbot0

   integer n,ith1,ith2,ith3,i0,j0,k0,i,j,k,i00,j00,iqk,iqj
   _REAL_ prod
   integer :: im,kq,jq,iq,iqi
   integer :: xgridmin,xgridmax,ygridmin,ygridmax,zgridmin, &
         zgridmax,jbot0,jbot,jtop
   
   !........Zero the Charge grids
   q(1:siz_q)=(ZERO,ZERO)

199 format(i4,a,10i5)

   call get_fft_limits(YZ_X_PARTITION,xgridmin,xgridmax,ygridmin,ygridmax, &
         zgridmin,zgridmax,mytaskid)

   sp_len1=nfftdim1
   sp_len2=1+ygridmax-ygridmin
   sp_len3=1+zgridmax-zgridmin
   kbot0 = zgridmin
   kbot = kbot0 + 1
   ktop = zgridmax+1
   jbot0 = ygridmin
   jbot = jbot0 + 1
   jtop = ygridmax+1

   do im = 1,nmine
      n = my_cg(im)
      k0 = int(fr3(im)) - order
      j00 = int(fr2(im)) - order
      i00 = int(fr1(im)) - order
      do ith3 = 1,order
         ! Z index into complete array
         k0 = k0 + 1
         if (k0 >= 0) then
            k = k0 + 1
         else
            k = k0 + 1 + nfft3
         end if

         if ( k >= kbot .and. k <= ktop ) then
            ! Z index and offset into local array
            kq = k - zgridmin     
            iqk = (kq-1)*2*sp_len1*sp_len2
            j0 = j00
            do ith2 = 1,order
               j0 = j0 + 1
               if (j0 >= 0) then
                  j = j0 + 1
               else
                  j = j0 + 1 + nfft2
               end if
               if ( j >= jbot .and. j <= jtop ) then
                  jq = j - ygridmin
                  iqj = iqk + (jq-1)*2*sp_len1
                  
                  prod = theta2(ith2,im)*theta3(ith3,im)*charge(n)
                  i0 = i00 + 1
                  do ith1 = 1,order
                     i0 = i0 + 1
                     if (i0 >= 1) then
                        iqi=i0+iqj                    
                     else
                        iqi=i0+(nfft1)+iqj
                     end if
                     iq=ishft(iqi+1,-1)
                     if(mod(iqi,2) == 0)then
                        q(iq) = q(iq) + (ZERO,ONE)*theta1(ith1,im)*prod
                     else
                        q(iq) = q(iq) + (ONE,ZERO)*theta1(ith1,im)*prod
                     endif
                  end do
               endif
            end do
         end if      
      end do  !  ith3 = 1,order
   end do  !  im = 1,nmine
   return
end subroutine spatial_fill_charge_grid 




!!-------------------------------------------------------------------
!!-------------------------------------------------------------------
!!  GRAD_SUM RC    **** REAL not complex ****
!!-----------------------------------------------------------------
!!-------------------------------------------------------------------
subroutine spatial_grad_sumrc( &
      natom,charge,recip,theta1,theta2,theta3, &
      dtheta1,dtheta2,dtheta3,frc,fr1,fr2,fr3, &
      order,nfft1,nfft2,nfft3,xdim,ydim,zdim, &
      q_spat,my_cg,nmine, &
      qm_pot_only )
  use ew_recip,only:frcx
  use constants, only:INV_AMBER_ELECTROSTATIC, zero, one
  use qmmm_module, only : qmewald, qmmm_struct,qmmm_nml
  use ew_bspline,only:kbot,ktop,jbot,jtop
  use fft,only: get_fft_limits,YZ_X_PARTITION
  
   implicit none
   
#  include "box.h"
#  include "ew_parallel.h"
#  include "parallel.h"
   include 'mpif.h'

   integer,intent(in) :: natom,order,nfft1,nfft2,nfft3,nmine,xdim,ydim,zdim
   integer,intent(in),dimension(nmine) :: my_cg
   logical,intent(in) :: qm_pot_only

   double complex,intent(in),dimension(*) :: q_spat

   _REAL_,intent(in)                       :: recip(3,3)
   _REAL_,intent(in),dimension(natom)      ::fr1,fr2,fr3,charge
   _REAL_,intent(inout)                    :: frc(3,natom)
   _REAL_,intent(in),dimension(order,natom) :: theta1,theta2,theta3,dtheta1, &
         dtheta2,dtheta3
   !--------- local vars ----------------------------------------
   integer :: iqk,iqj,j00,i00,kq
   integer :: n,nn,ith1,ith2,ith3,i0,j0,k0,i,j,jq,k,im
   integer :: qm_atom_index
   integer :: iq0,iqi
   _REAL_  :: recip11,recip22,recip33
   _REAL_  :: f1,f2,f3,term,chargen
   _REAL_  :: dfx,dfy,dfz,dnfft1,dnfft2,dnfft3
   _REAL_  :: f1fac,f2fac,f3fac,q
   integer :: xgridmin,xgridmax,ygridmin,ygridmax,zgridmin,zgridmax

   dnfft1 = real(nfft1,8)
   dnfft2 = real(nfft2,8)
   dnfft3 = real(nfft3,8)
   
   recip11 = recip(1,1)*dnfft1
   recip22 = recip(2,2)*dnfft2
   recip33 = recip(3,3)*dnfft3


   call get_fft_limits(YZ_X_PARTITION,xgridmin,xgridmax,ygridmin,ygridmax, &
         zgridmin,zgridmax,mytaskid)
   if(qm_pot_only)then
!  ===============  Below for qm grad sum yielding recip potential ======
      qmewald%mmpot(1:qmmm_struct%nquant)=ZERO      
      qm_atom_index = 0
      do im = 1,nmine
         n = my_cg(im)
         if (qmmm_struct%atom_mask(n)) then
           qm_atom_index = qm_atom_index+1
           f1 = ZERO
           i00 = int(fr1(im)) - order
           j00 = int(fr2(im)) - order
           k0 = int(fr3(im)) - order
           do ith3 = 1,order
              k0 = k0 + 1
              if (k0 >= 0) then
                 k = k0 + 1
              else
                 k = k0 + 1 + nfft3
              end if
            
              if ( k >= kbot .and. k <= ktop )then
                 kq = k - zgridmin
                 iqk = (kq-1)*2*xdim*ydim
                 j0 = j00
                 do ith2 = 1,order
                    j0 = j0 + 1
                    if (j0 >= 0) then
                       j = j0 + 1
                    else
                       j = j0 + 1 + nfft2
                    end if
                    if ( j >= jbot .and. j <= jtop ) then
                       jq = j - ygridmin
                       iqj = iqk + (jq-1)*2*xdim
                       i0 = i00 + 1
                       f1fac =  theta2(ith2,im) *  theta3(ith3,im)
                       do ith1 = 1,order
                          i0 = i0 + 1
                          if (i0 >= 1) then
                          else
                          end if
                          f1 = f1 + q * theta1(ith1,im) * f1fac
                       end do
                    endif
                 end do
              end if  ! ( k >= kbot .and. k <= ktop )
           end do  !  ith3 = 1,order
           do nn=1,qmmm_struct%nquant
              if(qmmm_struct%iqmatoms(nn) == n) &
                    qmewald%mmpot(nn) = f1*INV_AMBER_ELECTROSTATIC
           enddo
        end if !if (qmmm_struct%atom_mask(n)) then     
     end do  !  im = 1,nmine
     !  ===============  Below for nornmal grad sum yielding forces ======
     
     
  else
     do im = 1,nmine
        n = my_cg(im)
        f1 = ZERO
        f2 = ZERO
        f3 = ZERO
        i00 = int(fr1(im)) - order
        j00 = int(fr2(im)) - order
        k0 = int(fr3(im)) - order
        chargen = charge(n)
        do ith3 = 1,order
           k0 = k0 + 1
           k=mod(k0+nfft3,nfft3)+1
           if ( k >= kbot .and. k <= ktop )then
              kq = k - zgridmin
              iqk = (kq-1)*2*xdim*ydim
              j0 = j00
              do ith2 = 1,order
                 j0 = j0 + 1
                 if (j0 >= 0) then
                    j = j0 + 1
                 else
                    j = j0 + 1 + nfft2
                 end if
                 if ( j >= jbot .and. j <= jtop ) then
                    jq = j - ygridmin
                    iqj = iqk + (jq-1)*2*xdim
                    i0 = i00 + 1
                    f1fac =  theta2(ith2,im) *  theta3(ith3,im) * chargen
                    f2fac = dtheta2(ith2,im) *  theta3(ith3,im) * chargen
                    f3fac =  theta2(ith2,im) * dtheta3(ith3,im) * chargen
                    do ith1 = 1,order
                       i0 = i0 + 1
                       if (i0 >= 1) then
                          iq0 = i0 +iqj
                          iqi = ishft(iq0+1,-1)  ! divide 2 get index into

                       else
                          iq0 = i0 + nfft1 + iqj
                          iqi = ishft(iq0+1,-1)
                       end if
                       q = mod(iq0,2)*real(q_spat(iqi)) + &
                             (one - mod(iq0,2))*aimag(q_spat(iqi))
                       !---force is negative of grad
                       f1 = f1 - q * dtheta1(ith1,im) * f1fac
                       f2 = f2 - q *  theta1(ith1,im) * f2fac
                       f3 = f3 - q *  theta1(ith1,im) * f3fac
                    end do
                 endif
              end do
           end if  ! ( k >= kbot .and. k <= ktop )
        end do  !  ith3 = 1,order
        
        if ( ifbox == 1 ) then  ! orthogonal unit cell
           dfx = recip11 * f1
           dfy = recip22 * f2
           dfz = recip33 * f3
        else
           f1 = dnfft1*f1
           f2 = dnfft2*f2
           f3 = dnfft3*f3
           dfx=recip(1,1)*f1+recip(1,2)*f2+recip(1,3)*f3
           dfy=recip(2,1)*f1+recip(2,2)*f2+recip(2,3)*f3
           dfz=recip(3,1)*f1+recip(3,2)*f2+recip(3,3)*f3
        end if
        frc(1,n) = frc(1,n) + dfx
        frc(2,n) = frc(2,n) + dfy
        frc(3,n) = frc(3,n) + dfz
        frcx(1)=frcx(1)+dfx
        frcx(2)=frcx(2)+dfy
        frcx(3)=frcx(3)+dfz
        
     end do  !  im = 1,nmine
  endif   !  qm_pot_only
  return
end subroutine spatial_grad_sumrc 


!=================================================================
!     --- SPATIAL_SCALAR_SUM RC---
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dummy_scalar_sumrc here]
subroutine spatial_scalar_sumrc( &
      q_spat, &
      ewaldcof,volume,recip,prefac1,prefac2,prefac3, &
      nfft1,nfft2,nfft3,xdim,ydim, &
      eer,rec_vir)
  use fft,only:get_fft_limits,XY_Z_PARTITION
   use constants, only : pi, PI2, ZERO,HALF
   implicit none
#  include "extra.h"
#  include "box.h"
#  include "md.h"
#  include "ew_parallel.h"
#  include "parallel.h"
   integer ier,i
   integer nfft1,nfft2,nfft3,xdim,ydim
   _REAL_ prefac1(nfft1),prefac2(nfft2), &
         prefac3(nfft3),ewaldcof,volume
   _REAL_ eer,rec_vir(3,3)
   _REAL_ recip(3,3)
   integer k2q
   _REAL_ q(2,nfft3,xdim,ydim)
   double complex,intent(inout) :: q_spat(nfft3,xdim,ydim)
   _REAL_ fac,denom,eterm,vterm,energy
   integer k1,k2,k3,m1,m2,m3,nff,indtop
   integer nf1,nf2,nf3
   integer k10
   _REAL_ mhat1,mhat2,mhat3,msq,struc2,msq_inv,piv_inv
   integer k1s,k2s,k3s,m1s,m2s,m3s
   _REAL_ mhat1s,mhat2s,mhat3s,msqs
   _REAL_ struc2s,eterms,vterms,denoms,tmp1,tmp2
   integer :: xgmin,xgmax,ygmin,ygmax,zgmin,zgmax
   
   indtop = nfft1*nfft2*nfft3
   piv_inv = 1.d0/(pi*volume)
   fac = PI2/(ewaldcof*ewaldcof)
   nff = nfft1*nfft2
   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2 )nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3 )nf3 = nf3+1
   energy = 0.d0
   k10 = 1
#ifndef noVIRIAL
         rec_vir(m1,m2) = ZERO
#endif
   !........Insist that Q(1,1,1,1) is set to 0 (true already for neutral)
   
   if(master)then
      q(1,1,1,1) = 0.d0
      q(2,1,1,1) = 0.d0
   end if

   !======================================================================
   !        BIG LOOP
   !======================================================================
   call get_fft_limits(XY_Z_PARTITION,zgmin,zgmax,xgmin,xgmax,ygmin,ygmax, &
         mytaskid)
   do k2q = 1,ygmax-ygmin+1
      k2 = k2q + ygmin
      k2s=mod(nfft2-k2+1,nfft2)+1
      m2 = k2 - 1
      if ( k2 > nf2 )m2 = k2 - 1 - nfft2
      
      do k3 = 1,nfft3
         k3s=mod(nfft3-k3+1,nfft3)+1
         m3 = k3 - 1
         if ( k3 > nf3 )m3 = k3 - 1 - nfft3
         k10 = xgmin+1
         if(master)then
            if(k3+k2 == 2) k10 = 2
         end if
         do k1 = k10, xgmax+1
            k1s=nfft1-k1+2
            m1 = k1 - 1
            if ( k1 > nf1 ) m1 = k1 - 1 - nfft1
            mhat1 = recip(1,1)*m1+recip(1,2)*m2+recip(1,3)*m3
            mhat2 = recip(2,1)*m1+recip(2,2)*m2+recip(2,3)*m3
            mhat3 = recip(3,1)*m1+recip(3,2)*m2+recip(3,3)*m3
            msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
            msq_inv = 1.d0/msq
            eterm = exp(-fac*msq)*prefac1(k1)*prefac2(k2)*prefac3(k3) &
                  *piv_inv*msq_inv
            struc2 = q(1,k3,k1-xgmin,k2q)*q(1,k3,k1-xgmin,k2q) + &
                  q(2,k3,k1-xgmin,k2q)*q(2,k3,k1-xgmin,k2q)
            tmp1 = eterm*struc2
            energy = energy + tmp1
#ifndef noVIRIAL
            vterm = 2.d0*(fac*msq + 1.d0)*msq_inv
            tmp2 = tmp1*vterm
            rec_vir(1,1) = rec_vir(1,1) + tmp1 * (vterm*mhat1*mhat1 - 1.d0)
            rec_vir(1,2) = rec_vir(1,2) + tmp2*mhat1*mhat2
            rec_vir(1,3) = rec_vir(1,3) + tmp2*mhat1*mhat3
            rec_vir(2,1) = rec_vir(2,1) + tmp2*mhat2*mhat1
            rec_vir(2,2) = rec_vir(2,2) + tmp1 * (vterm*mhat2*mhat2 - 1.d0)
            rec_vir(2,3) = rec_vir(2,3) + tmp2*mhat2*mhat3
            rec_vir(3,1) = rec_vir(3,1) + tmp2*mhat3*mhat1
            rec_vir(3,2) = rec_vir(3,2) + tmp2*mhat3*mhat2
            rec_vir(3,3) = rec_vir(3,3) + tmp1 * (vterm*mhat3*mhat3 - 1.d0)
#endif
            
            
            if( k1 > 1 .and. k1 <= nfft1 )then
               m1s = k1s - 1
               if ( k1s > nf1 )m1s = k1s - 1 - nfft1
               m2s = k2s - 1
               if ( k2s > nf2 )m2s = k2s - 1 - nfft2
               m3s = k3s - 1
               if ( k3s > nf3 )m3s = k3s - 1 - nfft3
               mhat1s = recip(1,1)*m1s+recip(1,2)*m2s+recip(1,3)*m3s
               mhat2s = recip(2,1)*m1s+recip(2,2)*m2s+recip(2,3)*m3s
               mhat3s = recip(3,1)*m1s+recip(3,2)*m2s+recip(3,3)*m3s
               msqs = mhat1s*mhat1s+mhat2s*mhat2s+mhat3s*mhat3s
               msq_inv = 1.d0/msqs
               eterms = exp(-fac*msqs)*prefac1(k1s)*prefac2(k2s)* &
                     prefac3(k3s)*piv_inv*msq_inv
               tmp1 = eterms*struc2
               energy = energy + tmp1
#ifndef noVIRIAL
               vterms = 2.d0*(fac*msqs + 1.d0)*msq_inv
               tmp2 = tmp1*vterms
               rec_vir(1,1) = rec_vir(1,1) + tmp1*(vterms*mhat1s*mhat1s - 1.d0)
               rec_vir(1,2) = rec_vir(1,2) + tmp2*mhat1s*mhat2s
               rec_vir(1,3) = rec_vir(1,3) + tmp2*mhat1s*mhat3s
               rec_vir(2,1) = rec_vir(2,1) + tmp2*mhat2s*mhat1s
               rec_vir(2,2) = rec_vir(2,2) + tmp1*(vterms*mhat2s*mhat2s - 1.d0)
               rec_vir(2,3) = rec_vir(2,3) + tmp2*mhat2s*mhat3s
               rec_vir(3,1) = rec_vir(3,1) + tmp2*mhat3s*mhat1s
               rec_vir(3,2) = rec_vir(3,2) + tmp2*mhat3s*mhat2s
               rec_vir(3,3) = rec_vir(3,3) + tmp1*(vterms*mhat3s*mhat3s - 1.d0)
#endif
            end if

            q(1,k3,k1-xgmin,k2q) = eterm * q(1,k3,k1-xgmin,k2q)
            q(2,k3,k1-xgmin,k2q) = eterm * q(2,k3,k1-xgmin,k2q)
         end do  !  k1 = k10, nf1+1
      end do  !  k3 = 1,nfft3
   end do  !  k2q = 1, mxzslabs

   eer = HALF * energy

#ifndef noVIRIAL
   do m2 = 1,3
      do m1 = 1,3
         rec_vir(m1,m2) = 0.5d0*rec_vir(m1,m2)
      end do
   end do
#endif

   return
end subroutine spatial_scalar_sumrc 





!=================================================================
!     --- SPATIAL_SCALAR_SUM RC---  for ORTHOGONAL CELL
!-------------------------------------------------------------------

subroutine spatial_scalar_sumrc_orthog( &
      q,ewaldcof,volume,recip, &
      prefac1,prefac2,prefac3, &
      nfft1,nfft2,nfft3,xdim,ydim, &
      eer,rec_vir)
   use constants, only : PI, PI2,ZERO,HALF
   use ew_recip,only:first_pme
   use fft,only:get_xy_z_partition_limits,get_fft_limits,XY_Z_PARTITION
   
   implicit none
   
#  include "extra.h"
#  include "box.h"
#  include "md.h"

#  include "ew_parallel.h"
#  include "parallel.h"

   integer,intent(in) :: nfft1,nfft2,nfft3,xdim,ydim

   double complex,intent(inout) :: q(nfft3,xdim,ydim)
   
   _REAL_,intent(in) :: prefac1(nfft1),prefac2(nfft2), &
         prefac3(nfft3)
   _REAL_,intent(in) :: recip(3,3),ewaldcof,volume
   _REAL_,intent(out) :: eer,rec_vir(3,3)
   
   
   _REAL_ fac,denom,eterm,vterm,energy
   _REAL_ mhat1,mhat2,mhat3,msq,struc2,msq_inv,piv_inv
   _REAL_ mhat1s,mhat2s,mhat3s,msqs
   _REAL_ struc2s,eterms,vterms,denoms
   _REAL_ tmp1,tmp2,m2_m3_tbl,m2_m3_tbls
   _REAL_ recip11,recip22,recip33
   _REAL_ recip11sq,recip22sq,recip33sq
   integer nf1,nf2,nf3
   integer k10,k1s,k2s,k3s,m1s,m2s,m3s
   integer k2q,k1, k2, k3, m1, m2, m3, nff,indtop
   integer nzmax,nymax
   integer ier,i,itot
   integer :: xgmin,xgmax,ygmin,ygmax,zgmin,zgmax
   
   indtop = nfft1*nfft2*nfft3
   piv_inv = 1.d0/(pi*volume)
   fac = PI2/(ewaldcof*ewaldcof)
   nff = nfft1*nfft2
   nf1 = nfft1/2
   if ( 2*nf1 < nfft1 )nf1 = nf1+1
   nf2 = nfft2/2
   if ( 2*nf2 < nfft2 )nf2 = nf2+1
   nf3 = nfft3/2
   if ( 2*nf3 < nfft3 )nf3 = nf3+1
   recip11 = recip(1,1)
   recip22 = recip(2,2)
   recip33 = recip(3,3)

   if( first_pme) then
      allocate(m1_tbl(-(nfft1/2 + 1) : nfft1/2 + 1), &
            m2_tbl(-(nfft2/2 + 1) : nfft2/2 + 1), &
            m3_tbl(-(nfft3/2 + 1) : nfft3/2 + 1), &
            stat = ier)
      REQUIRE( ier==0 )
   end if
   
   if( ntp > 0 .or. first_pme) then
      recip11sq = recip11*recip11
      recip22sq = recip22*recip22
      recip33sq = recip33*recip33
      itot = 0
      do i = -(nfft1/2 + 1), nfft1/2 + 1
         itot = itot + 1
         m1_tbl(i) = -fac * dble(i) * dble(i) * recip11sq
      end do
      call vdexp(itot,m1_tbl(-nfft1/2-1),m1_tbl(-nfft1/2-1))

      itot = 0
      do i = -(nfft2/2 + 1), nfft2/2 + 1
         itot = itot + 1
         m2_tbl(i) = -fac * dble(i) * dble(i) * recip22sq
      end do
      call vdexp(itot,m2_tbl(-nfft2/2-1),m2_tbl(-nfft2/2-1))

      itot = 0
      do i = -(nfft3/2 + 1), nfft3/2 + 1
         itot = itot + 1
         m3_tbl(i) = -fac * dble(i) * dble(i) * recip33sq
      end do
      call vdexp(itot,m3_tbl(-nfft3/2-1),m3_tbl(-nfft3/2-1))

      first_pme = .false.
   end if
   energy = ZERO
   k10 = 1
#ifndef noVIRIAL
   rec_vir = ZERO
#endif
   
   !........Insist that Q(1,1,1) is set to 0 (true already for neutral)
   
   if(master)then
      q(1,1,1) = (ZERO,ZERO)
   end if

   !======================================================================
   !        BIG LOOP
   !======================================================================
   call get_fft_limits(XY_Z_PARTITION,zgmin,zgmax,xgmin,xgmax,ygmin,ygmax, &
         mytaskid)
   do k2q = 1, ygmax-ygmin+1
      if(master)then
         k2=k2q+ygmin
         k2s=mod(nfft2-k2+1,nfft2)+1
      else
         k2 = k2q + ygmin
         k2s=mod(nfft2-k2+1,nfft2)+1
      end if
      m2 = k2 - 1
      if ( k2 > nf2 )m2 = k2 - 1 - nfft2
      mhat2 = recip22*m2
      m2s = k2s - 1
      if ( k2s > nf2 )m2s = k2s - 1 - nfft2
      mhat2s = recip22*m2s
      
      do k3 = zgmin+1,zgmax+1
         k3s=mod(nfft3-k3+1,nfft3)+1
         m3 = k3 - 1
         if ( k3 > nf3 )m3 = k3 - 1 - nfft3
         mhat3 = recip33*m3
         m3s = k3s - 1
         if ( k3s > nf3 )m3s = k3s - 1 - nfft3
         mhat3s = recip33*m3s
         
         k10 = xgmin+1
         if(master .and. (k3+k2 == 2) )k10 = 2
         m2_m3_tbl = piv_inv*m2_tbl(m2)*m3_tbl(m3) &
               *prefac2(k2)*prefac3(k3)
         m2_m3_tbls = piv_inv*m2_tbl(m2s)*m3_tbl(m3s) &
               *prefac2(k2s)*prefac3(k3s)
         
         do k1 = k10, xgmax+1
            k1s=nfft1-k1+2
            m1 = k1 - 1
            if ( k1 > nf1 ) m1 = k1 - 1 - nfft1
            mhat1 = recip11*m1
            msq = mhat1*mhat1+mhat2*mhat2+mhat3*mhat3
            msq_inv = 1.d0/msq
            eterm = m1_tbl(m1) * m2_m3_tbl * &
                  prefac1(k1) * msq_inv
            struc2 = real(q(k3,k1-xgmin,k2q))*real(q(k3,k1-xgmin,k2q)) + &
                  aimag(q(k3,k1-xgmin,k2q))*aimag(q(k3,k1-xgmin,k2q))
            tmp1 = eterm*struc2
            energy = energy + tmp1
#ifndef noVIRIAL
            vterm = 2.d0*(fac*msq + 1.d0)*msq_inv
            tmp2 = tmp1*vterm
            rec_vir(1,1) = rec_vir(1,1) + tmp1 * (vterm*mhat1*mhat1 - 1.d0)
            rec_vir(1,2) = rec_vir(1,2) + tmp2*mhat1*mhat2
            rec_vir(1,3) = rec_vir(1,3) + tmp2*mhat1*mhat3
            rec_vir(2,1) = rec_vir(2,1) + tmp2*mhat2*mhat1
            rec_vir(2,2) = rec_vir(2,2) + tmp1 * (vterm*mhat2*mhat2 - 1.d0)
            rec_vir(2,3) = rec_vir(2,3) + tmp2*mhat2*mhat3
            rec_vir(3,1) = rec_vir(3,1) + tmp2*mhat3*mhat1
            rec_vir(3,2) = rec_vir(3,2) + tmp2*mhat3*mhat2
            rec_vir(3,3) = rec_vir(3,3) + tmp1 * (vterm*mhat3*mhat3 - 1.d0)
#endif
            
            
            if( k1 > 1 .and. k1 <= nfft1 )then
               m1s = k1s - 1
               if ( k1s > nf1 )m1s = k1s - 1 - nfft1
               mhat1s = recip11*m1s
               msqs = mhat1s*mhat1s+mhat2s*mhat2s+mhat3s*mhat3s
               msq_inv = 1.d0/msqs
               eterms = m1_tbl(m1s) * m2_m3_tbls * &
                     prefac1(k1s) * msq_inv
               tmp1 = eterms*struc2
               energy = energy + tmp1
#ifndef noVIRIAL
               vterms = 2.d0*(fac*msqs + 1.d0)*msq_inv
               tmp2 = tmp1*vterms
               rec_vir(1,1) = rec_vir(1,1) + tmp1*(vterms*mhat1s*mhat1s - 1.d0)
               rec_vir(1,2) = rec_vir(1,2) + tmp2*mhat1s*mhat2s
               rec_vir(1,3) = rec_vir(1,3) + tmp2*mhat1s*mhat3s
               rec_vir(2,1) = rec_vir(2,1) + tmp2*mhat2s*mhat1s
               rec_vir(2,2) = rec_vir(2,2) + tmp1*(vterms*mhat2s*mhat2s - 1.d0)
               rec_vir(2,3) = rec_vir(2,3) + tmp2*mhat2s*mhat3s
               rec_vir(3,1) = rec_vir(3,1) + tmp2*mhat3s*mhat1s
               rec_vir(3,2) = rec_vir(3,2) + tmp2*mhat3s*mhat2s
               rec_vir(3,3) = rec_vir(3,3) + tmp1*(vterms*mhat3s*mhat3s - 1.d0)
#endif
            end if

            q(k3,k1-xgmin,k2q) = eterm * q(k3,k1-xgmin,k2q)
         end do  !  k1 = k10, nf1+1
      end do  !  k3 = 1,nfft3
   end do  !  k2q = 1, nymax

   eer = HALF * energy
#ifndef noVIRIAL
   rec_vir = HALF*rec_vir
#endif

   return

end subroutine spatial_scalar_sumrc_orthog 

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!         SPATIAL GET_GRID_WEIGHTS
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine spatial_get_grid_weights( &
      natom,crd,recip,nfft1,nfft2,nfft3, &
      fr1,fr2,fr3,order,theta1,theta2,theta3, &
      dtheta1,dtheta2,dtheta3, &
      d2theta1,d2theta2,d2theta3, &
      my_cg,nmine,nderiv,num_ks)
   use fft,only:get_fft_limits,YZ_X_PARTITION
   use ew_bspline,only:kbot,ktop,jbot,jtop, &
         fill_bspline_0,fill_bspline_1,fill_bspline_2

   implicit none
   integer natom,nfft1,nfft2,nfft3,order
   _REAL_ crd(3,natom),recip(3,3)
   _REAL_ fr1(*),fr2(*),fr3(*)
   _REAL_ theta1(order,*),theta2(order,*), &
         theta3(order,*)
   _REAL_ dtheta1(order,*),dtheta2(order,*), &
         dtheta3(order,*)
   _REAL_ d2theta1(order,*),d2theta2(order,*), &
         d2theta3(order,*)
   integer my_cg(*),nmine,nderiv,num_ks
   _REAL_ fr3n,fr2n,fr1n,w,w1,w2,w3
   _REAL_ anint

   integer :: xgridmin,xgridmax,ygridmin,ygridmax,zgridmin,zgridmax

#  include "parallel.h"
#  include "ew_parallel.h"
   integer kbot0,ktop1,idoz
   integer j00,jbot0,jtop1,ido

   integer n,k00
   integer imine
   
   !     sanity check:
   
   if ( order < 2 + nderiv )then
      write(6,*)'too many B-spline derivs for order! '
      call mexit(6,1)
   end if
   if ( nderiv > 2 )then
      write(6,*)'More than 2 derivs of B-splines not implemented yet!'
      call mexit(6,1)
   end if

   call get_fft_limits(YZ_X_PARTITION,xgridmin,xgridmax,ygridmin,ygridmax, &
         zgridmin,zgridmax,mytaskid)

   kbot0 = zgridmin
   kbot = kbot0 + 1
   ktop = zgridmax+1
   ktop1 = ktop + order - 2

   jbot0 = ygridmin
   jbot = jbot0 + 1
   jtop = ygridmax+1
   jtop1 = jtop + order - 2
   imine = 0
   !  ------------------------------------------------
   !     First filter the atoms and make a list (my_cg)
   !          of atoms needed for generating this part of
   !          the grid
   !     If there are too many atoms for the allotted space,
   !          return when nmine exceedsn um_ks
   
   do n = 1,natom
      w = crd(1,n)*recip(1,3) &
            +crd(2,n)*recip(2,3)+crd(3,n)*recip(3,3)
      fr3n = nfft3*(w - (anint(w) - 0.5d0))
      k00 = int(fr3n)
      ! code for filtering atoms. In single proc mode, do all atoms
      ido  = 0
      idoz = 0
      if ( ktop1 >= nfft3)then
         if ( k00 >= kbot0 .or. k00 <= ktop1 - nfft3 )idoz = 1
      else
         if ( k00 >= kbot0 .and. k00 <= ktop1)idoz = 1
      end if
      if ( idoz == 1)then
         w = crd(1,n)*recip(1,2) &
               +crd(2,n)*recip(2,2)+crd(3,n)*recip(3,2)
         fr2n = nfft2*(w - (anint(w) - 0.5d0))
         j00 = int(fr2n)
         if ( jtop1 >= nfft2)then
            if( j00 >= jbot0 .or. j00 <= jtop1 - nfft2 )ido = 1
         else
            if( j00 >= jbot0 .and. j00 <= jtop1)ido = 1
         end if
      endif
      if ( ido == 1)then
         imine = imine + 1
         !           ----- do not fill my_cg past num_ks ------
         if( imine <= num_ks) my_cg(imine)=n
      end if
   end do
   nmine=imine
   !   ----- ERROR condition met --------------------
   if(imine > num_ks) return
   !   ----------------------------------------------
   imine = 0
   do imine = 1,nmine
      n=my_cg(imine)
      w = crd(1,n)*recip(1,3) &
            +crd(2,n)*recip(2,3)+crd(3,n)*recip(3,3)
      fr3n = nfft3*(w - (anint(w) - 0.5d0))
      k00 = int(fr3n)
      w = crd(1,n)*recip(1,1) &
            +crd(2,n)*recip(2,1)+crd(3,n)*recip(3,1)
      fr1n = nfft1*(w - (anint(w) - 0.5d0))
      w = crd(1,n)*recip(1,2) &
            +crd(2,n)*recip(2,2)+crd(3,n)*recip(3,2)
      fr2n = nfft2*(w - (anint(w) - 0.5d0))
      fr1(imine)=fr1n
      fr2(imine)=fr2n
      fr3(imine)=fr3n
      w1 = fr1n-int(fr1n)
      w2 = fr2n-int(fr2n)
      w3 = fr3n-int(fr3n)
      if ( nderiv == 0 )then
         call fill_bspline_0(w1,order,theta1(1,imine))
         call fill_bspline_0(w2,order,theta2(1,imine))
         call fill_bspline_0(w3,order,theta3(1,imine))
      else if ( nderiv == 1 )then
         call fill_bspline_1(w1,order,theta1(1,imine), &
               dtheta1(1,imine))
         call fill_bspline_1(w2,order,theta2(1,imine), &
               dtheta2(1,imine))
         call fill_bspline_1(w3,order,theta3(1,imine), &
               dtheta3(1,imine))
      else if ( nderiv == 2 )then
         call fill_bspline_2(w1,order,theta1(1,imine), &
               dtheta1(1,imine),d2theta1(1,imine))
         call fill_bspline_2(w2,order,theta2(1,imine), &
               dtheta2(1,imine),d2theta2(1,imine))
         call fill_bspline_2(w3,order,theta3(1,imine), &
               dtheta3(1,imine),d2theta3(1,imine))
      end if
   end do  !  imine = 1,nmine
   return
end subroutine spatial_get_grid_weights
#endif /* MPI */
end module ew_recip_spatial
