#include "dprec.fh"
#include "assert.fh"

!------------------------------------------------------------
module amoeba_vdw
implicit none
private

integer,save :: do_flag

integer,save,allocatable :: vdw_atom_parent(:),vdw_atom_type(:)
integer,save  :: num_vdw_atom_types
_REAL_,save,allocatable :: vdw_atom_parent_crd_wt(:)
_REAL_,save,allocatable :: vdw_epsilon(:,:),vdw_radius(:,:)
_REAL_,save :: vdw_buf_gamma,vdw_buf_delta  !halgren's vdw buffered by these
integer,save,allocatable :: vdw_type_count(:)
_REAL_, save :: ene_vdw_longrange_factor,vir_vdw_longrange_factor
_REAL_,save :: vdw_switch_on,vdw_switch_on_2, &
               vdw_switch_off,vdw_switch_off_2, &
               c0,c1,c2,c3,c4,c5


public AM_VDW_read_parm,AM_VDW_ADJUST_ene_frc, &
       AM_VDW_longrange_ene,AM_VDW_set_user_bit, &
       AM_VDW_DIRECT_ene_frc_i,AM_VDW_deallocate
#ifdef MPI 
public AM_VDW_bcast
#endif
contains
!------------------------------------------------------------
subroutine AM_VDW_set_user_bit(do_this)
  integer,intent(in) :: do_this
#include "do_flag.h"

  if ( do_this == 1 )then
    do_flag = ibset(do_flag,USER_BIT)
  else
    do_flag = ibclr(do_flag,USER_BIT)
  endif
end subroutine AM_VDW_set_user_bit
!------------------------------------------------------------

#ifdef MPI 
subroutine AM_VDW_bcast(natom)
   implicit none
   integer natom
   integer ierr

#  include "box.h"
   include 'mpif.h'
#  include "extra.h"
#  include "parallel.h"
   call mpi_bcast(do_flag,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(num_vdw_atom_types,1,MPI_INTEGER,0,commsander,ierr)
  
   if(.not.master) then
      allocate(vdw_atom_parent(natom),stat=ierr)
      REQUIRE(ierr==0)
      allocate(vdw_atom_type(natom),stat=ierr)
      REQUIRE(ierr==0)
      allocate(vdw_atom_parent_crd_wt(natom),stat=ierr)
      REQUIRE(ierr==0)
      allocate(vdw_epsilon(num_vdw_atom_types,num_vdw_atom_types),stat=ierr)
      REQUIRE(ierr==0)
      allocate(vdw_radius(num_vdw_atom_types,num_vdw_atom_types),stat=ierr)
      REQUIRE(ierr==0)
   end if

   call mpi_bcast(vdw_atom_parent,natom,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(vdw_atom_type  ,natom,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(vdw_atom_parent_crd_wt,natom,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(vdw_buf_gamma, 1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(vdw_buf_delta, 1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(vdw_epsilon,num_vdw_atom_types*num_vdw_atom_types,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(vdw_radius, num_vdw_atom_types*num_vdw_atom_types,MPI_DOUBLE_PRECISION,0,commsander,ierr)
  
   call mpi_bcast(ene_vdw_longrange_factor,1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(vir_vdw_longrange_factor,1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(vdw_switch_on,1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(vdw_switch_on_2, 1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(vdw_switch_off,1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(vdw_switch_off_2, 1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(c0,1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(c1,1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(c2,1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(c3,1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(c4,1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(c5,1, MPI_DOUBLE_PRECISION,0,commsander,ierr)
end subroutine AM_VDW_bcast
#endif

function AM_VDW_read_parm(nf)
  integer :: AM_VDW_read_parm
  integer, intent(in) :: nf
#include "do_flag.h"
  integer num_list,num_atoms,ier,dim1

  AM_VDW_read_parm = 0
  call AMOEBA_get_numlist('AMOEBA_VDW_ATOM_TYPES_',nf,num_list)
  if ( num_list <= 0 )then
    do_flag = ibclr(do_flag,VALID_BIT)
    return
  endif
  num_atoms = num_list
  allocate(vdw_atom_parent(num_atoms),stat=ier)
  REQUIRE(ier==0)
  allocate(vdw_atom_type(num_atoms),stat=ier)
  REQUIRE(ier==0)
  allocate(vdw_atom_parent_crd_wt(num_atoms),stat=ier)
  REQUIRE(ier==0)
  dim1 = 1
  call AMOEBA_read_list_data('AMOEBA_VDW_ATOM_TYPES_',nf, &
                  dim1,num_atoms,vdw_atom_type)
  call AMOEBA_read_list_data('AMOEBA_VDW_ATOM_PARENT_',nf, &
                  dim1,num_atoms,vdw_atom_parent)
  call AMOEBA_read_real_list_data('AMOEBA_VDW_PARENT_COORD_WEIGHT_',nf, &
                  dim1,num_atoms,vdw_atom_parent_crd_wt)
  call AMOEBA_read_real_scalar('AMOEBA_VDW_BUFFER_DELTA',nf,vdw_buf_delta)
  call AMOEBA_read_real_scalar('AMOEBA_VDW_BUFFER_GAMMA',nf,vdw_buf_gamma)
  call AMOEBA_get_numlist('AMOEBA_VDW_PARAMS_',nf,num_vdw_atom_types)
  if ( num_vdw_atom_types <= 0 )then
    write(6,*)'No VDW parameters present!'
    call mexit(6,1)
  endif
  allocate(vdw_epsilon(num_vdw_atom_types,num_vdw_atom_types),stat=ier)
  REQUIRE(ier==0)
  allocate(vdw_radius(num_vdw_atom_types,num_vdw_atom_types),stat=ier)
  REQUIRE(ier==0)
  dim1 = num_vdw_atom_types
  call AMOEBA_read_real_list_data('AMOEBA_VDW_MIXED_RADII_',nf, &
                  dim1,num_vdw_atom_types,vdw_radius)
  call AMOEBA_read_real_list_data('AMOEBA_VDW_MIXED_EPSILONS_',nf, &
                  dim1,num_vdw_atom_types,vdw_epsilon)
  ! calculate the vdw switch function
  call AM_VDW_switch()
  ! calculate long range factor
  call AM_VDW_longrange_factor(num_atoms)
  AM_VDW_read_parm = 1
  do_flag = ibset(do_flag,VALID_BIT)
end function AM_VDW_read_parm
!------------------------------------------------------------
subroutine AM_VDW_deallocate()
   if (allocated(vdw_atom_parent))deallocate(vdw_atom_parent)
   if (allocated(vdw_atom_parent_crd_wt))deallocate(vdw_atom_parent_crd_wt)
   if (allocated(vdw_epsilon))deallocate(vdw_epsilon)
   if (allocated(vdw_radius))deallocate(vdw_radius)
   if (allocated(vdw_type_count))deallocate(vdw_type_count)
end subroutine AM_VDW_deallocate
!------------------------------------------------------------
subroutine AM_VDW_switch()
use amoeba_mdin, only : vdw_taper
# include "box.h"
  _REAL_ denom

  if ( vdw_taper >= 1.d0 .or. vdw_taper <= 0.0 )then
     write(6,*)'Bad value for vdw_taper: ',vdw_taper
     write(6,*)'should be fraction'
     call mexit(6,1)
  endif
  vdw_switch_on = vdw_taper*cut
  vdw_switch_on_2 = vdw_switch_on*vdw_switch_on
  vdw_switch_off = cut
  vdw_switch_off_2 = cut*cut
  
  denom = (vdw_switch_off - vdw_switch_on)**5
  c0 = vdw_switch_off*vdw_switch_off_2 * &
         (vdw_switch_off_2-5.0d0*vdw_switch_off*vdw_switch_on+ &
          10.0d0*vdw_switch_on_2) / denom
  c1 = -30.0d0 * vdw_switch_on_2*vdw_switch_off_2 / denom
  c2 = 30.0d0 * (vdw_switch_off_2*vdw_switch_on+ &
                 vdw_switch_off*vdw_switch_on_2) / denom
  c3 = -10.0d0 * (vdw_switch_off_2+4.0d0*vdw_switch_off*vdw_switch_on+ &
                  vdw_switch_on_2) / denom
  c4 = 15.0d0 * (vdw_switch_off+vdw_switch_on) / denom
  c5 = -6.0d0 / denom

end subroutine AM_VDW_switch
!------------------------------------------------------------
subroutine AM_VDW_longrange_factor(num_atoms)
  use amoeba_mdin, only : do_vdw_taper, &
                         soft_lambda,soft_atom_range1,soft_atom_range2,&
                         soft_line,vdw_longrange_lambda,AMOEBA_read_soft
  use constants, only : zero,one,two,three,four,five,seven,pi
  integer,intent(in) :: num_atoms

  integer ier,n,nt,kdel,ndel,i,j,i_range,i_soft,i_soft_type
  integer,save,allocatable :: lig_type_ct(:)
  _REAL_ :: r,r1,r2,r3,r4,r5,req,eps,f,sume,sumv,t1,t2,rho,switch,delr
  _REAL_ :: dt1drho,dt2drho,drhodr,dfdr,dswitch_dr,g1,g2
  allocate(vdw_type_count(num_vdw_atom_types),stat=ier)
  REQUIRE(ier==0)
  allocate(lig_type_ct(num_vdw_atom_types),stat=ier)
  REQUIRE(ier==0)

  if (vdw_longrange_lambda .ne. 1.0 .or. soft_lambda .ne. 1.0) then
     call AMOEBA_read_soft()
  endif
  do n = 1,num_vdw_atom_types
    vdw_type_count(n) = 0
    lig_type_ct(n) = 0
  enddo
  do n = 1,num_atoms
    nt = vdw_atom_type(n)
    vdw_type_count(nt) = vdw_type_count(nt) + 1
    do i_range = 1, soft_line
      if (n .ge. soft_atom_range1(i_range).and. n.le.soft_atom_range2(i_range)) then
           lig_type_ct(nt) = lig_type_ct(nt) + 1
      endif
  enddo
  enddo
! choose r2 to be 60 Angstroms or ~20 sigma
! choose delr to be 0.05 angstroms, or ndel = 20*(r2-r1)
  if ( do_vdw_taper == 1 )then
     r1 = vdw_switch_on
  else
     r1 = vdw_switch_off
  endif
  r2 = 60.d0
  ndel = nint(100.d0*(r2 - r1))
  delr = (r2 - r1)/ndel
  ene_vdw_longrange_factor = zero
  vir_vdw_longrange_factor = zero
! use trapezoidal rule for each integral
  do j = 1,num_vdw_atom_types
    do i = 1,num_vdw_atom_types
      req = vdw_radius(i,j)
      eps = vdw_epsilon(i,j)
      sume = zero
      sumv = zero
      do kdel = 1,ndel
        r = r1 + (kdel-1)*delr
        rho = r / req
        t1 = ((one + vdw_buf_delta) / (rho + vdw_buf_delta))**7
        t2 = (one + vdw_buf_gamma) / (rho**7 + vdw_buf_gamma)
        dt1drho = -seven*t1 / (rho + vdw_buf_delta)
        dt2drho = -seven*t2 * (rho**6 / (rho**7 + vdw_buf_gamma))
        drhodr = one / req
        f = eps*t1*(t2 - two)
        dfdr = eps*(dt1drho*(t2 - two) + t1*dt2drho)*drhodr
        if ( r < vdw_switch_off )then
          g1 = f
          g2 = dfdr
          r2 = r*r
          r3 = r2*r
          r4 = r3*r
          r5 = r4*r
          switch = c5*r5 + c4*r4 + c3*r3 + c2*r2 + c1*r + c0
          dswitch_dr = five*c5*r4 + four*c4*r3 + three*c3*r2 + &
                       two*c2*r + c1
          dfdr = switch*dfdr + f*dswitch_dr
          f = switch*f
          f = g1 - f ! need what's not done explicitly
          dfdr = g2 - dfdr  ! need what's not done explicitly
        endif
        if ( kdel .eq. 1 .or. kdel .eq. ndel )then
          sume = sume + 0.5d0*r**2*f*delr
          sumv = sumv + 0.5d0*r**3*dfdr*delr
        else
          sume = sume + r**2*f*delr
          sumv = sumv + r**3*dfdr*delr
        endif
      enddo
      ! note the 2*pi below not 4*pi---since we do each i,j pair 2x
    ene_vdw_longrange_factor = ene_vdw_longrange_factor + two*pi*sume*(vdw_type_count(i)*vdw_type_count(j) &
      -(1-vdw_longrange_lambda)*(lig_type_ct(i)*vdw_type_count(j)+(vdw_type_count(i)-lig_type_ct(i))*lig_type_ct(j)))
    vir_vdw_longrange_factor = vir_vdw_longrange_factor + two*pi*sumv*(vdw_type_count(i)*vdw_type_count(j) &
      -(1-vdw_longrange_lambda)*(lig_type_ct(i)*vdw_type_count(j)+(vdw_type_count(i)-lig_type_ct(i))*lig_type_ct(j)))
    enddo
  enddo
  
end subroutine AM_VDW_longrange_factor
!------------------------------------------------------------
subroutine AM_VDW_longrange_ene(ene_vdw,vir_tensor)
  use nblist, only: volume
  use amoeba_mdin, only : do_vdw_longrange
  _REAL_,intent(out) :: ene_vdw
  _REAL_,intent(inout) :: vir_tensor(3,3)
#  include "do_flag.h"

   ene_vdw = 0.d0
   if ( do_flag /= PROCEED )return
   if ( do_vdw_longrange /= 1 )return

  ene_vdw = ene_vdw_longrange_factor / volume
  vir_tensor(1,1) = vir_tensor(1,1) + vir_vdw_longrange_factor / (3.d0*volume)
  vir_tensor(2,2) = vir_tensor(2,2) + vir_vdw_longrange_factor / (3.d0*volume)
  vir_tensor(3,3) = vir_tensor(3,3) + vir_vdw_longrange_factor / (3.d0*volume)

end subroutine AM_VDW_longrange_ene
!------------------------------------------------------------
subroutine AM_VDW_ADJUST_ene_frc(crd,num_adjust_list,adjust_list,vdw_weight, &
                                 ene_vdw,frc,virial)
   use constants, only : zero,one,two,three,four,five,seven
   use amoeba_mdin, only : do_vdw_taper
   _REAL_,intent(in) :: crd(3,*)
   integer,intent(in) :: adjust_list(3,*),num_adjust_list
   _REAL_,intent(in) :: vdw_weight(9)
   _REAL_,intent(inout) :: ene_vdw,frc(3,*),virial(3,3)

   integer :: i,j,k,ih,jh,it,jt,n
   _REAL_ :: wi,wj,xi,yi,zi,xj,yj,zj,delx,dely,delz,delr2,delr,eps,rad
   _REAL_ :: rho,t1,t2,dt1drho,dt2drho,drhodr,term,dfx,dfy,dfz
   _REAL_ :: vxx,vxy,vxz,vyy,vyz,vzz,prefac
   _REAL_ :: switch,dswitch_dr,f,dfdr,delr3,delr4,delr5
#  include "do_flag.h"

   ene_vdw = zero
   if ( do_flag /= PROCEED )return

   vxx = zero
   vxy = zero
   vxz = zero
   vyy = zero
   vyz = zero
   vzz = zero
   do n = 1,num_adjust_list
      i = adjust_list(1,n)
      j = adjust_list(2,n)
      k = adjust_list(3,n)
      ih = vdw_atom_parent(i)
      jh = vdw_atom_parent(j)
      it = vdw_atom_type(i)
      jt = vdw_atom_type(j)
      eps = vdw_epsilon(it,jt)
      rad = vdw_radius(it,jt)
      wi = vdw_atom_parent_crd_wt(i)
      wj = vdw_atom_parent_crd_wt(j)
      xi = wi*crd(1,ih) + (one-wi)*crd(1,i)
      yi = wi*crd(2,ih) + (one-wi)*crd(2,i)
      zi = wi*crd(3,ih) + (one-wi)*crd(3,i)
      xj = wj*crd(1,jh) + (one-wj)*crd(1,j)
      yj = wj*crd(2,jh) + (one-wj)*crd(2,j)
      zj = wj*crd(3,jh) + (one-wj)*crd(3,j)
      delx = xj - xi
      dely = yj - yi
      delz = zj - zi
      delr2 = delx*delx + dely*dely + delz*delz
      if ( delr2 < vdw_switch_off_2 )then
         delr = sqrt(delr2)
         rho = delr / rad
         t1 = ((one + vdw_buf_delta) / (rho + vdw_buf_delta))**7
         t2 = (one + vdw_buf_gamma) / (rho**7 + vdw_buf_gamma)
         dt1drho = -seven*t1 / (rho + vdw_buf_delta)
         dt2drho = -seven*t2 * (rho**6 / (rho**7 + vdw_buf_gamma))
         drhodr = one / rad
         prefac = vdw_weight(k)*eps
         f = prefac*t1*(t2 - two)
         dfdr = prefac*(dt1drho*(t2 - two) + t1*dt2drho)*drhodr
         if ( do_vdw_taper == 1 .and. delr2 > vdw_switch_on_2 )then
            delr3 = delr2*delr
            delr4 = delr3*delr
            delr5 = delr4*delr
            switch = c5*delr5 + c4*delr4 + c3*delr3 + c2*delr2 + &
                     c1*delr + c0
            dswitch_dr = five*c5*delr4 + four*c4*delr3 + three*c3*delr2 + &
                         two*c2*delr + c1
            f = switch*f
            dfdr = switch*dfdr + f*dswitch_dr
         endif
         ene_vdw = ene_vdw + f
         term = dfdr/delr
         dfx = term*delx
         dfy = term*dely
         dfz = term*delz
         frc(1,i) = frc(1,i) + (one-wi)*dfx
         frc(2,i) = frc(2,i) + (one-wi)*dfy
         frc(3,i) = frc(3,i) + (one-wi)*dfz
         frc(1,ih) = frc(1,ih) + wi*dfx
         frc(2,ih) = frc(2,ih) + wi*dfy
         frc(3,ih) = frc(3,ih) + wi*dfz
         frc(1,j) = frc(1,j) - (one-wj)*dfx
         frc(2,j) = frc(2,j) - (one-wj)*dfy
         frc(3,j) = frc(3,j) - (one-wj)*dfz
         frc(1,jh) = frc(1,jh) - wj*dfx
         frc(2,jh) = frc(2,jh) - wj*dfy
         frc(3,jh) = frc(3,jh) - wj*dfz
         vxx = vxx + dfx*delx
         vxy = vxy + dfx*dely
         vxz = vxz + dfx*delz
         vyy = vyy + dfy*dely
         vyz = vyz + dfy*delz
         vzz = vzz + dfz*delz
      endif
   enddo !n = 1,num_adjust_list
   virial(1,1) = virial(1,1) + vxx
   virial(1,2) = virial(1,2) + vxy
   virial(1,3) = virial(1,3) + vxz
   virial(2,1) = virial(2,1) + vxy
   virial(2,2) = virial(2,2) + vyy
   virial(2,3) = virial(2,3) + vyz
   virial(3,1) = virial(3,1) + vxz
   virial(3,2) = virial(3,2) + vyz
   virial(3,3) = virial(3,3) + vzz
end subroutine AM_VDW_ADJUST_ene_frc
!------------------------------------------------------------
subroutine AM_VDW_DIRECT_ene_frc_i(i,ipairs,numtot,xk,yk,zk, &
                                   crd,ene_vdw,frc,virial)
   use nblist, only: bckptr,imagcrds,tranvec
   use constants, only : zero,one,two,three,four,five,seven
   use amoeba_mdin, only : do_vdw_taper, &
                           soft_lambda,soft_alpha,soft_expo, &
                           soft_atom_range1,soft_atom_range2, soft_line
   integer,intent(in) :: i,ipairs(*),numtot
   _REAL_,intent(in) :: xk,yk,zk,crd(3,*)
   _REAL_,intent(inout) :: ene_vdw,frc(3,*),virial(3,3)

   integer :: itran,it,jt,ih,jh,j,m,np,mask27,i_range
   _REAL_ :: xktran(3,18)
   _REAL_ :: wi,wj,xi,yi,zi,delx,dely,delz,delr2,delr,eps,rad
   _REAL_ :: rho,t1,t2,dt1drho,dt2drho,drhodr,term,dfx,dfy,dfz
   _REAL_ :: rho6,rho7
   _REAL_ :: vxx,vxy,vxz,vyy,vyz,vzz
   _REAL_ :: switch,dswitch_dr,f,dfdr,delr3,delr4,delr5

#  include "do_flag.h"
#  include "parallel.h"

   if ( do_flag /= PROCEED )return

   mask27 = 2**27 - 1

   ih = vdw_atom_parent(i)
   wi = vdw_atom_parent_crd_wt(i)
   xi = wi*(crd(1,ih)-crd(1,i)) + xk 
   yi = wi*(crd(2,ih)-crd(2,i)) + yk
   zi = wi*(crd(3,ih)-crd(3,i)) + zk
   do m=1,18
      xktran(1,m) = tranvec(1,m) - xi
      xktran(2,m) = tranvec(2,m) - yi
      xktran(3,m) = tranvec(3,m) - zi
   end do
   it = vdw_atom_type(i)
   vxx = zero
   vxy = zero
   vxz = zero
   vyy = zero
   vyz = zero
   vzz = zero


   do m = 1,numtot
      np=ipairs(m)
      itran=ishft(np,-27)
      np = iand(np,mask27)
      j = bckptr(np)
      jh = vdw_atom_parent(j)
      wj = vdw_atom_parent_crd_wt(j)
      ! calculate delx etc between vdw points, not the atoms
      delx = wj*(crd(1,jh)-crd(1,j)) + imagcrds(1,np) + xktran(1,itran)
      dely = wj*(crd(2,jh)-crd(2,j)) + imagcrds(2,np) + xktran(2,itran)
      delz = wj*(crd(3,jh)-crd(3,j)) + imagcrds(3,np) + xktran(3,itran)
      delr2 = delx*delx + dely*dely+delz*delz
      if ( delr2 < vdw_switch_off_2 )then
         jt = vdw_atom_type(j)
         eps = vdw_epsilon(jt,it)
         rad = vdw_radius(jt,it)
         delr = sqrt(delr2)
         rho = delr / rad
         rho6 = rho**6
         rho7 = rho6*rho
         t1 = ((one + vdw_buf_delta) / (rho + vdw_buf_delta))**7
         t2 = (one + vdw_buf_gamma) / (rho7 + vdw_buf_gamma)
         dt1drho = -seven*t1 / (rho + vdw_buf_delta)
         dt2drho = -seven*t2 * (rho6 / (rho7 + vdw_buf_gamma))
         drhodr = one / rad
         do i_range = 1,soft_line
          if (       (i.ge.soft_atom_range1(i_range).and. &
                      i.le.soft_atom_range2(i_range)     )&
               .and.                                      &
                     (j.ge.soft_atom_range1(i_range).and. &
                      j.le.soft_atom_range2(i_range)     ))then
             continue
          else if ( .not. (i.ge.soft_atom_range1(i_range).and. &
                           i.le.soft_atom_range2(i_range)     )&
                    .and.                                      &
                    .not. (j.ge.soft_atom_range1(i_range).and. &
                           j.le.soft_atom_range2(i_range)     ))then
             continue
          else
                eps = eps * soft_lambda ** soft_expo
                t1 = (one + vdw_buf_delta)**7 / (soft_alpha * (1 - soft_lambda)**2 + (rho + vdw_buf_delta)**7)
                t2 = (one + vdw_buf_gamma) / (soft_alpha * (1 - soft_lambda)**2 + rho7 + vdw_buf_gamma)
                dt1drho = (-seven * (rho + vdw_buf_delta)**6 * t1) / (soft_alpha * (1 - soft_lambda)**2 & 
                         + (rho + vdw_buf_delta)**7)
                dt2drho = (-seven * rho6 * t2) / (soft_alpha * (1 - soft_lambda)**2 + rho7 + vdw_buf_gamma)
          endif
         end do
         f = eps*t1*(t2 - two)
         dfdr = eps*(dt1drho*(t2 - two) + t1*dt2drho)*drhodr
         if ( do_vdw_taper == 1 .and. delr2 > vdw_switch_on_2 )then
            delr3 = delr2*delr
            delr4 = delr3*delr
            delr5 = delr4*delr
            switch = c5*delr5 + c4*delr4 + c3*delr3 + c2*delr2 + &
                     c1*delr + c0
            dswitch_dr = five*c5*delr4 + four*c4*delr3 + three*c3*delr2 + &
                         two*c2*delr + c1
            dfdr = switch*dfdr + f*dswitch_dr
            f = switch*f
         endif
         ene_vdw = ene_vdw + f
         term = dfdr/delr
         dfx = term*delx
         dfy = term*dely
         dfz = term*delz
         frc(1,i) = frc(1,i) + (one-wi)*dfx
         frc(2,i) = frc(2,i) + (one-wi)*dfy
         frc(3,i) = frc(3,i) + (one-wi)*dfz
         frc(1,ih) = frc(1,ih) + wi*dfx
         frc(2,ih) = frc(2,ih) + wi*dfy
         frc(3,ih) = frc(3,ih) + wi*dfz
         frc(1,j) = frc(1,j) - (one-wj)*dfx
         frc(2,j) = frc(2,j) - (one-wj)*dfy
         frc(3,j) = frc(3,j) - (one-wj)*dfz
         frc(1,jh) = frc(1,jh) - wj*dfx
         frc(2,jh) = frc(2,jh) - wj*dfy
         frc(3,jh) = frc(3,jh) - wj*dfz
         vxx = vxx + dfx*delx
         vxy = vxy + dfx*dely
         vxz = vxz + dfx*delz
         vyy = vyy + dfy*dely
         vyz = vyz + dfy*delz
         vzz = vzz + dfz*delz
      endif
   end do
   virial(1,1) = virial(1,1) + vxx
   virial(1,2) = virial(1,2) + vxy
   virial(1,3) = virial(1,3) + vxz
   virial(2,1) = virial(2,1) + vxy
   virial(2,2) = virial(2,2) + vyy
   virial(2,3) = virial(2,3) + vyz
   virial(3,1) = virial(3,1) + vxz
   virial(3,2) = virial(3,2) + vyz
   virial(3,3) = virial(3,3) + vzz

end subroutine AM_VDW_DIRECT_ene_frc_i
!------------------------------------------------------------
end module amoeba_vdw
