#include "dprec.fh"
#include "assert.fh"

!---------------------------------------------------------
module amoeba_multipoles

  implicit none
  private

  integer,save :: do_flag
  _REAL_,parameter :: coulomb_const_kcal_per_mole = 332.05382d0 !Tinker
  integer,save :: num_chiral_frame_list
  integer,save :: start_chiral_frame_list,end_chiral_frame_list
type :: chiral_frame
  integer :: frame_index
  integer :: fourth_atom
  integer :: chirality
end type chiral_frame
  type(chiral_frame),save,allocatable :: chiral_frame_list(:)

type :: frame_def_list_entry
   integer :: frame_index ! which atomic frame does this refer to
   integer :: frame_point_number ! which frame point (1 or 2)
   integer :: vector_tail_index ! unit vector to add into point def
   integer :: vector_head_index
   integer :: num_vectors ! number of unit vec contribs to frame def point
end type frame_def_list_entry
integer, save :: num_frame_def_list
integer, save :: start_frame_def_list,end_frame_def_list
type(frame_def_list_entry),allocatable,save :: frame_def_list(:)

integer, save :: num_multipoles
integer, save :: start_multipoles,end_multipoles
_REAL_, save, allocatable :: local_multipole(:,:)
_REAL_, save, allocatable :: global_multipole(:,:)
_REAL_, save, allocatable :: torque_field(:,:)
integer,parameter :: MAXMP=35

type :: frame
   logical :: valid = .false.
   _REAL_,dimension(3) :: def_point1 = 0.d0,def_point2 = 0.d0
   _REAL_,dimension(3,3) :: axis
end type frame
type(frame),save, allocatable :: frame_list(:)
integer, save :: frame_axis_order(3) = (/3,1,2/)!gives axes with respect 
                                    !to def pts. default is z then x then y

public global_multipole,torque_field,AM_MPOLE_readparm, &
       AM_MPOLE_deallocate,AM_MPOLE_local_to_global, &
       AM_MPOLE_torque_to_force,coulomb_const_kcal_per_mole, &
       AM_MPOLE_bcast
contains
!---------------------------------------------------------
subroutine AM_MPOLE_get_start_end_lists()

  integer :: siztask,n
! first local mpoles boundaries
  call AMOEBA_get_startlist_endlist(num_multipoles,  &
                             start_multipoles,end_multipoles,siztask)
! next chiral
  start_chiral_frame_list = 0
  end_chiral_frame_list = 0
  do n = 1,num_chiral_frame_list
    if ( chiral_frame_list(n)%frame_index < start_multipoles )then
      start_chiral_frame_list = start_chiral_frame_list + 1
    endif
    if ( chiral_frame_list(n)%frame_index <= end_multipoles )then
      end_chiral_frame_list = end_chiral_frame_list + 1
      ! this will update through the last list item that works
    endif
  enddo
  ! one more to get 1st list item that works
  start_chiral_frame_list = start_chiral_frame_list + 1
  ! NOTE this list could be a problem with load balancing---since its the
  ! first few processors that will do the list--small list however
  ! note that typically all tests will come back negative--ie there is no 
  ! chirality flip---if its a backwards amino acid it will flip the first
  ! time and never again--maybe could split the list up and then broadcast 
  ! any changed local multipoles---to be fixed
  !write(6,*)'start_chiral_frame_list,end_chiral_frame_list = ', &
         !start_chiral_frame_list,end_chiral_frame_list
! next frame_def_list
  start_frame_def_list = 0
  end_frame_def_list = 0
  do n = 1,num_frame_def_list
    if ( frame_def_list(n)%frame_index < start_multipoles )then
      start_frame_def_list = start_frame_def_list + 1
    endif
    if ( frame_def_list(n)%frame_index <= end_multipoles )then
      end_frame_def_list = end_frame_def_list + 1
      ! this will update through the last list item that works
    endif
  enddo
  ! one more to get 1st list item that works
  start_frame_def_list = start_frame_def_list + 1
end subroutine AM_MPOLE_get_start_end_lists
!---------------------------------------------------------

subroutine AM_MPOLE_bcast()
   implicit none
#ifdef MPI
   integer ierr

   include 'mpif.h'
#  include "extra.h"
#  include "parallel.h"
   call mpi_bcast(do_flag,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(num_multipoles,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(num_chiral_frame_list,1,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(num_frame_def_list,1,MPI_INTEGER,0,commsander,ierr)

   if(.not.master) then
      allocate(local_multipole(10,num_multipoles),stat=ierr)
      REQUIRE(ierr==0)
      allocate(frame_list(num_multipoles),stat=ierr)
      REQUIRE(ierr==0)
      allocate(global_multipole(10,num_multipoles),stat=ierr)
      REQUIRE(ierr==0)
      allocate(torque_field(10,num_multipoles),stat=ierr)
      REQUIRE(ierr==0)
      allocate(chiral_frame_list(num_chiral_frame_list),stat=ierr)
      REQUIRE(ierr==0)
      allocate(frame_def_list(num_frame_def_list),stat=ierr)
      REQUIRE(ierr==0)
   end if
     
   call mpi_bcast(local_multipole,10*num_multipoles,MPI_DOUBLE_PRECISION,0,commsander,ierr)
   call mpi_bcast(chiral_frame_list,3*num_chiral_frame_list,MPI_INTEGER,0,commsander,ierr)
   call mpi_bcast(frame_def_list,5*num_frame_def_list,MPI_INTEGER,0,commsander,ierr)  
 
   call AM_MPOLE_get_start_end_lists()
#endif 
end subroutine AM_MPOLE_bcast

function AM_MPOLE_readparm(nf,num_atoms)
  integer :: AM_MPOLE_readparm
  integer, intent(in) :: nf,num_atoms
#include "do_flag.h"  
  integer :: n,ier,dim1
  integer, allocatable :: buf(:,:)

  AM_MPOLE_readparm = 0
  call AMOEBA_get_numlist('AMOEBA_LOCAL_FRAME_MULTIPOLES_',nf,num_multipoles)
  if ( num_multipoles <= 0 )then
    do_flag = ibclr(do_flag,VALID_BIT)
    return
  endif
  if ( num_multipoles /= num_atoms )then
    write(6,*)'mismatch between num_multipoles and num_atoms'
    call mexit(6,1)
  endif
  !allocate
  allocate(local_multipole(10,num_multipoles),stat=ier)
  REQUIRE(ier==0)
  dim1 = 10 ! symmetrized cartesian multipoles up to quadratic
  call AMOEBA_read_real_list_data('AMOEBA_LOCAL_FRAME_MULTIPOLES_',nf, &
                  dim1,num_multipoles,local_multipole)
  ! rescale them to units of angstroms---plus quadrupole correction for
  ! traceless cartesian
  call AM_MPOLE_rescale_multipoles()
  ! allocate frame list and global multipoles
  allocate(frame_list(num_multipoles),stat=ier)
  REQUIRE(ier==0)
  allocate(global_multipole(10,num_multipoles),stat=ier)
  REQUIRE(ier==0)
  allocate(torque_field(10,num_multipoles),stat=ier)
  REQUIRE(ier==0)

  call AMOEBA_get_numlist('AMOEBA_CHIRAL_FRAME_',nf,num_chiral_frame_list)
  !allocate
  if ( num_chiral_frame_list > 0 )then
    allocate(buf(3,num_chiral_frame_list),stat=ier)
    REQUIRE(ier==0)
    allocate(chiral_frame_list(num_chiral_frame_list),stat=ier)
    REQUIRE(ier==0)
    dim1 = 3
    call AMOEBA_read_list_data('AMOEBA_CHIRAL_FRAME_',nf, &
                  dim1,num_chiral_frame_list,buf)
    do n = 1,num_chiral_frame_list
      chiral_frame_list(n)%frame_index = buf(1,n)
      chiral_frame_list(n)%fourth_atom = buf(2,n)
      chiral_frame_list(n)%chirality = buf(3,n)
    enddo
    deallocate(buf)
  endif
  call AMOEBA_get_numlist('AMOEBA_FRAME_DEF_',nf,num_frame_def_list)
  if ( num_frame_def_list <= 0 )then
    do_flag = ibclr(do_flag,VALID_BIT)
    return
  endif
  !allocate
  allocate(buf(5,num_frame_def_list),stat=ier)
  REQUIRE(ier==0)
  allocate(frame_def_list(num_frame_def_list),stat=ier)
  REQUIRE(ier==0)
  dim1 = 5
  call AMOEBA_read_list_data('AMOEBA_FRAME_DEF_',nf, &
                  dim1,num_frame_def_list,buf)
  do n = 1,num_frame_def_list
    frame_def_list(n)%frame_index = buf(1,n)
    frame_def_list(n)%frame_point_number = buf(2,n)
    frame_def_list(n)%vector_tail_index = buf(3,n)
    frame_def_list(n)%vector_head_index = buf(4,n)
    frame_def_list(n)%num_vectors = buf(5,n)
  enddo
  deallocate(buf)

  call AM_MPOLE_get_start_end_lists()
  AM_MPOLE_readparm = 1
  do_flag = ibset(do_flag,VALID_BIT)
end function AM_MPOLE_readparm
!---------------------------------------------------------
subroutine AM_MPOLE_rescale_multipoles()

  integer j,n
  _REAL_ bohr,traced
! now change units from Bohr to Angstroms
   bohr = 0.5291772083d0
   do n = 1,num_multipoles
     do j = 2,4
       local_multipole(j,n) = bohr*local_multipole(j,n)
     enddo
     traced = 2.d0 / 3.d0  !spherical harm expansion diff from Taylor's
     ! see stone's book
     do j = 5,10
       local_multipole(j,n) = traced * bohr*bohr*local_multipole(j,n)
     enddo
  enddo
end subroutine AM_MPOLE_rescale_multipoles

!---------------------------------------------------------
subroutine AM_MPOLE_deallocate()

  if ( allocated(chiral_frame_list) )deallocate(chiral_frame_list)
  if ( allocated(frame_def_list) )deallocate(frame_def_list)
  if ( allocated(local_multipole) )deallocate(local_multipole)
end subroutine AM_MPOLE_deallocate
!---------------------------------------------------------
subroutine AM_MPOLE_local_to_global(crd)
  _REAL_,intent(in) :: crd(3,*)
#include "do_flag.h"
  
  integer :: n

  if ( .not. btest(do_flag,VALID_BIT) )return
  ! clear frames
  do n = start_multipoles,end_multipoles
    frame_list(n)%def_point1 = 0.d0
    frame_list(n)%def_point2 = 0.d0
    frame_list(n)%axis = 0.d0
  enddo
  call AM_MPOLE_build_frame_def_pts(crd, &
               start_frame_def_list,end_frame_def_list, &
               frame_def_list,frame_list)
  call AM_MPOLE_def_pts_to_axes(start_multipoles,end_multipoles, &
                                frame_axis_order,frame_list)
  call AM_MPOLE_check_chirality(crd, &
           start_chiral_frame_list,end_chiral_frame_list,frame_axis_order, &
           chiral_frame_list,frame_list,local_multipole)
  call AM_MPOLE_rotate_multipole(start_multipoles,end_multipoles, &
                   frame_list,local_multipole,global_multipole)
end subroutine AM_MPOLE_local_to_global
!---------------------------------------------------------
subroutine AM_MPOLE_torque_to_force(numatoms,crd,frc,virial)
  use stack
  character(kind=1,len=24) :: routine = "AM_MPOLE_torque_to_force"
  integer,intent(in) :: numatoms
  _REAL_,intent(in) :: crd(3,*)
  _REAL_,intent(inout) :: frc(3,*),virial(3,3)

  integer :: p_de_drotsite,p_de_drotframe,p_de_ddefpt

  call get_stack(p_de_drotsite,3*numatoms,routine)
  call get_stack(p_de_drotframe,3*numatoms,routine)
  call get_stack(p_de_ddefpt,3*2*numatoms,routine)
  if(.not. rstack_ok)then
      deallocate(r_stack)
      allocate(r_stack(1:lastrst),stat=alloc_ier)
      call reassign_rstack(routine)
  endif
  REQUIRE(rstack_ok)
  call AM_MPOLE_get_de_drot_mpoles(numatoms,global_multipole,torque_field, &
                                   r_stack(p_de_drotsite))
  call AM_MPOLE_accum_de_dframe_rot(numatoms,frame_axis_order, &
                                    r_stack(p_de_drotsite), &
                                    frame_list,r_stack(p_de_drotframe))
  call AM_MPOLE_accum_de_ddefpts(numatoms,frame_axis_order,   &  
                                 r_stack(p_de_drotframe),frame_list, &
                                 r_stack(p_de_ddefpt))
  call AM_MPOLE_de_ddefpts_to_atoms(num_frame_def_list,frame_def_list, &
                                    r_stack(p_de_ddefpt),crd,frc,virial)
  call free_stack(p_de_ddefpt,routine)
  call free_stack(p_de_drotframe,routine)
  call free_stack(p_de_drotsite,routine)
end subroutine AM_MPOLE_torque_to_force
!---------------------------------------------------------
! routine needed by AM_MPOLE_local_to_global,AM_MPOLE_torque_to_force
!---------------------------------------------------------
subroutine AM_MPOLE_build_frame_def_pts(crd,startlist,endlist, &
                                        fr_deflist,fr_list)
  _REAL_,intent(in) :: crd(3,*)
  integer,intent(in) :: startlist,endlist
  type(frame_def_list_entry),intent(in) :: fr_deflist(*)
  type(frame),intent(inout) :: fr_list(*)

  integer n,i1,i2,k,m,p
  _REAL_ dx,dy,dz,wt

  do n = startlist,endlist
    m = fr_deflist(n)%frame_index
    p = fr_deflist(n)%frame_point_number
    i1 = fr_deflist(n)%vector_tail_index
    i2 = fr_deflist(n)%vector_head_index
    k = fr_deflist(n)%num_vectors
    if ( i1 > 0 .and. i2 > 0 )then
      fr_list(m)%valid = .true.
      dx = crd(1,i2) - crd(1,i1)
      dy = crd(2,i2) - crd(2,i1)
      dz = crd(3,i2) - crd(3,i1)
      wt = k*sqrt(dx*dx+dy*dy+dz*dz) ! divide by length of i1i2 times num such
      if ( p == 1 )then
        fr_list(m)%def_point1(1) = fr_list(m)%def_point1(1) + dx / wt
        fr_list(m)%def_point1(2) = fr_list(m)%def_point1(2) + dy / wt
        fr_list(m)%def_point1(3) = fr_list(m)%def_point1(3) + dz / wt
      elseif ( p == 2 )then
        fr_list(m)%def_point2(1) = fr_list(m)%def_point2(1) + dx / wt
        fr_list(m)%def_point2(2) = fr_list(m)%def_point2(2) + dy / wt
        fr_list(m)%def_point2(3) = fr_list(m)%def_point2(3) + dz / wt
      else
        write(6,*)'AM_MPOLE_build_frame_def_pts: serious problem!'
        call mexit(6,1)
      endif
    else
      fr_list(m)%valid = .false.
    endif
  enddo
end subroutine AM_MPOLE_build_frame_def_pts
!---------------------------------------------------------
subroutine AM_MPOLE_def_pts_to_axes(startlist,endlist,axis_order,fr_list)
  integer,intent(in) :: startlist,endlist,axis_order(3)
  type(frame),intent(inout) :: fr_list(*)

  integer :: i,n,k1,k2,k3
  _REAL_ :: u(3),v(3),w(3),siz,dot

  k1 = axis_order(1)
  k2 = axis_order(2)
  k3 = axis_order(3)
  do n = startlist,endlist
    if ( fr_list(n)%valid )then
! u is unit vector in direction of primary def pt
      u(1) = fr_list(n)%def_point1(1)
      u(2) = fr_list(n)%def_point1(2)
      u(3) = fr_list(n)%def_point1(3)
      siz = sqrt(u(1)*u(1)+u(2)*u(2)+u(3)*u(3))
      u(1) = u(1) / siz
      u(2) = u(2) / siz
      u(3) = u(3) / siz
! v is unit vector given by component of secondary pt orthog to u
      v(1) = fr_list(n)%def_point2(1)
      v(2) = fr_list(n)%def_point2(2)
      v(3) = fr_list(n)%def_point2(3)
      dot = u(1)*v(1)+u(2)*v(2)+u(3)*v(3)
      v(1) = v(1) - dot*u(1)
      v(2) = v(2) - dot*u(2)
      v(3) = v(3) - dot*u(3)
      siz = sqrt(v(1)*v(1)+v(2)*v(2)+v(3)*v(3))
      v(1) = v(1) / siz
      v(2) = v(2) / siz
      v(3) = v(3) / siz
! w is u cross v
      w(1) = u(2)*v(3) - u(3)*v(2)
      w(2) = u(3)*v(1) - u(1)*v(3)
      w(3) = u(1)*v(2) - u(2)*v(1)
! build axes
      do i = 1,3
        fr_list(n)%axis(i,k1) = u(i)
        fr_list(n)%axis(i,k2) = v(i)
        fr_list(n)%axis(i,k3) = w(i)
      enddo
    endif ! fr_list(n)%valid
  enddo
end subroutine AM_MPOLE_def_pts_to_axes
!---------------------------------------------------------
subroutine AM_MPOLE_check_chirality(crd,startlist,endlist,axis_order, &
                                    chiral_frlist,fr_list,loc_mpole)
  _REAL_,intent(in) :: crd(3,*)
  integer,intent(in) :: startlist,endlist,axis_order(3)
  type(chiral_frame),intent(inout) :: chiral_frlist(*)
  type(frame),intent(in) :: fr_list(*)
  _REAL_,intent(inout) :: loc_mpole(10,*)

  integer :: n,j,k,k3
  _REAL_ dx,dy,dz,dot

  k3 = axis_order(3)
  do n = startlist,endlist
    j = chiral_frlist(n)%frame_index
    if ( .not. fr_list(j)%valid )then
      write(6,*)'chiral frame: serious problem in frame ',j
      call mexit(6,1)
    endif
    k = chiral_frlist(n)%fourth_atom
    dx = crd(1,k) - crd(1,j)
    dy = crd(2,k) - crd(2,j)
    dz = crd(3,k) - crd(3,j)
    dot = dx*fr_list(j)%axis(1,k3)  + &
          dy*fr_list(j)%axis(2,k3)  + &
          dz*fr_list(j)%axis(3,k3)
    ! this should be negative if chirality is positive and vice versa
    if ( ((chiral_frlist(n)%chirality == 1) .and. (dot > 0)) .or. &
         ((chiral_frlist(n)%chirality == -1) .and. (dot < 0)) )then
       chiral_frlist(n)%chirality = -chiral_frlist(n)%chirality
       if ( k3 == 1 )then
         loc_mpole(2,n) = -loc_mpole(2,n)
         loc_mpole(8,n) = -loc_mpole(8,n)
         loc_mpole(9,n) = -loc_mpole(9,n)
       elseif ( k3 == 2 )then
         loc_mpole(3,n) = -loc_mpole(3,n)
         loc_mpole(8,n) = -loc_mpole(8,n)
         loc_mpole(10,n) = -loc_mpole(10,n)
       else
         loc_mpole(4,n) = -loc_mpole(4,n)
         loc_mpole(9,n) = -loc_mpole(9,n)
         loc_mpole(10,n) = -loc_mpole(10,n)
       endif
    endif
  enddo
end subroutine AM_MPOLE_check_chirality
!---------------------------------------------------------
subroutine AM_MPOLE_rotate_multipole(startlist,endlist, &
                   fr_list,loc_mpole,glob_mpole)
  integer,intent(in) :: startlist,endlist
  type(frame),intent(inout) :: fr_list(*)
  _REAL_,intent(in) :: loc_mpole(10,*)
  _REAL_,intent(out) :: glob_mpole(10,*)

  integer order,dimxy,n,j
  _REAL_ Mpole_xy(MAXMP*MAXMP)
  order = 10
  dimxy = 10
  do n = startlist,endlist
    if ( fr_list(n)%valid  )then
      call XFORM_MPOLE_matrix(fr_list(n)%axis,Mpole_xy,order)
      call XFORM_MPOLE(Mpole_xy,dimxy,loc_mpole(:,n),glob_mpole(:,n),order)
    else
      do j = 1,10
        glob_mpole(j,n) = loc_mpole(j,n) !for charge only case e.g. ions
      enddo
    endif
  enddo

end subroutine AM_MPOLE_rotate_multipole
!---------------------------------------------------------
subroutine AM_MPOLE_get_de_drot_mpoles(numatoms,global_multipole, &
                                   torque_field,de_drotsite)
  integer ,intent(in) :: numatoms
  _REAL_,intent(in) :: global_multipole(10,*),torque_field(10,*)
  _REAL_,intent(out) :: de_drotsite(3,*)

!  get the derivative of electrostatic energy with respect to infinitesmal 
!  rotations of atomic frames about the x,y,z global axes
! Basic idea--electrostatic energy given by dot product of cartesian
! multipoles and the electrostatic potential and its cartesian derivatives
! i.e. ene = 1/2 *(q*phi + mux*dphidx + ...)
! Thus derivative obtained by rotating multipoles infinitesmally

  integer :: i,j,k,n,dimxy,order
  _REAL_ :: DMP_x(MAXMP*MAXMP),DMP_y(MAXMP*MAXMP),  &
            DMP_z(MAXMP*MAXMP),A_xy(3,3),DA_xy(3,3),  &
            Tmp_x(10),Tmp_y(10),Tmp_z(10)

  order = 10
! to get de_drot we calculate the deriv of mpole wrt infinitesmal
! rotations about x,y and z axis
  do i = 1,3
    do j = 1,3
      A_xy(i,j) = 0.d0
    enddo
    A_xy(i,i) = 1.d0
  enddo
     
! do the maximal order
! x-axis rotation of dtheta
  do i = 1,3
    do j = 1,3
      DA_xy(i,j) = 0.d0
    enddo
  enddo
  DA_xy(3,2) = 1.d0
  DA_xy(2,3) = -1.d0
  call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_x,order)
! y-axis
  do i = 1,3
    do j = 1,3
      DA_xy(i,j) = 0.d0
    enddo
  enddo
  DA_xy(3,1) = -1.d0
  DA_xy(1,3) = 1.d0
  call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_y,order)

! z-axis
  do i = 1,3
    do j = 1,3
      DA_xy(i,j) = 0.d0
    enddo
  enddo
  DA_xy(2,1) = 1.d0
  DA_xy(1,2) = -1.d0
  call XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMP_z,order)
 
  dimxy = order
  do n = 1,numatoms
    de_drotsite(1,n) = 0.d0
    de_drotsite(2,n) = 0.d0
    de_drotsite(3,n) = 0.d0
    call XFORM_MPOLE(DMP_x,dimxy,global_multipole(:,n),Tmp_x,order)
    call XFORM_MPOLE(DMP_y,dimxy,global_multipole(:,n),Tmp_y,order)
    call XFORM_MPOLE(DMP_z,dimxy,global_multipole(:,n),Tmp_z,order)
    do k = 1,order
      de_drotsite(1,n) = de_drotsite(1,n) +  &
                         coulomb_const_kcal_per_mole*Tmp_x(k)*torque_field(k,n)
      de_drotsite(2,n) = de_drotsite(2,n) +  &
                         coulomb_const_kcal_per_mole*Tmp_y(k)*torque_field(k,n)
      de_drotsite(3,n) = de_drotsite(3,n) +  &
                         coulomb_const_kcal_per_mole*Tmp_z(k)*torque_field(k,n)
    enddo
  enddo
end subroutine AM_MPOLE_get_de_drot_mpoles
!---------------------------------------------------------
subroutine AM_MPOLE_accum_de_dframe_rot(numatoms,axis_order,de_drotsite, &
                                        fr_list,de_drotframe)
  integer,intent(in) :: numatoms,axis_order(3)
  _REAL_,intent(in) :: de_drotsite(3,numatoms)
  type(frame),intent(in) :: fr_list(*)
  _REAL_,intent(out) :: de_drotframe(3,numatoms)

  integer :: n,k1,k2,k3
  _REAL_ :: p2unit(3),siz

  k1 = axis_order(1)
  k2 = axis_order(2)
  k3 = axis_order(3)
! deriv of energy with respect to rotation about unit vectors along
! p1 p2 and their cross product
! note that unit vector along p1 corresponds to k1st frame axis
! and unit vector in p1 x p2 direction corresponds to k3rd frame axis
! the energy derivative with respect to rotation about any unit vector
! is for each mpole given by the dot product of de_drotpole 
! (which is negative of torque due to that mpole) with the unit vector
        
  do n = 1,numatoms
    !initialize
    de_drotframe(1,n) = 0.d0
    de_drotframe(2,n) = 0.d0
    de_drotframe(3,n) = 0.d0
    if ( fr_list(n)%valid )then
      siz = sqrt(fr_list(n)%def_point2(1)*fr_list(n)%def_point2(1) + &
                 fr_list(n)%def_point2(2)*fr_list(n)%def_point2(2) + &
                 fr_list(n)%def_point2(3)*fr_list(n)%def_point2(3))
      p2unit(1) = fr_list(n)%def_point2(1) / siz
      p2unit(2) = fr_list(n)%def_point2(2) / siz
      p2unit(3) = fr_list(n)%def_point2(3) / siz
      de_drotframe(k1,n) = de_drotframe(k1,n) +   &
                           de_drotsite(1,n)*fr_list(n)%axis(1,k1) +  &
                           de_drotsite(2,n)*fr_list(n)%axis(2,k1) +  &
                           de_drotsite(3,n)*fr_list(n)%axis(3,k1)
      de_drotframe(k2,n) = de_drotframe(k2,n) +   &
                           de_drotsite(1,n)*p2unit(1) +   &
                           de_drotsite(2,n)*p2unit(2) +   &
                           de_drotsite(3,n)*p2unit(3)
      de_drotframe(k3,n) = de_drotframe(k3,n) +    &
                           de_drotsite(1,n)*fr_list(n)%axis(1,k3) +  &
                           de_drotsite(2,n)*fr_list(n)%axis(2,k3) +  &
                           de_drotsite(3,n)*fr_list(n)%axis(3,k3)
    endif ! fr_list(n)%valid
  enddo
      
end subroutine AM_MPOLE_accum_de_dframe_rot
!---------------------------------------------------------
subroutine AM_MPOLE_accum_de_ddefpts(numatoms,axis_order,   &  
                                     de_drotframe,fr_list,de_ddefpt)
  integer,intent(in) :: numatoms,axis_order(3)
  _REAL_,intent(in) :: de_drotframe(3,numatoms)
  type(frame),intent(in) :: fr_list(*)
  _REAL_,intent(out) :: de_ddefpt(3,2,numatoms)
! get the derivs of energy with respect to movement of defpoints
! expressed in the local frame coord system
  integer :: n,j,k1,k2,k3
  _REAL_ :: p1(3),p2(3),p2unit(3),p2perp1(3),p1perp2(3),  &
            u(3),v(3),w(3), &
            sizp1perp2,sizp2perp1,dot12,dot21,  &
            sizp1,sizp2,dedv,dedw, &
            de_drotu,de_drotv,de_drotw

  k1 = axis_order(1)
  k2 = axis_order(2)
  k3 = axis_order(3)
  do n = 1,numatoms
    if ( fr_list(n)%valid )then
      do j = 1,3
        p1(j) = fr_list(n)%def_point1(j)
        p2(j) = fr_list(n)%def_point2(j)
        u(j) = fr_list(n)%axis(j,k1)
        v(j) = fr_list(n)%axis(j,k2)
        w(j) = fr_list(n)%axis(j,k3)
      enddo
      de_drotu = de_drotframe(k1,n)
      de_drotv = de_drotframe(k2,n)
      de_drotw = de_drotframe(k3,n)

      sizp1 = sqrt( p1(1)**2 + p1(2)**2 + p1(3)**2 )
      sizp2 = sqrt( p2(1)**2 + p2(2)**2 + p2(3)**2 )
      do j = 1,3
        p2unit(j) = p2(j) / sizp2
      enddo
      dot21 = u(1)*p2(1) + u(2)*p2(2) + u(3)*p2(3)
      dot12 = p1(1)*p2unit(1) + p1(2)*p2unit(2) + p1(3)*p2unit(3)
      do j = 1,3
        p2perp1(j) = p2(j) - dot21*u(j)
        p1perp2(j) = p1(j) - dot12*p2unit(j)
      enddo
      sizp2perp1 = sqrt(p2perp1(1)**2 + p2perp1(2)**2 + p2perp1(3)**2)
      sizp1perp2 = sqrt(p1perp2(1)**2 + p1perp2(2)**2 + p1perp2(3)**2)
! def point one is along axis one. movement du parallel to that axis does
! not rotate the frame..so deriv is zero
!       dedu = 0.d0
! movement dv in v-axis direction corresponds to rotation about local w-axis
! of dtheta = dv/sizp1; thus a change in energy of 
!    dE = dedrotw*dtheta = dedrotw*dv/sizp1
!    dE/dv = dedrotw /sizp1
      dedv = de_drotw/sizp1
! movement dw in w-axis direction corresponds to rotation about p2unit
! of dtheta = -dw/sizp1perp2 (clockwise rotation) 
      dedw = -de_drotv/sizp1perp2
! Now convert to derivs wrt x,y,z. u = p.u = x*u(1)+y*u(2)+z*u(3)
! thus dudx = u(1). similaryl dvdx = v(1)
      de_ddefpt(1,1,n) = dedv*v(1) + dedw*w(1)
      de_ddefpt(2,1,n) = dedv*v(2) + dedw*w(2)
      de_ddefpt(3,1,n) = dedv*v(3) + dedw*w(3)
 
! for point 2..any movement in the local uv plane does not affect the frame
!       dedu = 0.d0
!       dedv = 0.d0
! movement dw in w direction corresponds to rotation about local u-axis 
! of dtheta = dw/sizp2perpu
      dedw = de_drotu/sizp2perp1
      de_ddefpt(1,2,n) = dedw*w(1)
      de_ddefpt(2,2,n) = dedw*w(2)
      de_ddefpt(3,2,n) = dedw*w(3)
    else ! fr_list(n)%valid .eq. .false.
      do j = 1,3
        de_ddefpt(j,1,n) = 0.d0
        de_ddefpt(j,2,n) = 0.d0
      enddo
    endif !fr_list(n)%valid
  enddo !n = 1,numatoms
end subroutine AM_MPOLE_accum_de_ddefpts
!---------------------------------------------------------
subroutine AM_MPOLE_de_ddefpts_to_atoms(num_fr_deflist,fr_deflist, &
                                        de_ddefpt,crd,frc,virial)
  use  constants, only : zero,half
  integer,intent(in) :: num_fr_deflist
  type(frame_def_list_entry),intent(in) :: fr_deflist(*)
  _REAL_,intent(in) :: de_ddefpt(3,2,*),crd(3,*)
  _REAL_,intent(inout) :: frc(3,*),virial(3,3)

  integer :: n,i,j,k,l,m
  _REAL_ :: siz,siz2,dx,dy,dz,dux_dx,dux_dy,dux_dz, &
            duy_dy,duy_dz,duz_dz,dedux,deduy,deduz,dedx,dedy,dedz,  &
            siz3inv,sizinv,vxx,vxy,vxz,vyx,vyy,vyz,vzx,vzy,vzz

  vxx = zero
  vxy = zero
  vxz = zero
  vyx = zero
  vyy = zero
  vyz = zero
  vzx = zero
  vzy = zero
  vzz = zero
  do n = 1,num_fr_deflist
    i = fr_deflist(n)%vector_tail_index
    j = fr_deflist(n)%vector_head_index
    if ( i > 0 .and. j > 0 )then
      k = fr_deflist(n)%frame_index
      l = fr_deflist(n)%frame_point_number
      m = fr_deflist(n)%num_vectors
      dx = crd(1,j) - crd(1,i)
      dy = crd(2,j) - crd(2,i)
      dz = crd(3,j) - crd(3,i)
      siz2 = dx*dx+dy*dy+dz*dz
      siz = sqrt(siz2)
      siz3inv = 1.d0/(siz2*siz)
      sizinv = 1.d0/siz
! ux, uy, uz are given by dx/siz, dy/siz, and dz/siz 
      dux_dx = sizinv - dx*dx*siz3inv
      dux_dy = -dx*dy*siz3inv   ! note duy_dx = dux_dy use this below
      dux_dz = -dx*dz*siz3inv
      duy_dy = sizinv - dy*dy*siz3inv
      duy_dz = -dy*dz*siz3inv
      duz_dz = sizinv - dz*dz*siz3inv
! the derivs of E wrt coordinates of unit vector in ij direction are given
! by (1/m) times the derivs of E wrt coords of def point (l,k)
! since the def point is the simple average of m of these unit vectors
      dedux = de_ddefpt(1,l,k) / m
      deduy = de_ddefpt(2,l,k) / m
      deduz = de_ddefpt(3,l,k) / m
! now apply chain rule, using symmetry e.g. dux_dy = duy_dx
      dedx = dedux*dux_dx + deduy*dux_dy + deduz*dux_dz
      dedy = dedux*dux_dy + deduy*duy_dy + deduz*duy_dz
      dedz = dedux*dux_dz + deduy*duy_dz + deduz*duz_dz
! finally apply forces. note force is negative of deriv of energy wrt position
! also note e.g. deriv of dx wrt x position of atoms i and j is -1,+1
      frc(1,i) = frc(1,i) + dedx
      frc(2,i) = frc(2,i) + dedy
      frc(3,i) = frc(3,i) + dedz
      frc(1,j) = frc(1,j) - dedx
      frc(2,j) = frc(2,j) - dedy
      frc(3,j) = frc(3,j) - dedz
      vxx = vxx + dedx*dx
      vxy = vxy + dedx*dy
      vxz = vxz + dedx*dz
      vyx = vyx + dedy*dx
      vyy = vyy + dedy*dy
      vyz = vyz + dedy*dz
      vzx = vzx + dedz*dx
      vzy = vzy + dedz*dy
      vzz = vzz + dedz*dz
    endif 
  enddo 
  virial(1,1) = virial(1,1) + vxx
  virial(1,2) = virial(1,2) + half*(vxy + vyx)
  virial(1,3) = virial(1,3) + half*(vxz + vzx)
  virial(2,1) = virial(2,1) + half*(vxy + vyx)
  virial(2,2) = virial(2,2) + vyy
  virial(2,3) = virial(2,3) + half*(vyz + vzy)
  virial(3,1) = virial(3,1) + half*(vxz + vzx)
  virial(3,2) = virial(3,2) + half*(vyz + vzy)
  virial(3,3) = virial(3,3) + vzz
end subroutine AM_MPOLE_de_ddefpts_to_atoms
!---------------------------------------------------------
end module amoeba_multipoles
!---------------------------------------------------------
! utility routines
!---------------------------------------------------------
subroutine XFORM_MPOLE_matrix(A_xy,Mpole_xy,order)
! A_xy is matrix of coefficients expanding y in terms of x
! i.e. y_i = A_i1 * x_1 + A_i2 * x_2 + A_i3 * x_3
! Mpole_xy is resulting matrix getting expansion of multipoles
! with basis y in terms of those with basis x
! order is order of multipoles 1 for charge, 4 for charge-dipoles
! 10 for charge,dipole,quadrupole, 20 for up to octupole and
! finally 35 for up to hexadecapole
! -------------------------------------------------------
! terms are obtained by examining taylor expansions of energies in terms
! of multipoles times field derivatives in coord system x or y,
! employing the chain rule to transform field derivatives as in
! subroutine XFORM_MPOLE_field_matrix below, reversing summation order and
! collecting terms. Note that symmetry applies, so that for example
! mpole_x(8) is the sum of the x1x2 and x2x1 multipole components
! this symmetry is applied in the matrix
!----------------------------------------------------------------
 
! Mpole 1 is charge... Mpole 10 is x_3,x_3 quadrupolar coefficients etc.
  implicit none
  integer order
  _REAL_ A_xy(3,3), Mpole_xy(order,order)
  integer i,j,k,l,m,n,p,ind1,ind2

  integer, parameter :: &
        qind1(6)  = (/1,2,3,1,1,2/), &
        qind2(6)  = (/1,2,3,2,3,3/), &
        oind1(10) = (/1,2,3,1,1,2,2,3,3,1/), &
        oind2(10) = (/1,2,3,1,1,2,2,3,3,2/), &
        oind3(10) = (/1,2,3,2,3,1,3,1,2,3/), &
        hind1(15) = (/1,2,3,1,1,2,2,3,3,1,1,2,1,2,3/), &
        hind2(15) = (/1,2,3,1,1,2,2,3,3,1,1,2,1,2,3/), &
        hind3(15) = (/1,2,3,1,1,2,2,3,3,2,3,3,2,1,1/), &
        hind4(15) = (/1,2,3,2,3,1,3,1,2,2,3,3,3,3,2/)

! CHARGE case
  Mpole_xy(1,1) = 1.d0
  if ( order .eq. 1 )return
  do j = 2,order
   Mpole_xy(1,j) = 0.d0
   Mpole_xy(j,1) = 0.d0
  enddo
! DIPOLES
  do j = 2,4
    do k = 2,4
      Mpole_xy(j,k) = A_xy(j-1,k-1)
    enddo
  enddo
  if ( order .eq. 4 )return
  do j = 2,4
    do k = 5,order
      Mpole_xy(j,k) = 0.d0
      Mpole_xy(k,j) = 0.d0
    enddo
  enddo
! QUADRUPOLES
  do ind1 = 1,3
    k = qind1(ind1)
    do ind2 = 1,6
      i = qind1(ind2)
      j = qind2(ind2)
      Mpole_xy(ind1+4,ind2+4) = A_xy(k,i)*A_xy(k,j)
    enddo
  enddo
! Q'_kl
  do ind1 = 4,6
    k = qind1(ind1)
    l = qind2(ind1)
    do ind2 = 1,6
      i = qind1(ind2)
      j = qind2(ind2)
      Mpole_xy(ind1+4,ind2+4) = A_xy(k,i)*A_xy(l,j) + A_xy(k,j)*A_xy(l,i)
    enddo
  enddo
  if ( order .eq. 10 )return
  do j = 5,10
    do k = 11,order
      Mpole_xy(k,j) = 0.d0
      Mpole_xy(j,k) = 0.d0
    enddo
  enddo
! OCTUPOLES
! O'_lll
  do ind1 = 1,3
    l = oind1(ind1)
    do ind2 = 1,10
      i = oind1(ind2)
      j = oind2(ind2)
      k = oind3(ind2)
      Mpole_xy(ind1+10,ind2+10) = A_xy(l,i)*A_xy(l,j)*A_xy(l,k)
    enddo
  enddo
! O'_llm
  do ind1 = 4,9
    l = oind1(ind1)
    m = oind3(ind1)
    do ind2 = 1,9
      i = oind1(ind2)
      k = oind3(ind2)
      Mpole_xy(ind1+10,ind2+10) =     &
           A_xy(l,i)*A_xy(l,i)*A_xy(m,k) +   &
           2.d0 * A_xy(l,i)*A_xy(m,i)*A_xy(l,k)
    enddo
    Mpole_xy(ind1+10,20) =   &
           A_xy(l,1)*A_xy(l,2)*A_xy(m,3) +   &
           A_xy(l,1)*A_xy(m,2)*A_xy(l,3) +   &
           A_xy(m,1)*A_xy(l,2)*A_xy(l,3)
  enddo
! O'_123
  Mpole_xy(20,11) = 6.d0*A_xy(1,1)*A_xy(2,1)*A_xy(3,1)
  Mpole_xy(20,12) = 6.d0*A_xy(1,2)*A_xy(2,2)*A_xy(3,2)
  Mpole_xy(20,13) = 6.d0*A_xy(1,3)*A_xy(2,3)*A_xy(3,3)
  do ind2 = 4,9
    i = oind1(ind2)
    k = oind3(ind2)
    Mpole_xy(20,10+ind2) = 2.d0*   &
           (  A_xy(1,i)*A_xy(2,i)*A_xy(3,k) +   &
              A_xy(1,i)*A_xy(2,k)*A_xy(3,i) +   &
              A_xy(1,k)*A_xy(2,i)*A_xy(3,i) )
  enddo
  Mpole_xy(20,20) =    &
              A_xy(1,1)*A_xy(2,2)*A_xy(3,3) +   &
              A_xy(1,1)*A_xy(3,2)*A_xy(2,3) +   &
              A_xy(2,1)*A_xy(1,2)*A_xy(3,3) +   &
              A_xy(2,1)*A_xy(3,2)*A_xy(1,3) +   &
              A_xy(3,1)*A_xy(1,2)*A_xy(2,3) +   &
              A_xy(3,1)*A_xy(2,2)*A_xy(1,3) 
  if ( order .eq. 20 )return
  do j = 11,20
    do k = 21,order
      Mpole_xy(k,j) = 0.d0
      Mpole_xy(j,k) = 0.d0
    enddo
  enddo
! HEXADECAPOLES
! H'_mmmm
  do ind1 = 1,3
    m = hind1(ind1)
    do ind2 = 1,15
      i = hind1(ind2)
      j = hind2(ind2)
      k = hind3(ind2)
      l = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) = A_xy(m,i)*A_xy(m,j)*A_xy(m,k)*A_xy(m,l)
    enddo
  enddo
! H'_mmmn
  do ind1 = 4,9
    m = hind1(ind1)
    n = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) =   &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,i)*A_xy(n,j) + 3.d0*  &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(m,j) 
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) =    &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,j)*A_xy(n,k) +        &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,k)*A_xy(n,j) + 2.d0*  &
               A_xy(m,i)*A_xy(m,j)*A_xy(m,k)*A_xy(n,i) 
    enddo
  enddo
! H'_mmnn
  do ind1 = 10,12
    m = hind1(ind1)
    n = hind3(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) = 3.d0 * (   &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(n,j) +    &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(n,i)  )
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) =   &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(n,k) +         &
               A_xy(m,j)*A_xy(m,k)*A_xy(n,i)*A_xy(n,i) +  2.d0*  &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(n,k) +  2.d0*  &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,i)*A_xy(n,j) 
    enddo
  enddo
! H'_mmnp
  do ind1 = 13,15
    m = hind1(ind1)
    n = hind3(ind1)
    p = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) =     &
              3.d0*A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(p,j) +     &
              3.d0*A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(p,i) +     &
              6.d0*A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(p,i) 
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Mpole_xy(ind1+20,ind2+20) =     &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(p,k) +     &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,k)*A_xy(p,j) + 2.d0 * (    &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(p,k) +     &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,k)*A_xy(p,i) +     &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,i)*A_xy(p,j) +     &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,j)*A_xy(p,i) +     &
               A_xy(m,j)*A_xy(m,k)*A_xy(n,i)*A_xy(p,i) )
    enddo
  enddo
  return
end subroutine XFORM_MPOLE_matrix
!------------------------------------------------------------
subroutine XFORM_MPOLE_deriv_matrix(A_xy,DA_xy,DMpole_xy,order)
! calculate derivative of Mpole_xy as DMpole_xy
! A_xy is matrix of coefficients expanding y in terms of x
! i.e. y_i = A_i1 * x_1 + A_i2 * x_2 + A_i3 * x_3
! DA_xy is deriv of A_xy wrt to some parameter
! DMpole_xy is deriv of Mpole_xy wrt to the same parameter
! order is order of multipoles 1 for charge, 4 for charge-dipoles
! 10 for charge,dipole,quadrupole, 20 for up to octupole and
! finally 35 for up to hexadecapole
 
! Mpole 1 is charge... Mpole 10 is x_3,x_3 quadrupolar coefficients etc.
  implicit none
  integer order
  _REAL_ A_xy(3,3)
  _REAL_ DA_xy(3,3), DMpole_xy(order,order)
  integer i,j,k,l,m,n,p,ind1,ind2

  integer, parameter :: &
        qind1(6)  = (/1,2,3,1,1,2/), &
        qind2(6)  = (/1,2,3,2,3,3/), &
        oind1(10) = (/1,2,3,1,1,2,2,3,3,1/), &
        oind2(10) = (/1,2,3,1,1,2,2,3,3,2/), &
        oind3(10) = (/1,2,3,2,3,1,3,1,2,3/), &
        hind1(15) = (/1,2,3,1,1,2,2,3,3,1,1,2,1,2,3/), &
        hind2(15) = (/1,2,3,1,1,2,2,3,3,1,1,2,1,2,3/), &
        hind3(15) = (/1,2,3,1,1,2,2,3,3,2,3,3,2,1,1/), &
        hind4(15) = (/1,2,3,2,3,1,3,1,2,2,3,3,3,3,2/)

! CHARGE case
!     Mpole_xy(1,1) = 1.d0
  DMpole_xy(1,1) = 0.d0
  if ( order .eq. 1 )return
  do j = 2,order
    DMpole_xy(1,j) = 0.d0
    DMpole_xy(j,1) = 0.d0
  enddo
! DIPOLES
  do j = 2,4
    do k = 2,4
      DMpole_xy(j,k) = DA_xy(j-1,k-1)
    enddo
  enddo
  if ( order .eq. 4 )return
  do j = 2,4
    do k = 5,order
      DMpole_xy(j,k) = 0.d0
      DMpole_xy(k,j) = 0.d0
    enddo
  enddo
! QUADRUPOLES
  do ind1 = 1,3
    k = qind1(ind1)
    do ind2 = 1,6
      i = qind1(ind2)
      j = qind2(ind2)
      DMpole_xy(ind1+4,ind2+4) = DA_xy(k,i)*A_xy(k,j) + A_xy(k,i)*DA_xy(k,j)
    enddo
  enddo
  do ind1 = 4,6
    k = qind1(ind1)
    l = qind2(ind1)
    do ind2 = 1,6
      i = qind1(ind2)
      j = qind2(ind2)
      DMpole_xy(ind1+4,ind2+4) =     &
            DA_xy(k,i)*A_xy(l,j) + DA_xy(k,j)*A_xy(l,i) +    &
             A_xy(k,i)*DA_xy(l,j) + A_xy(k,j)*DA_xy(l,i)
    enddo
  enddo
  if ( order .eq. 10 )return
  do j = 5,10
    do k = 11,order
      DMpole_xy(k,j) = 0.d0
      DMpole_xy(j,k) = 0.d0
    enddo
  enddo
! OCTUPOLES
! O'_lll
  do ind1 = 1,3
    l = oind1(ind1)
    do ind2 = 1,10
      i = oind1(ind2)
      j = oind2(ind2)
      k = oind3(ind2)
      DMpole_xy(ind1+10,ind2+10) =    &
                DA_xy(l,i)*A_xy(l,j)*A_xy(l,k) +    &
                A_xy(l,i)*DA_xy(l,j)*A_xy(l,k) +    &
                A_xy(l,i)*A_xy(l,j)*DA_xy(l,k) 
    enddo
  enddo
  do ind1 = 4,9
    l = oind1(ind1)
    m = oind3(ind1)
    do ind2 = 1,9
      i = oind1(ind2)
      k = oind3(ind2)
      DMpole_xy(ind1+10,ind2+10) =   &
           DA_xy(l,i)*A_xy(l,i)*A_xy(m,k) +     &
           2.d0 * DA_xy(l,i)*A_xy(m,i)*A_xy(l,k)
      DMpole_xy(ind1+10,ind2+10) = DMpole_xy(ind1+10,ind2+10) +   &
           A_xy(l,i)*DA_xy(l,i)*A_xy(m,k) +   &
           2.d0 * A_xy(l,i)*DA_xy(m,i)*A_xy(l,k)
      DMpole_xy(ind1+10,ind2+10) = DMpole_xy(ind1+10,ind2+10) +    &
           A_xy(l,i)*A_xy(l,i)*DA_xy(m,k) +    &
           2.d0 * A_xy(l,i)*A_xy(m,i)*DA_xy(l,k)
    enddo
    DMpole_xy(ind1+10,20) =    &
           DA_xy(l,1)*A_xy(l,2)*A_xy(m,3) +    &
           DA_xy(l,1)*A_xy(m,2)*A_xy(l,3) +    &
           DA_xy(m,1)*A_xy(l,2)*A_xy(l,3)
    DMpole_xy(ind1+10,20) = DMpole_xy(ind1+10,20) +   &
           A_xy(l,1)*DA_xy(l,2)*A_xy(m,3) +   &
           A_xy(l,1)*DA_xy(m,2)*A_xy(l,3) +   &
           A_xy(m,1)*DA_xy(l,2)*A_xy(l,3)
    DMpole_xy(ind1+10,20) = DMpole_xy(ind1+10,20) +   &
           A_xy(l,1)*A_xy(l,2)*DA_xy(m,3) +   &
           A_xy(l,1)*A_xy(m,2)*DA_xy(l,3) +   &
           A_xy(m,1)*A_xy(l,2)*DA_xy(l,3)
  enddo
  DMpole_xy(20,11) = 6.d0*DA_xy(1,1)*A_xy(2,1)*A_xy(3,1)
  DMpole_xy(20,11) = DMpole_xy(20,11) + 6.d0*A_xy(1,1)*DA_xy(2,1)*A_xy(3,1)
  DMpole_xy(20,11) = DMpole_xy(20,11) + 6.d0*A_xy(1,1)*A_xy(2,1)*DA_xy(3,1)
  DMpole_xy(20,12) = 6.d0*DA_xy(1,2)*A_xy(2,2)*A_xy(3,2)
  DMpole_xy(20,12) = DMpole_xy(20,12) + 6.d0*A_xy(1,2)*DA_xy(2,2)*A_xy(3,2)
  DMpole_xy(20,12) = DMpole_xy(20,12) + 6.d0*A_xy(1,2)*A_xy(2,2)*DA_xy(3,2)
  DMpole_xy(20,13) = 6.d0*DA_xy(1,3)*A_xy(2,3)*A_xy(3,3)
  DMpole_xy(20,13) = DMpole_xy(20,13) + 6.d0*A_xy(1,3)*DA_xy(2,3)*A_xy(3,3)
  DMpole_xy(20,13) = DMpole_xy(20,13) + 6.d0*A_xy(1,3)*A_xy(2,3)*DA_xy(3,3)
  do ind2 = 4,9
    i = oind1(ind2)
    k = oind3(ind2)
    DMpole_xy(20,10+ind2) = 2.d0*   &
           (  DA_xy(1,i)*A_xy(2,i)*A_xy(3,k) +    &
              DA_xy(1,i)*A_xy(2,k)*A_xy(3,i) +    &
              DA_xy(1,k)*A_xy(2,i)*A_xy(3,i) )
    DMpole_xy(20,10+ind2) = DMpole_xy(20,10+ind2) + 2.d0*   &
           (  A_xy(1,i)*DA_xy(2,i)*A_xy(3,k) +    &
              A_xy(1,i)*DA_xy(2,k)*A_xy(3,i) +    &
              A_xy(1,k)*DA_xy(2,i)*A_xy(3,i) )
    DMpole_xy(20,10+ind2) = DMpole_xy(20,10+ind2) + 2.d0*    &
           (  A_xy(1,i)*A_xy(2,i)*DA_xy(3,k) +    &
              A_xy(1,i)*A_xy(2,k)*DA_xy(3,i) +   &
              A_xy(1,k)*A_xy(2,i)*DA_xy(3,i) )
  enddo
  DMpole_xy(20,20) =    &
              DA_xy(1,1)*A_xy(2,2)*A_xy(3,3) +   &
              DA_xy(1,1)*A_xy(3,2)*A_xy(2,3) +   &
              DA_xy(2,1)*A_xy(1,2)*A_xy(3,3) +   &
              DA_xy(2,1)*A_xy(3,2)*A_xy(1,3) +   &
              DA_xy(3,1)*A_xy(1,2)*A_xy(2,3) +   &
              DA_xy(3,1)*A_xy(2,2)*A_xy(1,3)
  DMpole_xy(20,20) = DMpole_xy(20,20) +    &
              A_xy(1,1)*DA_xy(2,2)*A_xy(3,3) +   &
              A_xy(1,1)*DA_xy(3,2)*A_xy(2,3) +   &
              A_xy(2,1)*DA_xy(1,2)*A_xy(3,3) +   &
              A_xy(2,1)*DA_xy(3,2)*A_xy(1,3) +   &
              A_xy(3,1)*DA_xy(1,2)*A_xy(2,3) +   &
              A_xy(3,1)*DA_xy(2,2)*A_xy(1,3)
  DMpole_xy(20,20) = DMpole_xy(20,20) +   &
              A_xy(1,1)*A_xy(2,2)*DA_xy(3,3) +   &
              A_xy(1,1)*A_xy(3,2)*DA_xy(2,3) +   &
              A_xy(2,1)*A_xy(1,2)*DA_xy(3,3) +   &
              A_xy(2,1)*A_xy(3,2)*DA_xy(1,3) +   &
              A_xy(3,1)*A_xy(1,2)*DA_xy(2,3) +   &
              A_xy(3,1)*A_xy(2,2)*DA_xy(1,3)
  if ( order .eq. 20 )return
  do j = 11,20
    do k = 21,order
      DMpole_xy(k,j) = 0.d0
      DMpole_xy(j,k) = 0.d0
    enddo
  enddo

! HEXADECAPOLES
! H'_mmmm
  do ind1 = 1,3
    m = hind1(ind1)
    do ind2 = 1,15
      i = hind1(ind2)
      j = hind2(ind2)
      k = hind3(ind2)
      l = hind4(ind2)
      DMpole_xy(ind1+20,ind2+20) =   &
               DA_xy(m,i)*A_xy(m,j)*A_xy(m,k)*A_xy(m,l)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*DA_xy(m,j)*A_xy(m,k)*A_xy(m,l)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,j)*DA_xy(m,k)*A_xy(m,l)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,j)*A_xy(m,k)*DA_xy(m,l)
    enddo
  enddo
! H'_mmmn
  do ind1 = 4,9
    m = hind1(ind1)
    n = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
      DMpole_xy(ind1+20,ind2+20) =   &
               DA_xy(m,i)*A_xy(m,i)*A_xy(m,i)*A_xy(n,j) + 3.d0*   &
               DA_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(m,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*DA_xy(m,i)*A_xy(m,i)*A_xy(n,j) + 3.d0*  &
               A_xy(m,i)*DA_xy(m,i)*A_xy(n,i)*A_xy(m,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,i)*DA_xy(m,i)*A_xy(n,j) + 3.d0*  &
               A_xy(m,i)*A_xy(m,i)*DA_xy(n,i)*A_xy(m,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,i)*DA_xy(n,j) + 3.d0*  &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*DA_xy(m,j)
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      DMpole_xy(ind1+20,ind2+20) =   &
               DA_xy(m,i)*A_xy(m,i)*A_xy(m,j)*A_xy(n,k) +  &
               DA_xy(m,i)*A_xy(m,i)*A_xy(m,k)*A_xy(n,j) + 2.d0*  &
               DA_xy(m,i)*A_xy(m,j)*A_xy(m,k)*A_xy(n,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*DA_xy(m,i)*A_xy(m,j)*A_xy(n,k) +  &
               A_xy(m,i)*DA_xy(m,i)*A_xy(m,k)*A_xy(n,j) + 2.d0*  &
               A_xy(m,i)*DA_xy(m,j)*A_xy(m,k)*A_xy(n,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,i)*DA_xy(m,j)*A_xy(n,k) +  &
               A_xy(m,i)*A_xy(m,i)*DA_xy(m,k)*A_xy(n,j) + 2.d0*  &
               A_xy(m,i)*A_xy(m,j)*DA_xy(m,k)*A_xy(n,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +  &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,j)*DA_xy(n,k) +  &
               A_xy(m,i)*A_xy(m,i)*A_xy(m,k)*DA_xy(n,j) + 2.d0*  &
               A_xy(m,i)*A_xy(m,j)*A_xy(m,k)*DA_xy(n,i)
    enddo
  enddo

! H'_mmnn
  do ind1 = 10,12
    m = hind1(ind1)
    n = hind3(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
      DMpole_xy(ind1+20,ind2+20) = 3.d0 * (   &
               DA_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(n,j) +   &
               DA_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(n,i)  )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
          3.d0 * (                                                &
               A_xy(m,i)*DA_xy(m,i)*A_xy(n,i)*A_xy(n,j) +    &
               A_xy(m,i)*DA_xy(m,j)*A_xy(n,i)*A_xy(n,i)  )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
          3.d0 * (                                                &
               A_xy(m,i)*A_xy(m,i)*DA_xy(n,i)*A_xy(n,j) +    &
               A_xy(m,i)*A_xy(m,j)*DA_xy(n,i)*A_xy(n,i)  )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
          3.d0 * (                                                &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*DA_xy(n,j) +    &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*DA_xy(n,i)  )
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      DMpole_xy(ind1+20,ind2+20) =    &
               DA_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(n,k) +    &
               DA_xy(m,j)*A_xy(m,k)*A_xy(n,i)*A_xy(n,i) +  2.d0*    &
               DA_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(n,k) +  2.d0*    &
               DA_xy(m,i)*A_xy(m,k)*A_xy(n,i)*A_xy(n,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
               A_xy(m,i)*DA_xy(m,i)*A_xy(n,j)*A_xy(n,k) +    &
               A_xy(m,j)*DA_xy(m,k)*A_xy(n,i)*A_xy(n,i) +  2.d0*    &
               A_xy(m,i)*DA_xy(m,j)*A_xy(n,i)*A_xy(n,k) +  2.d0*    &
               A_xy(m,i)*DA_xy(m,k)*A_xy(n,i)*A_xy(n,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
               A_xy(m,i)*A_xy(m,i)*DA_xy(n,j)*A_xy(n,k) +    &
               A_xy(m,j)*A_xy(m,k)*DA_xy(n,i)*A_xy(n,i) +  2.d0*   &
               A_xy(m,i)*A_xy(m,j)*DA_xy(n,i)*A_xy(n,k) +  2.d0*   &
               A_xy(m,i)*A_xy(m,k)*DA_xy(n,i)*A_xy(n,j)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*DA_xy(n,k) +    &
               A_xy(m,j)*A_xy(m,k)*A_xy(n,i)*DA_xy(n,i) +  2.d0*   &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*DA_xy(n,k) +  2.d0*   &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,i)*DA_xy(n,j)
    enddo
  enddo

! H'_mmnp
  do ind1 = 13,15
    m = hind1(ind1)
    n = hind3(ind1)
    p = hind4(ind1)
! can put iiii & iiij together
    do ind2 = 1,9
      i = hind1(ind2)
      j = hind4(ind2)
      DMpole_xy(ind1+20,ind2+20) =    &
              3.d0*DA_xy(m,i)*A_xy(m,i)*A_xy(n,i)*A_xy(p,j) +    &
              3.d0*DA_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(p,i) +    &
              6.d0*DA_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(p,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
              3.d0*A_xy(m,i)*DA_xy(m,i)*A_xy(n,i)*A_xy(p,j) +    &
              3.d0*A_xy(m,i)*DA_xy(m,i)*A_xy(n,j)*A_xy(p,i) +    &
              6.d0*A_xy(m,i)*DA_xy(m,j)*A_xy(n,i)*A_xy(p,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
              3.d0*A_xy(m,i)*A_xy(m,i)*DA_xy(n,i)*A_xy(p,j) +    &
              3.d0*A_xy(m,i)*A_xy(m,i)*DA_xy(n,j)*A_xy(p,i) +    &
              6.d0*A_xy(m,i)*A_xy(m,j)*DA_xy(n,i)*A_xy(p,i)
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
              3.d0*A_xy(m,i)*A_xy(m,i)*A_xy(n,i)*DA_xy(p,j) +    &
              3.d0*A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*DA_xy(p,i) +    &
              6.d0*A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*DA_xy(p,i)
    enddo
! can put iijj & iijk together
    do ind2 = 10,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      DMpole_xy(ind1+20,ind2+20) =    &
               DA_xy(m,i)*A_xy(m,i)*A_xy(n,j)*A_xy(p,k) +    &
               DA_xy(m,i)*A_xy(m,i)*A_xy(n,k)*A_xy(p,j) + 2.d0 * (   &
               DA_xy(m,i)*A_xy(m,j)*A_xy(n,i)*A_xy(p,k) +   &
               DA_xy(m,i)*A_xy(m,j)*A_xy(n,k)*A_xy(p,i) +   &
               DA_xy(m,i)*A_xy(m,k)*A_xy(n,i)*A_xy(p,j) +   &
               DA_xy(m,i)*A_xy(m,k)*A_xy(n,j)*A_xy(p,i) +   &
               DA_xy(m,j)*A_xy(m,k)*A_xy(n,i)*A_xy(p,i) )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +   &
               A_xy(m,i)*DA_xy(m,i)*A_xy(n,j)*A_xy(p,k) +    &
               A_xy(m,i)*DA_xy(m,i)*A_xy(n,k)*A_xy(p,j) + 2.d0 * (    &
               A_xy(m,i)*DA_xy(m,j)*A_xy(n,i)*A_xy(p,k) +    &
               A_xy(m,i)*DA_xy(m,j)*A_xy(n,k)*A_xy(p,i) +    &
               A_xy(m,i)*DA_xy(m,k)*A_xy(n,i)*A_xy(p,j) +    &
               A_xy(m,i)*DA_xy(m,k)*A_xy(n,j)*A_xy(p,i) +    &
               A_xy(m,j)*DA_xy(m,k)*A_xy(n,i)*A_xy(p,i) )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
               A_xy(m,i)*A_xy(m,i)*DA_xy(n,j)*A_xy(p,k) +    &
               A_xy(m,i)*A_xy(m,i)*DA_xy(n,k)*A_xy(p,j) + 2.d0 * (    &
               A_xy(m,i)*A_xy(m,j)*DA_xy(n,i)*A_xy(p,k) +    &
               A_xy(m,i)*A_xy(m,j)*DA_xy(n,k)*A_xy(p,i) +    &
               A_xy(m,i)*A_xy(m,k)*DA_xy(n,i)*A_xy(p,j) +    &
               A_xy(m,i)*A_xy(m,k)*DA_xy(n,j)*A_xy(p,i) +    &
               A_xy(m,j)*A_xy(m,k)*DA_xy(n,i)*A_xy(p,i) )
      DMpole_xy(ind1+20,ind2+20) = DMpole_xy(ind1+20,ind2+20) +    &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,j)*DA_xy(p,k) +    &
               A_xy(m,i)*A_xy(m,i)*A_xy(n,k)*DA_xy(p,j) + 2.d0 * (    &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,i)*DA_xy(p,k) +    &
               A_xy(m,i)*A_xy(m,j)*A_xy(n,k)*DA_xy(p,i) +    &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,i)*DA_xy(p,j) +    &
               A_xy(m,i)*A_xy(m,k)*A_xy(n,j)*DA_xy(p,i) +    &
               A_xy(m,j)*A_xy(m,k)*A_xy(n,i)*DA_xy(p,i) )
    enddo
  enddo
  return
end subroutine XFORM_MPOLE_deriv_matrix
!------------------------------------------------------------
subroutine XFORM_MPOLE_field_matrix(A_xy,Field_xy,order)
! A_xy is matrix of first order terms connecting field in y-coord frame
! to those in the x-frame
! Field_xy expands this to get the higher order conversions between coord sys
! order is order of multipoles 1 for charge, 4 for charge-dipoles
! 10 for charge,dipole,quadrupole, 20 for up to octupole and
! finally 35 for up to hexadecapole
!--------------------------------------------------------------
! Field terms in y-frame are dphi_dy1, dphi_dy2, dphi_dy3
!       those in x-frame are dphi_dx1, dphi_dx2, dphi_dx3
! Thus:  A_xy(1,1) = dy1_dx1, A_xy(1,2) = dy2_dx1, A_xy(1,3) = dy3_dx1
!        A_xy(2,1) = dy1_dx2, A_xy(2,2) = dy2_dx2, A_xy(2,3) = dy3_dx2
!        A_xy(3,1) = dy1_dx3, A_xy(3,2) = dy2_dx3, A_xy(3,3) = dy3_dx3
! For example to get 2nd order derivs d^2_dx1_dx2 we expand
!  (d_dy1*dy1_dx1 + d_dy2*dy2_dx1 + d_dy3*dy3_dx1)*
!  (d_dy1*dy1_dx2 + d_dy2*dy2_dx2 + d_dy3*dy3_dx2)
! and collect terms involving each d^2_dyk_dyl
!------------------------------------------------------------ 
  implicit none
  integer order
  _REAL_ A_xy(3,3), Field_xy(order,order)
  integer i,j,k,l,m,n,p,ind1,ind2

  integer, parameter :: &
        qind1(6)  = (/1,2,3,1,1,2/), &
        qind2(6)  = (/1,2,3,2,3,3/), &
        oind1(10) = (/1,2,3,1,1,2,2,3,3,1/), &
        oind2(10) = (/1,2,3,1,1,2,2,3,3,2/), &
        oind3(10) = (/1,2,3,2,3,1,3,1,2,3/), &
        hind1(15) = (/1,2,3,1,1,2,2,3,3,1,1,2,1,2,3/), &
        hind2(15) = (/1,2,3,1,1,2,2,3,3,1,1,2,1,2,3/), &
        hind3(15) = (/1,2,3,1,1,2,2,3,3,2,3,3,2,1,1/), &
        hind4(15) = (/1,2,3,2,3,1,3,1,2,2,3,3,3,3,2/)

! CHARGE case
  Field_xy(1,1) = 1.d0
  if ( order .eq. 1 )return
  do j = 2,order
    Field_xy(1,j) = 0.d0
    Field_xy(j,1) = 0.d0
  enddo
! DIPOLES
  do j = 2,4
    do k = 2,4
      Field_xy(j,k) = A_xy(j-1,k-1)
    enddo
  enddo
  if ( order .eq. 4 )return
  do j = 2,4
    do k = 5,order
      Field_xy(j,k) = 0.d0
      Field_xy(k,j) = 0.d0
    enddo
  enddo
! QUADRUPOLES
  do ind1 = 1,3
    l = qind1(ind1)
    do ind2 = 1,3
      i = qind1(ind2)
      Field_xy(ind1+4,ind2+4) = A_xy(l,i)*A_xy(l,i)
    enddo
    do ind2 = 4,6
      i = qind1(ind2)
      j = qind2(ind2)
      Field_xy(ind1+4,ind2+4) = 2.d0*A_xy(l,i)*A_xy(l,j)
    enddo
  enddo
  do ind1 = 4,6
    l = qind1(ind1)
    m = qind2(ind1)
    do ind2 = 1,3
      i = qind1(ind2)
      Field_xy(ind1+4,ind2+4) = A_xy(l,i)*A_xy(m,i)
    enddo
    do ind2 = 4,6
      i = qind1(ind2)
      j = qind2(ind2)
      Field_xy(ind1+4,ind2+4) = (A_xy(l,i)*A_xy(m,j)+A_xy(m,i)*A_xy(l,j))
    enddo
  enddo
  if ( order .eq. 10 )return
  do j = 5,10
    do k = 11,order
      Field_xy(k,j) = 0.d0
      Field_xy(j,k) = 0.d0
    enddo
  enddo
! OCTUPOLES level field
! F'_lll
  do ind1 = 1,3
    l = oind1(ind1)
    do ind2 = 1,3
      i = oind1(ind2)
      Field_xy(ind1+10,ind2+10) = A_xy(l,i)**3
    enddo
    do ind2 = 4,9
      i = oind1(ind2)
      j = oind3(ind2)
      Field_xy(ind1+10,ind2+10) = 3.d0*A_xy(l,i)**2 * A_xy(l,j)
    enddo
    Field_xy(ind1+10,20) = 6.d0*A_xy(l,1)*A_xy(l,2)*A_xy(l,3) 
  enddo
! F'_llm
  do ind1 = 4,9
    l = oind1(ind1)
    m = oind3(ind1)
    do ind2 = 1,3
      i = oind1(ind2)
      Field_xy(ind1+10,ind2+10) = A_xy(l,i)**2 * A_xy(m,i)
    enddo
    do ind2 = 4,9
      i = oind1(ind2)
      j = oind3(ind2)
      Field_xy(ind1+10,ind2+10) = A_xy(l,i)**2 * A_xy(m,j) +   &
                                2.d0*A_xy(l,i)*A_xy(l,j)*A_xy(m,i)
    enddo
    Field_xy(ind1+10,20) = 2.d0*(A_xy(l,1)*A_xy(l,2)*A_xy(m,3)+   &
                                     A_xy(l,1)*A_xy(m,2)*A_xy(l,3)+   &
                                     A_xy(m,1)*A_xy(l,2)*A_xy(l,3) )
  enddo
! F'_123
  do ind2 = 1,3
    i = oind1(ind2)
    Field_xy(20,ind2+10) = A_xy(1,i)*A_xy(2,i)*A_xy(3,i)
  enddo
  do ind2 = 4,9
    i = oind1(ind2)
    j = oind3(ind2)
    Field_xy(20,ind2+10) = A_xy(1,i)*A_xy(2,i)*A_xy(3,j) +  &
                              A_xy(1,i)*A_xy(2,j)*A_xy(3,i) +  &
                              A_xy(1,j)*A_xy(2,i)*A_xy(3,i) 
  enddo
  Field_xy(20,20) = A_xy(1,1)*A_xy(2,2)*A_xy(3,3) +   &
                        A_xy(1,1)*A_xy(2,3)*A_xy(3,2) +   &
                        A_xy(1,2)*A_xy(2,1)*A_xy(3,3) +   &
                        A_xy(1,2)*A_xy(2,3)*A_xy(3,1) +   &
                        A_xy(1,3)*A_xy(2,1)*A_xy(3,2) +   &
                        A_xy(1,3)*A_xy(2,2)*A_xy(3,1)
  if ( order .eq. 20 )return
  do j = 11,20
    do k = 21,order
      Field_xy(k,j) = 0.d0
      Field_xy(j,k) = 0.d0
    enddo
  enddo
! HEXADECAPOLES
! F'_mmmm
  do ind1 = 1,3
    m = hind1(ind1)
    do ind2 = 1,3
      i = hind1(ind2)
      Field_xy(ind1+20,ind2+20) = A_xy(m,i)**4
    enddo
    do ind2 = 4,9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) = 4.d0 * A_xy(m,i)**3 * A_xy(m,j)
    enddo
    do ind2 = 10,12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) = 6.d0 * A_xy(m,i)**2 * A_xy(m,j)**2
    enddo
    do ind2 = 13,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1+20,ind2+20) = 12.d0*A_xy(m,i)**2 * A_xy(m,j)*A_xy(m,k)
    enddo
  enddo
  do ind1 = 4,9
    m = hind1(ind1)
    n = hind4(ind1)
    do ind2 = 1,3
      i = hind1(ind2)
      Field_xy(ind1+20,ind2+20) = A_xy(m,i)**3 * A_xy(n,i)
    enddo
    do ind2 = 4,9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) = A_xy(m,i)**3 * A_xy(n,j) +  &
                               3.d0*A_xy(m,i)**2 * A_xy(m,j)*A_xy(n,i)
    enddo
    do ind2 = 10,12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
              3.d0 * A_xy(m,i)**2 * A_xy(m,j) * A_xy(n,j) +    &
              3.d0 * A_xy(m,j)**2 * A_xy(m,i) * A_xy(n,i) 
    enddo
    do ind2 = 13,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
              6.d0 * A_xy(m,i) * A_xy(m,j) * A_xy(m,k) * A_xy(n,i) +  &
              3.d0 * A_xy(m,i)**2 * A_xy(m,j) * A_xy(n,k) +  &
              3.d0 * A_xy(m,i)**2 * A_xy(m,k) * A_xy(n,j)
    enddo
  enddo
  do ind1 = 10,12
    m = hind1(ind1)
    n = hind4(ind1)
    do ind2 = 1,3
      i = hind1(ind2)
      Field_xy(ind1+20,ind2+20) = A_xy(m,i)**2 * A_xy(n,i)**2
    enddo
    do ind2 = 4,9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
                 2.d0 * A_xy(m,i)**2 * A_xy(n,i) * A_xy(n,j) +   &
                 2.d0 * A_xy(n,i)**2 * A_xy(m,i) * A_xy(m,j)
    enddo
    do ind2 = 10,12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
                       A_xy(m,i)**2 * A_xy(n,j)**2 +   &
                       A_xy(m,j)**2 * A_xy(n,i)**2 +   &
               4.d0 *  A_xy(m,i) * A_xy(m,j) *A_xy(n,i) * A_xy(n,j)
    enddo
    do ind2 = 13,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =   &
                 2.d0 * A_xy(m,i)**2 * A_xy(n,j) * A_xy(n,k) +  &
                 2.d0 * A_xy(n,i)**2 * A_xy(m,j) * A_xy(m,k) +  &
           4.d0 * A_xy(m,i) * A_xy(m,j) * A_xy(n,i) * A_xy(n,k) +  &
           4.d0 * A_xy(m,i) * A_xy(m,k) * A_xy(n,i) * A_xy(n,j) 
    enddo
  enddo
  do ind1 = 13,15
    m = hind1(ind1)
    n = hind3(ind1)
    p = hind4(ind1)
    do ind2 = 1,3
      i = hind1(ind2)
      Field_xy(ind1+20,ind2+20) = A_xy(m,i)**2 * A_xy(n,i) * A_xy(p,i)
    enddo
    do ind2 = 4,9
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
            2.d0 * A_xy(m,i) * A_xy(m,j) * A_xy(n,i) * A_xy(p,i) +  &
            A_xy(m,i)**2 * (A_xy(n,i)*A_xy(p,j) + A_xy(n,j)*A_xy(p,i))
    enddo
    do ind2 = 10,12
      i = hind1(ind2)
      j = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =    &
               A_xy(m,i)**2 * A_xy(n,j) * A_xy(p,j) +   &
               A_xy(m,j)**2 * A_xy(n,i) * A_xy(p,i) +   &
               2.d0 * A_xy(m,i) * A_xy(m,j) *   &
                   (A_xy(n,i)*A_xy(p,j) + A_xy(n,j)*A_xy(p,i))
    enddo
    do ind2 = 13,15
      i = hind1(ind2)
      j = hind3(ind2)
      k = hind4(ind2)
      Field_xy(ind1+20,ind2+20) =   &
           A_xy(m,i)**2 * (A_xy(n,j)*A_xy(p,k) + A_xy(n,k)*A_xy(p,j)) +   &
                      2.d0 * A_xy(m,i)*A_xy(m,j) *    &
                      (A_xy(n,i)*A_xy(p,k) + A_xy(n,k)*A_xy(p,i)) +   &
                      2.d0 * A_xy(m,i)*A_xy(m,k) *    &
                      (A_xy(n,i)*A_xy(p,j) + A_xy(n,j)*A_xy(p,i)) +   &
                   2.d0 * A_xy(m,j) * A_xy(m,k) * A_xy(n,i) * A_xy(p,i)
    enddo
  enddo

  return
end subroutine XFORM_MPOLE_field_matrix
!------------------------------------------------------------
subroutine XFORM_MPOLE(Mpole_xy,dimxy,Mpole_in,Mpole_out,order)
  implicit none
  integer order,dimxy
  _REAL_ Mpole_xy(dimxy,dimxy)
  _REAL_ Mpole_in(*),Mpole_out(*)

  integer i,j

  if ( order .eq. 0 )return
  Mpole_out(1) = Mpole_xy(1,1)*Mpole_in(1)
  if ( order .eq. 1 )return
! DIPOLES
  do i = 2,4
    Mpole_out(i) = 0.d0
    do j = 2,4
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i,j)*Mpole_in(j)
    enddo
  enddo
  if ( order .eq. 4 )return
! QUADRUPOLES
  do i = 5,10
    Mpole_out(i) = 0.d0
    do j = 5,10
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i,j)*Mpole_in(j)
    enddo
  enddo
  if ( order .eq. 10 )return
! OCTUPOLES
  do i = 11,20
    Mpole_out(i) = 0.d0
    do j = 11,20
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i,j)*Mpole_in(j)
    enddo
  enddo
  if ( order .eq. 20 )return
! HEXADECAPOLES
  do i = 21,35
    Mpole_out(i) = 0.d0
    do j = 21,35
      Mpole_out(i) = Mpole_out(i) + Mpole_xy(i,j)*Mpole_in(j)
    enddo
  enddo
  if ( order .eq. 35 )return
  return
end subroutine XFORM_MPOLE
!------------------------------------------------------------
subroutine XFORM_FIELD(Field_xy,dimxy,Field_in,Field_out,order)
  implicit none
  integer order,dimxy
  _REAL_ Field_xy(dimxy,dimxy)
  _REAL_ Field_in(*),Field_out(*)

  integer i,j
  !Field_out is xform of Field_in

  if ( order .eq. 0 )return
  Field_out(1) = Field_xy(1,1)*Field_in(1)
  if ( order .eq. 1 )return
! DIPOLES
  do i = 2,4
    Field_out(i) = 0.d0
    do j = 2,4
      Field_out(i) = Field_out(i) + Field_xy(i,j)*Field_in(j)
    enddo
  enddo
  if ( order .eq. 4 )return
! QUADRUPOLES
  do i = 5,10
    Field_out(i) = 0.d0
    do j = 5,10
      Field_out(i) = Field_out(i) + Field_xy(i,j)*Field_in(j)
    enddo
  enddo
  if ( order .eq. 10 )return
! OCTUPOLES
  do i = 11,20
    Field_out(i) = 0.d0
    do j = 11,20
      Field_out(i) = Field_out(i) + Field_xy(i,j)*Field_in(j)
    enddo
  enddo
  if ( order .eq. 20 )return
! HEXADECAPOLES
  do i = 21,35
    Field_out(i) = 0.d0
    do j = 21,35
      Field_out(i) = Field_out(i) + Field_xy(i,j)*Field_in(j)
    enddo
  enddo
  if ( order .eq. 35 )return
  return
end subroutine XFORM_FIELD
