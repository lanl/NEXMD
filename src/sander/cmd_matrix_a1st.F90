
#include "dprec.fh"


!------------------------------------------------------------------------------
subroutine trans_full( a, b, a_to_b )

! Transform cartesian positions into normal mode positions.
  use pimd_vars, only: nbead
  use full_pimd_vars, only: mybeadid,xall
  !..................................................
  implicit none

#include "les.h"
#include "memory.h"
 
  integer :: iatm3, ibead,ierr
  integer :: jatm3, jbead
  integer :: myfirst,mylast

#ifdef MPI
  include 'mpif.h'
# include "parallel.h"
#endif

  _REAL_, intent(in) :: a(3*natom),a_to_b(nbead,nbead)
  _REAL_  :: b(3,natom)
  !..................................................
#ifdef MPI
  if( sanderrank.eq.0) then
      call mpi_allgather(a   ,3*natom,MPI_DOUBLE_PRECISION,&
                         xall,3*natom,MPI_DOUBLE_PRECISION,&
                         commmaster,ierr)
  endif
  call mpi_bcast(xall,3*natom*nbead,MPI_DOUBLE_PRECISION,0,commsander,ierr)
#endif
  b(1:3,1:natom) = 0.d0

  do jbead = 1,nbead
     b(1:3,1:natom)=b(1:3,1:natom)+a_to_b(mybeadid,jbead)*xall(1:3,1:natom,jbead)
  end do

end subroutine


!------------------------------------------------------------------------------
subroutine trans_pos_nmode_to_cart_a1st( pos_nmode, x )
!------------------------------------------------------------------------------
  use pimd_vars, only: nbead,natomCL
  use cmd_vars, only: nmode_to_cart
  !..................................................
  implicit none
#include "memory.h"
  _REAL_, intent(in) :: pos_nmode(3*natomCL,nbead)
  _REAL_ :: x(3*natomCL,nbead)
#include "les.h"

#ifdef MPI
   include 'mpif.h'
#  include "parallel.h"
#endif

  integer :: iatm3, ibead,iatom
  integer :: jatm3, jbead
  integer :: itask, myfirst, mylast,ierr
  integer :: base, bead_per_node

#ifdef MPI
  bead_per_node = nbead/numtasks
  myfirst = bead_per_node*mytaskid+1
  mylast  = bead_per_node*mytaskid+bead_per_node
#else
  myfirst = 1
  mylast  = nbead
#endif

#ifdef MPI
  if( numtasks>1) call xdist(pos_nmode, xx(lfrctmp), natom)
#endif

  do ibead = myfirst,mylast
     x(1:3*natomCL,ibead)=0.d0
     do jbead = 1,nbead
        x(1:3*natomCL,ibead)=x(1:3*natomCL,ibead)+nmode_to_cart(ibead,jbead)*pos_nmode(1:3*natomCL,jbead)
     end do
  end do

end subroutine


!------------------------------------------------------------------------------
subroutine trans_force_cart_to_nmode_a1st( f, force_nmode )

! Transform cartesian forces into normal mode forces
!------------------------------------------------------------------------------

  use pimd_vars, only: nbead,natomCL

  use cmd_vars, only: fcart_to_fnmode

  !..................................................

  implicit none

#include "memory.h"
#include "les.h"

#ifdef MPI
# include "parallel.h"
  include 'mpif.h'
  integer :: bead_per_node,base,itask,ierr
#endif

  _REAL_, intent(in) :: f(3*natomCL,nbead)
  _REAL_ :: force_nmode(3*natomCL,nbead)
  !...........................xx(lfrctmp).......................

  integer :: iatm3, ibead,iatom
  integer :: jatm3, jbead,myfirst,mylast
  !..................................................

#ifdef MPI
  if(numtasks>1) call xdist(f, xx(lfrctmp), natom)
#endif

#ifdef MPI
  bead_per_node = nbead/numtasks
  myfirst = bead_per_node*mytaskid+1
  mylast  = bead_per_node*mytaskid+bead_per_node
#else
  myfirst = 1
  mylast  = nbead
#endif

  do ibead = myfirst,mylast
     force_nmode(1:3*natomCL,ibead)=0.d0
     do jbead = 1,nbead
        force_nmode(1:3*natomCL,ibead)=force_nmode(1:3*natomCL,ibead)+fcart_to_fnmode(ibead,jbead)*f(1:3*natomCL,jbead)
     end do
  end do

end subroutine


!------------------------------------------------------------------------------
subroutine trans_vel_nmode_to_cart_a1st( vel_nmode, v )

! Transform normal mode positions into cartesian positions.
!------------------------------------------------------------------------------

  use pimd_vars, only: nbead,natomCL

  use cmd_vars, only: nmode_to_cart

  !..................................................

  implicit none
#include "memory.h"
#include "les.h"
 
  _REAL_, intent(in) :: vel_nmode(3*natom)

  _REAL_ :: v(3*natom)

  !..................................................

  integer :: iatm3, ibead
  integer :: jatm3, jbead

  v(1:3*natom) = 0.d0

  iatm3 = 0
  do ibead = 1,nbead
     jatm3 = 0
     do jbead = 1,nbead
        v(iatm3+1:iatm3+3*natomCL)=v(iatm3+1:iatm3+3*natomCL)+nmode_to_cart(ibead,jbead)*vel_nmode(jatm3+1:jatm3+3*natomCL)
        jatm3=jatm3+3*natomCL
     end do
     iatm3=iatm3+3*natomCL
  end do

end subroutine


!------------------------------------------------------------------------------
subroutine trans_vel_cart_to_nmode_a1st( v, vel_nmode )
                                                                                
! Transform normal mode positions into cartesian positions.
!------------------------------------------------------------------------------
                                                                                
  use pimd_vars, only: nbead,natomCL
                                                                                
  use cmd_vars, only: cart_to_nmode
                                                                                
  !..................................................
                                                                                
  implicit none
#include "memory.h"
#include "les.h"
 
  _REAL_, intent(in) :: v( 3*natom )
                                                                                
  _REAL_ :: vel_nmode( 3*natom )
                                                                                
  !..................................................
                                                                                
  integer :: iatm3, ibead
  integer :: jatm3, jbead
 
  vel_nmode(1:3*natom) = 0.d0

  iatm3 = 0
  do ibead = 1,nbead
     jatm3 = 0
     do jbead = 1,nbead
        vel_nmode(iatm3+1:iatm3+3*natomCL)=vel_nmode(iatm3+1:iatm3+3*natomCL)+cart_to_nmode(ibead,jbead)*v(jatm3+1:jatm3+3*natomCL)
        jatm3=jatm3+3*natomCL
     end do
     iatm3=iatm3+3*natomCL
  end do

end subroutine

/*

!------------------------------------------------------------------------------
subroutine scale_vel_centroid ( vel_nmode, istart, iend )
!
! scale centroid velocities and fix total momentum equal to zero
!------------------------------------------------------------------------------

  use constants, only: kB
  
  use pimd_vars, only: natomCL
  
  use cmd_vars, only: mass_nmode
  
  !..................................................

  implicit none

#include "md.h"
#include "les.h"
#include "memory.h"
#ifdef MPI
#  include "parallel.h"
   include 'mpif.h'
#endif
  _REAL_ :: vel_nmode(3,natom)
  
  integer   iatom,ibead, istart, iend, ierr

  !..................................................

  integer :: iseed = -2209

  _REAL_ :: scale_factor
  
  _REAL_ :: E_kin, temp_ini
  
  _REAL_ :: p(3), ptmp(3)
  
  !..................................................
  
  !! Set total momentum equal to zero for path centroid.
  p(1:3)=0.0
  do iatom = istart, iend
     ibead = (iatom-1)/natomCL+1
     if(ibead.eq.1) then
       p(1) = p(1) + mass_nmode(iatom) * vel_nmode( 1,iatom )
       p(2) = p(2) + mass_nmode(iatom) * vel_nmode( 2,iatom )
       p(3) = p(3) + mass_nmode(iatom) * vel_nmode( 3,iatom )
     end if
  end do

#ifdef MPI
  call mpi_allreduce(p,ptmp,3,MPI_DOUBLE_PRECISION,mpi_sum,commsander,ierr)
  p = ptmp
#endif

  p(1) = p(1) / dble( natomCL )
  p(2) = p(2) / dble( natomCL )
  p(3) = p(3) / dble( natomCL )

  E_kin = 0.d0

  do iatom = istart, iend
     ibead = (iatom-1)/natomCL+1
     if(ibead.eq.1) then
        vel_nmode( 1,iatom) = vel_nmode( 1, iatom )  &
                                - p(1) / mass_nmode( iatom )

        vel_nmode( 2,iatom) = vel_nmode( 2,iatom )  &
                                - p(2) / mass_nmode( iatom )

        vel_nmode( 3,iatom) = vel_nmode( 3,iatom )  &
                                - p(3) / mass_nmode( iatom )

        E_kin = E_kin + 0.5d0 * mass_nmode( iatom )  &
                           * ( vel_nmode( 1, iatom )**2  &
                             + vel_nmode( 2, iatom )**2  &
                             + vel_nmode( 3, iatom )**2 )
     end if
  end do

  temp_ini = 2.d0 * E_kin / 3.d0 / kB / dble( natomCL )

  scale_factor = sqrt( temp0 / temp_ini )

  !! Scale velocity according to target temperature.
  E_kin = 0.d0

  p(1:3) = 0.d0
  do iatom = istart, iend
     ibead = (iatom-1)/natomCL+1
     if(ibead.eq.1) then
        vel_nmode( 1, iatom ) = vel_nmode( 1, iatom ) * scale_factor
        vel_nmode( 2, iatom ) = vel_nmode( 2, iatom ) * scale_factor
        vel_nmode( 3, iatom ) = vel_nmode( 3, iatom ) * scale_factor
     end if
  end do


end subroutine

*/

!#else


!  iatm3 = 0
!  do ibead = 1,nbead
!     jatm3 = 0
!     do jbead = 1,nbead
!        force_nmode(iatm3+1:iatm3+3*natomCL)=force_nmode(iatm3+1:iatm3+3*natomCL)+fcart_to_fnmode(ibead,jbead)*f(1:3*natomCL,jbead)
!        jatm3=jatm3+3*natomCL
!     end do
!     iatm3=iatm3+3*natomCL
!  end do


