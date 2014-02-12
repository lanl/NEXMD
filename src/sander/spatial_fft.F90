! <compile=optimized>
#include "assert.fh"
#include "dprec.fh"
!+ Parallel 3D FFT 
!
      Module fft
!  
! Description:
! Parallel 3D FFT implemented via a 1D x 1D x 1D data distribution
! also known as a 2D domain decomposition.
! The yz grid plane is partitioned among the processors.  1D FFTs
! in the x direction are calculated locally.  The data is repartitioned
! with respect to the zx grid plane and the y FFTs calculated.
! Finally, the z FFTs are calculated in the xy grid partition.
! The data is left in this representation for use by the scalar sum.
! The public routines include:
!     Subroutine backward_rc_fft( data )
!          3D Backward FFT with Real to Complex transform in x direction.
!     Subroutine fft_init( nfft1, nfft2, nfft3 )
!          Initialize parallel FFT data structures
!     Subroutine forward_rc_fft( data )
!          3D Forward FFT with Real to Complex transform in x direction.
!
 
#ifdef DPREC
#  define _COMPLEX_ double complex
#else
#  define _COMPLEX_ complex
#endif
#define _DATA_ _COMPLEX_

Implicit none

logical,save,public :: column_fft_flag
#ifdef MPI
Private
Public :: backward_rc_fft
Public :: forward_rc_fft
Public :: fft_init
Public :: get_fft_limits
Public :: YZ_X_PARTITION,XY_Z_PARTITION
Public :: get_xy_z_partition_limits
! The representation of grid space is based completely on closed intervals.
! Each processor gets a parallelepiped of grid points.  One direction
! is completely local and the other two directions are decomposed across
! processors. 
! 1D FFTs are calculated on the local vector.
! For best FFT performance that local direction is stored contiguously,
! that is, with stride one, in memory.
! Consequently, the data must be globally transposed between calculating
! the directional FFTs.
! An advantage of the 1D x 1D x 1D distribution is that all to all
! communications are unnecessary for the transpose.
! In repartitioning the data the lengths of the edges of the parallelepipeds
! in one direction are preserved.
! Communication is only between processors that have common closed intervals
! in the preserved direction.
! The representation of processor space is also based on closed intervals
! with the usual MPI rank from 0 to the number of processors - 1.
! However, the stride through these intervals may not be 1.
!
!     Type ClosedInterval
!       Integer             :: infimum   ! greatest lower bound
!       Integer             :: supremum  ! least upper bound
!     End Type ClosedInterval
!
!     Type Rectangle
!       Type( ClosedInterval )    :: first_coordinate
!       Type( ClosedInterval )    :: second_coordinate
!     End Type Rectangle
!
!     Type GridPartition
!       Type( Rectangle ), dimension( 0:NUMBER_PROCESSORS - 1 )   :: decomposed
!       Type( ClosedInterval )                                    :: contiguous
!     End Type GridPartition
!
!     Type TransposeCommunications
!       Type( ClosedInterval ), dimension( 0:NUMBER_PROCESSORS -1):: comm_group
!       Integer                                                   :: stride
!     End Type TransposeCommunications
!
! Three grid partitions exist: yz_x, zx_y, and xy_z, where
! the variables before the underscore indicate which directions are
! decomposed across processor space and the variable after the
! underscore indicates which direction is contiguous on each processor.
! (GridPartition % decomposed and TransposeCommunications % comm_group
! would be dynamically allocated.)
! Since Fortran 77 does not support user defined types and because
! the partitioning and transposition algorithms are independent of
! the actual grid directions, these data structures are mutated
! into multidimensional arrays.

integer,parameter :: NDIRECTIONS = 3 ! Number of Grid Directions
integer,parameter :: X_DIRECTION = 1 ! X Grid Direction
integer,parameter :: Y_DIRECTION = 2 ! Y Grid Direction
integer,parameter :: Z_DIRECTION = 3 ! Z Grid Direction

integer,dimension(1:NDIRECTIONS) :: fftdim,fftdatalen
integer,parameter :: CONTIGUOUS_INFIMUM = 0 
! Interval infimum for contiguous direction

integer,parameter :: NRECT_COORS=2      ! Number of Rectangle Coordinates
integer,parameter :: FIRST_RECT_COOR=1  ! Rectangle % FIRST_RECT_COORdinate
integer,parameter :: SECOND_RECT_COOR=2 ! Rectangle % SECOND_RECT_COORdinate

integer,parameter ::   NPARTITIONS    = 3 ! Number of Partitions
integer,parameter ::   YZ_X_PARTITION = 1 ! yz decomposed, x contiguous
integer,parameter ::   ZX_Y_PARTITION = 2 ! zx decomposed, y contiguous
integer,parameter ::   XY_Z_PARTITION = 3 ! xy decomposed, z contiguous

! this will become MPI_MAX_PROCESSORS from parallel.h
Integer, parameter      :: MAX_PROCESSORS = 128

integer   decomposed_inf( NRECT_COORS, 0:MAX_PROCESSORS - 1, NPARTITIONS)
integer   decomposed_sup( NRECT_COORS, 0:MAX_PROCESSORS - 1, NPARTITIONS)
integer   contiguous_inf( NPARTITIONS )
integer   contiguous_sup( NPARTITIONS )
integer   contiguous_fftdim( NPARTITIONS )

integer,parameter :: NCOMMGROUPS    = 2  ! Number of Communication Groups
integer,parameter :: YZ_X_COMM_ZX_Y = 1  ! yz_x partition <--> zx_y partition
integer,parameter :: XY_Z_COMM_ZX_Y = 2  ! xy_z partition <--> zx_y partition

integer,dimension( 0:MAX_PROCESSORS - 1, NCOMMGROUPS) :: comm_group_begin,comm_group_end
integer,dimension( NCOMMGROUPS) :: comm_group_stride

! These map communication groups to rectangle coordinates with 
! respect to being preserved or altered during transpositions.

integer,dimension( NCOMMGROUPS ) ::  preserved_coor,altered_coor

! Underlying the above logical representation of grid space is the
! physical layout in memory.  
! The contiguous direction in each partition is stride one in memory, 
! i.e., the fastest varying dimension.
! By convention, the slowest varying dimensions are the variables
! immediately before the underscore in the partition name:
! Partition   Slowest  Middle  Fastest
!     yz_x       z       y        x
!     zx_y       x       z        y
!     xy_z       y       x        z
! All grid addressing is through absolute global addresses, that is,
! ordered triples of integers each in the range CONTIGUOUS_INFIMUM to
! fftdim(?_DIRECTION), that are mapped into relative processor-local indices.
!
! Because the general algorithms were designed from the specific case
! of the yz_x to zx_y backward transpose and are based on the concept
! of preserved coordinates instead of slowest varying dimensions,
! a problem exists between the last transpose of the backward FFT
! and the first transpose of the forward FFT: the original generic
! transpose algorithm has the slowest varying dimension hard coded
! for the backward transpose through specific loop nestings and
! argument orderings.  Several methods could be used to pass in the
! information that loop nestings and argument orderings need to be reversed
! for the forward FFT.  For ease and speed of development another transpose
! algorithm is duplicated and hard coded for the forward FFT.
! To return to a single transpose routine a data structure that maintains
! the direction of the slowest varying dimension, eg,
!      integer   slowest_coor( NRECT_COORS, NPARTITIONS, NCOMMGROUPS )
! could be used remap the pre_c and alt_c into slow_coor and middle_coor.

integer,parameter ::   TRANSFORM_BACKWARD = -1, TRANSFORM_FORWARD = 1

_DATA_ , dimension(:),allocatable :: transposed_data,recv_buffer,send_buffer
integer,save :: fft_alloc_ier, mpi_ierr
_REAL_, dimension(:),allocatable,save :: &
                                 fft_table_1, fft_table_2, fft_table_3, &
                                 alpha_rcfft, beta_rcfft
_REAL_, dimension(:,:),allocatable,save :: tmp_rc

character(len=40) :: fmt0='(i4,1x,a,10i5)', &
                     fmt1='(2i5,2e15.5,10i5)'

!interface print_data
!   module procedure print_data_c,print_data_i,print_data_r
!end interface
interface 
   subroutine dumpq(a,n1,n2,n3,lun,p_type,n10,n20,n30)
     implicit none
     integer,intent(in) :: n1,n2,n3,lun
     integer,optional :: p_type,n10,n20,n30
     double complex,intent(in) :: a(n1,n2,n3)
   end subroutine dumpq
end interface

Contains
! public routines in alphabetical order... mostly (mfc added some out of order)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine get_fft_limits(partition, &
      gridmin0,gridmax0,gridmin1,gridmax1,gridmin2,gridmax2, &
      mytaskid  )
  implicit none
  integer,intent(in) :: mytaskid,partition
  integer,intent(out) :: gridmin0,gridmax0,gridmin1,gridmax1,gridmin2,gridmax2

  !---- get max/min from the decomposed inf and sup

  gridmin0 = contiguous_inf( partition )
  gridmax0 = contiguous_sup( partition )
  gridmin1 = decomposed_inf(1,mytaskid,partition)
  gridmax1 = decomposed_sup(1,mytaskid,partition)
  gridmin2 = decomposed_inf(2,mytaskid,partition)
  gridmax2 = decomposed_sup(2,mytaskid,partition)

end Subroutine get_fft_limits


  
  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            BACKWARD_RC_FFT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

Subroutine backward_rc_fft( data )
  !
  ! 3D Backward FFT with Real to Complex transform in x direction.
  !
  !------------------------------------------------------------------
  use trace
  implicit none
#  include "def_time.h"
#include "parallel.h"
  _DATA_  data(0:*)          ! original data,              intent(in)
   integer xgridmin,xgridmax,ygridmin,ygridmax,zgridmin,zgridmax
   
  call x_fft( TRANSFORM_BACKWARD, data )
  
  call timer_start(TIME_FFTCOMM)
  call transpose( data,           YZ_X_PARTITION, &
                  transposed_data,ZX_Y_PARTITION, &
                  YZ_X_COMM_ZX_Y, recv_buffer, send_buffer )
   call timer_stop(TIME_FFTCOMM)
  
  call y_fft( TRANSFORM_BACKWARD, transposed_data )
  
   call timer_start(TIME_FFTCOMM)
                                   ! transpose routine
  call transpose( transposed_data, ZX_Y_PARTITION, data, &
        XY_Z_PARTITION, XY_Z_COMM_ZX_Y, recv_buffer, send_buffer )
    call timer_stop(TIME_FFTCOMM)
 
  call z_fft( TRANSFORM_BACKWARD, data )

End subroutine backward_rc_fft
  
  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            FFT_INIT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine fft_init( nfft1, nfft2, nfft3 )
  !
  ! Initialize parallel FFT data structures including the 
  ! grid partitioning information.
  ! This will be merged into ew_setup.f which will be reorganized.
  !
  !------------------------------------------------------------------
  use trace
  use constants, only: TWOPI
  Implicit none
  Integer,intent(in) :: nfft1,nfft2,nfft3
# include "parallel.h"
#include "ew_parallel.h"
  integer, save  :: focus = 0  ! Tracing output focus
  
  Integer      bisections
  Integer      maxfftdim
  Integer      size_fft_table
  Integer      size_fft_work
  
  Integer      xgrades
  Integer      ygrades
  Integer      zgrades
  Integer,dimension( MAX_PROCESSORS ) :: xinf,xsup,yinf,ysup,zinf,zsup
  integer i
  _REAL_ theta,pi2n
  
    call trace_enter( 'fft_init' )
 if ( mod( nfft1, 2 ) /= 0 ) then
     call sander_bomb("fft_init in ew_fft.f", &
           "For RealComplex FFT","nfft1 must be even")
  endif

  fftdim(X_DIRECTION) = nfft1/2
  fftdim(Y_DIRECTION) = nfft2
  fftdim(Z_DIRECTION) = nfft3
  fftdatalen(X_DIRECTION) = nfft1/2+1
  fftdatalen(Y_DIRECTION) = nfft2
  fftdatalen(Z_DIRECTION) = nfft3
  maxfftdim = max( fftdim(X_DIRECTION), fftdim(Y_DIRECTION), fftdim(Z_DIRECTION) )
  
  !  Using pubfft zfft1D for 1D ffts
  !  These do not require work space
  size_fft_table = 3*( 4*maxfftdim + 15 )  ! cffti.doc
  size_fft_work  = 3*( 4*maxfftdim + 15 )  ! guess
  
  ! Grid planes will be decomposed into rectangular areas formed from
  ! the Cartesian product of closed intervals.
  ! This simple algorithm to create the intervals divides the axes
  ! based on the prime factorization of the number of processors.
  ! Only supports powers of two.
  ! Generalization to include powers of 3 and 5 is expected.
  ! New algorithm outline: compute the number of multiplicative partitions
  ! of length two of numtasks; enumerate these partitions calculating the
  ! difference between the partitions; choose partition with minimal difference.
  ! Refer to a number theory or group theory text for the additive analog.
  
  bisections = log( real( numtasks  ) ) / log( 2.0 ) / 2
  call Trace_integer( 'bisections is ',bisections )
  ygrades = 2 ** bisections
  zgrades = 2 ** bisections
  if ( mod( log( real( numtasks ) ) / log( 2.0 ), 2.0 ) /= 0 ) then
     ! perform the extra bisection on the longer grid dimension
     if ( fftdim(Y_DIRECTION) > fftdim(Z_DIRECTION) ) then
        ygrades = 2 * ygrades
     else
        zgrades = 2 * zgrades
     endif
  endif
  call Trace_integer( 'ygrades is ',ygrades )
  call Trace_integer( 'zgrades is ',zgrades )

  ! The communication groups depend on only the number of interval divisions,
  ! here denoted grades to keep name lengths short, in the 2D data decomposition.
  call create_comm_groups( ygrades, zgrades )
  
  call divide_interval( fftdim(Y_DIRECTION), ygrades, yinf, ysup )
  call divide_interval( fftdim(Z_DIRECTION), zgrades, zinf, zsup )
  
  
  ! The first 1D FFT, which is real to complex, will be in the x direction.
  ! The choice of the direction to FFT first is arbitrary, but has been
  ! determined higher up in the reciprocal Ewald algorithm (fill_charge_grid?).
  ! The y 1D complex to complex FFT follows then the z 1D complex to complex FFT.
  ! The directional sequence of FFTs determines the optimal decomposition
  ! of grid space onto the processors to minimize communications.
  ! Hence the yz grid plane is partitioned among processors first.
  
  ! These are arbitrary and have been chosen to minimize storage requirements.
  preserved_coor( YZ_X_COMM_ZX_Y ) = SECOND_RECT_COOR
  altered_coor( YZ_X_COMM_ZX_Y )   = FIRST_RECT_COOR
  preserved_coor( XY_Z_COMM_ZX_Y ) = FIRST_RECT_COOR
  altered_coor( XY_Z_COMM_ZX_Y )   = SECOND_RECT_COOR
  
  call create_partition( YZ_X_PARTITION, &
        altered_coor( YZ_X_COMM_ZX_Y ), ygrades, yinf, ysup, &
        preserved_coor( YZ_X_COMM_ZX_Y ), zgrades, zinf, zsup, &
        CONTIGUOUS_INFIMUM, fftdatalen(X_DIRECTION)-1,fftdim(X_DIRECTION) )
  
  ! Optimal communications dictate that the z partitioning is preserved.
  ! Hence the x grades for the zx partitioning equal the y grades from
  ! the yz partitioning, and the processors are mapped to the
  ! x grades before the z grades, ie, the x direction is inner.
  
  xgrades = ygrades
  ! Since this was a real-complex fft in x direction, the size
  ! of the data is one bigger than the fftdim, use fftdatalen
  call divide_interval( fftdatalen(X_DIRECTION), xgrades, xinf, xsup )
  call create_partition( ZX_Y_PARTITION,  &
        preserved_coor( XY_Z_COMM_ZX_Y ), xgrades, xinf, xsup, &
        preserved_coor( YZ_X_COMM_ZX_Y ), zgrades, zinf, zsup, &
        CONTIGUOUS_INFIMUM, fftdim(Y_DIRECTION) - 1,fftdim(Y_DIRECTION) )
  
  ! Optimal communications dictate that the x partitioning is preserved.
  ! Hence the y grades for the xy partitioning equal the z grades from
  ! the zx partitioning, and the processors are mapped to the
  ! y grades after the x grades, ie, the x direction is inner.
  
  ygrades = zgrades
  call divide_interval( fftdim(Y_DIRECTION), ygrades, yinf, ysup )
  call create_partition( XY_Z_PARTITION, &
        preserved_coor( XY_Z_COMM_ZX_Y ), xgrades, xinf, xsup, &
        altered_coor( XY_Z_COMM_ZX_Y ), ygrades, yinf, ysup, &
        CONTIGUOUS_INFIMUM, fftdim(Z_DIRECTION) - 1,fftdim(Z_DIRECTION) )


  allocate (transposed_data(nfft1*nfft2*nfft3),stat=fft_alloc_ier)
  REQUIRE(fft_alloc_ier == 0)
  allocate (recv_buffer(nfft1*nfft2*nfft3), &
            send_buffer(nfft1*nfft2*nfft3), &
            stat=fft_alloc_ier)
  REQUIRE(fft_alloc_ier == 0)
  allocate ( &
        fft_table_1(fftdim(X_DIRECTION)*4+15), &
        fft_table_2(fftdim(Y_DIRECTION)*4+15), &
        fft_table_3(fftdim(Z_DIRECTION)*4+15), &
        stat=fft_alloc_ier)
  REQUIRE(fft_alloc_ier == 0)
  
  call cffti(fftdim(X_DIRECTION),fft_table_1)
  call cffti(fftdim(Y_DIRECTION),fft_table_2)
  call cffti(fftdim(Z_DIRECTION),fft_table_3)


  allocate (alpha_rcfft(0:nfft1),beta_rcfft(0:nfft1), &
        stat=fft_alloc_ier)
  REQUIRE(fft_alloc_ier == 0)
  
  allocate (tmp_rc(1:2, 0:fftdim(X_DIRECTION)-1), stat=fft_alloc_ier)
  REQUIRE(fft_alloc_ier == 0)
  
  pi2n=TWOPI/nfft1
  do i=0,fftdim(X_DIRECTION)-1
     theta=pi2n*i
     alpha_rcfft(i) = cos(theta)
     beta_rcfft(i)  = sin(theta)
  end do
  
  call trace_exit( 'fft_init' )
  
  
End subroutine fft_init


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!            FORWARD_RC_FFT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine forward_rc_fft( data, nfft1, nfft2, nfft3)
  !
  ! 3D Forward FFT with Real to Complex transform in x direction.
  !
  !------------------------------------------------------------------
  use trace
  implicit none
#  include "def_time.h"
  
  _DATA_ ,intent(inout) ::data(0:*)  ! original data
  integer,intent(in) :: nfft1,nfft2,nfft3
  
  call z_fft( TRANSFORM_FORWARD, data )
  
  call timer_start(TIME_FFTCOMM3)
  call ftranspose( data, XY_Z_PARTITION, transposed_data, &
        ZX_Y_PARTITION, XY_Z_COMM_ZX_Y, recv_buffer, send_buffer )
  call timer_stop(TIME_FFTCOMM3)
  
  call y_fft( TRANSFORM_FORWARD, transposed_data )
  
  call timer_start(TIME_FFTCOMM3)
  call ftranspose( transposed_data, ZX_Y_PARTITION, data, &
        YZ_X_COMM_ZX_Y, YZ_X_COMM_ZX_Y, recv_buffer, send_buffer )
  call timer_stop(TIME_FFTCOMM3)
 
  call x_fft( TRANSFORM_FORWARD, data )

  return
End subroutine forward_rc_fft


! private routines in alphabetical order


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine create_comm_groups( inner_grades, outer_grades )
  !
  ! Create the communication groups.
  ! Processors are mapped sequentially in one direction (denoted 
  ! the inner direction because of the algorithm in create_partition).
  ! In the other direction the processors are mapped with a stride
  ! of the inner direction's number of interval divisions.
  ! This routine initializes some of the module data structures.
  !
  !------------------------------------------------------------------
  Implicit none
#include "parallel.h"
  
  integer inner_grades  ! number of interval divisions, intent(in)
  integer outer_grades  ! number of interval divisions, intent(in)
  
  integer group_inf
  integer group_sup
  integer i
  integer j
  integer processor
  
  comm_group_stride(YZ_X_COMM_ZX_Y) = 1
  processor = 0
  do i = 1, outer_grades
     group_inf = processor
     group_sup = group_inf + comm_group_stride(YZ_X_COMM_ZX_Y) * &
           ( inner_grades - 1 )
     do j = 1, inner_grades
        comm_group_begin(processor, YZ_X_COMM_ZX_Y) = group_inf
        comm_group_end(processor, YZ_X_COMM_ZX_Y) = group_sup
        processor = processor + comm_group_stride(YZ_X_COMM_ZX_Y)
     enddo
  enddo
  ASSERT( processor == numtasks ) ! one stride past the end
  
  comm_group_stride(XY_Z_COMM_ZX_Y) = inner_grades
  do j = 1, inner_grades
     processor = j - 1
     group_inf = processor
     group_sup = group_inf + comm_group_stride(XY_Z_COMM_ZX_Y) * &
           ( outer_grades - 1 )
     do i = 1, outer_grades
        comm_group_begin(processor, XY_Z_COMM_ZX_Y) = group_inf
        comm_group_end(processor, XY_Z_COMM_ZX_Y) = group_sup
        processor = processor + comm_group_stride(XY_Z_COMM_ZX_Y)
     enddo
  enddo
  ASSERT( processor == numtasks - 1 + inner_grades )
  ! one stride past the end
  
End subroutine create_comm_groups
      

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine create_partition( partition, &
      inner_coor, inner_grades, inner_inf, inner_sup, &
      outer_coor, outer_grades, outer_inf, outer_sup, &
      local_inf, local_sup,local_fftdim )
  !
  ! Create a grid partition by distributing the grid points
  ! among the processors.
  ! This routine initializes some of the module data structures.
  !
  !------------------------------------------------------------------
  Implicit none
#include "parallel.h"
  
  integer partition          ! partition index,              intent(in)
  integer inner_coor         ! rectangle coordinate index,   intent(in)
  integer inner_grades       ! number of interval divisions, intent(in)
  integer inner_inf(MAX_PROCESSORS) ! interval infimums,     intent(in)
  integer inner_sup(MAX_PROCESSORS) ! interval supremums,    intent(in)
  integer outer_coor         ! rectangle coordinate index,   intent(in)
  integer outer_grades       ! number of interval divisions, intent(in)
  integer outer_inf(MAX_PROCESSORS) ! interval infimums,     intent(in)
  integer outer_sup(MAX_PROCESSORS) ! interval supremums,    intent(in)
  integer local_inf          ! contiguous interval infimum,  intent(in)
  integer local_sup          ! contiguous interval supremum, intent(in)
  integer,intent(in) :: local_fftdim
  
  integer i
  integer j
  integer processor


  processor = 0
  do i = 1, outer_grades
     do j = 1, inner_grades
        decomposed_inf(inner_coor, processor, partition) = inner_inf(j)
        decomposed_sup(inner_coor, processor, partition) = inner_sup(j)
        decomposed_inf(outer_coor, processor, partition) = outer_inf(i)
        decomposed_sup(outer_coor, processor, partition) = outer_sup(i)
        processor = processor + 1
     enddo
  enddo
  ASSERT( processor == numtasks ) ! one past the end
  
  contiguous_inf(partition) = local_inf
  contiguous_sup(partition) = local_sup
  contiguous_fftdim(partition) = local_fftdim
  
End subroutine create_partition


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine divide_interval( npoints, ngrades, inf, sup )
  !
  ! Apportion grid points in one direction into closed intervals.
  !
  !------------------------------------------------------------------
  use trace
  Implicit none
  
  integer npoints ! number of grid points in 1 direction, intent(in)
  integer ngrades ! number of divisions to create,        intent(in)
  integer inf(MAX_PROCESSORS) ! interval infimums,        intent(out)
  integer sup(MAX_PROCESSORS) ! interval supremums,       intent(out)
  
  integer i
  integer increment
  integer residue
  
  call Trace_enter( 'divide_interval' )
  
  increment = npoints / ngrades
  residue   = mod( npoints, ngrades )
  ASSERT( increment > 0 )  ! each grade must get at least 1 point
  
  ! [ inf(i), sup(i) ] is a closed interval of grid points which contains
  ! all the grid points for grade i.
  ! Grid indexing ranges from 0 to npoints - 1.
  ! Grade indexing ranges from 1 to ngrades.
  ! Space inefficient implementation since sup(i) = inf(i + 1) - 1.
  
  sup(ngrades) = npoints - 1
  do i = ngrades - 1, 1, -1
     sup(i) = sup(i + 1) - increment
     ! the residual grid points are divvied 1 by 1 
     ! starting with the last interval.
     ! these could have just as well started with the first interval.
     if ( residue > 0 ) then
        residue = residue - 1
        sup(i)  = sup(i) - 1
     endif
     inf(i + 1)  = sup(i) + 1
  enddo
  inf(1) = sup(1) - increment + 1
  ASSERT( inf(1) == 0 )  ! first grade starts at first grid point.
  
  call Trace_exit( 'divide_interval' )
End subroutine divide_interval


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine fft_1d_realcomplex( isign, length, data, table )
  !
  ! 1D real to complex FFT.
  ! The first argument specifies forward or backward transform.
  ! The real data is packed into consecutive real and imaginary
  ! slots in the complex array.
  !
  !------------------------------------------------------------------
  use constants
  Implicit none
  
  integer isign      ! forward or backward transform, intent(in)
  integer length     ! length of transform,           intent(in)
  _DATA_  data(0:*)  ! original data,               intent(inout)
  _REAL_, dimension(1:*) :: table
  _DATA_,parameter :: c_one=(one,zero),c_i=(zero,one)
  
  _REAL_ a
  _REAL_ b
  _REAL_ c
  _REAL_ d
  integer i
  
  if ( isign == TRANSFORM_BACKWARD ) then
     
     do i = 0, length - 1
        tmp_rc(1,i)=real(data(i))
        tmp_rc(2,i)=aimag(data(i))
     enddo
     call cfftf(length, tmp_rc(1,0), table)
     do i = 1, length - 1
        a =  Half*(tmp_rc(1,i)+tmp_rc(1,length-i)) ! Real F even
        b =  Half*(tmp_rc(2,i)-tmp_rc(2,length-i)) ! Imag F even
        c =  Half*(tmp_rc(2,i)+tmp_rc(2,length-i)) ! Real F odd
        d = -Half*(tmp_rc(1,i)-tmp_rc(1,length-i)) ! Imag F odd
        data(i) = (a + alpha_rcfft(i)*c + beta_rcfft(i)*d) * c_one &
              + (b + alpha_rcfft(i)*d - beta_rcfft(i)*c)   * c_i
     enddo
     ! DC and Nyquist terms
     data(0)  = c_one * (tmp_rc(1,0)+tmp_rc(2,0))
     data(length)=c_one * (tmp_rc(1,0)-tmp_rc(2,0))
  else if ( isign == TRANSFORM_FORWARD ) then
     
     do i = 1, length - 1
        a =   real(data(i)) +  real(data(length-i)) ! Real F even     
        b =  aimag(data(i)) - aimag(data(length-i)) ! Imag F even     
        c =  aimag(data(i)) + aimag(data(length-i)) ! F odd contrib   
        d =   real(data(i)) -  real(data(length-i)) ! F odd contrib   
        tmp_rc(1,i) = a - alpha_rcfft(i)*c - beta_rcfft(i)*d         
        tmp_rc(2,i) = b + alpha_rcfft(i)*d - beta_rcfft(i)*c         
     enddo
     tmp_rc(1,0) = real(data(0)) + real(data(length))
     tmp_rc(2,0) = real(data(0)) - real(data(length))

     call cfftb(length, tmp_rc(1,0), table)                 
     do i = 0, length - 1
        data(i) =   c_one * tmp_rc(1,i) &
                  + c_i   * tmp_rc(2,i)
     enddo
     
  else 
     ASSERT( .false. )  ! impossible
  endif
  
End subroutine fft_1d_realcomplex


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine fft_1d_cc( isign, length, fft_table, data )
  !
  ! 1D complex to complex FFT; wrapper for library routines.
  ! The first argument specifies forward or backward transform.
  !
  !------------------------------------------------------------------
  Implicit none

  _REAL_, dimension(1:*) :: fft_table
  integer isign      ! forward or backward transform, intent(in)
  _DATA_  data(*)    ! original data,                 intent(inout)
  integer length     ! length of transform,           intent(in)
  
  if ( isign == TRANSFORM_BACKWARD ) then
     call cfftf( length, data, fft_table )
  else if ( isign == TRANSFORM_FORWARD ) then
     call cfftb( length, data, fft_table )
  else 
     ASSERT( .false. )  ! impossible
  endif
  
End subroutine fft_1d_cc


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
integer function location( k_length, k_inf, k, &
                           j_length, j_inf, j, &
                           i_length, i_inf, i )
  !
  ! Return the grid location of the point i,j,k as an index
  ! into a three dimensional array stored in column major order
  ! and composed of the Cartesian product of three closed intervals.
  ! Essentially this routine converts absolute addresses in global
  ! grid space into relative addresses in processor local grid space.
  ! The kth coordinate is the slowest varying coordinate in memory, and
  ! the ith coordinate is the fastest varying.
  !
  !------------------------------------------------------------------
  implicit none
  integer  k_length  ! length  of the closed interval, intent(in)
  integer  k_inf     ! infimum of the closed interval, intent(in)
  integer  k         ! element of the closed interval, intent(in)
  integer  j_length  ! length  of the closed interval, intent(in)
  integer  j_inf     ! infimum of the closed interval, intent(in)
  integer  j         ! element of the closed interval, intent(in)
  integer  i_length  ! length  of the closed interval, intent(in)
  integer  i_inf     ! infimum of the closed interval, intent(in)
  integer  i         ! element of the closed interval, intent(in)

!  ASSERT( k_inf <= k                         )
!  ASSERT(          k <= k_inf + k_length - 1 )
!  ASSERT( j_inf <= j                         )
!  ASSERT(          j <= j_inf + j_length - 1 )
!  ASSERT( i_inf <= i                         )
!  ASSERT(          i <= i_inf + i_length - 1 )
  location = j_length * i_length * ( k - k_inf ) + &
                        i_length * ( j - j_inf ) + &
                                   ( i - i_inf )

End function location





!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     TRANSPOSE ------- for backward 3d FFT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine transpose( data, from_partition, transposed_data, &
      to_partition, comm_group, recv_buffer, send_buffer )
  !
  ! Global transpose of the grid data using MPI nonblocking calls
  ! The buffer areas are segmented in order of the MPI ranks of
  ! the remote processors in the appropriate communication group.
  ! This performs backwards tranposes and will be merged into a 
  ! single transpose subroutine.
  !
  !------------------------------------------------------------------
  use trace
  Implicit none
#include "ew_parallel.h"
   include 'mpif.h'
#include "parallel.h"
  
  _DATA_   data(0:*)          ! original data,              intent(in)
  integer  from_partition     ! partition index,            intent(in)
  _DATA_   transposed_data(0:*) ! transposed data,          intent(out)
  integer  to_partition       ! partition index,            intent(in)
  integer  comm_group         ! communication group index,  intent(in)
  _DATA_   recv_buffer(0:*)   ! scratch space for receives, intent()
  _DATA_   send_buffer(0:*)   ! scratch space for sends,    intent()
  
  
  integer   TRANSPOSE_TAG
  parameter(TRANSPOSE_TAG = 83)  ! closest prime to ascii( T )
  
  ! most of these are abbreviations, but initialization in declaration is not F77
  integer alt_c            
  integer alt_from_inf
  integer alt_from_sup
  integer alt_from_len
  integer alt_r_to_inf
  integer alt_r_to_sup
  integer alt_r_to_len
  integer alt_s_from_inf
  integer alt_s_from_sup
  integer alt_s_from_len
  integer alt_to_inf
  integer alt_to_sup
  integer alt_to_len
  integer begin
  integer cg
  integer con_from_inf
  integer con_from_sup
  integer con_from_len
  integer con_to_inf
  integer con_to_sup
  integer con_to_len
  integer end
  integer i,ito,ifrom
  integer j,jto,jfrom
  integer k,kto,kfrom
  integer p
  integer pre_c
  integer pre_from_inf
  integer pre_from_sup
  integer pre_from_len
  integer proc_offset    ! processor offset into buffers
  integer receiver
  integer recv_counter
  integer sender
  integer send_counter

  integer recv_offset(0:numtasks - 1) !saved recv buf proc offsets
  integer recv_request(numtasks)
  integer recv_status(MPI_STATUS_SIZE,numtasks)
  integer send_request(numtasks)
  integer send_status(MPI_STATUS_SIZE,numtasks)

!  integer recv_offset(0:MAX_PROCESSORS - 1) !saved recv buf proc offsets
!  integer recv_request(MAX_PROCESSORS)
!  integer recv_status(MPI_STATUS_SIZE,MAX_PROCESSORS)
!  integer send_request(MAX_PROCESSORS)
!  integer send_status(MPI_STATUS_SIZE,MAX_PROCESSORS)

  integer size
  integer step
  integer pe_list(numtasks),irecv
  
  integer mpi_data_type
  call mpi_type_contiguous( 1, MPI_DOUBLE_COMPLEX, mpi_data_type, mpi_ierr )
  call mpi_type_commit( mpi_data_type, mpi_ierr )
#define AMBER_MPI_REAL mpi_data_type
  
  ! abbreviations for shorter line lengths.
  alt_c        = altered_coor( comm_group )
  alt_from_inf = decomposed_inf(alt_c, mytaskid, from_partition)
  alt_from_sup = decomposed_sup(alt_c, mytaskid, from_partition)
  alt_from_len = alt_from_sup - alt_from_inf + 1
  alt_to_inf   = decomposed_inf(alt_c, mytaskid, to_partition)
  alt_to_sup   = decomposed_sup(alt_c, mytaskid, to_partition)
  alt_to_len   = alt_to_sup - alt_to_inf + 1
  begin        = comm_group_begin( mytaskid, comm_group )
  cg           = comm_group
  con_from_inf = contiguous_inf( from_partition )
  con_from_sup = contiguous_sup( from_partition )
  con_from_len = con_from_sup - con_from_inf + 1
  con_to_inf   = contiguous_inf( to_partition )
  con_to_sup   = contiguous_sup( to_partition )
  con_to_len   = con_to_sup - con_to_inf + 1
  end          = comm_group_end( mytaskid, comm_group )
  p            = mytaskid
  pre_c        = preserved_coor( comm_group )
  pre_from_inf = decomposed_inf(pre_c, mytaskid, from_partition)
  pre_from_sup = decomposed_sup(pre_c, mytaskid, from_partition)
  pre_from_len = pre_from_sup - pre_from_inf + 1
  step         = comm_group_stride( comm_group )

  !----------- Post receives --------------------------------------
  proc_offset = 0
  recv_counter = 0

  do sender = begin, end, step
  
     if ( sender .ne. mytaskid ) then
        recv_counter = recv_counter + 1
        pe_list(recv_counter)=sender
        recv_offset( sender ) = proc_offset
        size = alt_to_len * pre_from_len *  & ! alt_s_from_len 
              ( decomposed_sup(alt_c, sender, from_partition) - &
              decomposed_inf(alt_c, sender, from_partition) + 1 )
        call Trace_mpi( 'mpi_irecv', size, TRACE_DPREC, sender )
        call MPI_IRECV( recv_buffer( recv_offset( sender ) ), &
             size, &
             AMBER_MPI_REAL, &
             sender, &
             TRANSPOSE_TAG, &
             recip_comm, &
             recv_request( recv_counter ), &
             mpi_ierr )
        ASSERT(mpi_ierr == MPI_SUCCESS)
        proc_offset = proc_offset + size
     endif
  enddo
  
  !----------- Pack data in transposed order and send ----------------
  proc_offset = 0
  send_counter = 0
  do receiver = begin, end, step
     
     if ( receiver .ne. mytaskid ) then
        alt_r_to_inf = decomposed_inf(alt_c, receiver, to_partition)
        alt_r_to_sup = decomposed_sup(alt_c, receiver, to_partition)
        alt_r_to_len = alt_r_to_sup - alt_r_to_inf + 1
        do k = pre_from_inf, pre_from_sup
           kto  = (k - pre_from_inf) * alt_from_len               + proc_offset
           kfrom= (k - pre_from_inf) * alt_from_len * con_from_len

           do j = alt_from_inf, alt_from_sup
              jto   = kto   + (j - alt_from_inf) &
                   -alt_r_to_inf*pre_from_len*alt_from_len
              jfrom = kfrom + (j - alt_from_inf) * con_from_len - con_from_inf

              do i = alt_r_to_inf, alt_r_to_sup
                 ! the contiguous in the from-partition is the same direction
                 ! as the altered in the to-partition.

                 send_buffer(jto+ i*pre_from_len*alt_from_len ) &
                      =        data( jfrom  + i )

                 ! by convention, the altered coordinate in the
                 ! to-partition is the slowest varying.
                 ! by convention, the preserved coordinate in the
                 ! from-partition is the slowest varying.
              enddo
           enddo
        enddo
        send_counter = send_counter + 1
        size = pre_from_len * alt_from_len * alt_r_to_len
        call Trace_mpi( 'mpi_isend', size, TRACE_DPREC, receiver )
        call MPI_ISEND( send_buffer( proc_offset ), size, &
              AMBER_MPI_REAL, receiver, TRANSPOSE_TAG, recip_comm, &
              send_request( send_counter ), mpi_ierr )
        ASSERT(mpi_ierr == MPI_SUCCESS)
        proc_offset = proc_offset + size
     endif
  enddo
  call Trace_integer( 'Send Buffer on Processor ',p )
  
  ! Transpose local data
  call Trace_integer( 'Local Transpose on Processor ',p )
  do k = pre_from_inf, pre_from_sup
     kto   = (k - pre_from_inf)*con_to_len
     kfrom = (k - pre_from_inf)*alt_from_len*con_from_len
     do j = alt_from_inf, alt_from_sup
        jto   = kto + (j - con_to_inf)- alt_to_inf*pre_from_len*con_to_len
        jfrom = kfrom + (j - alt_from_inf) * con_from_len - con_from_inf
        do i = alt_to_inf, alt_to_sup
!           ito   = jto + i*pre_from_len * con_to_len
!           ifrom = jfrom + i

           ! the altered in the from-partition is the same direction
           ! as the contiguous in the to-partition.
           ! the contiguous in the from-partition is the same direction
           ! as the altered in the to-partition.

           transposed_data( jto+i*pre_from_len*con_to_len ) = &
                data( jfrom + i )

           ! by convention, the altered coordinate in the
           ! to-partition is the slowest varying.
           ! by convention, the preserved coordinate in the
           ! from-partition is the slowest varying.
        enddo
     enddo
  enddo
  
   do irecv=1,recv_counter 
      call mpi_waitany(recv_counter,recv_request,i,recv_status,mpi_ierr)
      ASSERT(mpi_ierr == MPI_SUCCESS)
      sender=pe_list(i)
      alt_s_from_inf = decomposed_inf(alt_c, sender, from_partition)
      alt_s_from_sup = decomposed_sup(alt_c, sender, from_partition)
      alt_s_from_len = alt_s_from_sup - alt_s_from_inf + 1

        do i = 0,alt_to_len-1
           ito   = i*pre_from_len*con_to_len
           ifrom = i*pre_from_len*alt_s_from_len + recv_offset(sender)

           do k = 0, pre_from_len-1
              kto   = ito + k*con_to_len
              kfrom = ifrom + k*alt_s_from_len

              do j = alt_s_from_inf, alt_s_from_sup
                 ! the altered in the from-partition is the same direction
                 ! as the contiguous in the to-partition.

                 transposed_data( kto + j - con_to_inf ) = &
                      recv_buffer( kfrom + j - alt_s_from_inf )

                 ! by convention, the altered coordinate in the
                 ! to-partition is the slowest varying.
              enddo
           enddo
        enddo
  enddo
  call MPI_Waitall( send_counter, send_request, send_status, mpi_ierr)
  call mpi_type_free( mpi_data_type, mpi_ierr )
  
End subroutine transpose



      
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     FTRANSPOSE ------- for foreward 3d FFT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine ftranspose( data, from_partition, transposed_data, &
      to_partition, comm_group, recv_buffer, send_buffer )
  !
  ! Global transpose of the grid data using MPI nonblocking calls
  ! The buffer areas are segmented in order of the MPI ranks of
  ! the remote processors in the appropriate communication group.
  ! This performs forward tranposes and will be merged into a 
  ! single transpose subroutine.
  !
  !------------------------------------------------------------------
  use trace
  Implicit none
#include "ew_parallel.h"
   include 'mpif.h'
#include "parallel.h"
#include "def_time.h"
  
  _DATA_   data(0:*)          ! original data,              intent(in)
  integer  from_partition     ! partition index,            intent(in)
  _DATA_   transposed_data(0:*) ! transposed data,          intent(out)
  integer  to_partition       ! partition index,            intent(in)
  integer  comm_group         ! communication group index,  intent(in)
  _DATA_   recv_buffer(0:*)   ! scratch space for receives, intent()
  _DATA_   send_buffer(0:*)   ! scratch space for sends,    intent()
  
  integer   TRANSPOSE_TAG
  parameter(TRANSPOSE_TAG = 83)  ! closest prime to ascii( T )
  
  ! most of these are abbreviations, 
  ! but initialization in declaration is not F77 (MFC what does this mean?)
  integer alt_c
  integer alt_from_inf
  integer alt_from_sup
  integer alt_from_len
  integer alt_r_to_inf
  integer alt_r_to_sup
  integer alt_r_to_len
  integer alt_s_from_inf
  integer alt_s_from_sup
  integer alt_s_from_len
  integer alt_to_inf
  integer alt_to_sup
  integer alt_to_len
  integer begin
  integer cg
  integer con_from_inf
  integer con_from_sup
  integer con_from_len
  integer con_to_inf
  integer con_to_sup
  integer con_to_len
  integer end
  integer i,ifrom,ito
  integer j,jfrom,jto
  integer k,kfrom,kto
  integer p
  integer pre_c
  integer pre_from_inf
  integer pre_from_sup
  integer pre_from_len
  integer proc_offset    ! processor offset into buffers
  integer receiver
  integer recv_counter
  integer recv_offset(0:numtasks - 1) !saved recv buf proc offsets
  integer recv_request(numtasks)
  integer recv_status(MPI_STATUS_SIZE,numtasks)
  integer sender
  integer send_counter
  integer send_request(numtasks)
  integer send_status(MPI_STATUS_SIZE,numtasks)
  integer size
  integer step

  integer pe_list(numtasks),irecv
  integer trans_timer,waitall_timer,waitany_timer,postrcv_timer,send_timer, &
       packsend_timer


  integer mpi_data_type
  
  trans_timer   = TIME_FFTTRANS3
  waitall_timer = TIME_FFTXTRA3
  postrcv_timer = TIME_FFTXTRA5
  waitany_timer = TIME_FFTXTRA4
  send_timer    = TIME_FFTXTRA6
  packsend_timer    = TIME_FFTCOMM4
  call timer_start(waitall_timer)
!  call timer_start(trans_timer)

  call mpi_type_contiguous( 1, MPI_DOUBLE_COMPLEX, mpi_data_type, mpi_ierr )
  call mpi_type_commit( mpi_data_type, mpi_ierr )
#define AMBER_MPI_REAL mpi_data_type
  
  call Trace_integer( 'Processor is ',p )  ! testing

  ! abbreviations for shorter line lengths.
  alt_c        = altered_coor( comm_group )
  alt_from_inf = decomposed_inf(alt_c, mytaskid, from_partition)
  alt_from_sup = decomposed_sup(alt_c, mytaskid, from_partition)
  alt_from_len = alt_from_sup - alt_from_inf + 1
  alt_to_inf   = decomposed_inf(alt_c, mytaskid, to_partition)
  alt_to_sup   = decomposed_sup(alt_c, mytaskid, to_partition)
  alt_to_len   = alt_to_sup - alt_to_inf + 1
  begin        = comm_group_begin( mytaskid, comm_group )
  cg           = comm_group
  con_from_inf = contiguous_inf( from_partition )
  con_from_sup = contiguous_sup( from_partition )
  con_from_len = con_from_sup - con_from_inf + 1
  con_to_inf   = contiguous_inf( to_partition )
  con_to_sup   = contiguous_sup( to_partition )
  con_to_len   = con_to_sup - con_to_inf + 1
  end          = comm_group_end( mytaskid, comm_group )
  p            = mytaskid
  pre_c        = preserved_coor( comm_group )
  pre_from_inf = decomposed_inf(pre_c, mytaskid, from_partition)
  pre_from_sup = decomposed_sup(pre_c, mytaskid, from_partition)
  pre_from_len = pre_from_sup - pre_from_inf + 1
  step         = comm_group_stride( comm_group )
  call timer_stop(waitall_timer)
  
  ! Post receives

  call timer_start(postrcv_timer)
  proc_offset = 0
  recv_counter = 0
  do sender = begin, end, step
     if ( sender .ne. mytaskid ) then
        recv_counter = recv_counter + 1
        pe_list(recv_counter)=sender
        recv_offset( sender ) = proc_offset
        size = alt_to_len * pre_from_len * & ! alt_s_from_len 
              ( decomposed_sup(alt_c, sender, from_partition) - &
              decomposed_inf(alt_c, sender, from_partition) + 1 )
        call Trace_mpi( 'mpi_irecv', size, TRACE_DPREC, sender )
        call MPI_IRECV( recv_buffer( recv_offset( sender ) ), size, &
              AMBER_MPI_REAL, sender, TRANSPOSE_TAG, recip_comm, &
              recv_request( recv_counter ), mpi_ierr )
        proc_offset = proc_offset + size
     endif
  enddo
  call timer_stop(postrcv_timer)
!  call timer_stop(trans_timer)

  ! Pack data in transposed order and send

  proc_offset = 0
  send_counter = 0
  do receiver = begin, end, step
     if ( receiver .ne. mytaskid ) then
        call timer_start(packsend_timer)
        alt_r_to_inf = decomposed_inf(alt_c, receiver, to_partition)
        alt_r_to_sup = decomposed_sup(alt_c, receiver, to_partition)
        alt_r_to_len = alt_r_to_sup - alt_r_to_inf + 1
        do j = alt_from_inf, alt_from_sup
           jto = j - alt_from_inf + proc_offset
           jfrom = (j - alt_from_inf)* pre_from_len * con_from_len
           do k = pre_from_inf, pre_from_sup
              kto = jto + (k - pre_from_inf) * alt_r_to_len * alt_from_len &
                   -alt_r_to_inf*alt_from_len
              kfrom = jfrom + (k - pre_from_inf) *  con_from_len-con_from_inf
              do i = alt_r_to_inf, alt_r_to_sup
!                 ito = kto + i*alt_from_len
!                 ifrom = kfrom + i

                 ! the contiguous in the from-partition is the same direction
                 ! as the altered in the to-partition.

                 send_buffer( kto + i*alt_from_len ) = data(  kfrom + i )

                 ! by convention, the preserved coordinate in the
                 ! to-partition is the slowest varying.
                 ! by convention, the altered coordinate in the
                 ! from-partition is the slowest varying.
              enddo
           enddo
        enddo
        send_counter = send_counter + 1
        size = pre_from_len * alt_from_len * alt_r_to_len
        call timer_stop(packsend_timer)
        call timer_start(send_timer)
        call Trace_mpi( 'mpi_isend', size, TRACE_DPREC, receiver )
        call MPI_ISEND( send_buffer( proc_offset ), size, &
              AMBER_MPI_REAL, receiver, TRANSPOSE_TAG, recip_comm, &
              send_request( send_counter ), mpi_ierr )
        proc_offset = proc_offset + size
        call timer_stop(send_timer)
     endif
  enddo
  call Trace_integer( 'Send Buffer on Processor ',p )
  
  ! Transpose local data
  call timer_start(trans_timer)
  call Trace_integer( 'Local Transpose on Processor ',p )
  do j = alt_from_inf, alt_from_sup
     jto   = j-con_to_inf
     jfrom = (j-alt_from_inf)*pre_from_len*con_from_len
     do k = pre_from_inf, pre_from_sup
        kto   = jto   + (k-pre_from_inf)*alt_to_len*con_to_len &
             -alt_to_inf*con_to_len
        kfrom = jfrom + (k-pre_from_inf)*con_from_len - con_from_inf
        do i = alt_to_inf, alt_to_sup
           ! the altered in the from-partition is the same direction
           ! as the contiguous in the to-partition.
           ! the contiguous in the from-partition is the same direction
           ! as the altered in the to-partition.

           transposed_data( kto   + i*con_to_len ) = data( kfrom + i )

           ! by convention, the preserved coordinate in the
           ! to-partition is the slowest varying.
           ! by convention, the altered coordinate in the
           ! from-partition is the slowest varying.
        enddo
     enddo
  enddo
  
   call Trace_integer( 'Global Transpose on Processor ',p )
   do irecv=1,recv_counter 
      call mpi_waitany(recv_counter,recv_request,i,recv_status,mpi_ierr)
      sender=pe_list(i)
        alt_s_from_inf = decomposed_inf(alt_c, sender, from_partition)
        alt_s_from_sup = decomposed_sup(alt_c, sender, from_partition)
        alt_s_from_len = alt_s_from_sup - alt_s_from_inf + 1

        do k = pre_from_inf, pre_from_sup
           kto   = (k-pre_from_inf)*alt_to_len*con_to_len
           kfrom = recv_offset( sender ) + &
                (k-pre_from_inf)*alt_to_len*alt_s_from_len

           do i = alt_to_inf, alt_to_sup
              ito   = kto   + (i -alt_to_inf)*con_to_len - con_to_inf
              ifrom = kfrom + (i -alt_to_inf)*alt_s_from_len - alt_s_from_inf

              do j = alt_s_from_inf, alt_s_from_sup
                 ! the altered in the from-partition is the same direction
                 ! as the contiguous in the to-partition.

                 transposed_data(ito+j) = recv_buffer(ifrom+j)

                 ! by convention, the preserved coordinate in the
                 ! to-partition is the slowest varying.
              enddo
           enddo
        enddo
  enddo
  call timer_stop(trans_timer)
  call timer_start(waitall_timer)
  call MPI_Waitall( send_counter, send_request, send_status, mpi_ierr)
  call timer_stop(waitall_timer)
  
  call timer_start(waitall_timer)
  call mpi_type_free( mpi_data_type, mpi_ierr )
  call timer_stop(waitall_timer)
  
End subroutine ftranspose
      

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine x_fft( isign, data )
  !
  ! 1D real to complex FFT in the x direction of grid space.
  ! The first argument specifies forward or backward transform.
  !
  !------------------------------------------------------------------
  Implicit none
#include "parallel.h"
  
  integer isign      ! forward or backward transform, intent(in)
  _DATA_  data(0:*)  ! original data,                 intent(in)
  
  integer alt_c
  integer alt_inf
  integer alt_sup
  integer alt_len
  integer con_inf
  integer con_sup
  integer con_len,con_fftdim
  integer i
  integer index
  integer j
  integer k
  integer p
  integer pre_c
  integer pre_inf
  integer pre_sup
  integer pre_len
  
  p = mytaskid
  
  ! 1D RC FFT in x direction
  alt_c   = altered_coor( YZ_X_COMM_ZX_Y )    ! y direction
  pre_c   = preserved_coor( YZ_X_COMM_ZX_Y )  ! z direction
  pre_inf = decomposed_inf( pre_c, p, YZ_X_PARTITION )
  pre_sup = decomposed_sup( pre_c, p, YZ_X_PARTITION )
  pre_len = pre_sup - pre_inf + 1
  alt_inf = decomposed_inf( alt_c, p, YZ_X_PARTITION )
  alt_sup = decomposed_sup( alt_c, p, YZ_X_PARTITION )
  alt_len = alt_sup - alt_inf + 1
  con_inf = contiguous_inf( YZ_X_PARTITION )  ! x direction
  con_sup = contiguous_sup( YZ_X_PARTITION )
  con_len = con_sup - con_inf + 1
  con_fftdim = contiguous_fftdim(YZ_X_PARTITION)

  do k = pre_inf, pre_sup
     do j = alt_inf, alt_sup
        ! by convention, the z direction is the slowest varying.
        index = location( pre_len, pre_inf, k, alt_len, alt_inf, j, &
              con_len, con_inf, 0 )
        call fft_1d_realcomplex( isign, con_fftdim, data( index ), &
              fft_table_1 )
     enddo
  enddo
End subroutine x_fft


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine y_fft( isign, data )
  !
  ! 1D complex to complex FFT in the y direction of grid space.
  ! The first argument specifies forward or backward transform.
  !
  !------------------------------------------------------------------
  Implicit none
#include "parallel.h"
  
  integer isign      ! forward or backward transform, intent(in)
  _DATA_  data(0:*)  ! original data,                 intent(in)
  
  integer alt_c
  integer alt_inf
  integer alt_sup
  integer alt_len
  integer con_inf
  integer con_sup
  integer con_len
  integer i
  integer index
  integer j
  integer k
  integer p
  integer pre_c
  integer pre_inf
  integer pre_sup
  integer pre_len
  
  p = mytaskid
  
  ! 1D CC FFT in y direction
  alt_c   = altered_coor( XY_Z_COMM_ZX_Y )    ! z direction
  pre_c   = preserved_coor( XY_Z_COMM_ZX_Y )  ! x direction
  pre_inf = decomposed_inf( pre_c, p, ZX_Y_PARTITION )
  pre_sup = decomposed_sup( pre_c, p, ZX_Y_PARTITION )
  pre_len = pre_sup - pre_inf + 1
  alt_inf = decomposed_inf( alt_c, p, ZX_Y_PARTITION )
  alt_sup = decomposed_sup( alt_c, p, ZX_Y_PARTITION )
  alt_len = alt_sup - alt_inf + 1
  con_inf = contiguous_inf( ZX_Y_PARTITION )  ! y direction
  con_sup = contiguous_sup( ZX_Y_PARTITION )
  con_len = con_sup - con_inf + 1
  do i = pre_inf, pre_sup
     do k = alt_inf, alt_sup
        ! by convention, the x direction is the slowest varying.
        index = location( pre_len, pre_inf, i, alt_len, alt_inf, k, &
              con_len, con_inf, 0 )
        call fft_1d_cc( isign, con_len, fft_table_2, data( index ) )
     enddo
  enddo
  
End subroutine y_fft


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Subroutine z_fft( isign, data )
  !
  ! 1D complex to complex FFT in the z direction of grid space.
  ! The first argument specifies forward or backward transform.
  !
  !------------------------------------------------------------------
  Implicit none
#include "parallel.h"
  
  integer isign      ! forward or backward transform, intent(in)
  _DATA_  data(0:*)  ! original data,                 intent(in)
  
  integer alt_c
  integer alt_inf
  integer alt_sup
  integer alt_len
  integer con_inf
  integer con_sup
  integer con_len
  integer i
  integer index
  integer j
  integer k
  integer p
  integer pre_c
  integer pre_inf
  integer pre_sup
  integer pre_len
  
  p = mytaskid
  
  pre_c   = preserved_coor( XY_Z_COMM_ZX_Y )  ! x direction
  alt_c   = altered_coor( XY_Z_COMM_ZX_Y )    ! y direction
  pre_inf = decomposed_inf( pre_c, p, XY_Z_PARTITION )
  pre_sup = decomposed_sup( pre_c, p, XY_Z_PARTITION )
  pre_len = pre_sup - pre_inf + 1
  alt_inf = decomposed_inf( alt_c, p, XY_Z_PARTITION )
  alt_sup = decomposed_sup( alt_c, p, XY_Z_PARTITION )
  alt_len = alt_sup - alt_inf + 1
  con_inf = contiguous_inf( XY_Z_PARTITION )  ! z direction
  con_sup = contiguous_sup( XY_Z_PARTITION )
  con_len = con_sup - con_inf + 1
  do j = alt_inf, alt_sup
     do i = pre_inf, pre_sup
        ! by convention, the y direction is the slowest varying.
        index = location( alt_len, alt_inf, j, pre_len, pre_inf, i, &
              con_len, con_inf, 0 )
        call fft_1d_cc( isign, con_len, fft_table_3, data( index ) )
     enddo
     
  enddo
  
End subroutine z_fft


subroutine get_xy_z_partition_limits(xgmin,xgmax,ygmin,ygmax, &
      zgmin,zgmax,mytaskid)
  implicit none
  integer,intent(in) :: mytaskid
  integer,intent(out) :: xgmin,xgmax,ygmin,ygmax,zgmin,zgmax

  integer :: alt_c,pre_c
  
  pre_c   = preserved_coor( XY_Z_COMM_ZX_Y )  ! x direction
  alt_c   = altered_coor( XY_Z_COMM_ZX_Y )    ! y direction
  xgmin = decomposed_inf( pre_c, mytaskid, XY_Z_PARTITION )
  xgmax = decomposed_sup( pre_c, mytaskid, XY_Z_PARTITION )
  ygmin = decomposed_inf( alt_c, mytaskid, XY_Z_PARTITION )
  ygmax = decomposed_sup( alt_c, mytaskid, XY_Z_PARTITION )
  zgmin = contiguous_inf( XY_Z_PARTITION )  ! z direction
  zgmax = contiguous_sup( XY_Z_PARTITION )

  return
end subroutine get_xy_z_partition_limits



#endif /* MPI */

End module FFT

