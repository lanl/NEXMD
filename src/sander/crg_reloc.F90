! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

!===============================================================================
! charge relocation module:
! add, subtract, multiply, divide and assign charges on atoms based on atomic
! distances and angles
!===============================================================================

module crg_reloc

implicit none

#include "parallel.h"
#ifdef MPI
   include 'mpif.h'
#endif

integer, parameter :: CR_MAX_NPTS = 500

integer, save :: ifcr, cropt, crprintcharges
_REAL_, save :: crcut, crcut2, crskin, crcut_skin, crcut_skin2, &
                half_crskin 
character(len=256) :: crin
! cr_charge holds the original charge array
_REAL_, allocatable, dimension(:), save :: cr_charge
_REAL_, allocatable, dimension(:,:), save :: cr_cub
integer, allocatable, dimension(:,:), save :: cr_info
integer, save :: cr_max_cub = 0, cr_max_info = 0

integer, allocatable, dimension(:,:), save :: cr_order
integer, save :: cr_max_order
! true if charge is updated
! size: cr_max_order
logical, allocatable, dimension(:), save :: cr_upcharge

integer, allocatable, dimension(:,:), save :: cr_info_i
integer, save :: cr_max_info_i
integer, allocatable, dimension(:,:), save :: cr_info_j
integer, save :: cr_max_info_j
integer, allocatable, dimension(:), save :: cr_info_ptr
integer, save :: cr_max_info_ptr

integer, save :: cr_max_info3_i, cr_max_info3_j, cr_max_info3_k, &
                 cr_max_info3_ptr
integer, allocatable, dimension(:,:), save :: cr_info3_i
integer, allocatable, dimension(:,:), save :: cr_info3_j
integer, allocatable, dimension(:,:), save :: cr_info3_k
integer, allocatable, dimension(:), save :: cr_info3_ptr

integer, allocatable, dimension(:), save :: cr_dcdr_tbl
integer, allocatable, dimension(:,:), save :: cr_dcdr_i
integer, save :: cr_max_dcdr_i
integer, allocatable, dimension(:), save :: cr_dcdr_j
integer, save :: cr_max_dcdr_j
_REAL_, allocatable, dimension(:,:), save :: cr_dcdr

! holds distance and delr of atom pairs
_REAL_, allocatable, dimension(:,:), save :: cr_geom
integer, save :: cr_max_geom
! holds distance and delr of atom pairs
_REAL_, allocatable, dimension(:,:), save :: cr_geom3
integer, save :: cr_max_geom3

logical, save :: cr_update_pair
_REAL_, allocatable, dimension(:), save :: cr_pair_distance
logical, allocatable, dimension(:), save :: cr_pair_eval
integer, allocatable, dimension(:,:), save :: cr_tranunit

logical, save :: cr_update_tb
_REAL_, allocatable, dimension(:,:), save :: cr_tb_distance
logical, allocatable, dimension(:), save :: cr_tb_eval
integer, allocatable, dimension(:,:), save :: cr_tb_tranunit

_REAL_, allocatable, dimension(:,:), save :: cr_cect

_REAL_, allocatable, dimension(:), save :: cr_dcdr_fac

contains

!===============================================================================
! read 'crset' from the input file and set 'cr_cub' and 'cr_info'
subroutine cr_read_input( natom )
   integer, intent(in) :: natom
   integer :: cur_cubi=0, cur_cubf=0
   character(len=80) :: line
   integer :: inunit=21, ierror=0, tag=0
   ! concatenated r and c array
   logical, dimension(natom) :: found

   ! namelist crset
   ! mode
   ! 0: atom set
   ! 1: cubic spline function set
   !
   ! type
   ! 1: charge addition
   ! 2: charge subtraction
   ! 3: charge multiplication
   ! 4: charge division
   ! 5: charge reassignment
   ! above 5, these operations are delayed until type <= 5 is found.
   ! 6: charge addition
   ! 7: charge subtraction
   ! 8: charge multiplication
   ! 9: charge division
   !
   ! npts: number of data points
   !
   ! at1: reference atom1
   ! at2: reference atom2
   ! at3: atom of which charge is actually affected
   integer :: mode=0, type=0, npts=0, at1=0, at2=0, at3=0, at4=0, pi1=2, pi2=2
   _REAL_ :: pr1=0.0d0, pr2=0.0d0
   integer :: cur_type=0
   _REAL_, dimension(2) :: cect

   ! namelist cr_data
   ! maximum size of r and c is CR_MAX_NPTS
   _REAL_, dimension(CR_MAX_NPTS) :: r, c

   namelist /crset/ mode, type, at1, at2, at3, at4, npts, r, c, pi1, pi2, pr1, &
                    pr2, cect

   cect = 0.0
   cr_max_geom = 0
   cr_max_geom3 = 0

   ! open input file
   open( unit=inunit, file=crin, status='old', iostat=ierror ) 
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to open ', crin
      call mexit(6,1)
   end if

   found = .false.
   do while ( ierror == 0 )
      read(unit=inunit, nml=crset, iostat=ierror)
      if ( ierror /= 0 ) exit
      tag = 1 ! There should be at least one 'crset' section
      if ( mode == 0 ) then
         call cr_check_crset_atom( cur_type, at1, at2, at3, at4, natom )
         if ( at4 /= 0 ) then
            found(at4) = .true.
         else
            found(at3) = .true.
         end if
         call cr_reallocate_info( at1, at2, at3, at4, cur_type, cur_cubi, &
                                  cur_cubf, cect )
         ! reset
         at1 = 0
         at2 = 0
         at3 = 0
         at4 = 0
      else if ( mode == 1 ) then
         call cr_check_crset_cubspl( type, npts, pi1, pi2 )
         call cr_reallocate_cub( npts, r(1:npts), c(1:npts), type, &
            cur_cubi, cur_cubf, pi1, pi2, pr1, pr2 )
         cur_type = type
         ! reset
         mode = 0
         npts = 0
      end if
      ! reset
      pi1 = 2
      pi2 = 2
      pr1 = 0.0d0
      pr2 = 0.0d0
      cect = 0.0d0
      type = 0
   end do

   if ( tag /= 1 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to read namelist crset'
      call mexit(6,1)
   end if

   call cr_prepare_info( natom )
   call cr_prepare_info3( natom )
   call cr_prepare_order( natom )
   call cr_prepare_dcdr( natom, found )

end subroutine cr_read_input
!===============================================================================


!===============================================================================
subroutine cr_prepare_order( natom )
   integer, intent(in) :: natom
   integer :: i, prev_i, ord_i, ierr, at, ati, at3, at4
   logical :: found(natom)

   found = .false.
   do i = 1, cr_max_info
      at4 = cr_info(4,i)
      if ( at4 == 0 ) then
         at3 = cr_info(3,i)
         found(at3) = .true.
      else
         found(at4) = .true.
      end if
   end do

   cr_max_order = count( found )
   allocate( cr_order(2,cr_max_order), stat=ierr ) 
   if ( ierr /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_order'
      call mexit(6,1)
   end if

   cr_info(9,:) = 0
   ord_i = 1
   do at = 1, natom
      if ( .not. found(at) ) cycle

      prev_i = 0
      do i = 1, cr_max_info
         at4 = cr_info(4,i)
         if ( at4 == 0 ) then
            ati = cr_info(3,i)
         else
            ati = at4
         end if
         if ( at .ne. ati ) cycle
         if ( prev_i > 0 ) then
            cr_info(9,prev_i) = i
         else
            cr_order(1,ord_i) = at 
            cr_order(2,ord_i) = i
         end if
         prev_i = i
      end do

      ord_i = ord_i + 1
   end do
end subroutine cr_prepare_order
!===============================================================================


!===============================================================================
subroutine cr_check_input( ips )
   integer, intent(in) :: ips

   if ( ips > 0 ) then
      write(6,'(x,a)') 'ifcr > 0 and ips > 0: This combination is not supported'
      call mexit(6,1)
   end if

end subroutine cr_check_input
!===============================================================================


!===============================================================================
! check if crset is correct when mode = 0
subroutine cr_check_crset_atom( cur_type, at1, at2, at3, at4, natom )
   integer, intent(in) :: cur_type, at1, at2, at3, at4, natom

   ! atom set
   if ( cur_type == 0 ) then
      write(6,'(x,a)') 'CRGRELOC: charge function should be set first'
      call mexit(6,1)
   end if
   if ( at1 <= 0 .or. at1 > natom ) then
      write(6,'(x,a,i8)') 'CRGRELOC: at1 is out of range ', at1
      call mexit(6,1)
   end if
   if ( at2 <= 0 .or. at2 > natom ) then
      write(6,'(x,a,i8)') 'CRGRELOC: at2 is out of range ', at2
      call mexit(6,1)
   end if
   if ( at3 <= 0 .or. at3 > natom ) then
      write(6,'(x,a,i8)') 'CRGRELOC: at3 is out of range ', at3
      call mexit(6,1)
   end if
   if ( at4 < 0 .or. at4 > natom ) then
      write(6,'(x,a,i8)') 'CRGRELOC: at4 is out of range ', at4
      call mexit(6,1)
   end if
   if ( at1 == at2 ) then
      write(6,'(x,a,i8)') 'CRGRELOC: at1 and at2 cannot be equal ', at1
      call mexit(6,1)
   end if
   if ( at4 /= 0 ) then
      if ( at2 == at3 ) then
         write(6,'(x,a,3i8)') &
            'CRGRELOC: at2 and at3 cannot be equal with at4 ', at2, at3, at4
         call mexit(6,1)
      end if
      if ( at1 == at3 ) then
         write(6,'(x,a,3i8)') &
            'CRGRELOC: at1 and at3 cannot be equal with at4 ', at1, at3, at4
         call mexit(6,1)
      end if
   end if
end subroutine cr_check_crset_atom
!===============================================================================


!===============================================================================
! check if crset is correct when mode = 1
subroutine cr_check_crset_cubspl( type, npts, pi1, pi2 )
   integer, intent(in) :: type, npts, pi1, pi2

   ! cubic spline set
   if ( type <= 0 .or. type > 9 ) then
      write(6,'(x,a)') 'CRGRELOC: type should be (1 <= type <= 9)'
      call mexit(6,1)
   end if
   if ( npts < 4 ) then
      write(6,'(x,a)') 'CRGRELOC: npts should be greater than 4'
      call mexit(6,1)
   end if
   if ( npts > CR_MAX_NPTS ) then
      write(6,'(x,a,i8)') 'CRGRELOC: npts exceeds CR_MAX_NPTS ', CR_MAX_NPTS
      call mexit(6,1)
   end if
   if ( pi1 < 1 .or. pi1 > 2 ) then
      write(6,'(x,a,i8)') 'CRGRELOC: pi1 should be 1 or 2 ', pi1
      call mexit(6,1)
   end if
   if ( pi2 < 1 .or. pi2 > 2 ) then
      write(6,'(x,a,i8)') 'CRGRELOC: pi2 should be 1 or 2 ', pi2
      call mexit(6,1)
   end if
end subroutine cr_check_crset_cubspl
!===============================================================================


!===============================================================================
! increase the size of cr_info
subroutine cr_reallocate_info( at1, at2, at3, at4, cur_type, cur_cubi, &
                               cur_cubf, cect )

   use constants, only : INV_AMBER_ELECTROSTATIC  

   integer, intent(in) :: at1, at2, at3, at4, cur_type, cur_cubi, cur_cubf
   _REAL_, dimension(2), intent(inout) :: cect
   integer, allocatable, dimension(:,:) :: info2
   _REAL_, allocatable, dimension(:,:) :: cect2
   integer :: i, ierror


   cect(1:2) = cect(1:2) * INV_AMBER_ELECTROSTATIC
   cect(2) = cect(2) * INV_AMBER_ELECTROSTATIC

   ! allocate_for info2
   cr_max_info = cr_max_info + 1
   allocate( info2(8,cr_max_info), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate info2'
      call mexit(6,1)
   end if
   allocate( cect2(cr_max_info,2), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cect2'
      call mexit(6,1)
   end if

   ! copy cr_info to info2
   if ( cr_max_info > 1 ) then
      do i = 1, cr_max_info-1
         info2(1:8,i) = cr_info(1:8,i)
         cect2(i,1:2) = cr_cect(i,1:2)
      end do
   end if

   ! add additional item to info2
   if ( at4 == 0 .and. at2 < at1 ) then
      info2(1,cr_max_info) = at2
      info2(2,cr_max_info) = at1
   else
      info2(1,cr_max_info) = at1
      info2(2,cr_max_info) = at2
   end if
   info2(3,cr_max_info) = at3
   info2(4,cr_max_info) = at4
   info2(5,cr_max_info) = cur_type
   info2(6,cr_max_info) = cur_cubi
   info2(7,cr_max_info) = cur_cubf
   if ( at4 == 0 ) then
      cr_max_geom = cr_max_geom + 1
      info2(8,cr_max_info) = cr_max_geom
   else
      cr_max_geom3 = cr_max_geom3 + 1
      info2(8,cr_max_info) = cr_max_geom3
   end if
   cect2(cr_max_info,1:2) = cect

   ! copy info2 back to cr_info
   if ( allocated(cr_info) ) deallocate( cr_info, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_info'
      call mexit(6,1)
   end if
   if ( allocated(cr_cect) ) deallocate( cr_cect, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_cect'
      call mexit(6,1)
   end if
   allocate( cr_info(9,cr_max_info), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info' 
      call mexit(6,1)
   end if
   cr_info(1:8,:) = info2(1:8,:)
   allocate( cr_cect(cr_max_info,2), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_cect' 
      call mexit(6,1)
   end if
   cr_cect = cect2
   deallocate( info2, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate info2' 
      call mexit(6,1)
   end if
   deallocate( cect2, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cect2' 
      call mexit(6,1)
   end if
end subroutine cr_reallocate_info
!===============================================================================


!===============================================================================
! increase the size of cr_cub
subroutine cr_reallocate_cub( npts, r, c, type, cur_cubi, cur_cubf, pi1, pi2, &
                              pr1, pr2 )
   integer, intent(in) :: npts, type
   _REAL_, dimension(npts), intent(in) :: r
   _REAL_, dimension(npts), intent(inout) :: c ! but c does not return values
   integer, intent(in) :: pi1, pi2
   _REAL_, intent(in) :: pr1, pr2 
   _REAL_ :: prev_r
   integer, intent(out) :: cur_cubi, cur_cubf
   integer :: i, ierror
   _REAL_, allocatable, dimension(:,:) :: cub2

   ! allocate for cub2
   allocate( cub2( 5, cr_max_cub+npts ), stat=ierror) 
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cub2'
      call mexit(6,1)
   end if

   ! copy cr_cub to cub2
   if ( cr_max_cub > 0 ) then
      cub2(:,1:cr_max_cub) = cr_cub(:,1:cr_max_cub)
   end if

   ! add additional cub to cub2
   prev_r = r(1)
   do i = 1, npts
      if ( prev_r > r(i) ) then
         write(6,'(x,a)') 'CRGRELOC: r should be in increasing order'
         write(6,'(x,a,f10.3,a,f10.3)' ) 'CRGRELOC: ', prev_r, ' > ', r(i)
      end if
      cub2(1,i+cr_max_cub) = r(i)
      cub2(2,i+cr_max_cub) = c(i)
      prev_r = r(i)
   end do

   ! copy cub2 back to cr_cub
   cur_cubi = cr_max_cub + 1
   cr_max_cub = cr_max_cub + npts
   cur_cubf = cr_max_cub
   if ( allocated(cr_cub) ) deallocate( cr_cub, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_cub'
      call mexit(6,1)
   end if
   allocate( cr_cub(5,cr_max_cub), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_cub' 
      call mexit(6,1)
   end if
   cr_cub = cub2
   deallocate( cub2, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cub2' 
      call mexit(6,1)
   end if

   call cr_calc_cub( cur_cubi, cur_cubf, pi1, pi2, pr1, pr2 )

end subroutine cr_reallocate_cub
!===============================================================================


!===============================================================================
! prepare search table for cr_info for three-body
subroutine cr_prepare_info3( natom )

   integer, intent(in) :: natom
   logical, dimension(natom) :: found
   integer :: ierror, i, j, k, l, si, sj, sk, sp

   ! allocate cr_info3_i
   found = .false.
   cr_max_info3_ptr = 0
   do i = 1, cr_max_info
      if ( cr_info(4,i) == 0 ) cycle
      cr_max_info3_ptr = cr_max_info3_ptr + 1 
      found( cr_info(1,i) ) = .true.
   end do
   cr_max_info3_i = count( found )
   allocate( cr_info3_i(cr_max_info3_i,3), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info3_i' 
      call mexit(6,1)
   end if
   cr_info3_i = 0
   j = 1
   do i = 1, natom
      if ( found(i) .eqv. .false. ) cycle
      cr_info3_i(j,1) = i
      j = j + 1
   end do
   ! allocate cr_info3_ptr
   allocate( cr_info3_ptr(cr_max_info3_ptr), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info3_ptr' 
      call mexit(6,1)
   end if

   ! allocate_cr_info3_j
   ! count the size of cr_max_info3_j
   cr_max_info3_j = 0
   do i = 1, cr_max_info3_i
      si = cr_info3_i(i,1)
      found = .false.
      do j = 1, cr_max_info
         if ( cr_info(4,j) == 0 ) cycle
         if ( cr_info(1,j) /= si ) cycle
         found( cr_info(2,j) ) = .true.
      end do
      cr_max_info3_j = cr_max_info3_j + count( found )
   end do
   allocate( cr_info3_j( cr_max_info3_j, 3 ), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info3_j' 
      call mexit(6,1)
   end if
   cr_info3_j = 0
   ! fill cr_info3_j
   sj = 1
   do i = 1, cr_max_info3_i
      si = cr_info3_i(i,1)
      found = .false.
      do j = 1, cr_max_info
         if ( cr_info(4,j) == 0 ) cycle
         if ( cr_info(1,j) /= si ) cycle
         found( cr_info(2,j) ) = .true.
      end do
      do j = 1, natom
         if ( found(j) .eqv. .false. ) cycle
         cr_info3_j( sj, 1 ) = j
         if ( cr_info3_i( i, 2 ) == 0 ) then
            cr_info3_i( i, 2 ) = sj
         end if
         cr_info3_i( i, 3 ) = sj
         sj = sj + 1
      end do
   end do

   ! allocate cr_info3_k
   ! count the size of cr_max_info3_k
   cr_max_info3_k = 0
   do i = 1, cr_max_info3_i
      si = cr_info3_i(i,1)
      do j = cr_info3_i(i,2), cr_info3_i(i,3)
         sj = cr_info3_j(j,1)
         found = .false.
         do k = 1, cr_max_info
            if ( cr_info(4,k) == 0 ) cycle
            if ( cr_info(1,k) /= si .or. cr_info(2,k) /= sj ) cycle
            found( cr_info(3,k) ) = .true.
         end do
         cr_max_info3_k = cr_max_info3_k + count( found )
      end do
   end do
   allocate( cr_info3_k( cr_max_info3_k, 3 ), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info3_k' 
      call mexit(6,1)
   end if
   cr_info3_k = 0
   ! fill cr_info3_k
   sk = 1 
   do i = 1, cr_max_info3_i
      si = cr_info3_i(i,1)
      do j = cr_info3_i(i,2), cr_info3_i(i,3)
         sj = cr_info3_j(j,1)
         found = .false.
         do k = 1, cr_max_info
            if ( cr_info(4,k) == 0 ) cycle
            if ( cr_info(1,k) /= si .or. cr_info(2,k) /= sj ) cycle
            found( cr_info(3,k) ) = .true.
         end do
         do k = 1, natom
            if ( found(k) .eqv. .false. ) cycle
            cr_info3_k(sk,1) = k
            if ( cr_info3_j( j, 2 ) == 0 ) then
               cr_info3_j( j, 2 ) = sk
            end if
            cr_info3_j( j, 3 ) = sk
            sk = sk + 1
         end do
      end do
   end do

   ! fill cr_info3_ptr
   sp = 1
   do i = 1, cr_max_info3_i
      si = cr_info3_i(i,1)
      do j = cr_info3_i(i,2), cr_info3_i(i,3)
         sj = cr_info3_j(j,1)
         do k = cr_info3_j(j,2), cr_info3_j(j,3)
            sk = cr_info3_k(k,1)
            do l = 1, cr_max_info
               if ( cr_info(4,l) == 0 ) cycle
               if (      cr_info(1,l) /= si &
                    .or. cr_info(2,l) /= sj &
                    .or. cr_info(3,l) /= sk &
                  ) cycle
               if ( cr_info3_k(k,2) == 0 ) then
                  cr_info3_k( k, 2 ) = sp
               end if
               cr_info3_k( k, 3 ) = sp
               cr_info3_ptr(sp) = l
               sp = sp + 1
            end do
         end do
      end do
   end do
            
end subroutine cr_prepare_info3
!===============================================================================


!===============================================================================
! prepare search table for cr_info
subroutine cr_prepare_info( natom )
   ! i < j
   integer, intent(in) :: natom
   logical, dimension(natom) :: found
   integer :: ierror, i, j, k, l, m, si, sj, sp
   integer, allocatable, dimension(:,:) :: j2

   ! allocate cr_info_i
   found = .false.
   cr_max_info_ptr = 0
   do i = 1, cr_max_info
      if ( cr_info(4,i) /= 0 ) cycle ! three-body
      cr_max_info_ptr = cr_max_info_ptr + 1
      found( cr_info(1,i) ) = .true.
   end do
   cr_max_info_i = count( found )
   allocate( cr_info_i(cr_max_info_i,3), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info_i' 
      call mexit(6,1)
   end if
   j = 1
   do i = 1, natom
      if ( found(i) .eqv. .false. ) cycle
      cr_info_i(j,1) = i
      j = j + 1
   end do

   ! allocate cr_info_j and cr_info_ptr
   allocate( cr_info_ptr(cr_max_info_ptr), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info_ptr' 
      call mexit(6,1)
   end if
   si = 1
   sp = 1
   cr_max_info_j = 0
   do i = 1, cr_max_info_i
      j = cr_info_i(i,1)
      sj = 0
      found = .false.
      do k = 1, cr_max_info
         if ( j /= cr_info(1,k) .or. cr_info(4,k) /= 0 ) cycle
         found( cr_info(2,k) ) = .true.
      end do
      sj = count( found )
      allocate( j2(cr_max_info_j + sj,3), stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to deallocate j2' 
         call mexit(6,1)
      end if
      ! copy cr_info_j to j2
      if ( cr_max_info_j > 0 ) then
         j2(1:cr_max_info_j,1:3) = cr_info_j(1:cr_max_info_j,1:3)
      end if
      k = cr_max_info_j + 1
      ! add additional info
      do l = 1, natom
         if ( found(l) .eqv. .false. ) cycle
         j2(k,1) = l
         j2(k,2) = 0
         do m = 1, cr_max_info
            if ( cr_info(4,m) /= 0 ) cycle
            if ( j /= cr_info(1,m) .or. l /= cr_info(2,m) ) cycle
            cr_info_ptr(sp) = m
            if ( j2(k,2) == 0 ) j2(k,2) = sp
            j2(k,3) = sp
            sp = sp + 1
         end do
         k = k + 1
      end do 
      cr_info_i(si,2) = cr_max_info_j + 1
      cr_max_info_j = cr_max_info_j + sj
      cr_info_i(si,3) = cr_max_info_j
      si = si + 1
      if ( allocated(cr_info_j) ) then
         deallocate( cr_info_j, stat=ierror )
         if ( ierror /= 0 ) then
            write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_info_j' 
            call mexit(6,1)
         end if
      end if
      allocate( cr_info_j(cr_max_info_j,3), stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info_j'
         call mexit(6,1)
      end if
      ! copy j2 back to cr_info_j
      cr_info_j = j2
      deallocate( j2, stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to deallocate j2'
         call mexit(6,1)
      end if
   end do

end subroutine cr_prepare_info
!===============================================================================


!===============================================================================
! calculate coefficients of cubic spline
subroutine cr_calc_cub( i, f, pi1, pi2, pr1, pr2 )
   integer, intent(in) :: i, f, pi1, pi2
   _REAL_, intent(in) :: pr1, pr2

   ! interface for cubspl subroutine (from ew_setup.f)
   interface
      subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
         integer ibcbeg,ibcend,n
         _REAL_ c(4,n),tau(n)
      end subroutine
   end interface

   ! calculate cub
   cr_cub( 3, i ) = pr1  ! either slope or 2nd-deriv at the beginning point
   cr_cub( 3, f ) = pr2  ! eigher slope or 2nd-deriv at the end point
   call cubspl( cr_cub(1,i:f), cr_cub(2:5,i:f), f-i+1, pi1, pi2 )
end subroutine cr_calc_cub
!===============================================================================


!===============================================================================
! print information of input to sander output file 
subroutine cr_print_info( outunit )
   integer, intent(in) :: outunit
   integer :: i, n

   write (outunit, '(/a)') 'Charge relocation setup:'
   do i=1,cr_max_info
      n = cr_info(7,i) - cr_info(6,i) + 1
      write (outunit,'(5x,3(a,i8))') 'at1     =', cr_info(1,i), &
         ', at2     =', cr_info(2,i), ', at3     = ', cr_info(3,i)
      write (outunit,'(5x,2(a,i8))') &
         'type    =', cr_info(5,i), ', npts    =', n
   end do
end subroutine cr_print_info
!===============================================================================


!===============================================================================
! Some arrays are only allocated for the master in cr_read_input.
! Here, they are allocated for the rest of the processors and data are copied
! from the master.
subroutine cr_allocate( master, natom )
   logical, intent(in) :: master
   integer, intent(in) :: natom
   integer :: ierror

! JMS: this subroutine should only be called when ifcr /= 0 (sander.f)
!  if ( ifcr == 0 ) return

#ifdef MPI
   call mpi_bcast(cropt,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(crcut,1,MPI_DOUBLE_PRECISION,0,commsander,ierror)
   call mpi_bcast(crskin,1,MPI_DOUBLE_PRECISION,0,commsander,ierror)
   call mpi_bcast(cr_max_cub,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_info,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_info_i,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_info_j,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_info_ptr,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_info3_i,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_info3_j,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_info3_k,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_info3_ptr,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_dcdr_i,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_dcdr_j,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_geom,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_geom3,1,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_max_order,1,MPI_INTEGER,0,commsander,ierror)
#endif

   allocate( cr_pair_distance(cr_max_info_j), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_pair_distance'
      call mexit(6,1)
   end if
   allocate( cr_pair_eval(cr_max_info_j), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_pair_eval' 
      call mexit(6,1)
   end if
   allocate( cr_tranunit(3,cr_max_info_j), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_tranunit'
      call mexit(6,1)
   end if

   allocate( cr_tb_distance(3,cr_max_info3_k), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_tb_distance'
      call mexit(6,1)
   end if
   allocate( cr_tb_eval(cr_max_info3_k), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_tb_eval' 
      call mexit(6,1)
   end if
   allocate( cr_tb_tranunit(9,cr_max_info3_k), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_tb_tranunit'
      call mexit(6,1)
   end if
   allocate( cr_dcdr_fac(cr_max_dcdr_i), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_max_dcdr_i'
      call mexit(6,1)
   end if
   allocate( cr_upcharge(cr_max_order), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_upcharge' 
      call mexit(6,1)
   end if

   if ( .not. master ) then
      allocate(cr_charge(natom), stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_charge' 
         call mexit(6,1)
      end if
      allocate(cr_cub(5,cr_max_cub), stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_cub' 
         call mexit(6,1)
      end if
      allocate(cr_info(9,cr_max_info), stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info' 
         call mexit(6,1)
      end if
      allocate(cr_info_i(cr_max_info_i,3), stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info_i' 
         call mexit(6,1)
      end if
      allocate(cr_info_j(cr_max_info_j,3), stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info_j' 
         call mexit(6,1)
      end if
      allocate(cr_info_ptr(cr_max_info_ptr), stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info_ptr' 
         call mexit(6,1)
      end if
      allocate(cr_info3_i(cr_max_info3_i,3), stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info3_i' 
         call mexit(6,1)
      end if
      allocate(cr_info3_j(cr_max_info3_j,3), stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info3_j' 
         call mexit(6,1)
      end if
      allocate(cr_info3_k(cr_max_info3_k,3), stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info3_k' 
         call mexit(6,1)
      end if
      allocate(cr_info3_ptr(cr_max_info3_ptr), stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_info3_ptr' 
         call mexit(6,1)
      end if
      allocate( cr_dcdr_tbl(natom), stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_dcdr_tbl' 
         call mexit(6,1)
      end if
      allocate( cr_dcdr_i(cr_max_dcdr_i,2), stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_dcdr_i' 
         call mexit(6,1)
      end if
      allocate( cr_dcdr_j(cr_max_dcdr_j), stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_dcdr_j' 
         call mexit(6,1)
      end if
      allocate( cr_dcdr(cr_max_dcdr_j,3), stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_dcdr' 
         call mexit(6,1)
      end if
      allocate( cr_order(2,cr_max_order), stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_order' 
         call mexit(6,1)
      end if
      allocate( cr_cect(cr_max_info,2), stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_cect' 
         call mexit(6,1)
      end if
   end if

#ifdef MPI
   call mpi_bcast(cr_charge,natom,MPI_DOUBLE_PRECISION,0,commsander,ierror)
   call mpi_bcast(cr_cub,5*cr_max_cub,MPI_DOUBLE_PRECISION,0,commsander,ierror)
   call mpi_bcast(cr_info,9*cr_max_info,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_info_i,3*cr_max_info_i,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_info_j,3*cr_max_info_j,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_info_ptr,cr_max_info_ptr,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_info3_i,3*cr_max_info3_i,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_info3_j,3*cr_max_info3_j,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_info3_k,3*cr_max_info3_k,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_info3_ptr,cr_max_info3_ptr,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_dcdr_tbl,natom,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_dcdr_i,cr_max_dcdr_i*2,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_dcdr_j,cr_max_dcdr_j,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_order,2*cr_max_order,MPI_INTEGER,0,commsander,ierror)
   call mpi_bcast(cr_cect,cr_max_info*2,MPI_DOUBLE_PRECISION,0,commsander,ierror)
#endif

   ! precalculation & initilization
   cr_update_pair = .true.
   cr_update_tb = .true.
   crcut2 = crcut * crcut
   crcut_skin = crcut + crskin
   crcut_skin2 = crcut_skin * crcut_skin
   half_crskin = crskin * 0.5

end subroutine cr_allocate
!===============================================================================


!===============================================================================
subroutine cr_prepare_dcdr( natom, foundi )
   integer, intent(in) :: natom
   logical, intent(in) :: foundi(natom)
   integer :: ierror
   integer :: i, j, k, m, si, sj, at
   integer :: at1n, at2n, at3n
   logical, dimension(natom) :: foundj
   integer, allocatable, dimension(:) :: jlist2

   ! cr_dcdr_tbl
   allocate( cr_dcdr_tbl(natom), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_dcdr_tbl' 
      call mexit(6,1)
   end if
   ! initialize
   cr_dcdr_tbl = 0

   ! cr_dcdr_i
   cr_max_dcdr_i = count( foundi )
   allocate( cr_dcdr_i(cr_max_dcdr_i,2), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_dcdr_i' 
      call mexit(6,1)
   end if
   j = 1
   do i=1,natom
      if ( foundi(i) .eqv. .false. ) cycle
      cr_dcdr_tbl(i) = j
      j = j + 1
   end do

   ! cr_dcdr_j
   si = 1
   cr_max_dcdr_j = 0
   do at=1,natom
      k = cr_dcdr_tbl(at)
      if ( k == 0 ) cycle
      sj = 0
      foundj = .false.
      do m=1,cr_max_info
         if ( cr_info(4,m) == 0 ) then
            if ( cr_info(3,m) /= at ) cycle
            at1n = cr_info(1,m)
            at2n = cr_info(2,m)
            foundj(at1n) = .true.
            foundj(at2n) = .true.
         else
            if ( cr_info(4,m) /= at ) cycle
            at1n = cr_info(1,m)
            at2n = cr_info(2,m)
            at3n = cr_info(3,m)
            foundj(at1n) = .true.
            foundj(at2n) = .true.
            foundj(at3n) = .true.
         end if
      end do
      sj = count( foundj )
      allocate( jlist2( cr_max_dcdr_j+sj ), stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate jlist2'
         call mexit(6,1)
      end if
      ! copy cr_dcdr_j to jlist2
      if ( cr_max_dcdr_j > 0 ) then
         jlist2(1:cr_max_dcdr_j) = cr_dcdr_j(1:cr_max_dcdr_j)
      end if
      j = cr_max_dcdr_j+1
      ! add additional info
      do i = 1, natom
         if ( foundj(i) .eqv. .false. ) cycle
         jlist2(j) = i
         j = j + 1
      end do
      cr_dcdr_i(si,1) = cr_max_dcdr_j + 1
      cr_max_dcdr_j = cr_max_dcdr_j + sj
      cr_dcdr_i(si,2) = cr_max_dcdr_j
      si = si + 1
      if ( allocated(cr_dcdr_j) ) deallocate(cr_dcdr_j, stat=ierror)
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_dcdr_j'
         call mexit(6,1)
      end if
      allocate( cr_dcdr_j(cr_max_dcdr_j), stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_dcdr_j'
         call mexit(6,1)
      end if
      ! copy jlist2 back to cr_dcdr_j
      cr_dcdr_j = jlist2
      deallocate( jlist2, stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to deallocate jlist2'
         call mexit(6,1)
      end if
   end do

   ! allocate cr_dcdr for later use
   allocate( cr_dcdr(cr_max_dcdr_j,3), stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_dcdr' 
      call mexit(6,1)
   end if

! following lines are for debugging only
end subroutine cr_prepare_dcdr
!===============================================================================


!===============================================================================
! cleanup routine
subroutine cr_cleanup()
   integer :: ierror=0

   deallocate( cr_charge, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_charge'
      call mexit(6,1)
   end if
   deallocate( cr_cub, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_cub' 
      call mexit(6,1)
   end if
   deallocate( cr_info, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_info'
      call mexit(6,1)
   end if
   deallocate( cr_info_i, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_info_i'
      call mexit(6,1)
   end if
   if ( allocated(cr_info_j) ) then
      deallocate( cr_info_j, stat=ierror )
      if ( ierror /= 0 ) then
         write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_info_j'
         call mexit(6,1)
      end if
   end if
   deallocate( cr_info_ptr, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_info_ptr'
      call mexit(6,1)
   end if
   deallocate( cr_info3_i, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_info3_i'
      call mexit(6,1)
   end if
   deallocate( cr_info3_j, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_info3_j'
      call mexit(6,1)
   end if
   deallocate( cr_info3_k, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_info3_k'
      call mexit(6,1)
   end if
   deallocate( cr_info3_ptr, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_info3_ptr'
      call mexit(6,1)
   end if
   deallocate( cr_dcdr_tbl, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_dcdr_tbl'
      call mexit(6,1)
   end if
   deallocate( cr_dcdr_i, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_dcdr_i'
      call mexit(6,1)
   end if
   deallocate( cr_dcdr_j, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_dcdr_j'
      call mexit(6,1)
   end if
   deallocate( cr_dcdr, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_dcdr'
      call mexit(6,1)
   end if
   deallocate( cr_pair_distance, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_pair_distance'
      call mexit(6,1)
   end if
   deallocate( cr_pair_eval, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_pair_eval' 
      call mexit(6,1)
   end if
   deallocate( cr_tranunit, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_tranunit'
      call mexit(6,1)
   end if
   deallocate( cr_tb_distance, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_tb_distance'
      call mexit(6,1)
   end if
   deallocate( cr_tb_eval, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_tb_eval' 
      call mexit(6,1)
   end if
   deallocate( cr_tb_tranunit, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_tb_tranunit'
      call mexit(6,1)
   end if
   deallocate( cr_order, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_order'
      call mexit(6,1)
   end if
   deallocate( cr_cect, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_cect'
      call mexit(6,1)
   end if
   deallocate( cr_dcdr_fac, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_dcdr_fac'
      call mexit(6,1)
   end if
   deallocate( cr_upcharge, stat=ierror )
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to deallocate cr_upcharge'
      call mexit(6,1)
   end if
end subroutine cr_cleanup
!==============================================================================


!==============================================================================
subroutine cr_acos( cos_angle, angle )
   use constants, only: PI

   _REAL_, intent(in) :: cos_angle
   _REAL_, intent(out) :: angle

   if ( cos_angle > 1.0 ) then
      angle = 0.0
   else if ( cos_angle < -1.0 ) then
      angle = PI
   else
      angle = acos(cos_angle)
   end if
end subroutine cr_acos
!==============================================================================


!==============================================================================
subroutine cr_calc_charge3( geom3, geom3_mask )

   _REAL_, intent(inout) :: geom3(13,cr_max_geom3)
   logical, intent(out) :: geom3_mask(cr_max_geom3)
   integer :: i, j, k, l, gi, a1, a2, a3, ierr
   _REAL_ :: d1, d2, d3, d1inv, d2inv
   _REAL_ :: delr1(3), delr2(3), delr3(3)
   _REAL_ :: angle, cos_angle, sin_angle, sin_inv, dfda, prod
   _REAL_ :: co1, co2, co3
   
#ifdef MPI
   integer :: jobcount, up_taskid

   jobcount = 1
#endif

   do i = 1, cr_max_info
      if ( cr_info(4,i) == 0 ) cycle
      gi = cr_info(8,i)
      d1 = geom3(2,gi)
      d2 = geom3(6,gi)
      d3 = geom3(10,gi)
      if (      d1 == 0.0 .or. d1 >= crcut &
         ) &
      then
         geom3_mask(gi) = .false.
         cycle
      end if
#ifdef MPI
      if ( mod(jobcount,numtasks) .ne. mytaskid ) then
         jobcount = jobcount + 1
         geom3_mask(gi) = .true.
         cycle
      else
         jobcount = jobcount + 1
      end if   
#endif
      geom3_mask(gi) = .true.
      cos_angle = geom3(1,gi)
      j = cr_info(6,i)
      k = cr_info(7,i)

      call cr_acos( cos_angle, angle )

      call cr_cspline( angle, geom3(1,gi), cr_cub(1,j:k), dfda, &
                       cr_cub(2:5,j:k), k-j+1, ierr )

      if ( ierr == 1 ) then
         write(6,'(x,a)')  &                 
            'CRGRELOC: Angle is less than the lower limit'
         write(6,'(x,a,3(a,i8),2(a,f4.2))')  &
            'CRGRELOC: ', 'at1 ', cr_info(1,i), ' at2 ', cr_info(2,i),&
            ' at3 ', cr_info(3,i), ' ANGLE ', angle, ' LIMIT ', &
            cr_cub(1,j)
         call mexit(6,1)
      else if ( ierr == 2 ) then
         write(6,'(x,a)')  &
            'CRGRELOC: Angle is greater than the upper limit'
         write(6,'(x,a,3(a,i8),2(a,f4.2))')  &
            'CRGRELOC: ', 'at1 ', cr_info(1,i), ' at2 ', cr_info(2,i),&
            ' at3 ', cr_info(3,i), ' ANGLE ', angle, ' LIMIT ', &
            cr_cub(1,k)
         call mexit(6,1)
      end if


      sin_angle = sin(angle)
      if ( sin_angle == 0.0 ) then
         sin_inv = 1.0d10  ! big number
      else
         sin_inv = 1.0 / sin_angle
      end if

      prod = -sin_inv * dfda
   
      delr1 = geom3(3:5,gi)
      delr2 = geom3(7:9,gi)
      delr3 = geom3(11:13,gi)

      d1inv = 1.0 / d1
      d2inv = 1.0 / d2

      co1 = prod * ( d2inv - cos_angle * d1inv ) * d1inv
      co2 = prod * ( d1inv - cos_angle * d2inv ) * d2inv
      co3 = prod * (-d1inv) * d2inv

      ! df/dx1, df/dy1, df/dz1
      geom3(2:4,gi) = co1 * (-delr1) + co3 * (-delr3)
      ! df/dx2, df/dy2, df/dz2
      geom3(5:7,gi) = co1 * delr1 + co2 * (-delr2)
      ! df/dx3, df/dy3, df/dz3
      geom3(8:10,gi) = co2 * delr2 + co3 * delr3
   end do

#ifdef MPI
   ! bcast geom3
   jobcount = 1
   do i = 1, cr_max_info
      if ( cr_info(4,i) == 0 ) cycle
      gi = cr_info(8,i)
      if (.not. geom3_mask(gi)) cycle

      up_taskid = mod( jobcount, numtasks )
      jobcount = jobcount + 1
      call mpi_bcast( geom3(1:10,gi), 10, MPI_DOUBLE_PRECISION, &
                      up_taskid, commsander, ierr )   
   end do
#endif

end subroutine cr_calc_charge3
!==============================================================================


!==============================================================================
subroutine cr_calc_charge( geom, geom_mask )

   _REAL_, intent(inout) :: geom(4,cr_max_geom)
   logical, intent(out) :: geom_mask(cr_max_geom)
   integer :: i, j, k, gi, ierr
   _REAL_ :: delr(3), dfdr
   _REAL_ :: distance
   
#ifdef MPI
   integer :: jobcount, up_taskid

   jobcount = 1
#endif

   do i = 1, cr_max_info
      if ( cr_info(4,i) /= 0 ) cycle
      gi = cr_info(8,i)
      distance = geom(1,gi)
      if ( distance == 0.0 .or. distance >= crcut ) then
         geom_mask(gi) = .false.
         cycle
      end if
#ifdef MPI
      if ( mod(jobcount,numtasks) .ne. mytaskid ) then
         jobcount = jobcount + 1
         geom_mask(gi) = .true.
         cycle
      else
         jobcount = jobcount + 1
      end if   
#endif

      geom_mask(gi) = .true.
      j = cr_info(6,i)
      k = cr_info(7,i)

      call cr_cspline( distance, geom(1,gi), cr_cub(1,j:k), dfdr, &
                       cr_cub(2:5,j:k), k-j+1, ierr )

      if ( ierr == 1 ) then
         write(6,'(x,a)')  &                 
            'CRGRELOC: Distance is shorter than the lower limit'
         write(6,'(x,a,3(a,i8),2(a,f4.2))')  &
            'CRGRELOC: ', 'at1 ', cr_info(1,i), ' at2 ', cr_info(2,i),&
            ' at3 ', cr_info(3,i), ' DISTANCE ', distance, ' LIMIT ', &
            cr_cub(1,j)
         call mexit(6,1)
      else if ( ierr == 2 ) then
         write(6,'(x,a)')  &
            'CRGRELOC: Distance is longer than the upper limit'
         write(6,'(x,a,3(a,i8),2(a,f4.2))')  &
            'CRGRELOC: ', 'at1 ', cr_info(1,i), ' at2 ', cr_info(2,i),&
            ' at3 ', cr_info(3,i), ' DISTANCE ', distance, ' LIMIT ', &
            cr_cub(1,k)
         call mexit(6,1)
      end if

      delr = geom(2:4,gi)
      geom(2:4,gi) = dfdr * delr / distance
   end do

#ifdef MPI
   ! bcast geom
   jobcount = 1
   do i = 1, cr_max_info
      if ( cr_info(4,i) /= 0 ) cycle
      gi = cr_info(8,i)
      if ( .not. geom_mask(gi)) cycle

      up_taskid = mod( jobcount, numtasks )
      jobcount = jobcount + 1
      call mpi_bcast( geom(1:4,gi), 4, MPI_DOUBLE_PRECISION, &
                      up_taskid, commsander, ierr )   
   end do
#endif

end subroutine cr_calc_charge
!==============================================================================


!==============================================================================
subroutine cr_cspline( x, y, absc, dydx, cub, n, stat )

   use constants, only : half, third

   integer, intent(in) :: n
   _REAL_, intent(in) :: x
   _REAL_, intent(out) :: y, dydx
   _REAL_, intent(in) :: absc(n), cub(4,n)
   integer, intent(out) :: stat
   _REAL_ :: h
   integer :: i, ind

   stat = 0
   if ( x < absc(1) ) then
      ! x is out-of-range (too small)
      stat = 1
      return
   end if

   ind = 0
   do i = 2, n
      if ( x <= absc(i) ) then
         ind = i - 1
         exit
      end if
   end do

   if ( ind == 0 ) then
      ! x is out-of-range (too big)
      stat = 2
      return
   end if

   h = x - absc(ind)
   y =   cub(1,ind)                              &
       + h * (   cub(2,ind)                      &
               + h * (   cub(3,ind)              &
                       + h * cub(4,ind) * third  &
                     ) * half                    &
             )

   dydx = (   cub(2,ind)                    &
            + h * (   cub(3,ind)            &
                    + h * cub(4,ind) * half &
                  )                         &
          )
end subroutine cr_cspline
!==============================================================================


!===============================================================================
subroutine cr_reassign_charge( crd, force, ect, charge, natom )

   integer, intent(in) :: natom
   _REAL_, intent(out) :: ect
   _REAL_, intent(in) :: crd(3,natom)
   _REAL_, intent(inout) :: force(3,natom)
   _REAL_, intent(out) :: charge(natom)
   logical :: geom_mask(cr_max_geom), geom3_mask(cr_max_geom3)
   _REAL_ :: geom(4,cr_max_geom), geom3(13,cr_max_geom3)

   ! initialize charge transfer energy
   ect = 0.0
   ! initialize dcdr factor
   cr_dcdr_fac = 0.0

   ! cr_fill_geom will fill geom with
   call cr_fill_geom( crd, geom, natom ) 

   ! cr_calc_charge will fill geom with
   call cr_calc_charge( geom, geom_mask )

   ! cr_fill_geom3 will fill geom3 with
   call cr_fill_geom3( crd, geom3, natom )

   ! cr_calc_charge3 will fill geom3 with
   call cr_calc_charge3( geom3, geom3_mask )

   call cr_update_charge( ect, force, &
      charge, natom, geom, geom3, geom_mask, geom3_mask, cr_upcharge )

end subroutine cr_reassign_charge 
!===============================================================================


!===============================================================================
subroutine cr_print_charge( charge, nstep )
   use constants, only : INV_AMBER_ELECTROSTATIC  

   integer, intent(in) :: nstep
   _REAL_, intent(in) :: charge(*) ! natom
   integer :: i, a_mod

   write(6,'(x,a,i8)') 'CRGRELOC: Modified Atomic Charges for Step: ', nstep
   write(6,'(x,a)') ' ATOM  CHARGE'
   do i = 1, cr_max_order
      if ( .not. cr_upcharge(i) ) cycle
      a_mod = cr_order(1,i)
      write(6,'(x,i5,a,f7.4)') a_mod, ' ', charge(a_mod)*INV_AMBER_ELECTROSTATIC
   end do
   write(6,'(x,a)') '-------------'
end subroutine cr_print_charge
!===============================================================================


!===============================================================================
subroutine cr_update_charge( ect, force, charge, natom, geom, geom3, &
                             geom_mask, geom3_mask, up_mask )
   integer, intent(in) :: natom
   _REAL_, intent(out) :: ect 
   _REAL_, intent(inout) :: force(3,natom)
   _REAL_, intent(out) :: charge(natom)
   _REAL_, intent(inout) :: geom(4,cr_max_geom), geom3(13,cr_max_geom3) 
   logical, intent(in) :: geom_mask(cr_max_geom), geom3_mask(cr_max_geom3)
   logical, intent(out) :: up_mask(cr_max_order)

   integer :: i, j, ierr, a4, gi, ind_i, cr_type, a_mod
   logical :: hold, release, up, start_hold
   _REAL_ :: hold_charge
   _REAL_, allocatable, dimension(:,:) :: hold_dcdr
   integer :: hold_type, j_beg, j_end
   
#ifdef MPI
   integer :: up_taskid
#endif

   ! initialize
   charge = cr_charge
   cr_dcdr = 0.0

   do i = 1, cr_max_order
      if ( mod(i,numtasks) .ne. mytaskid ) cycle

      hold = .false.
      up = .false.
      j = cr_order(2,i)

      ! dcdr_i
      a_mod = cr_order(1,i)
      ind_i = cr_dcdr_tbl(a_mod)
      j_beg = cr_dcdr_i(ind_i,1)
      j_end = cr_dcdr_i(ind_i,2)

      do while ( j > 0 )
         gi = cr_info(8,j)
         a4 = cr_info(4,j)

         if ( a4 == 0 ) then
            if ( geom_mask(gi) ) then

               cr_type = cr_info(5,j)
               call cr_get_hold( cr_type, hold, release, start_hold, &
                                 hold_type, hold_charge )
               if ( start_hold ) then
                  allocate( hold_dcdr(j_beg:j_end,3), stat=ierr )
                  hold_dcdr = 0.0
               end if

               up = .true.
               if ( hold ) then
                  call cr_update_distance_charge( hold_charge, hold_dcdr, &
                     geom(1,gi), geom(2:4,gi), cr_info(1,j), cr_info(2,j), &
                     cr_type, ind_i, j_beg, j_end, &
                     force, ect, cr_cect(j,1:2) )
               else
                  call cr_update_distance_charge( charge(a_mod), &
                     cr_dcdr(j_beg:j_end,1:3), &
                     geom(1,gi), geom(2:4,gi), cr_info(1,j), cr_info(2,j), &
                     cr_type, ind_i, j_beg, j_end, & 
                     force, ect, cr_cect(j,1:2) )
               end if
            end if
         else
            ! angle
            if ( geom3_mask(gi) ) then

               cr_type = cr_info(5,j)
               call cr_get_hold( cr_type, hold, release, start_hold, &
                                 hold_type, hold_charge )
               if ( start_hold ) then
                  allocate( hold_dcdr(j_beg:j_end,3), stat=ierr )
                  hold_dcdr = 0.0
               end if

               up = .true.
               if ( hold ) then
                  call cr_update_angle_charge( hold_charge, hold_dcdr, &
                     geom3(1,gi), geom3(2:4,gi), geom3(5:7,gi), geom3(8:10,gi),&
                     cr_info(1,j), cr_info(2,j), cr_info(3,j), cr_type, &
                     ind_i, j_beg, j_end, force, ect, cr_cect(j,1:2) )
               else
                  call cr_update_angle_charge( charge(a_mod), &
                     cr_dcdr(j_beg:j_end,1:3), &
                     geom3(1,gi), geom3(2:4,gi), geom3(5:7,gi), geom3(8:10,gi),&
                     cr_info(1,j), cr_info(2,j), cr_info(3,j), cr_type, &
                     ind_i, j_beg, j_end, force, ect, cr_cect(j,1:2) )
               end if
            end if
         end if

         if ( hold .and. release ) then
            call cr_combine_charge( charge(a_mod), cr_dcdr, hold_charge, &
               hold_dcdr, hold_type, j_beg, j_end )
            deallocate( hold_dcdr, stat=ierr )
            hold = .false.
         end if

         j = cr_info(9,j)
      end do

      if ( hold ) then
         call cr_combine_charge( charge(a_mod), cr_dcdr, hold_charge, &
            hold_dcdr, hold_type, j_beg, j_end )
         deallocate( hold_dcdr, stat=ierr )
      end if

      up_mask(i) = up
   end do

#ifdef MPI
   ind_i = 1
   do i = 1, cr_max_order
      up_taskid = mod( i, numtasks )
      call mpi_bcast( up_mask(i), 1, MPI_LOGICAL, up_taskid, &
                      commsander, ierr )
      if ( .not. up_mask(i)) cycle

      a_mod = cr_order(1,i)
      call mpi_bcast( charge(a_mod), 1, MPI_DOUBLE_PRECISION, up_taskid, &
                      commsander, ierr )
      ind_i = cr_dcdr_tbl(a_mod)
      j_beg = cr_dcdr_i(ind_i,1) 
      j_end = cr_dcdr_i(ind_i,2) 

      call mpi_bcast( cr_dcdr(j_beg:j_end,1:3), &
                      ( j_end - j_beg + 1 ) * 3, &
                      MPI_DOUBLE_PRECISION, up_taskid, commsander, ierr )
   end do
#endif

end subroutine cr_update_charge
!===============================================================================


!===============================================================================
subroutine cr_get_hold( cr_type, hold, release, start_hold, &
                        hold_type, hold_charge )
   integer, intent(inout) :: cr_type
   logical, intent(inout) :: hold
   logical, intent(out) :: release, start_hold
   integer, intent(out) :: hold_type
   _REAL_, intent(inout) :: hold_charge

   integer :: ierr

   start_hold = .false.
   if ( cr_type > 5 ) then
      release = .false.
      if ( hold )  then
         cr_type = cr_type - 5
      else
         hold = .true. 
         hold_type = cr_type - 5
         cr_type = 5
         hold_charge = 0.0
         start_hold = .true.
      end if
   else
      release = .true.
   end if
end subroutine cr_get_hold
!===============================================================================


!===============================================================================
subroutine cr_combine_charge( charge, dcdr, new_c, new_dcdr, cr_type, &
                              j_beg, j_end )
   integer, intent(in) :: cr_type, j_beg, j_end
   _REAL_, intent(inout) :: charge
   _REAL_, dimension(:,:), intent(inout) :: dcdr
   _REAL_, intent(in) :: new_c
   _REAL_, dimension(j_beg:j_end,3), intent(in) :: new_dcdr

   _REAL_ :: old_c, ncinv
   integer :: j

   select case ( cr_type )
      case (1) ! summation
         charge = charge + new_c
         do j = j_beg, j_end
            dcdr(j,1:3) = dcdr(j,1:3) + new_dcdr(j,1:3)
         end do
      case (2) ! subtraction
         charge = charge - new_c
         do j = j_beg, j_end
            dcdr(j,1:3) = dcdr(j,1:3) - new_dcdr(j,1:3)
         end do
      case (3) ! multiplication
         old_c = charge
         charge = old_c * new_c
         do j = j_beg, j_end
            dcdr(j,1:3) = dcdr(j,1:3) * new_c + old_c * new_dcdr(j,1:3)
         end do
      case (4) ! division
         ncinv = 1.0 / new_c
         old_c = charge
         charge = old_c * ncinv
         do j = j_beg, j_end
            dcdr(j,1:3) = ( dcdr(j,1:3) - charge * new_dcdr(j,1:3) ) * ncinv
         end do
      case default
         write(6,'(x,a,i8)') 'CRGRELOC: unknown cr_type', cr_type
         call mexit(6,1)
   end select

end subroutine cr_combine_charge
!===============================================================================


!===============================================================================
subroutine cr_update_angle_charge( charge, dcdr, new_c, new_dcdr1, new_dcdr2, &
                        new_dcdr3, a1, a2, a3, cr_type, ind_i, &
                        j_beg, j_end, force, ect, cect )

   use constants, only : AMBER_ELECTROSTATIC  

   integer, intent(in) :: a1, a2, a3, cr_type, ind_i, j_beg, j_end
   _REAL_, intent(inout) :: charge
   _REAL_, intent(inout) :: force(3,*) ! 3,natom
   _REAL_, intent(inout) :: ect
   _REAL_, intent(inout) :: cect(2)
   _REAL_, dimension(j_beg:j_end,3), intent(inout) :: dcdr
   _REAL_, intent(inout) :: new_c, new_dcdr1(3), new_dcdr2(3), new_dcdr3(3)
   _REAL_ :: old_c, ncinv, tmp, tmp2, dedr1(3), dedr2(3), dedr3(3)
   integer :: ind1, ind2, ind3

   select case ( cr_type )
      case (1) ! summation
         new_c = new_c * AMBER_ELECTROSTATIC
         new_dcdr1 = new_dcdr1 * AMBER_ELECTROSTATIC
         new_dcdr2 = new_dcdr2 * AMBER_ELECTROSTATIC
         new_dcdr3 = new_dcdr3 * AMBER_ELECTROSTATIC
         charge = charge + new_c
         call cr_get_dcdr_index_j( ind_i, a1, ind1 )
         call cr_get_dcdr_index_j( ind_i, a2, ind2 )
         call cr_get_dcdr_index_j( ind_i, a3, ind3 )
         dcdr(ind1,1:3) = dcdr(ind1,1:3) + new_dcdr1(1:3)
         dcdr(ind2,1:3) = dcdr(ind2,1:3) + new_dcdr2(1:3)
         dcdr(ind3,1:3) = dcdr(ind3,1:3) + new_dcdr3(1:3)
         ect = ect + new_c * ( cect(1) + new_c * cect(2) )
         tmp = cect(1) + 2.0 * cect(2) * new_c
         force(1:3,a1) = force(1:3,a1) - new_dcdr1(1:3) * tmp
         force(1:3,a2) = force(1:3,a2) - new_dcdr2(1:3) * tmp
         force(1:3,a3) = force(1:3,a3) - new_dcdr3(1:3) * tmp
      case (2) ! subtraction
         new_c = new_c * AMBER_ELECTROSTATIC
         new_dcdr1 = new_dcdr1 * AMBER_ELECTROSTATIC
         new_dcdr2 = new_dcdr2 * AMBER_ELECTROSTATIC
         new_dcdr3 = new_dcdr3 * AMBER_ELECTROSTATIC
         charge = charge - new_c
         call cr_get_dcdr_index_j( ind_i, a1, ind1 )
         call cr_get_dcdr_index_j( ind_i, a2, ind2 )
         call cr_get_dcdr_index_j( ind_i, a3, ind3 )
         dcdr(ind1,1:3) = dcdr(ind1,1:3) - new_dcdr1(1:3)
         dcdr(ind2,1:3) = dcdr(ind2,1:3) - new_dcdr2(1:3)
         dcdr(ind3,1:3) = dcdr(ind3,1:3) - new_dcdr3(1:3)
         ect = ect - new_c * ( cect(1) - new_c * cect(2) )
         tmp = - cect(1) + 2.0 * cect(2) * new_c
         force(1:3,a1) = force(1:3,a1) - new_dcdr1(1:3) * tmp
         force(1:3,a2) = force(1:3,a2) - new_dcdr2(1:3) * tmp
         force(1:3,a3) = force(1:3,a3) - new_dcdr3(1:3) * tmp
      case (3) ! multiplication
         old_c = charge
         charge = old_c * new_c
         call cr_get_dcdr_index_j( ind_i, a1, ind1 )
         call cr_get_dcdr_index_j( ind_i, a2, ind2 )
         call cr_get_dcdr_index_j( ind_i, a3, ind3 )
         tmp = new_c - 1.0
         dedr1(1:3) = dcdr(ind1,1:3) * tmp + old_c * new_dcdr1(1:3)
         dedr2(1:3) = dcdr(ind2,1:3) * tmp + old_c * new_dcdr2(1:3)
         dedr3(1:3) = dcdr(ind3,1:3) * tmp + old_c * new_dcdr3(1:3)
         dcdr(ind1,1:3) = dedr1(1:3) + dcdr(ind1,1:3)
         dcdr(ind2,1:3) = dedr2(1:3) + dcdr(ind2,1:3)
         dcdr(ind3,1:3) = dedr3(1:3) + dcdr(ind3,1:3)
         tmp = charge - old_c
         ect = ect + tmp * ( cect(1) + tmp * cect(2) )
         tmp = cect(1) + 2.0 * cect(2) * tmp
         dedr1(1:3) = dedr1(1:3) * tmp
         dedr2(1:3) = dedr2(1:3) * tmp
         dedr3(1:3) = dedr3(1:3) * tmp
         force(1:3,a1) = force(1:3,a1) - dedr1(1:3)
         force(1:3,a2) = force(1:3,a2) - dedr2(1:3)
         force(1:3,a3) = force(1:3,a3) - dedr3(1:3)
      case (4) ! division
         ncinv = 1.0 / new_c
         old_c = charge
         charge = old_c * ncinv
         call cr_get_dcdr_index_j( ind_i, a1, ind1 )
         call cr_get_dcdr_index_j( ind_i, a2, ind2 )
         call cr_get_dcdr_index_j( ind_i, a3, ind3 )
         tmp = ncinv - 1.0
         tmp2 = charge * ncinv
         dedr1(1:3) = dcdr(ind1,1:3) * tmp - tmp2 * new_dcdr1(1:3)
         dedr2(1:3) = dcdr(ind2,1:3) * tmp - tmp2 * new_dcdr2(1:3)
         dedr3(1:3) = dcdr(ind3,1:3) * tmp - tmp2 * new_dcdr3(1:3)
         dcdr(ind1,1:3) = dedr1(1:3) + dcdr(ind1,1:3)
         dcdr(ind2,1:3) = dedr2(1:3) + dcdr(ind2,1:3)
         dcdr(ind3,1:3) = dedr3(1:3) + dcdr(ind3,1:3)
         tmp = tmp * old_c
         ect = ect + tmp * ( cect(1) + tmp * cect(2) )
         tmp = cect(1) + 2.0 * cect(2) * tmp
         dedr1(1:3) = dedr1(1:3) * tmp
         dedr2(1:3) = dedr2(1:3) * tmp
         dedr3(1:3) = dedr3(1:3) * tmp
         force(1:3,a1) = force(1:3,a1) - dedr1(1:3)
         force(1:3,a2) = force(1:3,a2) - dedr2(1:3)
         force(1:3,a3) = force(1:3,a3) - dedr3(1:3)
      case (5) ! reassignment
         new_c = new_c * AMBER_ELECTROSTATIC
         new_dcdr1 = new_dcdr1 * AMBER_ELECTROSTATIC
         new_dcdr2 = new_dcdr2 * AMBER_ELECTROSTATIC
         new_dcdr3 = new_dcdr3 * AMBER_ELECTROSTATIC
         tmp = new_c - charge
         charge = new_c
         ect = ect + tmp * ( cect(1) + tmp * cect(2) )
         tmp = cect(1) + 2.0 * cect(2) * tmp
         call cr_get_dcdr_index_j( ind_i, a1, ind1 )
         call cr_get_dcdr_index_j( ind_i, a2, ind2 )
         call cr_get_dcdr_index_j( ind_i, a3, ind3 )
         force(1:3,a1) =   force(1:3,a1) &
                         - tmp * ( new_dcdr1(1:3) - dcdr(ind1,1:3) )
         force(1:3,a2) =   force(1:3,a2) &
                         - tmp * ( new_dcdr2(1:3) - dcdr(ind2,1:3) )
         force(1:3,a3) =   force(1:3,a3) &
                         - tmp * ( new_dcdr3(1:3) - dcdr(ind3,1:3) )
         dcdr(ind1,1:3) = new_dcdr1(1:3)
         dcdr(ind2,1:3) = new_dcdr2(1:3)
         dcdr(ind3,1:3) = new_dcdr3(1:3)
      case default
         write(6,'(x,a,i8)') 'CRGRELOC: unknown cr_type', cr_type
         call mexit(6,1)
   end select

end subroutine cr_update_angle_charge
!===============================================================================


!===============================================================================
subroutine cr_update_distance_charge( charge, dcdr, new_c, new_dcdr, a1, a2, &
                                      cr_type, ind_i, j_beg, j_end, &
                                      force, ect, cect )

   use constants, only : AMBER_ELECTROSTATIC  

   integer, intent(in) :: a1, a2, cr_type, ind_i, j_beg, j_end
   _REAL_, intent(inout) :: charge
   _REAL_, intent(inout), dimension(j_beg:j_end,3) :: dcdr
   _REAL_, intent(inout) :: new_c, new_dcdr(3)
   _REAL_, intent(inout) :: force(3,*) ! 3, natom
   _REAL_, intent(inout) :: ect
   _REAL_, intent(in) :: cect(2)
   _REAL_ :: prod(3), dedr1(3), dedr2(3), old_c, ncinv, tmp
   integer :: ind1, ind2, j, atr

   select case ( cr_type )
      ! r is increasing from a1 to a2
      ! so the sign of new_dcdr should be swiched for ind1
      case (1) ! summation
         new_c = new_c * AMBER_ELECTROSTATIC
         new_dcdr = new_dcdr * AMBER_ELECTROSTATIC
         charge = charge + new_c
         ! get index (of cr_dcdr_j) for atom a1
         call cr_get_dcdr_index_j( ind_i, a1, ind1 )
         ! get index (of cr_dcdr_j) for atom a2
         call cr_get_dcdr_index_j( ind_i, a2, ind2 )
         ! update dcdr for dq/dr1 and dq/dr2
         dcdr(ind1,1:3) = dcdr(ind1,1:3) - new_dcdr(1:3) ! dq3/dr1
         dcdr(ind2,1:3) = dcdr(ind2,1:3) + new_dcdr(1:3) ! dq3/dr2
         prod(1:3) = new_dcdr(1:3) * ( cect(1) + 2.0 * cect(2) * new_c )
         ect = ect + new_c * ( cect(1) + new_c * cect(2) )
         force(1:3,a1) = force(1:3,a1) + prod(1:3) ! add -dE/dr1
         force(1:3,a2) = force(1:3,a2) - prod(1:3) ! add -dE/dr2
      case (2) ! subtraction
         new_c = new_c * AMBER_ELECTROSTATIC
         new_dcdr = new_dcdr * AMBER_ELECTROSTATIC
         charge = charge - new_c
         call cr_get_dcdr_index_j( ind_i, a1, ind1 )
         call cr_get_dcdr_index_j( ind_i, a2, ind2 )
         dcdr(ind1,1:3) = dcdr(ind1,1:3) + new_dcdr(1:3)
         dcdr(ind2,1:3) = dcdr(ind2,1:3) - new_dcdr(1:3)
         prod(1:3) = new_dcdr(1:3) * ( - cect(1) + 2.0 * cect(2) * new_c )
         ect = ect - new_c * ( cect(1) - new_c * cect(2) )
         force(1:3,a1) = force(1:3,a1) + prod(1:3)
         force(1:3,a2) = force(1:3,a2) - prod(1:3)
      case (3) ! multiplication
         old_c = charge
         charge = old_c * new_c
         prod(1:3) = old_c * new_dcdr(1:3)
         call cr_get_dcdr_index_j( ind_i, a1, ind1 )
         call cr_get_dcdr_index_j( ind_i, a2, ind2 )
         tmp = new_c - 1.0
         dedr1(1:3) = dcdr(ind1,1:3) * tmp - prod(1:3)
         dedr2(1:3) = dcdr(ind2,1:3) * tmp + prod(1:3)
         dcdr(ind1,1:3) = dedr1(1:3) + dcdr(ind1,1:3)
         dcdr(ind2,1:3) = dedr2(1:3) + dcdr(ind2,1:3)
         tmp = charge - old_c
         ect = ect + tmp * ( cect(1) + tmp * cect(2) )
         tmp = cect(1) + 2.0 * cect(2) * tmp
         dedr1(1:3) = dedr1(1:3) * tmp
         dedr2(1:3) = dedr2(1:3) * tmp
         force(1:3,a1) = force(1:3,a1) - dedr1(1:3)
         force(1:3,a2) = force(1:3,a2) - dedr2(1:3)
      case (4) ! division
         ncinv = 1.0 / new_c
         old_c = charge
         charge = old_c * ncinv
         call cr_get_dcdr_index_j( ind_i, a1, ind1 )
         call cr_get_dcdr_index_j( ind_i, a2, ind2 )
         prod(1:3) = charge * ncinv * new_dcdr(1:3)
         tmp = ncinv - 1.0
         dedr1(1:3) = dcdr(ind1,1:3) * tmp + prod(1:3)
         dedr2(1:3) = dcdr(ind2,1:3) * tmp - prod(1:3)
         dcdr(ind1,1:3) = dedr1(1:3) + dcdr(ind1,1:3)
         dcdr(ind2,1:3) = dedr2(1:3) + dcdr(ind2,1:3)
         tmp = tmp * old_c
         ect = ect + tmp * ( cect(1) + tmp * cect(2) )
         tmp = cect(1) + 2.0 * cect(2) * tmp
         dedr1(1:3) = dedr1(1:3) * tmp
         dedr2(1:3) = dedr2(1:3) * tmp
         force(1:3,a1) = force(1:3,a1) - dedr1(1:3)
         force(1:3,a2) = force(1:3,a2) - dedr2(1:3)
      case (5) ! reassignment
         new_c = new_c * AMBER_ELECTROSTATIC
         new_dcdr = new_dcdr * AMBER_ELECTROSTATIC
         tmp = new_c - charge
         charge = new_c
         ect = ect + tmp * ( cect(1) + tmp * cect(2) )
         tmp = cect(1) + 2.0 * cect(2) * tmp
         call cr_get_dcdr_index_j( ind_i, a1, ind1 )
         call cr_get_dcdr_index_j( ind_i, a2, ind2 )
         force(1:3,a1) =   force(1:3,a1) &
                         - tmp * ( -new_dcdr(1:3) - dcdr(ind1,1:3) )
         force(1:3,a2) =   force(1:3,a2) &
                         - tmp * (  new_dcdr(1:3) - dcdr(ind2,1:3) )
         dcdr(ind1,1:3) = -new_dcdr(1:3)
         dcdr(ind2,1:3) = new_dcdr(1:3)
      case default
         write(6,'(x,a,i8)') 'CRGRELOC: unknown cr_type', cr_type
         call mexit(6,1)
   end select

end subroutine cr_update_distance_charge
!===============================================================================


!===============================================================================
subroutine cr_fill_geom3( crd, geom3, natom )

#include "box.h"

   _REAL_, intent(in) :: crd(3,natom)
   _REAL_, intent(out) :: geom3(13,cr_max_geom3)
   integer, intent(in) :: natom
   _REAL_ :: half_box_size2

   half_box_size2 =  sum( box * box ) * 0.25

   if ( cr_update_tb ) then
      call cr_fill_geom3_update( crd, geom3, natom, half_box_size2 )
   else
      call cr_fill_geom3_noupdate( crd, geom3, natom, half_box_size2 )
   end if
end subroutine cr_fill_geom3
!===============================================================================


!===============================================================================
subroutine cr_fill_geom( crd, geom, natom )

#include "box.h"

   _REAL_, intent(in) :: crd(3,natom)
   _REAL_, intent(out) :: geom(4,cr_max_geom)
   integer, intent(in) :: natom
   _REAL_ :: half_box_size2

   half_box_size2 =  sum( box * box ) * 0.25

   if ( cr_update_pair ) then
      call cr_fill_geom_update( crd, geom, natom, half_box_size2 )
   else
      call cr_fill_geom_noupdate( crd, geom, natom, half_box_size2 )
   end if
end subroutine cr_fill_geom
!===============================================================================


!===============================================================================
subroutine cr_fill_geom3_noupdate( crd, geom3, natom, half_box_size2 )
   use nblist, only: ucell

   _REAL_, intent(in) :: crd(3,natom)
   _REAL_, intent(out) :: geom3(13,cr_max_geom3)
   integer, intent(in) :: natom
   _REAL_, intent(in) :: half_box_size2
   _REAL_ :: half_box_size
   integer :: ai, aj, ak, al, am, i, j, k, ierr, gi
   _REAL_ :: dist1, delr12, delx1, dely1, delz1  ! j - i
   _REAL_ :: dist2, delr22, delx2, dely2, delz2  ! k - j
   _REAL_ :: dist3, delr32, delx3, dely3, delz3  ! k - i
   _REAL_ :: cos_angle
#ifdef MPI
   integer :: jobcount, up_taskid
   logical :: lor
   logical :: up_mask(cr_max_info3_k) 
#endif

   geom3 = 0.0
   half_box_size = sqrt(half_box_size2)

#ifdef MPI
   jobcount = 1
#endif

   do ai = 1,cr_max_info3_i
      i = cr_info3_i(ai,1)
      do aj = cr_info3_i(ai,2), cr_info3_i(ai,3)
         j = cr_info3_j(aj,1)
         do ak = cr_info3_j(aj,2), cr_info3_j(aj,3)
            k = cr_info3_k(ak,1)
            if ( .not. cr_tb_eval(ak) ) cycle
#ifdef MPI
            if ( mod(jobcount,numtasks) .ne. mytaskid ) then
               jobcount = jobcount + 1
               cycle
            else
               jobcount = jobcount + 1
            end if
            up_mask(ak) = .false.
#endif
            ! calculate distance
            delx1 =   crd(1,j) - crd(1,i) &
                    - sum( cr_tb_tranunit(1:3,ak) * ucell(1,1:3) ) 
            dely1 =   crd(2,j) - crd(2,i) &
                    - sum( cr_tb_tranunit(1:3,ak) * ucell(2,1:3) ) 
            delz1 =   crd(3,j) - crd(3,i) &
                    - sum( cr_tb_tranunit(1:3,ak) * ucell(3,1:3) ) 
            dist1 = sqrt( delx1 * delx1 + dely1 * dely1 + delz1 * delz1 )

            delx2 =   crd(1,k) - crd(1,j) &
                    - sum( cr_tb_tranunit(4:6,ak) * ucell(1,1:3) ) 
            dely2 =   crd(2,k) - crd(2,j) &
                    - sum( cr_tb_tranunit(4:6,ak) * ucell(2,1:3) ) 
            delz2 =   crd(3,k) - crd(3,j) &
                    - sum( cr_tb_tranunit(4:6,ak) * ucell(3,1:3) ) 
            dist2 = sqrt( delx2 * delx2 + dely2 * dely2 + delz2 * delz2 )

            delx3 =   crd(1,k) - crd(1,i) &
                    - sum( cr_tb_tranunit(7:9,ak) * ucell(1,1:3) ) 
            dely3 =   crd(2,k) - crd(2,i) &
                    - sum( cr_tb_tranunit(7:9,ak) * ucell(2,1:3) ) 
            delz3 =   crd(3,k) - crd(3,i) &
                    - sum( cr_tb_tranunit(7:9,ak) * ucell(3,1:3) ) 
            dist3 = sqrt( delx3 * delx3 + dely3 * dely3 + delz3 * delz3 )

            ! check if pair list update is necessary
            if ( abs( cr_tb_distance(1,ak) - dist1 ) > crskin ) then
               call cr_measure_min_distance( crd(1:3,i), crd(1:3,j), &
                  cr_tb_tranunit(1:3,ak), half_box_size2, delr12, delx1, dely1, &
                  delz1 )
               dist1 = sqrt(delr12)
               if ( abs( cr_tb_distance(1,ak) - dist1 ) > half_crskin ) then
                  cr_update_tb = .true.
               else
#ifdef MPI
                  up_mask(ak) = .true.
#endif
               end if
            else if ( abs( cr_tb_distance(1,ak) - dist1 ) > half_crskin ) then
               cr_update_tb = .true.
            end if
            if ( abs( cr_tb_distance(2,ak) - dist2 ) > crskin ) then
               call cr_measure_min_distance( crd(1:3,j), crd(1:3,k), &
                  cr_tb_tranunit(4:6,ak), half_box_size2, delr22, delx2, dely2, &
                  delz2 )
               dist2 = sqrt(delr22)
               cr_tb_distance(2,ak) = dist2
#ifdef MPI
               up_mask(ak) = .true.
#endif
            end if
            if ( abs( cr_tb_distance(3,ak) - dist3 ) > crskin ) then
               call cr_measure_min_distance( crd(1:3,i), crd(1:3,k), &
                  cr_tb_tranunit(7:9,ak), half_box_size2, delr32, delx3, dely3, &
                  delz3 )
               dist3 = sqrt(delr32)
               cr_tb_distance(3,ak) = dist3
#ifdef MPI
               up_mask(ak) = .true.
#endif
            end if

            call cr_calc_angle( dist1, dist2, dist3, cos_angle )

            do al = cr_info3_k(ak,2), cr_info3_k(ak,3)
               am = cr_info3_ptr(al)
               gi = cr_info(8,am)
               geom3(1,gi) = cos_angle
               geom3(2,gi) = dist1
               geom3(3,gi) = delx1 
               geom3(4,gi) = dely1
               geom3(5,gi) = delz1
               geom3(6,gi) = dist2
               geom3(7,gi) = delx2 
               geom3(8,gi) = dely2
               geom3(9,gi) = delz2
               geom3(10,gi) = dist3
               geom3(11,gi) = delx3 
               geom3(12,gi) = dely3
               geom3(13,gi) = delz3
            end do
         end do
      end do
   end do
#ifdef MPI
   call mpi_allreduce( cr_update_tb, lor, 1, MPI_LOGICAL, MPI_LOR, &
                       commsander, ierr )
   cr_update_tb = lor

   jobcount = 1
   do ai = 1,cr_max_info3_i
      do aj = cr_info3_i(ai,2), cr_info3_i(ai,3)
         do ak = cr_info3_j(aj,2), cr_info3_j(aj,3)
            if ( .not. cr_tb_eval(ak)) cycle
            up_taskid = mod(jobcount,numtasks)
            jobcount = jobcount + 1

            call mpi_bcast( up_mask(ak), 1, MPI_LOGICAL, up_taskid, &
                            commsander, ierr )
            do al = cr_info3_k(ak,2), cr_info3_k(ak,3)
               am = cr_info3_ptr(al)
               gi = cr_info(8,am)
               call mpi_bcast( geom3(1:13,gi), 13, MPI_DOUBLE_PRECISION, &
                               up_taskid, commsander, ierr )
            end do

            if ( .not. up_mask(ak)) cycle 
            call mpi_bcast( cr_tb_tranunit(1:9,ak), 9, MPI_INTEGER, up_taskid, &
                            commsander, ierr ) 
            call mpi_bcast( cr_tb_distance(2:3,ak), 2, MPI_DOUBLE_PRECISION, &
                            up_taskid, commsander, ierr ) 
         end do
      end do
   end do
#endif

end subroutine cr_fill_geom3_noupdate
!===============================================================================


!===============================================================================
subroutine cr_fill_geom_noupdate( crd, geom, natom, half_box_size2 )
   use nblist, only: ucell

   _REAL_, intent(in) :: crd(3,natom)
   _REAL_, intent(out) :: geom(4,cr_max_geom)
   integer, intent(in) :: natom
   _REAL_, intent(in) :: half_box_size2
   _REAL_ :: half_box_size
   integer :: ai, aj, ak, al, i, j, ierr, gi
   _REAL_ :: delr, delr2, delx, dely, delz
#ifdef MPI
   integer :: jobcount, up_taskid
   logical :: lor
   logical :: up_mask(cr_max_info_j) 
#endif

   geom = 0.0
   half_box_size = sqrt(half_box_size2)

#ifdef MPI
   jobcount = 1
#endif

   do ai = 1,cr_max_info_i
      i = cr_info_i( ai, 1 )
      do aj = cr_info_i( ai, 2 ), cr_info_i( ai, 3 )
         if ( .not. cr_pair_eval(aj) ) cycle
         j = cr_info_j( aj, 1 )
#ifdef MPI
         if ( mod(jobcount,numtasks) .ne. mytaskid ) then
            jobcount = jobcount + 1
            cycle
         else
            jobcount = jobcount + 1
         end if
         up_mask(aj) = .false.
#endif
         ! calculate distance
         delx =   crd(1,j) - crd(1,i) &
               - sum( cr_tranunit(1:3,aj) * ucell(1,1:3) ) 
         dely =   crd(2,j) - crd(2,i) &
               - sum( cr_tranunit(1:3,aj) * ucell(2,1:3) ) 
         delz =   crd(3,j) - crd(3,i) &
               - sum( cr_tranunit(1:3,aj) * ucell(3,1:3) ) 
         delr = sqrt( delx * delx + dely * dely + delz * delz )

         ! check if pair list update is necessary
         if ( abs( cr_pair_distance(aj) - delr ) > crskin ) then
            ! If distance difference is bigger than crskin (eg iwrap=1),
            ! this is possibly caused by wrapping molecules.  So first try
            ! to tune only tranunit.
            ! If successful, keep going on. If not, flag cr_update_pair.
            call cr_measure_min_distance( crd(1:3,i), crd(1:3,j), &
               cr_tranunit(1:3,aj), half_box_size2, delr2, delx, dely, delz )
            delr = sqrt(delr2)

            if ( abs( cr_pair_distance(aj) - delr ) > half_crskin ) then
               cr_update_pair = .true.
            else
#ifdef MPI
               up_mask(aj) = .true.
#endif
            end if
         else if ( abs( cr_pair_distance(aj) - delr ) > half_crskin ) then
            cr_update_pair = .true.
         end if

         do ak = cr_info_j( aj, 2 ), cr_info_j( aj, 3 )
            al = cr_info_ptr(ak)
            gi = cr_info(8,al) 
            geom(1,gi) = delr
            geom(2,gi) = delx
            geom(3,gi) = dely
            geom(4,gi) = delz
         end do

      end do
   end do
#ifdef MPI
   call mpi_allreduce( cr_update_pair, lor, 1, MPI_LOGICAL, MPI_LOR, &
                       commsander, ierr )
   cr_update_pair = lor

   jobcount = 1
   do ai = 1,cr_max_info_i
      do aj = cr_info_i( ai, 2 ), cr_info_i( ai, 3 )
         if ( .not.cr_pair_eval(aj)) cycle
         up_taskid = mod(jobcount,numtasks)
         jobcount = jobcount + 1

         do ak = cr_info_j(aj,2), cr_info_j(aj,3)
            al = cr_info_ptr(ak)
            gi = cr_info(8,al)
            call mpi_bcast( geom(1:4,gi), 4, MPI_DOUBLE_PRECISION, &
                            up_taskid, commsander, ierr )   
         end do

         call mpi_bcast( up_mask(aj), 1, MPI_LOGICAL, up_taskid, &
                         commsander, ierr )
         if ( .not. up_mask(aj)) cycle 
         call mpi_bcast( cr_tranunit(1:3,aj), 3, MPI_INTEGER, up_taskid, &
                         commsander, ierr ) 
      end do
   end do
#endif

end subroutine cr_fill_geom_noupdate
!===============================================================================


!===============================================================================
subroutine cr_calc_angle( r1, r2, r3, cos_angle )
   _REAL_, intent(in) :: r1, r2, r3
   _REAL_, intent(out) :: cos_angle

   cos_angle = ( r1*r1 + r2*r2 - r3*r3 ) / ( 2.0*r1*r2 )
end subroutine cr_calc_angle
!===============================================================================


!===============================================================================
subroutine cr_fill_geom3_update( crd, geom3, natom, half_box_size2 )

   _REAL_, intent(in) :: crd(3,natom)
   _REAL_, intent(out) :: geom3(13,cr_max_geom3)
   integer, intent(in) :: natom
   _REAL_, intent(in) :: half_box_size2
   integer :: ai, aj, ak, al, am, i, j, k, ierr, gi
   _REAL_ :: dist1, delr12, delx1, dely1, delz1  ! j - i
   _REAL_ :: dist2, delr22, delx2, dely2, delz2  ! k - j
   _REAL_ :: dist3, delr32, delx3, dely3, delz3  ! k - i
   _REAL_ :: cos_angle
#ifdef MPI
   integer :: up_taskid
   logical :: land
#endif

   ! initialize
   geom3 = 0.0

   ! update three-body list
   do ai = 1,cr_max_info3_i
      i = cr_info3_i(ai,1)
      do aj = cr_info3_i(ai,2), cr_info3_i(ai,3)
         j = cr_info3_j(aj,1)
         do ak = cr_info3_j(aj,2), cr_info3_j(aj,3)
            k = cr_info3_k(ak,1)
#ifdef MPI
            if ( mod(ak,numtasks) .ne. mytaskid ) cycle
#endif
            ! calculate distance
            call cr_measure_min_distance( crd(1:3,i), crd(1:3,j), &
               cr_tb_tranunit(1:3,ak), half_box_size2, delr12, &
               delx1, dely1, delz1 )
            call cr_measure_min_distance( crd(1:3,j), crd(1:3,k), &
               cr_tb_tranunit(4:6,ak), half_box_size2, delr22, &
               delx2, dely2, delz2 )
            call cr_measure_min_distance( crd(1:3,i), crd(1:3,k), &
               cr_tb_tranunit(7:9,ak), half_box_size2, delr32, &
               delx3, dely3, delz3 )

            if (       delr12 < crcut_skin2 &
                 .and. delr22 < crcut_skin2 &
                 .and. delr32 < crcut_skin2 &
               ) &
            then
               cr_tb_eval(ak) = .true.
               dist1 = sqrt(delr12)
               dist2 = sqrt(delr22)
               dist3 = sqrt(delr32)
               cr_tb_distance(1,ak) = dist1
               cr_tb_distance(2,ak) = dist2
               cr_tb_distance(3,ak) = dist3
               call cr_calc_angle( dist1, dist2, dist3, cos_angle )
               do al = cr_info3_k( ak, 2 ), cr_info3_k( ak, 3 )
                  am = cr_info3_ptr(al)
                  gi = cr_info(8,am)
                  geom3(1,gi) = cos_angle
                  geom3(2,gi) = dist1
                  geom3(3,gi) = delx1 
                  geom3(4,gi) = dely1
                  geom3(5,gi) = delz1
                  geom3(6,gi) = dist2
                  geom3(7,gi) = delx2 
                  geom3(8,gi) = dely2
                  geom3(9,gi) = delz2
                  geom3(10,gi) = dist3
                  geom3(11,gi) = delx3 
                  geom3(12,gi) = dely3
                  geom3(13,gi) = delz3
               end do
               cr_update_tb = .false.
            else
               cr_tb_eval(ak) = .false.
               cr_tb_distance(1:3,ak) = 0.0
            end if
         end do
      end do
   end do

#ifdef MPI
   ! bcast
   do ak = 1, cr_max_info3_k
      up_taskid = mod(ak,numtasks)
      call mpi_bcast( cr_tb_eval(ak), 1, MPI_LOGICAL, up_taskid, &
                      commsander, ierr )
      if ( .not. cr_tb_eval(ak)) cycle

      do al = cr_info3_k(ak,2), cr_info3_k(ak,3)
         am = cr_info3_ptr(al)
         gi = cr_info(8,am)
         call mpi_bcast( geom3(1:13,gi), 13, MPI_DOUBLE_PRECISION, &
                         up_taskid, commsander, ierr )
      end do
      call mpi_bcast( cr_tb_distance(1:3,ak), 3, MPI_DOUBLE_PRECISION,&
                      up_taskid, commsander, ierr )
      call mpi_bcast( cr_tb_tranunit(1:9,ak), 9, MPI_INTEGER, &
                      up_taskid, commsander, ierr )   
   end do
   call mpi_allreduce( cr_update_tb, land, 1, MPI_LOGICAL, MPI_LAND, &
                       commsander, ierr )
   cr_update_tb = land
#endif /* MPI */

end subroutine cr_fill_geom3_update
!===============================================================================


!===============================================================================
subroutine cr_fill_geom_update( crd, geom, natom, half_box_size2 )

   _REAL_, intent(in) :: crd(3,natom)
   _REAL_, intent(out) :: geom(4,cr_max_geom)
   integer, intent(in) :: natom
   _REAL_, intent(in) :: half_box_size2
   integer :: ai, aj, ak, al, i, j, ierr, gi
   _REAL_ :: delr, delr2, delx, dely, delz
#ifdef MPI
   integer :: up_taskid
   logical :: land
#endif

   ! initialize
   geom = 0.0

   ! update pair list
   do ai = 1,cr_max_info_i
      i = cr_info_i( ai, 1 )
      do aj = cr_info_i( ai, 2 ), cr_info_i( ai, 3 )
         j = cr_info_j( aj, 1 )
#ifdef MPI
         if ( mod(aj,numtasks) .ne. mytaskid ) cycle
#endif
         ! calculate distance
         call cr_measure_min_distance( crd(1:3,i), crd(1:3,j), &
            cr_tranunit(1:3,aj), half_box_size2, delr2, delx, dely, delz )

         if ( delr2 < crcut_skin2 ) then
            cr_pair_eval(aj) = .true.
            delr = sqrt(delr2)
            cr_pair_distance(aj) = delr
            do ak = cr_info_j( aj, 2 ), cr_info_j( aj, 3 )
               al = cr_info_ptr(ak)
               gi = cr_info(8,al)
               geom(1,gi) = delr
               geom(2,gi) = delx
               geom(3,gi) = dely
               geom(4,gi) = delz
            end do
            cr_update_pair = .false.
         else
            cr_pair_distance(aj) = 0.0
            cr_pair_eval(aj) = .false.
         end if
      end do
   end do

#ifdef MPI
   ! bcast
   do aj = 1, cr_max_info_j
      up_taskid = mod( aj, numtasks )
      call mpi_bcast( cr_pair_eval(aj), 1, MPI_LOGICAL, up_taskid, &
                      commsander, ierr )
      if ( .not. cr_pair_eval(aj)) cycle

      do ak = cr_info_j(aj,2), cr_info_j(aj,3)
         al = cr_info_ptr(ak)
         gi = cr_info(8,al)
         call mpi_bcast( geom(1:4,gi), 4, MPI_DOUBLE_PRECISION, &
                         up_taskid, commsander, ierr )   
      end do
      call mpi_bcast( cr_pair_distance(aj), 1, MPI_DOUBLE_PRECISION, &
                      up_taskid, commsander, ierr )
      call mpi_bcast( cr_tranunit(1:3,aj), 3, MPI_INTEGER, up_taskid, &
                      commsander, ierr )
   end do
   call mpi_allreduce( cr_update_pair, land, 1, MPI_LOGICAL, MPI_LAND, &
                       commsander, ierr )
   cr_update_pair = land
#endif /* MPI */

end subroutine cr_fill_geom_update
!===============================================================================


!===============================================================================
subroutine cr_measure_min_distance( ri, rj, tunit, half_box_size2, delr2, &
                                    delx, dely, delz )
   use nblist, only: recip, ucell

   _REAL_, intent(in) :: ri(3), rj(3)
   integer, intent(out) :: tunit(3)
   _REAL_, intent(in) :: half_box_size2
   _REAL_, intent(out) :: delr2, delx, dely, delz
   _REAL_ :: dr(3), drc(3)
   integer :: cor(3), ccor(3)
   _REAL_ :: mindelr2
   integer :: jx, jy, jz
   logical :: exit_loop

   dr = rj - ri

   mindelr2 = half_box_size2 * 4.001
   cor = floor( recip(1,1:3)*dr(1) + recip(2,1:3)*dr(2) + recip(3,1:3)*dr(3) )
   exit_loop = .false.

   do jx = -1, 1
      if ( exit_loop ) exit
      do jy = -1, 1
         if ( exit_loop ) exit
         do jz = -1, 1
            ccor(1) = cor(1) + jx
            ccor(2) = cor(2) + jy
            ccor(3) = cor(3) + jz

            drc = dr - (   ccor(1) * ucell(1:3,1) &
                         + ccor(2) * ucell(1:3,2) &
                         + ccor(3) * ucell(1:3,3) &
                       ) 
            delr2 = sum( drc * drc )
            if ( mindelr2 > delr2 ) then
               mindelr2 = delr2
               tunit = ccor 
               delx = drc(1)
               dely = drc(2)
               delz = drc(3) 
            end if
            if ( mindelr2 < half_box_size2 ) then
               exit_loop = .true.
               exit
            end if
         end do
      end do
   end do

end subroutine cr_measure_min_distance
!===============================================================================


!===============================================================================
! read the original charge and store them into cr_charge
subroutine cr_backup_charge( charge, natom )
   integer, intent(in) :: natom
   _REAL_, intent(in) :: charge(natom)
   integer :: ierror

   allocate( cr_charge(natom), stat=ierror)
   if ( ierror /= 0 ) then
      write(6,'(x,a)') 'CRGRELOC: failed to allocate cr_charge'
      call mexit(6,1)
   end if

   cr_charge = charge
end subroutine cr_backup_charge
!===============================================================================


!===============================================================================
subroutine cr_get_dcdr_index_j( i, a2, ind )
   integer, intent(in) :: i, a2
   integer, intent(out) :: ind
   integer :: j

   ind = 0
   do j = cr_dcdr_i(i,1), cr_dcdr_i(i,2)
      if ( cr_dcdr_j(j) < a2 ) cycle
      if ( cr_dcdr_j(j) > a2 ) exit
      ind = j 
      exit
   end do
end subroutine cr_get_dcdr_index_j
!===============================================================================


!===============================================================================
subroutine cr_add_dcdr_factor( at, factor )
   integer, intent(in) :: at
   _REAL_, intent(in) :: factor
   integer :: i

   i = cr_dcdr_tbl(at)
   if ( i /= 0 ) then
      cr_dcdr_fac(i) = cr_dcdr_fac(i) + factor
   end if
end subroutine cr_add_dcdr_factor
!===============================================================================


!===============================================================================
subroutine cr_calc_force( force )
   _REAL_, intent(inout) :: force(3,*)
   _REAL_ :: factor
   integer :: i, j, atr

   do i = 1, cr_max_dcdr_i
      factor = cr_dcdr_fac(i)
      if ( factor == 0.0 ) cycle
      do j = cr_dcdr_i(i,1), cr_dcdr_i(i,2)
         atr = cr_dcdr_j(j)
         force(1:3,atr) = force(1:3,atr) - cr_dcdr(j,1:3) * factor
      end do
   end do
end subroutine cr_calc_force
!===============================================================================


end module crg_reloc
