#include "copyright.i"

!*******************************************************************************
!
! Module: prfs_mod
!
! Description: Support for pseudo-residue fragments, and molecules described as
!              a list of prf's.  Only needed under MPI.
!*******************************************************************************

module prfs_mod

  use gbl_datatypes_mod

  implicit none

#if defined(MPI)

  ! Pseudo-residue fragment data, used for PME and GB:

  integer, save :: gbl_prf_cnt = 0
  integer, save :: gbl_max_prf_atms = 0

  type(listdata_rec), allocatable, save :: gbl_prf_listdata(:)
  integer, allocatable, save            :: gbl_prf_lists(:)

contains

!*******************************************************************************
!
! Subroutine:  setup_prfs
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine setup_prfs(num_ints, num_reals, atm_cnt, res_cnt, res_atms, &
                      bond_cnt, bonds, frame_cnt, frames)

  use extra_pnts_nb14_mod
  use gbl_datatypes_mod
  use parallel_dat_mod
  use pmemd_lib_mod
#ifdef PRF_DBG
  use prmtop_dat_mod
#endif /* PRF_DBG */

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)               :: num_ints, num_reals
  integer, intent(in)                   :: atm_cnt
  integer, intent(in)                   :: res_cnt
  integer, intent(in)                   :: res_atms(res_cnt + 1)
  integer, intent(in)                   :: bond_cnt
  type(bond_rec), intent(in)            :: bonds(*)
  integer, intent(in)                   :: frame_cnt
  type(ep_frame_rec), intent(in)        :: frames(*)

! Local variables:

  integer               :: alloc_failed
  integer               :: i
  integer               :: max_res_atms
  integer               :: atm_lst_cnt
  integer               :: new_prf_id
  integer               :: cur_prf_id
  integer               :: min_prf_id
  integer               :: max_prf_id
  integer               :: frame_id
  integer               :: bond_id
  integer               :: res_id
  integer               :: atm_id, last_atm_id
  integer               :: prf_id
  integer               :: offset, old_offset, new_offset
  integer               :: listcnt
  integer               :: max_atm_lst_cnt
  integer, allocatable  :: atm_lst(:)
  integer               :: prf_lst_head(atm_cnt)
  integer               :: prf_lst_tail(atm_cnt)
  integer               :: prf_lst(atm_cnt)
  integer               :: prf_map(atm_cnt)
  type(listdata_rec)    :: prf_listdata(atm_cnt)
  integer               :: prf_lists(atm_cnt)

! Pre-initialize as needed:

  prf_lst_head(:) = 0
  prf_lst_tail(:) = 0
  prf_lst(:) = 0        ! Marks empty prf list entries.
  prf_map(:) = 0        ! None of the atoms yet in pseudo-residue fragments.

  max_res_atms = get_max_res_atms()
  max_atm_lst_cnt = max(max_res_atms, 6)        ! 6 atoms in ep frame,
                                                ! 2 atoms in bond.

  allocate(atm_lst(max_atm_lst_cnt), stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  max_prf_id = atm_cnt  ! max_prf_id is the maximum possible prf_id, not the
                        ! real maximum prf_id here...

  new_prf_id = 0

  do frame_id = 1, frame_cnt

    call make_frame_atm_lst(frames(frame_id), atm_lst_cnt, atm_lst)

    min_prf_id = max_prf_id + 1 ! basically, min_prf_id is unassigned...

    ! Find any atoms in this frame already assigned to a prf:

    do i = 1, atm_lst_cnt
      cur_prf_id = prf_map(atm_lst(i))
      if (cur_prf_id .ne. 0) then
        if (cur_prf_id .lt. min_prf_id) min_prf_id = cur_prf_id
      end if
    end do

    if (min_prf_id .gt. max_prf_id) then
      ! No frame atoms previously assigned to another prf; create a new one:
      new_prf_id = new_prf_id + 1
      do i = 1, atm_lst_cnt
        call add_atom(atm_lst(i), new_prf_id)
      end do
    else
      ! Assign this frame's atoms to the lowest prf already in use for atoms
      ! in the frame.
      do i = 1, atm_lst_cnt
        cur_prf_id = prf_map(atm_lst(i))
        if (cur_prf_id .eq. 0) then
          call add_atom(atm_lst(i), min_prf_id)
        else if (cur_prf_id .ne. min_prf_id) then
          call change_prf_id(cur_prf_id, min_prf_id)
        end if
      end do
    end if

  end do

  ! This only covers bonds to H...

  do bond_id = 1, bond_cnt

    call make_bond_atm_lst(bonds(bond_id), atm_lst_cnt, atm_lst)

    min_prf_id = max_prf_id + 1 ! basically, min_prf_id is unassigned...

    ! Find any atoms in this bond already assigned to a prf:

    do i = 1, atm_lst_cnt
      cur_prf_id = prf_map(atm_lst(i))
      if (cur_prf_id .ne. 0) then
        if (cur_prf_id .lt. min_prf_id) min_prf_id = cur_prf_id
      end if
    end do

    if (min_prf_id .gt. max_prf_id) then
      ! No bond atoms previously assigned to another prf; create a new one:
      new_prf_id = new_prf_id + 1
      do i = 1, atm_lst_cnt
        call add_atom(atm_lst(i), new_prf_id)
      end do
    else
      ! Assign this bond's atoms to the lowest prf already in use for atoms
      ! in the bond.
      do i = 1, atm_lst_cnt
        cur_prf_id = prf_map(atm_lst(i))
        if (cur_prf_id .eq. 0) then
          call add_atom(atm_lst(i), min_prf_id)
        else if (cur_prf_id .ne. min_prf_id) then
          call change_prf_id(cur_prf_id, min_prf_id)
        end if
      end do
    end if

  end do

  ! Map any remaining atoms in a residue to new prf's, placing atoms with
  ! contiguous id's in the same prf.

  do res_id = 1, res_cnt

    call make_unassigned_res_atm_lst(res_atms(res_id), &
                                     res_atms(res_id+1) - 1, &
                                     atm_lst_cnt, atm_lst)

    if (atm_lst_cnt .gt. 0) then
      new_prf_id = new_prf_id + 1
      call add_atom(atm_lst(1), new_prf_id)
      last_atm_id = atm_lst(1)
      do i = 2, atm_lst_cnt
        if (atm_lst(i) - last_atm_id .ne. 1) new_prf_id = new_prf_id + 1
        call add_atom(atm_lst(i), new_prf_id)
        last_atm_id = atm_lst(i)
      end do
    end if

  end do

  ! Count up the atoms in each prf

  prf_listdata(:)%offset = 0
  prf_listdata(:)%cnt = 0

  do atm_id = 1, atm_cnt
    prf_id = prf_map(atm_id)
#ifdef PRF_DBG
    if (prf_id .eq. 0) then
      write(0, *) 'DBG: found atom not assigned to prf!', atm_id
      call mexit(6, 1)
    end if
#endif /* PRF_DBG */
    prf_listdata(prf_id)%cnt = prf_listdata(prf_id)%cnt + 1
  end do

  ! Use prf counts to get offsets.  NOTE that the maskdata might be nondense;
  ! therefore there may be unused prf_id's (well, probably not with the given
  ! algorithm, but we will eliminate any prf's with 0 atoms below, just to be
  ! sure).

  offset = 0

  do prf_id =  1, atm_cnt
    prf_listdata(prf_id)%offset = offset
    offset = offset + prf_listdata(prf_id)%cnt
    prf_listdata(prf_id)%cnt = 0 ! will be used as ctr in list construction.
  end do

  ! Now make the mask and sum up the listcnt data:

  prf_lists(:) = 0

  do atm_id = 1, atm_cnt
    prf_id = prf_map(atm_id)
    prf_listdata(prf_id)%cnt = prf_listdata(prf_id)%cnt + 1
    prf_lists(prf_listdata(prf_id)%offset+prf_listdata(prf_id)%cnt) = &
      atm_id
  end do

  ! Find the highest used prf_id; we are not certain that there are no prf's
  ! with 0 atoms here, so we may slightly overallocate...

  do prf_id = atm_cnt, 1, -1
    if (prf_listdata(prf_id)%cnt .gt. 0) then
      max_prf_id = prf_id
      exit
    end if
  end do

  allocate(gbl_prf_listdata(max_prf_id), &
           gbl_prf_lists(atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_prf_listdata) * listdata_rec_ints + &
                        size(gbl_prf_lists)

  ! Store the prf data structures in globally saved storage, reordering the
  ! prf's such that the lowest atm_id in prf(n-1) .lt. lowest atm_id in prf(n).
  ! The way this is done, any blank prf entries (prf's with 0 atoms) are
  ! also skipped.  Also determine size of largest prf for later use.

  new_prf_id = 0
  new_offset = 0
  gbl_max_prf_atms = 0

  do atm_id = 1, atm_cnt
    prf_id = prf_map(atm_id)

    if (prf_id .ne. 0) then

      new_prf_id = new_prf_id + 1
      listcnt = prf_listdata(prf_id)%cnt
      if (gbl_max_prf_atms .lt. listcnt) gbl_max_prf_atms = listcnt
      old_offset = prf_listdata(prf_id)%offset
      
      gbl_prf_listdata(new_prf_id)%cnt = listcnt
      gbl_prf_listdata(new_prf_id)%offset = new_offset

      gbl_prf_lists(new_offset+1:new_offset+listcnt) = &
        prf_lists(old_offset+1:old_offset+listcnt)

      do i = old_offset + 1, old_offset + listcnt
        prf_map(prf_lists(i)) = 0
      end do

      new_offset = new_offset + listcnt

    end if
  end do

  gbl_prf_cnt = new_prf_id     ! store global count of prf's.

#ifdef PRF_DBG
  if (master) then
  
    ! Have to rewrite the prf_map():
 
    do prf_id = 1, gbl_prf_cnt
      offset = gbl_prf_listdata(prf_id)%offset
      listcnt = gbl_prf_listdata(prf_id)%cnt
      do i = offset + 1, offset + listcnt
        prf_map(gbl_prf_lists(i)) = prf_id
      end do
    end do
 
    write(0, *)'Atom identifier, Atom symbol, Pseudo-residue fragment no.:'
    do atm_id = 1, atm_cnt
      prf_id = prf_map(atm_id)
      write(0, '(i8,4x,a4,4x,i8)') atm_id, atm_igraph(atm_id), prf_id
    end do
  
    write(0, *)'gbl_prf_listdata() size =', size(gbl_prf_listdata)
    write(0, *)'gbl_prf_cnt =', gbl_prf_cnt
    write(0, *)'max_prf_atms =', gbl_max_prf_atms
  
    do prf_id = 1, gbl_prf_cnt
      offset = gbl_prf_listdata(prf_id)%offset
      listcnt = gbl_prf_listdata(prf_id)%cnt
      write(0, *) 'prf id, offset, atm lst cnt =', &
        prf_id, offset, listcnt
      write(0, *) 'atm lst =', gbl_prf_lists(offset+1:offset+listcnt)
    end do
  
  end if
#endif /* PRF_DBG */

  deallocate(atm_lst)

  return

contains

!*******************************************************************************
!
! Internal function:  get_max_res_atms
!
! Description: <TBS>
!              
!*******************************************************************************

integer &
function get_max_res_atms()

  use prmtop_dat_mod

  implicit none

! Local variables:

  integer       :: i
  integer       :: res_atms

  get_max_res_atms = 0

  do i = 1, nres
    res_atms = gbl_res_atms(i+1) - gbl_res_atms(i)
    if (res_atms .gt. get_max_res_atms) get_max_res_atms = res_atms
  end do

end function get_max_res_atms

!*******************************************************************************
!
! Internal subroutine:  make_frame_atm_lst
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine make_frame_atm_lst(frame, atm_lst_cnt, atm_lst)

  implicit none

! Formal arguments:

  type(ep_frame_rec), intent(in)        :: frame
  integer, intent(out)                  :: atm_lst_cnt
  integer, intent(out)                  :: atm_lst(*)

! Local variables:

  integer               :: ep_idx

  atm_lst_cnt = 0

  do ep_idx = 1, frame%ep_cnt
    atm_lst_cnt = atm_lst_cnt + 1
    atm_lst(atm_lst_cnt) = frame%extra_pnt(ep_idx)
  end do

  atm_lst_cnt = atm_lst_cnt + 1
  atm_lst(atm_lst_cnt) = frame%parent_atm
  atm_lst_cnt = atm_lst_cnt + 1
  atm_lst(atm_lst_cnt) = frame%frame_atm1
  atm_lst_cnt = atm_lst_cnt + 1
  atm_lst(atm_lst_cnt) = frame%frame_atm2
  atm_lst_cnt = atm_lst_cnt + 1
  atm_lst(atm_lst_cnt) = frame%frame_atm3

  return

end subroutine make_frame_atm_lst

!*******************************************************************************
!
! Internal subroutine:  make_bond_atm_lst
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine make_bond_atm_lst(bond, atm_lst_cnt, atm_lst)

  implicit none

! Formal arguments:

  type(bond_rec), intent(in)    :: bond
  integer, intent(out)          :: atm_lst_cnt
  integer, intent(out)          :: atm_lst(*)

! Local variables:

  atm_lst(1) = bond%atm_i
  atm_lst(2) = bond%atm_j
  atm_lst_cnt = 2

  return

end subroutine make_bond_atm_lst

!*******************************************************************************
!
! Internal subroutine:  make_unassigned_res_atm_lst
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine make_unassigned_res_atm_lst(first_res_atm, last_res_atm, &
                                       atm_lst_cnt, atm_lst)

  implicit none

! Formal arguments:

  integer, intent(in)           :: first_res_atm
  integer, intent(in)           :: last_res_atm
  integer, intent(out)          :: atm_lst(*)
  integer, intent(out)          :: atm_lst_cnt

! Local variables:

  integer                       :: atm_id

  atm_lst_cnt = 0

  do atm_id = first_res_atm, last_res_atm
   if (prf_map(atm_id) .eq. 0) then
     atm_lst_cnt = atm_lst_cnt + 1
     atm_lst(atm_lst_cnt) = atm_id
   end if
  end do

  return

end subroutine make_unassigned_res_atm_lst

!*******************************************************************************
!
! Internal subroutine:  add_atom
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine add_atom(atm_id, prf_id)

  implicit none

! Formal arguments:

  integer, intent(in)   :: atm_id
  integer, intent(in)   :: prf_id

! Local variables:

  integer               :: last_atm_id

  ! Either start new prf list or add to existing prf list:

  if (prf_lst_head(prf_id) .eq. 0) then
    ! Starting a new prf list:
    prf_lst_head(prf_id) = atm_id
    prf_lst_tail(prf_id) = atm_id
  else
    ! Appending to an existing prf list:
    last_atm_id = prf_lst_tail(prf_id)
    prf_lst(last_atm_id) = atm_id
    prf_lst_tail(prf_id) = atm_id
  end if


  prf_lst(atm_id) = 0
  prf_map(atm_id) = prf_id

  return

end subroutine add_atom

!*******************************************************************************
!
! Internal subroutine:  change_prf_id
!
! Description: Move atoms from one prf_id to another.  This is done by moving
!              the list.  It is ASSUMED that the old list is not null!
!              
!*******************************************************************************

subroutine change_prf_id(old_prf_id, new_prf_id)

  implicit none

! Formal arguments:

  integer, intent(in)   :: old_prf_id
  integer, intent(in)   :: new_prf_id

! Local variables:

  integer       :: atm_id
  integer       :: last_atm_id

  ! Change the pref_id's in the prf_map:

  atm_id = prf_lst_head(old_prf_id)
  do while (atm_id .ne. 0)
    prf_map(atm_id) = new_prf_id
    atm_id = prf_lst(atm_id)
  end do

  ! Move the list head and tail.

  if (prf_lst_head(new_prf_id) .eq. 0) then

    ! The list is being moved to a currently null list:
    prf_lst_head(new_prf_id) = prf_lst_head(old_prf_id)
    prf_lst_tail(new_prf_id) = prf_lst_tail(old_prf_id)

  else

    ! The list is being appended to a list with entries:
    last_atm_id = prf_lst_tail(new_prf_id)
    prf_lst(last_atm_id) = prf_lst_head(old_prf_id)
    prf_lst_tail(new_prf_id) = last_atm_id

  end if

  prf_lst_head(old_prf_id) = 0
  prf_lst_tail(old_prf_id) = 0

  return

end subroutine change_prf_id

end subroutine setup_prfs

#endif /* MPI */
end module prfs_mod
