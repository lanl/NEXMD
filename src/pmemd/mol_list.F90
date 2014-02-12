#include "copyright.i"

!*******************************************************************************
!
! Module: mol_list_mod
!
! Description: Support for molecule lists.  This routine handles molecule
!              merging if no_intermolecular_bonds .eq. 1, and is used to
!              construct lists of atoms in molecules and prf's in molecules.
!              These larger more complex data structures are needed due to the
!              implications of molecule merging (to avoid intermolecular
!              bonds) and frames (which require the pseudo-residue fragment
!              based ownership scheme under MPI).
!*******************************************************************************

module mol_list_mod

  use gbl_datatypes_mod

  implicit none

  ! Pseudo-residue fragment data, used for PME and GB:

  integer, save :: gbl_mol_cnt = 0

  type(listdata_rec), allocatable, save :: gbl_mol_atms_listdata(:)
  integer, allocatable, save            :: gbl_mol_atms_lists(:)

#ifdef MPI
  type(listdata_rec), allocatable, save :: gbl_mol_prfs_listdata(:)
  integer, allocatable, save            :: gbl_mol_prfs_lists(:)
#endif /* MPI */

contains

!*******************************************************************************
!
! Subroutine:  setup_molecule_lists
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine setup_molecule_lists(num_ints, num_reals, atm_cnt, mol_cnt_in, &
                                mol_atm_cnts, bond_cnt, bonds, &
                                no_intermolecular_bonds)

  use gbl_datatypes_mod
  use parallel_dat_mod
  use pmemd_lib_mod
#ifdef DBG_MOL_LIST
  use prmtop_dat_mod
#endif /* DBG_MOL_LIST */
#ifdef CUDA
  use constraints_mod
#endif

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)               :: num_ints, num_reals
  integer, intent(in)                   :: atm_cnt
  integer, intent(in)                   :: mol_cnt_in
  integer, intent(in)                   :: mol_atm_cnts(mol_cnt_in)
  integer, intent(in)                   :: bond_cnt
  type(bond_rec), intent(in)            :: bonds(*)
  integer, intent(in)                   :: no_intermolecular_bonds

! Local variables:

  integer               :: mol_id
  integer               :: atm_id
  integer               :: first_atm_id, last_atm_id
  integer               :: bond_id
  integer               :: mol_id_i, mol_id_j
  integer               :: offset
  integer               :: listcnt
  integer               :: max_mol_id
  integer               :: alloc_failed
  integer               :: i
  integer               :: new_mol_id
  integer               :: old_offset, new_offset
  integer               :: mol_lst_head(mol_cnt_in)
  integer               :: mol_lst_tail(mol_cnt_in)
  integer               :: mol_lst(atm_cnt)
  integer               :: mol_map(atm_cnt)
  type(listdata_rec)    :: mol_listdata(mol_cnt_in)
  integer               :: mol_lists(atm_cnt)

! Pre-initialize as needed:

  mol_lst_head(:) = 0
  mol_lst_tail(:) = 0
  mol_lst(:) = 0        ! Marks empty mol list entries.

  ! Create the molecule map based on current knowledge.  It is not strictly
  ! needed if we are not merging molecules with an intermolecular bond, but
  ! some following code is simpler, possibly more efficient...

  last_atm_id = 0

  do mol_id = 1, mol_cnt_in
    first_atm_id = last_atm_id + 1
    last_atm_id = first_atm_id + mol_atm_cnts(mol_id) - 1
    mol_map(first_atm_id:last_atm_id) = mol_id
  end do

  ! Now create the molecule lists based on current knowledge...

  do atm_id = 1, atm_cnt
    call add_atom(atm_id, mol_map(atm_id))
  end do

  ! Now look for intermolecular bonds and merge molecules...  In theory, we
  ! should be able to ignore bonds containing H here; just to be sure we
  ! look at them all.

  if (no_intermolecular_bonds .ne. 0) then
    do bond_id = 1, bond_cnt
      mol_id_i = mol_map(bonds(bond_id)%atm_i)
      mol_id_j = mol_map(bonds(bond_id)%atm_j)
      if (mol_id_i .ne. mol_id_j) then
        if (mol_id_i .lt. mol_id_j) then
          call change_mol_id(mol_id_j, mol_id_i)
        else
          call change_mol_id(mol_id_i, mol_id_j)
        end if
      end if
    end do
  end if

  ! Count up the atoms in each molecule

  mol_listdata(:)%offset = 0
  mol_listdata(:)%cnt = 0

  do atm_id = 1, atm_cnt
    mol_id = mol_map(atm_id)
    mol_listdata(mol_id)%cnt = mol_listdata(mol_id)%cnt + 1
  end do

  ! Use mol counts to get offsets.  NOTE that the maskdata might be nondense;
  ! therefore there may be unused mol_id's (well, probably not with the given
  ! algorithm, but we will eliminate any mol's with 0 atoms below, just to be
  ! sure).

  offset = 0

  do mol_id =  1, mol_cnt_in
    mol_listdata(mol_id)%offset = offset
    offset = offset + mol_listdata(mol_id)%cnt
    mol_listdata(mol_id)%cnt = 0 ! will be used as ctr in list construction.
  end do

  ! Now make the mask and sum up the listcnt data:

  mol_lists(:) = 0

  do atm_id = 1, atm_cnt
    mol_id = mol_map(atm_id)
    mol_listdata(mol_id)%cnt = mol_listdata(mol_id)%cnt + 1
    mol_lists(mol_listdata(mol_id)%offset+mol_listdata(mol_id)%cnt) = atm_id
  end do

  ! Find the highest used mol_id; we are not certain that there are no
  ! molecules with 0 atoms here, so we may slightly overallocate...

  max_mol_id = 0

  do mol_id = mol_cnt_in, 1, -1
    if (mol_listdata(mol_id)%cnt .gt. 0) then
      max_mol_id = mol_id
      exit
    end if
  end do

  allocate(gbl_mol_atms_listdata(max_mol_id), &
           gbl_mol_atms_lists(atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_mol_atms_listdata) * listdata_rec_ints + &
                        size(gbl_mol_atms_lists)

  ! Store the mol data structures in globally saved storage, reordering the
  ! mol's such that the lowest atm_id in mol(n-1) .lt. lowest atm_id in mol(n).
  ! The way this is done, any blank mol entries (mol's with 0 atoms) are
  ! also skipped.

  new_mol_id = 0
  new_offset = 0

  do atm_id = 1, atm_cnt
    mol_id = mol_map(atm_id)

    if (mol_id .ne. 0) then

      new_mol_id = new_mol_id + 1
      listcnt = mol_listdata(mol_id)%cnt
      old_offset = mol_listdata(mol_id)%offset
      
      gbl_mol_atms_listdata(new_mol_id)%cnt = listcnt
      gbl_mol_atms_listdata(new_mol_id)%offset = new_offset

      gbl_mol_atms_lists(new_offset+1:new_offset+listcnt) = &
        mol_lists(old_offset+1:old_offset+listcnt)

      do i = old_offset + 1, old_offset + listcnt
        mol_map(mol_lists(i)) = 0
      end do

      new_offset = new_offset + listcnt

    end if
  end do

  gbl_mol_cnt = new_mol_id     ! store global count of mol's.
  
#ifdef CUDA
  call gpu_molecule_list_setup(gbl_mol_cnt, gbl_mol_atms_listdata)
  if (natc .gt. 0) then
    call gpu_constraint_molecule_list_setup(gbl_mol_cnt, gbl_mol_atms_listdata, atm_xc)
  end if
#endif

#ifdef DBG_MOL_LIST
  if (master) then
  
    ! Have to rewrite the mol_map():
 
    do mol_id = 1, gbl_mol_cnt
      offset = gbl_mol_atms_listdata(mol_id)%offset
      listcnt = gbl_mol_atms_listdata(mol_id)%cnt
      do i = offset + 1, offset + listcnt
        mol_map(gbl_mol_atms_lists(i)) = mol_id
      end do
    end do
 
    write(0, *)'Atom identifier, Atom symbol, mol id:'
    do atm_id = 1, atm_cnt
      mol_id = mol_map(atm_id)
      write(0, '(i8,4x,a4,4x,i8)') atm_id, atm_igraph(atm_id), mol_id
    end do
  
    write(0, *)'gbl_mol_atms_listdata() size =', size(gbl_mol_atms_listdata)
    write(0, *)'gbl_mol_cnt =', gbl_mol_cnt
  
    do mol_id = 1, gbl_mol_cnt
      offset = gbl_mol_atms_listdata(mol_id)%offset
      listcnt = gbl_mol_atms_listdata(mol_id)%cnt
      write(0, *) 'mol id, offset, atm lst cnt =', &
        mol_id, offset, listcnt
      write(0, *) 'atm lst =', gbl_mol_atms_lists(offset+1:offset+listcnt)
    end do
  
  end if
#endif /* DBG_MOL_LIST */

  return

contains

!*******************************************************************************
!
! Internal subroutine:  add_atom
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine add_atom(atm_id, mol_id)

  implicit none

! Formal arguments:

  integer, intent(in)   :: atm_id
  integer, intent(in)   :: mol_id

! Local variables:

  integer               :: last_atm_id

  ! Either start new molecule list or add to existing molecule list:

  if (mol_lst_head(mol_id) .eq. 0) then
    ! Starting a new molecule list:
    mol_lst_head(mol_id) = atm_id
    mol_lst_tail(mol_id) = atm_id
  else
    ! Appending to an existing molecule list:
    last_atm_id = mol_lst_tail(mol_id)
    mol_lst(last_atm_id) = atm_id
    mol_lst_tail(mol_id) = atm_id
  end if


  mol_lst(atm_id) = 0
  mol_map(atm_id) = mol_id

  return

end subroutine add_atom

!*******************************************************************************
!
! Internal subroutine:  change_mol_id
!
! Description: Move atoms from one mol_id to another.  This is done by moving
!              the list.  It is ASSUMED that the old list is not null!  Also
!              note that mol_map is updated!
!*******************************************************************************

subroutine change_mol_id(old_mol_id, new_mol_id)

  implicit none

! Formal arguments:

  integer, intent(in)   :: old_mol_id
  integer, intent(in)   :: new_mol_id

! Local variables:

  integer       :: atm_id
  integer       :: last_atm_id

  ! Change the pref_id's in the mol_map:

  atm_id = mol_lst_head(old_mol_id)
  do while (atm_id .ne. 0)
    mol_map(atm_id) = new_mol_id
    atm_id = mol_lst(atm_id)
  end do

  ! Move the list head and tail.

  if (mol_lst_head(new_mol_id) .eq. 0) then

    ! The list is being moved to a currently null list:
    mol_lst_head(new_mol_id) = mol_lst_head(old_mol_id)
    mol_lst_tail(new_mol_id) = mol_lst_tail(old_mol_id)

  else

    ! The list is being appended to a list with entries:
    last_atm_id = mol_lst_tail(new_mol_id)
    mol_lst(last_atm_id) = mol_lst_head(old_mol_id)
    mol_lst_tail(new_mol_id) = last_atm_id

  end if

  mol_lst_head(old_mol_id) = 0
  mol_lst_tail(old_mol_id) = 0

  return

end subroutine change_mol_id

end subroutine setup_molecule_lists

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  setup_mol_prf_data
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine setup_mol_prf_data(num_ints, num_reals, atm_cnt, prf_cnt, &
                              prf_listdata, prf_lists)

  use gbl_constants_mod
  use gbl_datatypes_mod
  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)               :: num_ints, num_reals
  integer, intent(in)                   :: atm_cnt
  integer, intent(in)                   :: prf_cnt
  type(listdata_rec), intent(in)        :: prf_listdata(prf_cnt)
  integer, intent(in)                   :: prf_lists(atm_cnt)

! Local variables:

  integer               :: alloc_failed
  integer               :: atm_id
  integer               :: i
  integer               :: first_mol_id
  integer               :: mol_id
  integer               :: listcnt
  integer               :: offset
  integer               :: prf_id
  integer               :: atm_mol_map(atm_cnt)
  type(listdata_rec)    :: mol_listdata(gbl_mol_cnt)
  integer               :: mol_lists(prf_cnt)
  integer               :: prf_molecules(prf_cnt)

  ! Make a atom to molecule map, to be used in finding the molecules in a prf.

  atm_mol_map(:) = -1   ! Mostly for debugging...

  do mol_id = 1, gbl_mol_cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    listcnt = gbl_mol_atms_listdata(mol_id)%cnt
    do i = offset + 1, offset + listcnt
      atm_mol_map(gbl_mol_atms_lists(i)) = mol_id
    end do
  end do

  mol_listdata(:)%offset = 0
  mol_listdata(:)%cnt = 0
  prf_molecules(:) = 0

  ! Map prf's to a molecule.  We do not allow intermolecular bonding if extra
  ! points are being used, so all the atoms in a prf should map to a single
  ! molecule.  We error exit if this is not the case.  This check may not
  ! be absolutely necessary, but it is cheap insurance against unexpected
  ! behaviour after setup.

  do prf_id = 1, prf_cnt

    listcnt = prf_listdata(prf_id)%cnt
    offset = prf_listdata(prf_id)%offset

    ! There should be no prf's with 0 atoms...

    first_mol_id = atm_mol_map(prf_lists(offset + 1))

    do i = offset + 2, offset + listcnt
      mol_id = atm_mol_map(prf_lists(i))
      if (mol_id .ne. first_mol_id) then
        write(mdout, '(a,a)') error_hdr, &
          'PMEMD does not support intermolecular PRFs!'
        call mexit(6, 1)
      end if
    end do

    mol_listdata(first_mol_id)%cnt = mol_listdata(first_mol_id)%cnt + 1
    prf_molecules(prf_id) = first_mol_id

  end do

  ! Use molecule counts (count of prf's in each molecule) to get offsets.

  offset = 0

  do mol_id = 1, gbl_mol_cnt
    mol_listdata(mol_id)%offset = offset
    offset = offset + mol_listdata(mol_id)%cnt
    mol_listdata(mol_id)%cnt = 0
  end do

  ! Now complete the molecule mask structures.

  mol_lists(:) = 0
  
  do prf_id = 1, prf_cnt

    mol_id = prf_molecules(prf_id)
    mol_listdata(mol_id)%cnt = mol_listdata(mol_id)%cnt + 1
    mol_lists(mol_listdata(mol_id)%offset+mol_listdata(mol_id)%cnt) = prf_id

  end do

  allocate(gbl_mol_prfs_listdata(gbl_mol_cnt), &
           gbl_mol_prfs_lists(prf_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_ints = num_ints + size(gbl_mol_prfs_listdata) * listdata_rec_ints + &
                        size(gbl_mol_prfs_lists)

  ! Store the mol mask data structures in globally saved storage.

  gbl_mol_prfs_listdata(:) = mol_listdata(:)
  gbl_mol_prfs_lists(1:prf_cnt) = mol_lists(1:prf_cnt)

#ifdef DBG_MOL_LIST
  if (master) then
 
    write(0, *) 'DBG: Dump of molecule prfs:'
     
    do mol_id = 1, mol_cnt
      write(0, *) 'Molecule id =', mol_id
      offset = gbl_mol_prfs_listdata(mol_id)%offset
      listcnt = gbl_mol_prfs_listdata(mol_id)%cnt
      write(0, *) '  offset, prf cnt =', offset, listcnt
      write(0, *) '  prf lst =', gbl_mol_prfs_lists(offset+1:offset+listcnt)
    end do
 
  end if
#endif /* DBG_MOL_LIST */

  return

end subroutine setup_mol_prf_data
#endif /* MPI */

end module mol_list_mod
