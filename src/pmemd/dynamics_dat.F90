#include "copyright.i"

!*******************************************************************************
!
! Module:  dynamics_dat_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module dynamics_dat_mod

  implicit none

! Global data definitions.

  ! The following storage is per-process common; ie., it SHOULD be
  ! broadcast from the master to the other processes!

  integer, parameter    :: dynamics_dat_dbl_cnt = 1

  ! tmass = total mass of atoms in system.
  
  double precision              tmass

  common / dynamics_dat_dbl /   tmass

  save  :: / dynamics_dat_dbl /

#ifdef MPI

  ! The following variables are not broadcast:
  !
  ! max_unfrag_mol_atms = maximum atoms allowed for a molecule that will be
  !                       handled by 1 processor without fragmentation.
  ! frag_mol_cnt = count of molecules that will be fragmented to enhance load
  !                balancing for constant pressure md.
  ! molfrags_cnt = the number of molecule fragments across all fragmented
  !                molecules.

  integer, save         :: max_unfrag_mol_atms = 0      ! catch bugs...
  integer, save         :: frag_mol_cnt = 0
  integer, save         :: molfrags_cnt = 0

  ! Type describes a fragmented molecule:

  type frag_mol_rec
    sequence
    integer     :: mol_idx           ! Points back into gbl_mol_atms_listdata().
    integer     :: first_molfrag_idx ! Where molecule fragment records start
    integer     :: molfrag_cnt       ! Ends up with exact value.
                                     ! in gbl_molfrags()
    integer     :: group             ! for mpi
    integer     :: communicator      ! for mpi
    integer     :: task_cnt          ! for mpi
  end type frag_mol_rec

  integer, parameter    :: frag_mol_rec_ints = 6 ! don't use for allocation!

  ! Type describes a molecule fragment:

  type molfrag_rec
    sequence
    integer     :: offset               ! into gbl_mol_list
    integer     :: cnt                  ! for this fragment.
    integer     :: owner                ! task id; -1 if unassigned.
  end type molfrag_rec

  integer, parameter    :: molfrag_rec_ints = 3 ! don't use for allocation!

  ! Task-specific fragmented molecule counts:

  integer, save         :: my_frag_mol_cnt = 0

#endif /* MPI */

  ! Task-specific molecule counts:

  integer, save         :: my_mol_cnt = 0

  ! atm_rel_crd = the atom xyz coordinates relative to the center of mass
  ! atm_mass_inv = the inverted atom mass array.
  ! gbl_mol_mass_inv = the inverted molecule mass array.
  ! gbl_mol_com = center of mass coordinates for molecules.
  ! gbl_my_mol_lst = indexes into gbl_mol_atms_listdata() for molecules
  !                  "owned" by this task.
  ! gbl_frag_mols = description of molecules that must be fragmented, plus
  !                 indexes into molecule fragment descriptors (gbl_frag_mols())
  ! gbl_molfrags = the fragmented molecule fragment descriptors.
  ! gbl_my_frag_mol_lst = list of indexes into gbl_frag_mols() for fragmented
  !                       molecules, one or more fragments of which is handled
  !                       by this task (it should normally be only one, but
  !                       more may be possible due to overflow algorithms).

  double precision,     allocatable, save       :: atm_rel_crd(:,:)
  double precision,     allocatable, save       :: atm_mass_inv(:)
  double precision,     allocatable, save       :: gbl_mol_mass_inv(:)
  double precision,     allocatable, save       :: gbl_mol_com(:,:)

  ! Data below this point is not broadcast:
  !
  ! The following molecule-related array is only available for CP MD
  integer,              allocatable, save       :: gbl_my_mol_lst(:)
#ifdef MPI
  type(frag_mol_rec),   allocatable, save       :: gbl_frag_mols(:)
  type(molfrag_rec),    allocatable, save       :: gbl_molfrags(:)
  integer,              allocatable, save       :: gbl_my_frag_mol_lst(:)
#endif /* MPI */

! Hide internal routines:

  private       alloc_dynamics_mem

contains

!*******************************************************************************
!
! Subroutine:  dynamics_dat_init
!
! Description: <TBS>
!
!*******************************************************************************

subroutine dynamics_dat_init(atm_cnt, ntp, imin, mass, crd, num_ints, num_reals)

  use file_io_dat_mod
  use gbl_constants_mod
  use mol_list_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in)           :: atm_cnt
  integer, intent(in)           :: ntp
  integer, intent(in)           :: imin
  double precision, intent(in)  :: mass(atm_cnt)          ! atom mass array
  double precision, intent(in)  :: crd(3, atm_cnt)        ! atom crd array
  integer, intent(in out)       :: num_ints, num_reals     

! Local variables:

  integer                       :: i                    ! mol atm list idx
  integer                       :: atm_id, mol_id
  integer                       :: mol_atm_cnt, offset
  double precision              :: mol_mass_inv(gbl_mol_cnt)    ! temp array

  ! Allocate base memory needed for dynamics.

  call alloc_dynamics_mem(atm_cnt, gbl_mol_cnt, ntp, imin, num_ints, num_reals)

  ! Calculate inverse, total masses:

  tmass = 0.d0

  if (ntp .gt. 0 .and. imin .eq. 0) then

    do mol_id = 1, gbl_mol_cnt
      mol_mass_inv(mol_id) = 0.d0
      mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
      offset = gbl_mol_atms_listdata(mol_id)%offset
      do i = offset + 1, offset + mol_atm_cnt
        atm_id = gbl_mol_atms_lists(i)
        mol_mass_inv(mol_id) = mol_mass_inv(mol_id) + mass(atm_id) ! Sum mols.
        if (mass(atm_id) .ne. 0.d0) then
          atm_mass_inv(atm_id) = 1.d0 / mass(atm_id) ! Save inverted mass.
        else
          atm_mass_inv(atm_id) = 0.d0 ! For extra points!
        end if
      end do
      tmass = tmass + mol_mass_inv(mol_id)
      mol_mass_inv(mol_id) = 1.d0 / mol_mass_inv(mol_id) ! Actually invert here.
    end do

    gbl_mol_mass_inv(:) = mol_mass_inv(:)

  ! Initialize molecule COM for all atoms if needed:

    call get_all_mol_com(crd, mass, gbl_mol_mass_inv, gbl_mol_com)

  else

    do atm_id = 1, atm_cnt
      tmass = tmass + mass(atm_id)
      if (mass(atm_id) .ne. 0.d0) then
        atm_mass_inv(atm_id) = 1.d0 / mass(atm_id) ! Save inverted mass.
      else
        atm_mass_inv(atm_id) = 0.d0 ! For extra points!
      end if
    end do

  end if

  return

end subroutine dynamics_dat_init

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  setup_fragmented_molecules
!
! Description: Only used to fragment molecules if (ntp .gt. 0 .and. imin .eq. 0)
!
!*******************************************************************************

subroutine setup_fragmented_molecules(atm_cnt, num_ints, num_reals)

  use file_io_dat_mod
  use gbl_constants_mod
  use parallel_dat_mod
  use pmemd_lib_mod
  use prfs_mod
  use mol_list_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in)           :: atm_cnt
  integer, intent(in out)       :: num_ints, num_reals     

! Local variables:

  integer                       :: alloc_failed
  integer                       :: molfrag_cnt
  integer                       :: gfm_idx
  integer                       :: gmf_idx
  integer                       :: i, j
  integer                       :: mol_id
  integer                       :: mol_atm_cnt
  integer                       :: mol_prf_cnt
  integer                       :: mol_prf_offset
  integer                       :: prf_cnt
  integer                       :: prf_cnt_total
  double precision              :: prfs_per_molfrag

  ! Note that this code is executed by all processes.

  ! Determine how much, if any memory is needed for fragmented molecule
  ! support (only for CP MD under MPI).


  ! The maximum unfragmented molecule atoms value is used to determine
  ! which atoms we fragment.

  max_unfrag_mol_atms = atm_cnt / numtasks + 1

  if (max_unfrag_mol_atms .lt. gbl_max_prf_atms * 2) &
    max_unfrag_mol_atms = gbl_max_prf_atms * 2

! BEGIN DBG
! max_unfrag_mol_atms = atm_cnt / 96 + 1        ! for debugging on small setup
! END DBG

  frag_mol_cnt = 0
  molfrags_cnt = 0

  do mol_id = 1, gbl_mol_cnt
    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    if (mol_atm_cnt .gt. max_unfrag_mol_atms) then
      frag_mol_cnt = frag_mol_cnt + 1
      molfrags_cnt = molfrags_cnt + mol_atm_cnt / max_unfrag_mol_atms
      if (mod(mol_atm_cnt, max_unfrag_mol_atms) .ne. 0) then
        molfrags_cnt = molfrags_cnt + 1
      end if
    end if
  end do

  if (frag_mol_cnt .gt. 0) then

    allocate(gbl_frag_mols(frag_mol_cnt), &
             gbl_molfrags(molfrags_cnt), &
             gbl_my_frag_mol_lst(frag_mol_cnt), &
             stat = alloc_failed)
               
    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(gbl_frag_mols) * frag_mol_rec_ints + &
                          size(gbl_molfrags) * molfrag_rec_ints + &
                          size(gbl_my_frag_mol_lst)

    ! Okay, now actually find the molecules you have to fragment and fragment
    ! them, storing information in the gbl_frag_mol and gbl_molfrags arrays.
    ! We do the actual fragmentation based on prf's. This may not produce
    ! perfect fragments sizes, but is certainly less error prone, and it
    ! is not going to really matter much...

    gbl_molfrags(1:molfrags_cnt)%offset = 0
    gbl_molfrags(1:molfrags_cnt)%cnt = 0
    gbl_molfrags(1:molfrags_cnt)%owner = -1

    gfm_idx = 1    ! Index into fragmented molecules.
    gmf_idx = 1    ! Index into molecule fragments.

    do mol_id = 1, gbl_mol_cnt

      mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt

      ! Find molecules that need to be fragmented:

      if (mol_atm_cnt .gt. max_unfrag_mol_atms) then

        molfrag_cnt = mol_atm_cnt / max_unfrag_mol_atms
        if (mod(mol_atm_cnt, max_unfrag_mol_atms) .ne. 0) then
          molfrag_cnt = molfrag_cnt + 1
        end if

        gbl_frag_mols(gfm_idx)%mol_idx = mol_id
        gbl_frag_mols(gfm_idx)%molfrag_cnt = molfrag_cnt
        gbl_frag_mols(gfm_idx)%first_molfrag_idx = gmf_idx
        gbl_frag_mols(gfm_idx)%group = MPI_GROUP_NULL
        gbl_frag_mols(gfm_idx)%communicator = MPI_COMM_NULL
        gbl_frag_mols(gfm_idx)%task_cnt = 0

        ! We now create exactly molfrag_cnt fragments on prf boundaries.

        mol_prf_offset = gbl_mol_prfs_listdata(mol_id)%offset
        mol_prf_cnt = gbl_mol_prfs_listdata(mol_id)%cnt

        prfs_per_molfrag = dble(mol_prf_cnt)/dble(molfrag_cnt)
        prf_cnt_total = 0

        do j = 1, molfrag_cnt
          gbl_molfrags(gmf_idx)%offset = mol_prf_offset + prf_cnt_total
          prf_cnt = int(prfs_per_molfrag * dble(j)) - prf_cnt_total
          gbl_molfrags(gmf_idx)%cnt = prf_cnt
          prf_cnt_total = prf_cnt_total + prf_cnt
          gmf_idx = gmf_idx + 1
        end do

        if (prf_cnt_total .ne. mol_prf_cnt) then
          gbl_molfrags(gmf_idx-1)%cnt = &
            gbl_molfrags(gmf_idx-1)%cnt + (mol_prf_cnt - prf_cnt_total)
! BEGIN DBG
          if (master) write(0,*)'DBG:corrected last cnt =', &
            gbl_molfrags(gmf_idx-1)%cnt
! END DBG
        end if

        gfm_idx = gfm_idx + 1

      end if

    end do

  end if

! BEGIN DBG
! if (master) then
!   write(0,*)'DBG: max_unfrag_mol_atms =', max_unfrag_mol_atms
!   write(0,*)'DBG: frag_mol_cnt =', frag_mol_cnt
!   write(0,*)'DBG: molfrags_cnt =', molfrags_cnt
!
!   do i = 1, frag_mol_cnt
!
!     write(0,*)'DBG: frag mol', i, ' describes mol', gbl_frag_mols(i)%mol_idx
!     write(0,*)'DBG: frag cnt =', gbl_frag_mols(i)%molfrag_cnt
!     write(0,*)'DBG: first frag idx =', gbl_frag_mols(i)%first_molfrag_idx
!
!     do j = gbl_frag_mols(i)%first_molfrag_idx, &
!            gbl_frag_mols(i)%first_molfrag_idx + &
!            gbl_frag_mols(i)%molfrag_cnt - 1
!
!       write(0,2000)' DBG: frag, offset, prf cnt =', j, &
!                    gbl_molfrags(j)%offset, &
!                    gbl_molfrags(j)%cnt
!     end do
!
!2000 format(a, 3i7)
!
!   end do
! end if
! END DBG

  return

end subroutine setup_fragmented_molecules
#endif /* MPI */

!*******************************************************************************
!
! Subroutine:  alloc_dynamics_mem
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine alloc_dynamics_mem(atm_cnt, mol_cnt, ntp, imin, num_ints, num_reals)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: atm_cnt
  integer, intent(in)           :: mol_cnt
  integer, intent(in)           :: ntp
  integer, intent(in)           :: imin

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     

! Local variables:

  integer                       :: alloc_failed

  ! Generic allocations required for just about everything:

  allocate(atm_mass_inv(atm_cnt), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(atm_mass_inv)

  ! Relative coordinates should only be needed for constant pressure MD:

  if (ntp .gt. 0 .and. imin .eq. 0) then

    allocate(atm_rel_crd(3, atm_cnt), &
#ifdef MPI
             gbl_my_mol_lst(mol_cnt), &
#endif
             gbl_mol_mass_inv(mol_cnt), &
             gbl_mol_com(3, mol_cnt), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

#ifdef MPI
    num_ints = num_ints + size(gbl_my_mol_lst)
#endif
    num_reals = num_reals + size(atm_rel_crd) + &
                            size(gbl_mol_mass_inv) + &
                            size(gbl_mol_com)
  end if

  ! No need to initialize atm_rel_crd.
  ! Other stuff initialized in other code.

  return

end subroutine alloc_dynamics_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  create_communicators
!
! Description: <TBS>
!
!*******************************************************************************

subroutine create_communicators

  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  ! None

! Local variables:

  integer       :: tasklist(numtasks)           ! final tasklist
  integer       :: task_cnt
  integer       :: mol_idx
  integer       :: frag_cnt
  integer       :: first_frag_idx
  integer       :: frag_idx
  integer       :: task_idx
  integer       :: taskid
  logical       :: nodup_found

  do mol_idx = 1, frag_mol_cnt

    frag_cnt = gbl_frag_mols(mol_idx)%molfrag_cnt
    first_frag_idx = gbl_frag_mols(mol_idx)%first_molfrag_idx

    tasklist(1) = gbl_molfrags(first_frag_idx)%owner
    task_cnt = 1

    ! Make a task list without any duplicates.  The algorithm is O(n**2), but
    ! n is typically pretty small, say 2-30.  All this hooplah is necessary
    ! because the target atom assignment may be larger or smaller than the
    ! default fragment size, which means that multiple atoms could get assigned
    ! to the same task.  Aside from having to do this duplicate search, that
    ! is not a bad thing.  The frequency of execution of this code should taper
    ! off quickly after initial loadbalancing.

    do frag_idx = first_frag_idx + 1, first_frag_idx + frag_cnt - 1
      taskid = gbl_molfrags(frag_idx)%owner
      nodup_found = .true.
      do task_idx = 1, task_cnt
        if (taskid .eq. tasklist(task_idx)) then
          nodup_found = .false.
          exit
        end if
      end do
      if (nodup_found) then
        task_cnt = task_cnt + 1
        tasklist(task_cnt) = taskid
      end if
    end do

    ! The group created below will have the task owning the first fragment
    ! of the molecule as the group master (task 0).  Very nice!

    if (task_cnt .gt. 1) then

      call mpi_group_incl(pmemd_group, task_cnt, tasklist, &
                          gbl_frag_mols(mol_idx)%group, err_code_mpi)

      call mpi_comm_create(pmemd_comm, gbl_frag_mols(mol_idx)%group, &
                           gbl_frag_mols(mol_idx)%communicator, err_code_mpi)

    end if
    gbl_frag_mols(mol_idx)%task_cnt = task_cnt

  end do

  return

end subroutine create_communicators
#endif /* MPI */

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  destroy_communicators
!
! Description:  Destroy mpi communications objects associated with fragmented
!               molecule processing.  Note that there is currently no mechanism
!               for calling this routine on either successful program exit or
!               error program exit; this is a circular dependency issue, and
!               we would have to factor things a bit differently.  However, I
!               have no reason to believe this is an issue; the mpi allocated
!               objects should be local to each task, and any mpi implementation
!               that does not protect itself from a failure to clean up allocs
!               is just plain broken.
!
!*******************************************************************************

subroutine destroy_communicators

  use parallel_dat_mod

  implicit none

! Formal arguments:

  ! None

! Local variables:

  integer       :: mol_idx

  do mol_idx = 1, frag_mol_cnt
    if (gbl_frag_mols(mol_idx)%communicator .ne. MPI_COMM_NULL) then
      call mpi_comm_free(gbl_frag_mols(mol_idx)%communicator, err_code_mpi)
      gbl_frag_mols(mol_idx)%communicator = MPI_COMM_NULL
    end if
    if (gbl_frag_mols(mol_idx)%group .ne. MPI_GROUP_NULL) then
      call mpi_group_free(gbl_frag_mols(mol_idx)%group, err_code_mpi)
      gbl_frag_mols(mol_idx)%group = MPI_GROUP_NULL
    end if
  end do

  return

end subroutine destroy_communicators
#endif /* MPI */

!*******************************************************************************
!
! Subroutine:   get_all_mol_com
!
! Description:  Pressure scaling routine for crds. ONLY used for constant
!               pressure scaling (ntp .gt. 0).  This version is used when all
!               coordinates are known, even in an mpi context.  It is a little
!               wasteful, but very infrequently used.
!
!*******************************************************************************

subroutine get_all_mol_com(crd, mass, mol_mass_inv, mol_com)

  use mol_list_mod

  implicit none

! Formal arguments:

  double precision              :: crd(3, *)            ! atom crd array
  double precision              :: mass(*)              ! atom mass array
  double precision              :: mol_mass_inv(*)
  double precision              :: mol_com(3, *)

! Local variables:

  double precision              :: com(3)
  integer                       :: i                    ! mol atm list idx
  integer                       :: atm_id, mol_id
  integer                       :: mol_atm_cnt, offset

! Get COM for all molecules:

  do mol_id = 1, gbl_mol_cnt
    com(:) = 0.d0
    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      com(:) = com(:) + mass(atm_id) * crd(:, atm_id)
    end do
    mol_com(:, mol_id) = com(:) * mol_mass_inv(mol_id)
  end do

  return

end subroutine get_all_mol_com

end module dynamics_dat_mod
