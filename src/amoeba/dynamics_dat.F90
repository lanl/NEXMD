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

  integer, parameter    :: dynamics_dat_int_cnt = 3

  ! max_unfrag_mol_len = maximum length allowed for a molecule that will be
  !                      handled by 1 processor without fragmentation.
  ! frag_mol_cnt = count of molecules that will be fragmented to enhance load
  !                balancing for constant pressure md.
  ! max_mol_frag_cnt = the maximum number of molecule fragments possible; this
  !                    is an estimate, and the actual number of molecule
  !                    fragments may be lower because the molecules must be
  !                    fragmented on residue boundaries.

  integer                       max_unfrag_mol_len, frag_mol_cnt, &
                                max_mol_frag_cnt

  common / dynamics_dat_int /   max_unfrag_mol_len, frag_mol_cnt, &
                                max_mol_frag_cnt

  save  :: / dynamics_dat_int /

  ! Type describes a fragmented molecule:

  type frag_mol_rec
    sequence
    integer     :: mol_idx              ! Points back into gbl_mol_atms().
    integer     :: frag_cnt             ! Ends up with exact value.
    integer     :: first_frag_idx       ! Where molecule fragment records start
                                        ! in gbl_mol_frags()
    integer     :: group                ! for mpi
    integer     :: communicator         ! for mpi
    integer     :: task_cnt             ! for mpi
  end type frag_mol_rec

  integer, parameter    :: frag_mol_rec_ints = 3 ! don't use for allocation!

  ! Type describes a molecule fragment:

  type mol_frag_rec
    sequence
    integer     :: first_atm_id         ! for this fragment.
    integer     :: atm_cnt              ! for this fragment.
    integer     :: owner                ! task id; -1 if unassigned.
  end type mol_frag_rec

  integer, parameter    :: mol_frag_rec_ints = 3 ! don't use for allocation!

#endif /* MPI */

  ! Task-specific molecule counts:

  integer, save         :: my_mol_cnt = 0
  integer, save         :: my_frag_mol_cnt = 0

  ! atm_rel_crd = the atom xyz coordinates relative to the center of mass
  ! atm_mass_inv = the inverted atom mass array.
  ! gbl_mol_mass_inv = the inverted molecule mass array.
  ! gbl_mol_atms = analogous to gbl_res_atms(), a list of first atm id's
  !                for all the molecules, capped by natom + 1.
  ! gbl_mol_com = center of mass coordinates for molecules.
  ! gbl_my_mol_lst = indexes into gbl_mol_atms() for molecules "owned" by
  !                  this task.
  ! gbl_frag_mols = description of molecules that must be fragmented, plus
  !                 indexes into molecule fragment descriptors (gbl_frag_mols())
  ! gbl_mol_frags = the fragmented molecule fragment descriptors.
  ! gbl_my_frag_mol_lst = list of indexes into gbl_frag_mols() for fragmented
  !                       molecules, one or more fragments of which is handled
  !                       by this task (it should normally be only one, but
  !                       more may be possible due to overflow algorithms).

  double precision,     allocatable, save       :: atm_rel_crd(:,:)
  double precision,     allocatable, save       :: atm_mass_inv(:)
  integer,              allocatable, save       :: gbl_mol_atms(:)
  double precision,     allocatable, save       :: gbl_mol_mass_inv(:)
  double precision,     allocatable, save       :: gbl_mol_com(:,:)
  ! The following molecule-related array is only available for CP MD:
  integer,              allocatable, save       :: gbl_my_mol_lst(:)
#ifdef MPI
  type(frag_mol_rec),   allocatable, save       :: gbl_frag_mols(:)
  type(mol_frag_rec),   allocatable, save       :: gbl_mol_frags(:)
  integer,              allocatable, save       :: gbl_my_frag_mol_lst(:)
#endif /* MPI */

! Hide internal routines:

  private       alloc_dynamics_mem

contains

!*******************************************************************************
!
! Subroutine:  init_dynamics_dat
!
! Description: <TBS>
!
!*******************************************************************************

#ifdef AMOEBA
subroutine init_dynamics_dat(natom, nres, nspm, ntp, imin, iamoeba, res_atms, &
                             atm_nsp, mass, crd, num_ints, num_reals)
#else
subroutine init_dynamics_dat(natom, nres, nspm, imin, ntp, res_atms, atm_nsp, &
                             mass, crd, num_ints, num_reals)
#endif /* AMOEBA */

  use parallel_dat_mod

  implicit none

! Formal arguments:

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in)           :: natom
  integer, intent(in)           :: nres
  integer, intent(in)           :: nspm
  integer, intent(in)           :: imin
  integer, intent(in)           :: ntp
#ifdef AMOEBA
  integer, intent(in)           :: iamoeba
#endif /* AMOEBA */
  integer, intent(in)           :: res_atms(nres + 1)   ! res atms list
  integer, intent(in)           :: atm_nsp(nspm)
  double precision, intent(in)  :: mass(natom)          ! atom mass array
  double precision, intent(in)  :: crd(3, natom)        ! atom crd array
  integer, intent(in out)       :: num_ints, num_reals     

! Local variables:

  integer                       :: atm_idx, last_atm_idx
  integer                       :: first_atm, last_atm, atm_cnt
  integer                       :: first_res, last_res
  integer                       :: last_res_first_atm
  integer                       :: last_res_atm_cnt
  integer                       :: last_res_atm_cnt_incl
  integer                       :: res_idx
  integer                       :: i, j, n
  integer                       :: mfrag_cnt            ! per molecule
  integer                       :: mfrag_idx
  integer                       :: cur_mfrag_idx
  integer                       :: last_mfrag_idx
  integer                       :: mol_atm_cnt

  double precision              :: com(3)
  double precision              :: mol_mass_inv(nspm)       ! temp array
  integer                       :: atm_res_map(natom)

  ! Note that this code is executed by the master.

#ifdef MPI

  ! Determine how much, if any memory is needed for fragmented molecule
  ! support (only for CP MD under MPI).

#ifdef AMOEBA
  if (ntp .gt. 0 .and. imin .eq. 0 .and. iamoeba .eq. 0) then
#else
  if (ntp .gt. 0 .and. imin .eq. 0) then
#endif /* AMOEBA */

    max_unfrag_mol_len = natom / numtasks + 1
!   max_unfrag_mol_len = natom / 64 + 1        ! for debugging on small setup
    if (max_unfrag_mol_len .lt. 3) max_unfrag_mol_len = 3
    frag_mol_cnt = 0
    max_mol_frag_cnt = 0

    do i = 1, nspm
      if (atm_nsp(i) .gt. max_unfrag_mol_len) then
        frag_mol_cnt = frag_mol_cnt + 1
        max_mol_frag_cnt = max_mol_frag_cnt + atm_nsp(i) / max_unfrag_mol_len
        if (mod(atm_nsp(i), max_unfrag_mol_len) .ne. 0) then
          max_mol_frag_cnt = max_mol_frag_cnt + 1
        end if
      end if
    end do

  else

    max_unfrag_mol_len = natom  ! for constant volume molecules don't matter.
    frag_mol_cnt = 0
    max_mol_frag_cnt = 0

  end if

#endif

  ! Allocate all memory needed for dynamics.

#ifdef AMOEBA
  call alloc_dynamics_mem(natom, nspm, ntp, imin, iamoeba,  num_ints, num_reals)
#else
  call alloc_dynamics_mem(natom, nspm, imin, ntp, num_ints, num_reals)
#endif /* AMOEBA */

  ! Create molecule atoms list if needed:

#ifdef AMOEBA
  if (ntp .gt. 0 .and. imin .eq. 0 .and. iamoeba .eq. 0) then
#else
  if (ntp .gt. 0 .and. imin .eq. 0) then
#endif /* AMOEBA */
    gbl_mol_atms(1) = 1
    do i = 2, nspm + 1
      gbl_mol_atms(i) = gbl_mol_atms(i - 1) + atm_nsp(i - 1)
    end do
  end if

  ! Calculate inverse, total masses:

  tmass = 0.d0

  j = 0

  do i = 1, nspm
    mol_mass_inv(i) = 0.d0
    mol_atm_cnt = atm_nsp(i)
    do n = 1, mol_atm_cnt
      j = j + 1
      mol_mass_inv(i) = mol_mass_inv(i) + mass(j)   ! Sum molecule.
      atm_mass_inv(j) = 1.d0 / mass(j)              ! Save inverted mass.
    end do
    tmass = tmass + mol_mass_inv(i)
    mol_mass_inv(i) = 1.d0 / mol_mass_inv(i)            ! Actually invert here.
  end do

  if (ntp .gt. 0 .and. imin .eq. 0) then
    gbl_mol_mass_inv(:) = mol_mass_inv(:)
  end if

  ! Initialize molecule COM for all atoms if needed:

#ifdef AMOEBA
  if (ntp .gt. 0 .and. imin .eq. 0 .and. iamoeba .eq. 0) then
#else
  if (ntp .gt. 0 .and. imin .eq. 0) then
#endif /* AMOEBA */
    call get_all_mol_com(nspm, crd, mass, gbl_mol_atms, gbl_mol_mass_inv, &
                         gbl_mol_com)
  end if

#ifdef MPI

  ! Okay, now actually find the molecules you have to fragment and fragment
  ! them, storing information in the gbl_frag_mol and gbl_mol_frags arrays.

  if (frag_mol_cnt .gt. 0) then

    ! First make the atom to residue map; this is a bit inefficient, but
    ! simplifies the heck out of things below, and only causes a small init
    ! cost in the master...

    do res_idx = 1, nres
      atm_idx = res_atms(res_idx)
      last_atm_idx = res_atms(res_idx + 1) - 1
      do while (atm_idx .le. last_atm_idx)
        atm_res_map(atm_idx) = res_idx
        atm_idx = atm_idx + 1
      end do
    end do

    gbl_mol_frags(1:max_mol_frag_cnt)%first_atm_id = 0
    gbl_mol_frags(1:max_mol_frag_cnt)%atm_cnt = 0
    gbl_mol_frags(1:max_mol_frag_cnt)%owner = -1

    frag_mol_cnt = 0    ! We redetermine the count of fragmented molecules.
    mfrag_idx = 1       ! This is the index into the whole fragments table.

    do i = 1, nspm

      ! Find molecules that need to be fragmented:

      if (atm_nsp(i) .gt. max_unfrag_mol_len) then

        frag_mol_cnt = frag_mol_cnt + 1
        mfrag_cnt = atm_nsp(i) / max_unfrag_mol_len
        if (mod(atm_nsp(i), max_unfrag_mol_len) .ne. 0) then
          mfrag_cnt = mfrag_cnt + 1
        end if
        gbl_frag_mols(frag_mol_cnt)%mol_idx = i
        gbl_frag_mols(frag_mol_cnt)%first_frag_idx = mfrag_idx
        gbl_frag_mols(frag_mol_cnt)%group = MPI_GROUP_NULL
        gbl_frag_mols(frag_mol_cnt)%communicator = MPI_COMM_NULL
        gbl_frag_mols(frag_mol_cnt)%task_cnt = 0
        
        first_atm = gbl_mol_atms(i)

        ! last_mfrag_idx is the last fragment in this fragmented molecule
        ! anticipated; it is possible there could be fewer fragments though,
        ! due to integral residue fragment requirements (worse case - imagine
        ! a huge molecule with 1 residue). The mfrag_cnt value here is the
        ! max value anticipated; it is reevaluated below.
        
        cur_mfrag_idx = mfrag_idx
        last_mfrag_idx = mfrag_idx + mfrag_cnt - 1
        mfrag_cnt = 0
        mol_atm_cnt = 0

        do while (cur_mfrag_idx .le. last_mfrag_idx)

          ! Make a guess about last_atm.  We know this is most likely
          ! incorrect, as we will round up or down to the nearest whole
          ! residue, but it helps us decide what is needed.
          
          last_atm = min(first_atm + max_unfrag_mol_len, gbl_mol_atms(i+1) - 1)

          first_res = atm_res_map(first_atm)
          last_res = atm_res_map(last_atm)

          ! If there is only one residue, or if we are dealing with the last
          ! fragment in this fragmented molecule, we just include all of the
          ! last residue; otherwise we check to see whether it is better to
          ! include it with this fragment or the next.

          if (first_res .ne. last_res .and. &
              cur_mfrag_idx .ne. last_mfrag_idx) then

            last_res_first_atm = res_atms(last_res)
            last_res_atm_cnt = res_atms(last_res + 1) - last_res_first_atm
            last_res_atm_cnt_incl = last_atm - last_res_first_atm + 1

            if (last_res_atm_cnt_incl * 2 .lt. last_res_atm_cnt) &
              last_res = last_res - 1

          end if

          atm_cnt = res_atms(last_res + 1) - first_atm
          gbl_mol_frags(cur_mfrag_idx)%first_atm_id = first_atm
          gbl_mol_frags(cur_mfrag_idx)%atm_cnt = atm_cnt
          first_atm = first_atm + atm_cnt
          mol_atm_cnt = mol_atm_cnt + atm_cnt

          mfrag_cnt = mfrag_cnt + 1
          cur_mfrag_idx = cur_mfrag_idx + 1

          ! BUGBUG - We need to write verification code that insures that
          ! atm_nsp() and gbl_res_atms() are consistent.  Then we should
          ! actually be able to just use mol_atm_cnt for loop control...

          if (mol_atm_cnt .ge. atm_nsp(i)) exit

        end do

        gbl_frag_mols(frag_mol_cnt)%frag_cnt = mfrag_cnt
        
        mfrag_idx = mfrag_idx + mfrag_cnt

      end if
    end do

  end if

#endif /* MPI */

! BEGIN DBG
#ifdef MPI
! write(0,*)'DBG: max_unfrag_mol_len =', max_unfrag_mol_len
! write(0,*)'DBG: frag_mol_cnt =', frag_mol_cnt
! write(0,*)'DBG: max_mol_frag_cnt =', max_mol_frag_cnt
!
! do i = 1, frag_mol_cnt
!
!   write(0,*)'DBG: frag mol', i, ' describes mol', gbl_frag_mols(i)%mol_idx
!   write(0,*)'DBG: frag cnt =', gbl_frag_mols(i)%frag_cnt
!   write(0,*)'DBG: first frag idx =', gbl_frag_mols(i)%first_frag_idx
!
!   do j = gbl_frag_mols(i)%first_frag_idx, &
!          gbl_frag_mols(i)%first_frag_idx + gbl_frag_mols(i)%frag_cnt - 1
!
!     write(0,2000)' DBG: frag, 1st atm, atm cnt =', j, &
!                   gbl_mol_frags(j)%first_atm_id, &
!                   gbl_mol_frags(j)%atm_cnt
!   end do
!
!2000 format(a, 3i7)
!
! end do
#endif /* MPI */
! END DBG

  return

end subroutine init_dynamics_dat

!*******************************************************************************
!
! Subroutine:  alloc_dynamics_mem
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef AMOEBA
subroutine alloc_dynamics_mem(natom, nspm, ntp, imin, iamoeba, &
                              num_ints, num_reals)
#else
subroutine alloc_dynamics_mem(natom, nspm, imin, ntp, &
                              num_ints, num_reals)
#endif /* AMOEBA */

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: natom
  integer, intent(in)           :: nspm
  integer, intent(in)           :: imin
  integer, intent(in)           :: ntp
#ifdef AMOEBA
  integer, intent(in)           :: iamoeba
#endif /* AMOEBA */

  ! num_ints and num_reals are used to return allocation counts. Don't zero.

  integer, intent(in out)       :: num_ints, num_reals     

! Local variables:

  integer                       :: alloc_failed

  ! Generic allocations required for just about everything:

  allocate(atm_mass_inv(natom), &
           stat = alloc_failed)

  if (alloc_failed .ne. 0) call setup_alloc_error

  num_reals = num_reals + size(atm_mass_inv)

  ! Relative coordinates should only be needed for constant pressure MD:

  if (ntp .gt. 0 .and. imin .eq. 0) then

    allocate(atm_rel_crd(3, natom), &
             gbl_mol_atms(nspm + 1), &
#ifdef MPI
             gbl_my_mol_lst(nspm), &
#endif
             gbl_mol_mass_inv(nspm), &
             gbl_mol_com(3, nspm), &
             stat = alloc_failed)

    if (alloc_failed .ne. 0) call setup_alloc_error

    num_ints = num_ints + size(gbl_mol_atms)
#ifdef MPI
    num_ints = num_ints + size(gbl_my_mol_lst)
#endif
    num_reals = num_reals + size(atm_rel_crd) + &
                            size(gbl_mol_mass_inv) + &
                            size(gbl_mol_com)

#ifdef MPI
    if (frag_mol_cnt .gt. 0) then

      allocate(gbl_frag_mols(frag_mol_cnt), &
               gbl_mol_frags(max_mol_frag_cnt), &
               gbl_my_frag_mol_lst(frag_mol_cnt), &
               stat = alloc_failed)
               
      num_ints = num_ints + size(gbl_frag_mols) * frag_mol_rec_ints + &
                            size(gbl_mol_frags) * mol_frag_rec_ints + &
                            size(gbl_my_frag_mol_lst)

    end if
#endif /* MPI */

  end if

  ! No need to initialize atm_rel_crd.
  ! Other stuff initialized in other code.

  return

end subroutine alloc_dynamics_mem

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_dynamics_dat
!
! Description: <TBS>
!              
!*******************************************************************************

#ifdef AMOEBA
subroutine bcast_dynamics_dat(natom, nspm, ntp, imin, iamoeba)
#else
subroutine bcast_dynamics_dat(natom, nspm, imin, ntp)
#endif /* AMOEBA */

  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer, intent(in)           :: natom
  integer, intent(in)           :: nspm
  integer, intent(in)           :: imin
  integer, intent(in)           :: ntp
#ifdef AMOEBA
  integer, intent(in)           :: iamoeba
#endif /* AMOEBA */

! Local variables:

  integer               :: num_ints, num_reals  ! returned values discarded
  integer               :: bytes_per_unit

  call mpi_bcast(tmass, dynamics_dat_dbl_cnt, mpi_double_precision, 0, &
                 mpi_comm_world, err_code_mpi)

  call mpi_bcast(max_unfrag_mol_len, dynamics_dat_int_cnt, mpi_integer, 0, &
                 mpi_comm_world, err_code_mpi)

  if (.not. master) then
    num_ints = 0
    num_reals = 0
#ifdef AMOEBA
    call alloc_dynamics_mem(natom, nspm, ntp, imin, iamoeba, &
                            num_ints, num_reals)
#else
    call alloc_dynamics_mem(natom, nspm, imin, ntp, num_ints, num_reals)
#endif /* AMOEBA */
  end if

  ! No need to broadcast atm_rel_crd.

  call mpi_bcast(atm_mass_inv, natom, mpi_double_precision, 0, &
                 mpi_comm_world, err_code_mpi)

  ! gbl_mol_atms and gbl_mol_mass_inv are only needed under CP MD:

#ifdef AMOEBA
  if (ntp .gt. 0 .and. imin .eq. 0 .and. iamoeba .eq. 0) then
#else
  if (ntp .gt. 0 .and. imin .eq. 0) then
#endif /* AMOEBA */

    call mpi_bcast(gbl_mol_atms, nspm + 1, mpi_integer, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_mol_mass_inv, nspm, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)
    call mpi_bcast(gbl_mol_com, nspm * 3, mpi_double_precision, 0, &
                   mpi_comm_world, err_code_mpi)

    if (frag_mol_cnt .gt. 0) then

      call get_bytesize(gbl_frag_mols(1), gbl_frag_mols(2), bytes_per_unit)
      call mpi_bcast(gbl_frag_mols, size(gbl_frag_mols) * bytes_per_unit, &
                     mpi_byte, 0, mpi_comm_world, err_code_mpi)

      call get_bytesize(gbl_mol_frags(1), gbl_mol_frags(2), bytes_per_unit)
      call mpi_bcast(gbl_mol_frags, size(gbl_mol_frags) * bytes_per_unit, &
                     mpi_byte, 0, mpi_comm_world, err_code_mpi)

    end if

  end if

  return

end subroutine bcast_dynamics_dat
#endif /* MPI */

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

    frag_cnt = gbl_frag_mols(mol_idx)%frag_cnt
    first_frag_idx = gbl_frag_mols(mol_idx)%first_frag_idx

    tasklist(1) = gbl_mol_frags(first_frag_idx)%owner
    task_cnt = 1

    ! Make a task list without any duplicates.  The algorithm is O(n**2), but
    ! n is typically pretty small, say 2-30.  All this hooplah is necessary
    ! because the target atom assignment may be larger or smaller than the
    ! default fragment size, which means that multiple atoms could get assigned
    ! to the same task.  Aside from having to do this duplicate search, that
    ! is not a bad thing.  The frequency of execution of this code should taper
    ! off quickly after initial loadbalancing.

    do frag_idx = first_frag_idx + 1, first_frag_idx + frag_cnt - 1
      taskid = gbl_mol_frags(frag_idx)%owner
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

      call mpi_group_incl(world_group, task_cnt, tasklist, &
                          gbl_frag_mols(mol_idx)%group, err_code_mpi)

      call mpi_comm_create(mpi_comm_world, gbl_frag_mols(mol_idx)%group, &
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

subroutine get_all_mol_com(mol_cnt, crd, mass, mol_atms, mol_mass_inv, mol_com)

  implicit none

! Formal arguments:

  integer                       :: mol_cnt
  double precision              :: crd(3, *)            ! atom crd array
  double precision              :: mass(*)              ! atom mass array
  integer                       :: mol_atms(*)
  double precision              :: mol_mass_inv(*)
  double precision              :: mol_com(3, *)

! Local variables:

  double precision              :: com(3)
  integer                       :: atm_idx, mol_idx

! Get COM for all molecules:

  do mol_idx = 1, mol_cnt
    com(:) = 0.d0
    do atm_idx = mol_atms(mol_idx), mol_atms(mol_idx + 1) - 1
      com(:) = com(:) + mass(atm_idx) * crd(:, atm_idx)
    end do
    mol_com(:, mol_idx) = com(:) * mol_mass_inv(mol_idx)
  end do

  return

end subroutine get_all_mol_com

end module dynamics_dat_mod
