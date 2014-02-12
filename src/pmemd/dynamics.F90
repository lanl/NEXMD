#include "copyright.i"

!*******************************************************************************
!
! Module:  dynamics_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module dynamics_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:   get_mol_com
!
! Description:  Pressure scaling routine for crds. ONLY used for constant
!               pressure scaling (ntp .gt. 0)!
!
!*******************************************************************************

#ifdef MPI
subroutine get_mol_com(mol_cnt, crd, mass, mol_mass_inv, my_mol_lst, &
                       mol_com)
#else
subroutine get_mol_com(mol_cnt, crd, mass, mol_mass_inv, mol_com)
#endif

  use dynamics_dat_mod
  use parallel_dat_mod
  use prfs_mod
  use mol_list_mod

  implicit none

! Formal arguments:

  integer, intent(in)                   :: mol_cnt
  double precision, intent(in)          :: crd(3, *)
  double precision, intent(in)          :: mass(*)      ! atom mass array
  double precision, intent(in)          :: mol_mass_inv(*)
#ifdef MPI
  integer, intent(in)                   :: my_mol_lst(*)
#endif
  double precision, intent(out)         :: mol_com(3, *)

! Local variables:

  double precision              :: com(3)
  integer                       :: i                    ! mol atm list idx
  integer                       :: atm_id, mol_id
  integer                       :: mol_atm_cnt, offset
#ifdef MPI
  integer                       :: j
  integer                       :: prf_id
  integer                       :: mol_offset, mol_listcnt
  integer                       :: prf_offset, prf_listcnt
  integer                       :: mylist_idx
  integer                       :: molfrag_idx
  integer                       :: frag_mol_idx
  integer                       :: task_cnt
  double precision, save        :: reduce_buf_in(3)
  double precision, save        :: reduce_buf_out(3)
#endif

! Get COM for molecules you own.

#ifdef MPI
  do mylist_idx = 1, mol_cnt
    mol_id = my_mol_lst(mylist_idx)
    mol_offset = gbl_mol_prfs_listdata(mol_id)%offset
    mol_listcnt = gbl_mol_prfs_listdata(mol_id)%cnt
    com(:) = 0.d0
    do i = mol_offset + 1, mol_offset + mol_listcnt
      prf_id = gbl_mol_prfs_lists(i)
      prf_offset = gbl_prf_listdata(prf_id)%offset
      prf_listcnt = gbl_prf_listdata(prf_id)%cnt
      do j = prf_offset + 1, prf_offset + prf_listcnt
        atm_id = gbl_prf_lists(j)
        com(:) = com(:) + mass(atm_id) * crd(:, atm_id)
      end do
    end do
    mol_com(:, mol_id) = com(:) * mol_mass_inv(mol_id)
  end do

  do mylist_idx = 1, my_frag_mol_cnt

    frag_mol_idx = gbl_my_frag_mol_lst(mylist_idx)
    task_cnt = gbl_frag_mols(frag_mol_idx)%task_cnt

    com(:) = 0.d0

    if (task_cnt .gt. 1) then

      do molfrag_idx = gbl_frag_mols(frag_mol_idx)%first_molfrag_idx, &
                       gbl_frag_mols(frag_mol_idx)%first_molfrag_idx + &
                       gbl_frag_mols(frag_mol_idx)%molfrag_cnt - 1

        ! This task id is in the "world" context.
        if (gbl_molfrags(molfrag_idx)%owner .eq. mytaskid) then
          mol_offset = gbl_molfrags(molfrag_idx)%offset
          mol_listcnt = gbl_molfrags(molfrag_idx)%cnt
          do i = mol_offset + 1, mol_offset + mol_listcnt
            prf_id = gbl_mol_prfs_lists(i)
            prf_offset = gbl_prf_listdata(prf_id)%offset
            prf_listcnt = gbl_prf_listdata(prf_id)%cnt
            do j = prf_offset + 1, prf_offset + prf_listcnt
              atm_id = gbl_prf_lists(j)
              com(:) = com(:) + mass(atm_id) * crd(:, atm_id)
            end do
          end do
        end if

      end do

      reduce_buf_in(:) = com(:)
      call mpi_allreduce(reduce_buf_in, reduce_buf_out, 3, &
                         mpi_double_precision, mpi_sum, &
                         gbl_frag_mols(frag_mol_idx)%communicator, err_code_mpi)
      com(:) = reduce_buf_out(:)

    else

      ! All the fragments are owned by this task...

      do molfrag_idx = gbl_frag_mols(frag_mol_idx)%first_molfrag_idx, &
                       gbl_frag_mols(frag_mol_idx)%first_molfrag_idx + &
                       gbl_frag_mols(frag_mol_idx)%molfrag_cnt - 1

        mol_offset = gbl_molfrags(molfrag_idx)%offset
        mol_listcnt = gbl_molfrags(molfrag_idx)%cnt
        do i = mol_offset + 1, mol_offset + mol_listcnt
          prf_id = gbl_mol_prfs_lists(i)
          prf_offset = gbl_prf_listdata(prf_id)%offset
          prf_listcnt = gbl_prf_listdata(prf_id)%cnt
          do j = prf_offset + 1, prf_offset + prf_listcnt
            atm_id = gbl_prf_lists(j)
            com(:) = com(:) + mass(atm_id) * crd(:, atm_id)
          end do
        end do

      end do

    end if

    mol_com(:, gbl_frag_mols(frag_mol_idx)%mol_idx) = &
      com(:) * mol_mass_inv(gbl_frag_mols(frag_mol_idx)%mol_idx)

  end do

#else
  do mol_id = 1, mol_cnt

    com(:) = 0.d0
    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      com(:) = com(:) + mass(atm_id) * crd(:, atm_id)
    end do

    mol_com(:, mol_id) = com(:) * mol_mass_inv(mol_id)

  end do
#endif

  return

end subroutine get_mol_com

!*******************************************************************************
!
! Subroutine:  get_ekcom
!
! Description:  Routine to calculate the total kinetic energy of the center of
!               mass of the sub-molecules and also the coordinates of the
!               molecules relative to the center of mass.
!*******************************************************************************

#ifdef MPI
subroutine get_ekcom(mol_cnt, tma_inv, ekcmt, vel, mass, my_mol_lst)
#else
subroutine get_ekcom(mol_cnt, tma_inv, ekcmt, vel, mass)
#endif

  use dynamics_dat_mod
  use parallel_dat_mod
  use prfs_mod
  use mol_list_mod

  implicit none

! Formal arguments:

  integer               :: mol_cnt
  double precision      :: tma_inv(*)
  double precision      :: ekcmt(3)
  double precision      :: vel(3, *)
  double precision      :: mass(*)              ! atom mass array
#ifdef MPI
  integer               :: my_mol_lst(*)
#endif

! Local variables:

  double precision              :: ekcml(3)
  double precision              :: vcm(3)

  integer                       :: i                    ! mol atm list idx
  integer                       :: atm_id, mol_id
  integer                       :: mol_atm_cnt, offset
#ifdef MPI
  integer                       :: j
  integer                       :: prf_id
  integer                       :: mol_offset, mol_listcnt
  integer                       :: prf_offset, prf_listcnt
  integer                       :: mylist_idx
  integer                       :: molfrag_idx
  integer                       :: frag_mol_idx
  integer                       :: task_cnt
  integer                       :: first_frag_idx
  double precision, save        :: reduce_buf_in(3)
  double precision, save        :: reduce_buf_out(3)
#endif

  ekcml(:) = 0.d0

#ifdef MPI
  do mylist_idx = 1, mol_cnt
    mol_id = my_mol_lst(mylist_idx)
    vcm(:) = 0.d0
    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      vcm(:) = vcm(:) + vel(:, atm_id) * mass(atm_id)
    end do
    ekcml(:) = ekcml(:) + tma_inv(mol_id) * vcm(:) * vcm(:)
  end do

  do mylist_idx = 1, my_frag_mol_cnt

    frag_mol_idx = gbl_my_frag_mol_lst(mylist_idx)
    task_cnt = gbl_frag_mols(frag_mol_idx)%task_cnt

    vcm(:) = 0.d0

    if (task_cnt .gt. 1) then

      do molfrag_idx = gbl_frag_mols(frag_mol_idx)%first_molfrag_idx, &
                       gbl_frag_mols(frag_mol_idx)%first_molfrag_idx + &
                       gbl_frag_mols(frag_mol_idx)%molfrag_cnt - 1

        ! This task id is in the "world" context.
        if (gbl_molfrags(molfrag_idx)%owner .eq. mytaskid) then
          mol_offset = gbl_molfrags(molfrag_idx)%offset
          mol_listcnt = gbl_molfrags(molfrag_idx)%cnt
          do i = mol_offset + 1, mol_offset + mol_listcnt
            prf_id = gbl_mol_prfs_lists(i)
            prf_offset = gbl_prf_listdata(prf_id)%offset
            prf_listcnt = gbl_prf_listdata(prf_id)%cnt
            do j = prf_offset + 1, prf_offset + prf_listcnt
              atm_id = gbl_prf_lists(j)
              vcm(:) = vcm(:) + vel(:, atm_id) * mass(atm_id)
            end do
          end do
        end if

      end do

      reduce_buf_in(:) = vcm(:)
      call mpi_reduce(reduce_buf_in, reduce_buf_out, 3, &
                         mpi_double_precision, mpi_sum, 0, &
                         gbl_frag_mols(frag_mol_idx)%communicator, err_code_mpi)

      first_frag_idx = gbl_frag_mols(frag_mol_idx)%first_molfrag_idx 
      if (mytaskid .eq. gbl_molfrags(first_frag_idx)%owner) then
        vcm(:) = reduce_buf_out(:)
        ekcml(:) = ekcml(:) + tma_inv(gbl_frag_mols(frag_mol_idx)%mol_idx) * &
                   vcm(:) * vcm(:)
      end if

    else

      ! All the fragments are owned by this task...

      do molfrag_idx = gbl_frag_mols(frag_mol_idx)%first_molfrag_idx, &
                       gbl_frag_mols(frag_mol_idx)%first_molfrag_idx + &
                       gbl_frag_mols(frag_mol_idx)%molfrag_cnt - 1


        mol_offset = gbl_molfrags(molfrag_idx)%offset
        mol_listcnt = gbl_molfrags(molfrag_idx)%cnt
        do i = mol_offset + 1, mol_offset + mol_listcnt
          prf_id = gbl_mol_prfs_lists(i)
          prf_offset = gbl_prf_listdata(prf_id)%offset
          prf_listcnt = gbl_prf_listdata(prf_id)%cnt
          do j = prf_offset + 1, prf_offset + prf_listcnt
            atm_id = gbl_prf_lists(j)
            vcm(:) = vcm(:) + vel(:, atm_id) * mass(atm_id)
          end do
        end do

      end do
      ekcml(:) = ekcml(:) + tma_inv(gbl_frag_mols(frag_mol_idx)%mol_idx) * &
                 vcm(:) * vcm(:)
    end if

  end do

#else
  do mol_id = 1, mol_cnt
    vcm(:) = 0.d0
    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      vcm(:) = vcm(:) + vel(:, atm_id) * mass(atm_id)
    end do

    ekcml(:) = ekcml(:) + tma_inv(mol_id) * vcm(:) * vcm(:)

  end do
#endif /* MPI */


  ekcmt(:) = ekcml(:)

  return

end subroutine get_ekcom

!*******************************************************************************
!
! Subroutine:  get_atm_rel_crd
!
! Description:  Routine to calculate the coordinate relative to the COM for
!               molecules.  This gets used in molecular virial calcs.
!*******************************************************************************

#ifdef MPI
subroutine get_atm_rel_crd(mol_cnt, mol_com, crd, rel_crd, my_mol_lst)
#else
subroutine get_atm_rel_crd(mol_cnt, mol_com, crd, rel_crd)
#endif

  use dynamics_dat_mod
  use parallel_dat_mod
  use prfs_mod
  use mol_list_mod

  implicit none

! Formal arguments:

  integer               :: mol_cnt
  double precision      :: mol_com(3, *)
  double precision      :: crd(3, *)
  double precision      :: rel_crd(3, *)
#ifdef MPI
  integer               :: my_mol_lst(*)
#endif

! Local variables:

  integer               :: i                    ! mol atm list idx
  integer               :: atm_id, mol_id
  integer               :: mol_atm_cnt, offset
#ifdef MPI
  integer               :: j
  integer               :: prf_id
  integer               :: mol_offset, mol_listcnt
  integer               :: prf_offset, prf_listcnt
  integer               :: mylist_idx
  integer               :: molfrag_idx
  integer               :: frag_mol_idx
#endif

#ifdef MPI
  do mylist_idx = 1, mol_cnt
    mol_id = my_mol_lst(mylist_idx)
    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      rel_crd(:, atm_id) = crd(:, atm_id) - mol_com(:, mol_id)
    end do
  end do

  do mylist_idx = 1, my_frag_mol_cnt

    frag_mol_idx = gbl_my_frag_mol_lst(mylist_idx)

    do molfrag_idx = gbl_frag_mols(frag_mol_idx)%first_molfrag_idx, &
                     gbl_frag_mols(frag_mol_idx)%first_molfrag_idx + &
                     gbl_frag_mols(frag_mol_idx)%molfrag_cnt - 1

      if (gbl_molfrags(molfrag_idx)%owner .eq. mytaskid) then
        mol_offset = gbl_molfrags(molfrag_idx)%offset
        mol_listcnt = gbl_molfrags(molfrag_idx)%cnt
        do i = mol_offset + 1, mol_offset + mol_listcnt
          prf_id = gbl_mol_prfs_lists(i)
          prf_offset = gbl_prf_listdata(prf_id)%offset
          prf_listcnt = gbl_prf_listdata(prf_id)%cnt
          do j = prf_offset + 1, prf_offset + prf_listcnt
            atm_id = gbl_prf_lists(j)
            rel_crd(:, atm_id) = crd(:, atm_id) - &
                                 mol_com(:, gbl_frag_mols(frag_mol_idx)%mol_idx)
          end do
        end do
      end if

    end do

  end do

#else
  do mol_id = 1, mol_cnt
    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      rel_crd(:, atm_id) = crd(:, atm_id) - mol_com(:, mol_id)
    end do
  end do
#endif /* MPI */

  return

end subroutine get_atm_rel_crd

!*******************************************************************************
!
! Subroutine:  langevin_setvel
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine langevin_setvel(atm_cnt, vel, frc, mass, mass_inv, &
                           dt, temp0, gamma_ln)

  use parallel_dat_mod
  use random_mod
  use mdin_ctrl_dat_mod, only : no_ntt3_sync

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: vel(3, atm_cnt)
  double precision      :: frc(3, atm_cnt)
  double precision      :: mass(atm_cnt)
  double precision      :: mass_inv(atm_cnt)
  double precision      :: dt
  double precision      :: temp0
  double precision      :: gamma_ln

! Local variables:

  double precision      :: aamass
  double precision      :: boltz2
  double precision      :: c_explic
  double precision      :: c_implic
  double precision      :: dtx
  double precision      :: fln1, fln2, fln3
  double precision      :: gammai
  double precision      :: half_dtx
  double precision      :: rsd
  double precision      :: sdfac
  double precision      :: wfac
  integer               :: j


  boltz2 = 8.31441d-3 * 0.5d0 / 4.184d0
  gammai = gamma_ln / 20.455d0
  dtx = dt * 20.455d+00
  half_dtx = dtx * 0.5d0
  c_implic = 1.d0 / (1.d0 + gammai * half_dtx)
  c_explic = 1.d0 - gammai * half_dtx
  sdfac = sqrt(4.d0 * gammai * boltz2 * temp0 / dtx)

! Split here depending on whether we are synching the random
! number stream across MPI tasks for ntt=3. If ig=-1 then we
! do not sync. This gives better scaling. We duplicate code here
! to avoid an if statement in the inner loop.
#if defined(MPI) && !defined(CUDA)
  if (no_ntt3_sync) then
    do j = 1, atm_cnt
      if (gbl_atm_owner_map(j) .eq. mytaskid) then
        wfac = mass_inv(j) * dtx
        aamass = mass(j)
        rsd = sdfac * sqrt(aamass)
        call gauss(0.d0, rsd, fln1)
        call gauss(0.d0, rsd, fln2)
        call gauss(0.d0, rsd, fln3)
        vel(1,j) = (vel(1,j) * c_explic + (frc(1,j) + fln1) * wfac) * c_implic
        vel(2,j) = (vel(2,j) * c_explic + (frc(2,j) + fln2) * wfac) * c_implic
        vel(3,j) = (vel(3,j) * c_explic + (frc(3,j) + fln3) * wfac) * c_implic
      end if
    end do
  else
    do j = 1, atm_cnt
      if (gbl_atm_owner_map(j) .eq. mytaskid) then
        wfac = mass_inv(j) * dtx
        aamass = mass(j)
        rsd = sdfac * sqrt(aamass)
        call gauss(0.d0, rsd, fln1)
        call gauss(0.d0, rsd, fln2)
        call gauss(0.d0, rsd, fln3)
        vel(1,j) = (vel(1,j) * c_explic + (frc(1,j) + fln1) * wfac) * c_implic
        vel(2,j) = (vel(2,j) * c_explic + (frc(2,j) + fln2) * wfac) * c_implic
        vel(3,j) = (vel(3,j) * c_explic + (frc(3,j) + fln3) * wfac) * c_implic
      else
        call gauss(0.d0, 1.d0)
        call gauss(0.d0, 1.d0)
        call gauss(0.d0, 1.d0)
      end if
    end do
  end if
#else
  do j = 1, atm_cnt
    wfac = mass_inv(j) * dtx
    aamass = mass(j)
    rsd = sdfac * sqrt(aamass)
    call gauss(0.d0, rsd, fln1)
    call gauss(0.d0, rsd, fln2)
    call gauss(0.d0, rsd, fln3)
    vel(1,j) = (vel(1,j) * c_explic + (frc(1,j) + fln1) * wfac) * c_implic
    vel(2,j) = (vel(2,j) * c_explic + (frc(2,j) + fln2) * wfac) * c_implic
    vel(3,j) = (vel(3,j) * c_explic + (frc(3,j) + fln3) * wfac) * c_implic
  end do
#endif

  return

end subroutine langevin_setvel

!*******************************************************************************
!
! Subroutine:   vrand_set_velocities
!
! Description:  Assign velocities from a Maxwellian distribution.
!              
!*******************************************************************************

subroutine vrand_set_velocities(atm_cnt, vel, mass_inv, temp)
   
  use parallel_dat_mod
  use random_mod

  implicit none

  integer               :: atm_cnt
  double precision      :: vel(3, atm_cnt)
  double precision      :: mass_inv(atm_cnt)
  double precision      :: temp

  double precision      :: boltz
  double precision      :: sd
  integer               :: j
   
  if (temp .lt. 1.d-6) then

    vel(:,:) = 0.d0

  else

    boltz = 8.31441d-3 * temp / 4.184d0

    do j = 1, atm_cnt

#if defined(MPI) && !defined(CUDA)
  ! In order to generate the same sequence of pseudorandom numbers that you
  ! would using a single processor or any other combo of multiple processors,
  ! you have to go through the atoms in order.  The unused results are not
  ! returned.

      if (gbl_atm_owner_map(j) .eq. mytaskid) then
        sd =  sqrt(boltz * mass_inv(j))
        call gauss(0.d0, sd, vel(1, j))
        call gauss(0.d0, sd, vel(2, j))
        call gauss(0.d0, sd, vel(3, j))
      else
        call gauss(0.d0, 1.d0)
        call gauss(0.d0, 1.d0)
        call gauss(0.d0, 1.d0)
      end if
#else
      sd =  sqrt(boltz * mass_inv(j))
      call gauss(0.d0, sd, vel(1, j))
      call gauss(0.d0, sd, vel(2, j))
      call gauss(0.d0, sd, vel(3, j))
#endif

    end do

  end if

  return

end subroutine vrand_set_velocities

!*******************************************************************************
!
! Subroutine:   all_atom_setvel
!
! Description:  Assign velocities from a Maxwellian distribution.
!              
!*******************************************************************************

subroutine all_atom_setvel(atm_cnt, vel, mass_inv, temp)
   
  use random_mod

  implicit none

  integer               :: atm_cnt
  double precision      :: vel(3, atm_cnt)
  double precision      :: mass_inv(atm_cnt)
  double precision      :: temp

  double precision      :: boltz
  double precision      :: sd
  integer               :: i
   
  if (temp .lt. 1.d-6) then

    vel(:,:) = 0.d0

  else
   
    boltz = 8.31441d-3 * temp / 4.184d0

    do i = 1, atm_cnt
      sd =  sqrt(boltz * mass_inv(i))
      call gauss(0.d0, sd, vel(1, i))
      call gauss(0.d0, sd, vel(2, i))
      call gauss(0.d0, sd, vel(3, i))
    end do

  end if
   
  return

end subroutine all_atom_setvel

end module dynamics_mod
