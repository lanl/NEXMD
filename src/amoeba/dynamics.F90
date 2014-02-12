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
subroutine get_mol_com(mol_cnt, crd, mass, mol_atms, mol_mass_inv, my_mol_lst, &
                       mol_com)
#else
subroutine get_mol_com(mol_cnt, crd, mass, mol_atms, mol_mass_inv, mol_com)
#endif

  use dynamics_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer                       :: mol_cnt
  double precision              :: crd(3, *)
  double precision              :: mass(*)              ! atom mass array
  integer                       :: mol_atms(*)
  double precision              :: mol_mass_inv(*)
#ifdef MPI
  integer                       :: my_mol_lst(*)
#endif
  double precision              :: mol_com(3, *)

! Local variables:

  double precision              :: com(3)
  integer                       :: atm_idx_lo, atm_idx_hi, atm_idx
  integer                       :: mol_idx
#ifdef MPI
  integer                       :: i, j, k
  integer                       :: mylist_idx
  integer                       :: frag_cnt
  integer                       :: task_cnt
  integer                       :: first_frag_idx
  integer                       :: frag_idx
  integer                       :: taskid
  double precision, save        :: reduce_buf_in(3)
  double precision, save        :: reduce_buf_out(3)
#endif

! Get COM for molecules you own.

#ifdef MPI
  do mylist_idx = 1, mol_cnt
    mol_idx = my_mol_lst(mylist_idx)
#else
  do mol_idx = 1, mol_cnt
#endif
    atm_idx_lo = mol_atms(mol_idx)
    atm_idx_hi = mol_atms(mol_idx + 1) - 1

    com(:) = 0.d0

    do atm_idx = atm_idx_lo, atm_idx_hi
      com(:) = com(:) + mass(atm_idx) * crd(:, atm_idx)
    end do

    mol_com(:, mol_idx) = com(:) * mol_mass_inv(mol_idx)

  end do

#ifdef MPI

  do mylist_idx = 1, my_frag_mol_cnt

    mol_idx = gbl_my_frag_mol_lst(mylist_idx)
    frag_cnt = gbl_frag_mols(mol_idx)%frag_cnt
    first_frag_idx = gbl_frag_mols(mol_idx)%first_frag_idx
    task_cnt = gbl_frag_mols(mol_idx)%task_cnt

    com(:) = 0.d0

    if (task_cnt .gt. 1) then

      do frag_idx = first_frag_idx, first_frag_idx + frag_cnt - 1
        ! This task id is in the "world" context.
        taskid = gbl_mol_frags(frag_idx)%owner

        if (taskid .eq. mytaskid) then

          atm_idx_lo = gbl_mol_frags(frag_idx)%first_atm_id
          atm_idx_hi = atm_idx_lo + gbl_mol_frags(frag_idx)%atm_cnt - 1

          do atm_idx = atm_idx_lo, atm_idx_hi
            com(:) = com(:) + mass(atm_idx) * crd(:, atm_idx)
          end do

        end if

      end do

      reduce_buf_in(:) = com(:)
      call mpi_allreduce(reduce_buf_in, reduce_buf_out, 3, &
                         mpi_double_precision, mpi_sum, &
                         gbl_frag_mols(mol_idx)%communicator, err_code_mpi)
      com(:) = reduce_buf_out(:)

    else

      ! All the fragments are owned by this task...

      do frag_idx = first_frag_idx, first_frag_idx + frag_cnt - 1

        atm_idx_lo = gbl_mol_frags(frag_idx)%first_atm_id
        atm_idx_hi = atm_idx_lo + gbl_mol_frags(frag_idx)%atm_cnt - 1

        do atm_idx = atm_idx_lo, atm_idx_hi
          com(:) = com(:) + mass(atm_idx) * crd(:, atm_idx)
        end do

      end do

    end if

    mol_com(:, gbl_frag_mols(mol_idx)%mol_idx) = &
      com(:) * mol_mass_inv(gbl_frag_mols(mol_idx)%mol_idx)

  end do

#endif /* MPI */

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
subroutine get_ekcom(mol_cnt, mol_atms, tma_inv, ekcmt, vel, mass, my_mol_lst)
#else
subroutine get_ekcom(mol_cnt, mol_atms, tma_inv, ekcmt, vel, mass)
#endif

  use dynamics_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: mol_cnt
  integer               :: mol_atms(*)
  double precision      :: tma_inv(*)
  double precision      :: ekcmt(3)
  double precision      :: vel(3, *)
  double precision      :: mass(*)              ! atom mass array
#ifdef MPI
  integer               :: my_mol_lst(*)
#endif

! Local variables:

  double precision      :: ekcml(3)
  double precision      :: vcm(3)

  integer               :: atm_idx, atm_idx_lo, atm_idx_hi
  integer               :: mol_idx
#ifdef MPI
  integer                       :: i, j, k      ! DBG
  integer                       :: mylist_idx
  integer                       :: frag_cnt
  integer                       :: first_frag_idx
  integer                       :: frag_idx
  integer                       :: task_cnt
  integer                       :: taskid
  double precision, save        :: reduce_buf_in(3)
  double precision, save        :: reduce_buf_out(3)
#endif

  ekcml(:) = 0.d0

#ifdef MPI
  do mylist_idx = 1, mol_cnt
    mol_idx = my_mol_lst(mylist_idx)
#else
  do mol_idx = 1, mol_cnt
#endif

    atm_idx_lo = mol_atms(mol_idx)
    atm_idx_hi = mol_atms(mol_idx + 1) - 1

    vcm(:) = 0.d0

    do atm_idx = atm_idx_lo, atm_idx_hi
      vcm(:) = vcm(:) + vel(:, atm_idx) * mass(atm_idx)
    end do

    ekcml(:) = ekcml(:) + tma_inv(mol_idx) * vcm(:) * vcm(:)

  end do

#ifdef MPI

  do mylist_idx = 1, my_frag_mol_cnt

    mol_idx = gbl_my_frag_mol_lst(mylist_idx)
    frag_cnt = gbl_frag_mols(mol_idx)%frag_cnt
    first_frag_idx = gbl_frag_mols(mol_idx)%first_frag_idx
    task_cnt = gbl_frag_mols(mol_idx)%task_cnt

    vcm(:) = 0.d0

    if (task_cnt .gt. 1) then

      do frag_idx = first_frag_idx, first_frag_idx + frag_cnt - 1
        ! This task id is in the "world" context.
        taskid = gbl_mol_frags(frag_idx)%owner

        if (taskid .eq. mytaskid) then

          atm_idx_lo = gbl_mol_frags(frag_idx)%first_atm_id
          atm_idx_hi = atm_idx_lo + gbl_mol_frags(frag_idx)%atm_cnt - 1

          do atm_idx = atm_idx_lo, atm_idx_hi
            vcm(:) = vcm(:) + vel(:, atm_idx) * mass(atm_idx)
          end do

        end if

      end do

      reduce_buf_in(:) = vcm(:)
      call mpi_reduce(reduce_buf_in, reduce_buf_out, 3, &
                         mpi_double_precision, mpi_sum, 0, &
                         gbl_frag_mols(mol_idx)%communicator, err_code_mpi)

      if (mytaskid .eq. gbl_mol_frags(first_frag_idx)%owner) then
        vcm(:) = reduce_buf_out(:)
        ekcml(:) = ekcml(:) + &
                   tma_inv(gbl_frag_mols(mol_idx)%mol_idx) * vcm(:) * vcm(:)
      end if

    else

      ! All the fragments are owned by this task...

      do frag_idx = first_frag_idx, first_frag_idx + frag_cnt - 1

        atm_idx_lo = gbl_mol_frags(frag_idx)%first_atm_id
        atm_idx_hi = atm_idx_lo + gbl_mol_frags(frag_idx)%atm_cnt - 1

        do atm_idx = atm_idx_lo, atm_idx_hi
          vcm(:) = vcm(:) + vel(:, atm_idx) * mass(atm_idx)
        end do

      end do

      ekcml(:) = ekcml(:) + &
                 tma_inv(gbl_frag_mols(mol_idx)%mol_idx) * vcm(:) * vcm(:)

    end if

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
subroutine get_atm_rel_crd(mol_cnt, mol_atms, mol_com, crd, rel_crd, my_mol_lst)
#else
subroutine get_atm_rel_crd(mol_cnt, mol_atms, mol_com, crd, rel_crd)
#endif

  use dynamics_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  integer               :: mol_cnt
  integer               :: mol_atms(*)
  double precision      :: mol_com(3, *)
  double precision      :: crd(3, *)
  double precision      :: rel_crd(3, *)
#ifdef MPI
  integer               :: my_mol_lst(*)
#endif

! Local variables:

  integer               :: atm_idx, mol_idx
#ifdef MPI
  integer               :: mylist_idx
  integer               :: frag_cnt
  integer               :: first_frag_idx
  integer               :: frag_idx
  integer               :: atm_idx_lo
  integer               :: atm_idx_hi
#endif

#ifdef MPI
  do mylist_idx = 1, mol_cnt
    mol_idx = my_mol_lst(mylist_idx)
#else
  do mol_idx = 1, mol_cnt
#endif

    do atm_idx = mol_atms(mol_idx), mol_atms(mol_idx + 1) - 1
      rel_crd(:, atm_idx) = crd(:, atm_idx) - mol_com(:, mol_idx)
    end do

  end do

#ifdef MPI
  do mylist_idx = 1, my_frag_mol_cnt

    mol_idx = gbl_my_frag_mol_lst(mylist_idx)
    frag_cnt = gbl_frag_mols(mol_idx)%frag_cnt
    first_frag_idx = gbl_frag_mols(mol_idx)%first_frag_idx

    do frag_idx = first_frag_idx, first_frag_idx + frag_cnt - 1
      if (gbl_mol_frags(frag_idx)%owner .eq. mytaskid) then
        atm_idx_lo = gbl_mol_frags(frag_idx)%first_atm_id
        atm_idx_hi = atm_idx_lo + gbl_mol_frags(frag_idx)%atm_cnt - 1
        do atm_idx = atm_idx_lo, atm_idx_hi
          rel_crd(:, atm_idx) = crd(:, atm_idx) - &
                                mol_com(:, gbl_frag_mols(mol_idx)%mol_idx)
        end do
      end if
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

  do j = 1, atm_cnt

#ifdef MPI
  ! In order to generate the same sequence of pseudorandom numbers that you
  ! would using a single processor or any other combo of multiple processors,
  ! you have to go through the atoms in order.  The unused results are not
  ! returned.

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
#else
    wfac = mass_inv(j) * dtx
    aamass = mass(j)
    rsd = sdfac * sqrt(aamass)
    call gauss(0.d0, rsd, fln1)
    call gauss(0.d0, rsd, fln2)
    call gauss(0.d0, rsd, fln3)
    vel(1,j) = (vel(1,j) * c_explic + (frc(1,j) + fln1) * wfac) * c_implic
    vel(2,j) = (vel(2,j) * c_explic + (frc(2,j) + fln2) * wfac) * c_implic
    vel(3,j) = (vel(3,j) * c_explic + (frc(3,j) + fln3) * wfac) * c_implic
#endif

  end do

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

#ifdef MPI
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
#ifdef AMOEBA
  use mdin_ctrl_dat_mod
  use mdin_amoeba_dat_mod
#endif /* AMOEBA */

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

#ifdef AMOEBA
    if (iamoeba .eq. 1) then
      if (beeman_integrator .eq. 1) then
        do i = 1, atm_cnt
          vel(:, i) = 20.455d0 * vel(:, i)
        end do
      end if
    end if
#endif /* AMOEBA */

  end if
   
  return

end subroutine all_atom_setvel

end module dynamics_mod
