!*******************************************************************************
!
! Module: pbc_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module pbc_mod

use file_io_dat_mod
use gbl_constants_mod

  implicit none

  ! The following storage is per-process common; ie., it SHOULD be
  ! broadcast from the master to the other processes!

  integer, parameter    :: pbc_int_cnt = 1

  integer            is_orthog

  common / pbc_int / is_orthog

  save  :: / pbc_int /

  integer, parameter    :: pbc_dbl_cnt = 32

  double precision   recip(3, 3), ucell(3, 3), pbc_box(3), cut_factor(3), &
                     reclng(3), pbc_alpha, pbc_beta, pbc_gamma, &
                     uc_volume, uc_sphere

  common / pbc_dbl / recip, ucell, pbc_box, cut_factor, &
                     reclng, pbc_alpha, pbc_beta, pbc_gamma, &
                     uc_volume, uc_sphere

  save  :: / pbc_dbl /

  ! Defining quantities for unit cell:
  !
  !   ucell is the 3x3 of direct lattice vectors.
  !   recip are the 3x3 of reciprocal lattice vectors.
  !   cut_factor handles "spherical cutoff protusion" in nonorthogonal unit cell

! Data that is not broadcast:

  double precision, save        :: last_recip(3, 3)

contains

!*******************************************************************************
!
! Subroutine:  init_pbc
!
! Description: This routine stores pbc unit cell box lengths, angles and
!              produces the direct and reciprocal lattice vectors from the unit
!              cell edge lengths and angles which are passed to it.  It is
!              assumed that the 1st vector (length a) lies along the cartesian
!              x-axis the 2nd vector (length b) is in the x-y plane with
!              positive y, and that the direct lattice vectors are a
!              non-degenerate right handed system.  Thus the 3rd vector has
!              positive z component.  Alpha is the angle (in degrees) between
!              2nd and 3rd vectors, beta is the angle (in degrees) between 1st
!              and 3rd vectors, and gamma is the angle (in degrees) between 1st
!              and 2nd vectors.  The lengths of the reciprocal vectors are
!              reclng(1),reclng(2) and reclng(3) (local to this routine).
!
! NOTE - If the unit cell is orthogonal, we do simpler calcs that ensure that
!        the off-diagonal elements are 0.d0.  This should improve accuracy.
!
!*******************************************************************************

subroutine init_pbc(a, b, c, alpha, beta, gamma, max_cutoff)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: a
  double precision, intent(in)  :: b
  double precision, intent(in)  :: c
  double precision, intent(in)  :: alpha
  double precision, intent(in)  :: beta
  double precision, intent(in)  :: gamma
  double precision, intent(in)  :: max_cutoff

! Local variables:

  double precision      :: distance
  double precision      :: factor
  double precision      :: u23(3), u31(3), u12(3)
  double precision      :: result
  integer               :: i

! Calculate/store box information:

  pbc_box(1) = a
  pbc_box(2) = b
  pbc_box(3) = c
  pbc_alpha = alpha
  pbc_beta = beta
  pbc_gamma = gamma

  if (alpha .eq. 90.d0 .and. beta  .eq. 90.d0 .and. gamma .eq. 90.d0) then
    is_orthog = 1
  else
    is_orthog = 0
  end if

  if (is_orthog .ne. 0) then
    ucell(:, :) = 0.d0
    ucell(1, 1) = a
    ucell(2, 2) = b
    ucell(3, 3) = c
    cut_factor(:) = 1.d0
  else
    factor = PI / 180.d0
    ucell(1, 1) = a
    ucell(2, 1) = 0.d0
    ucell(3, 1) = 0.d0
    ucell(1, 2) = b * cos(factor * gamma)
    ucell(2, 2) = b * sin(factor * gamma)
    ucell(3, 2) = 0.d0
    ucell(1, 3) = c * cos(factor * beta)
    ucell(2, 3) =  (b * c * cos(factor * alpha) - ucell(1, 3) * &
                   ucell(1, 2))/ucell(2, 2)
    ucell(3, 3) = sqrt(c * c - ucell(1, 3) * ucell(1, 3) - ucell(2, 3) * &
                  ucell(2, 3))

    ! Cut factors are used to correct for "spherical cutoff protrusion" into
    ! adjacent unit cells.  The problem is that the point at which a cutoff
    ! sphere is tangent to a unit cell side is not the contact point for 
    ! projection of an orthogonal x, y, or z vector in a nonorthogonal unit
    ! cell.  We thus have to increase the cutoff a bit to allow for the longer
    ! distance for the orthogonal projection.

    cut_factor(1) = 1.d0 / (sin(factor * beta) * sin(factor * gamma))
    cut_factor(2) = 1.d0 / (sin(factor * alpha) * sin(factor * gamma))
    cut_factor(3) = 1.d0 / (sin(factor * alpha) * sin(factor * beta))
  end if

! Now get reciprocal vectors:

  if (is_orthog .ne. 0) then
    recip(:, :) = 0.d0
    recip(1, 1) = 1.d0 / a
    recip(2, 2) = 1.d0 / b
    recip(3, 3) = 1.d0 / c
    reclng(1) = a
    reclng(2) = b
    reclng(3) = c
    uc_volume = a * b * c
  else
    call cross(ucell(1, 2), ucell(1, 3), u23)
    call cross(ucell(1, 3), ucell(1, 1), u31)
    call cross(ucell(1, 1), ucell(1, 2), u12)
    call dot(ucell(1, 1), u23, uc_volume)
    do i = 1, 3
      recip(i, 1) = u23(i)/uc_volume
      recip(i, 2) = u31(i)/uc_volume
      recip(i, 3) = u12(i)/uc_volume
    end do
    reclng(1) = 1.d0/sqrt(recip(1, 1) * recip(1, 1) + &
                recip(2, 1) * recip(2, 1) + recip(3, 1) * recip(3, 1))
    reclng(2) = 1.d0/sqrt(recip(1, 2) * recip(1, 2) + &
                recip(2, 2) * recip(2, 2) + recip(3, 2) * recip(3, 2))
    reclng(3) = 1.d0/sqrt(recip(1, 3) * recip(1, 3) + &
                recip(2, 3) * recip(2, 3) + recip(3, 3) * recip(3, 3))
  end if

! Interfacial distances given by dot of direct,recip
! sphere is radius of largest sphere inscribed in unit cell.
! The minimum image cutoff must be less than or equal to this.

  if (is_orthog .ne. 0) then
    uc_sphere = 0.5d0 * min(a, b, c)
  else
    uc_sphere = a + b + c
    do i = 1, 3
      call dot(recip(1, i), ucell(1, i), result)
      distance = result * reclng(i)
      if (distance .lt. uc_sphere) uc_sphere = distance
    end do
    uc_sphere = 0.5d0 * uc_sphere
  end if

  write(mdout, '(a,f9.3,/)') &
        '| Largest sphere to fit in unit cell has radius = ', uc_sphere

  ! Check to be sure pairlist cutoff is not too big:

  if (max_cutoff .gt. uc_sphere) then
    write(mdout, '(a,a)') error_hdr, &
      'max pairlist cutoff must be less than unit cell max sphere radius!'
      call mexit(6, 1)
  end if
  
#ifdef CUDA
  call gpu_init_pbc(a, b, c, alpha, beta, gamma, uc_volume, uc_sphere, max_cutoff, pbc_box, reclng, cut_factor, ucell, recip)
#endif

  return

end subroutine init_pbc

#ifdef MPI
!*******************************************************************************
!
! Subroutine:  bcast_pbc
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine bcast_pbc

  use parallel_dat_mod

  implicit none

  integer               :: num_ints, num_reals  ! returned values discarded
  
  call mpi_bcast(is_orthog, pbc_int_cnt, mpi_integer, 0, &
                 pmemd_comm, err_code_mpi)

  call mpi_bcast(recip, pbc_dbl_cnt, mpi_double_precision, 0, &
                 pmemd_comm, err_code_mpi)

  return

end subroutine bcast_pbc
#endif

!*******************************************************************************
!
! Subroutine:  get_fract_crds (in range of 0.0 - +.999...)
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine get_fract_crds(atm_cnt, crd, fraction)

  implicit none

! Formal arguments:

  integer               :: atm_cnt
  double precision      :: crd(3, atm_cnt)
  double precision      :: fraction(3, atm_cnt)

! Local variables:

  double precision      :: f1, f2, f3
  double precision      :: recip_11, recip_22, recip_33
  integer               :: i

! Get fractiontionals.  

  if (is_orthog .ne. 0) then

    recip_11 = recip(1, 1)
    recip_22 = recip(2, 2)
    recip_33 = recip(3, 3)

    do i = 1, atm_cnt
      f1 = recip_11 * crd(1, i)
      f2 = recip_22 * crd(2, i)
      f3 = recip_33 * crd(3, i)
      fraction(1, i) = f1 - dnint(f1) + 0.5d0
      fraction(2, i) = f2 - dnint(f2) + 0.5d0
      fraction(3, i) = f3 - dnint(f3) + 0.5d0
    end do

  else

    do i = 1, atm_cnt

      f1 = crd(1, i) * recip(1, 1) + crd(2, i) * recip(2, 1) + &
           crd(3, i) * recip(3, 1)

      f2 = crd(1, i) * recip(1, 2) + crd(2, i) * recip(2, 2) + &
           crd(3, i) * recip(3, 2)

      f3 = crd(1, i) * recip(1, 3) + crd(2, i) * recip(2, 3) + &
           crd(3, i) * recip(3, 3)

      fraction(1, i) = f1 - dnint(f1) + 0.5d0
      fraction(2, i) = f2 - dnint(f2) + 0.5d0
      fraction(3, i) = f3 - dnint(f3) + 0.5d0

    end do

  end if

  ! We must have fractional coordinates in the range of 0.0 - 0.999...
  ! The algorithm used above will produce fractionals in the range of 0.0 -
  ! 1.0, with some possibility of slight underflow (neg value) due to
  ! rounding error (confirmed). SO we force fractionals to be nonredundant
  ! and within the anticipated range here.

  do i = 1, atm_cnt
    if (fraction(1, i) .lt. 0.d0) fraction(1, i) = fraction(1, i) + 1.d0
    if (fraction(1, i) .ge. 1.d0) fraction(1, i) = fraction(1, i) - 1.d0
    if (fraction(2, i) .lt. 0.d0) fraction(2, i) = fraction(2, i) + 1.d0
    if (fraction(2, i) .ge. 1.d0) fraction(2, i) = fraction(2, i) - 1.d0
    if (fraction(3, i) .lt. 0.d0) fraction(3, i) = fraction(3, i) + 1.d0
    if (fraction(3, i) .ge. 1.d0) fraction(3, i) = fraction(3, i) - 1.d0
  end do

  return

end subroutine get_fract_crds

!*******************************************************************************
!
! All of the particle mesh Ewald code was written and contributed 
! by Tom Darden from the National Institute of Environmental Health 
! Sciences division of the NIH.  Originally written with a modified 
! version of AMBER 3A, the code was updated during the summer of 1994
! to be compatible with AMBER 4.1.
!
! The routines below are used in the particle mesh Ewald code 
! specifically for pressure scaling, volume calculation, imaging, 
! and other things related to the manipulating the periodic box.
!
!*******************************************************************************

!*******************************************************************************
!
! Subroutine:   pressure_scale_crds
!
! Description:  Pressure scaling routine for crds. ONLY used for constant
!               pressure scaling (ntp .gt. 0)!
!
!*******************************************************************************

#ifdef MPI
subroutine pressure_scale_crds(crd, mass, mol_mass_inv, my_mol_lst, mol_com)
#else
subroutine pressure_scale_crds(crd, mass, mol_mass_inv, mol_com)
#endif

  use dynamics_dat_mod
  use dynamics_mod
  use parallel_dat_mod
  use prfs_mod
  use mol_list_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  double precision      :: mass(*)              ! atom mass array.
  double precision      :: mol_mass_inv(*)
#ifdef MPI
  integer               :: my_mol_lst(*)
#endif
  double precision      :: mol_com(3, *)

! Local variables:

  double precision      :: last_com(3)
  double precision      :: new_com(3)
  double precision      :: frac1, frac2, frac3
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
#else
  integer               :: atm_id_lo
  integer               :: atm_id_hi
#endif

! Apply center of molecule based pressure scaling:

#ifdef MPI
  call get_mol_com(my_mol_cnt, crd, mass, mol_mass_inv, my_mol_lst, mol_com)

  do mylist_idx = 1, my_mol_cnt
    mol_id = my_mol_lst(mylist_idx)

! Now get fracs for c.o.m. using old cell params

    last_com(:) = mol_com(:, mol_id)

    frac1 = last_com(1) * last_recip(1, 1) + &
            last_com(2) * last_recip(2, 1) + &
            last_com(3) * last_recip(3, 1)

    frac2 = last_com(1) * last_recip(1, 2) + &
            last_com(2) * last_recip(2, 2) + &
            last_com(3) * last_recip(3, 2)

    frac3 = last_com(1) * last_recip(1, 3) + &
            last_com(2) * last_recip(2, 3) + &
            last_com(3) * last_recip(3, 3)

! Use these with new cell params to get new c.o.m. cartesians:

    new_com(1) = frac1 * ucell(1, 1) + frac2 * ucell(1, 2) + frac3 * ucell(1, 3)
    new_com(2) = frac1 * ucell(2, 1) + frac2 * ucell(2, 2) + frac3 * ucell(2, 3)
    new_com(3) = frac1 * ucell(3, 1) + frac2 * ucell(3, 2) + frac3 * ucell(3, 3)

! Now rigidly translate molecule:

    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      crd(:, atm_id) = crd(:, atm_id) + new_com(:) - last_com(:)
    end do

! Save the new COM:

    mol_com(:, mol_id) = new_com(:)

  end do

  do mylist_idx = 1, my_frag_mol_cnt

    frag_mol_idx = gbl_my_frag_mol_lst(mylist_idx)

    ! Now get fracs for c.o.m. using old cell params

    last_com(:) = mol_com(:, gbl_frag_mols(frag_mol_idx)%mol_idx)

    frac1 = last_com(1) * last_recip(1, 1) + &
            last_com(2) * last_recip(2, 1) + &
            last_com(3) * last_recip(3, 1)

    frac2 = last_com(1) * last_recip(1, 2) + &
            last_com(2) * last_recip(2, 2) + &
            last_com(3) * last_recip(3, 2)

    frac3 = last_com(1) * last_recip(1, 3) + &
            last_com(2) * last_recip(2, 3) + &
            last_com(3) * last_recip(3, 3)

! Use these with new cell params to get new c.o.m. cartesians:

    new_com(1) = frac1 * ucell(1, 1) + frac2 * ucell(1, 2) + frac3 * ucell(1, 3)
    new_com(2) = frac1 * ucell(2, 1) + frac2 * ucell(2, 2) + frac3 * ucell(2, 3)
    new_com(3) = frac1 * ucell(3, 1) + frac2 * ucell(3, 2) + frac3 * ucell(3, 3)

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
            crd(:, atm_id) = crd(:, atm_id) + new_com(:) - last_com(:)
          end do
        end do
      end if

    end do

    ! Save the new COM:

    mol_com(:, gbl_frag_mols(frag_mol_idx)%mol_idx) = new_com(:)

  end do

#else

  call get_mol_com(gbl_mol_cnt, crd, mass, mol_mass_inv, mol_com)

  do mol_id = 1, gbl_mol_cnt

! Now get fracs for c.o.m. using old cell params

    last_com(:) = mol_com(:, mol_id)

    frac1 = last_com(1) * last_recip(1, 1) + &
            last_com(2) * last_recip(2, 1) + &
            last_com(3) * last_recip(3, 1)

    frac2 = last_com(1) * last_recip(1, 2) + &
            last_com(2) * last_recip(2, 2) + &
            last_com(3) * last_recip(3, 2)

    frac3 = last_com(1) * last_recip(1, 3) + &
            last_com(2) * last_recip(2, 3) + &
            last_com(3) * last_recip(3, 3)

! Use these with new cell params to get new c.o.m. cartesians:

    new_com(1) = frac1 * ucell(1, 1) + frac2 * ucell(1, 2) + frac3 * ucell(1, 3)
    new_com(2) = frac1 * ucell(2, 1) + frac2 * ucell(2, 2) + frac3 * ucell(2, 3)
    new_com(3) = frac1 * ucell(3, 1) + frac2 * ucell(3, 2) + frac3 * ucell(3, 3)

! Now rigidly translate molecule:

    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      crd(:, atm_id) = crd(:, atm_id) + new_com(:) - last_com(:)
    end do

! Save the new COM:

    mol_com(:, mol_id) = new_com(:)

  end do

#endif /* MPI */

  return

end subroutine pressure_scale_crds

!*******************************************************************************
!
! Subroutine:   pressure_scale_restraint_crds
!
! Description:  Pressure scaling routine for restraint crds. ONLY used for
!               constant pressure scaling (ntp .gt. 0)!  This routine does not
!               need to keep COM data for the restraint crds, but we use a
!               temporary array in order to have one subroutine that does all
!               COM determination.  We are motivated here by the fact that
!               knowing the COM for all molecules owned is made much more
!               complicated by the new atom divison scheme which is on residue
!               boundaries regardless of whether we are running a constant
!               volume or constant pressure simulation.
!*******************************************************************************

#ifdef MPI
subroutine pressure_scale_restraint_crds(crd, mass, mol_mass_inv, my_mol_lst)
#else
subroutine pressure_scale_restraint_crds(crd, mass, mol_mass_inv)
#endif

  use dynamics_dat_mod
  use dynamics_mod
  use prfs_mod
  use mol_list_mod
  use prmtop_dat_mod
  use parallel_dat_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)
  double precision      :: mass(*)              ! atom mass array.
  double precision      :: mol_mass_inv(*)
#ifdef MPI
  integer               :: my_mol_lst(*)
#endif

! Local variables:

  double precision      :: mol_com(3, gbl_mol_cnt)      ! Temporary COM data.
  double precision      :: last_com(3)
  double precision      :: new_com(3)
  double precision      :: frac1, frac2, frac3
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
#else
  integer               :: atm_id_lo
  integer               :: atm_id_hi
#endif

! Apply center of molecule based pressure scaling:

#ifdef MPI
  call get_mol_com(my_mol_cnt, crd, mass, mol_mass_inv, my_mol_lst, mol_com)

  do mylist_idx = 1, my_mol_cnt
    mol_id = my_mol_lst(mylist_idx)

! Now get fracs for c.o.m. using old cell params

    last_com(:) = mol_com(:, mol_id)

    frac1 = last_com(1) * last_recip(1, 1) + &
            last_com(2) * last_recip(2, 1) + &
            last_com(3) * last_recip(3, 1)

    frac2 = last_com(1) * last_recip(1, 2) + &
            last_com(2) * last_recip(2, 2) + &
            last_com(3) * last_recip(3, 2)

    frac3 = last_com(1) * last_recip(1, 3) + &
            last_com(2) * last_recip(2, 3) + &
            last_com(3) * last_recip(3, 3)

! Use these with new cell params to get new c.o.m. cartesians:

    new_com(1) = frac1 * ucell(1, 1) + frac2 * ucell(1, 2) + frac3 * ucell(1, 3)
    new_com(2) = frac1 * ucell(2, 1) + frac2 * ucell(2, 2) + frac3 * ucell(2, 3)
    new_com(3) = frac1 * ucell(3, 1) + frac2 * ucell(3, 2) + frac3 * ucell(3, 3)

! Now rigidly translate molecule:

    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      crd(:, atm_id) = crd(:, atm_id) + new_com(:) - last_com(:)
    end do
  end do

  do mylist_idx = 1, my_frag_mol_cnt

    frag_mol_idx = gbl_my_frag_mol_lst(mylist_idx)

    ! Now get fracs for c.o.m. using old cell params

    last_com(:) = mol_com(:, gbl_frag_mols(frag_mol_idx)%mol_idx)

    frac1 = last_com(1) * last_recip(1, 1) + &
            last_com(2) * last_recip(2, 1) + &
            last_com(3) * last_recip(3, 1)

    frac2 = last_com(1) * last_recip(1, 2) + &
            last_com(2) * last_recip(2, 2) + &
            last_com(3) * last_recip(3, 2)

    frac3 = last_com(1) * last_recip(1, 3) + &
            last_com(2) * last_recip(2, 3) + &
            last_com(3) * last_recip(3, 3)

! Use these with new cell params to get new c.o.m. cartesians:

    new_com(1) = frac1 * ucell(1, 1) + frac2 * ucell(1, 2) + frac3 * ucell(1, 3)
    new_com(2) = frac1 * ucell(2, 1) + frac2 * ucell(2, 2) + frac3 * ucell(2, 3)
    new_com(3) = frac1 * ucell(3, 1) + frac2 * ucell(3, 2) + frac3 * ucell(3, 3)

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
            crd(:, atm_id) = crd(:, atm_id) + new_com(:) - last_com(:)
          end do
        end do
      end if
    end do

  end do

#else

  call get_mol_com(gbl_mol_cnt, crd, mass, mol_mass_inv, mol_com)

  do mol_id = 1, gbl_mol_cnt

! Now get fracs for c.o.m. using old cell params

    last_com(:) = mol_com(:, mol_id)

    frac1 = last_com(1) * last_recip(1, 1) + &
            last_com(2) * last_recip(2, 1) + &
            last_com(3) * last_recip(3, 1)

    frac2 = last_com(1) * last_recip(1, 2) + &
            last_com(2) * last_recip(2, 2) + &
            last_com(3) * last_recip(3, 2)

    frac3 = last_com(1) * last_recip(1, 3) + &
            last_com(2) * last_recip(2, 3) + &
            last_com(3) * last_recip(3, 3)

! Use these with new cell params to get new c.o.m. cartesians:

    new_com(1) = frac1 * ucell(1, 1) + frac2 * ucell(1, 2) + frac3 * ucell(1, 3)
    new_com(2) = frac1 * ucell(2, 1) + frac2 * ucell(2, 2) + frac3 * ucell(2, 3)
    new_com(3) = frac1 * ucell(3, 1) + frac2 * ucell(3, 2) + frac3 * ucell(3, 3)

! Now rigidly translate molecule:

    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      crd(:, atm_id) = crd(:, atm_id) + new_com(:) - last_com(:)
    end do

  end do

#endif /* MPI */

  return

end subroutine pressure_scale_restraint_crds

!*******************************************************************************
!
! Subroutine:  dot
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine dot(v1, v2, result)

  implicit none

  double precision      :: v1(3), v2(3), result

  result = v1(1) * v2(1) + v1(2) * v2(2) + v1(3) * v2(3)

  return

end subroutine dot

!*******************************************************************************
!
! Subroutine:  cross
!
! Description: <TBS>
!              
!*******************************************************************************

subroutine cross(v1, v2, v12)

  implicit none

! v12 is cross product of v1 and v2

  double precision      :: v1(3), v2(3), v12(3)

  v12(1) = v1(2) * v2(3) - v1(3) * v2(2)
  v12(2) = v1(3) * v2(1) - v1(1) * v2(3)
  v12(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return

end subroutine cross

!*******************************************************************************
!
! Subroutine:  pressure_scale_pbc_data
!
! Description: Scales unit cell uniformly by factor (isotrophically);
!              unit cell angles unchanged.
!              
!*******************************************************************************

subroutine pressure_scale_pbc_data(factor, max_cutoff, verbose)

  use parallel_dat_mod
  use pmemd_lib_mod

  implicit none

! Formal arguments:

  double precision, intent(in)  :: factor(3)
  double precision, intent(in)  :: max_cutoff
  integer, intent(in)           :: verbose

! Local variables:

  double precision      :: distance
  double precision      :: facinv(3)
  double precision      :: result
  integer               :: i, j


  last_recip(:,:) = recip(:,:)

  uc_volume = uc_volume * factor(1) * factor(2) * factor(3)
  pbc_box(:) = factor(:) * pbc_box(:)

  if (verbose .eq. 1 .and. master) then
     write(mdout, '(a, 4f12.3)') ' a,b,c,volume now equal to ', &
                                 pbc_box(1), pbc_box(2), pbc_box(3), uc_volume
  end if

  do i = 1, 3
    facinv(i) = 1.d0 / factor(i)
    reclng(i) = factor(i) * reclng(i)
  end do

  do j = 1, 3
    do i = 1, 3
      ucell(i, j) = factor(i) * ucell(i, j)
      recip(i, j) = facinv(i) * recip(i, j)
    end do
  end do

  uc_sphere = pbc_box(1) + pbc_box(2) + pbc_box(3)

  do i = 1, 3
    call dot(recip(1, i), ucell(1, i), result)
    distance = result * reclng(i)
    if (distance .lt. uc_sphere) uc_sphere = distance
  end do

  uc_sphere = 0.5d0 * uc_sphere

  ! Check to be sure pairlist cutoff is not too big:

  if (max_cutoff .gt. uc_sphere) then
    write(mdout, '(a,a)') error_hdr, &
      'max pairlist cutoff must be less than unit cell max sphere radius!'
      call mexit(6, 1)
  end if

  return

end subroutine pressure_scale_pbc_data

!*******************************************************************************
!
! Subroutine:  wrap_molecules
!
! Description:
!
! Wrap the molecules/coordinates across the periodic box.
! Geometric center of each molecule is checked to see if
! it is within the unit cell or not.
!              
!*******************************************************************************

subroutine wrap_molecules(crd)

  use mol_list_mod

  implicit none

  ! Formal arguments:

  double precision      :: crd(3, *)

  ! Local variables:

  double precision      :: tran(3), f1, f2, f3, g1, g2, g3
  integer               :: i    ! mol atm list idx
  integer               :: atm_id, mol_id
  integer               :: mol_atm_cnt, offset

  do mol_id = 1, gbl_mol_cnt
    f1 = 0.0d0
    f2 = 0.0d0
    f3 = 0.0d0
    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset
    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      f1 = f1 + crd(1, atm_id) * recip(1, 1) + crd(2, atm_id) * recip(2, 1) + &
                crd(3, atm_id) * recip(3, 1)

      f2 = f2 + crd(1, atm_id) * recip(1, 2) + crd(2, atm_id) * recip(2, 2) + &
                crd(3, atm_id) * recip(3, 2)

      f3 = f3 + crd(1, atm_id) * recip(1, 3) + crd(2, atm_id) * recip(2, 3) + &
                crd(3, atm_id) * recip(3, 3)
    end do

#if 0
    f1 = f1/mol_atm_cnt
    f2 = f2/mol_atm_cnt
    f3 = f3/mol_atm_cnt

    g1 = f1
    if (f1 .lt. 0.d0) g1 = f1 + 1.d0
    if (f1 .ge. 1.d0) g1 = f1 - 1.d0

    g2 = f2
    if (f2 .lt. 0.d0) g2 = f2 + 1.d0
    if (f2 .ge. 1.d0) g2 = f2 - 1.d0

    g3 = f3
    if (f3 .lt. 0.d0) g3 = f3 + 1.d0
    if (f3 .ge. 1.d0) g3 = f3 - 1.d0
#else
    f1 = f1/mol_atm_cnt - 0.5d0
    f2 = f2/mol_atm_cnt - 0.5d0
    f3 = f3/mol_atm_cnt - 0.5d0

    g1 = f1 - anint(f1)
    g2 = f2 - anint(f2)
    g3 = f3 - anint(f3)
#endif

    if (f1 .ne. g1 .or. f2 .ne. g2 .or. f3 .ne. g3) then
      tran(1) = (g1 - f1) * ucell(1, 1) + (g2 - f2) * ucell(1, 2) + &
                (g3 - f3) * ucell(1, 3)
      tran(2) = (g1 - f1) * ucell(2, 1) + (g2 - f2) * ucell(2, 2) + &
                (g3 - f3) * ucell(2, 3)
      tran(3) = (g1 - f1) * ucell(3, 1) + (g2 - f2) * ucell(3, 2) + &
                (g3 - f3) * ucell(3, 3)

      do i = offset + 1, offset + mol_atm_cnt
        atm_id = gbl_mol_atms_lists(i)
        crd(1, atm_id) = crd(1, atm_id) + tran(1)
        crd(2, atm_id) = crd(2, atm_id) + tran(2)
        crd(3, atm_id) = crd(3, atm_id) + tran(3)
      end do

      if (mol_id .eq. 1) write(mdout,'(a,3f15.5)') 'wrapping first mol.:', &
                          tran(1), tran(2), tran(3)
    end if

  end do

  return

end subroutine wrap_molecules

!*******************************************************************************
!
! Subroutine:  wrap_to
!
! Description:
!
! The trunf. oct. has:
!   * Center at (0,0,0) where the corner of the triclinic cell is. 
!   * One hex face with the x axis for a normal.
!                  Face is perpendicular to the x axis and 
!                  1/2 box(1) away from the origin.
!   * One hex face perp to xy plane, second edge vector of the
!                  triclinic cell is its normal, 109 degrees from
!                  the x axis in xy plane (-x,+y quadrant).
!
! Approach to reconstruct the t.o. is to rotate the coordinates
!    to put the hex faces in the (+-1,+-1,+-1) normal directions, and
!    the square (diamond) faces perp to xyz axes. This is 3 rotations:
!    We did +45 around z, +(90-tetra/2) around y, +90 around x to get
!       the to oriented for the triclinic cell.
!    Now we do the opposite: -90(x), -(90-tetra/2)(y), -45(z)
!       to reproduce the original orientation, map coords into 
!       a t.o. centered at origin, then do the rotations again
!       to put it back in the orientation that matches the restrt.
!              
!*******************************************************************************

subroutine wrap_to(crd)

  use mol_list_mod

  implicit none

! Formal arguments:

  double precision      :: crd(3, *)

! Local variables:

  double precision      :: cx, cy, cz, x0, y0, z0, x, y, z, xt, yt, zt
  double precision      :: phi, cos1, sin1, cos2, sin2
  double precision      :: tobox, tobinv, facecoord
  double precision      :: t11, t12, t13, t21, t22, t23, t31, t32, t33
  integer               :: i                    ! mol atm list idx
  integer               :: atm_id, mol_id
  integer               :: mol_atm_cnt, offset

  facecoord = pbc_box(1) / (2.d0 * sqrt(3.d0))
  tobox = 2.d0 * pbc_box(1) / sqrt(3.d0)
  tobinv = 1.d0 / tobox
  phi = PI / 4.d0
  cos1 = cos(phi)
  sin1 = sin(phi)
  cos2 = sqrt(2.d0)/sqrt(3.d0)
  sin2 = 1. / sqrt(3.d0)
  t11 = cos2 * cos1
  t12 = - cos2 * sin1
  t13 = - sin2
  t21 = - sin2 * cos1
  t22 = sin2 * sin1
  t23 = - cos2
  t31 = sin1
  t32 = cos1
  t33 = 0

  do mol_id = 1, gbl_mol_cnt

! Calculate center of geometry of molecule:

    cx = 0.d0
    cy = 0.d0
    cz = 0.d0
    mol_atm_cnt = gbl_mol_atms_listdata(mol_id)%cnt
    offset = gbl_mol_atms_listdata(mol_id)%offset

    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      cx = crd(1, atm_id) + cx
      cy = crd(2, atm_id) + cy
      cz = crd(3, atm_id) + cz
    end do

    cx = cx / mol_atm_cnt
    cy = cy / mol_atm_cnt
    cz = cz / mol_atm_cnt

! Rotate:

    x0 = cx * t11 + cy * t21 + cz * t31
    y0 = cx * t12 + cy * t22 + cz * t32
    z0 = cx * t13 + cy * t23 + cz * t33

! First map into cube of size 2 * pbc_box(1)/sqrt(3):

    xt = dnint(x0 * tobinv)
    x = x0 - xt * tobox
    yt = dnint(y0 * tobinv)
    y = y0 - yt * tobox
    zt = dnint(z0 * tobinv)
    z = z0 - zt * tobox

! Wrap molecules external to diag faces:

    xt = abs(x)
    yt = abs(y)
    zt = abs(z)
    if (xt + yt + zt .gt. 3. * facecoord) then

      if (x .gt. 0.d0) then
        x = x - 2. * facecoord
      else
        x = x + 2. * facecoord
      end if

      if (y .gt. 0.d0) then
        y = y - 2. * facecoord
      else
        y = y + 2. * facecoord
      end if

      if (z .gt. 0.d0) then
        z = z - 2. * facecoord
      else
        z = z + 2. * facecoord
      end if

    end if

! Get the translation in the rotated space for this molecules c-o-geom:

    xt = x - x0
    yt = y - y0
    zt = z - z0

! Rotate:

    cx = xt * t11 + yt * t12 + zt * t13
    cy = xt * t21 + yt * t22 + zt * t23
    cz = xt * t31 + yt * t32 + zt * t33

! Now move the molecule:

    do i = offset + 1, offset + mol_atm_cnt
      atm_id = gbl_mol_atms_lists(i)
      crd(1, atm_id) = crd(1, atm_id) + cx
      crd(2, atm_id) = crd(2, atm_id) + cy
      crd(3, atm_id) = crd(3, atm_id) + cz
    end do

  end do

  return

end subroutine wrap_to

end module pbc_mod
