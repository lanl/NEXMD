! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"

! ---------------------------------------------------------
! Ross Walker (SDSC) & Gustavo Seabra (UFL) 2008.
! with modifications by Andreas Goetz (SDSC) 2011.
!
! Routines for carrying out QM/MM simulations with a variable
! QM solvent selection based on the nearest N QM solvent molecules.
!
! 1) qmmm_vsolv_setup
!       This routine is called to setup the variable solvent
!       code. It does things like allocation and also
!       stores the ids of the solvent molecules etc.
!
! 2) qmmm_vsolv_update_iqmatoms
!       This routine updates iqmatoms to store the new atom
!       numbers for the new solvent.
!       Called from qmmm_vsolv_setup
!               and qmmm_vsolv_update
!
! 3) qmmm_vsolv_identify_nearest_solvent
!       This routine locates the nearest n solvent to the QM
!       region and stores the residue id of the solvent.
!
! 4) qmmm_vsolv_update
!       This routine is responsible for doing the actual update
!       of the QM region for the variable solvent. Essentially
!       it calls the routines necessary to find the nearest
!       solvent molecules, then calls the routine to update
!       iqmatoms, updates the charges, bonds, angles, dihedrals
!       etc.
! 
!
! Note about iqmatoms.
!
! Iqmatoms is normally sorted numerically, however, in the 
! general case for variable solvent this would mean that lots
! of additional arrays would need to be updated when the variable
! solvents change - this is true for the case where solvent residue
! numbers might bracket the fixed QM region residue numbers. In most
! cases this would not cause a problem since the solvent (typically WAT)
! is always last in the prmtop file and so would always be after the fixed
! QM region. However, in the general case this does not apply so in the case
! of variable QM solvent the iqmatoms array is as follows.
!
! fixed_QM_1
! fixed_QM_2
! fixed_QM_3
! ... to nquant - sorted numerically after reading qmmask.
! ...
! variable_QM
! variable_QM
! variable_QM
! ... to nearest_qm_solvent * natom_solv_res
! ... NOT sorted - but filled in order of distance by the identify
! ... routine - but may have lower numbers than the
! ... fixed_QM atoms.
! nlink_1
! nlink_2
! ... to nlink - contains atom number of MM link pair atom
! ...            I believe these are sorted numerically within
! ...            themselves.
!
! ---------------------------------------------------------

subroutine qmmm_vsolv_setup(nquant, max_quantum_atoms, &
                            iqmatoms, nres, lbres, res_pointers,&
                            natom_res, natom)

  use qmmm_module, only : qmmm_nml, qmmm_vsolv
  use molecule, only : mol_info

  implicit none

#ifdef MPI
  ! REB: Include parallel.h for adaptive QM/MM (DAS) with multisander
#include "parallel.h"
#endif

  ! Passed in
  integer, intent(inout) :: nquant
  integer, intent(in)    :: max_quantum_atoms, nres
  integer, intent(inout) :: iqmatoms(max_quantum_atoms)
  character(len=4), intent(in) :: lbres(nres)
  integer, intent(in) :: res_pointers(nres+1)
  integer, intent(in) :: natom_res(nres)
  integer, intent(in) :: natom
  
  ! Local
  integer :: i,j, ier, res_size, loop_count, iatom
  logical :: first_hit

#ifdef MPI
  if (qmmm_nml%vsolv > 1) then
     ! adaptive QM/MM based on vsolv code
     ! Change the number of nearest_qm_solvents in each multisander group
     qmmm_vsolv%nearest_qm_solvent = qmmm_vsolv%nearest_qm_solvent + masterrank 
  end if
#endif

  ! ATTENTION, ONLY REQUIRED FOR ADAPTIVE QM/MM CODE
  qmmm_vsolv%recalculate = .false.

  if ( qmmm_vsolv%debug ) then
     write(6,'(a)') '>>>>> entered qmmm_vsolv_setup'
     call flush(6)
  end if

  ! Store the user specified number of quantum atoms so we don't have
  ! to search for it again. Since we always sort the QM atoms so this
  ! represents the end of the non variable QM region. Also store 
  ! the original iqmatoms array since we need that because the main
  ! iqmatoms gets sorted.
  qmmm_vsolv%fixed_nquant = nquant
  allocate ( qmmm_vsolv%fixed_iqmatoms(nquant), stat=ier)
  REQUIRE(ier==0)
  qmmm_vsolv%fixed_iqmatoms(1:nquant) = iqmatoms(1:nquant)

  ! check if qm_center_atom_id is small 
  ! or equal to the number of fixed QM atoms
  if (qmmm_vsolv%qm_center_atom_id > nquant) then
     call sander_bomb('qm2_setup_vsolv', &
                      'qm_center_atom_id has to be smaller or ', &
                      'equal to the total number of fixed QM atoms.')
  end if

  !Next we need to work out how many atoms nearest_qm_solvent x atoms per solvent resname represents.

  !1) test 1 is to find out how many residues of the solvent resname specified are in the system.
  qmmm_vsolv%nsolv_res = 0
  first_hit = .true.
  do i = 1, nres
    if (trim(lbres(i))==trim(qmmm_vsolv%nearest_qm_solvent_resname)) then
       qmmm_vsolv%nsolv_res = qmmm_vsolv%nsolv_res + 1
       !store the number of atoms this residue corresponds to
       !Get the start and end location for this residue.
       res_size=res_pointers(i+1) - res_pointers(i)
       if (first_hit) then
         first_hit = .false.
         qmmm_vsolv%natom_solv_res = res_size
       else
         if (qmmm_vsolv%natom_solv_res /=  res_size) then
           write(6,'(a,I6,a,I6,a)') 'QMMM: Error - Residue ',i,' has a different number of atoms (',res_size,')'
           write(6,'(a,i6,a)') 'QMMM: Error - than the first selected solvent residue which has ', &
                          qmmm_vsolv%natom_solv_res ,' atoms.'
           call sander_bomb('qm2_setup_vsolv', &
                            'Selected solvent residues for variable QM region do not have constant atom count.',&
                            'It is a requirement that all selected solvent molecules have constant atom counts.')
         end if
       end if

       ! check if nearest_qm_solvent_center_id is small
       ! or equal to the number of 
       if (qmmm_vsolv%nearest_qm_solvent_center_id > res_size) then
          write(6,'(a)') 'QMMM: Error - nearest_qm_solvent_center_id'
          write(6,'(16x,a)') '> number of atoms in the solvent residue'
          call sander_bomb('qm2_setup_vsolv', &
                            'nearest_qm_solvent_center_id has to be small or', &
                            'equal to the number of atoms in the residue.')
       end if
    end if
  end do

  !TEST - Require at least as many solvent molecules in the system than the user asked to be QM.
  if (qmmm_vsolv%nsolv_res <  qmmm_vsolv%nearest_qm_solvent) then
    write(6,'(a,i8,a,a,/,a,i6,a)') 'QMMM: Error - Found only ',qmmm_vsolv%nsolv_res,' solvent residues of name ', &
                    qmmm_vsolv%nearest_qm_solvent_resname, 'QMMM: Error - but you requested ',qmmm_vsolv%nearest_qm_solvent, &
                    ' solvent residues to be in the QM region.'
    call sander_bomb('qm2_setup_vsolv','The requested number of solvent molecules in the QM region',&
                     'must be the same or smaller than the total number of solvent molecules in the simulation.')
  end if

  !Now we know how many residues there are to consider and the number of atoms per residue we will store the
  !starting atom of each solvent residue - this is effectively a rebasing of the residue pointers.
  allocate(qmmm_vsolv%solvent_pointers(qmmm_vsolv%nsolv_res),stat = ier)
  REQUIRE(ier==0)
  
  !We have to repeat the above loop.
  loop_count = 0
  do i = 1, nres
    if (trim(lbres(i))==trim(qmmm_vsolv%nearest_qm_solvent_resname)) then
      loop_count = loop_count+1
      qmmm_vsolv%solvent_pointers(loop_count) = res_pointers(i)
    end if
  end do

  REQUIRE(loop_count == qmmm_vsolv%nsolv_res)

  !Next test - check to make sure that the fixed QM region does not contain any of the
  !            atoms having the variable solvent residue id.
  do i=1,qmmm_vsolv%fixed_nquant
     iatom = qmmm_vsolv%fixed_iqmatoms(i)
     !check if this atoms residue ID matches the one we are considering to be solvent.
     if (trim(lbres(mol_info%atom_to_resid_map(iatom))) == trim(qmmm_vsolv%nearest_qm_solvent_resname)) then
        write(6,'(a,i8,a,i8,a,a)') 'QMMM: Error - QM atom id ',iatom, &
                                   '(residue ',mol_info%atom_to_resid_map(iatom), &
                                   ') has resname = ',lbres(mol_info%atom_to_resid_map(iatom))
        write(6,'(a)') 'QMMM: Error - which matches the specified nearest_qm_solvent_resname.'
        write(6,'(a)') 'QMMM: Error - This is NOT permitted.'
        call sander_bomb('qm2_setup_vsolv','Found solvent residue in fixed QM region.', &
                         'Variable QM solvent residues cannot be specified as permanently QM.')

     end if
  end do

  !Allocate an array that stores the nearest_solvent_pointers = the first atom number
  !of the nearest_qm_solvents to the QM region.
  !
  !We are allocating the array with size nearest_qm_solvent+1.
  !The last value in the array points the first MM solvent 
  !and it is going to be used in the qmmm_adaptive_module.
  allocate(qmmm_vsolv%nearest_solvent_pointers(qmmm_vsolv%nearest_qm_solvent+1), stat=ier)
  REQUIRE(ier==0)

  !It would now make sense to call the routine that identifies the nearest solvents
  !but we can't do this here because we don't have the coordinates read at this point.
  !so for the time being we just set iqmatoms to point to the first N solvent residues.
  do i = 1, qmmm_vsolv%nearest_qm_solvent + 1
    qmmm_vsolv%nearest_solvent_pointers(i)=qmmm_vsolv%solvent_pointers(i)
  end do

  !Allocate an array for the distances of nearest solvent residues 
  !from the QM center
  allocate(qmmm_vsolv%nearest_solvent_distances(qmmm_vsolv%nearest_qm_solvent+1),stat=ier)
  REQUIRE(ier==0)

  !Assign the initial distances as zeros
  qmmm_vsolv%nearest_solvent_distances(:)=0.0d0

  !Now we need to update iqmatoms and nquant.
  nquant = nquant + qmmm_vsolv%nearest_qm_solvent * qmmm_vsolv%natom_solv_res

  call qmmm_vsolv_update_iqmatoms(qmmm_vsolv%fixed_nquant, &
                                  qmmm_vsolv%nearest_qm_solvent, &
                                  qmmm_vsolv%nearest_solvent_pointers, &
                                  qmmm_vsolv%natom_solv_res, nquant, iqmatoms, &
                                  natom, qmmm_vsolv%verbosity, qmmm_vsolv%debug)

  !Print an information message to the output file if the nearest_qm_solvent routine is active.
  write(6,'(a)') ' '
  write(6,'(a)')      'QMMM:         Variable QM Solvent Region is Active'
  write(6,'(a)')      'QMMM: ------------------------------------------------------'
  write(6,'(a,a)')    'QMMM:             Residue name of solvent molecules : ',trim(qmmm_vsolv%nearest_qm_solvent_resname)
  write(6,'(a,i6)')   'QMMM:                    Atoms per solvent molecule : ',qmmm_vsolv%natom_solv_res
  write(6,'(a,i6)')   'QMMM: Total number of solvent molecules to consider : ',qmmm_vsolv%nsolv_res
  write(6,'(a,i6)')   'QMMM:                      Atoms in fixed QM region : ',qmmm_vsolv%fixed_nquant
  write(6,'(a,i6)')   'QMMM:           Atoms in variable QM solvent region : ',nquant - qmmm_vsolv%fixed_nquant
  write(6,'(a,i6)')   'QMMM:                      Total atoms in QM region : ',nquant
  write(6,'(a,i6,a)') 'QMMM:    QM Solvent region update frequency (steps) : ',qmmm_vsolv%nearest_qm_solvent_fq
  write(6,'(a)')      'QMMM: ------------------------------------------------------'
  write(6,'(a)') ' '

  if ( qmmm_vsolv%debug ) then
     write(6,'(a)') '<<<<< leaving qmmm_vsolv_setup'
     call flush(6)
  end if

end subroutine qmmm_vsolv_setup

! Update iqmatoms.
subroutine qmmm_vsolv_update_iqmatoms(fixed_nquant, nearest_qm_solvent, &
                         nearest_solvent_pointers, natom_solv_res, nquant, iqmatoms, &
                         natom, verbosity, debug)

  implicit none

  ! passed in
  integer, intent(in) :: fixed_nquant       ! Number of permanent QM atoms.
  integer, intent(in) :: nearest_qm_solvent ! Number of solvent molecules to include in QM region.
  integer, intent(in) :: nearest_solvent_pointers(nearest_qm_solvent+1) !First atom of each of the nearest
                                                                      !solvent residues.
  integer, intent(in) :: natom_solv_res ! Number of atoms per solvent residue.
  integer, intent(in) :: nquant         ! Total number of quantum atoms both fixed and variable.
  integer, intent(out) :: iqmatoms(nquant)
  integer, intent(in) :: natom
  integer, intent(in) :: verbosity
  logical, intent(in) :: debug

  ! local
  integer :: loop_count, iatom, i, j

  if ( debug ) then
     write(6,'(a)') '>>>>> entered qmmm_vsolv_update_iqmatoms'
     call flush(6)
  end if

  ! Note: iqmatoms is nquant_nlink here with the nlink entries containing the MM link pair
  !       atoms. However, we are only changing the solvent here so we can just do things over
  !       the nquant atoms.
  
  ! iqmatoms array is only sorted over the fixed_nquant - these always appear first so we
  ! only need to change the variable atoms which occur from fixed_nquant+1 upwards.

  loop_count = fixed_nquant
  do i = 1,nearest_qm_solvent
    iatom = nearest_solvent_pointers(i)
    do j = 1, natom_solv_res
      loop_count=loop_count+1
      iqmatoms(loop_count)= iatom
      iatom = iatom + 1
    end do
  end do

  if (verbosity>1) then
     write (6,'(a)') 'QMMM:'
     write (6,'(a)') 'QMMM: New iqmatoms listing:'
     do i = 1, nquant 
       write (6,'(a,i6,a,i6)') 'QMMM: ',i,'     ',iqmatoms(i)
     end do
     write (6,'(a)') 'QMMM:'
  end if

  if ( debug ) then
     write(6,'(a)') '<<<<< leaving qmmm_vsolv_update_iqmatoms'
     call flush(6)
  end if

end subroutine qmmm_vsolv_update_iqmatoms

subroutine qmmm_vsolv_identify_nearest_solvent(natom, nres, natom_res, unimaged_coords, atomic_masses, &
     qmmm_vsolv, natom_3_scratch)

!By Ross Walker (SDSC) and Gustavo Seabra (UFL), 2008.
!
!This routine finds the nearest solvent molecules to the QM
!region (excluding link atoms) and it returns an array
!containing the first atom number of the nearest_qm_solvent
!residues.
!
!Note: The algorithm used here is certainly not the most 
!      efficient and this could probably benefit later
!      on from being replaced with a better algorithm.
!
!Passed in
!                   natom = Total number of atoms in the simulation.
!                    nres = Total number of residues in the simulation.
!               natom_res = array containing the number of atoms in each residue
!                           (all residues) == ix(i70)
!         unimaged_coords = The unimaged MM coordinate array (natom long)
!  contained in qmmm_vsolv_type:
!            fixed_nquant = Number of fixed quantum atoms.
!          fixed_iqmatoms = Atom id list of fixed quantum atoms.
!               nsolv_res = Total number of solvent molecules to consider.
!        solvent_pointers = First atom id of each solvent residue we are considering.
!      nearest_qm_solvent = Number of residues to return that are closest to the QM region.
!nearest_solvent_pointers = This gets returned containing the first atom of each residue
!                           that is closest to the fixed QM region
!nearest_solvent_distances = Distances of nearest solvent residues from 
!                            the QM center
!          natom_solv_res = number of atoms in each solvent residue.
!         natom_3_scratch = 3*natom scratch array.

  use molecule, only : mol_info
  use qmmm_vsolv_module, only : qmmm_vsolv_type
  implicit none

  ! Passed in
  integer, intent(in)  :: natom, nres
  integer, intent(in)  :: natom_res(nres) 
  _REAL_,  intent(in)  :: unimaged_coords(3,natom)
  _REAL_,  intent(in)  :: atomic_masses(natom)
  _REAL_,  intent(out) :: natom_3_scratch(3,natom)
  type(qmmm_vsolv_type), intent(inout) :: qmmm_vsolv

  ! Local
  integer :: i,j,k,res_start,atomi,atomj,ier, current_resid
  integer :: n_qmc_points
  _REAL_  :: residue_distances(qmmm_vsolv%nsolv_res) !should be small enough not to blow the stack.
  _REAL_  :: min_R2, current_ceiling, current_distance
  _REAL_  :: vec1(3), vec2(3), vec3(3)
  _REAL_  :: box_center(3)
  _REAL_  :: vec_s_com(3)
  _REAL_  :: mass_atom, mass_solv_res
  _REAL_, parameter :: MAX_DISTANCE = 99999.0d0
  logical :: include_residue
 
  !Simple (unoptimized algorithm)
  !Loop 1 to nsolv_res - loop over solvent residues.
  !  Loop 1 to fixed_nquant - loop over all fixed quantum atoms.
  !    Loop 1 to natom res - loop over all atoms in the residue.
  !
  !      calculate min R2 between all atoms of this residue and
  !      all atoms of fixed quantum atoms.
  !
  !    end loop
  !  end loop
  !  store this minimum distance for this nsolv_res.

  if ( qmmm_vsolv%debug ) then
     write(6,'(a)') '>>>>> entered qmmm_vsolv_identify_nearest_solvent'
     call flush(6)
  end if

  if (qmmm_vsolv%verbosity > 0) write (6,'(a)') 'QMMM: Updating nearest solvent list.'

  ! ---------------------
  ! REIMAGING COORDINATES
  ! ---------------------
  ! 
  ! Reimage the system around the fixed qm zone.
  ! This is done only to get the minimum distance, so 
  ! we use a temporary coordinate array instead of the 
  ! real 'x'
  !
  ! NOTE: Could use a "minimum image" calculation instead?

  ! 1. Copy the coordinates to the temporary array
  natom_3_scratch(1:3,1:natom) = unimaged_coords(1:3,1:natom)

  ! 2. Image those coordinates  
  call iwrap2(qmmm_vsolv%fixed_nquant, qmmm_vsolv%fixed_iqmatoms, natom_3_scratch, &
              box_center)

  if ( qmmm_vsolv%debug ) then
     write(6,'(a)') ' unimaged coordinates'
     do i = 1, natom
        write(6,'(i7,3f12.6)') i, unimaged_coords(1:3,i)
     end do

     write(6,'(a)') ' re-imaged coordinates'
     do i = 1, natom
        write(6,'(i7,3f12.6)') i, natom_3_scratch(1:3,i)
     end do
  end if
  
  ! -------------------------
  ! END REIMAGING COORDINATES
  ! -------------------------

  if ( qmmm_vsolv%qm_center_atom_id == 0 ) then
     n_qmc_points = qmmm_vsolv%fixed_nquant
  else
     n_qmc_points = 1
  end if


  ! total mass of the solvent molecule
  mass_solv_res = 0.0d0
  res_start = qmmm_vsolv%solvent_pointers(1)

  do i = 1, qmmm_vsolv%natom_solv_res
     atomi = res_start + (i-1)
     mass_solv_res = mass_solv_res + mol_info%atom_mass(atomi)
  end do

  ! calculate distances or solvent molecules 
  do i = 1, qmmm_vsolv%nsolv_res

    res_start = qmmm_vsolv%solvent_pointers(i)
    min_R2 = MAX_DISTANCE

    atomi = res_start

    do j=1,n_qmc_points

      ! no qm center - closest solvent molecules from ANY qm atoms will be chosen
      if ( qmmm_vsolv%qm_center_atom_id == 0 ) then
         vec1(1:3) = natom_3_scratch(1:3,qmmm_vsolv%fixed_iqmatoms(j))

      ! qm center - center of mass of fixed QM region
      else if ( qmmm_vsolv%qm_center_atom_id == -1 ) then
         vec1(1:3) = box_center(1:3)

      ! qm center - atom id given in "qm_center_atom_id"
      else if ( qmmm_vsolv%qm_center_atom_id > 0 ) then
         vec1(1:3) = natom_3_scratch(1:3,qmmm_vsolv%qm_center_atom_id)
      end if

      ! no solvent center - ANY closest solvent atoms from qm center will be chosen
      if ( qmmm_vsolv%nearest_qm_solvent_center_id == 0 ) then

         do k = 1, qmmm_vsolv%natom_solv_res 
            atomi = res_start+(k-1)
            vec2(1:3) = natom_3_scratch(1:3,atomi)
            vec3(1:3) = vec2(1:3) - vec1(1:3)
            vec3(1:3) = vec3(1:3) * vec3(1:3)
            min_R2 = min(sum(vec3),min_R2)
         end do

      else

         ! solvent center - center of mass of the solvent molecule
         if ( qmmm_vsolv%nearest_qm_solvent_center_id == -1 ) then

            vec_s_com(1:3) = 0.0d0

            do k = 1, qmmm_vsolv%natom_solv_res
               atomi = res_start+(k-1)
               mass_atom = mol_info%atom_mass(atomi)
               vec_s_com(1:3) = vec_s_com(1:3) + natom_3_scratch(1:3,atomi)*mass_atom
            end do

            vec2(1:3) = vec_s_com(1:3) / mass_solv_res

         ! solvent center - 1st atom in solvent residue
         else if ( qmmm_vsolv%nearest_qm_solvent_center_id > 0 ) then
            atomi = res_start + qmmm_vsolv%nearest_qm_solvent_center_id - 1 
            vec2(1:3) = natom_3_scratch(1:3,atomi)
         end if

         vec3(1:3) = vec2(1:3) - vec1(1:3)
         vec3(1:3) = vec3(1:3) * vec3(1:3)
         min_R2 = min(sum(vec3),min_R2)

      end if

    end do

    ! We now have the minimum R2 distance between this residue and the
    ! fixed quantum region. Store this.
    residue_distances(i) = min_R2

  end do 

  ! Next find the N smallest values in the residue_distances array and
  ! store their residue id values.
  ! current_floor = 0.0d0
  do i = 1, qmmm_vsolv%nearest_qm_solvent + 1
    
    current_ceiling = MAX_DISTANCE

    do j = 1, qmmm_vsolv%nsolv_res

      current_distance = residue_distances(j)
      if (current_distance < current_ceiling) then

        ! This distance is the smallest we have found in this loop.
        ! Check that it is greater than the smallest we found
        ! last time. One way to do this is to use a floor as the
        ! commented out code implies but this does not work if two
        ! residues happen to be exactly the same distance apart. Hence
        ! we have to look over what we have found so far looking to make
        ! sure we haven't already recorded this residue.
        ! if (current_distance > current_floor) then
        !   current_ceiling = current_distance 
        !   current_resid = j
        ! end if
        include_residue = .true.
        already_included: do k = 1, i-1
          if (qmmm_vsolv%nearest_solvent_pointers(k) == qmmm_vsolv%solvent_pointers(j)) then
            !we have already included this residue.
            include_residue = .false.
            exit already_included
          end if
        end do already_included
        ! if we came out with include_residue still equal to true then
        ! we have not already included this residue
        if (include_residue) then
          current_ceiling = current_distance 
          current_resid = j
        end if
      end if

    end do

    ! When we get here we should have the resid corresponding to the 
    ! smallest R2 that is greater than all smaller ones we have found to date.
    qmmm_vsolv%nearest_solvent_pointers(i) = qmmm_vsolv%solvent_pointers(current_resid)
    qmmm_vsolv%nearest_solvent_distances(i) = residue_distances(current_resid)
    ! current_floor = current_ceiling 
  end do

  if (qmmm_vsolv%verbosity > 0) then
     write (6,'(a)') 'QMMM: New solvent residue list is:'
     write (6,'(a)') 'QMMM: Nearest ID    Res ID     Distance(A)'
     do i = 1, qmmm_vsolv%nearest_qm_solvent
       write (6,'(a,i6,a,i6,8x,f8.3)') 'QMMM:     ',i,'    ',&
              mol_info%atom_to_resid_map(qmmm_vsolv%nearest_solvent_pointers(i)),&
              DSQRT(qmmm_vsolv%nearest_solvent_distances(i)) 
     end do
     write (6,'(a)') 'QMMM:'
  end if
 

  if ( qmmm_vsolv%debug ) then
     write(6,'(a)') '<<<<< leaving qmmm_vsolv_identify_nearest_solvent'
     call flush(6)
  end if

end subroutine qmmm_vsolv_identify_nearest_solvent

subroutine qmmm_vsolv_update(nstep, natom, nres, natom_res, charge, unimaged_coords, atomic_masses, &
                             nbonh, nbona, ntheth, ntheta, nphih, nphia, &
                             ix_iibh, ix_ijbh, ix_iicbh, &
                             ix_iiba, ix_ijba, ix_iicba, &
                             ix_i24, ix_i26, ix_i28, ix_i30, &
                             ix_i32, ix_i34, ix_i36, ix_i38, &
                             ix_i40, ix_i42, ix_i44, ix_i46, ix_i48, &
                             ix_i50, ix_i52, ix_i54, ix_i56, ix_i58, ix_ibellygp)

  use qmmm_module, only : qmmm_struct, qmmm_nml, qmmm_vsolv, qmmm_scratch, qmmm_mpi
  use qmmm_vsolv_module, only : print
  use parms, only : numbnd

  implicit none

  ! passed in
  integer, intent(in) :: nstep     ! Current MD step
  integer, intent(in) :: natom     ! Number of total atoms in the system from prmtop (excluding link atoms)
  integer, intent(in) :: nres      ! Number of total residues in the system from prmtop
  integer, intent(in) :: natom_res(nres) ! array containing the number of atoms in each residue [ix(i70)]
  _REAL_, intent(inout) :: charge(natom) ! Main natom long amber charge array [xx(L15)]
  _REAL_, intent(in) ::  unimaged_coords(3,natom) !The unimaged MM coordinate array (3,natom long)
  _REAL_, intent(in) ::  atomic_masses(natom)     ! Mass of atoms.
  integer, intent(out) :: nbonh, nbona, ntheth, ntheta, nphih, nphia
  integer, intent(out) :: ix_iibh(qmmm_vsolv%nbonh), ix_ijbh(qmmm_vsolv%nbonh), ix_iicbh(qmmm_vsolv%nbonh)
  integer, intent(out) :: ix_iiba(qmmm_vsolv%nbona), ix_ijba(qmmm_vsolv%nbona), ix_iicba(qmmm_vsolv%nbona)
  integer, intent(out) :: ix_i24(qmmm_vsolv%ntheth), ix_i26(qmmm_vsolv%ntheth), ix_i28(qmmm_vsolv%ntheth), ix_i30(qmmm_vsolv%ntheth)
  integer, intent(out) :: ix_i32(qmmm_vsolv%ntheta), ix_i34(qmmm_vsolv%ntheta), ix_i36(qmmm_vsolv%ntheta), ix_i38(qmmm_vsolv%ntheta)
  integer, intent(out) :: ix_i40(qmmm_vsolv%nphih), ix_i42(qmmm_vsolv%nphih), ix_i44(qmmm_vsolv%nphih), &
                          ix_i46(qmmm_vsolv%nphih), ix_i48(qmmm_vsolv%nphih)
  integer, intent(out) :: ix_i50(qmmm_vsolv%nphia), ix_i52(qmmm_vsolv%nphia), ix_i54(qmmm_vsolv%nphia), &
                          ix_i56(qmmm_vsolv%nphia), ix_i58(qmmm_vsolv%nphia)
  integer, intent(out) :: ix_ibellygp(*)

! Local
  integer :: i, oldnumbnd

#ifdef MPI
   include 'mpif.h'
  integer :: ier
#endif

  if ( qmmm_vsolv%debug ) then
     write(6,'(a)') '>>>>> entered qmmm_vsolv_update_variable_solvent'
     call flush(6)
     call print(qmmm_vsolv)
     call flush(6)
  end if

  if (mod(nstep,qmmm_vsolv%nearest_qm_solvent_fq)/=0) return

! Notes about things that might need changing when we change the qm region.
!       Here we are just changing which solvent residues are in the QM region so a 
!       few things to note:
!
!        1) Nquant and nquant_nlink remain the same.
!        2) Solvent molecules are identical so atomic numbers etc do not change.
!
! Changes etc.
! 
! Could potentially need updating for variable QM but not needed for variable solvent
! -----------------------------------------------------------------------------------
! qmmm_struct%qm_resp_charges - contains resp charges for QM atoms as was read from the
!                               prmtop file. Since solvent residues are identical (or assumed to be)
!                               so this does NOT change if we change which solvent residues are in
!                               the QM region.
! qmmm_struct%qm_resp_charge_sum - does not change for same reasons as above.
! qmmm_struct%qm_coords - Does not need to be modified since the number of QM atoms and link atoms
!                         remains constant - this is refilled on every MD step inside QMMM.
! qmmm_struct%scaled_mm_charges - does not change for same reason as qm_resp_charges.
! qmmm_struct%qm_xcrd - Does not need to change since it is rebuilt on every qmmm call and resized as needed.
! qmmm_struct%iqm_atomic_numbers - Does not need to change since solvents are assumed to be identical and
!                                  always come after fixed_nquants - so order will not change.
! 
! Must be updated for variable solvent
! ------------------------------------
! qmmm_struct%iqmatoms - integer list of atom numbers making up QM region - fixed QM, then variable QM, then link.
! qmmm_struct%atom_mask - natom long, set to true if atom is QM - must be updated for new solvent. 
! 
! Bonds, angles and dihedrals - need to be updated.
! 


  ! Serial at the moment so only master should call and then relevant data is sent
  ! to other threads afterwards.
  if (qmmm_mpi%commqmmm_master) then
     ! note: if qmmm_vsolv%recalculate is .true. then we are now in a second call to
     !       force() and are reusing solvent pointers from the previous MD step
     !       THIS IS FOR ADAPTIVE QM/MM, see qmmm_adaptive_module for details
     if ( .not. qmmm_vsolv%recalculate ) then
        ! now call the routine that returns 
        ! nearest_solvent_pointers of size qmmm_vsolv%nearest_qm_solvent+1
        ! the last value points to the first MM solvent 
        call qmmm_vsolv_identify_nearest_solvent(natom, nres, natom_res, unimaged_coords, atomic_masses, &
                     qmmm_vsolv, qmmm_scratch%qm_real_scratch)
     end if
  end if

#ifdef MPI
  call mpi_bcast(qmmm_vsolv%nearest_solvent_pointers, &
                 qmmm_vsolv%nearest_qm_solvent+1,MPI_INTEGER, &
                 0,qmmm_mpi%commqmmm,ier)
  call mpi_bcast(qmmm_vsolv%nearest_solvent_distances, &
                 qmmm_vsolv%nearest_qm_solvent+1,MPI_DOUBLE_PRECISION, &
                 0,qmmm_mpi%commqmmm,ier)

#endif

  ! We should refill the main charge arrays before identifying new solvents for QM
  ! Don't need to restore the MMLink pair atoms since only the solvent is going to 
  ! change here and these can't involve link atoms when running with nearest qm solvent.
  ! in parallel all threads call this for now.
  call qmmm_restore_mm_charges(qmmm_struct%nquant,qmmm_struct%qm_resp_charges,charge, &
                               qmmm_struct%scaled_mm_charges, qmmm_struct%iqmatoms, &
                               qmmm_nml%chg_lambda,qmmm_struct%nlink, &
                               qmmm_struct%link_pairs, &
                               qmmm_struct%mm_link_pair_resp_charges, &
                               .false.)

  ! Update iqmatoms
  ! in parallel all threads call this for now.
  call qmmm_vsolv_update_iqmatoms(qmmm_vsolv%fixed_nquant, &
                                  qmmm_vsolv%nearest_qm_solvent, &
                                  qmmm_vsolv%nearest_solvent_pointers, &
                                  qmmm_vsolv%natom_solv_res, qmmm_struct%nquant, &
                                  qmmm_struct%iqmatoms, natom, &
                                  qmmm_nml%verbosity, qmmm_vsolv%debug)

  ! Also setup the qmmm atom mask array
  qmmm_struct%atom_mask(1:natom) = .false.
  do i = 1, qmmm_struct%nquant   !1 to fixed + variable QM.
    qmmm_struct%atom_mask(qmmm_struct%iqmatoms(i)) = .true.
  end do

  ! Now we have updated the QM atoms we need to zero the charge arrays again.
  ! We don't need to save the charges since they won't have changed since the
  ! very first time we did it.
  ! in parallel all threads call this for now.
  call qm_zero_charges(qmmm_struct, charge,qmmm_struct%scaled_mm_charges,.false.)

  ! Finally we must reinitialize the bonds, angles and dihedrals list
  ! for the new qm atoms.

  ! Note before calling setbon here we will reset the number of bond types to the prmtop
  ! value otherwise you can get problems with, for example runs with qmshake, where the number of
  ! bond types keeps increasing on each call.
  ! We also need to restore the original bond, angle and dihedral types
  oldnumbnd = numbnd
  numbnd = qmmm_vsolv%prmtop_numbnd
  nbonh  = qmmm_vsolv%nbonh
  nbona  = qmmm_vsolv%nbona
  ntheth = qmmm_vsolv%ntheth
  ntheta = qmmm_vsolv%ntheta
  nphih  = qmmm_vsolv%nphih
  nphia  = qmmm_vsolv%nphia

  ix_iibh  = qmmm_vsolv%iibh
  ix_ijbh  = qmmm_vsolv%ijbh
  ix_iicbh = qmmm_vsolv%icbh

  ix_iiba  = qmmm_vsolv%iiba
  ix_ijba  = qmmm_vsolv%ijba
  ix_iicba = qmmm_vsolv%icba

  ix_i24   = qmmm_vsolv%iith
  ix_i26   = qmmm_vsolv%ijth
  ix_i28   = qmmm_vsolv%ikth
  ix_i30   = qmmm_vsolv%icth

  ix_i32   = qmmm_vsolv%iita
  ix_i34   = qmmm_vsolv%ijta
  ix_i36   = qmmm_vsolv%ikta
  ix_i38   = qmmm_vsolv%icta

  ix_i40   = qmmm_vsolv%iiph
  ix_i42   = qmmm_vsolv%ijph
  ix_i44   = qmmm_vsolv%ikph
  ix_i46   = qmmm_vsolv%ilph
  ix_i48   = qmmm_vsolv%icph

  ix_i50   = qmmm_vsolv%iipa
  ix_i52   = qmmm_vsolv%ijpa
  ix_i54   = qmmm_vsolv%ikpa
  ix_i56   = qmmm_vsolv%ilpa
  ix_i58   = qmmm_vsolv%icpa

  ! Bonds with Hydrogen
  if(nbonh > 0) call setbon(nbonh,ix_iibh,ix_ijbh,ix_iicbh,ix_ibellygp)

  ! Bonds without Hydrogen
  if(nbona > 0) call setbon(nbona,ix_iiba,ix_ijba,ix_iicba, ix_ibellygp)

  ! Test that the new numbnd value matches our old value - otherwise something
  ! has changed in the bond topology since we started this run - i.e. 2 residues
  ! with the same name (and in the variable solvent residue) have different bond
  ! types.
  REQUIRE(numbnd==oldnumbnd)

  !Angles with Hydrogen   
  if(ntheth > 0) call setang(ntheth,ix_i24,ix_i26,ix_i28,ix_i30,ix_ibellygp) 
  
  !Angles without Hydrogen 
  if(ntheta > 0) call setang(ntheta,ix_i32,ix_i34,ix_i36,ix_i38,ix_ibellygp) 
   
  !Dihedrals with Hydrogen 
  if(nphih > 0) call setdih(nphih,ix_i40,ix_i42,ix_i44,ix_i46,ix_i48,ix_ibellygp)
   
  !Dihedrals without Hydrogen 
  if(nphia > 0) call setdih(nphia,ix_i50,ix_i52,ix_i54,ix_i56,ix_i58,ix_ibellygp)

  if ( qmmm_vsolv%debug ) then
     write(6,'(a)') '<<<<< leaving qmmm_vsolv_update'
     call flush(6)
  end if

end subroutine qmmm_vsolv_update
