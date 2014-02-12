! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"
module molecule
  implicit none
  
  ! If iwrap==2
  integer                          :: n_iwrap_mask_atoms ! number of atoms in iwrap_mask
  integer,allocatable,dimension(:) :: iwrap_mask_atoms  ! list of atoms in iwrap_mask

  type molecule_static_information
    ! Contains general molecule information that doesn't
    ! change during the course of the simulation
    integer :: natom ! Total number of atoms in the system
    integer :: nres  ! Total number of residues in the system
    integer, dimension(:), pointer :: natom_res ! (nres) number of atoms per residue
    integer, dimension(:), pointer :: atom_to_resid_map !(natom) residue number for each atom
    _REAL_ , dimension(:), pointer :: atom_mass ! (natom) atomic masses
  end type molecule_static_information
  type ( molecule_static_information ) mol_info

contains

subroutine allocate_molecule()
  implicit none
  integer :: ier

  allocate( mol_info%natom_res(mol_info%nres), STAT=ier )
  REQUIRE( ier == 0 )

  allocate( mol_info%atom_mass(mol_info%natom), STAT=ier )
  REQUIRE( ier == 0 )

  allocate( mol_info%atom_to_resid_map(mol_info%natom), STAT=ier )
  REQUIRE( ier == 0 )

end subroutine allocate_molecule

subroutine deallocate_molecule()
  implicit none
  integer :: ier

  deallocate( mol_info%natom_res, STAT=ier)
  REQUIRE( ier == 0 )
  
  deallocate( mol_info%atom_mass, STAT=ier)
  REQUIRE( ier == 0 )

  deallocate( mol_info%atom_to_resid_map, STAT=ier)
  REQUIRE( ier == 0 )

end subroutine deallocate_molecule
end module molecule
