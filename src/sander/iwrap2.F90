! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

! 2008-10-03
! Gustavo Seabra, University of Florida
!
! This routine wraps the system around the center of
! mass of a given set of atoms.

subroutine iwrap2(nwrap_atoms, wrap_atoms, coords, box_center)

    use molecule, only : mol_info
    implicit none

    ! Passed in
    integer, intent(in)    :: nwrap_atoms             ! Number of atoms to be kept in the center of box
    integer, intent(in)    :: wrap_atoms(nwrap_atoms) ! Atom numbers of the atoms to be kept in the center
    _REAL_, intent(inout) :: coords(3*mol_info%natom)! in:  The unimaged coordinates;
    _REAL_, intent(out) :: box_center(3)
    ! out: Imaged and wrapped coordinates
    ! Local
    _REAL_  :: iwrap_mask_com(3), shift(3)
    _REAL_  :: wrap_mass_inv, this_atom_mass
    integer :: m, atom, this_atom, this_atom_x, wrap_debug

! Included
#include "box.h"  /* Provides: ifbox  : box type */
    /*box(3):box dimensions*/

!--

    wrap_debug = 0
    !debug: write a restart file
    if (wrap_debug > 1) then
        write (35, '(a)') "ORIGINAL COORDINATES"
        write (35, '(i9)') mol_info%natom
        write (35, '(6f12.7)') (coords(m), m=1, mol_info%natom*3)
    end if
    ! 1. Find the COM of the atoms in the iwrap_mask
    iwrap_mask_com(1:3) = 0.d0
    wrap_mass_inv = 0.d0

    if (wrap_debug > 0) then
        write (6, '(a)') "IWRAP2 -- DETERMINING THE CENTER OF MASS OF IWRAP_MASK:"
        write (6, '(a)') "IWRAP2 -- "
        write (6, '(a10,a8,4a15)') "IWRAP2 -- ", "ATOM", "X", "Y", "Z", "MASS"
    end if
    do atom = 1, nwrap_atoms
        this_atom = wrap_atoms(atom)
        this_atom_mass = mol_info%atom_mass(this_atom)
        wrap_mass_inv = wrap_mass_inv + this_atom_mass
        this_atom_x = 1 + 3*(this_atom - 1)
        if (wrap_debug > 0) then
            write (6, '(a10,i9,4f15.6)') &
                "IWRAP2 -- ", this_atom, (coords(this_atom_x + m), m=0, 2), this_atom_mass
        end if
        iwrap_mask_com(1) = iwrap_mask_com(1) + coords(this_atom_x)*this_atom_mass
        iwrap_mask_com(2) = iwrap_mask_com(2) + coords(this_atom_x + 1)*this_atom_mass
        iwrap_mask_com(3) = iwrap_mask_com(3) + coords(this_atom_x + 2)*this_atom_mass
    end do

    wrap_mass_inv = 1.d0/wrap_mass_inv

    iwrap_mask_com(1) = iwrap_mask_com(1)*wrap_mass_inv
    iwrap_mask_com(2) = iwrap_mask_com(2)*wrap_mass_inv
    iwrap_mask_com(3) = iwrap_mask_com(3)*wrap_mass_inv

    if (wrap_debug > 0) then
        write (6, '(a10,a8,3f15.6)') &
            "IWRAP2 -- ", "COM", (iwrap_mask_com(m), m=1, 3)
    end if
    ! 2. Translate the system to put the COM of the
    !    iwrap_mask at the center of the box

    box_center(1:3) = box(1:3)*0.5d0
    shift(1:3) = iwrap_mask_com(1:3) - box_center(1:3)
    if (wrap_debug > 0) then
        iwrap_mask_com(1:3) = iwrap_mask_com(1:3) - shift(1:3)
        write (6, '(a,3f15.6)') "IWRAP2 -- SHIFT: ", (shift(m), m=1, 3)
        write (6, '(a)') "IWRAP2 -- BOX DIMENSIONS:"
        write (6, '(a,3f15.6)') "IWRAP2 -- ", (box(m), m=1, 3)
        write (6, '(a)') "IWRAP2 -- FINAL COORDINATES OF WRAPMASK ATOMS' COM:"
        write (6, '(a10,a6,2a10)') "IWRAP2 -- ", "CRD", "COM", "BOXCNTR"
        write (6, '(a,2f10.5)') "IWRAP2 --    X : ", iwrap_mask_com(1), box_center(1)
        write (6, '(a,2f10.5)') "IWRAP2 --    Y : ", iwrap_mask_com(2), box_center(2)
        write (6, '(a,2f10.5)') "IWRAP2 --    Z : ", iwrap_mask_com(3), box_center(3)
    end if

    ! Now, apply the shift to all atoms
    do atom = 1, mol_info%natom
        this_atom_x = 1 + 3*(atom - 1)
        coords(this_atom_x) = coords(this_atom_x) - shift(1)
        coords(this_atom_x + 1) = coords(this_atom_x + 1) - shift(2)
        coords(this_atom_x + 2) = coords(this_atom_x + 2) - shift(3)
    end do

    if (wrap_debug > 0) then
        write (6, '(a)') "IWRAP2 --  TRANSLATED COORDINATED FOR THE IWRAP_MASK:"
        write (6, '(a)') "IWRAP2 -- "
        write (6, '(a10,a8,4a15)') "IWRAP2 -- ", "ATOM", "X", "Y", "Z", "MASS"
        do atom = 1, nwrap_atoms
            this_atom = wrap_atoms(atom)
            this_atom_mass = mol_info%atom_mass(this_atom)
            this_atom_x = 1 + 3*(this_atom - 1)
            write (6, '(a10,i9,4f15.6)') &
                "IWRAP2 -- ", this_atom, (coords(this_atom_x + m), m=0, 2), this_atom_mass

        end do
    end if

    if (wrap_debug > 1) then
        write (36, '(a)') "TRANSLATED UNWRAPPED COORDINATES"
        write (36, '(i9)') mol_info%natom
        write (36, '(6f12.7)') (coords(m), m=1, mol_info%natom*3)
    end if

    ! 3. Call wrap molecules to wrap the whole system
    !    The COM of the iwrap_mask will be in the center
    !    of the box.
    call wrap_molecules(mol_info%nres, mol_info%natom_res, coords)
    if (ifbox == 2) call wrap_to(mol_info%nres, mol_info%natom_res, coords, box)

    !debug: write a restart file
    if (wrap_debug > 1) then
        write (37, '(a)') "FINAL WRAPPED COORDINATES"
        write (37, '(i9)') mol_info%natom
        write (37, '(6f12.7)') (coords(m), m=1, mol_info%natom*3)
        stop
    end if

    return
end subroutine iwrap2
