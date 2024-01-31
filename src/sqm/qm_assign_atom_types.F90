! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
subroutine qm_assign_atom_types(qmmm_struct)
!----------------------------------
!Written by Ross Walker (TSRI,2005)
!
! This routine will go through the atomic numbers
! for each of the QM atoms and work out how many
! different elements there are. It will then assign
! an atom type to each QM atom. This atom type is
! essentially a rebasing of the atomic number in
! order to save memory.
!----------------------------------

    use ElementOrbitalIndex, only : numberElements
    use qmmm_struct_module, only : qmmm_struct_type

    implicit none

    type(qmmm_struct_type), intent(inout) :: qmmm_struct
    integer ier
    integer i, j, natqmi
    logical assigned

    allocate (qmmm_struct%qm_atom_type(qmmm_struct%nquant_nlink), stat=ier)
    REQUIRE(ier == 0) !Deallocated in deallocate_qmmm

    qmmm_struct%qm_atom_type(1:qmmm_struct%nquant_nlink) = 0
    qmmm_struct%qm_type_id(1:numberElements) = 0
    qmmm_struct%qm_ntypes = 0
    do i = 1, qmmm_struct%nquant_nlink
        natqmi = qmmm_struct%iqm_atomic_numbers(i)
        assigned = .false.
        do j = 1, qmmm_struct%qm_ntypes
            !Loop over the number of types found so far and see if
            !this atomic number has already been assigned to a type.
            if (natqmi == qmmm_struct%qm_type_id(j)) then
                assigned = .true.
                qmmm_struct%qm_atom_type(i) = j
            end if
        end do
        !if assigned is true here then we have already assigned this type
        !we just simply move to the next atom.
        if (.NOT. assigned) then
            qmmm_struct%qm_ntypes = qmmm_struct%qm_ntypes + 1
            qmmm_struct%qm_type_id(qmmm_struct%qm_ntypes) = natqmi
            qmmm_struct%qm_atom_type(i) = qmmm_struct%qm_ntypes
        end if
    end do

end subroutine qm_assign_atom_types
