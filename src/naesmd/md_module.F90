#include "dprec.fh"
#include "assert.fh"

module md_module
    implicit none
    !mode and atom variables:
    private

    !data types
    public :: md_structure

    !subroutines
    public :: allocate_md_module_init
    public :: allocate_md_module

    !DATA Types
    type md_structure
        integer, allocatable:: atoms(:)        ! atomic numbers (currently max 1000)
        _REAL_, allocatable:: atmass(:)        ! atomic masses (currently max 1000)
        _REAL_, allocatable:: v(:, :)         ! unit matrix
        _REAL_, allocatable:: fm(:)            ! mode masses
        _REAL_, allocatable:: fo(:)            ! mode frequencies (used to determine friction)
        _REAL_, allocatable:: r0(:)            ! cartesian coordinates of the origin
        _REAL_ ttt                ! Initial Temperature, K
        integer imdtype           ! MD type: 0-ground,1-first exc. state
        integer ideriv            ! MD derivatives: 0-numerical, 1- analytic
        integer icart            ! 0 - along Cartesians, 1 - along vibrations
        integer ifric          ! Friction: 0-No, 1-yes
    end type md_structure

contains
    subroutine allocate_md_module_init(md_struct, natoms)
        implicit none
        type(md_structure), intent(inout) :: md_struct
        integer natoms

        write (6, *) 'Allocating md_module initial variables'
        allocate (md_struct%atoms(natoms), md_struct%atmass(natoms), md_struct%r0(3*natoms))
    end subroutine allocate_md_module_init
    subroutine allocate_md_module(md_struct, Na)
        implicit none
        type(md_structure), intent(inout) :: md_struct
        integer Na, Nm
        write (6, *) 'Allocating md_module variables', Na
        Nm = 3*Na
        if (.not. allocated(md_struct%v)) allocate (md_struct%v(Nm, Nm))
        if (.not. allocated(md_struct%fm)) allocate (md_struct%fm(Nm)) ! mode masses
        if (.not. allocated(md_struct%fo)) allocate (md_struct%fo(Nm))            ! mode frequencies (used to determine friction)
    end subroutine allocate_md_module
end module md_module
