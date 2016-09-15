#include "dprec.fh"
#include "assert.fh"

module md_module
    implicit none
    !mode and atom variables:
    integer,allocatable:: atoms(:)        ! atomic numbers (currently max 1000)
    _REAL_,allocatable:: atmass(:)        ! atomic masses (currently max 1000)
    _REAL_,allocatable:: v(:,:)         ! unit matrix
    _REAL_,allocatable:: fm(:)            ! mode masses
    _REAL_,allocatable:: fo(:)            ! mode frequencies (used to determine friction)
    _REAL_,allocatable:: r0(:)            ! cartesian coordinates of the origin
    _REAL_  ttt                ! Initial Temperature, K
    integer imdtype           ! MD type: 0-ground,1-first exc. state
    integer ideriv            ! MD derivatives: 0-numerical, 1- analytic
    integer icart            ! 0 - along Cartesians, 1 - along vibrations
    integer ifric          ! Friction: 0-No, 1-yes

contains
    subroutine allocate_md_module_init(natoms)
        implicit none
        integer natoms
        write(6,*)'Allocating md_module initial variables'
        allocate(atoms(natoms),atmass(natoms),r0(3*natoms))
    end
    subroutine allocate_md_module(Na)
        implicit none
        integer Na,Nm
        write(6,*)'Allocating md_module variables',Na
        Nm=3*Na
        if(.not.allocated(v)) allocate(v(Nm,Nm))
        if(.not.allocated(fm)) allocate(fm(Nm)) ! mode masses
        if(.not.allocated(fo)) allocate(fo(Nm))            ! mode frequencies (used to determine friction)
    end
end module
