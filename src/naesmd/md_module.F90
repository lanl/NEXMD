module md_module
    implicit none
    !mode and atom variables:
    !real(8),allocatable::xx(:),yy(:),zz(:)!,dij(:)
    !real(8),allocatable::xxp(:),yyp(:),zzp(:)
    !real(8),allocatable::xxm(:),yym(:),zzm(:)
    integer,allocatable:: atoms(:)        ! atomic numbers (currently max 1000)
    real*8,allocatable:: atmass(:)        ! atomic masses (currently max 1000)
    real*8,allocatable:: v(:,:)         ! unit matrix
    real*8,allocatable:: fm(:)            ! mode masses
    real*8,allocatable:: fo(:)            ! mode frequencies (used to determine friction)
    real*8,allocatable:: r0(:)            ! cartesian coordinates of the origin
    real*8  ttt                ! Initial Temperature, K
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
        !if(.not.allocated(dij)) allocate(dij(Na*3))
        !allocate(xx(Na),yy(Na),zz(Na))
        !allocate(xxp(Na),yyp(Na),zzp(Na))
        !allocate(xxm(Na),yym(Na),zzm(Na))
    end
end module
