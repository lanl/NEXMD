module md_module
      implicit none
      ! mode and atom variables:
      integer atoms(1000)        ! atomic numbers (currently max 1000)
      real*8  atmass(1000)        ! atomic masses (currently max 1000)
      real*8,allocatable:: v(:,:)         ! unit matrix
      real*8,allocatable:: fm(:)            ! mode masses
      real*8,allocatable:: fo(:)            ! mode frequencies (used to determine friction)
      real*8,allocatable:: r0(:)            ! cartesian coordinates of the origin
      real*8  ttt                ! Initial Temperature, K
      integer imdtype           ! MD type: 0-ground,1-first exc. state
      integer ideriv            ! MD derivatives: 0-numerical, 1- analytic
      integer icart            ! 0 - along Cartesians, 1 - along vibrations
      integer ifric             ! Friction: 0-No, 1-yes
      contains
      subroutine allocate_md_module(Na)
         implicit none
         integer Na,Nm !Number of atoms, Number of normal modes
         Nm=Na 
         allocate(v(3*Na,Na),fm(Nm)) ! mode masses
         allocate(fo(Nm))            ! mode frequencies (used to determine friction)
         allocate(r0(3*Na))          ! atomic numbers
      end
end module
