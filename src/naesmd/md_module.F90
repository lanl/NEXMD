module md_module
      implicit none
      ! mode and atom variables:
      real*8 v(N3_M,Nm_M)         ! unit matrix
      real*8 fm(Nm_M)            ! mode masses
      real*8 atmass(Na_M)        ! atomic masses
      real*8 fo(Nm_M)            ! mode frequencies (used to determine friction)
      real*8 r0(N3_M)            ! cartesian coordinates of the origin
      integer atoms(Na_M)        ! atomic numbers
      real*8 ttt                ! Initial Temperature, K
      integer imdtype           ! MD type: 0-ground,1-first exc. state
      integer ideriv            ! MD derivatives: 0-numerical, 1- analytic
      integer icart            ! 0 - along Cartesians, 1 - along vibrations
      integer ifric             ! Friction: 0-No, 1-yes

      function allocate_md_module(Na)
         implicit none
         integer Na,Nm !Number of atoms, Number of normal modes
         Nm=Na 
         allocate(v(3*Na,Na),fm(Nm)) ! mode masses
         allocate(atmass(Na))        ! atomic masses
         allocate(fo(Nm))            ! mode frequencies (used to determine friction)
         allocate(r0(3*Na))          ! atomic numbers
         allocate(atoms(Na))        ! atomic numbers
      end function
end module
