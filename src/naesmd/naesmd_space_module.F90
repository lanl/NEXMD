#include "dprec.fh"
#include "assert.fh"

module naesmd_space_module

implicit none

    type realp_t
        _REAL_, dimension(:), pointer :: p
    end type realp_t


    type xyz_t
        _REAL_, dimension(:), pointer  :: x
        _REAL_, dimension(:), pointer  :: y
        _REAL_, dimension(:), pointer  :: z

        type(realp_t), dimension(3)    :: v 

        _REAL_, dimension(:), pointer  :: xold
        _REAL_, dimension(:), pointer  :: yold
        _REAL_, dimension(:), pointer  :: zold

        type(realp_t), dimension(3)    :: vold
    end type xyz_t


    type xyz_and_derivs_t
        type(xyz_t),pointer :: r
        type(xyz_t),pointer :: v
        type(xyz_t),pointer :: a

        integer                        :: Na
       ! integer                        :: Nm
    end type xyz_and_derivs_t



    type  naesmd_data_t
        type(xyz_t)     :: r,v,a,deltaRp, deltaRm
        _REAL_, pointer :: mass(:)
        integer,pointer :: atomtype(:)

        _REAL_,pointer  :: Omega(:)
        _REAL_,pointer  :: E0

        integer         :: Na,Nm
        integer         :: npot

        _REAL_,pointer  :: dtnact
        _REAL_,pointer  :: dtmdqt

      ! cumulative times of sqm and davidson executaions from
      ! cpu_time - initially set to zero

    end type naesmd_data_t

    contains

    subroutine v_assoc(xyz)
        type(xyz_t), intent(inout) :: xyz

        xyz%v(1)%p => xyz%x
        xyz%v(2)%p => xyz%y
        xyz%v(3)%p => xyz%z

        xyz%vold(1)%p => xyz%xold
        xyz%vold(2)%p => xyz%yold
        xyz%vold(3)%p => xyz%zold
    end subroutine

    subroutine alloc_xyz(v,sz)
        type(xyz_t), intent(inout) :: v
        integer,     intent(inout) :: sz
        allocate(v%x(sz))
        allocate(v%y(sz))
        allocate(v%z(sz))

        allocate(v%xold(sz))
        allocate(v%yold(sz))
        allocate(v%zold(sz))

        v%x    = 0.d0
        v%x    = 0.d0
        v%x    = 0.d0
        v%xold = 0.d0
        v%yold = 0.d0
        v%zold = 0.d0
    end subroutine


    subroutine init_naesmd_space_globs(nsp,Na,Nm,pOmega, E0)
        type(naesmd_data_t), pointer :: nsp
        integer,intent(in)           :: Na, Nm
        _REAL_,              pointer :: pOmega(:), E0

        include 'sizes'
        include 'common'

        nsp%r%x    => rx
        nsp%r%y    => ry
        nsp%r%z    => rz
        nsp%r%xold => rxold
        nsp%r%yold => ryold
        nsp%r%zold => rzold
        call v_assoc(nsp%r)

        nsp%v%x    => vx
        nsp%v%y    => vy
        nsp%v%z    => vz
        nsp%v%xold => vxold
        nsp%v%yold => vyold
        nsp%v%zold => vzold
        call v_assoc(nsp%v)


        nsp%a%x    => ax
        nsp%a%y    => ay
        nsp%a%z    => az
        nsp%a%xold => axold
        nsp%a%yold => ayold
        nsp%a%zold => azold
        call v_assoc(nsp%a)

        nsp%deltaRp%x    => deltaxxpnew
        nsp%deltaRp%y    => deltayypnew
        nsp%deltaRp%z    => deltazzpnew
        nsp%deltaRp%xold => deltaxxpold
        nsp%deltaRp%yold => deltayypold
        nsp%deltaRp%zold => deltazzpold
        call v_assoc(nsp%deltaRp)

        nsp%deltaRm%x    => deltaxxmnew
        nsp%deltaRm%y    => deltayymnew
        nsp%deltaRm%z    => deltazzmnew
        nsp%deltaRm%xold => deltaxxmold
        nsp%deltaRm%yold => deltayymold
        nsp%deltaRm%zold => deltazzmold
        call v_assoc(nsp%deltaRm)

        nsp%mass     => massmdqt
        nsp%atomtype => atomtype

        nsp%Omega  => pOmega
        nsp%E0     => E0

        nsp%npot   = npot

        nsp%dtnact => dtnact
        nsp%dtmdqt => dtmdqt

    end subroutine


end module

!     subroutine init_xyz_and_derivs(s, r,v,a,Na)
!         type(xyz_and_derivs_t), intent(inout) :: s
!         type(xyz_t), target,    intent(in)    :: r,x,a
!         s%r => r
!         s%v => v
!         s%a => a
!     end subroutine
! 
