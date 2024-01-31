! <compile=optimized>
!------------------------------------------------------------------------------
! nose_hoover.f
!------------------------------------------------------------------------------
!
! Integrator of Nose'-Hoover chain of thermostat
!
! Ref: Jang & Voth, J. Chem. Phys. 107, 9514 (1997)
!
! Last update: 03/09/2006
!
!------------------------------------------------------------------------------
#include "dprec.fh"

!==============================================================================
module nose_hoover_module
!==============================================================================

    implicit none

    character(*), parameter :: module_name = "nose_hoover_module:"
    logical, save :: module_init = .false.
    integer, parameter :: dp = 8   !! double precision

    integer, parameter :: M = 9   !! max. chain length

    type SystemCoordinate_type
        real(dp) :: mass
        real(dp), pointer :: vel
    end type SystemCoordinate_type

    type Thermostat_type
        real(dp) :: eta(M)   !! position
        real(dp) :: eta_old(M)   !! old position for Leap-frog
        real(dp) :: v(M)   !! velocity
        real(dp) :: v_old(M)   !! old velocity for Leap-frog
        real(dp) :: a(M)   !! acceleration (without friction term)
        real(dp) :: Q(M)   !! mass
        real(dp) :: Q_inv(M)   !! inverse mass
        real(dp) :: kT, Ndof_kT   !! kB * temperature
        type(SystemCoordinate_type), pointer :: system_coords(:)
        integer :: num_system_coords
        integer :: system_coord_id = 0
        logical :: activated   !! flag to activate the thermostat
        logical :: initialized = .false.
    end type Thermostat_type

    !..................................................

    integer, save :: print_level = 0   !! output is suppressed

    !..................................................

!! Random number generator from Numerical Recipes.

    integer, save :: random_number_seed = 100
    integer, save :: idum_ran2

    !..................................................

    type(Thermostat_type), save :: thermo_lnv

    _REAL_, save :: c2_lnv, mass_lnv, v_lnv, f_lnv_v, f_lnv_p, x_lnv, x_lnv_old

contains

!------------------------------------------------------------------------------
    subroutine nose_hoover_module_init
!------------------------------------------------------------------------------

        implicit none

        if (module_init) return

        if (print_level > 0) then
            write (6, *)
            write (6, *) "Initializing ", module_name
        end if

        module_init = .true.

!!
!! Output.
!!

        if (print_level > 0) then
            write (6, *)
            write (6, *) module_name
            write (6, *)
            write (6, *) "   chain length = ", M
        end if

    end subroutine nose_hoover_module_init

!------------------------------------------------------------------------------
    subroutine Thermostat_init &
        (nchain, thermo, num_system_coords, num_constraints, kT, tau, &
        pos_init, vel_init, activate)

! initializes a Thermostat object.
!------------------------------------------------------------------------------

        implicit none

        type(Thermostat_type) :: thermo
        integer, intent(in) :: nchain  !! number of oscillators in each chain
        integer, intent(in) :: num_system_coords   !! number of system coordinates
        !! to be thermostatted
        integer, intent(in) :: num_constraints   !! number of geometrical constraints
        !! in the system
        real(dp), intent(in) :: kT   !! external temperature
        real(dp), intent(in) :: tau   !! characteristic timescale of the system
        real(dp), intent(in), optional :: pos_init, vel_init   !! initial condition
        !! of the thermostat

        logical, intent(in), optional :: activate   !! flag to activate the thermostat

        integer :: Ndof, j

        if (.not. module_init) call nose_hoover_module_init

        if (nchain > M - 1) then
            write (6, *) "Stop. Number of oscillators too large. nchain must be <= ", M - 1
            stop
        end if

        thermo%initialized = .true.

        allocate (thermo%system_coords(num_system_coords))
        thermo%num_system_coords = num_system_coords

        Ndof = num_system_coords - num_constraints

        thermo%Q(1) = kT*tau**2*dble(Ndof)
        thermo%Q(2:nchain) = kT*tau**2
        thermo%Q_inv(1:nchain) = 1.0d0/thermo%Q(1:nchain)

        thermo%Q(nchain + 1:M) = 0.d0
        thermo%Q_inv(nchain + 1:M) = 0.d0

        thermo%kT = kT
        thermo%Ndof_kT = dble(Ndof)*kT

!!
!! Initial position, velocity, and acceleration.
!!

        thermo%eta(:) = 0.0d0
        thermo%eta_old(:) = 0.0d0

        thermo%v(:) = 0.0d0
        thermo%v_old(:) = 0.0d0

        thermo%a(:) = 0.0d0

        do j = 1, nchain
            call gauss(0.d0, sqrt(kT/thermo%Q(j)), thermo%v(j))
            thermo%v_old(j) = thermo%v(j)
        end do

        if (present(pos_init)) thermo%eta(:) = pos_init
        if (present(vel_init)) thermo%v(:) = vel_init

        do j = 2, nchain
            thermo%a(j) = thermo%Q_inv(j) &
                *(thermo%Q(j - 1)*thermo%v(j - 1)**2 - thermo%kT)
            !! Note: thermo%a( 1 ) is computed in (subr.)Thermostat_integrate_1.
        end do

!! Flag to activate the thermostat.

        if (present(activate)) then
            thermo%activated = activate
        else
            thermo%activated = .true.
        end if

    end subroutine Thermostat_init

!------------------------------------------------------------------------------
    subroutine Thermostat_link(thermo, system_mass, system_velocity)

! establishes the link to system velocicities via pointer assignment.
!------------------------------------------------------------------------------

        implicit none

        type(Thermostat_type) :: thermo
        real(dp), intent(in) :: system_mass
        real(dp), intent(in), target :: system_velocity

        integer :: id

        if (.not. thermo%initialized) then
            write (6, *) module_name, " thermostat is not initialized. stop."
            stop
        end if

        thermo%system_coord_id = thermo%system_coord_id + 1
        id = thermo%system_coord_id
        if (id <= thermo%num_system_coords) then
            thermo%system_coords(id)%mass = system_mass
            thermo%system_coords(id)%vel => system_velocity   !! establish link
        else
            write (6, *) module_name, " thermo%num_system_coords is too small. stop."
            stop
        end if

    end subroutine Thermostat_link

!------------------------------------------------------------------------------
    subroutine Thermostat_switch(thermo, activate)
!
! Switch of thermostats
!------------------------------------------------------------------------------

        implicit none

        type(Thermostat_type) :: thermo
        logical, intent(in) :: activate

        thermo%activated = activate
        thermo%v(:) = 0.d0
        thermo%v_old(:) = 0.d0

    end subroutine Thermostat_switch

!------------------------------------------------------------------------------
    subroutine Thermostat_integrate_1 &
        (nchain, thermostats, num_thermostats, timestep, ntp)
!------------------------------------------------------------------------------

        implicit none

        integer, intent(in) :: nchain
        integer, intent(in) :: num_thermostats
        type(Thermostat_type), target :: thermostats(num_thermostats)
        real(dp), intent(in) :: timestep

        integer :: ithermo, j, k, ntp, istart3, iend3
        real(dp) :: dt, hdt, kT, Ndof_kT, Ekin2, exp1, exp2
        type(Thermostat_type), pointer :: thermo
        type(SystemCoordinate_type), pointer :: system(:)
#ifdef MPI
# include "parallel.h"
#endif

        dt = timestep
        hdt = 0.5d0*dt

        !------------------------------
#ifdef MPI
        istart3 = iparpt3(mytaskid) + 1
        iend3 = iparpt3(mytaskid + 1)
#else
        istart3 = 1
        iend3 = num_thermostats
#endif

        do ithermo = 1, num_thermostats
            thermo => thermostats(ithermo)
            if (.not. thermo%activated) cycle
            !! Error check.
            if (.not. thermo%initialized) goto 9000
            if (thermo%num_system_coords /= thermo%system_coord_id) goto 9100
            system => thermo%system_coords(:)
            kT = thermo%kT
            Ndof_kT = thermo%Ndof_kT
            !------------------------------
            !! Update coordinates for oscillators of Nose'-Hoover chains.
            !! - odd positions and forces.
            !! - even velocities.
            do j = 1, nchain, 2
                thermo%eta_old(j) = thermo%eta(j)
                thermo%eta(j) = thermo%eta(j) + dt*thermo%v(j)
            end do
            do j = 2, nchain, 2
                thermo%v_old(j) = thermo%v(j)
                exp1 = 1.d0
                if (j < nchain) exp1 = exp(-hdt*thermo%v(j + 1))
                exp2 = exp1*exp1
                thermo%v(j) = thermo%v(j)*exp2 &
                    + dt*thermo%a(j)*exp1
            end do
            do j = 3, nchain, 2
                thermo%a(j) = thermo%Q_inv(j) &
                    *(thermo%Q(j - 1)*thermo%v(j - 1)**2 - kT)
            end do
            !! Update force acting on the first thermostat oscillator.
            Ekin2 = 0.0d0

            do k = 1, thermo%num_system_coords
                Ekin2 = Ekin2 + system(k)%mass*system(k)%vel**2
            end do
            thermo%a(1) = thermo%Q_inv(1)*(Ekin2 - Ndof_kT)
        end do  !! loop over thermostats

        if (ntp > 0) then
            !------------------------------
            kT = thermo_lnv%kT
            do j = 1, nchain, 2
                thermo_lnv%eta_old(j) = thermo_lnv%eta(j)
                thermo_lnv%eta(j) = thermo_lnv%eta(j) + dt*thermo_lnv%v(j)
            end do

            do j = 2, nchain, 2
                thermo_lnv%v_old(j) = thermo_lnv%v(j)
                exp1 = 1.d0
                if (j < nchain) exp1 = exp(-hdt*thermo_lnv%v(j + 1))
                exp2 = exp1*exp1
                thermo_lnv%v(j) = thermo_lnv%v(j)*exp2 &
                    + dt*thermo_lnv%a(j)*exp1
            end do

            do j = 3, nchain, 2
                thermo_lnv%a(j) = thermo_lnv%Q_inv(j) &
                    *(thermo_lnv%Q(j - 1)*thermo_lnv%v(j - 1)**2 - kT)
            end do
            !! Update force acting on the first thermostat oscillator.
            Ekin2 = 0.0d0
            thermo_lnv%a(1) = thermo_lnv%Q_inv(1)*(mass_lnv*v_lnv*v_lnv - kT)
        end if

        return

!!
!! Error handling.
!!

9000    continue
        write (6, *) module_name, " thermostat is not initialized. stop."
        stop

9100    continue
        write (6, *) module_name, " link to system velocities is incomplete. stop."
        stop

    end subroutine Thermostat_integrate_1

!------------------------------------------------------------------------------
    subroutine Thermostat_integrate_2 &
        (nchain, thermostats, num_thermostats, timestep, ntp)
!------------------------------------------------------------------------------

        implicit none

        integer, intent(in) :: nchain, ntp
        integer, intent(in) :: num_thermostats
        type(Thermostat_type), target :: thermostats(num_thermostats)
        real(dp), intent(in) :: timestep

        integer :: ithermo, j, istart3, iend3
        real(dp) :: dt, hdt, kT, Ndof_kT, exp1, exp2
        type(Thermostat_type), pointer :: thermo
        type(SystemCoordinate_type), pointer :: system(:)

#ifdef MPI
# include "parallel.h"
#endif

        dt = timestep
        hdt = 0.5d0*dt

        !------------------------------
#ifdef MPI
        istart3 = iparpt3(mytaskid) + 1
        iend3 = iparpt3(mytaskid + 1)
#else
        istart3 = 1
        iend3 = num_thermostats
#endif

        do ithermo = 1, num_thermostats
            thermo => thermostats(ithermo)
            if (.not. thermo%activated) cycle
            system => thermo%system_coords(:)
            kT = thermo%kT
            Ndof_kT = thermo%Ndof_kT

            !------------------------------

            !! Update coordinates for oscillators of Nose'-Hoover chains.
            !! - even positions and forces.
            !! - odd velocities.

            do j = 2, nchain, 2
                thermo%eta(j) = thermo%eta(j) + dt*thermo%v(j)
                thermo%eta_old(j) = thermo%eta(j)
            end do

            do j = 1, nchain, 2
                exp1 = 1.d0
                if (j < nchain) exp1 = exp(-hdt*thermo%v(j + 1))
                exp2 = exp1**2
                thermo%v(j) = thermo%v(j)*exp2 &
                    + dt*thermo%a(j)*exp1
                thermo%v_old(j) = thermo%v(j)
            end do

            do j = 2, nchain, 2
                thermo%a(j) = thermo%Q_inv(j) &
                    *(thermo%Q(j - 1)*thermo%v(j - 1)**2 - kT)
            end do

        end do  !! loop over thermostats

        if (ntp > 0) then
            kT = thermo%kT
            do j = 2, nchain, 2
                thermo_lnv%eta(j) = thermo_lnv%eta(j) + dt*thermo_lnv%v(j)
                thermo_lnv%eta_old(j) = thermo_lnv%eta(j)
            end do

            do j = 1, nchain, 2
                exp1 = 1.d0
                if (j < nchain) exp1 = exp(-hdt*thermo_lnv%v(j + 1))
                exp2 = exp1**2

                thermo_lnv%v(j) = thermo_lnv%v(j)*exp2 &
                    + dt*thermo_lnv%a(j)*exp1

                thermo_lnv%v_old(j) = thermo_lnv%v(j)
            end do

            do j = 2, nchain, 2
                thermo_lnv%a(j) = thermo_lnv%Q_inv(j) &
                    *(thermo_lnv%Q(j - 1)*thermo_lnv%v(j - 1)**2 - kT)
            end do

        end if
        !------------------------------
    end subroutine Thermostat_integrate_2

!------------------------------------------------------------------------------
    function Thermostat_hamiltonian(nchain, thermostats, num_thermostats) result(E)

! computes the thermostat terms in the extended Hamiltonian.
!------------------------------------------------------------------------------

        implicit none

        integer, intent(in) :: nchain
        integer, intent(in) :: num_thermostats
        type(Thermostat_type), intent(in), target :: thermostats(num_thermostats)
        real(dp) :: E, av_v, av_eta
        integer :: j, ithermo
        type(Thermostat_type), pointer :: thermo

#ifdef MPI
        include 'mpif.h'
#  include "parallel.h"
        integer :: ierr
        _REAL_ :: Etmp
#endif
        integer :: istart3, iend3
#include "md.h"
#include "memory.h"

        E = 0.0d0

#ifdef MPI
        if (mpi_orig) then
            istart3 = 1
            iend3 = natom*3
        else
            istart3 = 3*iparpt(mytaskid) + 1
            iend3 = 3*iparpt(mytaskid + 1)
        end if
#else
        istart3 = 1
        iend3 = 3*natom
#endif

        do ithermo = istart3, iend3
            thermo => thermostats(ithermo)
            av_v = 0.5d0*(thermo%v_old(1) + thermo%v(1))
            av_eta = 0.5d0*(thermo%eta_old(1) + thermo%eta(1))
            E = E + 0.5d0*thermo%Q(1)*av_v**2 + thermo%Ndof_kT*av_eta
            do j = 2, nchain
                av_v = 0.5d0*(thermo%v_old(j) + thermo%v(j))
                av_eta = 0.5d0*(thermo%eta_old(j) + thermo%eta(j))
                E = E + 0.5d0*thermo%Q(j)*av_v**2 + thermo%kT*av_eta
            end do
        end do

#ifdef MPI
        call mpi_allreduce(E, Etmp, 1, MPI_DOUBLE_PRECISION, mpi_sum, commsander, ierr)
        E = Etmp
#endif

        if (ntp > 0) then
            do j = 1, nchain
                av_v = 0.5d0*(thermo_lnv%v_old(j) + thermo_lnv%v(j))
                av_eta = 0.5d0*(thermo_lnv%eta_old(j) + thermo_lnv%eta_old(j))
                E = E + 0.5*thermo_lnv%Q(j)*av_v*av_v + thermo_lnv%kT*av_eta
            end do
        end if

    end function Thermostat_hamiltonian

end module nose_hoover_module
