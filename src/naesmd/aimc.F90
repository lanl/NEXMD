#include "dprec.fh"
#include "assert.fh"

module AIMC
    use communism
    use AIMC_type_module, only : AIMC_type
    implicit none

    private

    public :: AIMC_clone_check, AIMC_clone

contains
    logical function AIMC_clone_check(sim)
        implicit none
        !This subroutined determines if a "cloning" of the simulation should take place
        !in accorfance with the AIMC18 cloning criterions
        !PCCP AIMC18 paper Section 2.6
        type(simulation_t), pointer :: sim
        integer :: i, j, k, l, n
        _REAL_ :: Wn, thetan, ratio, FM2, Fmax2, FE2

        AIMC_clone_check = .false.

        if (sim%aimc%nclones .ge. sim%aimc%nclones_max) then
            return
        end if
        Wn = 0.0
        do i = 1, sim%excN
            Wn = Wn + sim%naesmd%yg(i)**4
        end do
        Wn = 1.0/Wn

        ratio = 0.0
        FM2 = 0.0
        Fmax2 = 0.0
        do i = 1, sim%Na*3
            FM2 = FM2 + sim%aimc%FM(i)**2
            Fmax2 = Fmax2 + sim%aimc%Fmax(i)**2
            ratio = ratio + sim%aimc%FM(i)*sim%aimc%Fmax(i)
        end do
        ratio = ratio*2.0/(FM2 + Fmax2)
        thetan = acos(ratio)

        FE2 = 0.0
        do i = 1, sim%Na*3
            FE2 = FE2 + sim%aimc%FE(i)**2
        end do

        !old criterion 3
!        ratio=Sqrt(FE2/FM2)

        !new criterion 3
        ratio = 0.0
        do i = 1, sim%excN
            ratio = ratio + abs(2.0*sim%naesmd%yg(i)/sim%naesmd%yg(sim%aimc%imax)* &
                dcos(sim%naesmd%yg(i + sim%excN) - sim%naesmd%yg(sim%aimc%imax + sim%excN))*sim%naesmd%cadiab(sim%aimc%imax, i))
        end do

        if (Wn .lt. sim%aimc%delta_clone_1) then
            return
        end if

        if (thetan .lt. sim%aimc%delta_clone_2) then
            return
        end if

        if (ratio .gt. sim%aimc%delta_clone_3) then
            return
        end if

        AIMC_clone_check = .true.

889     format(10000(1X, F18.10))

    end function AIMC_clone_check

    subroutine AIMC_clone(old_sim, new_sim, nuclear, Nsim, ii)

!This subroutined modifies the initially identical old and new sim in accordance with the AIMC18 cloning algorithm
!PCCP AIMC18 paper Section 2.5
        implicit none
        type(simulation_t), pointer :: old_sim, new_sim
        type(MCE) :: nuclear
        integer :: i, j, k, l, m, n
        integer :: max_pop
        integer, intent(in) :: ii !Index of the cloning trajectory
        integer, intent(in) :: Nsim
        double complex, allocatable :: D0(:)

        max_pop = maxloc(abs(old_sim%naesmd%yg(1:old_sim%excN)), 1) !Cloned out state

!Redefining nuclear coefficients:
        allocate (D0(Nsim - 1))
        D0 = nuclear%D
        deallocate (nuclear%D)
        allocate (nuclear%D(Nsim))
        do i = 1, Nsim - 1
            nuclear%D(i) = D0(i)
        end do
        nuclear%D(Nsim) = D0(ii)*old_sim%naesmd%yg(max_pop)
        nuclear%D(ii) = D0(ii)*sqrt(1 - old_sim%naesmd%yg(max_pop)**2)

!Logical
        nuclear%cloned(ii) = .true.
        nuclear%cloned(Nsim) = .true.

!Redefining Heff
        nuclear%Heff(:, Nsim) = nuclear%Heff(:, ii)
        nuclear%Heff(Nsim, :) = nuclear%Heff(ii, :)
        nuclear%Heff(Nsim, Nsim) = nuclear%Heff(ii, ii)
        nuclear%Heff(ii, Nsim) = 0.0d0
        nuclear%Heff(Nsim, ii) = 0.0d0
        nuclear%Heff_old(:, Nsim) = nuclear%Heff_old(:, ii)
        nuclear%Heff_old(Nsim, :) = nuclear%Heff_old(ii, :)
        nuclear%Heff_old(Nsim, Nsim) = nuclear%Heff_old(ii, ii)
        nuclear%Heff_old(ii, Nsim) = 0.0d0
        nuclear%Heff_old(Nsim, ii) = 0.0d0

!Redefining electronic overlaps:
        nuclear%sE(Nsim, :, :, :) = nuclear%sE(ii, :, :, :)
        nuclear%sE(:, Nsim, :, :) = nuclear%sE(:, ii, :, :)
        nuclear%sE(Nsim, Nsim, :, :) = nuclear%sE(ii, ii, :, :)

!Redefining electronic coefficients:
        old_sim%naesmd%yg(1:old_sim%excN) = new_sim%naesmd%yg(1:old_sim%excN)/sqrt(1.0 - new_sim%naesmd%yg(max_pop)**2)
        old_sim%naesmd%yg(max_pop) = 0.0

        old_sim%naesmd%yg_new(1:old_sim%excN) = new_sim%naesmd%yg_new(1:old_sim%excN)/sqrt(1.0 - new_sim%naesmd%yg(max_pop)**2)
        old_sim%naesmd%yg_new(max_pop) = 0.0
        old_sim%naesmd%yg_new(old_sim%excN + 1:2*old_sim%excN) = old_sim%naesmd%yg_new(old_sim%excN + 1:2*old_sim%excN)/sqrt(1.0 - new_sim%naesmd%yg(max_pop)**2)
        old_sim%naesmd%yg_new(old_sim%excN + max_pop) = 0.0

        do i = 1, old_sim%excN
            if (i .eq. max_pop) then
                new_sim%naesmd%yg_new(i) = new_sim%naesmd%yg_new(i)/abs(new_sim%naesmd%yg(max_pop))
                new_sim%naesmd%yg_new(i + old_sim%excN) = new_sim%naesmd%yg_new(i + old_sim%excN)/abs(new_sim%naesmd%yg(max_pop))
            else
                new_sim%naesmd%yg_new(i) = 0.0
                new_sim%naesmd%yg_new(i + old_sim%excN) = 0.0
            end if
        end do

        new_sim%naesmd%yg(1:old_sim%excN) = 0.0
        new_sim%naesmd%yg(max_pop) = 1.0 !since new_sim yg is 1 for single state, it won't flag another clone

!print *, old_sim%naesmd%yg
!print *, new_sim%naesmd%yg

    end subroutine AIMC_clone

end module AIMC
