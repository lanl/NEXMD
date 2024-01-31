#include "dprec.fh"
#include "assert.fh"

!**********************************************************************!
!*This module contains auxiliary subroutines for dropping trajectories*!
!that reach S0/S1 conincal intersections                              *!
!**********************************************************************!
module dropout_module
    use naesmd_constants
    use communism
    use AIMC
    use clone_module
    use AIMC_type_module

    implicit none

contains

    subroutine dropout(sims, nuclear, S0_S1_threshold_tr, Nsim, drops, Nsim_max)
        use naesmd_constants
        use communism
        use AIMC
        use clone_module
        use AIMC_type_module

        implicit none
        type(sim_pointer_array), dimension(Nsim), intent(in) :: sims
        type(mce) :: nuclear
        integer, intent(in) :: S0_S1_threshold_tr !index of the dropout trajectory
        integer, intent(inout) :: Nsim!Number of trajectories
        integer, intent(inout), dimension(:) :: drops(Nsim_max)
        integer, intent(in) :: Nsim_max

        integer i, j

        if (Nsim .eq. 1) then
            print *, 'S0/S1 conical intersection reached'
            print *, 'S0/S1 energy gap: ', (sims(S0_S1_threshold_tr)%sim%naesmd%vgs - sims(S0_S1_threshold_tr)%sim%naesmd%vmdqt)*feVmdqt, ' eV'
            print *, 'less than: ', sims(S0_S1_threshold_tr)%sim%naesmd%S0_S1_threshold, 'eV threshold'
            stop
        else
            print *, 'S0/S1 conical intersection reached for trajectory ', S0_S1_threshold_tr
            print *, 'S0/S1 energy gap: ', (sims(S0_S1_threshold_tr)%sim%naesmd%vgs - sims(S0_S1_threshold_tr)%sim%naesmd%vmdqt)*feVmdqt, ' eV'
            print *, 'less than: ', sims(S0_S1_threshold_tr)%sim%naesmd%S0_S1_threshold, 'eV threshold'

            write (nuclear%outfile_6, *) sims(S0_S1_threshold_tr)%sim%naesmd%tfemto, sims(S0_S1_threshold_tr)%sim%id + sum(drops)
            drops(S0_S1_threshold_tr) = 1
!Renormalization
            do i = 1, Nsim
                if (i .ne. S0_S1_threshold_tr) then
                    nuclear%D(i) = nuclear%D(i)/sqrt(1 - abs(nuclear%D(S0_S1_threshold_tr))**2)
                end if
            end do
!Index rearrengement
            do i = S0_S1_threshold_tr, Nsim - 1
                sims(i)%sim%forces_allocated = .true.
                sims(i + 1)%sim%forces_allocated = .true.
                call clone_sim(sims(i + 1)%sim, sims(i)%sim)
                nuclear%D(i) = nuclear%D(i + 1)
                nuclear%cloned(i) = nuclear%cloned(i + 1)
            end do
            do i = S0_S1_threshold_tr, Nsim - 1
                nuclear%sE(i, :, :, :) = nuclear%sE(i + 1, :, :, :)
                nuclear%Heff(i, :) = nuclear%Heff(i + 1, :)
                nuclear%Heff_old(i, :) = nuclear%Heff_old(i + 1, :)
            end do
            do i = S0_S1_threshold_tr, Nsim - 1
                nuclear%sE(:, i, :, :) = nuclear%sE(:, i + 1, :, :)
                nuclear%Heff(:, i) = nuclear%Heff(:, i + 1)
                nuclear%Heff_old(:, i) = nuclear%Heff_old(:, i + 1)
            end do
            Nsim = Nsim - 1
        end if

    end subroutine dropout

end module dropout_module
