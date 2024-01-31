! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

!*******************************************************************************
!
! Module: amd_mod
!
! Description:
!
! Module for controlling accelerated molecular dynamics calculations
!
! Written by Romelia Salomon-Ferrer, 2/2012
!
!*******************************************************************************

module amd_mod

    use file_io_dat, only : MAX_FN_LEN

    implicit none

! Everything is private by default

! Variables
!

!AMD
    integer, save             :: iamd, iamdlag
    _REAL_, save              :: EthreshD, alphaD, EthreshP, alphaP
    _REAL_, save              :: totalenergy, tboostall, fwgt, totdih, tboost, dihsum
    _REAL_, save              :: amd_dih_noH, amd_dih_H, fwgtd
    integer, save                       :: num_amd_recs, num_amd_lag
    _REAL_, allocatable, save             :: amd_weights_and_energy(:, :)

! Files for AMD output:
!amdlog_unit        - records information for reweighting for each step.
    integer, parameter :: amdlog_unit = 77

! strings:
! amdlog          - AMD log filename
    character(len=MAX_FN_LEN), save :: amdlog

contains

!*******************************************************************************
!
! Subroutine: amd_setup
!
! Description: Sets up the AMD run. Opens the log file, sets up exchange tables
!              etc.
!
!*******************************************************************************

    subroutine amd_setup(ntwx)

        implicit none
#  include "parallel.h"

! Formal arguments:
        integer               :: ntwx

! Local variables:
        integer :: alloc_failed

        allocate (amd_weights_and_energy(6, ntwx), stat=alloc_failed)
        REQUIRE(alloc_failed == 0)

        fwgt = 1.0
        fwgtd = 1.0
        tboostall = 0.0
        tboost = 0.0
        num_amd_recs = 0
        num_amd_lag = 0

        ! Open and write out header to amdlog file. Only overall master
        !  deals with the amdlog.
        if (worldrank == 0) then
            call amopen(amdlog_unit, amdlog, 'U', 'F', 'W')
            write (amdlog_unit, '(a)') "# Accelerated Molecular Dynamics log file"
            write (amdlog_unit, '(a)') "# ntwx,total_nstep,istep+i,tot_potenergy, &
            &            totdih_ene,fwgt,fwgtd,tboostall,tboost"
        end if

        return
    end subroutine amd_setup

!*******************************************************************************
!
! Subroutine:  calculate_amd_dih_weights
!
! Description: <TBS>
!
!*******************************************************************************

    subroutine calculate_amd_dih_weights(totdih_ene, temp0)

        use constants, only : KB
        implicit none

! Formal arguments:
        _REAL_              :: totdih_ene
        _REAL_              :: temp0

! Calculate the boosting weight for amd

        totdih = 0.0d0
        tboost = 0.0d0
        fwgtd = 1.0d0
        totdih = totdih_ene

        if (totdih .le. EthreshD) then
            if (num_amd_lag .eq. 0) then
                tboost = ((EthreshD - totdih)**2)/ &
                    ((alphaD + (EthreshD - totdih))*temp0*KB)
                fwgtd = (alphaD**2)/((alphaD + EthreshD - totdih)**2)
            end if
        end if

        return

    end subroutine calculate_amd_dih_weights

!*******************************************************************************
!
! Subroutine:  calculate_amd_total_weights
!
! Description: <TBS>
!
!*******************************************************************************

    subroutine calculate_amd_total_weights(atm_cnt, tot_potenergy, totdih_ene, frc, temp0)

        use constants, only : KB
        implicit none
#  include "parallel.h"

! Formal arguments:
        integer             :: atm_cnt
        _REAL_              :: tot_potenergy
        _REAL_              :: totdih_ene
        _REAL_              :: frc(*)
        _REAL_              :: temp0

        integer             :: i

!AMD DUAL BOOST CALC START
        if (iamd .gt. 0) then
            totalenergy = 0.0d0
            tboostall = 0.0d0
            fwgt = 1.0d0
            totalenergy = tot_potenergy + (tboost*temp0*KB)
            if (((iamd == 1) .or. (iamd == 3)) .and. (totalenergy .le. EthreshP)) then
                if (num_amd_lag .eq. 0) then
                    tboostall = ((EthreshP - totalenergy)**2)/ &
                        ((alphaP + (EthreshP - totalenergy))*temp0*KB)
                    fwgt = (alphaP**2)/((alphaP + EthreshP - totalenergy)**2)

                    do i = 1, atm_cnt
                        frc(i*3 - 2) = frc(i*3 - 2)*fwgt
                        frc(i*3 - 1) = frc(i*3 - 1)*fwgt
                        frc(i*3) = frc(i*3)*fwgt
                    end do
                end if
            end if
            if (worldrank == 0 .and. (num_amd_recs .gt. 0)) then
                amd_weights_and_energy(1, num_amd_recs) = tot_potenergy
                amd_weights_and_energy(2, num_amd_recs) = totdih_ene
                amd_weights_and_energy(3, num_amd_recs) = fwgt
                amd_weights_and_energy(4, num_amd_recs) = fwgtd
                amd_weights_and_energy(5, num_amd_recs) = tboostall
                amd_weights_and_energy(6, num_amd_recs) = tboost
            end if
            num_amd_recs = num_amd_recs + 1
            if (num_amd_lag .eq. iamdlag) then
                num_amd_lag = 0
            else
                num_amd_lag = num_amd_lag + 1
            end if
        end if

        return

    end subroutine calculate_amd_total_weights

!*******************************************************************************
!
! Subroutine:  write_amd_weights
!
! Description: <TBS>
!
!*******************************************************************************

    subroutine write_amd_weights(ntwx, total_nstep)

        implicit none

! Formal arguments:

        integer               :: ntwx
        integer               :: total_nstep

! Local variables:

        integer               :: i, istep

        istep = total_nstep - ntwx
        do i = 1, ntwx
            write (amdlog_unit, '(2x,3i10,6f22.12)'), ntwx, total_nstep, istep + i, &
                amd_weights_and_energy(1, i), amd_weights_and_energy(2, i), &
                amd_weights_and_energy(3, i), amd_weights_and_energy(4, i), &
                amd_weights_and_energy(5, i), amd_weights_and_energy(6, i)

        end do

        num_amd_recs = 1

        return

    end subroutine write_amd_weights

!*********************************************************************
!               SUBROUTINE AMD_CLEANUP
!*********************************************************************
!  Close AMD files
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    subroutine amd_cleanup()

        implicit none
#  include "parallel.h"

        integer ier

        ! Close amdlog file
        if (worldrank == 0) then
            close (amdlog_unit)
        end if

! dealocate AMD vectors
        deallocate (amd_weights_and_energy, stat=ier)
        REQUIRE(ier == 0)

    end subroutine amd_cleanup

end module amd_mod
