#include "dprec.fh"
#include "assert.fh"

!
!     #############################################################
!     ##                                                         ##
!     ##  subroutine sdterm  --  frictional and random SD terms  ##
!     ##                                                         ##
!     #############################################################
!
!
!     "sdterm" gets frictional and random force terms needed to
!     update positions and velocities via stochastic dynamics
!
module langevin_temperature
    use naesmd_constants
    use random

contains

    subroutine sdterm(naesmd_struct)
        use naesmd_module, only:naesmd_structure
        implicit none
        type(naesmd_structure), intent(inout) :: naesmd_struct 

        integer i,j
        _REAL_  xgamma(naesmd_struct%natom) 
        _REAL_  gdt,egdt 
        _REAL_  pterm,vterm 
        _REAL_  rho,rhoc 
        _REAL_  gdt2,gdt3,gdt4,gdt5,gdt6,gdt7,gdt8,gdt9 
        _REAL_  ktm 
        _REAL_  psig,vsig 
        _REAL_  pnorm,vnorm 
        !
        !     set the atomic naesmd_struct%friction coefficients to the global value
        !
        do i = 1, naesmd_struct%natom
            xgamma(i) = naesmd_struct%friction
        end do
        !
        !     get the frictional and random terms for stochastic dynamics
        !
        do i = 1, naesmd_struct%natom
            gdt = xgamma(i) * naesmd_struct%dtmdqt
            !
            !     stochastic dynamics reduces to simple MD for zero naesmd_struct%friction
            !
            if (gdt .le. 0.0d0) then
                naesmd_struct%pfric(i) = 1.0d0
                naesmd_struct%vfric(i) = naesmd_struct%dtmdqt
                naesmd_struct%afric(i) = 0.5d0 * naesmd_struct%dtmdqt * naesmd_struct%dtmdqt
                do j = 1, 3
                    naesmd_struct%prand(j,i) = 0.0d0
                    naesmd_struct%vrand(j,i) = 0.0d0
                end do

            !
            !     analytical expressions when naesmd_struct%friction coefficient is large
            !
            else
                if (gdt .ge. 0.05d0) then
                    egdt = exp(-gdt)
                    naesmd_struct%pfric(i) = egdt
                    naesmd_struct%vfric(i) = (1.0d0-egdt) / xgamma(i)
                    naesmd_struct%afric(i) = (naesmd_struct%dtmdqt-naesmd_struct%vfric(i)) / xgamma(i)
                    pterm = 2.0d0*gdt - 3.0d0 + (4.0d0-egdt)*egdt
                    vterm = 1.0d0 - egdt**2
                    rho = (1.0d0-egdt)**2 / sqrt(pterm*vterm)
                !
                !     use series expansions when naesmd_struct%friction coefficient is small
                !
                else
                    gdt2 = gdt * gdt
                    gdt3 = gdt * gdt2
                    gdt4 = gdt2 * gdt2
                    gdt5 = gdt2 * gdt3
                    gdt6 = gdt3 * gdt3
                    gdt7 = gdt3 * gdt4
                    gdt8 = gdt4 * gdt4
                    gdt9 = gdt4 * gdt5
                    naesmd_struct%afric(i) = (gdt2/2.0d0 - gdt3/6.0d0 + gdt4/24.0d0 &
                        - gdt5/120.0d0 + gdt6/720.0d0 &
                        - gdt7/5040.0d0 + gdt8/40320.0d0 &
                        - gdt9/362880.0d0) / xgamma(i)**2
                    naesmd_struct%vfric(i) = naesmd_struct%dtmdqt - xgamma(i)*naesmd_struct%afric(i)
                    naesmd_struct%pfric(i) = 1.0d0 - xgamma(i)*naesmd_struct%vfric(i)
                    pterm = 2.0d0*gdt3/3.0d0 - gdt4/2.0d0 &
                        + 7.0d0*gdt5/30.0d0 - gdt6/12.0d0 &
                        + 31.0d0*gdt7/1260.0d0 - gdt8/160.0d0 &
                        + 127.0d0*gdt9/90720.0d0
                    vterm = 2.0d0*gdt - 2.0d0*gdt2 + 4.0d0*gdt3/3.0d0 &
                        - 2.0d0*gdt4/3.0d0 + 4.0d0*gdt5/15.0d0 &
                        - 4.0d0*gdt6/45.0d0 + 8.0d0*gdt7/315.0d0 &
                        - 2.0d0*gdt8/315.0d0 + 4.0d0*gdt9/2835.0d0
                    rho = sqrt(3.0d0) * (0.5d0 - 3.0d0*gdt/16.0d0 &
                        - 17.0d0*gdt2/1280.0d0 &
                        + 17.0d0*gdt3/6144.0d0 &
                        + 40967.0d0*gdt4/34406400.0d0 &
                        - 57203.0d0*gdt5/275251200.0d0 &
                        - 1429487.0d0*gdt6/13212057600.0d0)
                end if
                !
                !     compute random terms to thermostat the nonzero naesmd_struct%friction case
                !
                ktm = boltzman * naesmd_struct%temp0 / naesmd_struct%massmdqt(i)
                psig = sqrt(ktm*pterm) / xgamma(i)
                vsig = sqrt(ktm*vterm)
                rhoc = sqrt(1.0d0 - rho*rho)
                do j = 1, 3
                    pnorm = normal (naesmd_struct)
                    vnorm = normal (naesmd_struct)
                    naesmd_struct%prand(j,i) = psig * pnorm
                    naesmd_struct%vrand(j,i) = vsig * (rho*pnorm+rhoc*vnorm)
                end do
            end if
        end do
        return
    end subroutine

    !
    !
    !     ############################################################
    !     ##                                                        ##
    !     ##  function normal  --  random number from normal curve  ##
    !     ##                                                        ##
    !     ############################################################
    !
    !
    !     "normal" generates a random number from a normal Gaussian
    !     distribution with a mean of zero and a variance of one
    !
    function normal (naesmd_struct)
        use naesmd_module, only:naesmd_structure

        implicit none
        type(naesmd_structure), intent(inout) :: naesmd_struct 

        _REAL_ v1,v2,rsq
        _REAL_ factor,normal
        !
        !
        !     get a pair of random values from the distribution
        !
        if (naesmd_struct%compute) then
10      continue
        v1 = 2.0d0 * rranf1(naesmd_struct%iseedmdqt) - 1.0d0
        v2 = 2.0d0 * rranf1(naesmd_struct%iseedmdqt) - 1.0d0
        rsq = v1**2 + v2**2
        if (rsq .ge. 1.0d0)  goto 10
        factor = sqrt(-2.0d0*log(rsq)/rsq)
        naesmd_struct%store = v1 * factor
        normal = v2 * factor
        naesmd_struct%compute = .false.
    !
    !     use the second random value computed at the last call
    !
    else
        normal = naesmd_struct%store
        naesmd_struct%compute = .true.
    end if
    return
end function

! Subroutine to thermaize the velocities

SUBROUTINE temperature(i,naesmd_struct)
    use naesmd_module, only:naesmd_structure

    IMPLICIT NONE
    type(naesmd_structure), intent(inout) :: naesmd_struct 

    integer i,j
    _REAL_ scltmp

    call temper(naesmd_struct)


        if(naesmd_struct%ensemble.eq.'langev') naesmd_struct%tempf=naesmd_struct%temp0

 
    RETURN
END SUBROUTINE


! SUBROUTINE THAT COMPUTES THE TEMPERATURE OF THE SYSTEM FROM 
! THE TOTAL KINETIC ENERGY, USING EQUIPARTITION THEOREM:
! TOTAL KINETIC ENERGY=1/2*KBOLTZ*TEMP*(DEGREES OF FREEDOM)
! THE DEGREES OF FREEDOM OF OUR SYSTEM, COMPOSED OF N ATOMS,
! WITH CONSTANT TOTAL LINEAR MOMENTUM AND WITH TOTAL ANGULAR
! MOMENTUM, ARE 3*N

subroutine temper(naesmd_struct)
    use naesmd_module, only:naesmd_structure

    implicit none
    type(naesmd_structure), intent(inout) :: naesmd_struct 

    integer i
    _REAL_ xkboltz
    _REAL_ xkinsum
        

    xkboltz=3.166829662D-6

    ! KINETIC ENERGY OF ATOMS

    xkinsum=0.0d0
    do i=1,naesmd_struct%natom
        xkinsum=xkinsum+naesmd_struct%massmdqt(i)*(naesmd_struct%vx(i)**2+naesmd_struct%vy(i)**2+naesmd_struct%vz(i)**2)
    end do

    naesmd_struct%tempi=xkinsum/(3.0d0*dfloat(naesmd_struct%natom)*xkboltz)

    return
end subroutine
!
!********************************************************************
!
!********************************************************************
!
function rranf1(iseed)
    implicit real*8 (a-h,o-z)
    real*8 :: rranf1

    ! rranf() generates uniform in [0..1)
    ! NB: Seed = 0 is forbidden.
    ! Extremely fast.  232-1 cycle.
    ! Source taken from Numerical Recipes;

    parameter (ia = 16807)
    parameter (im = 2147483647)
    parameter (iq = 127773)
    parameter (ir = 2836)
    parameter (am = 1d0/im)
    parameter (am2 = 2*am)

    k = iseed / iq
    iseed = ia*(iseed-k*iq) - k*ir
    if (iseed.lt.0) iseed = iseed + im
    rranf1 = am*iseed
    return

end function
end module
