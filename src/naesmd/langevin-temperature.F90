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
!
!      subroutine sdterm (pfric,vfric,afric,prand,vrand)
module langevin_temperature
use naesmd_constants
use random
!implicit none

contains

      subroutine sdterm 

        implicit none

        include 'sizes'
        include 'common'

        integer i,j
        double precision  xgamma(nmax) 
        double precision  gdt,egdt 
        double precision  pterm,vterm 
        double precision  rho,rhoc 
        double precision  gdt2,gdt3,gdt4,gdt5,gdt6,gdt7,gdt8,gdt9 
        double precision  ktm 
        double precision  psig,vsig 
        double precision  pnorm,vnorm 
        ! double precision  normal 

!
!     set the atomic friction coefficients to the global value
!
      do i = 1, natom
         xgamma(i) = friction
      end do
!
!     get the frictional and random terms for stochastic dynamics
!
      do i = 1, natom
            gdt = xgamma(i) * dtmdqt
!
!     stochastic dynamics reduces to simple MD for zero friction
!
            if (gdt .le. 0.0d0) then
               pfric(i) = 1.0d0
               vfric(i) = dtmdqt
               afric(i) = 0.5d0 * dtmdqt * dtmdqt
               do j = 1, 3
                  prand(j,i) = 0.0d0
                  vrand(j,i) = 0.0d0
               end do

!
!     analytical expressions when friction coefficient is large
!
            else
               if (gdt .ge. 0.05d0) then
                  egdt = exp(-gdt)
                  pfric(i) = egdt
                  vfric(i) = (1.0d0-egdt) / xgamma(i)
                  afric(i) = (dtmdqt-vfric(i)) / xgamma(i)
                  pterm = 2.0d0*gdt - 3.0d0 + (4.0d0-egdt)*egdt
                  vterm = 1.0d0 - egdt**2
                  rho = (1.0d0-egdt)**2 / sqrt(pterm*vterm)
!
!     use series expansions when friction coefficient is small
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
                  afric(i) = (gdt2/2.0d0 - gdt3/6.0d0 + gdt4/24.0d0 &
                                - gdt5/120.0d0 + gdt6/720.0d0 &
                                - gdt7/5040.0d0 + gdt8/40320.0d0 &
                                - gdt9/362880.0d0) / xgamma(i)**2
                  vfric(i) = dtmdqt - xgamma(i)*afric(i)
                  pfric(i) = 1.0d0 - xgamma(i)*vfric(i)
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
!     compute random terms to thermostat the nonzero friction case
!
               ktm = boltzman * temp0 / massmdqt(i)
               psig = sqrt(ktm*pterm) / xgamma(i)
               vsig = sqrt(ktm*vterm)
               rhoc = sqrt(1.0d0 - rho*rho)
               do j = 1, 3
                  pnorm = normal ()
                  vnorm = normal ()
                  prand(j,i) = psig * pnorm
                  vrand(j,i) = vsig * (rho*pnorm+rhoc*vnorm)
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
!
      function normal ()

      implicit none

      include 'sizes'
      include 'common'

      double precision v1,v2,rsq,iseedhop
      double precision factor,store,normal
      ! double precision rranf1

      logical compute
      save compute,store
      data compute  / .true. /
!
!
!     get a pair of random values from the distribution
!
      if (compute) then
   10    continue
!         call random(iseedmdqt,iseedhop)
!         v1 = 2.0d0 * iseedhop - 1.0d0
!         call random(iseedmdqt,iseedhop)
!         v2 = 2.0d0 * iseedhop - 1.0d0
         v1 = 2.0d0 * rranf1(iseedmdqt) - 1.0d0
         v2 = 2.0d0 * rranf1(iseedmdqt) - 1.0d0
!         write(66,*) rranf1(iseedmdqt)
!         call flush(66)
         rsq = v1**2 + v2**2
         if (rsq .ge. 1.0d0)  goto 10
         factor = sqrt(-2.0d0*log(rsq)/rsq)
         store = v1 * factor
         normal = v2 * factor
         compute = .false.
!
!     use the second random value computed at the last call
!
      else
         normal = store
         compute = .true.
      end if
      return
      end function

! Subroutine to thermaize the velocities

        SUBROUTINE temperature(i) 

        IMPLICIT NONE

        integer i,j
        double precision scltmp
        include 'sizes'
        include 'common'


       call temper

! tempf is the target temperature at the specific heating step
! temp0 is the final target temperature
! tempi is the instantaneus temperature

       if(ensemble.eq.'temper') then
           if(prep.eq.'heat') then
               if(i.eq.1) then
                  iconttemperature=1
                  tempf=tempi
               endif
               if(iconttemperature.le.istepheat) then
                  iconttemperature=iconttemperature+1
                  if (iconttemperature.eq.istepheat) then
                      iconttemperature=1
                      tempf=tempf+1
                      if(tempf.gt.temp0) tempf=temp0
                  endif
               endif
           else 
               tempf=temp0
           endif

           if(prep.eq.'heat') then
               SCLTMP = DSQRT(1.0D0+dtmdqt*CONVT/TAO* &
      (TEMPf/TEMPI-1.0D0))
           else
               SCLTMP = DSQRT(1.0D0+dtmdqt*CONVT/TAO* &
      (TEMP0/TEMPI-1.0D0))
           endif
           do j=1,natom
               vx(j)=scltmp*vx(j)
               vy(j)=scltmp*vy(j)
               vz(j)=scltmp*vz(j)
           enddo
       else
           if(ensemble.eq.'langev') tempf=temp0 
       endif

 
      RETURN
      END SUBROUTINE


! SUBROUTINE THAT COMPUTES THE TEMPERATURE OF THE SYSTEM FROM 
! THE TOTAL KINETIC ENERGY, USING EQUIPARTITION THEOREM:
! TOTAL KINETIC ENERGY=1/2*KBOLTZ*TEMP*(DEGREES OF FREEDOM)
! THE DEGREES OF FREEDOM OF OUR SYSTEM, COMPOSED OF N ATOMS,
! WITH CONSTANT TOTAL LINEAR MOMENTUM AND WITH TOTAL ANGULAR
! MOMENTUM, ARE 3*N

   subroutine temper
   implicit none

   integer i
   _REAL_ xkboltz
   _REAL_ xkinsum
        
   include 'sizes'
   include 'common'

   xkboltz=3.166829662D-6

! KINETIC ENERGY OF ATOMS

   xkinsum=0.0d0
   do i=1,natom
      xkinsum=xkinsum+massmdqt(i)*(vx(i)**2+vy(i)**2+vz(i)**2)
   end do

   tempi=xkinsum/(3.0d0*dfloat(natom)*xkboltz)

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
      real*8 rrang
      logical flag

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
