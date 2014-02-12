! <compile=optimized>
module nmr

  implicit none

#include "assert.fh"
#include "copyright.h"
#include "dprec.fh"

public   chklin, echoin, impnum, modwt, ndvprt, nmrnrg, nmrprt, nmrrad, &
         nmrred, nmrsht, parsrest
private  angnrg, at2res, aveint, commass, disnrg, getnat, gendisnrg, jnrg, &
         nmrcmf, nmrcms, nmrgrp, parenfind, parscom, parsgen, plnnrg, plptnrg, &
         r6ave, r6drv, readmask, tornrg, vec_constr, smorse

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine angnrg here]
subroutine angnrg(x,f,dfr,thet,i1,i2,i3,e,r1,r2,r3,r4,k2,k3,ntb, &
      aave,aave0,nave,nexact,ipower,tauave,ravein, &
      dt,navint,iave,incflg,iravin,iflag)
   
   ! Subroutine ANGle eNeRGy
   
   ! This subroutine calculates the valence angle restraint energy.
   ! The energy function is a flat-bottomed well, with parabolic sides
   ! which become linear at greater distances.
   
   ! If r is the valence angle:
   
   !       r <= r1  : E = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
   ! r1 <= r <= r2  : E = k2*(r-r2)**2
   ! r2 <= r <= r3  : E = 0
   ! r3 <= r <= r4  : E = k3*(r-r3)**2
   ! r4 <= r        : E = 2*k3*(r4-r3)*(r-r4) + K3*(r4-r3)**2
   
   
   ! IAVE = 0: then the instantaneous (current) value of r is used
   !           in all calculations.
   ! IAVE = 1: The time-averaged value of r is used, and forces are calculated
   !           according to the above function, i.e.
   !               dE/dx = (dE/dr_ave)(dr_ave/dr)(dr/dx).
   ! IAVE = 2: The time-averaged value of r is used, and forces are calculated as
   !               dE/dx = (dE/dr_ave)(dr/dx).
   !           This expression, if integrated would not yield the energy function
   !           above, so energies reported are pseudo-energies. Also, it is a
   !           non-conservative force, so the system will tend to heat up. But it
   !           avoids potentially large force terms with (1+IPOWER)
   !           exponentiation, which can cause the simulation to become unstable.
   
   ! Other input:
   !    X(I)  : Coordinate array.
   !    F(I)  : Force array;
   !            Modified on output to include forces from this restraint.
   !  I1,I2,I3: Atom pointers for this angle (3*(I-1), where I is absolute
   !            atom number.
   !    NTB   : Periodic boundary conditions flag.
   !    IFLAG : =0, then calculate energy and derivatives
   !            =1, then only calculate current value of angle
   !            =2, the calculate current value of angle and energy, but
   !                do not accumulate derivatives
   !            =3, then calculate energy and derivatives, but do not
   !                accumulate forces into F() array. They are returned
   !                in the DFR array.
   !   INCFLG : Determines whether local saved pointers are updated in AVEINT
   
   ! (the following only used when time-averaging):
   ! AAVE(I) : Contains time-averaged values of angles, using TAUAVE damping.
   ! AAVE0(I): Contains time-averaged values of angles, using no damping.
   
   !           Both AAVE and AAVE0 are passed so that the first element of
   !           the array corresponds to the average value of the angle currently
   !           of interest, and the second element of the array contains
   !           the previous value of the angle of interest.
   
   !     NAVE : > 0 for time-averaging.
   !   NEXACT : Not currently used.
   !   IPOWER : The power used in averaging the angle data.
   !   TAUAVE : The exponential decay factor used in averaging the angle data
   !   IRAVIN : =0 do not modify the values of the angles by RAVEIN
   !            =1 modify theta by RAVEIN, as shown below.
   !   RAVEIN : At time 0, time-average integral is undefined.
   !            If -1000<RAVEIN<1000, use r_initial+ravein on this step.
   !            If RAVEIN<=-1000, use r_target+(ravein+1000) on this step.
   !            If RAVEIN>=1000, use r_target-(ravein-1000) on this step.
   !       DT : The integration time-step (ps)
   !   NAVINT : Time-averaged restraints are only updated every NAVINT
   !            steps. Typically, NAVINT=1.
   
   ! Output:
   !    E     : The energy contribution from this restraint.
   !    DFR   : If IFLAG=3, DFR(1->3) contains the d_E/d_distance forces
   !            for x,y, and z for position 1; those for position 2 are in
   !            DFR(4->6), and so on for each of the 5 atoms here
   !            Note that the IFLAG designation for this is different from
   !            that for disnrg, plnnrg, and plptnrg
   !  THET    : The calculated valence angle
   
   ! Author: David A. Pearlman
   ! Date: 7/89
   use constants, only : DEG_TO_RAD, half
   implicit none
   integer:: i1, i2, i3, iave, iflag, iflg, incflg, ipower, irav, &
        iravin, m, nave, navint, nexact, ntb
   _REAL_ :: aave, aave0, cii, cik, ckk, cst, denom, df, dfr, dif, &
        dif1, dravdr, dt, e, f, r1, r2, r3, r4, ravein, rdenom, rij, &
        rij2, rinc, rkj, rkj2, rnow, small, small2, snt, st, sth, tauave, &
        thet, x, xij, xkj
   _REAL_ k2,k3
   parameter (small=1.0d-14)
   parameter (small2=1.0d-5)
   dimension x(*),f(*),dfr(*),aave(2),aave0(2)
   dimension xij(3),xkj(3)
   rinc = 1000.0d0
   
   iflg = 1
   if (iflag == 1) iflg = 2
   
   ! Calculate the angle:
   
   rij2 = 0.0d0
   rkj2 = 0.0d0
   
   do m=1,3
      xij(m) = x(i1+m) - x(i2+m)
      xkj(m) = x(i3+m) - x(i2+m)
      rij2 = rij2 + xij(m)**2
      rkj2 = rkj2 + xkj(m)**2
   end do
   rij = sqrt(rij2)
   rkj = sqrt(rkj2)
   
   ! Calculate angle from a*b/|ab| = cos(theta)
   
   rdenom = rij*rkj
   cst = (xij(1)*xkj(1)+xij(2)*xkj(2)+xij(3)*xkj(3))/rdenom
   if (cst > 1.0d0) cst = 1.0d0
   if (cst < -1.0d0) cst = -1.0d0
   thet = acos(cst)
   
   ! Get averaged value, if requested
   
   dravdr = 1.0d0
   if (iave > 0) then
      rnow = thet
      ! Set the value on the first step (where integral not defined):
      if (iravin == 1) then
         irav = int(ravein)
         if (-1000 < irav .and. irav < 1000) then
            rnow = rnow + ravein*DEG_TO_RAD
         else
            rinc = -sign(rinc,ravein)
            rinc = (ravein+rinc) * DEG_TO_RAD
            if (r2 > small2) rnow =((r2+r3)*half) + rinc
            if (r2 <= small2) rnow = r3 + rinc
         end if
      end if
      call aveint(aave,aave0,rnow,tauave,dt,navint,ipower,2, &
            incflg,iflg,thet,denom)
      
      ! DRAVDR is the factor d(r_ave)/d(r_current)
      
      if (iave == 1 .and. iflag /= 1) &
            dravdr = ((thet/rnow)**(1+ipower))*dt/denom
   end if
   
   if (iflag == 1) return
   
   ! Calculate energy (E) and the derivative with respect to the valence
   ! angle (DF):
   
   if (thet < r1) then
      dif1 = r1-r2
      df = 2.0d0 * k2 * dif1
      e = df * (thet-r1) + k2*dif1*dif1
   else if (thet < r2) then
      dif = thet - r2
      df = 2.0d0 * k2 * dif
      e = k2*dif*dif
   else if (thet <= r3) then
      e = 0.0d0
      return
   else if (thet < r4) then
      dif = thet - r3
      df = 2.0d0 * k3 * dif
      e = k3*dif*dif
   else
      dif1 = r4-r3
      df = 2.0d0 * k3 * dif1
      e = df * (thet-r4) + k3*dif1*dif1
   end if
   if (iflag == 2) return
   
   ! Calculate the derivaties with respect to the coordinates, and add them
   ! into the force arrays. DRAVDR contains d(thet_ave)/d(thet(t)) if IAVE=1,
   ! 1.0 otherwise.
   
   df = dravdr * df
   snt = sin(thet)
   if (abs(snt) < small) snt = small
   st = -df/snt
   sth = st*cst
   cik = st/rdenom
   cii = sth/rij2
   ckk = sth/rkj2
   
   dfr(1) = cik*xkj(1) - cii*xij(1)
   dfr(2) = cik*xkj(2) - cii*xij(2)
   dfr(3) = cik*xkj(3) - cii*xij(3)
   dfr(7) = cik*xij(1) - ckk*xkj(1)
   dfr(8) = cik*xij(2) - ckk*xkj(2)
   dfr(9) = cik*xij(3) - ckk*xkj(3)
   dfr(4) = -dfr(1)-dfr(7)
   dfr(5) = -dfr(2)-dfr(8)
   dfr(6) = -dfr(3)-dfr(9)
   
   if (iflag == 3) return
   
   f(i1+1) = f(i1+1) - dfr(1)
   f(i1+2) = f(i1+2) - dfr(2)
   f(i1+3) = f(i1+3) - dfr(3)
   f(i2+1) = f(i2+1) - dfr(4)
   f(i2+2) = f(i2+2) - dfr(5)
   f(i2+3) = f(i2+3) - dfr(6)
   f(i3+1) = f(i3+1) - dfr(7)
   f(i3+2) = f(i3+2) - dfr(8)
   f(i3+3) = f(i3+3) - dfr(9)
   
   return
end subroutine angnrg 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine at2res here]
subroutine at2res(iat,ipres,nres,ires,iflag,iout)

   ! Subroutine ATom 2 RESidue number
   
   ! This routine takes the atom number IAT, and returns the number IRES of
   ! the residue containing it.
   
   ! Author: David A. Pearlman
   ! Date: 11/89
   
   ! IAT     : Input atom number
   ! IPRES(I): IPRES(I) is the first atom of residue I. IPRES(I+1)-1 is the
   !           last atom of residue I.
   ! NRES    : The number of residues in the system
   ! IRES    : The residue containing atom IAT (output).
   ! IFLAG   : =0 IAT is absolute atom number
   !           >1 IAT is (3*atom_number)-1
   ! IOUT    : Unit for informational prints
   
   implicit none
   integer:: i, iat, icheck, iflag, iout, ipres, ires, nres
   dimension ipres(*)
   
   icheck = iat
   if (iflag > 0) icheck = (iat/3) + 1
   do i=1,nres
      if (icheck >= ipres(i) .and. icheck < ipres(i+1)) then
         ires = i
         return
      end if
   end do
   
   write(iout,9001) icheck,ipres(nres+1)-1
   call mexit(iout, 1)
   9001 format(' Error from routine AT2RES: Atom number = ',i6, &
         ' cannot be translated',/,t29,'to an appropriate ', &
         'residue number. Last atom in',/,t29,'the system is ',i6)
end subroutine at2res 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Evaluate time-averaging integrals.
subroutine aveint(xave,xave0,xnow,tau,dt,navint,ipower,itype, &
      incflg,iflag,xaver,denom)
   
   
   ! Subroutine AVErage INTernals.
   
   ! This routine evaluates and updates the integrals needed to perform
   ! time-averaging
   
   ! Author: David A. Pearlman
   ! Date: 5/90
   ! Converted to Fortran 90; SRB Aug 2004
   
   ! Algorithm:
   !     Internals are averaged using inverse weighting to the
   ! IPOWER power (thus IPOWER = -1 gives straight linear weighting). In
   ! addition, an exponential decay function is included, so that the
   ! weighted average will continue to reflect the most recent values
   ! of the internals--even during long simulations.
   
   ! The formula for the time-averaged values is:
   
   !                  t
   !   rbar = 1/C * {int (exp(tp-t)/tau) r(tp)**-ipower dtp} **-1/ipower   (1)
   !                  0
   
   ! where int = integral (from 0-> t) over tp
   !      rbar = average value of internal
   !         t = current time
   !       tau = exponential decay constant
   !     r(tp) = value of internal at time tp
   !    ipower = average is over internals to the inverse ipower.
   !             Usually ipower = 3 for NOE distances, -1 for angles and
   !             torsions
   !         C = normalization integral.
   
   ! To minimize memory storage, we evaluate the integral as
   
   !    t
   !   int (exp(tp-t)/tau) r(tp)**-ipower dtp =
   !    0
   
   !                   t-dt
   !      exp(-dt/tau) int (exp(tp-t)/tau) r(tp)**-ipower dtp +
   !                    0
   
   !                    t
   !                   int (exp(tp-t)/tau) r(tp)**-ipower dtp
   !                   t-dt
   
   ! INPUT:
   !  XAVE(I)  : On input: XAVE(1) contains the value of the time-integral
   !                       shown above for the bounds 0->t-dt..
   !                       XAVE(2) contains r(t-dt)**-ipower
   !             On output (if IFLAG.NE.2)
   !                       XAVE(1) contains the value of the time-integral
   !                       shown above for the bounds 0->t.
   !                       XAVE(2) contains r(t)**-ipower
   
   ! XAVE0(I)  : Same as XAVE(1), except that the integral stored does not
   !             reflect exponential weighting (a real time-average).
   !  XNOW     : Current value of the internal on this call.
   !  TAU      : Characteristic time for the exponential decay.
   !             If TAU >= 1.D+6, then the non-exponential-damped values
   !             in XAVE0() will be used.
   !  DT       : The moleular dynamics timestep.
   !  NAVINT   : The current value of the internal is only stored every
   !             NAVINT steps (typically 1).
   !  IPOWER   : The exponent used in creating the weighted average. Standard
   !             values would be 3 for distances (resulting in **-3 weighting)
   !             and -1 for angles (resulting in **1/non-exponential weighting).
   !  ITYPE    : = 1 for bonds; =2 for angles; =3 for torsions (cos(tau));
   !             = 4 for torsions (sin(tau)); =5 for plane-point;
   !             = 6 for plane-plane; =7 for generalized distance 
   ! INCFLG    : Determines whether local saved pointers are updated
   !             on this call.
   !             INCFLG = 0 -- pointers are _not_ updated.
   !             INCFLG = 1 -- pointers are updated.
   !             Pointers should only be updated on the first call to AVEINT
   !             for each restraint type (bond, angle, torsion) in any call
   !             to NMRNRG.
   !  IFLAG    : Determines behavior on this call.
   !             IFLAG <= 1 -- Update the integral to reflect the value in XNOW,
   !                           set XAVE(2) to XNOW, and return
   !                           the calculated averaged value XAVE/DENOM.
   !             IFLAG  = 2 -- Just return the calculated values XAVE/DENOM
   
   ! OUTPUT:
   !  XAVER    : The calculated time-averaged value for the internal.
   !  DENOM    : The value of the exponential weight-factor divisor integral
   use constants, only : ZERO, one, half 
   implicit none

   _REAL_  xave(2), xave0(2)
   _REAL_  xnow, tau, dt
   integer navint, ipower, itype, incflg, iflag
   _REAL_  xaver, denom
   
   integer, parameter :: NUMBER_OF_TYPES = 7
   integer, save :: nsteps(NUMBER_OF_TYPES) = 0
   integer, save :: icalls(NUMBER_OF_TYPES) = 0
   ! nsteps(itype) : The number of calls with IFLAG.LE.1 (i.e. calls
   !                 where the integral sum was updated).
   ! icalls(itype) : Number of calls with IFLAG=0 or 1.

   integer i
   _REAL_  dtu
   _REAL_  tp
   _REAL_  xnowp
   _REAL_  xpon
   
   
   ASSERT(itype >= 1 .and. itype <= NUMBER_OF_TYPES)

   ! DTU is the time-step (modified to reflect NAVINT if it is not 1)
   
   dtu = dt * navint
   
   ! This section stores the value
   
   if (iflag <= 1) then
      icalls(itype) = icalls(itype) + incflg
      
      if (mod(icalls(itype)-1,navint) == 0) then
         nsteps(itype) = nsteps(itype) + incflg
         
         ! Update the integral summation using Simpsons rule to calculate the
         ! integral over the last interval.
         
         ! On the first step, set the integral to zero.
         ! XAVE(1) stores the integral with exponential decay. XAVE0(1) stores
         ! the integral without the exponential decay weighting factor.
         
         xnowp = xnow**(-ipower)
         if (nsteps(itype) <= 1) then
            xave(1) = ZERO
            xave0(1) = ZERO
         else
            xave(1) = exp(-dtu/tau)*xave(1) + &
                  (exp(-dtu/tau)*xave(2) + xnowp) * (dtu*half)
            xave0(1) = xave0(1) + (xave(2) + xnowp) * (dtu*half)
         end if
         xave(2) = xnowp
      end if
   end if
   
   ! This section calculates the average:
   
   ! TP is the current time, modified to reflect NAVINT
   
   tp = (nsteps(itype)-1)*dtu
   
   ! On the first step, the current value is the average:
   
   if (nsteps(itype) <= 1) then
      xaver = xnow
      denom = ONE
      return
   else
      xpon = -ONE/ipower
      if (tau < 1.d+6) then
         denom = tau*(ONE - exp(-tp/tau))
         xaver = sign((abs(xave(1)/denom))**xpon,xave(1)/denom)
      else
         denom = tp
         xaver = sign((abs(xave0(1)/denom))**xpon,xave0(1)/denom)
      end if
   end if
   
   return
end subroutine aveint 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine chklin here]
subroutine chklin(line,mtchst,ilnmat,imatch,istrt,ilengt,iout)
   
   
   ! Subroutine CHecK LINe
   
   ! This subroutine checks the character string LINE for the syntax
   ! MTCHST = FLNME, where MTCHST is a passed character string, and
   ! FLNME can be any character string. The equals sign ("=") is required.
   
   ! Author: David A. Pearlman
   ! Date: 3/90
   
   ! INPUT:
   
   ! LINE: The character line to be searched (character(len=80))
   ! MTCHST: The match string to be searched for in LINE. (Character(len=10))
   ! ILNMAT: The length of MTCHST.
   ! IMATCH: =0 if there is no match.
   !         =1 if there is a match.
   ! ISTRT: The starting position of FLNME
   ! ILENGT: The length of FLNME.
   ! IOUT: Unit for error prints
   
   ! MAXLEN is the length of LINE
   implicit none
   integer:: i, ilengt, ilnmat, imatch, iout, istrt, j, maxlen
   parameter (maxlen = 80)
   character(len=80) line
   character(len=10) mtchst
   
   ! See if the match string MTCHST occurs before the "=" sign:
   
   do i=1,maxlen-ilnmat+1
      if (line(i:i+ilnmat-1) == mtchst) goto 8
      if (line(i:i) == '=') goto 6
   end do
   
   ! Get here if MTCHST not found before the "=" sign:
   
   6 imatch = 0
   return
   
   ! Get here if MTCHST found:
   
   8 imatch = 1
   do i=1,maxlen
      if (line(i:i) == '=') then
         do j=i+1,maxlen
            if (line(j:j) /= ' ') then
               istrt = j
               goto 30
            end if
         end do
      end if
   end do
   
   ! If we get here (instead of skipping to line 30), keyword was found but
   ! either "=" was missing, or there was no string following the "=":
   
   write(iout,9001) line
   call mexit(iout, 1)
   
   30 continue
   do i=maxlen,1,-1
      if (line(i:i) /= ' ') then
         ilengt = i-istrt+1
         return
      end if
   end do
   
   9001 format(' Error: Missing "=" and/or filename after keyword in ', &
         'line:',/,a)
end subroutine chklin 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine commass here]
subroutine commass(comrimass,nmrat,nmrcom,rimass,i)
   
   
   ! Subroutine Center Of Mass MASS array maker
   
   ! This routine generates a miniature COMRIMASS array and populates it
   ! with the inverse masses of either COM groupings or single atoms
   ! This function is only necessary for use with plane-point restraints in the event
   ! that COMs are defined within one
   
   
   ! Author: Matthew Seetin
   ! Date: 9/2007
   
   implicit none
   integer:: i, iat, ip, j, nmrat, nmrcom
   _REAL_ :: comrimass, rimass
   dimension comrimass(4),nmrat(16,*),nmrcom(2,*),rimass(*)
   
   do iat = 1,4
      comrimass(iat) = 0.0d0
   end do
   do iat = 1,4
      if (nmrat(iat,i) < 0) then
         do ip = -nmrat(iat,i),-nmrat(8+iat,i)
            do j = nmrcom(1,ip),nmrcom(2,ip)
               if (comrimass(iat) > 0.0d0) then
                 comrimass(iat) = 1.0/(1.0/comrimass(iat) + 1.0/rimass(j))
               else
                 comrimass(iat) = rimass(j)
               end if
            end do
         end do
      else
         comrimass(iat) = rimass((nmrat(iat,i)/3) + 1)
      end if
   end do
   
   return
end subroutine commass 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine disnrg here]
subroutine disnrg(x,f,dfr,rij,i1,i2,e,r1,r2,r3,r4,k2,k3,ntb, &
      bave,bave0,nave,nexact,ipower,tauave,ravein, &
      dt,navint,iave,incflg,iravin,iflag,ialtd)
   
   
   ! Subroutine DIStance eNeRGy
   
   ! This subroutine calculates the distance restraint energy.
   ! The energy function is a flat-bottomed well, with parabolic sides
   ! which become linear at greater distances.
   
   ! If the values in nmrat(1) or nmrat(2) are negative, then we do
   ! r**-6 averaging on the appropriate distances between the groups;
   ! group information is kept in nmrcom.
   
   ! If r is the bond length, and IALTD=1, then:
   
   !       use a modification of a functional form suggested by
   !       Michael Nilges, Protein Eng. 2: 27-38 (1988):
   
   !           r <  r2  : E = k2*(r-r2)**2
   !     r2 <= r <= r3  : E = 0
   !     r3 <= r <= r4  : E = k3*(r-r3)**2
   !     r4 <  r        : E = k3*[c*(r-r3) + b/(r-r3) + a]
   !                        where: a = 3*(r4-r3)**2 - 2*c*(r4-r3)
   !                               b = -2*(r4-r3)**3 + c*(r4-r3)**2
   !                               c determines whether the potential flattens
   !                                 at large distances (c=0) or becomes linear.
   !                                 Current code assumes c=0, allowing "bad"
   !                                 restraints.
   
   ! otherwise, if IALTD=0, use the function form from Amber4/Amber5:
   
   !           r <= r1  : E = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
   !     r1 <= r <= r2  : E = k2*(r-r2)**2
   !     r2 <= r <= r3  : E = 0
   !     r3 <= r <= r4  : E = k3*(r-r3)**2
   !     r4 <= r        : E = 2*k3*(r4-r3)*(r-r4) + K3*(r4-r3)**2
   
   
   
   ! IAVE = 0: then the instantaneous (current) value of r is used
   !           in all calculations.
   ! IAVE = 1: The time-averaged value of r is used, and forces are calculated
   !           according to the above function, i.e.
   !               dE/dx = (dE/dr_ave)(dr_ave/dr)(dr/dx).
   ! IAVE = 2: The time-averaged value of r is used, and forces are calculated as
   !               dE/dx = (dE/dr_ave)(dr/dx).
   !           This expression, if integrated would not yield the energy function
   !           above, so energies reported are pseudo-energies. Also, it is a
   !           non-conservative force, so the system will tend to heat up. But it
   !           avoids potentially large force terms with (1+IPOWER)
   !           exponentiation, which can cause the simulation to become unstable.
   
   ! Other input:
   !    X(I)  : Coordinate array.
   !    F(I)  : Force array;
   !            Modified on output to include forces from this restraint.
   !    I1,I2 : Atom pointers for this bond (3*(I-1), where I is absolute
   !            atom number.
   !    NTB   : Periodic boundary conditions flag.
   !    IFLAG : =0, then calculate energy and forces.
   !            =1, then calculate only current value of bond.
   !            =2, then calculate energy and d_E/d_distance, but do not
   !                accumulate forces into F() array. They are returned
   !                in the DFR(24) array.
   !            =3, then calculate energy and dE/d_r. Assume the distance
   !                has already been calculated, and is passed in RIJ.
   !   INCFLG : Not currently used.
   
   ! (the following only used when time-averaging):
   ! BAVE(I) : Contains time-averaged values of bonds, using TAUAVE damping.
   ! BAVE0(I): Contains time-averaged values of bonds, using no damping.
   
   !           Both BAVE and BAVE0 are passed so that the first element of
   !           the array corresponds to the average value of the bond currently
   !           of interest, and the second element of the array contains
   !           the previous value of the bond of interest.
   
   !     NAVE : > 0 for time-averaging.
   !   NEXACT : Not currently used.
   !   IPOWER : The power used in averaging the distance data.
   !   TAUAVE : The exponential decay factor used in averaging the distance data
   !   IRAVIN : =0 do not modify the values of the distances by RAVEIN
   !            =1 modify r by RAVEIN, as shown below.
   !   RAVEIN : At time 0, time-average integral is undefined.
   !            If -1000<RAVEIN<1000, use r_initial+ravein on this step.
   !            If RAVEIN<=-1000, use r_target+(ravein+1000) on this step.
   !            If RAVEIN>=1000, use r_target-(ravein-1000) on this step.
   !       DT : The integration time-step (ps)
   !   NAVINT : Time-averaged restraints are only updated every NAVINT
   !            steps. Typically, NAVINT=1.
   
   ! Output:
   !    E     : The energy contribution from this restraint.
   !    DFR   : If IFLAG=2, DFR(1->3) contains the d_E/d_distance forces
   !            for x,y, and z for position 1; those for position 2 are in
   !            DFR(4->6).
   !    RIJ   : The actual bond distance
   
   ! Author: David A. Pearlman
   ! Date: 7/89
   
   implicit none
   integer:: i1, i2, ialtd, iave, iflag, iflg, incflg, ipower, irav, &
        iravin, m, nave, navint, nexact, ntb
   _REAL_ :: a, b, bave, bave0, denom, df, dfr, dif, dif1, dravdr, &
        dt, e, f, r1, r2, r3, r4, ravein, rij, rij2, rinc, rnow, small, &
        tauave, x, xij
   _REAL_ k2,k3
   parameter (small = 1.0d-5)
   dimension x(*),f(*),dfr(*),bave(2),bave0(2)
   dimension xij(3)
   rinc = 1000.0d0
   
   iflg = 1
   if (iflag == 1) iflg = 2
   
   ! Calculate distance:
   
   if (iflag /= 3) then
      rij2 = 0.0d0
      do m=1,3
         xij(m) = x(i1+m) - x(i2+m)
         rij2 = rij2 + xij(m)**2
      end do
      rij = sqrt(rij2)
   end if
   
   ! Get averaged value, if requested
   
   dravdr = 1.0d0
   if (iave > 0) then
      rnow = rij
      
      ! Set the value on the first step (where the integral is not defined):
      
      if (iravin == 1) then
         irav = int(ravein)
         if (-1000 < irav .and. irav < 1000) then
            rnow = rnow + ravein
         else
            rinc = -sign(rinc,ravein)
            if (r2 > small) rnow = ((r2+r3)/2.0d0) + ravein + rinc
            if (r2 <= small) rnow = r3 + ravein + rinc
         end if
      end if
      call aveint(bave,bave0,rnow,tauave,dt,navint,ipower,1, &
            incflg,iflg,rij,denom)
      
      ! DRAVDR is the factor d(r_ave)/d(r_current)
      
      if (iave == 1 .and. iflag /= 1) &
            dravdr = ((rij/rnow)**(1+ipower))*dt/denom
   end if
   
   if (iflag == 1) return
   
   ! Calculate energy (E) and the derivative with respect to the bond
   ! length (DF):
   
   if( ialtd == 1 ) then
      
      if (rij < r2) then
         dif = rij - r2
         df = 4.0d0 * k2 * dif**3
         e = k2*dif**4
      else if (rij <= r3) then
         e = 0.0d0
         do m=1,6
            dfr(m) = 0.0d0
         end do
         return
      else if (rij < r4) then
         dif = rij - r3
         df = 2.0d0 * k3 * dif
         e = k3*dif*dif
      else
         dif1 = r4 - r3
         
         !     --older version, with c = r1; dangerous if input files are not
         !           correctly modified!
         
         !         A = 3.D0*DIF1**2 - 2.D0*R1*DIF1
         !         B = -2.D0*DIF1**3 + R1*DIF1**2
         !         DIF = RIJ - R3
         !         E = K3*(R1*DIF + B/DIF + A)
         !         DF = K3*(R1 - B/DIF**2)
         
         !      -- following lines assume c=0, so potential becomes flat at
         !            large distances:
         
         a = 3.d0*dif1**2
         b = -2.d0*dif1**3
         dif = rij - r3
         e = k3*(b/dif + a)
         df = -k3*(b/dif**2)
         
      end if 
      
   else
      
      if (rij < r1) then
         dif1 = r1-r2
         df = 2.0d0 * k2 * dif1
         e = df * (rij-r1) + k2*dif1*dif1
      else if (rij < r2) then
         dif = rij - r2
         df = 2.0d0 * k2 * dif
         e = k2*dif*dif
      else if (rij <= r3) then
         e = 0.0d0
         do m=1,6
            dfr(m) = 0.0d0
         end do
         return
      else if (rij < r4) then
         dif = rij - r3
         df = 2.0d0 * k3 * dif
         e = k3*dif*dif
      else
         dif1 = r4-r3
         df = 2.0d0 * k3 * dif1
         e = df * (rij-r4) + k3*dif1*dif1
      end if
      
   end if  ! ( ialtd == 1 )
   
   ! Calculate the derivaties with respect to the coordinates, and add them
   ! into the force arrays. If IFLAG=2, do not add them into the force arrays
   ! (probably to be used in center-of-mass restraints). If IFLAG=3, return
   ! with just dE/dr force in DFR(1).
   
   ! DRAVDR contains dr_ave/dr(t) if IAVE=1, 1.0 otherwise.
   
   if (iflag == 3) then
      dfr(1) = df
      return
   end if
   
   df = dravdr*df/rij
   
   do m=1,3
      dfr(m) = xij(m)*df
      dfr(m+3) = -xij(m)*df
   end do
   
   if (iflag == 2) return
   
   do m=1,3
      f(i1+m) = f(i1+m) - dfr(m)
      f(i2+m) = f(i2+m) + dfr(m)
   end do
   
   return
end subroutine disnrg 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Lagrangian multiplier constrain routine
subroutine vec_constr(x,f,dfr,rimass,dt,rcurr,i,j,tgtvec,iflag)

   implicit none

   _REAL_ x(*), f(*), dfr(*), rimass(*), dt, rcurr, tgtvec(3)
   integer i, j, iflag

   _REAL_ dx, dy, dz, ti, tj, lagmul
   _REAL_ dfx, dfy, dfz, dtx
   integer ii

   dtx = dt*20.455d0
   ti = dtx*dtx*rimass(i/3+1)
   tj = dtx*dtx*rimass(j/3+1)
   dx = x(j+1) - x(i+1)
   dy = x(j+2) - x(i+2)
   dz = x(j+3) - x(i+3)
   lagmul = rcurr/(tj+ti)

   dfx = lagmul*(tgtvec(1)-dx/rcurr)
   dfy = lagmul*(tgtvec(2)-dy/rcurr)
   dfz = lagmul*(tgtvec(3)-dz/rcurr)

   if(iflag>0) then
      dfr(1) =  dfx
      dfr(2) =  dfy
      dfr(3) =  dfz
      dfr(4) = -dfx
      dfr(5) = -dfy
      dfr(6) = -dfz
      return
   end if

   f(j+1) = f(j+1) + dfx
   f(j+2) = f(j+2) + dfy
   f(j+3) = f(j+3) + dfz
   f(i+1) = f(i+1) - dfx
   f(i+2) = f(i+2) - dfy
   f(i+3) = f(i+3) - dfz
   return
end subroutine vec_constr

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine echoin here]
subroutine echoin(in,iout)
   
   
   ! Subroutine ECHO INput
   
   ! This routine is called by the MDREAD routine for an Amber suite program.
   ! It reads and echos the information in the input file to the user.
   
   ! Special provision is made for the lines written to the file by the
   ! Amber/Interface to log the Interface lines used to generate the input
   ! file.
   
   ! Before this routine returns, the file on unit IN will be rewound.
   
   ! Author: David A. Pearlman
   ! Date: 8/92
   
   ! IN: The unit the input file is read from.
   ! IOUT: The unit to which the echoed information is to be written.
   implicit none
   integer:: i, in, inter, intlin, iout
   
   character(len=80) aa
   character(len=14) bgflg,enflg
   data bgflg/'!! BEGIN Amber'/
   data enflg/'!! END   Amber'/
   
   ! First echo the Interface script, if any:
   
   inter = 0
   intlin = 0
   do i = 1,999999
      read(in,500,end=20) aa
      if (aa(1:14) == enflg) inter = 0
      if (inter == 1) write(iout,500) aa(1:79)
      if (aa(1:14) == bgflg) then
         inter = 1
         if (intlin == 0) then
            write(iout,*)
            write(iout,1000)
            write(iout,*)
         end if
         intlin = 1
      end if
   end do
   20 continue
   
   ! Now echo the standard Amber input lines:
   
   if (intlin > 0) write(iout,501)
   rewind(in)
   write(iout,*)
   write(iout,1010)
   write(iout,*)
   inter = 0
   do i = 1,999
      read(in,500,end=40) aa
      if (aa(1:14) == bgflg) inter = 1
      if (inter == 0) write(iout,500) aa(1:79)
      if (aa(1:14) == enflg) inter = 0
   end do
   40 continue
   
   rewind(in)
   return
   
   ! Format statements
   
   500 format(a)
   501 format(79('-'))
   1000 format(' The Interface script used to generate the input file:')
   1010 format(' Here is the input file:')
end subroutine echoin 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine disnrg here]
subroutine gendisnrg(x,f,dfr,coord,i1,i2,i3,i4,i5,i6,i7,i8,rstwt,&
      e,r1,r2,r3,r4,k2,k3,ntb, &
      gdave,gdave0,nave,nexact,ipower,tauave,ravein, &
      dt,navint,iave,incflg,iravin,iflag)
   
   
   ! Subroutine GENeralized DIStance reaction coordinate eNeRGy
   
   ! This subroutine calculates the generalized distance coordinate 
   ! restraint energy.
   ! The energy function is a flat-bottomed well, with parabolic sides
   ! which become linear at greater distances.
   
   !           r <= r1  : E = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
   !     r1 <= r <= r2  : E = k2*(r-r2)**2
   !     r2 <= r <= r3  : E = 0
   !     r3 <= r <= r4  : E = k3*(r-r3)**2
   !     r4 <= r        : E = 2*k3*(r4-r3)*(r-r4) + K3*(r4-r3)**2
   
   
   
   ! IAVE = 0: then the instantaneous (current) value of r is used
   !           in all calculations.
   ! IAVE = 1: The time-averaged value of r is used, and forces are calculated
   !           according to the above function, i.e.
   !               dE/dx = (dE/dr_ave)(dr_ave/dr)(dr/dx).
   ! IAVE = 2: The time-averaged value of r is used, and forces are calculated as
   !               dE/dx = (dE/dr_ave)(dr/dx).
   !           This expression, if integrated would not yield the energy function
   !           above, so energies reported are pseudo-energies. Also, it is a
   !           non-conservative force, so the system will tend to heat up. But it
   !           avoids potentially large force terms with (1+IPOWER)
   !           exponentiation, which can cause the simulation to become unstable.
   
   ! Other input:
   !    X(I)  : Coordinate array.
   !    F(I)  : Force array;
   !            Modified on output to include forces from this restraint.
   !    I1,I2 : Atom pointers for this bond (3*(I-1), where I is absolute
   !            atom number.
   !    NTB   : Periodic boundary conditions flag.
   !    IFLAG : =0, then calculate energy and forces.
   !            =1, then calculate only current value of bond.
   !            =2, then calculate energy and d_E/d_distance, but do not
   !                accumulate forces into F() array. They are returned
   !                in the DFR(24) array.
   !            =3, then calculate energy and dE/d_r. Assume the distance
   !                has already been calculated, and is passed in COORD.
   !   INCFLG : Not currently used.
   
   ! (the following only used when time-averaging):
   ! GDAVE(I) : Contains time-averaged values of bonds, using TAUAVE damping.
   ! GDAVE0(I): Contains time-averaged values of bonds, using no damping.
   
   !           Both GDAVE and GDAVE0 are passed so that the first element of
   !           the array corresponds to the average value of the coordinate currently
   !           of interest, and the second element of the array contains
   !           the previous value of the coordinate of interest.
   
   !     NAVE : > 0 for time-averaging.
   !   NEXACT : Not currently used.
   !   IPOWER : The power used in averaging the data.
   !   TAUAVE : The exponential decay factor used in averaging the distance data
   !   IRAVIN : =0 do not modify the values of the distances by RAVEIN
   !            =1 modify r by RAVEIN, as shown below.
   !   RAVEIN : At time 0, time-average integral is undefined.
   !            If -1000<RAVEIN<1000, use r_initial+ravein on this step.
   !            If RAVEIN<=-1000, use r_target+(ravein+1000) on this step.
   !            If RAVEIN>=1000, use r_target-(ravein-1000) on this step.
   !       DT : The integration time-step (ps)
   !   NAVINT : Time-averaged restraints are only updated every NAVINT
   !            steps. Typically, NAVINT=1.
   
   ! Output:
   !    E     : The energy contribution from this restraint.
   !    DFR   : If IFLAG=2, DFR(1->3) contains the d_E/d_distance forces
   !            for x,y, and z for position 1; those for position 2 are in
   !            DFR(4->6).
   !    COORD : The actual coordinate value
   
   ! Author:  Matthew Seetin 
   ! Date: 10/2007
   
   implicit none
   integer:: i1, i2, i3, i4, i5, i6, i7, i8, ialtd, iave, iflag, &
        iflg, incflg, ipower, irav, iravin, m, nave, navint, nexact, &
        ntb
   _REAL_ :: a, b, coord, denom, df, dfr, dif, dif1, dravdr, dt, e, &
        f, gdave, gdave0, r1, r2, r3, r4, ravein, rinc, rnow, rstwt, &
        small, tauave, x
   _REAL_ k2,k3,xij,rij,rij2
   INTEGER iat,used,i
   parameter (small = 1.0d-7)
   dimension x(*),f(*),dfr(*),rstwt(4),gdave(2),gdave0(2)
   dimension xij(3,4),iat(8),rij(4)
   rinc = 1000.0d0
   
   used = 0
   coord = 0.0d0
   
   iflg = 1
   if (iflag == 1) iflg = 2
   
   do i=1,4
     if ( abs(rstwt(i)) > small) used = used + 1
   end do
   
   iat(1) = i1
   iat(2) = i2
   iat(3) = i3
   iat(4) = i4
   iat(5) = i5
   iat(6) = i6
   iat(7) = i7
   iat(8) = i8
   
   do m=1,24
     dfr(m) = 0.0d0
   end do
   
   ! Calculate distances:
   
   if (iflag /= 3) then
     do i=1,used
       rij2 = 0.0d0
       do m=1,3
          xij(m,i) = x(iat(2*i-1)+m) - x(iat(2*i)+m)
          rij2 = rij2 + xij(m,i)**2
       end do
       rij(i) = sqrt(rij2)
       coord = coord + rstwt(i) * rij(i)
     end do
   end if
   
   ! Get averaged value, if requested
   
   dravdr = 1.0d0
   if (iave > 0) then
      rnow = coord
      
      ! Set the value on the first step (where the integral is not defined):
      
      if (iravin == 1) then
         irav = int(ravein)
         if (-1000 < irav .and. irav < 1000) then
            rnow = rnow + ravein
         else
            rinc = -sign(rinc,ravein)
            if (r2 > small) rnow = ((r2+r3)/2.0d0) + ravein + rinc
            if (r2 <= small) rnow = r3 + ravein + rinc
         end if
      end if
      call aveint(gdave,gdave0,rnow,tauave,dt,navint,ipower,7, &
            incflg,iflg,coord,denom)
      
      ! DRAVDR is the factor d(r_ave)/d(r_current)
      
      if (iave == 1 .and. iflag /= 1) &
            dravdr = ((coord/rnow)**(1+ipower))*dt/denom
   end if
   
   if (iflag == 1) return

   ! Calculate energy (E) and the derivative with respect to the bond
   ! length (DF):
   
   if( ialtd == 1 ) then
      if (coord < r2) then
         dif = coord - r2
         df = 4.0d0 * k2 * dif**3
         e = k2*dif**4
      else if (coord <= r3) then
         e = 0.0d0
         return
      else if (coord < r4) then
         dif = coord - r3
         df = 2.0d0 * k3 * dif
         e = k3*dif*dif
      else
         dif1 = r4 - r3
         
         !     --older version, with c = r1; dangerous if input files are not
         !           correctly modified!
         
         !         A = 3.D0*DIF1**2 - 2.D0*R1*DIF1
         !         B = -2.D0*DIF1**3 + R1*DIF1**2
         !         DIF = RIJ - R3
         !         E = K3*(R1*DIF + B/DIF + A)
         !         DF = K3*(R1 - B/DIF**2)
         
         !      -- following lines assume c=0, so potential becomes flat at
         !            large distances:
         
         a = 3.d0*dif1**2
         b = -2.d0*dif1**3
         dif = coord - r3
         e = k3*(b/dif + a)
         df = -k3*(b/dif**2)
         
      end if 
      
   else
      if (coord < r1) then
         dif1 = r1-r2
         df = 2.0d0 * k2 * dif1
         e = df * (coord-r1) + k2*dif1*dif1
      else if (coord < r2) then
         dif = coord - r2
         df = 2.0d0 * k2 * dif
         e = k2*dif*dif
      else if (coord <= r3) then
         e = 0.0d0
         return
      else if (coord < r4) then
         dif = coord - r3
         df = 2.0d0 * k3 * dif
         e = k3*dif*dif
      else
         dif1 = r4-r3
         df = 2.0d0 * k3 * dif1
         e = df * (coord-r4) + k3*dif1*dif1
      end if
      
   end if  ! ( ialtd == 1 )
   
   ! Calculate the derivaties with respect to the coordinates, and add them
   ! into the force arrays. If IFLAG=2, do not add them into the force arrays
   ! (probably to be used in center-of-mass restraints). If IFLAG=3, return
   ! with just dE/dr force in DFR(1).
   
   ! DRAVDR contains dr_ave/dr(t) if IAVE=1, 1.0 otherwise.
   
   if (iflag == 3) then
      dfr(1) = df
      return
   end if
   
   df = dravdr*df
   
   do i=1,used
     do m=1,3
       dfr(6*i-6+m) = df * rstwt(i) * xij(m,i) / rij(i)
       dfr(6*i-3+m) = -dfr(6*i-6+m)
     end do
   end do
   
   if (iflag == 2) return
   
   do i=1,2*used
     do m=1,3
       f(iat(i)+m) = f(iat(i)+m) - dfr(3*(i-1) + m)
     end do
   end do
   
   return
end subroutine gendisnrg 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine getnat here]
subroutine getnat(iatlocal,atnm,name,ipres,nres,iout,ierr)
   
   
   ! Subroutine GET Numbers of AToms.
   
   ! This routine takes a pointers to a residue numbers (IAT), and
   ! the name of the desired atom in this residue (ATNM), and
   ! returns the actual atom number pointers in IAT.
   
   ! Author: David A. Pearlman
   ! Date: 3/89
   
   implicit none
   integer:: i, iatlocal, ierr, iout, ipres, nres
   character(len=4) atnm, name
   dimension name(*),ipres(*)
   
   if (iatlocal <= 0) return
   if (iatlocal > nres) goto 8000
   
   do i = ipres(iatlocal),ipres(iatlocal+1)-1
      if (atnm == name(i)) then
         iatlocal = i
         return
      end if
   end do
   
   ! If we get here, the desired atom was not found. Report this to the user
   ! and return with IERR = 1
   
   write(iout,9000) atnm,iatlocal
   ierr = 1
   return
   
   8000 write(iout,9001) atnm,iatlocal
   ierr = 1
   return
   
   ! Error format statements
   
   9000 format(' Error: No atom ',a4,' in residue ',i5)
   9001 format(' Error: residue_number/atom_name reference specifies a', &
         /,'        residue number greater than the last residue in', &
         ' the system.',/,t9,'Res = ',i6,', Atom = "',a4,'"')
   
end subroutine getnat 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine impnum here]
subroutine impnum(lph,lpa,icph,icpa,nphih,nphia,nfpert,numphi, &
      nimprp)
   
   ! Subroutine IMProper Number
   
   ! This simple subroutine looks through the LP() array for improper torsions.
   ! The number of improper torsional parameters (not torsions themselves)
   ! is returned. This value is used in modifying the relative
   ! weight of regular torsion terms vs. improper torsion terms.
   
   ! Author: David A. Pearlman
   ! Date: 8/89
   
   ! INPUT
   
   ! LPH(I) : Array containing the fourth atom of torsion I. For torsions with
   !          hydrogens at one/both ends. If LP(I) < 0, then this is an
   !          improper torsion.
   ! LPA(I) : Array containing fourth atom of torsion I. For torsions with
   !          no hydrogens at ends. For perturbation calculations (NFPERT>0),
   !          these are followed by 2*NFPERT perturbed torsions.
   ! ICPH(I): ICP(I) is a pointer into the force constant array for torsion I.
   !          For torsions with hydrogens at one/both ends.
   ! ICPA(I): Force constant pointers corresponding to the torsions of LPA(I).
   ! NPHIH:   Number of torsion angles with hydrogens at one/both ends.
   ! NPHIA:   Number of torsion angles without hydrogens at the ends.
   ! NFPERT:  If NFPERT > 0, there are NFPERT perturbed torsions (GIBBS only).
   ! NUMPHI: Total number of torsional parameters (regular+improper)
   
   ! OUTPUT
   
   ! NIMPRP: The number of improper torsional parameters. The lowest value ICP(I)
   !         corresponding to an improper is used to calculate this value, assuming
   !         improper torsion parameters are packed after the regular torsions.
   
   implicit none
   integer:: i, icpa, icph, lpa, lph, nfpert, nimprp, nphia, nphih, numphi
   dimension lph(*),lpa(*),icph(*),icpa(*)
   
   nimprp = numphi + 1
   
   do i=1,nphih
      if (lph(i) < 0) then
         if (icph(i) < nimprp) nimprp = icph(i)
      end if
   end do
   
   do i=1,nphia + 2*nfpert
      if (lpa(i) < 0) then
         if (icpa(i) < nimprp) nimprp = icpa(i)
      end if
   end do
   
   nimprp = numphi - nimprp + 1
   
   return
end subroutine impnum 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine jnrg here]
subroutine jnrg  (x         ,f         ,ajval     ,i1        , &
      i2        ,i3        ,i4        ,e         ,r1        , &
      r2        ,r3        ,r4        ,k2        ,k3        , &
      ntb       ,tave1     ,tave01    ,tave2     ,tave02    , &
      nave      ,nexact    ,ipower    ,tauave    ,ravein    , &
      dt        ,navint    ,iave      ,incflg    ,iravin    , &
      ajcoef    ,iflag)
   
   
   ! Subroutine J coupling eNeRGy
   
   ! This routine calculates the restraint energy due to a J-Coupling
   ! restraint.
   
   ! J is calculated using a Karplus equation expressed in terms of
   ! of torsion angle cosines:
   
   !       J(tau) = A cos(tau)^2 + B cos(tau) + C
   
   ! The energy function is a flat-bottomed well, with parabolic sides
   ! which become linear at greater distances.
   
   ! Author: David A. Pearlman
   ! Date: 3/93
   
   ! If r is the torsion angle:
   
   !       r <= r1  : E = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
   ! r1 <= r <= r2  : E = k2*(r-r2)**2
   ! r2 <= r <= r3  : E = 0
   ! r3 <= r <= r4  : E = k3*(r-r3)**2
   ! r4 <= r        : E = 2*k3*(r4-r3)*(r-r4) + K3*(r4-r3)**2
   
#ifdef LES
   !                  if k2<0, then this torsion is assumed to be one of a
   !                    series to be averaged.  Use the absolute value of k2
   !                    as the constraint weight, and accumulate derivatives
   !                    in a local variable; when a positive k2 is found,
   !                    update global derivatives
#endif
   
   ! IAVE = 0: then the instantaneous (current) value of r is used
   !           in all calculations.
   ! IAVE = 1: The time-averaged value of r is used, and forces are calculated
   !           according to the above function, i.e.
   !               dE/dx = (dE/dr_ave)(dr_ave/dr)(dr/dx).
   ! IAVE = 2: The time-averaged value of r is used, and forces are calculated as
   !               dE/dx = (dE/dr_ave)(dr/dx).
   !           This expression, if integrated would not yield the energy function
   !           above, so energies reported are pseudo-energies. Also, it is a
   !           non-conservative force, so the system will tend to heat up. But it
   !           avoids potentially large force terms with (1+IPOWER)
   !           exponentiation, which can cause the simulation to become unstable.
   
   ! Other input:
   !    X(I)  : Coordinate array.
   !    F(I)  : Force array;
   !            Modified on output to include forces from this restraint.
   !  I1,I2,I3: Atom pointers for this angle (3*(I-1), where I is absolute
   !        I4  atom number.
   !    NTB   : Periodic boundary conditions flag.
   !    IFLAG : =0, then calculate energy and forces.
   !            =1, then only calculate value of torsion.
   !            =2, then calculate value and energy, but do not accumulate derivatives
   !   INCFLG : Determines whether local saved pointers are updated in AVEINT
   
   ! (following only used when time-averaging:)
   ! TAVE1(I)   Contain the time-averaged values of J(tau) using damped
   ! TAVE01(I): exponential weighting and no weighting, respectively.
   
   ! TAVE2(I):  Not currently used in this routine.
   ! TAVE02(I):  "     "       "    "   "    "
   
   !           Both TAVE and TAVE0 are passed so that the first element of
   !           the array corresponds to the average value of the J(tau) currently
   !           of interest, and the second element of the array contains
   !           the previous value of the J(tau) of interest.
   
   !     NAVE : > 0 for time-averaged restraints.
   !   NEXACT : Not currently used.
   !   IPOWER : The power used in averaging the J(tau) data.
   !   TAUAVE : The exponential decay factor used in averaging the J(tau) data
   !   IRAVIN : =0 do not modify the values of the initial J(tau) by RAVEIN
   !            =1 modify J(tau) by RAVEIN, as shown below.
   !   RAVEIN : At time 0, time-average integral is undefined.
   !            If -1000<RAVEIN<1000, use r_initial+ravein on this step.
   !            If RAVEIN<=-1000, use r_target+(ravein+1000) on this step.
   !            If RAVEIN>=1000, use r_target-(ravein-1000) on this step.
   !       DT : Integration time-step (ps).
   !   NAVINT : Only store restraints in TAVE/TAVE0 every NAVINT steps.
   
   ! AJCOEF(3): The 3 coefficients in the Karplus Equation for this torsion.
   
   ! Output:
   !    E     : The energy contribution from this restraint.
   !    AJVAL : The calculated J(tau)
   
   use constants, only : PI
   implicit none
   integer:: i1, i2, i3, i4, iave, iflag, iflg, incflg, ipower, &
        irav, iravin, m, nave, navint, nexact, nj, ntb
   _REAL_ :: ajave, ajcoef, ajnow, ajval, amlt, ap, bi, bk, cphi, &
        ct, dc, denom, df, dif, dif1, dr1, dr2, dr3, dr4, dr5, dr6, &
        dravdr, drx, dry, drz, dt, dx, dy, dz, e, f, gx, gy, gz, r1, r2, &
        r3, r4, ravein, rij2, rinc, rkj2, rkl2, s, small, small2, sphi, &
        t, tauave, tave01, tave02, tave1, tave2, x, xij, xkj, xkl, z1, &
        z11, z12, z2, z22
   _REAL_ k2,k3
   parameter (small=1.0d-14)
   parameter (small2=1.0d-5)
   dimension x(*),f(*),tave1(2),tave01(2),tave2(2),tave02(2)
   dimension xij(3),xkj(3),xkl(3),t(6),dc(6)
   dimension ajcoef(3)
#ifdef LES
   integer :: i
   integer, parameter :: maxatom3=15000
   _REAL_ drv(maxatom3)
   data drv /maxatom3 * 0.0/
   data ajave /0.0d0/
   data nj / 0 /
#  include "nmr.h"
#  include "memory.h"
#endif
   
   rinc = 1000.0d0
   
   iflg = 1
   if (iflag == 1) iflg = 2
   
   ! Calculate the underlying torsion:
   
   rij2 = 0.0d0
   rkj2 = 0.0d0
   rkl2 = 0.0d0
   
   do m=1,3
      xij(m) = x(i1+m) - x(i2+m)
      xkj(m) = x(i3+m) - x(i2+m)
      xkl(m) = x(i3+m) - x(i4+m)
      rij2 = rij2 + xij(m)**2
      rkj2 = rkj2 + xkj(m)**2
      rkl2 = rkl2 + xkl(m)**2
   end do
   
   ! Calculate ij X jk AND kl X jk
   
   dx = xij(2)*xkj(3) - xij(3)*xkj(2)
   dy = xij(3)*xkj(1) - xij(1)*xkj(3)
   dz = xij(1)*xkj(2) - xij(2)*xkj(1)
   
   gx = xkj(3)*xkl(2) - xkj(2)*xkl(3)
   gy = xkj(1)*xkl(3) - xkj(3)*xkl(1)
   gz = xkj(2)*xkl(1) - xkj(1)*xkl(2)
   
   
   ! Calculate the magnitudes of above vectors, and their dot product:
   
   bi = dx*dx + dy*dy + dz*dz
   bk = gx*gx + gy*gy + gz*gz
   ct = dx*gx + dy*gy + dz*gz
   
   ! If this is a linear dihedral, we cannot calculate a value for it,
   ! so set energy to zero and return:
   
   if (bk < 0.01d0 .or. bi < 0.01d0) then
      e = 0.0d0
      return
   end if
   
   bi = sqrt(bi)
   bk = sqrt(bk)
   z1 = 1.0d0/bi
   z2 = 1.0d0/bk
   ct = ct*z1*z2
   
   ! Cannot let ct be exactly 1. or -1., because this results
   ! in sin(tau) = 0., and we divide by sin(tau) [sphi] in calculating
   ! the derivatives.
   
   if (ct > 1.0d0-small) ct = 1.0d0-small
   if (ct < -1.0d0+small) ct = -1.0d0+small
   ap = acos(ct)
   
   s = xkj(1)*(dz*gy-dy*gz) + xkj(2)*(dx*gz-dz*gx) + &
         xkj(3)*(dy*gx-dx*gy)
   
   if (s < 0.0d0) ap = -ap
   ap = pi - ap
   
   cphi = cos(ap)
   sphi = sin(ap)
   ct = -ct
   
   ! Define the current value of the karplus-equation-calculated J value.
   ! If time-averaging is being performed, AJVAL will be replaced by its
   ! time-averaged value in routine AVEINT:
   
   ajnow = ajcoef(1) * cphi**2 + ajcoef(2)*cphi + ajcoef(3)
#ifdef LES
   ajave = ajave + ajnow
   nj = nj + 1
#endif
   ajval = ajnow
   
   ! Get averaged value, if requested
   
   dravdr = 1.0d0
   if (iave > 0) then
      
      ! Set the value on the first step (where integral not defined):
      
      if (iravin == 1) then
         irav = int(ravein)
         if (-1000 < irav .and. irav < 1000) then
            ajnow = ajnow + ravein
         else
            rinc = -sign(rinc,ravein)
            if (r2 > small2) ajnow = ((r2+r3)/2.0d0) + ravein+rinc
            if (r2 <= small2) ajnow = r3 + ravein + rinc
         end if
      end if
      
      call aveint(tave1,tave01,ajnow,tauave,dt,navint,ipower,3, &
            incflg,iflg,ajval,denom)
      
      ! DRAVDR is the factor d(<J>)/d(J(t))
      
      if (iave == 1 .and. iflag /= 1) &
            dravdr = ((ajval/ajnow)**(1+ipower))*dt/denom
   end if
   
   if (iflag == 1) then
      ajave = 0.0d0
      nj = 0
      return
   end if
#ifndef LES
   
   ! Calculate energy (E) and the derivative with respect to the torsion (DF):
   
   if (ajval < r1) then
      dif1 = r1-r2
      df = 2.0d0 * k2 * dif1
      e = df * (ajval-r1) + k2*dif1*dif1
   else if (ajval < r2) then
      dif = ajval - r2
      df = 2.0d0 * k2 * dif
      e = k2*dif*dif
   else if (ajval <= r3) then
      e = 0.0d0
      return
   else if (ajval < r4) then
      dif = ajval - r3
      df = 2.0d0 * k3 * dif
      e = k3*dif*dif
   else
      dif1 = r4-r3
      df = 2.0d0 * k3 * dif1
      e = df * (ajval-r4) + k3*dif1*dif1
   end if
   if (iflag == 2) return
   
   ! Calculate the total forces by the chain rule:
   
   !      dE/dx = dE/d<J> * d<J>/dJ(t) * dJ(t)/dcos(tau) * dcos(tau)/dx
   
   ! Multiply DF, which is dE/d<J> by DRAVDR, which is
   ! d<J>/d[J(t)] if IAVE=1, and 1.0 otherwise:
   
   df = dravdr*df
   
   ! Multiply DF by d[J(t)]/d[cos(tau)]
   
   amlt = 2.0d0*ajcoef(1)*cphi + ajcoef(2)
   df = df * amlt
   
   ! Now calculate derivatives; we have already calculated DF = DE/D{cos(tau)}.
   ! now we need D{cos(tau)}/D{X,Y,OR Z}
   
   t(1) = dx
   t(2) = dy
   t(3) = dz
   t(4) = -gx
   t(5) = -gy
   t(6) = -gz
   
   !  DC = First derivative of cos(phi) w/respect
   !  to the cartesian differences t
   
   z11 = z1*z1
   z12 = z1*z2
   z22 = z2*z2
   
   do m = 1,3
      dc(m) = t(m+3)*z12-cphi*t(m)*z11
      dc(m+3) = t(m)*z12-cphi*t(m+3)*z22
   end do
   
   dr1 = df*(dc(3)*xkj(2) - dc(2)*xkj(3))
   dr2 = df*(dc(1)*xkj(3) - dc(3)*xkj(1))
   dr3 = df*(dc(2)*xkj(1) - dc(1)*xkj(2))
   dr4 = df*(dc(6)*xkj(2) - dc(5)*xkj(3))
   dr5 = df*(dc(4)*xkj(3) - dc(6)*xkj(1))
   dr6 = df*(dc(5)*xkj(1) - dc(4)*xkj(2))
   drx = df*(-dc(2)*xij(3) + dc(3)*xij(2) + dc(5)*xkl(3) - &
         dc(6)*xkl(2))
   dry = df*( dc(1)*xij(3) - dc(3)*xij(1) - dc(4)*xkl(3) + &
         dc(6)*xkl(1))
   drz = df*(-dc(1)*xij(2) + dc(2)*xij(1) + dc(4)*xkl(2) - &
         dc(5)*xkl(1))
   
   f(i1+1) = f(i1+1) - dr1
   f(i1+2) = f(i1+2) - dr2
   f(i1+3) = f(i1+3) - dr3
   f(i2+1) = f(i2+1) - drx + dr1
   f(i2+2) = f(i2+2) - dry + dr2
   f(i2+3) = f(i2+3) - drz + dr3
   f(i3+1) = f(i3+1) + drx + dr4
   f(i3+2) = f(i3+2) + dry + dr5
   f(i3+3) = f(i3+3) + drz + dr6
   f(i4+1) = f(i4+1) - dr4
   f(i4+2) = f(i4+2) - dr5
   f(i4+3) = f(i4+3) - dr6
   
#else
   
   !           get D{cos(tau)}/D{X,Y,OR Z}
   
   t(1) = dx
   t(2) = dy
   t(3) = dz
   t(4) = -gx
   t(5) = -gy
   t(6) = -gz
   
   !              DC = First derivative of cos(phi) w/respect
   !                   to the cartesian differences t
   
   z11 = z1*z1
   z12 = z1*z2
   z22 = z2*z2
   
   do m = 1,3
      dc(m) = t(m+3)*z12-cphi*t(m)*z11
      dc(m+3) = t(m)*z12-cphi*t(m+3)*z22
   end do
   
   
   ! Multiply by AMLT = d[J(t)]/d[cos(tau)]
   
   amlt = 2.0d0*ajcoef(1)*cphi + ajcoef(2)
   
   dr1 = amlt*(dc(3)*xkj(2) - dc(2)*xkj(3))
   dr2 = amlt*(dc(1)*xkj(3) - dc(3)*xkj(1))
   dr3 = amlt*(dc(2)*xkj(1) - dc(1)*xkj(2))
   dr4 = amlt*(dc(6)*xkj(2) - dc(5)*xkj(3))
   dr5 = amlt*(dc(4)*xkj(3) - dc(6)*xkj(1))
   dr6 = amlt*(dc(5)*xkj(1) - dc(4)*xkj(2))
   drx = amlt*(-dc(2)*xij(3) + dc(3)*xij(2) + dc(5)*xkl(3) - &
         dc(6)*xkl(2))
   dry = amlt*( dc(1)*xij(3) - dc(3)*xij(1) - dc(4)*xkl(3) + &
         dc(6)*xkl(1))
   drz = amlt*(-dc(1)*xij(2) + dc(2)*xij(1) + dc(4)*xkl(2) - &
         dc(5)*xkl(1))
   
   drv(i1+1) = drv(i1+1) - dr1
   drv(i1+2) = drv(i1+2) - dr2
   drv(i1+3) = drv(i1+3) - dr3
   drv(i2+1) = drv(i2+1) - drx + dr1
   drv(i2+2) = drv(i2+2) - dry + dr2
   drv(i2+3) = drv(i2+3) - drz + dr3
   drv(i3+1) = drv(i3+1) + drx + dr4
   drv(i3+2) = drv(i3+2) + dry + dr5
   drv(i3+3) = drv(i3+3) + drz + dr6
   drv(i4+1) = drv(i4+1) - dr4
   drv(i4+2) = drv(i4+2) - dr5
   drv(i4+3) = drv(i4+3) - dr6
   
   if( iflag == 2 ) write( iuse, '(35x,f10.3,"(",f9.3,")")') &
         ajnow,57.296*ap
   if( k2 < 0.0 ) then
      e = 0.0
      return
   end if
   
   ! Calculate energy (E) and the derivative with respect to the torsion (DF):
   
   ajval = ajave/nj
   
   if (ajval < r1) then
      dif1 = r1-r2
      df = 2.0d0 * abs(k2) * dif1
      e = df * (ajval-r1) + abs(k2)*dif1*dif1
   else if (ajval < r2) then
      dif = ajval - r2
      df = 2.0d0 * abs(k2) * dif
      e = abs(k2)*dif*dif
   else if (ajval <= r3) then
      e = 0.0d0
      return
   else if (ajval < r4) then
      dif = ajval - r3
      df = 2.0d0 * k3 * dif
      e = k3*dif*dif
   else
      dif1 = r4-r3
      df = 2.0d0 * k3 * dif1
      e = df * (ajval-r4) + k3*dif1*dif1
   end if
   !     write(6,*) 'ajval,e: ',ajval,e
   if (iflag == 2) then
      ajave = 0.0
      nj = 0
      return
   end if
   
   ! Calculate the total forces by the chain rule:
   
   !      dE/dx = dE/d<J> * d<J>/dJ(t) * dJ(t)/dcos(tau) * dcos(tau)/dx
   
   ! Multiply DF, which is dE/d<J> by DRAVDR, which is
   ! d<J>/d[J(t)] if IAVE=1, and 1.0 otherwise:
   
   df = dravdr*df/nj
   
   ! Now update derivatives; we have already calculated DF = DE/D{cos(tau)}.
   !    D{cos(tau)}/D{X,Y,OR Z} is stored in drv()
   
   do i=1,3*natom
      f(i) = f(i) + df*drv(i)
      drv(i) = 0.0
   end do
   ajave = 0.0
   nj = 0
#endif
   return
end subroutine jnrg 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ MODify WeighTs.
subroutine modwt(wtnrg,iwtstp,iwttyp,ichgwt,ishrtb,nstep,temp0, &
      tautp,cut,rk,tk,pk,cn1,cn2,ag,bg,cg,rk2nmr, &
      rk3nmr,nmrshb,nmrst,numbnd,numang,numphi,nimprp, &
      nhb,ncharg,ntypes,rad,wel,radhb,welhb, &
      rwell,tauave,wtls,isftrp,tgtrmsd,nmrnum,temp0les)
   
   ! Subroutine MODify WeighTs
   
   ! This routine modifies the various energy coefficients as requested.
   
   ! Author: David A. Pearlman
   ! Date: 7/89
   
   ! INPUT
   ! -----
   !  WTNRG(2,I): Slope/intercept describing how the relative weights change
   !              with step number.
   ! IWTSTP(3,I): The beginning and ending steps over which a modification
   !              is to be active. If IWTSTP(3,I)>0, it gives the step increment
   !              between changes in the target value.
   !   IWTTYP(I): The type of modification I.
   !      ICHGWT: Number of relative weight change instructions.
   !      ISHRTB: If ISHRTB >0, short-range interactions have been defined;
   !              Relevant information is then stored in IWTSTP( ,ISHRTB) and
   !              WTNRG( ,ISHRTB).
   !       NSTEP: Current iteration/step number
   !       TEMP0: Target temperature (molecular dynamics only)
   !    TEMP0LES: Target temperature for LES region (molecular dynamics only)
   !       TAUTP: Berendsen temperature scaling parameter (md only)
   !       CUT  : (Non-bonded cutoff distance)**2
   !       RK(I): Force constants for bonds
   !       TK(I): Force constants for angles
   !       PK(I): Force constants for torsions. Regular torsions, followed by
   !              impropers.
   !      CN1(I): "A" (r**12) coefficient array for vdw term
   !      CN2(I): "B" (r**6) coefficient array for vdw term
   !       AG(I): "A" (r**12) coefficient for h-bond term
   !       BG(I): "B" (r**10) coefficient for h-bond term
   !       CG(I): Charge on atom I.
   ! RK2NMR(2,I): Slope/intercept of line describing how k2 force constants for NMR
   !              restraints vary with step number.
   ! RK3NMR(2,I): Slope/intercept of line describing how k3 force constants for NMR
   !              restraints vary with step number.
   !   NMRSHB(I): If NMRSHB(I) >0, then restraint I has been defined as a
   !              "short-range" restraint.
   !  NMRST(3,I): NMRST(3,I) > 0 if force constants for NMR restraint I are linearly
   !              interpolated. NMRST(3,I) < 0 if force constants for NMR restraint
   !              I are generated by constant factor multiplication. If the
   !              latter is true, do not modify RK?NMR(1,I) when restraint weights
   !              are modified.
   
   !      NUMBND: Number of bond parmameters
   !      NUMANG: Number of angle parameters
   !      NUMPHI: Number of torsional parameters
   !      NIMPRP: Number of "improper" torsional parameters (not included in NUMPHI
   !       NTTYP: Number of 6-12 (vdw) non-bond parameters
   !        NPHB: Number of hydrogen bond-parameters
   !      NCHARG: Number of charge parameters. Equal to the number of atoms
   !              for minimization and regular molecular dynamics. Equal to
   !              2*Number of atoms for perturbed systems (GIBBS).
   !      NTYPES: Number of non-bonded 6-12 parameter types.
   !      RAD(I): The ideal vdw radius (r_star) for atom-type i.
   !      WEL(I): The well depth (epsilon) for atom-type i.
   !    RADHB(I): The mixed interaction distance for h-bond interaction I.
   !    WELHB(I): The mixed interaction epsilon for h-bond interaction I.
   !       RWELL: The force constant used in the "soft-repulsion" function
   !              (only if ISFTRP > 0; see below).
   !   TAUAVE(3): The exponential decay constants used if time-averaged restraints
   !              are implemented
   !      ISFTRP: If > 0, then 6-12 potential replaced by a "soft-repulsion"
   !              function: E = K(r_star**2-r**2)**2; E=0 for r> r_star.
   !     TGTRMSD: current target RMSD value for targeted MD
   !      NMRNUM: Number of added NMR restraints.
   
   ! Output:
   !     WTLS(2): WTLS(1) is set = WEIGHT(10)  (used in resetting short-range)
   !              WTLS(2) is set = WEIGHT(11)  (used in resetting short-range)
   
   ! Local Common:
   ! Array WEIGHT--
   !  WEIGHT(1) : The current weight of the bonds
   !  WEIGHT(2) : The current weight of the angles
   !  WEIGHT(3) : The current weight of the torsions
   !  WEIGHT(4) : The current weight of the attractive 6 (vdw) term
   !  WEIGHT(5) : The current weight of the replusive 12 (vdw) term
   !  WEIGHT(6) : The current weight of the attractive 10 h-bond term
   !  WEIGHT(7) : The current weight of the replusive 12 h-bond term
   !  WEIGHT(8) : The current weight of the electrostatic terms
   !  WEIGHT(9) : The current weight of the "improper" torsions
   !  WEIGHT(10) : The current weight of the "short-range" NMR restraints
   !  WEIGHT(11) : The current weight of the non-"short-range" NMR restraints
   !  WEIGHT(12) : The current weight of the NOESY restraints
   !  WEIGHT(13) : The current weight of the chemical shift restraints
   !  WEIGHT(14) : Temperature scaling parameter on the first call
   !  WEIGHT(15): The "soft-repulsion" force constant as specified in the input
   !  WEIGHT(16) : The value of the target temperature on the first call
   !  WEIGHT(17) : The square of the non-bonded cutoff on the first call
   !  WEIGHT(18) : Exponential decay time constant for time-averaged bond restraints
   !  WEIGHT(19) : Exponential decay time constant for time-averaged angle restraints
   !  WEIGHT(20) : Exponential decay time constant for time-averaged torsion restraints
   !  WEIGHT(21) : The current target RMSD value for targeted MD (itgtmd=1)
   !  WEIGHT(22) : Exponential decay time constant for time-averaged plane-point restraints
   !  WEIGHT(23) : Exponential decay time constant for time-averaged plane-plane restraints
   !  WEIGHT(24) : The value of the LES target temperature on the first call
   use constants, only : one
   use parms, only : nttyp
#ifdef MPI
   ! initmodwt forces modwt() to re-read values such as temp0 - 
   ! otherwise they will be reset to the initial value after each 
   ! exchange.
   use remd, only : initmodwt
#endif
   implicit none
   
   integer jweit
   parameter (jweit=25)
#  include "nmr.h"

   _REAL_  ag
   _REAL_  bg
   _REAL_  cg
   _REAL_  cn1
   _REAL_  cn2
   _REAL_  cut
   integer ichgwt
   integer isftrp
   integer ishrtb
   integer iwtstp
   integer iwttyp
   integer ncharg
   integer nhb
   integer nimprp
   integer nmrnum
   integer nmrshb
   integer nmrst
   integer nstep
   integer ntypes
   integer numang
   integer numbnd
   integer numphi
   _REAL_  pk
   _REAL_  rad
   _REAL_  radhb
   _REAL_  rk
   _REAL_  rk2nmr
   _REAL_  rk3nmr
   _REAL_  rwell
   _REAL_  tauave
   _REAL_  tautp
   _REAL_  temp0
   _REAL_  temp0les
   _REAL_  tgtrmsd
   _REAL_  tk
   _REAL_  wel
   _REAL_  welhb
   _REAL_  wtls
   _REAL_  wtnrg
   dimension wtnrg(2,*),iwtstp(3,*),iwttyp(*)
   dimension rk(*),tk(*),pk(*),cn1(*),cn2(*),rad(*),wel(*),wtls(2)
   dimension ag(*),bg(*),cg(*),rk2nmr(2,*),rk3nmr(2,*),nmrshb(*)
   dimension nmrst(3,*),tauave(6),radhb(*),welhb(*)
   dimension ichang(jweit),ichold(jweit),ireset(jweit),weight(jweit)
   
   logical fixed
   integer i
   integer ichang
   integer ichold
   integer ihbp
   integer indx
   integer ionce
   integer ireset
   integer irstr
   integer irstyp
   integer itorfc
   integer itype
   integer ivdw
   integer j
   integer nstepu
   _REAL_  rmult
   _REAL_  weight,small,wt,wt2,wt3
   parameter (small = 1.0d-12)
   save weight,ichold
   data ionce/0/

   ! On the first call, initialize the ICHOLD indicator. On subsequent calls,
   ! ICHOLD(I)=0, if weight of term was set/remained at 1.0 on previous call.
   ! ICHOLD(I)=1, if weight of term was set to non-unity as of the previous call.
   ! ICHOLD(I)=2, if weight of term has been set to 0, or if the weight
   !              of the term is constant throughout the run (NSTEP1=NSTEP2=0).

#ifdef MPI
   if (ionce /= 99.or.initmodwt) then
#else   
   if (ionce /= 99) then
#endif

      ichold(1:jweit) = 0
      weight(1:13) = ONE
      weight(14) = tautp
      weight(15) = rwell
      weight(16) = temp0
      weight(17) = cut
      weight(18) = tauave(1)
      weight(19) = tauave(2)
      weight(20) = tauave(3)
      weight(21) = tgtrmsd
      weight(22) = tauave(4)
      weight(23) = tauave(5)
      weight(24) = tauave(6)
      weight(25) = temp0les
      ionce = 99
#ifdef MPI
      initmodwt = .false.
#endif

   end if
   ivdw = 0
   ihbp = 0
   itorfc = 0
   
   ! Initialize the change indicator. Any energy term not modified by one
   ! of the change instructions should be set to the default weight of 1.0
   ! at the end of this routine.
   
   ichang(1:jweit) = 0
   ireset(1:jweit) = 0
   
   ! Loop over the requested changes. Use an implied loop, rather than
   ! an explicit one, so that this code can be conveniently used to
   ! reset unchanged parameter weights to 1.0 later on.
   
   if (ichgwt <= 0) goto 75
   irstyp = 0
   i = 1
   20 continue
   fixed = (iwtstp(2,i) == 0)
   
   ! Skip, if the current step is outside of the desired range
   
   if (nstep < iwtstp(1,i) .or. (nstep > iwtstp(2,i).and. &
         iwtstp(2,i) /= 0)) goto 71
   itype = iwttyp(i)
   
   ! Modify the weight of the appropriate term. Changes are step-wise,
   ! every IWTSTP(3,I) steps (IWTSTP(3,I)=1, if user did not specify it).
   
   ! If IWTSTP(3,I) > 0, then the weight is linearly interpolated.
   ! If IWTSTP(3,I) < 0, then the weight is modified by a multplicative factor.
   
   if (iwtstp(3,i) > 0) then
      nstepu = nstep - mod(nstep-iwtstp(1,i),iwtstp(3,i))
      wt = nstepu*wtnrg(1,i) + wtnrg(2,i)
   else if (iwtstp(3,i) < 0) then
      nstepu = (nstep-iwtstp(1,i))/abs(iwtstp(3,i))
      wt = wtnrg(2,i) * wtnrg(1,i)**nstepu
   end if
   wt2 = wt
   
   ! TYPE = RSTAR:
   ! In this case, set WT to the appropriate R**12 weight; set WT2 to the
   ! appropriate attractive R**6 weight, and set ITYPE=4 (VDW). Also set
   ! IRSTR = 1 and WT3 = WT**10. This will ensure h-bond parameters will
   ! also be appropriately modified. The actual weight changes will be handled
   ! by the subsequent loops.
   
   irstr = 0
   if (itype == 14) then
      irstr = 1
      itype = 4
      wt2 = wt**6.0d0
      wt3 = wt**10.0d0
      wt = wt2*wt2
   end if
   
   
   ! TYPE = BOND/INTERN/ALL:
   
   25 if ((itype == 1 .or. itype == 10 .or. &
         itype == 11 .or. ireset(1) == 1) .and.ichold(1) /= 2) then
      rmult = wt/weight(1)
      rk(1:numbnd) = rk(1:numbnd)*rmult
      weight(1) = wt
      ichang(1) = 1
      if (abs(wt) <= small) ichold(1) = 2
      if (fixed) ichold(1) = 2
   end if
   
   ! TYPE = ANGLE/INTERN/ALL:
   
   if ((itype == 2 .or. itype == 10 .or. &
         itype == 11.or. ireset(2) == 1) .and. ichold(2) /= 2) then
      rmult = wt/weight(2)
      tk(1:numang) = tk(1:numang)*rmult
      weight(2) = wt
      ichang(2) = 1
      if (abs(wt) <= small) ichold(2) = 2
      if (fixed) ichold(2) = 2
   end if
   
   ! TYPE = TORSION/INTERN/ALL:
   
   if ((itype == 3 .or. itype == 10 .or. &
         itype == 11.or. ireset(3) == 1) .and. ichold(3) /= 2) then
      rmult = wt/weight(3)
      pk(1:numphi) = pk(1:numphi)*rmult
      weight(3) = wt
      ichang(3) = 1
      if (abs(wt) <= small) ichold(3) = 2
      if (fixed) ichold(3) = 2
      itorfc = 1
   end if
   
   ! TYPE = VDW/ATTRACT/NB/ALL:
   
   if ((itype == 4 .or. itype == 7 .or. itype == 8 .or. &
         itype == 11.or. ireset(4) == 1) .and. ichold(4) /= 2) then
      rmult = wt2/weight(4)
      do j=1,nttyp
         cn2(j) = cn2(j)*rmult
      end do
      weight(4) = wt2
      ichang(4) = 1
      if (abs(wt2) <= small) ichold(4) = 2
      if (fixed) ichold(4) = 2
      ivdw = 1
   end if
   
   ! TYPE = VDW/REPULSE/NB/ALL:
   
   if ((itype == 4 .or. itype == 7 .or. itype == 9 .or. &
         itype == 11.or. ireset(5) == 1) .and. ichold(5) /= 2) then
      rmult = wt/weight(5)
      do j=1,nttyp
         cn1(j) = cn1(j)*rmult
      end do
      weight(5) = wt
      ichang(5) = 1
      if (abs(wt) <= small) ichold(5) = 2
      if (fixed) ichold(5) = 2
      ivdw = 1
   end if
   
   ! TYPE = HB/ATTRACT/NB/ALL:
   
   if ((irstr == 1 .or. itype == 5 .or. itype == 7 .or. &
         itype == 8 .or. itype == 11 .or. ireset(6) == 1) .and. &
         ichold(6) /= 2) then
      
      rmult = wt/weight(6)
      if (irstr == 1) rmult = wt3/weight(6)
      do j=1,nhb
         bg(j) = bg(j)*rmult
      end do
      weight(6) = wt
      ichang(6) = 1
      if (abs(wt) <= small) ichold(6) = 2
      if (fixed) ichold(6) = 2
      ihbp = 1
   end if
   
   ! TYPE = HB/REPULSE/NB/ALL:
   
   if ((irstr == 1 .or. itype == 5 .or. itype == 7 .or. &
         itype == 9 .or. itype == 11 .or. ireset(7) == 1) .and. &
         ichold(7) /= 2) then
      
      rmult = wt/weight(7)
      do j=1,nhb
         ag(j) = ag(j)*rmult
      end do
      weight(7) = wt
      ichang(7) = 1
      if (abs(wt) <= small) ichold(7) = 2
      if (fixed) ichold(7) = 2
      ihbp = 1
   end if
   
   ! TYPE = ELEC/NB/ALL:
   
   if ((itype == 6 .or. itype == 7 .or. &
         itype == 11.or. ireset(8) == 1) .and. ichold(8) /= 2) then
      rmult = sqrt(wt/weight(8))
      do j=1,ncharg
         cg(j) = cg(j)*rmult
      end do
      weight(8) = wt
      ichang(8) = 1
      if (abs(wt) <= small) ichold(8) = 2
      if (fixed) ichold(8) = 2
   end if
   
   ! TYPE = IMPROP
   
   if ((itype == 15.or.ireset(9) == 1).and.ichold(9) /= 2) then
      rmult = wt/weight(9)
      do j=1,nimprp
         pk(numphi+j) = pk(numphi+j)*rmult
      end do
      weight(9) = wt
      ichang(9) = 1
      if (abs(wt) <= small) ichold(9) = 2
      if (fixed) ichold(9) = 2
      itorfc = 1
   end if
   
   ! TYPE = REST/RESTS:
   
   if ((itype == 12.or.itype == 19.or.ireset(10) == 1) .and. &
         ishrtb > 0 .and. ichold(10) /= 2) then
      rmult = wt/weight(10)
      do j=1,nmrnum
         if (nmrshb(j) > 0) then
            rk2nmr(2,j) = rk2nmr(2,j)*rmult
            rk3nmr(2,j) = rk3nmr(2,j)*rmult
            if (nmrst(3,j) > 0) then
               rk2nmr(1,j) = rk2nmr(1,j)*rmult
               rk3nmr(1,j) = rk3nmr(1,j)*rmult
            end if
         end if
      end do
      weight(10) = wt
      ichang(10) = 1
      if (abs(wt) <= small) ichold(10) = 2
      if (fixed) ichold(10) = 2
   end if
   
   ! TYPE = REST/RESTL:
   
   if ((itype == 12.or.itype == 20.or.ireset(11) == 1) .and. &
         ichold(11) /= 2) then
      rmult = wt/weight(11)
      do j=1,nmrnum
         if (nmrshb(j) <= 0) then
            rk2nmr(2,j) = rk2nmr(2,j)*rmult
            rk3nmr(2,j) = rk3nmr(2,j)*rmult
            if (nmrst(3,j) > 0) then
               rk2nmr(1,j) = rk2nmr(1,j)*rmult
               rk3nmr(1,j) = rk3nmr(1,j)*rmult
            end if
         end if
      end do
      
      weight(11) = wt
      ichang(11) = 1
      if (abs(wt) <= small) ichold(11) = 2
      if (fixed) ichold(11) = 2
   end if
   
   ! TYPE = NOESY
   
   if ((itype == 21.or.ireset(12) == 1).and.ichold(12) /= 2) then
      weight(12) = wt
      ichang(12) = 1
      if (abs(wt) < small .or. fixed) ichold(12) = 2
   end if
   
   ! TYPE = SHIFTS
   
   if ((itype == 22.or.ireset(13) == 1).and.ichold(13) /= 2) then
      weight(13) = wt
      ichang(13) = 1
      if (abs(wt) < small .or. fixed) ichold(13) = 2
   end if
   
   ! TYPE = SOFTR
   
   if (itype == 16) then
      rwell = wt
      ichang(15) = 1
   end if
   
   
   ! TYPE = TEMP0
   
   if (itype == 13) then
      temp0 = wt
      ichang(16) = 1
   end if
   
   ! TYPE = TAUTP
   
   if (itype == 23) then
      tautp = wt
      ichang(14) = 1
   end if
   
   ! TYPE = CUTOFF
   
   if (itype == 17) then
      cut = wt**2
      ichang(17) = 1
   end if
   
   ! TYPE = DISAVE, ANGAVE, TORAVE, PLPTAVE, PLNAVE, or GDISAVE
   
   if ((itype >= 24 .and. itype <= 26) .or. (itype >= 28 .and. itype <= 30)) then
      indx = itype-23
      if (indx > 3) indx = indx - 1
      tauave(indx) = wt
      if(itype >= 24 .and. itype <= 26) then
        ichang(17+indx) = 1
      else
        ichang(17+indx+1) = 1
      end if
   end if
   
   ! TYPE = TGTRMSD
   
   if (itype == 27) then
      tgtrmsd = wt
      ichang(21) = 1
   end if
   
   ! TYPE = TEMP0LES
   
   if (itype == 31) then
      temp0les = wt
      ichang(24) = 1
   end if
   
   ! End of implied loop over ICHGWT instructions:
   
   71 i = i + 1
   if (i <= ichgwt .and. irstyp == 0) goto 20
   
   ! Now see what weights were modified. Any that were not modified
   ! on this call, and that were not set to 1.0 previously, should be reset
   ! to 1.0.
   
   75 irstyp = 1
   do i = 1,14
      if (ireset(i) == 1) then
         cycle
      else if (ichang(i) == 0 .and. ichold(i) == 1) then
         ireset(i) = 1
         wt = ONE
         wt2 = ONE
         ichold(i) = 0
         goto 25
      else if (ichang(i) == 1 .and. ichold(i) == 0) then
         ichold(i) = 1
      end if
   end do
   
   wnoesy = weight(12)
   wshift = weight(13)
   if (ichang(14) /= 1) tautp = weight(14)
   if (ichang(15) /= 1) rwell = weight(15)
   if (ichang(16) /= 1) temp0 = weight(16)
   if (ichang(17) /= 1) cut = weight(17)
   if (ichang(18) /= 1) tauave(1) = weight(18)
   if (ichang(19) /= 1) tauave(2) = weight(19)
   if (ichang(20) /= 1) tauave(3) = weight(20)
   if (ichang(21) /= 1) tgtrmsd = weight(21)
   if (ichang(22) /= 1) tauave(4) = weight(22)
   if (ichang(23) /= 1) tauave(5) = weight(23)
   if (ichang(24) /= 1) tauave(6) = weight(24)
   if (ichang(25) /= 1) temp0les = weight(25)
   
   wtls(1) = weight(10)
   wtls(2) = weight(11)
   
   ! If the vdw parameters were changed, call nmrrad to re-set the values
   ! of the radii and well-depths
   
   if (ivdw == 1) call nmrrad(rad,wel,cn1,cn2,ntypes,0,0.0d0)
   
   ! DECNVH resets the h-bond r* and epsilon parameters when the h-bond
   ! coefficients have been changed.
   
   if (ihbp == 1) call decnvh(ag,bg,nhb,radhb,welhb)
   
   ! SETGMS sets a couple of arrays which depend on the torsional term
   ! force constants, if those force constants were changed here.
   
   if (itorfc == 1) call setgms(numphi+nimprp)
   
   return
end subroutine modwt 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ndvprt here]
subroutine ndvprt(x,f,name,irsnam,ipres,nres,iscrth,natom,rimass, &
      ntb,nmrnum,nstep,nmrat,resttype,rstwtarr,kxpk,nmrst,r1nmr,r2nmr,r3nmr, &
      r4nmr,rk2nmr,rk3nmr,nmrcom,nmrfty,igravt, &
      ialtdis,bave, &
      bave0,aave,aave0,tave1,tave01,tave2,tave02,ptave,ptave0, &
      plave,plave0,gdave,gdave0,ajcoef,nave,nexact,ipower,tauave,ravein, &
      navint,dt,iscopn,iout)
   
   
   ! Subroutine Nmr DeViations PRinT.
   
   ! This routine prints the deviations of each restraint from its
   ! target value, and the energy contribution for that restraint.
   
   ! Author: David A. Pearlman
   ! Date: 7/89
   
   ! INPUT:
   
   !  X(I)       : Coordinate array
   !  F(I)       : Force array; Updated during this call.
   !  NAME(I)    : Atom name array
   !  IRSNAM()   : Residue name array
   !  IPRES()    : Residue I contains atoms IPRES(I) -> IPRES(I+1)-1
   !  NRES       : Number of residues in the system
   !  ISCRTH()   : Scratch array of length at least NATOM.
   !  NATOM      : Number of atoms in the system.
   !  RIMASS(I)  : Array of inverse masses
   !  NTB        : Periodic boundary conditions flag
   !  NMRNUM     : The number of NMR restraints
   !  NSTEP      : The step/iteration number in the main calling program
   !  NMRAT(16,I) : The 2-8 atoms defining restraint I.
   !                NMRAT(9,I) through NMRAT(16,I) are used for atom groups
   !  KXPK(2,I)  : integer peak identifier, just for printing 1:peak 2:dataset
   !  NMRST(3,I) : Restraint I is imposed over the range of steps
   !               NMRST(1,I) -> NMRST(2,I). If NMRST(3,I) .NE. 0, change
   !               is implemented as a step function every ABS(NMRST(3,I)) steps.
   !               Otherwise, changes in target values occur continuously.
   !  R1NMR(2,I) : The slope/intercept of the dependence of r1 on step
   !               number (see routine DISNRG,ANGNRG or TORNRG for description
   !               of r1, r2, r3, r4, k2, and k3).
   !  R2NMR(2,I) : The slope/int. of the dependence of r2 on step number.
   !  R3NMR(2,I) : The slope/int. of the dependence of r3 on step number.
   !  R4NMR(2,I) : The slope/int. of the dependence of r4 on step number.
   !  RK2NMR(2,I) : The slope/int. of the dependence of k2 on step number.
   !  RK3NMR(2,I) : The slope/int. of the dependence of k3 on step number.
   !  NMRCOM(2,I) : The ranges of atoms defining the center-of-mass group,
   !                for applicable distance restraints
   !  NMRFTY(I)  : Flag for each restraint indicating what functional form should
   !               be used for that restraint. If time-averaged restraints
   !               are being applied to a particular class of internal
   !               (bond, angle, torsion) and NMRFTY(I)=1, then an instantaneous
   !               restraint will be applied to internal I.
   !  IGRAVT(I)  : For averaged group positions (if requested)
   !               = 0 for center of mass
   !               = 1 for r**-6 average
   !  IALTDIS(I) : flag to use alternate restraint functional form
   !    BAVE(I)
   !   BAVE0(I)
   !    AAVE(I)
   !   AAVE0(I)    For time-averaged restraints. _AVE(I) contains the time-averaged
   !   TAVE1(I)    value of the restraint, using the exponential decay value in
   !   TAVE01(I)   TAUAVE(I). _AVE0(I) contains the time-averaged value of the
   !   TAVE2(I)    restraint, using no exponential decay (real average). For
   !   TAVE02(I)   torsions, two values are required for average, corresponding
   !   PTAVE(I)    to the COS(tau) and SIN(tau)
   !  PTAVE0(I)
   !   PLAVE(I)
   !  PLAVE0(I)
   
   ! AJCOEF(3,I) : The coefficients to be used in the Karplus equation for
   !               any torsional restraints to be imposed using J-coupling.
   
   !     NAVE(3) : For time-averaged restraints. >0 for time-averaged restraints.
   !               Elements 1->3 correspond to values for bonds, angles, torsions.
   !   NEXACT(3) : Not currently used.
   !   IPOWER(3) : For time-averaged restraints. Gives the exponent used in
   !               developing the time average.
   !   TAUAVE(3) : For time-averaged restraints. Gives the time constant used
   !               in the exponential decay multiplier.
   !   RAVEIN(3) : For time-averaged restraints. On the first step, the value
   !               returned is the current value+RAVEIN(I) (I=1,3 for bonds,
   !               angles, torsions).
   !   NAVINT(3) : The averaged value is only recalculated every NAVINT() steps.
   !          DT : Integration timestep (ps). Used for time-averaged restraints.
   !  ISCOPN     : Scratch unit for redirected prints (if requested)
   !  IOUT       : Unit for prints (if requested to POUT)
   
   use constants, only : RAD_TO_DEG 
   use file_io_dat
   implicit none
   integer:: i, ialtdis, iat1, iat2, iat3, iat4, iat5, iat6, iat7, &
        iat8, iave, igravt, iout, ipower, ipres, iscopn, iscrth, itimes, &
        j, kxpk, natom, nave, navint, nexact, nmrat, nmrcom, nmrfty, &
        nmrnum, nmrst, nres, nstep, nstepu, ntb
   _REAL_ :: aave, aave0, ajcoef, bave, bave0, bound, convrt, dev, &
        dev2, dfr, dt, e, f, gdave, gdave0, plave, plave0, ptave, ptave0, &
        r1, r1nmr, r2, r2nmr, r3, r3nmr, r4, r4nmr, ravein, rimass, rint, &
        rk2, rk2nmr, rk3, rk3nmr, rmstot, rstwtarr, small, step, target, &
        tauave, tave01, tave02, tave1, tave2, x, xcom
   character(len=80) line
   logical jcoupl,usecom
   dimension x(*),f(*),rimass(*),nmrat(16,*),rstwtarr(4,*),nmrst(3,*),r1nmr(2,*)
   dimension r2nmr(2,*),r3nmr(2,*),r4nmr(2,*),rk2nmr(2,*),kxpk(2,*)
   character(len=14) tmpnam(4)
   character(len=4) name,irsnam
   integer resttypetmp
   integer resttype(*)
   
#  include "nmr.h"
   
   dimension rk3nmr(2,*),nmrcom(2,*),nmrfty(*),igravt(*),name(*)
   dimension irsnam(*),ipres(*),iscrth(*),ialtdis(*)
   dimension xcom(24),dfr(24)
   dimension bave(*),bave0(*),aave(*),aave0(*),tave1(*)
   dimension tave01(*),tave2(*),tave02(*),ptave(*),ptave0(*),plave(*),plave0(2),gdave(*),gdave0(*),&
             nave(6),rmstot(8)
   dimension nexact(6),ipower(6),tauave(6),ravein(6),navint(6)
   dimension ajcoef(3,*)
   
   save itimes,small
   data itimes/0/
   data small/1.0d-7/
   
   integer jjj
   _REAL_ ricmmass(8)
   _REAL_ comrimass, rstwt,edis,eang,etor,eplpt,epln,egendis
   dimension comrimass(4), rstwt(4)
   integer iat
   dimension iat(8)
   
   equivalence (iat(1),iat1),(iat(2),iat2)
   equivalence (iat(3),iat3),(iat(4),iat4)
   equivalence (iat(5),iat5),(iat(6),iat6)
   equivalence (iat(7),iat7),(iat(8),iat8)
   
   do jjj=1,4
     comrimass(jjj) = 0.0
   end do

   ! First, figure out where this information should be output (based on
   ! user-defined redirection).
   
   iuse = iout
   if (iredir(2) == 0) then
      return
   else
      if (redir(2)(1:iredir(2)) == 'POUT') then
         goto 1
         
         ! If LISTIN file = LISTOUT file, or if more than one LISTOUT type dump
         ! is being done, open the LISTOUT file as "old" 
         
      else if ((iredir(1) == iredir(2) .and. &
            redir(1)(1:iredir(1)) == redir(2)(1:iredir(2))) &
            .or. itimes > 0) then
         call amopen(iscopn, redir(2)(1:iredir(2)),'O','F','W')
      else
         call amopen(iscopn, redir(2)(1:iredir(2)), owrite, 'F','W')
      end if
   end if
   
   write(iout,9020) redir(2)(1:iredir(2))
   iuse = iscopn
   1 write(iuse,9032)
   write(iuse,9077) restrt(1:40)
   write(iuse,9030) pencut
   write(iuse,9032)
   write(iuse,9031)
   write(iuse,9032)
   
   ! Set up pointers to the residue number corresponding to each atom in
   ! the ISCRTH array.
   
   do i = 1,nres
      do j=ipres(i),ipres(i+1)-1
         iscrth(j) = i
      end do
   end do
   
   ! Main loop over all the restraints:
   
   edis = 0.0d0
   eang = 0.0d0
   etor = 0.0d0
   eplpt = 0.0d0
   epln = 0.0d0
   egendis = 0.0d0
   do i=1,nmrnum
      
      ! Skip restraints if the current step # is outside the applicable range:
      
      if (nstep < nmrst(1,i) .or. (nstep > nmrst(2,i) .and. &
            nmrst(2,i) > 0)) cycle
      
      ! Calculate the values of r1,r2,r3,r4,k2, and k3 for this step & restraint:
      
      ! Vary the step used in calculating the values only every NMRST(3,I) steps.
      ! If the user did not specify NINC in the input, NMRST(3,I)=1.
      
      nstepu = nstep - mod(nstep-nmrst(1,i),abs(nmrst(3,i)))
      step = nstepu
      r1 = r1nmr(1,i)*step + r1nmr(2,i)
      r2 = r2nmr(1,i)*step + r2nmr(2,i)
      r3 = r3nmr(1,i)*step + r3nmr(2,i)
      r4 = r4nmr(1,i)*step + r4nmr(2,i)
      
      ! If NMRST(3,I) > 0, then the weights are linearly interpolated.
      ! If NMRST(3,I) < 0, then the weights are modified by a multplicative factor.
      
      if (nmrst(3,i) > 0) then
         rk2 = rk2nmr(1,i)*step + rk2nmr(2,i)
         rk3 = rk3nmr(1,i)*step + rk3nmr(2,i)
      else
         nstepu = (nstep-nmrst(1,i))/abs(nmrst(3,i))
         rk2 = rk2nmr(2,i) * rk2nmr(1,i)**nstepu
         rk3 = rk3nmr(2,i) * rk3nmr(1,i)**nstepu
      end if
      
      ! Take the average of r2 & r3 as the average "target" value:
      
      target = (r2+r3)/2.0d0
      
      ! Determine what type of restraint (distance,angle,torsion) this
      ! is, and call the appropriate routine to calculate its contribution.
      
      usecom = .false.
      do jjj=9,16
        usecom = usecom .or. nmrat(jjj,i) < 0
      end do
      
      iave = 0
      jcoupl = .false.
      resttypetmp = resttype(i)
      select case (resttypetmp)
       case (1)
         if (nave(1) > 0 .and. nmrfty(i) == 0) iave = 1
         if (.not. usecom) then
            call disnrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i),e,r1, &
                  r2,r3,r4,rk2,rk3,ntb,bave(2*i-1), &
                  bave0(i),nave(1),nexact(1),ipower(1), &
                  tauave(1),ravein(1),dt,navint(1),iave,0, &
                  0,2,ialtdis(i))
         else if (igravt(i) == 1) then
            call r6ave(x,nmrat,nmrcom,rint,i)
            call disnrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i),e,r1, &
                  r2,r3,r4,rk2,rk3,ntb,bave(2*i-1), &
                  bave0(i),nave(1),nexact(1),ipower(1), &
                  tauave(1),ravein(1),dt,navint(1),iave,0, &
                  0,3,ialtdis(i))
         else
            call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
            call disnrg(xcom,f,dfr,rint,0,3,e,r1,r2,r3, &
                  r4,rk2,rk3,ntb,bave(2*i-1),bave0(i), &
                  nave(1),nexact(1),ipower(1),tauave(1), &
                  ravein(1),dt,navint(1),iave,0,0,2,ialtdis(i))
         end if
         edis = edis + e
       case (2)
         if (nave(2) > 0 .and. nmrfty(i) == 0) iave = 1
         if (usecom) then
           call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
           call angnrg(xcom,f,dfr,rint,0,3,6, &
               e,r1,r2,r3,r4,rk2,rk3,ntb,aave(2*i-1), &
               aave0(i),nave(2),nexact(2),ipower(2), &
               tauave(2),ravein(2),dt,navint(2),iave,0,0,2)
         else
           call angnrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i),nmrat(3,i), &
               e,r1,r2,r3,r4,rk2,rk3,ntb,aave(2*i-1), &
               aave0(i),nave(2),nexact(2),ipower(2), &
               tauave(2),ravein(2),dt,navint(2),iave,0,0,2)
         end if
         eang = eang + e
       case (3)
         if (nave(3) > 0 .and. nmrfty(i) == 0) iave = 1
         
         ! If J coefficients are 0 for this torsion, it is a normal torsional
         ! restraint. If they are .NE.0, then it is a J-coupling factor restraint.
         
         jcoupl = abs(ajcoef(1,i)) > small .or. &
               abs(ajcoef(2,i)) > small
         
         if (.not.jcoupl) then
           if (usecom) then
             call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
             call tornrg(x,f,dfr,rint,0,3,6,9, &
                  e,r1,r2,r3,r4,rk2,rk3,ntb, &
                  tave1(2*i-1),tave01(i),tave2(2*i-1), &
                  tave02(i),nave(3),nexact(3),ipower(3), &
                  tauave(3),ravein(3),dt,navint(3),iave,0,0,2)
           else
             call tornrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i),nmrat(3,i), &
                  nmrat(4,i),e,r1,r2,r3,r4,rk2,rk3,ntb, &
                  tave1(2*i-1),tave01(i),tave2(2*i-1), &
                  tave02(i),nave(3),nexact(3),ipower(3), &
                  tauave(3),ravein(3),dt,navint(3),iave,0,0,2)
           end if
         else
           if (usecom) then
             ! Since atom groups are not supported for j-coupling restraints, 
             ! print out an error message and terminate execution if the user 
             ! specifies atom groups AND j-coupling
             
             write(iout,9086)
             call mexit(iout, 1)
           else
             call jnrg(x,f,rint,nmrat(1,i),nmrat(2,i),nmrat(3,i), &
                  nmrat(4,i),e,r1,r2,r3,r4,rk2,rk3,ntb, &
                  tave1(2*i-1),tave01(i),tave2(2*i-1), &
                  tave02(i),nave(3),nexact(3),ipower(3), &
                  tauave(3),ravein(3),dt,navint(3),iave,0,0, &
                  ajcoef(1,i),2)
           end if
         end if
         etor = etor + e
       case (4)  ! This is a plane-point angle restraint
         if (nave(4) > 0 .and. nmrfty(i) == 0) iave = 1
         if (usecom) then
           call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
           call commass(comrimass,nmrat,nmrcom,rimass,i)
           call plptnrg(xcom,f,dfr,comrimass,rint,0,3,6,9,12, &
               e,r1,r2,r3,r4,rk2,rk3,ntb,ptave(2*i-1), &
               ptave0(i),nave(4),nexact(4),ipower(4), &
               tauave(4),ravein(4),dt,navint(4),iave, &
               0,0,2)
         else
           call plptnrg(x,f,dfr,rimass,rint,nmrat(1,i),nmrat(2,i), &
               nmrat(3,i),nmrat(4,i),nmrat(5,i), &
               e,r1,r2,r3,r4,rk2,rk3,ntb,ptave(2*i-1), &
               ptave0(i),nave(4),nexact(4),ipower(4), &
               tauave(4),ravein(4),dt,navint(4),iave, &
               0,0,2)
         end if
         eplpt = eplpt + e
       case (5)  ! This is a plane-plane angle restraint
         if (nave(5) > 0 .and. nmrfty(i) == 0) iave = 1
         if (usecom) then
           call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
           call plnnrg(xcom,f,dfr,rint,0,3,6,9,12,15,18,21, &
               e,r1,r2,r3,r4,rk2,rk3,ntb,plave(2*i-1), &
               plave0(i),nave(5),nexact(5),ipower(5), &
               tauave(5),ravein(5),dt,navint(5),iave, &
               0,0,2)
         else
           call plnnrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i), &
               nmrat(3,i),nmrat(4,i),nmrat(5,i),nmrat(6,i),nmrat(7,i),nmrat(8,i), &
               e,r1,r2,r3,r4,rk2,rk3,ntb,plave(2*i-1), &
               plave0(i),nave(5),nexact(5),ipower(5), &
               tauave(5),ravein(5),dt,navint(5),iave, &
               0,0,2)
         end if
         epln = epln + e
       case (6) ! This is a generalized distance coordinate restraint
        do jjj=1,4
          rstwt(jjj) = rstwtarr(jjj,i)
        end do
        if (nave(6) > 0 .and. nmrfty(i) == 0) iave = 1
        if (usecom) then
          call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
          call gendisnrg(xcom,f,dfr,rint,0,3,6,9,12,15,18,21, &
               rstwt,e,r1,r2,r3,r4,rk2,rk3,ntb,gdave(2*i-1), &
               gdave0(i),nave(6),nexact(6),ipower(6), &
               tauave(6),ravein(6),dt,navint(6),iave, &
               0,0,2)
        else
          call gendisnrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i), &
               nmrat(3,i),nmrat(4,i),nmrat(5,i),nmrat(6,i),nmrat(7,i),nmrat(8,i), &
               rstwt,e,r1,r2,r3,r4,rk2,rk3,ntb,gdave(2*i-1), &
               gdave0(i),nave(6),nexact(6),ipower(6), &
               tauave(6),ravein(6),dt,navint(6),iave, &
               0,0,2)
        end if
        egendis = egendis + e
       case default
        write(iout, 9087)
        call mexit(iout, 1)
      end select
      
      !  --- if the energy of this restraint is less than the cutoff,
      !      do not print anything
      
      if (e < pencut) cycle
      
      ! Calculate the deviation of the current value of the internal from
      ! its target
      
      dev = abs(target - rint)
      if (rint >= r2 .and. rint <= r3) then
         dev2 = 0.0d0
         !                              --if no penalty, print upper bound:
         bound = r3
      else if (rint < r2) then
         dev2 = r2 - rint
         bound = r2
      else
         dev2 = rint - r3
         bound = r3
      end if
      
      ! Print the internal, its value, its deviations from the target,
      ! and its energy contribution:
      
      
      do jjj=1,8
        iat(jjj) = nmrat(jjj,i)/3 + 1
        if(nmrat(jjj,i) < 0 .and. nmrat(jjj+8,i) < 0) iat(jjj) = nmrcom(2,-nmrat(jjj+8,i))
      end do
      
      usecom = .false.

      
      !  --- convention for printing: if a group is defined, identify it
      !        by its final atom, and put a "*" in the first column.
      !        This is usually sufficient, e.g. if
      !        the group is the three atoms of a methyl group, or two
      !        equivalent delta or epsilon protons on an aromatic ring,
      !        etc.  For more exotic groupings of atoms, this will not
      !        be as informative as it could be, but we will cross that
      !        bridge if it turns out to be a real problem.
      
      !        Also add a one-letter code at the end of each line:
      !              d       for distance constraint
      !              a       for angle constraint
      !              t       for torsion constraint
      !              j       for j-coupling constraint
      !              p       for plane-point restraint
      !              n       for plane-plane restraint
      !              g       for gemeralized distance coordinate restraint
      

      convrt = 1.0d0
      if (resttype(i) /= 1 .and. resttype(i) /= 6 .and. &
            .not.jcoupl) convrt = RAD_TO_DEG
      
      resttypetmp = resttype(i)
      select case (resttypetmp)
        case (1)
         write(tmpnam(1)(2:14),'(A4,1x,A4,I4)') &
               name(iat1),irsnam(iscrth(iat1)),iscrth(iat1)
         write(tmpnam(2)(2:14),'(A4,1x,A4,I4)') &
               name(iat2),irsnam(iscrth(iat2)),iscrth(iat2)
         tmpnam(1)(1:1) = ' '
         tmpnam(2)(1:1) = ' '
         if (nmrat(1,i) < 0) tmpnam(1)(1:1) = '*'
         if (nmrat(2,i) < 0) tmpnam(2)(1:1) = '*'

         write(iuse,9073) tmpnam(1),tmpnam(2), &
               rint,bound,abs(dev2),e,kxpk(1,i),kxpk(2,i)

       case (2)
         write(tmpnam(1)(2:14),'(A4,1x,A4,I4)') &
               name(iat1),irsnam(iscrth(iat1)),iscrth(iat1)
         write(tmpnam(3)(2:14),'(A4,1x,A4,I4)') &
               name(iat3),irsnam(iscrth(iat3)),iscrth(iat3)
         tmpnam(1)(1:1) = ' '
         tmpnam(3)(1:1) = ' '
         if (nmrat(1,i) < 0) tmpnam(1)(1:1) = '*'
         if (nmrat(3,i) < 0) tmpnam(3)(1:1) = '*'
         write(iuse,9041) tmpnam(1),tmpnam(3), &
               convrt*rint,convrt*bound,convrt*abs(dev2),e
       case (3)
         
         ! --- for torsions, print the two central atoms and the deviations:
         !       (for J-couplings, print the two end atoms, since this is generally
         !        more informative.)
         
         write(tmpnam(1)(2:14),'(A4,1x,A4,I4)') &
               name(iat1),irsnam(iscrth(iat1)),iscrth(iat1)
         write(tmpnam(2)(2:14),'(A4,1x,A4,I4)') &
               name(iat2),irsnam(iscrth(iat2)),iscrth(iat2)
         write(tmpnam(3)(2:14),'(A4,1x,A4,I4)') &
               name(iat3),irsnam(iscrth(iat3)),iscrth(iat3)
         write(tmpnam(4)(2:14),'(A4,1x,A4,I4)') &
               name(iat4),irsnam(iscrth(iat4)),iscrth(iat4)
         tmpnam(1)(1:1) = ' '
         tmpnam(2)(1:1) = ' '
         tmpnam(3)(1:1) = ' '
         tmpnam(4)(1:1) = ' '
         if (nmrat(1,i) < 0) tmpnam(1)(1:1) = '*'
         if (nmrat(2,i) < 0) tmpnam(2)(1:1) = '*'
         if (nmrat(3,i) < 0) tmpnam(3)(1:1) = '*'
         if (nmrat(4,i) < 0) tmpnam(4)(1:1) = '*'
         if( jcoupl ) then
            write(iuse,9082) tmpnam(1),tmpnam(4), &
                  convrt*rint,convrt*bound,convrt*abs(dev2),e
         else
            write(iuse,9083) tmpnam(2),tmpnam(3), &
                  convrt*rint,convrt*bound,convrt*abs(dev2),e
         end if
       case (4:5)
        ! For planes, identify each plane by the first atom in the plane
        ! For plane-point restraints, the second atom printed will correspond
        ! to the atom or group that makes up the point.
        
        write(tmpnam(1)(2:14),'(A4,1x,A4,I4)') &
               name(iat1),irsnam(iscrth(iat1)),iscrth(iat1)
        write(tmpnam(2)(2:14),'(A4,1x,A4,I4)') &
               name(iat5),irsnam(iscrth(iat5)),iscrth(iat5)
        tmpnam(1)(1:1) = ' '
        tmpnam(2)(1:1) = ' '
        if (nmrat(1,i) < 0) tmpnam(1)(1:1) = '*'
        if (nmrat(5,i) < 0) tmpnam(2)(1:1) = '*'
        if (resttype(i) == 4) then
          write(iuse,9089) tmpnam(1),tmpnam(2), &
                  convrt*rint,convrt*bound,convrt*abs(dev2),e
        else
          write(iuse,9090) tmpnam(1),tmpnam(2), &
                  convrt*rint,convrt*bound,convrt*abs(dev2),e
        end if
       case (6)
        ! For generalized distance coordinates, identify the restraint by the
        ! first atom of the first distance and the first atom of the second distance.
        ! This is pretty inadaquate, but there's not really a good way to do this
        ! in a way that resembles the output format of the other restraints.
        
        write(tmpnam(1)(2:14),'(A4,1x,A4,I4)') &
               name(iat1),irsnam(iscrth(iat1)),iscrth(iat1)
        write(tmpnam(2)(2:14),'(A4,1x,A4,I4)') &
               name(iat3),irsnam(iscrth(iat3)),iscrth(iat3)
        tmpnam(1)(1:1) = ' '
        tmpnam(2)(1:1) = ' '
        if (nmrat(1,i) < 0) tmpnam(1)(1:1) = '*'
        if (nmrat(3,i) < 0) tmpnam(2)(1:1) = '*'
        write(iuse,9091) tmpnam(1),tmpnam(2), &
                  convrt*rint,convrt*bound,convrt*abs(dev2),e
        
      end select
      
      ! If a group-type distance restraint was defined, echo the sorted group:
      
   end do

   if (edis > 0.0) write(iuse,9074) edis
   if (eang > 0.0) write(iuse,9075) eang
   if (etor > 0.0) write(iuse,9076) etor
   if (eplpt > 0.0) write(iuse,9078) eplpt
   if (epln > 0.0) write(iuse,9079) epln
   if (egendis > 0.0) write(iuse,9080) egendis
#ifndef MPI
   write(iuse,9084) ebdev
   write(iuse,9085) eadev
#endif
   write(iuse,9032)
   
   itimes = itimes + 1
#ifndef MPI
   
   !  ---- if NOESY volumes were calculated, print a summary here:
   
   if (iredir(4) /= 0) then
      
      write(iuse,9071)
      write(iuse,9032)
      rewind (81)
      11 read(81,'(a80)',end=12) line
      write(iuse,'(a80)') line
      goto 11
      12 rewind (81)
      write(iuse,9032)
   end if
#endif
   
   !  ---- if chemical shifts were calculated, print a summary here:
   
   if (iredir(5) /= 0) then
      write(iuse,9072)
      rewind (82)
      13 read(82,'(a80)',end=14) line
      write(iuse,'(a80)') line
      goto 13
      14 rewind (82)
      write(iuse,9032)
   end if
   
   !  ---- if paramagnetic shifts were calculated, print a summary here:

   if (iredir(7) /= 0) then
      write(iuse,9092)
      rewind (53)
      15 read(53,'(a80)',end=16) line
      write(iuse,'(a80)') line
      goto 15
      16 rewind (53)
      write(iuse,9032)
   end if
   
   !  ---- if RDC's were calculated, print a summary here:

   if (iredir(8) /= 0) then
      write(iuse,9093)
      rewind (57)
      17 read(57,'(a80)',end=18) line
      write(iuse,'(a80)') line
      goto 17
      18 rewind (57)
      write(iuse,9032)
   end if

   !  ---- if CSA restraints were calculated, print a summary here:

   if (iredir(9) /= 0) then
      write(iuse,9094)
      rewind (58)
      19 read(58,'(a80)',end=20) line
      write(iuse,'(a80)') line
      goto 19
      20 rewind (58)
      write(iuse,9032)
   end if

   if (iuse == iscopn) close(iuse)
   return
   9020 format(' Restraints/deviations being written to file: ',a)
   9030 format(' Restraints, deviations, and energy contributions:', &
         '    pencut = ',f7.2/)
   9031 format( &
         '     First atom        Last atom    curr. value target', &
         ' deviation  penalty')
   9032 format(' ',78('-'))
   9041 format(' ',a14,' -- ',a14,':',4f9.3,' a')
   9071 format(/ /' NOESY Volume analysis:'/ /)
   9072 format(/ /'  Chemical shifts analysis:'/ /)
   9092 format(/ /'  Paramagnetic shifts analysis:'/ /)
   9093 format(/ /'  Residual dipolar splittings:'/ /)
   9094 format(/ /'  Residual CSA splittings:'/ /)
   9073 format(' ',a14,' -- ',a14,':',4f9.3,' d',i5,':',i2)
   9074 format(39x,'Total distance penalty: ',f10.3)
   9075 format(39x,'Total angle    penalty: ',f10.3)
   9076 format(39x,'Total torsion  penalty: ',f10.3)
   9077 format(/ /' Final Restraint Analysis for coords: ',a40/ /)
   9078 format(39x,'Total plane-point angle penalty: ',f10.3)
   9079 format(39x,'Total plane-plane angle  penalty: ',f10.3)
   9080 format(39x,'Total generalized distance coord.  penalty: ',f10.3)
   9082 format(' ',a14,' -- ',a14,':',4f9.3,' j')
   9083 format(' ',a14,' -- ',a14,':',4f9.3,' t')
   9084 format('| ',30x,'RMS deviation from ideal bonds : ',f11.4)
   9085 format('| ',30x,'RMS deviation from ideal angles: ',f10.3)   
   9086 format('Error: Atom groupings not supported for j-coupling restraints.')
   9087 format('Error: Improper restraint specified.')
   ! 9088 format('Error: Invalid residue number or atom name in restraint.')
   9089 format(' ',a14,' -- ',a14,':',4f9.3,' p')
   9090 format(' ',a14,' -- ',a14,':',4f9.3,' n')
   9091 format(' ',a14,' -- ',a14,':',4f9.3,' g')
end subroutine ndvprt 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmrcmf here]


subroutine nmrcmf(f,dfr,nmrat,nmrcom,rmstot,rimass,i)
   
   ! Subroutine NMR Center of Mass Forces
   
   ! This routine calculates the correct -d_E/d_position forces for the atoms
   ! which are involved in a distance restraint between two center of mass
   ! positions. The force is given by d(E)/d(x_com) * d(x_com)/d(x_atom).
   ! The DF array contains the d(E)/d(x_com) contributions, as calculated in
   ! routine DISNRG.
   
   ! Author: David A. Pearlman
   ! Date: 8/89
   
   ! Modified by Matthew Seetin 9/2007
   ! Bug fixes 12/2008 by Carlos Simmerling
   ! Comments updated 12/2008 by Matthew Seetin
   
   ! INPUT:
   ! F(I):     The force array.
   !           Modified on output to reflect the forces from this restraint.
   ! DFR(24):   d(E)/d(x_com),d(E)/d(y_com) and d(E)/d(z_com) for position one
   !           (elements 1->3), position two (elements 4->6), etc.
   ! NMRAT(16,I)
   ! NMRCOM(2,I): Define the group of atoms used to calculate the center of mass.
   !              If NMRAT(1,I) < 0, then
   !                 (NMRCOM(1,K) -> NMRCOM(2,K)), K=-NMRAT(1,I),-NMRAT(3,I)
   !              gives the atoms in the group whose center of mass is position 1.
   
   !              If NMRAT(2,I) < 0, then
   !                 (NMRCOM(1,K) -> NMRCOM(2,K)), K=-NMRAT(2,I),-NMRAT(4,I)
   !              gives the atoms in the group whose center of mass is position 2.
   ! RMSTOT(8): The total masses for all the atoms defining the group for each
   !            atom or COM group of the restraint.
    ! RIMASS(I): The inverse mass for atom I.
   
   ! I: The number of the restraint in the list.
   
   implicit none
   integer:: i, iat, ip, j, jatm, jx, m, nmrat, nmrcom
   _REAL_ :: dcomdx, dfr, f, rimass, rmass, rmstot
   dimension f(*),dfr(*),nmrat(16,*),nmrcom(2,*),rimass(*)
   dimension rmstot(8)
   
   ! Loop over all atoms or COM groups in the restraint:
   
   do iat = 1,8
      jatm = 3*(iat-1)
      
      ! Calculate the d(c.o.m.)/d(coordinate) derivatives, multiply them
      ! by the derivatives in DF, and add them into the force array.
      
      if (nmrat(iat,i) < 0 .and. nmrat(iat+8,i) < 0) then
      
         
         ! See the comments in the NMRGRP subroutine for more information about 
         ! the non-intuitive scheme for looping over all atoms in a COM below.
         ! NMRAT stores which "atoms" are actually COM groups and which entries
         ! in NMRCOM are part of the group.  NMRCOM stores the list of atoms, 
         ! but in a manner that stores a set of atoms that are consecutive 
         ! in atom number as a two entries in the 2D array: the start atom and
         ! the end atom.
      
      
         do ip = -nmrat(iat,i),-nmrat(8+iat,i)
            do j = nmrcom(1,ip),nmrcom(2,ip)
               jx = 3*(j-1)
               rmass = 1.0d0/rimass(j)
               dcomdx = rmass/rmstot(iat)
               f(jx+1) = f(jx+1) - dfr(jatm+1)*dcomdx
               f(jx+2) = f(jx+2) - dfr(jatm+2)*dcomdx
               f(jx+3) = f(jx+3) - dfr(jatm+3)*dcomdx
            end do
         end do
      else if (nmrat(iat,i) >= 0) then
         
         ! Standard derivatives if no group defined:
         
         jx = nmrat(iat,i)
         do m = 1,3
            f(jx+m) = f(jx+m) - dfr(jatm+m)
         end do
      end if
   end do
   
   return
end subroutine nmrcmf 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmrcms here]
subroutine nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
   
   
   ! Subroutine NMR Center of MaSs
   
   ! This routine calculates the center of mass of the group(s) defined for
   ! distance restraint I, and places them in the XCOM array.
   
   
   ! Author: David A. Pearlman
   ! Date: 8/89
   
   ! Note: This routine assumes that all atoms of the group are in the
   !       same periodic "box" (if applicable).
   
   implicit none
   integer:: i, iat, ifill, ip, j, jx, nmrat, nmrcom
   _REAL_ :: ricmmass, rimass, rmass, rmstot, x, xcom, xtot, ytot, ztot
   dimension x(*),xcom(24),nmrat(16,*),nmrcom(2,*),rimass(*),ricmmass(8)
   dimension rmstot(8)
   
   ! Loop over all atoms that make up the restraint.  See NMRGRP or NMRCMF
   ! FMI about how NMRAT and NMRCOM combine to store the atoms involved in
   ! the group.
   
   do iat = 1,8
      xtot = 0.0d0
      ytot = 0.0d0
      ztot = 0.0d0
      rmstot(iat) = 0.0d0
      ricmmass(iat) = 0.0d0
      ifill = 3*(iat-1)
      
      if (nmrat(iat,i) < 0) then
         do ip = -nmrat(iat,i),-nmrat(8+iat,i)
            do j = nmrcom(1,ip),nmrcom(2,ip)
               
               jx = 3*(j-1)
               rmass = 1.0d0/rimass(j)
               
               xtot = xtot + x(jx+1)*rmass
               ytot = ytot + x(jx+2)*rmass
               ztot = ztot + x(jx+3)*rmass
               rmstot(iat) = rmstot(iat) + rmass
               ricmmass(iat) = 1.0d0/rimass(iat)
            end do
         end do
         xcom(ifill+1) = xtot/rmstot(iat)
         xcom(ifill+2) = ytot/rmstot(iat)
         xcom(ifill+3) = ztot/rmstot(iat)
      else
         xcom(ifill+1) = x(nmrat(iat,i)+1)
         xcom(ifill+2) = x(nmrat(iat,i)+2)
         xcom(ifill+3) = x(nmrat(iat,i)+3)
      end if
   end do
   
   return
end subroutine nmrcms 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmrgrp here]
subroutine nmrgrp(i,kat,lastpk,iscrth,natom,nmrat,nmrcom,in,iout, &
      maxgrp,igr,grnam,iresid,name,ipres,nres,maxigr)
   
   
   ! Subroutine NMR GRouP
   ! This routine is called by NMRRED if the user is defining groups of atoms
   ! whose center-of-mass is to be used in any restraints. Here, the
   ! atom definitions are read in, sequentially sorted, and stored in the
   ! NMRCOM/NMRAT arrays.
   
   ! Rather than always tracking a complete list of all the atoms in the group,
   ! the groups are stored as ranges of atoms in the array NMRCOM(2,I).  The
   ! first atom of a continuous block of atoms involved in the same restraint I
   ! is stored in NMRCOM(1,I), and the last atom in the block is stored in
   ! NMRCOM(2,I).  There may be many of these such blocks if many non-sequential
   ! atoms are used to make up the group.
   
   ! For the group taking the place of IAT1 (KAT = 0):
   
   ! (NMRCOM(1,K) -> NMRCOM(2,K)), K=-NMRAT(1,I), -NMRAT(9,I)
   
   ! -----
   
   ! For the group taking the place of IAT2 (KAT = 1):
   
   ! (NMRCOM(1,K) -> NMRCOM(2,K)), K=-NMRAT(2,I), -NMRAT(10,I)
   
   ! and so on up through NMRAT(8,I) and NMRAT(16,I) if that many groups are defined
   
   ! If IRESID = 1: IGR(I) is the residue # for atom I, and GRNAM(I) is the
   !                corresponding atom name.
   !           = 0: IGR(I) is the atom number for atom I.
   
   ! NAME() and IPRES() are the name and residue pointers for the system;
   !     NRES is the number of residues in the system; all are used by GETNAT
   !     when IRESID = 1.
   
   ! MAXIGR is the dimension of IGR.
   ! LASTPK is the first unfilled location in NMRAT
   ! ISCRTH is a scratch array of >= NATOM length
   ! NATOM is the number of atoms
   
   
   ! Author: David A. Pearlman
   ! Date: 8/89
   
   ! Updated by Matthew Seetin to be compatible with the use of up to 8 atom groups
   ! Date: 9/2007
   
   
   implicit none
   integer:: i, ierr, igr, ii, in, iout, ipack, ipres, iresid, &
        iscrth, k, kat, lastpk, maxgrp, maxigr, natom, nmrat, nmrcom, nres
   !     CHARACTER*80 ALINE
   character(len=4) grnam(*),name(*)
   dimension nmrat(16,*),nmrcom(2,*),iscrth(*),igr(200),ipres(*)
   do k = 1,natom
      iscrth(k) = 0
   end do
   
   !  -- old formatted input:
   !  15 READ(IN,9001,END=1001,ERR=1005) IGR
   
   ! Loop over the defining atoms. If IRESID = 1, then convert the
   ! residue #/atom-name reference into an actual atom number pointer, and store
   ! this pointer back in IGR(I).
   
   do k = 1,maxigr
      if (igr(k) <= 0) goto 30
      
      ! resolve pointer, if necessary:
      
      if (iresid == 1) then
         call getnat(igr(k),grnam(k),name,ipres,nres,iout,ierr)
         if (ierr == 1) then
            write(iout,9088)
            call mexit(6, 1)
         end if
      end if
      
      if (igr(k) > natom) goto 1007
      iscrth(igr(k)) = 1
   end do
   
   ! -- old formatted input:
   !     GO TO 15
   
   30 nmrat(1+kat,i) = -lastpk
   ipack = 0
   do k = 1,natom
      if (iscrth(k) == 1) then
         if (ipack == 0) then
            nmrcom(1,lastpk) = k
            ipack = 1
         end if
      else if (iscrth(k) == 0) then
         if (ipack == 1) then
            nmrcom(2,lastpk) = k-1
            ipack = 0
            lastpk = lastpk + 1
            if (lastpk > maxgrp) goto 1009
         end if
      end if
   end do
   
   if (ipack == 1) then
      nmrcom(2,lastpk) = natom
      ipack = 0
      lastpk = lastpk + 1
   end if
   if (lastpk > maxgrp) goto 1009
   
   nmrat(9+kat,i) = -(lastpk-1)
   
   return
   
   ! End-of-file read errors:
   
   !1001 WRITE(IOUT,2000)
   !     call mexit(iout, 1)
   
   ! Other read errors (probably format mismatch):
   
   !1005 BACKSPACE(IN)
   !     READ(IN,9004) ALINE
   !     WRITE(IOUT,2001) ALINE
   !     call mexit(iout, 1)
   
   ! A group-defining atom number was out-of-range
   
   1007 write(iout,2003) k
   write(iout,9001) (igr(ii),ii=1,k)
   call mexit(6, 1)
   
   ! Too many atom ranges need to be stored, relative to storage allocated
   ! (MAXGRP).
   
   1009 write(iout,2005) maxgrp
   call mexit(6, 1)
   
   
   ! 2000 format(' Error: End of file read.')
   ! 2001 format(' Error reading the line:',/,a80)
   2003 format(' Error: Atom ',i2,' in following group definition is', &
         ' greater than total # atoms')
   2005 format(' Error: Too many atom ranges need to be stored for ', &
         'center-of-mass distance',/,t9,'restraints. ', &
         'MAXGRP =',i5,'. This needs to be increased.')
   9001 format(16i5)
   ! 9004 format(a80)
   9088 format('Error: Invalid residue number or atom name in restraint.')
end subroutine nmrgrp 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmrnrg here]
subroutine nmrnrg(x,f,rimass,ntb,nmrnum,nstep,nmrat,resttype,rstwtarr,&
      nmrst,r1nmr,r2nmr,r3nmr,r4nmr,rk2nmr,rk3nmr,tgtvec,nmrcom,nmrfty, &
      igravt,ialtdis,ifconstr,bave,bave0,aave,aave0, &
      tave1,tave01,tave2,ptave,ptave0,plave,plave0,gdave,gdave0, &
      tave02,ajcoef,nave,nexact,ipower,tauave,ravein, &
      dt,navint,iavtyp,idmpav,idumpu,enmr, &
      devdis,devang,devtor,devplpt,devpln,devgendis,iout)
   
   ! Subroutine NMR eNeRGy.
   
   ! This routine calculates the contributions
   ! to the energy/forces due to distance/valence angle/torsion restraints.
   
   ! Author: David A. Pearlman
   ! Date: 7/89
   ! Modified 9/2007 by Matthew Seetin to enable group restraints on angles, 
   ! torsions, and my new planar restraints
   
   ! INPUT:
   
   !  X(I)       : Coordinate array
   !  F(I)       : Force array; Updated during this call.
   !  RIMASS(I)  : Array of inverse masses
   !  NTB        : Periodic boundary conditions flag
   !  NMRNUM     : The number of NMR restraints
   !  NSTEP      : The step/iteration number in the main calling program
   !  NMRAT(16,I) : The 2-6 atoms defining restraint I.
   !                NMRAT(9,I)-NMRAT(16,I) are used for atom groups 
   !  NMRST(3,I) : Restraint I is imposed over the range of steps
   !               NMRST(1,I) -> NMRST(2,I). If NMRST(3,I) .NE. 0, change
   !               is implemented as a step function every ABS(NMRST(3,I)) steps.
   !               Otherwise, changes in target values occur continuously.
   !  R1NMR(2,I) : The slope/intercept of the dependence of r1 on step
   !               number (see routine DISNRG,ANGNRG or TORNRG for description
   !               of r1, r2, r3, r4, k2, and k3).
   !  R2NMR(2,I) : The slope/int. of the dependence of r2 on step number.
   !  R3NMR(2,I) : The slope/int. of the dependence of r3 on step number.
   !  R4NMR(2,I) : The slope/int. of the dependence of r4 on step number.
   !  RK2NMR(2,I) : The slope/int. of the dependence of k2 on step number.
   !  RK3NMR(2,I) : The slope/int. of the dependence of k3 on step number.
   !  NMRCOM(2,I) : The ranges of atoms defining the center-of-mass group,
   !                for applicable distance restraints
   !  NMRFTY(I)  : Flag for each restraint indicating what functional form should
   !               be used for that restraint. If time-averaged restraints
   !               are being applied to a particular class of internal
   !               (bond, angle, torsion) and NMRFTY(I)=1, then an instantaneous
   !               restraint will be applied to internal I.
   !  IGRAVT(I)  : For averaged group restraint positions.
   !               = 0, use center of mass
   !               = 1, use r**-6 averaged position.
   !  IALTDIS(I) : Flag to use alternative restraint functional form
   !    BAVE(I)
   !   BAVE0(I)
   !    AAVE(I)
   !   AAVE0(I)    For time-averaged restraints. _AVE(I) contains the time-averaged
   !   TAVE1(I)    value of the restraint, using the exponential decay value in
   !   TAVE01(I)   TAUAVE(I). _AVE0(I) contains the time-averaged value of the
   !   TAVE2(I)    restraint, using no exponential decay (real average). For
   !   TAVE02(I)   torsions, two values are required for average, corresponding
   !   PTAVE(I)    to the COS(tau) and SIN(tau)
   !  PTAVE0(I)
   !   PLAVE(I)
   !  PLAVE0(I)
   ! AJCOEF(3,I) : The coefficients to be used in the Karplus equation for
   !               any torsional restraints to be imposed using J-coupling.
   !     NAVE(5) : If time-averaged restraint have been requested, NAVE(I)>0.
   !               Elements 1->3 correspond to values for bonds, angles, torsions.
   !   NEXACT(5) : Not currently used.
   !   IPOWER(5) : For time-averaged restraints. Gives the exponent used in
   !               developing the time average.
   !   TAUAVE(5) : For time-averaged restraints. Gives the time constant used
   !               in the exponential decay multiplier for the energies/forces
   !               reported back to the calling program.
   !   RAVEIN(5) : For time-averaged restraints. On the first step, the value
   !               returned is the current value+RAVEIN(I) (I=1,3 for bonds,
   !               angles, torsions).
   !          DT : Integration timestep (ps). Used for time-averaged restraints.
   !      NAVINT : Time-averaged restraints are only updated every NAVINT
   !               steps. Typically, NAVINT=1.
   !   IAVTYP(5) : Determines how forces will be calculated when time-averaged
   !               retraints are used. IAVTYP(I) = 1 --> Calculate forces
   !               according to the standard energy expression. IAVTYP(I) = 2 -->
   !               calculate forces as dE/d(r_ave) dr(t)/d(x). Latter form
   !               integrates into a non-intuitive pseudo-energy, but avoids
   !               large forces from (1+IPOWER) exponentiation.
   !      IDMPAV : Values of retraints will be dumped to the file REDIR(6) every
   !               IDMPAV steps, if IDMPAV > 0. A time-averaged value for a
   !               restraint will be reported, if time-averaging is being performed
   !               for that retraint.
   !      IDUMPU : Unit to be used for restraint dumps if IDMPAV > 0. This
   !               unit must be dedicated for the entire run.
   
   
   ! OUTPUT:
   
   ! ENMR(1)     : Total energy contribution from distance restraints.
   ! ENMR(2)     : Total energy contribution from angle restraints.
   ! ENMR(3)     : Total energy contribution from torsion restraints.
   ! ENMR(4)     : Total energy contribution from plane-point angle restraints.
   ! ENMR(5)     : Total energy contribution from plane-plane angle restraints.
   ! ENMR(6)     : Total energy contribution from generalized distance coordinate restraints.
   ! DEVDIS(1): The average deviation of the distance restraints. The
   !            average of r2 and r3 is used as the "target" distance.
   !       (2): The rms deviation of the distance restraints. The average
   !            of r2 and r3 is used as the "target" distance.
   !       (3): Average deviation, but a deviation of 0 is used for restraints
   !            falling between r2 and r3, and r2 or r3 (whichever is
   !            closer) is used outside this range.
   !       (4): The rms deviation, but same definitions of deviations as (3).
   ! DEVANG(4): Same as DEVDIS, but for angles
   ! DEVTOR(4): Same as DEVDIS, but for torsions
   ! DEVPLPT(4): Same as DEVDIS, but for plane-point angles
   ! DEVPLN(4): Same as DEVDIS, but for plane-plane angles
   ! DEVPLN(4): Same as DEVDIS, but for generalized distance coordinate angles
   
   
   use constants, only : half, DEG_TO_RAD, RAD_TO_DEG 
   use file_io_dat
   implicit none
   integer:: i, ialtdis, iave, iavtyp, idmpav, idumpu, ifirst, &
        igravt, incflg, inum, ionce, iout, ipower, itimes, j, natom, &
        nave, navint, nexact, nmrat, nmrcom, nmrfty, nmrnum, nmrst, &
        nstep, nstepu, ntb
   _REAL_ :: aave, aave0, ajcoef, bave, bave0, convrt, dev, dev2, &
        devang, devdis, devgendis, deviat, devpln, devplpt, devtor, dfr, &
        dt, e, enmr, f, fcurr, fold, gdave, gdave0, plave, plave0, ptave, &
        ptave0, r1, r1nmr, r2, r2nmr, r3, r3nmr, r4, r4nmr, ravein, &
        rimass, rint, rk2, rk2nmr, rk3, rk3nmr, rmstot, rnum, rstwtarr, &
        small, step, target, tauave, tave01, tave02, tave1, tave2, work, &
        x, xcom
   logical jcoupl, usecom
#  include "nmr.h"
#  include "extra.h"
   dimension x(*),f(*),rimass(*),nmrat(16,*),rstwtarr(4,*),nmrst(3,*),r1nmr(2,*)
   dimension r2nmr(2,*),r3nmr(2,*),r4nmr(2,*),rk2nmr(2,*)
   dimension rk3nmr(2,*),nmrcom(2,*),nmrfty(*),igravt(*)
   dimension ialtdis(*)
   dimension devdis(4),devang(4),devtor(4),devplpt(4),devpln(4),devgendis(4),deviat(6,4),inum(6)
   dimension enmr(6),xcom(24),dfr(24),rmstot(8),incflg(6)
   dimension bave(*),bave0(*),aave(*),aave0(*),tave1(*)
   dimension tave01(*),tave2(*),tave02(*),ptave(*),ptave0(*),plave(*),plave0(*),gdave(*),gdave0(*),&
             ajcoef(3,*),nave(6)
   dimension nexact(6),ipower(6),tauave(6),ravein(6)
   dimension navint(6),iavtyp(6)
   ! Modified 9/2007 by Matthew Seetin to enable group restraints on angles, 
   ! torsions, and my new planar restraints
   _REAL_    tgtvec(*)
   integer   ifconstr(*)
   integer :: jjj, resttypetmp
   integer resttype(*)
   common/nmrnlc/ionce,itimes
   
   save small,work,fold
   data small/1.0d-7/

   _REAL_ ricmmass(8)
   _REAL_ comrimass,rstwt
   dimension comrimass(4),rstwt(4)
   
   do jjj=1,4
     comrimass(jjj) = 0.0
   end do
   
   ifirst = 0
   if (ionce /= 9899) then
      ifirst = 1
      ionce = 9899
   end if
   
   ! If incremental "dumps" of restraint values have been requested,
   ! and this is the first call to NMRNRG, open the appropriate file:
   
   if (master .and. ifirst == 1 .and. iredir(6) > 0) &
         call amopen(idumpu,redir(6)(1:iredir(6)),owrite,'F','W')
   
   ! Zero the accumulators:
   
   do i=1,6
      do j=1,4
         deviat(i,j) = 0.0d0
      end do
      inum(i) = 0
      incflg(i) = 1
      enmr(i) = 0.0d0
   end do
   ! Modified 9/2007 by Matthew Seetin to enable group restraints on angles, 
   ! torsions, and my new planar restraints
   
   ! Main loop over all the restraints:
   
   do i=1,nmrnum
!        write(*,*) 'restr num ',i
        
        if(morse.eq.1) then
        
         if(jar.eq.0.and.i.eq.1)  then
!         call morse(natoms,De,ro,nat1,nat2,X,F)    
          call smorse(NATOM,rk2nmr(2,i),r2nmr(2,i),nmrat(1,i), &
                     nmrat(2,i),X,F)                     
         goto 12      
         elseif(jar.eq.1.and.i.eq.2) then
!         write(*,*) 'Force befoe',F(nmrat(1,i)+1),F(nmrat(1,i)+2),F(nmrat(1,i)+3)
!         write(*,*) 'MAIN: ',matom,(nmrat(1,i)+3)/3,(nmrat(2,i)+3)/3
         call smorse(matom,rk2nmr(2,i),r2nmr(2,i),nmrat(1,i),   &
                     nmrat(2,i),X,F)
!         write(*,*) 'Force post',F(nmrat(1,i)+1),F(nmrat(1,i)+2),F(nmrat(1,i)+3)            
         goto 12
         endif
   
   ! fin loop sobre Morse
        endif

      ! Skip restraints if the current step # is outside the applicable range:
      
      rint = 0.0d0
      if (nstep < nmrst(1,i) .or. (nstep > nmrst(2,i) .and. &
            nmrst(2,i) > 0)) goto 12
      
      ! Calculate the values of r1,r2,r3,r4,k2, and k3 for this step & restraint:
      
      ! Vary the step used in calculating the values only every ABS(NMRST(3,I)) steps.
      ! If the user did not specify NINC in the input, ABS(NMRST(3,I))=1.
      
      nstepu = nstep - mod(nstep-nmrst(1,i),abs(nmrst(3,i)))
      step = nstepu
      r1 = r1nmr(1,i)*step + r1nmr(2,i)
      r2 = r2nmr(1,i)*step + r2nmr(2,i)
      r3 = r3nmr(1,i)*step + r3nmr(2,i)
      r4 = r4nmr(1,i)*step + r4nmr(2,i)
      
      ! If NMRST(3,I) > 0, then the weights are linearly interpolated.
      ! If NMRST(3,I) < 0, then the weights are modified by a multplicative factor.
      
      if (nmrst(3,i) > 0) then
         rk2 = rk2nmr(1,i)*step + rk2nmr(2,i)
         rk3 = rk3nmr(1,i)*step + rk3nmr(2,i)
      else
         nstepu = (nstep-nmrst(1,i))/abs(nmrst(3,i))
         rk2 = rk2nmr(2,i) * rk2nmr(1,i)**nstepu
         rk3 = rk3nmr(2,i) * rk3nmr(1,i)**nstepu
      end if
      
      ! Take the average of r2 & r3 as the average "target" value:
      
      target = (r2+r3)*half
      
      ! Determine what type of restraint (distance,angle,torsion,plain-point, plane-plane) this
      ! is, and call the appropriate routine to calculate its contribution.
      ! Modified 9/2007 by Matthew Seetin to enable group restraints on angles, 
      ! torsions, and my new planar restraints
      
      usecom = .false.
      do jjj=9,16
        usecom = usecom .or. nmrat(jjj,i) < 0
      end do
      
      
      iave = 0
      jcoupl = .false.
      resttypetmp = resttype(i)
      select case (resttypetmp)
        case (1) ! Distance restraint
         if (nave(1) > 0 .and. nmrfty(i) == 0) iave = iavtyp(1)
         if (.not. usecom) then
            call disnrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i),e,r1, &
                  r2,r3,r4,rk2,rk3,ntb,bave(2*i-1), &
                  bave0(i),nave(1),nexact(1),ipower(1), &
                  tauave(1),ravein(1),dt,navint(1),iave, &
                  incflg(1),ifirst,0,ialtdis(i))
            if(ifconstr(i) > 0) call vec_constr(x,f,dfr,rimass,dt,rint, &
                   nmrat(1,i),nmrat(2,i),tgtvec(3*(i-1)+1),0)
         else if (igravt(i) == 1) then
            call r6ave(x,nmrat,nmrcom,rint,i)
            call disnrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i),e,r1, &
                  r2,r3,r4,rk2,rk3,ntb,bave(2*i-1), &
                  bave0(i),nave(1),nexact(1),ipower(1), &
                  tauave(1),ravein(1),dt,navint(1),iave, &
                  incflg(1),ifirst,3,ialtdis(i))
            call r6drv(f,dfr(1))
         else
            call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
            call disnrg(xcom,f,dfr,rint,0,3,e,r1,r2,r3, &
                  r4,rk2,rk3,ntb,bave(2*i-1),bave0(i), &
                  nave(1),nexact(1),ipower(1),tauave(1), &
                  ravein(1),dt,navint(1),iave,incflg(1), &
                  ifirst,2,ialtdis(i))
            if(ifconstr(i) > 0) call vec_constr(xcom,f,dfr,ricmmass,dt,rint, &
                  0,3,tgtvec(3*(i-1)+1),1)
            call nmrcmf(f,dfr,nmrat,nmrcom,rmstot,rimass,i)
         end if
      case (2) ! Angle restraint
         if (nave(2) > 0 .and. nmrfty(i) == 0) iave = iavtyp(2)
         if (usecom) then
           call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
           call angnrg(xcom,f,dfr,rint,0,3,6, &
               e,r1,r2,r3,r4,rk2,rk3,ntb,aave(2*i-1), &
               aave0(i),nave(2),nexact(2),ipower(2), &
               tauave(2),ravein(2),dt,navint(2),iave, &
               incflg(2),ifirst,3)
           call nmrcmf(f,dfr,nmrat,nmrcom,rmstot,rimass,i)
         else
           call angnrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i),nmrat(3,i), &
               e,r1,r2,r3,r4,rk2,rk3,ntb,aave(2*i-1), &
               aave0(i),nave(2),nexact(2),ipower(2), &
               tauave(2),ravein(2),dt,navint(2),iave, &
               incflg(2),ifirst,0)
         end if
      case (3) ! Torsion or J-coupling
         if (nave(3) > 0 .and. nmrfty(i) == 0) iave = iavtyp(3)
         
         ! If J coefficients are 0 for this torsion, it is a normal torsional
         ! restraint. If they are .NE.0, then it is a J-coupling factor restraint.
         
         jcoupl = abs(ajcoef(1,i)) > small .or. &
               abs(ajcoef(2,i)) > small
         
         if (.not.jcoupl) then
           if (usecom) then
             call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
             call tornrg(xcom,f,dfr,rint,0,3,6,9, &
                  e,r1,r2,r3,r4,rk2,rk3,ntb, &
                  tave1(2*i-1),tave01(i),tave2(2*i-1), &
                  tave02(i),nave(3),nexact(3),ipower(3), &
                  tauave(3),ravein(3),dt,navint(3),iave, &
                  incflg(3),ifirst,3)
             call nmrcmf(f,dfr,nmrat,nmrcom,rmstot,rimass,i)
           else
             call tornrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i),nmrat(3,i), &
                  nmrat(4,i),e,r1,r2,r3,r4,rk2,rk3,ntb, &
                  tave1(2*i-1),tave01(i),tave2(2*i-1), &
                  tave02(i),nave(3),nexact(3),ipower(3), &
                  tauave(3),ravein(3),dt,navint(3),iave, &
                  incflg(3),ifirst,0)
           end if
         else
           if (usecom) then
             ! Since atom groups are not supported for j-coupling restraints, 
             ! print out an error message and terminate execution if the user 
             ! specifies atom groups AND j-coupling
             ! This exit condition shouldn't be used, since it should have been 
             ! checked in nmrred.  However, I have it here just for safety.
             
             write(iout,9074)
             call mexit(iout, 1)
           else
             call jnrg(x,f,rint,nmrat(1,i),nmrat(2,i),nmrat(3,i), &
                  nmrat(4,i),e,r1,r2,r3,r4,rk2,rk3,ntb, &
                  tave1(2*i-1),tave01(i),tave2(2*i-1), &
                  tave02(i),nave(3),nexact(3),ipower(3), &
                  tauave(3),ravein(3),dt,navint(3),iave, &
                  incflg(3),ifirst,ajcoef(1,i),0)
           end if
         end if
      case (4)  ! This is a plane-point angle restraint
        if (nave(4) > 0 .and. nmrfty(i) == 0) iave = iavtyp(4)
        if (usecom) then
          call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
          call commass(comrimass,nmrat,nmrcom,rimass,i)
          call plptnrg(xcom,f,dfr,comrimass,rint,0,3,6,9,12, &
               e,r1,r2,r3,r4,rk2,rk3,ntb,ptave(2*i-1), &
               ptave0(i),nave(4),nexact(4),ipower(4), &
               tauave(4),ravein(4),dt,navint(4),iave, &
               incflg(4),ifirst,2)
          call nmrcmf(f,dfr,nmrat,nmrcom,rmstot,rimass,i)
        else
          call plptnrg(x,f,dfr,rimass,rint,nmrat(1,i),nmrat(2,i), &
               nmrat(3,i),nmrat(4,i),nmrat(5,i), &
               e,r1,r2,r3,r4,rk2,rk3,ntb,ptave(2*i-1), &
               ptave0(i),nave(4),nexact(4),ipower(4), &
               tauave(4),ravein(4),dt,navint(4),iave, &
               incflg(4),ifirst,0)
        end if
      case (5) ! This is a plane-plane angle restraint
        if (nave(5) > 0 .and. nmrfty(i) == 0) iave = iavtyp(5)
        if (usecom) then
          call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
          call plnnrg(xcom,f,dfr,rint,0,3,6,9,12,15,18,21, &
               e,r1,r2,r3,r4,rk2,rk3,ntb,plave(2*i-1), &
               plave0(i),nave(5),nexact(5),ipower(5), &
               tauave(5),ravein(5),dt,navint(5),iave, &
               incflg(5),ifirst,2)
          call nmrcmf(f,dfr,nmrat,nmrcom,rmstot,rimass,i)
        else
          call plnnrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i), &
               nmrat(3,i),nmrat(4,i),nmrat(5,i),nmrat(6,i),nmrat(7,i),nmrat(8,i), &
               e,r1,r2,r3,r4,rk2,rk3,ntb,plave(2*i-1), &
               plave0(i),nave(5),nexact(5),ipower(5), &
               tauave(5),ravein(5),dt,navint(5),iave, &
               incflg(5),ifirst,0)
        end if
      case (6) ! This is a generalized distance coordinate restraint
        do jjj=1,4
          rstwt(jjj) = rstwtarr(jjj,i)
        end do
        if (nave(6) > 0 .and. nmrfty(i) == 0) iave = 1
        if (usecom) then
          call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
          call gendisnrg(xcom,f,dfr,rint,0,3,6,9,12,15,18,21, &
               rstwt,e,r1,r2,r3,r4,rk2,rk3,ntb,gdave(2*i-1), &
               gdave0(i),nave(6),nexact(6),ipower(6), &
               tauave(6),ravein(6),dt,navint(6),iave, &
               0,0,2)
          call nmrcmf(f,dfr,nmrat,nmrcom,rmstot,rimass,i)
        else
          call gendisnrg(x,f,dfr,rint,nmrat(1,i),nmrat(2,i), &
               nmrat(3,i),nmrat(4,i),nmrat(5,i),nmrat(6,i),nmrat(7,i),nmrat(8,i), &
               rstwt,e,r1,r2,r3,r4,rk2,rk3,ntb,gdave(2*i-1), &
               gdave0(i),nave(6),nexact(6),ipower(6), &
               tauave(6),ravein(6),dt,navint(6),iave, &
               0,0,0)
        end if
      case default
        ! Improper restraint specified.  This should have also been checked by nmrred,
        ! but again I have the exit condition here just in case.
        write(iout,9075)
        call mexit(iout, 1)
      end select
            

      enmr(resttype(i)) = enmr(resttype(i)) + e
      

      ! Sum the appropriate quantities into the deviation accumulators.
      ! DEVIAT(1,1->4) stores bonds; (2,1->4) stores angles; (3,1->4) stores torsions
      ! If this was a coupling constant restraint (JCOUPL=.true.), scale the
      ! deviation by 3.14/180. here, as it will be scaled by 180./3.14 before
      ! printing...
      
      dev = abs(target - rint)
      if (rint >= r2 .and. rint <= r3) then
         dev2 = 0.0d0
      else if (rint < r2) then
         dev2 = r2 - rint
      else
         dev2 = rint - r3
      end if
      if (jcoupl) dev  = dev  * DEG_TO_RAD
      if (jcoupl) dev2 = dev2 * DEG_TO_RAD
      
      deviat(resttype(i),1) = deviat(resttype(i),1) + dev
      deviat(resttype(i),2) = deviat(resttype(i),2) + dev**2
      deviat(resttype(i),3) = deviat(resttype(i),3) + dev2
      deviat(resttype(i),4) = deviat(resttype(i),4) + dev2**2
      inum(resttype(i)) = inum(resttype(i)) + 1
      
      ! Set INCFLG to zero after the first time-averaged restraint of a
      ! particular type has been stored.
      
      if (iave > 0) incflg(resttype(i)) = 0
      
      ! Modified 9/2007 by Matthew Seetin to enable group restraints on angles, 
      ! torsions, and my new planar restraints
      
      ! If requested, write out restraint values, all on a single line:
      12 continue

      ! Adding jar option to write correctly (roit. 02/27/05)
      ! if jar=0, the behaviour is as usual, if jar=1,
      ! writting is done in the same file, but in the way:
      ! r2(t), rint(t), force(t), accumulated work(t)

        
      if(jar == 0) then
        if (master .and. idmpav > 0) then
         if (mod(nstep,idmpav) == 0) then
            if (i == 1) write(idumpu,'(i8,1x)',advance='no') nstep

            convrt = 1.0d0
            if (resttype(i) /= 1 .and. resttype(i) /= 6 .and. .not.jcoupl) convrt = RAD_TO_DEG
            write(idumpu,'(f9.3,1x)',advance='no') rint*convrt
         end if
        end if

      elseif(jar.eq.1.and.i.eq.1) then

        if (master .and. idmpav > 0) then
            convrt = 1.0d0
            if (resttype(i) /= 1 .and. resttype(i) /= 6 .and. .not.jcoupl) convrt = RAD_TO_DEG

      ! Set jar variables to the initial values (roit. 02/27/05)

      if(ifirst == 1) then
        work = 0.0
        fold = 2.0*rk2*(rint-r2)
      endif

            ! calculates the jar force and work (roit. 02/27/05)
            fcurr = -2.0*rk2*(rint-r2)
            work = work + (fcurr + fold) * drjar/2.0
            fold = fcurr

         if (mod(nstep,idmpav) == 0) then
            write(idumpu,'(f12.5,2x)',advance='no') r2*convrt
            write(idumpu,'(f12.5,2x)',advance='no') rint*convrt
            write(idumpu,'(f12.5,2x)',advance='no') fcurr/convrt
            write(idumpu,'(f12.5,2x)',advance='no') work/convrt
         end if
        end if !master
      end if !jar

! end loop over all restraints
   end do
   
   ! End the dump record
   if (master .and. idmpav > 0 ) then
      if (mod(nstep,idmpav) == 0) write(idumpu,'(a)')
   end if
   
   ! Determine the RMS averages and return
   
   do jjj=1,6
      if (inum(jjj) > 0) then
         rnum = inum(jjj)
         deviat(jjj,1) = deviat(jjj,1)/rnum
         deviat(jjj,2) = sqrt(deviat(jjj,2)/rnum)
         deviat(jjj,3) = deviat(jjj,3)/rnum
         deviat(jjj,4) = sqrt(deviat(jjj,4)/rnum)
      end if
   end do
   
   devdis(1) = deviat(1,1)
   devdis(2) = deviat(1,2)
   devdis(3) = deviat(1,3)
   devdis(4) = deviat(1,4)
   devang(1) = deviat(2,1)
   devang(2) = deviat(2,2)
   devang(3) = deviat(2,3)
   devang(4) = deviat(2,4)
   devtor(1) = deviat(3,1)
   devtor(2) = deviat(3,2)
   devtor(3) = deviat(3,3)
   devtor(4) = deviat(3,4)
   devplpt(1) = deviat(4,1)
   devplpt(2) = deviat(4,2)
   devplpt(3) = deviat(4,3)
   devplpt(4) = deviat(4,4)
   devpln(1) = deviat(5,1)
   devpln(2) = deviat(5,2)
   devpln(3) = deviat(5,3)
   devpln(4) = deviat(5,4)
   devgendis(1) = deviat(6,1)
   devgendis(2) = deviat(6,2)
   devgendis(3) = deviat(6,3)
   devgendis(4) = deviat(6,4)
   ! need gen. coord. stuff here
   ! Modified 9/2007 by Matthew Seetin to enable group restraints on angles, 
   ! torsions, and my new planar restraints
   
   return
   ! 9002 format(9f8.3)
   9074 format('Error: Atom groupings not supported for j-coupling restraints.')
   9075 format('Error: Improper restraint specified.')
end subroutine nmrnrg 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Print NMR deviation energies to the human-readable output
subroutine nmrprt(dvdis,dvang,dvtor,dvplpt,dvpln,dvgendis,eenmr,nstep,iout)
   
   
   ! Subroutine NMR PRinT
   
   ! This routine prints the energies/deviations due to the NMR restraints
   ! to the user-formatted file.
   
   ! Author: David A. Pearlman
   ! Date: 7/89
   
   ! DVDIS,DVANG,DVTOR,DVPLPT,DVPLN,DVGENDIS,EENMR(1,I) are the values on the last call to NMRNRG.
   ! DVDIS,DVANG,DVTOR,DVPLPT,DVPLN,DVGENDIS,EENMR(2,I) are the accumulated totals over the entire run.
   
   implicit none
   integer:: iout, j, nstep
   _REAL_ :: dvang, dvdis, dvgendis, dvpln, dvplpt, dvtor, eenmr, rstepu
#  include "nmr.h"
   dimension dvdis(2,4),dvang(2,4),dvtor(2,4),dvplpt(2,4),dvpln(2,4),dvgendis(2,4)
   dimension eenmr(2,6)
   logical, save :: printsecondline, printthirdline
   data printsecondline/.false./, printthirdline/.false./
   integer i
   
   ! The first solution to the issue of not printing the second line of restraint
   ! energies involved not printing them if the restraint energies were 0.  However,
   ! if a flat-bottomed well is used, the energies may very well be zero for some 
   ! steps even when a planar restraint is active.  I think it's important to 
   ! continue to print out the energies for those restraints in that case.  Checking
   ! the restraint type is a cleaner way to ensure that the energies for any active 
   ! restraint are used.  This checks what lines will be printed only on the first
   ! call to this subroutine.
   
   if (nmropt < 1) return
   if (nstep < 0) return
   
   if (.not. printsecondline) printsecondline = printsecondline .or. (eenmr(1,4) > 0.0d0 .or. &
                                                                      eenmr(1,5) > 0.0d0)
   if (.not. printthirdline) printthirdline = printthirdline .or. (eenmr(1,6) > 0.0d0)
   
   rstepu = max(nstep,1)
   write(iout,20) (eenmr(1,j),j=1,3)
   if (printsecondline) write(iout,21) (eenmr(1,j),j=4,5)
   if (printthirdline) write(iout,22) eenmr(1,6)       ! This will expand when we make gen. ang. 
                                                       ! and gen. tor. restraints.
   if (enoe > 0.0 .or. eshf > 0.0 .or. epcshf > 0.0) &
         write(iout,23) enoe,eshf,epcshf
   if( ealign > 0.d0 .or. ecsa > 0.d0 ) then
      write(iout,41) ealign,ecsa
      do i=1,num_datasets
         write(iout,42) s11(i), s12(i), s13(i)
         write(iout,43) s12(i), s22(i), s23(i)
         write(iout,43) s13(i), s23(i), -s11(i)-s22(i)
      end do
   end if
   write(iout,40)
   return
   
   20 format(' NMR restraints: Bond =',f9.3,3x,'Angle = ',f9.3,3x,'Torsion = ',f9.3)
   21 format('               : Plane-Point = ',f9.3,3x,'Plane-Plane = ',f9.3)
   22 format('               : Gen. Dist. Coord. = ',f9.3)
   23 format('               : Noesy=',f9.3,3x,'Shift = ',f9.3,3x, &
         'Pcshift = ',f9.3)
   40 format(79('='))
   41 format(' Energy (this step): Align=',f10.3, '  CSA=', f10.3)
   42 format('          Alignment tensor:',3f10.3)
   43 format('                           ',3f10.3)
end subroutine nmrprt 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmrrad here]
subroutine nmrrad(rad,wel,cn1,cn2,ntypes,ifill,rwell)
   
   
   ! Subroutine NMR RADius and well
   
   ! If IFILL = 0,
   ! this routine recalculates the values of the equilibrium radii and
   ! well depths if the A/B vdw coefficients have been changed.
   
   ! If IFILL = 1, this routine sets the well depths of all CN1,CN2 coefficients
   ! to RWELL, and recalculates the values of the equilibrium radii and
   ! well depths to reflect the A/B coefficients.
   
   ! Author: David A. Pearlman
   ! Date: 8/89
   
   ! INPUT:
   ! CN1(I) : Repulsive coeffient array
   ! CN2(I) : Attractive coefficient array
   ! NTYPES : The number of atom types
   ! IFILL  : Determines the action of this routine, as described above.
   ! RWELL  : If IFILL=1, the well depth for each pair of CN1,CN2 coefficients
   !          is set to RWELL.
   
   ! OUTPUT
   ! RAD(I): Equilibrium radius for atom-type I.
   ! WEL(I): Well depth for atom-type I.
   ! CN1(I)
   ! CN2(I): Modified to reflect RWELL, if IFILL=1.
   
   implicit  none
   _REAL_    cn1, cn2, rad, wel, rwell
   dimension cn1(*), cn2(*), rad(*), wel(*)
   integer   ntypes, ifill
   _REAL_    zero, two, half, fourth, sixth
   parameter (zero = 0.0d0, two = 2.0d0)
   parameter (half = 0.5d0, fourth = 0.25d0, sixth = 1.0d0/6.0d0)
   _REAL_    small
   parameter (small = 1.0d-7)
   _REAL_    a, b, welchg
   integer   i, ipt
   
   if (ifill == 0) then
      do i = 1,ntypes
         ipt = i*(i+1)/2
         a = cn1(ipt)
         b = cn2(ipt)
         rad(i) = zero
         wel(i) = zero
         if (a > small .and. b > small) then
            wel(i) = fourth * b*b/a
            rad(i) = half * (two * a/b)**sixth
         end if
      end do
   else
      do i=1,ntypes*(ntypes+1)/2
         a = cn1(i)
         b = cn2(i)
         if (a > small .and. b > small) then
            welchg = rwell/(fourth * b*b/a)
            cn1(i) = cn1(i)*welchg
            cn2(i) = cn2(i)*welchg
         end if
      end do
      
      do i=1,ntypes
         ipt = i*(i+1)/2
         a = cn1(ipt)
         b = cn2(ipt)
         rad(i) = zero
         wel(i) = zero
         if (a > small .and. b > small) then
            wel(i) = rwell
            rad(i) = half * (two * a/b)**sixth
         end if
      end do
   end if  ! (ifill == 0)
   
   return
end subroutine nmrrad 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ NMR REaD; read wt namelist containing the weight change info.
subroutine nmrred(x,name,ipres,rimass,r1nmr,r2nmr,r3nmr,r4nmr, &
      rk2nmr,rk3nmr,ajcoef,tgtvec,nmrat,resttype,rstwtarr,nmrst,nmrcom,nmrshb, &
      nmrfty,igravt,ialtdis,nmrnum,wtnrg,iwtstp, &
      iwttyp,kxpk,ifconstr, &
      iscrth,ichgwt,ishrtb,natom,nres,nstep0,stpmlt, &
      in,iout,iscop1,iscop2,maxrst,maxwt,maxgrp, &
      nave,nexact,navint,ipower,ravein,taufin, &
      iavtyp,idmpav,itimav)
   
   
   ! Subroutine NMR REaD.
   
   ! This subroutine reads 1) The NMR restraints;
   !                       2) Flags describing changes to the relative
   !                          weights of the energy terms.
   
   ! Author: David A. Pearlman
   ! Date: 7/89
   ! Modified by Matthew Seetin
   ! Date: 9/2007
   
   
   
   ! Information is read in three sections.
   ! SECTION 1 -- "WEIGHT CHANGE" cards:
   
   ! Information can either be read formatted, in the format described below,
   ! or as part of the "wt" namelist. To signal that formatted input is
   ! desired, you must preceed the weight-change card section by the flag
   ! "&formwt". This flag can be in any column, but must be the first not
   ! blank set of characters on the line.
   
   ! Namelist input follows standard namelist rules...
   
   ! Control lines of the format:
   !        TYPE , ISTEP1 , ISTEP2 , VALUE1 , VALUE2 , IINC , IMULT
   !            FORMAT(A8,2I5,2F12.6,2I5)
   
   
   ! -----------------------------------------------------------------------
   ! After "END" is read in the previous section,
   !    input/output redirection information can be read as described here.
   !    The inclusion of cards for redirection is optional.
   
   !    The format of the redirection cards is
   !         TYPE = filename
   !    where TYPE is any valid redirection keyword (see below), and filename
   !    is any character string. The equals sign ("=") is required.
   ! -----------------------------------------------------------------------------
   
   ! The input/output redirection cards (if any) are followed by the restraints.
   
   ! -----------------------------------------------------------------------------
   
   ! INPUT VARIABLES:
   ! ---------------
   
   !        X(I) : Coordinate Array
   !     NAME(I) : The name corresponding to atom I.
   !    IPRES(I) : Residue pointer array. IPRES(I) points to the first atom
   !               of residue I.
   !   RIMASS(I) : The inverse mass of atom I.
   !          IN : Unit from which the restraints/weight changes will be read.
   !        IOUT : Unit for informational prints.
   !      ISCOP1
   !      ISCOP2 : Scratch units for informational reads/writes. Used when
   !               redirection of input/output is requested.
   !       NATOM : The number of atoms.
   !        NRES : The number of residues.
   !      MAXRST : The maximum number of restraints for which array space is avail.
   !       MAXWT : The maximum number of weight changes
   !      MAXGRP : The maximum number of atom range definitions relating to
   !               group definitions (for center-of-mass distance restraints only)
   !               for which array space is available.
   
   ! OUTPUT:
   
   ! (re: restraints)
   !  R1NMR(2,I) : The slope and intercept, respectively, defining the dependency
   !               of r1 on step number for restraint I.
   !  R2NMR(2,I) : The slope and intercept for r2.
   !  R3NMR(2,I) : The slope and intercept for r3.
   !  R4NMR(2,I) : The slope and intercept for r4.
   ! RK2NMR(2,I) : The slope and intercept for k2.
   ! RK3NMR(2,I) : The slope and intercept for k3.
   ! AJCOEF(3,I) : For torsional restraints, if RJCOEF(3,I).NE.0, then
   !               a J-coupling restraint based on the Karplus relationship
   !               will be used.
   !  NMRAT(16,I) : The atoms defining restraint I. Stored as (3*(IAT-1)), where
   !               IAT is the absolute atom number.
   !  RESTTYPE(I) :  An integer identifying the type of restraint
   !                           = 1, this is a distance restraint
   !                           = 2, this is an angle restraint
   !                           = 3, this is a torsion or j-coupling restraint
   !                           = 4, this is a plane-point restraint
   !                           = 5, this is a plane-plane restraint
   !                           = 6, this is a generalized distance coordinate restraint
   !                               There are plans to implement generalized angle and torsion 
   !                               restraints, but not for this version of AMBER
   !  NMRST(3,I) : The beginning (1) and ending (2) steps over which this
   !               restraint is to be applied. NMRST(3,I) gives the increment
   !               between changes in the target values. If NMRST(3,I)=0,
   !               the target values are varied continuously.
   ! AJCOEF(3,I) : The coefficients to be used in the Karplus equation for
   !               any torsional restraints to be imposed using J-coupling.
   !   NMRSHB(I) : Flag for each restraint, which keeps track of whether that
   !               restraint is part of a defined "short-range" interactions
   !               group. NMRSHB(I) = 0 (not short-range or short range not
   !               defined); = 1 (short-range based on residue sequence);
   !               =2 (short-range based on distance criteria); =3
   !               (short-range based on both residue sequence and distance
   !               criteria).
   !   NMRFTY(I) : Flag for each restraint which indicates what functional
   !               form is used for the restraint. See description of IFNTYP
   !               above.
   !      NMRNUM : The total number of NMR restraints defined.
   
   ! (re: weight/temperature changes)
   !  WTNRG(2,I) : The slope and intercept, respectively, defining the dependency
   !               of weight/temperature/etc. change I on step number.
   ! IWTSTP(3,I) : The beginning and ending steps over which change I is to be
   !               applied (elements 1 and 2). If IWTSTP(3,I)>0, the change
   !               in target values is done as a step function every IWTSTP(3,I)
   !               steps. Otherwise, the change is done continuously.
   !   IWTTYP(I) : The type(s) of variable(s) being changed.
   !        KXPK : Peak identifier
   !      ICHGWT : Number of change parameters read.
   !      ISHRTB : If >0, then short-range interactions have been defined
   !               (using the SHORT command card). The definitions of the
   !               interactions are stored in IWTSTP( ,ISHRTB) and
   !               WTNRG( ,ISHRTB)>
   
   ! (re: other)
   !      NSTEP0 : If an NSTEP0 card is read, the initial value of the NSTEP0
   !               counter will be set to the given value. Otherwise, NSTEP0
   !               is returned with a value of 0.
   !      STPMLT : If a STPMLT card is read, the step counter multiplier will
   !               be set to the given value. Otherwise, STPMLT is returned
   !               with a value of 1.0.
   !     NAVE(5) : For time-averaged restraints. If time averaged restraints
   !               have been requested, NAVE(I) > 0.
   !               Elements 1->3 correspond to bonds, angles, torsions.
   !   NEXACT(5) : Not currently used.
   !   NAVINT(5) : If > 0, then averaged restraints are updated every
   !               NAVINT(I) steps.
   !   IPOWER(5) : For time-averaged restraints. Gives the exponent used in
   !               developing the time average.
   !   RAVEIN(5) : For time-averaged values, the value for the integral
   !               on the first step will be the instantaneous value+RAVEIN(I)
   !               (I=1,3 for bonds, angles, torsions).
   !   TAUFIN(5) : The value of the exponential decay constant used in calculating
   !               the final time-averaged values reported by NDVPRT.
   !   IAVTYP(5) : Determines how forces will be calculated when time-averaged
   !               retraints are used. IAVTYP(I) = 1 --> Calculate forces
   !               according to the standard energy expression. IAVTYP(I) = 2 -->
   !               calculate forces as dE/d(r_ave) dr(t)/d(x).
   !   IDMPAV    : If > 0, then "dump" the values of the restraints to a file
   !               every IDMPAV steps. The dump file is specified by a
   !               DUMPAVE redirection card.
   !   ITIMAV    : If = 0, then requests for time-averaged restraints will be
   !               ignored.

   use constants, only : zero, one, RAD_TO_DEG
   use file_io_dat
#ifdef MPI /* SOFT CORE */
  use softcore, only:  ifsc
#endif
   
   implicit none
   _REAL_  ajcoef
#ifdef MPI /* SOFT CORE */
   _REAL_ oneweight
#endif
   integer ialtdis
   integer iavtyp
   integer ichgwt
   integer idmpav
   integer igravt
   integer in
   integer iout
   integer ipower
   integer ipres
   integer iscop1
   integer iscop2
   integer iscrth
   integer ishrtb
   integer itimav
   integer iwtstp
   integer iwttyp
   integer kxpk
   integer maxgrp
   integer maxrst
   integer maxwt
   character(len=4) name(1)
   integer natom
   integer nave
   integer navint
   integer nexact
   integer nmrat
   integer resttype, resttypetmp
   integer nmrcom
   integer nmrfty
   integer nmrnum
   integer nmrshb
   integer nmrst
   integer nres
   integer nstep0
   integer ifconstr(1)       ! cuigl
   _REAL_  r1nmr
   _REAL_  r2nmr
   _REAL_  r3nmr
   _REAL_  r4nmr
   _REAL_  ravein
   _REAL_  rimass
   _REAL_  rk2nmr
   _REAL_  rk3nmr
   _REAL_  stpmlt
   _REAL_  taufin
   _REAL_  wtnrg
   _REAL_  x
   _REAL_  tgtvec(1)         ! cuigl
   _REAL_  rstwt, rstwtpass
   _REAL_  rstwtarr
   _REAL_  rstwttol
   parameter (rstwttol = 1.0d-7)

   character(len=80) aline
   character(256) restraint
   ! Added 9/2007 by Matthew Seetin to enable natural language input of restraints
   character(len=4) atnam(8),grnam1,grnam2,grnam3,grnam4,grnam5,grnam6,grnam7,grnam8, &
                    grnamarr, grnampass
   character(4) nametmp(8)
   ! Added 9/2007 by Matthew Seetin to enable group restraints on angles, 
   ! torsions, and my new planar restraints
   logical jcoupl
   logical usecom
   ! Added 9/2007 by Matthew Seetin to enable group restraints on angles, 
   ! torsions, and my new planar restraints
   
   ! ... IFLAG is the number of "TYPE" flags recognized (make sure "END" is last)
   
   integer iflag
#ifdef LES
   parameter (iflag = 32)
#  include "les.h"
#else
   parameter (iflag = 31)
#endif
   ! ... ZERNER is a value near zero used for weights set to 0.
   _REAL_ zerner
   parameter (zerner = 1.0d-7)
   character(len=8) type,flag(iflag)
   character(len=10) redirc(9)
#  include "nmr.h"
   ! Added to include 'nstlim' variable (roit. 02/27/05)
#  include "md.h"

   
   dimension x(*),ipres(*),rimass(*)
   dimension r1nmr(2,*),r2nmr(2,*),r3nmr(2,*),r4nmr(2,*),rk2nmr(2,*)
   dimension rk3nmr(2,*),nmrat(16,*),resttype(*),rstwt(4),rstwtpass(4),&
             rstwtarr(4,*),nmrst(3,*),nmrcom(2,*),nmrshb(*)
   ! Modified 9/2007 by Matthew Seetin to enable group restraints on angles, 
   ! torsions, and my new planar restraints
   dimension nmrfty(*),igravt(*),wtnrg(2,*),iwtstp(3,*),iwttyp(*)
   dimension ialtdis(*)
   dimension ajcoef(3,*),kxpk(2,*)
   dimension iscrth(*),xcom(24),dfr(24),rmstot(8),wtls(2)
   dimension ipower(6),nave(6),nexact(6),navint(6),ravein(6)
   dimension taufin(6),iavtyp(6),rjcoef(3)
   dimension iat(8)
   ! Modified 9/2007 by Matthew Seetin to enable new planar restraints
   
   ! ...  MAXIGR is the maximum number of atoms which can be used to define
   !      a center-of-mass group.
   
   integer ixpk,nxpk
   integer    maxigr
   parameter (maxigr=200)
   
   
   _REAL_  apmean
   _REAL_  convrt
   _REAL_  dev2
   _REAL_  dfr
   _REAL_  dumm
   _REAL_  dumma(6)
   _REAL_  e
   _REAL_  f(1)
   integer i
   integer ialtd
   integer iat
   integer iat1
   integer iat2
   integer iat3
   integer iat4
   integer iat5
   integer iat6
   integer iat7
   integer iat8
   ! Added 9/2007 by Matthew Seetin to enable new planar restraints
   integer ibeg
   integer idumm(6)
   integer ierr
   integer ifind
   integer ifntyp
   integer iformw
   integer ifvari
   integer iconstr     ! cuigl
   integer igr1
   integer igr2
   integer igr3
   integer igr4
   integer igr5
   integer igr6
   integer igr7
   integer igr8, igrpass
   integer igrarr
   ! Added 9/2007 by Matthew Seetin to enable group restraints on angles, 
   ! torsions, and my new planar restraints
   integer ihits
   integer ii
   integer iin
   integer iinc
   integer ilengt
   integer imat
   integer imult
   integer ir6
   integer iresid
   integer irstyp
   integer istep1
   integer istep2
   integer istrt
   integer j, jjj, iii
   integer jdx
   integer kk
   integer lastpk
   integer m
   integer n
   integer ninc
   integer nstep1
   integer nstep2
   integer ntu
   _REAL_  r1
   _REAL_  r1a
   _REAL_  r2
   _REAL_  r2a
   _REAL_  r3
   _REAL_  r3a
   _REAL_  r4
   _REAL_  r4a
   _REAL_  r0
   _REAL_  k0
   _REAL_  r0a
   _REAL_  k0a
   _REAL_  rcurr
   _REAL_  rint
   _REAL_  rjcoef
   _REAL_  rk2
   _REAL_  rk2a
   _REAL_  rk3
   _REAL_  rk3a
   _REAL_  rmstot
   _REAL_  rntu
   _REAL_  rprt
   _REAL_  slope
   _REAL_  step1
   _REAL_  step2
   _REAL_  value1
   _REAL_  value2
   _REAL_  wtls
   _REAL_  xcom
   _REAL_  ricmmass(8)   ! the inverse of total mass of a group -- cuigl
   
   _REAL_ comrimass
   dimension comrimass(4)
   
   dimension igr1(maxigr),igr2(maxigr),grnam1(maxigr),grnam2(maxigr), &
             igr3(maxigr),igr4(maxigr),grnam3(maxigr),grnam4(maxigr), &
             igr5(maxigr),igr6(maxigr),grnam5(maxigr),grnam6(maxigr), &
             igr7(maxigr),igr8(maxigr),grnam7(maxigr),grnam8(maxigr), &
             igrpass(maxigr), grnampass(maxigr) 
   dimension igrarr(maxigr,8), grnamarr(maxigr,8)
   ! Added 9/2007 by Matthew Seetin to enable group restraints on angles, 
   ! torsions, and my new planar restraints
   
   

   ! ----------------------------------------------------------------------------
   ! ----------------------------------------------------------------------------
   ! The following is all related to the use of NAMELIST reads.
   
   equivalence (iat(1),iat1),(iat(2),iat2)
   equivalence (iat(3),iat3),(iat(4),iat4)
   equivalence (iat(5),iat5),(iat(6),iat6)
   equivalence (iat(7),iat7),(iat(8),iat8)

   ! Added 9/2007 by Matthew Seetin to enable new planar restraints
   
   namelist /wt/ type,istep1,istep2,value1,value2,iinc,imult
   namelist /rst/ iat,restraint,rstwt,nstep1,nstep2,irstyp,ninc, &
         iresid,imult,atnam,igr1,igr2,igr3,igr4,igr5,igr6,igr7,igr8, &
         grnam1,grnam2,grnam3,grnam4,grnam5,grnam6,grnam7,grnam8,&
         r0,r1,r2,r3,r4,k0,rk2,rk3,r0a,r1a,r2a,r3a,r4a,k0a,rk2a,rk3a,ir6,ifntyp, &
         ifvari,rjcoef,ixpk,nxpk,ialtd,iconstr
   ! ----------------------------------------------------------------------------
   
   data flag/'BOND    ' , 'ANGLE   ' , 'TORSION ' , 'VDW     ', &
         'HB      ' , 'ELEC    ' , 'NB      ' , 'ATTRACT ', &
         'REPULSE ' , 'INTERN  ' , 'ALL     ' , 'REST    ', &
         'TEMP0   ' , 'RSTAR   ' , 'IMPROP  ' , 'SOFTR   ', &
         'CUT     ' , 'SHORT   ' , 'RESTS   ' , 'RESTL   ', &
         'NOESY   ' , 'SHIFTS  ' , 'TAUTP   ' , 'DISAVE  ', &
         'ANGAVE  ' , 'TORAVE  ' , 'TGTRMSD ' , 'PLPTAVE ', &
         'PLNAVE  ' , 'GDISAVE ' , &
#ifdef LES
         'TEMP0LES' , &
#endif
         'END     '/

   data redirc/'LISTIN    ' , 'LISTOUT   ' , 'DISANG    ', &
         'NOESY     ' , 'SHIFTS    ' , 'DUMPAVE   ', &
         'PCSHIFT   ' , 'DIPOLE    ' , 'CSA       ' /
   
   nstep0 = 0
   stpmlt = ONE
   lastpk = 1
   ishrtb = 0
   idmpav = 0
   ibeg = 1
   iat1 = 0
   wtls(1) = ONE
   wtls(2) = ONE
   ravein(1) = ZERO
   ravein(2) = ZERO
   ravein(3) = ZERO
   ravein(4) = ZERO
   ravein(5) = ZERO
   ravein(6) = ZERO
   ihits = 0
   
   do i=1,6
      iredir(i) = 0
      nave(i) = 0
      nexact(i) = 0
      taufin(i) = 1.0d+7
      navint(i) = 1
      iavtyp(i) = 1
   end do
   iredir(6) = 0
   iredir(7) = 0
   iredir(8) = 0
   iredir(9) = 0
   
   do jjj=1,4
     comrimass(jjj) = 0.0
   end do
   
   do jjj = 1,maxrst
     resttype(jjj) = 0
   end do
   
   write(iout,9000)
   
   ! Look for a card starting with the string "&formwt". If one is found, then
   ! the subsequent cards, up until one with TYPE="END" is encountered, are
   ! the weight cards, as formatted input. If this string is not found,
   ! assume input will be in namelist format.
   
   iformw = 0
   call nmlsrc('formwt',in,ifind)
   if (ifind == 1) then
      read(in,*)
      iformw = 1
   end if
   
   ! Read the weight changes
   
   9 do i = ibeg,99999
      
      ! Branch depending on whether formatted weight changes (IFORMW = 1) or
      ! namelist weight changes are being read:
      
      11 if (iformw /= 1) then
         
         ! NAMELIST-type read:
         ! ==================
         
         ! Set the default values of the parameters which can be set by the namelist
         
         istep1 = 0
         istep2 = 0
         value1 = ZERO
         value2 = ZERO
         iinc = 0
         imult = 0
         
         ! Look for "wt" namelist. If not found, stop with error.
         
         call nmlsrc('wt',in,ifind)
         if (ifind == 0) goto 1002
         
         read(in,nml=wt)
         !------------------------------------------------------------------
         
         ! Formatted-type weight card read:
         ! ==============================
         
      else
         read(in,9005,end=1001,err=1005) type,istep1,istep2, &
               value1,value2,iinc,imult
      end if
      
      ! End of weight card read. Now process given options:
      ! ==================================================
      
      ! ITYPE = "#" is a spacer card
      
      if (type(1:1) == '#') then
         write(iout,9006) type,istep1,istep2,value1,value2,iinc, &
               imult
         goto 11
      end if
      
      ! If ITYPE="END", then go to redirection reading section:
      
      if (type == flag(iflag) .or. type(1:3) == 'end') then
         if (i == 1) write(iout,9001)
         write(iout,9002)
         ichgwt = i-1
         goto 20
      end if
      
      ! If ITYPE="NSTEP0", then reset the initial value of the step counter
      ! and go to the next restraint.
      
      if (type == 'NSTEP0') then
         write(iout,9006) type,istep1,istep2,value1,value2,iinc, &
               imult
         nstep0 = istep1
         goto 11
      end if
      
      ! If ITYPE="STPMLT", then reset the value of the step counter multiplier.
      ! By default, this multiplier has a value of 1.0.
      
      if (type == 'STPMLT') then
         write(iout,9006) type,istep1,istep2,value1,value2,iinc, &
               imult
         if (value1 > zerner) stpmlt = value1
         goto 11
      end if
      
      ! If ITYPE="DISAVI", "ANGAVI", or "TORAVI", then set various values related
      ! to the time-integral. Also set the frequency of "dumping" to the restraint
      ! value dump file.
      
      select case (type)
       case('DISAVI') 
         navint(1) = 1
         if (idmpav == 0) idmpav = istep2
         ravein(1) = value1
         if (value2 > zerner) taufin(1) = value2
         if (iinc == 1) iavtyp(1) = 2
         write(iout,9006) type,istep1,istep2,value1,value2,iinc,imult
         goto 11
       case('ANGAVI')
         navint(2) = 1
         if (idmpav == 0) idmpav = istep2
         ravein(2) = value1
         if (value2 > zerner) taufin(2) = value2
         if (iinc == 1) iavtyp(2) = 2
         write(iout,9006) type,istep1,istep2,value1,value2,iinc,imult
         goto 11
       case ('TORAVI')
         navint(3) = 1
         if (idmpav == 0) idmpav = istep2
         ravein(3) = value1
         if (value2 > zerner) taufin(3) = value2
         if (iinc == 1) iavtyp(3) = 2
         write(iout,9006) type,istep1,istep2,value1,value2,iinc,imult
         goto 11
       case ('PLPTAVI')
         navint(4) = 1
         if (idmpav == 0) idmpav = istep2
         ravein(4) = value1
         if (value2 > zerner) taufin(4) = value2
         if (iinc == 1) iavtyp(4) = 2
         write(iout,9006) type,istep1,istep2,value1,value2,iinc,imult
         goto 11
       case ('PLNAVI') 
         navint(5) = 1
         if (idmpav == 0) idmpav = istep2
         ravein(5) = value1
         if (value2 > zerner) taufin(5) = value2
         if (iinc == 1) iavtyp(5) = 2
         write(iout,9006) type,istep1,istep2,value1,value2,iinc,imult
         goto 11
       case ('GDISAVI')
         navint(6) = 1
         if (idmpav == 0) idmpav = istep2
         ravein(6) = value1
         if (value2 > zerner) taufin(6) = value2
         if (iinc == 1) iavtyp(6) = 2
         write(iout,9006) type,istep1,istep2,value1,value2,iinc,imult
         goto 11
       case ('DUMPFREQ')
         idmpav = istep1
         write(iout,9006) type,istep1,istep2,value1,value2,iinc,imult
         goto 11
      end select
      
      ! Check to make sure maximum allowable number of weight cards not exceeded
      
      if (i > maxwt) then
         write(iout,9008) maxwt
         call mexit(iout, 1)
      end if
      
      
      ! Loop through valid TYPEs
      
      do j = 1,iflag
         if (type == flag(j)) goto 14
      end do
      write(iout,9010)
      write(iout,9006) type,istep1,istep2,value1,value2,iinc,imult
      call mexit(iout, 1)
      
      ! The TYPE is valid. Store an integer pointer to the type of flag
      ! in IWTTYP(I). Determine the slope and intercept of changes in VALUE
      ! with step number, and store these in WTNRG. Store the step range
      ! over which this change is valid in IWTSTP.
      
      14 write(iout,9006) type,istep1,istep2,value1,value2,iinc,imult
      iwttyp(i) = j
      if (istep1 > istep2) write(iout,9015)
      
      ! Special handling if averaging requested:
      
      if ( (j >= 24 .and. j <= 26) .or. (j >= 28 .and. j <= 30)) then
         if (j >= 24 .and. j <= 26) then
           jdx = j-23
         else
           jdx = j-24
         end if
         
         ! If ITIMAV=0, do not turn on averaging, even if requesting (e.g.
         ! if this is a MIN run).
         
         if (itimav == 0) then
            write(iout,9013)
            goto 11
         end if
         if (nave(jdx) == 0) nave(jdx) = max(iinc,1)
         !             IF (NEXACT(JDX).EQ.0) NEXACT(JDX) = MAX(IMULT,1)
         if (ipower(jdx) == 0) ipower(jdx) = nint(value2)
         imult = 0
         if (value1 < zerner) value1 = 1.0d+7
         value2 = value1
      end if
      
      if (value1 == ZERO ) value1 = zerner
      if (value2 == ZERO ) value2 = zerner
      if (value1 < ZERO ) value1 = ZERO
      if (value2 < ZERO ) value2 = ZERO
      
      if (istep1 <= istep2 .or. istep2 <= 0) then
         iwtstp(1,i) = istep1
         iwtstp(2,i) = istep2
      else
         iwtstp(1,i) = istep2
         iwtstp(2,i) = istep1
      end if
      iwtstp(3,i) = max(1,iinc)
      
      ! Store special information if TYPE = SHORT
      
      if (type(1:5) == 'SHORT') then
         slope = value1
         rint = value2
         iwtstp(1,i) = istep1
         iwtstp(2,i) = istep2
         if (iinc <= 0) iwtstp(3,i) = 0
         ishrtb = i
         
      else if (istep2 <= 0 .or. istep1 == istep2) then
         slope = ZERO
         rint = value1
         
         ! If IMULT > 0, then use a constant multipier every IINC steps, instead
         ! of a linear interpolation. Set IWTSTP(3,I) = -IINC, to signal this type
         ! of change, store the multiplicative factor in WTNRG(1,I), and store
         ! the initial value in WTNRG(2,I).
         
      else if (imult > 0) then
         if (abs(value1) < zerner) then
            write(iout,9011)
            call mexit(iout, 1)
         else if (value2/value1 < ZERO ) then
            write(iout,9012)
            call mexit(iout, 1)
         end if
         ntu = (iwtstp(2,i)-iwtstp(1,i))/iwtstp(3,i)
         rntu = max(1,ntu)
         slope = (value2/value1)**(ONE/rntu)
         rint = value1
         iwtstp(3,i) = -iwtstp(3,i)
      else
         step1 = istep1
         step2 = istep2
         slope = (value2 - value1) / (step2 - step1)
         rint = value2 - slope*step2
      end if  ! (type(1:5) == 'SHORT')
      wtnrg(1,i) = slope
      wtnrg(2,i) = rint
   end do
   
   ! Now read the redirection cards, if any. The first non-blank line encountered
   ! is assumed to be the start of the redirection section. If it is not a known
   ! option, assume no redirection cards where given.
   
   20 read(in,9004,end=25) aline
   
   ! skip blank lines
   
   do j = 1,len(aline)
      if (aline(j:j) /= ' ') goto 29
   end do
   goto 20
   
   ! Now see if this line starts with a valid redirection name...
   
   29 continue
   do m=1,9
      do n=10,1,-1
         if (redirc(m)(n:n) /= ' ') exit
      end do
      call chklin(aline,redirc(m),n,imat,istrt,ilengt,iout)
      if (imat /= 0) then
         redir(m) = aline(istrt:istrt+ilengt-1)
         iredir(m) = ilengt
         ihits = ihits + 1
         if (ihits == 1) write(iout,9069)
         write(iout,9070) redirc(m),redir(m)(1:ilengt)
         goto 20
      end if
   end do
   if (ihits == 0) write(iout,9071)
   
   backspace(in)
   
   ! Now read the restraints.
   
   25 continue
   !     ------ open scratch files for shifts and NOESY intensities,
   !            if needed
   
   if (iredir(4) /= 0) open (81,status='SCRATCH')
   if (iredir(5) /= 0) open (82,status='SCRATCH')
   if (iredir(7) /= 0) open (53,status='SCRATCH')
   if (iredir(8) /= 0) open (57,status='SCRATCH')
   if (iredir(9) /= 0) open (58,status='SCRATCH')
   
   ! If restraint input has been redirected, open the appropriate file
   
   if (iredir(3) /= 0) then
      call amopen(iscop1,redir(3)(1:iredir(3)),'O','F','R')
      iin = iscop1
      write(iout,9031) redir(3)(1:iredir(3))
      
      !  --- read and echo comment cards from this file:
      !  --- Mengjuei Hsieh: this will not print all the the comments
      !                      lines, because it starts parsing once
      !                      the first char is # and never goes back.
      
      write(6,'(a)') 'Here are comments from the DISANG input file:'
      42 read(iin,'(a)',end=43) aline
      if (aline(1:1) == '#') then
         write(6,'(a)') aline
         goto 42
      end if
      backspace (iin)
      write(6,'(a)')
      43 continue
   end if
   
   ! If echoing of restraints read has been redirected, open the appropriate file:
   
   iuse = 0
   if (iredir(1) /= 0) then
      if (redir(1)(1:iredir(1)) == 'POUT') then
         iuse = iout
      else
         call amopen(iscop2,redir(1)(1:iredir(1)),owrite,'F','W')
         iuse = iscop2
         write(iuse,9032)
         write(iout,9033) redir(1)(1:iredir(1))
      end if
   end if
   
   ! Set the 0 defaults for many variables which can be set by the
   ! "rst" namelist read. The variables which appear before the "do i=1,999999"
   ! loop take the value of the last &rst namelist that set them. The remainder
   ! revert to the defaults (0) in the "do i=1,999999" loop:
   
   ir6 = 0
   ialtd = 0
   ninc = 0
   imult = 0
   
   r1 = ZERO
   r2 = ZERO
   r3 = ZERO
   r4 = ZERO
   rk2 = ZERO
   rk3 = ZERO
   r1a = ZERO
   r2a = ZERO
   r3a = ZERO
   r4a = ZERO
   rk2a = ZERO
   rk3a = ZERO
 

! comienzo del loop lectura de los restraints indice i

   do i=1,999999
      do ii=1,8
         iat(ii) = 0
         atnam(ii) = '    '
      end do
      do ii = 1,4
        rstwt(ii) = 0.0
      end do
      restraint = ' '
      do ii = 1,3
         rjcoef(ii) = ZERO
      end do
      do ii=1,maxigr
         igr1(ii) = 0
         igr2(ii) = 0
         igr3(ii) = 0
         igr4(ii) = 0
         igr5(ii) = 0
         igr6(ii) = 0
         igr7(ii) = 0
         igr8(ii) = 0
         grnam1(ii) = '    '
         grnam2(ii) = '    '
         grnam3(ii) = '    '
         grnam4(ii) = '    '
         grnam5(ii) = '    '
         grnam6(ii) = '    '
         grnam7(ii) = '    '
         grnam8(ii) = '    '
      end do
      do ii=1,3
         tgtvec(3*(i-1)+ii) = ZERO
      end do
      iconstr = 0
      ifconstr(i) =0
      ixpk = 0
      nxpk = 0
      nstep1 = 0
      nstep2 = 0
      irstyp = 0
      ifvari = 0
      iresid = 0
      ifntyp = 0
      r0 = ZERO
      r0a = ZERO
      k0 = ZERO
      k0a = ZERO
      
      ! Look for "rst" namelist. If not found, assume we are done reading them
      
      if( iredir(3) /= 0 ) then
         call nmlsrc('rst',iin,ifind)
         if (ifind == 0) goto 227
         
         read (iin,nml=rst,end=227)
      end if
      
      do jjj=1,maxigr
        igrarr(jjj,1) = igr1(jjj)
        igrarr(jjj,2) = igr2(jjj)
        igrarr(jjj,3) = igr3(jjj)
        igrarr(jjj,4) = igr4(jjj)
        igrarr(jjj,5) = igr5(jjj)
        igrarr(jjj,6) = igr6(jjj)
        igrarr(jjj,7) = igr7(jjj)
        igrarr(jjj,8) = igr8(jjj)
        grnamarr(jjj,1) = grnam1(jjj)
        grnamarr(jjj,2) = grnam2(jjj)
        grnamarr(jjj,3) = grnam3(jjj)
        grnamarr(jjj,4) = grnam4(jjj)
        grnamarr(jjj,5) = grnam5(jjj)
        grnamarr(jjj,6) = grnam6(jjj)
        grnamarr(jjj,7) = grnam7(jjj)
        grnamarr(jjj,8) = grnam8(jjj)
      end do
      
      ! Added by Matthew Seetin 9/2007 to enable natural language input and planar
      ! restraints.
      
      ! If iat is specified in the old style, don't parse the restraint variable.

      227 if (restraint /= ' ' .and. iat1 ==0) then
            call parsrest(restraint,iat,igrarr,rstwt,name,ipres,nres,iout)
      else if (restraint /= ' ' .and. iat1 /=0) then
        write(iout,9079) ! Print a warning message if restraint and &
                         ! iat1 are both specified.  Only iat() will be used.
      end if
      
      !If IAT1=0, assume the last restraint has been read.
      
      if (iat1 == 0) then
         nmrnum = i-1
         if (nmrnum <= 0) write(iout,9025)
         if (nmrnum > 0) write(iout,9026) nmrnum
         write(iout,9027)
         goto 100
      end if
      
#ifdef MPI /* SOFT CORE */
      if (ifsc /= 0) then
         oneweight = 1.0d0 / ( 1.0d0 - clambda )
         write (6,'(a,f8.3,a,f8.3)') 'Scaling up weight of ',rk2,' by ',oneweight
         rk2 = rk2 * oneweight
         rk3 = rk3 * oneweight
         write (6,'(a,f8.3)') 'New weight ',rk2
      end if
#endif
      
      ! Figure out the type of this restraint
      
      resttype(i) = 0
      
      if (iat3 == 0) then        ! Distance restraint
        resttype(i) = 1  
        if (abs(rstwt(1)) > rstwttol) then
          rk2 = rk2 * rstwt(1)
          rk3 = rk3 * rstwt(1)
        end if
      else if (iat4 == 0) then   ! Angle restraint
        resttype(i) = 2
      else if (iat5 == 0) then
        if (abs(rstwt(1)) < rstwttol .and. abs(rstwt(2)) < rstwttol) then
          resttype(i) = 3         ! Torsion or j-coupling restraint
        else if (abs(rstwt(1)) < rstwttol) then  
        
        ! A generalized distance with only one non-zero weight, reduces to 
        ! a normal distance restraint
        
          iat1 = iat3
          iat2 = iat4
          iat3 = 0
          iat4 = 0
          rk2 = rk2 * rstwt(2)
          rk3 = rk3 * rstwt(2)
          rk2a = rk2a * rstwt(2)
          rk3a = rk3a * rstwt(2)
          resttype(i) = 1
        else if (abs(rstwt(2)) < rstwttol) then
          iat3 = 0
          iat4 = 0
          rk2 = rk2 * rstwt(1)
          rk3 = rk3 * rstwt(1)
          rk2a = rk2a * rstwt(1)
          rk3a = rk3a * rstwt(1)
          resttype(i) = 1
        else
          resttype(i) = 6
        end if
        rstwt(3) = 0.0
        rstwt(4) = 0.0
      else if (iat6 == 0) then
        resttype(i) = 4
      else if (iat7 == 0) then
        if (abs(rstwt(1)) > rstwttol .and. abs(rstwt(2)) > rstwttol .and. abs(rstwt(3)) > rstwttol) then
          resttype(i) = 6
        else if (abs(rstwt(1)) < rstwttol .and. abs(rstwt(2)) < rstwttol .and. abs(rstwt(3)) < rstwttol) then
          write(iout,'(a)') 'Error: At least one weight must be non-zero in a generalized restraint.'
          call mexit(iout, 1)
        else                  
        ! What I'm doing here is compacting the generalized restraint down so it only uses the parts that correspond to
        ! non-zero weights
          ii = 1
          do 
            if (ii > 3) then
              resttype(i) = 6
              exit
            end if
            if (abs(rstwt(ii)) < rstwttol) then
              if (ii < 3) then
                do jjj=2*ii-1,4
                  iat(jjj) = iat(jjj+2)
                end do
                if (ii == 2) then
                  rstwt(2) = rstwt(3)
                else
                  rstwt(1) = rstwt(2)
                  rstwt(2) = rstwt(3)
                end if
                rstwt(3) = 0.0
                if (iat5 == 0) then
                  iat3 = 0
                  iat4 = 0
                  rk2 = rk2 * rstwt(1)
                  rk3 = rk3 * rstwt(1)
                  rk2a = rk2a * rstwt(1)
                  rk3a = rk3a * rstwt(1)
                  resttype(i) = 1
                  exit
                else
                  iat5 = 0
                  iat6 = 0
                end if
              else
                iat5 = 0
                iat6 = 0
                resttype(i) = 6
                exit
              end if
            else
              ii = ii + 1
            end if
          end do
        end if
        rstwt(4) = 0.0
      else if (iat8 /= 0) then
        if (abs(rstwt(1)) < rstwttol .and. abs(rstwt(2)) < rstwttol .and. &
            abs(rstwt(3)) < rstwttol .and. abs(rstwt(4)) < rstwttol) then
          resttype(i) = 5
        else if (abs(rstwt(1)) > rstwttol .and. abs(rstwt(2)) > rstwttol .and. & 
                 abs(rstwt(3)) > rstwttol .and. abs(rstwt(4)) > rstwttol) then
          resttype(i) = 6
        else
        ! What I'm doing here is compacting the generalized restraint down so it only uses the parts that correspond to
        ! non-zero weights
          ii = 1
          do 
            if (ii > 4) then
              resttype(i) = 6
              exit
            end if
            if (abs(rstwt(ii)) < rstwttol) then
              if (ii < 4) then
                do jjj=2*ii-1,6
                  iat(jjj) = iat(jjj+2)
                end do
                do jjj = 1,3
                  rstwt(ii) = rstwt(ii+1)
                end do
                rstwt(4) = 0.0
                if (iat5 == 0) then
                  iat3 = 0
                  iat4 = 0
                  rk2 = rk2 * rstwt(1)
                  rk3 = rk3 * rstwt(1)
                  rk2a = rk2a * rstwt(1)
                  rk3a = rk3a * rstwt(1)
                  resttype(i) = 1
                  exit
                else if (iat7 == 0) then
                  iat5 = 0
                  iat6 = 0
                else
                  iat7 = 0
                  iat8 = 0
                end if
              else
                iat7 = 0
                iat8 = 0
                resttype(i) = 6
                exit
              end if
            else
              ii = ii + 1
            end if
          end do
        end if
      else
        write(iout,9075)
        call mexit(iout,1)
      end if
      
      !  Restraint input will now support more easily inputting a simple parabolic well
      !  This input will use r0, k0, r0a, and k0a.  A changing spring constant will not be 
      !  used if jar = 1.  The following code maps the above variables to r1-r4 and rk2-rk3.
      !  If the user specifies non-zero spring constants in the traditional style AND in this
      !  new style, the new-style input will overwrite the old.
      
      if ( (abs(k0) > ZERO) .or. ( (ifvari == 1) .and. (abs(k0a) > ZERO) ) ) then
        r1 = ZERO
        if(resttype(i) == 2 .or. resttype(i) == 4 .or. resttype(i) == 5) then
           r4 = 180.0 ! The acos function only returns values between 0 and 180
        else if ((resttype(i) == 3) .and. (abs(rjcoef(1)) < zerner) .and. (abs(rjcoef(2)) < zerner)) then
          r1 = r0 - 180.0
          r4 = r0 + 180.0
        else
          r4 = r0 + 500.0 
        end if
        r2 = r0
        r3 = r0
        rk2 = k0
        rk3 = k0
        if(ifvari == 1) then
          r1a = ZERO
          r2a = r0a
          r3a = r0a
          if(resttype(i) == 2 .or. resttype(i) == 4 .or. resttype(i) == 5) then
            r4a = 180.0 ! The acos function only returns values between 0 and 180
          else if ((resttype(i) == 3) .and. (abs(rjcoef(1)) < zerner) .and. (abs(rjcoef(2)) < zerner)) then
            r1a = r0a - 180.0
            r4a = r0a + 180.0
          else
            r4a = r0a + 500.0 
          end if
          if(jar == 1 .and. abs(k0a) > zerner) then
            write(iout,'(a)') ' Warning: Varying force constant not recommended when jar=1. Ignoring k0a.'
            write(iout,'(a)') '          User may override this by defining the energy well with r1-r4, ' 
            write(iout,'(a)') '          r1a-r4a, rk2, rk3, rk2a and rk3a only.'
          else
            rk2a = k0a
            rk3a = k0a
          end if
        end if
      end if


      ! Assignment of some variables if jar run (roit. 02/27/05)
      ! We are using the restraint potential as defined for an umbrella sampling
      ! calculation, where nstep2=nstlim, and the minimum of the well is
      ! linearly changed from r2=r3 to r2a=r3a. Cosntant is rk2=rk3 and r1 and r4
      ! are taken large enough [r1(r4)=r2(r3)-(+)100]

      if(i.eq.1.and.jar.eq.1) then
      write(6,*) 'jar option running '
!        if(nmrnum /= 1) then
!        write(*,*) 'Error asigando por roitberg si nmrnum > 1'
!        goto 1004
!        endif
        nstep2=nstlim
        ifvari=1

        if(r2.ne.r3) then
           if(r2.eq.0.) then
              r2=r3
           elseif(r3.eq.0.) then
              r3=r2
           else
              write(6,'(/2x,a)') 'r2 and r3 must be equal, or only one of them should be defined'
           endif
        endif

        if(r2a.ne.r3a) then
           if(r2a.eq.0.) then
              r2a=r3a
           elseif(r3a.eq.0.) then
              r3a=r2a
           else
              write(6,'(/2x,a)') 'r2a and r3a must be equal, or only one of them should be defined'
           endif
        endif

        if(rk2.ne.rk3) then
           if(rk2.eq.0.) then
              rk2=rk3
           elseif(rk3.eq.0.) then
              rk3=rk2
           else
              write(6,'(/2x,a)') 'rk2 and rk3 must be equal, or only one of them should be defined'
           endif
        endif
       
        if(rk2a.ne.rk2) then
           if(rk2a.eq.0.) then
              rk2a=rk2
           else if(rk2.eq.0.) then
              rk2=rk2a
           else
              write(6,'(/2x,a)') 'rk2 and rk2a must be equal, or only one of them should be defined'
           end if
        end if
        
        if(rk2a.ne.rk3a) then
           if(rk2a.eq.0.) then
              rk2a=rk3a
           elseif(rk3a.eq.0.) then
              rk3a=rk2a
           else
              write(6,'(/2x,a)') 'rk2a and rk3a must be equal, or only one of them should be defined'
           endif
        endif

        if(rk2.eq.0.) then
            write(6,'(/2x,a)') 'rk2 should not be 0'
        end if
        
        r1  = r2  - 100.0
        r1a = r2a - 100.0
        r4  = r3  + 100.0
        r4a = r3a + 100.0

        ! Modified by Matthew Seetin 11/2007 to be consistent with r0 and k0 input
        drjar=(r2a-r2)/(nstep2-nstep1)
        
      elseif (i.eq.1.and.jar.eq.0) then
      
! Asignacion de variables Morse 
! si jar = 0, Morse  es el primer restraint
! r2=r3=req  r1 y r4 suficientemente gdes
! rk2=rk3 its the binding energy
! iat(3) and iat(4) must be 0 (only distnace allowed)
        if (morse == 1 ) then
        write(6,*) '**********************************'
        write(6,*) 'Morse option running'
        if(iat(3).ne.0.or.iat(4).ne.0) & 
        write(6,*) 'Morse only allowed for distance'
        if(r2.eq.0.and.r3.eq.0) &
        write(6,*) 'equilibrium distance (r2 oR r3) cannot be 0' 
        if(rk2.eq.0.and.rk3.eq.0) &
        write(6,*) 'binding energy (rk2 oR rk3) cannot be 0'  
        r3=r2
        rk3=rk2
        r1=r2-100.0
        r4=r2+100.0
        r1a=r1
        r2a=r2
        r3a=r3
        r4a=r4
        write(6,*) 'Morse potential r0: ',r2,'De. ',rk2
        write(6,*) '**********************************'
        
! fin del loop asignacion Morse
        endif
      
      endif       

! Aignacion Morse si Jar = 1 y Morse =1 en el restr 2
        if(i.eq.2.and.jar.eq.1.and.Morse.eq.1) then
        write(6,*) '**********************************'
        write(6,*) 'Morse and Jar option running'
        if(iat(3).ne.0.or.iat(4).ne.0) &
        write(6,*) 'Morse only allowed for distance'
        if(r2.eq.0.and.r3.eq.0) &
        write(6,*) 'equilibrium distance (r2 oR r3) cannot be 0'
        if(rk2.eq.0.and.rk3.eq.0) &
        write(6,*) 'binding energy (rk2 oR rk3) cannot be 0'
        r3=r2
        rk3=rk2
        r1=r2-100.0
        r4=r2+100.0
        r1a=r1
        r2a=r2
        r3a=r3
        r4a=r4
        write(6,*) 'Morse potential r0: ',r2,'De. ',rk2                        
        write(6,*) '**********************************'
! fin del loop asignacion Morse (jar =1)
        endif


  ! Check to make sure maximum allowable number of weight cards not exceeded
      
      if (i > maxrst) then
         write(iout,9028) maxrst
         call mexit(iout, 1)
      end if
      

      
      ! If IRESID = 1, read another line and determine the absolute atom numbers
      ! of the atoms making up the restraint. The values in IAT1->IAT8 will
      ! be replaced by these absolute atom numbers. If the natural language input
      ! has been used, then this step has already been taken care of.  However, 
      ! iresid = 1 may still be necessary if the user has specified atom names in 
      ! atom groups.
      ! Updated 9/2007 by Matthew Seetin
      
      if (iresid == 1 .and. restraint == ' ') then
         ierr = 0
         do jjj=1,8
           if (iat(jjj) == 0) then
             exit
           else
             call getnat(iat(jjj),atnam(jjj),name,ipres,nres,iout,ierr)
           end if
         end do
         if (ierr == 1) then
            write(iout,9073)
            call mexit(iout, 1)
         end if
      end if
      
      ! Check and see if atoms in a plane restraint are the same.  If atom 1 and
      ! atom 3 are the same, that's allowed, but if atom 1 and atom 2 are the
      ! same, it will result in a divide by zero.

      if (iat5 /= 0) then
        do jjj=1,7,2
          if (iat(jjj) == iat(jjj+1) .and. iat(jjj) > 0 ) then
            write(iout,9078) jjj,jjj+1
            call mexit(iout,1)
          end if
        end do
      end if

      ! Set flag IGRAVT (GRoup AVerage Type) to IR6
      
      igravt(i) = ir6
      
      ! Set flag IALTDIS (Use alternative restraint function) to IALTD
      
      ialtdis(i) = ialtd
      
      ! Set the function-type flag to IFNTYP
      
      nmrfty(i) = ifntyp
      
      ! Store the atom pointers (3*(IAT-1)) and effective ranges. The atom
      ! pointers will be changed with the group info., if IAT1 or IAT2 < 0
      ! (see below).
      
      do jjj=1,8
        nmrat(jjj,i) = 3*(iat(jjj)-1)
      end do
      
      ! Initialize nmrat(9,I) - nmrat(16,i) to 0
      ! These will later be set to a negative number if a group is defined.
      
      do jjj=9,16
        nmrat(jjj,i) = 0
      end do
      
      do jjj = 1,4
        rstwtarr(jjj,i) = rstwt(jjj)
      end do
      
      if (nstep2 > nstep1 .or. nstep2 <= 0) then
         nmrst(1,i) = nstep1
         nmrst(2,i) = nstep2
      else
         nmrst(1,i) = nstep2
         nmrst(2,i) = nstep1
      end if
      nmrst(3,i) = max(1,ninc)
      kxpk(1,i) = ixpk
      kxpk(2,i) = nxpk
      
      ! Store the Karplus equation constants:
      
      ajcoef(1,i) = rjcoef(1)
      ajcoef(2,i) = rjcoef(2)
      ajcoef(3,i) = rjcoef(3)
      
      ! Added 9/2007 by Matthew Seetin to enable group restraints on angles, 
      ! torsions, and my new planar restraints
      ! If IAT is < 0, for any of the atoms specified in the restraint, then read
      ! and sort the group definition information for atom jjj of this restraint
      
      usecom = .false.
      
      do jjj=1,8
        if(iat(jjj) < 0) then
          do iii=1,maxigr
            igrpass(iii) = igrarr(iii,jjj)
            grnampass(iii) = grnamarr(iii,jjj)
          end do
          call nmrgrp(i,jjj-1,lastpk,iscrth,natom,nmrat,nmrcom,iin, &
            iout,maxgrp,igrpass,grnampass,iresid,name,ipres, &
            nres,maxigr)
        end if
        usecom = usecom .or. nmrat(jjj+8,i) < 0
        ! usecom tracks whether a group is involved in the restraint
      end do
      
      ! Calculate the current value of the internal. If IRSTYP=1, this
      ! will be used to calculate the target distance. If IRSTYP=0, it
      ! will be used simply to report the deviation between the current & target
      ! values.
      ! Modified by Matthew Seetin 9/2007.  This now includes atom group support
      ! for angles and torsions, and there are plane-plane and plane-point angle 
      ! restraints.
      
      jcoupl = .false.
      resttypetmp = resttype(i) 
      select case (resttypetmp)
        case (1)  ! distance
          if (.not. usecom) then
            call disnrg(x,f,dfr,rcurr,nmrat(1,i),nmrat(2,i),e, &
                  r1,r2,r3,r4,rk2,rk3,0,dumma,dumma,nave(1), &
                  nexact(1),ipower(1),dumm,ravein(1),dumm, &
                  0,0,0,0,1,ialtdis(i))
            ! calculate and normalize constraining vectors
            if(iconstr>0) then
               tgtvec(3*(i-1)+1) = x(nmrat(2,i)+1)-x(nmrat(1,i)+1)
               tgtvec(3*(i-1)+2) = x(nmrat(2,i)+2)-x(nmrat(1,i)+2)
               tgtvec(3*(i-1)+3) = x(nmrat(2,i)+3)-x(nmrat(1,i)+3)
               tgtvec(3*(i-1)+1) = tgtvec(3*(i-1)+1) / rcurr
               tgtvec(3*(i-1)+2) = tgtvec(3*(i-1)+2) / rcurr
               tgtvec(3*(i-1)+3) = tgtvec(3*(i-1)+3) / rcurr
            end if
            ifconstr(i) = iconstr
          else if (igravt(i) == 1) then
            call r6ave(x,nmrat,nmrcom,rcurr,i)
            call disnrg(x,f,dfr,rcurr,nmrat(1,i),nmrat(2,i),e, &
                  r1,r2,r3,r4,rk2,rk3,0,dumma,dumma,nave(1), &
                  nexact(1),ipower(1),dumm,ravein(1),dumm, &
                  0,0,0,0,3,ialtdis(i))
          else
            call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
            call disnrg(xcom,f,dfr,rcurr,0,3,e,r1,r2,r3,r4,rk2, &
                  rk3,0,dumma,dumma,nave(1),nexact(1), &
                  ipower(1),dumm,ravein(1),dumm,0,0,0,0,1, &
                  ialtdis(i))
            ! calculate and normalize constraining vectors
            if(iconstr>0) then
               tgtvec(3*(i-1)+1) = xcom(4) - xcom(1)
               tgtvec(3*(i-1)+2) = xcom(5) - xcom(2)
               tgtvec(3*(i-1)+3) = xcom(6) - xcom(3)
               tgtvec(3*(i-1)+1) = tgtvec(3*(i-1)+1) / rcurr
               tgtvec(3*(i-1)+2) = tgtvec(3*(i-1)+2) / rcurr
               tgtvec(3*(i-1)+3) = tgtvec(3*(i-1)+3) / rcurr
            end if
            ifconstr(i) = iconstr
         end if
      case (2) ! angle
        if (usecom) then
          call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
          call angnrg(xcom,f,dfr,rcurr,0,3,6, &
               e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma,nave(2), &
               nexact(2),ipower(2),dumm,ravein(2),dumm, &
               0,0,0,0,1)
        else
          call angnrg(x,f,dfr,rcurr,nmrat(1,i),nmrat(2,i),nmrat(3,i), &
               e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma,nave(2), &
               nexact(2),ipower(2),dumm,ravein(2),dumm, &
               0,0,0,0,1)
        end if
      case (3) ! Torsion or j-coupling
         
         ! If J coefficients are 0 for this torsion, it is a normal torsional
         ! restraint. If they are .NE.0, then it is a J-coupling factor restraint.
         
         jcoupl = abs(ajcoef(1,i)) > zerner .or. &
               abs(ajcoef(2,i)) > zerner
         
         if (.not.jcoupl) then
           if (usecom) then
             call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
             call tornrg(xcom,f,dfr,rcurr,0,3,6,9, &
                  e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma, &
                  dumma,dumma,nave(1),nexact(1),ipower(1),dumm, &
                  ravein(1),dumm,0,0,0,0,1)
           else
             call tornrg(x,f,dfr,rcurr,nmrat(1,i),nmrat(2,i),nmrat(3,i), &
                  nmrat(4,i),e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma, &
                  dumma,dumma,nave(1),nexact(1),ipower(1),dumm, &
                  ravein(1),dumm,0,0,0,0,1)
           end if
         else
           if (usecom) then
             ! Since atom groups are not supported for j-coupling restraints, 
             ! print out an error message and terminate execution if the user 
             ! specifies atom groups AND j-coupling
             
             write(iout,9074)
             call mexit(iout, 1)
           else
             call jnrg(x,f,rcurr,nmrat(1,i),nmrat(2,i),nmrat(3,i), &
                  nmrat(4,i),e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma, &
                  dumma,dumma,nave(1),nexact(1),ipower(1),dumm, &
                  ravein(1),dumm,0,0,0,0,ajcoef(1,i),1)
           end if
         end if
      case (4)  ! This is a plane-point angle restraint
        if (usecom) then
          call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
          call commass(comrimass,nmrat,nmrcom,rimass,i)
          call plptnrg(xcom,f,dfr,comrimass,rcurr,0,3,6,9,12, &
               e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma,nave(4), &
               nexact(4),ipower(4),dumm,ravein(4),dumm, &
               0,0,0,0,1)
        else
          call plptnrg(x,f,dfr,rimass,rcurr,nmrat(1,i),nmrat(2,i), &
               nmrat(3,i),nmrat(4,i),nmrat(5,i), &
               e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma,nave(4), &
               nexact(4),ipower(4),dumm,ravein(4),dumm, &
               0,0,0,0,1)
        end if
      case (5)  ! This is a plane-plane angle restraint
        if (usecom) then
          call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
          call plnnrg(xcom,f,dfr,rcurr,0,3,6,9,12,15,18,21, &
               e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma,nave(5), &
               nexact(5),ipower(5),dumm,ravein(5),dumm, &
               0,0,0,0,1)
        else
          call plnnrg(x,f,dfr,rcurr,nmrat(1,i),nmrat(2,i), &
               nmrat(3,i),nmrat(4,i),nmrat(5,i),nmrat(6,i),nmrat(7,i),nmrat(8,i),&
               e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma,nave(5), &
               nexact(5),ipower(5),dumm,ravein(5),dumm, &
               0,0,0,0,1)
        end if
      case (6) ! Generalized distance coordinate
        do jjj=1,4
          rstwtpass(jjj) = rstwtarr(jjj,i)
        end do
        if (usecom) then
          call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
          call gendisnrg(xcom,f,dfr,rcurr,0,3,6,9,12,15,18,21, &
               rstwtpass,e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma,nave(6), &
               nexact(6),ipower(6),dumm,ravein(6),dumm, &
               0,0,0,0,1)
        else
          call gendisnrg(x,f,dfr,rcurr,nmrat(1,i),nmrat(2,i), &
               nmrat(3,i),nmrat(4,i),nmrat(5,i),nmrat(6,i),nmrat(7,i),nmrat(8,i),&
               rstwtpass,e,r1,r2,r3,r4,rk2,rk3,0,dumma,dumma,nave(6), &
               nexact(6),ipower(6),dumm,ravein(6),dumm, &
               0,0,0,0,1)
        end if
      case default
        ! Improper restraint specified
        write(iout,9075)
        call mexit(iout, 1)
      end select
      
      rint = ZERO
      if (irstyp == 1) rint = rcurr
      
      ! Echo the restraint to the user, if the user has requested this
      
      do jjj=1,8
        if (iat(jjj) < 0) then
          nametmp(jjj) = 'COM '
        else if (iat(jjj) > 0) then
          nametmp(jjj) = name(iat(jjj))
        else
          nametmp(jjj) = '????'
        end if
      end do
      
      if (iuse /= 0) then
        resttypetmp = resttype(i)
        select case (resttypetmp)
          case (1)
            write(iuse,9040) nametmp(1),iat1,nametmp(2),iat2, &
                     nstep1,nstep2
          case (2)
            write(iuse,9041) nametmp(1),iat1,nametmp(2),iat2, &
                  nametmp(3),iat3,nstep1,nstep2
          case (3)
            write(iuse,9042) nametmp(1),iat1,nametmp(2),iat2, &
                  nametmp(3),iat3,nametmp(4),iat4, &
                  nstep1,nstep2
          case (4)
            write(iuse,9048) nametmp(1),iat1,nametmp(2),iat2, &
                  nametmp(3),iat3,nametmp(4),iat4,nametmp(5),iat5, &
                  nstep1,nstep2
          case (5)
            write(iuse,9049) nametmp(1),iat1,nametmp(2),iat2, &
                  nametmp(3),iat3,nametmp(4),iat4,nametmp(5),iat5, &
                  nametmp(6),iat6,nametmp(7),iat7,nametmp(8),iat8, &
                  nstep1,nstep2
          case (6)
            if (iat5 == 0) then
              write(iuse,9042) nametmp(1),iat1,nametmp(2),iat2, &
                    nametmp(3),iat3,nametmp(4),iat4, &
                    nstep1,nstep2
            else if (iat7 == 0) then
              write(iuse,9049) nametmp(1),iat1,nametmp(2),iat2, &
                    nametmp(3),iat3,nametmp(4),iat4,nametmp(5),iat5, &
                    nametmp(6),iat6
            else 
              write(iuse,9049) nametmp(1),iat1,nametmp(2),iat2, &
                    nametmp(3),iat3,nametmp(4),iat4,nametmp(5),iat5, &
                    nametmp(6),iat6,nametmp(7),iat7,nametmp(8),iat8, &
                    nstep1,nstep2
            end if
         end select
         if (ninc > 1) write(iuse,9043) ninc
         
         ! If a group-type distance restraint was defined, echo the sorted group:
         do jjj=1,8
           if (nmrat(jjj,i) < 0 .and. nmrat(jjj+8,i) < 0 .and. jjj == 1) then
             write(iuse,9065)
             write(iuse,9066) (nmrcom(1,kk),nmrcom(2,kk), &
                  kk=-nmrat(jjj,i),-nmrat(jjj+8,i))
           else if (nmrat(jjj,i) < 0 .and. nmrat(jjj+8,i) < 0 .and. jjj == 2) then
             write(iuse,9067)
             write(iuse,9066) (nmrcom(1,kk),nmrcom(2,kk), &
                  kk=-nmrat(jjj,i),-nmrat(jjj+8,i))
           else if (nmrat(jjj,i) < 0 .and. nmrat(jjj+8,i) < 0 .and. jjj == 3) then
             write(iuse,9076)
             write(iuse,9066) (nmrcom(1,kk),nmrcom(2,kk), &
                  kk=-nmrat(jjj,i),-nmrat(jjj+8,i))
           else if (nmrat(jjj,i) < 0 .and. nmrat(jjj+8,i) < 0) then
             write(iuse,9077) jjj
             write(iuse,9066) (nmrcom(1,kk),nmrcom(2,kk), &
                  kk=-nmrat(jjj,i),-nmrat(jjj+8,i))
           end if
         end do
      end if  ! (iuse /= 0)
      
      convrt = ONE
      if ((iat3 /= 0 .or. iat4 /= 0) .and. .not.jcoupl .and. resttype(i) /= 6) &
            convrt = RAD_TO_DEG
      rprt = rint*convrt
      rcurr = rcurr*convrt
      
      ! If this is a torsional restraint, translate the value of the torsion
      ! (by +- N*360) to bring it as close as possible to one of the two central
      ! "cutoff" points (r2,r3). Use this as the value of the torsion in the report.
      
      if (iat4 /= 0 .and. .not.jcoupl) then
         apmean = (r2+rprt+r3+rprt)*0.5
         18 if (rcurr-apmean > 180.0d0) then
            rcurr = rcurr - 360.0d0
            goto 18
         else if (apmean-rcurr > 180.0d0) then
            rcurr = rcurr + 360.0d0
            goto 18
         end if
      end if
      
      if (iuse /= 0) then
         write(iuse,9045) r1+rprt,r2+rprt,r3+rprt,r4+rprt,rk2,rk3
         if (ifvari > 0) &
               write(iuse,9046) r1a+rprt,r2a+rprt,r3a+rprt,r4a+rprt,rk2a, &
               rk3a
         
         dev2 = min(abs(rcurr-r2-rprt),abs(rcurr-r3-rprt))
         if (rcurr >= r2+rprt .and. rcurr <= r3+rprt) dev2 = 0
         write(iuse,9047) rcurr,abs(rcurr-(r2+rprt+r3+rprt)/2.0d0), &
               dev2
      end if
      
      ! The energy routines only work correctly if r1->r4 increase monotonically.
      ! If not, flag an error -- dap (8/91)
      
      if (r1 > r2 .or. r2 > r3 .or. r3 > r4) then
         goto 1003
      end if
      if (nstep2 /= 0 .and. ifvari > 0 .and. nstep1 /= nstep2) then
         if (r1a > r2a .or. r2a > r3a .or. r3a > r4a) then
            goto 1003
         end if
      end if
      
      
      ! Convert R1->R4 to radians, if it is an angle
      
      if ((iat3 /= 0 .or. iat4 /= 0) .and. .not.jcoupl .and. resttype(i) /= 6) then
         r1 = r1/convrt
         r2 = r2/convrt
         r3 = r3/convrt
         r4 = r4/convrt
         r1a = r1a/convrt
         r2a = r2a/convrt
         r3a = r3a/convrt
         r4a = r4a/convrt
      end if
      
      ! Now calculate the slope/intercepts for the r and k values with step
      ! number and store these in the appropriate arrays
      
      if (nstep2 == 0 .or. ifvari <= 0 .or. nstep1 == nstep2) then
         r1nmr(1,i) = ZERO
         r1nmr(2,i) = r1+rint
         r2nmr(1,i) = ZERO
         r2nmr(2,i) = r2+rint
         r3nmr(1,i) = ZERO
         r3nmr(2,i) = r3+rint
         r4nmr(1,i) = ZERO
         r4nmr(2,i) = r4+rint
         rk2nmr(1,i) = ZERO
         rk2nmr(2,i) = rk2
         rk3nmr(1,i) = ZERO
         rk3nmr(2,i) = rk3
      else
         step1 = nstep1
         step2 = nstep2
         r1nmr(1,i) = (r1a - r1) / (step2 - step1)
         r1nmr(2,i) = (r1a + rint) - r1nmr(1,i) * step2
         r2nmr(1,i) = (r2a - r2) / (step2 - step1)
         r2nmr(2,i) = (r2a + rint) - r2nmr(1,i) * step2
         r3nmr(1,i) = (r3a - r3) / (step2 - step1)
         r3nmr(2,i) = (r3a + rint) - r3nmr(1,i) * step2
         r4nmr(1,i) = (r4a - r4) / (step2 - step1)
         r4nmr(2,i) = (r4a + rint) - r4nmr(1,i) * step2
         
         ! If IMULT > 0, then use a constant multipier every IINC steps, instead
         ! of a linear interpolation for the force constants.
         
         if (imult > 0) then
            ntu = (nmrst(2,i)-nmrst(1,i))/nmrst(3,i)
            rntu = max(1,ntu)
            nmrst(3,i) = -nmrst(3,i)
            rk2nmr(1,i) = (rk2a/rk2)**(ONE/rntu)
            rk2nmr(2,i) = rk2
            rk3nmr(1,i) = (rk3a/rk3)**(ONE/rntu)
            rk3nmr(2,i) = rk3
         else
            rk2nmr(1,i) = (rk2a - rk2) / (step2 - step1)
            rk2nmr(2,i) = rk2a - rk2nmr(1,i) * step2
            rk3nmr(1,i) = (rk3a - rk3) / (step2 - step1)
            rk3nmr(2,i) = rk3a - rk3nmr(1,i) * step2
         end if
      end if  ! (nstep2 == 0 .or. ifvari <= 0 .or. nstep1 == nstep2)
      
      
      ! un-Convert R1->R4 to degrees, if it is an angle
      !  (this is done in case the next namelist does not reset these
      !   variables, with the user expecting the same values to be
      !   repeated for several constraints.)
      
      if ((iat3 /= 0 .or. iat4 /= 0) .and. .not.jcoupl) then
         r1 = r1*convrt
         r2 = r2*convrt
         r3 = r3*convrt
         r4 = r4*convrt
         r1a = r1a*convrt
         r2a = r2a*convrt
         r3a = r3a*convrt
         r4a = r4a*convrt
      end if
! fin del loop de lectira de restraints (indice es i)      
   end do
   
   ! Call NMRSHT to set up short-range definitions. Set NAVE(1) (IDUMM(1) in
   ! call list) to zero, so that time-averaged distances are not used on this
   ! (initial) call.
   
   idumm(1) = 0
   idumm(2) = 0
   100 call nmrsht(x,nmrnum,rimass,0,0,nmrat,resttype,nmrcom,nres,ipres,ishrtb, &
         iwtstp,wtnrg,rk2nmr,rk3nmr,nmrst,nmrshb,wtls,nmrfty, &
         igravt,idumm,dumma,dumma,idumm,idumm,idumm,dumma,idumm, &
         dumma,dumm,iout,0)
   
   if (iredir(3) /= 0) close(iin)
   if (iuse /= 0 .and. iuse /= iout) &
      close(iuse) 
   return
   
   ! End-of-file read errors:
   
   1001 write(iout,2000)
   call mexit(iout, 1)
   
   ! No "wt" namelist with type=END read:
   
   1002 write(iout,2002)
   call mexit(iout, 1)
   
   ! r1->r4 were not monotonic error:
   
   1003 write(iout,2003)
   do jjj=1,8
     if (iat(jjj) < 0) then
       nametmp(jjj) = 'COM '
     else
       nametmp(jjj) = name(iat(jjj))
     end if
   end do
   resttypetmp = resttype(i)
   select case (resttypetmp)
    case (1)
      write(iout,9040) nametmp(1),iat1,nametmp(2),iat2, &
               nstep1,nstep2
    case (2)
      write(iout,9041) nametmp(1),iat1,nametmp(2),iat2, &
            nametmp(3),iat3,nstep1,nstep2
    case (3)
      write(iout,9042) nametmp(1),iat1,nametmp(2),iat2, &
            nametmp(3),iat3,nametmp(4),iat4, &
            nstep1,nstep2
    case (4)
      write(iout,9048) nametmp(1),iat1,nametmp(2),iat2, &
            nametmp(3),iat3,nametmp(4),iat4,nametmp(5),iat5, &
            nstep1,nstep2
    case (5)
      write(iout,9049) nametmp(1),iat1,nametmp(2),iat2, &
            nametmp(3),iat3,nametmp(4),iat4,nametmp(5),iat5, &
            nametmp(6),iat6,nametmp(7),iat7,nametmp(8),iat8, &
            nstep1,nstep2
    case (6)
      if (iat5 == 0) then
        write(iout,9042) nametmp(1),iat1,nametmp(2),iat2, &
            nametmp(3),iat3,nametmp(4),iat4, &
            nstep1,nstep2
      else if (iat7 == 0) then
        write(iout,9049) nametmp(1),iat1,nametmp(2),iat2, &
            nametmp(3),iat3,nametmp(4),iat4,nametmp(5),iat5, &
            nametmp(6),iat6
      else 
        write(iout,9049) nametmp(1),iat1,nametmp(2),iat2, &
            nametmp(3),iat3,nametmp(4),iat4,nametmp(5),iat5, &
            nametmp(6),iat6,nametmp(7),iat7,nametmp(8),iat8, &
            nstep1,nstep2
      end if
   end select
   write(iout,9045) r1+rprt,r2+rprt,r3+rprt,r4+rprt,rk2,rk3
   if (ifvari > 0) &
         write(iout,9046) r1a+rprt,r2a+rprt,r3a+rprt,r4a+rprt,rk2a,rk3a
   call mexit(iout,1)
   
   ! Other read errors (probably format mismatch)
   ! (If "#" is in first column, this is a comment line):

   ! 1004 write(iout,9072) !jarzynski error exit (roit. 02/27/05)
   call mexit(iout, 1)

   
   1005 backspace(iin)
   read(iin,9004) aline
   if (aline(1:1) == '#') then
      write(iout,9068) aline
      ibeg = i+1
      goto 9
   else
      write(iout,2001) aline
      call mexit(iout, 1)
   end if
   
   2000 format(' ERROR: End of file read during formatted weight card ', &
         'read.')
   2001 format(' ERROR reading the line:',/,a80,/, &
         ' during formatted weight-card read.')
   2002 format(' ERROR: No "wt" namelist with TYPE=END found')
   2003 format(' ERROR: r1 -> r4 (and r1a -> r4a) must be monotonically ', &
         'increasing; Offending restraint:')
   9000 format(/ /,t12,'Begin reading energy term weight changes/NMR ', &
         'restraints',/,' WEIGHT CHANGES:')
   9001 format(t26,'** No weight changes given **')
   9002 format(/,' RESTRAINTS:')
   9004 format(a80)
   9005 format(a8,2i7,2f12.6,2i7)
   9006 format(1x,a8,2i7,2f12.6,2i7)
   9008 format(' Error: Maximum allowable number of weight change cards', &
         ' (',i5,') exceeded.')
   9010 format(' Error: Invalid TYPE flag in line:')
   9011 format(' Error: Multiplicative-type change (IMULT > 0) ', &
         'invalid when VALUE1 = 0.0.')
   9012 format(' Error: Multiplicative-type change (IMULT > 0) ', &
         'invalid when VALUE1 and',/,t8,'VALUE2 do not have the ', &
         'same sign.')
   9013 format(' Warning: Averaging request ignored.')
   9015 format(' Warning: ISTEP1 >= ISTEP2')
   !9020 FORMAT(11I5)
   9025 format(t27,'** No restraint defined **')
   9026 format(t23,' Number of restraints read = ',i5)
   9027 format(/,t19,'Done reading weight changes/NMR restraints',/ /)
   9028 format(' Error: Maximum allowable number of restraints', &
         ' (',i5,') exceeded.')
   !9029 FORMAT(4(1X,A4))
   !9030 FORMAT(6F12.6)
   9031 format(' Restraints will be read from file: ',a)
   9032 format(' INITIAL restraints/deviations/energy contributions:')
   9033 format(' Initial restraints echoed to file: ',a)
   9040 format('******',/, &
         1x,a4,'(',i5,')-',a4,'(',i5,')',t53,'NSTEP1=',i6, &
         ' NSTEP2=',i6)
   9041 format('******',/,1x,a4,'(',i5,')-',a4,'(',i5,')-',a4,'(',i5,')', &
         t53,'NSTEP1=',i6,' NSTEP2=',i6)
   9042 format('******',/,1x,a4,'(',i5,')-',a4,'(',i5,')-',a4,'(',i5,')-', &
         a4,'(',i5,')',t53,'NSTEP1=',i6,' NSTEP2=',i6)
   9043 format(1x,'NINC=',i6)
   9045 format('R1 =',f8.3,' R2 =',f8.3,' R3 =',f8.3,' R4 =',f8.3, &
         ' RK2 =',f8.3,' RK3 = ',f8.3)
   9046 format('R1A=',f8.3,' R2A=',f8.3,' R3A=',f8.3,' R4A=',f8.3, &
         ' RK2A=',f8.3,' RK3A= ',f8.3)

   9047 format(1x,'Rcurr: ',f8.3,'  Rcurr-(R2+R3)/2: ',f8.3, &
         '  MIN(Rcurr-R2,Rcurr-R3): ',f8.3)
   9048 format('******',/,1x,a4,'(',i5,')-',a4,'(',i5,')-',a4,'(',i5,')-', &
         a4,'(',i5,')-',/,1x,a4,'(',i5,')',t53,'NSTEP1=',i6,' NSTEP2=',i6)
   9049 format('******',/,1x,a4,'(',i5,')-',a4,'(',i5,')-',a4,'(',i5,')-', &
         a4,'(',i5,')-',/,1x,a4,'(',i5,')-',a4,'(',i5,')-',a4,'(',i5,')-', &
         a4,'(',i5,')',t53,'NSTEP1=',i6,' NSTEP2=',i6)
   ! 9050 format('******',/,1x,a4,'(',i5,')-',a4,'(',i5,')-',a4,'(',i5,')-', &
   !       a4,'(',i5,')-',/,1x,a4,'(',i5,')-',a4,'(',i5,')',t53,'NSTEP1=',i6,' NSTEP2=',i6)
   9065 format(' Atom ranges defining Center of Mass Group in the first position: ')
   9066 format(5(1x,i5,' -> ',i5,'/'))
   9067 format(' Atom ranges defining Center of Mass Group in the second position: ')
   9068 format(1x,a79)
   9069 format(' Requested file redirections:')
   9070 format(t3,a,t13,'= ',a)
   9071 format(t3,'No valid redirection requests found')
   ! 9072 format('In Jarzynsky runs, there must only one restraint, stopping program')
   9073 format('Error: Invalid residue number or atom name in restraint.')
   9074 format('Error: Atom groupings not supported for j-coupling restraints.')
   9075 format('Error: Improper restraint specified.')
   9076 format(' Atom ranges defining the Center of Mass Group in the third position: ')
   9077 format(' Atom ranges defining the Center of Mass Group in the ',i1,'th position: ')
   9078 format('Error: Atom ',i1,' and atom ',i1,&
               ' cannot be the same in a plane restraint.')
   9079 format('Warning: Both iat(1) and "restraint" variables specified.', &
               'Using only iat() and not parsing the contents of "restraint."')
   ! 9080 format('Warning: Both rk2/rk3/rk2a/rk3a and k0/k0a defined.', &
   !            'Using only original-style restraint definition and not r0/k0/r0a/k0a."')
end subroutine nmrred 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmrsht here]
subroutine nmrsht(x,nmrnum,rimass,ntb,nstep,nmrat,resttype,nmrcom,nres, &
      ipres,ishrtb,iwtstp,wtnrg,rk2nmr,rk3nmr,nmrst, &
      nmrshb,wtls,nmrfty,igravt,ialtdis, &
      bave,bave0,nave, &
      nexact,ipower,tauave,navint,ravein,dt,iout, &
      iflag)
   
   
   ! Subroutine NMR SHorT range interactions identification.
   
   ! This routine determines which restraints should be classified as
   ! "short-range", based on the information given by the user on a SHORT
   ! command line. The NMRSHB(I) array is set to indicate the short-range
   ! interactions:
   !       NMRSHB(I) = 0 : Restraint I is not short-range
   !                 = 1 : Restraint I is short-range based on residue # criterion
   !                 = 2 : Restraint I is short-range based on distance criterion
   !                 = 3 : Restraint I is short range based on both residue # and
   !                       distance criteria.
   
   ! Author: David A. Pearlman
   ! Date: 11/89
   
   ! Modified by Matthew Seetin
   ! Date: 9/2007
   ! This routine now supports short range interactions for angle and torsion restraints
   ! that involve atom groupings.  Plane-plane and plane-point restraints are not
   ! allowed to be short range.
   
   ! INPUT:
   
   !  X(I)       : Coordinate array
   !  NMRNUM     : The number of NMR restraints
   !  RIMASS(I)  : The inverse mass for atom I.
   !  NTB        : Periodic boundary conditions flag
   !  NSTEP      : The step/iteration number in the main calling program
   !  NMRAT(16,I) : The 2-8 atoms defining restraint I.
   !                I = 9,16 are < 0 if the atom is actually an atom group
   !  NMRCOM(2,I) : The ranges of atoms defining the center-of-mass group,
   !                for applicable distance restraints
   !  NRES       : The number of residues in the system.
   !  IPRES(I)   : IPRES(I) is the first atom in residue I. IPRES(I+1)-1
   !               is the last atom in residue I.
   !  ISHRTB     : If ISHRTB >0, short-range interactions have been defined;
   !               Atom-atom interactions are defined as short-range if either:
   
   !               i) IWTSTP(1,ISHRTB) <= ABS(residue_diff) <= IWTSTP(2,ISHRTB)
   
   !                  (where residue_diff is the difference in the numbers of
   !                   the residues containing atoms i and j)
   
   !               ii) WTNRG(1,ISHRTB) <= dist(i-j) <= WTNRG(2,ISHRTB)
   
   !                   (where dist(i-j) is the distance between atoms i and j).
   
   !               A restraint is considered short-range if all the bonds
   !               contributing to the restraint are short-range.
   
   !               The short-range distance classification is re-evaluated every
   !               IWTSTP(3,ISHRTB) steps, if IWTSTP(3,ISHRTB) is > 0.
   
   ! IWTSTP(3,I)
   ! WTNRG(2,I) : See above.
   ! RK2NMR(2,I) : The slope and intercept for k2.
   ! RK3NMR(2,I) : The slope and intercept for k3.
   ! NMRST(3,I)  : If NMRST(3,I) > 0, k2 & k3 are linearly interpolated with
   !               step number. If NMRST(3,I) < 0, then k2 & k3 are generated
   !               by multiplication by a constant factor.
   ! WTLS(2)     : WTLS(1) is current weight of short-range interactions (in MODWT)
   !               WTLS(2) is current weight of non "  "    interactions (in MODWT)
   !  NMRFTY(I)  : Flag for each restraint indicating what functional form should
   !               be used for that restraint. If time-averaged restraints
   !               are being applied to a particular class of internal and
   !               NMRFTY(I)=1, then an instantaneous restraint will be applied.
   !    BAVE(I)
   !   BAVE0(I) : For time-averaged restraints. BAVE(I) contains the time-averaged
   !              values of the distances, using exponential decay with
   !              a time-constant of TAUAVE(1). BAVE0(I) is the same thing,
   !              but without exponential decay (real average).
   !     NAVE(3) : NAVE(I) > 0 for time-averaged restraints.
   !   NEXACT(3) : Not currently used.
   !   IPOWER(3) : For time-averaged restraints. Gives the exponent used in
   !               developing the time average.
   !   TAUAVE(3) : For time-averaged restraints. Gives the time constant used
   !               in the exponential decay multiplier.
   !   NAVINT(3) : For time-averaged restraints. Averaged values only updated
   !               every NAVINT(I) steps.
   !   RAVEIN(3) : For time-averaged restraints. On the first step, the value
   !               returned is the current value+RAVEIN(I) (I=1 for bonds).
   !          DT : Integration timestep (ps). Used for time-averaged restraints.
   ! IOUT       : Unit for informational (error) prints.
   ! IFLAG      : = 0 First call to this routine. Determine short-range
   !                  restraints meeting either criterion (i) or (ii) above.
   !                  No check against IWTSTP(3,ISHRTB) is done.
   !              > 0 then only re-evaluate short-range restraints for
   !                  criterion (ii) above. A check against IWTSTP(3,ISHRTB)
   !                  is performed.
   
   ! OUTPUT:
   
   ! NMRSHB(I) : Short-range indicator flag, set for each restraint as
   !             described above.
   
   ! Notes: 1) When a distance restraint has been defined where one or
   !           both ends of the restraints refer to a center-of-mass coordinate,
   !           this routine assumes the atom ranges have been packed in NMRCOM
   !           in increasing sequential order, and that residue number increases
   !           monotonically with atom number. Thus, only the residues
   !           containing the first and last atoms of the center-of-mass definition
   !           need to be examined for the residue check (i).
   
   !         2) If time-averaged restraints are being monitored, the time
   !            averaged values of distances will be used in determining
   !            short-range interactions for bond restraints. HOWEVER, the
   !            short-range status of angle and torsion restraints will be
   !            determined using instantaneous values of the bond lengths
   !            comprising these internals.
   
   implicit none
   integer:: i, ialtdis, iat, iave, iflag, igravt, iout, ipower, &
        ipres, ires1, ires2, irs1, irs2, irs3, irs4, irsdif, ishort, &
        ishrtb, iwtstp, j, nave, navint, nexact, nmrat, nmrcom, nmrfty, &
        nmrnum, nmrshb, nmrst, nres, nstep, ntb
   _REAL_ :: bave, bave0, dt, dum, duma, ravein, rimass, rk2nmr, &
        rk3nmr, rlow, rmstot, rmult1, rmult2, rr, rr1, rr2, rr3, rup, &
        small, tauave, wtls, wtnrg, x, xcom
   parameter (small = 1.0d-12)
#  include "nmr.h"
   dimension x(*),nmrat(16,*),rimass(*),nmrcom(2,*),ipres(*),wtls(2)
   dimension iwtstp(3,*),wtnrg(2,*),nmrshb(*),rk2nmr(2,*),rk3nmr(2,*)
   dimension nmrst(3,*),nmrfty(*),igravt(*),ialtdis(*)
   dimension bave(*),bave0(*),nave(6),duma(2)
   dimension nexact(6),ipower(6),tauave(6),navint(6),ravein(6)
   dimension iat(8),xcom(24),rmstot(8)
   integer resttype(*)
   logical usecom
   integer atmarr
   dimension atmarr(16)
   integer jjj,resttypetmp

   _REAL_ ricmmass(8)
   
   if (iflag == 0) then
      do i=1,nmrnum
         nmrshb(i) = 0
      end do
   end if
   
   ! Return if no short-range interactions defined
   
   if (ishrtb == 0) return
   
   ! If IFLAG.NE.0, re-determine short-range interactions only every
   ! IWTSTP(3,ISHRTB) steps:
   
   if (iflag /= 0) then
      if (iwtstp(3,ishrtb) == 0) return
      if (mod(nstep,iwtstp(3,ishrtb)) /= 0 .or. nstep == 0) return
   end if
   
   ! Check for residue range defined short-range interactions if IFLAG=0:
   
   if (iflag == 0) then
      do i=1,nmrnum
         usecom = .false.
         do jjj=9,16
           usecom = usecom .or. nmrat(jjj,i) < 0
         end do
         resttypetmp = resttype(i)
         select case (resttypetmp)
           case (1)
            
            ! Center of mass-type distance restraints:
            
            if (usecom) then
               do j=1,2
                  if (nmrat(j,i) < 0) then
                     atmarr(2*(j-1)+1) = nmrcom(1,-nmrat(j,i))
                     atmarr(2*(j-1)+2) = nmrcom(2,-nmrat(j+8,i))
                  else
                     atmarr(2*(j-1)+1) = nmrat(j,i)
                     atmarr(2*(j-1)+2) = nmrat(j,i)
                  end if
               end do
               atmarr(1) = min(atmarr(1),atmarr(3))
               atmarr(2) = max(atmarr(2),atmarr(4))
               call at2res(atmarr(1),ipres,nres,ires1,1,iout)
               call at2res(atmarr(2),ipres,nres,ires2,1,iout)
            else
               
               ! Regular distance restraints:
               
               call at2res(nmrat(1,i),ipres,nres,ires1,1,iout)
               call at2res(nmrat(2,i),ipres,nres,ires2,1,iout)
            end if
            
            ! Angle restraints:
            
          case (2)
            if (usecom) then
              do j=1,3
                if (nmrat(j,i) < 0) then
                  atmarr(j) = nmrcom(1,-nmrat(j,i))
                  atmarr(j+8) = nmrcom(2,-nmrat(j+8,i))
                else
                  atmarr(j) = nmrat(j,i)
                  atmarr(j+8) = nmrat(j,i)
                end if
              end do
              atmarr(1) = min(atmarr(1),atmarr(2),atmarr(3))
              atmarr(2) = max(atmarr(9),atmarr(10),atmarr(11))
              call at2res(atmarr(1),ipres,nres,ires1,1,iout)
              call at2res(atmarr(2),ipres,nres,ires2,1,iout)
              
            else
              call at2res(nmrat(1,i),ipres,nres,irs1,1,iout)
              call at2res(nmrat(2,i),ipres,nres,irs2,1,iout)
              call at2res(nmrat(3,i),ipres,nres,irs3,1,iout)
              ires1 = min(irs1,irs2,irs3)
              ires2 = max(irs1,irs2,irs3)
            end if
            ! Torsional restraints:
            
         case (3)
            if (usecom) then
              do j=1,4
                if (nmrat(j,i) < 0) then
                  atmarr(j) = nmrcom(1,-nmrat(j,i))
                  atmarr(j+8) = nmrcom(2,-nmrat(j+8,i))
                else
                  atmarr(j) = nmrat(j,i)
                  atmarr(j+8) = nmrat(j,i)
                end if
              end do
              atmarr(1) = min(atmarr(1),atmarr(2),atmarr(3),atmarr(4))
              atmarr(2) = max(atmarr(9),atmarr(10),atmarr(11),atmarr(12))
              call at2res(atmarr(1),ipres,nres,ires1,1,iout)
              call at2res(atmarr(2),ipres,nres,ires2,1,iout)
              
            else
              call at2res(nmrat(1,i),ipres,nres,irs1,1,iout)
              call at2res(nmrat(2,i),ipres,nres,irs2,1,iout)
              call at2res(nmrat(3,i),ipres,nres,irs3,1,iout)
              call at2res(nmrat(4,i),ipres,nres,irs4,1,iout)
              ires1 = min(irs1,irs2)
              ires1 = min(ires1,irs3)
              ires1 = min(ires1,irs4)
              ires2 = max(irs1,irs2)
              ires2 = max(ires2,irs3)
              ires2 = max(ires2,irs4)
            end if
         case default
            ! Short range not supported for planar restraints
            nmrshb(i) = 0
            cycle
         end select
         
         ! Do the residue check:
         
         irsdif = abs(ires2-ires1)
         if (iwtstp(1,ishrtb) <= irsdif .and. &
               irsdif <= iwtstp(2,ishrtb)) nmrshb(i) = 1
      end do
   end if  ! iflag == 0
   
   ! Now do the distance check:
   
   ! Define conversion factors required to account for weight changes
   ! if short-range classification of a restraint changes:
   
   if (wtls(2) < small) then
      rmult1 = 0.0d0
   else
      rmult1 = wtls(1)/wtls(2)
   end if
   if (wtls(1) < small) then
      rmult2 = 0.0d0
   else
      rmult2 = wtls(2)/wtls(1)
   end if
   
   rlow = wtnrg(1,ishrtb)
   rup = wtnrg(2,ishrtb)
   do i=1,nmrnum
      ishort = 0
      iave = 0
      usecom = .false.
      do jjj=9,16
        usecom = usecom .or. nmrat(jjj,i) < 0
      end do
      resttypetmp = resttype(i)
      select case (resttypetmp)
         
       case(1)  ! Center of mass-type distance restraints:
         
         if (nave(1) > 0 .and. nmrfty(i) == 0) iave = 1
         if ((nmrat(1,i) < 0 .or. nmrat(2,i) < 0) .and. &
               igravt(i) == 0) then
            call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
            call disnrg(xcom,duma,duma,rr,0,3,dum,dum,dum,dum, &
                  dum,dum,dum,ntb,bave(2*i-1),bave0(i), &
                  nave(1),nexact(1),ipower(1),tauave(1), &
                  ravein(1),dt,navint(1),iave,0,0,1,ialtdis(i))
            if (rlow <= rr.and.rr <= rup) ishort = 1
            
            ! r**-6 averaged distance interaction restraints:
            
         else if ((nmrat(1,i) < 0 .or. nmrat(2,i) < 0) .and. &
               igravt(i) == 1) then
            call r6ave(x,nmrat,nmrcom,rr,i)
            
         else
            
            ! Regular distance restraints:
            
            call disnrg(x,duma,duma,rr,nmrat(1,i),nmrat(2,i),dum, &
                  dum,dum,dum,dum,dum,dum,ntb,bave(2*i-1), &
                  bave0(i),nave(1),nexact(1),ipower(1), &
                  tauave(1),ravein(1),dt,navint(1),iave, &
                  0,0,1,ialtdis(i))
            if (rlow <= rr.and.rr <= rup) ishort = 1
         end if
         
         ! Angle restraints:
         
      case (2)
         iave = 0
         if (usecom) then
           call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
           call disnrg(xcom,duma,duma,rr1,0,3,dum, &
               dum,dum,dum,dum,dum,dum,ntb,duma,duma, &
               nave(2),nexact(2),ipower(2),tauave(2), &
               ravein(2),dt,0,0,0,0,1,ialtdis(i))
           call disnrg(xcom,duma,duma,rr2,3,6,dum, &
               dum,dum,dum,dum,dum,dum,ntb,duma,duma, &
               nave(2),nexact(2),ipower(2),tauave(2), &
               ravein(2),dt,0,0,0,0,1,ialtdis(i))
         else
           call disnrg(x,duma,duma,rr1,nmrat(1,i),nmrat(2,i),dum, &
               dum,dum,dum,dum,dum,dum,ntb,duma,duma, &
               nave(2),nexact(2),ipower(2),tauave(2), &
               ravein(2),dt,0,0,0,0,1,ialtdis(i))
           call disnrg(x,duma,duma,rr2,nmrat(2,i),nmrat(3,i),dum, &
               dum,dum,dum,dum,dum,dum,ntb,duma,duma, &
               nave(2),nexact(2),ipower(2),tauave(2), &
               ravein(2),dt,0,0,0,0,1,ialtdis(i))
         end if
         if ((rlow <= rr1.and.rr1 <= rup) .and. &
               (rlow <= rr2.and.rr2 <= rup)) ishort = 1
         
         ! Torsional restraints:
         
      case (3)
         iave = 0
         if (usecom) then
           call nmrcms(x,xcom,nmrat,nmrcom,rmstot,rimass,ricmmass,i)
           call disnrg(x,duma,duma,rr1,0,3,dum, &
               dum,dum,dum,dum,dum,dum,ntb,duma,duma, &
               nave(3),nexact(3),ipower(3),tauave(3), &
               ravein(3),dt,0,iave,0,0,1,ialtdis(i))
           call disnrg(x,duma,duma,rr2,3,6,dum, &
               dum,dum,dum,dum,dum,dum,ntb,duma,duma, &
               nave(3),nexact(3),ipower(3),tauave(3), &
               ravein(3),dt,0,iave,0,0,1,ialtdis(i))
           call disnrg(x,duma,duma,rr3,6,9,dum, &
               dum,dum,dum,dum,dum,dum,ntb,duma,duma, &
               nave(3),nexact(3),ipower(3),tauave(3), &
               ravein(3),dt,0,iave,0,0,1,ialtdis(i))
         else
           call disnrg(x,duma,duma,rr1,nmrat(1,i),nmrat(2,i),dum, &
               dum,dum,dum,dum,dum,dum,ntb,duma,duma, &
               nave(3),nexact(3),ipower(3),tauave(3), &
               ravein(3),dt,0,iave,0,0,1,ialtdis(i))
           call disnrg(x,duma,duma,rr2,nmrat(2,i),nmrat(3,i),dum, &
               dum,dum,dum,dum,dum,dum,ntb,duma,duma, &
               nave(3),nexact(3),ipower(3),tauave(3), &
               ravein(3),dt,0,iave,0,0,1,ialtdis(i))
           call disnrg(x,duma,duma,rr3,nmrat(3,i),nmrat(4,i),dum, &
               dum,dum,dum,dum,dum,dum,ntb,duma,duma, &
               nave(3),nexact(3),ipower(3),tauave(3), &
               ravein(3),dt,0,iave,0,0,1,ialtdis(i))
         end if
         if ((rlow <= rr1.and.rr1 <= rup) .and. &
               (rlow <= rr2.and.rr2 <= rup) .and. &
               (rlow <= rr3.and.rr3 <= rup)) ishort = 1
      case default 
        ! Short range not supported for other restraints
        ishort = 0
      end select  
      
      ! If the short-range classification of the restraint has changed, change
      ! NMRSHB(I) to reflect this fact. Also, modify K2 and K3 so that they
      ! are multiplied by the weight factor appropriate for their new classification:
      
      if (nmrshb(i) <= 0 .and. ishort == 1) then
         nmrshb(i) = 2
         rk2nmr(2,i) = rk2nmr(2,i)*rmult1
         rk3nmr(2,i) = rk3nmr(2,i)*rmult1
         if (nmrst(3,i) > 0) then
            rk2nmr(1,i) = rk2nmr(1,i)*rmult1
            rk3nmr(1,i) = rk3nmr(1,i)*rmult1
         end if
      else if (nmrshb(i) == 2 .and. ishort == 0) then
         nmrshb(i) = 0
         rk2nmr(2,i) = rk2nmr(2,i)*rmult2
         rk3nmr(2,i) = rk3nmr(2,i)*rmult2
         if (nmrst(3,i) > 0) then
            rk2nmr(1,i) = rk2nmr(1,i)*rmult2
            rk3nmr(1,i) = rk3nmr(1,i)*rmult2
         end if
      else if (nmrshb(i) == 1 .and. ishort == 1) then
         nmrshb(i) = 3
      else if (nmrshb(i) == 3 .and. ishort == 0) then
         nmrshb(i) = 1
      end if
   end do
   
   return
end subroutine nmrsht 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PARENtheses FINDer locates the parentheses that denote a plane or COM
subroutine parenfind(restraint,k,restsize,openparen,closparen,iout)

  ! Subroutine PARENtheses FINDer
  
  ! This routine figures out the opening and closing parethesis for a 
  ! plane or COM grouping in the "restraint" variable.
  
  ! Authoer: Matthew G. Seetin
  ! Date: 10/2007
  
  !---------------------------
  ! INPUT VARIABLES
  !---------------------------
  
  ! restraint: a character variable of max 256 charaters that contains the
  !             natural-language input restraint
  ! restsize : equal to LEN_TRIM(restraint)
  ! openparen: the opening parenthesis of the group
  ! closparen: the closing parenthesis of the group
  !         k: the character being read in restraint
  !      iout: Unit number of the out file
  implicit none
  integer:: k
  CHARACTER(256) :: restraint
  INTEGER i, restsize,openparen,closparen,iout
  CHARACTER(1) :: opentype, clostype

  
  do i=k,restsize
    if(restraint(i:i) == ' ') then
      cycle
    else
      openparen = i
      closparen = i+1
      exit
    end if
  end do

  select case (restraint(i:i))
    case ('(')
      opentype = '('
      clostype = ')'
    case ('[')
      opentype = '['
      clostype = ']'
    case ('{')
      opentype = '{'
      clostype = '}'
    case default
      write(iout,9009)
      write(iout,9999) restraint(1:restsize)
      call mexit(iout, 1)
  end select
  
  do
    if (index(restraint(closparen:restsize), clostype) > &
    index(restraint(closparen:restsize), opentype) .and. &
    index(restraint(closparen:restsize), opentype) /= 0) then
      closparen = closparen + index(restraint(closparen:restsize), ')') 
    else if (index(restraint(closparen:restsize), clostype) == 0) then
      write(iout,9008)
      write(iout,9999) restraint(1:restsize)
      call mexit(iout, 1)
    else
      closparen = closparen + index(restraint(closparen:restsize), clostype) - 1
      exit
    end if
  end do
  
  return
  
  9008 format('Error: No closing parenthesis for plane or COM grouping in restraint.')
  9009 format('Error: No parentheses specified to define plane or COM grouping in restraint.')
  9999 format('restraint = "',(a),'"')

end subroutine parenfind

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PARSe COM grouping; read restraint char var, incorporate into igr
subroutine parscom(restraint,iat,igrarr,i,closparen,iatidx,restsize,name,ipres,nres,iout)

  ! Subroutine PARSe COM grouping
  
  ! This routine figures out the opening and closing parethesis for a 
  ! plane or COM grouping in the "restraint" variable.
  
  ! Authoer: Matthew G. Seetin
  ! Date: 10/2007
  
  !---------------------------
  ! INPUT VARIABLES
  !---------------------------
  
  ! restraint: a character variable of max 256 contains the
  !             natural-language input restraint
  ! restsize : equal to LEN_TRIM(restraint)
  ! closparen: the closing parenthesis of the group
  !         i: the character being read in restraint
  !igrarr(200,8): the arry that holds all the atoms for each COM grouping
  !    iat(8): Array of all atom numbers in the restraint
  !   iatidx : index of current atom being read
  !     name : Array of the names of all the atoms in the system.  It may be empty
  !            if this routine is called by restal
  !     ipres: Array of all the residue numbers in the system.  May also be empty if
  !            if this routine is called by restal
  !      nres: The number of residues in the system
  !      iout: Unit number of the out file
  implicit none
  integer:: i, iat, igrarr, iout, ipres, nres
  character(256) :: restraint
  integer comatms,numlen,namlen,igrpass,j,closparen,iatidx,restsize
  dimension igrarr(200,8), iat(8),ipres(*),igrpass(8)
  character(4) name(*)
  integer dumiatidx
  
  
  comatms = 0
  
  do j=1,8
    igrpass(j) = 0
  end do
  
  do 
    dumiatidx=1
    numlen = 0
    namlen = 0
    if(i>closparen-1) exit
    if (comatms > 200) then
      ! Print warning message
      exit
    end if
    if (restraint(i:i) == ' ' .or. restraint(i:i) == ',' .or. restraint(i:i) == '(' &
          .or. restraint(i:i) == ')' .or. restraint(i:i) == '['.or. restraint(i:i) == ']' &
          .or. restraint(i:i) == '{'.or. restraint(i:i) == '}') then
      i = i + 1
      cycle
    end if
    if (iachar(restraint(i:i)) >= 48 .AND. iachar(restraint(i:i)) <= 57) then ! Explicit atom number specified
      numlen = 1
      do 
        if (iachar(restraint(i+numlen:i+numlen)) >= 48 .AND. iachar(restraint(i+numlen:i+numlen)) <= 57) then
          numlen = numlen + 1
        else
          exit
        end if
      end do
      comatms = comatms+1
      read (restraint(i:i+numlen-1), *) igrarr(comatms,iatidx)
      i = i + numlen
    else if (restraint(i:i) == ':') then  !ptraj-style atom mask
      comatms = comatms + 1      
      call readmask(restraint,i,igrpass,dumiatidx,name,ipres,nres,iout)

      igrarr(comatms,iatidx) = igrpass(1)
        
    else
      write(iout,9006)
      write(iout,9999) restraint(1:restsize)
      call mexit(iout, 1)
    end if
    i = i + 1   
  end do             ! Done parsing COM grouping
  if (comatms == 0) then
    write(iout,9010)
    write(iout,9999) restraint(1:restsize)
    call mexit(iout, 1)
  else if (comatms == 1) then
    iat(iatidx) = igrarr(comatms,iatidx)
    iatidx = iatidx + 1
  else
    iat(iatidx) = -iatidx
    iatidx = iatidx + 1
  end if
  
  i = closparen + 1
  
  return
  
  9006 format('Error: Invalid atom or grouping specified in restraint.')
  9010 format('Error: No atoms specified within COM grouping in restraint.')
  9999 format('restraint = "',(a),'"')
   
end subroutine parscom

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PARSe GENeralzied coordinate restraint reads a coordinate subgrouping
subroutine parsgen(restraint,iat,igrarr,subgrpatoms,k,closparen,iatidx,comidx,currcom,&
                   restsize,name,ipres,nres,iout)
                   
  ! Subroutine GENeralized coordinate
  
  ! This subroutine parses the a distance, angle, or torsion (latter 2 not 
  ! yet active in AMBER) subgrouping from natural language restraint input
  
  ! Author: Matthew G. Seetin
  ! Date: 11/2007
   
  !---------------------------
  ! INPUT VARIABLES
  !---------------------------
  
  ! restraint : a character variable of max 200 charaters that contains the
  !             natural-language input restraint
  ! igrarr(200,8): The atoms in com grouping 1-8.  If igrarr(1,i) /= 0, iat(i) < 0.  
  !     name : Array of the names of all the atoms in the system.  It may be empty
  !            if this routine is called by restal
  !     ipres: Array of all the residue numbers in the system.  May also be empty if
  !            if this routine is called by restal
  !      nres: The number of residues in the system
  !      iout: Unit number of the out file
  !   iatidx : Index as to where in iat the next atom number read will be stored
  ! comidx(8): Locations of the "com" keyword in restraint
  !   currcom: Current COM grouping in use
  !  numatoms: Number of atoms expected in the grouping, two for distances, three for
  !            angles, and four for torsions.
  !         k: The index of the position of restraint that's currently being read
  ! closparen: The location of the parethesis that closes the subgrouping
  
  !---------------------------
  ! OUTPUT VARIABLES
  !---------------------------
  
  !    iat(8) : the atoms that make up the restraint.
  ! igrarr(200,8): The atoms in com grouping 1-8.  If igrarr(1,i) /= 0, iat(i) < 0.  
  implicit none
   integer:: iat, iattmp, igrarr, igrpass, iout, ipres, k, namlen, nres, numlen
   _REAL_ :: subgrpatms
  
  CHARACTER(256) :: restraint
  INTEGER subgrpatoms,numatoms,currcom,closparen,restsize,iatidx,comidx,openparen
  DIMENSION igrarr(200,8), iat(8),ipres(*),igrpass(200),comidx(8)
  CHARACTER(4) name(*)
   

  numatoms = 0
  
  do 
    if (k >= closparen) exit
    numlen = 0
    namlen = 0
    if (restraint(k:k) == ' ' .or. restraint(k:k) == ',' .or. restraint(k:k) == '(' &
        .or. restraint(k:k) == ')' .or. restraint(k:k) == '['.or. restraint(k:k) == ']' &
        .or. restraint(k:k) == '{'.or. restraint(k:k) == '}') then
        k = k + 1
        cycle
    end if
    if (numatoms > subgrpatoms) then
         write(iout,9000)
         write(iout,9001) subgrpatms
         write(iout,9999) restraint(1:restsize)
         exit
    end if
      
    if ((k < comidx(min(currcom,8))) .or. (comidx(min(currcom,8)) == 0)) then
       if (iachar(restraint(k:k)) >= 48 .AND. iachar(restraint(k:k)) <= 57) then ! User specifies explicit atom number
          do 
            if (iachar(restraint(k+numlen:k+numlen)) >= 48 .AND. iachar(restraint(k+numlen:k+numlen)) <= 57) then
              numlen = numlen + 1
            else
              exit
            end if
          end do
          read (restraint(k:k+numlen-1), *) iattmp
          iat(iatidx) = iattmp
          iatidx = iatidx + 1
          numatoms = numatoms + 1
          k = k + numlen
       else if (restraint(k:k) == ':') then  ! User specifies ptraj-style atom mask
         
         call readmask(restraint,k,iat,iatidx,name,ipres,nres,iout)
    
        else if (restraint(k:k) == '-') then   ! User specifies group of atoms in old way using igr
          if (iachar(restraint(k+1:k+1)) >= 48 .AND. iachar(restraint(k+1:k+1)) <= 57) then
            do
              if (iachar(restraint(k+numlen+1:k+numlen+1)) >= 48 .AND. iachar(restraint(k+numlen+1:k+numlen+1)) <= 57) then
                numlen = numlen + 1
              else
                exit
              end if
            end do
            read(restraint(k:k+numlen), *) iattmp
            iat(iatidx) = iattmp
            
            if (igrarr(1,iatidx) == 0) then
              write(iout,9004)
              write(iout,9999) restraint(1:LEN_TRIM(restraint))
              call mexit(iout, 1)
            end if
            
            iatidx = iatidx + 1
            numatoms = numatoms + 1
            k = k + numlen + 1
          else
            write(iout,9005)
            write(iout,9999) restraint(1:restsize)
            call mexit(iout, 1)
          end if
        else
          write(iout,9006)
          write(iout,9999) restraint(1:restsize)
          call mexit(iout, 1)
        end if
 
    else if (k == comidx(min(currcom,8))) then
      
        ! Find the parenthetical boundaries of the COM grouping
        k = k + 3
        call parenfind(restraint,k,restsize,openparen,closparen,iout)
        
        k = openparen + 1
        
        ! Parse a user-defined COM grouping
        
        call parscom(restraint,iat,igrarr,k,closparen,iatidx,restsize,name,ipres,nres,iout)
        
        k = closparen + 1
        currcom = currcom + 1
    else
         write(iout,9006)
         write(iout,9999) restraint(1:restsize)
         call mexit(iout, 1)
    end if
  end do
  
  9000 format('Warning: More characters in restraint beyond final atom necessary')
  9001 format('Using only first ',i1,' atoms.')
  ! 9002 format('Error: No valid restraint type specified.')
  ! 9003 format('Error: Invalid residue number or atom name in restraint.')
  9004 format('Error: Manual atom grouping specified, but none defined.')
  9005 format('Error: Invalid group number in restraint.')
  9006 format('Error: Invalid atom or grouping specified in restraint.')
  9999 format('restraint = "',(a),'"')
    
end subroutine parsgen

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ PARSe RESTraint: read restraint char var, incorporate into old vars
subroutine parsrest(restraint,iat,igrarr,rstwt,name,ipres,nres,iout)

  ! Subroutine PARSe RESTraint
  
  ! This subroutine parses the character variable restraint and places
  ! the extracted information into the variables IAT and IGR.
  
  ! Author: Matthew G. Seetin
  ! Date: 9/2007
  
  ! This function enables "natural language" input for restraints as an
  ! alternative to AMBER's traditional restraint input
  
  ! The general form for a restraint of this type is as follows
  ! restraint = "restrainttype(atom1,atom2,com(atom3,...,atomN),...)
  ! restrainttype will be either distance, angle, torsion, or coordinate
  ! "com" denotes a center of mass grouping
  ! Addtionally, the "plane" keyword may be used within the angle restraint
  ! type to invoke a new style of restraint, a restraint on the angle between
  ! a plane's normal vector and either another plane's normal vector, or the 
  ! vector between the planes com and another atom or com.
  
  ! Either (),{}, or [] may be used to partion off the restraint itself or a set 
  ! of atoms designated as a plane or as a com
  
  ! Atoms can be specified using either the explicit atom number or using a 
  ! notation similar to ptraj's atom masks, i.e., :residue_number@atom_name
  ! User may specify a negative atom number, in which case the progam will look for 
  ! a user-defined group of atoms in igr1-igr8.
  
  ! Atoms may be separated by either commas or spaces
  
  !---------------------------
  ! INPUT VARIABLES
  !---------------------------
  
  ! restraint : a character variable of max 200 charaters that contains the
  !             natural-language input restraint
  ! igrarr(200,8): The atoms in com grouping 1-8.  If igrarr(1,i) /= 0, iat(i) < 0.  
  !     name : Array of the names of all the atoms in the system.  It may be empty
  !            if this routine is called by restal
  !     ipres: Array of all the residue numbers in the system.  May also be empty if
  !            if this routine is called by restal
  !      nres: The number of residues in the system
  !      iout: Unit number of the out file
  
  !---------------------------
  ! OUTPUT VARIABLES
  !---------------------------
  
  !    iat(8) : the atoms that make up the restraint.
  ! igrarr(200,8): The atoms in com grouping 1-8.  If igrarr(1,i) /= 0, iat(i) < 0.  
  ! rstwt(4)  : The weights used in a generalized distance coordinate restraint
  
  !---------------------------
  ! INTERNAL VARIABLES
  !---------------------------
  
  !   iatidx : Index as to where in iat the next atom number read will be stored
  ! plnidx(2): Locations of the "plane" keyword in restraint
  ! comidx(8): Locations of the "com" keyword in restraint
  !distidx(4): Locations of the "distance" keyword in restraint (used only for gen. dist. coord.)
  
  !reststart : The position after the first open parenthesis.
  ! restatms : The number of atoms (or com groups) that is necessary for the specified 
  !            restraint type.
  !            distance:    restatms = 2
  !            angle:       restatms = 3
  !            torsion:     restatms = 4
  !            plane-point: restatms = 5
  !            plane-plane: restatms = 8
  !            generalized distance: restatms will be 2, 4, 6, or 8
  !openparen
  !closparen 
  !openparensub
  !closparensub : Positions of opening and closing parentheses when dealing with
  !               plane and com keywords
  ! plnatms : The number of atoms or com groups read in during the corse of dealing 
  !           with a "plane" keyword.  Should always be 4 when done with the plane.
  !           Otherwise execution will be terminated
  ! comatms : The nubmer of atoms in a com group. It should be that 2 <= comatms <= 200.
  !           If comatms = 1, then that atom is treated like a single atom.
  !           If comatms = 0, execution will be terminated
  !  numlen : The number of digits in either an atom number or a residue number
  !  namlen : The number of characters in a specified atom name
  ! nonblnk : The position of the first non-blank character in restraint
  ! restsize: LEN(restraint)
  ! itworked: Boolean for testing that everything is OK.
  ! emptytmp: Boolean for verifying that igrtmp is empty

  implicit none
  integer:: iat, igrarr, iout, ipres, nres
  INTEGER :: nonblnk,numlen,namlen,iatidx,reststart,restatms,closparen,openparen,plnatms, &
             closparensub,openparensub,restsize,currcom,currpln,currdist,i,j,k,m,n,iattmp
  INTEGER, DIMENSION(8) :: comidx
  INTEGER, DIMENSION(2) :: plnidx
  INTEGER, DIMENSION(4) :: distidx
  LOGICAL itworked, emptytmp, givemeane, givemeadot, itsnegative
  DIMENSION iat(8), igrarr(200,8), ipres(*)
  _REAL_ rstwt(4)
  CHARACTER(4) name(*)
  CHARACTER(256) :: restraint
  CHARACTER(32) resttmp

  comidx = (/ 0, 0, 0, 0, 0, 0, 0, 0 /)
  do i=1,8
    iat(i) = 0
  end do
  do i=1,4
    rstwt(i) = 0.0d0
  end do
  plnidx = (/ 0, 0 /)
  nonblnk = 1
  numlen = 0
  namlen = 0
  iatidx = 1
  reststart = 0
  restatms = 0
  closparen = 0
  openparen = 0
  plnatms =0
  closparensub = 0
  openparensub = 0
  restsize = 0
  iattmp=0
  currcom=1
  currpln=1
  currdist=1
  k=0

  restsize = LEN_TRIM(restraint)

  do i=1,restsize
    if (restraint(i:i) == ' ') Then
      nonblnk = nonblnk + 1
    else
      exit
    end if
  end do

  ! only parse restraint if it's not blank and if the restraint has not been specified in the old way

  if (nonblnk >= restsize - 5) return

  ! Locate instances where user specifies a center of mass or a plane instead of a single atom
  
  do j=1,8
     if (j > 1) then
        if (comidx(j-1) == 0) exit

        if (index(restraint(comidx(j-1)+3:), 'com') > 0) then
          comidx(j) = comidx(j-1) + 2 + index(restraint(comidx(j-1)+3:), 'com')
        else
          exit
        end if
     else
        comidx(j) = index(restraint, 'com')
     end if
  end do
  

! Identify the kind of restraint, the number of atoms necessary for that restrant, and 
! Locate where in restraint atoms start getting specified.

  if (restraint(nonblnk:nonblnk+7) == 'distance') then
    restatms = 2
    do i=nonblnk+8,restsize
      if (restraint(i:i) == ' ') cycle
      if (restraint(i:i) == '(' .OR. restraint(i:i) == '[' .OR. restraint(i:i) == '{') then
        reststart = i+1
        exit
      else
        exit
      end if
    end do
  else if (restraint(nonblnk:nonblnk+4) == 'angle') then
  
  ! Planes are only specified within angle restraints
  ! Locate where the user specifies planes of atoms
  
    do j=1,2
       if (j > 1) then
          if (plnidx(j-1) == 0) exit
          if ( index(restraint(plnidx(j-1)+5:), 'plane') > 0 ) &
            plnidx(j) = plnidx(j-1) + 4 + index(restraint(plnidx(j-1)+5:), 'plane')
       else
          plnidx(j) = index(restraint, 'plane')
       end if
    end do
  
    if(plnidx(1) > 0 .AND. plnidx(2) == 0) then
      restatms = 5
    else if (plnidx(1) > 0 .AND. plnidx(2) > 0) then
      restatms = 8
    else
      restatms = 3
    end if
    do i=nonblnk+5,restsize
      if (restraint(i:i) == ' ') cycle
      if (restraint(i:i) == '(' .OR. restraint(i:i) == '[' .OR. restraint(i:i) == '{') then
        reststart = i+1
        exit
      else
        exit
      end if
    end do
  else if (restraint(nonblnk:nonblnk+6) == 'torsion') then
    restatms = 4  
    do i=nonblnk+7,restsize
      if (restraint(i:i) == ' ') cycle
      if (restraint(i:i) == '(' .OR. restraint(i:i) == '[' .OR. restraint(i:i) == '{') then
        reststart = i+1
        exit
      else
        exit
      end if
    end do
  else if (restraint(nonblnk:nonblnk+9) == 'coordinate') then
    do j=1,4
       if (j > 1) then
          if (distidx(j-1) == 0) exit
          if ( index(restraint(distidx(j-1)+5:), 'distance') > 0 ) &
            distidx(j) = distidx(j-1) + 7 + index(restraint(distidx(j-1)+8:), 'distance')
       else
          distidx(j) = index(restraint, 'distance')
       end if
    end do
  
    if (distidx(4) > 0) then
      restatms = 8
    else if (distidx(3) > 0) then
      restatms = 6
    else if (distidx(2) > 0) then
      restatms = 4
    else if (distidx(1) > 0) then
      restatms = 2
    else
      write(iout,'(a)') 'Error: no distances specified inside a generalized distance coordinate restraint.'
      write(iout,9999) restraint(1:restsize)
      call mexit(iout,1)
    end if
    
    do i=nonblnk+10,restsize
      if (restraint(i:i) == ' ') cycle
      if (restraint(i:i) == '(' .OR. restraint(i:i) == '[' .OR. restraint(i:i) == '{') then
        reststart = i+1
        exit
      else
        exit
      end if
    end do
  else
    reststart = 0
  end if

  if (reststart == 0) then
    write(iout,9002)
    write(iout,9999) restraint(1:restsize)
    call mexit(iout, 1)
  else ! Begin reading the restraint
    k=reststart
    
    do 
      if (k > restsize) exit
      numlen = 0
      namlen = 0
      if (restraint(k:k) == ' ' .or. restraint(k:k) == ',' .or. restraint(k:k) == '(' &
          .or. restraint(k:k) == ')' .or. restraint(k:k) == '['.or. restraint(k:k) == ']' &
          .or. restraint(k:k) == '{'.or. restraint(k:k) == '}') then
          k = k + 1
          cycle
      end if
      if (iatidx > restatms) then
         write(iout,9000)
         write(iout,9001) restatms
         write(iout,9999) restraint(1:restsize)
         exit
      end if
        
      if ((k < comidx(min(currcom,8)) .or. (comidx(min(currcom,8)) == 0)) & ! Atoms specified not part  of a com or plane designation
      .and. (k < plnidx(min(currpln,2)) .or. plnidx(min(currpln,2)) == 0) &
      .and. (k < distidx(min(currdist,4)) .or. distidx(min(currdist,4)) == 0)) then
         if (iachar(restraint(k:k)) >= 48 .AND. iachar(restraint(k:k)) <= 57) then ! User specifies explicit atom number
            do 
              if (iachar(restraint(k+numlen:k+numlen)) >= 48 .AND. iachar(restraint(k+numlen:k+numlen)) <= 57) then
                numlen = numlen + 1
              else
                exit
              end if
            end do
            read (restraint(k:k+numlen-1), *) iattmp
            iat(iatidx) = iattmp
            iatidx = iatidx + 1
            k = k + numlen
         else if (restraint(k:k) == ':') then  ! User specifies ptraj-style atom mask
         
          call readmask(restraint,k,iat,iatidx,name,ipres,nres,iout)

         else if (restraint(k:k) == '-') then   ! User specifies group of atoms in old way using igr
            if (iachar(restraint(k+1:k+1)) >= 48 .AND. iachar(restraint(k+1:k+1)) <= 57) then
              do
                if (iachar(restraint(k+numlen+1:k+numlen+1)) >= 48 .AND. iachar(restraint(k+numlen:k+numlen)) <= 57) then
                  numlen = numlen + 1
                else
                  exit
                end if
              end do
              
              read(restraint(k:k+numlen), *) iattmp
              iat(iatidx) = iattmp
              
              if (igrarr(1,iatidx) == 0) then
                write(iout,9004)
                write(iout,9999) restraint(1:LEN_TRIM(restraint))
                call mexit(iout, 1)
              end if
              
              iatidx = iatidx + 1
              k = k + numlen + 1
            else
              write(iout,9005)
              write(iout,9999) restraint(1:restsize)
              call mexit(iout, 1)
            end if
         else
            write(iout,9006)
            write(iout,9999) restraint(1:restsize)
            call mexit(iout, 1)
         end if
      else if (k == plnidx(min(currpln,2))) then ! Atoms are within a plane designation
      
         if ((restatms - iatidx + 1) < 4 ) then  ! All planes contain exactly four atoms or COM groupings
           write(iout,9007)
           write(iout,9999) restraint(1:restsize)
           call mexit(iout, 1)
         else
           plnatms = 0
           k = k + 5
           call parenfind(restraint,k,restsize,openparen,closparen,iout)
   
           k=openparen+1
   
         ! The parser is smart enough to figure out the restraint if, in a plane-point restraint,
         ! the user specified the "point" part before the plane, even though the point, whether 
         ! atom or COM, will eventually be iat5.  If this happened, some juggling is necessary.
   
           if(restatms == 5 .and. iat(1) /= 0) then
             iatidx = 1
             iat(5) = iat(1)
             if (currcom > 1) then
               comidx(6) = comidx(1)
               do j=1,4
                 comidx(j) = comidx(j+1)
               end do
               comidx(5) = comidx(6)
               currcom = 1
             end if
           end if
           do
             if(k>closparen-1) exit
             numlen = 0
             namlen = 0
             if (restraint(k:k) == ' ' .or. restraint(k:k) == ',') then
               k = k + 1
               cycle
             end if
             if (k < comidx(min(currcom,8)) .or. (comidx(min(currcom,8)) == 0))  then
               if (iachar(restraint(k:k)) >= 48 .AND. iachar(restraint(k:k)) <= 57) then !explicit atom number
                 do 
                   if (iachar(restraint(k+numlen:k+numlen)) >= 48 .AND. iachar(restraint(k+numlen:k+numlen)) <= 57) then
                     numlen = numlen + 1
                   else
                     exit
                   end if
                 end do
                 read (restraint(k:k+numlen-1), *) iat(iatidx)
                 iatidx = iatidx + 1
                 plnatms = plnatms + 1
                 k = k + numlen
               else if (restraint(k:k) == ':') then  !ptraj-style atom mask
               
                 call readmask(restraint,k,iat,iatidx,name,ipres,nres,iout)

                 plnatms = plnatms + 1
                 
               else if (restraint(k:k) == '-') then   !user specifies group of atoms in old way
                 if (iachar(restraint(k+1:k+1)) >= 48 .AND. iachar(restraint(k+1:k+1)) <= 57) then
                   do
                     if (iachar(restraint(k+numlen+1:k+numlen+1)) >= 48 .AND. iachar(restraint(k+numlen:k+numlen)) <= 57) then
                       numlen = numlen + 1
                     else
                       exit
                     end if
                   end do
                   read(restraint(k:k+numlen), *) iat(iatidx)

                   if (igrarr(1,iatidx) == 0) then
                     write(iout,9004)
                     write(iout,9999) restraint(1:restsize)
                     call mexit(iout, 1)
                   end if
                   
                   iatidx = iatidx + 1
                   plnatms = plnatms + 1
                   k = k + numlen + 1
                 else
                   write(iout,9005)
                   write(iout,9999) restraint(1:restsize)
                   call mexit(iout, 1)
                 end if
               else
                 write(iout,9006)
                 write(iout,9999) restraint(1:restsize)
                 call mexit(iout, 1)
               end if
            else if (i == comidx(min(currcom,8))) then   ! Parse com info from w/in a plane designation
              k = k + 3
              call parenfind(restraint,k,restsize,openparensub,closparensub,iout)
       
              k=openparensub+1
      
              call parscom(restraint,iat,igrarr,k,closparensub,iatidx,restsize,name,ipres,nres,iout)
      
              plnatms = plnatms + 1
              currcom = currcom + 1
              k = closparensub + 1
            else
              write(iout,9006)
              write(iout,9999) restraint(1:restsize)
              call mexit(iout, 1)
            end if
          end do ! Done parsing plane grouping
         end if
         if(plnatms /= 4) then  ! Check that enough atoms were specified in the plane
           write(iout, 9011)
           write(iout,9999) restraint(1:restsize)
           call mexit(iout, 1)
         end if
 
         ! The parser is smart enough to figure out the restraint if, in a plane-point restraint,
         ! the user specified the "point" part before the plane, even though the point, whether 
         ! atom or COM, will eventually be iat5.  If this happened, some juggling is necessary.
 
         if(restatms == 5 .and. iat(5) /= 0) then
           iatidx = 6
         end if
         currpln = currpln + 1
         k = closparen + 1
 
      else if (k == comidx(min(currcom,8))) then
      
        ! Find the parenthetical boundaries of the COM grouping
        k = k + 3
        call parenfind(restraint,k,restsize,openparen,closparen,iout)
        
        k = openparen + 1
        
        ! Parse a user-defined COM grouping
        
        call parscom(restraint,iat,igrarr,k,closparen,iatidx,restsize,name,ipres,nres,iout)
        
        k = closparen + 1
        currcom = currcom + 1

      else if (k == distidx(min(currdist,4))) then
      
        ! Find the parenthetical boundaries of the distance grouping
        k = k + 8
        call parenfind(restraint,k,restsize,openparen,closparen,iout)
        
        k = openparen + 1
        
        ! Parse the distance grouping
        
        call parsgen(restraint,iat,igrarr,2,k,closparen,iatidx,comidx,currcom,restsize,name,ipres,nres,iout)
        
        k = closparen + 1
        
        ! Read the corresponding weight
        
        do
          if (restraint(k:k) == ' ' .or. restraint(k:k) == ',' .or. restraint(k:k) == '(' &
          .or. restraint(k:k) == ')' .or. restraint(k:k) == '['.or. restraint(k:k) == ']' &
          .or. restraint(k:k) == '{'.or. restraint(k:k) == '}') then
            k = k + 1
            cycle
          else 
            exit
          end if
        end do
        
        givemeane = .false.
        givemeadot = .false.
        itsnegative = .false.
        resttmp = ' '
        
        do 
          if (numlen > 32) then 
            exit
          else if ((iachar(restraint(k:k)) >= 48 .AND. iachar(restraint(k:k)) <= 57)) then
            resttmp = resttmp(:numlen) // restraint(k:k)
            numlen = numlen + 1
            k = k + 1
          else if (restraint(k:k) == '-' .and. .not. itsnegative) then
            itsnegative = .true.
            resttmp = resttmp(:numlen) // restraint(k:k)
            numlen = numlen + 1
            k = k + 1
          else if (restraint(k:k) == 'E' .and. .not. givemeane) then
            givemeane = .true.
            resttmp = resttmp(:numlen) // restraint(k:k)
            numlen = numlen + 1
            k = k + 1
            if (restraint(k:k) == '-') then
              k = k + 1
              resttmp = resttmp(:numlen) // restraint(k:k)
              numlen = numlen + 1
            end if
          else if (restraint(k:k) == '.' .and. .not. givemeadot) then
            givemeadot = .true.
            resttmp = resttmp(:numlen) // restraint(k:k)
            numlen = numlen + 1
            k = k + 1
          else if (restraint(k:k) == ' ' .or. restraint(k:k) == ',' .or. restraint(k:k) == '(' &
          .or. restraint(k:k) == ')' .or. restraint(k:k) == '['.or. restraint(k:k) == ']' &
          .or. restraint(k:k) == '{'.or. restraint(k:k) == '}') then
            exit
          else
            write(iout,'(a)') 'Error: Improper weight specified for generalized coordinate restraint.'
            write(iout,9999) restraint(1:restsize)
            call mexit(iout, 1)
          end if
        end do
        
        if (.not. givemeadot .and. .not. givemeane) resttmp = resttmp(:LEN_TRIM(resttmp)) // '.0'
        
        if (resttmp(1:1) == ' ') then
          read(resttmp(2:LEN_TRIM(resttmp)), *) rstwt(currdist)
        else
          read(resttmp(1:LEN_TRIM(resttmp)), *) rstwt(currdist)
        end if
        
        currdist = currdist + 1
        k = k + 1
        
      else
         write(iout,9006)
         write(iout,9999) restraint(1:restsize)
         call mexit(iout, 1)
      end if

    end do
  end if
  ! done parsing restraint
 
  itworked = .true.  ! Quick check to make sure everything is in order.
  do m=1,8
    if(m<=restatms) then
      itworked = itworked .and. (iat(m) /= 0)
      if(iat(m) < 0) itworked = itworked .and. (igrarr(1,m) > 0)
    else
      itworked = itworked .and. (iat(m) == 0)
    end if
  end do
  if(itworked) then
    return
  else
    write(iout,9012)
    write(iout,9999) restraint(1:restsize)
    call mexit(iout, 1)
  end if
  
  ! Error and warning messages
  9000 format('Warning: More characters in restraint beyond final atom necessary')
  9001 format('Using only first ',i1,' atoms.')
  9002 format('Error: No valid restraint type specified.')
  ! 9003 format('Error: Invalid residue number or atom name in restraint.')
  9004 format('Error: Manual atom grouping specified, but none defined.')
  9005 format('Error: Invalid group number in restraint.')
  9006 format('Error: Invalid atom or grouping specified in restraint.')
  9007 format('Error: Too many atoms specified to use plane restraint.')
  ! 9008 format('Error: No closing parenthesis for plane grouping in restraint.')
  ! 9014 format('Error: No closing parenthesis for COM grouping in restraint.')
  ! 9009 format('Error: No parentheses specified to define plane grouping in restraint.')
  ! 9010 format('Error: No atoms specified within COM grouping in restraint.')
  9011 format('Error: Must specify four atoms for a plane grouping.')
  9012 format('Error: Not enough atoms or groupings specified in restraint.')
  ! 9013 format('Warning: More atom groups defined than necessary in restraint.')
  9999 format('restraint = "',(a),'"')
  
  
end subroutine parsrest
  

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calc. the energy from the angle between the normals of two planes 
subroutine plnnrg(x,f,dfr,thet,i1,i2,i3,i4,i5,i6,i7,i8,e,r1,r2,r3,r4,k2,k3,ntb, &
      plave,plave0,nave,nexact,ipower,tauave,ravein, &
      dt,navint,iave,incflg,iravin,iflag)
   
   ! Subroutine PLaNe-plane eNeRGy
   
   ! This subroutine calculates the plane-plane angle restraint energy.
   ! The angles are taken between the planes' normal vectors, which are
   ! in turn defined by the cross products between i1-i2 and i3-i4, and 
   ! between i5-i6 and i7-i8.
   ! The energy function is a flat-bottomed well, with parabolic sides
   ! which become linear at greater distances.
   
   ! If r is the angle between the normals:
   
   !       r <= r1  : E = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
   ! r1 <= r <= r2  : E = k2*(r-r2)**2
   ! r2 <= r <= r3  : E = 0
   ! r3 <= r <= r4  : E = k3*(r-r3)**2
   ! r4 <= r        : E = 2*k3*(r4-r3)*(r-r4) + K3*(r4-r3)**2
   
   
   ! IAVE = 0: then the instantaneous (current) value of r is used
   !           in all calculations.
   ! IAVE = 1: The time-averaged value of r is used, and forces are calculated
   !           according to the above function, i.e.
   !               dE/dx = (dE/dr_ave)(dr_ave/dr)(dr/dx).
   ! IAVE = 2: The time-averaged value of r is used, and forces are calculated as
   !               dE/dx = (dE/dr_ave)(dr/dx).
   !           This expression, if integrated would not yield the energy function
   !           above, so energies reported are pseudo-energies. Also, it is a
   !           non-conservative force, so the system will tend to heat up. But it
   !           avoids potentially large force terms with (1+IPOWER)
   !           exponentiation, which can cause the simulation to become unstable.
   
   ! Other input:
   !    X(I)  : Coordinate array.
   !    F(I)  : Force array;
   !            Modified on output to include forces from this restraint.
   !I1,I2,I3,I4: Atom pointers for the first plane (3*(I-1), where I is absolute
   !            atom number.
   !I5,I6,I7,I8: Atom pointers for the second plane (3*(I-1), where I is
   !            absolute atom number
   !    NTB   : Periodic boundary conditions flag.
   !    IFLAG : =0, then calculate energy and derivatives
   !            =1, then only calculate current value of angle
   !            =2, then calculate energy and derivatives, but do not
   !                accumulate forces into F() array. They are returned
   !                in the DFR array.
   !            =3, the calculate current value of angle and energy, but
   !                do not accumulate derivatives
   !   INCFLG : Determines whether local saved pointers are updated in AVEINT
   
   ! (the following only used when time-averaging):
   ! PLAVE(I) : Contains time-averaged values of angles, using TAUAVE damping.
   ! PLAVE0(I): Contains time-averaged values of angles, using no damping.
   
   !           Both PLAVE and PLAVE0 are passed so that the first element of
   !           the array corresponds to the average value of the angle currently
   !           of interest, and the second element of the array contains
   !           the previous value of the angle of interest.
   
   !     NAVE : > 0 for time-averaging.
   !   NEXACT : Not currently used.
   !   IPOWER : The power used in averaging the angle data.
   !   TAUAVE : The exponential decay factor used in averaging the angle data
   !   IRAVIN : =0 do not modify the values of the angles by RAVEIN
   !            =1 modify theta by RAVEIN, as shown below.
   !   RAVEIN : At time 0, time-average integral is undefined.
   !            If -1000<RAVEIN<1000, use r_initial+ravein on this step.
   !            If RAVEIN<=-1000, use r_target+(ravein+1000) on this step.
   !            If RAVEIN>=1000, use r_target-(ravein-1000) on this step.
   !       DT : The integration time-step (ps)
   !   NAVINT : Time-averaged restraints are only updated every NAVINT
   !            steps. Typically, NAVINT=1.
   
   ! Output:
   !    E     : The energy contribution from this restraint.
   !    DFR   : If IFLAG=2, DFR(1->3) contains the d_E/d_distance forces
   !            for x,y, and z for position 1; those for position 2 are in
   !            DFR(4->6), and so on for each of the 8 atoms here
   !  THET    : The calculated valence angle
   
   ! Author: Matthew G. Seetin
   ! Date: 8/2007
   ! This routine is based heavily on the angnrg subroutine of David A. Pearlman
   ! from 7/89.
   use constants, only : DEG_TO_RAD, half
   implicit none
   integer:: i1, i2, i3, i4, i5, i6, i7, i8, iave, iflag, iflg, &
        incflg, ipower, irav, iravin, m, nave, navint, nexact, &
        ntb
   _REAL_ :: cst, denom, df, dfr, dif, dif1, dravdr, dt, e, f, &
        plave, plave0, r1, r2, r3, r4, ravein, rinc, rnow, small, small2, &
        snt, st, tauave, thet, x
   _REAL_ k2,k3,xab,xcd,xef,xgh,n1,n2,rdenom1,rdenom2
   parameter (small=1.0d-14)
   parameter (small2=1.0d-5)
   dimension x(*),f(*),dfr(*),plave(2),plave0(2)
   dimension xab(3),xcd(3),xef(3),xgh(3),n1(3),n2(3)
   rinc = 1000.0d0
   
   iflg = 1
   if (iflag == 1) iflg = 2
   
   ! Calculate the angle:
   
   do m=1,3
      xab(m) = x(i1+m) - x(i2+m)
      xcd(m) = x(i3+m) - x(i4+m)
      xef(m) = x(i5+m) - x(i6+m)
      xgh(m) = x(i7+m) - x(i8+m)
   end do
   
   ! Calculate normal vectors n1 and n2 from the cross products of xab and xcd,
   ! and xef and xgh

   n1(1) = xab(2)*xcd(3) - xab(3)*xcd(2)
   n1(2) = -xab(1)*xcd(3) + xab(3)*xcd(1)
   n1(3) = xab(1)*xcd(2) - xab(2)*xcd(1)

   n2(1) = xef(2)*xgh(3) - xef(3)*xgh(2)
   n2(2) = -xef(1)*xgh(3) + xef(3)*xgh(1)   
   n2(3) = xef(1)*xgh(2) - xef(2)*xgh(1)

   ! Calculate angle from n1*n2 = cos(theta)

   rdenom1 = sqrt(n1(1)**2 + n1(2)**2 + n1(3)**2)
   rdenom2 = sqrt(n2(1)**2 + n2(2)**2 + n2(3)**2)
   
   ! If the vectors defining either plane are parallel, we cannot calculate a value for 
   ! this restraint, so set energy to zero and return.  This does not conserve energy, 
   ! but otherwise this can result in a calculation blowing up for a relatively unclear reason.
   ! Users are encouraged to define their planes carefully.  A similar choice is made for 
   ! torsions if three consecutive atoms in the torsion are approx. colinear.  
   
   if (rdenom1 < 0.1d0 .or. rdenom2 < 0.1d0) then
      e = 0.0d0
      return
   end if
   
   ! Calculate angle from n1*n2 = cos(theta)
   
   
   cst = (n1(1)*n2(1) + n1(2)*n2(2) + n1(3)*n2(3))/rdenom1/rdenom2
   if (cst > 1.0d0) cst = 1.0d0
   if (cst < -1.0d0) cst = -1.0d0
   thet = acos(cst)
   
   ! Get averaged value, if requested
   
   dravdr = 1.0d0
   if (iave > 0) then
      rnow = thet
      ! Set the value on the first step (where integral not defined):
      if (iravin == 1) then
         irav = int(ravein)
         if (-1000 < irav .and. irav < 1000) then
            rnow = rnow + ravein*DEG_TO_RAD
         else
            rinc = -sign(rinc,ravein)
            rinc = (ravein+rinc) * DEG_TO_RAD
            if (r2 > small2) rnow =((r2+r3)*half) + rinc
            if (r2 <= small2) rnow = r3 + rinc
         end if
      end if
      call aveint(plave,plave0,rnow,tauave,dt,navint,ipower,6, &
            incflg,iflg,thet,denom)
      
      ! DRAVDR is the factor d(r_ave)/d(r_current)
      
      if (iave == 1 .and. iflag /= 1) &
            dravdr = ((thet/rnow)**(1+ipower))*dt/denom
   end if
   
   if (iflag == 1) return
   
   ! Calculate energy (E) and the derivative with respect to the valence
   ! angle (DF):
   
   if (thet < r1) then
      dif1 = r1-r2
      df = 2.0d0 * k2 * dif1
      e = df * (thet-r1) + k2*dif1*dif1
   else if (thet < r2) then
      dif = thet - r2
      df = 2.0d0 * k2 * dif
      e = k2*dif*dif
   else if (thet <= r3) then
      e = 0.0d0
      return
   else if (thet < r4) then
      dif = thet - r3
      df = 2.0d0 * k3 * dif
      e = k3*dif*dif
   else
      dif1 = r4-r3
      df = 2.0d0 * k3 * dif1
      e = df * (thet-r4) + k3*dif1*dif1
   end if
   if (iflag == 3) return
   
   ! Calculate the derivaties with respect to the coordinates, and add them
   ! into the force arrays. DRAVDR contains d(thet_ave)/d(thet(t)) if IAVE=1,
   ! 1.0 otherwise.
   
   df = dravdr * df
   snt = sin(thet)
   if (abs(snt) < small) snt = small
   st = -df/snt
   
   dfr(1) = st * ( (xcd(3)*xef(1)*xgh(3) - xcd(3)*xef(3)*xgh(1) &
       + xcd(2)*xef(1)*xgh(2) - xcd(2)*xef(2)*xgh(1))/rdenom1/rdenom2 &
       - cst/rdenom2/(rdenom1**3)*((xab(1)*xcd(3)-xab(3)*xcd(1))*xcd(3) &
       + (xab(1)*xcd(2)-xab(2)*xcd(1))*xcd(2)))
   dfr(2) = st * ( (xcd(3)*xef(2)*xgh(3) - xcd(3)*xef(3)*xgh(2) &
       - xcd(1)*xef(1)*xgh(2) + xcd(1)*xef(2)*xgh(1))/rdenom1/rdenom2 &
       - cst/rdenom2/(rdenom1**3)*((xab(2)*xcd(3)-xab(3)*xcd(2))*xcd(3) &
       - (xab(1)*xcd(2)-xab(2)*xcd(1))*xcd(1)))
   dfr(3) = st * ( (-xcd(2)*xef(2)*xgh(3) + xcd(2)*xef(3)*xgh(2) &
       - xcd(1)*xef(1)*xgh(3) + xcd(1)*xef(3)*xgh(1))/rdenom1/rdenom2 &
       - cst/rdenom2/(rdenom1**3)*((xab(1)*xcd(3)-xab(3)*xcd(1))*xcd(1) &
       - (xab(3)*xcd(2)-xab(2)*xcd(3))*xcd(2)))
   dfr(4) = -dfr(1)
   dfr(5) = -dfr(2)
   dfr(6) = -dfr(3)
   dfr(7) = st * ( (-xab(3)*xef(1)*xgh(3) + xab(3)*xef(3)*xgh(1) &
       - xab(2)*xef(1)*xgh(2) + xab(2)*xef(2)*xgh(1))/rdenom1/rdenom2 &
       - cst/rdenom2/(rdenom1**3)*(-(xab(1)*xcd(3)-xab(3)*xcd(1))*xab(3) &
       - (xab(1)*xcd(2)-xab(2)*xcd(1))*xab(2)))
   dfr(8) = st * ( (-xab(3)*xef(2)*xgh(3) + xab(3)*xef(3)*xgh(2) &
       + xab(1)*xef(1)*xgh(2) - xab(1)*xef(2)*xgh(1))/rdenom1/rdenom2 &
       - cst/rdenom2/(rdenom1**3)*((xab(1)*xcd(2)-xab(2)*xcd(1))*xcd(1) &
       - (xab(3)*xcd(2)-xab(2)*xcd(3))*xcd(3)))
   dfr(9) = st * ( (xab(2)*xef(2)*xgh(3) - xab(2)*xef(3)*xgh(2) &
       + xab(1)*xef(1)*xgh(3) - xab(1)*xef(3)*xgh(1))/rdenom1/rdenom2 &
       - cst/rdenom2/(rdenom1**3)*((xab(2)*xcd(3)-xab(3)*xcd(2))*xab(2) &
       + (xab(1)*xcd(3)-xab(3)*xcd(1))*xab(1)))
   dfr(10) = -dfr(7)
   dfr(11) = -dfr(8)
   dfr(12) = -dfr(9)
   dfr(13) = st * ( (xcd(3)*xab(1)*xgh(3) - xcd(1)*xab(3)*xgh(3) &
       + xcd(2)*xab(1)*xgh(2) - xcd(1)*xab(2)*xgh(2))/rdenom1/rdenom2 &
       - cst/rdenom1/(rdenom2**3)*((xef(1)*xgh(3)-xef(3)*xgh(1))*xgh(3) &
       + (xef(1)*xgh(2)-xef(2)*xgh(1))*xgh(2)))
   dfr(14) = st * ( (xcd(3)*xab(2)*xgh(3) - xcd(2)*xab(3)*xgh(3) &
       - xcd(2)*xab(1)*xgh(1) + xcd(1)*xab(2)*xgh(1))/rdenom1/rdenom2 &
       - cst/rdenom1/(rdenom2**3)*((xef(2)*xgh(3)-xef(3)*xgh(2))*xgh(3) &
       - (xef(1)*xgh(2)-xef(2)*xgh(1))*xgh(1)))
   dfr(15) = st * ( (-xcd(3)*xab(2)*xgh(2) + xcd(2)*xab(3)*xgh(2) &
       - xcd(3)*xab(1)*xgh(1) + xcd(1)*xab(3)*xgh(1))/rdenom1/rdenom2 &
       - cst/rdenom1/(rdenom2**3)*((xef(1)*xgh(3)-xef(3)*xgh(1))*xgh(1) &
       - (xef(3)*xgh(2)-xef(2)*xgh(3))*xgh(2)))
   dfr(16) = -dfr(13)
   dfr(17) = -dfr(14)
   dfr(18) = -dfr(15)
   dfr(19) = st * ( (-xcd(3)*xef(3)*xab(1) + xcd(1)*xef(3)*xab(3) &
       - xcd(2)*xef(2)*xab(1) + xcd(1)*xef(2)*xab(2))/rdenom1/rdenom2 &
       - cst/rdenom1/(rdenom2**3)*(-(xef(1)*xgh(3)-xef(3)*xgh(1))*xef(3) &
       - (xef(1)*xgh(2)-xef(2)*xgh(1))*xef(2)))
   dfr(20) = st * ( (-xcd(3)*xef(3)*xab(2) + xcd(2)*xef(3)*xab(3) &
       + xcd(2)*xef(2)*xab(1) - xcd(1)*xef(2)*xab(2))/rdenom1/rdenom2 &
       - cst/rdenom1/(rdenom2**3)*((xef(1)*xgh(2)-xef(2)*xgh(1))*xgh(1) &
       - (xef(3)*xgh(2)-xef(2)*xgh(3))*xgh(3)))
   dfr(21) = st * ( (xcd(3)*xef(2)*xab(2) - xcd(2)*xef(2)*xab(3) &
       + xcd(3)*xef(1)*xab(1) - xcd(1)*xef(1)*xab(3))/rdenom1/rdenom2 &
       - cst/rdenom1/(rdenom2**3)*((xef(2)*xgh(3)-xef(3)*xgh(2))*xef(2) &
       + (xef(1)*xgh(3)-xef(3)*xgh(1))*xef(1)))
   dfr(22) = -dfr(19)
   dfr(23) = -dfr(20)
   dfr(24) = -dfr(21)   
   
   if(iflag == 2) return
   
   f(i1+1) = f(i1+1) - dfr(1)
   f(i1+2) = f(i1+2) - dfr(2)
   f(i1+3) = f(i1+3) - dfr(3)
   f(i2+1) = f(i2+1) - dfr(4)
   f(i2+2) = f(i2+2) - dfr(5)
   f(i2+3) = f(i2+3) - dfr(6)
   f(i3+1) = f(i3+1) - dfr(7)
   f(i3+2) = f(i3+2) - dfr(8)
   f(i3+3) = f(i3+3) - dfr(9)
   f(i4+1) = f(i4+1) - dfr(10)
   f(i4+2) = f(i4+2) - dfr(11)
   f(i4+3) = f(i4+3) - dfr(12)
   f(i5+1) = f(i5+1) - dfr(13)
   f(i5+2) = f(i5+2) - dfr(14)
   f(i5+3) = f(i5+3) - dfr(15)
   f(i6+1) = f(i6+1) - dfr(16)
   f(i6+2) = f(i6+2) - dfr(17)
   f(i6+3) = f(i6+3) - dfr(18)
   f(i7+1) = f(i7+1) - dfr(19)
   f(i7+2) = f(i7+2) - dfr(20)
   f(i7+3) = f(i7+3) - dfr(21)
   f(i8+1) = f(i8+1) - dfr(22)
   f(i8+2) = f(i8+2) - dfr(23)
   f(i8+3) = f(i8+3) - dfr(24)

   
   return
end subroutine plnnrg 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calc. the energy from angle between the normal of a plane and an atom
subroutine plptnrg(x,f,dfr,rimassplpt,thet,i1,i2,i3,i4,i5,e,r1,r2,r3,r4,k2,k3,ntb, &
      ptave,ptave0,nave,nexact,ipower,tauave,ravein, &
      dt,navint,iave,incflg,iravin,iflag)
   
   ! Subroutine PLane-PoinT eNeRGy
   
   ! This subroutine calculates the plane-point angle restraint energy.
   ! The restraint is an angle restraint between a plane and a single atom.
   ! The angle is taken between the plane's normal vector, which are
   ! in turn defined by the cross products between i1-i2 and i3-i4, and 
   ! the vector from the center of mass of the plane atoms to the point,
   ! atom i5.  
   ! The energy function is a flat-bottomed well, with parabolic sides
   ! which become linear at greater distances.
   
   ! If r is the valence angle:
   
   !       r <= r1  : E = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
   ! r1 <= r <= r2  : E = k2*(r-r2)**2
   ! r2 <= r <= r3  : E = 0
   ! r3 <= r <= r4  : E = k3*(r-r3)**2
   ! r4 <= r        : E = 2*k3*(r4-r3)*(r-r4) + K3*(r4-r3)**2
   
   
   ! IAVE = 0: then the instantaneous (current) value of r is used
   !           in all calculations.
   ! IAVE = 1: The time-averaged value of r is used, and forces are calculated
   !           according to the above function, i.e.
   !               dE/dx = (dE/dr_ave)(dr_ave/dr)(dr/dx).
   ! IAVE = 2: The time-averaged value of r is used, and forces are calculated as
   !               dE/dx = (dE/dr_ave)(dr/dx).
   !           This expression, if integrated would not yield the energy function
   !           above, so energies reported are pseudo-energies. Also, it is a
   !           non-conservative force, so the system will tend to heat up. But it
   !           avoids potentially large force terms with (1+IPOWER)
   !           exponentiation, which can cause the simulation to become unstable.
   
   ! Other input:
   !    X(I)  : Coordinate array.
   !    F(I)  : Force array;
   !            Modified on output to include forces from this restraint.
   ! RIMASS(I): Array of inverse masses.
   !I1,I2,I3,I4: Atom pointers for the first plane (3*(I-1), where I is absolute
   !            atom number.
   !    I5    : Atom number for the "point" atom.
   !    NTB   : Periodic boundary conditions flag.
   !    IFLAG : =0, then calculate energy and derivatives
   !            =1, then only calculate current value of angle
   !            =2, then calculate energy and derivatives, but do not
   !                accumulate forces into F() array. They are returned
   !                in the DFR array.
   !            =3, the calculate current value of angle and energy, but
   !                do not accumulate derivatives
   !   INCFLG : Determines whether local saved pointers are updated in AVEINT
   
   ! (the following only used when time-averaging):
   ! PTAVE(I) : Contains time-averaged values of angles, using TAUAVE damping.
   ! PTAVE0(I): Contains time-averaged values of angles, using no damping.
   
   !           Both AAVE and AAVE0 are passed so that the first element of
   !           the array corresponds to the average value of the angle currently
   !           of interest, and the second element of the array contains
   !           the previous value of the angle of interest.
   
   !     NAVE : > 0 for time-averaging.
   !   NEXACT : Not currently used.
   !   IPOWER : The power used in averaging the angle data.
   !   TAUAVE : The exponential decay factor used in averaging the angle data
   !   IRAVIN : =0 do not modify the values of the angles by RAVEIN
   !            =1 modify theta by RAVEIN, as shown below.
   !   RAVEIN : At time 0, time-average integral is undefined.
   !            If -1000<RAVEIN<1000, use r_initial+ravein on this step.
   !            If RAVEIN<=-1000, use r_target+(ravein+1000) on this step.
   !            If RAVEIN>=1000, use r_target-(ravein-1000) on this step.
   !       DT : The integration time-step (ps)
   !   NAVINT : Time-averaged restraints are only updated every NAVINT
   !            steps. Typically, NAVINT=1.
   
   ! Output:
   !    E     : The energy contribution from this restraint.
   !    DFR   : If IFLAG=2, DFR(1->3) contains the d_E/d_distance forces
   !            for x,y, and z for position 1; those for position 2 are in
   !            DFR(4->6), and so on for each of the 5 atoms here
   !  THET    : The calculated valence angle
   
   ! Author: Matthew G. Seetin
   ! Date: 8/2007
   ! This routine is based heavily on the angnrg subroutine of David A. Pearlman
   ! from 7/89.
   use constants, only : DEG_TO_RAD, half
   implicit none
   integer:: i1, i2, i3, i4, i5, iave, iflag, iflg, incflg, ipower, &
        irav, iravin, l, m, nave, navint, nexact, ntb
   _REAL_ :: cst, denom, df, dfr, dif, dif1, dravdr, dt, e, f, &
        ptave, ptave0, r1, r2, r3, r4, ravein, rimassplpt, rinc, rnow, &
        small, small2, snt, st, tauave, thet, x
   _REAL_ k2,k3,xab,xcd,n1,rab2,rcd2,rab,rcd,rdenom1,rdenom2,rdenom2sq,bigm,plncom,vector2
   parameter (small=1.0d-14)
   parameter (small2=1.0d-5)
   dimension x(*),f(*),dfr(*),rimassplpt(*),ptave(2),ptave0(2)
   dimension xab(3),xcd(3),n1(3),vector2(3),plncom(3)
   rinc = 1000.0d0

   ! Calculate COM
 
   bigm = 1.0/rimassplpt(i1/3+1) + 1.0/rimassplpt(i2/3+1) + 1.0/rimassplpt(i3/3+1) + 1.0/rimassplpt(i4/3+1)
   do l=1,3
       plncom(l) = (x(i1+l)/rimassplpt(i1/3+1) + x(i2+l)/rimassplpt(i2/3+1) + &
       x(i3+l)/rimassplpt(i3/3+1) + x(i4+l)/rimassplpt(i4/3+1))/bigm
   end do
   
   iflg = 1
   if (iflag == 1) iflg = 2
   
   ! Calculate the relevant vectors:
   
   rab2 = 0.0d0
   rcd2 = 0.0d0
   rdenom2sq = 0.0d0
   
   do m=1,3
      xab(m) = x(i1+m) - x(i2+m)
      xcd(m) = x(i3+m) - x(i4+m)
      rab2 = rab2 + xab(m)**2
      rcd2 = rcd2 + xcd(m)**2

      vector2(m) = x(i5+m) - plncom(m)
      rdenom2sq = rdenom2sq + vector2(m)**2
      
   end do
   rab = sqrt(rab2)
   rcd = sqrt(rcd2)
   rdenom2 = sqrt(rdenom2sq)
   
   ! Calculate normal vector n1 from the cross product of xab and xcd,

   n1(1) = xab(2)*xcd(3) - xab(3)*xcd(2)
   n1(2) = -xab(1)*xcd(3) + xab(3)*xcd(1)
   n1(3) = xab(1)*xcd(2) - xab(2)*xcd(1)


   ! Calculate angle from n1*(x(i5)-plncom) = cos(theta)

   rdenom1 = sqrt(n1(1)**2 + n1(2)**2 + n1(3)**2)

   
   ! If the vectors defining the plane are parallel, we cannot calculate a value for 
   ! this restraint, so set energy to zero and return.  This does not conserve energy, 
   ! but otherwise this can result in a calculation blowing up for a relatively unclear reason.
   ! Users are encouraged to define their planes carefully.  A similar choice is made for 
   ! torsions if three consecutive atoms in the torsion are approx. colinear.  
   
   if (rdenom1 < 0.001d0) then
      e = 0.0d0
      return
   end if
   
   ! Calculate angle from n1*(x(i5)-plncom) = cos(theta)
   
   cst = (n1(1)*vector2(1) + n1(2)*vector2(2) + n1(3)*vector2(3))/rdenom1/rdenom2
   if (cst > 1.0d0) cst = 1.0d0
   if (cst < -1.0d0) cst = -1.0d0
   thet = acos(cst)
   
   ! Get averaged value, if requested
   
   dravdr = 1.0d0
   if (iave > 0) then
      rnow = thet
      ! Set the value on the first step (where integral not defined):
      if (iravin == 1) then
         irav = int(ravein)
         if (-1000 < irav .and. irav < 1000) then
            rnow = rnow + ravein*DEG_TO_RAD
         else
            rinc = -sign(rinc,ravein)
            rinc = (ravein+rinc) * DEG_TO_RAD
            if (r2 > small2) rnow =((r2+r3)*half) + rinc
            if (r2 <= small2) rnow = r3 + rinc
         end if
      end if
      call aveint(ptave,ptave0,rnow,tauave,dt,navint,ipower,5, &
            incflg,iflg,thet,denom)
      
      ! DRAVDR is the factor d(r_ave)/d(r_current)
      
      if (iave == 1 .and. iflag /= 1) &
            dravdr = ((thet/rnow)**(1+ipower))*dt/denom
   end if
   
   if (iflag == 1) return
   
   ! Calculate energy (E) and the derivative with respect to the valence
   ! angle (DF):
   
   if (thet < r1) then
      dif1 = r1-r2
      df = 2.0d0 * k2 * dif1
      e = df * (thet-r1) + k2*dif1*dif1
   else if (thet < r2) then
      dif = thet - r2
      df = 2.0d0 * k2 * dif
      e = k2*dif*dif
   else if (thet <= r3) then
      e = 0.0d0
      return
   else if (thet < r4) then
      dif = thet - r3
      df = 2.0d0 * k3 * dif
      e = k3*dif*dif
   else
      dif1 = r4-r3
      df = 2.0d0 * k3 * dif1
      e = df * (thet-r4) + k3*dif1*dif1
   end if
   if (iflag == 3) return
   
   ! Calculate the derivaties with respect to the coordinates, and add them
   ! into the force arrays. DRAVDR contains d(thet_ave)/d(thet(t)) if IAVE=1,
   ! 1.0 otherwise.
   
   df = dravdr * df
   snt = sin(thet)
   if (abs(snt) < small) snt = small
   st = -df/snt
   
   dfr(1) = st * ( ( -(xab(2)*xcd(3) - xab(3)*xcd(2))/rimassplpt(i1/3+1)/bigm &
       - xcd(3)*vector2(2) + xcd(2)*vector2(3))/rdenom1/rdenom2 &
       - ((xab(1)*xcd(3)-xab(3)*xcd(1))*xcd(3) + (xab(1)*xcd(2) - xab(2)*xcd(1))*xcd(2))*cst/(rdenom1**2) &
       + (vector2(1)/rimassplpt(i1/3+1)/bigm)*cst/rdenom2sq)
   dfr(2) = st * ( ( (xab(1)*xcd(3) - xab(3)*xcd(1))/rimassplpt(i1/3+1)/bigm &
       + xcd(3)*vector2(1) - xcd(1)*vector2(3))/rdenom1/rdenom2 &
       - ((xab(2)*xcd(3)-xab(3)*xcd(2))*xcd(3) - (xab(1)*xcd(2) - xab(2)*xcd(1))*xcd(1))*cst/(rdenom1**2) &
       + (vector2(2)/rimassplpt(i1/3+1)/bigm)*cst/rdenom2sq)
   dfr(3) = st * ( ( -(xab(1)*xcd(2) - xab(2)*xcd(1))/rimassplpt(i1/3+1)/bigm &
       - xcd(2)*vector2(1) + xcd(1)*vector2(2))/rdenom1/rdenom2 &
       + ((xab(1)*xcd(3)-xab(3)*xcd(1))*xcd(1) - (xab(2)*xcd(3) - xab(3)*xcd(2))*xcd(2))*cst/(rdenom1**2) &
       + (vector2(3)/rimassplpt(i1/3+1)/bigm)*cst/rdenom2sq)
   dfr(4) = st * ( ( -(xab(2)*xcd(3) - xab(3)*xcd(2))/rimassplpt(i2/3+1)/bigm &
       + xcd(3)*vector2(2) - xcd(2)*vector2(3))/rdenom1/rdenom2 &
       + ((xab(1)*xcd(3)-xab(3)*xcd(1))*xcd(3) + (xab(1)*xcd(2) - xab(2)*xcd(1))*xcd(2))*cst/(rdenom1**2) &
       + (vector2(1)/rimassplpt(i2/3+1)/bigm)*cst/rdenom2sq)
   dfr(5) = st * ( ( (xab(1)*xcd(3) - xab(3)*xcd(1))/rimassplpt(i2/3+1)/bigm &
       - xcd(3)*vector2(1) + xcd(1)*vector2(3))/rdenom1/rdenom2 &
       + ((xab(2)*xcd(3)-xab(3)*xcd(2))*xcd(3) - (xab(1)*xcd(2) - xab(2)*xcd(1))*xcd(1))*cst/(rdenom1**2) &
       + (vector2(2)/rimassplpt(i2/3+1)/bigm)*cst/rdenom2sq)
   dfr(6) = st * ( ( -(xab(1)*xcd(2) - xab(2)*xcd(1))/rimassplpt(i2/3+1)/bigm &
       + xcd(2)*vector2(1) - xcd(1)*vector2(2))/rdenom1/rdenom2 &
       - ((xab(1)*xcd(3)-xab(3)*xcd(1))*xcd(1) - (xab(2)*xcd(3) - xab(3)*xcd(2))*xcd(2))*cst/(rdenom1**2) &
       + (vector2(3)/rimassplpt(i2/3+1)/bigm)*cst/rdenom2sq)
   dfr(7) = st * ( ( -(xab(2)*xcd(3) - xab(3)*xcd(2))/rimassplpt(i3/3+1)/bigm &
       + xab(3)*vector2(2) - xab(2)*vector2(3))/rdenom1/rdenom2 &
       + ((xab(1)*xcd(3)-xab(3)*xcd(1))*xab(3) + (xab(1)*xcd(2) - xab(2)*xcd(1))*xab(2))*cst/(rdenom1**2) &
       + (vector2(1)/rimassplpt(i3/3+1)/bigm)*cst/rdenom2sq)
   dfr(8) = st * ( ( (xab(1)*xcd(3) - xab(3)*xcd(1))/rimassplpt(i3/3+1)/bigm &
       - xab(3)*vector2(1) + xab(1)*vector2(3))/rdenom1/rdenom2 &
       + ((xab(2)*xcd(3)-xab(3)*xcd(2))*xab(3) - (xab(1)*xcd(2) - xab(2)*xcd(1))*xab(1))*cst/(rdenom1**2) &
       + (vector2(2)/rimassplpt(i3/3+1)/bigm)*cst/rdenom2sq)
   dfr(9) = st * ( ( -(xab(1)*xcd(2) - xab(2)*xcd(1))/rimassplpt(i3/3+1)/bigm &
       + xab(2)*vector2(1) - xab(1)*vector2(2))/rdenom1/rdenom2 &
       - ((xab(1)*xcd(3)-xab(3)*xcd(1))*xab(1) + (xab(2)*xcd(3) - xab(3)*xcd(2))*xab(2))*cst/(rdenom1**2) &
       + (vector2(3)/rimassplpt(i3/3+1)/bigm)*cst/rdenom2sq)
   dfr(10) = st * ( ( -(xab(2)*xcd(3) - xab(3)*xcd(2))/rimassplpt(i4/3+1)/bigm &
       - xab(3)*vector2(2) + xab(2)*vector2(3))/rdenom1/rdenom2 &
       - ((xab(1)*xcd(3)-xab(3)*xcd(1))*xab(3) + (xab(1)*xcd(2) - xab(2)*xcd(1))*xab(2))*cst/(rdenom1**2) &
       + (vector2(1)/rimassplpt(i4/3+1)/bigm)*cst/rdenom2sq)
   dfr(11) = st * ( ( (xab(1)*xcd(3) - xab(3)*xcd(1))/rimassplpt(i4/3+1)/bigm &
       + xab(3)*vector2(1) - xab(1)*vector2(3))/rdenom1/rdenom2 &
       - ((xab(2)*xcd(3)-xab(3)*xcd(2))*xab(3) - (xab(1)*xcd(2) - xab(2)*xcd(1))*xab(1))*cst/(rdenom1**2) &
       + (vector2(2)/rimassplpt(i4/3+1)/bigm)*cst/rdenom2sq)
   dfr(12) = st * ( ( -(xab(1)*xcd(2) - xab(2)*xcd(1))/rimassplpt(i4/3+1)/bigm &
       - xab(2)*vector2(1) + xab(1)*vector2(2))/rdenom1/rdenom2 &
       + ((xab(1)*xcd(3)-xab(3)*xcd(1))*xab(1) + (xab(2)*xcd(3) - xab(3)*xcd(2))*xab(2))*cst/(rdenom1**2) &
       + (vector2(3)/rimassplpt(i4/3+1)/bigm)*cst/rdenom2sq)
   dfr(13) = st * ((xab(2)*xcd(3)-xab(3)*xcd(2))/rdenom1/rdenom2 + cst/rdenom2sq )
   dfr(14) = st * (-(xab(1)*xcd(3)-xab(3)*xcd(1))/rdenom1/rdenom2 + cst/rdenom2sq )
   dfr(15) = st * ((xab(1)*xcd(2)-xab(2)*xcd(1))/rdenom1/rdenom2 + cst/rdenom2sq )
   
   if(iflag == 2) return
   
   f(i1+1) = f(i1+1) - dfr(1)
   f(i1+2) = f(i1+2) - dfr(2)
   f(i1+3) = f(i1+3) - dfr(3)
   f(i2+1) = f(i2+1) - dfr(4)
   f(i2+2) = f(i2+2) - dfr(5)
   f(i2+3) = f(i2+3) - dfr(6)
   f(i3+1) = f(i3+1) - dfr(7)
   f(i3+2) = f(i3+2) - dfr(8)
   f(i3+3) = f(i3+3) - dfr(9)
   f(i4+1) = f(i4+1) - dfr(10)
   f(i4+2) = f(i4+2) - dfr(11)
   f(i4+3) = f(i4+3) - dfr(12)
   f(i5+1) = f(i5+1) - dfr(13)
   f(i5+2) = f(i5+2) - dfr(14)
   f(i5+3) = f(i5+3) - dfr(15)


   
   
   return
end subroutine plptnrg 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine r6ave here]
subroutine r6ave(x,nmrat,nmrcom,rave,indx)
   
   
   ! Subroutine R**-6 AVErage.
   
   ! This routine calculate the <r**-6>-6 average distance between two
   ! sets of atoms.
   
   ! The resulting averaged distance is placed in RAVE.
   ! A number of other quantities, required in calculating the derivatives,
   ! are calculated and stored in common block R6AVEC. This common block
   ! only appears in this routine and in R6DRV (where the derivatives are
   ! determined).
   
   ! Note: This routine assumes that all atoms of the group are in the
   !       same periodic "box" (if applicable).
   
   !       The maximum allowable number of pairs of atomic interactions used
   !       in calculating the r**-6 average is currently hardwired in this routine
   !       (MAXPR). If you change it here, you must also change it in R6DRV.
   
   
   implicit none
   integer:: i, i3, iat0, iat1, iat2, iat3, iat4, indx, ip, ipr, j, &
        j3, jat0, jp, maxpr, nmrat, nmrcom, nsum
   _REAL_ :: fsum, one6, rave, rij, rij2, rm2, rm6bar, rm8, x, xij
   parameter (maxpr = 5000)
   parameter (one6 = 1.0d0/6.0d0)
   
   ! r6avec common block. Shared between this routine and R6DRV.
   
   common/r6avec/xij(3,maxpr),iat0(maxpr), &
         jat0(maxpr),rm8(maxpr),fsum,rm6bar,rij,nsum
#  include "box.h"
#ifdef LES
#  include "les.h"
#endif
   
   dimension x(*),nmrat(16,*),nmrcom(2,*)
   
   iat1 = nmrat(1,indx)
   iat2 = nmrat(2,indx)
   iat3 = nmrat(3,indx)
   iat4 = nmrat(4,indx)
   
   ! Case where IAT1 > 0 (single atom) and IAT2 < 0 (r**-6 average pos.).
   
   ipr = 0
   rm6bar = 0.0d0
   if (iat1 >= 0) then
      i3 = iat1
      do jp=-iat2,-nmrat(10,indx)
      ! Modified 9/2007 by Matthew Seetin as a consequence of enabling restraints on
      ! group of atoms for angle, torsion, and my new planar restraints
         do j=nmrcom(1,jp),nmrcom(2,jp)
            j3 = 3*(j-1)
            ipr = ipr + 1
            if (ipr > maxpr) goto 500
            xij(1,ipr) = x(i3+1) - x(j3+1)
            xij(2,ipr) = x(i3+2) - x(j3+2)
            xij(3,ipr) = x(i3+3) - x(j3+3)
            iat0(ipr) = i3
            jat0(ipr) = j3
            rij2 = xij(1,ipr)**2 + xij(2,ipr)**2 + xij(3,ipr)**2
            rm2 = 1.0d0/rij2
            rm6bar = rm6bar + rm2**3
            rm8(ipr) = rm2**4
         end do
      end do
      nsum = ipr
      fsum = nsum
      rm6bar = rm6bar/fsum
      rij = rm6bar**(-one6)
      
      ! Case where IAT1 < 0 (r**-6 ave. pos.)  and IAT2 > 0
      
   else if (iat2 >= 0) then
      
      j3 = iat2
      do ip=-iat1,-nmrat(9,indx)
      ! Modified 9/2007 by Matthew Seetin as a consequence of enabling restraints on
      ! group of atoms for angle, torsion, and my new planar restraints
         do i=nmrcom(1,ip),nmrcom(2,ip)
            i3 = 3*(i-1)
            ipr = ipr + 1
            if (ipr > maxpr) goto 500
            xij(1,ipr) = x(i3+1) - x(j3+1)
            xij(2,ipr) = x(i3+2) - x(j3+2)
            xij(3,ipr) = x(i3+3) - x(j3+3)
            iat0(ipr) = i3
            jat0(ipr) = j3
            rij2 = xij(1,ipr)**2 + xij(2,ipr)**2 + xij(3,ipr)**2
            rm2 = 1.0d0/rij2
            rm6bar = rm6bar + rm2**3
            rm8(ipr) = rm2**4
         end do
      end do
      nsum = ipr
      fsum = nsum
      rm6bar = rm6bar/fsum
      rij = rm6bar**(-one6)
      
      ! Case where IAT1 < 0 and IAT2 < 0 (both r**-6 ave. pos.)
      
   else
      
      do jp=-iat2,-nmrat(10,indx)
      ! Modified 9/2007 by Matthew Seetin as a consequence of enabling restraints on
      ! group of atoms for angle, torsion, and my new planar restraints
         do j=nmrcom(1,jp),nmrcom(2,jp)
            j3 = 3*(j-1)
            do ip=-iat1,-nmrat(9,indx)
               do i=nmrcom(1,ip),nmrcom(2,ip)
#ifdef LES
                  if( subsp(i) == subsp(j)  .and. &
                        cnum(i) /= cnum(j) ) cycle
#endif
                  i3 = 3*(i-1)
                  ipr = ipr + 1
                  if (ipr > maxpr) goto 500
                  xij(1,ipr) = x(i3+1) - x(j3+1)
                  xij(2,ipr) = x(i3+2) - x(j3+2)
                  xij(3,ipr) = x(i3+3) - x(j3+3)
                  iat0(ipr) = i3
                  jat0(ipr) = j3
                  rij2 = xij(1,ipr)**2 + xij(2,ipr)**2 + xij(3,ipr)**2
                  rm2 = 1.0d0/rij2
                  rm6bar = rm6bar + rm2**3
                  rm8(ipr) = rm2**4
               end do
            end do
         end do
      end do
      nsum = ipr
      fsum = nsum
      rm6bar = rm6bar/fsum
      rij = rm6bar**(-one6)
      
   end if
   
   rave = rij
   
   return
   
   ! Errors:
   
   500 write(6,1000) maxpr
   call mexit(6, 1)
   
   ! Format statements:
   
   1000 format(' ERROR: Number of interaction distances in ', &
         'r**-6 averaged ', &
         'distance',/,t8,' exceeds maximum allowed by ', &
         'MAXPR = ',i5)

   
end subroutine r6ave 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine r6drv here]
subroutine r6drv(f,df)
   
   
   ! Subroutine R**-6 DeRiVatives
   
   ! This routine calcuates the derivatives associated with an r**-6
   ! average position defined group, and adds them into the force (F)
   ! array.
   
   ! The derivatives are calculated as
   
   !           dE/dr6 * dr6/dx
   
   ! where dE/dr6 is the derivative of the energy with respect to the
   ! averaged distance (calculated in routine DISNRG, and passed here as DF)
   ! and dr6/dx is the derivative of the averaged distance with repsect
   ! to each individual coordinate used in the distance.
   
   ! It is this latter distance which is calculated here. Note that
   ! most of the quantities required for calculating the derivatives are
   ! passed through common block R6AVEC. This block appears only here
   ! and in routine R6AVE.
   
   ! MAXPR is the maximum number of distances used in computing <r**-6>**-1/6.
   ! Must be changed here and in R6AVE, if desired.
   
   implicit none
   integer:: i3, iat0, ipr, j3, jat0, maxpr, nsum
   parameter (maxpr = 5000)
   _REAL_ :: df, f, fact, fact2, fsum, rij, rm6bar, rm8, xh, xij, yh, zh
   
   ! r6avec common block. Shared between this routine and R6DRV.
   
   common/r6avec/xij(3,maxpr),iat0(maxpr), &
         jat0(maxpr),rm8(maxpr),fsum,rm6bar,rij,nsum
   
   dimension f(*)
   
   fact = df*rij/(fsum*rm6bar)
   do ipr=1,nsum
      i3 = iat0(ipr)
      j3 = jat0(ipr)
      fact2 = fact*rm8(ipr)
      xh = fact2*xij(1,ipr)
      yh = fact2*xij(2,ipr)
      zh = fact2*xij(3,ipr)
      f(i3+1) = f(i3+1) - xh
      f(i3+2) = f(i3+2) - yh
      f(i3+3) = f(i3+3) - zh
      f(j3+1) = f(j3+1) + xh
      f(j3+2) = f(j3+2) + yh
      f(j3+3) = f(j3+3) + zh
   end do
   return
end subroutine r6drv 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ readname is called by the parser to read an atom name from restraint
subroutine readmask(restraint,k,iat,iatidx,name,ipres,nres,iout)

  ! Subroutine READ atom MASK
  
  ! When the parser encounters atom mask input in the restraint variable,
  ! it calls this subroutine to convert the mask into an atom number
  
  ! Authoer: Matthew G. Seetin
  ! Date: 10/2007
  
  !---------------------------
  ! INPUT VARIABLES
  !---------------------------
  
  ! restraint: a character variable of max 256 charaters that contains the
  !             natural-language input restraint
  !    iat(8) : the atoms that make up the restraint.
  !   iatidx : Index as to where in iat the next atom number read will be stored
  !     name : Array of the names of all the atoms in the system.  It may be empty
  !            if this routine is called by restal
  !     ipres: Array of all the residue numbers in the system.  May also be empty if
  !            if this routine is called by restal
  !      nres: The number of residues in the system
  !      iout: Unit number of the out file
  implicit none
  integer:: iat, iatidx, iout, ipres, k, nres
  CHARACTER(256) :: restraint
  CHARACTER(4) atname, name(*)
  INTEGER numlen,namlen,ierr,iattmp
  DIMENSION iat(8), ipres(*)
  
  numlen = 0
  namlen = 0
  iattmp = 0
  ierr = 0

  do
    if (iachar(restraint(k+numlen+1:k+numlen+1)) >= 48 &
    .AND. iachar(restraint(k+numlen+1:k+numlen+1)) <= 57) then
       numlen = numlen + 1
     else 
       exit
     end if
  end do
  if (numlen > 0 .AND. restraint(k+numlen+1:k+numlen+1) == '@') then
    do
       ! Read all characters up to the next delimiting character and assume that's the atom name
      
      if (namlen < 4 .AND. restraint(k+numlen+namlen+2:k+numlen+namlen+2) /= ',' &
      .and. restraint(k+numlen+namlen+2:k+numlen+namlen+2) /= '(' &
      .and. restraint(k+numlen+namlen+2:k+numlen+namlen+2) /= ')' &
      .and. restraint(k+numlen+namlen+2:k+numlen+namlen+2) /= '[' &
      .and. restraint(k+numlen+namlen+2:k+numlen+namlen+2) /= ']' &
      .and. restraint(k+numlen+namlen+2:k+numlen+namlen+2) /= '{' &
      .and. restraint(k+numlen+namlen+2:k+numlen+namlen+2) /= '}' &
      .and. restraint(k+numlen+namlen+2:k+numlen+namlen+2) /= ' ') then
        namlen = namlen + 1
      else
        exit
      end if
    end do
  end if
  if (numlen > 0 .AND. namlen > 0) then
    read(restraint(k+1:k+numlen), *) iat(iatidx)  ! Read residue number
    atname = restraint(k+numlen+2:k+numlen+namlen+1) ! Read atom name
     
      ! I've set restal to call this with dummy arrays for name and ipres since it doesn't
      ! have access to them.  It's not necessary to get a full, exact read of the rst file
      ! from restal (it didn't resolve atom names before), so having it receive a residue
      ! number instead of an atom number in iat is fine, as long as it counts something there.
      
    if(name(1) /= '    ' .and. ipres(1) > 0 .and. nres > 0) then
      call getnat(iat(iatidx),atname,name,ipres,nres,iout,ierr)
      if (ierr == 1) then
        write(iout,9003)
        write(iout,9999) restraint(1:LEN_TRIM(restraint))
        call mexit(iout, 1)
      end if
    end if
    iatidx = iatidx + 1
    k = k + numlen + namlen + 2
  else
    write(iout,9003)
    write(iout,9999) restraint(1:LEN_TRIM(restraint))
    call mexit(iout, 1)
  end if
  
  return
   
  9003 format('Error: Invalid residue number or atom name in restraint.')
  9999 format('restraint = "',(a),'"')
end subroutine readmask

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine tornrg here]
subroutine tornrg(x,f,dfr,ap,i1,i2,i3,i4,e,r1,r2,r3,r4,k2,k3,ntb, &
      tave1,tave01,tave2,tave02,nave,nexact,ipower, &
      tauave,ravein,dt,navint,iave,incflg,iravin, &
      iflag)
   
   
   ! Subroutine TORsion eNeRGy
   
   ! This subroutine calculates the torsion angle restraint energy.
   ! The energy function is a flat-bottomed well, with parabolic sides
   ! which become linear at greater distances.
   
   ! Author: David A. Pearlman
   ! Date: 7/89
   
   ! If r is the torsion angle:
   
   !       r <= r1  : E = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
   ! r1 <= r <= r2  : E = k2*(r-r2)**2
   ! r2 <= r <= r3  : E = 0
   ! r3 <= r <= r4  : E = k3*(r-r3)**2
   ! r4 <= r        : E = 2*k3*(r4-r3)*(r-r4) + K3*(r4-r3)**2
   
   
   ! IAVE = 0: then the instantaneous (current) value of r is used
   !           in all calculations.
   ! IAVE = 1: The time-averaged value of r is used, and forces are calculated
   !           according to the above function, i.e.
   !               dE/dx = (dE/dr_ave)(dr_ave/dr)(dr/dx).
   ! IAVE = 2: The time-averaged value of r is used, and forces are calculated as
   !               dE/dx = (dE/dr_ave)(dr/dx).
   !           This expression, if integrated would not yield the energy function
   !           above, so energies reported are pseudo-energies. Also, it is a
   !           non-conservative force, so the system will tend to heat up. But it
   !           avoids potentially large force terms with (1+IPOWER)
   !           exponentiation, which can cause the simulation to become unstable.
   
   ! When IAVE.GT.0, the averaged torsion is calculated as
   
   !        ATAN(YBAR/XBAR)    XBAR > 0
   !    180+ATAN(YBAR/XBAR)    XBAR < 0
   
   !    where XBAR = average value of cos(tau)
   !          YBAR = average value of sin(tau)
   
   ! Other input:
   !    X(I)  : Coordinate array.
   !    F(I)  : Force array;
   !            Modified on output to include forces from this restraint.
   !  I1,I2,I3: Atom pointers for this angle (3*(I-1), where I is absolute
   !        I4  atom number.
   !    NTB   : Periodic boundary conditions flag.
   !    IFLAG : =0, then calculate energy and forces.
   !            =1, then only calculate value of torsion.
   !            =2, then calculate value and energy, but do not accumulate derivatives
   !            =3, then calculate energy and derivatives, but do not
   !                accumulate forces into F() array. They are returned
   !                in the DFR array.
   !   INCFLG : Determines whether local saved pointers are updated in AVEINT
   
   ! (following only used when time-averaging:)
   ! TAVE1(I)   Contain the time-averaged values of sin(tau) using damped
   ! TAVE01(I): exponential weighting and no weighting, respectively.
   
   ! TAVE2(I):  Contain the time-averaged values of cos(tau) using damped
   ! TAVE02(I): exponential weighting and no weighting, respectively.
   
   !           Both TAVE and TAVE0 are passed so that the first element of
   !           the array corresponds to the average value of the torsion currently
   !           of interest, and the second element of the array contains
   !           the previous value of the torsion of interest.
   
   !     NAVE : > 0 for time-averaged restraings.
   !   NEXACT : Not currently used.
   !   IPOWER : The power used in averaging the torsion data.
   !   TAUAVE : The exponential decay factor used in averaging the torsion data
   !   IRAVIN : =0 do not modify the values of the torsions by RAVEIN
   !            =1 modify phi by RAVEIN, as shown below.
   !   RAVEIN : At time 0, time-average integral is undefined.
   !            If -1000<RAVEIN<1000, use r_initial+ravein on this step.
   !            If RAVEIN<=-1000, use r_target+(ravein+1000) on this step.
   !            If RAVEIN>=1000, use r_target-(ravein-1000) on this step.
   !       DT : Integration time-step (ps).
   !   NAVINT : Only store restraints in AAVE/AAVE0 every NAVINT steps.
   
   ! Output:
   !    E     : The energy contribution from this restraint.
   !    DFR   : If IFLAG=3, DFR(1->3) contains the d_E/d_distance forces
   !            for x,y, and z for position 1; those for position 2 are in
   !            DFR(4->6), and so on for each of the 5 atoms here
   !            Note that the IFLAG designation for this is different from
   !            that for disnrg, plnnrg, and plptnrg
   !   AP     : The calculated torsion angle
   
   ! Note: The torsion is translated by n*+-360, in order to bring
   !       it as close as possible to r2 or r3.
   
   use constants, only : PI, DEG_TO_RAD, half, TWOPI
   implicit none
   integer:: i1, i2, i3, i4, iave, iflag, iflg, incflg, ipower, &
        irav, iravin, m, nave, navint, nexact, ntb
   _REAL_ :: ap, apmean, bi, bk, cphi, ct, dc, denom, df, dfr, dif, &
        dif1, dr1, dr2, dr3, dr4, dr5, dr6, dravdr, drx, dry, drz, dt, &
        dx, dy, dz, e, f, gx, gy, gz, r1, r2, r3, r4, ravein, rij2, rinc, &
        rkj2, rkl2, rnow1, rnow2, s, small, small2, spc, sphi, t, tauave, &
        tave01, tave02, tave1, tave2, term1, term2, x, xij, xkj, xkl, z1, &
        z11, z12, z2, z22
   _REAL_ k2,k3
   parameter (small=1.0d-14)
   parameter (small2=1.0d-5)
   dimension x(*),f(*),dfr(*),tave1(2),tave01(2),tave2(2),tave02(2)
   dimension xij(3),xkj(3),xkl(3),t(6),dc(6)
   rinc = 1000.0d0
   
   iflg = 1
   if (iflag == 1) iflg = 2
   
   ! Calculate the torsion:
   
   rij2 = 0.0d0
   rkj2 = 0.0d0
   rkl2 = 0.0d0
   
   do m=1,3
      xij(m) = x(i1+m) - x(i2+m)
      xkj(m) = x(i3+m) - x(i2+m)
      xkl(m) = x(i3+m) - x(i4+m)
      rij2 = rij2 + xij(m)**2
      rkj2 = rkj2 + xkj(m)**2
      rkl2 = rkl2 + xkl(m)**2
   end do
   
   ! Calculate ij X jk AND kl X jk
   
   dx = xij(2)*xkj(3) - xij(3)*xkj(2)
   dy = xij(3)*xkj(1) - xij(1)*xkj(3)
   dz = xij(1)*xkj(2) - xij(2)*xkj(1)
   
   gx = xkj(3)*xkl(2) - xkj(2)*xkl(3)
   gy = xkj(1)*xkl(3) - xkj(3)*xkl(1)
   gz = xkj(2)*xkl(1) - xkj(1)*xkl(2)
   
   
   ! Calculate the magnitudes of above vectors, and their dot product:
   
   bi = dx*dx + dy*dy + dz*dz
   bk = gx*gx + gy*gy + gz*gz
   ct = dx*gx + dy*gy + dz*gz
   
   ! If this is a linear dihedral, we cannot calculate a value for it,
   ! so set energy to zero and return:
   
   if (bk < 0.01d0 .or. bi < 0.01d0) then
      e = 0.0d0
      return
   end if
   
   bi = sqrt(bi)
   bk = sqrt(bk)
   z1 = 1.0d0/bi
   z2 = 1.0d0/bk
   ct = ct*z1*z2
   
   ! Cannot let ct be exactly 1. or -1., because this results
   ! in sin(tau) = 0., and we divide by sin(tau) [sphi] in calculating
   ! the derivatives.
   
   if (ct > 1.0d0-small) ct = 1.0d0-small
   if (ct < -1.0d0+small) ct = -1.0d0+small
   ap = acos(ct)
   
   s = xkj(1)*(dz*gy-dy*gz) + xkj(2)*(dx*gz-dz*gx) + &
         xkj(3)*(dy*gx-dx*gy)
   
   if (s < 0.0d0) ap = -ap
   ap = pi - ap
   
   cphi = cos(ap)
   sphi = sin(ap)
   ct = -ct
   
   ! Get averaged value, if requested
   
   dravdr = 1.0d0
   if (iave > 0) then
      rnow1 = cphi
      rnow2 = sphi

      ! Set the value on the first step (where integral not defined):
      if (iravin == 1) then
         irav = int(ravein)
         if (-1000 < irav .and. irav < 1000) then
            rnow1 = cos(ap + ravein*DEG_TO_RAD)
            rnow2 = sin(ap + ravein*DEG_TO_RAD)
         else
            rinc = -sign(rinc,ravein)
            rinc = (ravein+rinc)*DEG_TO_RAD
            if (r2 > small2) rnow1 = cos(((r2+r3)*half) + rinc)
            if (r2 > small2) rnow2 = sin(((r2+r3)*half) + rinc)
            if (r2 <= small2) rnow1 = cos(r3 + rinc)
            if (r2 <= small2) rnow2 = sin(r3 + rinc)
         end if
      end if
      
      call aveint(tave1,tave01,rnow1,tauave,dt,navint,ipower,3, &
            incflg,iflg,cphi,denom)
      call aveint(tave2,tave02,rnow2,tauave,dt,navint,ipower,4, &
            incflg,iflg,sphi,denom)
      
      if (abs(cphi) < small .and. abs(sphi) < small) then
         e = 0.0d0
         return
      else if (abs(cphi) < small .and. sphi > 0) then
         ap = pi*half
      else if (abs(cphi) < small .and. sphi < 0) then
         ap = 3.0d0*pi*half
      else
         ap = atan(sphi/cphi)
         if (cphi < 0.0d0) ap = ap + pi
      end if
      
      ! If IAVE=1, calculate d(tau_ave)/d(tau(t)); break into different cases,
      ! depending on the value of IPOWER, to avoid singularities. If singularity
      ! found, leave DRAVDR at 1.0.
      
      if (iave == 1 .and. iflag /= 1) then
         spc = sphi**2 + cphi**2
         if (spc < small) goto 15
         if (ipower == -1) then
            dravdr = (rnow1*cphi+rnow2*sphi)*dt/(spc*denom)
         else if (ipower < -1) then
            if (abs(cphi) < small.or.abs(sphi) < small) goto 15
            term1 = rnow1*cphi*(rnow2/sphi)**(-1-ipower)
            term2 = rnow2*sphi*(rnow1/cphi)**(-1-ipower)
            dravdr = (term1+term2)*dt/(spc*denom)
         else
            if (abs(rnow1) < small.or.abs(rnow2) < small)goto 15
            term1 = rnow1*cphi*(sphi/rnow2)**(1+ipower)
            term2 = rnow2*sphi*(cphi/rnow1)**(1+ipower)
            dravdr = (term1+term2)*dt/(spc*denom)
         end if
      end if
   end if  ! (iave > 0)
   
   15 if (iflag == 1) return
   
   ! Translate the value of the torsion (by +- N*360) to bring it as close as
   ! possible to one of the two central "cutoff" points (r2,r3). Use this as
   ! the value of the torsion in the following comparison.
   
   apmean = (r2+r3)*0.5d0
   18 if (ap-apmean > pi) then
      ap = ap - TWOPI
      goto 18
   else if (apmean-ap > pi) then
      ap = ap + TWOPI
      goto 18
   end if

   
   ! Calculate energy (E) and the derivative with respect to the torsion (DF):
   
   if (ap < r1) then
      dif1 = r1-r2
      df = 2.0d0 * k2 * dif1
      e = df * (ap-r1) + k2*dif1*dif1
   else if (ap < r2) then
      dif = ap - r2
      df = 2.0d0 * k2 * dif
      e = k2*dif*dif
   else if (ap <= r3) then
      e = 0.0d0
      return
   else if (ap < r4) then
      dif = ap - r3
      df = 2.0d0 * k3 * dif
      e = k3*dif*dif
   else
      dif1 = r4-r3
      df = 2.0d0 * k3 * dif1
      e = df * (ap-r4) + k3*dif1*dif1
   end if
   if (iflag == 2) return
   
   ! Modify DF to be DE/D(cos(tau)); also multiply by DRAVDR, which is
   ! dtau_ave/dtau(t) if IAVE=1, and 1.0 otherwise:
   
   df = -dravdr*df/sphi
   
   ! Now calculate derivatives; we have already calculated DE/D{cos(tau)}.
   ! now we need D{cos(tau)}/D{X,Y,OR Z}
   
   t(1) = dx
   t(2) = dy
   t(3) = dz
   t(4) = -gx
   t(5) = -gy
   t(6) = -gz
   
   !  DC = First derivative of cos(phi) w/respect
   !  to the cartesian differences t
   
   z11 = z1*z1
   z12 = z1*z2
   z22 = z2*z2
   
   do m = 1,3
      dc(m) = t(m+3)*z12-cphi*t(m)*z11
      dc(m+3) = t(m)*z12-cphi*t(m+3)*z22
   end do
   
   dr1 = df*(dc(3)*xkj(2) - dc(2)*xkj(3))
   dr2 = df*(dc(1)*xkj(3) - dc(3)*xkj(1))
   dr3 = df*(dc(2)*xkj(1) - dc(1)*xkj(2))
   dr4 = df*(dc(6)*xkj(2) - dc(5)*xkj(3))
   dr5 = df*(dc(4)*xkj(3) - dc(6)*xkj(1))
   dr6 = df*(dc(5)*xkj(1) - dc(4)*xkj(2))
   drx = df*(-dc(2)*xij(3) + dc(3)*xij(2) + dc(5)*xkl(3) - &
         dc(6)*xkl(2))
   dry = df*( dc(1)*xij(3) - dc(3)*xij(1) - dc(4)*xkl(3) + &
         dc(6)*xkl(1))
   drz = df*(-dc(1)*xij(2) + dc(2)*xij(1) + dc(4)*xkl(2) - &
         dc(5)*xkl(1))
   
   dfr(1) =  dr1
   dfr(2) =  dr2
   dfr(3) =  dr3
   dfr(4) =  drx - dr1
   dfr(5) =  dry - dr2
   dfr(6) =  drz - dr3
   dfr(7) =  -drx - dr4
   dfr(8) =  -dry - dr5
   dfr(9) =  -drz - dr6
   dfr(10) =  dr4
   dfr(11) =  dr5
   dfr(12) =  dr6
   
   if(iflag == 3) return
   
   f(i1+1) = f(i1+1) - dr1
   f(i1+2) = f(i1+2) - dr2
   f(i1+3) = f(i1+3) - dr3
   f(i2+1) = f(i2+1) - drx + dr1
   f(i2+2) = f(i2+2) - dry + dr2
   f(i2+3) = f(i2+3) - drz + dr3
   f(i3+1) = f(i3+1) + drx + dr4
   f(i3+2) = f(i3+2) + dry + dr5
   f(i3+3) = f(i3+3) + drz + dr6
   f(i4+1) = f(i4+1) - dr4
   f(i4+2) = f(i4+2) - dr5
   f(i4+3) = f(i4+3) - dr6
   
   return
end subroutine tornrg 
!**********************************
! SUBROUTINE  MORSE by Anadra and MM (17-04-06)
subroutine smorse(natoms,De,r0,nori1,nori2,X,F)

!      call smorse(NATOM,rk2nmr(2,i),r2nmr(2,i),nmrat(1,i), &
          
        
        implicit none
        integer nori1,nori2,nat1,nat2,natoms
        double precision :: X(3*natoms),F(3*natoms),rc,dx,dy,dz, &
        Em,De,Beta,r0,Dem
        nat1=(nori1+3)/3
        nat2=(nori2+3)/3

!        write(*,*)  natoms 
!        write(*,*)  De                
!        write(*,*)  r0             
!        write(*,*)  nat1        
!        write(*,*)  nat2
        

!        write(6,*) 'MORSE SUBROUTINE WORKING'
        
!    De=21.383
    Beta=4.
!    r0=1.78


! calcula el valor de la coord de reaccion del SISTEMA
        rc=sqrt(                                &
          ( X(3*nat1-2)-X(3*nat2-2) )**2 +      &
          ( X(3*nat1-1)-X(3*nat2-1) )**2 +      &
          ( X(3*nat1)-X(3*nat2) )**2 )


! calcula la nueva energia
    
    Em=De*( exp(-2*Beta*(rc-r0)) - 2*exp(-1*Beta*(rc-r0)) )

!       if(MOD(nstep,NTWE).eq.0) then
!        write(10,'(a5,2f12.5)') 'MORSE',rc,Em
!        endif



    Dem=De*( exp(-1*Beta*(rc-r0)) - exp(-2*Beta*(rc-r0)) )
    Dem=Dem*2*Beta


          dx=(1.0/rc)*( X(3*nat1-2)-X(3*nat2-2) )
          dx=dx*Dem
          dy=(1.0/rc)*( X(3*nat1-1)-X(3*nat2-1) )
          dy=dy*Dem
          dz=(1.0/rc)*( X(3*nat1)-X(3*nat2) )
          dz=dz*Dem

! las nuevas fuerzas sobre 1 es -deriv
        F(3*nat1-2)=F(3*nat1-2)-dx
        F(3*nat2-2)=F(3*nat2-2)+dx
        F(3*nat1-1)=F(3*nat1-1)-dy
        F(3*nat2-1)=F(3*nat2-1)+dy
        F(3*nat1)=F(3*nat1)-dz
        F(3*nat2)=F(3*nat2)+dz
    
    return
end subroutine smorse

end module nmr
