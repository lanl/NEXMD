!*******************************************************************************
!
! Module: nmr_lib_mod
!
! Description: <TBS>
!              
!*******************************************************************************

module nmr_lib_mod

use file_io_dat_mod

  implicit none

contains

!*******************************************************************************
!
! Subroutine:  angnrg (ANGle eNeRGy)
!
! Description: 
!              
! This subroutine calculates the valence angle restraint energy.
! The energy function is a flat-bottomed well, with parabolic sides
! which become linear at greater distances.
!
! If r is the valence angle:
!
!       r <= r1  : energy = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
! r1 <= r <= r2  : energy = k2*(r-r2)**2
! r2 <= r <= r3  : energy = 0
! r3 <= r <= r4  : energy = k3*(r-r3)**2
! r4 <= r        : energy = 2*k3*(r4-r3)*(r-r4) + K3*(r4-r3)**2
!
! IAVE = 0: then the instantaneous (current) value of r is used
!           in all calculations. 
! IAVE = 1: The time-averaged value of r is used, and forces are calculated
!           according to the above function, i.e. 
!             denergy/dx = (denergy/dr_ave)(dr_ave/dr)(dr/dx). 
! IAVE = 2: The time-averaged value of r is used, and forces are calculated as
!             denergy/dx = (denergy/dr_ave)(dr/dx). 
!           This expression, if integrated would not yield the energy function
!           above, so energies reported are pseudo-energies. Also, it is a
!           non-conservative force, so the system will tend to heat up. But it
!           avoids potentially large force terms with (1+IPOWER)
!           exponentiation, which can cause the simulation to become unstable.
!
! Other input:
!
!    crd(i): Coordinate array.
!    frc(i): Force array;
!            Modified on output to include forces from this restraint.
!  I1,I2,I3: Atom pointers for this angle (3*(I-1), where I is absolute
!            atom number.
!    IFLAG : =0, then calculate energy and derivatives
!            =1, then only calculate current value of angle
!            =2, the calculate current value of angle and energy, but
!                do not accumulate derivatives
!   INCFLG : Determines whether local saved pointers are updated in AVEINT
!
! (the following only used when time-averaging):
! AAVE(I) : Contains time-averaged values of angles, using TAUAVE damping.
! AAVE0(I): Contains time-averaged values of angles, using no damping.
!
!           Both AAVE and AAVE0 are passed so that the first element of
!           the array corresponds to the average value of the angle currently
!           of interest, and the second element of the array contains
!           the previous value of the angle of interest.
!
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
!
! Output:
!  e       : The energy contribution from this restraint.
!  THET    : The calculated valence angle
!
! Author: David A. Pearlman
! Date: 7/89
!
!*******************************************************************************

subroutine angnrg(crd, frc, thet, i1, i2, i3, e, r1, r2, r3, r4, k2, k3, &
                  aave, aave0, ipower, tauave, ravein, dt, navint, iave, &
                  incflg, iravin, iflag)

  implicit none

! Formal arguments:

  double precision  crd(*)
  double precision  frc(*)
  double precision  thet
  integer           i1
  integer           i2
  integer           i3
  double precision  e
  double precision  r1
  double precision  r2
  double precision  r3
  double precision  r4
  double precision  k2
  double precision  k3
  double precision  aave(2)
  double precision  aave0(2)
  integer           ipower
  double precision  tauave
  double precision  ravein
  double precision  dt
  integer           navint
  integer           iave
  integer           incflg
  integer           iravin
  integer           iflag

! Local variables:

  double precision  cii
  double precision  cik
  double precision  ckk
  double precision  cst
  double precision  denom
  double precision  df
  double precision  dif
  double precision  dif1
  double precision  dravdr
  double precision  dt1
  double precision  dt2
  double precision  dt3
  double precision  dt4
  double precision  dt5
  double precision  dt6
  double precision  dt7
  double precision  dt8
  double precision  dt9
  integer           iflg
  integer           irav
  integer           m
  double precision  pi
  double precision  rdenom
  double precision  rinc
  double precision  rij
  double precision  rij2
  double precision  rkj
  double precision  rkj2
  double precision  rnow
  double precision  small
  double precision  small2
  double precision  snt
  double precision  st
  double precision  sth
  double precision  xij(3)
  double precision  xkj(3)

  parameter (pi = 3.141592653589793d0)
  parameter (small = 1.0d-14)
  parameter (small2 = 1.0d-5)

  rinc = 1000.0d0

  iflg = 1
  if (iflag .eq. 1) iflg = 2

! Calculate the angle:

  rij2 = 0.0d0
  rkj2 = 0.0d0

  do m = 1, 3
    xij(m) = crd(i1 + m) - crd(i2 + m)
    xkj(m) = crd(i3 + m) - crd(i2 + m)
    rij2 = rij2 + xij(m)**2
    rkj2 = rkj2 + xkj(m)**2
  end do

  rij = sqrt(rij2)
  rkj = sqrt(rkj2)

! Calculate angle from a*b/|ab| = cos(theta)

  rdenom = rij*rkj
  cst = (xij(1)*xkj(1) + xij(2)*xkj(2) + xij(3)*xkj(3)) / rdenom
  if (cst .gt. 1.0d0) cst = 1.0d0
  if (cst.lt. - 1.0d0) cst = - 1.0d0
  thet = acos(cst)

! Get averaged value, if requested

  dravdr = 1.0d0

  if (iave .gt. 0) then
    rnow = thet

! Set the value on the first step (where integral not defined):

    if (iravin .eq. 1) then
      irav = int(ravein)
      if (- 1000 .lt. irav .and. irav .lt. 1000) then
        rnow = rnow + ravein * pi / 180.0d0
      else
        rinc = -sign(rinc, ravein)
        rinc = (ravein + rinc) * pi / 180.0d0
        if (r2 .gt. small2) rnow = ((r2 + r3) / 2.0d0) + rinc
        if (r2 .le. small2) rnow = r3 + rinc
      end if
    end if

    call aveint(aave, aave0, rnow, tauave, dt, navint, ipower, 2, incflg, &
                iflg, thet, denom) 

! DRAVDR is the factor d(r_ave)/d(r_current)

    if (iave .eq. 1 .and. iflag .ne. 1) &
      dravdr = ((thet / rnow)**(1 + ipower))*dt / denom
      
  end if

  if (iflag .eq. 1) return

! Calculate energy (E) and the derivative with respect to the valence
! angle (DF):

  if (thet .lt. r1) then
    dif1 = r1 - r2
    df = 2.0d0 * k2 * dif1
    e = df * (thet - r1) + k2*dif1*dif1
  else if (thet .lt. r2) then
    dif = thet - r2
    df = 2.0d0 * k2 * dif
    e = k2*dif*dif
  else if (thet .le. r3) then
    e = 0.0d0
    return
  else if (thet .lt. r4) then
    dif = thet - r3
    df = 2.0d0 * k3 * dif
    e = k3*dif*dif
  else
    dif1 = r4 - r3
    df = 2.0d0 * k3 * dif1
    e = df * (thet - r4) + k3*dif1*dif1
  end if
  
  if (iflag .eq. 2) return

! Calculate the derivaties with respect to the coordinates, and add them
! into the force arrays. DRAVDR contains d(thet_ave)/d(thet(t)) if IAVE=1,
! 1.0 otherwise.

  df = dravdr * df
  snt = sin(thet)
  if (abs(snt) .lt. small) snt = small
  st = - df / snt
  sth = st*cst
  cik = st / rdenom
  cii = sth / rij2
  ckk = sth / rkj2

  dt1 = cik*xkj(1) - cii*xij(1)
  dt2 = cik*xkj(2) - cii*xij(2)
  dt3 = cik*xkj(3) - cii*xij(3)
  dt7 = cik*xij(1) - ckk*xkj(1)
  dt8 = cik*xij(2) - ckk*xkj(2)
  dt9 = cik*xij(3) - ckk*xkj(3)          
  dt4 = - dt1 - dt7
  dt5 = - dt2 - dt8
  dt6 = - dt3 - dt9

  frc(i1 + 1) = frc(i1 + 1) - dt1
  frc(i1 + 2) = frc(i1 + 2) - dt2
  frc(i1 + 3) = frc(i1 + 3) - dt3
  frc(i2 + 1) = frc(i2 + 1) - dt4
  frc(i2 + 2) = frc(i2 + 2) - dt5
  frc(i2 + 3) = frc(i2 + 3) - dt6
  frc(i3 + 1) = frc(i3 + 1) - dt7
  frc(i3 + 2) = frc(i3 + 2) - dt8
  frc(i3 + 3) = frc(i3 + 3) - dt9

  return

end subroutine angnrg

!*******************************************************************************
!
! Subroutine:  aveint (AVErage INTernals)
!
! Description:
!
! This routine evaluates and updates the integrals needed to perform 
! time-averaging
!
! Author: David A. Pearlman
! Date: 5/90
!
! Algorithm:
!     Internals are averaged using inverse weighting to the 
! ipower power (thus ipower = -1 gives straight linear weighting). In
! addition, an exponential decay function is included, so that the
! weighted average will continue to reflect the most recent values
! of the internals--even during long simulations.
!
! The formula for the time-averaged values is:
!
!                  t
!   rbar = 1/C * {int (exp(tp-t)/tau) r(tp)**-ipower dtp} **-1/ipower   (1)
!                  0
!
! where int = integral (from 0-> t)
!      rbar = average value of internal
!         t = current time
!       tau = exponential decay constant
!     r(tp) = value of internal at time tp
!    ipower = average is over internals to the inverse ipower.
!             Usually ipower = 3 for NOE distances, -1 for angles and
!             torsions
!         C = normalization integral.
!
! To minimize memory storage, we evaluate the integral as
!
!    t
!   int (exp(tp-t)/tau) r(tp)**-ipower dtp =
!    0
!
!                   t-dt
!      exp(-dt/tau) int (exp(tp-t)/tau) r(tp)**-ipower dtp +
!                    0
!
!                    t
!                   int (exp(tp-t)/tau) r(tp)**-ipower dtp
!                   t-dt
!
! INPUT:
!
!  xave(I)  : On input: xave(1) contains the value of the time-integral 
!                       shown above for the bounds 0->t-dt..
!                       xave(2) contains r(t-dt)**-ipower
!             On output (if iflag .ne. 2)
!                       xave(1) contains the value of the time-integral
!                       shown above for the bounds 0->t.
!                       xave(2) contains r(t)**-ipower
!
! xave0(I)  : Same as xave(1), except that the integral stored does not
!             reflect exponential weighting (a real time-average).
!  xnow     : Current value of the internal on this call.
!  tau      : Characteristic time for the exponential decay.
!             If tau >= 1.D+6, then the non-exponential-damped values
!             in xave0() will be used.
!  dt       : The moleular dynamics timestep.
!  navint   : The current value of the internal is only stored every
!             navint steps (typically 1).
!  ipower   : The exponent used in creating the weighted average. Standard
!             values would be 3 for distances (resulting in **-3 weighting)
!             and -1 for angles (resulting in **1/non-exponential weighting).
!  itype    : = 1 for bonds; =2 for angles; =3 for torsions (cos(tau));
!             = 4 for torsions (sin(tau)).
!  incflg   : Determines whether local saved pointers are updated
!             on this call.
!             incflg = 0 -- pointers are _not_ updated.
!             incflg = 1 -- pointers are updated.
!             Pointers should only be updated on the first call to AVEINT
!             for each restraint type (bond, angle, torsion) in any call
!             to NMRNRG.
!  iflag    : Determines behavior on this call.
!             iflag <= 1 -- Update the integral to reflect the value in xnow,
!                           set xave(2) to xnow, and return
!                           the calculated averaged value xave/denom.
!             iflag  = 2 -- Just return the calculated values xave/denom
!
! OUTPUT:
!
!  xaver    : The calculated time-averaged value for the internal.
!  denom    : The value of the exponential weight-factor divisor integral
!
! LOCAL STORAGE:
!
!  nsteps(4): The number of calls with iflag .le. 1 (i.e. calls where
!             the integral sum was updated)
!  icalls(4): Number of calls with iflag=0 or 1.
!              
!*******************************************************************************

subroutine aveint(xave, xave0, xnow, tau, dt, navint, ipower, itype, incflg, &
                  iflag, xaver, denom)

  implicit none

! Formal arguments:

  double precision  xave(2)
  double precision  xave0(2)
  double precision  xnow
  double precision  tau
  double precision  dt
  integer           navint
  integer           ipower
  integer           itype
  integer           incflg
  integer           iflag
  double precision  xaver
  double precision  denom

! Local variables:

  integer           icalls(4)
  integer           ionce
  integer           nsteps(4)

  common /avnloc/ icalls, ionce, nsteps

  double precision  dtu
  integer           i
  double precision  tp
  double precision  xnowp
  double precision  xpon

! Initialize values on first call

  if (ionce .ne. 4899) then
    do i = 1, 4
      nsteps(i) = 0
      icalls(i) = 0
    end do
    ionce = 4899
  end if

! DTU is the time-step (modified to reflect NAVINT if it is not 1)

  dtu = dt * dble(navint)

! This section stores the value

  if (iflag .le. 1) then
    icalls(itype) = icalls(itype) + incflg

    if (mod(icalls(itype) - 1, navint) .eq. 0) then
      nsteps(itype) = nsteps(itype) + incflg

! Update the integral summation using Simpsons rule to calculate the
! integral over the last interval.

! On the first step, set the integral to zero.
! XAVE(1) stores the integral with exponential decay. XAVE0(1) stores
! the integral without the exponential decay weighting factor.

      xnowp = xnow**(- ipower)
      if (nsteps(itype) .le. 1) then
        xave(1) = 0.0d0
        xave0(1) = 0.0d0
      else
        xave(1) = exp(- dtu / tau)*xave(1) +  &
                  (exp(- dtu / tau)*xave(2) + xnowp) * (dtu / 2.0d0)
        xave0(1) = xave0(1) + (xave(2) + xnowp) * (dtu / 2.0d0)
      end if
      xave(2) = xnowp
    end if
  end if

! This section calculates the average:

! TP is the current time, modified to reflect NAVINT 

  tp = dble(nsteps(itype) - 1)*dtu  

! On the first step, the current value is the average:

  if (nsteps(itype) .le. 1) then
    xaver = xnow
    denom = 1.0d0
    return
  else 
    xpon = - 1.0d0 / dble(ipower)
    if (tau .lt. 1.d+6) then
      denom = tau*(1.0d0 - exp(- tp / tau))
      xaver = sign((abs(xave(1) / denom))**xpon, xave(1) / denom)
    else
      denom = tp
      xaver = sign((abs(xave0(1) / denom))**xpon, xave0(1) / denom)
    end if
  end if

  return

end subroutine aveint

!*******************************************************************************
!
! Subroutine:  chklin (CHecK LINe)
!
! Description:
!
! This subroutine checks the character string LINE for the syntax
! MTCHST = FLNME, where MTCHST is a passed character string, and
! FLNME can be any character string. The equals sign ("=") is required.
!
! Author: David A. Pearlman
! Date: 3/90
!
! INPUT:
!
! LINE: The character line to be searched (character*80)
! MTCHST: The match string to be searched for in LINE. (Character*10).
! ILNMAT: The length of MTCHST.
! IMATCH: =0 if there is no match.
!         =1 if there is a match.
! ISTRT: The starting position of FLNME
! ILENGT: The length of FLNME.
! IOUT: Unit for error prints
!
! MAXLEN is the length of LINE
!              
!*******************************************************************************

subroutine chklin(line, mtchst, ilnmat, imatch, istrt, ilengt, iout)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  character(80)     line
  character(10)     mtchst
  integer           ilnmat
  integer           imatch
  integer           istrt
  integer           ilengt
  integer           iout

! Common variables in headers:

! Local variables:

  integer           i
  integer           j
  integer           maxlen

  parameter (maxlen = 80)

! See if the match string MTCHST occurs before the "=" sign:

  do i = 1, maxlen - ilnmat + 1
    if (line(i:i + ilnmat - 1) .eq. mtchst) go to 8
    if (line(i:i) .eq. '=') go to 6
  end do

! Get here if MTCHST not found before the "=" sign:

6 imatch = 0

  return

! Get here if MTCHST found:

8 imatch = 1

  do i = 1, maxlen
    if (line(i:i) .eq. '=') then
      do j = i + 1, maxlen
        if (line(j:j) .ne. ' ') then
          istrt = j
          go to 30
        end if
      end do
    end if
  end do

! If we get here (instead of skipping to line 30), keyword was found but
! either "=" was missing, or there was no string following the "=":

  write(iout, 9001) line
  call mexit(iout, 1)

30 do i = maxlen, 1, -1
     if (line(i:i) .ne. ' ') then
       ilengt = i - istrt + 1
       return
     end if
   end do

9001 format(' Error: Missing "=" and/or filename after keyword in line:', / , a)

end subroutine chklin

!*******************************************************************************
!
! Subroutine:   disnrg (DIStance eNeRGy)
!
! Description:
!
! This subroutine calculates the distance restraint energy.
! The energy function is a flat-bottomed well, with parabolic sides
! which become linear at greater distances.
!
! If the values in nmrat(1) or nmrat(2) are negative, then we do
! r**-6 averaging on the appropriate distances between the groups;
! group information is kept in nmrcom.
!
! If r is the bond length, and ialtd=1, then:
!
!       use a modification of a functional form suggested by 
!       Michael Nilges, Protein Eng. 2: 27-38 (1988):
!
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
!
! Otherwise, if ialtd=0, use the function form from Amber4/Amber5:
!
!           r <= r1  : E = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
!     r1 <= r <= r2  : E = k2*(r-r2)**2
!     r2 <= r <= r3  : E = 0
!     r3 <= r <= r4  : E = k3*(r-r3)**2
!     r4 <= r        : E = 2*k3*(r4-r3)*(r-r4) + K3*(r4-r3)**2
!
!
!
! Input:
!
! iave = 0: then the instantaneous (current) value of r is used
!           in all calculations. 
! iave = 1: The time-averaged value of r is used, and forces are calculated
!           according to the above function, i.e. 
!               dE/dx = (dE/dr_ave)(dr_ave/dr)(dr/dx). 
! iave = 2: The time-averaged value of r is used, and forces are calculated as
!               dE/dx = (dE/dr_ave)(dr/dx). 
!           This expression, if integrated would not yield the energy function
!           above, so energies reported are pseudo-energies. Also, it is a
!           non-conservative force, so the system will tend to heat up. But it
!           avoids potentially large force terms with (1+ipower)
!           exponentiation, which can cause the simulation to become unstable.
!
! Other input:
!
!    crd(i): Coordinate array.
!    frc(i): Force array;
!            Modified on output to include forces from this restraint.
!    i1,i2 : Atom pointers for this bond (3*(i-1), where I is absolute
!            atom number.
!    iflag : =0, then calculate energy and forces.
!            =1, then calculate only current value of bond.
!            =2, then calculate energy and d_E/d_distance, but do not
!                accumulate forces into frc() array. They are returned
!                in the dfr(6) array.
!            =3, then calculate energy and dE/d_r. Assume the distance
!                has already been calculated, and is passed in rij.
!   incflg : Not currently used.
!
! (the following only used when time-averaging):
!
! bave(i) : Contains time-averaged values of bonds, using tauave damping.
! bave0(i): Contains time-averaged values of bonds, using no damping.
!
!           Both bave and bave0 are passed so that the first element of
!           the array corresponds to the average value of the bond currently
!           of interest, and the second element of the array contains
!           the previous value of the bond of interest.
!
!   ipower : The power used in averaging the distance data.
!   tauave : The exponential decay factor used in averaging the distance data
!   iravin : =0 do not modify the values of the distances by ravein
!            =1 modify r by ravein, as shown below.
!   ravein : At time 0, time-average integral is undefined. 
!            If -1000<ravein<1000, use r_initial+ravein on this step.
!            If ravein<=-1000, use r_target+(ravein+1000) on this step.
!            If ravein>=1000, use r_target-(ravein-1000) on this step.
!       dt : The integration time-step (ps)
!   navint : Time-averaged restraints are only updated every navint
!            steps. Typically, navint=1.
!
! Output:
!
!    e     : The energy contribution from this restraint.
!    dfr   : If iflag=2, dfr(1->3) contains the d_E/d_distance forces
!            for x,y, and z for position 1; those for position 2 are in
!            dfr(4->6).
!    rij   : The actual bond distance
!
! Author: David A. Pearlman
! Date: 7/89
!
!
!*******************************************************************************

subroutine disnrg(crd, frc, dfr, rij, i1, i2, e, r1, r2, r3, r4, k2, k3, &
                  bave, bave0, ipower, tauave, ravein, dt, navint, iave, &
                  incflg, iravin, iflag, ialtd)

  implicit none

! Formal arguments:

  double precision  crd(*)
  double precision  frc(*)
  double precision  dfr(6)
  double precision  rij
  integer           i1
  integer           i2
  double precision  e
  double precision  r1
  double precision  r2
  double precision  r3
  double precision  r4
  double precision  k2
  double precision  k3
  double precision  bave(2)
  double precision  bave0(2)
  integer           ipower
  double precision  tauave
  double precision  ravein
  double precision  dt
  integer           navint
  integer           iave
  integer           incflg
  integer           iravin
  integer           iflag
  integer           ialtd

! Local variables:

  double precision  a
  double precision  b
  double precision  denom
  double precision  df
  double precision  dif
  double precision  dif1
  double precision  dravdr
  integer           iflg
  integer           irav
  integer           m
  double precision  rij2
  double precision  rinc
  double precision  rnow
  double precision  small
  double precision  xij(3)

  parameter (small = 1.0d-5)
  rinc = 1000.0d0

  iflg = 1
  if (iflag .eq. 1) iflg = 2

! Calculate distance:

  if (iflag .ne. 3) then
    rij2 = 0.0d0
    do m = 1, 3
      xij(m) = crd(i1 + m) - crd(i2 + m)
      rij2 = rij2 + xij(m)**2
    end do
    rij = sqrt(rij2)
  end if

! Get averaged value, if requested

  dravdr = 1.0d0

  if (iave .gt. 0) then
    rnow = rij

! Set the value on the first step (where the integral is not defined):

    if (iravin .eq. 1) then
      irav = int(ravein)
      if (- 1000 .lt. irav .and. irav .lt. 1000) then
        rnow = rnow + ravein
      else
        rinc = - sign(rinc, ravein)
        if (r2 .gt. small) rnow = ((r2 + r3) / 2.0d0) + ravein + rinc
        if (r2 .le. small) rnow = r3 + ravein + rinc
      end if
    end if
    call aveint(bave, bave0, rnow, tauave, dt, navint, ipower, 1, &
                  incflg, iflg, rij, denom)

! DRAVDR is the factor d(r_ave)/d(r_current)

    if (iave .eq. 1 .and. iflag .ne. 1)  &
      dravdr = ((rij / rnow)**(1 + ipower))*dt / denom
  end if

  if (iflag .eq. 1) return

! Calculate energy (E) and the derivative with respect to the bond
! length (DF):

  if (ialtd .eq. 1) then

    if (rij .lt. r2) then
      dif = rij - r2
      df = 4.0d0 * k2 * dif**3
      e = k2*dif**4
    else if (rij .le. r3) then
      e = 0.0d0
      do m = 1, 6
        dfr(m) = 0.0d0
      end do
      return
    else if (rij .lt. r4) then
      dif = rij - r3
      df = 2.0d0 * k3 * dif
      e = k3*dif*dif
    else
      dif1 = r4 - r3

! older version, with c = r1; dangerous if input files are not
!           correctly modified!

!         A = 3.D0*DIF1**2 - 2.D0*R1*DIF1
!         B = -2.D0*DIF1**3 + R1*DIF1**2
!         DIF = RIJ - R3
!         E = K3*(R1*DIF + B/DIF + A)
!         DF = K3*(R1 - B/DIF**2)

! following lines assume c=0, so potential becomes flat at
!            large distances:

      a = 3.d0*dif1**2
      b = - 2.d0*dif1**3
      dif = rij - r3
      e = k3*(b / dif + a)
      df = - k3*(b / dif**2)

    end if

  else

    if (rij .lt. r1) then
      dif1 = r1 - r2
      df = 2.0d0 * k2 * dif1
      e = df * (rij - r1) + k2*dif1*dif1
    else if (rij .lt. r2) then
      dif = rij - r2
      df = 2.0d0 * k2 * dif
      e = k2*dif*dif
    else if (rij .le. r3) then
      e = 0.0d0
      do m = 1, 6
        dfr(m) = 0.0d0
      end do
      return
    else if (rij .lt. r4) then
      dif = rij - r3
      df = 2.0d0 * k3 * dif
      e = k3*dif*dif
    else
      dif1 = r4 - r3
      df = 2.0d0 * k3 * dif1
      e = df * (rij - r4) + k3*dif1*dif1
    end if

  end if

! Calculate the derivaties with respect to the coordinates, and add them
! into the force arrays. If IFLAG=2, do not add them into the force arrays
! (probably to be used in center-of-mass restraints). If IFLAG=3, return
! with just dE/dr force in DFR(1).

! DRAVDR contains dr_ave/dr(t) if IAVE=1, 1.0 otherwise.

  if (iflag .eq. 3) then
    dfr(1) = df
    return
  end if

  df = dravdr*df / rij

  do m = 1, 3
    dfr(m) = xij(m) * df
    dfr(m + 3) = -xij(m) * df
  end do

  if (iflag .eq. 2) return

  do m = 1, 3
    frc(i1 + m) = frc(i1 + m) - dfr(m)
    frc(i2 + m) = frc(i2 + m) + dfr(m)
  end do

  return

end subroutine disnrg

!*******************************************************************************
!
! Subroutine:   jnrg    (J coupling eNeRGy)
!
! Description:
!
! This routine calculates the restraint energy due to a J-Coupling
! restraint.
!
! J is calculated using a Karplus equation expressed in terms of
! of torsion angle cosines:
!
!       J(tau) = A cos(tau)^2 + B cos(tau) + C
!
! The energy function is a flat-bottomed well, with parabolic sides
! which become linear at greater distances.
!
! Author: David A. Pearlman
! Date: 3/93
!
! If r is the torsion angle:
!
!       r <= r1  : E = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
! r1 <= r <= r2  : E = k2*(r-r2)**2
! r2 <= r <= r3  : E = 0
! r3 <= r <= r4  : E = k3*(r-r3)**2
! r4 <= r        : E = 2*k3*(r4-r3)*(r-r4) + K3*(r4-r3)**2
!
! Input:
!
! iave = 0: then the instantaneous (current) value of r is used
!           in all calculations. 
! iave = 1: The time-averaged value of r is used, and forces are calculated
!           according to the above function, i.e. 
!               dE/dx = (dE/dr_ave)(dr_ave/dr)(dr/dx). 
! iave = 2: The time-averaged value of r is used, and forces are calculated as
!               dE/dx = (dE/dr_ave)(dr/dx). 
!           This expression, if integrated would not yield the energy function
!           above, so energies reported are pseudo-energies. Also, it is a
!           non-conservative force, so the system will tend to heat up. But it
!           avoids potentially large force terms with (1+ipower)
!           exponentiation, which can cause the simulation to become unstable.
!
! Other input:
!
!    crd(i): Coordinate array.
!    frc(i): Force array;
!            Modified on output to include forces from this restraint.
!  i1,i2,i3: Atom pointers for this angle (3*(I-1), where I is absolute
!        i4  atom number.
!    iflag : =0, then calculate energy and forces.
!            =1, then only calculate value of torsion.
!            =2, then calculate value and energy, but do not accumulate
!                derivatives.
!   incflg : Determines whether local saved pointers are updated in aveint
!
! (following only used when time-averaging:)
!
! tave1(i)   Contain the time-averaged values of J(tau) using damped
! tave01(i): exponential weighting and no weighting, respectively.
!            
!           Both tave and tave0 are passed so that the first element of
!           the array corresponds to the average value of the J(tau) currently
!           of interest, and the second element of the array contains
!           the previous value of the J(tau) of interest.
!
!   ipower : The power used in averaging the J(tau) data.
!   tauave : The exponential decay factor used in averaging the J(tau) data
!   iravin : =0 do not modify the values of the initial J(tau) by ravein
!            =1 modify J(tau) by ravein, as shown below.
!   ravein : At time 0, time-average integral is undefined. 
!            If -1000<ravein<1000, use r_initial+ravein on this step.
!            If ravein<=-1000, use r_target+(ravein+1000) on this step.
!            If ravein>=1000, use r_target-(ravein-1000) on this step.
!       dt : Integration time-step (ps).
!   navint : Only store restraints in tave/tave0 every navint steps.
!
! ajcoef(3): The 3 coefficients in the Karplus Equation for this torsion.
!
! Output:
!
!    e     : The energy contribution from this restraint.
!    ajval : The calculated J(tau)
!
!*******************************************************************************

subroutine jnrg(crd, frc, ajval, i1, i2, i3, i4, e, r1, r2, r3, r4, k2, k3, &
                tave1, tave01, ipower, tauave, ravein, dt, navint, iave, &
                incflg, iravin, ajcoef, iflag)

  implicit none

! Formal arguments:

  double precision  crd(*)
  double precision  frc(*)
  double precision  ajval
  integer           i1
  integer           i2
  integer           i3
  integer           i4
  double precision  e
  double precision  r1
  double precision  r2
  double precision  r3
  double precision  r4
  double precision  k2
  double precision  k3
  double precision  tave1(2)
  double precision  tave01(2)
  integer           ipower
  double precision  tauave
  double precision  ravein
  double precision  dt
  integer           navint
  integer           iave
  integer           incflg
  integer           iravin
  double precision  ajcoef(3)
  integer           iflag

! Local variables:

  double precision  ajnow
  double precision  amlt
  double precision  ap
  double precision  bi
  double precision  bk
  double precision  cphi
  double precision  ct
  double precision  dc(6)
  double precision  denom
  double precision  df
  double precision  dif
  double precision  dif1
  double precision  dr1, dr2, dr3, dr4, dr5, dr6
  double precision  dravdr
  double precision  drx, dry, drz
  double precision  dx, dy, dz
  double precision  gx, gy, gz
  integer           iflg
  integer           irav
  integer           m
  double precision  pi
  double precision  rij2
  double precision  rkj2
  double precision  rkl2
  double precision  rinc
  double precision  s
  double precision  small
  double precision  small2
  double precision  t(6)
  double precision  xij(3)
  double precision  xkj(3)
  double precision  xkl(3)
  double precision  z1, z11, z12, z2, z22

  parameter (small = 1.0d-14)
  parameter (small2 = 1.0d-5)
  parameter (pi = 3.141592653589793d0)

  rinc = 1000.0d0

  iflg = 1
  if (iflag .eq. 1) iflg = 2

! Calculate the underlying torsion:

  rij2 = 0.0d0
  rkj2 = 0.0d0
  rkl2 = 0.0d0

  do m = 1, 3
    xij(m) = crd(i1 + m) - crd(i2 + m)
    xkj(m) = crd(i3 + m) - crd(i2 + m)
    xkl(m) = crd(i3 + m) - crd(i4 + m)
    rij2 = rij2 + xij(m)**2
    rkj2 = rkj2 + xkj(m)**2
    rkl2 = rkl2 + xkl(m)**2
  end do

! Calculate ij x jk and kl x jk:

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

  if (bk .lt. 0.01d0 .or. bi .lt. 0.01d0) then
    e = 0.0d0
    return
  end if

  bi = sqrt(bi)
  bk = sqrt(bk)
  z1 = 1.0d0/bi
  z2 = 1.0d0/bk
  ct = ct*z1*z2

! Cannot let ct be exactly 1. or -1., because this results in sin(tau) = 0.,
! and we divide by sin(tau) [sphi] in calculating the derivatives.

  if (ct .gt. 1.0d0 - small) ct = 1.0d0 - small
  if (ct .lt. - 1.0d0 + small) ct = - 1.0d0 + small
  ap = acos(ct)

  s = xkj(1)*(dz*gy - dy*gz) + xkj(2)*(dx*gz - dz*gx) + xkj(3)*(dy*gx - dx*gy)

  if (s .lt. 0.0d0) ap = -ap
  ap = pi - ap

  cphi = cos(ap)
  ct = - ct

! Define the current value of the karplus-equation-calculated J value.
! If time-averaging is being performed, ajval will be replaced by its
! time-averaged value in routine aveint:

  ajnow = ajcoef(1) * cphi**2 + ajcoef(2)*cphi + ajcoef(3)
  ajval = ajnow

! Get averaged value, if requested:

  dravdr = 1.0d0

  if (iave .gt. 0) then

! Set the value on the first step (where integral not defined):

    if (iravin .eq. 1) then
      irav = int(ravein)
      if (- 1000 .lt. irav .and. irav .lt. 1000) then
        ajnow = ajnow + ravein
      else
        rinc = - sign(rinc, ravein)
        if (r2 .gt. small2) ajnow = ((r2 + r3)/2.0d0) + ravein + rinc
        if (r2 .le. small2) ajnow = r3 + ravein + rinc
      end if
    end if

    call aveint(tave1, tave01, ajnow, tauave, dt, navint, ipower, 3, incflg, &
                iflg, ajval, denom) 

! dravdr is the factor d(<j>)/d(j(t)):

    if (iave .eq. 1 .and. iflag .ne. 1) &
        dravdr = ((ajval/ajnow)**(1 + ipower))*dt/denom
  end if

  if (iflag .eq. 1) return

! Calculate energy (e) and the derivative with respect to the torsion (df):

  if (ajval .lt. r1) then
    dif1 = r1 - r2
    df = 2.0d0 * k2 * dif1
    e = df * (ajval - r1) + k2*dif1*dif1
  else if (ajval .lt. r2) then
    dif = ajval - r2
    df = 2.0d0 * k2 * dif
    e = k2*dif*dif
  else if (ajval .le. r3) then
    e = 0.0d0
    return
  else if (ajval .lt. r4) then
    dif = ajval - r3
    df = 2.0d0 * k3 * dif
    e = k3*dif*dif
  else
    dif1 = r4 - r3
    df = 2.0d0 * k3 * dif1
    e = df * (ajval - r4) + k3*dif1*dif1
  end if

  if (iflag .eq. 2) return

! Calculate the total forces by the chain rule:
!
!      dE/dx = dE/d<J> * d<J>/dJ(t) * dJ(t)/dcos(tau) * dcos(tau)/dx
!
! Multiply df, which is de/d<j> by dravdr, which is
! d<J>/d[J(t)] if iave=1, and 1.0 otherwise:

  df = dravdr*df

! Multiply df by d[J(t)]/d[cos(tau)]

  amlt = 2.0d0*ajcoef(1)*cphi + ajcoef(2)
  df = df * amlt

! Now calculate derivatives; we have already calculated df = DE/D{cos(tau)};
! now we need D{cos(tau)}/D{X,Y,OR Z}

  t(1) = dx
  t(2) = dy
  t(3) = dz
  t(4) = -gx
  t(5) = -gy
  t(6) = -gz

! dc = First derivative of cos(phi) w/respect to the cartesian differences t.

  z11 = z1*z1
  z12 = z1*z2
  z22 = z2*z2

  do m = 1, 3
      dc(m) = t(m + 3)*z12 - cphi*t(m)*z11
      dc(m + 3) = t(m)*z12 - cphi*t(m + 3)*z22
  end do

  dr1 = df*(dc(3)*xkj(2) - dc(2)*xkj(3))
  dr2 = df*(dc(1)*xkj(3) - dc(3)*xkj(1))
  dr3 = df*(dc(2)*xkj(1) - dc(1)*xkj(2))
  dr4 = df*(dc(6)*xkj(2) - dc(5)*xkj(3))
  dr5 = df*(dc(4)*xkj(3) - dc(6)*xkj(1))
  dr6 = df*(dc(5)*xkj(1) - dc(4)*xkj(2))
  drx = df*(- dc(2)*xij(3) + dc(3)*xij(2) + dc(5)*xkl(3) - dc(6)*xkl(2))
  dry = df*(dc(1)*xij(3) - dc(3)*xij(1) - dc(4)*xkl(3) + dc(6)*xkl(1))
  drz = df*(- dc(1)*xij(2) + dc(2)*xij(1) + dc(4)*xkl(2) - dc(5)*xkl(1))

  frc(i1 + 1) = frc(i1 + 1) - dr1
  frc(i1 + 2) = frc(i1 + 2) - dr2
  frc(i1 + 3) = frc(i1 + 3) - dr3
  frc(i2 + 1) = frc(i2 + 1) - drx + dr1
  frc(i2 + 2) = frc(i2 + 2) - dry + dr2
  frc(i2 + 3) = frc(i2 + 3) - drz + dr3
  frc(i3 + 1) = frc(i3 + 1) + drx + dr4
  frc(i3 + 2) = frc(i3 + 2) + dry + dr5
  frc(i3 + 3) = frc(i3 + 3) + drz + dr6
  frc(i4 + 1) = frc(i4 + 1) - dr4
  frc(i4 + 2) = frc(i4 + 2) - dr5
  frc(i4 + 3) = frc(i4 + 3) - dr6

  return

end subroutine jnrg

!*******************************************************************************
!
! Subroutine:   tornrg (TORsion eNeRGy)
!
! Description:
!
! This subroutine calculates the torsion angle restraint energy.
! The energy function is a flat-bottomed well, with parabolic sides
! which become linear at greater distances.
!
! Author: David A. Pearlman
! Date: 7/89
!
! If r is the torsion angle:
!
!       r <= r1  : e = 2*k2*(r1-r2)*(r-r1) + k2*(r1-r2)**2
! r1 <= r <= r2  : e = k2*(r-r2)**2
! r2 <= r <= r3  : e = 0
! r3 <= r <= r4  : e = k3*(r-r3)**2
! r4 <= r        : e = 2*k3*(r4-r3)*(r-r4) + k3*(r4-r3)**2
!
!
! iave = 0: then the instantaneous (current) value of r is used
!           in all calculations. 
! iave = 1: The time-averaged value of r is used, and forces are calculated
!           according to the above function, i.e. 
!               dE/dx = (dE/dr_ave)(dr_ave/dr)(dr/dx). 
! iave = 2: The time-averaged value of r is used, and forces are calculated as
!               dE/dx = (dE/dr_ave)(dr/dx). 
!           This expression, if integrated would not yield the energy function
!           above, so energies reported are pseudo-energies. Also, it is a
!           non-conservative force, so the system will tend to heat up. But it
!           avoids potentially large force terms with (1+ipower)
!           exponentiation, which can cause the simulation to become unstable.
!
! when iave .gt. 0, the averaged torsion is calculated as
!
!        atan(ybar/xbar)    xbar > 0
!    180+atan(ybar/xbar)    xbar < 0
!
!    where xbar = average value of cos(tau)
!          ybar = average value of sin(tau)
!
! Other input:
!    crd(i): Coordinate array.
!    frc(i): Force array;
!            Modified on output to include forces from this restraint.
!  i1,i2,i3: Atom pointers for this angle (3*(i-1), where i is absolute
!        i4  atom number.
!    iflag : =0, then calculate energy and forces.
!            =1, then only calculate value of torsion.
!            =2, then calculate value and energy, but do not accumulate
!            derivatives.
!   incflg : Determines whether local saved pointers are updated in aveint
!
! (following only used when time-averaging:)
! tave1(i)   Contain the time-averaged values of sin(tau) using damped
! tave01(i): exponential weighting and no weighting, respectively.
!            
! tave2(i):  Contain the time-averaged values of cos(tau) using damped
! tave02(i): exponential weighting and no weighting, respectively.
!
!           Both tave and tave0 are passed so that the first element of
!           the array corresponds to the average value of the torsion currently
!           of interest, and the second element of the array contains
!           the previous value of the torsion of interest.
!
!   ipower : The power used in averaging the torsion data.
!   tauave : The exponential decay factor used in averaging the torsion data
!   iravin : =0 do not modify the values of the torsions by ravein
!            =1 modify phi by ravein, as shown below.
!   ravein : At time 0, time-average integral is undefined. 
!            If -1000<ravein<1000, use r_initial+ravein on this step.
!            If ravein<=-1000, use r_target+(ravein+1000) on this step.
!            If ravein>=1000, use r_target-(ravein-1000) on this step.
!       dt : Integration time-step (ps).
!   navint : Only store restraints in AAVE/AAVE0 every NAVINT steps.
!
! Output:
!    e     : The energy contribution from this restraint.
!   ap     : The calculated torsion angle
!
! Note: The torsion is translated by n*+-360, in order to bring
!       it as close as possible to r2 or r3.
!
!*******************************************************************************

subroutine tornrg(crd, frc, ap, i1, i2, i3, i4, e, r1, r2, r3, r4, k2, k3, &
                  tave1, tave01, tave2, tave02, ipower, tauave, ravein, dt, &
                  navint, iave, incflg, iravin, iflag)

  implicit none

! Formal arguments:

  double precision  crd(*)
  double precision  frc(*)
  double precision  ap
  integer           i1, i2, i3, i4
  double precision  e
  double precision  r1, r2, r3, r4
  double precision  k2
  double precision  k3
  double precision  tave1(2)
  double precision  tave01(2)
  double precision  tave2(2)
  double precision  tave02(2)
  integer           ipower
  double precision  tauave
  double precision  ravein
  double precision  dt
  integer           navint
  integer           iave
  integer           incflg
  integer           iravin
  integer           iflag

! Local variables:

  double precision  apmean
  double precision  bi
  double precision  bk
  double precision  cphi
  double precision  ct
  double precision  dc(6)
  double precision  denom
  double precision  df
  double precision  dif
  double precision  dif1
  double precision  dr1, dr2, dr3, dr4, dr5, dr6
  double precision  dravdr
  double precision  drx, dry, drz
  double precision  dx, dy, dz
  double precision  gx
  double precision  gy
  double precision  gz
  integer           iflg
  integer           irav
  integer           m
  double precision  pi
  double precision  pi2
  double precision  rinc
  double precision  rij2
  double precision  rkj2
  double precision  rkl2
  double precision  rnow1
  double precision  rnow2
  double precision  s
  double precision  small
  double precision  small2
  double precision  spc
  double precision  sphi
  double precision  t(6)
  double precision  term1
  double precision  term2
  double precision  xij(3)
  double precision  xkj(3)
  double precision  xkl(3)
  double precision  z1
  double precision  z11
  double precision  z12
  double precision  z2
  double precision  z22

  parameter (small = 1.0d-14)
  parameter (small2 = 1.0d-5)
  parameter (pi = 3.141592653589793d0)
  parameter (pi2 = 2.0d0*pi)
  rinc = 1000.0d0

  iflg = 1
  if (iflag .eq. 1) iflg = 2

! Calculate the torsion:

  rij2 = 0.0d0
  rkj2 = 0.0d0
  rkl2 = 0.0d0

  do m = 1, 3
    xij(m) = crd(i1 + m) - crd(i2 + m)
    xkj(m) = crd(i3 + m) - crd(i2 + m)
    xkl(m) = crd(i3 + m) - crd(i4 + m)
    rij2 = rij2 + xij(m)**2
    rkj2 = rkj2 + xkj(m)**2
    rkl2 = rkl2 + xkl(m)**2
  end do

! Calculate ij X jk AND kl X jk:

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

  if (bk .lt. 0.01d0 .or. bi .lt. 0.01d0) then
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

  if (ct .gt. 1.0d0 - small) ct = 1.0d0 - small
  if (ct .lt. - 1.0d0 + small) ct = - 1.0d0 + small
  ap = acos(ct)

  s = xkj(1)*(dz*gy - dy*gz) + xkj(2)*(dx*gz - dz*gx) + xkj(3)*(dy*gx - dx*gy)

  if (s .lt. 0.0d0) ap = -ap
  ap = pi - ap

  cphi = cos(ap)
  sphi = sin(ap)
  ct = - ct

! Get averaged value, if requested:

  dravdr = 1.0d0

  if (iave .gt. 0) then
    rnow1 = cphi
    rnow2 = sphi

! Set the value on the first step (where integral not defined):

    if (iravin .eq. 1) then
      irav = int(ravein)
      if (-1000 .lt. irav .and. irav .lt. 1000) then
        rnow1 = cos(ap + ravein*pi/180.0d0)
        rnow2 = sin(ap + ravein*pi/180.0d0)
      else
        rinc = - sign(rinc, ravein)
        rinc = (ravein + rinc)*pi/180.0d0
        if (r2 .gt. small2) rnow1 = cos(((r2 + r3)/2.0d0) + rinc)
        if (r2 .gt. small2) rnow2 = sin(((r2 + r3)/2.0d0) + rinc)
        if (r2 .le. small2) rnow1 = cos(r3 + rinc)
        if (r2 .le. small2) rnow2 = sin(r3 + rinc)
      end if
    end if

    call aveint(tave1, tave01, rnow1, tauave, dt, navint, ipower, 3,  &
                incflg, iflg, cphi, denom) 
    call aveint(tave2, tave02, rnow2, tauave, dt, navint, ipower, 4,  &
                incflg, iflg, sphi, denom) 

    if (abs(cphi) .lt. small .and. abs(sphi) .lt. small) then
      e = 0.0d0
      return
    else if (abs(cphi) .lt. small .and. sphi .gt. 0) then
      ap = pi/2.0d0
    else if (abs(cphi) .lt. small .and. sphi .lt. 0) then
      ap = 3.0d0*pi/2.0d0
    else
      ap = atan(sphi/cphi)
      if (cphi .lt. 0.0d0) ap = ap + pi
    end if

! If iave=1, calculate d(tau_ave)/d(tau(t)); break into different cases,
! depending on the value of ipower, to avoid singularities. If singularity
! found, leave dravdr at 1.0. 

    if (iave .eq. 1 .and. iflag .ne. 1) then

      spc = sphi**2 + cphi**2

      if (spc .lt. small) goto 15

      if (ipower.eq. - 1) then

        dravdr = (rnow1*cphi + rnow2*sphi)*dt/(spc*denom)

      else if (ipower.lt. - 1) then

        if (abs(cphi) .lt. small .or. abs(sphi) .lt. small) goto 15

        term1 = rnow1*cphi*(rnow2/sphi)**(- 1 - ipower)
        term2 = rnow2*sphi*(rnow1/cphi)**(- 1 - ipower)
        dravdr = (term1 + term2)*dt/(spc*denom)

      else

        if (abs(rnow1) .lt. small .or. abs(rnow2) .lt. small) goto 15

        term1 = rnow1*cphi*(sphi/rnow2)**(1 + ipower)
        term2 = rnow2*sphi*(cphi/rnow1)**(1 + ipower)
        dravdr = (term1 + term2)*dt/(spc*denom)
      end if

    end if
  end if

15 if (iflag .eq. 1) return

! Translate the value of the torsion (by +- n*360) to bring it as close as
! possible to one of the two central "cutoff" points (r2,r3). Use this as
! the value of the torsion in the following comparison.

  apmean = (r2 + r3)*0.5

18 if (ap - apmean .gt. pi) then
     ap = ap - pi2
     goto 18
   else if (apmean - ap .gt. pi) then
     ap = ap + pi2
     goto 18
   end if 

! Calculate energy (e) and the derivative with respect to the torsion (df):

  if (ap .lt. r1) then
    dif1 = r1 - r2
    df = 2.0d0 * k2 * dif1
    e = df * (ap - r1) + k2*dif1*dif1
  else if (ap .lt. r2) then
    dif = ap - r2
    df = 2.0d0 * k2 * dif
    e = k2*dif*dif
  else if (ap .le. r3) then
    e = 0.0d0
    return
  else if (ap .lt. r4) then
    dif = ap - r3
    df = 2.0d0 * k3 * dif
    e = k3*dif*dif
  else
    dif1 = r4 - r3
    df = 2.0d0 * k3 * dif1
    e = df * (ap - r4) + k3*dif1*dif1
  end if

  if (iflag .eq. 2) return

! Modify df to be de/d(cos(tau)); also multiply by dravdr, which is
! dtau_ave/dtau(t) if iave=1, and 1.0 otherwise:

  df = - dravdr*df/sphi

! Now calculate derivatives; we have already calculated de/d{cos(tau)}.
! now we need d{cos(tau)}/d{x,y,or z}

  t(1) = dx
  t(2) = dy
  t(3) = dz
  t(4) = -gx
  t(5) = -gy
  t(6) = -gz

! dc = First derivative of cos(phi) w/respect to the cartesian differences t.

  z11 = z1*z1
  z12 = z1*z2
  z22 = z2*z2

  do m = 1, 3
      dc(m) = t(m + 3)*z12 - cphi*t(m)*z11
      dc(m + 3) = t(m)*z12 - cphi*t(m + 3)*z22
  end do

  dr1 = df*(dc(3)*xkj(2) - dc(2)*xkj(3))
  dr2 = df*(dc(1)*xkj(3) - dc(3)*xkj(1))
  dr3 = df*(dc(2)*xkj(1) - dc(1)*xkj(2))
  dr4 = df*(dc(6)*xkj(2) - dc(5)*xkj(3))
  dr5 = df*(dc(4)*xkj(3) - dc(6)*xkj(1))
  dr6 = df*(dc(5)*xkj(1) - dc(4)*xkj(2))
  drx = df*(- dc(2)*xij(3) + dc(3)*xij(2) + dc(5)*xkl(3) -  dc(6)*xkl(2))
  dry = df*(dc(1)*xij(3) - dc(3)*xij(1) - dc(4)*xkl(3) +  dc(6)*xkl(1))
  drz = df*(- dc(1)*xij(2) + dc(2)*xij(1) + dc(4)*xkl(2) -  dc(5)*xkl(1))

  frc(i1 + 1) = frc(i1 + 1) - dr1
  frc(i1 + 2) = frc(i1 + 2) - dr2
  frc(i1 + 3) = frc(i1 + 3) - dr3
  frc(i2 + 1) = frc(i2 + 1) - drx + dr1
  frc(i2 + 2) = frc(i2 + 2) - dry + dr2
  frc(i2 + 3) = frc(i2 + 3) - drz + dr3
  frc(i3 + 1) = frc(i3 + 1) + drx + dr4
  frc(i3 + 2) = frc(i3 + 2) + dry + dr5
  frc(i3 + 3) = frc(i3 + 3) + drz + dr6
  frc(i4 + 1) = frc(i4 + 1) - dr4
  frc(i4 + 2) = frc(i4 + 2) - dr5
  frc(i4 + 3) = frc(i4 + 3) - dr6

  return

end subroutine tornrg

!*******************************************************************************
!
! Subroutine:   nmrcmf (NMR Center of Mass Forces)
!
! This routine calculates the correct -d_E/d_position forces for the atoms
! which are involved in a distance restraint between two center of mass
! positions. The force is given by d(E)/d(x_com) * d(x_com)/d(x_atom).
! The df array contains the d(E)/d(x_com) contributions, as calculated in
! routine disnrg.
!
! Description:
!              
! Author: David A. Pearlman
! Date: 8/89
!
! INPUT:
!
! frc(i):     The force array.
!           Modified on output to reflect the forces from this restraint.
! dfr(6):   d(E)/d(x_com),d(E)/d(y_com) and d(E)/d(z_com) for position one
!           (elements 1->3) and position two (elements 4->6).
! nmrat(4,i)
! nmrcom(2,i): Define the group of atoms used to calculate the center of mass.
!              If nmrat(1,i) < 0, then
!                 (nmrcom(1,k) -> nmrcom(2,k)), k=-nmrat(1,i),-nmrat(3,i)
!              gives the atoms in the group whose center of mass is position 1.
!
!              If nmrat(2,i) < 0, then
!                 (nmrcom(1,k) -> nmrcom(2,k)), k=-nmrat(2,i),-nmrat(4,i)
!              gives the atoms in the group whose center of mass is position 2.
! rmstot(2): The total masses for all the atoms defining the group at each
!            end of the distance restraint.
! rimass(i): The inverse mass for atom I.
!
! i: The number of the restraint in the list.
!
!*******************************************************************************

subroutine nmrcmf(frc, dfr, nmrat, nmrcom, rmstot, i)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision  frc(*)
  double precision  dfr(6)
  integer           nmrat(4, 1)
  integer           nmrcom(2, 1)
  double precision  rmstot(2)
  integer           i

! Local variables:

  double precision  dcomdx
  integer           iat
  integer           ip
  integer           j
  integer           jatm
  integer           jx
  integer           m
  double precision  rmass

! Loop over both ends of the distance restraint:

  do iat = 1, 2

    jatm = 3*(iat - 1)

! Calculate the d(c.o.m.)/d(coordinate) derivatives, multiply them
! by the derivatives in df, and add them into the force array.

    if (nmrat(iat, i) .lt. 0) then
      do ip = - nmrat(iat, i), - nmrat(2 + iat, i)
        do j = nmrcom(1, ip), nmrcom(2, ip)
          jx = 3*(j - 1)
          rmass = atm_mass(j)
          dcomdx = rmass/rmstot(iat)
          frc(jx + 1) = frc(jx + 1) - dfr(jatm + 1)*dcomdx
          frc(jx + 2) = frc(jx + 2) - dfr(jatm + 2)*dcomdx
          frc(jx + 3) = frc(jx + 3) - dfr(jatm + 3)*dcomdx
        end do
      end do
    else

! Standard derivatives if no group defined:

      jx = nmrat(iat, i)

      do m = 1, 3
        frc(jx + m) = frc(jx + m) - dfr(jatm + m)
      end do

    end if
  end do

  return

end subroutine nmrcmf

!*******************************************************************************
!
! Subroutine:   nmrcms (NMR Center of MaSs)
!
! Description:
!
! This routine calculates the center of mass of the group(s) defined for
! distance restraint i, and places them in the xcom array.
!
! Author: David A. Pearlman
! Date: 8/89
!
! Note: This routine assumes that all atoms of the group are in the
!       same periodic "box" (if applicable).
!              
!*******************************************************************************

subroutine nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision  crd(*)
  double precision  xcom(6)
  integer           nmrat(4, 1)
  integer           nmrcom(2, 1)
  double precision  rmstot(2)
  integer           i

! Local variables:

  integer           iat
  integer           ifill
  integer           ip
  integer           j
  integer           jx
  double precision  rmass
  double precision  xtot, ytot, ztot

  do iat = 1, 2
    xtot = 0.0d0
    ytot = 0.0d0
    ztot = 0.0d0
    rmstot(iat) = 0.0d0
    ifill = 3*(iat - 1)

    if (nmrat(iat, i) .lt. 0) then
      do ip = -nmrat(iat, i), -nmrat(2 + iat, i)
        do j = nmrcom(1, ip), nmrcom(2, ip)

          jx = 3*(j - 1)
          rmass = atm_mass(j)

          xtot = xtot + crd(jx + 1)*rmass
          ytot = ytot + crd(jx + 2)*rmass
          ztot = ztot + crd(jx + 3)*rmass
          rmstot(iat) = rmstot(iat) + rmass
        end do
      end do
        xcom(ifill + 1) = xtot/rmstot(iat)
        xcom(ifill + 2) = ytot/rmstot(iat)
        xcom(ifill + 3) = ztot/rmstot(iat)
      else
        xcom(ifill + 1) = crd(nmrat(iat, i) + 1)
        xcom(ifill + 2) = crd(nmrat(iat, i) + 2)
        xcom(ifill + 3) = crd(nmrat(iat, i) + 3)
      end if
  end do

  return

end subroutine nmrcms

!*******************************************************************************
!
! Subroutine:   nmrgrp  (NMR GRouP)
!
! Description:
!
! This routine is called by nmrred if the user is defining groups of atoms
! whose center-of-mass is to be used in distance restraints. Here, the
! atom definitions are read in, sequentially sorted, and stored in the
! nmrcom/nmrat arrays.
!
! The groups are stored as ranges of atoms in array nmrcom(2,i):
!
! For the group taking the place of iat1 (kat = 0):
!
! (nmrcom(1,k) -> nmrcom(2,k)), k=-nmrat(1,i), -nmrat(3,i)
!
! For the group taking the place of iat2 (kat = 1):
!
! (nmrcom(1,k) -> nmrcom(2,k)), k=-nmrat(2,i), -nmrat(4,i)
!
! If iresid = 1: igr(i) is the residue # for atom i, and grnam(i) is the
!                corresponding atom name.
!           = 0: igr(i) is the atom number for atom i.
!
! name() and ipres() are the name and residue pointers for the system;
!     nres is the number of residues in the system; all are used by getnat 
!     when iresid = 1.
!
! maxigr is the dimension of igr.
! lastpk is the first unfilled location in nmrat
! natom is the number of atoms
!
! iscrth is a scratch array of >= natom length
!
! Author: David A. Pearlman
! Date: 8/89
!
!*******************************************************************************

subroutine nmrgrp(i, kat, lastpk, natom, nmrat, nmrcom, iout, maxgrp, igr, &
                  grnam, iresid, name, ipres, nres, maxigr)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer           i
  integer           kat
  integer           lastpk
  integer           natom
  integer           nmrat(4, 1)
  integer           nmrcom(2, 1)
  integer           iout
  integer           maxgrp
  integer           igr(*)
  character(4)      grnam(*)
  integer           iresid
  character(4)      name(*)
  integer           ipres(*)
  integer           nres
  integer           maxigr

! Local variables:

  integer, dimension(natom) :: iscrth(natom)

  integer           ierr
  integer           ii
  character(4)      igrnam
  integer           ipack
  integer           k

  iscrth = 0    ! Array assignment

! Loop over the defining atoms. If iresid = 1, then convert the
! residue #/atom-name reference into an actual atom number pointer, and store
! this pointer back in igr(i).

  do k = 1, maxigr
    if (igr(k) .le. 0) goto 30

! Resolve pointer, if necessary:

    if (iresid .eq. 1) then
      read(grnam(k), '(a4)') igrnam
      call getnat(igr(k), igrnam, name, ipres, nres, iout, ierr)
      if (ierr .eq. 1) then
        call mexit(6, 1)
      end if
    end if

    if (igr(k) .gt. natom) goto 1007
    iscrth(igr(k)) = 1
  end do

30 nmrat(1 + kat, i) = - lastpk

  ipack = 0
  do k = 1, natom
    if (iscrth(k) .eq. 1) then
      if (ipack .eq. 0) then
        nmrcom(1, lastpk) = k
        ipack = 1
      end if
    else if (iscrth(k) .eq. 0) then
      if (ipack .eq. 1) then
        nmrcom(2, lastpk) = k - 1
        ipack = 0
        lastpk = lastpk + 1
        if (lastpk .gt. maxgrp) goto 1009
      end if
    end if
  end do

  if (ipack .eq. 1) then
    nmrcom(2, lastpk) = natom
    ipack = 0
    lastpk = lastpk + 1
  end if

  if (lastpk .gt. maxgrp) goto 1009

  nmrat(3 + kat, i) = -(lastpk - 1)

  return

! A group-defining atom number was out-of-range:

1007 write(iout, 2003) k

     write(iout, 9001) (igr(ii), ii = 1, k)
     call mexit(6, 1)

! Too many atom ranges to be stored, relative to storage allocated (MAXGRP):

1009 write(iout, 2005) maxgrp
     call mexit(6, 1)

2003 format(' Error: Atom ', i2, ' in following group definition is', &
            ' greater than total # atoms') 

2005 format(' Error: Too many atom ranges need to be stored for ', &
            'center-of-mass distance', /, t9, 'restraints. ', &
            'MAXGRP =', i5, '. This needs to be increased.')

9001 format(16i5)

end subroutine nmrgrp

!*******************************************************************************
!
! Subroutine:   getnat  (GET Numbers of AToms)
!
! Description:
!
! This routine takes pointers to residue numbers (iat), and
! the name of the desired atom in this residue (iatnm), and
! returns the actual atom number pointers in iat.
!
! Author: David A. Pearlman
! Date: 3/89
!
!*******************************************************************************

subroutine getnat(iat, iatnm, name, ipres, nres, iout, ierr)

  implicit none

! Formal arguments:

  integer           iat
  character(4)      iatnm
  character(4)      name(*)
  integer           ipres(*)
  integer           nres
  integer           iout
  integer           ierr

! Local variables:

  integer           i

  if (iat .le. 0) return
  if (iat .gt. nres) goto 8000

  do i = ipres(iat), ipres(iat + 1) - 1
    if (iatnm .eq. name(i)) then
      iat = i
      return
    end if
  end do

! If we get here, the desired atom was not found. Report this to the user
! and return with ierr = 1:

  write(iout, 9000) iatnm, iat
  ierr = 1
  return

8000 write(iout, 9001) iatnm, iat
     ierr = 1
     return

! Error format statements:

9000 format(' Error: No atom ', a4, ' in residue ', i5)

9001 format(' Error: residue_number/atom_name reference specifies a', &
            /, '        residue number greater than the last residue in', &
            ' the system.', /, t9, 'Res = ', i6, ', Atom = "', a4, '"')

end subroutine getnat

!*******************************************************************************
!
! Subroutine:   nmrsht (NMR SHorT range interactions identification)
!
! Description:
!
! This routine determines which restraints should be classified as
! "short-range", based on the information given by the user on a SHORT
! command line. The nmrshb(i) array is set to indicate the short-range
! interactions:
!       nmrshb(i) = 0 : Restraint i is not short-range
!                 = 1 : Restraint i is short-range based on residue # criterion
!                 = 2 : Restraint i is short-range based on distance criterion
!                 = 3 : Restraint i is short range based on both residue # and
!                       distance criteria.
!
! Author: David A. Pearlman
! Date: 11/89
!
! INPUT:
!
! crd(i)       : Coordinate array.
! nmrnum       : The number of NMR restraints.
! rimass(i)    : The inverse mass for atom i.
! nstep        : The step/iteration number in the main calling program.
! nmrat(4,i)   : The 2-4 atoms defining restraint i.
! nmrcom(2,i)  : The ranges of atoms defining the center-of-mass group,
!                for applicable distance restraints
! nres         : The number of residues in the system.
! ishrtb       : If ishrtb >0, short-range interactions have been defined;
!                Atom-atom interactions are defined as short-range if either:
!
!                i) iwtstp(1,ishrtb) <= abs(residue_diff) <= iwtstp(2,ishrtb)
!
!                   (where residue_diff is the difference in the numbers of
!                    the residues containing atoms i and j)
!
!                ii) wtnrg(1,ishrtb) <= dist(i-j) <= wtnrg(2,ishrtb)
!
!                    (where dist(i-j) is the distance between atoms i and j).
!
!                A restraint is considered short-range if all the bonds
!                contributing to the restraint are short-range.
! 
!                The short-range distance classification is re-evaluated every
!                iwtstp(3,ishrtb) steps, if iwtstp(3,ishrtb) is > 0.
!
! Input variables in global memory:
!
! gbl_res_atms(i) : gbl_res_atms(i) is the first atom in residue i.
!                   gbl_res_atms(i+1)-1 is the last atom in residue i.
!
! iwtstp(3,i)
! wtnrg(2,i)  : See above.
! rk2nmr(2,i) : The slope and intercept for k2.
! rk3nmr(2,i) : The slope and intercept for k3.
! nmrst(3,i)  : If nmrst(3,i) > 0, k2 & k3 are linearly interpolated with
!               step number. If nmrst(3,i) < 0, then k2 & k3 are generated
!               by multiplication by a constant factor.
! wtls(2)     : wtls(1) is current weight of short-range interactions (in modwt)
!               wtls(2) is current weight of non "  "    interactions (in modwt)
! nmrfty(i)   : Flag for each restraint indicating what functional form should
!               be used for that restraint. If time-averaged restraints
!               are being applied to a particular class of internal and
!               nmrfty(i)=1, then an instantaneous restraint will be applied.
! bave(i)
! bave0(i)    : For time-averaged restraints. bave(i) contains the time-averaged
!               values of the distances, using exponential decay with
!               a time-constant of tauave(1). bave0(i) is the same thing,
!               but without exponential decay (real average).
! nave(3)     : nave(i) > 0 for time-averaged restraints.
! ipower(3)   : For time-averaged restraints. Gives the exponent used in
!               developing the time average. 
! tauave(3)   : For time-averaged restraints. Gives the time constant used
!               in the exponential decay multiplier.
! navint(3)   : For time-averaged restraints. Averaged values only updated
!               every navint(i) steps.
! ravein(3)   : For time-averaged restraints. On the first step, the value
!               returned is the current value+ravein(i) (i=1 for bonds).
! dt          : Integration timestep (ps). Used for time-averaged restraints.
! iout        : Unit for informational (error) prints.
! iflag       : = 0 First call to this routine. Determine short-range
!                   restraints meeting either criterion (i) or (ii) above. 
!                   No check against iwtstp(3,ishrtb) is done.
!               > 0 then only re-evaluate short-range restraints for
!                   criterion (ii) above. A check against iwtstp(3,ishrtb)
!                   is performed.
!
! OUTPUT:
!
! nmrshb(i) : Short-range indicator flag, set for each restraint as
!             described above.
!
! Notes: 1) When a distance restraint has been defined where one or
!           both ends of the restraints refer to a center-of-mass coordinate,
!           this routine assumes the atom ranges have been packed in nmrcom
!           in increasing sequential order, and that residue number increases
!           monotonically with atom number. Thus, only the residues
!           containing the first and last atoms of the center-of-mass definition
!           need to be examined for the residue check (i).
!
!         2) If time-averaged restraints are being monitored, the time
!            averaged values of distances will be used in determining
!            short-range interactions for bond restraints. HOWEVER, the
!            short-range status of angle and torsion restraints will be
!            determined using instantaneous values of the bond lengths
!            comprising these internals.
!
!*******************************************************************************

subroutine nmrsht(crd, nmrnum, nstep, nmrat, nmrcom, ishrtb, iwtstp, wtnrg, &
                  rk2nmr, rk3nmr, nmrst, nmrshb, wtls, nmrfty, igravt, &
                  ialtdis, bave, bave0, nave, ipower, tauave, navint, &
                  ravein, dt, iout, iflag)

  use prmtop_dat_mod

  implicit none

! Formal arguments:

  double precision  crd(*)
  integer           nmrnum
  integer           nstep
  integer           nmrat(4, 1)
  integer           nmrcom(2, 1)
  integer           ishrtb
  integer           iwtstp(3, 1)
  double precision  wtnrg(2, 1)
  double precision  rk2nmr(2, 1)
  double precision  rk3nmr(2, 1)
  integer           nmrst(3, 1)
  integer           nmrshb(1)
  double precision  wtls(2)
  integer           nmrfty(1)
  integer           igravt(1)
  integer           ialtdis(1)
  double precision  bave(*)
  double precision  bave0(*)
  integer           nave(3)
  integer           ipower(3)
  double precision  tauave(3)
  integer           navint(3)
  double precision  ravein(3)
  double precision  dt
  integer           iout
  integer           iflag

! Local variables:

  double precision  dum
  double precision  dum_array(6)    ! maximum explicit size for dummy
  integer           i
  integer           iat1, iat2
  integer           iat(4)
  integer           iave
  integer           ires1, ires2
  integer           irs1, irs2, irs3, irs4
  integer           irsdif
  integer           ishort
  integer           j
  double precision  rlow
  double precision  rmult1
  double precision  rmult2
  double precision  rmstot(2)
  double precision  rr, rr1, rr2, rr3
  double precision  rup
  double precision  small
  double precision  xcom(6)

  parameter (small = 1.0d-12)

  if (iflag .eq. 0) then
    do i = 1, nmrnum
      nmrshb(i) = 0
    end do
  end if

! Return if no short-range interactions defined:

  if (ishrtb .eq. 0) return

! If iflag .ne. 0, re-determine short-range interactions only every
! iwtstp(3,ishrtb) steps:

  if (iflag .ne. 0) then
    if (iwtstp(3, ishrtb) .eq. 0) return
    if (mod(nstep, iwtstp(3, ishrtb)) .ne. 0 .or. nstep .eq. 0) return
  end if

! Check for residue range defined short-range interactions if iflag=0:

  if (iflag .eq. 0) then
    do i = 1, nmrnum
      if (nmrat(3, i) .lt. 0) then

! Center of mass-type distance restraints:

        if (nmrat(1, i) .lt. 0 .or. nmrat(2, i) .lt. 0) then
          do j = 1, 2
            if (nmrat(j, i) .lt. 0) then
              iat(2*(j - 1) + 1) = nmrcom(1, - nmrat(j, i))
              iat(2*(j - 1) + 2) = nmrcom(2, - nmrat(j + 2, i))
            else
              iat(2*(j - 1) + 1) = nmrat(j, i)
              iat(2*(j - 1) + 2) = nmrat(j, i)
            end if
          end do
          iat1 = min(iat(1), iat(3))
          iat2 = max(iat(2), iat(4))
          call at2res(iat1, gbl_res_atms, nres, ires1, 1, iout)
          call at2res(iat2, gbl_res_atms, nres, ires2, 1, iout)
        else

! Regular distance restraints:

          call at2res(nmrat(1, i), gbl_res_atms, nres, ires1, 1, iout)
          call at2res(nmrat(2, i), gbl_res_atms, nres, ires2, 1, iout)
        end if

! Angle restraints:

      else if (nmrat(4, i) .lt. 0) then
        call at2res(nmrat(1, i), gbl_res_atms, nres, irs1, 1, iout)
        call at2res(nmrat(2, i), gbl_res_atms, nres, irs2, 1, iout)
        call at2res(nmrat(3, i), gbl_res_atms, nres, irs3, 1, iout)
        ires1 = min(irs1, irs2)
        ires1 = min(ires1, irs3)
        ires2 = max(irs1, irs2)
        ires2 = max(ires2, irs3)

! Torsional restraints:

      else
        call at2res(nmrat(1, i), gbl_res_atms, nres, irs1, 1, iout)
        call at2res(nmrat(2, i), gbl_res_atms, nres, irs2, 1, iout)
        call at2res(nmrat(3, i), gbl_res_atms, nres, irs3, 1, iout)
        call at2res(nmrat(4, i), gbl_res_atms, nres, irs4, 1, iout)
        ires1 = min(irs1, irs2)
        ires1 = min(ires1, irs3)
        ires1 = min(ires1, irs4)
        ires2 = max(irs1, irs2)
        ires2 = max(ires2, irs3)
        ires2 = max(ires2, irs4)
      end if

! Do the residue check:

      irsdif = abs(ires2 - ires1)
      if (iwtstp(1, ishrtb) .le. irsdif .and. &
          irsdif .le. iwtstp(2, ishrtb)) nmrshb(i) = 1
    end do
  end if

! Now do the distance check:

! Define conversion factors required to account for weight changes
! if short-range classification of a restraint changes:

  if (wtls(2) .lt. small) then
    rmult1 = 0.0d0
  else
    rmult1 = wtls(1)/wtls(2)
  end if

  if (wtls(1) .lt. small) then
    rmult2 = 0.0d0
  else
    rmult2 = wtls(2)/wtls(1)
  end if

  rlow = wtnrg(1, ishrtb)
  rup = wtnrg(2, ishrtb)
  do i = 1, nmrnum
    ishort = 0
    iave = 0
    if (nmrat(3, i) .lt. 0) then

! Center of mass-type distance restraints:

      if (nave(1) .gt. 0 .and. nmrfty(i) .eq. 0) iave = 1

      if ((nmrat(1, i) .lt. 0 .or. nmrat(2, i) .lt. 0) .and. &
          igravt(i) .eq. 0) then

        call nmrcms(crd, xcom, nmrat, nmrcom, rmstot, i)

        call disnrg(xcom, dum_array, dum_array, rr, 0, 3, dum, dum, dum, &
                    dum, dum, dum, dum, bave(2*i - 1), bave0(i), ipower(1), &
                    tauave(1), ravein(1), dt, navint(1), iave, 0, 0, 1, &
                    ialtdis(i))

        if (rlow .le. rr .and. rr .le. rup) ishort = 1

! r**-6 averaged distance interaction restraints:

      else if ((nmrat(1, i) .lt. 0 .or. nmrat(2, i) .lt. 0) .and. &
                igravt(i) .eq. 1) then
        call r6ave(crd, nmrat, nmrcom, rr, i) 

      else

! Regular distance restraints:

        call disnrg(crd, dum_array, dum_array, rr, nmrat(1, i), nmrat(2, i), &
                    dum, dum, dum, dum, dum, dum, dum, bave(2*i - 1), &
                    bave0(i), ipower(1), tauave(1), ravein(1), dt, &
                    navint(1), iave, 0, 0, 1, ialtdis(i))

        if (rlow .le. rr .and. rr .le. rup) ishort = 1
      end if

! Angle restraints:

    else if (nmrat(4, i) .lt. 0) then
      iave = 0

      call disnrg(crd, dum_array, dum_array, rr1, nmrat(1, i), nmrat(2, i), &
                  dum, dum, dum, dum, dum, dum, dum, dum_array, dum_array, &
                  ipower(2), tauave(2), ravein(2), dt, 0, 0, 0, 0, 1, &
                  ialtdis(i))

      call disnrg(crd, dum_array, dum_array, rr2, nmrat(2, i), nmrat(3, i), &
                  dum, dum, dum, dum, dum, dum, dum, dum_array, dum_array, &
                  ipower(2), tauave(2), ravein(2), dt, 0, 0, 0, 0, 1, &
                  ialtdis(i))

      if ((rlow .le. rr1 .and. rr1 .le. rup) .and. &
           (rlow .le. rr2 .and. rr2 .le. rup)) ishort = 1

! Torsional restraints:

    else
      iave = 0

      call disnrg(crd, dum_array, dum_array, rr1, nmrat(1, i), nmrat(2, i), &
                  dum, dum, dum, dum, dum, dum, dum, dum_array, dum_array, &
                  ipower(3), tauave(3), ravein(3), dt, 0, iave, 0, 0, 1, &
                  ialtdis(i))

      call disnrg(crd, dum_array, dum_array, rr2, nmrat(2, i), nmrat(3, i), &
                  dum, dum, dum, dum, dum, dum, dum, dum_array, dum_array, &
                  ipower(3), tauave(3), ravein(3), dt, 0, iave, 0, 0, 1, &
                  ialtdis(i))

      call disnrg(crd, dum_array, dum_array, rr3, nmrat(3, i), nmrat(4, i), &
                  dum, dum, dum, dum, dum, dum, dum, dum_array, dum_array, &
                  ipower(3), tauave(3), ravein(3), dt, 0, iave, 0, 0, 1, &
                  ialtdis(i))

      if ((rlow .le. rr1 .and. rr1 .le. rup) .and. &
        (rlow .le. rr2 .and. rr2 .le. rup) .and. &
        (rlow .le. rr3 .and. rr3 .le. rup)) ishort = 1

    end if

! If the short-range classification of the restraint has changed, change
! nmrshb(i) to reflect this fact. Also, modify k2 and k3 so that they
! are multiplied by the weight factor appropriate for their new classification:

    if (nmrshb(i) .le. 0 .and. ishort .eq. 1) then
      nmrshb(i) = 2
      rk2nmr(2, i) = rk2nmr(2, i)*rmult1
      rk3nmr(2, i) = rk3nmr(2, i)*rmult1
      if (nmrst(3, i) .gt. 0) then
        rk2nmr(1, i) = rk2nmr(1, i)*rmult1
        rk3nmr(1, i) = rk3nmr(1, i)*rmult1
      end if
    else if (nmrshb(i) .eq. 2 .and. ishort .eq. 0) then
      nmrshb(i) = 0
      rk2nmr(2, i) = rk2nmr(2, i)*rmult2
      rk3nmr(2, i) = rk3nmr(2, i)*rmult2
      if (nmrst(3, i) .gt. 0) then
        rk2nmr(1, i) = rk2nmr(1, i)*rmult2
        rk3nmr(1, i) = rk3nmr(1, i)*rmult2
      end if
    else if (nmrshb(i) .eq. 1 .and. ishort .eq. 1) then
      nmrshb(i) = 3
    else if (nmrshb(i) .eq. 3 .and. ishort .eq. 0) then
      nmrshb(i) = 1
    end if
  end do

  return

contains

!*******************************************************************************
!
! Internal Subroutine:  at2res (ATom 2 RESidue number)
!
! Description:
!
! This routine takes the atom number iat, and returns the number ires of
! the residue containing it.
!
! Author: David A. Pearlman
! Date: 11/89
!
! iat     : Input atom number
! ipres(i): ipres(i) is the first atom of residue i. ipres(i+1)-1 is the
!           last atom of residue i.
! nres    : The number of residues in the system
! ires    : The residue containing atom iat (output).
! iflag   : =0 iat is absolute atom number
!           >1 iat is (3*atom_number)-1
! iout    : Unit for informational prints
!
!*******************************************************************************

subroutine at2res(iat, ipres, nres, ires, iflag, iout)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  integer           iat
  integer           ipres(*)
  integer           nres
  integer           ires
  integer           iflag
  integer           iout

! Local variables:

  integer           i
  integer           icheck

  icheck = iat
  if (iflag .gt. 0) icheck = (iat / 3) + 1
  do i = 1, nres
    if (icheck .ge. ipres(i) .and. icheck .lt. ipres(i + 1)) then
      ires = i
      return
    end if
  end do

  write(iout, 9001) icheck, ipres(nres + 1) - 1
  call mexit(iout, 1)

9001 format(' Error from routine AT2RES: Atom number = ', i6, &
            ' cannot be translated', / , t29, 'to an appropriate ', &
            'residue number. Last atom in', / , t29, 'the system is ', i6)

end subroutine at2res

end subroutine nmrsht

!*******************************************************************************
!
! Subroutine:   opnmrg (OPen NMR General)
!
! Description:
!
! This subroutine does a general-format open of the file toopen
! on the unit specified as iunit. 
!
! istat = 0 : Open a new file
! istat = 1 : Open an existing file and position at the end of the file.
! istat = 2 : Open an existing file and position at the beginning of file
! istat = 3 : Close the file.
!
! iout : Unit for informational error message prints
!
! ierror: = 0 if file opened successfully
!         = 1 if there was an error opening the file.
!
! Author: David A. Pearlman
! Date: 3/90
!              
!*******************************************************************************

subroutine opnmrg(toopen, iunit, istat, iout, ierror)

  implicit none

! Formal arguments:

  character(*)      toopen
  integer           iunit
  integer           istat
  integer           iout
  integer           ierror

! Local variables:

  ierror = 0

  if (istat .eq. 0) then
    open(unit = iunit, file = toopen, status = 'new', form = 'formatted', &
         err = 1001)
  else if (istat .eq. 1) then
    open(unit = iunit, file = toopen, status = 'old', form = 'formatted', &
         err = 1002)

10 read(iunit, 9000, end = 20, err = 20)
   goto 10

20 continue

  else if (istat .eq. 2) then
    open(unit = iunit, file = toopen, status = 'old', form = 'formatted', &
         err = 1002)
    rewind(iunit)
  else if (istat .eq. 3) then
    close(iunit)
  end if

  return

! Process errors on open:

1001 write(iout, 9001) toopen
     ierror = 1
     return

1002 write(iout, 9002) toopen
     ierror = 1
     return

9000 format(a)
9001 format(' Warning: Error opening "New" file from subroutine ', &
            'OPNMRG', /, ' File = ', a)
9002 format(' Warning: Error opening "Old" file from subroutine ', &
            'OPNMRG', /, ' File = ', a)

end subroutine opnmrg

!*******************************************************************************
!
! Subroutine:   r6ave (R**-6 AVErage)
!
! Description: 
!
! This routine calculate the <r**-6>-6 average distance between two
! sets of atoms.
! 
! The resulting averaged distance is placed in rave.
! A number of other quantities, required in calculating the derivatives,
! are calculated and stored in common block r6avec. This common block
! only appears in this routine and in r6drv (where the derivatives are
! determined).
!
! Note: This routine assumes that all atoms of the group are in the
!       same periodic "box" (if applicable).
!
!       The maximum allowable number of pairs of atomic interactions used
!       in calculating the r**-6 average is currently hardwired in this routine
!       (maxpr). If you change it here, you must also change it in r6drv.
!              
!*******************************************************************************

subroutine r6ave(crd, nmrat, nmrcom, rave, indx)

  use pmemd_lib_mod

  implicit none

! Formal arguments:

  double precision  crd(*)
  integer           nmrat(4, 1)
  integer           nmrcom(2, 1)
  double precision  rave
  integer           indx

! Local variables:

! BUGBUG - r6avec should be in a header!

! r6avec common block. Shared between this routine and r6drv:

  integer           maxpr

  parameter (maxpr = 5000)

  double precision  xij(3, maxpr)
  integer           iat0(maxpr)
  integer           jat0(maxpr)
  double precision  rm8(maxpr)
  double precision  fsum
  double precision  rm6bar
  double precision  rij
  integer           nsum

  common /r6avec/ xij, iat0, jat0, rm8, fsum, rm6bar, rij, nsum

  integer           i
  integer           i3
  integer           iat1, iat2, iat3, iat4
  integer           ip
  integer           ipr
  integer           j
  integer           j3
  integer           jp
  double precision  one6
  double precision  rij2
  double precision  rm2

  parameter (one6 = 1.0d0/6.0d0)

  iat1 = nmrat(1, indx)
  iat2 = nmrat(2, indx)
  iat3 = nmrat(3, indx)
  iat4 = nmrat(4, indx)

! Case where iat1 > 0 (single atom) and iat2 < 0 (r**-6 average pos.):

  ipr = 0
  rm6bar = 0.0d0
  if (iat1 .ge. 0) then
    i3 = iat1
    do jp = - iat2, - iat4
      do j = nmrcom(1, jp), nmrcom(2, jp)
        j3 = 3*(j - 1)
        ipr = ipr + 1
        if (ipr .gt. maxpr) goto 500
        xij(1, ipr) = crd(i3 + 1) - crd(j3 + 1)
        xij(2, ipr) = crd(i3 + 2) - crd(j3 + 2)
        xij(3, ipr) = crd(i3 + 3) - crd(j3 + 3)
        iat0(ipr) = i3
        jat0(ipr) = j3
        rij2 = xij(1, ipr)**2 + xij(2, ipr)**2 + xij(3, ipr)**2
        rm2 = 1.0/rij2
        rm6bar = rm6bar + rm2**3
        rm8(ipr) = rm2**4
      end do
    end do
    nsum = ipr
    fsum = dble(nsum)
    rm6bar = rm6bar/fsum
    rij = rm6bar**(- one6)

! Case where iat1 < 0 (r**-6 ave. pos.)  and iat2 > 0:

  else if (iat2 .ge. 0) then

    j3 = iat2
    do ip = - iat1, - iat3
      do i = nmrcom(1, ip), nmrcom(2, ip)
        i3 = 3*(i - 1)
        ipr = ipr + 1
        if (ipr .gt. maxpr) goto 500
        xij(1, ipr) = crd(i3 + 1) - crd(j3 + 1)
        xij(2, ipr) = crd(i3 + 2) - crd(j3 + 2)
        xij(3, ipr) = crd(i3 + 3) - crd(j3 + 3)
        iat0(ipr) = i3
        jat0(ipr) = j3
        rij2 = xij(1, ipr)**2 + xij(2, ipr)**2 + xij(3, ipr)**2
        rm2 = 1.0/rij2
        rm6bar = rm6bar + rm2**3
        rm8(ipr) = rm2**4
      end do
    end do

    nsum = ipr
    fsum = dble(nsum)
    rm6bar = rm6bar/fsum
    rij = rm6bar**(- one6)

! Case where iat1 < 0 and iat2 < 0 (both r**-6 ave. pos.):

  else

    do jp = - iat2, - iat4
      do j = nmrcom(1, jp), nmrcom(2, jp)
        j3 = 3*(j - 1)
        do ip = - iat1, - iat3
          do i = nmrcom(1, ip), nmrcom(2, ip)
            i3 = 3*(i - 1)
            ipr = ipr + 1
            if (ipr .gt. maxpr) goto 500
            xij(1, ipr) = crd(i3 + 1) - crd(j3 + 1)
            xij(2, ipr) = crd(i3 + 2) - crd(j3 + 2)
            xij(3, ipr) = crd(i3 + 3) - crd(j3 + 3)
            iat0(ipr) = i3
            jat0(ipr) = j3
            rij2 = xij(1, ipr)**2 + xij(2, ipr)**2 + xij(3, ipr)**2
            rm2 = 1.0/rij2
            rm6bar = rm6bar + rm2**3
            rm8(ipr) = rm2**4
          end do
        end do
      end do
    end do
    nsum = ipr
    fsum = dble(nsum)
    rm6bar = rm6bar/fsum
    rij = rm6bar**(- one6)

  end if

  rave = rij

  return

! Errors:

500 write(mdout, 1000) maxpr
    call mexit(6, 1)

! Format statements:

1000 format(' ERROR: Number of interaction distances in r**-6 averaged ', &
            'distance', /, t8, ' exceeds maximum allowed by MAXPR = ', i5)

end subroutine r6ave

!*******************************************************************************
!
! Subroutine:   r6drv   (R**-6 DeRiVatives)
!
! Description:
!
! This routine calcuates the derivatives associated with an r**-6
! average position defined group, and adds them into the force (f)
! array. 
!
! The derivatives are calculated as
!
!           dE/dr6 * dr6/dx
!
! where dE/dr6 is the derivative of the energy with respect to the
! averaged distance (calculated in routine disnrg, and passed here as df)
! and dr6/dx is the derivative of the averaged distance with respect
! to each individual coordinate used in the distance.
!
! It is this latter distance which is calculated here. Note that
! most of the quantities required for calculating the derivatives are
! passed through common block r6avec. this block appears only here
! and in routine r6ave.
!
! maxpr is the maximum number of distances used in computing <r**-6>**-1/6. 
! Must be changed here and in r6ave, if desired.
!              
!*******************************************************************************

subroutine r6drv(frc, df)

  implicit none

! Formal arguments:

  double precision  frc(*)
  double precision  df

! Local variables:

! BUGBUG - r6avec should be in a header!

! r6avec common block. Shared between this routine and r6drv:

  integer           maxpr

  parameter (maxpr = 5000)

  double precision  xij(3, maxpr)
  integer           iat0(maxpr)
  integer           jat0(maxpr)
  double precision  rm8(maxpr)
  double precision  fsum
  double precision  rm6bar
  double precision  rij
  integer           nsum

  common /r6avec/ xij, iat0, jat0, rm8, fsum, rm6bar, rij, nsum

  double precision  fact
  double precision  fact2
  integer           i3
  integer           ipr
  integer           j3
  double precision  xh, yh, zh

  fact = df*rij/(fsum*rm6bar)
  do ipr = 1, nsum
    i3 = iat0(ipr)
    j3 = jat0(ipr)
    fact2 = fact*rm8(ipr)
    xh = fact2*xij(1, ipr)
    yh = fact2*xij(2, ipr)
    zh = fact2*xij(3, ipr)
    frc(i3 + 1) = frc(i3 + 1) - xh
    frc(i3 + 2) = frc(i3 + 2) - yh
    frc(i3 + 3) = frc(i3 + 3) - zh
    frc(j3 + 1) = frc(j3 + 1) + xh
    frc(j3 + 2) = frc(j3 + 2) + yh
    frc(j3 + 3) = frc(j3 + 3) + zh
  end do

  return

end subroutine r6drv

!*******************************************************************************
!
! Subroutine:  setgms (SET GaMS)
!
! Description: <TBS>
!              
! This routine sets the values of the gamc() and gams() arrays, which
! are used in vectorized torsional energy routines. This routine is only
! called when the torsional force constants are changed in routine modwt.
! Otherwise, these arrays are set only once, by a call to dihpar from rdparm,
! at the start of the program. This is an abbreviated version of dihpar which
! passes most arguments by common, not call-list.
!
! Author: David A. Pearlman
! Date: 5/92
!
!*******************************************************************************

subroutine setgms(numphi)

  use gbl_constants_mod
  use prmtop_dat_mod

  implicit none

! Formal arguments:

  integer           numphi

! Local variables:

  double precision  dum
  double precision  dumc
  double precision  dums
  double precision  four
  integer           i
  double precision  one
  double precision  pim
  double precision  tenm3
  double precision  tenm10
  double precision  zero

  data zero, one, tenm3, tenm10/0.0d+00, 1.0d+00, 1.0d-03, 1.0d-06/
  data four /4.0d+00/

  pim = four*atan(one)
  do i = 1, numphi
    dum = gbl_phase(i)
    if (abs(dum - PI) .le. tenm3) dum = sign(pim, dum)
    dumc = cos(dum)
    dums = sin(dum)
    if (abs(dumc) .le. tenm10) dumc = zero
    if (abs(dums) .le. tenm10) dums = zero
    gbl_gamc(i) = dumc * gbl_pk(i)
    gbl_gams(i) = dums * gbl_pk(i)
  end do

  return

end subroutine setgms

end module nmr_lib_mod
