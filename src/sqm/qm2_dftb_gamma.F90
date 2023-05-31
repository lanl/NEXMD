! <compile=optimized>

!=============================================================================

! gamma resulting from exact Coulomb interaction of normalized
! exp(-a*r) charge distribution
!
! Attention: subroutine gamsub needed
!
! input:  r:     distance
!         uhub1: Hubbard parameter orbital 1
!         uhub2: Hubbard parameter orbital 2
!
! output: gval:  gamma12 function value
!=============================================================================
#include "dprec.fh"

subroutine gam12(r,uhub1,uhub2,gval)
  IMPLICIT NONE
  _REAL_, parameter :: zero=1.0d-4
  _REAL_ :: gval,a1,a2,src,avg,uhub1,uhub2,rrc,rrc3
  _REAL_ :: val12,val21,drv12,drv21,r,fac,fac2,efac
  ! neu fuer H-bonds
  _REAL_ :: fhbond,uhubh,cut1,cut2,rcuth
  _REAL_ :: kl1,kl2,kl3,kkk


  ! This will define an adjustment in case one (or both)
  ! of the atoms is a Hydrogen. 
  !
  !             WHERE IS IT DESCRIBED??
  !
  ! Hubbard parameter for H is 0.4195 from hh.spl (SK file)
  
  uhubh= 0.4195007d0
  fhbond = 1.0d0


  ! This will NEVER be true, because the Hubbard parameter
  ! for H is defined as 0.4195 in the parameter file.
  ! If we change uhubh as above, so that it can be true,
  ! the results are very strange.
  if ((uhub1 == uhubh) .OR. (uhub2 == uhubh)) then
     kl1 = 4.0d0
     fhbond = exp(-(((uhub1 + uhub2) / 2)**kl1)*r**2)
  endif


  ! gamma  besteht aus einem 1/r Term, und etwas,
  ! was fuer r=0 gegen Hubbard geht:
  ! multipliziere einfach den zweiten Term mit fhbond!

  gval= 0.0

  ! a1 = \tau_1 = 16/5 U_1 = 3.2 U_1
  a1= 3.2*uhub1
  a2= 3.2*uhub2

  IF (a1+a2 < zero) THEN
     RETURN
  ENDIF

  ! Takes an average of the 2 \taus. Why?
  ! this is only used in case the distance
  ! between the atoms is smaller than 1e-4, in 
  ! which case they are virtually the same atom, 
  ! OR if the hubbard parameters are within 
  ! 1e-5 of each other, in which case they 
  ! would be the same atom, but separated by
  ! some distance r (like two 'C' in the same 
  ! molecule.)

  ! If the 2 Hubbard parameters are the same,
  !         avg = \tau
  ! So that this probably is here to take care
  ! of small differences in the Hubbard pars.
 
  ! I CAN'T FIND THAT DESCRIBED ANYWHERE!

  src= 1.0/(a1+a2)
  fac= a1*a2*src
  avg= 1.6*(fac+fac*fac*src)

  IF (r < zero) THEN
     ! 'zero' means 1e-4. So, this is the r-->0 limit,
     ! and the atoms should be the same, and 
     !      gamma = 5/16 \tau = 0.3125 * avg
     gval= 0.3125*avg

  ELSE ! Atoms are not in the same position
     rrc= 1.0/r

     IF (abs(a1-a2) < 1.0d-5) THEN
        ! The atoms are (virtually) the same (e.g. two carbons).
        ! in this case, the equations are simplified.
        fac= avg*r
        fac2= fac*fac
        efac= exp(-fac)/48.0
        gval= (1.0-fhbond*(48.0+33.0*fac+fac2*(9.0+fac))*efac)*rrc
     ELSE
        ! Case of 2 different atoms
        call gamsub(a1,a2,r,rrc,val12,drv12)
        call gamsub(a2,a1,r,rrc,val21,drv21)
        gval= rrc-fhbond*val12-fhbond*val21
     ENDIF
  ENDIF
  RETURN
end subroutine gam12



!=============================================================================

! derivative of gamma resulting from exact Coulomb interaction
! of normalized  exp(-a*r) charge distribution
! Attention: subroutine gamsub needed

! input:  r:     distance
! uhub1: Hubbard parameter orbital 1
! uhub2: Hubbard parameter orbital 2
! output: gdrv:  ((d gamma12)/(d r)) * 1/r

!=============================================================================

subroutine gam121(r,uhub1,uhub2,gdrv)
  IMPLICIT NONE
  _REAL_ :: zero
  parameter(zero=1.0d-4)
  _REAL_ :: gdrv,a1,a2,src,avg,uhub1,uhub2,rrc,rrc3
  _REAL_ :: val12,val21,drv12,drv21,r,fac,fac2,efac
  ! neu fuer H-bonds

  _REAL_ :: fhbond,uhubh,cut1,cut2,rcuth,cutab(5,5)
  _REAL_ :: dhbond
  _REAL_ :: kl1,kl2,kl3,kkk
  kl1=4.0
  uhubh = 0.4195007d0
  fhbond = 1.0
  dhbond=0.0
  if((uhub1 == uhubh) .OR. (uhub2 == uhubh)) then
     fhbond= exp(-(((uhub1+uhub2)/2)**kl1)*r**2)
     dhbond= -2*fhbond*(((uhub1+uhub2)/2)**kl1)/r
  endif

  gdrv= 0.0
  a1= 3.2*uhub1
  a2= 3.2*uhub2
  IF (a1+a2 < zero) THEN
     RETURN
  ENDIF
  src= 1.0/(a1+a2)
  fac= a1*a2*src
  avg= 1.6*(fac+fac*fac*src)
  IF (r < zero) THEN
  ELSE
     rrc= 1.0/r
     rrc3= rrc*rrc*rrc
     IF (abs(a1-a2) < 1.0d-5) THEN
        fac= avg*r
        fac2= fac*fac
        efac= exp(-fac)/48.0
        gdrv= &
             -(1.0-fhbond*(48.0+48*fac+fac2*(24.0+7*fac+fac2))*efac)*rrc3 &
             + dhbond*(48.0+33*fac+fac2*(9.0+fac))*efac*rrc
     ELSE
        call gamsub(a1,a2,r,rrc,val12,drv12)
        call gamsub(a2,a1,r,rrc,val21,drv21)
        gdrv= -rrc3-fhbond*(drv12+drv21)*rrc &
             -dhbond*(val12+val21)
     ENDIF
  ENDIF
  RETURN
end subroutine gam121



!===========================================================================
! auxiliary routine needed by gam12 and gam121
! input   a:    alpha1
!         b:    alpha2
!         r:    distance
!         rrc:  1/distance
! output: gval: function value
!         gdrv: function derivative
!===========================================================================

subroutine gamsub(a,b,r,rrc,gval,gdrv)
  IMPLICIT NONE
  _REAL_ :: a,a2,b,b2,b4,b6,drc,drc2,r,efac,rrc,fac,gval,gdrv
  a2= a*a
  b2= b*b
  b4= b2*b2
  b6= b4*b2
  drc= 1.0/(a2-b2)
  drc2=drc*drc
  efac= exp(-a*r)
  fac= (b6-3*a2*b4)*drc2*drc*rrc
  gval= efac*(0.5*a*b4*drc2-fac)
  gdrv= -a*gval+efac*fac*rrc
  RETURN
end subroutine gamsub
