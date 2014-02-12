#include "dprec.fh"

!-------------BEGIN    nmr.h  ------------------------------------------------

!  ---Header file for the chemical shifts and NOESY intensity

!     Because of the complexity of the storage requirements for
!     these calculations, (and because one of the authors is very
!     lazy,) storage allocation for these is not done in the
!     LOCMEM routine, but rather arranged at compile time through
!     the information given below.

!     If you do not plan to make use of this section of the code,
!     set the parameters below to small numbers to avoid allocating
!     space unnecessarily.  When/if you change your mind, reset the
!     parameters and re-compile.

!     Input parameters (things you have to set):

!     MATOM = max # of atoms in the system
!     MXR = max # residues
!     MA = max # of protons in any sub-molecule
!     MXTAU = max number of mixing times
!     MXP = max number of input intensities (peaks) per mixing time
!     MTOT = max. # total peaks, all mixing times, all sub-molecules
!     MXVAR = max. # of "extra" dynamic variables
!     MRING = max. # of rings for chemical shift calculation
!     MSHF = max # of protons whose shifts are to be calculated
!     MAXDIP = max # of residueal dipolar couplings
!     MAXCSA = max # of residual csa measurements

integer mring,mshf,mxvar,matom,ma,ma2,lst,mxr,mxtau,mxp, &
      mtot,maxdip,maxdipsets,maxcsa

!  --- "standard" parameters for jobs like the test cases:

parameter (mring=50)
parameter (mshf=500)
parameter (mxvar=50)
parameter (matom=50000)
parameter (ma=100)
parameter (ma2=ma*ma)
parameter (lst=(ma2+ma)/2)
parameter (mxr=300)
parameter (mxtau=5)
parameter (mxp=100)
parameter (mtot=500)
parameter (maxdip=2000)
parameter (maxcsa=200)
parameter (maxdipsets=2)

integer isubi,isubr
parameter (isubi=8 + 3*ma + 2*matom + mxtau + 2*mxtau*mxp + mxp)
parameter (isubr=3*ma + mxtau + 3*mxtau*mxp + 4)

integer peakid(mxp)
_REAL_ tau(ma),pop(ma),popn(ma),emix(mxtau), &
      aexp(mxtau,mxp),awt(mxtau,mxp),arange(mxtau,mxp),oscale, &
      omega,taumet,taurot,invwt1,invwt2
integer nath,natmet,nummt,id2o,iroesy,ihet,nvect,ihsful,m2(ma), &
      inn(ma),ihyp(ma),ihyd(matom),inatom(matom),npeak(mxtau), &
      ihp(mxtau,mxp),jhp(mxtau,mxp)
common /methylr/ tau,pop,popn,emix,aexp,awt,arange,oscale, &
      omega,taumet,taurot,invwt1,invwt2
common /methyli/ nath,natmet,nummt,id2o,iroesy,ihet,nvect,ihsful,m2, &
      inn,ihyp,ihyd,inatom,npeak,ihp,jhp,peakid

!    Parameters for parallel broadcast

integer BC_METHYLR,BC_METHYLI,BC_ALIGNR,BC_ALIGNI
parameter(BC_METHYLR=3*ma+mxtau+3*mxtau*mxp+6)
parameter(BC_METHYLI=8 + 3*ma + 2*matom + mxtau + 2*mxtau*mxp)
parameter(BC_ALIGNR=5*maxdip + 1 + 6*maxdipsets)
parameter(BC_ALIGNI=3*maxdip + 4)

integer nmropt,iprint,noeskp,iscale,ipnlty,iuse,maxsub,jar,morse
_REAL_ scalm,pencut,ensave,tausw,ebdev,eadev,drjar
common/nmr1/scalm,pencut,ensave,tausw,ebdev,eadev,drjar, &
      nmropt,iprint,noeskp,iscale,ipnlty,iuse,maxsub,jar,morse

character(len=14) resat(matom)
common/nmr2/ resat

!   residual dipolar coupling (aka "alignment") restraint information:

_REAL_ dobsu(maxdip),dobsl(maxdip),dcut,gigj(maxdip), &
      dij(maxdip),dwt(maxdip), &
      s11(maxdipsets),s12(maxdipsets),s13(maxdipsets), &
      s22(maxdipsets),s23(maxdipsets),s33(maxdipsets)
integer ndip,num_datasets,id(maxdip),jd(maxdip),dataset(maxdip),ifreeze, &
        ifreezes
common/align/dobsu,dobsl,dcut,gigj,dij,dwt,s11,s12,s13,s22,s23,s33, &
      ndip,id,jd,dataset,num_datasets,ifreeze,ifreezes

!   residual csa shift restraint information:

_REAL_ cobsu(maxcsa),cobsl(maxcsa),ccut,cwt(maxcsa), &
       sigma(3,3,maxcsa),field(maxcsa)
integer ncsa,icsa(maxcsa),jcsa(maxcsa),kcsa(maxcsa),datasetc(maxcsa)
common/csa/cobsu,cobsl,ccut,cwt,sigma,field,ncsa,icsa,jcsa,kcsa,datasetc

#ifdef NMODE

   !   --- special variables only used for (unsupported-as-of-now)
   !         normal-mode/NMR code:

   integer mxvect,nmsnap
   logical bose,per_mode
   parameter (mxvect=70)
   _REAL_ vect(3*matom,mxvect),freq(mxvect),xdev,omegax,vtemp
   common/modes/ vect,freq,xdev,omegax,vtemp,bose,per_mode,nmsnap

   _REAL_ gamma_nmr, dgamma(3+mxvect),dratg(ma,ma,mxvect)
   integer iusev
   common/derivg/ gamma_nmr, dgamma,dratg,iusev

#else

   integer mxvect
   parameter (mxvect=1)

#endif

! Common block containing variables relating to nmr restraints.

integer       intreq,irlreq,lnmr01,inmr02,iprr,iprw
common/nmrstf/intreq,irlreq,lnmr01,inmr02,iprr,iprw

integer ntot,ntota,ipmix(mtot),ntotb
_REAL_ calc(mtot),exper(mtot),calca(mtot),expera(mtot),calcb(mtot),experb(mtot)
common/correl/ntot,ntota,calc,exper,calca,expera,calcb,experb,ipmix,ntotb

_REAL_        wnoesy,wshift,enoe,eshf,epcshf,ealign,ecsa
common/wremar/wnoesy,wshift,enoe,eshf,epcshf,ealign,ecsa

!-------------END      nmr.h  ------------------------------------------------

