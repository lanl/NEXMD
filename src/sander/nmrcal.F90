#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmrcal here]
subroutine nmrcal(x,f,name,irsnam,ipres,rimass,enmr,devdis,devang, &
      devtor,devplpt,devpln,devgendis,temp0,tautp,cut,ntb,work,iwork,iscrth, &
      in,iout,rk,tk,pk,cn1,cn2,ag,bg,cg,numbnd, &
      numang,numphi,nimprp,nhb,ncharg,natom, &
      ntypes,nres,rad,wel,radhb,welhb,rwell,isftrp, &
      tgtrmsd,temp0les,nprnt,caltyp)

   use nmr
#ifdef MPI
   use file_io_dat
#endif

   ! This is the top-level calling routine for the NMR/MD refinement
   ! suite of programs.
   
   ! This suite was written by
   
   !            David A. Pearlman,
   !            Department of Pharmaceutical Chemistry
   !            University of California, San Francisco, 94143-0446
   
   ! Subroutine NMR CALls. This routine makes the calls to the following
   ! NMR-md refinement related routines:
   !       RESTAL : Does a cursory read of the NMR restraint files to
   !                determine how much memory will be required. This
   !                can be used in dynamic memory allocation schemes.
   !       NMRRED : Reads information defining NMR restraints & changes
   !                in the relative weights of various energy terms.
   !       NMRNRG : Calculates the energy/derivatives resulting from the
   !                NMR restraints
   !       MODWT  : Modifies the relative weights of the various energy terms.
   
   !       NMRPRT : Prints the energy and average deviations for the restraints
   !       NDVPRT : Prints the deviation and energy contribution of each
   !                restraint
   
   ! Calls (except to RESTAL) are made thorugh a call to NMRCAL
   ! with the appropriate value of CALTYP (see below). Additionally,
   ! calls to NMRPRT can be made through the entry-point NMRPTX,
   ! and calls to NDVPRT can be made through the entry point NDVPTX.
   
   ! Calls to RESTAL are always made through the entry point RESTLX.
   
   ! The entry point NMRDCP decrements the step-number counter and returns.
   
   ! Author: David A. Pearlman
   ! Date: 7/89
   
   ! INPUT:
   !      CALTYP: (Character*4 variable)
   !              CALTYP = "READ": Call to NMRRED to read weight-change
   !                               and restraint info is made.
   !              CALTYP = "CALC": Call to NMRNRG is made.
   !              CALTYP = "WEIT": Call to MODWT is made.
   !              CALTYP = "PRNT": Call to NMRPRT is made (also done
   !                               every NPRNT steps).
   !              CALTYP = "DPRT": Call to NDVPRT is made.
   !        X(I): Coordinate array. Coordinates (x,y,z) for atom j are packed in
   !              elements 3*(j-1)+1, 3*(j-1)+2 and 3*(j-1)+3.
   !        F(I): Force array. Packed the same as the coordinate array.
   !     NAME(I): Name array. NAME(I) is the name of atom I.
   !   IRSNAM(I): Residue name array. RSNAM(I) is name of residue I.
   !    IPRES(I): Residue pointer array. IPRES(I) points to the first atom
   !              of residue I.
   !   RIMASS(I): The inverse mass of atom I.
   !       TEMP0: (MD only) On the first call to MODWT, should be the target
   !              temperature as specified in the normal input. On this and
   !              subsequent calls to MODWT, the appropriate new target
   !              temperature TEMP0 is returned.
   !    TEMP0LES: (MD only) On the first call to MODWT, should be the target
   !              temperature for LES region as specified in the normal input.
   !              On this and subsequent calls to MODWT, the appropriate
   !              new target temperature TEMP0LES is returned.
   !       TAUTP: same behavior as TEMP0, but for the Berendsen scaling
   !              parmameter.
   !         CUT: On the first call to MODWT, CUT should be the square of the
   !              non-bonded cutoff specified in the normal input. On this and
   !              subsequent calls to MODWT, the square of the appropriate new
   !              target cutoff distance is returned.
   !         NTB: Periodic boundary conditions flag.
   !     WORK(I): Real work array.
   !              If averaging of restraints has not been requested,
   !              WORK must be at least
   !              12*MAXRST + 12*MAXRST*ITIMAV + 2*MAXWT elements long.
   !              The WORK storage must not be modified between calls to NMRCAL.
   
   !              Note that  MAXRST, ITIMAV, MAXWT, and MAXGRP are all set
   !              during the initialization call to RESTAL. ITIMAV is a flag
   !              which is one if time-averaged restraints were requested.
   !    IWORK(I): Integer work array.
   !              IWORK must have at least 12*MAXRST+5*MAXWT+2*MAXGRP elements
   !              The IWORK storage must not be modified between calls to NMRCAL
   !   ISCRTH(I): Integer scratch array.
   !              ISCRTH must be at least NATOM (# of atoms) elements long.
   !              ISCRTH is used only by NMRRED, and can be reused outside of
   !              this routine.
   !      MAXRST: The maximum number of NMR restraints allowed.
   !       MAXWT: The maximum number of relative weight change definitions.
   !          IN: The unit from which the NMR restraints/weight change
   !              definitions will be read.
   !        IOUT: The unit for informational prints.
   !       NPRNT: Every NPRNT steps, info regarding the energies/deviations
   !              of the restraints will be printed to unit IOUT. No print-out
   !              if NPRNT < 0.
   
   !     The following are used only in MODWT calls. They are modified if the
   !     relative weights of the appropriate terms change:
   
   !       RK(I): Force constant array for bonds.
   !       TK(I): Force constant array for angles.
   !       PK(I): Force constant array for torsions. Regular torsions, followed
   !              by impropers.
   !      CN1(I): "A" (r**12) coefficient array for vdw term.
   !      CN2(I): "B" (r**6) coefficient array for vdw term.
   !       AG(I): "A" (r**12) coefficient array for h-bond term.
   !       BG(I): "B" (r**10) coefficient array for h-bond term.
   !       CG(I): Charge on atom I.
   !      NUMBND: Number of bond force constants.
   !      NUMANG: Number of angle force constants.
   !      NUMPHI: Number of torsional force constants.
   !      NIMPRP: Number of "improper" torsional force constants (not included
   !              in NUMPHI).
   !       NTTYP: Number of 6-12 (vdw) parameters
   !         NHB: Number of h-bond parameters.
   !      NCHARG: Number of charge parameters. Equal to the number of atoms for
   !              minimization and normal MD. Equal to 2*number_of_atoms for
   !              perturbed systems (GIBBS).
   !       NATOM: The number of atoms.
   !      NTYPES: The number of types of non-bonded (6-12) parameters.
   !        NRES: The number of residues in the system
   !      RAD(I): The van der Waals radius of atom-type I.
   !      WEL(I): The well depth (epsilon) for atom-type I.
   !    RADHB(I): The mixed radius for h-bond interaction set I.
   !    WELHB(I): The mixed epsilon for h-bond interaction set I.
   !       RWELL: The force constant used for non-bond interactions if
   !              ISFTRP > 0.
   !      ISFTRP: Soft-repulsion flag. If ISFTRP >0, the normal 6-12 vdw
   !              potential is replaced by a soft repulsion potential.
   !              If ISFTRP = 2, hydrogen bond 10-12s are also replaced
   !              by the soft repulsion potential.
   
   !     TGTRMSD: current target RMSD value for targeted MD
   
   
   ! The following are only passed in a call to RESTLX:
   
   !      ITOTST: If > 0, then time-averaged requests will be honored.
   !              Set to 0 to ignore time-averaged requests (e.g. for MIN).
   !         DTT: The integration time-step (ps). It is assumed that the
   !              time-step in constant throughout the run. The time-step
   !              is only required for time-averaged restraints.
   
   ! OUTPUT:
   ! ------
   ! (The following are set when a call to NMRNRG [CALTYP="CALC"] is made).
   
   !     ENMR(5): (1) is the total energy contribution from distance restraints.
   !              (2) is the total energy contribution from angle restraints.
   !              (3) is the total energy contribution from torsion restraints.
   !              (4) is the total energy contribution from plane-point restraints.
   !              (5) is the total energy contribution from plane-plane restraints.
   !   DEVDIS(4): (1) is the average deviation of the distance restraints from
   !                  the "target" value. The average of r2 and r3 (see routine
   !                  NMRNRG) is used as the "target distance"
   !              (2) is the rms deviation of the distance restraints from the
   !                  "target" value. The "target" value is defined as for (1).
   !              (3) is the average deviation of the distance restraints from
   !                  the "target". But in this case, a deviation of 0 is used
   !                  when restraints fall between r2 and r3. For restraints
   !                  outside this range, either r2 or r3 (whichever is closer)
   !                  is used as the target value.
   !              (4) is the rms deviation, calculated using the same deviations
   !                  defined for (3).
   !   DEVANG(4): Same as DEVDIS(4), but for angles. In radians.
   !   DEVTOR(4): Same as DEVDIS(4), but for torsions. In radians.
   !   DEVPLPT(4): Same as DEVDIS(4), but for plane-point restraints. In radians.
   !   DEVPLN(4): Same as DEVDIS(4), but for plane-plane restraints. In radians.
   
   ! The following are set when a call to RESTAL [entry point RESTLX] is made:
   
   !     NMRNUM:
   ! ITMAL(2,6):  Elements (I,1)->(I,5) contain the amount of storage allocated
   !              for time-averaging of bonds, angles, torsions, plane-point and plane-plane restraints (0 if
   !              time-averaging not requested). ITMAL(1,4) contains the total
   !              amount of space required. Generally calculated by routine
   !              RESTAL from the number of MD steps to be run, the numbers
   !              of restraints, and the values of NAVE and NEXACT (see the
   !              headers to routines NMRRED and RESTAL for more details).
   !    MAXRST :  The maximum allowable number of restraints.
   !    ITIMAV :  = 0 if no time-averaged restraints requested.
   !              = 1 if time-averaged restraints requested.
   !    MAXWT  :  The maximum allowable number of weight-change cards.
   !    MAXGRP  : The maximum number of ranges of atoms which can be stored
   !              as group information (only needed for group center-of-mass
   !              coordinates in distance restraints). This value is
   !              computed in the initial call to RESTAL.
   !              [If a non-zero value of MXGRP is passed
   !              to RESTLX, this value will be used as MAXGRP; this
   !              option is not used in the present code.]
   !    INTREQ :  The amount of integer work array (IWORK) storage required.
   !    IRLREQ :  The amount of real work arary (WORK) storage required.
   
   ! Partitioning the work arrays:
   !   The work storage is partitioned here. The following variables
   !   are set to indicate the storage locations of the various arrays
   !   in the work vectors:
   !         J1NMR:  R1NMR(2,I)
   !         J2NMR:  R2NMR(2,I)
   !         J3NMR:  R3NMR(2,I)
   !         J4NMR:  R4NMR(2,I)
   !        JK2NMR:  RK2NMR(2,I)
   !        JK3NMR:  RK3NMR(2,I)
   !        JNMRAT:  NMRAT(16,I)
   !       JRSTTYP:  RESTTYPE(I)
   !        JNMRST:  NMRST(3,I)
   !        JCOMDF:  NMRCOM(2,I)
   !        JSHRTD:  NMRSHB(I)
   !        JFNTYP:  NMRFTY(I)
   !        JGRAVT:  IGRAVT(I)
   !         JALTD:  IALTDIS(I)
   !        JWTNRG:  WTNRG(2,I)
   !         JBAVE:  BAVE(I)
   !         JAAVE:  AAVE(I)
   !        JTAVE1:  TAVE1(I)
   !        JTAVE2:  TAVE2(I)
   !         JBAV0:  BAVE0(I)
   !         JAAV0:  AAVE0(I)
   !        JTAV01:  TAVE01(I)
   !        JTAV02:  TAVE02(I)
   !        JPTAVE:  PTAVE(I)
   !       JPTAVE0:  PTAVE0(I)
   !        JPLAVE:  PLAVE(I)
   !       JPLAVE0:  PLAVE0(I)
   !        JJCOEF:  AJCOEF(3,I)
   !        JWTSTP:  IWTSTP(3,I)
   !        JWTTYP:  IWTTYP(I)
   !          JXPK:  KXPK(I)
   
   ! Other Local Common:
   !         NSTEPL: Local accumulator which keeps track of the number of
   !         calls to CALNRG. On each call, NSTEPL is incremented by 1.
   !         The effective value of NSTEP used is NINT(NSTEPL*STPMLT),
   !         where STPMLT=1.0 if not reset by the user (resulting in
   !         NSTEPL = NSTEP).
   !         Since NSTEPL is incrmented after CALNRG is called, MODWT
   !         should be called *before* CALNRG on any cycle.
   
   !         STPMLT modifies the effective increment in NSTEP per call.
   !         By default, STPMLT is 1.0.
   
   !         DVDIS(2,4),DVANG(2,4),DVTOR(2,4),DVPLPT(2,4),DVPLN(2,4) and EENMR(2,5) keep running totals
   !         of DEVDIS,DEVANG,DEVTOR,DEVPLPT,DEVPLN and ENMR, respectively, in the second
   !         elements (e.g. DVDIS(2,1), DVDIS(2,2),...). The first elements
   !         save the values of these quantities on the last call to NMRNRG.
   
   ! ISCOP1 and ISCOP2 are scratch unit numbers, used when input/output
   ! redirection is requested. ISCOP3 is a dedicated unit number, used for
   ! opening the file which contains periodic "dumps" of the restraint values.
   implicit none
   integer:: iavtyp, ichgwt, idmpav, ier, in, iout, ipower, ipres, &
        iscop1, iscop2, iscop3, iscrth, isftrp, ishrtb, itimav, itotst, &
        iwork, j1nmr, j2nmr, j3nmr, j3vec, j4nmr, jaav0, jaave, jaltd, &
        jbav0, jbave, jcomdf, jconstr, jfntyp, jgdave, jgdave0, jgravt, &
        jjcoef, jk2nmr, jk3nmr, jnmrat, jnmrst, jplave, jplave0, jptave, &
        jptave0, jrsttyp, jrstwt, jshrtd, jtav01, jtav02, jtave1, jtave2, &
        jwtnrg, jwtstp, jwttyp, jxpk, maxgrp, maxrst, maxwt, mxgrp, &
        natom, nave, navint, ncharg, nexact, nhb, nimprp, nmrnum, nprnt, &
        nres, nstep, nstepl, ntb, ntypes, numang, numbnd, numphi
   _REAL_ :: ag, bg, cg, cn1, cn2, cut, devang, devdis, devgendis, &
        devpln, devplpt, devtor, dt, dtt, dvang, dvdis, dvgendis, dvpln, &
        dvplpt, dvtor, eenmr, enmr, f, pk, rad, radhb, ravein, rimass, &
        rk, rwell, stpmlt, tauave, taufin, tautp, temp0, temp0les, &
        tgtrmsd, tk, wel, welhb, work, wtls, x
   
   parameter (iscop1 = 33)
   parameter (iscop2 = 34)
   parameter (iscop3 = 35)
   
   character(len=4) caltyp, name, irsnam
   dimension x(*),f(*),name(*),rimass(*),work(*),iwork(*),iscrth(*)
   dimension rk(*),tk(*),pk(*),cn1(*),cn2(*),ag(*),bg(*),cg(*)
   dimension rad(*),wel(*),irsnam(*),radhb(*),welhb(*),ipres(*)
   dimension enmr(6),devdis(4),devang(4),devtor(4),devplpt(4),devpln(4),devgendis(4)
#  include "nmr.h"
   
   integer B_NMRI,B_NMRR
   parameter (B_NMRI=76)    ! size of NMCLOCI common block
   parameter (B_NMRR=22)    ! size of NMCLOCR common block
   common/nmcloci/nstepl,nstep,j1nmr,j2nmr,j3nmr,j4nmr, &
         jk2nmr,jk3nmr,jnmrat,jrsttyp,jrstwt,jnmrst,jcomdf,jshrtd,jfntyp, &
         jgravt,jwtnrg,jbave,jaave,jtave1,jtave2,jptave,jplave,jgdave,jbav0,jaav0, &
         jtav01,jtav02,jptave0,jplave0,jgdave0,jjcoef,jwtstp,jwttyp,jxpk,ichgwt,nmrnum, &
         ishrtb,nave(6),nexact(6),ipower(6),navint(6), &
         iavtyp(6),idmpav,maxrst,itimav,maxwt,maxgrp,jaltd, &
         j3vec,jconstr
   common/nmclocr/stpmlt,wtls(2),dt,tauave(6),taufin(6),ravein(6)
   common/nmclc2/dvdis(2,4),dvang(2,4),dvtor(2,4),dvplpt(2,4),dvpln(2,4),dvgendis(2,4),eenmr(2,6)

   ! Local Variables
   integer i, j
   integer nstep0
   integer numgrp, itmav
   integer ionce
   save ionce
   data ionce/0/

#ifdef MPI
   integer :: ierr
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif
   
   ! On first non-initialization call to this routine, initialize NSTEPL and set
   ! the WORK/IWORK storage partition.
   
   if (ionce /= 99) then
      nstepl = 0
      nstep = 0
      ! WORK array:
      j1nmr = 1
      j2nmr = j1nmr + 2*maxrst
      j3nmr = j2nmr + 2*maxrst
      j4nmr = j3nmr + 2*maxrst
      jk2nmr = j4nmr + 2*maxrst
      jk3nmr = jk2nmr + 2*maxrst
      jrstwt = jk3nmr + 2*maxrst
      jwtnrg = jrstwt + 4*maxrst
      jbave = jwtnrg + 2*maxwt
      jbav0 = jbave + 2*maxrst*itimav
      jaave = jbav0 + maxrst*itimav
      jaav0 = jaave + 2*maxrst*itimav
      jtave1 = jaav0 + maxrst*itimav
      jtav01 = jtave1 + 2*maxrst*itimav
      jtave2 = jtav01 + maxrst*itimav
      jtav02 = jtave2 + 2*maxrst*itimav
      jptave = jtav02 + maxrst*itimav
      jptave0 = jptave + 2*maxrst*itimav
      jplave = jptave0 + maxrst*itimav
      jplave0 = jplave + 2*maxrst*itimav
      jgdave = jplave0 + maxrst*itimav
      jgdave0 = jgdave + 2*maxrst*itimav
      jjcoef = jgdave0 + maxrst*itimav
      j3vec  =  jjcoef + 3*maxrst       ! constraining vectors
      ! J3VEC require 3*MAXRST storage.
      ! IWORK array:
      jnmrat = 1
      jrsttyp = jnmrat + 16*maxrst
      jnmrst = jrsttyp + maxrst
      jcomdf = jnmrst + 3*maxrst
      jshrtd = jcomdf + 2*maxgrp
      jfntyp = jshrtd + maxrst
      jgravt = jfntyp + maxrst
      jaltd  = jgravt + maxrst
      jwtstp = jaltd  + maxrst
      jwttyp = jwtstp + 3*maxwt
      jxpk   = jwttyp + maxwt
      jconstr = jxpk + 2*maxrst         ! constrain flags
      ! JCONSTR requires MAXRST of storage
      
      do j=1,2
         do i = 1,4
            dvdis(j,i) = 0.0d0
            dvang(j,i) = 0.0d0
            dvtor(j,i) = 0.0d0
            dvplpt(j,i) = 0.0d0
            dvpln(j,i) = 0.0d0
            dvgendis(j,i) = 0.0d0
            eenmr(j,i) = 0.0d0
         end do
         eenmr(j,5) = 0.0d0
         eenmr(j,6) = 0.0d0
      end do
      wtls(1:2) = 1.0d0
      tauave(1:3) = 1.0d+7
      ionce = 99
   end if
   
   select case (caltyp)
    case ('READ')
      call nmrred(x,name,ipres,rimass,work(j1nmr),work(j2nmr), &
            work(j3nmr),work(j4nmr),work(jk2nmr),work(jk3nmr), &
            work(jjcoef),work(j3vec),iwork(jnmrat),iwork(jrsttyp), &
            work(jrstwt),iwork(jnmrst), &
            iwork(jcomdf),iwork(jshrtd),iwork(jfntyp), &
            iwork(jgravt),iwork(jaltd),nmrnum,work(jwtnrg), &
            iwork(jwtstp),iwork(jwttyp),iwork(jxpk),iwork(jconstr), &
            iscrth,ichgwt,ishrtb, &
            natom,nres,nstep0,stpmlt,in,iout,iscop1,iscop2, &
            maxrst,maxwt,maxgrp,nave,nexact,navint,ipower, &
            ravein,taufin,iavtyp,idmpav,itimav)
      if (nstep0 /= 0) nstepl = nstep0
      if (nstep0 /= 0) nstep = nint(nstep0*stpmlt)
    case ('CALC')
      call nmrnrg(x,f,rimass,ntb,nmrnum,nstep,iwork(jnmrat), &
            iwork(jrsttyp),work(jrstwt),iwork(jnmrst),work(j1nmr),&
            work(j2nmr),work(j3nmr), &
            work(j4nmr),work(jk2nmr),work(jk3nmr),work(j3vec), &
            iwork(jcomdf),iwork(jfntyp),iwork(jgravt), &
            iwork(jaltd),iwork(jconstr),work(jbave), &
            work(jbav0),work(jaave),work(jaav0), &
            work(jtave1),work(jtav01),work(jtave2), &
            work(jtav02),work(jptave),work(jptave0),work(jplave), &
            work(jplave0),work(jgdave),work(jgdave0),work(jjcoef),nave,nexact,ipower, &
            tauave,ravein,dt,navint,iavtyp,idmpav,iscop3, &
            enmr,devdis,devang,devtor,devplpt,devpln,devgendis,iout)
      do j = 1,4
         dvdis(1,j) = devdis(j)
         dvang(1,j) = devang(j)
         dvtor(1,j) = devtor(j)
         dvplpt(1,j) = devplpt(j)
         dvpln(1,j) = devpln(j)
         eenmr(1,j) = enmr(j)
         dvdis(2,j) = dvdis(2,j) + devdis(j)
         dvang(2,j) = dvang(2,j) + devang(j)
         dvtor(2,j) = dvtor(2,j) + devtor(j)
         dvplpt(2,j) = dvplpt(2,j) + devplpt(j)
         dvpln(2,j) = dvpln(2,j) + devpln(j)
         eenmr(2,j) = eenmr(2,j) + enmr(j)
      end do
      eenmr(1,5) = enmr(5)
      eenmr(1,6) = enmr(6)
      eenmr(2,5) = eenmr(2,5) + enmr(5)
      eenmr(2,6) = eenmr(2,6) + enmr(6)
      if (nprnt > 0 .and. mod(nstep,abs(nprnt)) == 1) then
        call nmrprt(dvdis,dvang,dvtor,dvplpt,dvpln,dvgendis,eenmr,nstep,iout)
      end if
      nstepl = nstepl + 1
      nstep = nint(nstepl*stpmlt)
    case ('WEIT')
      ! First check/reset short-range distance interactions flags
      call nmrsht(x,nmrnum,rimass,ntb,nstep,iwork(jnmrat), &
            iwork(jrsttyp),iwork(jcomdf),nres,ipres,ishrtb,iwork(jwtstp), &
            work(jwtnrg),work(jk2nmr),work(jk3nmr), &
            iwork(jnmrst),iwork(jshrtd),wtls,iwork(jfntyp), &
            iwork(jgravt),iwork(jaltd),work(jbave), &
            work(jbav0),nave, &
            nexact,ipower,tauave,navint,ravein,dt,iout,1)
      call modwt(work(jwtnrg),iwork(jwtstp),iwork(jwttyp),ichgwt, &
            ishrtb,nstep,temp0,tautp,cut,rk,tk,pk,cn1,cn2,ag,bg,cg, &
            work(jk2nmr),work(jk3nmr),iwork(jshrtd), &
            iwork(jnmrst),numbnd,numang,numphi,nimprp, &
            nhb,ncharg,ntypes,rad,wel,radhb,welhb,rwell,tauave, &
            wtls,isftrp,tgtrmsd,nmrnum, temp0les )

    case ('PRNT') 
      call nmrprt(dvdis,dvang,dvtor,dvplpt,dvpln,dvgendis,eenmr,nstep,iout)
    case ('DVPT') 
      call ndvprt(x,f,name,irsnam,ipres,nres,iscrth,natom, &
            rimass,ntb,nmrnum,nstep,iwork(jnmrat),iwork(jrsttyp), &
            work(jrstwt),iwork(jxpk), &
            iwork(jnmrst),work(j1nmr),work(j2nmr),work(j3nmr), &
            work(j4nmr),work(jk2nmr),work(jk3nmr), &
            iwork(jcomdf),iwork(jfntyp),iwork(jgravt),iwork(jaltd), &
            work(jbave),work(jbav0),work(jaave),work(jaav0), &
            work(jtave1),work(jtav01),work(jtave2), &
            work(jtav02),work(jptave),work(jptave0),work(jplave),work(jplave0), &
            work(jgdave),work(jgdave0),work(jjcoef),nave,nexact,ipower, &
            taufin,ravein,navint,dt,iscop1,iout)
#ifdef MPI
     case ('MPI ')
      
      !     ...when using current AMBER/MPI implementation, only the master
      !     processors reads and processes the input file.  This implies
      !     that only the master processor sets anything in the NMRCAL "READ"
      !     step or in the RESTAL entry point.  Therefore, we need to give
      !     the other processors the data.  This is done by adding the "MPI "
      !     CALTYP for NMRCAL() and calling this *after* startup (in order that
      !     the nmr and locmem pointers, etc be already properly broadcast).
      
      call mpi_bcast(nstepl,B_NMRI,mpi_integer,0,commsander,ierr)
      call mpi_bcast(stpmlt,B_NMRR,MPI_DOUBLE_PRECISION,0,commsander,ierr)
      call mpi_bcast(work(1),jjcoef+3*maxrst,MPI_DOUBLE_PRECISION, &
            0,commsander,ierr)
      call mpi_bcast(iwork(1),jxpk+2*maxrst,mpi_integer,0,commsander,ierr)
      call mpi_bcast(iredir(1),9,mpi_integer,0,commsander,ierr)
#endif
    case default
      write(iout,9000) caltyp
      call mexit(iout, 1)
   end select 
   
   return
   
   ! The following entry-point allows a call to RESTAL without having to
   ! specify lots of unnecessary stuff in the calling list:
   
   entry restlx(in,itotst,mxgrp,dtt,iout,ier)
   call restal(iscop1,in,itotst,ichgwt,nmrnum,itmav,iout,ier,numgrp)
   maxrst = nmrnum + 1 + 5*ier
   maxwt = ichgwt + 1 + 5*ier
   maxgrp = numgrp + 1
   if (mxgrp > 0) maxgrp = mxgrp
   itimav = itmav
   if (itotst == 0) itimav = 0
   intreq = 27*maxrst + 5*maxwt + 2*maxgrp
   irlreq = 22*maxrst + 21*maxrst*itimav + 2*maxwt
   dt = dtt
   return
   
   ! The following entry-point allows a call to NMRPRT without having to specify
   ! lots of unnecessary stuff in the calling list:
   
   entry nmrptx(iout)
   call nmrprt(dvdis,dvang,dvtor,dvplpt,dvpln,dvgendis,eenmr,nstep,iout)
   return
   
   ! The following entry-point allows a call to NDVPRT without having to
   ! specify lots of unnecessary stuff in the calling list:
   
   entry ndvptx(x,f,name,irsnam,ipres,nres,iscrth,natom,rimass, &
         ntb,work,iwork,iout)
   call ndvprt(x,f,name,irsnam,ipres,nres,iscrth,natom, &
            rimass,ntb,nmrnum,nstep,iwork(jnmrat),iwork(jrsttyp), &
            work(jrstwt),iwork(jxpk), &
            iwork(jnmrst),work(j1nmr),work(j2nmr),work(j3nmr), &
            work(j4nmr),work(jk2nmr),work(jk3nmr), &
            iwork(jcomdf),iwork(jfntyp),iwork(jgravt),iwork(jaltd), &
            work(jbave),work(jbav0),work(jaave),work(jaav0), &
            work(jtave1),work(jtav01),work(jtave2), &
            work(jtav02),work(jptave),work(jptave0),work(jplave),work(jplave0), &
            work(jgdave),work(jgdave0),work(jjcoef),nave,nexact,ipower, &
            taufin,ravein,navint,dt,iscop1,iout)
   return
   
   ! The following entry point allows the local step number counter,
   ! NSTEPL, to be decremented by one. This is necessary when a call to FORCE
   ! is made which does not correspond to a dynamics step:
   ! The accumuated deviation and energy counters are also decremented
   ! by the distances on the last call.
   
   entry nmrdcp
   nstepl = nstepl - 1
   nstep = nint(nstepl*stpmlt)
   do j = 1,4
      dvdis(2,j) = dvdis(2,j) - dvdis(1,j)
      dvang(2,j) = dvang(2,j) - dvang(1,j)
      dvtor(2,j) = dvtor(2,j) - dvtor(1,j)
      dvplpt(2,j) = dvplpt(2,j) - dvplpt(1,j)
      dvpln(2,j) = dvpln(2,j) - dvpln(1,j)
      dvgendis(2,j) = dvgendis(2,j) - dvgendis(1,j)
      eenmr(2,j) = eenmr(2,j) - eenmr(1,j)
   end do
   eenmr(2,5) = eenmr(2,5) - eenmr(2,5)
   eenmr(2,6) = eenmr(2,6) - eenmr(2,6)
   return
   
   9000 format(' Error: Unrecognized option passed to NMRCAL: CALTYP = ', &
         a4)
end subroutine nmrcal 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine restal here]
subroutine restal(iscop1,in,itotst,ichgwt,nmrnum,itimav, &
      iout,ierr,numgrp)

   ! Subroutine RESTraint Allocation
   
   ! This subroutine does a quick read of the weight-change/restraint
   ! cards used for NMR restraints, and returns the numbers of each (NMRNUM
   ! and ICHGWT, respectively).
   
   ! These quantities are required in setting up the dynamic memory
   ! allocation in routine LOCMEM
   
   ! It is expected that this routine will be called before any other of the
   ! NMR suite.
   
   ! Author: David A. Pearlman
   ! Date: 5/90
   
   ! INPUT:
   !      ISCOP1 : Scratch units for informational reads/writes. Used when
   !               redirection of input is requested.
   !          IN : Unit from which the restraints/weight changes will be read.
   !      ITOTST : Total number of dynamics steps to be run.
   !        IOUT : Unit for informational prints.
   
   ! OUTPUT:
   !      ICHGWT : Number of change parameters read.
   !      NMRNUM : The total number of NMR restraints defined.
   !      NUMGRP : The total number of NMR groups defined.
   !      ITIMAV : = 1 If any time-averaged restraints were requested.
   !               = 0 otherwise.
   !        IERR : Returned as 0 if no problems reading input. Returned as
   !               1 if there were problems. Actual informational messages to
   !               user about problems will be generated during the call to
   !               NMRRED
   
   use nmr
   implicit none

   integer:: i, iaaver, ialtd, iat, iat1, iat2, iat3, iat4, iat5, &
        iat6, iat7, iat8, ibaver, ibeg, ichgwt, iconstr, ierr, ifind, &
        ifntyp, iformw, ifvari, ig, igr1, igr2, igr3, igr4, igr5, igr6, &
        igr7, igr8, igrarr, iin, iinc, ilengt, imat, imult, in, iout, &
        ipaver, ir6, iredir, iresid, irstyp, iscop1, istep1, istep2, &
        istrt, itimav, itotst, ixpk, j, m, maxigr, n, ninc, nmrnum, nres, &
        nstep1, nstep2, numgrp, nxpk
   _REAL_ :: r1, r1a, r2, r2a, r3, r3a, r4, r4a, rjcoef, rk2, rk2a, &
        rk3, rk3a, value1, value2

   integer iii,jjj
   character(len=80) aline,redir
   character(len=10) redirc(9)
   character(len=8) type,flag(6)
   character(len=4) atnam(8),grnam1,grnam2,grnam3,grnam4,grnam5,grnam6,grnam7,grnam8, &
                    grnamarr
   dimension redir(9),iredir(9)
   dimension iat(8),rjcoef(3)
   character(256) restraint   
   character(4) dumnam(2)
   integer dumipres
   dimension dumipres(2)
   _REAL_ rstwt(4)
   _REAL_  r0
   _REAL_  k0
   _REAL_  r0a
   _REAL_  k0a

   parameter (maxigr=200)
   dimension igr1(maxigr),igr2(maxigr),grnam1(maxigr),grnam2(maxigr), &
             igr3(maxigr),igr4(maxigr),grnam3(maxigr),grnam4(maxigr), &
             igr5(maxigr),igr6(maxigr),grnam5(maxigr),grnam6(maxigr), &
             igr7(maxigr),igr8(maxigr),grnam7(maxigr),grnam8(maxigr) 
   dimension igrarr(maxigr,8), grnamarr(maxigr,8)
   !  --- old formatted input:
   !     DIMENSION IGR(16)
   
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
   
   data redirc/'LISTIN    ' , 'LISTOUT   ' , 'DISANG    ', &
         'NOESY     ' , 'SHIFTS    ' , 'DUMPAVE   ', &
         'PCSHIFT   ' , 'DIPOLE    ' , 'CSA       ' /
   
   data flag/'DISAVE  ' , 'ANGAVE  ' , 'TORAVE  ' , 'PLPTAVE ', 'PLNAVE  ', 'GDISAVE '/
   
   ichgwt = 0
   nmrnum = 0
   numgrp = 0
   ibaver = 0
   iaaver = 0
   ipaver = 0
   itimav = 0
   iat1 = 0
   restraint = ' '
   ierr = 0
   iin = in
   ibeg = 1
   do i = 1,8
      iredir(i) = 0
   end do
   
   dumnam(1) = '    '
   dumnam(2) = '    '
   dumipres = 0
   
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
   
   9 continue
   do i = ibeg,99999
      
      ! Branch depending on whether formatted weight changes (IFORMW = 1) or
      ! namelist weight changes are being read:
      
      if (iformw /= 1) then
         
         ! NAMELIST-type read:
         ! ==================
         
         ! Set the default values of the parameters which can be set by the namelist
         
         istep1 = 0
         istep2 = 0
         value1 = 0.0d0
         value2 = 0.0d0
         iinc = 0
         imult = 0
         
         ! Look for "wt" namelist. If not found, stop with error.
         
         call nmlsrc('wt',iin,ifind)
         if (ifind == 0) goto 1001
         
         read(iin,nml=wt)
         
         ! Formatted-type weight card read:
         ! ===============================
         
      else
         read(iin,9005,end=1001,err=1005) type,istep1,istep2, &
               value1,value2,iinc,imult
      end if
      
      ! End of weight card read. Now process given options:
      ! ==================================================
      
      ! Do not count comment lines:
      
      if (type(1:1) == '#') then
         cycle
         
         ! If ITYPE="END", then go to redirection reading section
         
      else if (type(1:3) == 'END' .or. type(1:3) == 'end') then
         exit
      else
         
         ! Check for cards requesting time-averaging
         
         do j = 1,6
            if (type == flag(j)) then
               itimav = 1
            end if
         end do
      end if
      
      ichgwt = ichgwt + 1
   end do
   
   ! Now read the redirection cards, if any. The first non-blank line encountered
   ! is assumed to be the start of the redirection section. If it is not a known
   ! option, assume no redirection cards were given.
   
   
   20 read(iin,9004,end=25) aline
   
   ! skip blank lines
   
   do j = 1,len(aline)
      if (aline(j:j) /= ' ') goto 29
   end do
   goto 20
   29 continue
   
   do m=1,9
      do n=10,1,-1
         if (redirc(m)(n:n) /= ' ') exit
      end do
      call chklin(aline,redirc(m),n,imat,istrt,ilengt,iout)
      if (imat /= 0) then
         redir(m) = aline(istrt:istrt+ilengt-1)
         iredir(m) = ilengt
         goto 20
      end if
   end do
   
   backspace(iin)
   
   ! Now read the restraints.
   
   25 continue
   
   ! If restraint input has been redirected, open the appropriate file
   
   if (iredir(3) /= 0) then
      call amopen(iscop1,redir(3)(1:iredir(3)),'O','F','R')
      iin = iscop1
   else
      goto 100    ! no DISANG restraints to read
   end if
   
   ! Count the restraints:
   
   do i=1,999999
   
      ! rezero iat
      do iii = 1,8
        iat(iii) = 0
      end do
      
      do jjj=1,maxigr
        igr1(jjj) = 0
        igr2(jjj) = 0
        igr3(jjj) = 0
        igr4(jjj) = 0
        igr5(jjj) = 0
        igr6(jjj) = 0
        igr7(jjj) = 0
        igr8(jjj) = 0
      end do
      
      ! Look for "rst" namelist. If not found, assume we are done reading them
      
      call nmlsrc('rst',iin,ifind)
      if (ifind == 0) goto 100
      
      read (iin,nml=rst,end=100)
      
      do jjj=1,maxigr
        igrarr(jjj,1) = igr1(jjj)
        igrarr(jjj,2) = igr2(jjj)
        igrarr(jjj,3) = igr3(jjj)
        igrarr(jjj,4) = igr4(jjj)
        igrarr(jjj,5) = igr5(jjj)
        igrarr(jjj,6) = igr6(jjj)
        igrarr(jjj,7) = igr7(jjj)
        igrarr(jjj,8) = igr8(jjj)
      end do
      
      ! If iat is specified in the old style, don't parse the restraint variable.
      
      if (restraint /= ' ' .and. iat1 ==0) then
            call parsrest(restraint,iat,igrarr,rstwt,dumnam,dumipres,nres,iout)
      end if
      
      do iii=1,8
        if (iat(iii) < 0) then
          numgrp = numgrp + 1
          do ig=2,maxigr
            if (igrarr(ig,iii) == 0) exit
            if (igrarr(ig,iii) /= igrarr(ig-1,iii)+1 .or. restraint /= ' ') numgrp = numgrp + 1  
            ! Since restraint is not parsed properly for atoms specified with atom names at this stage, 
            ! I'm going to be cautious and allocate one group for each entry in IGR, even if some will 
            ! not ultimately be used.
          end do
        end if
      end do
      34 nmrnum = nmrnum + 1

      ! If IAT1=0, assume this is the final restraint:
      if (iat1 == 0) goto 100
   end do
   
   100 continue
   
   if (iredir(3) /= 0) close(iin)
   return
   
   ! Errors:
   ! Do not report these here, just return with IERR=1 flag. The actual reporting
   ! and program termination will occur when routine NMRRED is called.
   
   ! End-of-file read errors:
   
   1001 ierr = 1
   return
   
   ! Other read errors (probably format mismatch)
   ! (If "#" is in first column, this is a comment line):
   
   1005 backspace(iin)
   read(iin,9004) aline
   if (aline(1:1) == '#') then
      ibeg = i+1
      goto 9
   else
      ierr = 1
      return
   end if

   !9001 FORMAT(16I5)
   9004 format(a80)
   9005 format(a8,2i7,2f12.6,2i7)
   !9020 FORMAT(11I5)
   !9029 FORMAT(4(1X,A4))
   !9030 FORMAT(6F12.6)
end subroutine restal 
