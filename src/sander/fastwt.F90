! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

module fastwt

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fastwat here]
subroutine fastwat(igraph,nres,ipres,lbres, &
      nbonh,nbona,ib,jb,ibelly,igrp, &
      iwtnm,iowtnm,ihwtnm,jfastw,ifstwt, &
      ifstwr,ibgwat,ienwat,ibgion,ienion,iorwat, &
      iout,natom)
   
   
   ! Subroutine FAST WaTer
   
   ! This routine determines which shaken bonds correspond to 3-point
   ! solvent molecules (typically TIP3P or SPC waters). For these bonds,
   ! IFSTWT(I) is set to 1. For all other bonds, IFSTWT(I) is set to 0.
   
   ! Waters are defined as residues that meet the following criteria:
   !   1) residue name = watnam
   !   2) atom names are owatnm, hwatnm(1) and hwatnm(2)
   
   ! Author: David A. Pearlman
   ! Date: 12/92
   
   ! Input:
   
   ! NATOM: Number of atoms
   ! IGRAPH(I): The name of atom I (I4)
   ! NRES: The number of residues
   ! IPRES(I): IPRES(I)->IPRES(I+1)-1 are the atoms of residue I
   ! LBRES(I): The name of residue I
   ! NBONH: The number of bonds to hydrogens
   ! NBONA: The number of bonds to non-hydrogens
   ! IB(I)
   ! JB(I): The two atoms of bond I. Stored as 3*(atom_#-1)
   ! IBELLY: > 0 if belly run is being performed.
   ! IGRP(I): > 0 if atom I is part of moving belly (only if IBELLY >0).
   ! IWTNM: The name of residues to be considered water (I)
   ! IOWTNM: The name of oxygen atom in the water residues (I)
   ! IHWTNM(2): The name of 2 hydrogen atoms in the water residues (I)
   ! JFASTW: = 0: Use fast water routine (possibly with modified names)
   !           4: Do not use fast waters anywhere.
   
   ! Output:
   ! IFSTWT(I): = 0 if bond should be constrained by standard routine.
   !            = 1 if bond should be constrained by fast 3-point const. routine
   ! IFSTWR(I): = 0 if bonds of residue are constrained by standard routine.
   !            = 1 if bonds of residue are constrained by fast 3-point const.rout.
   
   ! IBGWAT: The first water residue found.
   ! IENWAT: The last water residue found, if all residues IBGWAT->IENWAT are
   !         waters. If the waters are not contiguous in the list, IENWAT
   !         will be 0.
   ! IORWAT: The position of the oxygen in each residue (1, 2, or 3).
   use qmmm_module, only : qmmm_nml, qmmm_struct 
   implicit none
   character(len=4) igraph,lbres,iwtnm,iowtnm,ihwtnm
   integer nres,ipres
   integer nbonh,nbona,ib,jb,ibelly,igrp
   integer jfastw,ifstwt
   integer ifstwr,ibgwat,ienwat,ibgion,ienion,iorwat
   integer iout
   integer natom,nrbi,jj
   dimension igraph(*),ipres(*),lbres(*),ib(*),jb(*)
   dimension igrp(*)
   dimension ifstwt(*),ifstwr(*),ihwtnm(*)
   integer numfst
   common/fastw/numfst
   logical all3point
   integer i, j, i1, i2
   integer iat, ihf1, ihf2
   integer iof, iofm
   integer ifind
   integer ncnt
   
   integer :: rb(3,nres), atrn(natom), nrb(nres)

   ! Initialize IFSTWT(I),IFSTWR(I)
      
   do i = 1,nbonh+nbona
      ifstwt(i) = 0
   end do
      
   do i = 1,nres
      ifstwr(i) = 0
   end do
      
   ! Search the list of residues for those a) with 3 atoms; 2) with the
   ! appropriate residue name; 3) with the appropriate atom names
   
   ibgwat = 0
   ienwat = 0
   
   numfst = 0
   
   ! If JFASTW=4 user has requested we not use fast water routine in any case.
   ! simply report no fast waters defined and then return.
      
   if (jfastw == 4) then
      write(iout,9001) numfst
      return
   end if
      
   ! If using belly (IBELLY >0), skip any residue where all atoms of residue
   ! cannot move (where all atoms do not have IGRP(I)>0).
   
   !    setting up bond list for waters
   !       nrb(i)=number of bonds inside this residue i
   !       rb(j,i)=bond entry number for the first three bonds
   !               j=1,2,3  for residue i
   !       atrn(j)=residue number of atom j
      
   do i=1,nres
      nrb(i)=0
      rb(1,i)=0
      rb(2,i)=0
      rb(3,i)=0
   end do
      
   do i=1,nres-1
      do j=ipres(i),ipres(i+1)-1
         atrn(j)=i
      end do
   end do
      
   do j=ipres(nres),natom
      atrn(j)=nres
   end do
      
   do  j = 1,nbonh+nbona
      nrbi=atrn(ib(j)/3+1)
      if(nrbi == atrn(jb(j)/3+1))then
         nrb(nrbi)=nrb(nrbi)+1
         if(nrb(nrbi) < 4)rb(nrb(nrbi),nrbi)=j
      end if
   end do
      
   all3point = .true.
   nres_loop: do i = 1,nres
      
      !       -- screen out residue if wrong resname
      
      if (lbres(i) /= iwtnm) cycle nres_loop

      !       -- check to see if waters have exactly three atoms: if this
      !          is true for all waters, and if the waters are all
      !          contiguous, we can use a single call to setlep() to shake
      !          all of the waters.  If not, we have to generate a call to
      !          setlep() for each water molecule.  For example, TIP4P
      !          and TIP5P waters can use setlep(), but only in the
      !          one-molecule-at-a-time mode.

      if (ipres(i+1)-ipres(i) /= 3) all3point = .false.
      
      !       -- screen out residue if an atom is fixed (belly) or
      !          if an atom name does not match or is duplicated
      !          - also note atom #s
      
      ihf1 = 0
      ihf2 = 0
      iof = 0
      do  j = 1,3
         iat = ipres(i)+j-1
         if (ibelly > 0 .and. igrp(iat) <= 0 ) cycle nres_loop !Skip this residue - move to next
         !If this water residue is currently in the QM region
         !do not flag it for fast water (quick3) treatment
         if (qmmm_nml%ifqnt) then
            if (qmmm_struct%atom_mask(iat)) cycle nres_loop
         end if
         if (igraph(iat) == iowtnm) then
            if (iof /= 0) cycle nres_loop
            iof = 3*(iat-1)
            iofm = iat
         else if (igraph(iat) == ihwtnm(1)) then
            if (ihf1 /= 0) cycle nres_loop
            ihf1 = 3*(iat-1)
         else if (igraph(iat) == ihwtnm(2)) then
            if (ihf2 /= 0) cycle nres_loop
            ihf2 = 3*(iat-1)
         else
            cycle nres_loop
         end if
      end do
      
      ! Set IBGWAT & IENWAT. At the end, IBGWAT will be the first water
      ! residue found. If all water residues are contiguous in the list,
      ! IENWAT will = the number of the last water, otherwise IENWAT will
      ! be 0. Also set IORWAT when the first water is found
      
      if (ibgwat == 0) then
         ibgwat = i
         if (iofm == ipres(i)) then
            iorwat = 1
         else if (iofm == ipres(i)+1) then
            iorwat = 2
         else
            iorwat = 3
         end if
      end if
      
      if (ienwat == i-1 .or. ibgwat == i) then
         ienwat = i
      else
         ienwat = 0
      end if
      
      ! If we get here, everything is OK. Now we need to pick the bonds out
      ! of the list which comprise this water molecule and flag them. Note
      ! that if one atom of a bond is in one water molecule and one is in
      ! another water molecule, this bond is not flagged. Assuming water
      ! molecules have been properly defined, this case should never occur...:
      
      ifstwr(i) = 1
      ifind = 0
      do  jj = 1,3
         j=rb(jj,i)
         i1 = 0
         i2 = 0
         if (ib(j) == iof .or. ib(j) == ihf1 .or. ib(j) == ihf2) i1 = 1
         if (jb(j) == iof .or. jb(j) == ihf1 .or. jb(j) == ihf2) i2 = 1
         if (i1 == 1 .and. i2 == 1) then
            ifstwt(j) = 1
            ifind = ifind + 1
         end if
      end do
      
      ! If IFIND.NE.3, there is a problem (there should be 3 and only 3 bonds
      ! defining the constrained triangle). Stop with error message
      
      if (ifind /= 3) then
         write(iout,9000) i,ifind
         call mexit(6,1)
      end if
      numfst = numfst + 1
   end do nres_loop
   if ( .not. all3point ) ienwat = 0
   write(iout,9001) numfst
   
   if(ibgwat /= 0 .and. ienwat /= 0) then
      ! HG Try to locate monoatomic ions that come _before_ the water molecules
      !   (and after the "solute")
      !    Assumes that all monoatomic ions "come in one stretch"
      !    Stores the results in ibgion and ienion
      ibgion = 0; ienion = 0
      do i=ibgwat-1, 1, -1
         ncnt = ipres(i+1) - ipres(i)
         if(ncnt == 1) then
            ! Found residue with only one atom -> assuming monoatomic ion
            if(ienion == 0) ienion = ibgwat - 1
            ibgion = i
         else
            ! Found residue with more than one atom -> exiting search
            exit
         end if
      end do
      if(ibgion /= 0 .and. ienion /= 0) then
         ! Output result
         !!!write(iout,9002) ienion - ibgion + 1
      end if
   end if
   
   return
   
   ! Format statements:
   
   9000 format(' Error: A residue defined as a "fast 3-point water"', &
         /,8x,'is not defined by a triangle of three bonds.' &
         /,8x,'Residue ',i8,' contains ',i8,' bonds.')
   9001 format(' Number of triangulated 3-point waters found: ',i8)
 9002   format(' Number of monoatomic ions found: ',i8)
end subroutine fastwat 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ determines how to call setlep(), the 3-point water shake routine
subroutine quick3(x0, xh, ifstwr, natom, nres, ipres)
   
   ! Subroutine QUICK 3 point water
   
   ! This routine makes the calls to carry out the fast, analytic 3-point
   ! water constraints.
   
   ! Author: David A. Pearlman
   ! Date: 12/92
   
   ! INPUT:
   ! X0(I): Coordinate array corresponding to time t-dt/2. Coordinates
   !        are packed such that x(atom i), y(atom i) and z(atom i) are
   !        found at X0(3*(i-1)+1), X0(3*(i-1)+2) and X0(3*(i-1)+3), respectively.
   
   ! XH(I): Coordinate array corresponding to time t+dt/2, but not corrected
   !        (on input) for internal constraints, i.e. an unconstrained move.
   !        Packing is the same as for X0.
   
   ! IBGWAT
   ! IENWAT: If IENWAT>0, then all waters to be constrained by the
   !        quick formula are contiguous, and are in the residue
   !        range IBGWAT->IENWAT. If IBGWAT=0, no waters exist for
   !        quick constraints. If IENWAT=0, the waters are not contiguous
   !        in the residue array.
   
   ! IFSTWR(I): If = 0, the bonds of this residue are constrained by
   !        SHAKE/TORCON. If IFSTWR(I)=1, the bonds of the residue are constrained
   !        by the fast, analytic 3-point routine.
   
   ! NATOM: Number of atoms in the system
   ! NRES:  Number of residues in the system.
   
   ! IPRES(I): IPRES(I): IPRES(I+1)-1 are the atoms in residue I.
   
   ! IORWAT: THe position of the oxygen atom in each water (1,2, or 3).
   
   ! RBTARG(9): The five geometry parameters used by the fast water routine
   !         to impose constraints, the masses of the atoms in the water,
   !         and the total mass of the water.
   
   implicit none
   integer ifstwr, natom, nres, ipres
   _REAL_  x0, xh
   dimension x0(*),xh(*),ifstwr(*),ipres(*)
   integer i, iend, istart
   
#ifdef MPI
   
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
   
#endif
#  include "md.h"
   
   ! If IBGWAT=0, no constraints to be done by fast routine. Return now.
   
   if (ibgwat <= 0) return
   
   ! If IENWAT > 0, then all the waters are contiguous in the residue array.
   ! We can use one call to SETLEP for all of them. This will be faster.
   
   if (ienwat > 0) then
      
#ifdef MPI
      
      if ( mpi_orig ) then
         istart = ipres(ibgwat)
         iend = ipres(ienwat)+2
      else
         istart = max(ipres(ibgwat),iparpt(mytaskid) + 1)
         iend = min(ipres(ienwat)+2,iparpt(mytaskid+1))
      end if
      ! write(0,*) 'calling setlep',mytaskid,istart,iend,ibgwat,ienwat
      
#else
      istart = ipres(ibgwat)
      iend = ipres(ienwat)+2
#endif
      
      call setlep (istart,iend,x0,xh,natom,iorwat,rbtarg)
      
      
      ! If IENWAT=0, then the waters are not contiguous. We need to loop over
      ! all the residues and call SETLEP to perform the constraints on the 
      ! waters one at a time. On vector machines, this may be slower than the
      ! single call above.
      
   else
      do i = 1,nres
         if (ifstwr(i) > 0) then
            istart = ipres(i)
            iend = ipres(i) + 2
#ifdef MPI
            if( mpi_orig ) then
               call setlep (istart,iend,x0,xh,natom,iorwat,rbtarg)
            else
               if (istart >= iparpt(mytaskid)+1 .and. &
                    iend <= iparpt(mytaskid+1)) then
                  call setlep (istart,iend,x0,xh,natom,iorwat,rbtarg)
               end if
            end if
#else
            call setlep (istart,iend,x0,xh,natom,iorwat,rbtarg)
#endif
         end if
      end do
   end if  ! (ienwat > 0)
   return
end subroutine quick3 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Apply SHAKE to a three-point water model
subroutine setlep (iwatat, natom, x0, x1, mxatm, iorwat, rbtarg)
   !*****************************************************************
   !                                                               **
   !    Subroutine : setlep - reset positions of TIP3P waters      **
   !    Author : Shuichi Miyamoto                                  **
   !                                                               **
   !    Reference for the SETTLE algorithm                         **
   !      S. Miyamoto et al., J. Comp. Chem.,  13, 952 (1992)      **
   !                                                               **
   !                                                               **
   !*****************************************************************
   implicit none
   integer iwatat, natom, mxatm, iorwat
   _REAL_  x0, x1, rbtarg
   dimension  x0(3,mxatm),x1(3,mxatm),rbtarg(9)

   !*****************************************************************
   ! Most of the variable names correspond to those in the reference.
   
   !     iwatat   : pointer for the first atom of waters
   !     natom    : No. of total atoms in system
   !     x0(3,i)  : position at previous time step (t0) (input)
   !     x1(3,i)  : position at present  time step (t0 + dt)
   !              : position before applying constraints (input)
   !              : position after  applying constraints (output)
   !     mxatm    : maximum No. of atoms in system
   !     iorwat   : position of oxygen in water (1, 2, or 3)
   
   !     x,y,zcom       : center of mass
   !     x,y,zaksX(YZ)d : axis vectors of the alternative orthogonal
   !                      coordinate system Xprime Yprime Zprime
   !     trns..         : matrix of orthogonal transformation
   
   !     wo,wh : mass of oxygen and hydrogen
   !     wohh  : mass of water (H2O)
   !     hhhh  : rc2*rc2 (square of HH distance)
   
   !*****************************************************************
   _REAL_  ra, onera
   _REAL_  rb
   _REAL_  rc
   _REAL_  rc2
   _REAL_  hhhh
   _REAL_  wo
   _REAL_  wh
   _REAL_ onewohh
   !     parameter (wo = 16.000d0 , wh = 1.008d0 , wohh = 18.016d0)
   !     parameter (ra = 0.0655822665508295d0,
   !    .           rb = 0.5204941789748370d0,
   !    .           rc = 0.75680d0,
   !    .          rc2 = 1.51360d0)
   !     parameter (hhhh = 2.29098496d0)
   
   integer i
   integer ind1, ind2, ind3
   _REAL_  axlng, aylng, azlng
   _REAL_  alpha, beta, gamma, al2be2
   _REAL_  cosphi, cospsi, sinphi, sinpsi, costhe, sinthe
   _REAL_  deltx, hh2
   _REAL_  trns11, trns21, trns31
   _REAL_  trns12, trns22, trns32
   _REAL_  trns13, trns23, trns33
   _REAL_  xb0, yb0, zb0, xc0, yc0, zc0
   _REAL_  xb0d, yb0d, zb0d, xc0d, yc0d, zc0d
   _REAL_  xa1, ya1, za1, xb1, yb1, zb1, xc1, yc1, zc1
   _REAL_  xa1d,ya1d,za1d,xb1d,yb1d,zb1d,xc1d,yc1d,zc1d
   _REAL_  xa2d,ya2d,za2d,xb2d,yb2d,zb2d,xc2d,yc2d,zc2d
   _REAL_  xa3, ya3, za3, xb3, yb3, zb3, xc3, yc3, zc3
   _REAL_  xa3d,ya3d,za3d,xb3d,yb3d,zb3d,xc3d,yc3d,zc3d
   _REAL_  xb2d2
   _REAL_  xcom, ycom, zcom
   _REAL_  xakszd, yakszd, zakszd, xaksxd, yaksxd, zaksxd
   _REAL_  xaksyd, yaksyd, zaksyd
   
   !        if ( abs(za1d) .ge. ra) then

   ra = rbtarg(1)
   rb = rbtarg(2)
   rc = rbtarg(3)
   rc2 = rbtarg(4)
   hhhh = rbtarg(5)
   wo = rbtarg(6)
   wh = rbtarg(7)
   onewohh = rbtarg(8)
   onera = rbtarg(9)
   
   if (iorwat == 1) then
      ind1 = 0
      ind2 = 1
      ind3 = 2
   else if (iorwat == 2) then
      ind1 = 1
      ind2 = 2
      ind3 = 0
   else
      ind1 = 2
      ind2 = 0
      ind3 = 1
   end if
   
   !forcevector
   do i=iwatat,natom,3
      !                                                --- Step1  A1_prime ---
      
      xb0 = x0(1,i+ind2) - x0(1,i+ind1)
      yb0 = x0(2,i+ind2) - x0(2,i+ind1)
      zb0 = x0(3,i+ind2) - x0(3,i+ind1)
      xc0 = x0(1,i+ind3) - x0(1,i+ind1)
      yc0 = x0(2,i+ind3) - x0(2,i+ind1)
      zc0 = x0(3,i+ind3) - x0(3,i+ind1)
      
      xcom = ( x1(1,i+ind1)*wo + (x1(1,i+ind2) + x1(1,i+ind3)) * wh ) &
            * onewohh
      ycom = ( x1(2,i+ind1)*wo + (x1(2,i+ind2) + x1(2,i+ind3)) * wh ) &
            * onewohh
      zcom = ( x1(3,i+ind1)*wo + (x1(3,i+ind2) + x1(3,i+ind3)) * wh ) &
            * onewohh
      
      xa1 = x1(1,i+ind1) - xcom
      ya1 = x1(2,i+ind1) - ycom
      za1 = x1(3,i+ind1) - zcom
      xb1 = x1(1,i+ind2) - xcom
      yb1 = x1(2,i+ind2) - ycom
      zb1 = x1(3,i+ind2) - zcom
      xc1 = x1(1,i+ind3) - xcom
      yc1 = x1(2,i+ind3) - ycom
      zc1 = x1(3,i+ind3) - zcom
      
      xakszd = yb0*zc0 - zb0*yc0
      yakszd = zb0*xc0 - xb0*zc0
      zakszd = xb0*yc0 - yb0*xc0
      xaksxd = ya1*zakszd - za1*yakszd
      yaksxd = za1*xakszd - xa1*zakszd
      zaksxd = xa1*yakszd - ya1*xakszd
      xaksyd = yakszd*zaksxd - zakszd*yaksxd
      yaksyd = zakszd*xaksxd - xakszd*zaksxd
      zaksyd = xakszd*yaksxd - yakszd*xaksxd
      
      axlng = 1.0d0 / sqrt ( xaksxd * xaksxd + yaksxd * yaksxd &
            + zaksxd * zaksxd )
      aylng = 1.0d0 / sqrt ( xaksyd * xaksyd + yaksyd * yaksyd &
            + zaksyd * zaksyd )
      azlng = 1.0d0 / sqrt ( xakszd * xakszd + yakszd * yakszd &
            + zakszd * zakszd )
      trns11 = xaksxd * axlng
      trns21 = yaksxd * axlng
      trns31 = zaksxd * axlng
      trns12 = xaksyd * aylng
      trns22 = yaksyd * aylng
      trns32 = zaksyd * aylng
      trns13 = xakszd * azlng
      trns23 = yakszd * azlng
      trns33 = zakszd * azlng
      
      xb0d = trns11*xb0 + trns21*yb0 + trns31*zb0
      yb0d = trns12*xb0 + trns22*yb0 + trns32*zb0
      xc0d = trns11*xc0 + trns21*yc0 + trns31*zc0
      yc0d = trns12*xc0 + trns22*yc0 + trns32*zc0
      xa1d = trns11*xa1 + trns21*ya1 + trns31*za1
      ya1d = trns12*xa1 + trns22*ya1 + trns32*za1
      za1d = trns13*xa1 + trns23*ya1 + trns33*za1
      xb1d = trns11*xb1 + trns21*yb1 + trns31*zb1
      yb1d = trns12*xb1 + trns22*yb1 + trns32*zb1
      zb1d = trns13*xb1 + trns23*yb1 + trns33*zb1
      xc1d = trns11*xc1 + trns21*yc1 + trns31*zc1
      yc1d = trns12*xc1 + trns22*yc1 + trns32*zc1
      zc1d = trns13*xc1 + trns23*yc1 + trns33*zc1
      
      !        if ( abs(za1d) .ge. ra) then
      !          write  (6,699)
      !  699     format (5x,' ### SETLEP : deviation is too big !! ')
      !          call exit
      !        endif
      !                                                --- Step2  A2_prime ---
      
      sinphi = za1d * onera
      cosphi = sqrt (1.d0 - sinphi*sinphi)
      sinpsi = ( zb1d - zc1d ) / (rc2 * cosphi)
      cospsi = sqrt (1.d0 - sinpsi*sinpsi)

      ya2d =   ra * cosphi
      xb2d = - rc * cospsi
      !        xc2d =   rc * cospsi
      yb2d = - rb * cosphi - rc *sinpsi * sinphi
      yc2d = - rb * cosphi + rc *sinpsi * sinphi
      xb2d2 = xb2d * xb2d
      hh2 = 4.d0 * xb2d2 + (yb2d-yc2d) * (yb2d-yc2d) &
            + (zb1d-zc1d) * (zb1d-zc1d)
      deltx = 2.d0 * xb2d + sqrt ( 4.d0 * xb2d2 - hh2 + hhhh )
      xb2d = xb2d - deltx * 0.5d0
      !        xc2d = xc2d + deltx * 0.5d0
      
      !                                                --- Step3  al,be,ga ---
      
      alpha = ( xb2d * (xb0d-xc0d) + yb0d * yb2d + yc0d * yc2d )
      beta = ( xb2d * (yc0d-yb0d) + xb0d * yb2d + xc0d * yc2d )
      gamma = xb0d * yb1d - xb1d * yb0d + xc0d * yc1d - xc1d * yc0d
      
      al2be2 = alpha * alpha + beta * beta
      sinthe = ( alpha*gamma - beta * sqrt ( al2be2 - gamma * gamma ) ) &
            / al2be2
      
      !                                                --- Step4  A3_prime ---
      
      costhe = sqrt (1.d0 - sinthe * sinthe )
      xa3d = - ya2d * sinthe
      ya3d =   ya2d * costhe
      za3d = za1d
      xb3d =   xb2d * costhe - yb2d * sinthe
      yb3d =   xb2d * sinthe + yb2d * costhe
      zb3d = zb1d
      xc3d = - xb2d * costhe - yc2d * sinthe
      yc3d = - xb2d * sinthe + yc2d * costhe
      zc3d = zc1d
      !                                                --- Step5  A3 ---
      
      
      xa3 = trns11*xa3d + trns12*ya3d + trns13*za3d
      ya3 = trns21*xa3d + trns22*ya3d + trns23*za3d
      za3 = trns31*xa3d + trns32*ya3d + trns33*za3d
      xb3 = trns11*xb3d + trns12*yb3d + trns13*zb3d
      yb3 = trns21*xb3d + trns22*yb3d + trns23*zb3d
      zb3 = trns31*xb3d + trns32*yb3d + trns33*zb3d
      xc3 = trns11*xc3d + trns12*yc3d + trns13*zc3d
      yc3 = trns21*xc3d + trns22*yc3d + trns23*zc3d
      zc3 = trns31*xc3d + trns32*yc3d + trns33*zc3d
      
      x1(1,i+ind1)   = xcom + xa3
      x1(2,i+ind1)   = ycom + ya3
      x1(3,i+ind1)   = zcom + za3
      x1(1,i+ind2) = xcom + xb3
      x1(2,i+ind2) = ycom + yb3
      x1(3,i+ind2) = zcom + zb3
      x1(1,i+ind3) = xcom + xc3
      x1(2,i+ind3) = ycom + yc3
      x1(3,i+ind3) = zcom + zc3
      
   end do
   return
end subroutine setlep 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ shake velocities of 3-point waters; unused right now, but please keep!
subroutine setlev (iwatat, natom, x, y, z, vx, vy, vz)
   !Note from Ross Walker : If you decide to use this routine at some point
   !you should invert all the divisions in the same fashion as setlep above.

   !*****************************************************************
   !                                                               **
   !    Subroutine : setlev - reset velocities of TIP3P waters     **
   !    Author : Shuichi Miyamoto                                  **
   !    Date of last update : Dec. 3, 1991                         **
   !                                                               **
   !    Reference for the SETTLE algorithm                         **
   !      S. Miyamoto et al., J. Comp. Chem.,  13, 952 (1992)      **
   !                                                               **
   !*****************************************************************
   implicit none
   integer iwatat, natom
   _REAL_  x, y, z, vx, vy, vz
   dimension  x(*),y(*),z(*),vx(*),vy(*),vz(*)
   !*****************************************************************
   ! Most of the variable names correspond to those in the reference.
   
   !     iwatat : pointer for the first atom of waters
   !     natom  : No. of total atoms in system
   !     x,y,z  : position after  applying constraints (input)
   !     vx,y,z : velocity before applying constraints (output)
   !            : velocity after  applying constraints (output)
   
   !     ab,bc,ca : vectors AB, BC and CA
   !     vabab,.. : eab*vab, ebc*vbc, .
   !     tab,..  = tabd / deno ,  tbcd / deno , .
   
   !     woh,whh : mass of O+H, H+H
   !     wowh2   : O x H x 2
   !     wohwoh  : woh x woh
   
   !*****************************************************************
   _REAL_  wh, wo, wohh, woh, whh, woh2, wowh2, whwh
   _REAL_  wohwoh, wohwo, whhwh
   !     parameter (wo = 16.000d0 , wh = 1.008d0 , wohh = 18.016d0)
   data wo   /16.000d0/ , wh  /1.008d0/ &
         , woh  /17.008d0/ , whh /2.016d0/ &
         , woh2 /34.016d0/
   data wowh2  /32.256d0/     , whwh  /1.016064d0/
   data wohwoh /289.272064d0/ , wohwo /272.128d0/ &
         , whhwh  /2.032128d0/

   _REAL_  ablng, bclng, calng
   _REAL_  abmc, bcma, camb
   _REAL_  cosa, cosb, cosc
   _REAL_  deno
   integer i
   _REAL_  tabd, tbcd, tcad
   _REAL_  vabab, vbcbc, vcaca
   _REAL_  xab, yab, zab, xbc, ybc, zbc, xca, yca, zca
   _REAL_  xeab, yeab, zeab
   _REAL_  xebc, yebc, zebc
   _REAL_  xeca, yeca, zeca
   _REAL_  xvab, yvab, zvab
   _REAL_  xvbc, yvbc, zvbc
   _REAL_  xvca, yvca, zvca
   
   
   do i = iwatat, natom, 3
      !                                                  --- Step1  AB,VAB ---
      
      xab = x(i+1) - x(i)
      yab = y(i+1) - y(i)
      zab = z(i+1) - z(i)
      xbc = x(i+2) - x(i+1)
      ybc = y(i+2) - y(i+1)
      zbc = z(i+2) - z(i+1)
      xca = x(i)   - x(i+2)
      yca = y(i)   - y(i+2)
      zca = z(i)   - z(i+2)
      xvab = vx(i+1) - vx(i)
      yvab = vy(i+1) - vy(i)
      zvab = vz(i+1) - vz(i)
      xvbc = vx(i+2) - vx(i+1)
      yvbc = vy(i+2) - vy(i+1)
      zvbc = vz(i+2) - vz(i+1)
      xvca = vx(i)   - vx(i+2)
      yvca = vy(i)   - vy(i+2)
      zvca = vz(i)   - vz(i+2)
      
      !                                                  --- Step2  eab ---
      
      ablng = sqrt ( xab * xab + yab * yab + zab * zab )
      bclng = sqrt ( xbc * xbc + ybc * ybc + zbc * zbc )
      calng = sqrt ( xca * xca + yca * yca + zca * zca )
      xeab = xab / ablng
      yeab = yab / ablng
      zeab = zab / ablng
      xebc = xbc / bclng
      yebc = ybc / bclng
      zebc = zbc / bclng
      xeca = xca / calng
      yeca = yca / calng
      zeca = zca / calng
      
      !                                                  --- Step3  vabab ---
      
      vabab = xvab * xeab + yvab * yeab + zvab * zeab
      vbcbc = xvbc * xebc + yvbc * yebc + zvbc * zebc
      vcaca = xvca * xeca + yvca * yeca + zvca * zeca
      
      !                                                  --- Step4  tab ---
      
      cosa = - xeab * xeca - yeab * yeca - zeab * zeca
      cosb = - xebc * xeab - yebc * yeab - zebc * zeab
      cosc = - xeca * xebc - yeca * yebc - zeca * zebc
      abmc = wh * cosa * cosb - woh * cosc
      bcma = wo * cosb * cosc - whh * cosa
      camb = wh * cosc * cosa - woh * cosb
      tabd = vabab * ( woh2 - wo * cosc * cosc ) &
            + vbcbc * camb + vcaca * bcma
      tbcd = vbcbc * ( wohwoh - whwh * cosa * cosa ) &
            + vcaca * abmc * wo + vabab * camb * wo
      tcad = vcaca * ( woh2 - wo * cosb * cosb ) &
            + vabab * bcma + vbcbc * abmc
      deno = 2.d0 * wohwoh + wowh2 * cosa * cosb * cosc &
            - whhwh * cosa * cosa &
            - wohwo * ( cosb * cosb + cosc * cosc )
      
      !                                                  --- Step5  V ---
      
      vx(i)   = vx(i)   + ( xeab * tabd - xeca * tcad ) * wh / deno
      vy(i)   = vy(i)   + ( yeab * tabd - yeca * tcad ) * wh / deno
      vz(i)   = vz(i)   + ( zeab * tabd - zeca * tcad ) * wh / deno
      vx(i+1) = vx(i+1) + ( xebc * tbcd - xeab * tabd * wo ) / deno
      vy(i+1) = vy(i+1) + ( yebc * tbcd - yeab * tabd * wo ) / deno
      vz(i+1) = vz(i+1) + ( zebc * tbcd - zeab * tabd * wo ) / deno
      vx(i+2) = vx(i+2) + ( xeca * tcad * wo - xebc * tbcd ) / deno
      vy(i+2) = vy(i+2) + ( yeca * tcad * wo - yebc * tbcd ) / deno
      vz(i+2) = vz(i+2) + ( zeca * tcad * wo - zebc * tbcd ) / deno
   end do   
   return
end subroutine setlev 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine getwds here]
subroutine getwds(igraph    ,nres      ,ipres     ,lbres     , &
      nbonh     ,nbona     ,nbper     ,ib        ,jb        , &
      iwtnm     ,iowtnm    ,ihwtnm    ,jfastw    ,icb       , &
      req       ,winv      ,rbtarg    ,ibelly    ,igrp      , &
      iout)
   
   ! Subroutine GET Water DiStances
   
   ! This routine determines the target lengths for the water molecules and
   ! derives from them the paramters required by fast 3-point water routine
   ! SETLEP.
   
   ! Author: David A. Pearlman
   ! Date: 10/93
   
   ! Input:
   
   ! IGRAPH(I): The name of atom I (I4)
   ! NRES: The number of residues
   ! IPRES(I): IPRES(I)->IPRES(I+1)-1 are the atoms of residue I
   ! LBRES(I): The name of residue I
   ! NBONH: The number of bonds to hydrogens
   ! NBONA: The number of bonds to non-hydrogens
   ! NBPER: The number of bonds in the perturbed group (Gibbs only)
   ! IB(I)
   ! JB(I): The two atoms of bond I. Stored as 3*(atom_#-1)
   ! IWTNM: The name of residues to be considered water (I)
   ! IOWTNM: The name of oxygen atom in the water residues (I)
   ! IHWTNM(2): The name of 2 hydrogen atoms in the water residues (I)
   ! JFASTW: = 0: Use fast water routine (possibly with modified names)
   !           4: Do not use fast waters anywhere.
   ! ICB(I): Parameter pointer for bond number I.
   ! REQ(I): REQ(ICB(I)) is the target bond length of bond I.
   ! WINV(I): The inverse mass of atom I.
   ! IOUT: Unit for warning writes.
   
   ! Output:
   
   ! RBTARG(9): ra, rb, rc, rc2, hhhh mass(O), mass(H), one/mass(HOH) and one/ra
   ! as required by the SETLEP routine.
   
   
   
   implicit none
   character(len=4) igraph, lbres, iwtnm, iowtnm, ihwtnm
   integer nres      ,ipres
   integer nbonh     ,nbona     ,nbper     ,ib        ,jb
   integer jfastw    ,icb
   integer ibelly, igrp, iout
   _REAL_  req       ,winv      ,rbtarg
   dimension igraph(*),ipres(*),lbres(*),ib(*),jb(*)
   dimension ihwtnm(*),icb(*),req(*),rbtarg(9),winv(*)
   dimension igrp(*)
   
   integer i, j
   integer iat, ihf1, ihf2
   integer iof
   integer nofast
   _REAL_  roh1, roh2, rhh, roh

   _REAL_  small
   save small
   data small/1.0d-4/
   
   if (jfastw == 4) return
   
   ! Search through the list for a water molecule. If one is found,
   ! determine the O-H and H-H distances ascribed to it.
   
   nofast = 0
   nres_loop: do i = 1,nres
      iof = -1   ! -1 is not a valid index into IB and JB
      ihf1 = -1
      ihf2 = -1
      if (lbres(i) /= iwtnm) cycle nres_loop
      iof = 0
      ihf1 = 0
      ihf2 = 0
      do  j = 1,3
         iat = ipres(i)+j-1
         
         !         -- any waters not part of moving belly (IGRP(I)=0)
         !            will have had their bond parameters removed from
         !            list already and these must be skipped here
         
         if (ibelly > 0 .and. igrp(iat) <= 0) cycle nres_loop
         
         if (igraph(iat) == iowtnm) then
            iof = 3*(iat-1)
            rbtarg(6) = 1.0d0/winv(iat)
         else if (igraph(iat) == ihwtnm(1)) then
            ihf1 = 3*(iat-1)
            rbtarg(7) = 1.0d0/winv(iat)
         else if (igraph(iat) == ihwtnm(2)) then
            ihf2 = 3*(iat-1)
            rbtarg(7) = 1.0d0/winv(iat)
         else
            cycle nres_loop
         end if
      end do
      if (iof > 0 .and. ihf1 > 0 .and. ihf2 > 0) goto 50
   end do nres_loop
   nofast = 1
   
   ! We have found the three atoms. Now search the bonds list for the
   ! corresponding bonds. At that point, we can determine the target distances
   ! from REQ.
   
   50 roh1 = -10.0d0
   roh2 = -10.0d0
   rhh  = -10.0d0
   do  j = 1,nbonh+nbona+nbper
      if (ib(j) == iof .and. jb(j) == ihf1) then
         roh1 = req(icb(j))
      else if (jb(j) == iof .and. ib(j) == ihf1) then
         roh1 = req(icb(j))
      else if (ib(j) == iof .and. jb(j) == ihf2) then
         roh2 = req(icb(j))
      else if (jb(j) == iof .and. ib(j) == ihf2) then
         roh2 = req(icb(j))
      else if (ib(j) == ihf1.and. jb(j) == ihf2) then
         rhh  = req(icb(j))
      else if (jb(j) == ihf1.and. ib(j) == ihf2) then
         rhh  = req(icb(j))
      end if
   end do
   
   ! If all three bond lengths were not assigned in list, or if ROH1 and
   ! ROH2 are not the same length, assume something is a bit cockeye. In
   ! this case, issue a warning and use the default TIP3P values.
   
   if (roh1 < 0.0d0 .or. roh2 < 0.0d0 .or. rhh < 0.0d0) then
      if (nofast == 0) write(iout,1001)
      rhh = 1.5136d0
      roh = 0.9572d0
      rbtarg(6) = 16.0000d0
      rbtarg(7) = 1.008d0
   else if (abs(roh1-roh2) > small) then
      write(iout,1002)
      rhh = 1.5136d0
      roh = 0.9572d0
      rbtarg(6) = 16.0000d0
      rbtarg(7) = 1.008d0
   else
      roh = roh1
   end if
   rbtarg(8) = 1.0d0/(rbtarg(6)+2.0d0*rbtarg(7))
   
   ! B4SETL sets the appropriate constants based on the bond lengths
   
   call b4setl(rhh,roh,rbtarg)
   
   return
   
   ! WARNINGS:
   
   1001 format('WARNING: Bond lengths params not found in PARM file', &
         ' for fast water model;',/,t10,'using TIP3P defaults')
   1002 format('WARNING: R(O-H) bond lengths params found in PARM file', &
         ' for fast water model',/,t10, &
         'not identical; using TIP3P defaults')
end subroutine getwds 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine b4setl here]
subroutine b4setl(rhh,roh,rbtarg)
  
   ! Subroutine Before SETLep
   
   ! This routine determines the values of ra, rb, rc, rc2, and hhhh
   ! required by the SETLEP program. These are calculated from the
   ! target rigid water distances RHH (H...H) and ROH (O...H).
   
   ! The calulated values are returned in RBTARG(1-5):
   
   !    RBTARG(1) ... ra
   !    RBTARG(2) ... rb
   !    RBTARG(3) ... rc
   !    RBTARG(4) ... rc2
   !    RBTARG(5) ... hhhh
   !    (RBTARG(6->8) are masses already set in calling routine)
   !    RBTARG(9) ... one/ra (Ross Walker - Mod for speed - avoid excessive divisions)
   
   ! Author: David A. Pearlman
   ! Date: 10/93
   use constants, only : zero, half
   implicit none
   _REAL_  rhh,roh,rbtarg(9), x(2,3)

   _REAL_  wh, wo, onewohh
   !     PARAMETER (WO = 16.000D0 , WH = 1.008D0 , WOHH = 18.016D0)
   !     PARAMETER (RHH=1.5136D0 , ROH=0.9572D0)

   _REAL_  comx, comy
   _REAL_  dis
   _REAL_  height
   
   wo = rbtarg(6)
   wh = rbtarg(7)
   onewohh = rbtarg(8)
   
   ! create triangle in x,y plane with base parallel to x.
   
   height = sqrt(roh**2-(rhh*half)**2)
   x(1,1) = -rhh*half
   x(2,1) = -height
   x(1,2) = +rhh*half
   x(2,2) = -height
   x(1,3) = zero
   x(2,3) = zero
   
   ! calculate the center of mass of the triangle:
   
   comx = (x(1,1)*wh + x(1,2)*wh + x(1,3)*wo)*onewohh
   comy = (x(2,1)*wh + x(2,2)*wh + x(2,3)*wo)*onewohh
   
   ! the distance between the center of mass and the apex is ra. rb is
   ! the height - ra.
   
   dis = sqrt(comx*comx+comy*comy)
   rbtarg(1) = dis
   rbtarg(2) = height-rbtarg(1)
   rbtarg(3) = rhh*half
   rbtarg(4) = rhh
   rbtarg(5) = rbtarg(4)*rbtarg(4)
   rbtarg(9) = 1.0d0/dis
   
   !     print*,'ra=  ',rbtarg(1)
   !     print*,'rb=  ',rbtarg(2)
   !     print*,'rc=  ',rbtarg(3)
   !     print*,'rc2= ',rbtarg(4)
   !     print*,'hhhh=',rbtarg(5)

   return
end subroutine b4setl 

end module fastwt
