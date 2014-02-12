! <compile=optimized>
#include "copyright.h"
#define _REAL_ double precision
#include "pb_def.h"
!#ifndef SANDER
#include "timer.h"
!#endif

module dispersion_cavity

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Chuck Tan, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   implicit none

   logical donpsa
   integer npopt
   integer decompopt
   integer use_rmin
   integer use_sav
   _REAL_ rhow_effect
   _REAL_ cavity_surften
   _REAL_ cavity_offset
 
contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver of nonelectrostatic solvation energy and forces
subroutine np_force( natom,nres,ntypes,ipres,iac,ico,cn1,cn2,x,f,enbrfcav,enbrfdis )

   use poisson_boltzmann, only : cutnb, acrd, iprshrt, iar1pb, nex, iex, &
                                 ligand, outflag
   use solvent_accessibility
!#ifndef SANDER
   use pbtimer_module
!#endif

   ! Common variables
    
#  include "pb_md.h"
#ifdef SANDER
#  undef _REAL_
#  include "md.h"
#  include "box.h"
#else
#  include "md.h"
#  include "box.h"
#endif
#  include "pb_constants.h"
#  include "flocntrl.h"
#define _REAL_ double precision

   ! Passed variables
    
   integer natom, nres, ntypes, ipres(*), iac(*), ico(*)
   _REAL_ cn1(*), cn2(*), x(3,natom), f(3,natom)
   _REAL_ enbrfcav, enbrfdis
    
   ! Local variables
    
   integer ic, iatm
!   _REAL_ area(natom), darea(3,natom)
!   _REAL_, allocatable:: xarea(:,:), yarea(:,:), zarea(:,:)
   _REAL_ mdsig_ow, mdsig_iatm, epsln_iatm, sigow2, sigow4, sigow6
   _REAL_ a(natom), b(natom), epsow(natom), sigow(natom), rminow(natom)

   if ( do_pbnp == 0 ) return

   enbrfcav= ZERO; enbrfdis = ZERO
    
   ! compute sa surface and arcs for dispersion energy and forces
 
!#ifndef SANDER
   call pbtimer_start(PBTIME_NPSAS)
!#endif
   if ( donpsa ) then
      if ( use_rmin == 0 ) call sa_init(pbverbose,pbprint,natom,natom,ifcap,sprob,mdsig,radip,radip2,outflag)
      if ( use_rmin == 1 ) call sa_init(pbverbose,pbprint,natom,natom,ifcap,sprob,rmin,radip,radip2,outflag)
      call sa_driver(pbverbose,pbprint,ipb,inp,natom,natom,dosas,ndosas,npbstep,nsaslag,&
              ligand, outflag,&
              acrd(1,1),iar1pb(1,0),iprshrt,nex,iex,.true.)
   end if
!#ifndef SANDER
   call pbtimer_stop(PBTIME_NPSAS)
!#endif

   ! compute cavity solvation energy and forces

!#ifndef SANDER
   call pbtimer_start(PBTIME_NPCAV)
!#endif
!   allocate(xarea(natom,natom))
!   allocate(yarea(natom,natom))
!   allocate(zarea(natom,natom))
!   call np_cavity(natom,cavity_surften,radip,x,enbrfcav,area,darea,xarea,yarea,zarea) 
   if ( use_sav == 0 ) enbrfcav = cavity_surften * prtsas + cavity_offset
   if ( use_sav == 1 ) enbrfcav = cavity_surften * prtsav + cavity_offset
   if ( pbprint ) then
      if ( pbverbose ) write(6, '(1x,a,2f12.4)') 'Nonpolar SAS/SAV', prtsas, prtsav
      write(6, '(1x,a,f12.4)') 'Cavity solvation energy', enbrfcav
   end if
!#ifndef SANDER
   call pbtimer_stop(PBTIME_NPCAV)
!#endif
 
   ! compute dispersion solvation energy forces
   ! first obtain van der Waals A and B parameters between iatm and TIP3P/OW

!#ifndef SANDER
   call pbtimer_start(PBTIME_NPDIS)
!#endif
   if ( npopt == 2 ) then
      mdsig_ow = 1.7683d0*(TWO**(-SIXTH))
      do iatm = 1, natom
         ic = ico(ntypes*(iac(iatm)-1) + iac(iatm))
         if (cn2(ic) /= ZERO) then
            mdsig_iatm = (cn1(ic)/cn2(ic))**SIXTH/2
            epsln_iatm = cn2(ic)/(256.0d0*mdsig(iatm)**6)   
         else
            mdsig_iatm = ZERO
            epsln_iatm = ZERO
         endif
         sigow(iatm) = mdsig_iatm + mdsig_ow
         rminow(iatm) = sigow(iatm)*(TWO**SIXTH)
         epsow(iatm) = sqrt(epsln_iatm*0.1520d0)
         sigow2 = sigow(iatm)*sigow(iatm)
         sigow4 = sigow2*sigow2
         sigow6 = sigow2*sigow4
         b(iatm) = FOUR*epsow(iatm)*sigow6
         a(iatm) = b(iatm)*sigow6
      enddo
   
      call np_dispersion( )
      if ( pbprint ) then
         write(6, '(1x,a,f12.4)') 'Dispersion solvation energy', enbrfdis
      end if
   end if
!#ifndef SANDER
   call pbtimer_stop(PBTIME_NPDIS)
!#endif

!#ifndef SANDER
   call pbtimer_start(PBTIME_NPSAS)
!#endif
   if ( donpsa ) then
      call sa_free( dosas, ndosas, .true. )
   end if
!#ifndef SANDER
   call pbtimer_stop(PBTIME_NPSAS)
!#endif

   ! returning:

contains
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dispersion solvation energy and force
subroutine np_dispersion( )

#  include "pb_constants.h"

   ! Passed variables
    
   ! integer natom              ! no of atoms
   ! integer nres               ! no of residues
   ! integer ipres(*)           ! first atom of a residue
   ! integer maxsph             ! max no of surface dots per atom
   ! integer narc               ! no of arcs per molecule
   ! integer nsrfdot            ! no of surface _dots_ per molecule
   ! integer narcdot            ! no of arc _dots_ per molecule
   ! integer fstsdot(natom)     ! the first surface dot of an atom
   ! integer lstsdot(natom)     ! the last surface dot of an atom
   ! integer fstadot(natom)     ! the first arc dot of an atom
   ! integer lstadot(natom)     ! the last arc dot of an atom
   ! integer arcatm(2,*)        ! the two atoms of a given arc
   ! integer dotarc(*)          ! the arc that generates the current arc dot
   ! _REAL_ cutnb             ! cutoff for vdw interactions
   ! _REAL_ savactr(3,*)      ! coordinates of the centers of arcs
   ! _REAL_ srfcrd(3,*)       ! coordinates of all surface dots
   ! _REAL_ arccrd(3,*)       ! coordinates of all arc dots
   ! _REAL_ a(natom), b(natom)! atomic vdw coefficients
   ! _REAL_ radip(natom)      ! atomic vdw radii
   ! _REAL_ x(3,natom)        ! atomic coordinates
   ! _REAL_ f(3,natom)        ! returned solvation dispersion forces
   ! _REAL_ enbrfdis          ! returned solvation dispersion energy

   ! Local variables

   integer iatm, jatm, katm
   ! integer iarc, idot
   integer jdot
   integer i, j, ii, ip1, ip2
   ! integer firstdot, lastarc
   ! integer narcik
   integer nbindex
   integer nblist(natom)

   _REAL_ dx, dy, dz, d2
   _REAL_ ris(3), risnorm(3), ris1, ris2, ris4, ris_1
   _REAL_ snorm(3), r_1, dotprot
   ! _REAL_ rls(3), rls2, rls4
   _REAL_ adis
   _REAL_ cutoff, cutoff4, cutlng, rcut2    
   ! _REAL_ rarc_1, rarcx, rarcy, rarcz
   ! _REAL_ avegamma, tmpgamma
   _REAL_ costheta, cross_sect
   _REAL_ tmpg, tmpf, gcorrec
   _REAL_ tmpfx, tmpfy, tmpfz !, fxcorrec, fycorrec, fzcorrec
   _REAL_ delts(natom)
   _REAL_ nbd, nbd2(natom), nbdx(natom), nbdy(natom), nbdz(natom)

   _REAL_, parameter :: rhow = 3.3330d-02

   ! for InsightII display

   !print *, 'Number of remaining surface points', nsrfdot
   !open (unit=55, file='sasrf.dot')
   !write (55, '("DOTS")')
   !do jatm = 1, natom
   !   do idot = fstsdot(jatm), lstsdot(jatm)
   !      write (55,'(4(f8.3,2x))') srfcrd(1,idot)+x(1,jatm), srfcrd(2,idot)+x(2,jatm), srfcrd(3,idot)+x(3,jatm), 300.
   !   enddo
   !enddo
   !close(55)
 
   adis = ZERO

   if (cutnb == ZERO) then
      cutoff = 999.9d0
   else
      cutoff = sqrt(cutnb)
   endif
   cutlng = (cutoff + 6.0d0)*(cutoff + 6.0d0)
   cutoff4 = (cutoff*cutoff)**2

   ! convert b's to those for bulk water with a prefactor of 1/3 from the surface integration.

   b = rhow_effect*rhow*THIRD*b
   a = rhow_effect*rhow*THIRD*THIRD*a
   epsow = rhow_effect*rhow*THIRD*epsow

   ! detlaS per SAS dot

   delts = FOURPI*(radip**2)/maxsph

   do iatm = 1, natom

      if ( radip(iatm) == ZERO ) cycle

      !
      ! find iatm's neighbors
      !

      nbindex = 0; nbd2(1:natom) = ZERO; nblist(1:natom) = 0
      nbdx(1:natom) = ZERO; nbdy(1:natom) = ZERO; nbdz(1:natom) = ZERO

      do ii = 1, nres

         ip1 = ipres(ii)
         ip2 = ipres(ii+1)-1

         dx  = (x(1,ip1) + x(1,ip2))/TWO - x(1,iatm)
         dy  = (x(2,ip1) + x(2,ip2))/TWO - x(2,iatm)
         dz  = (x(3,ip1) + x(3,ip2))/TWO - x(3,iatm)
         d2  = dx*dx + dy*dy +dz*dz
         if ( d2 > cutlng ) cycle

         do j = ip1, ip2

            if (radip(j) == ZERO ) cycle

            dx  = x(1,j) - x(1,iatm)
            dy  = x(2,j) - x(2,iatm)
            dz  = x(3,j) - x(3,iatm)
            d2 = dx*dx + dy*dy + dz*dz
            rcut2  = (radip(j) + cutoff)*(radip(j) + cutoff)

            if ( d2 < rcut2 ) then
               nbindex = nbindex + 1
               nbd2(nbindex) = d2
               nblist(nbindex) = j
               nbdx(nbindex) = dx; nbdy(nbindex) = dy; nbdz(nbindex) = dz 
            endif 

         enddo ! do j = ip1, ip2

      enddo ! do ii = 1, nres

      !
      ! ----- Gamma Part -----
      ! For every pair of atom (iatm) and tessera (on jatm), calculate its
      ! contribution to the total dispersion energy and to the force on iatm and jatm.
      !

      do j = 1, nbindex

         jatm  = nblist(j)

         tmpg = ZERO
         gcorrec = ZERO
         tmpfx = ZERO !; fxcorrec = ZERO
         tmpfy = ZERO !; fycorrec = ZERO
         tmpfz = ZERO !; fzcorrec = ZERO
         r_1   = ONE/radip(jatm)
         nbd = sqrt(nbd2(j))

         ! if jatm is exposed

         if ( fstsdot(jatm) <= lstsdot(jatm) ) then
            do jdot = fstsdot(jatm), lstsdot(jatm)
               ris(1) = nbdx(j) + srfcrd(1,jdot)
               ris(2) = nbdy(j) + srfcrd(2,jdot)
               ris(3) = nbdz(j) + srfcrd(3,jdot)
               snorm(1) = srfcrd(1,jdot)*r_1
               snorm(2) = srfcrd(2,jdot)*r_1
               snorm(3) = srfcrd(3,jdot)*r_1
               ris2 = ris(1)*ris(1)+ris(2)*ris(2)+ris(3)*ris(3)
               ris4 = ris2*ris2
               ris1 = sqrt(ris2)
               ris_1 = ONE/ris1
               dotprot = ris(1)*snorm(1)+ris(2)*snorm(2)+ris(3)*snorm(3)

               ! dispersion energy

               ! 6/12 decomposition
                
               if ( decompopt == 1 ) then
                
                  adis = - b(iatm)/(ris2*ris4)
                
               ! sigma decomposition
                
               else if  ( decompopt == 2 ) then
                   
                  if ( ris1 > sigow(iatm) ) then
                     adis = - b(iatm)/( ris2*ris4     ) &
                            + a(iatm)/( ris4*ris4*ris4)
                  else
                     adis = - b(iatm)/((sigow(iatm)   *ris1 )**3) &
                            + a(iatm)/( sigow(iatm)**3*ris1 )**3
                  endif
                
               ! WCA decomposition
                
               else if ( decompopt == 3 ) then
                  if ( ris1 > rminow(iatm) ) then
                     adis = - b(iatm)/(ris2*ris4) &
                            + a(iatm)/(ris4*ris4*ris4)
                  else
                     adis = epsow(iatm)*(ONE - (rminow(iatm)/ris1)**3) &
                         - b(iatm)/((rminow(iatm)*ris1)**3) &
                         + a(iatm)/(rminow(iatm)**3*ris1)**3              
                  endif
                
               end if
                
               tmpg = tmpg + adis*dotprot

               ! derivative of dispersion energy

               risnorm(1) = ris(1)*ris_1
               risnorm(2) = ris(2)*ris_1
               risnorm(3) = ris(3)*ris_1
               tmpf = TWO*b(iatm)/(ris4*ris2*ris1)*dotprot
               tmpfx = tmpfx - adis*snorm(1) - tmpf*risnorm(1)
               tmpfy = tmpfy - adis*snorm(2) - tmpf*risnorm(2)
               tmpfz = tmpfz - adis*snorm(3) - tmpf*risnorm(3)

            enddo ! do jdot = fstsdot(jatm), lstsdot(jatm) 

         ! if the pair of atoms are too far away from each other
         ! add energy correction. However, if jatm is exposed, do not
         ! add correction.
         ! force correction can be ignored

         elseif ( nbd > (cutoff - radip(jatm)) ) then

            costheta = ( nbd2(j) + cutnb - radip(jatm)*radip(jatm) )/( TWO*nbd*cutoff )
            cross_sect = TWOPI*cutoff*(ONE - costheta) ! cross-section area
            if ( decompopt == 1 ) then
               adis = -b(iatm)/(cutoff*cutoff4)
            else if ( decompopt == 2 .or. decompopt == 3 ) then
               adis = -b(iatm)/(cutoff*cutoff4) + a(iatm)*cutoff/(cutoff4*cutoff4*cutoff4)
            end if
            gcorrec = gcorrec + adis*cross_sect

!           tmpf  = -FIVE*b(iatm)/(cutoff4*cutnb)
!           fxcorrec  = fxcorrec - tmpf*(nbdx(j)/nbd)*cross_sect
!           fycorrec  = fycorrec - tmpf*(nbdy(j)/nbd)*cross_sect
!           fzcorrec  = fzcorrec - tmpf*(nbdz(j)/nbd)*cross_sect

         endif
                
         enbrfdis = enbrfdis + tmpg*delts(jatm) + gcorrec
         f(1, iatm) = f(1, iatm) + tmpfx*delts(jatm) !+ fxcorrec
         f(2, iatm) = f(2, iatm) + tmpfy*delts(jatm) !+ fycorrec
         f(3, iatm) = f(3, iatm) + tmpfz*delts(jatm) !+ fzcorrec
         f(1, jatm) = f(1, jatm) - tmpfx*delts(jatm) !- fxcorrec
         f(2, jatm) = f(2, jatm) - tmpfy*delts(jatm) !- fycorrec
         f(3, jatm) = f(3, jatm) - tmpfz*delts(jatm) !- fzcorrec

      enddo ! do j = 1, nbindex

!      !
!      ! ----- S Part -----
!      ! For every arc, which is indexed by iatm, calculate its contributions to the
!      ! two atoms that generate it. To do this, the average of integrand over all the arcdot
!      ! on the arc is calculated.
!      !
! 
!      firstdot = fstadot(iatm)
!      lastarc = dotarc(firstdot)
!      narcik = 0
!      avegamma = ZERO
!
!      do idot = fstadot(iatm), lstadot(iatm)
!
!         ! Retrieve the arc that generates the dot
!
!         iarc = dotarc(idot)
!          
!         ! Get the norm of this dot
!
!         rarcx = arccrd(1, idot) - savactr(1, iarc)
!         rarcy = arccrd(2, idot) - savactr(2, iarc)
!         rarcz = arccrd(3, idot) - savactr(3, iarc)
!         rarc_1 = ONE/sqrt(rarcx*rarcx+rarcy*rarcy+rarcz*rarcz)
!         snorm(1) = rarcx*rarc_1
!         snorm(2) = rarcy*rarc_1
!         snorm(3) = rarcz*rarc_1
!
!         ! Compute Gamma over iatm's neighbors
!         ! tmpgamma is for this dot only
!         ! Gamma is computed with cutoff and with no correction
!         ! because the correction is very small just like above
!
!         tmpgamma = ZERO
! 
!         do j = 1, nbindex
!
!            jatm = nblist(j)
!
!            if ( radip(jatm) == ZERO  ) cycle
!
!            rls(1) = arccrd(1,idot) - x(1,jatm)
!            rls(2) = arccrd(2,idot) - x(2,jatm)
!            rls(3) = arccrd(3,idot) - x(3,jatm)
!            rls2 = rls(1)*rls(1) + rls(2)*rls(2) + rls(3)*rls(3)
!            rls4 = rls2*rls2
!            dotprot = rls(1)*snorm(1) + rls(2)*snorm(2) + rls(3)*snorm(3)
!            adis = -b(jatm)/(rls2*rls4)
!            tmpgamma = tmpgamma + adis*dotprot
!
!         enddo ! do j = 1, nbindex
!
!         ! Are we still on the same arc .AND. NOT the last arcdot of this atom?
!
!         if ( (iarc == lastarc) .and. (idot /= lstadot(iatm)) ) then
! 
!            ! We are on the same arc and not the last arcdot of iatm.
!
!            narcik = narcik + 1
!            avegamma = avegamma + tmpgamma  
!
!         elseif ( iarc /= lastarc ) then
!
!            ! We are on the next arc, so accumulate the force for the previous pair of atoms.
!            ! retrieve katm that generates the arc with iatm
!              
!            if ( iatm == arcatm(1,lastarc) ) then
!               katm = arcatm(2,lastarc)
!            else
!               katm = arcatm(1,lastarc)
!            endif
!
!            ! Convert to forces and accumulate
!
!            avegamma = avegamma/narcik
!            tmpfx = xarea(iatm, katm)*avegamma
!            tmpfy = yarea(iatm, katm)*avegamma
!            tmpfz = zarea(iatm, katm)*avegamma
!            f(1,iatm) = f(1,iatm) + tmpfx
!            f(2,iatm) = f(2,iatm) + tmpfy
!            f(3,iatm) = f(3,iatm) + tmpfz
!            f(1,katm) = f(1,katm) - tmpfx
!            f(2,katm) = f(2,katm) - tmpfy
!            f(3,katm) = f(3,katm) - tmpfz
!
!            ! Reinitialization with the contribution of this dot to the new arc.
!
!            lastarc = iarc; narcik = 1; avegamma = tmpgamma
!
!         else
!
!            ! We're now on the last arcdot of iatm, so do the force accumulation
!            ! for the previous pair of atoms. Accumulate this dot's contribution first.
!
!            narcik = narcik + 1
!            avegamma = avegamma + tmpgamma 
!
!            ! Retrieve katm that generates the arc
!
!            if ( iatm == arcatm(1, lastarc) ) then
!               katm = arcatm(2, lastarc)
!            else
!               katm = arcatm(1, lastarc)
!            endif
!
!            ! Convert to forces and accumulate
!
!            avegamma = avegamma/narcik
!            tmpfx = xarea(iatm, katm)*avegamma
!            tmpfy = yarea(iatm, katm)*avegamma
!            tmpfz = zarea(iatm, katm)*avegamma
!            f(1, iatm) = f(1, iatm) + tmpfx
!            f(2, iatm) = f(2, iatm) + tmpfy
!            f(3, iatm) = f(3, iatm) + tmpfz
!            f(1, katm) = f(1, katm) - tmpfx
!            f(2, katm) = f(2, katm) - tmpfy
!            f(3, katm) = f(3, katm) - tmpfz
!
!         endif
!
!      enddo ! do idot = fstadot(iatm), lstadot(iatm)

   enddo ! do iatm = 1, natom                

end subroutine np_dispersion

end subroutine np_force

end module dispersion_cavity
