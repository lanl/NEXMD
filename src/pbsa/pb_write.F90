#include "copyright.h"
#define _REAL_ double precision

!  Author Meng-Juei Hsieh
subroutine croak(strwhich, iwhere)
   implicit none
   character(*) strwhich
   integer      iwhere

   write(6,'(a,1x,a,1x,a,i6)') 'Error at ',strwhich,', line ',iwhere
   call mexit(6,1)
end subroutine croak

#ifndef SANDER
! This is a refactorized version of prntmd
! Energy output for md, in human-readable form.
subroutine prntmd(nstep,nitp,nits,time,ener,fac,iout7,rms)
   implicit none
#  include "files.h"
#  include "md.h"

   integer nstep, nitp, nits, iout7
   _REAL_  time, ener, fac
   dimension ener(*),fac(*)
   logical rms

   _REAL_  etot,   ektot,   rms_pbs
   _REAL_  volume, densit 
   _REAL_  press
   _REAL_  ekcmt,   virt
   _REAL_  epot,   enonb,   eel,    ehbond, ebond,  eangle 
   _REAL_  edihed, enb14,   eel14,  econst, epol
   _REAL_  esurf,  e3bod,  dvdl,   edisp

   if ( imin == 1 .or. imin == 6 ) then ! reused to print minimization results
      etot   = ener(1)
      enonb  = ener(2)
      eel    = ener(3)
      ehbond = ener(4)
      ebond  = ener(5)
      eangle = ener(6)
      edihed = ener(7)
      enb14  = ener(8)
      eel14  = ener(9)
      econst = ener(10)
      epol   = ener(11)
      !aveper = ener(12)
      !aveind = ener(13)
      !avetot = ener(14)
      esurf  = ener(15)
      e3bod  = ener(16)
      edisp  = ener(18)
      !diprms = ener(22)
      !dipiter = ener(23)
      !dipole_temp = ener(24)
      dvdl   = 0d0
      volume = 0d0
      epot   = 0d0
      ektot  = 0d0
      temp   = 0d0
      press  = 0d0
   else
      etot   = ener(1)
      ektot  = ener(2)
      temp   = ektot/fac(1)
      !eksolt = ener(3)/fac(2)
      
      if(ntt == 5) then
         !eksolv = ener(4)/fac(3)
      else
         rms_pbs = ener(4)
      end if
      
      !scaltp = ener(5)
      
      !boxx   = ener(7)
      !boxy   = ener(8)
      !boxz   = ener(9)
      volume = ener(10)
      densit = ener(42)
      
      !presx  = ener(11)
      !presy  = ener(12)
      !presz  = ener(13)
      press  = ener(14)
      
      !ekcmx  = ener(15)
      !ekcmy  = ener(16)
      !ekcmz  = ener(17)
      ekcmt  = ener(18)
      
      !virx   = ener(19)
      !viry   = ener(20)
      !virz   = ener(21)
      virt   = ener(22)
      
      epot   = ener(23)
      enonb  = ener(24)
      eel    = ener(25)
      ehbond = ener(26)
      ebond  = ener(27)
      eangle = ener(28)
      edihed = ener(29)
      enb14  = ener(30)
      eel14  = ener(31)
      econst = ener(32)
      epol   = ener(33)
      !aveper   = ener(34)
      !aveind   = ener(35)
      !avetot   = ener(36)
      esurf = ener(37)
      e3bod = ener(38)
      dvdl = ener(39)
      edisp = ener(40)
      !virvsene = ener(43)
      !diprms = ener(44)
      !dipiter = ener(45)
      !dipole_temp = ener(46)
   endif
   
   write(6,9018) nstep,time,temp,press
   write(6,9028) etot,ektot,epot
   write(6,9038) ebond,eangle,edihed
   write(6,9048) enb14,eel14,enonb
   if (ipb < 1) then
      write(6,9059) eel,ehbond,econst
   else
      write(6,9060) eel,ehbond,econst
   endif
   if ( gbsa > 0         ) write(6,9077) esurf
   if ( ipb >= 1         ) write(6,9074) esurf,edisp
   if ( econst /= 0.0    ) write(6,9076) epot-econst
   if ( dvdl /= 0.d0     ) write(6,9089) dvdl
   if ( rms .and. ntt==0 ) write(6,9075) rms_pbs
   if ( volume /= 0.0    ) write(6,9078) ekcmt,virt,volume
   if ( epol /= 0.0 .or. e3bod /= 0.0 ) &
        write(6,9070) epol,e3bod
   if ( volume /= 0.0    ) write(6,9079) densit

   write(6,8088)
   
   !     --- flush i/o buffer ---
   
   call amflsh(6)
   if (iout7 == 0) return
   
   !       ----- OUTPUT THE INFO FILE if requested -----
   
   write(7,9018) nstep,time,temp,press
   write(7,9028) etot,ektot,epot
   write(7,9038) ebond,eangle,edihed
   write(7,9048) enb14,eel14,enonb
   write(7,9059) eel,ehbond,econst
   if (gbsa > 0) write(7,9077) esurf
   if (ipb >= 1) write(7,9074) esurf,edisp
   if (econst /= 0.0) write(7,9076) epot-econst
   if ( dvdl /= 0.d0) write(7,9089) dvdl
   if (rms .and. ntt==0 ) write(6,9075) rms_pbs
   if (volume /= 0.0) write(7,9078) ekcmt,virt,volume
   if (epol /= 0.0 .or. e3bod /= 0.0) &
         write(7,9070) epol,e3bod
   if (volume /= 0.0) write(7,9079) densit
   return
   
   8088 format(t2,78('-'),/)
   9018 format(/1x, 'NSTEP =',i9,3x,'TIME(PS) =',f12.3,2x, &
         'TEMP(K) =',f9.2,2x,'PRESS =',f8.1)
   9028 format (1x,'Etot   = ',f14.4,2x,'EKtot   = ',f14.4,2x, &
         'EPtot      = ',f14.4)
   9038 format (1x,'BOND   = ',f14.4,2x,'ANGLE   = ',f14.4,2x, &
         'DIHED      = ',f14.4)
   9048 format (1x,'1-4 NB = ',f14.4,2x,'1-4 EEL = ',f14.4,2x, &
         'VDWAALS    = ',f14.4)
   9059 format (1x,'EELEC  = ',f14.4,2x,'EGB     = ',f14.4,2x, &
         'RESTRAINT  = ',f14.4)
   9060 format (1x,'EELEC  = ',f14.4,2x,'EPB     = ',f14.4,2x, &
         'RESTRAINT  = ',f14.4)
   9075 format ('|E(PBS) = ',f14.4)
   9076 format (1x,'EAMBER (non-restraint)  = ',f14.4)
   9077 format (1x,'ESURF= ',f14.4)
   9074 format (1x,'ECAVITY= ',f14.4,2x,'EDISPER = ',f14.4)
   9078 format (1x,'EKCMT  = ',f14.4,2x,'VIRIAL  = ',f14.4,2x, &
         'VOLUME     = ',f14.4)
   9079 format (52x,'Density    = ',f14.4)
   
   9070 format (1x,'EPOLZ  = ',f14.4,2x,'E3BODY  = ',f14.4)
   
   9089 format (1x,'DV/DL  = ',f14.4)

end subroutine prntmd 
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  Printing the final minimization report 
!  It becomes clear that many of the arguments in this subroutine are
!  now dummies, we still preserve the interface for future usage.
subroutine report_min_results( nstep, gradient_rms, coordinates, &
      forces, energies, igraph, xx, ix, ih )

   use decomp, only: checkdec, printdec, printpdec
   implicit none

   integer, intent(in)    :: nstep
   _REAL_,  intent(in)    :: gradient_rms
   _REAL_,  intent(in)    :: coordinates(*)
   _REAL_,  intent(in)    :: forces(*)
   _REAL_,  intent(in)    :: energies(51)
   character(len=4), intent(in)    :: igraph(*) ! atom name map
   _REAL_,  intent(inout) :: xx(*)              ! real dynamic memory
   integer, intent(inout) :: ix(*)              ! integer dynamic memory
   character(len=4), intent(inout) :: ih(*)     ! hollerith dynamic memory

#  include "extra.h"
#  include "md.h"

   _REAL_  :: dummyfac(3)

   dummyfac=1

   if (master) then
      write(6, '(/ /20x,a,/)' ) 'FINAL RESULTS'
      call prntmd( nstep, 0, 0, 0d0, energies, dummyfac, 0, .false. )
      if (idecomp > 0) call checkdec(idecomp)
      if (idecomp == 1 .or. idecomp == 2) call printdec(ix)
      if (idecomp == 3 .or. idecomp == 4) call printpdec(ix)
   end if

   return
end subroutine report_min_results


!  Print out a "minimization progress" report 
!  It becomes clear that many of the arguments in this subroutine are
!  now dummies, we still preserve the interface for future usage.
subroutine report_min_progress( nstep, gradient_rms, forces, energies, igraph )
   implicit none
   integer, intent(in)    :: nstep
   _REAL_,  intent(in)    :: gradient_rms
   _REAL_,  intent(in)    :: forces(*)
   _REAL_,  intent(in)    :: energies(51)
   character(len=4), intent(in)    :: igraph(*)    ! atom name map

#  include "extra.h"

   _REAL_  :: dummyfac(3)

   dummyfac=1

   if (master) then
      call prntmd( nstep, 0, 0, 0d0, energies, dummyfac, 0, .false. )
   end if

   return
end subroutine report_min_progress
#endif /*ifndef SANDER*/

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ print out information on current fdpb calculation
subroutine pb_print ( ifcap, ipb, natom ) 
   
   ! Module variables
   
   use poisson_boltzmann
   use solvent_accessibility
   implicit none
   
   ! Common variables
   
#  include "pb_def.h"
#  include "pb_md.h"
   
   ! Passed vriables

   integer ifcap, ipb, natom
   
   ! Local variables
   
   integer ip, iatm
   
   ! begin code
   
   write(6, '(a)'        )
   write(6, '(a)'        ) ' ======== FDPB Summary ========'
   write(6, '(a)'        )
   write(6, '(a,i12,x,a)') ' Do FDPB every', ndofd, 'steps'
   write(6, '(a)'        ) ' Nonbonded Update'
   write(6, '(a, f19.13)') '   residue cutoff is set to', sqrt(cutres)
   write(6, '(a, f19.13)') '   fdpb cutoff is set to', sqrt(cutfd)
   write(6, '(a, f19.13)') '   sas cutoff is set to', sqrt(cutsa)
   write(6, '(a, f19.13)') '   nonbonded cutoff is set to', sqrt(cutnb)
   write(6, '(a)'        ) ' Grid Constants'
   write(6, '(a, 3i8)') '   Grid dimension', xm, ym, zm
   write(6, '(a, f10.3)' ) '   Grid spacing set to', h
   write(6, '(a)'        ) '   Grid boundary'
   write(6, '(4x, 2f10.3)'   ) gox, gox + (xm+1) * h
   write(6, '(4x, 2f10.3)'   ) goy, goy + (ym+1) * h
   write(6, '(4x, 2f10.3)'   ) goz, goz + (zm+1) * h
   write(6, *)
   write(6, *) 'Dielectric Map'
   if ( ifcap /= 0 ) then
      write(6, *) '  A single spherical dielectric boundary is used'
   else
      if ( radiopt == 0 ) then
         write(6, *) '  Use cavity radii in the prmtop file'
      else if ( radiopt == 1 ) then
         write(6, *) '  Use Tan, Yang, and Luo optimized cavity radii definition'
      else if ( radiopt == 2 ) then
         write(6, *) '  Use Pymol radii definition for potential surface display'
      else
         write(6, *) 'PB Bomb in pb_init(): Unknown radiopt for cavity radii', radiopt
         call mexit(6, 1)
      end if
   end if
   if ( .not. srsas ) then
      write(6, *) '  The following dynamic PB radii are used at this step:'
      do iatm = 1, natom
         write(6, *) iatm, radi(iatm), radip3(iatm)
      end do
   end if
   write(6, *)
   if ( ipb == 1 ) then
      write(6, '(x,a)'       ) '  Use probe-based SES definition'
      write(6, '(x,a,i12,x,a)')'    Compute SAS/SAR every', ndosas, 'steps'
      write(6, '(x,a,f19.14)') '    Solvent probe radius', dprob
      write(6, '(x,a,i8)'    ) '    Surface dots per atom', maxsph
      write(6, '(x,a,f19.14)') '    Buried atom radii increment', radinc
      write(6, '(x,a,f19.14)') '    Threshhold for exposed atom', expthresh
      ! f19.14 for sas area is too small and you may see *********
      ! Haven't changed the format cause it may affect test cases
      ! by Qin
      write(6, '(x,a,f19.14)') '    Current SAS', prtsas
      !print *,'    Current SAS', prtsas
   else if ( ipb == 2 .or. ipb == 4 ) then
      write(6, *) '  Use level-set-based SES definition'
      write(6, *) '    Compute SAS/SAR every', ndosas, 'steps'
      write(6, *) '    Solvent probe radius', dprob
      write(6, *) '    Surface dots per atom', maxsph
   else
      write(6, *) '  Use undocumented definition'
   endif
   write(6, *)
   write(6, *) 'Boundary conditions'
   if ( bcopt == 1 ) then
      write(6, *) '  zero potential on boundary'
   else if ( bcopt == 2 ) then
      write(6, *) '  solute as a single DH sphere'
   else if ( bcopt == 3 ) then
      write(6, *) '  sum of residues as independent DH spheres'
   else if ( bcopt == 4 ) then
      write(6, *) '  sum of atoms as independent DH spheres'
   else if ( bcopt == 5 ) then
      write(6, *) '  sum of grid charges as independent DH spheres'
   else if ( bcopt == 0 ) then
      write(6, *) '  electrostatic focus boundary condition'
   end if
   write(6, '(a)')
   write(6, '(x,a)') 'Physical constants'
   write(6, '(x,a,f19.14)') '  Solute dielectric constant  :', epsin/eps0
   write(6, '(x,a,f19.13)') '  Solvent dielectric constant :', epsout/eps0
   write(6, '(x,a,f19.12)') '  Temperature (K)             :', pbtemp
   write(6, '(x,a,f19.14)') '  Ionic strength (mM)         :', fiono*istrng
   write(6, '(x,a,f19.14)') '  Debye-Huckel parameter (1/A):', pbkappa
   write(6, '(a)')
   write(6, '(a)') 'FD Solver Option'
   if (solvopt == 1) then
      write(6, *) '  Use Modified ICCG solver'
   else if (solvopt == 2) then
      write(6, *) '  Use MG solver'
   else if (solvopt == 3) then
      write(6, *) '  Use CG solver'
   else if (solvopt == 4) then
      write(6, *) '  Use SOR solver'
   else if (solvopt == 5) then
      write(6, *) '  Use Adaptive SOR solver'
   else if (solvopt == 6) then
      write(6, *) '  Use Damped SOR solver'
   end if
   write(6, '(a)')
   write(6, '(x,a)'        ) 'Iteration data'
   write(6, '(x,a,i12)'    ) '  Maximum iterations  :', maxitn
   write(6, '(x,a,es23.15)') '  Convergence criteria:', accept
   
   return
end subroutine pb_print

!  Wrapper for i/o buffer flushing routine
!  Author: George Seibel
!  Rewritten by: Meng-Juei Hsieh
!  Working for most Unix (BSD, Convex, Sun, Stellar, SGI Iris...)
subroutine amflsh(filenum)
   implicit none
   integer filenum ! unit file number to flush
   integer istat   ! return status from flush
#if defined(AIX) || defined( XLF90)
   call flush_(filenum) !page 222 in the V2.3 Language Reference of XLF
#else
#  ifdef SGI
      call flush(filenum,istat)
#  else
      call flush(filenum)
#  endif
#endif
   return
end subroutine amflsh
