#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Count number of *active* atoms in belly
subroutine setatm(nat,natb,igrp)
   implicit none
   ! Passed variables
   integer, intent(in)  :: nat, igrp(*)
   integer, intent(out) :: natb
   ! Internal variables
   integer :: i,idum

   idum = 0
   do i = 1,nat
      if(igrp(i) > 0) idum = idum+1
   end do
   natb = idum
   
   return
end subroutine setatm 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Get rid of bond force constant for QM/MM or non-belly atom pairs.
subroutine setbon(nb,ib,jb,icb,igrp)
   use qmmm_module, only : qmmm_nml,qmmm_struct
   use parms, only: numbnd,max_bond_type,rk,req
   use constants, only : zero
   implicit none
   ! Passed arguments
   integer, intent(inout) :: nb
   integer, intent(inout) :: ib(nb), jb(nb), icb(nb)
   integer, intent(in)    :: igrp(*)

   ! Internal variables
   integer nba, i, iat, jat
   logical iq, jq
   integer j,iatnum,jatnum
  
!numbnd = number of bond parameters
 
   nba = 0
   do i = 1,nb
      iat = ib(i)/3+1
      jat = jb(i)/3+1
      iq = .FALSE.
      jq = .FALSE.
      if (qmmm_nml%ifqnt) then
        ! remove quantum - quantum bond pairs from the bond list. This is done by rebuilding
        ! the list without these pairs.
        ! Note if lnk_method == 2 then we will treat the MMLink pair atom as being a QM atom so
        ! we also delete any MML-QM bonds.
        do j=1, qmmm_struct%nquant
           if (iat==qmmm_struct%iqmatoms(j)) then
              iq = .true.
              iatnum = qmmm_struct%iqm_atomic_numbers(j)
           elseif (jat==qmmm_struct%iqmatoms(j)) then
              jq = .true.
              jatnum = qmmm_struct%iqm_atomic_numbers(j)
           end if
        end do
        if (qmmm_nml%lnk_method == 2) then
          ! MMlink atoms treated as QM atoms.
          ! Need to check if either iat or jat is an MM link pair atom.
          if (qmmm_struct%mm_link_mask(iat)) then
            iq = .true.
            iatnum = -1  !set iatnum to -1 for MMlink atoms since we don't want these bonds getting shaken.
          elseif (qmmm_struct%mm_link_mask(jat)) then
            jq = .true.
            jatnum = -1  !set jatnum to -1 for MMlink atoms since we don't want these bonds getting shaken.
          end if 
        end if
        iq = iq .and. jq
        !iq will be true if BOTH iat and jat are in the QM region
      end if
      if((igrp(iat) > 0 .or. igrp(jat) > 0)) then
        if (iq) then !Both are QM atoms (or also QM-MMLink if lnk_method==2)
          !In order to shake QM atoms we need to add the bonds to the bond list
          !but we need to make the force constant zero while preserving the eqm bond length.
          !We will do this by creating a new type for each of the QM-QM bonds.
          !For the moment only do QM-H bonds for shake NTC=2
          !Only do this if qmshake is set.
          if ((iatnum == 1 .OR. jatnum == 1) .AND. qmmm_nml%qmshake > 0) then
            !at least one of them is a hydrogen
            !We need to make a new bond type here.
            numbnd = numbnd+1
            if (numbnd > MAX_BOND_TYPE) then
               !exceeded the maximum number of bond types
               call sander_bomb('setbon (set.f)', &
               'Exceeded MAX_BOND_TYPE while setting up bond types for QM atom shake', &
               'Increase MAX_BOND_TYPE in parms.f and recompile.')
            end if
            rk(numbnd) = zero  !Set force constant to zero
            req(numbnd) = req(icb(i)) !Preserve the eqm distance
            nba = nba+1
            ib(nba) = ib(i)
            jb(nba) = jb(i)
            icb(nba) = numbnd
          !else we do nothing and the bond doesn't get added to the bond list.
          end if 
        else
          !We add the current pair to the list if both atoms are not bellied AND/OR
          !both in the QM region
          nba = nba+1
          ib(nba) = ib(i)
          jb(nba) = jb(i)
          icb(nba) = icb(i)
        end if
      end if
   end do
   nb = nba
   return
end subroutine setbon 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Remove QM/MM or non-belly angle parameters
subroutine setang(nt,it,jt,kt,ict,igrp)
   use qmmm_module, only : qmmm_nml,qmmm_struct
   implicit NONE
   integer nt
   integer it(nt),jt(nt),kt(nt),ict(nt),igrp(*)
   integer nta, i, iat, jat, kat
   logical iq, jq, kq
   integer j
   
   nta = 0
   do i = 1,nt
      iat = it(i)/3+1
      jat = jt(i)/3+1
      kat = kt(i)/3+1
      iq = .FALSE.
      jq = .FALSE.
      kq = .FALSE.
      if (qmmm_nml%ifqnt) then
        ! remove quantum - quantum - quantum angle triplets from the angle list. This is done by rebuilding
        ! the list without these triplets.
        ! if lnk_method=2 we treat the MMlink pair atoms as being QM so we also remove ML-QM-QM angle triplets.
        do j=1, qmmm_struct%nquant
          iq = iq.or.iat==qmmm_struct%iqmatoms(j)
          jq = jq.or.jat==qmmm_struct%iqmatoms(j)
          kq = kq.or.kat==qmmm_struct%iqmatoms(j)
        end do
        if (qmmm_nml%lnk_method == 2) then
          ! MMlink atoms treated as QM atoms.
          ! Need to check if iat, jat or kat are an MM link pair atom.
          iq = iq .or. qmmm_struct%mm_link_mask(iat)
          jq = jq .or. qmmm_struct%mm_link_mask(jat)
          kq = kq .or. qmmm_struct%mm_link_mask(kat)
        end if
        iq = iq .and. jq .and. kq
        !iq will be true if iat,jat and kat are all in the QM region
      end if
      if((igrp(iat) > 0 .or. igrp(jat) > 0 .or. igrp(kat) > 0) .and. (.not. iq) ) then
        !We add the current triplet to the list if all atoms are not bellied AND/OR
        !in the QM region
        nta = nta+1
        it(nta) = it(i)
        jt(nta) = jt(i)
        kt(nta) = kt(i)
        ict(nta) = ict(i)
      end if
   end do
   nt = nta
   return
end subroutine setang 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setdih here]
subroutine setdih(np,ip,jp,kp,lp,icp,igrp)
   use qmmm_module, only : qmmm_nml,qmmm_struct
   implicit NONE
   integer np
   integer ip(np),jp(np),kp(np),lp(np),icp(np),igrp(*)
   integer npa, i, iat, jat, kat, lat
   logical iq, jq, kq, lq
   integer j
   
   npa = 0
   do i = 1,np
      iat = ip(i)/3+1
      jat = jp(i)/3+1
      kat = iabs(kp(i))/3+1
      lat = iabs(lp(i))/3+1
      iq = .FALSE.
      jq = .FALSE.
      kq = .FALSE.
      lq = .FALSE.
      if (qmmm_nml%ifqnt) then
        ! remove quantum - quantum - quantum - quantum dihedrals from the dihedra; list. This is done by rebuilding
        ! the list without these dihedrals.
        do j=1, qmmm_struct%nquant
          iq = iq.or.iat==qmmm_struct%iqmatoms(j)
          jq = jq.or.jat==qmmm_struct%iqmatoms(j)
          kq = kq.or.kat==qmmm_struct%iqmatoms(j)
          lq = lq.or.lat==qmmm_struct%iqmatoms(j)
        end do
        if (qmmm_nml%lnk_method == 2) then
          ! MMlink atoms treated as QM atoms.
          ! Need to check if iat, jat, kat or lat are an MM link pair atom.
          iq = iq .or. qmmm_struct%mm_link_mask(iat)
          jq = jq .or. qmmm_struct%mm_link_mask(jat)
          kq = kq .or. qmmm_struct%mm_link_mask(kat)
          lq = lq .or. qmmm_struct%mm_link_mask(lat)
        end if
        iq = iq .and. jq .and. kq .and. lq
        !iq will be true if iat,jat, kat and lat are all in the QM region
      end if
      if((igrp(iat) > 0 .or. igrp(jat) > 0 .or. igrp(kat) > 0 .or. igrp(lat) > 0) .and. (.not. iq)) then
        !add current dihedral set to the list since all atoms are not bellied or in the QM region
        npa = npa+1
        ip(npa) = ip(i)
        jp(npa) = jp(i)
        kp(npa) = kp(i)
        lp(npa) = lp(i)
        icp(npa) = icp(i)
      end if
   end do
   np = npa
   return
end subroutine setdih 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read coordinates from a restart file
!-------------------------------------------------------------------
!     --- RDREST ---
!-------------------------------------------------------------------
!     Read coordinates from a restart file, the format of which is 
!     determined by ntrx:
!       0 - Unformatted (binary) restart
!       1 - Formatted (ascii) restart
!       2 - Netcdf restart
!     nr is the # atoms, filename is the path to the file, xc is
!     the array of coordinates to be read to.
subroutine rdrest(nr,ntrx,filename,xc)
   
   use file_io_dat
   use binrestart, only: read_nc_restart,check_nc_restart
   implicit none
   
   !     ----- ROUTINE TO READ THE REFERENCE POSITIONS FOR RESTRAINING ----
   integer, intent(in) :: nr,ntrx
   character(len=*), intent(in) :: filename
   _REAL_, intent(out) :: xc(*)
   !--- local -------
   character(len=80) :: line
   integer :: nr3,natom,i3,i,reflun
   ! Dummy variables for reading temp/time from netcdf restart
   _REAL_ :: ltemp0,restarttime

   ! reflun is the logical unit that will be associated with filename.
   ! It is assigned the value of 10 since this routine was ALWAYS called
   ! with reflun = 10, but it really should have its definition in files.h
   ! for consitency
   reflun = 10
   
   nr3 = 3*nr
   write(6,191)
   
   if (check_nc_restart(filename)) then
     ! Netcdf file - only reading coordinates, like ntx = 1, so the second
     ! instance of xc (normally velocities) will not be written to and its
     ! safe to pass it in as a dummy arg.
     call read_nc_restart(filename,title1,1,nr,xc,xc,ltemp0,restarttime)
   else if(ntrx == 0) then
     ! Binary file
     call amopen(reflun,filename,'O','U','R') 
     read(reflun) title1
     read(reflun) natom
     if(natom /= nr) then
       write(6,1000)
       call mexit(6, 1)
     endif
     read(reflun) (xc(i3),i3 = 1,nr3)
     close(reflun)
   else
     ! Formatted file
     call amopen(reflun,filename,'O','F','R')
     read(reflun,40) title1
     read(reflun,'(a)') line

     if( line(6:6) == ' ' ) then   ! this is an old, i5 file
        read(line,'(i5)') natom
     elseif( line(7:7) == ' ' ) then ! sander 7/8/9/10 large system format...
        read(line,'(i6)') natom
     elseif( line(8:8) == ' ' ) then ! Sander 11 - 1 mil+ format
        read(line,'(i7)') natom
     else                          ! assume amber 11 VERY large system format. 10 mil+
        read(line,'(i8)') natom
     end if
     if(natom /= nr) then
       write(6,1000)
       call mexit(6, 1)
     endif
     read(reflun,92) (xc(i),i = 1,nr3)
     close(reflun)
   endif
   write(6,31) title1
   return
   191 format(/,'   5.  REFERENCE ATOM COORDINATES',/)
   31 format(2x,20a4)
   40 format(20a4)
   92 format(6f12.7)
   1000 format(/2x,'FATAL: NATOM mismatch in constraint coord', &
         ' and topology files')
end subroutine rdrest 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dihdup here]
subroutine dihdup(nphia,ip,jp,kp,lp,icp,pn)
   
   !     duplicates pointers to multi-term dihedrals for vector ephi.
   !     H-atom diheds are duplicated at the end of the h-atom lists,
   !     but heavy atom diheds must be inserted between the end of the
   !     heavy-atom diheds, but before the constraint diheds if they are
   !     present.  In order to use this code, extra space MUST be allocated
   !     in LOCMEM
   
   !     (Note: this is only for ancient prmtop files: anything created by
   !     LEaP should report no duplications)
   
   !     Author:  George Seibel
   
   implicit none
   
   !     COMMON:
   
#  include "md.h"
   
   !     INPUT:
   
   integer :: nphia
   !        ... number of dihedrals
   integer ip(*), jp(*), kp(*), lp(*)
   !        ... pointers to atoms of dihedrals
   integer icp(*)
   !        ... pointers to dihedral parameters
   _REAL_ :: pn(*)
   !        ... periodicity; negative if multi-term, read until + encountered
   
   !     INTERNAL:
   
   integer ndup
   !        ... number of diheds duplicated
   integer ic
   !        ... working parameter pointer
   integer i,ierr
   !        ... do loop index
   
   ndup = 0
   ierr = 0
   
   do i = 1, nphia
      ic = icp(i)
      100 continue
      if (pn(ic) < 0) then
         ndup = ndup + 1
         if (ndup > maxdup) then
            ierr = 1
         else
            
            !            --- duplicate pointers of multi-term dihedrals ---
            
            ip(nphia+ndup)  = ip(i)
            jp(nphia+ndup)  = jp(i)
            kp(nphia+ndup)  = kp(i)
            lp(nphia+ndup)  = lp(i)
            
            !            --- but use the NEXT parameter pointer ---
            
            icp(nphia+ndup) = ic + 1
         end if
         
         !            --- check for a third or higher term ---
         
         ic = ic + 1
         goto 100
      end if
   end do
   
   if (ierr == 1) then
      write(6,'(/,5x,a,i5,a)') 'MAXDUP =',maxdup,' exceeded'
      write(6,'(/,5x,a,i5,a)') 'set MAXDUP =',ndup,' in locmem.f'
      call mexit(6, 1)
   end if
   
   nphia = nphia + ndup
   write(6,'(a,i5,a)') '| Duplicated',ndup,' dihedrals'
   return
end subroutine dihdup 
!-----------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine dihpar here]
subroutine dihpar(numphi,pk,pn,phase,gamc,gams,ipn)
   use constants, only : zero, one, four, TEN_TO_MINUS3, TEN_TO_MINUS6, PI
   implicit none
   
   !     ----- ROUTINE TO GET ADDITIONAL PARAMETERS FOR THE
   !           VECTORISED EPHI -----
   
   integer :: numphi, ipn(*), kp(numphi), lp(numphi)
   _REAL_ :: pk(*),pn(*),phase(*),gamc(*),gams(*)

!----- local -----------------   
   _REAL_ :: pim,dum,dumc,dums
   integer i
   
   pim = four*atan(one)
   do i = 1,numphi
      dum = phase(i)
      if(abs(dum-pi) <= TEN_TO_MINUS3) dum = sign(pim,dum)
      dumc = cos(dum)
      dums = sin(dum)
      if(abs(dumc) <= TEN_TO_MINUS6) dumc = zero
      if(abs(dums) <= TEN_TO_MINUS6) dums = zero
      gamc(i) = dumc*pk(i)
      gams(i) = dums*pk(i)

      !Build the 1-4 mask for all the dihedral types
      !pn=-ve = duplicate dihedral
      !if 4th atom of dihedral is -ve then this is an improper
      !if 3rd atom of dihedral is -ve then this implies that it is part
      !of a ring or similar system where the 1-4's should be ignored.
      !This array essentially contains a true if it is not a duplicate.

      !RCW: Note pn is always positive according to the prmtop file
      !     definition since duplicate 1-4's are dealt with by
      !     setting entry 3 or 4 of each dihedral type -ve. PN is
      !     only set to be negative in the frcmod / parm file to tell
      !     leap to look for duplicate dihedrals. Beyond leap it gets
      !     ignored. Earlier versions of AMBER used to set FMULN to zero
      !     if pn was negative but I think this is redundant - perhaps
      !     some VERY OLD prmtop's still have negative pn values but
      !     none of the test cases show problems if pn being negative
      !     is ignored. It is still set as the absolute value here
      !     just to be sure.
      !    if (pn(i) .le. zero) fmn(i) = zero
      pn(i) = abs(pn(i))
      ipn(i) = int(pn(i)+TEN_TO_MINUS3)
   end do

   

   return
end subroutine dihpar 
!-----------------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setgms here]
subroutine setgms(numphi)
   
   ! Subroutine SET GaMS.
   
   ! This routine sets the values of the GAMC() and GAMS() arrays, which
   ! are used in vectorized torsional energy routines. This routine is only
   ! called when the torsional force constants are changed in routine MODWT.
   ! Otherwise, these arrays are set only once, by a call to DIHPAR from RDPARM,
   ! at the start of the program. This is an abbreviated version of DIHPAR which
   ! passes most arguments by common, not call-list.
   
   ! Author: David A. Pearlman
   ! Date: 5/92
   use parms, only: gams, gamc, phase, pk
   use constants, only : zero, one, four, TEN_TO_MINUS3, TEN_TO_MINUS6, PI
   implicit none

   integer :: numphi
!---- local -----------------   
   integer i
   _REAL_ pim
   _REAL_ :: dumc,dums,dum
   
   pim = four*atan(one)
   do i = 1,numphi
      dum = phase(i)
      if(abs(dum-pi) <= TEN_TO_MINUS3) dum = sign(pim,dum)
      dumc = cos(dum)
      dums = sin(dum)
      if(abs(dumc) <= TEN_TO_MINUS6) dumc = zero
      if(abs(dums) <= TEN_TO_MINUS6) dums = zero
      gamc(i) = dumc*pk(i)
      gams(i) = dums*pk(i)
   end do
   return
end subroutine setgms 
#ifdef MPI

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Distribute atoms to processors using molecule boundaries
subroutine setpar(nspm, nsp, ntp, ipres, amass)
   
   implicit none
   integer nspm, nsp(*), ntp, ipres(*)
   _REAL_ amass(*)

#  include "memory.h"
#  include "parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#  include "extra.h"
#  include "nmr.h"
   integer target,i,iat,imol,ipmol(nspm),node,j,ires,fraction


   if (mpi_orig) then
      
      !   when only master will be running the runmd,runmin routines:
      
      iparpt(0) = 0
      iparpt(1) = natom

   else
      
      !   all nodes will divide the work in runmd,runmin:
      
      if( ntp > 0 ) then   !  constant pressure run, divide on molecules:

         if (nspm < numtasks) then
            write(6,*) 'Must have more molecules than processors!'
            call mexit(6,1)
         end if

         !   set up an ipmol array, giving the first atom in each molecule;
         !     (by analogy to the ipres array)

         ipmol(1) = 1
         do imol=2,nspm
            ipmol(imol) = ipmol(imol-1) + nsp(imol-1)
         end do

         iparpt(0) = 0
         iparpt(numtasks) = natom
         do node=1,numtasks-1
            target = iparpt(node-1) + (natom-iparpt(node-1))/(numtasks-(node-1))
            do imol=2,nspm
               if (ipmol(imol) > target) exit
            end do

            !   --- this molecule starts past the target atom; hence, end the 
            !       list at the end of the previous molecule:

            iparpt(node) = ipmol(imol) - 1

         end do
      
      else   !  this is not a constant pressure run; divide on residues:

         if (nres < numtasks) then
            write(6,*) 'Must have more residues than processors!'
            call mexit(6,1)
         end if
         fraction = natom/numtasks
         target = 0
         iparpt(0) = 0
         iparpt(numtasks) = natom
         nodes: do node=1,numtasks-1
            target = target + fraction
            residues: do ires=1,nres

               !  don't stop before a residue whose single atom is a hydrogen:
            
               if ((ipres(ires+1) - ipres(ires)) < 2) then
                  iat = 0
                  do i=1,ires
                     iat = iat + ipres(i+1) - ipres(i)
                  enddo
                  if (amass(iat) < 3.0) cycle residues
               endif

               !  don't stop after a residue whose single atom is a hydrogen:

               if ( ires > 1 )then
                  if( (ipres(ires) - ipres(ires-1)) < 2 ) then
                     iat = 0
                     do i=1,ires-1
                        iat = iat + ipres(i+1) - ipres(i)
                     enddo
                     if (amass(iat) < 3.0) cycle residues
                  endif
               endif
 
               !   --- if this residue starts past the target atom, end 
               !       the atom list at the end of the previous residue:
            
               if (ipres(ires) > target) then
                  iparpt(node) = ipres(ires) - 1
                  cycle nodes
               end if
            end do residues

            !  should never get here, unless single atom residues are really
            !  messed up somehow....
            write(6,*) 'ERROR IN SETPAR() upon atom distribution'
            call mexit(6,1)

         end do nodes
      
      end if  !  (ntp > 0 )

      !       --- following used when forces, coords. communicated:
      
      do node=0,numtasks
         iparpt3(node) = 3*iparpt(node)
      end do
      
      !  --- include "extra" dynamical variables for final processor:
      
      iparpt3(numtasks) = iparpt3(numtasks) + iscale
      
      do i=0,numtasks-1
         rcvcnt(i) = iparpt(i+1) - iparpt(i)
         rcvcnt3(i) = iparpt3(i+1) - iparpt3(i)
      end do
      if (master) then
         write(6,'(a)') '|  Atom division among processors:'
         write(6,'("|  ", 8i8)') (iparpt(j),j=0,numtasks)
      end if
   end if  
   return
end subroutine setpar 
#endif

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setvar here]
subroutine setvar(ix,belly)
   
   implicit none
   logical belly
   integer i
   
   !     ----- ROUTINE TO DO THE NECESSARY ACCOMODATIONS FOR PROTEIN
   !           BELLY MINIMISATIONS -----
   
#  include "memory.h"
#  include "box.h"
   
   integer ix(*)
   
   !     --- SETUP THE BELLY GROUP ARRAY AND RESIDUE ARRAY for resnba ---
   
   if(belly) goto 190
   do i = 1,natom
      ix(ibellygp+i-1) = 1
   end do
   return
   
   !     ----- DELETE BONDS WHICH ARE IN THE BELLY ALONE -----
   
   190 continue
   call setbon(nbonh,ix(iibh),ix(ijbh),ix(iicbh),ix(ibellygp))
   call setbon(nbona,ix(iiba),ix(ijba),ix(iicba),ix(ibellygp))
   
   !     ----- MAKE THE BOND ARRAYS SEQUENTIAL FOR SHAKE AND FORCE
   !           ROUTINES -----
   
   do i = 1,nbona
      ix(iibh+nbonh+i-1)  = ix(iiba+i-1)
      ix(ijbh+nbonh+i-1)  = ix(ijba+i-1)
      ix(iicbh+nbonh+i-1) = ix(iicba+i-1)
   end do
   iiba = iibh+nbonh
   ijba = ijbh+nbonh
   iicba = iicbh+nbonh
   
   !     ----- DELETE THE BONDS WHICH ARE IN THE BELLY ALONE -----
   
   call setang(ntheth,ix(i24),ix(i26),ix(i28),ix(i30), &
         ix(ibellygp))
   call setang(ntheta,ix(i32),ix(i34),ix(i36),ix(i38), &
         ix(ibellygp))
   
   !     ----- MAKE THE ANGLE ARRAYS SEQUENTIAL -----
   
   do i = 1,ntheta
      ix(i24+ntheth+i-1) = ix(i32+i-1)
      ix(i26+ntheth+i-1) = ix(i34+i-1)
      ix(i28+ntheth+i-1) = ix(i36+i-1)
      ix(i30+ntheth+i-1) = ix(i38+i-1)
   end do
   i32 = i24+ntheth
   i34 = i26+ntheth
   i36 = i28+ntheth
   i38 = i30+ntheth
   
   !     ----- DELETE THE DIHEDRALS -----
   
   call setdih(nphih,ix(i40),ix(i42),ix(i44),ix(i46), &
         ix(i48),ix(ibellygp))
   call setdih(nphia,ix(i50),ix(i52),ix(i54),ix(i56), &
         ix(i58),ix(ibellygp))
   
   !     ----- MAKE THE DIHEDRALS SEQUENTIAL ------
   
   do i = 1,nphia
      ix(i40+nphih+i-1) = ix(i50+i-1)
      ix(i42+nphih+i-1) = ix(i52+i-1)
      ix(i44+nphih+i-1) = ix(i54+i-1)
      ix(i46+nphih+i-1) = ix(i56+i-1)
      ix(i48+nphih+i-1) = ix(i58+i-1)
   end do
   i50 = i40+nphih
   i52 = i42+nphih
   i54 = i44+nphih
   i56 = i46+nphih
   i58 = i48+nphih
   
   !     ----- FIND THE TOTAL NUMBER OF ACTIVE ATOMS AND RESIDUES -----
   
   call setatm(natom,natbel,ix(ibellygp))
   
   !     ----- ALL DONE RETURN -----
   
   return
end subroutine setvar 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine setnoshake(ix,noshakegp,ntc,num_noshake)
   
   use qmmm_module, only : qmmm_nml, qmmm_struct
   implicit none
   integer ix(*),noshakegp(*),ntc,num_noshake
   
   !     ----- don't shake any bonds with atoms in the noshakegp list
   
#include "memory.h"
#include "box.h"
#include "nmr.h"
   
   integer iano,jano,i,j,maxbond
   logical ilap, jlap

   ! write(6,'(20i3)') (noshakegp(i), i=1,natom)

   num_noshake = 0
   ix(noshake:noshake+nbonh-nbona-1) = 0

   if( ntc == 1 ) then
      return
   else if( ntc == 2 ) then
      maxbond = nbonh
   else if( ntc == 3 ) then
      maxbond = nbonh + nbona
   end if

   if (qmmm_nml%ifqnt) then
      qmmm_struct%noshake_overlap = 0
   end if

   !  loop over all of the bonds:
   do i = 1,maxbond
      iano = ix(iibh+i-1)/3 + 1
      jano = ix(ijbh+i-1)/3 + 1
      ilap = .FALSE.
      jlap = .FALSE.
      if( noshakegp( iano ) == 1  .or. noshakegp( jano ) == 1 ) then
         ix(noshake+i-1) = 1
         write(6,'(a,a,a,a)') '   Removing shake constraints from ', &
              resat(iano)(1:13),' -- ',resat(jano)(1:13)
         num_noshake = num_noshake + 1
         if (qmmm_nml%ifqnt) then
            do j = 1, qmmm_struct%nquant
               if (iano == qmmm_nml%iqmatoms(j)) then
                  if (qmmm_nml%qmshake == 0) then
                     ilap = .TRUE.
                  else ! qmshake > 0
                     call sander_bomb('setnoshake (in set.f)', &
                          'there are QM atoms in noshakemask', &
                          'in conflict with qmshake /= 0')
                  end if
               elseif (jano == qmmm_nml%iqmatoms(j)) then
                  if (qmmm_nml%qmshake == 0) then
                     jlap = .TRUE.
                  else
                     call sander_bomb('setnoshake (in set.f)', &
                          'there are QM atoms in noshakemask', &
                          'in conflict with qmshake /= 0')
                  end if
               end if
            end do
            if (ilap .and. jlap) then
               qmmm_struct%noshake_overlap = qmmm_struct%noshake_overlap + 1
            end if
         end if
      end if
   end do
   return

end subroutine setnoshake

#ifdef MPI /* SOFT CORE */
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine setnoshake_sc(ix,ntc,num_noshake,master)

   use softcore, only : nsc
   
   implicit none
   integer ix(*),num_noshake,ntc
   logical master
   
   !     ----- don't shake any bonds that cross between SC and common atoms
   
#include "memory.h"
#include "box.h"
#include "nmr.h"
   
   integer iano,jano,i,maxbond

   if( ntc == 1 ) then
      return  
   else if( ntc == 2 ) then
      maxbond = nbonh
   else if( ntc == 3 ) then
      maxbond = nbonh + nbona
   end if

   if ( master ) write (6,'(a)') '     Checking for SHAKE constraints on bonds crossing into the SC region'

   !  loop over all of the bonds:
   do i = 1,maxbond
      iano = ix(iibh+i-1)/3 + 1
      jano = ix(ijbh+i-1)/3 + 1
      if( ( nsc( iano ) == 1 .and. nsc( jano ) == 0 ) .or. &
          ( nsc( iano ) == 0 .and. nsc( jano ) == 1 )) then
         ix(noshake+i-1) = 1
         if ( master ) write(6,'(a,a,a,a)') '   Removing shake constraints from ', &
              resat(iano)(1:13),' -- ',resat(jano)(1:13)
         num_noshake = num_noshake + 1
      end if
   end do
   return

end subroutine setnoshake_sc
#endif
