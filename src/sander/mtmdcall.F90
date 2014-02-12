#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine nmrcal here]
subroutine mtmdcall(emtmd,mtmd_reals,mtmd_ints,x,f,name,irsnam,ipres,&
                    isymbl,mass,natom,nres,calltype)

  use multitmd
  use file_io_dat
  
   ! This is the top-level calling routine for the Multiple Target
   ! Targeted MD subroutines.  
   
   ! This was written by
   
   !            Matthew G. Seetin,
   !            Department of Biochemistry and Biophysics
   !            University of Rochester Medical Center, 14624
   !
   !     It is modeled after the nmr/nmrcal module by David Pearlman.
   
   
   !     In AMBER 10 and earlier, one could only do TMD relative to a single target.
   !     With this, the user will be able to do TMD with multiple reference structures,
   !     each with their own force constant and mask to which the RMSD is considered.
   !     That is, if one wants to force one group of atoms (say, one domain of a protein)
   !     to look like one reference structure, and a different group of atoms to look 
   !     like a different reference structure, that will now be possible.
   !
   !     The way this interfaces with the rest of the code is modeled after the NMR 
   !     restraint module: nmr.f along with nmrcal.f and nmr.h.
   
   !
   ! Subroutine Multi-TMD Call. This routine makes the calls to the following
   ! Multi-TMD related routines:
   
   !    MTMDALLOC : Does a cursory read of the mtmd input files to
   !                determine how much memory will be required. This
   !                will be used by locmem to dynamically allocate memory
   !     MTMDREAD : Reads information defining each target for MTMD,
   !                including the reference files, force constants, target
   !                RMSDs, and time-varying information
   !   MTMDENERGY : Calculates the energy/derivatives resulting from each
   !                target.  It evaluates the current values of the force
   !                constant and target RMSD and uses the RMSFIT and XTGTMD
   !                subroutines from the original TMD implimentation to 
   !                calculate the current RMSD, forces and energy.
   !    MTMDPRINT : Prints the current RMSD and the current target RMSD in 
   !                the out file.  Energy is reported in the RESTRAINT value
   !                of the out file.
   
   ! Calls (except to MTMDALLOC) are made thorugh a call to MTMDCALL
   ! with the appropriate value of CALLTYPE (see below). 
   
   ! Calls to MTMDALLOC are always made through the entry point MTMDLX.
   
   ! The entry point MTMDUNSTEP decrements the step-number counter and returns.
   
   ! INPUT:
   !    CALLTYPE: (Character*4 variable)
   !              CALLTYPE = "READ": Call to MTMDREAD to read mtmd file
   !              CALLTYPE = "CALC": Call to MTMDENERGY is made.
   !              CALLTYPE = "PRNT": Call to MTMDPRINT is made 
   !              CALLTYPE = "MPI ": Broadcasts MTMD common blocks and dynamic arrays
   !        X(I): Coordinate array. Coordinates (x,y,z) for atom j are packed in
   !              elements 3*(j-1)+1, 3*(j-1)+2 and 3*(j-1)+3.
   !        F(I): Force array. Packed the same as the coordinate array.
   !     NAME(I): Name array. NAME(I) is the name of atom I.
   !   IRSNAM(I): Residue name array. RSNAM(I) is name of residue I.
   !    IPRES(I): Residue pointer array. IPRES(I) points to the first atom
   !              of residue I.
   ! MTMD_REALS(I): Dynamically allocated array for reals.
   !                Requires  maxtgt * (natom + 4) storage.
   ! MTMD_INTS(I): Dynamically allocated array for integers.
   !               Requires maxtgt * (3*natom + 6) storage.
   !      MAXTGT: The maximum number of targets allowed for MTMD
   !        IOUT: The unit for informational prints
   
   ! OUTPUT:
   ! ------
   ! (The following are set when a call to NMRNRG [CALTYP="CALC"] is made).
   
   !     EMTMD:  The total energy contribution from each target
   
   ! The following are set when a call to RESTAL [entry point RESTLX] is made:
   
   !     TGTNUM :  The number of targets used.
   !    MAXTGT  :  The maximum number of targets allowed for MTMD.
   ! MTMDINTREQ :  The amount of integer work array (MTMD_INTS) storage required.
   ! MTMDIRLREQ :  The amount of real work arary (MTMD_REALS) storage required.
   
   ! Partitioning the work arrays:
   !   The work storage is partitioned here. The following variables
   !   are set to indicate the storage locations of the various arrays
   !   in the work vectors:
   !
   !  Reals:
   !         refcrds_ptr:  REFCRDS(*)
   !         mtmdrms_ptr:  MTMDRMSDARR(2,I)
   !         mtmdfrc_ptr:  MTMDFRC(2,I)
   !     tgtrmsdcurr_ptr:  TGTRMSDCURR(I)
   !        rmsdcurr_ptr:  RMSDCURR(I)

   ! Ints:
   !        mtmdmask_ptr:  MTMDMASK(*)
   !       mtmdsteps_ptr:  MTMDSTEPS(3,I)
   !     mtmdtgtatms_ptr:  MTMDTGTATMS(I)

   
   ! Other Local Common:
   !         NSTEP: Local accumulator which keeps track of the number of
   !         calls to MTMDENERGY. On each call, NSTEP is incremented by 1.
   
   implicit none
   
   
   integer B_MTMDI,mtmdlun,reflun,iout
   parameter (B_MTMDI=11)   ! size of mtmdloci common block
   !parameter (B_MTMDR=)    ! size of mtmdlocr common block (not yet used)
   parameter (mtmdlun = 33) ! same logical unit number as nmr restraint file.  
                            ! should be closed by the time this routine is called
   parameter (reflun = 10) ! same logical unit number as used for refc.  It 
                           ! should also have been closed by the time this is read
                           
   character(len=4) calltype, name, irsnam, isymbl
   dimension name(*), irsnam(*), isymbl(*)
   
   _REAL_ emtmd
   _REAL_ x, f, mtmd_reals, mass
   dimension x(*),f(*),mtmd_reals(*), mass(*)
   
   integer natom,nres
   parameter (iout = 6)
   
   integer mtmd_ints, ipres
   dimension mtmd_ints(*),ipres(*)
   
   integer tgtnum,maxtgt,refcrds_ptr,mtmdrms_ptr,mtmdfrc_ptr,&
           tgtrmsdcurr_ptr,rmsdcurr_ptr,mtmdmask_ptr,mtmdsteps_ptr,&
           mtmdtgtatms_ptr,nstep
   
   common/mtmdloci/tgtnum,maxtgt,refcrds_ptr,mtmdrms_ptr,mtmdfrc_ptr,&
                   tgtrmsdcurr_ptr,rmsdcurr_ptr,mtmdmask_ptr,mtmdsteps_ptr,&
                   mtmdtgtatms_ptr,nstep
   
   !common/mtmdlocr/
   
   logical runonce
   save runonce
   data runonce/.false./
   
#ifdef MPI
#  include "parallel.h"
   integer :: ierr
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif

#include "multitmd.h"
   
   if(.not.runonce) then
     nstep = 0
     
   ! mtmd_reals array
     refcrds_ptr = 1
     mtmdrms_ptr = refcrds_ptr + 3*natom*maxtgt
     mtmdfrc_ptr = mtmdrms_ptr + 2 * maxtgt
     tgtrmsdcurr_ptr = mtmdfrc_ptr + 2 * maxtgt
     rmsdcurr_ptr = tgtrmsdcurr_ptr + maxtgt
     ! rmsdcurr_ptr requires maxtgt storage
     
     ! mtmd_ints array
     mtmdmask_ptr = 1
     mtmdsteps_ptr = mtmdmask_ptr + natom*maxtgt
     mtmdtgtatms_ptr = mtmdsteps_ptr + 3 * maxtgt
     ! mtmdtgtatms_ptr requires maxtgt storage
     
     runonce = .true.
   end if
   
   select case(calltype)
     case('READ')
       call mtmdread(mtmd_reals(refcrds_ptr),mtmd_reals(mtmdrms_ptr),mtmd_reals(mtmdfrc_ptr),&
                     mtmd_ints(mtmdmask_ptr),mtmd_ints(mtmdsteps_ptr),mtmd_ints(mtmdtgtatms_ptr),&
                     x,name,irsnam,ipres,isymbl,&
                     natom,nres,tgtnum,maxtgt,mtmdlun,mtmd,reflun)
     case('CALC')
       call mtmdenergy(emtmd,mtmd_reals(refcrds_ptr),mtmd_reals(mtmdrms_ptr),&
                     mtmd_reals(mtmdfrc_ptr),mtmd_ints(mtmdmask_ptr),mtmd_ints(mtmdsteps_ptr),&
                     mtmd_ints(mtmdtgtatms_ptr),x,f,&
                     mtmd_reals(tgtrmsdcurr_ptr),mtmd_reals(rmsdcurr_ptr),name,irsnam,&
                     mass,natom,nres,tgtnum,maxtgt,nstep)
       nstep = nstep + 1
     case('PRNT')
       call mtmdprint(tgtnum, mtmd_reals(rmsdcurr_ptr), mtmd_reals(tgtrmsdcurr_ptr), iout)
#ifdef MPI
     case ('MPI ')
      
      !     ...when using current AMBER/MPI implementation, only the master
      !     processors reads and processes the input file.  This implies
      !     that only the master processor sets anything in the MTMDCALL "READ"
      !     step or in the MTMDALLOC entry point.  Therefore, we need to give
      !     the other processors the data.  This is done by adding the "MPI "
      !     CALLTYPE for MTMDCALL() and calling this *after* startup (in order that
      !     the mtmd and locmem pointers, etc be already properly broadcast).
      
      call mpi_bcast(tgtnum,B_MTMDI,mpi_integer,0,commsander,ierr)
      !call mpi_bcast(??,B_MTMDR,MPI_DOUBLE_PRECISION,0,commsander,ierr)
      call mpi_bcast(mtmd_reals(1),rmsdcurr_ptr + maxtgt,MPI_DOUBLE_PRECISION, &
            0,commsander,ierr)
      call mpi_bcast(mtmd_ints(1),mtmdtgtatms_ptr + maxtgt,mpi_integer,0,commsander,ierr)
#endif
     case default
      write(iout,9000) calltype
      call mexit(iout, 1)
   end select
   
   return
   
   entry mtmdlx(natom)
   
   call mtmdalloc(mtmd,reflun,tgtnum)

   maxtgt = tgtnum
   
   mtmdintreq = maxtgt * (natom + 4) ! Mask + start and end steps of variation
   mtmdirlreq = maxtgt * (3*natom + 6) ! Coords + force const + target RMSD
   
   return 
   
   entry mtmdunstep
   nstep = nstep - 1
   return
   
   9000 format(' Error: Unrecognized option passed to MTMDCALL: CALLTYPE = ', a4)
end subroutine mtmdcall

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine restal here]
subroutine mtmdalloc(tgtin,mtmdlun,tgtnum)

  ! Subroutine Multiple Targeted MD memory ALLOCation
  
  ! This subroutine does a quick read of the multiple-target input file
  ! and returns the number of targets TGTNUM
   
   ! These quantities are required in setting up the dynamic memory
   ! allocation in routine LOCMEM
   
   ! It is expected that this routine will be called before any other of the
   ! NMR suite.
   
   ! Author: Matthew G. Seetin
   ! Date: 6/2008
   
   ! INPUT:
   !       TGTIN : The file name of the mtmd file
   !     MTMDLUN : Logical unit number for mtmd file
   !        IOUT : Unit for informational prints.
   
   ! OUTPUT:
   !      TGTNUM : The total number of targets defined.

  use file_io_dat, only : MAX_FN_LEN
  
  implicit none
  
  integer mtmdlun,tgtnum
  character(MAX_FN_LEN) tgtin
  
  character(MAX_FN_LEN) refin
  character(256) mtmdmask
  integer mtmdstep1
  integer mtmdstep2
  _REAL_ mtmdforce
  _REAL_ mtmdforce2
  _REAL_ mtmdrmsd
  _REAL_ mtmdrmsd2
  integer mtmdform
  integer mtmdninc
  integer mtmdmult
  integer mtmdvari
  integer ifind
  integer i
  
  namelist /tgt/ refin,mtmdmask,mtmdstep1,mtmdstep2,mtmdforce,mtmdforce2,&
                mtmdrmsd,mtmdrmsd2,mtmdform,mtmdninc,mtmdmult,mtmdvari
                
  refin = ' '
  mtmdmask = ' '
  mtmdstep1 = 0
  mtmdstep2 = 0
  mtmdforce = 0.0d0
  mtmdforce2 = 0.0d0
  mtmdrmsd = 0.0d0
  mtmdrmsd2 = 0.0d0
  mtmdform = 1
  mtmdninc = 0.0
  mtmdvari = 0
  mtmdmult = 0
  ifind = 0
  tgtnum = 0
  
  call amopen(mtmdlun,tgtin,'O','F','R')
  
  do i=1,999999
    refin = ' '
    
    call nmlsrc('tgt',mtmdlun,ifind)
    if (ifind == 0) exit
    
    ! Read the tgt namelist
    read (mtmdlun,nml=tgt,end=100)
    
    ! Only count the restraint if there was an input reference file defined.
    ! If there was no reference file given in the namelist, assume that there
    ! will be no further input.
    if(refin == ' ') then
      exit
    else
      tgtnum = tgtnum + 1
    end if
  end do
  
  100 continue
  
  close(mtmdlun)

end subroutine mtmdalloc
