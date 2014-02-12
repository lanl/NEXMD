!+ Specification and control of Amber's Input/Output

! File names
integer, parameter :: MAX_FN_LEN = 256
integer, parameter :: MAX_LINE_BUF_LEN = 256

character(len=4096) groupbuffer
character(len=MAX_FN_LEN) mdin, mdout, inpcrd, parm, restrt, &
      refc, mdvel, mden, mdcrd, mdinfo, mtmd, nmrf, mincor, &
      vecs, radii, freqe,redir(9),rstdip,mddip,inpdip,groups,gpes, &
      cpin, cpout, cprestrt, evbin, evbout, mmtsb_setup_file,pimdout, &
      inptraj

character owrite, facc
common /files/ groupbuffer, mdin, mdout, inpcrd, parm, restrt, &
      refc, mdvel, mden, mdcrd, mdinfo, mtmd, nmrf, mincor, &
      vecs, radii, freqe, owrite, facc,rstdip,mddip,inpdip,groups,gpes, &
      cpin, cpout, cprestrt, evbin, evbout, mmtsb_setup_file,pimdout, &
      inptraj

#ifdef RISM
#  include "../rism/files.h"
#endif

! put this in a separate common block to stop the compiler from
! complaining about misalignment
integer numgroup, nslice
common/nmgrp/ numgroup, nslice

! File units
! An I/O Unit resource manager does not exist.
integer     MDCRD_UNIT
integer     INPTRAJ_UNIT
integer     MDEN_UNIT
integer     MDINFO_UNIT
integer     MDVEL_UNIT
parameter ( MDINFO_UNIT =  7 )
parameter ( MDCRD_UNIT  = 12 )
parameter ( INPTRAJ_UNIT = 24 )
parameter ( MDEN_UNIT   = 15 )
parameter ( MDVEL_UNIT  = 13 )
integer, parameter :: CNSTPH_UNIT = 18, CPOUT_UNIT = 19
integer, parameter :: INPCRD_UNIT = 9;

! 18 was picked because CNSTPH uses it; conflicts are not expected.
integer     MMTSB_UNIT
parameter ( MMTSB_UNIT = 18 )

!!
!! EVB I/O unit
!!
   integer, parameter :: evb_unit = 75
   integer, parameter :: schlegel_unit = 80

!! FULL PIMD I/O unit
   integer, parameter :: pimd_unit = 277
! File related controls and options
character(len=80) title,title1
common/runhed/ title, title1

logical mdin_ewald,mdin_pb,mdin_amoeba

#ifdef APBS
logical mdin_apbs, sp_apbs
common/mdin_flags/mdin_ewald,mdin_pb,mdin_amoeba,mdin_apbs,sp_apbs
#else
common/mdin_flags/mdin_ewald,mdin_pb,mdin_amoeba
#endif /* APBS */

integer BC_HULP  ! size in integers of common HULP
parameter ( BC_HULP = 10 )

integer     ntpr,ntwr,ntwx,ntwv,ntwe,ntpp,ioutfm,ntwprt,&
            ntwrism, ntave
common/hulp/ntpr,ntwr,ntwx,ntwv,ntwe,ntpp,ioutfm,ntwprt,&
            ntwrism, ntave

!      NMRRDR : Contains information about input/output file redirection
!               REDIR and IREDIR contain information regarding
!               LISTIN, LISTOUT, READNMR, NOESY, SHIFTS, DUMPAVE,
!               PCSHIFT and DIPOLE respectively. If IREDIR(I) > 0,
!               then that input/output has been redirected.

integer iredir(9)
common/nmrrdr/redir,iredir
