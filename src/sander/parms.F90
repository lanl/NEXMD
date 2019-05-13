 !<compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

module parms


!    BC_PARMR is the number of reals in common RPARMS; BC_PARMI is the
!    number of ints in IPARMS.  (Change these if you change the sizes below).

#define BC_PARMR 24560
#define BC_PARMI 1200
#define BC_PARML 1200

!RCW: WARNING - DO NOT DELETE THE COMMON BLOCKS HERE - THEY ARE A HACK TO FORCE
!               THE MEMORY LAYOUT TO BE LINEAR FOR DOING MPI BROADCASTS. ULTIMATELY
!               WE SHOULD GET RID OF THESE, MAKE EVERYTHING DYNAMIC AND DO THE BROADCASTS
!               IN A BETTER WAY.

integer, parameter :: num_bc_parmr = BC_PARMR, num_bc_parmi = BC_PARMI
integer, parameter :: num_bc_parml = BC_PARML
integer, parameter :: MAX_BOND_TYPE = 5000 !NUMBND
integer, parameter :: MAX_ATOM_TYPE = 100  !NATYP

!one_scnb is 1.0d0/scnb as an array of dihedral types since you can now
!set scnb and scee values for individual dihedrals.

_REAL_ rk(MAX_BOND_TYPE),req(MAX_BOND_TYPE),tk(900),teq(900),pk(1200), &
      pn(1200),phase(1200),cn1(1830),cn2(1830), one_scnb(1200), &
      one_scee(1200), &
      solty(MAX_ATOM_TYPE), &
      gamc(1200),gams(1200), &
      asol(200),bsol(200),hbcut(200), &
      cn3(1830),cn4(1830),cn5(1830) ! mjhsieh: for another vdwmodel
common/rparms/rk,req,tk,teq,pk, &
      pn,phase,cn1,cn2,one_scnb, one_scee, solty, &
      gamc,gams,asol,bsol,hbcut, &
      cn3,cn4,cn5 ! mjhsieh: PLEASE, there is no need for common block here.

integer ipn(1200)
common/iparms/ipn

! NPHB is the number of h-bond parameters. NIMPRP is the number of
! improper torsional parameters (NPTRA-NIMPRP is the number of regular
! torsional parameters).


#define BC_PRMLIM 6
integer, parameter :: num_bc_prmlim = BC_PRMLIM
integer       numbnd,numang,nptra,nphb,nimprp,nttyp
common/prmlim/numbnd,numang,nptra,nphb,nimprp,nttyp

end module parms
