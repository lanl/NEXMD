#include "dprec.fh"
! This module sets up pointer arrays for allocated segments of the LOCMEM
! stack arrays. Arrays from the PRMTOP have names matching the %FLAG
! names.
module memory_module
   implicit none
   public

#include "memory.h"

! Old stack arrays:
   _REAL_,  dimension(:), allocatable, target :: x
   integer, dimension(:), allocatable, target :: ix
   character(len=4), dimension(:), allocatable, target :: ih

! New allocatable arrays, using pointers as a transitional form:
   character(len=4), dimension(:), pointer :: atom_name, residue_label, &
      amber_atom_type, tree_chain_classification

   integer, dimension(:), pointer :: &
      bonds_inc_hydrogen_1, &
      bonds_inc_hydrogen_2, &
      bonds_inc_hydrogen_type, &
      bonds_without_hydrogen_1, &
      bonds_without_hydrogen_2, &
      bonds_without_hydrogen_type, &
      angles_inc_hydrogen_1, &
      angles_inc_hydrogen_2, &
      angles_inc_hydrogen_3, &
      angles_inc_hydrogen_type, &
      angles_without_hydrogen_1, &
      angles_without_hydrogen_2, &
      angles_without_hydrogen_3, &
      angles_without_hydrogen_type, &
      dihedrals_inc_hydrogen_1, &
      dihedrals_inc_hydrogen_2, &
      dihedrals_inc_hydrogen_3, &
      dihedrals_inc_hydrogen_4, &
      dihedrals_inc_hydrogen_type, &
      dihedrals_without_hydrogen_1, &
      dihedrals_without_hydrogen_2, &
      dihedrals_without_hydrogen_3, &
      dihedrals_without_hydrogen_4, &
      dihedrals_without_hydrogen_type


   integer, dimension(:), pointer :: &
      atom_type_index, number_excluded_atoms, &
      nonbonded_parm_index, residue_pointer, &
      excluded_atoms_list, join_array, atoms_per_molecule, &
      restraint_group, tgt_fit_group, tgt_rms_group, &
      belly_group, atom_noshake, num_bonds, nmr_iwork

   _REAL_, dimension(:), pointer :: charge, massinv, mass, &
      bond_bcoef, hbcut, box_dimensions, radii, screen, &
      polarizability, nmr_work, reff, onereff, group_weight, &
      gb_vdw_radii, gb_p1, gb_p2, gb_p3, gb_p4, &
!      rborn_max, rborn_min, rborn_ave, rborn_fluct
! Modified by WJM, YD
      rborn_max, rborn_min, rborn_ave, rborn_fluct, dampfactor, pol2
!!

   integer, dimension(:), pointer :: egb_neighbor_list

   _REAL_, dimension(:,:), pointer :: ref_coordinate, coordinate, &
      velocity, velocity_old, force

   integer, pointer, dimension(:) :: nmr_imet
   _REAL_, pointer, dimension(:) :: &
       nmr_xmet, nmr_force, nmr_ddep, nmr_dddep, &
       nmr_dorat, nmr_ddrat, nmr_rate, nmr_trp, nmr_dint

   _REAL_, pointer, dimension(:) :: cnstph_dcharge

contains
   ! This routine intitializes all of the pointer arrays after they
   ! have been allocated by LOCMEM. Currently, some of the arrays
   ! that have different sizes based on input options are just
   ! set to the maximum size, but are valid for the range that
   ! was allocated by LOCMEM.
   subroutine memory_init()
      implicit none
#include "md.h"
#include "nmr.h"
#include "dynph.h"

      residue_label => ih(m02:m02+nres-1)
      atom_name => ih(m04:m04+natom-1)
      amber_atom_type => ih(m06:m06+natom-1)
      tree_chain_classification => ih(m08:m08+natom-1)

  ! FIXME: n12, m14, m16   DEAD but still allocated??

      atom_type_index => ix(i04:i04+natom-1)
      number_excluded_atoms => ix(i08:i08+natom-1)
      nonbonded_parm_index => ix(i06:i06+ntypes**2-1)
      residue_pointer => ix(i02:i02+nres)

      ! Should these atom index arrays be converted
      ! from multiple arrays to multi-dimensional arrays?

      bonds_inc_hydrogen_1 => ix(iibh:iibh+nbonh-1)
      bonds_inc_hydrogen_2 => ix(ijbh:ijbh+nbonh-1)
      bonds_inc_hydrogen_type => ix(iicbh:iicbh+nbonh-1)

      bonds_without_hydrogen_1 => ix(iiba:iiba+nbona-1)
      bonds_without_hydrogen_2 => ix(ijba:ijba+nbona-1)
      bonds_without_hydrogen_type => ix(iicba:iicba+nbona-1)

      angles_inc_hydrogen_1 => ix(i24:i24+ntheth-1)
      angles_inc_hydrogen_2 => ix(i26:i26+ntheth-1)
      angles_inc_hydrogen_3 => ix(i28:i28+ntheth-1)
      angles_inc_hydrogen_type => ix(i30:i30+ntheth-1)

      angles_without_hydrogen_1 => ix(i32:i32+ntheta-1)
      angles_without_hydrogen_2 => ix(i34:i34+ntheta-1)
      angles_without_hydrogen_3 => ix(i36:i36+ntheta-1)
      angles_without_hydrogen_type => ix(i38:i38+ntheta-1)

      dihedrals_inc_hydrogen_1 => ix(i40:i40+nphih-1)
      dihedrals_inc_hydrogen_2 => ix(i42:i42+nphih-1)
      dihedrals_inc_hydrogen_3 => ix(i44:i44+nphih-1)
      dihedrals_inc_hydrogen_4 => ix(i46:i46+nphih-1)
      dihedrals_inc_hydrogen_type => ix(i48:i48+nphih-1)

      dihedrals_without_hydrogen_1 => ix(i50:i50+nphia-1)
      dihedrals_without_hydrogen_2 => ix(i52:i52+nphia-1)
      dihedrals_without_hydrogen_3 => ix(i54:i54+nphia-1)
      dihedrals_without_hydrogen_4 => ix(i56:i56+nphia-1)
      dihedrals_without_hydrogen_type => ix(i58:i58+nphia-1)

      excluded_atoms_list => ix(i10:i10+nnb-1) !FIXME: locmem says nnb*2
      join_array => ix(i64:i64+natom-1)

      atoms_per_molecule => ix(i70:i70+nspm-1)

      ! i78: UNUSED

      num_bonds => ix(i80:i80+natom-1)
      egb_neighbor_list => ix(i82:i82+natom*merge(80,40,gbsa==2))

      charge => x(l15:l15+natom-1)
      massinv => x(lwinv:lwinv+natom-1)
      mass => x(lmass:lmass+natom-1)
      radii => x(l97:l97+natom-1)
      screen => x(l96:l96+natom-1)
      polarizability => x(lpol:lpol+natom-1)
! by WJM, YD
      pol2 => x(lpol2:lpol2+natom-1)
      dampfactor => x(ldf:ldf+natom-1)
!!

      restraint_group => ix(icnstrgp:icnstrgp+ndper-1)
      tgt_fit_group => ix(itgtfitgp:itgtfitgp+natom-1)
      tgt_rms_group => ix(itgtrmsgp:itgtrmsgp+natom-1)
      belly_group => ix(ibelly:ibelly+natom-1)
      atom_noshake => ix(noshake:noshake+natom-1)

      call set_rank2_pointer(coordinate,x(lcrd),3,natom)
      call set_rank2_pointer(ref_coordinate,x(lcrdr),3,natom)
      call set_rank2_pointer(velocity,x(lvel),3,natom) ! 6* when imin/=0 ???
      call set_rank2_pointer(velocity_old,x(lvel2),3,natom) ! 6* when imin/=0 ???
      call set_rank2_pointer(force,x(lforce),3,natom)

      !coor_ref?   x(l45:l45+3*natom*am_nbead+mxvar-1)
      group_weight => x(l60:l60+natom-1)

      ! l65: polarization  DEAD??

      nmr_work => x(lnmr01:lnmr01+irlreq-1)
      nmr_iwork => ix(inmr02:inmr02+intreq-1)

#ifdef LES
      reff => x(l98:l98+natom*ncopy-1)
      onereff => x(l99:l99+natom*ncopy-1)
#else
      reff => x(l98:l98+natom-1)
      onereff => x(l99:l99+natom-1)
#endif

      gb_vdw_radii => x(l165:l165+natom-1)
      gb_p1 => x(l170:l170+natom-1)
      gb_p2 => x(l175:l175+natom-1)
      gb_p3 => x(l180:l180+natom-1)
      gb_p4 => x(l185:l185+natom-1)
      rborn_max => x(l186:l186+natom-1)
      rborn_min => x(l187:l187+natom-1)
      rborn_ave => x(l188:l188+natom-1)
      rborn_fluct => x(l189:l189+natom-1)


      if (nmropt >= 2) then
         nmr_imet => ix(i65:i65+mxsub*isubi-1)
         nmr_xmet => x(l105:l105+mxsub*isubr-1)
         nmr_force => x(l110:l110+3*natom+mxvar-1)
         nmr_ddep => x(l115:l115+ma*ma-1)
         nmr_dddep => x(l120:l120+3*ma*ma-1)
         nmr_dorat => x(l125:l125+3*ma*ma-1)
         nmr_ddrat => x(l130:l130+3*ma*ma-1)
         nmr_rate => x(l135:l135+ma*ma-1)
         nmr_trp => x(l140:l140+ma*ma-1)
         nmr_dint => x(l145:l145+3*ma+mxvar-1)
      end if

      if ( icnstph /= 0 ) then
         cnstph_dcharge => x(l190:l190+natom-1)
      end if

      return
   end subroutine memory_init

   subroutine set_rank2_pointer(ptr,array,dim1,dim2)
      _REAL_, pointer :: ptr(:,:)
      integer, intent(in) :: dim1, dim2
      _REAL_, target :: array(dim1,dim2)
      ptr => array 
   end subroutine set_rank2_pointer

end module memory_module
