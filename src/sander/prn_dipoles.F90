! <compile=optimized>
!+ Routines by Ross Walker for dipole manipulation and printing
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ []
#include "copyright.h"
#include "dprec.fh"
subroutine printdip(numgroups,ipres,coord,q,inddip,mass,natom)
   
   !subroutine to print out current dipoles for atoms in specified groups
   !Example: Called at end of each MD step in a polarization run
   !it will print out the current dipoles for the chosen groups

   !Since the groups will typically be charged so the calculated dipole
   !will be ambiguous. Hence the centre of mass of the specified group
   !is first calculated and then the coordinates are reset such that
   !the dipole moment is calculated with respect to the group centre
   !of mass.

   use constants, only : INV_AMBER_ELECTROSTATIC, AMBER_ELECTROSTATIC, zero
   use qmmm_module, only : qmmm_nml, qm2_struct, qmmm_struct
   implicit none
   integer numgroups
   integer ipres(*)
   integer natom
   _REAL_ coord(3,*)
   _REAL_ q(*)
   _REAL_ inddip(3,*)
   _REAL_ mass(*)

   integer grpcount   
   integer atomcount
   integer qm_temp_count
   _REAL_ centmass(3)
   _REAL_ permdipx, permdipy, permdipz
   _REAL_ inddipx, inddipy, inddipz
   _REAL_ permdiptotal, inddiptotal, diptotal

   ! Convert electron angstrom to Debye
   ! 1 A = 1.0 / 0.529177249 BOHR
   ! 1 Electron BOHR = 2.54158 Debye
   ! Therefore 1 Electron Angstom = 4.802889778 Debye
   ! Note to convert sander internal unit of charge to electrons
   ! divide by 18.2223 (1/18.2223) = 0.05487781455
   _REAL_ elecang_to_debye
   elecang_to_debye = 4.802889778d0

   !if QMMM is in use then use the mulliken charges to calculate dipoles.
   if (qmmm_nml%ifqnt) then
     do qm_temp_count = 1, qmmm_struct%nquant_nlink
         q(qmmm_struct%iqmatoms(qm_temp_count)) = &
                 qm2_struct%scf_mchg(qm_temp_count)*AMBER_ELECTROSTATIC
     end do
   end if

   ! Loop over groups 1 to numgroups - 1
   ! minus 1 since the last group contains everything else
   ! so we don't want it

   !except this is a problem if the only group required 
   !contains all atoms. Therefore check if numgroups=1 and
   !if it does increment to 2
   if (numgroups == 1) then
      numgroups = 2
   end if

   do grpcount = 1,(numgroups-1)
      !Step 1 is to find the centre of mass of the chosen group
      !since the groups can be charged and so we need a common
      !reference frame to obtain a dipole moment that is not
      !ambiguous.

      !Call centre of mass routine - this puts the centre of mass
      !of group 'grpcount' in centmass(3), in order X,Y,Z
      call calc_grp_centre_of_mass(grpcount, centmass, ipres, coord, &
                                   mass, natom, q)

      write(6,'(37x,''x'',9x,''y'',9x,''z'',7x,''Total'')')
      write(6,'(''DIPGRP -'',i3, '': centre  of  mass :'',3f10.3)') &
           grpcount, centmass(1), centmass(2), centmass(3)

      !Step 2 is to find the permanent and induced dipoles for the group
      permdipx = 0.0d0
      permdipy = 0.0d0
      permdipz = 0.0d0
      inddipx = 0.0d0
      inddipy = 0.0d0
      inddipz = 0.0d0
      
      do atomcount = 1,natom
         if (ipres(atomcount) == grpcount) then
            permdipx = permdipx + (q(atomcount)*(coord(1,atomcount)-centmass(1)))
            permdipy = permdipy + (q(atomcount)*(coord(2,atomcount)-centmass(2)))
            permdipz = permdipz + (q(atomcount)*(coord(3,atomcount)-centmass(3)))
            inddipx = inddipx + inddip(1,atomcount) 
            inddipy = inddipy + inddip(2,atomcount)
            inddipz = inddipz + inddip(3,atomcount)
         end if
      end do    !  atomcount = 1,natom
      permdipx = permdipx * elecang_to_debye * INV_AMBER_ELECTROSTATIC
      permdipy = permdipy * elecang_to_debye * INV_AMBER_ELECTROSTATIC
      permdipz = permdipz * elecang_to_debye * INV_AMBER_ELECTROSTATIC
      inddipx = inddipx * elecang_to_debye * INV_AMBER_ELECTROSTATIC
      inddipy = inddipy * elecang_to_debye * INV_AMBER_ELECTROSTATIC
      inddipz = inddipz * elecang_to_debye * INV_AMBER_ELECTROSTATIC

      ! Calculate totals
      ! Permanent dipole total = Sqrt(px^2+py^2+pz^2)
      ! Induced dipole total = Sqrt(ix^2+iy^2+iz^2)
      ! Total Dipole = Sqrt([px+ix]^2+[py+iy]^2+[pz+iz]^2)
      permdiptotal = sqrt((permdipx*permdipx)+(permdipy*permdipy) &
                         +(permdipz*permdipz))
      inddiptotal  = sqrt((inddipx*inddipx)+(inddipy*inddipy) &
                         +(inddipz*inddipz))
      diptotal     = sqrt((permdipx+inddipx)**2+(permdipy+inddipy)**2 &
                         +(permdipz+inddipz)**2)

     !Now we write out the results for this group
      write(6,'(''DIPGRP -'',i3, '': permanent dipole :'',4f10.3,'' D'')') &
           grpcount, permdipx, permdipy, permdipz, permdiptotal
      write(6,'(''DIPGRP -'',i3, '': induced   dipole :'',4f10.3,'' D'')') &
           grpcount, inddipx, inddipy, inddipz, inddiptotal
      write(6,'(''DIPGRP -'',i3, '': total     dipole :'',4f10.3,'' D'')') &
           grpcount, permdipx+inddipx, permdipy+inddipy, permdipz+inddipz, &
           diptotal
     !Move to next group
   end do  !  grpcount = 1,(numgroups-1)

   if (qmmm_nml%ifqnt) then
     do qm_temp_count = 1, qmmm_struct%nquant_nlink
        q(qmmm_struct%iqmatoms(qm_temp_count)) = zero
     end do
   end if

   return

end subroutine printdip

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Finds centre of mass of a chosen group]
subroutine calc_grp_centre_of_mass(grpcount,centmass,ipres,coord,mass,natom,q)

!Calculates the centre of mass for all atoms in the group 'grpcount'
!Using the following:
!      1 N
! Cx = -Sum(m[i]x[i])
!      M i
!
!For X, Y and Z independently
!

   implicit none
   integer grpcount
   integer natom
   integer ipres(*)
   _REAL_ coord(3,*)
   _REAL_ mass(*)
   _REAL_ centmass(*)
   _REAL_ q(*)
   _REAL_ masssum

   integer atomcount

   centmass(1) = 0.0d0
   centmass(2) = 0.0d0
   centmass(3) = 0.0d0
   masssum = 0.0d0

   do atomcount = 1,natom
      if (ipres(atomcount) == grpcount) then

         ! Atom is part of group 'grpcount'
         ! add it to the centre of mass array 'centmass(1 to 3)'
         centmass(1) = centmass(1) + (mass(atomcount)*coord(1,atomcount))
         centmass(2) = centmass(2) + (mass(atomcount)*coord(2,atomcount))
         centmass(3) = centmass(3) + (mass(atomcount)*coord(3,atomcount))
         masssum = masssum + mass(atomcount)

      end if
   end do  !  atomcount = 1, natom

   ! Divide sum of mass(i)*coord(i) by total mass of group
   centmass(1) = centmass(1) / masssum
   centmass(2) = centmass(2) / masssum
   centmass(3) = centmass(3) / masssum

   return

end subroutine calc_grp_centre_of_mass

