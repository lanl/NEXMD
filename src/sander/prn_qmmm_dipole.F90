! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

subroutine qmmm_dipole(coord,mass,ipres,lbres,nres)
 use qmmm_module, only : qm2_params, qm2_struct, qmmm_struct, qmmm_nml
 use constants, only : light_speed, bohr_radius, charge_on_elec

! use findmask
      implicit none
!
! Calculates the quantum dipole moment by using NDDO mehodology
! from the atomic charges and the lone-pairs 

      integer, intent(in) :: nres, ipres(*)
      _REAL_, intent(inout) :: mass(*)
      _REAL_, intent(inout) :: coord(*)
      character(len=4), intent(in) :: lbres(*)
      
! Determination of constants are  needed for the computation of the dipole

      integer :: imm,i,j,iqm,ilm,atomdipo

! Number of atoms without solvent

      _REAL_ :: mmdipole(3), totaldipol
      _REAL_ :: centermass(3)
      _REAL_ :: qm_cm_coord(3*qmmm_struct%nquant_nlink)
      _REAL_ :: axis(3), densisp, summasp
      _REAL_ :: sc_const


      sc_const=light_speed*charge_on_elec

      centermass(1:3)=0.0d0

! Computation dipole moment of atoms

      mmdipole(1:3)=0.0d0
      
      qm_cm_coord=0.0d0

         
      select case (qmmm_nml%printdipole)
         case (1)
            atomdipo = qmmm_struct%nquant
            ! Calculation of center of mass
            call calc_center_mass(atomdipo,coord,mass,centermass)       
            !Building array of coordinates to compute the QM dipole moment
            do i=1,qmmm_struct%nquant
                imm = qmmm_struct%iqmatoms(i)
               do j=1,3
                  qm_cm_coord((i-1)*3+j)=coord((imm-1)*3+j)-centermass(j)
               end do
            end do 

            iqm=qmmm_struct%nquant
         
            if (qmmm_struct%nlink.ne.0) then
               do i=1,qmmm_struct%nlink
                  iqm=iqm+1
                  call positionla(i,coord,axis)
                  do j=1,3
                     qm_cm_coord((iqm-1)*3+j)=axis(j)-centermass(j)
                  end do
               end do
            end if

            call qm2_calc_dipole(qm2ds,qmmm_struct,qm_cm_coord)
 
         ! Asking if MM + QM dipole moment is required
         case (2)
            
            ! Stripping the solvent assuming that it is WAT
            call strip_solvent_atoms(atomdipo,ipres,lbres,nres)

            ! Calculation of center of mass
            call calc_center_mass(atomdipo,coord,mass,centermass)
            iqm=0
            ilm=0
            do imm=1,atomdipo

               ! Asking if the atom is quantum atom 
               if(.not.qmmm_struct%atom_mask(imm)) then

                  ! Asking if the classical atom is a link classical atom 
                  if(.not.qmmm_struct%mm_link_mask(imm)) then
                     do i=1,3 
                        mmdipole(i)=mmdipole(i)+sc_const*(coord((imm-1)*3+i)-centermass(i))*&
                        qmmm_struct%scaled_mm_charges(imm)*1.0d-10
                     end do
                  else
                     ilm=ilm+1
                     do i=1,3 
                        mmdipole(i)=mmdipole(i)+sc_const*(coord((imm-1)*3+i)-centermass(i))*&
                        qmmm_struct%mm_link_pair_resp_charges(ilm)/18.2223*1.0d-10
                     end do
                  end if
               else
                  !Building the coordinate array for QM atoms
                  iqm=iqm+1
                  do i=1,3
                     qm_cm_coord((iqm-1)*3+i)=coord((imm-1)*3+i)-centermass(i)
                  end do 
               end if
            end do
            !Quantum Link atom
            if (qmmm_struct%nlink.ne.0) then
               do i=1,qmmm_struct%nlink
                  iqm=iqm+1
                  call positionla(i,coord,axis)
                  do j=1,3
                     qm_cm_coord((iqm-1)*3+j)=axis(j)-centermass(j)
                  end do
               end do
            end if
            
            !Calculation of MM dipole in Debye
            
            mmdipole(1:3)=mmdipole(1:3)/1.0d-21
            totaldipol=sqrt(mmdipole(1)**2+mmdipole(2)**2+mmdipole(3)**2)

            write(6,'(" ","          ","       X    ","    Y    ","    Z    "," TOTAL  ")')
            write(6,'(" "," MM DIPOLE ",F9.3,F9.3,F9.3,F9.3)'), mmdipole(1:3), totaldipol
            
            !Computing the QM dipole moment
            call qm2_calc_dipole(qm2ds,qmmm_struct, qm_cm_coord)

          case default
 
       end select 
      
       return
end subroutine qmmm_dipole

subroutine strip_solvent_atoms(atomdipo,ipres,lbres,nres)

! This subroutine strip water molecules
 
    integer atomdipo,i,j,nres,ipres(*)
    character(len=4) lbres(*)

    atomdipo=0
      do i=1,nres
         if (lbres(i).ne."WAT") then
            do j=ipres(i),ipres(i+1)-1
               atomdipo=atomdipo+1
            enddo
         endif
      enddo
    return
end subroutine strip_solvent_atoms

subroutine calc_center_mass(atomdipo,coord,mass,centermass)

! This subroutine compute the center of mass over both parts the quantum and the 
! classical part.

 use qmmm_module, only : qmmm_struct, qmmm_nml
 implicit none
   integer atomdipo
   _REAL_ mass(*)
   _REAL_ coord(3*atomdipo)
   _REAL_ centermass(3)
   _REAL_ masssum, axis(3)

   integer j,i,iqm

   centermass(1) = 0.0d0
   centermass(2) = 0.0d0
   centermass(3) = 0.0d0
   masssum = 0.0d0

   do j=1,atomdipo
      if (qmmm_nml%printdipole == 1) then
         iqm=qmmm_struct%iqmatoms(j)
         centermass(1) = centermass(1) + (mass(iqm)*coord((iqm-1)*3+1))
         centermass(2) = centermass(2) + (mass(iqm)*coord((iqm-1)*3+2))
         centermass(3) = centermass(3) + (mass(iqm)*coord((iqm)*3))
         masssum = masssum + mass(iqm)
      else
         centermass(1) = centermass(1) + (mass(j)*coord((j-1)*3+1))
         centermass(2) = centermass(2) + (mass(j)*coord((j-1)*3+2))
         centermass(3) = centermass(3) + (mass(j)*coord((j)*3))
         masssum = masssum + mass(j)
      end if 
   end do

! The link atom must be Hydrogen 

   if (qmmm_struct%nlink.ne.0) then 
       do i=1,qmmm_struct%nlink
           call positionla(i,coord,axis)
           centermass(1) = centermass(1) + (1.007d0*axis(1))
           centermass(2) = centermass(2) + (1.007d0*axis(2))
           centermass(3) = centermass(3) + (1.007d0*axis(3))
           masssum = masssum + 1.007d0
       end do
   end if

   centermass(1) = centermass(1) / masssum
   centermass(2) = centermass(2) / masssum
   centermass(3) = centermass(3) / masssum

   return

end subroutine calc_center_mass

subroutine positionla(linkatom,coord,axis)
! This subroutine compute the coordinates for the link atom
! in the crd system of coordinates

 use qmmm_module, only : qmmm_struct

   _REAL_ coord(*),axis(3)
   integer it,iq,linkatom,linkpos,x,i

   linkpos=qmmm_struct%nquant+linkatom
   it = qmmm_struct%link_pairs(2,linkatom)
   iq = qmmm_struct%iqmatoms(it)
   axis=0.0d0
   do i=1,3
      axis(i)=coord((iq-1)*3+i)+(qmmm_struct%qm_coords(i,it)-qmmm_struct%qm_coords(i,linkpos))
   end do

   return
end subroutine 
