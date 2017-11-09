! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

subroutine qm2_calc_dipole(qmmm_struct, coord)
 use qmmm_module, only : qm2_params, qm2_struct,  qmmm_nml
 use constants, only : light_speed, bohr_radius, charge_on_elec, CODATA08_AU_TO_DEBYE, BOHRS_TO_A
 use qm2_davidson_module
 use qmmm_struct_module, only : qmmm_struct_type

! use findmask
      implicit none
 
      type(qmmm_struct_type), intent(inout) :: qmmm_struct

!
! Calculates the quantum dipole moment by using NDDO mehodology
! from the atomic charges and the lone-pairs 

      _REAL_, intent(inout) :: coord(*)
      
! Determination of constants are  needed for the computation of the dipole

      integer :: loop_count,orb_beg,orb_end,orb_size,i,nquant,j

! Number of atoms without solvent

      _REAL_ :: qmdipole(3), totaldipol
      _REAL_ :: densisp, summasp
      integer :: ndubl

      _REAL_ :: sc_const, c2_const, iqm_charge

	! CML TEST
	_REAL_ :: dip(3, qm2ds%Nb*(qm2ds%Nb+1)/2)

      ! AWG: Currently not supported with d orbitals for NDDO methods
      if ( .not. qmmm_nml%qmtheory%DFTB ) then
         do i = 1,qmmm_struct%nquant_nlink
            if ( qm2_params%natomic_orbs(i) > 5 ) return
         end do
      end if

      sc_const=light_speed*charge_on_elec

! Computation dipole moment of quatum atoms

      qmdipole(1:3)=0.0d0

      ndubl=qm2_struct%nclosed

      ! Increasing the number of quantum atoms and determining orbital numbers
      do nquant=1,qmmm_struct%nquant_nlink

         ! Asking for the Method, if DFTB is used it does not compute the changes in charge due to orbitals.             
         if (.not. qmmm_nml%qmtheory%DFTB) then
            orb_beg=qm2_params%orb_loc(1,nquant)
            orb_end=qm2_params%orb_loc(2,nquant)

            ! Determination of the values difference between the orbitals
            orb_size=orb_end-orb_beg

            ! Trying to get de density value for P and S interaction over the same atom 

           c2_const=2.0D0*bohr_radius*sc_const*qm2_params%multip_2c_elec_params(1,nquant)


            if (orb_size>0) then 
               i=0
               ! do loop_count=orb_beg+1,orb_end
               do loop_count=orb_beg+1,orb_beg+3 ! only for p-orbitals
                  i=i+1
                  densisp=0.0d0
                  do j=1,ndubl
                    densisp = densisp+(2.0d0*qm2_struct%eigen_vectors(orb_beg,j)*qm2_struct%eigen_vectors(loop_count,j))
                  end do
                   qmdipole(i)=qmdipole(i)-c2_const*densisp;
               end do

            end if 
 
         end if

         ! Getting the mulliken charge 
         iqm_charge=qm2_struct%scf_mchg(nquant)
	 !iqm_charge=qm2_params%core_chg(nquant) !FOR TESTING PURPUSES ONLY, COUNTING CORES CHARGES
         do i=1,3
            qmdipole(i)=qmdipole(i)+coord((nquant-1)*3+i)*iqm_charge*sc_const*1.0d-10
         end do
      end do

!Calculation of the dipole in Debye
            
      qmdipole(1:3)=qmdipole(1:3)/1.0d-21
      totaldipol=sqrt(qmdipole(1)**2+qmdipole(2)**2+qmdipole(3)**2)

      write(6,'(" ","          ","       X    ","    Y    ","    Z    "," TOTAL  ")')
      write(6,'(" "," QM DIPOLE (Debye)",4F9.3)') &
            qmdipole(1), qmdipole(2), qmdipole(3), totaldipol
      !write(6,*) &
      !      qmdipole(1), qmdipole(2), qmdipole(3), totaldipol

! 	  ! CML TEST
 		!write(6,*) 'sc_const', sc_const
 		!write(6,*) 'c2_const', c2_const
 		!write(6,*) 'AU_TO_DEBYE', CODATA08_AU_TO_DEBYE
! 		call unpacking(qm2ds%Nb,qm2_struct%den_matrix,qm2ds%eta_scratch_2,'s')
! 		qm2ds%eta_scratch_2 = qm2ds%eta_scratch_2 * 2.d0
! 		call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%eta_scratch_2,qm2ds%eta_scratch,qm2ds%xi_scratch)
! 		call get_dipole_matrix(qmmm_struct%qm_coords, dip)
! 		dip = dip * (light_speed * charge_on_elec) * 1.d11   ! 4.80... D / 1 eA
! 		call unpacking(qm2ds%Nb,dip(1,:),qm2ds%xi_scratch_2,'s')
! 		call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%xi_scratch_2,qm2ds%xi_scratch,qm2ds%eta_scratch_2)
! 		densisp = 0.d0
! 		do i = 1,qm2ds%Nb
! 			do j = 1,qm2ds%Nb
! 				densisp = densisp + 2.d0*qm2ds%eta_scratch((i-1)*qm2ds%Nb+j)*qm2ds%xi_scratch((j-1)*qm2ds%Nb+i)
! 			end do
! 		end do
! 
! 		write(6,*) 'X dipole', densisp
! 		flush(6)
! 		qmdipole(1) = densisp
! 		densisp = 0.d0
! 		call unpacking(qm2ds%Nb,dip(2,:),qm2ds%xi_scratch_2,'s')
! 		call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%xi_scratch_2,qm2ds%xi_scratch,qm2ds%eta_scratch_2)
! 		do i = 1,qm2ds%Nb
! 			do j = 1,qm2ds%Nb
! 				densisp = densisp + 2.d0*qm2ds%eta_scratch((i-1)*qm2ds%Nb+j)*qm2ds%xi_scratch((j-1)*qm2ds%Nb+i)
! 			end do
! 		end do
! 
! 		write(6,*) 'Y dipole', densisp
! 		flush(6)
! 		qmdipole(2) = densisp
! 		densisp = 0.d0
! 		call unpacking(qm2ds%Nb,dip(3,:),qm2ds%xi_scratch_2,'s')
! 		call mo2sitef(qm2ds%Nb,qm2ds%vhf,qm2ds%xi_scratch_2,qm2ds%xi_scratch,qm2ds%eta_scratch_2)
! 		do i = 1,qm2ds%Nb
! 			do j = 1,qm2ds%Nb
! 				densisp = densisp + 2.d0*qm2ds%eta_scratch((i-1)*qm2ds%Nb+j)*qm2ds%xi_scratch((j-1)*qm2ds%Nb+i)
! 			end do
! 		end do
! 
! 		write(6,*) 'Z dipole', densisp
! 		flush(6)
! 		qmdipole(3) = densisp
      return
end subroutine qm2_calc_dipole

