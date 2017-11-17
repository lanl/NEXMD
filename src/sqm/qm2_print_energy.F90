#include "copyright.h"
#include "dprec.fh"

subroutine qm2_print_energy(qm2_struct, verbosity, qmtheory, escf, qmmm_struct)

  ! Print QM energy contributions
  ! For historic reasons the nomenclature for the energy terms is misleading:
  ! SCF Energy:
  !     This is the total heat of formation
  !     It contains the Total Energy from below and plus 
  !     the dispersion and hydrogen bond correction (if in use) plus
  !     the additional semiempirical term, derived from experiment, 
  !     to obtain the heat of formation
  ! Total Energy:
  !     This is the total SCF energy, that is, the electronic energy
  !     plus the core repulsion energy
  ! Dispersion Energy:
  !     Empirical dispersion correction
  ! H-Bond energy:
  !     Empirical correction for hydrogen bonding

  use constants, only : J_PER_CAL, EV_TO_KCAL, KCAL_TO_EV
  use qmmm_struct_module, only : qmmm_struct_type
  use qmmm_module, only: qmmm_scratch,qm2_structure
  use qmmm_qmtheorymodule, only : qmTheoryType

   use cosmo_C, only: ediel,onsagE,solvent_model,potential_type
   use qm2_davidson_module

  implicit none
  
  type(qm2_structure),intent(inout) :: qm2_struct
  integer, intent(in) :: verbosity
  _REAL_, intent(in) :: escf
  type(qmmm_struct_type), intent(in) :: qmmm_struct
  type(qmTheoryType) :: qmtheory

  _REAL_ :: total_energy

!Printing molecular orbital energies
	if (verbosity > 4) then
		write (6,'("QMMM: Occupied MO Energies (eV):")')
		write (6,'(5f18.8)') qmmm_scratch%mat_diag_workspace(1:qm2_struct%nclosed,1)
		write (6,'("QMMM: Virtual MO Energies (eV):")') 
		write(6,'(5f18.8)') qmmm_scratch%mat_diag_workspace(qm2_struct%nclosed+1:qm2_struct%norbs,1)
	endif

!Printing other energies
  if (verbosity > 0) then
     !Verbosity level of 1 or more = print more accurate SCF energy
     write (6,'("QMMM:")')
     write (6,'("QMMM: SCF Energy =",f18.8," KCal/mol, ",f18.8," KJ/mol")') escf, escf*J_PER_CAL
     write (6,'("QMMM: SCF Energy = Heat of formation")')
     !If verbosity level is greater than 1 we also print the nuclear and electronic energies.
     if (verbosity > 1) then
        write (6,'("QMMM:")')
        write (6,'("QMMM:        Electronic energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
             qmmm_struct%elec_eng, qmmm_struct%elec_eng*EV_TO_KCAL
        if ((solvent_model.gt.0).and.(potential_type.eq.3)) then   ! Dielectric Permittivity from COSMO module
           write (6,'("QMMM: Dielectric(COSMO) energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
              ediel,ediel*EV_TO_KCAL
        elseif ((solvent_model.gt.0).and.(potential_type.eq.2)) then !Onsager Solvent Model
           write (6,'("QMMM: Dielectric(Onsager) energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
              onsagE,onsagE*EV_TO_KCAL
        end if

        if  ( qmtheory%DFTB ) then
           write (6,'("QMMM:         Repulsive energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                qmmm_struct%enuclr_qmqm,qmmm_struct%enuclr_qmqm*EV_TO_KCAL
           !! write (6,'("QMMM:        Careful: Dispersion Energy Already Included.")')
           !! write (6,'("QMMM:        Dispersion Energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
           !!           dftb_edisp*AU_TO_KCAL
           total_energy = qmmm_struct%elec_eng + qmmm_struct%enuclr_qmqm
        else
           write (6,'("QMMM: QM core - QM core energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                qmmm_struct%enuclr_qmqm,qmmm_struct%enuclr_qmqm*EV_TO_KCAL
           write (6,'("QMMM: QM core - MM atom energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                qmmm_struct%enuclr_qmmm,qmmm_struct%enuclr_qmmm*EV_TO_KCAL
           write (6,'("QMMM: Total core - core energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                qmmm_struct%enuclr_qmmm+qmmm_struct%enuclr_qmqm, &
                (qmmm_struct%enuclr_qmmm+qmmm_struct%enuclr_qmqm)*EV_TO_KCAL
           total_energy = qmmm_struct%elec_eng + qmmm_struct%enuclr_qmmm + qmmm_struct%enuclr_qmqm
        end if
        write (6,'("QMMM:             Total energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
             total_energy, total_energy*EV_TO_KCAL
        ! Print Dispersion energy if in use
        if (qmtheory%DISPERSION .or. qmtheory%DISPERSION_HYDROGENPLUS) then
           write (6,'("QMMM:")')
           write (6,'("QMMM:        Dispersion energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                qmmm_struct%dCorrection*KCAL_TO_EV, qmmm_struct%dCorrection
        end if
        ! Print H-bond energy if in use
        if (qmtheory%DISPERSION_HYDROGENPLUS) then
           write (6,'("QMMM:            H-bond energy = ",f18.8," eV (",f18.8," KCal/mol)")') &
                qmmm_struct%hCorrection*KCAL_TO_EV, qmmm_struct%hCorrection
        end if
     end if
  end if

end subroutine qm2_print_energy
