! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "assert.fh"
subroutine qm2_calc_rij_and_eqns(qm2_params,qmmm_nml,qmmm_mpi,qm2_rij_eqns,qmmm_struct, coords, nquant_nlink, crdsmm, natom, npairs)

!-----------------------------------------------------------------
! Written by Ross Walker (TSRI, 2005)
! 
! This routine should be called on each call to QM_MM. It's purpose
! is to calculate RIJ for each QM-QM pair. It stores this is in
! the qm2_rij_eqns structure. It also calculated a number of RIJ
! related equations that involve exponentials, sqrts etc. In this
! way they only need to be calculated once rather than several times
! during the energy calculation and then in the derivative code.
!
! AWG: Note this is currently only done for atom pairs with at
! maximum p orbitals since different subroutines are used to
! compute integrals involving semiempirical d orbitals.
!
!-----------------------------------------------------------------

  use qmmm_module, only : qm2_rij_eqns_structure, &
                          qmmm_mpi_structure,alph_mm
  use constants, only : A_TO_BOHRS, A2_TO_BOHRS2, one, two
  use qmmm_struct_module, only : qmmm_struct_type
  use qm2_params_module,  only : qm2_params_type
  use qmmm_nml_module   , only : qmmm_nml_type

  implicit none

! Passed in
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qm2_rij_eqns_structure),intent(inout) :: qm2_rij_eqns
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   integer, intent(in) :: nquant_nlink, natom, npairs
   _REAL_, intent(in) :: coords(3,nquant_nlink), crdsmm(4,npairs)
                         !crdsmm array is laid out as x,y,z,chg,x,y,z,chg...
! Local
   integer :: i,j, loop_count
   logical :: sp_atom, spd_atom
   _REAL_ :: r2, rr2, rr, vec(3), onerij, rij
   _REAL_ :: qmi_oneBDD1, qmi_oneBDD2, qmi_oneBDD3
   _REAL_ :: qmj_oneBDD1, pddge1i, pddge2i,qmi_alpa
   _REAL_ :: ijBDD1, qmi_DD, qmi_QQ, qmi_QQ2
   _REAL_ :: RRADD, RRMDD, RRAQQ, RRMQQ
   _REAL_ :: RIJ_temp, qmx(3)
   integer :: ier=0
   integer :: num_per_thread, jstart, jend

#include "qm2_array_locations.h"
 
   !Only full QM-MM interaction (multipole) benefits from qmmmrij_incore
   if (qmmm_nml%qmmmrij_incore) then
     !Allocate memory here for the qmmmrijdata array
     !If this is the first call it will be allocated as either
     !min(npairs*qmmm_struct%nquant+npairs,qmmm_struct%nquant * (natom - qmmm_struct%nquant))
     !Since the number of pairs can change on each call we should check each time to see
     !if we need to re-allocate.
     num_per_thread = qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1
     if (qmmm_struct%qm2_calc_rij_eqns_first_call) then
       call qm2_allocate_qm2_qmmm_rij_eqns(qmmm_nml,qm2_rij_eqns, qmmm_mpi, qmmm_struct, natom,npairs)
       qmmm_struct%qm2_calc_rij_eqns_first_call = .false.
     else if (min(npairs*(num_per_thread)+npairs, &
              (num_per_thread) * (natom - qmmm_struct%nquant)) > qm2_rij_eqns%qmmmrij_allocated) then
       !We need to de-allocate the array and then call the allocation routine again so it gets
       !allocated larger.
       deallocate(qm2_rij_eqns%qmmmrijdata,stat=ier)
       REQUIRE(ier==0)
       call qm2_allocate_qm2_qmmm_rij_eqns(qmmm_nml,qm2_rij_eqns, qmmm_mpi,qmmm_struct, natom,npairs)
     end if

     loop_count = 0

     do i=qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end 

       qmi_alpa = qm2_params%cc_exp_params(i)

       sp_atom     = (qm2_params%natomic_orbs(i) == 4)
       spd_atom    = (qm2_params%natomic_orbs(i) == 9)
       qmi_DD      = qm2_params%multip_2c_elec_params(1,i)
       qmi_QQ      = qm2_params%multip_2c_elec_params(2,i)*two
       qmi_oneBDD1 = qm2_params%multip_2c_elec_params(3,i)**2
       qmi_oneBDD2 = qm2_params%multip_2c_elec_params(4,i)**2
       qmi_oneBDD3 = qm2_params%multip_2c_elec_params(5,i)**2
       qmi_QQ2 = qmi_QQ*qmi_QQ+qmi_oneBDD3

       do j=1,npairs

         loop_count = loop_count+1
         vec(1:3) = coords(1:3,i)-crdsmm(1:3,j)
         r2  = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
         rr2 = r2*A2_TO_BOHRS2
         onerij=one/sqrt(r2)
         qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count) = onerij
         rij=r2*onerij !one/onerij
         qm2_rij_eqns%qmmmrijdata(QMMMRIJ,    loop_count) = rij
         qm2_rij_eqns%qmmmrijdata(QMMMEXP1,   loop_count) = exp(-qmi_alpa*rij)
         qm2_rij_eqns%qmmmrijdata(QMMMEXP2,   loop_count) = exp(-ALPH_MM*rij)
         qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count) = one/sqrt(RR2+qmi_oneBDD1)

         if (sp_atom) then
           ! SP-atom specific stuff
           rr=rij*A_TO_BOHRS
           RRADD = RR+qmi_DD
           RRADD = RRADD*RRADD+qmi_oneBDD2
           qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRADDADE,loop_count) = one/sqrt(RRADD)
           RRMDD = RR-qmi_DD
           RRMDD = RRMDD*RRMDD+qmi_oneBDD2
           qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMDDADE,loop_count) = one/sqrt(RRMDD)
           qm2_rij_eqns%qmmmrijdata(QMMMSQRTRR2AQE,  loop_count) = one/SQRT(RR2+qmi_oneBDD3)
           RRAQQ=RR+qmi_QQ
           RRAQQ=RRAQQ*RRAQQ+qmi_oneBDD3
           qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQAQE,loop_count) = one/SQRT(RRAQQ)
           RRMQQ=RR-qmi_QQ
           RRMQQ=RRMQQ*RRMQQ+qmi_oneBDD3
           qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMQQAQE, loop_count) = one/SQRT(RRMQQ)
           qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQ2AQE,loop_count) = one/SQRT(RR2+qmi_QQ2)
         end if !(sp_atom)

       end do

     end do

  end if

end subroutine qm2_calc_rij_and_eqns

subroutine qm2_allocate_qm2_qmmm_rij_eqns(qmmm_nml,qm2_rij_eqns, qmmm_mpi, qmmm_struct,natom,npairs)

   use qmmm_nml_module   , only : qmmm_nml_type
   use qmmm_module, only : qm2_rij_eqns_structure, qmmm_mpi_structure
   use qmmm_struct_module, only : qmmm_struct_type

   implicit none
!Passed in
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qm2_rij_eqns_structure),intent(inout) :: qm2_rij_eqns
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
  type(qmmm_struct_type), intent(inout) :: qmmm_struct
  integer, intent(in) :: natom, npairs

!Local
  integer :: ier=0
  integer array_size, num_per_thread

  !We will allocate the qmmmrijdata array enough to hold the npairs+a bit
  !or qmmm_struct%nquant * (natom - qmmm_struct%nquant)
  !whichever is smaller. We then store the amount we allocated
  !in qm2_rij_eqns%qmmmrij_allocated. In this way when the number
  !of pairs changes this value can be checked and the array reallocated
  !large if necessary.

  !In parallel it only needs to be (qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1)*qmmm_struct%qm_mm_pairs
  !in size

  if (qmmm_nml%qmmmrij_incore) then
     num_per_thread = qmmm_mpi%nquant_nlink_end-qmmm_mpi%nquant_nlink_start+1
     array_size = npairs*(num_per_thread)+npairs
     array_size = min(array_size,(num_per_thread) &
                               * (natom - qmmm_struct%nquant))
     allocate(qm2_rij_eqns%qmmmrijdata(QMMMNORIJ,array_size),stat=ier)
     REQUIRE(ier==0)
     qm2_rij_eqns%qmmmrij_allocated=array_size
     if (qmmm_nml%verbosity>1 .and. qmmm_mpi%commqmmm_master) then
       write(6,'(a,i6,a,i6)') 'QMMM: Allocating qmmmrijdata array as ', &
           QMMMNORIJ,' x ',array_size
       write(6,'(a,i10,a,i8)') 'QMMM: to hold min npairs of ', &
          npairs,' and nquant per thread of ',num_per_thread
     end if
  end if
  return
end subroutine qm2_allocate_qm2_qmmm_rij_eqns

