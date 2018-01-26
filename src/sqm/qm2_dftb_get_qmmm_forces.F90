! <compile=optimized>

#include "dprec.fh"

subroutine qm2_dftb_get_qmmm_forces(qmmm_mpi,qm2_struct, qmmm_struct, dxyzcl,dxyzqm, vectmp1,vectmp2,vectmp3,vectmp4)
!Calculate interaction between REAL QM and MM atoms in the list. Exclude MM atoms.
! Due to the different charge magnitude in quantum and classical parts,
! a scaling may be necessary in this interaction. See eq. (4) in
! Elstner et al, J. Mol. Struc. (Theochem), 632, 24--41 (2003).

! WHEN CHANGING HERE, REMEMBER TO MAKE THE SAME CHANGES IN 
! qm2_dftb_externalshift.

! Hubbard Parameter (for the scaling factor):
! uhub=mcharge%uhubb(izp(j))

! \zeta = 4.0 in Eq. (4):
! gamma = 1.0/sqrt(r2 + (0.5/uhub + 0.5/uhub)**2)

! \zeta = 0.4:
! gamma = 1.0/sqrt(r2 + 0.1 * (1.0/uhub)**2)
! In DFTB, \zeta in Eq. (4) is set to zero. 
! This is similar to no scaling at all. (p.33)

! Vector version written by Ross Walker (SDSC 2006)

  use qmmm_module, only : qm2_structure, qmmm_mpi_structure
  use constants, only: AU_TO_KCAL, BOHRS_TO_A
  use qmmm_struct_module, only : qmmm_struct_type

  implicit none

  !Passed in
   type(qm2_structure),intent(inout) :: qm2_struct
   type(qmmm_struct_type), intent(in) :: qmmm_struct
   type(qmmm_mpi_structure), intent(in) :: qmmm_mpi
  _REAL_ , intent(out) :: dxyzcl(3,qmmm_struct%qm_mm_pairs)
  _REAL_ , intent(inout) :: dxyzqm(3,qmmm_struct%nquant_nlink)
  _REAL_ , intent(out) :: vectmp1(*), vectmp2(*), &
                          vectmp3(*), vectmp4(*)
                         !These should have been allocated to at least qm_mm_npair long.
  
  !Local
   _REAL_  :: qmx, qmy, qmz, scf_mchgi
  integer :: i

!  do i = 1,qmmm_struct%nquant_nlink
  do i = qmmm_mpi%nquant_nlink_start,qmmm_mpi%nquant_nlink_end
     scf_mchgi = qm2_struct%scf_mchg(i)*AU_TO_KCAL*BOHRS_TO_A !Unit conversion
     qmx = qmmm_struct%qm_coords(1,i)
     qmy = qmmm_struct%qm_coords(2,i)
     qmz = qmmm_struct%qm_coords(3,i)
!Vector Code
     vectmp1(1:qmmm_struct%qm_mm_pairs) = qmx - qmmm_struct%qm_xcrd(1,1:qmmm_struct%qm_mm_pairs)
     vectmp2(1:qmmm_struct%qm_mm_pairs) = qmy - qmmm_struct%qm_xcrd(2,1:qmmm_struct%qm_mm_pairs)
     vectmp3(1:qmmm_struct%qm_mm_pairs) = qmz - qmmm_struct%qm_xcrd(3,1:qmmm_struct%qm_mm_pairs)
     vectmp4(1:qmmm_struct%qm_mm_pairs) = vectmp1(1:qmmm_struct%qm_mm_pairs)*vectmp1(1:qmmm_struct%qm_mm_pairs) &
                                        + vectmp2(1:qmmm_struct%qm_mm_pairs)*vectmp2(1:qmmm_struct%qm_mm_pairs) &
                                        + vectmp3(1:qmmm_struct%qm_mm_pairs)*vectmp3(1:qmmm_struct%qm_mm_pairs)
     !vectmp4 contains r2 in angstroms^2
     call vdinvsqrt(qmmm_struct%qm_mm_pairs,vectmp4,vectmp4)
     !vectmp4 now contains 1/r
     vectmp4(1:qmmm_struct%qm_mm_pairs) = vectmp4(1:qmmm_struct%qm_mm_pairs) &
                                         *vectmp4(1:qmmm_struct%qm_mm_pairs) &
                                         *vectmp4(1:qmmm_struct%qm_mm_pairs)
     !vectmp4 now contains 1/r^3 = gamma^3
     !qm_xcrd(4,...) = mm charge in electrons
     vectmp4(1:qmmm_struct%qm_mm_pairs) = qmmm_struct%qm_xcrd(4,1:qmmm_struct%qm_mm_pairs) &
                                          *vectmp4(1:qmmm_struct%qm_mm_pairs) &
                                          *scf_mchgi
     vectmp1(1:qmmm_struct%qm_mm_pairs) = vectmp1(1:qmmm_struct%qm_mm_pairs)*vectmp4(1:qmmm_struct%qm_mm_pairs)
     vectmp2(1:qmmm_struct%qm_mm_pairs) = vectmp2(1:qmmm_struct%qm_mm_pairs)*vectmp4(1:qmmm_struct%qm_mm_pairs)
     vectmp3(1:qmmm_struct%qm_mm_pairs) = vectmp3(1:qmmm_struct%qm_mm_pairs)*vectmp4(1:qmmm_struct%qm_mm_pairs)
     !vectmp1,2,3 now contain the gradients in the x,y, and z direction in kcal/mol/a^2
     dxyzcl(1,1:qmmm_struct%qm_mm_pairs) = dxyzcl(1,1:qmmm_struct%qm_mm_pairs)+vectmp1(1:qmmm_struct%qm_mm_pairs)
     dxyzcl(2,1:qmmm_struct%qm_mm_pairs) = dxyzcl(2,1:qmmm_struct%qm_mm_pairs)+vectmp2(1:qmmm_struct%qm_mm_pairs)
     dxyzcl(3,1:qmmm_struct%qm_mm_pairs) = dxyzcl(3,1:qmmm_struct%qm_mm_pairs)+vectmp3(1:qmmm_struct%qm_mm_pairs)
     dxyzqm(1,i) = dxyzqm(1,i) - sum(vectmp1(1:qmmm_struct%qm_mm_pairs))
     dxyzqm(2,i) = dxyzqm(2,i) - sum(vectmp2(1:qmmm_struct%qm_mm_pairs))
     dxyzqm(3,i) = dxyzqm(3,i) - sum(vectmp3(1:qmmm_struct%qm_mm_pairs))
  end do 
    
!Original Scalar code 
!     do j = 1,qmmm_struct%qm_mm_pairs
!        vec(1:3) = qmx(1:3) - qmmm_struct%qm_xcrd(1:3,j)
!        r2 = vec(1)*vec(1)+vec(2)*vec(2)+vec(3)*vec(3)
!        gammar = 1.0d0/sqrt(r2) !1/r
!        gammar = gammar*gammar*gammar !Gamma = Gamma^3
!        qdiff = qmmm_struct%qm_xcrd(4,j)*scf_mchgi*gammar
!        dgr(1:3) = vec(1:3)*qdiff*AU_TO_KCAL*BOHRS_TO_A
!        dxyzmm(1:3,j) = dxyzmm(1:3,j) + dgr(1:3)
!        qmsum(1:3) = qmsum(1:3) + dgr(1:3)
!     enddo
!     qmmm_struct%dxyzqm(1:3,i) = qmmm_struct%dxyzqm(1:3,i) - qmsum(1:3)
!  enddo

  return 
  
end subroutine qm2_dftb_get_qmmm_forces

