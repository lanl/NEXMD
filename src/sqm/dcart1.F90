! 
#include "copyright.h"
#include "dprec.fh"

! CML Keep an eye on qm2_get_qm_forces() in future versions to make sure we
! CML implement any changes made to them so we can make the same changes in
! CML this subroutine, as well as DCART1() and DCART2(), since they are based
! CML on the same subroutine. 7/13/12

subroutine dcart1(dxyzqm, gs_dm, ex_dm, xyz_in) ! CML add coordinates passed in 7/13/12
!Current code maintained by: Ross Walker (TSRI 2004)

!This routine calculates the derivatives of the energy for QM-QM
!interactions.
!
! xyz_in - QM Atom Coordinates
! dxyzqm  - Returned with the forces in for each QM atom.

      use constants          , only : EV_TO_KCAL
      use ElementOrbitalIndex, only: MaxValenceOrbitals
      use qmmm_module        , only : qmmm_nml,qmmm_struct, qm2_struct, qm2_params, qmmm_mpi
      use qm2_pm6_hof_module
      use dh_correction_module, only : dh_correction_grad
      use dcart_xpm_module !BTN 8/7/17 use dhc1 and h1elec functions from here
	  use qm2_davidson_module ! CML 7/13/12

 
       implicit none     
      _REAL_, parameter :: change=2.0D-6, halfChange=change/2.0D0, oneChange=1.0D0/change
      _REAL_, parameter :: delAdj =1.0D-8, twoOnedelAdj= 0.5D0/delAdj    

!Passed in
      _REAL_, intent(inout) :: dxyzqm(3,qmmm_struct%nquant_nlink)                                  
	  _REAL_, intent(in) :: xyz_in(3,qmmm_struct%nquant_nlink)	! CML Just in case we don't update coords in qmmm_struct 7/13/12
	  _REAL_, intent(in) :: gs_dm(qm2ds%Nb*(qm2ds%Nb+1)/2) ! CML 7/13/12
	  _REAL_, intent(in) :: ex_dm(qm2ds%Nb*(qm2ds%Nb+1)/2) ! CML 7/13/12
!Local


      _REAL_ e_repul(22) !Used when qmqm_erep_incore = false
      _REAL_ pair_force(3)
      integer loop_count !Keeps track of number of times through nquant * (nquant-1)/2 loop
!      _REAL_ psum(36) !36 = max combinations with heavy and heavy = 4 orbs * 4 orbs (Note, no d orb support)
      _REAL_ psum(MaxValenceOrbitals**2*3) 
      _REAL_ ptzsum(MaxValenceOrbitals**2*3) ! CML for excited state DM 7/13/12
      _REAL_ xyz_qmi(3), xyz_qmj(3), vec_qm_qm1, vec_qm_qm2, vec_qm_qm3
      _REAL_ DENER ! BTN 08/10/2017 done to store energy derivative
      integer natqmi, natqmj, qmitype, qmjtype
      integer ii, iif, iil, jj, jjf, jjl, ij
      integer i,j,k,l
      integer n_atomic_orbi, n_atomic_orbj
      integer jstart, jend
      _REAL_ aa,ee,deriv,angle,refh,heat,sum
      _REAL_ corei, corej
      _REAL_ betasas, betasap, betapas, betapap
      _REAL_ bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, bdd3j
      _REAL_ qqi, qqi2, qqj, qqj2, ddi,ddj
      _REAL_ htype, fqmii(3)
      _REAL_, target :: F(MaxValenceOrbitals*(MaxValenceOrbitals*2+1)) !BTN 10/08/2017 place to store fock matrix
      integer natom
      
!#define change 1.D-4
!#define halfChange 5.D-5
!!one/change = 10000
!#define onechange 10000
!#define delAdj 1.0D-8
!#define TWOONEdelAdj 50000000

   if (qmmm_nml%qmqm_exc_analyt) then !We do analytical derivatives
		! CML as of right now, fully analytical derivatives for the excited
		! CML state are not possible. So, we must exit.
		write(6,*)  'QMMM: Analytical excited state derivatives are not implemented! Exiting.'
		call mexit(6,1)
   else !We will do (pseudo numerical derivatives)
!************** PSEUDO NUMERICAL DERIVATIVES **************
#ifdef MPI
      do ii = qmmm_mpi%nquant_nlink_istart, qmmm_mpi%nquant_nlink_iend
         jstart =  qmmm_mpi%nquant_nlink_jrange(1,ii)
         jend = qmmm_mpi%nquant_nlink_jrange(2,ii)
#else
         do II=2,qmmm_struct%nquant_nlink
           jstart = 1
           jend = ii-1
#endif
          !Loop over all pairs of quantum atoms
          iif=qm2_params%orb_loc(1,II)
          iil=qm2_params%orb_loc(2,II)
          qmitype = qmmm_struct%qm_atom_type(ii)
          natqmi=qmmm_struct%iqm_atomic_numbers(II)
          xyz_qmi(1)=xyz_in(1,II)
          xyz_qmi(2)=xyz_in(2,II)
          xyz_qmi(3)=xyz_in(3,II)
          do JJ=jstart,jend !jj=1,ii-1
!  FORM DIATOMIC MATRICES
            jjf=qm2_params%orb_loc(1,JJ)
            jjl=qm2_params%orb_loc(2,JJ)
!   GET FIRST ATOM
            qmjtype = qmmm_struct%qm_atom_type(jj)
            natqmj=qmmm_struct%iqm_atomic_numbers(JJ)
            xyz_qmj(1)=xyz_in(1,JJ)
            xyz_qmj(2)=xyz_in(2,JJ)
            xyz_qmj(3)=xyz_in(3,JJ)
            IJ=0
            do I=jjf,jjl
              K=qm2_params%pascal_tri1(i)+jjf-1
              do J=jjf,I
                IJ=IJ+1
                K=K+1
                psum(IJ)=gs_dm(K)
                ptzsum(IJ)=ex_dm(K)
              end do
            end do
! GET SECOND ATOM FIRST ATOM INTERSECTION
            do I=iif,iil
               L=qm2_params%pascal_tri1(i)
               K=L+jjf-1 
               do J=jjf,jjl
                  IJ=IJ+1
                  K=K+1
                  psum(IJ)=gs_dm(K)
                  ptzsum(IJ)=ex_dm(K)
               end do
               K=L+iif-1
               do L=iif,I
                   K=K+1
                   IJ=IJ+1
                   psum(IJ)=gs_dm(K)
                   ptzsum(IJ)=ex_dm(K)
               end do
            end do
            do K=1,3
              xyz_qmi(K)=xyz_qmi(K)+halfChange
              call qm2_dhc1(psum,ii,jj,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi,natqmj,iif,iil,jjf, &
                       jjl,F)
              
              DENER=qm2_helect1(iil-iif+jjl-jjf+1,ptzsum,F)   ! CML 7/13/12 BTN 08/10/2017
              AA=DENER*2.d0
                       
              xyz_qmi(K)=xyz_qmi(K)-change
              call qm2_dhc1(psum,ii,jj,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi,natqmj,iif,iil,jjf, &
                       jjl,F)
                       
              DENER=qm2_helect1(iil-iif+jjl-jjf+1,ptzsum,F)   ! CML 7/13/12 BTN 08/10/2017
              EE=DENER*2.d0
                       
                       
              xyz_qmi(K)=xyz_qmi(K)+halfChange
                   
              DERIV=(EE-AA)*EV_TO_KCAL*onechange
              dxyzqm(K,II)=dxyzqm(K,II)-DERIV
              dxyzqm(K,JJ)=dxyzqm(K,JJ)+DERIV

            end do
          end do
       end do
!************** end PSEUDO NUMERICAL DERIVATIVES **************
   end if                                           

end subroutine dcart1
