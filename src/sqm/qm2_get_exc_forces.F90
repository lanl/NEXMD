! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"

! CML Keep an eye on qm2_get_qm_forces() in future versions to make sure we
! CML implement any changes made to them so we can make the same changes in
! CML this subroutine, as well as DCART1() and DCART2(), since they are based
! CML on the same subroutine. 7/13/12

subroutine qm2_get_exc_forces(qm2_params,qmmm_nml, qm2_rij_eqns, qmmm_mpi,qm2_struct,qmmm_struct,dxyzqm, xyz_in) ! CML add coordinates passed in 7/13/12
!Current code maintained by: Ross Walker (TSRI 2004)

!This routine calculates the derivatives of the energy for QM-QM
!interactions.

! xyz_in - QM Atom Coordinates
! dxyzqm  - Returned with the forces in for each QM atom.

      use constants          , only : EV_TO_KCAL
      use ElementOrbitalIndex, only: MaxValenceOrbitals
      use qmmm_module        , only : qm2_structure, qmmm_mpi_structure, qm2_rij_eqns_structure
      use qm2_pm6_hof_module
      use dh_correction_module, only : dh_correction_grad
      use qmmm_struct_module, only : qmmm_struct_type
      use qm2_params_module,  only : qm2_params_type
      use qmmm_nml_module   , only : qmmm_nml_type

 
       implicit none     
      _REAL_, parameter :: change=2.0D-6, halfChange=change/2.0D0, oneChange=1.0D0/change
      _REAL_, parameter :: delAdj =1.0D-8, twoOnedelAdj= 0.5D0/delAdj    

!Passed in
   type(qm2_rij_eqns_structure),intent(inout) :: qm2_rij_eqns
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
      type(qm2_structure),intent(inout) :: qm2_struct
      type(qmmm_struct_type), intent(inout) :: qmmm_struct
      _REAL_, intent(out) :: dxyzqm(3,qmmm_struct%nquant_nlink)                                  
      _REAL_, intent(inout) :: xyz_in(3,qmmm_struct%nquant_nlink)	! CML Just in case we don't update coords in qmmm_struct 7/13/12

!Local


      _REAL_ e_repul(22) !Used when qmqm_erep_incore = false
      _REAL_ pair_force(3)
      integer loop_count !Keeps track of number of times through nquant * (nquant-1)/2 loop
      _REAL_ psum(MaxValenceOrbitals**2*3) 
      _REAL_ xyz_qmi(3), xyz_qmj(3), vec_qm_qm1, vec_qm_qm2, vec_qm_qm3
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
      integer natom
      

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
                psum(IJ)=qm2_struct%den_matrix(K)
              end do
            end do
! GET SECOND ATOM FIRST ATOM INTERSECTION
            do I=iif,iil
               L=qm2_params%pascal_tri1(i)
               K=L+jjf-1
               do J=jjf,jjl
                  IJ=IJ+1
                  K=K+1
                  psum(IJ)=qm2_struct%den_matrix(K)
               end do
               K=L+iif-1
               do L=iif,I
                   K=K+1
                   IJ=IJ+1
                   psum(IJ)=qm2_struct%den_matrix(K)
               end do
            end do
            do K=1,3
              xyz_qmi(K)=xyz_qmi(K)+halfChange
              call qm2_exc_dhc(qm2_params, qmmm_nml, qm2_rij_eqns, qm2_struct, qmmm_struct, &
                       psum,ii,jj,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi,natqmj,iif,iil,jjf, &
                       jjl,AA)
              xyz_qmi(K)=xyz_qmi(K)-change
              call qm2_exc_dhc(qm2_params, qmmm_nml, qm2_rij_eqns, qm2_struct, qmmm_struct, &
                       psum,ii,jj,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi,natqmj,iif,iil,jjf, &
                       jjl,EE)
              xyz_qmi(K)=xyz_qmi(K)+halfChange
                   
              DERIV=(EE-AA)*EV_TO_KCAL*onechange
              dxyzqm(K,II)=dxyzqm(K,II)-DERIV
              dxyzqm(K,JJ)=dxyzqm(K,JJ)+DERIV

            end do
          end do
       end do
!************** end PSEUDO NUMERICAL DERIVATIVES **************
   end if
   ! --------------------------------------------
   ! PM6: Gradient of PM6 correction to HOF
   !      (nitrogen non-planarity energy penalty)
   ! --------------------------------------------
   if (qmmm_mpi%commqmmm_master) then
      ! this is not parallelized - do only on the master
      if (qmmm_nml%qmtheory%PM6) then
         natom = qmmm_struct%nquant_nlink
         call hofCorrectionGradient(qmmm_struct, natom, dxyzqm)
      end if
      if (qmmm_nml%qmtheory%DISPERSION .or. qmmm_nml%qmtheory%DISPERSION_HYDROGENPLUS) then
         call dh_correction_grad(qm2_params,qm2_struct, qmmm_struct%nquant_nlink,xyz_in, &
                                 qmmm_struct%iqm_atomic_numbers,qmmm_nml%qmtheory,dxyzqm)
      endif
   end if

   if(qmmm_nml%peptide_corr) then
!  NOW ADD IN MOLECULAR-MECHANICS CORRECTION TO THE H-N-C=O TORSION            
     if (qmmm_nml%qmtheory%PM3 .OR. qmmm_nml%qmtheory%PDDGPM3 .OR. qmmm_nml%qmtheory%PM3CARB1 &
         .OR. qmmm_nml%qmtheory%PM3ZNB .OR. qmmm_nml%qmtheory%PDDGPM3_08) then
       htype = 7.1853D0                                                      
     elseif (qmmm_nml%qmtheory%AM1 .OR. qmmm_nml%qmtheory%RM1) then
       htype = 3.3191D0                                                      
     else !Assume MNDO
       htype = 6.1737D0                                                      
     end if
!Parallel
     do I=qmmm_mpi%mytaskid+1,qm2_struct%n_peptide_links,qmmm_mpi%numthreads !1,n_peptide_links
       do J=1,4                                    
         do K=1,3                                
           xyz_in(K,qm2_struct%peptide_links(J,I))= &
                          xyz_in(K,qm2_struct%peptide_links(J,I))-delAdj
           call qm2_dihed(xyz_in,qm2_struct%peptide_links(1,I), &
                          qm2_struct%peptide_links(2,I),qm2_struct%peptide_links(3,I), &
                          qm2_struct%peptide_links(4,I),ANGLE)
           REFH=HTYPE*SIN(ANGLE)**2         
           xyz_in(K,qm2_struct%peptide_links(J,I))= &
                          xyz_in(K,qm2_struct%peptide_links(J,I))+delAdj*2.D0
           call qm2_dihed(xyz_in,qm2_struct%peptide_links(1,I),qm2_struct%peptide_links(2,I), &
                          qm2_struct%peptide_links(3,I),qm2_struct%peptide_links(4,I),ANGLE)
           xyz_in(K,qm2_struct%peptide_links(J,I))= &
                          xyz_in(K,qm2_struct%peptide_links(J,I))-delAdj
           HEAT=HTYPE*SIN(ANGLE)**2         
           SUM=(REFH-HEAT)*TWOONEdelAdj
           dxyzqm(K,qm2_struct%peptide_links(J,I))=dxyzqm(K,qm2_struct%peptide_links(J,I))-SUM 
         end do
       end do                                   
     end do                                    
   end if                                           

end subroutine qm2_get_exc_forces

subroutine qm2_exc_dhc(qm2_params, qmmm_nml, qm2_rij_eqns, qm2_struct, qmmm_struct, &
                   P,iqm, jqm,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi, &
                   natqmj, iif, iil, jjf, jjl, DENER)
!***********************************************************************
!
! Ross Walker (SDSC, 2006) : Do 'Pseudo' Numerical Derivatives for QM
!
! d-orbital extension,  (Taisung Lee, 2011)  
!
!***********************************************************************

      use constants          , only: ONE, A_TO_BOHRS, A2_TO_BOHRS2, EV_TO_KCAL
      use ElementOrbitalIndex, only: MaxValenceOrbitals,MaxValenceDimension 
      use qmmm_module        , only: OVERLAP_CUTOFF, qm2_structure, qm2_rij_eqns_structure
      use Rotation           , only: GetRotationMatrix, Rotate2Center2Electron, RotateCore   
      use qm2_fock_d         , only: W2Fock_atompair
      use qmmm_struct_module, only : qmmm_struct_type
      use qm2_params_module,  only : qm2_params_type
      use qmmm_nml_module   , only : qmmm_nml_type
 
      implicit none

!Passed in
      type(qm2_rij_eqns_structure),intent(inout) :: qm2_rij_eqns
      type(qm2_params_type),intent(inout) :: qm2_params
      type(qmmm_nml_type),intent(inout) :: qmmm_nml
      type(qm2_structure),intent(inout) :: qm2_struct
      type(qmmm_struct_type), intent(inout) :: qmmm_struct
      _REAL_ P(*)
      _REAL_, intent(in)  :: xyz_qmi(3),xyz_qmj(3)
      integer, intent(in) :: iqm, jqm, natqmi, natqmj, qmitype, qmjtype
      integer, intent(in) :: iif, iil, jjf, jjl
      _REAL_, intent(out) :: DENER

! Local
      integer :: n_atomic_orbi, n_atomic_orbj
      integer :: i_dimension, j_dimension
      integer :: i,j,k,j1,jj,i1, linear, i2, ii, j2
      integer :: firstIndexAO_i, lastIndexAO_i, firstIndexAO_j, lastIndexAO_j      

      _REAL_ :: H(MaxValenceOrbitals*(MaxValenceOrbitals*2+1))
      _REAL_ :: F(MaxValenceOrbitals*(MaxValenceOrbitals*2+1))
      _REAL_ :: SHMAT(MaxValenceOrbitals,MaxValenceOrbitals)
      _REAL_ :: W(MaxValenceDimension**2)
      _REAL_ :: enuclr, ee, temp
      _REAL_ :: r2, rij, r2InAu, rijInAu, oneOverRij
      
      _REAL_ :: RI(22), CORE(10,2)
      _REAL_ :: rotationMatrix(15,45)
      _REAL_, allocatable:: WW(:,:)
      integer ::orb_loc(2,2),KR
      logical::hasDOrbital

      !qm2_Helect is a function
      _REAL_ qm2_helect
      _REAL_ vasya;

      if (iif < jjf) then
        i=iif-1
        j=jjf-iil+iif-2
      else
        i=iif-jjl+jjf-2
        j=jjf-1
      end if
      
      firstIndexAO_i=iif-i
      lastIndexAO_i=iil-i
      firstIndexAO_j=jjf-j
      lastIndexAO_j=jjl-j     
      
      n_atomic_orbi = lastIndexAO_i-firstIndexAO_i+1
      n_atomic_orbj = lastIndexAO_j-firstIndexAO_j+1
      
      linear=qm2_params%pascal_tri2(n_atomic_orbi+n_atomic_orbj)
      F(1:linear)=0.0D0
      H(1:linear)=0.0D0

! RCW: Caution, i and j reversed here.

        r2 = sum((xyz_qmj-xyz_qmi)**2) 
        rij=sqrt(r2)
        rijInAu=rij*A_TO_BOHRS
        oneOverRij=one/rij
        r2InAu=r2*A2_TO_BOHRS2
        
      if (r2InAu < OVERLAP_CUTOFF) then
      
        if ((n_atomic_orbi.lt.9) .and. (n_atomic_orbj.lt.9)) then  ! SP case
            SHMAT=0.0d0
             call qm2_h1elec(qm2_params,r2InAu,xyz_qmi(1),                                  &
                   xyz_qmj(1),n_atomic_orbi, n_atomic_orbj, SHMAT,           &
                   qmitype,qmjtype)
            
              I2=0
              do I1=firstIndexAO_i,lastIndexAO_i
                 II=qm2_params%pascal_tri1(i1)+firstIndexAO_j-1
                 I2=I2+1
                 J2=0
                 JJ=MIN(I1,lastIndexAO_j)
                 do J1=firstIndexAO_j,JJ                                                   
                    II=II+1                                                       
                    J2=J2+1                                                       
                    H(II)=H(II)+SHMAT(I2,J2) 
                 end do
              end do
            
        else  ! for atoms with d orbitals
        
            call qm2_h1elec_d(qm2_params,qm2_struct,r2InAu,xyz_qmi(1:3), xyz_qmj(1:3),  &
                        n_atomic_orbi,n_atomic_orbj,                &
                        firstIndexAO_i, firstIndexAO_j, qmitype, qmjtype,  &
                        n_atomic_orbi+n_atomic_orbj, H)                

        end if ! ((n_atomic_orbi.lt.9) .and. (n_atomic_orbj.lt.9))
      end if !(R2 < OVERLAP_CUTOFF)

      KR=1
      hasDOrbital=((n_atomic_orbi.ge.9) .or. (n_atomic_orbj.ge.9))
      call GetRotationMatrix(qm2_params,xyz_qmj-xyz_qmi, rotationMatrix, hasDOrbital)        
      call qm2_rotate_qmqm(qmmm_nml, qm2_params, qm2_rij_eqns,qm2_struct, qmmm_struct,&
                  -1,iqm,jqm,natqmi,natqmj,xyz_qmi,xyz_qmj,            &
                  W(KR),KR, RI, core)

      if (hasDOrbital) then   ! spd case

        i_dimension=n_atomic_orbi*(n_atomic_orbi+1)/2
        j_dimension=n_atomic_orbj*(n_atomic_orbj+1)/2
       
        allocate(ww(1:j_dimension, 1:i_dimension)) 
        WW=0.0D0

        ! calculate the 2-center integrals and core-core interaction integrals
        call qm2_repp_d(qm2_params,qmmm_struct, qmitype,qmjtype,rijInAu,RI,CORE,WW,i_dimension,j_dimension,1)
 
        ! put 2-center 2-electron integrals to the linearized matrix W

        k=0
        do ii=1,i_dimension
          do jj=1, j_dimension
            k=k+1
            W(k)=WW(jj,ii)
          end do
        end do
                   
        deallocate(ww)
        call Rotate2Center2Electron(W, i_dimension, j_dimension,rotationMatrix)
        
   end if ! ((n_atomic_orbi.ge.9) .and. (n_atomic_orbj.ge.9))
   
   
   ! calculate the core-core contribution to the H matrix
    ii=qm2_params%pascal_tri2(firstIndexAO_i)
    jj=qm2_params%pascal_tri2(firstIndexAO_j)   
      
   call RotateCore(firstIndexAO_i,firstIndexAO_j,              &
    n_atomic_orbi,n_atomic_orbj,  &
    ii,jj,core,rotationMatrix,H)
   
   call qm2_core_core_repulsion(qmmm_nml, qm2_params,qmmm_struct, iqm, jqm, rij, oneOverRij, RI, enuclr)         
        
    ! put what we have now to the Fock matrix
    F(1:linear)=H(1:linear)
       
    ! 2-center 2-electron contribution to the Fock matrix      
    call W2Fock_atompair(W, F, P, n_atomic_orbj, n_atomic_orbi,  &
      firstIndexAO_j, firstIndexAO_i,qmmm_struct%W2Fock_atompair_initialized, qmmm_struct%w_index)
    call W2Fock_atompair(W, F, P, n_atomic_orbi, n_atomic_orbj,  &
      firstIndexAO_i, firstIndexAO_j,qmmm_struct%W2Fock_atompair_initialized, qmmm_struct%w_index)   

    EE=qm2_HELECT(n_atomic_orbi+n_atomic_orbj-1,P,H,F,vasya)   
    DENER=EE+ENUCLR
   
   
end subroutine qm2_exc_dhc




