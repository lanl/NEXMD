! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
subroutine qm2_get_qm_forces(qm2_struct, qmmm_struct, dxyzqm)

!Current code maintained by: Ross Walker (TSRI 2004)

!This routine calculates the derivatives of the energy for QM-QM
!interactions.

! qmmm_struct%qm_coords - QM Atom Coordinates
! dxyzqm  - Returned with the forces in for each QM atom.

      use constants          , only : EV_TO_KCAL
      use ElementOrbitalIndex, only: MaxValenceOrbitals
      use qmmm_module        , only : qmmm_nml, qm2_structure, qm2_params, qmmm_mpi
      use qm2_pm6_hof_module
      use dh_correction_module, only : dh_correction_grad
      use qmmm_struct_module, only : qmmm_struct_type

 
       implicit none     

      type(qmmm_struct_type), intent(inout) :: qmmm_struct
      type(qm2_structure), intent(inout) :: qm2_struct
 
      _REAL_, parameter :: change=2.0D-6, halfChange=change/2.0D0, oneChange=1.0D0/change
      _REAL_, parameter :: delAdj =1.0D-8, twoOnedelAdj= 0.5D0/delAdj    

!Passed in
      _REAL_, intent(out) :: dxyzqm(3,qmmm_struct%nquant_nlink)                                  

!Local


      _REAL_ e_repul(22) !Used when qmqm_erep_incore = false
      _REAL_ pair_force(3)
      integer loop_count !Keeps track of number of times through nquant * (nquant-1)/2 loop
!      _REAL_ psum(36) !36 = max combinations with heavy and heavy = 4 orbs * 4 orbs (Note, no d orb support)
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
      
!#define change 1.D-4
!#define halfChange 5.D-5
!!one/change = 10000
!#define onechange 10000
!#define delAdj 1.0D-8
!#define TWOONEdelAdj 50000000

   if (qmmm_nml%qmqm_analyt) then !We do analytical derivatives
   !RCW: Note: there is a lot of repeated code in the two options below
   !but this is necessary to factor this if statement out of the inner loop.
      loop_count = 0
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
!   GET FIRST ATOM INFO                                                             
         iif=qm2_params%orb_loc(1,II)                                                          
         iil=qm2_params%orb_loc(2,II) 
!             n_atomic_orbi = iil - iif + 1
         n_atomic_orbi = qm2_params%natomic_orbs(ii)
         natqmi=qmmm_struct%iqm_atomic_numbers(II)               
         corei=qm2_params%core_chg(ii)

         ddi = qm2_params%multip_2c_elec_params(1,ii)
         qqi = qm2_params%multip_2c_elec_params(2,ii)
         qqi2 = qqi*qqi
         bdd1i = qm2_params%multip_2c_elec_params(3,ii)
         bdd2i = qm2_params%multip_2c_elec_params(4,ii)
         bdd3i = qm2_params%multip_2c_elec_params(5,ii)
         xyz_qmi(1:3)=qmmm_struct%qm_coords(1:3,II)        
         qmitype = qmmm_struct%qm_atom_type(ii)
         fqmii(1:3) = 0.0d0
         do JJ=jstart,jend  !jj=1,ii-1
!   GET SECOND ATOM INFO                                 
           jjf=qm2_params%orb_loc(1,JJ)                                                       
           jjl=qm2_params%orb_loc(2,JJ)                                                        
!           n_atomic_orbj = jjl - jjf + 1
           n_atomic_orbj = qm2_params%natomic_orbs(jj)
           natqmj=qmmm_struct%iqm_atomic_numbers(JJ)                                                      
           corej=qm2_params%core_chg(jj)
           ddj = qm2_params%multip_2c_elec_params(1,jj)
           qqj = qm2_params%multip_2c_elec_params(2,jj)
           qqj2 = qqj*qqj
           bdd1j = qm2_params%multip_2c_elec_params(3,jj)
           bdd2j = qm2_params%multip_2c_elec_params(4,jj)
           bdd3j = qm2_params%multip_2c_elec_params(5,jj)
           xyz_qmj(1:3)=qmmm_struct%qm_coords(1:3,JJ)    
           vec_qm_qm1=xyz_qmi(1) - xyz_qmj(1)
           vec_qm_qm2=xyz_qmi(2) - xyz_qmj(2)
           vec_qm_qm3=xyz_qmi(3) - xyz_qmj(3)
           qmjtype = qmmm_struct%qm_atom_type(jj)
           betasas = qm2_params%betasas(qmitype,qmjtype)
           betasap = qm2_params%betasap(qmitype,qmjtype)
           betapas = qm2_params%betasap(qmjtype,qmitype)
           betapap = qm2_params%betapap(qmitype,qmjtype)

           IJ=0                                                       
!  FORM DIATOMIC MATRICES                                                       
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

           loop_count=loop_count+1
           if (qmmm_nml%qmqm_erep_incore) then
             call qm2_deriv_qm_analyt(qmmm_struct, ii,jj,loop_count,qm2_struct%qm_qm_e_repul(1:22,loop_count), &
                       psum,n_atomic_orbj,n_atomic_orbi, &
                       corei,corej,betasas,betasap,betapas,betapap,vec_qm_qm1,vec_qm_qm2,  &
                       vec_qm_qm3, pair_force, qqi, qqi2, qqj, qqj2, ddi, ddj, &
                       bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, bdd3j, &
                       qm2_params%atom_orb_zz_sxs_over_sas(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_ss_eqn_adb(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xx1(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xx2(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xx1(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xx2(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xy(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xy(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_zz_pxp_over_pap(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_pp_eqn_xxy1(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_pp_eqn_xxy2(1,1,qmitype,qmjtype))
           else
!Same call as above, just qm2_struct%qm_qm_e_repul(1,loop_count) replaced with local e_repul
             call qm2_deriv_qm_analyt(qmmm_struct,ii,jj,loop_count,e_repul, &
                       psum,n_atomic_orbj,n_atomic_orbi, &
                       corei,corej,betasas,betasap,betapas,betapap,vec_qm_qm1,vec_qm_qm2,  &
                       vec_qm_qm3, pair_force, qqi, qqi2, qqj, qqj2, ddi, ddj, &
                       bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, bdd3j, &
                       qm2_params%atom_orb_zz_sxs_over_sas(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_ss_eqn_adb(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_zz_sxp_over_sap(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xx1(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xx2(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xx1(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xx2(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_sp_eqn_xy(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_sp_eqn_xy(1,1,qmjtype,qmitype), &
                       qm2_params%atom_orb_zz_pxp_over_pap(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_pp_eqn_xxy1(1,1,qmitype,qmjtype), &
                       qm2_params%atom_orb_pp_eqn_xxy2(1,1,qmitype,qmjtype))
           end if
           fqmii(1:3) = fqmii(1:3) + pair_force(1:3)
           dxyzqm(1:3,JJ)=dxyzqm(1:3,JJ)+pair_force(1:3)                          
         end do
         dxyzqm(1:3,II)=dxyzqm(1:3,II)-fqmii(1:3)
      end do

!************** end ANALYTICAL DERIVATIVES **************
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
          xyz_qmi(1)=qmmm_struct%qm_coords(1,II)
          xyz_qmi(2)=qmmm_struct%qm_coords(2,II)
          xyz_qmi(3)=qmmm_struct%qm_coords(3,II)
          do JJ=jstart,jend !jj=1,ii-1
!  FORM DIATOMIC MATRICES
            jjf=qm2_params%orb_loc(1,JJ)
            jjl=qm2_params%orb_loc(2,JJ)
!   GET FIRST ATOM
            qmjtype = qmmm_struct%qm_atom_type(jj)
            natqmj=qmmm_struct%iqm_atomic_numbers(JJ)
            xyz_qmj(1)=qmmm_struct%qm_coords(1,JJ)
            xyz_qmj(2)=qmmm_struct%qm_coords(2,JJ)
            xyz_qmj(3)=qmmm_struct%qm_coords(3,JJ)
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
              call qm2_dhc(qm2_struct, qmmm_struct, psum,ii,jj,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi,natqmj,iif,iil,jjf, &
                       jjl,AA)
              xyz_qmi(K)=xyz_qmi(K)-change
              call qm2_dhc(qm2_struct, qmmm_struct, psum,ii,jj,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi,natqmj,iif,iil,jjf, &
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
         call hofCorrectionGradient(qmmm_struct,natom, dxyzqm)
      end if
      if (qmmm_nml%qmtheory%DISPERSION .or. qmmm_nml%qmtheory%DISPERSION_HYDROGENPLUS) then
         call dh_correction_grad(qm2_struct, qmmm_struct%nquant_nlink,qmmm_struct%qm_coords, &
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
           qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))= &
                          qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))-delAdj
           call qm2_dihed(qmmm_struct%qm_coords,qm2_struct%peptide_links(1,I), &
                          qm2_struct%peptide_links(2,I),qm2_struct%peptide_links(3,I), &
                          qm2_struct%peptide_links(4,I),ANGLE)
           REFH=HTYPE*SIN(ANGLE)**2         
           qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))= &
                          qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))+delAdj*2.D0
           call qm2_dihed(qmmm_struct%qm_coords,qm2_struct%peptide_links(1,I),qm2_struct%peptide_links(2,I), &
                          qm2_struct%peptide_links(3,I),qm2_struct%peptide_links(4,I),ANGLE)
           qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))= &
                          qmmm_struct%qm_coords(K,qm2_struct%peptide_links(J,I))-delAdj
           HEAT=HTYPE*SIN(ANGLE)**2         
           SUM=(REFH-HEAT)*TWOONEdelAdj
           dxyzqm(K,qm2_struct%peptide_links(J,I))=dxyzqm(K,qm2_struct%peptide_links(J,I))-SUM 
         end do
       end do                                   
     end do                                    
   end if                                           

end subroutine qm2_get_qm_forces

subroutine qm2_deriv_qm_analyt(qmmm_struct, iqm,jqm,loop_count,qm_qm_e_repul,PSUM, &
           n_atomic_orbj,n_atomic_orbi,corei,corej,betasas,betasap,betapas,betapap, &
           vec_qm_qm1,vec_qm_qm2, vec_qm_qm3, pair_force, &
           qqi, qqi2, qqj, qqj2, ddi, ddj, bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, &
           bdd3j,zz_sxs_over_sas,ss_eqn_adb,zz_sxp_over_sapij, &
           zz_sxp_over_sapji,sp_eqn_xx1ij,sp_eqn_xx2ij,sp_eqn_xx1ji,sp_eqn_xx2ji, &
           sp_eqn_xyij, sp_eqn_xyji,zz_pxp_over_pap,pp_eqn_xxy1,pp_eqn_xxy2)

!************************************************************************        
!*                                                                      *        
!*         CALCULATION OF ANALYTICAL DERIVATIVES                        *      
!*                                                                      *
!* Current routine maintained by Ross Walker (TSRI, 2004)               *
!* Inlining and optimisation by Ross Walker (TSRI, 2004)                *  
!*                                                                      *        
!************************************************************************        


!*** Variable Definitions:
!
!    qm_qm_e_repul = QM-QM electron repulsion integrals for this QM-QM pair
!             psum = Density matrix elements for orbitals centered on this
!                    QM-QM pair. - For coulomb term
!    n_atomic_orbj = number of atomic orbitals on atom j
!    n_atomic_orbi = number of atomic orbitals on atom i
!       vec_qm_qmx = xyz_qmi(x) - xyz_qmj(x) 
!       pair_force = Returned with cartesian force on this QM-QM pair

   use constants          , only: A_TO_BOHRS, A3_TO_BOHRS3, A2_TO_BOHRS2, AU_TO_EV, zero, &
                                  half, one, two, four, fourth, eighth, sixteenth, thirtysecond, &
                                  EV_TO_KCAL
   use ElementOrbitalIndex, only: NumberElements, MaxValenceOrbitals, MaxValenceDimension
   use qmmm_module        , only: qmmm_nml, qm2_rij_eqns, qm2_params, &
                                  AXIS_TOL, OVERLAP_CUTOFF, EXPONENTIAL_CUTOFF
   use qmmm_struct_module, only : qmmm_struct_type


   implicit none

!Passed in
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   integer, intent(in) :: iqm,jqm,loop_count
   _REAL_, intent(inout) :: qm_qm_e_repul(22) !for qmqm_erep_incore=true it is in, for =false it is out.
   _REAL_, intent(in) :: psum(MaxValenceOrbitals**2*3)
   integer, intent(in) :: n_atomic_orbj, n_atomic_orbi
   _REAL_, intent(in) :: corei,corej,betasas,betasap,betapas,betapap
   _REAL_, intent(in) :: vec_qm_qm1, vec_qm_qm2, vec_qm_qm3
   _REAL_, intent(in) :: qqi, qqi2, qqj, qqj2, ddi, ddj
   _REAL_, intent(in) :: bdd1i, bdd2i, bdd3i, bdd1j, bdd2j, bdd3j
   _REAL_, intent(out):: pair_force(3)
   _REAL_, intent(in) :: zz_sxs_over_sas(6,6),ss_eqn_adb(6,6)
   _REAL_, intent(in) :: zz_sxp_over_sapij(6,6), zz_sxp_over_sapji(6,6)
   _REAL_, intent(in) :: sp_eqn_xx1ij(6,6), sp_eqn_xx2ij(6,6)
   _REAL_, intent(in) :: sp_eqn_xx1ji(6,6), sp_eqn_xx2ji(6,6)
   _REAL_, intent(in) :: sp_eqn_xyij(6,6), sp_eqn_xyji(6,6)
   _REAL_, intent(in) :: zz_pxp_over_pap(6,6)
   _REAL_, intent(in) :: pp_eqn_xxy1(6,6), pp_eqn_xxy2(6,6)

!Local
   _REAL_ FAAX, FAAY, FAAZ,FABX, FABY, FABZ,FNUCX, FNUCY, FNUCZ
   _REAL_ dsx,dsy,dsz,DSPX(3), DSPY(3), DSPZ3
   _REAL_ DPXPX(3), DPYPY(3), DPZPZ(3), DPCROSS
   _REAL_ dgx(22), dgy(22), dgz(22)
   _REAL_ drx(MaxValenceDimension**2), dry(MaxValenceDimension**2), drz(MaxValenceDimension**2)
   _REAL_ r2, rij, onerij, rr, rr2
   _REAL_ bb, aa
   integer total_atomic_orb, qm_atomi_orb_start, orb_offset
   integer isp, k, l
   integer m, n, mn, kk, ll, kl, mk, nk, ml, nl

   _REAL_ vec_qm_qm_onerij1, vec_qm_qm_onerij2, vec_qm_qm_onerij3
   _REAL_ vec2_qm_qm1, vec2_qm_qm2, vec2_qm_qm3
   _REAL_ vec_qm_qm123, vec_qm_qm12, vec_qm_qm23, vec_qm_qm13
   integer i, j
   _REAL_ ABN, ADBR2, SS
   _REAL_ temp_real, temp_real2, temp_real3, temp_real4
   _REAL_  SQRTAEE, bdd1ij

!Originally in delri
   _REAL_ termx, termy, termz, ee, ade, aqe, dze, qzze, qxxe
   _REAL_ aed, aeq, edz, eqzz, eqxx
   _REAL_ add, adq, aqd, aqq, dxdx, dzdz,dzqxx,qxxdz,dzqzz,qzzdz
   _REAL_ qxxqxx, qxxqyy, qxxqzz, qzzqxx, qzzqzz, dxqxz, qxzdx, qxzqxz

!Originally in delmol
   integer mm, nn
   _REAL_ ::xTDX(MaxValenceOrbitals-1),xTDY(MaxValenceOrbitals-1),xTDZ(MaxValenceOrbitals-1)
   _REAL_ ::yTDX(MaxValenceOrbitals-1),yTDY(MaxValenceOrbitals-1),yTDZ(MaxValenceOrbitals-1)
   _REAL_ ::zTDX(MaxValenceOrbitals-1),zTDY(MaxValenceOrbitals-1),zTDZ(MaxValenceOrbitals-1)
   _REAL_ ::TX(MaxValenceOrbitals-1),TY(MaxValenceOrbitals-1),TZ(MaxValenceOrbitals-1)
   _REAL_ ::TZ3i, TZ3i2
   _REAL_ ::xtemp1, ytemp1, ztemp1, xtemp2, ytemp2, ztemp2

!Originally for rotat
   _REAL_ RXY2, RYZ2, RZX2, oneRXY
   logical LRXY2, LRYZ2, LRZX2

!For calls to core repulsion gradient routine
   _REAL_ :: xyzij(3), dxyz(3)
   _REAL_ :: gij, dgij(3)

!AWG: Now needed for call to qm2_repp
   _REAL_ :: core(10,2)

   FAAX=zero; FAAY=zero; FAAZ=zero
   FABX=zero; FABY=zero; FABZ=zero

   FNUCX = zero
   FNUCY = zero
   FNUCZ = zero

   if (n_atomic_orbi > 1 .OR. n_atomic_orbj >1) then      
     vec_qm_qm123 = vec_qm_qm1*vec_qm_qm2*vec_qm_qm3
     vec_qm_qm12 = vec_qm_qm1*vec_qm_qm2
     vec_qm_qm23 = vec_qm_qm2*vec_qm_qm3
     vec_qm_qm13 = vec_qm_qm1*vec_qm_qm3
   end if

   vec2_qm_qm1 = vec_qm_qm1*vec_qm_qm1 ! (xi - xj)**2
   vec2_qm_qm2 = vec_qm_qm2*vec_qm_qm2 ! (yi - yj)**2
   vec2_qm_qm3 = vec_qm_qm3*vec_qm_qm3 ! (zi - zj)**2

   r2=vec2_qm_qm1+vec2_qm_qm2+vec2_qm_qm3 ! (Ri-Rj)**2
   rr2=r2 * A2_TO_BOHRS2
   oneRIJ = one/sqrt(r2) !1/sqrt is faster than doing sqrt.
   RIJ = r2*oneRIJ !one/oneRIJ
   RR=RIJ * A_TO_BOHRS                                                                
   bdd1ij = bdd1i+bdd1j
   bdd1ij = bdd1ij*bdd1ij
   SQRTAEE=1.0d0/sqrt(RR2+bdd1ij)

   vec_qm_qm_onerij1 = vec_qm_qm1 * oneRIJ ! (xi - xj) / Rij
   vec_qm_qm_onerij2 = vec_qm_qm2 * oneRIJ ! (yi - yj) / Rij
   vec_qm_qm_onerij3 = vec_qm_qm3 * oneRIJ ! (zi - zj) / Rij

   total_atomic_orb=n_atomic_orbi+n_atomic_orbj

   qm_atomi_orb_start=n_atomic_orbj+1

!If we don't have the 1-e repul integrals for this QM-QM pair in memory we need
!to calculate them now.
   if (.NOT. qmmm_nml%qmqm_erep_incore) then
     call qm2_repp(qmmm_struct, iqm,jqm,rr,rr2,qm_qm_e_repul,core,SQRTAEE)
   end if

!   THE FIRST DERIVATIVES OF OVERLAP INTEGRALS                                  
      !Loop over atomic orbitals for QM atom i
      !for PX, PY, PZ - So K=0 = S orbital, K=1 = PX, 2 = PY, 3 = PZ
      ! CALCULATE OVERLAP DERIVATIVES, STORE RESULTS IN DS                     

   if (RR2 < OVERLAP_CUTOFF) then  !10A cutoff on overlaps
!ALL Atom pairs must have at least an S-S interaction - Min 1 orbital per atom.
! (S/S) Terms
     dsx=zero ; dsy=zero; dsz = zero
     do I=1,6 !1 to NGAUSS
       do J=1,6
         ADBR2=zz_sxs_over_sas(i,j)*RR2
!Check overlap is non-zero before starting
         if(ADBR2 < EXPONENTIAL_CUTOFF) then
           SS=ss_eqn_adb(i,j)*EXP(-ADBR2)
           DSx=DSx+vec_qm_qm1*SS
           DSy=DSy+vec_qm_qm2*SS
           DSz=DSz+vec_qm_qm3*SS
         end if
       end do
     end do

!    COMBINE TOGETHER THE OVERLAP DERIVATIVE PARTS
     orb_offset=1+qm2_params%pascal_tri1(qm_atomi_orb_start)
     temp_real=betasas*PSUM(orb_offset)
     FABx=FABx+temp_real*DSx
     FABy=FABy+temp_real*DSy
     FABz=FABz+temp_real*DSz

! end OF S/S terms, if both QM atoms have only 1 orbital then we don't do the loops below.

!NOW DO S with P orbitals if necessary (K=0)
     if (n_atomic_orbj>1) then
       DSPX=zero ; DSPY=zero !Zeros entire arrays of 3
       DSPZ3 = zero
       do j=1,6
         do i=1,6
           ADBR2=zz_sxp_over_sapij(i,j)*RR2
           !Check overlap is non-zero - otherwise we can skip this bit.
           if (ADBR2<EXPONENTIAL_CUTOFF) then
             ADBR2 = EXP(-ADBR2)
!   (S/PX-x) TERM 
             ABN=sp_eqn_xx1ij(i,j) - (sp_eqn_xx2ij(i,j)*vec2_qm_qm1)
             DSPX(1)=DSPX(1)+ABN*ADBR2
!   (S/PY-y) TERM
             ABN=sp_eqn_xx1ij(i,j) - (sp_eqn_xx2ij(i,j)*vec2_qm_qm2)
             DSPY(2)=DSPY(2)+ABN*ADBR2
!   (S/PZ-z) TERM
             ABN=sp_eqn_xx1ij(i,j) - (sp_eqn_xx2ij(i,j)*vec2_qm_qm3)
             DSPZ3=DSPZ3+ABN*ADBR2

             ABN=ADBR2*sp_eqn_xyij(i,j)
!   (PX-y/S) TERM
             DSPX(2)=DSPX(2)+ABN*vec_qm_qm12
!   (PX-z/S) TERM
             DSPX(3)=DSPX(3)+ABN*vec_qm_qm13
!   (PY-z/S) TERM
             DSPY(3)=DSPY(3)+ABN*vec_qm_qm23
!   (PY-x/S) TERM
             !DSPY(1)=DSPY(1)-ABN*vec_qm_qm12 DSPY(1)=DSPX(2)
!   (PZ-x/S) TERM
             !DSPZ(1)=DSPZ(1)-ABN*vec_qm_qm13 DSPZ(1)=DSPX(3)
!   (PZ-y/S) TERM
             !DSPZ(2)=DSPZ(2)-ABN*vec_qm_qm23 DSPZ(2)=DSPY(3)
           end if !(ADBR2<EXPONENTIAL_CUTOFF)
         end do
       end do
       !   COMBINE TOGETHER THE OVERLAP DERIVATIVE PARTS
       temp_real = betasap*PSUM(orb_offset+1)
       FABx=FABx+temp_real*DSPX(1)
       FABy=FABy+temp_real*DSPX(2)
       FABz=FABz+temp_real*DSPX(3)
       temp_real = betasap*PSUM(orb_offset+2)
       FABx=FABx+temp_real*DSPX(2) !DSPY(1)=DSPX(2)
       FABy=FABy+temp_real*DSPY(2)
       FABz=FABz+temp_real*DSPY(3)
       temp_real = betasap*PSUM(orb_offset+3)
       FABx=FABx+temp_real*DSPX(3) !DSPZ(1)=DSPX(3)
       FABy=FABy+temp_real*DSPY(3) !DSPZ(2)=DSPY(3)
       FABz=FABz+temp_real*DSPZ3
     end if ! if (n_atomic_orbj>1)

!Now do P orbitals of atom 1 with S of atom 2 (K>0 and L=0)
     if (n_atomic_orbi>1) then
       DSPX=zero ; DSPY=zero !Zeros entire arrays of 3
       DSPZ3=zero
       do I=1,6
         do J=1,6
           ADBR2=zz_sxp_over_sapji(j,i)*RR2
           if (ADBR2<EXPONENTIAL_CUTOFF) then !Only do the rest if the overlap is not zero
             ADBR2 = EXP(-ADBR2)
!   (PX-x/S) TERM
             ABN=-sp_eqn_xx1ji(j,i) + (sp_eqn_xx2ji(j,i)*vec2_qm_qm1)
             DSPX(1)=DSPX(1)+ABN*ADBR2
!   (PY-y/S) TERM
             ABN=-sp_eqn_xx1ji(j,i) + (sp_eqn_xx2ji(j,i)*vec2_qm_qm2)
             DSPY(2)=DSPY(2)+ABN*ADBR2
!   (PZ-z/S) TERM
             ABN=-sp_eqn_xx1ji(j,i) + (sp_eqn_xx2ji(j,i)*vec2_qm_qm3)
             DSPZ3=DSPZ3+ABN*ADBR2

             ABN=ADBR2*sp_eqn_xyji(j,i)
!   (PX-y/S) TERM
             DSPX(2)=DSPX(2)-ABN*vec_qm_qm12
!   (PX-z/S) TERM
             DSPX(3)=DSPX(3)-ABN*vec_qm_qm13
!   (PY-z/S) TERM
             DSPY(3)=DSPY(3)-ABN*vec_qm_qm23
!   (PY-x/S) TERM
             !DSPY(1)=DSPY(1)-ABN*vec_qm_qm12 DSPY(1)=DSPX(2)
!   (PZ-x/S) TERM
             !DSPZ(1)=DSPZ(1)-ABN*vec_qm_qm13 DSPZ(1)=DSPX(3)
!   (PZ-y/S) TERM
             !DSPZ(2)=DSPZ(2)-ABN*vec_qm_qm23 DSPZ(2)=DSPY(3)
           end if !(ADBR2<EXPONENTIAL_CUTOFF)
         end do
       end do
       temp_real = betapas*PSUM(1+qm2_params%pascal_tri1(qm_atomi_orb_start+1))
       FABx=FABx+temp_real*DSPX(1)
       FABy=FABy+temp_real*DSPX(2)
       FABz=FABz+temp_real*DSPX(3)
       temp_real = betapas*PSUM(1+qm2_params%pascal_tri1(qm_atomi_orb_start+2))
       FABx=FABx+temp_real*DSPX(2) !DSPY(1)=DSPX(2)
       FABy=FABy+temp_real*DSPY(2)
       FABz=FABz+temp_real*DSPY(3)
       temp_real = betapas*PSUM(1+qm2_params%pascal_tri1(qm_atomi_orb_start+3))
       FABx=FABx+temp_real*DSPX(3) !DSPZ(1)=DSPX(3)
       FABy=FABy+temp_real*DSPY(3) !DSPZ(2)=DSPY(3)
       FABz=FABz+temp_real*DSPZ3
     end if !if (n_atomic_orbi>1)

!Ross Walker - PP combinations that are the same have been condensed
!to one variable for speed and to save memory.

!Finally do P's of atom 1 with P's of atom 2 (if necessary)
     if (n_atomic_orbi > 1 .and. n_atomic_orbj > 1) then
       DPXPX=zero; DPYPY=zero; DPZPZ=zero !Zeros entire arrays of 3
       DPCROSS=zero
       do I=1,6
          do J=1,6
            ADBR2=zz_pxp_over_pap(i,j)*RR2
            if (ADBR2<EXPONENTIAL_CUTOFF) then !Only do the next bit if the overlap is not zero
              ADBR2=EXP(-ADBR2)
!    (PX / PX) - x
              ABN=((3.0D0*pp_eqn_xxy1(i,j)) + pp_eqn_xxy2(i,j) * vec2_qm_qm1)*ADBR2
              DPXPX(1)=DPXPX(1)+ABN*vec_qm_qm1
!    (PY / PY) - y
              ABN=((3.0D0*pp_eqn_xxy1(i,j)) + pp_eqn_xxy2(i,j) * vec2_qm_qm2)*ADBR2
              DPYPY(2)=DPYPY(2)+ABN*vec_qm_qm2
!    (PZ / PZ) - z
              ABN=((3.0D0*pp_eqn_xxy1(i,j)) + pp_eqn_xxy2(i,j) * vec2_qm_qm3)*ADBR2
              DPZPZ(3)=DPZPZ(3)+ABN*vec_qm_qm3

              ABN = (pp_eqn_xxy1(i,j) + pp_eqn_xxy2(i,j) * vec2_qm_qm1)*ADBR2
!    (PX / PX) - y
              DPXPX(2)=DPXPX(2)+ABN*vec_qm_qm2
!    (PX / PX) - z
              DPXPX(3)=DPXPX(3)+ABN*vec_qm_qm3
!    (PX / PY) - x
              !DPXPY(1)=DPXPY(1)+ABN*vec_qm_qm2 PXPY(1)=PXPX(2)
!    (PX / PZ) - x
              !DPXPZ(1)=DPXPZ(1)+ABN*vec_qm_qm3 PXPZ(1)=PXPX(3)
!    (PY / PX) - x
              !DPYPX(1)=DPYPX(1)+ABN*vec_qm_qm2 PXPY=PYPX
!    (PZ / PX) - x
              !DPZPX(1)=DPZPX(1)+ABN*vec_qm_qm3 PXPZ=PZPX

              ABN=(pp_eqn_xxy1(i,j) + pp_eqn_xxy2(i,j) * vec2_qm_qm2)*ADBR2
!    (PY / PY) - x
              DPYPY(1)=DPYPY(1)+ABN*vec_qm_qm1
!    (PY / PY) - z
              DPYPY(3)=DPYPY(3)+ABN*vec_qm_qm3
!    (PX / PY) - y
!              DPXPY(2)=DPXPY(2)+ABN*vec_qm_qm1 PXPY(2)=PYPX(2)=PYPY(1)
!    (PY / PZ) - y
!              DPYPZ(2)=DPYPZ(2)+ABN*vec_qm_qm3 PYPZ(2)=PZPY(2)=PYPY(3)
!    (PY / PX) - y
              !DPYPX(2)=DPYPX(2)+ABN*vec_qm_qm1 PXPY=PYPX
!    (PZ / PY) - y
              !DPZPY(2)=DPZPY(2)+ABN*vec_qm_qm3 PYPZ=PZPY

              ABN=(pp_eqn_xxy1(i,j) + pp_eqn_xxy2(i,j) * vec2_qm_qm3)*ADBR2
!    (PZ / PZ) - x
              DPZPZ(1)=DPZPZ(1)+ABN*vec_qm_qm1
!    (PZ / PZ) - y
              DPZPZ(2)=DPZPZ(2)+ABN*vec_qm_qm2
!    (PX / PZ) - z
!              DPXPZ(3)=DPXPZ(3)+ABN*vec_qm_qm1 PXPZ(3)=PZPX(3)=PZPZ(1)
!    (PY / PZ) - z
!              DPYPZ(3)=DPYPZ(3)+ABN*vec_qm_qm2 PYPZ(3)=PZPY(3)=PZPZ(2)
!    (PZ / PY) - z
              !DPZPY(3)=DPZPY(3)+ABN*vec_qm_qm2 PYPZ=PZPY
!    (PZ / PX) - z
              !DPZPX(3)=DPZPX(3)+ABN*vec_qm_qm1 PXPZ=PZPX

              ABN=pp_eqn_xxy2(i,j) * ADBR2*vec_qm_qm123
              DPCROSS=DPCROSS+ABN
!    (PX / PY) - z
              !DPXPY(3)=DPXPY(3)+ABN !PXPY(3)=PYPX(3)=DPCROSS
!    (PX / PZ) - y
              !DPXPZ(2)=DPXPZ(2)+ABN !PXPZ(2)=PZPX(2)=DPCROSS
!    (PY / PZ) - x
              !DPYPZ(1)=DPYPZ(1)+ABN !PYPZ(1)=PZPY(1)=DPCROSS
!    (PZ / PY) - x
              !DPZPY(1)=DPZPY(1)+ABN PYPZ=PZPY
!    (PY / PX) - z
              !DPYPX(3)=DPYPX(3)+ABN PXPY=PYPX
!    (PZ / PX) - y
              !DPZPX(2)=DPZPX(2)+ABN PXPZ=PZPX
           end if !ADBR2>EXPONENTIAL_CUTOFF
         end do
       end do
!      COMBINE TOGETHER THE OVERLAP DERIVATIVE PARTS
       orb_offset=2+qm2_params%pascal_tri1(qm_atomi_orb_start+1)
       temp_real = betapap*PSUM(orb_offset)
       FABx=FABx+temp_real*DPXPX(1)
       FABy=FABy+temp_real*DPXPX(2)
       FABz=FABz+temp_real*DPXPX(3)
       temp_real = betapap*PSUM(orb_offset+1)
       FABx=FABx+temp_real*DPXPX(2) !PYPX(1)=PXPY(1)=PXPX(2)
       FABy=FABy+temp_real*DPYPY(1) !PXPY(2)=PYPX(2)=PYPY(1)
       FABz=FABz+temp_real*DPCROSS  !PXPY(3)=PYPX(3)=DPCROSS
       temp_real = betapap*PSUM(orb_offset+2)
       FABx=FABx+temp_real*DPXPX(3) !PXPZ(1)=PZPX(1)=PXPX(3)
       FABy=FABy+temp_real*DPCROSS  !PXPZ(2)=PZPX(2)=DPCROSS
       FABz=FABz+temp_real*DPZPZ(1) !PXPZ(3)=PZPX(3)=PZPZ(1)
       orb_offset=2+qm2_params%pascal_tri1(qm_atomi_orb_start+2)
       temp_real=betapap*PSUM(orb_offset)
       FABx=FABx+temp_real*DPXPX(2) !PYPX(1)=PXPY(1)=PXPX(2)
       FABy=FABy+temp_real*DPYPY(1) !PXPY(2)=PYPX(2)=PYPY(1)
       FABz=FABz+temp_real*DPCROSS  !PXPY(3)=PYPX(3)=DPCROSS 
       temp_real=betapap*PSUM(orb_offset+1)
       FABx=FABx+temp_real*DPYPY(1)
       FABy=FABy+temp_real*DPYPY(2)
       FABz=FABz+temp_real*DPYPY(3)
       temp_real=betapap*PSUM(orb_offset+2)
       FABx=FABx+temp_real*DPCROSS  !PYPZ(1)=PZPY(1)=DPCROSS
       FABy=FABy+temp_real*DPYPY(3) !PYPZ(2)=PZPY(2)=PYPY(3)
       FABz=FABz+temp_real*DPZPZ(2) !PYPZ(3)=PZPY(3)=PZPZ(2)
       orb_offset=2+qm2_params%pascal_tri1(qm_atomi_orb_start+3)
       temp_real = betapap*PSUM(orb_offset)
       FABx=FABx+temp_real*DPXPX(3) !PXPZ(1)=PZPX(1)=PXPX(3)
       FABy=FABy+temp_real*DPCROSS  !PXPZ(2)=PZPX(2)=DPCROSS
       FABz=FABz+temp_real*DPZPZ(1) !PXPZ(3)=PZPX(3)=PZPZ(1)
       temp_real = betapap*PSUM(orb_offset+1)
       FABx=FABx+temp_real*DPCROSS  !PYPZ(1)=PZPY(1)=DPCROSS
       FABy=FABy+temp_real*DPYPY(3) !PYPZ(2)=PZPY(2)=PYPY(3)
       FABz=FABz+temp_real*DPZPZ(2) !PYPZ(3)=PZPY(3)=PZPZ(2)
       temp_real = betapap*PSUM(orb_offset+2)
       FABx=FABx+temp_real*DPZPZ(1)
       FABy=FABy+temp_real*DPZPZ(2)
       FABz=FABz+temp_real*DPZPZ(3)
     end if !n_atomic_orbi >1 .and. n_atomic_orbj > 1
   end if !if (RR2 < 100.D0 * A2_TO_BOHRS2) then  !10A cutoff on overlaps
!Code is linear from this point on for a given pair.

!NOTE: Ross Walker - In the equations below the only thing that varies for this pair
!during a simulation is the value of RR so we may want to consider pre-computing
!a lot of this. E.g. For a given pair AEE is constant.

!  S-orbital of QM-QM pairs - all QM pairs have s-s interactions.
   TERMX=vec_qm_qm_onerij1*AU_TO_EV*A_TO_BOHRS
   TERMY=vec_qm_qm_onerij2*AU_TO_EV*A_TO_BOHRS
   TERMZ=vec_qm_qm_onerij3*AU_TO_EV*A_TO_BOHRS
   EE    =-RR*(SQRTAEE)**3
   DGX(1)=TERMX*EE
   DGY(1)=TERMY*EE
   DGZ(1)=TERMZ*EE
   if(n_atomic_orbi > 2) then
!   SP-ATOM-HYDROGEN
      ADE=(bdd2i+bdd1j)**2

      AQE=(bdd3i+bdd1j)**2

      DZE   = (RR+ddi)/(SQRT((RR+ddi)**2+ADE))**3 &
              -(RR-ddi)/(SQRT((RR-ddi)**2+ADE))**3

      temp_real = (two*RR)/(SQRT(RR2+AQE))**3

      QZZE  =-(RR+two*qqi)/(SQRT((RR+two*qqi)**2+AQE))**3 &
             -(RR-two*qqi)/(SQRT((RR-two*qqi)**2+AQE))**3 &
             +temp_real

      QXXE  =-(two*RR)/(SQRT(RR2+four*qqi2+AQE))**3 &
             +temp_real

      temp_real = -DZE*half

      DGX(2)=TERMX*temp_real
      DGY(2)=TERMY*temp_real
      DGZ(2)=TERMZ*temp_real
      temp_real = EE+QZZE*fourth
      DGX(3)=TERMX*temp_real
      DGY(3)=TERMY*temp_real
      DGZ(3)=TERMZ*temp_real
      temp_real = EE+QXXE*fourth
      DGX(4)=TERMX*temp_real
      DGY(4)=TERMY*temp_real
      DGZ(4)=TERMZ*temp_real
   end if
   if (n_atomic_orbj > 2) then
!   HYDROGEN-SP-ATOM
      AED=(bdd1i+bdd2j)**2

      AEQ=(bdd1i+bdd3j)**2

      EDZ   = (RR-ddj)/(SQRT((RR-ddj)**2+AED))**3 &
             -(RR+ddj)/(SQRT((RR+ddj)**2+AED))**3

      temp_real = (two*RR)/(SQRT(RR2+AEQ))**3

      EQZZ  =-(RR-two*qqj)/(SQRT((RR-two*qqj)**2+AEQ))**3 &
             -(RR+two*qqj)/(SQRT((RR+two*qqj)**2+AEQ))**3 &
             +temp_real

      EQXX  =-(two*RR)/(SQRT(RR2+four*qqj2+AEQ))**3 &
             +temp_real

      temp_real = -EDZ*half
      DGX(5)=TERMX*temp_real
      DGY(5)=TERMY*temp_real
      DGZ(5)=TERMZ*temp_real
      temp_real = EE+EQZZ*fourth
      DGX(11)=TERMX*temp_real
      DGY(11)=TERMY*temp_real
      DGZ(11)=TERMZ*temp_real
      temp_real = EE+EQXX*fourth
      DGX(12)=TERMX*temp_real
      DGY(12)=TERMY*temp_real
      DGZ(12)=TERMZ*temp_real
   end if
   if (n_atomic_orbi > 2 .and. n_atomic_orbj > 2) then
!   SP-ATOM-SP-ATOM
      ADD=(bdd2i+bdd2j)**2
      ADQ=(bdd2i+bdd3j)**2
      AQD=(bdd3i+bdd2j)**2
      AQQ=(bdd3i+bdd3j)**2

      DXDX  =-(two*RR)/(SQRT(RR2+(ddi-ddj)**2+ADD))**3 &
             +(two*RR)/(SQRT(RR2+(ddi+ddj)**2+ADD))**3

      DZDZ  =-(RR+ddi-ddj)/(SQRT((RR+ddi-ddj)**2+ADD))**3 &
             -(RR-ddi+ddj)/(SQRT((RR-ddi+ddj)**2+ADD))**3 &
             +(RR-ddi-ddj)/(SQRT((RR-ddi-ddj)**2+ADD))**3 &
             +(RR+ddi+ddj)/(SQRT((RR+ddi+ddj)**2+ADD))**3

      DZQXX = two*(RR+ddi)/(SQRT((RR+ddi)**2+four*qqj2+ADQ))**3 &
             -two*(RR-ddi)/(SQRT((RR-ddi)**2+four*qqj2+ADQ))**3 &
             -two*(RR+ddi)/(SQRT((RR+ddi)**2+ADQ))**3 &
             +two*(RR-ddi)/(SQRT((RR-ddi)**2+ADQ))**3

      QXXDZ = two*(RR-ddj)/(SQRT((RR-ddj)**2+four*qqi2+AQD))**3 &
             -two*(RR+ddj)/(SQRT((RR+ddj)**2+four*qqi2+AQD))**3 &
             -two*(RR-ddj)/(SQRT((RR-ddj)**2+AQD))**3 &
             +two*(RR+ddj)/(SQRT((RR+ddj)**2+AQD))**3

      DZQZZ = (RR+ddi-two*qqj)/(SQRT((RR+ddi-two*qqj)**2+ADQ))**3 &
             -(RR-ddi-two*qqj)/(SQRT((RR-ddi-two*qqj)**2+ADQ))**3 &
             +(RR+ddi+two*qqj)/(SQRT((RR+ddi+two*qqj)**2+ADQ))**3 &
             -(RR-ddi+two*qqj)/(SQRT((RR-ddi+two*qqj)**2+ADQ))**3 &
             +two*(RR-ddi)/(SQRT((RR-ddi)**2+ADQ))**3 &
             -two*(RR+ddi)/(SQRT((RR+ddi)**2+ADQ))**3

      QZZDZ = (RR+two*qqi-ddj)/(SQRT((RR+two*qqi-ddj)**2+AQD))**3 &
             -(RR+two*qqi+ddj)/(SQRT((RR+two*qqi+ddj)**2+AQD))**3 &
             +(RR-two*qqi-ddj)/(SQRT((RR-two*qqi-ddj)**2+AQD))**3 &
             -(RR-two*qqi+ddj)/(SQRT((RR-two*qqi+ddj)**2+AQD))**3 &
             -two*(RR-ddj)/(SQRT((RR-ddj)**2+AQD))**3 &
             +two*(RR+ddj)/(SQRT((RR+ddj)**2+AQD))**3

      QXXQXX=-(two*RR)/(SQRT(RR2+four*(qqi-qqj)**2+AQQ))**3 &
             -(two*RR)/(SQRT(RR2+four*(qqi+qqj)**2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqi2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqj2+AQQ))**3 &
             -(four*RR)/(SQRT(RR2+AQQ))**3

      QXXQYY=-(four*RR)/(SQRT(RR2+four*qqi2+four*qqj2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqi2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqj2+AQQ))**3 &
             -(four*RR)/(SQRT(RR2+AQQ))**3

      QXXQZZ= &
             -two*(RR-two*qqj)/(SQRT((RR-two*qqj)**2+four*qqi2+AQQ))**3 &
             -two*(RR+two*qqj)/(SQRT((RR+two*qqj)**2+four*qqi2+AQQ))**3 &
             +two*(RR-two*qqj)/(SQRT((RR-two*qqj)**2+AQQ))**3 &
             +two*(RR+two*qqj)/(SQRT((RR+two*qqj)**2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqi2+AQQ))**3 &
             -(four*RR)/(SQRT(RR2+AQQ))**3

      QZZQXX= &
             -two*(RR+two*qqi)/(SQRT((RR+two*qqi)**2+four*qqj2+AQQ))**3 &
             -two*(RR-two*qqi)/(SQRT((RR-two*qqi)**2+four*qqj2+AQQ))**3 &
             +two*(RR+two*qqi)/(SQRT((RR+two*qqi)**2+AQQ))**3 &
             +two*(RR-two*qqi)/(SQRT((RR-two*qqi)**2+AQQ))**3 &
             +(four*RR)/(SQRT(RR2+four*qqj2+AQQ))**3 &
             -(four*RR)/(SQRT(RR2+AQQ))**3

      QZZQZZ= &
           -(RR+two*qqi-two*qqj)/(SQRT((RR+two*qqi-two*qqj)**2+AQQ))**3 &
           -(RR+two*qqi+two*qqj)/(SQRT((RR+two*qqi+two*qqj)**2+AQQ))**3 &
           -(RR-two*qqi-two*qqj)/(SQRT((RR-two*qqi-two*qqj)**2+AQQ))**3 &
           -(RR-two*qqi+two*qqj)/(SQRT((RR-two*qqi+two*qqj)**2+AQQ))**3 &
             +two*(RR-two*qqi)/(SQRT((RR-two*qqi)**2+AQQ))**3 &
             +two*(RR+two*qqi)/(SQRT((RR+two*qqi)**2+AQQ))**3 &
             +two*(RR-2.D0*qqj)/(SQRT((RR-2.D0*qqj)**2+AQQ))**3 &
             +two*(RR+2.D0*qqj)/(SQRT((RR+2.D0*qqj)**2+AQQ))**3 &
             -(four*RR)/(SQRT(RR2+AQQ))**3

      DXQXZ = two*(RR-qqj)/(SQRT((RR-qqj)**2+(ddi-qqj)**2+ADQ))**3 &
             -two*(RR+qqj)/(SQRT((RR+qqj)**2+(ddi-qqj)**2+ADQ))**3 &
             -two*(RR-qqj)/(SQRT((RR-qqj)**2+(ddi+qqj)**2+ADQ))**3 &
             +two*(RR+qqj)/(SQRT((RR+qqj)**2+(ddi+qqj)**2+ADQ))**3

      QXZDX = two*(RR+qqi)/(SQRT((RR+qqi)**2+(qqi-ddj)**2+AQD))**3 &
             -two*(RR-qqi)/(SQRT((RR-qqi)**2+(qqi-ddj)**2+AQD))**3 &
             -two*(RR+qqi)/(SQRT((RR+qqi)**2+(qqi+ddj)**2+AQD))**3 &
             +two*(RR-qqi)/(SQRT((RR-qqi)**2+(qqi+ddj)**2+AQD))**3

      QXZQXZ=-two*(RR+qqi-qqj)/(SQRT((RR+qqi-qqj)**2+(qqi-qqj)**2+AQQ))**3 &
             +two*(RR+qqi+qqj)/(SQRT((RR+qqi+qqj)**2+(qqi-qqj)**2+AQQ))**3 &
             +two*(RR-qqi-qqj)/(SQRT((RR-qqi-qqj)**2+(qqi-qqj)**2+AQQ))**3 &
             -two*(RR-qqi+qqj)/(SQRT((RR-qqi+qqj)**2+(qqi-qqj)**2+AQQ))**3 &
             +two*(RR+qqi-qqj)/(SQRT((RR+qqi-qqj)**2+(qqi+qqj)**2+AQQ))**3 &
             -two*(RR+qqi+qqj)/(SQRT((RR+qqi+qqj)**2+(qqi+qqj)**2+AQQ))**3 &
             -two*(RR-qqi-qqj)/(SQRT((RR-qqi-qqj)**2+(qqi+qqj)**2+AQQ))**3 &
             +two*(RR-qqi+qqj)/(SQRT((RR-qqi+qqj)**2+(qqi+qqj)**2+AQQ))**3

      temp_real = DZDZ*fourth
      DGX(6)=TERMX*temp_real
      DGY(6)=TERMY*temp_real
      DGZ(6)=TERMZ*temp_real
      temp_real = DXDX*fourth
      DGX(7)=TERMX*temp_real
      DGY(7)=TERMY*temp_real
      DGZ(7)=TERMZ*temp_real
      temp_real = -(EDZ*half+QZZDZ*eighth)
      DGX(8)=TERMX*temp_real
      DGY(8)=TERMY*temp_real
      DGZ(8)=TERMZ*temp_real
      temp_real = -(EDZ*half+QXXDZ*eighth)
      DGX(9)=TERMX*temp_real
      DGY(9)=TERMY*temp_real
      DGZ(9)=TERMZ*temp_real
      temp_real = -QXZDX*eighth
      DGX(10)=TERMX*temp_real
      DGY(10)=TERMY*temp_real
      DGZ(10)=TERMZ*temp_real
      temp_real = -(DZE*half+DZQZZ*eighth) 
      DGX(13)=TERMX*temp_real
      DGY(13)=TERMY*temp_real
      DGZ(13)=TERMZ*temp_real
      temp_real = -(DZE*half+DZQXX*eighth)
      DGX(14)=TERMX*temp_real
      DGY(14)=TERMY*temp_real
      DGZ(14)=TERMZ*temp_real
      temp_real = -DXQXZ*eighth
      DGX(15)=TERMX*temp_real
      DGY(15)=TERMY*temp_real
      DGZ(15)=TERMZ*temp_real
      temp_real2 = EE+EQZZ*fourth
      temp_real = temp_real2+QZZE*fourth+QZZQZZ*sixteenth
      DGX(16)=TERMX*temp_real
      DGY(16)=TERMY*temp_real
      DGZ(16)=TERMZ*temp_real
      temp_real = temp_real2+QXXE*fourth+QXXQZZ*sixteenth
      DGX(17)=TERMX*temp_real
      DGY(17)=TERMY*temp_real
      DGZ(17)=TERMZ*temp_real
      temp_real2 = EE+EQXX*fourth
      temp_real = temp_real2+QZZE*fourth+QZZQXX*sixteenth
      DGX(18)=TERMX*temp_real
      DGY(18)=TERMY*temp_real
      DGZ(18)=TERMZ*temp_real
      temp_real = temp_real2+QXXE*fourth+QXXQXX*sixteenth
      DGX(19)=TERMX*temp_real
      DGY(19)=TERMY*temp_real
      DGZ(19)=TERMZ*temp_real
      temp_real = QXZQXZ*sixteenth
      DGX(20)=TERMX*temp_real
      DGY(20)=TERMY*temp_real
      DGZ(20)=TERMZ*temp_real
      temp_real = temp_real2+QXXE*fourth+QXXQYY*sixteenth
      DGX(21)=TERMX*temp_real
      DGY(21)=TERMY*temp_real
      DGZ(21)=TERMZ*temp_real
      temp_real = (QXXQXX-QXXQYY)*thirtysecond
      DGX(22)=TERMX*temp_real
      DGY(22)=TERMY*temp_real
      DGZ(22)=TERMZ*temp_real
   end if

   if(n_atomic_orbi>1.OR.n_atomic_orbj>1) then
      RXY2 = vec2_qm_qm1 + vec2_qm_qm2
      RYZ2 = vec2_qm_qm2 + vec2_qm_qm3
      RZX2 = vec2_qm_qm3 + vec2_qm_qm1
      LRXY2 = RXY2 < AXIS_TOL
      LRYZ2 = RYZ2 < AXIS_TOL
      LRZX2 = RZX2 < AXIS_TOL

      XTDX=zero; YTDX=zero; ZTDX=zero
      XTDY=zero; YTDY=zero; ZTDY=zero
      XTDZ=zero; YTDZ=zero; ZTDZ=zero

      if (.NOT.(LRXY2 .OR. LRYZ2 .OR. LRZX2)) then
        TERMX = vec_qm_qm_onerij1 * oneRIJ !vec_qm_qm1/(RIJ*RIJ)
        TERMY = vec_qm_qm_onerij2 * oneRIJ
        TERMZ = vec_qm_qm_onerij3 * oneRIJ
        oneRXY  = one/sqrt(RXY2)
        !Ross Walker - rearranged order here slightly for speed.
        TZ(3) = oneRIJ/oneRXY
        TZ3i = RIJ * oneRXY !Inverse of TZ(3) to avoid other divisions.
        TZ3i2 = TZ3i*TZ3i   !Square of 1/TZ(3)

        TX(1) = vec_qm_qm_onerij1
        TX(2) = vec_qm_qm_onerij2
        TX(3) = vec_qm_qm_onerij3

        TY(1) = -TX(2)*SIGN(one,TX(1)) * TZ3i
        TY(2) = ABS(TX(1) * TZ3i)
        TY(3) = zero

        TZ(1) = -TX(1)*TX(3)*TZ3i
        TZ(2) = -TX(2)*TX(3)*TZ3i
        !TZ(3) = see above

        XTDX(1)=oneRIJ-TX(1)*TERMX
        YTDX(1)=-TX(1)*TERMY
        ZTDX(1)=-TX(1)*TERMZ

        XTDX(2)=-TX(2)*TERMX
        YTDX(2)=oneRIJ-TX(2)*TERMY
        ZTDX(2)=-TX(2)*TERMZ

        XTDX(3)=-TX(3)*TERMX
        YTDX(3)=-TX(3)*TERMY
        ZTDX(3)=oneRIJ-TX(3)*TERMZ

        XTDZ(3)=TX(1)*oneRXY-TZ(3)*TERMX
        YTDZ(3)=TX(2)*oneRXY-TZ(3)*TERMY
        ZTDZ(3)=-TZ(3)*TERMZ

        XTDY(1)=-XTDX(2)*TZ3i+TX(2)*XTDZ(3)*TZ3i2
        YTDY(1)=-YTDX(2)*TZ3i+TX(2)*YTDZ(3)*TZ3i2
        ZTDY(1)=-ZTDX(2)*TZ3i+TX(2)*ZTDZ(3)*TZ3i2

        XTDY(1)=XTDY(1)*sign(one,TX(1))
        YTDY(1)=YTDY(1)*sign(one,TX(1))
        ZTDY(1)=ZTDY(1)*sign(one,TX(1))

        XTDY(2)=XTDX(1)*TZ3I-TX(1)*XTDZ(3)*TZ3I2
        YTDY(2)=YTDX(1)*TZ3I-TX(1)*YTDZ(3)*TZ3I2
        ZTDY(2)=ZTDX(1)*TZ3I-TX(1)*ZTDZ(3)*TZ3I2

        XTDY(2)=XTDY(2)*sign(one,TX(1))
        YTDY(2)=YTDY(2)*sign(one,TX(1))
        ZTDY(2)=ZTDY(2)*sign(one,TX(1))

! Don't need to zero these again as they were zeroed above
!       XTDY(3)=0.0D0
!       YTDY(3)=0.0D0
!       ZTDY(3)=0.0D0

!Note: Ross Walker and Mike Crowley, we could factor out TZ3I here or we could
!pre-compute -TX(3)*TZ3I etc etc. But for the moment we will leave it as is since
!this is really just doing the compiler's work.
         XTDZ(1)=-TX(3)*XTDX(1)*TZ3I-TX(1)*XTDX(3)*TZ3I &
                +TX(1)*TX(3)*XTDZ(3)*TZ3I2
         YTDZ(1)=-TX(3)*YTDX(1)*TZ3I-TX(1)*YTDX(3)*TZ3I &
                +TX(1)*TX(3)*YTDZ(3)*TZ3I2
         ZTDZ(1)=-TX(3)*ZTDX(1)*TZ3I-TX(1)*ZTDX(3)*TZ3I &
                +TX(1)*TX(3)*ZTDZ(3)*TZ3I2

         XTDZ(2)=-TX(3)*XTDX(2)*TZ3I-TX(2)*XTDX(3)*TZ3I &
                +TX(2)*TX(3)*XTDZ(3)*TZ3I2
         YTDZ(2)=-TX(3)*YTDX(2)*TZ3I-TX(2)*YTDX(3)*TZ3I &
                +TX(2)*TX(3)*YTDZ(3)*TZ3I2
         ZTDZ(2)=-TX(3)*ZTDX(2)*TZ3I-TX(2)*ZTDX(3)*TZ3I &
                +TX(2)*TX(3)*ZTDZ(3)*TZ3I2
      elseif (LRXY2) then
!      MOLECULAR Z AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(1)=zero
         TX(2)=zero
         TX(3)=sign(one,vec_qm_qm3)
         TY(1)=zero
         TY(2)=one
         TY(3)=zero
         TZ(1)=TX(3)
         TZ(2)=zero
         TZ(3)=zero
       !X Axis
         XTDX(1)=oneRIJ
         XTDZ(3)=-oneRIJ
       !Y Axis
         YTDX(2)=oneRIJ
         YTDY(3)=-TX(3)*oneRIJ
       !Z Axis
      elseif (LRYZ2) then
!      MOLECULAR X AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(1)=sign(one,vec_qm_qm1)
         TX(2)=zero
         TX(3)=zero
         TY(1)=zero
         TY(2)=TX(1)
         TY(3)=zero
         TZ(1)=zero
         TZ(2)=zero
         TZ(3)=one
         !X Axis
         !Y Axis
         YTDX(2)=oneRIJ
         YTDY(1)=-oneRIJ
      !Z Axis
         ZTDX(3)=oneRIJ
         ZTDZ(1)=-TX(1)*oneRIJ
      else !if (LRZX2) then
!      MOLECULAR Y AXIS IS PARALLEL TO DIATOMIC Z AXIS
         TX(1)=zero
         TX(2)=sign(one,vec_qm_qm2)
         TX(3)=zero
         TY(1)=-TX(2)
         TY(2)=zero
         TY(3)=zero
         TZ(1)=zero
         TZ(2)=zero
         TZ(3)=one
       !X Axis
         XTDX(1)=oneRIJ
         XTDY(2)=oneRIJ
       !Y Axis
       !Z Axis
         ZTDX(3)=oneRIJ
         ZTDZ(2)=-TX(2)*oneRIJ
      end if
      !At least one atom has more than one atomic orbital so need to consider S and P interactions.
      isp = 0
      do K=qm_atomi_orb_start,total_atomic_orb
         KK=K-qm_atomi_orb_start
         do L=K,total_atomic_orb
            LL=L-qm_atomi_orb_start
            do M=1,n_atomic_orbj
               MM=M-1
               do N=M,n_atomic_orbj
                  NN=N-1
                  ISP=ISP+1
                  if(NN == 0)then
                     if(LL == 0) then
!   (SS/SS)
                        DRX(ISP)=DGX(1)
                        DRY(ISP)=DGY(1)
                        DRZ(ISP)=DGZ(1)
                    elseif(KK == 0) then
!   (SP/SS)
                        DRX(ISP)=DGX(2)*TX(LL)+qm_qm_e_repul(2)*xTDX(LL)
                        DRY(ISP)=DGY(2)*TX(LL)+qm_qm_e_repul(2)*yTDX(LL)
                        DRZ(ISP)=DGZ(2)*TX(LL)+qm_qm_e_repul(2)*zTDX(LL)
                     else
!   (PP/SS)
                        DRX(ISP)=DGX(3)*TX(KK)*TX(LL)+qm_qm_e_repul(3)*(xTDX(KK)    &
                                *TX(LL)+TX(KK)*xTDX(LL))+DGX(4)*(TY(KK)*TY(LL)     &
                                +TZ(KK)*TZ(LL))+qm_qm_e_repul(4)*(xTDY(KK)*TY(LL) &
                                +TY(KK)*xTDY(LL)+xTDZ(KK)*TZ(LL)+TZ(KK)*xTDZ(LL))
                        DRY(ISP)=DGY(3)*TX(KK)*TX(LL)+qm_qm_e_repul(3)*(yTDX(KK)    &
                                *TX(LL)+TX(KK)*yTDX(LL))+DGY(4)*(TY(KK)*TY(LL)     &
                                +TZ(KK)*TZ(LL))+qm_qm_e_repul(4)*(yTDY(KK)*TY(LL) &
                                +TY(KK)*yTDY(LL)+yTDZ(KK)*TZ(LL)+TZ(KK)*yTDZ(LL))
                        DRZ(ISP)=DGZ(3)*TX(KK)*TX(LL)+qm_qm_e_repul(3)*(zTDX(KK)    &
                                *TX(LL)+TX(KK)*zTDX(LL))+DGZ(4)*(TY(KK)*TY(LL)     &
                                +TZ(KK)*TZ(LL))+qm_qm_e_repul(4)*(zTDY(KK)*TY(LL) &
                                +TY(KK)*zTDY(LL)+zTDZ(KK)*TZ(LL)+TZ(KK)*zTDZ(LL))
                     endif
                 elseif(MM == 0) then
                     if(LL == 0) then
!   (SS/SP)
                        DRX(ISP)=DGX(5)*TX(NN)+qm_qm_e_repul(5)*xTDX(NN)
                        DRY(ISP)=DGY(5)*TX(NN)+qm_qm_e_repul(5)*yTDX(NN)
                        DRZ(ISP)=DGZ(5)*TX(NN)+qm_qm_e_repul(5)*zTDX(NN)
                    elseif(KK == 0) then
!   (SP/SP)
                        DRX(ISP)=DGX(6)*TX(LL)*TX(NN)+qm_qm_e_repul(6)*(xTDX(LL)*TX(NN) &
                                +TX(LL)*xTDX(NN))+DGX(7)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN)) &
                                +qm_qm_e_repul(7)*(xTDY(LL)*TY(NN)+TY(LL)*xTDY(NN)     &
                                +xTDZ(LL)*TZ(NN)+TZ(LL)*xTDZ(NN))
                        DRY(ISP)=DGY(6)*TX(LL)*TX(NN)+qm_qm_e_repul(6)*(yTDX(LL)*TX(NN) &
                                +TX(LL)*yTDX(NN))+DGY(7)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN)) &
                                +qm_qm_e_repul(7)*(yTDY(LL)*TY(NN)+TY(LL)*yTDY(NN)     &
                                +yTDZ(LL)*TZ(NN)+TZ(LL)*yTDZ(NN))
                        DRZ(ISP)=DGZ(6)*TX(LL)*TX(NN)+qm_qm_e_repul(6)*(zTDX(LL)*TX(NN) &
                                +TX(LL)*zTDX(NN))+DGZ(7)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN)) &
                                +qm_qm_e_repul(7)*(zTDY(LL)*TY(NN)+TY(LL)*zTDY(NN)     &
                                +zTDZ(LL)*TZ(NN)+TZ(LL)*zTDZ(NN))
                     else
!   (PP/SP)
                        DRX(ISP)=DGX(8)*TX(KK)*TX(LL)*TX(NN)+qm_qm_e_repul(8)*(xTDX(KK)  &
                                *TX(LL)*TX(NN)+TX(KK)*xTDX(LL)*TX(NN)+TX(KK)*TX(LL)    &
                                *xTDX(NN))+DGX(9)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(NN)  &
                                +qm_qm_e_repul(9)*((xTDY(KK)*TY(LL)+TY(KK)*xTDY(LL)     &
                                +xTDZ(KK)*TZ(LL)+TZ(KK)*xTDZ(LL))*TX(NN)+(TY(KK)*TY(LL) &
                                +TZ(KK)*TZ(LL))*xTDX(NN))+DGX(10)*(TX(KK)*(TY(LL)*TY(NN)&
                                +TZ(LL)*TZ(NN))+TX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))) &
                                +qm_qm_e_repul(10)*(xTDX(KK)*(TY(LL)*TY(NN)+TZ(LL)     &
                                *TZ(NN))+xTDX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(KK) &
                                *(xTDY(LL)*TY(NN)+TY(LL)*xTDY(NN)+xTDZ(LL)*TZ(NN)+TZ(LL) &
                                *xTDZ(NN))+TX(LL)*(xTDY(KK)*TY(NN)+TY(KK)*xTDY(NN)       &
                                +xTDZ(KK)*TZ(NN)+TZ(KK)*xTDZ(NN)))
                        DRY(ISP)=DGY(8)*TX(KK)*TX(LL)*TX(NN)+qm_qm_e_repul(8)*(yTDX(KK)  &
                                *TX(LL)*TX(NN)+TX(KK)*yTDX(LL)*TX(NN)+TX(KK)*TX(LL)    &
                                *yTDX(NN))+DGY(9)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(NN)  &
                                +qm_qm_e_repul(9)*((yTDY(KK)*TY(LL)+TY(KK)*yTDY(LL)     &
                                +yTDZ(KK)*TZ(LL)+TZ(KK)*yTDZ(LL))*TX(NN)+(TY(KK)*TY(LL) &
                                +TZ(KK)*TZ(LL))*yTDX(NN))+DGY(10)*(TX(KK)*(TY(LL)*TY(NN)&
                                +TZ(LL)*TZ(NN))+TX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))) &
                                +qm_qm_e_repul(10)*(yTDX(KK)*(TY(LL)*TY(NN)+TZ(LL)     &
                                *TZ(NN))+yTDX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(KK) &
                                *(yTDY(LL)*TY(NN)+TY(LL)*yTDY(NN)+yTDZ(LL)*TZ(NN)+TZ(LL) &
                                *yTDZ(NN))+TX(LL)*(yTDY(KK)*TY(NN)+TY(KK)*yTDY(NN)       &
                                +yTDZ(KK)*TZ(NN)+TZ(KK)*yTDZ(NN)))
                        DRZ(ISP)=DGZ(8)*TX(KK)*TX(LL)*TX(NN)+qm_qm_e_repul(8)*(zTDX(KK)  &
                                *TX(LL)*TX(NN)+TX(KK)*zTDX(LL)*TX(NN)+TX(KK)*TX(LL)    &
                                *zTDX(NN))+DGZ(9)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(NN)  &
                                +qm_qm_e_repul(9)*((zTDY(KK)*TY(LL)+TY(KK)*zTDY(LL)     &
                                +zTDZ(KK)*TZ(LL)+TZ(KK)*zTDZ(LL))*TX(NN)+(TY(KK)*TY(LL) &
                                +TZ(KK)*TZ(LL))*zTDX(NN))+DGZ(10)*(TX(KK)*(TY(LL)*TY(NN)&
                                +TZ(LL)*TZ(NN))+TX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))) &
                                +qm_qm_e_repul(10)*(zTDX(KK)*(TY(LL)*TY(NN)+TZ(LL)     &
                                *TZ(NN))+zTDX(LL)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(KK) &
                                *(zTDY(LL)*TY(NN)+TY(LL)*zTDY(NN)+zTDZ(LL)*TZ(NN)+TZ(LL) &
                                *zTDZ(NN))+TX(LL)*(zTDY(KK)*TY(NN)+TY(KK)*zTDY(NN)       &
                                +zTDZ(KK)*TZ(NN)+TZ(KK)*zTDZ(NN)))
                     endif
                 elseif(LL == 0) then
!   (SS/PP)
                     DRX(ISP)=DGX(11)*TX(MM)*TX(NN)+qm_qm_e_repul(11)*(xTDX(MM)*TX(NN)+TX(MM) &
                             *xTDX(NN))+DGX(12)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))                &
                             +qm_qm_e_repul(12)*(xTDY(MM)*TY(NN)+TY(MM)*xTDY(NN)+xTDZ(MM)     &
                             *TZ(NN)+TZ(MM)*xTDZ(NN))
                     DRY(ISP)=DGY(11)*TX(MM)*TX(NN)+qm_qm_e_repul(11)*(yTDX(MM)*TX(NN)+TX(MM) &
                             *yTDX(NN))+DGY(12)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))                &
                             +qm_qm_e_repul(12)*(yTDY(MM)*TY(NN)+TY(MM)*yTDY(NN)+yTDZ(MM)     &
                             *TZ(NN)+TZ(MM)*yTDZ(NN))
                     DRZ(ISP)=DGZ(11)*TX(MM)*TX(NN)+qm_qm_e_repul(11)*(zTDX(MM)*TX(NN)+TX(MM) &
                             *zTDX(NN))+DGZ(12)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))                &
                             +qm_qm_e_repul(12)*(zTDY(MM)*TY(NN)+TY(MM)*zTDY(NN)+zTDZ(MM)     &
                             *TZ(NN)+TZ(MM)*zTDZ(NN))
                 elseif(KK == 0) then
!   (SP/PP)
                     DRX(ISP)=DGX(13)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(13)*(xTDX(LL)*TX(MM) &
                             *TX(NN)+TX(LL)*xTDX(MM)*TX(NN)+TX(LL)*TX(MM)*xTDX(NN))+DGX(14)   &
                             *TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(14)       &
                             *(xTDX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(LL)*(xTDY(MM)*TY(NN)&
                             +TY(MM)*xTDY(NN)+xTDZ(MM)*TZ(NN)+TZ(MM)*xTDZ(NN)))+DGX(15)*(TY(LL)&
                             *(TY(MM)*TX(NN)+TY(NN)*TX(MM))+TZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)   &
                             *TX(MM)))+qm_qm_e_repul(15)*(xTDY(LL)*(TY(MM)*TX(NN)+TY(NN)    &
                             *TX(MM))+xTDZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM))+TY(LL)*(xTDY(MM)&
                             *TX(NN)+TY(MM)*xTDX(NN)+xTDY(NN)*TX(MM)+TY(NN)*xTDX(MM))+TZ(LL)  &
                             *(xTDZ(MM)*TX(NN)+TZ(MM)*xTDX(NN)+xTDZ(NN)*TX(MM)+TZ(NN)*xTDX(MM)))
                     DRY(ISP)=DGY(13)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(13)*(yTDX(LL)*TX(MM) &
                             *TX(NN)+TX(LL)*yTDX(MM)*TX(NN)+TX(LL)*TX(MM)*yTDX(NN))+DGY(14)   &
                             *TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(14)       &
                             *(yTDX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(LL)*(yTDY(MM)*TY(NN)&
                             +TY(MM)*yTDY(NN)+yTDZ(MM)*TZ(NN)+TZ(MM)*yTDZ(NN)))+DGY(15)*(TY(LL)&
                             *(TY(MM)*TX(NN)+TY(NN)*TX(MM))+TZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)   &
                             *TX(MM)))+qm_qm_e_repul(15)*(yTDY(LL)*(TY(MM)*TX(NN)+TY(NN)    &
                             *TX(MM))+yTDZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM))+TY(LL)*(yTDY(MM)&
                             *TX(NN)+TY(MM)*yTDX(NN)+yTDY(NN)*TX(MM)+TY(NN)*yTDX(MM))+TZ(LL)  &
                             *(yTDZ(MM)*TX(NN)+TZ(MM)*yTDX(NN)+yTDZ(NN)*TX(MM)+TZ(NN)*yTDX(MM)))
                     DRZ(ISP)=DGZ(13)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(13)*(zTDX(LL)*TX(MM) &
                             *TX(NN)+TX(LL)*zTDX(MM)*TX(NN)+TX(LL)*TX(MM)*zTDX(NN))+DGZ(14)   &
                             *TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(14)       &
                             *(zTDX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(LL)*(zTDY(MM)*TY(NN)&
                             +TY(MM)*zTDY(NN)+zTDZ(MM)*TZ(NN)+TZ(MM)*zTDZ(NN)))+DGZ(15)*(TY(LL)&
                             *(TY(MM)*TX(NN)+TY(NN)*TX(MM))+TZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)   &
                             *TX(MM)))+qm_qm_e_repul(15)*(zTDY(LL)*(TY(MM)*TX(NN)+TY(NN)    &
                             *TX(MM))+zTDZ(LL)*(TZ(MM)*TX(NN)+TZ(NN)*TX(MM))+TY(LL)*(zTDY(MM)&
                             *TX(NN)+TY(MM)*zTDX(NN)+zTDY(NN)*TX(MM)+TY(NN)*zTDX(MM))+TZ(LL)  &
                             *(zTDZ(MM)*TX(NN)+TZ(MM)*zTDX(NN)+zTDZ(NN)*TX(MM)+TZ(NN)*zTDX(MM)))
                  else
!   (PP/PP)
                     DRX(ISP)=DGX(16)*TX(KK)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(16)*(xTDX(KK)*TX(LL)*TX(MM)     &
                             *TX(NN)+TX(KK)*xTDX(LL)*TX(MM)*TX(NN)+TX(KK)*TX(LL)*xTDX(MM)*TX(NN)+TX(KK)*TX(LL) &
                             *TX(MM)*xTDX(NN))+DGX(17)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(MM)*TX(NN)             &
                             +qm_qm_e_repul(17)*((xTDY(KK)*TY(LL)+TY(KK)*xTDY(LL)+xTDZ(KK)*TZ(LL)+TZ(KK)        &
                             *xTDZ(LL))*TX(MM)*TX(NN)+(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*(xTDX(MM)*TX(NN)+TX(MM)    &
                             *xTDX(NN)))+DGX(18)*TX(KK)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(18) &
                             *((xTDX(KK)*TX(LL)+TX(KK)*xTDX(LL))*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(KK)*TX(LL)   &
                             *(xTDY(MM)*TY(NN)+TY(MM)*xTDY(NN)+xTDZ(MM)*TZ(NN)+TZ(MM)*xTDZ(NN)))
                     DRY(ISP)=DGY(16)*TX(KK)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(16)*(yTDX(KK)*TX(LL)*TX(MM)     &
                             *TX(NN)+TX(KK)*yTDX(LL)*TX(MM)*TX(NN)+TX(KK)*TX(LL)*yTDX(MM)*TX(NN)+TX(KK)*TX(LL) &
                             *TX(MM)*yTDX(NN))+DGY(17)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(MM)*TX(NN)             &
                             +qm_qm_e_repul(17)*((yTDY(KK)*TY(LL)+TY(KK)*yTDY(LL)+yTDZ(KK)*TZ(LL)+TZ(KK)        &
                             *yTDZ(LL))*TX(MM)*TX(NN)+(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*(yTDX(MM)*TX(NN)+TX(MM)    &
                             *yTDX(NN)))+DGY(18)*TX(KK)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(18) &
                             *((yTDX(KK)*TX(LL)+TX(KK)*yTDX(LL))*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(KK)*TX(LL)   &
                             *(yTDY(MM)*TY(NN)+TY(MM)*yTDY(NN)+yTDZ(MM)*TZ(NN)+TZ(MM)*yTDZ(NN)))
                     DRZ(ISP)=DGZ(16)*TX(KK)*TX(LL)*TX(MM)*TX(NN)+qm_qm_e_repul(16)*(zTDX(KK)*TX(LL)*TX(MM)     &
                             *TX(NN)+TX(KK)*zTDX(LL)*TX(MM)*TX(NN)+TX(KK)*TX(LL)*zTDX(MM)*TX(NN)+TX(KK)*TX(LL) &
                             *TX(MM)*zTDX(NN))+DGZ(17)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*TX(MM)*TX(NN)             &
                             +qm_qm_e_repul(17)*((zTDY(KK)*TY(LL)+TY(KK)*zTDY(LL)+zTDZ(KK)*TZ(LL)+TZ(KK)        &
                             *zTDZ(LL))*TX(MM)*TX(NN)+(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))*(zTDX(MM)*TX(NN)+TX(MM)    &
                             *zTDX(NN)))+DGZ(18)*TX(KK)*TX(LL)*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+qm_qm_e_repul(18) &
                             *((zTDX(KK)*TX(LL)+TX(KK)*zTDX(LL))*(TY(MM)*TY(NN)+TZ(MM)*TZ(NN))+TX(KK)*TX(LL)   &
                             *(zTDY(MM)*TY(NN)+TY(MM)*zTDY(NN)+zTDZ(MM)*TZ(NN)+TZ(MM)*zTDZ(NN)))
                     DRX(ISP)=DRX(ISP)+DGX(19)*(TY(KK)*TY(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*TZ(MM)*TZ(NN))        &
                             +qm_qm_e_repul(19)*(xTDY(KK)*TY(LL)*TY(MM)*TY(NN)+TY(KK)*xTDY(LL)*TY(MM)*TY(NN)   &
                             +TY(KK)*TY(LL)*xTDY(MM)*TY(NN)+TY(KK)*TY(LL)*TY(MM)*xTDY(NN)+xTDZ(KK)*TZ(LL)*TZ(MM)&
                             *TZ(NN)+TZ(KK)*xTDZ(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*xTDZ(MM)*TZ(NN)+TZ(KK)*TZ(LL) &
                             *TZ(MM)*xTDZ(NN))+DGX(20)*(TX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)    &
                             *(TY(LL)*TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))    &
                             +TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM))))
                     DRY(ISP)=DRY(ISP)+DGY(19)*(TY(KK)*TY(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*TZ(MM)*TZ(NN))        &
                             +qm_qm_e_repul(19)*(yTDY(KK)*TY(LL)*TY(MM)*TY(NN)+TY(KK)*yTDY(LL)*TY(MM)*TY(NN)   &
                             +TY(KK)*TY(LL)*yTDY(MM)*TY(NN)+TY(KK)*TY(LL)*TY(MM)*yTDY(NN)+yTDZ(KK)*TZ(LL)*TZ(MM)&
                             *TZ(NN)+TZ(KK)*yTDZ(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*yTDZ(MM)*TZ(NN)+TZ(KK)*TZ(LL) &
                             *TZ(MM)*yTDZ(NN))+DGY(20)*(TX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)    &
                             *(TY(LL)*TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))    &
                             +TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM))))
                     DRZ(ISP)=DRZ(ISP)+DGZ(19)*(TY(KK)*TY(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*TZ(MM)*TZ(NN))        &
                             +qm_qm_e_repul(19)*(zTDY(KK)*TY(LL)*TY(MM)*TY(NN)+TY(KK)*zTDY(LL)*TY(MM)*TY(NN)   &
                             +TY(KK)*TY(LL)*zTDY(MM)*TY(NN)+TY(KK)*TY(LL)*TY(MM)*zTDY(NN)+zTDZ(KK)*TZ(LL)*TZ(MM)&
                             *TZ(NN)+TZ(KK)*zTDZ(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*zTDZ(MM)*TZ(NN)+TZ(KK)*TZ(LL) &
                             *TZ(MM)*zTDZ(NN))+DGZ(20)*(TX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)    &
                             *(TY(LL)*TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))    &
                             +TX(NN)*(TY(KK)*TY(MM)+TZ(KK)*TZ(MM))))
!      TO AVOID COMPILER DIFFICULTIES THIS IS DIVIDED
                     xTEMP1=  xTDX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)*(TY(LL)*TY(MM)+TZ(LL)      &
                             *TZ(MM)))+xTDX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(NN)*(TY(KK)*TY(MM)   &
                             +TZ(KK)*TZ(MM)))+TX(KK)*(xTDX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+xTDX(NN)*(TY(LL)  &
                             *TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(xTDX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+xTDX(NN)   &
                             *(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))
                     yTEMP1=  yTDX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)*(TY(LL)*TY(MM)+TZ(LL)      &
                             *TZ(MM)))+yTDX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(NN)*(TY(KK)*TY(MM)   &
                             +TZ(KK)*TZ(MM)))+TX(KK)*(yTDX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+yTDX(NN)*(TY(LL)  &
                             *TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(yTDX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+yTDX(NN)   &
                             *(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))
                     zTEMP1=  zTDX(KK)*(TX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+TX(NN)*(TY(LL)*TY(MM)+TZ(LL)      &
                             *TZ(MM)))+zTDX(LL)*(TX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+TX(NN)*(TY(KK)*TY(MM)   &
                             +TZ(KK)*TZ(MM)))+TX(KK)*(zTDX(MM)*(TY(LL)*TY(NN)+TZ(LL)*TZ(NN))+zTDX(NN)*(TY(LL)  &
                             *TY(MM)+TZ(LL)*TZ(MM)))+TX(LL)*(zTDX(MM)*(TY(KK)*TY(NN)+TZ(KK)*TZ(NN))+zTDX(NN)   &
                             *(TY(KK)*TY(MM)+TZ(KK)*TZ(MM)))

                     xTEMP2=  TX(KK)*(TX(MM)*(xTDY(LL)*TY(NN)+TY(LL)*xTDY(NN)+xTDZ(LL)*TZ(NN)+TZ(LL)*xTDZ(NN))    &
                             +TX(NN)*(xTDY(LL)*TY(MM)+TY(LL)*xTDY(MM)+xTDZ(LL)*TZ(MM)+TZ(LL)*xTDZ(MM)))+TX(LL)   &
                             *(TX(MM)*(xTDY(KK)*TY(NN)+TY(KK)*xTDY(NN)+xTDZ(KK)*TZ(NN)+TZ(KK)*xTDZ(NN))+TX(NN)   &
                             *(xTDY(KK)*TY(MM)+TY(KK)*xTDY(MM)+xTDZ(KK)*TZ(MM)+TZ(KK)*xTDZ(MM)))
                     yTEMP2=  TX(KK)*(TX(MM)*(yTDY(LL)*TY(NN)+TY(LL)*yTDY(NN)+yTDZ(LL)*TZ(NN)+TZ(LL)*yTDZ(NN))    &
                             +TX(NN)*(yTDY(LL)*TY(MM)+TY(LL)*yTDY(MM)+yTDZ(LL)*TZ(MM)+TZ(LL)*yTDZ(MM)))+TX(LL)   &
                             *(TX(MM)*(yTDY(KK)*TY(NN)+TY(KK)*yTDY(NN)+yTDZ(KK)*TZ(NN)+TZ(KK)*yTDZ(NN))+TX(NN)   &
                             *(yTDY(KK)*TY(MM)+TY(KK)*yTDY(MM)+yTDZ(KK)*TZ(MM)+TZ(KK)*yTDZ(MM)))
                     zTEMP2=  TX(KK)*(TX(MM)*(zTDY(LL)*TY(NN)+TY(LL)*zTDY(NN)+zTDZ(LL)*TZ(NN)+TZ(LL)*zTDZ(NN))    &
                             +TX(NN)*(zTDY(LL)*TY(MM)+TY(LL)*zTDY(MM)+zTDZ(LL)*TZ(MM)+TZ(LL)*zTDZ(MM)))+TX(LL)   &
                             *(TX(MM)*(zTDY(KK)*TY(NN)+TY(KK)*zTDY(NN)+zTDZ(KK)*TZ(NN)+TZ(KK)*zTDZ(NN))+TX(NN)   &
                             *(zTDY(KK)*TY(MM)+TY(KK)*zTDY(MM)+zTDZ(KK)*TZ(MM)+TZ(KK)*zTDZ(MM)))

                     DRX(ISP)=DRX(ISP)+qm_qm_e_repul(20)*(xTEMP1+xTEMP2)
                     DRY(ISP)=DRY(ISP)+qm_qm_e_repul(20)*(yTEMP1+yTEMP2)
                     DRZ(ISP)=DRZ(ISP)+qm_qm_e_repul(20)*(zTEMP1+zTEMP2)
                     DRX(ISP)=DRX(ISP)+DGX(21)*(TY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*TY(MM)*TY(NN))        &
                             +qm_qm_e_repul(21)*(xTDY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TY(KK)*xTDY(LL)*TZ(MM)*TZ(NN)   &
                             +TY(KK)*TY(LL)*xTDZ(MM)*TZ(NN)+TY(KK)*TY(LL)*TZ(MM)*xTDZ(NN)+xTDZ(KK)*TZ(LL)*TY(MM)&
                             *TY(NN)+TZ(KK)*xTDZ(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*xTDY(MM)*TY(NN)+TZ(KK)*TZ(LL) &
                             *TY(MM)*xTDY(NN))
                     DRY(ISP)=DRY(ISP)+DGY(21)*(TY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*TY(MM)*TY(NN))        &
                             +qm_qm_e_repul(21)*(yTDY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TY(KK)*yTDY(LL)*TZ(MM)*TZ(NN)   &
                             +TY(KK)*TY(LL)*yTDZ(MM)*TZ(NN)+TY(KK)*TY(LL)*TZ(MM)*yTDZ(NN)+yTDZ(KK)*TZ(LL)*TY(MM)&
                             *TY(NN)+TZ(KK)*yTDZ(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*yTDY(MM)*TY(NN)+TZ(KK)*TZ(LL) &
                             *TY(MM)*yTDY(NN))
                     DRZ(ISP)=DRZ(ISP)+DGZ(21)*(TY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TZ(KK)*TZ(LL)*TY(MM)*TY(NN))        &
                             +qm_qm_e_repul(21)*(zTDY(KK)*TY(LL)*TZ(MM)*TZ(NN)+TY(KK)*zTDY(LL)*TZ(MM)*TZ(NN)   &
                             +TY(KK)*TY(LL)*zTDZ(MM)*TZ(NN)+TY(KK)*TY(LL)*TZ(MM)*zTDZ(NN)+zTDZ(KK)*TZ(LL)*TY(MM)&
                             *TY(NN)+TZ(KK)*zTDZ(LL)*TY(MM)*TY(NN)+TZ(KK)*TZ(LL)*zTDY(MM)*TY(NN)+TZ(KK)*TZ(LL) &
                             *TY(MM)*zTDY(NN))

                     DRX(ISP)=DRX(ISP)+DGX(22)*(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))      &
                             +qm_qm_e_repul(22)*((xTDY(KK)*TZ(LL)+TY(KK)*xTDZ(LL)+xTDZ(KK)*TY(LL)+TZ(KK)        &
                             *xTDY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))+(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(xTDY(MM)  &
                             *TZ(NN)+TY(MM)*xTDZ(NN)+xTDZ(MM)*TY(NN)+TZ(MM)*xTDY(NN)))
                     DRY(ISP)=DRY(ISP)+DGY(22)*(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))      &
                             +qm_qm_e_repul(22)*((yTDY(KK)*TZ(LL)+TY(KK)*yTDZ(LL)+yTDZ(KK)*TY(LL)+TZ(KK)        &
                             *yTDY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))+(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(yTDY(MM)  &
                             *TZ(NN)+TY(MM)*yTDZ(NN)+yTDZ(MM)*TY(NN)+TZ(MM)*yTDY(NN)))
                     DRZ(ISP)=DRZ(ISP)+DGZ(22)*(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))      &
                             +qm_qm_e_repul(22)*((zTDY(KK)*TZ(LL)+TY(KK)*zTDZ(LL)+zTDZ(KK)*TY(LL)+TZ(KK)        &
                             *zTDY(LL))*(TY(MM)*TZ(NN)+TZ(MM)*TY(NN))+(TY(KK)*TZ(LL)+TZ(KK)*TY(LL))*(zTDY(MM)  &
                             *TZ(NN)+TY(MM)*zTDZ(NN)+zTDZ(MM)*TY(NN)+TZ(MM)*zTDY(NN)))
                  endif
               end do !N=M,n_atomic_orbj
            end do !M=1,n_atomic_orbj
         end do !L=K,total_atomic_orb
      end do !K=qm_atomi_orb_start,total_atomic_orb
   else
   ! Just have 1 atomic orbital interacting with 1 atomic orbital so only need
   ! to do (SS/SS)
      DRX(1)=DGX(1)
      DRY(1)=DGY(1)
      DRZ(1)=DGZ(1)
   end if !if(n_atomic_orbi>1.OR.n_atomic_orbj>1)

   ! ---------------------------------------
   ! Gradient of the core repulsion function
   ! ---------------------------------------
   xyzij(1) = vec_qm_qm1
   xyzij(2) = vec_qm_qm2
   xyzij(3) = vec_qm_qm3
   gij = qm_qm_e_repul(1)
   dgij(1) = dgx(1)
   dgij(2) = dgy(1)
   dgij(3) = dgz(1)
   dxyz = zero
   call qm2_core_core_repulsion_dxyz(qmmm_struct, iqm, jqm, rij, onerij, xyzij, gij, dgij, dxyz)
   FNUCX = FNUCX + dxyz(1)
   FNUCY = FNUCY + dxyz(2)
   FNUCZ = FNUCZ + dxyz(3)

!  FIRST, CORE-ELECTRON ATTRACTION DERIVATIVES (MNDO, AM1 and PM3)             
!  ATOM CORE I AFFECTING A.O.S ON J                                     
   ISP=0                                                               
   do M=1,n_atomic_orbj                                                      
      BB=one                                                        
      do N=M,n_atomic_orbj                                                    
         MN=M+qm2_params%pascal_tri1(N)
         ISP=ISP+1      
         temp_real = BB*corei*PSUM(MN)                                               
         FABx=FABx-temp_real*DRX(ISP)                    
         FABy=FABy-temp_real*DRY(ISP)                    
         FABz=FABz-temp_real*DRZ(ISP)                    
         BB=two                                                       
      end do
   end do
!  ATOM CORE J AFFECTING A.O.S ON I                                     
   K=n_atomic_orbj
   K=qm2_params%pascal_tri2(n_atomic_orbj)
   ISP=-K+1                                                            
   do M=qm_atomi_orb_start,total_atomic_orb                                                      
      BB=one
      do N=M,total_atomic_orb                                                    
         MN=M+qm2_params%pascal_tri1(N)
         ISP=ISP+K              
         temp_real = BB*corej*PSUM(MN)                                       
         FABx=FABx-temp_real*DRX(ISP)                    
         FABy=FABy-temp_real*DRY(ISP)                    
         FABz=FABz-temp_real*DRZ(ISP)                    
         BB=two
      end do
   end do

!   NOW FOR COULOMB AND EXCHANGE TERMS (MNDO, AM1 and PM3)                           
   ISP=0                                                               
   do K=qm_atomi_orb_start,total_atomic_orb                                                      
      AA=one                                                     
      KK=qm2_params%pascal_tri1(K)
      do L=K,total_atomic_orb                                                    
         LL=qm2_params%pascal_tri1(L)
         KL=K+LL                                                 
         do M=1,n_atomic_orbj                                                
            BB=one                                                   
            MK=M+KK                                                 
            ML=M+LL                                                 
            do N=M,n_atomic_orbj                                              
               ISP=ISP+1                                               
               MN=M+qm2_params%pascal_tri1(N)
!    COULOMB TERM
               temp_real=AA*BB*PSUM(KL)*PSUM(MN)                                                              
               FAAx=FAAx+temp_real*DRX(ISP)           
               FAAy=FAAy+temp_real*DRY(ISP)           
               FAAz=FAAz+temp_real*DRZ(ISP)           
!  EXCHANGE TERM                                                              
               NK=N+KK                                                 
               NL=N+LL 
               temp_real = AA*BB*0.25D0*(PSUM(MK)*PSUM(NL)+PSUM(NK)*PSUM(ML))                                 
               FAAx = FAAx-temp_real*DRX(ISP) 
               FAAy = FAAy-temp_real*DRY(ISP) 
               FAAz = FAAz-temp_real*DRZ(ISP) 
               BB=two
            end do !N=M,n_atomic_orbj
         end do !M=1,n_atomic_orbj
         AA=two
      end do ! L=K,total_atomic_orb
   end do ! K=qm_atomi_orb_start,total_atomic_orb

   pair_force(1)=FAAx+FABx+FNUCX
   pair_force(1) = -pair_force(1)*EV_TO_KCAL                                 
   pair_force(2)=FAAy+FABy+FNUCY                                           
   pair_force(2) = -pair_force(2)*EV_TO_KCAL                                 
   pair_force(3)=FAAz+FABz+FNUCZ
   pair_force(3) = -pair_force(3)*EV_TO_KCAL                                 

end subroutine qm2_deriv_qm_analyt

subroutine qm2_dhc(qm2_struct, qmmm_struct, P,iqm, jqm,qmitype,qmjtype,xyz_qmi,xyz_qmj,natqmi, &
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
      use qmmm_module        , only: qm2_params, OVERLAP_CUTOFF, qmmm_nml, qm2_structure
      use Rotation           , only: GetRotationMatrix, Rotate2Center2Electron, RotateCore   
      use qm2_fock_d         , only: W2Fock_atompair
      use qmmm_struct_module, only : qmmm_struct_type
 
      implicit none

!Passed in
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
             call qm2_h1elec(r2InAu,xyz_qmi(1),                                  &
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
        
            call qm2_h1elec_d(qm2_struct,r2InAu,xyz_qmi(1:3), xyz_qmj(1:3),  &
                        n_atomic_orbi,n_atomic_orbj,                &
                        firstIndexAO_i, firstIndexAO_j, qmitype, qmjtype,  &
                        n_atomic_orbi+n_atomic_orbj, H)                

        end if ! ((n_atomic_orbi.lt.9) .and. (n_atomic_orbj.lt.9))
      end if !(R2 < OVERLAP_CUTOFF)

      KR=1
      hasDOrbital=((n_atomic_orbi.ge.9) .or. (n_atomic_orbj.ge.9))
      call GetRotationMatrix(xyz_qmj-xyz_qmi, rotationMatrix, hasDOrbital)        
      call qm2_rotate_qmqm(qm2_struct,qmmm_struct,-1,iqm,jqm,natqmi,natqmj,xyz_qmi,xyz_qmj,            &
                  W(KR),KR, RI, core)

      if (hasDOrbital) then   ! spd case

        i_dimension=n_atomic_orbi*(n_atomic_orbi+1)/2
        j_dimension=n_atomic_orbj*(n_atomic_orbj+1)/2
       
        allocate(ww(1:j_dimension, 1:i_dimension)) 
        WW=0.0D0

        ! calculate the 2-center integrals and core-core interaction integrals
        call qm2_repp_d(qmmm_struct, qmitype,qmjtype,rijInAu,RI,CORE,WW,i_dimension,j_dimension,1)
 
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
   
   call qm2_core_core_repulsion(qmmm_struct, iqm, jqm, rij, oneOverRij, RI, enuclr)         
        
    ! put what we have now to the Fock matrix
    F(1:linear)=H(1:linear)
       
    ! 2-center 2-electron contribution to the Fock matrix      
    call W2Fock_atompair(W, F, P, n_atomic_orbj, n_atomic_orbi,  &
      firstIndexAO_j, firstIndexAO_i)
    call W2Fock_atompair(W, F, P, n_atomic_orbi, n_atomic_orbj,  &
      firstIndexAO_i, firstIndexAO_j)   

    EE=qm2_HELECT(n_atomic_orbi+n_atomic_orbj-1,P,H,F)   
    DENER=EE+ENUCLR
   
   
end subroutine qm2_dhc




