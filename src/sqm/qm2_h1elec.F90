! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
subroutine qm2_h1elec(qm2_params,R2,XI,XJ, n_atomic_orbi,n_atomic_orbj,SHMAT, qmitype,qmjtype)

!***********************************************************************
!
!  qm2_h1elec forms the one-electron matrix between two atoms and
!  calculates the overlaps.
!
!  Current code optimised by Ross Walker (TSRI, 2005)
!
!   ON INPUT
!               XI   = COORDINATES OF FIRST ATOM.
!               XJ   = COORDINATES OF SECOND ATOM.
!  n_atomic_orbi,j   = number of atomic orbitals on i and j.
!                                                 
!   ON OUTPUT   SHMAT = MATRIX OF ONE-ELECTRON INTERACTIONS.
!                                                           
!***********************************************************************
      use constants          , only : A2_TO_BOHRS2, A_TO_BOHRS, half
      use ElementOrbitalIndex, only : MaxValenceOrbitals,MaxGaussianExpansion
      use SlaterOverlap      , only : GetSlaterOverlap
      use qmmm_module        , only : EXPONENTIAL_CUTOFF
      use qm2_params_module,  only : qm2_params_type

      implicit none
      type(qm2_params_type), intent(inout) :: qm2_params
      _REAL_ , intent(in) :: R2,XI(3),XJ(3)
      integer, intent(in) :: n_atomic_orbi, n_atomic_orbj
      integer, intent(in) :: qmitype,qmjtype      
      _REAL_, intent(out) :: SHMAT(MaxValenceOrbitals, MaxValenceOrbitals)
      
      ! Pointers to avoid having to do a look up on a 4 dimensional array within a loop.      
      _REAL_, pointer :: sxs_over_sas(:,:), ss_eqn(:,:)
      _REAL_, pointer :: sxp_over_sap(:,:), pxs_over_pas(:,:), sp_ovlp(:,:), ps_ovlp(:,:)
      _REAL_, pointer :: sxd_over_sad(:,:), dxs_over_das(:,:), sd_ovlp(:,:), ds_ovlp(:,:)
      _REAL_, pointer :: pxd_over_pad(:,:), dxp_over_dap(:,:), pd_ovlp(:,:), dp_ovlp(:,:)
          
      _REAL_, pointer :: pxp_over_pap(:,:), pp_ovlp_ieqj1(:,:), pp_ovlp_ieqj2(:,:), pp_ovlp_inj(:,:)    
      _REAL_, pointer :: dxd_over_pap(:,:), dd_ovlp_ieqj1(:,:), dd_ovlp_ieqj2(:,:), dd_ovlp_inj(:,:) 
              
      _REAL_          :: betasas,betasap,betapas,betasad, betadas, betapap, betapad, betadap, betadad

      _REAL_ :: bij_temp
      _REAL_ ::SH, vec_qm_qm(3)
      _REAL_ ::ADBR2, TOMB
      _REAL_ ::rab
      _REAL_ ::zeta_si, zeta_sj, zeta_pi, zeta_pj, zeta_di, zeta_dj
      _REAL_ :: zeta_a_old=-1.0D99, zeta_b_old=-1.0D99, rab_old=-1.0D99
      _REAL_ :: local_A(15), local_B(15)
      integer:: ntotal=-1
      integer:: i,j,k,l, ii,jj
      integer:: ng, ni, nj

      logical::isISAtom=.false., isISPAtom=.false.,isISPDAtom=.false.
      logical::isJSAtom=.false., isJSPAtom=.false.,isJSPDAtom=.false.
      
      !  Assoicate the pointers

      ng = MaxGaussianExpansion
        
      sxs_over_sas => qm2_params%atom_orb_zz_sxs_over_sas(1:ng, 1:ng,qmitype,qmjtype)
      ss_eqn       => qm2_params%atom_orb_ss_eqn(1:ng, 1:ng,qmitype,qmjtype)
        
      sxp_over_sap => qm2_params%atom_orb_zz_sxp_over_sap(1:ng, 1:ng,qmitype,qmjtype)
      pxs_over_pas => qm2_params%atom_orb_zz_sxp_over_sap(1:ng, 1:ng,qmjtype,qmitype)
      sp_ovlp      => qm2_params%atom_orb_sp_ovlp(1:ng, 1:ng,qmitype,qmjtype)
      ps_ovlp      => qm2_params%atom_orb_sp_ovlp(1:ng, 1:ng,qmjtype,qmitype)
 
      sxd_over_sad => qm2_params%atom_orb_zz_sxd_over_sad(1:ng, 1:ng,qmitype,qmjtype)
      dxs_over_das => qm2_params%atom_orb_zz_sxd_over_sad(1:ng, 1:ng,qmjtype,qmitype)
      sd_ovlp      => qm2_params%atom_orb_sd_ovlp(1:ng, 1:ng,qmitype,qmjtype)
      ds_ovlp      => qm2_params%atom_orb_sd_ovlp(1:ng, 1:ng,qmjtype,qmitype)
        
      pxp_over_pap => qm2_params%atom_orb_zz_pxp_over_pap(1:ng, 1:ng,qmitype,qmjtype)
      pp_ovlp_ieqj1=> qm2_params%atom_orb_pp_ovlp_ieqj1(1:ng, 1:ng,qmitype,qmjtype)
      pp_ovlp_ieqj2=> qm2_params%atom_orb_pp_ovlp_ieqj2(1:ng, 1:ng,qmitype,qmjtype)
      pp_ovlp_inj  => qm2_params%atom_orb_pp_ovlp_inj(1:ng, 1:ng,qmitype,qmjtype)

      pxd_over_pad => qm2_params%atom_orb_zz_pxd_over_pad(1:ng, 1:ng,qmitype,qmjtype)
      dxp_over_dap => qm2_params%atom_orb_zz_pxd_over_pad(1:ng, 1:ng,qmjtype,qmitype)
      pd_ovlp      => qm2_params%atom_orb_pd_ovlp(1:ng, 1:ng,qmitype,qmjtype)
      dp_ovlp      => qm2_params%atom_orb_pd_ovlp(1:ng, 1:ng,qmjtype,qmitype)
 
      dxd_over_pap => qm2_params%atom_orb_zz_dxd_over_dad(1:ng, 1:ng,qmitype,qmjtype)
      dd_ovlp_ieqj1=> qm2_params%atom_orb_dd_ovlp_ieqj1(1:ng, 1:ng,qmitype,qmjtype)
      dd_ovlp_ieqj2=> qm2_params%atom_orb_dd_ovlp_ieqj2(1:ng, 1:ng,qmitype,qmjtype)
      dd_ovlp_inj  => qm2_params%atom_orb_dd_ovlp_inj(1:ng, 1:ng,qmitype,qmjtype)


        
      betasas = qm2_params%betasas(qmitype,qmjtype)
      betasap = qm2_params%betasap(qmitype,qmjtype)
      betapas = qm2_params%betasap(qmjtype,qmitype)
      betasad = qm2_params%betasad(qmitype,qmjtype)
      betadas = qm2_params%betasad(qmjtype,qmitype)        
      betapap = qm2_params%betapap(qmitype,qmjtype)
      betapad = qm2_params%betapad(qmitype,qmjtype)
      betadap = qm2_params%betapad(qmjtype,qmitype)        
      betadad = qm2_params%betadad(qmitype,qmjtype) 
        
      ni = qm2_params%sp_quantum_number(qmitype)
      nj = qm2_params%sp_quantum_number(qmjtype) 
      zeta_si = qm2_params%s_orb_exp_by_type(qmitype)
      zeta_sj = qm2_params%s_orb_exp_by_type(qmjtype)  
      zeta_pi = qm2_params%p_orb_exp_by_type(qmitype)
      zeta_pj = qm2_params%p_orb_exp_by_type(qmjtype) 
      zeta_di = qm2_params%d_orb_exp_by_type(qmitype)
      zeta_dj = qm2_params%d_orb_exp_by_type(qmjtype)
      rab = sqrt(R2)            
             
             
 ! logicals
                      
      isISAtom   = (n_atomic_orbi==1)
      isISPAtom  = (n_atomic_orbi==4)
      isISPDAtom = (n_atomic_orbi==9)   
      
      isJSAtom   = (n_atomic_orbj==1)
      isJSPAtom  = (n_atomic_orbj==4)
      isJSPDAtom = (n_atomic_orbj==9)         

!***********************************************************************        
!   CALCULATE THE OVERLAP INTEGRALS USING A GAUSSIAN EXPANSION         *        
!         STO-6G BY R.F. STEWART, J. CHEM. PHYS., 52 431-438, 1970     *        
!                                                                      *
!         FILL SHMAT=  4X4 ARRAY OF OVERLAPS, IN ORDER S,PX,PY,PZ      *        
!***********************************************************************        

!------------------------------
!Current Code and optimisation:
!       Ross Walker (TSRI, 2004)
!------------------------------

!     R2   =  INTERATOMIC DISTANCE^2 IN BOHRS2

      vec_qm_qm(1) = (XI(1)-XJ(1))*A_TO_BOHRS
      vec_qm_qm(2) = (XI(2)-XJ(2))*A_TO_BOHRS
      vec_qm_qm(3) = (XI(3)-XJ(3))*A_TO_BOHRS

!All atoms have S-orbitals
!    S-S
      do K=1,6 !1 to NGAUSS
         do L=1,6
           ADBR2=sxs_over_sas(k,l)*R2
!          CHECK OF OVERLAP IS NON-ZERO BEFORE DOING EXPONENTIAL
           if(ADBR2 < EXPONENTIAL_CUTOFF) then
              SHMAT(1,1)=SHMAT(1,1)+ss_eqn(k,l)*EXP(-ADBR2)
           endif
         end do
      end do
      !Multiply by S-S beta factor
      SHMAT(1,1)=SHMAT(1,1)*half*betasas

      !I:S and J:P
      if (.not. isJSAtom) then
         do K=1,6 !1 to NGAUSS
            do L=1,6
              ADBR2=sxp_over_sap(k,l)*R2
              ! CHECK OF OVERLAP IS NON-ZERO BEFORE DOING EXPONENTIAL
              if(ADBR2 < EXPONENTIAL_CUTOFF) then
                SH=sp_ovlp(k,l)*EXP(-ADBR2)
                SHMAT(1,2)=SHMAT(1,2)+SH*vec_qm_qm(1)
                SHMAT(1,3)=SHMAT(1,3)+SH*vec_qm_qm(2)
                SHMAT(1,4)=SHMAT(1,4)+SH*vec_qm_qm(3)
              endif
            end do
         end do
         !Multiply by S-P beta factor
         bij_temp=GetSlaterOverlap(ni, 0, nj, 1, 0, zeta_si, zeta_a_old, zeta_pj, zeta_b_old, rab, rab_old,ntotal,local_A,local_B)
         bij_temp=GetSlaterOverlap(ni, 1, nj, 0, 0, zeta_pi, zeta_a_old, zeta_sj, zeta_b_old, rab, rab_old,ntotal,local_A,local_B)

         bij_temp=half*betasap
         SHMAT(1,2)=SHMAT(1,2)*bij_temp
         SHMAT(1,3)=SHMAT(1,3)*bij_temp
         SHMAT(1,4)=SHMAT(1,4)*bij_temp
      end if

      !I:S and J:D
      if (isJSPDAtom) then
         do K=1,6 !1 to NGAUSS
            do L=1,6
              ADBR2=sxd_over_sad(k,l)*R2
              ! CHECK OF OVERLAP IS NON-ZERO BEFORE DOING EXPONENTIAL
              if(ADBR2 < EXPONENTIAL_CUTOFF) then
                SH=sd_ovlp(k,l)*EXP(-ADBR2)
                SHMAT(1,5)=SHMAT(1,5)+SH*vec_qm_qm(1)
                SHMAT(1,6)=SHMAT(1,6)+SH*vec_qm_qm(2)
                SHMAT(1,7)=SHMAT(1,7)+SH*vec_qm_qm(3)
                SHMAT(1,8)=SHMAT(1,8)+SH*vec_qm_qm(3)
                SHMAT(1,9)=SHMAT(1,9)+SH*vec_qm_qm(3)
              endif
            end do
         end do
         !Multiply by S-D beta factor
         bij_temp=half*betasad
         SHMAT(1,5)=SHMAT(1,5)*bij_temp
         SHMAT(1,6)=SHMAT(1,6)*bij_temp
         SHMAT(1,7)=SHMAT(1,7)*bij_temp
         SHMAT(1,8)=SHMAT(1,8)*bij_temp
         SHMAT(1,9)=SHMAT(1,9)*bij_temp                  
      end if

      !I:P and J:S
      if (.not. isISAtom) then
         do K=1,6 !1 to NGAUSS
            do L=1,6
              ADBR2=pxs_over_pas(l,k)*R2
              ! CHECK OF OVERLAP IS NON-ZERO BEFORE DOING EXPONENTIAL
              if(ADBR2 < EXPONENTIAL_CUTOFF) then
                SH=-ps_ovlp(l,k)*EXP(-ADBR2)
                SHMAT(2,1)=SHMAT(2,1)+SH*vec_qm_qm(1)
                SHMAT(3,1)=SHMAT(3,1)+SH*vec_qm_qm(2)
                SHMAT(4,1)=SHMAT(4,1)+SH*vec_qm_qm(3)
              endif
            end do
         end do
         !Multiply by P-S beta factor
         bij_temp=half*betapas
         SHMAT(2,1)=SHMAT(2,1)*bij_temp
         SHMAT(3,1)=SHMAT(3,1)*bij_temp
         SHMAT(4,1)=SHMAT(4,1)*bij_temp
      end if

      !I:P and J:P
      if ( (.not. isISAtom) .AND. (.not. isJSAtom)) then
        bij_temp=half*betapap
        !  the "3" here is for (x,y,z} of p-orbitals
        do I=1,3
           ii=i+1
           do J=1,3
           jj=j+1
           ! P-P
            TOMB=vec_qm_qm(i)*vec_qm_qm(j)
            do K=1,6 !1 to NGAUSS
               do L=1,6
                 ADBR2=pxp_over_pap(k,l)*R2
                 ! CHECK OF OVERLAP IS NON-ZERO BEFORE DOING THE EXPONENTIAL
                 if(ADBR2 < EXPONENTIAL_CUTOFF) then
                   if(ii.EQ.jj) then
                     SH=EXP(-ADBR2)*(pp_ovlp_ieqj1(k,l)*TOMB+pp_ovlp_ieqj2(k,l))
                   else
                     SH=EXP(-ADBR2)*TOMB*pp_ovlp_inj(k,l)
                   end if
                   SHMAT(ii,jj)=SHMAT(ii,jj)+SH
                 endif
               end do !I=1,6
            end do !K=1,6
            !Multiply by P-P beta factor
            SHMAT(ii,jj)=SHMAT(ii,jj)*bij_temp
          end do !j=1,3: x,y,z for j
        end do !i=1,3: x,y,z for i
      end if
 
      return
end subroutine qm2_h1elec


