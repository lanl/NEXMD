! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
subroutine qm2_get_qmmm_forces(qmmm_struct, dxyzqm,qm_xcrd,dxyzmm,scf_mchg)

! This routine calculates force on the MM atoms due to the QM atoms
! and vice versa. The forces are added to dxyzmm and dxyzqm respectively.

! Currently only analytical derivatives are supported (for S and P atoms).
!
! only numerical derivatives are available for D atoms
! TL (Rutgers, 2011)
!
!     Variable Definitions:
!
!     qm_coords - Cartesian coordinates of QM atoms.
!     dxyzqm - Coupled potential energy derivatives with respect to
!              movement of QM atoms.
!     qm_xcrd - Cartesian coordinates of QM and MM atoms. In same order as cutoff list.
!               also contains scaled charge in 4th element.
!     dxyzmm - Coupled potential energy derivatives with respect to
!              movement of MM atoms.
!     Note this routine does not need iqmatoms since it gets mm
!     atoms from the pair list and puts QM forces in it's own
!     array which gets moved to the main force array later.

! Current Code by Ross Walker (TSRI, 2004)
! Derivative optimisations by Ross Walker and Mike Crowley (TSRI, 2004)

      use constants  , only : zero, one, A_TO_BOHRS, A2_TO_BOHRS2, EV_TO_KCAL, &
                              AU_TO_EV, BOHRS_TO_A
      use Rotation   , only : GetRotationMatrix, RotateCore         
      use qmmm_module, only : qmmm_nml, qm2_struct, qm2_params, &
        qmmm_mpi, alph_MM, qmmm_opnq
      use opnq, only: Opnq_fock_atom_pair, Opnq_LJ_atom_pair  
      use qmmm_struct_module, only : qmmm_struct_type

      implicit none

!Passed in
      type(qmmm_struct_type), intent(inout) :: qmmm_struct
      _REAL_ , intent(out) :: dxyzqm(3,qmmm_struct%nquant_nlink)
      _REAL_ , intent(in) :: qm_xcrd(4,qmmm_struct%qm_mm_pairs)
      _REAL_ , intent(out) :: dxyzmm(*)
      _REAL_ , intent(in) :: scf_mchg(*)

!Local
      _REAL_ psum(45), psum_light
      _REAL_ pair_force(3)
      _REAL_ qm_atom_coord(3), fqm(3)
      _REAL_ qm_atom_core, qm_atom_alpa
       integer jj, jf, jl, ij, i, k, ii, j
      integer n_atomic_orb !Number of atomic orbitals on qm atom
      integer loop_count !Keeps track of number of times through nquant * ni_mm loop
      integer inner_loop_count !xyz offset for loop over pairs.  
 
 !qm2_Helect is a function
      _REAL_ :: qm2_HELECT     
         
      _REAL_, parameter :: change=2.0D-6, halfChange=change/2.0D0, oneChange=1.0D0/change
      _REAL_, parameter :: changeScale=1.0d3
      _REAL_ :: temp(3)
      _REAL_ :: w(45), RI(22), core(10,2), rotationMatrix(15,45), H(45) 
      _REAL_ :: r2, r2InAu, rij, rijInAu, oneOverRij, AA, EE, scale
      _REAL_ :: total_density, opnq_pair, fock_opnq_pair, LJ_pair      
      _REAL_ :: qm2_switch_func, qm2_switch_derv  ! functions
      _REAL_ :: f_switch, df_switch, pot_tmp
      _REAL_ :: vasya
      logical :: s_atom, sp_atom, spd_atom
      
      if (qmmm_struct%qm_mm_pairs == 0) return

      loop_count=0

      !do jj=1,qmmm_struct%nquant_nlink
      do jj=qmmm_mpi%nquant_nlink_start,qmmm_mpi%nquant_nlink_end
         jf=qm2_params%orb_loc(1,jj)
         jl=qm2_params%orb_loc(2,jj)
         qm_atom_coord(1:3) = qmmm_struct%qm_coords(1:3,jj)
         qm_atom_core = qm2_params%core_chg(jj)

         qm_atom_alpa = qm2_params%cc_exp_params(jj)
         n_atomic_orb = qm2_params%natomic_orbs(jj)
         s_atom   = (n_atomic_orb==1)         
         sp_atom  = (n_atomic_orb==4)
         spd_atom = (n_atomic_orb==9)
                  
         fqm(1:3)=0.0d0
! Split into heavy and light atoms here - means code duplication but is good for speed.
         if (sp_atom .or. spd_atom) then
!        Get elements of density matrix involving only orbitals centered
!        on current qm atom.
           ij = 0
           do i=jf,jl
              k = qm2_params%pascal_tri1(i) + jf - 1
              do j=jf,i
                ij = ij + 1
                k = k + 1
                psum(ij) = qm2_struct%den_matrix(k)
              end do
           end do
           
         end if
         
         if (spd_atom) then
          inner_loop_count=1
          do j=1,qmmm_struct%qm_mm_pairs
              loop_count=loop_count+1
              do k=1,3
                qm_atom_coord(k)=qm_atom_coord(k)+halfChange
                    !AA               
                    H=0.0D0
                    temp(1:3)=qm_xcrd(1:3,j)-qm_atom_coord(1:3)
                    r2=sum( temp(1:3) **2) 
                    rij=sqrt(r2)
                    rijInAu=rij*A_TO_BOHRS
                    oneOverRij=one/rij
                    r2InAu=r2*A2_TO_BOHRS2
                    call GetRotationMatrix(temp(1:3), rotationMatrix, spd_atom)   
                    call qm2_repp_d(qmmm_struct, qmmm_struct%qm_atom_type(jj),0,rijInAu,RI,CORE,W,45,1,0) 
                    core=core*qm_xcrd(4,j)
                    call RotateCore(1,0, n_atomic_orb,0,1, 0            &
                        ,core,rotationMatrix,H)
                         
                    AA=- core(1,1)*qm_atom_core
                    scale=abs(AA*exp(-qm2_params%cc_exp_params(jj)*RIJ)+EXP(-ALPH_MM*RIJ))
                    AA=qm2_HELECT(8, psum, H, H, vasya)+AA+scale
                    
                qm_atom_coord(k)=qm_atom_coord(k)-change
                    !EE
                    H=0.0D0
                    temp(1:3)=qm_xcrd(1:3,j)-qm_atom_coord(1:3)
                    r2=sum( temp(1:3) **2) 
                    rij=sqrt(r2)
                    rijInAu=rij*A_TO_BOHRS
                    oneOverRij=one/rij
                    r2InAu=r2*A2_TO_BOHRS2                   
                    call GetRotationMatrix(temp(1:3), rotationMatrix, spd_atom)   
                    call qm2_repp_d(qmmm_struct, qmmm_struct%qm_atom_type(jj),0,rijInAu,RI,CORE,W,45,1,0) 
                    core=core*qm_xcrd(4,j)
                    call RotateCore(1,0, n_atomic_orb,0, 1, 0            &
                        ,core,rotationMatrix,H)
     
                    EE=- core(1,1)*qm_atom_core
                    scale=abs(EE*exp(-qm2_params%cc_exp_params(jj)*RIJ)+EXP(-ALPH_MM*RIJ))
                    EE=qm2_HELECT(8, psum, H, H, vasya)+EE+scale
                    
                    
                qm_atom_coord(k)=qm_atom_coord(k)+halfChange
                
                pair_force(k)=(AA-EE)*EV_TO_KCAL*onechange
              
              end do  !k
  
              if ( qmmm_nml%qmmm_switch ) then
                 H=0.0D0
                 temp(1:3)=qm_xcrd(1:3,j)-qm_atom_coord(1:3)
                 r2=sum( temp(1:3) **2) 
                 rij=sqrt(r2)
                 rijInAu=rij*A_TO_BOHRS
                 oneOverRij=one/rij
                 r2InAu=r2*A2_TO_BOHRS2
                 call GetRotationMatrix(temp(1:3), rotationMatrix, spd_atom)   
                 call qm2_repp_d(qmmm_struct, qmmm_struct%qm_atom_type(jj),0,rijInAu,RI,CORE,W,45,1,0) 
                 core=core*qm_xcrd(4,j)
                 call RotateCore(1,0, n_atomic_orb,0,1, 0            &
                     ,core,rotationMatrix,H)
                      
                 AA=- core(1,1)*qm_atom_core
                 scale=abs(AA*exp(-qm2_params%cc_exp_params(jj)*RIJ)+EXP(-ALPH_MM*RIJ))
                 AA=qm2_HELECT(8, psum, H, H,vasya)+AA+scale

                 ! switch function
                 if (rij <= qmmm_nml%r_switch_lo) then
                   f_switch = one
                   df_switch= zero
                 else if (rij < qmmm_nml%r_switch_hi) then  
                   f_switch = qm2_switch_func(RIJ, qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi)
                   df_switch= qm2_switch_derv(RIJ, qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi)
                 else
                   f_switch = zero
                   df_switch= zero
                 end if

                 AA=oneOverRij*df_switch*AA*EV_TO_KCAL
                                           !converts EV to kcal

                 pair_force(1) = pair_force(1)*f_switch - temp(1)*AA 
                 pair_force(2) = pair_force(2)*f_switch - temp(2)*AA 
                 pair_force(3) = pair_force(3)*f_switch - temp(3)*AA 

                 pot_tmp = scf_mchg(jj) * qm_xcrd(4,j) * oneOverRij &
                       & * AU_TO_EV * EV_TO_KCAL * BOHRS_TO_A
                          !converts AU to kcal
                 AA=oneOverRij*df_switch*pot_tmp

                 pair_force(1) = pair_force(1) + temp(1)*AA 
                 pair_force(2) = pair_force(2) + temp(2)*AA 
                 pair_force(3) = pair_force(3) + temp(3)*AA 

                 AA = pot_tmp*(one-f_switch)*oneOverRij*oneOverRij

                 pair_force(1) = pair_force(1) + temp(1)*AA
                 pair_force(2) = pair_force(2) + temp(2)*AA
                 pair_force(3) = pair_force(3) + temp(3)*AA
              end if
  
              dxyzmm(inner_loop_count) = dxyzmm(inner_loop_count) - pair_force(1)
              dxyzmm(inner_loop_count+1) = dxyzmm(inner_loop_count+1) - pair_force(2)
              dxyzmm(inner_loop_count+2) = dxyzmm(inner_loop_count+2) - pair_force(3)
              fqm(1:3) = fqm(1:3) + pair_force(1:3)
              inner_loop_count=inner_loop_count+3

           end do         

         else if (sp_atom) then
!RCW: We could maybe move this outside of the loop for optimisation purposes.
!        Loop over MM atoms that are in the interaction list for current
!        QM atom. We have the loop repeated twice here. Once for when we
!        have the qm_mm 1 elec repulsion integrals in memory and once for
!        when we need to calculate them on the fly.
!        This is excess code but it avoids the if(...) being inside the
!        inner loop.
!We don't have QM-MM 1 e repul in memory so calc on the fly.
           inner_loop_count=1
           do ii=1,qmmm_struct%qm_mm_pairs
              loop_count=loop_count+1
!             Get analytical derivatives with respect to x, y, and z
!             directions for the interaction of the current QM-MM pair.
              call qm2_deriv_qmmm_heavy(qmmm_struct, jj,loop_count, &
                                  psum,qm_atom_coord,qm_xcrd(1,ii),n_atomic_orb, &
                                  pair_force,qm_atom_core,qm_atom_alpa,scf_mchg(jj))
  
              dxyzmm(inner_loop_count) = dxyzmm(inner_loop_count) - pair_force(1)
              dxyzmm(inner_loop_count+1) = dxyzmm(inner_loop_count+1) - pair_force(2)
              dxyzmm(inner_loop_count+2) = dxyzmm(inner_loop_count+2) - pair_force(3)
              fqm(1:3) = fqm(1:3) + pair_force(1:3)
              inner_loop_count=inner_loop_count+3
             
           end do  

        else !If (sp_atom)
          !Light atom - see above for comments
           k = qm2_params%pascal_tri1(jf) + jf 
           psum_light = qm2_struct%den_matrix(k)
           inner_loop_count=1
           do ii=1,qmmm_struct%qm_mm_pairs
              loop_count=loop_count+1
              call qm2_deriv_qmmm_light(qmmm_struct, jj,ii,loop_count,psum_light,qm_atom_coord,qm_xcrd(1,ii), &
                                        pair_force,qm_atom_core,qm_atom_alpa,scf_mchg(jj))
              dxyzmm(inner_loop_count:inner_loop_count+2) = &
                dxyzmm(inner_loop_count:inner_loop_count+2) - pair_force(1:3)
              fqm(1:3) = fqm(1:3) + pair_force(1:3)
              inner_loop_count=inner_loop_count+3
           end do
        end if !if (sp_atom)
        
        ! OPNQ correction
        if (qmmm_opnq%useOPNQ) then
          
            total_density=0.0D0

            do i=qm2_params%orb_loc(1,jj), qm2_params%orb_loc(2,jj)
                total_density=total_density+qm2_struct%den_matrix(qm2_params%pascal_tri2(i))
            end do
            
          inner_loop_count=1
          do j=1,qmmm_struct%qm_mm_pairs
              
              pair_force=0.d0
             
              ! electronic part; use finite difference to calculate the derivatives
              ! the OPNQ subroutines will not work properly for very small difference
              ! hence the step is multiplied by changeScale
!              do k=1,3
!                qmmm_struct%qm_coords(k,jj)=qmmm_struct%qm_coords(k,jj)+halfChange*changeScale
!                call Opnq_fock_atom_pair(jj, j, opnq_pair, fock_opnq_pair)
!                AA=total_density*fock_opnq_pair*.5d0
!                   
!                qmmm_struct%qm_coords(k,jj)=qmmm_struct%qm_coords(k,jj)-change*changeScale
!                call Opnq_fock_atom_pair(jj, j, opnq_pair, fock_opnq_pair)
!                EE=total_density*fock_opnq_pair*.5d0
! 
!                qmmm_struct%qm_coords(k,jj)=qmmm_struct%qm_coords(k,jj)+halfChange*changeScale
!                pair_force(k)=(AA-EE)*EV_TO_KCAL*onechange/changeScale
!              end do  !k
              
              ! opnq 
              pair_force=0.d0
              call Opnq_fock_atom_pair(qmmm_struct,jj, j, opnq_pair, fock_opnq_pair, temp(1), temp(2), temp(3))
              temp=temp*EV_TO_KCAL
              pair_force=pair_force + temp
               
              ! LJ
              call Opnq_LJ_atom_pair(qmmm_struct, jj, j, LJ_pair, temp(1), temp(2), temp(3))
              temp= - temp*EV_TO_KCAL
              pair_force=pair_force + temp       
               
              dxyzmm(inner_loop_count:inner_loop_count+2) = &
                dxyzmm(inner_loop_count:inner_loop_count+2) - pair_force(1:3)
              fqm(1:3) = fqm(1:3) + pair_force(1:3)
              inner_loop_count=inner_loop_count+3
           end do         
        end if ! (qmmm_opnq%useOPNQ) 
        
        
        dxyzqm(1:3,jj) = dxyzqm(1:3,jj) + fqm(1:3) 
      end do !jj=1,qmmm_struct%nquant
      
      return

end subroutine qm2_get_qmmm_forces

subroutine qm2_deriv_qmmm_light(qmmm_struct, iqm,jpair,loop_count,psum_light,xyz_qm,xyz_mm,pair_force, &
                         qm_atom_core,alpa,qm_charge)
!For light atoms
!See heavy version of routine for comments
!  Current Version: Ross Walker (TSRI, 2005)

      use qmmm_module, only : qmmm_nml,qm2_params, qm2_struct, qm2_rij_eqns, &
                              alph_mm, EXPONENTIAL_CUTOFF
      use constants  , only : A_TO_BOHRS, A2_TO_BOHRS2, AU_TO_EV, A2_TO_BOHRS2xAU_TO_EV, &
                              EV_TO_KCAL, AU_TO_KCAL, BOHRS_TO_A, zero, one, two
      use qmmm_struct_module, only : qmmm_struct_type

      implicit none

! ON return, pair_force HOLDS ANALYTICAL DERIVATIVES                                   

!Passed in
      type(qmmm_struct_type), intent(in) :: qmmm_struct
      _REAL_, intent(in) :: xyz_qm(3),xyz_mm(4),psum_light
      _REAL_, intent(out) :: pair_force(3)
      integer, intent(in) :: iqm, loop_count, jpair
      _REAL_, intent(in) ::  qm_atom_core, alpa
      _REAL_, intent(in) :: qm_charge

!Local
      _REAL_ FABX, FABY, FABZ ,FNUCX, FNUCY, FNUCZ
      _REAL_ r2, rij, onerij, rr2
      _REAL_ vec_qm_mm1, vec_qm_mm2, vec_qm_mm3
      _REAL_ c1, mm_charge
      _REAL_ ee, sqrtaee
      _REAL_ DGX, DGY, DGZ
      _REAL_ EXP1, EXP2, EXP3, EXP4
      _REAL_ qmmm_erep, anam1, oner2, temp_real, temp_real2
      _REAL_ corrfx, corrfy, corrfz, temp_exp3, temp_exp4
      _REAL_ sf1, sf2, rho_pm3mmx
      _REAL_ qm2_switch_func, qm2_switch_derv  ! functions
      _REAL_ f_switch, df_switch, ELEC, ENUC, scale, AA, pot_tmp
      integer :: qmitype, i

      vec_qm_mm1=xyz_qm(1) - xyz_mm(1)
      vec_qm_mm2=xyz_qm(2) - xyz_mm(2)
      vec_qm_mm3=xyz_qm(3) - xyz_mm(3)
      mm_charge = xyz_mm(4)
#include "qm2_array_locations.h"

! BEGIN --- THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM
!           Note: This could technically be done at the same time we do the
!                 energy in hcore_qmmm.

    !if ((qmmm_nml%qmmm_int == 3 .or. qmmm_nml%qmmm_int == 4) .and. qmmm_nml%qmtheory%PM3) then  ! PM3/MM*
    if (qmmm_struct%PM3MMX_INTERFACE) then
      ! PM3/MM* MODIFIED QM-MM INTERFACE

      qmitype = qmmm_struct%qm_atom_type(iqm)

      ! Different s factors are applied to H-H and H-heavy
      if (qmmm_struct%qm_mm_pair_atom_numbers(jpair) == 1) then
         sf1 = qm2_params%scale_factor1_pm3mmx(1,qmitype)
         sf2 = qm2_params%scale_factor2_pm3mmx(1,qmitype)
      else
         sf1 = qm2_params%scale_factor1_pm3mmx(2,qmitype)
         sf2 = qm2_params%scale_factor2_pm3mmx(2,qmitype)
      end if
      ! rho_pm3mmx in qmmm_int==3 is always zero
      rho_pm3mmx = qm2_params%rho_pm3mmx(qmitype)

      r2 = vec_qm_mm1*vec_qm_mm1+vec_qm_mm2*vec_qm_mm2+vec_qm_mm3*vec_qm_mm3
      RR2=r2*A2_TO_BOHRS2 !Conversion to bohrs^2
      oneRIJ=one/sqrt(r2) !1/sqrt is faster than doing sqrt.
!      RIJ=one/oneRIJ
      RIJ=R2*oneRIJ
      SQRTAEE=one/sqrt(RR2+(qm2_params%multip_2c_elec_params(3,iqm)+rho_pm3mmx)**2)

      qmmm_erep = AU_TO_EV*SQRTAEE

      EE=-A2_TO_BOHRS2xAU_TO_EV*SQRTAEE*SQRTAEE*SQRTAEE
!        S-orbital of QM atom:
      DGX=vec_qm_mm1*EE
      DGY=vec_qm_mm2*EE
      DGZ=vec_qm_mm3*EE

      C1=qm_atom_core*mm_charge

      ENUC = qmmm_erep*C1
      ENUC = ENUC + ENUC*SIGN(1.0D0,mm_charge)*(-exp(-RIJ*sf1) + exp(-RIJ*sf2))

      FNUCx = dgx*C1
      FNUCy = dgy*C1
      FNUCz = dgz*C1

      temp_exp3 = SIGN(1.0D0,mm_charge)*exp(-RIJ*sf1)*C1
      temp_exp4 = SIGN(1.0D0,mm_charge)*exp(-RIJ*sf2)*C1
      corrfx =  dgx*(-temp_exp3+temp_exp4) + vec_qm_mm1*onerij*qmmm_erep &
                                             *(temp_exp3*sf1-temp_exp4*sf2)
      corrfy =  dgy*(-temp_exp3+temp_exp4) + vec_qm_mm2*onerij*qmmm_erep &
                                             *(temp_exp3*sf1-temp_exp4*sf2)
      corrfz =  dgz*(-temp_exp3+temp_exp4) + vec_qm_mm3*onerij*qmmm_erep &
                                             *(temp_exp3*sf1-temp_exp4*sf2)

      FNUCx = FNUCx + corrfx
      FNUCy = FNUCy + corrfy
      FNUCz = FNUCz + corrfz

    else  ! Follow normal QM/MM procedure.

      if (qmmm_nml%qmmmrij_incore) then
        oneRIJ=qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)
        RIJ=qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)
        EXP1=qm2_rij_eqns%qmmmrijdata(QMMMEXP1,loop_count)
        EXP2=qm2_rij_eqns%qmmmrijdata(QMMMEXP2,loop_count)
        SQRTAEE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count)
      else
        r2 = vec_qm_mm1*vec_qm_mm1+vec_qm_mm2*vec_qm_mm2+vec_qm_mm3*vec_qm_mm3
        RR2=r2*A2_TO_BOHRS2 !Conversion to bohrs^2
        oneRIJ=one/sqrt(r2) !1/sqrt is faster than doing sqrt.
!        RIJ=one/oneRIJ
        RIJ=R2*oneRIJ
        EXP1 = exp(-alpa*rij)
        EXP2 = exp(-ALPH_MM*rij)
        SQRTAEE=one/sqrt(RR2+qm2_params%multip_2c_elec_params(3,iqm)**2)
      end if

      qmmm_erep = AU_TO_EV*SQRTAEE

      EE=-A2_TO_BOHRS2xAU_TO_EV*SQRTAEE*SQRTAEE*SQRTAEE
!        S-orbital of QM atom:
      DGX=vec_qm_mm1*EE
      DGY=vec_qm_mm2*EE
      DGZ=vec_qm_mm3*EE

      C1=qm_atom_core*mm_charge

      ENUC = qmmm_erep*C1
      scale= ABS((EXP1+EXP2)*ENUC)

!     ****   START OF THE AM1 and PM3 RM1 etc SPECIFIC DERIVATIVE CODE   ***
!     ANALYT=-A*(1/(R*R)+2.D0*B*(R-C)/R)*EXP(-B*(R-C)**2)
!     ANALYTICAL DERIVATIVES

      if ( qmmm_nml%qmmm_int==2 .and. (qmmm_struct%AM1_OR_PM3 .or. qmmm_nml%qmtheory%PM6) ) then
        qmitype = qmmm_struct%qm_atom_type(iqm)
        ANAM1=zero
        oner2=oneRIJ*oneRIJ
        do I=1,qm2_params%num_fn(qmitype)
          temp_real=RIJ-qm2_params%FN3(i,qmitype)
          temp_real2=qm2_params%FN2(i,qmitype)*temp_real*temp_real
          if (temp_real2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
            ANAM1=ANAM1+qm2_params%FN1(i,qmitype)* &
                 (oner2+two*qm2_params%FN2(i,qmitype)*temp_real*oneRIJ)*EXP(-temp_real2)
            scale=scale+qm2_params%FN1(i,qmitype)*EXP(-temp_real2)*C1*oneRIJ
          end if
        end do
        ANAM1=-ANAM1*c1*onerij
        FNUCX=ANAM1*vec_qm_mm1
        FNUCY=ANAM1*vec_qm_mm2
        FNUCZ=ANAM1*vec_qm_mm3
      else
         FNUCX=zero; FNUCY=zero; FNUCZ=zero
      endif

      ENUC = ENUC + scale

!     ****   end OF THE AM1 and PM3 SPECIFIC DERIVATIVE CODE   ***

      EXP3 = (EXP1+EXP2)*abs(c1)
      EXP4 = qmmm_erep*onerij*(alpa*EXP1 + ALPH_MM*EXP2)*abs(c1)
      FNUCX = FNUCX+dgx*c1-vec_qm_mm1*EXP4+dgX*EXP3
      FNUCY = FNUCY+dgy*c1-vec_qm_mm2*EXP4+dgy*EXP3
      FNUCZ = FNUCZ+dgz*c1-vec_qm_mm3*EXP4+dgz*EXP3

    endif  ! qtw - end OF if(qmmm_int == 3) 

! end --- THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM

!     MM CORE AFFECTING AO'S ON QM ATOM.
      mm_charge=-mm_charge*psum_light
      FABX=mm_charge*DGX
      FABY=mm_charge*DGY
      FABZ=mm_charge*DGZ

      pair_force(1) = (FABX+FNUCX)*EV_TO_KCAL
      pair_force(2) = (FABY+FNUCY)*EV_TO_KCAL 
      pair_force(3) = (FABZ+FNUCZ)*EV_TO_KCAL 

      if (qmmm_nml%qmmm_switch) then
         if (rij <= qmmm_nml%r_switch_lo) then
            f_switch = one
            df_switch= zero
         else if (rij < qmmm_nml%r_switch_hi) then  ! qtw
            f_switch = qm2_switch_func(RIJ, qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi)
            df_switch= qm2_switch_derv(RIJ, qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi)
         else
            f_switch = zero
            df_switch= zero
         end if

         ! electronic energy in EV
         ELEC = -xyz_mm(4)*qmmm_erep*psum_light + ENUC
         ! converts EV to kcal
         ELEC = ELEC * EV_TO_KCAL

         AA = ELEC * oneRIJ * df_switch
         pair_force(1) = pair_force(1)*f_switch + vec_qm_mm1*AA
         pair_force(2) = pair_force(2)*f_switch + vec_qm_mm2*AA
         pair_force(3) = pair_force(3)*f_switch + vec_qm_mm3*AA

         pot_tmp = qm_charge*xyz_mm(4)*oneRIJ*AU_TO_EV*EV_TO_KCAL*BOHRS_TO_A
                                             !converts AU to kcal
         AA=oneRIJ*df_switch*pot_tmp

         pair_force(1) = pair_force(1) - vec_qm_mm1*AA 
         pair_force(2) = pair_force(2) - vec_qm_mm2*AA 
         pair_force(3) = pair_force(3) - vec_qm_mm3*AA 

         AA = pot_tmp*(one-f_switch)*oneRIJ*oneRIJ

         pair_force(1) = pair_force(1) - vec_qm_mm1*AA
         pair_force(2) = pair_force(2) - vec_qm_mm2*AA
         pair_force(3) = pair_force(3) - vec_qm_mm3*AA

      end if

end subroutine qm2_deriv_qmmm_light

subroutine qm2_deriv_qmmm_heavy(qmmm_struct, iqm,loop_count,psum,xyz_qm,xyz_mm,n_atomic_orb,pair_force, &
                         qm_atom_core,alpa,qm_charge)

!     This routine computes the analytical energy derivatives for the QM-MM
!     interaction energy arising from a single QM-MM pair.  The contibutions
!     to the derivatives come from the electron-core and core-core interactions.

!     Variable Definitions:
!
!     psum   - Density matrix elements for orbitals centered on QM atom.
!     xyz_qm - Cartesian coordinates of QM atom.
!     xyz_mm - Cartesian coordinates of MM atom. and charge
! n_atomic_orb - Number of atomic orbitals
! pair_force - Energy derivatives in the x, y, and z directions for
!              the interaction of the QM-MM atom pair.  The algebraic
!              sign of pair_force corresponds to dE/dR for the MM atom and
!              -dE/dR for the QM atom.
!
!  Current Version: Ross Walker (TSRI, 2004)
    use ElementOrbitalIndex, only : ElementSymbol, MaxValenceOrbitals, MaxValenceDimension
    use qmmm_module        , only : qmmm_nml,qm2_params, qm2_struct, qm2_rij_eqns, &
                                    alph_mm, AXIS_TOL, EXPONENTIAL_CUTOFF
    use constants          , only : A_TO_BOHRS, A2_TO_BOHRS2, AU_TO_EV, HALF_AU_TO_EV, FOURTH_AU_TO_EV, &
                                    A2_TO_BOHRS2xAU_TO_EV, one, zero, four, two, half, EV_TO_KCAL, &
                                    BOHRS_TO_A
    use qmmm_struct_module, only : qmmm_struct_type

    implicit none

! ON return, pair_force HOLDS ANALYTICAL DERIVATIVES                                   

!Passed in
      type(qmmm_struct_type), intent(in) :: qmmm_struct
      _REAL_, intent(in) :: xyz_qm(3),xyz_mm(4),psum(45)
      _REAL_, intent(out) :: pair_force(3)
      integer, intent(in) :: iqm, loop_count, n_atomic_orb
      _REAL_, intent(in) ::  qm_atom_core, alpa
      _REAL_, intent(in) :: qm_charge

!Local
      _REAL_ FABX, FABY, FABZ ,FNUCX, FNUCY, FNUCZ
      _REAL_ r2, rij, onerij, rr2, rr, one_rija0ev
      _REAL_ vec_qm_mm1, vec_qm_mm2, vec_qm_mm3
      _REAL_ c1, bb, mm_charge
      _REAL_ sqrtaee, dze, qzze, qxxe
      _REAL_ SQRTRRMDDADE, SQRTRRADDADE, SQRTRR2AQE
      _REAL_ SQRTRRAQQAQE, SQRTRRMQQAQE, SQRTRRAQQ2AQE
      _REAL_ Xtdx(MaxValenceOrbitals), Xtdy(MaxValenceOrbitals), Xtdz(MaxValenceOrbitals)
      _REAL_ Ytdx(MaxValenceOrbitals), Ytdy(MaxValenceOrbitals), Ytdz(MaxValenceOrbitals)
      _REAL_ Ztdx(MaxValenceOrbitals), Ztdy(MaxValenceOrbitals), Ztdz(MaxValenceOrbitals)
      _REAL_ TX(MaxValenceOrbitals),TY(MaxValenceOrbitals),TZ(MaxValenceOrbitals), TZ3i, TZ3i2
      _REAL_ RXY2, RYZ2, RZX2, oneRXY
      logical LRXY2, LRYZ2, LRZX2
      _REAL_ TERMX, TERMY, TERMZ
      _REAL_ DGX(MaxValenceOrbitals), DGY(MaxValenceOrbitals), DGZ(MaxValenceOrbitals)
      _REAL_ DRX(MaxValenceDimension)
      _REAL_ DRY(MaxValenceDimension)
      _REAL_ DRZ(MaxValenceDimension)
      _REAL_ EXP1, EXP2, EXP3, EXP4
      _REAL_ DD, QQ
      _REAL_ qm_mm_e_repul(4)
      _REAL_ :: anam1, oner2, temp_real, temp_real2, tmp1, tmp2
      _REAL_ :: corrfx, corrfy, corrfz, temp_exp3, temp_exp4
      _REAL_ :: sf1, sf2, rho_pm3mmx
      _REAL_ :: qm2_switch_func, qm2_switch_derv, qm2_HELECT  ! functions
      _REAL_ :: E1B(10)
      _REAL_ :: f_switch, df_switch, ELEC, ENUC, scale, AA, pot_tmp
      _REAL_ :: CHGMM_RI2, CHGMM_RI3, CHGMM_RI4
      _REAL_ :: vasya;
      integer:: qmitype
      integer isp, m, n, mn
      integer k, kk, l, ll, i

!****************************************************************************
!*                                                                          *
!*             CALCULATION OF ANALYTICAL DERIVATIVES                        *
!*                                                                          *
!* Routine inlined and optimised by Ross Walker and Mike Crowley (TSRI 2004)*
!****************************************************************************
                                                                               

      vec_qm_mm1=xyz_qm(1) - xyz_mm(1)
      vec_qm_mm2=xyz_qm(2) - xyz_mm(2)
      vec_qm_mm3=xyz_qm(3) - xyz_mm(3)
      mm_charge = xyz_mm(4)
#include "qm2_array_locations.h"
      !if ((qmmm_nml%qmmm_int == 3 .or. qmmm_nml%qmmm_int == 4) .and. qmmm_nml%qmtheory%PM3) then  ! PM3/MM*
      if (qmmm_struct%PM3MMX_INTERFACE) then
        ! PM3/MM* - MODIFIED QM-MM INTERFACE

        qmitype = qmmm_struct%qm_atom_type(iqm)
        ! rho_pm3mmx in qmmm_int==3 is always zero
        rho_pm3mmx = qm2_params%rho_pm3mmx(qmitype)

        r2 = vec_qm_mm1*vec_qm_mm1+vec_qm_mm2*vec_qm_mm2+vec_qm_mm3*vec_qm_mm3
        RR2=r2*A2_TO_BOHRS2 !Conversion to bohrs^2
        oneRIJ=one/sqrt(r2) !1/sqrt is faster than doing sqrt.
        RIJ=one/oneRIJ
        RR=RIJ*A_TO_BOHRS
        one_rija0ev = A_TO_BOHRS * oneRIJ*AU_TO_EV
        EXP1 = exp(-alpa*rij)
        EXP2 = exp(-ALPH_MM*rij)
        SQRTAEE=one/sqrt(RR2+(qm2_params%multip_2c_elec_params(3,iqm)+rho_pm3mmx)**2)
!SP-ATOM specific
          DD=qm2_params%multip_2c_elec_params(1,iqm)
          QQ=qm2_params%multip_2c_elec_params(2,iqm)
          tmp1 = (qm2_params%multip_2c_elec_params(4,iqm)+rho_pm3mmx)**2
          tmp2 = (qm2_params%multip_2c_elec_params(5,iqm)+rho_pm3mmx)**2
          SQRTRRADDADE=one/SQRT((RR+DD)**2+tmp1)
          SQRTRRMDDADE=one/SQRT((RR-DD)**2+tmp1)
          SQRTRR2AQE=one/SQRT(RR2+tmp2)
          SQRTRRAQQAQE=one/SQRT((RR+two*QQ)**2+tmp2)
          SQRTRRMQQAQE=one/SQRT((RR-two*QQ)**2+tmp2)
          SQRTRRAQQ2AQE=one/SQRT(RR2+(four*(QQ**2)+tmp2))
!end sp_atom specific
      else if (qmmm_nml%qmmmrij_incore) then
        oneRIJ=qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)
        RIJ=qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)
        one_rija0ev = oneRIJ*A_TO_BOHRS*AU_TO_EV
        EXP1=qm2_rij_eqns%qmmmrijdata(QMMMEXP1,loop_count)
        EXP2=qm2_rij_eqns%qmmmrijdata(QMMMEXP2,loop_count)
        SQRTAEE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count)
!SP-ATOM Specific
          SQRTRRADDADE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRADDADE,loop_count)
          SQRTRRMDDADE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMDDADE,loop_count)
          SQRTRR2AQE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRR2AQE,loop_count)
          SQRTRRAQQAQE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQAQE,loop_count)
          SQRTRRMQQAQE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMQQAQE,loop_count)
          SQRTRRAQQ2AQE=qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQ2AQE,loop_count)
!End SP-ATOM specific
      else
        r2 = vec_qm_mm1*vec_qm_mm1+vec_qm_mm2*vec_qm_mm2+vec_qm_mm3*vec_qm_mm3
        RR2=r2*A2_TO_BOHRS2 !Conversion to bohrs^2
        oneRIJ=one/sqrt(r2) !1/sqrt is faster than doing sqrt.
        RIJ=one/oneRIJ
        RR=RIJ*A_TO_BOHRS
        one_rija0ev = A_TO_BOHRS * oneRIJ*AU_TO_EV
        EXP1 = exp(-alpa*rij)
        EXP2 = exp(-ALPH_MM*rij)
        SQRTAEE=one/sqrt(RR2+qm2_params%multip_2c_elec_params(3,iqm)**2)
!SP-ATOM specific
          DD=qm2_params%multip_2c_elec_params(1,iqm)
          QQ=qm2_params%multip_2c_elec_params(2,iqm)
          tmp1 = qm2_params%multip_2c_elec_params(4,iqm)**2
          tmp2 = qm2_params%multip_2c_elec_params(5,iqm)**2
          SQRTRRADDADE=one/SQRT((RR+DD)**2+tmp1)
          SQRTRRMDDADE=one/SQRT((RR-DD)**2+tmp1)
          SQRTRR2AQE=one/SQRT(RR2+tmp2)
          SQRTRRAQQAQE=one/SQRT((RR+two*QQ)**2+tmp2)
          SQRTRRMQQAQE=one/SQRT((RR-two*QQ)**2+tmp2)
          SQRTRRAQQ2AQE=one/SQRT(RR2+(four*(QQ**2)+tmp2))
!end sp_atom specific
      end if

      qm_mm_e_repul(1) = AU_TO_EV*SQRTAEE
      qm_mm_e_repul(2) = HALF_AU_TO_EV*(SQRTRRADDADE - SQRTRRMDDADE)
      qm_mm_e_repul(3) = qm_mm_e_repul(1) + FOURTH_AU_TO_EV*(SQRTRRAQQAQE + SQRTRRMQQAQE)-HALF_AU_TO_EV*SQRTRR2AQE
      qm_mm_e_repul(4) = qm_mm_e_repul(1) + HALF_AU_TO_EV*(SQRTRRAQQ2AQE - SQRTRR2AQE)

!     Returns the
!     derivatives of the electron-core interaction energies in a local
!     diatomic frame.  At most only four terms are computed because there
!     are at most four unique electron-core interactions for the QM-MM
!     pair: (ss|  ), (so|  ), (oo|  ), (pp|  ).  Note that it is not
!     necessary to specify whether an x, y, or z derivative is being
!     evaluated because the formula is the same for all three directions.
!      one_rija0 = one_rija0

      SQRTAEE = -SQRTAEE*SQRTAEE*SQRTAEE*A2_TO_BOHRS2xAU_TO_EV
!        S-orbital of QM atom:
      DGX(1) = vec_qm_mm1*SQRTAEE
      DGY(1) = vec_qm_mm2*SQRTAEE
      DGZ(1) = vec_qm_mm3*SQRTAEE

      DD=qm2_params%multip_2c_elec_params(1,iqm)*one_rija0ev
      SQRTRRADDADE=SQRTRRADDADE*SQRTRRADDADE*SQRTRRADDADE
      SQRTRRMDDADE=SQRTRRMDDADE*SQRTRRMDDADE*SQRTRRMDDADE

      DZE = half*(A2_TO_BOHRS2xAU_TO_EV*(SQRTRRMDDADE-SQRTRRADDADE)-DD*(SQRTRRMDDADE+SQRTRRADDADE))

      SQRTRR2AQE=SQRTRR2AQE*SQRTRR2AQE*SQRTRR2AQE
      SQRTRRAQQ2AQE=SQRTRRAQQ2AQE*SQRTRRAQQ2AQE*SQRTRRAQQ2AQE
      QXXE  = SQRTAEE+half*A2_TO_BOHRS2xAU_TO_EV*(SQRTRR2AQE-SQRTRRAQQ2AQE)

      SQRTRRAQQAQE=SQRTRRAQQAQE*SQRTRRAQQAQE*SQRTRRAQQAQE
      SQRTRRMQQAQE=SQRTRRMQQAQE*SQRTRRMQQAQE*SQRTRRMQQAQE

      QZZE  = SQRTAEE+half*(one_rija0ev*qm2_params%multip_2c_elec_params(2,iqm)*(SQRTRRMQQAQE-SQRTRRAQQAQE) - &
                            half*A2_TO_BOHRS2xAU_TO_EV*(SQRTRRMQQAQE+SQRTRRAQQAQE-two*SQRTRR2AQE))

      DGX(2)=vec_qm_mm1*DZE
      DGX(3)=vec_qm_mm1*QZZE
      DGX(4)=vec_qm_mm1*QXXE

      DGY(2)=vec_qm_mm2*DZE
      DGY(3)=vec_qm_mm2*QZZE
      DGY(4)=vec_qm_mm2*QXXE

      DGZ(2)=vec_qm_mm3*DZE
      DGZ(3)=vec_qm_mm3*QZZE
      DGZ(4)=vec_qm_mm3*QXXE

!     This routine
!     takes the electron-core integral derivatives DG which have been
!     computed for a local frame and rotates them to the molecular frame.
!     The rotated derivatives are returned in DR

!     Determines the transformation (TX,TY,TZ) and derivatives of the transformation
!     (TDX,TDY,TDZ) involved in rotating from a local frame to the molecular frame 
!     for calculation of electron-core interactions.
      RXY2=vec_qm_mm1*vec_qm_mm1+vec_qm_mm2*vec_qm_mm2   
      RYZ2=vec_qm_mm2*vec_qm_mm2+vec_qm_mm3*vec_qm_mm3
      RZX2=vec_qm_mm3*vec_qm_mm3+vec_qm_mm1*vec_qm_mm1
      LRXY2 = RXY2 < AXIS_TOL
      LRYZ2 = RYZ2 < AXIS_TOL
      LRZX2 = RZX2 < AXIS_TOL

!     Zeros entire array of 3
      XTDX=zero
      YTDX=zero
      ZTDX=zero
      XTDY=zero
      YTDY=zero
      ZTDY=zero
      XTDZ=zero
      YTDZ=zero
      ZTDZ=zero

      if(.NOT.(LRXY2 .OR. LRYZ2 .OR. LRZX2)) then 
        oneRXY = one/sqrt(RXY2)
      !Ross Walker + Mike Crowley - rearranged order here slightly for speed.
        TZ(3)=oneRIJ/oneRXY
        TZ3i = one/TZ(3)  !Inverse of TZ(3) to avoid other divisions.
        TZ3i2 = TZ3i*TZ3i  !Square of 1/TZ(3)

        TX(1)=vec_qm_mm1*oneRIJ
        TX(2)=vec_qm_mm2*oneRIJ
        TX(3)=vec_qm_mm3*oneRIJ

        TY(1)=-TX(2)*SIGN(one,TX(1))*TZ3i
        TY(2)=ABS(TX(1)*TZ3i)
        TY(3)=zero

        TZ(1)=-TX(1)*TX(3)*TZ3i
        TZ(2)=-TX(2)*TX(3)*TZ3i

        TERMX = TX(1)*oneRIJ
        TERMY = TX(2)*oneRIJ
        TERMZ = TX(3)*oneRIJ

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
!        XTDY(3)=0.0D0
!        YTDY(3)=0.0D0
!        ZTDY(3)=0.0D0

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
!     MOLECULAR Z AXIS IS PARALLEL TO DIATOMIC Z AXIS
        TX(1)=zero
        TX(2)=zero
        TX(3)=sign(one,vec_qm_mm3)
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
!     MOLECULAR X AXIS IS PARALLEL TO DIATOMIC Z AXIS
        TX(1)=sign(one,vec_qm_mm1)
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
!     MOLECULAR Y AXIS IS PARALLEL TO DIATOMIC Z AXIS
        TX(1)=zero
        TX(2)=sign(one,vec_qm_mm2)
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


      ISP=0
      do K=1,n_atomic_orb
         KK=K-1
         do L=K,n_atomic_orb
            LL=L-1
            ISP=ISP+1
            if(LL == 0) then
               DRX(ISP)=DGX(1)                       !(SS/SS)
               DRY(ISP)=DGY(1)                       !(SS/SS)
               DRZ(ISP)=DGZ(1)                       !(SS/SS)
           elseif(KK == 0) then
               DRX(ISP)=DGX(2)*TX(LL)+qm_mm_e_repul(2)*XTDX(LL)   !(SP/SS)
               DRY(ISP)=DGY(2)*TX(LL)+qm_mm_e_repul(2)*YTDX(LL)   !(SP/SS)
               DRZ(ISP)=DGZ(2)*TX(LL)+qm_mm_e_repul(2)*ZTDX(LL)   !(SP/SS)
            else
               DRX(ISP)=DGX(3)*TX(KK)*TX(LL)        &  !(PP/SS)
              +qm_mm_e_repul(3)*(XTDX(KK)*TX(LL)+TX(KK)*XTDX(LL))  &
              +DGX(4)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))   &
              +qm_mm_e_repul(4)*(XTDY(KK)*TY(LL)+TY(KK)*XTDY(LL)   &
              +XTDZ(KK)*TZ(LL)+TZ(KK)*XTDZ(LL))

               DRY(ISP)=DGY(3)*TX(KK)*TX(LL)        &  !(PP/SS)
              +qm_mm_e_repul(3)*(YTDX(KK)*TX(LL)+TX(KK)*YTDX(LL))  &
              +DGY(4)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))   &
              +qm_mm_e_repul(4)*(YTDY(KK)*TY(LL)+TY(KK)*YTDY(LL)   &
              +YTDZ(KK)*TZ(LL)+TZ(KK)*YTDZ(LL))

               DRZ(ISP)=DGZ(3)*TX(KK)*TX(LL)        &  !(PP/SS)
              +qm_mm_e_repul(3)*(ZTDX(KK)*TX(LL)+TX(KK)*ZTDX(LL))  &
              +DGZ(4)*(TY(KK)*TY(LL)+TZ(KK)*TZ(LL))   &
              +qm_mm_e_repul(4)*(ZTDY(KK)*TY(LL)+TY(KK)*ZTDY(LL)   &
              +ZTDZ(KK)*TZ(LL)+TZ(KK)*ZTDZ(LL))
            endif
         end do
      end do

! BEGIN --- THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM
!           Note: This could technically be done at the same time we do the
!                 energy in hcore_qmmm.

    !if ((qmmm_nml%qmmm_int == 3 .or. qmmm_nml%qmmm_int == 4) .and. qmmm_nml%qmtheory%PM3) then  ! PM3/MM*
    if (qmmm_struct%PM3MMX_INTERFACE) then
      ! PM3/MM* - MODIFIED QM-MM INTERFACE

      qmitype = qmmm_struct%qm_atom_type(iqm)

      sf1 = qm2_params%scale_factor1_pm3mmx(1,qmitype)
      sf2 = qm2_params%scale_factor2_pm3mmx(1,qmitype)

      C1=qm_atom_core*mm_charge

      ENUC = qm_mm_e_repul(1)*C1
      ENUC = ENUC + ENUC*SIGN(1.0D0,mm_charge)*(-exp(-RIJ*sf1) + exp(-RIJ*sf2))

      FNUCx = dgx(1)*C1
      FNUCy = dgy(1)*C1
      FNUCz = dgz(1)*C1

      temp_exp3 = SIGN(1.0D0,mm_charge)*exp(-RIJ*sf1)*C1
      temp_exp4 = SIGN(1.0D0,mm_charge)*exp(-RIJ*sf2)*C1
      corrfx =  dgx(1)*(-temp_exp3+temp_exp4) + vec_qm_mm1*onerij*qm_mm_e_repul(1) &
                                                *(temp_exp3*sf1-temp_exp4*sf2)
      corrfy =  dgy(1)*(-temp_exp3+temp_exp4) + vec_qm_mm2*onerij*qm_mm_e_repul(1) &
                                                *(temp_exp3*sf1-temp_exp4*sf2)
      corrfz =  dgz(1)*(-temp_exp3+temp_exp4) + vec_qm_mm3*onerij*qm_mm_e_repul(1) &
                                                *(temp_exp3*sf1-temp_exp4*sf2)

      FNUCx = FNUCx + corrfx
      FNUCy = FNUCy + corrfy
      FNUCz = FNUCz + corrfz

    else  ! then FOLLOW THE NORMAL QM/MM PROCEDURE

      C1=qm_atom_core*mm_charge

      ENUC = qm_mm_e_repul(1)*C1
      scale= ABS((EXP1+EXP2)*ENUC)

!     ****   START OF THE AM1 and PM3 SPECIFIC DERIVATIVE CODE   ***
!     ANALYT=-A*(1/(R*R)+2.D0*B*(R-C)/R)*EXP(-B*(R-C)**2)
!     ANALYTICAL DERIVATIVES

      if ( qmmm_nml%qmmm_int==2 .and. (qmmm_struct%AM1_OR_PM3 .or. qmmm_nml%qmtheory%PM6) ) then
        qmitype = qmmm_struct%qm_atom_type(iqm)
        ANAM1=zero
        oner2=oneRIJ*oneRIJ
        do I=1,qm2_params%num_fn(qmitype)
          temp_real=RIJ-qm2_params%FN3(i,qmitype)
          temp_real2=qm2_params%FN2(i,qmitype)*temp_real*temp_real
          if (temp_real2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
            ANAM1=ANAM1+qm2_params%FN1(i,qmitype)* &
                 (oner2+two*qm2_params%FN2(i,qmitype)*temp_real*oneRIJ)*EXP(-temp_real2)
            scale=scale+qm2_params%FN1(i,qmitype)*EXP(-temp_real2)*C1*oneRIJ
          end if
        end do
        ANAM1=-ANAM1*c1*onerij
        FNUCX=ANAM1*vec_qm_mm1
        FNUCY=ANAM1*vec_qm_mm2
        FNUCZ=ANAM1*vec_qm_mm3
     else
        FNUCX=zero; FNUCY=zero; FNUCZ=zero
     endif

     ENUC = ENUC + scale

!     ****   end OF THE AM1 and PM3 SPECIFIC DERIVATIVE CODE   ***

      FNUCX = FNUCX+dgX(1)*c1
      FNUCY = FNUCY+dgy(1)*c1
      FNUCZ = FNUCZ+dgz(1)*c1
      EXP3 = (EXP1+EXP2)*abs(c1)
      EXP4 = qm_mm_e_repul(1)*onerij*(alpa*EXP1 + ALPH_MM*EXP2)*abs(c1)
      FNUCX = FNUCX-vec_qm_mm1*EXP4
      FNUCY = FNUCY-vec_qm_mm2*EXP4
      FNUCZ = FNUCZ-vec_qm_mm3*EXP4

      FNUCX = FNUCX+dgX(1)*EXP3
      FNUCY = FNUCY+dgy(1)*EXP3
      FNUCZ = FNUCZ+dgz(1)*EXP3

    endif  ! qtw - end OF if(qmmm_int == 3) 

! end --- THE FIRST DERIVATIVE OF NUCLEAR REPULSION TERM

      FABX=zero
      FABY=zero
      FABZ=zero
!     MM CORE AFFECTING AO'S ON QM ATOM.
      ISP=0
      do M=1,n_atomic_orb
         BB=one
         do N=M,n_atomic_orb
            MN=M+qm2_params%pascal_tri1(N)
            ISP=ISP+1
            FABX=FABX-BB*mm_charge*PSUM(MN)*DRX(ISP)
            FABY=FABY-BB*mm_charge*PSUM(MN)*DRY(ISP)
            FABZ=FABZ-BB*mm_charge*PSUM(MN)*DRZ(ISP)                    
            BB=two
         end do
      end do

      pair_force(1) = (FABX+FNUCX)*EV_TO_KCAL                                          
      pair_force(2) = (FABY+FNUCY)*EV_TO_KCAL                                           
      pair_force(3) = (FABZ+FNUCZ)*EV_TO_KCAL                                           

      if (qmmm_nml%qmmm_switch) then
         if (rij <= qmmm_nml%r_switch_lo) then  ! qtw
            f_switch = one
            df_switch= zero
         else if (rij < qmmm_nml%r_switch_hi) then
            f_switch = qm2_switch_func(RIJ, qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi)
            df_switch= qm2_switch_derv(RIJ, qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi)
         else
            f_switch = zero
            df_switch= zero
         end if

         ! Calculate energy
         E1B(1) = -mm_charge*qm_mm_e_repul(1)

         CHGMM_RI2 = -mm_charge*qm_mm_e_repul(2)
         CHGMM_RI3 = -mm_charge*qm_mm_e_repul(3)
         CHGMM_RI4 = -mm_charge*qm_mm_e_repul(4)

         E1B(2) = CHGMM_RI2*TX(1)
         E1B(3) = CHGMM_RI3*TX(1)*TX(1)+CHGMM_RI4*((TY(1)*TY(1))+(TZ(1)*TZ(1)))
         E1B(4) = CHGMM_RI2*TX(2)
         E1B(5) = CHGMM_RI3*TX(2)*TX(1)+CHGMM_RI4*((TY(2)*TY(1))+(TZ(2)*TZ(1)))
         E1B(6) = CHGMM_RI3*TX(2)*TX(2)+CHGMM_RI4*((TY(2)*TY(2))+(TZ(2)*TZ(2)))
         E1B(7) = CHGMM_RI2*TX(3)
         E1B(8) = CHGMM_RI3*TX(3)*TX(1)+CHGMM_RI4*TZ(3)*TZ(1)
         E1B(9) = CHGMM_RI3*TX(3)*TX(2)+CHGMM_RI4*TZ(3)*TZ(2)
         E1B(10)= CHGMM_RI3*TX(3)*TX(3)+CHGMM_RI4*TZ(3)*TZ(3)

         ! electronic energy in EV
         ELEC = qm2_HELECT(3, psum, E1B, E1B,vasya) + ENUC
         ! converts EV to kcal
         ELEC = ELEC * EV_TO_KCAL

         AA = ELEC * oneRIJ * df_switch
         pair_force(1) = pair_force(1)*f_switch + vec_qm_mm1*AA
         pair_force(2) = pair_force(2)*f_switch + vec_qm_mm2*AA
         pair_force(3) = pair_force(3)*f_switch + vec_qm_mm3*AA

         pot_tmp = qm_charge*mm_charge*oneRIJ*AU_TO_EV*EV_TO_KCAL*BOHRS_TO_A
                                             !converts AU to kcal
         AA=oneRIJ*df_switch*pot_tmp

         pair_force(1) = pair_force(1) - vec_qm_mm1*AA 
         pair_force(2) = pair_force(2) - vec_qm_mm2*AA 
         pair_force(3) = pair_force(3) - vec_qm_mm3*AA 

         AA = pot_tmp*(one-f_switch)*oneRIJ*oneRIJ

         pair_force(1) = pair_force(1) - vec_qm_mm1*AA
         pair_force(2) = pair_force(2) - vec_qm_mm2*AA
         pair_force(3) = pair_force(3) - vec_qm_mm3*AA

      end if

end subroutine qm2_deriv_qmmm_heavy

