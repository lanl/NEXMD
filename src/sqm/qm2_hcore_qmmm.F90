! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
subroutine  qm2_hcore_qmmm(qm2_params,qmmm_nml,qmmm_mpi, qm2_rij_eqns, qm2_struct, qmmm_struct, H,ENUCLR,qm_xcrd)

! This routine is repsonsible for calculating the QM-MM interactions between
! QM-MM pairs that are in the QMMM pair list. This routine calls qm2_rotate_qmmm.
! It is responsible for updating the one-electron matix H with the QM electron
! MM core interaction. It also calculates the QM-MM core-core interactions that
! are returned as an addition to ENUCLR.
!
! QM core charge interacts with MM resp charge. The MM resp charge, in electron
! units instead of AMBER's internal units is stored as the 4th index of the qm_COORD array.
! Hence the qm_coord array runs: x1, y1, z1, q1, x2, y2, z2, q2...

! Current code: Ross Walker (TSRI, 2005)
! d-orbital extension,  (Taisung Lee, 2011)  

      use constants  , only : zero, one, A_TO_BOHRS, A2_TO_BOHRS2, BOHRS_TO_A, AU_TO_EV
      use qmmm_module, only : qm2_structure, alph_MM, qmmm_mpi_structure, qm2_rij_eqns_structure
      use Rotation   , only : GetRotationMatrix, RotateCore
      use qmmm_struct_module, only : qmmm_struct_type
      use qm2_params_module,  only : qm2_params_type
      use qmmm_nml_module   , only : qmmm_nml_type

      implicit none
	   type(qm2_params_type),intent(inout) :: qm2_params
	   type(qmmm_nml_type),intent(inout) :: qmmm_nml
	   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
           type(qm2_rij_eqns_structure),intent(inout) :: qm2_rij_eqns
 
      type(qmmm_struct_type), intent(inout) :: qmmm_struct
      type(qm2_structure),intent(inout) :: qm2_struct
      _REAL_, intent(inout) :: H(qm2_struct%matsize)
      _REAL_, intent(out)   :: enuclr
      _REAL_, intent(in)    :: qm_xcrd(4,qmmm_struct%qm_mm_pairs) 

! Keeps track of the nquant * ni_mm loop iterations
      _REAL_ qm_atom_coord(3)
      integer loop_count, i, j, ia, ib, i2, i1, j1, ii,jj
      logical s_atom, sp_atom, spd_atom
      _REAL_ enucij
      _REAL_ E1B(10), E1B_light
      
      ! d-orbital related
      integer :: firstIndexAO_i, n_orbs
      _REAL_  :: w(45), RI(22), core(10,2), rotationMatrix(15,45) 
      _REAL_  :: r2, r2InAu, rij, rijInAu, oneOverRij, scale, temp     
      _REAL_  :: qm2_switch_func  ! function
      _REAL_  :: f_switch, enuc_switch, pot_tmp


      ! Initializing array
      if (qmmm_nml%qmmm_switch) then
         qmmm_struct%switched_mmpot = zero
      end if

!     Loop over REAL QM atoms. - no need to worry about link atoms here
      loop_count=0
!      do i=1,qmmm_struct%nquant_nlink
      do i=qmmm_mpi%nquant_nlink_start,qmmm_mpi%nquant_nlink_end
!
!       Definitions:
!
!       orb_loc(1,i) = first atomic orbital on atom i
!       orb_loc(2,i) =  last atomic orbital on atom i
! iqm_atomic_numbers(i) =    atomic number for atom i
!
        s_atom   = (qm2_params%natomic_orbs(i) ==1)
        sp_atom  = (qm2_params%natomic_orbs(i) ==4)        
        spd_atom = (qm2_params%natomic_orbs(i) ==9)
        
        ia = qm2_params%orb_loc(1,i)
        qm_atom_coord(1) = qmmm_struct%qm_coords(1,i)
        qm_atom_coord(2) = qmmm_struct%qm_coords(2,i)
        qm_atom_coord(3) = qmmm_struct%qm_coords(3,i)

!       Loop over MM atoms that interact with QM atom i, update the
!       one-electron matrix H, and accumulate core-core repulsions
!       ENUCLR.

! Loop in steps of since qm_xcrd array is x,y,z,chg,x,y,z,chg...
! Duplicated code here but it avoids the if(sp_atom) inside the loop
!          Get the QM-MM interactions:
!          QM electrons - MM core ---> E1B
!          QM core - MM core ----> enucij 

        if (s_atom) then

          !It is a light atom (Hydrogen) - same notes as below for SP-ATOM
          do j=1,qmmm_struct%qm_mm_pairs
            loop_count = loop_count+1
            call qm2_rotate_qmmm_light(qm2_params, qmmm_nml, qm2_rij_eqns, qm2_struct,qmmm_struct, &
                                       loop_count,i,j,qm_atom_coord,qm_xcrd(1,j),E1B_light,enucij)
            ii = qm2_params%pascal_tri1(ia) + ia 
            H(ii) = H(ii) + E1B_light
            ENUCLR = ENUCLR + enucij
          end do
        
        else if (sp_atom) then        
 
          ib = qm2_params%orb_loc(2,i) !Will be same as ia for light atom
          do j=1,qmmm_struct%qm_mm_pairs
            loop_count = loop_count+1

            call qm2_rotate_qmmm_heavy(qm2_params,qmmm_nml, qm2_rij_eqns, qm2_struct,qmmm_struct, &
                                       loop_count,i,qm_atom_coord,qm_xcrd(1,j),E1B,enucij)
!          Add E1B to the corresponding diagonal elements of the one-electron
!          matrix H.
            i2 = 0
                        
            do i1=ia,ib
              ii = qm2_params%pascal_tri1(i1) + ia - 1
              do j1=ia,i1
                ii = ii + 1
                i2 = i2 + 1
                H(ii) = H(ii) + E1B(i2)
              end do 
            end do 
!         Add on QM core - MM core interactions.
            ENUCLR = ENUCLR + enucij

          end do

        else if (spd_atom) then
        
            firstIndexAO_i=qm2_params%orb_loc(1,i)
            ii=qm2_params%pascal_tri2(firstIndexAO_i)
            n_orbs=qm2_params%natomic_orbs(i)

            do j=1,qmmm_struct%qm_mm_pairs
                loop_count = loop_count+1

                r2=sum( (qm_xcrd(1:3,j)-qm_atom_coord(1:3)) **2) 
                rij=sqrt(r2)
                rijInAu=rij*A_TO_BOHRS
                oneOverRij=one/rij
                r2InAu=r2*A2_TO_BOHRS2

                call GetRotationMatrix(qm2_params,qm_xcrd(1:3,j)-qm_atom_coord(1:3), rotationMatrix, spd_atom)   
                call qm2_repp_d(qm2_params,qmmm_struct, qmmm_struct%qm_atom_type(i),0,rijInAu,RI,CORE,W,45,1,0) 
                core(1,2)=0.d0
                core=core*qm_xcrd(4,j)
                
                if (qmmm_nml%qmmm_switch) then
                   if (rij <= qmmm_nml%r_switch_lo) then
                     f_switch = one
                   else if (rij < qmmm_nml%r_switch_hi) then
                     f_switch = qm2_switch_func(RIJ, qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi)
                   else
                     f_switch = zero
                   end if
                   core=core*f_switch

                   pot_tmp = (one - f_switch)*qm_xcrd(4,j)*oneOverRij*AU_TO_EV*BOHRS_TO_A
                                                                     !converts (AU/A) to (ev/Bohrs)
                   qmmm_struct%switched_mmpot(i) = qmmm_struct%switched_mmpot(i) + pot_tmp

                   enuc_switch = qm2_params%core_chg(i)*pot_tmp
                else
                   enuc_switch = zero
                end if

                call RotateCore(firstIndexAO_i,0, n_orbs,0, ii, 0            &
                    ,core,rotationMatrix,H)
                 
                temp = -core(1,1)*qm2_params%core_chg(i)
                scale = abs(temp*exp(-qm2_params%cc_exp_params(i)*RIJ)+EXP(-ALPH_MM*RIJ))
                
                ENUCLR = ENUCLR + temp + scale + enuc_switch

           end do
          
        end if        

      end do

end subroutine qm2_hcore_qmmm

SUBROUTINE qm2_rotate_qmmm_light(qm2_params, qmmm_nml, qm2_rij_eqns, qm2_struct,qmmm_struct, &
                                 loop_count,IQM,jpair,xyz_qm,xyz_mm,E1B,ENUC)

!For light atoms
!-----------------------------------
!Written by Ross Walker (TSRI, 2005)
!-----------------------------------
!See heavy routine below for notes.

      use qmmm_module, only : qm2_structure, qm2_rij_eqns_structure,alph_mm, &
                              EXPONENTIAL_CUTOFF
      use constants  , only : A2_TO_BOHRS2, AU_TO_EV, BOHRS_TO_A, one, zero
      use qmmm_struct_module, only : qmmm_struct_type
      use qm2_params_module,  only : qm2_params_type
      use qmmm_nml_module   , only : qmmm_nml_type

      implicit none
! Passed in
      type(qmmm_struct_type), intent(inout) :: qmmm_struct
      type(qm2_structure),intent(inout) :: qm2_struct
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qm2_rij_eqns_structure),intent(inout) :: qm2_rij_eqns

      integer, intent(in) :: loop_count, iqm, jpair
      _REAL_, intent(in) :: xyz_qm(3), xyz_mm(4)
      _REAL_, intent(out) :: e1b, enuc

! Local
      _REAL_ :: RI
      _REAL_ :: oneBDDi1
      _REAL_ :: X1, X2, X3
      _REAL_ :: RR2, oneRIJ, RIJ, RIJ2
      _REAL_ :: scale, anam1, temp_real, temp_real2
      _REAL_ :: sf1, sf2, rho_pm3mmx
      _REAL_ :: qm2_switch_func  ! function
      _REAL_ :: f_switch, pot_tmp
      integer :: i, qmitype, mm_natomic
    
#include "qm2_array_locations.h"

    !PM3/MM* - MODIFIED QM-MM INTERFACE
    !if ((qmmm_nml%qmmm_int == 3 .or. qmmm_nml%qmmm_int == 4) .and. qmmm_nml%qmtheory%PM3) then  ! PM3/MM*
    if (qmmm_struct%PM3MMX_INTERFACE) then

      qmitype = qmmm_struct%qm_atom_type(iqm)

      ! Different scale factors are applied to H(QM)-H(MM) and H(QM)-heavy(MM)
      if (qmmm_struct%qm_mm_pair_atom_numbers(jpair) == 1) then
         ! H MM atom
         sf1 = qm2_params%scale_factor1_pm3mmx(1,qmitype)
         sf2 = qm2_params%scale_factor2_pm3mmx(1,qmitype)
      else
         ! Heavy MM atom
         sf1 = qm2_params%scale_factor1_pm3mmx(2,qmitype)
         sf2 = qm2_params%scale_factor2_pm3mmx(2,qmitype)
      end if
      ! rho_pm3mmx in qmmm_int==3 is always zero
      rho_pm3mmx = qm2_params%rho_pm3mmx(qmitype)

      X1 = xyz_qm(1) - xyz_mm(1)
      X2 = xyz_qm(2) - xyz_mm(2)
      X3 = xyz_qm(3) - xyz_mm(3)
      RIJ2=X1*X1+X2*X2+X3*X3
      RR2=RIJ2*A2_TO_BOHRS2
      oneRIJ=1.0d0/sqrt(RIJ2)
      RIJ =RIJ2*oneRIJ
      !oneBDDi1=qm2_params%multip_2c_elec_params(3,iqm)**2
      oneBDDi1=(qm2_params%multip_2c_elec_params(3,iqm)+rho_pm3mmx)**2
      RI = AU_TO_EV/sqrt(RR2+oneBDDi1)
      E1B=-xyz_mm(4)*RI
      ENUC = -qm2_params%core_chg(IQM)*E1B
      ! Added ENUC correction
      scale = -exp(-RIJ*sf1) + exp(-RIJ*sf2)
      ENUC = ENUC + SIGN(1.0D0,xyz_mm(4))*ENUC*scale

    else  ! Follow normal QM/MM procedure

      if (qmmm_nml%qmmmrij_incore) then
        scale=qm2_rij_eqns%qmmmrijdata(QMMMEXP1,loop_count)+qm2_rij_eqns%qmmmrijdata(QMMMEXP2,loop_count)
        RIJ = qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)
        oneRIJ = qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)
        RI = AU_TO_EV*qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count)
      else
        X1 = xyz_qm(1) - xyz_mm(1)
        X2 = xyz_qm(2) - xyz_mm(2)
        X3 = xyz_qm(3) - xyz_mm(3)
        RIJ2=X1*X1+X2*X2+X3*X3
        RR2=RIJ2*A2_TO_BOHRS2
        oneRIJ=1.0d0/sqrt(RIJ2)
        RIJ =RIJ2*oneRIJ
        scale=EXP(-qm2_params%cc_exp_params(iqm)*RIJ)+EXP(-ALPH_MM*RIJ)
        oneBDDi1=qm2_params%multip_2c_elec_params(3,iqm)**2
        RI = AU_TO_EV/sqrt(RR2+oneBDDi1)
      end if

      E1B=-xyz_mm(4)*RI
      ENUC = -qm2_params%core_chg(IQM)*E1B
      scale = ABS(scale*ENUC)

!Add in the extra core core repulsion terms if qmmm_int==2. This is equivalent of the
!way CHARMM and DYNAMO do QMMM interactions.
      if(qmmm_nml%qmmm_int == 2) then
         if (qmmm_struct%AM1_OR_PM3 .or. qmmm_nml%qmtheory%PM6) then
!          Add gaussians.
           if (qmmm_nml%qmmmrij_incore) then
             rij = qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)
             onerij = qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)
           end if
           qmitype = qmmm_struct%qm_atom_type(iqm)
           anam1=0.0d0
           do I=1,qm2_params%num_fn(qmitype)
              temp_real=RIJ-qm2_params%FN3(i,qmitype)
              temp_real2=qm2_params%FN2(i,qmitype)*temp_real*temp_real
              if (temp_real2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
                anam1=anam1+qm2_params%FN1(i,qmitype)*EXP(-temp_real2)
              end if
           end do
           anam1=anam1*qm2_params%core_chg(IQM)*xyz_mm(4)*oneRIJ
           scale = scale + anam1
         end if
      endif

      ENUC = ENUC+scale !E(AB) = ZaZb(SASA,SBSB)(1+exp(-alphaA*RAB)+exp(-alphaB*RAB))

    endif  ! qtw - end OF if(qmmm_int == 3) 

    if (qmmm_nml%qmmm_switch) then
       if (rij <= qmmm_nml%r_switch_lo) then
         f_switch = one
       else if (rij < qmmm_nml%r_switch_hi) then 
         f_switch = qm2_switch_func(RIJ, qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi)
       else 
         f_switch = zero
       end if
       E1B=E1B*f_switch
       ENUC=ENUC*f_switch

       pot_tmp = (one - f_switch)*xyz_mm(4)*oneRij*AU_TO_EV*BOHRS_TO_A
                                                  !converts (AU/A) to (ev/Bohrs)
       qmmm_struct%switched_mmpot(iqm) = qmmm_struct%switched_mmpot(iqm) + pot_tmp

       ENUC = ENUC + qm2_params%core_chg(iqm)*pot_tmp
    end if

end subroutine qm2_rotate_qmmm_light

SUBROUTINE qm2_rotate_qmmm_heavy(qm2_params,qmmm_nml, qm2_rij_eqns, qm2_struct, qmmm_struct, loop_count,IQM,xyz_qm,xyz_mm,E1B,ENUC)

!For SP-ATOMs
!-----------------------------------
!Written by Ross Walker (TSRI, 2005)
!-----------------------------------

!-----------------------------------
!This routine calculates the QMelectron - MMcore interactions and the
!QMcore-MMcore interaction for a give QM-MM pair. Initially calculated in
!a local frame by qm2_repp_qmmm they are subsequently rotated into the
!molecular frame by this routine. The core charge for the MM atom is the
!RESP charge from the prmtop file in electron units.

!On Input:
! loop_count - Offset into certain arrays
!        iqm - Current qm atom in our loop from 1 to nquant
!     xyz_qm - Cartesian coordinates of the QM atom.
!     xyz_mm - Cartesian coordinates and charge (4) of MM atom.
!        E1B - QM-MM electron core interaction in KCal/mol for each
!              of the 10 multipoles of the QM atom. (s,s), (s,px), (s,py)
!              (s,pz), (px,px), (px,py), (px,pz), (py,py), (py,pz), (pz,pz)
!       ENUC - QM-MM core core interaction in KCal/mol.
!-----------------------------------

      use qmmm_module, only : qm2_structure, qm2_rij_eqns_structure, &
                              alph_mm, AXIS_TOL, EXPONENTIAL_CUTOFF
      use constants  , only : A_TO_BOHRS, A2_TO_BOHRS2, AU_TO_EV, HALF_AU_TO_EV, FOURTH_AU_TO_EV, &
                              BOHRS_TO_A, one, zero
      use qmmm_struct_module, only : qmmm_struct_type
   use qm2_params_module,  only : qm2_params_type
   use qmmm_nml_module   , only : qmmm_nml_type
                           
      implicit none
! Passed in
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qm2_rij_eqns_structure),intent(inout) :: qm2_rij_eqns
      type(qm2_structure),intent(inout) :: qm2_struct
      type(qmmm_struct_type), intent(inout) :: qmmm_struct
      integer, intent(in) :: loop_count, iqm
      _REAL_, intent(in) :: xyz_qm(3), xyz_mm(4)
      _REAL_, intent(out) :: e1b(10), enuc

! Local
      _REAL_ RI(4)
      _REAL_ :: oneBDDi1, oneBDDi2, oneBDDi3
      _REAL_ :: sqrtrr2aqe, qmi_DD, qmi_QQ
      _REAL_ :: X1, X2, X3,Y1, Y2, Y3,Z1, Z2,Z3, oneZ3, X3_2
      _REAL_ :: chrgmm, RR2, oneRIJ, RIJ, RR, RIJ2, CHGMM_RI2, CHGMM_RI3, CHGMM_RI4
      _REAL_ :: scale, anam1, temp_real, temp_real2
      _REAL_ :: sf1, sf2, rho_pm3mmx
      _REAL_ :: qm2_switch_func  ! function
      _REAL_ :: f_switch, pot_tmp
      integer :: i, qmitype
      
      X1 = xyz_qm(1) - xyz_mm(1)
      X2 = xyz_qm(2) - xyz_mm(2)
      X3 = xyz_qm(3) - xyz_mm(3)
      chrgmm = xyz_mm(4)
#include "qm2_array_locations.h"

! qtw - PM3/MM* - MODIFIED QM-MM INTERFACE
    !if ((qmmm_nml%qmmm_int == 3 .or. qmmm_nml%qmmm_int == 4) .AND. qmmm_nml%qmtheory%PM3) then  ! PM3/MM*
    if (qmmm_struct%PM3MMX_INTERFACE) then

      qmitype = qmmm_struct%qm_atom_type(iqm)

      sf1 = qm2_params%scale_factor1_pm3mmx(1,qmitype)
      sf2 = qm2_params%scale_factor2_pm3mmx(1,qmitype)
      ! rho_pm3mmx in qmmm_int==3 is always zero
      rho_pm3mmx = qm2_params%rho_pm3mmx(qmitype)

      RIJ2=X1*X1+X2*X2+X3*X3
      RR2=RIJ2*A2_TO_BOHRS2
      oneRIJ=1.0d0/sqrt(RIJ2)
      RIJ =RIJ2*oneRIJ
      RR = RIJ*A_TO_BOHRS
      !oneBDDi1=qm2_params%multip_2c_elec_params(3,iqm)**2
      oneBDDi1=(qm2_params%multip_2c_elec_params(3,iqm)+rho_pm3mmx)**2
      RI(1) = AU_TO_EV/sqrt(RR2+oneBDDi1)
      qmi_DD=qm2_params%multip_2c_elec_params(1,iqm)
      qmi_QQ=2.0d0*qm2_params%multip_2c_elec_params(2,iqm)
      !oneBDDi2=qm2_params%multip_2c_elec_params(4,iqm)**2
      !oneBDDi3=qm2_params%multip_2c_elec_params(5,iqm)**2
      oneBDDi2=(qm2_params%multip_2c_elec_params(4,iqm)+rho_pm3mmx)**2
      oneBDDi3=(qm2_params%multip_2c_elec_params(5,iqm)+rho_pm3mmx)**2
      RI(2)= HALF_AU_TO_EV*(1.0d0/SQRT((RR+qmi_DD)**2+oneBDDi2) - 1.0d0/SQRT((RR-qmi_DD)**2+oneBDDi2))

      SQRTRR2AQE=1.0d0/SQRT(RR2+oneBDDi3)
      RI(3)= RI(1) + FOURTH_AU_TO_EV*(1.0d0/SQRT((RR+qmi_QQ)**2+oneBDDi3) + 1.0d0/SQRT((RR-qmi_QQ)**2+oneBDDi3)) &
             - HALF_AU_TO_EV*SQRTRR2AQE

      RI(4) = RI(1) + HALF_AU_TO_EV*(1.0d0/SQRT(RR2+(qmi_QQ*qmi_QQ)+oneBDDi3) - SQRTRR2AQE)
      E1B(1)=-chrgmm*RI(1)
      ENUC = -qm2_params%core_chg(IQM)*E1B(1)
      ! Added ENUC correction
      scale = -exp(-RIJ*sf1) + exp(-RIJ*sf2)
      ENUC = ENUC + SIGN(1.0D0,chrgmm)*ENUC*scale

    else  ! qtw - then FOLLOW THE NORMAL QM/MM PROCEDURE

!---RIJ and related equation storage specifics---
      if (qmmm_nml%qmmmrij_incore) then
        scale=qm2_rij_eqns%qmmmrijdata(QMMMEXP1,loop_count)+qm2_rij_eqns%qmmmrijdata(QMMMEXP2,loop_count)
        RIJ = qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)
        oneRIJ=qm2_rij_eqns%qmmmrijdata(QMMMONERIJ,loop_count)
        RI(1) = AU_TO_EV*qm2_rij_eqns%qmmmrijdata(QMMMSQRTAEE,loop_count)
        RI(2) = HALF_AU_TO_EV*(qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRADDADE,loop_count)- &
                               qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMDDADE,loop_count))

        RI(3) = RI(1) + FOURTH_AU_TO_EV*(qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQAQE,loop_count)+ &
                qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRMQQAQE,loop_count))- &
                HALF_AU_TO_EV*qm2_rij_eqns%qmmmrijdata(QMMMSQRTRR2AQE,loop_count)

        RI(4) = RI(1) + HALF_AU_TO_EV*(qm2_rij_eqns%qmmmrijdata(QMMMSQRTRRAQQ2AQE,loop_count)- &
                                       qm2_rij_eqns%qmmmrijdata(QMMMSQRTRR2AQE,loop_count))
      else
        RIJ2=X1*X1+X2*X2+X3*X3
        RR2=RIJ2*A2_TO_BOHRS2
        oneRIJ=1.0d0/sqrt(RIJ2)
        RIJ =RIJ2*oneRIJ
        RR = RIJ*A_TO_BOHRS
        scale=EXP(-qm2_params%cc_exp_params(iqm)*RIJ)+EXP(-ALPH_MM*RIJ)
        oneBDDi1=qm2_params%multip_2c_elec_params(3,iqm)**2
        RI(1) = AU_TO_EV/sqrt(RR2+oneBDDi1)
        qmi_DD=qm2_params%multip_2c_elec_params(1,iqm)
        qmi_QQ=2.0d0*qm2_params%multip_2c_elec_params(2,iqm)
        oneBDDi2=qm2_params%multip_2c_elec_params(4,iqm)**2
        oneBDDi3=qm2_params%multip_2c_elec_params(5,iqm)**2
        RI(2)= HALF_AU_TO_EV*(1.0d0/SQRT((RR+qmi_DD)**2+oneBDDi2) - 1.0d0/SQRT((RR-qmi_DD)**2+oneBDDi2))

        SQRTRR2AQE=1.0d0/SQRT(RR2+oneBDDi3)
        RI(3)= RI(1) + FOURTH_AU_TO_EV*(1.0d0/SQRT((RR+qmi_QQ)**2+oneBDDi3) + 1.0d0/SQRT((RR-qmi_QQ)**2+oneBDDi3)) &
               - HALF_AU_TO_EV*SQRTRR2AQE

        RI(4) = RI(1) + HALF_AU_TO_EV*(1.0d0/SQRT(RR2+(qmi_QQ*qmi_QQ)+oneBDDi3) - SQRTRR2AQE)
      end if
!---end RIJ and related equation storage specifics---

! Step 1 is to calculate the interaction between the QM electrons and the MM
! atomic charge for this QM-MM pair in the local frame. MM charge is the RESP
! charge which is currently stored in the 4th index of the MM cartesian coordinate
! array. The reason the charge is packed into the 4th index of the coordinate array
! is so that we move linearly through memory = more cache hits. (RCW). - DONE ABOVE FOR SPEED

!All atoms have S orbitals
      E1B(1)=-chrgmm*RI(1)

!Calculate QM-QM nuclear term
!This is calculates as the QM core charge * chrgmm * RI(1)
!Reuse E1B term to save doing the multiplication of chrgmm*RI(1) again.
      ENUC = -qm2_params%core_chg(IQM)*E1B(1)
      scale = ABS(scale*ENUC)

!Add in the extra core core repulsion terms if qmmm_int==2. This is equivalent of the
!way CHARMM and DYNAMO do QMMM interactions.
      if(qmmm_nml%qmmm_int == 2 .and. (qmmm_struct%AM1_OR_PM3 .or. qmmm_nml%qmtheory%PM6) ) then
!        Add gaussians.
         if (qmmm_nml%qmmmrij_incore) then
           rij = qm2_rij_eqns%qmmmrijdata(QMMMRIJ,loop_count)
         end if
         qmitype = qmmm_struct%qm_atom_type(iqm)
         anam1=0.0d0
         do I=1,qm2_params%num_fn(qmitype)
            temp_real=RIJ-qm2_params%FN3(i,qmitype)
            temp_real2=qm2_params%FN2(i,qmitype)*temp_real*temp_real
            if (temp_real2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
              anam1=anam1+qm2_params%FN1(i,qmitype)*EXP(-temp_real2)
            end if
         end do
         anam1=anam1*qm2_params%core_chg(IQM)*xyz_mm(4)*oneRIJ
         scale = scale + anam1
      end if
      ENUC = ENUC + scale

    endif  ! qtw - end OF if(qmmm_int == 3) 

      X1 = X1*oneRIJ
      X2 = X2*oneRIJ
      X3 = X3*oneRIJ
      X3_2 = X3*X3

      if (abs(X3) > (1.0D0-AXIS_TOL)) then
         X3 = SIGN(1.0D0,X3)
         Y1 = 0.0D0
         Y2 = 1.0D0
         Y3 = 0.0D0
         Z1 = 1.0D0
         Z2 = 0.0D0
         Z3 = 0.0D0
      else
         Z3=SQRT(1.0D0-X3_2)
         oneZ3=1.0D0/Z3
         Y1=-oneZ3*X2*SIGN(1.D0,X1)
         Y2=ABS(oneZ3*X1)
         Y3=0.0D0
         Z1=-oneZ3*X1*X3
         Z2=-oneZ3*X2*X3
      endif

      CHGMM_RI2=-chrgmm*RI(2)
      CHGMM_RI3=-chrgmm*RI(3)
      CHGMM_RI4=-chrgmm*RI(4)

      E1B(2) = CHGMM_RI2*X1
      E1B(3) = CHGMM_RI3*X1*X1+CHGMM_RI4*((Y1*Y1)+(Z1*Z1))
      E1B(4) = CHGMM_RI2*X2
      E1B(5) = CHGMM_RI3*X2*X1+CHGMM_RI4*((Y2*Y1)+(Z2*Z1))
      E1B(6) = CHGMM_RI3*X2*X2+CHGMM_RI4*((Y2*Y2)+(Z2*Z2))
      E1B(7) = CHGMM_RI2*X3
      E1B(8) = CHGMM_RI3*X3*X1+CHGMM_RI4*Z3*Z1
      E1B(9) = CHGMM_RI3*X3*X2+CHGMM_RI4*Z3*Z2
      E1B(10)= CHGMM_RI3*X3_2+CHGMM_RI4*Z3*Z3

      if (qmmm_nml%qmmm_switch) then
         if (rij <= qmmm_nml%r_switch_lo) then
           f_switch = one
         else if (rij < qmmm_nml%r_switch_hi) then 
           f_switch = qm2_switch_func(RIJ, qmmm_nml%r_switch_lo, qmmm_nml%r_switch_hi)
         else
           f_switch = zero
         end if
         E1B(1:10)=E1B(1:10)*f_switch
         ENUC=ENUC*f_switch

         pot_tmp = (one - f_switch)*xyz_mm(4)*oneRij*AU_TO_EV*BOHRS_TO_A
                                                    !converts (AU/A) to (ev/Bohrs)
         qmmm_struct%switched_mmpot(iqm) = qmmm_struct%switched_mmpot(iqm) + pot_tmp

         ENUC = ENUC + qm2_params%core_chg(iqm)*pot_tmp
      end if

end subroutine qm2_rotate_qmmm_heavy


subroutine qm2_hcore_add_switched(qm2_params,qmmm_mpi, H, switched_mmpot)

!--------------------------------------------------------------------
! This routine adds the switched_mmpot to the diagonal elements 
! of the H matrix.
!--------------------------------------------------------------------

   use qmmm_module, only : qmmm_mpi_structure
   use qm2_params_module,  only : qm2_params_type

   implicit none

   ! Passed in
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
   _REAL_, intent(inout) :: H(*)
   _REAL_, intent(in) :: switched_mmpot(*)

   ! Local
   _REAL_ :: temp_pot
   integer :: i, ia, ib, i1, i2, k

   ! add the switched_mmpot to the diagonal elements of the H matrix
   !do i = 1, qmmm_struct%nquant_nlink
   do i=qmmm_mpi%nquant_nlink_start,qmmm_mpi%nquant_nlink_end
      ia = qm2_params%orb_loc(1,i)
      ib = qm2_params%orb_loc(2,i)
      temp_pot = switched_mmpot(i)
      do i1 = ia, ib
         i2 = qm2_params%pascal_tri2(i1)
         H(i2) = H(i2) - temp_pot
      end do
   end do

   return
end subroutine qm2_hcore_add_switched


function qm2_switch_func(r, r_lower, r_upper)

   implicit none

   _REAL_, intent(in) :: r, r_lower, r_upper

   _REAL_ :: qm2_switch_func, f_tmp

   f_tmp = (r_upper - r)/(r_upper - r_lower)

   qm2_switch_func =    10.0D0 * f_tmp**3 &
                    & - 15.0D0 * f_tmp**4 &
                    & +  6.0D0 * f_tmp**5

end function qm2_switch_func


function qm2_switch_derv(r, r_lower, r_upper)

   implicit none

   _REAL_, intent(in) :: r, r_lower, r_upper

   _REAL_ :: qm2_switch_derv, f_tmp, df_tmp

   f_tmp = (r_upper - r)/(r_upper - r_lower)
   df_tmp = -1.0D0 / (r_upper - r_lower)

   qm2_switch_derv =    30.0D0 * f_tmp**2 &
                    & - 60.0D0 * f_tmp**3 &
                    & + 30.0D0 * f_tmp**4

   qm2_switch_derv = qm2_switch_derv * df_tmp

end function qm2_switch_derv

