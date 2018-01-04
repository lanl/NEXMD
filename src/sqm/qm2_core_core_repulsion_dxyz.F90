! <compile=optimized>
#include "dprec.fh"
subroutine qm2_core_core_repulsion_dxyz(qmmm_nml,qm2_params,qmmm_struct, iat, jat, rij, onerij, xyzij, gij, dgij, dxyz)
! ----------------------------------------------------------------------
! PURPOSE: Calculate the derivative of the core core repulsion energy
!          between atoms iat, jat
!          and add it to the gradient dxyz
!
! E^{\rm core}_{i,j} = Z_i Z_j (s_i s_i|s_j s_j) [1 + f_{i,j}] + g_{i,j}
!
! MNDO   : f_{i,j} = e^{-\alpha_i R_{i,j}} + e^{-alpha_j R_{i,j}}
!          g_{i,j} = 0
!
! AM1/PM3: f_{i,j} = f^{\rm MNDO}_{i,j} 
!          g_{i,j} = (Z_i Z_j)/R_{i,j}
!                              \sum_k a_k^i e^{-b_k^i (R_{i,j}-c_k^i)^2}
!
! PM6    : f_{i,j} = x_{ij}e^{-alpha_{ij}(R_ij+0.0003R_{ij}^6)}
!          with special cases for atom pairs of the type
!          OH, NH, CH, CC, SiO
!          g_{i,j} = g^{\rm AM1}_{i,j}
!                      + 10^{-8}[(Z_i^{1/3}+Z_j^{1/3})/R_{ij}]^{12}
!
! PM3-MAIS: f_{i,j} = f^{\rm MNDO}_{i,j}
!           g_{i,j} = sum_k a_k(ij) e^{-b_K(ij) [(R_ij-c_k(ij)^2]}
!
! The internuclear distance and its inverse are passed down for 
! efficiency reasons. They have been calculated before for the
! evaluation of the ERI (s_i s_i|s_j s_j)
!
! Author: Andreas W. Goetz
!         <agoetz@sdsc.edu>
! Date  : April 2010
! ----------------------------------------------------------------------
  use qm2_params_module,  only : qm2_params_type
  use qmmm_nml_module   , only : qmmm_nml_type
  use qmmm_module, only :EXPONENTIAL_CUTOFF
  use constants, only : zero, one, two, ten, twelve, twenty, ten_to_minus8, third
  use qmmm_struct_module, only : qmmm_struct_type
  implicit none

  type(qm2_params_type), intent(in) :: qm2_params
  type(qmmm_nml_type), intent(in) :: qmmm_nml
  type(qmmm_struct_type), intent(in) :: qmmm_struct
  integer, intent(in)  :: iat, jat
  _REAL_,  intent(in)  :: rij, onerij, xyzij(3)
  _REAL_,  intent(in)  :: gij, dgij(3)
  _REAL_,  intent(out) :: dxyz(3)

  integer :: i, j
  integer :: iatyp, jatyp
  integer :: iatnum, jatnum, ijmin, ijmax
  _REAL_  :: alphai, alphaj, alpi, alpj, ccij
  _REAL_  :: anam1, scale, dscale(3), tmp, tmp2
  _REAL_  :: pddg_exp1, pddg_exp2, pddg_exp3, pddg_exp4, pddg_correction
  _REAL_  :: xyzonerij(3)
  ! PM6 specific stuff
  _REAL_  :: upc
  _REAL_, parameter :: threemin4 = 3.0d-04
  _REAL_, parameter :: opemin3   = 1.8d-3
  _REAL_, parameter :: ccfac = 9.28d0
  _REAL_, parameter :: ccexp = 5.98d0
  _REAL_, parameter :: siofac = 0.7d-03
  _REAL_, parameter :: sioexp = 2.9d0
  _REAL_ :: alpab, xab

  dxyz(1:3) = zero

  xyzonerij(1:3) = xyzij(1:3) * onerij

  iatyp  = qmmm_struct%qm_atom_type(iat) 
  jatyp  = qmmm_struct%qm_atom_type(jat) 
  iatnum = qmmm_struct%iqm_atomic_numbers(iat)
  jatnum = qmmm_struct%iqm_atomic_numbers(jat)
  ijmin = min(iatnum, jatnum)
  ijmax = max(iatnum, jatnum)
  
  ccij = qm2_params%core_chg(iat) * qm2_params%core_chg(jat)

  ! ------------------------------------------------
  ! Scaling factor and derivatives for monopole term
  ! ------------------------------------------------
  SWITCHMETHOD: if (qmmm_nml%qmtheory%PM6) then

     ! -------------------------------------
     ! PM6 correction for unpolarizable core
     ! -------------------------------------
     upc = onerij * (dble(iatnum)**third + dble(jatnum)**third)
     upc = -twelve * ten_to_minus8 * onerij * upc**12
     dxyz(1:3) = dxyz(1:3) + upc*xyzonerij(1:3)

     ! ---------------------------------
     ! PM6 scaling factor and derivative
     ! ---------------------------------
     alpab = qm2_params%pm6_alpab(iatyp, jatyp)
     xab   = qm2_params%pm6_xab  (iatyp, jatyp)

     if ( (ijmin == 1) .and. ( (ijmax == 6) .or. (ijmax == 7) .or. (ijmax == 8) ) ) then
        ! CH, NH and OH special case
        tmp = alpab*rij*rij
        tmp2 = two*xab*exp(-tmp)
        dscale(1:3) = -two*alpab*tmp2*xyzij(1:3)
        scale = one + tmp2
     else
        ! general case
        tmp = alpab*( rij + threemin4*(rij**6) )
        tmp2 = two*xab*exp(-tmp)
        tmp = one + opemin3*rij**5
        dscale(1:3) = -alpab*tmp*tmp2*xyzonerij(1:3)
        scale = one + tmp2
     end if

     if ( (iatnum == 6) .and. (jatnum == 6) ) then
        ! CC scaling factor correction, equation (10)
        tmp = ccexp * rij
        tmp2 = ccfac*exp(-tmp)
        scale = scale + tmp2
        dscale(1:3) = dscale(1:3) - ccexp*tmp2*xyzonerij(1:3)
     else if ( (ijmax == 14) .and. (ijmin == 8) ) then
        ! SiO scaling factor correction, equation (11)
        tmp = (rij - sioexp)
        tmp2 = -siofac*exp(-tmp**2)
        scale = scale + tmp2
        dscale(1:3) = dscale(1:3) + two*tmp*tmp2*xyzonerij(1:3)
     end if

  else SWITCHMETHOD

     ! -------------------
     ! MNDO scaling factor
     ! -------------------
     alphai = qm2_params%cc_exp_params(iat)
     alpi = exp( -alphai * rij )
     if (iatyp == jatyp) then
        alphaj = alphai
        alpj = alpi
     else
        alphaj = qm2_params%cc_exp_params(jat)
        alpj = exp( -alphaj * rij )
     end if

     if ( (ijmin == 1) .and. ( (ijmax == 7) .or. (ijmax == 8) ) ) then
        ! NH and OH atom pairs
        if (iatnum == 1) then
           scale = one + alpi + rij*alpj
           dscale(1:3) = -( alphai*alpi + (rij*alphaj-one)*alpj )*xyzonerij(1:3)
        else
           scale = one + rij*alpi + alpj
           dscale(1:3) = -( (rij*alphai-one)*alpi + alphaj*alpj)*xyzonerij(1:3)
        end if
     else
        ! default case
        scale = one + alpi + alpj
        dscale(1:3) = -(alphai*alpi + alphaj*alpj)*xyzonerij(1:3)
     end if

  end if SWITCHMETHOD

  ! ---------------------------------
  ! Monopole contribution to gradient
  ! ---------------------------------
  dxyz(1:3) = dxyz(1:3) + ccij*scale*dgij(1:3)
  dxyz(1:3) = dxyz(1:3) + ccij*gij*dscale(1:3)

  ! -------------------------------------
  ! AM1, PM3 and PM6 Gaussian corrections
  ! dxyz=-A*[1/(R*R)+2.D0*B*(R-C)/R)*EXP(-B*(R-C)**2]*ccij*xyzij/R
  ! -------------------------------------
  if (qmmm_struct%AM1_OR_PM3 .or. qmmm_nml%qmtheory%PM6) then
     anam1 = zero
     do i = 1, qm2_params%num_fn(iatyp)
        tmp = rij - qm2_params%FN3(i,iatyp)
        tmp2 = qm2_params%FN2(i,iatyp)*tmp*tmp
        if (tmp2 < EXPONENTIAL_CUTOFF) then ! Skip doing the exponential if it is essentially zero
           anam1 = anam1 + qm2_params%FN1(i,iatyp)* &
                (onerij*onerij + two*qm2_params%FN2(i,iatyp)*tmp*onerij)*exp(-tmp2)
        end if
     end do
     if (iatyp == jatyp) then
        anam1 = anam1 + anam1
     else
        do i = 1, qm2_params%num_fn(jatyp)
           tmp = rij - qm2_params%FN3(i,jatyp)
           tmp2 = qm2_params%FN2(i,jatyp)*tmp*tmp
           if (tmp < EXPONENTIAL_CUTOFF) then ! Skip doing the exponential if it is essentially zero
              anam1 = anam1 + qm2_params%FN1(i,jatyp)* &
                   (onerij*onerij + two*qm2_params%FN2(i,jatyp)*tmp*onerij)*exp(-tmp2)
           end if
        end do
     end if
     dxyz(1:3) = dxyz(1:3) - ccij*anam1*xyzonerij(1:3)
  end if

  ! ---------------------------------
  ! PM3-MAIS SPECIFIC DERIVATIVE CODE
  ! dxyz = -2.D0 * A * B * (R-C) * EXP(-B*(R-C)**2)
  ! ---------------------------------
  if (qmmm_nml%qmtheory%PM3MAIS) then
     anam1 = zero
     do i = 1, 3
        tmp = rij - qm2_params%pm3mais_gamab(iatyp,jatyp,i)
        tmp2 = qm2_params%pm3mais_betab(iatyp,jatyp,i) * tmp * tmp
        if (tmp2 < EXPONENTIAL_CUTOFF) then!Skip doing the exponential if it is essentially zero
           anam1 = anam1 + 2.d0*qm2_params%pm3mais_alpab(iatyp,jatyp,i) &
                *qm2_params%pm3mais_betab(iatyp,jatyp,i)*tmp*exp(-tmp2)
        endif
     enddo
     dxyz(1:3) = dxyz(1:3) - anam1*xyzonerij(1:3)
  end if

  ! ---------------
  ! PDDG correction
  ! ---------------
  if (qmmm_struct%PDDG_IN_USE) then

     pddg_exp1 = exp(-ten * (rij - qm2_params%pddge1(iat) - qm2_params%pddge1(jat))**2)
     pddg_exp2 = exp(-ten * (rij - qm2_params%pddge1(iat) - qm2_params%pddge2(jat))**2)
     pddg_exp3 = exp(-ten * (rij - qm2_params%pddge2(iat) - qm2_params%pddge1(jat))**2)
     pddg_exp4 = exp(-ten * (rij - qm2_params%pddge2(iat) - qm2_params%pddge2(jat))**2)
     
     pddg_exp1 = -twenty * (rij - qm2_params%pddge1(iat) - qm2_params%pddge1(jat)) * pddg_exp1
     pddg_exp2 = -twenty * (rij - qm2_params%pddge1(iat) - qm2_params%pddge2(jat)) * pddg_exp2
     pddg_exp3 = -twenty * (rij - qm2_params%pddge2(iat) - qm2_params%pddge1(jat)) * pddg_exp3
     pddg_exp4 = -twenty * (rij - qm2_params%pddge2(iat) - qm2_params%pddge2(jat)) * pddg_exp4

     pddg_correction = qm2_params%pddg_term1(iatyp,jatyp)*pddg_exp1 + &
          qm2_params%pddg_term2(iatyp,jatyp)*pddg_exp2 + &
          qm2_params%pddg_term3(iatyp,jatyp)*pddg_exp3 + &
          qm2_params%pddg_term4(iatyp,jatyp)*pddg_exp4

     dxyz(1:3) = dxyz(1:3) + pddg_correction*xyzonerij(1:3)

    end if
  
end subroutine qm2_core_core_repulsion_dxyz
