#include "dprec.fh"
#include "assert.fh"

!------------------------------------------------------------
module amoeba_self
implicit none
private

#  include "amoeba_mpole_index.h"

integer,save :: do_flag

  public AM_SELF_permfield,AM_SELF_dipole_field,AM_SELF_ene_torque, &
         AM_SELF_set_user_bit
#ifdef MPI
  public AM_SELF_bcast
#endif
contains
!-------------------------------------------------------
#ifdef MPI
subroutine AM_SELF_bcast
  implicit none
  integer ierr

  include 'mpif.h'
# include "parallel.h"
  call mpi_bcast(do_flag,1,MPI_INTEGER,0,commsander,ierr)
end subroutine AM_SELF_bcast
#endif

subroutine AM_SELF_set_user_bit(do_this)
  integer,intent(in) :: do_this
#include "do_flag.h"

  ! set the valid bit---this part always since no parmread needed
  do_flag = ibset(do_flag,VALID_BIT)

  if ( do_this == 1 )then ! do in all cases
    do_flag = ibset(do_flag,USER_INDUCE_BIT)
    do_flag = ibset(do_flag,USER_POSTINDUCE_BIT)
  elseif ( do_this == 2 )then ! do the induction, not the post-induction
    do_flag = ibset(do_flag,USER_INDUCE_BIT)
    do_flag = ibclr(do_flag,USER_POSTINDUCE_BIT)
  elseif ( do_this == 3 )then ! do the post-induction, not the induction
    do_flag = ibclr(do_flag,USER_INDUCE_BIT)
    do_flag = ibset(do_flag,USER_POSTINDUCE_BIT)
  elseif ( do_this == 0 )then 
    do_flag = ibclr(do_flag,USER_INDUCE_BIT)
    do_flag = ibclr(do_flag,USER_POSTINDUCE_BIT)
  else
    write(6,*)'AM_SELF_set_user_bit: bad value of user do_this'
    call mexit(6,1)
  endif
end subroutine AM_SELF_set_user_bit
!-------------------------------------------------------
subroutine AM_SELF_permfield(numatoms,direct_field,polar_field)
   use amoeba_multipoles, only : global_multipole
   use constants, only : three,four,pi
   integer,intent(in) :: numatoms
   _REAL_,intent(inout) :: direct_field(3,*),polar_field(3,*)

#  include "ew_pme_recip.h"
#include "do_flag.h"
   _REAL_ :: factor
   integer n

   if ( iand(do_flag,PROCEED_INDUCE) /= PROCEED_INDUCE )return

   factor = four * ew_coeff**3 / (three*sqrt(PI))
   do n = 1,numatoms
      direct_field(1,n) = direct_field(1,n)-factor*global_multipole(Ind_100,n)
      direct_field(2,n) = direct_field(2,n)-factor*global_multipole(Ind_010,n)
      direct_field(3,n) = direct_field(3,n)-factor*global_multipole(Ind_001,n)
      polar_field(1,n) = polar_field(1,n)-factor*global_multipole(Ind_100,n)
      polar_field(2,n) = polar_field(2,n)-factor*global_multipole(Ind_010,n)
      polar_field(3,n) = polar_field(3,n)-factor*global_multipole(Ind_001,n)
   enddo
end subroutine AM_SELF_permfield
!-------------------------------------------------------
subroutine AM_SELF_dipole_field(numatoms, &
                  ind_dip1,ind_dip2,dip_field1,dip_field2)
   use constants, only : three,four,pi
   integer,intent(in) :: numatoms
   _REAL_,intent(in) :: ind_dip1(3,*),ind_dip2(3,*)
   _REAL_,intent(inout) :: dip_field1(3,*),dip_field2(3,*)

#  include "ew_pme_recip.h"
#include "do_flag.h"
   _REAL_ :: factor
   integer n

   if ( iand(do_flag,PROCEED_INDUCE) /= PROCEED_INDUCE )return

   factor = four * ew_coeff**3 / (three*sqrt(PI))
   do n = 1,numatoms
      dip_field1(1,n) = dip_field1(1,n) - factor*ind_dip1(1,n)
      dip_field1(2,n) = dip_field1(2,n) - factor*ind_dip1(2,n)
      dip_field1(3,n) = dip_field1(3,n) - factor*ind_dip1(3,n)
      dip_field2(1,n) = dip_field2(1,n) - factor*ind_dip2(1,n)
      dip_field2(2,n) = dip_field2(2,n) - factor*ind_dip2(2,n)
      dip_field2(3,n) = dip_field2(3,n) - factor*ind_dip2(3,n)
   enddo
end subroutine AM_SELF_dipole_field
!-------------------------------------------------------
subroutine AM_SELF_ene_torque(numatoms,  &
                      ind_dip_d,ind_dip_p,ene_perm,ene_ind)
   use amoeba_multipoles, only : global_multipole, &
                                 coulomb_const_kcal_per_mole,torque_field
   use constants, only : zero,half,two,pi
   integer,intent(in) :: numatoms
   _REAL_,intent(in) :: ind_dip_d(3,*),ind_dip_p(3,*)
   _REAL_,intent(inout) :: ene_perm,ene_ind

   _REAL_ :: delx,dely,delz,B(0:4),fac,fact,gmi(10),phi(10),i_di(3), &
             i_mi(3),e_pp,e_ind
   _REAL_  :: Rn(1),Rn_1(4),Rn_2(10),Rn_3(20),Rn_4(35)
   integer i,j,n
#include "do_flag.h"
#  include "ew_pme_recip.h"

   ene_perm = zero
   ene_ind = zero

   if ( iand(do_flag,PROCEED_POSTINDUCE) /= PROCEED_POSTINDUCE )return

   fact = two*ew_coeff / sqrt(PI)
   fac = -two*ew_coeff*ew_coeff
   do j = 0,4
     B(j) = fact/(two*j+1)
     fact = fac*fact
   enddo
   delx = zero
   dely = zero
   delz = zero

   n = 4
   Rn(Ind_000) = B(n)
   Rn_1(Ind_000) = B(n-1)
   Rn_1(Ind_100) = delx*Rn(Ind_000)
   Rn_1(Ind_010) = dely*Rn(Ind_000)
   Rn_1(Ind_001) = delz*Rn(Ind_000)
   Rn_2(Ind_000) = B(n-2)
   Rn_2(Ind_100) = delx*Rn_1(Ind_000)
   Rn_2(Ind_010) = dely*Rn_1(Ind_000)
   Rn_2(Ind_001) = delz*Rn_1(Ind_000)
   Rn_2(Ind_200) = Rn_1(Ind_000) + delx*Rn_1(Ind_100)
   Rn_2(Ind_020) = Rn_1(Ind_000) + dely*Rn_1(Ind_010)
   Rn_2(Ind_002) = Rn_1(Ind_000) + delz*Rn_1(Ind_001)
   Rn_2(Ind_110) = delx*Rn_1(Ind_010)
   Rn_2(Ind_101) = delx*Rn_1(Ind_001)
   Rn_2(Ind_011) = dely*Rn_1(Ind_001)
   Rn_3(Ind_000) = B(n-3) 
   Rn_3(Ind_100) = delx*Rn_2(Ind_000)
   Rn_3(Ind_010) = dely*Rn_2(Ind_000)
   Rn_3(Ind_001) = delz*Rn_2(Ind_000)
   Rn_3(Ind_200) = Rn_2(Ind_000) + delx*Rn_2(Ind_100)
   Rn_3(Ind_020) = Rn_2(Ind_000) + dely*Rn_2(Ind_010)
   Rn_3(Ind_002) = Rn_2(Ind_000) + delz*Rn_2(Ind_001)
   Rn_3(Ind_110) = delx*Rn_2(Ind_010)
   Rn_3(Ind_101) = delx*Rn_2(Ind_001)
   Rn_3(Ind_011) = dely*Rn_2(Ind_001)
   Rn_3(Ind_300) = 2.d0*Rn_2(Ind_100) + delx*Rn_2(Ind_200)
   Rn_3(Ind_030) = 2.d0*Rn_2(Ind_010) + dely*Rn_2(Ind_020)
   Rn_3(Ind_003) = 2.d0*Rn_2(Ind_001) + delz*Rn_2(Ind_002)
   Rn_3(Ind_210) = dely*Rn_2(Ind_200)
   Rn_3(Ind_201) = delz*Rn_2(Ind_200)
   Rn_3(Ind_120) = delx*Rn_2(Ind_020)
   Rn_3(Ind_021) = delz*Rn_2(Ind_020)
   Rn_3(Ind_102) = delx*Rn_2(Ind_002)
   Rn_3(Ind_012) = dely*Rn_2(Ind_002)
   Rn_3(Ind_111) = delx*Rn_2(Ind_011)
   Rn_4(Ind_000) = B(n-4) 
   Rn_4(Ind_100) = delx*Rn_3(Ind_000)
   Rn_4(Ind_010) = dely*Rn_3(Ind_000)
   Rn_4(Ind_001) = delz*Rn_3(Ind_000)
   Rn_4(Ind_200) = Rn_3(Ind_000) + delx*Rn_3(Ind_100)
   Rn_4(Ind_020) = Rn_3(Ind_000) + dely*Rn_3(Ind_010)
   Rn_4(Ind_002) = Rn_3(Ind_000) + delz*Rn_3(Ind_001)
   Rn_4(Ind_110) = delx*Rn_3(Ind_010)
   Rn_4(Ind_101) = delx*Rn_3(Ind_001)
   Rn_4(Ind_011) = dely*Rn_3(Ind_001)
   Rn_4(Ind_300) = 2.d0*Rn_3(Ind_100) + delx*Rn_3(Ind_200)
   Rn_4(Ind_030) = 2.d0*Rn_3(Ind_010) + dely*Rn_3(Ind_020)
   Rn_4(Ind_003) = 2.d0*Rn_3(Ind_001) + delz*Rn_3(Ind_002)
   Rn_4(Ind_210) = dely*Rn_3(Ind_200)
   Rn_4(Ind_201) = delz*Rn_3(Ind_200)
   Rn_4(Ind_120) = delx*Rn_3(Ind_020)
   Rn_4(Ind_021) = delz*Rn_3(Ind_020)
   Rn_4(Ind_102) = delx*Rn_3(Ind_002)
   Rn_4(Ind_012) = dely*Rn_3(Ind_002)
   Rn_4(Ind_111) = delx*Rn_3(Ind_011)
   Rn_4(Ind_400) = 3.d0*Rn_3(Ind_200) + delx*Rn_3(Ind_300)
   Rn_4(Ind_040) = 3.d0*Rn_3(Ind_020) + dely*Rn_3(Ind_030)
   Rn_4(Ind_004) = 3.d0*Rn_3(Ind_002) + delz*Rn_3(Ind_003)
   Rn_4(Ind_310) = dely*Rn_3(Ind_300)
   Rn_4(Ind_301) = delz*Rn_3(Ind_300)
   Rn_4(Ind_130) = delx*Rn_3(Ind_030)
   Rn_4(Ind_031) = delz*Rn_3(Ind_030)
   Rn_4(Ind_103) = delx*Rn_3(Ind_003)
   Rn_4(Ind_013) = dely*Rn_3(Ind_003)
   Rn_4(Ind_220) = Rn_3(Ind_020) + delx*Rn_3(Ind_120)
   Rn_4(Ind_202) = Rn_3(Ind_002) + delx*Rn_3(Ind_102)
   Rn_4(Ind_022) = Rn_3(Ind_002) + dely*Rn_3(Ind_012)
   Rn_4(Ind_211) = dely*Rn_3(Ind_201)
   Rn_4(Ind_121) = delx*Rn_3(Ind_021)
   Rn_4(Ind_112) = delx*Rn_3(Ind_012)
   do i = 1,numatoms
      do j = 1,10
         gmi(j) = global_multipole(j,i)
      enddo
      do j = 1,3
         i_di(j) = ind_dip_d(j,i)
         i_mi(j) = ind_dip_d(j,i) + ind_dip_p(j,i)
      enddo
      ! self-field due to permanent mpoles at i and derivs wrt r_j-r_i (at 0)
      phi(Ind_000)=Rn_4(Ind_000)*gmi(Ind_000)+Rn_4(Ind_100)*gmi(Ind_100)+ &
                   Rn_4(Ind_010)*gmi(Ind_010)+Rn_4(Ind_001)*gmi(Ind_001)+ &
                   Rn_4(Ind_200)*gmi(Ind_200)+Rn_4(Ind_020)*gmi(Ind_020)+ &
                   Rn_4(Ind_002)*gmi(Ind_002)+Rn_4(Ind_110)*gmi(Ind_110)+ &
                   Rn_4(Ind_101)*gmi(Ind_101)+Rn_4(Ind_011)*gmi(Ind_011)
      phi(Ind_100)=-(Rn_4(Ind_100)*gmi(Ind_000)+Rn_4(Ind_200)*gmi(Ind_100)+ &
                     Rn_4(Ind_110)*gmi(Ind_010)+Rn_4(Ind_101)*gmi(Ind_001)+ &
                     Rn_4(Ind_300)*gmi(Ind_200)+Rn_4(Ind_120)*gmi(Ind_020)+ &
                     Rn_4(Ind_102)*gmi(Ind_002)+Rn_4(Ind_210)*gmi(Ind_110)+ &
                     Rn_4(Ind_201)*gmi(Ind_101)+Rn_4(Ind_111)*gmi(Ind_011))
      phi(Ind_010)=-(Rn_4(Ind_010)*gmi(Ind_000)+Rn_4(Ind_110)*gmi(Ind_100)+ &
                     Rn_4(Ind_020)*gmi(Ind_010)+Rn_4(Ind_011)*gmi(Ind_001)+ &
                     Rn_4(Ind_210)*gmi(Ind_200)+Rn_4(Ind_030)*gmi(Ind_020)+ &
                     Rn_4(Ind_012)*gmi(Ind_002)+Rn_4(Ind_120)*gmi(Ind_110)+ &
                     Rn_4(Ind_111)*gmi(Ind_101)+Rn_4(Ind_021)*gmi(Ind_011))
      phi(Ind_001)=-(Rn_4(Ind_001)*gmi(Ind_000)+Rn_4(Ind_101)*gmi(Ind_100)+ &
                     Rn_4(Ind_011)*gmi(Ind_010)+Rn_4(Ind_002)*gmi(Ind_001)+ &
                     Rn_4(Ind_201)*gmi(Ind_200)+Rn_4(Ind_021)*gmi(Ind_020)+ &
                     Rn_4(Ind_003)*gmi(Ind_002)+Rn_4(Ind_111)*gmi(Ind_110)+ &
                     Rn_4(Ind_102)*gmi(Ind_101)+Rn_4(Ind_012)*gmi(Ind_011))
      phi(Ind_200)=Rn_4(Ind_200)*gmi(Ind_000)+Rn_4(Ind_300)*gmi(Ind_100)+ &
                   Rn_4(Ind_210)*gmi(Ind_010)+Rn_4(Ind_201)*gmi(Ind_001)+ &
                   Rn_4(Ind_400)*gmi(Ind_200)+Rn_4(Ind_220)*gmi(Ind_020)+ &
                   Rn_4(Ind_202)*gmi(Ind_002)+Rn_4(Ind_310)*gmi(Ind_110)+ &
                   Rn_4(Ind_301)*gmi(Ind_101)+Rn_4(Ind_211)*gmi(Ind_011)
      phi(Ind_020)=Rn_4(Ind_020)*gmi(Ind_000)+Rn_4(Ind_120)*gmi(Ind_100)+ &
                   Rn_4(Ind_030)*gmi(Ind_010)+Rn_4(Ind_021)*gmi(Ind_001)+ &
                   Rn_4(Ind_220)*gmi(Ind_200)+Rn_4(Ind_040)*gmi(Ind_020)+ &
                   Rn_4(Ind_022)*gmi(Ind_002)+Rn_4(Ind_130)*gmi(Ind_110)+ &
                   Rn_4(Ind_121)*gmi(Ind_101)+Rn_4(Ind_031)*gmi(Ind_011)
      phi(Ind_002)=Rn_4(Ind_002)*gmi(Ind_000)+Rn_4(Ind_102)*gmi(Ind_100)+ &
                   Rn_4(Ind_012)*gmi(Ind_010)+Rn_4(Ind_003)*gmi(Ind_001)+ &
                   Rn_4(Ind_202)*gmi(Ind_200)+Rn_4(Ind_022)*gmi(Ind_020)+ &
                   Rn_4(Ind_004)*gmi(Ind_002)+Rn_4(Ind_112)*gmi(Ind_110)+ &
                   Rn_4(Ind_103)*gmi(Ind_101)+Rn_4(Ind_013)*gmi(Ind_011)
      phi(Ind_110)=Rn_4(Ind_110)*gmi(Ind_000)+Rn_4(Ind_210)*gmi(Ind_100)+ &
                   Rn_4(Ind_120)*gmi(Ind_010)+Rn_4(Ind_111)*gmi(Ind_001)+ &
                   Rn_4(Ind_310)*gmi(Ind_200)+Rn_4(Ind_130)*gmi(Ind_020)+ &
                   Rn_4(Ind_112)*gmi(Ind_002)+Rn_4(Ind_220)*gmi(Ind_110)+ &
                   Rn_4(Ind_211)*gmi(Ind_101)+Rn_4(Ind_121)*gmi(Ind_011)
      phi(Ind_101)=Rn_4(Ind_101)*gmi(Ind_000)+Rn_4(Ind_201)*gmi(Ind_100)+ &
                   Rn_4(Ind_111)*gmi(Ind_010)+Rn_4(Ind_102)*gmi(Ind_001)+ &
                   Rn_4(Ind_301)*gmi(Ind_200)+Rn_4(Ind_121)*gmi(Ind_020)+ &
                   Rn_4(Ind_103)*gmi(Ind_002)+Rn_4(Ind_211)*gmi(Ind_110)+ &
                   Rn_4(Ind_202)*gmi(Ind_101)+Rn_4(Ind_112)*gmi(Ind_011)
      phi(Ind_011)=Rn_4(Ind_011)*gmi(Ind_000)+Rn_4(Ind_111)*gmi(Ind_100)+ &
                   Rn_4(Ind_021)*gmi(Ind_010)+Rn_4(Ind_012)*gmi(Ind_001)+ &
                   Rn_4(Ind_211)*gmi(Ind_200)+Rn_4(Ind_031)*gmi(Ind_020)+ &
                   Rn_4(Ind_013)*gmi(Ind_002)+Rn_4(Ind_121)*gmi(Ind_110)+ &
                   Rn_4(Ind_112)*gmi(Ind_101)+Rn_4(Ind_022)*gmi(Ind_011)
      e_pp = phi(Ind_000)*gmi(Ind_000) + phi(Ind_100)*gmi(Ind_100) + &
             phi(Ind_010)*gmi(Ind_010) + phi(Ind_001)*gmi(Ind_001) + &
             phi(Ind_200)*gmi(Ind_200) + phi(Ind_020)*gmi(Ind_020) + &
             phi(Ind_002)*gmi(Ind_002) + phi(Ind_110)*gmi(Ind_110) + &
             phi(Ind_101)*gmi(Ind_101) + phi(Ind_011)*gmi(Ind_011)
      e_ind = phi(Ind_100)*i_di(1) + phi(Ind_010)*i_di(2) + &
              phi(Ind_001)*i_di(3)
      ene_perm = ene_perm - e_pp ! subtract self-term
      ene_ind = ene_ind - e_ind ! subtract self-term
      ! torque due to permanent mpoles (should be zero)
      do j = 1,10
         torque_field(j,i) = torque_field(j,i) - phi(j)
      enddo
      ! field due to induced dipoles
      phi(Ind_000) = Rn_4(Ind_100)*i_mi(1) + Rn_4(Ind_010)*i_mi(2) + &
                     Rn_4(Ind_001)*i_mi(3)
      phi(Ind_100) = -(Rn_4(Ind_200)*i_mi(1) + Rn_4(Ind_110)*i_mi(2) + &
                       Rn_4(Ind_101)*i_mi(3))
      phi(Ind_010) = -(Rn_4(Ind_110)*i_mi(1) + Rn_4(Ind_020)*i_mi(2) + &
                       Rn_4(Ind_011)*i_mi(3))
      phi(Ind_001) = -(Rn_4(Ind_101)*i_mi(1) + Rn_4(Ind_011)*i_mi(2) + &
                       Rn_4(Ind_002)*i_mi(3))
      phi(Ind_200) = Rn_4(Ind_300)*i_mi(1) + Rn_4(Ind_210)*i_mi(2) + &
                     Rn_4(Ind_201)*i_mi(3)
      phi(Ind_020) = Rn_4(Ind_120)*i_mi(1) + Rn_4(Ind_030)*i_mi(2) + &
                     Rn_4(Ind_021)*i_mi(3)
      phi(Ind_002) = Rn_4(Ind_102)*i_mi(1) + Rn_4(Ind_012)*i_mi(2) + &
                     Rn_4(Ind_003)*i_mi(3)
      phi(Ind_110) = Rn_4(Ind_210)*i_mi(1) + Rn_4(Ind_120)*i_mi(2) + &
                     Rn_4(Ind_111)*i_mi(3)
      phi(Ind_101) = Rn_4(Ind_201)*i_mi(1) + Rn_4(Ind_111)*i_mi(2) + &
                     Rn_4(Ind_102)*i_mi(3)
      phi(Ind_011) = Rn_4(Ind_111)*i_mi(1) + Rn_4(Ind_021)*i_mi(2) + &
                     Rn_4(Ind_012)*i_mi(3)
      do j = 1,10
         torque_field(j,i) = torque_field(j,i) - half*phi(j)
      enddo
   enddo
   ene_perm = half*coulomb_const_kcal_per_mole*ene_perm
   ene_ind = half*coulomb_const_kcal_per_mole*ene_ind

end subroutine AM_SELF_ene_torque
!-------------------------------------------------------
end module amoeba_self
