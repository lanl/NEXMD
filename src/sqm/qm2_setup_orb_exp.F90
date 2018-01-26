! <compile=optimized>

! Add d-orbitals, correct coefficients for n=4,5; by Taisung Lee (Rutgers, 2011) 
!
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"
subroutine qm2_setup_orb_exp(qm2_params,qmmm_struct)

      use constants, only : A_TO_BOHRS, A2_TO_BOHRS2, A3_TO_BOHRS3
      use ElementOrbitalIndex, only: MaxAngularQuantumNumber, MaxGaussianExpansion
      use qmmm_struct_module, only : qmmm_struct_type
      use qm2_params_module,  only : qm2_params_type

      implicit none
      type(qmmm_struct_type), intent(inout) :: qmmm_struct
      type(qm2_params_type),intent(inout) :: qm2_params

!Local
      _REAL_:: ALLC(MaxGaussianExpansion,MaxGaussianExpansion,MaxAngularQuantumNumber+1)
      _REAL_:: ALLZ(MaxGaussianExpansion,MaxGaussianExpansion,MaxAngularQuantumNumber+1)
      integer:: ni,nqn
      integer:: i,j,k,l,ier
      integer::ng, np
      _REAL_:: xi(1:MaxAngularQuantumNumber+1)

!For pre-computing the overlap equations used in energy and derivative calculation.
      
      _REAL_, pointer, dimension(:,:) :: atom_orb_cc_s_by_type
      _REAL_, pointer, dimension(:,:) :: atom_orb_cc_p_by_type      
      _REAL_, pointer, dimension(:,:) :: atom_orb_cc_d_by_type
      _REAL_, pointer, dimension(:,:) :: atom_orb_zz_s_by_type
      _REAL_, pointer, dimension(:,:) :: atom_orb_zz_p_by_type
      _REAL_, pointer, dimension(:,:) :: atom_orb_zz_d_by_type
            
      _REAL_:: atom_orb_cc_s_x_s, atom_orb_zz_s_x_s
      _REAL_:: atom_orb_cc_s_x_p, atom_orb_zz_s_x_p
      _REAL_:: atom_orb_cc_s_x_d, atom_orb_zz_s_x_d
      _REAL_:: atom_orb_cc_p_x_p, atom_orb_zz_p_x_p 
      _REAL_:: atom_orb_cc_p_x_d, atom_orb_zz_p_x_d
      _REAL_:: atom_orb_cc_d_x_d, atom_orb_zz_d_x_d               
      _REAL_:: atom_orb_zz_one_s_a_s, atom_orb_zz_one_s_a_p, atom_orb_zz_one_s_a_d
      _REAL_:: atom_orb_zz_one_p_a_p, atom_orb_zz_one_p_a_d, atom_orb_zz_one_d_a_d     
      _REAL_:: atom_orb_sp_eqn, atom_orb_sd_eqn, atom_orb_pp_eqn, atom_orb_pd_eqn, atom_orb_dd_eqn    

      _REAL_ :: temp_real=0.0

!This routine fills atom_orb_cc_s, atom_orb_cc_p, atom_orb_zz_s, atom_orb_zz_p
!with the STO-6G orbital expansion data. It then de-allocates the no longer
!needed qm2_params%s_orb_exp and p_orb_exp. For this reason it should be called
!once and only once per run.

!  The coefficients are from 
!Calculate Overlap Integrals using a Gaussian Expansion
!STO-6G BY R.F. STEWART, J. CHEM. PHYS., 52 431-438, 1970

!     SET-UP THE STEWART'S STO-6G EXPANSIONS
!                                      1S
     ALLZ(1,1,1) =2.310303149D01
     ALLZ(2,1,1) =4.235915534D00
     ALLZ(3,1,1) =1.185056519D00
     ALLZ(4,1,1) =4.070988982D-01
     ALLZ(5,1,1) =1.580884151D-01
     ALLZ(6,1,1) =6.510953954D-02
 
     ALLC(1,1,1) =9.163596280D-03
     ALLC(2,1,1) =4.936149294D-02
     ALLC(3,1,1) =1.685383049D-01
     ALLC(4,1,1) =3.705627997D-01
     ALLC(5,1,1) =4.164915298D-01
     ALLC(6,1,1) =1.303340841D-01
!                                      2S 
     ALLZ(1,2,1) =2.768496241D01
     ALLZ(2,2,1) =5.077140627D00
     ALLZ(3,2,1) =1.426786050D00
     ALLZ(4,2,1) =2.040335729D-01
     ALLZ(5,2,1) =9.260298399D-02
     ALLZ(6,2,1) =4.416183978D-02
 
     ALLC(1,2,1) =-4.151277819D-03
     ALLC(2,2,1) =-2.067024148D-02
     ALLC(3,2,1) =-5.150303337D-02
     ALLC(4,2,1) =3.346271174D-01
     ALLC(5,2,1) =5.621061301D-01
     ALLC(6,2,1) =1.712994697D-01
!                                     2P
     ALLZ(1,2,2) =5.868285913D00
     ALLZ(2,2,2) =1.530329631D00
     ALLZ(3,2,2) =5.475665231D-01
     ALLZ(4,2,2) =2.288932733D-01
     ALLZ(5,2,2) =1.046655969D-01
     ALLZ(6,2,2) =4.948220127D-02
 
     ALLC(1,2,2) =7.924233646D-03
     ALLC(2,2,2) =5.144104825D-02
     ALLC(3,2,2) =1.898400060D-01
     ALLC(4,2,2) =4.049863191D-01
     ALLC(5,2,2) =4.012362861D-01
     ALLC(6,2,2) =1.051855189D-01
!                                      3S
     ALLZ(1,3,1) =3.273031938D00
     ALLZ(2,3,1) =9.200611311D-01
     ALLZ(3,3,1) =3.593349765D-01
     ALLZ(4,3,1) =8.636686991D-02
     ALLZ(5,3,1) =4.797373812D-02
     ALLZ(6,3,1) =2.724741144D-02
 
     ALLC(1,3,1) =-6.775596947D-03
     ALLC(2,3,1) =-5.639325779D-02
     ALLC(3,3,1) =-1.587856086D-01
     ALLC(4,3,1) =5.534527651D-01
     ALLC(5,3,1) =5.015351020D-01
     ALLC(6,3,1) =7.223633674D-02
!                                     3P
     ALLZ(1,3,2) =5.077973607D00
     ALLZ(2,3,2) =1.340786940D00
     ALLZ(3,3,2) =2.248434849D-01
     ALLZ(4,3,2) =1.131741848D-01
     ALLZ(5,3,2) =6.076408893D-02
     ALLZ(6,3,2) =3.315424265D-02
 
     ALLC(1,3,2) =-3.329929840D-03
     ALLC(2,3,2) =-1.419488340D-02
     ALLC(3,3,2) =1.639395770D-01
     ALLC(4,3,2) =4.485358256D-01
     ALLC(5,3,2) =3.908813050D-01
     ALLC(6,3,2) =7.411456232D-02
!                                     3D
     ALLZ(1,3,3) =2.488296923D00
     ALLZ(2,3,3) =7.981487853D-01
     ALLZ(3,3,3) =3.311327490D-01
     ALLZ(4,3,3) =1.559114463D-01
     ALLZ(5,3,3) =7.817734732D-02
     ALLZ(6,3,3) =4.058484363D-02
 
     ALLC(1,3,3) =7.283828112D-03
     ALLC(2,3,3) =5.386799363D-02
     ALLC(3,3,3) =2.072139149D-01
     ALLC(4,3,3) =4.266269092D-01
     ALLC(5,3,3) =3.843100204D-01
     ALLC(6,3,3) =8.902827546D-02     
!                                     4S
     ALLZ(1,4,1) = 3.232838646D+00
     ALLZ(2,4,1) = 3.605788802D-01
     ALLZ(3,4,1) = 1.717905487D-01
     ALLZ(4,4,1) = 5.277666487D-02
     ALLZ(5,4,1) = 3.163400284D-02
     ALLZ(6,4,1) = 1.874093091D-02
 
     ALLC(1,4,1) = 1.374817488D-03
     ALLC(2,4,1) =-8.666390043D-02
     ALLC(3,4,1) =-3.130627309D-01
     ALLC(4,4,1) = 7.812787397D-01
     ALLC(5,4,1) = 4.389247988D-01
     ALLC(6,4,1) = 2.487178756D-02
!                                     4P
     ALLZ(1,4,2) = 2.389722618D+00
     ALLZ(2,4,2) = 7.960947826D-01
     ALLZ(3,4,2) = 3.415541380D-01
     ALLZ(4,4,2) = 8.847434525D-02
     ALLZ(5,4,2) = 4.958248334D-02
     ALLZ(6,4,2) = 2.816929784D-02
 
     ALLC(1,4,2) =-1.665913575D-03
     ALLC(2,4,2) =-1.657464971D-02
     ALLC(3,4,2) =-5.958513378D-02
     ALLC(4,4,2) = 4.053115554D-01
     ALLC(5,4,2) = 5.433958189D-01
     ALLC(6,4,2) = 1.204970491D-01

!                                     4D
     ALLZ(1,4,3) = 4.634239420D+00
     ALLZ(2,4,3) = 1.341648295D+00
     ALLZ(3,4,3) = 2.209593028D-01
     ALLZ(4,4,3) = 1.101467943D-01
     ALLZ(5,4,3) = 5.904190370D-02
     ALLZ(6,4,3) = 3.232628887D-02  
 
     ALLC(1,4,3) =-4.749842876D-04
     ALLC(2,4,3) =-3.566777891D-03
     ALLC(3,4,3) = 1.111867048D-01
     ALLC(4,4,3) = 4.159646930D-01
     ALLC(5,4,3) = 4.621672517D-01
     ALLC(6,4,3) = 1.081250196D-01
 
!                                     5S
     ALLZ(1,5,1) = 1.410128298D+00
     ALLZ(2,5,1) = 5.077878915D-01
     ALLZ(3,5,1) = 1.847926858D-01
     ALLZ(4,5,1) = 1.061070594D-01
     ALLZ(5,5,1) = 3.669584901D-02
     ALLZ(6,5,1) = 2.213558430D-02
 
     ALLC(1,5,1) = 2.695439582D-03
     ALLC(2,5,1) = 1.850157487D-02
     ALLC(3,5,1) =-9.588628125D-02
     ALLC(4,5,1) =-5.200673560D-01
     ALLC(5,5,1) = 1.087619490D+00
     ALLC(6,5,1) = 3.103964343D-01
!                                      5P
     ALLZ(1,5,2) = 3.778623374D+00
     ALLZ(2,5,2) = 3.499121109D-01
     ALLZ(3,5,2) = 1.683175469D-01
     ALLZ(4,5,2) = 5.404070736D-02
     ALLZ(5,5,2) = 3.328911801D-02
     ALLZ(6,5,2) = 2.063815019D-02
 
     ALLC(1,5,2) = 1.163246387D-04
     ALLC(2,5,2) =-2.920771322D-02
     ALLC(3,5,2) =-1.381051233D-01
     ALLC(4,5,2) = 5.706134877D-01
     ALLC(5,5,2) = 4.768808140D-01
     ALLC(6,5,2) = 6.021665516D-02
!                                      5D
     ALLZ(1,5,3) = 8.820520428D-01
     ALLZ(2,5,3) = 3.410838409D-01
     ALLZ(3,5,3) = 9.204308840D-02
     ALLZ(4,5,3) = 5.472831774D-02
     ALLZ(5,5,3) = 3.391202830D-02
     ALLZ(6,5,3) = 2.108227374D-02
 
     ALLC(1,5,3) =-4.097311019D-03
     ALLC(2,5,3) =-2.508271857D-02
     ALLC(3,5,3) = 2.648458555D-01
     ALLC(4,5,3) = 5.097437054D-01
     ALLC(5,5,3) = 2.654483467D-01
     ALLC(6,5,3) = 2.623132212D-02     
!                                      6S
     ALLZ(1,6,1) = 5.800292686D-01
     ALLZ(2,6,1) = 2.718262251D-01
     ALLZ(3,6,1) = 7.938523262D-02
     ALLZ(4,6,1) = 4.975088254D-02
     ALLZ(5,6,1) = 2.983643556D-02
     ALLZ(6,6,1) = 1.886067216D-02
 
     ALLC(1,6,1) = 4.554359511D-03
     ALLC(2,6,1) = 5.286443143D-02
     ALLC(3,6,1) =-7.561016358D-01
     ALLC(4,6,1) =-2.269803820D-01
     ALLC(5,6,1) = 1.332494651D+00
     ALLC(6,6,1) = 3.622518293D-01
!                                      6P
     ALLZ(1,6,2) = 6.696537714D-01
     ALLZ(2,6,2) = 1.395089793D-01
     ALLZ(3,6,2) = 8.163894960D-02
     ALLZ(4,6,2) = 4.586329272D-02
     ALLZ(5,6,2) = 2.961305556D-02
     ALLZ(6,6,2) = 1.882221321D-02
 
     ALLC(1,6,2) = 2.782723680D-03
     ALLC(2,6,2) =-1.282887780D-01
     ALLC(3,6,2) =-2.266255943D-01
     ALLC(4,6,2) = 4.682259383D-01
     ALLC(5,6,2) = 6.752048848D-01
     ALLC(6,6,2) = 1.091534212D-01

!Allocate the local arrays

     ng=MaxGaussianExpansion
     np=qmmm_struct%qm_ntypes
     
     allocate (atom_orb_cc_s_by_type(ng,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (atom_orb_cc_p_by_type(ng,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (atom_orb_cc_d_by_type(ng,np), stat=ier )
     REQUIRE(ier == 0)  
        
     allocate (atom_orb_zz_s_by_type(ng,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (atom_orb_zz_p_by_type(ng,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (atom_orb_zz_d_by_type(ng,np), stat=ier )
     REQUIRE(ier == 0)
     
     do i=1,np
       ni = qmmm_struct%qm_type_id(i)
       if(NI.LT.2) then
         NQN=1
      elseif(NI.LT.10)then
         NQN=2
      elseif(NI.LT.18)then
         NQN=3
      elseif(NI.LT.36)then
         NQN=4
      elseif(NI.LT.54)then
         NQN=5
       else
         NQN=6
       endif


       !  Fill out all orbitals even if atom type doesn't have p or d orbitals, they just won't be used.
       xi(1)= qm2_params%s_orb_exp_by_type(i)*qm2_params%s_orb_exp_by_type(i)
       xi(2)= qm2_params%p_orb_exp_by_type(i)*qm2_params%p_orb_exp_by_type(i)
       xi(3)= qm2_params%d_orb_exp_by_type(i)*qm2_params%d_orb_exp_by_type(i)              

       do j=1,ng
         atom_orb_cc_s_by_type(j,i)=ALLC(j,NQN,1)
         atom_orb_zz_s_by_type(j,i)=ALLZ(j,NQN,1)*xi(1)

         atom_orb_cc_p_by_type(j,i)=ALLC(j,NQN,2)
         atom_orb_zz_p_by_type(j,i)=ALLZ(j,NQN,2)*xi(2)

         atom_orb_cc_d_by_type(j,i)=ALLC(j,NQN,3)
         atom_orb_zz_d_by_type(j,i)=ALLZ(j,NQN,3)*xi(3)
       end do
       
     end do

!Deallocate the original 2 parameter arrays here as they are no longer needed from this point on.
!  commented by Taisung Lee: they will be needed for analytic calculations of slater overlap
!
!     deallocate (qm2_params%s_orb_exp_by_type, stat = ier)
!     REQUIRE(ier == 0)
!     deallocate (qm2_params%p_orb_exp_by_type, stat = ier)
!     REQUIRE(ier == 0)
!     deallocate (qm2_params%d_orb_exp_by_type, stat = ier)
!     REQUIRE(ier == 0)
     
!Now, lets pre-compute a lot of orbital interaction solutions to save time later.
!Allocate the memory required
     allocate (qm2_params%atom_orb_zz_sxs_over_sas(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_zz_sxp_over_sap(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_zz_sxd_over_sad(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)     
     allocate (qm2_params%atom_orb_zz_pxp_over_pap(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_zz_pxd_over_pad(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)  
     allocate (qm2_params%atom_orb_zz_dxd_over_dad(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)         
     allocate (qm2_params%atom_orb_ss_eqn(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_sp_ovlp(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_sd_ovlp(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0) 
     allocate (qm2_params%atom_orb_pd_ovlp(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)         
     allocate (qm2_params%atom_orb_pp_ovlp_inj(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_pp_ovlp_ieqj1(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_pp_ovlp_ieqj2(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_dd_ovlp_inj(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_dd_ovlp_ieqj1(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_dd_ovlp_ieqj2(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)     
     
     allocate (qm2_params%atom_orb_ss_eqn_adb(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_sp_eqn_xy(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_sp_eqn_xx1(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_sp_eqn_xx2(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_pp_eqn_xxy1(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
     allocate (qm2_params%atom_orb_pp_eqn_xxy2(ng,ng,np,np), stat=ier )
     REQUIRE(ier == 0)
  
     !Ross Walker - The following pre-computation saves a LOT of time in calculating the overlap energy
     !              and derivatives. This routine should only be called once so this loop here doesn't need to
     !              be amazingly efficient.

     !We do p orbital expansions even when the system doesn't have p orbitals. This is just to make things
     !simpler, the result of the p orbital calculations on systems without p orbitals will be garbage but we
     !won't reference the value so it shouldn't matter. Note we have to put in checks in this situation
     !to avoid divide by zeros.
     do k=1,np
       do l=1,np
         do i=1,6
           do j=1,6
             atom_orb_cc_s_x_s = atom_orb_cc_s_by_type(i,k)*atom_orb_cc_s_by_type(j,l)
             atom_orb_zz_s_x_s = atom_orb_zz_s_by_type(i,k)*atom_orb_zz_s_by_type(j,l)

             atom_orb_cc_s_x_p = atom_orb_cc_s_by_type(i,k)*atom_orb_cc_p_by_type(j,l)
             atom_orb_zz_s_x_p = atom_orb_zz_s_by_type(i,k)*atom_orb_zz_p_by_type(j,l)

             atom_orb_cc_s_x_d = atom_orb_cc_s_by_type(i,k)*atom_orb_cc_d_by_type(j,l)
             atom_orb_zz_s_x_d = atom_orb_zz_s_by_type(i,k)*atom_orb_zz_d_by_type(j,l)
                          
             atom_orb_cc_p_x_p = atom_orb_cc_p_by_type(i,k)*atom_orb_cc_p_by_type(j,l)
             atom_orb_zz_p_x_p = atom_orb_zz_p_by_type(i,k)*atom_orb_zz_p_by_type(j,l)
             
             atom_orb_cc_p_x_d = atom_orb_cc_p_by_type(i,k)*atom_orb_cc_d_by_type(j,l)
             atom_orb_zz_p_x_d = atom_orb_zz_p_by_type(i,k)*atom_orb_zz_d_by_type(j,l)             

             atom_orb_cc_d_x_d = atom_orb_cc_d_by_type(i,k)*atom_orb_cc_d_by_type(j,l)
             atom_orb_zz_d_x_d = atom_orb_zz_d_by_type(i,k)*atom_orb_zz_d_by_type(j,l)
             
             temp_real = atom_orb_zz_s_by_type(i,k)+atom_orb_zz_s_by_type(j,l)
             if (temp_real /= 0.0d0) atom_orb_zz_one_s_a_s = 1.0d0/temp_real
             
             temp_real = atom_orb_zz_s_by_type(i,k)+atom_orb_zz_p_by_type(j,l)
             if (temp_real /= 0.0d0) atom_orb_zz_one_s_a_p = 1.0d0/temp_real
             
             temp_real = atom_orb_zz_s_by_type(i,k)+atom_orb_zz_d_by_type(j,l)
             if (temp_real /= 0.0d0) atom_orb_zz_one_s_a_d = 1.0d0/temp_real
                        
             temp_real = atom_orb_zz_p_by_type(i,k)+atom_orb_zz_p_by_type(j,l)
             if (temp_real /= 0.0d0) atom_orb_zz_one_p_a_p = 1.0d0/temp_real

             temp_real = atom_orb_zz_p_by_type(i,k)+atom_orb_zz_d_by_type(j,l)
             if (temp_real /= 0.0d0) atom_orb_zz_one_p_a_d = 1.0d0/temp_real

             temp_real = atom_orb_zz_d_by_type(i,k)+atom_orb_zz_d_by_type(j,l)
             if (temp_real /= 0.0d0) atom_orb_zz_one_d_a_d = 1.0d0/temp_real

             qm2_params%atom_orb_zz_sxs_over_sas(i,j,k,l)=atom_orb_zz_s_x_s*atom_orb_zz_one_s_a_s
             qm2_params%atom_orb_zz_sxp_over_sap(i,j,k,l)=atom_orb_zz_s_x_p*atom_orb_zz_one_s_a_p
             qm2_params%atom_orb_zz_sxd_over_sad(i,j,k,l)=atom_orb_zz_s_x_d*atom_orb_zz_one_s_a_d
             
             qm2_params%atom_orb_zz_pxp_over_pap(i,j,k,l)=atom_orb_zz_p_x_p*atom_orb_zz_one_p_a_p
             qm2_params%atom_orb_zz_pxd_over_pad(i,j,k,l)=atom_orb_zz_p_x_d*atom_orb_zz_one_p_a_d
             qm2_params%atom_orb_zz_dxd_over_dad(i,j,k,l)=atom_orb_zz_d_x_d*atom_orb_zz_one_d_a_d
                          
             qm2_params%atom_orb_ss_eqn(i,j,k,l)=sqrt((2.0D0*sqrt(atom_orb_zz_s_x_s) &
                                                 *atom_orb_zz_one_s_a_s)**3) &
                                                 *atom_orb_cc_s_x_s
             atom_orb_sp_eqn=sqrt((2.0D0*sqrt(atom_orb_zz_s_x_p) &
                                                 *atom_orb_zz_one_s_a_p)**3) &
                                                 *atom_orb_cc_s_x_p
             atom_orb_sd_eqn=sqrt((2.0D0*sqrt(atom_orb_zz_s_x_d) &
                                                 *atom_orb_zz_one_s_a_d)**3) &
                                                 *atom_orb_cc_s_x_d


             !SQRT((two*SQRT(APB)*AMB)**3)*atom_orb_cc_p_x_p
             atom_orb_pp_eqn=sqrt((2.0D0*sqrt(atom_orb_zz_p_x_p) &
                                                 *atom_orb_zz_one_p_a_p)**3) &
                                                 *atom_orb_cc_p_x_p
             atom_orb_pd_eqn=sqrt((2.0D0*sqrt(atom_orb_zz_p_x_d) &
                                                 *atom_orb_zz_one_p_a_d)**3) &
                                                 *atom_orb_cc_p_x_d
             atom_orb_dd_eqn=sqrt((2.0D0*sqrt(atom_orb_zz_d_x_d) &
                                                 *atom_orb_zz_one_d_a_d)**3) &
                                                 *atom_orb_cc_d_x_d

             !2.0D0 * atom_orb_zz_s_by_type(K,qmitype)*SQRT(atom_orb_zz_p_by_type(L,qmjtype))
             !*atom_orb_zz_one_s_a_p*atom_orb_sp_eqn
             !Used in gover for Si-Pj overlap energy and -Sj-Pi overlap.
             qm2_params%atom_orb_sp_ovlp(i,j,k,l)=2.0D0*atom_orb_zz_s_by_type(i,k)* &
                                                  sqrt(atom_orb_zz_p_by_type(j,l))* &
                                                  atom_orb_zz_one_s_a_p*atom_orb_sp_eqn

             !-4.0D0*sqrt(atom_orb_zz_p_x_p)* 
             !atom_orb_zz_one_p_a_p*qm2_params%atom_orb_zz_pxp_over_pap(k,l,qmitype,qmjtype)* &
             !atom_orb_pp_eqn
             !Used in gover for Pi-Pj overlap energy when i!=j
             qm2_params%atom_orb_pp_ovlp_inj(i,j,k,l)=-4.0D0*sqrt(atom_orb_zz_p_x_p)* &
                                                      atom_orb_zz_one_p_a_p* &
                                                      qm2_params%atom_orb_zz_pxp_over_pap(i,j,k,l)* &
                                                      atom_orb_pp_eqn

             !-4.0D0*atom_orb_pp_eqn*
             !sqrt(atom_orb_zz_p_x_p)*
             !atom_orb_zz_one_p_a_p* 
             !qm2_params%atom_orb_zz_pxp_over_pap(k,l,qmitype,qmjtype)
             !Used in gover for Pi-Pj overlap energy when i==j
             qm2_params%atom_orb_pp_ovlp_ieqj1(i,j,k,l)=-4.0D0*atom_orb_pp_eqn* &
                                                        sqrt(atom_orb_zz_p_x_p)*atom_orb_zz_one_p_a_p* &
                                                        qm2_params%atom_orb_zz_pxp_over_pap(i,j,k,l)
             !2.0D0*atom_orb_pp_eqn* 
             !sqrt(atom_orb_zz_p_x_p)*atom_orb_zz_one_p_a_p
             !Used in gover for Pi-Pj overlap energy when i==j
             qm2_params%atom_orb_pp_ovlp_ieqj2(i,j,k,l)=2.0D0*atom_orb_pp_eqn* &
                                                        sqrt(atom_orb_zz_p_x_p)*atom_orb_zz_one_p_a_p

              !---specifics for QM-QM derivatives
              !-2.0D0*A2_TO_BOHRS2*qm2_params%atom_orb_ss_eqn_adb(i,j,qmitype,qmjtype)*qm2_params%atom_orb_zz_sxs_over_sas(i,j,qmitype,qmjtype)
              !Used for S-S overlap in QM-QM derivatives
              qm2_params%atom_orb_ss_eqn_adb(i,j,k,l)=-2.0D0*A2_TO_BOHRS2*qm2_params%atom_orb_ss_eqn(i,j,k,l)* &
                                                      qm2_params%atom_orb_zz_sxs_over_sas(i,j,k,l)

              !-four*A3_TO_BOHRS3*qm2_params%atom_orb_zz_sxp_over_sap(i,j,qmitype,qmjtype)**2* &
              !(one/(SQRT(atom_orb_zz_p_by_type(J,qmjtype))))*atom_orb_sp_eqn 
              !Used for S-P overlap in QM-QM derivatives where P... /= axis
              if (atom_orb_zz_p_by_type(j,l)/=0.0d0) then
                temp_real = 1.0d0/SQRT(atom_orb_zz_p_by_type(j,l))
              else
                temp_real = 0.0d0
              endif 
              qm2_params%atom_orb_sp_eqn_xy(i,j,k,l)=-4.0D0*A3_TO_BOHRS3* &
                                                     qm2_params%atom_orb_zz_sxp_over_sap(i,j,k,l)**2* &
                                                     (temp_real*atom_orb_sp_eqn) 
              
              !Used for S-P overlap in QM-QM derivatives where P... == axis
              if (atom_orb_zz_p_by_type(j,l)/=0.0d0) then
                temp_real =1.0D0/SQRT(atom_orb_zz_p_by_type(j,l))
              else 
                temp_real = 0.0d0
              end if
              qm2_params%atom_orb_sp_eqn_xx1(i,j,k,l)=2.0D0*A_TO_BOHRS*qm2_params%atom_orb_zz_sxp_over_sap(i,j,k,l)* &
                                                      temp_real * atom_orb_sp_eqn
              qm2_params%atom_orb_sp_eqn_xx2(i,j,k,l)=4.0D0*A3_TO_BOHRS3*qm2_params%atom_orb_zz_sxp_over_sap(i,j,k,l)**2* &
                                                      temp_real*atom_orb_sp_eqn

              !-four*A2_TO_BOHRS2*(ADB_array(inner_index)**2)*(one/(SQRT(atom_orb_zz_p_x_p)))*
              !atom_orb_pp_eqn
              !Used for P-P overlap in QM-QM derivatives where P... = P... /= axis as 3.0D*xxx when P... = P... == axis
              if (atom_orb_zz_p_x_p/=0.0d0) then
                 temp_real = 1.0d0/SQRT(atom_orb_zz_p_x_p)
              else
                 temp_real = 0.0d0
              end if
              qm2_params%atom_orb_pp_eqn_xxy1(i,j,k,l)=-4.0D0*A2_TO_BOHRS2* &
                                                       (qm2_params%atom_orb_zz_pxp_over_pap(i,j,k,l)**2)* &
                                                       temp_real*atom_orb_pp_eqn
              !eight*A2_TO_BOHRS2*A2_TO_BOHRS2*(ADB_array(inner_index)**3)* 
              !(one/(SQRT(atom_orb_zz_p_x_p)))*atom_orb_pp_eqn
              qm2_params%atom_orb_pp_eqn_xxy2(i,j,k,l)=8.0D0*A2_TO_BOHRS2*A2_TO_BOHRS2* &
                                                       (qm2_params%atom_orb_zz_pxp_over_pap(i,j,k,l)**3)* &
                                                       temp_real*atom_orb_pp_eqn
           end do
         end do
       end do
     end do

!Finally deallocate no longer needed arrays
     deallocate (atom_orb_zz_s_by_type, stat = ier)
     REQUIRE(ier == 0)
     deallocate (atom_orb_zz_p_by_type, stat = ier)
     REQUIRE(ier == 0)
     deallocate (atom_orb_zz_d_by_type, stat = ier)
     REQUIRE(ier == 0)
     deallocate (atom_orb_cc_s_by_type, stat = ier)
     REQUIRE(ier == 0)
     deallocate (atom_orb_cc_p_by_type, stat = ier)
     REQUIRE(ier == 0)
     deallocate (atom_orb_cc_d_by_type, stat = ier)
     REQUIRE(ier == 0)

     return
end subroutine qm2_setup_orb_exp

