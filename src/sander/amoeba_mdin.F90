#include "dprec.fh"
#include "assert.fh"

module amoeba_mdin
  implicit none
  private
  integer,save :: iamoeba=0,do_amoeba_valence=1,do_amoeba_nonbond=1, &
                  do_vdw_taper=1,do_vdw_longrange=1,am_nbead=1
  _REAL_,save :: sor_coefficient = 0.75d0
  _REAL_,save :: dipole_scf_tol = 0.01d0
  integer,save :: dipole_scf_iter_max = 50
#ifndef DISABLE_AMOEBA_CG
  integer,save :: dipole_scf_use_cg = 0
  integer,save :: dipole_scf_cg_niter = 4
#endif /* DISABLE_AMOEBA_CG */
  _REAL_,save :: ee_dsum_cut=7.d0
  _REAL_,save :: ee_damped_cut=4.5d0
  _REAL_,save :: vdw_taper = 0.9d0
  _REAL_,save :: thole_expon_coeff=0.39d0
  _REAL_,save :: compress = 0.000046d0
  _REAL_,save :: soft_lambda = 1.0d0
  _REAL_,save :: soft_alpha = 0.5d0
  _REAL_,save :: soft_expo = 4
  _REAL_,save :: vdw_longrange_lambda = 1.0d0
  integer, save :: amoeba_verbose = 0
  integer,save :: beeman_integrator = 0

  !variables and arrays for softcore file
  integer,save :: soft_atom_range1(20),soft_atom_range2(20),soft_line

  public AMOEBA_read_mdin,iamoeba,do_amoeba_valence, &
         do_amoeba_nonbond,amoeba_verbose,beeman_integrator, &
         sor_coefficient,dipole_scf_tol,dipole_scf_iter_max, &
#ifndef DISABLE_AMOEBA_CG
         dipole_scf_use_cg, dipole_scf_cg_niter, &
#endif /* DISABLE_AMOEBA_CG */
         ee_dsum_cut,ee_damped_cut,thole_expon_coeff,vdw_taper, &
         compress,do_vdw_taper,do_vdw_longrange,am_nbead, &
         soft_lambda,soft_alpha,soft_expo, &
         AMOEBA_read_soft, soft_atom_range1,soft_atom_range2, soft_line, &
         vdw_longrange_lambda
  contains
!-------------------------------------------------------------------------------
subroutine AMOEBA_read_mdin(nf)
  use file_io_dat

  implicit none
  integer,intent(in) :: nf

  integer            :: do_bond=1,do_ureyb=1,do_reg_angle=1,  &
                        do_trig_angle=1,do_opbend=1,do_torsion=1, &
                        do_pi_torsion=1,do_strbend=1,do_torsion_torsion=1, &
                        do_str_torsion=1, &
                        do_recip=1,do_adjust=1,do_direct=1,do_self=1, &
                        do_vdw=1,do_induced=1
  namelist/amoeba/do_amoeba_valence,do_amoeba_nonbond, &
                  do_bond,do_ureyb,do_reg_angle,  &
                  do_trig_angle,do_opbend,do_torsion,do_str_torsion, &
                  do_pi_torsion,do_strbend,do_torsion_torsion, &
                  do_recip,do_adjust,do_direct,do_self, &
                  do_vdw,do_induced,amoeba_verbose,beeman_integrator, & 
                  sor_coefficient,dipole_scf_tol, &
                  dipole_scf_iter_max,&
#ifndef DISABLE_AMOEBA_CG
                  dipole_scf_use_cg, dipole_scf_cg_niter, &
#endif /* DISABLE_AMOEBA_CG */
                  ee_dsum_cut,ee_damped_cut, &
                  thole_expon_coeff,vdw_taper,do_vdw_taper,do_vdw_longrange, &
                  compress,am_nbead,&
                  soft_lambda,soft_alpha,soft_expo,vdw_longrange_lambda

  read(nf,nml=amoeba)
  call AM_VAL_set_user_bit(do_bond,do_ureyb,do_reg_angle,do_trig_angle, &
                          do_opbend,do_torsion,do_str_torsion, &
                          do_pi_torsion,do_strbend, &
                          do_torsion_torsion)
  call AM_NONBOND_set_user_bit(do_recip,do_adjust,do_direct,do_self, &
                               do_vdw,do_induced)
  
end subroutine AMOEBA_read_mdin
!----------------------------------------------------------
subroutine AMOEBA_read_soft()
  integer pos_dash,length,i_range,i
  character(len=20) atm_range
  character(len=6) temp
  
  !flag whether soft_atm.txt exists
  logical alive
  inquire(file='soft_atm.txt',exist=alive)
  if(alive) then
     do i=1,20
        soft_atom_range1(i)=0;
        soft_atom_range2(i)=0;
     end do
   
     soft_line=0
     open (11,FILE='soft_atm.txt')
     do while (.true.)
     
        read(11,*,end=99) atm_range
        !write(*,*) atm_range
        pos_dash=scan(atm_range,'-')
        !write(*,*) "position",pos_dash
        length=len_trim(atm_range)
        if (length > 0) then
           soft_line=soft_line+1
        endif
        !write(*,*) "length",length
        if(pos_dash.gt.0) then
           temp=atm_range(1:pos_dash-1)
           read(temp,*) soft_atom_range1(soft_line)
           temp=atm_range(pos_dash+1:length)
           read(temp,*) soft_atom_range2(soft_line)
  
        else
           read(atm_range,*) soft_atom_range1(soft_line)
           read(atm_range,*) soft_atom_range2(soft_line)
        endif
      end do
   else
     write(6,*) 'soft_atm.txt not found'
     call mexit(6,1)
   endif
   99 continue
   close(11)
 
end subroutine AMOEBA_read_soft
!---------------------------------------------------------
end module amoeba_mdin
!-------------------------------------------------------------------------------
