#include "dprec.fh"
!------------------------------------------------------------------------------
module lscivr_vars
!
! internal variables used in LSC-IVR
!------------------------------------------------------------------------------
!
  implicit none
  save
!    
  integer :: file_pos_lsc = 3001
  integer :: file_vel_lsc = 3002
  integer :: file_rhoa_lsc = 3003
  integer :: file_cor_lsc = 3004  
  integer :: file_Tmat_lsc = 3005
! the sign to start the LSC-IVR calculation
! ilscivr = 0,     no LSC-IVR 
!          = 1,     LSC-IVR
  integer :: ilscivr  
! the number of LSC-IVR trajectories 
  integer:: ntraj_lsc
! the sign to determine which type of the correlation function
! icorf_lsc = 0,       scalar function
!             1,       position <x(0)x(t)>
!             2,       momentum <v(0)v(t)>
!             3,       nonlinear function <f(0)f(t)>         
!             4,       Kubo-transformed <v(0)v(t)>
! Also see lsc_init.f & lsc_xp.f
  integer :: icorf_lsc
! an auxiliary variable
  integer :: ilsc
! # of atoms and degrees of freedom
  integer :: ndof_lsc, natom_lsc
! # of time steps
  integer :: nstep_lsc
! time along the trajectory
  _REAL_ :: t_lsc
! define the inverse temperature
  _REAL_ :: beta_lsc
! the small displacement used to calculate the derivative  
  _REAL_ :: dx_lsc
! a variable for the Gaussian distribution of each mode
  _REAL_ :: sigma_lsc
! frequency of the mode
  _REAL_ :: omega_lsc
! mass array
  _REAL_, allocatable :: mass_lsc(:) 
! a 1-dim auxiliary array to help compute the 2nd derivative 
!   of the mass-weighted potential
  _REAL_, allocatable ::  x_lsc(:)
! a 1-dim auxiliary array for the force
  _REAL_, allocatable ::  f_lsc(:)
! a 1-dim auxiliary array for the frequency 
  _REAL_, allocatable :: v2_freq(:)
! a 2-dim array for the 2nd derivative of the mass-weighted potential
  _REAL_, allocatable :: v2_lsc(:,:)
! a 1-dim array for the mode
  _REAL_, allocatable :: alpha_lsc(:)
! a 1-dim array for the momentum generated from the Wigner distribution
  _REAL_, allocatable :: p_lsc(:)
! a 1-dim auxiliary array for diagonalization
  _REAL_, allocatable :: work_lsc(:) 
!
! ---------------------------------------------------------
! start to define the (exp(-beta*H)*A) part of 
!  the correlation function 
! this variable will need to be written to an output file 
! ---------------------------------------------------------
  _REAL_, allocatable :: rho_a(:) 

  _REAL_:: Ek_lsc
! for gaussian variable

end module lscivr_vars
