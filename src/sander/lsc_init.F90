!#include "../include/assert.fh"
#include "dprec.fh"

subroutine lsc_init
!
! --------------------------------------------------------------
! In this subroutine, one has to call PIMD subroutine to 
!  prepare the configurations of beads
!
! beta, mass, and # of atoms
! set up the value for the displacement dx
! 
   use lscivr_vars
   use constants, only: kB
   implicit none
#include "md.h"
#include "extra.h"
!   
   integer:: ierr
!
   open( file_pos_lsc, file = 'LSC_position.dat')
   open( file_vel_lsc, file = 'LSC_velocity.dat')
   open( file_rhoa_lsc, file = 'LSC_rhoa.dat')
   open( file_cor_lsc, file = 'LSC_beta_A.dat')
   open( file_Tmat_lsc, file = 'LSC_Tmat.dat')

   beta_lsc = 1.0d0/(kB * temp0)
   
!
   dx_lsc = 1.0d-5
!
!   
   allocate(mass_lsc(ndof_lsc),stat=ierr)
!
   allocate(v2_freq(ndof_lsc),stat=ierr)
!
   allocate(x_lsc(ndof_lsc),stat=ierr)
   allocate(f_lsc(ndof_lsc),stat=ierr)
!
   allocate(v2_lsc(ndof_lsc, ndof_lsc), stat=ierr)
!
   allocate(alpha_lsc(ndof_lsc), stat=ierr) 
!
   allocate(p_lsc(ndof_lsc), stat=ierr)
!
   allocate(work_lsc(10*ndof_lsc), stat=ierr)
! --------------------------------------------------------------
! the sign to determine which type of the correlation function
! icorf_lsc = 0,       scalar function
!             1,       position <x(0)x(t)>
!             2,       momentum <v(0)v(t)>
!             3,       nonlinear vector function <f(0)f(t)>
!             4,       Kubo-transformed <v(0)v(t)>
! Also see lscivr_vars.f & lsc_xp.f
! --------------------------------------------------------------
end subroutine lsc_init
!
!-------------------------------------------------------------------
! lsc_final is the subroutine to deallocate arrays used by the LSC
!-------------------------------------------------------------------
subroutine lsc_final
!
   use lscivr_vars
   implicit none
   integer:: alloc_error
!
   deallocate(mass_lsc, stat = alloc_error )
!
   deallocate(x_lsc, stat = alloc_error )
!
   deallocate(f_lsc, stat = alloc_error )
!
   deallocate(v2_lsc, stat = alloc_error )
!
   deallocate(alpha_lsc, stat = alloc_error )
!
   deallocate(p_lsc, stat = alloc_error )
!
   deallocate(work_lsc, stat = alloc_error )
!
   deallocate(rho_a, stat = alloc_error)
!
   close( file_pos_lsc)
   close( file_vel_lsc)
   close( file_rhoa_lsc)
   close( file_cor_lsc )
   close( file_Tmat_lsc )

!
end subroutine lsc_final
