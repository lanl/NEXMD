! <compile=optimized> 
!  -*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.fh"


subroutine check_dftb_3rd_order(dftb_3rd_order)

   ! Called during the reading of the qmmm namelist, this
   ! routine checks the validity of the 3rd order keyword.
   !
   ! Allowed options are:
   !  'PA'     -->   Use the Proton Affinities parameterization
   !  'PR'     -->   Use the Phosphate Reactions parameterization
   !  'READ'   -->   Reads the parameters from namelist 'dftb_3rd_order'
   !                 in the input file
   !  '<file>' -->   Reads the namelist from file <file>. (NOT IMPLEMENTED YET)

   use qm2_dftb_module, only: DFTB_3rd_order_str
   implicit none
   character(len=4), intent(in) :: dftb_3rd_order

   dftb_3rd_order_str%do_3rd_order = .true.
   return
   
end subroutine check_dftb_3rd_order

subroutine qm2_dftb_read_3rd_order(qmmm_struct, qm_natoms,qm_ntypes,izp,skroot) 

   use qmmm_module, only : qmmm_nml
   use ElementOrbitalIndex, only : elementSymbol
   use qm2_dftb_module, only: mol, DFTB_3rd_order_str
   use qmmm_struct_module, only : qmmm_struct_type
   implicit none
#  include "files.h"


!! Passed in:
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   integer, intent(in ) :: qm_natoms ! qmmm_struct%nquant_nlink
   integer, intent(in ) :: qm_ntypes ! qmmm_struct%qm_ntypes
   integer, intent(in ) :: izp(*)    ! izp_str%izp
   character(len=*), intent(in ) :: skroot ! The path to the sk-files directory
   
!! Locals
   character(len=1024) :: dftb_3rd_order_file
   integer this_at_number, j, dftb_3rd_order_unit, ier

!! Namelist input
   _REAL_ :: Gaussian_d0, Gaussian_g0, Gaussian_q0
   _REAL_, dimension(100) :: Hubbard_deriv
   logical :: debug_print
   namelist /dftb_3rd_order/ Gaussian_d0, Gaussian_g0, Gaussian_q0, &
                             Hubbard_deriv, debug_print
  
   if ( (qmmm_nml%dftb_3rd_order == 'PA').or. &
        (qmmm_nml%dftb_3rd_order == 'PR')) then
      dftb_3rd_order_file = &
         TRIM(skroot)//'DFTB_3RD_ORDER_'//TRIM(qmmm_nml%dftb_3rd_order)//'.DAT'
      dftb_3rd_order_unit = 54
   else if (qmmm_nml%dftb_3rd_order == 'READ') then
      ! Read namelist from mdin file
      dftb_3rd_order_file = mdin
      dftb_3rd_order_unit = 5 
   else
      dftb_3rd_order_unit = 54
      dftb_3rd_order_file = TRIM(qmmm_nml%dftb_3rd_order)
      ! Read namelist from external file
      ! (Not implemented yet)
      ! dftb_3rd_order_file = qmmm_nml%dftb_3rd_order_file
   endif
   write(6,*)
   write(6,'("QMMM: 3rd Order SCC-DFTB")')
   write(6,'("QMMM: ------------------")')
   write(6,'("QMMM: Reading 3rd Order parameters from file:")')
   write(6,'("| ",5X,A)'), TRIM(dftb_3rd_order_file)
   write(6,*)

   Hubbard_deriv(1:100) = 0.0d0
   debug_print = .false.
  
   if (dftb_3rd_order_unit == 54) then
      open(54,file=TRIM(dftb_3rd_order_file), status="OLD", IOSTAT=ier)
      if (ier /=0) call sander_bomb('qm2_dftb_read_3rd_order',&
                        'Could not find the file:', TRIM(dftb_3rd_order_file))
   end if

   ! If we got here, the file exists. Now we check for 
   ! the existence of the 'dftb_3rd_order' namelist

!   WRITE(6,*)"qm2_dftb_3rd_order:",dftb_3rd_order_unit,dftb_3rd_order_str%file_unit 

   ! gfortran bombs on /Run.ala8_PA without this with
   ! At line 202 of file _qm2_dftb_3rd_order.f
   ! Fortran runtime error: Illegal seek
   dftb_3rd_order_str%file_unit = dftb_3rd_order_unit

   rewind dftb_3rd_order_str%file_unit 
   call nmlsrc('dftb_3rd_order',dftb_3rd_order_unit,ier)

   if (ier /= 0) then
      !Read qmmm namelist
      rewind dftb_3rd_order_unit
      read(dftb_3rd_order_unit, nml=dftb_3rd_order)
   else
      call sander_bomb('qm2_dftb_read_3rd_order',&
                       'Could not find the "dftb_3rd_order" namelist in the file:',&
                       TRIM(dftb_3rd_order_file))
   endif

   ! Puts the parameters in the proper places
   DFTB_3rd_order_str%Gaussian_d0 = Gaussian_d0
   DFTB_3rd_order_str%Gaussian_g0 = Gaussian_g0
   DFTB_3rd_order_str%Gaussian_q0 = Gaussian_q0
   DFTB_3rd_order_str%Hubbard_deriv = Hubbard_deriv
   DFTB_3rd_order_str%debug_print = debug_print
  

   !Prints information 
   write(6,*)
   write(6,'(A)')"QMMM: Gaussian Parameters:"
   write(6,'(A,F8.3)')"QMMM:          D0 =    ", DFTB_3rd_order_str%Gaussian_d0
   write(6,'(A,F8.3)')"QMMM:          g0 =    ", DFTB_3rd_order_str%Gaussian_g0
   write(6,'(A,F8.3)')"QMMM:          q0 =    ", DFTB_3rd_order_str%Gaussian_q0
   write(6,*)
   write(6,'(A)')"QMMM: Hubbard Derivatives:"
   
   do j=1,qm_natoms
      this_at_number = qmmm_struct%iqm_atomic_numbers(j)
      write(6,'(A,I5,2X,"(",I2,")",A3,2X,F8.3)')"QMMM:  ",&
            j, &
            this_at_number, &
            elementSymbol( this_at_number ),&
            DFTB_3rd_order_str%Hubbard_deriv( this_at_number )
   end do

   
   return
end subroutine qm2_dftb_read_3rd_order




