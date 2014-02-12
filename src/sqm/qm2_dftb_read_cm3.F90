! <compile=optimized>
!  -*- mode: f90; coding: iso-8859-15; -*-
#include "dprec.fh"

subroutine read_cm3(natoms,ntype,izp,cm3_file)
   !
   ! Subroutine to read the CM3 parameters from the file
   ! $AMBERHOME/dat/slko/CM3_PARAMETERS.DAT
   !
   use qmmm_module, only: qmmm_nml, qmmm_struct
   use qm2_dftb_module, only: mol, cm3

!! Passed in:
   integer, intent(in) :: natoms
   integer, intent(in) :: ntype
   integer, intent(in) :: izp(*)  ! (NNDIM)
   character(len=*), intent(in ) :: cm3_file

!! Locals
   logical :: found
   character(len=2) :: atom1, atom2,atypi,atypj
   integer line, i, j
   _REAL_ :: temp_cm3_c, temp_cm3_d

   !--

   found = .false.
   inquire(FILE=TRIM(cm3_file), EXIST=found)
   open(55,file=TRIM(cm3_file),status="unknown")
   if (.not.found) then
      write(6,*)
      write(6,*)"****************************************************"
      write(6,*)"*     !! A FILE NEEDED BY DFTB WAS NOT FOUND !!    *"
      write(6,*)"****************************************************"
      write(6,*)
      write(6,*)" Sander could not find the file containing the "
      write(6,*)" CM3 parameter table needed for this calculation."
      write(6,*)
      write(6,*)" This file was supposed to be located in the "
      write(6,*)" $(AMBERHOME)/dat/slko/"
      write(6,*)" directory."
      write(6,*)
      call sander_bomb("qm2_dftb_read_cm3 <qm2_dftb_read_cm3.f>", &
            "File not found.", &
            "Exiting.")
   else
      write(6,'("|")')
      write(6,'("|******************************************************************************")')
      write(6,'("|")')
      write(6,'("|  The DFTB-CM3 parametrization is described at:")')
      write(6,'("|  Kalinowski et al., J. Phys. Chem. A 2004, 108, 2545 ")')
      write(6,'("|")')
      write(6,'("|  The parameters used here were kindly provided by Dr. Jaroslaw Kalinowski by")')
      write(6,'("|  private communication. Note that the C parameter for the S-P interaction")')
      write(6,'("|  used here differs from the published value.")')
      write(6,'("|")')
      write(6,'("|  Parameters read from file: ")')
      write(6,'("|  ",A)')cm3_file
      write(6,'("|")')
      write(6,'("|   1  (TYP)   2  (TYP)         D            C   ")')
      !         "|   H  (  1)   H  (  1)      0.00000      0.00000"
   end if

   do i=1,ntype
      write(6,'("|")')
      atypi = mol%atyp(i)
      atom1 = atypi
      do j = 1, ntype
         atypj = mol%atyp(j)
         atom2 = atypj

         cm3%c(i,j) = 0.0d0
         cm3%d(i,j) = 0.0d0

         if (i /= j) then
            found = .false.
            line = 0
            temp_cm3_d = 0.0d0
            temp_cm3_c = 0.0d0
            rewind(55)
            do while ( (.not.found) .and. (line <= cm3%num_params) )
               line = line + 1
               read(UNIT=55,FMT='(2(A2,2X),2F11.6)',end=10) atom1,atom2,temp_cm3_d, temp_cm3_c
               if ( TRIM(atom1) == TRIM(atypi) .and. TRIM(atom2) == TRIM(atypj) ) then
                  found = .true.
                  cm3%d(i,j) = temp_cm3_d
                  cm3%c(i,j) = temp_cm3_c
               else if ( TRIM(atom1) == TRIM(atypj) .and. TRIM(atom2) == TRIM(atypi) ) then
                  !c|d(j,i) = -c|d(i,j)
                  found = .true.
                  cm3%d(i,j) = -temp_cm3_d
                  cm3%c(i,j) = -temp_cm3_c
               end if
            end do
            ! If line > cm3%num_params it means we couldn't find this atom in the parameter file. Bomb.
            if ( line > cm3%num_params ) then
               write(6, '(" CM3 parameters for bond between atoms ",A," and ",A," not found.")') TRIM(atypi), TRIM(atypj)
               call sander_bomb("<qm2_dftb_cm3_read> qm2_dftb_cm3_read.f",&
                     "Missing CM3 parameter.", "Exiting.")
            end if
         end if

         ! If we got here, it is because we did find the parameter we are looking for.
         write(6,'("|",3X,A2," (",I3,")",3X,A2," (",I3,")",3X,f10.6,3X,f10.6)') &
               atypi,i,atypj,j,cm3%d(i,j),cm3%c(i,j)
      end do
   end do

   !<--
   write(6,'("|")')
   write(6,'("|******************************************************************************")')
   write(6,'("|")')
   !-->

   goto 20
10 call sander_bomb("qm2_dftb_read_cm3 <qm2_dftb_read_cm3.f>", &
            "End of file!!!", &
            "Exiting.")
   continue

20 return

end subroutine read_cm3
