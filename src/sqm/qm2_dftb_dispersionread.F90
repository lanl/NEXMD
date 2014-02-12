! <compile=optimized> 
!  -*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.fh"
subroutine dispersionread(nn,ntype,izp,disp_file) 
! This is in the dispfile structure: ,read_DISP_DOT_INP,h1,h2,Ni0,scale)

   use qmmm_module, only : qmmm_nml, qmmm_struct
   use ElementOrbitalIndex, only : elementSymbol
   use qm2_dftb_module, only: mol, dispertmp,dispfile
   implicit none
   integer:: i, j

!! Passed in:
   integer, intent(in ) :: nn
   integer, intent(in ) :: ntype
   integer, intent(in ) :: izp(*)  ! (NNDIM)
   character(len=*), intent(in ) :: disp_file

!! Locals
   logical :: found
   character(len=2) atom
   integer ier

!! Pointers
   character(len=2), dimension(:), pointer :: atyp

!! Temporary:
   logical :: ex

!! Set pointers
   atyp => mol%atyp

   ! Start by assuming the information is in DISPERSION.INP_ONCHSP
   dispfile%read_DISP_DOT_INP = .false.
   open(54,file=TRIM(disp_file),status="unknown")

   inquire(FILE="DISP.INP",EXIST=dispfile%read_DISP_DOT_INP) 

   ! parameters as in Elstner et al., JCP 2000
   dispfile%scale=1.0

   write(6,*)
   write(6,*) "******************************************************************************"
   if (dispfile%scale <= 0.0) then
      write(6,*) '  DFTB DISPERSION: London type dispersion, pol and dis in eV,scale=',dispfile%scale
   else
      write(6,*) '  DFTB DISPERSION: Slater-Kirkwood dispersion switched on.'
   endif
   write(6,*)
   write(6,'(A)') '|  Dispersion parameters read from file:'
   if (dispfile%read_DISP_DOT_INP) then
      write(6,'(A,3X,A)') '|','DISP.INP'
   else
      write(6,'(A,3X,A)') '|',TRIM(disp_file)
   end if
   write(6,*)

   if (dispfile%read_DISP_DOT_INP) then
      open(16,FILE="DISP.INP")
      write(6,*)  "   Parameters read from DISP.INP file:"
      write(6,*)  "       Atom        hh1       hh2        Ni"
      do i=1,nn
         read(16,*)  dispfile%hh1(i),dispfile%hh2(i),dispfile%Ni(i)
         write(6,'(3X,I6,": ",A2,3F10.6)') i, elementSymbol( qmmm_struct%iqm_atomic_numbers(i) ), &
               dispfile%hh1(i),  dispfile%hh2(i), dispfile%Ni(i) 
      end do
      close(16)

!!      return
!!
!!      write(6,*) "Did not return!!! why?"
!!   end if

   else
      do i = 1, nn
         atyp( izp(i) ) = elementSymbol( qmmm_struct%iqm_atomic_numbers(i) )
      end do
      
      do i=1,ntype
         found = .false.
         do while ( (.not.found) )
            if ( dispfile%scale >= 0.0 ) then
               read(54,*,iostat=ier) atom, (dispfile%h1(i,j),j=1,4),&
                     (dispfile%h2(i,j),j=1,4),&
                     dispfile%Ni0(i)
            else
               read(54,*,iostat=ier) atom, (dispfile%h1(i,j),j=1,4),&
                     (dispfile%h2(i,j),j=1,4)
            endif
            !If we hit the end of the file then we didn't find the parameters
            if ( ier<0 ) then
              write(6, '(" Dispersion parameters for atom ",A," not found.")') TRIM(atyp(i))
              call sander_bomb("<qm2_dftb_dispersionread> qm2_dftb_dispersionread.f",&
                    "Missing dispersion parameter.", "Exiting.")
            end if
            if (TRIM(atom) == TRIM(atyp(i))) then
               found = .true.
               rewind(54)
            end if
         end do
         
         ! If we got here, it is because we did find the atom we are looking for.
         write(6,'(1X,"Atom: ",A)')  atom
         write(6,'(2X,"h1  =",4f7.3)') (dispfile%h1(i,j),j=1,4)
         write(6,'(2X,"h2  =",4f7.3)') (dispfile%h2(i,j),j=1,4)
         if  ( dispfile%scale >= 0.0 ) write(6,'(2X,"Ni0 =",f7.3)')dispfile%Ni0(i)
      enddo
      
   end if !! if (dispfile%read_DISP_DOT_INP) then ... else ...

   !-----------------!
   ! Warning message !
   !-----------------!
   !<--
   write(6,*)
   write(6,*) "  WARNING: This is still being developed. Check your results carefully."
   write(6,*) "  Careful: The specific parameter to be used will be determined solely from"
   write(6,*) "           the number of atoms within 1.2 Angstroms from the central"
   write(6,*)
   write(6,*) "           ALSO, note that these paramters were developed for DNA bases only."
   write(6,*)
   write(6,*) "  For details about the theory, see: "
   write(6,*) '      "Hydrogen bonding and stacking interactions of nucleic acid base pairs:'
   write(6,*) '       a density-functional-theory based treatment",'
   write(6,*) "       M. Elstner, P. Hobza, T. Frauenheim,  S. Suhai, and E. Kaxiras, "
   write(6,*) "       J. Chem.  Phys. , 114, 5149 (2001)."
   write(6,*) "******************************************************************************"
   write(6,*)
   !-->

   return
end subroutine dispersionread





