! <compile=optimized> !  -*- mode: f90; coding: iso-8859-15; -*-
#include "dprec.fh"

subroutine qm2_dftb_cm3(qmmm_nml,qmmm_mpi,qmmm_struct, scf_mchg,scf_cm3)

!! Calculates the CM3 charges, according to 
!! Kalinowsli et al, JPC-A, 108, 2545 (2004).
!!
!!
!! After returning, the cm3%qcm3(j) array will contain
!! CM3 charges.
!!
!! -------------------------------------------------
!! Implementation in Amber/SCC-DFTB by
!!
!! GUSTAVO SEABRA
!! Quantum Theory Project, University of Florida
!! e-mail: seabra@qtp.ufl.edu
!! url:    http://www.qtp.ufl.edu/~seabra
!! -------------------------------------------------

   use qmmm_module, only: qmmm_mpi_structure
   use qm2_dftb_module, only: cm3, mol, ks_struct, izp_str, mcharge, NDIM, dacc
   use qmmm_struct_module, only : qmmm_struct_type
   use qmmm_nml_module   , only : qmmm_nml_type

   implicit none

!! Passed in
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   _REAL_, intent(in) :: scf_mchg(qmmm_struct%nquant_nlink) !Mulliken charge per atom.
   _REAL_, intent(out) :: scf_cm3(qmmm_struct%nquant_nlink)

!! Pointers
   integer, pointer :: ind(:)
   integer, pointer :: izp(:)
   _REAL_ , pointer :: qzero(:)     ! (NNDIM)     Atomic electron population of the neutral atom
   _REAL_ , pointer :: qmat(:)      ! (NNDIM)     Atomic electron population
   _REAL_ , pointer :: overl(:,:)   ! (MDIM,MDIM) Overlap Matrix
   _REAL_ , pointer :: density(:,:) ! (NDIM,NDIM) Density Matrix

!! Locals
   integer :: i,j, izpi, izpj, lambda, omega, adim
   integer :: lambda_begin, lambda_end, omega_begin, omega_end, l
   _REAL_  :: PS1, PS2, temp_b, temp_c, temp_d,temp_b2

!! For Matrix Printing   
   integer :: tot_n
   integer :: n_blocks_to_print, n_extra_columns,begin
   _REAL_  :: tr

!--

!! Set pointers
   ind     => ks_struct%ind
   overl   => ks_struct%overl
   density => ks_struct%density
   izp     => izp_str%izp
   qmat    => mol%qmat
   qzero   => mcharge%qzero
!!

   cm3%b = 0.0d0
   cm3%t = 0.0d0

!! Find the number of OCCUPIED MO's (adim)
!! MOVED HERE BY ROSS WALKER. It used to be stored in the module but this was
!! causing strange performance behaviour the mulliken calculation. So for the
!! time being we just recalculate it here and pass it to qm2_dftb_dens_matrix.
   do  i = 1,NDIM
      if(ks_struct%occ(i) < dacc) exit
   end do

   ! The number of occupied orbitals will get stored into adim
   adim = i-1

   ! Calculates the density matrix (DFTB doesn't by default)
   call qm2_dftb_dens_matrix(adim, ndim)

!! -------------------------------------------------------------
!! ---  Build Mayer bond order matrix, B: (the CM3 B-Matrix) ---
!! -------------------------------------------------------------
   do i = 1, qmmm_struct%nquant_nlink      ! (k)
      lambda_begin = ind(i)+1
      lambda_end   = ind(i+1)
      do j = 1, qmmm_struct%nquant_nlink
         omega_begin = ind(j)+1
         omega_end   = ind(j+1)
         temp_b = 0.0d0
         do lambda = lambda_begin, lambda_end
            do omega = omega_begin, omega_end
               PS1 = 0.0d0
               PS2 = 0.0d0
               do l = 1, NDIM
                  PS1 = PS1 + Density(omega ,l) * Overl(lambda,l)
                  PS2 = PS2 + Density(lambda,l) * Overl(omega ,l)
               end do
               temp_b = temp_b +  PS1 * PS2
            end do ! omega
         end do ! lambda
         cm3%b(i,j) = temp_b
      end do ! j
   end do ! i




!! -------------------------------------
!! ----  Calculates the CM3 T-Matrix ---
!! ------------------------------------- 
   do i = 1, qmmm_struct%nquant_nlink
      izpi = izp(i)
      do j = 1, qmmm_struct%nquant_nlink
         izpj = izp(j)
         temp_b  = cm3%b(i,j)
         temp_c  = cm3%c(izpi,izpj)
         temp_d  = cm3%d(izpi,izpj)
         temp_b2 = temp_b * temp_b
         cm3%t(i,j) =  temp_d * temp_b + temp_c * temp_b2 
      end do
   end do



!! ----------------------------------
!! ---  Calculate the CM3 charges ---
!! ----------------------------------
   do i = 1, qmmm_struct%nquant_nlink
      cm3%qcm3(i) = scf_mchg(i)
      do j = 1, qmmm_struct%nquant_nlink
         if (i /= j) then
            cm3%qcm3(i) = cm3%qcm3(i) + cm3%t(i,j)
         end if
      end do
   end do


   ! DEBUG: WRITE OUTPUT
   if (qmmm_nml%verbosity > 3 .and. qmmm_mpi%commqmmm_master) then
      write(6,*)
      write(6,'(" AT  (TYP)  Symb      Qzero        Qmat     QMullik        QCM3    new Qmat")')
      !         " 12  (  1)   H      1.00000     0.91905     0.08095     0.11405     0.88595
      do j = 1, qmmm_struct%nquant_nlink
         write(6,'(I3,2X,"(",I3,")",3X,A2,5(2X,F10.5))') &
               j, izp(j),mol%atyp(izp(j)),qzero( izp(j) ), qmat(j), scf_mchg(j),cm3%qcm3(j),(qzero( izp(j) ) - cm3%qcm3(j))
      end do
   end if


   do j = 1, qmmm_struct%nquant_nlink
      scf_cm3(j) = cm3%qcm3(j)
   end do

   return

!end subroutine qm2_dftb_cm3



end subroutine qm2_dftb_cm3

subroutine print_orbitals(nstep,eigenvalues,eigenvectors, m, n, ldm, first, last)
   implicit none
   
   !! Passed in
   integer, intent(in) :: nstep ! MD step number  
   integer, intent(in) :: m     ! Number of colums
   integer, intent(in) :: n     ! Number of lines
   integer, intent(in) :: ldm   ! Linear dimension of the matrix
                                ! (number of atomic orbitals)
   integer, intent(in) :: first ! First column to print
   integer, intent(in) :: last  ! Last column to print
   _REAL_ , intent(in) :: eigenvalues(ldm)
   _REAL_ , intent(in) :: eigenvectors(ldm,n)


   !! Local
   integer :: n_blocks_to_print
   integer :: n_extra_columns
   integer :: n_columns_to_print
   integer :: n_columns_by_block
   integer :: begin, i, j, l, iout


   iout = 35 

   n_columns_to_print = last - first + 1
   n_blocks_to_print  = int(n_columns_to_print / 5)
   n_extra_columns = n_columns_to_print - (5 * n_blocks_to_print)
   if (n_extra_columns < 0) n_extra_columns = 0

   begin = first
   n_columns_by_block = 5
   if (n_columns_to_print < 5) n_columns_by_block = n_columns_to_print

   write(iout,*) " "
   write(iout,*) " "
   write(iout,*) " ===================================="
   write(iout,*) " ORBITALS FROM MD STEP NUMBER", nstep
   write(iout,*) " ===================================="

   do i = 1, n_blocks_to_print
      write(iout,'("")')
      !Matrix header
      write(iout,'(14X,5(3X,I7))') &
      (l, l = begin, begin + n_columns_by_block -1)
      write(iout,'("")')
      !Eigenvalues
      write(iout,'(" Eigenvalues: ",5(F10.5))') &
      (eigenvalues(l), l = begin, begin + n_columns_by_block -1)
      write(iout,'("")')
      !Eigenvectors
      do j = 1, n
         write(iout,'(I14,5(F10.5))')j,&
         (eigenvectors(j,l),l = begin, begin + n_columns_by_block - 1)
      end do
      begin = begin + 5 
      write(iout,'("")')
   end do

   if (n_extra_columns.gt.0) then
      write(iout,'("")')
      !Matrix header 
      write(iout,'(14X,5(3X,I7))') &
      (l, l = begin, begin + n_extra_columns - 1)
      write(iout,'("")')
      !Eigenvalues
      write(iout,'(" Eigenvalues: ",5(F10.5))') &
      (eigenvalues(l), l = begin, begin + n_extra_columns -1)
      write(iout,'("")')
      !Eigenvectors
      do j = 1, n
         write(iout,'(I14,5(F10.5))')j,&
         (eigenvectors(j,l),l=begin,begin+n_extra_columns-1)
      end do
   end if

end subroutine print_orbitals
