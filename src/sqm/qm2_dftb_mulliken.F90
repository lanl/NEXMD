! <compile=optimized>
!-*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.fh"

!! =============================================================================
!!                         MULLIKEN POPULATION ANALYSIS
!! =============================================================================
subroutine mulliken(qm2_struct,qmmm_struct, nquant_nlink,NDIM,izp,lmax,dacc,qmat,qzero,scf_mchg)

   use qmmm_module, only : qm2_structure, qmmm_nml, qmmm_mpi
   use qm2_dftb_module, only: MDIM, ks_struct, izp_str, cm3, mol
   use qmmm_struct_module, only : qmmm_struct_type

   implicit none

!!Passed in:
   type(qm2_structure),intent(inout) :: qm2_struct
   type(qmmm_struct_type), intent(in) :: qmmm_struct
   integer, intent(in)  :: nquant_nlink
   integer, intent(in)  :: NDIM
   integer, intent(in)  :: izp(*)
   integer, intent(in)  :: lmax(*)
   _REAL_ , intent(in)  :: dacc
   _REAL_ , intent(in)  :: qzero(*)   ! Electron numbers on neutral atom
   _REAL_ , intent(out) :: qmat(*)    ! Electron numbers
   _REAL_ , intent(out) :: scf_mchg(qmmm_struct%nquant_nlink) !Mulliken charge per atom.

!!Locals
   integer :: i, adim
   _REAL_ :: sum,qhelp,conv, temp, tr
   integer :: k
   logical :: calc_dens

!! Testing
   integer :: n, j, ljmax, indj, lj, jofn, mj, mjmax
   _REAL_  :: qtot

!! Find the number of OCCUPIED MO's (adim)
   do  i = 1,NDIM
      if(ks_struct%occ(i) < dacc) exit
   end do

   ! The number of occupied orbitals will get stored into adim
   adim = i-1 

!! =============================================================================
!!                  CALCULATION OF the MULLIKEN POPULATIONS
!! =============================================================================

   call qm2_dftb_qmulli(NDIM, MDIM, adim, calc_dens) ! Orbital-based Mulliken pop. (qmulli)
   call qm2_dftb_mull_pop(nquant_nlink,NDIM,izp,lmax,qmat)     ! Atom-based Mulliken pop.    (qmat)

   ! MULLIKEN CHARGES
   ! ----------------
   ! sander needs the Mulliken charges stored in scf_mchg(natoms)
   !
   ! Note: What is stored in qmat are NOT the Mulliken charges, but 
   !       the electron population from Mulliken analysis.
   !
   scf_mchg(1:qmmm_struct%nquant_nlink) =         &
         qzero( izp(1:qmmm_struct%nquant_nlink) ) &
         - qmat(1:qmmm_struct%nquant_nlink)

!! ===========================================================
!!               EVERYTHING BELOW THIS POINT 
!!                IS NOT USED IN AMBER (yet)
!! ===========================================================
!!   ! Total molecular charge
!!   do j = 1,nquant_nlink
!!      qtot = qtot+qmat(j)
!!   end do
!!
!!   ! ===============
!!   ! Electric dipole
!!   ! ===============
!!
!!   ! Dipole vector components
!!   do i = 1,3
!!      dipol(i) = 0.0d0
!!
!!      do j = 1,nquant_nlink
!!         izpj = izp(j)
!!         qhelp = qzero(izpj) - qmat(j)
!!         dipol(i) = dipol(i) + qhelp*qm_coords(i,j)*A_TO_BOHRS
!!      end do
!!      dipol(i) = dipol(i)*2.541765d0 !to Debye?
!!   end do
!!
!!   ! Norm of dipole vector
!!   dipabs = 0.0d0
!!   do i = 1,3
!!      dipabs = dipabs + dipol(i)**2
!!   end do
!!   dipabs = sqrt(dipabs)

end subroutine mulliken

!===============================================================================

subroutine qm2_dftb_dsymm(NDIM, MDIM, ADIM)

!! Just calls the LAPACK routine DSYMM
!! to calculate the density (or afoo) matrix.
!! The result is stored in the
!! ks_struct%density(MDIM,MDIM) matrix.

   use qm2_dftb_module, only: ks_struct

   implicit none

   integer, intent(in) :: NDIM
   integer, intent(in) :: MDIM
   integer, intent(in) :: ADIM

!!---------------------------------------------------
!! This is just a matrix multiplication!!
!! After this call, 
!! afoo(NDIM, adim) = overl(MDIM,NDIM) * a(MDIM,adim)
!!---------------------------------------------------
!!
!!      afoo = S x C
!! So,
!!      afoo x C = S x C x C = S x D'
!!
!! Where D' = C x C would be the density matrix IF the 
!! basis where ortogonal. (The REAL density matrix, in this
!! nonorthogonal basis is D = S x D' x S.)
!!
!!----------------------------------------------------

   call dsymm('L','U',NDIM,adim,1.0d0,ks_struct%overl,MDIM,ks_struct%a,MDIM,0.0d0,ks_struct%density,NDIM)

   return

end subroutine qm2_dftb_dsymm

!===============================================================================

subroutine qm2_dftb_qmulli(NDIM, MDIM, adim, calc_dens)
!!
!! Mulliken analysis:
!! ==================
!!
!! This subroutine builds the 'ks_struct%qmulli(NDIM)' 
!! vector, which is a Mulliken population per orbital.
!!
!! 'qmulli(n)' will be the charge localized on ATOMIC orbital 'n',
!! where the 'n' are organized by atom.
!!
!! OBS: If using the 'old' form, "density" is a misleading name.
!!      It is NOT the density matrix. Rather, it is the 'afoo'
!!      matrix fromthe original DFTB code:
!!                density ( = afoo ) = S x C
!!

   use qm2_dftb_module, only: ks_struct

   implicit none

!! Passed in:
   integer, intent(in) :: NDIM      ! Total number of orbitals
   integer, intent(in) :: MDIM      ! Maximum number of orbitals
   integer, intent(in) :: adim      ! Number of occupied orbitals
   logical, intent(in) :: calc_dens ! Do we have density matrix? 

!! Locals:
   integer :: i, n
   _REAL_  :: occi

   ks_struct%qmulli(1:NDIM) = 0.0d0

   if (calc_dens) then
      !   Use this if we have a density matrix in storage.
      !   (DFTB by default does not calculate it.)
      !   density = C x C, the real density matrix 
      !   The other option is actually faster. 
      !   Probably because it uses LAPACK, and because
      !   my routine to calcualate the dens matrix is too slow.

      call qm2_dftb_dens_matrix(adim, ndim)         ! Build density matrix      
      do i = 1, NDIM
         do n = 1, NDIM
            ks_struct%qmulli(i) = ks_struct%qmulli(i) + ks_struct%density(i,n) * ks_struct%overl(i,n)
         end do
      end do

   else
      !== Original version:
      !   (No density matrix, but much faster)
      !   Notice that, although it is called 'density', this is NOT the real
      !   density matrix, but the 'afoo' matrix from the original DFTB implementation. 
      !   density ( = afoo) = S x C
      !   'n' loops through atomic orbitals
      !   'i' loops through occupied (molecular) KS orbitals

      call qm2_dftb_dsymm(NDIM, MDIM, adim)         ! Builds 'density' (actually 'afoo = S x C')
      do i = 1,adim
         occi = ks_struct%occ(i)
         do n = 1,NDIM
            ks_struct%qmulli(n) = ks_struct%qmulli(n) + occi * ks_struct%density(n,i) * ks_struct%a(n,i)
         end do
      end do
      !== End old version
   end if

   return

end subroutine qm2_dftb_qmulli

!===============================================================================

subroutine qm2_dftb_mull_pop(nquant_nlink,NDIM,izp,lmax,qmat)
!! 
!! MULLIKEN POPULATIONS PER ATOM
!! =============================
!! Sum up the charge contributions from each atomic
!! orbital into the atoms, to find out the electron 
!! population for each atom.
!! 
!! The atomic Mulliken populations are returned
!! in 'qmulli(NATOMS)'.
!!

   use qm2_dftb_module, only: ks_struct
   implicit none

!! Passed in
   integer, intent(in)  :: nquant_nlink         ! Number of atoms
   integer, intent(in)  :: NDIM       ! Total number of orbitals
   integer, intent(in)  :: izp(*)     ! Atom types
   integer, intent(in)  :: lmax(*)    ! Maximum 'l' (per atom type)
   _REAL_ , intent(out) :: qmat(*)    ! Atomic electron populations

!! Locals
   integer :: i, j, lj, mj, jofn, ljmax, indj, mjmax
   _REAL_  :: qtot, qhelp



   ! Loop over all atoms
   do j = 1,nquant_nlink
      ! For each atom, loop over shells (s,p,d)
      ! lmax=1 for H, =2 for O,C,N
      qtot  = 0.0d0
      ljmax = lmax(izp(j))
      indj  = ks_struct%ind(j)

      do lj = 1,ljmax
         ! For each shell, loop over magnetic quantum numbers
         ! (s, px,py,pz ...)
         ! qhelp accumulates the charge contribution from each AO

         jofn  = indj + (lj-1) * (lj-1)
         qhelp = 0.0d0
         mjmax = 2*lj-1

         do mj = 1,mjmax
            qhelp = qhelp + ks_struct%qmulli(jofn+mj)
         end do

         ! qtot accumulates the charge contribution from each shell
         qtot = qtot + qhelp
      end do
      qmat(j) = qtot
   end do


end subroutine qm2_dftb_mull_pop

!===============================================================================

subroutine qm2_dftb_dens_matrix(adim, ndim)
!! Builds density matrix into ks_struct%density
   !ks_struct%density(i,j) = ks_struct%density(i,j) 
   !                       + ks_struct%occ(k) * ks_struct%a(i,k) * ks_struct%a(j,k)

   use qm2_dftb_module, only: ks_struct

   implicit none

!Passed in
   integer, intent(in) :: adim, ndim

!! Locals
   integer :: i,j,k
   _REAL_  :: temp, tempi

   ks_struct%density(1:NDIM,1:NDIM) = 0.0d0

   do k = 1, adim
      temp = ks_struct%occ(k)
      do i = 1, NDIM
         tempi = temp * ks_struct%a(i,k)
         do j = 1, NDIM
            ks_struct%density(i,j) = ks_struct%density(i,j) + tempi * ks_struct%a(j,k)
         end do
      end do
   end do

end subroutine qm2_dftb_dens_matrix

!===============================================================================
