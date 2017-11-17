! <compile=optimized>
!-*- mode: f90; coding: iso-8859-15; -*-

!
!     Loads parameters for DFTB calculations
!     **************************************
!
#include "copyright.h"
#include "dprec.fh"

subroutine qm2_dftb_load_params(qm2_struct, qmmm_struct)

!In parallel all threads call this routine

   use qmmm_module, only : qmmm_nml, qm2_structure, &
                           qmmm_mpi
   use ElementOrbitalIndex, only : elementSymbol
   use constants, only: EV_TO_KCAL, AU_TO_EV, AU_TO_KCAL, A_TO_BOHRS
   use qm2_dftb_module, only: NDIM, LDIM, MDIM, &     ! Fixed parameters
         MAX_BRD_ITER, IMATSZ,&                       ! Broyden mixing
         NNDIM, MAXSIZ, mol, lmax, mcharge, izp_str, &
         sktab, log_racc, dacc, dispfile, ks_struct, fermi_str, DFTB_3rd_order_str
   use qmmm_struct_module, only : qmmm_struct_type
   
   implicit none
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   type(qm2_structure),intent(inout) :: qm2_struct


   !Locals
   logical do_repeat
   integer :: pair_count

   _REAL_ :: yhlp ! Used to set sktab%slkcutoff
   _REAL_ :: racc ! Used for calculating machine precision limit.

   integer k, i, j

   !===
   character(LEN=2) :: atom_uc, atom_lc
   character(LEN=2) :: lowercase_atom
   !===

   character(len=1024) :: skroot
   character(*), parameter :: skext =".skf"

   ! Dispersion
   character(len=255) ::  disp_file

   ! CM3 charges
   character(len=255) ::  cm3_file

if (qmmm_mpi%commqmmm_master) then
   !===================
   ! Parameters
   !===================

   ! Broyden
   MAX_BRD_ITER = qmmm_nml%dftb_maxiter
   IMATSZ       = qmmm_nml%dftb_maxiter
   NNDIM  = qmmm_struct%nquant_nlink
   !MDIM   = 400! 9 * NNDIM
   MDIM   = LDIM * NNDIM
   MAXSIZ = MDIM

   !===================
   ! Memory Allocation
   !===================

   call qm2_dftb_allocate(qmmm_struct)

   !=======================
   !     Initializaton
   !=======================

   call getenv('AMBERHOME',skroot)
   skroot = TRIM(skroot) // '/dat/slko/'

   ! izp (atom types)
   izp_str%izp(1:NNDIM)=0

   ! S-K files
   sktab%skfiles(1:qmmm_struct%qm_ntypes,1:qmmm_struct%qm_ntypes)=""

   ! Gets this machine accuracy --> racc
   racc = 1.0d0
   DO
      if((1.0d0+racc) <= 1.0d0) EXIT
      racc = racc*0.5d0
   END DO
   racc = 2.0d0*racc
   log_racc = log(racc)
   dacc = 4.0d0*racc

   !=========================
   !     Number of Atoms 
   !=========================
   ! qmmm_struct%nquant_nlink


   !======================================
   !     Sets the control paramenters
   !======================================
   izp_str%nel     = 0         ! Total charge

   !=====================================================
   !     Find the number of (independent) atom types
   !=====================================================
   qmmm_struct%qm_ntypes = 0
   do i = 1,qmmm_struct%nquant_nlink 
      do_repeat = .false.
      do j = 1, i
         if ( i /= j .and. &
               qmmm_struct%iqm_atomic_numbers(i) == &
                     qmmm_struct%iqm_atomic_numbers(j)) then
            do_repeat = .true.
            izp_str%izp(i) =  izp_str%izp(j)
         end if
      end do
      if (.not.do_repeat) then
         qmmm_struct%qm_ntypes = qmmm_struct%qm_ntypes+1
         !=========================================
         !     Set the maximum l for each atom
         !=========================================
         ! There is NO check to wether an atom is supported any longer.
         ! The reason is that the user can just add any parameter
         ! sets he sees fit.
         !
         ! CURRENTLY: H,C,N,O,S,Zn
         !
         if (qmmm_struct%iqm_atomic_numbers(i) .le. 2) then  
            ! 1st row
            lmax(qmmm_struct%qm_ntypes) = 1
         else if ( qmmm_struct%iqm_atomic_numbers(i) .ge. 3 .and. &
                   qmmm_struct%iqm_atomic_numbers(i) .le. 10      ) then
            ! 2nd row 
            lmax(qmmm_struct%qm_ntypes) = 2
         else if ( qmmm_struct%iqm_atomic_numbers(i) .gt. 10 ) then
            ! 3rd row on
            lmax(qmmm_struct%qm_ntypes) = 3
         end if
         izp_str%izp(i) = qmmm_struct%qm_ntypes
      end if
   end do

   ! Fill mol%atyp with atomic symbols
   do i = 1, qmmm_struct%nquant_nlink
      mol%atyp( izp_str%izp(i) ) = elementSymbol( qmmm_struct%iqm_atomic_numbers(i) )
   end do

   ! Calculate the number of orbitals
   qm2_struct%norbs=0
   do i=1, qmmm_struct%nquant_nlink
      if (lmax(izp_str%izp(i)).eq.1) then
         qm2_struct%norbs=qm2_struct%norbs+1
      else if (lmax(izp_str%izp(i)).eq.2) then
         qm2_struct%norbs=qm2_struct%norbs+4
      else 
         qm2_struct%norbs=qm2_struct%norbs+9
      end if        
   end do

   !===============================================
   !     Fill the Slater-Kosher integral tables
   !===============================================
   ! Note: There are (qmmm_struct%qm_ntypes)**2 tables to fill.
   !
   ! Note: interaction between atom1-atom2 is DIFFERENT
   !       than atom2-atom1, so for each pair there's
   !       really 2 files to load.

   do i = 1,qmmm_struct%nquant_nlink 
      atom_uc=elementSymbol(qmmm_struct%iqm_atomic_numbers(i))
      atom_lc=lowercase_atom(atom_uc)
      sktab%latyp( izp_str%izp(i)) = TRIM(atom_lc)
   end do
   write(6,'(a,i4)') ' DFTB: Number of atom types = ', qmmm_struct%qm_ntypes
   write(6,*)
   write(6, '(a)') " Parameter files:"
   write(6,'( 5X,"TYP (AT)",2X,"TYP (AT)",5X, "SK integral FILE")')
   pair_count = 0
   do i = 1, qmmm_struct%qm_ntypes
      do j = 1, qmmm_struct%qm_ntypes
         pair_count = pair_count + 1
         ! New file naming convention from www.dftb.org
         sktab%skfiles(i,j)=TRIM(skroot)//TRIM(mol%atyp(i))//"-"//TRIM(mol%atyp(j))//skext 
         write (6,'("|",i3,1x,i2,2x,"(",a2,")",2x,i2,2x,"(",a2,")",5x,A)') &
               pair_count, &
               i, &         
               mol%atyp(i), &
               j, & 
               mol%atyp(j), &
               TRIM(sktab%skfiles(i,j))
         call qm2_dftb_check_slko_file(sktab%skfiles(i,j))
      end do
   end do

   call gettab(qmmm_struct%qm_ntypes)!,sktab%skfiles)

   ! Set the SK cutoff
   sktab%slkcutoff=0.0d0
   do i = 1,qmmm_struct%qm_ntypes
      do j = 1,qmmm_struct%qm_ntypes
         yhlp=sktab%sr(i,j)*sktab%dimens(i,j)+0.3d0
         if(yhlp > sktab%slkcutoff)then
            sktab%slkcutoff=yhlp
         endif
      end do
   end do

!!=================================
!!    Read the dispersion file
!!=================================
   if (qmmm_nml%dftb_disper == 1) then
      disp_file = TRIM(skroot)//"DISPERSION.INP_ONCHSP"
      call dispersionread(qmmm_struct, qmmm_struct%nquant_nlink,qmmm_struct%qm_ntypes, &
            izp_str%izp, &
            disp_file)
   end if

!!=================================
!!    Third order
!!=================================
   if (DFTB_3rd_order_str%do_3rd_order) then
      call qm2_dftb_read_3rd_order(qmmm_struct, qmmm_struct%nquant_nlink,qmmm_struct%qm_ntypes, & 
            izp_str%izp,skroot)
   end if


!!=================================
!!    CM3 Charges
!!=================================
   if (qmmm_nml%dftb_chg == 1) then

      write(6,*) "READING CM3 CHARGES PARAMETERS..."

      cm3_file = TRIM(skroot)//"CM3_PARAMETERS.DAT"
      call read_cm3(qmmm_struct%nquant_nlink, qmmm_struct%qm_ntypes, izp_str%izp, cm3_file)
   end if

!!=================================
!! Fermi Distribution
!!=================================
   fermi_str%telec = qmmm_nml%dftb_telec
   fermi_str%telec_step = qmmm_nml%dftb_telec_step


!! Initialize the charges

   ! Initial electron population per atom (# valence electrons.)
   mol%qmat(1:qmmm_struct%nquant_nlink) = mcharge%qzero( izp_str%izp(1:qmmm_struct%nquant_nlink) )

   ! Mulliken Charges
   !   WRITE(22,*)"DBG PRT qm2_dftb_load_params",qmmm_struct%nquant_nlink,ASSOCIATED(qm2_struct%scf_mchg)
   IF ( .NOT. ASSOCIATED( qm2_struct%scf_mchg ) ) THEN
      WRITE(6,'(A)')"sqm/qm2_dftb_load_params: scf_mchg is not associated"
      WRITE(6,'(A)')"Cannot currently use DFTB"
      STOP
   END IF
   qm2_struct%scf_mchg(1:qmmm_struct%nquant_nlink) = 0.0d0

   ! Calculates the number of electrons in the system.
   izp_str%nel = -qmmm_nml%qmcharge
   do i = 1,qmmm_struct%nquant_nlink
      izp_str%nel =  izp_str%nel+mcharge%qzero( izp_str%izp(i))
   end do

   ! Indices for H and S matrices
   ks_struct%ind(1) = 0
   do j = 1, qmmm_struct%nquant_nlink
      ks_struct%ind(j+1) = ks_struct%ind(j) + lmax( izp_str%izp(j) )**2
   end do

   ! actual dimension of matrix
   ndim = ks_struct%ind(qmmm_struct%nquant_nlink+1)

   ! Allocate matrices which size depend on slko parameters
   call qm2_dftb_allocate_slko_depn

   ! Static memory limitation
   if (ndim > MDIM) then
      write(6,*) ' eglcao: ndim > ', MDIM
      stop
   endif
end if

   return
end subroutine qm2_dftb_load_params


subroutine qm2_dftb_check_slko_file(sk_file)

   implicit none
   
   character (len=*) :: sk_file
   logical :: found

   found = .false.
   inquire(FILE=TRIM(sk_file), EXIST=found)
   
   if (.not.found) then
      write(6,*)"****************************************************"
      write(6,*)"*     !! A FILE NEEDED BY DFTB WAS NOT FOUND !!    *"
      write(6,*)"****************************************************"
      write(6,*)
      write(6,*)" Missing file:"
      write(6,*)" ",TRIM(sk_file)
      write(6,*)
      write(6,*)" Sander could not find the file containing the "
      write(6,*)" Slater-Koster integral tables needed. This means "
      write(6,*)" either that the file was placed in the wrong "
      write(6,*)" directory, is using a different naming scheme, "
      write(6,*)" or that you just don't have it." 
      write(6,*)
      write(6,*)" These files are supposd to be located in the "
      write(6,*)" $(AMBERHOME)/dat/slko/"
      write(6,*)" directory, and be named as <Atom1>-<Atom2>.skf,"
      write(6,*)" where <Atom1> and <Atom2> are the atomic symbols."
      write(6,*)
      write(6,*)" Note that the integral table files needed for"
      write(6,*)" a DFTB calculation are not distributed with Amber."
      write(6,*)" To obtain those files, you must point your browser to"
      write(6,*)" http://www.dftb.org and follow the instructions"
      write(6,*)" to obtain the parameter files."
      write(6,*)
      write(6,*)" Also, you must make sure that the parameters for the"
      write(6,*)" atoms you are interested in exist. Due to the continuous "
      write(6,*)" development of new parameters, SANDER NO LONGER DOES"
      write(6,*)" THIS CHECK. Rather, Sander only checks for the "
      write(6,*)" presence of the parameter file."
      write(6,*)

      call sander_bomb("qm2_dftb_check_slko_file <qm2_dftb_load_params.f>", &
                       "File not found.", &
                       "Exiting.")
   endif

   return
end subroutine qm2_dftb_check_slko_file


! #############################################################################
! #                                                                           #
! #                    MEMORY ALLOCATION / DEALLOCATION                       #
! #                                                                           #
! #############################################################################

! Allocation error ('REQUIRE...')
#include "assert.fh"


subroutine qm2_dftb_allocate_slko_depn
   use qm2_dftb_module
   implicit none

   integer :: ier = 0 ! Allocation status

   allocate(ks_struct%density(NDIM,NDIM) , stat=ier)
   REQUIRE(ier == 0)


end subroutine qm2_dftb_allocate_slko_depn


subroutine qm2_dftb_allocate(qmmm_struct)
!!============================================
!!           Memory Allocation
!!============================================

   use qm2_dftb_module
   use qmmm_struct_module, only : qmmm_struct_type

   implicit none

   integer :: ier=0 ! Allocation status
   type(qmmm_struct_type), intent(inout) :: qmmm_struct

!! Begin Memory Allocations
!!==========================

   ! ----------------------------------------
   ! Arrays that do NOT belong to a structure
   ! ----------------------------------------
   allocate( dummy(LDIM,LDIM), stat=ier ) ! Scratch space used by slkode <qm2_dftb_slkode.f>
   REQUIRE(ier == 0)

   allocate( espin(qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier == 0)

   allocate( lmax(qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier ==0)

   allocate(scr_space(3*MDIM), stat=ier)
   REQUIRE(ier ==0)


   ! --------------------------------
   ! Arrays that belong to structures
   ! --------------------------------
   
   ! dispertmp_structure
   allocate(dispertmp%C6(NNDIM,NNDIM), stat=ier )
   REQUIRE(ier ==0)

   allocate(dispertmp%Rvdw(NNDIM,NNDIM), stat=ier )
   REQUIRE(ier ==0)

   ! dispfile_structure
   allocate(dispfile%Ni0(qmmm_struct%qm_ntypes), stat=ier )
   REQUIRE(ier ==0)

   allocate(dispfile%h1(qmmm_struct%qm_ntypes,4), stat=ier )
   REQUIRE(ier ==0)

   allocate(dispfile%h2(qmmm_struct%qm_ntypes,4), stat=ier )
   REQUIRE(ier ==0)

   allocate(dispfile%nei(NNDIM), stat=ier )
   REQUIRE(ier ==0)

   allocate(dispfile%Ni(NNDIM), stat=ier )
   REQUIRE(ier ==0)

   allocate(dispfile%hh1(NNDIM), stat=ier )
   REQUIRE(ier ==0)

   allocate(dispfile%hh2(NNDIM), stat=ier )
   REQUIRE(ier ==0)

   allocate(dispfile%R0vdw(qmmm_struct%qm_ntypes,4), stat=ier )
   REQUIRE(ier ==0)

   ! DFTB_3rd_order_str
   allocate(DFTB_3rd_order_str%Hubbard_deriv(100), stat=ier )
   REQUIRE(ier ==0)
   
   ! mol_structure
   allocate(mol%qmat(NNDIM), stat=ier)
   REQUIRE(ier == 0)

   allocate(mol%atyp(qmmm_struct%qm_ntypes), stat=ier )
   REQUIRE(ier ==0)

   ! sktab_structure
   allocate(sktab%skfiles(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier == 0)

   allocate(sktab%latyp(qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier == 0)

   allocate(sktab%sr(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes),stat=ier)
   REQUIRE(ier == 0)
   
   allocate(sktab%dimens(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes),stat=ier)
   REQUIRE(ier == 0)

   allocate(sktab%skself(3,qmmm_struct%qm_ntypes),stat=ier)
   REQUIRE(ier == 0)

   allocate(sktab%skhtab(10,MAXTAB,qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes),stat=ier)
   REQUIRE(ier == 0)

   allocate(sktab%skstab(10,MAXTAB,qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes),stat=ier)
   REQUIRE(ier == 0)

   ! spltab_structure
   allocate(spltab%numint(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier == 0)
   
   allocate(spltab%cutoff(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier == 0)
   
   allocate(spltab%efkt(3,qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier == 0)
   
   allocate(spltab%xr(2,MAXINT,qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier == 0)
   
   allocate(spltab%coeff(6,MAXINT,qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier == 0)

   ! mcharge_structure
   allocate(mcharge%qzero(qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier == 0)

   allocate(mcharge%uhubb(qmmm_struct%qm_ntypes), stat=ier)
   REQUIRE(ier == 0)

   ! izp_structure
   allocate(izp_str%izp(NNDIM), stat=ier)
   REQUIRE(ier == 0)

   ! ks_dftb_structure
   allocate(ks_struct%ind(NNDIM+1) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%hgrad(3,NNDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%au(LDIM,LDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%bu(LDIM,LDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%auh(LDIM,LDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%buh(LDIM,LDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%a(MDIM,MDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%b(MDIM,MDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%gammamat(NNDIM,NNDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%derivx(NNDIM,NNDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%derivy(NNDIM,NNDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%derivz(NNDIM,NNDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%ev(MDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%occ(MDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%qmulli(MDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%qmold(NNDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%hamil(MDIM,MDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%overl(MDIM,MDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%shift(NNDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%shiftE(NNDIM) , stat=ier)
   REQUIRE(ier == 0)


   !==
   
   allocate(ks_struct%scr1(MDIM,MDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%scr2(MDIM,MDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(ks_struct%xtrans(MDIM,MDIM) , stat=ier)
   REQUIRE(ier == 0)

   !==
   
   ! broyden_structure
   allocate(brd_struct%f(MAXSIZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%ui(MAXSIZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%vti(MAXSIZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%t1(MAXSIZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%dumvi(MAXSIZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%df(MAXSIZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%a(IMATSZ,IMATSZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%b(IMATSZ,IMATSZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%cm(IMATSZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%d(IMATSZ,IMATSZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%w(IMATSZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%unit31(MAXSIZ,2) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%unit32(MAXSIZ,2,MAX_BRD_ITER) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%td(IMATSZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%ad(IMATSZ) , stat=ier)
   REQUIRE(ier == 0)

   allocate(brd_struct%bd(IMATSZ) , stat=ier)
   REQUIRE(ier == 0)

   ! cm3_structure
   allocate(cm3%qcm3(NNDIM), stat=ier)
   REQUIRE(ier == 0)

   allocate(cm3%d(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes) , stat=ier)
   REQUIRE(ier == 0)

   allocate(cm3%c(qmmm_struct%qm_ntypes,qmmm_struct%qm_ntypes) , stat=ier)
   REQUIRE(ier == 0)

   allocate(cm3%t(NNDIM,NNDIM) , stat=ier)
   REQUIRE(ier == 0)

   allocate(cm3%b(NNDIM,NNDIM) , stat=ier)
   REQUIRE(ier == 0)

!!==========================
!! End Memory Allocations
!!==========================
end subroutine qm2_dftb_allocate


!! ##############################################################

subroutine qm2_dftb_deallocate
!!============================================
!!           Memory Deallocation
!!============================================
!! Called from deallocate_qmmm < qmmm_module.f >

   use qm2_dftb_module

   implicit none

   integer :: ier=0 ! Dellocation status

   ! ----------------------------------------
   ! Arrays that do NOT belong to a structure
   ! ----------------------------------------
   deallocate( dummy, stat=ier ) ! Scratch space used by slkode <qm2_dftb_slkode.f>
   REQUIRE(ier == 0)

   deallocate( espin, stat=ier)
   REQUIRE(ier == 0)

   deallocate( lmax, stat=ier)
   REQUIRE(ier ==0)

   deallocate(scr_space, stat=ier)
   REQUIRE(ier ==0)


   ! --------------------------------
   ! Arrays that belong to structures
   ! --------------------------------
   
   ! dispertmp_structure
   deallocate(dispertmp%C6, stat=ier )
   REQUIRE(ier ==0)

   deallocate(dispertmp%Rvdw, stat=ier )
   REQUIRE(ier ==0)

   ! dispfile_structure
   deallocate(dispfile%Ni0, stat=ier )
   REQUIRE(ier ==0)

   deallocate(dispfile%h1, stat=ier )
   REQUIRE(ier ==0)

   deallocate(dispfile%h2, stat=ier )
   REQUIRE(ier ==0)

   deallocate(dispfile%nei, stat=ier )
   REQUIRE(ier ==0)

   deallocate(dispfile%Ni, stat=ier )
   REQUIRE(ier ==0)

   deallocate(dispfile%hh1, stat=ier )
   REQUIRE(ier ==0)

   deallocate(dispfile%hh2, stat=ier )
   REQUIRE(ier ==0)

   deallocate(dispfile%R0vdw, stat=ier )
   REQUIRE(ier ==0)


   
   ! mol_structure
   deallocate(mol%qmat, stat=ier)
   REQUIRE(ier == 0)

   deallocate(mol%atyp, stat=ier )
   REQUIRE(ier ==0)

   ! sktab_structure
   deallocate(sktab%skfiles, stat=ier)
   REQUIRE(ier == 0)

   deallocate(sktab%latyp, stat=ier)
   REQUIRE(ier == 0)

   deallocate(sktab%sr,stat=ier)
   REQUIRE(ier == 0)
   
   deallocate(sktab%dimens,stat=ier)
   REQUIRE(ier == 0)

   deallocate(sktab%skself,stat=ier)
   REQUIRE(ier == 0)

   deallocate(sktab%skhtab,stat=ier)
   REQUIRE(ier == 0)

   deallocate(sktab%skstab,stat=ier)
   REQUIRE(ier == 0)

   ! spltab_structure
   deallocate(spltab%numint, stat=ier)
   REQUIRE(ier == 0)
   
   deallocate(spltab%cutoff, stat=ier)
   REQUIRE(ier == 0)
   
   deallocate(spltab%efkt, stat=ier)
   REQUIRE(ier == 0)
   
   deallocate(spltab%xr, stat=ier)
   REQUIRE(ier == 0)
   
   deallocate(spltab%coeff, stat=ier)
   REQUIRE(ier == 0)

   ! mcharge_structure
   deallocate(mcharge%qzero, stat=ier)
   REQUIRE(ier == 0)

   deallocate(mcharge%uhubb, stat=ier)
   REQUIRE(ier == 0)

   ! ks_dftb_structure
   deallocate(ks_struct%ind, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%hgrad, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%au, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%bu, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%auh, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%buh, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%a, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%b, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%gammamat, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%derivx, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%derivy, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%derivz, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%ev, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%occ, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%qmulli, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%qmold, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%hamil, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%overl, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%density, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%shift, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%shiftE, stat=ier)
   REQUIRE(ier == 0)

   !==

   deallocate(ks_struct%scr1, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%scr2, stat=ier)
   REQUIRE(ier == 0)

   deallocate(ks_struct%xtrans, stat=ier)
   REQUIRE(ier == 0)

   !==
   
   ! broyden_structure
   deallocate(brd_struct%f, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%ui, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%vti, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%t1, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%dumvi, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%df, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%a, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%b, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%cm, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%d, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%w, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%unit31, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%unit32, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%td, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%ad, stat=ier)
   REQUIRE(ier == 0)

   deallocate(brd_struct%bd, stat=ier)
   REQUIRE(ier == 0)

   ! cm3_structure
   deallocate(cm3%qcm3, stat=ier)
   REQUIRE(ier == 0)

   deallocate(cm3%d, stat=ier)
   REQUIRE(ier == 0)

   deallocate(cm3%c, stat=ier)
   REQUIRE(ier == 0)

   deallocate(cm3%t, stat=ier)
   REQUIRE(ier == 0)

   deallocate(cm3%b, stat=ier)
   REQUIRE(ier == 0)

!!==========================
!! End Memory Deallocations
!!==========================
end subroutine qm2_dftb_deallocate

character(LEN=2) function lowercase_atom(atom)
! Returns the lowercase of string 'str'
  implicit none
  character(LEN=2), intent(in) :: atom
  character(LEN=1) :: lowercase_char, char_uc, char_lc
  integer i

  !write(6,*) " Converting to lowercase. LEN=", len(atom)
  !write(6,*) " atom      = ", atom

  !write(6,*) " atom(1:1) = ", "'",atom(1:1),"'"
  
  !write(6,*) " atom(2:2) = ", "'",atom(2:2),"'"
  lowercase_atom=''
  do i=1, len(atom)
    char_uc = atom(i:i)
    char_lc = lowercase_char(char_uc)
    !write(6,*) " --> Lowercase: ", "'",char_lc,"'"
    lowercase_atom = TRIM(lowercase_atom) // char_lc
    !write(6,*) "'",lowercase_atom,"'"
  end do
  return

end function lowercase_atom

character(len=1) function lowercase_char(char)
! Returns the lowercase character

  implicit none

  !! Passed in
  character (LEN=1), intent (in) :: char

  !! Locals
  !write(6,*) "Routine lowercase char :: CHAR = ", "'",char,"'"
  lowercase_char = char
  if (char == 'A') lowercase_char = 'a'
  if (char == 'B') lowercase_char = 'b'
  if (char == 'C') lowercase_char = 'c'
  if (char == 'D') lowercase_char = 'd'
  if (char == 'E') lowercase_char = 'e'
  if (char == 'F') lowercase_char = 'f'
  if (char == 'G') lowercase_char = 'g'
  if (char == 'H') lowercase_char = 'h'
  if (char == 'I') lowercase_char = 'i'
  if (char == 'J') lowercase_char = 'j'
  if (char == 'K') lowercase_char = 'k'
  if (char == 'L') lowercase_char = 'l'
  if (char == 'M') lowercase_char = 'm'
  if (char == 'N') lowercase_char = 'n'
  if (char == 'O') lowercase_char = 'o'
  if (char == 'P') lowercase_char = 'p'
  if (char == 'Q') lowercase_char = 'q'
  if (char == 'R') lowercase_char = 'r'
  if (char == 'S') lowercase_char = 's'
  if (char == 'T') lowercase_char = 't'
  if (char == 'U') lowercase_char = 'u'
  if (char == 'V') lowercase_char = 'v'
  if (char == 'W') lowercase_char = 'w'
  if (char == 'X') lowercase_char = 'x'
  if (char == 'Y') lowercase_char = 'y'
  if (char == 'Z') lowercase_char = 'z'
  !write(6,*) "Lowecase char = ", "'",lowercase_char,"'"
  return
end function lowercase_char

