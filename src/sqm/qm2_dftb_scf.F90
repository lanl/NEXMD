! <compile=optimized>
#include "copyright.h"
#include "dprec.fh"
#include "def_time.h"
subroutine qm2_dftb_scf(escf, elec_eng,enuclr_qmqm,scf_mchg)

   use qmmm_module, only : qmmm_nml, qmmm_mpi, qmmm_struct
   use qm2_dftb_module, only : mcharge, izp_str, fermi_str, mol, disper, espin,dftb_3rd_order_str
   use constants, only : AU_TO_EV, AU_TO_KCAL

   implicit none

!For irespa
#include "md.h" 

!Passed in
   _REAL_, intent(out) :: elec_eng, enuclr_qmqm, escf
   _REAL_, intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)

!Local
   _REAL_ :: geseatom, total_e, ethird
   integer :: outer_scf_count, inner_scf_count, i
   logical :: scf_convgd

   total_e      = 0.0d0
   
   scf_convgd=.false.
   outer_scf_count = 0
   inner_scf_count = 0

   !---------------------
   !BEGIN OUTER SCF LOOP
   !---------------------
   do_scf_outer: do while ( (.not.scf_convgd) .and. (outer_scf_count < qmmm_nml%itrmax) )
     
      call eglcao(qmmm_struct%qm_coords,total_e,elec_eng,ethird,enuclr_qmqm,&
                  inner_scf_count, outer_scf_count, scf_convgd,mol%qmat,scf_mchg )

      if ( .not.scf_convgd ) then 
        ! Take steps to improve convergence

        ! Reset Broyden mixing by restarting the scf process in eglcao.
        if (qmmm_nml%verbosity > 0 .and. qmmm_mpi%commqmmm_master) then
           write(6,'(" QMMM SCC-DFTB: SCC-DFTB FOR STEP ",i5," DID NOT CONVERGE AFTER ",i3," cycles.")') &
                 irespa,outer_scf_count
           if (outer_scf_count < qmmm_nml%itrmax) write(6,'(" QMMM SCC-DFTB: Resetting Broyden mixing.")')
        end if

        ! In the first time, maybe the initial charges were way too wrong
        ! for this iteration, re-set qmat:
        if (outer_scf_count == qmmm_nml%dftb_maxiter .and. qmmm_mpi%commqmmm_master) then
           if (qmmm_nml%verbosity > 0) &
                 write(6,'(" QMMM SCC-DFTB: Resetting initial charges.")')
           mol%qmat(1:qmmm_struct%nquant_nlink) = mcharge%qzero( izp_str%izp(1:qmmm_struct%nquant_nlink) )
        end if

        ! Increase telec...
        ! WARNING: This is optional, and should be used only with great care. I
        !          have no idea what effects could happen.
        if (fermi_str%telec_step > 0.0d0 .and. qmmm_mpi%commqmmm_master) then
           ! Increase tlec and give a warning
           fermi_str%telec = fermi_str%telec + fermi_str%telec_step
           if (qmmm_nml%verbosity > 0) &
                 write(6,'(" QMMM SCC-DFTB: Increasing the electronic temperature of step ",i5," to ",f10.3," K")')&
                 irespa, fermi_str%telec
        end if
      end if

   end do do_scf_outer ! do while ( (.not.scf_convgd) .and. (outer_scf_count < qmmm_nml%itrmax) )

   !---------------------
   !END OF OUTER SCF LOOP
   !---------------------

   ! If we had to increase telec, now we restore it to the original value.
   if ( fermi_str%telec /= qmmm_nml%dftb_telec ) then
     if (qmmm_mpi%commqmmm_master) then
       write(6,'(" QMMM SCC-DFTB: **WARNING** : The energy for this step was calculated with a higher")')
       write(6,'(" QMMM SCC-DFTB:               elctronic temperature. You should check the results.")')
       write(6,'(" QMMM SCC-DFTB:               STEP=",i5,5X," TELEC=",f10.3)') irespa,fermi_str%telec
     end if
     ! NOT IMPLEMENTED YET:
     ! CALL PRINT_DFTB_CHARGES_AND_DIPOLE_FOR_CHECKING
     fermi_str%telec = qmmm_nml%dftb_telec
     if (qmmm_mpi%commqmmm_master) &
       write(6,'(" QMMM SCC-DFTB: Electronic temperature restored to ",f10.3," K.")') fermi_str%telec
   end if

   ! If it still didn't converge, something must be really wrong. Bomb the calculation.
   if ( outer_scf_count >= qmmm_nml%itrmax .or. .not.scf_convgd) call dftb_conv_failure("dylcao <qm2_dftb_main.f> : ", &
          "SCC Convergence failure - ITRMAX exceeded.", "Exiting")


   if (qmmm_nml%verbosity > 0 .and. qmmm_mpi%commqmmm_master) then
     if (scf_convgd) then
        write(6,'(" QMMM SCC-DFTB: SCC-DFTB for step ",i5," converged in ",i3," cycles.")') irespa,outer_scf_count
     else
        write(6,'(" QMMM SCC-DFTB: **WARNING** SCC-DFTB FOR STEP ",i5," DID NOT CONVERGE AFTER ",i3," cycles.")') &
              irespa,inner_scf_count
     end if
        write(6,*)
   end if


   if (qmmm_mpi%commqmmm_master) then

      ! Atomization energy
      geseatom=0.0d0
      do i=1,qmmm_struct%nquant_nlink
         geseatom = geseatom+espin(izp_str%izp(i))
      end do

     ! Binding Energy. Will be put into Amber's escf energy.
     ! (This is done here mainly to shift the zero of energy,
     !  so that the DFTB energy scale is closer to the MM energy scale.)
     ! Convert results to units used by Amber
     ! (The results from dylcao come in a.u.)
     escf = (total_e-geseatom) * AU_TO_KCAL
     elec_eng  = elec_eng * AU_TO_EV
     enuclr_qmqm = enuclr_qmqm * AU_TO_EV

     !==============================
     !     Prints DFTB results
     !==============================
     if (qmmm_nml%verbosity > 1) then
        write(6,'(" QMMM SCC-DFTB:")')
        write(6,'(" QMMM SCC-DFTB:    Atomization Energy (eV)       = ",f20.12)') geseatom  * AU_TO_EV
        write(6,'(" QMMM SCC-DFTB:    Electronic Energy  (eV)       = ",f20.12)') elec_eng    
        write(6,'(" QMMM SCC-DFTB:    Repulsive Energy   (eV)       = ",f20.12)') enuclr_qmqm 
        if (qmmm_nml%dftb_disper ==1 )&
              write(6,'(" QMMM SCC-DFTB:    Dispersion Energy  (eV)       = ",f20.12)') disper%edis * AU_TO_EV
        if (dftb_3rd_order_str%do_3rd_order)&
              write(6,'(" QMMM SCC-DFTB:    Third Order Energy (eV)       = ",f20.12)') ethird * AU_TO_EV
        write(6,'(" QMMM SCC-DFTB:    Total Energy       (eV)       = ",f20.12)') total_e*AU_TO_EV
        write(6,'(" QMMM SCC-DFTB:    SCF Energy         (eV)       = ",f20.12)') escf * AU_TO_EV / AU_TO_KCAL
        write(6,'(" QMMM SCC-DFTB:")')

        if ( qmmm_nml%verbosity > 3) then
           write(6,'(" QMMM SCC-DFTB:")')
           write(6,'(" QMMM SCC-DFTB:    Atomization Energy (a.u.)     = ",f20.12)') geseatom 
           write(6,'(" QMMM SCC-DFTB:    Electronic Energy  (a.u.)     = ",f20.12)') elec_eng / AU_TO_EV
           write(6,'(" QMMM SCC-DFTB:    Repulsive Energy   (a.u.)     = ",f20.12)') enuclr_qmqm / AU_TO_EV
           if (qmmm_nml%dftb_disper ==1 )&
                 write(6,'(" QMMM SCC-DFTB:    Dispersion Energy  (a.u.)     = ",f20.12)') disper%edis
           if (dftb_3rd_order_str%do_3rd_order)&
                 write(6,'(" QMMM SCC-DFTB:    Third Order Energy (a.u.)     = ",f20.12)') ethird
           write(6,'(" QMMM SCC-DFTB:    Total Energy       (a.u.)     = ",f20.12)') total_e
           write(6,'(" QMMM SCC-DFTB:    SCF Energy         (a.u.)     = ",f20.12)') escf / AU_TO_KCAL
           write(6,'(" QMMM SCC-DFTB:")')
        end if
     end if

   end if ! (qmmm_mpi%commqmmm_master)

   return
end subroutine qm2_dftb_scf

subroutine eglcao(qm_coords,total_e,elec_eng,ethird,enuclr_qmqm, &
      inner_scf_count, outer_scf_count, scc_converged,qmat,scf_mchg)

! SUBROUTINE EGLCAO
! =================
!
! Copyright 1997 by Peter Blaudeck, Dirk Porezag, Michael Haugk,
! Joachim Elsner
!
! The original routine has been largely modified by Gustavo Seabra
! for inclusion in the Amber package. (2005)
!
! *********************************************************************
!
! PROGRAM CHARACTERISTICS
! -----------------------
!
! eglcao calculates energy and gradient for dylcao as shown by Seifert.
! The determination of the occupation numbers has been changed to be
! also valid for metallic systems.
!


!In parallel all threads enter here.

   use qm2_dftb_module, only: MDIM,LDIM,NDIM,disper, lmax, dacc, mcharge, &
         izp_str, ks_struct, fermi_str, dftb_3rd_order_str
   use qmmm_module, only: qmmm_nml, qmmm_struct, qm2_struct, qm_gb, qmewald, qmmm_mpi,qm2_params
   use ElementOrbitalIndex, only : elementSymbol
   use constants, only : BOHRS_TO_A, AU_TO_KCAL, AU_TO_EV

   implicit none


   ! Parameters passed in:
   ! =====================
   _REAL_ , intent(in )   :: qm_coords(3,qmmm_struct%nquant_nlink)       ! QM atoms coordinates
   _REAL_ , intent(out)   :: total_e         ! Total energy
   _REAL_ , intent(out)   :: elec_eng        ! Electronic energy
   _REAL_ , intent(out)   :: ethird          ! Third order energy contribution
   _REAL_ , intent(out)   :: enuclr_qmqm     ! Repulsive energy
   integer, intent(out)   :: inner_scf_count ! Number of SCC iterations performed in this trial
   integer, intent(out)   :: outer_scf_count ! Total number of SCC iterations performed
   logical, intent(out)   :: scc_converged   ! SCC procedure converged?
   _REAL_ , intent(inout) :: qmat(*)         ! Electron population per atom
   _REAL_ , intent(out)   :: scf_mchg(qmmm_struct%nquant_nlink) ! Mulliken charges per atom

   ! Locals
   ! ======
   integer :: indj1,indk1,lumo
   integer :: cai, liend, ljend
   integer :: j, izpj, k, li, lj, i
   integer :: n, m
   integer :: nstart, nend, mstart, mend
   integer :: indkn, indjm, indi, indj,indk, indili, indjlj
   _REAL_  :: shifti, shiftj

   !Ewald stuff
   _REAL_  :: ew_corr ! Ewald correction to energy in ev. Needed only if qm_ewald > 0.
   _REAL_  :: temp_pot
   _REAL_  :: gb_escf_corr !GB correction for escf.
   _REAL_  :: elec_eng_old, qtot, eext, egb
   _REAL_  :: ecoul, efermi, spro, dipabs
   integer :: IA, IB, i1, i2

   ! SCC Charge convergency
   _REAL_  :: ediff
   _REAL_  :: chdiff, chdiff_tmp
   integer :: atdiff

   ! Error code from diagonalizer
   integer :: ier

   ! Third order
   _REAL_  :: gaussian
   !  correct dimensions of Uhub and DUhub are the number of atoms
   _REAL_, dimension(100) :: Uhub
   _REAL_, dimension(100) :: DUhub
   _REAL_ :: third_order_h_contrib


#ifdef MPI
   include 'mpif.h'
#endif

   if (qmmm_mpi%commqmmm_master) then

     ! Setup of charge-independent part of H and S
     ! -------------------------------------------
     !
     ! Those matrices depend on the geometry, but do 
     ! not change during the SCF
     do j = 1,qmmm_struct%nquant_nlink
        do k = 1,j   !--> Calculates only the lower triangle

           ! Gets the hamiltonian and overlap terms
           ! referring to this specific pair (j,k)

           call slkmatrices(j,k,qm_coords,ks_struct%au,ks_struct%bu,LDIM)

           ! Puts the calculated block inside the big
           ! hamiltonian and ovelap matrices
           nstart = 1
           nend = ks_struct%ind(k+1)-ks_struct%ind(k)

           do n = nstart, nend

              indkn = ks_struct%ind(k) + n
              mstart = 1
              mend = ks_struct%ind(j+1)-ks_struct%ind(j)

              do m = mstart, mend

                 indjm = ks_struct%ind(j) + m
                 ks_struct%hamil(indjm,indkn) = ks_struct%au(m,n)
                 ks_struct%overl(indjm,indkn) = ks_struct%bu(m,n)

                 ! The actual matrices are symmetric
                 ks_struct%hamil(indkn,indjm) = ks_struct%au(m,n)
                 ks_struct%overl(indkn,indjm) = ks_struct%bu(m,n)
              end do
           end do
        end do
     end do

     ! EXTERNAL FIELD OF MM POINT CHARGES
     ! ----------------------------------
     ! The effect of the MM atoms is taken as an extra term to 
     ! the Hamiltonian, i.e., as an external energy shift to H.
     !
     ! This is calculated only once, outside the SCC process
     ! because the charges don't move.
     !
     if ( qmmm_nml%qmmm_int > 0 .and. (qmmm_nml%qmmm_int /= 5) ) then
        call externalshift(qm_coords,izp_str%izp,ks_struct%shiftE)
     else
        ks_struct%shiftE(1:qmmm_struct%nquant_nlink)=0.0d0
     endif

     elec_eng_old = 0.0d0
     ks_struct%shift  = 0.0d0
     ks_struct%qmold(1:qmmm_struct%nquant_nlink) = qmat(1:qmmm_struct%nquant_nlink)

   end if !(qmmm_mpi%commqmmm_master)

   call timer_start(TIME_QMMMENERGYSCF)

   scc_converged = .false.

!! =============================================================================
!!                           SCF Loop Starts here.
!! =============================================================================
   if ((qmmm_nml%verbosity > 2) .and. qmmm_mpi%commqmmm_master) then
      write (6,*)
      write (6,'(" QMMM SCC-DFTB:  SCC convergence criteria: ")')
      write (6,'(" QMMM SCC-DFTB:      Energy:", 1P,e10.1)') qmmm_nml%scfconv
      write (6,'(" QMMM SCC-DFTB:      Charge:", 1P,e10.1)') qmmm_nml%density_conv
      write (6,'(" QMMM SCC-DFTB:      Telec :", f10.3,"K")') fermi_str%telec
      write (6,'(" QMMM SCC-DFTB:  It#",20X,"Energy",19X,"Diff",14X,"MaxChDiff",3X,"At#",3X,"Index",3X,"Symb",3X,"MullikCh")')
   end if

   !For the moment the master thread does the full scf loop while the non-master
   !threads do a reduced version of it that basically just calls the qmewald code
   !in parallel and then reduces the qmpot array to the master.
   if (qmmm_mpi%commqmmm_master) then
      
      ! In its original format, DFTB uses a general 
      ! eigenvalue diagonalizer to solve HC=SCE.

      ! We here use Canonical Orthogonalization 
      ! (See Szabo section 3.4.5, eq. 3.169) 
      ! to convert it to the more simple H'C'=C'E

      ! BEFORE ENTERING THE SCF LOOP:
      ! 1. Diagonalize the overlap matrix to obtain
      !    the transformation matrix X (ks_struct%xtrans).
      !    X(i,j)=U(i,j).s(j)^-1/2
      
      ! IN EACH SCF ITERATION:
      ! 2. Build the Hamiltonian matrix H (ks_struct%a).
      ! 3. Calculate the transformed hamiltonian H'=Xt.H.X
      ! 4. Diagonalize H' to obtain C' and e
      ! 5. Obtain the eigenvectors of H by: C=X.C'

      ! Copy the overlap matrix into ks_struct%b
      
      ks_struct%b(1:qm2_struct%norbs,1:qm2_struct%norbs) = 0.0d0
      do j=1,qm2_struct%norbs
         ks_struct%b(1:qm2_struct%norbs,j) = ks_struct%overl(1:qm2_struct%norbs,j)
      end do
      
      ! Diagonalize the overlap matrix

      call dftb_matrix_diag(MDIM,ndim,ks_struct%b,ks_struct%ev,ier)

      ! Return:
      ! b  --> 'u' (eigenvectors of S)
      ! ev --> Eigenvalues of 'S' (s)

      ! Invert and sqrt the eigenvalues
      ks_struct%scr1(1:qm2_struct%norbs,1:qm2_struct%norbs)=0.0d0
      call vdinvsqrt( qm2_struct%norbs, ks_struct%ev, ks_struct%ev )
      do i=1,qm2_struct%norbs
         ks_struct%scr1(i,i)=ks_struct%ev(i)
      end do
       
      ! Canonical orthogonalization 
      ! X(i,j)=U(i,j).s(j)^-1/2
      ks_struct%xtrans(1:qm2_struct%norbs,1:qm2_struct%norbs) = 0.0d0
      do j=1,qm2_struct%norbs
         ks_struct%xtrans(1:qm2_struct%norbs,j) &
            = ks_struct%b(1:qm2_struct%norbs,j) * ks_struct%ev(j)
      end do

      



      scc_loop: do inner_scf_count = 1,qmmm_nml%dftb_maxiter

        outer_scf_count = outer_scf_count + 1

        ! Add charge dependent terms (Hubbard, etc. )

        ! Start by defining shift. 
        ! The 'shift' is the difference between the
        ! SCC and non-SCC energies, which is the H^1_{\mu\nu}
        ! in the secular equations.

        ! HAMILTONIAN SHIFT
        ! -----------------
        ! Calculate hamiltonian shift due to SCC - also calculate scf_mchg.

        ! zero the whole ks_struct%shift vector (qmmm_struct%nquant_nlink long)
        ks_struct%shift(1:qmmm_struct%nquant_nlink) = 0.0d0

        call HAMILSHIFT(qm_coords, izp_str%izp,&
              mcharge%uhubb,inner_scf_count,ks_struct%gammamat, &
              ks_struct%shift, qm2_struct%scf_mchg)

        ! EXTERNAL CHARGES SHIFT
        ! ----------------------
        ! Add external charges (QM/MM)
        ! This is added as another energy shift in the Hamiltonian
        ! (Notice the different sign)

         ks_struct%shift(1:qmmm_struct%nquant_nlink) =    &
              ks_struct%shift(1:qmmm_struct%nquant_nlink) &
            - ks_struct%shiftE(1:qmmm_struct%nquant_nlink)

        ! QM EWALD SHIFT OR GB SHIFT
        ! --------------------------
        if ( qmmm_nml%qm_ewald>0 .or. qmmm_nml%qmgb == 2) then
#ifdef MPI
          !At present only the master has up to date scf_mchg. Broadcast it
          !to all threads
          call mpi_bcast(scf_mchg, qmmm_struct%nquant_nlink, &
                         MPI_DOUBLE_PRECISION, 0, qmmm_mpi%commqmmm, ier)
#endif

          if (qmmm_nml%qm_ewald>0) then
!Semi-Parallel
#ifndef SQM
            call qm2_dftb_ewald_shift(scf_mchg)
#endif
          else
            !Must be qmgb==2
!Semi-Parallel
            call qm2_dftb_gb_shift(scf_mchg)
          end if
        end if
  
        ! FINAL HAMILTONIAN
        ! -----------------
        ! Update hamiltonian matrix (here in "a"),
        ! with the charge-dependent part
        
        ! Restore the Hamiltonian matrix in ks_struct%a
        ks_struct%a(1:qm2_struct%norbs,1:qm2_struct%norbs) = 0.0d0
        
        do i = 1,qmmm_struct%nquant_nlink
           indi = ks_struct%ind(i)
           liend = lmax( izp_str%izp(i) )**2
           shifti = ks_struct%shift(i)
           do li = 1,liend
              indili = indi + li
               do j = 1,i          
               ! --> Note: Only the lower triangle is used
                 indj = ks_struct%ind(j)
                 ljend = lmax( izp_str%izp(j) )**2
                 shiftj = ks_struct%shift(j)
                 do lj = 1,ljend
                    indjlj = indj + lj
                    ks_struct%a(indili,indjlj) = ks_struct%hamil(indili,indjlj) &
                          + 0.5*ks_struct%overl(indili,indjlj)*(shifti+shiftj)
                 end do
              end do
           end do
        end do

      ! ==========================
      !   Third order SCC-DFTB
      ! ==========================
      if (dftb_3rd_order_str%do_3rd_order) then

         !Prints debug information 
         if ( dftb_3rd_order_str%debug_print ) then
            write(6,*) 
            if (qmmm_nml%verbosity > 3) then
               write(6,'(A,F8.3)')"DEBUG ==> 3RD ORDER SCC-DFTB: Gaussian_D0 = ", dftb_3rd_order_str%Gaussian_D0
               write(6,'(A,F8.3)')"DEBUG ==> 3RD ORDER SCC-DFTB: Gaussian_g0 = ", dftb_3rd_order_str%Gaussian_G0
               write(6,'(A,F8.3)')"DEBUG ==> 3RD ORDER SCC-DFTB: Gaussian_q0 = ", dftb_3rd_order_str%Gaussian_q0
               write(6,'(A)')"DEBUG ==> 3RD ORDER SCC-DFTB: Hubbard Derivatives:"
               do j=1,qmmm_struct%nquant_nlink
                  write(6,'(A,I5,2X,"(",I2,")",A3,2X,F8.3)')"DEBUG ==> 3RD ORDER SCC-DFTB:  ",&
                        j, &
                        qmmm_struct%iqm_atomic_numbers(j), &
                        elementSymbol( qmmm_struct%iqm_atomic_numbers(j) ),&
                        dftb_3rd_order_str%Hubbard_deriv( qmmm_struct%iqm_atomic_numbers(j) )
               end do
            end if
            write (6,'(A,2X,A3,2X,A15,2X,A15,2x,A15)') "DEBUG ==> 3RD ORDER SCC-DFTB: ", "(i)", "Uhub(i)","Duhub(i)"
         end if ! debug_print

         ! build the atomic intermediates needed for third order that depend on mulliken charges (scf_mchg) for the current SCF cycle
         do i = 1, qmmm_struct%nquant_nlink
            gaussian = -1.0d0 * dftb_3rd_order_str%Gaussian_D0 &
                     * exp( - dftb_3rd_order_str%Gaussian_G0 * (-scf_mchg(i) - dftb_3rd_order_str%Gaussian_q0 )**2 )
            Uhub(i)  = dftb_3rd_order_str%Hubbard_deriv( qmmm_struct%iqm_atomic_numbers(i) ) + gaussian
            DUhub(i) = -2.0d0*dftb_3rd_order_str%Gaussian_G0 * (-scf_mchg(i) - dftb_3rd_order_str%Gaussian_q0 )*gaussian

            if ( dftb_3rd_order_str%debug_print ) then
               write (6,'(A,2X,I3,2X,E15.8,$)') "DEBUG ==> 3RD ORDER SCC-DFTB: ", i, Uhub(i)
               write (6,'(2X,E15.8)')DUhub(i)
            endif
         end do
         ! add the third order Fock matrix contribution
         do i = 1,qmmm_struct%nquant_nlink
            indi = ks_struct%ind(i)
            liend = lmax( izp_str%izp(i) )**2
            shifti = ks_struct%shift(i)
            do li = 1,liend
               indili = indi + li
               do j = 1,i          
                  ! --> Note: Only the lower triangle is used
                  indj = ks_struct%ind(j)
                  ljend = lmax( izp_str%izp(j) )**2
                  shiftj = ks_struct%shift(j)
                  do lj = 1,ljend
                     indjlj = indj + lj
                     third_order_h_contrib = 0.25d0 * ks_struct%overl(indili,indjlj)    * &
                                                    ( Uhub(i)  * scf_mchg(i)**2         + &
                                                      Uhub(j)  * scf_mchg(j)**2         - &
                                                      DUhub(i) * scf_mchg(i)**3 / 3.0d0 - &
                                                      DUhub(j) * scf_mchg(j)**3 / 3.0d0 )
                     ks_struct%a(indili,indjlj) = ks_struct%a(indili,indjlj) + third_order_h_contrib
                  end do
               end do
            end do
         end do
      end if !if (dftb_3rd_order_str%do_3rd_order) then
      ! End 3rd Order contribution, part 1

        ! Now, build H''=Xt . H' .X
        ! 1. Call dsymm to make the SECOND product first: scr1 = H'.  X
        call dsymm('L','L', qm2_struct%norbs, qm2_struct%norbs,&
                    1.0d0 , ks_struct%a     ,MDIM,&
                            ks_struct%xtrans,MDIM,&
                    0.0d0 , ks_struct%scr1  ,MDIM)
        
        ! 2. Call dgemm for the FIRST product: H'' = Xt . scr1
        call dgemm('t','n', qm2_struct%norbs,qm2_struct%norbs,qm2_struct%norbs,&
                    1.0d0 , ks_struct%xtrans,MDIM,&
                            ks_struct%scr1  ,MDIM,&
                    0.0d0 , ks_struct%a     ,MDIM)
        
        call timer_start(TIME_QMMMENERGYSCFDIAG) 
        
        ! ----------------------------
        ! SOLVE the EIGENVALUE PROBLEM
        ! ----------------------------
        ier = 0
        call dftb_matrix_diag(MDIM,ndim,ks_struct%a,ks_struct%ev,ier)
        
        call timer_stop(TIME_QMMMENERGYSCFDIAG) 
  
        ! Convergence failure? 
        if (ier /= 0) then
           write(6,*)
           write(6,*)" QMMM SCC-DFTB: ***************************************************"
           write(6,*)" QMMM SCC-DFTB: ERROR ON EWEVGE (Eigenvalue solver). "
           write(6,*)" QMMM SCC-DFTB: ewevge: ier =",ier,"inner_scf_count=",inner_scf_count
           write(6,*)" QMMM SCC-DFTB: ***************************************************"
           call dftb_conv_failure("eglcao <qm2_dftb_eglcao.f> : ", &
                             "Convergence failure on EWEVGE (Eigenvalue solver).", "Exiting")
        endif

        ! At this point, ks_struct%a has the eigenvectors of F' (C')
        ! we now need to back-transform them into the 
        ! eigenvectors of F: C = X C'
        ! 1. scr1 = X C'
        call dgemm('n','n', qm2_struct%norbs,qm2_struct%norbs,qm2_struct%norbs,     &
                    1.0d0 , ks_struct%xtrans ,MDIM,&
                            ks_struct%a      ,MDIM,&
                    0.0d0 , ks_struct%scr1   ,MDIM)

        ! 2. Transfer the resulting eigenvectors (C) to ks_struct%a:
        ks_struct%a(1:qm2_struct%norbs,1:qm2_struct%norbs)=0.0d0
        do j=1, qm2_struct%norbs
            ks_struct%a(1:qm2_struct%norbs,j) = ks_struct%scr1(1:qm2_struct%norbs,j)
        end do

 
        ! Calculate occupation (occ) and fermi energy (efermi)
        ! using fermi-distribution
        call FERMI(ndim,izp_str%nel,ks_struct%ev,ks_struct%occ,efermi)

        ! Electronic Energy
        ! -----------------
        ! (sum of occupied eigenvalues)       
        elec_eng = 0.0d0
        qm2_struct%nclosed = 0
        qm2_struct%nopenclosed = 0

        do i = 1,ndim
           if (ks_struct%occ(i) < dacc) exit
           elec_eng = elec_eng + ks_struct%occ(i)*ks_struct%ev(i)
           if ( ks_struct%occ(i) == 2.0d0) &
                         qm2_struct%nclosed = qm2_struct%nclosed + 1
           qm2_struct%nopenclosed = qm2_struct%nopenclosed + 1
        end do

        ! Lowest unoccupied level
        lumo = i

        ! MULLIKEN CHARGES
        ! ----------------
        ! qmat will contain the electron populations per atom
        call MULLIKEN(qmmm_struct%nquant_nlink,NDIM,izp_str%izp, &
                      lmax,dacc,qmat,mcharge%qzero,scf_mchg) 
!!!           ! OUTPUT EIGENVECTORS
!        call outeigenvectors(lumo,qm_coords,qmmm_struct%nquant_nlink)


        ! complete calculation of electronic energy:
        ! charge-dependent energy contribution
        ! warning: this will only lead to the right result if convergence
        ! has been reached
        ecoul = 0.0d0
        eext  = 0.0d0

        ! COULOMBIC ELECTRONIC REPULSION ENERGY
        ! -------------------------------------
        do i = 1, qmmm_struct%nquant_nlink
           ecoul = ecoul + &
                   ks_struct%shift(i) &
                     *(qmat(i)+mcharge%qzero( izp_str%izp(i)))
        end do

        ! EXTERNAL CHARGES (QM/MM) COULOMB ENERGY
        ! ---------------------------------------
        ! This term is a simple coulomb interaction between 
        ! the external charge (in ks_struct%shiftE) and the mulliken charge
        ! on the atom. 
!           do i= qmmm_mpi%nquant_nlink_start, qmmm_mpi%nquant_nlink_end
        do i= 1, qmmm_struct%nquant_nlink
           eext = eext + ks_struct%shiftE(i)*(mcharge%qzero( izp_str%izp(i))-qmat(i))
        end do

        ! Third order energy
        ethird=0.0d0
        if (dftb_3rd_order_str%do_3rd_order) then
           do i = 1,qmmm_struct%nquant_nlink
                ethird=ethird-Uhub(i)*scf_mchg(i)**3/6.0d0 - &
                      Uhub(i)*scf_mchg(i)**2*qmat(i)/2.0d0 + &
                     DUhub(i)*scf_mchg(i)**3*qmat(i)/6.0d0
           end do
           ! adjust shift to include the third order part 
           do i = 1,qmmm_struct%nquant_nlink
              ks_struct%shift(i)=ks_struct%shift(i)+Uhub(i)*scf_mchg(i)**2/2.0d0 - &
                               DUhub(i)*scf_mchg(i)**3/6.0d0
           end do
        end if


        ! Electronic Energy
        ! -----------------
        ! remark: elec_eng contains ks_struct%shiftE aready via ev,
        ! ks_struct%shift also contains -ks_struct%shiftE, i.e. ecoul also
        ! contains contributions from EXT
        elec_eng = elec_eng-0.5d0*ecoul + 0.5d0*eext + ethird

        ! =================
        !  SCC CONVERGENCE
        ! =================

        ! Energy difference
        ediff = elec_eng-elec_eng_old

        ! Maximum charge difference
        chdiff = 0.0d0
        do i = 1, qmmm_struct%nquant_nlink
           chdiff_tmp = abs(ks_struct%qmold(i) - qmat(i))
           if (chdiff_tmp > chdiff) then
              chdiff = chdiff_tmp
              atdiff = i
           end if
        end do

        ! Print convergence progress
        if (qmmm_nml%verbosity > 2) then
           if (inner_scf_count > 1) then
              write (6,'(" QMMM SCC-DFTB: ",i4,3X,3(3X, f20.15),3X,I3,3X,I5,3X,A,3X,F10.5)') &
                    outer_scf_count, elec_eng, ediff, chdiff, atdiff, qmmm_nml%iqmatoms(atdiff),&
                    elementSymbol(qmmm_struct%iqm_atomic_numbers(atdiff)), scf_mchg(atdiff)
           else
              write (6,'(" QMMM SCC-DFTB: ", 1X,i3,6X, f20.15)') outer_scf_count, elec_eng
           end if

        end if

        ! Convergency Test
        if ( (abs(ediff) < qmmm_nml%scfconv) &
             .and. (chdiff < qmmm_nml%density_conv) ) then
           scc_converged = .true.
#ifdef MPI
           !Tell other threads that we have converged.
           call mpi_bcast(scc_converged, 1, MPI_LOGICAL, 0, &
                          qmmm_mpi%commqmmm, ier)
#endif

           exit scc_loop
#ifdef MPI
        else
           !Tell other threads that we have not converged.
           call mpi_bcast(scc_converged, 1, MPI_LOGICAL, 0, &
                          qmmm_mpi%commqmmm, ier)
#endif
        end if

!! ------------------------------
!! Here begins the next iteration
!! ------------------------------

        ! Save the last energy
        elec_eng_old = elec_eng

        ! BROYDEN MIXING
        ! --------------
        call broyden(inner_scf_count,qmmm_struct%nquant_nlink,ks_struct%qmold,qmat)
        qmat(1:qmmm_struct%nquant_nlink) &
            = ks_struct%qmold(1:qmmm_struct%nquant_nlink)

     end do scc_loop ! (end SCC)
   else !This is a slave thread, NOT the master
      scc_slave_loop: do inner_scf_count = 1,qmmm_nml%dftb_maxiter
        outer_scf_count = outer_scf_count + 1
        if ( qmmm_nml%qm_ewald>0 .or. qmmm_nml%qmgb == 2) then
#ifdef MPI
          !At present only the master has up to date scf_mchg. Broadcast it
          !to all threads
          call mpi_bcast(scf_mchg, qmmm_struct%nquant_nlink, &
                         MPI_DOUBLE_PRECISION, 0, qmmm_mpi%commqmmm, ier)
#endif              
        
          if (qmmm_nml%qm_ewald>0) then
!Semi-Parallel
#ifndef SQM
            call qm2_dftb_ewald_shift(scf_mchg)
#endif
          else
            !Must be qmgb==2
!Semi-Parallel
            call qm2_dftb_gb_shift(scf_mchg)
          end if             
        end if
        !Wait to be told by the master if we have converged - implicit barrier in the bcast.
#ifdef MPI
        call mpi_bcast(scc_converged, 1, MPI_LOGICAL, 0, qmmm_mpi%commqmmm, ier)
#endif
        if (scc_converged) exit scc_slave_loop
      end do scc_slave_loop
   end if

   call timer_stop(TIME_QMMMENERGYSCF)

!! =============================================================================
!!                               End of SCF loop
!! =============================================================================

   ! Calculate the Mulliken charges from the final electron population.
   if (qmmm_mpi%commqmmm_master .and. .not. scc_converged) then
     !We need to calculate scf_mchg again because we called broyden on the way out
     !of scc_loop and this changed qmat.
     do j = 1, qmmm_struct%nquant_nlink
        scf_mchg(j) = mcharge%qzero( izp_str%izp(j) ) - qmat(j)
     end do
   endif 
#ifdef MPI
   !At present only the master has up to date scf_mchg. Broadcast it
   !to all threads
   if (.not. scc_converged .or. qmmm_nml%qm_ewald==0) &
        call mpi_bcast(scf_mchg, qmmm_struct%nquant_nlink, &
                       MPI_DOUBLE_PRECISION, 0, qmmm_mpi%commqmmm, ier)
#endif

   if (qmmm_mpi%commqmmm_master) then

     ! REPULSIVE ENERGY
     ! ----------------
     ! Repulsive Energy is added after the SCC process is complete.
     !
     call dftb_repulsive(qmmm_struct%nquant_nlink,izp_str%izp,qm_coords,enuclr_qmqm)
     
     
     ! EWALD ENERGY CORRECTION
     ! -----------------------
     ! This correction is necessary because, up to this point,
     ! the energy for the Ewald sum only included half of the term from the MM atoms
     ! and half from the QM atoms. The QM atoms should contribute only half but the
     ! MM atoms should contribute full. The corrects qm_ewald_correct_ee routine 
     ! supposedly corrects for this.
     !
     ! That routine uses the density matrix, and expects it to be indexed according
     ! to the contents of qm2_params%orb_loc. 

#ifndef SQM
     if ( qmmm_nml%qm_ewald>0 ) then

        ! Calculate the correction in Hartree/Bohrs

        call qm2_dftb_ewald_corr(qmmm_struct%nquant_nlink, ew_corr, &
                                 qmewald%mmpot, scf_mchg)!qmat)
        ! Put correction into elec_eng.
        elec_eng = elec_eng + ew_corr

     end if

     if ( qmmm_nml%qmgb == 2) then
        !Calculate the correction for GB - removes double counting of
        !GB energy in both escf and egb.
        gb_escf_corr = 0.0d0
        do i = 1,qmmm_struct%nquant_nlink
           gb_escf_corr = gb_escf_corr + (qm_gb%gb_mmpot(i)+qm_gb%gb_qmpot(i))*scf_mchg(i)
        end do
        elec_eng = elec_eng + (0.5d0*gb_escf_corr*BOHRS_TO_A)
    end if
#endif

     ! total energy
     total_e = elec_eng + enuclr_qmqm
     !Calculate Dispersion Energy
     if (qmmm_nml%dftb_disper == 1) then
        call timer_start(TIME_QMMMDFTBDISPE)
        call  dispersion_energy(qmmm_struct%nquant_nlink,qm_coords)
        total_e = total_e + disper%Edis
        call timer_stop(TIME_QMMMDFTBDISPE)
     endif
   end if !(qmmm_mpi%commqmmm_master)
end subroutine eglcao


! Routine to give a message if convergence fails. It has the same syntax as
! sander_bomb.
subroutine dftb_conv_failure(string1,string2,string3)

   use qmmm_module, only:  qmmm_nml, qmmm_struct, qm2_struct, qmmm_mpi
   use ElementOrbitalIndex, only : elementSymbol
   implicit none
   integer :: i, j

!! Passed in:
   ! The only arguments are the strings which are passed to
   ! sander_bomb in case the execution must be stopped.
   character(len=*) :: string1 ! Usually expected to be the name of the subroutine / file
   character(len=*) :: string2 ! Usually expected to be a message
   character(len=*) :: string3 ! Usually expected to be another message.

!! Locals
   logical :: continue_job
   
!! --
   continue_job = .true.
!    continue_job = .false.
   if (qmmm_mpi%commqmmm_master) then
      write(6,'(/," QMMM SCC-DFTB: !!!! ============= WARNING ============= !!!!")')
      write(6,'(  " QMMM SCC-DFTB: Convergence could not be achieved in this step.")')
      if (continue_job) then
         write(6,'(  " QMMM SCC-DFTB: The calculation will continue, but energies and ")')
         write(6,'(  " QMMM SCC-DFTB: forces for this step will not be accurate. ")')
      else
         write(6,'(  " QMMM SCC-DFTB: The calculation will stop.")')
         write(6,'(/," QMMM SCC-DFTB: Last QM Region Cartesian Coordinates ")')
         write(6,'(" QMMM SCC-DFTB: ","SYM",10X,"X",16X,"Y",16X,"Z",13X,"Charge")')
         do i = 1, qmmm_struct%nquant_nlink
            write(6,'(" QMMM Mullik: ",A2,2X,4F16.10)')  &
                  elementSymbol(qmmm_struct%iqm_atomic_numbers(i)), &
                  (qmmm_struct%qm_coords(j,i), j=1,3), qm2_struct%scf_mchg(i)
         end do
         
         if (qmmm_struct%qm_mm_pairs > 0) then
            write(6,*) ' QMMM SCC-DFTB: number of external charges', qmmm_struct%qm_mm_pairs
            write(6,*) ' QMMM SCC-DFTB: Coordinates of external charges (XYZ)'
            write(6,*)
            write(6,'(" QMMM SCC-DFTB: ",i3)') qmmm_struct%qm_mm_pairs
            do i=1,qmmm_struct%qm_mm_pairs
               write(6,'(" QMMM SCC-DFTB: ",4(2x,f10.6))') (qmmm_struct%qm_xcrd(j,i),j=1,3)
            end do
         end if
         write(6,*) "***************************************************"
         call sander_bomb(string1,string2,string3)
      end if ! if (continue_job) then ...
   end if

   
   
end subroutine dftb_conv_failure


!!!=========================================
!!! Currently, this subroutine is not used.
!!! But I kept it here for the future.
!!!=========================================
!subroutine outeigenvectors(lumo,qm_coords,nn)
!!
!   use qm2_dftb_module, only: MDIM,LDIM,NDIM,disper, lmax, dacc, mcharge, &
!         izp_str, ks_struct, fermi_str
!   use qmmm_module, only : qmmm_nml, qmmm_mpi, qmmm_struct
!!
!   implicit none
!   _REAL_  :: qm_coords(3,nn)
!   integer :: i,j,k,l,nn,lumo,norbs
!   character :: filename*10,tmpstr*3
!
!!   Do k=1,mdim
!!    If (k.lt.10) write(tmpstr,'(I1)')k
!!    If (k.lt.100) write(tmpstr,'(I2)')k
!!    If (k.lt.1000) write(tmpstr,'(I3)')k
!!    filename='CUBE.DAT'//tmpstr
!   open (1,file='CUBE.DAT',status='unknown')
!   rewind 1
!   write(1,*)nn
!   Do l=1,nn
!       write(1,'(I3,3f15.10)') &
!       qmmm_struct%iqm_atomic_numbers(l), &
!       qm_coords(1,l),qm_coords(2,l),qm_coords(3,l)
!   End do
!   norbs=ks_struct%ind(nn+1)
!   write(1,*)norbs
!      do i=1,norbs
!         write(1,'(f15.10)') ks_struct%a(i,6)
!      end do
!   write(1,*)
!   endfile 1
!   close (1)

!end subroutine outeigenvectors



