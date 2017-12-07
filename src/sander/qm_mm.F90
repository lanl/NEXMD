#include "copyright.h"
#include "dprec.fh"
#include "def_time.h"
!-----------------------------------
!Principal code for calculating
!QM potential for QMMM simulations.
!Principal Authors of current code:
!           Ross Walker
!           Mike Crowley
!
! Please send all comments or
! queries to: ross@rosswalker.co.uk
!-----------------------------------

subroutine qm_mm(coords,natom,scaled_mm_charges,f,escf, &
                 periodic,born_radii,one_born_radii, &
                 intdiel,extdiel,Arad,mmcut2, scf_mchg, &
                 ntype, atom_type, atomic_mass, atom_type_index, nstep) 
!
!     Argument list variables:
!
!     coords(natom*3)                 - Cartesian coordinates for all atoms.
!                                       Amber array
!     natom                           - Total number of REAL atoms.
!     qmmm_struct%nquant              - Number of REAL quantum atoms as specified in mdin.
!     iqmatoms(qmmm_struct%nquant)
!                                     - Atom numbers for quantum atoms link atoms given values of -1
!     scaled_mm_charges(natom)          - Atomic charges for mm atoms. (Scaled to elec units)
!     qmmm_struct%iqm_atomic_numbers(qmmm_struct%nquant) - Atomic numbers for qm atoms.
!     qmmm_struct%nlink               - Number of link atoms.
!     f((natom)*3)                    - Atomic forces.
!     qmmm_nml%qmcut               - cutoff in angstroms.
!     qmmm_nml%qmcut2              - cutoff^2 in angstroms^2.
!     escf                            - Heat of formation from QM.
!     qm_mm_pair_list(*)              - Atom-based list for quantum atoms:
!                                       i1,i2,i3. Each QM atom shares same list.
!     qmmm_struct%qm_coords(3,qmmm_struct%nquant+qmmm_struct%nlink)  
!                                     - Cartesian coordinates of quantum atoms.
!                                       (Extracted from coords by qm_extract_coords)

!     Locally defined arrays:
!     dxyzqm(3,qmmm_struct%nquant+qmmm_struct%nlink)     
!                                - Quantum mechanical derivatives from qm-mm
!                                                        interactions.
!     dxyzcl(3,natom)          - Classical derivatives from qm-mm interaction.
!     qm_xcrd(4,natom)         - imaged coords array from amber's array - note is 
!                                ordered in the same order as the pair_list. This
!                                is a DIFFERENT ORDER to the coords array.
!                                1->3 = coordinates, 4 = scaled_mm_charge
!    born_radii(1->natom)      - Effective GB radii - only used when doing qm with gb (and qm_gb==2)
!                                Calculated via an initial call to egb.
!    one_born_radii(1->natom)  - 1.0d0/born_radii(i)
!    scf_mchg                  - nquant long, gets filled with the mulliken charges during scf.
!    mmcut2 - cut^2 in angstroms^2 - cut from cntrl namelist

   use qmmm_module, only : qmmm_nml,qmmm_struct, qm2_struct, qm2_rij_eqns, qmewald,  &
                           qm_gb, qmmm_mpi, qmmm_scratch, qmmm_opnq, qmmm_vsolv, get_atomic_number
   use constants, only : EV_TO_KCAL, KCAL_TO_EV, zero, one, alpb_alpha
   
   use parms, only : cn1, cn2, nttyp 
   use nblist,only: alpha,beta,gamma
   use qm2_extern_module, only: qm2_extern_get_qm_forces
  
   implicit none

#include "assert.fh"

#ifdef MPI
   include 'mpif.h'
#endif

   ! Passed in
   integer, intent(in) :: natom,periodic,nstep, ntype
   _REAL_ , intent(inout)  :: coords(natom*3) !Amber array - adjusted for link atoms
   _REAL_ , intent(in)  :: scaled_mm_charges(natom)
   _REAL_ , intent(out) :: f(natom*3)
   _REAL_ , intent(out) :: escf
   _REAL_ , intent(in) :: born_radii(natom), one_born_radii(natom)
   _REAL_ , intent(in) :: intdiel, extdiel, Arad, mmcut2
   _REAL_ , intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)
   character(len=4), intent(in) :: atom_type(natom)
   integer, intent(in) :: atom_type_index(natom)  
   _REAL_,intent(in):: atomic_mass(natom) 

   ! Locals
   _REAL_ :: alpb_beta, temp
   logical::somethingWrong

   integer :: ier=0
   integer i, j, n, m, offset, qm_no

   ! Locals for link atoms
   _REAL_ :: forcemod(3)
   integer :: lnk_no, mm_no
   integer :: nclatoms

!=============================================================================
!                   START OF QMMM SETUP: allocate list memory
!=============================================================================

      call timer_start(TIME_QMMMSETUP)

   ! Increment the counter of how many time qm_mm routine has been called.
   qmmm_struct%num_qmmm_calls = qmmm_struct%num_qmmm_calls + 1

   ! If this is the first call to the routine, do some initial allocation
   ! that has not been done elsewhere.
   if (qmmm_struct%qm_mm_first_call) then

      allocate ( qmmm_struct%qm_coords(3,qmmm_struct%nquant+qmmm_struct%nlink), stat=ier )
                !Stores the REAL and link atom qm coordinates
      REQUIRE(ier == 0)

     ! Do some initial setup for qm_ewald if in use
     if (qmmm_nml%qm_ewald>0) then
        call timer_start(TIME_QMMMEWALDSETUP)
        qmewald%ewald_startup = .true. !Specify that we haven't done any qm_ewald stuff before this point.
                                        !Essentially that this is the first MD step.
        qmewald%natom = natom !QM ewald needs access to natom in some deep QM routines where it would
                               !not normally be available
        call qm_ewald_setup(qmewald%totkq,qmmm_nml%kappa,qmmm_nml%kmaxqx,qmmm_nml%kmaxqy,qmmm_nml%kmaxqz, &
                            qmmm_nml%ksqmaxq,natom,qmmm_struct%nquant,qmmm_struct%nlink) 
                                               !Also allocates kvec memory (totkq reals)
                                               !,KVec table which is 6 lots of (natom,totkq) reals.
                                               !This is also responsible for dividing up kvectors between
                                               !cpus.
        !If diagnostics are on we can print some info here.
        qmewald%kappa=qmmm_nml%kappa
        if (qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 2) then
           write(6,*) ''
           write(6,*) 'QMMM: Ewald - kmaxq(x,y,z) = ', qmmm_nml%kmaxqx, ',', qmmm_nml%kmaxqy, ',', &
                      qmmm_nml%kmaxqz
           write(6,*) 'QMMM: Ewald -      ksqmaxq = ', qmmm_nml%ksqmaxq
           write(6,*) 'QMMM: Ewald - Total number of k vectors = ', qmewald%totkq
           write(6,*) 'QMMM: Ewald - Kappa = ',  qmewald%kappa
           write(6,*) ''
        end if

        !If we are NOT updating the qm image atom charges on every SCF step for qm_ewald
        !then we initially need to zero the scf_mchg array. Since in this situation on the first step
        !we allow the QM image charges to vary with the SCF. It is only on steps 2 -> N that we
        !keep them fixed during the SCF. This avoids the need for an explicit test of first
        !call in the qm_ewald_calc_mm_pot code.
        if (qmmm_nml%qm_ewald==2) scf_mchg(1:qmmm_struct%nquant_nlink) = zero

        call timer_stop(TIME_QMMMEWALDSETUP)
     else
       qmewald%ewald_startup = .false.
     end if
     
     ! Allocation for QM_GB (qmgb==2)
     if (qmmm_nml%qmgb == 2) then
       ! Calculate dielectric factor
       if (qm_gb%alpb_on) then
         alpb_beta=alpb_alpha*(intdiel/extdiel)
         qm_gb%intdieli = one/(intdiel*(one + alpb_beta))
         qm_gb%extdieli = one/(extdiel*(one + alpb_beta))
         qm_gb%one_Arad_beta = alpb_beta/Arad
       else
         qm_gb%intdieli = 1.0d0/intdiel
         qm_gb%extdieli = 1.0d0/extdiel
       end if
       qm_gb%mmcut2 = mmcut2
     end if
   end if ! ---- first call endif ----------
   call timer_stop(TIME_QMMMSETUP)

!=============================================================================
!                   BUILD NONBOND LIST for QM region
!=============================================================================
   if(periodic ==1) then
      !---- check periodic status and run some checks on the system
      call timer_start(TIME_QMMMCOORDSX)
      !Note also builds list and starts / stops relevant timers.
      call qm_fill_qm_xcrd_periodic(coords,natom, &
           qmmm_struct%iqmatoms,scaled_mm_charges,qmmm_scratch%qm_real_scratch)

      call timer_stop(TIME_QMMMLISTBUILD) !Since QMMMCOORDSX was stopped in qm_fill_qm_xcrd_periodic
                                          !and QMMMLISTBUILD was started.
   else
      !---------------Not Periodic ------------------------------------
      call timer_start(TIME_QMMMLISTBUILD)
      !---- cutoff based on distance from whole QM region
      !nb_list also fills qm_xcrd array and extracts qm atoms.
      call qm_fill_qm_xcrd( coords, natom, scaled_mm_charges) 
      call timer_stop(TIME_QMMMLISTBUILD)
   endif

   call timer_start(TIME_QMMMCOORDSX)

   ! If verbosity is on print the number of QMMM pairs
   if (qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 1) then
      write(6,*) 'QMMM: No. QMMM Pairs per QM atom: ', qmmm_struct%qm_mm_pairs
   end if

   !Finally we need to position the link atoms for this step
   !We base this on the imaged coordinates so we don't need
   !to worry about periodic boundaries. Writes the link atoms
   !to the end of the qmmm_struct%qm_coords array.
   if ( qmmm_struct%nlink > 0 ) then
     call position_link_atoms(qmmm_struct, coords)
   end if

   if (qmmm_mpi%commqmmm_master .AND. (qmmm_nml%verbosity>1 .OR. qmmm_struct%qm_mm_first_call)) then
      call print_link_atom_info(qmmm_struct, qmmm_struct%qm_coords,atom_type )
   end if
       
!=============================================================================
!                   START OF REST OF QMMM SETUP
!=============================================================================
   call timer_stop_start(TIME_QMMMCOORDSX,TIME_QMMMSETUP)
   if(qmmm_struct%qm_mm_first_call) then 

     if (qmmm_mpi%commqmmm_master) then
       ! ---------------------------------------
       ! Print the initial QM region coordinates
       ! ---------------------------------------
       call qm_print_coords(nstep,.true.)
        write(6,'(/80("-")/"  3.1 QM CALCULATION INFO",/80("-"))')
     end if

   end if !if (qmmm_struct%qm_mm_first_call)

   !Also if we are doing nearest solvent and the QM region may have changed reprint
   !things.
   if (qmmm_mpi%commqmmm_master      .and. &
       qmmm_nml%vsolv > 0            .and. &
       qmmm_nml%verbosity > 0        .and. &
       .not. qmmm_struct%qm_mm_first_call) then
     if (mod(nstep,qmmm_vsolv%nearest_qm_solvent_fq)==0) then
        call qm_print_coords(nstep,.false.)
     end if
   end if

   !======================
   !  Setup for QM EWALD
   !======================
   !If we are doing QM Ewald then we need to calculate the Kvectors. We do this on
   !every step since the box dimensions may change and these values depend on the
   !box dimensions.
   if (qmmm_nml%qm_ewald>0) then
     call timer_stop_start(TIME_QMMMSETUP,TIME_QMMMEWALDKTABLE)
!Parallel
     call qm_ewald_calc_kvec(qmewald%kvec, qmewald%dkvec, qmewald%dmkv, qmmm_nml%kmaxqx, &
                             qmmm_nml%kmaxqy, qmmm_nml%kmaxqz,qmmm_nml%ksqmaxq)
     !Next we calculate the KTABLE which is an array of exponentials (in complex
     !Sin, cos notation) in the form of 1->6 by 1->NKvectors wide by 1->Natom (all atoms) long.
     !Note, this routine is very memory intensive and very time consuming. Although if
     !we are using Walker and Crowley PME implementation then we only do the QM-QM
     !table which is fast.
     !Note this routine uses the unimaged coordinates from AMBER coords array so here the
     !MMlink atom must have had its coordinates replaced with the link atom.
     if (qmmm_nml%ifqnt) call adj_mm_link_pair_crd(qmmm_struct,coords)
!Parallel
     call qm_ewald_calc_ktable(natom, qmmm_struct%nquant, qmmm_struct%nlink, coords, qmewald%dmkv)
                                 !Will fill the kvector tables with the exponential values
                                 !memory for ktable_x_cos... should already have been allocated.

     if (qmmm_nml%ifqnt) call rst_mm_link_pair_crd(qmmm_struct, coords)

     call timer_stop_start(TIME_QMMMEWALDKTABLE,TIME_QMMMSETUP)

   end if
   !==========================
   !  End Setup for QM EWALD
   !==========================
   
   !====================================
   ! Set up OPNQ stuff
   ! Taisung Lee, Rutgers, 2011
   !====================================

   if (qmmm_struct%qm_mm_first_call .and. qmmm_opnq%useOPNQ) then   
       qmmm_opnq%OPNQCorrection=zero
       qmmm_opnq%vdWCorrection=zero       
       qmmm_opnq%NB_cutoff=qmmm_nml%qmcut
       qmmm_opnq%switch_cutoff1=min(5.0D0, qmmm_nml%qmcut*0.6)
       qmmm_opnq%switch_cutoff2=min(6.5D0, qmmm_nml%qmcut*0.8)

       if (qmmm_mpi%commqmmm_master .AND. qmmm_nml%verbosity > 1) then
          write(6,'("QMMM: OPNQ correction is turned on")')
          write(6,'("QMMM: OPNQ Switching cutoff1=",f18.8)') &
                 qmmm_opnq%switch_cutoff1
          write(6,'("QMMM: OPNQ Switching cutoff2=",f18.8)') &
                 qmmm_opnq%switch_cutoff2
       end if    
          
       !store the mm type info       
       allocate(qmmm_opnq%MM_atomType(natom), stat=ier )
       REQUIRE(ier == 0)      
       qmmm_opnq%MM_atomType=atom_type_index
       
       n=floor(sqrt(nttyp*2.D0))
       REQUIRE(n==ntype)
       if (associated(qmmm_opnq%supported))  nullify(qmmm_opnq%supported)
       allocate(qmmm_opnq%supported(n), stat=ier )
       REQUIRE(ier == 0)

       if (associated(qmmm_opnq%LJ_r))  nullify(qmmm_opnq%LJ_r)                     
       allocate(qmmm_opnq%LJ_r(n), stat=ier )
       REQUIRE(ier == 0)
       if (associated(qmmm_opnq%LJ_epsilon))  nullify(qmmm_opnq%LJ_epsilon)       
       allocate(qmmm_opnq%LJ_epsilon(n), stat=ier )
       REQUIRE(ier == 0)

       !Get the atomic number for each atom
       if (associated(qmmm_opnq%atomic_number))  nullify(qmmm_opnq%atomic_number)       
       allocate(qmmm_opnq%atomic_number(n), stat=ier )
       REQUIRE(ier == 0)   
       do i=1, natom
         qmmm_opnq%atomic_number(atom_type_index(i))=0
         if (atomic_mass(i) >= 0.01 )  then 
            somethingWrong=.false.
            call get_atomic_number(atom_type(i), atomic_mass(i), j, somethingWrong)
            if (.not. somethingWrong) then
               qmmm_opnq%atomic_number(atom_type_index(i))=j
            end if
         else
            qmmm_opnq%atomic_number(atom_type_index(i))=0
         endif
       end do 
       
       ! store the mm LJ parameters
       temp=2.d0**(-5.d0/6.d0)
       do i=1,n
          j=i*(i+1)/2
          if ( (abs(cn1(j)) <1.0D-5) .and. (abs(cn2(j))<1.0D-5) ) then
             qmmm_opnq%supported(i)=.false.
             qmmm_opnq%LJ_r(i)=0.d0
             qmmm_opnq%LJ_epsilon(i)=0.d0
          else   
            qmmm_opnq%supported(i)=.true. 
            qmmm_opnq%LJ_r(i)=temp*(  ( (cn1(j)/cn2(j)) )**(1.d0/6.d0) )
            qmmm_opnq%LJ_epsilon(i)=0.25d0*cn2(j)*cn2(j)/cn1(j)
          end if
       end do
   
   end if !(qmmm_struct%qm_mm_first_call .and. qmmm_opnq%useOPNQ)
   
   !=============================
   ! Setup for PM3MMX_INTERFACE
   !=============================
   if (qmmm_struct%PM3MMX_INTERFACE) then
      ! at this point, we know qm_mm_pairs, so can allocate memory for qm_mm_pair_atom_numbers array, 
      ! which is qm_mm_pairs long, and stores atomic number of MM atoms in the pair list
      ! This uses less memory since qm_mm_pairs <= (natom-nquant+1)
      if (associated(qmmm_struct%qm_mm_pair_atom_numbers)) &
         nullify(qmmm_struct%qm_mm_pair_atom_numbers)
      allocate ( qmmm_struct%qm_mm_pair_atom_numbers(qmmm_struct%qm_mm_pairs), stat=ier )
      REQUIRE(ier == 0)

      do i=1, qmmm_struct%qm_mm_pairs
         ! atom_type(1:natom)
         ! atomic_mass(1:natom)
         ! qm_mm_pair_atom_numbers(1:qm_mm_pairs) 
         ! so we need to convert index i back to index of (1:natom)
         j=qmmm_struct%qm_mm_pair_list(i)
         call get_atomic_number(atom_type(j), atomic_mass(j), &
                                qmmm_struct%qm_mm_pair_atom_numbers(i))
      end do 
   end if 
   !================================
   ! End Setup for PM3MMX_INTERFACE
   !================================
   
   call timer_stop(TIME_QMMMSETUP)
!======================END OF QMMM SETUP ======================================

   !==================
   ! Calculate Forces
   !==================
   if (qmmm_nml%qmtheory%EXTERN) then
      nclatoms = qmmm_struct%qm_mm_pairs
      if (qmmm_nml%qmmm_int == 5) then
         nclatoms = 0
      end if
      call qm2_extern_get_qm_forces(nstep, qmmm_struct%nquant_nlink, qmmm_struct%qm_coords, &
           nclatoms, qmmm_struct%qm_xcrd, &
           escf, qmmm_struct%dxyzqm, qmmm_struct%dxyzcl)
   else
      call get_qm2_forces(qmmm_mpi%commqmmm_master, qm2_struct%calc_mchg_scf, natom, &
           born_radii, one_born_radii, coords, scaled_mm_charges, qmmm_nml, &
           qmmm_struct, qmmm_scratch, scf_mchg, escf)
   end if

   ! If we are doing qm_ewald then we need to calculate the gradients due to the ewald potential here.
   ! Only available analytically, no numerical gradients are available.
   ! With qm_pme the forces are calculated in ew_recip using the mulliken charges
   if (qmmm_nml%qm_ewald>0) then
     call timer_start(TIME_QMMMFQMEWALD)
     ! Parallel
     call qm_ewald_get_forces(qmmm_struct%qm_xcrd, qmmm_struct%qm_coords,&
                              natom, scf_mchg, &
                              qmmm_nml%qmmmrij_incore, &
                              qmmm_struct%dxyzqm, qmmm_struct%dxyzcl, &
                              qmewald%dkvec, scaled_mm_charges)
     call timer_stop(TIME_QMMMFQMEWALD)
   end if

   !=============================
   ! End Calculation of Forces
   !=============================

   ! NOW WE NEED TO PUT THE CALCULATED FORCES INTO THE SANDER FORCE ARRAY
   call timer_start(TIME_QMMMCOLLATEF)
   do i=1,qmmm_struct%nquant
     m = qmmm_struct%iqmatoms(i)
     m = (m-1)*3
     f(m+1) = f(m+1) - qmmm_struct%dxyzqm(1,i)
     f(m+2) = f(m+2) - qmmm_struct%dxyzqm(2,i)
     f(m+3) = f(m+3) - qmmm_struct%dxyzqm(3,i)
   enddo
!Only need to do MM atoms that are in the list. Only need to do this if the QMMM interaction
!is being calculated using the full orbital interaction.
   if (qmmm_nml%qmmm_int >0 .and. (qmmm_nml%qmmm_int /= 5) ) then
     do i = 1,qmmm_struct%qm_mm_pairs
        m = (qmmm_struct%qm_mm_pair_list(i)-1)*3
        f(m+1) = f(m+1) - qmmm_struct%dxyzcl(1,i)
        f(m+2) = f(m+2) - qmmm_struct%dxyzcl(2,i)
        f(m+3) = f(m+3) - qmmm_struct%dxyzcl(3,i)
     end do
   end if


   if (qmmm_nml%qm_ewald>0 .and. .not. qmmm_nml%qm_pme) then
!If we are doing QM ewald then we have a set of forces on all MM atoms that we need to put
!into the main force array. For qm_pme this is done in ew_recip later on in force.
!Parallel division in the ewald code is over kvectors so all threads have some forces on all
!atoms.

     do i = 1, natom
       offset = (3*i) - 2
       f(offset)   = f(offset)   - qmewald%d_ewald_mm(1,i)
       f(offset+1) = f(offset+1) - qmewald%d_ewald_mm(2,i)
       f(offset+2) = f(offset+2) - qmewald%d_ewald_mm(3,i)
     end do
   end if

   !We need to divide the force on the link atom up between the
   !QM and MM atom.
   do i=1,qmmm_struct%nlink
     mm_no = 3*qmmm_struct%link_pairs(1,i)-2  !location of atom in x array
     lnk_no = qmmm_struct%link_pairs(2,i) !Nquant number of QM atom bound to link atom
     qm_no = 3*qmmm_struct%iqmatoms(lnk_no)-2
     !Note this routine uses the flink in the form -flink. 
     call distribute_lnk_f(forcemod,qmmm_struct%dxyzqm(1,qmmm_struct%nquant+i),coords(mm_no), &
                           coords(qm_no),qmmm_nml%lnk_dis)

     !NOTE: forces are reversed in QM calc with respect to amber force array
     !so we subtract forcemod from MM atom and add it to QM atom.
     !MM atom's new force = FMM(x,y,z) - FORCEMOD(x,y,z)
     f(mm_no) = f(mm_no) - forcemod(1)
     f(mm_no+1) = f(mm_no+1) - forcemod(2)
     f(mm_no+2) = f(mm_no+2) - forcemod(3)

     !QM atom's new force = FQM(x,y,z) - Flink(x,y,z) + FORCEMOD(x,y,z)
     !Note QM forces should be subtracted from sander F array to leave total force.
     f(qm_no) = f(qm_no) - qmmm_struct%dxyzqm(1,qmmm_struct%nquant+i) + forcemod(1)
     f(qm_no+1) = f(qm_no+1) - qmmm_struct%dxyzqm(2,qmmm_struct%nquant+i) + forcemod(2)
     f(qm_no+2) = f(qm_no+2) - qmmm_struct%dxyzqm(3,qmmm_struct%nquant+i) + forcemod(3)

   end do

   call timer_stop(TIME_QMMMCOLLATEF)

   qmmm_struct%qm_mm_first_call = .false.
   qmewald%ewald_startup = .false.

end subroutine qm_mm

!======================END OF QM_MM ======================================

!============== NON PERIODIC QM_XCRD, QM_COORDS and LIST =================
!Build non-periodic list - fill qm_xcrd and extract qm_coords.
subroutine qm_fill_qm_xcrd( x, natom, scaled_mm_charges ) 
   use qmmm_module, only : qmmm_nml,qmmm_struct, qmmm_scratch
   implicit none
!#  include "memory.h" 
   !This routine calculates a qm_mm_pair_list which is of the form
   !listtype atom atom atom ... 
   !for qmmm_struct%nquant atoms (I.e each QM atom excluding any link atoms

   !The list of MM atoms that interact with each QM atom is based
   !on the cut off distance. It then fills qm_xcrd.

   !It also extracts the qm coordinates from the amber x array.

   !Each QM atom will use the identical list since to be included an MM atom
   !only needs to be within cut of any QM atom

   !qm_mm_pair_list comes from qmmm_module.

!Passed in
   integer natom
   _REAL_ , intent(in) ,dimension(3,natom) :: x
   _REAL_, intent(in), dimension(natom) :: scaled_mm_charges

   !Local Variables!
   integer j,m,i,n1
   _REAL_ , dimension(6) :: bxbnd
   _REAL_ x_qm, y_qm, z_qm, dx2, xbnd0, xbnd1, ybnd0, ybnd1, zbnd0, zbnd1
   logical include_atom

!     Find the bounding box limits of the QM region, then find all atoms
!     inside the box + cutoff before calculating or testing distances.
   m=qmmm_struct%iqmatoms(1)
   xbnd0=x(1,m)
   xbnd1=x(1,m)
   ybnd0=x(2,m)
   ybnd1=x(2,m)
   zbnd0=x(3,m)
   zbnd1=x(3,m)
   do j=2,qmmm_struct%nquant
      m = qmmm_struct%iqmatoms(j)
      xbnd0=min(xbnd0,x(1,m))
      xbnd1=max(xbnd1,x(1,m))
      ybnd0=min(ybnd0,x(2,m))
      ybnd1=max(ybnd1,x(2,m))
      zbnd0=min(zbnd0,x(3,m))
      zbnd1=max(zbnd1,x(3,m))
   enddo

   bxbnd(1)=xbnd0-qmmm_nml%qmcut
   bxbnd(2)=xbnd1+qmmm_nml%qmcut
   bxbnd(3)=ybnd0-qmmm_nml%qmcut
   bxbnd(4)=ybnd1+qmmm_nml%qmcut
   bxbnd(5)=zbnd0-qmmm_nml%qmcut
   bxbnd(6)=zbnd1+qmmm_nml%qmcut

   include_atom = .false.

   n1 = 0  ! index of qm_mm_pair_list to which we will put the current MM atom
 !---------- FIRST PASS ----------------------------------------------
   qmmm_scratch%qm_int_scratch(1:natom)=0 !Used for a mask
   do m=1,natom 
      !No short circuit evaluation in fortran so having a series of
      !sperate if statements here should be faster.
      if ( x(1,m) <= bxbnd(1) ) cycle
      if ( x(1,m) >= bxbnd(2) ) cycle
      if ( x(2,m) <= bxbnd(3) ) cycle
      if ( x(2,m) >= bxbnd(4) ) cycle
      if ( x(3,m) <= bxbnd(5) ) cycle
      if ( x(3,m) >= bxbnd(6) ) cycle
          
      qmmm_scratch%qm_int_scratch(m)=1

   enddo
   !Set QM atoms in qm_int_scratch to zero so they don't get counted as pairs.
   qmmm_scratch%qm_int_scratch(qmmm_struct%iqmatoms(1:qmmm_struct%nquant)) = 0

 !---------- Second PASS ----------------------------------------------
   do m=1,natom          ! we loop over all atoms - excluding QM atoms
       if ( qmmm_scratch%qm_int_scratch(m) == 1 ) then
          check_cut: do j=1,qmmm_struct%nquant
             i = qmmm_struct%iqmatoms(j) ! the atom number of the first QM atom
                              ! Find all MM atoms that are within CUT of this QM atom
             x_qm = x(1,i)    ! Get coordinates of QM atom
             y_qm = x(2,i)
             z_qm = x(3,i)
                              ! calculate the distance and see if it is within the cut off
             dx2 = ( ( x_qm - x(1,m) ) * ( x_qm - x(1,m) ) &
                   + ( y_qm - x(2,m) ) * ( y_qm - x(2,m) ) &
                   + ( z_qm - x(3,m) ) * ( z_qm - x(3,m) ) &
                   ) 
             if ( dx2 < qmmm_nml%qmcut2 ) then
                !We include the atom.
                !however, if this is a mm link pair
                !atom then we don't include it.
                if (qmmm_struct%mm_link_mask(m)) then
                  !We don't include it since it is a link atom
                  exit check_cut
                end if
                include_atom = .true.
                exit check_cut
             end if
          end do check_cut

          if ( include_atom ) then
                              !include this mm atom in the list
             n1 = n1+1
             qmmm_struct%qm_mm_pair_list( n1 ) = m
             qmmm_struct%qm_xcrd(1:3,n1) = x(1:3,m)
             qmmm_struct%qm_xcrd(4,n1) = scaled_mm_charges(m) 
             include_atom=.false.
          end if
       end if
    end do
    qmmm_struct%qm_mm_pairs = n1

    !Extract QM atoms from x into qmcoords array.
    do m=1,qmmm_struct%nquant
        i=qmmm_struct%iqmatoms(m)
        qmmm_struct%qm_coords(1:3,m) = x(1:3,i)
    enddo

end subroutine qm_fill_qm_xcrd

!=============================================================================
!             QM_FILL_QM_XCRD PERIODIC
!=============================================================================
subroutine qm_fill_qm_xcrd_periodic(x,natom, &
                                    iqmatoms,scaled_mm_chrgs,real_scratch)

   use qmmm_module, only: qmmm_nml,qmmm_struct, qmmm_scratch
   use nblist, only: a,b,c,alpha,beta,gamma,ucell,recip,sphere
   use constants, only : zero, one, half, two

   implicit none
   integer , intent(in) :: natom
   _REAL_ , intent(in), dimension(3,natom) :: x
   integer, intent(in), dimension(qmmm_struct%nquant) :: iqmatoms
   _REAL_, intent(in), dimension(natom) :: scaled_mm_chrgs
   _REAL_, intent(out), dimension(3,natom) :: real_scratch 

   integer j, jqmatom, ier
   _REAL_ :: xbnd0,xbnd1,ybnd0,ybnd1,zbnd0,zbnd1
   _REAL_ :: offset(3), frac(3), xx,yy,zz,xtmp,ytmp,ztmp,one_nquant
   _REAL_ , dimension(6) :: bxbnd

   integer :: m,n, n1
   _REAL_ :: fbndx0, fbndx1, fbndy0, fbndy1, fbndz0, fbndz1, dx2
   logical :: include_atom

!Move the first QM atom to be at the origin.
!   iqm_one = iqmatoms(1)
!   frac(1:3) = x(1,iqm_one)*recip(1,1:3)+x(2,iqm_one)*recip(2,1:3)+ &
!               x(3,iqm_one)*recip(3,1:3)

!Moving the first QM atom can be an issue if the box size is small
!since it could end up with atoms that can stick outside the box.
!Thus we can change this to calculate the center of the coordinates of
!the QM region and move this instead. This could cause problems if the QM
!region diffuses significantly but if that happens we are in a whole
!world of problems anyway and it should trigger the cutoff too large error.
   xtmp=zero; ytmp=zero; ztmp=zero
   do j = 1, qmmm_struct%nquant
     jqmatom = iqmatoms(j)
     xtmp = xtmp + x(1,jqmatom)
     ytmp = ytmp + x(2,jqmatom)
     ztmp = ztmp + x(3,jqmatom)
   end do
   one_nquant = 1.0d0/real(qmmm_struct%nquant)
   xtmp = xtmp * one_nquant
   ytmp = ytmp * one_nquant
   ztmp = ztmp * one_nquant

   frac(1:3) = xtmp*recip(1,1:3)+ytmp*recip(2,1:3)+ &
               ztmp*recip(3,1:3)
   frac(1:3) = frac(1:3) - anint(frac(1:3))
   offset(1) = frac(1)*ucell(1,1) + frac(2)*ucell(1,2) + frac(3)*ucell(1,3)
   offset(2) = frac(1)*ucell(2,1) + frac(2)*ucell(2,2) + frac(3)*ucell(2,3)
   offset(3) = frac(1)*ucell(3,1) + frac(2)*ucell(3,2) + frac(3)*ucell(3,3)

   do j = 1, qmmm_struct%nquant
      jqmatom = iqmatoms(j)
      xx = x(1,jqmatom) - offset(1)
      yy = x(2,jqmatom) - offset(2)
      zz = x(3,jqmatom) - offset(3)
      frac(1:3) = xx*recip(1,1:3) + yy*recip(2,1:3) + zz*recip(3,1:3)
      frac(1:3) = frac(1:3) - anint(frac(1:3))
      real_scratch(1,j)=frac(1)*ucell(1,1) + frac(2)*ucell(1,2) + frac(3)*ucell(1,3)
      real_scratch(2,j)=frac(1)*ucell(2,1) + frac(2)*ucell(2,2) + frac(3)*ucell(2,3)
      real_scratch(3,j)=frac(1)*ucell(3,1) + frac(2)*ucell(3,2) + frac(3)*ucell(3,3)
   end do

   xbnd0=zero; xbnd1=zero; ybnd0=zero; ybnd1=zero; zbnd0=zero; zbnd1=zero
   
   do j=1,qmmm_struct%nquant
      xbnd0=min(xbnd0,real_scratch(1,j))
      xbnd1=max(xbnd1,real_scratch(1,j))
      ybnd0=min(ybnd0,real_scratch(2,j))
      ybnd1=max(ybnd1,real_scratch(2,j))
      zbnd0=min(zbnd0,real_scratch(3,j))
      zbnd1=max(zbnd1,real_scratch(3,j))
   enddo
   xbnd0=xbnd0-qmmm_nml%qmcut
   ybnd0=ybnd0-qmmm_nml%qmcut
   zbnd0=zbnd0-qmmm_nml%qmcut
   xbnd1=xbnd1+qmmm_nml%qmcut
   ybnd1=ybnd1+qmmm_nml%qmcut
   zbnd1=zbnd1+qmmm_nml%qmcut

  !------ Check if QM region plus cutoff around it is too large for this box
  !---       The sphere method is the simplest check, but we can get more
  !---       sophisticated using the distance between parallel faces later if we need to...mfc
  !---       sphere is calculated in ew_box.f and used here.
   if( (xbnd1-xbnd0 > two*sphere) .or. (ybnd1-ybnd0 > two*sphere) .or. (zbnd1-zbnd0 > two*sphere) )then
      write(6,*) " ****************************************************"
      write(6,*) " ERROR: QM region + cutoff larger than box dimension:"
      write(6,'(2X,"QM-MM Cutoff = ",f8.4)') qmmm_nml%qmcut
      write(6,*) "  Coord   Lower     Upper    Size    Radius of largest sphere inside unit cell"
      write(6,'(5X,"X",4(2X,f8.3))') xbnd0, xbnd1, xbnd1-xbnd0, sphere
      write(6,'(5X,"Y",4(2X,f8.3))') ybnd0, ybnd1, ybnd1-ybnd0, sphere
      write(6,'(5X,"Z",4(2X,f8.3))') zbnd0, zbnd1, zbnd1-zbnd0, sphere
      write(6,*) " ****************************************************"
      call sander_bomb("QM_CHECK_PERIODIC<qm_mm.f>", &
        "QM region + cutoff larger than box", &
        "cannot continue, need larger box.")
   endif

   bxbnd(1)=xbnd0+offset(1)
   bxbnd(2)=xbnd1+offset(1)
   bxbnd(3)=ybnd0+offset(2)
   bxbnd(4)=ybnd1+offset(2)
   bxbnd(5)=zbnd0+offset(3)
   bxbnd(6)=zbnd1+offset(3)

!FILL QM_XCRD Periodic
  !---- First move center of QM region to origin: find offset
   offset(1)=(bxbnd(2)+bxbnd(1))*half
   offset(2)=(bxbnd(4)+bxbnd(3))*half
   offset(3)=(bxbnd(6)+bxbnd(5))*half

!  !---- Create new fractional coordinates with new origin
!    !--- find bounds and cut in fracs
!   fbndx0 = (bxbnd(1)-offset(1))*bxinv(1)
!   fbndx1 = (bxbnd(2)-offset(1))*bxinv(1)
!   fbndy0 = (bxbnd(3)-offset(2))*bxinv(2)
!   fbndy1 = (bxbnd(4)-offset(2))*bxinv(2)
!   fbndz0 = (bxbnd(5)-offset(3))*bxinv(3)
!   fbndz1 = (bxbnd(6)-offset(3))*bxinv(3)

    !---- run through coords to get fracs, select those within
    !     bounding box + cutoff

   qmmm_scratch%qm_int_scratch(1:natom) = 0  !Used for a mask
   do m=1,natom
      xx=x(1,m)-offset(1)
      yy=x(2,m)-offset(2)
      zz=x(3,m)-offset(3)
      frac(1:3) = xx*recip(1,1:3) + yy*recip(2,1:3) + zz*recip(3,1:3)
      frac(1:3) = frac(1:3) - anint(frac(1:3))
      xx = frac(1)*ucell(1,1) + frac(2)*ucell(1,2) + frac(3)*ucell(1,3)
      if( (xx > bxbnd(1)-offset(1)) .and. (xx < bxbnd(2)-offset(1)) ) then
         yy=frac(1)*ucell(2,1) + frac(2)*ucell(2,2) + frac(3)*ucell(2,3)
         if( ( yy > bxbnd(3)-offset(2)) .and. (yy < bxbnd(4)-offset(2)) ) then
            zz=frac(1)*ucell(3,1) + frac(2)*ucell(3,2) + frac(3)*ucell(3,3)
            if( (zz > bxbnd(5)-offset(3)) .and. (zz < bxbnd(6) -offset(3)) ) then
               !------- This one inside box ------------------
               qmmm_scratch%qm_int_scratch(m)=1
               real_scratch(1,m) = xx
               real_scratch(2,m) = yy
               real_scratch(3,m) = zz
            endif
         endif
      endif
   enddo

!Fill the qm coordinate array with the imaged QM atoms
   do n=1,qmmm_struct%nquant
      m=iqmatoms(n)
      qmmm_struct%qm_coords(1:3,n) = real_scratch(1:3,m)
      qmmm_scratch%qm_int_scratch(iqmatoms(n)) = 0
   enddo
   call timer_stop_start(TIME_QMMMCOORDSX,TIME_QMMMLISTBUILD)
!Now calculate the pair list and fill qm_xcrd
   include_atom = .false.
   n1 = 0
 !---------- Second PASS ----------------------------------------------
   do m=1,natom          ! we loop over all atoms - excluding QM atoms
      if ( qmmm_scratch%qm_int_scratch(m) > 0 ) then
         check_cut: do j=1,qmmm_struct%nquant
            ! Find all MM atoms that are within CUT of this QM atom
            xx = qmmm_struct%qm_coords(1,j) ! Get coordinates of QM atom
            yy = qmmm_struct%qm_coords(2,j)
            zz = qmmm_struct%qm_coords(3,j)
            ! calculate the distance and see if it is within the cut off
            dx2 = ( ( xx - real_scratch(1,m) ) * ( xx - real_scratch(1,m) ) &
               + ( yy - real_scratch(2,m) ) * ( yy - real_scratch(2,m) ) &
               + ( zz - real_scratch(3,m) ) * ( zz - real_scratch(3,m) ) &
               )
            if ( dx2 < qmmm_nml%qmcut2 ) then
              !We include the atom.
              !however, if this is a mm link pair
              !atom then we don't include it.
              if (qmmm_struct%mm_link_mask(m)) then
                !We don't include it since it is a link atom
                exit check_cut
              end if
              include_atom = .true.
              exit check_cut
            end if
         end do check_cut
         if ( include_atom ) then
            n1 = n1+1
            qmmm_struct%qm_mm_pair_list( n1 ) = m
            qmmm_struct%qm_xcrd(1,n1)=real_scratch(1,m)
            qmmm_struct%qm_xcrd(2,n1)=real_scratch(2,m)
            qmmm_struct%qm_xcrd(3,n1)=real_scratch(3,m)
            qmmm_struct%qm_xcrd(4,n1)=scaled_mm_chrgs(m)
            include_atom=.false.
         end if
      end if                    ! endif of skip_m
   end do
   qmmm_struct%qm_mm_pairs = n1

end subroutine qm_fill_qm_xcrd_periodic

! -----------------------------------------------------------------------
! Calculate energy and forces with built-in semiempirical or DFTB methods
! -----------------------------------------------------------------------
subroutine get_qm2_forces(master, calc_mchg_scf, natom, &
     born_radii, one_born_radii, coords, scaled_mm_charges, qmmm_nml, &
     qmmm_struct, qmmm_scratch, scf_mchg, escf)

  use qmmm_nml_module, only : qmmm_nml_type
  use qmmm_struct_module, only : qmmm_struct_type
  use qmmm_module, only : qmmm_scratch_structure, qmmm_opnq
  use constants, only : zero, EV_TO_KCAL

  implicit none

  logical, intent(in) :: master
  logical, intent(in) :: calc_mchg_scf
  integer, intent(in) :: natom
  _REAL_ , intent(in) :: born_radii(natom), one_born_radii(natom)
  _REAL_ , intent(in) :: coords(natom*3)
  _REAL_ , intent(in) :: scaled_mm_charges(natom)
  type(qmmm_nml_type)         , intent(inout) :: qmmm_nml
  type(qmmm_struct_type)      , intent(inout) :: qmmm_struct
  type(qmmm_scratch_structure), intent(inout) :: qmmm_scratch
  _REAL_                      , intent(inout) :: scf_mchg(qmmm_struct%nquant_nlink)
  _REAL_ , intent(out) :: escf
  
  integer :: i, j
  _REAL_  :: total_energy
  logical, save :: first_call = .true.

  if (first_call) then

     first_call = .false.

     call timer_start(TIME_QMMMSETUP)
     ! Load semiempirical parameters
     ! Also does a lot of memory allocation and pre-calculates all the STO-6G orbital expansions.
     call qm2_load_params_and_allocate(qm2_struct,qmmm_struct)
     call timer_stop(TIME_QMMMSETUP)

     if (master) then

       ! Print a summary about memory usage
       ! WARNING - FOR THE NUMBERS PRODUCED BY THE PRINT DYN MEM ROUTINE TO BE ACCURATE ALL
       ! MEMORY ALLOCATION MUST HAVE BEEN DONE BY THIS STAGE.
        call qm_print_dyn_mem(qm2_struct,qmmm_struct, natom,qmmm_struct%qm_mm_pairs)

       ! Finally print the result header that was skipped in sander.
        write(6,'(/80("-")/"   4.  RESULTS",/80("-")/)')

     end if

   end if

   call timer_start(TIME_QMMMRIJEQNS)
   ! Parallel: Calculate RIJ and many related equations here.
   ! Necessary memory allocation is done inside the routine.
   call qm2_calc_rij_and_eqns(qmmm_struct, qmmm_struct%qm_coords, qmmm_struct%nquant_nlink,qmmm_struct%qm_xcrd, &
                              natom, qmmm_struct%qm_mm_pairs)
                                !and store them in memory to save time later.
   call timer_stop_start(TIME_QMMMRIJEQNS,TIME_QMMMENERGY)

   ! Parallel: Calculate SCF Energy
   call qm2_energy(cosmo_c_struct,qm2_struct,qm2ds, qmmm_struct, escf, scf_mchg, natom, born_radii, one_born_radii, coords, scaled_mm_charges)

   call timer_stop(TIME_QMMMENERGY)

   ! Calculate forces
   call timer_start(TIME_QMMMFQM)

   qmmm_struct%dxyzqm=zero
   if (qmmm_nml%qmtheory%DFTB) then
      ! Partially Parallel
      call qm2_dftb_get_qm_forces(qm2_sturct,qmmm_struct, qmmm_struct%dxyzqm)
   else
      ! Parallel
      call qm2_get_qm_forces(qm2_struct,qmmm_struct, qmmm_struct%dxyzqm)
   end if

   call timer_stop(TIME_QMMMFQM)

   call timer_start(TIME_QMMMFQMMM)

   if ( qmmm_nml%qmmm_int>0 .and. (qmmm_nml%qmmm_int /= 5) ) then
      qmmm_struct%dxyzcl=zero
      if (qmmm_nml%qmtheory%DFTB) then
         ! Parallel
         call qm2_dftb_get_qmmm_forces(qm2_struct,qmmm_struct, qmmm_struct%dxyzcl,qmmm_struct%dxyzqm, qmmm_scratch%qm_real_scratch(1), &
              qmmm_scratch%qm_real_scratch(natom+1), &
              qmmm_scratch%qm_real_scratch(2*natom+1), &
              qmmm_scratch%qm_real_scratch(3*natom+1))
      else
         ! Parallel
         if (qmmm_nml%qmmm_switch) then
            !Calculate Mulliken charges that will be used in QM/MM switching function
            do i=1,qmmm_struct%nquant_nlink
              call qm2_calc_mulliken(qm2_struct,i,scf_mchg(i))
            end do
         end if
         call qm2_get_qmmm_forces(qm2_struct,qmmm_struct, qmmm_struct%dxyzqm,qmmm_struct%qm_xcrd,qmmm_struct%dxyzcl,scf_mchg)
      end if
   end if

   call timer_stop(TIME_QMMMFQMMM)

   ! -----------------------------------------------------------
   ! Calculate Mulliken Charges for Later Printing if Needed
   ! Note: at present we calculate the mulliken charges even
   !       if we don't need them for printing since one might
   !       want to use them for dipole calculations etc. It is
   !       not very expensive to calculate them so might as well
   !       do it on every MD step.
   ! -----------------------------------------------------------
   if (master) then
      if ( qmmm_nml%qmtheory%DFTB .or. calc_mchg_scf &
           .or. qmmm_nml%qm_ewald == 2 .or. qmmm_nml%qmmm_switch ) then
         ! Mulliken charges have already been calculated and stored.
      else
         do i=1,qmmm_struct%nquant_nlink
            !Need to calculate Mulliken charges here.
            call qm2_calc_mulliken(qm2_struct,i,scf_mchg(i))
         end do
      end if
   end if

   ! Print some extra information if verbosity level is > 0
   if (master) then
      call qm2_print_energy(cosmo_c_struct,qm2_struct, qmmm_nml%verbosity, qmmm_nml%qmtheory, escf, qmmm_struct)
      if (qmmm_nml%verbosity > 3) then
      
        if (qmmm_nml%qm_ewald>0) then 
           write (6,'("QMMM: QM - MM Coulombic Interaction = ",f18.8," eV (",f18.8," KCal/mol)")') &
                     qmmm_struct%coulombic_eng, qmmm_struct%coulombic_eng *EV_TO_KCAL  
        end if
        if (qmmm_opnq%useOPNQ) then     
                write (6,'("QMMM: QM - MM OPNQ correction = ",f18.8," eV (",f18.8," KCal/mol)")') &
                     qmmm_opnq%OPNQCorrection, qmmm_opnq%OPNQCorrection*EV_TO_KCAL                 
                write (6,'("QMMM: QM - MM vdW correction = ",f18.8," eV (",f18.8," KCal/mol)")') &
                     qmmm_opnq%vdWCorrection, qmmm_opnq%vdWCorrection*EV_TO_KCAL                 
                write (6,'("QMMM: QM - MM Total OPNQ correction = ",f18.8," eV (",f18.8," KCal/mol)")') &
                     qmmm_opnq%vdWcorrection+qmmm_opnq%OPNQCorrection,  &
                     (qmmm_opnq%vdWCorrection+qmmm_opnq%OPNQCorrection)*EV_TO_KCAL                 

        end if
      
         write (6,'("QMMM:")')
         write (6,'("QMMM: Forces on QM atoms from SCF calculation")')
         write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,qmmm_struct%dxyzqm(1,j), qmmm_struct%dxyzqm(2,j), &
              qmmm_struct%dxyzqm(3,j), j=1,qmmm_struct%nquant_nlink)
         if (qmmm_nml%verbosity > 4) then
            !Also print info in KJ/mol
            write (6,'("QMMM:")')
            write (6,'("QMMM: Forces on QM atoms from SCF calculation (KJ/mol)")')
            write (6,'("QMMM: Atm ",i6,": ",3f20.14)') (j,qmmm_struct%dxyzqm(1,j)*4.184d0, &
                 qmmm_struct%dxyzqm(2,j)*4.184d0, qmmm_struct%dxyzqm(3,j)*4.184d0, &
                 j=1,qmmm_struct%nquant_nlink)
         end if
      end if
   end if

end subroutine get_qm2_forces
