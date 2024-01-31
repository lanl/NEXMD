! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"
module qmmm_vsolv_module
! ----------------------------------------------------------------------
! PURPOSE: Data type holding variable solvent information
!
! Author: Andreas W. Goetz
!         <agoetz@sdsc.edu>
! Date  : October 2010
!         March 2011
!
! Based on previous data type and subroutines contained in qmmm_module
! as written by Ross Walker
!
! Methods:
!   qmmm_vsolv_store_parameters
!      initializes data type and stores parameters
!
!   new
!     allocate data type
!
!   delete
!      deallocate data type
!
!   broadcast
!      broadcast data type
!      and allocate arrays if this has not been done before
!          (for example, slave threads)
!
! ----------------------------------------------------------------------

    implicit none

    private

    public :: qmmm_vsolv_type

    public :: read_vsolv_nml
    public :: qmmm_vsolv_store_parameters
    public :: new, delete
    public :: print
#ifdef MPI
    public :: broadcast
#endif

    type qmmm_vsolv_type

        logical :: debug

        integer :: verbosity

        ! NOTE:
        ! THIS VARIABLE IS USED FOR ADAPTIVE QM/MM
        ! IT HAS NOTHING TO DO WITH THE VARIABLE SOLVENT SCHEME
        ! IT IS SET IN adaptive_qmmmm()[qmmm_adaptive_module] BECAUSE WE NEED TO
        ! ACCESS DATA OUTSIDE OF THE qmmm_adaptive_module
        ! qmmm_adaptive_module ASSUMES THAT THIS VARIABLE IS ACCESSIBLE ON ALL THREADS,
        ! BOUT GROUP MASTERS AND SLAVES FOR MULTISANDER ADAPTIVE QM/MM RUNS
        logical :: recalculate

        ! Number of nearest QM solvent molecules to the QM region
        ! These will be treated as QM, too.
        integer :: nearest_qm_solvent

        ! Frequency with which to check for nearest solvent
        ! 1 = (Default) on every step
        integer :: nearest_qm_solvent_fq

        ! The residue name of the solvent molecules you want considered as QM, e.g. = WAT
        character(len=4) :: nearest_qm_solvent_resname

        ! Choose which atom of the solvent residue should determine the distance to the QM region
        ! 0 = the closest atom
        integer :: nearest_qm_solvent_center_id

        ! Choose which atom in the fixed QM region should be the center of QM region
        ! 0 = center of mass of fixed QM region
        integer :: qm_center_atom_id

        ! Number of fixed quantum atoms (excluding link atoms)
        integer :: fixed_nquant

        ! Total number of solvent molecules in simulation
        integer :: nsolv_res

        ! Number of atoms per solvent residue
        integer :: natom_solv_res

        ! Original sorted iqmatoms array for the fixed QM region, orig_nquant long
        integer, dimension(:), pointer :: fixed_iqmatoms => null()

        ! First atom of solvent residue 1 to nsolv_res
        integer, dimension(:), pointer :: solvent_pointers => null()

        ! First atom of solvent residue 1 to nearest_qm_solvent
        ! I.e. the first atom number of the nearest residues to the QM region
        ! that should be treated as QM.
        integer, dimension(:), pointer :: nearest_solvent_pointers => null()

        ! Distances of QM solvents from the QM center
        ! These values are calculated in qmmm_vsolv_identify_nearest_solvent
        ! in order to determine the solvent residues which needs to be included in
        ! QM region
        ! These distances are going to be used in qmmm_adaptive_module for
        ! calculating the weights of each partition
        _REAL_, dimension(:), pointer :: nearest_solvent_distances => null()

        ! Number of bond types as read from the prmtop before any QM modifications are made - needed
        ! to be able to call setbon multiple times and not have new QM-H atom types keep being created.
        integer :: prmtop_numbnd

        ! Bond, angle and dihedral types with bonds, angles, dihedrals
        ! involving hydrogen and other atoms
        ! Need to be stored since setbon, setang and setdih will simply
        ! delete list parameters that are in the quantum region and we will need
        ! to rebuild the list after each re-determination of the QM region

        ! bonds with hydrogens
        integer :: nbonh
        integer, dimension(:), pointer :: iibh => null()
        integer, dimension(:), pointer :: ijbh => null()
        integer, dimension(:), pointer :: icbh => null()

        ! bonds without hydrogens
        integer :: nbona
        integer, dimension(:), pointer :: iiba => null()
        integer, dimension(:), pointer :: ijba => null()
        integer, dimension(:), pointer :: icba => null()

        ! angles with hydrogens
        integer :: ntheth
        integer, dimension(:), pointer :: iith => null() ! i24
        integer, dimension(:), pointer :: ijth => null() ! i26
        integer, dimension(:), pointer :: ikth => null() ! i28
        integer, dimension(:), pointer :: icth => null() ! i30

        ! angles without hydrogens
        integer :: ntheta
        integer, dimension(:), pointer :: iita => null() ! i32
        integer, dimension(:), pointer :: ijta => null() ! i34
        integer, dimension(:), pointer :: ikta => null() ! i36
        integer, dimension(:), pointer :: icta => null() ! i38

        ! dihedrals with hydrogens
        integer :: nphih
        integer, dimension(:), pointer :: iiph => null() ! i40
        integer, dimension(:), pointer :: ijph => null() ! i42
        integer, dimension(:), pointer :: ikph => null() ! i44
        integer, dimension(:), pointer :: ilph => null() ! i46
        integer, dimension(:), pointer :: icph => null() ! i48

        ! dihedrals without hydrogens
        integer :: nphia
        integer, dimension(:), pointer :: iipa => null() ! i50
        integer, dimension(:), pointer :: ijpa => null() ! i52
        integer, dimension(:), pointer :: ikpa => null() ! i54
        integer, dimension(:), pointer :: ilpa => null() ! i56
        integer, dimension(:), pointer :: icpa => null() ! i58

    end type qmmm_vsolv_type

    interface delete
        module procedure delete_qmmm_vsolv_type
    end interface delete

    interface new
        module procedure new_qmmm_vsolv_type
    end interface new

    interface print
        module procedure print_qmmm_vsolv_type
    end interface print

#ifdef MPI
    interface broadcast
        module procedure broadcast_qmmm_vsolv_type
    end interface broadcast
#endif

contains

    subroutine read_vsolv_nml(self, nres)

        use file_io_dat

        implicit none
        type(qmmm_vsolv_type), intent(out) :: self
        integer, intent(in) :: nres

        integer :: ierr

        integer :: nearest_qm_solvent
        integer :: nearest_qm_solvent_fq
        integer :: nearest_qm_solvent_center_id
        integer :: qm_center_atom_id
        character(len=4) :: nearest_qm_solvent_resname
        integer :: verbosity
        integer :: debug
        namelist /vsolv/ nearest_qm_solvent, nearest_qm_solvent_fq, &
            nearest_qm_solvent_center_id, qm_center_atom_id, &
            nearest_qm_solvent_resname, verbosity, debug

        call flush (6)
        ! Default values
        debug = 0
        verbosity = 0
        nearest_qm_solvent = 0    ! Do not change the QM region to keep the nearest solvents.
        nearest_qm_solvent_fq = 1 ! Do update of QM solvent on every step if nearest_qm_solvent > 0
        nearest_qm_solvent_resname = 'WAT ' !by default assume solvent is resname WAT
        nearest_qm_solvent_center_id = 0 ! By default use the closest atom of the residue
        qm_center_atom_id = 0 ! By default use the center of mass of fixed QM region

        ! Read namelist
        rewind 5
        read (5, nml=vsolv, iostat=ierr)

        if (ierr > 0) then
            call sander_bomb('read_vsolv_nml (qmmm_vsolv_module)', &
                '&vsolv namelist read error', &
                'Please check your input.')
        end if

        self%debug = .false.
        if (debug > 0) then
            self%debug = .true.
        end if

        call int_legal_range('QMMM: (QM-MM nearest_qm_solvent) ', nearest_qm_solvent, 0, nres)
        call int_legal_range('QMMM: (QM-MM nearest_qm_solvent_fq) ', nearest_qm_solvent_fq, 1, 99999999)
        self%nearest_qm_solvent = nearest_qm_solvent
        self%nearest_qm_solvent_fq = nearest_qm_solvent_fq
        self%nearest_qm_solvent_resname = nearest_qm_solvent_resname
        self%nearest_qm_solvent_center_id = nearest_qm_solvent_center_id
        self%qm_center_atom_id = qm_center_atom_id
        self%verbosity = verbosity

        ! call print(self)

    end subroutine read_vsolv_nml

    ! The pristine MM bond/angle/dihedral status needs to be stored
    ! to be able to reset the bond/angle/dihedral list for changing
    ! QM/MM regions
    subroutine qmmm_vsolv_store_parameters(self, numbnd, &
        iibh, ijbh, icbh, &
        iiba, ijba, icba, &
        iith, ijth, ikth, icth, &
        iita, ijta, ikta, icta, &
        iiph, ijph, ikph, ilph, icph, &
        iipa, ijpa, ikpa, ilpa, icpa)

        implicit none

        type(qmmm_vsolv_type), intent(inout) :: self
        integer, intent(in) :: numbnd
        integer, intent(in) :: iibh(self%nbonh), ijbh(self%nbonh), icbh(self%nbonh)
        integer, intent(in) :: iiba(self%nbona), ijba(self%nbona), icba(self%nbona)
        integer, intent(in) :: iith(self%ntheth), ijth(self%ntheth), ikth(self%ntheth), icth(self%ntheth)
        integer, intent(in) :: iita(self%ntheta), ijta(self%ntheta), ikta(self%ntheta), icta(self%ntheta)
        integer, intent(in) :: iiph(self%nphih), ijph(self%nphih), ikph(self%nphih), ilph(self%nphih), icph(self%nphih)
        integer, intent(in) :: iipa(self%nphia), ijpa(self%nphia), ikpa(self%nphia), ilpa(self%nphia), icpa(self%nphia)

        integer :: ierr

        if (self%debug) then
            write (6, '(a)') '>>>>> entered qmmm_vsolv_store_parameters'
            call flush (6)
        end if

        self%prmtop_numbnd = numbnd

        self%iibh = iibh
        self%ijbh = ijbh
        self%icbh = icbh

        self%iiba = iiba
        self%ijba = ijba
        self%icba = icba

        self%iith = iith
        self%ijth = ijth
        self%ikth = ikth
        self%icth = icth

        self%iita = iita
        self%ijta = ijta
        self%ikta = ikta
        self%icta = icta

        self%iiph = iiph
        self%ijph = ijph
        self%ikph = ikph
        self%ilph = ilph
        self%icph = icph

        self%iipa = iipa
        self%ijpa = ijpa
        self%ikpa = ikpa
        self%ilpa = ilpa
        self%icpa = icpa

        if (self%debug) then
            write (6, '(a)') '<<<<< leaving qmmm_vsolv_store_parameters'
            call flush (6)
        end if

    end subroutine qmmm_vsolv_store_parameters

    subroutine new_qmmm_vsolv_type(self, nbonh, nbona, ntheth, ntheta, nphih, nphia)

        ! Note: fixed_iqmatoms, solvent_pointers, nearest_solvent_pointers
        !       and nearest_solvent_distances
        !       are all allocated manually in subroutine qmmm_vsolv_setup
        !       which is called from read_qmmm_nm_and_alloc()
        !       or in broadcast_qmmm_vsolv_type() below (for slave threads)
        implicit none

        type(qmmm_vsolv_type), intent(inout) :: self
        integer, intent(in) :: nbonh
        integer, intent(in) :: nbona
        integer, intent(in) :: ntheth
        integer, intent(in) :: ntheta
        integer, intent(in) :: nphih
        integer, intent(in) :: nphia

        integer :: ierr

        if (self%debug) then
            write (6, '(a)') '>>>>> entered new_qmmm_vsolv_type'
            call flush (6)
        end if

        self%nbonh = nbonh
        self%nbona = nbona
        self%ntheth = ntheth
        self%ntheta = ntheta
        self%nphih = nphih
        self%nphia = nphia

        ! bonds with Hydrogen
        if (.not. associated(self%iibh)) then
            allocate (self%iibh(nbonh), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ijbh)) then
            allocate (self%ijbh(nbonh), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%icbh)) then
            allocate (self%icbh(nbonh), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! bonds without hydrogen atoms
        if (.not. associated(self%iiba)) then
            allocate (self%iiba(nbona), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ijba)) then
            allocate (self%ijba(nbona), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%icba)) then
            allocate (self%icba(nbona), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! angles with hydrogen atoms
        if (.not. associated(self%iith)) then
            allocate (self%iith(ntheth), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ijth)) then
            allocate (self%ijth(ntheth), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ikth)) then
            allocate (self%ikth(ntheth), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%icth)) then
            allocate (self%icth(ntheth), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! angles without hydrogen atoms
        if (.not. associated(self%iita)) then
            allocate (self%iita(ntheta), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ijta)) then
            allocate (self%ijta(ntheta), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ikta)) then
            allocate (self%ikta(ntheta), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%icta)) then
            allocate (self%icta(ntheta), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! dihedrals with hydrogen atoms
        if (.not. associated(self%iiph)) then
            allocate (self%iiph(nphih), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ijph)) then
            allocate (self%ijph(nphih), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ikph)) then
            allocate (self%ikph(nphih), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ilph)) then
            allocate (self%ilph(nphih), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%icph)) then
            allocate (self%icph(nphih), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! dihedrals without hydrogen atoms
        if (.not. associated(self%iipa)) then
            allocate (self%iipa(nphia), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ijpa)) then
            allocate (self%ijpa(nphia), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ikpa)) then
            allocate (self%ikpa(nphia), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%ilpa)) then
            allocate (self%ilpa(nphia), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (.not. associated(self%icpa)) then
            allocate (self%icpa(nphia), stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (self%debug) then
            write (6, '(a)') '<<<<< leaving new_qmmm_vsolv_type'
            call flush (6)
        end if

    end subroutine new_qmmm_vsolv_type

    subroutine delete_qmmm_vsolv_type(self)

        implicit none

        type(qmmm_vsolv_type), intent(inout) :: self

        integer :: ierr

        if (self%debug) then
            write (6, '(a)') '>>>>> entered delete_qmmm_vsolv_type'
            call flush (6)
        end if

        if (associated(self%fixed_iqmatoms)) then
            deallocate (self%fixed_iqmatoms, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%solvent_pointers)) then
            deallocate (self%solvent_pointers, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%nearest_solvent_pointers)) then
            deallocate (self%nearest_solvent_pointers, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%nearest_solvent_distances)) then
            deallocate (self%nearest_solvent_distances, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! bonds with hydrogens
        if (associated(self%iibh)) then
            deallocate (self%iibh, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ijbh)) then
            deallocate (self%ijbh, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%icbh)) then
            deallocate (self%icbh, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! bonds without hydrogens
        if (associated(self%iiba)) then
            deallocate (self%iiba, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ijba)) then
            deallocate (self%ijba, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%icba)) then
            deallocate (self%icba, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! angles with hydrogens
        if (associated(self%iith)) then
            deallocate (self%iith, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ijth)) then
            deallocate (self%ijth, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ikth)) then
            deallocate (self%ikth, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%icth)) then
            deallocate (self%icth, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! angles without hydrogens
        if (associated(self%iita)) then
            deallocate (self%iita, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ijta)) then
            deallocate (self%ijta, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ikta)) then
            deallocate (self%ikta, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%icta)) then
            deallocate (self%icta, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! dihedrals with hydrogens
        if (associated(self%iiph)) then
            deallocate (self%iiph, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ijph)) then
            deallocate (self%ijph, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ikph)) then
            deallocate (self%ikph, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ilph)) then
            deallocate (self%ilph, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%icph)) then
            deallocate (self%icph, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        ! dihedrals without hydrogens
        if (associated(self%iipa)) then
            deallocate (self%iipa, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ijpa)) then
            deallocate (self%ijpa, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ikpa)) then
            deallocate (self%ikpa, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%ilpa)) then
            deallocate (self%ilpa, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (associated(self%icpa)) then
            deallocate (self%icpa, stat=ierr)
            REQUIRE(ierr == 0)
        end if

        if (self%debug) then
            write (6, '(a)') '<<<<< leaving delete_qmmm_vsolv_type'
            call flush (6)
        end if

    end subroutine delete_qmmm_vsolv_type

    subroutine print_qmmm_vsolv_type(self)

        implicit none

        type(qmmm_vsolv_type), intent(in) :: self

        if (self%debug) then
            write (6, '(a)') '>>>>> entered print_qmmm_vsolv_type'
            call flush (6)
        end if

        write (6, '(/,a)') 'QMMM VSOLV options:'

        write (6, '(3x,a,l5)') 'debug                        = ', self%debug
        write (6, '(3x,a,i5)') 'verbosity                    = ', self%verbosity
        write (6, '(3x,a,i5)') 'nearest_qm_solvent           = ', self%nearest_qm_solvent
        write (6, '(3x,a,i5)') 'nearest_qm_solvent_fq        = ', self%nearest_qm_solvent_fq
        write (6, '(3x,2a)') 'nearest_qm_solvent_resname   = ', trim(self%nearest_qm_solvent_resname)
        write (6, '(3x,a,i5)') 'nearest_qm_solvent_center_id = ', self%nearest_qm_solvent_center_id
        write (6, '(3x,a,i5)') 'qm_center_atom_id            = ', self%qm_center_atom_id

        ! Don't print the rest unless debug is on
        ! FIXIT: MOST OF THE REMAINING VARIABLES SHOULD BE MOVED TO ANOTHER MODULE ANYWAYS
        !        THAT IS A MODULE HOLDING THE TOPOLOGY DATA AND PARAMETERS
        if (.not. self%debug) return

        write (6, '(a,i5)') 'fixed_nquant                 = ', self%fixed_nquant
        write (6, '(a,i5)') 'nsolv_res                    = ', self%nsolv_res
        write (6, '(a,i5)') 'natom_solv_res               = ', self%natom_solv_res

        call print_integer_array('fixed_iqmatoms', self%fixed_iqmatoms)
        call print_integer_array('solvent_pointers', self%solvent_pointers)
        call print_integer_array('nearest_solvent_pointers', self%nearest_solvent_pointers)

        ! This is a hack:
        ! If iibh etc are not allocated, this print routine is called before these arrays
        ! are allocated. So return without printing if these are not allocated yet.
        if (.not. associated(self%iibh)) return

        write (6, '(a,i5)') 'prmtop_numbnd     = ', self%prmtop_numbnd
        write (6, '(a,i5)') 'nbonh             = ', self%nbonh
        write (6, '(a,i5)') 'nbona             = ', self%nbona
        write (6, '(a,i5)') 'ntheth            = ', self%ntheth
        write (6, '(a,i5)') 'ntheta            = ', self%ntheta
        write (6, '(a,i5)') 'nphih             = ', self%nphih
        write (6, '(a,i5)') 'nphia             = ', self%nphia

        call flush (6)
        call print_integer_array('iibh', self%iibh)
        call print_integer_array('ijbh', self%ijbh)
        call print_integer_array('icbh', self%icbh)
        call print_integer_array('iiba', self%iiba)
        call print_integer_array('ijba', self%ijba)
        call print_integer_array('icba', self%icba)
        call print_integer_array('iith', self%iith)
        call print_integer_array('ijth', self%ijth)
        call print_integer_array('ikth', self%ikth)
        call print_integer_array('icth', self%icth)
        call print_integer_array('iita', self%iita)
        call print_integer_array('ijta', self%ijta)
        call print_integer_array('ikta', self%ikta)
        call print_integer_array('icta', self%icta)
        call print_integer_array('iiph', self%iiph)
        call print_integer_array('ijph', self%ijph)
        call print_integer_array('ikph', self%ikph)
        call print_integer_array('ilph', self%ilph)
        call print_integer_array('icph', self%icph)
        call print_integer_array('iipa', self%iipa)
        call print_integer_array('ijpa', self%ijpa)
        call print_integer_array('ikpa', self%ikpa)
        call print_integer_array('ilpa', self%ilpa)
        call print_integer_array('icpa', self%icpa)

        if (self%debug) then
            write (6, '(a)') '<<<<< leavng print_qmmm_vsolv_type'
            call flush (6)
        end if

    end subroutine print_qmmm_vsolv_type

    subroutine print_integer_array(name, array)

        implicit none

        character(len=*), intent(in) :: name
        integer, dimension(:), intent(in) :: array

        integer :: i, j, jstart, jend
        integer :: jstep = 10

        write (6, '(a)') trim(name)//':'
        jstart = 1
        do i = 1, size(array)/jstep + 1
            jend = min((jstart + jstep - 1), size(array))
            write (6, '(10(i8))') (array(j), j=jstart, jend)
            jstart = jstart + jstep
        end do

    end subroutine print_integer_array

#ifdef MPI
    subroutine broadcast_qmmm_vsolv_type(self)

        implicit none
#include "parallel.h"
        include 'mpif.h'

        type(qmmm_vsolv_type), intent(inout) :: self

        integer :: ierr

        call mpi_bcast(self%debug, 1, mpi_logical, 0, commsander, ierr)

        if (self%debug) then
            write (6, '(a)') '>>>> entered broadcast'
            call flush (6)
        end if

        ! broadcast data that was set / allocated in qmmm_vsolv_setup()
        call mpi_bcast(self%debug, 1, mpi_logical, 0, commsander, ierr)
        call mpi_bcast(self%verbosity, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%nearest_qm_solvent, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%nearest_qm_solvent_fq, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%nearest_qm_solvent_resname, 4, mpi_character, 0, commsander, ierr)
        call mpi_bcast(self%nearest_qm_solvent_center_id, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%fixed_nquant, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%nsolv_res, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%natom_solv_res, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%qm_center_atom_id, 1, mpi_integer, 0, commsander, ierr)

        if (.not. associated(self%fixed_iqmatoms)) then
            allocate (self%fixed_iqmatoms(self%fixed_nquant), stat=ierr)
            REQUIRE(ierr == 0)
        end if
        if (.not. associated(self%solvent_pointers)) then
            allocate (self%solvent_pointers(self%nsolv_res), stat=ierr)
            REQUIRE(ierr == 0)
        end if
        if (.not. associated(self%nearest_solvent_pointers)) then
            allocate (self%nearest_solvent_pointers(self%nearest_qm_solvent + 1), stat=ierr)
            REQUIRE(ierr == 0)
        end if
        if (.not. associated(self%nearest_solvent_distances)) then
            allocate (self%nearest_solvent_distances(self%nearest_qm_solvent + 1), stat=ierr)
            REQUIRE(ierr == 0)
        end if
        call mpi_bcast(self%fixed_iqmatoms, self%fixed_nquant, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%solvent_pointers, self%nsolv_res, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%nearest_solvent_pointers, self%nearest_qm_solvent + 1, mpi_integer, 0, commsander, ierr)

        call mpi_bcast(self%nearest_solvent_distances, self%nearest_qm_solvent + 1, mpi_double_precision, 0, commsander, ierr)

        ! broadcast data that was allocated in new_qmmm_vsolv_type()
        ! and initialized in qmmm_vsolv_store_parameters()
        call mpi_bcast(self%nbonh, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%nbona, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ntheth, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ntheta, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%nphih, 1, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%nphia, 1, mpi_integer, 0, commsander, ierr)

        call new(self, self%nbonh, self%nbona, self%ntheth, self%ntheta, self%nphih, self%nphia)

        call mpi_bcast(self%prmtop_numbnd, 1, mpi_integer, 0, commsander, ierr)

        call mpi_bcast(self%iibh, self%nbonh, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ijbh, self%nbonh, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%icbh, self%nbonh, mpi_integer, 0, commsander, ierr)

        call mpi_bcast(self%iiba, self%nbona, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ijba, self%nbona, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%icba, self%nbona, mpi_integer, 0, commsander, ierr)

        call mpi_bcast(self%iith, self%ntheth, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ijth, self%ntheth, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ikth, self%ntheth, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%icth, self%ntheth, mpi_integer, 0, commsander, ierr)

        call mpi_bcast(self%iita, self%ntheta, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ijta, self%ntheta, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ikta, self%ntheta, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%icta, self%ntheta, mpi_integer, 0, commsander, ierr)

        call mpi_bcast(self%iiph, self%nphih, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ijph, self%nphih, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ikph, self%nphih, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ilph, self%nphih, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%icph, self%nphih, mpi_integer, 0, commsander, ierr)

        call mpi_bcast(self%iipa, self%nphia, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ijpa, self%nphia, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ikpa, self%nphia, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%ilpa, self%nphia, mpi_integer, 0, commsander, ierr)
        call mpi_bcast(self%icpa, self%nphia, mpi_integer, 0, commsander, ierr)

        if (self%debug) then
            call print_qmmm_vsolv_type(self)
            write (6, '(a)') '<<<<< leaving broadcast'
            call flush (6)
        end if

    end subroutine broadcast_qmmm_vsolv_type
#endif

end module qmmm_vsolv_module
