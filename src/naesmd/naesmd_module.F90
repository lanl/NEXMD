#include "dprec.fh"
#include "assert.fh"

module naesmd_module
   
    implicit none
    
    private

    ! data types
    public :: naesmd_structure, realp_t, apc_common_struct
    
    ! functions and subroutines
    public :: allocate_naesmd_module_init
    public :: allocate_naesmd_module
    public :: allocate_naesmd_module2
    ! DATA Types
    type realp_t
        _REAL_, dimension(:), pointer :: p
    end type realp_t
    type apc_common_struct
	integer U(260),V(260),FB(260),RC(260),P(260) 
    end type
    
    type naesmd_structure
	    !data from main.f90
            integer crosstot
            integer, allocatable :: cross(:)
            _REAL_, dimension(:), allocatable :: yg, ygprime
            _REAL_, dimension(:), allocatable :: yg_new, ygprime_new
            _REAL_ :: pha !Nuclear phase for MCE (gamma)
            _REAL_, dimension(:,:), allocatable :: sgn
            integer lprint
	    
            !common data from apc.f90
	    type(apc_common_struct) :: apc_common 
	    !save data from function normal
	    _REAL_ store
            logical :: compute = .true.
	    !stolen from naesmd_data_t
	    _REAL_,allocatable:: Omega(:)
            _REAL_ E0
            integer Na,Nm
	    !sizes variables
	    integer natom,npot,ihop,nquantumstep,nstep,nstep0,nquantumreal
	    integer nstepcross,iredpot,nstates
	    integer icontw,nstepw,icontini
	    integer istepheat,iconttemperature,icontpdb,icont,nstepcoord
	    integer txtinicoord,txtendcoord
	    integer nbasis,nmaxbasis,nao
	    integer ihopprev
	    integer uumdqtflag
	    integer iview,jend,decorhop
	    integer iseedmdqt,conthop,conthop2
            integer fix
            integer nmc !number of normal modes to freeze
            integer, allocatable ::  mc(:) !indexes for freezing normal modes
            _REAL_ masatotal !total mass for freezing normal modes
            _REAL_,allocatable :: enm(:,:) !Normal modes read from fort.34 file 
            _REAL_,allocatable :: xbf0(:), ybf0(:), zbf0(:) !Equilibrium positions for freezing normal modes, read from tempA
            _REAL_,allocatable :: massqrt(:)
            _REAL_,allocatable :: w(:) !Gaussian widths for MCE
            _REAL_ temp0,tempf,tempi,tao
            integer npc !number of pairs of atoms to freeze the distance between them
            integer, allocatable :: dtc(:,:) !indexes of pairs of atoms to freeze distances between them
            integer printTdipole !to print transition dipole moments in separate files
            integer printTDM !to print the complete TDM
            _REAL_ :: S0_S1_threshold !to drop trajectory when S0/S1 gap is lower than 'S0_S1_threshold' (eV)
            integer,allocatable :: iordenhop(:),iorden(:),torden(:)
	    integer,allocatable:: atomtype(:) !atom types currently max 1000
	    integer,allocatable:: lowvaluestep(:)
	    _REAL_,allocatable:: lowvalue(:)
	    _REAL_ tini0,deltared
	    _REAL_,allocatable:: rx(:),ry(:),rz(:)
	    _REAL_,allocatable:: rxold(:),ryold(:),rzold(:)
	    _REAL_,allocatable:: deltaxxpold(:),deltayypold(:)
	    _REAL_,allocatable:: deltazzpold(:),deltaxxmold(:)
	    _REAL_,allocatable:: deltayymold(:),deltazzmold(:)
	    _REAL_,allocatable:: deltaxxpnew(:),deltayypnew(:)
	    _REAL_,allocatable:: deltazzpnew(:),deltaxxmnew(:)
	    _REAL_,allocatable:: deltayymnew(:),deltazzmnew(:)
	    _REAL_,allocatable:: vx(:),vy(:),vz(:)
	    _REAL_,allocatable:: vxold(:),vyold(:),vzold(:)
	    _REAL_,allocatable:: ax(:),ay(:),az(:)
	    _REAL_,allocatable:: axold(:),ayold(:),azold(:)
	    _REAL_,allocatable:: fxmdqt(:),fymdqt(:),fzmdqt(:)
	    _REAL_,allocatable::massmdqt(:)
	    _REAL_ dtquantum,kin
	    _REAL_ dtnact, dtmdqt
	    _REAL_,allocatable:: vmdqt(:)
	    _REAL_,allocatable:: bcoeffvmdqt(:)
	    _REAL_,allocatable:: vmdqtnew(:)
	    _REAL_,allocatable:: vmdqtmiddle(:)
	    _REAL_,allocatable:: vmdqtmiddleold(:)
	    _REAL_,allocatable:: vmdqtold(:)
	    _REAL_,allocatable:: vnqcorrhoptot(:,:)
	    _REAL_,allocatable:: vnqcorrhop(:,:)
	    _REAL_,allocatable:: cadiab(:,:)
	    _REAL_,allocatable:: cadiab_analt(:,:)
	    _REAL_,allocatable:: cadiabnew(:,:)
	    _REAL_,allocatable:: cadiabold(:,:)
	    _REAL_,allocatable:: cadiabmiddle(:,:)
	    _REAL_,allocatable:: cadiabmiddleold(:,:)
	    _REAL_,allocatable:: bcoeffcadiab(:,:)
	    _REAL_,allocatable:: cmdqt(:,:)
	    _REAL_,allocatable:: cmdqtold(:,:)
	    _REAL_,allocatable:: cmdqtmiddleold(:,:)
	    _REAL_,allocatable:: cmdqtmiddle(:,:)
	    _REAL_,allocatable:: cmdqtnew(:,:)
	    _REAL_,allocatable:: scpr(:,:)
	    _REAL_,allocatable:: cicoeffao2(:,:)
	    _REAL_,allocatable:: cicoeffao3(:)

	    _REAL_,allocatable:: uuold(:,:)
	    _REAL_ xcmini,ycmini,zcmini,masstot
	    _REAL_ tfemto
	    _REAL_ tfemtoquantum
	    _REAL_ nqold
	    _REAL_ kinini,vini,etotini,vmfini
	    _REAL_ vgs
	    _REAL_ friction
	    _REAL_,allocatable:: pfric(:),vfric(:),afric(:)
	    _REAL_,allocatable:: prand(:,:),vrand(:,:)
	    character*2,allocatable:: atomtype2(:)
	    character*4 state,prep
	    character*6 ensemble
	    character*200 cardini
	    character*10 dynam_type 
	    character*4 ktbig(0:9999)
	    character*200 txtinput(1000)
	    _REAL_ cadiabhop
	    _REAL_,allocatable:: scprreal(:,:)
            integer ifxd(1:1000)
  

	end type naesmd_structure
contains
    subroutine allocate_naesmd_module_init(naesmd_struct,natoms)
        implicit none
        integer natoms
        type(naesmd_structure), intent(inout) :: naesmd_struct 
        write(6,*)'Allocating initial naesmd_module variables'
        allocate(naesmd_struct%atomtype(natoms),naesmd_struct%atomtype2(natoms))
        allocate(naesmd_struct%massmdqt(natoms))
    end subroutine
    subroutine allocate_naesmd_module(naesmd_struct,Na,Nexc)
        implicit none
        integer Na,Nexc
        type(naesmd_structure), intent(inout) :: naesmd_struct
        write(6,*)'Allocating naesmd_module variables',Na,Nexc
        if(Nexc.ne.0) then
                allocate(naesmd_struct%iordenhop(Nexc),naesmd_struct%iorden(Nexc),naesmd_struct%torden(Nexc))
                allocate(naesmd_struct%lowvaluestep(Nexc))
                allocate(naesmd_struct%lowvalue(Nexc))
                allocate(naesmd_struct%vmdqt(Nexc))
                allocate(naesmd_struct%bcoeffvmdqt(Nexc))
                allocate(naesmd_struct%vmdqtnew(Nexc))
                allocate(naesmd_struct%vmdqtmiddle(Nexc))
                allocate(naesmd_struct%vmdqtmiddleold(Nexc))
                allocate(naesmd_struct%vmdqtold(Nexc))
                allocate(naesmd_struct%vnqcorrhoptot(Nexc,Nexc))
                allocate(naesmd_struct%vnqcorrhop(Nexc,Nexc))
                allocate(naesmd_struct%cadiab(Nexc,Nexc))
                allocate(naesmd_struct%cadiab_analt(Nexc,Nexc))
                allocate(naesmd_struct%cadiabnew(Nexc,Nexc))
                allocate(naesmd_struct%cadiabold(Nexc,Nexc))
                allocate(naesmd_struct%cadiabmiddle(Nexc,Nexc))
                allocate(naesmd_struct%cadiabmiddleold(Nexc,Nexc))
                allocate(naesmd_struct%bcoeffcadiab(Nexc,Nexc))
                allocate(naesmd_struct%scpr(Nexc,Nexc))
                allocate(naesmd_struct%scprreal(Nexc,Nexc))
        endif
        if(naesmd_struct%nmc.gt.0) then
            allocate(naesmd_struct%mc(naesmd_struct%nmc)) !Allocating array for indexes of normal modes to freeze
            allocate(naesmd_struct%enm(3*Na,3*Na)) !Allocating normal modes
            allocate(naesmd_struct%xbf0(Na),naesmd_struct%ybf0(Na),naesmd_struct%zbf0(Na))
            allocate(naesmd_struct%massqrt(Na))
        elseif(naesmd_struct%npc.gt.0) then
            allocate(naesmd_struct%dtc(naesmd_struct%npc,2))
        endif
        allocate(naesmd_struct%rx(Na),naesmd_struct%ry(Na),naesmd_struct%rz(Na))
        allocate(naesmd_struct%pfric(Na),naesmd_struct%vfric(Na),naesmd_struct%afric(Na))
        allocate(naesmd_struct%prand(3,Na),naesmd_struct%vrand(3,Na))
        allocate(naesmd_struct%rxold(Na),naesmd_struct%ryold(Na),naesmd_struct%rzold(Na))
        allocate(naesmd_struct%deltaxxpold(Na),naesmd_struct%deltayypold(Na))
        allocate(naesmd_struct%deltazzpold(Na),naesmd_struct%deltaxxmold(Na))
        allocate(naesmd_struct%deltayymold(Na),naesmd_struct%deltazzmold(Na))
        allocate(naesmd_struct%deltaxxpnew(Na),naesmd_struct%deltayypnew(Na))
        allocate(naesmd_struct%deltazzpnew(Na),naesmd_struct%deltaxxmnew(Na))
        allocate(naesmd_struct%deltayymnew(Na),naesmd_struct%deltazzmnew(Na))
        allocate(naesmd_struct%vx(Na),naesmd_struct%vy(Na),naesmd_struct%vz(Na))
        allocate(naesmd_struct%vxold(Na),naesmd_struct%vyold(Na),naesmd_struct%vzold(Na))
        allocate(naesmd_struct%ax(Na),naesmd_struct%ay(Na),naesmd_struct%az(Na))
        allocate(naesmd_struct%axold(Na),naesmd_struct%ayold(Na),naesmd_struct%azold(Na))
        allocate(naesmd_struct%fxmdqt(Na),naesmd_struct%fymdqt(Na),naesmd_struct%fzmdqt(Na))
        allocate(naesmd_struct%cicoeffao3(Na))
    end subroutine
    subroutine allocate_naesmd_module2(naesmd_struct,Ncis,Nmo,Nexc)
        implicit none
        type(naesmd_structure), intent(inout) :: naesmd_struct 
        integer Ncis,Nmo,Nexc
        allocate(naesmd_struct%cmdqt(Ncis,Nexc))
        allocate(naesmd_struct%cmdqtold(Ncis,Nexc))
        allocate(naesmd_struct%cmdqtmiddleold(Ncis,Nexc))
        allocate(naesmd_struct%cmdqtmiddle(Ncis,Nexc))
        allocate(naesmd_struct%cmdqtnew(Ncis,Nexc))
        allocate(naesmd_struct%cicoeffao2(Ncis,Nexc))
        allocate(naesmd_struct%uuold(Nmo,Nmo))
    end subroutine
end module
