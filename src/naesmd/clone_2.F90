subroutine clone_naesmd(naesmd2, naesmd1, sim1)
use naesmd_module

implicit none
type(naesmd_structure),pointer            :: naesmd1, naesmd2
type(simulation_t), pointer  :: sim1

call allocate_naesmd_module_init(naesmd2, sim1%Na)  
call allocate_naesmd_module(naesmd2, sim1%Na, sim1%excN)  
if(sim1%excN.ne.0) then
	call allocate_naesmd_module2(naesmd2, sim1%dav%Ncis, sim1%nbasis, sim1%excN)  
endif

naesmd2%npot = naesmd1%npot 
if(naesmd2%npot.gt.0) then
	allocate(naesmd2%yg(2*sim1%excN))
    allocate(naesmd2%yg_new(3*sim1%excN))
	allocate(naesmd2%ygprime(2*sim1%excN))
    allocate(naesmd2%ygprime_new(3*sim1%excN))
	allocate(naesmd2%cross(sim1%excN))
	allocate(naesmd2%Omega(sim1%excN))
endif
naesmd2%crosstot=naesmd1%crosstot
naesmd2%cross(:) = naesmd1%cross(:) 
naesmd2%yg(:) = naesmd1%yg(:) 
naesmd2%yg_new(:) = naesmd1%yg_new(:)
naesmd2%sgn(:,:) = naesmd1%sgn(:,:)
naesmd2%ygprime(:) = naesmd1%ygprime(:) 
naesmd2%ygprime_new(:) = naesmd1%ygprime_new(:)
naesmd2%lprint = naesmd1%lprint
naesmd2%apc_common%U(:) = naesmd1%apc_common%U(:) 
naesmd2%apc_common%V(:) = naesmd1%apc_common%V(:) 
naesmd2%apc_common%FB(:) = naesmd1%apc_common%FB(:) 
naesmd2%apc_common%RC(:) = naesmd1%apc_common%RC(:) 
naesmd2%apc_common%P(:) = naesmd1%apc_common%P(:) 
naesmd2%store = naesmd1%store 
naesmd2%compute = naesmd1%compute 
naesmd2%Omega(:) = naesmd1%Omega(:) 
naesmd2%E0 = naesmd1%E0 
naesmd2%Na = naesmd1%Na 
naesmd2%Nm = naesmd1%Nm 
naesmd2%natom = naesmd1%natom 
naesmd2%ihop = naesmd1%ihop 
naesmd2%nquantumstep = naesmd1%nquantumstep 
naesmd2%nstep = naesmd1%nquantumreal 
naesmd2%nquantumreal = naesmd1%nquantumreal 
naesmd2%nstepcross = naesmd1%nstepcross 
naesmd2%iredpot = naesmd1%iredpot 
naesmd2%nstates = naesmd1%nstates 
naesmd2%icontw = naesmd1%icontw 
naesmd2%nstepw = naesmd1%nstepw 
naesmd2%icontini= naesmd1%icontini 
naesmd2%istepheat = naesmd1%istepheat 
naesmd2%iconttemperature = naesmd1%iconttemperature 
naesmd2%icontpdb = naesmd1%icontpdb 
naesmd2%icont = naesmd1%icont 
naesmd2%nstepcoord = naesmd1%nstepcoord 
naesmd2%txtinicoord = naesmd1%txtinicoord 
naesmd2%txtendcoord = naesmd1%txtendcoord 
naesmd2%nbasis = naesmd1%nbasis 
naesmd2%nmaxbasis = naesmd1%nmaxbasis 
naesmd2%nao = naesmd1%nao 
naesmd2%ihopprev = naesmd1%ihopprev 
naesmd2%uumdqtflag = naesmd1%uumdqtflag 
naesmd2%iview = naesmd1%iview 
naesmd2%jend = naesmd1%jend 
naesmd2%decorhop= naesmd1%decorhop 
naesmd2%iseedmdqt = naesmd1%iseedmdqt 
naesmd2%conthop = naesmd1%conthop 
naesmd2%conthop2 = naesmd1%conthop2 
naesmd2%temp0 = naesmd1%temp0 
naesmd2%tempf = naesmd1%tempf 
naesmd2%tempi = naesmd1%tempi 
naesmd2%tao = naesmd1%tao 
naesmd2%nmc = naesmd1%nmc
naesmd2%npc = naesmd1%npc
naesmd2%mc(:) = naesmd1%mc(:)
naesmd2%enm(:,:) = naesmd1%enm(:,:)
naesmd2%printTdipole = naesmd1%printTdipole
naesmd2%xbf0(:) = naesmd1%xbf0(:)
naesmd2%ybf0(:) = naesmd1%ybf0(:)
naesmd2%zbf0(:) = naesmd1%zbf0(:)
naesmd2%massqrt(:) = naesmd1%massqrt(:)
naesmd2%iordenhop(:) = naesmd1%iordenhop(:) 
naesmd2%iorden(:) = naesmd1%iorden(:) 
naesmd2%atomtype(:) = naesmd1%atomtype(:) 
naesmd2%lowvaluestep(:) = naesmd1%lowvaluestep(:) 
naesmd2%lowvalue(:) = naesmd1%lowvalue(:) 
naesmd2%tini0 = naesmd1%tini0 
naesmd2%deltared = naesmd1%deltared 
naesmd2%rx(:) = naesmd1%rx(:) 
naesmd2%ry(:) = naesmd1%ry(:) 
naesmd2%rz(:) = naesmd1%rz(:) 
naesmd2%rxold(:) = naesmd1%rxold(:) 
naesmd2%ryold(:) = naesmd1%ryold(:) 
naesmd2%rzold(:) = naesmd1%rzold(:) 
naesmd2%deltaxxpold(:) = naesmd1%deltaxxpold(:) 
naesmd2%deltayypold(:) = naesmd1%deltayypold(:) 
naesmd2%deltazzpold(:) = naesmd1%deltazzpold(:)
naesmd2%deltaxxmold(:) = naesmd1%deltaxxmold(:) 
naesmd2%deltayymold(:) = naesmd1%deltayymold(:)
naesmd2%deltazzmold(:) = naesmd1%deltazzmold(:) 
naesmd2%deltaxxpnew(:) = naesmd1%deltaxxpnew(:)
naesmd2%deltayypnew(:) = naesmd1%deltayypnew(:)
naesmd2%deltazzpnew(:) = naesmd1%deltazzpnew(:)
naesmd2%deltaxxmnew(:) = naesmd1%deltaxxmnew(:)
naesmd2%deltayymnew(:) = naesmd1%deltayymnew(:) 
naesmd2%deltazzmnew(:) = naesmd1%deltazzmnew(:)
naesmd2%vx(:) = naesmd1%vx(:)
naesmd2%vy(:) = naesmd1%vy(:)
naesmd2%vz(:) = naesmd1%vz(:)
naesmd2%vxold(:) = naesmd1%vxold(:)
naesmd2%vyold(:) = naesmd1%vyold(:)
naesmd2%vzold(:) = naesmd1%vzold(:)
naesmd2%ax(:)= naesmd1%ax(:)
naesmd2%ay(:) = naesmd1%ay(:)
naesmd2%az(:) = naesmd1%az(:)
naesmd2%axold(:) = naesmd1%axold(:)
naesmd2%ayold(:) = naesmd1%ayold(:)
naesmd2%azold(:) = naesmd1%azold(:)
naesmd2%fxmdqt(:) = naesmd1%fxmdqt(:)
naesmd2%fymdqt(:) = naesmd1%fymdqt(:)
naesmd2%fzmdqt(:) = naesmd1%fzmdqt(:)
naesmd2%massmdqt(:) = naesmd1%massmdqt(:)
naesmd2%dtquantum = naesmd1%dtquantum
naesmd2%kin = naesmd1%kin
naesmd2%dtnact = naesmd1%dtnact 
naesmd2%dtmdqt = naesmd1%dtmdqt
naesmd2%vmdqt(:) = naesmd1%vmdqt(:)
naesmd2%bcoeffvmdqt(:) = naesmd1%bcoeffvmdqt(:)
naesmd2%vmdqtnew(:) = naesmd1%vmdqtnew(:)
naesmd2%vmdqtmiddle(:) = naesmd1%vmdqtmiddle(:)
naesmd2%vmdqtmiddleold(:) = naesmd1%vmdqtmiddleold(:)
naesmd2%vmdqtold(:) = naesmd1%vmdqtold(:)
naesmd2%vnqcorrhoptot(:,:) = naesmd1%vnqcorrhoptot(:,:)
naesmd2%vnqcorrhop(:,:) = naesmd1%vnqcorrhop(:,:)
naesmd2%cadiab(:,:) = naesmd1%cadiab(:,:)
naesmd2%cadiab_analt(:,:) = naesmd1%cadiab_analt(:,:)
naesmd2%cadiabnew(:,:) = naesmd1%cadiabnew(:,:)
naesmd2%cadiabold(:,:) = naesmd1%cadiabold(:,:)
naesmd2%adiabmiddle(:,:) = naesmd1%adiabmiddle(:,:)
naesmd2%cadiabmiddleold(:,:) = naesmd1%cadiabmiddleold(:,:)
naesmd2%bcoeffcadiab(:,:) = naesmd1%bcoeffcadiab(:,:)
naesmd2%cmdqt(:,:) = naesmd1%cmdqt(:,:)
naesmd2%cmdqtold(:,:) = naesmd1%cmdqtold(:,:)
naesmd2%cmdqtmiddleold(:,:) = naesmd1%cmdqtmiddleold(:,:)
naesmd2%cmdqtmiddle(:,:) = naesmd1%cmdqtmiddle(:,:)
naesmd2%cmdqtnew(:,:) = naesmd1%cmdqtnew(:,:)
naesmd2%scpr(:,:) = naesmd1%scpr(:,:)
naesmd2%cicoeffao2(:,:) = naesmd1%cicoeffao2(:,:)
naesmd2%cicoeffao3(:)  = naesmd1%cicoeffao3(:)
naesmd2%uuold(:,:) = naesmd1%uuold(:,:)
naesmd2%xcmini = naesmd1%xcmini
naesmd2%ycmini = naesmd1%ycmini 
naesmd2%zcmini = naesmd1%zcmini
naesmd2%masstot = naesmd1%masstot
naesmd2%tfemto = naesmd1%tfemto
naesmd2%tfemtoquantum = naesmd1%tfemtoquantum 
naesmd2%nqold = naesmd1%nqold
naesmd2%kinini = naesmd1%kinini
naesmd2%vini = naesmd1%vini
naesmd2%etotini = naesmd1%etotini
naesmd2%vmfini = naesmd1%vmfini
naesmd2%vgs = naesmd1%vgs
naesmd2%friction = naesmd1%friction 
naesmd2%pfric(:) = naesmd1%pfric(:)
naesmd2%vfric(:) = naesmd1%vfric(:)
naesmd2%afric(:) = naesmd1%afric(:)
naesmd2%prand(:,:) = naesmd1%prand(:,:)
naesmd2%vrand(:,:) = naesmd1%vrand(:,:)
naesmd2%atomtype2(:) = naesmd1%atomtype2(:)
naesmd2%state = naesmd1%state
naesmd2%prep = naesmd1%prep
naesmd2%ensemble = naesmd1%ensemble 
naesmd2%cardini = naesmd1%cardini
naesmd2%dynam_type = naesmd1%dynam_type
naesmd2%ktbig(:) = naesmd1%ktbig(:)
naesmd2%txtinput(:) = naesmd1%txtinput(:)
naesmd2%cadiabhop = naesmd1%cadiabhop
naesmd2%scprreal(:,:) = naesmd1%scprreal(:,:)
naesmd2%ifxd(:) = naesmd1%ifxd(:)
!naesmd2% = naesmd1% 

end subroutine clone_naesmd

subroutine clone_md(md2, md1, sim1)
use md_module

implicit none
type(naesmd_structure),pointer            :: md1, md2
type(simulation_t), pointer  :: sim1

call allocate_md_module_init(md2,sim1%Na)
call allocate_md_module(md2,sim1%Na)

md2%atoms(:) = md1%atoms(:)
md2%atmass(:) = md1%atmass(:)
md2%v(:,:) = md1%v(:,:)
md2%fm(:) = md1%fm(:)
md2%fo(:) = md1%fo(:)
md2%r0(:) = md1%r0(:)
md2%ttt = md1%ttt
md2%imdtype = md1%imdtype
md2%ideriv = md1%ideriv
md2%icart = md1%icart
md2%ifric = md1%ifric

end subroutine clone_md



subroutine clone_dav(dav2, dav1, sim1)
use qm2_davidson_module

!must be called after clone of qmmm_struct qmmm_scratch, qmmm_nml, qm2_struct

implicit none
type(naesmd_structure),pointer            :: dav1, dav2
integer natoms
type(simulation_t), pointer  :: sim1

call allocate_davidson(sim1%qmmm_scratch, sim

end subroutine clone_dav


end module

