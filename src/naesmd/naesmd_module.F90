#include "dprec.fh"
#include "assert.fh"

module naesmd_module
    !sizes variables
    integer natom,npot,ihop,nquantumstep,nstep,nquantumreal
    integer nstepcross
    integer icontw,nstepw,icontini
    integer istepheat,iconttemperature,icontpdb,icont,nstepcoord
    integer txtinicoord,txtendcoord
    integer nbasis,nmaxbasis,nao
    integer ihopprev
    integer uumdqtflag
    integer iview,jend,decorhop
    integer iseedmdqt,conthop,conthop2
    _REAL_ temp0,tempf,tempi,tao
    integer,allocatable :: iordenhop(:),iorden(:)
    integer,target,allocatable:: atomtype(:) !atom types currently max 1000
    integer,allocatable:: lowvaluestep(:)
    _REAL_,allocatable:: lowvalue(:)
    _REAL_ tini0
    _REAL_,target,allocatable:: rx(:),ry(:),rz(:)
    _REAL_,target,allocatable:: rxold(:),ryold(:),rzold(:)
    _REAL_,target,allocatable:: deltaxxpold(:),deltayypold(:)
    _REAL_,target,allocatable:: deltazzpold(:),deltaxxmold(:)
    _REAL_,target,allocatable:: deltayymold(:),deltazzmold(:)
    _REAL_,target,allocatable:: deltaxxpnew(:),deltayypnew(:)
    _REAL_,target,allocatable:: deltazzpnew(:),deltaxxmnew(:)
    _REAL_,target,allocatable:: deltayymnew(:),deltazzmnew(:)
    _REAL_,target,allocatable:: vx(:),vy(:),vz(:)
    _REAL_,target,allocatable:: vxold(:),vyold(:),vzold(:)
    _REAL_,target,allocatable:: ax(:),ay(:),az(:)
    _REAL_,target,allocatable:: axold(:),ayold(:),azold(:)
    _REAL_,allocatable:: fxmdqt(:),fymdqt(:),fzmdqt(:)
    _REAL_,target,allocatable::massmdqt(:)
    _REAL_ dtquantum,kin
    _REAL_,target,allocatable::dtnact, dtmdqt
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
    _REAL_,allocatable:: uuold(:,:)
    _REAL_ xcmini,ycmini,zcmini,masstot
    _REAL_ tfemto
    _REAL_ tfemtoquantum
    _REAL_ nqold
    _REAL_ kinini,vini,etotini
    _REAL_ vgs
    _REAL_ friction
    _REAL_,allocatable:: pfric(:),vfric(:),afric(:)
    _REAL_,allocatable:: prand(:,:),vrand(:,:)
    _REAL_ deltax
    parameter(deltax=1.0d-4)
    character*2,allocatable:: atomtype2(:)
    character*4 state,prep
    character*6 ensemble
    character*200 cardini
    character*4 ktbig(0:9999)
    character*200 txtinput(1000)
    _REAL_ cadiabhop
    _REAL_,allocatable:: scprreal(:,:)

contains
    subroutine allocate_naesmd_module_init(natoms)
        implicit none
        integer natoms
        write(6,*)'Allocating initial naesmd_module variables'
        allocate(atomtype(natoms),atomtype2(natoms))
        allocate(massmdqt(natoms))
    end
    subroutine allocate_naesmd_module(Na,Nexc)
        implicit none
        integer Na,Nexc
        write(6,*)'Allocating naesmd_module variables',Na,Nexc
        if(Nexc.ne.0) then
                allocate(iordenhop(Nexc),iorden(Nexc))
                allocate(lowvaluestep(Nexc))
                allocate(lowvalue(Nexc))
                allocate(vmdqt(Nexc))
                allocate(bcoeffvmdqt(Nexc))
                allocate(vmdqtnew(Nexc))
                allocate(vmdqtmiddle(Nexc))
                allocate(vmdqtmiddleold(Nexc))
                allocate(vmdqtold(Nexc))
                allocate(vnqcorrhoptot(Nexc,Nexc))
                allocate(vnqcorrhop(Nexc,Nexc))
                allocate(cadiab(Nexc,Nexc))
                allocate(cadiab_analt(Nexc,Nexc))
                allocate(cadiabnew(Nexc,Nexc))
                allocate(cadiabold(Nexc,Nexc))
                allocate(cadiabmiddle(Nexc,Nexc))
                allocate(cadiabmiddleold(Nexc,Nexc))
                allocate(bcoeffcadiab(Nexc,Nexc))
                allocate(scpr(Nexc,Nexc))
                allocate(scprreal(Nexc,Nexc))
        endif
        allocate(rx(Na),ry(Na),rz(Na))
        allocate(pfric(Na),vfric(Na),afric(Na))
        allocate(prand(3,Na),vrand(3,Na))
        allocate(rxold(Na),ryold(Na),rzold(Na))
        allocate(deltaxxpold(Na),deltayypold(Na))
        allocate(deltazzpold(Na),deltaxxmold(Na))
        allocate(deltayymold(Na),deltazzmold(Na))
        allocate(deltaxxpnew(Na),deltayypnew(Na))
        allocate(deltazzpnew(Na),deltaxxmnew(Na))
        allocate(deltayymnew(Na),deltazzmnew(Na))
        allocate(vx(Na),vy(Na),vz(Na))
        allocate(vxold(Na),vyold(Na),vzold(Na))
        allocate(ax(Na),ay(Na),az(Na))
        allocate(axold(Na),ayold(Na),azold(Na))
        allocate(fxmdqt(Na),fymdqt(Na),fzmdqt(Na))
    end
    subroutine allocate_naesmd_module2(Ncis,Nmo,Nexc)
        implicit none
        integer Ncis,Nmo,Nexc
        allocate(cmdqt(Ncis,Nexc))
        allocate(cmdqtold(Ncis,Nexc))
        allocate(cmdqtmiddleold(Ncis,Nexc))
        allocate(cmdqtmiddle(Ncis,Nexc))
        allocate(cmdqtnew(Ncis,Nexc))
        allocate(cicoeffao2(Ncis,Nexc))
        allocate(uuold(Nmo,Nmo))
    end
end module
