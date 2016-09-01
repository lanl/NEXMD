module naesmd_module
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
      double precision temp0,tempf,tempi,tao
      integer,allocatable :: iordenhop(:)
      integer iorden(260)
      integer,target,allocatable:: atomtype(:)
      integer,allocatable:: lowvaluestep(:)
      double precision,allocatable:: lowvalue(:)
      double precision tini0
      double precision,target,allocatable:: rx(:),ry(:),rz(:)
      double precision,target,allocatable:: rxold(:),ryold(:),rzold(:)
      double precision,target,allocatable:: deltaxxpold(:),deltayypold(:)
      double precision,target,allocatable:: deltazzpold(:),deltaxxmold(:)
      double precision,target,allocatable:: deltayymold(:),deltazzmold(:)
      double precision,target,allocatable:: deltaxxpnew(:),deltayypnew(:)
      double precision,target,allocatable:: deltazzpnew(:),deltaxxmnew(:)
      double precision,target,allocatable:: deltayymnew(:),deltazzmnew(:)
      double precision,target,allocatable:: vx(:),vy(:),vz(:)
      double precision,target,allocatable:: vxold(:),vyold(:),vzold(:)
      double precision,target,allocatable:: ax(:),ay(:),az(:)
      double precision,target,allocatable:: axold(:),ayold(:),azold(:)
      double precision,allocatable:: fxmdqt(:),fymdqt(:),fzmdqt(:)
      double precision,target,allocatable::massmdqt(:)
      double precision dtquantum,kin
      double precision,target,allocatable::dtnact, dtmdqt
      double precision,allocatable:: vmdqt(:)
      double precision,allocatable:: bcoeffvmdqt(:)
      double precision,allocatable:: vmdqtnew(:)
      double precision,allocatable:: vmdqtmiddle(:)
      double precision,allocatable:: vmdqtmiddleold(:)
      double precision,allocatable:: vmdqtold(:)
      double precision,allocatable:: vnqcorrhoptot(:,:)
      double precision,allocatable:: vnqcorrhop(:,:)
      double precision,allocatable:: cadiab(:,:)
      double precision,allocatable:: cadiab_analt(:,:)
      double precision,allocatable:: cadiabnew(:,:)
      double precision,allocatable:: cadiabold(:,:)
      double precision,allocatable:: cadiabmiddle(:,:)
      double precision,allocatable:: cadiabmiddleold(:,:)
      double precision,allocatable:: bcoeffcadiab(:,:)
      double precision,allocatable:: cmdqt(:,:)
      double precision,allocatable:: cmdqtold(:,:)
      double precision,allocatable:: cmdqtmiddleold(:,:)
      double precision,allocatable:: cmdqtmiddle(:,:)
      double precision,allocatable:: cmdqtnew(:,:)
      double precision,allocatable:: scpr(:,:)
      double precision,allocatable:: cicoeffao2(:,:)
      double precision,allocatable:: uuold(:,:)
      double precision xcmini,ycmini,zcmini,masstot
      double precision tfemto
      double precision tfemtoquantum
      double precision nqold
      double precision kinini,vini,etotini
      double precision vgs
      double precision friction
      double precision,allocatable:: pfric(:),vfric(:),afric(:)
      double precision,allocatable:: prand(:,:),vrand(:,:)
      double precision deltax
      parameter(deltax=1.0d-4)
      character*2,allocatable:: atomtype2(:)
      character*4 state,prep
      character*6 ensemble
      character*200 cardini
      character*4 ktbig(0:9999)
      character*200 txtinput(1000)
      double precision cadiabhop
      double precision,allocatable:: scprreal(:,:)

      contains 
      subroutine allocate_naesmd_module(Na,Nexc,Nmo,Nbasis)
        allocate(iordenhop(Nexc))
        allocate(atomtype(Na))
        allocate(lowvaluestep(Nexc))
        allocate(lowvalue(Nexc))
        allocate(rx(Na),ry(Na),rz(Na))
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
        allocate(massmdqt(Na))
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
        allocate(cmdqt(Nbasis,Nexc))
        allocate(cmdqtold(Nbasis,Nexc))
        allocate(cmdqtmiddleold(Nbasis,Nexc))
        allocate(cmdqtmiddle(Nbasis,Nexc))
        allocate(cmdqtnew(Nbasis,Nexc))
        allocate(scpr(Nexc,Nexc))
        allocate(cicoeffao2(Nbasis,Nexc))
        allocate(uuold(Nmo,Nmo))
        allocate(pfric(Na),vfric(Na),afric(Na))
        allocate(prand(3,Na),vrand(3,Na))
        allocate(atomtype2(Na))
        allocate(scprreal(Nexc,Nexc))
      end

end module
