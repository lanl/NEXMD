integer, parameter :: mshfd=12000, maxfe=5, &
      mpar=5, mp=maxfe*mpar+1, np=maxfe*mpar

_REAL_ :: &
      coox(mshfd),cooy(mshfd),cooz(mshfd), &
      cmx(maxfe),cmy(maxfe),cmz(maxfe), &
      tolpro(mshfd),wt(mshfd),shift(mshfd), &
      obs(mshfd), obsp(mshfd),shavp(mshfd), &
      phi(maxfe),teta(maxfe),omg(maxfe), &
      a1dip(maxfe),a2dip(maxfe), &
      optphi(maxfe),opttet(maxfe),optomg(maxfe), &
      opta1(maxfe),opta2(maxfe), &
      toldip,resid,oldresid,optkon

integer :: &
      mltpro(mshfd),iprot(mshfd), &
      ippmc(maxfe),ioldvio,iviolation, &
      nhp,nfe,nprot,ipear,nstampa,nstampaten

common /simplex/ &
      coox,cooy,cooz, &
      cmx,cmy,cmz, &
      tolpro,wt,shift, &
      obs,mltpro,iprot, &
      obsp,shavp, &
      phi,teta,omg, &
      a1dip,a2dip, &
      optphi,opttet,optomg, &
      opta1,opta2, &
      toldip,resid,oldresid,optkon, &
      ippmc,ioldvio,iviolation, &
      nhp,nfe,nprot,ipear,nstampa,nstampaten

_REAL_ :: d(3,mshfd)
common/derive/ d

