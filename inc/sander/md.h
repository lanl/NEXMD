#include "dprec.fh"
!-------------BEGIN    md.h  ------------------------------------------------

! NOTE: if you change this file, make sure you change the corresponding file
! in both src/sander/md.h and AmberTools/src/sqm/md.h

integer BC_MDI  ! size in integers of common block mdi
integer BC_MDR  ! size in Reals of common block mdr

! ... integer variables:

integer nrp,nspm,ig,ntx,ntcx,            &!5
      ntxo,ntt,ntp,ntr,init,             &!10
      ntcm,nscm,isolvp,nsolut,klambda,   &!15
      ntc,ntcc,ntf,ntid,ntn,             &!20
      ntnb,nsnb,ndfmin,nstlim,nrc,       &!25
      ntrx,npscal,imin,maxcyc,ncyc,      &!30
      ntmin,irest,jfastw,                &!33
      ibgwat,ienwat,iorwat,              &!36
      iwatpr,nsolw,igb,alpb,iyammp,      &!41
      gbsa,vrand,iwrap,nrespa,irespa,nrespai,icfe,       &!48
      rbornstat,ivcap,iconstreff,                        &!51
      neb,vv,tmode,ipol,iesp,ievb,nodeid,num_noshake,    &!59
      idecomp,icnstph,ntcnstph,maxdup,numexchg,repcrd,numwatkeep,hybridgb, &!67
      ibgion,ienion,profile_mpi,                         &!70
      ipb,inp,ntrelax,relaxing,dec_verbose,vdwmodel,     &!76
      csurften, ninterface, no_ntt3_sync                  !79

common/mdi/nrp,nspm,ig, &                                               !3
      ntx,ntcx,ntxo,ntt,ntp,ntr,init,ntcm,nscm, &                       !12
      isolvp,nsolut,ntc,ntcc,ntf,ntid,ntn,ntnb,nsnb,ndfmin, &           !22
      nstlim,nrc,ntrx,npscal,imin,maxcyc,ncyc,ntmin, &                  !30
      irest,jfastw,ibgwat,ienwat,iorwat, &                              !35
      iwatpr,nsolw,igb,alpb,iyammp,gbsa,vrand,numexchg,repcrd, &        !44
      numwatkeep,hybridgb, &                                            !46
      iwrap,nrespa,irespa,nrespai,icfe,rbornstat, &                     !52
      ivcap,iconstreff,idecomp,klambda,icnstph,ntcnstph,maxdup,neb,vv, &!61
      tmode,ipol,iesp,ievb,nodeid,num_noshake,ibgion,ienion, &          !69
      profile_mpi,ipb,inp,ntrelax,relaxing,dec_verbose,vdwmodel,  &     !76
      csurften, ninterface, no_ntt3_sync                                !79

parameter (BC_MDI=79) ! Number of elements in the common block;
                      ! Be sure to update if you change things

! ... floats:

_REAL_ t,dt,temp0,tautp,pres0,comp,taup,temp,tempi, & !9
      tol,taur,dx0,drms,vlimit,rbtarg(9),tmass,tmassinv,  & !25
      kappa,offset,surften,gamma_ln,extdiel,intdiel,rdt,  & !32
      gbalpha,gbbeta,gbgamma,cut_inner,clambda,saltcon,  & !38
      solvph,rgbmax,fsmax,restraint_wt, &  !42
      skmin,skmax,vfac,gbneckscale,v11,v12,v22,kevb,evbt,Arad, & !52
      gbalphaH,gbbetaH,gbgammaH, & !55 Hai Nguyen
      gbalphaC,gbbetaC,gbgammaC, & !58
      gbalphaN,gbbetaN,gbgammaN, & !61
      gbalphaOS,gbbetaOS,gbgammaOS, & !64 Hai Nguyen
      gbalphaP,gbbetaP,gbgammaP, &    !67
      Sh,Sc,Sn,So,Ss,Sp, &  ! 73 Hai Nguyen
      gamma_ten !74

common/mdr/t,dt,temp0,tautp,pres0,comp,taup,temp,tempi, &             !9
      tol,taur,dx0,drms,vlimit,rbtarg,tmass,tmassinv, &               !25
      kappa,offset,surften,gamma_ln,extdiel,intdiel,rdt, &            !32
      gbalpha,gbbeta,gbgamma,cut_inner,clambda,saltcon, &             !38
      solvph,rgbmax,fsmax,restraint_wt,skmin,skmax,vfac,gbneckscale, &!46
      v11,v12,v22,kevb,evbt,Arad, & !52
      gbalphaH,gbbetaH,gbgammaH,   & !55    !Hai Nguyen: add igb 8 parameters     
      gbalphaC,gbbetaC,gbgammaC, &   !58
      gbalphaN,gbbetaN,gbgammaN, &   !61
      gbalphaOS,gbbetaOS,gbgammaOS, & !64
      gbalphaP,gbbetaP,gbgammaP, & !67
      Sh,Sc,Sn,So,Ss,Sp, & !73
      gamma_ten !74

parameter (BC_MDR=74) ! Number of elements in the common block;
                      ! Be sure to update if you change things

! ... strings:

character(len=4) iwtnm,iowtnm,ihwtnm
  character(len=256) restraintmask,bellymask,tgtfitmask,&
            tgtrmsmask,noshakemask,crgmask,iwrap_mask
common/mds/ restraintmask,bellymask,tgtfitmask,tgtrmsmask,noshakemask,crgmask,  &
            iwtnm,iowtnm,ihwtnm(2),iwrap_mask

!-------------END    md.h  ------------------------------------------------

