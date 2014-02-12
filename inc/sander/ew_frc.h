!---------------------- ew_frc.h --------------------
#include "dprec.fh"

_REAL_ eer,eed,evdw,evdwr,ehb, &
      eedvir,eea,ees,epold,epola,epols
_REAL_ dipself,dipkine,diprms,dipndf
_REAL_ rec_vir(3,3),dir_vir(3,3), &
      adj_vir(3,3),rec_vird(3,3),self_vir(3,3)
_REAL_ atvir(3,3),molvir(3,3),subvir(3,3)
_REAL_ c1,c2,c3,xr1,xr2,xr3
_REAL_ ee14,enb14,epol14
_REAL_ e14vir(3,3),framevir(3,3)

#define BC_EW_COMM3 114
common/ew_comm3/eer,eed,evdw,evdwr,ehb,eedvir,eea,ees, &
      epold,epola,epols, &
      dipself,dipkine,diprms,dipndf, &
      ee14,enb14,epol14, &
      e14vir,framevir, &
      rec_vir,dir_vir,adj_vir,rec_vird,self_vir, &
      c1,c2,c3,xr1,xr2,xr3, &
      atvir,molvir,subvir
!-------------------END ew_frc.h --------------------

