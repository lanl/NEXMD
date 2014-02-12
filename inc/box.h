!     common block sizes:

#define BC_BOXI 7
#define BC_BOXR 621

! ... floats:

#include "dprec.fh"
_REAL_ box,cut,dielc,rad,wel,radhb,welhb, &
      cutcap,xcap,ycap,zcap,fcap, &
      xlorth,ylorth,zlorth,xorth,yorth,zorth,forth, &
      rwell,xbox0,ybox0,zbox0
common/boxr/box(3),cut,dielc,xbox0,ybox0,zbox0, &       !8
      cutcap,xcap,ycap,zcap,fcap, &                     !13
      xlorth,ylorth,zlorth,xorth,yorth,zorth,forth, &   !20
      rwell,rad(100),wel(100),radhb(200),welhb(200)     !621

! ... integers:

integer ntb,ifbox,numpk,nbit,ifcap,natcap,isftrp
common/boxi/ntb,ifbox,numpk,nbit,ifcap,natcap,isftrp

_REAL_ extraboxdim
parameter (extraboxdim=30.d0)
