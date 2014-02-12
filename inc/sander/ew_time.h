#define EW_LST_TIM_NUM 10
#define EW_FRC_TIM_NUM 30
#include "dprec.fh"

! timer info

_REAL_ listtime,forctime
character(len=30) liststr
character(len=30) forcestr
_REAL_ time1,time2
common /ew_time/ listtime(EW_LST_TIM_NUM),forctime(EW_FRC_TIM_NUM) &
      ,time1,time2, &
      liststr(EW_LST_TIM_NUM),forcestr(EW_FRC_TIM_NUM)
integer &
      iself,iadjst,ibspl,ifillq,ifft,iscsum,igrsum, &
      iclear,iadjmap,irecsum,idirsum,iaccfrc,inbvir,ifrctim, &
      imap,isetgrd,igrduc,igrdim,ibldlst,ilsttim, &
      mfct1,mfct2,stack_time,mfct4,f_wait
common/indtime/ &
      iself,iadjst,ibspl,ifillq,ifft,iscsum,igrsum, &
      iclear,iadjmap,irecsum,idirsum,iaccfrc,inbvir,ifrctim, &
      imap,isetgrd,igrduc,igrdim,ibldlst,ilsttim, &
      mfct1,mfct2,stack_time,mfct4,f_wait
