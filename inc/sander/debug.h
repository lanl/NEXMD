integer natomn
parameter (natomn = 25)
integer do_debugf,atomn(natomn),neglgdel,zerochg,zerovdw, &
      zerodip,nranatm,ranseed,chkvir,dumpfrc,rmsfrc
#define BC_DEBUG 35
common/debug/do_debugf,neglgdel,zerochg,zerovdw,atomn, &
      zerodip,nranatm,ranseed,chkvir,dumpfrc,rmsfrc
! memory offsets
#define BC_DEB_HEAP 12
integer lscg,lsrms,lscn1,lscn2,lsdip,lsind, &
      lf1,lf2,lf3,lf4,lf5,lf6
common/deb_heap/lscg,lsrms,lscn1,lscn2,lsdip, &
      lsind,lf1,lf2,lf3,lf4,lf5,lf6
