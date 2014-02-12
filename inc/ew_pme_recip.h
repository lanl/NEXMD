#include "dprec.fh"

! SIZES

#define BC_PME_PARS_INT 18
integer sizfftab,sizffwrk,siztheta,siz_q,sizheap,sizstack,sizscr
integer lfftable,lprefac1,lprefac2,lprefac3
integer order,nfft1,nfft2,nfft3,mlimit(3)
common/pme_size/ sizfftab,sizffwrk,siztheta,siz_q,sizheap,sizstack,sizscr, &
         lfftable,lprefac1,lprefac2,lprefac3, &
         order,nfft1,nfft2,nfft3,mlimit

#define BC_PME_PARS_REAL 4
_REAL_ dsum_tol,rsum_tol,maxexp,ew_coeff
common/pme_pars_real/dsum_tol,rsum_tol,maxexp,ew_coeff
