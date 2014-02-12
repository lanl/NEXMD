!  ------------ begin ew_cntrl.h -----------------------------------------
!  control parameters:
!   verbose controls level of output. look in force_info
!   checkacc allows rigouous RMS force error checks
!   netfrc       = 1 if remove average force (due to analytic forces in pme)

!   ew_type      = 0 for pme,
!                  1 for reg ewald
!   vdwmeth      = 0 for cutoff of vdw,
!                  1 for analytic fix,
!                  2 for pme with geo mean mixing,
!                  3 for pme lorentz-berthelot mixing
!   use_pme      = 0 to skip reciprocal part of PME,
!                  1 to include it (default)

! M-WJ
!!   ipol is now specified in prmtop. YD
!   induced      = ipol (can be 0,1,2,3,4)
!   mpoltype     = induced
!                = 0 for no dipoles
!                = 1 for Applequist model
!                = 2 for exponential Thole model 1 ( as in AMOEBA )
!                = 3 for exponential Thole model 2
!                = 4 for linear Thole model
!



!   ***********************************************************
!   BC_EWCTNRL needs to be set to the size of the common block:

#define BC_EWCNTRL 15

!   ***********************************************************
integer verbose,netfrc,     ew_type,    vdwmeth, &
      periodic,  use_pme,    opt_infl,   ischrgd, fix_dip, &
      fix_quad,  mpoltype,   induced,    frameon, chngmask, &
      scaldip

common/ewcntrl/ &
      verbose,   netfrc,     ew_type,    vdwmeth,    &! 4
      periodic,  use_pme,    opt_infl,   ischrgd, fix_dip,    &!9
      fix_quad,  mpoltype,   induced,    frameon, chngmask,   &!14
      scaldip

logical nogrdptrs,nocutoff,boxbad
common/nogrd_flags/nogrdptrs,nocutoff,boxbad
#define BC_EWCNTRL_NP 2
integer inogrdptrs,inocutoff
common/nogrd_flags/inogrdptrs,inocutoff
!  ------------ end   ew_cntrl.h -----------------------------------------
