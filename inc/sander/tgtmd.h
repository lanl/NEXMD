!     note: if this common block if changed, dont forget to
!     update the initial broadcast in subroutine startup()
!     in parallel.f

!     some ntr stuff is broadcast as part of memory.h in NATOM block
!     ntr itself is in md.h in the NRP block

!       itgtmd    0 is default
!                 1 implies targeted md is to be used
!                  xc() will be the refc for RMSD calculation

!       dotgtmd    like "konst" logical for restrained (ntr=1) md

!       tgtrmsd    target rmsd value

!       tgtmdfrc    force constant for tmd


!        cannot be used with ntr=1 since refc shared
!        some info shared with variables normally used with ntr=1


integer itgtmd
common/tmd_int/ itgtmd

logical dotgtmd,rmsok

_REAL_ tgtrmsd,tgtmdfrc
common/tmd_real/ tgtrmsd,tgtmdfrc


! this does not need to be broadcast by parallel.f, it is done in runmd
! will need to be used in ene.f

_REAL_ rmsdvalue
common/tmd_real2/ rmsdvalue



