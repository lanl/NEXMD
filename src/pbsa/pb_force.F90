! <compile=optimized>
#include "copyright.h"
#define _REAL_ double precision
#define REQUIRE(e) if(.not.(e)) call croak(__FILE__,__LINE__)
#include "pb_def.h"
#include "timer.h"

module poisson_boltzmann

   implicit none

#  include "pb_constants.h"

   ! PBMD parameters

   _REAL_, parameter :: pbkb   = 1.3807D-23 / 1.6606D-27 / (1.00D+12)**2 * (1.00D+10)**2
   _REAL_, parameter :: fioni  = 6.0220D+23 / 1.00D+30
   _REAL_, parameter :: fiono  = ONE / fioni
   _REAL_, parameter :: eps0   = 8.8542D-12 / (1.6022D-19)**2 / (1.00D+10)**3 * (1.00D+12)**2 * 1.6606D-27
   _REAL_, parameter :: frcfac = 0.01D0 / 4.1840D0

   ! PBMD FD control variables

   logical :: outphi
   logical :: srsas
   logical :: scalerf
   logical :: outlvlset    !WMBS output total lvlset 
   logical :: outmlvlset   !WMBS output membrane lvlset

  
   integer :: phiform 
   integer :: saopt
   integer :: sasopt
   integer :: dbfopt
   integer :: eneopt
   integer :: npbopt
   integer :: solvopt
   integer :: frcopt
   integer :: intopt
   integer :: bcopt
   integer :: smoothopt
   integer :: xm
   integer :: ym
   integer :: zm
   integer :: xmym
   integer :: xmymzm
   integer :: nbuffer
   integer :: level
   integer :: nfocus
   integer :: fscale
   integer :: maxitn
   integer :: itn
   integer :: m, n
   integer :: savbcopt(MAXLEVEL)
   integer :: levelblock(MAXLEVEL)
   integer :: savxm(MAXLEVEL)
   integer :: savym(MAXLEVEL)
   integer :: savzm(MAXLEVEL)
   integer :: savxo(MAXLEVEL) !Origin offset
   integer :: savyo(MAXLEVEL) !
   integer :: savzo(MAXLEVEL)
   integer :: savxmym(MAXLEVEL)
   integer :: savxmymzm(MAXLEVEL)
   integer :: isurfchg           ! print surface charges
   integer :: membraneopt        !WMBS - membrane type 0 for none 1 for slab
   !integer :: mdielectricopt    !WMBS - place holder for future use
   integer :: poretype   !WMBS- membrane exclusion region flag
   integer :: augtoltype         !WMBS - Aug IIM gmres Tolerance type
                                 !       0 for absolute, 1 for relative

   _REAL_ :: h
   _REAL_ :: gox
   _REAL_ :: goy
   _REAL_ :: goz
   _REAL_ :: fmiccg
   _REAL_ :: accept
   _REAL_ :: laccept
   _REAL_ :: wsor
   _REAL_ :: lwsor
   _REAL_ :: norm
   _REAL_ :: inorm
   _REAL_ :: xmax
   _REAL_ :: xmin
   _REAL_ :: ymax
   _REAL_ :: ymin
   _REAL_ :: zmax
   _REAL_ :: zmin
   _REAL_ :: gxmax
   _REAL_ :: gxmin
   _REAL_ :: gymax
   _REAL_ :: gymin
   _REAL_ :: gzmax
   _REAL_ :: gzmin
   _REAL_ :: savxbox(MAXLEVEL)
   _REAL_ :: savybox(MAXLEVEL)
   _REAL_ :: savzbox(MAXLEVEL)
   _REAL_ :: cxbox(MAXLEVEL)
   _REAL_ :: cybox(MAXLEVEL)
   _REAL_ :: czbox(MAXLEVEL)
   _REAL_ :: savh(MAXLEVEL)
   _REAL_ :: savgox(MAXLEVEL)
   _REAL_ :: savgoy(MAXLEVEL)
   _REAL_ :: savgoz(MAXLEVEL)
   _REAL_ :: offx
   _REAL_ :: offy
   _REAL_ :: offz
   _REAL_ :: fillratio

   _REAL_ :: epsin
   _REAL_ :: epsout
   _REAL_ :: pbkappa
   _REAL_ :: istrng
   _REAL_ :: ivalence
   _REAL_ :: pbtemp
   _REAL_ :: totcrg
   _REAL_ :: totcrgp
   _REAL_ :: totcrgn

   _REAL_ :: pbgamma_int
   _REAL_ :: pbgamma_ext
   _REAL_ :: qef(3)
   _REAL_ :: ref(3)

   _REAL_ :: mthick     !WMBS- membrane thickness
   _REAL_ :: mctrdz     !WMBS- membrane z offset
   !_REAL_ :: mdielectric !WMBS placeholder for later
   _REAL_ :: poreradius !WMBS- radius for cylindrical region

   _REAL_ :: augctf     !WMBS- Aug IIM cluster radius (0 for h*h)
   _REAL_ :: augtol     !WMBS- Aug IIM gmres tolerance

   ! PBMD topology information

   integer              :: lastp
   integer              :: ngrdcrg

   integer, allocatable ::    icrd(:,:)
   integer, allocatable ::  grdcrg(:,:)
   _REAL_, allocatable :: qgrdcrg(:)
   _REAL_, allocatable ::    gcrd(:,:)
   _REAL_, allocatable ::    acrd(:,:)
   _REAL_, allocatable ::    acrg(:)
   _REAL_, allocatable ::    gcrg(:,:)
 
   ! PBMD nblist information

   integer              :: maxnbr
   integer              :: maxnba
   _REAL_              :: cutres, cutnb, cutfd, cutsa
 
   integer, allocatable ::   nshrt(:)
   integer, allocatable ::     nex(:)
   integer, allocatable ::     iex(:,:)
   integer, allocatable :: iprlong(:)
   integer, allocatable :: iprshrt(:)
   integer, allocatable ::  iar1pb(:,:)
   _REAL_, allocatable ::   cn1pb(:)
   _REAL_, allocatable ::   cn2pb(:)
   _REAL_, allocatable ::   cn3pb(:)

   ! PBMD cap water simulation information

   integer              :: mpopt
   integer              :: lmax
   integer              :: inatm
   integer              :: outwat
   integer              :: oution
   integer, allocatable :: outflag(:)
   integer, allocatable :: outflagorig(:)
   integer, allocatable :: mapout(:)
   integer, allocatable :: ibelly(:)
   _REAL_              :: sepbuf

   ! physical variables for energy and force calculations

   integer:: nbnd
   integer:: nbndx
   integer:: nbndy
   integer:: nbndz
   _REAL_, allocatable :: pos_crg(:,:,:)
   _REAL_, allocatable :: surf_crg(:,:)
   integer, allocatable :: ipos_crg(:,:,:)
   integer, allocatable :: crg_num(:)

   ! physical variable maps for numerical solutions

   _REAL_, allocatable ::     phi(:)
   _REAL_, allocatable ::      bv(:)
   _REAL_, allocatable ::   chgrd(:)
   _REAL_, allocatable ::    epsx(:)
   _REAL_, allocatable ::    epsy(:)
   _REAL_, allocatable ::    epsz(:)
   _REAL_, allocatable :: saltgrd(:)

   ! geometry maps for dielectric interface
 
   integer, allocatable ::  insas(:)
   integer, allocatable :: atmsas(:)
   _REAL_, allocatable ::  lvlset(:)
   _REAL_, allocatable ::      zv(:)

   ! physical variable maps for force calculations

   _REAL_, allocatable ::     cphi(:)
   integer, allocatable ::  iepsav(:,:)
   integer, allocatable :: iepsavx(:,:)
   integer, allocatable :: iepsavz(:,:)
   integer, allocatable :: iepsavy(:,:)
   _REAL_, allocatable ::   fedgex(:)
   _REAL_, allocatable ::   fedgey(:)
   _REAL_, allocatable ::   fedgez(:)

   ! saved phi array for pbmd
 
   _REAL_, allocatable :: xs(:)
   integer :: xsoffset

   ! ligand focusing options

   logical :: ligand
   character(len=256) ligandmask
   integer, allocatable :: liveflag(:)
   integer, allocatable :: realflag(:)
   integer :: ntrajmol
   _REAL_ :: buffer

   ! Multiple distributive fine grid geometry / Multiblock focusing

   logical :: multiblock                          !TRUE if specified multiblock
   logical :: firstleveldone                      !TRUE if level 1 is done once
   integer, allocatable :: blkxo(:)               !block origin in x dir
   integer, allocatable :: blkyo(:)               !block origin in y dir
   integer, allocatable :: blkzo(:)               !block origin in z dir
   integer, allocatable :: blkxlo(:)              !block lower bound in x
   integer, allocatable :: blkylo(:)              !block lower bound in y
   integer, allocatable :: blkzlo(:)              !block lower bound in z
   integer, allocatable :: blkxup(:)              !block upper bound in x
   integer, allocatable :: blkyup(:)              !block upper bound in y
   integer, allocatable :: blkzup(:)              !block upper bound in z
   integer, allocatable :: blknx(:)               !block grid length in x
   integer, allocatable :: blkny(:)               !block grid length in y
   integer, allocatable :: blknz(:)               !block grid length in z
   integer, allocatable :: blknxny(:)             !blknx . blkny
   integer, allocatable :: blknynz(:)             !blkny . blknz
   integer, allocatable :: blknxnz(:)             !blknx . blknz
   integer, allocatable :: blknxnynz(:)           !blknx . blkny . blknz
   integer              :: ngrdblkx               !# of grids per block in x
   integer              :: ngrdblky               !# of grids per block in y
   integer              :: ngrdblkz               !# of grids per block in z
   integer              :: xmblk                  !total blocks in x
   integer              :: ymblk                  !total blocks in y
   integer              :: zmblk                  !total blocks in z
   integer              :: ilower                 !
   integer              :: iupper                 !
   integer              :: jlower                 !
   integer              :: jupper                 !
   integer              :: klower                 !
   integer              :: kupper                 !
   integer              :: nfirstwpad(MAXBLOCK+1) !Described in the text
   integer              :: nfirst0pad(MAXBLOCK+1) !Described in the text
   integer              :: nfirst4s3f(MAXBLOCK+1) !Described in the text
   integer, allocatable :: lfocuswpad(:)          !Described in the text
   integer, allocatable :: lfocus0pad(:)          !Described in the text
   integer, allocatable :: lfocus4s3f(:)          !Described in the text
   integer, allocatable :: fineblockindex(:)      !for load balancing
   integer              :: iblkdebug              !
   _REAL_,  allocatable :: blkgox(:)              !block origin in x
   _REAL_,  allocatable :: blkgoy(:)              !block origin in y
   _REAL_,  allocatable :: blkgoz(:)              !block origin in z
   _REAL_,  allocatable :: coarsephi(:)           !saved 1st level solution
!  _REAL_, allocatable  :: spv(:)                 !Ion placement
   integer              :: saltout                !Ion placement
   logical              :: outsalt                !Ion placement
   _REAL_               :: stern                  !Ion Placement

contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver of PBMD energy and forces
!     call pb_force( natom,nres,ntypes,npdec,ix(i02),ix(i04),ix(i06),&
!                    ix(i10),cn1,cn2,xx(l15),x,f,evdw,eelt,epol)
subroutine pb_force( natom,nres,ntypes,npdec,ipres,iac,ico,natex,cn1,cn2,cg,x,f,enb,eel,eelrf )
    
   use solvent_accessibility, only : dprob, radi, radip, radip2, radip3, nzratm, &
#if defined SANDER || defined LIBPBSA
                                     sa_init, sa_driver, sa_free, sa_free_mb
#else
#ifdef MPI
                                     sa_init, sa_driver, sa_free, sa_free_mb,  &
                                     saslave_init
#else
                                     sa_init, sa_driver, sa_free, sa_free_mb
#endif /* MPI */
#endif /*SANDER or LIBPBSA*/
   use decomp, only : irespw, jgroup
   use pbtimer_module

   ! Common variables
    
#  include "pb_md.h"
#ifdef SANDER
#  undef _REAL_
#  include "../../../src/sander/md.h"
#  include "../../../src/sander/box.h"
   integer, parameter :: mytaskid = 0
   integer, parameter :: numtasks = 1
#elif defined LIBPBSA
#  include "md.h"
#  include "box.h"
   integer, parameter :: mytaskid = 0
   integer, parameter :: numtasks = 1
#else
#  include "md.h"
#  include "box.h"
#ifdef MPI
  include "mpif.h"
#  include "parallel.h"
#else  /*MPI*/
   integer, parameter :: mytaskid = 0
   integer, parameter :: numtasks = 1
#endif /*MPI*/
#endif /* SANDER */
#  include "extra.h"
    
   ! Passed variables
    
   integer natom, nres, ntypes, npdec, ipres(*), iac(*), ico(*), natex(*)
   _REAL_ cn1(*), cn2(*), cg(natom), x(3,natom), f(3,natom)
   _REAL_ enb, eel, eelrf
 
   ! Local variables

   integer i, j, k
   integer iatm, proatm, atmfirst, atmlast
   integer atmind(natom)
   _REAL_ acg(natom)
   _REAL_ pbcutcap, pbxcap, pbycap, pbzcap
   _REAL_ eelrffd, eelrfmp
   _REAL_ pbfrc(3,natom)
   _REAL_ ionene

   ! Local multiblock variables

   integer ipermute     !
   integer iblock       !
!  integer boolsorted
   integer mynblock     !
   integer lvl2begin    !
   integer orphanblkn   !
   integer guess_int    !
!#endif
   integer ierr         !
   integer taskpiece    !
   integer myblkstart   !
   integer myblkend     !
   integer tasknatom    !
   integer ihavedone    !
   ! readin info verification
   integer tmp_nsatm    !
   integer mingrdblk    !
   integer myldim       !
   logical indexmatched !

   !This is not an efficient way to use the memory, need to estimate
   !the size of the array.
   integer, allocatable :: tmpindex(:) !

   _REAL_ myh           !
   _REAL_ guess_float   !
   _REAL_ blk_eelrffd   ! temporary storage for eelrffd
   _REAL_ blk_eel       ! temporary storage for eel
   _REAL_ blk_enb       ! temporary storage for enb
#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
   _REAL_ , allocatable :: recvbuf1(:), recvbuf2(:), recvbuf3(:)
#endif /*MPI*/
#endif /*ndef SANDER or LIBPBSA*/

   ! end of multi-block

   logical localpbgrid

   ! Variables initialization
 
   enb = ZERO; eel = ZERO; eelrf = ZERO
   eelrffd = ZERO; eelrfmp = ZERO
   pbfrc = ZERO; ionene = ZERO
   atmind = 0

   ! Variables initialization, multi-block

   myblkstart = 0
   myblkend   = 0
   ierr = 0
   acg = ZERO !making sure clien(s) have clean acg

   firstleveldone = .false.

#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
   call pbslave_init(natom)
   call saslave_init(natom)
#endif /* MPI */
#endif /*ndef SANDER or LIBPBSA*/

   !End of multi-block initializatioin

   if ( ifcap /= 0 .and. (ifcap < 3 .or. ifcap > 5) ) then
      pbcutcap = cutcap+TWO; pbxcap = xcap; pbycap = ycap; pbzcap = zcap 
      radi(0) = pbcutcap; acrd(1,0) = pbxcap; acrd(2,0) = pbycap; acrd(3,0) = pbzcap
      radip3(1) = radi(0); nzratm(1) = 0
   else
      pbcutcap = ZERO; pbxcap = ZERO; pbycap = ZERO; pbzcap = ZERO
   end if

   ! split atoms into internal/external and update nblist
 
   call pbtimer_start(PBTIME_PBLIST)
   if ( pbgrid ) then
      if( mpopt == 1 ) then
         !multipole expansion (commented by mj)
         call pb_atmpart(pbverbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                         inatm,outwat,oution,ipres,outflag, &
                         pbxcap,pbycap,pbzcap,pbcutcap,sepbuf,x,ifcap)
      else if ( ifcap == 2 ) then 
         ! Use cutcap here, not pbcutcap because the latter is augmented by TWO
         call pb_atmpart(pbverbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                         inatm,outwat,oution,ipres,outflag, &
                         pbxcap,pbycap,pbzcap,cutcap,0.0d0,x,ifcap)
      else if ( ifcap == 5 ) then 
         ! Use cutcap here, not pbcutcap because the latter is augmented by TWO
         call pb_atmpart2(pbverbose,pbprint,natom,ibgwat,ienwat,ibgion,ienion, &
                          inatm,outwat,oution,ipres,outflag, &
                          cutcap,x)
      else if ( ligand ) then
         ! This option will be visited if it is not the first time pb_force got
         ! called, however there would be something to do here once we've done
         ! MD.
         continue
      else
         ! Multiblock will go here as well as other conditions
         outflag = 0
      end if
   else 
      if ( mpopt == 1 .or. ifcap == 2 .or. ifcap == 5 ) outflag = outflagorig
   end if

   call pb_atmconv(mpopt,ifcap,natom,ibgwat,ienwat,ibgion,ienion,atmind,ipres,x,cg,acg)

   ! This is for the global run for the coarse grid
   ! If ligand/multiple block is used, these will be updated later in docklist

   if ( ntnba == 1 .and. max(cutnb,cutsa,cutfd) > ZERO ) call pb_atmlist(pbverbose,pbprint,&
      maxnba,natom,ntypes,iac,ico,natex,nshrt,nex,iex,iar1pb,iprshrt,&
      cutnb,cutsa,cutfd,cn1,cn2,cn1pb,cn2pb,cn3pb,cg,acrd(1,1))
   if ( ntnbr == 1 ) ntnbr = 0
   if ( ntnba == 1 ) ntnba = 0
   call pbtimer_stop(PBTIME_PBLIST)
 
   call pbtimer_start(PBTIME_PBSETUP)
   if ( mpopt /=2 .and. pbgrid ) then
      if ( ligand ) &
         call pb_atmpart3(pbverbose,pbprint,natom,buffer,xmin,xmax,ymin,ymax,&
              zmin,zmax,liveflag,realflag,outflag,x)
      if (ifcap == 2 .or. ifcap == 5) then
         call pb_setgrd(ipb,pbverbose,pbprint,pbinit,pbgrid,ifcap,1,inatm,pbxcap,pbycap,pbzcap,pbcutcap)
      else
         call pb_setgrd(ipb,pbverbose,pbprint,pbinit,pbgrid,ifcap,1,natom,pbxcap,pbycap,pbzcap,pbcutcap)
      end if
   end if
   call pbtimer_stop(PBTIME_PBSETUP)

   ! compute grid-independent sas calculations for dielectric assignment
   ! when ifcap /= 0, no need to comptue sas

   call pbtimer_start(PBTIME_PBSAS)
!write(600+mytaskid,*)natom,natom,ifcap,dprob,radi,radip,radip2,outflag;call mexit(0,0)
   if ( srsas .and. ( ifcap == 0 .or. ifcap == 5 ) ) then
      if( ifcap == 5 ) then
         call sa_init(pbverbose,pbprint,natom,inatm,ifcap,dprob,radi,radip,radip2,outflag)
         ! the call here requires verification if we need to take care of multiblock as well
         call sa_driver(pbverbose,pbprint,ipb,inp,natom,inatm,dosas,ndosas,npbstep,nsaslag,&
                        ligand, outflag,&
                        acrd(1,1),iar1pb(1,0),iprshrt,nex,iex,.false.)
      else
         call sa_init(pbverbose,pbprint,natom,natom,ifcap,dprob,radi,radip,radip2,outflag)
!write(600+mytaskid,*)ipb,inp,natom,natom,dosas,ndosas,npbstep,nsaslag,ligand,multiblock,outflag;call mexit(0,0)
         call sa_driver(pbverbose,pbprint,ipb,inp,natom,natom,dosas,ndosas,npbstep,nsaslag,&
                        (ligand .or. multiblock), outflag,&
                        acrd(1,1),iar1pb(1,0),iprshrt,nex,iex,.false.)
      end if
   end if
   call pbtimer_stop(PBTIME_PBSAS)

   call pbtimer_start(PBTIME_PBLIST)
   ! Multi-block:
   ! Following part is for level 2 focusing task partition
   if ( multiblock ) then
      ! Current task partitioning: Predetermined distribution
      ! o o o o    For 4 simultaneous tasks doing 13 blocks, mytaskid
      ! o o o o    0, 1, 2 all get 3 blocks respectively. Mytaskid 3 
      ! o o o o    gets 4 blocks. Note: Queuing the blocks inside MPI
      !       o    is a pain in the neck, so it's not in our future
      !            plan.
      mynblock = levelblock(2)
      taskpiece  = levelblock(2)/numtasks
      orphanblkn = levelblock(2)-taskpiece*numtasks
      ! doing serial
      if ( numtasks < 2 ) then
         myblkstart = 1
         myblkend   = mynblock
      ! doing parallel
      else if ( orphanblkn == 0 .or. mytaskid < numtasks-orphanblkn ) then
         myblkstart =     mytaskid*taskpiece + 1
         myblkend   = (mytaskid+1)*taskpiece
      else if ( mytaskid == numtasks-orphanblkn ) then
         ! when orphanblkn == 1, it adds one block to the last task
         myblkstart =     mytaskid*taskpiece+1
         myblkend   = (mytaskid+1)*taskpiece+1
      else if ( mytaskid > numtasks-orphanblkn ) then
         myblkstart =     mytaskid*taskpiece+1+mytaskid+orphanblkn-numtasks
         myblkend   = (mytaskid+1)*taskpiece+1+mytaskid+orphanblkn-numtasks
      else
         write(6,*) "exception caught, block division error for MPI"
         call mexit(6,1)
      end if

      if ( master ) then
         do i = 1, numtasks
            if ( orphanblkn == 0 .or. i < numtasks-orphanblkn ) then
               write(6,*) "thread:",i,"is for",taskpiece  ,"fine blocks."
            else if ( i == numtasks-orphanblkn ) then
               write(6,*) "thread:",i,"is for",taskpiece+1,"fine blocks."
            else if ( i > numtasks-orphanblkn ) then
               write(6,*) "thread:",i,"is for",taskpiece+1,"fine blocks."
            endif
         enddo
      endif
      taskpiece = myblkend - myblkstart + 1

      myh = savh(2)

      mingrdblk=min(ngrdblkx,ngrdblky,ngrdblkz)
      if ( buffer/myh > mingrdblk ) then
         write(6,*) "Big padding",int(buffer/myh), &
                    "is unrealistic (the least grdblk is", &
                    mingrdblk,"), please give up."
      endif
      if (mingrdblk == 0) then
         write(*,*) "grdblkx or y, z contains zero."; stop
      endif
      ! Possible problem: mjhsieh forgot why there is a "+ 2"
      myldim=ceiling(2 + (buffer*2/myh/(mingrdblk-1)+3)**3) * natom
      myldim=min(myldim,mynblock*natom)
      ! Multithread situation is not considered so far, thus no benifit
      ! from multithread to reducing array sizes "for now".
      allocate(lfocuswpad(myldim),    stat=ierr); if (ierr /= 0) call mexit(6,1)
      allocate(lfocus0pad(2*natom+1), stat=ierr); if (ierr /= 0) call mexit(6,1)
      allocate(lfocus4s3f(myldim),    stat=ierr); if (ierr /= 0) call mexit(6,1)

      lfocus4s3f = 0 !initialization
      lfocuswpad = 0 !initialization
      lfocus0pad = 0 !initialization
      nfirst4s3f = 0 !initialization
      nfirstwpad = 0 !initialization
      nfirst0pad = 0 !initialization

      ! Following part initializes for level one
      nfirst4s3f(1) = 1
      nfirstwpad(1) = 1
      nfirst0pad(1) = 1

      ! Grid Partition for Level 2
      ! lfocuswpad, lfocus0pad, nfirstwpad and nfirst0pad got updated
      call pb_atmpart_mb(pbverbose, pbprint, natom, 1, mynblock, myh,&
           blknx,  blkny,  blknz,   blkxlo, blkylo, blkzlo,          &
           blkxup, blkyup, blkzup,  blkgox, blkgoy, blkgoz,          &
           x, lfocuswpad, lfocus0pad, lfocus4s3f,                    &
              nfirstwpad, nfirst0pad, nfirst4s3f, master)
      ! So the information needed for the solver are:
      !    myblkstart: the first block for focus finegrid in this thread
      !    myblkend:   the last block for focus finegrid in this thread
      !    lfocuswpad: serial number list of the atoms covered by the
      !                focus fine grid with the padding zone
      !    lfocus0pad: serial number list of the atoms covered by the
      !                focus fine grid without the padding zone
      !    lfocus4s3f: atom list for doing surface for blocks
      !    nfirstwpad: the first position in lfocuswpad that belongs to
      !                current block (level 2 blocks start from 2)
      !    nfirst0pad: the first position in lfocus0pad that belongs to
      !                current block (level 2 blocks start from 2)
      !    nfirst4s3f: the first position in lfocus4s3f
      tasknatom = nfirst0pad(mynblock+1)-nfirst0pad(1)
      if ( tasknatom .ne. natom ) then
         write(6, *) 'pb_force(): Atom partition error',tasknatom,'/=',natom
         call mexit(6, 1)
      endif
   end if
   call pbtimer_stop(PBTIME_PBLIST)

   ! for focussing run, liveflag, outflag, realflag are updated
   ! atom list is updated next
   ! surface area is then updated

   if ( ligand ) then
      call pb_atmlist(pbverbose,pbprint,maxnba,natom,ntypes,iac,ico,natex, &
              nshrt,nex,iex,iar1pb,iprshrt,cutnb,cutsa,cutfd,cn1,cn2,cn1pb, &
              cn2pb,cn3pb,cg,acrd)
      call sa_driver(pbverbose,pbprint,ipb,inp,natom,natom,dosas,ndosas,    &
              npbstep,nsaslag,ligand,outflag,acrd(1,1),iar1pb(1,0),iprshrt, &
              nex,iex,.false.)

   ! BEGIN OF MULTIPLE BLOCK LOOP: LOOPING OVER BLOCKS

   else if ( multiblock ) then

      ! fineblockindex: storing the order of the multiblock for
      !                 load balancing or just sequential order

      allocate( fineblockindex(levelblock(2)), stat = ierr )
      do i = 1, levelblock(2)
         fineblockindex(i)=i
      enddo

      ! This will shuffle the order, which is a fake load balancing.

      call amrset(8249562)
      if ( master ) then

         allocate(       tmpindex(levelblock(2)), stat = ierr )
         do i = 1, levelblock(2)
            call amrand(guess_float) ! from 0 to 1
            guess_float = guess_float * (levelblock(2)+1-i)
            guess_int   = ceiling(guess_float)
            ipermute    = fineblockindex(guess_int)
            fineblockindex(guess_int) = fineblockindex(levelblock(2)+1-i)
            fineblockindex(levelblock(2)+1-i) = ipermute
         enddo
         tmpindex=fineblockindex

         ! Sort the list according to block gridsize (disabled)
         ! To enable it, you need to verify the index with the new scheme.
!        do while (boolsorted == 0)
!           boolsorted = 1
!           do i = 2, levelblock(2)
!              if (blknxnynz(tmpindex(i-1)) < blknxnynz(tmpindex(i))) then
!                 ipermute     =tmpindex(i-1)
!                 tmpindex(i-1)=tmpindex(i  )
!                 tmpindex(i  )=ipermute
!                 boolsorted = 0
!              endif
!           enddo
!        enddo

         !distributing blocks among available nodes

         k = 1
         do i = 0, numtasks-1
            do j = 1, levelblock(2)
               if ( mod(j,numtasks) == i ) then
                  fineblockindex(k)=tmpindex(j)
                  k = k + 1
               endif
            enddo
         enddo
         deallocate( tmpindex, stat = ierr )
         ihavedone = 0
      endif

#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
      ! broadcast the ownership of blocks to each thread
      call MPI_BCAST(fineblockindex(1),levelblock(2),MPI_INTEGER,0,CommSANDER,ierr)
      REQUIRE(ierr==0)
      !call MPI_BARRIER( CommSANDER, ierr );REQUIRE(ierr==0)
#endif /*MPI*/
#endif /*ndef SANDER or LIBPBSA*/

!#ifdef MPI
!      call mpi_bcast(fineblockindex,mynblock,MPI_INTEGER,0,commsander,ierr)
!#endif

      ! add FD reaction field energy/force

      call pbtimer_start(PBTIME_PBFDFRC)
      blk_eelrffd = ZERO
      blk_eel     = ZERO
      blk_enb     = ZERO
      allocate(coarsephi(savxmymzm(1)), stat = ierr)
      blkloop: do iblock=1, levelblock(2)
        iblkdebug = iblock
        indexmatched = .false.
        matchingindex: do i = myblkstart, myblkend
           if ( iblock == fineblockindex(i) ) then
              indexmatched = .true.
              cycle matchingindex
           end if
        end do matchingindex
        if ( .not. indexmatched ) cycle blkloop
        liveflag = 0 ! 1 if inside the block  w/o pad
        realflag = 0 ! 1 if inside the block with pad
        outflag  = 1 ! 1 if outside the geometry for surface
        do i = nfirst0pad(iblock), nfirst0pad(iblock+1)-1
           liveflag( lfocus0pad(i) ) = 1
        end do
        do i = nfirstwpad(iblock), nfirstwpad(iblock+1)-1
           realflag( lfocuswpad(i) ) = 1
        end do
        do i = nfirst4s3f(iblock), nfirst4s3f(iblock+1)-1
           outflag ( lfocus4s3f(i) ) = 0
        end do
        savxm(2)=blknx(iblock)
        savym(2)=blkny(iblock)
        savzm(2)=blknz(iblock)
        savxmym(2)=blknxny(iblock)
        savxmymzm(2)=blknxnynz(iblock)
        savgox(2)=blkgox(iblock)
        savgoy(2)=blkgoy(iblock)
        savgoz(2)=blkgoz(iblock)

!write(600+mytaskid,*)maxnba,natom,ntypes,iac(1:natom),ico(1:ntypes*ntypes),natex(1:392);call mexit(0,0)
!write(600+mytaskid,*)nshrt(0:natom),nex(natom),iex(1:64,1:natom),iar1pb(1:6,0:natom);call mexit(0,0)
!write(600+mytaskid,*)iprshrt(maxnba),cutnb,cutsa,cutfd,cn1(1:1830),cn2(1:1830),cn1pb(1:maxnba);call mexit(0,0)
!write(600+mytaskid,*)cn2pb(1:maxnba),cn3pb(1:maxnba),cg,acrd(1:3,1:natom);call mexit(0,0)
        call pb_atmlist(pbverbose,pbprint,maxnba,natom,ntypes,iac,ico,natex, &
                nshrt,nex,iex,iar1pb,iprshrt,cutnb,cutsa,cutfd,cn1,cn2,cn1pb, &
                cn2pb,cn3pb,cg,acrd)
!write(600+mytaskid,*)cn3pb(1:maxnba),cg,acrd(1:3,1:natom);call mexit(0,0)
!write(600+mytaskid,*)cn1(1:1830),cn2(1:1830),cn1pb(1:maxnba),cn2pb(1:maxnba);call mexit(0,0)
!write(600+mytaskid,*)iar1pb(1:6,0:natom),iprshrt(maxnba),cutnb,cutsa,cutfd;call mexit(0,0)
!write(600+mytaskid,*)ntypes,iac(1:natom),ico(1:ntypes*ntypes),natex(1:392),nshrt(0:natom),nex(natom),iex(1:64,1:natom);call mexit(0,0)
!write(600+mytaskid,*)ntypes,iac(1:natom),ico(1:ntypes*ntypes),natex(1:392),nshrt(0:natom),nex(natom),iex(1:64,1:natom);call mexit(0,0)

        ! In this implementation of multiblock focusing, each block is treated
        ! as exactly like a ligand box, so we hijack the ligand flag for now.
        ! Here is the interface before the hijack:
        ! call sa_driver(pbverbose,pbprint,ipb,inp,natom,natom,dosas,ndosas,  &
        !       npbstep,nsaslag,ligand,outflag,acrd(1,1),iar1pb(1,0),iprshrt, &
        !       nex,iex)

        call sa_driver(pbverbose,pbprint,ipb,inp,natom,natom,dosas,ndosas,    &
                npbstep,nsaslag,.true.,outflag,acrd(1,1),iar1pb(1,0),iprshrt, &
                nex,iex,.false.)
        !write(6,*) '-debug: pb_force: done with sa_driver';flush(6)
!write(600+mytaskid,*)dosas,ndosas,npbstep,nsaslag,outflag(1:natom),acrd(1:3,1:natom);call mexit(0,0)
!write(600+mytaskid,*)iar1pb(1:6,0:natom),iprshrt(1:maxnba),nex(1:natom),iex(1:64,1:natom);call mexit(0,0)
        ! End of sa_driver hijack.

        call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,natom,    &
                npdec,idecomp,irespw,ipres,jgroup,ibgwat,ibgion,pbfrc,        &
                eelrffd,ionene,npbstep,npbgrid,nstlim)
        blk_eelrffd = blk_eelrffd + eelrffd
        ihavedone = ihavedone + 1
        if ( srsas .and. ihavedone < taskpiece ) then
           call sa_free_mb( dosas,ndosas )
        end if

        ! cutnb > 0 is required in multiblock

        !if ( cutnb == ZERO ) then
        !  call pb_directnocut(natom,proatm,inatm,ipres,ibgwat,ienwat,ibgion,ienion,ntypes,eneopt,idecomp,ifcap, &
        !                      iac,ico,nex,iex,cn1,cn2,acg,acrd(1,1),pbfrc,eel,enb)
        !else
           if( ifcap == 2 .or. ifcap == 5) then
              call pb_directwtcut(natom,inatm,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1), &
                                  pbfrc,eel,enb)
           else
              call pb_directwtcut(natom,natom,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1), &
                                  pbfrc,eel,enb)
           end if 
         !endif
         blk_eel = blk_eel + eel
         blk_enb = blk_enb + enb
      end do blkloop
      deallocate (coarsephi, stat = ierr)
      outflag = 0
      realflag = 1
      eelrffd = blk_eelrffd
      eel = blk_eel
      enb = blk_enb
      ! PHI MAP SHOULD BE EXPORTED IN FDFRC, NOT HERE.
      ! To do so xs is required for assembling the phi
      call pbtimer_stop(PBTIME_PBFDFRC)
#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
      allocate(recvbuf1(numtasks), stat = ierr); recvbuf1 = ZERO
      allocate(recvbuf2(numtasks), stat = ierr); recvbuf2 = ZERO
      allocate(recvbuf3(numtasks), stat = ierr); recvbuf3 = ZERO
      call MPI_GATHER(eelrffd,  1, MPI_DOUBLE_PRECISION, recvbuf1, 1, &
         MPI_DOUBLE_PRECISION, 0, CommSANDER, ierr)
      call MPI_GATHER(eel    ,  1, MPI_DOUBLE_PRECISION, recvbuf2, 1, &
         MPI_DOUBLE_PRECISION, 0, CommSANDER, ierr)
      call MPI_GATHER(enb    ,  1, MPI_DOUBLE_PRECISION, recvbuf3, 1, &
         MPI_DOUBLE_PRECISION, 0, CommSANDER, ierr)
      if (master) then
         eelrffd = SUM(recvbuf1)
         eel     = SUM(recvbuf2)
         enb     = SUM(recvbuf3)
      endif
#endif /* MPI */
#endif /*ndef SANDER or LIBPBSA*/
   end if
   ! HERE GOES THE ENERGY SUMMATION (parallel)
   ! END OF MULTIPLE BLOCK LOOP

   ! add FD reaction field energy/force

   call pbtimer_start(PBTIME_PBFDFRC)
   if ( multiblock ) then
      ! For multiblock focusing, FD force is already done in the previous
      ! section. But that's not the case for ligand.
      continue
   else if ( epsout /= epsin .and. mpopt /= 2 ) then
      ! In the case of ifcap == 2,5, only map crg within cap to grid (atmlast == inatm),
      ! else map all crg (atmlast == natom)
      if( ifcap == 2 .or. ifcap == 5 ) then
         call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,inatm,npdec,idecomp,irespw, &
                       ipres,jgroup,ibgwat,ibgion,pbfrc,eelrffd,ionene,npbstep,npbgrid,nstlim)
      else
         call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,natom,npdec,idecomp,irespw, &
                       ipres,jgroup,ibgwat,ibgion,pbfrc,eelrffd,ionene,npbstep,npbgrid,nstlim)
      end if
   else if (epsout == epsin) then
        call pb_fdfrc(pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,natom,npdec,idecomp,irespw, &
                       ipres,jgroup,ibgwat,ibgion,pbfrc,eelrffd,ionene,npbstep,npbgrid,nstlim) 
  ! else
  !    write(6,*) "Exception caught in pb_force."; call mexit(6,1)
   end if
   call pbtimer_stop(PBTIME_PBFDFRC)

   ! clean up for sas calculations

   call pbtimer_start(PBTIME_PBSETUP)
   if ( srsas .and. (ifcap == 0 .or. ifcap == 5) ) then
      call sa_free( dosas,ndosas,.false. )
   end if
   call pbtimer_stop(PBTIME_PBSETUP)

   if ( saopt < 0 ) return

   ! add MP reaction field energy/forces when ifcap /= 0
 
   call pbtimer_start(PBTIME_PBMP)
   if ( mpopt /= 0 .and. epsout /= epsin ) then
      if ( mpopt == 1 ) then      ! multipole expansion for boundary atoms
         atmfirst = inatm + 1
         atmlast  = natom
      else if ( mpopt == 2 ) then ! multipole expansion for all atoms
         atmfirst = 1
         atmlast  = natom
      end if
      call pb_mpfrc(natom,atmfirst,atmlast,lmax,pbcutcap,pbxcap,pbycap,pbzcap,&
              epsin,epsout,acrg,acrd(1,1),pbfrc,eelrfmp)
   end if
   call pbtimer_stop(PBTIME_PBMP)
 
   ! add direct coulombic and nonbonded forces
    
   call pbtimer_start(PBTIME_PBDIRECT)
   if ( multiblock ) then
      ! For multiblock focusing, pb_direct SHOULD BE DONE inside the loop
      continue
   else if ( cutnb == ZERO ) then
      call pb_directnocut(natom,proatm,inatm,ipres,ibgwat,ienwat,ibgion,ienion,ntypes,eneopt,idecomp,ifcap, &
                          iac,ico,nex,iex,cn1,cn2,acg,acrd(1,1),pbfrc,eel,enb)
   else
      if( ifcap == 2 .or. ifcap == 5) then
         call pb_directwtcut(natom,inatm,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1), &
                             pbfrc,eel,enb)
      else
         call pb_directwtcut(natom,natom,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,acrd(1,1), &
                             pbfrc,eel,enb)
      end if 
   end if
   call pbtimer_stop(PBTIME_PBDIRECT)

   ! returning:
   ! adding the nonbonded forces to the MD forces

   eel = eel * eps0/epsin
   if ( eneopt == 1 .and. (bcopt < 6 .or. bcopt > 9) ) then 
      eel = eel + eelrffd + eelrfmp
      eelrf = ZERO
   else
      eelrf = eelrffd + eelrfmp
   end if

   if( ifcap == 2 .or. ifcap == 5 ) then
     atmlast = inatm
   else
     atmlast = natom
   endif

   do iatm = 1, natom
      f(1,iatm) = f(1,iatm) + pbfrc(1,mapout(iatm))
      f(2,iatm) = f(2,iatm) + pbfrc(2,mapout(iatm))
      f(3,iatm) = f(3,iatm) + pbfrc(3,mapout(iatm))
   end do

!  open (unit = 120, file = 'total.dat')

!  write(120,*) ' :::: Atomic forces ::::'
!  do iatm = 1, natom
!     write(120,'(3e20.6)') f(1,iatm),f(2,iatm),f(3,iatm)
!  end do

!  open (unit = 121, file = 'total2.dat')

!  write(121,*) ' :::: Atomic forces ::::'
!  do iatm = 1, natom
!     write(121,'(3e20.6)') pbfrc(1,iatm),pbfrc(2,iatm),pbfrc(3,iatm)
!  end do

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ convert passed coordinates and charges to the internal format
subroutine pb_atmconv( mpopt,ifcap,natom,ibgwat,ienwat,ibgion,ienion,atmind,ipres,x,cg,acg )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   ! Passed variables
   
   integer mpopt, ifcap, natom, ibgwat, ienwat, ibgion, ienion
   integer atmind(natom), ipres(*) 
   _REAL_ x(3,natom), cg(natom), acg(natom)
    
   ! Local variables
    
   integer i, j, ifirst, ilast, iatm, ires, num
    
   if ( mpopt == 1 .or. ifcap == 2 .or. ifcap == 5 ) then
       
      ! copy reordered coord/charge to private arrays for pb/mp or ifcap == 2,5
      ! protein atoms go into internal portion (so do IONS!)
       
      ifirst = 1; ilast = ipres(ibgwat)-1

      do iatm = ifirst, ilast
         acrd(1,iatm) = x(1,iatm); acrd(2,iatm) = x(2,iatm); acrd(3,iatm) = x(3,iatm)
         acrg(iatm) = cg(iatm)/18.2223d0; acg(iatm) = cg(iatm); atmind(iatm) = iatm
         mapout(iatm) = iatm
      end do
       
      ! water atoms go into internal/external portion
       
      ifirst = ipres(ibgwat); ilast = natom

      i = ifirst; j =  inatm + 1
      do iatm = ifirst, ilast
         if ( outflag(iatm) == 0 ) then
            acrd(1,i   ) = x(1,iatm); acrd(2,i   ) = x(2,iatm); acrd(3,i   ) = x(3,iatm)
            acrg(i   ) = cg(iatm)/18.2223d0; acg(i   ) = cg(iatm); atmind(i   ) = iatm
            mapout(iatm) = i
            i = i + 1
         else
            acrd(1,j   ) = x(1,iatm); acrd(2,j   ) = x(2,iatm); acrd(3,j   ) = x(3,iatm);
            acrg(j   ) = cg(iatm)/18.2223d0; acg(j   ) = cg(iatm); atmind(j   ) = iatm
            mapout(iatm) = j
            j = j + 1
         end if
      end do

      ! store original outflag array and prepare an updated one for water atoms

      outflagorig = outflag
      do iatm = ifirst, ilast
         if( iatm <= inatm ) then
            outflag(iatm) = 0
         else
            outflag(iatm) = 1
         end if
      end do
       
   else
       
      do iatm = 1, natom
         acrd(1,iatm) = x(1,iatm); acrd(2,iatm) = x(2,iatm); acrd(3,iatm) = x(3,iatm)
         acrg(iatm) = cg(iatm)/18.2223d0; acg(iatm) = cg(iatm)
         mapout(iatm) = iatm
      end do
   end if
 
end subroutine pb_atmconv

end subroutine pb_force
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ driver for FDPB forces and energy
subroutine pb_fdfrc( pbverbose,pbprint,pbgrid,ifcap,ipb,imin,natom,atmlast,npdec,idecomp,irespw, &
                     ipres,jgroup,ibgwat,ibgion,pbfrc,eelrf,ionene,npbstep,npbgrid,nstlim )

   use solvent_accessibility, only : dprob,iprob,radi,radip3,nzratm, &
       narcdot,maxarc,marc,m2narc,fstarc,arcatm,arccrd,savarc,dotarc
!#ifndef SANDER
   use pbtimer_module
!#endif

#  include "flocntrl.h"
#include "parallel.h"

   ! passed variables
    
   logical pbverbose, pbprint, pbgrid  
   integer ifcap, ipb, imin, natom, atmlast, npdec, idecomp, ibgwat, ibgion
   integer npbstep, npbgrid, nstlim
   integer irespw(*), ipres(*), jgroup(*)
   _REAL_ ionene, eelrf, pbfrc(3,natom)!, fnet(3)
    
   ! local variables
    
   integer atmfirst
   integer iatm, lastatm, mpdec, i
   integer cnt, j, k ! phiform = 2
   _REAL_ eelself, eelcoul, rh, fcrd(3,atmlast)
   _REAL_ aa, bb, cc, aa1, bb1, cc1
   _REAL_ bb1cc1, bb_cc1, bb1cc, bb_cc
   _REAL_ acrgtmp, factor
   _REAL_ grdreac, drcreac 
 
   !XP: test variables

   _REAL_  px0,py0,pz0,pu,pex,pey,pez,pcirreg(1:137440,1:6)
   integer ipp,iip
   character(len=20):: filename = "pcrd"
   logical alive
   integer :: status = 0

   character(len=12) phifilename
   character(len=23) phidataname
   integer           phifilenum

   if (phiform == 0 .OR. phiform == 1) then
      phifilename="pbsa_phi.phi"
   else
      phifilename="pbsa_phi.dx "
   end if
   phidataname="Electrostatic Potential"
   phifilenum=64

   if ( do_pbfd == 0 ) return

   mpdec = 1
   if ( idecomp > 2 ) mpdec = npdec

   !-- PB pairwise decomp <<1,mpdec>>
 
   atmfirst = 1
   ! atmfirst is supposedly to be passed

   do m = 1, mpdec
    
      ! do fdpb calculations upto nfocus

      eelself = ZERO; eelcoul = ZERO
      if ( .not. firstleveldone ) xsoffset = 1
      if ( ligand ) then
         icrd = 0
         phi = ZERO; bv = ZERO; xs = ZERO!; fcrd = ZERO; chgrd = ZERO
      end if
      !xs shouldn't be reset in MD, to be fixed

      do level = 1, nfocus
       
         !  re-feed the phi with first level phi

         if ( firstleveldone .and. level == 1 .and. multiblock ) then
            phi(1:savxmymzm(1))=coarsephi(1:savxmymzm(1))
            cycle
         end if

         call pbtimer_start(PBTIME_PBBUILDSYS)

         ! retrieving saved grid data into working variables

         bcopt = savbcopt(level)
         xm = savxm(level); ym = savym(level); zm = savzm(level)
         xmym = xm*ym; xmymzm = xmym*zm
         h = savh(level)
         gox = savgox(level); goy = savgoy(level); goz = savgoz(level)
       
         ! grid-unit version of the atom data array at the working level
       
         rh = ONE/h
         if ( ifcap == 1 ) then
            gcrd(1,0) = (acrd(1,0) - gox)*rh
            gcrd(2,0) = (acrd(2,0) - goy)*rh
            gcrd(3,0) = (acrd(3,0) - goz)*rh
         end if
         do iatm = atmfirst, atmlast
            gcrd(1,iatm) = (acrd(1,iatm) - gox)*rh
            gcrd(2,iatm) = (acrd(2,iatm) - goy)*rh
            gcrd(3,iatm) = (acrd(3,iatm) - goz)*rh
         end do
!write(600,*) xm,ym,zm,gox,goy,goz
!write(600,*) gcrd(1:3,1:natom); close(600)
         !icrd = 0; fcrd = 0d0 ! shouldn't be necessary
         do iatm = atmfirst, atmlast
            if ( ligand .or. multiblock ) then
               if ( level > 1 .and. realflag(iatm) == 0 ) cycle
            endif
            !if ( level > 1 .and. ligand     .and. realflag(iatm) == 0 ) cycle
            !if ( level > 1 .and. multiblock .and. realflag(iatm) == 0 ) cycle
            icrd(1,iatm) = floor(gcrd(1,iatm))
            icrd(2,iatm) = floor(gcrd(2,iatm))
            icrd(3,iatm) = floor(gcrd(3,iatm))
            fcrd(1,iatm) = dble(icrd(1,iatm))
            fcrd(2,iatm) = dble(icrd(2,iatm))
            fcrd(3,iatm) = dble(icrd(3,iatm))
            if ( icrd(1,iatm) > xm-1 .or. icrd(2,iatm) > ym-1 .or. icrd(3,iatm) > zm-1 ) then
               write(6,*) acrd(1:3,iatm)
               write(6,*) "pb_fdfrc(): Atom out of focusing box",icrd(1:3,iatm)
               call mexit(6,1)
            end if
         end do
         gcrg=0d0!shouldn't be necessary
         do iatm = atmfirst, atmlast
            if ( ligand .or. multiblock ) then
               if ( level > 1 .and. realflag(iatm) == 0 ) cycle
            endif
            !if ( level > 1 .and. ligand     .and. realflag(iatm) == 0 ) cycle
            !if ( level > 1 .and. multiblock .and. realflag(iatm) == 0 ) cycle
            aa = gcrd(1,iatm) - fcrd(1,iatm)
            bb = gcrd(2,iatm) - fcrd(2,iatm)
            cc = gcrd(3,iatm) - fcrd(3,iatm)
            bb1 = ONE - bb; cc1 = ONE - cc
            !-- PB decomp
            if (idecomp < 3 ) then
               acrgtmp = acrg(iatm)
            else if (iatm >= ipres(irespw(m)) .and. iatm < ipres(irespw(m)+1)) then
               acrgtmp = acrg(iatm)
            else
               acrgtmp = ZERO
            end if
            aa  = acrgtmp*aa; aa1 = acrgtmp - aa
            bb1cc1 = bb1*cc1; bb_cc1 = bb *cc1
            bb1cc  = bb1*cc ; bb_cc  = bb *cc
            if ( (ifcap == 2 .or. ifcap == 5) .and. outflag(iatm) == 1 ) then
               gcrg(1,iatm) = ZERO; gcrg(2,iatm) = ZERO
               gcrg(3,iatm) = ZERO; gcrg(4,iatm) = ZERO
               gcrg(5,iatm) = ZERO; gcrg(6,iatm) = ZERO
               gcrg(7,iatm) = ZERO; gcrg(8,iatm) = ZERO
            else
               gcrg(1,iatm) = aa1*bb1cc1; gcrg(2,iatm) = aa *bb1cc1
               gcrg(3,iatm) = aa1*bb_cc1; gcrg(4,iatm) = aa *bb_cc1
               gcrg(5,iatm) = aa1*bb1cc ; gcrg(6,iatm) = aa *bb1cc
               gcrg(7,iatm) = aa1*bb_cc ; gcrg(8,iatm) = aa *bb_cc
            end if
         end do
!write(600+mytaskid,*)gcrg(1:8,1:natom);call mexit(0,0)
!write(600+mytaskid,*)ligand,firstleveldone,istrng,iprob,h,gox,goy,goz,xm,ym,xmymzm,saltgrd;call mexit(0,0)

         ! set up dielectric map
         ! when ifcap == 0, do dielectric map assignment everystep
         ! when ifcap /= 0, do dielectric map assignment once only when grid is set up

         if ( ifcap /= 1 ) then

            ! part I,
            ! when solving systems with salt, set up stern layer map, on the grid
            ! points

            if ( istrng /= ZERO ) call pb_ionmap( pbverbose,ifcap,natom,&
                 iprob,h,gox,goy,goz,xm,ym,zm,xmymzm,&
                 outflag,gcrd,radi,atmsas,insas,zv,saltgrd)
!           if(level>1)then
!             write(6666,*) atmsas(1:xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8), &
!                           insas(1:xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8), &
!                              zv(1:xmymzm + xm*ym*2 + ym*zm*2 + xm*zm*2 + xm*4 + ym*4 + zm*4 + 8);stop
!write(600+mytaskid,*)insas(1:xmymzm+xm*ym*2+ym*zm*2+xm*zm*2+xm*4+ym*4+zm*4+8),&
!zv(1:xmymzm+xm*ym*2+ym*zm*2+xm*zm*2+xm*4+ym*4+zm*4+8);call mexit(0,0)
!           endif

            ! part II,
            ! here comes the dielectric map on the grid edges: x, y, and z

!write(600+mytaskid,*)pbverbose,ifcap,ipb,natom,smoothopt,dprob,epsin,epsout,h,gox,goy,goz,xm,ym,zm,xmymzm;call mexit(0,0)
!write(600+mytaskid,*)level,nfocus,narcdot,maxarc,nbnd,nbndx,nbndy,nbndz;call mexit(0,0)
            call pb_exmol_ses( pbverbose,ifcap,ipb,savbcopt,saopt,sasopt,natom,&
                 smoothopt,dprob,epsin,epsout,&
                 h,gox,goy,goz,xm,ym,zm,xmymzm,&
                 level,nfocus,&
                 narcdot,maxarc,nbnd,nbndx,nbndy,nbndz,&
                 outflag,gcrd,acrd,radi,radip3,&
                 marc,m2narc,fstarc,arcatm,dotarc,arccrd,savarc,&
                 atmsas,insas,lvlset,zv,epsx,epsy,epsz,&
                 iepsav,iepsavx,iepsavy,iepsavz,fedgex,fedgey,fedgez)
         elseif ( pbgrid ) then
            call pb_exmol_cap( pbverbose,ifcap )
         end if
!print *,gox,gox+xm*h,goy+h,goy+ym*h,goz+h,goz+zm*h
!print *,gox+h,gox+xm*h,goy,goy+ym*h,goz+h,goz+zm*h
!print *,gox+h,gox+xm*h,goy+h,goy+ym*h,goz,goz+zm*h
         
         ! set up grid charges and store in chgrd
          
         chgrd = ZERO
         call pb_crggrd( xm, ym, zm, 1, atmlast, icrd, gcrg, chgrd(1) )
!write(600+mytaskid,*)gcrg(1:8,1:natom),chgrd(1:xmymzm);call mexit(0,0)

         call pbtimer_stop(PBTIME_PBBUILDSYS)

         if ( saopt < 0 ) then
            if ( level == nfocus ) then
               return
            else
               cycle
            end if
         end if

         ! now call fdpb driver

         call pbtimer_start(PBTIME_PBSOLV)
         if ( (ipb == 1 .or. ipb == 2 .or. ipb ==3 .or. (ipb==4 .and. nfocus /= level)) .and. solvopt /= 7 ) then
!print *,mytaskid,level
!write(600+mytaskid,*)npbstep,npbgrid,nstlim,npbopt,solvopt,level,nfocus,bcopt,h,savh,gox,goy,goz;call mexit(0,0)
!write(600+mytaskid,*)savgox,savgoy,savgoz,xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm;call mexit(0,0)
!write(600+mytaskid,*)pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,ngrdcrg,grdcrg,qgrdcrg;call mexit(0,0)
!write(600+mytaskid,*)nbnd,iepsav,insas,epsx,epsy,epsz,chgrd,saltgrd,phi,bv,cphi;call mexit(0,0)
            call pb_fddrv( npbstep,npbgrid,nstlim,npbopt,solvopt,level,nfocus,bcopt,&
                           h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                           xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                           maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
                           pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
                           ngrdcrg,grdcrg,qgrdcrg,&
                           nbnd,iepsav,insas,epsx,epsy,epsz,&
                           chgrd,saltgrd,phi,&
                           bv,cphi,xs(xsoffset)&
                         )
     !       call printphidx(phi(1:xmymzm)*eps0,xm,ym,zm,gox,goy,goz,h)
            if ( level == 1 .and. multiblock ) coarsephi(1:xmymzm)=phi(1:xmymzm)
         else if ( (ipb == 1 .or. ipb == 2 .or. ipb ==3 .or. (ipb==4 .and. nfocus /= level)) .and. &
          solvopt == 7 ) then
#ifndef SANDER
            ! let's do fft test here
            ! atomic coordinates in Angstrom
            ! atomic charges in electron (au)
            ! electrostatic energy in e^2/A-mol, which can be converted to kcal/mol by
            ! multiplying with a constant of AMBER_ELECTROSTATIC2
            ! electrostatic force in e^2/A^2-mol, which can also be converted to
            ! kcal/mol-A by multiplying the same constant of AMBER_ELECTROSTATIC2
            ! call pb_fft(acrd(1,1),acg,eel,pbfrc,savxm(1),savym(1),savzm(1),savh(1))
            ! acrd --> atom coordinates
            ! acrg --> atom charges
            !! chgrd --> charge grid (will be added later to make use of current charge
            !           grid set up routines, will use internal ngp method for now
            ! xm, ym, zm, h, gox,goy,goz--> grid setup parameters
	    ! eel --> output for calculated long range coulombic energy
            ! pbfrc --> output place holder for forces (not implemented yet)
            ! outphi is a flag (true or false) if set to true, pg_fft will write
            ! the calculated potential grid to the file pb_fft.phi
            ! the output potential map file will be in amber format.
            !-- This implementation currently only performs partical mesh ewald 
            !   (pme) method 
	    ! arguments: coordinates, charges, number atoms, x nodes, y nodes, z nodes, gridspacing
	    ! 		  forces (output), kappa, epsilon inside, epsilon outside, phi (output)
     write(6,*) "Calling pb_fftsolv fftw3"
     call pb_fftsolv(chgrd(1),2,xm,ym,zm,phi,h,64+512)
      !     call printphidx(phi(1:xmymzm),xm,ym,zm,gox,goy,goz,h)
            phi = phi/eps0
            !stop
#else
            write(6,*) "PB_FFT is not supported in SANDER";call mexit(6,1)
#endif /*#ifndef SANDER*/
         else if ( ipb == 4 .and. nfocus == level ) then
            if ( ligand .or. multiblock ) then
               write(6,*) "IIM doesn't work with ligandmask."; call mexit(6,1)
            end if
            call pb_iimdrv( npbstep,npbgrid,nstlim,1,natom,npbopt,solvopt,level,nfocus,bcopt,&
                           natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                           xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                           maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
                           pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
                           gcrd,acrg,&
                           nbnd,iepsav,insas,lvlset,&
                           chgrd,saltgrd,phi,&
                           bv,cphi,xs(xsoffset)&
                         )
         else if ( ipb == 5 .and. nfocus == level ) then
            if ( ligand ) then
               write(6,*) "AUG doesn't work with ligandmask."; call mexit(6,1)
            end if
#ifndef SANDER
            call pb_augdrv( npbstep,npbgrid,nstlim,1,natom,npbopt,solvopt,level,nfocus,bcopt,&
                           natom,h,savh,gox,goy,goz,savgox,savgoy,savgoz,&
                           xm,ym,zm,xmym,xmymzm,savxm,savym,savzm,&
                           maxitn,itn,fmiccg,accept,laccept,wsor,lwsor,inorm,norm,&
                           pbkappa,pbkb,pbtemp,ivalence,istrng,eps0,epsin,epsout,ionene,&
                           gcrd,acrg,&
                           nbnd,iepsav,insas,lvlset,&
                           chgrd,saltgrd,phi,&
                           bv,cphi,xs(xsoffset),&
                           pbverbose, &
                           augtoltype, augctf, augtol &
                         )
#else
            write(6,*) "PB_AUGFFT is not supported in SANDER";call mexit(6,1)
#endif
         else
            write(6,*) "PB Bomb: pb_fdfrc: unknown ipb";call mexit(6,1)
         end if
!        if ( level == nfocus ) then
!          if ( ligand ) then 
!             call mapoutviaidx(phi(1),9999,      1,    xm,      1,     ym,     1,    zm)
!          else
!             call mapoutviacrd(phi(1),9999,-13.0d0,11.0d0,-44.0d0,-21.0d0,11.0d0,35.0d0)
!          endif
!          stop
!        else
!          write(9999,*)phi(1:xmymzm)
!        end if

!#ifndef SANDER
         call pbtimer_stop(PBTIME_PBSOLV)
!#endif
    ! Here to begin the task. XP
    ! What we are doing here is to be fed in the projection crds and calculate
    ! the related grid points to calculate the potential and field of the
    ! projection points of interests. FYI, for sphere system, we can compare
    ! them with the analytical results we have in exactu package. XP
    
       !        do k=1,zm; do j=1,ym; do i=1,xm
       !        write(14,*) phi(i,j,k)
       !        end do; end do; end do
!    inquire ( file='irreg',exist=alive)
!    if(.not. alive) then
!    write(6,*) 'file does not exist' 
!    call mexit(6,1)
!    end if
!    open ( unit=010,file='irreg')
!     do ipp=1,2131
!     read( unit = 010,fmt=*,iostat=status) pcirreg(ipp,1:3)
!     end do
!     do  ipp=1,2131
!     px0 = (pcirreg(ipp,1)-gox)*rh
!     py0 = (pcirreg(ipp,2)-goy)*rh
!     pz0 = (pcirreg(ipp,3)-goz)*rh;!write(111,*) pcirreg(ipp,1:3),goz,rh
!    iip= floor(px0)+(floor(py0)-1)*xm+ (floor(pz0)-1)*xmym 
!     write(13,*) pcirreg(ipp,1:3),phi(iip)*eps0*FOURPI
!  end do
!    pcirreg=0.0d0
!    inquire ( file='pcrd',exist=alive)
!    if(.not. alive) then
!    write(6,*) 'file does not exist' 
!    call mexit(6,1)
!    end if
!    open ( unit=011,file='pcrd')
!     do ipp=1,137440
!     read( unit = 011,fmt=*,iostat=status) pcirreg(ipp,1:6)
!     end do
!     do  ipp=1,137440
!     px0 = (pcirreg(ipp,1)-gox)*rh
!     py0 = (pcirreg(ipp,2)-goy)*rh
!     pz0 = (pcirreg(ipp,3)-goz)*rh;!write(111,*) pcirreg(ipp,1:3),goz,rh
!     call gradu(xm,ym,zm,-ONE,27,10,px0,py0,pz0,pu,pex,pey,pez,phi,zv)
!     write(16,*) pu*eps0*FOURPI,(pex*pcirreg(ipp,4)+pey*pcirreg(ipp,5) & 
!                +pez*pcirreg(ipp,6))*eps0*FOURPI
!     end do

         ! if requested, print a summary when the grid is set up

         if ( pbverbose .and. pbprint ) then
            call pb_print( ifcap, ipb, natom )
            write(6, '(a,I6     )') '   Iterations required        :', itn
            write(6, '(a,F21.10 )') '   Norm of the constant vector:', inorm
            write(6, '(a,F21.13 )') '   Norm of the residual vector:', norm
            write(6, '(a,ES24.16)') '   Convergence achieved       :', norm/inorm
            write(6, '()')
         end if  !  pbverbose .and. pbprint
       
         ! if requested, output phi map
       
         if ( outphi .and. level == nfocus ) then
            if ( phiform == 0 ) then
               write(6,*) 'writing potential map in delphi format'
               open(64,file='pbsa.phi',form="unformatted")
               write(64) ' start of phimap    '
               write(64) ' potential', ' ------ AMBER PBSA ------ phimap in kT/e (0.593kcal/mol-e)  '
               write(64) real((frcfac/0.593d0)*phi(1:xmymzm))
               write(64) ' end of phimap  '
               write(64) real(1.0d0/h),real(cxbox(level)),real(cybox(level)),real(czbox(level)),xm
               close(64)
            elseif ( phiform == 1 ) then
               write(6,*) 'writing potential map in amber format'
               open(64,file='pbsa.phi',form="formatted")
               write(64,*) '# the following data is provided:'
               write(64,*) '# h, gox, goy, goz'
               write(64,*) '# xm, ym, zm'
               write(64,*) '# phi(1:xmymzm) in kcal/mol-e'
               write(64,*) '# mapping between (i,j,k) and phi index:'
               write(64,*) '# i + xm * ( j-1 + ym * ( k-1 ) )'
               write(64,*) '# grid coordinates: xg = gox + h*i; '
               write(64,*) '# yg = goy + h*j; zg = goz + h*k'
               write(64,*) h, gox, goy, goz
               write(64,*) xm, ym, zm
               write(64,*) frcfac*phi(1:xmymzm)
               close(64)
            elseif ( phiform == 2 ) then
            !write dx format phi
               call gen_dx_file(xm,ym,zm,h,gox,goy,goz,phi(1:xmymzm),&
                           phifilename,phifilenum,phidataname)
!                write(6,*) 'writing potential map in dx format'
!                open(67,file='pbsa_phi.dx')
            !   write(67,100) xm,ym,zm
            !   write(67,110) gox+h,goy+h,goz+h
            !   write(67,120) sngl(h),0,0
            !   write(67,121) 0,sngl(h),0
            !   write(67,122) 0,0,sngl(h)
            !   write(67,130) xm,ym,zm
            !   write(67,140) xm*ym*zm
            !   cnt=0
            !   do k=1,zm; do j=1,ym;
            !       do i=1,xm
            !           cnt = cnt+1
            !           write(67,150,advance='no') frcfac*phi(&
            !                                      (i-1)*xm*ym+(j-1)*ym+k)
            !           if (cnt >= 3) then
            !               cnt = 0
            !               write(67,*)''
            !           end if
            !       end do
            !   end do; end do
            !   write(67,*) 'attribute "dep" string "positions"'
            !   write(67,*) 'object "Potential" class field'
            !   write(67,*) 'component "positions" value 1'
            !   write(67,*) 'component "connections" value 2'
            !   write(67,*) 'component "data" value 3'
!                close(67)
!            100 format('object 1 class gridpositions counts',3I5)
!            110 format('origin',3F9.4)
!            120 format('delta',F10.7,2I2)
!            121 format('delta',I2,F10.7,I2)
!            122 format('delta',2I2,F10.7)
!            130 format('object 2 class gridconnections counts',3I5)
!            140 format('object 3 class array type double rank 0 items',I10,' data follows')
!            150 format(ES17.10,' ')
            else
               write(6,*) 'Warning, unrecognizable phimap output format.'
               write(6,*) 'No phimap will be written out.'
            endif
         end if  ! outphi .and. level == nfocus

         if ( firstleveldone ) cycle
         xsoffset = xsoffset + savxmymzm(level) + 2*savxmym(level)
         firstleveldone = .true.


      end do  !  level = 1, nfocus

      !  Useful start if you want a working xs
      !if ( firstleveldone ) &
      !   xsoffset = xsoffset + savxmymzm(level) + 2*savxmym(level)

      ! compute fd energy and force by the qE option
      ! note that self forces are zero
      ! delete fd coulomb energy and forces for all close pairs
      ! dbf is computed by Gilson et al


      if ( eneopt == 1 ) then

         ! compute total qE energy and forces and delete self energy, note that self forces are zero

         if ( intopt == 1 ) call pb_qefrc( natom, 1, atmlast, idecomp, eelrf, eelself, outflag, pbfrc, phi )
         zv = -dble(insas) ! pseudo signed distance function
         if ( intopt == 2 ) then
            if ( bcopt < 6 .or. bcopt > 9 ) then
               write(6,*) 'PB Bomb in pb_force(): intop = 2 does not support bcopt < 6'
               call mexit(6, 1)
            end if
            call pb_qefrc2( natom, 1, atmlast, eelrf, outflag, pbfrc, zv(1), phi )
         end if
!write(6666,*) eelrf

         ! delete fd grid energy and forces for all close pairs
         ! when bcopt >= 6, we only have reaction field energy in eelrf
         ! no need to any corrections, coulombic energy should come from
         ! pb_direct ...

         if ( bcopt < 6 .or. bcopt > 9 ) then
            call pb_fdcoulomb( natom, 1, atmlast, idecomp, eelcoul, outflag, pbfrc )
            eelrf = eelrf - eelself - eelcoul
         end if
!write(7777,*) eelrf,eelself+eelcoul

!        write(6,*) 'final eelrf/ionene in kcal/mol', frcfac*eelrf, frcfac*ionene
!        write(6,*) 'final eelself in kcal/mol', frcfac*eelself
!        write(6,*) 'final eelcoul in kcal/mol', frcfac*eelcoul

         ! add ion contributions for nonlinear PBE ...

         if ( npbopt /= 0 ) then
            eelrf = eelrf + ionene
         end if

         eelrf = frcfac*eelrf

         if ( frcopt == 1 ) then

            ! first printing of qE forces, in kcal/mol-Angstrom
         qef=0.0d0 !XP: calculate qE force
            open (unit = 102, file = 'force.dat')
            write(102,*) ' :::: Atomic qE forces ::::'
            do iatm = 1, natom
               write(102,'(3e20.6)') pbfrc(1:3,iatm)*frcfac
               qef(1)=qef(1)+pbfrc(1,iatm)
               qef(2)=qef(2)+pbfrc(2,iatm)
               qef(3)=qef(3)+pbfrc(3,iatm)
            end do
               qef=qef*frcfac
            write(6,*) 'Vector qE forces is',qef(1:3)
            write(6,*) 'Total qE forces is',sqrt(qef(1)**2+ qef(2)**2+ qef(3)**2)
!           pbfrc(1:3,iatm) = 0.d0

            ! Ray: This should be turned on for the printing purpose.
!           pbfrc = ZERO ! resetting to zero for printing only
            call pb_dbfrc_fld(pbverbose,pbprint,natom,pbfrc,epsx,epsy,epsz,phi,cphi)

            ! second printing of DB forces, in kcal/mol-Angstrom

            pbfrc = frcfac*pbfrc
!           write(102,*) ' :::: Atomic DB forces ::::'
!           do iatm = 1, natom
!              write(102,'(3e20.6)') pbfrc(1:3,iatm)
!           end do
         else
            pbfrc = frcfac*pbfrc

            !fnet = ZERO
            !do iatm = 1, atmlast
            !   fnet(1) = fnet(1) + pbfrc(1, iatm)
            !   fnet(2) = fnet(2) + pbfrc(2, iatm)
            !   fnet(3) = fnet(3) + pbfrc(3, iatm)
            !end do
            !fnet = fnet/dble(atmlast)
            !do iatm = 1, atmlast
            !   pbfrc(1,iatm) = frcfac*(pbfrc(1,iatm) - fnet(1))
            !   pbfrc(2,iatm) = frcfac*(pbfrc(2,iatm) - fnet(2))
            !   pbfrc(3,iatm) = frcfac*(pbfrc(3,iatm) - fnet(3))
            !end do

         end if

      ! compute fdfrc by the charge option
      ! dbf is computed by Cai and Ye et al

      else if ( eneopt == 2 .and. epsin /= epsout ) then

         if (ifcap == 2 .or. ifcap == 5) then
            ! Only consider protein atoms for calc of total qE energy
            if(ibgion /= 0) then
               lastatm = ipres(ibgion) - 1
            else if(ibgwat /= 0) then
               lastatm = ipres(ibgwat) - 1
            else
               lastatm = atmlast
            end if
         else
            ! Consider all atoms
            lastatm = atmlast
         end if

         factor = AMBER_ELECTROSTATIC2/frcfac/(epsin/eps0)
         chgrd(1:xmymzm) = chgrd(1:xmymzm)*factor

         if ( frcopt == 2 ) then
            zv = -dble(insas) ! pseudo signed distance function
            call pb_dbfrc_crg ( pbverbose,pbprint,natom,eelrf,pbfrc, &
                                epsx,epsy,epsz,zv,phi,chgrd,cphi )

         else if ( frcopt == 3 ) then
            zv = -dble(insas) ! pseudo signed distance function
            call pb_dbfrc_fld2( pbverbose,pbprint,natom,eelrf,pbfrc, &
                                epsx,epsy,epsz,zv,phi,chgrd,cphi )
!           call pb_dbfrc_fld3( pbverbose,pbprint,natom,eelrf,pbfrc, &
!                               epsx,epsy,epsz,zv,phi,chgrd,cphi )
         else if ( frcopt == 4 ) then
            zv = -dble(insas) ! pseudo signed distance function
            call pb_dbfrc_crg2( pbverbose,pbprint,natom,eelrf,pbfrc, &
                                epsx,epsy,epsz,zv,phi,chgrd,cphi )
         else 
            call pb_dbene( pbverbose,pbprint,natom,lastatm,ifcap,npdec,idecomp,m,irespw, &
                           ipres,eelrf,pbfrc,insas,phi,chgrd,cphi )

         end if

      ! compute fd energy and force by the P3M option
      ! note that self forces are zero
      ! delete fd coulomb energy and forces for all close pairs
      ! delete fd reaction energy and forces for all close pairs

      else if ( eneopt == 3 ) then

!        print *, 'OK0', intopt
         if (ifcap == 2 .or. ifcap == 5) then
            ! Only consider protein atoms for calc of total qE energy
            if(ibgion /= 0) then
               lastatm = ipres(ibgion) - 1
            else if(ibgwat /= 0) then
               lastatm = ipres(ibgwat) - 1
            else
               lastatm = atmlast
            end if
         else
            ! Consider all atoms
            lastatm = atmlast
         end if

         factor = AMBER_ELECTROSTATIC2/frcfac/(epsin/eps0)
         chgrd(1:xmymzm) = chgrd(1:xmymzm)*factor

         ! compute total qE energy and forces and delete self energy, note that self forces are zero

         if ( intopt == 1 ) call pb_qefrc( natom, 1, atmlast, idecomp, eelrf, eelself, outflag, pbfrc, phi )
         zv = -dble(insas) ! pseudo signed distance function
         if ( intopt == 2 ) then
            if ( bcopt < 6 .or. bcopt > 9 ) then
               write(6,*) 'PB Bomb in pb_force(): intop = 2 does not support bcopt < 6'
               call mexit(6, 1)
            end if

            call pb_qefrc2( natom, 1, atmlast, eelrf, outflag, pbfrc, zv(1), phi )
         end if

         ! delete fd coulomb energy and forces for all close pairs
         ! when bcopt >= 6, we only have reaction field energy in eelrf
         ! no need to any corrections, coulombic energy should come from
         ! pb_direct ...

         !if ( bcopt < 6 ) then
         !   print *, 'call pb_fdcoulomb'
         !   call pb_fdcoulomb( natom, 1, atmlast, idecomp, eelcoul, outflag, pbfrc )
         !   eelrf = eelrf - eelself - eelcoul
         !end if
         !print *, 'fd react1 ', eelrf*frcfac, bcopt

         ! delete fd reaction energy and forces for all close pairs
         ! when bcopt >= 6, we only have reaction field energy in eelrf
         ! no need to any corrections, coulombic energy should come from
         ! pb_dbfrc_crg ...

         grdreac = ZERO
         drcreac = ZERO
         allocate (pos_crg(3,10000,natom))
         allocate (ipos_crg(3,10000,natom))
         allocate (surf_crg(10000,natom))
         allocate (crg_num(1:natom))
         if ( epsin /= epsout ) then
!           print *, 'OK1'
            call pb_crgview ( pbverbose,pbprint,natom,pbfrc,epsx,epsy,epsz,zv(1),phi,chgrd,cphi )
!           print *, 'OK2'
            call pb_fdreaction( natom, 1, atmlast, idecomp, grdreac, outflag, pbfrc )
            call pb_direct_reaction( natom, 1, atmlast, idecomp, drcreac, outflag, pbfrc )
         end if
         deallocate (pos_crg)
         deallocate (ipos_crg)
         deallocate (surf_crg)
         deallocate (crg_num)

         ! add ion contributions for nonlinear PBE ...

         if ( npbopt /= 0 ) then
            eelrf = eelrf + ionene
         end if

         eelrf = frcfac*eelrf
         write(6,'(a,f20.10)') 'fld', eelrf
         eelrf = eelrf - grdreac + drcreac
         write(6,'(a,f20.10)') 'p3m', eelrf 
!        write(5001,'(3f20.10)') frcfac*eelrf, grdreac, drcreac
!        write(5002,'(f20.10)') frcfac*eelrf - grdreac + drcreac

!        write(6,*) 'final eelrf/ionene in kcal/mol', frcfac*eelrf, frcfac*ionene
!        write(6,*) 'final eelself in kcal/mol', frcfac*eelself
!        write(6,*) 'final eelcoul in kcal/mol', frcfac*eelcoul

         if ( frcopt == 2 ) then
            zv = -dble(insas) ! pseudo signed distance function
            call pb_dbfrc_crg ( pbverbose,pbprint,natom,eelrf,pbfrc, &
                                epsx,epsy,epsz,zv(1),phi,chgrd,cphi )
         end if

      ! probably a uniform dielectric run for reference state calculation and
      ! development

      else
 
         eelrf = ZERO

      end if

   end do !-- PB pairwise decomp <<1,mpdec>>
 
end subroutine pb_fdfrc
subroutine mapoutviaidx(mymap,fn,lowi,highi,lowj,highj,lowk,highk)
   implicit none
   _REAL_ mymap(1:xm,1:ym,1:zm)
   integer i,j,k,fn,lowi,highi,lowj,highj,lowk,highk
   do i=lowi,highi
      do j=lowj,highj
         do k=lowk,highk
            write(fn,'(f20.6)') mymap(i,j,k)
            !write(fn,*) i*h+gox,j*h+goy,k*h+goz
   enddo; enddo; enddo
end subroutine mapoutviaidx
subroutine mapoutviacrd(mymap,fn,minx,maxx,miny,maxy,minz,maxz)
   implicit none
   _REAL_ mymap(1:xm,1:ym,1:zm)
   integer i,j,k,fn
   _REAL_ minx,maxx,miny,maxy,minz,maxz
   do i=1,xm
      if ( i*h+gox < minx .or. i*h+gox > maxx ) cycle
      do j=1,ym
         if ( j*h+goy < miny .or. j*h+goy > maxy ) cycle
         do k=1,zm
            if ( k*h+goz < minz .or. k*h+goz > maxz ) cycle
            !write(fn,*) mymap(i,j,k)
            write(fn,'(f20.6)') mymap(i,j,k)
            !write(fn,*) i*h+gox,j*h+goy,k*h+goz
   enddo; enddo; enddo
end subroutine mapoutviacrd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Map atomic charges onto grid
subroutine pb_crggrd ( xm, ym, zm, atmfirst, atmlast, icrd, gcrg, chgrd )

   integer xm, ym, zm, atmfirst, atmlast
   integer icrd(3,*)
   _REAL_ gcrg(8,*)
   _REAL_ chgrd(xm,ym,zm)

   integer iatm, i, j, k
   do iatm = atmfirst, atmlast
      if ( ligand .or. multiblock ) then
         if ( level == nfocus .and. realflag(iatm) == 0 ) cycle
      endif
      !if ( level == nfocus .and. ligand     .and. realflag(iatm) == 0 ) cycle
      !if ( level == nfocus .and. multiblock .and. realflag(iatm) == 0 ) cycle
      i = icrd(1, iatm); j = icrd(2, iatm); k = icrd(3, iatm)
      if ( i<1 .or. j<1 .or. k<1 .or. i+1>xm .or. j+1>ym .or. k+1>zm ) then
         write(6,*) i,j,k,xm,ym,zm
         write(6,*) gox,goy,goz
         write(6,*) "pb_crggrd: index exceeds box range, bail."
         call mexit(6,1)
      endif
      chgrd(i  ,j  ,k  ) = chgrd(i  ,j  ,k  ) + gcrg(1, iatm)
      chgrd(i+1,j  ,k  ) = chgrd(i+1,j  ,k  ) + gcrg(2, iatm)
      chgrd(i  ,j+1,k  ) = chgrd(i  ,j+1,k  ) + gcrg(3, iatm)
      chgrd(i+1,j+1,k  ) = chgrd(i+1,j+1,k  ) + gcrg(4, iatm)
      chgrd(i  ,j  ,k+1) = chgrd(i  ,j  ,k+1) + gcrg(5, iatm)
      chgrd(i+1,j  ,k+1) = chgrd(i+1,j  ,k+1) + gcrg(6, iatm)
      chgrd(i  ,j+1,k+1) = chgrd(i  ,j+1,k+1) + gcrg(7, iatm)
      chgrd(i+1,j+1,k+1) = chgrd(i+1,j+1,k+1) + gcrg(8, iatm)
   end do

end subroutine pb_crggrd
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ total finite difference es energy and forces for intopt = 1
subroutine pb_qefrc( natom, atmfirst, atmlast, idecomp, grdnrg, grdself, outflag, pbfrc, phi )
   
   use decomp, only: decpair

   ! Common variables
   
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
   
   ! Passed variables
   
   integer natom, atmfirst, atmlast, idecomp
   integer outflag(natom)
   _REAL_ grdnrg, grdself
   _REAL_ pbfrc(3, natom)
   _REAL_ phi(xm,ym,zm)
   
   ! Local variables
   
   integer iatm
   integer i, j, k
   _REAL_ :: g000, g100, g110, g111
   _REAL_ :: gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ :: grdnrgtmp, grdselftmp, hpsnrv
   
   ! begin code
   
   g000 = INV_FOURPI*green(0,0,0)
   g100 = INV_FOURPI*green(1,0,0)
   g110 = INV_FOURPI*green(1,1,0)
   g111 = INV_FOURPI*green(1,1,1)
    
   hpsnrv = ONE/(h*epsin)
    
   grdnrg = ZERO
   grdself = ZERO
    
   ! split each atoms charge over the eight surrounding
   ! grid points according to the trilinear weighting
   ! function and add up each of the contributions.
!  if ( ligand ) then
!    do iatm = 1,natom
!       write(1234,'(i4)') realflag(iatm)
!    enddo
!    continue
!  else
!    realflag = 0
!    do iatm = 1 , natom
!      read(1234,'(i4)') realflag(iatm)
!    enddo
!    rewind(1234)
!    ligand = .true.
!  endif
      
   do iatm = atmfirst, atmlast
      if ( ligand .or. multiblock ) then
         if ( liveflag(iatm) == 0 ) cycle
      !if ( ligand .and. liveflag(iatm) == 0 ) then
      !   cycle
      !else if ( multiblock .and. liveflag(iatm) == 0 ) then
      !   cycle
      else if ( outflag(iatm) == 1 ) then 
         cycle
      end if
      i = icrd(1,iatm); j = icrd(2,iatm); k = icrd(3,iatm)
      gci1 = gcrg(1,iatm); gci2 = gcrg(2,iatm)
      gci3 = gcrg(3,iatm); gci4 = gcrg(4,iatm)
      gci5 = gcrg(5,iatm); gci6 = gcrg(6,iatm)
      gci7 = gcrg(7,iatm); gci8 = gcrg(8,iatm)
      
      grdnrgtmp = &
         gci1*phi(i  ,j  ,k  ) + gci2*phi(i+1,j  ,k  ) + &
         gci3*phi(i  ,j+1,k  ) + gci4*phi(i+1,j+1,k  ) + &
         gci5*phi(i  ,j  ,k+1) + gci6*phi(i+1,j  ,k+1) + &
         gci7*phi(i  ,j+1,k+1) + gci8*phi(i+1,j+1,k+1)

      grdnrg = grdnrg + grdnrgtmp
      
      grdselftmp = &
         g000 * (gci1*gci1 + gci2*gci2 + gci3*gci3 + gci4*gci4 + &
                 gci5*gci5 + gci6*gci6 + gci7*gci7 + gci8*gci8 )*HALF + &
         g100 * (gci1*gci2 + gci1*gci3 + gci1*gci5 + gci2*gci4 + &
                 gci2*gci6 + gci4*gci3 + gci4*gci8 + gci3*gci7 + &
                 gci5*gci6 + gci5*gci7 + gci6*gci8 + gci8*gci7 ) + &
         g110 * (gci1*gci4 + gci1*gci6 + gci1*gci7 + gci2*gci3 + &
                 gci2*gci5 + gci2*gci8 + gci4*gci6 + gci4*gci7 + &
                 gci3*gci5 + gci3*gci8 + gci5*gci8 + gci6*gci7 ) + &
         g111 * (gci1*gci8 + gci2*gci7 + gci4*gci5 + gci6*gci3 )

      grdself = grdself + grdselftmp

      !-- PB decomp
      if(idecomp == 1 .or. idecomp == 2) then
         grdnrgtmp = frcfac*(HALF*grdnrgtmp - grdselftmp*hpsnrv)
         call decpair(1,iatm,iatm,grdnrgtmp)
      end if

      if (frcopt /= 0) then
         pbfrc(1,iatm) = &
         gci1*ex(i  ,j  ,k  ) + gci2*ex(i+1,j  ,k  ) + &
         gci3*ex(i  ,j+1,k  ) + gci4*ex(i+1,j+1,k  ) + &
         gci5*ex(i  ,j  ,k+1) + gci6*ex(i+1,j  ,k+1) + &
         gci7*ex(i  ,j+1,k+1) + gci8*ex(i+1,j+1,k+1)
         pbfrc(2,iatm) = &
         gci1*ey(i  ,j  ,k  ) + gci2*ey(i+1,j  ,k  ) + &
         gci3*ey(i  ,j+1,k  ) + gci4*ey(i+1,j+1,k  ) + &
         gci5*ey(i  ,j  ,k+1) + gci6*ey(i+1,j  ,k+1) + &
         gci7*ey(i  ,j+1,k+1) + gci8*ey(i+1,j+1,k+1) 
         pbfrc(3,iatm) = &
         gci1*ez(i  ,j  ,k  ) + gci2*ez(i+1,j  ,k  ) + &
         gci3*ez(i  ,j+1,k  ) + gci4*ez(i+1,j+1,k  ) + &
         gci5*ez(i  ,j  ,k+1) + gci6*ez(i+1,j  ,k+1) + &
         gci7*ez(i  ,j+1,k+1) + gci8*ez(i+1,j+1,k+1) 
      end if
   end do
  
   grdnrg  = HALF*grdnrg
   grdself = grdself/( h*epsin )

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ e field vector x
   _REAL_ function ex(i, j, k)

   integer, intent(in) :: i, j, k

   ex = ( phi(i-1,j  ,k  ) - phi(i+1,j  ,k  ) )/(2*h)

end function ex
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ e field vector y
_REAL_ function ey(i, j, k)

   integer, intent(in) :: i, j, k

   ey = ( phi(i  ,j-1,k  ) - phi(i  ,j+1,k  ) )/(2*h)

end function ey
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ e field vector z
_REAL_ function ez(i, j, k)

   integer, intent(in) :: i, j, k

   ez = ( phi(i  ,j  ,k-1) - phi(i  ,j  ,k+1) )/(2*h)

end function ez

end subroutine pb_qefrc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ total finite difference es energy and forces for intopt = 2
subroutine pb_qefrc2( natom, atmfirst, atmlast, grdnrg, outflag, pbfrc, insas, phi )
   
   ! Passed variables
   
   integer natom, atmfirst, atmlast
   integer outflag(natom)
   _REAL_ grdnrg
   _REAL_ pbfrc(3, natom)
   _REAL_ insas(0:xm+1, 0:ym+1, 0:zm+1), phi(xm,ym,zm)
   
   ! Local variables
   
   integer iatm
   integer i, j, k
   _REAL_ :: gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ :: fx0, fy0, fz0, charge
   _REAL_ :: up, dudxi0, dudyi0, dudzi0
   integer, parameter :: n_point = 8
    
   ! begin code
    
   grdnrg = ZERO
    
   ! split each atoms charge over the eight surrounding
   ! grid points according to the trilinear weighting
   ! function and add up each of the contributions.
    
   do iatm = atmfirst, atmlast
      if ( ligand .or. multiblock ) then
         if ( liveflag(iatm) == 0 ) cycle
      !if ( ligand .and. liveflag(iatm) == 0 ) then
      !   cycle
      !else if ( multiblock .and. liveflag(iatm) == 0 ) then
      !   cycle
      else if ( outflag(iatm) == 1 ) then
         cycle
      end if
      i = icrd(1,iatm); j = icrd(2,iatm); k = icrd(3,iatm)
      fx0 = gcrd(1,iatm); fy0 = gcrd(2,iatm); fz0 = gcrd(3,iatm)
      charge = acrg(iatm)
       
      call gradu(xm,ym,zm,-ONE,n_point,4,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,insas)
       
      grdnrg = grdnrg + charge*up
       
      pbfrc(1,iatm) = pbfrc(1,iatm) - charge*dudxi0
      pbfrc(2,iatm) = pbfrc(2,iatm) - charge*dudyi0
      pbfrc(3,iatm) = pbfrc(3,iatm) - charge*dudzi0
   end do
    
   grdnrg = HALF*grdnrg

end subroutine pb_qefrc2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ FD coulombic energy and forces.
subroutine pb_fdcoulomb( natom, atmfirst, atmlast, idecomp, grdcoul, outflag, pbfrc )
   
   use decomp, only: decpair

   ! Common variables
   
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
   
   ! Passed variables
   
   integer natom, atmfirst, atmlast, idecomp
   integer outflag(natom)
   _REAL_ grdcoul
   _REAL_ pbfrc(3, natom)
   
   ! Local variables
   
   integer iatm, jatm, jp, ilast, jfirst, jlast
   integer dijx, dijx0, dijx1, dijx2, dijx3, dijy, dijy0, dijy1, dijy2, dijy3, dijz, dijz0, dijz1, dijz2, dijz3
   integer ix, iy, iz, jx, jy, jz
   _REAL_ gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ gcj1, gcj2, gcj3, gcj4, gcj5, gcj6, gcj7, gcj8
   _REAL_ gcij(27)
   _REAL_ factor, factor1, decfac
   _REAL_ ffx, ffy, ffz
   _REAL_ dumx, dumy, dumz
   _REAL_ frc(3, natom)
   _REAL_ grdcoultmp
   _REAL_ pair_correct
   
   ! begin code
   
   factor  = ONE/( FOURPI*epsin*h )
   factor1 = HALF/( FOURPI*epsin*h*h )
   decfac  = -factor*frcfac
   
   grdcoul = ZERO
   frc = ZERO
   do iatm = atmfirst, atmlast
      if ( ligand .or. multiblock ) then
         if ( liveflag(iatm) == 0 ) cycle
      !if ( ligand .and. liveflag(iatm) == 0 ) then
      !   cycle
      !!else if ( multiblock .and. realflag(iatm) == 0 ) then
      !!   cycle ! this is only for debugging
      !else if ( multiblock .and. liveflag(iatm) == 0 ) then
      !   cycle
      else if ( outflag(iatm) == 1 ) then
         cycle
      end if
      ix = icrd(1,iatm); iy = icrd(2,iatm); iz = icrd(3,iatm)
      
      gci1 = gcrg(1,iatm); gci2 = gcrg(2,iatm)
      gci3 = gcrg(3,iatm); gci4 = gcrg(4,iatm)
      gci5 = gcrg(5,iatm); gci6 = gcrg(6,iatm)
      gci7 = gcrg(7,iatm); gci8 = gcrg(8,iatm)
      
      jfirst = iar1pb(4, iatm-1) + 1
      jlast  = iar1pb(2, iatm)
      dumx = ZERO; dumy = ZERO; dumz = ZERO
      do jp = jfirst, jlast
         jatm = iprshrt(jp)
         pair_correct = ONE
!         if ( multiblock .and. liveflag(jatm) == 0 ) then
!            pair_correct = HALF
!write(6666,*) iatm,jatm
!         endif
         jx = icrd(1,jatm); jy = icrd(2,jatm); jz = icrd(3,jatm)
         
         gcj1 = gcrg(1,jatm); gcj2 = gcrg(2,jatm)
         gcj3 = gcrg(3,jatm); gcj4 = gcrg(4,jatm)
         gcj5 = gcrg(5,jatm); gcj6 = gcrg(6,jatm)
         gcj7 = gcrg(7,jatm); gcj8 = gcrg(8,jatm)
         
         dijx  =       ix-jx; dijy  =       iy-jy; dijz  =       iz-jz
         dijx0 = abs(dijx-2); dijy0 = abs(dijy-2); dijz0 = abs(dijz-2)
         dijx1 = abs(dijx-1); dijy1 = abs(dijy-1); dijz1 = abs(dijz-1)
         dijx2 = abs(dijx+1); dijy2 = abs(dijy+1); dijz2 = abs(dijz+1)
         dijx3 = abs(dijx+2); dijy3 = abs(dijy+2); dijz3 = abs(dijz+2)
         dijx  = abs(dijx  ); dijy  = abs(dijy  ); dijz  = abs(dijz  )
         
         gcij( 1) = gci1*gcj1 + gci2*gcj2 + gci3*gcj3 + gci4*gcj4 + gci5*gcj5 + gci6*gcj6 + gci7*gcj7 + gci8*gcj8
         gcij( 2) = gci1*gcj2 + gci3*gcj4 + gci5*gcj6 + gci7*gcj8
         gcij( 3) = gci1*gcj3 + gci2*gcj4 + gci5*gcj7 + gci6*gcj8
         gcij( 4) = gci1*gcj4 + gci5*gcj8
         gcij( 5) = gci1*gcj5 + gci2*gcj6 + gci3*gcj7 + gci4*gcj8
         gcij( 6) = gci1*gcj6 + gci3*gcj8
         gcij( 7) = gci1*gcj7 + gci2*gcj8
         gcij( 8) = gci1*gcj8
         gcij( 9) = gci2*gcj1 + gci4*gcj3 + gci6*gcj5 + gci8*gcj7
         gcij(10) = gci2*gcj3 + gci6*gcj7
         gcij(11) = gci2*gcj5 + gci4*gcj7
         gcij(12) = gci2*gcj7
         gcij(13) = gci3*gcj1 + gci4*gcj2 + gci7*gcj5 + gci8*gcj6
         gcij(14) = gci3*gcj2 + gci7*gcj6
         gcij(15) = gci3*gcj5 + gci4*gcj6
         gcij(16) = gci3*gcj6
         gcij(17) = gci4*gcj1 + gci8*gcj5
         gcij(18) = gci4*gcj5
         gcij(19) = gci5*gcj1 + gci6*gcj2 + gci7*gcj3 + gci8*gcj4
         gcij(20) = gci5*gcj2 + gci7*gcj4
         gcij(21) = gci5*gcj3 + gci6*gcj4
         gcij(22) = gci5*gcj4
         gcij(23) = gci6*gcj1 + gci8*gcj3
         gcij(24) = gci6*gcj3
         gcij(25) = gci7*gcj1 + gci8*gcj2
         gcij(26) = gci7*gcj2
         gcij(27) = gci8*gcj1

         grdcoultmp = dble( &
              gci1*( gcj1*l_green(dijx ,dijy ,dijz ) + gcj2*l_green(dijx1,dijy ,dijz ) + &
                     gcj3*l_green(dijx ,dijy1,dijz ) + gcj4*l_green(dijx1,dijy1,dijz ) + &
                     gcj5*l_green(dijx ,dijy ,dijz1) + gcj6*l_green(dijx1,dijy ,dijz1) + &
                     gcj7*l_green(dijx ,dijy1,dijz1) + gcj8*l_green(dijx1,dijy1,dijz1) ) + &
              gci2*( gcj1*l_green(dijx2,dijy ,dijz ) + gcj2*l_green(dijx ,dijy ,dijz ) + &
                     gcj3*l_green(dijx2,dijy1,dijz ) + gcj4*l_green(dijx ,dijy1,dijz ) + &
                     gcj5*l_green(dijx2,dijy ,dijz1) + gcj6*l_green(dijx ,dijy ,dijz1) + &
                     gcj7*l_green(dijx2,dijy1,dijz1) + gcj8*l_green(dijx ,dijy1,dijz1) ) + &
              gci3*( gcj1*l_green(dijx ,dijy2,dijz ) + gcj2*l_green(dijx1,dijy2,dijz ) + &
                     gcj3*l_green(dijx ,dijy ,dijz ) + gcj4*l_green(dijx1,dijy ,dijz ) + &
                     gcj5*l_green(dijx ,dijy2,dijz1) + gcj6*l_green(dijx1,dijy2,dijz1) + &
                     gcj7*l_green(dijx ,dijy ,dijz1) + gcj8*l_green(dijx1,dijy ,dijz1) ) + &
              gci4*( gcj1*l_green(dijx2,dijy2,dijz ) + gcj2*l_green(dijx ,dijy2,dijz ) + &
                     gcj3*l_green(dijx2,dijy ,dijz ) + gcj4*l_green(dijx ,dijy ,dijz ) + &
                     gcj5*l_green(dijx2,dijy2,dijz1) + gcj6*l_green(dijx ,dijy2,dijz1) + &
                     gcj7*l_green(dijx2,dijy ,dijz1) + gcj8*l_green(dijx ,dijy ,dijz1) ) + &
              gci5*( gcj1*l_green(dijx ,dijy ,dijz2) + gcj2*l_green(dijx1,dijy ,dijz2) + &
                     gcj3*l_green(dijx ,dijy1,dijz2) + gcj4*l_green(dijx1,dijy1,dijz2) + &
                     gcj5*l_green(dijx ,dijy ,dijz ) + gcj6*l_green(dijx1,dijy ,dijz ) + &
                     gcj7*l_green(dijx ,dijy1,dijz ) + gcj8*l_green(dijx1,dijy1,dijz ) ) + &
              gci6*( gcj1*l_green(dijx2,dijy ,dijz2) + gcj2*l_green(dijx ,dijy ,dijz2) + &
                     gcj3*l_green(dijx2,dijy1,dijz2) + gcj4*l_green(dijx ,dijy1,dijz2) + &
                     gcj5*l_green(dijx2,dijy ,dijz ) + gcj6*l_green(dijx ,dijy ,dijz ) + &
                     gcj7*l_green(dijx2,dijy1,dijz ) + gcj8*l_green(dijx ,dijy1,dijz ) ) + &
              gci7*( gcj1*l_green(dijx ,dijy2,dijz2) + gcj2*l_green(dijx1,dijy2,dijz2) + &
                     gcj3*l_green(dijx ,dijy ,dijz2) + gcj4*l_green(dijx1,dijy ,dijz2) + &
                     gcj5*l_green(dijx ,dijy2,dijz ) + gcj6*l_green(dijx1,dijy2,dijz ) + &
                     gcj7*l_green(dijx ,dijy ,dijz ) + gcj8*l_green(dijx1,dijy ,dijz ) ) + &
              gci8*( gcj1*l_green(dijx2,dijy2,dijz2) + gcj2*l_green(dijx ,dijy2,dijz2) + &
                     gcj3*l_green(dijx2,dijy ,dijz2) + gcj4*l_green(dijx ,dijy ,dijz2) + &
                     gcj5*l_green(dijx2,dijy2,dijz ) + gcj6*l_green(dijx ,dijy2,dijz ) + &
                     gcj7*l_green(dijx2,dijy ,dijz ) + gcj8*l_green(dijx ,dijy ,dijz ) ) )

         grdcoul = grdcoul + grdcoultmp*pair_correct

         !-- PB decomp
         if(idecomp == 1 .or. idecomp == 2) call decpair(1,iatm,jatm,decfac*grdcoultmp)
        
         ffx = dble( gcij( 1)*(l_green(dijx1,dijy ,dijz ) - l_green(dijx2,dijy ,dijz )) + &
                     gcij( 2)*(l_green(dijx0,dijy ,dijz ) - l_green(dijx ,dijy ,dijz )) + &
                     gcij( 3)*(l_green(dijx1,dijy1,dijz ) - l_green(dijx2,dijy1,dijz )) + &
                     gcij( 4)*(l_green(dijx0,dijy1,dijz ) - l_green(dijx ,dijy1,dijz )) + &
                     gcij( 5)*(l_green(dijx1,dijy ,dijz1) - l_green(dijx2,dijy ,dijz1)) + &
                     gcij( 6)*(l_green(dijx0,dijy ,dijz1) - l_green(dijx ,dijy ,dijz1)) + &
                     gcij( 7)*(l_green(dijx1,dijy1,dijz1) - l_green(dijx2,dijy1,dijz1)) + &
                     gcij( 8)*(l_green(dijx0,dijy1,dijz1) - l_green(dijx ,dijy1,dijz1)) + &
                     gcij( 9)*(l_green(dijx ,dijy ,dijz ) - l_green(dijx3,dijy ,dijz )) + &
                     gcij(10)*(l_green(dijx ,dijy1,dijz ) - l_green(dijx3,dijy1,dijz )) + &
                     gcij(11)*(l_green(dijx ,dijy ,dijz1) - l_green(dijx3,dijy ,dijz1)) + &
                     gcij(12)*(l_green(dijx ,dijy1,dijz1) - l_green(dijx3,dijy1,dijz1)) + &
                     gcij(13)*(l_green(dijx1,dijy2,dijz ) - l_green(dijx2,dijy2,dijz )) + &
                     gcij(14)*(l_green(dijx0,dijy2,dijz ) - l_green(dijx ,dijy2,dijz )) + &
                     gcij(15)*(l_green(dijx1,dijy2,dijz1) - l_green(dijx2,dijy2,dijz1)) + &
                     gcij(16)*(l_green(dijx0,dijy2,dijz1) - l_green(dijx ,dijy2,dijz1)) + &
                     gcij(17)*(l_green(dijx ,dijy2,dijz ) - l_green(dijx3,dijy2,dijz )) + &
                     gcij(18)*(l_green(dijx ,dijy2,dijz1) - l_green(dijx3,dijy2,dijz1)) + &
                     gcij(19)*(l_green(dijx1,dijy ,dijz2) - l_green(dijx2,dijy ,dijz2)) + &
                     gcij(20)*(l_green(dijx0,dijy ,dijz2) - l_green(dijx ,dijy ,dijz2)) + &
                     gcij(21)*(l_green(dijx1,dijy1,dijz2) - l_green(dijx2,dijy1,dijz2)) + &
                     gcij(22)*(l_green(dijx0,dijy1,dijz2) - l_green(dijx ,dijy1,dijz2)) + &
                     gcij(23)*(l_green(dijx ,dijy ,dijz2) - l_green(dijx3,dijy ,dijz2)) + &
                     gcij(24)*(l_green(dijx ,dijy1,dijz2) - l_green(dijx3,dijy1,dijz2)) + &
                     gcij(25)*(l_green(dijx1,dijy2,dijz2) - l_green(dijx2,dijy2,dijz2)) + &
                     gcij(26)*(l_green(dijx0,dijy2,dijz2) - l_green(dijx ,dijy2,dijz2)) + &
                     gcij(27)*(l_green(dijx ,dijy2,dijz2) - l_green(dijx3,dijy2,dijz2)) )
       
         ffy = dble( gcij( 1)*(l_green(dijx ,dijy1,dijz ) - l_green(dijx ,dijy2,dijz )) + &
                     gcij( 2)*(l_green(dijx1,dijy1,dijz ) - l_green(dijx1,dijy2,dijz )) + &
                     gcij( 3)*(l_green(dijx ,dijy0,dijz ) - l_green(dijx ,dijy ,dijz )) + &
                     gcij( 4)*(l_green(dijx1,dijy0,dijz ) - l_green(dijx1,dijy ,dijz )) + &
                     gcij( 5)*(l_green(dijx ,dijy1,dijz1) - l_green(dijx ,dijy2,dijz1)) + &
                     gcij( 6)*(l_green(dijx1,dijy1,dijz1) - l_green(dijx1,dijy2,dijz1)) + &
                     gcij( 7)*(l_green(dijx ,dijy0,dijz1) - l_green(dijx ,dijy ,dijz1)) + &
                     gcij( 8)*(l_green(dijx1,dijy0,dijz1) - l_green(dijx1,dijy ,dijz1)) + &
                     gcij( 9)*(l_green(dijx2,dijy1,dijz ) - l_green(dijx2,dijy2,dijz )) + &
                     gcij(10)*(l_green(dijx2,dijy0,dijz ) - l_green(dijx2,dijy ,dijz )) + &
                     gcij(11)*(l_green(dijx2,dijy1,dijz1) - l_green(dijx2,dijy2,dijz1)) + &
                     gcij(12)*(l_green(dijx2,dijy0,dijz1) - l_green(dijx2,dijy ,dijz1)) + &
                     gcij(13)*(l_green(dijx ,dijy ,dijz ) - l_green(dijx ,dijy3,dijz )) + &
                     gcij(14)*(l_green(dijx1,dijy ,dijz ) - l_green(dijx1,dijy3,dijz )) + &
                     gcij(15)*(l_green(dijx ,dijy ,dijz1) - l_green(dijx ,dijy3,dijz1)) + &
                     gcij(16)*(l_green(dijx1,dijy ,dijz1) - l_green(dijx1,dijy3,dijz1)) + &
                     gcij(17)*(l_green(dijx2,dijy ,dijz ) - l_green(dijx2,dijy3,dijz )) + &
                     gcij(18)*(l_green(dijx2,dijy ,dijz1) - l_green(dijx2,dijy3,dijz1)) + &
                     gcij(19)*(l_green(dijx ,dijy1,dijz2) - l_green(dijx ,dijy2,dijz2)) + &
                     gcij(20)*(l_green(dijx1,dijy1,dijz2) - l_green(dijx1,dijy2,dijz2)) + &
                     gcij(21)*(l_green(dijx ,dijy0,dijz2) - l_green(dijx ,dijy ,dijz2)) + &
                     gcij(22)*(l_green(dijx1,dijy0,dijz2) - l_green(dijx1,dijy ,dijz2)) + &
                     gcij(23)*(l_green(dijx2,dijy1,dijz2) - l_green(dijx2,dijy2,dijz2)) + &
                     gcij(24)*(l_green(dijx2,dijy0,dijz2) - l_green(dijx2,dijy ,dijz2)) + &
                     gcij(25)*(l_green(dijx ,dijy ,dijz2) - l_green(dijx ,dijy3,dijz2)) + &
                     gcij(26)*(l_green(dijx1,dijy ,dijz2) - l_green(dijx1,dijy3,dijz2)) + &
                     gcij(27)*(l_green(dijx2,dijy ,dijz2) - l_green(dijx2,dijy3,dijz2)) )
      
         ffz = dble( gcij( 1)*(l_green(dijx ,dijy ,dijz1) - l_green(dijx ,dijy ,dijz2)) + &
                     gcij( 2)*(l_green(dijx1,dijy ,dijz1) - l_green(dijx1,dijy ,dijz2)) + &
                     gcij( 3)*(l_green(dijx ,dijy1,dijz1) - l_green(dijx ,dijy1,dijz2)) + &
                     gcij( 4)*(l_green(dijx1,dijy1,dijz1) - l_green(dijx1,dijy1,dijz2)) + &
                     gcij( 5)*(l_green(dijx ,dijy ,dijz0) - l_green(dijx ,dijy ,dijz )) + &
                     gcij( 6)*(l_green(dijx1,dijy ,dijz0) - l_green(dijx1,dijy ,dijz )) + &
                     gcij( 7)*(l_green(dijx ,dijy1,dijz0) - l_green(dijx ,dijy1,dijz )) + &
                     gcij( 8)*(l_green(dijx1,dijy1,dijz0) - l_green(dijx1,dijy1,dijz )) + &
                     gcij( 9)*(l_green(dijx2,dijy ,dijz1) - l_green(dijx2,dijy ,dijz2)) + &
                     gcij(10)*(l_green(dijx2,dijy1,dijz1) - l_green(dijx2,dijy1,dijz2)) + &
                     gcij(11)*(l_green(dijx2,dijy ,dijz0) - l_green(dijx2,dijy ,dijz )) + &
                     gcij(12)*(l_green(dijx2,dijy1,dijz0) - l_green(dijx2,dijy1,dijz )) + &
                     gcij(13)*(l_green(dijx ,dijy2,dijz1) - l_green(dijx ,dijy2,dijz2)) + &
                     gcij(14)*(l_green(dijx1,dijy2,dijz1) - l_green(dijx1,dijy2,dijz2)) + &
                     gcij(15)*(l_green(dijx ,dijy2,dijz0) - l_green(dijx ,dijy2,dijz )) + &
                     gcij(16)*(l_green(dijx1,dijy2,dijz0) - l_green(dijx1,dijy2,dijz )) + &
                     gcij(17)*(l_green(dijx2,dijy2,dijz1) - l_green(dijx2,dijy2,dijz2)) + &
                     gcij(18)*(l_green(dijx2,dijy2,dijz0) - l_green(dijx2,dijy2,dijz )) + &
                     gcij(19)*(l_green(dijx ,dijy ,dijz ) - l_green(dijx ,dijy ,dijz3)) + &
                     gcij(20)*(l_green(dijx1,dijy ,dijz ) - l_green(dijx1,dijy ,dijz3)) + &
                     gcij(21)*(l_green(dijx ,dijy1,dijz ) - l_green(dijx ,dijy1,dijz3)) + &
                     gcij(22)*(l_green(dijx1,dijy1,dijz ) - l_green(dijx1,dijy1,dijz3)) + &
                     gcij(23)*(l_green(dijx2,dijy ,dijz ) - l_green(dijx2,dijy ,dijz3)) + &
                     gcij(24)*(l_green(dijx2,dijy1,dijz ) - l_green(dijx2,dijy1,dijz3)) + &
                     gcij(25)*(l_green(dijx ,dijy2,dijz ) - l_green(dijx ,dijy2,dijz3)) + &
                     gcij(26)*(l_green(dijx1,dijy2,dijz ) - l_green(dijx1,dijy2,dijz3)) + &
                     gcij(27)*(l_green(dijx2,dijy2,dijz ) - l_green(dijx2,dijy2,dijz3)) )
      
         dumx = dumx + ffx; dumy = dumy + ffy; dumz = dumz + ffz
         frc(1,jatm) = frc(1,jatm) - ffx
         frc(2,jatm) = frc(2,jatm) - ffy
         frc(3,jatm) = frc(3,jatm) - ffz
      end do  !  jp = jfirst, jlast
      
      frc(1,iatm) = frc(1,iatm) + dumx
      frc(2,iatm) = frc(2,iatm) + dumy
      frc(3,iatm) = frc(3,iatm) + dumz
   end do  !  iatm = 1, ilast
   
   grdcoul = factor*grdcoul
   pbfrc  = pbfrc - factor1*frc
   
contains

function l_green (i,j,k)

   implicit none
   integer i,j,k
   _REAL_ l_green

   if ( i <= 20  .and. j <= 20 .and. k <= 20 ) then
      l_green = green(i,j,k)
   else
      l_green = ONE/sqrt(dble(i*i+j*j+k*k))
   end if

end function l_green

end subroutine pb_fdcoulomb 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ FD coulombic energy and forces.
subroutine pb_fdreaction( natom, atmfirst, atmlast, idecomp, grdreac, outflag, pbfrc )
   
   use decomp, only: decpair

   ! Common variables
   
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
   
   ! Passed variables
   
   integer natom, atmfirst, atmlast, idecomp
   integer outflag(natom)
   _REAL_ grdreac
   _REAL_ pbfrc(3, natom)
   
   ! Local variables
   
   integer iatm, jatm, jp, ilast, jfirst, jlast
   integer dijx, dijx0, dijx1, dijx2, dijx3, dijy, dijy0, dijy1, dijy2, dijy3, dijz, dijz0, dijz1, dijz2, dijz3
   integer ix, iy, iz, jx, jy, jz
   integer iix, iiy, iiz, jjx, jjy, jjz, i, j
   _REAL_ gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ gcj1, gcj2, gcj3, gcj4, gcj5, gcj6, gcj7, gcj8
   _REAL_ gcij(27)
   _REAL_ factor, factor1, decfac
   _REAL_ ffx, ffy, ffz
   _REAL_ dumx, dumy, dumz
   _REAL_ frc(3, natom)
   _REAL_ grdreactmp
   
   integer atmflag, crgflag
   integer isrf, jsrf, icrg, jcrg
   !allocate (pos_crg(3,10000,natom))
   !allocate (ipos_crg(3,10000,natom))
   !allocate (surf_crg(10000,natom))
   !_REAL_ pos_crg(3,10000,natom)
   !integer ipos_crg(3,10000,natom)
   !_REAL_ surf_crg(10000,natom)

   ! begin code
   
   factor  = ONE/( FOURPI*epsin*h )
   factor1 = HALF/( FOURPI*epsin*h*h )
   decfac  = -factor*frcfac
   
   grdreac = ZERO
   frc = ZERO

   ! first loop over iatm's grid charges with iatm's surface charges
   do iatm = atmfirst, atmlast
!     if ( ligand .and. realflag(iatm) == 0 ) then
!        cycle
!     else if ( multiblock .and. realflag(iatm) == 0 ) then
!        cycle ! this is only for debugging
!     else if ( multiblock .and. liveflag(iatm) == 0 ) then
!        cycle
!     else if ( outflag(iatm) == 1 ) then
!        cycle
!     end if
      ix = icrd(1,iatm); iy = icrd(2,iatm); iz = icrd(3,iatm)

      ! loop over iatm's grid charges over iatm's surface charges 

      grdreactmp = ZERO
      do isrf = 1, crg_num(iatm)
            jx = ipos_crg(1, isrf, iatm)
            jy = ipos_crg(2, isrf, iatm)
            jz = ipos_crg(3, isrf, iatm)
            dijx  =       jx-ix; dijy  =       jy-iy; dijz  =       jz-iz
            dijx1 = abs(dijx-1); dijy1 = abs(dijy-1); dijz1 = abs(dijz-1)
            dijx  = abs(dijx  ); dijy  = abs(dijy  ); dijz  = abs(dijz  )
            grdreactmp = grdreactmp + &
            surf_crg(isrf,iatm)*( gcrg(1,iatm)*l_green(dijx ,dijy ,dijz ) + gcrg(2,iatm)*l_green(dijx1,dijy ,dijz ) + &
                                  gcrg(3,iatm)*l_green(dijx ,dijy1,dijz ) + gcrg(4,iatm)*l_green(dijx1,dijy1,dijz ) + &
                                  gcrg(5,iatm)*l_green(dijx ,dijy ,dijz1) + gcrg(6,iatm)*l_green(dijx1,dijy ,dijz1) + &
                                  gcrg(7,iatm)*l_green(dijx ,dijy1,dijz1) + gcrg(8,iatm)*l_green(dijx1,dijy1,dijz1) ) 

      end do
         grdreac = grdreac + grdreactmp

   end do
   ! second loop over iatm's grid charges with iatm's surface charges

   do iatm = atmfirst, atmlast
!     if ( ligand .and. realflag(iatm) == 0 ) then
!        cycle
!     else if ( multiblock .and. realflag(iatm) == 0 ) then
!        cycle ! this is only for debugging
!     else if ( multiblock .and. liveflag(iatm) == 0 ) then
!        cycle
!     else if ( outflag(iatm) == 1 ) then
!        cycle
!     end if
      ix = icrd(1,iatm); iy = icrd(2,iatm); iz = icrd(3,iatm)

      jfirst = iar1pb(4, iatm-1) + 1
      jlast  = iar1pb(2, iatm)
      dumx = ZERO; dumy = ZERO; dumz = ZERO
      do jp = jfirst, jlast
         jatm = iprshrt(jp)
         if( jatm > atmlast ) cycle
         jx = icrd(1,jatm); jy = icrd(2,jatm); jz = icrd(3,jatm)
        
         ! loop over jatm's grid charges over iatm's surface charges 

         grdreactmp = ZERO
         do isrf = 1, crg_num(iatm)
               iix = ipos_crg(1, isrf, iatm)
               iiy = ipos_crg(2, isrf, iatm)
               iiz = ipos_crg(3, isrf, iatm)
            dijx  =       iix-jx; dijy  =       iiy-jy; dijz  =       iiz-jz
            dijx1 = abs(dijx-1); dijy1 = abs(dijy-1); dijz1 = abs(dijz-1)
            dijx  = abs(dijx  ); dijy  = abs(dijy  ); dijz  = abs(dijz  )
            grdreactmp = grdreactmp + &
            surf_crg(isrf,iatm)*( gcrg(1,jatm)*l_green(dijx ,dijy ,dijz ) + gcrg(2,jatm)*l_green(dijx1,dijy ,dijz ) + &
                                  gcrg(3,jatm)*l_green(dijx ,dijy1,dijz ) + gcrg(4,jatm)*l_green(dijx1,dijy1,dijz ) + &
                                  gcrg(5,jatm)*l_green(dijx ,dijy ,dijz1) + gcrg(6,jatm)*l_green(dijx1,dijy ,dijz1) + &
                                  gcrg(7,jatm)*l_green(dijx ,dijy1,dijz1) + gcrg(8,jatm)*l_green(dijx1,dijy1,dijz1) ) 

         end do
         grdreac = grdreac + grdreactmp

         ! loop over iatm's grid charges over jatm's surface charges 

         grdreactmp = ZERO
         do jsrf = 1, crg_num(jatm)
               jjx = ipos_crg(1, jsrf, jatm)
               jjy = ipos_crg(2, jsrf, jatm)
               jjz = ipos_crg(3, jsrf, jatm)
            dijx  =       jjx-ix; dijy  =       jjy-iy; dijz  =       jjz-iz
            dijx1 = abs(dijx-1); dijy1 = abs(dijy-1); dijz1 = abs(dijz-1)
            dijx  = abs(dijx  ); dijy  = abs(dijy  ); dijz  = abs(dijz  )
            grdreactmp = grdreactmp + &
            surf_crg(jsrf,jatm)*( gcrg(1,iatm)*l_green(dijx ,dijy ,dijz ) + gcrg(2,iatm)*l_green(dijx1,dijy ,dijz ) + &
                                  gcrg(3,iatm)*l_green(dijx ,dijy1,dijz ) + gcrg(4,iatm)*l_green(dijx1,dijy1,dijz ) + &
                                  gcrg(5,iatm)*l_green(dijx ,dijy ,dijz1) + gcrg(6,iatm)*l_green(dijx1,dijy ,dijz1) + &
                                  gcrg(7,iatm)*l_green(dijx ,dijy1,dijz1) + gcrg(8,iatm)*l_green(dijx1,dijy1,dijz1) ) 

         end do
         grdreac = grdreac + grdreactmp

      end do  !  jp = jfirst, jlast
      
   end do  !  iatm = 1, ilast
   
   grdreac = grdreac*AMBER_ELECTROSTATIC*HALF/h
!   print *, grdreac, 'FD reaction'
contains

function l_green (i,j,k)

   implicit none
   integer i,j,k
   _REAL_ l_green

   if ( i <= 20  .and. j <= 20 .and. k <= 20 ) then
      l_green = green(i,j,k)
   else
      l_green = ONE/sqrt(dble(i*i+j*j+k*k))
   end if

end function l_green

end subroutine pb_fdreaction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Direct coulombic energy and forces.
subroutine pb_direct_reaction( natom, atmfirst, atmlast, idecomp, drcreac, outflag, pbfrc )
   
   use decomp, only: decpair

   ! Common variables
   
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
   
   ! Passed variables
   
   integer natom, atmfirst, atmlast, idecomp
   integer outflag(natom)
   _REAL_ drcreac
   _REAL_ pbfrc(3, natom)
   
   ! Local variables
   
   integer iatm, jatm, jp, ilast, jfirst, jlast
   integer dijx, dijx0, dijx1, dijx2, dijx3, dijy, dijy0, dijy1, dijy2, dijy3, dijz, dijz0, dijz1, dijz2, dijz3
   integer ix, iy, iz, jx, jy, jz
   integer iix, iiy, iiz, jjx, jjy, jjz, i, j
   _REAL_ gci1, gci2, gci3, gci4, gci5, gci6, gci7, gci8
   _REAL_ gcj1, gcj2, gcj3, gcj4, gcj5, gcj6, gcj7, gcj8
   _REAL_ gcij(27)
   _REAL_ factor, factor1, decfac
   _REAL_ ffx, ffy, ffz
   _REAL_ dumx, dumy, dumz
   _REAL_ frc(3, natom)
   _REAL_ drcreactmp 
   _REAL_ xx_i, yy_i, zz_i, xx_j, yy_j, zz_j, dx(1:3)
   
   integer atmflag, crgflag, ncount
   integer isrf, jsrf, icrg, jcrg

   ! begin code
   
   factor  = ONE/( FOURPI*epsin*h )
   factor1 = HALF/( FOURPI*epsin*h*h )
   decfac  = -factor*frcfac
   
   drcreac = ZERO
   frc = ZERO

   ! first loop over iatm's grid charges with iatm's surface charges
   do iatm = atmfirst, atmlast
!     if ( ligand .and. realflag(iatm) == 0 ) then
!        cycle
!     else if ( multiblock .and. realflag(iatm) == 0 ) then
!        cycle ! this is only for debugging
!     else if ( multiblock .and. liveflag(iatm) == 0 ) then
!        cycle
!     else if ( outflag(iatm) == 1 ) then
!        cycle
!     end if
      xx_i = acrd(1,iatm); yy_i = acrd(2,iatm); zz_i = acrd(3,iatm)

      ! loop over iatm's grid charges over iatm's surface charges 

      drcreactmp = ZERO
      do isrf = 1, crg_num(iatm)
            dx(1) = pos_crg(1, isrf, iatm) - xx_i 
            dx(2) = pos_crg(2, isrf, iatm) - yy_i
            dx(3) = pos_crg(3, isrf, iatm) - zz_i
            drcreactmp = drcreactmp + surf_crg(isrf,iatm)*acrg(iatm)/sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
      end do
         drcreac = drcreac + drcreactmp

   end do
   ! second loop over iatm's grid charges with iatm's surface charges

   do iatm = atmfirst, atmlast
!     if ( ligand .and. realflag(iatm) == 0 ) then
!        cycle
!     else if ( multiblock .and. realflag(iatm) == 0 ) then
!        cycle ! this is only for debugging
!     else if ( multiblock .and. liveflag(iatm) == 0 ) then
!        cycle
!     else if ( outflag(iatm) == 1 ) then
!        cycle
!     end if
      xx_i = acrd(1,iatm); yy_i = acrd(2,iatm); zz_i = acrd(3,iatm)

      jfirst = iar1pb(4, iatm-1) + 1
      jlast  = iar1pb(2, iatm)
      dumx = ZERO; dumy = ZERO; dumz = ZERO
      do jp = jfirst, jlast
         jatm = iprshrt(jp)
         if( jatm > atmlast ) cycle
         xx_j = acrd(1,jatm); yy_j = acrd(2,jatm); zz_j = acrd(3,jatm)
        
         ! loop over jatm's grid charges over iatm's surface charges 

         drcreactmp = ZERO
         do isrf = 1, crg_num(iatm)
               dx(1) = pos_crg(1, isrf, iatm) - xx_j 
               dx(2) = pos_crg(2, isrf, iatm) - yy_j
               dx(3) = pos_crg(3, isrf, iatm) - zz_j
            drcreactmp = drcreactmp + surf_crg(isrf,iatm)*acrg(jatm)/sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         end do
         drcreac = drcreac + drcreactmp

         ! loop over iatm's grid charges over jatm's surface charges 

         drcreactmp = ZERO
         do jsrf = 1, crg_num(jatm)
               dx(1) = pos_crg(1, jsrf, jatm) - xx_i 
               dx(2) = pos_crg(2, jsrf, jatm) - yy_i
               dx(3) = pos_crg(3, jsrf, jatm) - zz_i
            drcreactmp = drcreactmp + surf_crg(jsrf,jatm)*acrg(iatm)/sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
         end do
         drcreac = drcreac + drcreactmp

      end do  !  jp = jfirst, jlast
      
   end do  !  iatm = 1, ilast
   
   drcreac = drcreac*AMBER_ELECTROSTATIC*HALF
!  print *, drcreac, 'direct reaction'
end subroutine pb_direct_reaction
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric boundary energy and forces
subroutine pb_dbene( verbose,eneout,natom,atmlast,ifcap,npdec,idecomp,m,irespw,ipres, &
                     eel,f,insas,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd
   use decomp, only: decpair
 
   ! Passed variables
 
   logical verbose, eneout
   integer natom,atmlast,ifcap,npdec,idecomp,m,irespw(*),ipres(*)
   integer insas(0:xm+1,0:ym+1,0:zm+1)
   _REAL_ phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)
   _REAL_ eel, f(3,natom)
 
   ! Local variables
 
   integer  i, j, k, iatm, jatm, matm, natm, iarc, ip
!  character*10 str
   _REAL_ srfcrg, factor, scalfact
   _REAL_ g(3), x(3), dx(3), crd(3), dx2, dist, rdist, acg
   _REAL_ eelrf
   _REAL_ d2inv, dinv, de, dff, df(3), dum(3)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ fx(natom), fy(natom), fz(natom)
   _REAL_ pcolomb

   _REAL_, parameter :: smallcrg = 0.5d0

   ! initialization
   factor = THREE*h/(TWOPI)
   scalfact = ONE
 
   ax = acrd(1,1:natom)
   ay = acrd(2,1:natom)
   az = acrd(3,1:natom)

   srfcrg = ZERO; eel = ZERO
   fx = ZERO; fy = ZERO; fz = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO

   ! calculate the total potential inside solute for singular-free PB (bcopt = 6/7)

   if ( bcopt == 6 .or. bcopt == 7 ) then
      do k = 1, zm; do j = 1, ym; do i = 1, xm
         if ( insas(i,j,k) > 0 ) phi(i,j,k) = phi(i,j,k) + cphi(i,j,k)
      end do; end do; end do
   end if

   ! for InsightII display
   !open (unit=55, file='ms.dot')
   !write (55, '("DOTS")')
 
   do ip = 1, nbnd
      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k

      ! project the surface grid point on to the molecular surface, crd() is the new coord

      if      ( iatm > 0 ) then
         ! the contact boundary grid points are projected to the atom spheres
         x(1:3) = acrd(1:3,iatm)
         dist = radi(iatm)
      else if ( iatm < 0 ) then
         ! the reentry boundary grid points are projected to the solvent spheres
         x(1:3) = arccrd(1:3,-iatm)
         dist = dprob
      else
         ! otherwise do not project. The setting should make crd = g
         x(1:3) = g(1:3)
         dist = ONE
      end if

      dx = g - x; dx2 = dx(1)**2 + dx(2)**2 + dx(3)**2
      if ( dx2 == ZERO ) then
         rdist = ONE
      else
         rdist = dist*ONE/sqrt(dx2)
      end if
      crd = x + dx*rdist

!if ( (.not.ligand).and.&
!    g(1)>-14.0.and.g(1)<12.00.and.&
!    g(2)>-10.5.and.g(2)< 9.50.and.&
!    g(3)>-12.5.and.g(3)< 3.50) then
!write(str,'(i10)') ip
!str=adjustl(str)
!write(3000,'(2x,a,a6,f9.5,2f15.5)') "H",str,crd(1:3)
!else if ( ligand ) then
!write(str,'(i10)') ip
!str=adjustl(str)
!write(3000,'(2x,a,a6,f9.5,2f15.5)') "H",str,crd(1:3)
!if ( ip == 2982 ) then
!  write(6,*)"mjhsieh: ",insas(i,j,k),i,j,k,g
!  write(6,*)"mjhsieh: ",iepsav(4,ip),arccrd(1:3,iepsav(4,ip))
! write(6,*)"mjhsieh: ",dot_product(crd(1:3)-arccrd(1:3,iepsav(4,ip)),crd(1:3)-arccrd(1:3,iepsav(4,ip)))
! write(6,*)"mjhsieh: ",g
! iarc = dotarc(iepsav(4,ip)) ! this leads to the circle/arc
! write(6,*) 'the two atoms forming the arc', iarc, arcatm(1:2,iarc)
!end if
!end if

      ! for InsightII display
      !write (55,'(4(f8.3,2x))') crd(1:3), 300.
 
      ! compute induced charge on the molecular surface

      acg = factor*(phi(i,j,k)-&
      SIXTH*( phi(i-1,j,k)+phi(i+1,j,k)+phi(i,j-1,k)+phi(i,j+1,k)+phi(i,j,k-1)+phi(i,j,k+1) ))

      ! if only the reaction field potential is computed (bcopt = 8/9), there is no atomic grid
      ! charge in acg.
      ! otherwise, atomic grid charges should be removed from acg.

      if ( bcopt < 8 ) acg = acg - chgrd(i,j,k)
      srfcrg = srfcrg + acg*frcfac*INV_AMBER_ELECTROSTATIC2
 
      ! compute reaction field energy and forces due to this boundary grid point

      eelrf = ZERO
      dum = ZERO
      if(idecomp < 3) then
         do jatm = 1, atmlast
            if((ifcap == 2 .or. ifcap == 5) .and. outflag(jatm) == 1) cycle
            if ( ligand .or. multiblock ) then
               if (liveflag(jatm) == 0) cycle
            endif
            !if(ligand .and. liveflag(jatm) == 0) cycle
            !if(multiblock .and. liveflag(jatm) == 0) cycle

            dx(1) = crd(1) - ax(jatm)
            dx(2) = crd(2) - ay(jatm)
            dx(3) = crd(3) - az(jatm)
            dinv = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2); d2inv = dinv**2

            de = acg*acrg(jatm)*dinv
            eelrf = eelrf + de

            !-- PB decomp
            if(idecomp == 1 .or. idecomp == 2) call decpair(1,jatm,jatm,de*frcfac*HALF)

            if ( frcopt == 1 ) then
               dff = de*d2inv; df(1) = dx(1)*dff; df(2) = dx(2)*dff; df(3) = dx(3)*dff
               fx(jatm) = fx(jatm) + df(1)
               fy(jatm) = fy(jatm) + df(2)
               fz(jatm) = fz(jatm) + df(3)
               dum(1) = dum(1) - df(1)
               dum(2) = dum(2) - df(2)
               dum(3) = dum(3) - df(3)
            end if
         end do
      !-- PB pairwise decomp
      else if(idecomp >= 3) then
         do n = 1, npdec
            do jatm = ipres(irespw(n)), ipres(irespw(n)+1)-1
               if ( ligand .or. multiblock ) then
                  if (liveflag(jatm) == 0) cycle
               endif
               !if(ligand .and. liveflag(jatm) == 0) cycle
               !if(multiblock .and. liveflag(jatm) == 0) cycle
               dx(1) = crd(1) - ax(jatm)
               dx(2) = crd(2) - ay(jatm)
               dx(3) = crd(3) - az(jatm)
               dinv = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
               de = acg*acrg(jatm)*dinv
               eelrf = eelrf + de
               call decpair(-1,ipres(irespw(m)),ipres(irespw(n)),de*frcfac*HALF)
            end do
         end do
      end if

      ! collecting energy

      eel = eel + eelrf
 
      if (sasopt /= 2 .and. frcopt == 1) call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)
   end do ! end of ip = 1, nbnd

   ! for InsightII display
   !close(55)
   !stop

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz
   fx = fx+dbx; fy = fy+dby; fz = fz+dbz
   if ( scalerf .and. abs(totcrg) > smallcrg ) then
      scalfact = abs( totcrg/srfcrg*(ONE/epsin - ONE/epsout)*eps0 )
      srfcrg = scalfact*srfcrg
      eel = scalfact*eel
      fx = scalfact*fx
      fy = scalfact*fy
      fz = scalfact*fz
   end if
    
   eel = HALF*frcfac*eel
   do iatm = 1, atmlast
      if((ifcap == 2 .or. ifcap == 5) .and. outflag(iatm) == 1) cycle
      f(1,iatm) = f(1,iatm) - frcfac*fx(iatm)
      f(2,iatm) = f(2,iatm) - frcfac*fy(iatm)
      f(3,iatm) = f(3,iatm) - frcfac*fz(iatm)
   end do
    
   if ( eneout ) then
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eel
   end if


end subroutine pb_dbene
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_dbfrc_fld(verbose,eneout,natom,f,epsx,epsy,epsz,phi,cphi)
  
   use solvent_accessibility, only : dprob, radi, arccrd

   implicit none

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ f(3,natom)
   _REAL_ epsx(0:xm,1:ym,1:zm)
   _REAL_ epsy(1:xm,0:ym,1:zm)
   _REAL_ epsz(1:xm,1:ym,0:zm)
   _REAL_ phi(xm,ym,zm)
   _REAL_ cphi(xm,ym,zm)

   ! local variables

   integer i, j, k, iatm
   integer xp, yp, zp
   _REAL_ df, factor, factor1, sfactor
   _REAL_ lEx, lEy, lEz
   _REAL_ dx(3), dum(3)
   _REAL_ dfx, dfy, dfz 
   _REAL_ adx, ady, adz
   _REAL_ x(3), crd(3)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_ rsphere, cut, r, r2, h2

   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO

   ! initialization

   h2 = h*h
   sfactor = -HALF*(epsout-epsin)/(epsin*epsout)

   ! contributions from epsx boundary edges

!  write (58, '("contact")')
!  write (59, '("reentry")')
   do xp = 1, nbndx
      i = iepsavx(1,xp); j = iepsavx(2,xp); k = iepsavx(3,xp); iatm = iepsavx(4,xp)
      lEx = phi(i+1,j,k)+cphi(i+1,j,k) - phi(i,j,k)-cphi(i,j,k)
      crd(1) = gox + h*(i+fedgex(xp)); crd(2) = goy + h*j; crd(3) = goz + h*k
      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         dx = crd - x
!        write(58,'(a,3(f20.15),2x)') "H", crd(1:3)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
      else
         x(1:3) = arccrd(1:3,-iatm)
         dx = x - crd
!        write(59,'(a,3(f20.15),2x)') "H", crd(1:3)
         rsphere = dprob
      end if
      adx = abs(dx(1))

      df = sfactor*lEx*lEx*epsx(i,j,k)**2
      factor1 = df/adx
      dum = factor1*dx

      ref(1)=ref(1)+dum(1)
      ref(2)=ref(2)+dum(2)
      ref(3)=ref(3)+dum(3)
!     df = sfactor*lEx*lEx*epsx(i,j,k)**2
!     cut = gox+h*(i+HALF)-x(1)
!     r = ONE/sqrt(cut**2+dx(2)**2+dx(3)**2)
!     r2 = r*r
!     factor1 = df*(rsphere*r)**3*abs(cut)/adx**2
!     factor1 = factor1*(ONE+EIGHTH*h2*r2*(THREE-FIVE*cut**2*r2))
!     dum = factor1*dx

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)
!     write(31,*) i,j,k,epsx(i,j,k)*lEx/h*FOURPI

   enddo

   ! contributions from epsy boundary grid edges
             
   do yp = 1, nbndy
      i = iepsavy(1,yp); j = iepsavy(2,yp); k = iepsavy(3,yp); iatm = iepsavy(4,yp)
      lEy = phi(i,j+1,k)+cphi(i,j+1,k) - phi(i,j,k)-cphi(i,j,k)
      crd(1) = gox + h*i; crd(2) = goy + h*(j+fedgey(yp)); crd(3) = goz + h*k

      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         dx = crd - x
!        write(58,'(a,3(f20.15),2x)') "H", crd(1:3)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
      else
         x(1:3) = arccrd(1:3,-iatm)
         dx = x - crd
!        write(59,'(a,3(f20.15),2x)') "H", crd(1:3)
         rsphere = dprob
      end if
      ady = abs(dx(2))

      df = sfactor*lEy*lEy*epsy(i,j,k)**2
      factor1= df/ady
      dum = factor1*dx

      ref(1)=ref(1)+dum(1)
      ref(2)=ref(2)+dum(2)
      ref(3)=ref(3)+dum(3)
!     df = sfactor*lEy*lEy*epsy(i,j,k)**2
!     cut = goy+h*(j+HALF)-x(2)
!     r = ONE/sqrt(cut**2+dx(1)**2+dx(3)**2)
!     r2 = r*r
!     factor1 = df*(rsphere*r)**3*abs(cut)/ady**2
!     factor1 = factor1*(ONE+EIGHTH*h2*r2*(THREE-FIVE*cut**2*r2))
!     dum = factor1*dx

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)

   enddo
           
   ! contributions from epsz boundary grid edges
           
   do zp = 1, nbndz
      i = iepsavz(1,zp); j = iepsavz(2,zp); k = iepsavz(3,zp); iatm = iepsavz(4,zp)
      lEz = phi(i,j,k+1)+cphi(i,j,k+1) - phi(i,j,k)-cphi(i,j,k)
      crd(1) = gox + h*i; crd(2) = goy + h*j; crd(3) = goz + h*(k+fedgez(zp))
      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         dx = crd - x
!        write(58,'(a,3(f20.15),2x)') "H", crd(1:3)
         if ( sasopt > 0 ) then
            rsphere = radi(iatm)+dprob
         else
            rsphere = radi(iatm)
         end if
      else
         x(1:3) = arccrd(1:3,-iatm)
         dx = x - crd
!        write(59,'(a,3(f20.15),2x)') "H", crd(1:3)
         rsphere = dprob
      end if
      adz = abs(dx(3))

      df = sfactor*lEz*lEz*epsz(i,j,k)**2
      factor1= df/adz
      dum = factor1*dx

      ref(1)=ref(1)+dum(1)
      ref(2)=ref(2)+dum(2)
      ref(3)=ref(3)+dum(3)
!     df = sfactor*lEz*lEz*epsz(i,j,k)**2
!     cut = goz+h*(k+HALF)-x(3)
!     r = ONE/sqrt(cut**2+dx(1)**2+dx(2)**2)
!     r2 = r*r
!     factor1 = df*(rsphere*r)**3*abs(cut)/adz**2
!     factor1 = factor1*(ONE+EIGHTH*h2*r2*(THREE-FIVE*cut**2*r2))
!     dum = factor1*dx

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)

   end do

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz

!  write(102,*) ' :::: Atomic contact forces ::::'
!  do iatm = 1, natom
!     write(102,'(3e20.6)') cnx(iatm)*frcfac,cny(iatm)*frcfac,cnz(iatm)*frcfac
!  end do

!  write(102,*) ' :::: Atomic reentry forces ::::'
!  do iatm = 1, natom
!     write(102,'(3e20.6)') rnx(iatm)*frcfac,rny(iatm)*frcfac,rnz(iatm)*frcfac
!  end do
   ref=0.0d0 !XP: calculate DB force
   write(102,*) ' :::: Atomic DB forces ::::'
   do iatm = 1, natom
      write(102,'(3e20.6)') dbx(iatm)*frcfac,dby(iatm)*frcfac,dbz(iatm)*frcfac
   end do
   ref=ref*frcfac
   write(6,*) 'Vector DB forces is',ref(1:3)
   write(6,*) 'Total DB forces is',sqrt(ref(1)**2+ref(2)**2+ref(3)**2)

   f(1,1:natom) = f(1,1:natom) + dbx(1:natom)
   f(2,1:natom) = f(2,1:natom) + dby(1:natom)
   f(3,1:natom) = f(3,1:natom) + dbz(1:natom)

end subroutine pb_dbfrc_fld
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_dbfrc_crg( verbose,eneout,natom,eelrf,f,epsx,epsy,epsz,zv,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd

   implicit none

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ eelrf
   _REAL_ f(3,natom)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)
   _REAL_ zv(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)

   ! local variables

   integer, parameter :: n_point = 27
   integer i, j, k, iatm, ip
   _REAL_, parameter :: small_dn = 1.d-3 * AMBER_ELECTROSTATIC
   _REAL_ g(3), x(3), dx(3), crd(3), rn(3), sgn, rdx, rsphere, crg
   _REAL_ srfcrg, coulomb(3)
   _REAL_ fx0, fy0, fz0, rh
   _REAL_ up, dudxi0, dudyi0, dudzi0
   _REAL_ dn, dt2, dum(3)
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_ Ex, Ey, Ez, w1, w2, w3, cnnt, factor, reps0
   integer sgnx, sgny, sgnz
   _REAL_, allocatable :: pol_charge(:)

   allocate (pol_charge(1:nbnd))

   ! initialization

   rh = ONE/h; reps0 = ONE/eps0
   factor = FOURPI*eps0*AMBER_ELECTROSTATIC
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO
   x = ZERO; rsphere = ZERO; sgn = ONE; coulomb = ZERO

   ! compute polarization charges on the boundary grid points and
   ! report total in srfcrg
   ! pol_charge is in kcal/mol

   srfcrg = ZERO
   call get_charge_pol(nbnd,phi,cphi,chgrd,pol_charge,srfcrg)

   ! main double loops over polarization charges and atom charges

   eelrf = ZERO
   if( isurfchg > 0 ) open(131,file='srfchg.pos',form="formatted")
   do ip = 1, nbnd

      ! collect boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
      crg = pol_charge(ip)

      ! project the boundary grid point onto the molecular surface, crd() is the
      ! projected coord, x() is the atom/probe center coord, and fx/y/z0 is the grid
      ! version of crd()

      if ( iatm == 0 ) then
         write(6,*) 'PB Bomb in pb_dbfrc_fld(): can not find projection atom/probe'
         call mexit(6, 1)
      end if

      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         sgn = ONE
      else if ( iatm < 0 ) then
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         sgn = -ONE
      end if
      dx = g - x
      rdx = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
      crd = x + dx*rdx*rsphere

      !  write surface charges to output
      if (isurfchg > 0 ) then
        write(131,'(3f15.5,f10.6)') crd(1:3),crg*INV_AMBER_ELECTROSTATIC
      endif

      ! normal vector, it should be sgn*dx*rdx

      rn = sgn*dx*rdx

      ! compute reaction field energy, forces, and coulomb field
      ! between the current surface charge and natom atomic charges
      ! inner loop over atoms...

      call get_coulomb(natom,crg,crd,acrg,ax,ay,az,qex,qey,qez,eelrf,coulomb)

      ! now it is time to compute dbf
      dn = ZERO
      if ( iatm > 0 ) then

         ! interpolate E

         ! grid unit of the projected boundary point

         fx0 = (crd(1)-gox)*rh
         fy0 = (crd(2)-goy)*rh
         fz0 = (crd(3)-goz)*rh

         ! compute E on the inner side of surface position crd(1:3)

         call gradu(xm,ym,zm,-ONE,n_point,10,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,zv)

         dudxi0 = -dudxi0*factor
         dudyi0 = -dudyi0*factor
         dudzi0 = -dudzi0*factor
       
         ! add the coulomb field to get the total E of inner side  
         ! converted to electric displacement

         dudxi0 = (dudxi0 + coulomb(1))*epsin*reps0
         dudyi0 = (dudyi0 + coulomb(2))*epsin*reps0
         dudzi0 = (dudzi0 + coulomb(3))*epsin*reps0

         ! compute the DBF as HALF*Q*D^2/Dn
         ! when denominator dn is tiny, use the limiting law.

         dt2 = dudxi0*dudxi0 + dudyi0*dudyi0 + dudzi0*dudzi0
         dn = dudxi0*rn(1)+dudyi0*rn(2)+dudzi0*rn(3)

      elseif ( iatm < 0 ) then

         ! weighted sum of D
         ! bcopt = 9 is better?

         cnnt = ZERO

         !  contact
         !  consider inside, positive rn
         !  should be 1, but inside -2, outside 2
         !  reentry
         !  consider inside, positive rn
         !  should be -1, but there is sgn in rn, and inside -1, outside 1 
         sgnx = -1; sgny = -1; sgnz = -1

         Ex = ZERO; Ey = ZERO; Ez = ZERO
         w1 = ZERO; w2 = ZERO; w3 = ZERO
         if ( rn(1) < ZERO ) sgnx = -sgnx
         if ( rn(2) < ZERO ) sgny = -sgny
         if ( rn(3) < ZERO ) sgnz = -sgnz
         if ( zv(i,j,k) < ZERO ) then
            sgnx = -sgnx
            sgny = -sgny
            sgnz = -sgnz
         end if
         if ( zv(i,j,k)*zv(i+sgnx,j,k) < ZERO .and. abs(rn(1)) > 1.d-2 ) then
            !  Ex = dble(sign(1,sgnx))*(phi(i,j,k)-phi(i+sgnx,j,k))
            Ex = dble(sign(1,sgnx))*(phi(i,j,k)+cphi(i,j,k)-phi(i+sgnx,j,k)-cphi(i+sgnx,j,k))
            w1 = ONE/rn(1)
            cnnt = cnnt + ONE
         end if
         if ( zv(i,j,k)*zv(i,j+sgny,k) < ZERO .and. abs(rn(2)) > 1.d-2 ) then
            !  Ey = dble(sign(1,sgny))*(phi(i,j,k)-phi(i,j+sgny,k))
            Ey = dble(sign(1,sgny))*(phi(i,j,k)+cphi(i,j,k)-phi(i,j+sgny,k)-cphi(i,j+sgny,k))
            w2 = ONE/rn(2)
            cnnt = cnnt + ONE
         end if
         if ( zv(i,j,k)*zv(i,j,k+sgnz) < ZERO .and. abs(rn(3)) > 1.d-2 ) then
            !  Ez = dble(sign(1,sgnz))*(phi(i,j,k)-phi(i,j,k+sgnz))
            Ez = dble(sign(1,sgnz))*(phi(i,j,k)+cphi(i,j,k)-phi(i,j,k+sgnz)-cphi(i,j,k+sgnz))
            w3 = ONE/rn(3)
            cnnt = cnnt + ONE
         end if
         if ( cnnt > ZERO ) then
            if ( sgnx == 1 ) sgnx = 0
            if ( sgny == 1 ) sgny = 0
            if ( sgnz == 1 ) sgnz = 0
            Ex = Ex*epsx(i+sgnx,j,k)*reps0
            Ey = Ey*epsy(i,j+sgny,k)*reps0
            Ez = Ez*epsz(i,j,k+sgnz)*reps0
            dn = (Ex*w1+Ey*w2+Ez*w3)/cnnt*rh*factor
         else 
!           write(6,*) "PB warning: cnnt = 0, at", i, j, k
            continue
         end if

      end if

!     if ( abs(dn) > small_dn ) then
!        dum = HALF*crg/dn*dt2*rn
!     else
!        dum = ZERO
!        dum = HALF*INV_FOURPI*(epsin/epsout-ONE)*h*h*dt2*dx
!     end if
      dum = HALF*crg*dn*rn

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)
   end do

   deallocate (pol_charge)

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz

   open (unit = 103, file = 'force.dat')

   write(103,*) ' :::: Atomic qE forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') qex(iatm),qey(iatm),qez(iatm)
   end do

   write(103,*) ' :::: Atomic contact forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') cnx(iatm),cny(iatm),cnz(iatm)
   end do

   write(103,*) ' :::: Atomic reentry forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') rnx(iatm),rny(iatm),rnz(iatm)
   end do

   write(103,*) ' :::: Atomic DB forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') dbx(iatm),dby(iatm),dbz(iatm)
   end do

   f(1,1:natom) = f(1,1:natom) + dbx(1:natom) + qex(1:natom)
   f(2,1:natom) = f(2,1:natom) + dby(1:natom) + qey(1:natom)
   f(3,1:natom) = f(3,1:natom) + dbz(1:natom) + qez(1:natom)

   eelrf = HALF*eelrf
   !if ( eneout ) then
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg*INV_AMBER_ELECTROSTATIC
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eelrf
      write(6, '(1x,a,3f12.4)') 'Total solvation force',&
                                sum(dbx(1:natom))+sum(qex(1:natom)),&
                                sum(dby(1:natom))+sum(qey(1:natom)),&
                                sum(dbz(1:natom))+sum(qez(1:natom))
   !end if

end subroutine pb_dbfrc_crg
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric boundary forces, second field-based strategy
!+ loop over boundary grid edges, three loops
subroutine pb_dbfrc_fld2( verbose,eneout,natom,eelrf,f,epsx,epsy,epsz,zv,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd

   ! Passed variables

   logical verbose, eneout
   integer natom
   _REAL_ zv(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)
   _REAL_ eelrf, f(3,natom)

   ! Local variables

   integer i, j, k, iatm, ip
   _REAL_ x(3), crd(3)
   _REAL_ fbnd, dum(3)
   _REAL_ fx0, fy0, fz0
   _REAL_ rn(1:3), rsphere, dr
   _REAL_ ds, srfarea
   _REAL_ crg, h2, hh, r1, r2, r3
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_ up, dudxi0, dudyi0, dudzi0
   _REAL_ dudxo0, dudyo0, dudzo0
   _REAL_ dudni, dudno, coulomb(3) 
   _REAL_ E2

   ! mjhsieh: warning eliminator
   rsphere = -1d0
   ! initialization for DBF

   h2 = h*h; hh = HALF*h
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO
   coulomb = ZERO

   eelrf = ZERO
   srfarea = ZERO
   do ip = 1, nbndx
      i = iepsavx(1,ip); j = iepsavx(2,ip); k = iepsavx(3,ip); iatm = iepsavx(4,ip)
      fx0 = i + fedgex(ip); fy0 = j; fz0 = k
!     crd(1) = gox + h*i + fedgex(ip)*h; crd(2) = goy + h*j; crd(3) = goz + h*k
      crd(1) = gox + h*i + hh; crd(2) = goy + h*j; crd(3) = goz + h*k

      if ( iatm == 0 ) then
         x(1:3) = ZERO
         write(6,*) 'PBMD FATAL ERROR: can not find projection atom/probe' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         rn = crd - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn = x - crd
      end if

      dr = abs(rn(1))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1

      ! surface element
      ds = dr*rsphere**2*r3*h2
!     ds = ds*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(1)**2*r2))
      srfarea = srfarea + ds

#     include "pb_dbfrc_fld2.h"

   end do !nbndx

   do ip = 1, nbndy
      i = iepsavy(1,ip); j = iepsavy(2,ip); k = iepsavy(3,ip); iatm = iepsavy(4,ip)
      fx0 = i; fy0 = j + fedgey(ip); fz0 = k
!     crd(1) = gox + h*i; crd(2) = goy + h*j + fedgey(ip)*h; crd(3) = goz + h*k
      crd(1) = gox + h*i; crd(2) = goy + h*j + hh; crd(3) = goz + h*k

      if ( iatm == 0 ) then
         x(1:3) = ZERO
         write(6,*) 'PBMD FATAL ERROR: can not find projection atom/probe' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         rn = crd - x
      else
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn = x - crd
      end if

      dr = abs(rn(2))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1

      ! surface element
      ds = dr*rsphere**2*r3*h2
!     ds = ds*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(2)**2*r2))
      srfarea = srfarea + ds

#     include "pb_dbfrc_fld2.h"

   end do !nbndy

   do ip = 1, nbndz
      i = iepsavz(1,ip); j = iepsavz(2,ip); k = iepsavz(3,ip); iatm = iepsavz(4,ip)
      fx0 = i; fy0 = j; fz0 = k + fedgez(ip)
!     crd(1) = gox + h*i ; crd(2) = goy + h*j; crd(3) = goz + h*k+fedgez(ip)*h
      crd(1) = gox + h*i ; crd(2) = goy + h*j; crd(3) = goz + h*k+hh

      if ( iatm == 0 ) then
         x(1:3) = ZERO
         write(6,*) 'PBMD FATAL ERROR: can not find projection atom/probe' 
         call mexit(6, 1)
      elseif ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         rn = crd - x
      else 
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         rn = x - crd
      end if

      dr = abs(rn(3))
      r2 = ONE/(rn(1)**2+rn(2)**2+rn(3)**2)
      r1 = sqrt(r2)
      r3 = r2*r1

      ! surface element
      ds = dr*rsphere**2*r3*h2
!     ds = ds*(ONE+EIGHTH*h2*r2*(THREE-FIVE*rn(3)**2*r2))
      srfarea = srfarea + ds

#     include "pb_dbfrc_fld2.h"

   end do !nbndz

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz

   open (unit = 103, file = 'force.dat')

!  write(103,*) ' :::: Atomic qE forces ::::'
!  do iatm = 1, natom
!     write(103,'(3e20.6)') qex(iatm),qey(iatm),qez(iatm)
!  end do

!  write(103,*) ' :::: Atomic contact forces ::::'
!  do iatm = 1, natom
!     write(103,'(3e20.6)') cnx(iatm),cny(iatm),cnz(iatm)
!  end do

!  write(103,*) ' :::: Atomic reentry forces ::::'
!  do iatm = 1, natom
!     write(103,'(3e20.6)') rnx(iatm),rny(iatm),rnz(iatm)
!  end do

   write(103,*) ' :::: Atomic DB forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') dbx(iatm),dby(iatm),dbz(iatm)
   end do

   f(1,1:natom) = f(1,1:natom) + dbx(1:natom) + qex(1:natom)
   f(2,1:natom) = f(2,1:natom) + dby(1:natom) + qey(1:natom)
   f(3,1:natom) = f(3,1:natom) + dbz(1:natom) + qez(1:natom)

   !if ( eneout ) then
      write(6, '(1x,a,f12.4)') 'Total molecular surface', srfarea
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eelrf
      write(6, '(1x,a,3f12.4)') 'Total solvation force',&
                                sum(dbx(1:natom))+sum(qex(1:natom)),&
                                sum(dby(1:natom))+sum(qey(1:natom)),&
                                sum(dbz(1:natom))+sum(qez(1:natom))
   !end if

end subroutine pb_dbfrc_fld2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric boundary forces, second charge-based strategy
subroutine pb_dbfrc_crg2( verbose,eneout,natom,eelrf,f,epsx,epsy,epsz,zv,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd
   implicit none

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ eelrf
   _REAL_ f(3,natom)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)
   _REAL_ zv(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)

   ! local variables

   integer, parameter :: n_point = 7
   integer i, j, k, iatm, jatm, matm, natm, iarc, ip
   _REAL_ srfcrg, crg
   _REAL_ g(3), x(3), dx(3), crd(3), rn(3), rdx, rsphere, sgn
   _REAL_ dum(3), coulomb(3)
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_ fx0, fy0, fz0, rh
   _REAL_ up,dudxi0, dudyi0, dudzi0, dudni, E2
   _REAL_, allocatable :: pol_charge(:)

   allocate (pol_charge(1:nbnd))

   ! mjhsieh: warning eliminator
   sgn = ZERO

   ! initialization

   rh = ONE/h
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO
   x = ZERO; rsphere = ZERO; coulomb = ZERO

   ! compute polarization charges on the boundary grid points and
   ! report total in srfcrg
   ! pol_charge is in kcal/mol

   srfcrg = ZERO
   call get_charge_pol(nbnd,phi,cphi,chgrd,pol_charge,srfcrg)

   ! main double loops over polarization charges and atom charges

   eelrf = ZERO
   do ip = 1, nbnd

      ! collect boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
      crg = pol_charge(ip)

      ! project the surface grid point on to the molecular surface, crd() is the
      ! new coord, and x() is the atom/probe coord, dx() is the projection
      ! direction vector, fx/y/z0 is the grid version of crd()

      if ( iatm == 0 ) then
         write(6,*) 'PB Bomb in pb_dbfrc_crg2(): can not find projection atom/probe'
         call mexit(6, 1)
      end if

      ! collect data

      if ( iatm > 0 ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         sgn = ONE
      else if ( iatm < 0 ) then
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         sgn = -ONE
      end if

      ! project onto the sphere

      dx = g - x
      rdx = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
      crd = x + dx*rdx*rsphere

      ! normal direction vector

      rn = sgn*dx*rdx

      ! grid unit of the projected boundary point

      fx0 = (crd(1)-gox)*rh
      fy0 = (crd(2)-goy)*rh
      fz0 = (crd(3)-goz)*rh

      ! compute reaction field energy, forces, and coulomb field
      ! between the current surface charge and natom atomic charges
      ! inner loop over atoms...

      call get_coulomb(natom,crg,crd,acrg,ax,ay,az,qex,qey,qez,eelrf,coulomb)

      ! now it is time to compute dbf
      ! compute E on the inner side of surface position crd(1:3)

      call gradu(xm,ym,zm,-ONE,n_point,4,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,zv)

      dudxi0 = -dudxi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudyi0 = -dudyi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudzi0 = -dudzi0*FOURPI*eps0*AMBER_ELECTROSTATIC
    
      ! add the coulomb field to get the total E of inner side
      ! converted to electric displacement

      dudxi0 = (dudxi0 + coulomb(1))*epsin/eps0
      dudyi0 = (dudyi0 + coulomb(2))*epsin/eps0
      dudzi0 = (dudzi0 + coulomb(3))*epsin/eps0

      !dt2 = dudxi0*dudxi0 + dudyi0*dudyi0 + dudzi0*dudzi0
      dudni = dudxi0*rn(1)+dudyi0*rn(2)+dudzi0*rn(3)

      !dudno = dudni
      !dudxo0 = (dudxi0 - dudni*rn(1))/epsin*epsout + dudno*rn(1)
      !dudyo0 = (dudyi0 - dudni*rn(2))/epsin*epsout + dudno*rn(2)
      !dudzo0 = (dudzi0 - dudni*rn(3))/epsin*epsout + dudno*rn(3)

      ! part b: apply the normal field approximation
      !         or use the total field

      E2 = dudni
      !E2 = dudxi0*dudxo0 + dudyi0*dudyo0 + dudzi0*dudzo0

      ! compute the DBF as HALF*Q*Dn with the normal field approximation
      ! compute the DBF as HALF*Q*Di*Do/Dn with the total field

      dum = HALF*crg*E2*rn
      !if ( abs(dn) > small_dn ) then
      !   dum = HALF*crg/dn*E2*rn
      !else
      !   dum = HALF*crg*E2*rn
      !end if

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)

   end do

   deallocate (pol_charge)

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz

   open (unit = 103, file = 'force.dat')

   write(103,*) ' :::: Atomic qE forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') qex(iatm),qey(iatm),qez(iatm)
   end do

   write(103,*) ' :::: Atomic DB forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') dbx(iatm),dby(iatm),dbz(iatm)
   end do

   f(1,1:natom) = f(1,1:natom) + dbx(1:natom) + qex(1:natom)
   f(2,1:natom) = f(2,1:natom) + dby(1:natom) + qey(1:natom)
   f(3,1:natom) = f(3,1:natom) + dbz(1:natom) + qez(1:natom)

   eelrf = HALF*eelrf
   !if ( eneout ) then
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg*INV_AMBER_ELECTROSTATIC
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eelrf
      write(6, '(1x,a,3f12.4)') 'Total solvation force',&
                                sum(dbx(1:natom))+sum(qex(1:natom)),&
                                sum(dby(1:natom))+sum(qey(1:natom)),&
                                sum(dbz(1:natom))+sum(qez(1:natom))
   !end if

end subroutine pb_dbfrc_crg2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ dielectric boundary forces, third field-based strategy,
!+ loop over boundary grid points, one loop
subroutine pb_dbfrc_fld3( verbose,eneout,natom,eelrf,f,epsx,epsy,epsz,zv,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd

   implicit none

#  include "pb_constants.h"

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ eelrf
   _REAL_ f(3,natom)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)
   _REAL_ zv(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)

   ! local variables

   integer, parameter :: n_point = 7
   integer i, j, k, iatm, ip
   _REAL_ g(3), x(3), crd(3), dx(3), rn(3), sgn, rdx, rsphere
   _REAL_ fx0, fy0, fz0, rh
   _REAL_ ds, srfarea, epsth, repsin, repsout, repsp(3), repsm(3)
   _REAL_ E2, fbnd, dum(3)
   _REAL_ srfcrg, crg
   _REAL_ up, dudxi0, dudyi0, dudzi0
   _REAL_ dudxo0, dudyo0, dudzo0
   _REAL_ dudni, dudno, coulomb(3)
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)
   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_, allocatable :: pol_charge(:)

   allocate (pol_charge(1:nbnd))

   ! initialization

   rh = ONE/h
   repsin = ONE/epsin; repsout = ONE/epsout
   epsth = TWO/(repsin+repsout)
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO
   x = ZERO; sgn = ONE; rsphere = ZERO; crg = ZERO; coulomb = ZERO

   ! compute polarization charges on the boundary grid points and
   ! report total in srfcrg
   ! pol_charge is in kcal/mol

   srfcrg = ZERO
   call get_charge_pol(nbnd,phi,cphi,chgrd,pol_charge,srfcrg)

   eelrf = ZERO
   srfarea = ZERO
   do ip = 1, nbnd

      ! collect boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
      crg = pol_charge(ip)

      ! project the boundary grid point onto the molecular surface, crd() is the
      ! projected coord, x() is the atom/probe center coord, and fx/y/z0 is the grid
      ! version of crd()

      if ( iatm == 0 ) then
         write(6,*) 'PB Bomb in pb_dbfrc_fld(): can not find projection atom/probe'
         call mexit(6, 1)
      end if
      if      ( iatm > 0 ) then
         ! the contact boundary grid points are projected to the atom spheres
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         sgn = ONE
      else if ( iatm < 0 ) then
         ! the reentry boundary grid points are projected to the solvent spheres
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob
         sgn = -ONE
      else
         ! otherwise do not project. The setting should make crd = g
         x(1:3) = g(1:3)
         rsphere = ONE
         sgn = ONE
      end if
      dx = g - x; rdx = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
      crd = x + dx*rdx*rsphere

      ! normal direction vector

      rn = sgn*dx*rdx

      ! grid unit of the projected boundary point

      fx0 = (crd(1)-gox)*rh
      fy0 = (crd(2)-goy)*rh
      fz0 = (crd(3)-goz)*rh

      ! compute the effective surface element according to
      ! Cai, et al. JCC, in prep.
      ! first manipulate dielectric constants if smoothing is used

      if ( smoothopt == 0 ) then
         ds = dx(1)/epsx(i,j,k) - dx(1)/epsx(i-1,j,k) + &
              dx(2)/epsy(i,j,k) - dx(2)/epsy(i,j-1,k) + &
              dx(3)/epsz(i,j,k) - dx(3)/epsz(i,j,k-1)
      else
         repsp = repsin
         repsm = repsin
         if ( epsx(i-1,j,k) > epsth ) repsm(1) = repsout
         if ( epsx(i  ,j,k) > epsth ) repsp(1) = repsout
         if ( epsy(i,j-1,k) > epsth ) repsm(2) = repsout
         if ( epsy(i,j ,k ) > epsth ) repsp(2) = repsout
         if ( epsz(i,j,k-1) > epsth ) repsm(3) = repsout
         if ( epsz(i,j,k  ) > epsth ) repsp(3) = repsout
         ds = dx(1)*repsp(1) - dx(1)*repsm(1) + &
              dx(2)*repsp(2) - dx(2)*repsm(2) + &
              dx(3)*repsp(3) - dx(3)*repsm(3)
      end if
      ds = sgn*ds*rdx*h*h/(repsout-repsin)

      srfarea = srfarea + ds

      ! compute reaction field energy, forces, and coulomb field
      ! between the current surface charge and all atomic charges
      ! inner loop over "natom" natoms...

      call get_coulomb(natom,crg,crd,acrg,ax,ay,az,qex,qey,qez,eelrf,coulomb)

      ! now it is time to compute dbf
      ! part a: compute E_in and E_out ...
      !    compute E on the inner side of surface position crd(1:3)

      call gradu(xm,ym,zm,-ONE,n_point,4,fx0,fy0,fz0,up,dudxi0,dudyi0,dudzi0,phi,zv)

      dudxi0 = -dudxi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudyi0 = -dudyi0*FOURPI*eps0*AMBER_ELECTROSTATIC
      dudzi0 = -dudzi0*FOURPI*eps0*AMBER_ELECTROSTATIC

      !    add the coulomb field to get the total E on the inner side
      !    convert to displacement, D, for consistency with other methods

      dudxi0 = (dudxi0 + coulomb(1))*epsin/eps0
      dudyi0 = (dudyi0 + coulomb(2))*epsin/eps0
      dudzi0 = (dudzi0 + coulomb(3))*epsin/eps0

      !    get normal field component, which is continuous
      !    get D of the outer side based on the jump conditions

      dudni = dudxi0*rn(1)+dudyi0*rn(2)+dudzi0*rn(3)
      dudno = dudni
      !dudxo0 = (dudxi0 - dudni*rn(1))/epsin*epsout + dudno*rn(1)
      !dudyo0 = (dudyi0 - dudni*rn(2))/epsin*epsout + dudno*rn(2)
      !dudzo0 = (dudzi0 - dudni*rn(3))/epsin*epsout + dudno*rn(3)

      ! part b: apply the normal field approximation
      !         or use the total field

      E2 = dudni*dudni
      !E2 = dudxi0*dudxo0 + dudyi0*dudyo0 + dudzi0*dudzo0

      ! part c: compute the surface force element

      fbnd = INV_EIGHTPI*(eps0*(epsin-epsout)/(epsin*epsout))*E2*ds
      dum(1) = fbnd*rn(1)
      dum(2) = fbnd*rn(2)
      dum(3) = fbnd*rn(3)

      call dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)
   end do

   deallocate (pol_charge)

   dbx = cnx+rnx; dby = cny+rny; dbz = cnz+rnz

   open (unit = 103, file = 'force.dat')

   write(103,*) ' :::: Atomic qE forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') qex(iatm),qey(iatm),qez(iatm)
   end do

   write(103,*) ' :::: Atomic DB forces ::::'
   do iatm = 1, natom
      write(103,'(3e20.6)') dbx(iatm),dby(iatm),dbz(iatm)
   end do

   f(1,1:natom) = f(1,1:natom) + dbx(1:natom) + qex(1:natom)
   f(2,1:natom) = f(2,1:natom) + dby(1:natom) + qey(1:natom)
   f(3,1:natom) = f(3,1:natom) + dbz(1:natom) + qez(1:natom)

   eelrf = HALF*eelrf
   !if ( eneout ) then
      write(6, '(1x,a,f12.4)') 'Total molecular surface', srfarea
      write(6, '(1x,a,f12.4)') 'Total surface charge ', srfcrg*INV_AMBER_ELECTROSTATIC
      write(6, '(1x,a,f12.4)') 'Reaction field energy', eelrf
      write(6, '(1x,a,3f12.4)') 'Total solvation force',&
                                sum(dbx(1:natom))+sum(qex(1:natom)),&
                                sum(dby(1:natom))+sum(qey(1:natom)),&
                                sum(dbz(1:natom))+sum(qez(1:natom))
   !end if

end subroutine pb_dbfrc_fld3
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine pb_crgview( verbose,eneout,natom,f,epsx,epsy,epsz,insas,phi,chgrd,cphi )

   use solvent_accessibility, only : dprob, radi, arccrd, dotarc, arcatm, narcdot, ntri, triatm, triopt

   implicit none

   ! passed variables

   logical verbose, eneout
   integer natom
   _REAL_ f(3,natom)
   _REAL_ epsx(0:xm,1:ym,1:zm), epsy(1:xm,0:ym,1:zm), epsz(1:xm,1:ym,0:zm)
   _REAL_ insas(0:xm+1,0:ym+1,0:zm+1), phi(xm,ym,zm), chgrd(xm,ym,zm), cphi(xm,ym,zm)

   ! local variables

   integer, parameter :: n_point = 27
   integer i, j, k, iatm, jatm, matm, natm, iarc, ip
   _REAL_, parameter :: small_dn = 1.d-3 * AMBER_ELECTROSTATIC
   _REAL_ g(3), x(3), dx(3), crd(3), rn(3), sgn, rdx, rsphere, crg
   _REAL_ srfcrg, coulomb(3)
   _REAL_ fx0, fy0, fz0, rh
   _REAL_ up, dudxi0, dudyi0, dudzi0
   _REAL_ dn, dt2, dum(3)
   _REAL_ mvec(3), nvec(3), mdotn, mxnv(3), rmdotn2, fdotm, fdotn
   _REAL_ dfm, dfn, dum_norm(3), dum_tang(3), dumnorm
   _REAL_ ax(natom), ay(natom), az(natom)
   _REAL_ qex(natom), qey(natom), qez(natom)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)

   _REAL_ dbx(natom), dby(natom), dbz(natom)
   _REAL_, allocatable :: pol_charge(:)
   integer   atmflag, crgflag

   integer reentopt, m, n, mp, np, patm, qatm
  ! _REAL_, allocatable :: a(:,:), u(:,:), w(:), v(:,:), b(:), t(:)
   _REAL_ d1, d2, wmax, thresh, TOL
   _REAL_ pvec(3), qvec(3), mdist, ndist, pdist, qdist
   _REAL_ scalar
   mp = 3; np = 3;
   TOL = 1.d-5


   allocate (pol_charge(1:nbnd))
   ! initialization

   rh = ONE/h
   !  a = ZERO; u = ZERO; w = ZERO; v = ZERO; b = ZERO; t = ZERO
   ax = acrd(1,1:natom); ay = acrd(2,1:natom); az = acrd(3,1:natom)
   qex = ZERO; qey = ZERO; qez = ZERO
   cnx = ZERO; cny = ZERO; cnz = ZERO
   rnx = ZERO; rny = ZERO; rnz = ZERO
   dbx = ZERO; dby = ZERO; dbz = ZERO
   x = ZERO; rsphere = ZERO; sgn = ONE; coulomb = ZERO
   crg_num = ZERO

   ! compute polarization charges on the boundary grid points and
   ! report total in srfcrg
   ! pol_charge is in kcal/mol

   srfcrg = ZERO
   call get_charge_pol(nbnd,phi,cphi,chgrd,pol_charge,srfcrg)
   ! main double loops over polarization charges and atom charges

   do ip = 1, nbnd

      ! collect boundary grid point info ...

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip); iatm = iepsav(4,ip)
      g(1) = gox + h*i; g(2) = goy + h*j; g(3) = goz + h*k
      crg = pol_charge(ip)
!     write(400,'(3I7,f15.6)')i, j, k, crg

      ! project the boundary grid point onto the molecular surface, crd() is the
      ! projected coord, x() is the atom/probe center coord, and fx/y/z0 is the grid
      ! version of crd()

      if ( iatm == 0 ) then
         write(6,*) 'PB Bomb in pb_dbfrc_fld(): can not find projection atom/probe'
         call mexit(6, 1)
      end if

      if ( abs(insas(i,j,k)) == TWO ) then
         x(1:3) = acrd(1:3,iatm)
         rsphere = radi(iatm)
         sgn = ONE
   
         atmflag = iatm
         crg_num(atmflag) = crg_num(atmflag) + 1
         crgflag = crg_num(atmflag)
         surf_crg(crgflag,atmflag) = pol_charge(ip)

      else if ( abs(insas(i,j,k)) == ONE ) then
         x(1:3) = arccrd(1:3,-iatm)
         rsphere = dprob 
         sgn = -ONE
   
         if ( triopt == 1 ) then
            if ( -iatm > narcdot-ntri ) then
               matm = triatm(1,ntri-iatm-narcdot)
               natm = triatm(2,ntri-iatm-narcdot)
               patm = triatm(3,ntri-iatm-narcdot)
            else
               iarc = dotarc(-iatm)
               matm = arcatm(1,iarc)
               natm = arcatm(2,iarc)
            end if
         else
            iarc = dotarc(-iatm)
            matm = arcatm(1,iarc)
            natm = arcatm(2,iarc)
         end if
         atmflag = matm
         crg_num(atmflag) = crg_num(atmflag) + 1
         crgflag = crg_num(atmflag)
         surf_crg(crgflag,atmflag) = pol_charge(ip)
   
      end if
      dx = g - x
      rdx = ONE/sqrt(dx(1)**2 + dx(2)**2 + dx(3)**2)
      crd = x + dx*rdx*rsphere
!     crd = g
      pos_crg(1,crgflag,atmflag) = crd(1)
      pos_crg(2,crgflag,atmflag) = crd(2)
      pos_crg(3,crgflag,atmflag) = crd(3)
      ipos_crg(1,crgflag,atmflag) = i
      ipos_crg(2,crgflag,atmflag) = j
      ipos_crg(3,crgflag,atmflag) = k
  
   end do
   deallocate (pol_charge)
!    do i = 1, natom
!         do j = 1, crg_num(i)
!            write(401,'(4f15.6)')ipos_crg(1:3,j,i), surf_crg(j,i)
!         end do
!    end do

end subroutine pb_crgview
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine dbfrc(natom,iatm,x,dum,cnx,cny,cnz,rnx,rny,rnz)

   use solvent_accessibility, only : dotarc, arcatm, narcdot, ntri, triatm, triopt

   integer natom, iatm
   _REAL_ dum(3), x(3)
   _REAL_ cnx(natom), cny(natom), cnz(natom)
   _REAL_ rnx(natom), rny(natom), rnz(natom)

   integer, parameter :: mp = 3, np = 3
   integer i, j, iarc
   integer m, n, patm, matm, natm
   _REAL_ mvec(3), nvec(3), mdotn, mxnv(3), rmdotn2, fdotm, fdotn
   _REAL_ dfm, dfn, dum_norm(3), dum_tang(3), dumnorm
   _REAL_ a(mp,np), u(mp,np), w(np), v(mp,np), b(mp), t(np)
   _REAL_ wmax, thresh, TOL
   _REAL_ pvec(3), qvec(3), mdist, ndist, pdist, qdist

   TOL = 1.d-5
   patm = -1; natm = -1; matm = -1
   a = ZERO; u = ZERO; w = ZERO; v = ZERO; b = ZERO; t = ZERO

   ! try different force decomposition scheme
#  include "pb_dbfrc.h"
!#  include "pb_dbfrc1.h"

end subroutine dbfrc
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine gradu(l,m,n,phip,n2,n3,xp,yp,zp,up,dudx,dudy,dudz,u,phi)
!
! passed variables:
! phip: flag to choose inside or outside points
!       -1, inside
!       +1, outside
! n2: no. of grid points to be used for least square fitting
! n3: no. of unknowns to be fitted
!        4, first-order linear fitting
!       10, second-roder quadratic fitting
! (Mengjuei strongly opposes variable name phi.)

   integer l,m,n
   integer n2,n3
   _REAL_  xp,yp,zp,phip
   _REAL_  u(1:l,1:m,1:n), phi(0:l+1,0:m+1,0:n+1)
   _REAL_  up, dudx,dudy,dudz

   ! local variable

   integer job,info,isvd
   integer ix0,iy0,iz0,ix1,iy1,iz1,ix,iy,iz
   integer nsub
   integer i1, i2, j2, k2, i, j, k, ij
   integer, parameter :: nn = 100
   _REAL_  w1(1:n3,1:nn)
   !_REAL_  w2(1:n3,1:nn)
   _REAL_  b(1:nn)
   _REAL_  sd(1:n3+1)
   _REAL_  ew(1:nn),w4(1:nn)
   _REAL_  w3u(1:nn),w3ux(1:nn),w3uy(1:nn),w3uz(1:nn)
   _REAL_  ewu(1:nn),ewux(1:n3),ewuy(1:n3),ewuz(1:n3)
   _REAL_  uw(1:n3,1:n3)
   _REAL_  vl(1:nn,1:nn)
   !_REAL_  uxw(1:n3,1:n2)
   !_REAL_  uxwxv(1:n2,1:n2)
   _REAL_  dist, dist0
   _REAL_  dx,dy,dz
   _REAL_  sd2(n3),work(1000)

   ix = 0; iy = 0; iz = 0
   isvd = 2
   !w1(1:n3,1:nn) = 0.0d0; w2(1:n3,1:nn) = 0.0d0

   ! select the grid (ix,iy,iz) closest to the surface charge (xp,yp,zp)
   
   ix0 = nint(xp)
   iy0 = nint(yp)
   iz0 = nint(zp)
   dist0 = 9999.0d0
   do i = -1 ,1
      do j = -1 ,1
         do k = -1 ,1
            ix1 = ix0 + i
            iy1 = iy0 + j
            iz1 = iz0 + k
            if ( phi(ix1,iy1,iz1) > ZERO ) cycle ! only choose interior points
            dist = (xp-ix1)**2+(yp-iy1)**2+(zp-iz1)**2
            if ( dist < dist0 ) then
               ix = ix1
               iy = iy1
               iz = iz1
               dist0 = dist
            end if
         end do
      end do
   end do

   ! select the closest inside grid points to interplate E_in  
   ! and construcut the matrix for SVD

   nsub = 0
   i1 = 0
   do while ( nsub < n2 )
      do i2 = -i1, i1
         do j2 = -i1, i1
            do k2 = -i1, i1
               i = i2 + ix
               j = j2 + iy
               k = k2 + iz
               if ( i > l .or. i < 1 ) cycle
               if ( j > m .or. j < 1 ) cycle
               if ( k > n .or. k < 1 ) cycle
               if ( i2*i2 + j2*j2 + k2*k2 == i1 ) then
                  if ( phi(i,j,k)*phip >= ZERO ) then ! phip = -1 for interior points
                     nsub = nsub + 1
                     dx = (i - xp)*h
                     dy = (j - yp)*h
                     dz = (k - zp)*h
!                    if ( nsub <= n2 ) then
                        w1(1,nsub) = ONE
                        w1(2,nsub) = dx
                        w1(3,nsub) = dy
                        w1(4,nsub) = dz
                        if ( n3 == 10 ) then
                           w1(5,nsub) = 0.5d0*dx*dx
                           w1(6,nsub) = 0.5d0*dy*dy
                           w1(7,nsub) = 0.5d0*dz*dz
                           w1(8,nsub) = dx*dy
                           w1(9,nsub) = dx*dz
                           w1(10,nsub) = dy*dz
                        end if
                        b(nsub) = u(i,j,k)
!                    end if       
                 end if
              end if
            end do
         end do
      end do
      i1 = i1 + 1
   end do
!  if ( nsub > n2 ) nsub = n2

!  do i = 1, nsub
!     write(120,'(5x,5e20.6)') w1(1:4, i), b(i)
!  end do
!  do i =1,n2
!     w3ux(i) = 0.0
!     w3uy(i) = 0.0
!     w3uz(i) = 0.0
!  end do
!  do i =1,n3
!     ewux(i) = 0.0
!     ewuy(i) = 0.0
!     ewuz(i) = 0.0
!  end do

   ! call dsvdc for the singular value decomposition

!  print *,'n2=',n2
!  print *,'n3=',n3
!  print *,'b'
!  print *,v
!  w2 = w1
   job = 11

   if ( isvd == 1 ) then
      call dsvdc(w1,n3,n3,nsub,sd,ew,uw,n3,vl,nn,w4,job,info)

   else
      call dgesvd('A','A',n3,nsub,w1,n3,sd2,uw,n3,vl,nn, &
                  work,1000,info)
      do ij = 1, n3
         sd(ij) = sd2(ij)
      end do
   end if

!  print *,v
!  print *,'SVD is completed'

!  print *,'~~~~~~~~~~~~ u ~~~~~~~~~~~'
!  do i = 1, n3
!     do j = 1, n3
!        print *,sum(uw(1:n3,i)*uw(1:n3,j))
!     end do
!  end do
!  print *,'~~~~~~~~~~~~ u ~~~~~~~~~~~'

!  print *,'~~~~~~~~~~~~ v ~~~~~~~~~~~'
!  do i = 1, n2
!     do j = 1, n2
!        print *,sum(vl(1:n2,i)*vl(1:n2,j))
!     end do
!  end do
!  print *,'~~~~~~~~~~~~ v ~~~~~~~~~~~'

!  do i = 1, n3
!     do j = 1, n2
!        uxw(i,j)=sum(uw(1:n3,i)*w2(1:n3,j))
!     end do
!  end do
!  print *,'~~~~~~~~~~~~~~~~~~~~~~~~~'
!  do i = 1, n3
!     do j = 1, n2
!        uxwxv(i,j) = sum(uxw(i,1:n2)*vl(1:n2,j))
!     end do
!     print "(100f10.3)",uxwxv(i,1:n2)
!  end do
!  print *,'~~~~~~~~~~~~~~~~~~~~~~~~~'

   ! calculate E_in using the returned least-squared coefficients 

   do i = 1, n3
      if ( sd(i) > 1.0d-12 ) then
         ewu(i)  = uw(1,i)/sd(i)
         ewux(i) = uw(2,i)/sd(i)
         ewuy(i) = uw(3,i)/sd(i)
         ewuz(i) = uw(4,i)/sd(i)
      else
         ewu(i)  = ZERO
         ewux(i) = ZERO
         ewuy(i) = ZERO
         ewuz(i) = ZERO
      endif
   enddo

!  print *,v

   if ( isvd == 1 ) then
      do i = 1, nsub
         w3u(i)  = ZERO
         w3ux(i) = ZERO
         w3uy(i) = ZERO
         w3uz(i) = ZERO
         do j = 1, n3
            w3u(i)  = w3u(i)  + vl(i,j)*ewu(j)
            w3ux(i) = w3ux(i) + vl(i,j)*ewux(j)
            w3uy(i) = w3uy(i) + vl(i,j)*ewuy(j)
            w3uz(i) = w3uz(i) + vl(i,j)*ewuz(j)
         enddo
!        print *,i,v
      enddo
   else
      do i = 1, nsub
         w3u(i)  = ZERO
         w3ux(i) = ZERO
         w3uy(i) = ZERO
         w3uz(i) = ZERO
         do j = 1, n3
            w3u(i)  = w3u(i)  + vl(j,i)*ewu(j)
            w3ux(i) = w3ux(i) + vl(j,i)*ewux(j)
            w3uy(i) = w3uy(i) + vl(j,i)*ewuy(j)
            w3uz(i) = w3uz(i) + vl(j,i)*ewuz(j)
         enddo
!        print *,i,v
      enddo
   end if
!  print *,v

   up   = ZERO
   dudx = ZERO
   dudy = ZERO
   dudz = ZERO
   do i = 1, nsub
      up   = up   +  w3u(i) *b(i)
      dudx = dudx +  w3ux(i)*b(i)
      dudy = dudy +  w3uy(i)*b(i)
      dudz = dudz +  w3uz(i)*b(i)
   end do
!  print *,v


end subroutine gradu
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_charge_pol(nbnd,phi,cphi,chgrd,pol_charge,srfcrg)

   ! passed variables
   
   integer nbnd  
   _REAL_ phi(xm,ym,zm), cphi(xm,ym,zm), chgrd(xm,ym,zm)
   _REAL_ pol_charge(nbnd), srfcrg

   ! local variables

   integer i, j, k, ip
   _REAL_ total_phi(7), factor, total

   factor = INV_FOURPI*SIX*h

   total = ZERO
   do ip = 1, nbnd

      i = iepsav(1,ip); j = iepsav(2,ip); k = iepsav(3,ip)

      total_phi(1) = phi(i,  j,k) + cphi(i,  j,k)  
      total_phi(2) = phi(i-1,j,k) + cphi(i-1,j,k)  
      total_phi(3) = phi(i+1,j,k) + cphi(i+1,j,k)  
      total_phi(4) = phi(i,j-1,k) + cphi(i,j-1,k)  
      total_phi(5) = phi(i,j+1,k) + cphi(i,j+1,k)  
      total_phi(6) = phi(i,j,k-1) + cphi(i,j,k-1)  
      total_phi(7) = phi(i,j,k+1) + cphi(i,j,k+1)  

      pol_charge(ip) = factor*( total_phi(1)-&
           SIXTH*( total_phi(2)+total_phi(3)+&
                   total_phi(4)+total_phi(5)+&
                   total_phi(6)+total_phi(7) ) )
      pol_charge(ip) = ( pol_charge(ip) - chgrd(i,j,k) )*frcfac*INV_AMBER_ELECTROSTATIC
      total = total + pol_charge(ip) 
   end do
   srfcrg = total


end subroutine get_charge_pol

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine get_coulomb(natom,crg,crd,ac,ax,ay,az,fx,fy,fz,eelrf,coulomb)

   ! passed variables

   integer natom
   _REAL_ crg, crd(1:3), coulomb(1:3), eelrf
   _REAL_ ac(natom), ax(natom), ay(natom), az(natom)
   _REAL_ fx(natom), fy(natom), fz(natom)

   ! local variables

   integer jatm
   _REAL_ dinv, d2inv, de, dx(1:3), dff, dcc

   coulomb(1:3) = ZERO
   do jatm = 1, natom
      if ( ligand .or. multiblock ) then
         if (liveflag(jatm) == 0) cycle
      endif
      !if ( ligand .and. liveflag(jatm) == 0 ) cycle
      !if ( multiblock .and. liveflag(jatm) == 0 ) cycle
      dx(1) = crd(1) - ax(jatm)
      dx(2) = crd(2) - ay(jatm)
      dx(3) = crd(3) - az(jatm)
      d2inv = ONE/(dx(1)**2 + dx(2)**2 + dx(3)**2); dinv = sqrt(d2inv)

      de = ac(jatm)*dinv*AMBER_ELECTROSTATIC
      dcc = de*d2inv
      dff = crg*dcc

      ! calculate reaction field energy

      eelrf = eelrf + de*crg

      ! calculate the QE force on the atom

      fx(jatm) = fx(jatm) - dx(1)*dff
      fy(jatm) = fy(jatm) - dx(2)*dff
      fz(jatm) = fz(jatm) - dx(3)*dff

      ! calculate the coulomb field on the surface charge

      coulomb(1) = coulomb(1) + dx(1)*dcc
      coulomb(2) = coulomb(2) + dx(2)*dcc
      coulomb(3) = coulomb(3) + dx(3)*dcc
   end do


end subroutine get_coulomb

#if !defined SANDER && !defined LIBPBSA
#ifdef MPI
! Mengjuei Hsieh, University of California Irvine
subroutine pbslave_init(natom)
   implicit none
#  include "flocntrl.h"
#  include "pb_md.h"
  include "mpif.h"
#  include "parallel.h"
#  include "extra.h"
   integer natom
   _REAL_ green(0:20, 0:20, 0:20)
   common /blk_green/ green
   !LOCAL
   integer ierr
   !MPI initialization VERY BAD IMPLEMENTATION
   if ( .not. master ) then
      pbprint = .false.
      allocate(liveflag(    natom),STAT=ierr)
      allocate(realflag(    natom),STAT=ierr)
      allocate(outflag (    natom),STAT=ierr)
      allocate(outflagorig( natom),STAT=ierr)
      allocate(  iar1pb(4,0:natom),STAT=ierr)
      allocate(   nshrt(  0:natom),STAT=ierr)
      allocate(    acrd(3,  natom),STAT=ierr)
      allocate(    gcrd(3,  natom),STAT=ierr)
      allocate(    icrd(3,  natom),STAT=ierr)
      allocate(    acrg(    natom),STAT=ierr)
      allocate(    gcrg(8,  natom),STAT=ierr)
      allocate(  grdcrg(3,8*natom),STAT=ierr)
      allocate( qgrdcrg(  8*natom),STAT=ierr)
      allocate(  mapout(    natom),STAT=ierr)
      allocate(     nex(    natom),STAT=ierr)
      allocate(     iex(64, natom),STAT=ierr)
   endif
   call MPI_BCAST(    ligand,          1,MPI_LOGICAL,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(multiblock,          1,MPI_LOGICAL,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     srsas,          1,MPI_LOGICAL,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    maxnbr,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    maxnba,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    nfocus,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    fscale,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(   nbuffer,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(   solvopt,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     bcopt,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    eneopt,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    intopt,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    maxitn,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  ngrdblkx,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  ngrdblky,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  ngrdblkz,          1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  nshrt(0),    natom+1,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    nex(1),      natom,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  iex(1,1),      natom,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    do_dir,BC_FLOCNTRL,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    buffer,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    cutres,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     cutnb,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     cutsa,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     cutfd,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(      offx,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(      offy,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(      offz,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST( fillratio,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    fmiccg,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    accept,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    pbtemp,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(  ivalence,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(      wsor,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     lwsor,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(     epsin,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(    epsout,      1,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST( acrd(1,1),natom*3,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(green(0,0,0), 9261,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BARRIER( CommSANDER, ierr );REQUIRE(ierr==0)
   call MPI_BCAST(savbcopt(1),nfocus,MPI_INTEGER,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BCAST(   savh(1), nfocus,MPI_DOUBLE_PRECISION,0,CommSANDER,ierr); REQUIRE(ierr==0)
   call MPI_BARRIER( CommSANDER, ierr );REQUIRE(ierr==0)
   if ( .not. master ) then
      !allocate(iprlong (maxnbr),stat=ierr)
      allocate(iprshrt (maxnba),stat=ierr)
      allocate(cn1pb   (maxnba),stat=ierr)
      allocate(cn2pb   (maxnba),stat=ierr)
      allocate(cn3pb   (maxnba),stat=ierr)
   end if
end subroutine pbslave_init
#endif /*def MPI*/
#endif /*ndef SANDER or LIBPBSA*/

end module poisson_boltzmann
