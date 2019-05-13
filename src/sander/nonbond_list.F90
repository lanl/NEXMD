 !<compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

module nblist

_REAL_, dimension(:,:), allocatable, save :: imagcrds,fraction,savfrac,dfrac, &
                                             savcrd
integer, dimension(:),  allocatable, save :: indatg,atmcell,numatg,indoff, &
                                             bckptr,nlogrid,nhigrid,my_grids, &
                                             lstmask,nummask,maskptr,iwa,iwh,iws, &
                                             numvdw,numhbnd,numsc,atmlist, &
                                             mygrdlist, numimg
                                             
integer, dimension(:), allocatable, save :: indexlo,indexhi
integer, dimension(:), allocatable, save :: nvdwcls
integer, dimension(:,:),  allocatable, save :: nghbptr,nghtran
integer, dimension(:),    allocatable, save :: itran
integer, save :: nucgrd1_0,nucgrd2_0,nucgrd3_0
integer, dimension(1:7,1:10), save :: xtran
integer, save :: myindexlo,myindexhi,inddelta,nblist_allint,nblist_allreal
integer,private,allocatable,dimension(:),save :: exclude

logical :: first_list_flag
integer, save :: num_calls_nblist=0

_REAL_, dimension(1:3,1:18), save :: tranvec

!-------------  ew_unitcell.h --------------------------
integer, parameter :: BC_EWUCR=54, BC_EWUCI=3

_REAL_ :: a,b,c,alpha,beta,gamma,volume
_REAL_ :: ucell,recip,dirlng,reclng,sphere,cutoffnb
_REAL_ :: olducell,oldrecip
_REAL_ :: skinnb,cutlist,nbfilter
common/unitcell/ &
      ucell(3,3),recip(3,3),dirlng(3),reclng(3),   &! 24
      olducell(3,3),oldrecip(3,3),                 &! 42
      a,b,c,alpha,beta,gamma, volume,              &! 49
      sphere,cutoffnb,skinnb,cutlist,nbfilter     ! 54

integer :: nbflag,nbtell,steps_since_list_build
common/ew_upd/nbflag,nbtell,steps_since_list_build

!END----------  ew_unitcell.h --------------------------
!-------------  ew_localnb.h --------------------------
! UNIT CELL  GRID of mapped unit cell coords for preimaging
! NUCGRD1 is number of cells along 1st direction
! NUCGRD2 is number of cells along 2nd direction
! NUCGRD3 is number of cells along 3rd direction
! The subcell neighborhood is obtained by considering
! cells within +- NGHB1 in the first direction
! cells within +- NGHB2 in the second direction, and
! cells within +- NGHB3 in the third direction

! The distance between parallel faces of a subcell is then
!  reclng(1)/NUCGRD1, reclng(2)/NUCGRD2 or
!  reclng(3)/NUCGRD3
! Thus the short range cutoff is the minimum of
!  NGHB1*reclng(1)/NUCGRD1,NGHB2*reclng(2)/NUCGRD2 and
!  NGHB3*reclng(3)/NUCGRD3
! MXATCELL is the maximum number of atoms per subcell
! NOTE!!!!!! YOU MUST HAVE
!           NGHB1 < NUCGRD1, NGHB2 < NUCGRD2, and NGHB3 < NUCGRD3
! imagptr(i) is the subcell number in the imaged grid corresponding to
! subcell i in the unit cell grid
! nghbptr is array from which the neighbor atoms
! in the image grid can be rapidly retrieved.


!----------------------------------------------------------
!  neighbor cell #
integer nghb
parameter (nghb = 3)
!ccccccccccccccccccccccccccccccccccccccccccccccccc
integer, parameter :: BC_DIRPARS=13
integer :: numnptrs, nucgrd1,nucgrd2,nucgrd3,nucgmax,nucgrd
integer :: nghb1,nghb2,nghb3
integer :: maxnblst,maxnptrs,maximage,mxlstmsk

common/dirpars/ &
      numnptrs,nucgrd1,nucgrd2, &
      nucgrd3,nghb1,nghb2,nghb3,nucgmax,nucgrd, &
      maxnblst,maxnptrs,maximage,mxlstmsk


!ccccccccccccccccccccccccccccccccccccccccccccccccc
integer,parameter :: BC_NB_GINDEX = 2
integer gindexlo,gindexhi
common/nb_gindex/gindexlo,gindexhi
!END----------  ew_localnb.h --------------------------



contains

!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Nonbond list management.

!-----------------------------------------------------------------------
!     --- NBLIST_ALLOCATE ---
!-----------------------------------------------------------------------
subroutine nblist_allocate(natom,ntypes,num_direct,numtasks)
   implicit none
   integer, intent(in) :: natom,ntypes,num_direct,numtasks

   integer :: allocate_error,allint,allreal

   allreal=0   !report back how many reals allocated
   allint=0    !                     ints

   allocate(imagcrds(1:3,1:natom), stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint+3*natom

   allocate(fraction(1:3,1:natom), savfrac(1:3,1:natom), dfrac(1:3,1:natom), &
            savcrd(1:3,1:natom), &
            stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allreal=allreal+12*natom


   allocate (indatg(1:natom),atmcell(1:natom),atmlist(1:natom), &
             stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint+3*natom

   allocate (nlogrid(1:nucgmax), nhigrid(1:nucgmax), my_grids(1:nucgmax), &
            indoff(1:nucgmax), numimg(1:nucgmax), &
            numatg(1:nucgmax), &
            stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint+6*nucgmax


   allocate (nghbptr(0:maxnptrs,1:nucgmax),bckptr(1:natom), &
            stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint+natom + (maxnptrs+1)*nucgmax

   allocate (nummask(1:natom),maskptr(1:natom), lstmask(1:mxlstmsk), &
            stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint = allint + 2*natom + mxlstmsk

   allocate (iwa(1:natom),iwh(1:natom),numvdw(1:natom),numhbnd(1:natom), &
            stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint+4*natom

#ifdef MPI /* SOFT CORE */
   allocate (iws(1:natom),numsc(1:natom),stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint+2*natom
#endif

   allocate (itran(1:natom), &
            stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint+natom

   allocate(indexlo(0:numtasks-1),indexhi(0:numtasks-1), stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint+2*numtasks

 !---- private arrays ---------------------------------
   allocate (exclude(1:natom),mygrdlist(1:natom), &
            stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint+2*natom
   
#ifndef MPI
   allocate (   nghtran( 1:maxnptrs , 1:nucgmax+1 ), &
            stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint + maxnptrs*(nucgmax+1)
#else
   allocate (   nghtran( 1:maxnptrs , 1+nucgmax/num_direct ), &
            stat=allocate_error)
      REQUIRE (allocate_error == 0)
      allint=allint + maxnptrs*(nucgmax/num_direct+1)
#endif 


   allocate (nvdwcls(1:ntypes), stat=allocate_error)
   REQUIRE (allocate_error == 0)
   allint = allint + ntypes

   nblist_allint  = allint
   nblist_allreal = allreal

   return
end subroutine nblist_allocate

!-----------------------------------------------------------------------
!     --- NBLIST_DEALLOCATE ---
!-----------------------------------------------------------------------
subroutine nblist_deallocate()
   implicit none

   integer :: allocate_error

   if(.not. allocated(imagcrds) ) return ! None of this memory has been
                                         !  allocated, just return

   deallocate(imagcrds, stat=allocate_error)
      REQUIRE (allocate_error == 0)

   deallocate(fraction, savfrac, dfrac, savcrd, stat=allocate_error)
      REQUIRE (allocate_error == 0)

   deallocate (indatg,atmcell,atmlist, stat=allocate_error)
      REQUIRE (allocate_error == 0)

   deallocate(nlogrid,nhigrid,my_grids,indoff,numimg,numatg, &
          stat=allocate_error)
      REQUIRE (allocate_error == 0)

   deallocate (nghbptr,bckptr,stat=allocate_error)
      REQUIRE (allocate_error == 0)

   deallocate (nummask,maskptr, lstmask,stat=allocate_error)
      REQUIRE (allocate_error == 0)

   deallocate (iwa,iwh,numvdw,numhbnd,stat=allocate_error)
      REQUIRE (allocate_error == 0)

#ifdef MPI /* SOFT CORE */
   deallocate (iws, numsc, stat=allocate_error)
      REQUIRE (allocate_error == 0)
#endif

   deallocate (itran,stat=allocate_error)
      REQUIRE (allocate_error == 0)

   deallocate(indexlo,indexhi,stat=allocate_error)
      REQUIRE (allocate_error == 0)

 !---- private arrays ---------------------------------
   deallocate (exclude,mygrdlist,stat=allocate_error)
      REQUIRE (allocate_error == 0)
   
   deallocate (   nghtran,stat=allocate_error)
      REQUIRE (allocate_error == 0)

   deallocate (nvdwcls,stat=allocate_error)
   REQUIRE (allocate_error == 0)

   return
end subroutine nblist_deallocate

!
!-----------------------------------------------------------------------
!     --- NONBOND_LIST ---
!-----------------------------------------------------------------------
!     Handles set-up and error checking for calling of
!     get_nb_list which creates the nonbond list.

subroutine nonbond_list(crd,iac,ico,iblo,inb,ntypes, &
      natom,x,ix, &
      ipairs,ntnb,ibelly,belly,newbalance,cn1, &
      v,vold,ntp,xr,qsetup,&
      do_list_update)
   use amoeba_mdin,only:iamoeba
   use trace
   implicit none

! ----- INPUT variables ------------------------   
   integer :: natom,ntnb,newbalance
   integer :: ipairs(*),ix(*),ibelly(*)
   _REAL_  :: crd(3,natom),cn1(*)
   integer :: iac(natom),ico(*),iblo(*),inb(*),ntypes
   _REAL_  :: x(*)
   logical qsetup



   integer :: ngrd1,ngrd2,ngrd3,sizgrdprs
   integer last_numlist,isiz
   logical belly,trial,balance
#  include "extra.h"
#  include "ew_cntrl.h"
#  include "ew_mpole.h"
   
#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
   integer listdiff, listdiffmax
   integer indexlo(0:numtasks-1),indexhi(0:numtasks-1)
# ifdef MPI_DOUBLE_PRECISION
#  undef MPI_DOUBLE_PRECISION
# endif
   include 'mpif.h'
# ifdef CRAY_PVP
#  define MPI_DOUBLE_PRECISION MPI_REAL8
# endif
   integer tmplist(0:MPI_MAX_PROCESSORS), &
         alllist(0:MPI_MAX_PROCESSORS)
#else   /* not parallel needs numtasks and mytaskid */
   integer numtasks,mytaskid,indexlo,indexhi
#endif
   
!#  include "ew_tran.h"
#  include "def_time.h"
   
   integer ier
   integer, dimension(:), allocatable :: nnghbptr 
   integer, dimension(:), allocatable :: nghtranptr
   integer, dimension(:), allocatable :: gridpairs
   
   !     ARGUMENTS: (all are input)
   !     CRD is the array holding atomic coords.
   !     IAC and ICO are used to look up vdw and hbond interactions.
   !     ICO is an NTYPES by NTYPES array giving lookup into VDW and HBOND
   !        coefficients.
   !     IAC(i) is the atom type of atom i. The coefficients
   !        for 6-12 interactions for atoms i and j are indexed by
   !        ICO(iac(i),iac(j)). In practice ICO is unrolled to a 1
   !        dimensional array.  They are needed here to split nonbond
   !        list into vdw,hbonds (see pack_nb_list())
   !     IBLO and INB are used to mask out some nonbond interactions
   !        (bonded pairs,1-3,1-4 pairs) etc.  IBLO(i) stores the number
   !        of atoms masked by atom i.  INB stores their indices.  The
   !        masked pairs are stored back to back in INB.
   !        For our  purposes we need to produce our own variants.
   !        See "load_mask".
   !     NTYPES is used with IAC and ICO
   !     X and IX are real and integer arrays which comprise the total "dynamic"
   !        memory  for amber; coords,forces,bond lists etc are accessed as
   !        offsets in them.
   
   integer ifail,listtot,listtotall,task
#ifdef MPI
   integer inddel, i, ierr
#endif
   logical, intent(out) :: do_list_update
   _REAL_ v(*),xr(*),vold(*)
   integer ntp
   save trial,balance
   save last_numlist
#ifdef MPI
   save listdiffmax
#endif
   
   !---------------------------- code starts here ---------------------------
   call trace_enter( 'nonbond_list' )
#ifndef MPI
   mytaskid=0
   numtasks=1
#endif


   !  ---------------------------------------------------------------------
   !  ------        FIRST TIME SETUP STUFF --------------------------------
   !  -------- If qsetup is true then this is the
   !   ------- first time through, do the setup stuff
   !    ------ do not set the qsetup flag false till after the
   !     ----- if(.not. qsetup ) block
   
   if ( qsetup ) then
      call ew_startup(natom,iblo,inb,x,ix)
      nucgrd1_0=-1
      nucgrd2_0=-1
      nucgrd3_0=-1
      call save_crds(natom,crd)

      balance=.false.
      newbalance=0
      if(periodic == 0) balance= (numtasks > 1)
      if(balance)newbalance=2
      last_numlist=0
      steps_since_list_build = 0
   end if
   !  ------   DONE FIRST TIME SETUP STUFF --------------------------------
   
   trial=newbalance > 0
   
   if ( .not. qsetup )then
     
      ! The skin check is done in runmd for spatial, do_list_update is passed
      !   in as optional arg
         !----------------------------------------------------------------
         !   ---- SKIN (BUFFER) CHECK    see if new list is required
         !        Do this block except on first pass. 
         !
         !   ---- if the nocutoff flag is set, then there is no need to
         !        update the list, all pairs are already determined
      
         if (nocutoff) then
            call trace_exit( 'nonbond_list' )
            return
         end if
      
         !   if nbflag=0, use old logic; i.e. only update if ntnb .ne. 0
         !   if nbflag=1, the skin check determines it
      
         if ( nbflag == 0 )then
            if ( nbtell /= 0 ) then ! list update info requested
               if(master)write(6,'(1x,A,I2)') 'OLD LIST LOGIC; ntnb = ',ntnb
            end if
            if ( ntnb == 0 ) then
               call trace_exit( 'nonbond_list' )
               return
            end if
         else
            call check_skin(crd,do_list_update)
            if ( do_list_update ) then
               call save_crds(natom,crd)
            else
               call trace_exit( 'nonbond_list' )
               return
            end if
         end if
   end if  ! ( .not. qsetup )
   qsetup = .false.

   !------------------------------------------------------------------------
   !            START THE LIST BUILD, FIRST THE LIST GRID SETUP
   !---------------------------------------------------------------------

   num_calls_nblist = num_calls_nblist + 1

   call timer_start(TIME_BLDLST)
   if ( master .and. nbtell /= 0 ) then
      ! list update info requested
      write(6, '(1x,A,I7,A)') 'Building list: ', &
            steps_since_list_build, ' steps since previous list build.'
   end if
   steps_since_list_build = 0
#ifdef MPI
   if(i_do_direct)then
      call mpi_comm_size(direct_comm,numtasks,ierr)
      call mpi_comm_rank(direct_comm,mytaskid,ierr)
#endif
   !  ---- Get nb grid information
   !       This is necessary in the case where the system dimensions are
   !       changing: nonperiodic and const pressure
      call map_coords(crd,natom,recip)
      call save_frac_crds(natom)
      call setup_grids(periodic,nogrdptrs,verbose)
      if ( nucgrd1*nucgrd2*nucgrd3 > nucgmax ) &
            call sander_bomb('nonbond_list', &
            ' volume of ucell too big, too many subcells', &
            ' list grid memory needs to be reallocated, restart sander')
      
      !  --------- Nonperiodic systems: --------------------------------------
      if(periodic == 0)then
         !       ----- If system has too few cells for the pointer method
         !       -----   to be efficient, use the no-grid-pointer system
         !       -----   where all forward cells are checked for pair atoms
         !       -----   rather than just the forward cells on the grid-pointer-list
         
         nogrdptrs=nogrdptrs .or.( (nucgrd1 <= 2*nghb+2) .or. &
               (nucgrd2 <= 2*nghb+2) .or. &
               (nucgrd3 <= 2*nghb+2) )
         if(verbose >= 3)write(6,*)" List using nogrdptrs ",nogrdptrs
         !       ----- Make the subcells about 3 A in size now so that there
         !       ----- are more subcells, easier to load balance.
         if(nogrdptrs)then
            if(dirlng(1)/nucgrd1 > 3)then
               ngrd1 = dirlng(1)/3
               ngrd2 = dirlng(2)/3
               ngrd3 = dirlng(3)/3
               if ( verbose >= 1)then
                  write(6,'("|    New Grids set up for nogrdptrs ")')
                  write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)') &
                        'Number of grids per unit cell in x,y,z:', &
                        ngrd1, ngrd2, ngrd3
               end if
               nucgrd1=ngrd1
               nucgrd2=ngrd2
               nucgrd3=ngrd3
            end if
         end if
      end if
      ! ---- end of nonperiodic grid stuff ------------------------
      
      call fill_tranvec()
      nucgrd = nucgrd1*nucgrd2*nucgrd3

#ifdef MPI
      
      !---------------------------------------------------------------------
      !      For trial run of a parallel run, an initial guess at
      !         the distribution of subcells must be made. (trial=.true.)
      !      For a periodic system, the balance is assumed to be
      !         good for evenly dividing the subcells. (balance = .false.)
      if(iamoeba.eq.1) then
         myindexlo=1
         myindexhi=nucgrd
      else
         if(trial .or. .not.balance)then
            inddel = (nucgrd-1) / numtasks + 1
            if ( inddel == 0 )inddel = 1
            myindexlo = 1+(mytaskid)*inddel
            myindexhi = myindexlo+inddel-1
            if ( mytaskid == numtasks-1 ) myindexhi = nucgrd
            last_numlist = 0
         end if
      end if
#else
      !---------------------------------------------------------------------
      !     ----- for non parallel runs, do all the subcells
      myindexlo=1
      myindexhi=nucgrd
#endif

      !---------------------------------------------------------------------
      !----- assign atoms to cells (atmcell) and make cell atom list (indatg)
      call grid_ucell(natom, periodic)
      !---------------------------------------------------------------------
      if(periodic == 0)then
         !----- Nonperiodic systems:
         isiz=max(2,(myindexhi-myindexlo+1)*(maxnptrs+1))
         if(nogrdptrs)isiz=1
         allocate( nnghbptr(isiz), stat=ier )
         REQUIRE( ier==0 )

         isiz=max(2,(myindexhi-myindexlo+1)*(maxnptrs))
         if(nogrdptrs)isiz=1
         allocate( nghtranptr(isiz), stat=ier )
         REQUIRE( ier==0 )

         call grid_pointers( &
               nnghbptr, &
               nghtranptr, &
               periodic,nogrdptrs)
      else
         !-----  Periodic systems:
         !            only call grid_pointers  if the unit cell has changed
         !            so much that nucgrd[123] have changed. Otherwise the 
         !            grid pointers do not change.
         if(nucgrd1 /= nucgrd1_0 &
               .or. nucgrd2 /= nucgrd2_0 &
               .or. nucgrd3 /= nucgrd3_0)then
            nucgrd1_0=nucgrd1
            nucgrd2_0=nucgrd2
            nucgrd3_0=nucgrd3
            call grid_pointers( &
                  nghbptr, &
                  nghtran, &
                  periodic,nogrdptrs)
         end if
      end if  ! (periodic == 0)
      
      !  ---------------------------------------------------------------------
      
      call grid_image( natom,verbose,periodic )

      if(balance)then
         sizgrdprs=nucgrd*2
      else
         sizgrdprs=2
      end if
      allocate( gridpairs(sizgrdprs), stat=ier )
      REQUIRE( ier==0 )
      if(periodic == 0)then
         call get_nb_list(iac,ico,ntypes,ifail, &
               listtot,natom, &
               ipairs, &
               nnghbptr, &
               verbose, &
               nghtranptr,tranvec, &
               belly,ibelly,balance,gridpairs, &
               periodic,nogrdptrs,cn1)
      else
         call get_nb_list(iac,ico,ntypes,ifail, &
               listtot,natom, &
               ipairs, &
               nghbptr, &
               verbose, &
               nghtran,tranvec, &
               belly,ibelly,balance,gridpairs, &
               periodic,nogrdptrs,cn1)
      end if  ! (periodic == 0)

      if ( ifail == 1) then
         write(6, '(5x,a,i10)') 'SIZE OF NONBOND LIST = ', listtot
         call sander_bomb('nonbond_list','Non bond list overflow!', &
               'check MAXPR in locmem.f')
      end if

#ifdef MPI
      if(trial)then
         ASSERT ( periodic == 0 )
         call trace_mpi('mpi_allreduce', nucgrd, 'MPI_INTEGER', mpi_sum)
         call mpi_allreduce(gridpairs, &
               gridpairs(1+nucgrd),nucgrd, &
               mpi_integer, mpi_sum,commsander,ierr)
         call fix_grid_balance( &
               gridpairs(1+nucgrd), &
               nucgrd,numtasks,mytaskid, &
               listdiffmax)
         
         deallocate( gridpairs, nghtranptr, nnghbptr )
         
         isiz=max(2,(myindexhi-myindexlo+1)*(maxnptrs+1))
         if(nogrdptrs)isiz=1
         allocate( nnghbptr(isiz), stat=ier )
         REQUIRE( ier==0 )

         isiz=max(2,(myindexhi-myindexlo+1)*(maxnptrs))
         if(nogrdptrs)isiz=1
         allocate( nghtranptr(isiz), stat=ier )
         REQUIRE( ier==0 )

         allocate( gridpairs(2*nucgrd), stat=ier )
         REQUIRE( ier==0 )
         
         call grid_pointers( &
               nnghbptr, &
               nghtranptr, &
               periodic,nogrdptrs)
         call grid_image( natom,verbose,periodic)

         call get_nb_list(iac,ico,ntypes,ifail, &
               listtot,natom, &
               ipairs, &
               nnghbptr, &
               verbose, &
               nghtranptr,tranvec, &
               belly,ibelly,balance,gridpairs, &
               periodic,nogrdptrs,cn1)
         trial=.false.
         last_numlist=listtot
         newbalance=0
      else if (balance) then
         
         !        --- If this was not a rebalance step, then check
         !            whether the list has gotten more than a tolerance
         !            away from the last balance. If so, trigger
         !            a new balance by setting trial and newbalance
         !            newbalance will need to be communicated to
         !            all processors before the next list build.
         
         listdiff=iabs(last_numlist-listtot)
         if(listdiff > listdiffmax)then
            newbalance=2
         end if
      end if  ! (trial)
      
      if( (nbtell > 1 ) .or.  ( balance .and. (verbose > 0))  )then
         if(master) write(6,105)
         if(master) write(6,106)
         do i=0,numtasks-1
            tmplist(i)=0
         end do
         tmplist(mytaskid)=listtot
         call trace_mpi('mpi_allreduce', numtasks, 'MPI_INTEGER', mpi_sum)
         call mpi_allreduce(tmplist(0),alllist(0),numtasks, &
               mpi_integer, mpi_sum,commsander,ierr)
         if(master) then
            write(6,110)(i,alllist(i),i=0,numtasks-1)
            do i=1,numtasks-1
               tmplist(0)=tmplist(0)+alllist(i)
            end do
            write(6,103)tmplist(0)
         end if
         110 format("|       ",2i12)
         103 format("|            Total: ",i12)
         105 format("| ",'----------------List Breakdown----------------')
         106 format("|       ",'list processor   listtot')
      end if
#else
      if(master)then
         if ( verbose > 0)write(6,*)'listtot = ',listtot
      end if
#endif
      
      deallocate( gridpairs )
      if( periodic == 0 )then
         deallocate( nghtranptr, nnghbptr )
      end if
      listtotall=listtot
#ifdef MPI
      if(first_list_flag .or. (nbtell >= 1))then
         call trace_mpi('mpi_allreduce', 1, 'MPI_INTEGER', mpi_sum)
         call mpi_allreduce(listtot,listtotall,1, &
               mpi_integer, mpi_sum,commsander,ierr)
      end if
      call mpi_comm_rank(commsander,mytaskid,ierr)
      call mpi_comm_size(commsander,numtasks,ierr)
   end if  ! (i_do_direct)
   call mpi_barrier(commsander,ierr)
#endif
!            '| Local SIZE OF NONBOND LIST = ', listtot,'   PE',mytaskid
   if( (first_list_flag .or. (nbtell == 1))  .and. master ) then
      write(6, '(a,i10)') &
            '| Local SIZE OF NONBOND LIST = ', listtot
      write(6, '(a,i10)') &
            '| TOTAL SIZE OF NONBOND LIST = ', listtotall
   end if
   call amflsh(6)
   first_list_flag = .false.
   call timer_stop(TIME_BLDLST)
   call trace_exit( 'nonbond_list' )
   return
end subroutine nonbond_list 

!----------------------------------------------------------------------------
!     --- GET_NB_LIST ---    *********
!     ...this routine carries out the short range part of nonbond
!     calculation.  Note that I am using a gridding routine, or
!     geometric hashing with renumbering of atoms to speed list
!     generation and to promote locality of reference, important
!     for cache memory usage.  The pre-imaging is to avoid expensive
!     coordinate transforms in the minimum image pair calculations.
!     I set up pointers between the subcells to speed this process.
!     NOTE: this comment also refers to routine pack_nb_list().
!----------------------------------------------------------------------------
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine get_nb_list here]
subroutine get_nb_list(iac,ico,ntypes,ifail, &
      listtot,natom, &
      ipairs, &
      nghbptr, &
      verbose, &
      nghtran, &
      tranvec, &
      belly,ibelly,balance,gridpairs,periodic,nogrdptrs,cn1)

   !     ---main routine for list building
   !     note it is structured for parallelism
   implicit none
   
#  include "extra.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
   include 'mpif.h'
#else   /* not parallel needs numtasks and mytaskid */
   integer numtasks,mytaskid
#endif
   
   _REAL_ cn1(*)
   _REAL_ tranvec(*)
   integer iac(*),ico(*),ntypes,ifail,listtot
   integer natom
   integer ipairs(maxnblst)
   integer ibelly(*), gridpairs(nucgrd,2)
   logical belly,balance
#ifdef MPI
   integer np0
#endif

   integer nghbptr(0:maxnptrs,*), verbose
   
   integer nghtran(maxnptrs,*)
   integer numlist
   integer index,index0,index1
   integer j,jj,j0,k,ncell_lo,ncell_hi,m,m1,m2
   integer i,kk,numpack
   integer nstart,tranyz,xtindex,jtran,indi,nstop,numindex
   _REAL_ xk,yk,zk
   _REAL_ cutoffsq
   integer periodic
   
   logical nogrdptrs
   
#ifndef MPI
   numtasks=1
   mytaskid=0
#endif
   
   mygrdlist(1:natom)=0
   ifail = 0
   cutoffsq = cutlist*cutlist 
   exclude(1:natom) = 0  

   numpack = 0
# ifdef MPI
   if(balance)then
      np0=0
      do i=1,nucgrd
         gridpairs(i,1)=0
      end do
   end if
# endif
   do index = myindexlo,myindexhi
      index0=index-myindexlo+1
      if ( numimg(index) > 0 )then
         ncell_lo = nlogrid(index)
         ncell_hi = nhigrid(index)
         numlist = 0
         !==========================================================
         if(nogrdptrs)then
         !==========================================================
            if(periodic /= 0) &
                  call sander_bomb('get_nb_list', &
                  'Cannot run nogrdptrs with ntb=1', &
                  'turn off nogrdptrs ')

            do j0 = index,nucgrd
               m1 = nlogrid(j0)
               m2 = nhigrid(j0)
               if ( m2 >= m1 )then
                  do m = m1,m2
                     numlist = numlist+1
                     atmlist(numlist) = m
                  end do
               end if
            end do
            if ( numlist > 0 )then
               do k = ncell_lo,ncell_hi
                  kk = k - ncell_lo + 1
                  i = bckptr(k)
                  xk = imagcrds(1,k)
                  yk = imagcrds(2,k)
                  zk = imagcrds(3,k)
                  mygrdlist(k)=1
                  call pack_nb_nogrdptrs(kk,i,xk,yk,zk,imagcrds,cutoffsq, &
                        numlist,numpack, &
                        iac,ico,ntypes,ipairs,ifail, &
                        belly,ibelly)
                  if ( ifail == 1)then
                     listtot = maxnblst
                     return
                  end if
               end do
            end if
            
            !==========================================================
         else   ! nogrdptrs if statement
            !==========================================================
            !  -- get list of atoms in the neighborhood of cell(index0)
            
            numindex= numnptrs
            if(periodic == 0) numindex = nghbptr(0,index0)
            do j0 = 1,numindex
               index1 = nghbptr(j0,index0)
               if ( j0 == 1 )then
                  nstart=0
               else
                  nstart=-nghb1
               end if
               nstop=nghb1
               if(periodic == 0)then
                  indi=mod(index1-1,nucgrd1)+1
                  if(indi > nucgrd1-nghb1)nstop=nucgrd1-indi
                  if(indi <= nghb1)nstart=1-indi
                  if ( j0 == 1 ) nstart=0
               end if
               xtindex=ishft(nghtran(j0,index0),-8)
                  REQUIRE(xtindex <= 10 .and. xtindex >= 0)
               tranyz = nghtran(j0,index0)-ishft(xtindex,8)

               do j = nstart,nstop
                  jtran = tranyz+xtran(j-nstart+1,xtindex)
                  jj = index1+j-xtran(j-nstart+1,xtindex)*nucgrd1
                  m1 = nlogrid(jj)
                  m2 = nhigrid(jj)
                  if ( m2 >= m1 )then
                     do m = m1,m2
                        numlist = numlist+1
                        atmlist(numlist) = m
                        itran(numlist)=jtran
                     end do
                  end if
               end do
            end do  !  j0 = 1,numindex

            if ( numlist > 0 )then
               do k = ncell_lo,ncell_hi
                  kk = k - ncell_lo + 1
                  i = bckptr(k)
                  xk = imagcrds(1,k)
                  yk = imagcrds(2,k)
                  zk = imagcrds(3,k)
                  mygrdlist(k)=1
                  call pack_nb_list(kk,i,xk,yk,zk,imagcrds,cutoffsq, &
                        numlist,numpack, &
                        iac,ico,ntypes,ipairs,ifail, &
                        tranvec, &
                        belly,ibelly,cn1)
                  if ( ifail == 1)then
                     listtot = maxnblst
                     return
                  end if
               end do
            end if
         !==========================================================
         end if   ! nogrdptrs if statement
         !==========================================================
      end if  ! ( numimg(index) > 0 )
#  ifdef MPI
      if(balance)then
         gridpairs(index,1)=numpack-np0
         np0=numpack
      end if
#  endif
   end do  !  index = myindexlo,myindexhi
   listtot = numpack
   return
end subroutine get_nb_list 

!-------------------------------------------------------------------
!     --- GRID_IMAGE ---
!     ...this routine grids the imaged atoms in the unit cell plus
!     neighbor cells while building the array of sorted image atoms.
!     The sorting should improve locality of reference...
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine grid_image here]
subroutine grid_image(natom, verbose, periodic)

   use constants, only : zero, half
   implicit none
   integer, intent(in) ::  natom,verbose,periodic
   _REAL_  f1,f2,f3,shft
   integer index,i,j,n
   integer i0,i1,i2,numimage
#  include "extra.h"

   shft=half
   if(periodic == 0)shft=zero
   numimage = 0
   do index = 1,nucgrd
      if(my_grids(index) == 1)then

         n = numatg(index)
         if (numimage + n > natom ) then

            write(6,*)'natom = ',natom
            write(6,*)'numimage = ',numimage
            call sander_bomb('grid_image', &
                  'num image atoms exceeds natom!!', &
                  '<ew_direct.f>grid_image() ')
         end if
         i0 = indoff(index)
         i1 = i0 + 1
         i2 = i0 + n
         nlogrid(index) = numimage + 1
         nhigrid(index) = numimage + n
         numimg(index) = n
            do i = i1,i2
               j = indatg(i)
               f1 = fraction(1,j)+shft
               f2 = fraction(2,j)+shft
               f3 = fraction(3,j)+shft
               numimage = numimage + 1
               imagcrds(1,numimage) = f1*ucell(1,1) + f2*ucell(1,2) + &
                  f3*ucell(1,3)
               imagcrds(2,numimage) = f1*ucell(2,1) + f2*ucell(2,2) + &
                  f3*ucell(2,3)
               imagcrds(3,numimage) = f1*ucell(3,1) + f2*ucell(3,2) + &
                  f3*ucell(3,3)
               bckptr(numimage) = j
            end do
      end if  
   end do 
   if(master) then
      if (verbose == 1)write(6,*)'grid_image: numimage = ',numimage
   end if
   return
end subroutine grid_image 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                GRID POINTERS
!-----------------------------------------------------------------------
!+ list of nearby cell line centers that need to be searched for nb list
!  example:   2d grid of cells, X is cell of interest.Here we are using 3
!             cells to the cutoff, so need 3 cells distant in each direction
!             We only point to the center of each string of 7 cells in the
!             x direction, and the nb list routine knows to check this cell
!             and three on each side of it in the +x and -x directions.
!              This routine will list the center of cell strings:
!         A                          B
!    ..................    ..................   and three more like layer B
!    ..................    ..................
!    .....ooo#ooo......    .....ooo#ooo......
!    .....ooo#ooo......    .....ooo#ooo......
!    .....ooo#ooo......    .....ooo#ooo......
!    ........Xooo......    .....ooo#ooo......
!    ..................    .....ooo#ooo......
!    ..................    .....ooo#ooo......
!    ..................    .....ooo#ooo......
!    ..................    ..................
!         cell X and its      Next layer ahead
!         neighbors o
!         center cell #
!         (ahead only search)
!
!    The pointer list contains the identity of all the # cells and the X cell
!
subroutine grid_pointers(nghbptr,nghtran,periodic,nogrdptrs)
   implicit none
   integer nghbptr(0:maxnptrs,*),nghtran(maxnptrs,*)
   integer index0
   integer periodic

   integer i1,i2,i3,j2,j3,i2grd,i3grd,j3grd,j2grd,jj2
   integer j2strt,j2stop,jj2strt,jj2stop,j3stop
   integer index,num
   integer k3off
   integer xtindex_1,xtindex_2
   integer sizgrd12
   integer index0_min,index0_max

   logical nogrdptrs

   index0_min = maxnptrs
   index0_max = 0
   sizgrd12=nucgrd1*nucgrd2
   my_grids=0               ! whole array zeroed

   !     ---get forward pointers for unit cell grid
   !     neighbor cells must be later in the subcell list
   !     to get unique atom - imageatom pairing
   !-------------------------------------------------------------------------
   !          Non-Periodic cell list pointers
   !            
   if(periodic == 0)then
      if(nogrdptrs)then
         do i1=myindexlo,nucgrd1*nucgrd2*nucgrd3
            my_grids(i1)=1
         end do
      else
         do i3 = 1,nucgrd3
            i3grd=sizgrd12*(i3-1)
            do i2 = 1,nucgrd2
               i2grd=i3grd+nucgrd1*(i2-1)
               do i1 = 1,nucgrd1
                  index = i2grd+i1
                  if(index < myindexlo.or.index > myindexhi)goto 180
                  my_grids(index)=1
                  index0=index-myindexlo+1
                  index0_max = max (index0_max,index0)
                  index0_min = min (index0_min,index0)
                  !
                  !     Set up xtran index for this cell
                  !
                  !     UNSTABLE till fixed:
                  !       Need to make sure that atoms do not drift to edges
                  !       Otherwise this simplified xtindex will not work.
                  xtindex_1 = 1
                  xtindex_2 = 1
                  j2strt=-nghb2
                  if(i2+j2strt < 1)j2strt=1-i2
                  j2stop=nghb2
                  if(i2+j2stop > nucgrd2)j2stop=nucgrd2-i2
                  jj2strt=-nghb1
                  if(i1+jj2strt < 1)jj2strt=1-i1
                  jj2stop=nghb1
                  if(i1+jj2stop > nucgrd1)jj2stop=nucgrd1-i1
                  !
                  ! first do the plane that this cell is in.
                  ! first row of neighbors starts with the cell itself (Cell0)
                  !
                  num = 1
                  nghbptr(num,index0) = index
                  !
                  !     no y or z translate, this cell HAS to be in uc
                  !
                  nghtran(num,index0) = 5+ishft(1,8)
                  do jj2=1,jj2stop
                     my_grids(index+jj2)=1
                  end do
                  !
                  !next rows of neighbors are the nghb2 rows of 2*nghb1+1 cells
                  !  above.  Use the center cell of the row as the pointer
                  !  (with the same x position as Cell0)
                  !This way the pointer cell must be in the uc wrt the x 
                  !   direction.
                  do j2 = 1,j2stop
                     j2grd=index+j2*nucgrd1
                     if(j2+i2 <= nucgrd2)then
                        num=num+1
                        nghtran(num,index0)=5+ishft(xtindex_2,8)
                        nghbptr(num,index0) = j2grd
                        do jj2=jj2strt,jj2stop
                           my_grids(j2grd+jj2 &
                                 -nucgrd1*xtran(jj2+1+nghb1,xtindex_2))=1
                        end do
                     end if
                  end do
                  !
                  ! Now there are nghb3 planes of (2*nghb1+1)(2*nghb2+1) cells
                  !     ahead that are good neighbors.
                  !
                  j3stop=nghb3
                  if(i3+j3stop > nucgrd3)j3stop=nucgrd3-i3
                  do j3 = 1,j3stop
                     j3grd=j3*sizgrd12
                     do j2 = j2strt,j2stop
                        num=num+1
                        j2grd=index+j3grd+j2*nucgrd1
                        nghtran(num,index0)= &
                              5+ishft(xtindex_2,8)
                        nghbptr(num,index0) = j2grd
                        do jj2=jj2strt,jj2stop
                           my_grids(j2grd+jj2)=1
                        end do
                     end do
                  end do
                  nghbptr(0,index0) = num
                  180 continue
               end do  
            end do  
         end do 
      end if  ! (nogrdptrs)
      numnptrs = maxnptrs
      return
   end if 
   !-------------------------------------------------------------------------
   do i3 = 1,nucgrd3
      i3grd=sizgrd12*(i3-1)
      do i2 = 1,nucgrd2
         i2grd=i3grd+nucgrd1*(i2-1)
         do i1 = 1,nucgrd1
            index = i2grd+i1
            if( (index >= myindexlo) .and. (index <= myindexhi) ) then
               my_grids(index)=1
               index0=index-myindexlo+1
               index0_max = max (index0_max,index0)
               index0_min = min (index0_min,index0)
            !
            !  Set up xtran index for this cell
            !
               if(i1+nghb1 > nucgrd1)then
                  xtindex_1 = 1+i1+nghb1-nucgrd1
                  xtindex_2 = xtindex_1 + 2*nghb1
               else if(i1-nghb1 < 1)then
                  xtindex_1 = 1
                  xtindex_2 = 2*nghb1+2-i1
               else
                  xtindex_1 = 1
                  xtindex_2 = 1
               end if
            !
            !  first do the plane that this cell is in.
            !  first row of neighbors starts with the cell itself (Cell0)
            !
               num = 1
               nghbptr(num,index0) = index
            !
            !   no y or z translate, this cell HAS to be in uc
            !
               nghtran(num,index0) = 5+ishft(xtindex_1,8)
               do jj2=1,nghb1
                  my_grids(index+jj2 &
                     -nucgrd1*xtran(jj2+1,xtindex_1))=1
               end do
            !
            ! next rows of neighbors are the nghb2 rows of 2*nghb1+1 cells
            !   above.  Use the center cell of the row as the pointer
            !   (with the same x position as Cell0)
            !  This way the pointer cell must be in the uc wrt the x direction.
            !
               do j2 = 1,nghb2
                  num=num+1
                  j2grd=index+j2*nucgrd1
                  if(j2+i2 <= nucgrd2)then
                     nghtran(num,index0)=5+ishft(xtindex_2,8)
                  else
                     j2grd=j2grd-sizgrd12
                     nghtran(num,index0)=8+ishft(xtindex_2,8)
                  end if
                  nghbptr(num,index0) = j2grd
                  do jj2=-nghb1,nghb1
                     my_grids(j2grd+jj2 &
                        -nucgrd1*xtran(jj2+1+nghb1,xtindex_2))=1
                  end do
               end do
            !
            ! Now there are nghb3 planes of (2*nghb1+1)(2*nghb2+1) cells
            !     ahead that are good neighbors.
            !
               do j3 = 1,nghb3
                  j3grd=j3*sizgrd12
                  if(j3+i3 > nucgrd3)then
                     j3grd=j3grd-sizgrd12*nucgrd3
                     k3off=9
                  else
                     k3off=0
                  end if
                  do j2 = -nghb2,nghb2
                     num=num+1
                     j2grd=index+j3grd+j2*nucgrd1
                     if(j2+i2 > nucgrd2)then
                        nghtran(num,index0)=8+k3off+ishft(xtindex_2,8)
                        j2grd=j2grd-sizgrd12
                     else if(j2+i2 < 1)then
                        nghtran(num,index0)=2+k3off+ishft(xtindex_2,8)
                        j2grd=j2grd+sizgrd12
                     else
                        nghtran(num,index0)=5+k3off+ishft(xtindex_2,8)
                     end if
                     nghbptr(num,index0) = j2grd
                     do jj2=-nghb1,nghb1
                        my_grids(j2grd+jj2 &
                            -nucgrd1*xtran(jj2+1+nghb1,xtindex_2))=1
                     end do
                  end do
               end do
            endif ! index within lo and hi
         end do   
      end do 
   end do

   numnptrs = num

   return
end subroutine grid_pointers 
!-------------------------------------------------------------------
!     --- SETUP_GRIDS ---
!     ...this routine checks on the necessary resources for the unit
!     cell and image cell grids used for short range particle pair
!     calculations.  It is assumed that unit cell setup has already
!     occurred.  This routine then checks the short range cutoff.
!     The unit cell will be split into NUCGRD1 x NUCGRD2 x NUCGRD3
!     geometrically similar subcells of size dirlng(1)/NUCGRD1  by
!     dirlng(2)/NUCGRD2  by  dirlng(3)/NUCGRD3.
!     The short range interactions will involve pairs in the subcell
!     neighborhood of  +- NGHB1  by  +- NGHB2  by  +- NGHB3  subcells.,
!     about any given subcell.
!     The distances between  parallel faces of the unit cell are ,
!     respectively reclng(1), reclng(2) and reclng(3).
!     Thus these subcell neighborhoods are guaranteed to contain all points
!     within the minimum of (NGHB1/NUCGRD1)*reclng(1),
!     (NGHB2/NUCGRD2)*reclng(2),and (NGHB3/NUCGRD3)*1.0/reclng(3).
!     This minimum is taken to be the short range cutoff...


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine setup_grids here]
subroutine setup_grids(periodic,nogrdptrs,verbose)
   implicit none
#  include "extra.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#else
   integer mytaskid,numtasks
   parameter (mytaskid=0,numtasks=1)
#endif
   _REAL_ sizmaxhb
   integer  periodic,verbose
#ifdef MPI
   integer i,indtop,inddel
#endif
   integer nghb0
   integer ngrd1,ngrd2,ngrd3
   logical             :: nogrdptrs

   _REAL_ dc1,dc2,dc3,cut

   parameter(sizmaxhb=1.34d0)    ! max bond length for h to heavy atom

   nghb0=nghb
   nghb1 = nghb0
   nghb2 = nghb0
   nghb3 = nghb0

   dc1 = cutlist / nghb1
   dc2 = cutlist / nghb2
   dc3 = cutlist / nghb3
   nucgrd1 = max(1,int(reclng(1) / dc1) )
   nucgrd2 = max(1,int(reclng(2) / dc2) )
   nucgrd3 = max(1,int(reclng(3) / dc3) )
   if(periodic == 0)then
      !       ----- If system has too few cells for the pointer method
      !       -----   to be efficient, use the no-grid-pointer system
      !       -----   where all forward cells are checked for pair atoms
      !       -----   rather than just the forward cells on the grid-pointer-list
      
      nogrdptrs=nogrdptrs .or.( (nucgrd1 <= 2*nghb+2) .or. &
            (nucgrd2 <= 2*nghb+2) .or. &
            (nucgrd3 <= 2*nghb+2) )
      if(verbose >= 3)write(6,*)" List using nogrdptrs ",nogrdptrs
      !       ----- Make the subcells about 3 A in size now so that there
      !       ----- are more subcells, easier to load balance.
      if(nogrdptrs)then
         if(dirlng(1)/nucgrd1 > 3)then
            ngrd1 = dirlng(1)/3
            ngrd2 = dirlng(2)/3
            ngrd3 = dirlng(3)/3
            if ( verbose >= 1)then
               write(6,'("|    New Grids set up for nogrdptrs ")')
               write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)') &
                     'Number of grids per unit cell x,y,z:', &
                     ngrd1, ngrd2, ngrd3
            end if
            nucgrd1=ngrd1
            nucgrd2=ngrd2
            nucgrd3=ngrd3
         end if
      end if
   end if
   
   ! check the short range cutoff:
   
   dc1 = reclng(1)/nucgrd1
   dc2 = reclng(2)/nucgrd2
   dc3 = reclng(3)/nucgrd3
   cut = nghb1*dc1
   if ( nghb2*dc2 < cut )cut = nghb2*dc2
   if ( nghb3*dc3 < cut )cut = nghb3*dc3
   if(nogrdptrs) cut=cutlist
#ifdef MPI
   if(master) then
#endif
      if ( verbose >= 1)then
         write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)') &
               'Number of grids per unit cell in each dimension:', &
               nucgrd1, nucgrd2, nucgrd3
         write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
               'Unit cell edge lengths in each dimension:', &
               dirlng(1), dirlng(2), dirlng(3)
         write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
               'Distance between parallel faces of unit cell:', &
               reclng(1), reclng(2), reclng(3)
         write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
               'Distance between faces of short range grid subcells:', &
               dc1, dc2, dc3
         write(6, '(5X,a,F9.3)') &
               'Resulting cutoff from subcell neighborhoods is ', cut
      end if
#ifdef MPI
   end if
#endif
   if ( cut < cutlist )then
      call sander_bomb('setup_grids', &
            'Resulting cutoff is too small for your lower limit',' ')
   end if
   return
end subroutine setup_grids 

subroutine setup_grid_sizes(periodic,nogrdptrs,verbose)

   implicit none
#  include "extra.h"
#ifdef MPI
#  include "ew_parallel.h"
#  include "parallel.h"
#endif
   _REAL_ sizmaxhb
   integer periodic,verbose
#ifdef MPI
   integer myindexlo,myindexhi
   integer i,indtop,inddel
#endif
   integer nghb0
   integer ngrd1,ngrd2,ngrd3
   logical nogrdptrs

   _REAL_ dc1,dc2,dc3,cut

   parameter(sizmaxhb=1.34d0)    ! max bond length for h to heavy atom

   nghb0=nghb
   nghb1 = nghb0
   nghb2 = nghb0
   nghb3 = nghb0

   dc1 = (cutlist + sizmaxhb) / nghb1
   dc2 = (cutlist + sizmaxhb) / nghb2
   dc3 = (cutlist + sizmaxhb) / nghb3
   nucgrd1 = max(1,int(reclng(1) / dc1) )
   nucgrd2 = max(1,int(reclng(2) / dc2) )
   nucgrd3 = max(1,int(reclng(3) / dc3) )
   if(periodic == 0)then
      !       ----- If system has too few cells for the pointer method
      !       -----   to be efficient, use the no-grid-pointer system
      !       -----   where all forward cells are checked for pair atoms
      !       -----   rather than just the forward cells on the grid-pointer-list
      
      nogrdptrs=nogrdptrs .or.( (nucgrd1 <= 2*nghb+2) .or. &
            (nucgrd1 <= 2*nghb+2) .or. &
            (nucgrd1 <= 2*nghb+2) )
      if(verbose >= 3)write(6,*)" List using nogrdptrs ",nogrdptrs
      !       ----- Make the subcells about 3 A in size now so that there
      !       ----- are more subcells, easier to load balance.
      if(nogrdptrs)then
         if(dirlng(1)/nucgrd1 > 3)then
            ngrd1 = dirlng(1)/3
            ngrd2 = dirlng(2)/3
            ngrd3 = dirlng(3)/3
            if ( verbose >= 1)then
               write(6,'("|    New Grids set up for nogrdptrs ")')
               write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)') &
                     'Number of grids per unit cell x,y,z:', &
                     ngrd1, ngrd2, ngrd3
            end if
            nucgrd1=ngrd1
            nucgrd2=ngrd2
            nucgrd3=ngrd3
         end if
      end if
   end if
   
   ! check the short range cutoff:
   
   dc1 = reclng(1)/nucgrd1
   dc2 = reclng(2)/nucgrd2
   dc3 = reclng(3)/nucgrd3
   cut = nghb1*dc1
   if ( nghb2*dc2 < cut )cut = nghb2*dc2
   if ( nghb3*dc3 < cut )cut = nghb3*dc3
   if(nogrdptrs) cut=cutlist
#ifdef MPI
   if(master) then
#endif
      if ( verbose >= 1)then
         write(6, '(5X,a,/,5X,i9,1X,i9,1X,i9)') &
               'Number of grids per unit cell in each dimension:', &
               nucgrd1, nucgrd2, nucgrd3
         write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
               'Unit cell edge lengths in each dimension:', &
               dirlng(1), dirlng(2), dirlng(3)
         write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
               'Distance between parallel faces of unit cell:', &
               reclng(1), reclng(2), reclng(3)
         write(6, '(5X,a,/,5X,F9.3,1X,F9.3,1X,F9.3)') &
               'Distance between faces of short range grid subcells:', &
               dc1, dc2, dc3
         write(6, '(5X,a,F9.3)') &
               'Resulting cutoff from subcell neighborhoods is ', cut
      end if
#ifdef MPI
   end if
#endif
   if ( cut < cutlist )then
      call sander_bomb('setup_grids', &
            'Resulting cutoff is too small for your lower limit',' ')
   end if
   return
end subroutine setup_grid_sizes

!-------------------------------------------------------------------
!     --- GRID_UCELL ---

!     ...this routine grids the mapped atoms in the unit cell
!     into the NUCGRD1 x NUCGRD2 x NUCGRD3 subcells according
!     to their fractional coordinates.

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine grid_ucell here]
subroutine grid_ucell(natom,periodic)
   use constants, only : half
   implicit none
   integer natom,periodic
   integer i,j,i1,i2,i3,index
!   integer ichk
   _REAL_ shft

   shft=half
   if(periodic == 0)shft=0.0

   do index = 1,nucgrd
      numatg(index) = 0
   end do
   !-----------------------------------------------------------------
   !     ---find out which ucgrd subcell each atom is in.
   !          atmcell(i) is the ucgrd subcell that contains atom i
   !          numatg(index) is the number of atoms in the index ucgrd
   do i = 1,natom
         i1 = (fraction(1,i) + shft) * nucgrd1 + 1
         i2 = (fraction(2,i) + shft) * nucgrd2 + 1
         i3 = (fraction(3,i) + shft) * nucgrd3 + 1
         index = nucgrd1*nucgrd2*(i3-1)+nucgrd1*(i2-1)+i1
         atmcell(i) = index
         numatg(index) = numatg(index) + 1
   end do
   !   
   !     ---find the offset of the starting atoms for each ucgrd subcell.
   !         zero the numatg()entries as you go
   indoff(1) = 0
   do i = 2,nucgrd
      indoff(i) = indoff(i-1) + numatg(i-1)
      numatg(i-1) = 0
   end do
   
   ! ichk = natom - indoff(nucgrd) - numatg(nucgrd)
   ! this crash also only seems to occur when atoms can not be gridded due to
   ! truly awful coordinates e.g. NaN (We get with really bad input data
   ! sometimes)
   ! if ( ichk .ne. 0 )then
   !   write(6,*)'BIG PROBLEM in grid_ucell!!!'
   !   call mexit(6,1)
   ! endif
   !-----------------------------------------------------------------
   !  Fill indatg() as a list of atoms in each subcell such that
   !       the list of atoms in subcell 1 are at the beginning,
   !       subcell 2 list is right after that (starting at indoff(2)+1)
   !
   numatg(nucgrd) = 0
   do i = 1,natom
         index = atmcell(i)
         numatg(index) = numatg(index) + 1
         j = numatg(index)+indoff(index)
         indatg(j) = i
   end do
   return
end subroutine grid_ucell 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!             PACK_NB_LIST
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
subroutine pack_nb_list(kk,i,xk,yk,zk,imagcrds,cutoffsq, &
      numlist,numpack, &
      iac,ico,ntypes, &
      ipairs, &
      ifail, &
      tranvec, &
      belly,ibelly,cn1)
#ifdef LES
   use pimd_vars, only: ipimd
#endif
#ifdef MPI /* SOFT CORE */
   use softcore, only : nsc, ifsc
#endif
   use qmmm_module, only : qmmm_nml,qmmm_struct
   use constants, only : zero
   implicit none
   integer numpack
   integer kk,i,numlist, &
         iac(*),ico(*),ntypes, &
         ifail 
   _REAL_ xk,yk,zk,imagcrds(3,*),cutoffsq,cn1(*)
   integer ipairs(*),ibelly(*)
   logical belly,deadi,deadik
   integer num
!Local variables for QMMM
   integer qm_temp_count2

#ifdef MPI
#  include "parallel.h"
#endif

#ifdef LES
#  include "les.h"
#endif

   _REAL_ dx,dy,dz,r2
   integer iaci,ic,index,k,lpr,lps,lhb,m,n,npr

   integer jtran
   _REAL_ x_tran,y_tran,z_tran,tranvec(3,*)
#ifdef CRAY_X1
   integer m2,mm,mlist(numlist)
#endif

   num = 0
   m = maskptr(i)
   lpr = 0
   lps = 0
   lhb = 0
   iaci = ntypes*(iac(i)-1)
   
   do n = 1,nummask(i)
      k = lstmask(m+n)
      exclude(k) = i
   end do

#ifdef LES
   if(ipimd>0) then
      do m = kk+1,numlist
         n = atmlist(m)
         k = bckptr(n)
         if( lestyp(i).eq.lestyp(k).and.cnum(i).ne.cnum(k) ) then
            exclude(k)=i
         end if
      end do
   endif
#endif
 
        

   if ( qmmm_nml%ifqnt ) then
     if ( qmmm_struct%atom_mask(i) ) then ! then current atom is a QM atom
       ! skip interaction with all other QM atoms:
       do qm_temp_count2=1, qmmm_struct%nquant
         exclude(qmmm_struct%iqmatoms(qm_temp_count2))=i
       end do
     end if
   end if

   deadi=.false.
   deadi=(ibelly(i) == 0).and.belly

#ifdef CRAY_X1
   m2=0
   do m = kk+1,numlist
      jtran=itran(m)
      x_tran=tranvec(1,itran(m))
      y_tran=tranvec(2,itran(m))
      z_tran=tranvec(3,itran(m))
      n = atmlist(m)
      k = bckptr(n)
      deadik=(ibelly(k) == 0).and.deadi
      dx = imagcrds(1,n) - xk + x_tran
      dy = imagcrds(2,n) - yk + y_tran
      dz = imagcrds(3,n) - zk + z_tran
      r2 = dx*dx + dy*dy + dz*dz
      if ( r2 < cutoffsq ) then
#ifdef TEST_CLOSEATOMS
         if(r2 < .5) write(6,190)i,k,r2
 190     format("<pack_nb_list> Atoms close:",2i8,"  r**2 ",f10.6)
#endif
         if ( (exclude(k) /= i) .and. .not.deadik) then
            mygrdlist(n)=ior(mygrdlist(n),1)
            num = num + 1
            m2=m2+1
            mlist(m2)=m
         endif
      endif
   enddo
   do mm=1,m2
      m=mlist(mm)
      jtran=itran(m)
      x_tran=tranvec(1,itran(m))
      y_tran=tranvec(2,itran(m))
      z_tran=tranvec(3,itran(m))
      n = atmlist(m)
      k = bckptr(n)
      index = iaci+iac(k)
      ic = ico(index)
      if (ic >= 0) then
         lpr = lpr+1
         iwa(lpr) = ior(n,ishft(itran(m),27)) ! bitwise optimization
      else
         lhb = lhb+1
         iwh(lhb) = ior(n,ishft(itran(m),27)) ! bitwise optimization
      end if
   enddo

#else /* Generic (Not cray_x1) follows */

#ifdef MPI /* SOFT CORE */
   softcore_on: if(ifsc /= 0) then ! softcore potential in use, check for each atom on which list it goes
      check_softcore: if (nsc(i) == 0) then ! atom i is not a softcore atom
         do m = kk+1,numlist            
            
            jtran=itran(m)
            x_tran=tranvec(1,itran(m))
            y_tran=tranvec(2,itran(m))
            z_tran=tranvec(3,itran(m))
            n = atmlist(m)
            k = bckptr(n)
            deadik=(ibelly(k) == 0).and.deadi
            if ( (exclude(k) /= i) .and. .not.deadik) then
               dx = imagcrds(1,n) - xk + x_tran
               dy = imagcrds(2,n) - yk + y_tran
               dz = imagcrds(3,n) - zk + z_tran
               r2 = dx*dx + dy*dy + dz*dz
               if ( r2 < cutoffsq ) then
#  ifdef TEST_CLOSEATOMS
                  if(r2 < .5) &
                       write(6,190)i,k,r2
190               format("<pack_nb_list> Atoms close:",2i8, &
                       "  r**2 ",f10.6)
#  endif
                  mygrdlist(n)=ior(mygrdlist(n),1)
                  num = num + 1
                  
                  index = iaci+iac(k)
                  ic = ico(index)
                  if (ic > 0) then
                     if (nsc(k) == 0) then
                        lpr = lpr+1
                        iwa(lpr) = ior(n,ishft(itran(m),27)) ! bitwise optimization
                     else
                        lps = lps+1
                        iws(lps) = ior(n,ishft(itran(m),27)) ! bitwise optimization
                     end if
                  else
                     if (nsc(k) == 0) then
                        lhb = lhb+1
                        iwh(lhb) = ior(n,ishft(itran(m),27)) ! bitwise optimization
                     else
                        lps = lps+1
                        iws(lps) = ior(n,ishft(itran(m),27)) ! bitwise optimization
                     end if
                  end if
               
               end if
            end if
         end do
      else check_softcore ! atom i IS a softcore atom
         do m = kk+1,numlist            
         
            jtran=itran(m)
            x_tran=tranvec(1,itran(m))
            y_tran=tranvec(2,itran(m))
            z_tran=tranvec(3,itran(m))
            n = atmlist(m)
            k = bckptr(n)
            deadik=(ibelly(k) == 0).and.deadi
            if ( (exclude(k) /= i) .and. .not.deadik) then
               dx = imagcrds(1,n) - xk + x_tran
               dy = imagcrds(2,n) - yk + y_tran
               dz = imagcrds(3,n) - zk + z_tran
               r2 = dx*dx + dy*dy + dz*dz
               if ( r2 < cutoffsq ) then
                  mygrdlist(n)=ior(mygrdlist(n),1)
                  num = num + 1
                  index = iaci+iac(k)
                  ic = ico(index)
                  lps = lps+1
                  iws(lps) = ior(n,ishft(itran(m),27)) ! bitwise optimization
               end if
            end if
         end do
      end if check_softcore
   else softcore_on
# endif
      ! softcore potential not on, build list the normal way
      do m = kk+1,numlist            
            
         jtran=itran(m)
         x_tran=tranvec(1,itran(m))
         y_tran=tranvec(2,itran(m))
         z_tran=tranvec(3,itran(m))
         n = atmlist(m)
         k = bckptr(n)
         deadik=(ibelly(k) == 0).and.deadi
         if ( (exclude(k) /= i) .and. .not.deadik) then
            dx = imagcrds(1,n) - xk + x_tran
            dy = imagcrds(2,n) - yk + y_tran
            dz = imagcrds(3,n) - zk + z_tran
            r2 = dx*dx + dy*dy + dz*dz
            if ( r2 < cutoffsq ) then
#  ifdef TEST_CLOSEATOMS
               if(r2 < .5) &
                    write(6,190)i,k,r2
190            format("<pack_nb_list> Atoms close:",2i8, &
                    "  r**2 ",f10.6)
#  endif
               mygrdlist(n)=ior(mygrdlist(n),1)
               num = num + 1
               
               index = iaci+iac(k)
               ic = ico(index)
               if (ic >= 0) then
                  lpr = lpr+1
                  iwa(lpr) = ior(n,ishft(itran(m),27)) ! bitwise optimization
               else
                  lhb = lhb+1
                  iwh(lhb) = ior(n,ishft(itran(m),27)) ! bitwise optimization
               end if
            end if
         end if
      end do
# ifdef MPI /* SOFT CORE */
   end if softcore_on
# endif 
#endif /* cray X1 */
   if ( num + numpack > maxnblst )then
#ifdef MPI
      write(6, '(/,a,i12,i12,a,i12,a,i4)') &
            ' * NB pairs ', num,numpack, &
            ' exceeds capacity (',maxnblst,')', mytaskid
#else
      write(6, '(/,a,i12,i12,a,i12,a)') &
            ' * NB pairs ', num,numpack, &
            ' exceeds capacity (',maxnblst,')'
#endif
      ifail=1
      return
   end if

   !     ---now put all pairs into iwa

   do k = 1,lhb
      iwa(lpr+k) = iwh(k)
   end do
   npr = lpr+lhb

#ifdef MPI /* SOFT CORE */
   iwa(npr+1:npr+lps) = iws(1:lps)
   npr = npr + lps
   numsc(i) = lps
#endif

   numvdw(i) = lpr
   numhbnd(i) = lhb
   do k = 1,npr
      numpack = numpack+1
      ipairs(numpack) = iwa(k)
   end do
   return
end subroutine pack_nb_list 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pack_nb_nogrdptrs here]
subroutine pack_nb_nogrdptrs(kk,i,xk,yk,zk,imagcrds,cutoffsq, &
      numlist,numpack, &
      iac,ico,ntypes, &
      ipairs, &
      ifail, &
      belly,ibelly)
   use qmmm_module, only : qmmm_nml,qmmm_struct
#ifdef MPI /* SOFT CORE */
   use softcore, only : nsc, ifsc
#endif
   implicit none
   integer numpack
   integer kk,i,numlist, &
         iac(*),ico(*),ntypes,ifail 
   _REAL_ xk,yk,zk,imagcrds(3,*),cutoffsq
   integer ipairs(*),ibelly(*)
   logical belly,deadi,deadik
   integer jtran,num

!Local variables for QMMM
   integer qm_temp_count2

#ifdef MPI
#  include "parallel.h"
#endif

   _REAL_ dx,dy,dz,r2
   integer iaci,ic,index,k,lpr,lps,lhb,m,n,npr


   num = 0
   m = maskptr(i)
   lpr = 0
   lps = 0
   lhb = 0
   iaci = ntypes*(iac(i)-1)

   do n = 1,nummask(i)
      k = lstmask(m+n)
      exclude(k) = i
   end do


 

!QMMM 
   if ( qmmm_nml%ifqnt ) then
      if ( qmmm_struct%atom_mask(i) ) then ! then current atom is a QM atom
        do qm_temp_count2=1, qmmm_struct%nquant
          exclude(qmmm_struct%iqmatoms(qm_temp_count2))=i ! skip interaction with all other QM atoms
        end do
      end if
   end if
!END QMMM


   deadi=.false.
   deadi=(ibelly(i) == 0).and. belly
#ifdef MPI /* SOFT CORE */
   softcore_on: if (ifsc /= 0) then
      check_softcore: if (nsc(i) == 0) then ! atom i is NOT a softcore atom
         do m = kk+1,numlist
            jtran=5
            n = atmlist(m)
            k = bckptr(n)
            deadik=(ibelly(k) == 0).and.deadi
            if ( (exclude(k) /= i) .and. .not.deadik) then
               dx = imagcrds(1,n) - xk
               dy = imagcrds(2,n) - yk
               dz = imagcrds(3,n) - zk
               r2 = dx*dx + dy*dy + dz*dz
               if ( r2 < cutoffsq ) then
# ifdef TEST_CLOSEATOMS
                  if(r2 < .5) &
                       write(6,190)i,k,r2
190               format("<pack_nb_list> Atoms close:",2i8, &
                       "  r**2 ",f10.6)
                  
# endif
                  if (mygrdlist(n) == 0) then
                     mygrdlist(n)=1
                  end if
                  num = num + 1

                  index = iaci+iac(k)
                  ic = ico(index)
                  if(ic >= 0) then
                     if (nsc(k) == 0) then
                        lpr = lpr+1
                        iwa(lpr) = ior(n,ishft(5,27)) ! bitwise optimization
                     else
                        lps = lps+1
                        iws(lps) = ior(n,ishft(itran(m),27)) ! bitwise optimization
                     end if
                  else
                     lhb = lhb+1
                     iwh(lhb) = ior(n,ishft(5,27)) ! bitwise optimization
                  end if
               end if
            end if
         end do  
      else check_softcore ! atom i IS a softcore atom
         do m = kk+1,numlist
            jtran=5
            n = atmlist(m)
            k = bckptr(n)
            deadik=(ibelly(k) == 0).and.deadi
            if ( (exclude(k) /= i) .and. .not.deadik) then
               dx = imagcrds(1,n) - xk
               dy = imagcrds(2,n) - yk
               dz = imagcrds(3,n) - zk
               r2 = dx*dx + dy*dy + dz*dz
               if ( r2 < cutoffsq ) then
                  if (mygrdlist(n) == 0) then
                     mygrdlist(n)=1
                  end if
                  num = num + 1

                  index = iaci+iac(k)
                  ic = ico(index)
                  if(ic >= 0) then
                     lps = lps+1
                     iws(lps) = ior(n,ishft(5,27)) ! bitwise optimization
                  else ! unlikely case atom is softcore and 10-12 flagged
                     lhb = lhb+1
                     iwh(lhb) = ior(n,ishft(5,27)) ! bitwise optimization
                  end if
               end if
            end if
         end do
      end if check_softcore

   else softcore_on ! softcore is not on, build list the usual way
#endif /* SOFT CORE  */
      do m = kk+1,numlist
         jtran=5
         n = atmlist(m)
         k = bckptr(n)
         deadik=(ibelly(k) == 0).and.deadi
         if ( (exclude(k) /= i) .and. .not.deadik) then
            dx = imagcrds(1,n) - xk
            dy = imagcrds(2,n) - yk
            dz = imagcrds(3,n) - zk
            r2 = dx*dx + dy*dy + dz*dz
            if ( r2 < cutoffsq ) then
# ifdef TEST_CLOSEATOMS
               if(r2 < .5) &
                    write(6,190)i,k,r2
190            format("<pack_nb_list> Atoms close:",2i8, &
                    "  r**2 ",f10.6)
# endif
               if (mygrdlist(n) == 0) then
                  mygrdlist(n)=1
               end if
               num = num + 1

               index = iaci+iac(k)
               ic = ico(index)
               if(ic >= 0) then
                  lpr = lpr+1
                  iwa(lpr) = ior(n,ishft(5,27)) ! bitwise optimization
               else
                  lhb = lhb+1
                  iwh(lhb) = ior(n,ishft(5,27)) ! bitwise optimization
               end if

            end if
         end if
      end do
#ifdef MPI /* SOFT CORE */
   end if softcore_on
#endif

   if ( num + numpack > maxnblst )then
#ifdef MPI
      write(6, '(/,a,i12,i12,a,i12,a,i4)') &
            ' * NB pairs ', num,numpack, &
            ' exceeds capacity (',maxnblst,')', mytaskid
#else
      write(6, '(/,a,i12,i12,a,i12,a)') &
            ' * NB pairs ', num,numpack, &
            ' exceeds capacity (',maxnblst,')'
#endif
      ifail=1
      return
   end if

   !     ---now put all pairs into iwa

   do k = 1,lhb
      iwa(lpr+k) = iwh(k)
   end do
   npr = lpr+lhb

#ifdef MPI /* SOFT CORE */
   iwa(npr+1:npr+lps) = iws(1:lps)
   npr = npr + lps
   numsc(i) = lps
#endif

   numvdw(i) = lpr
   numhbnd(i) = lhb
   do k = 1,npr
      numpack = numpack+1
      ipairs(numpack) = iwa(k)
   end do
   return
end subroutine pack_nb_nogrdptrs 


!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Re-balance the grid by divvying cells among processors.
!
!-----------------------------------------------------------------------
!     --- FIX_GRID_BALANCE ---
!-----------------------------------------------------------------------

subroutine fix_grid_balance(gridpairs, &
      nucgrd,numtasks,mytaskid,listdiffmax)
   implicit none

   integer gridpairs(*)
   integer nucgrd,listtot,numtasks,mytaskid,listdiffmax
   integer i,n,lsum,i0,ishare

   n=0
   lsum=0
   myindexlo=0
   myindexhi=0
   listtot=0
   do i=1,nucgrd
      listtot=listtot+gridpairs(i)
   end do

   i0=1
   ishare=listtot/numtasks
   
   !   ------ set trigger for new balance as 10 percent of list size
   
   listdiffmax=max(1000,ishare/10)
   
   !   ------ give cells to each processor .lt.mytaskid till
   !          they have ishare, then give the next
   !          cells to this pe till it has its share, then return
   
   do i=1,nucgrd
      lsum=lsum+gridpairs(i)
      if(lsum >= ishare)then
         if(n == mytaskid)then
            myindexlo=i0
            myindexhi=i
            return
         end if
         lsum=0
         i0=i+1
         n=n+1
      end if
   end do
   if(n == mytaskid)then
      myindexlo=i0
      myindexhi=nucgrd
      return
   end if
   myindexlo=nucgrd
   return
end subroutine fix_grid_balance 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Determine imaged coordinates from the unit cell and fractional crds.
!-----------------------------------------------------------------------
!     --- ADJUST_IMAGCRDS ---
!-----------------------------------------------------------------------
!     ...this is needed in case you are wrapping coords in the box.
!     The code can be run without wrapping (since it is easier to
!     analyze the results) but it should work equally well,
!     either way...

!---------------------------------------------------------------------

subroutine adjust_imagcrds(crd,natom)

   use trace
   use constants, only : half
   implicit none
   integer, intent(in) :: natom
   _REAL_ , dimension(3,natom), intent(in) :: crd

   integer i,j
   _REAL_ f1,f2,f3
   _REAL_ anint
#  include "parallel.h"
#  include "ew_cntrl.h"
#  include "box.h"
   !   ---- For nonperiodic, there is no imaging so frac coords are not
   !   ---- used and no wrapping is done, but the imgcrds need to be
   !   ---- in the right order
   if(periodic == 0)then
      do i=1,natom
         if (mygrdlist(i) == 1) then
           j = bckptr(i)
           imagcrds(1,i) = crd(1,j)+xbox0
           imagcrds(2,i) = crd(2,j)+ybox0
           imagcrds(3,i) = crd(3,j)+zbox0
         end if
      end do
      return
   end if
   
   !   ---- Periodic systems
   !        find change since last step:
   
   i=0
   do while ( i < natom )
      i=i+1
      if(mygrdlist(i) == 1)then
         j = bckptr(i)
         dfrac(1,j) = fraction(1,j) - savfrac(1,j)
         dfrac(2,j) = fraction(2,j) - savfrac(2,j)
         dfrac(3,j) = fraction(3,j) - savfrac(3,j)
         dfrac(1,j) = dfrac(1,j) - anint(dfrac(1,j))
         dfrac(2,j) = dfrac(2,j) - anint(dfrac(2,j))
         dfrac(3,j) = dfrac(3,j) - anint(dfrac(3,j))
         savfrac(1,j) = savfrac(1,j) + dfrac(1,j)
         savfrac(2,j) = savfrac(2,j) + dfrac(2,j)
         savfrac(3,j) = savfrac(3,j) + dfrac(3,j)
         f1 = savfrac(1,j)+half
         f2 = savfrac(2,j)+half
         f3 = savfrac(3,j)+half
         if( ifbox == 1 ) then
            imagcrds(1,i) = f1*ucell(1,1)
            imagcrds(2,i) = f2*ucell(2,2)
            imagcrds(3,i) = f3*ucell(3,3)
         else
            imagcrds(1,i) = f1*ucell(1,1) + f2*ucell(1,2) + f3*ucell(1,3)
            imagcrds(2,i) = f1*ucell(2,1) + f2*ucell(2,2) + f3*ucell(2,3)
            imagcrds(3,i) = f1*ucell(3,1) + f2*ucell(3,2) + f3*ucell(3,3)
         end if
      end if                 ! mygrdlist == 1
   end do                    ! while i<natom
   
   return
end subroutine adjust_imagcrds 

!-------------------------------------------------------------------
!     --- MAP_COORDS ---
!     ...this routine takes the cartesian coordinates of the atoms
!     and maps them to fractional coords. Note that everything is
!     done atom by atom, i.e. we are not doing residue based or
!     molecule based translations.  The fractional coordinates for
!     atoms, are obtained from the dot products with the reciprocal
!     lattice vectors...


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine map_coords here]
subroutine map_coords(crd,natom,recip)
   
   use constants, only : half, zero, one
   implicit none
#  include "ew_cntrl.h"
#include "extra.h"
#include "box.h"
   integer natom
   _REAL_ crd(3,natom)
   _REAL_ recip(3,3)
   integer i
   _REAL_ anint,fracmax,fracmin
   
   fracmax=half
   fracmin=half
   
   if(periodic == 0)then
      do i = 1,natom
         fraction(1,i) = (crd(1,i)+xbox0)*recip(1,1)
         fraction(2,i) = (crd(2,i)+ybox0)*recip(2,2)
         fraction(3,i) = (crd(3,i)+zbox0)*recip(3,3)
         fracmax=max(fracmax,fraction(1,i), &
               fraction(2,i),fraction(3,i))
         fracmin=min(fracmin,fraction(1,i), &
               fraction(2,i),fraction(3,i))
      end do
      if(fracmax > one .or. fracmin < zero)then
         write(6,*)"Frac coord min, max: ",fracmin,fracmax
         write(6,*)"The system has extended beyond "
         write(6,*)"    the extent of the virtual box."
         write(6,*)"Restarting sander will recalculate"
         write(6,*)"   a new virtual box with 30 Angstroms"
         write(6,*)"   extra on each side, if there is a"
         write(6,*)"   restart file for this configuration."
         call sander_bomb('Routine: map_coords (ew_force.f)', &
               'Atom out of bounds. If a restart has been written,', &
               'restarting should resolve the error')
      end if
      return
   end if

   !     ---get fractiontionals
   
   if( ifbox == 1 ) then  !  orthogonal unit cell
   
      do i = 1,natom
         fraction(1,i) = crd(1,i)*recip(1,1)
         fraction(2,i) = crd(2,i)*recip(2,2)
         fraction(3,i) = crd(3,i)*recip(3,3)
      end do

   else
      do i = 1,natom
         fraction(1,i) = crd(1,i)*recip(1,1)+crd(2,i)*recip(2,1)+ &
               crd(3,i)*recip(3,1)
         fraction(2,i) = crd(1,i)*recip(1,2)+crd(2,i)*recip(2,2)+ &
               crd(3,i)*recip(3,2)
         fraction(3,i) = crd(1,i)*recip(1,3)+crd(2,i)*recip(2,3)+ &
               crd(3,i)*recip(3,3)
      end do
   end if
   
   !      --- Check if system has gone out of box for nonperiodic
   !      --- simulations with finite cutoff
   
   if(periodic == 0)then
      if( .not. nocutoff)then
         boxbad=.false.
         do i = 1,natom
            if(anint(fraction(1,i)-half) /= zero)boxbad=.true.
            if(anint(fraction(2,i)-half) /= zero)boxbad=.true.
            if(anint(fraction(3,i)-half) /= zero)boxbad=.true.
         end do
      end if
      if(boxbad) write(6,*) "**********BOX IS BAD****************"
   else
      
      !      --- map them inside, if periodic:
      
      do i = 1,natom
         fraction(1,i) = fraction(1,i) - anint(fraction(1,i))
         fraction(2,i) = fraction(2,i) - anint(fraction(2,i))
         fraction(3,i) = fraction(3,i) - anint(fraction(3,i))
      end do
   end if
   return
end subroutine map_coords 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Store the fractional coordinates.
!
!-----------------------------------------------------------------------
!     --- SAVE_FRAC_CRDS ---
!-----------------------------------------------------------------------
!     This is needed in case you are wrapping coords in the box;
!     note, td actually runs ewald without wrapping since it is easier
!     to analyze but it should work equally well, either way.

subroutine save_frac_crds(natom)
   implicit none
   integer natom
   !     These are actually two-dimensional (3,natom), but to enable
   !     vectorization on IA32 SSE platforms they are treated as
   !     one-dimensional; this may also improve software pipelining !

   integer i
   do i = 1,natom
      savfrac(1,i) = fraction(1,i)
      savfrac(2,i) = fraction(2,i)
      savfrac(3,i) = fraction(3,i)
   end do
   return
end subroutine save_frac_crds 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Store the atomic coordinates.
!
!-----------------------------------------------------------------------
!     --- SAVE_CRDS ---
!-----------------------------------------------------------------------
!   Save the atomic coordinates that were used to build the list.
!   This is needed for the skin test for buffered pairlists,
!   where the current coordinates are compared with the saved
!   coordinates relative to the skin criterion; see check_skin.
!--------------------------------------------------------------------

subroutine save_crds(natom,crd)
   implicit none
   integer natom
   _REAL_  crd(3,natom)
   !     These are actually two-dimensional (3,natom), but to enable
   !     vectorization on IA32 SSE platforms they are treated as
   !     one-dimensional; this may also improve software pipelining !

   integer i
   do i = 1,natom
      savcrd(1,i) = crd(1,i)
      savcrd(2,i) = crd(2,i)
      savcrd(3,i) = crd(3,i)
   end do
   return
end subroutine save_crds 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Determine if the list needs rebuilding.
!-----------------------------------------------------------------------
!     --- CHECK_SKIN ---
!---------------------------------------------------------------------
!  Check if any atom has moved more than half the skin distance,
!  which is half the nbskin added to the cutoff in generating the
!  verlet list; in which case a list build is flagged.
!  An obvious parallel implementation, in which each processor checks its
!  atoms and communicates the results to all other processors, produced
!  on Linux clusters large losses with large numbers of processors.
!  The separate sections based on nbtell exist for improved performance;
!  computing the list update info has a small performance cost.

subroutine check_skin(crd,do_list_update)

   use trace
   use constants, only : zero, fourth
   implicit none
   _REAL_, intent(in) :: crd(3,natom)    ! current atom coordinates, intent(in)
   logical, intent(out) :: do_list_update  ! true if a list update is needed, intent(out)

#  include "extra.h"
#  include "memory.h"

   integer first_atom
   integer last_atom
   integer i
   integer nmoved_atoms     ! total number of atoms triggering a list update
   _REAL_  dx,dy,dz,dis2
   _REAL_  maxdis2

   call trace_enter( 'check_skin' )
   steps_since_list_build = steps_since_list_build + 1
   first_atom = 1
   last_atom  = natom
   maxdis2    = zero

   if ( nbtell == 0 ) then
      ! list update info not requested
      do i = first_atom, last_atom
         dx = crd(1,i) - savcrd(1,i)
         dy = crd(2,i) - savcrd(2,i)
         dz = crd(3,i) - savcrd(3,i)
         dis2 = dx**2 + dy**2 + dz**2
         maxdis2 = max(dis2,maxdis2)
      end do
      do_list_update = maxdis2 > fourth * skinnb * skinnb
   else
      ! list update info requested
      nmoved_atoms = 0
      do i = first_atom, last_atom
         dx = crd(1,i) - savcrd(1,i)
         dy = crd(2,i) - savcrd(2,i)
         dz = crd(3,i) - savcrd(3,i)
         dis2 = dx**2 + dy**2 + dz**2
         if ( dis2 > fourth * skinnb * skinnb) then
            maxdis2 = max(dis2,maxdis2)
            nmoved_atoms = nmoved_atoms + 1
         end if
      end do
      do_list_update = nmoved_atoms > 0
      if ( do_list_update ) then
         if ( master ) then
            write(6, '(1x,A,I7,/,1x,A,F8.4,I7)') &
                  'List Build Triggered: Number of atoms triggering = ', &
                  nmoved_atoms, ' Maximum distance moved = ',sqrt(maxdis2)
         end if
      end if
   end if  ! ( nbtell == 0 )

   call trace_logical('do_list_update = ', do_list_update)
   call trace_exit( 'check_skin' )
   return
end subroutine check_skin

!-------------------------------------------------------------------
!        FILL_TRANVEC    fill the translate vector array
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_tranvec here]
subroutine fill_tranvec()
   implicit none
   integer iv,i0,i1,i2,i3
   
   iv=0
   do i3=0,1
      do i2=-1,1
         do i1= -1,1
            iv=iv+1
            do i0=1,3
               tranvec(i0,iv)= &
                     i1*ucell(i0,1) + &
                     i2*ucell(i0,2) + &
                     i3*ucell(i0,3)
            end do
         end do
      end do
   end do
   return
end subroutine fill_tranvec 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_xtran here]
subroutine fill_xtran()
   implicit none
!!!#include "extra.h"
   integer i,j
! Hard Wired for nghb1=3 thus the hard dimensions (7,10)
   REQUIRE ( nghb1 == 3)
!
!  The neighbor cells of a cell of interest (call it A) that are only
!      "forward" of that cell reside in the plane of the cell,
!      and three planes "forward" in the z direction.
!
!         A = cell for which we want to get neighbor list
!         x = forward neighbor cell within 3 cells
!         o = same as x except this cell has same x index as A
!
!            i3          i3+1            i3+2         i3+3
!        ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!        ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!        ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
! ---->       Axxx...  ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!                      ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!                      ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!                      ..xxxoxxx...  ..xxxoxxx...  ..xxxoxxx...
!
!
!
! A cell and its neighbors span 
!   the x direction over 7 cells (3 on each side and the cell iself)
! XTRAN array contains a 0, 1, or -1 for whether the neighbor cell
!    has x index outside the unit cell and needs to be translated
!    along the x uc vector positive one cell, negative, or not at all.
! There are 10 cases of neighbor cells in x direction
!     xtran(*,1) (0000000) All have x indices within the unit cell
!  cases 2,3,4 are special for the row of neighbors containing the 
!        cell A itself (see arrow). This is the only case where neighbors with
!        x index are not included in the list search since those
!        cells are "before" the cell of interest, and this is
!        a "look forward only" method. So for this row of cells,
!        only 4 xtran values are needed: cell A, and the three cells to right.
!        The cases represent whether the neighbors extend out of the 
!        unit cell by one, two, or three cells. Entry 1 is for cell A and
!        must be 0 since it must be in the unit cell. (last 4 entries are 
!                                                        ignored for this set)
!          (*,2) (0001000)
!          (*,3) (0011000)
!          (*,4) (0111000)
!  cases 5,6,7 are same as 2,3,4 except that there are 7 cells in all other
!        rows
!          (*,5) (0000001)
!          (*,6) (0000011)
!          (*,7) (0000111)
!  cases 8,9,10 are for neighbors that extend to the left out of the UC.
!          (*,8) (-1000000)
!          (*,9) (-1-100000)
!          (*,10)(-1-1-10000)
   !----- most cells will not be translated, so set xtran to 0
   !----- for all possibilities, and then fill in the 1 and -1
   !----- entries for neighbor cells that are beyond uc edges
   !
   !----- zero out entire array ----- CASE 1
   do i=1, 2*nghb1+1
      do j=1,nghb1*3+1
         xtran(i,j)=0
      end do
   end do
   !
   !-----CASES 2,3,4
   do i = 0,nghb1-1
      do j=nghb1+1-i,nghb1+1
         xtran(j,i+2)=1
      end do
   end do
   !
   !------CASES 5,6,7
   do i=1,nghb1
      do j=1,i
         xtran(j,nghb1+1+i)=-1
      end do
   end do
   !
   !------CASES 8,9,10
   do i = 0,nghb1-1
      do j = 0,i
         xtran(2*nghb1+1-j,2*nghb1+2+i)=1
      end do
   end do

   return
end subroutine fill_xtran 



end module nblist
