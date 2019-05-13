! <compile=optimized>
module emap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  EMAP module  for flexible fitting macromolecules
!  into experimental maps
!         By Xiongwu Wu:  wuxw@nhlbi.nih.gov
!            updated 10/1/2011
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  implicit none

#include "assert.fh"
#include "copyright.h"
#include "dprec.fh"

public   

    TYPE EMAPOBJECT           ! map object type
      integer :: mapid              ! map object identity
      character(len=80) mapfile
      integer :: type               ! map object type
      integer :: MODE               ! map object mode
      character(LEN=80) :: name     ! map object name
      character(LEN=516) :: info    ! information about the map object
      _REAL_,pointer :: xatom(:)=> null(),yatom(:)=> null(),zatom(:)=> null()        !  pointer to atom coordinates
      integer :: lx,ly,lz     !  map size in grid numbers (voxels)
      integer :: mx,my,mz     !  starting indexes (voxels)
      _REAL_    :: avg,std,max,min   ! statistics
      _REAL_    :: alpha,beta,gamma   ! lattice angles (degrees)
      _REAL_    :: cx,cy,cz  !  map center in the map coordinate (vexols)
      _REAL_    :: dx,dy,dz  !  grid sizes (angstroms/vexol)
      _REAL_    :: ux,uy,uz           !  unit cell sizes (angstroms)
      _REAL_    :: ox,oy,oz           !  absolute coordinates of map center (angstroms)
      _REAL_    :: sx,sy,sz           !  parameters for reduced coordinates (angstroms)
      integer(kind=1),pointer :: bdata(:)=> null()     !   integer data
      integer(kind=2),pointer :: idata(:)=> null()     !   integer data
      real(kind=4),pointer:: rdata(:)=> null()        !   real data
      complex(kind=8), pointer :: Cdata(:)=> null()  !   double precision data
      type (emapobject), pointer :: next   !  pointer to next emap object
      type (emapobject), pointer :: ref    !  pointer to a reference emap object
    END TYPE EMAPOBJECT

    TYPE EMAPRIGID        ! rigid domain type
      integer :: rigid      ! rigid domain identity
      integer :: mapid      ! rigid domain identity
      _REAL_ fcons
      _REAL_ :: RESOLUTION        !  resolution of the map object
      character(len=256) atmask
      character(LEN=80) :: name     ! rigid domain name
      character(LEN=80) :: mapfit     ! output map file name
      character(LEN=80) :: molfit     ! output molecular file name
      integer :: natc             !  number of atoms represented by the map
      integer, pointer, dimension(:) :: idxatom => null()  !  pointer to atom index array
      logical :: movable,fitting        ! flag indicating a fixed rigid domain
      _REAL_ :: mass,winv,tcm         ! mass, invert
      _REAL_ :: energy         ! energy
      INTEGER,dimension(6) :: grids   ! grid numbers in x,y,z,phi,psi,theta
      integer :: nminim             !  number of minimization steps
      _REAL_,dimension(3) :: prinv   ! inverse of principle inertia
      _REAL_,dimension(3) :: center         ! Position center coordinates 
      _REAL_,dimension(3) :: eular ! Eular angles:phi,theta,psi  
      _REAL_,dimension(4) :: quad   ! quarterions  
      _REAL_,dimension(3) :: v         ! velocities 
      _REAL_,dimension(3) :: omiga        ! angular velocities 
      _REAL_,dimension(3) :: angm        ! angular momentum 
      _REAL_,dimension(3) :: TRAN,TRANREF ! Translation vectors
      _REAL_,dimension(3,3) :: ROT,ROTREF ! rotational matrix
      _REAL_,dimension(3) :: force ! force on center of mass
      _REAL_,dimension(3) :: torq ! torque
      type (emapobject), pointer :: emap ! pointer to the associated emap object
      type (emaprigid), pointer :: next ! pointer to next rigid domain
      type (emaprigid), pointer :: ref  ! pointer to a reference rigid domain
      integer :: nc      ! number of constraint points
      real, pointer :: xc,yc,zc  ! coordinate array of constraint points
    END TYPE EMAPRIGID


    LOGICAL TEMAP
    INTEGER NEMAP,NRIGID,BORDER
    _REAL_ GAMMAMAP,EKMAP,MAPDT,XAVG0,XAVG1,EMSIGMA
    type(EMAPRIGID), save :: rigmin
    type(EMAPOBJECT), allocatable,dimension(:),save :: EMAPS
    type(EMAPRIGID), allocatable,dimension(:),save  :: EMRIGS

contains



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE EMAP_OPTIONS(UNIT)
!_________________________________________________________________
!  read in map constraint information
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      INTEGER UNIT
      character(len=80) mapfile
      character(len=80) mapfit,molfit
      character(len=256) atmask
      _REAL_ fcons,resolution
      INTEGER move,ifit,minim,grids(6)
      namelist /emap/mapfile,mapfit,molfit,atmask,fcons,move,ifit,grids,minim,resolution
      INTEGER I,J,IFIND,mapid,rigid,alloc_err
      nrigid=0
      nemap=0
!  search for namelist EMAP
      rewind UNIT
      DO
        call nmlsrc('emap',UNIT,ifind)
        if (ifind == 0) exit
        read(UNIT,nml=emap)
        nrigid=nrigid+1
      ENDDO
      if(allocated(EMRIGS))deallocate(EMRIGS)
      allocate(EMRIGS(NRIGID),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate EMRIGS"
      rewind UNIT
      DO I=1,NRIGID
        mapfile=''
        atmask=':*'
        mapfit=''
        molfit=''
        move=0
        ifit=0
        minim=100
        grids=1
        fcons=0.05d0
        resolution=0.0d0
        rigid=i
        read(UNIT,nml=emap)
        if(len_trim(mapfile)==0)then
          nemap=nemap+1
          mapid=nemap
        else
          DO J=1,I-1
            if(emrigs(j)%name==mapfile)then
              mapid=emrigs(j)%mapid
              exit
            endif
          ENDDO
          if(j==i)then
            nemap=nemap+1
            mapid=nemap
          endif
        endif
        call rigid_init(emrigs(i))
        emrigs(i)%name=mapfile
        emrigs(i)%mapfit=mapfit
        emrigs(i)%molfit=molfit
        emrigs(i)%rigid=rigid
        emrigs(i)%atmask=atmask
        emrigs(i)%mapid=mapid
        emrigs(i)%fcons=fcons
        emrigs(i)%resolution=resolution
        emrigs(i)%movable=move>0
        emrigs(i)%fitting=ifit>0
        emrigs(i)%grids=grids
        emrigs(i)%nminim=minim
      ENDDO
      if(allocated(EMAPS))deallocate(EMAPS)
      allocate(EMAPS(NEMAP),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate EMAPS"
      J=0
      DO I=1,NRIGID
        if(emrigs(i)%mapid>j)then
          j=emrigs(i)%mapid
          emaps(j)%mapfile=emrigs(i)%name
        endif
      ENDDO
999   return
      end subroutine emap_options


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE PEMAP(dt,temp,x,ix,ih)
!_________________________________________________________________
!  Preprocessing EMAP module
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use constants, only: KB
      use findmask,only:atommask
#ifdef MPI
#  include "parallel.h"
   include 'mpif.h'
#endif
#  include "memory.h"
      _REAL_ dt,temp
      _REAL_ x(*)
      integer ix(*)
      character(len=4) :: ih(*)

      INTEGER I,J,nm,iat,i3,j3,idemp
      INTEGER alloc_err,ierr,UNITOUT
      _REAL_ RESO
      INTEGER, allocatable, dimension(:)::mask
      _REAL_, allocatable, dimension(:)::TMP_CRD,TMP_AMASS
!
      BORDER=6
      UNITOUT=6
      mapdt=20.455d0*20.455d0*dt/GAMMAMAP
      EKMAP=3.0d0*KB*300.0d0*GAMMAMAP*GAMMAMAP/(20.455d0*20.455d0)
      IF(TEMP>0.0d0)EKMAP=3.0d0*KB*TEMP*GAMMAMAP*GAMMAMAP/(20.455d0*20.455d0)
      XAVG1=0.01d0
      XAVG0=1.0d0-XAVG1
#ifdef MPI
      if ( sandersize > 1 ) then
        if ( sanderrank > 0 ) then
          UNITOUT=-1
          if(allocated(EMRIGS))deallocate(EMRIGS)
          allocate(EMRIGS(NRIGID),stat=alloc_err)
          if(alloc_err /= 0 ) write(6,*)"unable to allocate EMRIGS"
          if(allocated(EMAPS))deallocate(EMAPS)
          allocate(EMAPS(NEMAP),stat=alloc_err)
          if(alloc_err /= 0 ) write(6,*)"unable to allocate EMAPS"
        end if
        DO I=1,NRIGID
          if ( sanderrank > 0 )call rigid_init(emrigs(i))
          call mpi_bcast(emrigs(i)%mapid,1,MPI_INTEGER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%rigid,1,MPI_INTEGER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%name,80,MPI_CHARACTER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%mapfit,80,MPI_CHARACTER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%atmask,256,MPI_CHARACTER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%fcons,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%resolution,1,MPI_DOUBLE_PRECISION,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%movable,1,MPI_LOGICAL,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%fitting,1,MPI_LOGICAL,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%nminim,1,MPI_INTEGER,0,commsander,ierr)
          call mpi_bcast(emrigs(i)%grids,6,MPI_INTEGER,0,commsander,ierr)
        ENDDO
        DO I=1,NEMAP
          call mpi_bcast(emaps(i)%mapfile,80,MPI_CHARACTER,0,commsander,ierr)
        ENDDO
      end if
#endif

      DO I=1,NEMAP
        CALL RDEMAP(EMAPS(I)%MAPFILE,I,UNITOUT)
      ENDDO
      if(allocated(mask))deallocate(mask)
      allocate(mask(NATOM),stat=alloc_err)
      DO I=1,NRIGID
        call atommask( natom, nres, 0, ih(m04), ih(m06), &
                ix(i02), ih(m02), x(lcrd), emrigs(i)%atmask, mask )
        nm=0
        DO J=1,NATOM
          if(mask(J)>0)nm=nm+1
        ENDDO
        emrigs(i)%natc=nm
        if(associated(EMRIGS(I)%IDXATOM))deallocate(EMRIGS(I)%IDXATOM)
        allocate(EMRIGS(I)%IDXATOM(NM),stat=alloc_err)
        if(alloc_err /= 0 ) write(6,*)"unable to allocate IDXATOM"
        nm=0
        DO J=1,NATOM
          if(mask(J)>0)then
            nm=nm+1
            emrigs(i)%idxatom(nm)=j
          endif
        ENDDO
        call rigid_prop(x(lmass),x(lcrd),emrigs(i))
        if(len_trim(emrigs(i)%name)==0)then
          ! create map from the coordinates of masked atoms
          if(allocated(TMP_CRD))deallocate(TMP_CRD)
          allocate(TMP_CRD(3*NM),stat=alloc_err)
          if(alloc_err /= 0 ) write(6,*)"unable to allocate TMP_CRD"
          if(allocated(TMP_AMASS))deallocate(TMP_AMASS)
          allocate(TMP_AMASS(NM),stat=alloc_err)
          if(alloc_err /= 0 ) write(6,*)"unable to allocate TMP_AMASS"
          DO IAT=1,NM
            J=emrigs(i)%idxatom(iat)
            tmp_amass(iat)=x(lmass+j-1)
            i3=3*iat-2
            j3=lcrd+3*j-3
            tmp_crd(i3)=x(j3)
            tmp_crd(i3+1)=x(j3+1)
            tmp_crd(i3+2)=x(j3+2)
          enddo
          IDEMP=emrigs(i)%mapid
          RESO=emrigs(i)%resolution
          if(RESO<0.1D0)RESO=2.0D0
          if(unitout>0)then
            write(unitout,'("Map ",I4," is created from ",i6," constrained atoms with resolution: ",F4.1)') &
              IDEMP,nm, RESO
          endif
          CALL COR2MAP(IDEMP,RESO,nm,tmp_amass,tmp_crd,unitout)      
        endif
        if(unitout>0)then
           write(unitout,'("Rigid ",I4," has ",i6," constrained atoms with mask: ",A)') &
              I,nm, emrigs(i)%atmask
           write(unitout,'("Fitting: ",L1," movability ",L1)')emrigs(i)%fitting, emrigs(i)%movable
        endif
        ! perform map to atom fit if required
        if(emrigs(i)%fitting)then
           j=emrigs(i)%mapid
           call emapfit(emrigs(i),emaps(j),UNITOUT)
           if(unitout>0)then
              write(unitout,'("Rigid ",I4," has been fitted to: ",A," with energy:",E14.6)') &
              I, emaps(j)%mapfile,emrigs(i)%energy
           endif
        endif
      ENDDO
      if(allocated(MASK))deallocate(mask)
      if(allocated(TMP_AMASS))deallocate(tmp_amass)
      if(allocated(TMP_CRD))deallocate(tmp_crd)
      return
      end subroutine pemap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE QEMAP()
!_________________________________________________________________
!  Clean and quit EMAP module
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
#ifdef MPI
#  include "parallel.h"
   include 'mpif.h'
#endif
#  include "memory.h"

      INTEGER I,J,nm,dealloc_err,ierr,UNITOUT
!
      UNITOUT=6
#ifdef MPI
      if ( sanderrank > 0 ) then
        UNITOUT=-1
      end if
#endif

      DO I=1,NRIGID
        if(UNITOUT>0)call RIGMAP(emrigs(i),UNITOUT)
        if(associated(EMRIGS(I)%IDXATOM))deallocate(EMRIGS(I)%IDXATOM,stat=dealloc_err)
        if(dealloc_err /= 0 ) write(6,*)"unable to deallocate IDXATOM for RIGID ",I
      ENDDO
      if(allocated(EMRIGS))deallocate(EMRIGS)
      DO I=1,0
        if(associated(EMAPS(I)%BDATA))deallocate(EMAPS(I)%BDATA,stat=dealloc_err)
        if(dealloc_err /= 0 ) write(6,*)"unable to deallocate BDATA for EMAP ",I
        if(associated(EMAPS(I)%CDATA))deallocate(EMAPS(I)%CDATA,stat=dealloc_err)
        if(dealloc_err /= 0 ) write(6,*)"unable to deallocate CDATA for EMAP ",I
        if(associated(EMAPS(I)%IDATA))deallocate(EMAPS(I)%IDATA,stat=dealloc_err)
        if(dealloc_err /= 0 ) write(6,*)"unable to deallocate IDATA for EMAP ",I
        if(associated(EMAPS(I)%RDATA))deallocate(EMAPS(I)%RDATA,stat=dealloc_err)
        if(dealloc_err /= 0 ) write(6,*)"unable to deallocate RDATA for EMAP ",I
      ENDDO
      if(allocated(EMAPS))deallocate(EMAPS)
      return
      end subroutine qemap

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE EMAPFIT(RIGOBJ,MAPOBJ,OUTU)
!_________________________________________________________________
!  fit map to a rigid with grid-threading minimization
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
   use memory_module, only: xx =>x
   use constants, only: PI,TWOPI
      type(EMAPRIGID) :: rigobj
      type(EMAPOBJECT) :: mapobj
      INTEGER OUTU

#  include "memory.h"
#  include "md.h"
      INTEGER I,J,NDATA,alloc_err
      INTEGER IDEMP,UNIT,ktx,kty,ktz,krx,kry,krz
      INTEGER it,jt,kt,ir,jr,kr,m,iter
      _REAL_ dtx,gx,gy,gz,phi,psi,tht,dphi,dpsi,dtht
      _REAL_ xl,yl,zl,enrig,pos0(6),pos(6),posi(6),gtol
      _REAL_ p(7,6),engs(7),engi
!
      dtx = dt*20.455d+00
      gtol=1.0d-20
!
      ktx=rigobj%grids(1)
      kty=rigobj%grids(2)
      ktz=rigobj%grids(3)
      krx=rigobj%grids(4)
      kry=rigobj%grids(5)
      krz=rigobj%grids(6)
      xl=mapobj%dx*mapobj%mx
      yl=mapobj%dy*mapobj%my
      zl=mapobj%dz*mapobj%mz
      gx=mapobj%dx*mapobj%lx/(ktx)
      gy=mapobj%dy*mapobj%ly/(kty)
      gz=mapobj%dz*mapobj%lz/(ktz)
      dphi=TWOPI/krx
      dpsi=TWOPI/kry
      dtht=PI/krz
      rigmin=rigobj
      pos0(1)=rigobj%tran(1)+anint((rigobj%center(1)-rigobj%tran(1)-xl)/gx-0.5d0)*gx
      pos0(2)=rigobj%tran(2)+anint((rigobj%center(2)-rigobj%tran(2)-yl)/gy-0.5d0)*gy
      pos0(3)=rigobj%tran(3)+anint((rigobj%center(3)-rigobj%tran(3)-zl)/gz-0.5d0)*gz
      pos0(4:6)=rigobj%eular
      do it=1,ktx
        pos(1)=pos0(1)-(it-1)*gx
        do jt=1,kty
          pos(2)=pos0(2)-(jt-1)*gy
          do kt=1,ktz
            pos(3)=pos0(3)-(kt-1)*gz
            do ir=1,krx
              pos(4)=pos0(4)+(ir-1)*dphi
              do jr=1,kry
                pos(5)=pos0(5)+(jr-1)*dpsi
                do kr=1,krz
                  pos(6)=pos0(6)+(kr-1)*dtht
                  do i=1,7
                    p(i,1:6)=pos
                  enddo                  
                  p(1,1)=p(1,1)+gx/2
                  p(2,2)=p(2,2)+gy/2
                  p(3,3)=p(3,3)+gz/2
                  p(4,4)=p(4,4)+dphi/2
                  p(5,5)=p(5,5)+dpsi/2
                  p(6,6)=p(6,6)+dtht/2
                  do i=1,7
                    posi=p(i,1:6)
                    engs(i)=rigeng(posi)
                  enddo
                  if(outu>0)write(outu,'("initial:",6f8.2,e14.6)')(p(7,i),i=1,6),engs(7)
                  call simplex(p,engs,7,6,6,gtol,rigeng,iter)  
                  if(outu>0)write(outu,'("result:",6f8.2,e14.6)')(p(1,i),i=1,6),engs(1)
                  if(rigmin%energy<rigobj%energy)then
                     rigobj=rigmin
                 endif
                enddo  !kr
              enddo  !jr
            enddo  !ir
          enddo  !kt
        enddo  !jt
      enddo  !it
      RETURN
      END SUBROUTINE EMAPFIT


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE rigid_init(rigobj)
!_________________________________________________________________
!  initialize rigid domain
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      type(EMAPRIGID) :: rigobj

      INTEGER I,J

      rigobj%energy=1.0D9
      DO I=1,3
        DO J=1,3
          rigobj%ROT(I,J)=0.0D0
        ENDDO
        rigobj%ROT(I,I)=1.0D0
        rigobj%TRAN(I)=0.0D0
        rigobj%center(I)=0.0D0
        rigobj%FORCE(I)=0.0D0
        rigobj%TORQ(I)=0.0D0
        rigobj%QUAD(I+1)=0.0D0
      ENDDO
      rigobj%QUAD(1)=1.0D0
      return
      end subroutine rigid_init

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE rigid_prop(amass,crd,rigobj)
!_________________________________________________________________
!  initialize rigid domain
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      _REAL_ amass(*),CRD(*)
      type(EMAPRIGID) :: rigobj

      INTEGER I,J,N,IA
      _REAL_ AMASSI,X,Y,Z,XC,YC,ZC,TCM,TMASS

      N=rigobj%natc

      XC=0.0D0
      YC=0.0D0
      ZC=0.0D0
      DO I=1,N
        IA=rigobj%idxatom(i)
        AMASSI=AMASS(IA)
        X=crd(3*I-2)
        Y=crd(3*I-1)
        Z=crd(3*I)
        XC=XC+X
        YC=YC+Y
        ZC=ZC+Z
      ENDDO
      XC=XC/N
      YC=YC/N
      ZC=ZC/N
      TCM=0.0D0
      TMASS=0.0D0
      DO I=1,N
        IA=rigobj%idxatom(i)
        AMASSI=AMASS(IA)
        X=crd(3*I-2)-XC
        Y=crd(3*I-1)-YC
        Z=crd(3*I)-ZC
        TCM=TCM+AMASSI*(X*X+Y*Y+Z*Z)
        TMASS=TMASS+AMASSI
      ENDDO
      rigobj%mass=TMASS
      rigobj%tcm=TCM
      rigobj%center(1)=xc
      rigobj%center(2)=yc
      rigobj%center(3)=zc
      return
      end subroutine rigid_prop

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE EMAP_move()
!_________________________________________________________________
!  Perform motion step for all rigid domains
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      _REAL_ forcetmp(3),torqtmp(3),enemap

      INTEGER I,J,OUTU,IERR
      enemap=0.0d0
      DO I=1,NRIGID
        enemap=enemap+emrigs(i)%energy
        if(emrigs(i)%movable)call rigid_move(emrigs(i))
        DO J=1,3
          emrigs(i)%force(j)=0.0d0
          emrigs(i)%torq(j)=0.0d0
        ENDDO
        emrigs(i)%energy=0.0d0
      ENDDO
      return
      end subroutine EMAP_move

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE rigid_move(rigobj)
!_________________________________________________________________
!  moving rigid domain
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!
!     dEmap/dr=-f    r=(x,y,z)
!     dEmap/da=-Ta=(dEmap/drc)(drc/da)=-T(d(A*r0)/da)=-T(dA/da)(A-)rc
!
!     dr=f*dt2/(2m)
!     da=Ta*dt2/(2I)
!
!_______________________________________________________________________
!
      use constants,only: TWOPI
      
      type(EMAPRIGID) :: rigobj

      INTEGER I
      _REAL_  ff,tt,scalf,scalt
      _REAL_  DRF,DGT,phi,psi,tht,dphi,dpsi,dtht
      _REAL_ a(3,3),daf(3,3),day(3,3),das(3,3),b(3),t(3),c(3)

      ff=dot_product(rigobj%force,rigobj%force)
      tt=dot_product(rigobj%torq,rigobj%torq)
      scalf=sqrt(EKMAP*rigobj%mass/ff)
      scalt=sqrt(EKMAP*rigobj%tcm/tt)
      DRF=min(1.0d0,scalf)*mapdt/rigobj%mass
      DGT=min(1.0d0,scalt)*mapdt/rigobj%tcm
      phi=rigobj%eular(1)
      psi=rigobj%eular(2)
      tht=rigobj%eular(3)

      call rotdeular(a,phi,psi,tht,daf,day,das)
      b(1)=a(1,1)+a(2,1)+a(3,1)
      b(2)=a(1,2)+a(2,2)+a(3,2)
      b(3)=a(1,3)+a(2,3)+a(3,3)
      t=rigobj%torq
      c=matmul(daf,b)
      dphi=dot_product(t,c)
      c=matmul(day,b)
      dpsi=dot_product(t,c)
      c=matmul(das,b)
      dtht=dot_product(t,c)
      rigobj%eular(1)=phi+dphi*dgt
      rigobj%eular(2)=psi+dpsi*dgt
      rigobj%eular(3)=tht+dtht*dgt
      rigobj%eular(1)=rigobj%eular(1)-TWOPI*ANINT(rigobj%eular(1)/TWOPI)
      rigobj%eular(2)=rigobj%eular(2)-TWOPI*ANINT(rigobj%eular(2)/TWOPI)
      rigobj%eular(3)=rigobj%eular(3)-TWOPI*ANINT(rigobj%eular(3)/TWOPI)
      ! Update rotation matrix
      call roteular(rigobj%rot,rigobj%eular(1),rigobj%eular(2),rigobj%eular(3))
      DO I=1,3
        rigobj%TRAN(I)=rigobj%TRAN(I)+rigobj%FORCE(I)*DRF
      ENDDO
      return
      end subroutine rigid_move

      SUBROUTINE roteular(U,PHI,PSI,THT)
!_______________________________________________________________________
!     This routine convert Eular angles
!     to rotational matrix
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!     emap energy and derivatives to x,y,z, and phi,psi,theta
!_______________________________________________________________________
    !     axis: ax(3)  
    !     angle: phir in rad
    !     matrix: u
    !     lOK: validate
    !
    !       { cosFcosY-sinFsinYcosS   sinFcosY+cosFsinYcosS   sinSsinY  }
    !     u={ -cosFsinY-sinFcosYcosS  -sinFsinY+cosFcosYcosS   sinScosY  }
    !       { sinFsinS                 -cosFsinS                cosS     }
    !
    !
    !
      implicit none
    !
      _REAL_ U(3,3),PHI,PSI,THT
    !
      _REAL_ sinF,cosF,sinY,cosY,sinS,cosS
    !
      sinF=sin(phi)
      cosF=cos(phi)
      sinY=sin(psi)
      cosY=cos(psi)
      sinS=sin(tht)
      cosS=cos(tht)
      u(1,1)=cosF*cosY-sinF*sinY*cosS 
      u(1,2)=sinF*cosY+cosF*sinY*cosS
      u(1,3)=sinS*sinY
      u(2,1)=-cosF*sinY-sinF*cosY*cosS
      u(2,2)=-sinF*sinY+cosF*cosY*cosS
      u(2,3)=sinS*cosY 
      u(3,1)=sinF*sinS  
      u(3,2)=-cosF*sinS
      u(3,3)=cosS
      RETURN
      END SUBROUTINE ROTEULAR

      SUBROUTINE rotdeular(A,PHI,PSI,THT,daf,day,das)
!_______________________________________________________________________
!     This routine convert Eular angles
!     to rotational matrix and it derivative to Eular angles
!
!       ( cosf cosy-sinf siny coss,  sinf cosy+cosf siny coss,   siny sins  )
!     A=(-cosf siny-sinf cosy coss, -sinf siny+cosf cosy coss,   cosy sins  )
!       (   sinf sins             ,    -cosf sins            ,    coss      )
!
!       (-sinf cosy-cosf siny coss,  cosf cosy-sinf siny coss,   0  ) (-a12,a11,0)
! dA/df=( sinf siny-cosf cosy coss, -cosf siny-sinf cosy coss,   0  )=(-a22,a21,0)
!       (   cosf sins             ,    sinf sins            ,    0  ) (-a32,a31,0)
!
!       (-cosf siny-sinf cosy coss,-sinf siny+cosf cosy coss, cosy sins) ( a21, a22, a23)
! dA/dy=(-cosf cosy+sinf siny coss,-sinf cosy-cosf siny coss,-siny sins)=(-a11,-a12,-a31)
!       (   0                     ,    0                    ,   0      )=(  0 ,  0 ,  0 )
!
!       ( sinf siny sins,  -cosf siny sins,   siny coss  )    (a13a31,a32a13,a13a33   )
! dA/ds=( sinf cosy sins, -cosf cosy sins,    cosy coss  )=   (a31a23,a32a23,a23a33   )/b33
!       (   sinf coss   ,    -cosf coss   ,   -sins      )    (a31a33,a32a33,a33*a33-1)
!    here b33=sqrt(1-a33*a33)
!_________________________________________________________________
    !
    !
      implicit none
    !
      _REAL_ PHI,PSI,THT
    !
      _REAL_ a(3,3),daf(3,3),day(3,3),das(3,3)
      _REAL_ sinf,cosf,siny,cosy,sins,coss
    !
      sinf=sin(phi)
      siny=sin(psi)
      sins=sin(tht)
      cosf=cos(phi)
      cosy=cos(psi)
      coss=cos(tht)
      a(1,1)=cosF*cosY-sinF*sinY*cosS 
      a(1,2)=sinF*cosY+cosF*sinY*cosS
      a(1,3)=sinS*sinY
      a(2,1)=-cosF*sinY-sinF*cosY*cosS
      a(2,2)=-sinF*sinY+cosF*cosY*cosS
      a(2,3)=sinS*cosY 
      a(3,1)=sinF*sinS  
      a(3,2)=-cosF*sinS
      a(3,3)=cosS
      daf=0.0d0
      daf(1,1)=-a(1,2)
      daf(1,2)=a(1,1)
      daf(2,1)=-a(2,2)
      daf(2,2)=a(2,1)
      daf(3,1)=-a(3,2)
      daf(3,2)=a(3,1)
      day=0.0d0
      day(1,1)=a(2,1)
      day(1,2)=a(2,2)
      day(1,3)=a(2,3)
      day(2,1)=-a(1,1)
      day(2,2)=-a(1,2)
      day(2,3)=-a(1,3)
      das(1,1)=sinf*siny*sins
      das(1,2)=-cosf*siny*sins
      das(1,3)=siny*coss 
      das(2,1)=sinf*cosy*sins 
      das(2,2)=-cosf*cosy*sins 
      das(2,3)=cosy*coss 
      das(3,1)=sinf*coss 
      das(3,2)=-cosf*coss 
      das(3,3)= -sins 
      RETURN
      END SUBROUTINE ROTDEULAR


      SUBROUTINE rotaxisphi(U,RN,PHIR)
    !-------------  --------------------------------------------------------
    !     This routine convert rotation axis and angle
    !     to rotational matrix
    !     axis: ax(3)  
    !     angle: phir in rad
    !     matrix: u
    !     lOK: validate
    !
    !       { cos+ux2(1-cos)   uxuy(1-cos)-uzsin   uxuz(1-cos)+uysin   }
    !     u={uyux(1-cos)+uzsin   cos+uy2(1-cos)   uyuz(1-cos)-uxsin    }
    !       {uzux(1-cos)-uysin   uzuy(1-cos)+uxsin   cos+uz2(1-cos)    }
    !
    !     Xiongwu Wu
    !
    !
      implicit none
    !
      _REAL_ U(3,3),RN(3),PHIR
    !
      _REAL_ rnorm,rx,ry,rz,cosp,cosp1,sinp
      INTEGER I,J
    !
      cosp=cos(phir)
      cosp1=1.0D0-cosp
      sinp=sin(phir)
      rnorm=sqrt(rn(1)*rn(1)+rn(2)*rn(2)+rn(3)*rn(3))
      if(rnorm<1.0d-8)then
        u(1,1)=1.0d0
        u(1,2)=0.0d0
        u(1,3)=0.0d0
        u(2,1)=0.0d0
        u(2,2)=1.0d0
        u(2,3)=0.0d0
        u(3,1)=0.0d0
        u(3,2)=0.0d0
        u(3,3)=1.0d0
        return
      endif
      rx=rn(1)/rnorm
      ry=rn(2)/rnorm
      rz=rn(3)/rnorm
      u(1,1)=cosp+rx*rx*cosp1
      u(1,2)=rx*ry*cosp1-rz*sinp
      u(1,3)=rx*rz*cosp1+ry*sinp
      u(2,1)=rx*ry*cosp1+rz*sinp
      u(2,2)=cosp+ry*ry*cosp1
      u(2,3)=ry*rz*cosp1-rz*sinp
      u(3,1)=rx*rz*cosp1-ry*sinp
      u(3,2)=ry*rz*cosp1+rx*sinp
      u(3,3)=cosp+rz*rz*cosp1
      RETURN
      END SUBROUTINE rotaxisphi



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE RIGMAP(RIGOBJ,OUTU)
!_________________________________________________________________
!  write out map and molecule of a rigid 
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      type(EMAPRIGID) :: rigobj
      type(EMAPOBJECT) :: mapobj
      INTEGER OUTU

      INTEGER I,J,NDATA,alloc_err
      INTEGER IDEMP,UNIT
      CHARACTER*80 FNAME
!
      IF(rigobj%mapfit/='')then
        FNAME=rigobj%mapfit
        IDEMP=rigobj%mapid
        NDATA=EMAPS(IDEMP)%LX*EMAPS(IDEMP)%LY*EMAPS(IDEMP)%LZ
        if(associated(MAPOBJ%RDATA))deallocate(MAPOBJ%RDATA)
        allocate(MAPOBJ%RDATA(NDATA),stat=alloc_err)
        if(alloc_err /= 0 ) write(6,*)"unable to allocate RDATA"
        !  Transform map to rigid position
        CALL MAPCAST(MAPOBJ,EMAPS(IDEMP),RIGOBJ%TRAN,RIGOBJ%ROT,0.0D0)
        UNIT=99
        OPEN(UNIT,FILE=FNAME,ACCESS='STREAM',STATUS='UNKNOWN')
        IF(INDEX(FNAME,'ccp4')>0.or.INDEX(FNAME,'CCP4')>0 .or. &
           INDEX(FNAME,'map')>0.or.INDEX(FNAME,'MAP')>0)THEN
          CALL WRTCCP4(UNIT,MAPOBJ,OUTU)      
        ELSE
          write(outu,*)' Map format is not supported. Skip output:',fname
        ENDIF
        CLOSE(UNIT)
        deallocate(MAPOBJ%RDATA)
      ENDIF
      IF(rigobj%molfit/='')then
        FNAME=rigobj%molfit
        IDEMP=rigobj%mapid
        UNIT=99
        OPEN(UNIT,FILE=FNAME,STATUS='UNKNOWN')
        IF(INDEX(FNAME,'pdb')>0.or.INDEX(FNAME,'PDB')>0 )THEN
          CALL WRTFITPDB(UNIT,RIGOBJ%NATC,RIGOBJ%IDXATOM, &
       EMAPS(IDEMP)%CX,EMAPS(IDEMP)%CY,EMAPS(IDEMP)%CZ,RIGOBJ%TRAN,RIGOBJ%ROT,OUTU)      
        ELSE
          write(outu,*)' Map format is not supported. Skip output:',fname
        ENDIF
        CLOSE(UNIT)
      ENDIF
      RETURN
      END SUBROUTINE RIGMAP

      SUBROUTINE MAPCAST(MAPN,MAP,T,U,RESO)
!_________________________________________________________________
!  cast map according transform
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use ew_bspline,only:kbot,ktop,fill_bspline_0
!
      type(EMAPOBJECT) :: mapn,map
      INTEGER UNIT,IDEMP,OUTU
      _REAL_  T(3),U(3,3),RESO
!
      INTEGER*4 LX,LY,LZ,MODE,MNX,MNY,MNZ,NX,NY,NZ
      _REAL_ CX,CY,CZ,DX,DY,DZ,ALPHA,BETA,GAMMA
      INTEGER*4 MAPC,MAPR,MAPS
      REAL*4 AMAX,AMIN,AMEAN,ARMS
      INTEGER*4 ISPG,NSYMBT,LSKFLG,NNOTE
      REAL*4 SKWMAT(9),SKWTRN(3),EXTRA(15)
      CHARACTER*4 MAPLABLE,MACHST
      CHARACTER*80 NOTES(10)
      INTEGER*1 BDATAI
      INTEGER*2 IDATAI
      REAL*4 RDATAI
      COMPLEX(KIND=4) CDATAI
      INTEGER*4 LXN,LYN,LZN,MXN,MYN,MZN,NDATA
      _REAL_ CXN,CYN,CZN,DXN,DYN,DZN,XL,YL,Zl,XH,YH,ZH
!
      INTEGER I,J,K,IN,JN,KN,M1,M2,M3,alloc_err
      INTEGER IPT1,IPT2,IPT3,ITH1,ITH2,ITH3
      _REAL_ X,Y,Z,XN,YN,ZN
      _REAL_ FR1,FR2,FR3,W,RHOI
      REAL*8 BS1(8),BS2(8),BS3(8)
      REAL*8 VAL0,VAL0A
!
      LX=MAP%LX
      LY=MAP%LY
      LZ=MAP%LZ
      MODE=MAP%MODE
      MNX=MAP%MX
      MNY=MAP%MY
      MNZ=MAP%MZ
      DX=MAP%DX
      DY=MAP%DY
      DZ=MAP%DZ
      CX=MAP%CX
      CY=MAP%CY
      CZ=MAP%CZ
      XL=MNX*DX
      YL=MNY*DY
      ZL=MNZ*DZ
      XH=XL+LX*DX
      YH=YL+LY*DY
      ZH=ZL+LZ*DZ
      IF(RESO>0.0D0)THEN
        DXN=2.0D0*RESO/BORDER
        DYN=2.0D0*RESO/BORDER
        DZN=2.0D0*RESO/BORDER
      ELSE
        DXN=DX
        DYN=DY
        DZN=DZ
      ENDIF
      MXN=NINT((XL+T(1))/DXN)
      MYN=NINT((YL+T(2))/DYN)
      MZN=NINT((ZL+T(3))/DZN)
      LXN=NINT((XH-XL)/DXN)
      LYN=NINT((YH-YL)/DYN)
      LZN=NINT((ZH-ZL)/DZN)
      CXN=CX+T(1)
      CYN=CY+T(2)
      CZN=CZ+T(3)
      MAPN=MAP
      NDATA=LXN*LYN*LZN
      allocate(MAPN%RDATA(NDATA),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate RDATA"
      MAPN%LX=LXN
      MAPN%LY=LYN
      MAPN%LZ=LZN
      MAPN%MX=MXN
      MAPN%MY=MYN
      MAPN%MZ=MZN
      MAPN%CX=CXN
      MAPN%CY=CYN
      MAPN%CZ=CZN
      MAPN%DX=DXN
      MAPN%DY=DYN
      MAPN%DZ=DZN
      MAPN%MODE=2
      ALPHA=MAP%ALPHA
      BETA=MAP%BETA
      GAMMA=MAP%GAMMA
      MAPC=1
      MAPR=2
      MAPS=3
      AMIN=MAP%MAX
      AMAX=MAP%MIN
      AMEAN=0.0D0
      ARMS=0.0D0
      ISPG=0
      NSYMBT=0
      LSKFLG=0
      MAPLABLE='MAP '
      MACHST='ALL '
      NNOTE=3
      NOTES(1)=" This map is created with the emap module "
      NOTES(2)=" Report questions to Dr. Xiongwu Wu  "
      NOTES(3)="             Email: wuxw@nhlbi.nih.gov "
      DO KN=1,LZN
        ZN=(MZN+KN-1)*DZN-CZN
        DO JN=1,LYN
          YN=(MYN+JN-1)*DYN-CYN
          DO IN=1,LXN
            XN=(MXN+IN-1)*DXN-CXN
            X=U(1,1)*XN+U(2,1)*YN+U(3,1)*ZN+CX
            Y=U(1,2)*XN+U(2,2)*YN+U(3,2)*ZN+CY
            Z=U(1,3)*XN+U(2,3)*YN+U(3,3)*ZN+CZ
! calculate B-spline parameters
            FR1=(X)/DX
            FR2=(Y)/DY
            FR3=(Z)/DZ
            M1=ANINT(FR1-0.5D0)
            W=FR1-M1
            call fill_bspline_0(w,BORDER,BS1)
            M2=ANINT(FR2-0.5D0)
            W=FR2-M2
            call fill_bspline_0(w,BORDER,BS2)
            M3=ANINT(FR3-0.5D0)
            W=FR3-M3
            call fill_bspline_0(w,BORDER,BS3)
! Calculate B-spline interception
            RHOI=0.0D0
            K = M3-MNZ+1 - BORDER/2
            DO 100 ITH3 = 1,BORDER
               K=K+1
               IF(K<1)THEN
                 IPT1=0
               ELSE IF(K>LZ)THEN
                 IPT1=LX*LY*(LZ-1)
               ELSE
                 IPT1=(K-1)*LY*LX
               ENDIF
               VAL0A =  BS3(ITH3)
!
               J = M2-MNY+1 - BORDER/2
!
               DO 200 ITH2 = 1,BORDER
                 J=J+1
                 IF(J<1)THEN
                   IPT2=IPT1
                 ELSE IF(J>LY)THEN
                   IPT2=IPT1+LX*(LY-1)
                 ELSE
                   IPT2=IPT1+(J-1)*LX
                 ENDIF
!
                 VAL0= VAL0A * BS2(ITH2)
!
                 I = M1-MNX+1 - BORDER/2
!
                 DO 300 ITH1 = 1,BORDER
                   I=I+1
                   IF(I<1)THEN
                     IPT3=IPT2+1
                   ELSE IF(I>LX)THEN
                     IPT3=IPT2+LX
                   ELSE
                     IPT3=IPT2+I
                   ENDIF
                   RHOI=RHOI+ VAL0 * BS1(ITH1)*MAP%RDATA(IPT3)
!
300              CONTINUE
200            CONTINUE
100         CONTINUE
!
            MAPN%RDATA(IN+LXN*(JN-1+LYN*(KN-1)))=RHOI
            IF(RHOI<AMIN)AMIN=RHOI
            IF(RHOI>AMAX)AMAX=RHOI
            AMEAN=AMEAN+RHOI
            ARMS=ARMS+RHOI*RHOI
          ENDDO
        ENDDO
      ENDDO
      AMEAN=AMEAN/NDATA
      ARMS=ARMS/NDATA
      ARMS=ARMS-AMEAN*AMEAN
      IF(ARMS>0.0D0)ARMS=SQRT(ARMS)
      MAPN%AVG=AMEAN
      MAPN%STD=ARMS
      MAPN%MIN=AMIN
      MAPN%MAX=AMAX
      RETURN
      END SUBROUTINE MAPCAST


      SUBROUTINE WRTFITPDB(UNIT,NATRIG,IDXATM,CX,CY,CZ,T,U,OUTU)
!_________________________________________________________________
!  transform atom coordinates to their fitting position
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use memory_module, only: memory_init,  &
            residue_pointer,residue_label,atom_name,coordinate
!
#  include "memory.h"
      INTEGER UNIT,NATRIG,OUTU,IDXATM(*)
      _REAL_  CX,CY,CZ,T(3),U(3,3)
      _REAL_ CRDN(3)
!
      INTEGER IRES,IATOM,I,IA
      INTEGER IPT1,IPT2,IPT3,ITH1,ITH2,ITH3
      _REAL_ XN,YN,ZN
      LOGICAL ATSKIP
      character(len=4) :: name
      character(len=3) :: resName
!
      call memory_init()
      write(unit,'(A)') 'REMARK  Written by Amber 12, SANDER, EMAP'
      do ires = 1,nres
         do iatom = residue_pointer(ires), residue_pointer(ires+1)-1
           ATSKIP=.TRUE.
           DO I=1,NATRIG
             IA=IDXATM(I)
             IF(IATOM==IA)THEN
               ATSKIP=.FALSE.
               EXIT
             ENDIF
           ENDDO
           IF(ATSKIP)CYCLE
           XN=coordinate(1,ia)-T(1)-CX
           YN=coordinate(2,ia)-T(2)-CY
           ZN=coordinate(3,ia)-T(3)-CZ
           CRDN(1)=U(1,1)*XN+U(2,1)*YN+U(3,1)*ZN+CX
           CRDN(2)=U(1,2)*XN+U(2,2)*YN+U(3,2)*ZN+CY
           CRDN(3)=U(1,3)*XN+U(2,3)*YN+U(3,3)*ZN+CZ
            name = atom_name(iatom)
            resName=residue_label(ires)
            resName=adjustr(resName)
            write(unit,'(A6,I5,1X,A4,A1,A3,1X,A1,I4,A1,3X,3F8.3)')&
                  'ATOM  ', &
                  iatom,name,' ', &
                  resName,' ', &
                  ires,' ', &
                  CRDN(1:3)
         end do
      end do
      write(unit,'(A)') 'END'
      RETURN
      END SUBROUTINE WRTFITPDB


      SUBROUTINE WRTCCP4(UNIT,MAP,OUTU)
!_________________________________________________________________
!  write out CCP4 map information and its density array
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      INTEGER UNIT,OUTU
      type(EMAPOBJECT) :: map
!
      INTEGER*4 LX,LY,LZ,MODE,MNX,MNY,MNZ,NX,NY,NZ
      REAL*4 DX,DY,DZ,XL,YL,ZL,ALPHA,BETA,GAMMA
      INTEGER*4 MAPC,MAPR,MAPS
      REAL*4 AMAX,AMIN,AMEAN,ARMS
      INTEGER*4 ISPG,NSYMBT,LSKFLG,NNOTE
      REAL*4 SKWMAT(9),SKWTRN(3),EXTRA(15)
      CHARACTER*4 MAPLABLE,MACHST
      CHARACTER*80 NOTES(10)
      INTEGER*1 BDATAI
      INTEGER*2 IDATAI
      REAL*4 RDATAI
      COMPLEX(KIND=4) CDATAI
!
      INTEGER I,NDATA
!
      LX=MAP%LX
      LY=MAP%LY
      LZ=MAP%LZ
      NDATA=LX*LY*LZ
      MODE=MAP%MODE
      MNX=MAP%MX
      MNY=MAP%MY
      MNZ=MAP%MZ
      NX=MAP%LX
      NY=MAP%LY
      NZ=MAP%LZ
      DX=MAP%DX
      DY=MAP%DY
      DZ=MAP%DZ
      XL=NX*DX
      YL=NY*DY
      ZL=NZ*DZ
      ALPHA=MAP%ALPHA
      BETA=MAP%BETA
      GAMMA=MAP%GAMMA
      MAPC=1
      MAPR=2
      MAPS=3
      AMIN=MAP%MIN
      AMAX=MAP%MAX
      AMEAN=MAP%AVG
      ARMS=MAP%STD
      ISPG=0
      NSYMBT=0
      LSKFLG=0
      SKWMAT=0.0d0
      SKWMAT(1)=1.0D0
      SKWMAT(5)=1.0D0
      SKWMAT(9)=1.0D0
      SKWTRN=0.0d0
      EXTRA=0.0D0
      MAPLABLE='EMAP '
      MACHST='ALL '
      NNOTE=3
      NOTES=""
      NOTES(1)=" This map is created with the emap module "
      NOTES(2)=" Report questions to Dr. Xiongwu Wu  "
      NOTES(3)="             Email: wuxw@nhlbi.nih.gov "
      WRITE(UNIT)LX,LY,LZ               !1,2,3
      WRITE(UNIT)MODE                   ! 4
      WRITE(UNIT)MNX,MNY,MNZ            ! 5,6,7
      WRITE(UNIT)NX,NY,NZ               ! 8,9,10
      WRITE(UNIT)XL,YL,ZL               ! 11,12,13
      WRITE(UNIT)ALPHA,BETA,GAMMA       ! 14,15,16
      WRITE(UNIT)MAPC,MAPR,MAPS         ! 17,18,19
      WRITE(UNIT)AMIN,AMAX,AMEAN        ! 20,21,22
      WRITE(UNIT)ISPG,NSYMBT,LSKFLG     ! 23,24,25
      WRITE(UNIT)SKWMAT                 ! 26-34
      WRITE(UNIT)SKWTRN                 ! 35-37
      WRITE(UNIT)EXTRA                  ! 38-52
      WRITE(UNIT)MAPLABLE               ! 53
      WRITE(UNIT)MACHST                 ! 54
      WRITE(UNIT)ARMS                   ! 55
      WRITE(UNIT)NNOTE                  ! 56
      WRITE(UNIT)NOTES                  ! 57-256
! write data
      IF(MODE==0)THEN
        DO I=1,NDATA
          BDATAI=MAP%BDATA(I)
          WRITE(UNIT)BDATAI
        ENDDO
      ELSE IF(MODE==1)THEN
        DO I=1,NDATA
          IDATAI=MAP%IDATA(I)
          WRITE(UNIT)IDATAI
        ENDDO
      ELSE IF(MODE==2)THEN
        DO I=1,NDATA
          RDATAI=MAP%RDATA(I)
          WRITE(UNIT)RDATAI
        ENDDO
      ELSE IF(MODE==3)THEN
        DO I=1,NDATA
          CDATAI=MAP%CDATA(I)
          WRITE(UNIT)CDATAI
        ENDDO
      ELSE IF(MODE==4)THEN
        DO I=1,NDATA
          CDATAI=MAP%CDATA(I)
          WRITE(UNIT)CDATAI
        ENDDO
      ELSE IF(MODE==5)THEN
        DO I=1,NDATA
          BDATAI=MAP%BDATA(I)
          WRITE(UNIT)BDATAI
        ENDDO
      ENDIF
! print map information
      IF(OUTU > 0)THEN
      WRITE(OUTU,1000)MAP%NAME, UNIT  !map name to be written
      WRITE(OUTU,1010)LX,LY,LZ               !1,2,3
      WRITE(OUTU,1020)MODE                   ! 4
      WRITE(OUTU,1030)MNX,MNY,MNZ            ! 5,6,7
      WRITE(OUTU,1040)NX,NY,NZ               ! 8,9,10
      WRITE(OUTU,1050)XL,YL,ZL               ! 11,12,13
      WRITE(OUTU,1060)ALPHA,BETA,GAMMA       ! 14,15,16
      WRITE(OUTU,1070)MAPC,MAPR,MAPS         ! 17,18,19
      WRITE(OUTU,1080)AMIN,AMAX,AMEAN,ARMS       ! 20,21,22,55
      WRITE(OUTU,1090)ISPG,NSYMBT,LSKFLG,NNOTE     ! 23,24,25,56
      WRITE(OUTU,1110)SKWMAT                 ! 26-34
      WRITE(OUTU,1120)SKWTRN                 ! 35-37
      WRITE(OUTU,1140)MAPLABLE               ! 53
      WRITE(OUTU,1150)MACHST                 ! 54
      WRITE(OUTU,1160)(I,NOTES(I),I=1,NNOTE)   ! 57-256
      WRITE(OUTU,1210)NDATA                  ! DATA NUMBER
      ENDIF
1000  FORMAT(" map object: "/A/" is written to unit: ",I4)
1010  FORMAT(" LX, LY, LZ              = ",3I8)
1020  FORMAT(" MODE                    = ",I8)
1030  FORMAT(" MX, MY, MZ              = ",3I8)
1040  FORMAT(" NX, NY, NZ              = ",3I8)
1050  FORMAT(" XL, YL, ZL              = ",3F8.2)
1060  FORMAT(" ALPHA,BETA,GAMMA        = ",3F8.2)
1070  FORMAT(" MAPC, MAPR, MAPS        = ",3I8)
1080  FORMAT(" MIN,MAX,MEAN,STD        = ",4E12.4)
1090  FORMAT(" ISPG,NSYMBT,LSKFLG,NNOTE= ",4I8)
1110  FORMAT(" SKWMAT                  = ",3F8.2)
1120  FORMAT(" SKWTRN                  = ",3F8.2)
1130  FORMAT(" EXTRA                   = ",15F8.2)
1140  FORMAT(" MAPLABLE                = ",A)
1150  FORMAT(" MACHST                  = ",A)
1160  FORMAT(" NOTES ",I2,":"A80)
1210  FORMAT(" DATA POINT NUMBER       = ",I8)
      RETURN
      END SUBROUTINE WRTCCP4

      SUBROUTINE RDEMAP(FNAME,IDEMP,OUTU)
!-----------------------------------------------------------------------
!  read in map information and its density array
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!                        
      INTEGER IDEMP,OUTU
      CHARACTER*(*)FNAME
      type(EMAPOBJECT)  EMPOBJ
!
      INTEGER I
      _REAL_ RESO
!
      EMAPS(IDEMP)%name=FNAME
      IF(LEN_TRIM(FNAME)==0)THEN
        ! map will be created from input coordinates based on rigid mask
        RETURN      
      ELSE IF(INDEX(FNAME,'.ccp4')>0.or.INDEX(FNAME,'.CCP4')>0 .or.  &
         INDEX(FNAME,'.map')>0.or.INDEX(FNAME,'.MAP')>0)THEN
        CALL READCCP4(FNAME,IDEMP,OUTU)      
      ELSE IF(INDEX(FNAME,'.pdb')>0.or.INDEX(FNAME,'.PDB')>0)THEN
        DO I=1,NRIGID
          IF(EMRIGS(I)%MAPID==IDEMP)RESO=EMRIGS(I)%resolution
        ENDDO
        IF(RESO<1.0D-1)RESO=2.0D0
        CALL PDB2MAP(FNAME,IDEMP,RESO,OUTU)      
      ELSE
        write(outu,*)' Constraint map format is not supported:',fname
        stop
      ENDIF
      RETURN
      END SUBROUTINE RDEMAP

      SUBROUTINE READCCP4(FNAME,IDEMP,OUTU)
!_________________________________________________________________
!  read in CCP4 map information and its density array
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
     CHARACTER*(*)FNAME
     INTEGER UNIT,IDEMP,OUTU
!
      INTEGER*4 LX,LY,LZ,MODE,MNX,MNY,MNZ,NX,NY,NZ
      REAL*4 XL,YL,ZL,ALPHA,BETA,GAMMA
      INTEGER*4 MAPC,MAPR,MAPS
      REAL*4 AMAX,AMIN,AMEAN,ARMS
      INTEGER*4 ISPG,NSYMBT,LSKFLG,NNOTE
      REAL*4 SKWMAT(9),SKWTRN(3),EXTRA(15)
      CHARACTER*4 MAPLABLE,MACHST
      CHARACTER*80 NOTES(10)
      INTEGER*1 BDATAI
      INTEGER*2 IDATAI
      REAL*4 RDATAI
      COMPLEX(KIND=4) CDATAI
!
      INTEGER I,NDATA,alloc_err

!
      UNIT=99
      OPEN(UNIT,FILE=FNAME,ACCESS='STREAM',STATUS='OLD')
      !
      READ(UNIT)LX,LY,LZ               !1,2,3
      READ(UNIT)MODE                   ! 4
      READ(UNIT)MNX,MNY,MNZ            ! 5,6,7
      READ(UNIT)NX,NY,NZ               ! 8,9,10
      READ(UNIT)XL,YL,ZL               ! 11,12,13
      READ(UNIT)ALPHA,BETA,GAMMA       ! 14,15,16
      READ(UNIT)MAPC,MAPR,MAPS         ! 17,18,19
      READ(UNIT)AMIN,AMAX,AMEAN        ! 20,21,22
      READ(UNIT)ISPG,NSYMBT,LSKFLG     ! 23,24,25
      READ(UNIT)SKWMAT                 ! 26-34
      READ(UNIT)SKWTRN                 ! 35-37
      READ(UNIT)EXTRA                  ! 38-52
      READ(UNIT)MAPLABLE               ! 53
      READ(UNIT)MACHST                 ! 54
      READ(UNIT)ARMS                   ! 55
      READ(UNIT)NNOTE                  ! 56
      READ(UNIT)(NOTES(I),I=1,10)      ! 57-256
! write data
      NDATA=LX*LY*LZ
      IF(MODE==0)THEN
        if(associated(EMAPS(IDEMP)%BDATA))deallocate(EMAPS(IDEMP)%BDATA)
        allocate(EMAPS(IDEMP)%BDATA(NDATA),stat=alloc_err)
        if(alloc_err /= 0 ) write(6,*)"unable to allocate BDATA"
        DO I=1,NDATA
            READ(UNIT)BDATAI
            EMAPS(IDEMP)%BDATA(I)=BDATAI
        ENDDO
      ELSE IF(MODE==1)THEN
        if(associated(EMAPS(IDEMP)%IDATA))deallocate(EMAPS(IDEMP)%IDATA)
        allocate(EMAPS(IDEMP)%IDATA(NDATA),stat=alloc_err)
        if(alloc_err /= 0 ) write(6,*)"unable to allocate IDATA"
        DO I=1,NDATA
          READ(UNIT)IDATAI
          EMAPS(IDEMP)%IDATA(I)=IDATAI
        ENDDO
      ELSE IF(MODE==2)THEN
        if(associated(EMAPS(IDEMP)%RDATA))deallocate(EMAPS(IDEMP)%RDATA)
        allocate(EMAPS(IDEMP)%RDATA(NDATA),stat=alloc_err)
        if(alloc_err /= 0 ) write(6,*)"unable to allocate RDATA"
        DO I=1,NDATA
          READ(UNIT)RDATAI
          EMAPS(IDEMP)%RDATA(I)=RDATAI
        ENDDO
      ELSE IF(MODE==3)THEN
        if(associated(EMAPS(IDEMP)%CDATA))deallocate(EMAPS(IDEMP)%CDATA)
        allocate(EMAPS(IDEMP)%CDATA(NDATA),stat=alloc_err)
        if(alloc_err /= 0 ) write(6,*)"unable to allocate CDATA"
        DO I=1,NDATA
          READ(UNIT)CDATAI
          EMAPS(IDEMP)%CDATA(I)=CDATAI
        ENDDO
      ELSE IF(MODE==4)THEN
        if(associated(EMAPS(IDEMP)%CDATA))deallocate(EMAPS(IDEMP)%CDATA)
        allocate(EMAPS(IDEMP)%CDATA(NDATA),stat=alloc_err)
        if(alloc_err /= 0 ) write(6,*)"unable to allocate CDATA"
        DO I=1,NDATA
          READ(UNIT)CDATAI
          EMAPS(IDEMP)%CDATA(I)=CDATAI
        ENDDO
      ELSE IF(MODE==5)THEN
        if(associated(EMAPS(IDEMP)%BDATA))deallocate(EMAPS(IDEMP)%BDATA)
        allocate(EMAPS(IDEMP)%BDATA(NDATA),stat=alloc_err)
        if(alloc_err /= 0 ) write(6,*)"unable to allocate BDATA"
        DO I=1,NDATA
          READ(UNIT)BDATAI
          EMAPS(IDEMP)%BDATA(I)=BDATAI
        ENDDO
      ENDIF
      EMAPS(IDEMP)%LX=LX
      EMAPS(IDEMP)%LY=LY
      EMAPS(IDEMP)%LZ=LZ
      EMAPS(IDEMP)%MODE=MODE
      EMAPS(IDEMP)%MX=MNX
      EMAPS(IDEMP)%MY=MNY
      EMAPS(IDEMP)%MZ=MNZ
      EMAPS(IDEMP)%LX=NX
      EMAPS(IDEMP)%LY=NY
      EMAPS(IDEMP)%LZ=NZ
      EMAPS(IDEMP)%DX=XL/NX
      EMAPS(IDEMP)%DY=YL/NY
      EMAPS(IDEMP)%DZ=ZL/NZ
      EMAPS(IDEMP)%ALPHA=ALPHA
      EMAPS(IDEMP)%BETA=BETA
      EMAPS(IDEMP)%GAMMA=GAMMA
      EMAPS(IDEMP)%CX=(MNX+LX/2.0D0)*EMAPS(IDEMP)%DX
      EMAPS(IDEMP)%CY=(MNY+LY/2.0D0)*EMAPS(IDEMP)%DY
      EMAPS(IDEMP)%CZ=(MNZ+LZ/2.0D0)*EMAPS(IDEMP)%DZ
!     statistics
      CALL MAPSTAT(IDEMP,OUTU)
! print map information
      IF(OUTU > 0 )THEN
      WRITE(OUTU,1000)IDEMP,EMAPS(IDEMP)%NAME, UNIT  !map name to be written
      WRITE(OUTU,1010)LX,LY,LZ               !1,2,3
      WRITE(OUTU,1020)MODE                   ! 4
      WRITE(OUTU,1030)MNX,MNY,MNZ            ! 5,6,7
      WRITE(OUTU,1040)NX,NY,NZ               ! 8,9,10
      WRITE(OUTU,1050)XL,YL,ZL               ! 11,12,13
      WRITE(OUTU,1060)ALPHA,BETA,GAMMA       ! 14,15,16
      WRITE(OUTU,1070)MAPC,MAPR,MAPS         ! 17,18,19
      WRITE(OUTU,1080)AMIN,AMAX,AMEAN,ARMS       ! 20,21,22,55
      WRITE(OUTU,1090)ISPG,NSYMBT,LSKFLG,NNOTE     ! 23,24,25,56
      WRITE(OUTU,1110)SKWMAT                 ! 26-34
      WRITE(OUTU,1120)SKWTRN                 ! 35-37
      WRITE(OUTU,1130)EXTRA                  ! 38-52
      WRITE(OUTU,1210)NDATA                  ! DATA NUMBER
      WRITE(OUTU,1220)IDEMP                    ! DATA NUMBER
      ENDIF
1000  FORMAT(" ------------------EMAP IMAGE ",I4," INPUT ---------------" / &
            "map file: "/A/" is read from unit: ",I4)
1010  FORMAT(" LX, LY, LZ              = ",3I8)
1020  FORMAT(" MODE                    = ",I8)
1030  FORMAT(" MX, MY, MZ              = ",3I8)
1040  FORMAT(" NX, NY, NZ              = ",3I8)
1050  FORMAT(" XL, YL, ZL              = ",3F8.2)
1060  FORMAT(" ALPHA,BETA,GAMMA        = ",3F8.2)
1070  FORMAT(" MAPC, MAPR, MAPS        = ",3I8)
1080  FORMAT(" MIN,MAX,MEAN,STD        = ",4E12.4)
1090  FORMAT(" ISPG,NSYMBT,LSKFLG,NNOTE= ",4I8)
1110  FORMAT(" SKWMAT                  = ",3F8.2)
1120  FORMAT(" SKWTRN                  = ",3F8.2)
1130  FORMAT(" EXTRA                   = ",5F8.2)
1140  FORMAT(" MAPLABLE                = ",A)
1150  FORMAT(" MACHST                  = ",A)
1160  FORMAT(" NOTES ",I2,":"/A)
1210  FORMAT(" DATA POINT NUMBER       = ",I8)
1220  FORMAT(" ----------------------- END OF EMAP IMAGE ",I4, &
             "  -------------------------- ")
        CLOSE(UNIT)
      RETURN
      END SUBROUTINE READCCP4

      SUBROUTINE PDB2MAP(FNAME,IDEMP,RESO,OUTU)
!_________________________________________________________________
!  read in PDB coordinates to create a map object
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use ew_bspline,only:kbot,ktop,fill_bspline_0
!
      CHARACTER*(*)FNAME
      INTEGER UNIT,IDEMP,OUTU
      _REAL_ RESO
!
      CHARACTER*80 pdblin
      CHARACTER*4 CHAINI,CHAINP,RESIN,ATOMIN,SID,SIDP,RID,RIDP,ATOMNM
      _REAL_ XIN,YIN,ZIN,AMASSI,RHOI,DG,WIN
      INTEGER NAT,IAT,ISEQ
!
      INTEGER alloc_err
      _REAL_, allocatable, dimension(:)::tmp_crd,tmp_amass
!
      UNIT=99
      OPEN(UNIT,FILE=FNAME,ACCESS='SEQUENTIAL',STATUS='OLD')
      NAT=0
99    READ(UNIT,'(A)',END=120,ERR=120) PDBLIN
      PDBLIN=TRIM(PDBLIN)
      IF(LEN_TRIM(PDBLIN)<44)GOTO 99
      IF (PDBLIN(1:6).NE.'ATOM  '.AND.PDBLIN(1:6).NE.'HETATM') GOTO 99
      NAT=NAT+1
      GOTO 99
120   CONTINUE
      if(allocated(tmp_crd))deallocate(tmp_crd)
      allocate(tmp_crd(3*NAT),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate CRD !"
      if(allocated(tmp_amass))deallocate(tmp_amass)
      allocate(tmp_amass(NAT),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate AMASS !"
      CLOSE(UNIT)
      OPEN(UNIT,FILE=FNAME,ACCESS='SEQUENTIAL',STATUS='OLD')
      IAT=0
299   READ(UNIT,'(A)',END=320,ERR=320) PDBLIN
      PDBLIN=TRIM(PDBLIN)
      IF(LEN(PDBLIN)<44)GOTO 299
      IF (PDBLIN(1:6).NE.'ATOM  '.AND.PDBLIN(1:6).NE.'HETATM') GOTO 299
      READ(PDBLIN, &
          '(6X,I5,1X,A4,1X,A3,1X,A1,A5,3X,3F8.3)') &
               ISEQ,ATOMIN,RESIN,CHAINI,RID,XIN,YIN,ZIN     
!  determine atom mass
      ATOMNM=ADJUSTL(ATOMIN)
      IF(INDEX(ATOMNM,'CAL')==1)THEN
        AMASSI=40.0D0
      ELSE IF(INDEX(ATOMNM,'CL')==1)THEN
        AMASSI=35.5D0
      ELSE IF(INDEX(ATOMNM,'MG')==1)THEN
        AMASSI=24.3D0
      ELSE IF(INDEX(ATOMNM,'S')==1)THEN
        AMASSI=32.06D0
      ELSE IF(INDEX(ATOMNM,'O')==1)THEN
        AMASSI=16.0D0
      ELSE IF(INDEX(ATOMNM,'N')==1)THEN
        AMASSI=14.01D0
      ELSE IF(INDEX(ATOMNM,'P')==1)THEN
        AMASSI=31.0D0
      ELSE IF(INDEX(ATOMNM,'F')==1)THEN
        AMASSI=19.0D0
      ELSE IF(INDEX(ATOMNM,'C')==1)THEN
        AMASSI=12.01D0
      ELSE IF(INDEX(ATOMNM,'H')==1)THEN
        AMASSI=1.008D0
      ELSE 
        AMASSI=1.008D0
      ENDIF
      IAT=IAT+1
      tmp_amass(IAT)=AMASSI
      tmp_crd(3*IAT-2)=XIN
      tmp_crd(3*IAT-1)=YIN
      tmp_crd(3*IAT)=ZIN
      goto 299
320   CONTINUE
!     assign map object
      CALL COR2MAP(IDEMP,RESO,NAT,tmp_amass,tmp_crd,OUTU)
      if(allocated(tmp_crd))deallocate(tmp_crd)
      if(allocated(tmp_amass))deallocate(tmp_amass)
      CLOSE(UNIT)
      RETURN
      END SUBROUTINE PDB2MAP


      SUBROUTINE COR2MAP(IDEMP,RESO,NATC,AMASS,CRD,OUTU)
!_________________________________________________________________
!  Blur input coordinates to a map object
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      use ew_bspline,only:kbot,ktop,fill_bspline_0
!
      INTEGER IDEMP,NATC,OUTU
      _REAL_ AMASS(*)
      _REAL_ CRD(*)
      _REAL_ RESO
!
      real(kind=4),pointer:: rho(:)=> null()        !   real data
!
      INTEGER*4 LX,LY,LZ,MODE,MNX,MNY,MNZ,NX,NY,NZ
      REAL*4 XC,YC,ZC,ALPHA,BETA,GAMMA
      INTEGER*4 MAPC,MAPR,MAPS
      INTEGER*4 ISPG,NSYMBT,LSKFLG,NNOTE
      REAL*4 SKWMAT(9),SKWTRN(3),EXTRA(15)
      CHARACTER*4 MAPLABLE,MACHST
      CHARACTER*80 NOTES(10)
      REAL*4 RDATAI,XIN,YIN,ZIN,AMASSI,RHOI,DG,WIN
      INTEGER ISEQ,IAT
!
      INTEGER I,J,K,IN,JN,KN,M1,M2,M3,NDATA
      INTEGER IPT1,IPT2,IPT3,ITH1,ITH2,ITH3
      _REAL_ XMIN,YMIN,ZMIN,XMAX,YMAX,ZMAX
      _REAL_ FR1,FR2,FR3,W
      REAL*8 BS1(8),BS2(8),BS3(8)
      REAL*8 VAL0,VAL0A
      INTEGER alloc_err

!
      MODE=2
      XC=0.0D0
      YC=0.0D0
      ZC=0.0D0
      XMIN=1.0D9
      YMIN=1.0D9
      ZMIN=1.0D9
      XMAX=-1.0D9
      YMAX=-1.0D9
      ZMAX=-1.0D9
      DO IAT=1,NATC
        XIN=CRD(3*IAT-2)
        YIN=CRD(3*IAT-1)
        ZIN=CRD(3*IAT)
        XC=XC+XIN
        YC=YC+YIN
        ZC=ZC+ZIN
        IF(XMIN>XIN)XMIN=XIN
        IF(YMIN>YIN)YMIN=YIN
        IF(ZMIN>ZIN)ZMIN=ZIN
        IF(XMAX<XIN)XMAX=XIN
        IF(YMAX<YIN)YMAX=YIN
        IF(ZMAX<ZIN)ZMAX=ZIN
      ENDDO
      XC=XC/NATC
      YC=YC/NATC
      ZC=ZC/NATC
      DG=2.0D0*RESO/BORDER
      MNX=nint((XMIN-RESO)/DG-0.5d0)
      MNY=nint((YMIN-RESO)/DG-0.5d0)
      MNZ=nint((ZMIN-RESO)/DG-0.5d0)
      LX=nint((XMAX-XMIN+2.0D0*RESO)/DG+0.5d0)+1
      LY=nint((YMAX-YMIN+2.0D0*RESO)/DG+0.5d0)+1
      LZ=nint((ZMAX-ZMIN+2.0D0*RESO)/DG+0.5d0)+1
      NDATA=LX*LY*LZ
      if(associated(EMAPS(IDEMP)%RDATA))deallocate(EMAPS(IDEMP)%RDATA)
      allocate(EMAPS(IDEMP)%RDATA(NDATA),stat=alloc_err)
      if(alloc_err /= 0 ) write(6,*)"unable to allocate RDATA"
      RHO=>EMAPS(IDEMP)%RDATA
      RHO=0.0D0
!   generate map with b-spline
      DO IAT=1,NATC
!  determine atom mass
        AMASSI=AMASS(IAT)
        XIN=CRD(3*IAT-2)
        YIN=CRD(3*IAT-1)
        ZIN=CRD(3*IAT)
! calculate B-spline parameters
        FR1=XIN/DG
        FR2=YIN/DG
        FR3=ZIN/DG
        M1=ANINT(FR1-0.5D0)
        W=FR1-M1
        call fill_bspline_0(w,BORDER,BS1)
        M2=ANINT(FR2-0.5D0)
        W=FR2-M2
        call fill_bspline_0(w,BORDER,BS2)
        M3=ANINT(FR3-0.5D0)
        W=FR3-M3
        call fill_bspline_0(w,BORDER,BS3)
! Calculate B-spline interception
         K = M3-MNZ+1 - BORDER/2
         DO 410 ITH3 = 1,BORDER
           K=K+1
           IF(K<1)THEN
             IPT1=0
             ELSE IF(K>LZ)THEN
               IPT1=LX*LY*(LZ-1)
             ELSE
               IPT1=(K-1)*LY*LX
             ENDIF
             VAL0A =  AMASSI*BS3(ITH3)
!
             J = M2-MNY+1 - BORDER/2
!
             DO 420 ITH2 = 1,BORDER
               J=J+1
               IF(J<1)THEN
                 IPT2=IPT1
               ELSE IF(J>LY)THEN
                 IPT2=IPT1+LX*(LY-1)
               ELSE
                 IPT2=IPT1+(J-1)*LX
               ENDIF
!
               VAL0= VAL0A * BS2(ITH2)
!
               I = M1-MNX+1 - BORDER/2
!
               DO 430 ITH1 = 1,BORDER
                 I=I+1
                 IF(I<1)THEN
                   IPT3=IPT2+1
                 ELSE IF(I>LX)THEN
                   IPT3=IPT2+LX
                 ELSE
                   IPT3=IPT2+I
                 ENDIF
                 RHO(IPT3)=RHO(IPT3)+ VAL0 * BS1(ITH1)
!
430          CONTINUE
420        CONTINUE
410      CONTINUE
!
      ENDDO
!     assign map object
      
      EMAPS(IDEMP)%LX=LX
      EMAPS(IDEMP)%LY=LY
      EMAPS(IDEMP)%LZ=LZ
      EMAPS(IDEMP)%MODE=MODE
      EMAPS(IDEMP)%MX=MNX
      EMAPS(IDEMP)%MY=MNY
      EMAPS(IDEMP)%MZ=MNZ
      EMAPS(IDEMP)%DX=DG
      EMAPS(IDEMP)%DY=DG
      EMAPS(IDEMP)%DZ=DG
      EMAPS(IDEMP)%ALPHA=90.0D0
      EMAPS(IDEMP)%BETA=90.0D0
      EMAPS(IDEMP)%GAMMA=90.0D0
      NOTES(1)="          " &
         //" This map is created with the emap module of charmm "
      NOTES(2)="          " &
         //" Report questions to Dr. Xiongwu Wu  "
      NOTES(3)="          " &
        //"             Email: wuxw@nhlbi.nih.gov "
      EMAPS(IDEMP)%CX=XC
      EMAPS(IDEMP)%CY=YC
      EMAPS(IDEMP)%CZ=ZC
!     statistics
      CALL MAPSTAT(IDEMP,OUTU)
! print map information
      IF(OUTU > 0 )THEN
      WRITE(OUTU,1000)IDEMP  !map id
      WRITE(OUTU,1010)LX,LY,LZ               !1,2,3
      WRITE(OUTU,1020)MODE                   ! 4
      WRITE(OUTU,1030)MNX,MNY,MNZ            ! 5,6,7
      WRITE(OUTU,1050)LX*DG,LY*DG,LZ*DG            ! 5,6,7
      WRITE(OUTU,1080)EMAPS(IDEMP)%MIN,EMAPS(IDEMP)%MAX,EMAPS(IDEMP)%AVG,EMAPS(IDEMP)%STD       ! 20,21,22,55
      WRITE(OUTU,1210)NDATA                  ! DATA NUMBER
      WRITE(OUTU,1220)IDEMP                    ! DATA NUMBER
      ENDIF
1000  FORMAT(" ------------------EMAP ID ",I4," CREATED ---------------" )
1010  FORMAT(" LX, LY, LZ              = ",3I8)
1020  FORMAT(" MODE                    = ",I4)
1030  FORMAT(" MX, MY, MZ              = ",3I8)
1040  FORMAT(" NX, NY, NZ              = ",3I8)
1050  FORMAT(" XL, YL, ZL              = ",3F8.2)
1060  FORMAT(" ALPHA,BETA,GAMMA        = ",3F8.2)
1070  FORMAT(" MAPC, MAPR, MAPS        = ",3I8)
1080  FORMAT(" MIN,MAX,MEAN,STD        = ",4E12.4)
1090  FORMAT(" ISPG,NSYMBT,LSKFLG,NNOTE= ",4I8)
1110  FORMAT(" SKWMAT                  = ",3F8.2)
1120  FORMAT(" SKWTRN                  = ",3F8.2)
1130  FORMAT(" EXTRA                   = ",5F8.2)
1140  FORMAT(" MAPLABLE                = ",A)
1150  FORMAT(" MACHST                  = ",A)
1160  FORMAT(" NOTES ",I2,":"/A)
1210  FORMAT(" DATA POINT NUMBER       = ",I8)
1220  FORMAT(" ----------------------- END OF EMAP IMAGE ",I4, &
             "  -------------------------- ")
      RETURN
      END SUBROUTINE COR2MAP


      SUBROUTINE MAPSTAT(IDEMP,OUTU)
!_________________________________________________________________
!  calculate statistics of a map
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
!
      INTEGER IDEMP,OUTU
      _REAL_ RESO
      real(kind=4),pointer:: rho(:)=> null()        !   real data
!
      INTEGER NDATA,I
      REAL*8 AMAX,AMIN,AMEAN,ARMS,RHOI
!
      NDATA=EMAPS(IDEMP)%LX*EMAPS(IDEMP)%LY*EMAPS(IDEMP)%LZ
      RHO=>EMAPS(IDEMP)%RDATA
!     map statistics
      AMIN=1.0D9
      AMAX=-1.0D9
      AMEAN=0.0D0
      ARMS=0.0D0
      DO I=1,NDATA
        RHOI=RHO(I)
        IF(AMIN>RHOI)AMIN=RHOI
        IF(AMAX<RHOI)AMAX=RHOI
        AMEAN=AMEAN+RHOI
        ARMS=ARMS+RHOI*RHOI
      ENDDO
      AMEAN=AMEAN/NDATA
      ARMS=ARMS/NDATA-AMEAN*AMEAN
      IF(ARMS>0.0D0)ARMS=SQRT(ARMS)
      EMAPS(IDEMP)%MIN=AMIN
      EMAPS(IDEMP)%MAX=AMAX
      EMAPS(IDEMP)%AVG=AMEAN
      EMAPS(IDEMP)%STD=ARMS
      RETURN
      END SUBROUTINE MAPSTAT

!+ Read maps for constraint
subroutine emapforce(natom,enemap,amass,x,f)
!________________________________________________________________
!  Calculate map constraint forces
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
      _REAL_ ENEMAP,AMASS(*),X(*),F(*)
      INTEGER NATOM,I
      _REAL_ ENRIG
      ENEMAP=0.0D0
!  Loop over all active rigid domains
      DO I=1,NRIGID
        CALL RIGIDFORCE(NATOM,ENRIG,AMASS,X,F,EMRIGS(I))
        ENEMAP=ENEMAP+ENRIG
      ENDDO
      return
      end subroutine emapforce

      SUBROUTINE RIGIDENG(NATOM,ENRIG,AMASS,X,rigobj)
!_________________________________________________________________
!   This routine transform reference atoms to map frame and
!     calculate B-splie parameters and the energy and force
!     then transform the force back to simulation frame
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
!
      use ew_bspline,only:kbot,ktop,fill_bspline_1
!
#ifdef MPI
#  include "parallel.h"
   include 'mpif.h'
#endif
      INTEGER NATOM
      REAL*8 ENRIG,AMASS(*),X(*)
      type(EMAPRIGID) :: rigobj
      type(EMAPOBJECT) :: mapobj
      real(kind=4),pointer:: rho(:)=> null()        !   real data
      integer,pointer:: idx(:)=> null()        !   real data
      INTEGER NATC,MNX,MNY,MNZ,LX,LY,LZ
      REAL*8 TR(3),U(3,3),PGUID,RIG_F(3),RIG_T(3)
      REAL*8 RAVG,RHO0,CX,CY,CZ,GX,GY,GZ
      REAL*8 BS1(8),BS2(8),BS3(8)
      REAL*8 DBS1(8),DBS2(8),DBS3(8)
      REAL*8 XI,YI,ZI,XJ,YJ,ZJ
      REAL*8 AMASSI,RHOI,ENRIGM
      REAL*8 CFACT,CFACTX,CFACTY,CFACTZ,FR1,FR2,FR3,W
      REAL*8 VAL0,VAL1,VAL2,VAL3,VAL0A,VAL1A,VAL2A,VAL3A
      INTEGER ID,IAT,IATX,IATY,IATZ,I,J,K,M1,M2,M3,NSTA,NEND
      INTEGER ITH1,ITH2,ITH3,IPT1,IPT2,IPT3,IERR
      _REAL_ TEMP1(10),TEMP2(10)
!
#ifdef MPI
   if ( mpi_orig ) then
      nsta = 1
      nend = natom
   else
      nsta = iparpt(mytaskid) + 1
      nend = iparpt(mytaskid+1)
   end if
#else
   nsta = 1
   nend = natom
#endif
      mapobj=emaps(rigobj%mapid)
      MNX=MAPOBJ%MX
      MNY=MAPOBJ%MY
      MNZ=MAPOBJ%MZ
      LX=MAPOBJ%LX
      LY=MAPOBJ%LY
      LZ=MAPOBJ%LZ
      GX=MAPOBJ%DX
      GY=MAPOBJ%DY
      GZ=MAPOBJ%DZ
      CX=MAPOBJ%CX
      CY=MAPOBJ%CY
      CZ=MAPOBJ%CZ
      RHO=>MAPOBJ%RDATA
      RHO0=mapobj%AVG
      CFACT=rigobj%fcons/mapobj%std
      CFACTX=CFACT/GX
      CFACTY=CFACT/GY
      CFACTZ=CFACT/GZ
      NATC=RIGOBJ%NATC
      U=RIGOBJ%ROT
      TR=RIGOBJ%TRAN
      IDX=>RIGOBJ%IDXATOM
      ENRIG=0.0D0
!  Calculate selected atoms
      DO ID=1,NATC
        IAT=IDX(ID)
        IF(IAT<NSTA.OR.IAT>NEND)CYCLE
        IATX=3*IAT-2
        IATY=IATX+1
        IATZ=IATY+1
        AMASSI=AMASS(IAT)
        XI=X(IATX)-CX-TR(1)
        YI=X(IATY)-CY-TR(2)
        ZI=X(IATZ)-CZ-TR(3)
! transform to map frame
        XJ=XI*U(1,1)+YI*U(2,1)+ZI*U(3,1)+CX
        YJ=XI*U(1,2)+YI*U(2,2)+ZI*U(3,2)+CY
        ZJ=XI*U(1,3)+YI*U(2,3)+ZI*U(3,3)+CZ
! calculate B-spline parameters
        FR1=XJ/GX
        FR2=YJ/GY
        FR3=ZJ/GZ
        M1=ANINT(FR1-0.5D0)
        W=FR1-M1
        call fill_bspline_1(w,BORDER,BS1,dBS1)
        M2=ANINT(FR2-0.5D0)
        W=FR2-M2
        call fill_bspline_1(w,BORDER,BS2,dBS2)
        M3=ANINT(FR3-0.5D0)
        W=FR3-M3
        call fill_bspline_1(w,BORDER,BS3,dBS3)
! Calculate B-spline interception
        K = M3-MNZ+1 - BORDER/2
        DO 100 ITH3 = 1,BORDER
           K=K+1
           IF(K*(LZ-K)<=0)GOTO 100
           IPT1=(K-1)*LY*LX
            VAL0A = AMASSI  * BS3(ITH3)
!
            J = M2-MNY+1 - BORDER/2
!
            DO 200 ITH2 = 1,BORDER
               J=J+1
               IF(J*(LY-J)<=0)GOTO 200
               IPT2=IPT1+(J-1)*LX
!
               VAL0= VAL0A * BS2(ITH2)
!
               I = M1-MNX+1 - BORDER/2
!
               DO 300 ITH1 = 1,BORDER
                  I=I+1
                  IF(I*(LX-I)<=0)GOTO 300
                  IPT3=IPT2+I
                  RHOI=RHO(IPT3)-RHO0
                  ENRIG = ENRIG + VAL0 * RHOI * BS1(ITH1)
!
300          CONTINUE
200        CONTINUE
100      CONTINUE
!
      ENDDO
      ENRIG=-CFACT*ENRIG
      ENRIGM=ENRIG
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=ENRIGM
# ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,temp1,1,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          ENRIGM=TEMP1(1)
#else
          CALL MPI_ALLREDUCE(TEMP1,TEMP2,1, &
          MPI_DOUBLE_PRECISION,MPI_SUM,COMMSANDER,IERR)
          ENRIGM=TEMP2(1)
# endif
        ENDIF
#endif
      rigobj%energy=enrigm
      RETURN
      END SUBROUTINE RIGIDENG

      SUBROUTINE RIGIDFORCE(NATOM,ENRIG,AMASS,X,F,rigobj)
!_________________________________________________________________
!   This routine transform reference atoms to map frame and
!     calculate B-splie parameters and the energy and force
!     then transform the force back to simulation frame
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
!
      use ew_bspline,only:kbot,ktop,fill_bspline_1
!
#ifdef MPI
#  include "parallel.h"
   include 'mpif.h'
#endif
      INTEGER NATOM
      REAL*8 ENRIG,AMASS(*),X(*),F(*)
      type(EMAPRIGID) :: rigobj
      type(EMAPOBJECT) :: mapobj
      real(kind=4),pointer:: rho(:)=> null()        !   real data
      integer,pointer:: idx(:)=> null()        !   real data
      INTEGER NATC,MNX,MNY,MNZ,LX,LY,LZ
      REAL*8 TR(3),U(3,3),PGUID,RIG_F(3),RIG_T(3)
      REAL*8 RAVG,RHO0,CX,CY,CZ,GX,GY,GZ
      REAL*8 BS1(8),BS2(8),BS3(8)
      REAL*8 DBS1(8),DBS2(8),DBS3(8)
      REAL*8 XI,YI,ZI,XJ,YJ,ZJ
      REAL*8 AMASSI,RHOI,ENRIGM,F1,F2,F3,FX,FY,FZ
      REAL*8 CFACT,CFACTX,CFACTY,CFACTZ,FR1,FR2,FR3,W
      REAL*8 rfx,rfy,rfz,rtx,rty,rtz
      REAL*8 VAL0,VAL1,VAL2,VAL3,VAL0A,VAL1A,VAL2A,VAL3A
      INTEGER ID,IAT,IATX,IATY,IATZ,I,J,K,M1,M2,M3,NSTA,NEND
      INTEGER ITH1,ITH2,ITH3,IPT1,IPT2,IPT3,IERR
      _REAL_ TEMP1(10),TEMP2(10)
!
#ifdef MPI
   if ( mpi_orig ) then
      nsta = 1
      nend = natom
   else
      nsta = iparpt(mytaskid) + 1
      nend = iparpt(mytaskid+1)
   end if
#else
   nsta = 1
   nend = natom
#endif
      mapobj=emaps(rigobj%mapid)
      MNX=MAPOBJ%MX
      MNY=MAPOBJ%MY
      MNZ=MAPOBJ%MZ
      LX=MAPOBJ%LX
      LY=MAPOBJ%LY
      LZ=MAPOBJ%LZ
      GX=MAPOBJ%DX
      GY=MAPOBJ%DY
      GZ=MAPOBJ%DZ
      CX=MAPOBJ%CX
      CY=MAPOBJ%CY
      CZ=MAPOBJ%CZ
      RHO=>MAPOBJ%RDATA
      RHO0=mapobj%AVG
      CFACT=rigobj%fcons/mapobj%std
      CFACTX=CFACT/GX
      CFACTY=CFACT/GY
      CFACTZ=CFACT/GZ
      NATC=RIGOBJ%NATC
      U=RIGOBJ%ROT
      TR=RIGOBJ%TRAN
      IDX=>RIGOBJ%IDXATOM
      ENRIG=0.0D0
      rfx=0.0D0
      rfy=0.0D0
      rfz=0.0D0
      rtx=0.0D0
      rty=0.0D0
      rtz=0.0D0
!  Calculate selected atoms
      DO ID=1,NATC
        IAT=IDX(ID)
        IF(IAT<NSTA.OR.IAT>NEND)CYCLE
        IATX=3*IAT-2
        IATY=IATX+1
        IATZ=IATY+1
        AMASSI=AMASS(IAT)
        XI=X(IATX)-CX-TR(1)
        YI=X(IATY)-CY-TR(2)
        ZI=X(IATZ)-CZ-TR(3)
! transform to map frame
        XJ=XI*U(1,1)+YI*U(2,1)+ZI*U(3,1)+CX
        YJ=XI*U(1,2)+YI*U(2,2)+ZI*U(3,2)+CY
        ZJ=XI*U(1,3)+YI*U(2,3)+ZI*U(3,3)+CZ
! calculate B-spline parameters
        FR1=XJ/GX
        FR2=YJ/GY
        FR3=ZJ/GZ
        M1=ANINT(FR1-0.5D0)
        W=FR1-M1
        call fill_bspline_1(w,BORDER,BS1,dBS1)
        M2=ANINT(FR2-0.5D0)
        W=FR2-M2
        call fill_bspline_1(w,BORDER,BS2,dBS2)
        M3=ANINT(FR3-0.5D0)
        W=FR3-M3
        call fill_bspline_1(w,BORDER,BS3,dBS3)
! Calculate B-spline interception
        F1 = 0.0D0 
        F2 = 0.0D0 
        F3 = 0.0D0 
        K = M3-MNZ+1 - BORDER/2
        DO 100 ITH3 = 1,BORDER
           K=K+1
           IF(K*(LZ-K)<=0)GOTO 100
           IPT1=(K-1)*LY*LX
            VAL0A = AMASSI  * BS3(ITH3)
            VAL1A = AMASSI * CFACTX * BS3(ITH3)
            VAL2A = AMASSI * CFACTY * BS3(ITH3)
            VAL3A = AMASSI * CFACTZ * DBS3(ITH3)
!
            J = M2-MNY+1 - BORDER/2
!
            DO 200 ITH2 = 1,BORDER
               J=J+1
               IF(J*(LY-J)<=0)GOTO 200
               IPT2=IPT1+(J-1)*LX
!
               VAL0= VAL0A * BS2(ITH2)
               VAL1= VAL1A * BS2(ITH2)
               VAL2= VAL2A * DBS2(ITH2)
               VAL3= VAL3A * BS2(ITH2)

!
               I = M1-MNX+1 - BORDER/2
!
               DO 300 ITH1 = 1,BORDER
                  I=I+1
                  IF(I*(LX-I)<=0)GOTO 300
                  IPT3=IPT2+I
                  RHOI=RHO(IPT3)-RHO0
                  ENRIG = ENRIG + VAL0 * RHOI * BS1(ITH1)
!                                   force is negative of grad
                  F1 = F1 - VAL1 * RHOI * DBS1(ITH1)
                  F2 = F2 - VAL2 * RHOI * BS1(ITH1)
                  F3 = F3 - VAL3 * RHOI * BS1(ITH1)
!
300          CONTINUE
200        CONTINUE
100      CONTINUE
!
         FX = (U(1,1)*F1+U(1,2)*F2+U(1,3)*F3)
         FY = (U(2,1)*F1+U(2,2)*F2+U(2,3)*F3)
         FZ = (U(3,1)*F1+U(3,2)*F2+U(3,3)*F3)
         F(IATX) = F(IATX)-FX
         F(IATY) = F(IATY)-FY
         F(IATZ) = F(IATZ)-FZ
        ! Calculate force and torque on the rigid
         rfx=rfx+fx
         rfy=rfy+fy
         rfz=rfz+fz
         rtx=rtx+yi*fz-zi*fy
         rty=rty+zi*fx-xi*fz
         rtz=rtz+xi*fy-yi*fx
      ENDDO
      ENRIG=-CFACT*ENRIG
      ENRIGM=ENRIG
#ifdef MPI
        IF(SANDERSIZE > 1)THEN
!  Combining all node results
!
          TEMP1(1)=ENRIGM
          TEMP1(2)=RFX
          TEMP1(3)=RFY
          TEMP1(4)=RFZ
          TEMP1(5)=RTX
          TEMP1(6)=RTY
          TEMP1(7)=RTZ
# ifdef USE_MPI_IN_PLACE
          call mpi_allreduce(MPI_IN_PLACE,temp1,7,MPI_DOUBLE_PRECISION,MPI_SUM,commsander,ierr)
          ENRIGM=TEMP1(1)
          RFX=TEMP1(2)
          RFY=TEMP1(3)
          RFZ=TEMP1(4)
          RTX=TEMP1(5)
          RTY=TEMP1(6)
          RTZ=TEMP1(7)
#else
          CALL MPI_ALLREDUCE(TEMP1,TEMP2,7, &
          MPI_DOUBLE_PRECISION,MPI_SUM,COMMSANDER,IERR)
          ENRIGM=TEMP2(1)
          RFX=TEMP2(2)
          RFY=TEMP2(3)
          RFZ=TEMP2(4)
          RTX=TEMP2(5)
          RTY=TEMP2(6)
          RTZ=TEMP2(7)
# endif
        ENDIF
        !write(*,*)'rank,rfx,rtx=',sanderrank,rigobj%rigid,nsta,nend,rfx,rtx
#endif
      rigobj%energy=enrigm
      rigobj%force(1)=rfx
      rigobj%force(2)=rfy
      rigobj%force(3)=rfz
      rigobj%torq(1)=rtx
      rigobj%torq(2)=rty
      rigobj%torq(3)=rtz
      ! If the map object is movable, remove net force to constraint atoms
      if(rigobj%movable)then
        rfx=-rfx/rigobj%mass
        rfy=-rfy/rigobj%mass
        rfz=-rfz/rigobj%mass
        DO ID=1,NATC
          IAT=IDX(ID)
          IF(IAT<NSTA.OR.IAT>NEND)CYCLE
          IATX=3*IAT-2
          IATY=IATX+1
          IATZ=IATY+1
          AMASSI=AMASS(IAT)
          F(IATX) = F(IATX)-AMASSI*RFX
          F(IATY) = F(IATY)-AMASSI*RFY
          F(IATZ) = F(IATZ)-AMASSI*RFZ
        enddo
      endif
      RETURN
      END SUBROUTINE RIGIDFORCE

      SUBROUTINE EMAPGUID(IDRIG,IDEMP,AMASS,NAT,IDX)
!_________________________________________________________________
!     Initialize constraint parameters
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
!
!
!
      INTEGER IDRIG,IDEMP,NAT,IDX(*)
      _REAL_ AMASS(*)
      INTEGER IAT,I,NDATA
      _REAL_ TMASS,DELT
!
!
      TMASS=0.0D0
!  Loop over all active rigid domains
      DO I=1,NAT
        IAT=IDX(I)
        TMASS=TMASS+AMASS(IAT)
      ENDDO
!  Statistics of map object
      NDATA=EMAPS(IDEMP)%LX*EMAPS(IDEMP)%LY*EMAPS(IDEMP)%LZ
      IF(DELT>0.0D0)DELT=SQRT(DELT)
      RETURN
      END SUBROUTINE EMAPGUID



      FUNCTION rigeng(pos)
!_______________________________________________________________________
!     emap energy for rigid fitting
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!_________________________________________________________________
!
      use memory_module, only: xx =>x
#  include "memory.h"
      _REAL_ rigeng,rigengi
      _REAL_ pos(6)
!
      integer i
      
       do i=1,3
         rigmin%tran(i)=pos(i)
         rigmin%eular(i)=pos(i+3)
       enddo
       call ROTEULAR(rigmin%rot,pos(4),pos(5),pos(6))
       CALL RIGIDENG(NATOM,RIGENGI,xx(lmass),xx(lcrd),RIGMIN)
       RIGENG=RIGMIN%ENERGY
      return
      end function rigeng

      FUNCTION rigengd(pos,dpos)
!_______________________________________________________________________
!     emap energy and derivatives to x,y,z, and phi,psi,theta
!                    By Xiongwu Wu, wuxw@nhlbi.nih.gov
!
!     dEmap/dr=-f    r=(x,y,z)
!     dEmap/da=(dEmap/drc)(drc/da)=-T(d(A*r0)/da)=-T(dA/da)(A-)rc
!     a=(phi,psi,theta)=(f,y,s), rc=(1,1,1), T=(Tx,Ty, Tz)
!_________________________________________________________________
!
!

      use memory_module, only: xx =>x
#  include "memory.h"
      _REAL_ rigengd,rigengi
      _REAL_ pos(6),dpos(6)
!
      integer i
      _REAL_ a(3,3),daf(3,3),day(3,3),das(3,3),b(3),t(3),c(3)
      
       do i=1,3
         rigmin%tran(i)=pos(i)
         rigmin%eular(i)=pos(i+3)
       enddo
       call rotdeular(rigmin%rot,pos(4),pos(5),pos(6),daf,day,das)
       CALL RIGIDFORCE(NATOM,RIGENGI,xx(lmass),xx(lcrd),xx(lforce),RIGMIN)
       RIGENGD=RIGMIN%ENERGY
       do i=1,3
         dpos(i)=-rigmin%force(i)
       enddo
       a=rigmin%rot
       t=rigmin%torq
      b(1)=a(1,1)+a(2,1)+a(3,1)
      b(2)=a(1,2)+a(2,2)+a(3,2)
      b(3)=a(1,3)+a(2,3)+a(3,3)
      c=matmul(daf,b)
      dpos(4)=-dot_product(t,c)
      c=matmul(day,b)
      dpos(5)=-dot_product(t,c)
      c=matmul(das,b)
      dpos(6)=-dot_product(t,c)
      return
      end function rigengd


      subroutine amotry(simtry,p,y,psum,mp,np,ndim,funk,ihi,fac)
!_____________________________________________
!     Trypoint converted from F77 numerical recipes: amotry.for
!______________________________________________
      INTEGER ihi,mp,ndim,np,NMAX  
      _REAL_ simtry,fac,p(mp,np),psum(np),y(mp),funk
      PARAMETER (NMAX=20)  
      EXTERNAL funk  
!    USES funk  
      INTEGER j  
      _REAL_ fac1,fac2,ytry,ptry(NMAX)  
      fac1=(1.-fac)/ndim  
      fac2=fac1-fac  
      do 11 j=1,ndim  
        ptry(j)=psum(j)*fac1-p(ihi,j)*fac2  
11    continue  
      ytry=funk(ptry)  
      if (ytry.lt.y(ihi)) then  
        y(ihi)=ytry  
        do 12 j=1,ndim  
          psum(j)=psum(j)-p(ihi,j)+ptry(j)  
          p(ihi,j)=ptry(j)  
12      continue  
      endif  
      simtry=ytry  
      return  
      END subroutine amotry
 

      SUBROUTINE simplex(p,y,mp,np,ndim,ftol,funk,iter)  
!_____________________________________________
!     Downhill simplex method converted from F77 numerical recipes: amoeba.for
!______________________________________________
      INTEGER iter,mp,ndim,np,NMAX,ITMAX  
      _REAL_ ftol,p(mp,np),y(mp),funk
      PARAMETER (NMAX=20,ITMAX=10000)  
      EXTERNAL funk
!    USES simtry,funk  
      INTEGER i,ihi,ilo,inhi,j,m,n  
      _REAL_ rtol,sum,swap,ysave,ytry,psum(NMAX)
      iter=0  
1     do 12 n=1,ndim  
        sum=0.  
        do 11 m=1,ndim+1  
          sum=sum+p(m,n)  
11      continue  
        psum(n)=sum  
12    continue  
2     ilo=1  
      if (y(1).gt.y(2)) then  
        ihi=1  
        inhi=2  
      else  
        ihi=2  
        inhi=1  
      endif  
      do 13 i=1,ndim+1  
        if(y(i).le.y(ilo)) ilo=i  
        if(y(i).gt.y(ihi)) then  
          inhi=ihi  
          ihi=i  
        else if(y(i).gt.y(inhi)) then  
          if(i.ne.ihi) inhi=i  
        endif  
13    continue  
      rtol=2.*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))  
      if (rtol.lt.ftol) then  
        swap=y(1)  
        y(1)=y(ilo)  
        y(ilo)=swap  
        do 14 n=1,ndim  
          swap=p(1,n)  
          p(1,n)=p(ilo,n)  
          p(ilo,n)=swap  
14      continue  
        return  
      endif  
      if (iter.ge.ITMAX) then
         return
      endif
      iter=iter+2  
      call amotry(ytry,p,y,psum,mp,np,ndim,funk,ihi,-1.0d0)  
      if (ytry.le.y(ilo)) then  
        call amotry(ytry,p,y,psum,mp,np,ndim,funk,ihi,2.0d0)  
      else if (ytry.ge.y(inhi)) then  
        ysave=y(ihi)  
        call amotry(ytry,p,y,psum,mp,np,ndim,funk,ihi,0.5d0)  
        if (ytry.ge.ysave) then  
          do 16 i=1,ndim+1  
            if(i.ne.ilo)then  
              do 15 j=1,ndim  
                psum(j)=0.5*(p(i,j)+p(ilo,j))  
                p(i,j)=psum(j)  
15            continue  
              y(i)=funk(psum)  
            endif  
16        continue  
          iter=iter+ndim  
          goto 1  
        endif  
      else  
        iter=iter-1  
      endif  
      goto 2  
      RETURN
      END SUBROUTINE simplex

  

end module emap
