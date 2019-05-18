#include "dprec.fh"
!
!********************************************************************
!
!  Davidson diagonalization
!
!  Converted to fortran 90 (from fortran 77)
!  and interfaced with sqm(AmberTools) by
!  Kirill A. Velizhanin (kirill@lanl.gov)
!
!********************************************************************	
!
!********************************************************************
!   
   subroutine davidson(qm2_params,qmmm_nml,qmmm_mpi,cosmo_c_struct,qm2_struct,qm2ds, qmmm_struct)
   use qmmm_module,only:qm2_structure, qmmm_mpi_structure !cml-test
   use qm2_davidson_module
   use qmmm_struct_module, only : qmmm_struct_type
   use cosmo_C, only : cosmo_C_structure
   use qm2_params_module,  only : qm2_params_type
   use qmmm_nml_module   , only : qmmm_nml_type

   implicit none
   type(qmmm_nml_type),intent(inout) :: qmmm_nml
   type(qm2_params_type),intent(inout) :: qm2_params
   type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi 
   type(cosmo_C_structure), intent (inout) :: cosmo_c_struct
   type(qm2_structure),intent(inout) :: qm2_struct
   type(qmmm_struct_type), intent(inout) :: qmmm_struct
   type(qm2_davidson_structure_type), intent(inout) :: qm2ds

   _REAL_ random,rranset,rranf
   _REAL_ ferr1,ferr2,fn
   _REAL_ f,f1,f2,ddot,f0,f11

   integer i,j,lprint
   integer iseed,ip,ih,j0,Mx0,imin,one,iloops
   integer kflag,nd,nd1,nd1_old,j1

   integer Lt,Mb
!
!--------------------------------------------------------------------
!
!  Initialize Davidson
!
!--------------------------------------------------------------------
!
   Lt=qm2ds%Nb*(qm2ds%Nb+1)/2
   Mb=qm2ds%Nb**2
   one=1 
   lprint=qm2ds%verbosity       

   nd=qm2ds%nd 

   iloops=0

   if(qm2ds%dav_guess==0) then
      qm2ds%istore=0 ! overwriting qm2ds%istore to not use guess
   end if

! Split Mx into batches of j1 size
   j1=nd/qm2ds%idavfrac
   j1=min(j1,qm2ds%Mx)
   nd1=min(j1+2,nd/2)
   if (qm2ds%mdflag.ge.0) j1=qm2ds%Mx ! No shift is allowed for MD points
      

   if (qm2ds%irflag.gt.1) then
!     Calculations will start from irflag state for ceo
      j0=qm2ds%irflag-1
      !write(6,*)"begin at ", j0,qm2ds%Mx,qm2ds%irflag
      if (j0.gt.qm2ds%Mx) then
         write(6,*) 'Looks like irflag= ', qm2ds%irflag
         write(6,*) 'Show that states Mx= ',qm2ds%Mx
         write(6,*) 'are already calculated, exiting'

         stop 
      end if

!     Read irflag-1 states 
      open (qm2ds%modes_unit,file=trim(qm2ds%modes_b),form='unformatted',status='old')
      open (qm2ds%ee_unit,file=trim(qm2ds%ee_b),status='old')

      do j=1,j0
         read (qm2ds%modes_unit) (qm2ds%v0(i,j),i=1,qm2ds%Nrpa)
         read (qm2ds%ee_unit,*) qm2ds%e0(j)
      enddo
   
      close(qm2ds%modes_unit)
      close(qm2ds%ee_unit)

      if (j0.eq.qm2ds%Mx) goto 70

   else    ! Assume start from the beginning 
      j0=0
   end if

   qm2ds%Mj=j0

   if (lprint.gt.1) then
      write(6,*) 'Davidson parameters'
      write(6,*) 'mdflag=',qm2ds%mdflag
      write(6,*) 'irflag=',qm2ds%irflag
      write(6,*) 'M4=',qm2ds%Ncis
      write(6,*) 'Mx=',qm2ds%Mx
      write(6,*) 'Mj=',qm2ds%Mj
      write(6,*) 'nd=',nd 
      write(6,*) 'nd1=',nd1
      write(6,*) 'j1=',j1
      write(6,*)'qm2ds%istore=',qm2ds%istore
      write(6,*)'qm2ds%istore_M=',qm2ds%istore_M
   end if


   if(qm2ds%istore.gt.0) then ! MD point only!!!!	 
!     recover excited state vectors from AO representation in v2
!     recovered state vectors are put to v0
      do j=1,qm2ds%istore
         call site2mo(qm2ds,qm2ds%rrwork,qm2ds%v2(1,j),qm2ds%v0(1,j))
      end do
   endif

! kav: the lines below are commented out since they are not present in original davidson
! 
! CML TEST 7/15/12
!
!--------------------------------------------------------------------
!  Begin big loop
!--------------------------------------------------------------------
! find many vectors:
10 continue
   if (j0.lt.qm2ds%Mx) then 

   iloops=iloops+1
   j1=min(j1,(qm2ds%Mx-j0))
   nd1=min(j1+2,nd/2)

   if (lprint.gt.3) then
      write(6,*)
      write(6,*) 'Davidson batch ', iloops
      write(6,*) 'So far found ', j0, ' states'
      write(6,*) 'out of requested ',qm2ds%Mx, ' states'
      write(6,*) 'This batch will seek',j1,' vectors'

      if(j0.gt.0) write(6,*) 'Shift is',qm2ds%fs+qm2ds%e0(j0),' eV'
   end if 

! order quasidiagonal:
   i=0
   do ip=1,qm2ds%Np
    do ih=1,qm2ds%Nh
       i=i+1
       qm2ds%rrwork(i)=qm2ds%ehf(ih+qm2ds%Np)-qm2ds%ehf(ip)
    end do
   end do

   call rrdpsort(qm2ds%rrwork,qm2ds%Ncis,qm2ds%ix,1) 
!  Account for found vectors
   do j=1,j0
    do i=1,qm2ds%Ncis
       qm2ds%rrwork(i)=qm2ds%rrwork(i)+(qm2ds%fs+qm2ds%e0(j0)) &
          *abs(qm2ds%v0(i,j)**2-qm2ds%v0(qm2ds%Ncis+i,j)**2)
    end do
   end do

! try to find vector new vectors in the batch:
   call davidson0(qm2_params,qmmm_nml,qmmm_mpi, cosmo_c_struct, qm2_struct,qm2ds,&
      qmmm_struct,qm2ds%Ncis,lprint,qm2ds%ftol0,qm2ds%ftol1,qm2ds%ferr, &
      qm2ds%Np,qm2ds%Nh,j0,j1, &
      qm2ds%e0,qm2ds%v0,kflag,qm2ds%ix, &
      qm2ds%rrwork(1),qm2ds%rrwork(2*qm2ds%Ncis+1), &
      qm2ds%rrwork(4*qm2ds%Ncis+1), &
      nd,nd1,qm2ds%vexp1,qm2ds%vexp,qm2ds%ray,qm2ds%rayv,qm2ds%rayvL, &
      qm2ds%rayvR,qm2ds%raye,qm2ds%raye1, &
      qm2ds%ray1,qm2ds%ray1a,qm2ds%ray2,qm2ds%idav,qm2ds%istore)

! Printing out found eigenvalues, error and tolerance
   if(lprint>0) then
      write(6,*)' i, e0(i), ferr(i), ftol0'
      do i=1,j0
         write(6,111) i,' +++ ',qm2ds%e0(i),qm2ds%ferr(i),qm2ds%ftol0
      end do
      write(6,*)'-------------------------------------------------'
      write(6,*) 
   end if

   if (qm2ds%mdflag.le.-3) qm2ds%Mj=j0
111   format (i3,a,g24.16,2(' ',e8.2))
   call flush(6)

   
! Write vectors only for BIG sizes in the case of crash/restart       	  
   
   if(qm2ds%mdflag.lt.0.and.qm2ds%Nb.gt.100) then
      open (qm2ds%modes_unit,file=trim(qm2ds%modes_b),form='unformatted')
      open (qm2ds%ee_unit,file=trim(qm2ds%ee_b))

      do j=1,j0
         write (qm2ds%modes_unit) (qm2ds%v0(i,j),i=1,qm2ds%Nrpa)
         write (qm2ds%ee_unit,*)  qm2ds%e0(j)
      end do

      close(qm2ds%modes_unit)
      close(qm2ds%ee_unit)
   end if

   goto 10
   end if
       
!
!--------------------------------------------------------------------
!
! End big loop
!
!--------------------------------------------------------------------
!
 70   continue

! end find many vectors
!--------------------------------------------------------------------
!
   if(qm2ds%mdflag.ne.0) then   ! initial MD point, store vector
!  normalize vectors
      do j=1,qm2ds%Mx
         fn=0
! CIS
         if(qm2ds%idav==1) then
            do i=1,qm2ds%Ncis
               qm2ds%v0(i+qm2ds%Ncis,j)=0.0
            end do
         end if
! end CIS

         fn=ddot(qm2ds%Ncis,qm2ds%v0(1,j),one,qm2ds%v0(1,j),one) &
            -ddot(qm2ds%Ncis,qm2ds%v0(qm2ds%Ncis+1,j),one, &
            qm2ds%v0(qm2ds%Ncis+1,j),one)


         f=1/sqrt(abs(fn))
         call dscal(qm2ds%Nrpa,f,qm2ds%v0(1,j),one)

      end do

!     write to the hard disk
      call rrdpsort(qm2ds%e0,qm2ds%Mx,qm2ds%kx,2)

      if (qm2ds%mdflag.lt.0) then
         open (qm2ds%modes_unit,file=trim(qm2ds%modes_b),form='unformatted')
         open (qm2ds%ee_unit,file=trim(qm2ds%ee_b))
!
         do j=1,qm2ds%Mx
            write (qm2ds%modes_unit) (qm2ds%v0(i,qm2ds%kx(j)),i=1,qm2ds%Nrpa)
            write (qm2ds%ee_unit,*)  qm2ds%e0(j)
         end do
 
         close(qm2ds%modes_unit)
         close(qm2ds%ee_unit)
      end if
   end if

   if (qm2ds%mdflag.ge.0) then ! MD point only!!!!
!     store excited state vectors in AO representation in v2
!     qm2ds%istore=min(Mx,Mx_ss)
      if (qm2ds%istore.le.qm2ds%Mx) qm2ds%istore=qm2ds%Mx
      if (qm2ds%istore_M.le.qm2ds%Mx) qm2ds%istore_M=qm2ds%Mx
      if (qm2ds%istore.le.qm2ds%istore_M) qm2ds%istore=qm2ds%istore_M

      do j=1,qm2ds%Mx
         call mo2site(qm2ds,qm2ds%v0(1,j),qm2ds%v2(1,j),qm2ds%rrwork)
      end do
   end if

   return     
   end

!
!********************************************************************
!
!********************************************************************
!
   subroutine davidson0(qm2_params,qmmm_nml,qmmm_mpi,cosmo_c_struct,qm2_struct,qm2ds,qmmm_struct,M4,lprint,ftol0,ftol1,ferr, &
      Np,Nh,j0,j1,e0,v0,kflag,iee2,ee2,eta,xi, &
      nd,nd1,vexp1,vexp,ray,rayv,rayvL,rayvR,raye,raye1, &
      ray1,ray1a,ray2,idav,istore)
   use qm2_davidson_module   
   use cosmo_C,only:cosmo_C_structure
   use qmmm_struct_module, only : qmmm_struct_type
   use qmmm_module,only:qm2_structure, qmmm_mpi_structure
   use qm2_params_module,  only : qm2_params_type
    use qmmm_nml_module   , only : qmmm_nml_type

   implicit none
     type(qmmm_nml_type),intent(inout) :: qmmm_nml
     type(cosmo_C_structure), intent (inout) :: cosmo_c_struct
     type(qmmm_struct_type), intent(inout) :: qmmm_struct
     type(qm2_davidson_structure_type), intent(inout) :: qm2ds
     type(qm2_structure),intent(inout) :: qm2_struct
     type(qmmm_mpi_structure),intent(inout) :: qmmm_mpi
     type(qm2_params_type),intent(inout) :: qm2_params

   !logical check_symmetry; !!JAB Testing

   integer Np,Nh,M4,lprint,kflag,nd,nd1,nd1_old,j0,j1
   integer one,i,j,k,n,m,icount,idav,istore,iloop,itarget
   integer info,iee2(qm2ds%Ncis)
   !integer l,u,c
   _REAL_ ftol0,ftol1,ferr(j1+j0),ddot,fn,f2m
   _REAL_ f0,f1,f2,f3,f4,tresh2
   _REAL_ e0(M4),v0(qm2ds%Nrpa,qm2ds%Mx),ee2(qm2ds%Ncis)
   _REAL_ eta(qm2ds%Nrpa),xi(qm2ds%Nrpa)
!  Davidson expansion vectors:
   _REAL_ vexp1(qm2ds%Ncis,nd),vexp(qm2ds%Nrpa,nd)

! Davidson Rayleigh matrices:
   _REAL_ ray1(nd,nd),ray2(nd,nd),ray1a(nd,nd)
   _REAL_ ray(nd,nd),raye(nd),raye1(nd)
   _REAL_ rayv(nd,nd),rayvL(nd,nd),rayvR(nd,nd)
   _REAL_, allocatable :: dtmp(:)
   !_REAL_ f11;

   if (lprint.gt.4) write(6,*)' Entering davidson0'

   call clearing(nd*nd,ray1(1,1))
   call clearing(nd*nd,ray1a(1,1))
   call clearing(nd*nd,ray2(1,1))
   call clearing(nd*nd,ray(1,1))
   call clearing(nd*nd,rayv(1,1))
   call clearing(nd*nd,rayvL(1,1))
   call clearing(nd*nd,rayvR(1,1))
   call clearing(nd,raye)
   call clearing(nd,raye1)

   one=1

   icount=0
   nd1_old=0
   m=0
   n=0
   iloop=0

   if (istore.gt.0) then    ! Use vectors from the previous step
      j1=min(j1,istore)
      goto 70 
   end if 

! **** Assign Davidson trial vectors

85 continue

   tresh2=qm2ds%tresh*0.9
   do j=1,nd   
      call clearing (2*M4,vexp(1,j))
      call clearing (M4,vexp1(1,j))
   end do 

! vexp and vexp1 are now zero !!JAB
! assign trial vectors based on the MO
   itarget =0
   do j=1,nd1
80    continue

      itarget=itarget+1
      if (itarget.ge.Np*Nh) goto 85 ! Restart vectors
      f1=0.0

      do i=1,j0
         f1=f1+v0(iee2(itarget),i)**2+v0(iee2(itarget)+M4,i)**2
      end do

      if (f1.ge.tresh2) goto 80  ! MO pair is not accepted
      if (lprint.gt.2) write(6,*)'++ START ',itarget,iee2(itarget)

      vexp1(iee2(itarget),j) = vexp1(iee2(itarget),j)+1.0
   end do

!  Orthogonolize trial vectors (vexp1) to found eigenvectors, some strange
!  function XI(V0) !!JAB
 
  !!ORTHONORMALIZE WITH XI=X+Y??
   
   do i=1,j0
      !!MAKE XI
      do k=1,M4
         xi(k)=v0(k,i)+v0(k+M4,i) !xi=X+Y??
      end do
      !!NORMALIZE
      f2=ddot(M4,xi,one,xi,one) !xi.xi
      f2=1.0D0/sqrt(abs(f2)) !normalization constant

      call dscal(M4,f2,xi,one) !normalize xi
     
      do j=1,nd1
         !!ORTHOGONALIZE
         f1=-ddot(M4,xi,one,vexp1(1,j),one) !xi.vexp1
         call daxpy(M4,f1,xi,one,vexp1(1,j),one) !vexp1=vexp1+xi.vexp1*xi
      end do

      f2=ddot(M4,vexp(1,j),one,vexp1(1,j),one) !vexp.vexp1 but vexp=0 if started from 85??
      f2=1.0D0/sqrt(abs(f2)) !normalization constant
      call dscal(M4,f2,vexp1(1,j),one) !normalize

      !!SAME AS ABOVE BUT WITH XI=X-Y
      do k=1,M4
         xi(k)=v0(k,i)-v0(k+M4,i) !xi=X-Y??
      end do

      f2=ddot(M4,xi,one,xi,one) !xi.xi
      f2=1.0D0/sqrt(abs(f2)) 
      call dscal(M4,f2,xi,one) !normalize xi

      do j=1,nd1   
         f1=-ddot(M4,xi,one,vexp1(1,j),one) 
         call daxpy(M4,f1,xi,one,vexp1(1,j),one)
      end do

      f2=ddot(M4,vexp(1,j),one,vexp1(1,j),one)
      f2=1.0D0/sqrt(abs(f2))
      call dscal(M4,f2,vexp1(1,j),one)
   end do

!  ORTHOGONALIZE AND NORMALIZE TRIAL VECTORS
90     continue
   do j=1,nd1

      !!NORMALIZE
      f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(M4,f2,vexp1(1,j),one)

      do i=1,j-1
         
         !!ORTHOGONALIZE
         f1=ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
         if(abs(f1).gt.qm2ds%tresh1) then !!REMOVE VERY NONORTHOGONAL VECTORS
            call dcopy(M4,vexp1(1,nd1),one,vexp1(1,j),one) 
            nd1=nd1-1 !!NUMBER OF VECTORS
            if(lprint.gt.2) write(6,*) 'Trial removed #',j,i,f1,nd1
            goto 90
         end if
         
         call daxpy(M4,-f1,vexp1(1,i),one,vexp1(1,j),one) !!ORTHOGONALIZE
      end do

      f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one) !!NORMALIZE AGAIN
      f2=1/sqrt(abs(f2))
      call dscal(M4,f2,vexp1(1,j),one)
   end do 

!  CHECK ORTHOGONALITY OF EXPANSIONS
   if(lprint.gt.3) write(6,*) 'Check expansion to expansion'
   do j=1,nd1
      do i=1,nd1
         f1= ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
         if(f1.gt.1.0E-15.and.i.ne.j) goto 90
      end do
   end do

   if (nd1.le.j1) j1=nd1
! **** End assign Davidson trial vectors
10 continue
! **** Write some things, check the number of iterations, etc.
   icount=icount+1 !iteration counter
   if(lprint.gt.4) write(6,*) 'COUNT=',icount,'Exp=',nd1
   if(icount.gt.qm2ds%icount_M) then
         write(6,*) "Number of davidson iterations exceeded, exiting"
         stop
   endif
   if(lprint.gt.4) write(6,*) 'nd1,nd1_old,j0,j1',nd1,nd1_old,j0,j1

! **** Davidson Restart
   if(nd1.gt.nd) then !if the number of Krylov subspace expansions is  
      iloop=iloop+1
      if(lprint.gt.4) write(6,*)
      write(6,*) 'Davidson not converged, expansion=',nd
   end if

70     continue

   if(nd1.gt.nd.or.istore.gt.0) then  ! Use vectors from the previous step
      istore=0
      if (lprint.gt.4) write(6,*) 'Restart Davidson with guesses from the previous step'
      if(lprint.gt.4) write(6,*) 'Restart loop = ',iloop
      icount=0
      nd1_old=0
      m=0
      n=0

      if(idav.eq.2) nd1=min((nd-2),2*j1)  ! RPA
      if(idav.eq.1) nd1=min((nd-2),j1)  ! CIS
      if(lprint.gt.1) write(6,*) 'Currently have ',j1,' states'
      if(lprint.gt.1) write(6,*) 'With ',nd1,' initial guesses'

      !NEW EXPANSION VECTORS
      if(idav.eq.2) then ! RPA 
         do j=1,j1
          do i=1,M4
             if(j.le.nd1) vexp1(i,j)=v0(i,j+j0)+v0(i+M4,j+j0)
             if((j+j1).le.nd1) vexp1(i,j+j1)=v0(i,j+j0)-v0(i+M4,j+j0)
          end do
         end do
      end if

      if(idav.eq.1) then ! CIS 
         do j=1,j1
          do i=1,M4
             if(j.le.nd1) vexp1(i,j)=v0(i,j+j0) 
          end do
         end do
      end if


!  Orthogonalize and normalize new expansion vectors
12    continue
      do j=1,nd1
         f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
         f2=1/sqrt(abs(f2))
         call dscal(M4,f2,vexp1(1,j),one)

         do i=1,j-1
            f1=ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
            if (abs(f1).gt.qm2ds%tresh1) then
               call dcopy(M4,vexp1(1,nd1),one,vexp1(1,j),one) 
               nd1=nd1-1
               if(lprint.gt.2) write(6,*) 'New expansion removed #',j,i,f1,nd1
               goto 12
            end if

            call daxpy(M4,-f1,vexp1(1,i),one,vexp1(1,j),one)
         end do

         f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
         f2=1/sqrt(abs(f2))
         call dscal(M4,f2,vexp1(1,j),one)
      end do


!  Check
      do j=1,nd1
       do i=1,nd1
          f1=ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
          if (f1.gt.1.0E-15.and.i.ne.j.and.lprint.gt.0) write(6,*) i,j,f1
          if(f1.gt.1.0E-15.and.i.ne.j) goto 12
       end do
      end do

      goto 10
   end if

! **** End Davidson Restart	
! Next is calculated vexp=L(vexp1), so vexp becomes the Lxi of the guess xi !!JAB
! Vexp1 was the orthonormalized guess expansion vectors from the previous
! section. Throughout comments using A and B are related to  L(xi)=L([X;Y]) similarly 
! to [A,B;B,A][X;Y] from RPA eigenvalue problem, but definately not the band 
! ABBA as one might mistakenly assume !!JAB

   do i=nd1_old+1,nd1
      call clearing (2*M4,eta)
      call dcopy(M4,vexp1(1,i),one,eta,one)
      call Lxi_testing(qm2_params,qmmm_nml,qmmm_mpi,cosmo_c_struct,qm2_struct, &
                       qm2ds,qmmm_struct,eta,vexp(1,i),cosmo_c_struct%solvent_model)
       
! CIS - set Y=0
      if(idav.eq.1) then
         do j=1,M4
            vexp(M4+j,i)=0.0
         enddo
      endif

! form vexp(1,i) having Ab and Bb to (A+B)b and (A-B)b	  
! thus vexp was L([X;Y], now it is [LX+LY;LX-LY] !!JAB
      do j=1,M4
         f1=vexp(j,i)
         vexp(j,i)=vexp(j,i)+vexp(M4+j,i)
         vexp(M4+j,i)=f1-vexp(M4+j,i)
      end do
   end do

!       stop
! **** Operations in Krylov space
! form ray1=b(A-B)b and ray2=b(A+B)b	
!That is, these are formed from vexp=L(xi) and the new guess expansions Vexp1
   do i=1,nd1
      do j=nd1_old+1,nd1
         ray1(i,j)=ddot(M4,vexp1(1,i),one,vexp(1,j),one)
         ray2(i,j)=ddot(M4,vexp1(1,i),one,vexp(M4+1,j),one)
      end do
   end do
  
   if(nd1_old.ne.0) then
      do i=nd1_old+1,nd1
       do j=1,nd1
          ray1(i,j)=ddot(M4,vexp1(1,i),one,vexp(1,j),one)
          ray2(i,j)=ddot(M4,vexp1(1,i),one,vexp(M4+1,j),one)
       end do
      end do
   end if

   nd1_old=nd1


! form ray1a=sqrt(b(A-B)b)
   f1=1.0
   f0=0.0

   call dcopy(nd*nd,ray1,one,rayvR,one)


   write(6,*)'info00',info 

   allocate(dtmp(4*nd1))
   call dsyev ('v','u',nd1,rayvR,nd,raye,dtmp,4*nd1,info) 
        !Eigenvalues of ray1 in raye and eigenvectors in rayvR
   write(6,*)'info000',info
   write(6,*)'info Nrpa,nd,nd1',qm2ds%Nrpa,nd,nd1
   write(6,*)'info shapes',shape(rayvR),shape(raye),shape(xi)
 
   do j=1,nd1
      raye1(j)=Sign(Sqrt(Abs(raye(j))),raye(j)) !make raye1 sqrt(eigenvalues) with same sign as eigenvalues
   end do

   do j=1,nd1
    do i=1,nd1
       rayv(i,j)=rayvR(j,i)*raye1(i)
    end do
   end do
 
   call dgemm('N','N',nd1,nd1,nd1,f1,rayvR,nd,rayv,nd,f0,ray1a,nd) !


   call dgemm('N','N',nd1,nd1,nd1,f1,ray2,nd,ray1a,nd,f0,rayv,nd)
   call dgemm('N','N',nd1,nd1,nd1,f1,ray1a,nd,rayv,nd,f0,ray,nd)
   call dcopy(nd*nd,ray,one,rayv,one)
   call symmetr(nd,rayv)


! find eigenvalues and eigenvectors of ray
   write(6,*)'info0',info
   call dsyev ('v','u',nd1,rayv,nd,raye,xi,qm2ds%Nrpa,info)
   do j=1,nd1
      raye(j) = Sign(Sqrt(Abs(raye(j))),raye(j))
   enddo
   deallocate(dtmp)

! Solve for Right EigenVector  rayvR = |X+Y>=ray1a*rayv
   call dgemm('N','N',nd1,nd1,nd1,f1,ray1a,nd,rayv,nd,f0,rayvR,nd)
! Solve for Left EigenVector  
! rayvL = |X-Y> = 1/E b(A+B)b|X+Y> =1/raye * ray2 * rayvR
   call dgemm('N','N',nd1,nd1,nd1,f1,ray2,nd,rayvR,nd,f0,rayvL,nd)

   do j=1,nd1
    do i=1,nd1
       rayvL(i,j)=rayvL(i,j)/raye(j)
    end do
   end do

   if(lprint.gt.2) then 
      write(6,*)'info',info
      write(6,910) (raye(i),i=1,j1)
   end if 

   if(raye(1).le.0.1) goto 100           ! RPA bad behavior
! **** End operations in Krylov space

!  Form approximate eigenvectors
   do i=j0+1,j0+j1
      call clearing(2*M4,v0(1,i))
   end do

   do k=j0+1,j0+j1

 !new eigenvectors from converged onward (j0 is the number of converged vectors)
 !vexp1 is calculated at the beginning of the routine. The previous Krylov space
 !calculations result in f1 and f2 below which are used to mix vexp1 with
 !unconverged vectors in v0
      do i=1,nd1
         f1=rayvR(i,(k-j0))
         f2=rayvL(i,(k-j0))
         call daxpy(M4,f1,vexp1(1,i),one,v0(1,k),one) !V0=V0+f1*Vexp1
         call daxpy(M4,f2,vexp1(1,i),one,v0(1+M4,k),one)
      end do

      do i=1,M4 !Whatever is going on here, it looks like some X-Y, X+Y stuff
         f3=v0(i,k)
         v0(i,k)=-v0(i,k)-v0(M4+i,k)
         v0(M4+i,k)=v0(M4+i,k)-f3
      end do
       
! CIS : Y=0
      if(idav.eq.1) then
         do i=1,M4
            v0(i+M4,k)=0.0
         end do
      end if
   end do

!  Orthogonalize and normalize approximate eigenvectors
   do j=j0+1,j0+j1
      do i=j0+1,j-1
         f1=-ddot(M4,v0(1,i),one,v0(1,j),one) &
            +ddot(M4,v0(1+M4,i),one,v0(1+M4,j),one)
         call daxpy(2*M4,f1,v0(1,i),one,v0(1,j),one)
      end do

      f2=ddot(M4,v0(1,j),one,v0(1,j),one) &
         -ddot(M4,v0(1+M4,j),one,v0(1+M4,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(2*M4,f2,v0(1,j),one)

   end do
 
! FIND eigenvalue, residual vectors and residual norm:
! print results of this iteration and check for converged vectors

   if (lprint.gt.1) write(6,*)'eigenvalues and residual norm'
   n=0
   m=0

   do j=j0+1,j0+j1
      
      call Lxi_testing(qm2_params,qmmm_nml,qmmm_mpi,cosmo_c_struct, qm2_struct,qm2ds,qmmm_struct, &
                       v0(1,j),eta,cosmo_c_struct%solvent_model) !L(xi) output in eta !!JAB
      f1=ddot(M4,v0(1,j),one,eta(1),one) &
         -ddot(M4,v0(1+M4,j),one,eta(1+M4),one)
! CIS
      if (idav.eq.1) f1=ddot(M4,v0(1,j),one,eta(1),one) !E=<xi|L(xi)> now in f1 !!JAB
      call dcopy(2*M4,eta,one,xi,one)
      call daxpy(2*M4,-f1,v0(1,j),one,xi,one)
      f2=ddot(M4,xi(1),one,xi(1),one)
      f3=ddot(M4,xi(1+M4),one,xi(1+M4),one)
! CIS
      if(idav.eq.1) then !Clear Y after Liouville equation
         f3=0.0 
         call clearing(M4,xi(1+M4))
         call clearing(M4,eta(1+M4))
      endif

      f2=f2+f3 !f(x) + f(y)
      if(f2.le.ftol0) then ! Converged vector
         n=n+1
         call dcopy(2*M4,v0(1,j0+n),one,v0(1,j),one) !Move converged vector to beginning
         e0(j0+n)= f1 !move converged energy to beginning
         ferr(j0+n)=abs(f2)+abs(f3) !

         if(lprint.gt.3) then
            write(6,"(3i5,1x,2f14.9,2g10.3,2x,A)") &
               j,n,m,raye(j-j0),f1,f2,f3,' Converged!'
         end if
! CIS
      else if(idav.eq.1) then       

 ! FORM perturbed residual vectors for new initial guesses
         if((nd1+m).eq.nd) goto 45
         m=m+1
        
         !Not obvious how this is done, but f1/(f1-ee2(i) is some sort of energy
         !weighting and eta is L(xi) v0
         do i=1,M4
            vexp1(i,nd1+m)=(eta(i)-f1*v0(i,j))/(f1-ee2(i))
         end do

         f4=ddot(M4,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
         f4=1/sqrt(abs(f4))
         call dscal(M4,f4,vexp1(1,nd1+m),one)

         if(lprint.gt.3) then 
                write(6,"(3i5,1x,2f14.9,2g10.3)") &
               j,n,m,raye(j-j0),f1,f2,f3
         end if
      else
         if((nd1+m).eq.nd) goto 45
         m=m+1

         do i=1,M4
            vexp1(i,nd1+m)=(eta(i)-eta(i+M4)-f1*(v0(i,j)-v0(i+M4,j))) &
               /(f1-ee2(i))
         end do

         f4=ddot(M4,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
         f4=1/sqrt(abs(f4))
         call dscal(M4,f4,vexp1(1,nd1+m),one)

         if((nd1+m).eq.nd) goto 45
         m=m+1

         do i=1,M4
            vexp1(i,nd1+m)=(eta(i)+eta(i+M4)-f1*(v0(i,j)+v0(i+M4,j))) &
               /(f1-ee2(i))
         end do

         f4=ddot(M4,vexp1(1,nd1+m),one,vexp1(1,nd1+m),one)
         f4=1/sqrt(abs(f4))
         call dscal(M4,f4,vexp1(1,nd1+m),one)
         if(lprint.gt.1) write(6,"(3i5,1x,2f14.9,2g10.3)") &
            j,n,m,raye(j-j0),f1,f2,f3
      end if          
   end do
!          
45 continue
        
   if((j1-n).eq.0) then
      if(lprint.gt.4) write(6,*) 'All vectors found after loop' &
         ,iloop, ', Expansion ', nd1
      goto 100
   end if

   if(m.eq.0) then ! Restart Davidson
      nd1=nd+1
      goto 10
   end if

   if(lprint.gt.4) write(6,*) 'New perturbed m=',m
  
15 continue

   do j=1+nd1,nd1+m
      f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(M4,f2,vexp1(1,j),one)

      do i=1,j-1
         f1=ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
         if(abs(f1).gt.qm2ds%tresh) then
            call dcopy(M4,vexp1(1,nd1+m),one,vexp1(1,j),one) 
            m=m-1     
            if(lprint.gt.2) write(6,*) 'exp removed #',j,i,f1,nd1+m
            goto 15
         end if

         call daxpy(M4,-f1,vexp1(1,i),one,vexp1(1,j),one)
      end do

      f2=ddot(M4,vexp1(1,j),one,vexp1(1,j),one)
      f2=1/sqrt(abs(f2))
      call dscal(M4,f2,vexp1(1,j),one)
   end do
!  Check
   do j=1+nd1,nd1+m
    do i=1,j-1
       f1=ddot(M4,vexp1(1,i),one,vexp1(1,j),one)
       if(f1.gt.1.0E-15.and.i.ne.j) goto 15          
    end do
   end do


   if(lprint.gt.1) write(6,*) 'Perturbed left m=',m

   if(m.eq.0) then      
      write(6,*) 'Run out of expansion vectors'
      goto 100 
   end if

   if(iloop.gt.qm2ds%iloop_M) then
      write(6,*)'Davidson not converged after ',iloop,' loops'
      goto 100 
   end if

   if((nd1+m).gt.(nd-1).and.nd1.ne.nd) then
      nd1=nd
   else
      nd1=nd1+m
   endif

   goto 10

100   continue
   j0=j0+n
   if (n.eq.0) then
      write(6,*)'Could not go further ',j0,' vectors'
      kflag=1
   end if
    
   if(lprint.gt.0) then
      write(6,*) '@@@@ Davidson subroutine Found vectors',j0
   end if

910   format(' ', 100g11.3)

   return
   end
!
