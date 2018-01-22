#include "dprec.fh"                                    
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! This program calculates the transition dipole moments
! between excited states using CEO formalism                                         
!   
!    MODIFIED FOR SQM_NAESMD by JAKB LANL 2015
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!   
    
subroutine polarizab(qm2_params,qmmm_nml,qm2_struct,qm2ds,qmmm_struct)                         
      use qmmm_module,only:qm2_structure
      use qm2_davidson_module
      use constants, only : BOHRS_TO_A, SQRT2, ONE_AU
!TODO compare units between CEO and SQM_NAESMD
      use qmmm_struct_module, only : qmmm_struct_type
      use qm2_params_module,  only : qm2_params_type
      use qmmm_nml_module   , only : qmmm_nml_type

    
    
      implicit none                            
        type(qmmm_nml_type),intent(inout) :: qmmm_nml
        type(qm2_params_type),intent(inout) :: qm2_params
        type(qm2_structure),intent(inout) :: qm2_struct
        type(qm2_davidson_structure_type), intent(inout) :: qm2ds
        type(qmmm_struct_type), intent(inout) :: qmmm_struct
        _REAL_ f0,f1,f2,f3,ddot,freq
        integer zero,one,two,three
        !integer nfreq_M                       

        parameter (f0 = 0.0)                   
        parameter (f1 = 1.0)                   
        parameter (f2 = -2.0)                  
        parameter (zero=0)                     
        parameter (one=1)                      
        parameter (two=2)                      
        parameter (three=3)                    
        !parameter (nfreq_M=10000)
    
        integer modf(5),modpr(qm2ds%Mx),Mpp,Mhh
        integer nfreq1,nfreq2,nfreq3           
        _REAL_ freq01,freq02,freq03            
        _REAL_ freq11,freq12,freq13
        _REAL_ rdamp1,rdamp2,rdamp3
        _REAL_ omega12,omega13,omega23         
        character*3 type2,type3                
        _REAL_ tr,f                            
        _REAL_ muax(qm2ds%Mx),muay(qm2ds%Mx),muaz(qm2ds%Mx)
        _REAL_ muabx(-qm2ds%Mx:qm2ds%Mx,-qm2ds%Mx:qm2ds%Mx)
        _REAL_ muaby(-qm2ds%Mx:qm2ds%Mx,-qm2ds%Mx:qm2ds%Mx)
        _REAL_ muabz(-qm2ds%Mx:qm2ds%Mx,-qm2ds%Mx:qm2ds%Mx)
    
    
        _REAL_ tmp1(qm2ds%Nrpa),tmp2(qm2ds%Nrpa)
        _REAL_ temp1(qm2ds%Nb**2),temp2(qm2ds%Nb**2),temp3(qm2ds%Nb**2),temp4(qm2ds%Nb**2)
        _REAL_ temp5(qm2ds%Nb**2),temp6(qm2ds%Nb**2)

      integer itime1,itime2,itime3,itime11,get_time
      integer i,j                              
      real time11                              
      character*20 datetime
      character*(150) txt, txt1*15, machname*36    
      _REAL_ dip(3,qm2ds%Lt) !added            
    
      call get_date(datetime)
      call get_machine(machname)
      itime1=get_time()                        


      write (6,*)
      write (6,*) '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
      write (6,8) '| POL execution started at ',         datetime,' |'
      write (6,7) '| Computer: ',                         machname,'|'
      write (6,*) '|________________________________________________|'
      write (6,*)

! read input parameters and Cartesian coordinates modes frequencies      

210   format(A)

      Mpp = qm2ds%Np*qm2ds%Np
      Mhh = qm2ds%Nh*qm2ds%Nh


      print *, '!!!!!!-----INPUT-----!!!!!!'
      print *
      print *, 'INDO basis set size          ', qm2ds%Nb
      print *, 'Number of electrons (pairs)  ', qm2ds%Np
      print *, 'Number of holes (pairs)      ', qm2ds%Nh
      print *, 'Modes to be sorted           ', qm2ds%Mx
      print *, 'Total number of atoms        ', qmmm_struct%nquant_nlink
      print *, 'Calculated number of states  ', qm2ds%Mx
      print *
      do i=1,qmmm_struct%nquant_nlink
         write(6,100) qmmm_struct%iqm_atomic_numbers(i),qmmm_struct%qm_coords(:,i)
      enddo
      print *
      print *, '!!!!!!-----End INPUT-----!!!!!!'

      print *, 'WARNING, the units output from this POL calculation are in &
                unknown units and are only for relative magnitude (FIXME)'
100   format(I5,'     ',3F12.6,2i4)
110   format(I5,'     ',F12.6,i4)
120   format(' Modes factors Av,S,Px,Py,Pz     ',5i4)

if(1==1) then

! Note after mo2siteph the mode eta is ksi(-alpha)* sqrt(2) different 
! from usual normalization condition, i.e Tr(rho[eta,eta^T])=2
!     Compute mu_alpha (transition dipoles)
        call get_dipole_matrix(qm2_params,qmmm_nml,qm2_struct,qmmm_struct,qmmm_struct%qm_coords, dip)
        call unpacking(qm2ds%Nb,dip(1,:),temp1,'s')
        do j=1,qm2ds%Mx
           muax(j)=ddot(qm2ds%Nb**2,temp1,one,qm2ds%v2(:,j),one)*sqrt(2.0)
        enddo
        call unpacking(qm2ds%Nb,dip(2,:),temp1,'s')
        do j=1,qm2ds%Mx
           muay(j)=ddot(qm2ds%Nb**2,temp1,one,qm2ds%v2(:,j),one)*sqrt(2.0)
        enddo
        call unpacking(qm2ds%Nb,dip(3,:),temp1,'s')
        do j=1,qm2ds%Mx
           muaz(j)=ddot(qm2ds%Nb**2,temp1,one,qm2ds%v2(:,j),one)*sqrt(2.0)
        enddo
!  Print mu_alpha
      print *
      print *, 'mu_alpha'
      print *, 'i e(i) fx(i) fy(i) fz(i) f(i)'
      do i=1,qm2ds%Mx
         f=sqrt(muax(i)**2+muay(i)**2+muaz(i)**2)
         write(6,130) i, qm2ds%e0(i), 2*qm2ds%e0(i)*muax(i)**2/21.0, 2*qm2ds%e0(i)*muay(i)**2/21.0, &
        2*qm2ds%e0(i)*muaz(i)**2/21.0, 2*qm2ds%e0(i)*f**2/21.0
      enddo
130   format(I5,5F12.4)

!     Zero all nonlinear couplings matrices
      do i=-qm2ds%Mx,qm2ds%Mx
         do j=-qm2ds%Mx,qm2ds%Mx
            muabx(i,j)=0.0
            muaby(i,j)=0.0
            muabz(i,j)=0.0
         enddo
      enddo
      call unpacking(qm2ds%Nb,dip(1,:),temp1,'s')
      call unpacking(qm2ds%Nb,dip(1,:),temp3,'s')  !Necessary?
      call unpacking(qm2ds%Nb,qm2_struct%den_matrix,temp2,'s')
      call dgemm ('N','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,f2,temp1,qm2ds%Nb,temp2,qm2ds%Nb,f1,temp3,qm2ds%Nb)
      do j=-qm2ds%Mx,qm2ds%Mx
         if (j.ne.0) then
            temp1=qm2ds%v2(:,abs(j))
            if (j.gt.0) call transp1(qm2ds%Nb,temp1)
            call dgemm ('N','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,f1,temp3,qm2ds%Nb,temp1,qm2ds%Nb,f0,temp4,qm2ds%Nb)
            do i=-qm2ds%Mx,qm2ds%Mx
               if (i.ne.0) then
                  temp2=qm2ds%v2(:,abs(i))
                  if (i.lt.0) call transp1(qm2ds%Nb,temp2)
                  tr=ddot(qm2ds%Nb**2,temp4,one,temp2,one)
                  muabx(j,i)=muabx(j,i)+tr
                  muabx(i,j)=muabx(i,j)+tr
               endif
            enddo
         endif
       enddo

       call unpacking(qm2ds%Nb,dip(2,:),temp1,'s')
       call unpacking(qm2ds%Nb,dip(2,:),temp3,'s')
       call unpacking(qm2ds%Nb,qm2_struct%den_matrix,temp2,'s')
       call dgemm ('N','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,f2,temp1,qm2ds%Nb,temp2,qm2ds%Nb,f1,temp3,qm2ds%Nb)
       do j=-qm2ds%Mx,qm2ds%Mx
          if (j.ne.0) then
             temp1=qm2ds%v2(:,abs(j))
             if (j.gt.0) call transp1(qm2ds%Nb,temp1)
             call dgemm ('N','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,f1,temp3,qm2ds%Nb,temp1,qm2ds%Nb,f0,temp4,qm2ds%Nb)
             do i=-qm2ds%Mx,qm2ds%Mx
                if (i.ne.0) then
                   temp1=qm2ds%v2(:,abs(i))
                   if (i.lt.0) call transp1(qm2ds%Nb,temp2)
                   tr=ddot(qm2ds%Nb**2,temp4,one,temp2,one)
                   muaby(j,i)=muaby(j,i)+tr
                   muaby(i,j)=muaby(i,j)+tr
                endif
             enddo
          endif
       enddo
      call unpacking(qm2ds%Nb,dip(3,:),temp1,'s')
      call unpacking(qm2ds%Nb,dip(3,:),temp3,'s')
      call unpacking(qm2ds%Nb,qm2_struct%den_matrix,temp2,'s')
      call dgemm ('N','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,f2,temp1,qm2ds%Nb,temp2,qm2ds%Nb,f1,temp3,qm2ds%Nb)
      do j=-qm2ds%Mx,qm2ds%Mx
         if (j.ne.0) then
            temp1=qm2ds%v2(:,abs(j))
            if (j.gt.0) call transp1(qm2ds%Nb,temp1)
            call dgemm ('N','N',qm2ds%Nb,qm2ds%Nb,qm2ds%Nb,f1,temp3,qm2ds%Nb,temp1,qm2ds%Nb,f0,temp4,qm2ds%Nb)
            do i=-qm2ds%Mx,qm2ds%Mx
               if (i.ne.0) then
                  temp2=qm2ds%v2(:,abs(i))
                  if (i.lt.0) call transp1(qm2ds%Nb,temp2)
                  tr=ddot(qm2ds%Nb**2,temp4,one,temp2,one)
                  muabz(j,i)=muabz(j,i)+tr
                  muabz(i,j)=muabz(i,j)+tr
               endif
             enddo
          endif
       enddo

       itime2=get_time()
       time11=real(itime2-itime1)/100
       print *
       print *, 'Computed mu_alpha_beta, time',time11, 'sec'
       open (qm2ds%muab_unit,file=trim(qm2ds%muab_out))
      print *
      print *, 'mu_alpha_beta'
      print *, 'j i e(ji) fabx(ji) faby(ji) fabz(ji) fab(ji)'
      do j=-1,-qm2ds%Mx,-1
         do i=abs(j),qm2ds%Mx
            f=sqrt(muabx(j,i)**2+muaby(j,i)**2+muabz(j,i)**2)
            write(6,140) j,i,qm2ds%e0(i)-qm2ds%e0(-j),2*(qm2ds%e0(i)-qm2ds%e0(-j))*muabx(j,i)**2/21.0, &
            2*(qm2ds%e0(i)-qm2ds%e0(-j))*muaby(j,i)**2/21.0, &
            2*(qm2ds%e0(i)-qm2ds%e0(-j))*muabz(j,i)**2/21.0,2*(qm2ds%e0(i)-qm2ds%e0(-j))*f**2/21.0
            write(qm2ds%muab_unit,140) abs(j),i,qm2ds%e0(i)-qm2ds%e0(-j),2*(qm2ds%e0(i)-qm2ds%e0(-j))*muabx(j,i)**2/21.0, &
            2*(qm2ds%e0(i)-qm2ds%e0(-j))*muaby(j,i)**2/21.0, &
            2*(qm2ds%e0(i)-qm2ds%e0(-j))*muabz(j,i)**2/21.0,2*(qm2ds%e0(i)-qm2ds%e0(-j))*f**2/21.0
         enddo
      enddo
      close(qm2ds%muab_unit)

140   format(2I5,5F12.4)

      open (qm2ds%ceo_unit,file=trim(qm2ds%ceo_out))
      j=-qmmm_struct%state_of_interest
      write(qm2ds%ceo_unit,*)
      write(qm2ds%ceo_unit,*) 'Energies (eV) and oscillator strengths for transitions'
      write(qm2ds%ceo_unit,*) 'from state ',abs(j),' to all other states'
      do i=1,qm2ds%Mx
         f=sqrt(muabx(j,i)**2+muaby(j,i)**2+muabz(j,i)**2)
         write(qm2ds%ceo_unit,130) i,qm2ds%e0(i)-qm2ds%e0(-j),2*(qm2ds%e0(i)-qm2ds%e0(-j))*f**2/21.0
      enddo
      flush(qm2ds%ceo_unit)
      close(qm2ds%ceo_unit)

      itime11=get_time()
      time11=real(itime11-itime1)/100
      itime11=time11
      itime1=MOD(itime11,60)
      itime11=itime11/60
      itime2=MOD(itime11,60)
      itime11=itime11/60
      itime3=MOD(itime11,24)
      itime11=itime11/24
      call get_date(datetime)

endif
      write (6,*)
      write (6,*) '|^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|'
      write (6,8) '| POL normal termination at',         datetime,' |'
      write (6,9) '| POL total CPU time    ',   time11,   ' seconds |'
      write(6,301) itime11,itime3,itime2,itime1
      write (6,*) '|________________________________________________|'
      write (6,*)

 7    format (' ',A12,' ',A35,A2)
 8    format (' ',A27,'  ',A19,A2)
 9    format (' ',A20,'    ',g16.5,A10)
301   format(' |     ',i2,' days ',i2,' hours ',i2,' minutes ',i2,' seconds     |')
      call flush(6)

      stop

29    print *, txt
      stop 'Bad input file'

end subroutine

