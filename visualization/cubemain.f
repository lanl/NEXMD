        Program cubemain

        Implicit none
      Integer numatoms,numorbs
      OPEN(8, FILE='EVEC.DATA')
      READ(8,*)
      READ(8,*)numatoms
      READ(8,*)
      READ(8,*)numorbs
      rewind 8
      Close (8)
      Call readcube(numorbs,numatoms)
	  End

        Subroutine readcube(numorbs,numatoms)

        Implicit none
      Integer numatoms,numorbs
      Double Precision xyz(numatoms,3), readcoeff(numorbs)
      Integer i,atype(numatoms),ireaderror,IERR,
     &   print_count,ivec,numoccorbs,itally,inummd
      real etime          ! Declare the type of etime()
      real elapsed(2)     ! For receiving user and system time
      real total          ! For receiving total time



      itally=0
      ireaderror=1
      OPEN(8, FILE='EVEC.DATA')
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)
      READ(8,*)numoccorbs
      READ(8,*)
      READ(8,*)print_count
C          Ben, looking at current variable sizes
C      WRITE(*,*) "Inside readcube"
C      WRITE(*,*) "numatoms: ",numatoms
C      WRITE(*,*) "numorbs: ",numorbs
C      WRITE(*,*) "numoccorbs: ",numoccorbs
C      WRITE(*,*) "print_count: ",print_count
C      WRITE(*,*) "Done"
C           Do inummd=1,10
      Do While (ireaderror.eq.1)
C      This read statement reads eigenvector line	  
C       WRITE(*,*) "Reading Atomic coordinates line"
       READ(8,*,IOSTAT=IERR,ERR=100)
C check to see if we are at the end of the file
       If (IERR.NE.0) GOTO 100
       Do i=1,numatoms
C       READ(8,*,IOSTAT=IERR,ERR=100)atype(i),xyz(i,1),xyz(i,2),xyz(i,3)
C       This reads the geometry
C        WRITE(*,*) "Reading Geometry line"
        READ(8,*)atype(i),xyz(i,1),xyz(i,2),xyz(i,3)
       EndDo
       Do ivec=1,print_count
C       This reads the eigenvector header
C        WRITE(*,*) "Reading Eigenvector headewr"
        READ(8,*)
        Do i=1,numorbs
C        This reads the eigenvector coeffienct
C         Write(*,*) "Reading Coeffiencient"
         READ(8,*)readcoeff(i)
        EndDo
C only o the homo for now
        If (ivec.eq.numoccorbs) Then
         itally=itally+1
         write(*,*)itally
         total = etime(elapsed)
         print *, 'End: total=', total, ' user=', elapsed(1),
     &        ' system=', elapsed(2)
         call buildcube(numorbs,numatoms,atype,xyz,readcoeff,itally)
        EndIf
C end if numoccorbs=ivec
       EndDo
C end loop over ivec
       READ(8,*)
C read trailing blank line
       EndDo
C end loop over inummd
       Close (8)
100    return
	  End

        double precision Function sto(itype,xi,r)
         Double Precision xi,r,Pi
         Integer itype
           pi=dacos(0.0d0)*2.D0
           If (itype.eq.1) Then
             sto=DSQRT(xi**3/PI)*DEXP(-xi*r)
           EndIf
           If (itype.eq.2) Then
             sto=DSQRT(xi**5/(3.0d0*PI))*DEXP(-xi*r)
           EndIf
           If (itype.eq.3) Then
             sto=DSQRT(xi**5/PI)*DEXP(-xi*r)
           EndIf
        Return
        End


      subroutine buildcube(numorbs,numatoms,atype,xyz,readcoeff,itally)
        Implicit none
	Integer numatoms,numorbs
	Double Precision xyz(numatoms,3),cubedim(3)
        Double Precision sto,voxel,CM(3),coeff1,coeff2
        Double Precision local_x,local_y,local_z,local_r,
     &                   current_x,current_y,current_z,
     &                   mocoeff(numatoms,5),orbcoeff(numatoms,2),
     &                   readcoeff(numorbs),cubit(6),padding,
     &                   cmax,cmin,cormax(3)
	Integer i,j,atype(numatoms),orbcount,
     &   ix,iy,iz,iatom,numvox(3),jj,itally
        Character filename*14
	Character tmpstr*4

        voxel=0.25D0
        padding=10.0D0
		
C       Ben: Printing values to verify I'm not working with zeros
C        write(*,*) "Ben: printing Coefficients"
C        do i=1,numorbs
C         write(*,"(f12.4)") readcoeff(i)
C		enddo
C		write(*,*) "Ben: writing geometry"
C		do i=1,numatoms
C		 write(*,"(f12.8)",advance='no') ( xyz(i,j), j=1,3 )
C		 write(*,*) ""
C		enddo

C read the number of orbitals
C read the coeff for the MO of interest
C assumed that Hydrogen has one 1s orbital 
C and all others have zero 1s and one 2s and three 2p

C To simplify things later on we are going to repack the Coeff
C So there are 5 AO coeff for each atom
         orbcount=1
         Do i=1,numatoms
           If (atype(i).eq.1) Then
              mocoeff(i,1)=readcoeff(orbcount)
              orbcount=orbcount+1
              Do j=2,5
                mocoeff(i,j)=0.0D0
              EndDo
           Else
              mocoeff(i,1)=0.0D0
              Do j=2,5
                mocoeff(i,j)=readcoeff(orbcount)
                orbcount=orbcount+1
              EndDo
           EndIf
         EndDo
		 
C		write(*,*) "Ben: printing MO coefficients"
C		do i=1,numatoms
C         write(*,"(f12.4)",advance='no') ( mocoeff(i,j), j=1,5 )
C		 write(*,*) ""
C		enddo

C set the orbital type and paramters for each atom
C reference Szabo and Ostlund p. 186
        Do i=1,numatoms
          If (atype(i).eq.1) Then
            orbcoeff(i,1)=1.24D0
            orbcoeff(i,2)=0.0D0
          ElseIf (atype(i).eq.3) Then
            orbcoeff(i,1)=2.69D0
            orbcoeff(i,2)=0.75D0
          ElseIf (atype(i).eq.4) Then
            orbcoeff(i,1)=3.68D0
            orbcoeff(i,2)=1.1D0
          ElseIf (atype(i).eq.5) Then
            orbcoeff(i,1)=4.68D0
            orbcoeff(i,2)=1.45D0
          ElseIf (atype(i).eq.6) Then
            orbcoeff(i,1)=5.67D0
            orbcoeff(i,2)=1.72D0
          ElseIf (atype(i).eq.7) Then
            orbcoeff(i,1)=6.67D0
            orbcoeff(i,2)=1.95D0
          ElseIf (atype(i).eq.8) Then
            orbcoeff(i,1)=7.66D0
            orbcoeff(i,2)=2.25D0
          ElseIf (atype(i).eq.9) Then
            orbcoeff(i,1)=8.65D0
            orbcoeff(i,2)=2.55D0
          Else
C I am approximating these coefficients
            orbcoeff(i,1)=atype(i)*1.0D0
            orbcoeff(i,2)=2.8D0
          EndIf
        EndDo

C scale coordinate system
           Do i=1,3
           Do j=1,numatoms
             xyz(j,i)=xyz(j,i)/0.529177249D0
           EndDo
           EndDo

           Do i=1,3
           cubedim(i)=0.0D0
             cmin=xyz(1,i)
             cmax=cmin
           Do j=2,numatoms
             If (xyz(j,i).lt.cmin) cmin=xyz(j,i)
             If (xyz(j,i).gt.cmax) cmax=xyz(j,i)
           EndDo
             cormax(i)=cmax
             CM(i)=cmax-cmin
C Add a small region so that the molecule is well inside the cube 
             If (CM(i).gt.cubedim(i)) cubedim(i) = CM(i) + padding
C determine the number of voxel elements
           numvox(i) = Int(cubedim(i)/voxel)
C make numvox even or don't
           numvox(i) = numvox(i) + Mod(numvox(i),2)
           EndDo
            write(tmpstr,'(I4.4)')itally
            filename='EVEC.'//tmpstr//'.cube'
          OPEN(7, FILE=filename)
C first two lines are title
         write(7,*)'cube file'
         write(7,*)
C third line is number of atoms and origin
         write(7,114) numatoms,
     &  cormax(1)-(CM(1)+cubedim(1))/2.0D0,
     &  cormax(2)-(CM(2)+cubedim(2))/2.0D0,
     &  cormax(3)-(CM(3)+cubedim(3))/2.0D0
C lines 4,5,6   give the number of voxels and size of each voxel?
C allegedly, if the number of voxels is negative Angstroms are used
        write(7,113)numvox(1),voxel,0.0D0,0.0D0
        write(7,113)numvox(2),0.0D0,voxel,0.0D0
        write(7,113)numvox(3),0.0D0,0.0D0,voxel
        Do i=1,numatoms
          write(7,112)atype(i),xyz(i,1),xyz(i,2),xyz(i,3)
        EndDo
C to call 1s slater:  sto(1,xi,r)
C to call 2s slater:  r*sto(2,xi,r)
C to call 2p[x,y,z] slater:  [x,y,z]*sto(3,xi,r)
  
C initallize current position
       current_x=cormax(1)-(CM(1)+cubedim(1))/2.0D0
       Do ix=1,numvox(1)
       current_y=cormax(2)-(CM(2)+cubedim(2))/2.0D0
         Do iy=1,numvox(2)
           current_z=cormax(3)-(CM(3)+cubedim(3))/2.0D0
           Do iz=1,numvox(3),6
             Do jj=1,6
               cubit(jj)=0.0D0
               Do iatom=1,numatoms
                 local_x=current_x-xyz(iatom,1)
                 local_y=current_y-xyz(iatom,2)
                 local_z=current_z-xyz(iatom,3)
                 local_r=DSQRT((local_x**2+local_y**2+local_z**2))
                 coeff1=orbcoeff(iatom,1)
                 coeff2=orbcoeff(iatom,2)
                 cubit(jj)=cubit(jj)
     &                +mocoeff(iatom,1)*sto(1,coeff1,local_r)
     &                +mocoeff(iatom,2)*local_r*sto(2,coeff2,local_r)
     &                +mocoeff(iatom,3)*local_x*sto(3,coeff2,local_r)
     &                +mocoeff(iatom,4)*local_y*sto(3,coeff2,local_r)
     &                +mocoeff(iatom,5)*local_z*sto(3,coeff2,local_r)
                  EndDo
                  current_z=current_z+voxel
                EndDo
                If (iz.le.(numvox(3)-6)) write(7,111)(cubit(jj),jj=1,6)
              EndDo
              write(7,111)(cubit(jj),jj=1,numvox(3)-(iz-7))
              current_y=current_y+voxel
            EndDo
            current_x=current_x+voxel
          EndDo
          Close (7)
        Return
111     Format(6(E13.5))
112     Format('   ',I2,'    0.000000  ',F10.6,'  ',F10.6,'  ',F10.6)
113     Format(2x,I3,4x,F8.6,4x,F8.6,4x,F8.6)
114     Format('  ',I6,'  ',F12.8,'  ',F12.8,'  ',F12.8)
        End


