!!********************************************************************
!
!  Program to calculate electronic coupling from NDDO style charges calculated
!  with subroutine printCfitNM() using input files coupling.in molecule1.in and
!  molecule2.in and outputs coupling.out
!
!********************************************************************
!  1.) Read coupling.in for number of states 1 and 2, displacement of 1 and 2,
!  and rotation of 1 and 2.
!  2.) Read molecule1 and 2.in for coordinates and charges
!--------------------------------------------------------------------
!
!  Josiah A. Bjorgaard NDSU 2013
!
!--------------------------------------------------------------------

PROGRAM calccoup
implicit none

!VARIABLES FOR INPUT AND READING FILES
      real, allocatable, dimension(:,:) :: xyz1, charge1, xyz2, charge2
      integer, allocatable, dimension(:) :: nstates, ncharges
      integer :: nmolec
      character(LEN=2), allocatable,dimension(:) :: poletype1,poletype2
      character(LEN=5) :: fname1,junk
      character(LEN=1) :: fname2
      character(LEN=3) :: fname3
      character(LEN=9) :: filename,filename2
      real :: garbage

!VARIABLES FOR CALCULATING THE COULOMB POTENTIAL
      integer :: n1,n2,l1,l2
      real :: displace1(3),displace2(3),R,elem_chg,perm,J_to_eV,A_to_Bohr
      real,allocatable,dimension(:,:) :: V,VM,VD,VLQ,VSQ     
      
!VARIABLES TO COMMUNICATE WITH SUBROUTINE TO READ INPUT FILES
INTERFACE
subroutine readfile(nstates,ncharges,xyz,charge,poletype,filename)
implicit none
      real, intent(out),dimension(:,:) :: xyz , charge
      integer, intent(in) :: nstates, ncharges
      integer ::  n, m, atomnumber
      character(2), dimension(:), intent(out) :: poletype
      character(9), intent(in) :: filename
      character(200) :: header
end subroutine readfile
end INTERFACE

nmolec=2
fname1='molec'
fname3='.in'
filename2='couple.in'

Elem_chg=1.60217657E-19
perm=8.854187817E-12
A_to_Bohr=1.8897
J_to_eV=6.24E18

allocate(nstates(nmolec))
allocate(ncharges(nmolec))

!Open instructions file couple.in
open(71,file=filename2,status='old')
read(71,'(I5,I5)') nstates(1),nstates(2)
read(71,'(F5.1,1X,F5.1,1X,F5.1)') displace1(1),displace1(2),displace1(3)
read(71,'(F5.1,1X,F5.1,1X,F5.1)') displace2(1),displace2(2),displace2(3)
close(71)

!Writing instructions that were read
write(6,*)'Calculating V for',nstates(1),'and',nstates(2),'states.'
write(6,*)'Displacement for molecule 1:',displace1
write(6,*)'Displacement for molecule 2:',displace2

!Read files to find the number of charges
1000 format(F10.7,2X,F10.7,2X,F10.7,2X,F10.7,2X,A2)!molec file
OPEN(71,file='molec1.in')
do n1=1,5
READ(71,"(A)",ERR=5,END=20) junk
enddo
5 CONTINUE
ncharges(1)=1
do while(1==1)
READ(71,1000,ERR=10,END=20) garbage
ncharges(1)=ncharges(1)+1
enddo
10 CONTINUE
20 CONTINUE
CLOSE(71)


write(6,*) 'Ncharges',ncharges

allocate(xyz1(3,ncharges(1)))
allocate(xyz2(3,ncharges(2)))
allocate(charge1(ncharges(1),nstates(1)))
allocate(charge2(ncharges(2),nstates(2)))
allocate(poletype2(ncharges(2)))
allocate(poletype1(ncharges(1)))
allocate(VM(nstates(1),nstates(2)))
allocate(VD(nstates(1),nstates(2)))
allocate(VLQ(nstates(1),nstates(2)))
allocate(VSQ(nstates(1),nstates(2)))

!!READ FILES
nmolec=1
    write(6,*) 'Opening file ',filename
    call readfile(nstates(1),ncharges(1),xyz1,charge1,poletype1,'molec1.in')
write(6,*) size(charge1)
nmolec=2
    write(fname2,'(I1)') nmolec !molecule number for filename
    call readfile(nstates(2),ncharges(2),xyz2,charge2,poletype2,'molec2.in')
write(6,*) size(charge2)

!Calculate Coulomb Coupling
!displace coordinates
do n1=1,3
xyz1(n1,:)=xyz1(n1,:)+displace1(n1)
xyz2(n1,:)=xyz2(n1,:)+displace2(n1)
enddo

write(6,*) 'Summed Charges 1:',sum(charge1,DIM=1)
write(6,*) 'Summed Charges 2:',sum(charge2,DIM=1)
write(6,*) '*******************************'
write(6,*) '***   Electronic Coupling   ***'
write(6,*) '*******************************'
write(6,*) 'Mom V(a.u)    V(eV)     St1 St2'

do l1=1,nstates(1)
do l2=1,nstates(2)
do n1=1,ncharges(1)
do n2=1,ncharges(2)
    
    R=sqrt((xyz1(1,n1)-xyz2(1,n2))**2+(xyz1(2,n1)-xyz2(2,n2))**2+(xyz1(3,n1)-xyz2(3,n2))**2)*A_to_Bohr
    if (R<1) then
        write(6,*)'**Charges',n1,n2,'are close at',R,'Bohr'
        STOP
    end if

    V(l1,l2)=V(l1,l2)+charge1(n1,l1)*charge2(n2,l2)/R
    if (poletype1(n1).EQ.'M'.AND.poletype2(n2).EQ.'M') then
    VM(l1,l2)=VM(l1,l2)+charge1(n1,l1)*charge2(n2,l2)/R
    endif
    if (poletype1(n1).EQ.'D'.AND.poletype2(n2).EQ.'D') then
    VD(l1,l2)=VD(l1,l2)+charge1(n1,l1)*charge2(n2,l2)/R
    endif
    if (poletype1(n1).EQ.'LQ'.AND.poletype2(n2).EQ.'LQ') then
    VLQ(l1,l2)=VLQ(l1,l2)+charge1(n1,l1)*charge2(n2,l2)/R
    endif
    if (poletype1(n1).EQ.'SQ'.AND.poletype2(n2).EQ.'SQ') then
    VSQ(l1,l2)=VSQ(l1,l2)+charge1(n1,l1)*charge2(n2,l2)/R
    endif
   
end do
end do

write(6,"(A2,2X,F8.5,2X,F8.5,2X,I2,2X,I2)") 'T',V(l1,l2),V(l1,l2)*27.2114,l1,l2
write(6,"(A2,2X,F8.5,2X,F8.5,2X,I2,2X,I2)") 'M',VM(l1,l2),VM(l1,l2)*27.2114,l1,l2
write(6,"(A2,2X,F8.5,2X,F8.5,2X,I2,2X,I2)") 'D',VD(l1,l2),VD(l1,l2)*27.2114,l1,l2
write(6,"(A2,2X,F8.5,2X,F8.5,2X,I2,2X,I2)") 'LQ',VLQ(l1,l2),VLQ(l1,l2)*27.2114,l1,l2
write(6,"(A2,2X,F8.5,2X,F8.5,2X,I2,2X,I2)") 'SQ',VSQ(l1,l2),VSQ(l1,l2)*27.2114,l1,l2

enddo
enddo

STOP

end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  SUBROUTINE TO READ FILES CONTAINING CHARGES FROM NDDO TYPE CALCULATION WITH
!  NORMAL MODES
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Josiah A. Bjorgaard, 2013, North Dakota State University
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine readfile(nstates,ncharges,xyz,charge,poletype,filename)
implicit none
      real, dimension(:,:), intent(inout) :: xyz , charge
      integer, intent(in) :: nstates, ncharges
      integer ::  n, m, atomnumber
      character(2), dimension(:), intent(inout) :: poletype
      character(9), intent(in) :: filename
      character(100) :: header

      1000 format(F10.7,2X,F10.7,2X,F10.7,2X,F10.7,2X,A2)!molec file
      
      open(71,file=filename,status='old') !Open molec file
    
      m=0 !set number of states to zero
      n=ncharges

      100 Continue

      m=m+1

      if (n.NE.ncharges) then
        write(6,*) 'Wrong number of charges specified'
        RETURN
      end if
      n=0

  do while (1==1)
      
      READ(71,'(A)',END=200) header !Check for the header to start reading
    write(6,*) 'reading header',header
    if (index(header,'Atom') .NE. 0 ) then !Check for start to read
    write(6,*) 'Reading transition ',m
      do while (1==1) !n=1,ncharges !capture data for normal mode nstates
          n=n+1
          if (n.GT.ncharges) GOTO 100
          !the read statement gives an error at the end of the normal mode
          !because the formatting is incorrect, so return to the top of the loop
          READ(71,1000,ERR=100,END=200) xyz(1,n),xyz(2,n),&
                                xyz(3,n),charge(n,m),&
                                poletype(n)
write(6,*) xyz(1,n)
       enddo

    endif

  enddo

  200 Continue
  write(6,*) 'File succefully read'
  write(6,*) charge
  close(71)

RETURN
end subroutine
