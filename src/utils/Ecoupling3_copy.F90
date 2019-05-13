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
      real, allocatable, dimension(:,:):: displace, rotate
      integer, allocatable, dimension(:) :: nstates, ncharges, mu
      integer :: nmolec,n,m
      character(LEN=2), allocatable,dimension(:) :: poletype1,poletype2
      character(LEN=5) :: fname1,junk
      character(LEN=1) :: fname2
      character(LEN=3) :: fname3
      character(LEN=9) :: filename,filename2
      character(LEN=50) :: header
      real :: garbage

!VARIABLES FOR CALCULATING THE COULOMB POTENTIAL
      integer :: n1,n2,l1,l2
      real :: displace1(3),displace2(3),R,test1,test2,elem_chg,perm,J_to_eV,A_to_Bohr
      real,allocatable,dimension(:,:) :: V     
      
!VARIABLES TO COMMUNICATE WITH SUBROUTINE TO READ INPUT FILES
INTERFACE
 subroutine readfile(filename,nstates,ncharges,xyz,charge,poletype)
      real, allocatable, dimension(:,:) :: xyz , charge
      integer :: nstates, ncharges, n, m, atomnumber
      character(2), allocatable, dimension(:) :: poletype
      character(9) :: filename
      character(4) :: header
 end subroutine readfile
end INTERFACE

nmolec=2

allocate(nstates(nmolec))
allocate(ncharges(nmolec))
allocate(mu(nmolec))

fname1='molec'
fname3='.in'
filename2='couple.in'

!Open instructions file couple.in
open(71,file=filename2,status='old')
read(71,'(I5,I5)') nstates(1),nstates(2)
read(71,'(F5.1,1X,F5.1,1X,F5.1)') displace1(1),displace1(2),displace1(3)
read(71,'(F5.1,1X,F5.1,1X,F5.1)') displace2(1),displace2(2),displace2(3)
read(71,'(F5,.1,1X,F5.1,1X,F5.1)') mu(1),mu(2)
close(71)

write(6,*)'Calculating V for',nstates(1),'and',nstates(2),'states.'
write(6,*)'Displacement for molecule 1:',displace1
write(6,*)'Displacement for molecule 2:',displace2

1000 format(F10.7,2X,F10.7,2X,F10.7,2X,F10.7,2X,A2)!molec file
!Read files to find the number of charges
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
CLOSE(71)
OPEN(71,file='molec2.in')
do n1=1,5
READ(71,"(A)",ERR=15,END=20) junk
enddo
15 CONTINUE
ncharges(2)=1
do while(1==1)
READ(71,1000,ERR=20,END=20) garbage
ncharges(2)=ncharges(2)+1
enddo
20 CONTINUE
CLOSE(71)



elem_chg=1.60217657E-19
perm=8.854187817E-12
A_to_Bohr=1.8897
J_to_eV=6.24E18


allocate(V(nstates(1),nstates(2)))

!!READ FILES
nmolec=1
    allocate(xyz1(3,ncharges(nmolec)))
    allocate(charge1(ncharges(nmolec),nstates(nmolec)))
    allocate(poletype1(ncharges(nmolec)))
    write(fname2,'(I1)') nmolec !molecule number for filename
    filename=fname1//fname2//fname3 !piece together filename
    write(6,*) 'Opening file ',filename
nmolec=2
    allocate(xyz2(3,ncharges(nmolec)))
    allocate(charge2(ncharges(nmolec),nstates(nmolec)))
    allocate(poletype2(ncharges(nmolec)))
    write(fname2,'(I1)') nmolec !molecule number for filename
    filename=fname1//fname2//fname3 !piece together filename
    write(6,*) 'Opening file ',filename

!Calculate Coulomb Coupling
!displace coordinates
do n1=1,3
xyz1(n1,:)=xyz1(n1,:)+displace1(n1)
xyz2(n1,:)=xyz2(n1,:)+displace2(n1)
enddo
write(6,*) 'Summed Charges 1:',sum(charge1,DIM=1)
write(6,*) 'Summed Charges 2:',sum(charge2,DIM=1)
write(6,*) '***************************'
write(6,*) '*** Electronic Coupling ***'
write(6,*) '***************************'
write(6,*) 'V(a.u)    V(eV)     St1 St2'

do l1=1,nstates(1)
do l2=1,nstates(2)

test1=0
test2=0

do n1=1,ncharges(1)
  do n2=1,ncharges(2)
    
    R=sqrt((xyz1(1,n1)-xyz2(1,n2))**2+(xyz1(2,n1)-xyz2(2,n2))**2+(xyz1(3,n1)-xyz2(3,n2))**2)*A_to_Bohr
    if (R<1) then
        write(6,*)'**Charges',n1,n2,'are close at',R,'Bohr'
        STOP
    end if

    V(l1,l2)=V(l1,l2)+charge1(n1,l1)*charge2(n2,l2)/R
   
    test2=test2+charge2(n2,l2)
    test1=test1+charge1(n1,l1)
  end do
end do

write(6,"(F8.5,2X,F8.5,2X,I2,2X,I2)") V(l1,l2),V(l1,l2)*27.2114,l1,l2
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

subroutine readfile(filename,nstates,ncharges,xyz,charge,poletype)
implicit none

      real, allocatable, dimension(:,:) :: xyz , charge
      integer :: nstates, ncharges, n, m, atomnumber
      character(2), allocatable, dimension(:) :: poletype
      character(9) :: filename
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
       enddo

    endif

  enddo

  200 Continue
  write(6,*) 'File succefully read'
  close(71)

RETURN
end subroutine
