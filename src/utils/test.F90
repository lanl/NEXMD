
PROGRAM calccoup
implicit none
!VARIABLES FOR INPUT AND READING FILES
      real, allocatable, dimension(:,:) :: xyz1, xyz2,charge1,charge2
     ! real, allocatable, dimension(:,:):: displace, rotate
      integer, allocatable, dimension(:) :: nstates, ncharges
      !real, allocatable, dimension(:) :: mu
      integer :: nmolec, mcount(4)
      character(LEN=2), allocatable,dimension(:) ::poletype1,poletype2
      character(LEN=5) :: fname1,junk
      character(LEN=1) :: fname2
      character(LEN=3) :: fname3
      character(LEN=9) :: filename,filename2
      real :: garbage

!VARIABLES FOR CALCULATING THE COULOMB POTENTIAL
      integer :: n1,n2,l1,l2
      real,dimension(3) ::displace1(3),displace2(3)
      real :: R,elem_chg,perm,J_to_eV,A_to_Bohr
      real,allocatable,dimension(:,:) :: V,VM,VD,VLQ,VSQ

nmolec=2

allocate(nstates(nmolec))
allocate(ncharges(nmolec))
nstates(1)=5
nstates(2)=5
!Open instructions file couple.in
open(10,file='couple.in',status='old')
read(10,'(I5,I5)') nstates(1),nstates(2)
read(10,'(F8.1,F8.1,F8.1)') displace1(1),displace1(2),displace1(3)
read(10,'(F8.1,F8.1,F8.1)') displace2(1),displace2(2),displace2(3)
!read(71,'(F5.1,1X,F5.1,1X,F5.1)') mu(1),mu(2)



!allocate variables to calculate
allocate(VLQ(nstates(1),nstates(2)))
allocate(V(nstates(1),nstates(2)))
allocate(VM(nstates(1),nstates(2)))
allocate(VD(nstates(1),nstates(2)))
allocate(VSQ(nstates(1),nstates(2)))
write(6,*) size(V,1),size(V,2)

close(10)
write (6,*) V,VM,VD,VSQ,VLQ
write (6,*) displace1,displace2
end program
