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

program calccoup
    implicit none

!VARIABLES FOR INPUT AND READING FILES
    real, allocatable, dimension(:, :) :: xyz1, charge1, xyz2, charge2
    real, allocatable, dimension(:, :):: displace, rotate
    integer, allocatable, dimension(:) :: nstates, ncharges
    integer :: nmolec
    character(LEN=2), allocatable, dimension(:) :: poletype1, poletype2
    character(LEN=5) :: fname1, junk
    character(LEN=1) :: fname2
    character(LEN=3) :: fname3
    character(LEN=9) :: filename, filename2
    real :: garbage

!VARIABLES FOR CALCULATING THE COULOMB POTENTIAL
    integer :: n1, n2, l1, l2
    real :: displace1(3), displace2(3), R, elem_chg, perm, J_to_eV, A_to_Bohr
    real, allocatable, dimension(:, :) :: V, VM, VD, VLQ, VSQ

!VARIABLES TO COMMUNICATE WITH SUBROUTINE TO READ INPUT FILES
    interface
        subroutine readfile(filename, nstates, ncharges, xyz, charge, poletype)
            real, dimension(:, :), intent(out) :: xyz, charge
            integer, intent(in) :: nstates, ncharges
            character(2), dimension(:), intent(out) :: poletype
            character(9), intent(in) :: filename
        end subroutine readfile
    end interface

    nmolec = 2

    fname1 = 'molec'
    fname3 = '.in'
    filename2 = 'couple.in'

    Elem_chg = 1.60217657E-19
    perm = 8.854187817E-12
    A_to_Bohr = 1.8897
    J_to_eV = 6.24E18

    allocate (nstates(nmolec))
    allocate (ncharges(nmolec))

!Open instructions file couple.in
    open (71, file=filename2, status='old')
    read (71, '(I5,I5)') nstates(1), nstates(2)
    read (71, '(F5.1,1X,F5.1,1X,F5.1)') displace1(1), displace1(2), displace1(3)
    read (71, '(F5.1,1X,F5.1,1X,F5.1)') displace2(1), displace2(2), displace2(3)
    close (71)

!Writing instructions that were read
    write (6, *) 'Calculating V for', nstates(1), 'and', nstates(2), 'states.'
    write (6, *) 'Displacement for molecule 1:', displace1
    write (6, *) 'Displacement for molecule 2:', displace2

!Allocating space for couplings
    allocate (V(nstates(1), nstates(2)))
    allocate (VM(nstates(1), nstates(2)))
    allocate (VD(nstates(1), nstates(2)))
    allocate (VLQ(nstates(1), nstates(2)))
    allocate (VSQ(nstates(1), nstates(2)))

!Read files to find the number of charges
1000 format(F10.7, 2X, F10.7, 2X, F10.7, 2X, F10.7, 2X, A2)!molec file
    open (71, file='molec1.in')
    do n1 = 1, 5
        read (71, "(A)", ERR=5, end=20) junk
    end do
5   continue
    ncharges(1) = 1
    do while (1 == 1)
        read (71, 1000, ERR=10, end=20) garbage
        ncharges(1) = ncharges(1) + 1
    end do
10  continue
    close (71)
    open (71, file='molec2.in')
    do n1 = 1, 5
        read (71, "(A)", ERR=15, end=20) junk
    end do
15  continue
    ncharges(2) = 1
    do while (1 == 1)
        read (71, 1000, ERR=20, end=20) garbage
        ncharges(2) = ncharges(2) + 1
    end do
20  continue
    close (71)

    write (6, *) 'Ncharges', ncharges

!Allocating space for variables read from files
    allocate (xyz1(3, ncharges(1)))
    allocate (charge1(ncharges(1), nstates(1)))
    allocate (poletype1(ncharges(1)))
    allocate (xyz2(3, ncharges(2)))
    allocate (charge2(ncharges(2), nstates(2)))
    allocate (poletype2(ncharges(2))) EAD FILES

    nmolec = 1
    write (fname2, '(I1)') nmolec !molecule number for filename
    filename = fname1//fname2//fname3 !piece together filename
    write (6, *) 'Opening file ', filename
    call readfile(filename, nstates(nmolec), ncharges(nmolec), xyz1, charge1, poletype1)
    nmolec = 2
    write (fname2, '(I1)') nmolec !molecule number for filename
    filename = fname1//fname2//fname3 !piece together filename
    write (6, *) 'Opening file ', filename
    call readfile(filename, nstates(nmolec), ncharges(nmolec), xyz2, charge2, poletype2)

!Calculate Coulomb Coupling
!displace coordinates
    do n1 = 1, 3
        xyz1(n1, :) = xyz1(n1, :) + displace1(n1)
        xyz2(n1, :) = xyz2(n1, :) + displace2(n1)
    end do

    write (6, *) 'Summed Charges 1:', sum(charge1, DIM=1)
    write (6, *) 'Summed Charges 2:', sum(charge2, DIM=1)
    write (6, *) '*******************************'
    write (6, *) '***   Electronic Coupling   ***'
    write (6, *) '*******************************'
    write (6, *) 'Mom V(a.u)    V(eV)     St1 St2'

    do l1 = 1, nstates(1)
        do l2 = 1, nstates(2)
            do n1 = 1, ncharges(1)
                do n2 = 1, ncharges(2)

                    R = sqrt((xyz1(1, n1) - xyz2(1, n2))**2 + (xyz1(2, n1) - xyz2(2, n2))**2 + (xyz1(3, n1) - xyz2(3, n2))**2)*A_to_Bohr
                    if (R < 1) then
                        write (6, *) '**Charges', n1, n2, 'are close at', R, 'Bohr'
                        stop
                    end if

                    V(l1, l2) = V(l1, l2) + charge1(n1, l1)*charge2(n2, l2)/R
                    if (poletype1(n1) .EQ. 'M' .AND. poletype2(n2) .EQ. 'M') then
                        VM(l1, l2) = VM(l1, l2) + charge1(n1, l1)*charge2(n2, l2)/R
                    end if
                    if (poletype1(n1) .EQ. 'D' .AND. poletype2(n2) .EQ. 'D') then
                        VD(l1, l2) = VD(l1, l2) + charge1(n1, l1)*charge2(n2, l2)/R
                    end if
                    if (poletype1(n1) .EQ. 'LQ' .AND. poletype2(n2) .EQ. 'LQ') then
                        VLQ(l1, l2) = VLQ(l1, l2) + charge1(n1, l1)*charge2(n2, l2)/R
                    end if
                    if (poletype1(n1) .EQ. 'SQ' .AND. poletype2(n2) .EQ. 'SQ') then
                        VSQ(l1, l2) = VSQ(l1, l2) + charge1(n1, l1)*charge2(n2, l2)/R
                    end if

                end do
            end do

            write (6, "(A2,2X,F8.5,2X,F8.5,2X,I2,2X,I2)") 'T', V(l1, l2), V(l1, l2)*27.2114, l1, l2
            write (6, "(A2,2X,F8.5,2X,F8.5,2X,I2,2X,I2)") 'M', VM(l1, l2), VM(l1, l2)*27.2114, l1, l2
            write (6, "(A2,2X,F8.5,2X,F8.5,2X,I2,2X,I2)") 'D', VD(l1, l2), VD(l1, l2)*27.2114, l1, l2
            write (6, "(A2,2X,F8.5,2X,F8.5,2X,I2,2X,I2)") 'LQ', VLQ(l1, l2), VLQ(l1, l2)*27.2114, l1, l2
            write (6, "(A2,2X,F8.5,2X,F8.5,2X,I2,2X,I2)") 'SQ', VSQ(l1, l2), VSQ(l1, l2)*27.2114, l1, l2

        end do
    end do

    stop

end program calccoup

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

subroutine readfile(filename, nstates, ncharges, xyz, charge, poletype)
    implicit none

    real, dimension(:, :), intent(out) :: xyz, charge
    integer, intent(in) :: nstates, ncharges
    integer ::  n, m, atomnumber
    character(2), dimension(:), intent(out) :: poletype
    character(9), intent(in) :: filename
    character(100) :: header

1000 format(F10.7, 2X, F10.7, 2X, F10.7, 2X, F10.7, 2X, A2)!molec file

    open (71, file=filename, status='old') !Open molec file

    m = 0 !set number of states to zero
    n = ncharges

100 continue

    m = m + 1

    if (n .NE. ncharges) then
        write (6, *) 'Wrong number of charges specified'
        return
    end if
    n = 0

    do while (1 == 1)

        read (71, '(A)', end=200) header !Check for the header to start reading
        if (index(header, 'Atom') .NE. 0) then !Check for start to read
            write (6, *) 'Reading transition ', m
            do while (1 == 1) !n=1,ncharges !capture data for normal mode nstates
                n = n + 1
                if (n .GT. ncharges) goto 100
                !the read statement gives an error at the end of the normal mode
                !because the formatting is incorrect, so return to the top of the loop
                read (71, 1000, ERR=100, end=200) xyz(1, n), xyz(2, n), &
                    xyz(3, n), charge(n, m), &
                    poletype(n)
            end do

        end if

    end do

200 continue
    write (6, *) 'File succefully read'
    close (71)

    return
end subroutine readfile
