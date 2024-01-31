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

    real, dimension(1:, 1:) :: xyz, charge
    real, dimension(1:3, 1:):: displace, rotate
    integer, dimension(1:) :: nstates, ncharges
    integer :: nstates1, ncharges1, nmolec, n, m, l, atomnumber
    character(2), dimension(1:) :: poletype
    character :: fname1, fname2, fname3, filename, header

    nmolec = 2

!!READ FILES
    do n = 1, nmolec !loop over number of molecules

        write (fname2, '(I1)') nmolec !molecule number for filename
        filename = fname1//fname2//fname3 !piece together filename
        open (71, file=filename, status='old') !Open molec file
        nstates1 = 0 !set number of states to zero

        do while (1 == 1)

100         continue !Restart point when the end of a normal mode is reached
            ncharges1 = 0 !set number of charges to zero
            nstates1 = nstates1 + 1

            read (71, '(A)', end=200) header !Check for the header to start reading a section

            if (index(header, 'Atom') .NE. 0) then !Check for start to read

                do while (1 == 1) !capture data for normal mode nstates

                    ncharges1 = ncharges1 + 1
                    !the read statement gives an error at the end of the normal m
                    !because the formatting is incorrect, so return to the top of the loop
                    read (71, 1000, ERR=100) atomnumber, xyz(1, ncharges, n), xyz(2, ncharges, n), &
                        xyz(3, ncharges, n), charge(ncharges, nstates, n), &
                        poletype(ncharges, n)
                end do

            end if

        end do

200     continue
        close (71)
        ncharges(n) = ncharges1
        nstates(n) = nstates1
    end do

!!WRITE FILES
    open (71, file='output.test')

    do n = 1, nmolec
        do m = 1, nstates(n)
            do l = 1, ncharges(n)

                write (71, *) xyz(1, m, n), xyz(2, m, n), xyz(3, m, n), charge(l, m, n), poletype(m, n)

            end do
        end do
    end do

    close (71)

    stop

end program calccoup

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine readfile(filename, nstates, ncharges, xyz, charge, poletype)
    implicit none

    real, allocatable(:, :) :: xyz, charge
    integer :: nstates, ncharges, n, m, atomnumber
    character(2), allocatable, dimension(:) :: poletype
    character :: filename, header

1000 format(I4, 3X, F10.7, 3X, F10.7, 3X, F10.7, 3X, F10.7, A10, 5X, A5, A5)!molec file

    allocate (xyz(3, ncharges))
    allocate (charge(ncharges, nstates))
    allocate (poletype(ncharges))

    open (71, file=filename, status='old') !Open molec file

    m = 0 !set number of states to zero
    n = ncharges

    do while (1 == 1)

100     continue !Restart point when the end of a normal mode is reached
        if (n .NE. ncharges) print(6, *) 'Error, number of charges is incorrect'

        n = 0 !set number of charges to zero
        m = m + 1

        read (71, '(A)', end=200) header !Check for the header to start reading

        if (index(header, 'Atom') .NE. 0) then !Check for start to read

            do while (1 == 1) !capture data for normal mode nstates

                m = m + 1
                !the read statement gives an error at the end of the normal mode
                !because the formatting is incorrect, so return to the top of the loop
                read (71, 1000, ERR=100) atomnumber, xyz(1, n), xyz(2, n), &
                    xyz(3, n), charge(n, m), &
                    poletype(n)
            end do

        end if

    end do

200 continue
    close (71)

    return
end subroutine readfile
