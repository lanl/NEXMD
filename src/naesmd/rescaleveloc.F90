#include "dprec.fh"
!*****************************
!Subroutine to eliminate rotation and translation in dynamics by rescaling
!velocities
!
!****************************
subroutine rescaleveloc(rx,ry,rz,vx,vy,vz,massmdqt,natom)
    _REAL_ :: rxcm,rycm,rzcm
        _REAL_ :: rx(natom),ry(natom),rz(natom),massmdqt(natom)
        _REAL_ :: vxcm,vycm,vzcm,masstot,kincm
        _REAL_ :: vx(natom),vy(natom),vz(natom)
        _REAL_ :: mang(3),vang(3)
        integer :: natom,i,j
        _REAL_ :: xxi,xyi,xzi,yyi,yzi,zzi
        _REAL_ :: xdel,ydel,zdel
        _REAL_ :: tensor(3,3)
        _REAL_ :: erot
        rxcm=0.0d0
        rycm=0.0d0
        rzcm=0.0d0
        masstot=0.0d0
        do j=1,natom
            masstot=masstot+massmdqt(j)
        enddo
        do j=1,natom
            rxcm=rxcm+rx(j)*massmdqt(j)/masstot
            rycm=rycm+ry(j)*massmdqt(j)/masstot
            rzcm=rzcm+rz(j)*massmdqt(j)/masstot
        enddo
        ! calculate the velocity of the center of mass
        vxcm=0.0d0
        vycm=0.0d0
        vzcm=0.0d0
        do j=1,natom
            vxcm=vxcm+vx(j)*massmdqt(j)/masstot
            vycm=vycm+vy(j)*massmdqt(j)/masstot
            vzcm=vzcm+vz(j)*massmdqt(j)/masstot
        enddo
        ! eliminate any rotation about the system center of mass
        ! calculate translational kinetic energy of overall system
        kincm=0.5d0*masstot*(vxcm**2+vycm**2+vzcm**2)
        ! calculate the angular momentum of the overall system
        do j = 1, 3
            mang(j) = 0.0d0
        end do
        do j=1,natom
            mang(1)=mang(1)+(ry(j)*vz(j)-rz(j)*vy(j))*massmdqt(j)
            mang(2)=mang(2)+(rz(j)*vx(j)-rx(j)*vz(j))*massmdqt(j)
            mang(3)=mang(3)+(rx(j)*vy(j)-ry(j)*vx(j))*massmdqt(j)
        enddo
        mang(1)=mang(1)-(rycm*vzcm-rzcm*vycm)*masstot
        mang(2)=mang(2)-(rzcm*vxcm-rxcm*vzcm)*masstot
        mang(3)=mang(3)-(rxcm*vycm-rycm*vxcm)*masstot
        !
        ! calculate and then invert the inertia tensor
        !
        xxi = 0.0d0
        xyi = 0.0d0
        xzi = 0.0d0
        yyi = 0.0d0
        yzi = 0.0d0
        zzi = 0.0d0
        do i = 1, natom
            xdel = rx(i) - rxcm
            ydel = ry(i) - rycm
            zdel = rz(i) - rzcm
            xxi = xxi + xdel*xdel*massmdqt(i)
            xyi = xyi + xdel*ydel*massmdqt(i)
            xzi = xzi + xdel*zdel*massmdqt(i)
            yyi = yyi + ydel*ydel*massmdqt(i)
            yzi = yzi + ydel*zdel*massmdqt(i)
            zzi = zzi + zdel*zdel*massmdqt(i)
        end do
        tensor(1,1) = yyi + zzi
        tensor(2,1) = -xyi
        tensor(3,1) = -xzi
        tensor(1,2) = -xyi
        tensor(2,2) = xxi + zzi
        tensor(3,2) = -yzi
        tensor(1,3) = -xzi
        tensor(2,3) = -yzi
        tensor(3,3) = xxi + yyi
        call invert(3,3,tensor)
        !
        ! compute angular velocity and rotational kinetic energy
        ! using L=I*Omega
        !
        erot = 0.0d0
        do i = 1, 3
            vang(i) = 0.0d0
            do j = 1, 3
                vang(i) = vang(i) + tensor(i,j)*mang(j)
            end do
            erot = erot + vang(i)*mang(i)
        end do
        erot = 0.5d0 * erot
        !
        ! eliminate any translation of the overall system
        !
        do i = 1, natom
            vx(i)=vx(i)-vxcm
            vy(i)=vy(i)-vycm
            vz(i)=vz(i)-vzcm
        end do
        !
        ! eliminate any rotation about the system center of mass
        !
        do i = 1, natom
            xdel = rx(i) - rxcm
            ydel = ry(i) - rycm
            zdel = rz(i) - rzcm
            vx(i) = vx(i) - vang(2)*zdel + vang(3)*ydel
            vy(i) = vy(i) - vang(3)*xdel + vang(1)*zdel
            vz(i) = vz(i) - vang(1)*ydel + vang(2)*xdel
        end do
    end subroutine
    ! ############################################################
    ! ## ##
    ! ## subroutine invert -- gauss-jordan matrix inversion ##
    ! ## ##
    ! ############################################################
    !
    !
    ! "invert" inverts a matrix using the Gauss-Jordan method
    !
    ! variables and parameters:
    !
    ! n logical dimension of the matrix to be inverted
    ! np physical dimension of the matrix storage area
    ! a matrix to invert; contains inverse on exit
    !
    !
    subroutine invert (n,np,a)
        implicit none
        integer maxinv
        parameter (maxinv=100)
        integer i,j,k,n,np
        integer icol,irow
        integer ipivot(maxinv)
        integer indxc(maxinv)
        integer indxr(maxinv)
        _REAL_ big,temp,pivot
        _REAL_ a(np,np)
        !
        ! perform matrix inversion via the Gauss-Jordan algorithm
        !
        do i = 1, n
            ipivot(i) = 0
        end do
        do i = 1, n
            big = 0.0d0
            do j = 1, n
                if (ipivot(j) .ne. 1) then
                    do k = 1, n
                        if (ipivot(k) .eq. 0) then
                            if (abs(a(j,k)) .ge. big) then
                                big = abs(a(j,k))
                                irow = j
                                icol = k
                            end if
                        end if
                    end do
                end if
            end do
            ipivot(icol) = ipivot(icol) + 1
            if (irow .ne. icol) then
                do j = 1, n
                    temp = a(irow,j)
                    a(irow,j) = a(icol,j)
                    a(icol,j) = temp
                end do
            end if
            indxr(i) = irow
            indxc(i) = icol
            pivot = a(icol,icol)
            a(icol,icol) = 1.0d0
            do j = 1, n
                a(icol,j) = a(icol,j) / pivot
            end do
            do j = 1, n
                if (j .ne. icol) then
                    temp = a(j,icol)
                    a(j,icol) = 0.0d0
                    do k = 1, n
                        a(j,k) = a(j,k) - a(icol,k)*temp
                    end do
                end if
            end do
        end do
        do i = n, 1, -1
            if (indxr(i) .ne. indxc(i)) then
                do k = 1, n
                    temp = a(k,indxr(i))
                    a(k,indxr(i)) = a(k,indxc(i))
                    a(k,indxc(i)) = temp
                end do
            end if
        end do
        return
    end subroutine
