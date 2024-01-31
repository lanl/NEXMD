!**********************************************************************!
!*This module contains auxiliary subroutines for freezing normal modes*!
!**********************************************************************!
module freezing_module
    use communism
    implicit none

contains

    subroutine rmsd(rx1, ry1, rz1, rx2, ry2, rz2, natom, r)

        implicit none

        type(simulation_t), pointer :: sim

        integer i, j
        integer natom
        real rx1(natom), ry1(natom), rz1(natom)
        real rx2(natom), ry2(natom), rz2(natom)
        real r

        r = 0.0d0
        do i = 1, natom
            r = r + ((rx1(i) - rx2(i))**2 + (ry1(i) - ry2(i))**2 + (rz1(i) - rz2(i))**2)/natom
        end do
        r = sqrt(r)

    end subroutine rmsd

    subroutine freezePairs1(sim)!Subroutine for verlet1 for freezing distances between pairs of atoms

        implicit none

        type(simulation_t), pointer :: sim

        integer i, j, l
        real :: mx(sim%naesmd%natom), my(sim%naesmd%natom), mz(sim%naesmd%natom)
        real :: mxo(sim%naesmd%natom), myo(sim%naesmd%natom), mzo(sim%naesmd%natom)
        real :: sij(sim%naesmd%natom, sim%naesmd%natom), gg(sim%naesmd%natom, sim%naesmd%natom)
        real :: dij(sim%naesmd%natom, sim%naesmd%natom), dis(sim%naesmd%natom, sim%naesmd%natom)
        real dt

        dt = sim%naesmd%dtmdqt

        do i = 1, sim%naesmd%natom
            mx(i) = sim%naesmd%rx(i)
            my(i) = sim%naesmd%ry(i)
            mz(i) = sim%naesmd%rz(i)
            mxo(i) = sim%naesmd%rxold(i)
            myo(i) = sim%naesmd%ryold(i)
            mzo(i) = sim%naesmd%rzold(i)
        end do

        do l = 1, sim%naesmd%npc
            i = sim%naesmd%dtc(l, 1)
            j = sim%naesmd%dtc(l, 2)
            dis(i, j) = (mxo(i) - mxo(j))**2 + (myo(i) - myo(j))**2 + (mzo(i) - mzo(j))**2
            sij(i, j) = (mx(i) - mx(j))**2 + (my(i) - my(j))**2 + (mz(i) - mz(j))**2
            if (abs(sij(i, j) - dis(i, j)) .gt. 10.0d-17) then
                gg(i, j) = (sij(i, j) - dis(i, j))/(2.0d0*dt*((mx(i) - mx(j))*(mxo(i) - mxo(j)) + &
                    (my(i) - my(j))*(myo(i) - myo(j)) + (mz(i) - mz(j))*(mzo(i) - mzo(j)))*(1.0d0/sim%naesmd%massmdqt(i) + &
                    1.0d0/sim%naesmd%massmdqt(j)))
                mx(i) = mx(i) - dt*gg(i, j)*(mx(i) - mx(j))/sim%naesmd%massmdqt(i)
                mx(j) = mx(j) + dt*gg(i, j)*(mx(i) - mx(j))/sim%naesmd%massmdqt(j)
                my(i) = my(i) - dt*gg(i, j)*(my(i) - my(j))/sim%naesmd%massmdqt(i)
                my(j) = my(j) + dt*gg(i, j)*(my(i) - my(j))/sim%naesmd%massmdqt(j)
                mz(i) = mz(i) - dt*gg(i, j)*(mz(i) - mz(j))/sim%naesmd%massmdqt(i)
                mz(j) = mz(j) + dt*gg(i, j)*(mz(i) - mz(j))/sim%naesmd%massmdqt(j)
                mxo(i) = mx(i)
                mxo(j) = mx(j)
                myo(i) = my(i)
                myo(j) = my(j)
                mzo(i) = mz(i)
                mzo(j) = mz(j)
                sim%naesmd%vx(i) = sim%naesmd%vx(i) - gg(i, j)*(mx(i) - mx(j))/sim%naesmd%massmdqt(i)
                sim%naesmd%vy(i) = sim%naesmd%vy(i) - gg(i, j)*(my(i) - my(j))/sim%naesmd%massmdqt(i)
                sim%naesmd%vz(i) = sim%naesmd%vz(i) - gg(i, j)*(mz(i) - mz(j))/sim%naesmd%massmdqt(i)
                sim%naesmd%vx(j) = sim%naesmd%vx(j) + gg(i, j)*(mx(i) - mx(j))/sim%naesmd%massmdqt(j)
                sim%naesmd%vy(j) = sim%naesmd%vy(j) + gg(i, j)*(my(i) - my(j))/sim%naesmd%massmdqt(j)
                sim%naesmd%vz(j) = sim%naesmd%vz(j) + gg(i, j)*(mz(i) - mz(j))/sim%naesmd%massmdqt(j)
            end if
            sim%naesmd%rx(i) = mx(i)
            sim%naesmd%rx(j) = mx(j)
            sim%naesmd%ry(i) = my(i)
            sim%naesmd%ry(j) = my(j)
            sim%naesmd%rz(i) = mz(i)
            sim%naesmd%rz(j) = mz(j)
        end do

    end subroutine freezePairs1

    subroutine freezePairs2(sim)!Subroutine for verlet2 for freezing distances between pairs of atoms

        implicit none

        type(simulation_t), pointer :: sim

        integer i, j, l
        real :: mx(sim%naesmd%natom), my(sim%naesmd%natom), mz(sim%naesmd%natom)
        real :: mxo(sim%naesmd%natom), myo(sim%naesmd%natom), mzo(sim%naesmd%natom)
        real :: sij(sim%naesmd%natom, sim%naesmd%natom), gg(sim%naesmd%natom, sim%naesmd%natom)
        real :: dij(sim%naesmd%natom, sim%naesmd%natom), kk(sim%naesmd%natom, sim%naesmd%natom)
        real :: dis(sim%naesmd%natom, sim%naesmd%natom)
        real dt

        dt = sim%naesmd%dtmdqt

        do i = 1, sim%naesmd%natom
            mx(i) = sim%naesmd%rx(i)
            my(i) = sim%naesmd%ry(i)
            mz(i) = sim%naesmd%rz(i)
        end do
        do l = 1, sim%naesmd%npc
            i = sim%naesmd%dtc(l, 1)
            j = sim%naesmd%dtc(l, 2)
            dis(i, j) = (mx(i) - mx(j))**2 + (my(i) - my(j))**2 + (mz(i) - mz(j))**2
        end do

        do l = 1, sim%naesmd%npc
            i = sim%naesmd%dtc(l, 1)
            j = sim%naesmd%dtc(l, 2)
200         continue
            if (abs((mx(i) - mx(j))*(sim%naesmd%vx(i) - sim%naesmd%vx(j)) + (my(i) - my(j))*(sim%naesmd%vy(i) - sim%naesmd%vy(j)) + &
                (mz(i) - mz(j))*(sim%naesmd%vz(i) - sim%naesmd%vz(j))) .gt. 10d-17) then
                kk(i, j) = ((mx(i) - mx(j))*(sim%naesmd%vx(i) - sim%naesmd%vx(j)) + (my(i) - my(j))*(sim%naesmd%vy(i) - sim%naesmd%vy(j)) + &
                    (mz(i) - mz(j))*(sim%naesmd%vz(i) - sim%naesmd%vz(j)))/(dis(i, j)*(1.0d0/sim%naesmd%massmdqt(i) + 1.0d0/sim%naesmd%massmdqt(j)))
                sim%naesmd%vx(i) = sim%naesmd%vx(i) - kk(i, j)*(mx(i) - mx(j))/sim%naesmd%massmdqt(i)
                sim%naesmd%vy(i) = sim%naesmd%vy(i) - kk(i, j)*(my(i) - my(j))/sim%naesmd%massmdqt(i)
                sim%naesmd%vz(i) = sim%naesmd%vz(i) - kk(i, j)*(mz(i) - mz(j))/sim%naesmd%massmdqt(i)
                sim%naesmd%vx(j) = sim%naesmd%vx(j) + kk(i, j)*(mx(i) - mx(j))/sim%naesmd%massmdqt(j)
                sim%naesmd%vy(j) = sim%naesmd%vy(j) + kk(i, j)*(my(i) - my(j))/sim%naesmd%massmdqt(j)
                sim%naesmd%vz(j) = sim%naesmd%vz(j) + kk(i, j)*(mz(i) - mz(j))/sim%naesmd%massmdqt(j)
                go to 200
            end if
        end do

    end subroutine freezePairs2

    subroutine freeze(sim) !For freezing normal modes

        implicit none

        type(simulation_t), pointer :: sim

        integer i, j, l, k, kk
        real r1ini(sim%naesmd%natom, 3), r2ini(sim%naesmd%natom, 3)
        real r3ini(sim%naesmd%natom, 3), r4ini(sim%naesmd%natom, 3)
        real r5ini(sim%naesmd%natom, 3)
        real mx(sim%naesmd%natom), my(sim%naesmd%natom), mz(sim%naesmd%natom)
        real cm1(3), cm2(3), vcm(3), acm(3), cm(3)
        real r1(sim%naesmd%natom, 3), r2(sim%naesmd%natom, 3)
        real r3(sim%naesmd%natom, 3), r4(sim%naesmd%natom, 3)
        real r5(sim%naesmd%natom, 3)
        real xxyx, xxyy, xxyz, xyyx, xyyy, xzyx, xzyy, xzyz, xyyz
        real*8 c(4, 4), v(4, 4), work1(4), work2(4), d(4), q(4)
        real*8 rot(3, 3)
        real xrot, yrot, zrot
        real ldmode(sim%naesmd%natom*3), gmmode(sim%naesmd%natom*3)
        real dt2_2, dt_2
        real r

        dt_2 = sim%naesmd%dtmdqt/2.d0 ! dt/2
        dt2_2 = sim%naesmd%dtmdqt**2/2 ! dt^2/2

!Equilibrium positions
        do i = 1, sim%naesmd%natom
            r1ini(i, 1) = sim%naesmd%xbf0(i)
            r1ini(i, 2) = sim%naesmd%ybf0(i)
            r1ini(i, 3) = sim%naesmd%zbf0(i)
        end do
!Positions
        do i = 1, sim%naesmd%natom
            r2ini(i, 1) = sim%naesmd%rx(i)
            r2ini(i, 2) = sim%naesmd%ry(i)
            r2ini(i, 3) = sim%naesmd%rz(i)
        end do
!Updated positions without freezing
        do i = 1, sim%naesmd%natom
            mx(i) = sim%naesmd%vxold(i)*sim%naesmd%dtmdqt + sim%naesmd%ax(i)*dt2_2
            my(i) = sim%naesmd%vyold(i)*sim%naesmd%dtmdqt + sim%naesmd%ay(i)*dt2_2
            mz(i) = sim%naesmd%vzold(i)*sim%naesmd%dtmdqt + sim%naesmd%az(i)*dt2_2
        end do
!*****************************************************
! Calculo del CM de los atomos de la estructura 1
!*****************************************************
        do j = 1, 3
            cm1(j) = 0.0d0
        end do
        do i = 1, sim%naesmd%natom
            do j = 1, 3
                cm1(j) = cm1(j) + sim%naesmd%massmdqt(i)*r1ini(i, j)/sim%naesmd%masatotal
            end do
        end do
!*****************************************************
! Traslacion de los atomos de la estructura 1
!*****************************************************
        do j = 1, sim%naesmd%natom
            do kk = 1, 3
                r1(j, kk) = r1ini(j, kk) - cm1(kk)
            end do
        end do
!*****************************************************
! Calculo del CM de los atomos de la estructura 2
!*****************************************************
        do j = 1, 3
            cm2(j) = 0.0d0
        end do
        do i = 1, sim%naesmd%natom
            do j = 1, 3
                cm2(j) = cm2(j) + sim%naesmd%massmdqt(i)*r2ini(i, j)/sim%naesmd%masatotal
            end do
        end do
!*****************************************************
! Traslacion de los atomos de la estructura 2
!*****************************************************
        do j = 1, sim%naesmd%natom
            do kk = 1, 3
                r2(j, kk) = r2ini(j, kk) - cm2(kk)
            end do
        end do
!        call rmsd(r1(:,1),r1(:,2),r1(:,3),r2(:,1),r2(:,2),r2(:,3),sim%naesmd%natom,r)
!        write(6,*) 'RMSD before rotation: ', r
!Rotation:
        xxyx = 0.0d0
        xxyy = 0.0d0
        xxyz = 0.0d0
        xyyx = 0.0d0
        xyyy = 0.0d0
        xyyz = 0.0d0
        xzyx = 0.0d0
        xzyy = 0.0d0
        xzyz = 0.0d0
        do j = 1, sim%naesmd%natom
            xxyx = xxyx + r1(j, 1)*r2(j, 1)
            xxyy = xxyy + r1(j, 2)*r2(j, 1)
            xxyz = xxyz + r1(j, 3)*r2(j, 1)
            xyyx = xyyx + r1(j, 1)*r2(j, 2)
            xyyy = xyyy + r1(j, 2)*r2(j, 2)
            xyyz = xyyz + r1(j, 3)*r2(j, 2)
            xzyx = xzyx + r1(j, 1)*r2(j, 3)
            xzyy = xzyy + r1(j, 2)*r2(j, 3)
            xzyz = xzyz + r1(j, 3)*r2(j, 3)
        end do
        c(1, 1) = xxyx + xyyy + xzyz
        c(1, 2) = xzyy - xyyz
        c(2, 2) = xxyx - xyyy - xzyz
        c(1, 3) = xxyz - xzyx
        c(2, 3) = xxyy + xyyx
        c(3, 3) = xyyy - xzyz - xxyx
        c(1, 4) = xyyx - xxyy
        c(2, 4) = xzyx + xxyz
        c(3, 4) = xyyz + xzyy
        c(4, 4) = xzyz - xxyx - xyyy

!     diagonalize the quadratic form matrix

        call jacobi(4, 4, c, d, v, work1, work2)

!     extract the desired quaternion
!
        q(1) = v(1, 4)
        q(2) = v(2, 4)
        q(3) = v(3, 4)
        q(4) = v(4, 4)

!     assemble rotation matrix that superimposes the molecules
!
        rot(1, 1) = q(1)**2 + q(2)**2 - q(3)**2 - q(4)**2
        rot(2, 1) = 2.0d0*(q(2)*q(3) - q(1)*q(4))
        rot(3, 1) = 2.0d0*(q(2)*q(4) + q(1)*q(3))
        rot(1, 2) = 2.0d0*(q(3)*q(2) + q(1)*q(4))
        rot(2, 2) = q(1)**2 - q(2)**2 + q(3)**2 - q(4)**2
        rot(3, 2) = 2.0d0*(q(3)*q(4) - q(1)*q(2))
        rot(1, 3) = 2.0d0*(q(4)*q(2) - q(1)*q(3))
        rot(2, 3) = 2.0d0*(q(4)*q(3) + q(1)*q(2))
        rot(3, 3) = q(1)**2 - q(2)**2 - q(3)**2 + q(4)**2

!
!     rotate second molecule to best fit with first molecule
!

        do j = 1, sim%naesmd%natom
            xrot = r2(j, 1)*rot(1, 1) + r2(j, 2)*rot(1, 2) + r2(j, 3)*rot(1, 3)
            yrot = r2(j, 1)*rot(2, 1) + r2(j, 2)*rot(2, 2) + r2(j, 3)*rot(2, 3)
            zrot = r2(j, 1)*rot(3, 1) + r2(j, 2)*rot(3, 2) + r2(j, 3)*rot(3, 3)
            sim%naesmd%rx(j) = xrot
            sim%naesmd%ry(j) = yrot
            sim%naesmd%rz(j) = zrot
            r2(j, 1) = xrot
            r2(j, 2) = yrot
            r2(j, 3) = zrot
        end do
!        call rmsd(r1(:,1),r1(:,2),r1(:,3),r2(:,1),r2(:,2),r2(:,3),sim%naesmd%natom,r)
!        write(6,*) 'RMSD after rotation: ', r

!******************************************************
! RMSD VELOCIDADES, ACELERACION, X-Y-Z
!******************************************************
        do i = 1, sim%naesmd%natom
            r3ini(i, 1) = sim%naesmd%vx(i)
            r3ini(i, 2) = sim%naesmd%vy(i)
            r3ini(i, 3) = sim%naesmd%vz(i)
        end do

        do i = 1, sim%naesmd%natom
            r4ini(i, 1) = sim%naesmd%ax(i)
            r4ini(i, 2) = sim%naesmd%ay(i)
            r4ini(i, 3) = sim%naesmd%az(i)
        end do

        do i = 1, sim%naesmd%natom
            r5ini(i, 1) = mx(i)
            r5ini(i, 2) = my(i)
            r5ini(i, 3) = mz(i)
        end do

!*****************************************************
! Calculo de VCM y ACM de los atomos de la estructura 2
!*****************************************************
        do j = 1, 3
            vcm(j) = 0.0d0
        end do
        do i = 1, sim%naesmd%natom
            do j = 1, 3
                vcm(j) = vcm(j) + sim%naesmd%massmdqt(i)*r3ini(i, j)/sim%naesmd%masatotal
            end do
        end do

        do j = 1, 3
            acm(j) = 0.0d0
        end do
        do i = 1, sim%naesmd%natom
            do j = 1, 3
                acm(j) = acm(j) + sim%naesmd%massmdqt(i)*r4ini(i, j)/sim%naesmd%masatotal
            end do
        end do

        do j = 1, 3
            cm(j) = 0.0d0
        end do
        do i = 1, sim%naesmd%natom
            do j = 1, 3
                cm(j) = cm(j) + sim%naesmd%massmdqt(i)*r5ini(i, j)/sim%naesmd%masatotal
            end do
        end do

!*****************************************************
! Traslacion de los atomos de la estructura 2
!*****************************************************
        do j = 1, sim%naesmd%natom
            do kk = 1, 3
                r3(j, kk) = r3ini(j, kk) - vcm(kk)
            end do
        end do

        do j = 1, sim%naesmd%natom
            do kk = 1, 3
                r4(j, kk) = r4ini(j, kk) - acm(kk)
            end do
        end do

        do j = 1, sim%naesmd%natom
            do kk = 1, 3
                r5(j, kk) = r5ini(j, kk) - cm(kk)
            end do
        end do
        do j = 1, sim%naesmd%natom
            xrot = r3(j, 1)*rot(1, 1) + r3(j, 2)*rot(1, 2) + r3(j, 3)*rot(1, 3)
            yrot = r3(j, 1)*rot(2, 1) + r3(j, 2)*rot(2, 2) + r3(j, 3)*rot(2, 3)
            zrot = r3(j, 1)*rot(3, 1) + r3(j, 2)*rot(3, 2) + r3(j, 3)*rot(3, 3)
            sim%naesmd%vx(j) = xrot
            sim%naesmd%vy(j) = yrot
            sim%naesmd%vz(j) = zrot
            r3(j, 1) = xrot
            r3(j, 2) = yrot
            r3(j, 3) = zrot

            xrot = r4(j, 1)*rot(1, 1) + r4(j, 2)*rot(1, 2) + r4(j, 3)*rot(1, 3)
            yrot = r4(j, 1)*rot(2, 1) + r4(j, 2)*rot(2, 2) + r4(j, 3)*rot(2, 3)
            zrot = r4(j, 1)*rot(3, 1) + r4(j, 2)*rot(3, 2) + r4(j, 3)*rot(3, 3)
            sim%naesmd%ax(j) = xrot
            sim%naesmd%ay(j) = yrot
            sim%naesmd%az(j) = zrot
            r4(j, 1) = xrot
            r4(j, 2) = yrot
            r4(j, 3) = zrot

            xrot = r5(j, 1)*rot(1, 1) + r5(j, 2)*rot(1, 2) + r5(j, 3)*rot(1, 3)
            yrot = r5(j, 1)*rot(2, 1) + r5(j, 2)*rot(2, 2) + r5(j, 3)*rot(2, 3)
            zrot = r5(j, 1)*rot(3, 1) + r5(j, 2)*rot(3, 2) + r5(j, 3)*rot(3, 3)
            mx(j) = xrot
            my(j) = yrot
            mz(j) = zrot
            r5(j, 1) = xrot
            r5(j, 2) = yrot
            r5(j, 3) = zrot
        end do

!Calculate Lambdas
        ldmode = 0.0d0
        do l = 1, sim%naesmd%nmc
            kk = sim%naesmd%mc(l)

!               ldmode(kk)=0.0d0
            do k = 1, sim%naesmd%natom
                ldmode(kk) = ldmode(kk) + &
                    sim%naesmd%massqrt(k)*sim%naesmd%enm(3*k - 2, kk)*mx(k) + &
                    sim%naesmd%massqrt(k)*sim%naesmd%enm(3*k - 1, kk)*my(k) + &
                    sim%naesmd%massqrt(k)*sim%naesmd%enm(3*k, kk)*mz(k)
            end do
            ldmode(kk) = ldmode(kk)/dt2_2
        end do

!            write(6,*) 'Lambda:', ldmode(144), sim%naesmd%massqrt(1), sim%naesmd%enm(1,1), mx(1)

!Calculate Gammas
        gmmode = 0.0d0
        do l = 1, sim%naesmd%nmc
            kk = sim%naesmd%mc(l)
            gmmode(kk) = 0.0d0
            do k = 1, sim%naesmd%natom
                gmmode(kk) = gmmode(kk) + &
                    sim%naesmd%massqrt(k)*sim%naesmd%enm(3*k - 2, kk)*sim%naesmd%vx(k) + &
                    sim%naesmd%massqrt(k)*sim%naesmd%enm(3*k - 1, kk)*sim%naesmd%vy(k) + &
                    sim%naesmd%massqrt(k)*sim%naesmd%enm(3*k, kk)*sim%naesmd%vz(k)
            end do
            gmmode(kk) = gmmode(kk)/dt_2   !!!!!!!!!!!!!!!Gamma+Lambda
        end do

!            write(6,*) 'Gamma:', gmmode

!Freezing
        do k = 1, sim%naesmd%natom
            do l = 1, sim%naesmd%nmc
                j = sim%naesmd%mc(l)
                sim%naesmd%rx(k) = sim%naesmd%rx(k) - &
                    dt2_2*ldmode(j)*sim%naesmd%enm(3*k - 2, j)/sim%naesmd%massqrt(k)
                sim%naesmd%ry(k) = sim%naesmd%ry(k) - &
                    dt2_2*ldmode(j)*sim%naesmd%enm(3*k - 1, j)/sim%naesmd%massqrt(k)
                sim%naesmd%rz(k) = sim%naesmd%rz(k) - &
                    dt2_2*ldmode(j)*sim%naesmd%enm(3*k, j)/sim%naesmd%massqrt(k)
            end do
        end do

        do k = 1, sim%naesmd%natom
            do l = 1, sim%naesmd%nmc
                j = sim%naesmd%mc(l)
                sim%naesmd%vx(k) = sim%naesmd%vx(k) - &
                    dt_2*gmmode(j)*sim%naesmd%enm(3*k - 2, j)/sim%naesmd%massqrt(k)
                sim%naesmd%vy(k) = sim%naesmd%vy(k) - &
                    dt_2*gmmode(j)*sim%naesmd%enm(3*k - 1, j)/sim%naesmd%massqrt(k)
                sim%naesmd%vz(k) = sim%naesmd%vz(k) - &
                    dt_2*gmmode(j)*sim%naesmd%enm(3*k, j)/sim%naesmd%massqrt(k)
            end do
        end do

!Backward rotation
        do j = 1, sim%naesmd%natom
            r3(j, 1) = sim%naesmd%vx(j)
            r3(j, 2) = sim%naesmd%vy(j)
            r3(j, 3) = sim%naesmd%vz(j)
            xrot = r3(j, 1)*rot(1, 1) + r3(j, 2)*rot(2, 1) + r3(j, 3)*rot(3, 1)
            yrot = r3(j, 1)*rot(1, 2) + r3(j, 2)*rot(2, 2) + r3(j, 3)*rot(3, 2)
            zrot = r3(j, 1)*rot(1, 3) + r3(j, 2)*rot(2, 3) + r3(j, 3)*rot(3, 3)
            sim%naesmd%vx(j) = xrot + vcm(1)
            sim%naesmd%vy(j) = yrot + vcm(2)
            sim%naesmd%vz(j) = zrot + vcm(3)

            r4(j, 1) = sim%naesmd%ax(j)
            r4(j, 2) = sim%naesmd%ay(j)
            r4(j, 3) = sim%naesmd%az(j)
            xrot = r4(j, 1)*rot(1, 1) + r4(j, 2)*rot(2, 1) + r4(j, 3)*rot(3, 1)
            yrot = r4(j, 1)*rot(1, 2) + r4(j, 2)*rot(2, 2) + r4(j, 3)*rot(3, 2)
            zrot = r4(j, 1)*rot(1, 3) + r4(j, 2)*rot(2, 3) + r4(j, 3)*rot(3, 3)
            sim%naesmd%ax(j) = xrot + acm(1)
            sim%naesmd%ay(j) = yrot + acm(2)
            sim%naesmd%az(j) = zrot + acm(3)

            r2(j, 1) = sim%naesmd%rx(j)
            r2(j, 2) = sim%naesmd%ry(j)
            r2(j, 3) = sim%naesmd%rz(j)
            xrot = r2(j, 1)*rot(1, 1) + r2(j, 2)*rot(2, 1) + r2(j, 3)*rot(3, 1)
            yrot = r2(j, 1)*rot(1, 2) + r2(j, 2)*rot(2, 2) + r2(j, 3)*rot(3, 2)
            zrot = r2(j, 1)*rot(1, 3) + r2(j, 2)*rot(2, 3) + r2(j, 3)*rot(3, 3)
            sim%naesmd%rx(j) = xrot + cm2(1)
            sim%naesmd%ry(j) = yrot + cm2(2)
            sim%naesmd%rz(j) = zrot + cm2(3)
        end do

    end subroutine freeze

!
!     ###################################################
!     ##  COPYRIGHT (C)  1990  by  Jay William Ponder  ##
!     ##              All Rights Reserved              ##
!     ###################################################
!
!     ############################################################
!     ##                                                        ##
!     ##  subroutine jacobi  --  jacobi matrix diagonalization  ##
!     ##                                                        ##
!     ############################################################
!
!
!     "jacobi" performs a matrix diagonalization of a real
!     symmetric matrix by the method of Jacobi rotations
!
!     variables and parameters:
!
!     n     logical dimension of the matrix to be diagonalized
!     np    physical dimension of the matrix storage area
!     a     input with the matrix to be diagonalized; only
!              the upper triangle and diagonal are required
!     d     returned with the eigenvalues in ascending order
!     v     returned with the eigenvectors of the matrix
!     b     temporary work vector
!     z     temporary work vector
!
!
    subroutine jacobi(n, np, a, d, v, b, z)
        implicit none
        integer i, j, k
        integer n, np, ip, iq
        integer nrot, maxrot
        real*8 sm, tresh, s, c, t
        real*8 theta, tau, h, g, p
        real*8 d(np), b(np), z(np)
        real*8 a(np, np), v(np, np)
!
!
!     setup and initialization
!
        maxrot = 100
        nrot = 0
        do ip = 1, n
            do iq = 1, n
                v(ip, iq) = 0.0d0
            end do
            v(ip, ip) = 1.0d0
        end do
        do ip = 1, n
            b(ip) = a(ip, ip)
            d(ip) = b(ip)
            z(ip) = 0.0d0
        end do
!
!     perform the jacobi rotations
!
        do i = 1, maxrot
            sm = 0.0d0
            do ip = 1, n - 1
                do iq = ip + 1, n
                    sm = sm + abs(a(ip, iq))
                end do
            end do
            if (sm .eq. 0.0d0) goto 10
            if (i .lt. 4) then
                tresh = 0.2d0*sm/n**2
            else
                tresh = 0.0d0
            end if
            do ip = 1, n - 1
                do iq = ip + 1, n
                    g = 100.0d0*abs(a(ip, iq))
                    if (i .gt. 4 .and. abs(d(ip)) + g .eq. abs(d(ip)) &
                        .and. abs(d(iq)) + g .eq. abs(d(iq))) then
                        a(ip, iq) = 0.0d0
                    else if (abs(a(ip, iq)) .gt. tresh) then
                        h = d(iq) - d(ip)
                        if (abs(h) + g .eq. abs(h)) then
                            t = a(ip, iq)/h
                        else
                            theta = 0.5d0*h/a(ip, iq)
                            t = 1.0d0/(abs(theta) + sqrt(1.0d0 + theta**2))
                            if (theta .lt. 0.0d0) t = -t
                        end if
                        c = 1.0d0/sqrt(1.0d0 + t**2)
                        s = t*c
                        tau = s/(1.0d0 + c)
                        h = t*a(ip, iq)
                        z(ip) = z(ip) - h
                        z(iq) = z(iq) + h
                        d(ip) = d(ip) - h
                        d(iq) = d(iq) + h
                        a(ip, iq) = 0.0d0
                        do j = 1, ip - 1
                            g = a(j, ip)
                            h = a(j, iq)
                            a(j, ip) = g - s*(h + g*tau)
                            a(j, iq) = h + s*(g - h*tau)
                        end do
                        do j = ip + 1, iq - 1
                            g = a(ip, j)
                            h = a(j, iq)
                            a(ip, j) = g - s*(h + g*tau)
                            a(j, iq) = h + s*(g - h*tau)
                        end do
                        do j = iq + 1, n
                            g = a(ip, j)
                            h = a(iq, j)
                            a(ip, j) = g - s*(h + g*tau)
                            a(iq, j) = h + s*(g - h*tau)
                        end do
                        do j = 1, n
                            g = v(j, ip)
                            h = v(j, iq)
                            v(j, ip) = g - s*(h + g*tau)
                            v(j, iq) = h + s*(g - h*tau)
                        end do
                        nrot = nrot + 1
                    end if
                end do
            end do
            do ip = 1, n
                b(ip) = b(ip) + z(ip)
                d(ip) = b(ip)
                z(ip) = 0.0d0
            end do
        end do
!
!     print warning if not converged
!
10      continue
        if (nrot .eq. maxrot) then
            write (99, 20)
20          format(/, ' JACOBI  --  Matrix Diagonalization not Converged')
        end if
!
!     sort the eigenvalues and vectors
!
        do i = 1, n - 1
            k = i
            p = d(i)
            do j = i + 1, n
                if (d(j) .lt. p) then
                    k = j
                    p = d(j)
                end if
            end do
            if (k .ne. i) then
                d(k) = d(i)
                d(i) = p
                do j = 1, n
                    p = v(j, i)
                    v(j, i) = v(j, k)
                    v(j, k) = p
                end do
            end if
        end do
        return
    end subroutine jacobi

end module freezing_module
