#include "dprec.fh"
#include "assert.fh"

!-------------------------------------------------------------!
!This subroutine checks that files are prepared for restarting!
!-------------------------------------------------------------!
subroutine check_files(Nsim,nstep0,id,id0)

    use communism

    implicit none

!Passed in
    integer, intent(in) :: Nsim
    integer, intent(out) :: nstep0
    integer, intent(in) :: id
    integer, intent(in) :: id0

!Internal
    integer inputSteps
    integer restartSteps
    integer steps
    integer aimcSteps
    integer l,l1,l2,i
    character*150 filename
    character*150 dynam_type
    character*10 strl
    character*9 lastLine
    integer excN
    integer cs
    integer outData,outDataCoords,natoms,printTdipole
    integer stepClone
    double precision timeClone
    integer dbl
    _REAL_ dt
    integer, allocatable, dimension(:) :: clone_steps
    integer, allocatable, dimension(:) :: dropout_steps
    integer nlines
    integer lprint
    integer nq
    integer nclones0
    integer ired
    integer dps
    logical file_exists2

    dps = 0
    inquire(FILE="dropped.out", EXIST=file_exists2)
    if(file_exists2) then
        call system('wc dropped.out > dropped.tmp')
        open(1,file='dropped.tmp')
            read(1,*) dps
        close(1)
        call system('rm dropped.tmp')
    endif

    if(dps.gt.1) then
        allocate(dropout_steps(dps))
        allocate(clone_steps(Nsim + dps))
        dropout_steps = -1
    else
        allocate(dropout_steps(1))
        allocate(clone_steps(Nsim))
        dropout_steps = -1
    endif

    filename = ''

!Reading steps from input
    call read_input_moldyn_block("input.ceon",inputSteps,dynam_type,excN,cs,&
        outData,outDataCoords,natoms,printTdipole,dt,lprint,nq,nclones0,ired) 
!Checking synchronization
    if(trim(adjustl(dynam_type)).eq.'aimc'.and.Nsim.gt.1) call check_synchronization(Nsim)
    nstep0=inputSteps
    stepClone = 0
    timeClone = 0.0d0
!Checking restart file
    !Calculating time of cloning
    if(trim(adjustl(dynam_type)).eq.'aimc'.and.Nsim.gt.1) then
        filename = ''
        write (filename, "(a8,i4.4,a4)") "restart_", id, ".out"
        call system("tail -1 "//filename(1:16)//" > file1.tmp")
    else
        call system("tail -1 restart.out > file1.tmp")
    endif
    open(1,file="file1.tmp")
        read(1,*) lastLine
    close(1)
    call system("rm file1.tmp")
    if(lastLine.ne.'&endcoeff') then
        print *, "Restart file is incomplete!"
!        call remake_restart(Nsim) !All restart files has to be remaked
!        print *, 'Restart files remaked, try launching again!'
        STOP;
    else
        if(Nsim.eq.1) then
            print *, 'Restart file complete'
        else
            print *, trim(adjustl(filename))//' file complete'
        endif
    endif
!Reading steps from restart
    if(Nsim.eq.1) then
        call read_input_moldyn_block("restart.out",restartSteps,dynam_type,excN,cs,&
            outData,outDataCoords,natoms,printTdipole,dt,lprint,nq,nclones0,ired)
    else
        filename = ''
        write (filename, "(a8,i4.4,a4)") "restart_", id, ".out"
        call read_input_moldyn_block(filename(1:16),restartSteps,dynam_type,excN,cs,&
            outData,outDataCoords,natoms,printTdipole,dt,lprint,nq,nclones0,ired)
        write (filename, "(a8,i4.4,a4)") "coeff-n_", id, ".out"
        call system('head -1 '//trim(adjustl(filename))//' > file.tmp')
        open(1,file='file.tmp')
        read(1,*) dbl, timeClone
        close(1)
        call system('rm file.tmp')
        stepClone = nint(timeClone/dt/outData*outDataCoords)
    endif
    steps = (inputSteps-restartSteps)
!        print *, 'Checking files for restart'
!        print *, 'Total number of classical steps:', inputSteps
!        print *, 'Number of steps already finished:', steps
!        print *, 'Remaining steps:', restartSteps
    if(restartSteps.le.0) then
        print *, 'MD already finished :)'
        STOP;
    endif
    aimcSteps=steps/outData*outDataCoords+1
    steps=steps/outData*outDataCoords+1-stepClone
 
    if(excN.gt.0) then
!Checking coefficient.out file
        if(cs.gt.0) then
        if(ired.eq.0) then
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a12,i4.4,a4)") "coefficient_", id, ".out"
            call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
        else
            call system("grep -c $ coefficient.out > file.tmp")
        endif
        open(1,file='file.tmp')
        read(1,*) l
        close(1)
        call system("rm file.tmp")
        write(strl, '(I10)') steps*(excN+3)
        if(l.eq.steps*(excN+3)) then
!            print *, 'coefficient.out is ready for restart!'
        elseif(l.lt.steps*(excN+3)) then
            print *, 'coefficient.out file is incomplete, edit restart file and relaunch'
            STOP;
        elseif(l.gt.steps*(excN+3)-stepClone) then
!            print *, 'Cutting exceding lines from coefficient.out'
            if(trim(adjustl(dynam_type)).eq.'aimc') then
                write (filename, "(a12,i4.4,a4)") "coefficient_", id, ".out"
                call system('head -n '//trim(adjustl(strl))//' '//filename//' > coefficient.out_tmp')
                call system('mv coefficient.out_tmp '//filename)
            else
                call system('head -n '//trim(adjustl(strl))//' coefficient.out > coefficient.out_tmp')
                call system('mv coefficient.out_tmp coefficient.out')
            endif
!            print *, 'coefficient.out is ready for restart!'
        endif   
        endif
        endif

!Checking gamma.out file
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a6,i4.4,a4)") "gamma_", id, ".out"
            call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
            open(1,file='file.tmp')
            read(1,*) l
            close(1)
            call system("rm file.tmp")
            write(strl, '(I10)') steps
            if(l.eq.steps) then
!                print *, 'gamma.out is ready for restart!'
            elseif(l.lt.steps) then
                print *, 'gamma.out file is incomplete, edit restart file and relaunch'
                STOP;
            elseif(l.gt.steps) then
!                print *, 'Cutting exceding lines from gamma.out'
                write (filename, "(a6,i4.4,a4)") "gamma_", id, ".out"
                call system('head -n '//trim(adjustl(strl))//' '//filename//' > gamma.out_tmp')
                call system('mv gamma.out_tmp '//filename)
!                print *, 'gamma.out is ready for restart!'
            endif
        endif

!Checking coeff-n.out file
        if(cs.gt.0) then
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a8,i4.4,a4)") "coeff-n_", id, ".out"
            call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
        else
            call system("grep -c $ coeff-n.out > file.tmp")
        endif
        open(1,file='file.tmp')
        read(1,*) l
        close(1)
        call system("rm file.tmp")
        write(strl, '(I10)') steps
        if(l.eq.steps) then
!            print *, 'coeff-n.out is ready for restart!'
        elseif(l.lt.steps) then
            print *, 'coeff-n.out file is incomplete, edit restart file and relaunch'
            STOP;
        elseif(l.gt.steps) then
!            print *, 'Cutting exceding lines from coeff-n.out'
            if(trim(adjustl(dynam_type)).eq.'aimc') then
                write (filename, "(a8,i4.4,a4)") "coeff-n_", id, ".out"
                call system('head -n '//trim(adjustl(strl))//' '//filename//' > coeff-n.out_tmp')
                call system('mv coeff-n.out_tmp '//filename)
            else
                call system('head -n '//trim(adjustl(strl))//' coeff-n.out > coeff-n.out_tmp')
                call system('mv coeff-n.out_tmp coeff-n.out')
            endif
!            print *, 'coeff-n.out is ready for restart!'
        endif
        endif

!Checking coeff-q.out file
        if(cs.gt.0) then
        if(lprint.gt.3) then
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a8,i4.4,a4)") "coeff-q_", id, ".out"
            call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
        else
            call system("grep -c $ coeff-q.out > file.tmp")
        endif
        open(1,file='file.tmp')
        read(1,*) l
        close(1)
        call system("rm file.tmp")
        write(strl, '(I10)') steps
        if(l.eq.steps) then
!            print *, 'coeff-q.out is ready for restart!'
        elseif(l.lt.steps) then
            print *, 'coeff-q.out file is incomplete, edit restart file and relaunch'
            STOP;
        elseif(l.gt.steps) then
!            print *, 'Cutting exceding lines from coeff-q.out'
            if(trim(adjustl(dynam_type)).eq.'aimc') then
                write (filename, "(a8,i4.4,a4)") "coeff-q_", id, ".out"
                call system('head -n '//trim(adjustl(strl))//' '//filename//' > coeff-q.out_tmp')
                call system('mv coeff-q.out_tmp '//filename)
            else
                call system('head -n '//trim(adjustl(strl))//' coeff-q.out > coeff-q.out_tmp')
                call system('mv coeff-q.out_tmp coeff-q.out')
            endif
!            print *, 'coeff-q.out is ready for restart!'
        endif
        endif
    endif
    endif

!Checking coords.xyz
    if(trim(adjustl(dynam_type)).eq.'aimc') then
        write (filename, "(a7,i4.4,a4)") "coords_", id, ".xyz"
        call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
    else
        call system("grep -c $ coords.xyz > file.tmp")
    endif
    open(1,file='file.tmp')
    read(1,*) l
    close(1)
    call system("rm file.tmp")
    write(strl, '(I10)') (natoms+2)*steps
    if(l.eq.(natoms+2)*steps) then
!            print *, 'coords.xyz is ready for restart!'
    elseif(l.lt.(natoms+2)*steps) then
        print *, 'coords.xyz file is incomplete, edit restart file and relaunch'
        STOP;
    elseif(l.gt.(natoms+2)*steps) then
!            print *, 'Cutting exceding lines from coords.xyz'
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a7,i4.4,a4)") "coords_", id, ".xyz"
            call system('head -n '//trim(adjustl(strl))//' '//filename//' > coords.xyz_tmp')
            call system('mv coords.xyz_tmp '//filename)
        else
            call system('head -n '//trim(adjustl(strl))//' coords.xyz > coords.xyz_tmp')
            call system('mv coords.xyz_tmp coords.xyz')
        endif
!            print *, 'coords.xyz is ready for restart!'
    endif

!Checking energy-ev.out
    if(lprint.gt.0) then
    if(trim(adjustl(dynam_type)).eq.'aimc') then
        write (filename, "(a10,i4.4,a4)") "energy-ev_", id, ".out"
        call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
    else
        call system("grep -c $ energy-ev.out > file.tmp")
    endif
    open(1,file='file.tmp')
    read(1,*) l
    close(1)
    call system("rm file.tmp")
    if(trim(adjustl(dynam_type)).eq.'aimc'.and.id.ne.0) then
        l1 = steps
    else
        l1 = steps+1
    endif
    write(strl, '(I10)') l1
    if(l.eq.l1) then
!            print *, 'energy-ev.out is ready for restart!'
    elseif(l.lt.l1) then
        print *, 'energy-ev.out file is incomplete, edit restart file and relaunch'
        STOP;
    elseif(l.gt.l1) then
!            print *, 'Cutting exceding lines from energy-ev.out'
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a10,i4.4,a4)") "energy-ev_", id, ".out"
            call system('head -n '//trim(adjustl(strl))//' '//filename//' > energy-ev.out_tmp')
            call system('mv energy-ev.out_tmp '//filename)
        else
            call system('head -n '//trim(adjustl(strl))//' energy-ev.out > energy-ev.out_tmp')
            call system('mv energy-ev.out_tmp energy-ev.out')
        endif
!            print *, 'energy-ev.out is ready for restart!'
    endif
    endif

!Checking forces.out
    if(lprint.gt.3) then
    if(trim(adjustl(dynam_type)).eq.'aimc') then
        write (filename, "(a7,i4.4,a4)") "forces_", id, ".out"
        call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
    else
        call system("grep -c $ forces.out > file.tmp")
    endif
    open(1,file='file.tmp')
    read(1,*) l
    close(1)
    call system("rm file.tmp")
    write(strl, '(I10)') (natoms+1)*steps
    if(l.eq.(natoms+1)*steps) then
!            print *, 'forces.out is ready for restart!'
    elseif(l.lt.(natoms+1)*steps) then
        print *, 'forces.out file is incomplete, edit restart file and relaunch'
        STOP;
    elseif(l.gt.(natoms+1)*steps) then
!            print *, 'Cutting exceding lines from forces.out'
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a7,i4.4,a4)") "forces_", id, ".out"
            call system('head -n '//trim(adjustl(strl))//' '//filename//' > forces.out_tmp')
            call system('mv forces.out_tmp '//filename)
        else
            call system('head -n '//trim(adjustl(strl))//' forces.out > forces.out_tmp')
            call system('mv forces.out_tmp forces.out')
        endif
!            print *, 'forces.out is ready for restart!'
    endif
    endif

    if(excN.gt.0) then
!Checking nact.out
        if(cs.gt.0) then
        if(lprint.gt.1) then
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a5,i4.4,a4)") "nact_", id, ".out"
            call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
        else
            call system("grep -c $ nact.out > file.tmp")
        endif
        open(1,file='file.tmp')
        read(1,*) l
        close(1)
        call system("rm file.tmp")
        write(strl, '(I10)') steps
        if(l.eq.steps) then
!            print *, 'nact.out is ready for restart!'
        elseif(l.lt.steps) then
            print *, 'nact.out file is incomplete, edit restart file and relaunch'
            STOP;
        elseif(l.gt.steps) then
!            print *, 'Cutting exceding lines from nact.out'
            if(trim(adjustl(dynam_type)).eq.'aimc') then
                write (filename, "(a5,i4.4,a4)") "nact_", id, ".out"
                call system('head -n '//trim(adjustl(strl))//' '//filename//' > nact.out_tmp')
                call system('mv nact.out_tmp '//filename)
            else
                call system('head -n '//trim(adjustl(strl))//' nact.out > nact.out_tmp')
                call system('mv nact.out_tmp nact.out')
            endif
!            print *, 'nact.out is ready for restart!'
        endif
        endif
        endif

!Checking pes.out file
        if(lprint.gt.1) then
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a4,i4.4,a4)") "pes_", id, ".out"
            call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
        else
            call system("grep -c $ pes.out > file.tmp")
        endif
        open(1,file='file.tmp')
        read(1,*) l
        close(1)
        call system("rm file.tmp")
        write(strl, '(I10)') steps
        if(l.eq.steps) then
!            print *, 'pes.out is ready for restart!'
        elseif(l.lt.steps) then
            print *, 'pes.out file is incomplete, edit restart file and relaunch'
            STOP;
        elseif(l.gt.steps) then
!            print *, 'Cutting exceding lines from pes.out'
            if(trim(adjustl(dynam_type)).eq.'aimc') then
                write (filename, "(a4,i4.4,a4)") "pes_", id, ".out"
                call system('head -n '//trim(adjustl(strl))//' '//filename//' > pes.out_tmp')
                call system('mv pes.out_tmp '//filename)
            else
                call system('head -n '//trim(adjustl(strl))//' pes.out > pes.out_tmp')
                call system('mv pes.out_tmp pes.out')
            endif
!            print *, 'pes.out is ready for restart!'
        endif
        endif

!Checking order.out file
        if(cs.gt.0) then
        if(lprint.gt.2) then
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a6,i4.4,a4)") "order_", id, ".out"
            call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
        else
            call system("grep -c $ order.out > file.tmp")
        endif
        open(1,file='file.tmp')
        read(1,*) l
        close(1)
        call system("rm file.tmp")
        write(strl, '(I10)') steps
        if(l.eq.steps) then
!            print *, 'order.out is ready for restart!'
        elseif(l.lt.steps) then
            print *, 'order.out file is incomplete, edit restart file and relaunch'
            STOP;
        elseif(l.gt.steps) then
!            print *, 'Cutting exceding lines from order.out'
            if(trim(adjustl(dynam_type)).eq.'aimc') then
                write (filename, "(a6,i4.4,a4)") "order_", id, ".out"
                call system('head -n '//trim(adjustl(strl))//' '//filename//' > order.out_tmp')
                call system('mv order.out_tmp '//filename)
            else
                call system('head -n '//trim(adjustl(strl))//' order.out > order.out_tmp')
                call system('mv order.out_tmp order.out')
            endif
!            print *, 'order.out is ready for restart!'
        endif
        endif
        endif

!Checking tdipole.out file
        if(printTdipole.eq.1.and.ired.eq.0) then
            if(trim(adjustl(dynam_type)).eq.'aimc') then
                write (filename, "(a8,i4.4,a4)") "tdipole_", id, ".out"
                call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
            else
                call system("grep -c $ tdipole.out > file.tmp")
            endif
            open(1,file='file.tmp')
            read(1,*) l
            close(1)
            call system("rm file.tmp")
            write(strl, '(I10)') steps*excN
            if(l.eq.steps*excN) then
!                print *, 'tdipole.out is ready for restart!'
            elseif(l.lt.steps*excN) then
                print *, 'tdipole.out file is incomplete, edit restart file and relaunch'
                STOP;
            elseif(l.gt.steps*excN) then
!                print *, 'Cutting exceding lines from tdipole.out'
                if(trim(adjustl(dynam_type)).eq.'aimc') then
                    write (filename, "(a8,i4.4,a4)") "tdipole_", id, ".out"
                    call system('head -n '//trim(adjustl(strl))//' '//filename//' > tdipole.out_tmp')
                    call system('mv tdipole.out_tmp '//filename)
                else
                    call system('head -n '//trim(adjustl(strl))//' tdipole.out > tdipole.out_tmp')
                    call system('mv tdipole.out_tmp tdipole.out')
                endif
!                print *, 'tdipole.out is ready for restart!'
            endif
        endif
    endif

!Checking temperature.out file
    if(lprint.gt.0) then
    if(trim(adjustl(dynam_type)).eq.'aimc') then
        write (filename, "(a12,i4.4,a4)") "temperature_", id, ".out"
        call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
    else
        call system("grep -c $ temperature.out > file.tmp")
    endif
    open(1,file='file.tmp')
    read(1,*) l
    close(1)
    call system("rm file.tmp")
    if(trim(adjustl(dynam_type)).eq.'aimc'.and.id.ne.0) then
        l1 = steps
    else
        l1 = steps+1
    endif
    write(strl, '(I10)') l1
    if(l.eq.l1) then
    elseif(l.lt.steps) then
        print *, 'temperature.out file is incomplete, edit restart file and relaunch'
        STOP;
    elseif(l.gt.l1) then
!            print *, 'Cutting exceding lines from temperature.out'
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a12,i4.4,a4)") "temperature_", id, ".out"
            call system('head -n '//trim(adjustl(strl))//' '//filename//' > temperature.out_tmp')
            call system('mv temperature.out_tmp '//filename)
        else
            call system('head -n '//trim(adjustl(strl))//' temperature.out > temperature.out_tmp')
            call system('mv temperature.out_tmp temperature.out')
        endif
!            print *, 'temperature.out is ready for restart!'
    endif
    endif

    if(excN.gt.0.and.lprint.gt.1) then
!Checking transition-densities.out file
        if(cs.gt.0) then
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a21,i4.4,a4)") "transition-densities_", id, ".out"
            call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
        else
            call system("grep -c $ transition-densities.out > file.tmp")
        endif
        open(1,file='file.tmp')
        read(1,*) l
        close(1)
        call system("rm file.tmp")
        if(trim(adjustl(dynam_type)).eq.'tully') then
            l2 = steps
        else
            l2 = steps*excN
        endif
        write(strl, '(I10)') l2
        if(l.eq.l2) then
!            print *, 'transition-densities.out is ready for restart!'
        elseif(l.lt.l2) then
            print *, 'transition-densities.out file is incomplete, edit restart file and relaunch'
            STOP;
        elseif(l.gt.l2) then
!            print *, 'Cutting exceding lines from transition-densities.out'
            if(trim(adjustl(dynam_type)).eq.'aimc') then
                write (filename, "(a21,i4.4,a4)") "transition-densities_", id, ".out"
                call system('head -n '//trim(adjustl(strl))//' '//filename//' > transition-densities.out_tmp')
                call system('mv transition-densities.out_tmp '//filename)
            else
                call system('head -n '//trim(adjustl(strl))//' transition-densities.out > transition-densities.out_tmp')
                call system('mv transition-densities.out_tmp transition-densities.out')
            endif
!            print *, 'transition-densities.out is ready for restart!'
        endif
        endif
    endif

!Checking velocity.out
    if(lprint.gt.0) then
    if(trim(adjustl(dynam_type)).eq.'aimc') then
        write (filename, "(a9,i4.4,a4)") "velocity_", id, ".out"
        call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
    else
        call system("grep -c $ velocity.out > file.tmp")
    endif
    open(1,file='file.tmp')
    read(1,*) l
    close(1)
    call system("rm file.tmp")
    write(strl, '(I10)') (natoms+3)*steps
    if(l.eq.(natoms+3)*steps) then
!            print *, 'velocity.out is ready for restart!'
    elseif(l.lt.(natoms+3)*steps) then
        print *, 'velocity.out file is incomplete, edit restart file and relaunch'
        STOP;
    elseif(l.gt.(natoms+3)*steps) then
!            print *, 'Cutting exceding lines from velocity.out'
        if(trim(adjustl(dynam_type)).eq.'aimc') then
            write (filename, "(a9,i4.4,a4)") "velocity_", id, ".out"
            call system('head -n '//trim(adjustl(strl))//' '//filename//' > velocity.out_tmp')
            call system('mv velocity.out_tmp '//filename)
        else
            call system('head -n '//trim(adjustl(strl))//' velocity.out > velocity.out_tmp')
            call system('mv velocity.out_tmp velocity.out')
        endif
!            print *, 'velocity.out is ready for restart!'
    endif
    endif
    
    if(excN.gt.0) then
        if(trim(adjustl(dynam_type)).eq.'aimc'.or.trim(adjustl(dynam_type)).eq.'mf') then
!Checking nacr.out file
            if(cs.gt.0) then
            if(lprint.gt.1) then
            if(trim(adjustl(dynam_type)).eq.'aimc') then
                write (filename, "(a5,i4.4,a4)") "nacr_", id, ".out"
                call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
            else
                call system("grep -c $ nacr.out > file.tmp")
            endif
            open(1,file='file.tmp')
            read(1,*) l
            close(1)
            call system("rm file.tmp")
            write(strl, '(I10)') (excN*(excN-1))/2*steps
            if(l.eq.(excN*(excN-1))/2*steps) then
!                print *, 'nacr.out is ready for restart!'
            elseif(l.lt.(excN*(excN-1))/2*steps) then
                print *, 'nacr.out file is incomplete, edit restart file and relaunch'
                STOP;
            elseif(l.gt.(excN*(excN-1))/2*steps) then
!                print *, 'Cutting exceding lines from nacr.out'
                if(trim(adjustl(dynam_type)).eq.'aimc') then
                    write (filename, "(a5,i4.4,a4)") "nacr_", id, ".out"
                    call system('head -n '//trim(adjustl(strl))//' '//filename//' > nacr.out_tmp')
                    call system('mv nacr.out_tmp '//filename)
                else
                    call system('head -n '//trim(adjustl(strl))//' nacr.out > nacr.out_tmp')
                    call system('mv nacr.out_tmp nacr.out')
                endif
!                print *, 'nacr.out is ready for restart!'
            endif
            endif
!Checking state_forces.out file
            if(trim(adjustl(dynam_type)).eq.'aimc') then
                write (filename, "(a13,i4.4,a4)") "state_forces_", id, ".out"
                call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
            else
                call system("grep -c $ state_forces.out > file.tmp")
            endif
            open(1,file='file.tmp')
            read(1,*) l
            close(1)
            call system("rm file.tmp")
            write(strl, '(I10)') (natoms*excN+1)*steps
            if(l.eq.(natoms*excN+1)*steps) then
!                print *, 'state_forces.out is ready for restart!'
            elseif(l.lt.(natoms*excN+1)*steps) then
                print *, 'state_forces.out file is incomplete, edit restart file and relaunch'
                STOP;
            elseif(l.gt.(natoms*excN+1)*steps) then
!                print *, 'Cutting exceding lines from state_forces.out'
                if(trim(adjustl(dynam_type)).eq.'aimc') then
                    write (filename, "(a13,i4.4,a4)") "state_forces_", id, ".out"
                    call system('head -n '//trim(adjustl(strl))//' '//filename//' > state_forces.out_tmp')
                    call system('mv state_forces.out_tmp '//filename)
                else
                    call system('head -n '//trim(adjustl(strl))//' state_forces.out > state_forces.out_tmp')
                    call system('mv state_forces.out_tmp state_forces.out')
                endif
!                print *, 'state_forces.out is ready for restart!'
            endif
            endif
        endif
    endif

!Checking MCE files:
    if(trim(adjustl(dynam_type)).eq.'aimc'.and.id0.eq.0) then
!Checking pop.dat file
        write (filename, "(a7)") 'pop.dat'
        call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
        open(1,file='file.tmp')
        read(1,*) l
        close(1)
        call system("rm file.tmp")
        write(strl, '(I10)') aimcSteps
        if(l.eq.aimcSteps) then
!           print *, 'pop.dat is ready for restart!'
        elseif(l.lt.aimcSteps) then
            print *, 'pop.dat file is incomplete, edit restart file and relaunch'
            STOP;
        elseif(l.gt.aimcSteps) then
!           print *, 'Cutting exceding lines from pop.dat'
            call system('head -n '//trim(adjustl(strl))//' '//filename//' > pop.dat_tmp')
            call system('mv pop.dat_tmp '//filename)
!           print *, 'pop.dat is ready for restart!'
        endif
!Checking nuclear_coeff.dat file
        write (filename, "(a17)") 'nuclear_coeff.dat'
        call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
        open(1,file='file.tmp')
        read(1,*) l
        close(1)
        call system("rm file.tmp")
        write(strl, '(I10)') aimcSteps
        if(l.eq.aimcSteps) then
!           print *, 'nuclear_coeff.dat is ready for restart!'
        elseif(l.lt.aimcSteps) then
            print *, 'nuclear_coeff.dat file is incomplete, edit restart file and relaunch'
            STOP;
        elseif(l.gt.aimcSteps) then
!           print *, 'Cutting exceding lines from nuclear_coeff.dat'
            call system('head -n '//trim(adjustl(strl))//' '//filename//' > nuclear_coeff.dat_tmp')
            call system('mv nuclear_coeff.dat_tmp '//filename)
           print *, 'nuclear_coeff.dat is ready for restart!'
        endif
!Checking electronic_overlaps.dat file
        call get_all_clone_steps(clone_steps, Nsim, dps, dt, outData, outDataCoords)
        call get_all_dropouts(dropout_steps,dps,dt)
        call get_electronic_overlaps_n_lines(clone_steps, dropout_steps, Nsim, dps, aimcSteps, excN, nlines)
        write (filename, "(a24)") 'electronic_overlaps.dat'
        call system("grep -c $ "//trim(adjustl(filename))//" > file.tmp")
        open(1,file='file.tmp')
        read(1,*) l
          close(1)
        call system("rm file.tmp")
        write(strl, '(I10)') nlines
        if(l.eq.nlines) then
            print *, 'electronic_overlaps.dat is ready for restart!'
        elseif(l.lt.nlines) then
            print *, 'electronic_overlaps.dat file is incomplete, edit restart file and relaunch'
            STOP;
        elseif(l.gt.nlines) then
            print *, 'Cutting exceding lines from electronic_overlaps.dat'
            call system('head -n '//trim(adjustl(strl))//' '//filename//' > electronic_overlaps.dat_tmp')
            call system('mv electronic_overlaps.dat_tmp '//filename)
            print *, 'nuclear_coeff.dat is ready for restart!'
        endif
    endif

endsubroutine check_files

!----------------------------------------------------------------!
!This subroutine calculates the number of steps before each clone!
!----------------------------------------------------------------!
subroutine get_all_clone_steps(clone_steps, Nsim, Ndps, dt, outData, outDataCoords)
    
    implicit none

    integer, intent(out), dimension(Nsim + Ndps) :: clone_steps
    integer, intent(in) :: Nsim
    integer, intent(in) :: Ndps
    _REAL_, intent(in) :: dt
    integer, intent(in) :: outData
    integer, intent(in) :: outDataCoords

    integer :: i, j
    character*150 filename
    _REAL_ :: time
    
    clone_steps = 1
    do i=1,Nsim-1+Ndps
        filename = ''
        write (filename, "(a8,i4.4,a4)") "coeff-n_", i, ".out"
        open(1,file=filename)
            read(1,*) j, time
            clone_steps(i+1) = nint(time/dt/outData*outDataCoords) + 1
        close(1)
    enddo

end subroutine get_all_clone_steps

!----------------------------------------------------!
!This subroutine calculates the step for each dropout!
!----------------------------------------------------!
subroutine get_all_dropouts(dropout_steps, Ndps, dt)

    implicit none

    integer, intent(out), dimension(Ndps) :: dropout_steps
    integer, intent(in) :: Ndps
    _REAL_, intent(in) :: dt
    logical file_exists2
    integer i
    _REAL_ :: dropout_time

    file_exists2 = .false.

    inquire(FILE="dropped.out", EXIST=file_exists2)
    if(file_exists2) then
        open(1, file='dropped.out')
        do i = 1, Ndps
            read(1,*) dropout_time
            dropout_steps(i) = dropout_time/dt + 1
        enddo
        close(1)
    endif
    
end subroutine get_all_dropouts

!------------------------------------------------------------------------------!
!This subroutine calculates the number of lines in the 'electronic_overlap.dat'!
!file                                                                          !
!------------------------------------------------------------------------------!
subroutine get_electronic_overlaps_n_lines(clone_steps,dropout_steps,NsimLast,Ndps,nsteps,excN,nlines)

    implicit none

    integer, dimension(NsimLast + Ndps), intent(in) :: clone_steps
    integer, dimension(Ndps), intent(in) :: dropout_steps
    integer, intent(in) :: NsimLast
    integer, intent(in) :: Ndps
    integer, intent(in) :: nsteps
    integer, intent(in) :: excN
    integer, intent(out) :: nlines

    integer i,j,k
    integer nsim

    nlines = 0
    nsim = 1
    j = 2
    k = 1

    do i = 1, nsteps
        if(i.eq.clone_steps(j)) then
            j = j + 1
            if(j.eq.NsimLast + Ndps + 1) j = 2
            nsim = nsim + 1
        endif
        if(i.eq.dropout_steps(k) + 1) then
            k = k + 1
            if (k.eq.Ndps + 1) k = 1
            nsim = nsim - 1
        endif
        nlines = nlines + (excN**2 * nsim * (nsim-1)) / 2
    enddo

end subroutine get_electronic_overlaps_n_lines

!--------------------------------------------------------------------!
!This subroutine reads the number of classical steps in an input.ceon!
!--------------------------------------------------------------------!
subroutine read_input_moldyn_block(file_name,nstep,dynam_type,excN,cs,outData,outDataCoords,natom,p,dt,lprint,nq,nclones00,ired)

    implicit none

!Passed in
    character*150, intent(in) :: file_name
    integer, intent(out) :: nstep
    character*150, intent(out) :: dynam_type
    integer, intent(out) :: excN
    integer, intent(out) :: cs
    integer, intent(out) :: outData
    integer, intent(out) :: outDataCoords
    integer, intent(out) :: natom
    integer, intent(out) :: p
    _REAL_, intent(out) ::  dt
    integer, intent(out) :: lprint
    integer, intent(out) :: nq
    integer, intent(out) :: nclones00
    integer, intent(out) :: ired

!Variables of the moldyn namelist
    integer natoms
    integer bo_dynamics_flag,exc_state_init,n_exc_states_propagate
    integer out_count_init
    _REAL_ time_init,time_step
    integer n_class_steps,n_quant_steps,quant_coeffs_reinit
    _REAL_ num_deriv_step,therm_temperature
    integer therm_type
    _REAL_ berendsen_relax_const
    integer heating,heating_steps_per_degree,out_data_steps
    integer out_coords_steps
    _REAL_ therm_friction
    integer rnd_seed,out_data_cube,verbosity,moldyn_deriv_flag
    _REAL_ quant_step_reduction_factor
    _REAL_ decoher_e0,decoher_c
    integer decoher_type,dotrivial
    _REAL_ deltared
    integer iredpot,nstates
    integer ifixed
    integer :: nmc, npc
    integer :: printTdipole
    _REAL_ :: AIMC_dclone_1,AIMC_dclone_2,AIMC_dclone_3,AIMC_max_clone,AIMC_force_pop_min
    integer nclones0
    _REAL_ :: dpt 
    character(100) :: NAMD_type
    namelist /moldyn/ natoms,bo_dynamics_flag,NAMD_type,exc_state_init, &
        n_exc_states_propagate,out_count_init,time_init, &
        time_step,n_class_steps,n_quant_steps, &
        num_deriv_step, &
        therm_temperature,therm_type, &
        berendsen_relax_const,heating, &
        heating_steps_per_degree,out_data_steps,out_coords_steps, &
        therm_friction,rnd_seed,out_data_cube,verbosity,moldyn_deriv_flag, &
        quant_step_reduction_factor,decoher_e0,decoher_c,decoher_type,dotrivial,&
        iredpot,nstates,deltared,ifixed,AIMC_dclone_1,AIMC_dclone_2,AIMC_dclone_3,&
        AIMC_max_clone,AIMC_force_pop_min,nmc,npc,printTdipole,nclones0,dpt

!Default values
    n_class_steps = 0
    NAMD_type='tully'
    n_exc_states_propagate = 0
    exc_state_init = 0
    out_coords_steps = 1
    out_data_steps = 1
    printTdipole = 0
    time_step = 0.1d0
    verbosity = -1
    n_quant_steps = 1
    nclones0=0
    ired=0

!Reading &moldyn
    open(1,file=trim(adjustl(file_name)),status='old')
    read(1,nml=moldyn)
    close(1)
    nstep=n_class_steps
    dynam_type=NAMD_type
    excN=n_exc_states_propagate
    cs=exc_state_init
    outData=out_coords_steps*out_data_steps
    outDataCoords=out_coords_steps
    natom=natoms
    p=printTdipole
    dt=time_step
    lprint=verbosity
    nq=n_quant_steps
    nclones00=nclones0
    ired=iredpot

endsubroutine read_input_moldyn_block

!----------------------------------------------------------!
!This subroutine reads the initial data for restarting aimc!
!----------------------------------------------------------!
subroutine read_aimc_for_restart(nuclear,Nsim,excN)

    use communism

    implicit none

    type(mce) :: nuclear
    integer, intent(in) :: Nsim
    integer, intent(in) :: excN

    integer i,j,n,m,k
    _REAL_ time
    _REAL_, dimension(2*Nsim) :: D_aux
    character*10 strl
    _REAL_, dimension(Nsim,2*Nsim) :: Heff_aux

!Reading initial values for the nuclear coefficients
    call system("tail -1 nuclear_coeff.dat > file.tmp")
    open(1,file='file.tmp')
        read(1,*) i, time, D_aux(1:2*Nsim)
    close(1)
    do i = 1, Nsim
        nuclear%D(i) = D_aux(2*i-1)+(0.0d0,1.0d0)*D_aux(2*i)
    enddo
!Reading initial values for the electronic overlaps
    nuclear%sE(1:Nsim,1:Nsim,1:excN,1:excN) = 0.0d0
    write(strl, '(I10)') (excN**2 * Nsim * (Nsim-1)) / 2
    call system('tail -'//trim(adjustl(strl))//' electronic_overlaps.dat > file.tmp')
    open(1,file='file.tmp')
        do n = 1, Nsim
            do m = n + 1, Nsim
                do i = 1, excN
                    do j = 1, excN
                        read(1,*) k, time, k, k, k, k, nuclear%sE(n,m,i,j)
                        nuclear%sE(m,n,j,i) = nuclear%sE(n,m,i,j)
                    enddo
                enddo
            enddo
        enddo
    close(1)
    do n = 1, Nsim
        do i = 1, excN
            nuclear%sE(n,n,i,i) = 1.0d0
        enddo
    enddo
    call system('rm file.tmp')
!Reading the initial Heff
    open(1,file='Heff_last.out')
        read(1,*)
        do i = 1, Nsim
            read(1,*) Heff_aux(i,:)
        enddo
    close(1)
    do i = 1, Nsim
        do j = 1, Nsim
            nuclear%Heff(i,j) = Heff_aux(i,2*j-1) + (0.0d0,1.0d0) * Heff_aux(i,2*j)
        enddo
    enddo

end subroutine read_aimc_for_restart

!---------------------------------------------------------!
!This subroutines check syncrhronization for AIMC dynamics!
!---------------------------------------------------------!
subroutine check_synchronization(Nsim)

    implicit none

    integer, intent(in) :: Nsim
    
    integer i,j,k
    character*150 filename
    integer, dimension(Nsim) :: steps
    integer step
    integer maxSteps
    character*150 dynam_type
    integer excN, cs, outData, outDataCoords,natoms,printTdipole
    _REAL_ dt
    character*150, dimension(200) :: b1
    integer b1lines
    integer, allocatable, dimension(:) :: atoms
    _REAL_, allocatable, dimension(:,:) :: coords
    _REAL_, allocatable, dimension(:,:) :: veloc
    _REAL_, allocatable, dimension(:,:) :: coeff
    integer lprint
    integer nq
    integer nclones0
    integer ired
    integer dps
    _REAL_ junk_time
    integer, allocatable :: indexes(:)
    logical file_exists2

    dps = 0
    inquire(FILE="dropped.out", EXIST=file_exists2)
    if(file_exists2) then
        call system('wc dropped.out > dropped.tmp')
        open(1,file='dropped.tmp')
            read(1,*) dps
        close(1)
        call system('rm dropped.tmp')
        if(dps.gt.0) then
            allocate(indexes(dps))
            open(1,file='dropped.out')
                do i = 1, dps
                    read(1,*) junk_time, indexes(i)
                enddo
            close(1)
        else
            allocate(indexes(1))
            indexes(1) = -1
        endif
    endif

    call read_input_moldyn_block('input.ceon',maxSteps,dynam_type,excN,cs,&
                    outData,outDataCoords,natoms,printTdipole,dt,lprint,nq,nclones0,ired)
    steps = 0
    j = 0
    do i = 1, Nsim + dps
        if(.not.any(i - 1 == indexes)) then
            j = j + 1
            write (filename, "(a8,i4.4,a4)") "restart_", i-1, ".out"
            call read_input_moldyn_block(filename,step,dynam_type,excN,cs,&
                outData,outDataCoords,natoms,printTdipole,dt,lprint,nq,nclones0,ired)
            steps(j) = step
        endif
    enddo
    j = 0
    if(minval(steps).ne.maxval(steps)) then
        allocate(atoms(natoms))
        allocate(coords(natoms,3),veloc(natoms,3))
        allocate(coeff(excN,2))
        do i = 1, Nsim + dps
            if(file_exists2) then
                if(.not.any(i - 1 == indexes)) then
                    j = j + 1
                    if(maxval(steps).ne.steps(i)) then
                        print *, 'Remaking restart for tr: ', i - 1
                        write (filename, "(a8,i4.4,a4)") "restart_", i-1, ".out"
                        call get_b1(filename,b1,b1lines)
                        call read_input_moldyn_block(filename,steps(i),dynam_type,excN,cs,&
                            outData,outDataCoords,natoms,printTdipole,dt,lprint,nq,nclones0,ired)
                        write (filename, "(a7,i4.4,a4)") "coords_", i-1, ".xyz"
                        call get_coords(filename,maxSteps-maxval(steps),atoms,coords,natoms,dt)
                        write (filename, "(a9,i4.4,a4)") "velocity_", i-1, ".out"
                        call get_veloc(filename,maxSteps-maxval(steps),veloc,natoms,dt)
                        write (filename, "(a12,i4.4,a4)") "coefficient_", i-1, ".out"
                        call get_coeff(filename,maxSteps-maxval(steps),coeff,excN,dt)
                        write (filename, "(a8,i4.4,a4)") "restart_", i-1, ".out"
                        call write_new_restart(filename,b1,b1lines,atoms,coords,veloc,natoms,coeff,excN,maxval(steps),maxSteps-maxval(steps),dt,nclones0)
                    endif
                endif
            else
                if(maxval(steps).ne.steps(i)) then
                    print *, 'Remaking restart for tr: ', i - 1
                    write (filename, "(a8,i4.4,a4)") "restart_", i-1, ".out"
                    call get_b1(filename,b1,b1lines)
                    call read_input_moldyn_block(filename,steps(i),dynam_type,excN,cs,&
                        outData,outDataCoords,natoms,printTdipole,dt,lprint,nq,nclones0,ired)
                    write (filename, "(a7,i4.4,a4)") "coords_", i-1, ".xyz"
                    call get_coords(filename,maxSteps-maxval(steps),atoms,coords,natoms,dt)
                    write (filename, "(a9,i4.4,a4)") "velocity_", i-1, ".out" 
                    call get_veloc(filename,maxSteps-maxval(steps),veloc,natoms,dt)
                    write (filename, "(a12,i4.4,a4)") "coefficient_", i-1, ".out"
                    call get_coeff(filename,maxSteps-maxval(steps),coeff,excN,dt)
                    write (filename, "(a8,i4.4,a4)") "restart_", i-1, ".out"
                    call write_new_restart(filename,b1,b1lines,atoms,coords,veloc,natoms,coeff,excN,maxval(steps),maxSteps-maxval(steps),dt,nclones0)
                endif
            endif
        enddo
    endif
end subroutine check_synchronization

!---------------------------------------------------------------!
!This subroutine gets the first block of a restart_????.out file!
!---------------------------------------------------------------!
subroutine get_b1(filename,b1,b1lines)

    implicit none

    character*150, intent(in) :: filename
    character*150, dimension(200), intent(out) :: b1
    integer, intent(out) :: b1lines

    integer ios

    b1lines = 1
    open(1,file=trim(adjustl(filename)))
        do
            read(1, '(A)', iostat=ios) b1(b1lines)
            if (ios /= 0 .or. b1(b1lines)(2:6).eq.'coord') exit
            b1lines = b1lines + 1 
        enddo
    close(1)

end subroutine get_b1

!----------------------------------------------------------------------------!
!This subroutine gets the coords from a coords_????.xyz file for a given step!
!----------------------------------------------------------------------------!
subroutine get_coords(filename,step,atoms,coords,natoms,time_step)

    implicit none

    character*150, intent(in) :: filename
    integer, intent(in) :: step
    integer, intent(out), dimension(natoms) :: atoms
    _REAL_, intent(out), dimension(natoms,3) :: coords
    integer, intent(in) :: natoms
    _REAL_, intent(in) :: time_step

    integer i, j, ios
    character*28 str1
    _REAL_ time
    character*9 str2
    _REAL_ energy
    character*3 element
    logical found

    found = .false.
    open(1,file=filename)
    do
        read(1,*,iostat=ios)
        if (ios /= 0) then
            print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
            STOP;
        endif
        read(1,449,iostat=ios) str1, energy, str2, time
        if (ios /= 0) then
            print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
            STOP;
        endif
!        print *, time, step*time_step, nint(time*100000).eq.nint(step*time_step*100000)
        if (nint(time*100000).eq.nint(step*time_step*100000)) then
            found = .true.
            do i = 1, natoms
                read(1,302,iostat=ios) element, (coords(i,j), j = 1, 3)
                if (ios /= 0) then
                    print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
                endif
                if(trim(adjustl(element)).eq.'H') then
                    atoms(i) = 1
                elseif(trim(adjustl(element)).eq.'He') then
                    atoms(i) = 2
                elseif(trim(adjustl(element)).eq.'Li') then
                    atoms(i) = 3
                elseif(trim(adjustl(element)).eq.'Be') then
                    atoms(i) = 4
                elseif(trim(adjustl(element)).eq.'B') then
                    atoms(i) = 5
                elseif(trim(adjustl(element)).eq.'C') then
                    atoms(i) = 6
                elseif(trim(adjustl(element)).eq.'N') then
                    atoms(i) = 7
                elseif(trim(adjustl(element)).eq.'O') then
                    atoms(i) = 8
                elseif(trim(adjustl(element)).eq.'F') then
                    atoms(i) = 9
                elseif(trim(adjustl(element)).eq.'Ne') then
                    atoms(i) = 10
                elseif(trim(adjustl(element)).eq.'Si') then
                    atoms(i) = 14
                elseif(trim(adjustl(element)).eq.'P') then
                    atoms(i) = 15
                elseif(trim(adjustl(element)).eq.'S') then
                    atoms(i) = 16
                elseif(trim(adjustl(element)).eq.'Cl') then
                    atoms(i) = 17
                elseif(trim(adjustl(element)).eq.'Zn') then
                    atoms(i) = 30
                else
                    print *, 'Fatal error: element not supported for remaking restart:', i, trim(adjustl(element))
                    STOP;
                endif
            enddo
        else
            do i = 1, natoms
                read(1,*,iostat=ios)
                if (ios /= 0) then
                    print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
                endif
            enddo
        endif
        if(found) then 
            close(1)
            exit
        endif
    enddo

449     format(a28,F18.10,a9,F18.10)
302     format(A3,3f16.10)

end subroutine get_coords

!----------------------------------------------------------------------------------!
!This subroutine gets the velocities from a velocity_????.out file for a given step!
!----------------------------------------------------------------------------------!
subroutine get_veloc(filename,step,veloc,natoms,time_step)

    implicit none

    character*150, intent(in) :: filename
    integer, intent(in) :: step
    _REAL_, intent(out), dimension(natoms,3) :: veloc
    integer, intent(in) :: natoms
    _REAL_, intent(in) :: time_step

    integer i, j, ios
    character*28 str1
    _REAL_ time
    character*9 str2
    _REAL_ energy
    logical found

    found = .false.
    open(1,file=filename)
    do
        read(1,449,iostat=ios) str1, energy, str2, time
        if (ios /= 0) then
            print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
            STOP;
        endif
        read(1,*,iostat=ios)
        if (ios /= 0) then
            print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
            STOP;
        endif
!        print *, time, step*time_step, nint(time*100000).eq.nint(step*time_step*100000)
        if (nint(time*100000).eq.nint(step*time_step*100000)) then
            found = .true.
            do i = 1, natoms
                read(1,223,iostat=ios) (veloc(i,j), j = 1, 3)
                if (ios /= 0) then
                    print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
                endif
            enddo
            read(1,*,iostat=ios)
            if (ios /= 0) then
                print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
                STOP;
            endif
        else
            do i = 1, natoms
                read(1,*,iostat=ios)
                if (ios /= 0) then
                    print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
                endif
            enddo
            read(1,*,iostat=ios)
            if (ios /= 0) then
                print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
                STOP;
            endif
        endif
        if(found) then
            close(1)
            exit
        endif
    enddo

449     format(a28,F18.10,a9,F18.10)
223     FORMAT(3(1x,F16.10))

end subroutine get_veloc

!--------------------------------------------------------------------------------!
!This subroutine gets the coeff from a coefficient_????.out file for a given step!
!--------------------------------------------------------------------------------!
subroutine get_coeff(filename,step,coeff,excN,time_step)

    implicit none 

    character*150, intent(in) :: filename
    integer, intent(in) :: step 
    integer, intent(in) :: excN
    _REAL_, intent(out), dimension(excN,2) :: coeff
    _REAL_, intent(in) :: time_step

    integer i, j, ios
    character*28 str1 
    _REAL_ time 
    character*9 str2 
    _REAL_ energy
    logical found
    _REAL_, dimension(5) :: act

    found = .false.
    open(1,file=trim(adjustl(filename)))
    do   
        read(1,449,iostat=ios) str1, energy, str2, time 
        if (ios /= 0) then 
            print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
            STOP;
        endif
        read(1,*,iostat=ios)
        if (ios /= 0) then 
            print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
            STOP;
        endif
        if (nint(time*100000).eq.nint(step*time_step*100000)) then 
            found = .true.
            do i = 1, excN
                read(1,223,iostat=ios) (act(j), j = 1, 5)
                if (ios /= 0) then 
                    print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
                endif
                coeff(i,1) = sqrt(act(1)**2 + act(2)**2)
                coeff(i,2) = atan2(act(2),act(1))
            enddo
            read(1,*,iostat=ios)
            if (ios /= 0) then 
                print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
                STOP;
            endif
        else 
            do i = 1, excN
                read(1,*,iostat=ios)
                if (ios /= 0) then 
                    print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
                endif
            enddo
            read(1,*,iostat=ios)
            if (ios /= 0) then 
                print *, 'Fatal error, file ', trim(adjustl(filename)), 'incomplete, unable to remake the corresponding restart'
                STOP;
            endif
        endif
        if(found) then 
            close(1)
            exit 
        endif
    enddo

449     format(a28,F18.10,a9,F18.10)
223     FORMAT(5(1x,F16.10))

end subroutine get_coeff

!-------------------------------------------------------------------!
!This subroutines write new restart files to achieve synchronization!
!-------------------------------------------------------------------!
subroutine write_new_restart(filename,b1,b1lines,atoms,coords,veloc,natoms,coeff,excN,steps1,steps2,dt,nclones0)

    implicit none

    character*150, intent(in) :: filename
    character*150, dimension(200), intent(in) :: b1
    integer, intent(in) :: b1lines
    integer, intent(in) :: natoms
    integer, intent(in), dimension(natoms) :: atoms
    _REAL_, intent(in), dimension(natoms,3) :: coords
    _REAL_, intent(in), dimension(natoms,3) :: veloc
    integer, intent(in) :: excN
    _REAL_, intent(in), dimension(excN,2) :: coeff
    integer, intent(in) :: steps1
    integer, intent(in) :: steps2
    _REAL_, intent(in) :: dt
    integer, intent(in) :: nclones0

    integer i,j,k

    open(1,file=filename,action='write',status='replace')
    do i=1, b1lines
        if(b1(i)(1:13).eq.'   time_init=') then
            write(1,'(A13,F9.2,A26)') '   time_init=',dt*steps2,', ! Initial time, fs [0.0]'
        else if(b1(i)(1:17).eq.'   n_class_steps=') then
            write(1,'(A17,I7,A33)') '  n_class_steps=', steps1,', ! Number of classical steps [1]' 
        else if(b1(i)(1:12).eq.'   nclones0=') then
            write(1,'(A12,I2,A76)') '   nclones0=', nclones0, ", ! Clones count for 'aimc' (must be declared here for restarting 'aimc') [0]"
        else 
            write(1,'(A)') b1(i)(1:len_trim(b1(i)))
        endif
    enddo
    do i = 1, natoms
        write(1,'(I4,3F16.10)') atoms(i), (coords(i,j),j=1,3)
    enddo
    write(1,'(A)') '&endcoord'
    write(1,'(A)') '&veloc'
    do i = 1, natoms
        write(1,'(3F16.10)') (veloc(i,j),j=1,3)
    enddo
    write(1,'(A)') '&endveloc'
    write(1,'(A)') '&coeff'
    do i = 1, excN
        write(1,'(2F16.10)') coeff(i,1), coeff(i,2)
    enddo
    write(1,'(A)') '&endcoeff'
    close(1)

end subroutine write_new_restart

