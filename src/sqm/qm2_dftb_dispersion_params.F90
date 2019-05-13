! <compile=optimized>
!  -*- mode: f90; coding: iso-8859-15; -*-

#include "dprec.fh"
subroutine dispersion_params(qmmm_struct, nn,izp)

  use qm2_dftb_module, only: dispertmp, dispfile
  use qmmm_struct_module, only : qmmm_struct_type

  implicit none

!! Passed in:
  type(qmmm_struct_type), intent(in) :: qmmm_struct
  integer, intent(in) :: nn
  integer, intent(in) :: izp(*) ! izp(NNDIM)
! Local
  integer :: i,j,n
  _REAL_  :: r2

! if we read from DISPERSION.INP:
! determine hybridization from number of Hydrogens
! wir starten mit 1 Nachbarn. Dann zaehlen wir die
! Anzahl der Wasserstoffe durch: fuer jeden Wasserstoff gibt es
! einen Nachbarn mehr. Die werte fuer C und N unterscheiden sich nur in den
! Waserstoffen. N mit 2H ist anders als nur N oder N mit 1H
! C ist nur fuer 3H unterschiedlich

  open(15,file='DISP.CHEK')
  rewind(15)


! parameters as in Elstner et al., JCP 2000
    dispertmp%A=7
    dispertmp%B=4
    dispertmp%C=-3.0d0
    dispertmp%rv=0.0d0
    dispertmp%r0=0.0d0
    dispfile%scale=1.0d0

    ! determine parameters for all atoms
    ! If the "DISP.INP" file is present, it will already contain the hh1, hh2 and Ni parameters for every atom.
    ! If it is not present, then those values need to be calculated from data in DISPERSION.INP
    if (.not.dispfile%read_DISP_DOT_INP)  then 
       do i=1,nn
          dispfile%nei(i) = 1
          do j=1,nn
             if ( j /= i) then
                r2 = (qmmm_struct%qm_coords(1,i) - qmmm_struct%qm_coords(1,j))**2 &
                    +(qmmm_struct%qm_coords(2,i) - qmmm_struct%qm_coords(2,j))**2 &
                    +(qmmm_struct%qm_coords(3,i) - qmmm_struct%qm_coords(3,j))**2 
                if (r2 <= 1.2d0**2) then
                   dispfile%nei(i) = dispfile%nei(i) + 1
                endif
             endif
          enddo !j=1,nn
          dispfile%hh1(i)=dispfile%h1(izp(i),dispfile%nei(i))
          dispfile%hh2(i)=dispfile%h2(izp(i),dispfile%nei(i))
          dispfile%Ni(i) =dispfile%Ni0(izp(i))
          write(15,'(3F12.6)')  dispfile%hh1(i),dispfile%hh2(i),dispfile%Ni(i)
       end do
    endif

    do i=1,nn
       ! check values
       if ((dispfile%hh1(i) == 0.0d0) .OR. (dispfile%hh2(i) == 0.0d0) .OR. (dispfile%Ni(i) == 0.0d0)) then
          ! The calculation must stop!!
          write(6,'(10X,"*******************************************")')
          write(6,'(10X," WARNING: a parameter is 0.0 for atom ",i3)')i
          write(6,'(10X,"          IZP = ",i3)') izp(i)
          write(6,'(10X,"          nei = ",i3," (",i2," neighbors.)")') dispfile%nei(i), dispfile%nei(i)-1
          write(6,'(10X,"          hh1 = ",f6.2)') dispfile%hh1(i)
          write(6,'(10X,"          hh2 = ",f6.2)') dispfile%hh2(i)
          write(6,'(10X,"          Ni  = ",f6.2)') dispfile%Ni(i)
          write(6,'(10X,"*******************************************")')
          call sander_bomb("<dispersion_params>","qm2_dftb_dispersion_params.f","Exiting.")
       endif
           
    enddo !i=1,nn

! set up mixed coefficients
! mixing from Halgren JACS 1992 114 p7827
    if ( .NOT. dispfile%read_DISP_DOT_INP)  then
        write(15,*) ' --------------'
        write(15,*) ' I J  typeI typeJ C6 R NeiI NeiJ'
    endif


    do i=1,nn
        do j=1,i
            if (dispfile%scale <= 0.0d0) then
                dispertmp%C6(i,j) = -dispfile%scale*1.5d0*dispfile%hh1(i)*dispfile%hh1(j)* &
                dispfile%hh2(i)*dispfile%hh2(j)/ &
                ( dispfile%hh2(i)+dispfile%hh2(j) )
            else
            ! cc  17.532 conversion from eV in a.u. for polarizability
            ! 0.5975 conversion from [H au**6] to [eV A**6]
            ! total 17.532*0.5975 = 10.476
            ! * 1.5d0 = 15.714
               dispertmp%C6(i,j) = dispfile%scale*15.714d0*dispfile%hh1(i)*dispfile%hh1(j)/ &
                    ( sqrt(dispfile%hh1(i)/dispfile%Ni(i) ) + &
                    sqrt( dispfile%hh1(j)/dispfile%Ni(j) ) )
                dispertmp%Rvdw(i,j)=(dispfile%hh2(i)**3 + dispfile%hh2(j)**3)/ &
                     (dispfile%hh2(i)**2 + dispfile%hh2(j)**2 )
                dispertmp%Rvdw(j,i) =  dispertmp%Rvdw(i,j)
            endif
            dispertmp%C6(j,i) = dispertmp%C6(i,j)
            if ( .NOT. dispfile%read_DISP_DOT_INP)  then
                write(15,'(4I4,2F12.6,2I4)') i,j,izp(i),izp(j), &
                dispertmp%C6(i,j),dispertmp%Rvdw(i,j),dispfile%nei(i),dispfile%nei(j)
            endif
        enddo
    enddo

    close(15)
end subroutine dispersion_params
