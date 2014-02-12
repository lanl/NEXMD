!Standardized dx file writting routine.
!Prints out the data contained in voldata to file 'filename'
!voldata is assumed to be a grid of values of shape: voldata(1:xm,1:ym,1:zm)
!xm, ym, zm are number of nodes along x,y,z
!gox, goy, goz are grid origin in x,y,z, h is grid spacing (assumed
!uniform for now)
#  define _REAL_ double precision
subroutine gen_dx_file(xm,ym,zm,h,gox,goy,goz,voldata,filename,fn,dataname)

    implicit none
    !passed variables
    integer xm,ym,zm,fn
    _REAL_  h, gox, goy, goz
    _REAL_  voldata(1:xm,1:ym,1:zm)
    character*(*), intent(in) :: filename, dataname

    !internal variables
    integer cnt, i, j, k   
     
    !          write(6,*) '-debug:writing volumetric map in dx format';flush(6)
                open(fn,file=filename)
                write(fn,100) xm,ym,zm
                write(fn,110) gox+h,goy+h,goz+h
                write(fn,120) sngl(h),0,0
                write(fn,121) 0,sngl(h),0
                write(fn,122) 0,0,sngl(h)
                write(fn,130) xm,ym,zm
                write(fn,140) xm*ym*zm
                cnt=0
                do i=1,xm; do j=1,ym; 
                    do k=1,zm
                        cnt = cnt+1 
                        write(fn,150,advance='no') voldata(i,j,k)
                        if (cnt >= 3) then
                            cnt = 0
                            write(fn,*)''
                        end if
                    end do
                end do; end do
                write(fn,*) 'attribute "dep" string "positions"'
                write(fn,*) 'object "',dataname,'" class field'
                write(fn,*) 'component "positions" value 1'
                write(fn,*) 'component "connections" value 2'
                write(fn,*) 'component "data" value 3'
                close(fn)
    !           write(6,*) "-debug: done writting volumetric dx file";flush(6)
            100 format('object 1 class gridpositions counts',3I5)
            110 format('origin',3F9.4)
            120 format('delta',F10.7,2I2)
            121 format('delta',I2,F10.7,I2)
            122 format('delta',2I2,F10.7)
            130 format('object 2 class gridconnections counts',3I5)
            140 format('object 3 class array type double rank 0 items',I10,' data follows')
            150 format(ES19.10,' ') 
end subroutine gen_dx_file

subroutine gen_integer_dx_file(xm,ym,zm,h,gox,goy,goz,voldata,filename,fn,dataname)

    implicit none
    !passed variables
    integer xm,ym,zm,fn
    _REAL_  h, gox, goy, goz
    integer  voldata(1:xm,1:ym,1:zm)
    character*(*), intent(in) :: filename, dataname

    !internal variables
    integer cnt, i, j, k   
     
!               write(6,*) 'writing potential map in dx format'
                open(fn,file=filename)
                write(fn,200) xm,ym,zm
                write(fn,210) gox+h,goy+h,goz+h
                write(fn,220) sngl(h),0,0
                write(fn,221) 0,sngl(h),0
                write(fn,222) 0,0,sngl(h)
                write(fn,230) xm,ym,zm
                write(fn,240) xm*ym*zm
                cnt=0
                do i=1,xm; do j=1,ym; 
                    do k=1,zm
                        cnt = cnt+1 
                        write(fn,250,advance='no') voldata(i,j,k)
                        if (cnt >= 3) then
                            cnt = 0
                            write(fn,*)''
                        end if
                    end do
                end do; end do
                write(fn,*) 'attribute "dep" string "positions"'
                write(fn,*) 'object "',dataname,'" class field'
                write(fn,*) 'component "positions" value 1'
                write(fn,*) 'component "connections" value 2'
                write(fn,*) 'component "data" value 3'
                close(fn)
            200 format('object 1 class gridpositions counts',3I5)
            210 format('origin',3F9.4)
            220 format('delta',F10.7,2I2)
            221 format('delta',I2,F10.7,I2)
            222 format('delta',2I2,F10.7)
            230 format('object 2 class gridconnections counts',3I5)
            240 format('object 3 class array type double rank 0 items',I10,' data follows')
            250 format(I12,' ') 
end subroutine gen_integer_dx_file
subroutine printphidx(phigrid,xm,ym,zm,gox,goy,goz,h)
    use poisson_boltzmann, only : frcfac, eps0
    !Passed Variables
    integer xm,ym,zm
    _REAL_ phigrid(1:xm,1:ym,1:zm)
    _REAL_ gox,goy,goz,h
    !Internal Variables
    _REAL_ phitemp(xm,ym,zm)
    character*14 phi_filename
    character*23 phi_setname
    integer phi_filenum
    
    phi_filename = "pbsa_phi.dx"
    phi_setname = "Electrostatic Potential"
    phi_filenum = 8888

    phitemp=phigrid/eps0!*frcfac
    call gen_dx_file(xm,ym,zm,h,gox,goy,goz,phitemp,&
                    phi_filename,phi_filenum,phi_setname)

end subroutine printphidx
