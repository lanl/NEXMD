	function get_time()
	integer tms(4), get_time
	real*4 etime,etime1,tarry(2)
	
!        call times(tms)           ! Linux setting
!        get_time=tms(1)+tms(2)    ! Linux setting
!        get_time= mclock()	   ! SGI IRIX setting
	 get_time=100*etime(tarry) ! DEC Alpha setting
        return
	end    

	subroutine get_date(datetime)
	
	character*20 datetime
	character*10 chdate	
	character*8 chtime

      character*8 date
      character*10 time
      character*5 zone
      integer*4 values(8)

      datetime='                   '
	chtime='        '
	chdate='          '	
	call date(chdate)
	call time(chtime)
!      call date_and_time(date,time,zone,values)

!      datetime(1:2)=date(5:6)
!      datetime(3:3)='/'
!      datetime(4:5)=date(7:8)
!      datetime(6:6)='/' 
!      datetime(7:10)=date(1:4)
!      datetime(11:11)=' '
!      datetime(12:13)=time(1:2)
!      datetime(14:14)= 'h'
!      datetime(15:16)=time(3:4)
!      datetime(17:17)= 'm'
!      datetime(18:19)=time(5:6)     
!      datetime(20:20)= 's'
				
!	print *,datetime 			
	datetime(1:10)=chdate
	datetime(11:11)=' '
	datetime(12:19)=chtime

        return
	end


	subroutine get_machine(machname)
	
	character*36 machname, str*20
	integer system, i
	
	str='uname -nms > syst'
      i=system(str)
	
	if (i.eq.0) then
	  open (2,file='syst')
	  read (2,10) machname 
	  close (2,status='delete')
	endif 	
 10   format (A36)

        return
	end
