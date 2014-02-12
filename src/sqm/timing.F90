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
   implicit none

   character*30 datetime

   character*20 str
   integer system,res

   str='date > syst'
   res=system(str)

   if(res==0) then ! everything is ok
      open(12,file='syst',status='old')
      read(12,'(A30)') datetime
      close(12,status='delete')
   else
      datetime='system time problem'
   end if

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
