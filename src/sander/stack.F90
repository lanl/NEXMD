! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

module stack

!--------------------------------------------------------------------------
!                      USAGE of stack
!------------------------------------------------------------------------- 
! 
!    Stack usage is through this module for both r_stack (real stack)
!                                           and  i_stack (integer stack)
!
!    The restrictions are that all calls to either stack must start with
!    a clean stack and finish with the stack freed, and all within the same
!    subroutine. The stack pointers can be passed to other routines but
!    the stack calls (get and free) must not originate from another routine
!    until the stack is cleared.
!
!    this restriction comes from the new freedom to increase the stack by
!    allocating a new stack and copying the contents of the old one in. If 
!    that occurs in a routine that is downstream from a routine that first 
!    called the get_stack or get_istack, the new pointer will be fine, but 
!    all the other space will that was passed in will be obliterated 
!    (deallocated) and lost....bombs away.
!
!    Each set of calls to get_stack and get_istack MUST have following it and
!    preceding the use of the stack arrays, a check that the stacks are OK.
!    Either "REQUIRE(rstack_ok)", "REQUIRE(istack_ok)" to cause a bomb if the
!    stack was too small, 
!    OR,
!      call get_stack(l_th3,num_ks*order,routine)
!      if(.not. rstack_ok)then
!         deallocate(r_stack)
!         allocate(r_stack(1:lastrst),stat=alloc_ier)
!         call reassign_rstack(routine)
!      endif
!      REQUIRE(rstack_ok)
!
!      call get_istack(imy_cg,num_ks,routine)
!      if(.not. istack_ok)then
!         deallocate(i_stack)
!         allocate(i_stack(1:lastist),stat=alloc_ier)
!         call reassign_istack(routine)
!      endif
!      REQUIRE(istack_ok)
!    to reallocate new stack and copy in anything relevant from the old stack
!     
!
!    To check the validity of the programming using stack, turn on 
!    DEBUG_STACK in AMBER_BUILDFLAGS. Verbose messages about allocation and
!    get_stack will be printed and calling routines will be checked.
!
!---------------------------------------------------------------------------

integer, parameter :: max_stack_ptrs=100
integer :: alloc_ier

!=========== RSTACK ===============================================

character (kind=1,len=80), private, save :: stk_call_routine, &
                                           ist_call_routine, &
                                           blank_routine=" "
integer,save :: len_stk_callrtn,len_ist_callrtn

integer,save :: bot_stk,top_stk,num_stkptrs, &
                last_stkptr,highest_stk,lastrst,top_stk_last
integer, save, allocatable, dimension(:) :: rstk_ptr
_REAL_, save, allocatable, target, dimension(:) :: r_stack_0, r_stack_1
_REAL_, save, allocatable, dimension(:) :: r_stack, r_stack_old
logical, save :: rstack_ok, rstack_resized

!=========== ISTACK ===============================================
  
integer,save :: ibot_stk,itop_stk,inum_stkptrs, &
                ilast_stkptr,ihighest_stk,lastist,itop_stk_last
integer, save, allocatable, dimension(:) :: istk_ptr

integer, save, allocatable, target, dimension(:) :: i_stack_0, i_stack_1
integer, save, allocatable, dimension(:) :: i_stack, i_stack_old
logical, save :: istack_ok, istack_resized

integer, parameter, private :: stk_verbose=1


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!------------------------------------------------------
!         GET_ISTACK
!------------------------------------------------------
subroutine get_istack(ipointer,isize,routine)
   implicit none
   integer,intent(out) :: ipointer
   integer,intent(in)  :: isize
   character(kind=1,len=*),intent(in) :: routine
#ifdef DEBUG_STACK
   integer :: len_routine

   if(stk_verbose > 1) &
        write(6,'("| get_istack called from ",a,1x,a,2i12)') &
        routine, ist_call_routine(1:len_ist_callrtn), isize,itop_stk
   ! --- Check for name of calling routine
   !     All calls must originate from same routine unless
   !      stack is clear

   len_routine = len_trim(routine)
   if(len_routine == 0)then
      write(6,*)"ERROR get_istack passed no routine name"
      call mexit(6,1)
   endif
   if(itop_stk == 0) then
      len_ist_callrtn = len_routine
      ist_call_routine(1:len_routine) = routine(1:len_routine)
   else
      if( (len_routine /= len_ist_callrtn) .or. &
         (index(ist_call_routine,routine) /= 1) ) then
         write(6,*)"ERROR get_istack called from several routines"
         write(6,*)"Legal routine       : ", &
               ist_call_routine(1:len_ist_callrtn),len_ist_callrtn
         write(6,*)"Illegal call routine: ",routine(1:len_routine)
         call mexit(6,1)
      endif
   endif
#endif

   ipointer=itop_stk+1

   if(itop_stk+isize > lastist )then
      istack_ok = .false.

#ifdef DEBUG_STACK
      write(6,'("| Exceeding lastist in get_istack called from ",a)')routine
      write(6,'("|  lastist = ",i12)')lastist
      write(6,'("|  itop_stk= ",i12)')itop_stk
      write(6,'("|  isize   = ",i12)')isize
      write(6,'("|  request = ",i12)')itop_stk+isize
#endif

      call resize_istack(itop_stk+isize,lastist)
      if(alloc_ier /= 0)then
         write(6,*) " REALLOCATION of int stack FAILED... "
         write(6,*) " Request that could not be allocated was size :",isize
         call mexit(6,1)
      endif
      istack_resized=.true.
!      write(6,'("|  lastist = ",i12)')lastist
   end if
   itop_stk = itop_stk+isize
   inum_stkptrs=inum_stkptrs+1
   ilast_stkptr=ipointer
   if(inum_stkptrs > max_stack_ptrs)then
      write(6,*) "Exceeding MAX_STACK_PTRS in get_istack() "
      write(6,*) "increase max_stack_ptrs in stack.f"
      call mexit(6,1)
   end if
   istk_ptr(inum_stkptrs)=ipointer
   ihighest_stk=max(ihighest_stk,itop_stk)
!  write(6,'(a,2i10)') '  get_istack: ', itop_stk, isize
   return
 end subroutine get_istack

!------------------------------------------------------
!         GET_STACK
!------------------------------------------------------
subroutine get_stack(ipointer,isize,routine)
   implicit none
   integer, intent(out) :: ipointer
   integer, intent(in)  :: isize
   character(kind=1,len=*),intent(in) :: routine
#ifdef DEBUG_STACK
   integer :: len_routine
   
   ! --- Check for name of calling routine
   !     All calls must originate from same routine unless
   !      stack is clear

   !write(6,*)"Get_stack called from",routine,isize
   len_routine = len_trim(routine)
   if(len_routine == 0)then
      write(6,*)"ERROR get_stack passed no routine name"
      call mexit(6,1)
   endif
   if(top_stk == 0) then
      len_stk_callrtn = len_routine
      stk_call_routine(1:len_routine) = routine(1:len_routine)
   else
      if( (len_routine /= len_stk_callrtn) .or. &
         (index(stk_call_routine,routine) /= 1) ) then
         write(6,*)"ERROR get_stack called from several routines"
         write(6,*)"Legal routine  : ",stk_call_routine(1:len_stk_callrtn), &
               &   "Illegal routine: :",routine(1:len_routine)
         call mexit(6,1)
      endif
   endif
#endif
         
            
   ipointer=top_stk+1
   if(top_stk+isize > lastrst )then
      rstack_ok=.false.
#ifdef DEBUG_STACK
      write(6,'("| Exceeding lastrst in get_stack from",a)')routine
      write(6,'("|  lastrst = ",i12)')lastrst
      write(6,'("|  top_stk= ",i12)')top_stk
      write(6,'("|  isize   = ",i12)')isize
      write(6,'("|  request = ",i12)')top_stk+isize
#endif
      call resize_stack(top_stk+isize,lastrst)
      if(alloc_ier /= 0)then
         write(6,*) " REALLOCATION of real stack FAILED... "
         write(6,*) " Request that could not be allocated was size :",isize
         call mexit(6,1)
      endif
      rstack_resized=.true.
   end if
   top_stk = top_stk+isize
   num_stkptrs=num_stkptrs+1
   last_stkptr=ipointer
   if(num_stkptrs > max_stack_ptrs)then
      write(6,*) "Exceeding MAX_STACK_PTRS in get_stack() "
      write(6,*) "increase MAX_STACK_PTRS in stack.f"
      call mexit(6,1)
   end if
   rstk_ptr(num_stkptrs)=ipointer
   highest_stk=max(highest_stk,top_stk)
   return
end subroutine get_stack 
 
!------------------------------------------------------
!         FREE_ISTACK
!------------------------------------------------------
subroutine free_istack(ipointer,routine)
   implicit none
   integer ipointer
   character(kind=1,len=*),intent(in) :: routine
#ifdef DEBUG_STACK
   integer :: len_routine

   if(stk_verbose > 1) &
   write(6,'("| free_istack called from ",2(a,1x),i12)')routine, &
      ist_call_routine(1:len_ist_callrtn),ipointer 
   len_routine = len_trim(routine)
   if(len_routine == 0)then
      write(6,*)"ERROR free_istack passed no routine name"
      call mexit(6,1)
   endif
   if(itop_stk == 0) then
      ist_call_routine = blank_routine
      len_ist_callrtn = 0
   else
      if( (len_routine /= len_ist_callrtn) .or. &
         (index(ist_call_routine,routine) /= 1) ) then
         write(6,*)"ERROR free_istack called from several routines"
         write(6,*)"Legal routine  : ",ist_call_routine(1:len_ist_callrtn), &
               &   "Illegal routine: :",routine(1:len_routine)
         call mexit(6,1)
      endif
   endif
#endif
   
   if(ipointer /= ilast_stkptr)then
      write(6,*)"Trying to free int stack from the middle"
      write(6,*)"highest pointer is ",ilast_stkptr
      write(6,*)"free req ptr is    ",ipointer
      write(6,*)"   Must free from top (last one made freed first"
      call mexit(6,1)
   end if
   inum_stkptrs=inum_stkptrs-1
   if(inum_stkptrs < 0)then
      write(6,*)"Problem in int stack, freed too much."
      write(6,*)" There are zero pointers left on the stack"
      write(6,*)" cannot free anymore"
      call mexit(6,1)
   end if
   if (inum_stkptrs > 0)ilast_stkptr=istk_ptr(inum_stkptrs)
   itop_stk = ipointer-1
   if(itop_stk < ibot_stk)then
      write(6,*)"Problem in integer stack, freed too much."
      write(6,*)"top of stack is below the top of heap"
      call mexit(6,1)
   end if
#ifdef DEBUG_STACK
   if(stk_verbose > 1) then
      if(itop_stk == 0)write(6,'(3(a,1x))') &
            '| free_istack clear ', routine,ist_call_routine(1:len_ist_callrtn)
   endif
#endif
   return
end subroutine free_istack 
!------------------------------------------------------
!         FREE_STACK
!------------------------------------------------------
subroutine free_stack(ipointer,routine)
   implicit none
   integer, intent(in) :: ipointer
   character(kind=1,len=*),intent(in) :: routine
#ifdef DEBUG_STACK
   integer :: len_routine
   
   len_routine = len_trim(routine)
   if(len_routine == 0)then
      write(6,*)"ERROR free_stack passed no routine name"
      call mexit(6,1)
   endif
   if(top_stk == 0) then
      stk_call_routine = blank_routine
      len_stk_callrtn = 0
   else
      if( (len_routine /= len_stk_callrtn) .or. &
         (index(stk_call_routine,routine) /= 1) ) then
         write(6,*)"ERROR free_stack called from several routines"
         write(6,*)"Legal routine  : ",stk_call_routine(1:len_stk_callrtn), &
               &   "Illegal routine: :",routine(1:len_routine)
         call mexit(6,1)
      endif
   endif
#endif
   
   if(ipointer /= last_stkptr)then
      write(6,*) "Trying to free stack from the middle"
      write(6,*) "highest pointer is ",last_stkptr
      write(6,*) "free req ptr is    ",ipointer
      write(6,*) "   Must free from top (last one freed first)"
      call mexit(6,1)
   end if
   num_stkptrs=num_stkptrs-1
   if(num_stkptrs < 0)then
      write(6,*) "Problem in stack, freed too much."
      write(6,*) " There are zero pointers left on the stack"
      write(6,*) " cannot free anymore"
      call mexit(6,1)
   end if
   if(num_stkptrs > 0)last_stkptr=rstk_ptr(num_stkptrs)
   top_stk = ipointer-1
   if(top_stk < bot_stk)then
      write(6,*) "Problem in stack, freed too much, corrupt pointer"
      write(6,*) "top of stack is below the top of heap"
      call mexit(6,1)
   end if
   return
end subroutine free_stack 

!------------------------------------------------------
!         STACK_SETUP
!------------------------------------------------------
subroutine stack_setup()
   implicit none
   integer ier
   itop_stk=0
   ibot_stk=0
   ihighest_stk=0
   inum_stkptrs=0
   allocate(istk_ptr(max_stack_ptrs), stat = ier )
   REQUIRE( ier == 0 )
   allocate( i_stack(1:lastist), stat = ier )
   REQUIRE( ier == 0 )
   istack_ok = .true.
   istack_resized = .false.
   itop_stk_last=0

   stk_call_routine="get_stack"
   len_stk_callrtn = len_trim(stk_call_routine)
   ist_call_routine="get_istack"
   len_ist_callrtn = len_trim(ist_call_routine)

   top_stk=0
   bot_stk=0
   highest_stk=0
   num_stkptrs=0
   allocate(rstk_ptr(max_stack_ptrs), stat = ier )
   REQUIRE( ier == 0 )
   allocate( r_stack(1:lastrst), stat = ier )
   REQUIRE( ier == 0 )
   rstack_ok = .true.
   rstack_resized = .false.
   top_stk_last=0
   return
end subroutine stack_setup

 
!------------------------------------------------------
!         DEALLOCATE_STACKS
!------------------------------------------------------
subroutine deallocate_stacks()
   implicit none
   integer ier

   if(top_stk >0 )write(6,'("| ERROR top_stk not zero!!!",i12)')top_stk
   if(itop_stk >0 )write(6,'("| ERROR itop_stk not zero!!!",i12)')itop_stk
   deallocate(istk_ptr, stat = ier )
   REQUIRE( ier == 0 )
   deallocate(rstk_ptr, stat = ier )
   REQUIRE( ier == 0 )
   
   if(allocated(i_stack))then
      deallocate( i_stack, stat = ier )
      REQUIRE( ier == 0 )
   endif
   if(allocated(r_stack))then
      deallocate( r_stack, stat = ier )
      REQUIRE( ier == 0 )
   endif
   return
end subroutine deallocate_stacks

!------------------------------------------------------
!         RESIZE_ISTACK
!------------------------------------------------------
subroutine resize_istack(isize,lastist)
  integer, intent(in) :: isize
  integer, intent(inout) :: lastist

!  write(6,'("| Reallocating istack from ", i12, " to ",i12)')lastist,isize
 
   alloc_ier=0
   if(.not. istack_resized) then
       allocate(i_stack_old(1:itop_stk),stat=alloc_ier)
       if(itop_stk > 0 ) then
          i_stack_old(1:itop_stk) = i_stack(1:itop_stk)
       endif
       itop_stk_last=itop_stk
    endif
   lastist=isize
   return
end subroutine resize_istack


!------------------------------------------------------
!         REASSIGN_ISTACK
!------------------------------------------------------
subroutine reassign_istack(routine)
  implicit none
  character(kind=1,len=*), intent(in) :: routine

  if(alloc_ier /= 0 ) then
     write(6,*)"Unable to reallocate more stack"
     write(6,*)"Routine ",routine,"  size needed: ",lastist
     call mexit(6,1)
  endif
  if(itop_stk_last > 0 ) then
     i_stack(1:itop_stk_last) = i_stack_old(1:itop_stk_last)
  endif
  deallocate(i_stack_old)
  istack_resized=.false.
  istack_ok = .true.
  
  return
end subroutine reassign_istack

!------------------------------------------------------
!         RESIZE_STACK
!------------------------------------------------------
subroutine resize_stack(isize,lastrst)
  integer, intent(in) :: isize
  integer, intent(out) :: lastrst
  
  
!  write(6,'("| Reallocating r_stack from ", i12, " to ",i12)')lastrst,isize
 
  alloc_ier=0
  if(.not. rstack_resized) then
     top_stk_last = top_stk
     if(top_stk > 0 ) then
        allocate(r_stack_old(1:top_stk),stat=alloc_ier)
        if(alloc_ier /= 0 ) then
           write(6,*)"Unable to reallocate more stack"
           return
        endif
        r_stack_old(1:top_stk) = r_stack(1:top_stk)
     endif
     rstack_resized=.true.
  endif
  lastrst=isize
  return
end subroutine resize_stack

!------------------------------------------------------
!         REASSIGN_RSTACK
!------------------------------------------------------
subroutine reassign_rstack(routine)
  implicit none
  character(kind=1,len=*), intent(in) :: routine

  if(alloc_ier /= 0 ) then
     write(6,*)"Unable to reallocate more stack"
     write(6,*)"Routine ",routine,"  size needed: ",lastrst
     call mexit(6,1)
  endif
  if(top_stk_last > 0 ) then
     r_stack(1:top_stk_last) = r_stack_old(1:top_stk_last)
     deallocate(r_stack_old) 
  endif
  rstack_resized=.false.
  rstack_ok = .true.
  
  return
end subroutine reassign_rstack

end module stack
