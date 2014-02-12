#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Instrumentation facility to log diagnostic information.

module trace
   
   ! Description:
   ! Module for tracing execution flow.
   ! Activation of tracing is currently controlled by the
   ! preprocessor variable TRACE.
   ! This module is accessed via a use statement.
   ! The use statement must occur before any implicit
   ! statements or variable declarations.
   ! parallel.f contains examples of the usage of this module.
   ! This is a traditional instrumentation facility.
   ! Presently, emitted information is sent to the standard output
   ! and is prepended with the token "Trace:".
   ! Following this token is a possibly null string of blanks that
   ! indicates the nesting depth of subprogram calls.
   
   ! The current set of public tracing routines includes:
   !   Subroutine Trace_debug( filename, line, message )
   !          log a message, prepending the filename and line number
   !   Subroutine Trace_enter( subprogram_name )
   !          log entry into a function or subroutine
   !   Subroutine Trace_exit( subprogram_name )
   !          log exit from a function or subroutine
   !   Subroutine Trace_integer( label, value )
   !          log an integer
   !   Subroutine Trace_logical( label, value )
   !          log an logical
   !   Subroutine Trace_mpi( name, size, type, rank )
   !          log an MPI call
   !   Subroutine Trace_note( message )
   !          log a message
   !   Subroutine Trace_output_mpi_tally( )
   !          log the mpi tallies and restart the per step counters
   !   Subroutine Trace_real( label, value )
   !          log a real
   !   Subroutine Trace_set_ouput_focus( mpi_rank )
   !          assign the output focus to the processor of rank mpi_rank
   !           default rank is 0.
   !   Subroutine Trace_table( i1, r1, i2, r2, ... )
   !          log a table row consisting of column 1 (i1 xor r1),
   !           column 2 (i2 xor r2), etc.
   !           all arguments are optional.
   !           Caution g77 emits syntax errors for optional arguments;
   !           calls to Trace_table must be commented out.
   !           Example usage:
   !             call Trace_note( "     i     j     k  data" )
   !             do i = 1, IMAX ; do j = 1, JMAX ; do k = 1 to KMAX
   !               call Trace_table( i1=i, i2=j, i3=k, r4=data(i,j,k) )
   !             enddo; enddo; enddo
   
   ! The current set of public tracing parameters includes:
   !     TRACE_DPREC
   !     TRACE_INT
   !     TRACE_LOG
   ! which are abbreviations for the corresponding MPI data types as
   ! strings.
   
   ! History:
   ! $Id: trace.f,v 10.1 2009/08/22 02:18:16 case Exp $
   
   ! Code Description:
   !   Language:           Fortran 90.
   !   Software Standards: "European Standards for Writing and
   !     Documenting Exchangeable Fortran 90 Code":
   !     http://www.physics.reading.ac.uk/CompPhys/fortran90/F90Style.htm
   

   ! **********************************************************************
   !  Copyright 2003                                                      *
   !                                                                      *
   !   Modified BSD license                                               *
   !                                                                      *
   !   Redistribution and use in source and binary forms, with or without *
   !   modification, are permitted provided that the following conditions *
   !   are met:                                                           *
   !                                                                      *
   !    1.Redistributions of source code must retain the above copyright  *
   !      notice, this list of conditions and the following disclaimer.   *
   !    2.Redistributions in binary form must reproduce the above         *
   !      copyright notice, this list of conditions and the following     *
   !      disclaimer in the documentation and/or other materials provided *
   !      with the distribution.                                          *
   !    3.The name of the author may not be used to endorse or promote    *
   !      products derived from this software without specific prior      *
   !      written permission.                                             *
   !                                                                      *
   !   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS''                  *
   !   AND ANY EXPRESS OR IMPLIED WARRANTIES,                             *
   !   INCLUDING, BUT NOT LIMITED TO, THE IMPLIED                         *
   !   WARRANTIES OF MERCHANTABILITY AND FITNESS FOR                      *
   !   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT                   *
   !   SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,                         *
   !   INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR                       *
   !   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT                          *
   !   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR                     *
   !   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR                        *
   !   BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON                       *
   !   ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,                      *
   !   STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE                    *
   !   OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE                    *
   !   OF THIS SOFTWARE, EVEN IF ADVISED OF THE                           *
   !   POSSIBILITY OF SUCH DAMAGE.                                        *
   !                                                                      *
   !  To report bugs, suggest enhancements, etc., contact:                *
   !    Scott Brozell                                                     *
   !                  send email to sbrozell@chemistry.ohio-state.edu     *
   !                                                                      *
   ! **********************************************************************

   implicit none

   private
   public :: trace_dprec
   public :: trace_int
   public :: trace_log
   public :: trace_debug
   public :: trace_enter
   public :: trace_exit
   public :: trace_integer
   public :: trace_logical
   public :: trace_mpi
   public :: trace_note
   public :: trace_output_mpi_tally
   public :: trace_real
   public :: trace_set_ouput_focus
   public :: trace_table

   character(*), parameter :: trace_dprec = 'MPI_DOUBLE_PRECISION'
   character(*), parameter :: trace_int   = 'MPI_INTEGER'
   character(*), parameter :: trace_log   = 'MPI_LOGICAL'


   type tally
      integer             :: per_step
      integer             :: total
   end type tally

   integer                 :: indentation_level = 0
   character(*), parameter :: indentation_spacing = " "
   type( tally )           :: mpi_calls = tally( 0, 0 )
   type( tally )           :: mpi_sizes = tally( 0, 0 )
   character(*), parameter :: output_token = "Trace:"
   integer                 :: processor_with_output_focus = 0
#ifdef TRACE
   logical                 :: is_tracing_on = .true.
#else
   logical                 :: is_tracing_on = .false.
#endif


   contains
   ! public routines in alphabetical order

   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_debug( filename, line, message )
      !
      ! log a message, prepending the filename and line number
      !
      !------------------------------------------------------------------
      implicit none
      character(*), intent(in) :: filename
      integer,      intent(in) :: line
      character(*), intent(in) :: message
      if ( is_trace_active() .and. may_node_emit() ) then
         call indent()
         write ( *, '(2A,I6,2A)' ) filename, ": ", line, ": ", message
      end if
   end subroutine Trace_debug


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_enter( subprogram_name )
      !
      ! log entry into a function or subroutine
      !
      !------------------------------------------------------------------
      implicit none
      character(*), intent(in) :: subprogram_name
      _REAL_                   :: time
      if ( is_trace_active() .and. may_node_emit() ) then
         call indent()
         call wallclock( time )
         write ( *, '(3A,G19.14,A)' ) "Enter " &
               , subprogram_name, "   at ", time, " seconds"
         call indentation_increment()
      end if
   end subroutine Trace_enter


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_exit( subprogram_name )
      !
      ! log exit from a function or subroutine
      !
      !------------------------------------------------------------------
      implicit none
      character(*), intent(in) :: subprogram_name
      _REAL_                   :: time
      if ( is_trace_active() .and. may_node_emit() ) then
         call indentation_decrement()
         call indent()
         call wallclock( time )
         write ( *, '(3A,G19.14,A)' ) "Exit  " &
               , subprogram_name, "   at ", time, " seconds"
      end if
   end subroutine Trace_exit


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_integer( label, value )
      !
      ! log an integer
      !
      !------------------------------------------------------------------
      implicit none
      character(*), intent(in) :: label
      integer,      intent(in) :: value
      if ( is_trace_active() .and. may_node_emit() ) then
         call indent()
         write ( *, '(A,I11)' ) label, value
      end if
   end subroutine Trace_integer


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_logical( label, value )
      !
      ! log a logical
      !
      !------------------------------------------------------------------
      implicit none
      character(*), intent(in) :: label
      logical,      intent(in) :: value
      if ( is_trace_active() .and. may_node_emit() ) then
         call indent()
         write ( *, '(A,L11)' ) label, value
      end if
   end subroutine Trace_logical


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_mpi( name, size, type, rank )
      !
      ! log an MPI call
      !
      !------------------------------------------------------------------
      implicit none
      character(*), intent(in) :: name
      integer,      intent(in) :: size
      character(*), intent(in) :: type
      integer,      intent(in) :: rank
      if ( is_trace_active() .and. may_node_emit() ) then
         call indent()
         write ( *, '(2A,I7,3A,I3)' ) name, ": size = ", size &
               , ": type = ", type, ": rank = ", rank
         call increase_mpi_tally( 1, size )
      end if
   end subroutine Trace_mpi


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_note( message )
      !
      ! log a message
      !
      !------------------------------------------------------------------
      implicit none
      character(*), intent(in) :: message
      if ( is_trace_active() .and. may_node_emit() ) then
         call indent()
         write ( *, '(A)' ) message
      end if
   end subroutine Trace_note


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_output_mpi_tally( )
      !
      ! log the mpi tallies and restart the per step counters
      !
      !------------------------------------------------------------------
      implicit none
      if ( is_trace_active() .and. may_node_emit() ) then
         call indent()
         write ( *, '(A,I11)' ) "MPI calls this step " &
               , mpi_calls%per_step
         call indent()
         write ( *, '(A,I11)' ) "MPI sizes this step " &
               , mpi_sizes%per_step
         call indent()
         write ( *, '(A,I11)' ) "MPI total calls " &
               , mpi_calls%total
         call indent()
         write ( *, '(A,I11)' ) "MPI total sizes " &
               , mpi_sizes%total
         call restart_mpi_per_step_tally()
      end if
   end subroutine Trace_output_mpi_tally


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_real( label, value )
      !
      ! log a real
      !
      !------------------------------------------------------------------
      implicit none
      character(*), intent(in) :: label
      _REAL_ ,      intent(in) :: value
      if ( is_trace_active() .and. may_node_emit() ) then
         call indent()
         write ( *, '(A,G12.5)' ) label, value
      end if
   end subroutine Trace_real


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_set_ouput_focus( mpi_rank )
      !
      ! assign the output focus to the processor of rank mpi_rank
      !
      !------------------------------------------------------------------
      implicit none
      integer,      intent(in) :: mpi_rank
      if ( is_trace_active() ) then
         processor_with_output_focus = mpi_rank
      end if
   end subroutine Trace_set_ouput_focus


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine Trace_table( i1,r1, i2,r2, i3,r3, i4,r4, i5,r5, &
         i6,r6, i7,r7 )
      !
      ! log a table row consisting of column 1 (i1 xor r1),
      ! column 2 (i2 xor r2), etc.
      ! all arguments are optional.
      ! Example invocation:
      !   call Trace_table( i1=i, i2=j, i3=k, r4=data(i,j,k) )
      !
      !------------------------------------------------------------------
      implicit none
      integer,      optional, intent(in) :: i1
      _REAL_ ,      optional, intent(in) :: r1
      integer,      optional, intent(in) :: i2
      _REAL_ ,      optional, intent(in) :: r2
      integer,      optional, intent(in) :: i3
      _REAL_ ,      optional, intent(in) :: r3
      integer,      optional, intent(in) :: i4
      _REAL_ ,      optional, intent(in) :: r4
      integer,      optional, intent(in) :: i5
      _REAL_ ,      optional, intent(in) :: r5
      integer,      optional, intent(in) :: i6
      _REAL_ ,      optional, intent(in) :: r6
      integer,      optional, intent(in) :: i7
      _REAL_ ,      optional, intent(in) :: r7
      if ( is_trace_active() .and. may_node_emit() ) then
         call indent()
         if ( present( i1 ) ) then
            write ( *, '(I9)', advance = "no" ) i1
         end if
         if ( present( r1 ) ) then
            write ( *, '(G12.5)', advance = "no" ) r1
         end if
         if ( present( i2 ) ) then
            write ( *, '(I9)', advance = "no" ) i2
         end if
         if ( present( r2 ) ) then
            write ( *, '(G12.5)', advance = "no" ) r2
         end if
         if ( present( i3 ) ) then
            write ( *, '(I9)', advance = "no" ) i3
         end if
         if ( present( r3 ) ) then
            write ( *, '(G12.5)', advance = "no" ) r3
         end if
         if ( present( i4 ) ) then
            write ( *, '(I9)', advance = "no" ) i4
         end if
         if ( present( r4 ) ) then
            write ( *, '(G12.5)', advance = "no" ) r4
         end if
         if ( present( i5 ) ) then
            write ( *, '(I9)', advance = "no" ) i5
         end if
         if ( present( r5 ) ) then
            write ( *, '(G12.5)', advance = "no" ) r5
         end if
         if ( present( i6 ) ) then
            write ( *, '(I9)', advance = "no" ) i6
         end if
         if ( present( r6 ) ) then
            write ( *, '(G12.5)', advance = "no" ) r6
         end if
         if ( present( i7 ) ) then
            write ( *, '(I9)', advance = "no" ) i7
         end if
         if ( present( r7 ) ) then
            write ( *, '(G12.5)', advance = "no" ) r7
         end if
         write ( *, '(A)' ) " "
      end if  ! ( is_trace_active() .and. may_node_emit() )
   end subroutine Trace_table


   ! private routines in alphabetical order


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine indent()
      !
      ! prepends nesting level indentation to the log output
      !
      !------------------------------------------------------------------
      implicit none

      integer i
      write( *, '(A,I3,A)', advance = "no" ) output_token, &
            processor_with_output_focus, ':'
      do i = 1, indentation_level
         write( *, '(A)', advance = "no" ) indentation_spacing
      end do
   end subroutine indent


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine indentation_decrement()
      !
      ! decreases the indentation level for the log output
      !
      !------------------------------------------------------------------
      implicit none

      indentation_level = indentation_level - 1
      ASSERT( indentation_level >= 0 )
   end subroutine indentation_decrement


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine indentation_increment()
      !
      ! increases the indentation level for the log output
      !
      !------------------------------------------------------------------
      implicit none

      indentation_level = indentation_level + 1
      ASSERT( indentation_level > 0 )
   end subroutine indentation_increment


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine increase_mpi_tally( calls, size )
      !
      ! increases the various mpi tallies
      !
      !------------------------------------------------------------------
      implicit none
      integer,      intent(in) :: calls
      integer,      intent(in) :: size

      ASSERT( calls >= 0 )
      ASSERT( size >= 0 )
      mpi_calls%per_step = mpi_calls%per_step + calls
      mpi_sizes%per_step = mpi_sizes%per_step + size
      mpi_calls%total    = mpi_calls%total    + calls
      mpi_sizes%total    = mpi_sizes%total    + size
   end subroutine increase_mpi_tally


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   function may_node_emit()
      !
      ! returns true if the node is allowed to produce output
      !
      !------------------------------------------------------------------
      implicit none
#include "parallel.h"
      logical :: may_node_emit

      may_node_emit = mytaskid == processor_with_output_focus
   end function may_node_emit


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   subroutine restart_mpi_per_step_tally()
      !
      ! zeros the various per step mpi tallies
      !
      !------------------------------------------------------------------
      implicit none

      mpi_calls%per_step = 0
      mpi_sizes%per_step = 0
   end subroutine restart_mpi_per_step_tally


   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   ! This subprogram will likely become public when runtime control of
   ! tracing is implemented.

   function is_trace_active()
      !
      ! returns true if tracing is active, false if tracing is unactive
      !
      !------------------------------------------------------------------
      implicit none
      logical :: is_trace_active

      is_trace_active = is_tracing_on
   end function is_trace_active


end module trace

