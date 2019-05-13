! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"

!-------------------------------------------------------------------
!     --- ASSIGN_IND ---
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine assign_ind here]
subroutine assign_ind(ind0,ind)
   integer ind0,ind
   ind0 = ind
   ind = ind+1
   return
end subroutine assign_ind 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute and emit net charge, neutralizing that from roundoff error.
!-----------------------------------------------------------------------
!     --- CHECK_NEUTRAL ---
!-----------------------------------------------------------------------
!     In general, the Ewald method is only truly applicable
!     under conditions of charge neutrality.  When the system is not
!     net neutral, the direct sum and reciprocal sums are not "beta"
!     independent.  Regardless, the Ewald method can be applied
!     with the fictitious assumption that there is a
!     "uniform net neutralizing plasma".
!
!     This routine will remove any net charge resulting from
!     conversion of the low precision parm topology charges.

subroutine check_neutral(charge,natom)

   use constants, only : INV_AMBER_ELECTROSTATIC, TEN_TO_MINUS2, zero
   implicit none
   _REAL_  charge(*), sum
   integer natom, i

#  include "box.h"
#  include "ew_cntrl.h"
#  include "extra.h"

   sum = ZERO
   do i = 1,natom
      sum = sum + charge(i)
   end do

   if (master) then
      write(6, '(/,5x,a,f12.8)') 'Sum of charges from parm topology file = ', &
         sum * INV_AMBER_ELECTROSTATIC
   end if
   
   ! Amber6 employed ischrgd to control charge neutralization.
   ! But this seems too dangerous, and users could easily forget
   ! to set ischgrd=1.  Moreover, the only time it makes sense to
   ! force neutrality is when it is just roundoff error.  So now
   ! we ignore ischrgd, and do the right thing:
   
   if ( abs( sum*INV_AMBER_ELECTROSTATIC ) > TEN_TO_MINUS2 ) then
      ! Significant nonzero net charge
      if ( ntb /= 0 ) then 
         ! Periodic simulation
         if ( use_pme /= 0 ) then
            ! Both direct and reciprocal Ewald
            if (master) write(6, '(5x,a)') &
               'Assuming uniform neutralizing plasma'
         else
            ! No reciprocal Ewald
            if (master) write(6, '(5x,a,a)') 'Unusual simulation - ', &
               'Periodic and No Reciprocal Ewald and Nonzero Net Charge'
         end if
      end if
   else
      ! Insignificant nonzero net charge, assuming due to roundoff error
      if (master) write(6, '(5x,a)') 'Forcing neutrality...'
      sum = sum/natom
      do i = 1,natom
         charge(i) = charge(i) - sum
      end do
   end if

   return
end subroutine check_neutral 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize direct Ewald Coulomb switch function.
!-----------------------------------------------------------------------
!     --- INIT_COULOMB_SWITCH ---
!-----------------------------------------------------------------------
subroutine init_coulomb_switch(cutoffnb,dxdr,eed_cub,eed_lin, &
      eedtbdns,eedmeth,ee_type)

   use constants, only : zero, one, half, third, tenth, TEN_TO_MINUS4
   implicit none
   _REAL_  eed_lin(2,*),eed_cub(4,*), cutoffnb,dxdr,eedtbdns
   integer eedmeth,ee_type

   _REAL_ x,xx,dx,switch,d_switch_dx
   _REAL_ err,derr,maxerr,maxderr,del,small
   _REAL_ aswitch,ad_switch_dx,xerr,dxerr
   integer ind
   integer i,j,max,min

#  include "extra.h"

   if ( master )then
      write(6,*) &
            '---------------------------------------------------'
   end if
   if ( eedmeth == 3 )then
      if ( master )then
         write(6, '(2(/,5x,a))' ) &
               'Using subroutine call for high precision switch', &
               'FOR FASTER DIRECT SUM try eedmeth =1 (cubic) or 2 (linear)'
         write(6,*) &
               '---------------------------------------------------'
      end if
      return
   end if
   if ( eedmeth == 4 )then
      if ( master )then
         write(6, '(2(/,5x,a))' ) &
               'eedmeth=4: Setting switch to one everywhere'
         write(6,*) &
               '---------------------------------------------------'
      end if
      return
   end if
   if ( eedmeth == 5 )then
      if ( master )then
         write(6, '(2(/,5x,a))' ) &
               'eedmeth=5: Using 1/r dielectric'
         write(6,*) &
               '---------------------------------------------------'
      end if
      return
   end if
   if ( eedmeth == 6 )then
      if ( master )then
         write(6, '(2(/,5x,a))' ) &
               'eedmeth=6: Using IPS method for electrostatic energy'
         write(6,*) &
               '---------------------------------------------------'
      end if
      return
   end if
   small = TEN_TO_MINUS4
   max = dxdr*cutoffnb*eedtbdns
   min = dxdr*eedtbdns
   min = 1
   maxerr = ZERO
   maxderr = ZERO
   xerr = -ONE
   dxerr = -ONE
   del = ONE / eedtbdns

   if ( master )then
      if ( eedmeth == 1 )then
         write(6,1000)
      else if ( eedmeth == 2 ) then 
         write(6,1100)
      else
         ASSERT( .false. )  ! eedmeth input validation occurs in mdread.f
      end if
      write(6,1001) eedtbdns
      write(6,1002)
   end if
   1000 format(1x,'APPROXIMATING switch and d/dx switch', &
         ' using CUBIC SPLINE INTERPOLATION')
   1100 format(1x,'APPROXIMATING switch and d/dx switch', &
         ' using LINEAR INTERPOLATION')
   1001 format(1x,'using ',f8.1,' points per unit in tabled values')
   1002 format(1x,'TESTING RELATIVE ERROR over r ranging from ', &
         '0.0 to cutoff')

   do i = min,max
      do j = 1,10
         x = del*(i-1) + j*TENTH*del
         call get_ee_func(x,switch,d_switch_dx,ee_type)
         if ( eedmeth == 1 )then
            
            !           -- cubic spline on switch:
            
            ind = eedtbdns*x + 1
            dx = x - (ind-ONE)*del
            aswitch = eed_cub(1,ind)+dx*(eed_cub(2,ind)+ &
                  dx*(eed_cub(3,ind)+dx*eed_cub(4,ind)*THIRD)*HALF)
            ad_switch_dx = eed_cub(2,ind)+dx*(eed_cub(3,ind)+ &
                  dx*eed_cub(4,ind)*HALF)
            
         else if ( eedmeth == 2 )then
            
            !           -- linear lookup on switch:
            
            xx = eedtbdns*x + 1
            ind = xx
            dx = xx - ind
            aswitch = (ONE - dx)*eed_lin(1,ind) + &
                  dx*eed_lin(1,ind+1)
            ad_switch_dx = (ONE - dx)*eed_lin(2,ind) + &
                  dx*eed_lin(2,ind+1)
         end if
         err = abs(aswitch - switch) / (abs(switch)+small)
         derr = abs(d_switch_dx - ad_switch_dx) / &
               (abs(d_switch_dx)+small)
         
         
         if ( err > maxerr )then
            maxerr = err
            xerr = x
         end if
         if ( derr > maxderr )then
            maxderr = derr
            dxerr = x
         end if
      end do  !  j = 1,10
   end do  !  i = min,max

   if ( master )then
      write(6,1003)maxerr,xerr
      write(6,1004)maxderr,dxerr
      write(6,*) &
            '---------------------------------------------------'
   end if
   1003 format('| CHECK switch(x): max rel err = ',e12.4, &
         '   at ',f10.6)
   1004 format('| CHECK d/dx switch(x): max rel err = ',e12.4, &
         '   at ',f10.6)
   return
end subroutine init_coulomb_switch 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine test_prime_factors here]
!-----------------------------------------------------------------------
!     --- PRIME FACTORS ---
!-----------------------------------------------------------------------

subroutine test_prime_factors(string,n)
   implicit none
   character(len=*) string
   integer n,result
   call check_prime_factors(n,result)
   if ( result == 0 )then
      write(6,*)'Checking for valid grid dimensions'
      write(6,60)string,n
      60 format(1x,a,'= ',i5,' not a product of 2s, 3s and 5s')
      call mexit(6,1)
   end if
   return
end subroutine test_prime_factors 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ result is 1 if positive integer n is a product of powers of 2,3,5.
subroutine check_prime_factors(n,result)
   implicit none
   integer, intent(in)   :: n
   ! result is 0 if n is not a product of powers of 2,3,5
   integer, intent(out)  :: result
   integer nq,nl

   ASSERT ( n > 0 )

   ! don't hurt n
   nl = n
   result = 1
   !     first divide down by 2
   10 continue
   if ( nl == 1 )return
   nq = nl / 2
   if ( nl /= nq * 2 )goto 20
   nl = nq
   goto 10
   20 continue
   !     now try 3
   30 continue
   if ( nl == 1 )return
   nq = nl / 3
   if ( nl /= nq * 3 )goto 40
   nl = nq
   goto 30
   40 continue
   !     now try 5
   50 continue
   if ( nl == 1 )return
   nq = nl / 5
   if ( nl /= nq * 5 )goto 60
   nl = nq
   goto 50
   60 continue
   result = 0
   return
end subroutine check_prime_factors 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute the ceil of x that is also a product of powers of 2,3,5.
!-------------------------------------------------------------------
!     --- COMPUTE_NFFT ---
!-------------------------------------------------------------------
!     this routine uses check_prime_factors to get the smallest
!     integer greater or equal than x which is decomposable into
!     powers of 2,3,5
!     not so efficient... bfd

subroutine compute_nfft(x,n)
   implicit none
   _REAL_ x
   integer n

   integer m,result,i
   m = int(x) - 1
   do i = 1,100
      m = m + 1
      call check_prime_factors(m,result)
      if ( result == 1 )then
         n = m
         return
      end if
   end do
   write(6,*)'compute_nfft: failed to get good fft array size'
   write(6,*)'x = ',x
   call mexit(6,1)
end subroutine compute_nfft 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Cubic spline interpolation.
!-------------------------------------------------------------------
!     --- CUBSPL ---
!-------------------------------------------------------------------
!     This code is from netlib.

subroutine cubspl ( tau, c, n, ibcbeg, ibcend )
   
   !     from  * a practical guide to splines *  by c. de boor
   
   !     ************************  input  ***************************
   !     n = number of data points. assumed to be .ge. 2.
   !     (tau(i), c(1,i), i=1,...,n) = abscissae and ordinates of the
   !        data points. tau is assumed to be strictly increasing.
   !     ibcbeg, ibcend = boundary condition indicators, and
   !     c(2,1), c(2,n) = boundary condition information. specifically,
   !        ibcbeg = 0  means no boundary condition at tau(1) is given.
   !           in this case, the not-a-knot condition is used, i.e. the
   !           jump in the third derivative across tau(2) is forced to
   !           zero, thus the first and the second cubic polynomial pieces
   !           are made to coincide.)
   !        ibcbeg = 1  means that the slope at tau(1) is made to equal
   !           c(2,1), supplied by input.
   !        ibcbeg = 2  means that the second derivative at tau(1) is
   !           made to equal c(2,1), supplied by input.
   !        ibcend = 0, 1, or 2 has analogous meaning concerning the
   !           boundary condition at tau(n), with the additional infor-
   !           mation taken from c(2,n).
   
   !     ***********************  output  **************************
   !     c(j,i), j=1,...,4; i=1,...,l (= n-1) = the polynomial coefficients
   !        of the cubic interpolating spline with interior knots (or
   !        joints) tau(2), ..., tau(n-1). precisely, in the interval
   !        (tau(i), tau(i+1)), the spline f is given by
   !           f(x) = c(1,i)+h*(c(2,i)+h*(c(3,i)+h*c(4,i)/3.)/2.)
   !        where h = x - tau(i). the function program *ppvalu* may be
   !        used to evaluate f or its derivatives from tau,c, l = n-1,
   !        and k=4.
   
   implicit none
   integer ibcbeg,ibcend,n,   i,j,l,m
   _REAL_ c(4,n),tau(n),   divdf1,divdf3,dtau,g
   
   !****** a tridiagonal linear system for the unknown slopes s(i) of
   !  f  at tau(i), i=1,...,n, is generated and then solved by gauss elim-
   !  ination, with s(i) ending up in c(2,i), all i.
   !     c(3,.) and c(4,.) are used initially for temporary storage.
   
   l = n-1
   
   !  compute first differences of tau sequence and store in c(3,.). also,
   !  compute first divided difference of data and store in c(4,.).
   
   do m=2,n
      c(3,m) = tau(m) - tau(m-1)
      c(4,m) = (c(1,m) - c(1,m-1))/c(3,m)
   end do
   
   !onstruct first equation from the boundary condition, of the form
   !             c(4,1)*s(1) + c(3,1)*s(2) = c(2,1)
   
   if (ibcbeg-1)                     11,15,16
   11 if (n > 2)                     goto 12
   
   !     no condition at left end and n = 2.
   
   c(4,1) = 1.d0
   c(3,1) = 1.d0
   c(2,1) = 2.d0*c(4,2)
   goto 25
   
   !     not-a-knot condition at left end and n .gt. 2.
   
   12 c(4,1) = c(3,3)
   c(3,1) = c(3,2) + c(3,3)
   c(2,1) =((c(3,2)+2.d0*c(3,1))*c(4,2)*c(3,3)+ &
         c(3,2)*c(3,2)*c(4,3))/c(3,1)
   goto 19
   
   !     slope prescribed at left end.
   
   15 c(4,1) = 1.d0
   c(3,1) = 0.d0
   goto 18
   
   !     second derivative prescribed at left end.
   
   16 c(4,1) = 2.d0
   c(3,1) = 1.d0
   c(2,1) = 3.d0*c(4,2) - c(3,2)/2.d0*c(2,1)
   18 if(n == 2)                      goto 25
   
   !  if there are interior knots, generate the corresp. equations and car-
   !  ry out the forward pass of gauss elimination, after which the m-th
   !  equation reads    c(4,m)*s(m) + c(3,m)*s(m+1) = c(2,m).
   
   19 continue
   do m=2,l
      g = -c(3,m+1)/c(4,m-1)
      c(2,m) = g*c(2,m-1) + 3.d0 * &
            (c(3,m)*c(4,m+1)+c(3,m+1)*c(4,m))
      c(4,m) = g*c(3,m-1) + 2.d0*(c(3,m) + c(3,m+1))
   end do
   
   ! construct last equation from the second boundary condition, of the form
   !           (-g*c(4,n-1))*s(n-1) + c(4,n)*s(n) = c(2,n)
   !     if slope is prescribed at right end, one can go directly to back-
   !     substitution, since c array happens to be set up just right for it
   !     at this point.
   
   if (ibcend-1)                     21,30,24
   21 if (n == 3 .and. ibcbeg == 0) goto 22
   
   !     not-a-knot and n .ge. 3, and either n.gt.3 or  also not-a-knot at
   !     left end point.
   
   g = c(3,n-1) + c(3,n)
   c(2,n) = ((c(3,n)+2.d0*g)*c(4,n)*c(3,n-1) &
         + c(3,n)*c(3,n)*(c(1,n-1)-c(1,n-2))/c(3,n-1))/g
   g = -g/c(4,n-1)
   c(4,n) = c(3,n-1)
   goto 29
   
   !     either (n=3 and not-a-knot also at left) or (n=2 and not not-a-
   !     knot at left end point).
   
   22 c(2,n) = 2.d0*c(4,n)
   c(4,n) = 1.d0
   goto 28
   
   !     second derivative prescribed at right endpoint.
   
   24 c(2,n) = 3.d0*c(4,n) + c(3,n)/2.d0*c(2,n)
   c(4,n) = 2.d0
   goto 28
   25 if (ibcend-1)                     26,30,24
   26 if (ibcbeg > 0)                goto 22
   
   !     not-a-knot at right endpoint and at left endpoint and n = 2.
   
   c(2,n) = c(4,n)
   goto 30
   28 g = -1.d0/c(4,n-1)
   
   ! complete forward pass of gauss elimination.
   
   29 c(4,n) = g*c(3,n-1) + c(4,n)
   c(2,n) = (g*c(2,n-1) + c(2,n))/c(4,n)
   
   ! carry out back substitution
   
   30 j = l
   40 c(2,j) = (c(2,j) - c(3,j)*c(2,j+1))/c(4,j)
   j = j - 1
   if (j > 0)                  goto 40
   
   !****** generate cubic coefficients in each interval, i.e., the deriv.s
   !  at its left endpoint, from value and slope at its endpoints.
   
   do i=2,n
      dtau = c(3,i)
      divdf1 = (c(1,i) - c(1,i-1))/dtau
      divdf3 = c(2,i-1) + c(2,i) - 2.d0*divdf1
      c(3,i-1) = 2.d0*(divdf1 - c(2,i-1) - divdf3)/dtau
      c(4,i-1) = (divdf3/dtau)*(6.d0/dtau)
   end do
   return
end subroutine cubspl 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Allocate memory for Ewald switch function lookup table.
!-------------------------------------------------------------------
!     --- EED_TABLE_MEM ---
!-------------------------------------------------------------------
!     Patterned after locmem in amber code.

subroutine eed_table_mem(startreal,endreal, &
      mxeedtab,leed_cub,leed_lin, &
      dxdr,cutoffnb,eedtbdns, eedmeth,verbose)

   implicit none
#  include "extra.h"
   integer verbose
   integer startreal,endreal
   integer mxeedtab,leed_cub,leed_lin,eedmeth
   _REAL_ dxdr,cutoffnb,eedtbdns

   integer mem_ptr

   !     eed table: assume all nonbond distances are less than 1.5 X cutoff
   !     between nonbond updates; i.e. no excess motion
   !     this is enforced by nbfilter

   mxeedtab = int(dxdr*eedtbdns*cutoffnb*1.5d0)
   if(master.and. verbose > 0) &
         write(6, '(a,6X,a,i10)') &
         '|','Size of EEDTABLE                 = ', mxeedtab

   !     do real array offsets
   mem_ptr = startreal

   !     permanent or heap REAL storage
   if ( eedmeth == 1 )then
      call adj_mem_ptr(mem_ptr,leed_lin,0)
      call adj_mem_ptr(mem_ptr,leed_cub,4*mxeedtab)
   else if ( eedmeth == 2 )then
      call adj_mem_ptr(mem_ptr,leed_lin,2*mxeedtab)
      call adj_mem_ptr(mem_ptr,leed_cub,0)
   else
      call adj_mem_ptr(mem_ptr,leed_lin,0)
      call adj_mem_ptr(mem_ptr,leed_cub,0)
   end if

   endreal =  mem_ptr

   return
end subroutine eed_table_mem 



!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine ewald_mem here]
!-------------------------------------------------------------------
!     --- EWALD_MEM ---
!-------------------------------------------------------------------
!     This routine is patterned after locmem in amber code,
!     in order to allow seamless integration.

subroutine ewald_mem(maxpr,natom_local,numbad_nnb, &
      startreal,endreal,startint,endint )
   
   !     ARGUMENT LIST:
   !     MAXPR (INPUT), the size of the pairlist array
   !       NOTE: currently the nonbond list is not packed into smaller words!
   !       This is because the nonbond routines use preimaging which may
   !       generate larger indices in imaged atoms, or in other words 16 bits
   !       may not suffice even if natom_local < 32K!!
   !     natom_local (INPUT) is the number of atoms in the system
   !     startreal (INPUT) is the starting offset in X of the ewald specific
   !       real memory  e.g. in amber many other things such as coords and
   !       forces are stored in X(1) thru X(startreal-1)
   !     endreal (OUTPUT) is computed as the end of that memory
   use nblist, only: reclng,dirlng,cutoffnb, &
                     nghb,maxnptrs, &
                     maximage, &
                     maxnblst,nblist_allocate

   implicit none
   integer natom_local,numbad_nnb
   integer startreal,endreal,startint,endint,maxpr
#  include "extra.h"
#  include "ew_pme_recip.h"
#  include "memory.h"
#  include "ew_erfc_spline.h"
#  include "ew_cntrl.h"
#  include "ew_mpole.h"
#  include "md.h"

   integer startr,starti
#ifdef MPI
#  include "parallel.h"
#else
   integer :: numtasks=1
#endif
#  include "ew_parallel.h"
   logical do_output

   do_output=master.and. (verbose >= 1)
   if (do_output) write(6, '(/,"|  ",a,/)') 'EWALD MEMORY USE:'
   
   !     --- FIRST GET SOME SIZES
   
   !     -- pme sizes
   
   call pmesh_kspace_get_sizes(natom_local,ncopy,verbose)
   
   !     -- adjust ewald code: mask size
   
   mxadjmsk = 2*numbad_nnb
   if (do_output) write(6, '(a,4x,a,i10)') &
         '|','Adjacent nonbond minimum mask    = ', mxadjmsk
   
   !     -- local nonbond code sizes
   
   call local_nb_get_sizes(numbad_nnb, &
         periodic,nogrdptrs,verbose, &
         natom_local,ncopy)

   
   !     --- NOW DO POINTER OFFSETS
   
   if (do_output) write(6, '("|  ",a)') &
         'EWALD LOCMEM POINTER OFFSETS'
   startr = startreal
   starti = startint
   ! point multipole stuff
   if ( mpoltype >= 1 )then
      call mpole_mem(natom_local,startr,endreal, &
            starti,endint, &
            ifirst,imiddle,ithird, &
            lfixdip,ldipole,linddip,lquad,lfield,ltorque, &
            leold1,leold2,leold3,ldipvel,indmeth)
   else
      linddip=1
      lfield=1
   end if
   starti = endint
   
   !     -- reals
   
   startr = endreal
   call pme_mem(natom_local,startr,endreal, &
         sizfftab,sizffwrk,siztheta,siz_q,sizscr, &
         nfft1,nfft2,nfft3,order, &
         lfftable,lprefac1,lprefac2,lprefac3)
   if (do_output) write(6, '(a,6X,a,i10)') &
         '|','Real memory needed by PME        = ', endreal-startr
   startr = endreal
   !     end if

   call eed_table_mem(startr,endreal, &
         mxeedtab,leed_cub,leed_lin,dxdr,cutoffnb,eedtbdns, &
         eedmeth,verbose)
   if (do_output) write(6, '(a,6X,a,i10)') &
         '|','Real memory needed by EEDTABLE   = ', endreal-startr
   startr = endreal
   
   !     -- integers
   
   endint = starti
   call adj_mem_ptr(endint,imask1,mxadjmsk)
   call adj_mem_ptr(endint,imask2,mxadjmsk)

   if (do_output) write(6, '(a,6X,a,i10)') &
         '|','Integer memory needed by ADJ     = ', endint-starti
   starti = endint
   
   !     --- THE LOCAL NONBOND NEEDS TO BE MAPPED LAST DUE TO NB LIST
   call nblist_allocate(natom_local,ntypes,num_direct,numtasks)
   
   call loc_nb_mem(natom_local, &
         startr,starti,endreal,endint,cutoffnb, &
         ntypes,periodic )
   if (do_output) write(6, '(a,6X,a,i10)') &
         '|','Integer memory used by local nonb= ', endint-starti
   if (do_output) write(6, '(a,6X,a,i10)') &
         '|','Real memory used by local nonb   = ', endreal-startr
   if (do_output) write(6, '(/,a,4X,a,i10)') &
         '|','MAX NONBOND PAIRS = ', maxpr
   maxnblst = maxpr
   return
end subroutine ewald_mem 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize Ewald calculations.
!-------------------------------------------------------------------
!     --- EW_STARTUP ---
!-------------------------------------------------------------------
!     called after load_ewald_info() and locmem() to fill in some
!     arrays and perform other initial chores.

subroutine ew_startup(natom_local,iblo,inb,x,ix)

   use nblist, only:cutoffnb,nvdwcls,mxlstmsk,fill_xtran,fill_tranvec
   use stack
   use qmmm_module, only : qmmm_nml, qmmm_struct
   implicit none
   character(kind=1,len=10) :: routine="ew_startup"
#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "ew_erfc_spline.h"
#  include "extra.h"
#  include "def_time.h"
#  include "md.h"
#  include "memory.h"

   integer natom_local,iblo(*),inb(*),l_tau
   _REAL_ x(*)
   integer ix(*)
   _REAL_ tim1,tim2
#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif
#endif
   integer ierr

   call fill_xtran()
   call fill_tranvec()
   if (periodic /= 0) then !Do not need this memory for non-periodic sim.
      call get_stack(l_tau,mxeedtab,routine)
      if(.not. rstack_ok)then
         deallocate(r_stack)
         allocate(r_stack(1:lastrst),stat=alloc_ier)
         call reassign_rstack(routine)
      endif
      REQUIRE(rstack_ok)
      
      call fill_eed_table(eedtbdns,mxeedtab, &
           x(leed_cub),r_stack(l_tau),x(leed_lin),eedmeth,ee_type)
      call free_stack(l_tau,routine)
   end if
   call init_coulomb_switch(cutoffnb,dxdr, &
         x(leed_cub),x(leed_lin),eedtbdns,eedmeth,ee_type)
   call vdw_correct_setup(natom_local,ix(i04),ntypes,nvdwcls)
   call load_adj_mask(iblo,inb,natom_local, &
         mxadjmsk,ix(imask1),ix(imask2),numadjst,verbose)
   if ( qmmm_nml%ifqnt .and. (qmmm_nml%qmmm_int == 5) ) then
      ! Mechanical QM/MM embedding
      ! => need to exclude also all QM-QM atom pairs
      ! NOTES: We assume natom_local to be the total number of atoms
      !        in the system
      call load_adj_mask_qmqm(iblo,inb,natom_local,mxadjmsk,ix(imask1),ix(imask2), &
           numadjst,qmmm_struct%atom_mask,verbose)
   end if
   call load_list_mask(iblo,inb,natom_local, &
         mxlstmsk)
#ifdef MPI
   if(i_do_recip)then
      call mpi_comm_size(recip_comm,numtasks,ierr)
      call mpi_comm_rank(recip_comm,mytaskid,ierr)
#endif
      if (periodic /= 0) then !do not need to do this for non-periodic
        call pmesh_kspace_setup( &
              x(lprefac1),x(lprefac2),x(lprefac3),x(lfftable), &
              nfft1,nfft2,nfft3,order,sizfftab,sizffwrk,opt_infl,ew_type)
      end if
#ifdef MPI
      call mpi_comm_size(world_comm,numtasks,ierr)
      call mpi_comm_rank(world_comm,mytaskid,ierr)
   end if
#endif
   return
end subroutine ew_startup 

!-------------------------------------------------------------------

!     --- FILL_EED_TABLE ---

!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fill_eed_table here]
subroutine fill_eed_table(eedtbdns,mxeedtab, &
      eed_cub,tau,eed_lin,eedmeth,ee_type)

   implicit none
   integer mxeedtab,eedmeth,ee_type
   _REAL_ eedtbdns,eed_cub(4,*),tau(*),eed_lin(2,*)

   _REAL_ del,x,switch,d_switch_dx
   integer i,ibcbeg,ibcend

   del = 1.d0 / eedtbdns

   if ( eedmeth >= 3 )return
   if ( eedmeth == 1 )then
      x = 0.d0
      call get_ee_func(x,switch,d_switch_dx,ee_type)
      eed_cub(2,1) = d_switch_dx
      x = (mxeedtab-1)*del
      call get_ee_func(x,switch,d_switch_dx,ee_type)
      eed_cub(2,mxeedtab) = d_switch_dx
      do i = 1,mxeedtab
         x = del*(i-1)
         call get_ee_func(x,switch,d_switch_dx,ee_type)
         tau(i) = x
         eed_cub(1,i) = switch
      end do
      ibcbeg = 1
      ibcend = 1
      call cubspl ( tau, eed_cub, mxeedtab, ibcbeg, ibcend )
   else if ( eedmeth == 2 )then
      do i = 1,mxeedtab
         x = del*(i-1)
         call get_ee_func(x,switch,d_switch_dx,ee_type)
         eed_lin(1,i) = switch
         eed_lin(2,i) = d_switch_dx
      end do
   end if
   return
end subroutine fill_eed_table 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine find_ewaldcof here]
!-------------------------------------------------------------------
!     --- FIND_EWALDCOF ---
!-------------------------------------------------------------------

subroutine find_ewaldcof(cutoff,dtol,ewaldcof)
   implicit none
   _REAL_ cutoff,dtol,ewaldcof
   integer i,n
   _REAL_ term,x,xlo,xhi,y,erfc

   ! first get direct sum tolerance. How big must ewaldcof be to get
   ! terms outside the cutoff below tol

   x = 0.5d0
   i = 0
   10 continue
   x = 2.d0 * x
   i = i + 1
   y = x * cutoff
   call erfcfun(y,erfc)
   term = erfc/cutoff
   if ( term >= dtol) goto 10

   ! binary search tolerance is 2 to the -50th

   n = i + 50
   xlo = 0.d0
   xhi = x
   do i = 1,n
      x = (xlo+xhi)/2
      y = x * cutoff
      call erfcfun(y,erfc)
      term = erfc/cutoff
      if ( term >= dtol )then
         xlo = x
      else
         xhi = x
      end if
   end do
   ewaldcof = x

   return
end subroutine find_ewaldcof 

!-------------------------------------------------------------------

!     --- LOAD_ADJ_MASK ---

!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine load_adj_mask here]
subroutine load_adj_mask(iblo,inb,natom, &
      mxadjmsk,mask1,mask2,numadjst,verbose)
#ifdef LES
   use pimd_vars, only: ipimd
#endif
   implicit none
   integer iblo(*),inb(*),natom
   integer mask1(*),mask2(*),numadjst,mxadjmsk

   integer ji,i,j,nx,k,verbose
#  include "extra.h"
#ifdef MPI
#  include "parallel.h"
#endif
#ifdef LES
#  include "les.h"
#endif
   ! PASS 1 check sizes

   ji = 0
   numadjst = 0
   do i = 1,natom-1
      nx = iblo(i)
      do j = 1,nx
         k = inb(ji+j)

         ! use extra arrays mask1,mask2 to speed the nb_adjust routine

         if ( k > i )then
            numadjst = numadjst + 1
         end if
      end do
      ji = ji + nx
   end do
   if (master .and. verbose > 0) &
         write(6, '(5x,a,i10)') 'Total number of mask terms = ', numadjst
   if ( numadjst > mxadjmsk) &
         call sander_bomb('load_mask','MXADJMSK not big enough!!',' ')

   ! PASS 2 fill mask array

   ji = 0
   numadjst = 0
   do i = 1,natom-1
      nx = iblo(i)
      do j = 1,nx
         k = inb(ji+j)

         !         use extra arrays mask1,mask2 to speed the nb_adjust routine
         if ( k > i )then
#ifdef LES
         if( ipimd==0 .or. lestyp(i).ne.lestyp(k) .or. cnum(i).eq.cnum(k) ) then
#endif             
             numadjst = numadjst + 1
             mask1(numadjst) = i
             mask2(numadjst) = k
         end if
#ifdef LES
         end if
#endif
      end do
      ji = ji + nx
   end do

   return
end subroutine load_adj_mask 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [load_adj_mask_qmqm: excludes QM-QM pairs for mechanical embedding]
subroutine load_adj_mask_qmqm(iblo,inb,natom, &
      mxadjmsk,mask1,mask2,numadjst,qm_mask,verbose)
#ifdef LES
   use pimd_vars, only: ipimd
#endif
   implicit none
   integer, intent(in) :: iblo(*), inb(*), natom, mxadjmsk, verbose
   logical, intent(in) :: qm_mask(*)
   integer, intent(inout) ::  mask1(*), mask2(*), numadjst

   integer :: i,j,k,ji,nx,ni,nexcl,numqmadjst
   logical :: excluded
#  include "extra.h"
#ifdef MPI
#  include "parallel.h"
#endif
#ifdef LES
#  include "les.h"
#endif
   ! PASS 1 check sizes

   ji=0
   numqmadjst = 0
   do i = 1, natom-1
      nx = iblo(i)
      if (qm_mask(i)) then
         do j = i+1, natom
            if (qm_mask(j)) then
               ! both atoms are QM
               ! check whether they have already been excluded
               ! because of being in the list of 1-2, 1-3 etc interactions
               excluded = .false.
               do ni = 1, nx
                  nexcl = inb(ji+ni)
                  if (nexcl == j) then
                     excluded = .true.
                  end if
               end do

               if (.not.excluded) then
                  ! atom pair has not yet been excluded
                  ! exclude it now
                  numqmadjst = numqmadjst + 1
               end if
            end if
         end do
      end if
      ji=ji+nx
   end do
   if (master .and. verbose > 0) then
      write(6, '(5x,a,i10)') 'Total number of QM mask terms = ', numqmadjst
      call flush(6)
   end if

   if ( numadjst + numqmadjst > mxadjmsk) &
         call sander_bomb('load_mask','MXADJMSK not big enough!!',' ')

   ! PASS 2 fill mask array
   ji=0
   do i = 1, natom-1
      nx = iblo(i)
      if (qm_mask(i)) then
         do j = i+1, natom
            if (qm_mask(j)) then
               ! both atoms are QM
               ! check whether they have already been excluded
               ! because of being in the list of 1-2, 1-3 etc interactions
               excluded = .false.
               do ni = 1, nx
                  nexcl = inb(ji+ni)
                  if (nexcl == j) then
                     excluded = .true.
                  end if
               end do

               if (.not.excluded) then
                  ! atom pair has not yet been excluded
                  ! exclude it now
#ifdef LES
                  if( ipimd==0 .or. lestyp(i).ne.lestyp(k) .or. cnum(i).eq.cnum(k) ) then
#endif             
                     numadjst = numadjst + 1
                     mask1(numadjst) = i
                     mask2(numadjst) = j
#ifdef LES
                  end if
#endif
               end if
            end if
         end do
      end if
      ji=ji+nx
   end do
 
end subroutine load_adj_mask_qmqm


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Ewald namelist read postprocessing.
!-------------------------------------------------------------------
!     --- LOAD_EWALD_INFO ---
!-------------------------------------------------------------------
!     Routine which sets up some sizes and parameters necessary
!     for locmem() and must be called prior to calling locmem().
!     Currently this is called inside mdread().

subroutine load_ewald_info(parm,inpcrd,ntp)
   use nblist, only: a, b, c, alpha, beta, gamma, &
                     skinnb, nbflag, nbfilter,cutoffnb,cutlist
   use amoeba_runmd, only: AM_RUNMD_get_ucell_info
   use binrestart, only: check_nc_restart, read_nc_restart_box
   use file_io_dat, only: MAX_FN_LEN
   implicit none
#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "ew_erfc_spline.h"
#  include "box.h"
#  include "ew_mpole.h"
#ifdef MPI
#  include "ew_parallel.h"
#endif
   _REAL_ ax,bx,cx,alphax,betax,gammax
   integer ipol,ntp
   character(len=MAX_FN_LEN) parm,inpcrd
   logical newstyle

!   induced = ipol
!   mpoltype= ipol
   
   ! get the values for ucell:
   
   if( ntb > 0 ) then
      ! Check for new Netcdf restart format
      if ( check_nc_restart(inpcrd) ) then
        write(6,*) 'getting box info from netcdf restart file'
        call read_nc_restart_box(inpcrd,ax,bx,cx,alphax,betax,gammax)
      else
         write(6,*)'getting new box info from bottom of inpcrd'
         call AMOEBA_check_newstyle_inpcrd(inpcrd,newstyle)
         if ( newstyle )then
            call AM_RUNMD_get_ucell_info(inpcrd,ax,bx,cx,alphax,betax,gammax)
         else
            call peek_ewald_inpcrd(inpcrd,ax,bx,cx,alphax,betax,gammax)
         end if
      endif

   end if
   call read_ewald(ax,bx,cx,alphax,betax,gammax)
   
   !     ---- also set the older parameters in box.h to these values:
   
   box(1) = a
   box(2) = b
   box(3) = c
   
   !  Check if non-isotropic scaling is requested with non-orthorhombic cell:
   
   if ( ntp > 1 ) then
      if ( abs(alpha-90.0d0 ) > 1.d-5 .or. &
            abs(beta - 90.0d0) > 1.d-5 .or. &
            abs(gamma - 90.0d0) > 1.d-5 ) then
         call sander_bomb('read_ewald', &
               'Cannot do non-isotropic scaling unless orthorhombic' &
               ,'use ntp=1 if angles are not 90 degrees. ')
      end if
   end if

   ! implement extended list
   
   if ( nbflag == 1 .or. skinnb > 1.d-10 ) then
      nbfilter = cutoffnb
      cutlist = cutoffnb + skinnb
   else
      
      ! don't exceed the erfc table maximum
      
      cutlist = cutoffnb
      nbfilter = 1.4 * cutoffnb
   end if
   
   ! setup dxdr map from r to x in table lookup of eed
   
   if ( ee_type == 1 )then
      dxdr = ew_coeff
   else if ( ee_type == 2 )then
      dxdr = 1.d0/cutoffnb
   end if

   return
end subroutine load_ewald_info 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!-------------------------------------------------------------------
!     --- LOAD_EWALD_POL_INFO ---
!-------------------------------------------------------------------
!     Routine which sets up multipole-related parameters.
!
subroutine load_ewald_pol_info(ipol, pol, pol2, df, natom, iok)
   implicit none
#  include "ew_cntrl.h"
#  include "ew_mpole.h"
   integer ipol, i, natom, iok
   _REAL_ pol(*), pol2(*), df(*), sixth, small
   
   sixth = 1.d0/6.d0
   small = 1.d-6

!   induced = ipol

   ! setup point multipoles
   
!   mpoltype = induced
   if ( induced == 0 )then
      irstdip =0
      indmeth = 0
   end if
!
   if ( ipol > 1 ) then
   if ( dipdamp <= 0.d0 ) then
           if ( ipol == 2 ) then
         dipdamp = 0.572d0
      else if ( ipol == 3 ) then
         dipdamp = 2.089d0
      else if ( ipol == 4 ) then
         dipdamp = 1.662d0
      end if
      write(6,'(5X,a/5X,a,f12.4/)') &
     '|  Thole coefficient dipdamp was not found in ewald control. ', &
     '|  Use default value = ',dipdamp
   else
      write(6,'(5X,a,f12.4/)') &
     '|  Thole coefficient dipdamp read from ewald control = ',dipdamp
   end if

      if ( iok .eq. 0 ) then
         write(6,'(5X,A)') &
         '| Thole damping coefficients read from prmtop.'
      else
         write(6,'(5X,A/A,F12.5)') &
         '| Thole damping coefficients are unavailable in prmtop.', &
         '| Use default value dipdamp = ',dipdamp
      end if
      do i = 1, natom
         if (df(i) .le. small ) df(i) = dipdamp
      end do

   end if
!
   if ( ipol == 2 ) then
! YD pre-calculate the factors for performance
! Thole exponential
       do i = 1,natom
         pol2(i) = sqrt(pol(i)/df(i))
       end do
   else if ( ipol == 3) then
! Thole exponential
       do i = 1,natom
         pol2(i) = pol(i)**sixth/sqrt(df(i))
       end do
   else if ( ipol == 4) then
! Thole-linear
       do i = 1,natom
         pol2(i) = pol(i)**sixth*sqrt(df(i))
       end do
   end if

   return
end subroutine load_ewald_pol_info 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine load_list_mask here]
!-------------------------------------------------------------------
!     --- LOAD_LIST_MASK ---
!-------------------------------------------------------------------

subroutine load_list_mask(iblo,inb,natom, &
      mxlstmsk)
   use nblist, only:lstmask,nummask,maskptr
   implicit none
   integer iblo(*),inb(*),natom
   integer mxlstmsk

   ! double the mask to deal with our list generator

   integer ji,i,j,m,nx,off,k,tot

#  include "extra.h"
#  include "ew_cntrl.h"
#ifdef MPI
#  include "parallel.h"
#endif

   ! PASS 1 get pointers, check size

   ji = 0
   do i = 1,natom
      nummask(i) = 0
   end do
   do i = 1,natom - 1
      nx = iblo(i)
      do j = 1,nx
         k = inb(ji+j)
         if ( k > 0 )then
            nummask(k) = nummask(k) + 1
            nummask(i) = nummask(i) + 1
         end if
      end do
      ji = ji + nx
   end do
   tot = 0
   do i = 1,natom
      tot = tot + nummask(i)
   end do
   if (master .and. verbose > 0) &
         write(6, '(5x,a,i10)') 'Total number of mask terms = ', tot
   if ( tot > mxlstmsk) &
         call sander_bomb('load_mask','MXLSTMSK not big enough!!',' ')
   off = 0
   do i = 1,natom
      maskptr(i) = off
      off = off + nummask(i)
   end do

   ! PASS 2 fill mask array

   ji = 0
   do i = 1,natom
      nummask(i) = 0
   end do
   do i = 1,natom-1
      nx = iblo(i)
      do j = 1,nx
         k = inb(ji+j)
         if ( k > 0 )then
            nummask(k) = nummask(k) + 1
            m = maskptr(k) + nummask(k)
            lstmask(m) = i
            nummask(i) = nummask(i) + 1
            m = maskptr(i) + nummask(i)
            lstmask(m) = k
         end if
      end do
      ji = ji + nx
   end do

   return
end subroutine load_list_mask 

!-------------------------------------------------------------------

!     --- LOC_NB_MEM ---

!-------------------------------------------------------------------
!     ...patterned after locmem in amber code...


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine loc_nb_mem here]
subroutine loc_nb_mem(natom, &
      startreal,startint,endreal,endint,cutoffnb, &
      ntypes, periodic )
   implicit none
   _REAL_ cutoffnb
   integer natom,startreal,endreal,startint,endint
   ! INTEGER POINTERS
   integer ntypes, periodic 
         
#  include "extra.h"
#  include "ew_parallel.h"
#ifdef MPI
#  include "parallel.h"
#else
   integer numtasks,mytaskid
   parameter (numtasks=1,mytaskid=0)
#endif

   integer mem_ptr


   ! do integer array  offsets
#ifndef MPI
   num_direct=1
#endif

   mem_ptr = startint
   endint = mem_ptr

   return
end subroutine loc_nb_mem 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine print_sz here]
subroutine print_sz(label,isize)
#include "extra.h"
   integer isize
   character(len=*) label
   if(master)then
      write(6, 100)isize,label
      100 format("Ewald loc_nb_mem allocating: ",i9,a)
   end if
   return
end subroutine print_sz 

!-------------------------------------------------------------------

!     --- LOCAL_NB_GET_SIZES ---


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine local_nb_get_sizes here]
subroutine local_nb_get_sizes(numbad_nnb,periodic,nogrdptrs,verbose,natom,ncopy)

   use nblist, only: setup_grids,setup_grid_sizes,nghb,maxnptrs, &
                     nucgrd1,nucgrd2,nucgrd3,nucgmax, &
                     maximage,mxlstmsk
   implicit none
   integer numbad_nnb,periodic,verbose, &
         natom, ncopy

   _REAL_ density
   integer nucgmin
   logical nogrdptrs
#  include "extra.h"
#ifdef MPI
#  include "parallel.h"
#endif
   
   maxnptrs = ((2*nghb+1)*(2*nghb+1)+1)/2
   if (master.and.verbose > 0) &
         write(6, '(a,4X,a,i10)') &
         '|','Max number of pointers           = ', maxnptrs
   mxlstmsk = 4*numbad_nnb
   if (master.and.verbose > 0) &
         write(6, '(a,4X,a,i10)') &
         '|','List build mxlstmsk               = ', mxlstmsk

   !     --- use task 0 so only sizes are calculated for memory allocation
   call setup_grids(periodic,nogrdptrs,verbose)

   ! set upper bounds for ucell grid dimensions to allow for volume fluctuations
   
   nucgmax = nucgrd1*nucgrd2*nucgrd3 * 1.67
   nucgmin = nucgrd1*nucgrd2*nucgrd3 / 1.67

   if(periodic == 0)then
      nucgmax=max(nucgmax,1000)
      nucgmin=max(1,nucgmin)
   end if
   
   ! density is average number of atoms per subcell
   
   density = natom
   density = density / (nucgmin)
   
   ! the maximum image is density times nimgrd1*nimgrd2*nimgrd3
   ! times a fudge factor of 1.1
   
   maximage = density * nucgrd1*nucgrd2*nucgrd3 * 1.1
   if (master.and.verbose > 0) &
         write(6, '(a,4X,a,i10,/)') &
         '|','Maximage  = ', maximage
   return
end subroutine local_nb_get_sizes 

!-------------------------------------------------------------------
!    --- mpole_mem
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine mpole_mem here]
subroutine mpole_mem(natom,startreal,endreal, &
      startint,endint, &
      ifirst,imiddle,ithird, &
      lfixdip,ldipole,linddip,lquad,lfield,ltorque, &
      leold1,leold2,leold3,ldipvel,indmeth)
   implicit none
   integer natom
   integer startreal,endreal
   integer startint,endint
   integer ifirst,imiddle,ithird, &
         lfixdip,ldipole,linddip,lquad,lfield,ltorque, &
         leold1,leold2,leold3,ldipvel,indmeth
   integer mem_ptr

   mem_ptr = startreal
   ! comment out the fixed dipole and quad code
   !     call adj_mem_ptr(mem_ptr,lfixdip,3*natom)
   !     call adj_mem_ptr(mem_ptr,ldipole,3*natom)
   call adj_mem_ptr(mem_ptr,linddip,3*natom)
   !     call adj_mem_ptr(mem_ptr,lquad,9*natom)
   call adj_mem_ptr(mem_ptr,lfield,3*natom)
   ! initialize to default beginning array index
   ! to enable compiler based array bounds checking
   leold1 = 1
   leold2 = 1
   leold3 = 1
   ldipvel = 1
   if ( indmeth >= 0 .and. indmeth <= 2 )then
      call adj_mem_ptr(mem_ptr,leold1,3*natom)
   end if
   if ( indmeth >= 1 .and. indmeth <= 2 )then
      call adj_mem_ptr(mem_ptr,leold2,3*natom)
   end if
   if ( indmeth == 2 )then
      call adj_mem_ptr(mem_ptr,leold3,3*natom)
   end if
   if ( indmeth == 3 )then
      call adj_mem_ptr(mem_ptr,ldipvel,3*natom)
   end if

   ! comment out the fixed dipole and quad code
   !     call adj_mem_ptr(mem_ptr,ltorque,3*natom)
   endreal =  mem_ptr
   mem_ptr = startint
   ! comment out the fixed dipole and quad code
   !     call adj_mem_ptr(mem_ptr,ifirst,natom)
   !     call adj_mem_ptr(mem_ptr,imiddle,natom)
   !     call adj_mem_ptr(mem_ptr,ithird,natom)
   endint = mem_ptr
   return
end subroutine mpole_mem 
!-------------------------------------------------------------------
!     --- PME_MEM ---

!     ...patterned after locmem in amber code...


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pme_mem here]
subroutine pme_mem(natom,startreal,endreal, &
      sizfftab,sizffwrk,siztheta,siz_q,sizscr, &
      nfft1,nfft2,nfft3,order, &
      lfftable,lprefac1,lprefac2,lprefac3)

   implicit none
   integer startreal,endreal
   integer nfft1,nfft2,nfft3,natom,order, &
         sizfftab,sizffwrk,siztheta,siz_q,sizscr, &
         lfftable,lprefac1,lprefac2,lprefac3

   integer mem_ptr
   
   mem_ptr = startreal

   ! permanent or heap REAL storage
   call adj_mem_ptr(mem_ptr,lfftable,sizfftab)
   call adj_mem_ptr(mem_ptr,lprefac1,nfft1)
   call adj_mem_ptr(mem_ptr,lprefac2,nfft2)
   call adj_mem_ptr(mem_ptr,lprefac3,nfft3)

   endreal =  mem_ptr

   return
end subroutine pme_mem 

!-------------------------------------------------------------------
!     --- PMESH_KSPACE_GET_SIZES ---
!     ...this routine computes the parameters needed for heap or
!     stack allocation specified below...


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pmesh_kspace_get_sizes here]
subroutine pmesh_kspace_get_sizes(natom,ncopy,verbose)
   implicit none
   !      integer nfft1,nfft2,nfft3,natom,order,
   !     $     sizfftab,sizffwrk,siztheta,siz_Q,sizheap,sizstack,sizscr
   integer natom,ncopy,verbose

#  include "ew_pme_recip.h"

   ! INPUT
   !      nfft1,nfft2,nfft3,natom,order
   !      nfft1,nfft2,nfft3 are the dimensions of the charge grid array
   !      natom is number of atoms
   !      order is the order of B-spline interpolation

   ! OUTPUT
   !      sizfftab,sizffwrk,siztheta,siz_Q
   !      sizfftab is permanent 3d fft table storage
   !      sizffwrk is temporary 3d fft work storage
   !      siztheta is size of arrays theta1-3 dtheta1-3
   !      sizscr is size of temporary arrays needed by vector code
   !      sizheap is total size of permanent storage
   !      sizstack is total size of temporary storage

   integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork
#  include "extra.h"
#ifdef MPI
#  include "parallel.h"
#  include "ew_parallel.h"
#endif



   sizscr = 1
   call get_fftdims(nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork, &
         sizfftab,sizffwrk)
#ifndef MPI
   siztheta = natom*order
   siz_q = 2*nfftdim1*nfftdim2*nfftdim3
#else
   siztheta = natom*order*(nxyslab(0)+6)/nfft3
   siz_q = max( ntxyslab*(nxyslab(0)), ntxzslab*nxzslab(0))
#endif
   sizheap = nfft1+nfft2+nfft3+sizfftab
   sizstack = siz_q+6*siztheta+sizffwrk+3*natom+sizscr

   if (master.and. verbose >= 1) &
         write(6, '(a,4X,a,i10)') &
         '|','Total heap storage needed        = ', sizheap
   
   return
end subroutine pmesh_kspace_get_sizes 

!-------------------------------------------------------------------
!     --- PMESH_KSPACE_SETUP ---

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine pmesh_kspace_setup here]
subroutine pmesh_kspace_setup( &
      prefac1,prefac2,prefac3,fftable, &
      nfft1,nfft2,nfft3,order,sizfftab,sizffwrk,opt_infl,ew_type)

  use ew_bspline
#ifdef MPI
  use fft,only:fft_init,column_fft_flag
#endif
   implicit none

   !  see DO_PMESH_KSPACE for explanation of arguments

   integer nfft1,nfft2,nfft3,order,sizfftab,sizffwrk,opt_infl, &
         ew_type
   _REAL_ prefac1(*),prefac2(*),prefac3(*)
   _REAL_ fftable(sizfftab)
   _REAL_ ffwork_dummy

   _REAL_ dummy
   integer nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw

   call get_fftdims(nfft1,nfft2,nfft3, &
         nfftdim1,nfftdim2,nfftdim3,nfftable,nffwork,sfft,sffw)
   call load_prefacs(prefac1,prefac2,prefac3, &
         nfft1,nfft2,nfft3,order,opt_infl)
   call fft_setup(dummy,fftable,ffwork_dummy, &
         nfft1,nfft2,nfft3,nfftdim1,nfftdim2,nfftdim3, &
         nfftable,nffwork)
#ifdef MPI
   if(column_fft_flag)then
      call fft_init(nfft1,nfft2,nfft3)
   endif
#endif
   
   return
end subroutine pmesh_kspace_setup 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Read box info from inpcrd.  Abort if box info is not found.
subroutine peek_ewald_inpcrd(inpcrd,a,b,c,alpha,beta,gamma)
   use file_io_dat, only : MAX_FN_LEN, INPCRD_UNIT, MAX_LINE_BUF_LEN
   use sander_lib, only  : check_inpcrd_overflow
   implicit none
   character(len=MAX_FN_LEN) inpcrd
   character(len=MAX_LINE_BUF_LEN) line
   _REAL_ a,b,c,alpha,beta,gamma

   integer natom
   _REAL_ tt,x1,x2,x3,x4,x5,x6
   integer i,ic,nrow,nextra,justcrd,vel
   
   call amopen(INPCRD_UNIT,inpcrd,'O','F','R')
   read(INPCRD_UNIT,'(a)') line
   
   read(INPCRD_UNIT,'(a)') line
   if( line(6:6) == ' ' ) then ! this is an old, i5 file
     read(line,'(i5,e15.7)', err=666) natom,tt
   elseif( line(7:7) == ' ' ) then ! sander 7/8/9/10 large system format...
     read(line,'(i6,e15.7)', err=666) natom,tt
   elseif( line(8:8) == ' ' ) then ! Sander 11 - 1 mil+ format
     read(line,'(i7,e15.7)', err=666) natom,tt
   else                   ! assume amber 11 VERY large system format. 10 mil+
     read(line,'(i8,e15.7)', err=666) natom,tt
   end if

   if ( natom <= 2 )then
      write(6,*)'peek_ewald_inpcrd: ', &
            ' Cannot Deduce box info from inpcrd. Too few atoms'
      call mexit(6,1)
      return
   end if
   ic = 0
   do i = 1,9999999
      read(INPCRD_UNIT,9028,end=81,err=667)x1,x2,x3,x4,x5,x6
      ic = ic+1
   end do
   81 continue
   close(INPCRD_UNIT)
   nrow = natom/2
   nextra = mod(natom,2)
   justcrd = nrow+nextra
   vel = 2*justcrd
   if ( ic == justcrd .or. ic == vel )then
      write(6,'(a)') &
            '| peek_ewald_inpcrd: Box info not found in inpcrd'
      call mexit(6,1)
      return
   end if
   if ( ic == justcrd+1 .or. ic == vel+1 )then
      write(6,'(a)') '| peek_ewald_inpcrd: Box info found'
      a = x1
      b = x2
      c = x3
      if ( x4 > 0.d0 .and. x5 == 0.d0 .and. x6 == 0.d0)then
         
         !         ---only has beta
         alpha = 90.d0
         beta = x4
         gamma = 90.d0
         
      else if (x4 == 0.d0 .and. x5 == 0.d0 .and. x6 == 0.d0 ) then
         
         !         --- no angles in input: assume they are all 90:
         alpha = 90.d0
         beta  = 90.d0
         gamma = 90.d0

      else

         !         --- found the angles:
         alpha = x4
         beta = x5
         gamma = x6
      end if
      return
   end if
   
   write(6,*) 'peek_ewald_inpcrd: SHOULD NOT BE HERE'
   call mexit(6,1)
   
   9028 format(6f12.7)
   return

666 write(6, '(a)') 'ERROR: I could not find the number of atoms or the time on'
    write(6, '(3a)') '       the second line of your inpcrd file [', trim(inpcrd), &
                    ']. Bad INPCRD file!'
    close(INPCRD_UNIT)
    call mexit(6,1)
667 write(6, '(2a)') 'ERROR: Problem reading coordinates or velocities from ', &
                    trim(inpcrd)
    write(6, '()')
    write(6, '(a,i5,a)') 'I could not understand line ', ic + 3, ' :'

! I want to print the offending line. However, there's no convenient way to
! do this in Fortran. So I will rewind the whole file, read ic+2 lines, then
! print the ic+3'th line to unit 6
    rewind(INPCRD_UNIT)
    do i = 1, ic + 2
      read(INPCRD_UNIT, '(a80)') line
    end do

    ! Now we're there. Read the line and print it
    read(INPCRD_UNIT, '(a80)') line
    write(6, '(a)')  line
    write(6, '()')
    call check_inpcrd_overflow(line, .true.) ! .true. because this is PBC
    close(INPCRD_UNIT)

    call mexit(6,1)

end subroutine peek_ewald_inpcrd 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initialize to default values and read Ewald namelist.
!-----------------------------------------------------------------------
!     --- READ_EWALD ---
!-----------------------------------------------------------------------
!         called from load_ewald_info from mdread

subroutine read_ewald(ax,bx,cx,alphax,betax,gammax)
   use qmmm_module , only : qmmm_struct
   use constants, only : RETIRED_INPUT_OPTION, one, half
   use nblist, only: amod=>a, bmod=>b, cmod=>c, &
                     alphamod=>alpha, betamod=>beta, gammamod=>gamma, &
                     skinnbmod=>skinnb, nbflagmod=>nbflag, nbtellmod=>nbtell, &
                     nbfilter,cutoffnb, &
                     ucell,recip,dirlng,reclng,sphere,volume
   use fft,only: column_fft_flag
   use constants, only: NO_INPUT_VALUE
   use file_io_dat
                     
   implicit none
   _REAL_ a,b,c,alpha,beta,gamma,skinnb
   integer nbflag,nbtell
   _REAL_ ax,bx,cx,alphax,betax,gammax
   _REAL_ y,erfc,eigmin
   integer maxmlim

#  include "ew_cntrl.h"
#  include "ew_pme_recip.h"
#  include "ew_erfc_spline.h"
#  include "ew_legal.h"
#  include "ew_mpole.h"
#  include "box.h"
#  include "md.h"
#ifdef MPI
#  include "parallel.h"
#endif
   integer gridpointers,column_fft
   
   !  Amber 8 (and earlier) codes allowed the user to enter a,b,c,alpha,
   !    beta, gamma, but there is no ordinary reason to do so, and it can
   !    be dangerous.  If you feel a need for this flexibility , uncomment
   !    out the following line.  Otherwise, these will not be legal namelist
   !    variables.
   !    namelist /ewald/a,b,c,alpha,beta,gamma,dsum_tol,ew_coeff, &

   namelist /ewald/dsum_tol,ew_coeff, &
         skinnb,diptol,dipmass,diptau, &
         nfft1,nfft2,nfft3,order,opt_infl, &
         ischrgd,verbose,nbflag,nbtell,netfrc, &
         ew_type,vdwmeth,eedmeth,ee_type, &
         eedtbdns,rsum_tol,maxexp,mlimit,use_pme, &
         maxiter,indmeth,irstdip,nquench, &
         frameon,chngmask,scaldip, &
         gridpointers,column_fft
   
   !  ---- Determine if nonperiodic system, set big box to satisfy sanity checks.
   !  ---- Box will be determined after &ewald is read.
   if(ntb == 0)then
      if( igb == 0 .and. ipb == 0 ) then
         a=90.d0
         b=90.d0
         c=90.d0
         alpha = 90.d0
         beta = 90.d0
         gamma = 90.d0
         use_pme = 0
         eedmeth = 4
         periodic=0
         write(6,'(a,a)' ) "|    NONPERIODIC ", &
               '  ntb=0 and igb=0: Setting up nonperiodic simulation'
      end if
   else
      !  ---- we already got the unitcell params from inpcrd if ntx=7
      !  ---- set the default namelist values before reading in
      a = ax
      b = bx
      c = cx
      alpha = alphax
      beta = betax
      gamma = gammax
      use_pme = 1
      eedmeth = 1
      periodic=1
   end if
   
   cutoffnb = cut
   dsum_tol = 1.d-5
   ew_coeff = 0.d0
   diptol = 1.d-4
   dipmass = 0.33d0
   diptau = 11.0d0
   maxiter = 20
   nquench = 0
   indmeth = 3
   irstdip = 0
   if( imin == 0 .or. ntb > 0 ) then
      nbflag = 1
      skinnb = 2.0d0
   else
      nbflag = 0
      skinnb = 0.0d0
   end if
   ischrgd = RETIRED_INPUT_OPTION
   nbtell = 0
   nfft1 = 0
   nfft2 = 0
   nfft3 = 0
   order = 4
   opt_infl = 1
   verbose = 0
   netfrc = NO_INPUT_VALUE
   ew_type = 0
   vdwmeth = 1
   ee_type = 1
   eedtbdns = 5000.d0
   rsum_tol = 5.d-5
   maxexp = 0.d0
   mlimit(1) = 0
   mlimit(2) = 0
   mlimit(3) = 0
   frameon = 1
   chngmask = 1
   scaldip = 1
   gridpointers=1
   nogrdptrs=.false.
   column_fft = 0
! M-WJ
   dipdamp = -1.0d0
   
   if ( mdin_ewald ) then
      rewind 5
      read(5,nml=ewald,err=668)
   end if
   
   if ( ischrgd /= RETIRED_INPUT_OPTION ) then
      write(6,'(/,a,/,a,/,a)') 'Warning: ischrgd has been retired.', &
            'A nonzero net charge due to roundoff error is always neutralized.'
   end if

   !     --- set temp. regulation of dipoles based on diptau:
   
   if( diptau < 10d0) then
      nttdip = 1
   else
      nttdip = 0
   end if
   
   !  ---- adjust the netfrc value, different for minimization than dynamics:
   if (netfrc .eq. NO_INPUT_VALUE) then
      if (imin .eq. 0) then             
         netfrc = 1         
      else        
         netfrc = 0
      end if      
   else
      if (imin .ne. 0 .and. netfrc .ne. 0) then
         write(6,'(a)') '| Setting netfrc to 0 for minimization'
         netfrc = 0                                           
      end if
   end if

   !   ---- for nonperiodic user can request no grid pointers method
   !   ----   of list build by setting gridpointers=0 in &ewald
   nogrdptrs=gridpointers /= 1
   
   ! Assume unit cell params are complete
   !     get some related quantities:
   
   !    For non-periodic we need to make dummy cell parameters
   !      We need a cell that is bigger than the system, so
   !      the coords need to be read and analyzed. Allocate
   !      space for coord on the real stack.
   
   if(periodic == 0)then
      !        --- For nonperiodic systems, figure out the box size
      !        --- so that proper allocation of memory can be performed.
      call firstbox()
      a=amod
      b=bmod
      c=cmod
   end if
   
   ! Assume unit cell params are complete
   !     get some related quantities:
   

   call get_ucell(a,b,c,alpha,beta,gamma, &
         ucell,recip,dirlng,reclng,sphere,volume)
   
   ! Assign ewald coefficient;
   ! necessary even when the reciprocal Ewald is not computed.
   ! ew_coeff is sometimes used to get dxdr; dxdr determines mxeedtab.
   
   if ( ew_coeff < 1.d-6 )then
      call float_legal_range('dsum_tol: (direct sum tol) ', &
            dsum_tol,tollo,tolhi)
      call find_ewaldcof(cutoffnb,dsum_tol,ew_coeff)
   else
      call float_legal_range('ew_coeff: (ewald conv. coeff.) ', &
            ew_coeff,ew_coefflo,ew_coeffhi)
      y = ew_coeff*cutoffnb
      call erfcfun(y,erfc)
      dsum_tol = erfc/cutoffnb
   end if
   
   ! grid sizes
   
   if ( nfft1 == 0 ) then
      !        RCFFT needs an even dimension for x direction
      call compute_nfft((a + one)*half   ,nfft1)
      nfft1=nfft1*2
   end if
   if ( nfft2 == 0 ) call compute_nfft(b,nfft2)
   if ( nfft3 == 0 ) call compute_nfft(c,nfft3)
   call int_legal_range('nfft1: (grid size) ', &
         nfft1,gridlo,gridhi)
   call int_legal_range('nfft2: (grid size) ', &
         nfft2,gridlo,gridhi)
   call int_legal_range('nfft3: (grid size) ', &
         nfft3,gridlo,gridhi)
   call test_prime_factors('nfft1',nfft1)
   call test_prime_factors('nfft2',nfft2)
   call test_prime_factors('nfft3',nfft3)
   if ( ew_type /= 0 )then
      write(6,*)'-----------------------------------------'
      write(6,*)'====== Running Regular Ewald code'
      maxmlim = max(mlimit(1),mlimit(2),mlimit(3))
      if ( maxmlim > 0 )then
         call maxexp_from_mlim(maxexp,mlimit,reclng,recip)
         write(6,100)maxexp
         100 format(1x,'maxexp calculated from mlimit: ',e8.3)
      else
         call float_legal_range('rsum_tol: (Ewald recip sum tol) ', &
               rsum_tol,tollo,tolhi)
         if ( maxexp < 1.d-6 )then
            call find_maxexp(ew_coeff,rsum_tol,maxexp)
            write(6,101)maxexp
            101 format(1x,'maxexp calculated from rsum_tol: ',e8.3)
         end if
         
         ! eigmin typically bigger than this (unless badly distorted cell)
         
         eigmin = 0.5d0
         call get_mlim(maxexp,mlimit,eigmin,reclng,recip)
      end if
   end if
   if ( vdwmeth > 2 )then
      write(6,*)'pme vdw not supported yet'
      call mexit(6,1)
   end if
   if ( ee_type /= 1 )then
      write(6,*)'only switch supported is erfc'
      call mexit(6,1)
   end if
#ifndef MPI
   if(column_fft == 1)then
      write(6,'("| Column fft is only for parallel, setting to false")')
      column_fft = 0
   endif
#endif
   if(column_fft == 1) then
      if ( ntp > 1 ) then
         if ( abs(alpha-90.0d0 ) > 1.d-5 .or. &
              abs(beta - 90.0d0) > 1.d-5 .or. &
              abs(gamma - 90.0d0) > 1.d-5 ) then
            write(6,'(a,i3/a/a)') &
                 "| Column fft: ifbox is :",ifbox, &
                 "| Column fft: only orthorhombic periodic boxes only at present.", & 
                 "|             setting to false"
            column_fft = 0
         endif
      endif
   endif
   column_fft_flag = (column_fft == 1)
   if(column_fft_flag ) &
        write(6,'(a/a)') &
        "| Column fft: orthorhombic periodic box, using COLUMN_FFT.", & 
        "|             "
   
   amod=a
   bmod=b
   cmod=c
   alphamod=alpha
   betamod=beta
   gammamod=gamma
   skinnbmod=skinnb
   nbtellmod=nbtell
   nbflagmod=nbflag

   return

668 write(6, '(a)') 'ERROR: Problem reading &ewald section of MDIN ', trim(mdin)
end subroutine read_ewald 
!-------------------------------------------------------------------
!     --- VDW_CORRECT_SETUP ----
!   sets up the numbers of atoms in each vdw type. used for
!   analytic pressure, energy, correction to vdw dispersion


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine vdw_correct_setup here]
subroutine vdw_correct_setup(natom,iac,ntypes,nvdwclas)
   implicit none
   integer natom,iac(*),nvdwclas(*),ntypes

   integer n,j
   do n = 1,ntypes
      nvdwclas(n) = 0
   end do
   do n = 1,natom
      j = iac(n)
      nvdwclas(j) = nvdwclas(j) + 1
   end do
   
   
   return
end subroutine vdw_correct_setup 
!--------------------------------------------------------------------
!     --- VDW_CORRECTION ----
! get analytic estimate of energy and virial corrections
! due to dispersion interactions beyond the cutoff
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine vdw_correction here]
subroutine vdw_correction(ico,ntypes,nvdwclas, &
      volume,evdwr,rec_vird,cn2,cutoffnb)
   use constants, only : TWOPI
   implicit none
   integer ico(*),ntypes,nvdwclas(*)
   _REAL_ volume,evdwr,rec_vird(3,3),cn2(*),cutoffnb

   integer i,j,ic,iaci,index
   _REAL_ term,pi,prefac
   term = 0.d0
   prefac = TWOPI/(3.d0*volume*cutoffnb*cutoffnb*cutoffnb)
   do i = 1,ntypes
      iaci = ntypes*(i-1)
      do j = 1,ntypes
         index = iaci + j
         ic = ico(index)
         if ( ic > 0 )then
            term = term + nvdwclas(i)*nvdwclas(j)*cn2(ic)
         end if
      end do
   end do
   evdwr = -prefac*term
   do j = 1,3
      do i = 1,3
         rec_vird(i,j) = 0.d0
      end do
      rec_vird(j,j) = -2.d0*evdwr
   end do

   return
end subroutine vdw_correction 

!-------------------------------------------------------------------

! -- ti decomp
subroutine vdwdec_correction(natom,ntypes,iac,ico,nvdwclas,cn2,volume,cutoffnb)
   use constants, only: TWOPI
   use decomp,    only: decpair
   use file_io_dat

   implicit none
#  include "md.h"
   integer i,j,natom,ntypes,iac(*),ico(*),nvdwclas(*)
   _REAL_  cn2(*),volume,cutoffnb

   integer iaci,ic
   _REAL_ prefac,term

   prefac = -TWOPI/(3.d0*volume*(cutoffnb**3))

   do i = 1, natom
      iaci = ntypes*(iac(i)-1)
      if(idecomp < 3) then
         do j = 1, ntypes
            ic = ico(iaci + j)
            if(ic > 0) then
               term = prefac*nvdwclas(j)*cn2(ic)
               call decpair(3,i,i,term/(nstlim/ntpr))
            endif
         end do ! j = 1, ntypes
      else
         do j = 1, natom
            ic = ico(iaci + iac(j))
            if(ic > 0) then
               term = prefac*cn2(ic)
               call decpair(3,i,j,term/(nstlim/ntpr))
            endif
         end do ! j = 1, natom
      endif
   end do ! i = 1, natom

   return
end subroutine vdwdec_correction
! -- ti decomp end

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine zero_array here]
subroutine zero_array(array,num)
   implicit none
   _REAL_ array(*)
   integer num

   integer i
   do i = 1,num
      array(i) = 0.d0
   end do
   return
end subroutine zero_array 
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine array_copy here]
subroutine array_copy(a,b,num)
   implicit none
   _REAL_ a(*),b(*)
   integer num

   integer i
   do i = 1,num
      b(i) = a(i)
   end do
   return
end subroutine array_copy 

#ifdef MPI
!-------------------------------------------------------------------
!     --- STARTUP_GROUPS ---
!-------------------------------------------------------------------

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine startup_groups here]
subroutine startup_groups(err)
   implicit none
   integer i,j,err
#  include "parallel.h"
#  include "ew_parallel.h"
#  include "extra.h"
#  include "memory.h"

#ifdef MPI_DOUBLE_PRECISION
#undef MPI_DOUBLE_PRECISION
#endif
   include 'mpif.h'
   integer ierr
#ifdef CRAY_PVP
#define MPI_DOUBLE_PRECISION MPI_REAL8
#endif


   call mpi_bcast(num_recip,1,mpi_integer,0,commsander,ierr)
# ifndef EWGRPS
   if(num_recip < numtasks .and. master) &
         write(6,*) "Groups not compiled, cannot use nrecip"
   world_comm=commsander
   recip_comm=world_comm
   num_recip = numtasks
   i_do_recip = .true.

   direct_comm=world_comm
   i_do_direct= .true.
   num_direct = numtasks

   do i=1,numtasks
      ranks(i)=i-1
   end do
# else
   world_comm=commsander
   call mpi_comm_group(commsander,world_group,ierr)

   if(num_recip > numtasks)then
      write(6,*) num_recip, &
            'NUM_RECIP exceeds number of available PEs',numtasks
      num_recip = numtasks
   end if

   do i=1,num_recip
      ranks(i)=i-1
   end do
   if(mytaskid < num_recip)then
      i_do_recip = .true.
   else
      i_do_recip = .false.
   end if
   call mpi_group_incl(world_group,num_recip,ranks,recip_group,err)
   call mpi_barrier(commsander,ierr)
   call mpi_comm_create(commsander,recip_group,recip_comm,err)
   call mpi_barrier(commsander,ierr)

   if(num_recip == numtasks)then
      i_do_direct = .true.
      num_direct = numtasks
      do i=1,numtasks
         ranks(i)=i-1
      end do
   else
      i_do_direct = .false.
      num_direct = numtasks-num_recip
      do i=1,numtasks-num_recip
         ranks(i)=num_recip-1+i
         if(mytaskid == num_recip-1+i)i_do_direct = .true.
      end do
   end if
   call mpi_group_incl(world_group,num_direct,ranks,direct_group,err)
   call mpi_barrier(commsander,ierr)
   call mpi_comm_create(commsander,direct_group,direct_comm,err)
   call mpi_barrier(commsander,ierr)

# endif
   return
end subroutine startup_groups 

#endif
!-------------------------------------------------

