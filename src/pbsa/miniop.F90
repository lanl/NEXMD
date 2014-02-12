
#define _REAL_ double precision

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine distan here]
subroutine distan(x1,y1,z1,x2,y2,z2, dis)
 ! implicit none
   implicit none
  _REAL_ x1,y1,z1,x2,y2,z2, dis
   dis = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)
   dis = sqrt(dis)

   return
end subroutine distan 



!-------print out a vector x


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine privec here]
subroutine privec(n,x)
  !implicit none
   implicit none
   _REAL_  x(n)
   integer i,n
   do i=1,n
      write(*,11) x(i)
      11 format(6x, f14.6)
   end do
   write(*,*) '    '
   return
end subroutine privec 



!-------print out the first 8 columns of a matrix a

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine primat here]
subroutine primat(l,m,n,a)
  !implicit none
   implicit none
  integer l,m,n,mtemp,i,j
   _REAL_ a(l,m,n)

   mtemp = 8
   if (mtemp > n) mtemp = n
   do i=1,m
      write(*,11) (a(i,j,6), j=1,mtemp)
      11 format(6x, 8f14.6)
   end do
   write(*,*) '    '
   return
end subroutine primat 




!-------fmax finds the larger between a and b


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function fmax here]
function fmax(a,b)
   implicit none
   _REAL_ a,b,fmax
   fmax=a
   if (a < b) fmax=b
   return
end function fmax 



!-------fmin finds the smaller between a and b


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function fmin here]
function fmin(a,b)
   implicit none
  _REAL_ a,b,fmin
   fmin=a
   if (a > b) fmin=b
   return
end function fmin 



!-------sort a list in increasing order


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine sortin here]
subroutine sortin(n,x,io)
   implicit none
   _REAL_  x(n),io(n),temp,itemp
   integer n,i,j

   io(1) = 1
   do i=2, n
      io(i)=i
      do j=i,2,-1
         if (x(j) <= x(j-1)) then
            temp = x(j)
            x(j) = x(j-1)
            x(j-1) = temp
            itemp = io(j)
            io(j) = io(j-1)
            io(j-1) = itemp
         end if
      end do
   end do

   return
end subroutine sortin 




!-------cpyvec copys the vector a into b.


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cpyvec here]
subroutine cpyvec(n,a,b)
   implicit none
   integer i,n
   _REAL_ a(n),b(n)

   do  i=1,n
      b(i) = a(i)
   end do

   return
end subroutine cpyvec 



!-------cpymatcopys the matrix a into b


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine cpymat here]
subroutine cpymat(m,n,a,b)
   implicit none
   integer m,n,i,j
   _REAL_ a(m,n),b(m,n)

   do i=1,m
      do j=1,n
         b(i,j) = a(i,j)
      end do
   end do

   return
end subroutine cpymat 



!------- couputes the dot product of two vectors: a and b


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function dotpro here]
function dotpro(n,a,b)
   implicit none
   integer n,i
   _REAL_ a(n),b(n),dotpro

   dotpro = 0.0
   do i=1,n
      dotpro = dotpro + a(i)*b(i)
   end do

   return
end function dotpro 



!-------cosin couputes cos of the angle between vectors: a and b


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of function cosin here]
function cosin(n,a,b)
   implicit none
   integer n,i
   _REAL_ a(n),b(n),prod0,prod1,prod2,cosin

   prod0 = 0.0
   prod1 = 0.0
   prod2 = 0.0
   do i=1,n
      prod0 = prod0 + a(i)*b(i)
      prod1 = prod1 + a(i)*a(i)
      prod2 = prod2 + b(i)*b(i)
   end do

   cosin = sqrt(prod0)/sqrt(prod1)/sqrt(prod2)
   return
end function cosin 




!-------matvec finds the product a*x of a matrix a and a vector x


!!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!+ [Enter a one-line description of subroutine matvec here]
!subroutine matvec(m,n,a,x,y)
!   implicit none
!   _REAL_ a(m,n), x(n),y(m)
!   do i=1,m
!      y(i) = 0.0
!      do j=1,n
!         y(i) = y(i) + a(i,j)*x(j)
!      end do
!   end do
!
!   return
!end subroutine matvec 
!
!
!!-------matmat finds the product a*b of two matrices a and b


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine matmat here]
subroutine matmat(l,m,n,a,b,c)
   implicit none
   integer l,m,n,i,j,k
   _REAL_ a(l,m), b(m,n), c(l,n)

   do i=1,l
      do j=1,n
         c(i,j)=0.0
         do k=1,m
            c(i,j)=c(i,j)+a(i,k)*b(k,j)
         end do
      end do
   end do
   return
end subroutine matmat 

!-------transp finds the transpose of a 2D matrix


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine transp here]
subroutine transp(m,n,a,b)
   implicit none
   integer m,n,i,j
   _REAL_ a(m,n), b(n,m)
   do i=1,m
      do j=1,n
         b(j,i)=a(i,j)
      end do
   end do
   return
end subroutine transp 






!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine rootp2 here]
subroutine rootp2(a0,a1,a2,r1,r2,info)
   implicit none
   _REAL_ a0,a1,a2,r1,r2,t
   integer info
   
   !       **************************************************************
   !       *                                                            *
   !       *  rootp2 find the two roots of the quadratic  equations:    *
   !       *             a2*x^2+a1*x+a0 = 0                             *
   !       *                                                            *
   !       *  input:  a2, a1, a0                                        *
   !       *                                                            *
   !       *  output: r1, r2: the two roots if they exist               *
   !       *          info: =  0   if a0=a1=a2=0                        *
   !       *                = -1   if a1=a2=0, a1!=0                    *
   !       *                =  1   if a2=0 or one double root           *
   !       *                =  2   if two distinct real roots           *
   !       *                = -2   if two different complex roots       *
   !       *                                                            *
   !       **************************************************************

   if (abs(a2) <= 1.0d-10) then
      if (abs(a1) <= 1.0d-10) then
         if (abs(a0) <= 1.0d-10) then
            info = 0
         else
            info = -1
         end if
         return
      else
         info = 1
         r1 = -a0/a1
         r2 = r1
         return
      end if
   else
      t = a1*a1 - 4.0*a0*a2
      if (t >= 0.0) then
         info = 2
         t = sqrt(t)
         if (a1 >= 0.0) then
            r1 = ( -a1 - t)/(2.0*a2)
         else
            r1 = ( -a1 + t)/(2.0*a2)
         end if

         if (t == 0.0) then
            info = 1
            r2 = r1
         else
            r2 = a0/(a2*r1)
         end if
         return
      else
         info = -2
         return
      end if
   end if  ! (abs(a2) <= 1.0d-10)

end subroutine rootp2 

