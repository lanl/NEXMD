#define _REAL_ double precision

! The IIM method solves the singular free PB equation.
! bcopt = 4 or 5: inside is total field and outside is total field.
! bcopt = 6 or 7: inside is reaction field and outside is total field.
! bcopt = 8 or 9: inside is reaction field and outside is reaction field.
!
! the right side of equation inside
!
function ff_in(x,y,z)

   implicit none

   _REAL_ ff_in
   _REAL_ x,y,z

   ff_in = 0.d0


end function ff_in
!
! the right side of equation outside
!
function ff_out(x,y,z)

   implicit none

   _REAL_ ff_out
   _REAL_ x,y,z

   ff_out = 0.d0


end function ff_out
!
! the jump condition of potential at projection point (x,y,z)
! bcopt = 4 or 5: it is zero.
! bcopt = 6 or 7: it is the coulomb potential.
! bcopt = 8 or 9: it is zero.
!
function fw(bcopt,bi,bo,atmfirst,atmlast,x,y,z)

   use poisson_boltzmann, only : acrg,acrd
   implicit none

   _REAL_ :: bi,bo
   !common /para/bi,bo

   integer bcopt,atmfirst,atmlast
   _REAL_ fw
   _REAL_ x,y,z

   integer k
   _REAL_ r(1:3),d
   _REAL_ phi0

   fw = 0.d0
   if ( bcopt == 6 .or. bcopt == 7 ) then
         if(x==-2.0d0 .and. y==0.0d0 .and. z==0.0d0) then
         end if
         
      do k = atmfirst, atmlast
         r(1) = x - acrd(1,k)
         r(2) = y - acrd(2,k)
         r(3) = z - acrd(3,k)
         d = sqrt(sum(r*r))
         if(x==-2.0d0 .and. y==0.0d0 .and. z==0.0d0) then

         end if
         fw = fw + acrg(k) / d
      end do
   end if


end function fw
!
! the jump condition of field (b*E_n) at projection point (x,y,z)
! bcopt = 4 or 5: it is zero. 
! bcopt = 6 or 7: it is  bi       * coulomb field.
! bcopt = 8 or 9: it is (bi - bo) * coulomb field.
!
function fq(bcopt,bi,bo,atmfirst,atmlast,x,y,z,t)

   use poisson_boltzmann, only : acrg,acrd
   use iim_use
   implicit none

   _REAL_ :: bi,bo
   !common /para/bi,bo

   integer bcopt, atmfirst, atmlast
   _REAL_ fq
   _REAL_ x,y,z,t(3,3)

   integer k
   _REAL_ r(1:3),rn(1:3),d,dn,cst,d2,f(3),ft(3)

   fq = 0.d0
   if ( bcopt > 5 .and. bcopt < 10 )  then
      do k = atmfirst, atmlast
         r(1) = x - acrd(1,k)
         r(2) = y - acrd(2,k)
         r(3) = z - acrd(3,k)
         d2 = sum(r*r)
         d = sqrt(d2)
         r = r / d

         f = - acrg(k) / d2 * r
         call matvec(3,3,t,f,ft)
         fq = fq + ft(1)
!        rn(1) = x - xctr(1)
!        rn(2) = y - xctr(2)
!        rn(3) = z - xctr(3)
!        dn = sum(rn*rn)
!        cst = sum(r*rn)/sqrt(d*dn)
!        fq = fq - ecg(k)/d*cst
      end do
      if ( bcopt == 6 .or. bcopt == 7 ) fq =  bi      *fq
      if ( bcopt == 8 .or. bcopt == 9 ) fq = (bi - bo)*fq
   end if


end function fq
!
! the dielectric constant inside
!
function fb_in(bi,bo,x,y,z)

   implicit none
      
   _REAL_ :: bi,bo
   !common /para/bi,bo
         
   _REAL_ fb_in
   _REAL_ x,y,z

   fb_in = bi

end function fb_in
!
! the dielectric constant outside
!
function fb_out(bi,bo,x,y,z)

   implicit none

   _REAL_ :: bi,bo
   !common /para/bi,bo

   _REAL_ fb_out
   _REAL_ x,y,z

   fb_out = bo

end function fb_out
!
! the sigma term, i.e. the linear Boltzmann term
! set to zero for poisson equation right now
!
! inside term
!
function fk_in(x,y,z)
   implicit none

   _REAL_ fk_in,x,y,z

   fk_in = 0.0d0


end function fk_in
!
! outside term
!
function fk_out(x,y,z)
   implicit none

   _REAL_ fk_out,x,y,z

   fk_out = 0.0d0


end function fk_out
