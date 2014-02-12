
#define _REAL_ double precision
! ----- find jump in f(x,y,z)


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fjmps here]
subroutine fjmps(x1,y1,z1, fjmp)
   implicit none
   _REAL_ fjmp,x1,y1,z1,ff_out,ff_in

   fjmp  = ff_out(x1,y1,z1)-ff_in(x1,y1,z1)

   return
end subroutine fjmps 


! ----- find jump in k(x,y,z)


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine fkjmps here]
subroutine fkjmps(x1,y1,z1, fkjmp)
   implicit none
 _REAL_ x1,y1,z1,fkjmp,fk_out,fk_in

   fkjmp = fk_out(x1,y1,z1) - fk_in(x1,y1,z1)

   return
end subroutine fkjmps 


! ----- find beta_in, beta_out at the projection (x1, y1, z1)


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine betas here]
subroutine betas(hx,hy,hz,x1,y1,z1,t,bl_in,bl_out,b_in,b_out)
  use iim_use
   implicit none
   _REAL_ hx,hy,hz,x1,y1,z1,b_in,b_out,fb_in,fb_out

   _REAL_ t(3,3)
   _REAL_ bg_in(3),bg_out(3),bl_in(3),bl_out(3)

   !b_in  = fb_in(x1,y1,z1)
   !b_out = fb_out(x1,y1,z1)
   
    
   bg_in(1)  = (fb_in(b_in,b_out,x1+hx,y1,z1)-fb_in(b_in,b_out,x1-hx,y1,z1))/(2.0*hx)
   bg_in(2)  = (fb_in(b_in,b_out,x1,y1+hy,z1)-fb_in(b_in,b_out,x1,y1-hy,z1))/(2.0*hy)
   bg_in(3)  = (fb_in(b_in,b_out,x1,y1,z1+hz)-fb_in(b_in,b_out,x1,y1,z1-hz))/(2.0*hz)
   bg_out(1) = (fb_out(b_in,b_out,x1+hx,y1,z1)-fb_out(b_in,b_out,x1-hx,y1,z1))/(2.0*hx)
   bg_out(2) = (fb_out(b_in,b_out,x1,y1+hy,z1)-fb_out(b_in,b_out,x1,y1-hy,z1))/(2.0*hy)
   bg_out(3) = (fb_out(b_in,b_out,x1,y1,z1+hz)-fb_out(b_in,b_out,x1,y1,z1-hz))/(2.0*hz)

   call matvec(3,3,t,bg_in, bl_in)         !# derivative in local coordinates
   call matvec(3,3,t,bg_out,bl_out)

   return
end subroutine betas 

