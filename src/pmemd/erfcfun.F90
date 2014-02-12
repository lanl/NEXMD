!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Compute the complementary error function.
!-----------------------------------------------------------------------
!     --- ERFCFUN ---
!---------------------------------------------------------------------
!  The algorithm is ??

subroutine derfcfun(x,erfc)

   implicit none
   double precision, intent(in)  :: x
   double precision, intent(out) :: erfc
   double precision :: absx, c, p, q, nonexperfc, erf

   absx=abs(x)
   if (x > 26.d0) then
      erfc = 0.d0

   else if (x < -5.5d0) then
      erfc = 2.0d0

   else if (absx <= 0.5d0) then
      c = x * x
      p=((-0.356098437018154d-1*c+0.699638348861914d1)*c+   &
            0.219792616182942d2)*c+0.242667955230532d3
      q=((c+0.150827976304078d2)*c+0.911649054045149d2)*c+  &
            0.215058875869861d3
      erf = x*p/q
      erfc = 1.0-erf

   else if (absx < 4.0) then
      c=absx
      p=((((((-0.136864857382717d-6*c+0.564195517478974d0)*c+    &
         0.721175825088309d1)*c+0.431622272220567d2)*c+    &
         0.152989285046940d3)*c+0.339320816734344d3)*c+    &
         0.451918953711873d3)*c+0.300459261020162d3
      q=((((((c+0.127827273196294d2)*c+0.770001529352295d2)*c+    &
         0.277585444743988d3)*c+0.638980264465631d3)*c+    &
         0.931354094850610d3)*c+0.790950925327898d3)*c+    &
         0.300459260956983d3
      if ( x > 0.d0 ) then
         nonexperfc = p/q
      else
         nonexperfc = 2.d0*exp(x*x) - p/q
      end if
      erfc = exp(-absx*absx)*nonexperfc
      if (x < 0.d0) erfc = 2.0- erfc

   else
      c=1.0/(x*x)
      p=(((0.223192459734185d-1*c+0.278661308609648d0)*c+      &
         0.226956593539687d0)*c+0.494730910623251d-1)*c+      &
         0.299610707703542d-2
      q=(((c+0.198733201817135d1)*c+0.105167510706793d1)*c+      &
         0.191308926107830d0)*c+0.106209230528468d-1
      c=(-c*p/q + 0.564189583547756d0)/absx
      if( x > 0.d0 ) then
         nonexperfc = c
      else
         nonexperfc = 2.d0*exp(x*x) - c
      end if
      erfc = exp(-absx*absx)*nonexperfc
      if (x < 0.d0) erfc = 2.0- erfc
   end if

   return
   
end subroutine derfcfun 
