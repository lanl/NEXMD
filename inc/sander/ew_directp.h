
! prologue: gather the data and put it in temporary arrays.

   n=ipairs(m)
   itran=ishft(n,-27)
   n = iand(n,mask27)
   j = bckptr(n)
   delx = imagcrds(1,n) + xktran(1,itran)
   dely = imagcrds(2,n) + xktran(2,itran)
   delz = imagcrds(3,n) + xktran(3,itran)
   delr2 = delx*delx + dely*dely+delz*delz
   
   if ( delr2 < filter_cut2 )then
      icount = icount + 1
      cache_x(icount) = delx
      cache_y(icount) = dely
      cache_z(icount) = delz
      cache_r2(icount) = delr2
      cache_bckptr(icount) = j
   end if
