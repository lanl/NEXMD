
! epilogue: 12-6 LF terms

do im_new = 1,icount
   j = cache_bckptr(im_new)

   dfee = cache_df(im_new)
   delx = cache_x(im_new)
   dely = cache_y(im_new)
   delz = cache_z(im_new)
   delr2inv = cache_r2(im_new)

   ic = ico(iaci+iac(j))
   r6 = delr2inv*delr2inv*delr2inv
#ifdef LES 
   lfac=lesfac(lestmp+lestyp(j))
   f6 = cn2(ic)*r6*lfac
   f12 = cn1(ic)*(r6*r6)*lfac

   if(ipimd>0) then
      if(cnum(i).eq.0.and.cnum(j).eq.0) then
         nrg_all(1:nbead)=nrg_all(1:nbead) + (f12-f6)/nbead
      else
         if(cnum(i).ne.0) then
            nrg_all(cnum(i)) = nrg_all(cnum(i)) + f12-f6
         else
            nrg_all(cnum(j)) = nrg_all(cnum(j)) + f12-f6
         endif
      endif
   endif
#else
#  ifdef TVDW
   !         "truncated" LJ repulsion:
   r6pinv = 1.d0/(1.d0/r6 + 0.01d0)
   f12 = cn1(ic)*r6pinv*r6pinv
   f6 = cn2(ic)*r6pinv
#  else
   if ( vdwmodel == 0 ) then
      f6 = cn2(ic)*r6
      f12 = cn1(ic)*(r6*r6)
   else
      f6 = cn5(ic)*r6
      f12 = cn4(ic)/exp(cn3(ic)/sqrt(delr2inv))
   endif
#  endif
#endif
   ! -- ti decomp
   if(decpr .and. idecomp > 0) call decpair(3,i,j,(f12 - f6)/(nstlim/ntpr))
   evdw = evdw + f12 - f6
#ifdef TVDW
   r4 = 1.d0/(delr2inv*delr2inv)
   df = dfee + (12.d0*f12 - 6.d0*f6)*r4*r6pinv
#else
   df = dfee + (12.d0*f12 - 6.d0*f6)*delr2inv
#endif

   dfx = delx*df
   dfy = dely*df
   dfz = delz*df
#ifndef noVIRIAL
   vxx = vxx - dfx*delx
   vxy = vxy - dfx*dely
   vxz = vxz - dfx*delz
   vyy = vyy - dfy*dely
   vyz = vyz - dfy*delz
   vzz = vzz - dfz*delz
#endif
   dumx = dumx + dfx
   dumy = dumy + dfy
   dumz = dumz + dfz
   force(1,j) = force(1,j) + dfx
   force(2,j) = force(2,j) + dfy
   force(3,j) = force(3,j) + dfz
end do  !  im_new = 1,icount
