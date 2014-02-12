
! epilogue: 12-10 LJ terms.


!  --- this code allows 10-12 terms; in many (most?) (all?) cases, the
!       only "nominal" 10-12 terms are on waters, where the asol and bsol
!       parameters are always zero; hence we can skip this part.

do im_new = 1,icount
   j = cache_bckptr(im_new)

   df =   cache_df(im_new)
   delx = cache_x(im_new)
   dely = cache_y(im_new)
   delz = cache_z(im_new)

#ifdef HAS_10_12
   delr2inv = cache_r2(im_new)
   ic = -ico(iaci+iac(j))
   r10 = delr2inv*delr2inv*delr2inv*delr2inv*delr2inv
#  ifdef LES
   lfac=lesfac(lestmp+lestyp(j))
   f10 = bsol(ic)*r10*lfac
   f12 = asol(ic)*(r10*delr2inv)*lfac
#  else
   f10 = bsol(ic)*r10
   f12 = asol(ic)*(r10*delr2inv)
#  endif
   ehb = ehb + f12 - f10
   ! -- ti decomp: add hbond-terms to vdw-energy in energy decomposition
   if(decpr .and. idecomp > 0) call decpair(3,i,j,(f12 - f10)/(nstlim/ntpr))
   df = df + (12.d0*f12 - 10.d0*f10)*delr2inv
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
