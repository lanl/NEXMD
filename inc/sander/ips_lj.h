
! epilogue: 12-6 LF IPS terms  by Xiongwu Wu

do im_new = 1,icount
   j = cache_bckptr(im_new)

   dfee = cache_df(im_new)
   delx = cache_x(im_new)
   dely = cache_y(im_new)
   delz = cache_z(im_new)
   delr2inv = cache_r2(im_new)

   ic = ico(iaci+iac(j))
   f6 = cn2(ic)*rips6r
   f12 = cn1(ic)*rips12r
!  L-J r6 term
!   etr6=1/r6+a0+r2*(a1+r2*(a2+a3*r2))
!   detr6/dr*r1=-6/r6+r2*(d1+r2*(d2+d3*r2))
!
        UIPS2R=(DELR2INV*RIPS2)
        UIPS2=ONE/UIPS2R
        UIPS4=UIPS2*UIPS2
        UIPS6R=UIPS2R*UIPS2R*UIPS2R
        UIPS12R=UIPS6R*UIPS6R
            PVC=UIPS6R+AIPSVC(0)+UIPS2*(AIPSVC(1)+UIPS2*(AIPSVC(2)+UIPS2*AIPSVC(3)))
            DVCU=-SIX*UIPS6R+UIPS2*(BIPSVC(1)+UIPS2*(BIPSVC(2)+UIPS2*BIPSVC(3)))
!  L-J r12 term 
!   etr12=1/r12+a0+r2*(a1+r2*(a2+a3*r2))
!   detr12/dr*r1=-12/r12+r4*(d1+r4*(d2+d3*r4))
!
            PVA=UIPS12R+AIPSVA(0)+UIPS4*(AIPSVA(1)+UIPS4*(AIPSVA(2)+UIPS4*AIPSVA(3)))
            DVAU=-TWELVE*UIPS12R+UIPS4*(BIPSVA(1)+UIPS4*(BIPSVA(2)+UIPS4*BIPSVA(3) ))
            evdw = evdw +F12*(PVA-PIPSVAC)-F6*(PVC-PIPSVCC)
            df = dfee -(F12*DVAU-F6*DVCU)*DELR2INV
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
