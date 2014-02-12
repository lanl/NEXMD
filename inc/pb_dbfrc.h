
      ! collect contact forces

      if ( iatm > 0 ) then
         cnx(iatm) = cnx(iatm) + dum(1)
         cny(iatm) = cny(iatm) + dum(2)
         cnz(iatm) = cnz(iatm) + dum(3)

      ! collect reentry forces

      else 
         if ( triopt == 1 ) then
            if ( -iatm > narcdot-ntri ) then
!              write(58,'(a,3f20.15,i,3i5,f)') "H",crd(1),crd(2),crd(3)
               do i = 1, ntri
                  if ( -iatm == narcdot-ntri+i ) then
                     matm = triatm(1,i)
                     natm = triatm(2,i)
                     patm = triatm(3,i)
                     exit
                  end if
               end do
               mvec(1:3) = acrd(1:3,matm) - x(1:3)
               nvec(1:3) = acrd(1:3,natm) - x(1:3)
               pvec(1:3) = acrd(1:3,patm) - x(1:3)
               mdist = sqrt(mvec(1)**2 + mvec(2)**2 + mvec(3)**2)
               ndist = sqrt(nvec(1)**2 + nvec(2)**2 + nvec(3)**2)
               pdist = sqrt(pvec(1)**2 + pvec(2)**2 + pvec(3)**2)
               mvec = mvec/mdist
               nvec = nvec/ndist
               pvec = pvec/pdist

               m = 3; n = 3;
               a(1:3,1) = mvec(1:3)
               a(1:3,2) = nvec(1:3)
               a(1:3,3) = pvec(1:3)
               b(1:3) = dum(1:3)
               u = a;
               call svdcmp(u,m,n,mp,np,w,v)
               wmax = ZERO
               do j = 1, n
                 if ( w(j) > wmax ) wmax = w(j)
               end do
               thresh = TOL * wmax
               do j = 1, n
                 if ( w(j) < thresh ) w(j) = ZERO
               end do
               call svbksb(u,w,v,m,n,mp,np,b,t)
               rnx(matm) = rnx(matm) + t(1)*mvec(1)
               rny(matm) = rny(matm) + t(1)*mvec(2)
               rnz(matm) = rnz(matm) + t(1)*mvec(3)
               rnx(natm) = rnx(natm) + t(2)*nvec(1)
               rny(natm) = rny(natm) + t(2)*nvec(2)
               rnz(natm) = rnz(natm) + t(2)*nvec(3)
               rnx(patm) = rnx(patm) + t(3)*pvec(1)
               rny(patm) = rny(patm) + t(3)*pvec(2)
               rnz(patm) = rnz(patm) + t(3)*pvec(3)
            else
               iarc = dotarc(-iatm)
               matm = arcatm(1,iarc); natm = arcatm(2,iarc)
               mvec(1:3) = acrd(1:3,matm) - x(1:3)
               nvec(1:3) = acrd(1:3,natm) - x(1:3)
               mdist = sqrt(mvec(1)**2 + mvec(2)**2 + mvec(3)**2)
               ndist = sqrt(nvec(1)**2 + nvec(2)**2 + nvec(3)**2)
               mvec = mvec/mdist
               nvec = nvec/ndist

               m = 3; n = 2
               a(1:3,1) = mvec(1:3)
               a(1:3,2) = nvec(1:3)
               b(1:3) = dum(1:3)
               u = a;
               call svdcmp(u,m,n,mp,np,w,v)
               wmax = ZERO
               do j = 1, n
                 if ( w(j) > wmax ) wmax = w(j)
               end do
               thresh = TOL * wmax
               do j = 1, n
                 if ( w(j) < thresh ) w(j) = ZERO
               end do
               call svbksb(u,w,v,m,n,mp,np,b,t)
               rnx(matm) = rnx(matm) + t(1)*mvec(1)
               rny(matm) = rny(matm) + t(1)*mvec(2)
               rnz(matm) = rnz(matm) + t(1)*mvec(3)
               rnx(natm) = rnx(natm) + t(2)*nvec(1)
               rny(natm) = rny(natm) + t(2)*nvec(2)
               rnz(natm) = rnz(natm) + t(2)*nvec(3)
            end if
         else
            iarc = dotarc(-iatm)
            matm = arcatm(1,iarc); natm = arcatm(2,iarc)

            mvec(1:3) = x(1:3) - acrd(1:3,matm)
            nvec(1:3) = x(1:3) - acrd(1:3,natm)
            mvec = mvec/sqrt(mvec(1)**2 + mvec(2)**2 + mvec(3)**2)
            nvec = nvec/sqrt(nvec(1)**2 + nvec(2)**2 + nvec(3)**2)

            mxnv(1) = mvec(2)*nvec(3) - nvec(2)*mvec(3)
            mxnv(2) = nvec(1)*mvec(3) - mvec(1)*nvec(3)
            mxnv(3) = mvec(1)*nvec(2) - nvec(1)*mvec(2)
            mxnv = mxnv/sqrt(mxnv(1)**2 + mxnv(2)**2 + mxnv(3)**2)

            ! split dum() into tangent and normal directions wrt the plan of mvec/nvec

            dumnorm = dum(1)*mxnv(1) + dum(2)*mxnv(2) + dum(3)*mxnv(3)
            dum_norm = dumnorm*mxnv; dum_tang = dum - dum_norm

            ! further split dum_tangent into mvec and nvec directions

            mdotn = mvec(1)*nvec(1) + mvec(2)*nvec(2) + mvec(3)*nvec(3)
            rmdotn2 = ONE/(ONE - mdotn**2)
            fdotm = dum_tang(1)*mvec(1) + dum_tang(2)*mvec(2) + dum_tang(3)*mvec(3)
            fdotn = dum_tang(1)*nvec(1) + dum_tang(2)*nvec(2) + dum_tang(3)*nvec(3)
            if ( fdotm < ZERO .and. fdotn < ZERO) then
               mvec = -mvec; nvec = -nvec
            else if ( fdotm < ZERO ) then
               mvec = -mvec
               mdotn = -mdotn
            else if ( fdotn < ZERO ) then
               nvec = -nvec
               mdotn = -mdotn
            end if
            fdotm = abs(fdotm); fdotn = abs(fdotn)

            dfm = (fdotm - fdotn*mdotn)*rmdotn2
            dfn = (fdotn - fdotm*mdotn)*rmdotn2
            rnx(matm) = rnx(matm) + dfm*mvec(1) + HALF*dum_norm(1)
            rny(matm) = rny(matm) + dfm*mvec(2) + HALF*dum_norm(2)
            rnz(matm) = rnz(matm) + dfm*mvec(3) + HALF*dum_norm(3)
            rnx(natm) = rnx(natm) + dfn*nvec(1) + HALF*dum_norm(1)
            rny(natm) = rny(natm) + dfn*nvec(2) + HALF*dum_norm(2)
            rnz(natm) = rnz(natm) + dfn*nvec(3) + HALF*dum_norm(3)
         end if
      end if
