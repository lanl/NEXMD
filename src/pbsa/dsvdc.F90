#define _REAL_ double precision

      subroutine dsvdc(x,ldx,n,p,s,e,u,ldu,v,ldv,work,job,info)

      implicit none

      integer ldx,n,p,ldu,ldv,job,info
      _REAL_ x(ldx,ldv),s(ldx+1),e(ldv),u(ldu,ldu),v(ldv,ldv),work(ldv)
!
!
!     dsvdc is a subroutine to reduce a double precision nxp matrix x
!     by orthogonal transformations u and v to diagonal form.  the
!     diagonal elements s(i) are the singular values of x.  the
!     columns of u are the corresponding left singular vectors,
!     and the columns of v the right singular vectors.
!
!     on entry
!
!         x         double precision(ldx,p), where ldx.ge.n.
!                   x contains the matrix whose singular value
!                   decomposition is to be computed.  x is
!                   destroyed by dsvdc.
!
!         ldx       integer.
!                   ldx is the leading dimension of the array x.
!
!         n         integer.
!                   n is the number of rows of the matrix x.
!
!         p         integer.
!                   p is the number of columns of the matrix x.
!
!         ldu       integer.
!                   ldu is the leading dimension of the array u.
!                   (see below).
!
!         ldv       integer.
!                   ldv is the leading dimension of the array v.
!                   (see below).
!
!         work      double precision(n).
!                   work is a scratch array.
!
!         job       integer.
!                   job controls the computation of the singular
!                   vectors.  it has the decimal expansion ab
!                   with the following meaning
!
!                        a.eq.0    do not compute the left singular
!                                  vectors.
!                        a.eq.1    return the n left singular vectors
!                                  in u.
!                        a.ge.2    return the first min(n,p) singular
!                                  vectors in u.
!                        b.eq.0    do not compute the right singular
!                                  vectors.
!                        b.eq.1    return the right singular vectors
!                                  in v.
!
!     on return
!
!         s         double precision(mm), where mm=min(n+1,p).
!                   the first min(n,p) entries of s contain the
!                   singular values of x arranged in descending
!                   order of magnitude.
!
!         e         double precision(p), 
!                   e ordinarily contains zeros.  however see the
!                   discussion of info for exceptions.
!
!         u         double precision(ldu,k), where ldu.ge.n.  if
!                                   joba.eq.1 then k.eq.n, if joba.ge.2
!                                   then k.eq.min(n,p).
!                   u contains the matrix of left singular vectors.
!                   u is not referenced if joba.eq.0.  if n.le.p
!                   or if joba.eq.2, then u may be identified with x
!                   in the subroutine call.
!
!         v         double precision(ldv,p), where ldv.ge.p.
!                   v contains the matrix of right singular vectors.
!                   v is not referenced if job.eq.0.  if p.le.n,
!                   then v may be identified with x in the
!                   subroutine call.
!
!         info      integer.
!                   the singular values (and their corresponding
!                   singular vectors) s(info+1),s(info+2),...,s(m)
!                   are correct (here m=min(n,p)).  thus if
!                   info.eq.0, all the singular values and their
!                   vectors are correct.  in any event, the matrix
!                   b = trans(u)*x*v is the bidiagonal matrix
!                   with the elements of s on its diagonal and the
!                   elements of e on its super-diagonal (trans(u)
!                   is the transpose of u).  thus the singular
!                   values of x and b are the same.
!
!     linpack. this version dated 08/14/78 .
!              correction made to shift 2/84.
!     g.w. stewart, university of maryland, argonne national lab.
!
!     dsvdc uses the following functions and subprograms.
!
!     external drot
!     blas daxpy,ddot,dscal,dswap,dnrm2,drotg
!     fortran dabs,dmax1,max0,min0,mod,dsqrt
!
!     internal variables
!
      integer i,iter,j,jobu,k,kase,kk,l,ll,lls,lm1,lp1,ls,lu,m,maxit, &
              mm,mm1,mp1,nct,nctp1,ncu,nrt,nrtp1
      _REAL_ ddot,dnrm2,t
      _REAL_ b,c,cs,el,emm1,f,g,scale,shift,sl,sm,sn, &
                       smm1,t1,test,ztest
      logical wantu,wantv
!
!     mjhsieh: warning eliminator
      l = -1; ls = -1
!
!     set the maximum number of iterations.
!
      maxit = 30
!
!     determine what is to be computed.
!
      wantu = .false.
      wantv = .false.
      jobu = mod(job,100)/10
      ncu = n
      if (jobu .gt. 1) ncu = min0(n,p)
      if (jobu .ne. 0) wantu = .true.
      if (mod(job,10) .ne. 0) wantv = .true.
!
!     reduce x to bidiagonal form, storing the diagonal elements
!     in s and the super-diagonal elements in e.
!
      info = 0
      nct = min0(n-1,p)
      nrt = max0(0,min0(p-2,n))
      lu = max0(nct,nrt)
      if (lu .lt. 1) go to 170
      do 160 l = 1, lu
         lp1 = l + 1
         if (l .gt. nct) go to 20
!
!           compute the transformation for the l-th column and
!           place the l-th diagonal in s(l).
!
            s(l) = dnrm2(n-l+1,x(l,l),1)
            if (s(l) .eq. 0.0d0) go to 10
               if (x(l,l) .ne. 0.0d0) s(l) = dsign(s(l),x(l,l))
               call dscal(n-l+1,1.0d0/s(l),x(l,l),1)
               x(l,l) = 1.0d0 + x(l,l)
   10       continue
            s(l) = -s(l)
   20    continue
         if (p .lt. lp1) go to 50
         do 40 j = lp1, p
            if (l .gt. nct) go to 30
            if (s(l) .eq. 0.0d0) go to 30
!
!              apply the transformation.
!
               t = -ddot(n-l+1,x(l,l),1,x(l,j),1)/x(l,l)
               call daxpy(n-l+1,t,x(l,l),1,x(l,j),1)
   30       continue
!
!           place the l-th row of x into  e for the
!           subsequent calculation of the row transformation.
!
            e(j) = x(l,j)
   40    continue
   50    continue
         if (.not.wantu .or. l .gt. nct) go to 70
!
!           place the transformation in u for subsequent back
!           multiplication.
!
            do 60 i = l, n
               u(i,l) = x(i,l)
   60       continue
   70    continue
         if (l .gt. nrt) go to 150
!
!           compute the l-th row transformation and place the
!           l-th super-diagonal in e(l).
!
            e(l) = dnrm2(p-l,e(lp1),1)
            if (e(l) .eq. 0.0d0) go to 80
               if (e(lp1) .ne. 0.0d0) e(l) = dsign(e(l),e(lp1))
               call dscal(p-l,1.0d0/e(l),e(lp1),1)
               e(lp1) = 1.0d0 + e(lp1)
   80       continue
            e(l) = -e(l)
            if (lp1 .gt. n .or. e(l) .eq. 0.0d0) go to 120
!
!              apply the transformation.
!
               do 90 i = lp1, n
                  work(i) = 0.0d0
   90          continue
               do 100 j = lp1, p
                  call daxpy(n-l,e(j),x(lp1,j),1,work(lp1),1)
  100          continue
               do 110 j = lp1, p
                  call daxpy(n-l,-e(j)/e(lp1),work(lp1),1,x(lp1,j),1)
  110          continue
  120       continue
            if (.not.wantv) go to 140
!
!              place the transformation in v for subsequent
!              back multiplication.
!
               do 130 i = lp1, p
                  v(i,l) = e(i)
  130          continue
  140       continue
  150    continue
  160 continue
  170 continue
!
!     set up the final bidiagonal matrix or order m.
!
      m = min0(p,n+1)
      nctp1 = nct + 1
      nrtp1 = nrt + 1
      if (nct .lt. p) s(nctp1) = x(nctp1,nctp1)
      if (n .lt. m) s(m) = 0.0d0
      if (nrtp1 .lt. m) e(nrtp1) = x(nrtp1,m)
      e(m) = 0.0d0
!
!     if required, generate u.
!
      if (.not.wantu) go to 300
         if (ncu .lt. nctp1) go to 200
         do 190 j = nctp1, ncu
            do 180 i = 1, n
               u(i,j) = 0.0d0
  180       continue
            u(j,j) = 1.0d0
  190    continue
  200    continue
         if (nct .lt. 1) go to 290
         do 280 ll = 1, nct
            l = nct - ll + 1
            if (s(l) .eq. 0.0d0) go to 250
               lp1 = l + 1
               if (ncu .lt. lp1) go to 220
               do 210 j = lp1, ncu
                  t = -ddot(n-l+1,u(l,l),1,u(l,j),1)/u(l,l)
                  call daxpy(n-l+1,t,u(l,l),1,u(l,j),1)
  210          continue
  220          continue
               call dscal(n-l+1,-1.0d0,u(l,l),1)
               u(l,l) = 1.0d0 + u(l,l)
               lm1 = l - 1
               if (lm1 .lt. 1) go to 240
               do 230 i = 1, lm1
                  u(i,l) = 0.0d0
  230          continue
  240          continue
            go to 270
  250       continue
               do 260 i = 1, n
                  u(i,l) = 0.0d0
  260          continue
               u(l,l) = 1.0d0
  270       continue
  280    continue
  290    continue
  300 continue
!
!     if it is required, generate v.
!
      if (.not.wantv) go to 350
         do 340 ll = 1, p
            l = p - ll + 1
            lp1 = l + 1
            if (l .gt. nrt) go to 320
            if (e(l) .eq. 0.0d0) go to 320
               do 310 j = lp1, p
                  t = -ddot(p-l,v(lp1,l),1,v(lp1,j),1)/v(lp1,l)
                  call daxpy(p-l,t,v(lp1,l),1,v(lp1,j),1)
  310          continue
  320       continue
            do 330 i = 1, p
               v(i,l) = 0.0d0
  330       continue
            v(l,l) = 1.0d0
  340    continue
  350 continue
!
!     main iteration loop for the singular values.
!
      mm = m
      iter = 0
  360 continue
!
!        quit if all the singular values have been found.
!
!     ...exit
         if (m .eq. 0) go to 620
!
!        if too many iterations have been performed, set
!        flag and return.
!
         if (iter .lt. maxit) go to 370
            info = m
!     ......exit
            go to 620
  370    continue
!
!        this section of the program inspects for
!        negligible elements in the s and e arrays.  on
!        completion the variables kase and l are set as follows.
!
!           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m
!           kase = 2     if s(l) is negligible and l.lt.m
!           kase = 3     if e(l-1) is negligible, l.lt.m, and
!                        s(l), ..., s(m) are not negligible (qr step).
!           kase = 4     if e(m-1) is negligible (convergence).
!
         do 390 ll = 1, m
            l = m - ll
!        ...exit
            if (l .eq. 0) go to 400
            test = dabs(s(l)) + dabs(s(l+1))
            ztest = test + dabs(e(l))
            if (ztest .ne. test) go to 380
               e(l) = 0.0d0
!        ......exit
               go to 400
  380       continue
  390    continue
  400    continue
         if (l .ne. m - 1) go to 410
            kase = 4
         go to 480
  410    continue
            lp1 = l + 1
            mp1 = m + 1
            do 430 lls = lp1, mp1
               ls = m - lls + lp1
!           ...exit
               if (ls .eq. l) go to 440
               test = 0.0d0
               if (ls .ne. m) test = test + dabs(e(ls))
               if (ls .ne. l + 1) test = test + dabs(e(ls-1))
               ztest = test + dabs(s(ls))
               if (ztest .ne. test) go to 420
                  s(ls) = 0.0d0
!           ......exit
                  go to 440
  420          continue
  430       continue
  440       continue
            if (ls .ne. l) go to 450
               kase = 3
            go to 470
  450       continue
            if (ls .ne. m) go to 460
               kase = 1
            go to 470
  460       continue
               kase = 2
               l = ls
  470       continue
  480    continue
         l = l + 1
!
!        perform the task indicated by kase.
!
         go to (490,520,540,570), kase
!
!        deflate negligible s(m).
!
  490    continue
            mm1 = m - 1
            f = e(m-1)
            e(m-1) = 0.0d0
            do 510 kk = l, mm1
               k = mm1 - kk + l
               t1 = s(k)
               call drotg(t1,f,cs,sn)
               s(k) = t1
               if (k .eq. l) go to 500
                  f = -sn*e(k-1)
                  e(k-1) = cs*e(k-1)
  500          continue
               if (wantv) call drot(p,v(1,k),1,v(1,m),1,cs,sn)
  510       continue
         go to 610
!
!        split at negligible s(l).
!
  520    continue
            f = e(l-1)
            e(l-1) = 0.0d0
            do 530 k = l, m
               t1 = s(k)
               call drotg(t1,f,cs,sn)
               s(k) = t1
               f = -sn*e(k)
               e(k) = cs*e(k)
               if (wantu) call drot(n,u(1,k),1,u(1,l-1),1,cs,sn)
  530       continue
         go to 610
!
!        perform one qr step.
!
  540    continue
!
!           calculate the shift.
!
            scale = dmax1(dabs(s(m)),dabs(s(m-1)),dabs(e(m-1)), &
                          dabs(s(l)),dabs(e(l)))
            sm = s(m)/scale
            smm1 = s(m-1)/scale
            emm1 = e(m-1)/scale
            sl = s(l)/scale
            el = e(l)/scale
            b = ((smm1 + sm)*(smm1 - sm) + emm1**2)/2.0d0
            c = (sm*emm1)**2
            shift = 0.0d0
            if (b .eq. 0.0d0 .and. c .eq. 0.0d0) go to 550
               shift = dsqrt(b**2+c)
               if (b .lt. 0.0d0) shift = -shift
               shift = c/(b + shift)
  550       continue
            f = (sl + sm)*(sl - sm) + shift
            g = sl*el
!
!           chase zeros.
!
            mm1 = m - 1
            do 560 k = l, mm1
               call drotg(f,g,cs,sn)
               if (k .ne. l) e(k-1) = f
               f = cs*s(k) + sn*e(k)
               e(k) = cs*e(k) - sn*s(k)
               g = sn*s(k+1)
               s(k+1) = cs*s(k+1)
               if (wantv) call drot(p,v(1,k),1,v(1,k+1),1,cs,sn)
               call drotg(f,g,cs,sn)
               s(k) = f
               f = cs*e(k) + sn*s(k+1)
               s(k+1) = -sn*e(k) + cs*s(k+1)
               g = sn*e(k+1)
               e(k+1) = cs*e(k+1)
               if (wantu .and. k .lt. n) call drot(n,u(1,k),1,u(1,k+1),1,cs,sn)
  560       continue
            e(m-1) = f
            iter = iter + 1
         go to 610
!
!        convergence.
!
  570    continue
!
!           make the singular value  positive.
!
            if (s(l) .ge. 0.0d0) go to 580
               s(l) = -s(l)
               if (wantv) call dscal(p,-1.0d0,v(1,l),1)
  580       continue
!
!           order the singular value.
!
  590       if (l .eq. mm) go to 600
!           ...exit
               if (s(l) .ge. s(l+1)) go to 600
               t = s(l)
               s(l) = s(l+1)
               s(l+1) = t
               if (wantv .and. l .lt. p) call dswap(p,v(1,l),1,v(1,l+1),1)
               if (wantu .and. l .lt. n) call dswap(n,u(1,l),1,u(1,l+1),1)
               l = l + 1
            go to 590
  600       continue
            iter = 0
            m = m - 1
  610    continue
      go to 360
  620 continue
      return
      end
