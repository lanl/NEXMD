c------------------------------------------------------------------
      subroutine direct(numwats,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial)
      implicit none
      double precision x(*),y(*),z(*),cg(*),cutoff,ewaldcof,box,
     $      ene,fx(*),fy(*),fz(*),virial(6)
      integer numwats

      integer i,n,j,ilo,jlo,numatoms,minimg_flag

      numatoms = 3*numwats
      ene = 0.d0
      do 10 i = 1,6
        virial(i) = 0.d0
10    continue
      
c first handle atoms with themselves. Get all but minimum image
      minimg_flag = 0
      do 100 n = 1,numwats
       ilo = 3*(n-1)+1
       call get_dir_pair(ilo,ilo,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial,minimg_flag)
       call get_dir_pair(ilo+1,ilo+1,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial,minimg_flag)
       call get_dir_pair(ilo+2,ilo+2,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial,minimg_flag)
100   continue
      ene = 0.5d0*ene
      do 150 i = 1,6
        virial(i) = 0.5d0*virial(i)
150   continue
c Next handle pairs on same water. Get all but minimum image
      do 200 n = 1,numwats
       ilo = 3*(n-1)+1
       call get_dir_pair(ilo,ilo+1,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial,minimg_flag)
       call get_dir_pair(ilo,ilo+2,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial,minimg_flag)
       call get_dir_pair(ilo+1,ilo+2,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial,minimg_flag)
200   continue
       
c Next handle all other pairs. Get all images.
      minimg_flag = 1
      do 300 n = 1,numwats 
       jlo = 3*n+1
       ilo = 3*(n-1)+1
       do 250 j = jlo,numatoms
         call get_dir_pair(ilo,j,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial,minimg_flag)
         call get_dir_pair(ilo+1,j,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial,minimg_flag)
         call get_dir_pair(ilo+2,j,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial,minimg_flag)
250    continue
300   continue
      return
      end
c------------------------------------------------------------------
      subroutine get_dir_pair(i,j,x,y,z,cg,cutoff,ewaldcof,box,
     $      ene,fx,fy,fz,virial,minimg_flag)
      implicit none
      integer i,j,minimg_flag
      double precision x(*),y(*),z(*),cg(*),cutoff,ewaldcof,box,
     $      ene,fx(*),fy(*),fz(*),virial(6)

      double precision xij,yij,zij,r,crd(3),grad(3),pot
      double precision xij0,yij0,zij0,boxinv
      integer k,l,m

      boxinv = 1.d0 / box
      xij0 = x(i) - x(j)
      yij0 = y(i) - y(j)
      zij0 = z(i) - z(j)
      xij0 = xij0 - box*anint(xij0*boxinv)
      yij0 = yij0 - box*anint(yij0*boxinv)
      zij0 = zij0 - box*anint(zij0*boxinv)
      do 50 m = 1,3
       zij = zij0 + (m-2)*box
       do 45 l = 1,3
        yij = yij0 + (l-2)*box
        do 40 k = 1,3
c check for doing the minimium image
         if ( k.eq.2 .and. l.eq.2 .and. m.eq.2 
     $     .and. minimg_flag .eq. 0)goto 39
         xij = xij0 + (k-2)*box
         r = sqrt(xij**2+yij**2+zij**2) 
         if ( r .le. cutoff )then
          crd(1) = xij
          crd(2) = yij
          crd(3) = zij
          call ew_direct(ewaldcof,crd,pot,grad)
          ene = ene + cg(i)*cg(j)*pot
          fx(i) = fx(i) - cg(i)*cg(j)*grad(1)
          fx(j) = fx(j) + cg(i)*cg(j)*grad(1)
          fy(i) = fy(i) - cg(i)*cg(j)*grad(2)
          fy(j) = fy(j) + cg(i)*cg(j)*grad(2)
          fz(i) = fz(i) - cg(i)*cg(j)*grad(3)
          fz(j) = fz(j) + cg(i)*cg(j)*grad(3)
          virial(1) = virial(1) + cg(i)*cg(j)*xij*grad(1)
          virial(2) = virial(2) + cg(i)*cg(j)*xij*grad(2)
          virial(3) = virial(3) + cg(i)*cg(j)*xij*grad(3)
          virial(4) = virial(4) + cg(i)*cg(j)*yij*grad(2)
          virial(5) = virial(5) + cg(i)*cg(j)*yij*grad(3)
          virial(6) = virial(6) + cg(i)*cg(j)*zij*grad(3)
         endif
39       continue
40      continue
45     continue
50    continue
      return
      end
c---------------------------------------------------------
      subroutine ew_direct(ewaldcof,crd,pot,grad)
      implicit none
      double precision ewaldcof,crd(3),pot,grad(3)

      double precision pi,r,x,term,merfc
      integer i

      pi = 3.14159265358979323846

c assume x,y,z are minimum image
      
      r = sqrt(crd(1)**2+crd(2)**2+crd(3)**2)
      x = r*ewaldcof
      call erfcfun(x,merfc)
      pot = merfc/r
      term = (2.d0/sqrt(pi))*ewaldcof*exp(-x**2)/r +
     $       merfc/r**2
      do 20 i = 1,3
        grad(i) = -term*crd(i)/r
20    continue
      return
      end
c------------------------------------------------------------------------
      subroutine  adjust(numwats,x,y,z,cg,ewaldcof,ene,fx,fy,fz,vir)
      implicit none
      double precision x(*),y(*),z(*),cg(*),ewaldcof,
     $      ene,fx(*),fy(*),fz(*),vir(6)
      integer numwats

      integer i,n,ilo,numatoms

      numatoms = 3*numwats
      ene = 0.d0
      do 10 i = 1,6
        vir(i) = 0.d0
10    continue

      do 100 n = 1,numwats
       ilo = 3*(n-1)+1
       call get_adj_pair(ilo,ilo+1,x,y,z,cg,ewaldcof,ene,fx,fy,fz,vir)
       call get_adj_pair(ilo,ilo+2,x,y,z,cg,ewaldcof,ene,fx,fy,fz,vir)
       call get_adj_pair(ilo+1,ilo+2,x,y,z,cg,ewaldcof,ene,fx,fy,fz,vir)
100   continue
      return
      end
c------------------------------------------------------------------
      subroutine get_adj_pair(i,j,x,y,z,cg,ewaldcof,ene,fx,fy,fz,vir)
      implicit none
      integer i,j
      double precision x(*),y(*),z(*),cg(*),ewaldcof,
     $      ene,fx(*),fy(*),fz(*),vir(6)

      double precision crd(3),grad(3),pot
      crd(1) = x(i)-x(j)
      crd(2) = y(i)-y(j)
      crd(3) = z(i)-z(j)
      call ew_adjust(ewaldcof,crd,pot,grad)
      ene = ene + cg(i)*cg(j)*pot
      fx(i) = fx(i) - cg(i)*cg(j)*grad(1)
      fx(j) = fx(j) + cg(i)*cg(j)*grad(1)
      fy(i) = fy(i) - cg(i)*cg(j)*grad(2)
      fy(j) = fy(j) + cg(i)*cg(j)*grad(2)
      fz(i) = fz(i) - cg(i)*cg(j)*grad(3)
      fz(j) = fz(j) + cg(i)*cg(j)*grad(3)
      vir(1) = vir(1) + cg(i)*cg(j)*crd(1)*grad(1)
      vir(2) = vir(2) + cg(i)*cg(j)*crd(1)*grad(2)
      vir(3) = vir(3) + cg(i)*cg(j)*crd(1)*grad(3)
      vir(4) = vir(4) + cg(i)*cg(j)*crd(2)*grad(2)
      vir(5) = vir(5) + cg(i)*cg(j)*crd(2)*grad(3)
      vir(6) = vir(6) + cg(i)*cg(j)*crd(3)*grad(3)
      return
      end
c------------------------------------------------------------------
      subroutine ew_adjust(ewaldcof,crd,pot,grad)
      implicit none
      double precision ewaldcof,crd(3),pot,grad(3)

      double precision pi,r,x,term,merfc
      integer i

      pi = 3.14159265358979323846

c assume x,y,z are minimum image
      
      r = sqrt(crd(1)**2+crd(2)**2+crd(3)**2)
      x = r*ewaldcof
      call erfcfun(x,merfc)
      pot = -(1.d0-merfc)/r
      term = -(2.d0/sqrt(pi))*ewaldcof*exp(-x**2)/r +
     $       (1.d0-merfc)/r**2
      do 20 i = 1,3
        grad(i) = term*crd(i)/r
20    continue
      return
      end
c------------------------------------------------------------------------
      subroutine self(cg,numatoms,ene,ewaldcof)
      implicit none
      integer numatoms
      double precision cg(*),ene,ewaldcof

      integer i
      double precision ee,pi
      ee = 0.d0
      pi = 3.14159265358979323846
      do 10 i = 1,numatoms
        ee = ee + cg(i)**2
10    continue
      ene = -ee*ewaldcof/sqrt(pi)
      return
      end
c------------------------------------------------------------------
      subroutine find_ewaldcof(cutoff,dtol,ewaldcof)
      implicit none
      double precision cutoff,dtol,ewaldcof

      integer i,n
      double precision pi,term,x,xlo,xhi,y,erfc

c first get direct sum tolerance. How big must ewaldcof be to get
c terms outside the cutoff below tol
      pi = 3.14159265358979323846

      x = 0.5d0
      i = 0
10    x = 2.d0 * x
      i = i + 1
      y = x*cutoff
      call erfcfun(y,erfc)
      term = erfc/cutoff
      if ( term .ge. dtol)goto 10
c binary search tolerance is 2 to the -60th
      n = i + 60
      xlo = 0.d0
      xhi = x
      do 20 i = 1,n
        x = (xlo+xhi)/2
        y = x*cutoff
        call erfcfun(y,erfc)
        term = erfc/cutoff
        if ( term .ge. dtol )then
           xlo = x
        else 
           xhi = x
        endif
20    continue
      ewaldcof = x

      return
      end
c----------------------------------------------
      subroutine find_cutoff(cutoff,dtol,ewaldcof)
      implicit none
      double precision cutoff,dtol,ewaldcof

      integer i,n
      double precision pi,term,x,xlo,xhi,y,erfc

c first get direct sum tolerance. How big must ewaldcof be to get
c terms outside the cutoff below tol
      pi = 3.14159265358979323846

      x = 0.5d0
      i = 0
10    x = 2.d0 * x
      i = i + 1
      y = x*ewaldcof
      call erfcfun(y,erfc)
      term = erfc/x
      if ( term .ge. dtol)goto 10
c binary search tolerance is 2 to the -60th
      n = i + 60
      xlo = 0.d0
      xhi = x
      do 20 i = 1,n
        x = (xlo+xhi)/2
        y = x*ewaldcof
        call erfcfun(y,erfc)
        term = erfc/x
        if ( term .ge. dtol )then
           xlo = x
        else 
           xhi = x
        endif
20    continue
      cutoff = x

      return
      end
c------------------------------------------------------
      subroutine find_maxexp(ewaldcof,rtol,maxexp)
      implicit none
      double precision ewaldcof,rtol,maxexp

      integer i,n
      double precision pi,term,x,xlo,xhi,y,erfc

      pi = 3.14159265358979323846
      x = 0.5d0
      i = 0
30    x = 2.d0 * x
      i = i + 1
      y = pi*x/ewaldcof
      call erfcfun(y,erfc)
      term = 2.d0*ewaldcof*erfc/sqrt(pi)
      if ( term .ge. rtol)goto 30
c binary search tolerance is 2 to the -60th
      n = i + 60
      xlo = 0.d0
      xhi = x
      do 40 i = 1,n
        x = (xlo+xhi)/2
        y = pi*x/ewaldcof
        call erfcfun(y,erfc)
        term = 2.d0*ewaldcof*erfc/sqrt(pi)
        if ( term .gt. rtol )then
           xlo = x
        else 
           xhi = x
        endif
40    continue
      maxexp = x

      return
      end
c------------------------------------------------------
      subroutine get_mlim(maxexp,mlimit,eigmin,reclng,recip)
      implicit none
      double precision maxexp
      double precision eigmin
      integer mlimit(3)
      double precision reclng(3),recip(3,3)

c get coefficients for reciprocal space ewald sum

      integer mtop1,mtop2,mtop3,mlim1,mlim2,mlim3
      integer m1,m2,m3,nrecvecs
      double precision z1,z2,z3,expo
      double precision pi

      pi = 3.14159265358979323846
      mtop1 = reclng(1)*maxexp/sqrt(eigmin)
      mtop2 = reclng(2)*maxexp/sqrt(eigmin)
      mtop3 = reclng(3)*maxexp/sqrt(eigmin)

      nrecvecs = 0
      mlim1 = 0
      mlim2 = 0
      mlim3 = 0
      do 100 m1 = -mtop1,mtop1
      do 100 m2 = -mtop2,mtop2
      do 100 m3 = -mtop3,mtop3
        z1 = m1*recip(1,1)+m2*recip(1,2)+m3*recip(1,3) 
        z2 = m1*recip(2,1)+m2*recip(2,2)+m3*recip(2,3) 
        z3 = m1*recip(3,1)+m2*recip(3,2)+m3*recip(3,3) 
        expo = z1**2 + z2**2 + z3**2
        if ( expo .le. maxexp**2 )then
          nrecvecs = nrecvecs + 1
          if ( abs(m1) .gt. mlim1 )mlim1 = abs(m1) 
          if ( abs(m2) .gt. mlim2 )mlim2 = abs(m2) 
          if ( abs(m3) .gt. mlim3 )mlim3 = abs(m3) 
        endif
100   continue
      write(6,*)'number of reciprocal vecs = ',nrecvecs
      write(6,*)'mlim1,2,3 = ',mlim1,mlim2,mlim3
      mlimit(1) = mlim1
      mlimit(2) = mlim2
      mlimit(3) = mlim3
      return
      end
c-------------------------------------------------------------
c------------------------------------------------------
      SUBROUTINE SECOND(TIM)
      implicit none
      real tim

c     1 f77 version
c     use etime() for elapsed time from start of execution
c     use dtime() for delta time from last call to dtime()
      real dumm(2),etime
      tim = etime(dumm)
      return
      end
c-------------------------------------------------------------
      subroutine check_virial(self_ene,adj_ene,dir_ene,rec_ene,
     $       adj_vir,rec_vir,dir_vir)
      implicit none
      double precision self_ene,adj_ene,dir_ene,rec_ene
      double precision adj_vir(6),rec_vir(6),dir_vir(6)

      double precision etot,svir,relerr
      etot = self_ene+adj_ene+dir_ene+rec_ene
      svir = adj_vir(1)+rec_vir(1)+dir_vir(1)+
     $       adj_vir(4)+rec_vir(4)+dir_vir(4)+
     $       adj_vir(6)+rec_vir(6)+dir_vir(6)
      relerr = 2.d0*abs(etot+svir)/(abs(etot)+abs(svir))
      write(6,*)'tot ene =   ',etot
      write(6,*)'trace vir = ',svir
      write(6,*)'rel error = ',relerr
      return
      end
c-------------------------------------------------------------
      subroutine check_force(numatoms,
     $      fx1,fy1,fz1,fx2,fy2,fz2,fdx,fdy,fdz)
      integer numatoms
      double precision fx1(*),fy1(*),fz1(*),fx2(*),fy2(*),fz2(*),
     $      fdx(*),fdy(*),fdz(*)

      double precision rms_num,rms_den,rms
      integer i
      rms_num = 0.d0
      rms_den = 0.d0
      do 100 i = 1,numatoms
       rms_num = rms_num + (fx1(i)-fx2(i))**2 + (fy1(i)-fy2(i))**2 +
     $          (fz1(i)-fz2(i))**2
       rms_den = rms_den + fdx(i)**2 + fdy(i)**2 + fdz(i)**2
100   continue
      rms = sqrt(rms_num/rms_den)
      write(6,*)'rms force err = ',rms
      return
      end
c-------------------------------------------------------------
