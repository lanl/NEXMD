
!       ---LCPO code: to be included in egb.f

call timer_start(TIME_GBSA)
max_count = count
count=1
#ifdef MPI
do i=mytaskid+1,natom,numtasks
#else
do i=1,natom
#endif
   if (ineighbor(count) == 0) then
      count=count+1
   else
      si = FOURPI*(vdwrad(i)*vdwrad(i))
      if(idecomp == 1 .or. idecomp == 2) then
         call decsasa(1,i,0,0,surften*p1(i)*si)
      else if(idecomp == 3 .or. idecomp == 4) then
         call decsasa(-1,i,0,0,surften*p1(i)*si)
      end if
      sumaij            = 0.d0
      sumajk            = 0.d0
      sumaijajk         = 0.d0
      sumdaijddijdxi    = 0.d0
      sumdaijddijdyi    = 0.d0
      sumdaijddijdzi    = 0.d0
      sumdaijddijdxiajk = 0.d0
      sumdaijddijdyiajk = 0.d0
      sumdaijddijdziajk = 0.d0
      
      !       --- loop over j ---
      
      icount=count
      num_j_vals         = 0
      70 j                  = ineighbor(count)
      num_j_vals         = num_j_vals + 1
      j_vals(num_j_vals) = j
      count              = count + 1
      if (ineighbor(count) /= 0) then
         goto 70
      else
         count            = count + 1
      end if
      xi = x(3*i-2)
      yi = x(3*i-1)
      zi = x(3*i  )
      do jjj = 1, num_j_vals
         j = j_vals(jjj)
         xj = x(3*j-2)
         yj = x(3*j-1)
         zj = x(3*j  )
         r2 = (xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
         dij1i = 1.d0/sqrt(r2)
         rij = r2*dij1i
         tmpaij = vdwrad(i) - rij*0.5d0 &
               - (vdwrad(i)**2 - vdwrad(j)**2)*0.5d0*dij1i
         aij = 2.d0*pi*vdwrad(i)*tmpaij
         
         !       --- first derivatives ---
         
         daijddij = pi*vdwrad(i)* &
               (dij1i*dij1i*(vdwrad(i)**2 - vdwrad(j)**2) - 1)
         daijddijdxj=daijddij*(xj-xi)*dij1i
         daijddijdyj=daijddij*(yj-yi)*dij1i
         daijddijdzj=daijddij*(zj-zi)*dij1i
         
         !       ----------------------------
         
         sumaij = sumaij+aij
         if(idecomp == 1 .or. idecomp == 2) then
            call decsasa(2,i,j,0,surften*p2(i)*aij)
         else if(idecomp == 3 .or. idecomp == 4) then
            call decsasa(-2,i,j,0,surften*p2(i)*aij)
         end if
         
         count2=icount
         sumajk2=0.d0
         sumdajkddjkdxj=0.d0
         sumdajkddjkdyj=0.d0
         sumdajkddjkdzj=0.d0
         p3p4aij = -surften*(p3(i) + p4(i)*aij)*frespa
         
         !       --- loop over k ---
         
         do kk=count2,max_count
            if(ineighbor(kk) == 0)then
               count2_fin=kk-1
               goto 81
            end if
         end do
         81 continue
         num_k_vals = 0
         do kk=count2,count2_fin
            k=ineighbor(kk)
            if(k /= j)then
               num_k_vals= num_k_vals+1
               k_vals(num_k_vals) = k
            end if
         end do
         count2=icount
         
!!!         !dir$ ivdep
         do kk=1,num_k_vals
            k=k_vals(kk)
            xk = x(3*k-2)
            yk = x(3*k-1)
            zk = x(3*k  )
            
            rjk2 = (xj-xk)**2 + (yj-yk)**2 + (zj-zk)**2
            djk1i = 1.d0/sqrt(rjk2)
            rjk = rjk2*djk1i
            if ((vdwrad(j) + vdwrad(k)) > rjk) then
               vdw2dif = vdwrad(j)**2 - vdwrad(k)**2
               tmpajk = 2.d0*vdwrad(j) - rjk - vdw2dif*djk1i
               
               ajk = pi*vdwrad(j)*tmpajk
               sumajk = sumajk+ajk
               sumajk2 = sumajk2+ajk
               if(idecomp == 1 .or. idecomp == 2) then
                  call decsasa(3,i,j,k, &
                        surften*(p3(i)*ajk + p4(i)*aij*ajk))
               else if(idecomp == 3 .or. idecomp == 4) then
                  call decsasa(-3,i,j,k, &
                        surften*(p3(i)*ajk + p4(i)*aij*ajk))
               end if
               
               !             --- first derivatives ---
               
               dajkddjk=pi*vdwrad(j)*djk1i* &
                     (djk1i*djk1i*vdw2dif - 1.d0)
               
               dajkddjkdxj=dajkddjk*(xj-xk)
               dajkddjkdyj=dajkddjk*(yj-yk)
               dajkddjkdzj=dajkddjk*(zj-zk)
               
               f(3*k-2) = f(3*k-2) - dajkddjkdxj*p3p4aij
               f(3*k-1) = f(3*k-1) - dajkddjkdyj*p3p4aij
               f(3*k  ) = f(3*k  ) - dajkddjkdzj*p3p4aij
               
               sumdajkddjkdxj=sumdajkddjkdxj+dajkddjkdxj
               sumdajkddjkdyj=sumdajkddjkdyj+dajkddjkdyj
               sumdajkddjkdzj=sumdajkddjkdzj+dajkddjkdzj
            end if  ! ((vdwrad(j) + vdwrad(k)) > rjk)
         end do  !  kk=1,num_k_vals
         
         
         !         --- finished looping over k ---
         
         sumaijajk = sumaijajk+aij*sumajk2
         
         !         --- first derivatives ---
         
         sumdaijddijdxi = sumdaijddijdxi - daijddijdxj
         sumdaijddijdyi = sumdaijddijdyi - daijddijdyj
         sumdaijddijdzi = sumdaijddijdzi - daijddijdzj
         sumdaijddijdxiajk = sumdaijddijdxiajk - &
               daijddijdxj*sumajk2
         sumdaijddijdyiajk = sumdaijddijdyiajk - &
               daijddijdyj*sumajk2
         sumdaijddijdziajk = sumdaijddijdziajk - &
               daijddijdzj*sumajk2
         
         lastxj=daijddijdxj*sumajk2+aij*sumdajkddjkdxj
         lastyj=daijddijdyj*sumajk2+aij*sumdajkddjkdyj
         lastzj=daijddijdzj*sumajk2+aij*sumdajkddjkdzj
         
         daidxj=surften*(p2(i)*daijddijdxj &
               +p3(i)*sumdajkddjkdxj+p4(i)*lastxj)
         daidyj=surften*(p2(i)*daijddijdyj &
               +p3(i)*sumdajkddjkdyj+p4(i)*lastyj)
         daidzj=surften*(p2(i)*daijddijdzj &
               +p3(i)*sumdajkddjkdzj+p4(i)*lastzj)
         
         f(3*j-2) = f(3*j-2) - daidxj*frespa
         f(3*j-1) = f(3*j-1) - daidyj*frespa
         f(3*j  ) = f(3*j  ) - daidzj*frespa
         
         !             write(6,*) 'Aij,sumAjk,sumAijAjk',Aij,sumAjk,sumAijAjk
         
      end do  !  jjj = 1, num_j_vals
      
      !         --- finished looping over j ---
      
      !           write(6,'(3f12.5)') sumAij,sumAjk,sumAijAjk
      
      ai = p1(i)*si+p2(i)*sumaij+p3(i)*sumajk &
            +p4(i)*sumaijajk
      
      daidxi = surften*(p2(i)*sumdaijddijdxi &
            +p4(i)*sumdaijddijdxiajk)
      daidyi = surften*(p2(i)*sumdaijddijdyi &
            +p4(i)*sumdaijddijdyiajk)
      daidzi = surften*(p2(i)*sumdaijddijdzi &
            +p4(i)*sumdaijddijdziajk)
      
      f(3*i-2) = f(3*i-2) - daidxi*frespa
      f(3*i-1) = f(3*i-1) - daidyi*frespa
      f(3*i  ) = f(3*i  ) - daidzi*frespa
      
      totsasa = totsasa + ai
      
   end if  ! (ineighbor(count) == 0)
end do  !  i=mytaskid+1,natom,numtasks

!       --- finished looping over i ---

call timer_stop(TIME_GBSA)
