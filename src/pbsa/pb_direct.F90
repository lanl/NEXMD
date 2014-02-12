
! <compile=optimized>
#include "copyright.h"
#define _REAL_ double precision
#include "pb_def.h"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Driver of direct coulombic and vdw energy and force
subroutine pb_directnocut( natom,proatm,inatm,ipres,ibgwat,ienwat,ibgion,ienion,ntypes,eneopt, &
                           idecomp,ifcap,iac,ico,nex,iex,cn1,cn2,cg,x,f,eel,enb )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   implicit none

#  include "pb_constants.h"
#  include "flocntrl.h"
    
   ! Passed variables
    
   integer natom, proatm, inatm, ibgwat, ienwat, ibgion, ienion, ntypes, eneopt, idecomp, ifcap
   integer ipres(*), iac(*), ico(*), nex(natom), iex(64,natom)
   _REAL_ cn1(*), cn2(*), cg(natom)
   _REAL_ x(3,natom)
   _REAL_ enb, eel, f(3,natom)
    
   ! Local variables
    
   integer iaci, ic, watres, watfirst, ionfirst
   _REAL_ cn1oo, cn2oo, cgoo, cgoh, cghh
   _REAL_ eelwat, enbwat, eelpro, enbpro, eelprowat, enbprowat
    
   if ( do_pbdir == 0 ) return

   eelwat = ZERO; enbwat = ZERO
   eelpro = ZERO; enbpro = ZERO
   eelprowat = ZERO; enbprowat = ZERO

   if ( ibgwat /= 0 ) then
      proatm = ipres(ibgwat) - 1
   else
      proatm = natom
   end if

   if((ifcap == 2 .or. ifcap == 5) .and. proatm > inatm) then
      write(6, *) 'PB bomb in pb_directnocut(): proatm must not be larger than inatm here'
      call mexit(6, 1)
   end if

   watfirst = 0
   watres = 0
   ionfirst = 0
   if ( ibgwat /= 0 ) then
      watfirst = proatm + 1
      if(ifcap == 2 .or. ifcap == 5) then
         watres = (inatm - watfirst + 1) / 3
      else
         watres = ienwat - ibgwat + 1
      end if
      if( ibgion /= 0) then
        ionfirst = ipres(ibgion)
      endif
   end if
    
   ! get protein/protein and protein/water interactions
    
   call pb_dirpro(proatm,natom,watres,watfirst,ionfirst,ntypes,ifcap,nex,iex,iac,ico,cn1,cn2,cg,x,f, &
                  eelpro,enbpro,eelprowat,enbprowat)

   if(ifcap == 2 .or. ifcap == 5) then
      write(6,'(a,2f10.3)') ' Protein-solvent interactions: ', eelprowat, enbprowat
   end if
    
   ! get water/water interactions
   ! first get van der Waals and charge-charge pairs
    
   if ( ibgwat /= 0 ) then
      iaci = ntypes*(iac(watfirst)-1)
      ic = ico(iaci+iac(watfirst))
      cn1oo = cn1(ic); cn2oo = cn2(ic)
      cgoo = cg(watfirst)**2
      cgoh = cg(watfirst)*cg(watfirst+1)
      cghh = cg(watfirst+1)**2
      call pb_dirwat(watfirst,watres,ibgwat,ienwat,cn1oo,cn2oo,cgoo,cgoh,cghh,x,f,eelwat,enbwat)
   end if
    
   eel = eelwat + eelpro + eelprowat
   enb = enbwat + enbpro + enbprowat
    
contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ compute direct non-water/all nonbonded interactions
subroutine pb_dirpro( proatm,natom,watres,watfirst,ionfirst,ntypes,ifcap,nex,iex,iac,ico,cn1,cn2,cg,x,f, &
                      eelpro,enbpro,eelprowat,enbprowat )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   use poisson_boltzmann, only: outflag
   use decomp, only: decpair
    
   implicit none

#  include "pb_constants.h"
    
   ! Passed variables
    
   integer proatm, natom, watres, watfirst, ionfirst, ntypes, ifcap
   integer iac(*), ico(*), nex(natom), iex(64, natom)
   _REAL_ cn1(*), cn2(*), cg(natom)
   _REAL_ x(3,natom)
   _REAL_ eelpro, enbpro, eelprowat, enbprowat, f(3,natom)
    
   ! Local variables
    
   integer jex,i,j,ilast,jfirst,jlast,iaci,ic,jp,iatm
   _REAL_ cgi,cgj,cn1ij,cn2ij,cgio,cgih,cn1io,cn2io
   _REAL_ dumx,dumy,dumz,xi,yi,zi,dx,dy,dz,d2inv,r6
   _REAL_ f1,f2,df,df2,fw1,fw2,fw3
   _REAL_ eelprotmp, enbprotmp, eelprowattmp, enbprowattmp
   _REAL_ xwij1(3*watres), xwij2(3*watres), xwij3(3*watres), rw1(3*watres)
   _REAL_ fx(3*watres), fy(3*watres), fz(3*watres)
   _REAL_ coulx(natom),couly(natom),coulz(natom)
   _REAL_ vdwx(natom),vdwy(natom),vdwz(natom)
   _REAL_ ofx(natom),ofy(natom),ofz(natom)

    
   ! Compute interaction energy and forces for non-water atoms 
    
   vdwx = ZERO; vdwy=ZERO; vdwz=ZERO
   coulx = ZERO; couly=ZERO; coulz=ZERO
   ofx = f(1,1:natom); ofy = f(2,1:natom); ofz = f(3,1:natom)

   ilast = proatm
   do i = 1, ilast
      if((ifcap == 2 .or. ifcap == 5) .and. outflag(i) == 1) cycle
      iaci = ntypes*(iac(i)-1); cgi = cg(i)
      xi = x(1,i); yi = x(2,i); zi = x(3,i)
      dumx = ZERO; dumy = ZERO; dumz = ZERO
      eelprotmp = ZERO; enbprotmp = ZERO
      eelprowattmp = ZERO; enbprowattmp = ZERO
       
      ! pro/pro interactions
       
      do j = i+1, ilast
         if((ifcap == 2 .or. ifcap == 5) .and. outflag(j) == 1) cycle
         do jp = 1, nex(i)
            jex = iex(jp,i)
            if (j == jex) goto 10
         end do
         dx = xi - x(1,j); dy = yi - x(2,j); dz = zi - x(3,j); cgj=cg(j)
         d2inv = ONE/(dx**2+dy**2+dz**2)
         df2 = -cgi*cgj*sqrt(d2inv)
         if(i < ionfirst .and. j >= ionfirst) then
            ! Count protein <-> ion interactions separately
            eelprowattmp = eelprowattmp + df2
         else
            eelprotmp = eelprotmp + df2
         end if
         ic = ico(iaci+iac(j)); cn1ij = cn1(ic); cn2ij = cn2(ic)
         r6 = d2inv**3; f2 = cn2ij*r6; f1 = cn1ij*(r6*r6)
         if(i < ionfirst .and. j >= ionfirst) then
            ! Count protein <-> ion interactions separately
            enbprowattmp = enbprowattmp + (f2-f1)
         else
            enbprotmp = enbprotmp + (f2-f1)
         end if

         !-- PB decomp
         if(idecomp == 1 .or. idecomp == 2) then
            if(eneopt == 1) then
               call decpair(1,i,j,df2)
            end if
            call decpair(2,i,j,-df2)
            call decpair(3,i,j,(f1 - f2))
         else if(idecomp == 3 .or. idecomp == 4) then
            call decpair(-2,i,j,-df2)
            call decpair(-3,i,j,(f1 - f2))
         end if

         df = ( df2 + SIX*( (f2-f1)-f1 ) )*d2inv
         fw1 = dx*df; fw2 = dy*df; fw3 = dz*df
         dumx = dumx + fw1
         dumy = dumy + fw2
         dumz = dumz + fw3
         f(1,j) = f(1,j) + fw1
         f(2,j) = f(2,j) + fw2
         f(3,j) = f(3,j) + fw3

         ! for "waterless" only
         df = df2*d2inv
         fw1 = dx*df; fw2 = dy*df; fw3 = dz*df
         coulx(j) = coulx(j) + fw1
         couly(j) = couly(j) + fw2
         coulz(j) = coulz(j) + fw3
         coulx(i) = coulx(i) - fw1
         couly(i) = couly(i) - fw2
         coulz(i) = coulz(i) - fw3
         df = SIX*( (f2-f1)-f1 )*d2inv
         fw1 = dx*df; fw2 = dy*df; fw3 = dz*df
         vdwx(j) = vdwx(j) + fw1
         vdwy(j) = vdwy(j) + fw2
         vdwz(j) = vdwz(j) + fw3
         vdwx(i) = vdwx(i) - fw1
         vdwy(i) = vdwy(i) - fw2
         vdwz(i) = vdwz(i) - fw3
 
10       continue
      end do
       
      ! pro/wat interactions
       
      if ( watres /= 0 ) then
       
      ic = ico(iaci+iac(watfirst)); cn1io = cn1(ic); cn2io = cn2(ic)
      cgio = cgi*cg(watfirst); cgih = cgi*cg(watfirst+1)
       
      ! stacking i <-> ow
       
      j = watfirst
      jfirst = 1; jlast = watres
      do jp = jfirst, jlast
         xwij1(jp) = xi - x(1,j); xwij2(jp) = yi - x(2,j); xwij3(jp) = zi - x(3,j)
         j = j + 3
      end do
       
      ! stacking i <-> hw1
       
      j = watfirst + 1
      jfirst = watres+1; jlast = 2*watres
      do jp = jfirst, jlast
         xwij1(jp) = xi - x(1,j); xwij2(jp) = yi - x(2,j); xwij3(jp) = zi - x(3,j)
         j = j + 3
      end do
       
      ! stacking i <-> hw2
       
      j = watfirst + 2
      jfirst = 2*watres+1; jlast = 3*watres
      do jp = jfirst, jlast
         xwij1(jp) = xi - x(1,j); xwij2(jp) = yi - x(2,j); xwij3(jp) = zi - x(3,j)
         j = j + 3
      end do
       
      do jp = 1, 3*watres
         rw1(jp) = ONE/( xwij1(jp)**2 + xwij2(jp)**2 + xwij3(jp)**2 )
      end do
       
      ! compute forces for i <-> ow
       
      j = watfirst
      jfirst = 1; jlast = watres
!c !DIR$ IVDEP
      do jp = jfirst, jlast
         df2 = -cgio*sqrt(rw1(jp))
         eelprowattmp = eelprowattmp + df2
         r6 = rw1(jp)**3; f2 = cn2io*r6; f1 = cn1io*(r6*r6)
         enbprowattmp = enbprowattmp + (f2-f1)

         !-- PB decomp
         if(idecomp == 1 .or. idecomp == 2) then
            if(eneopt == 1) then
               call decpair(1,i,j,df2)
            end if
            call decpair(2,i,j,-df2)
            call decpair(3,i,j,(f1 - f2))
         else if(idecomp == 3 .or. idecomp == 4) then
            call decpair(-2,i,j,-df2)
            call decpair(-3,i,j,(f1 - f2))
         end if

         df = ( df2+SIX*( (f2-f1)-f1 ) )*rw1(jp)
         fw1 = xwij1(jp)*df; fw2 = xwij2(jp)*df; fw3 = xwij3(jp)*df
         dumx = dumx + fw1
         dumy = dumy + fw2
         dumz = dumz + fw3
         fx(jp) = fw1
         fy(jp) = fw2
         fz(jp) = fw3

         j = j + 3
      end do
       
      ! compute forces for i <-> hw1, hw2
       
      j = watfirst + 1
      jfirst = watres+1; jlast = 3*watres
!c !DIR$ IVDEP
      do jp = jfirst, jlast
         df2 = -cgih*sqrt(rw1(jp))
         eelprowattmp = eelprowattmp + df2

         !-- PB decomp
         if(idecomp == 1 .or. idecomp == 2) then
            if(eneopt == 1) then
               call decpair(1,i,j,df2)
            end if
            call decpair(2,i,j,-df2)
         else if(idecomp == 3 .or. idecomp == 4) then
            call decpair(-2,i,j,-df2)
         end if

         df = ( df2 )*rw1(jp)
         fw1 = xwij1(jp)*df; fw2 = xwij2(jp)*df; fw3 = xwij3(jp)*df
         dumx = dumx + fw1
         dumy = dumy + fw2
         dumz = dumz + fw3
         fx(jp) = fw1
         fy(jp) = fw2
         fz(jp) = fw3

         j = j + 3
         if ( jp == 2*watres ) j = watfirst + 2 
      end do
       
      ! collecting forces
       
      j = watfirst
!c !DIR$ IVDEP
      do jp = 1, watres
         f(1,j  ) = f(1,j  ) + fx(         jp)
         f(2,j  ) = f(2,j  ) + fy(         jp)
         f(3,j  ) = f(3,j  ) + fz(         jp)
         f(1,j+1) = f(1,j+1) + fx(  watres+jp)
         f(2,j+1) = f(2,j+1) + fy(  watres+jp)
         f(3,j+1) = f(3,j+1) + fz(  watres+jp)
         f(1,j+2) = f(1,j+2) + fx(2*watres+jp)
         f(2,j+2) = f(2,j+2) + fy(2*watres+jp)
         f(3,j+2) = f(3,j+2) + fz(2*watres+jp)
         j = j + 3
      end do
       
      end if
       
      eelpro = eelpro - eelprotmp; enbpro = enbpro - enbprotmp
      eelprowat = eelprowat - eelprowattmp; enbprowat = enbprowat - enbprowattmp
       
      f(1,i) = f(1,i) - dumx
      f(2,i) = f(2,i) - dumy
      f(3,i) = f(3,i) - dumz
   end do  ! i = 1, ilast

!  open (unit = 107, file = 'coul.dat')

!  write(107,*) ' :::: Coulombic forces ::::'
!  do iatm = 1, natom
!     write(107,'(3e20.6)') coulx(iatm),couly(iatm),coulz(iatm)
!  end do

!  open (unit = 108, file = 'vdw.dat')

!  write(108,*) ' :::: Van der Waals forces ::::'
!  do iatm = 1, natom
!     write(108,'(3e20.6)') vdwx(iatm),vdwy(iatm),vdwz(iatm)
!  end do
!   
!  open (unit = 109, file = 'dir.dat')

!  write(109,*) ' :::: Direct forces ::::' 
!  do iatm = 1, natom
!     write(109,'(3e20.6)') f(1,iatm)-ofx(iatm),f(2,iatm)-ofy(iatm),f(3,iatm)-ofz(iatm)
!  end do

end subroutine pb_dirpro
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ compute direct water/water nonbonded forces 
subroutine pb_dirwat( watfirst,watres,ibgwat,ienwat,cn1oo,cn2oo,cgoo,cgoh,cghh,x,f,eelwat,enbwat )
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Authors:
   ! Lijiang Yang, Luo Research Group, UC-Irvine
   !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
   implicit none
    
#  include "pb_constants.h"

   ! Passed variables
    
   integer watfirst, watres, ibgwat, ienwat
   _REAL_ cn1oo, cn2oo, cgoo, cgoh, cghh
   _REAL_ x(3,*)
   _REAL_ eelwat, enbwat, f(3,*)
    
   ! Local variables
    
   integer i,j,ii,jj,jn,ni,ilast,natmo,jo,jh1,jh2
   _REAL_ dfa,dfb,dfc,r6,f1,f2,fw1,fw2,fw3,fw4,fw5,fw6,fw7,fw8,fw9
   _REAL_ dumxa,dumya,dumza,dumxb,dumyb,dumzb,dumxc,dumyc,dumzc
   _REAL_ xo1, xo2, xo3, xa1, xa2, xa3, xb1, xb2, xb3
   _REAL_ eelwattmp,enbwattmp
   _REAL_ xwij1(3*watres),xwij2(3*watres),xwij3(3*watres)
   _REAL_ xwij4(3*watres),xwij5(3*watres),xwij6(3*watres)
   _REAL_ xwij7(3*watres),xwij8(3*watres),xwij9(3*watres)
   _REAL_ rw1(3*watres),rw2(3*watres),rw3(3*watres)
   _REAL_ fx(3*watres),fy(3*watres),fz(3*watres)
    
   ilast = ibgwat + watres - 1
   natmo = watres   ! the number of water's oxygen atoms
   ni    = 1        ! ni is the counter of waters that have been considered
                    ! for every ii increment, it is incremented by 1
   i = watfirst - 3 ! the first water oxygen - 3
   do ii = ibgwat, ilast-1
       
      ! get the distances between each water(i)'s (O,H1,H2) and O(j). xwij(1:3,*): for pairs
      ! of O(i) and O(j); xwij(4:6,*): for H1(i) and O(j); xwji(7:9,*): for H2(i) and O(j).
      ! do the same for H1(j) and H2(j). the number of water(i)-O(j) pairs is natmo-ni.
      ! this is also the case for water(i)-H1(j) and for water(i)-H2(j).
       
      i = i + 3; jn = 1  ! jn is index of atom pairs
      xo1 = x(1,i  ); xo2 = x(2,i  ); xo3 = x(3,i  )
      xa1 = x(1,i+1); xa2 = x(2,i+1); xa3 = x(3,i+1)
      xb1 = x(1,i+2); xb2 = x(2,i+2); xb3 = x(3,i+2)
       
      jo = i; jh1 = i + 1; jh2 = i + 2
      do jj = ii+1, ilast
         jo = jo + 3
         xwij1(jn) = xo1 - x(1,jo)
         xwij2(jn) = xo2 - x(2,jo)
         xwij3(jn) = xo3 - x(3,jo)
         xwij4(jn) = xa1 - x(1,jo)
         xwij5(jn) = xa2 - x(2,jo)
         xwij6(jn) = xa3 - x(3,jo)
         xwij7(jn) = xb1 - x(1,jo)
         xwij8(jn) = xb2 - x(2,jo)
         xwij9(jn) = xb3 - x(3,jo)
         jn = jn + 1
      end do
      do jj = ii+1, ilast
         jh1 = jh1 + 3
         xwij1(jn) = xo1 - x(1,jh1)
         xwij2(jn) = xo2 - x(2,jh1)
         xwij3(jn) = xo3 - x(3,jh1)
         xwij4(jn) = xa1 - x(1,jh1)
         xwij5(jn) = xa2 - x(2,jh1)
         xwij6(jn) = xa3 - x(3,jh1)
         xwij7(jn) = xb1 - x(1,jh1)
         xwij8(jn) = xb2 - x(2,jh1)
         xwij9(jn) = xb3 - x(3,jh1)
         jn = jn + 1
      end do
      do jj = ii+1, ilast
         jh2 = jh2 + 3
         xwij1(jn) = xo1 - x(1,jh2)
         xwij2(jn) = xo2 - x(2,jh2)
         xwij3(jn) = xo3 - x(3,jh2)
         xwij4(jn) = xa1 - x(1,jh2)
         xwij5(jn) = xa2 - x(2,jh2)
         xwij6(jn) = xa3 - x(3,jh2)
         xwij7(jn) = xb1 - x(1,jh2)
         xwij8(jn) = xb2 - x(2,jh2)
         xwij9(jn) = xb3 - x(3,jh2)
         jn = jn + 1
      end do
      do j = 1, jn
         rw1(j) = ONE/( xwij1(j)**2 + xwij2(j)**2 + xwij3(j)**2 )
         rw2(j) = ONE/( xwij4(j)**2 + xwij5(j)**2 + xwij6(j)**2 )
         rw3(j) = ONE/( xwij7(j)**2 + xwij8(j)**2 + xwij9(j)**2 )
      end do
       
      ! compute the energy and force of the oxygen atom for water(i)
      ! if the ni'th water is under consideration, water after it will be natmo-ni
       
      eelwattmp = ZERO; enbwattmp = ZERO

      dumxa = ZERO; dumya = ZERO; dumza = ZERO
      dumxb = ZERO; dumyb = ZERO; dumzb = ZERO
      dumxc = ZERO; dumyc = ZERO; dumzc = ZERO
!c !DIR$ IVDEP
      do jn = 1, natmo-ni
         dfa = -cgoo*sqrt(rw1(jn))
         dfb = -cgoh*sqrt(rw2(jn))
         dfc = -cgoh*sqrt(rw3(jn))
         eelwattmp = eelwattmp + dfa+dfb+dfc
          
         r6 = rw1(jn)**3
         f2 = cn2oo*r6
         f1 = cn1oo*(r6*r6)
         enbwattmp = enbwattmp + (f2-f1)
          
         dfa = (dfa+SIX*((f2-f1)-f1))*rw1(jn)
         dfb = dfb*rw2(jn)
         dfc = dfc*rw3(jn)
          
         fw1 = xwij1(jn)*dfa
         fw2 = xwij2(jn)*dfa
         fw3 = xwij3(jn)*dfa
         fw4 = xwij4(jn)*dfb
         fw5 = xwij5(jn)*dfb
         fw6 = xwij6(jn)*dfb
         fw7 = xwij7(jn)*dfc
         fw8 = xwij8(jn)*dfc
         fw9 = xwij9(jn)*dfc
            
         dumxa = dumxa + fw1
         dumya = dumya + fw2
         dumza = dumza + fw3
         dumxb = dumxb + fw4
         dumyb = dumyb + fw5
         dumzb = dumzb + fw6
         dumxc = dumxc + fw7
         dumyc = dumyc + fw8
         dumzc = dumzc + fw9
            
         fx(jn) = fw1+fw4+fw7
         fy(jn) = fw2+fw5+fw8
         fz(jn) = fw3+fw6+fw9
      end do  !  jn = 1, natmo-ni 
        
      ! compute the energy and force of H1 and H2 atoms of water(i)
      ! there is natmo-ni water(i)-o(j) pairs, so H1 pairs will be from natmo-ni+1 to 2*(natmo-ni)
        
!c !DIR$ IVDEP
      do jn = natmo-ni+1, 2*(natmo-ni)
         dfa = -cgoh*sqrt(rw1(jn))
         dfb = -cghh*sqrt(rw2(jn))
         dfc = -cghh*sqrt(rw3(jn))
         eelwattmp = eelwattmp + dfa+dfb+dfc
            
         dfa = dfa*rw1(jn)
         dfb = dfb*rw2(jn)
         dfc = dfc*rw3(jn)
            
         fw1 = xwij1(jn)*dfa
         fw2 = xwij2(jn)*dfa
         fw3 = xwij3(jn)*dfa
         fw4 = xwij4(jn)*dfb
         fw5 = xwij5(jn)*dfb
         fw6 = xwij6(jn)*dfb
         fw7 = xwij7(jn)*dfc
         fw8 = xwij8(jn)*dfc
         fw9 = xwij9(jn)*dfc
            
         dumxa = dumxa + fw1
         dumya = dumya + fw2
         dumza = dumza + fw3
         dumxb = dumxb + fw4
         dumyb = dumyb + fw5
         dumzb = dumzb + fw6
         dumxc = dumxc + fw7
         dumyc = dumyc + fw8
         dumzc = dumzc + fw9
            
         fx(jn) = fw1+fw4+fw7
         fy(jn) = fw2+fw5+fw8
         fz(jn) = fw3+fw6+fw9
      end do  !  jn = natmo-ni+1, 2*(natmo-ni) 
        
!c !DIR$ IVDEP
      do jn = 2*(natmo-ni)+1, 3*(natmo-ni)
         dfa = -cgoh*sqrt(rw1(jn))
         dfb = -cghh*sqrt(rw2(jn))
         dfc = -cghh*sqrt(rw3(jn))
         eelwattmp = eelwattmp + dfa+dfb+dfc
            
         dfa = dfa*rw1(jn)
         dfb = dfb*rw2(jn)
         dfc = dfc*rw3(jn)
            
         fw1 = xwij1(jn)*dfa
         fw2 = xwij2(jn)*dfa
         fw3 = xwij3(jn)*dfa
         fw4 = xwij4(jn)*dfb
         fw5 = xwij5(jn)*dfb
         fw6 = xwij6(jn)*dfb
         fw7 = xwij7(jn)*dfc
         fw8 = xwij8(jn)*dfc
         fw9 = xwij9(jn)*dfc
            
         dumxa = dumxa + fw1
         dumya = dumya + fw2
         dumza = dumza + fw3
         dumxb = dumxb + fw4
         dumyb = dumyb + fw5
         dumzb = dumzb + fw6
         dumxc = dumxc + fw7
         dumyc = dumyc + fw8
         dumzc = dumzc + fw9
            
         fx(jn) = fw1+fw4+fw7
         fy(jn) = fw2+fw5+fw8
         fz(jn) = fw3+fw6+fw9
      end do  !  jn = 2*(natmo-ni)+1, 3*(natmo-ni)
       
      ! Now collecting forces ... for j's and i's
       
      jn = 1; j = i
!c !DIR$ IVDEP
      do jj = ii+1, ilast
         j = j + 3
         f(1,j  ) = f(1,j  ) + fx(             jn)
         f(2,j  ) = f(2,j  ) + fy(             jn)
         f(3,j  ) = f(3,j  ) + fz(             jn)
         f(1,j+1) = f(1,j+1) + fx(   natmo-ni +jn)
         f(2,j+1) = f(2,j+1) + fy(   natmo-ni +jn)
         f(3,j+1) = f(3,j+1) + fz(   natmo-ni +jn)
         f(1,j+2) = f(1,j+2) + fx(2*(natmo-ni)+jn)
         f(2,j+2) = f(2,j+2) + fy(2*(natmo-ni)+jn)
         f(3,j+2) = f(3,j+2) + fz(2*(natmo-ni)+jn)
         jn = jn + 1
      end do
       
      eelwat = eelwat - eelwattmp; enbwat = enbwat - enbwattmp
       
      f(1,i  ) = f(1,i  ) - dumxa
      f(2,i  ) = f(2,i  ) - dumya
      f(3,i  ) = f(3,i  ) - dumza
      f(1,i+1) = f(1,i+1) - dumxb
      f(2,i+1) = f(2,i+1) - dumyb
      f(3,i+1) = f(3,i+1) - dumzb
      f(1,i+2) = f(1,i+2) - dumxc
      f(2,i+2) = f(2,i+2) - dumyc
      f(3,i+2) = f(3,i+2) - dumzc
        
      ni = ni + 1
        
   end do  ! ii = ibgwat, ilast-1
    
    
end subroutine pb_dirwat

end subroutine pb_directnocut
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ pairwise nonbonded forces
subroutine pb_directwtcut( natom,atmlast,ifcap,idecomp,iprshrt,iar1pb,cn1pb,cn2pb,cn3pb,x,f,eel,enb )
    
   use poisson_boltzmann, only: outflag, ligand, multiblock, liveflag
   use decomp, only: decpair

   implicit none

#  include "pb_constants.h"
#  include "flocntrl.h"
    
   ! Passed variables
    
   integer natom, atmlast, ifcap, idecomp, iprshrt(*), iar1pb(4,0:natom)
   _REAL_ x(3,natom)
   _REAL_ cn1pb(*), cn2pb(*), cn3pb(*)
   _REAL_ f(3,natom)
   _REAL_ eel, enb
    
   ! Local variables
    
   integer i, j, jp, ilast, jfirst1, jlast1, jfirst2, jlast2
   _REAL_ xi, yi, zi
   _REAL_ dumx, dumy, dumz
   _REAL_ dx, dy, dz, d2inv, r6
   _REAL_ df2, f2, f1, df, fw1, fw2, fw3

   if ( do_pbdir == 0 ) return

   ! initialization
    
   eel     = ZERO; enb     = ZERO
   ilast   = atmlast - 1
   do i = 1, ilast
      if((ifcap == 2 .or. ifcap == 5) .and. outflag(i) == 1) cycle
      xi   = x(1,i); yi   = x(2,i); zi   = x(3,i)
       
      dumx = ZERO; dumy = ZERO; dumz = ZERO
       
      jfirst1 = iar1pb(1, i) + 1; jlast1 = iar1pb(2, i)
      jfirst2 = iar1pb(2, i) + 1; jlast2 = iar1pb(4, i)

      ! loop over direct coulombic pairs
       
!print *,"debug: jfirst1,jlast1,i",jfirst1,jlast1,i
      do jp = jfirst1, jlast1
         j = iprshrt(jp)
         if( j > atmlast ) cycle   ! Skip waters that should not be considered
         if((ifcap == 2 .or. ifcap == 5) .and. outflag(j) == 1 ) cycle   ! Skip ions that should not be considered
         if( multiblock ) then
             if( liveflag(j) == 0 ) cycle
         endif
         !if(multiblock .and. liveflag(j) == 0) cycle
         ! mjhsieh: in theory ligand should do the same thing, but it
         !          changes the result. please discuss.
         dx = xi - x(1,j); dy = yi - x(2,j); dz = zi - x(3,j)
         d2inv = ONE/(dx**2+dy**2+dz**2)
         df2 = cn3pb(jp)*sqrt(d2inv)
         eel = eel+df2

         !-- PB decomp
         if(idecomp == 1 .or. idecomp == 2) then
            call decpair(2,i,j,-df2)
         else if(idecomp == 3 .or. idecomp == 4) then
            call decpair(-2,i,j,-df2)
         end if

         r6 = d2inv**3
         f2 = cn2pb(jp)*r6
         f1 = cn1pb(jp)*(r6*r6)
         enb = enb + (f2-f1)
         df = (df2+SIX*((f2-f1)-f1))*d2inv

         !-- PB decomp
         if(idecomp == 1 .or. idecomp == 2) then
            call decpair(3,i,j,(f1 - f2))
         else if(idecomp == 3 .or. idecomp == 4) then
            call decpair(-3,i,j,(f1 - f2))
         end if

         fw1 = dx*df; fw2 = dy*df; fw3 = dz*df
         dumx = dumx + fw1
         dumy = dumy + fw2
         dumz = dumz + fw3
         f(1,j) = f(1,j) + fw1
         f(2,j) = f(2,j) + fw2
         f(3,j) = f(3,j) + fw3
      end do
       
      ! loop over nb pairs
       
!print *,"debug: jfirst1,jlast1,i",jfirst2,jlast2,i
      do jp = jfirst2, jlast2
         j = iprshrt(jp)
         if( j > atmlast ) cycle   ! Skip waters that should not be considered
         if((ifcap == 2 .or. ifcap == 5) .and. outflag(j) == 1) cycle   ! Skip ions that should not be considered
         if ( multiblock ) then
            if ( liveflag(j) == 0 ) cycle
         endif
         !if(multiblock .and. liveflag(j) == 0) cycle
         ! mjhsieh: in theory ligand should do the same thing, but it
         !          changes the result. please discuss.
         dx = xi - x(1,j); dy = yi - x(2,j); dz = zi - x(3,j)
         d2inv = ONE/(dx**2+dy**2+dz**2)
         r6 = d2inv**3
         f2 = cn2pb(jp)*r6
         f1 = cn1pb(jp)*(r6*r6)
         enb = enb + (f2-f1)
!print *,"debug: ENB =",enb,i,j
         df = SIX*((f2-f1)-f1)*d2inv

         !-- PB decomp
         if(idecomp == 1 .or. idecomp == 2) then
            call decpair(3,i,j,(f1 - f2))
         else if(idecomp == 3 .or. idecomp == 4) then
            call decpair(-3,i,j,(f1 - f2))
         end if

         fw1 = dx*df; fw2 = dy*df; fw3 = dz*df
         dumx = dumx + fw1
         dumy = dumy + fw2
         dumz = dumz + fw3
         f(1,j) = f(1,j) + fw1
         f(2,j) = f(2,j) + fw2
         f(3,j) = f(3,j) + fw3
      end do
!print *,"debug: EEL =",eel,i
!print *,"debug: ENB =",enb,i
      f(1,i) = f(1,i) - dumx
      f(2,i) = f(2,i) - dumy
      f(3,i) = f(3,i) - dumz
   end do  !  i = 1, ilast
 
   ! --- END OF PAIRS INVOLVING ATOM I ---
 
   eel = -eel; enb = -enb
 
end subroutine pb_directwtcut
