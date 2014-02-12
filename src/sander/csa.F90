#include "copyright.h"
#include "dprec.fh"

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ compute penalty function for residual csa (or pseudo-csa)
subroutine csa1( natom, x, f )
   
   
   use file_io_dat
   implicit none
   integer   natom
   _REAL_    x(*),f(*)
#  include "nmr.h"
#  include "extra.h"
   integer   iof, i, m, j, iset, k, p, q, ii,jj,kk
   _REAL_    dx, dy, dz, dr, temp, denom, temp1
   _REAL_    ccalc
   _REAL_    ecsai
   _REAL_    target

   _REAL_  dhn,dhn15,dhn2,dcn,dcn15,dcn2
   _REAL_  hn1,hn2,hn3,cn1,cn2,cn3,dc2n,dc3n
   _REAL_  dc(3,3), ddc(3,3,6), fc(9), s(3,3), fs(3,3)
   
   if (master .and. iprint /= 0) then
      rewind (58)
      write(58,9032)
      write(58,9031)
      write(58,9032)
   end if
   9031 format( &
         '     atom         curr. value target  deviation  penalty')
   9032 format(' ',78('-'))
   
   iof = 0
   do i=1,num_datasets
      s11(i)  = x(3*natom + iof + 1)
      s12(i)  = x(3*natom + iof + 2)
      s22(i)  = x(3*natom + iof + 3)
      s13(i)  = x(3*natom + iof + 4)
      s23(i)  = x(3*natom + iof + 5)
      s33(i)  = -s11(i) - s22(i)
      iof = iof + 5
   end do
   
   !    Loop over obvserved residual CSA couplings:
   
   do m=1,ncsa
      
      i = icsa(m)
      j = jcsa(m)
      k = kcsa(m)
      iset = dataset(m)  ! each "set" indexes a different alignment tensor
      iof = 5*(iset-1)
      !======================================================================
      
      !   The anisotropy of the alignment tensor is defined by
      !   S11,  S22, S33, S12, S23, S13, in the  coordinate
      !   system X={1,0,0}, Y={0,1,0}, Z={0,0,1}
      
      !   In order to have the order of magnitude of the S values be
      !   commensurate with coordinates, we will mutliply them by 10**5.
      
      !   ccalc: calculated resudual csa for spin i.
      !   cobsu: upper bound for the observed resudual csa for atom i;
      !         if calc. value is larger than this, a penalty will be imposed
      !   cobsl: lower bound for the observed resudual csa for atom i;
      !         if calc. value is smaller than this, a penalty will be imposed
      !   (No penalty if cobsl < ccalc < cobsu).
      !   ecsa: Resudual csa restrant energy term is E=cwt*(DJ-OJ)**2
      
      !======================================================================
      
      s(1,1) = s11(iset)*1.0d-5
      s(1,2) = s12(iset)*1.0d-5
      s(1,3) = s13(iset)*1.0d-5
      s(2,1) = s12(iset)*1.0d-5
      s(2,2) = s22(iset)*1.0d-5
      s(2,3) = s23(iset)*1.0d-5
      s(3,1) = s13(iset)*1.0d-5
      s(3,2) = s23(iset)*1.0d-5
      s(3,3) = s33(iset)*1.0d-5

      hn1 = x(3*i-2)-x(3*j-2)
      hn2 = x(3*i-1)-x(3*j-1)
      hn3 = x(3*i  )-x(3*j  )
      cn1 = x(3*k-2)-x(3*j-2)
      cn2 = x(3*k-1)-x(3*j-1)
      cn3 = x(3*k  )-x(3*j  )

      dhn2= 1.d0/(hn1**2 + hn2**2 + hn3**2)
      dhn = sqrt(dhn2)
      dhn15 = dhn*dhn2
      dcn2= 1.d0/(cn1**2 + cn2**2 + cn3**2)
      dcn = sqrt(dcn2)
      dcn15 = dcn*dcn2

      ! dc(,) will hold direction cosines, first index is local frame,
      !       second index is molecular frame

      dc(1,1) = hn1*dhn
      dc(1,2) = hn2*dhn
      dc(1,3) = hn3*dhn

      dc2n = 1.d0/Sqrt((hn1**2 + hn2**2 + hn3**2)* &
          (cn3**2*(hn1**2 + hn2**2) -  &
            2*cn1*cn3*hn1*hn3 -  &
            2*cn2*hn2*(cn1*hn1 + cn3*hn3) +  &
            cn2**2*(hn1**2 + hn3**2) +  &
            cn1**2*(hn2**2 + hn3**2)))
      dc(2,1) = (-(cn2*hn1*hn2) - cn3*hn1*hn3 + cn1*(hn2**2 + hn3**2))*dc2n
      dc(2,2) = (-(hn2*(cn1*hn1 + cn3*hn3)) + cn2*(hn1**2 + hn3**2))*dc2n
      dc(2,3) = (cn3*(hn1**2 + hn2**2) - (cn1*hn1 + cn2*hn2)*hn3)*dc2n

      dc3n = 1.d0/sqrt((cn2*hn1 - cn1*hn2)**2 + (-(cn3*hn1) + cn1*hn3)**2 + &
                       (cn3*hn2 - cn2*hn3)**2 )
      dc(3,1) = (cn3*hn2 - cn2*hn3)*dc3n
      dc(3,2) = (-(cn3*hn1) + cn1*hn3)*dc3n
      dc(3,3) = (cn2*hn1 - cn1*hn2)*dc3n

#if 0
      ! check orthonormality of dc:

      do ii=1,3
        do jj=1,3
          fs(ii,jj) = 0.d0
          do kk=1,3
             fs(ii,jj) = fs(ii,jj) + dc(ii,kk)*dc(jj,kk)
          end do
        end do
      end do
      write(6,*) 'dc orthonorm: ', i,j
      write(6,'(3f20.10)') fs(1,1), fs(1,2), fs(1,3)
      write(6,'(3f20.10)') fs(2,1), fs(2,2), fs(2,3)
      write(6,'(3f20.10)') fs(3,1), fs(3,2), fs(3,3)
      write(6,*) 'dc itself: '
      write(6,'(3f20.10)') dc(1,1), dc(1,2), dc(1,3)
      write(6,'(3f20.10)') dc(2,1), dc(2,2), dc(2,3)
      write(6,'(3f20.10)') dc(3,1), dc(3,2), dc(3,3)
#endif


      ! ddc(i,j,k) will hold D[ dc[i,j], k], where k=1,2,3 is hn1,hn2,hn3
      !            and k=4,5,6 is cn1,cn2,cn3

      ddc(1,1,1) = (hn2**2 + hn3**2)*dhn15
      ddc(1,1,2) = -((hn1*hn2)*dhn15)
      ddc(1,1,3) = -((hn1*hn3)*dhn15)
      ddc(1,1,4) = 0.d0
      ddc(1,1,5) = 0.d0
      ddc(1,1,6) = 0.d0

      ddc(1,2,1) = -((hn1*hn2)*dhn15)
      ddc(1,2,2) = (hn1**2 + hn3**2)*dhn15
      ddc(1,2,3) = -((hn2*hn3)*dhn15)
      ddc(1,2,4) = 0.d0
      ddc(1,2,5) = 0.d0
      ddc(1,2,6) = 0.d0

      ddc(1,3,1) = -((hn1*hn3)*dhn15)
      ddc(1,3,2) = -((hn2*hn3)*dhn15)
      ddc(1,3,3) = (hn1**2 + hn2**2)*dhn15
      ddc(1,3,4) = 0.d0
      ddc(1,3,5) = 0.d0
      ddc(1,3,6) = 0.d0

      denom = (hn1**2 + hn2**2 + hn3**2)*  &
             (cn3**2*(hn1**2 + hn2**2) - 2*cn1*cn3*hn1*hn3 - &
               2*cn2*hn2*(cn1*hn1 + cn3*hn3) + &
               cn2**2*(hn1**2 + hn3**2) + cn1**2*(hn2**2 + hn3**2))
      denom = 1.d0/(denom*sqrt(denom))

      ddc(2,1,1) = 0.5d0*(2*(-(cn2*hn2) - cn3*hn3)* &
           (hn1**2 + hn2**2 + hn3**2)* (cn3**2*(hn1**2 + hn2**2) -  &
             2*cn1*cn3*hn1*hn3 - 2*cn2*hn2*(cn1*hn1 + cn3*hn3) +  &
             cn2**2*(hn1**2 + hn3**2) + cn1**2*(hn2**2 + hn3**2)) -  &
          2*(-(cn2*hn1*hn2) - cn3*hn1*hn3 + cn1*(hn2**2 + hn3**2))* &
           ((2*cn3*hn1 - cn1*hn3)* (cn3*(hn1**2 + hn2**2) -  &
                (cn1*hn1 + cn2*hn2)*hn3) + (2*cn2*hn1 - cn1*hn2)* &
              (-(hn2*(cn1*hn1 + cn3*hn3)) + cn2*(hn1**2 + hn3**2)) +  &
             (cn2*hn2 + cn3*hn3)* (cn2*hn1*hn2 + cn3*hn1*hn3 -  &
                cn1*(hn2**2 + hn3**2)))) * denom

      ddc(2,1,2) = (-(cn2**3*hn1*(hn1**2 + hn3**2)**2) +  &
          hn2*(cn1*hn1 + cn3*hn3)* &
           (cn1**2*hn1*(hn2**2 + hn3**2) - cn1*cn3*hn3* &
              (3*hn1**2 + hn2**2 + hn3**2) +  &
             cn3**2*hn1*(2*hn1**2 + 2*hn2**2 + hn3**2)) + cn2**2*hn2* &
           (cn3*hn1*hn3* (2*hn1**2 - hn2**2 + 2*hn3**2) +  &
             cn1*(3*hn1**4 + 4*hn1**2*hn3**2 + hn2**2*hn3**2 + hn3**4)) -  &
          cn2*(3*cn1**2*hn1**3*hn2**2 + cn3**2*hn1* &
            (hn1**4 - hn2**4 + 2*hn1**2*hn3**2 + 3*hn2**2*hn3**2 + hn3**4) - &
             cn1*cn3*hn3* (hn1**4 - hn2**4 + hn3**4 +  &
                hn1**2*(-6*hn2**2 + 2*hn3**2)))) * denom

      ddc(2,1,3) = 0.5d0*(2*(-(cn3*hn1) + 2*cn1*hn3)* &
           (hn1**2 + hn2**2 + hn3**2)* (cn3**2*(hn1**2 + hn2**2) -  &
             2*cn1*cn3*hn1*hn3 - 2*cn2*hn2*(cn1*hn1 + cn3*hn3) +  &
             cn2**2*(hn1**2 + hn3**2) + cn1**2*(hn2**2 + hn3**2)) -  &
          2*(-(cn2*hn1*hn2) - cn3*hn1*hn3 + cn1*(hn2**2 + hn3**2))* &
           ((cn1*hn1 + cn2*hn2)* (-(cn3*(hn1**2 + hn2**2)) +  &
                (cn1*hn1 + cn2*hn2)*hn3) + (cn3*hn2 - 2*cn2*hn3)* &
              (hn2*(cn1*hn1 + cn3*hn3) - cn2*(hn1**2 + hn3**2)) +  &
             (cn3*hn1 - 2*cn1*hn3)* (cn2*hn1*hn2 + cn3*hn1*hn3 -  &
                cn1*(hn2**2 + hn3**2)))) * denom

      ddc(2,1,4) = ((cn3*hn2 - cn2*hn3)**2*  &
          (hn1**2 + hn2**2 + hn3**2)**2) * denom
      ddc(2,1,5) = -((cn3*hn1 - cn1*hn3)*(cn3*hn2 - cn2*hn3)* &
            (hn1**2 + hn2**2 + hn3**2)**2) * denom
      ddc(2,1,6) = -((cn2*hn1 - cn1*hn2)*(-(cn3*hn2) + cn2*hn3)*  &
            (hn1**2 + hn2**2 + hn3**2)**2) * denom

      ddc(2,2,1) = (-(cn1**3*hn2*(hn2**2 + hn3**2)**2) +  &
          hn1*(cn2*hn2 + cn3*hn3)* (cn2**2*hn2*(hn1**2 + hn3**2) +  &
             cn3**2*hn2* (2*hn1**2 + 2*hn2**2 + hn3**2) -  &
             cn2*cn3*hn3*(hn1**2 + 3*hn2**2 + hn3**2)) + cn1**2*hn1* &
           (cn2*(3*hn2**4 + hn1**2*hn3**2 + 4*hn2**2*hn3**2 + hn3**4) +  &
             cn3*hn2*hn3* (-hn1**2 + 2*(hn2**2 + hn3**2))) +  &
          cn1*(-3*cn2**2*hn1**2*hn2**3 + cn2*cn3*hn3* &
              (-hn1**4 - 6*hn1**2*hn2**2 + (hn2**2 + hn3**2)**2) -  &
             cn3**2*hn2* (-hn1**4 + 3*hn1**2*hn3**2 +  &
                (hn2**2 + hn3**2)**2))) * denom

      ddc(2,2,2) = 0.5d0*(2*(-(cn1*hn1) - cn3*hn3)* &
           (hn1**2 + hn2**2 + hn3**2)* (cn3**2*(hn1**2 + hn2**2) -  &
             2*cn1*cn3*hn1*hn3 - 2*cn2*hn2*(cn1*hn1 + cn3*hn3) +  &
             cn2**2*(hn1**2 + hn3**2) + cn1**2*(hn2**2 + hn3**2)) -  &
          2*(-(hn2*(cn1*hn1 + cn3*hn3)) + cn2*(hn1**2 + hn3**2))* &
           ((2*cn3*hn2 - cn2*hn3)* (cn3*(hn1**2 + hn2**2) -  &
                (cn1*hn1 + cn2*hn2)*hn3) + (cn1*hn1 + cn3*hn3)* &
              (hn2*(cn1*hn1 + cn3*hn3) - cn2*(hn1**2 + hn3**2)) +  &
             (cn2*hn1 - 2*cn1*hn2)* (cn2*hn1*hn2 + cn3*hn1*hn3 -  &
                cn1*(hn2**2 + hn3**2)))) * denom

      ddc(2,2,3) = (-(cn3**3*hn2*(hn1**2 + hn2**2)**2) -  &
          3*cn2**2*cn3*hn2**3*hn3**2 + cn2**3*hn2**2*hn3*(hn1**2 + hn3**2) + &
          cn1**3*hn1*hn2*hn3* (hn1**2 + 2*(hn2**2 + hn3**2)) +  &
          cn2*cn3**2*hn3* (hn1**4 + 3*hn2**4 +  &
             hn1**2*(4*hn2**2 + hn3**2)) + cn1*hn1*(-3*cn2**2*hn2**3*hn3 +  &
             cn3**2*hn2*hn3* (2*hn1**2 + 2*hn2**2 - hn3**2) +  &
             cn2*cn3*(hn1**4 + 2*hn1**2*hn2**2 +  &
                hn2**4 - 6*hn2**2*hn3**2 - hn3**4)) - cn1**2*(cn2*hn3* &
              (hn1**4 - 2*hn2**2*(hn2**2 + hn3**2) +  &
                hn1**2*(2*hn2**2 + hn3**2)) +  &
             cn3*hn2*(hn1**4 + hn2**4 - hn3**4 +  &
                hn1**2*(2*hn2**2 + 3*hn3**2)))) * denom

      ddc(2,2,4) = -((cn3*hn1 - cn1*hn3)*(cn3*hn2 - cn2*hn3)* &
            (hn1**2 + hn2**2 + hn3**2)**2) * denom
      ddc(2,2,5) = ((cn3*hn1 - cn1*hn3)**2* &
          (hn1**2 + hn2**2 + hn3**2)**2) * denom
      ddc(2,2,6) = -((cn2*hn1 - cn1*hn2)*(cn3*hn1 - cn1*hn3)* &
            (hn1**2 + hn2**2 + hn3**2)**2) * denom

      ddc(2,3,1) = (-(cn1**3*hn3*(hn2**2 + hn3**2)**2) +  &
          hn1*(cn2*hn2 + cn3*hn3)* (cn3**2*(hn1**2 + hn2**2)*hn3 +  &
             cn2**2*hn3* (2*hn1**2 + hn2**2 + 2*hn3**2) -  &
             cn2*cn3*hn2*(hn1**2 + hn2**2 + 3*hn3**2)) + cn1**2*hn1* &
           (cn3*(hn1**2*hn2**2 + hn2**4 + 4*hn2**2*hn3**2 + 3*hn3**4) +  &
             cn2*hn2*hn3* (-hn1**2 + 2*(hn2**2 + hn3**2))) +  &
          cn1*(-3*cn3**2*hn1**2*hn3**3 - cn2**2*hn3* &
              (-hn1**4 + 3*hn1**2*hn2**2 + (hn2**2 + hn3**2)**2) +  &
             cn2*cn3*hn2* (-hn1**4 - 6*hn1**2*hn3**2 +  &
                (hn2**2 + hn3**2)**2))) * denom

      ddc(2,3,2) = (cn3**3*hn2*(hn1**2 + hn2**2)*hn3**2 -  &
          3*cn2*cn3**2*hn2**2*hn3**3 - cn2**3*hn3*(hn1**2 + hn3**2)**2 +  &
          cn1**3*hn1*hn2*hn3* (hn1**2 + 2*(hn2**2 + hn3**2)) +  &
          cn2**2*cn3*hn2* (hn1**4 + 3*hn3**4 +  &
             hn1**2*(hn2**2 + 4*hn3**2)) + cn1*hn1*(-3*cn3**2*hn2*hn3**3 +  &
             cn2**2*hn2*hn3* (2*hn1**2 - hn2**2 + 2*hn3**2) +  &
             cn2*cn3*(hn1**4 - hn2**4 + 2*hn1**2*hn3**2 - 6*hn2**2*hn3**2 +  &
                hn3**4)) - cn1**2*(cn3*hn2* &
              (hn1**4 - 2*hn3**2*(hn2**2 + hn3**2) +  &
                hn1**2*(hn2**2 + 2*hn3**2)) +  &
             cn2*hn3*(hn1**4 - hn2**4 + hn3**4 +  &
                hn1**2*(3*hn2**2 + 2*hn3**2)))) * denom

      ddc(2,3,3) = 0.5d0*(2*(-(cn1*hn1) - cn2*hn2)* &
           (hn1**2 + hn2**2 + hn3**2)* (cn3**2*(hn1**2 + hn2**2) -  &
             2*cn1*cn3*hn1*hn3 - 2*cn2*hn2*(cn1*hn1 + cn3*hn3) +  &
             cn2**2*(hn1**2 + hn3**2) + cn1**2*(hn2**2 + hn3**2)) -  &
          2*(cn3*(hn1**2 + hn2**2) - (cn1*hn1 + cn2*hn2)*hn3)* &
           ((cn1*hn1 + cn2*hn2)* (-(cn3*(hn1**2 + hn2**2)) +  &
                (cn1*hn1 + cn2*hn2)*hn3) + (cn3*hn2 - 2*cn2*hn3)* &
              (hn2*(cn1*hn1 + cn3*hn3) - cn2*(hn1**2 + hn3**2)) +  &
             (cn3*hn1 - 2*cn1*hn3)* (cn2*hn1*hn2 + cn3*hn1*hn3 -  &
                cn1*(hn2**2 + hn3**2)))) * denom

      ddc(2,3,4) = -((cn2*hn1 - cn1*hn2)*(-(cn3*hn2) + cn2*hn3)* &
            (hn1**2 + hn2**2 + hn3**2)**2) * denom
      ddc(2,3,5) = -((cn2*hn1 - cn1*hn2)*(cn3*hn1 - cn1*hn3)* &
            (hn1**2 + hn2**2 + hn3**2)**2) * denom
      ddc(2,3,6) = ((cn2*hn1 - cn1*hn2)**2* &
          (hn1**2 + hn2**2 + hn3**2)**2) * denom

      denom = (cn2*hn1 - cn1*hn2)**2 + (cn3*hn1 - cn1*hn3)**2 + &
              (cn3*hn2 - cn2*hn3)**2
      denom = 1.d0/(denom*sqrt(denom))

      ddc(3,1,1) = -((cn3*hn2 - cn2*hn3)* (cn2**2*hn1 - cn1*cn2*hn2 + &
              cn3*(cn3*hn1 - cn1*hn3))) * denom
      ddc(3,1,2) = ((cn3*hn1 - cn1*hn3)* (cn2**2*hn1 - cn1*cn2*hn2 + &
            cn3*(cn3*hn1 - cn1*hn3))) * denom
      ddc(3,1,3) = -((cn2*hn1 - cn1*hn2)* (cn2**2*hn1 - cn1*cn2*hn2 + &
              cn3*(cn3*hn1 - cn1*hn3))) * denom
      ddc(3,1,4) = ((cn3*hn2 - cn2*hn3)* (cn2*hn1*hn2 + cn3*hn1*hn3 - &
            cn1*(hn2**2 + hn3**2))) * denom
      ddc(3,1,5) = -((cn3*hn1 - cn1*hn3)* (cn2*hn1*hn2 + cn3*hn1*hn3 - &
              cn1*(hn2**2 + hn3**2))) * denom
      ddc(3,1,6) = ((cn2*hn1 - cn1*hn2)* (cn2*hn1*hn2 + cn3*hn1*hn3 - &
            cn1*(hn2**2 + hn3**2))) * denom

      ddc(3,2,1) = ((-(cn3*hn2) + cn2*hn3)* (-(cn1*cn2*hn1) + cn1**2*hn2 + &
            cn3*(cn3*hn2 - cn2*hn3))) * denom
      ddc(3,2,2) = ((cn3*hn1 - cn1*hn3)* (-(cn1*cn2*hn1) + cn1**2*hn2 + &
            cn3*(cn3*hn2 - cn2*hn3))) * denom
      ddc(3,2,3) = ((cn2*hn1 - cn1*hn2)* (cn1*cn2*hn1 - cn1**2*hn2 + &
            cn3*(-(cn3*hn2) + cn2*hn3))) * denom
      ddc(3,2,4) = ((cn3*hn2 - cn2*hn3)* (hn2*(cn1*hn1 + cn3*hn3) - &
            cn2*(hn1**2 + hn3**2))) * denom
      ddc(3,2,5) = -((cn3*hn1 - cn1*hn3)* (hn2*(cn1*hn1 + cn3*hn3) - &
              cn2*(hn1**2 + hn3**2))) * denom
      ddc(3,2,6) = -((cn2*hn1 - cn1*hn2)* (-(hn2*(cn1*hn1 + cn3*hn3)) + &
              cn2*(hn1**2 + hn3**2))) * denom

      ddc(3,3,1) = ((cn3*hn2 - cn2*hn3)* (cn1*cn3*hn1 - cn1**2*hn3 + &
            cn2*(cn3*hn2 - cn2*hn3))) * denom
      ddc(3,3,2) = -((cn3*hn1 - cn1*hn3)* (cn1*cn3*hn1 - cn1**2*hn3 + &
              cn2*(cn3*hn2 - cn2*hn3))) * denom
      ddc(3,3,3) = -((cn2*hn1 - cn1*hn2)* (-(cn1*cn3*hn1) + cn1**2*hn3 + &
              cn2*(-(cn3*hn2) + cn2*hn3))) * denom
      ddc(3,3,4) = -((cn3*hn2 - cn2*hn3)* (cn3*(hn1**2 + hn2**2) - &
              (cn1*hn1 + cn2*hn2)*hn3)) * denom
      ddc(3,3,5) = ((cn3*hn1 - cn1*hn3)* (cn3*(hn1**2 + hn2**2) - &
            (cn1*hn1 + cn2*hn2)*hn3)) * denom
      ddc(3,3,6) = ((cn2*hn1 - cn1*hn2)* (-(cn3*(hn1**2 + hn2**2)) + &
            (cn1*hn1 + cn2*hn2)*hn3)) * denom

      ccalc = 0.d0
      do kk=1,9
         fc(kk) = 0.d0
      end do
      do p=1,3
         do q = 1,3
            fs(p,q) = 0.d0
            do ii=1,3
               do jj=1,3
                  ccalc = ccalc + s(p,q)*dc(ii,p)*sigma(ii,jj,m)*dc(jj,q)
                  fs(p,q) = fs(p,q) + dc(ii,p)*sigma(ii,jj,m)*dc(jj,q)
                  do kk=1,3
                     temp = s(p,q)*sigma(ii,jj,m)*(dc(ii,p)*ddc(jj,q,kk) + &
                           ddc(ii,p,kk)*dc(jj,q))
                     fc(kk) = fc(kk) + temp
                     fc(kk+3) = fc(kk+3) - temp
                     temp = s(p,q)*sigma(ii,jj,m)*(dc(ii,p)*ddc(jj,q,kk+3) + &
                           ddc(ii,p,kk+3)*dc(jj,q))
                     fc(kk+6) = fc(kk+6) + temp
                     fc(kk+3) = fc(kk+3) - temp
                  end do
               end do
            end do
         end do
      end do
      ccalc = (2.d0/3.d0)*ccalc   !  value in ppm
      ccalc = ccalc*field(m)      !  convert from ppm to Hz
      
      if( ccalc > cobsu(m) ) then
         ecsai = cwt(m)*(ccalc - cobsu(m))**2
      else if( ccalc < cobsl(m) ) then
         ecsai = cwt(m)*(ccalc - cobsl(m))**2
      else
         ecsai = 0.d0
      end if
      ecsa = ecsa + ecsai
      
      if( ccalc > cobsu(m) ) then
         temp1=-4.0*cwt(m)*(ccalc-cobsu(m))*field(m)/3.0
      else if( ccalc < cobsl(m) ) then
         temp1=-4.0*cwt(m)*(ccalc-cobsl(m))*field(m)/3.0
      else
         temp1 = 0.d0
      end if
      
      !---   Derivatives with respect to Sij

      if (ifreezes .eq. 0 ) then
         f(3*natom+iof+1) = f(3*natom+iof+1) + temp1*(fs(1,1)-fs(3,3))*1.0d-5
         f(3*natom+iof+3) = f(3*natom+iof+3) + temp1*(fs(2,2)-fs(3,3))*1.0d-5
      
         f(3*natom+iof+2) = f(3*natom+iof+2) + temp1*(fs(1,2)+fs(2,1))*1.0d-5
         f(3*natom+iof+4) = f(3*natom+iof+4) + temp1*(fs(1,3)+fs(3,1))*1.0d-5
         f(3*natom+iof+5) = f(3*natom+iof+5) + temp1*(fs(2,3)+fs(3,2))*1.0d-5
      end if
      
      !---   Derivatives with respect to x,y,z
      
      f(3*i-2) = f(3*i-2) + temp1*fc(1)
      f(3*i-1) = f(3*i-1) + temp1*fc(2)
      f(3*i  ) = f(3*i  ) + temp1*fc(3)
      f(3*j-2) = f(3*j-2) + temp1*fc(4)
      f(3*j-1) = f(3*j-1) + temp1*fc(5)
      f(3*j  ) = f(3*j  ) + temp1*fc(6)
      f(3*k-2) = f(3*k-2) + temp1*fc(7)
      f(3*k-1) = f(3*k-1) + temp1*fc(8)
      f(3*k  ) = f(3*k  ) + temp1*fc(9)

      if ( master .and. iprint /= 0 ) then
         if ( ecsai > dcut ) then
            target = 0.5d0*(cobsu(m)+cobsl(m))
            write(58,9073) resat(j)(1:13), &
                ccalc,target,ccalc-target,ecsai
         end if
      end if
      9073 format(' ',a13,' :',4f9.3)
      
   end do  !  m=1,ncsa

   !     ---printing of results:
   
   if (master .and. iprint /= 0) then
      write(58,44) ecsa
      44 format(40x,'Total csa    constraint:',f8.2)
   end if
   
   return
end subroutine csa1 


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ [Enter a one-line description of subroutine csaread here]
subroutine csaread
   
   use file_io_dat
   implicit none
#  include "nmr.h"
#  include "md.h"
   character(len=80) line
   integer   iin, ifind, i, iof
   _REAL_ sigma11(maxcsa), sigma22(maxcsa), sigma12(maxcsa), &
          sigma13(maxcsa), sigma23(maxcsa)

   namelist /csa/ ncsa, icsa, jcsa, kcsa, cobsu, cobsl, datasetc, &
         sigma11, sigma22, sigma12, sigma13, sigma23, field, cwt, ccut
   
   ! If restraint input has been redirected, open the appropriate file
   
   call amopen(38,redir(9)(1:iredir(9)),'O','F','R')
   iin = 38
   write(6,10) redir(9)(1:iredir(9))
   10 format(' CSA info will be read from file: ',a)
   
   !  --- read and echo title from csa file:
   
   write(6,*) 'Here are comments from the csa input file:'
   42 read(iin,'(a)') line
   if (line(1:1) == '#') then
      write(6,*) line
      goto 42
   end if
   backspace (iin)
   write(6,*)
   
   ! read the namelist csa, first setting up defaults:
   
   ccut = 0.1
   do i=1,maxcsa
      cwt(i) = 1.0d0
      datasetc(i) = 1
      cobsu(i) = 0.d0
      cobsl(i) = 0.d0
   end do
   read(iin,nml=csa,err=30)
   goto 40
   30 write(6,*) 'namelist reports error reading &csa'
   call mexit(6,1)
   40 continue
   if( ncsa > maxcsa) then
      write(6,*) 'ncsa is too big: ',ncsa,maxcsa
      call mexit(6,1)
   end if

   do i=1,ncsa
      sigma(1,1,i) = sigma11(i)
      sigma(2,2,i) = sigma22(i)
      sigma(3,3,i) = -sigma11(i) - sigma22(i)
      sigma(1,2,i) = sigma12(i)
      sigma(2,1,i) = sigma12(i)
      sigma(1,3,i) = sigma13(i)
      sigma(3,1,i) = sigma13(i)
      sigma(2,3,i) = sigma23(i)
      sigma(3,2,i) = sigma23(i)
   end do
   
   return
end subroutine csaread 
