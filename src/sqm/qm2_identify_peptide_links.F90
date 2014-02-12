! <compile=optimized>
#include "copyright.h"
#include "assert.fh"
#include "dprec.fh"
subroutine qm2_identify_peptide_links(n_peptide_links,coord)

!Ross Walker (TSRI, 2005)
!Identifies peptide linkages based on distance. Also allocates
!the necessary memory for storing the identities of the
!peptide linkages. Note, at some point it would probably
!be better to use sander's bond information for this but
!for the moment this will suffice.

      use qmmm_module, only : qmmm_struct, qm2_struct
      implicit none

!Passed in
      integer, intent(out) :: n_peptide_links
      _REAL_, intent(in) :: coord(3,qmmm_struct%nquant_nlink)

!Local
      _REAL_, dimension(:), pointer :: RXYZ !Distance matrix used in finding peptide linkages
      integer :: l,i,j,k,jk,kj,ij,ji,kl,lk,m,mk,km
      integer :: ier=0

      ! Distance matrix is currently only needed for calculating peptide linkages
      ! Note: Ross Walker originally had a sqrt in this but I removed this for speed.
      ! Also - ultimately it might be better to use ambers bonding info for this.
      allocate (RXYZ(ishft(qmmm_struct%nquant_nlink*(qmmm_struct%nquant_nlink+1),-1)), stat=ier)
      REQUIRE(ier == 0)
     
      L=0
      do I=1,qmmm_struct%nquant_nlink
        do J=1,I
          L=L+1
          RXYZ(L)=((COORD(1,I)-COORD(1,J))**2+(COORD(2,I)-COORD(2,J))**2+(COORD(3,I)-COORD(3,J))**2)
        end do
      end do

!     IDENTIFY HOW MANY O=C-N-H SYSTEMS VIA THE INTERATOMIC DISTANCES MATRIX
      do i=1,qmmm_struct%nquant_nlink
        if(qmmm_struct%iqm_atomic_numbers(i) == 8) then
          do j=1,qmmm_struct%nquant_nlink
            if(qmmm_struct%iqm_atomic_numbers(j) == 6) then
              IJ=MAX(I,J)
              JI=I+J-IJ
!                                    RCW  1.69=1.3^2                                   
              if (RXYZ((IJ*(IJ-1))/2+JI) <= 1.69) then
                do k=1,qmmm_struct%nquant_nlink
                  if(qmmm_struct%iqm_atomic_numbers(k) == 7) then
                    JK=MAX(J,K)
                    KJ=J+K-JK
!                                               2.56=1.6^2
                    if (RXYZ((JK*(JK-1))/2+KJ) <= 2.56) then
                      do l=1,qmmm_struct%nquant_nlink
                        if(qmmm_struct%iqm_atomic_numbers(L) == 1) then
                          KL=MAX(K,L)
                          LK=K+L-KL
                          if (RXYZ((KL*(KL-1))/2+LK) <= 1.69) then
!   WE HAVE A H-N-C=O SYSTEM.  THE ATOM NUMBERS ARE L-K-J-I
!   NOW SEARCH OUT ATOM ATTACHED TO NITROGEN, THIS SPECIFIES
!   THE SYSTEM X-N-C=O
                            do M=1,qmmm_struct%nquant_nlink
                              if (M /= K .AND. M /= L .AND. M /= J) then
                                MK=MAX(M,K)
                                KM=M+K-MK
!                                                   2.89=1.7^2                          
                                if(RXYZ((MK*(MK-1))/2+KM) <= 2.89) then
                                  n_peptide_links=n_peptide_links+2
                                end if
                              end if
                            end do
                          end if
                        end if
                      end do
                    end if
                  end if
                end do
              end if
            end if
          end do
        end if
      end do

!Stage 2 allocate the identity array and fill it.
      allocate (qm2_struct%peptide_links(4,n_peptide_links), stat=ier)
      REQUIRE(ier==0)

      n_peptide_links=0
      do i=1,qmmm_struct%nquant_nlink
        if(qmmm_struct%iqm_atomic_numbers(i) == 8) then
          do j=1,qmmm_struct%nquant_nlink
            if(qmmm_struct%iqm_atomic_numbers(j) == 6) then
              IJ=MAX(I,J)
              JI=I+J-IJ
!                                    RCW  1.69=1.3^2                                   
              if (RXYZ((IJ*(IJ-1))/2+JI) <= 1.69) then
                do k=1,qmmm_struct%nquant_nlink
                  if(qmmm_struct%iqm_atomic_numbers(k) == 7) then
                    JK=MAX(J,K)
                    KJ=J+K-JK
!                                               2.56=1.6^2
                    if (RXYZ((JK*(JK-1))/2+KJ) <= 2.56) then
                      do l=1,qmmm_struct%nquant_nlink
                        if(qmmm_struct%iqm_atomic_numbers(L) == 1) then
                          KL=MAX(K,L)
                          LK=K+L-KL
                          if (RXYZ((KL*(KL-1))/2+LK) <= 1.69) then
!   WE HAVE A H-N-C=O SYSTEM.  THE ATOM NUMBERS ARE L-K-J-I
!   NOW SEARCH OUT ATOM ATTACHED TO NITROGEN, THIS SPECIFIES
!   THE SYSTEM X-N-C=O
                            do M=1,qmmm_struct%nquant_nlink
                              if (M /= K .AND. M /= L .AND. M /= J) then
                                MK=MAX(M,K)
                                KM=M+K-MK
!                                                   2.89=1.7^2                          
                                if(RXYZ((MK*(MK-1))/2+KM) <= 2.89) then
                                  n_peptide_links=n_peptide_links+1
                                  qm2_struct%peptide_links(1,n_peptide_links)=I
                                  qm2_struct%peptide_links(2,n_peptide_links)=J
                                  qm2_struct%peptide_links(3,n_peptide_links)=K
                                  qm2_struct%peptide_links(4,n_peptide_links)=M
                                  n_peptide_links=n_peptide_links+1
                                  qm2_struct%peptide_links(1,n_peptide_links)=I
                                  qm2_struct%peptide_links(2,n_peptide_links)=J
                                  qm2_struct%peptide_links(3,n_peptide_links)=K
                                  qm2_struct%peptide_links(4,n_peptide_links)=L
                                end if
                              end if
                            end do
                          end if
                        end if
                      end do
                    end if
                  end if
                end do
              end if
            end if
          end do
        end if
      end do

      deallocate ( RXYZ, stat=ier )
      REQUIRE(ier == 0)

      return

end subroutine qm2_identify_peptide_links


