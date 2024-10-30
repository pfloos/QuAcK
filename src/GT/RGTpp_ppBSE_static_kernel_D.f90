subroutine RGTpp_ppBSE_static_kernel_D(ispin,eta,nBas,nC,nO,nV,nR,nOO,nVV,lambda,eGF,Taaaa,Tabab,Tbaab,KD_sta)

! Compute the OOOO block of the static T-matrix 

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: Taaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Tabab(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Tbaab(nBas,nBas,nBas,nBas)

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,k,l,ij,kl,ef,mn,m,e

! Output variables

  double precision,intent(out)  :: KD_sta(nOO,nOO)

! Initialization

  KD_sta(:,:) = 0d0
  
!===============!
! singlet block !
!===============!

  if(ispin == 1) then

    ij = 0
    do i=nC+1,nO
      do j=i,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
          do l=k,nO
            kl = kl + 1
 
            do m=nC+1,nO
              do e=nO+1,nBas-nR
                 dem = eGF(m) - eGF(e)
                 ! Wabab_{ijkl}
                 num = Taaaa(i,m,k,e)*Tabab(e,j,m,l) + Tabab(i,m,k,e)*Taaaa(e,j,m,l)
                 KD_sta(ij,kl) = KD_sta(ij,kl) + num*dem/(dem**2 + eta**2)

                 num = Taaaa(i,e,k,m)*Tabab(m,j,e,l) + Tabab(i,e,k,m)*Taaaa(m,j,e,l)
                 KD_sta(ij,kl) = KD_sta(ij,kl) + num*dem/(dem**2 + eta**2)
                 
                 num = Taaaa(j,m,k,e)*Tabab(e,i,m,l) + Tabab(j,m,k,e)*Taaaa(e,i,m,l)
                 KD_sta(ij,kl) = KD_sta(ij,kl) + num*dem/(dem**2 + eta**2)

                 num = Taaaa(j,e,k,m)*Tabab(m,i,e,l) + Tabab(j,e,k,m)*Taaaa(m,i,e,l)
                 KD_sta(ij,kl) = KD_sta(ij,kl) + num*dem/(dem**2 + eta**2)

                 num = Tabab(i,m,k,e)*Tbaab(e,j,m,l) + Tbaab(i,e,k,m)*Tabab(m,j,e,l)
                 KD_sta(ij,kl) = KD_sta(ij,kl) - num*dem/(dem**2 + eta**2)
        
                 num = Tabab(j,m,k,e)*Tbaab(e,i,m,l) + Tbaab(j,e,k,m)*Tabab(m,i,e,l)
                 KD_sta(ij,kl) = KD_sta(ij,kl) - num*dem/(dem**2 + eta**2)
              end do
            end do

            KD_sta(ij,kl) = KD_sta(ij,kl) / sqrt((1d0 + Kronecker_delta(i,j))*(1d0 + Kronecker_delta(k,l)))
          end do
        end do

      end do
    end do

  end if

!===============!
! triplet block !
!===============!

  if(ispin == 2) then

    ij = 0
    do i=nC+1,nO
      do j=i+1,nO
        ij = ij + 1

        kl = 0
        do k=nC+1,nO
          do l=k+1,nO
            kl = kl + 1
 
            do m=nC+1,nO
              do e=nO+1,nBas-nR
                 dem = eGF(m) - eGF(e)
                 num = Taaaa(i,m,k,e)*Taaaa(e,j,m,l) + Tabab(i,m,k,e)*Tabab(e,j,m,l)
                           
                 KD_sta(ij,kl) = KD_sta(ij,kl) + num*dem/(dem**2 + eta**2)

                 num = Taaaa(i,e,k,m)*Taaaa(m,j,e,l) + Tabab(i,e,k,m)*Tabab(m,j,e,l)
                           
                 KD_sta(ij,kl) = KD_sta(ij,kl) + num*dem/(dem**2 + eta**2)
                 
                 num = Taaaa(j,m,k,e)*Taaaa(e,i,m,l) + Tabab(j,m,k,e)*Tabab(e,i,m,l)
                           
                 KD_sta(ij,kl) = KD_sta(ij,kl) - num*dem/(dem**2 + eta**2)

                 num = Taaaa(j,e,k,m)*Taaaa(m,i,e,l) + Tabab(j,e,k,m)*Tabab(m,i,e,l)
                           
                 KD_sta(ij,kl) = KD_sta(ij,kl) - num*dem/(dem**2 + eta**2)
                 
              end do
            end do

          end do
        end do

      end do
    end do

  end if

end subroutine 
