subroutine GGF2_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,eGF,KA_sta)

! Compute the resonant part of the static BSE@GF2 matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  
! Local variables

  double precision              :: dem,num
  integer                       :: i,j,k,l
  integer                       :: a,b,c,d
  integer                       :: ia,jb

! Output variables

  double precision,intent(out)  :: KA_sta(nS,nS)

! Initialization

  KA_sta(:,:) = 0d0

! Second-order correlation kernel for the block A of the spinorbital manifold

    jb = 0

!$omp parallel do default(private) shared(KA_sta,ERI,num,dem,eGF,nO,nBas,eta,nC,nR)
  do j=nC+1,nO
    do b=nO+1,nBas-nR
      jb = (b-nO) + (j-1)*(nBas-nO)

      ia = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          ia = (a-nO) + (i-1)*(nBas-nO)

          do k=nC+1,nO
            do c=nO+1,nBas-nR
   
              dem = - (eGF(c) - eGF(k))
              num = ERI(j,k,i,c)*ERI(a,c,b,k) - ERI(j,k,i,c)*ERI(a,c,k,b) & 
                  - ERI(j,k,c,i)*ERI(a,c,b,k) + ERI(j,k,c,i)*ERI(a,c,k,b)

              KA_sta(ia,jb) = KA_sta(ia,jb) - num*dem/(dem**2 + eta**2)
          
              dem = + (eGF(c) - eGF(k))
              num = ERI(j,c,i,k)*ERI(a,k,b,c) - ERI(j,c,i,k)*ERI(a,k,c,b) & 
                  - ERI(j,c,k,i)*ERI(a,k,b,c) + ERI(j,c,k,i)*ERI(a,k,c,b)

              KA_sta(ia,jb) = KA_sta(ia,jb) + num*dem/(dem**2 + eta**2)
          
            end do
          end do

          do c=nO+1,nBas-nR
            do d=nO+1,nBas-nR
   
              dem = - (eGF(c) + eGF(d))
              num = ERI(a,j,c,d)*ERI(c,d,i,b) - ERI(a,j,c,d)*ERI(c,d,b,i) & 
                  - ERI(a,j,d,c)*ERI(c,d,i,b) + ERI(a,j,d,c)*ERI(c,d,b,i)

              KA_sta(ia,jb) = KA_sta(ia,jb) + 0.5d0*num*dem/(dem**2 + eta**2)
          
            end do
          end do

          do k=nC+1,nO
            do l=nC+1,nO

              dem = - (eGF(k) + eGF(l))
              num = ERI(a,j,k,l)*ERI(k,l,i,b) - ERI(a,j,k,l)*ERI(k,l,b,i) & 
                  - ERI(a,j,l,k)*ERI(k,l,i,b) + ERI(a,j,l,k)*ERI(k,l,b,i)

              KA_sta(ia,jb) = KA_sta(ia,jb) - 0.5d0*num*dem/(dem**2 + eta**2)
          
            end do
          end do

        end do
      end do

    end do
  end do
!$omp end parallel do

end subroutine 
