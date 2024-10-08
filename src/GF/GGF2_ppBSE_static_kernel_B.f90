subroutine GGF2_ppBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nOO,nVV,lambda,ERI,eGF,KB_sta)

! Compute the resonant part of the dynamic BSE@GF2 matrix

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGF(nBas)
  
! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  integer                       :: i,j,a,b,m,e
  integer                       :: ab,ij

! Output variables

  double precision,intent(out)  :: KB_sta(nVV,nOO)

! Initialization

  KB_sta(:,:) = 0d0

! Second-order correlation kernel for the block B of the spinorbital manifold

  ab = 0
  do a=nO+1,nBas-nR
    do b=a+1,nBas-nR
      ab = ab + 1

      ij = 0
      do i=nC+1,nO
        do j=i+1,nO
          ij = ij + 1

          do m=nC+1,nO
            do e=nO+1,nBas-nR
   
              dem = eGF(m) - eGF(e)
              num =       (ERI(a,m,i,e) - ERI(a,m,e,i)) * (ERI(e,b,m,j) - ERI(e,b,j,m))
              num = num + (ERI(a,e,i,m) - ERI(a,e,m,i)) * (ERI(m,b,e,j) - ERI(m,b,j,e))
              num = num - (ERI(b,m,i,e) - ERI(b,m,e,i)) * (ERI(e,a,m,j) - ERI(e,a,j,m))
              num = num - (ERI(b,e,i,m) - ERI(b,e,m,i)) * (ERI(m,a,e,j) - ERI(m,a,j,e))                                                                                
              KB_sta(ab,ij) = KB_sta(ab,ij) + num*dem/(dem**2 + eta**2)
              
            end do
          end do

        end do
      end do

    end do
  end do

end subroutine 
