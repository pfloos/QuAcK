subroutine GGW_ppBSE_dynamic_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,lambda,eGW,Om,rho,KB_dyn)

! Compute the dynamic part of the Bethe-Salpeter equation matrices

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  
! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: dem,num
  integer                       :: m
  integer                       :: a,b,i,j
  integer                       :: ab,ij

! Output variables

  double precision,intent(out)  :: KB_dyn(nVV,nOO)

! Initialization

  KB_dyn(:,:) = 0d0

! Build dynamic B matrix

  ab = 0
  do a=nO+1,nBas-nR
    do b=a+1,nBas-nR
      ab = ab + 1
 
      ij = 0
      do i=nC+1,nO
        do j=i+1,nO     
          ij = ij + 1
  
          do m=1,nS
             
               num = (rho(a,i,m)*rho(b,j,m) - rho(b,i,m)*rho(a,j,m))/2
               dem = - Om(m) - eGW(b) + eGW(j)
               KB_dyn(ab,ij) = KB_dyn(ab,ij) + num*dem/(dem**2 + eta**2)

               dem = - Om(m) - eGW(a) + eGW(i)
               KB_dyn(ab,ij) = KB_dyn(ab,ij) + num*dem/(dem**2 + eta**2)

               dem = - Om(m) - eGW(a) + eGW(j)
               KB_dyn(ab,ij) = KB_dyn(ab,ij) + num*dem/(dem**2 + eta**2)
       
               dem = - Om(m) - eGW(b) + eGW(i)
               KB_dyn(ab,ij) = KB_dyn(ab,ij) + num*dem/(dem**2 + eta**2)

          end do            

          KB_dyn(ab,ij) = 0.5d0*KB_dyn(ab,ij)

        end do
      end do

    end do
  end do

end subroutine 
