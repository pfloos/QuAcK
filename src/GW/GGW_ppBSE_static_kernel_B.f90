subroutine GGW_ppBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,lambda,ERI,Om,rho,KB)

! Compute the VVOO block of the static screening W for the pp-BSE

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
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: chi
  double precision              :: eps
  integer                       :: a,b,i,j,ab,ij,m

! Output variables

  double precision,intent(out)  :: KB(nVV,nOO)

! Initialization

  KB(:,:) = 0d0

!---------------!
! SpinOrb block !
!---------------!

  ab = 0
  do a=nO+1,nBas-nR
   do b=a+1,nBas-nR
      ab = ab + 1
      ij = 0
      do i=nC+1,nO
       do j=i+1,nO
          ij = ij + 1

          chi = 0d0
          do m=1,nS
            eps = Om(m)**2 + eta**2
            chi = chi - rho(i,a,m)*rho(j,b,m)*Om(m)/eps &
                      + rho(i,b,m)*rho(a,j,m)*Om(m)/eps
          end do

          KB(ab,ij) = 2d0*lambda*chi

        end do
      end do
    end do
  end do

end subroutine 
