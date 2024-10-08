subroutine GGTpp_ppBSE_static_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nOO,nVV,nOOx,nVVx,lambda,Om1,rho1,Om2,rho2,TC)

! Compute the VVVV block of the static T-matrix

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
  integer,intent(in)            :: nOOx
  integer,intent(in)            :: nVVx
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: Om1(nVV)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: Om2(nOO)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: chi
  double precision              :: eps
  integer                       :: a,b,c,d,ab,cd,ef,mn

! Output variables

  double precision,intent(out)  :: TC(nVVx,nVVx)

  ab = 0
  do a=nO+1,nBas-nR
    do b=a+1,nBas-nR
      ab = ab + 1

      cd = 0
      do c=nO+1,nBas-nR
        do d=c+1,nBas-nR
          cd = cd + 1

          chi = 0d0

          do ef=1,nVV
            eps = + Om1(ef)
            chi = chi + rho1(a,b,ef)*rho1(c,d,ef)*eps/(eps**2 + eta**2)
          end do

          do mn=1,nOO
            eps = - Om2(mn)
            chi = chi + rho2(a,b,mn)*rho2(c,d,mn)*eps/(eps**2 + eta**2)
          end do

          TC(ab,cd) = lambda*chi

        end do
      end do

    end do
  end do

end subroutine 
