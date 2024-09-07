subroutine GGW_ppBSE_static_kernel_C(eta,nBas,nC,nO,nV,nR,nS,nVV,lambda,ERI,Om,rho,KC)

! Compute the VVVV block of the static screening W for the pp-BSE

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
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
  double precision              :: tmp_ab, lambda4, eta2
  integer                       :: a,b,c,d,ab,cd,m
  integer                       :: a0, aa

  double precision, allocatable :: Om_tmp(:)
  double precision, allocatable :: tmp_m(:,:,:)
  double precision, allocatable :: tmp(:,:,:,:)

! Output variables

  double precision,intent(out)  :: KC(nVV,nVV)

!---------------!
! SpinOrb block !
!---------------!

  ab = 0
  do a=nO+1,nBas-nR
    do b=a+1,nBas-nR
      ab = ab + 1
      cd = 0
      do c=nO+1,nBas-nR
        do d=c+1,nBas-nR
          cd = cd + 1

          chi = 0d0
          do m=1,nS
            eps = Om(m)**2 + eta**2
            chi = chi - rho(a,c,m)*rho(b,d,m)*Om(m)/eps &
                      + rho(a,d,m)*rho(b,c,m)*Om(m)/eps
          end do
         
          KC(ab,cd) = 2d0*lambda*chi

        end do
      end do
    end do
  end do

end subroutine 
