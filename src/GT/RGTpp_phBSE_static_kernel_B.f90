subroutine RGTpp_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,lambda,Omega1,rho1,Omega2,rho2,KB)

! Compute the OVVO block of the static T-matrix

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: Omega1(nVV)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: Omega2(nOO)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kl,cd

! Output variables

  double precision,intent(out)  :: KB(nS,nS)

  KB(:,:) = 0d0

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

          chi = 0d0

          do cd=1,nVV
            eps = + Omega1(cd)
            chi = chi + rho1(i,j,cd)*rho1(a,b,cd)*eps/(eps**2 + eta**2)
          end do

          do kl=1,nOO
            eps = - Omega2(kl)
            chi = chi + rho2(i,j,kl)*rho2(a,b,kl)*eps/(eps**2 + eta**2)
          end do

          KB(ia,jb) = lambda*chi

        end do
      end do
    end do
  end do

end subroutine 
