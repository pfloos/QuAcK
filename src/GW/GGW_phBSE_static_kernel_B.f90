subroutine GGW_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,Om,rho,KB)

! Compute the BSE static kernel for the coupling block

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: KB(nS,nS)

! Compute static kernel

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1

          chi = 0d0
          do kc=1,nS
            eps = Om(kc)**2 + eta**2
            chi = chi + rho(i,b,kc)*rho(a,j,kc)*Om(kc)/eps
          end do

          KB(ia,jb) = 2d0*lambda*chi

        end do
      end do
    end do
  end do

end subroutine 
