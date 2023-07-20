subroutine GTpp_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,lambda,Omega1,rho1,Omega2,rho2,KA)

! Compute the OOVV block of the static T-matrix

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
  integer                       :: i,j,a,b,ia,jb,kl,cd,c,d

! Output variables

  double precision,intent(out)  :: KA(nS,nS)

  KA(:,:) = 0d0

  jb = 0
!$omp parallel do default(private) shared(KA,Omega1,Omega2,rho1,rho2,nO,nBas,nVV,nOO,chi,eps,eta,nC,nR,lambda)
  do j=nC+1,nO
    do b=nO+1,nBas-nR
      jb = (b-nO) + (j-1)*(nBas-nO) 

      ia = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          ia = (a-nO) + (i-1)*(nBas-nO) 

          chi = 0d0

          do cd=1,nVV
            eps = + Omega1(cd)
            chi = chi + rho1(i,b,cd)*rho1(a,j,cd)*eps/(eps**2 + eta**2)

          enddo

          do kl=1,nOO
            eps = - Omega2(kl)
            chi = chi + rho2(i,b,kl)*rho2(a,j,kl)*eps/(eps**2 + eta**2)
          enddo

          KA(ia,jb) = lambda*chi

        enddo
      enddo
    enddo
  enddo

!$omp end parallel do

end subroutine 
