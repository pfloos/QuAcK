subroutine dynamic_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOO,nVV,lambda,eGT,Omega1,Omega2,rho1,rho2,OmBSE,TA,ZA)

! Compute the dynamic part of the Bethe-Salpeter equation matrices for GT

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
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: OmBSE
  

  double precision,intent(in)   :: Omega1(nVV)
  double precision,intent(in)   :: Omega2(nOO)
  double precision,intent(in)   :: rho1(nBas,nBas,nVV)
  double precision,intent(in)   :: rho2(nBas,nBas,nOO)

! Local variables

  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,cd,kl

! Output variables

  double precision,intent(out)  :: TA(nS,nS)
  double precision,intent(out)  :: ZA(nS,nS)

! Initialization

  TA(:,:) = 0d0
  ZA(:,:) = 0d0

! Build dynamic A matrix

  jb = 0
!$omp parallel do default(private) shared(TA,ZA,Omega1,Omega2,OmBSE,eGT,rho1,rho2,nO,nBas,nVV,nOO,chi,eps,eta,nC,nR,lambda)
  do j=nC+1,nO
    do b=nO+1,nBas-nR
      jb = (b-nO) + (j-1)*(nBas-nO) 

      ia = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          ia = (a-nO) + (i-1)*(nBas-nO) 
 
          chi = 0d0

          do cd=1,nVV
            eps = + OmBSE - Omega1(cd) + (eGT(i) + eGT(j))
            chi = chi + rho1(i,b,cd)*rho1(a,j,cd)*eps/(eps**2 + eta**2) 
          end do

          do kl=1,nOO
            eps = + OmBSE + Omega2(kl) - (eGT(a) + eGT(b))
            chi = chi + rho2(i,b,kl)*rho2(a,j,kl)*eps/(eps**2 + eta**2)
          end do

          TA(ia,jb) = TA(ia,jb) - lambda*chi

          chi = 0d0

          do cd=1,nVV
            eps = + OmBSE - Omega1(cd) + (eGT(i) + eGT(j))
            chi = chi + rho1(i,b,cd)*rho1(a,j,cd)*(eps**2 - eta**2)/(eps**2 + eta**2)**2
          end do

          do kl=1,nOO
            eps = + OmBSE + Omega2(kl) - (eGT(a) + eGT(b))
            chi = chi + rho2(i,b,kl)*rho2(a,j,kl)*(eps**2 - eta**2)/(eps**2 + eta**2)**2
          end do

          ZA(ia,jb) = ZA(ia,jb) + lambda*chi

        end do
      end do
    end do
  end do

!$omp end parallel do

end subroutine dynamic_Tmatrix_A
