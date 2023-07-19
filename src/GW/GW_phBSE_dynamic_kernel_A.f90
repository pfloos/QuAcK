subroutine GW_phBSE_dynamic_kernel_A(eta,nBas,nC,nO,nV,nR,nS,lambda,eGW,OmRPA,rho_RPA,OmBSE,A_dyn,ZA_dyn)

! Compute the dynamic part of the Bethe-Salpeter equation matrices

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: rho_RPA(nBas,nBas,nS)
  double precision,intent(in)   :: OmBSE
  
! Local variables

  integer                       :: maxS
  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: A_dyn(nS,nS)
  double precision,intent(out)  :: ZA_dyn(nS,nS)

! Initialization

   A_dyn(:,:) = 0d0
  ZA_dyn(:,:) = 0d0

! Number of poles taken into account 

  maxS = nS

! Build dynamic A matrix

  jb = 0
!$omp parallel do default(private) shared(A_dyn,ZA_dyn,OmRPA,OmBSE,eGW,rho_RPA,nO,nBas,maxS,chi,eps,eta,nC,nR,lambda)
  do j=nC+1,nO
    do b=nO+1,nBas-nR
      jb = (b-nO) + (j-1)*(nBas-nO) 

      ia = 0
      do i=nC+1,nO
        do a=nO+1,nBas-nR
          ia = (a-nO) + (i-1)*(nBas-nO) 
 
          chi = 0d0

          do kc=1,maxS
            eps = OmRPA(kc)
            chi = chi + rho_RPA(i,j,kc)*rho_RPA(a,b,kc)*eps/(eps**2 + eta**2)
          enddo

          A_dyn(ia,jb) = A_dyn(ia,jb) - 4d0*lambda*chi

          chi = 0d0

          do kc=1,maxS

            eps = + OmBSE - OmRPA(kc) - (eGW(a) - eGW(j))
            chi = chi + rho_RPA(i,j,kc)*rho_RPA(a,b,kc)*eps/(eps**2 + eta**2)

            eps = + OmBSE - OmRPA(kc) - (eGW(b) - eGW(i))
            chi = chi + rho_RPA(i,j,kc)*rho_RPA(a,b,kc)*eps/(eps**2 + eta**2)

          enddo

          A_dyn(ia,jb) = A_dyn(ia,jb) - 2d0*lambda*chi

          chi = 0d0
          do kc=1,maxS

            eps = + OmBSE - OmRPA(kc) - (eGW(a) - eGW(j))
            chi = chi + rho_RPA(i,j,kc)*rho_RPA(a,b,kc)*(eps**2 - eta**2)/(eps**2 + eta**2)**2

            eps = + OmBSE - OmRPA(kc) - (eGW(b) - eGW(i))
            chi = chi + rho_RPA(i,j,kc)*rho_RPA(a,b,kc)*(eps**2 - eta**2)/(eps**2 + eta**2)**2

          enddo

          ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + 2d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

!$omp end parallel do

end subroutine 
