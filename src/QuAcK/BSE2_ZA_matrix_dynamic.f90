subroutine BSE2_ZA_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,lambda,eGW,OmRPA,OmBSE,rho,ZA_dyn)

! Compute the dynamic part of the Bethe-Salpeter equation matrices

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: OmBSE
  double precision,intent(in)   :: rho(nBas,nBas,nS)
  
! Local variables

  integer                       :: maxS
  double precision              :: chi
  double precision              :: eps
  integer                       :: i,j,a,b,ia,jb,kc

! Output variables

  double precision,intent(out)  :: ZA_dyn(nS,nS)

! Initialization

  ZA_dyn(:,:) = 0d0

! Number of poles taken into account 

  maxS = nS

! Build dynamic A matrix

  ia = 0
  do i=nC+1,nO
    do a=nO+1,nBas-nR
      ia = ia + 1
      jb = 0
      do j=nC+1,nO
        do b=nO+1,nBas-nR
          jb = jb + 1
 
          chi = 0d0
          do kc=1,maxS

            eps = (OmBSE - OmRPA(kc) - (eGW(a) - eGW(j)))**2 + eta**2
            chi = chi + rho(i,j,kc)*rho(a,b,kc)*((OmBSE - OmRPA(kc) - (eGW(a) - eGW(j)))/eps)**2

            eps = (OmBSE - OmRPA(kc) - (eGW(b) - eGW(i)))**2 + eta**2
            chi = chi + rho(i,j,kc)*rho(a,b,kc)*((OmBSE - OmRPA(kc) - (eGW(b) - eGW(i)))/eps)**2

          enddo

          ZA_dyn(ia,jb) = ZA_dyn(ia,jb) + 2d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

end subroutine BSE2_ZA_matrix_dynamic
