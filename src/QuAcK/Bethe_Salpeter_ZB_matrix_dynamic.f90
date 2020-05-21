subroutine Bethe_Salpeter_ZB_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,lambda,eGW,OmRPA,OmBSE,rho,ZB_dyn)

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

  double precision,intent(out)  :: ZB_dyn(nS,nS)

! Initialization

  ZB_dyn(:,:) = 0d0

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

            eps = (OmBSE - OmRPA(kc) - (eGW(j) - eGW(i)))**2 + eta**2
            chi = chi + rho(i,b,kc)*rho(a,j,kc)*((OmBSE - OmRPA(kc) - (eGW(j) - eGW(i)))/eps)**2

            eps = (OmBSE - OmRPA(kc) - (eGW(a) - eGW(b)))**2 + eta**2
            chi = chi + rho(i,b,kc)*rho(a,j,kc)*((OmBSE - OmRPA(kc) - (eGW(a) - eGW(b)))/eps)**2

          enddo

          ZB_dyn(ia,jb) = ZB_dyn(ia,jb) + 2d0*lambda*chi

        enddo
      enddo
    enddo
  enddo

end subroutine Bethe_Salpeter_ZB_matrix_dynamic
