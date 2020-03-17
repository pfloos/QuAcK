subroutine RMFL20_lda_exchange_potential(nEns,wEns,nGrid,weight,nBas,AO,rho,Fx)

! Compute the restricted version of the weight-dependent MFL20 exchange potential

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: Cx0
  double precision              :: Cx1
  double precision              :: CxLDA
  double precision              :: Cxw
  double precision              :: r,vAO

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Weight-dependent Cx coefficient for RMFL20 exchange functional

  Cx0   = -(4d0/3d0)*(1d0/pi)**(1d0/3d0)
  Cx1   = -(176d0/105d0)*(1d0/pi)**(1d0/3d0)
  CxLDA = -(3d0/4d0)*(3d0/pi)**(1d0/3d0)

  Cxw   = CxLDA + wEns(2)*(Cx1 - Cx0)

! Compute LDA exchange matrix in the AO basis

  Fx(:,:) = 0d0
  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        r = max(0d0,rho(iG))

        if(r > threshold) then

          vAO = weight(iG)*AO(mu,iG)*AO(nu,iG)
          Fx(mu,nu) = Fx(mu,nu) + vAO*4d0/3d0*Cxw*r**(1d0/3d0)

        endif

      enddo
    enddo
  enddo

end subroutine RMFL20_lda_exchange_potential
