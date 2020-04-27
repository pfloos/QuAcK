subroutine RGIC_lda_exchange_potential(nEns,wEns,nGrid,weight,nBas,AO,rho,Fx)

! Compute the restricted version of the GIC exchange potential

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
  double precision              :: r,vAO

  double precision              :: a,b,c,w
  double precision              :: CxGIC

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Weight-dependent Cx coefficient for RMFL20 exchange functional

! Parameters for H2 at equilibrium

  a = + 0.5739189000851961d0
  b = - 0.0003469882157336496d0
  c = - 0.2676338054343272d0

! Parameters for stretch H2

! a = + 0.01918229168254928d0
! b = - 0.01545313842512261d0
! c = - 0.012720073519142448d0

! Parameters for He

! a = 1.9015719148496788d0
! b = 2.5236598782764412d0
! c = 1.6652282199359842d0

  w = 0.5d0*wEns(2) + wEns(3)
  CxGIC = 1d0 - w*(1d0 - w)*(a + b*(w - 0.5d0) + c*(w - 0.5d0)**2)
  CxGIC = CxLDA*CxGIC

! Compute LDA exchange matrix in the AO basis

  Fx(:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        r = max(0d0,rho(iG))

        if(r > threshold) then

          vAO = weight(iG)*AO(mu,iG)*AO(nu,iG)
          Fx(mu,nu) = Fx(mu,nu) + vAO*4d0/3d0*CxGIC*r**(1d0/3d0)

        endif

      enddo
    enddo
  enddo

end subroutine RGIC_lda_exchange_potential
