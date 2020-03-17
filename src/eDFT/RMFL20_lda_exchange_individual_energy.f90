subroutine RMFL20_lda_exchange_individual_energy(nEns,wEns,nGrid,weight,rhow,rho,Ex)

! Compute the restricted version of the Marut-Fromager-Loos 2020 weight-dependent exchange functional

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: Cx0
  double precision              :: Cx1
  double precision              :: CxLDA
  double precision              :: Cxw
  double precision              :: r,rI
  double precision              :: e,dedr

! Output variables

  double precision,intent(out)  :: Ex

! Weight-dependent Cx coefficient for RMFL20 exchange functional

  Cx0   = -(4d0/3d0)*(1d0/pi)**(1d0/3d0)
  Cx1   = -(176d0/105d0)*(1d0/pi)**(1d0/3d0)
  CxLDA = -(3d0/4d0)*(3d0/pi)**(1d0/3d0)

  Cxw   = CxLDA + wEns(2)*(Cx1 - Cx0)

! Compute LDA exchange matrix in the AO basis

  Ex = 0d0
  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold .and. rI > threshold) then

      e    =         Cxw*r**(1d0/3d0)
      dedr = 1d0/3d0*Cxw*r**(-2d0/3d0)
      Ex = Ex + weight(iG)*(e*rI + dedr*r*rI - dedr*r*r)

    endif

  enddo

end subroutine RMFL20_lda_exchange_individual_energy
