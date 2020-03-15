subroutine RMFL20_lda_exchange_energy(nEns,wEns,nGrid,weight,rho,Ex)

! Compute restricted version of Marut-Fromager-Loos weight-dependent LDA exchange energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: CxLDA
  double precision              :: Cx2
  double precision              :: Cxw
  double precision              :: r

! Output variables

  double precision              :: Ex

! Cx coefficient for Slater LDA exchange

  CxLDA = -(4d0/3d0)*(1d0/pi)**(1d0/3d0)
  Cx2   = -(176d0/105d0)*(1d0/pi)**(1d0/3d0)

  Cxw   = (1d0 - wEns(2))*CxLDA + wEns(2)*(Cx2 - CxLDA)

! Compute LDA exchange energy

  Ex = 0d0
  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      Ex = Ex + weight(iG)*Cxw*r**(4d0/3d0)

    endif

  enddo

end subroutine RMFL20_lda_exchange_energy
