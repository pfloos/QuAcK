subroutine RMFL20_lda_exchange_energy(nEns,wEns,nGrid,weight,rho,Ex)

! Compute the restricted version of the Marut-Fromager-Loos weight-dependent exchange functional
! The RMFL20 is a two-state, single-weight exchange functional

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
  double precision              :: Cx0
  double precision              :: Cx1
  double precision              :: CxLDA
  double precision              :: Cxw
  double precision              :: r

! Output variables

  double precision              :: Ex

! Cx coefficient for Slater LDA exchange

  Cx0   = -(4d0/3d0)*(1d0/pi)**(1d0/3d0)
  Cx1   = -(176d0/105d0)*(1d0/pi)**(1d0/3d0)
  CxLDA = -(3d0/4d0)*(3d0/pi)**(1d0/3d0)

  Cxw   = CxLDA + wEns(1)*(Cx1 - Cx0)

! Compute LDA exchange energy

  Ex = 0d0
  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      Ex = Ex + weight(iG)*Cxw*r**(4d0/3d0)

    endif

  enddo

end subroutine RMFL20_lda_exchange_energy
