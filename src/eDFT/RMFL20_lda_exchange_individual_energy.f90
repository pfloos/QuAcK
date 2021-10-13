subroutine RMFL20_lda_exchange_individual_energy(LDA_centered,nEns,wEns,nGrid,weight,rhow,rho,Ex)

! Compute the restricted version of the Marut-Fromager-Loos 2020 weight-dependent exchange functional

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: LDA_centered
  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: Cxw
  double precision              :: r,rI
  double precision              :: e_p,dedr

  double precision,parameter    :: Cx0   = - (4d0/3d0)*(1d0/pi)**(1d0/3d0)
  double precision,parameter    :: Cx1   = - (176d0/105d0)*(1d0/pi)**(1d0/3d0)

! Output variables

  double precision,intent(out)  :: Ex

! Weight-dependent Cx coefficient for RMFL20 exchange functional

  if(LDA_centered) then
    Cxw = CxLDA + (Cx1 - Cx0)*wEns(2)
  else
    Cxw = wEns(1)*Cx0 + wEns(2)*Cx1
  end if

! Compute LDA exchange matrix in the AO basis

  Ex = 0d0
  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold .or. rI > threshold) then

      e_p  =         Cxw*r**(1d0/3d0)
      dedr = 1d0/3d0*Cxw*r**(-2d0/3d0)
      Ex = Ex + weight(iG)*(e_p*rI + dedr*r*rI - dedr*r*r)

    endif

  enddo

end subroutine RMFL20_lda_exchange_individual_energy
