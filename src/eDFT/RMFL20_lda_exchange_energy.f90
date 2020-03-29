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

  logical                       :: LDA_centered = .true.
  integer                       :: iG
  double precision              :: Cxw
  double precision              :: r

! Output variables

  double precision              :: Ex

! Weight-denepdent Cx coefficient

  if(LDA_centered) then 
    Cxw = CxLDA + (Cx1 - Cx0)*wEns(2)
  else 
    Cxw = wEns(1)*Cx0 + wEns(2)*Cx1
  end if

! Cxw = CxLDA + (Cx1 - Cx0)*wEns(2)*(cos(2d0*pi*wEns(2)) + 1d0)

! Compute LDA exchange energy

  Ex = 0d0
  do iG=1,nGrid
    r = max(0d0,rho(iG))
    if(r > threshold) then
      Ex = Ex + weight(iG)*Cxw*r**(4d0/3d0)
    endif
  enddo

end subroutine RMFL20_lda_exchange_energy
