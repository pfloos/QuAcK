subroutine RGIC_lda_exchange_energy(nEns,wEns,nGrid,weight,rho,Ex)

! Compute the restricted version of the GIC exchange functional

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
  double precision              :: r

  double precision              :: a,b,c,w
  double precision              :: CxGIC

! Output variables

  double precision              :: Ex

! Weight-denepdent Cx coefficient

  a = + 0.5751782560799208d0
  b = - 0.021108186591137282d0
  c = - 0.36718902716347124d0

  w = wEns(2)
  CxGIC = CxLDA*w*(1d0 - w)*(a + b*(w - 0.5d0) + c*(w - 0.5d0)**2)

! Compute GIC-LDA exchange energy

  Ex = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then
      Ex = Ex + weight(iG)*CxLDA*r**(4d0/3d0) 
      Ex = Ex + weight(iG)*CxGIC*r**(4d0/3d0) 
    endif

  enddo

end subroutine RGIC_lda_exchange_energy
