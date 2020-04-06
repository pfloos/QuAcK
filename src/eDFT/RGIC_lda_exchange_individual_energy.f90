subroutine RGIC_lda_exchange_individual_energy(nEns,wEns,nGrid,weight,rhow,rho,Ex)

! Compute the restricted version of the GIC exchange functional

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
  double precision              :: CxGIC
  double precision              :: r,rI
  double precision              :: e_p,dedr

  double precision              :: a,b,c,w

! Output variables

  double precision,intent(out)  :: Ex

! Weight-dependent Cx coefficient for RMFL20 exchange functional

  a = + 0.5751782560799208d0
  b = - 0.021108186591137282d0
  c = - 0.36718902716347124d0

  w = wEns(2)
  CxGIC = CxLDA*w*(1d0 - w)*(a + b*(w - 0.5d0) + c*(w - 0.5d0)**2)
 
! Compute LDA exchange matrix in the AO basis

  Ex = 0d0
  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold .and. rI > threshold) then

      e_p  =         CxLDA*r**(1d0/3d0)
      dedr = 1d0/3d0*CxLDA*r**(-2d0/3d0)
      Ex = Ex + weight(iG)*(e_p*rI + dedr*r*rI - dedr*r*r)

      e_p  =         CxGIC*r**(1d0/3d0)
      dedr = 1d0/3d0*CxGIC*r**(-2d0/3d0)
      Ex = Ex + weight(iG)*(e_p*rI + dedr*r*rI - dedr*r*r)

    endif

  enddo

end subroutine RGIC_lda_exchange_individual_energy
