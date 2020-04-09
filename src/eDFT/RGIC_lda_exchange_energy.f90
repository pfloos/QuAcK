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

! Weight-dependent Cx coefficient

! Parameters for H2 at equilibrium

! a = + 0.5751782560799208d0
! b = - 0.021108186591137282d0
! c = - 0.36718902716347124d0

! Parameters for stretch H2

! a = + 0.01922622507087411d0
! b = - 0.01799647558018601d0
! c = - 0.022945430666782573d0

! Parameters for He

  a = 1.9125735895875828d0
  b = 2.715266992840757d0
  c = 2.1634223380633086d0

  w = wEns(2)
  CxGIC = 1d0 - w*(1d0 - w)*(a + b*(w - 0.5d0) + c*(w - 0.5d0)**2)
  CxGIC = CxLDA*CxGIC

! Compute GIC-LDA exchange energy

  Ex = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then
      Ex = Ex + weight(iG)*CxGIC*r**(4d0/3d0) 
    endif

  enddo

end subroutine RGIC_lda_exchange_energy
