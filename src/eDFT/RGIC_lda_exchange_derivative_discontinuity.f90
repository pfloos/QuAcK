subroutine RGIC_lda_exchange_derivative_discontinuity(nEns,wEns,nGrid,weight,rhow,ExDD)

! Compute the restricted version of the GIC exchange individual energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)

! Local variables

  integer                       :: iEns,jEns
  integer                       :: iG
  double precision              :: r
  double precision,allocatable  :: dExdw(:)
  double precision,external     :: Kronecker_delta

  double precision              :: a,b,c,w
  double precision              :: dCxGICdw

! Output variables

  double precision,intent(out)  :: ExDD(nEns)

! Memory allocation

  allocate(dExdw(nEns))

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
  dCxGICdw = (0.5d0*b + (2d0*a + 0.5d0*c)*(w - 0.5d0) - (1d0 - w)*w*(3d0*b + 4d0*c*(w - 0.5d0)))
  dCxGICdw = CxLDA*dCxGICdw

  dExdw(:) = 0d0

  do iG=1,nGrid
    
    r = max(0d0,rhow(iG))
    
    if(r > threshold) then
 
      dExdw(1) = 0d0
      dExdw(2) = dExdw(2) + 0.5d0*weight(iG)*dCxGICdw*r**(4d0/3d0)
      dExdw(3) = dExdw(3) + 1.0d0*weight(iG)*dCxGICdw*r**(4d0/3d0)

    end if
     
  end do 

  ExDD(:) = 0d0
 
  do iEns=1,nEns
    do jEns=2,nEns

      ExDD(iEns) = ExDD(iEns) + (Kronecker_delta(iEns,jEns) - wEns(jEns))*dExdw(jEns)

    end do
  end do

end subroutine RGIC_lda_exchange_derivative_discontinuity
