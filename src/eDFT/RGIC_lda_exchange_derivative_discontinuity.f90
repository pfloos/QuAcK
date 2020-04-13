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

! Compute correlation energy for ground- and doubly-excited states


! Parameters for H2 at equilibrium

! a = + 0.5751782560799208d0
! b = - 0.021108186591137282d0
! c = - 0.36718902716347124d0

! Parameters for stretch H2

! a = + 0.01922622507087411d0
! b = - 0.01799647558018601d0
! c = - 0.022945430666782573d0

! Parameters for He

! a = 1.9125735895875828d0
! b = 2.715266992840757d0
! c = 2.1634223380633086d0

! Parameters for HNO

  a = 0.0061158387543040335d0
  b = -0.00005968703047293955d0
  c = -0.00001692245714408755d0

  w = wEns(2)
  dCxGICdw = (0.5d0*b + (2d0*a + 0.5d0*c)*(w - 0.5d0) - (1d0 - w)*w*(3d0*b + 4d0*c*(w - 0.5d0)))
  dCxGICdw = CxLDA*dCxGICdw

  dExdw(:) = 0d0

  do iG=1,nGrid
    
    r = max(0d0,rhow(iG))
    
    if(r > threshold) then
 
      dExdw(1) = 0d0
      dExdw(2) = dExdw(2) + weight(iG)*dCxGICdw*r**(4d0/3d0)

    end if
     
  end do 

  ExDD(:) = 0d0
 
  do iEns=1,nEns
    do jEns=2,nEns

      ExDD(iEns) = ExDD(iEns) + (Kronecker_delta(iEns,jEns) - wEns(jEns))*dExdw(jEns)

    end do
  end do

end subroutine RGIC_lda_exchange_derivative_discontinuity
