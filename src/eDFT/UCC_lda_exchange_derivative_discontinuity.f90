subroutine UCC_lda_exchange_derivative_discontinuity(nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rhow,ExDD)

! Compute the unrestricted version of the curvature-corrected exchange ensemble derivative

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  double precision,intent(in)   :: aCC_w1(3)
  double precision,intent(in)   :: aCC_w2(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)

! Local variables

  integer                       :: iEns,jEns
  integer                       :: iG
  double precision              :: r,alpha
  double precision,allocatable  :: dExdw(:)
  double precision,external     :: Kronecker_delta

  double precision              :: a1,b1,c1,w1
  double precision              :: a2,b2,c2,w2
  double precision              :: dCxdw1,dCxdw2

! Output variables

  double precision,intent(out)  :: ExDD(nEns)

! Memory allocation

  allocate(dExdw(nEns))

! Single excitation parameters

!  a1 = 0.0d0
!  b1 = 0.0d0
!  c1 = 0.0d0

! Parameters for H2 at equilibrium

! a2 = +0.5751782560799208d0
! b2 = -0.021108186591137282d0
! c2 = -0.36718902716347124d0

! Parameters for stretch H2

!  a2 = + 0.01922622507087411d0
!  b2 = - 0.01799647558018601d0
!  c2 = - 0.022945430666782573d0

! Parameters for He

! a2 = 1.9125735895875828d0
! b2 = 2.715266992840757d0
! c2 = 2.1634223380633086d0

! Parameters for He N -> N-1

  a1 = aCC_w1(1)
  b1 = aCC_w1(2)
  c1 = aCC_w1(3)


! Parameters for He N -> N+1

  a2 = aCC_w2(1)
  b2 = aCC_w2(2)
  c2 = aCC_w2(3)
 

! Cx coefficient for unrestricted Slater LDA exchange

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)

  w1 = wEns(2)
  w2 = wEns(3)

! Double weight-dependency

!  dCxdw1 = (0.5d0*b1 + (2d0*a1 + 0.5d0*c1)*(w1 - 0.5d0) - (1d0 - w1)*w1*(3d0*b1 + 4d0*c1*(w1 - 0.5d0))) &
!         * (1d0 - w2*(1d0 - w2)*(a2 + b2*(w2 - 0.5d0) + c2*(w2 - 0.5d0)**2))

!  dCxdw2 = (1d0 - w1*(1d0 - w1)*(a1 + b1*(w1 - 0.5d0) + c1*(w1 - 0.5d0)**2))                            &
!         * (0.5d0*b2 + (2d0*a2 + 0.5d0*c2)*(w2 - 0.5d0) - (1d0 - w2)*w2*(3d0*b2 + 4d0*c2*(w2 - 0.5d0)))  

! left single-weight-dependency
!  dCxdw1 = (0.5d0*b1 + (2d0*a1 + 0.5d0*c1)*(w1 - 0.5d0) - (1d0 - w1)*w1*(3d0*b1 + 4d0*c1*(w1 - 0.5d0)))
!  dCxdw2 = 0.d0

! right single-weight-dependency
  dCxdw1 = 0.d0
  dCxdw2 =(0.5d0*b2 + (2d0*a2 + 0.5d0*c2)*(w2 - 0.5d0) - (1d0 - w2)*w2*(3d0*b2 + 4d0*c2*(w2 - 0.5d0)))



  dCxdw1 = alpha*dCxdw1
  dCxdw2 = alpha*dCxdw2

  dExdw(:) = 0d0

  do iG=1,nGrid
    
    r = max(0d0,rhow(iG))
    
    if(r > threshold) then
 
      dExdw(1) = 0d0
      dExdw(2) = dExdw(2) + weight(iG)*dCxdw1*r**(4d0/3d0)
      dExdw(3) = dExdw(3) + weight(iG)*dCxdw2*r**(4d0/3d0)

    end if
     
  end do 

  ExDD(:) = 0d0
 
  do iEns=1,nEns
    do jEns=2,nEns

      ExDD(iEns) = ExDD(iEns) + (Kronecker_delta(iEns,jEns) - wEns(jEns))*dExdw(jEns)

    end do
  end do

end subroutine UCC_lda_exchange_derivative_discontinuity
