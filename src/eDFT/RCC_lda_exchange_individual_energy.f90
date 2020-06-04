subroutine RCC_lda_exchange_individual_energy(nEns,wEns,nGrid,weight,rhow,rho,Ex)

! Compute the restricted version of the curvature-corrected exchange functional

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
  double precision              :: r,rI
  double precision              :: e_p,dedr

  double precision              :: a1,b1,c1,w1
  double precision              :: a2,b2,c2,w2
  double precision              :: Fx1,Fx2,Cx

! Output variables

  double precision,intent(out)  :: Ex

! Single excitation parameter

  a1 = 0.0d0
  b1 = 0.0d0
  c1 = 0.0d0

! Parameters for H2 at equilibrium

! a2 = +0.5751782560799208d0
! b2 = -0.021108186591137282d0
! c2 = -0.36718902716347124d0

! Parameters for stretch H2

  a2 = + 0.01922622507087411d0
  b2 = - 0.01799647558018601d0
  c2 = - 0.022945430666782573d0

! Parameters for He

! a2 = 1.9125735895875828d0
! b2 = 2.715266992840757d0
! c2 = 2.1634223380633086d0

  w1 = wEns(2)
  Fx1 = 1d0 - w1*(1d0 - w1)*(a1 + b1*(w1 - 0.5d0) + c1*(w1 - 0.5d0)**2)

  w2 = wEns(3)
  Fx2 = 1d0 - w2*(1d0 - w2)*(a2 + b2*(w2 - 0.5d0) + c2*(w2 - 0.5d0)**2)

  Cx = CxLDA*Fx1*Fx2

! Compute LDA exchange matrix in the AO basis

  Ex = 0d0
  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold .and. rI > threshold) then

      e_p  =         Cx*r**(1d0/3d0)
      dedr = 1d0/3d0*Cx*r**(-2d0/3d0)
      Ex = Ex + weight(iG)*(e_p*rI + dedr*r*rI - dedr*r*r)

    endif

  enddo

end subroutine RCC_lda_exchange_individual_energy
