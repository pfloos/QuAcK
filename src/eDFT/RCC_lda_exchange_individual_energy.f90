subroutine RCC_lda_exchange_individual_energy(nEns,wEns,aCC_w1,aCC_w2,nGrid,weight,rhow,rho,Ex,Cx_choice)

! Compute the restricted version of the curvature-corrected exchange functional

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
  double precision,intent(in)   :: rho(nGrid)
  integer,intent(in)            :: Cx_choice

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



  a1 = aCC_w1(1)
  b1 = aCC_w1(2)
  c1 = aCC_w1(3)
 


  a2 = aCC_w2(1)
  b2 = aCC_w2(2)
  c2 = aCC_w2(3)


  w1 = wEns(2)
  Fx1 = 1d0 - w1*(1d0 - w1)*(a1 + b1*(w1 - 0.5d0) + c1*(w1 - 0.5d0)**2)

  w2 = wEns(3)
  Fx2 = 1d0 - w2*(1d0 - w2)*(a2 + b2*(w2 - 0.5d0) + c2*(w2 - 0.5d0)**2)

  select case (Cx_choice)

    case(1)
      Cx = CxLDA*Fx1

    case(2)
      Cx = CxLDA*Fx2

    case(3)
      Cx = CxLDA*Fx2*Fx1

    case default
      Cx = CxLDA

  end select

! Compute LDA exchange matrix in the AO basis

  Ex = 0d0
  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold .or. rI > threshold) then

      e_p  =         Cx*r**(1d0/3d0)
      dedr = 1d0/3d0*Cx*r**(-2d0/3d0)
      Ex = Ex + weight(iG)*(e_p*rI + dedr*r*rI - dedr*r*r)

    endif

  enddo

end subroutine RCC_lda_exchange_individual_energy
