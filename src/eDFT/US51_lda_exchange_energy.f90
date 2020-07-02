subroutine US51_lda_exchange_energy(nGrid,weight,rho,Ex,wEns,nEns)

! Compute Slater's LDA exchange energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)

! Local variables

  integer                       :: iG
  double precision              :: alpha,r,alphaw,a2,b2,c2,a1,b1,c1

! Output variables

  double precision              :: Ex

! Cxw2 parameters for He N->N+1
!  a2 = 0.135068d0
!  b2 = -0.00774769d0
!  c2 = -0.0278205d0

! Cxw1 parameters for He N->N-1
!  a1 = 0.420243d0
!  b1 = 0.0700561d0
!  c1 = -0.288301d0

! Cx coefficient for Slater LDA exchange

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)

!  alphaw = alpha*(1d0 - wEns(2)*(1d0 - wEns(2))*(a1 + b1*(wEns(2) - 0.5d0) + c1*(wEns(2) - 0.5d0)**2))
! Compute LDA exchange energy

  Ex = 0d0
  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      Ex = Ex  + weight(iG)*alpha*r**(4d0/3d0)

    endif

  enddo

end subroutine US51_lda_exchange_energy
