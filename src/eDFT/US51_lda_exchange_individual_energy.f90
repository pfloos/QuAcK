subroutine US51_lda_exchange_individual_energy(nGrid,weight,rhow,rho,Ex)

! Compute the restricted version of Slater's LDA exchange individual energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r,rI,alpha
  double precision              :: e,dedr

! Output variables

  double precision,intent(out)  :: Ex

! Compute LDA exchange matrix in the AO basis

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)


  Ex = 0d0
  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold .or. rI > threshold) then

      e    =         alpha*r**(1d0/3d0)
      dedr = 1d0/3d0*alpha*r**(-2d0/3d0)

      Ex = Ex + weight(iG)*(e*rI + dedr*r*rI - dedr*r*r)

    endif

  enddo

end subroutine US51_lda_exchange_individual_energy
