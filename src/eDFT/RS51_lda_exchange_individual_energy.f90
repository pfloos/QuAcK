subroutine RS51_lda_exchange_individual_energy(nGrid,weight,rhow,rho,Ex)

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
  double precision              :: r,rI
  double precision              :: e,dedr

! Output variables

  double precision,intent(out)  :: Ex

! Compute LDA exchange matrix in the AO basis

  Ex = 0d0
  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))

    if(r > threshold .or. rI > threshold) then

      e    =         CxLDA*r**(1d0/3d0)
      dedr = 1d0/3d0*CxLDA*r**(-2d0/3d0)

      Ex = Ex + weight(iG)*(e*rI + dedr*r*rI - dedr*r*r)

    endif

  enddo

end subroutine RS51_lda_exchange_individual_energy
