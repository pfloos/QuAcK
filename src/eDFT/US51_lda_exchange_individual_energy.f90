subroutine US51_lda_exchange_individual_energy(nGrid,weight,rhow,Ex)

! Compute the restricted version of Slater's LDA exchange individual energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r
  double precision              :: dedr

! Output variables

  double precision,intent(out)  :: Ex

! Compute LDA exchange matrix in the AO basis

  Ex = 0d0

  do iG=1,nGrid

    r  = max(0d0,rhow(iG))
    if(r > threshold) then

      dedr = 1d0/3d0*CxLSDA*r**(-2d0/3d0)

      Ex = Ex - weight(iG)*dedr*r*r      

    endif

  enddo

end subroutine US51_lda_exchange_individual_energy
