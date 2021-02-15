subroutine US51_lda_exchange_energy(nGrid,weight,rho,Ex)

  use xc_f90_lib_m

! Compute Slater's LDA exchange energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r

! Output variables

  double precision              :: Ex

! Compute LDA exchange energy

  Ex = 0d0
  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      Ex = Ex  + weight(iG)*CxLSDA*r**(4d0/3d0)

    endif

  enddo

end subroutine US51_lda_exchange_energy
