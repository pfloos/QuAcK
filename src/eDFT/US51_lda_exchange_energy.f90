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
  double precision              :: alpha,r,alphaw,a2,b2,c2,a1,b1,c1

! Output variables

  double precision              :: Ex

! Cx coefficient for Slater LDA exchange

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)

! Compute LDA exchange energy

  Ex = 0d0
  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      Ex = Ex  + weight(iG)*alpha*r**(4d0/3d0)

    endif

  enddo

end subroutine US51_lda_exchange_energy
