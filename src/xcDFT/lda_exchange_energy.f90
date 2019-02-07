subroutine lda_exchange_energy(nGrid,weight,rho,Ex)

! Compute LDA exchange energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: alpha 

! Output variables

  double precision              :: Ex

! Cx coefficient for Slater LDA exchange

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)

! Compute LDA exchange energy

  Ex = 0d0
  do iG=1,nGrid
    Ex = Ex  + weight(iG)*alpha*rho(iG)**(4d0/3d0)
  enddo

  Ex = 2d0*Ex

end subroutine lda_exchange_energy
