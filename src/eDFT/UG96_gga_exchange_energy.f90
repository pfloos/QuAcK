subroutine UG96_gga_exchange_energy(nGrid,weight,rho,drho,Ex)

! Compute Gill's 96 GGA exchange energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Local variables

  integer                       :: iG
  double precision              :: alpha,beta
  double precision              :: r,g

! Output variables

  double precision              :: Ex

! Coefficients for G96 GGA exchange functional

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)
  beta  = 1d0/137d0

! Compute GGA exchange energy

  Ex = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then 
 
      g = drho(1,iG)**2 + drho(2,iG)**2 + drho(3,iG)**2

      Ex = Ex + weight(iG)*r**(4d0/3d0)*(alpha - beta*g**(3d0/4d0)/r**2)

    end if

  end do

end subroutine UG96_gga_exchange_energy
