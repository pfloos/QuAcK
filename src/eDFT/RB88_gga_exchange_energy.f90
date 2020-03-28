subroutine RB88_gga_exchange_energy(nGrid,weight,rho,drho,Ex)

! Compute restricted version of Becke's 88 GGA exchange energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Local variables

  integer                       :: iG
  double precision              :: alpha
  double precision              :: beta
  double precision              :: r,g,x

! Output variables

  double precision              :: Ex

! Coefficients for B88 GGA exchange functional

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)
  beta  = 0.0042d0

! Compute GGA exchange energy

  Ex = 0d0

  do iG=1,nGrid

    r = max(0d0,0.5d0*rho(iG))

    if(r > threshold) then 
      g = drho(1,iG)**2 + drho(2,iG)**2 + drho(3,iG)**2
      x = sqrt(g)/r**(4d0/3d0)

      Ex = Ex + weight(iG)*alpha*r**(4d0/3d0) & 
              - weight(iG)*beta*x**2*r**(4d0/3d0)/(1d0 + 6d0*beta*x*asinh(x))

    end if

  end do

  Ex = 2d0*Ex

end subroutine RB88_gga_exchange_energy
