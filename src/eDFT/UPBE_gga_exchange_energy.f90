subroutine UPBE_gga_exchange_energy(nGrid,weight,rho,drho,Ex)

! Compute PBE GGA exchange energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)

! Local variables

  integer                       :: iG
  double precision              :: alpha,mupbe,kappa
  double precision              :: r,g,s2

! Output variables

  double precision              :: Ex

! Coefficients for PBE exchange functional

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)
  mupbe = ((1d0/2d0**(1d0/3d0))/(2d0*(3d0*pi**2)**(1d0/3d0)))**2*0.21951d0
  kappa = 0.804d0

! Compute GGA exchange energy

  Ex = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then 
      g = drho(1,iG)**2 + drho(2,iG)**2 + drho(3,iG)**2
      s2 = g/r**(8d0/3d0)

      Ex = Ex + weight(iG)*alpha*r**(4d0/3d0)*(1d0 + kappa - kappa/(1d0 + mupbe*s2/kappa))

    end if

  end do

end subroutine UPBE_gga_exchange_energy
