subroutine UB88_gga_exchange_energy(nGrid,weight,rho,drho,Ex)

! Compute Becke's 88 GGA exchange energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)

! Local variables

  integer                       :: iG
  double precision              :: b
  double precision              :: r,g,x

! Output variables

  double precision              :: Ex

! Coefficients for B88 GGA exchange functional

  b = 0.0042d0

! Compute GGA exchange energy

  Ex = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then 
      g = drho(1,iG)**2 + drho(2,iG)**2 + drho(3,iG)**2
      x = sqrt(g)/r**(4d0/3d0)

      Ex = Ex + weight(iG)*r**(4d0/3d0)*(CxLSDA - b*x**2/(1d0 + 6d0*b*x*asinh(x)))

    end if

  end do

end subroutine UB88_gga_exchange_energy
