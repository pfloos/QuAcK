subroutine RB88_gga_exchange_individual_energy(nGrid,weight,rhow,drhow,rho,drho,Ex)

! Compute restricted Becke's GGA indivudal energy

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: drhow(ncart,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Local variables

  integer                       :: iG
  double precision              :: alpha
  double precision              :: beta
  double precision              :: r,rI,g,x
  double precision              :: ex_p,dexdr_p

! Output variables

  double precision,intent(out)  :: Ex

! Coefficients for B88 GGA exchange functional

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)
  beta  = 0.0042d0

! Compute GGA exchange matrix in the AO basis

  Ex = 0d0

  do iG=1,nGrid

    r  = max(0d0,0.5d0*rhow(iG))
    rI = max(0d0,0.5d0*rho(iG))

    if(r > threshold .and. rI > threshold) then

      g = 0.25d0*(drho(1,iG)**2 + drho(2,iG)**2 + drho(3,iG)**2)
      x = sqrt(g)/r**(4d0/3d0)

      dexdr_p = 4d0/3d0*r**(1d0/3d0)*(alpha - beta*g**(3d0/4d0)/r**2) &
              + 2d0*beta*g**(3d0/4d0)/r**(5d0/3d0) &
              - 2d0*3d0/4d0*beta*g**(-1d0/4d0)/r**(2d0/3d0)

      ex_p = alpha*r**(4d0/3d0) &
           - weight(iG)*beta*x**2*r**(4d0/3d0)/(1d0 + 6d0*beta*x*asinh(x))

      Ex = Ex + weight(iG)*(ex_p*rI + dexdr_p*r*rI - dexdr_p*r*r)

    end if

  end do

end subroutine RB88_gga_exchange_individual_energy
