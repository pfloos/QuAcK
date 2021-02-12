subroutine ULYP_gga_correlation_energy(nGrid,weight,rho,drho,Ec)

! Compute unrestricted LYP GGA correlation energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin)

! Local variables

  integer                       :: iG
  double precision              :: ra,rb,r
  double precision              :: ga,gab,gb,g

  double precision              :: a,b,c,d
  double precision              :: Cf,omega,delta

! Output variables

  double precision              :: Ec(nsp)

! Parameters of the functional

  a = 0.04918d0
  b = 0.132d0
  c = 0.2533d0
  d = 0.349d0

  Cf = 3d0/10d0*(3d0*pi**2)**(2d0/3d0)

! Initialization

  Ec(:) = 0d0

  do iG=1,nGrid

    ra = max(0d0,rho(iG,1))
    rb = max(0d0,rho(iG,2))
    r  = ra + rb

    if(r > threshold) then

      ga  = drho(1,iG,1)**2 + drho(2,iG,1)**2 + drho(3,iG,1)**2
      gb  = drho(1,iG,2)**2 + drho(2,iG,2)**2 + drho(3,iG,2)**2
      gab = drho(1,iG,1)*drho(1,iG,2) + drho(2,iG,1)*drho(2,iG,2) + drho(3,iG,1)*drho(3,iG,2)
      g   = ga + gab + gb

      omega = exp(-c*r**(-1d0/3d0))/(1d0 + d*r**(-1d0/3d0))*r**(-11d0/3d0)
      delta = c*r**(-1d0/3d0) + d*r**(-1d0/3d0)/(1d0 + d*r**(-1d0/3d0))

      Ec(2) = Ec(2) - weight(iG)*4d0*a/(1d0 + d*r**(-1d0/3d0))*ra*rb/r     &
                    - weight(iG)*a*b*omega*ra*rb*(                         &
                        2d0**(11d0/3d0)*Cf*(ra**(8d0/3d0) + rb**(8d0/3d0)) &
                      + (47d0/18d0 - 7d0*delta/18d0)*g                     &
                      - (5d0/2d0 - delta/18d0)*(ga + gb)                   & 
                      - (delta - 11d0)/9d0*(ra/r*ga + rb/r*gb) )           &
                    - weight(iG)*a*b*omega*(                               &
                        - 2d0*r**2/3d0*g                                   &
                        + (2d0*r**2/3d0 - ra**2)*gb                        &
                        + (2d0*r**2/3d0 - rb**2)*ga )            

    end if

  end do

end subroutine ULYP_gga_correlation_energy
