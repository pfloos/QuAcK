subroutine PBE_gga_correlation_energy(nGrid,weight,rho,drho,Ec)

! Compute unrestricted PBE GGA correlation energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin)

! Local variables

  integer                       :: iG
  double precision              :: ra,rb,r,rs,z
  double precision              :: ga,gab,gb,g

  double precision              :: a,b,c,d
  double precision              :: gam,beta

  double precision              :: A_p,a1_p,b1_p,b2_p,b3_p,b4_p
  double precision              :: A_f,a1_f,b1_f,b2_f,b3_f,b4_f
  double precision              :: A_a,a1_a,b1_a,b2_a,b3_a,b4_a

  double precision              :: ec_z,ec_p,ec_f,ec_a
  double precision              :: fz,d2fz

  double precision              :: H,kf,ks,t,phi

! Output variables

  double precision              :: Ec(nsp)

! Parameters for PW92

  A_p  = 0.031091d0
  a1_p = 0.21370d0
  b1_p = 7.5957d0
  b2_p = 3.5876d0
  b3_p = 1.6382d0
  b4_p = 0.49294d0

  A_f  = 0.015545d0
  a1_f = 0.20548d0
  b1_f = 14.1189d0
  b2_f = 6.1977d0
  b3_f = 3.3662d0
  b4_f = 0.62517d0

  A_a  = 0.016887d0
  a1_a = 0.11125d0
  b1_a = 10.357d0
  b2_a = 3.6231d0
  b3_a = 0.88026d0
  b4_a = 0.49671d0

! Parameters PBE

  gam  = (1d0 - log(2d0))/pi**2
  beta = 0.066725d0

! Initialization

  Ec(:) = 0d0

  do iG=1,nGrid

    ra = max(0d0,rho(iG,1))
    rb = max(0d0,rho(iG,2))
    r  = ra + rb
    z  = (ra - rb)/r

!   alpha-alpha contribution

    if(ra > threshold) then

      rs = (4d0*pi*ra/3d0)**(-1d0/3d0)

      ec_f = b1_f*sqrt(rs) + b2_f*rs + b3_f*rs**(3d0/2d0) + b4_f*rs**2
      ec_f = -2d0*A_f*(1d0 + a1_f*rs)*log(1d0 + 1d0/(2d0*A_f*ec_f))

      ga  = drho(1,iG,1)*drho(1,iG,1) + drho(2,iG,1)*drho(2,iG,1) + drho(3,iG,1)*drho(3,iG,1)

      kf = (3d0*pi**2*ra)**(1d0/3d0)
      ks = sqrt(4d0*kf/pi)
      phi = 1d0
      t = sqrt(ga)/(2d0*phi*ks*ra)

      A = beta/gam/(exp(-ec_f/(gam*phi**3)) - 1d0)

      H = gam*phi**3*log(1d0 + beta/gam*t**2*((1d0 + A*t**2)/(1d0 + A*t**2 + A**2*t**4)))

      Ec(1) = Ec(1) + weight(iG)*(ec_f + H)*ra

    end if

    r  = ra + rb

!   alpha-beta contribution

    if(r > threshold) then

      rs = (4d0*pi*r/3d0)**(-1d0/3d0)

      fz = (1d0 + z)**(4d0/3d0) + (1d0 - z)**(4d0/3d0) - 2d0
      fz = fz/(2d0*(2d0**(1d0/3d0) - 1d0))
      d2fz = 4d0/(9d0*(2**(1d0/3d0) - 1d0))

      ec_p = b1_p*sqrt(rs) + b2_p*rs + b3_p*rs**(3d0/2d0) + b4_p*rs**2
      ec_p = -2d0*A_p*(1d0 + a1_p*rs)*log(1d0 + 1d0/(2d0*A_p*ec_p))

      ec_f = b1_f*sqrt(rs) + b2_f*rs + b3_f*rs**(3d0/2d0) + b4_f*rs**2
      ec_f = -2d0*A_f*(1d0 + a1_f*rs)*log(1d0 + 1d0/(2d0*A_f*ec_f))

      ec_a = b1_a*sqrt(rs) + b2_a*rs + b3_a*rs**(3d0/2d0) + b4_a*rs**2
      ec_a = +2d0*A_a*(1d0 + a1_a*rs)*log(1d0 + 1d0/(2d0*A_a*ec_a))

      ec_z = ec_p + ec_a*fz/d2fz*(1d0-z**4) + (ec_f - ec_p)*fz*z**4

      ga  = drho(1,iG,1)*drho(1,iG,1) + drho(2,iG,1)*drho(2,iG,1) + drho(3,iG,1)*drho(3,iG,1)
      gb  = drho(1,iG,2)*drho(1,iG,2) + drho(2,iG,2)*drho(2,iG,2) + drho(3,iG,2)*drho(3,iG,2)
      gab = drho(1,iG,1)*drho(1,iG,2) + drho(2,iG,1)*drho(2,iG,2) + drho(3,iG,1)*drho(3,iG,2)
      g   = ga + 2d0*gab + gb
   
      rs = (4d0*pi*r/3d0)**(-1d0/3d0)
      kf = (3d0*pi**2*r)**(1d0/3d0)
      ks = sqrt(4d0*kf/pi)
      phi = ((1d0 + z)**(2d0/3d0) + (1d0 - z)**(2d0/3d0))/2d0
      t = sqrt(g)/(2d0*phi*ks*r)

      A = beta/gam/(exp(-ec_p/(gam*phi**3)) - 1d0)

      H = gam*phi**3*log(1d0 + beta/gam*t**2*((1d0 + A*t**2)/(1d0 + A*t**2 + A**2*t**4)))

      Ec(2) = Ec(2) - weight(iG)*(ec_p + H)*r

    end if

!   beta-beta contribution

    if(rb > threshold) then

      rs = (4d0*pi*rb/3d0)**(-1d0/3d0)

      ec_f = b1_f*sqrt(rs) + b2_f*rs + b3_f*rs**(3d0/2d0) + b4_f*rs**2
      ec_f = -2d0*A_f*(1d0 + a1_f*rs)*log(1d0 + 1d0/(2d0*A_f*ec_f))

      gb  = drho(1,iG,2)*drho(1,iG,2) + drho(2,iG,2)*drho(2,iG,2) + drho(3,iG,2)*drho(3,iG,2)

      kf = (3d0*pi**2*rb)**(1d0/3d0)
      ks = sqrt(4d0*kf/pi)
      phi = 1d0
      t = sqrt(gb)/(2d0*phi*ks*rb)

      A = beta/gam/(exp(-ec_f/(gam*phi**3)) - 1d0)

      H = gam*phi**3*log(1d0 + beta/gam*t**2*((1d0 + A*t**2)/(1d0 + A*t**2 + A**2*t**4)))

      Ec(3) = Ec(3) + weight(iG)*(ec_f + H)*rb

    end if

  end do

  Ec(2) = Ec(2) - Ec(1) - Ec(3)


end subroutine PBE_gga_correlation_energy
