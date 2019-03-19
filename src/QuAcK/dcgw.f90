function SigC_dcgw(x,y) result(SigC)

! Degeneracy-corrected GW

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: x,y

! Local variables

  double precision,parameter    :: eta = 0.1d0
  double precision              :: r

! Output variables

  double precision              :: SigC

! Compute the divergence-free term

  r = y/x

! Bare style

  SigC = y*r

! DCPT style

! SigC = -0.5d0*x*(1d0-sqrt(1d0+4d0*r*r))

! QDPT style

! SigC = y*r/sqrt(1d0+4d0*r*r)

! Infinitesimal

! SigC = y*y*x/(x*x+eta*eta)

end function SigC_dcgw

function Z_dcgw(x,y) result(Z)


! Derivative of the degeneracy-corrected GW

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: x,y

! Local variables

  double precision,parameter    :: eta = 1d0
  double precision              :: r

! Output variables

  double precision              :: Z

! Compute the divergence-free term

  r = y/x 

! Bare style

  Z = r*r

! DCPT style

! Z = 0.5d0*(1d0-1d0/sqrt(1d0+4d0*r*r))

! QDPT style

! Z = r/sqrt(1d0+4d0*r*r)/(1d0+4d0*r*r)

! Infinitesimal

! Z = y*y*(x*x-eta*eta)/(x*x+eta*eta)**2

end function Z_dcgw
