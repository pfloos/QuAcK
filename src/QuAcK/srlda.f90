subroutine srlda(rs,mu,ecsr)

! Correlation energy of a spin unpolarized uniform electron gas 
! with short-range interaction erfc(mu*r)/r
! See Zecca et al. PRB 70, 205127 (2004)

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)  :: rs
  double precision,intent(in)  :: mu

! Local variables

  double precision             :: ec
  double precision             :: cf
  double precision             :: b1
  double precision             :: b2
  double precision             :: b3
  double precision             :: b4
  double precision             :: a0
  double precision             :: bb
  double precision             :: m1

! Ouput variables

  double precision,intent(out) :: ecsr

! Compute PW LDA correlation energy

  call ecPW(rs,0d0,ec)

! Define various stuff

  cf = (9d0*pi/4d0)**(1d0/3d0)
  bb = 1.27329d0    
  m1 = 0.0357866d0  
  a0 = ec
  b3 = bb*rs**(7d0/2d0)
  b2 = -3d0/2d0/pi/cf*rs/a0
  b1 = (b3-1d0/sqrt(3d0*pi)*rs**(3d0/2d0)/a0)/b2
  b4 = -a0*b1*rs**3/m1

! Compute short-range correlation energy

  ecsr = a0*(1d0 + b1*mu)/(1d0 + b1*mu+b2*mu**2 + b3*mu**3 + b4*mu**4)

end subroutine srlda

!==========================================================================================

subroutine ecPW(x,y,ec)

! Correlation energy of the 3D electron gas of density rs and spin polarization z
! Perdew & Wang, PRB 45, 13244 (1992)

  implicit none
  include 'parameters.h'

! Input variables 

  double precision,intent(in)  :: x
  double precision,intent(in)  :: y

! Local variables 

  double precision             :: f02
  double precision             :: ff
  double precision             :: aaa
  double precision             :: G
  double precision             :: ec0
  double precision             :: ec1
  double precision             :: alfac
  
! Ouput variables

  double precision,intent(out) :: ec

  f02 = 4d0/(9d0*(2d0**(1d0/3d0) - 1d0))

  ff = ((1d0+y)**(4d0/3d0) + (1d0-y)**(4d0/3d0)-2d0)/(2d0**(4d0/3d0) - 2d0)

  aaa = (1d0 - log(2d0))/pi**2

  call GPW(x,aaa,0.21370d0,7.5957d0,3.5876d0,1.6382d0,0.49294d0,G)
  ec0 = G

  aaa=aaa/2d0
  call GPW(x,aaa,0.20548d0,14.1189d0,6.1977d0,3.3662d0,0.62517d0,G)
  ec1 = G

  call GPW(x,0.016887d0,0.11125d0,10.357d0,3.6231d0,0.88026d0,0.49671d0,G)
  alfac = -G

  ec = ec0 + alfac*ff/f02*(1d0 - y**4) + (ec1 - ec0)*ff*y**4

end subroutine ecPW

!==========================================================================================

subroutine GPW(x,Ac,alfa1,beta1,beta2,beta3,beta4,G)

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)  :: Ac
  double precision,intent(in)  :: alfa1
  double precision,intent(in)  :: beta1
  double precision,intent(in)  :: beta2
  double precision,intent(in)  :: beta3
  double precision,intent(in)  :: beta4
  double precision,intent(in)  :: x

! Ouput variables

  double precision,intent(out) :: G

  G = -2d0*Ac*(1d0 + alfa1*x)*log(1d0 & 
    + 1d0/(2d0*Ac*(beta1*sqrt(x) + beta2*x + beta3*x*sqrt(x) + beta4*x**2)))

end subroutine GPW
