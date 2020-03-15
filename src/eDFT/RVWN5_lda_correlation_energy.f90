subroutine RVWN5_lda_correlation_energy(nGrid,weight,rho,Ec)

! Compute the restricted VWN5 LDA correlation energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r,rs,x
  double precision              :: a_p,x0_p,xx0_p,b_p,c_p,x_p,q_p

  double precision              :: ec_p

! Output variables

  double precision              :: Ec(nsp)

! Parameters of the functional

  a_p  = +0.0621814D0/2D0
  x0_p = -0.10498d0
  b_p  = +3.72744d0
  c_p  = +12.9352d0

! Initialization

  Ec = 0d0

  do iG=1,nGrid

    r = max(0d0,rho(iG))

    if(r > threshold) then

      rs = (4d0*pi*r/3d0)**(-1d0/3d0)
      x = sqrt(rs)

      x_p = x*x + b_p*x + c_p

      xx0_p = x0_p*x0_p + b_p*x0_p + c_p

      q_p = sqrt(4d0*c_p - b_p*b_p)

      ec_p = a_p*( log(x**2/x_p) + 2d0*b_p/q_p*atan(q_p/(2d0*x + b_p)) & 
                 - b_p*x0_p/xx0_p*( log((x - x0_p)**2/x_p) + 2d0*(b_p + 2d0*x0_p)/q_p*atan(q_p/(2d0*x + b_p)) ) )

      Ec = Ec + weight(iG)*ec_p*r

    end if

  end do

end subroutine RVWN5_lda_correlation_energy
