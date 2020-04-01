subroutine RVWN5_lda_correlation_individual_energy(nGrid,weight,rhow,rho,Ec)

! Compute the restricted VWN5 LDA correlation individual energies

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: iG
  double precision              :: r,rI,rs,x
  double precision              :: a_p,x0_p,xx0_p,b_p,c_p,x_p,q_p
  double precision              :: dxdrs,dxdx_p,decdx_p
  double precision              :: drsdr,decdr_p
  double precision              :: ec_p

! Output variables

  double precision              :: Ec

! Parameters of the functional

  a_p  = +0.0621814D0/2D0
  x0_p = -0.10498d0
  b_p  = +3.72744d0
  c_p  = +12.9352d0

! Initialization

  Ec = 0d0

  do iG=1,nGrid

    r = max(0d0,rhow(iG))
    rI = max(0d0,rho(iG))
   
    if(r > threshold .and. rI > threshold) then
   
      rs = (4d0*pi*r/3d0)**(-1d0/3d0)
      x = sqrt(rs)
   
      x_p = x*x + b_p*x + c_p
      xx0_p = x0_p*x0_p + b_p*x0_p + c_p
      q_p = sqrt(4d0*c_p - b_p*b_p)

      ec_p = a_p*( log(x**2/x_p) + 2d0*b_p/q_p*atan(q_p/(2d0*x + b_p)) & 
                 - b_p*x0_p/xx0_p*( log((x - x0_p)**2/x_p) + 2d0*(b_p + 2d0*x0_p)/q_p*atan(q_p/(2d0*x + b_p)) ) )
   
      drsdr = - (36d0*pi)**(-1d0/3d0)*r**(-4d0/3d0)
      dxdrs = 0.5d0/sqrt(rs)

      dxdx_p = 2d0*x + b_p

      decdx_p = a_p*( 2d0/x - 4d0*b_p/( (b_p+2d0*x)**2 + q_p**2) - dxdx_p/x_p &
                      - b_p*x0_p/xx0_p*( 2/(x-x0_p) - 4d0*(b_p+2d0*x0_p)/( (b_p+2d0*x)**2 + q_p**2) - dxdx_p/x_p ) )

      decdr_p = drsdr*dxdrs*decdx_p

      Ec = Ec + weight(iG)*(ec_p*rI + decdr_p*r*rI - decdr_p*r*r)
   
    end if

  end do

end subroutine RVWN5_lda_correlation_individual_energy
