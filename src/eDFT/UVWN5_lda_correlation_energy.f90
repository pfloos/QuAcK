subroutine UVWN5_lda_correlation_energy(nGrid,weight,rho,Ec)

! Compute unrestricted VWN5 LDA correlation energy

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  integer                       :: iG
  double precision              :: ra,rb,r,rs,x,z
  double precision              :: a_p,x0_p,xx0_p,b_p,c_p,x_p,q_p
  double precision              :: a_f,x0_f,xx0_f,b_f,c_f,x_f,q_f
  double precision              :: a_a,x0_a,xx0_a,b_a,c_a,x_a,q_a

  double precision              :: ec_z,ec_p,ec_f,ec_a
  double precision              :: fz,d2fz

! Output variables

  double precision              :: Ec(nsp)

! Parameters of the functional

  a_p  = +0.0621814D0/2D0
  x0_p = -0.10498d0
  b_p  = +3.72744d0
  c_p  = +12.9352d0

  a_f  = +0.0621814D0/4D0
  x0_f = -0.325d0
  b_f  = +7.06042d0
  c_f  = +18.0578d0

  a_a  = -1d0/(6d0*pi**2)
  x0_a = -0.0047584D0
  b_a  = 1.13107d0
  c_a  = 13.0045d0

! Initialization

  Ec(:) = 0d0

  do iG=1,nGrid

    ra = max(0d0,rho(iG,1))
    rb = max(0d0,rho(iG,2))
    r = ra + rb
    z = (ra - rb)/r

!   alpha-alpha contribution

    if(ra > threshold) then

      rs = (4d0*pi*ra/3d0)**(-1d0/3d0)
      x = sqrt(rs)

      x_f   = x*x + b_f*x + c_f
      xx0_f = x0_f*x0_f + b_f*x0_f + c_f
      q_f   = sqrt(4d0*c_f - b_f*b_f)

      ec_f = a_f*( log(x**2/x_f) + 2d0*b_f/q_f*atan(q_f/(2d0*x + b_f)) &
                 - b_f*x0_f/xx0_f*( log((x - x0_f)**2/x_f) + 2d0*(b_f + 2d0*x0_f)/q_f*atan(q_f/(2d0*x + b_f)) ) )

      Ec(1) = Ec(1) + weight(iG)*ec_f*ra

    end if

!   alpha-beta contribution

    if(r > threshold) then

      rs = (4d0*pi*r/3d0)**(-1d0/3d0)
      x = sqrt(rs)

      fz = (1d0 + z)**(4d0/3d0) + (1d0 - z)**(4d0/3d0) - 2d0
      fz = fz/(2d0*(2d0**(1d0/3d0) - 1d0))

      d2fz = 4d0/(9d0*(2**(1d0/3d0) - 1d0))

      x_p = x*x + b_p*x + c_p
      x_f = x*x + b_f*x + c_f
      x_a = x*x + b_a*x + c_a

      xx0_p = x0_p*x0_p + b_p*x0_p + c_p
      xx0_f = x0_f*x0_f + b_f*x0_f + c_f
      xx0_a = x0_a*x0_a + b_a*x0_a + c_a

      q_p = sqrt(4d0*c_p - b_p*b_p)
      q_f = sqrt(4d0*c_f - b_f*b_f)
      q_a = sqrt(4d0*c_a - b_a*b_a)

      ec_p = a_p*( log(x**2/x_p) + 2d0*b_p/q_p*atan(q_p/(2d0*x + b_p)) & 
                 - b_p*x0_p/xx0_p*( log((x - x0_p)**2/x_p) + 2d0*(b_p + 2d0*x0_p)/q_p*atan(q_p/(2d0*x + b_p)) ) )

      ec_f = a_f*( log(x**2/x_f) + 2d0*b_f/q_f*atan(q_f/(2d0*x + b_f)) & 
                 - b_f*x0_f/xx0_f*( log((x - x0_f)**2/x_f) + 2d0*(b_f + 2d0*x0_f)/q_f*atan(q_f/(2d0*x + b_f)) ) )

      ec_a = a_a*( log(x**2/x_a) + 2d0*b_a/q_a*atan(q_a/(2d0*x + b_a)) & 
                 - b_a*x0_a/xx0_a*( log((x - x0_a)**2/x_a) + 2d0*(b_a + 2d0*x0_a)/q_a*atan(q_a/(2d0*x + b_a)) ) )

      ec_z = ec_p + ec_a*fz/d2fz*(1d0-z**4) + (ec_f - ec_p)*fz*z**4

      Ec(2) = Ec(2) + weight(iG)*ec_z*r

    end if

!   beta-beta contribution

    if(rb > threshold) then

      rs = (4d0*pi*rb/3d0)**(-1d0/3d0)
      x = sqrt(rs)

      x_f   = x*x + b_f*x + c_f
      xx0_f = x0_f*x0_f + b_f*x0_f + c_f
      q_f   = sqrt(4d0*c_f - b_f*b_f)

      ec_f = a_f*( log(x**2/x_f) + 2d0*b_f/q_f*atan(q_f/(2d0*x + b_f)) &
                 - b_f*x0_f/xx0_f*( log((x - x0_f)**2/x_f) + 2d0*(b_f + 2d0*x0_f)/q_f*atan(q_f/(2d0*x + b_f)) ) )

      Ec(3) = Ec(3) + weight(iG)*ec_f*rb

    end if

  end do

  Ec(2) = Ec(2) - Ec(1) - Ec(3)

end subroutine UVWN5_lda_correlation_energy
