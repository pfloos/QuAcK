subroutine UVWN5_lda_correlation_potential(nGrid,weight,nBas,AO,rho,Fc)

! Compute unrestricted VWN5 LDA correlation potential

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: ra,rb,r,rs,x,z
  double precision              :: a_p,x0_p,xx0_p,b_p,c_p,x_p,q_p
  double precision              :: a_f,x0_f,xx0_f,b_f,c_f,x_f,q_f
  double precision              :: a_a,x0_a,xx0_a,b_a,c_a,x_a,q_a
  double precision              :: dfzdz,dxdrs,dxdx_p,dxdx_f,dxdx_a,decdx_p,decdx_f,decdx_a
  double precision              :: dzdra,dfzdra,drsdra,decdra_p,decdra_f,decdra_a,decdra
  double precision              :: dzdrb,dfzdrb,drsdrb,decdrb_p,decdrb_f,decdrb_a,decdrb

  double precision              :: ec_z,ec_p,ec_f,ec_a
  double precision              :: fz,d2fz

! Output variables

  double precision              :: Fc(nBas,nBas,nspin)

! Parameters of the functional

  a_p  = +0.0621814D0/2d0
  x0_p = -0.10498d0
  b_p  = +3.72744d0
  c_p  = +12.9352d0

  a_f  = +0.0621814D0/4D0
  x0_f = -0.325d0
  b_f  = +7.06042d0
  c_f  = +18.0578d0

  a_a  = -1d0/(6d0*pi**2)
  x0_a = -0.0047584D0
  b_a  = +1.13107d0
  c_a  = +13.0045d0

! Initialization

  Fc(:,:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        ra = max(0d0,rho(iG,1))
        rb = max(0d0,rho(iG,2))
        r  = ra + rb
        z  = (ra - rb)/r

        rs = (4d0*pi*r/3d0)**(-1d0/3d0)
        x  = sqrt(rs)
       
        fz = (1d0 + z)**(4d0/3d0) + (1d0 - z)**(4d0/3d0) - 2d0
        fz = fz/(2d0*(2d0**(1d0/3d0) - 1d0))
       
        d2fz = 4d0/(9d0*(2**(1d0/3d0) - 1d0))
       
!       spin-up contribution
       
        if(ra > threshold) then
       
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

          dzdra = + (1d0 - z)/r
          dfzdz = (4d0/3d0)*((1d0 + z)**(1d0/3d0) - (1d0 - z)**(1d0/3d0))/(2d0*(2d0**(1d0/3d0) - 1d0))
          dfzdra = dzdra*dfzdz

          drsdra = - (36d0*pi)**(-1d0/3d0)*r**(-4d0/3d0)
          dxdrs  = 0.5d0/sqrt(rs)

          dxdx_p = 2d0*x + b_p
          dxdx_f = 2d0*x + b_f
          dxdx_a = 2d0*x + b_a

          decdx_p = a_p*( 2d0/x - 4d0*b_p/( (b_p+2d0*x)**2 + q_p**2) - dxdx_p/x_p &
                          - b_p*x0_p/xx0_p*( 2d0/(x-x0_p) - 4d0*(b_p+2d0*x0_p)/( (b_p+2d0*x)**2 + q_p**2) - dxdx_p/x_p ) )

          decdx_f = a_f*( 2d0/x - 4d0*b_f/( (b_f+2d0*x)**2 + q_f**2) - dxdx_f/x_f &
                          - b_f*x0_f/xx0_f*( 2d0/(x-x0_f) - 4d0*(b_f+2d0*x0_f)/( (b_f+2d0*x)**2 + q_f**2) - dxdx_f/x_f ) )

          decdx_a = a_a*( 2d0/x - 4d0*b_a/( (b_a+2d0*x)**2 + q_a**2) - dxdx_a/x_a &
                          - b_a*x0_a/xx0_a*( 2d0/(x-x0_a) - 4d0*(b_a+2d0*x0_a)/( (b_a+2d0*x)**2 + q_a**2) - dxdx_a/x_a ) )

          decdra_p = drsdra*dxdrs*decdx_p
          decdra_f = drsdra*dxdrs*decdx_f
          decdra_a = drsdra*dxdrs*decdx_a

          decdra = decdra_p + decdra_a*fz/d2fz*(1d0-z**4) + ec_a*dfzdra/d2fz*(1d0-z**4) - 4d0*ec_a*fz/d2fz*dzdra*z**3 &
                 + (decdra_f - decdra_p)*fz*z**4 + (ec_f - ec_p)*dfzdra*z**4 + 4d0*(ec_f - ec_p)*fz*dzdra*z**3

          Fc(mu,nu,1) = Fc(mu,nu,1) + weight(iG)*AO(mu,iG)*AO(nu,iG)*(ec_z + decdra*r)
       
        end if

!       spin-down contribution
       
        if(rb > threshold) then
       
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
       
          dzdrb = - (1d0 + z)/r
          dfzdz = (4d0/3d0)*((1d0 + z)**(1d0/3d0) - (1d0 - z)**(1d0/3d0))/(2d0*(2d0**(1d0/3d0) - 1d0))
          dfzdrb = dzdrb*dfzdz

          drsdrb = - (36d0*pi)**(-1d0/3d0)*r**(-4d0/3d0)
          dxdrs  = 0.5d0/sqrt(rs)

          dxdx_p = 2d0*x + b_p
          dxdx_f = 2d0*x + b_f
          dxdx_a = 2d0*x + b_a

          decdx_p = a_p*( 2d0/x - 4d0*b_p/( (b_p+2d0*x)**2 + q_p**2) - dxdx_p/x_p &
                          - b_p*x0_p/xx0_p*( 2d0/(x-x0_p) - 4d0*(b_p+2d0*x0_p)/( (b_p+2d0*x)**2 + q_p**2) - dxdx_p/x_p ) )

          decdx_f = a_f*( 2d0/x - 4d0*b_f/( (b_f+2d0*x)**2 + q_f**2) - dxdx_f/x_f &
                          - b_f*x0_f/xx0_f*( 2d0/(x-x0_f) - 4d0*(b_f+2d0*x0_f)/( (b_f+2d0*x)**2 + q_f**2) - dxdx_f/x_f ) )

          decdx_a = a_a*( 2d0/x - 4d0*b_a/( (b_a+2d0*x)**2 + q_a**2) - dxdx_a/x_a &
                          - b_a*x0_a/xx0_a*( 2d0/(x-x0_a) - 4d0*(b_a+2d0*x0_a)/( (b_a+2d0*x)**2 + q_a**2) - dxdx_a/x_a ) )

          decdrb_p = drsdrb*dxdrs*decdx_p
          decdrb_f = drsdrb*dxdrs*decdx_f
          decdrb_a = drsdrb*dxdrs*decdx_a

          decdrb = decdrb_p + decdrb_a*fz/d2fz*(1d0-z**4) + ec_a*dfzdrb/d2fz*(1d0-z**4) - 4d0*ec_a*fz/d2fz*dzdrb*z**3 &
                 + (decdrb_f - decdrb_p)*fz*z**4 + (ec_f - ec_p)*dfzdrb*z**4 + 4d0*(ec_f - ec_p)*fz*dzdrb*z**3

          Fc(mu,nu,2) = Fc(mu,nu,2) + weight(iG)*AO(mu,iG)*AO(nu,iG)*(ec_z + decdrb*r)
       
        end if

      end do
    end do
  end do

end subroutine UVWN5_lda_correlation_potential
