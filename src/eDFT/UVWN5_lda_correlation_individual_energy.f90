subroutine UVWN5_lda_correlation_individual_energy(nGrid,weight,rhow,rho,doNcentered,kappa,Ec)

! Compute VWN5 LDA correlation potential

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: rho(nGrid,nspin)
  logical,intent(in)            :: doNcentered
  double precision,intent(in)   :: kappa

! Local variables

  integer                       :: iG
  double precision              :: ra,rb,r,raI,rbI,rI,rs,x,z
  double precision              :: a_p,x0_p,xx0_p,b_p,c_p,x_p,q_p
  double precision              :: a_f,x0_f,xx0_f,b_f,c_f,x_f,q_f
  double precision              :: a_a,x0_a,xx0_a,b_a,c_a,x_a,q_a
  double precision              :: dfzdz,dxdrs,dxdx_p,dxdx_f,dxdx_a,decdx_p,decdx_f,decdx_a
  double precision              :: dzdra,dzdrb,dfzdra,dfzdrb,drsdr,decdr_p,decdr_f,decdr_a,decdra,decdrb,decdr
  double precision              :: ec_z,ec_p,ec_f,ec_a
  double precision              :: fz,d2fz

  double precision              :: Ecrr(nsp)
  double precision              :: EcrI(nsp)
  double precision              :: EcrrI(nsp)

! Output variables

  double precision              :: Ec(nsp)

! Parameters of the functional

  a_p  = +0.0621814d0/2d0
  x0_p = -0.10498d0
  b_p  = +3.72744d0
  c_p  = +12.9352d0

  a_f  = +0.0621814d0/4d0
  x0_f = -0.325d0
  b_f  = +7.06042d0
  c_f  = +18.0578d0

  a_a  = -1d0/(6d0*pi**2)
  x0_a = -0.0047584d0
  b_a  = +1.13107d0
  c_a  = +13.0045d0

! Initialization

  Ec(:)    = 0d0
  Ecrr(:)  = 0d0
  EcrI(:)  = 0d0
  EcrrI(:) = 0d0

  do iG=1,nGrid

    ra = max(0d0,rhow(iG,1))
    rb = max(0d0,rhow(iG,2))

    raI = max(0d0,rho(iG,1))
    rbI = max(0d0,rho(iG,2))
   
    r  = ra + rb
    rI = raI + rbI

!   spin-up contribution

 !  if(r > threshold) then

 !    rs = (4d0*pi*r/3d0)**(-1d0/3d0)
 !    x = sqrt(rs)
 ! 
 !    x_f   = x*x + b_f*x + c_f
 !    xx0_f = x0_f*x0_f + b_f*x0_f + c_f
 !    q_f   = sqrt(4d0*c_f - b_f*b_f)
 ! 
 !    ec_f = a_f*( log(x**2/x_f) + 2d0*b_f/q_f*atan(q_f/(2d0*x + b_f)) & 
 !         - b_f*x0_f/xx0_f*( log((x - x0_f)**2/x_f) + 2d0*(b_f + 2d0*x0_f)/q_f*atan(q_f/(2d0*x + b_f)) ) )
 ! 
 !    drsdr = - (36d0*pi)**(-1d0/3d0)*r**(-4d0/3d0)
 !    dxdrs = 0.5d0/sqrt(rs)

 !    dxdx_f = 2d0*x + b_f

 !    decdx_f = a_f*( 2d0/x - 4d0*b_f/( (b_f+2d0*x)**2 + q_f**2) - dxdx_f/x_f &
 !            - b_f*x0_f/xx0_f*( 2/(x-x0_f) - 4d0*(b_f+2d0*x0_f)/( (b_f+2d0*x)**2 + q_f**2) - dxdx_f/x_f ) )

 !    decdr_f = drsdr*dxdrs*decdx_f

 !    Ecrr(1)  = Ecrr(1)  - weight(iG)*decdr_f*r*r

 !    if(rI > threshold) then 

 !      EcrI(1)  = EcrI(1)  + weight(iG)*ec_f*rI
 !      EcrrI(1) = EcrrI(1) + weight(iG)*decdr_f*r*rI

 !    end if
 ! 
 !  end if

!   up-down contribution

    if(r > threshold) then

      rs = (4d0*pi*r/3d0)**(-1d0/3d0)
      z = (ra - rb)/r
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
   
      ec_z = ec_p + ec_a*fz/d2fz*(1d0 - z**4) + (ec_f - ec_p)*fz*z**4

      dfzdz = (4d0/3d0)*((1d0 + z)**(1d0/3d0) - (1d0 - z)**(1d0/3d0))/(2d0*(2d0**(1d0/3d0) - 1d0))

      dzdra = + (1d0 - z)/r
      dfzdra = dzdra*dfzdz

      dzdrb = - (1d0 + z)/r
      dfzdrb = dzdrb*dfzdz

      drsdr = - (36d0*pi)**(-1d0/3d0)*r**(-4d0/3d0)
      dxdrs = 0.5d0/sqrt(rs)

      dxdx_p = 2d0*x + b_p
      dxdx_f = 2d0*x + b_f
      dxdx_a = 2d0*x + b_a

      decdx_p = a_p*( 2d0/x - 4d0*b_p/( (b_p+2d0*x)**2 + q_p**2) - dxdx_p/x_p &
                      - b_p*x0_p/xx0_p*( 2d0/(x-x0_p) - 4d0*(b_p+2d0*x0_p)/( (b_p+2d0*x)**2 + q_p**2) - dxdx_p/x_p ) )

      decdx_f = a_f*( 2d0/x - 4d0*b_f/( (b_f+2d0*x)**2 + q_f**2) - dxdx_f/x_f &
                      - b_f*x0_f/xx0_f*( 2d0/(x-x0_f) - 4d0*(b_f+2d0*x0_f)/( (b_f+2d0*x)**2 + q_f**2) - dxdx_f/x_f ) )

      decdx_a = a_a*( 2d0/x - 4d0*b_a/( (b_a+2d0*x)**2 + q_a**2) - dxdx_a/x_a &
                      - b_a*x0_a/xx0_a*( 2d0/(x-x0_a) - 4d0*(b_a+2d0*x0_a)/( (b_a+2d0*x)**2 + q_a**2) - dxdx_a/x_a ) )

      decdr_p = drsdr*dxdrs*decdx_p
      decdr_f = drsdr*dxdrs*decdx_f
      decdr_a = drsdr*dxdrs*decdx_a

      decdra = decdr_p + decdr_a*fz/d2fz*(1d0-z**4) + ec_a*dfzdra/d2fz*(1d0-z**4) - 4d0*ec_a*fz/d2fz*dzdra*z**3 &
             + (decdr_f - decdr_p)*fz*z**4 + (ec_f - ec_p)*dfzdra*z**4 + 4d0*(ec_f - ec_p)*fz*dzdra*z**3

      decdrb = decdr_p + decdr_a*fz/d2fz*(1d0-z**4) + ec_a*dfzdrb/d2fz*(1d0-z**4) - 4d0*ec_a*fz/d2fz*dzdrb*z**3 &
             + (decdr_f - decdr_p)*fz*z**4 + (ec_f - ec_p)*dfzdrb*z**4 + 4d0*(ec_f - ec_p)*fz*dzdrb*z**3
      
      decdr = 0d0 
      if(ra > threshold) decdr = decdr + decdra
      if(rb > threshold) decdr = decdr + decdrb
 
      Ecrr(2)  = Ecrr(2) - weight(iG)*decdr*r*r

      if(rI > threshold) then

        EcrI(2)  = EcrI(2)  + weight(iG)*ec_z*rI
        EcrrI(2) = EcrrI(2) + weight(iG)*decdr*r*rI

      end if
   
    end if

!   spin-down contribution
   
 !  if(r > threshold) then

 !    rs = (4d0*pi*r/3d0)**(-1d0/3d0)
 !    x  = sqrt(rs)

 !    x_f   = x*x + b_f*x + c_f
 !    xx0_f = x0_f*x0_f + b_f*x0_f + c_f
 !    q_f   = sqrt(4d0*c_f - b_f*b_f)

 !    ec_f = a_f*( log(x**2/x_f) + 2d0*b_f/q_f*atan(q_f/(2d0*x + b_f)) &
 !               - b_f*x0_f/xx0_f*( log((x - x0_f)**2/x_f) + 2d0*(b_f + 2d0*x0_f)/q_f*atan(q_f/(2d0*x + b_f)) ) )

 !    drsdr = - (36d0*pi)**(-1d0/3d0)*r**(-4d0/3d0)
 !    dxdrs = 0.5d0/sqrt(rs)

 !    dxdx_f = 2d0*x + b_f

 !    decdx_f = a_f*( 2d0/x - 4d0*b_f/( (b_f+2d0*x)**2 + q_f**2) - dxdx_f/x_f &
 !            - b_f*x0_f/xx0_f*( 2/(x-x0_f) - 4d0*(b_f+2d0*x0_f)/( (b_f+2d0*x)**2 + q_f**2) - dxdx_f/x_f ) )

 !    decdr_f = drsdr*dxdrs*decdx_f

 !    Ecrr(3)  = Ecrr(3)  - weight(iG)*decdr_f*r*r

 !    if(rI > threshold) then

 !      EcrI(3)  = EcrI(3)  + weight(iG)*ec_f*rI
 !      EcrrI(3) = EcrrI(3) + weight(iG)*decdr_f*r*rI

 !    end if

 !  end if

  end do

  Ecrr(2)  = Ecrr(2)  - Ecrr(1)  - Ecrr(3)
  EcrI(2)  = EcrI(2)  - EcrI(1)  - EcrI(3)
  EcrrI(2) = EcrrI(2) - EcrrI(1) - EcrrI(3)

! De-scaling for N-centered ensemble

  if(doNcentered) then

    Ecrr(:)  = kappa*Ecrr(:)
    EcrI(:)  = kappa*EcrI(:)

  endif

  Ec(:) = Ecrr(:) + EcrI(:) + EcrrI(:)

end subroutine UVWN5_lda_correlation_individual_energy
