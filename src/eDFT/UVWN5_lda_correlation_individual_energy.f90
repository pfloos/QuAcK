subroutine UVWN5_lda_correlation_individual_energy(nEns,nGrid,weight,rhow,rho,doNcentered,LZc,Ec)

! Compute VWN5 LDA correlation potential

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)
  double precision,intent(in)   :: rho(nGrid,nspin,nEns)
  logical,intent(in)            :: doNcentered

! Local variables

  integer                       :: iG
  integer                       :: iEns
  double precision              :: ra,rb,r,raI,rbI,rI,rs,x,z
  double precision              :: a_p,x0_p,xx0_p,b_p,c_p,x_p,q_p
  double precision              :: a_f,x0_f,xx0_f,b_f,c_f,x_f,q_f
  double precision              :: a_a,x0_a,xx0_a,b_a,c_a,x_a,q_a
  double precision              :: dfzdz,dxdrs,dxdx_p,dxdx_f,dxdx_a,decdx_p,decdx_f,decdx_a
  double precision              :: dzdra,dzdrb,dfzdra,dfzdrb,drsdr,decdr_p,decdr_f,decdr_a,decdra,decdrb
  double precision              :: ec_z,ec_p,ec_f,ec_a
  double precision              :: fz,d2fz

! Output variables

  double precision,intent(out)  :: LZc(nsp)
  double precision,intent(out)  :: Ec(nsp,nEns)

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

  LZc(:)  = 0d0
  Ec(:,:) = 0d0

  do iG=1,nGrid

    ra = max(0d0,rhow(iG,1))
    rb = max(0d0,rhow(iG,2))

    r  = ra + rb

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

      dzdra = + (1d0 - z)/r
      dfzdra = dzdra*dfzdz

      decdra = decdr_p + decdr_a*fz/d2fz*(1d0-z**4) + ec_a*dfzdra/d2fz*(1d0-z**4) - 4d0*ec_a*fz/d2fz*dzdra*z**3 &
             + (decdr_f - decdr_p)*fz*z**4 + (ec_f - ec_p)*dfzdra*z**4 + 4d0*(ec_f - ec_p)*fz*dzdra*z**3

      dzdrb = - (1d0 + z)/r
      dfzdrb = dzdrb*dfzdz

      decdrb = decdr_p + decdr_a*fz/d2fz*(1d0-z**4) + ec_a*dfzdrb/d2fz*(1d0-z**4) - 4d0*ec_a*fz/d2fz*dzdrb*z**3 &
             + (decdr_f - decdr_p)*fz*z**4 + (ec_f - ec_p)*dfzdrb*z**4 + 4d0*(ec_f - ec_p)*fz*dzdrb*z**3
      
      ! spin-up contribution

      if(ra > threshold) then 
 
        LZc(1) = LZc(1) - weight(iG)*decdra*ra*ra
        
        do iEns=1,nEns

          raI = max(0d0,rho(iG,1,iEns))

          if(raI > threshold) then 

            Ec(1,iEns) = Ec(1,iEns) + weight(iG)*(ec_z + decdra*ra)*raI
            Ec(2,iEns) = Ec(2,iEns) + weight(iG)*decdra*rb*raI

          end if

        end do

      end if

      ! spin-down contribution

      if(rb > threshold) then 

        LZc(3) = LZc(3) - weight(iG)*decdrb*rb*rb
        
        do iEns=1,nEns

          rbI = max(0d0,rho(iG,2,iEns))

          if(rbI > threshold) then 

            Ec(3,iEns) = Ec(3,iEns) + weight(iG)*(ec_z + decdrb*rb)*rbI
            Ec(2,iEns) = Ec(2,iEns) + weight(iG)*decdrb*ra*rbI

          end if

        end do

      end if

    end if

  end do

end subroutine UVWN5_lda_correlation_individual_energy
