subroutine UPW92_lda_correlation_potential(nGrid,weight,nBas,AO,rho,Fc)

! Compute unrestricted PW92 LDA correlation potential

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
  double precision              :: ra,rb,r,rs,z,t,dt
  double precision              :: A_p,a1_p,b1_p,b2_p,b3_p,b4_p
  double precision              :: A_f,a1_f,b1_f,b2_f,b3_f,b4_f
  double precision              :: A_a,a1_a,b1_a,b2_a,b3_a,b4_a
  double precision              :: dfzdz,decdrs_p,decdrs_f,decdrs_a
  double precision              :: dzdra,dfzdra,drsdra,decdra_p,decdra_f,decdra_a,decdra
  double precision              :: dzdrb,dfzdrb,drsdrb,decdrb_p,decdrb_f,decdrb_a,decdrb

  double precision              :: ec_z,ec_p,ec_f,ec_a
  double precision              :: fz,d2fz

! Output variables

  double precision              :: Fc(nBas,nBas,nspin)

! Parameters of the functional

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

! Initialization

  Fc(:,:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        ra = max(0d0,rho(iG,1))
        rb = max(0d0,rho(iG,2))
       
!       spin-up contribution
       
        if(ra > threshold) then
       
          r = ra + rb
          rs = (4d0*pi*r/3d0)**(-1d0/3d0)
          z = (ra - rb)/r
       
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

          dzdra = (1d0 - z)/r
          dfzdz = (4d0/3d0)*((1d0 + z)**(1d0/3d0) - (1d0 - z)**(1d0/3d0))/(2d0*(2d0**(1d0/3d0) - 1d0))
          dfzdra = dzdra*dfzdz
          drsdra = - (36d0*pi)**(-1d0/3d0)*r**(-4d0/3d0)

          t  = b1_p*sqrt(rs) + b2_p*rs + b3_p*rs**(3d0/2d0) + b4_p*rs**2
          dt = 0.5d0*b1_p*sqrt(rs) + b2_p*rs + 1.5d0*b3_p*rs**(3d0/2d0) + 2d0*b4_p*rs**2
          decdrs_p = (1d0 + a1_p*rs)/rs*dt/(t**2*(1d0 + 1d0/(2d0*A_p*t))) &
                   - 2d0*A_p*a1_p*log(1d0 + 1d0/(2d0*A_p*t))

          t  = b1_f*sqrt(rs) + b2_f*rs + b3_f*rs**(3d0/2d0) + b4_f*rs**2
          dt = 0.5d0*b1_f*sqrt(rs) + b2_f*rs + 1.5d0*b3_f*rs**(3d0/2d0) + 2d0*b4_f*rs**2
          decdrs_f = (1d0 + a1_f*rs)/rs*dt/(t**2*(1d0 + 1d0/(2d0*A_f*t))) &
                   - 2d0*A_f*a1_f*log(1d0 + 1d0/(2d0*A_f*t))

          t  = b1_a*sqrt(rs) + b2_a*rs + b3_a*rs**(3d0/2d0) + b4_a*rs**2
          dt = 0.5d0*b1_a*sqrt(rs) + b2_a*rs + 1.5d0*b3_a*rs**(3d0/2d0) + 2d0*b4_a*rs**2
          decdrs_a = (1d0 + a1_a*rs)/rs*dt/(t**2*(1d0 + 1d0/(2d0*A_a*t))) &
                   - 2d0*A_a*a1_a*log(1d0 + 1d0/(2d0*A_a*t))

          decdra_p = drsdra*decdrs_p
          decdra_f = drsdra*decdrs_f
          decdra_a = drsdra*decdrs_a

          decdra = decdra_p + decdra_a*fz/d2fz*(1d0-z**4) + ec_a*dfzdra/d2fz*(1d0-z**4) - 4d0*ec_a*fz/d2fz*dzdra*z**3 &
                 + (decdra_f - decdra_p)*fz*z**4 + (ec_f - ec_p)*dfzdra*z**4 + 4d0*(ec_f - ec_p)*fz*dzdra*z**3

          Fc(mu,nu,1) = Fc(mu,nu,1) + weight(iG)*AO(mu,iG)*AO(nu,iG)*(ec_z + decdra*r)
       
        end if

!       spin-down contribution
       
        if(rb > threshold) then
       
          r = ra + rb
          rs = (4d0*pi*r/3d0)**(-1d0/3d0)
          z = (ra - rb)/r
       
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
       
          dzdrb = - (1d0 + z)/r
          dfzdz = (4d0/3d0)*((1d0 + z)**(1d0/3d0) - (1d0 - z)**(1d0/3d0))/(2d0*(2d0**(1d0/3d0) - 1d0))
          dfzdrb = dzdrb*dfzdz

          drsdrb = - (36d0*pi)**(-1d0/3d0)*r**(-4d0/3d0)


          t  = b1_p*sqrt(rs) + b2_p*rs + b3_p*rs**(3d0/2d0) + b4_p*rs**2
          dt = 0.5d0*b1_p*sqrt(rs) + b2_p*rs + 1.5d0*b3_p*rs**(3d0/2d0) + 2d0*b4_p*rs**2
          decdrs_p = (1d0 + a1_p*rs)/rs*dt/(t**2*(1d0 + 1d0/(2d0*A_p*t))) &
                   - 2d0*A_p*a1_p*log(1d0 + 1d0/(2d0*A_p*t))

          t  = b1_f*sqrt(rs) + b2_f*rs + b3_f*rs**(3d0/2d0) + b4_f*rs**2
          dt = 0.5d0*b1_f*sqrt(rs) + b2_f*rs + 1.5d0*b3_f*rs**(3d0/2d0) + 2d0*b4_f*rs**2
          decdrs_f = (1d0 + a1_f*rs)/rs*dt/(t**2*(1d0 + 1d0/(2d0*A_f*t))) &
                   - 2d0*A_f*a1_f*log(1d0 + 1d0/(2d0*A_f*t))

          t  = b1_a*sqrt(rs) + b2_a*rs + b3_a*rs**(3d0/2d0) + b4_a*rs**2
          dt = 0.5d0*b1_a*sqrt(rs) + b2_a*rs + 1.5d0*b3_a*rs**(3d0/2d0) + 2d0*b4_a*rs**2
          decdrs_a = (1d0 + a1_a*rs)/rs*dt/(t**2*(1d0 + 1d0/(2d0*A_a*t))) &
                   - 2d0*A_a*a1_a*log(1d0 + 1d0/(2d0*A_a*t))

          decdrb_p = drsdrb*decdrs_p
          decdrb_f = drsdrb*decdrs_f
          decdrb_a = drsdrb*decdrs_a

          decdrb = decdrb_p + decdrb_a*fz/d2fz*(1d0-z**4) + ec_a*dfzdrb/d2fz*(1d0-z**4) - 4d0*ec_a*fz/d2fz*dzdrb*z**3 &
                 + (decdrb_f - decdrb_p)*fz*z**4 + (ec_f - ec_p)*dfzdrb*z**4 + 4d0*(ec_f - ec_p)*fz*dzdrb*z**3

          Fc(mu,nu,2) = Fc(mu,nu,2) + weight(iG)*AO(mu,iG)*AO(nu,iG)*(ec_z + decdrb*r)
       
        end if

      end do
    end do
  end do

end subroutine UPW92_lda_correlation_potential
