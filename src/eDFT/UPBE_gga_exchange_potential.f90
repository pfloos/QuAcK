subroutine UPBE_gga_exchange_potential(nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

! Compute PBE GGA exchange potential

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(3,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(3,nGrid)

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: mupbe,kappa
  double precision              :: r,g,s2,vAO,gAO

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Coefficients for PBE exchange functional

  mupbe = ((1d0/2d0**(1d0/3d0))/(2d0*(3d0*pi**2)**(1d0/3d0)))**2*0.21951d0
  kappa = 0.804d0

! Compute GGA exchange matrix in the AO basis

  Fx(:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        r = max(0d0,rho(iG))
 
        if(r > threshold) then

          g = drho(1,iG)**2 + drho(2,iG)**2 + drho(3,iG)**2
          s2 = g/r**(8d0/3d0)

          vAO = weight(iG)*AO(mu,iG)*AO(nu,iG)

          Fx(mu,nu) = Fx(mu,nu) &
                    + vAO*4d0/3d0*CxLSDA*r**(1d0/3d0)*(1d0 + kappa - kappa/(1d0 + mupbe*s2/kappa)) &
                    - vAO*8d0/3d0*CxLSDA*r**(1d0/3d0)*mupbe*s2/(1d0 + mupbe*s2/kappa)**2
          
          gAO = drho(1,iG)*(dAO(1,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(1,nu,iG)) & 
              + drho(2,iG)*(dAO(2,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(2,nu,iG)) & 
              + drho(3,iG)*(dAO(3,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(3,nu,iG))
          gAO = weight(iG)*gAO
          
          Fx(mu,nu) = Fx(mu,nu) + 2d0*gAO*CxLSDA*r**(-4d0/3d0)*mupbe/(1d0 + mupbe*s2/kappa)**2

        end if

      end do
    end do
  end do

end subroutine UPBE_gga_exchange_potential
