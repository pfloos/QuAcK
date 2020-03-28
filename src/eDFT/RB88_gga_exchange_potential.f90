subroutine RB88_gga_exchange_potential(nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

! Compute restricted Becke's GGA exchange potential

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)
  double precision,intent(in)   :: drho(ncart,nGrid)

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: alpha
  double precision              :: beta
  double precision              :: r,g,vAO,gAO

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Coefficients for B88 GGA exchange functional

  alpha = -(3d0/2d0)*(3d0/(4d0*pi))**(1d0/3d0)
  beta  = 0.0042d0

! Compute GGA exchange matrix in the AO basis

  Fx(:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        r = max(0d0,0.5d0*rho(iG))
 
        if(r > threshold) then

          g = drho(1,iG)**2 + drho(2,iG)**2 + drho(3,iG)**2
          vAO = weight(iG)*AO(mu,iG)*AO(nu,iG)
          Fx(mu,nu) = Fx(mu,nu) &
                    + vAO*(4d0/3d0*r**(1d0/3d0)*(alpha - beta*g**(3d0/4d0)/r**2) &
                         + 2d0*beta*g**(3d0/4d0)/r**(5d0/3d0))
          
          gAO = drho(1,iG)*(dAO(1,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(1,nu,iG)) & 
              + drho(2,iG)*(dAO(2,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(2,nu,iG)) & 
              + drho(3,iG)*(dAO(3,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(3,nu,iG))

          gAO = weight(iG)*gAO
          
          Fx(mu,nu) = Fx(mu,nu) - 2d0*gAO*3d0/4d0*beta*g**(-1d0/4d0)/r**(2d0/3d0)

        end if

      end do
    end do
  end do

end subroutine RB88_gga_exchange_potential
