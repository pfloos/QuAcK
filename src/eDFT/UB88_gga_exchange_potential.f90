subroutine UB88_gga_exchange_potential(nGrid,weight,nBas,AO,dAO,rho,drho,Fx)

! Compute Becke's GGA exchange potential

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
  double precision              :: b
  double precision              :: vAO,gAO
  double precision              :: r,g,x,dxdr,dxdg,f

! Output variables

  double precision,intent(out)  :: Fx(nBas,nBas)

! Coefficients for B88 GGA exchange functional

  b = 0.0042d0

! Compute GGA exchange matrix in the AO basis

  Fx(:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        r = max(0d0,rho(iG))
 
        if(r > threshold) then

          vAO = weight(iG)*AO(mu,iG)*AO(nu,iG)

          g = drho(1,iG)**2 + drho(2,iG)**2 + drho(3,iG)**2
          x = sqrt(g)/r**(4d0/3d0)
          dxdr = - 4d0*sqrt(g)/(3d0*r**(7d0/3d0))/x
          dxdg = + 1d0/(2d0*sqrt(g)*r**(4d0/3d0))/x

          f = b*x**2/(1d0 + 6d0*b*x*asinh(x))

          Fx(mu,nu) = Fx(mu,nu) + vAO*(                 &
                      4d0/3d0*r**(1d0/3d0)*(CxLSDA - f) &
                    - 2d0*r**(4d0/3d0)*dxdr*f           &
                    + r**(4d0/3d0)*dxdr*(6d0*b*x*asinh(x) + 6d0*b*x**2/sqrt(1d0+x**2))*f/(1d0 + 6d0*b*x*asinh(x)) )
          
          gAO = drho(1,iG)*(dAO(1,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(1,nu,iG)) & 
              + drho(2,iG)*(dAO(2,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(2,nu,iG)) & 
              + drho(3,iG)*(dAO(3,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(3,nu,iG))
          gAO = weight(iG)*gAO
          
          Fx(mu,nu) = Fx(mu,nu) + 2d0*gAO*r**(4d0/3d0)*dxdg*(   &
                    - 2d0*f + (6d0*b*x*asinh(x) + 6d0*b*x**2/sqrt(1d0+x**2))*f/(1d0 + 6d0*b*x*asinh(x)) )

        end if

      end do
    end do
  end do

end subroutine UB88_gga_exchange_potential
