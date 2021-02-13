subroutine ULYP_gga_correlation_potential(nGrid,weight,nBas,AO,dAO,rho,drho,Fc)

! Compute LYP correlation potential

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)
  double precision,intent(in)   :: drho(ncart,nGrid,nspin)

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: vAO,gaAO,gbAO
  double precision              :: ra,rb,r
  double precision              :: ga,gab,gb,g
  double precision              :: dfdra,dfdrb
  double precision              :: dfdga,dfdgab,dfdgb
  double precision              :: dodra,dodrb,dddra,dddrb

  double precision              :: a,b,c,d
  double precision              :: Cf,omega,delta

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas,nspin)

! Prameter of the functional

  a = 0.04918d0
  b = 0.132d0
  c = 0.2533d0
  d = 0.349d0

  Cf = 3d0/10d0*(3d0*pi**2)**(2d0/3d0)

! Compute matrix elements in the AO basis

  Fc(:,:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        ra = max(0d0,rho(iG,1))
        rb = max(0d0,rho(iG,2))
        r  = ra + rb
 
        if(r > threshold) then

          ga  = drho(1,iG,1)*drho(1,iG,1) + drho(2,iG,1)*drho(2,iG,1) + drho(3,iG,1)*drho(3,iG,1)
          gb  = drho(1,iG,2)*drho(1,iG,2) + drho(2,iG,2)*drho(2,iG,2) + drho(3,iG,2)*drho(3,iG,2)
          gab = drho(1,iG,1)*drho(1,iG,2) + drho(2,iG,1)*drho(2,iG,2) + drho(3,iG,1)*drho(3,iG,2)
          g   = ga + 2d0*gab + gb

          omega = exp(-c*r**(-1d0/3d0))/(1d0 + d*r**(-1d0/3d0))*r**(-11d0/3d0)
          delta = c*r**(-1d0/3d0) + d*r**(-1d0/3d0)/(1d0 + d*r**(-1d0/3d0))

          vAO = weight(iG)*AO(mu,iG)*AO(nu,iG)

          dodra = (d/(3d0*r**(4d0/3d0)*(1d0 + d*r**(-1d0/3d0))) + c/(3d0*r**(4d0/3d0)) - 11d0/(3d0*r))*omega
          dodrb = dodra
         
          dddra = - c/3d0*r**(-4d0/3d0)                                 & 
                  + d**2/(3d0*(1d0 + d*r**(-1d0/3d0))**2)*r**(-5d0/3d0) &
                  - d/(3d0*(1d0 + d*r**(-1d0/3d0)))*r**(-4d0/3d0)
          dddrb = dddra

          dfdra = - 4d0*a/(1d0 + d*r**(-1d0/3d0))*rb/r                        & 
                  - 4d0/3d0*a*d/(1d0 + d*r**(-1d0/3d0))**2*ra*rb/r**(7d0/3d0) &
                  + 4d0*a/(1d0 + d*r**(-1d0/3d0))*ra*rb/r**2                  &
                  - a*b*omega*rb*(                                            &
                    + 2d0**(11d0/3d0)*Cf*(ra**(8d0/3d0) + rb**(8d0/3d0))      &
                    + (47d0/18d0 - 7d0*delta/18d0)*g                          &
                    - (5d0/2d0 - delta/18d0)*(ga + gb)                        &
                    - (delta - 11d0)/9d0*(ra/r*ga + rb/r*gb)                  &
                    - 4d0/3d0*r/rb*g                                          &
                    + (4d0/3d0*r/rb - 2d0*ra/rb)*gb                           &
                    + 4d0/3d0*r/rb*ga )                                       &
                  - a*b*omega*ra*rb*(                                         &
                    + 8d0/3d0*2d0**(11d0/3d0)*Cf*ra**(5d0/3d0)                &
                    - 7d0*dddra/18d0*g                                        &
                    + dddra/18d0*(ga + gb)                                    &
                    - dddra/9d0*(ra/r*ga + rb/r*gb)                           &
                    - (delta - 11d0)/(9d0*r)*(-ra/r*ga - rb/r*gb + ga) )      &
                  - a*b*dodra*ra*rb*(                                         &
                    + 2d0**(11d0/3d0)*Cf*(ra**(8d0/3d0) + rb**(8d0/3d0))      &
                    + (47d0/18d0 - 7d0*delta/18d0)*g                          &
                    - (5d0/2d0 - delta/18d0)*(ga + gb)                        &
                    - (delta - 11d0)/9d0*(ra/r*ga + rb/r*gb) )                &
                  - a*b*dodra*(                                               &
                    - 2d0*r**2/3d0*g                                          &
                    + (2d0*r**2/3d0 - ra**2)*gb                               &
                    + (2d0*r**2/3d0 - rb**2)*ga )

          dfdrb = - 4d0*a/(1d0 + d*r**(-1d0/3d0))*ra/r                        & 
                  - 4d0/3d0*a*d/(1d0 + d*r**(-1d0/3d0))**2*ra*rb/r**(7d0/3d0) &
                  + 4d0*a/(1d0 + d*r**(-1d0/3d0))*ra*rb/r**2                  &
                  - a*b*omega*ra*(                                            &
                    + 2d0**(11d0/3d0)*Cf*(ra**(8d0/3d0) + rb**(8d0/3d0))      &
                    + (47d0/18d0 - 7d0*delta/18d0)*g                          &
                    - (5d0/2d0 - delta/18d0)*(ga + gb)                        &
                    - (delta - 11d0)/9d0*(ra/r*ga + rb/r*gb)                  &
                    - 4d0/3d0*r/ra*g                                          &
                    + (4d0/3d0*r/ra - 2d0*rb/ra)*ga                           &
                    + 4d0/3d0*r/ra*gb )                                       &
                  - a*b*omega*ra*rb*(                                         &
                    + 8d0/3d0*2d0**(11d0/3d0)*Cf*rb**(5d0/3d0)                &
                    - 7d0*dddrb/18d0*g                                        &
                    + dddrb/18d0*(ga + gb)                                    &
                    - dddrb/9d0*(ra/r*ga + rb/r*gb)                           &
                    - (delta - 11d0)/(9d0*r)*(-ra/r*ga - rb/r*gb + gb) )      &
                  - a*b*dodrb*ra*rb*(                                         &
                    + 2d0**(11d0/3d0)*Cf*(ra**(8d0/3d0) + rb**(8d0/3d0))      &
                    + (47d0/18d0 - 7d0*delta/18d0)*g                          &
                    - (5d0/2d0 - delta/18d0)*(ga + gb)                        &
                    - (delta - 11d0)/9d0*(ra/r*ga + rb/r*gb) )                &
                  - a*b*dodrb*(                                               &
                    - 2d0*r**2/3d0*g                                          &
                    + (2d0*r**2/3d0 - ra**2)*gb                               &
                    + (2d0*r**2/3d0 - rb**2)*ga )

          Fc(mu,nu,1) = Fc(mu,nu,1) + vAO*dfdra
          Fc(mu,nu,2) = Fc(mu,nu,2) + vAO*dfdrb
          
          gaAO = drho(1,iG,1)*(dAO(1,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(1,nu,iG)) & 
               + drho(2,iG,1)*(dAO(2,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(2,nu,iG)) & 
               + drho(3,iG,1)*(dAO(3,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(3,nu,iG))
          gaAO = weight(iG)*gaAO

          gbAO = drho(1,iG,2)*(dAO(1,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(1,nu,iG)) & 
               + drho(2,iG,2)*(dAO(2,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(2,nu,iG)) & 
               + drho(3,iG,2)*(dAO(3,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(3,nu,iG))
          gbAO = weight(iG)*gbAO

          dfdga  = -a*b*omega*(-rb**2 + ra*rb*(1d0/9d0 - (delta-11d0)/9d0*ra/r - delta/3d0))
          dfdgab = -a*b*omega*(-4d0/3d0*r**2 + 2d0*ra*rb*(47d0/18d0 - 7d0*delta/18d0))
          dfdgb  = -a*b*omega*(-ra**2 + ra*rb*(1d0/9d0 - (delta-11d0)/9d0*rb/r - delta/3d0))
          
          Fc(mu,nu,1) = Fc(mu,nu,1) + 2d0*gaAO*dfdga + gbAO*dfdgab
          Fc(mu,nu,2) = Fc(mu,nu,2) + 2d0*gbAO*dfdgb + gaAO*dfdgab

        end if

      end do
    end do
  end do

end subroutine ULYP_gga_correlation_potential
