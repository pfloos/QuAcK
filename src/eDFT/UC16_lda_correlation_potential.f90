subroutine UC16_lda_correlation_potential(nGrid,weight,nBas,AO,rho,Fc)

! Compute unrestricted LDA correlation potential

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
  double precision              :: ra,rb,r,rs
  double precision              :: a_p,b_p,ec_p,decdrs_p,decdra_p,decdrb_p
  double precision              :: a_f,b_f,ec_f,decdrs_f,decdra_f,decdrb_f
  double precision              :: ec_z,decdra_z,decdrb_z
  double precision              :: z,dzdra,dzdrb,fz,dfzdz,dfzdra,dfzdrb
  double precision              :: drsdra,drsdrb,dFcdra,dFcdrb

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas,nspin)

! Coefficients for Chachiyo's LDA correlation

  a_p = (log(2d0) - 1d0)/(2d0*pi**2)
  b_p = 20.4562557d0

  a_f = (log(2d0) - 1d0)/(4d0*pi**2)
  b_f = 27.4203609d0

! Compute LDA correlation matrix in the AO basis

  Fc(:,:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

!       Spin-up and spin-down densities

        ra = max(0d0,rho(iG,1))
        rb = max(0d0,rho(iG,2))
      
!       Total density

        r = ra + rb

!       Spin-up part contribution

        if(ra > threshold) then

          rs  = (4d0*pi*r/3d0)**(-1d0/3d0)

          ec_p = a_p*log(1d0 + b_p/rs + b_p/rs**2)
          ec_f = a_f*log(1d0 + b_f/rs + b_f/rs**2)

          z = (ra-rb)/r

          fz = (1d0 + z)**(4d0/3d0) + (1d0 - z)**(4d0/3d0) - 2d0
          fz = fz/(2d0*(2d0**(1d0/3d0) - 1d0))

          ec_z = ec_p + (ec_f - ec_p)*fz

          dzdra  = (1d0 - z)/r
          dfzdz  = (4d0/3d0)*((1d0 + z)**(1d0/3d0) - (1d0 - z)**(1d0/3d0))/(2d0*(2d0**(1d0/3d0) - 1d0))
          dfzdra = dzdra*dfzdz

          drsdra = - (36d0*pi)**(-1d0/3d0)*r**(-4d0/3d0)

          decdrs_p = - a_p/rs**2*(b_p + 2d0*b_p/rs)/(1d0 + b_p/rs + b_p/rs**2)
          decdrs_f = - a_f/rs**2*(b_f + 2d0*b_f/rs)/(1d0 + b_f/rs + b_f/rs**2)

          decdra_p = drsdra*decdrs_p
          decdra_f = drsdra*decdrs_f

          decdra_z = decdra_p + (decdra_f - decdra_p)*fz + (ec_f - ec_p)*dfzdra

          dFcdra = decdra_z*r + ec_z

          Fc(mu,nu,1) = Fc(mu,nu,1) + weight(iG)*AO(mu,iG)*AO(nu,iG)*dFcdra
       
        endif

!       Spin-down part contribution

        if(rb > threshold) then

          rs  = (4d0*pi*r/3d0)**(-1d0/3d0)

          ec_p = a_p*log(1d0 + b_p/rs + b_p/rs**2)
          ec_f = a_f*log(1d0 + b_f/rs + b_f/rs**2)

          z = (ra-rb)/r

          fz = (1d0 + z)**(4d0/3d0) + (1d0 - z)**(4d0/3d0) - 2d0
          fz = fz/(2d0*(2d0**(1d0/3d0) - 1d0))

          ec_z = ec_p + (ec_f - ec_p)*fz

          dzdrb  = - (1d0 + z)/r
          dfzdz  = (4d0/3d0)*((1d0 + z)**(1d0/3d0) - (1d0 - z)**(1d0/3d0))/(2d0*(2d0**(1d0/3d0) - 1d0))
          dfzdrb = dzdrb*dfzdz

          drsdrb = - (36d0*pi)**(-1d0/3d0)*r**(-4d0/3d0)

          decdrs_p = - a_p/rs**2*(b_p + 2d0*b_p/rs)/(1d0 + b_p/rs + b_p/rs**2)
          decdrs_f = - a_f/rs**2*(b_f + 2d0*b_f/rs)/(1d0 + b_f/rs + b_f/rs**2)

          decdrb_p = drsdrb*decdrs_p
          decdrb_f = drsdrb*decdrs_f

          decdrb_z = decdrb_p + (decdrb_f - decdrb_p)*fz + (ec_f - ec_p)*dfzdrb

          dFcdrb = decdrb_z*r + ec_z

          Fc(mu,nu,2) = Fc(mu,nu,2) + weight(iG)*AO(mu,iG)*AO(nu,iG)*dFcdrb

        endif

      enddo
    enddo
  enddo

end subroutine UC16_lda_correlation_potential
