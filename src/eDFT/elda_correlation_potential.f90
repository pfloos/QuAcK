subroutine elda_correlation_potential(nEns,aLF,nGrid,weight,nBas,AO,rho,Fc)

! Compute LDA correlation energy of 2-glomium for various states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: aLF(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: ra,rb,r,ec_p
  double precision              :: dFcdra,dFcdrb

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas,nspin)

! Compute eLDA correlation potential

  Fc(:,:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        ra = max(0d0,rho(iG,1))
        rb = max(0d0,rho(iG,2))

        if(ra > threshold) then

          r = ra + rb

          ec_p = aLF(1)/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
          dFcdra = aLF(2)*r**(-1d0/6d0) + 2d0*aLF(3)*r**(-1d0/3d0)
          dFcdra = dFcdra/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
          dFcdra = ec_p*dFcdra/(6d0*r)
          dFcdra = ec_p + dFcdra*r

          Fc(mu,nu,1) = Fc(mu,nu,1) + weight(iG)*AO(mu,iG)*AO(nu,iG)*dFcdra

         endif

        if(rb > threshold) then

          r = ra + rb

          ec_p = aLF(1)/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
          dFcdrb = aLF(2)*r**(-1d0/6d0) + 2d0*aLF(3)*r**(-1d0/3d0)
          dFcdrb = dFcdrb/(1d0 + aLF(2)*r**(-1d0/6d0) + aLF(3)*r**(-1d0/3d0))
          dFcdrb = ec_p*dFcdrb/(6d0*r)
          dFcdrb = ec_p + dFcdrb*r

          Fc(mu,nu,2) = Fc(mu,nu,2) + weight(iG)*AO(mu,iG)*AO(nu,iG)*dFcdrb

         endif

      end do
    end do
  end do

end subroutine elda_correlation_potential
