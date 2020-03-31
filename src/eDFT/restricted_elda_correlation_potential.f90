subroutine restricted_elda_correlation_potential(aMFL,nGrid,weight,nBas,AO,rho,Fc)

! Compute LDA correlation energy of 2-glomium for various states

  implicit none
  include 'parameters.h'

! Input variables

  double precision,intent(in)   :: aMFL(3)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  integer                       :: mu,nu,iG
  double precision              :: r,ec_p
  double precision              :: dFcdr

! Output variables

  double precision,intent(out)  :: Fc(nBas,nBas)

! Compute eLDA correlation potential

  Fc(:,:) = 0d0

  do mu=1,nBas
    do nu=1,nBas
      do iG=1,nGrid

        r = max(0d0,rho(iG))

        if(r > threshold) then

          ec_p  = aMFL(1)/(1d0 + aMFL(2)*r**(-1d0/6d0) + aMFL(3)*r**(-1d0/3d0))
          dFcdr = aMFL(2)*r**(-1d0/6d0) + 2d0*aMFL(3)*r**(-1d0/3d0)
          dFcdr = dFcdr/(1d0 + aMFL(2)*r**(-1d0/6d0) + aMFL(3)*r**(-1d0/3d0))
          dFcdr = ec_p*dFcdr/(6d0*r)
          dFcdr = ec_p + dFcdr*r

          Fc(mu,nu) = Fc(mu,nu) + weight(iG)*AO(mu,iG)*AO(nu,iG)*dFcdr

         endif

      end do
    end do
  end do

end subroutine restricted_elda_correlation_potential
