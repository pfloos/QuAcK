subroutine density(nGrid,nBas,P,AO,rho)

! Calculate one-electron density

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: AO(nBas,nGrid)

! Local variables

  integer                       :: iG,mu,nu

! Output variables

  double precision,intent(out)  :: rho(nGrid) 

  rho(:) = 0d0

  do iG=1,nGrid
    do mu=1,nBas
      do nu=1,nBas

        rho(iG) = rho(iG) + AO(mu,iG)*P(mu,nu)*AO(nu,iG)

      enddo
    enddo
  enddo

end subroutine 
