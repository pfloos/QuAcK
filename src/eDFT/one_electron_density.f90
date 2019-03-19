subroutine one_electron_density(nGrid,nBas,P,AO,dAO,rho,drho)

! Calculate one-electron density

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nGrid
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(3,nBas,nGrid)

! Local variables

  integer                       :: ixyz,iG,mu,nu
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: rho(nGrid) 
  double precision,intent(out)  :: drho(3,nGrid)

  rho(:) = 0d0
  do iG=1,nGrid
    do mu=1,nBas
      do nu=1,nBas
        rho(iG) = rho(iG) + AO(mu,iG)*P(mu,nu)*AO(nu,iG)
      enddo
    enddo
  enddo
 
  drho(:,:) = 0d0
  do ixyz=1,3
    do iG=1,nGrid
      do mu=1,nBas
        do nu=1,nBas
          drho(ixyz,iG) = drho(ixyz,iG) & 
                        + P(mu,nu)*(dAO(ixyz,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(ixyz,nu,iG))
        enddo
      enddo
    enddo
  enddo

end subroutine one_electron_density
