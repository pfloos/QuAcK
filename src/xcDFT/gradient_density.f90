subroutine gradient_density(nGrid,nBas,P,AO,dAO,drho)

! Calculate gradient of the one-electron density

  implicit none
  include 'parameters.h'

! Input variables

  double precision,parameter    :: thresh = 1d-15

  integer,intent(in)            :: nGrid
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: P(nBas,nBas)
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(3,nBas,nGrid)

! Local variables

  integer                       :: ixyz,iG,mu,nu
  double precision,external     :: trace_matrix

! Output variables

  double precision,intent(out)  :: drho(3,nGrid)

  drho(:,:) = 0d0
  do iG=1,nGrid
    do mu=1,nBas
      do nu=1,nBas
        do ixyz=1,3
          drho(ixyz,iG) = drho(ixyz,iG) & 
                        + P(mu,nu)*(dAO(ixyz,mu,iG)*AO(nu,iG) + AO(mu,iG)*dAO(ixyz,nu,iG))
        enddo
      enddo
    enddo
  enddo

  do iG=1,nGrid
    do ixyz=1,3
      if(abs(drho(ixyz,iG)) < thresh) drho(ixyz,iG) = thresh
    enddo
  enddo

end subroutine gradient_density
