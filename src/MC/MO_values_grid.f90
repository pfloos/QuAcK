subroutine MO_values_grid(nBas,nGrid,c,AO,dAO,MO,dMO)

! Compute values of the MOs and their derivatives with respect to the cartesian coordinates

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: AO(nBas,nGrid)
  double precision,intent(in)   :: dAO(ncart,nBas,nGrid)

! Local variables

  integer                       :: p,mu,iG

! Output variables

  double precision,intent(out)  :: MO(nBas,nGrid)
  double precision,intent(out)  :: dMO(ncart,nBas,nGrid)

! Initialization

  MO(:,:)    = 0d0
  dMO(:,:,:) = 0d0

  do p=1,nBas
    do mu=1,nBas
      do iG=1,ngrid

        MO(p,iG) = MO(p,iG) + c(mu,p)*AO(mu,iG)

        dMO(1,p,iG) = dMO(1,p,iG) + c(mu,p)*dAO(1,mu,iG)
        dMO(2,p,iG) = dMO(2,p,iG) + c(mu,p)*dAO(2,mu,iG)
        dMO(3,p,iG) = dMO(3,p,iG) + c(mu,p)*dAO(3,mu,iG)

      end do
    end do
  end do

end subroutine MO_values_grid
