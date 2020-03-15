subroutine restricted_density_matrix(nBas,nEns,nO,c,P)

! Calculate density matrices

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  integer,intent(in)            :: nO
  double precision,intent(in)   :: c(nBas,nBas)

! Local variables

  integer                       :: iEns

! Output variables

  double precision,intent(out)  :: P(nBas,nBas,nEns)

! Ground state density matrix

  iEns = 1
  P(:,:,iEns) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

! Doubly-excited state density matrix

  iEns = 2 
  P(:,:,iEns) = 2d0*matmul(c(:,1:nO-1),transpose(c(:,1:nO-1))) &
              + 2d0*matmul(c(:,nO+1:nO+1),transpose(c(:,nO+1:nO+1)))

end subroutine restricted_density_matrix
