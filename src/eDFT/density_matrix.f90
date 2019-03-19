subroutine density_matrix(nBas,nEns,nO,c,P)

! Calculate density matrices

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nEns
  integer,intent(in)            :: nO(nspin)
  double precision,intent(in)   :: c(nBas,nBas,nspin)

! Local variables

  integer                       :: ispin

! Output variables

  double precision,intent(out)  :: P(nBas,nBas,nspin,nEns)

! Ground state density matrix

  do ispin=1,nspin
    P(:,:,ispin,1) = matmul(c(:,1:nO(ispin),ispin),transpose(c(:,1:nO(ispin),ispin)))
  end do

! Singly-excited state density matrix

  P(:,:,1,2) = matmul(c(:,1:nO(1)-1,1),transpose(c(:,1:nO(1)-1,1))) &
             + matmul(c(:,nO(1)+1:nO(1)+1,1),transpose(c(:,nO(1)+1:nO(1)+1,1)))
  P(:,:,2,2) = matmul(c(:,1:nO(2),2),transpose(c(:,1:nO(2),2)))

! Doubly-excited state density matrix

  do ispin=1,nspin
    P(:,:,ispin,3) = matmul(c(:,1:nO(ispin)-1,ispin),transpose(c(:,1:nO(ispin)-1,ispin))) &
                   + matmul(c(:,nO(ispin)+1:nO(ispin)+1,ispin),transpose(c(:,nO(ispin)+1:nO(ispin)+1,ispin)))
  end do

end subroutine density_matrix
