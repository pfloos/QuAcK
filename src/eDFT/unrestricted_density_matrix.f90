subroutine unrestricted_density_matrix(nBas,nEns,nO,c,P)

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
  integer                       :: iEns

! Output variables

  double precision,intent(out)  :: P(nBas,nBas,nspin,nEns)

! N-electron ground state

  iEns = 1
  do ispin=1,nspin
    P(:,:,ispin,iEns) = matmul(c(:,1:nO(ispin),ispin),transpose(c(:,1:nO(ispin),ispin)))
  end do

! (N-1)-electron state: remove spin-down electrons

  iEns = 2
  P(:,:,1,iEns) = matmul(c(:,1:nO(1)  ,1),transpose(c(:,1:nO(1)  ,1)))  
  if (nO(2) > 1) then
    P(:,:,2,iEns) = matmul(c(:,1:nO(2)-1,2),transpose(c(:,1:nO(2)-1,2)))
  else
    P(:,:,2,iEns) = 0.d0
  end if
! (N+1)-electron state: remove spin-up electrons

  iEns = 3 
  P(:,:,1,iEns) = matmul(c(:,1:nO(1)+1,1),transpose(c(:,1:nO(1)+1,1)))  
  P(:,:,2,iEns) = matmul(c(:,1:nO(2)  ,2),transpose(c(:,1:nO(2)  ,2)))

end subroutine unrestricted_density_matrix
