subroutine natural_orbital(nBas,nO,cHF,c)

! Compute natural orbitals and natural occupancies

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nO
  double precision,intent(in)   :: cHF(nBas,nBas),c(nBas,nBas)

! Local variables

  integer :: i,j,k
  double precision,allocatable  :: eNO(:),cNO(:,:),P(:,:)

! Allocate

  allocate(eNO(nBas),cNO(nBas,nBas),P(nBas,nBas))

! Compute density matrix

  P = matmul(transpose(cHF),cHF)

  call matout(nBas,nBas,P)

  cNO = 0d0

  do i=1,nBas
    do j=1,nBas
      do k=1,1
        cNO(i,j) = cNO(i,j) + 2d0*P(i,k)*P(j,k)
      enddo
    enddo
  enddo

! cNO(:,:) = matmul(c(:,1:nO),transpose(c(:,1:nO)))

! cNO = matmul(transpose(cHF),matmul(cNO,cHF))

  call diagonalize_matrix(nBas,cNO,eNO)

! Print results

  write(*,'(A50)') '---------------------------------------'
  write(*,'(A32)') ' Natural orbitals '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,nBas,cNO)
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A32)') ' Natural occupancies'
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,eNO)
  write(*,*)

end subroutine natural_orbital
