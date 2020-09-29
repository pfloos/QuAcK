subroutine overlap(nBas,bra,ket)

! Compute the overlap between two sets of coefficients

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: bra(nBas,nBas),ket(nBas,nBas)

! Local variables

  double precision,allocatable  :: s(:),Ov(:,:)

! Allocate

  allocate(s(nBas),Ov(nBas,nBas))

! Compute overlap

  Ov = matmul(transpose(bra),ket)

  call diagonalize_matrix(nBas,Ov,s)

! Print results

  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' Overlap '
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,nBas,Ov)
  write(*,*)
  write(*,'(A50)') '---------------------------------------'
  write(*,'(A50)') ' Eigenvalues of overlap matrix'
  write(*,'(A50)') '---------------------------------------'
  call matout(nBas,1,s)
  write(*,*)

end subroutine overlap
