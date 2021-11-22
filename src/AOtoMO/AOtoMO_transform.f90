subroutine AOtoMO_transform(nBas,c,A)

! Perform AO to MO transformation of a matrix A for given coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)

  integer :: i,j

! Output variables

  double precision,intent(inout):: A(nBas,nBas)

  A = matmul(transpose(c),matmul(A,c))

! do j=1,nBas
!   do i=1,nBas
!     write(10,'(I5,I5,F16.10)') i,j,A(i,j)
!   enddo
! enddo


end subroutine AOtoMO_transform
