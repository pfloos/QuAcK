subroutine AOtoMO_transform(nBas,c,A)

! Perform AO to MO transformation of a matrix A for given coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)

! Output variables

  double precision,intent(inout):: A(nBas,nBas)

  A = matmul(transpose(c),matmul(A,c))

end subroutine AOtoMO_transform
