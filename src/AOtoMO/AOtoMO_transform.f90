subroutine AOtoMO_transform(nBas,c,A,B)

! Perform AO to MO transformation of a matrix A for given coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: A(nBas,nBas)

! Output variables

  double precision,intent(out)  :: B(nBas,nBas)

  B = matmul(transpose(c),matmul(A,c))

end subroutine 
