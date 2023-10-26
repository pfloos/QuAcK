subroutine AOtoMO_transform_GHF(nBas,nBas2,Ca,Cb,A,B)

! Perform AO to MO transformation of a matrix A for given coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  double precision,intent(in)   :: Ca(nBas,nBas2)
  double precision,intent(in)   :: Cb(nBas,nBas2)
  double precision,intent(inout):: A(nBas,nBas)

! Output variables

  double precision,intent(inout):: B(nBas2,nBas2)

  B = matmul(transpose(Ca),matmul(A,Ca)) &
    + matmul(transpose(Cb),matmul(A,Cb))

end subroutine 
