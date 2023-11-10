subroutine AOtoMO(nBas,C,A,B)

! Perform AO to MO transformation of a matrix A for given coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: C(nBas,nBas)
  double precision,intent(in)   :: A(nBas,nBas)

! Local variables

  double precision,allocatable  :: AC(:,:)

! Output variables

  double precision,intent(out)  :: B(nBas,nBas)

  allocate(AC(nBas,nBas))

  AC = matmul(A,C)
  B  = matmul(transpose(C),AC)

end subroutine 
