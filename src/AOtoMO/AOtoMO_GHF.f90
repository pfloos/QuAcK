subroutine AOtoMO_GHF(nBas,nBas2,Ca,Cb,A,B)

! Perform AO to MO transformation of a matrix A for given coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  double precision,intent(in)   :: Ca(nBas,nBas2)
  double precision,intent(in)   :: Cb(nBas,nBas2)
  double precision,intent(in)   :: A(nBas,nBas)

! Local variables

  double precision,allocatable  :: AC(:,:)

! Output variables

  double precision,intent(out)  :: B(nBas2,nBas2)

  allocate(AC(nBas,nBas2))

  AC = matmul(A,Ca)
  B  = matmul(transpose(Ca),AC)

  AC = matmul(A,Cb)
  B  = B + matmul(transpose(Cb),AC)

end subroutine 
