subroutine MOtoAO_transform(nBas,S,c,A)

! Perform MO to AO transformation of a matrix A for a given metric S
! and coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: S(nBas,nBas),c(nBas,nBas)

! Local variables

  double precision,allocatable  :: Sc(:,:)

! Output variables

  double precision,intent(inout):: A(nBas,nBas)

! Memory allocation
  allocate(Sc(nBas,nBas))

  Sc = matmul(S,c)
  A = matmul(Sc,matmul(A,transpose(Sc)))

end subroutine 
