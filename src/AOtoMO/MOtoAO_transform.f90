subroutine MOtoAO_transform(nBas,S,c,B,A)

! Perform MO to AO transformation of a matrix A for a given metric S
! and coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: c(nBas,nBas)
  double precision,intent(in)   :: B(nBas,nBas)

! Local variables

  double precision,allocatable  :: Sc(:,:)

! Output variables

  double precision,intent(out)  :: A(nBas,nBas)

! Memory allocation

  allocate(Sc(nBas,nBas))

  Sc = matmul(S,c)
  A = matmul(Sc,matmul(B,transpose(Sc)))

end subroutine 
