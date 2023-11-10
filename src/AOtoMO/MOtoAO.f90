subroutine MOtoAO(nBas,S,C,B,A)

! Perform MO to AO transformation of a matrix A for a given metric S
! and coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: C(nBas,nBas)
  double precision,intent(in)   :: B(nBas,nBas)

! Local variables

  double precision,allocatable  :: SC(:,:),BSC(:,:)

! Output variables

  double precision,intent(out)  :: A(nBas,nBas)

! Memory allocation

  allocate(SC(nBas,nBas),BSC(nBas,nBas))

  SC  = matmul(S,C)
  BSC = matmul(B,transpose(SC))
  A   = matmul(SC,BSC)

end subroutine 
