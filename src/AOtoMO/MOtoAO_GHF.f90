subroutine MOtoAO_GHF(nBas2,nBas,S,Ca,Cb,B,A)

! Perform MO to AO transformation of a matrix A for a given metric S
! and coefficients c

  implicit none

! Input variables

  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nBas
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: Ca(nBas,nBas2)
  double precision,intent(in)   :: Cb(nBas,nBas2)
  double precision,intent(in)   :: B(nBas2,nBas2)

! Local variables

  double precision,allocatable  :: SC(:,:)
  double precision,allocatable  :: BSC(:,:)

! Output variables

  double precision,intent(inout):: A(nBas,nBas)

! Memory allocation

  allocate(SC(nBas,nBas2),BSC(nBas2,nBas))

  SC  = matmul(S,Ca)
  BSC = matmul(B,transpose(SC))
  A   = matmul(SC,BSC)

  SC  = matmul(S,Cb)
  BSC = matmul(B,transpose(SC))
  A   = A + matmul(SC,BSc)

end subroutine 
