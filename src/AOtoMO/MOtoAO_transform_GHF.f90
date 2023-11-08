subroutine MOtoAO_transform_GHF(nBas2,nBas,S,Ca,Cb,B,A)

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

  double precision,allocatable  :: Sc(:,:)

! Output variables

  double precision,intent(inout):: A(nBas,nBas)

! Memory allocation

  allocate(Sc(nBas,nBas2))

  Sc = matmul(S,Ca)
  A = matmul(Sc,matmul(B,transpose(Sc)))

  Sc = matmul(S,Cb)
  A = A + matmul(Sc,matmul(B,transpose(Sc)))

end subroutine 
