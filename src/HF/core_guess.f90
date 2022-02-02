subroutine core_guess(nBas,Hc,X,c)

!  Core guess of the molecular orbitals for HF calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)

! Local variables

  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: e(:)


! Output variables

  double precision,intent(out)  :: c(nBas,nBas)

! Memory allocation

  allocate(cp(nBas,nBas),e(nBas))

! Core guess

  cp(:,:) = matmul(transpose(X(:,:)),matmul(Hc(:,:),X(:,:)))
  call diagonalize_matrix(nBas,cp,e)
  c(:,:) = matmul(X(:,:),cp(:,:))

end subroutine
