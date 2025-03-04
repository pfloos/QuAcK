subroutine complex_core_guess(nBas, nOrb, Hc, X, c)

!  Core guess of the molecular orbitals for HF calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas, nOrb
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)

! Local variables

  complex*16,allocatable  :: cp(:,:)
  complex*16,allocatable  :: e(:)


! Output variables

  complex*16,intent(out)  :: c(nBas,nOrb)

! Memory allocation

  allocate(cp(nOrb,nOrb), e(nOrb))

! Core guess

  cp(:,:) = matmul(transpose(X(:,:)), matmul(Hc(:,:), X(:,:)))

  call diagonalize_matrix(nOrb, cp, e)
  c(:,:) = matmul(X(:,:), cp(:,:))

  deallocate(cp, e)

end subroutine
