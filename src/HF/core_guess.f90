subroutine core_guess(nBas_AOs, nBas_MOs, Hc, X, c)

!  Core guess of the molecular orbitals for HF calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas_AOs, nBas_MOs
  double precision,intent(in)   :: Hc(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: X(nBas_AOs,nBas_MOs)

! Local variables

  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: e(:)


! Output variables

  double precision,intent(out)  :: c(nBas_AOs,nBas_MOs)

! Memory allocation

  allocate(cp(nBas_MOs,nBas_MOs), e(nBas_MOs))

! Core guess

  cp(:,:) = matmul(transpose(X(:,:)), matmul(Hc(:,:), X(:,:)))

  call diagonalize_matrix(nBas_MOs, cp, e)
  c(:,:) = matmul(X(:,:), cp(:,:))

  deallocate(cp, e)

end subroutine
