subroutine core_guess(nBas, nOrb, Hc, X, c)

!  Core guess of the molecular orbitals for HF calculation

  implicit none

! Input variables

  integer,intent(in)            :: nBas, nOrb
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)

! Local variables

  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: e(:)


! Output variables

  double precision,intent(out)  :: c(nBas,nOrb)

! Memory allocation

  allocate(cp(nOrb,nOrb), e(nOrb))

! Core guess

  cp(:,:) = matmul(transpose(X(:,:)), matmul(Hc(:,:), X(:,:)))

  call diagonalize_matrix(nOrb, cp, e)
  write(*,*) 'cp'
  call matout(nOrb,nOrb,cp)
  write(*,*) 'Eigenvalues e'
  call vecout(nOrb,e)
  c(:,:) = matmul(X(:,:), cp(:,:))

  deallocate(cp, e)

end subroutine
