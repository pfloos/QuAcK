subroutine complex_MOM_density_matrix(nBas, nOrb, nO, S, c, cGuess, occupations, occupationsGuess, P)

! Calculate Overlap projected to occupied space of guess for one spin component.

  implicit none
  include 'parameters.h'

! Input variables
  integer, intent(in)             :: nBas, nOrb, nO
  double precision,intent(in)     :: S(nBas,nBas)
  complex*16, intent(in)          :: c(nBas, nOrb), cGuess(nBas,nOrb)
  integer,intent(in)              :: occupationsGuess(nO)

! Output variables
  integer, intent(inout)          :: occupations(nO)
  complex*16, intent(out)         :: P(nBas,nBas)

! Local variables
  integer                         :: i
  double precision,allocatable    :: projO(:) 

  allocate(projO(nOrb))
  
  call complex_projected_overlap(nBas, nOrb,nO, S, c, cGuess, occupationsGuess, projO)
  ! Select orbitals with maximum overlap
  call MOM_idx(nO,nOrb,projO,occupations(1:nO))
  ! Compute density with selected orbitals
  P(:,:) = matmul(c(:,occupations(1:nO)),&
        transpose(c(:,occupations(1:nO))))

  deallocate(projO)
end subroutine
