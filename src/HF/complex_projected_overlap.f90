subroutine complex_projected_overlap(nBas, nOrb, nO, S, c, cGuess, occupationsGuess, projO)

! Calculate Overlap projected to occupied space of guess

  implicit none
  include 'parameters.h'

! Input variables
  integer, intent(in)             :: nBas, nOrb, nO
  double precision,intent(in)     :: S(nBas,nBas)
  complex*16, intent(in)          :: c(nBas, nOrb), cGuess(nBas,nOrb)
  integer,intent(in)              :: occupationsGuess(nO)

! Output variables
  double precision, intent(out) :: projO(nOrb) 

! Local variables
  integer                         :: i
  complex*16,allocatable          :: O(:,:)

  allocate(O(nOrb,nOrb))

  projO = 0d0
  O = matmul(matmul(conjg(transpose(cGuess)),S),c)
  do i=1,nO
    projO(:) = projO(:) + real(O(occupationsGuess(i),:))**2 + aimag(O(occupationsGuess(i),:))**2
  end do
  
  deallocate(O)
end subroutine
