subroutine level_shifting(level_shift, nBas, nOrb, nO, S, c, F)

! Perform level-shifting on the Fock matrix

  implicit none

! Input variables

  double precision,intent(in)   :: level_shift
  integer,intent(in)            :: nBas, nOrb
  integer,intent(in)            :: nO
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: c(nBas,nOrb)

! Local variables

  double precision,allocatable  :: F_MO(:,:)
  double precision,allocatable  :: Sc(:,:)

  integer                       :: a

! Output variables

  double precision,intent(inout):: F(nBas,nBas)

  allocate(F_MO(nOrb,nOrb), Sc(nBas,nOrb))

  F_MO(:,:) = matmul(transpose(c), matmul(F, c))

  do a = nO+1, nOrb
    F_MO(a,a) = F_MO(a,a) + level_shift
  end do

  Sc(:,:) = matmul(S, c)
  F(:,:) = matmul(Sc, matmul(F_MO, transpose(Sc)))

  deallocate(F_MO, Sc)

end subroutine 
