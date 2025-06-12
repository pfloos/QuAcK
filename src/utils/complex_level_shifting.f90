subroutine complex_level_shifting(level_shift, nBas, nOrb, nO, S, c, F)

! Perform level-shifting on the Fock matrix

  implicit none

! Input variables

  double precision,intent(in)   :: level_shift
  integer,intent(in)            :: nBas, nOrb
  integer,intent(in)            :: nO
  double precision,intent(in)   :: S(nBas,nBas)
  complex*16,intent(in)         :: c(nBas,nOrb)

! Local variables

  complex*16,allocatable        :: F_MO(:,:)
  complex*16,allocatable        :: Sc(:,:)
  complex*16                    :: complex_level_shift
  integer                       :: a

! Output variables

  complex*16,intent(inout)      :: F(nBas,nBas)

  allocate(F_MO(nOrb,nOrb), Sc(nBas,nOrb))
  complex_level_shift = cmplx(level_shift, 0.0,kind=8)
  F_MO(:,:) = matmul(transpose(c), matmul(F, c))

  do a = nO+1, nOrb
    F_MO(a,a) = F_MO(a,a) + level_shift
  end do

  Sc(:,:) = matmul(S, c)
  F(:,:) = matmul(Sc, matmul(F_MO, transpose(Sc)))

  deallocate(F_MO, Sc)

end subroutine 
