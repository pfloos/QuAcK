subroutine level_shifting(level_shift, nBas_AOs, nBas_MOs, nO, S, c, F)

! Perform level-shifting on the Fock matrix

  implicit none

! Input variables

  double precision,intent(in)   :: level_shift
  integer,intent(in)            :: nBas_AOs, nBas_MOs
  integer,intent(in)            :: nO
  double precision,intent(in)   :: S(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: c(nBas_AOs,nBas_MOs)

! Local variables

  double precision,allocatable  :: F_MO(:,:)
  double precision,allocatable  :: Sc(:,:)

  integer                       :: a

! Output variables

  double precision,intent(inout):: F(nBas_AOs,nBas_AOs)

  allocate(F_MO(nBas_MOs,nBas_MOs), Sc(nBas_AOs,nBas_MOs))

  F_MO(:,:) = matmul(transpose(c), matmul(F, c))

  do a = nO+1, nBas_MOs
    F_MO(a,a) = F_MO(a,a) + level_shift
  end do

  Sc(:,:) = matmul(S, c)
  F(:,:) = matmul(Sc, matmul(F_MO, transpose(Sc)))

  deallocate(F_MO, Sc)

end subroutine 
