subroutine level_shifting(nBas,nO,S,c,F)

! Perform level-shifting on the Fock matrix

  implicit none

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: c(nBas,nBas)

! Local variables

  double precision,allocatable  :: F_MO(:,:)
  double precision,allocatable  :: Sc(:,:)
  double precision,parameter    :: LS = 0.1d0

  integer                       :: a

! Output variables

  double precision,intent(inout):: F(nBas,nBas)

  allocate(F_MO(nBas,nBas),Sc(nBas,nBas))

  F_MO(:,:) = matmul(transpose(c),matmul(F,c))

  do a=nO+1,nBas
    F_MO(a,a) = F_MO(a,a) + LS
  end do

  Sc(:,:) = matmul(S,c)
  F(:,:) = matmul(Sc,matmul(F_MO,transpose(Sc)))

end subroutine level_shifting
