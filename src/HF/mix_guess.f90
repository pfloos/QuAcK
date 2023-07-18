subroutine mix_guess(nBas,nO,c)

!  Guess mixing for UHF calculations on singlet states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO(nspin)

! Local variables

  double precision,parameter   :: th = 0.1d0

  double precision,allocatable :: HOMO(:), LUMO(:)
  double precision,allocatable :: HOMOa(:),LUMOa(:)
  double precision,allocatable :: HOMOb(:),LUMOb(:)

! Output variables

  double precision,intent(inout):: c(nBas,nBas,nspin)

! Memory allocation

  allocate(HOMO(nBas), LUMO(nBas),  &
           HOMOa(nBas),LUMOa(nBas), &
           HOMOb(nBas),LUMOb(nBas))

! Perform HOMO and LUMO rotation for guess mixing

  HOMO(:) = c(:,nO(1),1)
  LUMO(:) = c(:,nO(1)+1,1)

  HOMOa(:) =  cos(th)*HOMO(:) + sin(th)*LUMO(:)
  HOMOb(:) =  cos(th)*HOMO(:) - sin(th)*LUMO(:)

  LUMOa(:) = -sin(th)*HOMO(:) + cos(th)*LUMO(:)
  LUMOb(:) =  sin(th)*HOMO(:) + cos(th)*LUMO(:)

  c(:,nO(1),1)   = HOMOa(:)
  c(:,nO(1)+1,1) = LUMOa(:)

  c(:,nO(2),2)   = HOMOb(:)
  c(:,nO(2)+1,2) = LUMOb(:)

end subroutine 
