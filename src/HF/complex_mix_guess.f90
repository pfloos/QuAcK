subroutine complex_mix_guess(nBas,nO,mix,c)

!  Guess mixing for UHF calculations on singlet states

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO(nspin)
  double precision,intent(in)   :: mix

! Local variables

  complex*16,allocatable        :: HOMO(:), LUMO(:)
  complex*16,allocatable        :: HOMOa(:),LUMOa(:)
  complex*16,allocatable        :: HOMOb(:),LUMOb(:)

  double precision              :: th

! Output variables

  complex*16,intent(inout)      :: c(nBas,nBas,nspin)

! Memory allocation

  allocate(HOMO(nBas), LUMO(nBas),  &
           HOMOa(nBas),LUMOa(nBas), &
           HOMOb(nBas),LUMOb(nBas))

! Perform HOMO and LUMO rotation for guess mixing

  th = 2d0*pi*mix

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
