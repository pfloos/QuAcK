subroutine RHF_mix_guess(nBas,nOrb,nO,mix,c)

!  Mix core orbitals in solution of RHF

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nOrb
  integer,intent(in)            :: nO(nspin)
  double precision,intent(in)   :: mix

! Local variables

  double precision,allocatable :: core1(:), core2(:)
  double precision,allocatable :: Core1rot(:),Core2rot(:)
  integer                      :: core1idx=1, core2idx=2

  double precision              :: th

! Output variables

  double precision,intent(inout):: c(nBas,nOrb)

! Memory allocation

  allocate(core1(nBas), core2(nBas),  &
           core1rot(nBas),core2rot(nBas))

! Perform HOMO and LUMO rotation for guess mixing

  th = 2d0*pi*mix

  core1(:) = c(:,core1idx)
  core2(:) = c(:,core2idx)
  
  core1rot(:) =  cos(th)*core1(:)  - sin(th)*core2(:)
  core2rot(:) =  cos(th)*core2(:)  + sin(th)*core1(:)

  c(:,core1idx) = core1rot(:)
  c(:,core2idx) = core2rot(:)

  deallocate(core1,core2,core1rot,core2rot)
end subroutine 
