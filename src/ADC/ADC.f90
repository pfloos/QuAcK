subroutine ADC(singlet_manifold,triplet_manifold,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,e,ERI)

! ADC main routine

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: singlet_manifold,triplet_manifold
  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: e(nBas),ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin


! Hello world

  write(*,*)
  write(*,*)'**********************'
  write(*,*)'|    ADC(n) module   |'
  write(*,*)'**********************'
  write(*,*)

! ADC(2) calculation for singlet manifold

  if(singlet_manifold) then

    ispin = 1
    call ADC2(ispin,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,e,ERI)

  endif

! ADC(2) calculation for triplet manifold

  if(triplet_manifold) then

    ispin = 2
    call ADC2(ispin,maxSCF,thresh,max_diis,nBas,nC,nO,nV,nR,e,ERI)

  endif

end subroutine ADC
