subroutine FCI(nBas,nC,nO,nV,nR,ERI,e)

! Perform a full configuration interaction calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

! Hello world

  write(*,*)
  write(*,*)'**********************************'
  write(*,*)'| Full Configuration Interaction |'
  write(*,*)'**********************************'
  write(*,*)

! Form FCI vector

! Form FCI matrix

! Diagonalize FCI matrix


end subroutine 
