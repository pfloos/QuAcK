subroutine read_MOs(nBas,C,e,EJ)

! Read normalization factor and MOs (coefficients and eigenvalues)

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas

! Local variables

  integer                       :: i,j

! Output variables

  double precision,intent(out)  :: EJ
  double precision,intent(out)  :: C(nBas,nBas),e(nBas)

!------------------------------------------------------------------------
! Primary basis set information
!------------------------------------------------------------------------

! Open file with basis set specification

  open(unit=3,file='input/MOs')

! Read MO information

  read(3,*) EJ

  do i=1,nBas
    read(3,*) (C(i,j),j=1,nBas)
  enddo

  do i=1,nBas
    read(3,*) e(i)
  enddo

! Print results

  write(*,'(A28)') '----------------------'
  write(*,'(A28)') 'MO coefficients'
  write(*,'(A28)') '----------------------'
  call matout(nBas,nBas,C)
  write(*,*)
  write(*,'(A28)') '----------------------'
  write(*,'(A28)') 'MO energies'
  write(*,'(A28)') '----------------------'
  call matout(nBas,1,e)
  write(*,*)

! Close file

  close(unit=3)

end subroutine read_MOs
