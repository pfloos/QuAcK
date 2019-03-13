subroutine read_molecule(nAt,nEl,nO)

! Read number of atoms nAt and number of electrons nEl

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(out)           :: nAt,nEl(nspin),nO(nspin)

! Local variables

  integer                       :: n

! Open file with geometry specification

  open(unit=1,file='input/molecule')

! Read number of atoms and number of electrons

  read(1,*) 
  read(1,*) nAt,n,nEl(1),nEl(2)

! Check imput consistency

  if(n /= sum(nEl(:))) then
    write(*,*) 'number of electrons inconsistent'
    stop
  endif

  nO(:) = nEl(:)

! Print results

  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of atoms',nAt
  write(*,'(A28)') '----------------------'
  write(*,*)
  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of spin-up   electrons',nEl(1)
  write(*,'(A28,1X,I16)') 'Number of spin-down electrons',nEl(2)
  write(*,'(A28,1X,I16)') 'Total number of     electrons',sum(nEl(:))
  write(*,'(A28)') '----------------------'
  write(*,*)

! Close file with geometry specification

  close(unit=1)

end subroutine read_molecule
