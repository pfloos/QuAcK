subroutine read_molecule(nAt,nEl,nO)

! Read number of atoms nAt and number of electrons nEl

  implicit none

! Input variables
  integer,intent(out)           :: nAt,nEl,nO

! Open file with geometry specification

  open(unit=1,file='input/molecule')

! Read number of atoms and number of electrons

  read(1,*) 
  read(1,*) nAt,nEl

! Number of occupied orbitals

  if(mod(nEl,2) /= 0) then
    write(*,*) 'closed-shell system required!'
    stop
  endif
  nO = nEl/2

! Print results

  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of atoms',nAt
  write(*,'(A28)') '----------------------'
  write(*,*)
  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of electrons',nEl
  write(*,'(A28)') '----------------------'
  write(*,*)

! Close file with geometry specification

  close(unit=1)

end subroutine read_molecule
