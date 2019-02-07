subroutine read_molecule(nAt,nEl,nC,nO,nR)

! Read number of atoms nAt, 
!      number of electrons nEl, 
!      number of core electrons nC,
!      number of Rydberg orbitals nR

  implicit none

! Input variables
  integer,intent(out)           :: nAt,nEl,nC,nO,nR

! Open file with geometry specification

  open(unit=1,file='input/molecule')

! Read number of atoms and number of electrons

  read(1,*) 
  read(1,*) nAt,nEl,nC,nR

! Number of occupied orbitals

  if(mod(nEl,2) /= 0) then
    write(*,*) 'closed-shell system required!'
!   stop
  endif
  nO = nEl/2

! Number of core orbitals

  if(mod(nC,2) /= 0) then
    write(*,*) 'Number of core electrons not even!'
    stop
  endif
  nC = nC/2

  if(nC > nO) then
    write(*,*) 'Number of core electrons greater than number of electrons!'
    stop
  endif

! Print results

  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of atoms',nAt
  write(*,'(A28)') '----------------------'
  write(*,*)
  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of electrons',nEl
  write(*,'(A28)') '----------------------'
  write(*,*)
  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of core electrons',2*nC
  write(*,'(A28)') '----------------------'
  write(*,*)
  write(*,*)
  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of Rydberg orbitals',nR
  write(*,'(A28)') '----------------------'
  write(*,*)

! Close file with geometry specification

  close(unit=1)

end subroutine read_molecule
