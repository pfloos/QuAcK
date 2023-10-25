subroutine read_molecule(nNuc,nO,nC,nR)

! Read number of atoms and number of electrons 

  implicit none

  include 'parameters.h'

! Local variables

  integer                       :: nCore
  integer                       :: nRyd

! Output variables

  integer,intent(out)           :: nNuc
  integer,intent(out)           :: nO(nspin)
  integer,intent(out)           :: nC(nspin)
  integer,intent(out)           :: nR(nspin)

! Open file with geometry specification

  open(unit=1,file='input/molecule')

! Read number of atoms and number of electrons

  read(1,*) 
  read(1,*) nNuc,nO(1),nO(2),nCore,nRyd

  if(mod(nCore,2) /= 0 .or. mod(nRyd,2) /= 0) then

    print*, 'The number of core and Rydberg electrons must be even!'
    stop

  end if

  nC(:) = nCore/2
  nR(:) = nRyd/2

! Print results

  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of atoms',nNuc
  write(*,'(A28)') '----------------------'
  write(*,*)
  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of spin-up   electrons',nO(1)
  write(*,'(A28,1X,I16)') 'Number of spin-down electrons',nO(2)
  write(*,'(A28,1X,I16)') '    Total number of electrons',sum(nO(:))
  write(*,'(A28)') '----------------------'
  write(*,*)
  write(*,'(A28)') '----------------------'
  write(*,'(A28,1X,I16)') 'Number of core      electrons',sum(nC(:))
  write(*,'(A28,1X,I16)') 'Number of Rydberg   electrons',sum(nR(:))
  write(*,'(A28)') '----------------------'
  write(*,*)

! Close file with geometry specification

  close(unit=1)

end subroutine read_molecule
