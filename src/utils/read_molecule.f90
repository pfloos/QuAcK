subroutine read_molecule(working_dir,nNuc,nO,nC,nR)

! Read number of atoms and number of electrons 

  implicit none

  include 'parameters.h'

! Local variables

  integer                       :: nCore
  integer                       :: nRyd

! Input variables

  character(len=256),intent(in) :: working_dir

! Output variables

  integer,intent(out)           :: nNuc
  integer,intent(out)           :: nO(nspin)
  integer,intent(out)           :: nC(nspin)
  integer,intent(out)           :: nR(nspin)

  integer                       :: status
  character(len=256)            :: file_path

! Open file with geometry specification

  file_path = trim(working_dir) // '/input/molecule'
  open(unit=1, file=file_path, status='old', action='read', iostat=status)

    if(status /= 0) then

      print *, "Error opening file: ", file_path
      stop

    else

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

    endif

  ! Close file with geometry specification
  close(unit=1)

end subroutine 
