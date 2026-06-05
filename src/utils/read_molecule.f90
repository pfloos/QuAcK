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

      write(*,*) '============================================================'
      write(*,*) '                     Molecular System                       '
      write(*,*) '============================================================'
      write(*,'(A35,I8)') ' Number of atoms                = ', nNuc
      write(*,*)
      write(*,'(A35,I8)') ' Spin-up electrons              = ', nO(1)
      write(*,'(A35,I8)') ' Spin-down electrons            = ', nO(2)
      write(*,'(A35,I8)') ' Total electrons                = ', sum(nO(:))
      write(*,*)
      write(*,'(A35,I8)') ' Core electrons                 = ', sum(nC(:))
      write(*,'(A35,I8)') ' Rydberg electrons              = ', sum(nR(:))
      write(*,*) '============================================================'
      write(*,*)
      
    endif
      
  ! Close file with geometry specification
  close(unit=1)

end subroutine 
