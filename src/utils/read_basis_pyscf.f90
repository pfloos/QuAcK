subroutine read_basis_pyscf(working_dir,nBas,nO,nV)

! Read basis set information from PySCF

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nO(nspin)
  character(len=256),intent(in) :: working_dir

! Output variables

  integer,intent(out)            :: nV(nspin)
  integer,intent(out)            :: nBas

! Local variables

  integer                       :: status
  character(len=256)            :: file_path

!------------------------------------------------------------------------
! Primary basis set information
!------------------------------------------------------------------------

  file_path = trim(working_dir) // '/int/nBas.dat'
  open(unit=3, file=file_path, status='old', action='read', iostat=status)

    if(status /= 0) then

      print *, "Error opening file: ", file_path
      stop

    else

      read(3, *) nBas

    endif

  close(unit=3)

!  write(*,'(A38)') '--------------------------------------'
!  write(*,'(A38,1X,I16)') 'Number of basis functions (AOs)', nBas
!  write(*,'(A38)') '--------------------------------------'
!  write(*,*)

! Number of virtual orbitals

  nV(:) = nBas - nO(:)

end subroutine 
