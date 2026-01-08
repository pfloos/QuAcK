subroutine read_mom_occupations(working_dir,nO,occupations)

! Reads the occupations of spin up and down for the MOM guess from the input dir.

  implicit none
  include 'parameters.h'

  character(len=256),intent(in) :: working_dir
  integer,intent(in)            :: nO(nspin)
  integer,intent(out)           :: occupations(maxval(nO),nspin)

  integer                       :: status
  character(len=256)            :: file_path
  character(len=1024)           :: line
  integer                       :: ispin, nread
  integer                       :: tmp(maxval(nO))

  file_path = trim(working_dir) // '/input/mom_occupations'
  occupations(:,:) = 0

  open(unit=1, file=file_path, status='old', action='read', iostat=status)
  if (status /= 0) then
          write(*,*) "Error in opening mom_occupations file !!!"
          write(*,*) "If this file is not needed, please just provide an empty file with the name mom_occupations in the input dir."
          return
  end if

  do ispin = 1, nspin

    read(1,*,iostat=status)    ! comment line
    if (status /= 0) exit

    read(1,'(A)',iostat=status) line
    if (status /= 0) exit

    tmp = 0
    read(line,*,iostat=status) tmp

    occupations(:,ispin) = tmp(:)
    print *, occupations(:,ispin)

  end do
  close(1)

end subroutine
