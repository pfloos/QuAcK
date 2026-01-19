subroutine read_mom_occupations(working_dir,nO,doMOM,occupations)

! Reads the occupations of spin up and down for the MOM guess from the input dir, if doMOM is true and otherwise initializes
! the occupation following the aufbau-principle

  implicit none
  include 'parameters.h'

  character(len=256),intent(in) :: working_dir
  integer,intent(in)            :: nO(nspin)
  logical,intent(in)            :: doMOM
  integer,intent(out)           :: occupations(maxval(nO),nspin)

  integer                       :: status
  character(len=256)            :: file_path
  character(len=1024)           :: line
  integer                       :: ispin,i,nread
  integer                       :: tmp(maxval(nO))

  
  occupations(:,:) = 0

  if(doMOM) then
    file_path = trim(working_dir) // '/input/mom_occupations'

    open(unit=1, file=file_path, status='old', action='read', iostat=status)
    if (status /= 0) then
      write(*,*) "Error in opening mom_occupations file !!!"
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

    end do
  
    if(occupations(nO(1),1)==0) then
      print *, "Not enough alpha occupations provided !"
      print *, "Keep in mind that number of alpha electrons has to be >= number of beta electrons"
      print *, "Alpha occupation:"
      print *, occupations(:,1)
      error stop
    end if

    if(occupations(nO(2),2)==0) then
      print *, "Not enough beta occupations provided !"
      print *, "Beta occupation:"
      print *, occupations(:,2)
      error stop
    end if
 
    close(1)
  else

    do ispin=1,nspin
      do i=1,nO(ispin)
        occupations(i,ispin) = i
      end do
    end do

  end if
end subroutine
