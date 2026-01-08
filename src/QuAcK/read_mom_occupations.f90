subroutine read_mom_occupations(working_dir,eta_cap)

! Read mom occupations for guess

  implicit none

! Input variables

  character(len=256),intent(in) :: working_dir

! Output variables
  
  double precision,intent(out)  :: eta_cap  

! Local variables

  integer                       :: status
  character(len=256)            :: file_path


  file_path = trim(working_dir) // '/input/eta_opt.dat'
  open(unit=1, file=file_path, status='old', action='read', iostat=status)

    if(status /= 0) then

      print *, "Error opening file: ", file_path
      stop

    else

      read(1,*) eta_cap
    endif

  ! Close file
  close(unit=1)

end subroutine 
