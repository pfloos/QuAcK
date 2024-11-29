subroutine read_hardware(working_dir,use_gpu)

! Read desired methods 

  implicit none

! Input variables

  character(len=256),intent(in) :: working_dir

! Output variables

  logical,intent(out)           :: use_gpu

! Local variables

  character(len=1)              :: ans
  integer                       :: ios
  character(len=256)            :: file_path

! Open file with method specification

  file_path = trim(working_dir) // '/input/hardware'
  open(unit=1, file=file_path, status='old', action='read', iostat=ios)

    if(ios /= 0) then

      use_gpu = .False.

    else

      read(1,*) 
      read(1,*) ans
      if(ans == 'T') then
        use_gpu = .true.
      else
        use_gpu = .False.
      endif

    endif

  ! Close file with options
  close(unit=1)

end subroutine 
