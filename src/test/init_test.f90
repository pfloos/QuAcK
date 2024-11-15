subroutine init_test(working_dir,doRtest,doUtest,doGtest)

  implicit none

! Input variables

  character(len=256),intent(in) :: working_dir
  logical,intent(in)            :: doRtest
  logical,intent(in)            :: doUtest
  logical,intent(in)            :: doGtest
 
! Local variables

  integer                       :: status
  character(len=256)            :: file_path

! Output variables

  if(doRtest) then
    file_path = trim(working_dir) // '/test/Rtest.dat'
    open(unit=1231597, file=file_path, iostat=status)
    if(status /= 0) then
      print *, "Error opening file: ", file_path
      stop
    else
      print *, "opening file successed"
    endif
  endif

  if(doUtest) then
    file_path = trim(working_dir) // '/test/Utest.dat'
    open(unit=1232584, file=file_path)
  endif

  if(doGtest) then
    file_path = trim(working_dir) // '/test/Gtest.dat'
    open(unit=1234181, file=file_path)
  endif

end subroutine

