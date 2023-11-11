subroutine dump_test_value(branch,description,value)

  implicit none

! Input variables

  character(len=1),intent(in)        :: branch
  character(len=*),intent(in)        :: description

  double precision,intent(in)        :: value
 
! Local variables

! Output variables

  if(branch == 'R') then 
    open(unit=11,file='test/Rtest.dat')
  elseif(branch == 'U') then 
    open(unit=11,file='test/Utest.dat')
  elseif(branch == 'G') then
    open(unit=11,file='test/Gtest.dat')
  else
    write(*,*) 'Wrong branch name in dump_test_value'
  end if

  write(11,*) '# ',trim(description)
  write(11,'(F20.15)') value

  close(unit=11)

end subroutine
