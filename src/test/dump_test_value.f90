subroutine dump_test_value(branch,description,value)

  implicit none

! Input variables

  character(len=1),intent(in)        :: branch
  character(len=*),intent(in)        :: description

  double precision,intent(in)        :: value
 
! Local variables

! Output variables

  if(branch == 'R') then 

    write(11,*) trim(description)
    write(11,'(F20.15)') value

  elseif(branch == 'U') then 

    write(12,*) trim(description)
    write(12,'(F20.15)') value

  elseif(branch == 'G') then

    write(13,*) trim(description)
    write(13,'(F20.15)') value

  else

    write(*,*) 'Wrong branch name in dump_test_value'

  end if

end subroutine
