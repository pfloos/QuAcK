subroutine dump_test_value(branch, description, val)

  implicit none

! Input variables

  character(len=1),intent(in)        :: branch
  character(len=*),intent(in)        :: description

  double precision,intent(in)        :: val
 
! Local variables

! Output variables

  if(branch == 'R') then 

    !write(1231597, '(A, ": ", F20.15)') '"' // trim(description) // '"', val
    write(1231597, *) trim(description)
    write(1231597, '(F20.15)') val

  elseif(branch == 'U') then 

    write(1232584,*) trim(description)
    write(1232584,'(F20.15)') val

  elseif(branch == 'G') then

    write(1234181,*) trim(description)
    write(1234181,'(F20.15)') val

  else

    write(*,*) 'Wrong branch name in dump_test_value'

  end if

end subroutine
