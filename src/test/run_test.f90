subroutine run_test(doRtest,doUtest,doGtest)

  implicit none

! Input variables

  logical,intent(in)                 :: doRtest
  logical,intent(in)                 :: doUtest
  logical,intent(in)                 :: doGtest
 
! Local variables

! Output variables

  if(doRtest) then 
  
    write(*,*) '****************************************'
    write(*,*) '* Testing Restricted Branch of QuAcK...*'
    write(*,*) '****************************************'
    write(*,*) 

    call check_test_value('R')

    write(*,*) 
    write(*,*) '**************************'
    write(*,*) '* End of Restricted Test *'
    write(*,*) '**************************'
    write(*,*) 

  end if

  if(doUtest) then

    write(*,*) '******************************************'
    write(*,*) '* Testing Unrestricted Branch of QuAcK...*'
    write(*,*) '******************************************'
    write(*,*) 

    call check_test_value('U')

    write(*,*) 
    write(*,*) '****************************'
    write(*,*) '* End of Unrestricted Test *'
    write(*,*) '****************************'
    write(*,*) 

  end if

  if(doGtest) then 

    write(*,*) '*****************************************'
    write(*,*) '* Testing Generalized Branch of QuAcK...*'
    write(*,*) '*****************************************'
    write(*,*) 

    call check_test_value('G')

    write(*,*) 
    write(*,*) '***************************'
    write(*,*) '* End of Generalized Test *'
    write(*,*) '***************************'
    write(*,*) 

  end if

end subroutine
