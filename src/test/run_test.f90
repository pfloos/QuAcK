subroutine run_test(doRtest,doUtest,doGtest)

  implicit none

! Input variables

  logical,intent(in)                 :: doRtest
  logical,intent(in)                 :: doUtest
  logical,intent(in)                 :: doGtest
 
! Local variables

  double precision              :: start_test     ,end_test       ,t_test

! Output variables

  if(doRtest) then 
  
    write(*,*) '****************************************'
    write(*,*) '* Testing Restricted Branch of QuAcK...*'
    write(*,*) '****************************************'
    write(*,*) 

    call wall_time(start_test)
    call check_test_value('R')
    call wall_time(end_test)

    t_test = end_test - start_test
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for test of restricted branch = ',t_test,' seconds'

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

    call wall_time(start_test)
    call check_test_value('U')
    call wall_time(end_test)

    t_test = end_test - start_test
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for test of unrestricted branch = ',t_test,' seconds'

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

    call wall_time(start_test)
    call check_test_value('G')
    call wall_time(end_test)

    t_test = end_test - start_test
    write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for test of generalized branch = ',t_test,' seconds'

    write(*,*) 
    write(*,*) '***************************'
    write(*,*) '* End of Generalized Test *'
    write(*,*) '***************************'
    write(*,*) 

  end if

end subroutine
