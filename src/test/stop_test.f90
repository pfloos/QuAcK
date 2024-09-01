subroutine stop_test(doRtest,doUtest,doGtest)

  implicit none

! Input variables

  logical,intent(in)                 :: doRtest
  logical,intent(in)                 :: doUtest
  logical,intent(in)                 :: doGtest
 
! Local variables

! Output variables

  if(doRtest) close(unit=1231597)

  if(doUtest) close(unit=1231597)

  if(doGtest) close(unit=1234181)

end subroutine
