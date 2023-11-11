subroutine stop_test(doRtest,doUtest,doGtest)

  implicit none

! Input variables

  logical,intent(in)                 :: doRtest
  logical,intent(in)                 :: doUtest
  logical,intent(in)                 :: doGtest
 
! Local variables

! Output variables

  if(doRtest) close(unit=11)

  if(doUtest) close(unit=12)

  if(doGtest) close(unit=13)

end subroutine
