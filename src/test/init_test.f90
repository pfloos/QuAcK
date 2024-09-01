subroutine init_test(doRtest,doUtest,doGtest)

  implicit none

! Input variables

  logical,intent(in)                 :: doRtest
  logical,intent(in)                 :: doUtest
  logical,intent(in)                 :: doGtest
 
! Local variables

! Output variables

  if(doRtest) open(unit=1231597, file='test/Rtest.dat')

  if(doUtest) open(unit=1232584, file='test/Utest.dat')

  if(doGtest) open(unit=1234181, file='test/Gtest.dat')

end subroutine
