subroutine init_test(doRtest,doUtest,doGtest)

  implicit none

! Input variables

  logical,intent(in)                 :: doRtest
  logical,intent(in)                 :: doUtest
  logical,intent(in)                 :: doGtest
 
! Local variables

! Output variables

  if(doRtest) open(unit=11,file='test/Rtest.dat')

  if(doUtest) open(unit=12,file='test/Utest.dat')

  if(doGtest) open(unit=13,file='test/Gtest.dat')

end subroutine
