subroutine check_test_value(branch)

  implicit none

! Input variables

  character(len=1),intent(in)        :: branch
 
! Local variables

  character(len=30)                  :: description
  double precision                   :: value
  double precision                   :: reference
  character(len=15)                  :: answer

  logical                            :: failed
  double precision,parameter         :: cutoff = 1d-10

! Output variables

  if(branch == 'R') then 

    open(unit=11,file='test/Rtest.dat')
    open(unit=12,file='test/Rtest_ref.dat')

  elseif(branch == 'U') then 

    open(unit=11,file='test/Utest.dat')
    open(unit=12,file='test/Utest_ref.dat')

  elseif(branch == 'G') then

    open(unit=11,file='test/Gtest.dat')
    open(unit=12,file='test/Gtest_ref.dat')

  else

    write(*,*) 'Wrong branch name in check_test_value'

  end if

  failed = .false.

  write(*,*) '----------------------------------------------------------------------------------------------------'
  do

    read(11,'(A30)',end=11) description
    read(11,'(F20.15)',end=11) value

    read(12,*,end=12) 
    read(12,'(F20.15)',end=12) reference

    if(abs(value-reference) < cutoff) then
      answer = '.......... :-)'
    else
      answer = '.......... :-( '
      failed = .true.
    end if
    write(*,'(1X,A1,1X,A30,1X,A1,1X,3F15.10,1X,A1,1X,A15,1X,A1)') & 
      '|',description,'|',value,reference,abs(value-reference),'|',answer,'|'

  end do

  11 close(unit=11)
  12 close(unit=12)

  write(*,*) '----------------------------------------------------------------------------------------------------'
  if(failed) then 
    write(*,'(3X,A30)') ' :-( Tests have failed    :-( '
  else 
    write(*,'(3X,A30)') ' :-) Tests have succeeded :-) '
  end if
  write(*,*) '----------------------------------------------------------------------------------------------------'

end subroutine
