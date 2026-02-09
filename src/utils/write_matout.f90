subroutine write_matout(m,n,A,fpath)

! Print the MxN array A

  implicit none

  integer,intent(in)            :: m,n
  double precision,intent(in)   :: A(m,n)
  character(len=*),intent(in) :: fpath
  integer                       :: i,j
  integer                       :: unit,ios

  unit = 69

  ! Open the file for writing
  open(newunit=unit, file=fpath, status='replace', action='write', form='formatted',iostat=ios)
  if(ios/=0) then
    print *, "Error File could not be opened:", ios
    return
  endif
  do i = 1, m
    do j = 1, n
      write(unit, '(F20.12)', advance='no') A(i,j)
      if (j < n) write(unit, '(" ")', advance='no')
    end do
    write(unit, *)
  end do

  ! Close the file
  close(unit)

end subroutine
