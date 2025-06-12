subroutine write_matout(m,n,A,fpath)

! Print the MxN array A

  implicit none

  integer,intent(in)            :: m,n
  double precision,intent(in)   :: A(m,n)
  character(len=256),intent(in) :: fpath
  integer                       :: i,j
  integer                       :: unit

  unit = 69

  ! Open the file for writing
  open(unit=unit, file=fpath, status='unknown', action='write', form='formatted')

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
