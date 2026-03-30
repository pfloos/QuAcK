subroutine write_occupations(nOa,nOb,occupationsa,occupationsb,fpath)

! Print the MxN array A

  implicit none

  integer,intent(in)            :: nOa,nOb
  integer,intent(in)            :: occupationsa(nOa),occupationsb(nOb)
  character(len=*),intent(in)   :: fpath
  integer                       :: unit,ios,i

  unit = 69

  ! Open the file for writing
  open(newunit=unit, file=fpath, status='replace', action='write', form='formatted',iostat=ios)
  if(ios/=0) then
    print *, "Error File could not be opened:", ios
    return
  endif
  
  write(unit, *) "# Alpha occupations"
  write(unit, *) occupationsa(1:nOa)
  write(unit, *) "# Beta occupations"
  write(unit, *) occupationsb(1:nOb)

  ! Close the file
  close(unit)

end subroutine

