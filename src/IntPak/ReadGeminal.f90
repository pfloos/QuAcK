subroutine ReadGeminal(ExpS)

! Read the geminal information

  implicit none

! Input variables
  double precision,intent(out)  :: ExpS

! Open file with geometry specification
  open(unit=4,file='input/geminal')

! Read exponent of Slater geminal
  read(4,*) ExpS


  write(*,'(A28)') '------------------'
  write(*,'(A28,1X,F16.10)') 'Slater geminal exponent',ExpS
  write(*,'(A28)') '------------------'
  write(*,*)

! Close file with geminal information
  close(unit=4)

end subroutine ReadGeminal
