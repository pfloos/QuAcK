subroutine read_options(rung,SGn)

! Read DFT options

  implicit none

! Input variables

  integer,intent(out)           :: rung
  integer,intent(out)           :: SGn

! Open file with method specification

  open(unit=1,file='input/options')

! Default values

  rung = 1
  SGn  = 0

! Read rung of Jacob's ladder

  read(1,*)
  read(1,*) rung

! Read SG-n grid

  read(1,*)
  read(1,*) SGn

end subroutine read_options
