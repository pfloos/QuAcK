subroutine read_options(debug,chemist_notation,ExpS,doOv,doKin,doNuc,doERI,doF12,doYuk,doErf,do3eInt,do4eInt)

! Read desired methods 

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(out)           :: debug
  logical,intent(out)           :: chemist_notation

  double precision,intent(out)  :: ExpS

  logical,intent(out)           :: doOv
  logical,intent(out)           :: doKin
  logical,intent(out)           :: doNuc

  logical,intent(out)           :: doERI
  logical,intent(out)           :: doF12
  logical,intent(out)           :: doYuk
  logical,intent(out)           :: doErf

  logical,intent(out)           :: do3eInt(n3eInt)
  logical,intent(out)           :: do4eInt(n4eInt)
  
! Local variables

  character(len=1)              :: answer1,answer2,answer3,answer4

! Open file with method specification

  open(unit=1,file='input/int')

! Read HF options

  debug = .false.
  chemist_notation = .false.
  ExpS = 1d0

  doOv  = .false.
  doKin = .false.
  doNuc = .false.

  doERI = .false.
  doF12 = .false.
  doYuk = .false.
  doErf = .false.

  do3eInt(1) = .false.
  do3eInt(2) = .false.
  do3eInt(3) = .false.

  do4eInt(1) = .false.
  do4eInt(2) = .false.
  do4eInt(3) = .false.
 
! Debugging mode 

  read(1,*) 
  read(1,*) answer1

  if(answer1 == 'T') debug = .true.

! Chem or Phys notations?

  read(1,*) 
  read(1,*) answer1

  if(answer1 == 'T') chemist_notation = .true.

! Read exponent of Slater geminal
  read(1,*) 
  read(1,*) ExpS

  write(*,'(A28)') '------------------'
  write(*,'(A28,1X,F16.10)') 'Slater geminal exponent',ExpS
  write(*,'(A28)') '------------------'
  write(*,*)

! One-electron integrals

  read(1,*) 
  read(1,*) answer1,answer2,answer3

  if(answer1 == 'T') doOv  = .true.
  if(answer2 == 'T') doKin = .true.
  if(answer3 == 'T') doNuc  = .true.

! Two-electron integrals

  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4

  if(answer1 == 'T') doERI = .true.
  if(answer2 == 'T') doF12 = .true.
  if(answer3 == 'T') doYuk  = .true.
  if(answer4 == 'T') doErf = .true.

! Three-electron integrals

  read(1,*) 
  read(1,*) answer1,answer2,answer3

  if(answer1 == 'T') do3eInt(1) = .true.
  if(answer2 == 'T') do3eInt(2) = .true.
  if(answer3 == 'T') do3eInt(3) = .true.

! Four-electron integrals

  read(1,*) 
  read(1,*) answer1,answer2,answer3

  if(answer1 == 'T') do4eInt(1) = .true.
  if(answer2 == 'T') do4eInt(2) = .true.
  if(answer3 == 'T') do4eInt(3) = .true.

! Close file with options

  close(unit=1)

end subroutine read_options
