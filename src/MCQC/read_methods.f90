subroutine read_methods(doHF,doMOM,           & 
                        doMP2,doMP3,doMP2F12, & 
                        doCC,                 & 
                        doCIS,doTDHF,doADC,   & 
                        doGF2,doGF3,          & 
                        doG0W0,doevGW,doqsGW, & 
                        doMCMP2)

! Read desired methods 

  implicit none

! Input variables

  logical,intent(out)           :: doHF,doMOM
  logical,intent(out)           :: doMP2,doMP3,doMP2F12
  logical,intent(out)           :: doCC
  logical,intent(out)           :: doCIS,doTDHF,doADC
  logical,intent(out)           :: doGF2,doGF3  
  logical,intent(out)           :: doG0W0,doevGW,doqsGW
  logical,intent(out)           :: doMCMP2

! Local variables

  character(len=1)              :: answer1,answer2,answer3

! Open file with method specification

  open(unit=1,file='input/methods')

! Set all the booleans to false

  doHF  = .false.
  doMOM = .false.

  doMP2    = .false.
  doMP3    = .false.
  doMP2F12 = .false.
 
  doCC = .false.

  doCIS  = .false.
  doTDHF = .false.
  doADC  = .false.

  doGF2 = .false.
  doGF3 = .false.

  doG0W0 = .false.
  doevGW = .false.
  doqsGW = .false.

  doMCMP2 = .false.

! Read mean-field methods

  read(1,*) 
  read(1,*) answer1,answer2
  if(answer1 == 'T') doHF  = .true.
  if(answer2 == 'T') doMOM = .true.

! Read MPn methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3
  if(answer1 == 'T') doMP2    = .true.
  if(answer2 == 'T') doMP3    = .true.
  if(answer3 == 'T') doMP2F12 = .true.

! Read CC methods

  read(1,*) 
  read(1,*) answer1
  if(answer1 == 'T') doCC    = .true.

! Read excited state methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3
  if(answer1 == 'T') doCIS  = .true.
  if(answer2 == 'T') doTDHF = .true.
  if(answer3 == 'T') doADC  = .true.

! Read Green function methods

  read(1,*) 
  read(1,*) answer1,answer2
  if(answer1 == 'T') doGF2 = .true.
  if(answer2 == 'T') doGF3 = .true.

! Read GW methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3
  if(answer1 == 'T') doG0W0 = .true.
  if(answer2 == 'T') doevGW = .true.
  if(answer3 == 'T') doqsGW = .true.

! Read stochastic methods

  read(1,*) 
  read(1,*) answer1
  if(answer1 == 'T') doMCMP2 = .true.

! Close file with geometry specification

  close(unit=1)

end subroutine read_methods
