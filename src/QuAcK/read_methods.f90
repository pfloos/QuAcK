subroutine read_methods(doRHF,doUHF,doKS,doMOM,           & 
                        doMP2,doMP3,doMP2F12,             & 
                        doCCD,doDCD,doCCSD,doCCSDT,       & 
                        do_drCCD,do_rCCD,do_lCCD,do_pCCD, &
                        doCIS,doCIS_D,doCID,doCISD,doFCI, & 
                        doRPA,doRPAx,doppRPA,             & 
                        doG0F2,doevGF2,doqsGF2,           &
                        doG0F3,doevGF3,                   & 
                        doG0W0,doevGW,doqsGW,             & 
                        doG0T0,doevGT,doqsGT,             & 
                        doMCMP2)

! Read desired methods 

  implicit none

! Input variables

  logical,intent(out)           :: doRHF,doUHF,doKS,doMOM
  logical,intent(out)           :: doMP2,doMP3,doMP2F12
  logical,intent(out)           :: doCCD,doDCD,doCCSD,doCCSDT
  logical,intent(out)           :: do_drCCD,do_rCCD,do_lCCD,do_pCCD
  logical,intent(out)           :: doCIS,doCIS_D,doCID,doCISD,doFCI
  logical,intent(out)           :: doRPA,doRPAx,doppRPA
  logical,intent(out)           :: doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3  
  logical,intent(out)           :: doG0W0,doevGW,doqsGW
  logical,intent(out)           :: doG0T0,doevGT,doqsGT
  logical,intent(out)           :: doMCMP2

! Local variables

  character(len=1)              :: answer1,answer2,answer3,answer4,answer5,answer6,answer7

! Open file with method specification

  open(unit=1,file='input/methods')

! Set all the booleans to false

  doRHF = .false.
  doUHF = .false.
  doKS  = .false.
  doMOM = .false.

  doMP2    = .false.
  doMP3    = .false.
  doMP2F12 = .false.
 
  doCCD   = .false.
  doDCD   = .false.
  doCCSD  = .false.
  doCCSDT = .false.

  do_drCCD = .false.
  do_rCCD  = .false.
  do_lCCD  = .false.
  do_pCCD  = .false.

  doCIS   = .false.
  doCIS_D = .false.
  doCID   = .false.
  doCISD  = .false.
  doFCI   = .false.

  doRPA   = .false.
  doRPAx  = .false.
  doppRPA = .false.

  doG0F2  = .false.
  doevGF2 = .false.
  doqsGF2 = .false.
  doG0F3  = .false.
  doevGF3 = .false.

  doG0W0 = .false.
  doevGT = .false.
  doqsGT = .false.

  doG0T0 = .false.
  doevGW = .false.
  doqsGW = .false.

  doMCMP2 = .false.

! Read mean-field methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4
  if(answer1 == 'T') doRHF = .true.
  if(answer2 == 'T') doUHF = .true.
  if(answer3 == 'T') doKS  = .true.
  if(answer4 == 'T') doMOM = .true.

! Read MPn methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3
  if(answer1 == 'T') doMP2    = .true.
  if(answer2 == 'T') doMP3    = .true.
  if(answer3 == 'T') doMP2F12 = .true.

! Read CC methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4
  if(answer1 == 'T') doCCD    = .true.
  if(answer2 == 'T') doDCD    = .true.
  if(answer3 == 'T') doCCSD   = .true.
  if(answer4 == 'T') doCCSDT  = .true.

! Read weird CC methods
  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4
  if(answer1 == 'T') do_drCCD = .true.
  if(answer2 == 'T') do_rCCD  = .true.
  if(answer3 == 'T') do_lCCD  = .true.
  if(answer4 == 'T') do_pCCD  = .true.

! Read excited state methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4,answer5
  if(answer1 == 'T') doCIS   = .true.
  if(answer2 == 'T') doCIS_D = .true.
  if(answer3 == 'T') doCID   = .true.
  if(answer4 == 'T') doCISD  = .true.
  if(answer5 == 'T') doFCI   = .true.
  if(doCIS_D)        doCIS   = .true.

  read(1,*) 
  read(1,*) answer1,answer2,answer3
  if(answer1 == 'T') doRPA   = .true.
  if(answer2 == 'T') doRPAx  = .true.
  if(answer3 == 'T') doppRPA = .true.

! Read Green function methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4,answer5
  if(answer1 == 'T') doG0F2  = .true.
  if(answer2 == 'T') doevGF2 = .true.
  if(answer3 == 'T') doqsGF2 = .true.
  if(answer4 == 'T') doG0F3  = .true.
  if(answer5 == 'T') doevGF3 = .true.

! Read GW methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3
  if(answer1 == 'T') doG0W0 = .true.
  if(answer2 == 'T') doevGW = .true.
  if(answer3 == 'T') doqsGW = .true.

! Read GT methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3
  if(answer1 == 'T') doG0T0 = .true.
  if(answer2 == 'T') doevGT = .true.
  if(answer3 == 'T') doqsGT = .true.

! Read stochastic methods

  read(1,*) 
  read(1,*) answer1
  if(answer1 == 'T') doMCMP2 = .true.

! Close file with geometry specification

  close(unit=1)

end subroutine read_methods
