subroutine read_methods(doRHF,doUHF,doKS,doMOM,            & 
                        doMP2,doMP3,                       & 
                        doCCD,dopCCD,doDCD,doCCSD,doCCSDT, & 
                        do_drCCD,do_rCCD,do_crCCD,do_lCCD, &
                        doCIS,doCIS_D,doCID,doCISD,doFCI,  & 
                        dophRPA,dophRPAx,docrRPA,doppRPA,  & 
                        doG0F2,doevGF2,doqsGF2,            &
                        doG0F3,doevGF3,                    & 
                        doG0W0,doevGW,doqsGW,doSRGqsGW,    & 
                        doufG0W0,doufGW,                   & 
                        doG0T0pp,doevGTpp,doqsGTpp,        &
                        doG0T0eh,doevGTeh,doqsGTeh)

! Read desired methods 

  implicit none

! Input variables

  logical,intent(out)           :: doRHF,doUHF,doKS,doMOM
  logical,intent(out)           :: doMP2,doMP3
  logical,intent(out)           :: doCCD,dopCCD,doDCD,doCCSD,doCCSDT
  logical,intent(out)           :: do_drCCD,do_rCCD,do_crCCD,do_lCCD
  logical,intent(out)           :: doCIS,doCIS_D,doCID,doCISD,doFCI
  logical,intent(out)           :: dophRPA,dophRPAx,docrRPA,doppRPA
  logical,intent(out)           :: doG0F2,doevGF2,doqsGF2,doG0F3,doevGF3  
  logical,intent(out)           :: doG0W0,doevGW,doqsGW,doSRGqsGW,doufG0W0,doufGW
  logical,intent(out)           :: doG0T0pp,doevGTpp,doqsGTpp
  logical,intent(out)           :: doG0T0eh,doevGTeh,doqsGTeh

! Local variables

  character(len=1)              :: answer1,answer2,answer3,answer4,answer5,answer6

! Open file with method specification

  open(unit=1,file='input/methods')

! Set all the booleans to false

  doRHF = .false.
  doUHF = .false.
  doKS  = .false.
  doMOM = .false.

  doMP2    = .false.
  doMP3    = .false.
 
  doCCD   = .false.
  dopCCD  = .false.
  doDCD   = .false.
  doCCSD  = .false.
  doCCSDT = .false.

  do_drCCD = .false.
  do_rCCD  = .false.
  do_crCCD = .false.
  do_lCCD  = .false.

  doCIS   = .false.
  doCIS_D = .false.
  doCID   = .false.
  doCISD  = .false.
  doFCI   = .false.

  dophRPA  = .false.
  dophRPAx = .false.
  docrRPA  = .false.
  doppRPA  = .false.

  doG0F2  = .false.
  doevGF2 = .false.
  doqsGF2 = .false.
  doG0F3  = .false.
  doevGF3 = .false.

  doG0W0    = .false.
  doevGW    = .false.
  doqsGW    = .false.
  doSRGqsGW = .false.
  doufG0W0  = .false.
  doufGW    = .false.

  doG0T0pp = .false.
  doevGTpp = .false.
  doqsGTpp = .false.
  doG0T0eh = .false.
  doevGTeh = .false.
  doqsGTeh = .false.

! Read mean-field methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4
  if(answer1 == 'T') doRHF = .true.
  if(answer2 == 'T') doUHF = .true.
  if(answer3 == 'T') doKS  = .true.
  if(answer4 == 'T') doMOM = .true.

! Read MPn methods

  read(1,*) 
  read(1,*) answer1,answer2
  if(answer1 == 'T') doMP2    = .true.
  if(answer2 == 'T') doMP3    = .true.

! Read CC methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4,answer5
  if(answer1 == 'T') doCCD   = .true.
  if(answer2 == 'T') dopCCD  = .true.
  if(answer3 == 'T') doDCD   = .true.
  if(answer4 == 'T') doCCSD  = .true.
  if(answer5 == 'T') doCCSDT = .true.

! Read weird CC methods
  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4
  if(answer1 == 'T') do_drCCD = .true.
  if(answer2 == 'T') do_rCCD  = .true.
  if(answer3 == 'T') do_crCCD = .true.
  if(answer4 == 'T') do_lCCD  = .true.

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
  read(1,*) answer1,answer2,answer3,answer4
  if(answer1 == 'T') dophRPA  = .true.
  if(answer2 == 'T') dophRPAx = .true.
  if(answer3 == 'T') docrRPA  = .true.
  if(answer4 == 'T') doppRPA  = .true.

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
  read(1,*) answer1,answer2,answer3,answer4,answer5,answer6
  if(answer1 == 'T') doG0W0    = .true.
  if(answer2 == 'T') doevGW    = .true.
  if(answer3 == 'T') doqsGW    = .true.
  if(answer4 == 'T') doSRGqsGW = .true.
  if(answer5 == 'T') doufG0W0  = .true.
  if(answer6 == 'T') doufGW    = .true.

! Read GT methods

  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4,answer5,answer6
  if(answer1 == 'T') doG0T0pp   = .true.
  if(answer2 == 'T') doevGTpp   = .true.
  if(answer3 == 'T') doqsGTpp   = .true.
  if(answer4 == 'T') doG0T0eh   = .true.
  if(answer5 == 'T') doevGTeh   = .true.
  if(answer6 == 'T') doqsGTeh   = .true.

! Close file with geometry specification

  close(unit=1)

end subroutine read_methods
