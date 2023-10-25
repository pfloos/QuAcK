subroutine read_methods(doRHF,doUHF,doGHF,doROHF,          &
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

  logical,intent(out)           :: doRHF,doUHF,doGHF,doROHF
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

  character(len=1)              :: ans1,ans2,ans3,ans4,ans5,ans6

! Open file with method specification

  open(unit=1,file='input/methods')

! Read mean-field methods

  doRHF  = .false.
  doUHF  = .false.
  doGHF  = .false.
  doROHF = .false.

  read(1,*) 
  read(1,*) ans1,ans2,ans3,ans4
  if(ans1 == 'T') doRHF  = .true.
  if(ans2 == 'T') doUHF  = .true.
  if(ans3 == 'T') doGHF  = .true.
  if(ans4 == 'T') doROHF = .true.

! Read MPn methods

  doMP2    = .false.
  doMP3    = .false.

  read(1,*) 
  read(1,*) ans1,ans2
  if(ans1 == 'T') doMP2    = .true.
  if(ans2 == 'T') doMP3    = .true.

! Read CC methods

  doCCD   = .false.
  dopCCD  = .false.
  doDCD   = .false.
  doCCSD  = .false.
  doCCSDT = .false.

  read(1,*) 
  read(1,*) ans1,ans2,ans3,ans4,ans5
  if(ans1 == 'T') doCCD   = .true.
  if(ans2 == 'T') dopCCD  = .true.
  if(ans3 == 'T') doDCD   = .true.
  if(ans4 == 'T') doCCSD  = .true.
  if(ans5 == 'T') doCCSDT = .true.

! Read weird CC methods

  do_drCCD = .false.
  do_rCCD  = .false.
  do_crCCD = .false.
  do_lCCD  = .false.

  read(1,*) 
  read(1,*) ans1,ans2,ans3,ans4
  if(ans1 == 'T') do_drCCD = .true.
  if(ans2 == 'T') do_rCCD  = .true.
  if(ans3 == 'T') do_crCCD = .true.
  if(ans4 == 'T') do_lCCD  = .true.

! Read CI methods

  doCIS   = .false.
  doCIS_D = .false.
  doCID   = .false.
  doCISD  = .false.
  doFCI   = .false.

  read(1,*) 
  read(1,*) ans1,ans2,ans3,ans4,ans5
  if(ans1 == 'T') doCIS   = .true.
  if(ans2 == 'T') doCIS_D = .true.
  if(ans3 == 'T') doCID   = .true.
  if(ans4 == 'T') doCISD  = .true.
  if(ans5 == 'T') doFCI   = .true.
  if(doCIS_D)        doCIS   = .true.

! Read RPA methods

  dophRPA  = .false.
  dophRPAx = .false.
  docrRPA  = .false.
  doppRPA  = .false.

  read(1,*) 
  read(1,*) ans1,ans2,ans3,ans4
  if(ans1 == 'T') dophRPA  = .true.
  if(ans2 == 'T') dophRPAx = .true.
  if(ans3 == 'T') docrRPA  = .true.
  if(ans4 == 'T') doppRPA  = .true.

! Read Green's function methods

  doG0F2  = .false.
  doevGF2 = .false.
  doqsGF2 = .false.
  doG0F3  = .false.
  doevGF3 = .false.

  read(1,*) 
  read(1,*) ans1,ans2,ans3,ans4,ans5
  if(ans1 == 'T') doG0F2  = .true.
  if(ans2 == 'T') doevGF2 = .true.
  if(ans3 == 'T') doqsGF2 = .true.
  if(ans4 == 'T') doG0F3  = .true.
  if(ans5 == 'T') doevGF3 = .true.

! Read GW methods

  doG0W0    = .false.
  doevGW    = .false.
  doqsGW    = .false.
  doSRGqsGW = .false.
  doufG0W0  = .false.
  doufGW    = .false.

  read(1,*) 
  read(1,*) ans1,ans2,ans3,ans4,ans5,ans6
  if(ans1 == 'T') doG0W0    = .true.
  if(ans2 == 'T') doevGW    = .true.
  if(ans3 == 'T') doqsGW    = .true.
  if(ans4 == 'T') doSRGqsGW = .true.
  if(ans5 == 'T') doufG0W0  = .true.
  if(ans6 == 'T') doufGW    = .true.

! Read GT methods

  doG0T0pp = .false.
  doevGTpp = .false.
  doqsGTpp = .false.
  doG0T0eh = .false.
  doevGTeh = .false.
  doqsGTeh = .false.

  read(1,*) 
  read(1,*) ans1,ans2,ans3,ans4,ans5,ans6
  if(ans1 == 'T') doG0T0pp   = .true.
  if(ans2 == 'T') doevGTpp   = .true.
  if(ans3 == 'T') doqsGTpp   = .true.
  if(ans4 == 'T') doG0T0eh   = .true.
  if(ans5 == 'T') doevGTeh   = .true.
  if(ans6 == 'T') doqsGTeh   = .true.

! Close file with geometry specification

  close(unit=1)

end subroutine 
