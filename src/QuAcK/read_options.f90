subroutine read_options(maxSCF_HF,thresh_HF,DIIS_HF,max_diis_HF,guess_type,ortho_type,mix,level_shift,dostab, &
                        reg_MP,                                                                               &
                        maxSCF_CC,thresh_CC,DIIS_CC,max_diis_CC,                                              &
                        TDA,singlet,triplet,spin_conserved,spin_flip,                                         &
                        maxSCF_GF,thresh_GF,DIIS_GF,max_diis_GF,lin_GF,eta_GF,renorm_GF,reg_GF,               &
                        maxSCF_GW,thresh_GW,DIIS_GW,max_diis_GW,lin_GW,eta_GW,reg_GW,TDA_W,                   &
                        maxSCF_GT,thresh_GT,DIIS_GT,max_diis_GT,lin_GT,eta_GT,reg_GT,TDA_T,                   &
                        doACFDT,exchange_kernel,doXBS,                                                        &
                        dophBSE,dophBSE2,doppBSE,dBSE,dTDA)

! Read desired methods 

  implicit none

! Input variables

  integer,intent(out)           :: maxSCF_HF
  double precision,intent(out)  :: thresh_HF
  logical,intent(out)           :: DIIS_HF
  integer,intent(out)           :: max_diis_HF
  integer,intent(out)           :: guess_type
  integer,intent(out)           :: ortho_type
  logical,intent(out)           :: mix
  double precision,intent(out)  :: level_shift
  logical,intent(out)           :: dostab

  logical,intent(out)           :: reg_MP

  integer,intent(out)           :: maxSCF_CC
  double precision,intent(out)  :: thresh_CC
  logical,intent(out)           :: DIIS_CC
  integer,intent(out)           :: max_diis_CC

  logical,intent(out)           :: TDA
  logical,intent(out)           :: singlet
  logical,intent(out)           :: triplet
  logical,intent(out)           :: spin_conserved
  logical,intent(out)           :: spin_flip

  integer,intent(out)           :: maxSCF_GF
  double precision,intent(out)  :: thresh_GF
  logical,intent(out)           :: DIIS_GF
  integer,intent(out)           :: max_diis_GF
  logical,intent(out)           :: lin_GF
  integer,intent(out)           :: renorm_GF
  double precision,intent(out)  :: eta_GF
  logical,intent(out)           :: reg_GF

  integer,intent(out)           :: maxSCF_GW
  double precision,intent(out)  :: thresh_GW
  logical,intent(out)           :: DIIS_GW
  integer,intent(out)           :: max_diis_GW
  logical,intent(out)           :: TDA_W
  logical,intent(out)           :: lin_GW
  double precision,intent(out)  :: eta_GW
  logical,intent(out)           :: reg_GW

  integer,intent(out)           :: maxSCF_GT
  double precision,intent(out)  :: thresh_GT
  logical,intent(out)           :: DIIS_GT
  integer,intent(out)           :: max_diis_GT
  logical,intent(out)           :: TDA_T
  logical,intent(out)           :: lin_GT
  double precision,intent(out)  :: eta_GT
  logical,intent(out)           :: reg_GT

  logical,intent(out)           :: doACFDT
  logical,intent(out)           :: exchange_kernel
  logical,intent(out)           :: doXBS

  logical,intent(out)           :: dophBSE
  logical,intent(out)           :: dophBSE2
  logical,intent(out)           :: doppBSE
  logical,intent(out)           :: dBSE
  logical,intent(out)           :: dTDA

! Local variables

  character(len=1)              :: ans1,ans2,ans3,ans4,ans5

! Open file with method specification

  open(unit=1,file='input/options')

! Read HF options

  maxSCF_HF    = 64
  thresh_HF    = 1d-6
  DIIS_HF      = .false.
  max_diis_HF    = 5
  guess_type   = 1
  ortho_type   = 1
  mix          = .false.
  level_shift  = 0d0
  dostab       = .false.

  read(1,*) 
  read(1,*) maxSCF_HF,thresh_HF,ans1,max_diis_HF,guess_type,ortho_type,ans2,level_shift,ans3

  if(ans1 == 'T') DIIS_HF     = .true.
  if(ans2 == 'T') mix         = .true.
  if(ans3 == 'T') dostab      = .true.

  if(.not.DIIS_HF) max_diis_HF = 1

! Read MPn options

  reg_MP = .false.
  read(1,*)
  read(1,*) ans1

  if(ans1 == 'T') reg_MP = .true.

! Read CC options

  maxSCF_CC    = 64
  thresh_CC    = 1d-5
  DIIS_CC      = .false.
  max_diis_CC    = 5

  read(1,*)
  read(1,*) maxSCF_CC,thresh_CC,ans1,max_diis_CC

  if(ans1 == 'T') DIIS_CC    = .true.

  if(.not.DIIS_CC) max_diis_CC = 1

! Read excited state options

  TDA            = .false.
  singlet        = .false.
  triplet        = .false.
  spin_conserved = .false.
  spin_flip      = .false.

  read(1,*) 
  read(1,*) ans1,ans2,ans3,ans4,ans5

  if(ans1 == 'T') TDA            = .true.
  if(ans2 == 'T') singlet        = .true.
  if(ans3 == 'T') triplet        = .true.
  if(ans4 == 'T') spin_conserved = .true.
  if(ans5 == 'T') spin_flip      = .true.

! Read GF options

  maxSCF_GF   = 64
  thresh_GF   = 1d-5
  DIIS_GF     = .false.
  max_diis_GF = 5
  lin_GF      = .false.
  eta_GF      = 0d0
  renorm_GF   = 0
  reg_GF      = .false.

  read(1,*) 
  read(1,*) maxSCF_GF,thresh_GF,ans1,max_diis_GF,ans2,eta_GF,renorm_GF,ans3

  if(ans1 == 'T') DIIS_GF = .true.
  if(ans2 == 'T') lin_GF  = .true.
  if(ans3 == 'T') reg_GF  = .true.
  if(.not.DIIS_GF) max_diis_GF = 1

! Read GW options

  maxSCF_GW   = 64
  thresh_GW   = 1d-5
  DIIS_GW     = .false.
  max_diis_GW = 5
  lin_GW      = .false.
  eta_GW      = 0d0
  reg_GW      = .false.
  TDA_W       = .false.

  read(1,*) 
  read(1,*) maxSCF_GW,thresh_GW,ans1,max_diis_GW,ans2,eta_GW,ans3,ans4

  if(ans1 == 'T') DIIS_GW = .true.
  if(ans2 == 'T') lin_GW  = .true.
  if(ans3 == 'T') TDA_W   = .true.
  if(ans4 == 'T') reg_GW  = .true.
  if(.not.DIIS_GW) max_diis_GW = 1

! Read GT options

  maxSCF_GT   = 64
  thresh_GT   = 1d-5
  DIIS_GT     = .false.
  max_diis_GT = 5
  lin_GT      = .false.
  eta_GT      = 0d0
  reg_GT      = .false.
  TDA_T       = .false.

  read(1,*) 
  read(1,*) maxSCF_GT,thresh_GT,ans1,max_diis_GT,ans2,eta_GT,ans3,ans4

  if(ans1 == 'T') DIIS_GT = .true.
  if(ans2 == 'T') lin_GT  = .true.
  if(ans3 == 'T') TDA_T   = .true.
  if(ans4 == 'T') reg_GT  = .true.
  if(.not.DIIS_GT) max_diis_GT = 1

! Options for adiabatic connection

  doACFDT         = .false.
  exchange_kernel = .false.
  doXBS           = .false.

  read(1,*) 
  read(1,*) ans1,ans2,ans3

  if(ans1 == 'T') doACFDT = .true.
  if(ans2 == 'T') exchange_kernel = .true.
  if(ans3 == 'T') doXBS   = .true.

! Options for dynamical BSE calculations

  dophBSE  = .false.
  dophBSE2 = .false.
  doppBSE  = .false.
  dBSE     = .false.
  dTDA     = .true.

  read(1,*) 
  read(1,*) ans1,ans2,ans3,ans4,ans5

  if(ans1 == 'T') dophBSE  = .true.
  if(ans2 == 'T') dophBSE2 = .true.
  if(ans3 == 'T') doppBSE  = .true.
  if(ans4 == 'T') dBSE     = .true.
  if(ans5 == 'F') dTDA     = .false.

! Close file with options

  close(unit=1)

end subroutine 
