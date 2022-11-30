subroutine read_options(maxSCF_HF,thresh_HF,DIIS_HF,n_diis_HF,guess_type,ortho_type,mix,level_shift,dostab, &
                        regMP,                                                                              &
                        maxSCF_CC,thresh_CC,DIIS_CC,n_diis_CC,                                              &
                        TDA,singlet,triplet,spin_conserved,spin_flip,                                       &
                        maxSCF_GF,thresh_GF,DIIS_GF,n_diis_GF,linGF,eta_GF,renormGF,regGF,                  &
                        maxSCF_GW,thresh_GW,DIIS_GW,n_diis_GW,linGW,eta_GW,regGW,                           &
                        COHSEX,SOSEX,TDA_W,                                                                 &
                        maxSCF_GT,thresh_GT,DIIS_GT,n_diis_GT,linGT,eta_GT,regGT,TDA_T,                     &
                        doACFDT,exchange_kernel,doXBS,                                                      &
                        BSE,dBSE,dTDA,evDyn,ppBSE,BSE2)

! Read desired methods 

  implicit none

! Input variables

  integer,intent(out)           :: maxSCF_HF
  double precision,intent(out)  :: thresh_HF
  logical,intent(out)           :: DIIS_HF
  integer,intent(out)           :: n_diis_HF
  integer,intent(out)           :: guess_type
  integer,intent(out)           :: ortho_type
  logical,intent(out)           :: mix
  double precision,intent(out)  :: level_shift
  logical,intent(out)           :: dostab

  logical,intent(out)           :: regMP

  integer,intent(out)           :: maxSCF_CC
  double precision,intent(out)  :: thresh_CC
  logical,intent(out)           :: DIIS_CC
  integer,intent(out)           :: n_diis_CC

  logical,intent(out)           :: TDA
  logical,intent(out)           :: singlet
  logical,intent(out)           :: triplet
  logical,intent(out)           :: spin_conserved
  logical,intent(out)           :: spin_flip

  integer,intent(out)           :: maxSCF_GF
  double precision,intent(out)  :: thresh_GF
  logical,intent(out)           :: DIIS_GF
  integer,intent(out)           :: n_diis_GF
  logical,intent(out)           :: linGF
  integer,intent(out)           :: renormGF
  double precision,intent(out)  :: eta_GF
  logical,intent(out)           :: regGF

  integer,intent(out)           :: maxSCF_GW
  double precision,intent(out)  :: thresh_GW
  logical,intent(out)           :: DIIS_GW
  integer,intent(out)           :: n_diis_GW
  logical,intent(out)           :: COHSEX
  logical,intent(out)           :: SOSEX
  logical,intent(out)           :: TDA_W
  logical,intent(out)           :: linGW
  double precision,intent(out)  :: eta_GW
  logical,intent(out)           :: regGW

  integer,intent(out)           :: maxSCF_GT
  double precision,intent(out)  :: thresh_GT
  logical,intent(out)           :: DIIS_GT
  integer,intent(out)           :: n_diis_GT
  logical,intent(out)           :: TDA_T
  logical,intent(out)           :: linGT
  double precision,intent(out)  :: eta_GT
  logical,intent(out)           :: regGT

  logical,intent(out)           :: doACFDT
  logical,intent(out)           :: exchange_kernel
  logical,intent(out)           :: doXBS

  logical,intent(out)           :: BSE
  logical,intent(out)           :: dBSE
  logical,intent(out)           :: dTDA
  logical,intent(out)           :: evDyn
  logical,intent(out)           :: ppBSE
  logical,intent(out)           :: BSE2

! Local variables

  character(len=1)              :: answer1,answer2,answer3,answer4,answer5,answer6,answer7,answer8

! Open file with method specification

  open(unit=1,file='input/options')

! Read HF options

  maxSCF_HF    = 64
  thresh_HF    = 1d-6
  DIIS_HF      = .false.
  n_diis_HF    = 5
  guess_type   = 1
  ortho_type   = 1
  mix          = .false.
  level_shift  = 0d0
  dostab       = .false.

  read(1,*) 
  read(1,*) maxSCF_HF,thresh_HF,answer1,n_diis_HF,guess_type,ortho_type,answer2,level_shift,answer3

  if(answer1 == 'T') DIIS_HF     = .true.
  if(answer2 == 'T') mix         = .true.
  if(answer3 == 'T') dostab      = .true.

  if(.not.DIIS_HF) n_diis_HF = 1

! Read MPn options

  regMP = .false.
  read(1,*)
  read(1,*) answer1

  if(answer1 == 'T') regMP = .true.

! Read CC options

  maxSCF_CC    = 64
  thresh_CC    = 1d-5
  DIIS_CC      = .false.
  n_diis_CC    = 5

  read(1,*)
  read(1,*) maxSCF_CC,thresh_CC,answer1,n_diis_CC

  if(answer1 == 'T') DIIS_CC    = .true.

  if(.not.DIIS_CC) n_diis_CC = 1

! Read excited state options

  TDA            = .false.
  singlet        = .false.
  triplet        = .false.
  spin_conserved = .false.
  spin_flip      = .false.

  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4,answer5

  if(answer1 == 'T') TDA            = .true.
  if(answer2 == 'T') singlet        = .true.
  if(answer3 == 'T') triplet        = .true.
  if(answer4 == 'T') spin_conserved = .true.
  if(answer5 == 'T') spin_flip      = .true.

! Read GF options

  maxSCF_GF = 64
  thresh_GF = 1d-5
  DIIS_GF   = .false.
  n_diis_GF = 5
  linGF     = .false.
  eta_GF    = 0d0
  renormGF  = 0
  regGF     = .false.

  read(1,*) 
  read(1,*) maxSCF_GF,thresh_GF,answer1,n_diis_GF,answer2,eta_GF,renormGF,answer3

  if(answer1 == 'T') DIIS_GF = .true.
  if(answer2 == 'T') linGF   = .true.
  if(answer3 == 'T') regGF   = .true.
  if(.not.DIIS_GF) n_diis_GF = 1

! Read GW options

  maxSCF_GW = 64
  thresh_GW = 1d-5
  DIIS_GW   = .false.
  n_diis_GW = 5
  linGW     = .false.
  eta_GW    = 0d0
  regGW     = .false.
  COHSEX    = .false.
  SOSEX     = .false.
  TDA_W     = .false.

  read(1,*) 
  read(1,*) maxSCF_GW,thresh_GW,answer1,n_diis_GW,answer2,eta_GW, &
            answer3,answer4,answer5,answer6

  if(answer1 == 'T') DIIS_GW = .true.
  if(answer2 == 'T') linGW   = .true.
  if(answer3 == 'T') COHSEX  = .true.
  if(answer4 == 'T') SOSEX   = .true.
  if(answer5 == 'T') TDA_W   = .true.
  if(answer6 == 'T') regGW   = .true.
  if(.not.DIIS_GW) n_diis_GW = 1

! Read GT options

  maxSCF_GT = 64
  thresh_GT = 1d-5
  DIIS_GT   = .false.
  n_diis_GT = 5
  linGT     = .false.
  eta_GT    = 0d0
  regGT     = .false.
  TDA_T     = .false.

  read(1,*) 
  read(1,*) maxSCF_GT,thresh_GT,answer1,n_diis_GT,answer2,eta_GT, &
            answer3,answer4

  if(answer1 == 'T') DIIS_GT = .true.
  if(answer2 == 'T') linGT   = .true.
  if(answer3 == 'T') TDA_T   = .true.
  if(answer4 == 'T') regGT   = .true.
  if(.not.DIIS_GT) n_diis_GT = 1

! Options for adiabatic connection

  doACFDT         = .false.
  exchange_kernel = .false.
  doXBS           = .false.

  read(1,*) 
  read(1,*) answer1,answer2,answer3

  if(answer1 == 'T') doACFDT = .true.
  if(answer2 == 'T') exchange_kernel = .true.
  if(answer3 == 'T') doXBS   = .true.

! Options for dynamical BSE calculations

  BSE   = .false.
  dBSE  = .false.
  dTDA  = .true.
  evDyn = .false.
  ppBSE = .false.
  BSE2  = .false.

  read(1,*) 
  read(1,*) answer1,answer2,answer3,answer4,answer5,answer6

  if(answer1 == 'T') BSE   = .true.
  if(answer2 == 'T') dBSE  = .true.
  if(answer3 == 'F') dTDA  = .false.
  if(answer4 == 'T') evDyn = .true.
  if(answer5 == 'T') ppBSE = .true.
  if(answer6 == 'T') BSE2  = .true.

! Close file with options

  close(unit=1)

end subroutine read_options
