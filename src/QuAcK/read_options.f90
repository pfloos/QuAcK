subroutine read_options(maxSCF_HF,thresh_HF,DIIS_HF,n_diis_HF,guess_type,ortho_type,   &
                        maxSCF_CC,thresh_CC,DIIS_CC,n_diis_CC,                         &
                        singlet_manifold,triplet_manifold,                             &
                        maxSCF_GF,thresh_GF,DIIS_GF,n_diis_GF,linGF,renormGF,          &
                        BSE_GF,TDA_GF,eta_GF,                                          &
                        maxSCF_GW,thresh_GW,DIIS_GW,n_diis_GW,                         &
                        COHSEX,SOSEX,BSE_GW,TDA_W,TDA_GW,G0W,GW0,linGW,eta_GW,         &
                        doACFDT,exchange_kernel,doXBS,                                 &
                        nMC,nEq,nWalk,dt,nPrint,iSeed,doDrift)

! Read desired methods 

  implicit none

! Input variables

  integer,intent(out)           :: maxSCF_HF
  double precision,intent(out)  :: thresh_HF
  logical,intent(out)           :: DIIS_HF
  integer,intent(out)           :: n_diis_HF
  integer,intent(out)           :: guess_type
  integer,intent(out)           :: ortho_type

  integer,intent(out)           :: maxSCF_CC
  double precision,intent(out)  :: thresh_CC
  logical,intent(out)           :: DIIS_CC
  integer,intent(out)           :: n_diis_CC

  logical,intent(out)           :: singlet_manifold
  logical,intent(out)           :: triplet_manifold

  integer,intent(out)           :: maxSCF_GF
  double precision,intent(out)  :: thresh_GF
  logical,intent(out)           :: DIIS_GF
  integer,intent(out)           :: n_diis_GF
  logical,intent(out)           :: linGF
  integer,intent(out)           :: renormGF
  logical,intent(out)           :: BSE_GF
  logical,intent(out)           :: TDA_GF
  double precision,intent(out)  :: eta_GF

  integer,intent(out)           :: maxSCF_GW
  double precision,intent(out)  :: thresh_GW
  logical,intent(out)           :: DIIS_GW
  integer,intent(out)           :: n_diis_GW
  logical,intent(out)           :: COHSEX
  logical,intent(out)           :: SOSEX
  logical,intent(out)           :: BSE_GW
  logical,intent(out)           :: TDA_W
  logical,intent(out)           :: TDA_GW
  logical,intent(out)           :: G0W
  logical,intent(out)           :: GW0
  logical,intent(out)           :: linGW
  double precision,intent(out)  :: eta_GW

  logical,intent(out)           :: doACFDT
  logical,intent(out)           :: exchange_kernel
  logical,intent(out)           :: doXBS

  integer,intent(out)           :: nMC
  integer,intent(out)           :: nEq
  integer,intent(out)           :: nWalk
  double precision,intent(out)  :: dt
  integer,intent(out)           :: nPrint
  integer,intent(out)           :: iSeed
  logical,intent(out)           :: doDrift
  
! Local variables

  character(len=1)              :: answer1,answer2,answer3,answer4,answer5,answer6,answer7,answer8,answer9

! Open file with method specification

  open(unit=1,file='input/options')

! Read HF options

  maxSCF_HF    = 64
  thresh_HF    = 1d-6
  DIIS_HF      = .false.
  n_diis_HF    = 5
  guess_type   = 1
  ortho_type   = 1

  read(1,*) 
  read(1,*) maxSCF_HF,thresh_HF,answer1,n_diis_HF,guess_type,ortho_type

  if(answer1 == 'T') DIIS_HF    = .true.

  if(.not.DIIS_HF) n_diis_HF = 1

! Read MPn options

  read(1,*) 
  read(1,*)

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

  singlet_manifold = .false.
  triplet_manifold = .false.

  read(1,*) 
  read(1,*) answer1,answer2

  if(answer1 == 'T') singlet_manifold = .true.
  if(answer2 == 'T') triplet_manifold = .true.

! Read Green function options

  maxSCF_GF = 64
  thresh_GF = 1d-5
  DIIS_GF   = .false.
  n_diis_GF = 5
  linGF     = .false.
  renormGF  = 0
  BSE_GF    = .false.
  TDA_GF    = .false.
  eta_GF    = 0d0

  read(1,*) 
  read(1,*) maxSCF_GF,thresh_GW,answer1,n_diis_GF,answer2,renormGF,answer3,answer4,eta_GF

  if(answer1 == 'T') DIIS_GF = .true.
  if(answer2 == 'T') linGF   = .true.
  if(answer3 == 'T') BSE_GF  = .true.
  if(answer4 == 'T') TDA_GF  = .true.
  if(.not.DIIS_GF) n_diis_GF = 1

! Read GW options

  maxSCF_GW = 64
  thresh_GW = 1d-5
  DIIS_GW   = .false.
  n_diis_GW = 5
  COHSEX    = .false.
  SOSEX     = .false.
  BSE_GW    = .false.
  TDA_W     = .false.
  TDA_GW    = .false.
  G0W       = .false.
  GW0       = .false.
  linGW     = .false.
  eta_GW    = 0d0

  read(1,*) 
  read(1,*) maxSCF_GW,thresh_GW,answer1,n_diis_GW,answer2, &
            answer3,answer4,answer5,answer6,answer7,answer8,answer9,eta_GW

  if(answer1 == 'T') DIIS_GW = .true.
  if(answer2 == 'T') COHSEX  = .true.
  if(answer3 == 'T') SOSEX   = .true.
  if(answer4 == 'T') BSE_GW  = .true.
  if(answer5 == 'T') TDA_W   = .true.
  if(answer6 == 'T') TDA_GW  = .true.
  if(answer7 == 'T') G0W     = .true.
  if(answer8 == 'T') GW0     = .true.
  if(answer9 == 'T') linGW   = .true.
  if(.not.DIIS_GW) n_diis_GW = 1

! Options for adiabatic connection

  doACFDT = .false.
  exchange_kernel = .false.
  doXBS   = .false.

  read(1,*) 
  read(1,*) answer1,answer2,answer3

  if(answer1 == 'T') doACFDT = .true.
  if(answer2 == 'T') exchange_kernel = .true.
  if(answer3 == 'T') doXBS   = .true.

! Read options for MC-MP2: Monte Carlo steps, number of equilibration steps, number of walkers,
! Monte Carlo time step, frequency of output results, and seed for random number generator

  nMC    = 100000
  nEq    = 10000
  nWalk  = 10
  dt     = 0.3d0
  nPrint = 1000
  iSeed  = 0
  doDrift = .false.

  read(1,*)
  read(1,*) nMC,nEq,nWalk,dt,nPrint,iSeed,answer1

  if(answer1 == 'T') doDrift = .true.

! Close file with options

  close(unit=1)

end subroutine read_options
