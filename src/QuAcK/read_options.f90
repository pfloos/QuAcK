subroutine read_options(working_dir,                                                                        &
                        maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,dostab,dosearch,         &
                        reg_MP,                                                                             &
                        maxSCF_CC,thresh_CC,max_diis_CC,                                                    &
                        TDA,spin_conserved,spin_flip,                                                       &
                        maxSCF_GF,thresh_GF,max_diis_GF,lin_GF,eta_GF,renorm_GF,reg_GF,                     &
                        maxSCF_GW,thresh_GW,max_diis_GW,lin_GW,eta_GW,reg_GW,doOO,mu,                       &
                        im_freqs,nfreqs,ntimes,TDA_W,                                                       &
                        maxSCF_GT,thresh_GT,max_diis_GT,lin_GT,eta_GT,reg_GT,TDA_T,                         &
                        doACFDT,exchange_kernel,doXBS,                                                      &
                        dophBSE,dophBSE2,doppBSE,dBSE,dTDA,                                                 &
                        temperature,sigma,chem_pot_hf,restart_hfb,                                          &
                        TDAeh,TDApp,max_diis_1b,max_diis_2b,max_it_1b,conv_1b,max_it_2b,conv_2b,lin_parquet,reg_1b,reg_2b)

! Read desired methods 

  implicit none

! Input variables

  character(len=256),intent(in) :: working_dir

! Output variables

  integer,intent(out)           :: maxSCF_HF
  double precision,intent(out)  :: thresh_HF
  integer,intent(out)           :: max_diis_HF
  integer,intent(out)           :: guess_type
  double precision,intent(out)  :: mix
  double precision,intent(out)  :: level_shift
  logical,intent(out)           :: dostab
  logical,intent(out)           :: dosearch

  logical,intent(out)           :: reg_MP

  integer,intent(out)           :: maxSCF_CC
  double precision,intent(out)  :: thresh_CC
  integer,intent(out)           :: max_diis_CC

  logical,intent(out)           :: TDA
  logical,intent(out)           :: spin_conserved
  logical,intent(out)           :: spin_flip

  integer,intent(out)           :: maxSCF_GF
  double precision,intent(out)  :: thresh_GF
  integer,intent(out)           :: max_diis_GF
  logical,intent(out)           :: lin_GF
  integer,intent(out)           :: renorm_GF
  double precision,intent(out)  :: eta_GF
  logical,intent(out)           :: reg_GF

  integer,intent(out)           :: maxSCF_GW
  double precision,intent(out)  :: thresh_GW
  integer,intent(out)           :: max_diis_GW
  logical,intent(out)           :: TDA_W
  logical,intent(out)           :: lin_GW
  double precision,intent(out)  :: eta_GW
  logical,intent(out)           :: reg_GW
  logical,intent(out)           :: doOO
  integer,intent(out)           :: mu
  logical,intent(out)           :: im_freqs
  integer,intent(out)           :: nfreqs
  integer,intent(out)           :: ntimes

  integer,intent(out)           :: maxSCF_GT
  double precision,intent(out)  :: thresh_GT
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

  logical,intent(out)           :: chem_pot_hf
  logical,intent(out)           :: restart_hfb
  double precision,intent(out)  :: temperature
  double precision,intent(out)  :: sigma

  integer,intent(out)           :: max_it_1b,max_it_2b
  double precision,intent(out)  :: conv_1b,conv_2b
  integer,intent(out)           :: max_diis_1b,max_diis_2b
  logical,intent(out)           :: TDAeh,TDApp
  double precision,intent(out)  :: reg_1b,reg_2b
  logical,intent(out)           :: lin_parquet

! Local variables

  character(len=1)              :: ans1,ans2,ans3,ans4,ans5
  integer                       :: status
  character(len=256)            :: file_path

! Open file with method specification

  file_path = trim(working_dir) // '/input/options'
  open(unit=1, file=file_path, status='old', action='read', iostat=status)

    if(status /= 0) then

      print *, "Error opening file: ", file_path
      stop

    else

      ! Read HF options
    
      maxSCF_HF    = 64
      thresh_HF    = 1d-6
      max_diis_HF  = 1
      guess_type   = 1
      mix          = 0d0
      level_shift  = 0d0
      dostab       = .false.
      dosearch     = .false.
    
      read(1,*) 
      read(1,*) maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,ans1,ans2
    
      if(ans1 == 'T') dostab   = .true.
      if(ans2 == 'T') dosearch = .true.
    
      ! Read MPn options
    
      reg_MP = .false.
      read(1,*)
      read(1,*) ans1
    
      if(ans1 == 'T') reg_MP = .true.
    
      ! Read CC options
    
      maxSCF_CC    = 64
      thresh_CC    = 1d-5
      max_diis_CC  = 1
    
      read(1,*)
      read(1,*) maxSCF_CC,thresh_CC,max_diis_CC
    
      ! Read excited state options
    
      TDA            = .false.
      spin_conserved = .false.
      spin_flip      = .false.
    
      read(1,*) 
      read(1,*) ans1,ans2,ans3
    
      if(ans1 == 'T') TDA            = .true.
      if(ans2 == 'T') spin_conserved = .true.
      if(ans3 == 'T') spin_flip      = .true.
    
      ! Read GF options
    
      maxSCF_GF   = 64
      thresh_GF   = 1d-5
      max_diis_GF = 1
      lin_GF      = .false.
      eta_GF      = 0d0
      renorm_GF   = 0
      reg_GF      = .false.
    
      read(1,*) 
      read(1,*) maxSCF_GF,thresh_GF,max_diis_GF,ans1,eta_GF,renorm_GF,ans2
    
      if(ans1 == 'T') lin_GF  = .true.
      if(ans2 == 'T') reg_GF  = .true.
    
      ! Read GW options
    
      maxSCF_GW   = 64
      thresh_GW   = 1d-5
      max_diis_GW = 1
      lin_GW      = .false.
      eta_GW      = 0d0
      reg_GW      = .false.
      doOO        = .false.
      mu          = 0
      TDA_W       = .false.
      im_freqs    = .false.
    
      read(1,*) 
      read(1,*) maxSCF_GW,thresh_GW,max_diis_GW,ans1,eta_GW,ans2,ans3,ans4,mu,ans5,nfreqs,ntimes
    
      if(ans1 == 'T') lin_GW   = .true.
      if(ans2 == 'T') TDA_W    = .true.
      if(ans3 == 'T') reg_GW   = .true.
      if(ans4 == 'T') doOO     = .true. 
      if(ans5 == 'T') im_freqs = .true.
   
      ! Read GT options
    
      maxSCF_GT   = 64
      thresh_GT   = 1d-5
      max_diis_GT = 1
      lin_GT      = .false.
      eta_GT      = 0d0
      reg_GT      = .false.
      TDA_T       = .false.
    
      read(1,*) 
      read(1,*) maxSCF_GT,thresh_GT,max_diis_GT,ans1,eta_GT,ans2,ans3
    
      if(ans1 == 'T') lin_GT  = .true.
      if(ans2 == 'T') TDA_T   = .true.
      if(ans3 == 'T') reg_GT  = .true.
    
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

      ! Options for Hartree-Fock Bogoliubov
    
      temperature  = 0d0
      sigma        = 1d0
      chem_pot_hf  = .false.
      restart_hfb  = .false.
    
      read(1,*) 
      read(1,*) temperature,sigma,ans1,ans2

      if(ans1 == 'T') chem_pot_hf  = .true.
      if(ans2 == 'T') restart_hfb  = .true.

      ! Options for Parquet module

      TDAeh       = .false.
      TDApp       = .false.
      max_diis_1b = 1
      max_diis_2b = 1
      max_it_1b   = 1
      conv_1b     = 1d-2
      max_it_2b   = 1
      conv_2b     = 1d-2
      lin_parquet = .false.
      reg_1b      = 0d0
      reg_2b      = 0d0
    
      read(1,*) 
      read(1,*) ans1,ans2,max_it_1b,conv_1b,max_it_2b,conv_2b,max_diis_1b,max_diis_2b,ans3,reg_1b,reg_2b

      if(ans1 == 'T') TDAeh = .true.
      if(ans2 == 'T') TDApp = .true.
      if(ans3 == 'T') lin_parquet = .true.
 
    endif

  ! Close file with options

  close(unit=1)

end subroutine 
