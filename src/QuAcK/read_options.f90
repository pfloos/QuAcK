subroutine read_options(working_dir,                                                                         &
                        maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,dostab,dosearch,doaordm,  &
                        readFCIDUMP,doMOM,writeMOs,reg_MP,                                                   &
                        maxSCF_CC,thresh_CC,max_diis_CC,                                                     &
                        TDA,spin_conserved,spin_flip,nCVS,FC,                                                &
                        max_iter_OO,thresh_OO,dRPA_OO,mu_OO,diagHess_OO,                                     &
                        maxSCF_GF,thresh_GF,max_diis_GF,lin_GF,eta_GF,renorm_GF,reg_GF,do_linDM_GF2,         &
                        maxSCF_GW,thresh_GW,max_diis_GW,lin_GW,eta_GW,shift_GW,reg_GW,do_linDM_GW,dosign_XoB,&
                        nfreqs,TDA_W,restart_scGW,restart_scGF2,verbose_scGW,verbose_scGF2,chem_pot_scG,     &
                        maxSCF_GT,thresh_GT,max_diis_GT,lin_GT,eta_GT,reg_GT,TDA_T,do_linDM_GT,              &
                        doACFDT,exchange_kernel,doXBS,                                                       &
                        dophBSE,dophBSE2,doppBSE,dBSE,dTDA,                                                  &
                        temperature,sigma,chem_pot_hf,restart_hfb,error_P,                                   &
                        TDAeh,TDApp,max_diis_1b,max_diis_2b,max_it_1b,conv_1b,max_it_2b,conv_2b,lin_parquet, &
                        reg_1b,reg_2b,reg_PA,                                                                &
                  parquet_method,faddeev_coupling,faddeev_davidson,n_IP_roots,faddeev_static3,gw_faddeev,gw_external, &
                        do_dyson,diag_approx,sig_inf,lin_ADC,reg_ADC,eta_ADC,                                &
                        eweight,eforward)

! Read desired methods 

  implicit none
  include 'parameters.h'

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
  logical,intent(out)           :: doaordm
  logical,intent(out)           :: readFCIDUMP
  logical,intent(out)           :: doMOM
  logical,intent(out)           :: writeMOs
  logical,intent(out)           :: error_P

  logical,intent(out)           :: reg_MP

  integer,intent(out)           :: maxSCF_CC
  double precision,intent(out)  :: thresh_CC
  integer,intent(out)           :: max_diis_CC

  logical,intent(out)           :: TDA
  logical,intent(out)           :: spin_conserved
  logical,intent(out)           :: spin_flip

  integer,intent(out)           :: nCVS(nspin)
  integer,intent(out)           :: FC(nspin)

  integer,intent(out)           :: max_iter_OO
  double precision,intent(out)  :: thresh_OO
  logical,intent(out)           :: dRPA_OO
  integer,intent(out)           :: mu_OO
  logical,intent(out)           :: diagHess_OO
  
  integer,intent(out)           :: maxSCF_GF
  double precision,intent(out)  :: thresh_GF
  integer,intent(out)           :: max_diis_GF
  logical,intent(out)           :: lin_GF
  integer,intent(out)           :: renorm_GF
  double precision,intent(out)  :: eta_GF
  logical,intent(out)           :: reg_GF
  logical,intent(out)           :: restart_scGF2
  logical,intent(out)           :: verbose_scGF2
  logical,intent(out)           :: do_linDM_GF2

  integer,intent(out)           :: maxSCF_GW
  double precision,intent(out)  :: thresh_GW
  integer,intent(out)           :: max_diis_GW
  logical,intent(out)           :: TDA_W
  logical,intent(out)           :: lin_GW
  double precision,intent(out)  :: eta_GW
  double precision,intent(out)  :: shift_GW
  logical,intent(out)           :: restart_scGW
  logical,intent(out)           :: verbose_scGW
  logical,intent(out)           :: chem_pot_scG
  logical,intent(out)           :: reg_GW
  integer,intent(out)           :: nfreqs
  logical,intent(out)           :: do_linDM_GW
  logical,intent(out)           :: dosign_XoB

  integer,intent(out)           :: maxSCF_GT
  double precision,intent(out)  :: thresh_GT
  integer,intent(out)           :: max_diis_GT
  logical,intent(out)           :: TDA_T
  logical,intent(out)           :: lin_GT
  double precision,intent(out)  :: eta_GT
  logical,intent(out)           :: reg_GT
  logical,intent(out)           :: do_linDM_GT

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
  logical,intent(out)           :: lin_parquet,reg_PA

  integer,intent(out) :: parquet_method,faddeev_coupling,n_IP_roots
  logical,intent(out) :: faddeev_davidson,faddeev_static3,gw_faddeev,gw_external

  logical,intent(out)           :: do_dyson,diag_approx,sig_inf,lin_ADC,reg_ADC
  double precision,intent(out)  :: eta_ADC

  logical,intent(out)           :: eforward
  double precision,intent(out)  :: eweight

! Local variables

  character(len=1)              :: ans1,ans2,ans3,ans4,ans5,ans6,ans7,ans8
  integer                       :: status
  character(len=256)            :: file_path,option_line

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
      doaordm      = .false.
      readFCIDUMP  = .false.
      doMOM        = .false.
      writeMOs     = .false.
    
      read(1,*) 
      read(1,*) maxSCF_HF,thresh_HF,max_diis_HF,guess_type,mix,level_shift,ans1,ans2,ans3,ans4,ans5,ans6
    
      if(ans1 == 'T') dostab   = .true.
      if(ans2 == 'T') dosearch = .true.
      if(ans3 == 'T') doaordm  = .true.
      if(ans4 == 'T') readFCIDUMP  = .true.
      if(ans5 == 'T') doMOM    = .true.
      if(ans6 == 'T') writeMOs = .true.
    
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
    
      ! Number of inactive virtuals (CVS)
    
      read(1,*) 
      read(1,*) nCVS(1),nCVS(2), FC(1),FC(2)
    
      ! Read Orbital optimization options

      max_iter_OO = 256
      thresh_OO   = 1e-5
      dRPA_OO  = .false.
      mu_OO    = 0
      diagHess_OO = .false.

      read(1,*) 
      read(1,*) max_iter_OO,thresh_OO,ans1,mu_OO,ans2
      if(ans1 == 'T') dRPA_OO     = .true.
      if(ans2 == 'T') diagHess_OO = .true.
      

      ! Read GF options
    
      maxSCF_GF   = 64
      thresh_GF   = 1d-5
      max_diis_GF = 1
      lin_GF      = .false.
      eta_GF      = 0d0
      renorm_GF   = 0
      reg_GF      = .false.
      do_linDM_GF2 = .false.
      restart_scGF2 = .false.
      verbose_scGF2 = .false.
    
      read(1,*) 
      read(1,*) maxSCF_GF,thresh_GF,max_diis_GF,ans1,eta_GF,renorm_GF,ans2,ans3,ans4,ans5
    
      if(ans1 == 'T') lin_GF  = .true.
      if(ans2 == 'T') reg_GF  = .true.
      if(ans3 == 'T') do_linDM_GF2  = .true.
      if(ans4 == 'T') restart_scGF2  = .true.
      if(ans5 == 'T') verbose_scGF2  = .true.
    
      ! Read GW options
    
      maxSCF_GW   = 64
      thresh_GW   = 1d-5
      max_diis_GW = 1
      lin_GW      = .false.
      eta_GW      = 0d0
      shift_GW    = 1d-3
      reg_GW      = .false.
      TDA_W       = .false.
      do_linDM_GW = .false.
      dosign_XoB  = .false.
      restart_scGW = .false.
      verbose_scGW = .false.
      chem_pot_scG = .false.
    
      read(1,*) 
      read(1,*) maxSCF_GW,thresh_GW,max_diis_GW,ans1,eta_GW,ans2,ans3,nfreqs,shift_GW,ans4,ans5,ans6,ans7,ans8
    
      if(ans1 == 'T') lin_GW       = .true.
      if(ans2 == 'T') TDA_W        = .true.
      if(ans3 == 'T') reg_GW       = .true.
      if(ans4 == 'T') do_linDM_GW  = .true.
      if(ans5 == 'T') restart_scGW = .true.
      if(ans6 == 'T') verbose_scGW = .true.
      if(ans7 == 'T') chem_pot_scG = .true.
      if(ans8 == 'T') dosign_XoB   = .true.
   
      ! Read GT options
    
      maxSCF_GT   = 64
      thresh_GT   = 1d-5
      max_diis_GT = 1
      lin_GT      = .false.
      eta_GT      = 0d0
      reg_GT      = .false.
      TDA_T       = .false.
      do_linDM_GT = .false.
    
      read(1,*) 
      read(1,*) maxSCF_GT,thresh_GT,max_diis_GT,ans1,eta_GT,ans2,ans3,ans4
    
      if(ans1 == 'T') lin_GT      = .true.
      if(ans2 == 'T') TDA_T       = .true.
      if(ans3 == 'T') reg_GT      = .true.
      if(ans4 == 'T') do_linDM_GT = .true.
    
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
      error_P      = .false.
    
      read(1,*) 
      read(1,*) temperature,sigma,ans1,ans2,ans3

      if(ans1 == 'T') chem_pot_hf  = .true.
      if(ans2 == 'T') restart_hfb  = .true.
      if(ans3 == 'T') error_P      = .true.

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
      reg_PA      = .false.
    
      read(1,*) 
      read(1,*) ans1,ans2,max_it_1b,conv_1b,max_it_2b,conv_2b,max_diis_1b,max_diis_2b,ans3,reg_1b,reg_2b,ans4

      if(ans1 == 'T') TDAeh = .true.
      if(ans2 == 'T') TDApp = .true.
      if(ans3 == 'T') lin_parquet = .true.
      if(ans4 == 'T') reg_PA = .true.

      ! Options for ADC module

      do_dyson    = .false.
      diag_approx = .false.
      sig_inf     = .false.
      lin_ADC     = .false.
      reg_ADC     = .false.
      eta_ADC     = 0d0

      read(1,*)
      read(1,*) ans1,ans2,ans3,ans4,ans5,eta_ADC

      if(ans1 == 'T') do_dyson    = .true.
      if(ans2 == 'T') diag_approx = .true.
      if(ans3 == 'T') sig_inf     = .true.
      if(ans4 == 'T') lin_ADC     = .true.
      if(ans5 == 'T') reg_ADC     = .true.
      
      ! Read ensemble HF options
    
      eforward    = .false.
      eweight     = 0d0

      read(1,*)
      read(1,*) eweight,ans1

      if(ans1 == 'T') eforward = .true.
 
      ! Optional trailing section for evParquet self-energy selection.
      parquet_method = 2
      faddeev_coupling = 0
      faddeev_davidson = .true.
      n_IP_roots = 1
      faddeev_static3 = .false.
      gw_faddeev = .false.
      gw_external = .false.
      do
        read(1,'(A)',iostat=status) option_line
        if(status < 0) exit
        if(status /= 0) error stop 'Error reading Parquet self-energy options'
        if(index(adjustl(option_line),'# GWParquet Faddeev external GW:') == 1) then
          read(1,'(A)',iostat=status) option_line
          if(status /= 0) error stop 'Missing GWParquet Faddeev external GW flag'
          ans1 = ' '
          read(option_line,*,iostat=status) ans1
          if(status /= 0.or.(ans1 /= 'T'.and.ans1 /= 'F')) error stop 'External GW flag must be T or F'
          gw_external = ans1 == 'T'
          cycle
        end if
        if(index(adjustl(option_line),'# GWParquet Faddeev:') == 1) then
          read(1,'(A)',iostat=status) option_line
          if(status /= 0) error stop 'Missing GWParquet Faddeev flag'
          ans1 = ' '
          read(option_line,*,iostat=status) ans1
          if(status /= 0.or.(ans1 /= 'T'.and.ans1 /= 'F')) error stop 'GWParquet Faddeev flag must be T or F'
          gw_faddeev = ans1 == 'T'
          cycle
        end if
        if(index(adjustl(option_line),'# Faddeev static self-energy:') == 1) then
          read(1,'(A)',iostat=status) option_line
          if(status /= 0) error stop 'Missing Faddeev static self-energy flag'
          ans1 = ' '
          read(option_line,*,iostat=status) ans1
          if(status /= 0.or.(ans1 /= 'T'.and.ans1 /= 'F')) error stop 'Faddeev static self-energy flag must be T or F'
          faddeev_static3 = ans1 == 'T'
          cycle
        end if
        if(index(adjustl(option_line),'# Parquet self-energy:') /= 1) cycle
        read(1,'(A)',iostat=status) option_line
        if(status /= 0) error stop 'Missing Parquet self-energy option values'
        ans1 = " "
        read(option_line,*,iostat=status) parquet_method,faddeev_coupling,ans1,n_IP_roots
        if(status /= 0) error stop 'Expected Parquet method, Faddeev coupling, Davidson flag and IP root count'
        if(parquet_method < 0.or.parquet_method > 2) error stop 'Parquet method must be 0, 1 or 2'
        if(faddeev_coupling /= 0.and.(faddeev_coupling < 2.or.faddeev_coupling > 4)) &
          error stop 'Faddeev coupling must be 0, 2, 3 or 4'
        if(ans1 /= 'T'.and.ans1 /= 'F') error stop 'Faddeev Davidson flag must be T or F'
        faddeev_davidson = ans1 == 'T'
        if(n_IP_roots < 1) error stop 'Faddeev IP root count must be positive'
      end do

    endif

  ! Close file with options

  close(unit=1)

end subroutine 
