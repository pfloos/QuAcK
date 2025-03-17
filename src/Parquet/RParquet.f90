subroutine RParquet(max_it_macro,conv_one_body,max_it_micro,conv_two_body,nOrb,nC,nO,nV,nR,nS,eHF,ERI)

! Spatial orbital Parquet implementation
  implicit none
  include 'parameters.h'

! Hard-coded parameters
  logical                       :: linearize = .true.
  logical                       :: TDA = .true.
  logical                       :: print_phLR = .false.
  logical                       :: print_ppLR = .false.
  
! Input variables
  integer,intent(in)            :: max_it_macro,max_it_micro
  double precision,intent(in)   :: conv_one_body,conv_two_body
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  
! Local variables
  integer                       :: n_it_macro,n_it_micro
  double precision              :: err_one_body,err_two_body
  double precision              :: err_eh_sing,err_eh_trip
  double precision              :: err_hh_sing,err_hh_trip
  double precision              :: err_ee_sing,err_ee_trip
  double precision              :: start_t, end_t, t

  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  
  double precision              :: EcRPA
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: sing_XpY(:,:),trip_XpY(:,:)
  double precision,allocatable  :: sing_XmY(:,:),trip_XmY(:,:)
  double precision,allocatable  :: eh_sing_Om(:), old_eh_sing_Om(:)
  double precision,allocatable  :: eh_trip_Om(:), old_eh_trip_Om(:)

  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: ee_sing_Om(:), old_ee_sing_Om(:)
  double precision,allocatable  :: ee_trip_Om(:), old_ee_trip_Om(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: hh_sing_Om(:), old_hh_sing_Om(:)
  double precision,allocatable  :: hh_trip_Om(:), old_hh_trip_Om(:)

  double precision,allocatable  :: eh_sing_rho(:,:,:),eh_trip_rho(:,:,:)
  double precision,allocatable  :: ee_sing_rho(:,:,:),hh_sing_rho(:,:,:)
  double precision,allocatable  :: ee_trip_rho(:,:,:),hh_trip_rho(:,:,:)

  double precision,allocatable  :: eh_sing_Gam_A(:,:),eh_sing_Gam_B(:,:)
  double precision,allocatable  :: eh_trip_Gam_A(:,:),eh_trip_Gam_B(:,:)
  double precision,allocatable  :: pp_sing_Gam_B(:,:),pp_sing_Gam_C(:,:),pp_sing_Gam_D(:,:)
  double precision,allocatable  :: pp_trip_Gam_B(:,:),pp_trip_Gam_C(:,:),pp_trip_Gam_D(:,:)
  double precision,allocatable  :: eh_sing_Gam(:,:,:,:),eh_trip_Gam(:,:,:,:)
  double precision,allocatable  :: pp_sing_Gam(:,:,:,:),pp_trip_Gam(:,:,:,:)

  double precision,allocatable  :: eParquetlin(:),eParquet(:),old_eParquet(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision              :: EcGM
  
! Output variables
  nOOs = nO*(nO + 1)/2
  nVVs = nV*(nV + 1)/2
  nOOt = nO*(nO - 1)/2
  nVVt = nV*(nV - 1)/2

  allocate(eParquet(nOrb),old_eParquet(nOrb))
  
  write(*,*)
  write(*,*)'**********************************'
  write(*,*)'* Restricted Parquet Calculation *'
  write(*,*)'**********************************'
  write(*,*)

! Print parameters
  write(*,*)'Parameters for this run:'
  write(*,*)'Maximum number of macro iteration:', max_it_macro
  write(*,*)'Convergence threshold for one-body energies:', conv_one_body
  write(*,*)'Maximum number of micro iteration:', max_it_micro
  write(*,*)'Convergence threshold for two-body energies:', conv_two_body
  
  if (linearize) then 
      write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
      write(*,*)
   else 
      write(*,*) ' *** Quasiparticle energies obtained by root search *** '
      write(*,*)
  endif

! Initialization
  n_it_macro      = 1
  err_one_body    = 1d0
  n_it_micro      = 1
  err_two_body    = 1d0
  old_eParquet(:) = eHF(:)
  
  write(*,*)
  write(*,*)'************ Solving initial linear-response problems ************'
  write(*,*)'------------------------------------------------------------------'

  !-------------------
  ! Density channel
  !-------------------
  allocate(Aph(nS,nS),Bph(nS,nS),eh_sing_Om(nS),sing_XpY(nS,nS),sing_XmY(nS,nS),old_eh_sing_Om(nS))
  Aph(:,:) = 0d0
  Bph(:,:) = 0d0
  call wall_time(start_t)
  if(.not.TDA) call phRLR_B(1,.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)
               call phRLR_A(1,.false.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  call phRLR(TDA,nS,Aph,Bph,EcRPA,eh_sing_Om,sing_XpY,sing_XmY)
  call wall_time(end_t)
  t = end_t - start_t
  write(*,*)
  write(*,'(A51,1X,F9.3,A8)') 'Total wall time for initial singlet phRPA problem =',t,' seconds'
  !if (print_phLR) call print_excitation_energies('phRPA@RHF','singlet',nS,eh_sing_Om)
  call print_excitation_energies('phRPA@RHF','singlet',nS,eh_sing_Om)
  
  deallocate(Aph,Bph)
  
  !-------------------
  ! Magnetic channel
  !-------------------
  allocate(Aph(nS,nS),Bph(nS,nS),eh_trip_Om(nS),trip_XpY(nS,nS),trip_XmY(nS,nS),old_eh_trip_Om(nS))
  Aph(:,:) = 0d0
  Bph(:,:) = 0d0
  call wall_time(start_t)
  if(.not.TDA) call phRLR_B(2,.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)
               call phRLR_A(2,.false.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  call phRLR(TDA,nS,Aph,Bph,EcRPA,eh_trip_Om,trip_XpY,trip_XmY)
  call wall_time(end_t)
  t = end_t - start_t
  write(*,*)
  write(*,'(A51,1X,F9.3,A8)') 'Total wall time for initial triplet phRPA problem =',t,' seconds'
  !if (print_phLR) call print_excitation_energies('phRPA@RHF','triplet',nS,eh_trip_Om)
  call print_excitation_energies('phRPA@RHF','triplet',nS,eh_trip_Om)
  
  deallocate(Aph,Bph)

  !-------------------
  ! Singlet channel
  !-------------------
  allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs),                        &
           ee_sing_Om(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),old_ee_sing_Om(nVVs), & 
           hh_sing_Om(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs),old_hh_sing_Om(nOOs))
  Bpp(:,:)=0d0
  Cpp(:,:)=0d0
  Dpp(:,:)=0d0
  call wall_time(start_t)
  if(.not.TDA) call ppRLR_B(1,nOrb,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
               call ppRLR_C(1,nOrb,nC,nO,nV,nR,nVVs,1d0,eHF,ERI,Cpp)
               call ppRLR_D(1,nOrb,nC,nO,nV,nR,nOOs,1d0,eHF,ERI,Dpp)
  call ppRLR(TDA,nOOs,nVVs,Bpp,Cpp,Dpp,ee_sing_Om,X1s,Y1s,hh_sing_Om,X2s,Y2s,EcRPA)
  call wall_time(end_t)
  t = end_t - start_t
  write(*,*)
  write(*,'(A51,1X,F9.3,A8)') 'Total wall time for initial singlet ppRPA problem =',t,' seconds'
  if (print_ppLR) call print_excitation_energies('ppRPA@RHF','2p (singlets)',nVVs,ee_sing_Om)
  if (print_ppLR) call print_excitation_energies('ppRPA@RHF','2h (singlets)',nOOs,hh_sing_Om)
  
  deallocate(Bpp,Cpp,Dpp)

  !-------------------
  ! Triplet channel
  !-------------------
  allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt),                        &
           ee_trip_Om(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),old_ee_trip_Om(nVVt), & 
           hh_trip_Om(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt),old_hh_trip_Om(nOOt))
  Bpp(:,:)=0d0
  Cpp(:,:)=0d0
  Dpp(:,:)=0d0
  call wall_time(start_t)
  if(.not.TDA) call ppRLR_B(2,nOrb,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
               call ppRLR_C(2,nOrb,nC,nO,nV,nR,nVVt,1d0,eHF,ERI,Cpp)
               call ppRLR_D(2,nOrb,nC,nO,nV,nR,nOOt,1d0,eHF,ERI,Dpp)
  call ppRLR(TDA,nOOt,nVVt,Bpp,Cpp,Dpp,ee_trip_Om,X1t,Y1t,hh_trip_Om,X2t,Y2t,EcRPA)
  call wall_time(end_t)
  t = end_t - start_t
  write(*,*)
  write(*,'(A51,1X,F9.3,A8)') 'Total wall time for initial triplet ppRPA problem =',t,' seconds'
  write(*,*)
  if(print_ppLR) call print_excitation_energies('ppRPA@RHF','2p (triplets)',nVVt,ee_trip_Om)
  if(print_ppLR) call print_excitation_energies('ppRPA@RHF','2h (triplets)',nOOt,hh_trip_Om)

  deallocate(Bpp,Cpp,Dpp)

  !-------------------
  !    Updating
  !-------------------
  old_eh_sing_Om(:) = eh_sing_Om(:)
  old_eh_trip_Om(:) = eh_trip_Om(:)
  old_ee_sing_Om(:) = ee_sing_Om(:)
  old_hh_sing_Om(:) = hh_sing_Om(:)
  old_ee_trip_Om(:) = ee_trip_Om(:)
  old_hh_trip_Om(:) = hh_trip_Om(:)
  
  deallocate(eh_sing_Om,eh_trip_Om,ee_sing_Om,hh_sing_Om,ee_trip_Om,hh_trip_Om)

  allocate(eh_sing_rho(nOrb,nOrb,nS))
  call R_eh_singlet_screened_integral(nOrb,nC,nO,nR,nS,ERI,sing_XpY,eh_sing_rho)
  deallocate(sing_XpY,sing_XmY)

  allocate(eh_trip_rho(nOrb,nOrb,nS))
  call R_eh_triplet_screened_integral(nOrb,nC,nO,nR,nS,ERI,trip_XpY,eh_trip_rho)
  deallocate(trip_XpY,trip_XmY)
  
  allocate(ee_sing_rho(nOrb,nOrb,nVVs),hh_sing_rho(nOrb,nOrb,nOOs))
  allocate(pp_sing_Gam(nOrb,nOrb,nOrb,nOrb))
  pp_sing_Gam(:,:,:,:) = 0d0
  call R_pp_singlet_screened_integral(nOrb,nC,nO,nV,nR,nOOs,nVVs,ERI,pp_sing_Gam,X1s,Y1s,ee_sing_rho,X2s,Y2s,hh_sing_rho)
  deallocate(X1s,Y1s,X2s,Y2s)
  deallocate(pp_sing_Gam)
  
  allocate(ee_trip_rho(nOrb,nOrb,nVVt),hh_trip_rho(nOrb,nOrb,nOOt))
  allocate(pp_trip_Gam(nOrb,nOrb,nOrb,nOrb))
  pp_trip_Gam(:,:,:,:) = 0d0
  call R_pp_triplet_screened_integral(nOrb,nC,nO,nV,nR,nOOt,nVVt,ERI,pp_trip_Gam,X1t,Y1t,ee_trip_rho,X2t,Y2t,hh_trip_rho)
  deallocate(X1t,Y1t,X2t,Y2t)
  deallocate(pp_trip_Gam)
  
!------------------------------------------------------------------------
! Loop on one-body energies
!------------------------------------------------------------------------

  do while(err_one_body > conv_one_body .and. n_it_macro <= max_it_macro)

     write(*,*)
     write(*,*)'************ Macro iteration number ',n_it_macro,' ************'
     write(*,*)'---------------------------------------------------------------'
     write(*,*)

     
     do while(err_two_body > conv_two_body .and. n_it_micro <= max_it_micro)

        !TODO add some timers everywhere
        write(*,*)
        write(*,*)'          Micro iteration number ',n_it_micro
        write(*,*)'          -----------------------------------'
        write(*,*)

        !-------------------
        ! Density channel
        !-------------------
        write(*,*)'Diagonalizing singlet ehBSE:'
        write(*,*)'----------------------------'
        allocate(Aph(nS,nS),Bph(nS,nS),eh_sing_Om(nS),sing_XpY(nS,nS),sing_XmY(nS,nS),eh_sing_Gam_A(nS,nS),eh_sing_Gam_B(nS,nS))
        Aph(:,:) = 0d0
        Bph(:,:) = 0d0
        call wall_time(start_t)
        if(.not.TDA) call phRLR_B(1,.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)
                     call phRLR_A(1,.false.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)

        call R_eh_singlet_Gamma_A(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                    old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, &
                    old_ee_sing_Om,ee_sing_rho,old_ee_trip_Om,ee_trip_rho, &
                    old_hh_sing_Om,hh_sing_rho,old_hh_trip_Om,hh_trip_rho, eh_sing_Gam_A)
        call R_eh_singlet_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                    old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, &
                    old_ee_sing_Om,ee_sing_rho,old_ee_trip_Om,ee_trip_rho, &
                    old_hh_sing_Om,hh_sing_rho,old_hh_trip_Om,hh_trip_rho, eh_sing_Gam_B)
        
        Aph(:,:) = Aph(:,:) + eh_sing_Gam_A(:,:)
        Bph(:,:) = Bph(:,:) + eh_sing_Gam_B(:,:)             
        call phRLR(TDA,nS,Aph,Bph,EcRPA,eh_sing_Om,sing_XpY,sing_XmY)
        call wall_time(end_t)
        t = end_t - start_t
        write(*,'(A52,1X,F9.3,A8)') 'Total wall time for singlet phBSE problem =',t,' seconds'
        write(*,*)
        if (print_phLR) call print_excitation_energies('phBSE@Parquet','singlet',nS,eh_sing_Om)
        err_eh_sing = maxval(abs(old_eh_sing_Om - eh_sing_Om))
        deallocate(Aph,Bph,eh_sing_Gam_A,eh_sing_Gam_B)

        !-------------------
        ! Magnetic channel
        !-------------------
        write(*,*)'Diagonalizing triplet ehBSE:'
        write(*,*)'----------------------------'
        allocate(Aph(nS,nS),Bph(nS,nS),eh_trip_Om(nS),trip_XpY(nS,nS),trip_XmY(nS,nS),eh_trip_Gam_A(nS,nS),eh_trip_Gam_B(nS,nS))
        Aph(:,:) = 0d0
        Bph(:,:) = 0d0
        call wall_time(start_t)
        if(.not.TDA) call phRLR_B(2,.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)
                     call phRLR_A(2,.false.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)

        call R_eh_triplet_Gamma_A(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                    old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, &
                    old_ee_sing_Om,ee_sing_rho,old_ee_trip_Om,ee_trip_rho, &
                    old_hh_sing_Om,hh_sing_rho,old_hh_trip_Om,hh_trip_rho, eh_trip_Gam_A)
        call R_eh_triplet_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                    old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, &
                    old_ee_sing_Om,ee_sing_rho,old_ee_trip_Om,ee_trip_rho, &
                    old_hh_sing_Om,hh_sing_rho,old_hh_trip_Om,hh_trip_rho, eh_trip_Gam_B)
        
        Aph(:,:) = Aph(:,:) + eh_trip_Gam_A(:,:)
        Bph(:,:) = Bph(:,:) + eh_trip_Gam_B(:,:)
        call phRLR(TDA,nS,Aph,Bph,EcRPA,eh_trip_Om,trip_XpY,trip_XmY)
        call wall_time(end_t)
        t = end_t - start_t
        write(*,'(A52,1X,F9.3,A8)') 'Total wall time for triplet phBSE problem =',t,' seconds'
        write(*,*)
        if (print_phLR) call print_excitation_energies('phBSE@Parquet','triplet',nS,eh_trip_Om)
        err_eh_trip = maxval(abs(old_eh_trip_Om - eh_trip_Om))
        deallocate(Aph,Bph,eh_trip_Gam_A,eh_trip_Gam_B)
        
        !-------------------
        ! Singlet channel
        !-------------------
        write(*,*)'Diagonalizing singlet ppBSE:'
        write(*,*)'----------------------------'
        allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs),   &
                 ee_sing_Om(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs), & 
                 hh_sing_Om(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs), &
                 pp_sing_Gam_B(nVVs,nOOs),pp_sing_Gam_C(nVVs,nVVs),pp_sing_Gam_D(nOOs,nOOs))
        Bpp(:,:) = 0d0
        Cpp(:,:) = 0d0
        Dpp(:,:) = 0d0
        call wall_time(start_t)
        if(.not.TDA) call ppRLR_B(1,nOrb,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                     call ppRLR_C(1,nOrb,nC,nO,nV,nR,nVVs,1d0,eHF,ERI,Cpp)
                     call ppRLR_D(1,nOrb,nC,nO,nV,nR,nOOs,1d0,eHF,ERI,Dpp)

        call R_pp_singlet_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,&
                                                           old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, pp_sing_Gam_B)
        call R_pp_singlet_Gamma_C(nOrb,nC,nO,nV,nR,nS,nVVs,old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, pp_sing_Gam_C)
        call R_pp_singlet_Gamma_D(nOrb,nC,nO,nV,nR,nS,nOOs,old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, pp_sing_Gam_D)
                     
        Bpp(:,:) = Bpp(:,:) + pp_sing_Gam_B(:,:)
        Cpp(:,:) = Cpp(:,:) + pp_sing_Gam_C(:,:)
        Dpp(:,:) = Dpp(:,:) + pp_sing_Gam_D(:,:)
        
        call ppRLR(TDA,nOOs,nVVs,Bpp,Cpp,Dpp,ee_sing_Om,X1s,Y1s,hh_sing_Om,X2s,Y2s,EcRPA)
        call wall_time(end_t)
        t = end_t - start_t
        write(*,'(A52,1X,F9.3,A8)') 'Total wall time for singlet ppBSE problem =',t,' seconds'
        write(*,*)
        if (print_ppLR) call print_excitation_energies('ppBSE@Parquet','2p (singlets)',nVVs,ee_sing_Om)
        if (print_ppLR) call print_excitation_energies('ppBSE@Parquet','2h (singlets)',nOOs,hh_sing_Om)
        err_ee_sing = maxval(abs(old_ee_sing_Om - ee_sing_Om))
        err_hh_sing = maxval(abs(old_hh_sing_Om - hh_sing_Om))
        deallocate(Bpp,Cpp,Dpp,pp_sing_Gam_B,pp_sing_Gam_C,pp_sing_Gam_D)
        
        !-------------------
        ! Triplet channel
        !-------------------
        write(*,*)'Diagonalizing triplet ppBSE:'
        write(*,*)'----------------------------'
        allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt),   &
                 ee_trip_Om(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt), & 
                 hh_trip_Om(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt), &
                 pp_trip_Gam_B(nVVt,nOOt),pp_trip_Gam_C(nVVt,nVVt),pp_trip_Gam_D(nOOt,nOOt))
        Bpp(:,:)=0d0
        Cpp(:,:)=0d0
        Dpp(:,:)=0d0
        call wall_time(start_t)
        if(.not.TDA) call ppRLR_B(2,nOrb,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                     call ppRLR_C(2,nOrb,nC,nO,nV,nR,nVVt,1d0,eHF,ERI,Cpp)
                     call ppRLR_D(2,nOrb,nC,nO,nV,nR,nOOt,1d0,eHF,ERI,Dpp)

        call R_pp_triplet_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOOt,nVVt,&
                                                           old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, pp_trip_Gam_B)
        call R_pp_triplet_Gamma_C(nOrb,nC,nO,nV,nR,nS,nVVt,old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, pp_trip_Gam_C)
        call R_pp_triplet_Gamma_D(nOrb,nC,nO,nV,nR,nS,nOOt,old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, pp_trip_Gam_D)
                     
        Bpp(:,:) = Bpp(:,:) + pp_trip_Gam_B(:,:)
        Cpp(:,:) = Cpp(:,:) + pp_trip_Gam_C(:,:)
        Dpp(:,:) = Dpp(:,:) + pp_trip_Gam_D(:,:)
        
        call ppRLR(TDA,nOOt,nVVt,Bpp,Cpp,Dpp,ee_trip_Om,X1t,Y1t,hh_trip_Om,X2t,Y2t,EcRPA)
        call wall_time(end_t)
        t = end_t - start_t
        write(*,'(A52,1X,F9.3,A8)') 'Total wall time for triplet ppBSE problem =',t,' seconds'
        write(*,*)
        if(print_ppLR) call print_excitation_energies('ppBSE@Parquet','2p (triplets)',nVVt,ee_trip_Om)
        if(print_ppLR) call print_excitation_energies('ppBSE@Parquet','2h (triplets)',nOOt,hh_trip_Om)
        err_ee_trip = maxval(abs(old_ee_trip_Om - ee_trip_Om))
        err_hh_trip = maxval(abs(old_hh_trip_Om - hh_trip_Om))
        deallocate(Bpp,Cpp,Dpp,pp_trip_Gam_B,pp_trip_Gam_C,pp_trip_Gam_D)
        
        write(*,*)
        write(*,*)'Error for density channel:', err_eh_sing
        write(*,*)'Error for magnetic channel:',err_eh_trip
        write(*,*)'Error for singlet channel:', max(err_ee_sing,err_hh_sing)
        write(*,*)'Error for triplet channel:', max(err_ee_trip,err_hh_trip)
        write(*,*)

        !-------------------
        !    Updating
        !-------------------
        old_eh_sing_Om(:) = eh_sing_Om(:)
        old_eh_trip_Om(:) = eh_trip_Om(:)
        old_ee_sing_Om(:) = ee_sing_Om(:)
        old_hh_sing_Om(:) = hh_sing_Om(:)
        old_ee_trip_Om(:) = ee_trip_Om(:)
        old_hh_trip_Om(:) = hh_trip_Om(:)

        allocate(pp_trip_Gam(nOrb,nOrb,nOrb,nOrb))
        call wall_time(start_t)
        call R_pp_triplet_Gamma(nOrb,nC,nR,nS,eh_sing_Om,eh_sing_rho,eh_trip_Om,eh_trip_rho,pp_trip_Gam)
        call wall_time(end_t)
        t = end_t - start_t
        write(*,'(A52,1X,F9.3,A8)') 'Total wall time for building pp triplet Gamma =',t,' seconds'
        write(*,*)
        allocate(pp_sing_Gam(nOrb,nOrb,nOrb,nOrb))
        call wall_time(start_t)
        call R_pp_singlet_Gamma(nOrb,nC,nR,nS,eh_sing_Om,eh_sing_rho,eh_trip_Om,eh_trip_rho,pp_sing_Gam)
        call wall_time(end_t)
        t = end_t - start_t
        write(*,'(A52,1X,F9.3,A8)') 'Total wall time for building pp singlet Gamma =',t,' seconds'
        write(*,*)
        
        deallocate(eh_sing_Om,eh_trip_Om,ee_sing_Om,hh_sing_Om,ee_trip_Om,hh_trip_Om)
        deallocate(eh_sing_rho,eh_trip_rho,ee_sing_rho,ee_trip_rho,hh_sing_rho,hh_trip_rho)
        
        allocate(eh_sing_rho(nOrb,nOrb,nS))
        call R_eh_singlet_screened_integral(nOrb,nC,nO,nR,nS,ERI,sing_XpY,eh_sing_rho)
        deallocate(sing_XpY,sing_XmY)
  
        allocate(eh_trip_rho(nOrb,nOrb,nS))
        call R_eh_triplet_screened_integral(nOrb,nC,nO,nR,nS,ERI,trip_XpY,eh_trip_rho)
        deallocate(trip_XpY,trip_XmY)
        
        allocate(ee_sing_rho(nOrb,nOrb,nVVs),hh_sing_rho(nOrb,nOrb,nOOs))
        call R_pp_singlet_screened_integral(nOrb,nC,nO,nV,nR,nOOs,nVVs,ERI,pp_sing_Gam,X1s,Y1s,ee_sing_rho,X2s,Y2s,hh_sing_rho)
        deallocate(X1s,Y1s,X2s,Y2s,pp_sing_Gam)

        allocate(ee_trip_rho(nOrb,nOrb,nVVt),hh_trip_rho(nOrb,nOrb,nOOt))
        call R_pp_triplet_screened_integral(nOrb,nC,nO,nV,nR,nOOt,nVVt,ERI,pp_trip_Gam,X1t,Y1t,ee_trip_rho,X2t,Y2t,hh_trip_rho)
        deallocate(X1t,Y1t,X2t,Y2t,pp_trip_Gam)
        
        ! Convergence criteria
        err_two_body = max(err_eh_sing,err_eh_trip,err_ee_sing,err_ee_trip,err_hh_sing,err_hh_trip)
        
        n_it_micro = n_it_micro + 1
     end do
     !------------------------------------------------------------------------
     ! End main loop
     !------------------------------------------------------------------------

     ! Did it actually converge?
     if(n_it_micro == max_it_micro+1) then
        write(*,*)
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'             Two-body convergence failed            '
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
        stop
     else
        write(*,*)
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'             Two-body convergence success           '
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
        call print_excitation_energies('phBSE@Parquet','singlet',nS,old_eh_sing_Om)
        call print_excitation_energies('phBSE@Parquet','triplet',nS,old_eh_trip_Om)
        call print_excitation_energies('ppBSE@Parquet','2p (singlets)',nVVs,old_ee_sing_Om)
        call print_excitation_energies('ppBSE@Parquet','2h (singlets)',nOOs,old_hh_sing_Om)
        call print_excitation_energies('ppBSE@Parquet','2p (triplets)',nVVt,old_ee_trip_Om)
        call print_excitation_energies('ppBSE@Parquet','2h (triplets)',nOOt,old_hh_trip_Om)
     end if

     allocate(eParquetlin(nOrb),Z(nOrb),SigC(nOrb)) 
     write(*,*) 'Building self-energy'
     call wall_time(start_t)
     call R_irred_Parquet_self_energy(nOrb,nC,nO,nV,nR,old_eParquet,EcGM,SigC,Z)
     call wall_time(end_t)
     t = end_t - start_t
     write(*,'(A52,1X,F9.3,A8)') 'Total wall time for building self energy =',t,' seconds'
     write(*,*)

     eParquetlin(:) = eHF(:) !+ Z(:)*SigC(:)
     ! Solve the quasi-particle equation
     if(linearize) then
        eParquet(:) = eParquetlin(:)
     else
        write(*,*)
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'  Newton-Raphson for Dyson equation not implemented  '
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
        stop
     end if
     deallocate(eParquetlin,Z,SigC)
     
     err_one_body =  maxval(abs(old_eParquet - eParquet))
     old_eParquet(:) = eParquet(:)
     
     n_it_macro = n_it_macro + 1
  end do ! End the macro loop

       ! Did it actually converge?
     if(n_it_macro == max_it_macro+1) then
        write(*,*)
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'             One-body convergence failed            '
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
        stop
     else
        write(*,*)
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)'             One-body convergence success           '
        write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(*,*)
     end if
     
end subroutine
