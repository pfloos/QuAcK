subroutine RParquet(max_it_1b,conv_1b,max_it_2b,conv_2b,nOrb,nC,nO,nV,nR,nS,eHF,ERI)

! Parquet approximation based on restricted orbitals

  implicit none
  include 'parameters.h'

! Hard-coded parameters

  logical                       :: linearize = .true.
  logical                       :: TDA = .true.
  logical                       :: print_phLR = .true.
  logical                       :: print_ppLR = .true.
  
! Input variables

  integer,intent(in)            :: max_it_1b,max_it_2b
  double precision,intent(in)   :: conv_1b,conv_2b
  integer,intent(in)            :: nOrb,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  
! Local variables

  integer                       :: ispin

  integer                       :: n_it_1b,n_it_2b
  double precision              :: err_1b,err_2b
  double precision              :: err_eh_sing,err_eh_trip
  double precision              :: err_hh_sing,err_hh_trip
  double precision              :: err_ee_sing,err_ee_trip
  double precision              :: start_t, end_t, t
  double precision              :: start_1b, end_1b, t_1b
  double precision              :: start_2b, end_2b, t_2b

  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  
  double precision              :: EcRPA
  double precision,allocatable  :: Aph(:,:), Bph(:,:)
  double precision,allocatable  :: sing_XpY(:,:),trip_XpY(:,:)
  double precision,allocatable  :: sing_XmY(:,:),trip_XmY(:,:)
  double precision,allocatable  :: eh_sing_Om(:), old_eh_sing_Om(:)
  double precision,allocatable  :: eh_trip_Om(:), old_eh_trip_Om(:)

  double precision,allocatable  :: Bpp(:,:), Cpp(:,:), Dpp(:,:)
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

  double precision,allocatable  :: eQPlin(:),eQP(:),eOld(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision              :: EcGM
  
! Output variables
! None
  
  nOOs = nO*(nO + 1)/2
  nVVs = nV*(nV + 1)/2
  nOOt = nO*(nO - 1)/2
  nVVt = nV*(nV - 1)/2

  allocate(eQP(nOrb),eOld(nOrb))
  
  write(*,*)
  write(*,*)'**********************************'
  write(*,*)'* Restricted Parquet Calculation *'
  write(*,*)'**********************************'
  write(*,*)

! Print parameters

  write(*,*)'---------------------------------------------------------------'
  write(*,*)' Parquet parameters for one-body and two-body self-consistency '
  write(*,*)'---------------------------------------------------------------'
  write(*,'(1X,A50,1X,I5)')    'Maximum number for one-body self-consistency:', max_it_1b
  write(*,'(1X,A50,1X,E10.5)') 'Convergence threshold for one-body energies:', conv_1b
  write(*,*)'---------------------------------------------------------------'
  write(*,'(1X,A50,1X,I5)')    'Maximum number for two-body self-consistency:', max_it_2b
  write(*,'(1X,A50,1X,E10.5)') 'Convergence threshold for two-body energies:', conv_2b
  write(*,*)'---------------------------------------------------------------'
  write(*,*)
  
  if(linearize) then 
      write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
      write(*,*)
   else 
      write(*,*) ' *** Quasiparticle energies obtained by root search *** '
      write(*,*)
  endif
  
! Memory allocation 

  allocate(old_eh_sing_Om(nS),old_eh_trip_Om(nS))
  allocate(old_ee_sing_Om(nVVs),old_hh_sing_Om(nOOs))
  allocate(old_ee_trip_Om(nVVt),old_hh_trip_Om(nOOt))
  allocate(eh_sing_rho(nOrb,nOrb,nS),eh_trip_rho(nOrb,nOrb,nS))
  allocate(ee_sing_rho(nOrb,nOrb,nVVs),hh_sing_rho(nOrb,nOrb,nOOs))
  allocate(ee_trip_rho(nOrb,nOrb,nVVt),hh_trip_rho(nOrb,nOrb,nOOt))

! Initialization

  n_it_1b = 0
  err_1b  = 1d0

  n_it_2b = 0 
  err_2b  = 1d0

  eQP(:)  = eHF(:)
  eOld(:) = eHF(:)

  eh_sing_rho(:,:,:) = 0d0
  eh_trip_rho(:,:,:) = 0d0
  ee_sing_rho(:,:,:) = 0d0
  ee_trip_rho(:,:,:) = 0d0
  hh_sing_rho(:,:,:) = 0d0
  hh_trip_rho(:,:,:) = 0d0

  old_eh_sing_Om(:) = 1d0
  old_eh_trip_Om(:) = 1d0
  old_ee_sing_Om(:) = 1d0
  old_ee_trip_Om(:) = 1d0
  old_hh_sing_Om(:) = 1d0
  old_hh_trip_Om(:) = 1d0

  !-----------------------------------------!
  ! Main loop for one-body self-consistency !
  !-----------------------------------------!

  do while(err_1b > conv_1b .and. n_it_1b < max_it_1b)

    n_it_1b = n_it_1b + 1
    call wall_time(start_1b)

    write(*,*)
    write(*,*)'====================================='
    write(*,'(1X,A30,1X,I4)') 'One-body iteration #',n_it_1b
    write(*,*)'====================================='
    write(*,*)

    !-----------------------------------------!
    ! Main loop for two-body self-consistency !
    !-----------------------------------------!
    
    do while(err_2b > conv_2b .and. n_it_2b < max_it_2b)

      n_it_2b = n_it_2b + 1
      call wall_time(start_2b)

      write(*,*)' ***********************************'
      write(*,'(1X,A30,1X,I4)') 'Two-body iteration #',n_it_2b
      write(*,*)' ***********************************'
      write(*,*)

      !--------------------------------!
      ! Compute effective interactions !
      !--------------------------------!

      ! Memory allocation
      allocate(eh_sing_Gam(nOrb,nOrb,nOrb,nOrb))
      allocate(eh_trip_Gam(nOrb,nOrb,nOrb,nOrb))
      allocate(pp_sing_Gam(nOrb,nOrb,nOrb,nOrb))
      allocate(pp_trip_Gam(nOrb,nOrb,nOrb,nOrb))
      
      ! Build singlet eh effective interaction
      write(*,*) 'Computing singlet eh effective interaction...'

      call wall_time(start_t)
      call R_eh_singlet_Gamma(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                  old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, &
                  old_ee_sing_Om,ee_sing_rho,old_ee_trip_Om,ee_trip_rho, &
                  old_hh_sing_Om,hh_sing_rho,old_hh_trip_Om,hh_trip_rho, eh_sing_Gam)
      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for eh singlet Gamma =',t,' seconds'
      write(*,*)

     ! Build triplet eh effective interaction
      write(*,*) 'Computing triplet eh effective interaction...'

      call wall_time(start_t)
      call R_eh_triplet_Gamma(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                  old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, &
                  old_ee_sing_Om,ee_sing_rho,old_ee_trip_Om,ee_trip_rho, &
                  old_hh_sing_Om,hh_sing_rho,old_hh_trip_Om,hh_trip_rho, eh_trip_Gam)
      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for eh triplet Gamma =',t,' seconds'
      write(*,*)

     ! Build singlet pp effective interaction
      write(*,*) 'Computing singlet pp effective interaction...'

      call wall_time(start_t)
      call R_pp_singlet_Gamma(nOrb,nC,nR,nS,old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho,pp_sing_Gam)
      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for pp singlet Gamma =',t,' seconds'
      write(*,*)

     ! Build triplet pp effective interaction
      write(*,*) 'Computing triplet pp effective interaction...'

      call wall_time(start_t)
      call R_pp_triplet_Gamma(nOrb,nC,nR,nS,old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho,pp_trip_Gam)
      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for pp triplet Gamma =',t,' seconds'
      write(*,*)

      !-----------------!
      ! Density channel !
      !-----------------!

      write(*,*)'   -------------------------------'
      write(*,*)'   | Diagonalizing singlet ehBSE |'
      write(*,*)'   -------------------------------'
      write(*,*)

      allocate(Aph(nS,nS),Bph(nS,nS),eh_sing_Om(nS),sing_XpY(nS,nS),sing_XmY(nS,nS),eh_sing_Gam_A(nS,nS),eh_sing_Gam_B(nS,nS))

      ispin = 1
      Aph(:,:) = 0d0
      Bph(:,:) = 0d0

      call wall_time(start_t)

                   call phRLR_A(ispin,.false.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
      if(.not.TDA) call phRLR_B(ispin,.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

      if(n_it_2b == 1) then

        eh_sing_Gam_A(:,:) = 0d0
        eh_sing_Gam_B(:,:) = 0d0

      else

        call R_eh_singlet_Gamma_A(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                    old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, &
                    old_ee_sing_Om,ee_sing_rho,old_ee_trip_Om,ee_trip_rho, &
                    old_hh_sing_Om,hh_sing_rho,old_hh_trip_Om,hh_trip_rho, & 
                    eh_sing_Gam_A)
       
        if(.not.TDA) call R_eh_singlet_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                                 old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, &
                                 old_ee_sing_Om,ee_sing_rho,old_ee_trip_Om,ee_trip_rho, &
                                 old_hh_sing_Om,hh_sing_rho,old_hh_trip_Om,hh_trip_rho, & 
                                 eh_sing_Gam_B)

      end if      
      Aph(:,:) = Aph(:,:) + eh_sing_Gam_A(:,:)
      Bph(:,:) = Bph(:,:) + eh_sing_Gam_B(:,:)             

      call phRLR(TDA,nS,Aph,Bph,EcRPA,eh_sing_Om,sing_XpY,sing_XmY)

      call wall_time(end_t)

      t = end_t - start_t
      write(*,'(A50,1X,F9.3,A8)') 'Wall time for singlet phBSE =',t,' seconds'
      write(*,*)

      if(print_phLR) call print_excitation_energies('phBSE@Parquet','singlet',nS,eh_sing_Om)

      err_eh_sing = maxval(abs(old_eh_sing_Om - eh_sing_Om))

      deallocate(Aph,Bph,eh_sing_Gam_A,eh_sing_Gam_B)

      !------------------!
      ! Magnetic channel !
      !------------------!

      write(*,*)'   -------------------------------'
      write(*,*)'   | Diagonalizing triplet ehBSE |'
      write(*,*)'   -------------------------------'
      write(*,*)

      allocate(Aph(nS,nS),Bph(nS,nS),eh_trip_Om(nS),trip_XpY(nS,nS),trip_XmY(nS,nS),eh_trip_Gam_A(nS,nS),eh_trip_Gam_B(nS,nS))

      ispin = 2
      Aph(:,:) = 0d0
      Bph(:,:) = 0d0

      call wall_time(start_t)

                   call phRLR_A(ispin,.false.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
      if(.not.TDA) call phRLR_B(ispin,.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

      if(n_it_2b == 1) then

        eh_trip_Gam_A(:,:) = 0d0
        eh_trip_Gam_B(:,:) = 0d0

      else

        call R_eh_triplet_Gamma_A(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                    old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, &
                    old_ee_sing_Om,ee_sing_rho,old_ee_trip_Om,ee_trip_rho, &
                    old_hh_sing_Om,hh_sing_rho,old_hh_trip_Om,hh_trip_rho, & 
                    eh_trip_Gam_A)
       
        if(.not.TDA) call R_eh_triplet_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                                 old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, &
                                 old_ee_sing_Om,ee_sing_rho,old_ee_trip_Om,ee_trip_rho, &
                                 old_hh_sing_Om,hh_sing_rho,old_hh_trip_Om,hh_trip_rho, & 
                                 eh_trip_Gam_B)

      end if
      
      Aph(:,:) = Aph(:,:) + eh_trip_Gam_A(:,:)
      Bph(:,:) = Bph(:,:) + eh_trip_Gam_B(:,:)

      call phRLR(TDA,nS,Aph,Bph,EcRPA,eh_trip_Om,trip_XpY,trip_XmY)

      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for triplet phBSE =',t,' seconds'
      write(*,*)

      if(print_phLR) call print_excitation_energies('phBSE@Parquet','triplet',nS,eh_trip_Om)

      err_eh_trip = maxval(abs(old_eh_trip_Om - eh_trip_Om))

      deallocate(Aph,Bph,eh_trip_Gam_A,eh_trip_Gam_B)
      
      !-----------------!
      ! Singlet channel !
      !-----------------!

      write(*,*)'   -------------------------------'
      write(*,*)'   | Diagonalizing singlet ppBSE |'
      write(*,*)'   -------------------------------'
      write(*,*)

      allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs),   &
               ee_sing_Om(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs), & 
               hh_sing_Om(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs), &
               pp_sing_Gam_B(nVVs,nOOs),pp_sing_Gam_C(nVVs,nVVs),pp_sing_Gam_D(nOOs,nOOs))

      ispin = 1
      Bpp(:,:) = 0d0
      Cpp(:,:) = 0d0
      Dpp(:,:) = 0d0

      call wall_time(start_t)
      if(.not.TDA) call ppRLR_B(ispin,nOrb,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                   call ppRLR_C(ispin,nOrb,nC,nO,nV,nR,nVVs,1d0,eHF,ERI,Cpp)
                   call ppRLR_D(ispin,nOrb,nC,nO,nV,nR,nOOs,1d0,eHF,ERI,Dpp)

      if(n_it_2b == 1) then

        pp_sing_Gam_B(:,:) = 0d0
        pp_sing_Gam_C(:,:) = 0d0
        pp_sing_Gam_D(:,:) = 0d0

      else

        if(.not.TDA) call R_pp_singlet_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,&
                                                           old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho,pp_sing_Gam_B)
        call R_pp_singlet_Gamma_C(nOrb,nC,nO,nV,nR,nS,nVVs,old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho,pp_sing_Gam_C)
        call R_pp_singlet_Gamma_D(nOrb,nC,nO,nV,nR,nS,nOOs,old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho,pp_sing_Gam_D)

      end if
                   
      Bpp(:,:) = Bpp(:,:) + pp_sing_Gam_B(:,:)
      Cpp(:,:) = Cpp(:,:) + pp_sing_Gam_C(:,:)
      Dpp(:,:) = Dpp(:,:) + pp_sing_Gam_D(:,:)
      
      call ppRLR(TDA,nOOs,nVVs,Bpp,Cpp,Dpp,ee_sing_Om,X1s,Y1s,hh_sing_Om,X2s,Y2s,EcRPA)
      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for singlet ppBSE =',t,' seconds'
      write(*,*)

      if(print_ppLR) call print_excitation_energies('ppBSE@Parquet','2p (singlets)',nVVs,ee_sing_Om)
      if(print_ppLR) call print_excitation_energies('ppBSE@Parquet','2h (singlets)',nOOs,hh_sing_Om)

      err_ee_sing = maxval(abs(old_ee_sing_Om - ee_sing_Om))
      err_hh_sing = maxval(abs(old_hh_sing_Om - hh_sing_Om))

      deallocate(Bpp,Cpp,Dpp,pp_sing_Gam_B,pp_sing_Gam_C,pp_sing_Gam_D)
      
      !-----------------!
      ! Triplet channel !
      !-----------------!

      write(*,*)'   -------------------------------'
      write(*,*)'   | Diagonalizing triplet ppBSE |'
      write(*,*)'   -------------------------------'
      write(*,*)

      allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt),   &
               ee_trip_Om(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt), & 
               hh_trip_Om(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt), &
               pp_trip_Gam_B(nVVt,nOOt),pp_trip_Gam_C(nVVt,nVVt),pp_trip_Gam_D(nOOt,nOOt))

      ispin = 2
      Bpp(:,:) = 0d0
      Cpp(:,:) = 0d0
      Dpp(:,:) = 0d0

      call wall_time(start_t)
      if(.not.TDA) call ppRLR_B(ispin,nOrb,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                   call ppRLR_C(ispin,nOrb,nC,nO,nV,nR,nVVt,1d0,eHF,ERI,Cpp)
                   call ppRLR_D(ispin,nOrb,nC,nO,nV,nR,nOOt,1d0,eHF,ERI,Dpp)

      if(n_it_2b == 1) then

        pp_trip_Gam_B(:,:) = 0d0
        pp_trip_Gam_C(:,:) = 0d0
        pp_trip_Gam_D(:,:) = 0d0

      else

        if(.not.TDA) call R_pp_triplet_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOOt,nVVt,&
                                                           old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, pp_trip_Gam_B)
        call R_pp_triplet_Gamma_C(nOrb,nC,nO,nV,nR,nS,nVVt,old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, pp_trip_Gam_C)
        call R_pp_triplet_Gamma_D(nOrb,nC,nO,nV,nR,nS,nOOt,old_eh_sing_Om,eh_sing_rho,old_eh_trip_Om,eh_trip_rho, pp_trip_Gam_D)

      end if
                   
      Bpp(:,:) = Bpp(:,:) + pp_trip_Gam_B(:,:)
      Cpp(:,:) = Cpp(:,:) + pp_trip_Gam_C(:,:)
      Dpp(:,:) = Dpp(:,:) + pp_trip_Gam_D(:,:)
      
      call ppRLR(TDA,nOOt,nVVt,Bpp,Cpp,Dpp,ee_trip_Om,X1t,Y1t,hh_trip_Om,X2t,Y2t,EcRPA)

      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for triplet ppBSE =',t,' seconds'
      write(*,*)

      if(print_ppLR) call print_excitation_energies('ppBSE@Parquet','2p (triplets)',nVVt,ee_trip_Om)
      if(print_ppLR) call print_excitation_energies('ppBSE@Parquet','2h (triplets)',nOOt,hh_trip_Om)

      err_ee_trip = maxval(abs(old_ee_trip_Om - ee_trip_Om))
      err_hh_trip = maxval(abs(old_hh_trip_Om - hh_trip_Om))

      deallocate(Bpp,Cpp,Dpp,pp_trip_Gam_B,pp_trip_Gam_C,pp_trip_Gam_D)
      
      write(*,*) '----------------------------------------'
      write(*,*) ' Two-body convergence '
      write(*,*) '----------------------------------------'
      write(*,'(1X,A30,F10.6)')'Error for density  channel = ',err_eh_sing
      write(*,'(1X,A30,F10.6)')'Error for magnetic channel = ',err_eh_trip
      write(*,'(1X,A30,F10.6)')'Error for singlet  channel = ',max(err_ee_sing,err_hh_sing)
      write(*,'(1X,A30,F10.6)')'Error for triplet  channel = ',max(err_ee_trip,err_hh_trip)
      write(*,*) '----------------------------------------'
      write(*,*)

      !----------!
      ! Updating !
      !----------!

      old_eh_sing_Om(:) = eh_sing_Om(:)
      old_eh_trip_Om(:) = eh_trip_Om(:)
      old_ee_sing_Om(:) = ee_sing_Om(:)
      old_hh_sing_Om(:) = hh_sing_Om(:)
      old_ee_trip_Om(:) = ee_trip_Om(:)
      old_hh_trip_Om(:) = hh_trip_Om(:)

      deallocate(eh_sing_Om,eh_trip_Om,ee_sing_Om,hh_sing_Om,ee_trip_Om,hh_trip_Om)
      
      !----------------------------!
      ! Compute screened integrals !
      !----------------------------!

      ! Free memory
      deallocate(eh_sing_rho,eh_trip_rho,ee_sing_rho,ee_trip_rho,hh_sing_rho,hh_trip_rho)
      ! TODO Once we will compute the blocks of kernel starting from the 4-tensors we can move the freeing up
      ! Memory allocation
      allocate(eh_sing_rho(nOrb,nOrb,nS))
      allocate(eh_trip_rho(nOrb,nOrb,nS))
      allocate(ee_sing_rho(nOrb,nOrb,nVVs),hh_sing_rho(nOrb,nOrb,nOOs))
      allocate(ee_trip_rho(nOrb,nOrb,nVVt),hh_trip_rho(nOrb,nOrb,nOOt))
      
      
      ! Build singlet eh screened integrals
      write(*,*) 'Computing singlet eh screened integrals...'

      call wall_time(start_t)
      call R_eh_singlet_screened_integral(nOrb,nC,nO,nR,nS,ERI,eh_sing_Gam,sing_XpY,sing_XmY,eh_sing_rho)
      call wall_time(end_t)
      t = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for singlet eh integrals =',t,' seconds'
      write(*,*)
      ! Done with eigenvectors and kernel
      deallocate(sing_XpY,sing_XmY,eh_sing_Gam)
  
      ! Build triplet eh screened integrals
      write(*,*) 'Computing triplet eh screened integrals...'

      call wall_time(start_t)
      call R_eh_triplet_screened_integral(nOrb,nC,nO,nR,nS,ERI,eh_trip_Gam,trip_XpY,trip_XmY,eh_trip_rho)
      call wall_time(end_t)
      t = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet eh integrals =',t,' seconds'
      write(*,*)
      ! Done with eigenvectors and kernel
      deallocate(trip_XpY,trip_XmY,eh_trip_Gam)
      
      ! Build singlet pp screened integrals
      write(*,*) 'Computing singlet pp screened integrals...'

      call wall_time(start_t)
      call R_pp_singlet_screened_integral(nOrb,nC,nO,nV,nR,nOOs,nVVs,ERI,pp_sing_Gam,X1s,Y1s,ee_sing_rho,X2s,Y2s,hh_sing_rho)
      call wall_time(end_t)
      t = end_t - start_t
      ! Done with eigenvectors and kernel
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for singlet pp integrals =',t,' seconds'
      write(*,*)

      deallocate(X1s,Y1s,X2s,Y2s,pp_sing_Gam)

      ! Build triplet pp screened integrals
      write(*,*) 'Computing triplet pp screened integrals...'

      call wall_time(start_t)
      call R_pp_triplet_screened_integral(nOrb,nC,nO,nV,nR,nOOt,nVVt,ERI,pp_trip_Gam,X1t,Y1t,ee_trip_rho,X2t,Y2t,hh_trip_rho)
      call wall_time(end_t)
      t = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for triplet pp integrals =',t,' seconds'
      write(*,*)
      ! Done with eigenvectors and kernel
      deallocate(X1t,Y1t,X2t,Y2t,pp_trip_Gam)

      ! Convergence criteria
      err_2b = max(err_eh_sing,err_eh_trip,err_ee_sing,err_ee_trip,err_hh_sing,err_hh_trip)
       
      call wall_time(end_2b)
      t_2b = end_2b - start_2b
      write(*,'(A50,1X,F9.3,A8)') 'Wall time for two-body iteration =',t_2b,' seconds'
      write(*,*)

    end do
    !---------------------------------------------!
    ! End main loop for two-body self-consistency !
    !---------------------------------------------!

    ! Did it actually converge?

    if(n_it_2b == max_it_2b) then

      write(*,*)
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)'             Two-body convergence failed            '
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      stop

    else

      write(*,*)
      write(*,*)'****************************************************'
      write(*,*)'             Two-body convergence success           '
      write(*,*)'****************************************************'
      write(*,*)

      call print_excitation_energies('phBSE@Parquet','singlet',nS,old_eh_sing_Om)
      call print_excitation_energies('phBSE@Parquet','triplet',nS,old_eh_trip_Om)
      call print_excitation_energies('ppBSE@Parquet','2p (singlets)',nVVs,old_ee_sing_Om)
      call print_excitation_energies('ppBSE@Parquet','2h (singlets)',nOOs,old_hh_sing_Om)
      call print_excitation_energies('ppBSE@Parquet','2p (triplets)',nVVt,old_ee_trip_Om)
      call print_excitation_energies('ppBSE@Parquet','2h (triplets)',nOOt,old_hh_trip_Om)

    end if

    allocate(eQPlin(nOrb),Z(nOrb),SigC(nOrb)) 

    write(*,*) 'Building self-energy'

    call wall_time(start_t)
    call R_irred_Parquet_self_energy(nOrb,nC,nO,nV,nR,eOld,EcGM,SigC,Z)
    call wall_time(end_t)
    t = end_t - start_t
    write(*,'(A50,1X,F9.3,A8)') 'Wall time for self energy =',t,' seconds'
    write(*,*)

    eQPlin(:) = eHF(:) !+ Z(:)*SigC(:)

    ! Solve the quasi-particle equation

    if(linearize) then

       eQP(:) = eQPlin(:)

    else

      write(*,*)
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)'  Newton-Raphson for Dyson equation not implemented  '
      write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
      write(*,*)
      stop

    end if

    deallocate(eQPlin,Z,SigC)

    ! Check one-body converge

    err_1b =  maxval(abs(eOld - eQP))
    eOld(:) = eQP(:)
 
    call wall_time(end_1b)
    t_1b = end_1b - start_1b
    write(*,'(A50,1X,F9.3,A8)') 'Wall time for one-body iteration =',t_1b,' seconds'
    
  end do 
  !---------------------------------------------!
  ! End main loop for one-body self-consistency !
  !---------------------------------------------!

  ! Did it actually converge?
  if(n_it_1b == max_it_1b) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'             One-body convergence failed            '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)
    stop

  else

    write(*,*)
    write(*,*)'****************************************************'
    write(*,*)'             One-body convergence success           '
    write(*,*)'****************************************************'
    write(*,*)

  end if
     
end subroutine RParquet
