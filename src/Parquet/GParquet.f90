subroutine GParquet(max_it_1b,conv_1b,max_it_2b,conv_2b,nOrb,nC,nO,nV,nR,nS,eHF,ERI)

! Parquet approximation based on restricted orbitals

  implicit none
  include 'parameters.h'

! Hard-coded parameters

  logical                       :: linearize = .true.
  logical                       :: TDA = .true.
  logical                       :: print_phLR = .false.
  logical                       :: print_ppLR = .false.

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
  double precision              :: err_eh, err_hh, err_ee
  double precision              :: start_t, end_t, t
  double precision              :: start_1b, end_1b, t_1b
  double precision              :: start_2b, end_2b, t_2b

  integer                       :: nOO,nVV

  double precision              :: EcRPA
  double precision,allocatable  :: Aph(:,:), Bph(:,:)
  double precision,allocatable  :: XpY(:,:),XmY(:,:)
  double precision,allocatable  :: eh_Om(:), old_eh_Om(:)

  double precision,allocatable  :: Bpp(:,:), Cpp(:,:), Dpp(:,:)
  double precision,allocatable  :: X1(:,:),Y1(:,:)
  double precision,allocatable  :: ee_Om(:), old_ee_Om(:)
  double precision,allocatable  :: X2(:,:),Y2(:,:)
  double precision,allocatable  :: hh_Om(:), old_hh_Om(:)

  double precision,allocatable  :: eh_rho(:,:,:), ee_rho(:,:,:), hh_rho(:,:,:)

  double precision,allocatable  :: eh_Gam_A(:,:),eh_Gam_B(:,:)
  double precision,allocatable  :: pp_Gam_B(:,:),pp_Gam_C(:,:),pp_Gam_D(:,:)
  double precision,allocatable  :: eh_Gam(:,:,:,:),pp_Gam(:,:,:,:)

  double precision,allocatable  :: eQPlin(:),eQP(:),eOld(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision              :: EcGM
  
! Output variables
! None
  
  nOO = nO*(nO - 1)/2
  nVV = nV*(nV - 1)/2

  allocate(eQP(nOrb),eOld(nOrb))
    
  write(*,*)
  write(*,*)'***********************************'
  write(*,*)'* Generalized Parquet Calculation *'
  write(*,*)'***********************************'
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

  allocate(old_eh_Om(nS),old_ee_Om(nVV),old_hh_Om(nOO))
  allocate(eh_rho(nOrb,nOrb,nS),ee_rho(nOrb,nOrb,nVV),hh_rho(nOrb,nOrb,nOO))

! Initialization

  n_it_1b = 0
  err_1b  = 1d0

  n_it_2b = 0 
  err_2b  = 1d0

  eQP(:)  = eHF(:)
  eOld(:) = eHF(:)

  eh_rho(:,:,:) = 0d0
  ee_rho(:,:,:) = 0d0
  hh_rho(:,:,:) = 0d0

  old_eh_Om(:) = 1d0
  old_ee_Om(:) = 1d0
  old_hh_Om(:) = 1d0

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
      allocate(eh_Gam(nOrb,nOrb,nOrb,nOrb))
      allocate(pp_Gam(nOrb,nOrb,nOrb,nOrb))

      ! Build eh effective interaction
      write(*,*) 'Computing eh effective interaction...'

      call wall_time(start_t)
      !call R_eh_Gamma(nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
      !     old_eh_Om,eh_rho,old_ee_Om,ee_rho,old_hh_Om,hh_rho, &
      !     eh_trip_Gam)
      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for eh Gamma =',t,' seconds'
      write(*,*)

     ! Build singlet pp effective interaction
      write(*,*) 'Computing pp effective interaction...'

      call wall_time(start_t)
      !call R_pp_Gamma(nOrb,nC,nR,nS,old_eh_Om,eh_rho,pp_Gam)
      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for pp Gamma =',t,' seconds'
      write(*,*)


      
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

      call print_excitation_energies('phBSE@Parquet','1h1p',nS,old_eh_Om)
      call print_excitation_energies('ppBSE@Parquet','2p',nVV,old_ee_Om)
      call print_excitation_energies('ppBSE@Parquet','2h',nOO,old_hh_Om)

    end if

    allocate(eQPlin(nOrb),Z(nOrb),SigC(nOrb)) 

    write(*,*) 'Building self-energy'

    call wall_time(start_t)
    !call G_irred_Parquet_self_energy(nOrb,nC,nO,nV,nR,eOld,EcGM,SigC,Z)
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
  
end subroutine GParquet
