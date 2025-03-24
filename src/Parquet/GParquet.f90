subroutine GParquet(max_it_1b,conv_1b,max_it_2b,conv_2b,nOrb,nC,nO,nV,nR,nS,eHF,ERI)

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
      !call R_eh_Gamma(nOrb,nC,nO,nV,nR,nS,nOO,nVV, &
      !     old_eh_Om,eh_rho,old_ee_Om,ee_rho,old_hh_Om,hh_rho, &
      !     eh_Gam)
      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for eh Gamma =',t,' seconds'
      write(*,*)

     ! Build singlet pp effective interaction
      write(*,*) 'Computing pp effective interaction...'

      call wall_time(start_t)
      call G_pp_Gamma(nOrb,nC,nO,nV,nR,nS,old_eh_Om,eh_rho,pp_Gam)
      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for pp Gamma =',t,' seconds'
      write(*,*)

      !-----------------!
      !    eh channel   !
      !-----------------!

      write(*,*)'   ------------------------------'
      write(*,*)'   |     Diagonalizing ehBSE    |'
      write(*,*)'   ------------------------------'
      write(*,*)

      allocate(Aph(nS,nS),Bph(nS,nS),eh_Om(nS),XpY(nS,nS),XmY(nS,nS),eh_Gam_A(nS,nS),eh_Gam_B(nS,nS))


      Aph(:,:) = 0d0
      Bph(:,:) = 0d0

      call wall_time(start_t)

                   call phGLR_A(.false.,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
      if(.not.TDA) call phGLR_B(.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)
      
      if(n_it_2b == 1) then

        eh_Gam_A(:,:) = 0d0
        eh_Gam_B(:,:) = 0d0

      else

        call G_eh_Gamma_A(nOrb,nC,nO,nV,nR,nS,nOO,nVV,           &
             old_eh_Om,eh_rho,old_ee_Om,ee_rho,old_hh_Om,hh_rho, & 
             eh_Gam_A)
       
        if(.not.TDA) call G_eh_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOO,nVV,       &
                      old_eh_Om,eh_rho,old_ee_Om,ee_rho,old_hh_Om,hh_rho, & 
                      eh_Gam_B)

      end if
      
      Aph(:,:) = Aph(:,:) + eh_Gam_A(:,:)
      Bph(:,:) = Bph(:,:) + eh_Gam_B(:,:)             

      call phGLR(TDA,nS,Aph,Bph,EcRPA,eh_Om,XpY,XmY)

      call wall_time(end_t)

      t = end_t - start_t
      write(*,'(A50,1X,F9.3,A8)') 'Wall time for phBSE =',t,' seconds'
      write(*,*)

      if(print_phLR) call print_excitation_energies('phBSE@Parquet','eh generalized',nS,eh_Om)

      err_eh = maxval(abs(old_eh_Om - eh_Om))

      deallocate(Aph,Bph,eh_Gam_A,eh_Gam_B)

      !-----------------!
      !    pp channel   !
      !-----------------!

      write(*,*)'   ------------------------------'
      write(*,*)'   |     Diagonalizing ppBSE    |'
      write(*,*)'   ------------------------------'
      write(*,*)

      allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO),   &
               ee_Om(nVV),X1(nVV,nVV),Y1(nOO,nVV), & 
               hh_Om(nOO),X2(nVV,nOO),Y2(nOO,nOO), &
               pp_Gam_B(nVV,nOO),pp_Gam_C(nVV,nVV),pp_Gam_D(nOO,nOO))

      Bpp(:,:) = 0d0
      Cpp(:,:) = 0d0
      Dpp(:,:) = 0d0

      call wall_time(start_t)
      if(.not.TDA) call ppGLR_B(nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)
                   call ppGLR_C(nOrb,nC,nO,nV,nR,nVV,1d0,eHF,ERI,Cpp)
                   call ppGLR_D(nOrb,nC,nO,nV,nR,nOO,1d0,eHF,ERI,Dpp)

      if(n_it_2b == 1) then

        pp_Gam_B(:,:) = 0d0
        pp_Gam_C(:,:) = 0d0
        pp_Gam_D(:,:) = 0d0

      else

        if(.not.TDA) call G_pp_Gamma_B(nOrb,nC,nO,nV,nR,nS,nOO,nVV,old_eh_Om,eh_rho,pp_Gam_B)
        call G_pp_Gamma_C(nOrb,nC,nO,nV,nR,nS,nVV,old_eh_Om,eh_rho,pp_Gam_C)
        call G_pp_Gamma_D(nOrb,nC,nO,nV,nR,nS,nOO,old_eh_Om,eh_rho,pp_Gam_D)

      end if
                   
      Bpp(:,:) = Bpp(:,:) + pp_Gam_B(:,:)
      Cpp(:,:) = Cpp(:,:) + pp_Gam_C(:,:)
      Dpp(:,:) = Dpp(:,:) + pp_Gam_D(:,:)
      
      call ppGLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,ee_Om,X1,Y1,hh_Om,X2,Y2,EcRPA)
      call wall_time(end_t)
      t = end_t - start_t

      write(*,'(A50,1X,F9.3,A8)') 'Wall time for ppBSE =',t,' seconds'
      write(*,*)

      if(print_ppLR) call print_excitation_energies('ppBSE@Parquet','2p generalized',nVV,ee_Om)
      if(print_ppLR) call print_excitation_energies('ppBSE@Parquet','2h generalized',nOO,hh_Om)

      err_ee = maxval(abs(old_ee_Om - ee_Om))
      err_hh = maxval(abs(old_hh_Om - hh_Om))

      deallocate(Bpp,Cpp,Dpp,pp_Gam_B,pp_Gam_C,pp_Gam_D)

            
      write(*,*) '----------------------------------------'
      write(*,*) ' Two-body convergence '
      write(*,*) '----------------------------------------'
      write(*,'(1X,A30,F10.6)')'Error for eh channel = ',err_eh
      write(*,'(1X,A30,F10.6)')'Error for pp channel = ',max(err_ee,err_hh)
      write(*,*) '----------------------------------------'
      write(*,*)

      !----------!
      ! Updating !
      !----------!

      old_eh_Om(:) = eh_Om(:)
      old_ee_Om(:) = ee_Om(:)
      old_hh_Om(:) = hh_Om(:)

      deallocate(eh_Om,ee_Om,hh_Om)
      
      !----------------------------!
      ! Compute screened integrals !
      !----------------------------!

      ! Free memory
      deallocate(eh_rho,ee_rho,hh_rho)
      ! TODO Once we will compute the blocks of kernel starting from the 4-tensors we can move the freeing up
      ! Memory allocation
      allocate(eh_rho(nOrb,nOrb,nS))
      allocate(ee_rho(nOrb,nOrb,nVV),hh_rho(nOrb,nOrb,nOO))

      ! Build singlet eh integrals
      write(*,*) 'Computing eh screened integrals...'

      call wall_time(start_t)
      call G_eh_screened_integral(nOrb,nC,nO,nR,nS,ERI,eh_Gam,XpY,XmY,eh_rho)
      call wall_time(end_t)
      t = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for eh integrals =',t,' seconds'
      write(*,*)
      ! Done with eigenvectors and kernel
      deallocate(XpY,XmY,eh_Gam)

      ! Build singlet pp integrals
      write(*,*) 'Computing pp screened integrals...'

      call wall_time(start_t)
      call G_pp_screened_integral(nOrb,nC,nO,nV,nR,nOO,nVV,ERI,pp_Gam,X1,Y1,ee_rho,X2,Y2,hh_rho)
      call wall_time(end_t)
      t = end_t - start_t
      ! Done with eigenvectors and kernel
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for pp integrals =',t,' seconds'
      write(*,*)

      deallocate(X1,Y1,X2,Y2,pp_Gam)

      ! Convergence criteria
      err_2b = max(err_eh,err_ee,err_hh)
      
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
