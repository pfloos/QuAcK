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

  integer                       :: n_it_1b,n_it_2b
  double precision              :: err_1b,err_2b
  double precision              :: err_eh, err_pp
  double precision              :: err_eig_eh, err_eig_hh, err_eig_ee
  double precision              :: start_t, end_t, t
  double precision              :: start_1b, end_1b, t_1b
  double precision              :: start_2b, end_2b, t_2b

  integer                       :: nOO,nVV

  ! eh BSE
  double precision              :: EcRPA
  double precision,allocatable  :: Aph(:,:), Bph(:,:)
  double precision,allocatable  :: XpY(:,:), XmY(:,:)
  double precision,allocatable  :: eh_Om(:), old_eh_Om(:)
  double precision,allocatable  :: eh_Gam_A(:,:),eh_Gam_B(:,:)
  ! pp BSE
  double precision,allocatable  :: Bpp(:,:), Cpp(:,:), Dpp(:,:)
  double precision,allocatable  :: X1(:,:),Y1(:,:)
  double precision,allocatable  :: ee_Om(:), old_ee_Om(:)
  double precision,allocatable  :: X2(:,:),Y2(:,:)
  double precision,allocatable  :: hh_Om(:), old_hh_Om(:)
  double precision,allocatable  :: pp_Gam_B(:,:),pp_Gam_C(:,:),pp_Gam_D(:,:)
  ! Effective integrals
  double precision,allocatable  :: eh_rho(:,:,:), ee_rho(:,:,:), hh_rho(:,:,:)
  ! Reducible kernels
  double precision,allocatable  :: eh_Phi(:,:,:,:), pp_Phi(:,:,:,:)
  double precision,allocatable  :: old_eh_Phi(:,:,:,:), old_pp_Phi(:,:,:,:)
  ! One-body quantities
  double precision,allocatable  :: eQPlin(:),eQP(:),eOld(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision              :: EcGM
  ! DIIS
  integer                       :: max_diis,n_diis
  double precision              :: rcond
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: Om_diis(:,:)
  double precision,allocatable  :: err(:)
  double precision,allocatable  :: Om(:)
  double precision              :: alpha

! Output variables
! None
  
! Useful parameters
  nOO = nO*(nO - 1)/2
  nVV = nV*(nV - 1)/2

  allocate(eQP(nOrb),eOld(nOrb))
    
! DIIS parameters

  ! max_diis = 10
  ! n_diis   = 0
  ! rcond    = 0d0

  ! allocate(err_diis(nS+nOO+nVV,max_diis),Om_diis(nS+nOO+nVV,max_diis))
  ! allocate(err(nS+nOO+nVV),Om(nS+nOO+nVV))

  ! err_diis(:,:) = 0d0
  ! Om_diis(:,:)  = 0d0

! Start
 
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
  allocate(eh_rho(nOrb,nOrb,nS+nS),ee_rho(nOrb,nOrb,nVV),hh_rho(nOrb,nOrb,nOO))
  allocate(old_eh_Phi(nOrb,nOrb,nOrb,nOrb),old_pp_Phi(nOrb,nOrb,nOrb,nOrb))

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

  old_eh_Om(:) = 0d0
  old_ee_Om(:) = 0d0
  old_hh_Om(:) = 0d0

  old_eh_Phi(:,:,:,:) = 0d0
  old_pp_Phi(:,:,:,:) = 0d0
    
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

                     call G_eh_Gamma_A(nOrb,nC,nO,nR,nS,old_eh_Phi,old_pp_Phi,eh_Gam_A)
        if(.not.TDA) call G_eh_Gamma_B(nOrb,nC,nO,nR,nS,old_eh_Phi,old_pp_Phi,eh_Gam_B)
        
      end if
      
       Aph(:,:) = Aph(:,:) + eh_Gam_A(:,:)
       Bph(:,:) = Bph(:,:) + eh_Gam_B(:,:)    
      

      call phGLR(TDA,nS,Aph,Bph,EcRPA,eh_Om,XpY,XmY)

      call wall_time(end_t)

      t = end_t - start_t
      write(*,'(A50,1X,F9.3,A8)') 'Wall time for phBSE =',t,' seconds'
      write(*,*)

      if(print_phLR) call print_excitation_energies('phBSE@Parquet','eh generalized',nS,eh_Om)

      err_eig_eh = maxval(abs(old_eh_Om - eh_Om))

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

        if(.not.TDA) call G_pp_Gamma_B(nOrb,nC,nO,nR,nOO,nVV,old_eh_Phi,pp_Gam_B)
        call G_pp_Gamma_C(nOrb,nO,nR,nVV,old_eh_Phi,pp_Gam_C)
        call G_pp_Gamma_D(nOrb,nC,nO,nOO,old_eh_Phi,pp_Gam_D)

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

      err_eig_ee = maxval(abs(old_ee_Om - ee_Om))
      err_eig_hh = maxval(abs(old_hh_Om - hh_Om))

      deallocate(Bpp,Cpp,Dpp,pp_Gam_B,pp_Gam_C,pp_Gam_D)

            
      write(*,*) '----------------------------------------'
      write(*,*) ' Two-body (eigenvalue) convergence '
      write(*,*) '----------------------------------------'
      write(*,'(1X,A30,F10.6)')'Error for eh channel = ',err_eig_eh
      write(*,'(1X,A30,F10.6)')'Error for pp channel = ',max(err_eig_ee,err_eig_hh)
      write(*,*) '----------------------------------------'
      write(*,*)

      !--------------------!
      ! DIIS extrapolation !
      !--------------------!

      ! err(       1:nS        ) = eh_Om(:) - old_eh_Om(:)
      ! err(    nS+1:nS+nVV    ) = ee_Om(:) - old_ee_Om(:)
      ! err(nVV+nS+1:nS+nVV+nOO) = hh_Om(:) - old_hh_Om(:) 

      ! Om(       1:nS        ) = eh_Om(:)
      ! Om(    nS+1:nS+nVV    ) = ee_Om(:)
      ! Om(nVV+nS+1:nS+nVV+nOO) = hh_Om(:)

      ! if(max_diis > 1) then
     
      !   n_diis = min(n_diis+1,max_diis)
      !   call DIIS_extrapolation(rcond,nS+nOO+nVV,nS+nOO+nVV,n_diis,err_diis,Om_diis,err,Om)
     
      ! end if

      ! eh_Om(:) = Om(       1:nS        )
      ! ee_Om(:) = Om(    nS+1:nS+nVV    )
      ! hh_Om(:) = Om(nVV+nS+1:nS+nVV+nOO)

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
      allocate(eh_rho(nOrb,nOrb,nS+nS))
      allocate(ee_rho(nOrb,nOrb,nVV),hh_rho(nOrb,nOrb,nOO))

      ! Build singlet eh integrals
      write(*,*) 'Computing eh screened integrals...'

      call wall_time(start_t)
      call G_eh_screened_integral(nOrb,nC,nO,nR,nS,ERI,old_eh_Phi,old_pp_Phi,XpY,XmY,eh_rho)
      call wall_time(end_t)
      t = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for eh integrals =',t,' seconds'
      write(*,*)
      ! Done with eigenvectors and kernel
      deallocate(XpY,XmY)

      ! Build singlet pp integrals
      write(*,*) 'Computing pp screened integrals...'

      call wall_time(start_t)
      call G_pp_screened_integral(nOrb,nC,nO,nR,nOO,nVV,ERI,old_eh_Phi,X1,Y1,ee_rho,X2,Y2,hh_rho)
      call wall_time(end_t)
      t = end_t - start_t
      ! Done with eigenvectors and kernel
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for pp integrals =',t,' seconds'
      write(*,*)

      deallocate(X1,Y1,X2,Y2)

      !----------------------------!
      ! Compute reducible kernels  !
      !----------------------------!

      ! Memory allocation
      allocate(eh_Phi(nOrb,nOrb,nOrb,nOrb))
      allocate(pp_Phi(nOrb,nOrb,nOrb,nOrb))

      ! Build eh reducible kernels
      write(*,*) 'Computing eh reducible kernel ...'

      call wall_time(start_t)
      call G_eh_Phi(nOrb,nC,nR,nS,old_eh_Om,eh_rho,eh_Phi)
      call wall_time(end_t)
      t = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for eh reducible kernel =',t,' seconds'
      write(*,*)

      ! Build pp reducible kernels
      write(*,*) 'Computing pp reducible kernel ...'

      call wall_time(start_t)
      call G_pp_Phi(nOrb,nC,nR,nOO,nVV,old_ee_Om,ee_rho,old_hh_Om,hh_rho,pp_Phi)
      call wall_time(end_t)
      t = end_t - start_t
      write(*,'(1X,A50,1X,F9.3,A8)') 'Wall time for pp reducible kernel =',t,' seconds'
      write(*,*)

      ! alpha = 0.01d0
      ! eh_Phi(:,:,:,:) = alpha * eh_Phi(:,:,:,:) + (1d0 - alpha) * old_eh_Phi(:,:,:,:)
      ! pp_Phi(:,:,:,:) = alpha * pp_Phi(:,:,:,:) + (1d0 - alpha) * old_pp_Phi(:,:,:,:)

      err_eh = maxval(abs(old_eh_Phi - eh_Phi))
      err_pp = maxval(abs(old_pp_Phi - pp_Phi))

      old_eh_Phi(:,:,:,:) = eh_Phi(:,:,:,:)
      old_pp_Phi(:,:,:,:) = pp_Phi(:,:,:,:)
      
      ! Free memory
      deallocate(eh_Phi,pp_Phi)
      
      !--------------------!
      ! DIIS extrapolation !
      !--------------------!


      write(*,*) '----------------------------------------'
      write(*,*) ' Two-body (kernel) convergence '
      write(*,*) '----------------------------------------'
      write(*,'(1X,A30,F10.6)')'Error for eh channel = ',err_eh
      write(*,'(1X,A30,F10.6)')'Error for pp channel = ',err_pp
      write(*,*) '----------------------------------------'
      write(*,*)

      
      ! Convergence criteria
      err_2b = max(err_eh,err_pp)
      
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
