subroutine evGW(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,COHSEX,SOSEX,BSE,TDA_W,TDA,  & 
                G0W,GW0,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,        &
                ERI_AO,ERI_MO,dipole_int,PHF,cHF,eHF,Vxc,eG0W0)

! Perform self-consistent eigenvalue-only GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: SOSEX
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: G0W
  logical,intent(in)            :: GW0
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: Vxc(nBas,nspin)
  double precision,intent(in)   :: eG0W0(nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_MO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: linear_mixing
  integer                       :: ispin
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: rcond
  double precision              :: Conv
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM
  double precision              :: alpha
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: SigX(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|       Self-consistent evGW calculation       |'
  write(*,*)'************************************************'
  write(*,*)

! SOSEX correction

  if(SOSEX) then 
    write(*,*) 'SOSEX correction activated but BUG!'
    stop
  end if

! COHSEX approximation

  if(COHSEX) then 
    write(*,*) 'COHSEX approximation activated!'
    write(*,*)
  end if

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
    write(*,*)
  end if

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! GW0

  if(GW0) then 
    write(*,*) 'GW0 scheme activated!'
    write(*,*)
  end if

! G0W

  if(G0W) then 
    write(*,*) 'G0W scheme activated!'
    write(*,*)
  end if

! Linear mixing

  linear_mixing = .false.
  alpha = 0.2d0

! Memory allocation

  allocate(eGW(nBas),eOld(nBas),Z(nBas),SigX(nBas),SigC(nBas),OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS), & 
           rho_RPA(nBas,nBas,nS),error_diis(nBas,max_diis),e_diis(nBas,max_diis))

! Compute the exchange part of the self-energy

  call self_energy_exchange_diag(nBas,cHF,PHF,ERI_AO,SigX)

! Initialization

  nSCF            = 0
  ispin           = 1
  n_diis          = 0
  Conv            = 1d0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGW(:)          = eG0W0(:)
  eOld(:)         = eGW(:)
  Z(:)            = 1d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

   ! Compute screening

    if(.not. GW0 .or. nSCF == 0) then

      call linear_response(ispin,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI_MO,OmRPA, & 
                           rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

    endif

!   Compute spectral weights

    call excitation_density(nBas,nC,nO,nR,nS,ERI_MO,XpY_RPA,rho_RPA)

!   Compute correlation part of the self-energy 

    if(G0W) then

      call self_energy_correlation_diag(COHSEX,eta,nBas,nC,nO,nV,nR,nS,eHF,OmRPA,rho_RPA,EcGM,SigC)
      call renormalization_factor(COHSEX,eta,nBas,nC,nO,nV,nR,nS,eHF,OmRPA,rho_RPA,Z)

    else 

      call self_energy_correlation_diag(COHSEX,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA,rho_RPA,EcGM,SigC)
      call renormalization_factor(COHSEX,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA,rho_RPA,Z)

    endif

    ! Solve the quasi-particle equation

    eGW(:) = eHF(:) + SigC(:)

    ! Convergence criteria

    Conv = maxval(abs(eGW - eOld))

    ! Print results

    call print_evGW(nBas,nO,nSCF,Conv,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA,EcGM)

    ! Linear mixing or DIIS extrapolation

    if(linear_mixing) then
 
      eGW(:) = alpha*eGW(:) + (1d0 - alpha)*eOld(:)
 
    else

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBas,nBas,n_diis,error_diis,e_diis,eGW-eOld,eGW)

!    Reset DIIS if required

      if(abs(rcond) < 1d-15) n_diis = 0

    endif

    ! Save quasiparticles energy for next cycle

    eOld(:) = eGW(:)

    ! Increment

    nSCF = nSCF + 1

  enddo
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Plot stuff

! call plot_GW(nBas,nC,nO,nV,nR,nS,eHF,eGW,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin))

! Did it actually converge?

  if(nSCF == maxSCF+1) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  endif

! Deallocate memory

  deallocate(eOld,Z,SigC,OmRPA,XpY_RPA,XmY_RPA,rho_RPA,error_diis,e_diis)

! Perform BSE calculation

  if(BSE) then

    call Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eGW,eGW,EcBSE)

    if(exchange_kernel) then

      EcBSE(1) = 0.5d0*EcBSE(1)
      EcBSE(2) = 1.5d0*EcBSE(2)

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW correlation energy           =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '------------------------------------------------------'
      write(*,*) 'Adiabatic connection version of BSE correlation energy'
      write(*,*) '------------------------------------------------------'
      write(*,*)

      if(doXBS) then

        write(*,*) '*** scaled screening version (XBS) ***'
        write(*,*)

      end if

      call ACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,eGW,eGW,EcAC)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy (singlet) =',EcAC(1)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy (triplet) =',EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy           =',EcAC(1) + EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  endif

end subroutine evGW
