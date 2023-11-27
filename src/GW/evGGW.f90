subroutine evGGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE, & 
                 linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)

! Perform self-consistent eigenvalue-only GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: linear_mixing
  logical                       :: dRPA = .true.
  integer                       :: ispin
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: rcond
  double precision              :: Conv
  double precision              :: EcRPA
  double precision              :: EcBSE
  double precision              :: EcGM
  double precision              :: alpha
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)

! Hello world

  write(*,*)
  write(*,*)'********************************'
  write(*,*)'* Generalized evGW Calculation *'
  write(*,*)'********************************'
  write(*,*)

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

! Linear mixing

  linear_mixing = .false.
  alpha = 0.2d0

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),eGW(nBas),eOld(nBas),Z(nBas),SigC(nBas), & 
           Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nBas,nBas,nS),error_diis(nBas,max_diis),e_diis(nBas,max_diis))

! Initialization

  nSCF            = 0
  ispin           = 3
  n_diis          = 0
  Conv            = 1d0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGW(:)          = eHF(:)
  eOld(:)         = eGW(:)
  Z(:)            = 1d0
  rcond           = 0d0
  
!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

   ! Compute screening

    call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
    if(.not.TDA_W) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

    call phLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

    ! Compute spectral weights

    call GW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY,rho)

    ! Compute correlation part of the self-energy 

    if(regularize) call GW_regularization(nBas,nC,nO,nV,nR,nS,eGW,Om,rho)

    call GGW_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS,eGW,Om,rho,EcGM,SigC,Z)

    ! Solve the quasi-particle equation

    if(linearize) then 
 
      write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
      write(*,*)

      eGW(:) = eHF(:) + SigC(:)

    else 

      write(*,*) ' *** Quasiparticle energies obtained by root search *** '
      write(*,*)
  
      call GGW_QP_graph(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,eOld,eOld,eGW,Z)
 
    end if

    ! Convergence criteria

    Conv = maxval(abs(eGW - eOld))

    ! Print results

    call print_evGGW(nBas,nO,nSCF,Conv,eHF,ENuc,EGHF,SigC,Z,eGW,EcRPA,EcGM)

    ! Linear mixing or DIIS extrapolation

    if(linear_mixing) then
 
      eGW(:) = alpha*eGW(:) + (1d0 - alpha)*eOld(:)
 
    else

      n_diis = min(n_diis+1,max_diis)
      if(abs(rcond) > 1d-7) then
        call DIIS_extrapolation(rcond,nBas,nBas,n_diis,error_diis,e_diis,eGW-eOld,eGW)
      else
        n_diis = 0
      end if

    end if

    ! Save quasiparticles energy for next cycle

    eOld(:) = eGW(:)

    ! Increment

    nSCF = nSCF + 1

  end do
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF+1) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    stop

  end if

! Deallocate memory

  deallocate(Aph,Bph,eOld,Z,SigC,Om,XpY,XmY,rho,error_diis,e_diis)

! Perform BSE calculation

  if(dophBSE) then

    call GGW_phBSE(dophBSE2,TDA_W,TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGW,eGW,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@GHF correlation energy = ',EcBSE,' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@GHF total       energy = ',ENuc + EGHF + EcBSE,' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

!   if(doACFDT) then

!     write(*,*) '------------------------------------------------------'
!     write(*,*) 'Adiabatic connection version of BSE correlation energy'
!     write(*,*) '------------------------------------------------------'
!     write(*,*)

!     if(doXBS) then

!       write(*,*) '*** scaled screening version (XBS) ***'
!       write(*,*)

!     end if

!     call GW_phACFDT(exchange_kernel,doXBS,dRPA,TDA_W,TDA,dophBSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,eGW,eGW,EcBSE)

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy (singlet) =',EcBSE(1)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy (triplet) =',EcBSE(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy           =',EcBSE(1) + EcBSE(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW total energy                 =',ENuc + EGHF + EcBSE(1) + EcBSE(2)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

  end if

! if(doppBSE) then

!   call GW_ppBSE(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE)

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@evGW correlation energy (singlet) =',EcBSE(1)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@evGW correlation energy (triplet) =',3d0*EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@evGW correlation energy =',EcBSE(1) + 3d0*EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@evGW total energy =',ENuc + EGHF + EcBSE(1) + 3d0*EcBSE(2)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

! end if

! Testing zone

  if(dotest) then

    call dump_test_value('G','evGW correlation energy',EcRPA)
    call dump_test_value('G','evGW HOMO energy',eGW(nO))
    call dump_test_value('G','evGW LUMO energy',eGW(nO+1))

  end if

end subroutine 
