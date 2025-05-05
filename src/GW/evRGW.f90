subroutine evRGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE, & 
                 singlet,triplet,linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform self-consistent eigenvalue-only GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
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
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  logical                       :: dRPA = .true.
  integer                       :: ispin
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: flow
  double precision              :: rcond
  double precision              :: Conv
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
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
  write(*,*)'*******************************'
  write(*,*)'* Restricted evGW Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
    write(*,*)
  end if

! SRG regularization

  flow = 500d0

  if(doSRG) then

    write(*,*) '*** SRG regularized evGW scheme ***'
    write(*,*)

  end if

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),eGW(nOrb),eOld(nOrb),Z(nOrb),SigC(nOrb), & 
           Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nOrb,nOrb,nS),error_diis(nOrb,max_diis),e_diis(nOrb,max_diis))

! Initialization

  nSCF            = 0
  ispin           = 1
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

                   call phRLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
    if(.not.TDA_W) call phRLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

    call phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

    ! Compute spectral weights

    call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

    ! Compute correlation part of the self-energy 

    if(doSRG) then
      call RGW_SRG_self_energy_diag(flow,nBas,nOrb,nC,nO,nV,nR,nS,eGW,Om,rho,EcGM,SigC,Z)
    else
      call RGW_self_energy_diag(eta,nBas,nOrb,nC,nO,nV,nR,nS,eGW,Om,rho,EcGM,SigC,Z)
    end if

    ! Solve the quasi-particle equation

    if(linearize) then 
 
      write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
      write(*,*)

      eGW(:) = eHF(:) + SigC(:)

    else 

      write(*,*) ' *** Quasiparticle energies obtained by root search *** '
      write(*,*)
  
      call RGW_QP_graph(doSRG,eta,flow,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,eOld,eOld,eGW,Z)
 
    end if

    ! Convergence criteria

    Conv = maxval(abs(eGW - eOld))

    ! Print results

    call print_evRGW(nOrb,nO,nSCF,Conv,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA,EcGM)

    ! DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nOrb,nOrb,n_diis,error_diis,e_diis,eGW-eOld,eGW)

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

!--------------------!
! Cumulant expansion !
!--------------------!

! call RGWC(dotest,eta,nOrb,nC,nO,nV,nR,nS,Om,rho,eHF,eGW,eGW,Z)

! Deallocate memory

  deallocate(Aph,Bph,eOld,Z,SigC,Om,XpY,XmY,rho,error_diis,e_diis)

! Perform BSE calculation

  if(dophBSE) then

    call RGW_phBSE(dophBSE2,exchange_kernel,TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta, & 
                   nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eGW,eGW,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@RHF correlation energy           = ',sum(EcBSE),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      call RGW_phACFDT(exchange_kernel,doXBS,TDA_W,TDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,eGW,eGW,EcBSE)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@RHF correlation energy           = ',sum(EcBSE),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

  if(doppBSE) then

    call RGW_ppBSE(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGW@RHF correlation energy           = ',sum(EcBSE),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGW@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Testing zone

  if(dotest) then

    call dump_test_value('R','evGW correlation energy',EcRPA)
    call dump_test_value('R','evGW HOMO energy',eGW(nO))
    call dump_test_value('R','evGW LUMO energy',eGW(nO+1))

  end if

end subroutine 
