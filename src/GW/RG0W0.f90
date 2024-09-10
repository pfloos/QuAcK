subroutine RG0W0(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE,singlet,triplet, & 
                 linearize,eta,regularize,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform G0W0 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)
  double precision,intent(in)   :: eHF(nOrb)

! Local variables

  logical                       :: print_W   = .true.
  logical                       :: plot_self = .false.
  logical                       :: dRPA_W
  integer                       :: isp_W
  double precision              :: lambda
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)

  double precision,allocatable  :: eGWlin(:)
  double precision,allocatable  :: eGW(:)


! Output variables

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted G0W0 Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Initialization

  lambda = 1d0

! Spin manifold and TDA for dynamical screening

  isp_W = 1
  dRPA_W = .true.

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamical screening!'
    write(*,*)
  end if

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),SigC(nOrb),Z(nOrb),Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nOrb,nOrb,nS), & 
           eGW(nOrb),eGWlin(nOrb))

!-------------------!
! Compute screening !
!-------------------!

                 call phLR_A(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,lambda,eHF,ERI,Aph)
  if(.not.TDA_W) call phLR_B(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,lambda,ERI,Bph)

  call phLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

  if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)

!--------------------------!
! Compute spectral weights !
!--------------------------!

  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

!------------------------!
! Compute GW self-energy !
!------------------------!

  if(regularize) call GW_regularization(nOrb,nC,nO,nV,nR,nS,eHF,Om,rho)

  call RGW_self_energy_diag(eta,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigC,Z)

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  ! Linearized or graphical solution?

  eGWlin(:) = eHF(:) + Z(:)*SigC(:)

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGW(:) = eGWlin(:)

  else 

    write(*,*) ' *** Quasiparticle energies obtained by root search *** '
    write(*,*)
  
    call RGW_QP_graph(eta,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,eGWlin,eHF,eGW,Z)

  end if

! Plot self-energy, renormalization factor, and spectral function

  if(plot_self) call RGW_plot_self_energy(nOrb,eta,nC,nO,nV,nR,nS,eHF,eHF,Om,rho)

!--------------------!
! Cumulant expansion !
!--------------------!

! call RGWC(dotest,eta,nOrb,nC,nO,nV,nR,nS,Om,rho,eHF,eHF,eGW,Z)

! Compute the RPA correlation energy

                 call phLR_A(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,lambda,eGW,ERI,Aph)
  if(.not.TDA_W) call phLR_B(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,lambda,ERI,Bph)

  call phLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

!--------------!
! Dump results !
!--------------!

  call print_RG0W0(nOrb,nO,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA,EcGM)

! Perform BSE calculation

  if(dophBSE) then

    call RGW_phBSE(dophBSE2,exchange_kernel,TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta, & 
                   nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0@RHF correlation energy (singlet) = ',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0@RHF correlation energy (triplet) = ',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0@RHF correlation energy           = ',sum(EcBSE),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '-------------------------------------------------------------'
      write(*,*) ' Adiabatic connection version of BSE@G0W0 correlation energy '
      write(*,*) '-------------------------------------------------------------'
      write(*,*) 

      call RGW_phACFDT(exchange_kernel,doXBS,TDA_W,TDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,eHF,eGW,EcBSE)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0W0@RHF correlation energy (singlet) = ',EcBSE(1),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0W0@RHF correlation energy (triplet) = ',EcBSE(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0W0@RHF correlation energy           = ',sum(EcBSE),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0W0@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

  if(doppBSE) then

    call RGW_ppBSE(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0@RHF correlation energy (singlet) = ',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0@RHF correlation energy (triplet) = ',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0@RHF correlation energy           = ',sum(EcBSE),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if
  
! Testing zone

  if(dotest) then

    call dump_test_value('R','G0W0 correlation energy',EcRPA)
    call dump_test_value('R','G0W0 HOMO energy',eGW(nO))
    call dump_test_value('R','G0W0 LUMO energy',eGW(nO+1))

  end if

end subroutine 
