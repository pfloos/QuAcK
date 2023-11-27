subroutine RG0T0eh(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE, & 
                   singlet,triplet,linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform ehG0T0 calculation

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
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: eHF(nBas)

! Local variables

  logical                       :: print_T = .true.
  logical                       :: dRPA = .false.
  integer                       :: ispin
  integer                       :: isp_W
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Sig(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rhoL(:,:,:)
  double precision,allocatable  :: rhoR(:,:,:)

  double precision,allocatable  :: eGTlin(:)
  double precision,allocatable  :: eGT(:)

  double precision,allocatable  :: KA_sta(:,:)
  double precision,allocatable  :: KB_sta(:,:)
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Restricted G0T0eh Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! Initialization

  EcRPA = 0d0

! TDA for T

  if(TDA_T) then 
    write(*,*) 'Tamm-Dancoff approximation for eh T-matrix!'
    write(*,*)
  end if

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),Sig(nBas),Z(nBas),Om(nS),XpY(nS,nS),XmY(nS,nS), & 
           rhoL(nBas,nBas,nS),rhoR(nBas,nBas,nS),eGT(nBas),eGTlin(nBas))

!---------------------------------
! Compute (triplet) RPA screening 
!---------------------------------

  ispin = 2

  call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  if(.not.TDA_T) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phLR(TDA_T,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

  if(print_T) call print_excitation_energies('phRPA@RHF','triplet',nS,Om)

!--------------------------!
! Compute spectral weights !
!--------------------------!

  call GTeh_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY,XmY,rhoL,rhoR)

!------------------------!
! Compute GW self-energy !
!------------------------!

  if(regularize) call GTeh_regularization(nBas,nC,nO,nV,nR,nS,eHF,Om,rhoL,rhoR)

  call GTeh_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rhoL,rhoR,EcGM,Sig,Z)

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  ! Linearized or graphical solution?

  eGTlin(:) = eHF(:) + Z(:)*Sig(:)

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGT(:) = eGTlin(:) 

  else 

    write(*,*) ' *** Quasiparticle energies obtained by root search *** '
    write(*,*)

    call GTeh_QP_graph(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rhoL,rhoR,eGTlin,eHF,eGT,Z)

  end if

! call GTeh_plot_self_energy(nBas,nC,nO,nV,nR,nS,eHF,eHF,Om,rhoL,rhoR)

! Compute the RPA correlation energy based on the G0T0eh quasiparticle energies

  call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI,Aph)
  if(.not.TDA_T) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phLR(TDA_T,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

!--------------!
! Dump results !
!--------------!

  call print_RG0T0eh(nBas,nO,eHF,ENuc,ERHF,Sig,Z,eGT,EcRPA,EcGM)

! Testing zone

  if(dotest) then

    call dump_test_value('R','G0T0eh correlation energy',EcRPA)
    call dump_test_value('R','G0T0eh HOMO energy',eGT(nO))
    call dump_test_value('R','G0T0eh LUMO energy',eGT(nO+1))

  end if

end subroutine 
