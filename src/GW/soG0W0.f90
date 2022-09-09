subroutine soG0W0(doACFDT,exchange_kernel,doXBS,COHSEX,BSE,TDA_W,TDA,dBSE,dTDA,evDyn,ppBSE, & 
                  singlet,triplet,linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,   & 
                  ERI_AO,ERI_MO,dipole_int,PHF,cHF,eHF,Vxc,eGW)

! Perform G0W0 calculation in the spin-orbital basis

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: BSE
  logical,intent(in)            :: ppBSE
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
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
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_MO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: Vxc(nBas)
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: PHF(nBas,nBas)

! Local variables

  logical                       :: print_W = .true.
  integer                       :: ispin
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcppBSE(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

  double precision,allocatable  :: eGWlin(:)

  integer                       :: nBas2
  integer                       :: nC2
  integer                       :: nO2
  integer                       :: nV2
  integer                       :: nR2
  integer                       :: nS2

  double precision,allocatable  :: seHF(:),seGW(:),sERI(:,:,:,:)

! Output variables

  double precision              :: eGW(nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|        One-shot soG0W0 calculation           |'
  write(*,*)'************************************************'
  write(*,*)

! Initialization

  EcRPA = 0d0

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

! spatial to spin transformation

  nBas2 = 2*nBas
  nO2   = 2*nO
  nV2   = 2*nV
  nC2   = 2*nC
  nR2   = 2*nR
  nS2   = nO2*nV2

  allocate(seHF(nBas2),seGW(nBas2),sERI(nBas2,nBas2,nBas2,nBas2))

  call spatial_to_spin_MO_energy(nBas,eHF,nBas2,seHF)
  call spatial_to_spin_MO_energy(nBas,eGW,nBas2,seGW)
  call spatial_to_spin_ERI(nBas,ERI_MO,nBas2,sERI)

! Spin manifold 

  ispin = 3

! Memory allocation

  allocate(SigC(nBas2),Z(nBas2),OmRPA(nS2),XpY_RPA(nS2,nS2),XmY_RPA(nS2,nS2),rho_RPA(nBas2,nBas2,nS2),eGWlin(nBas2))

!-------------------!
! Compute screening !
!-------------------!

  call linear_response(ispin,.true.,TDA_W,eta,nBas2,nC2,nO2,nV2,nR2,nS2,1d0, & 
                       seHF,sERI,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

  if(print_W) call print_excitation('RPA@HF      ',ispin,nS2,OmRPA)

!--------------------------!
! Compute spectral weights !
!--------------------------!

  call excitation_density(nBas2,nC2,nO2,nR2,nS2,sERI,XpY_RPA,rho_RPA)

!------------------------!
! Compute GW self-energy !
!------------------------!

  if(regularize) then 

    call regularized_self_energy_correlation_diag(COHSEX,eta,nBas,nC,nO,nV,nR,nS,eHF,OmRPA,rho_RPA,EcGM,SigC)
    call regularized_renormalization_factor(COHSEX,eta,nBas,nC,nO,nV,nR,nS,eHF,OmRPA,rho_RPA,Z)

  else

    call self_energy_correlation_diag_so(COHSEX,eta,nBas2,nC2,nO2,nV2,nR2,nS2,seHF,OmRPA,rho_RPA,EcGM,SigC)
    call renormalization_factor_so(COHSEX,eta,nBas2,nC2,nO2,nV2,nR2,nS2,seHF,OmRPA,rho_RPA,Z)

  end if

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  eGWlin(:) = seHF(:) + Z(:)*SigC(:)

  ! Linearized or graphical solution?

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    seGW(:) = eGWlin(:)

  end if

! Compute the RPA correlation energy

  call linear_response(ispin,.true.,TDA_W,eta,nBas2,nC2,nO2,nV2,nR2,nS2,1d0,seGW,sERI, & 
                       EcRPA,OmRPA,XpY_RPA,XmY_RPA)

!--------------!
! Dump results !
!--------------!

  call print_G0W0(nBas2,nO2,seHF,ENuc,ERHF,SigC,Z,seGW,EcRPA,EcGM)

! Deallocate memory

  deallocate(SigC,Z,OmRPA,XpY_RPA,XmY_RPA,rho_RPA,eGWlin)

! Perform BSE calculation

  if(ppBSE) then

    call  Bethe_Salpeter_pp_so(TDA_W,TDA,singlet,triplet,eta,nBas2,nC2,nO2,nV2,nR2,nS2,sERI,dipole_int,seHF,seGW,EcppBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 correlation energy (singlet) =',EcppBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 correlation energy (triplet) =',3d0*EcppBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 correlation energy           =',EcppBSE(1) + 3d0*EcppBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 total energy                 =',ENuc + ERHF + EcppBSE(1) + 3d0*EcppBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine soG0W0
