subroutine crRPA(TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,EHF,ERI,dipole_int,e)

! Crossed-ring channel of the random phase approximation 

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: TDA
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin
  logical                       :: dRPA
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)

  double precision              :: EcTr(nspin)
  double precision              :: EcAC(nspin)

! Hello world

  write(*,*)
  write(*,*)'***********************************************************'
  write(*,*)'|  Random phase approximation calculation: cr channel     |'
  write(*,*)'***********************************************************'
  write(*,*)

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Initialization

  dRPA = .false.

  EcTr(:) = 0d0
  EcAC(:)   = 0d0

! Memory allocation

  allocate(Om(nS),XpY(nS,nS),XmY(nS,nS),Aph(nS,nS))
  if(.not.TDA) allocate(Bph(nS,nS))

! Singlet manifold

  if(singlet) then 

    ispin = 1

    call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,-1d0,e,ERI,Aph)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,-1d0,ERI,Bph)

    call phLR(TDA,nS,Aph,Bph,EcTr(ispin),Om,XpY,XmY)
    call print_excitation_energies('crRPA@HF',ispin,nS,Om)
    call phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,XpY,XmY)

  endif

! Triplet manifold 

  if(triplet) then 

    ispin = 2

    call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,-1d0,e,ERI,Aph)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,-1d0,ERI,Bph)

    call phLR(TDA,nS,Aph,Bph,EcTr(ispin),Om,XpY,XmY)
    call print_excitation_energies('crRPA@HF',ispin,nS,Om)
    call phLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,XpY,XmY)

  endif

  if(exchange_kernel) then

    EcTr(1) = 0.5d0*EcTr(1)
    EcTr(2) = 1.5d0*EcTr(2)

  end if

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@crRPA correlation energy (singlet) =',EcTr(1)
  write(*,'(2X,A50,F20.10)') 'Tr@crRPA correlation energy (triplet) =',EcTr(2)
  write(*,'(2X,A50,F20.10)') 'Tr@crRPA correlation energy           =',EcTr(1) + EcTr(2)
  write(*,'(2X,A50,F20.10)') 'Tr@crRPA total energy                 =',ENuc + EHF + EcTr(1) + EcTr(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then

    write(*,*) '-------------------------------------------------------'
    write(*,*) 'Adiabatic connection version of crRPA correlation energy'
    write(*,*) '-------------------------------------------------------'
    write(*,*)

    call crACFDT(exchange_kernel,dRPA,TDA,singlet,triplet,nBas,nC,nO,nV,nR,nS,ERI,e,EcAC)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'AC@crRPA correlation energy (singlet) =',EcAC(1)
    write(*,'(2X,A50,F20.10)') 'AC@crRPA correlation energy (triplet) =',EcAC(2)
    write(*,'(2X,A50,F20.10)') 'AC@crRPA correlation energy           =',EcAC(1) + EcAC(2)
    write(*,'(2X,A50,F20.10)') 'AC@crRPA total energy                 =',ENuc + EHF + EcAC(1) + EcAC(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine 
