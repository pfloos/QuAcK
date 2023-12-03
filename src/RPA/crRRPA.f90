subroutine crRRPA(dotest,TDA,doACFDT,exchange_kernel,singlet,triplet,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Crossed-ring channel of the random phase approximation 

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: dotest

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
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
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

  double precision              :: EcRPA(nspin)

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Restricted cr-RPA Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Initialization

  dRPA     = .false.
  EcRPA(:) = 0d0

! Memory allocation

  allocate(Om(nS),XpY(nS,nS),XmY(nS,nS),Aph(nS,nS),Bph(nS,nS))

! Singlet manifold

  if(singlet) then 

    ispin = 1

    call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,-1d0,eHF,ERI,Aph)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,-1d0,ERI,Bph)

    call phLR(TDA,nS,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call print_excitation_energies('crRPA@RHF','singlet',nS,Om)
    call phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,XpY,XmY)

  end if

! Triplet manifold 

  if(triplet) then 

    ispin = 2

    call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,-1d0,eHF,ERI,Aph)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,-1d0,ERI,Bph)

    call phLR(TDA,nS,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)
    call print_excitation_energies('crRPA@RHF','triplet',nS,Om)
    call phLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,XpY,XmY)

  end if

  if(exchange_kernel) then

    EcRPA(1) = 0.5d0*EcRPA(1)
    EcRPA(2) = 1.5d0*EcRPA(2)

  end if

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@crRRPA correlation energy (singlet) =',EcRPA(1),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@crRRPA correlation energy (triplet) =',EcRPA(2),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@crRRPA correlation energy           =',sum(EcRPA),' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@crRRPA total energy                 =',ENuc + ERHF + sum(EcRPA),' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then

    write(*,*) '-------------------------------------------------------'
    write(*,*) 'Adiabatic connection version of crRPA correlation energy'
    write(*,*) '-------------------------------------------------------'
    write(*,*)

    call crACFDT(exchange_kernel,dRPA,TDA,singlet,triplet,nBas,nC,nO,nV,nR,nS,ERI,eHF,EcRPA)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'AC@crRRPA correlation energy (singlet) =',EcRPA(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@crRRPA correlation energy (triplet) =',EcRPA(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@crRRPA correlation energy           =',sum(EcRPA),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@crRRPA total energy                 =',ENuc + ERHF + sum(EcRPA),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

  if(dotest) then
  
    call dump_test_value('R','crRPA correlation energy',sum(EcRPA))
  
  end if

end subroutine 
