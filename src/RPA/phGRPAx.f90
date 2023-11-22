subroutine phGRPAx(dotest,TDA,nBas,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)

! Perform random phase approximation calculation with exchange (aka TDHF)

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
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

  double precision              :: EcRPA

! Hello world

  write(*,*)
  write(*,*)'************************************'
  write(*,*)'* Generalized ph-RPAx Calculation  *'
  write(*,*)'************************************'
  write(*,*)

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Initialization

  dRPA  = .false.
  EcRPA = 0d0

! Memory allocation

  allocate(Om(nS),XpY(nS,nS),XmY(nS,nS),Aph(nS,nS),Bph(nS,nS))

  ispin = 3

  call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phLR(TDA,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
  call print_excitation_energies('phRPAx@GHF','spinorbital',nS,Om)
  call phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,Om,XpY,XmY)

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phRPAx correlation energy           = ',EcRPA,' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@phRPAx total energy                 = ',ENuc + EGHF + EcRPA,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

  if(dotest) then

    call dump_test_value('G','phRPAx correlation energy',EcRPA)

  end if

end subroutine 
