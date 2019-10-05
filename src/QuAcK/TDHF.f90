subroutine TDHF(singlet_manifold,triplet_manifold,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,e)

! Perform random phase approximation calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dRPA
  logical                       :: TDA
  logical                       :: BSE
  integer                       :: ispin
  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)

  double precision              :: rho
  double precision              :: EcRPA(nspin)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|  Time-dependent Hartree-Fock calculation     |'
  write(*,*)'************************************************'
  write(*,*)

! Initialization

  EcRPA(:) = 0d0

! Switch on exchange for TDHF

  dRPA = .false.
 
! Switch off Tamm-Dancoff approximation for TDHF

  TDA = .false.
 
! Switch off Bethe-Salpeter equation for TDHF

  BSE = .false. 

! Memory allocation

  allocate(Omega(nS,nspin),XpY(nS,nS,nspin))

! Singlet manifold

  if(singlet_manifold) then 

    ispin = 1

    call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI,rho, &
                         EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin))
    call print_excitation('TDHF ',ispin,nS,Omega(:,ispin))

    call linear_response_ph(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI,rho, &
                         EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin))
    call print_excitation('p-h TDHF ',ispin,nS,Omega(:,ispin))

    call linear_response_pp(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI,rho)

  endif

! Triplet manifold 

  if(triplet_manifold) then 

    ispin = 2

    call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI,rho, &
                         EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin))
    call print_excitation('TDHF ',ispin,nS,Omega(:,ispin))

  endif

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A40,F15.6)') 'RPA@TDHF correlation energy (singlet) =',EcRPA(1)
  write(*,'(2X,A40,F15.6)') 'RPA@TDHF correlation energy (triplet) =',EcRPA(2)
  write(*,'(2X,A40,F15.6)') 'RPA@TDHF correlation energy           =',EcRPA(1) + EcRPA(2)
  write(*,'(2X,A40,F15.6)') 'RPA@TDHF total energy                 =',ENuc + ERHF + EcRPA(1) + EcRPA(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine TDHF
