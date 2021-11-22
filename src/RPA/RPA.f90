subroutine RPA(TDA,doACFDT,exchange_kernel,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform a direct random phase approximation calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: TDA
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta
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
  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)
  double precision,allocatable  :: XmY(:,:,:)

  double precision              :: rho
  double precision              :: EcRPA(nspin)
  double precision              :: EcAC(nspin)

! Hello world

  write(*,*)
  write(*,*)'***********************************************'
  write(*,*)'|  Random-phase approximation calculation     |'
  write(*,*)'***********************************************'
  write(*,*)

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Initialization

  EcRPA(:) = 0d0
  EcAC(:)  = 0d0

! Memory allocation

  allocate(Omega(nS,nspin),XpY(nS,nS,nspin),XmY(nS,nS,nspin))

! Singlet manifold

  if(singlet) then 

    ispin = 1

    call linear_response(ispin,.true.,TDA,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI,rho,Omega(:,ispin), &
                         EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('RPA@HF       ',ispin,nS,Omega(:,ispin))
    call print_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

  endif

! Triplet manifold 

  if(triplet) then 

    ispin = 2

    call linear_response(ispin,.true.,TDA,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI,rho,Omega(:,ispin), &
                         EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('RPA@HF      ',ispin,nS,Omega(:,ispin))
    call print_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

  endif

! if(exchange_kernel) then

!   EcRPA(1) = 0.5d0*EcRPA(1)
!   EcRPA(2) = 1.5d0*EcRPA(2)

! end if

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@RPA  correlation energy (singlet) =',EcRPA(1)
  write(*,'(2X,A50,F20.10)') 'Tr@RPA  correlation energy (triplet) =',EcRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@RPA  correlation energy           =',EcRPA(1) + EcRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@RPA  total energy                 =',ENuc + ERHF + EcRPA(1) + EcRPA(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 
! Switch off ACFDT for RPA as the trace formula is equivalent 

  if(doACFDT) then

    write(*,*) '------------------------------------------------------'
    write(*,*) 'Adiabatic connection version of RPA correlation energy'
    write(*,*) '------------------------------------------------------'
    write(*,*) 

    call ACFDT(exchange_kernel,.false.,.true.,.false.,TDA,.false.,singlet,triplet,eta, &
               nBas,nC,nO,nV,nR,nS,ERI,eHF,eHF,EcAC)

    if(exchange_kernel) then
    
      EcAC(1) = 0.5d0*EcAC(1)
      EcAC(2) = 1.5d0*EcAC(2)
    
    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'AC@RPA  correlation energy (singlet) =',EcAC(1)
    write(*,'(2X,A50,F20.10)') 'AC@RPA  correlation energy (triplet) =',EcAC(2)
    write(*,'(2X,A50,F20.10)') 'AC@RPA  correlation energy           =',EcAC(1) + EcAC(2)
    write(*,'(2X,A50,F20.10)') 'AC@RPA  total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)


  end if

end subroutine RPA
