subroutine RPA(doACFDT,exchange_kernel,singlet_manifold,triplet_manifold,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,e)

! Perform a direct random phase approximation calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
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
  write(*,*)'|  random-phase approximation calculation     |'
  write(*,*)'***********************************************'
  write(*,*)

! Initialization

  EcRPA(:) = 0d0
  EcAC(:)  = 0d0

! Memory allocation

  allocate(Omega(nS,nspin),XpY(nS,nS,nspin),XmY(nS,nS,nspin))

! Singlet manifold

  if(singlet_manifold) then 

    ispin = 1

    call linear_response(ispin,.true.,.false.,.false.,nBas,nC,nO,nV,nR,nS,1d0,e,ERI,rho, &
                         EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('RPA   ',ispin,nS,Omega(:,ispin))

  endif

! Triplet manifold 

  if(triplet_manifold) then 

    ispin = 2

    call linear_response(ispin,.true.,.false.,.false.,nBas,nC,nO,nV,nR,nS,1d0,e,ERI,rho, &
                         EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('RPA   ',ispin,nS,Omega(:,ispin))

  endif

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F15.6)') 'Tr@RPA  correlation energy (singlet) =',EcRPA(1)
  write(*,'(2X,A50,F15.6)') 'Tr@RPA  correlation energy (triplet) =',EcRPA(2)
  write(*,'(2X,A50,F15.6)') 'Tr@RPA  correlation energy           =',EcRPA(1) + EcRPA(2)
  write(*,'(2X,A50,F15.6)') 'Tr@RPA  total energy                 =',ENuc + ERHF + EcRPA(1) + EcRPA(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then

    write(*,*) '------------------------------------------------------'
    write(*,*) 'Adiabatic connection version of RPA correlation energy'
    write(*,*) '------------------------------------------------------'
    write(*,*) 

    call ACFDT(exchange_kernel,.false.,.true.,.false.,.false.,singlet_manifold,triplet_manifold, &
               nBas,nC,nO,nV,nR,nS,ERI,e,Omega,XpY,XmY,rho,EcAC)


  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F15.6)') 'AC@RPA  correlation energy (singlet) =',EcAC(1)
  write(*,'(2X,A50,F15.6)') 'AC@RPA  correlation energy (triplet) =',EcAC(2)
  write(*,'(2X,A50,F15.6)') 'AC@RPA  correlation energy           =',EcAC(1) + EcAC(2)
  write(*,'(2X,A50,F15.6)') 'AC@RPA  total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)


  end if

end subroutine RPA
