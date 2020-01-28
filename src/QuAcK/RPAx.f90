subroutine RPAx(doACFDT,exchange_kernel,singlet_manifold,triplet_manifold,eta, & 
                nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,e)

! Perform random phase approximation calculation with exchange (aka TDHF)

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: singlet_manifold
  double precision,intent(in)   :: eta
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
  double precision              :: EcRPAx(nspin)
  double precision              :: EcAC(nspin)

! Hello world

  write(*,*)
  write(*,*)'***********************************************************'
  write(*,*)'|  Random phase approximation calculation with exchange   |'
  write(*,*)'***********************************************************'
  write(*,*)

! Initialization

  EcRPAx(:) = 0d0
  EcAC(:)   = 0d0

! Memory allocation

  allocate(Omega(nS,nspin),XpY(nS,nS,nspin),XmY(nS,nS,nspin))

! Singlet manifold

  if(singlet_manifold) then 

    ispin = 1

    call linear_response(ispin,.false.,.false.,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,e,ERI,rho, &
                         EcRPAx(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('RPAx  ',ispin,nS,Omega(:,ispin))

  endif

! Triplet manifold 

  if(triplet_manifold) then 

    ispin = 2

    call linear_response(ispin,.false.,.false.,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,e,ERI,rho, &
                         EcRPAx(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('RPAx  ',ispin,nS,Omega(:,ispin))

  endif

  if(exchange_kernel) then

    EcRPAx(1) = 0.5d0*EcRPAx(1)
    EcRPAx(2) = 1.5d0*EcRPAx(2)

  end if

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@RPAx correlation energy (singlet) =',EcRPAx(1)
  write(*,'(2X,A50,F20.10)') 'Tr@RPAx correlation energy (triplet) =',EcRPAx(2)
  write(*,'(2X,A50,F20.10)') 'Tr@RPAx correlation energy           =',EcRPAx(1) + EcRPAx(2)
  write(*,'(2X,A50,F20.10)') 'Tr@RPAx total energy                 =',ENuc + ERHF + EcRPAx(1) + EcRPAx(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then

    write(*,*) '-------------------------------------------------------'
    write(*,*) 'Adiabatic connection version of RPAx correlation energy'
    write(*,*) '-------------------------------------------------------'
    write(*,*)

    call ACFDT(exchange_kernel,.false.,.false.,.false.,.false.,singlet_manifold,triplet_manifold,eta, &
               nBas,nC,nO,nV,nR,nS,ERI,e,Omega,XpY,XmY,rho,EcAC)

    if(exchange_kernel) then

      EcAC(1) = 0.5d0*EcAC(1)
      EcAC(2) = 1.5d0*EcAC(2)

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'AC@RPAx correlation energy (singlet) =',EcAC(1)
    write(*,'(2X,A50,F20.10)') 'AC@RPAx correlation energy (triplet) =',EcAC(2)
    write(*,'(2X,A50,F20.10)') 'AC@RPAx correlation energy           =',EcAC(1) + EcAC(2)
    write(*,'(2X,A50,F20.10)') 'AC@RPAx total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine RPAx
