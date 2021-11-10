subroutine ppRPA(TDA,doACFDT,exchange_kernel,singlet,triplet,eta,nBas,nC,nO,nV,nR,ENuc,ERHF,ERI,e)

! Perform pp-RPA calculation

  implicit none
  include 'parameters.h'

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
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin
  integer                       :: nS
  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  double precision,allocatable  :: Omega1s(:),Omega1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: Omega2s(:),Omega2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)

  double precision              :: Ec_ppRPA(nspin)
  double precision              :: EcAC(nspin)

! Hello world

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'|  particle-particle RPA calculation   |'
  write(*,*)'****************************************'
  write(*,*)

! Initialization

  Ec_ppRPA(:) = 0d0
  EcAC(:)   = 0d0

! Useful quantities

  nS   = nO*nV

  nOOs = nO*(nO+1)/2
  nVVs = nV*(nV+1)/2

  nOOt = nO*(nO-1)/2
  nVVt = nV*(nV-1)/2

 ! Memory allocation

  allocate(Omega1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs), & 
           Omega2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs))

  allocate(Omega1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt), & 
           Omega2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt))

! Singlet manifold

  if(singlet) then 

    ispin = 1

    call linear_response_pp(ispin,TDA,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,e,ERI, & 
                            Omega1s,X1s,Y1s,Omega2s,X2s,Y2s,Ec_ppRPA(ispin))

    call print_excitation('pp-RPA (N+2)',ispin,nVVs,Omega1s)
    call print_excitation('pp-RPA (N-2)',ispin,nOOs,Omega2s)

  endif

! Triplet manifold 

  if(triplet) then 

    ispin = 2

    call linear_response_pp(ispin,TDA,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,e,ERI, &
                            Omega1t,X1t,Y1t,Omega2t,X2t,Y2t,Ec_ppRPA(ispin))

    call print_excitation('pp-RPA (N+2)',ispin,nVVt,Omega1t)
    call print_excitation('pp-RPA (N-2)',ispin,nOOt,Omega2t)

  endif

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy (singlet) =',Ec_ppRPA(1)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy (triplet) =',3d0*Ec_ppRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy           =',Ec_ppRPA(1) + 3d0*Ec_ppRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA total energy                 =',ENuc + ERHF + Ec_ppRPA(1) + 3d0*Ec_ppRPA(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then

    write(*,*) '---------------------------------------------------------'
    write(*,*) 'Adiabatic connection version of pp-RPA correlation energy'
    write(*,*) '---------------------------------------------------------'
    write(*,*)

    call ACFDT_Tmatrix(exchange_kernel,.false.,.false.,.false.,TDA,.false.,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS, &
                       ERI,e,e,EcAC)

    if(exchange_kernel) then

      EcAC(1) = 0.5d0*EcAC(1)
      EcAC(2) = 1.5d0*EcAC(1)

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA correlation energy (singlet) =',EcAC(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA correlation energy (triplet) =',EcAC(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA correlation energy           =',EcAC(1) + EcAC(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if


end subroutine ppRPA
