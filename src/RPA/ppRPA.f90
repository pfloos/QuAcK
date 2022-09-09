subroutine ppRPA(TDA,doACFDT,singlet,triplet,nBas,nC,nO,nV,nR,ENuc,ERHF,ERI,dipole_int,e)

! Perform pp-RPA calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin
  integer                       :: nS
  integer                       :: nOO
  integer                       :: nVV
  double precision,allocatable  :: Omega1(:)
  double precision,allocatable  :: X1(:,:)
  double precision,allocatable  :: Y1(:,:)
  double precision,allocatable  :: Omega2(:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: Y2(:,:)

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
  EcAC(:)     = 0d0

! Useful quantities

  nS = nO*nV

! Singlet manifold

  if(singlet) then 

    write(*,*) '****************'
    write(*,*) '*** Singlets ***'
    write(*,*) '****************'
    write(*,*)

    ispin = 1

    nOO = nO*(nO+1)/2
    nVV = nV*(nV+1)/2

    allocate(Omega1(nVV),X1(nVV,nVV),Y1(nOO,nVV),Omega2(nOO),X2(nVV,nOO),Y2(nOO,nOO))

    call linear_response_pp(ispin,TDA,nBas,nC,nO,nV,nR,nOO,nVV,1d0,e,ERI, & 
                            Omega1,X1,Y1,Omega2,X2,Y2,Ec_ppRPA(ispin))

    call print_transition_vectors_pp(.true.,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Omega1,X1,Y1,Omega2,X2,Y2)

!   call print_excitation('pp-BSE (N+2)',ispin,nVV,Omega1)
!   call print_excitation('pp-BSE (N-2)',ispin,nOO,Omega2)

    deallocate(Omega1,X1,Y1,Omega2,X2,Y2)

  endif

! Triplet manifold 

  if(triplet) then 

    write(*,*) '****************'
    write(*,*) '*** Triplets ***'
    write(*,*) '****************'
    write(*,*)

    ispin = 2

    nOO = nO*(nO-1)/2
    nVV = nV*(nV-1)/2

    allocate(Omega1(nVV),X1(nVV,nVV),Y1(nOO,nVV),Omega2(nOO),X2(nVV,nOO),Y2(nOO,nOO))

    call linear_response_pp(ispin,TDA,nBas,nC,nO,nV,nR,nOO,nVV,1d0,e,ERI, &
                            Omega1,X1,Y1,Omega2,X2,Y2,Ec_ppRPA(ispin))

    call print_transition_vectors_pp(.false.,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Omega1,X1,Y1,Omega2,X2,Y2)

!   call print_excitation('pp-BSE (N+2)',ispin,nVV,Omega1)
!   call print_excitation('pp-BSE (N-2)',ispin,nOO,Omega2)

    deallocate(Omega1,X1,Y1,Omega2,X2,Y2)

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

    call ACFDT_pp(TDA,singlet,triplet,nBas,nC,nO,nV,nR,nS,ERI,e,EcAC)

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
