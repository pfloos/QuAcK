subroutine ppRPA(TDA,doACFDT,singlet,triplet,nBas,nC,nO,nV,nR,ENuc,EHF,ERI,dipole_int,e)

! Perform ppRPA calculation

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
  double precision,intent(in)   :: EHF
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin
  integer                       :: nOO
  integer                       :: nVV
  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)
  double precision,allocatable  :: Om1(:)
  double precision,allocatable  :: X1(:,:)
  double precision,allocatable  :: Y1(:,:)
  double precision,allocatable  :: Om2(:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: Y2(:,:)

  double precision              :: EcRPA(nspin)
  double precision              :: EcAC(nspin)

! Hello world

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'|  particle-particle RPA calculation   |'
  write(*,*)'****************************************'
  write(*,*)

! Initialization

  EcRPA(:) = 0d0
  EcAC(:)  = 0d0

! Singlet manifold

  if(singlet) then 

    write(*,*) '****************'
    write(*,*) '*** Singlets ***'
    write(*,*) '****************'
    write(*,*)

    ispin = 1

    nOO = nO*(nO+1)/2
    nVV = nV*(nV+1)/2

    allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO), &
             Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))

    if(.not.TDA) call ppLR_B(ispin,nBas,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)
                 call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVV,1d0,e,ERI,Cpp)
                 call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOO,1d0,e,ERI,Dpp)

    call ppLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA(ispin))

!   call print_transition_vectors_pp(.true.,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

    call print_excitation_energies('ppRPA@HF (N+2)',ispin,nVV,Om1)
    call print_excitation_energies('ppRPA@HF (N-2)',ispin,nOO,Om2)

    deallocate(Om1,X1,Y1,Om2,X2,Y2,Bpp,Cpp,Dpp)

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

    allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO), &
             Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))

    if(.not.TDA) call ppLR_B(ispin,nBas,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)
                 call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVV,1d0,e,ERI,Cpp)
                 call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOO,1d0,e,ERI,Dpp)

    call ppLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA(ispin))

!   call print_transition_vectors_pp(.false.,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

    call print_excitation_energies('ppRPA@HF (N+2)',ispin,nVV,Om1)
    call print_excitation_energies('ppRPA@HF (N-2)',ispin,nOO,Om2)

    deallocate(Om1,X1,Y1,Om2,X2,Y2,Bpp,Cpp,Dpp)

  endif

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy (singlet) =',EcRPA(1)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy (triplet) =',3d0*EcRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA correlation energy           =',EcRPA(1) + 3d0*EcRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA total energy                 =',ENuc + EHF + EcRPA(1) + 3d0*EcRPA(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then

    write(*,*) '--------------------------------------------------------'
    write(*,*) 'Adiabatic connection version of ppRPA correlation energy'
    write(*,*) '--------------------------------------------------------'
    write(*,*)

    call ppACFDT(TDA,singlet,triplet,nBas,nC,nO,nV,nR,ERI,e,EcAC)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA correlation energy (singlet) =',EcAC(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA correlation energy (triplet) =',EcAC(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA correlation energy           =',EcAC(1) + EcAC(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA total energy                 =',ENuc + EHF + EcAC(1) + EcAC(2),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine
