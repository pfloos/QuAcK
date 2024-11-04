subroutine ppRRPA(dotest,TDA,doACFDT,singlet,triplet,nBas,nC,nO,nV,nR,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform ppRPA calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

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
  double precision,intent(in)   :: eHF(nBas)
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

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Restricted pp-RPA Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! Initialization

  EcRPA(:) = 0d0

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
                 call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVV,1d0,eHF,ERI,Cpp)
                 call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOO,1d0,eHF,ERI,Dpp)

    call ppLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA(ispin))

    call ppLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

!    call print_excitation_energies('ppRPA@RHF','2p (singlet)',nVV,Om1)
!    call print_excitation_energies('ppRPA@RHF','2h (singlet)',nOO,Om2)

    deallocate(Om1,X1,Y1,Om2,X2,Y2,Bpp,Cpp,Dpp)

  end if

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
                 call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVV,1d0,eHF,ERI,Cpp)
                 call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOO,1d0,eHF,ERI,Dpp)

    call ppLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA(ispin))

    call ppLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

!    call print_excitation_energies('ppRPA@RHF','2p (triplet)',nVV,Om1)
!    call print_excitation_energies('ppRPA@RHF','2h (triplet)',nOO,Om2)

    deallocate(Om1,X1,Y1,Om2,X2,Y2,Bpp,Cpp,Dpp)

  end if

  EcRPA(2) = 3d0*EcRPA(2)

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF correlation energy (singlet) = ',EcRPA(1),'au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF correlation energy (triplet) = ',EcRPA(2),'au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF correlation energy           = ',sum(EcRPA),'au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF total       energy           = ',ENuc + ERHF + sum(EcRPA),'au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Compute the correlation energy via the adiabatic connection 

  if(doACFDT) then

    call ppACFDT(TDA,singlet,triplet,nBas,nC,nO,nV,nR,ERI,eHF,EcRPA)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA@RHF correlation energy (singlet) = ',EcRPA(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA@RHF correlation energy (triplet) = ',EcRPA(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA@RHF correlation energy           = ',EcRPA(1) + EcRPA(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'AC@ppRPA@RHF total       energy           = ',ENuc + ERHF + EcRPA(1) + EcRPA(2),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

  if(dotest) then

    call dump_test_value('R','ppRPA correlation energy',sum(EcRPA))

  end if

end subroutine
