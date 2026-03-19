subroutine complex_ppRRPA(dotest,TDA,doACFDT,singlet,triplet,nBas,nC,nO,nV,nR,ENuc,ERHF,ERI,dipole_int,eHF)

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
  complex*16,intent(in)         :: ERHF
  complex*16,intent(in)         :: eHF(nBas)
  complex*16,intent(in)         :: ERI(nBas,nBas,nBas,nBas)
  complex*16,intent(in)         :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin
  integer                       :: nOO
  integer                       :: nVV
  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)
  complex*16,allocatable        :: complex_Bpp(:,:)
  complex*16,allocatable        :: complex_Cpp(:,:)
  complex*16,allocatable        :: complex_Dpp(:,:)
  complex*16,allocatable        :: Om1(:)
  complex*16,allocatable        :: X1(:,:)
  complex*16,allocatable        :: Y1(:,:)
  complex*16,allocatable        :: Om2(:)
  complex*16,allocatable        :: X2(:,:)
  complex*16,allocatable        :: Y2(:,:)

  complex*16                    :: EcRPA(nspin)

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

    if(.not.TDA) call ppRLR_B(ispin,nBas,nC,nO,nV,nR,nOO,nVV,1d0,real(ERI),Bpp)
                 call ppRLR_C(ispin,nBas,nC,nO,nV,nR,nVV,1d0,real(eHF),real(ERI),Cpp)
                 call ppRLR_D(ispin,nBas,nC,nO,nV,nR,nOO,1d0,real(eHF),real(ERI),Dpp)
    
    allocate(complex_Bpp(nVV,nOO),complex_Cpp(nVV,nVV),complex_Dpp(nOO,nOO))
    complex_Bpp = cmplx(Bpp,0d0,kind=8)
    complex_Cpp = cmplx(Cpp,0d0,kind=8)
    complex_Dpp = cmplx(Dpp,0d0,kind=8)
    deallocate(Bpp,Cpp,Dpp)
     
    call complex_ppRLR(TDA,nOO,nVV,complex_Bpp,complex_Cpp,complex_Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA(ispin))

    call print_excitation_energies('Real ppRPA@RHF','2p (singlet)',nVV,real(Om1))
    call print_excitation_energies('Imag ppRPA@RHF','2p (singlet)',nVV,aimag(Om1))
    call print_excitation_energies('Real ppRPA@RHF','2h (singlet)',nOO,real(Om2))
    call print_excitation_energies('Imag ppRPA@RHF','2h (singlet)',nOO,aimag(Om2))

    deallocate(Om1,X1,Y1,Om2,X2,Y2,complex_Bpp,complex_Cpp,complex_Dpp)

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

    if(.not.TDA) call ppRLR_B(ispin,nBas,nC,nO,nV,nR,nOO,nVV,1d0,real(ERI),Bpp)
                 call ppRLR_C(ispin,nBas,nC,nO,nV,nR,nVV,1d0,real(eHF),real(ERI),Cpp)
                 call ppRLR_D(ispin,nBas,nC,nO,nV,nR,nOO,1d0,real(eHF),real(ERI),Dpp)

    allocate(complex_Bpp(nVV,nOO),complex_Cpp(nVV,nVV),complex_Dpp(nOO,nOO))
    complex_Bpp = cmplx(Bpp,0d0,kind=8)
    complex_Cpp = cmplx(Cpp,0d0,kind=8)
    complex_Dpp = cmplx(Dpp,0d0,kind=8)
    deallocate(Bpp,Cpp,Dpp)
    
    call complex_ppRLR(TDA,nOO,nVV,complex_Bpp,complex_Cpp,complex_Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA(ispin))

    call print_excitation_energies('Real ppRPA@RHF','2p (triplet)',nVV,real(Om1))
    call print_excitation_energies('Imag ppRPA@RHF','2p (triplet)',nVV,aimag(Om1))
    call print_excitation_energies('Real ppRPA@RHF','2h (triplet)',nOO,real(Om2))
    call print_excitation_energies('Imag ppRPA@RHF','2h (triplet)',nOO,aimag(Om2))

    deallocate(Om1,X1,Y1,Om2,X2,Y2,complex_Bpp,complex_Cpp,complex_Dpp)

  end if

  EcRPA(2) = 3d0*EcRPA(2)

  write(*,*)
  write(*,*)'Real contribution'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF correlation energy (singlet) = ',real(EcRPA(1)),'au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF correlation energy (triplet) = ',real(EcRPA(2)),'au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF correlation energy           = ',real(sum(EcRPA)),'au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF total       energy           = ',real(ENuc + ERHF + sum(EcRPA)),'au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

  write(*,*)
  write(*,*)'Imaginary contribution'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF correlation energy (singlet) = ',aimag(EcRPA(1)),'au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF correlation energy (triplet) = ',aimag(EcRPA(2)),'au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF correlation energy           = ',aimag(sum(EcRPA)),'au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppRPA@RHF total       energy           = ',aimag(ENuc + ERHF + sum(EcRPA)),'au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)
end subroutine
