subroutine ppGRPA(dotest,TDA,nBas,nC,nO,nV,nR,ENuc,EGHF,ERI,dipole_int,eHF)

! Perform ppGRPA calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

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

  double precision              :: EcRPA

! Hello world

  write(*,*)
  write(*,*)'**********************************'
  write(*,*)'* Generalized pp-RPA Calculation *'
  write(*,*)'**********************************'
  write(*,*)

! Initialization

  EcRPA = 0d0

  nOO = nO*(nO-1)/2
  nVV = nV*(nV-1)/2

  allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO), &
           Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))

  if(.not.TDA) call ppGLR_B(nBas,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)
               call ppGLR_C(nBas,nC,nO,nV,nR,nVV,1d0,eHF,ERI,Cpp)
               call ppGLR_D(nBas,nC,nO,nV,nR,nOO,1d0,eHF,ERI,Dpp)

  call ppGLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA)

!   call print_transition_vectors_pp(.true.,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

  call print_excitation_energies('ppRPA@GHF','2p (spinorbital)',nVV,Om1)
  call print_excitation_energies('ppRPA@GHF','2h (spinorbital)',nOO,Om2)

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppGRPA correlation energy           = ',EcRPA,' au'
  write(*,'(2X,A50,F20.10,A3)') 'Tr@ppGRPA total energy                 = ',ENuc + EGHF + EcRPA,' au'
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

  if(dotest) then
  
    call dump_test_value('G','ppRPA correlation energy',EcRPA)

  end if

end subroutine
