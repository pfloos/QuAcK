subroutine ppRPA(TDA,singlet,triplet,nBas,nC,nO,nV,nR,ENuc,ERHF,ERI,e)

! Perform pp-RPA calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
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

! Local variables

  integer                       :: ispin
  integer                       :: nOO
  integer                       :: nVV
  double precision,allocatable  :: Omega1(:,:)
  double precision,allocatable  :: X1(:,:,:)
  double precision,allocatable  :: Y1(:,:,:)
  double precision,allocatable  :: Omega2(:,:)
  double precision,allocatable  :: X2(:,:,:)
  double precision,allocatable  :: Y2(:,:,:)

  double precision              :: Ec_ppRPA(nspin)

! Hello world

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'|  particle-particle RPA calculation   |'
  write(*,*)'****************************************'
  write(*,*)

! Initialization

  Ec_ppRPA(:) = 0d0

! Singlet manifold

  if(singlet) then 

    ispin = 1

  ! Useful quantities

    nOO = nO*(nO+1)/2
    nVV = nV*(nV+1)/2

   ! Memory allocation

    allocate(Omega1(nVV,nspin),X1(nVV,nVV,nspin),Y1(nOO,nVV,nspin), & 
             Omega2(nOO,nspin),X2(nVV,nOO,nspin),Y2(nOO,nOO,nspin))

    call linear_response_pp(ispin,TDA,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,    & 
                            Omega1(:,ispin),X1(:,:,ispin),Y1(:,:,ispin), & 
                            Omega2(:,ispin),X2(:,:,ispin),Y2(:,:,ispin), & 
                            Ec_ppRPA(ispin))

    call print_excitation('pp-RPA (N+2)',ispin,nVV,Omega1(:,ispin))
    call print_excitation('pp-RPA (N-2)',ispin,nOO,Omega2(:,ispin))

    deallocate(Omega1,X1,Y1,Omega2,X2,Y2)

  endif

! Triplet manifold 

  if(triplet) then 

    ispin = 2

  ! Useful quantities

    nOO = nO*(nO-1)/2
    nVV = nV*(nV-1)/2

  ! Memory allocation

    allocate(Omega1(nVV,nspin),X1(nVV,nVV,nspin),Y1(nOO,nVV,nspin), & 
             Omega2(nOO,nspin),X2(nVV,nOO,nspin),Y2(nOO,nOO,nspin))


    call linear_response_pp(ispin,TDA,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI,    &
                            Omega1(:,ispin),X1(:,:,ispin),Y1(:,:,ispin), & 
                            Omega2(:,ispin),X2(:,:,ispin),Y2(:,:,ispin), & 
                            Ec_ppRPA(ispin))

    call print_excitation('pp-RPA (N+2)',ispin,nVV,Omega1(:,ispin))
    call print_excitation('pp-RPA (N-2)',ispin,nOO,Omega2(:,ispin))

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

end subroutine ppRPA
