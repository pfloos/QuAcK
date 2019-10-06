subroutine ppRPA(singlet_manifold,triplet_manifold,nBas,nC,nO,nV,nR,ENuc,ERHF,ERI,e)

! Perform pp-RPA calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
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

  logical                       :: dRPA
  logical                       :: TDA
  logical                       :: BSE
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

! Useful quantities

 nOO = nO*(nO-1)/2
 nVV = nV*(nV-1)/2

! Initialization

  Ec_ppRPA(:) = 0d0

! Switch on exchange for TDHF

  dRPA = .false.
 
! Switch off Tamm-Dancoff approximation for TDHF

  TDA = .false.
 
! Switch off Bethe-Salpeter equation for TDHF

  BSE = .false. 

! Memory allocation

  allocate(Omega1(nVV,nspin),X1(nVV,nVV,nspin),Y1(nOO,nVV,nspin), & 
           Omega2(nOO,nspin),X2(nVV,nOO,nspin),Y2(nOO,nOO,nspin))

! Singlet manifold

  if(singlet_manifold) then 

    ispin = 1

    call linear_response_pp(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI, & 
                            Omega1(:,ispin),X1(:,:,ispin),Y1(:,:,ispin),       & 
                            Omega2(:,ispin),X2(:,:,ispin),Y2(:,:,ispin),       & 
                            Ec_ppRPA(ispin))
    call print_excitation('pp-RPA (N+2)',ispin,nVV,Omega1(:,ispin))
    call print_excitation('pp-RPA (N-2)',ispin,nOO,Omega2(:,ispin))

  endif

! Triplet manifold 

  if(triplet_manifold) then 

    ispin = 2

    call linear_response_pp(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nOO,nVV,e,ERI, &
                            Omega1(:,ispin),X1(:,:,ispin),Y1(:,:,ispin),       & 
                            Omega2(:,ispin),X2(:,:,ispin),Y2(:,:,ispin),       & 
                            Ec_ppRPA(ispin))
    call print_excitation('pp-RPA (N+2)',ispin,nVV,Omega1(:,ispin))
    call print_excitation('pp-RPA (N-2)',ispin,nOO,Omega2(:,ispin))

  endif

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A40,F15.6)') 'pp-RPA   correlation energy (singlet) =',Ec_ppRPA(1)
  write(*,'(2X,A40,F15.6)') 'pp-RPA   correlation energy (triplet) =',Ec_ppRPA(2)
  write(*,'(2X,A40,F15.6)') 'pp-RPA   correlation energy           =',Ec_ppRPA(1) + Ec_ppRPA(2)
  write(*,'(2X,A40,F15.6)') 'pp-RPA   total energy                 =',ENuc + ERHF + Ec_ppRPA(1) + Ec_ppRPA(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

end subroutine ppRPA