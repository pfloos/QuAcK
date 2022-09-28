subroutine soRPAx(TDA,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Perform random phase approximation calculation with exchange (aka TDHF) in the
! spinorbital basis

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: TDA
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin
  double precision,allocatable  :: Omega(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: Y(:,:)
  double precision,allocatable  :: Xinv(:,:)
  double precision,allocatable  :: t(:,:,:,:)

  double precision              :: EcRPAx

  integer                       ::i,a,j,b,k,c,ia,jb,kc

! Hello world

  write(*,*)
  write(*,*)'***********************************************************'
  write(*,*)'|  Random phase approximation calculation with exchange   |'
  write(*,*)'***********************************************************'
  write(*,*)

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*) ' => RPAx + TDA = CIS '
    write(*,*)
  end if

! Initialization

  EcRPAx = 0d0

! Memory allocation

  allocate(Omega(nS),XpY(nS,nS),XmY(nS,nS))

  ispin = 3

  call linear_response(ispin,.false.,TDA,0d0,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI,EcRPAx,Omega,XpY,XmY)
  call print_excitation('soRPAx@HF   ',ispin,nS,Omega)

  EcRPAx = 0.5d0*EcRPAx

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@RPAx correlation energy           =',EcRPAx
  write(*,'(2X,A50,F20.10)') 'Tr@RPAx total energy                 =',ENuc + ERHF + EcRPAx
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! allocate(X(nS,nS),Y(nS,nS),Xinv(nS,nS),t(nO,nO,nV,nV))

! X(:,:) = transpose(0.5d0*(XpY(:,:) + XmY(:,:)))
! Y(:,:) = transpose(0.5d0*(XpY(:,:) - XmY(:,:)))

! call matout(nS,nS,matmul(transpose(X),X)-matmul(transpose(Y),Y))

! call inverse_matrix(nS,X,Xinv)

! t = 0d0
! ia = 0
! do i=1,nO
!   do a=1,nV
!     ia = ia + 1

!     jb = 0
!     do j=1,nO
!       do b=1,nV
!       jb = jb + 1

!       kc = 0
!       do k=1,nO
!         do c=1,nV
!         kc = kc + 1

!           t(i,j,a,b) = t(i,j,a,b) + Y(ia,kc)*Xinv(kc,jb)

!           end do
!         end do
!       end do
!     end do
!   end do
! end do

! call matout(nO*nO,nV*nV,t)

end subroutine soRPAx
