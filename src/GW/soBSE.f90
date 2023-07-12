subroutine soBSE(TDA_W,TDA,eta,nBas,nC,nO,nV,nR,nS,ERI,eW,eGW)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin
  integer                       :: isp_W

  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY_BSE(:,:)
  double precision,allocatable  :: XmY_BSE(:,:)

  double precision,allocatable  :: WA_sta(:,:)
  double precision,allocatable  :: WB_sta(:,:)

  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: Y(:,:)
  double precision,allocatable  :: Xinv(:,:)
  double precision,allocatable  :: t(:,:,:,:)

  integer                       ::i,a,j,b,k,c,ia,jb,kc

  double precision              :: EcBSE

! Memory allocation

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nBas,nBas,nS), &
           WA_sta(nS,nS),WB_sta(nS,nS),OmBSE(nS),XpY_BSE(nS,nS),XmY_BSE(nS,nS))

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  isp_W = 3
  EcRPA = 0d0

  call phLR(isp_W,.true.,TDA_W,eta,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call GW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  call static_screening_WA_so(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,WA_sta)
  call static_screening_WB_so(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,WB_sta)


  ispin = 3
  EcBSE = 0d0

  ! Compute BSE excitation energies

  call linear_response_BSE(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI,WA_sta,WB_sta, &
                           EcBSE,OmBSE,XpY_BSE,XmY_BSE)
  call print_excitation('soBSE@GW    ',ispin,nS,OmBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0 correlation energy =',0.5d0*EcBSE,' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

! allocate(X(nS,nS),Y(nS,nS),Xinv(nS,nS),t(nO,nO,nV,nV))

! X(:,:) = transpose(0.5d0*(XpY_BSE(:,:) + XmY_BSE(:,:)))
! Y(:,:) = transpose(0.5d0*(XpY_BSE(:,:) - XmY_BSE(:,:)))

! call matout(nS,nS,matmul(X,transpose(X))-matmul(Y,transpose(Y)))

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

end subroutine 
