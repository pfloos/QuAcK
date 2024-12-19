subroutine GGF2_ppBSE(TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,ERI,dipole_int,eGF,EcBSE)

! Compute the Bethe-Salpeter excitation energies at the pp level

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  double precision,intent(in)   :: eGF(nBas)
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

  double precision,allocatable  :: KB_sta(:,:)
  double precision,allocatable  :: KC_sta(:,:)
  double precision,allocatable  :: KD_sta(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE

! Initialization

  EcBSE = 0d0

  nOO = nO*(nO-1)/2
  nVV = nV*(nV-1)/2

  allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),       &
           Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO),       &
           Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO), &
           KB_sta(nVV,nOO),KC_sta(nVV,nVV),KD_sta(nOO,nOO))

  ! Compute BSE excitation energies

  if(.not.TDA) call GGF2_ppBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nOO,nVV,1d0,ERI,eGF,KB_sta)
               call GGF2_ppBSE_static_kernel_C(eta,nBas,nC,nO,nV,nR,nVV,1d0,ERI,eGF,KC_sta)
               call GGF2_ppBSE_static_kernel_D(eta,nBas,nC,nO,nV,nR,nOO,1d0,ERI,eGF,KD_sta)

  if(.not.TDA) call ppGLR_B(nBas,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)
               call ppGLR_C(nBas,nC,nO,nV,nR,nVV,1d0,eGF,ERI,Cpp)
               call ppGLR_D(nBas,nC,nO,nV,nR,nOO,1d0,eGF,ERI,Dpp)

  Bpp(:,:) = Bpp(:,:) + KB_sta(:,:)
  Cpp(:,:) = Cpp(:,:) + KC_sta(:,:)
  Dpp(:,:) = Dpp(:,:) + KD_sta(:,:)

  call ppGLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcBSE)

  call ppLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

  !----------------------------------------------------!
  ! Compute the dynamical screening at the ppBSE level !
  !----------------------------------------------------!

! if(dBSE) &
!     call GGF2_ppBSE_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nOO,nVV,eGF,ERI,dipole_int, &
!                                          Om1,X1,Y1,Om2,X2,Y2,KB_sta,KC_sta,KD_sta)

  deallocate(Om1,X1,Y1,Om2,X2,Y2,Bpp,Cpp,Dpp,KB_sta,KC_sta,KD_sta)

end subroutine 
