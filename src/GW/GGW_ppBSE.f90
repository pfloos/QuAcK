subroutine GGW_ppBSE(TDA_W,TDA,dBSE,dTDA,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies at the pp level based on a GHF reference

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eW(nOrb)
  double precision,intent(in)   :: eGW(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  logical                       :: dRPA   = .false.
  logical                       :: dRPA_W = .true.

  integer                       :: nOO
  integer                       :: nVV

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)

  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

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

!-----------------------!
! Compute RPA screening !
!-----------------------!

  EcRPA = 0d0

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nOrb,nOrb,nS), &
           Aph(nS,nS),Bph(nS,nS))
 
                 call phGLR_A(dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,eW,ERI,Aph)
  if(.not.TDA_W) call phGLR_B(dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  deallocate(XpY_RPA,XmY_RPA,Aph,Bph)

  EcBSE = 0d0

  nOO = nO*(nO-1)/2
  nVV = nV*(nV-1)/2

  allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),       &
           Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO),       &
           Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO), &
           KB_sta(nVV,nOO),KC_sta(nVV,nVV),KD_sta(nOO,nOO))

  ! Compute BSE excitation energies

               call GGW_ppBSE_static_kernel_C(eta,nOrb,nC,nO,nV,nR,nS,nVV,1d0,ERI,OmRPA,rho_RPA,KC_sta)
               call GGW_ppBSE_static_kernel_D(eta,nOrb,nC,nO,nV,nR,nS,nOO,1d0,ERI,OmRPA,rho_RPA,KD_sta)
  if(.not.TDA) call GGW_ppBSE_static_kernel_B(eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,1d0,ERI,OmRPA,rho_RPA,KB_sta)
  
               call ppGLR_C(nOrb,nC,nO,nV,nR,nVV,1d0,eGW,ERI,Cpp)
               call ppGLR_D(nOrb,nC,nO,nV,nR,nOO,1d0,eGW,ERI,Dpp)
  if(.not.TDA) call ppGLR_B(nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)

  Bpp(:,:) = Bpp(:,:) + KB_sta(:,:)
  Cpp(:,:) = Cpp(:,:) + KC_sta(:,:)
  Dpp(:,:) = Dpp(:,:) + KD_sta(:,:)

  call ppLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcBSE)

  call ppLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

  !----------------------------------------------------!
  ! Compute the dynamical screening at the ppBSE level !
  !----------------------------------------------------!

! if(dBSE) &
!     call GGW_ppBSE_dynamic_perturbation(dTDA,eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,eW,eGW,ERI,dipole_int,OmRPA,rho_RPA, &
!                                        Om1,X1,Y1,Om2,X2,Y2,KB_sta,KC_sta,KD_sta)


  deallocate(Om1,X1,Y1,Om2,X2,Y2,Bpp,Cpp,Dpp,KB_sta,KC_sta,KD_sta)

end subroutine 
