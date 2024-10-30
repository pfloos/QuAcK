subroutine GGTpp_ppBSE(TDA_T,TDA,dBSE,dTDA,eta,nOrb,nC,nO,nV,nR,nOO,nVV, &
                      ERI,dipole_int,eT,eGT,EcBSE)

! Compute the Bethe-Salpeter excitation energies with the T-matrix kernel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR

  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV

  double precision,intent(in)   :: eT(nOrb)
  double precision,intent(in)   :: eGT(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: Bpp(:,:),Cpp(:,:),Dpp(:,:)
  double precision,allocatable  :: Om1(:), Om2(:)
  double precision,allocatable  :: X1(:,:), X2(:,:)
  double precision,allocatable  :: Y1(:,:), Y2(:,:)
  double precision,allocatable  :: rho1(:,:,:), rho2(:,:,:)
  double precision,allocatable  :: KB_sta(:,:),KC_sta(:,:),KD_sta(:,:)
  double precision,allocatable  :: T(:,:,:,:)
  
! Output variables

  double precision,intent(out)  :: EcBSE

!----------------------------------------------
! Compute linear response
!----------------------------------------------
  allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO))
  allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))
  Bpp(:,:) = 0d0
  Cpp(:,:) = 0d0
  Dpp(:,:) = 0d0
  call ppGLR_C(nOrb,nC,nO,nV,nR,nVV,1d0,eT,ERI,Cpp)
  call ppGLR_D(nOrb,nC,nO,nV,nR,nOO,1d0,eT,ERI,Dpp)
  if(.not.TDA_T) call ppGLR_B(nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)

  call ppGLR(TDA_T,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA)

  deallocate(Bpp,Cpp,Dpp)

!----------------------------------------------
! Compute excitation densities
!----------------------------------------------
  allocate(rho1(nOrb,nOrb,nVV),rho2(nOrb,nOrb,nOO))
  call GGTpp_excitation_density(nOrb,nC,nO,nV,nR,nOO,nVV,ERI,X1,Y1,rho1,X2,Y2,rho2)

  deallocate(X1,Y1,X2,Y2)
  
!----------------------------------------------
! Compute new ppRPA block
!----------------------------------------------

  allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))
  
  call ppGLR_C(nOrb,nC,nO,nV,nR,nVV,1d0,eGT,ERI,Cpp)
  call ppGLR_D(nOrb,nC,nO,nV,nR,nOO,1d0,eGT,ERI,Dpp)
  if(.not.TDA_T) call ppGLR_B(nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)

!----------------------------------------------
! Compute T matrix tensor
!----------------------------------------------

  allocate(T(nOrb,nOrb,nOrb,nOrb))

  call GGT_Tmatrix(nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,eGT,Om1,rho1,Om2,rho2,T)  
  
!----------------------------------------------
! Compute kernels
!----------------------------------------------

  allocate(KB_sta(nVV,nOO),KC_sta(nVV,nVV),KD_sta(nOO,nOO))
  
  call GGTpp_ppBSE_static_kernel_C(eta,nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,eGT, & 
                                  Om1,rho1,Om2,rho2,T,KC_sta)
  call GGTpp_ppBSE_static_kernel_D(eta,nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,eGT, & 
                                  Om1,rho1,Om2,rho2,T,KD_sta)
  if(.not.TDA_T) call GGTpp_ppBSE_static_kernel_B(eta,nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,eGT, & 
                                                  Om1,rho1,Om2,rho2,T,KB_sta)

  deallocate(Om1,Om2,rho1,rho2)
! Deallocate the 4-tensor T 
  deallocate(T)

!----------------------------------------------
! Diagonalize ppBSE
!----------------------------------------------

  Bpp(:,:) = Bpp(:,:) + KB_sta(:,:) 
  Cpp(:,:) = Cpp(:,:) + KC_sta(:,:) 
  Dpp(:,:) = Dpp(:,:) + KD_sta(:,:) 

  allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO))
  
  call ppGLR(TDA_T,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcBSE)

  call ppLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

  !----------------------------------------------------!
  ! Compute the dynamical screening at the ppBSE level !
  !----------------------------------------------------!

  ! TODO
  
  deallocate(Om1,X1,Y1,Om2,X2,Y2,Bpp,Cpp,Dpp,KB_sta,KC_sta,KD_sta) 

end subroutine
