subroutine GGTpp_ehBSE(TDA_T,TDA,dBSE,dTDA,eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,ERI,dipole_int,eT,eGT,EcBSE)

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
  integer,intent(in)            :: nS

  integer,intent(in)            :: nOO
  integer,intent(in)            :: nVV

  double precision,intent(in)   :: eT(nOrb)
  double precision,intent(in)   :: eGT(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  double precision              :: EcRPA
  double precision,allocatable  :: Bpp(:,:),Cpp(:,:),Dpp(:,:)
  double precision,allocatable  :: Aph(:,:),Bph(:,:)
  double precision,allocatable  :: Om1(:), Om2(:), Om(:)
  double precision,allocatable  :: X1(:,:), X2(:,:)
  double precision,allocatable  :: Y1(:,:), Y2(:,:)
  double precision,allocatable  :: XpY(:,:), XmY(:,:)
  double precision,allocatable  :: rho1(:,:,:), rho2(:,:,:)
  double precision,allocatable  :: KA_sta(:,:),KB_sta(:,:)
  
! Output variables

  double precision,intent(out)  :: EcBSE

!----------------------------------------------
! Compute linear response
!----------------------------------------------

  allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO))
  allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))

  if(.not.TDA_T) call ppGLR_B(nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)
                 call ppGLR_C(nOrb,nC,nO,nV,nR,nVV,1d0,eT,ERI,Cpp)
                 call ppGLR_D(nOrb,nC,nO,nV,nR,nOO,1d0,eT,ERI,Dpp)

  call ppGLR(TDA_T,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA)

  deallocate(Bpp,Cpp,Dpp)

!----------------------------------------------
! Compute excitation densities
!----------------------------------------------

  allocate(rho1(nOrb,nOrb,nVV),rho2(nOrb,nOrb,nOO))

  call GGTpp_excitation_density(nOrb,nC,nO,nV,nR,nOO,nVV,ERI,X1,Y1,rho1,X2,Y2,rho2)

  deallocate(X1,Y1,X2,Y2)
  
!----------------------------------------------
! Compute ehRPA block
!----------------------------------------------

  allocate(Aph(nS,nS),Bph(nS,nS))

               call phGLR_A(.false.,nOrb,nC,nO,nV,nR,nS,1d0,eGT,ERI,Aph)
  if(.not.TDA) call phGLR_B(.false.,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)
  
!----------------------------------------------
! Compute kernels
!----------------------------------------------

  allocate(KA_sta(nS,nS),KB_sta(nS,nS))
  
               call GGTpp_phBSE_static_kernel_A(eta,nOrb,nC,nO,nV,nR,nOO,nVV,nS,1d0,ERI,eGT,Om1,rho1,Om2,rho2,KA_sta)
  if(.not.TDA) call GGTpp_phBSE_static_kernel_B(eta,nOrb,nC,nO,nV,nR,nOO,nVV,nS,1d0,ERI,eGT,Om1,rho1,Om2,rho2,KB_sta)

  deallocate(Om1,Om2,rho1,rho2)

!----------------------------------------------
! Diagonalize ehBSE
!----------------------------------------------
 
  Aph(:,:) = Aph(:,:) + KA_sta(:,:)
  Bph(:,:) = Bph(:,:) + KB_sta(:,:)

  allocate(Om(nS),XpY(nS,nS),XmY(nS,nS))

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated in phBSE!'
    write(*,*)
  end if
  
  call phGLR(TDA,nS,Aph,Bph,EcBSE,Om,XpY,XmY)

  call phLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nS,dipole_int,Om,XpY,XmY)
  
  deallocate(Om,XpY,XmY,Aph,Bph,KA_sta,KB_sta) 

end subroutine
