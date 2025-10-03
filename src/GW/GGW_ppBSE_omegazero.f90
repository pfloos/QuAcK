subroutine GGW_ppBSE_omegazero(TDA_W,TDA,dBSE,dTDA,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

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

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nOrb,nOrb,nS), &
           Aph(nS,nS),Bph(nS,nS))
 
                 call phGLR_A(dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,eW,ERI,Aph)
  if(.not.TDA_W) call phGLR_B(dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phGLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
!  call phLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nS,dipole_int,OmRPA,XpY_RPA,XmY_RPA)
  
  call GGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  deallocate(XpY_RPA,XmY_RPA,Aph,Bph)

  nOO = nO*(nO-1)/2
  nVV = nV*(nV-1)/2

  allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),       &
           Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO),       &
           Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO), &
           KB_sta(nVV,nOO),KC_sta(nVV,nVV),KD_sta(nOO,nOO))

  KB_sta(:,:) = 0d0
  KC_sta(:,:) = 0d0
  
  ! Compute BSE excitation energies

!               call GGW_ppBSE_static_kernel_C(eta,nOrb,nC,nO,nV,nR,nS,nVV,1d0,ERI,OmRPA,rho_RPA,KC_sta)
               call GGW_ppBSE_static_kernel_D_omegazero(eta,nOrb,nC,nO,nV,nR,nS,nOO,1d0,ERI,eGW,OmRPA,rho_RPA,KD_sta)
!  if(.not.TDA) call GGW_ppBSE_static_kernel_B(eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV,1d0,ERI,OmRPA,rho_RPA,KB_sta)
  
               call ppGLR_C(nOrb,nC,nO,nV,nR,nVV,1d0,eGW,ERI,Cpp)
               call ppGLR_D(nOrb,nC,nO,nV,nR,nOO,1d0,eGW,ERI,Dpp)
  if(.not.TDA) call ppGLR_B(nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)

  Bpp(:,:) = Bpp(:,:) + KB_sta(:,:)
  Cpp(:,:) = Cpp(:,:) + KC_sta(:,:)
  Dpp(:,:) = Dpp(:,:) + KD_sta(:,:)

  call ppGLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcBSE)

  call ppLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nOO,nVV,dipole_int,Om1,X1,Y1,Om2,X2,Y2)

  deallocate(Om1,X1,Y1,Om2,X2,Y2,Bpp,Cpp,Dpp,KB_sta,KC_sta,KD_sta)

end subroutine GGW_ppBSE_omegazero

subroutine GGW_ppBSE_static_kernel_D_omegazero(eta,nBas,nC,nO,nV,nR,nS,nOO,lambda,ERI,eGW,Om,rho,KD)

! Compute the OOOO block of the static screening W for the pp-BSE

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  integer,intent(in)            :: nOO
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: lambda
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: Om(nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  double precision,external     :: Kronecker_delta
  double precision              :: num
  double precision              :: dem
  integer                       :: i,j,k,l,ij,kl,m

! Output variables

  double precision,intent(out)  :: KD(nOO,nOO)

! Initialization

  KD(:,:) = 0d0

!---------------!
! SpinOrb block !
!---------------!

  ij = 0
  do i=nC+1,nO
    do j=i+1,nO
      ij = ij + 1
      kl = 0
      do k=nC+1,nO
        do l=k+1,nO
          kl = kl + 1
  
          do m=1,nS

               num = (rho(i,k,m)*rho(j,l,m) - rho(j,k,m)*rho(i,l,m))/2d0
               dem = - Om(m) + eGW(j) + eGW(l)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)

               dem = - Om(m) + eGW(i) + eGW(k)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)

               dem = - Om(m) + eGW(i) + eGW(l)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)
       
               dem = - Om(m) + eGW(j) + eGW(k)
               KD(ij,kl) = KD(ij,kl) + num*dem/(dem**2 + eta**2)

          end do

        end do
      end do
    end do
  end do

end subroutine GGW_ppBSE_static_kernel_D_omegazero
