subroutine RGTpp_ppBSE(TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt, &
                      ERI,dipole_int,eT,eGT,EcBSE)

! Compute the Bethe-Salpeter excitation energies with the T-matrix kernel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR

  integer,intent(in)            :: nOOs
  integer,intent(in)            :: nOOt
  integer,intent(in)            :: nVVs
  integer,intent(in)            :: nVVt

  double precision,intent(in)   :: eT(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin, isp_T

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: Bpp(:,:),Cpp(:,:),Dpp(:,:)
  double precision,allocatable  :: Om1s(:),Om1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Om2s(:),Om2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)
  double precision,allocatable  :: KB_sta(:,:)
  double precision,allocatable  :: KC_sta(:,:)
  double precision,allocatable  :: KD_sta(:,:)
  double precision,allocatable  :: Taaaa(:,:,:,:),Tabab(:,:,:,:),Tbaab(:,:,:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

!---------------------------------
! Compute ppRPA excitation density 
!---------------------------------

  ! Singlet contribution

  isp_T  = 1

  allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

  if(.not.TDA_T) call ppLR_B(isp_T,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                 call ppLR_C(isp_T,nBas,nC,nO,nV,nR,nVVs,1d0,eT,ERI,Cpp)
                 call ppLR_D(isp_T,nBas,nC,nO,nV,nR,nOOs,1d0,eT,ERI,Dpp)

  allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs))
  allocate(Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs))

  call ppLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(isp_T))
!  call ppLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nOOs,nVVs,dipole_int,Om1s,X1s,Y1s,Om2s,X2s,Y2s)
  
  allocate(rho1s(nBas,nBas,nVVs),rho2s(nBas,nBas,nOOs))

  call RGTpp_excitation_density(isp_T,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

  deallocate(X1s,Y1s,X2s,Y2s,Bpp,Cpp,Dpp) 

  ! Triplet contribution

  isp_T  = 2

  allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))

  if(.not.TDA_T) call ppLR_B(isp_T,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                 call ppLR_C(isp_T,nBas,nC,nO,nV,nR,nVVt,1d0,eT,ERI,Cpp)
                 call ppLR_D(isp_T,nBas,nC,nO,nV,nR,nOOt,1d0,eT,ERI,Dpp)

  allocate(Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt))
  allocate(Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt))

  call ppLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(isp_T))
!  call ppLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nOOt,nVVt,dipole_int,Om1t,X1t,Y1t,Om2t,X2t,Y2t)

  allocate(rho1t(nBas,nBas,nVVt),rho2t(nBas,nBas,nOOt))

  call RGTpp_excitation_density(isp_T,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

  deallocate(X1t,Y1t,X2t,Y2t,Bpp,Cpp,Dpp) 
  
!---------------------------------
! Compute T matrix elements
!---------------------------------
  
  ! Elements aaaa

  isp_T = 1
  allocate(Taaaa(nBas,nBas,nBas,nBas))

  call RGT_Tmatrix(isp_T,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,1d0,ERI,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,Taaaa)

  ! Elements abab

  isp_T = 2
  allocate(Tabab(nBas,nBas,nBas,nBas))

  call RGT_Tmatrix(isp_T,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,1d0,ERI,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,Tabab)
  
  ! Elements baab

  isp_T = 3
  allocate(Tbaab(nBas,nBas,nBas,nBas))

  call RGT_Tmatrix(isp_T,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,1d0,ERI,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,Tbaab)

  deallocate(Om1s,Om2s,Om1t,Om2t,rho1s,rho2s,rho1t,rho2t) 
  
!------------------!
! Singlet manifold !
!------------------!

  if(singlet) then

    ispin = 1
    
    allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs))
    allocate(Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs))

    ! Compute BSE excitation energies

    allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))
    allocate(KB_sta(nVVs,nOOs),KC_sta(nVVs,nVVs),KD_sta(nOOs,nOOs))

    if(.not.TDA) call ppLR_B(ispin,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                 call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVVs,1d0,eGT,ERI,Cpp)
                 call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOOs,1d0,eGT,ERI,Dpp)

    if(.not.TDA) call RGTpp_ppBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,eGT,Taaaa,Tabab,Tbaab,KB_sta)
                 call RGTpp_ppBSE_static_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,eGT,Taaaa,Tabab,Tbaab,KC_sta)
                 call RGTpp_ppBSE_static_kernel_D(ispin,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,eGT,Taaaa,Tabab,Tbaab,KD_sta)
                 
    Bpp(:,:) = Bpp(:,:) + KB_sta(:,:)
    Cpp(:,:) = Cpp(:,:) + KC_sta(:,:)
    Dpp(:,:) = Dpp(:,:) + KD_sta(:,:)

    call ppLR(TDA,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcBSE(ispin))

    call ppLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nOOs,nVVs,dipole_int,Om1s,X1s,Y1s,Om2s,X2s,Y2s)

    deallocate(Om1s,X1s,Y1s,Om2s,X2s,Y2s,Bpp,Cpp,Dpp,KB_sta,KC_sta,KD_sta) 

  end if

!------------------!
! Triplet manifold !
!------------------!

  if(triplet) then

    ispin  = 2
    
    EcBSE(ispin) = 0d0

    allocate(Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt))
    allocate(Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt))

    ! Compute BSE excitation energies

    allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))
    allocate(KB_sta(nVVt,nOOt),KC_sta(nVVt,nVVt),KD_sta(nOOt,nOOt))
  
    if(.not.TDA) call ppLR_B(ispin,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                 call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVVt,1d0,eGT,ERI,Cpp)
                 call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOOt,1d0,eGT,ERI,Dpp)

    if(.not.TDA) call RGTpp_ppBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,eGT,Taaaa,Tabab,Tbaab,KB_sta)
                 call RGTpp_ppBSE_static_kernel_C(ispin,eta,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,eGT,Taaaa,Tabab,Tbaab,KC_sta)
                 call RGTpp_ppBSE_static_kernel_D(ispin,eta,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,eGT,Taaaa,Tabab,Tbaab,KD_sta)
                 
    Bpp(:,:) = Bpp(:,:) + KB_sta(:,:)
    Cpp(:,:) = Cpp(:,:) + KC_sta(:,:)
    Dpp(:,:) = Dpp(:,:) + KD_sta(:,:)

    call ppLR(TDA,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcBSE(ispin))

    call ppLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nOOt,nVVt,dipole_int,Om1t,X1t,Y1t,Om2t,X2t,Y2t)

    deallocate(Om1t,X1t,Y1t,Om2t,X2t,Y2t,Bpp,Cpp,Dpp)

  end if

end subroutine 
