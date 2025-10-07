subroutine RGTpp_phBSE(exchange_kernel,TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                       Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,Om2t,X2t,Y2t,rho1t,rho2t,ERI,dipole_int,eT,eGT,EcBSE)

! Compute the Bethe-Salpeter excitation energies with the T-matrix kernel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  integer,intent(in)            :: nOOs
  integer,intent(in)            :: nOOt
  integer,intent(in)            :: nVVs
  integer,intent(in)            :: nVVt

  double precision,intent(in)   :: eT(nOrb)
  double precision,intent(in)   :: eGT(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

  double precision,intent(in)   :: Om1s(nVVs)
  double precision,intent(in)   :: X1s(nVVs,nVVs)
  double precision,intent(in)   :: Y1s(nOOs,nVVs)
  double precision,intent(in)   :: Om2s(nOOs)
  double precision,intent(in)   :: X2s(nVVs,nOOs)
  double precision,intent(in)   :: Y2s(nOOs,nOOs)
  double precision,intent(in)   :: rho1s(nOrb,nOrb,nVVs)
  double precision,intent(in)   :: rho2s(nOrb,nOrb,nOOs)
  double precision,intent(in)   :: Om1t(nVVt)
  double precision,intent(in)   :: X1t(nVVt,nVVt)
  double precision,intent(in)   :: Y1t(nOOt,nVVt)
  double precision,intent(in)   :: Om2t(nOOt)
  double precision,intent(in)   :: X2t(nVVt,nOOt)
  double precision,intent(in)   :: Y2t(nOOt,nOOt) 
  double precision,intent(in)   :: rho1t(nOrb,nOrb,nVVt)
  double precision,intent(in)   :: rho2t(nOrb,nOrb,nOOt)

! Local variables

  logical                       :: dRPA = .false.

  integer                       :: ispin

  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: KA_sta(:,:),KB_sta(:,:)
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY_BSE(:,:)
  double precision,allocatable  :: XmY_BSE(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),KA_sta(nS,nS),KB_sta(nS,nS), & 
           OmBSE(nS),XpY_BSE(nS,nS),XmY_BSE(nS,nS))

!-----!
! TDA !
!-----!

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated in phBSE!'
    write(*,*)
  end if

!------------------------------------!
! Compute T-matrix for singlet block !
!------------------------------------!

  ispin  = 1

  allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

  if(.not.TDA_T) call ppRLR_B(ispin,nOrb,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                 call ppRLR_C(ispin,nOrb,nC,nO,nV,nR,nVVs,1d0,eT,ERI,Cpp)
                 call ppRLR_D(ispin,nOrb,nC,nO,nV,nR,nOOs,1d0,eT,ERI,Dpp)

  call ppRLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

!------------------------------------!
! Compute T-matrix for triplet block !
!------------------------------------!

  ispin  = 2

  allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))

  if(.not.TDA_T) call ppRLR_B(ispin,nOrb,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                 call ppRLR_C(ispin,nOrb,nC,nO,nV,nR,nVVt,1d0,eT,ERI,Cpp)
                 call ppRLR_D(ispin,nOrb,nC,nO,nV,nR,nOOt,1d0,eT,ERI,Dpp)

  call ppRLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

!----------------------------------------------
! Compute excitation densities
!----------------------------------------------

  ispin = 1

  call RGTpp_excitation_density(ispin,nOrb,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

  ispin = 2

  call RGTpp_excitation_density(ispin,nOrb,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)
  
!------------------!
! Singlet manifold !
!------------------!

 if(singlet) then

    ispin = 1

    ! Compute BSE singlet excitation energies

                 call phRLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eGT,ERI,Aph)
    if(.not.TDA) call phRLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

                 call RGTpp_phBSE_static_kernel_A(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,1d0,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,KA_sta)
    if(.not.TDA) call RGTpp_phBSE_static_kernel_B(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,1d0,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,KB_sta)

                 Aph(:,:) = Aph(:,:) + KA_sta(:,:) 
    if(.not.TDA) Bph(:,:) = Bph(:,:) + KB_sta(:,:)

    call phRLR(TDA,nS,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GTpp','singlet',nS,OmBSE)
    call phLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

    ! TODO This is old code and should be properly spin adapted now
    ! Compute dynamic correction for BSE via renormalized perturbation theory 
    ! if(dBSE) call RGTpp_phBSE_dynamic_perturbation(ispin,dTDA,eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
    !                                                Om1s,Om2s,Om1t,Om2t,rho1s,rho2s,rho1t,rho2t,eT,eGT, & 
    !                                                dipole_int,OmBSE,XpY_BSE,XmY_BSE,TAs,TAt)
    
  end if

!------------------!
! Triplet manifold !
!------------------!

  if(triplet) then

    ispin = 2

    ! Compute BSE triplet excitation energies

                 call phRLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eGT,ERI,Aph)
    if(.not.TDA) call phRLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

                 call RGTpp_phBSE_static_kernel_A(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,1d0,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,KA_sta)
    if(.not.TDA) call RGTpp_phBSE_static_kernel_B(ispin,eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,1d0,Om1s,rho1s,Om2s,rho2s,Om1t,rho1t,Om2t,rho2t,KB_sta)
    
                 Aph(:,:) = Aph(:,:) + KA_sta(:,:)
    if(.not.TDA) Bph(:,:) = Bph(:,:) + KB_sta(:,:)

    call phRLR(TDA,nS,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GTpp','triplet',nS,OmBSE)
    call phLR_transition_vectors(.false.,nOrb,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

    ! TODO This is old code and should be properly spin adapted now
    ! Compute dynamic correction for BSE via renormalized perturbation theory 
    ! if(dBSE) call RGTpp_phBSE_dynamic_perturbation(ispin,dTDA,eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
    !                                                Om1s,Om2s,Om1t,Om2t,rho1s,rho2s,rho1t,rho2t,eT,eGT, & 
    !                                                dipole_int,OmBSE,XpY_BSE,XmY_BSE,TAs,TAt)

  end if

  if(exchange_kernel) then

    EcBSE(1) = 0.5d0*EcBSE(1)
    EcBSE(2) = 1.5d0*EcBSE(1)

  end if

end subroutine 
