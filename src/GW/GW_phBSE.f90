subroutine GW_phBSE(dophBSE2,TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: TDA_W
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
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: dRPA   = .false.
  logical                       :: dRPA_W = .true.

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

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)

  double precision,allocatable  :: KA_sta(:,:)
  double precision,allocatable  :: KB_sta(:,:)

  double precision,allocatable  :: W(:,:,:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nBas,nBas,nS), &
           Aph(nS,nS),Bph(nS,nS),KA_sta(nS,nS),KB_sta(nS,nS), &
           OmBSE(nS),XpY_BSE(nS,nS),XmY_BSE(nS,nS))

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  isp_W = 1
  EcRPA = 0d0

                 call phLR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI,Aph)
  if(.not.TDA_W) call phLR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call GW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  call GW_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KA_sta)
  call GW_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KB_sta)

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    ispin = 1
    EcBSE(ispin) = 0d0

    ! Compute BSE excitation energies

                 call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

    ! Second-order BSE static kernel
  
    if(dophBSE2) then 

      allocate(W(nBas,nBas,nBas,nBas))

      write(*,*) 
      write(*,*) '*** Second-order BSE static kernel activated! ***'
      write(*,*) 

      call GW_phBSE_static_kernel(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,W)
      call GW_phBSE2_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,1d0,eW,W,KA_sta)

      if(.not.TDA) call GW_phBSE2_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,1d0,eW,W,KB_sta)

      deallocate(W)

    end if

                 Aph(:,:) = Aph(:,:) + KA_sta(:,:)
    if(.not.TDA) Bph(:,:) = Bph(:,:) + KB_sta(:,:)

    call phLR(TDA,nS,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GW@RHF','singlet',nS,OmBSE)
    call phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

    !--------------------!
    ! Cumulant expansion !
    !--------------------!

    call RGWC(.false.,nBas,nC,nO,nR,nS,OmBSE,rho_RPA,eGW)

    !----------------------------------------------------!
    ! Compute the dynamical screening at the phBSE level !
    !----------------------------------------------------!

    if(dBSE) &
        call GW_phBSE_dynamic_perturbation(dophBSE2,dTDA,eta,nBas,nC,nO,nV,nR,nS,eW,eGW,ERI,dipole_int,OmRPA,rho_RPA, &
                                           OmBSE,XpY_BSE,XmY_BSE,KA_sta,KB_sta)

  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    ispin = 2
    EcBSE(ispin) = 0d0

    ! Compute BSE excitation energies

                 call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

                 Aph(:,:) = Aph(:,:) + KA_sta(:,:)
    if(.not.TDA) Bph(:,:) = Bph(:,:) + KB_sta(:,:)

    call phLR(TDA,nS,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GW@RHF','triplet',nS,OmBSE)
    call phLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) &
        call GW_phBSE_dynamic_perturbation(dophBSE2,dTDA,eta,nBas,nC,nO,nV,nR,nS,eW,eGW,ERI,dipole_int,OmRPA,rho_RPA, &
                                           OmBSE,XpY_BSE,XmY_BSE,KA_sta,KB_sta)

  end if

end subroutine 
