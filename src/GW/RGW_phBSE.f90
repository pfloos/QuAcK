subroutine RGW_phBSE(dophBSE2,exchange_kernel,TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta, & 
                     nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: TDA_W
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
  double precision,intent(in)   :: eW(nOrb)
  double precision,intent(in)   :: eGW(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

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

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nOrb,nOrb,nS), &
           Aph(nS,nS),Bph(nS,nS),KA_sta(nS,nS),KB_sta(nS,nS), &
           OmBSE(nS),XpY_BSE(nS,nS),XmY_BSE(nS,nS))

!-----!
! TDA !
!-----!

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated in phBSE!'
    write(*,*)
  end if

! Initialization

  EcBSE(:) = 0d0

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  isp_W = 1
  EcRPA = 0d0

                 call phLR_A(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,eW,ERI,Aph)
  if(.not.TDA_W) call phLR_B(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  call RGW_phBSE_static_kernel_A(eta,nOrb,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KA_sta)
  call RGW_phBSE_static_kernel_B(eta,nOrb,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KB_sta)

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    ispin = 1

    ! Compute BSE excitation energies

                 call phLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
    if(.not.TDA) call phLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

    ! Second-order BSE static kernel
  
    if(dophBSE2) then 

      allocate(W(nOrb,nOrb,nOrb,nOrb))

      write(*,*) 
      write(*,*) '*** Second-order BSE static kernel activated! ***'
      write(*,*) 

      call RGW_phBSE_static_kernel(eta,nOrb,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,W)
      call RGW_phBSE2_static_kernel_A(eta,nOrb,nC,nO,nV,nR,nS,1d0,eW,W,KA_sta)

      if(.not.TDA) call RGW_phBSE2_static_kernel_B(eta,nOrb,nC,nO,nV,nR,nS,1d0,eW,W,KB_sta)

      deallocate(W)

    end if

                 Aph(:,:) = Aph(:,:) + KA_sta(:,:)
    if(.not.TDA) Bph(:,:) = Bph(:,:) + KB_sta(:,:)

    call phLR(TDA,nS,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GW@RHF','singlet',nS,OmBSE)
    call phLR_transition_vectors(.true.,nOrb,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

    !----------------------------------------------------!
    ! Compute the dynamical screening at the phBSE level !
    !----------------------------------------------------!

    if(dBSE) &
        call RGW_phBSE_dynamic_perturbation(dophBSE2,dTDA,eta,nOrb,nC,nO,nV,nR,nS,eW,eGW,ERI,dipole_int,OmRPA,rho_RPA, &
                                           OmBSE,XpY_BSE,XmY_BSE,KA_sta,KB_sta)

    !----------------!
    ! Upfolded phBSE !
    !----------------!

!   call RGW_phBSE_upfolded_sym(ispin,nOrb,nOrb,nC,nO,nV,nR,nS,ERI,rho_RPA,OmRPA,eW)

!   call RGW_phBSE_upfolded(ispin,nOrb,nOrb,nC,nO,nV,nR,nS,ERI,rho_RPA,OmRPA,eGW)

  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    ispin = 2

    ! Compute BSE excitation energies

                 call phLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
    if(.not.TDA) call phLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

                 Aph(:,:) = Aph(:,:) + KA_sta(:,:)
    if(.not.TDA) Bph(:,:) = Bph(:,:) + KB_sta(:,:)

    call phLR(TDA,nS,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GW@RHF','triplet',nS,OmBSE)
    call phLR_transition_vectors(.false.,nOrb,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) &
        call RGW_phBSE_dynamic_perturbation(dophBSE2,dTDA,eta,nOrb,nC,nO,nV,nR,nS,eW,eGW,ERI,dipole_int,OmRPA,rho_RPA, &
                                           OmBSE,XpY_BSE,XmY_BSE,KA_sta,KB_sta)

    !----------------!
    ! Upfolded phBSE !
    !----------------!

!   call RGW_phBSE_upfolded(ispin,nOrb,nOrb,nC,nO,nV,nR,nS,ERI,rho_RPA,OmRPA,eGW)

  end if

! Scale properly correlation energy if exchange is included in interaction kernel

  if(exchange_kernel) then

    EcBSE(1) = 0.5d0*EcBSE(1)
    EcBSE(2) = 1.5d0*EcBSE(2)

  end if

end subroutine 
