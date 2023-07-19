subroutine GW_phBSE(BSE2,TDA_W,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: BSE2
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
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

  double precision,allocatable  :: A_sta(:,:)
  double precision,allocatable  :: B_sta(:,:)

  double precision,allocatable  :: KA_sta(:,:)
  double precision,allocatable  :: KB_sta(:,:)

  double precision,allocatable  :: W(:,:,:,:)
  double precision,allocatable  :: KA2_sta(:,:)
  double precision,allocatable  :: KB2_sta(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nBas,nBas,nS), &
           A_sta(nS,nS),KA_sta(nS,nS),KA2_sta(nS,nS),OmBSE(nS),XpY_BSE(nS,nS),XmY_BSE(nS,nS))
  allocate(B_sta(nS,nS),KB_sta(nS,nS),KB2_sta(nS,nS))

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  isp_W = 1
  EcRPA = 0d0

                 call phLR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI,A_sta)
  if(.not.TDA_W) call phLR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,1d0,ERI,B_sta)

  call phLR(TDA_W,nS,A_sta,B_sta,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
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

                 call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI,A_sta)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,B_sta)

                 A_sta(:,:) = A_sta(:,:) + KA_sta(:,:)
    if(.not.TDA) B_sta(:,:) = B_sta(:,:) + KB_sta(:,:)

    ! Second-order BSE static kernel
  
    if(BSE2) then 

      allocate(W(nBas,nBas,nBas,nBas))

      write(*,*) 
      write(*,*) '*** Second-order BSE static kernel activated! ***'
      write(*,*) 

      call static_kernel_W(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,W)
      call GW_phBSE2_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,1d0,eW,W,KA2_sta)

      if(.not.TDA) call GW_phBSE2_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,1d0,eW,W,KB2_sta)

      deallocate(W)

                   A_sta(:,:) = A_sta(:,:) + KA2_sta(:,:)
      if(.not.TDA) B_sta(:,:) = B_sta(:,:) + KB2_sta(:,:)

    end if

    call phLR(TDA,nS,A_sta,B_sta,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation('phBSE@GW    ',ispin,nS,OmBSE)
    call print_transition_vectors_ph(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) then

      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)
 
      if(evDyn) then
 
        call GW_phBSE_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,dipole_int,OmRPA,rho_RPA, &
                                                     OmBSE,XpY_BSE,XmY_BSE)
      else

        call GW_phBSE_dynamic_perturbation(BSE2,dTDA,eta,nBas,nC,nO,nV,nR,nS,eW,eGW,dipole_int,OmRPA,rho_RPA, &
                                           OmBSE,XpY_BSE,XmY_BSE,W,KA2_sta)
      end if

    end if
 
  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    ispin = 2
    EcBSE(ispin) = 0d0

    ! Compute BSE excitation energies

                 call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI,A_sta)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,B_sta)

                 A_sta(:,:) = A_sta(:,:) + KA_sta(:,:)
    if(.not.TDA) B_sta(:,:) = B_sta(:,:) + KB_sta(:,:)

    call phLR(TDA,nS,A_sta,B_sta,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation('phBSE@GW    ',ispin,nS,OmBSE)
    call print_transition_vectors_ph(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY_BSE,XmY_BSE)

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) then

      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)

      if(evDyn) then
     
        call GW_phBSE_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,dipole_int,OmRPA,rho_RPA, &
                                                     OmBSE,XpY_BSE,XmY_BSE)
      else
     
        call GW_phBSE_dynamic_perturbation(BSE2,dTDA,eta,nBas,nC,nO,nV,nR,nS,eW,eGW,dipole_int,OmRPA,rho_RPA, &
                                           OmBSE,XpY_BSE,XmY_BSE,W,KA2_sta)
      end if

    end if

  end if

end subroutine 
