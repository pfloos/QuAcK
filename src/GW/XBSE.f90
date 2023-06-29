subroutine XBSE(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

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

  integer                       :: ispin
  integer                       :: isp_W

  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

  double precision,allocatable  :: OmBSE(:,:)
  double precision,allocatable  :: XpY_BSE(:,:,:)
  double precision,allocatable  :: XmY_BSE(:,:,:)

  double precision,allocatable  :: KA_sta(:,:)
  double precision,allocatable  :: KB_sta(:,:)

  double precision,allocatable  :: W(:,:,:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nBas,nBas,nS), &
           KA_sta(nS,nS),KB_sta(nS,nS),OmBSE(nS,nspin),XpY_BSE(nS,nS,nspin),XmY_BSE(nS,nS,nspin))

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  isp_W = 1
  EcRPA = 0d0

  call linear_response(isp_W,.true.,TDA_W,eta,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI, &
                       EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

  call XBSE_static_kernel_KA(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KA_sta,eW,eGW)
  call XBSE_static_kernel_KB(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KB_sta,eW)

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    ispin = 1
    EcBSE(ispin) = 0d0

    ! Compute BSE excitation energies

    call linear_response_BSE(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI,KA_sta,KB_sta, &
                             EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
    call print_excitation('BSE@GW      ',ispin,nS,OmBSE(:,ispin))
    call print_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int, & 
                                  OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) then

      call Bethe_Salpeter_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,eW,eW,dipole_int,OmRPA,rho_RPA, &
                                               OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
    end if
 
  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    ispin = 2
    EcBSE(ispin) = 0d0

    ! Compute BSE excitation energies

    call linear_response_BSE(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI,KA_sta,KB_sta, &
                             EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
    call print_excitation('BSE@GW      ',ispin,nS,OmBSE(:,ispin))
    call print_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int, & 
                                  OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) then

      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)

      call Bethe_Salpeter_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,eW,eW,dipole_int,OmRPA,rho_RPA, &
                                               OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    end if

  end if

end subroutine 
