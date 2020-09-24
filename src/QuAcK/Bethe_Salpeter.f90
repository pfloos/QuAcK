subroutine Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

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

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nBas,nBas,nS))
  allocate(OmBSE(nS,nspin),XpY_BSE(nS,nS,nspin),XmY_BSE(nS,nS,nspin))

!---------------------------------
! Compute (singlet) RPA screening 
!---------------------------------

  isp_W = 1
  EcRPA = 0d0

  call linear_response(isp_W,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI,OmRPA, &
                       rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    ispin = 1
    EcBSE(ispin) = 0d0

    ! Compute BSE excitation energies

    call linear_response(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI,OmRPA, &
                         rho_RPA,EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
    call print_excitation('BSE@GW      ',ispin,nS,OmBSE(:,ispin))

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) then

      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)
 
      if(evDyn) then
 
        call Bethe_Salpeter_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA,rho_RPA, &
                                                           OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      else
 
        call Bethe_Salpeter_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA,rho_RPA, &
                                                 OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
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

    call linear_response(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI,OmRPA, &
                         rho_RPA,EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
    call print_excitation('BSE@GW      ',ispin,nS,OmBSE(:,ispin))

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) then

      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)

      if(evDyn) then
     
        call Bethe_Salpeter_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA,rho_RPA, &
                                                           OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      else
     
        call Bethe_Salpeter_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA,rho_RPA, &
                                                 OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      end if

    end if

  end if

end subroutine Bethe_Salpeter
