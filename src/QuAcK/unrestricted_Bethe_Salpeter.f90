subroutine unrestricted_Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,spin_conserved,spin_flip,eta,  & 
                                       nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab, & 
                                       eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: spin_conserved 
  logical,intent(in)            :: spin_flip

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: eW(nBas,nspin)
  double precision,intent(in)   :: eGW(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_abab(nBas,nBas,nBas,nBas)


! Local variables

  integer                       :: ispin
  integer                       :: isp_W

  integer                       :: nS_aa,nS_bb,nS_sc
  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:,:)

  double precision,allocatable  :: OmBSE_sc(:)
  double precision,allocatable  :: XpY_BSE_sc(:,:)
  double precision,allocatable  :: XmY_BSE_sc(:,:)

  integer                       :: nS_ab,nS_ba,nS_sf
  double precision,allocatable  :: OmBSE_sf(:)
  double precision,allocatable  :: XpY_BSE_sf(:,:)
  double precision,allocatable  :: XmY_BSE_sf(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

  ! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb
  
  allocate(OmRPA(nS_sc),XpY_RPA(nS_sc,nS_sc),XmY_RPA(nS_sc,nS_sc),rho_RPA(nBas,nBas,nS_sc,nspin))

!--------------------------!
! Spin-conserved screening !
!--------------------------!

  isp_W = 1
  EcRPA = 0d0

  ! Compute spin-conserved RPA screening 

  call unrestricted_linear_response(isp_W,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
                                    eW,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

  call unrestricted_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)


!----------------------------!
! Spin-conserved excitations !
!----------------------------!

  if(spin_conserved) then

    ispin = 1
    EcBSE(ispin) = 0d0

    allocate(OmBSE_sc(nS_sc),XpY_BSE_sc(nS_sc,nS_sc),XmY_BSE_sc(nS_sc,nS_sc))

    ! Compute spin-conserved BSE excitation energies

    call unrestricted_linear_response(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
                                      eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,OmRPA,rho_RPA,EcBSE(ispin),       & 
                                      OmBSE_sc,XpY_BSE_sc,XmY_BSE_sc)

    call print_excitation('BSE@UG0W0',5,nS_sc,OmBSE_sc)

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

!   if(dBSE) then

!     ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)
!
!     if(evDyn) then
!
!       call Bethe_Salpeter_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW(:),OmRPA(:,ispin),OmBSE(:,ispin), & 
!                                                          XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin),rho_RPA(:,:,:,ispin))
!     else
!
!       call Bethe_Salpeter_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW(:),OmRPA(:,ispin),OmBSE(:,ispin), & 
!                                                XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin),rho_RPA(:,:,:,ispin))
!     end if

!   end if
 
  end if

!-----------------------!
! Spin-flip excitations !
!-----------------------!

 if(spin_flip) then

    ispin = 2
    EcBSE(ispin) = 0d0

    ! Memory allocation

    nS_ab = (nO(1) - nC(1))*(nV(2) - nR(2))
    nS_ba = (nO(2) - nC(2))*(nV(1) - nR(1))
    nS_sf = nS_ab + nS_ba
  
    allocate(OmBSE_sf(nS_sf),XpY_BSE_sf(nS_sf,nS_sf),XmY_BSE_sf(nS_sf,nS_sf))

    call unrestricted_linear_response(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS_ab,nS_ba,nS_sf,nS_sc,1d0, &
                                      eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,OmRPA,rho_RPA,EcBSE(ispin),       & 
                                      OmBSE_sf,XpY_BSE_sf,XmY_BSE_sf)

    call print_excitation('BSE@UG0W0',6,nS_sf,OmBSE_sf)

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

!   if(dBSE) then

!     ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)

!     if(evDyn) then
!    
!       call Bethe_Salpeter_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA(:,ispin),OmBSE(:,ispin), & 
!                                                          XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin),rho_RPA(:,:,:,ispin))
!     else
!    
!       call Bethe_Salpeter_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA(:,ispin),OmBSE(:,ispin), & 
!                                                XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin),rho_RPA(:,:,:,ispin))
!     end if

!   end if

  end if

end subroutine unrestricted_Bethe_Salpeter
