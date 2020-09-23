subroutine unrestricted_Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,spin_conserved,spin_flip,eta,  & 
                                       nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab, & 
                                       eW,eGW,EcRPA,EcBSE)

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
  integer                       :: nS_ab,nS_ba,nS_sf
  double precision,allocatable  :: OmRPA_sc(:)
  double precision,allocatable  :: XpY_RPA_sc(:,:)
  double precision,allocatable  :: XmY_RPA_sc(:,:)
  double precision,allocatable  :: rho_RPA_sc(:,:,:,:)
  double precision,allocatable  :: OmBSE_sc(:)
  double precision,allocatable  :: XpY_BSE_sc(:,:)
  double precision,allocatable  :: XmY_BSE_sc(:,:)

! Output variables

  double precision,intent(out)  :: EcRPA(nspin)
  double precision,intent(out)  :: EcBSE(nspin)

!----------------------------!
! Spin-conserved excitations !
!----------------------------!

 if(spin_conserved) then

    ispin = 1
    isp_W = 1
    EcBSE(ispin) = 0d0

    ! Memory allocation

    nS_aa = nS(1)
    nS_bb = nS(2)
    nS_sc = nS_aa + nS_bb
  
    allocate(OmRPA_sc(nS_sc),XpY_RPA_sc(nS_sc,nS_sc),XmY_RPA_sc(nS_sc,nS_sc),rho_RPA_sc(nBas,nBas,nS_sc,nspin))
    allocate(OmBSE_sc(nS_sc),XpY_BSE_sc(nS_sc,nS_sc),XmY_BSE_sc(nS_sc,nS_sc))

    ! Compute spin-conserved RPA screening 

    call unrestricted_linear_response(isp_W,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,1d0, &
                                      eW,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,rho_RPA_sc,EcRPA(ispin),        & 
                                      OmRPA_sc,XpY_RPA_sc,XmY_RPA_sc)

!   call print_excitation('RPA@UG0W0',5,nS_sc,OmRPA_sc)

    call unrestricted_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb, & 
                                         XpY_RPA_sc,rho_RPA_sc)

    ! Compute spin-conserved BSE excitation energies

    OmBSE_sc(:) = OmRPA_sc(:)

    call unrestricted_linear_response(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,1d0, &
                                      eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,rho_RPA_sc,EcBSE(ispin),    & 
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

!if(spin_flip) then

!   ispin = 2
!   isp_W = 1
!   EcBSE(ispin) = 0d0

!  ! Compute (singlet) RPA screening 

!   call linear_response(isp_W,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI, &
!                        rho_RPA(:,:,:,ispin),EcRPA(ispin),OmRPA(:,ispin),XpY_RPA(:,:,ispin),XmY_RPA(:,:,ispin))
!   call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA(:,:,ispin),rho_RPA(:,:,:,ispin))

!   ! Compute BSE excitation energies

!   OmBSE(:,ispin) = OmRPA(:,ispin)

!   call linear_response(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI, &
!                        rho_RPA(:,:,:,ispin),EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
!   call print_excitation('BSE         ',ispin,nS,OmBSE(:,ispin))

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

! end if

end subroutine unrestricted_Bethe_Salpeter
