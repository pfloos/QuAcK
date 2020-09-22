subroutine unrestricted_Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,spin_conserved,spin_flip,eta,  & 
                                       nBas,nC,nO,nV,nR,nSa,nSb,nSt,ERI_aaaa,ERI_aabb,ERI_bbbb, & 
                                       eW,eGW,OmRPA,XpY_RPA,XmY_RPA,rho_RPA,EcRPA,EcBSE)

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
  integer,intent(in)            :: nSa
  integer,intent(in)            :: nSb
  integer,intent(in)            :: nSt
  double precision,intent(in)   :: eW(nBas,nspin)
  double precision,intent(in)   :: eGW(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)

  double precision              :: OmRPA(nSt)
  double precision              :: XpY_RPA(nSt,nSt)
  double precision              :: XmY_RPA(nSt,nSt)
  double precision              :: rho_RPA(nBas,nBas,nSt,nspin)

! Local variables

  integer                       :: ispin
  integer                       :: isp_W
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY_BSE(:,:)
  double precision,allocatable  :: XmY_BSE(:,:)

! Output variables

  double precision,intent(out)  :: EcRPA
  double precision,intent(out)  :: EcBSE

! Memory allocation

  allocate(OmBSE(nSt),XpY_BSE(nSt,nSt),XmY_BSE(nSt,nSt))

!----------------------------!
! Spin-conserved excitations !
!----------------------------!

 if(spin_conserved) then

    ispin = 1
    isp_W = 1
    EcBSE = 0d0

   ! Compute spin-conserved RPA screening 

    call unrestricted_linear_response(isp_W,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,1d0, &
                                      eW,ERI_aaaa,ERI_aabb,ERI_bbbb,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

    call unrestricted_excitation_density(nBas,nC,nO,nR,nSa,nSb,nSt,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

    ! Compute BSE excitation energies

    OmBSE(:) = OmRPA(:)

    call unrestricted_linear_response(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,1d0, &
                                      eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,rho_RPA,EcBSE,OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation('BSE@UG0W0',5,nSt,OmBSE)

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
