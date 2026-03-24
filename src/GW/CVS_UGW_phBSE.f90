subroutine CVS_UGW_phBSE(exchange_kernel,TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS, &
                     nCVS,nFC,occupations,virtuals,S,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,cW,eW,eGW,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: spin_conserved 
  logical,intent(in)            :: spin_flip

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  integer,intent(in)            :: nCVS(nspin)
  integer,intent(in)            :: nFC(nspin)
  integer,intent(in)            :: occupations(maxval(nO-nFC),nspin)
  integer,intent(in)            :: virtuals(nBas-minval(nO),nspin)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: cW(nBas,nBas,nspin)
  double precision,intent(in)   :: eW(nBas,nspin)
  double precision,intent(in)   :: eGW(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin
  integer                       :: isp_W
  logical                       :: dRPA   = .false.
  logical                       :: dRPA_W = .true.

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: KA(:,:)
  double precision,allocatable  :: KB(:,:)

  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:,:)

  integer                       :: nS_aa,nS_bb,nS_sc
  integer                       :: nS_ab,nS_ba,nS_sf
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY_BSE(:,:)
  double precision,allocatable  :: XmY_BSE(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  nS_aa = (nO(1) - nFC(1))*(nV(1)-nCVS(1))
  nS_bb = (nO(2) - nFC(2))*(nV(2)-nCVS(2))
  nS_sc = nS_aa + nS_bb

  nS_ab = (nO(1) - nFC(1))*(nV(2) - nCVS(2))
  nS_ba = (nO(2) - nFC(2))*(nV(1) - nCVS(1))
  nS_sf = nS_ab + nS_ba
  
! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated in phBSE!'
    write(*,*)
  end if

! Initialization

  EcBSE(:) = 0d0

!--------------------------!
! Spin-conserved screening !
!--------------------------!

  isp_W = 1

  ! Compute spin-conserved RPA screening 

  allocate(Aph(nS_sc,nS_sc),Bph(nS_sc,nS_sc))
  allocate(OmRPA(nS_sc),XpY_RPA(nS_sc,nS_sc),XmY_RPA(nS_sc,nS_sc),rho_RPA(nBas,nBas,nS_sc,nspin))

                 call CVS_phULR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nCVS,nFC,occupations,virtuals,  &
                                    1d0,eW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
  if(.not.TDA_W) call CVS_phULR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nCVS,nFC,occupations,virtuals,  &
                                    1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)
  call CVS_phULR(TDA_W,nS_aa,nS_bb,nS_sc,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call CVS_UGW_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,nCVS,nFC,occupations,virtuals,  &
                                          ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)
    
  deallocate(Aph,Bph)

!----------------------------!
! Spin-conserved excitations !
!----------------------------!

  if(spin_conserved) then

    ispin = 1

    allocate(Aph(nS_sc,nS_sc),Bph(nS_sc,nS_sc),KA(nS_sc,nS_sc),KB(nS_sc,nS_sc))
    allocate(OmBSE(nS_sc),XpY_BSE(nS_sc,nS_sc),XmY_BSE(nS_sc,nS_sc))

    ! Compute spin-conserved BSE excitation energies

                 call CVS_phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nCVS,nFC,occupations,virtuals,  &
                                   1d0,eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
    if(.not.TDA) call CVS_phULR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nCVS,nFC,occupations,virtuals,  &
                                   1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

                 call CVS_UGW_phBSE_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,         &
                                   nCVS,nFC,occupations,virtuals,                                               &
                                   1d0,OmRPA,rho_RPA,KA)
    if(.not.TDA) call CVS_UGW_phBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,         &
                                   nCVS,nFC,occupations,virtuals,                                               &
                                   1d0,OmRPA,rho_RPA,KB)

                 Aph(:,:) = Aph(:,:) + KA(:,:)
    if(.not.TDA) Bph(:,:) = Bph(:,:) + KB(:,:)

    call CVS_phULR(TDA,nS_aa,nS_bb,nS_sc,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GW@UHF','spin-conserved',nS_sc,OmBSE)
    call CVS_phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nS_aa,nS_bb,nS_sc,nCVS,nFC,occupations,virtuals, &
                                  dipole_int_aa,dipole_int_bb,                                                   &
                                  cW,S,OmBSE,XpY_BSE,XmY_BSE)

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) &
      call CVS_UGW_phBSE_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nS_aa,nS_bb,nS_sc,nS_sc,    &
                                               nCVS,nFC,occupations,virtuals,                                &
                                          eW,eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,     &
                                          OmRPA,rho_RPA,OmBSE,XpY_BSE,XmY_BSE)

    deallocate(Aph,Bph,KA,KB)
    deallocate(OmBSE,XpY_BSE,XmY_BSE)

  end if

!-----------------------!
! Spin-flip excitations !
!-----------------------!

 if(spin_flip) then

    ispin = 2

    allocate(Aph(nS_sf,nS_sf),Bph(nS_sf,nS_sf),KA(nS_sf,nS_sf),KB(nS_sf,nS_sf))
    allocate(OmBSE(nS_sf),XpY_BSE(nS_sf,nS_sf),XmY_BSE(nS_sf,nS_sf))

    ! Compute spin-conserved BSE excitation energies

                 call CVS_phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS_ab,nS_ba,nS_sf,nCVS,nFC,occupations,virtuals,  &
                                1d0,eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
    if(.not.TDA) call CVS_phULR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS_ab,nS_ba,nS_sf,nCVS,nFC,occupations,virtuals,  &
                                1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

                 call CVS_UGW_phBSE_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS_ab,nS_ba,nS_sf,nS_sc,         &
                                                     nCVS,nFC,occupations,virtuals,                             &
                                                      1d0,OmRPA,rho_RPA,KA)
    if(.not.TDA) call CVS_UGW_phBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS_ab,nS_ba,nS_sf,nS_sc,         &
                                                      nCVS,nFC,occupations,virtuals,                            &
                                                      1d0,OmRPA,rho_RPA,KB)

                 Aph(:,:) = Aph(:,:) + KA(:,:)
    if(.not.TDA) Bph(:,:) = Bph(:,:) + KB(:,:)

    call CVS_phULR(TDA,nS_ab,nS_ba,nS_sf,Aph,Bph,EcBSE(ispin),OmBSE,XpY_BSE,XmY_BSE)

    call print_excitation_energies('phBSE@GW@UHF','spin-flip',nS_sf,OmBSE)
    call CVS_phULR_transition_vectors(ispin,nBas,nC,nO,nV,nR,nS,nS_ab,nS_ba,nS_sf,nCVS,nFC,occupations,virtuals,&
                                      dipole_int_aa,dipole_int_bb,                                              &
                                      cW,S,OmBSE,XpY_BSE,XmY_BSE)

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) &
      call CVS_UGW_phBSE_dynamic_perturbation(ispin,dTDA,eta,nBas,nC,nO,nV,nR,nS,nS_ab,nS_ba,nS_sf,nS_sc,    &
                                               nCVS,nFC,occupations,virtuals,                                &
                                          eW,eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,     &
                                          OmRPA,rho_RPA,OmBSE,XpY_BSE,XmY_BSE)

    deallocate(Aph,Bph,KA,KB)
    deallocate(OmBSE,XpY_BSE,XmY_BSE)
 
  end if
  
! Scale properly correlation energy if exchange is included in interaction kernel

  if(exchange_kernel) then

    EcBSE(:) = 0.5d0*EcBSE(:)

  end if

end subroutine 
