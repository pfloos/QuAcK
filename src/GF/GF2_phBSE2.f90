subroutine GF2_phBSE2(TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,EcBSE)

! Compute the second-order Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

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
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: dRPA = .false.
  integer                       :: ispin
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: A_sta(:,:)
  double precision,allocatable  :: B_sta(:,:)
  double precision,allocatable  :: KA_sta(:,:)
  double precision,allocatable  :: KB_sta(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(OmBSE(nS),XpY(nS,nS),XmY(nS,nS),A_sta(nS,nS),KA_sta(nS,nS))
  allocate(B_sta(nS,nS),KB_sta(nS,nS))

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    ispin = 1
    EcBSE(ispin) = 0d0

                 call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGF,ERI,A_sta)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,B_sta)

    ! Compute static kernel

                 call GF2_phBSE2_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,KA_sta)
    if(.not.TDA) call GF2_phBSE2_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,KB_sta)

                 A_sta(:,:) = A_sta(:,:) + KA_sta(:,:)
    if(.not.TDA) B_sta(:,:) = B_sta(:,:) + KB_sta(:,:)

    ! Compute phBSE2@GF2 excitation energies

    call phLR(TDA,nS,A_sta,B_sta,EcBSE(ispin),OmBSE,XpY,XmY)
    call print_excitation_energies('phBSE2@GF2',ispin,nS,OmBSE)
    call phLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY,XmY)

    ! Compute dynamic correction for BSE via perturbation theory

    if(dBSE) &
      call GF2_phBSE2_dynamic_perturbation(dTDA,ispin,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,KA_sta,KB_sta,OmBSE,XpY,XmY)

  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    ispin = 2
    EcBSE(ispin) = 0d0

                 call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGF,ERI,A_sta)
    if(.not.TDA) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,B_sta)

    ! Compute static kernel

                 call GF2_phBSE2_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,KA_sta)
    if(.not.TDA) call GF2_phBSE2_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,KB_sta)

                 A_sta(:,:) = A_sta(:,:) + KA_sta(:,:)
    if(.not.TDA) B_sta(:,:) = B_sta(:,:) + KB_sta(:,:)

    ! Compute phBSE2@GF2 excitation energies

    call phLR(TDA,nS,A_sta,B_sta,EcBSE(ispin),OmBSE,XpY,XmY)
    call print_excitation_energies('phBSE2@GF2',ispin,nS,OmBSE)
    call phLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY,XmY)

    ! Compute dynamic correction for BSE via perturbation theory

    if(dBSE) &
      call GF2_phBSE2_dynamic_perturbation(dTDA,ispin,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,KA_sta,KB_sta,OmBSE,XpY,XmY)

  end if

end subroutine 
