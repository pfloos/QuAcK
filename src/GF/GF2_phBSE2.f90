subroutine GF2_phBSE2(TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGF,EcBSE)

! Compute the second-order Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

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
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin
  double precision,allocatable  :: OmBSE(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision              :: rho
  double precision,allocatable  :: A_sta(:,:)
  double precision,allocatable  :: B_sta(:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Memory allocation

  allocate(OmBSE(nS),XpY(nS,nS),XmY(nS,nS),A_sta(nS,nS),B_sta(nS,nS))

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    ispin = 1
    EcBSE(ispin) = 0d0

    ! Compute static kernel

                 call GF2_phBSE2_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,A_sta)
    if(.not.TDA) call GF2_phBSE2_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,B_sta)

    ! Compute BSE2 excitation energies

    call linear_response_BSE(ispin,.false.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGF,ERI,-A_sta,-B_sta,EcBSE(ispin),OmBSE,XpY,XmY)
    call print_excitation('BSE2        ',ispin,nS,OmBSE)
    call print_transition_vectors_ph(.true.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY,XmY)


    ! Compute dynamic correction for BSE via perturbation theory

    if(dBSE) then

     if(evDyn) then

      call GF2_phBSE2_dynamic_perturbation_iterative(dTDA,ispin,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGF, & 
                                                     A_sta,B_sta,OmBSE,XpY,XmY)
     else
 
      call GF2_phBSE2_dynamic_perturbation(dTDA,ispin,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGF,A_sta,B_sta,OmBSE,XpY,XmY)
 
      end if

    end if

  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    ispin = 2
    EcBSE(ispin) = 0d0

    ! Compute static kernel

                 call GF2_phBSE2_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,A_sta)
    if(.not.TDA) call GF2_phBSE2_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,eGF,B_sta)

    ! Compute BSE2 excitation energies

    call linear_response_BSE(ispin,.false.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGF,ERI,-A_sta,-B_sta,EcBSE(ispin),OmBSE,XpY,XmY)
    call print_excitation('BSE2        ',ispin,nS,OmBSE)
    call print_transition_vectors_ph(.false.,nBas,nC,nO,nV,nR,nS,dipole_int,OmBSE,XpY,XmY)

    ! Compute dynamic correction for BSE via perturbation theory

    if(dBSE) then

     if(evDyn) then

      call GF2_phBSE2_dynamic_perturbation_iterative(dTDA,ispin,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGF, & 
                                                     A_sta,B_sta,OmBSE,XpY,XmY)
     else
 
      call GF2_phBSE2_dynamic_perturbation(dTDA,ispin,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGF,A_sta,B_sta,OmBSE,XpY,XmY)
 
      end if

    end if

  end if

end subroutine 