subroutine Bethe_Salpeter_Tmatrix(TDA_T,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eT,eGT,EcBSE)

! Compute the Bethe-Salpeter excitation energies with the T-matrix kernel

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA_T
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

  double precision,intent(in)   :: eT(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: nOOs
  integer                       :: nOOt
  integer                       :: nVVs
  integer                       :: nVVt
  integer                       :: ispin
  integer                       :: iblock
  integer                       :: dERI
  integer                       :: xERI

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: Omega1s(:),Omega1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Omega2s(:),Omega2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)

  double precision,allocatable  :: TA(:,:),TB(:,:)
  double precision,allocatable  :: OmBSE(:,:)
  double precision,allocatable  :: XpY_BSE(:,:,:)
  double precision,allocatable  :: XmY_BSE(:,:,:)

! Output variables

  double precision,intent(out)  :: EcBSE(nspin)

! Dimensions of the pp-RPA linear reponse matrices

  nOOs = nO*nO
  nVVs = nV*nV

  nOOt = nO*(nO - 1)/2
  nVVt = nV*(nV - 1)/2

! Memory allocation

  allocate(Omega1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs), &
           Omega2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs), &
           rho1s(nBas,nO,nVVs),rho2s(nBas,nV,nOOs), &
           Omega1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt), &
           Omega2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt), &
           rho1t(nBas,nO,nVVt),rho2t(nBas,nV,nOOt))

  allocate(TA(nS,nS),TB(nS,nS),OmBSE(nS,nspin),XpY_BSE(nS,nS,nspin),XmY_BSE(nS,nS,nspin))

! Initialize T matrix

  TA(:,:) = 0d0
  TB(:,:) = 0d0

!----------------------------------------------
! Compute T-matrix for alpha-beta block
!----------------------------------------------

  ispin  = 1
  iblock = 3
  dERI   = +1d0
  xERI   = +0d0

  call linear_response_pp(iblock,.true.,.false.,nBas,nC,nO,nV,nR,nOOs,nVVs,eT,ERI,  &
                          Omega1s,X1s,Y1s,Omega2s,X2s,Y2s,EcRPA(ispin))

  call excitation_density_Tmatrix(iblock,dERI,xERI,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI, &
                                  X1s,Y1s,rho1s,X2s,Y2s,rho2s)

               call static_Tmatrix_TA(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,1d0,ERI,Omega1s,rho1s,Omega2s,rho2s,TA)
  if(.not.TDA) call static_Tmatrix_TB(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,1d0,ERI,Omega1s,rho1s,Omega2s,rho2s,TB)

!----------------------------------------------
! Compute T-matrix for alpha-alpha block
!----------------------------------------------

  ispin  = 2
  iblock = 4
  dERI   = +1d0
  xERI   = -1d0

  call linear_response_pp(iblock,.true.,.false.,nBas,nC,nO,nV,nR,nOOt,nVVt,eT,ERI,  &
                          Omega1t,X1t,Y1t,Omega2t,X2t,Y2t,EcRPA(ispin))

  call excitation_density_Tmatrix(iblock,dERI,xERI,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI, &
                                  X1t,Y1t,rho1t,X2t,Y2t,rho2t)

               call static_Tmatrix_TA(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,1d0,ERI,Omega1t,rho1t,Omega2t,rho2t,TA)
  if(.not.TDA) call static_Tmatrix_TB(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,1d0,ERI,Omega1t,rho1t,Omega2t,rho2t,TB)

!-------------------
! Singlet manifold
!-------------------

 if(singlet) then

    ispin = 1
    EcBSE(ispin) = 0d0

    ! Compute BSE singlet excitation energies

    call linear_response_Tmatrix(ispin,.true.,TDA,eta,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI,TA,TB, &
                                 EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
    call print_excitation('BSE@GT      ',ispin,nS,OmBSE(:,ispin))
    call print_transition_vectors(.true.,nBas,nC,nO,nV,nR,nS,dipole_int, & 
                                  OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) then

      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)

      if(evDyn) then

!       call Bethe_Salpeter_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,dipole_int,OmRPA,rho_RPA, &
!                                                          OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      else

!       call Bethe_Salpeter_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,eW,eGW,dipole_int,OmRPA,rho_RPA, &
!                                                OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      end if

    end if

  end if

!-------------------
! Triplet manifold
!-------------------

 if(triplet) then

    ispin = 2
    EcBSE(ispin) = 0d0

    ! Compute BSE triplet excitation energies

    call linear_response_Tmatrix(ispin,.true.,TDA,eta,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI,TA,TB, &
                                 EcBSE(ispin),OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
    call print_excitation('BSE@GT      ',ispin,nS,OmBSE(:,ispin))
    call print_transition_vectors(.false.,nBas,nC,nO,nV,nR,nS,dipole_int, & 
                                  OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))

    !-------------------------------------------------
    ! Compute the dynamical screening at the BSE level
    !-------------------------------------------------

    if(dBSE) then

      ! Compute dynamic correction for BSE via perturbation theory (iterative or renormalized)

      if(evDyn) then

!       call Bethe_Salpeter_dynamic_perturbation_iterative(dTDA,eta,nBas,nC,nO,nV,nR,nS,eGW,dipole_int,OmRPA,rho_RPA, &
!                                                          OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      else

!       call Bethe_Salpeter_dynamic_perturbation(dTDA,eta,nBas,nC,nO,nV,nR,nS,eW,eGW,dipole_int,OmRPA,rho_RPA, &
!                                                OmBSE(:,ispin),XpY_BSE(:,:,ispin),XmY_BSE(:,:,ispin))
      end if

    end if

  end if

end subroutine Bethe_Salpeter_Tmatrix
