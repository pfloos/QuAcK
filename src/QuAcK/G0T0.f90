subroutine G0T0(BSE,singlet_manifold,triplet_manifold,eta,nBas,nC,nO,nV,nR,ENuc,ERHF,ERI,eHF,eG0T0)

! Perform G0W0 calculation with a T-matrix self-energy (G0T0)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: BSE
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  double precision,intent(in)   :: eta

  integer,intent(in)            :: nBas,nC,nO,nV,nR
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin
  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision,allocatable  :: Omega1s(:),Omega1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Omega2s(:),Omega2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)
  double precision,allocatable  :: SigT(:)
  double precision,allocatable  :: Z(:)

! Output variables

  double precision,intent(out)  :: eG0T0(nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0T0 calculation           |'
  write(*,*)'************************************************'
  write(*,*)

! Dimensions of the rr-RPA linear reponse matrices

  nOOs = nO*(nO + 1)/2
  nVVs = nV*(nV + 1)/2

  nOOt = nO*(nO - 1)/2
  nVVt = nV*(nV - 1)/2

! Memory allocation

  allocate(Omega1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs), & 
           Omega2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs), & 
           rho1s(nBas,nO-nC,nVVs),rho2s(nBas,nV-nR,nOOs), & 
           Omega1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt), & 
           Omega2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt), & 
           rho1t(nBas,nO-nC,nVVt),rho2t(nBas,nV-nR,nOOt), & 
           SigT(nBas),Z(nBas))

!----------------------------------------------
! Singlet manifold
!----------------------------------------------

 ispin = 1

! Compute linear response

  call linear_response_pp(ispin,.false.,nBas,nC,nO,nV,nR, & 
                          nOOs,nVVs,eHF(:),ERI(:,:,:,:),  & 
                          Omega1s(:),X1s(:,:),Y1s(:,:),   & 
                          Omega2s(:),X2s(:,:),Y2s(:,:),   & 
                          EcRPA(ispin))

  call print_excitation('pp-RPA (N+2)',ispin,nVVs,Omega1s(:))
  call print_excitation('pp-RPA (N-2)',ispin,nOOs,Omega2s(:))

! Compute excitation densities for the T-matrix

  call excitation_density_Tmatrix(ispin,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI(:,:,:,:), & 
                                  X1s(:,:),Y1s(:,:),rho1s(:,:,:),                & 
                                  X2s(:,:),Y2s(:,:),rho2s(:,:,:))

!----------------------------------------------
! Triplet manifold
!----------------------------------------------

 ispin = 2

! Compute linear response

  call linear_response_pp(ispin,.false.,nBas,nC,nO,nV,nR, & 
                          nOOt,nVVt,eHF(:),ERI(:,:,:,:),  & 
                          Omega1t(:),X1t(:,:),Y1t(:,:),   & 
                          Omega2t(:),X2t(:,:),Y2t(:,:),   & 
                          EcRPA(ispin))

  call print_excitation('pp-RPA (N+2)',ispin,nVVt,Omega1t(:))
  call print_excitation('pp-RPA (N-2)',ispin,nOOt,Omega2t(:))

! Compute excitation densities for the T-matrix

  call excitation_density_Tmatrix(ispin,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI(:,:,:,:), & 
                                  X1t(:,:),Y1t(:,:),rho1t(:,:,:),                & 
                                  X2t(:,:),Y2t(:,:),rho2t(:,:,:))

!----------------------------------------------
! Compute T-matrix version of the self-energy 
!----------------------------------------------

  call self_energy_Tmatrix_diag(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF(:), & 
                                Omega1s(:),rho1s(:,:,:),Omega2s(:),rho2s(:,:,:), & 
                                Omega1t(:),rho1t(:,:,:),Omega2t(:),rho2t(:,:,:), & 
                                SigT(:))

! Compute renormalization factor for T-matrix self-energy

  call renormalization_factor_Tmatrix(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF(:), & 
                                      Omega1s(:),rho1s(:,:,:),Omega2s(:),rho2s(:,:,:), & 
                                      Omega1t(:),rho1t(:,:,:),Omega2t(:),rho2t(:,:,:), & 
                                      Z(:))

!----------------------------------------------
! Solve the quasi-particle equation
!----------------------------------------------

  eG0T0(:) = eHF(:) + Z(:)*SigT(:)

!----------------------------------------------
! Dump results
!----------------------------------------------

  call print_G0T0(nBas,nO,eHF(:),ENuc,ERHF,SigT(:),Z(:),eG0T0(:),EcRPA(:))

! Perform BSE calculation

! if(BSE) then

!  ! Singlet manifold

!  if(singlet_manifold) then

!     ispin = 1
!     EcBSE(ispin) = 0d0

!     call linear_response(ispin,.false.,.false.,.false.,nBas,nC,nO,nV,nR,nS,eHF,ERI, &
!                          rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin))
!     call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

!     call linear_response(ispin,.false.,.false.,BSE,nBas,nC,nO,nV,nR,nS,eG0T0,ERI, &
!                          rho(:,:,:,ispin),EcBSE(ispin),Omega(:,ispin),XpY(:,:,ispin))
!     call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

!   endif

!  ! Triplet manifold

!  if(triplet_manifold) then

!     ispin = 2
!     EcBSE(ispin) = 0d0

!     call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,eHF,ERI, &
!                          rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin))
!     call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

!     call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,eG0T0,ERI, &
!                          rho(:,:,:,ispin),EcBSE(ispin),Omega(:,ispin),XpY(:,:,ispin))
!     call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

!   endif

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A40,F15.6)') 'BSE@G0T0 correlation energy (singlet) =',EcBSE(1)
!   write(*,'(2X,A40,F15.6)') 'BSE@G0T0 correlation energy (triplet) =',EcBSE(2)
!   write(*,'(2X,A40,F15.6)') 'BSE@G0T0 correlation energy           =',EcBSE(1) + EcBSE(2)
!   write(*,'(2X,A40,F15.6)') 'BSE@G0T0 total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

! endif

end subroutine G0T0
