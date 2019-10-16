subroutine G0T0(BSE,singlet_manifold,triplet_manifold,eta,nBas,nC,nO,nV,nR,nOO,nVV,ENuc,ERHF,Hc,H,ERI,PHF,cHF,eHF,eG0T0)

! Perform G0W0 calculation with a T-matrix self-energy (G0T0)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: BSE
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  double precision,intent(in)   :: eta

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nOO,nVV
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: H(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dRPA
  integer                       :: ispin,jspin
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: SigT(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Omega1(:,:)
  double precision,allocatable  :: X1(:,:,:)
  double precision,allocatable  :: Y1(:,:,:)
  double precision,allocatable  :: rho1(:,:,:,:)
  double precision,allocatable  :: Omega2(:,:)
  double precision,allocatable  :: X2(:,:,:)
  double precision,allocatable  :: Y2(:,:,:)
  double precision,allocatable  :: rho2(:,:,:,:)

! Output variables

  double precision              :: eG0T0(nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0T0 calculation           |'
  write(*,*)'************************************************'
  write(*,*)

! Spin manifold 

  ispin = 1

! Memory allocation

  allocate(Omega1(nVV,nspin),X1(nVV,nVV,nspin),Y1(nOO,nVV,nspin), & 
           Omega2(nOO,nspin),X2(nVV,nOO,nspin),Y2(nOO,nOO,nspin), & 
           rho1(nBas,nBas,nVV,nspin),rho2(nBas,nBas,nOO,nspin),   & 
           SigT(nBas),Z(nBas))

! Compute linear response

  call linear_response_pp(ispin,.false.,nBas,nC,nO,nV,nR,nOO,nVV,eHF,ERI, & 
                          Omega1(:,ispin),X1(:,:,ispin),Y1(:,:,ispin),    & 
                          Omega2(:,ispin),X2(:,:,ispin),Y2(:,:,ispin),    & 
                          EcRPA(ispin))

! Compute excitation densities for the T-matrix

  call excitation_density_Tmatrix(nBas,nC,nO,nR,nOO,nVV,ERI,                     & 
                                  X1(:,:,ispin),Y1(:,:,ispin),rho1(:,:,:,ispin), & 
                                  X2(:,:,ispin),Y2(:,:,ispin),rho2(:,:,:,ispin))

! Compute T-matrix version of the self-energy 

  call self_energy_Tmatrix_diag(eta,nBas,nC,nO,nV,nR,nOO,nVV,eHF,  & 
                                Omega1(:,ispin),rho1(:,:,:,ispin), &  
                                Omega2(:,ispin),rho2(:,:,:,ispin), & 
                                SigT)

! Compute renormalization factor for T-matrix self-energy

  call renormalization_factor_Tmatrix(eta,nBas,nC,nO,nV,nR,nOO,nVV,eHF,  & 
                                      Omega1(:,ispin),rho1(:,:,:,ispin), & 
                                      Omega2(:,ispin),rho2(:,:,:,ispin), & 
                                      Z(:))

! Solve the quasi-particle equation

  eG0T0(:) = eHF(:) + Z(:)*SigT(:)

! Dump results

  call print_excitation('pp-RPA (N+2)',ispin,nVV,Omega1(:,ispin))
  call print_excitation('pp-RPA (N-2)',ispin,nOO,Omega2(:,ispin))

  call print_G0T0(nBas,nO,eHF,ENuc,ERHF,SigT,Z,eG0T0,EcRPA(ispin))

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
