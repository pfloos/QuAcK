subroutine Bethe_Salpeter(TDA,singlet_manifold,triplet_manifold,eta, & 
                 nBas,nC,nO,nV,nR,nS,ERI,eW,eGW,Omega,XpY,XmY,rho,EcRPA,EcBSE)

! Compute the Bethe-Salpeter excitation energies

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

  double precision              :: Omega(nS,nspin)
  double precision              :: XpY(nS,nS,nspin)
  double precision              :: XmY(nS,nS,nspin)
  double precision              :: rho(nBas,nBas,nS,nspin)

! Local variables

  integer                       :: ispin

! Output variables

  double precision,intent(out)  :: EcRPA(nspin)
  double precision,intent(out)  :: EcBSE(nspin)

! Singlet manifold

 if(singlet_manifold) then

    ispin = 1
    EcBSE(ispin) = 0d0

    call linear_response(ispin,.true.,TDA,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI, &
                         rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

    call linear_response(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI, &
                         rho(:,:,:,ispin),EcBSE(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('BSE   ',ispin,nS,Omega(:,ispin))

  end if

! Triplet manifold

 if(triplet_manifold) then

    ispin = 2
    EcBSE(ispin) = 0d0

    call linear_response(ispin,.true.,TDA,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI, &
                         rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

    call linear_response(ispin,.true.,TDA,.true.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI, &
                         rho(:,:,:,ispin),EcBSE(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('BSE   ',ispin,nS,Omega(:,ispin))

  end if

end subroutine Bethe_Salpeter
