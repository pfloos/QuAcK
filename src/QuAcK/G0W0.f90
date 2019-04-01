subroutine G0W0(COHSEX,SOSEX,BSE,TDA,singlet_manifold,triplet_manifold, & 
                nBas,nC,nO,nV,nR,nS,ENuc,ERHF,Hc,ERI_AO_basis,ERI_MO_basis,P,cHF,eHF,eG0W0)

! Perform G0W0 calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: COHSEX,SOSEX,BSE,TDA,singlet_manifold,triplet_manifold
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ENuc,ERHF
  double precision,intent(in)   :: cHF(nBas,nBas),eHF(nBas),Hc(nBas,nBas),P(nBas,nBas)
  double precision,intent(in)   :: ERI_AO_basis(nBas,nBas,nBas,nBas),ERI_MO_basis(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dRPA
  integer                       :: ispin
  double precision              :: EcRPA,EcGM
  double precision,allocatable  :: H(:,:),SigmaC(:),Z(:)
  double precision,allocatable  :: Omega(:,:),XpY(:,:,:),rho(:,:,:,:),rhox(:,:,:,:)

! Output variables

  double precision              :: eG0W0(nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0W0 calculation           |'
  write(*,*)'************************************************'
  write(*,*)

! SOSEX correction

  if(SOSEX) write(*,*) 'SOSEX correction activated!'
  write(*,*)

! Switch off exchange for G0W0

  dRPA = .true.

! Spin manifold 

  ispin = 1

! Memory allocation

  allocate(H(nBas,nBas),SigmaC(nBas),Z(nBas), &
           Omega(nS,nspin),XpY(nS,nS,nspin),  & 
           rho(nBas,nBas,nS,nspin),rhox(nBas,nBas,nS,nspin))

! Compute Hartree Hamiltonian in the MO basis

  call Hartree_matrix_MO_basis(nBas,cHF,P,Hc,ERI_AO_basis,H)

! Compute linear response

  call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,eHF,ERI_MO_basis, & 
                       rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))

! Compute correlation part of the self-energy 

  call excitation_density_from_MO(nBas,nC,nO,nR,nS,ERI_MO_basis,XpY(:,:,ispin),rho(:,:,:,ispin))

  if(SOSEX) call excitation_density_SOSEX_from_MO(nBas,nC,nO,nR,nS,ERI_MO_basis,XpY(:,:,ispin),rhox(:,:,:,ispin))

  call self_energy_correlation_diag(COHSEX,SOSEX,nBas,nC,nO,nV,nR,nS,eHF, & 
                                    Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),EcGM,SigmaC)

! COHSEX static approximation

  if(COHSEX) then

    Z(:) = 1d0

  else

    call renormalization_factor(SOSEX,nBas,nC,nO,nV,nR,nS,eHF,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),Z)

  endif

! Solve the quasi-particle equation

  eG0W0(:) = eHF(:) + Z(:)*SigmaC(:)

! Dump results

  call print_excitation('RPA  ',ispin,nS,Omega(:,ispin))
  call print_G0W0(nBas,nO,eHF,ENuc,ERHF,SigmaC,Z,eG0W0,EcRPA,EcGM)

! Plot stuff

  call plot_GW(nBas,nC,nO,nV,nR,nS,eHF,eG0W0,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin))
 
! Perform BSE calculation

  if(BSE) then

   ! Singlet manifold

   if(singlet_manifold) then

      ispin = 1
      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,eG0W0,ERI_MO_basis, &
                           rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    endif

   ! Triplet manifold

   if(triplet_manifold) then

      ispin = 2
      call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,eHF,ERI_MO_basis, &
                           rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))
      call excitation_density(nBas,nC,nO,nR,nS,cHF,ERI_AO_basis,XpY(:,:,ispin),rho(:,:,:,ispin))

      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,eG0W0,ERI_MO_basis, &
                           rho(:,:,:,1),EcRPA,Omega(:,ispin),XpY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    endif

  endif


end subroutine G0W0
