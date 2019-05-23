subroutine G0W0(COHSEX,SOSEX,BSE,TDA,singlet_manifold,triplet_manifold, & 
                nBas,nC,nO,nV,nR,nS,ENuc,ERHF,Hc,H,ERI,PHF,cHF,eHF,eG0W0)

! Perform G0W0 calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: SOSEX
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
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
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)
  double precision,allocatable  :: rho(:,:,:,:)
  double precision,allocatable  :: rhox(:,:,:,:)

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

  allocate(SigC(nBas),Z(nBas),Omega(nS,nspin),XpY(nS,nS,nspin),  & 
           rho(nBas,nBas,nS,nspin),rhox(nBas,nBas,nS,nspin))

! Compute linear response

  call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,eHF,ERI, & 
                       rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin))

! Compute correlation part of the self-energy 

  call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

  if(SOSEX) call excitation_density_SOSEX(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rhox(:,:,:,ispin))

  call self_energy_correlation_diag(COHSEX,SOSEX,nBas,nC,nO,nV,nR,nS,eHF, & 
                                    Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),EcGM,SigC)

! COHSEX static approximation

  if(COHSEX) then

    Z(:) = 1d0

  else

    call renormalization_factor(SOSEX,nBas,nC,nO,nV,nR,nS,eHF,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),Z)

  endif

! Solve the quasi-particle equation

  eG0W0(:) = eHF(:) + Z(:)*SigC(:)

! Dump results

  call print_excitation('RPA  ',ispin,nS,Omega(:,ispin))
  call print_G0W0(nBas,nO,eHF,ENuc,ERHF,SigC,Z,eG0W0,EcRPA(ispin),EcGM)

! Plot stuff

  call plot_GW(nBas,nC,nO,nV,nR,nS,eHF,eG0W0,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin))
 
! Perform BSE calculation

  if(BSE) then

   ! Singlet manifold

   if(singlet_manifold) then

      ispin = 1
      EcBSE(ispin) = 0d0

      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,eG0W0,ERI, &
                           rho(:,:,:,ispin),EcBSE(ispin),Omega(:,ispin),XpY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    endif

   ! Triplet manifold

   if(triplet_manifold) then

      ispin = 2
      EcBSE(ispin) = 0d0

      call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,eHF,ERI, &
                           rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin))
      call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,eG0W0,ERI, &
                           rho(:,:,:,1),EcBSE(ispin),Omega(:,ispin),XpY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    endif

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A40,F15.6)') 'BSE@G0W0 correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A40,F15.6)') 'BSE@G0W0 correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A40,F15.6)') 'BSE@G0W0 correlation energy           =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A40,F15.6)') 'BSE@G0W0 total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  endif


end subroutine G0W0
