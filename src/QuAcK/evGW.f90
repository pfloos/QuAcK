subroutine evGW(maxSCF,thresh,max_diis,COHSEX,SOSEX,BSE,TDA,G0W,GW0,singlet_manifold,triplet_manifold,linearize, & 
                nBas,nC,nO,nV,nR,nS,ENuc,ERHF,Hc,H,ERI_MO_basis,PHF,cHF,eHF,eG0W0)

! Perform self-consistent eigenvalue-only GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: SOSEX
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: G0W
  logical,intent(in)            :: GW0
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  logical,intent(in)            :: linearize
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: PHF(nBas,nBas) 
  double precision,intent(in)   :: eG0W0(nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: H(nBas,nBas)
  double precision,intent(in)   :: ERI_MO_basis(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dRPA
  logical                       :: linear_mixing
  integer                       :: ispin
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: rcond
  double precision              :: Conv
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcGM
  double precision              :: lambda
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)
  double precision,allocatable  :: rho(:,:,:,:)
  double precision,allocatable  :: rhox(:,:,:,:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|       Self-consistent evGW calculation       |'
  write(*,*)'************************************************'
  write(*,*)

! SOSEX correction

  if(SOSEX) write(*,*) 'SOSEX correction activated!'
  write(*,*)

! COHSEX approximation

  if(COHSEX) write(*,*) 'COHSEX approximation activated!'
  write(*,*)

! Switch off exchange for G0W0

  dRPA = .true.

! Linear mixing

  linear_mixing = .false.
  lambda = 0.2d0

! Memory allocation

  allocate(eGW(nBas),eOld(nBas),Z(nBas),SigC(nBas),Omega(nS,nspin),           & 
           XpY(nS,nS,nspin),rho(nBas,nBas,nS,nspin),rhox(nBas,nBas,nS,nspin), &
           error_diis(nBas,max_diis),e_diis(nBas,max_diis))

! Initialization

  nSCF            = 0
  ispin           = 1
  n_diis          = 0
  Conv            = 1d0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGW(:)          = eG0W0(:)
  eOld(:)         = eGW(:)
  Z(:)            = 1d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

   ! Compute linear response

    if(.not. GW0 .or. nSCF == 0) then

      call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,eGW,ERI_MO_basis, & 
                           rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin))

    endif

!   Compute correlation part of the self-energy 

    call excitation_density(nBas,nC,nO,nR,nS,ERI_MO_basis,XpY(:,:,ispin),rho(:,:,:,ispin))

    if(SOSEX) call excitation_density_SOSEX(nBas,nC,nO,nR,nS,ERI_MO_basis,XpY(:,:,ispin),rhox(:,:,:,ispin))

    ! Correlation self-energy

    if(G0W) then

      call self_energy_correlation_diag(COHSEX,SOSEX,nBas,nC,nO,nV,nR,nS,eHF, & 
                                        Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),EcGM,SigC)
      call renormalization_factor(SOSEX,nBas,nC,nO,nV,nR,nS,eHF,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),Z)

    else 

      call self_energy_correlation_diag(COHSEX,SOSEX,nBas,nC,nO,nV,nR,nS,eGW, & 
                                        Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),EcGM,SigC)
      call renormalization_factor(SOSEX,nBas,nC,nO,nV,nR,nS,eGW,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),Z)

    endif

    ! Solve the quasi-particle equation (linearized or not)

    if(linearize) then 

      eGW(:) = eHF(:) + Z(:)*SigC(:)

    else

      eGW(:) = eHF(:) + SigC(:)

    endif

    ! Convergence criteria

    Conv = maxval(abs(eGW - eOld))

    ! Print results

    call print_excitation('RPA  ',ispin,nS,Omega(:,ispin))
    call print_evGW(nBas,nO,nSCF,Conv,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA(ispin),EcGM)

    ! Linear mixing or DIIS extrapolation

    if(linear_mixing) then
 
      eGW(:) = lambda*eGW(:) + (1d0 - lambda)*eOld(:)
 
    else

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBas,nBas,n_diis,error_diis,e_diis,eGW-eOld,eGW)

!    Reset DIIS if required

      if(abs(rcond) < 1d-15) n_diis = 0

    endif

    ! Save quasiparticles energy for next cycle

    eOld(:) = eGW(:)

    ! Increment

    nSCF = nSCF + 1

  enddo
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Plot stuff

  call plot_GW(nBas,nC,nO,nV,nR,nS,eHF,eGW,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin))

! Did it actually converge?

  if(nSCF == maxSCF+1) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    if(BSE) stop

  endif

! Perform BSE calculation

  if(BSE) then

    ! Singlet manifold
    if(singlet_manifold) then

      ispin = 1
      EcBSE(ispin) = 0d0

      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,eGW,ERI_MO_basis, &
                           rho(:,:,:,ispin),EcBSE(ispin),Omega(:,ispin),XpY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    endif

    ! Triplet manifold
    if(triplet_manifold) then

      ispin = 2
      EcBSE(ispin) = 0d0

      call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,eGW,ERI_MO_basis, &
                           rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin))
      call excitation_density(nBas,nC,nO,nR,nS,ERI_MO_basis,XpY(:,:,ispin),rho(:,:,:,ispin))

      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,eGW,ERI_MO_basis, &
                           rho(:,:,:,ispin),EcBSE(ispin),Omega(:,ispin),XpY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    endif

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A40,F15.6)') 'BSE@evGW correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A40,F15.6)') 'BSE@evGW correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A40,F15.6)') 'BSE@evGW correlation energy           =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A40,F15.6)') 'BSE@evGW total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  endif

end subroutine evGW
