subroutine evGW(maxSCF,thresh,max_diis,COHSEX,SOSEX,BSE,TDA,G0W,GW0,singlet_manifold,triplet_manifold,linearize, & 
                nBas,nC,nO,nV,nR,nS,ENuc,ERHF,Hc,ERI_AO_basis,ERI_MO_basis,PHF,cHF,eHF,eG0W0)

! Perform self-consistent eigenvalue-only GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF,max_diis
  double precision,intent(in)   :: thresh,ENuc,ERHF
  logical,intent(in)            :: COHSEX,SOSEX,BSE,TDA,G0W,GW0
  logical,intent(in)            :: singlet_manifold,triplet_manifold,linearize
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: cHF(nBas,nBas),eHF(nBas),eG0W0(nBas),Hc(nBas,nBas),PHF(nBas,nBas) 
  double precision,intent(in)   :: ERI_AO_basis(nBas,nBas,nBas,nBas),ERI_MO_basis(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dRPA,linear_mixing
  integer                       :: ispin,nSCF,n_diis
  double precision              :: rcond
  double precision              :: Conv,EcRPA,EcGM,lambda
  double precision,allocatable  :: error_diis(:,:),e_diis(:,:)
  double precision,allocatable  :: eGW(:),eOld(:),Z(:)
  double precision,allocatable  :: H(:,:),SigmaC(:)
  double precision,allocatable  :: Omega(:,:),XpY(:,:,:),rho(:,:,:,:),rhox(:,:,:,:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|       Self-consistent evGW calculation       |'
  write(*,*)'************************************************'
  write(*,*)

! SOSEX correction

  if(SOSEX) write(*,*) 'SOSEX correction activated!'
  write(*,*)

! Switch off exchange for G0W0

  dRPA = .true.

! Linear mixing

  linear_mixing = .false.
  lambda = 0.2d0

! Memory allocation

  allocate(eGW(nBas),eOld(nBas),Z(nBas),                       &
           H(nBas,nBas),SigmaC(nBas),                          &
           Omega(nS,nspin),XpY(nS,nS,nspin),                   &
           rho(nBas,nBas,nS,nspin),rhox(nBas,nBas,nS,nspin),   &
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

! Compute Hartree Hamiltonian in the MO basis

  call Hartree_matrix_MO_basis(nBas,cHF,PHF,Hc,ERI_AO_basis,H)

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

   ! Compute linear response

    if(.not. GW0 .or. nSCF == 0) then

      call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,eGW,ERI_MO_basis, & 
                           rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))

    endif

!   Compute correlation part of the self-energy 

    call excitation_density(nBas,nC,nO,nR,nS,cHF,ERI_AO_basis,XpY(:,:,ispin),rho(:,:,:,ispin))

    if(SOSEX) call excitation_density_SOSEX(nBas,nC,nO,nR,nS,cHF,ERI_AO_basis,XpY(:,:,ispin),rhox(:,:,:,ispin))

    ! Correlation self-energy

    if(G0W) then

      call self_energy_correlation_diag(COHSEX,SOSEX,nBas,nC,nO,nV,nR,nS,eHF, & 
                                        Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),EcGM,SigmaC)
      call renormalization_factor(SOSEX,nBas,nC,nO,nV,nR,nS,eHF,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),Z)

    else 

      call self_energy_correlation_diag(COHSEX,SOSEX,nBas,nC,nO,nV,nR,nS,eGW, & 
                                        Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),EcGM,SigmaC)
      call renormalization_factor(SOSEX,nBas,nC,nO,nV,nR,nS,eGW,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),Z)

    endif

    ! Solve the quasi-particle equation (linearized or not)

    if(linearize) then 

      eGW(:) = eHF(:) + Z(:)*SigmaC(:)

    else

      eGW(:) = eHF(:) + SigmaC(:)

    endif

    ! Convergence criteria

    Conv = maxval(abs(eGW - eOld))

    ! Print results

    call print_excitation('RPA  ',ispin,nS,Omega(:,ispin))
    call print_evGW(nBas,nO,nSCF,Conv,eHF,ENuc,ERHF,SigmaC,Z,eGW,EcRPA)

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
      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,eGW,ERI_MO_basis, &
                           rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    endif

    ! Triplet manifold
    if(triplet_manifold) then

      ispin = 2
      call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,eGW,ERI_MO_basis, &
                           rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))
      call excitation_density(nBas,nC,nO,nR,nS,cHF,ERI_AO_basis,XpY(:,:,ispin),rho(:,:,:,ispin))

      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,eGW,ERI_MO_basis, &
                           rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    endif

  endif

end subroutine evGW
