subroutine qsGW(maxSCF,thresh,max_diis,COHSEX,BSE,TDA,G0W,GW0,singlet_manifold,triplet_manifold, &
                nBas,nC,nO,nV,nR,nS,ENuc,S,X,T,V,Hc,ERI_AO_basis,PHF,cHF,eHF)

! Compute linear response

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF,max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: COHSEX,BSE,TDA,G0W,GW0,singlet_manifold,triplet_manifold
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: PHF(nBas,nBas),cHF(nBas,nBas),eHF(nBas)
  double precision,intent(in)   :: S(nBas,nBas),T(nBas,nBAs),V(nBas,nBas),Hc(nBas,nBas),X(nBas,nBas)
  double precision,intent(in)   :: ERI_AO_basis(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dRPA
  integer                       :: nSCF,nBasSq,ispin,n_diis
  double precision              :: EcRPA,Conv
  double precision,external     :: trace_matrix
  double precision,allocatable  :: error_diis(:,:),F_diis(:,:)
  double precision,allocatable  :: Omega(:,:),XpY(:,:,:),rho(:,:,:,:)
  double precision,allocatable  :: c(:,:),cp(:,:),e(:),P(:,:)
  double precision,allocatable  :: F(:,:),Fp(:,:),J(:,:),K(:,:),SigmaC(:,:)
  double precision,allocatable  :: error(:,:),ERI_MO_basis(:,:,:,:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|       Self-consistent qsGW calculation       |'
  write(*,*)'************************************************'
  write(*,*)

! Stuff 

  nBasSq = nBas*nBas

! Switch off exchange for G0W0

  dRPA = .true.

! Memory allocation

  allocate(e(nBas),c(nBas,nBas),cp(nBas,nBas),P(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas), &
           J(nBas,nBas),K(nBas,nBas),SigmaC(nBas,nBas),              & 
           ERI_MO_basis(nBas,nBas,nBas,nBas),error(nBas,nBas),       &
           Omega(nS,nspin),XpY(nS,nS,nspin),rho(nBas,nBas,nS,nspin), &
           error_diis(nBasSq,max_diis),F_diis(nBasSq,max_diis))

! Initialization
  
  nSCF            = 0
  n_diis          = 0
  ispin           = 1
  Conv            = 1d0
  P(:,:)          = PHF(:,:)
  e(:)            = eHF(:)
  c(:,:)          = cHF(:,:)
  F_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

    ! Buid Coulomb matrix

    call Coulomb_matrix_AO_basis(nBas,P,ERI_AO_basis,J)

    ! Compute exchange part of the self-energy 

    call exchange_matrix_AO_basis(nBas,P,ERI_AO_basis,K)

    ! AO to MO transformation of two-electron integrals

    call AOtoMO_integral_transform(nBas,c,ERI_AO_basis,ERI_MO_basis)

    ! Compute linear response

    if(.not. GW0 .or. nSCF == 0) then

      call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,e,ERI_MO_basis, &
                           rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))

    endif

    ! Compute correlation part of the self-energy 

    call excitation_density(nBas,nC,nO,nR,nS,c,ERI_AO_basis,XpY(:,:,ispin),rho(:,:,:,ispin))

    if(G0W) then

      call self_energy_correlation(COHSEX,nBas,nC,nO,nV,nR,nS,eHF,Omega(:,ispin),rho(:,:,:,ispin),SigmaC)

     else

      call self_energy_correlation(COHSEX,nBas,nC,nO,nV,nR,nS,e,Omega(:,ispin),rho(:,:,:,ispin),SigmaC)

     endif

    ! Make correlation self-energy Hermitian and transform it back to AO basis
   
    SigmaC =  0.5d0*(SigmaC + transpose(SigmaC))

    call MOtoAO_transform(nBas,S,c,SigmaC)
 
    ! Solve the quasi-particle equation

    F(:,:) = Hc(:,:) + J(:,:) + K(:,:) + SigmaC(:,:)

    ! Compute commutator and convergence criteria

    error = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
    Conv = maxval(abs(error))

    ! DIIS extrapolation 

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(nBasSq,nBasSq,n_diis,error_diis,F_diis,error,F)

    ! Diagonalize Hamiltonian in AO basis

    Fp = matmul(transpose(X),matmul(F,X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp,e)
    c = matmul(X,cp)

    ! Compute new density matrix in the AO basis

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

    ! Print results

    call print_qsGW(nBas,nO,nSCF,Conv,thresh,eHF,e,c,ENuc,P,T,V,Hc,J,K,F,EcRPA)

    ! Increment

    nSCF = nSCF + 1

  enddo
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

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
      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI_MO_basis, & 
                           rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    endif

    ! Triplet manifold
    if(triplet_manifold) then

      ispin = 2
      call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,e,ERI_MO_basis, &
                             rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))
      call excitation_density(nBas,nC,nO,nR,nS,c,ERI_AO_basis,XpY(:,:,ispin),rho(:,:,:,ispin))
     
      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,e,ERI_MO_basis, &
                           rho(:,:,:,ispin),EcRPA,Omega(:,ispin),XpY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    endif

  endif

end subroutine qsGW
