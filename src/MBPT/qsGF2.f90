subroutine qsGF2(maxSCF,thresh,max_diis,BSE,TDA,dBSE,dTDA,evDyn,singlet,triplet, & 
                eta,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,ERHF, &
                S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)

! Perform a quasiparticle self-consistent GF2 calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_MO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_MO(nBas,nBas,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: ispin
  integer                       :: n_diis
  double precision              :: EqsGF2
  double precision              :: Conv
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision              :: dipole(ncart)
  double precision              :: Ec
  double precision              :: EcBSE(nspin)

  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: eGF2(:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: SigCp(:,:)
  double precision,allocatable  :: SigCm(:,:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: error(:,:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|       Self-consistent qsGF2 calculation      |'
  write(*,*)'************************************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsGF2 !!'
  write(*,*)

! Stuff 

  nBasSq = nBas*nBas

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  allocate(eGF2(nBas),c(nBas,nBas),cp(nBas,nBas),P(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas),        &
           J(nBas,nBas),K(nBas,nBas),SigC(nBas,nBas),SigCp(nBas,nBas),SigCm(nBas,nBas),Z(nBas), & 
           error(nBas,nBas),error_diis(nBasSq,max_diis),F_diis(nBasSq,max_diis))

! Initialization
  
  nSCF            = -1
  n_diis          = 0
  ispin           = 1
  Conv            = 1d0
  P(:,:)          = PHF(:,:)
  eGF2(:)         = eHF(:)
  c(:,:)          = cHF(:,:)
  F_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

    ! Increment

    nSCF = nSCF + 1

    ! Buid Coulomb matrix

    call Coulomb_matrix_AO_basis(nBas,P,ERI_AO,J)

    ! Compute exchange part of the self-energy 

    call exchange_matrix_AO_basis(nBas,P,ERI_AO,K)

    ! AO to MO transformation of two-electron integrals

    call AOtoMO_integral_transform(1,1,1,1,nBas,c,ERI_AO,ERI_MO)

    ! Compute self-energy and renormalization factor

    call self_energy_GF2(eta,nBas,nC,nO,nV,nR,nS,eHF,eGF2,ERI_MO,SigC,Z,Ec)

    ! Make correlation self-energy Hermitian and transform it back to AO basis
   
    SigCp = 0.5d0*(SigC + transpose(SigC))
    SigCm = 0.5d0*(SigC - transpose(SigC))

    call MOtoAO_transform(nBas,S,c,SigCp)
 
    ! Solve the quasi-particle equation

    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + SigCp(:,:)

    ! Compute commutator and convergence criteria

    error = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
    Conv = maxval(abs(error))

    ! DIIS extrapolation 

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nBasSq,nBasSq,n_diis,error_diis,F_diis,error,F)

!   Reset DIIS if required

    if(abs(rcond) < 1d-15) n_diis = 0

    ! Diagonalize Hamiltonian in AO basis

    Fp = matmul(transpose(X),matmul(F,X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp,eGF2)
    c = matmul(X,cp)

    ! Compute new density matrix in the AO basis

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

    ! Print results

    call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsGF2(nBas,nO,nSCF,Conv,thresh,eHF,eGF2,c,ENuc,P,T,V,J,K,F,SigCp,Z,EqsGF2,Ec,dipole)

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

    stop

  endif

! Deallocate memory

  deallocate(c,cp,P,F,Fp,J,K,SigC,SigCp,SigCm,Z,error,error_diis,F_diis)

! Perform BSE calculation

  if(BSE) then

    call BSE2(TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO,eHF,eGF2,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGF2 correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGF2 correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGF2 correlation energy           =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGF2 total energy                 =',ENuc + EqsGF2 + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine qsGF2
