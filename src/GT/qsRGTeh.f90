
! ---

subroutine qsRGTeh(dotest, maxSCF, thresh, max_diis, doACFDT, exchange_kernel, doXBS, dophBSE, &
                   dophBSE2, TDA_T, TDA, dBSE, dTDA, singlet, triplet, eta, regularize, nNuc,  &
                   ZNuc, rNuc, ENuc, nBas_AOs, nBas_MOs, nC, nO, nV, nR, nS, ERHF, S, X, T, V, &
                   Hc, ERI_AO, ERI_MO, dipole_int_AO, dipole_int_MO, PHF, cHF, eHF)

! Perform a quasiparticle self-consistent GTeh calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas_AOs, nBas_MOs
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas_MOs)
  double precision,intent(in)   :: cHF(nBas_AOs,nBas_MOs)
  double precision,intent(in)   :: PHF(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: S(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: T(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: V(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: Hc(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: X(nBas_AOs,nBas_MOs)
  double precision,intent(in)   :: ERI_AO(nBas_AOs,nBas_AOs,nBas_AOs,nBas_AOs)
  double precision,intent(inout):: ERI_MO(nBas_MOs,nBas_MOs,nBas_MOs,nBas_MOs)
  double precision,intent(in)   :: dipole_int_AO(nBas_AOs,nBas_AOs,ncart)
  double precision,intent(in)   :: dipole_int_MO(nBas_MOs,nBas_MOs,ncart)

! Local variables

  logical                       :: dRPA = .false.
  integer                       :: nSCF
  integer                       :: nBas_AOs_Sq
  integer                       :: ispin
  integer                       :: n_diis
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: Ex
  double precision              :: EqsGT
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM
  double precision              :: Conv
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision              :: dipole(ncart)

  logical                       :: print_T = .true.
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rhoL(:,:,:)
  double precision,allocatable  :: rhoR(:,:,:)
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: Sig(:,:)
  double precision,allocatable  :: Sigp(:,:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: err(:,:)

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Restricted qsGTeh Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsGTeh !!'
  write(*,*)

! Stuff 

  nBas_AOs_Sq = nBas_AOs*nBas_AOs

! TDA for T

  if(TDA_T) then 
    write(*,*) 'Tamm-Dancoff approximation for eh T-matrix!'
    write(*,*)
  end if

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  allocate(Aph(nS,nS), Bph(nS,nS), Om(nS), XpY(nS,nS), XmY(nS,nS))

  allocate(eGT(nBas_MOs))
  allocate(eOld(nBas_MOs))
  allocate(Z(nBas_MOs))

  allocate(c(nBas_AOs,nBas_MOs))

  allocate(cp(nBas_MOs,nBas_MOs))
  allocate(Fp(nBas_MOs,nBas_MOs))
  allocate(Sig(nBas_MOs,nBas_MOs))

  allocate(P(nBas_AOs,nBas_AOs))
  allocate(F(nBas_AOs,nBas_AOs))
  allocate(J(nBas_AOs,nBas_AOs))
  allocate(K(nBas_AOs,nBas_AOs))
  allocate(Sigp(nBas_AOs,nBas_AOs))
  allocate(err(nBas_AOs,nBas_AOs))

  allocate(err_diis(nBas_AOs_Sq,max_diis), F_diis(nBas_AOs_Sq,max_diis))

  allocate(rhoL(nBas_MOs,nBas_MOs,nS), rhoR(nBas_MOs,nBas_MOs,nS))


! Initialization
  
  nSCF          = -1
  n_diis        = 0
  ispin         = 2
  Conv          = 1d0
  P(:,:)        = PHF(:,:)
  eGT(:)        = eHF(:)
  eOld(:)       = eHF(:)
  c(:,:)        = cHF(:,:)
  F_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0
  rcond         = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

    ! Increment

    nSCF = nSCF + 1

    ! Buid Hartree matrix

    call Hartree_matrix_AO_basis(nBas_AOs,P,ERI_AO,J)

    ! Compute exchange part of the self-energy 

    call exchange_matrix_AO_basis(nBas_AOs,P,ERI_AO,K)

    ! AO to MO transformation of two-electron integrals

    call AOtoMO_ERI_RHF(nBas_AOs, nBas_MOs, c, ERI_AO, ERI_MO)

    ! Compute linear response

    call phLR_A(ispin,dRPA,nBas_MOs,nC,nO,nV,nR,nS,1d0,eGT,ERI_MO,Aph)
    if(.not.TDA_T) call phLR_B(ispin,dRPA,nBas_MOs,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)

    call phLR(TDA_T,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

    if(print_T) call print_excitation_energies('phRPA@RHF','triplet',nS,Om)

    ! Compute correlation part of the self-energy 

    call GTeh_excitation_density(nBas_MOs,nC,nO,nR,nS,ERI_MO,XpY,XmY,rhoL,rhoR)

    if(regularize) call GTeh_regularization(nBas_MOs,nC,nO,nV,nR,nS,eGT,Om,rhoL,rhoR)

    call GTeh_self_energy(eta,nBas_MOs,nC,nO,nV,nR,nS,eGT,Om,rhoL,rhoR,EcGM,Sig,Z)

    ! Make correlation self-energy Hermitian and transform it back to AO basis
   
    Sig = 0.5d0*(Sig + transpose(Sig))

    call MOtoAO(nBas_AOs, nBas_MOs, S, c, Sig, Sigp)
 
    ! Solve the quasi-particle equation

    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + Sigp(:,:)

    ! Compute commutator and convergence criteria

    err = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)

    ! DIIS extrapolation 

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBas_AOs_Sq,nBas_AOs_Sq,n_diis,err_diis,F_diis,err,F)

    end if

    ! Diagonalize Hamiltonian in AO basis

    Fp = matmul(transpose(X),matmul(F,X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas_MOs, cp, eGT)
    c = matmul(X,cp)

    ! Compute new density matrix in the AO basis

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

    ! Save quasiparticles energy for next cycle

    Conv = maxval(abs(err))
    eOld(:) = eGT(:)

    !------------------------------------------------------------------------
    !   Compute total energy
    !------------------------------------------------------------------------

    ! Kinetic energy

    ET = trace_matrix(nBas_AOs,matmul(P,T))

    ! Potential energy

    EV = trace_matrix(nBas_AOs,matmul(P,V))

    ! Hartree energy

    EJ = 0.5d0*trace_matrix(nBas_AOs,matmul(P,J))

    ! Exchange energy

    Ex = 0.25d0*trace_matrix(nBas_AOs,matmul(P,K))

    ! Total energy

    EqsGT = ET + EV + EJ + Ex 

    ! Print results

    call dipole_moment(nBas_AOs,P,nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsRGTeh(nBas_AOs, nBas_MOs, nO, nSCF, Conv, thresh, eHF, eGT, c, Sig, &
                       Z, ENuc, ET, EV, EJ, Ex, EcGM, EcRPA, EqsGT, dipole)

  end do
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

    deallocate(c, cp, P, F, Fp, J, K, Sig, Sigp, Z, Om, XpY, XmY, rhoL, rhoR, err, err_diis, F_diis)

    stop

  end if

! Deallocate memory

  deallocate(c, cp, P, F, Fp, J, K, Sig, Sigp, Z, Om, XpY, XmY, rhoL, rhoR, err, err_diis, F_diis)

! Perform BSE calculation

! if(BSE) then

!   call Bethe_Salpeter(BSE2,TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nBas_AOs,nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO, & 
!                       eGT,eGT,EcBSE)

!   if(exchange_kernel) then

!     EcBSE(1) = 0.5d0*EcBSE(1)
!     EcBSE(2) = 1.5d0*EcBSE(2)

!   end if

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW correlation energy (singlet) =',EcBSE(1)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW correlation energy (triplet) =',EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW correlation energy           =',EcBSE(1) + EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW total energy                 =',ENuc + EqsGW + EcBSE(1) + EcBSE(2)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

!   if(doACFDT) then

!     write(*,*) '------------------------------------------------------'
!     write(*,*) 'Adiabatic connection version of BSE correlation energy'
!     write(*,*) '------------------------------------------------------'
!     write(*,*)

!     if(doXBS) then

!       write(*,*) '*** scaled screening version (XBS) ***'
!       write(*,*)

!     end if

!     call ACFDT(exchange_kernel,doXBS,.true.,TDA_T,TDA,BSE,singlet,triplet,eta,nBas_AOs,nC,nO,nV,nR,nS,ERI_MO,eGW,eGW,EcAC)

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy (singlet) =',EcAC(1)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy (triplet) =',EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy           =',EcAC(1) + EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW total energy                 =',ENuc + EqsGW + EcAC(1) + EcAC(2)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

! end if

! Testing zone

  if(dotest) then

    call dump_test_value('R','qsGTeh correlation energy',EcRPA)
    call dump_test_value('R','qsGTeh HOMO energy',eGT(nO))
    call dump_test_value('R','qsGTeh LUMO energy',eGT(nO+1))

  end if

end subroutine
