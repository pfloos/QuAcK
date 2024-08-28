
! ---

subroutine ROHF(dotest, maxSCF, thresh, max_diis, guess_type, mix, level_shift, nNuc, ZNuc, rNuc, ENuc, & 
                nBas_AOs, nBas_MOs, nO, S, T, V, Hc, ERI, dipole_int, X, EROHF, eHF, c, Ptot)

! Perform restricted open-shell Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: mix 
  double precision,intent(in)   :: level_shift
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: nBas_AOs, nBas_MOs

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nO(nspin)
  double precision,intent(in)   :: S(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: T(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: V(nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: Hc(nBas_AOs,nBas_AOs) 
  double precision,intent(in)   :: X(nBas_AOs,nBas_MOs) 
  double precision,intent(in)   :: ERI(nBas_AOs,nBas_AOs,nBas_AOs,nBas_AOs)
  double precision,intent(in)   :: dipole_int(nBas_AOs,nBas_AOs,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: nBas_AOs_Sq
  integer                       :: n_diis
  double precision              :: Conv
  double precision              :: rcond
  double precision              :: ET(nspin)
  double precision              :: EV(nspin)
  double precision              :: EJ(nsp)
  double precision              :: EK(nspin)
  double precision              :: dipole(ncart)

  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: J(:,:,:)
  double precision,allocatable  :: F(:,:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: Ftot(:,:)
  double precision,allocatable  :: P(:,:,:)
  double precision,allocatable  :: K(:,:,:)
  double precision,allocatable  :: err(:,:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,external     :: trace_matrix

  integer                       :: ispin

! Output variables

  double precision,intent(out)  :: EROHF
  double precision,intent(out)  :: eHF(nBas_MOs)
  double precision,intent(inout):: c(nBas_AOs,nBas_MOs)
  double precision,intent(out)  :: Ptot(nBas_AOs,nBas_AOs)

! Hello world

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'* Restricted Open-Shell HF Calculation *'
  write(*,*)'****************************************'
  write(*,*)

! Useful stuff

  nBas_AOs_Sq = nBas_AOs*nBas_AOs

! Memory allocation

  allocate(J(nBas_AOs,nBas_AOs,nspin))
  allocate(K(nBas_AOs,nBas_AOs,nspin))
  allocate(F(nBas_AOs,nBas_AOs,nspin))
  allocate(Ftot(nBas_AOs,nBas_AOs))
  allocate(P(nBas_AOs,nBas_AOs,nspin))
  allocate(err(nBas_AOs,nBas_AOs))

  allocate(Fp(nBas_MOs,nBas_MOs))
  allocate(cp(nBas_MOs,nBas_MOs))

  allocate(err_diis(nBas_AOs_Sq,max_diis))
  allocate(F_diis(nBas_AOs_Sq,max_diis))

! Guess coefficients and demsity matrices

  call mo_guess(nBas_AOs, nBas_MOs, guess_type, S, Hc, X, c)

  do ispin = 1, nspin
    !P(:,:,ispin) = matmul(c(:,1:nO(ispin)), transpose(c(:,1:nO(ispin))))
    call dgemm('N', 'T', nBas_AOs, nBas_AOs, nO(ispin), 1.d0, c, nBas_AOs, c, nBas_AOs, 0.d0, P(1,1,ispin), nBas_AOs)
  end do
  Ptot(:,:) = P(:,:,1) + P(:,:,2)

! Initialization

  n_diis        = 0
  F_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0
  rcond         = 0d0

  nSCF = 0
  Conv = 1d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'-----------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(ROHF)','|','EJ(ROHF)','|','EK(ROHF)','|','Conv','|'
  write(*,*)'-----------------------------------------------------------------------------'
  
  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Build Hartree repulsion

    do ispin = 1, nspin
      call Hartree_matrix_AO_basis(nBas_AOs, P(:,:,ispin), ERI(:,:,:,:), J(:,:,ispin))
    end do

!   Compute exchange potential

    do ispin = 1, nspin
      call exchange_matrix_AO_basis(nBas_AOs, P(:,:,ispin), ERI(:,:,:,:), K(:,:,ispin))
    end do
 
!   Build Fock operator

    do ispin = 1, nspin
      F(:,:,ispin) = Hc(:,:) + J(:,:,ispin) + J(:,:,mod(ispin,2)+1) + K(:,:,ispin)
    end do

    call ROHF_fock_matrix(nBas_AOs, nBas_MOs, nO(1), nO(2), S, c, F(:,:,1), F(:,:,2), Ftot)

!   Check convergence 

    err(:,:) = matmul(Ftot, matmul(Ptot, S)) - matmul(matmul(S, Ptot), Ftot)
    if(nSCF > 1) Conv = maxval(abs(err(:,:)))
    
!   Kinetic energy

    do ispin = 1, nspin
      ET(ispin) = trace_matrix(nBas_AOs, matmul(P(:,:,ispin), T(:,:)))
    end do

!   Potential energy

    do ispin = 1, nspin
      EV(ispin) = trace_matrix(nBas_AOs, matmul(P(:,:,ispin), V(:,:)))
    end do

!   Hartree energy

    EJ(1) = 0.5d0*trace_matrix(nBas_AOs, matmul(P(:,:,1), J(:,:,1)))
    EJ(2) = trace_matrix(nBas_AOs, matmul(P(:,:,1), J(:,:,2)))
    EJ(3) = 0.5d0*trace_matrix(nBas_AOs, matmul(P(:,:,2), J(:,:,2)))

!   Exchange energy

    do ispin = 1, nspin
      EK(ispin) = 0.5d0*trace_matrix(nBas_AOs, matmul(P(:,:,ispin), K(:,:,ispin)))
    end do

!   Total energy

    EROHF = sum(ET) + sum(EV) + sum(EJ) + sum(EK)

!   DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBas_AOs_Sq,nBas_AOs_Sq,n_diis,err_diis,F_diis,err,Ftot)

    end if

!   Level-shifting

    if(level_shift > 0d0 .and. Conv > thresh) then

      do ispin=1,nspin
        call level_shifting(level_shift, nBas_AOs, nBas_MOs, maxval(nO), S, c, Ftot)
      end do

    end if

!  Transform Fock matrix in orthogonal basis

    Fp(:,:) = matmul(transpose(X(:,:)), matmul(Ftot(:,:), X(:,:)))

!  Diagonalize Fock matrix to get eigenvectors and eigenvalues

    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas_MOs, cp, eHF)
    
!   Back-transform eigenvectors in non-orthogonal basis

    c(:,:) = matmul(X(:,:), cp(:,:))

!   Compute density matrix 

    do ispin = 1, nspin
      !P(:,:,ispin) = matmul(c(:,1:nO(ispin)), transpose(c(:,1:nO(ispin))))
      call dgemm('N', 'T', nBas_AOs, nBas_AOs, nO(ispin), 1.d0, c, nBas_AOs, c, nBas_AOs, 0.d0, P(1,1,ispin), nBas_AOs)
    end do
    Ptot(:,:) = P(:,:,1) + P(:,:,2) 

!   Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,E10.2,1X,A1,1X)') &
      '|',nSCF,'|',EROHF + ENuc,'|',sum(EJ),'|',sum(EK),'|',Conv,'|'
 
  end do
  write(*,*)'-----------------------------------------------------------------------------'
!------------------------------------------------------------------------
! End of SCF loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    deallocate(J, K, F, Ftot, P, err, Fp, cp, err_diis, F_diis)

    stop

  end if

! Compute final UHF energy

  call dipole_moment(nBas_AOs,Ptot,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_ROHF(nBas_AOs, nBas_MOs, nO, eHF, c, ENuc, ET, EV, EJ, EK, EROHF, dipole)

! Print test values

  if(dotest) then

    call dump_test_value('R','ROHF energy',EROHF)

  end if

  deallocate(J, K, F, Ftot, P, err, Fp, cp, err_diis, F_diis)

end subroutine 
