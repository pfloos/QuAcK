
! ---

subroutine RHF(dotest, maxSCF, thresh, max_diis, guess_type, level_shift, nNuc, ZNuc, rNuc, ENuc, & 
               nBas_AOs, nBas_MOs, nO, S, T, V, Hc, ERI, dipole_int, X, ERHF, eHF, c, P)

! Perform restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift

  integer,intent(in)            :: nBas_AOs, nBas_MOs
  integer,intent(in)            :: nO
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
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
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: EK
  double precision              :: dipole(ncart)

  double precision              :: Conv
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision,allocatable  :: err(:,:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Fp(:,:)

! Output variables

  double precision,intent(out)  :: ERHF
  double precision,intent(out)  :: eHF(nBas_MOs)
  double precision,intent(inout):: c(nBas_AOs,nBas_MOs)
  double precision,intent(out)  :: P(nBas_AOs,nBas_AOs)

! Hello world

  write(*,*)
  write(*,*)'*****************************'
  write(*,*)'* Restricted HF Calculation *'
  write(*,*)'*****************************'
  write(*,*)

! Useful quantities

  nBas_AOs_Sq = nBas_AOs*nBas_AOs

! Memory allocation

  allocate(J(nBas_AOs,nBas_AOs))
  allocate(K(nBas_AOs,nBas_AOs))

  allocate(err(nBas_AOs,nBas_AOs))
  allocate(F(nBas_AOs,nBas_AOs))

  allocate(cp(nBas_MOs,nBas_MOs))
  allocate(Fp(nBas_MOs,nBas_MOs))

  allocate(err_diis(nBas_AOs_Sq,max_diis))
  allocate(F_diis(nBas_AOs_Sq,max_diis))

! Guess coefficients and density matrix

  call mo_guess(nBas_AOs, nBas_MOs, guess_type, S, Hc, X, c)

  !P(:,:) = 2d0 * matmul(c(:,1:nO), transpose(c(:,1:nO)))
  call dgemm('N', 'T', nBas_AOs, nBas_AOs, nO, 2.d0, c, nBas_AOs, c, nBas_AOs, 0.d0, P, nBas_AOs)

! Initialization

  n_diis        = 0
  F_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0
  rcond = 0d0

  Conv   = 1d0
  nSCF   = 0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'-----------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(RHF)','|','EJ(RHF)','|','EK(RHF)','|','Conv','|'
  write(*,*)'-----------------------------------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Increment 

    nSCF = nSCF + 1

    ! Build Fock matrix
    
    call Hartree_matrix_AO_basis(nBas_AOs, P, ERI, J)
    call exchange_matrix_AO_basis(nBas_AOs, P, ERI, K)
    
    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)

    ! Check convergence 

    err = matmul(F, matmul(P, S)) - matmul(matmul(S, P), F)
    if(nSCF > 1) Conv = maxval(abs(err))

    ! Kinetic energy

    ET = trace_matrix(nBas_AOs, matmul(P, T))

    ! Potential energy

    EV = trace_matrix(nBas_AOs, matmul(P, V))

    ! Hartree energy

    EJ = 0.5d0*trace_matrix(nBas_AOs, matmul(P, J))

    ! Exchange energy

    EK = 0.25d0*trace_matrix(nBas_AOs, matmul(P, K))

    ! Total energy

    ERHF = ET + EV + EJ + EK

    ! DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1, max_diis)
      call DIIS_extrapolation(rcond, nBas_AOs_Sq, nBas_AOs_Sq, n_diis, err_diis, F_diis, err, F)

    end if

    ! Level shift

    if(level_shift > 0d0 .and. Conv > thresh) then
      call level_shifting(level_shift, nBas_AOs, nBas_MOs, nO, S, c, F)
    endif

    ! Diagonalize Fock matrix

    Fp = matmul(transpose(X), matmul(F, X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas_MOs, cp, eHF)
    c = matmul(X, cp)

    ! Density matrix

    !P(:,:) = 2d0*matmul(c(:,1:nO), transpose(c(:,1:nO)))
    call dgemm('N', 'T', nBas_AOs, nBas_AOs, nO, 2.d0, c, nBas_AOs, c, nBas_AOs, 0.d0, P, nBas_AOs)

    ! Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,E10.2,1X,A1,1X)') &
      '|',nSCF,'|',ERHF + ENuc,'|',EJ,'|',EK,'|',Conv,'|'

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

    deallocate(J, K, err, F, cp, Fp, err_diis, F_diis)

    stop

  end if

! Compute dipole moments

  call dipole_moment(nBas_AOs, P, nNuc, ZNuc, rNuc, dipole_int, dipole)
  call print_RHF(nBas_AOs, nBas_MOs, nO, eHF, c, ENuc, ET, EV, EJ, EK, ERHF, dipole)

! Testing zone

  if(dotest) then
 
    call dump_test_value('R','RHF energy',ERHF)
    call dump_test_value('R','RHF HOMO energy',eHF(nO))
    call dump_test_value('R','RHF LUMO energy',eHF(nO+1))
    call dump_test_value('R','RHF dipole moment',norm2(dipole))

  end if

  deallocate(J, K, err, F, cp, Fp, err_diis, F_diis)

end subroutine 
