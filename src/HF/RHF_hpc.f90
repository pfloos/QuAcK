subroutine RHF_hpc(working_dir,dotest,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
                    nBas,nOrb,nO,S,T,V,Hc,dipole_int,X,ERHF,eHF,c,P,F)

! Perform restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  character(len=256),intent(in) :: working_dir

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ii, jj
  integer                       :: nSCF
  integer                       :: nBas_Sq
  integer                       :: n_diis
  integer*8                     :: ERI_size
  double precision              :: t1, t2
  double precision              :: diff, diff_loc
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
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: ERI_chem(:)
  double precision,allocatable  :: ERI_phys(:,:,:,:), J_deb(:,:), K_deb(:,:)
  double precision,allocatable  :: tmp1(:,:), FX(:,:)


! Output variables

  double precision,intent(out)  :: ERHF
  double precision,intent(out)  :: eHF(nOrb)
  double precision,intent(inout):: c(nBas,nOrb)
  double precision,intent(out)  :: P(nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'****************************************'
  write(*,*)'* Restricted HF Calculation (HPC mode) *'
  write(*,*)'****************************************'
  write(*,*)

! Useful quantities

  nBas_Sq = nBas*nBas

! Memory allocation

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))

  allocate(err(nBas,nBas))

  allocate(cp(nOrb,nOrb))
  allocate(Fp(nOrb,nOrb))

  allocate(err_diis(nBas_Sq,max_diis))
  allocate(F_diis(nBas_Sq,max_diis))

  allocate(tmp1(nBas,nBas))
  allocate(FX(nBas,nOrb))

! Guess coefficients and density matrix
  call mo_guess(nBas,nOrb,guess_type,S,Hc,X,c)

  call dgemm('N', 'T', nBas, nBas, nO, 2.d0, &
             c(1,1), nBas, c(1,1), nBas,     &
             0.d0, P(1,1), nBas)

  ERI_size = shiftr(nBas * (nBas + 1), 1)
  ERI_size = shiftr(ERI_size * (ERI_size + 1), 1)
  allocate(ERI_chem(ERI_size))
  call read_2e_integrals_hpc(working_dir, ERI_size, ERI_chem)

  !call wall_time(t1)
  !call Hartree_matrix_AO_basis_hpc(nBas, ERI_size, P, ERI_chem, J)
  !call wall_time(t2)
  !print*, " J built in (sec):", (t2-t1)
  !call wall_time(t1)
  !call exchange_matrix_AO_basis_hpc(nBas, ERI_size, P, ERI_chem, K)
  !call wall_time(t2)
  !print*, " K built in (sec):", (t2-t1)
  !allocate(ERI_phys(nBas,nBas,nBas,nBas))
  !allocate(J_deb(nBas,nBas))
  !allocate(K_deb(nBas,nBas))
  !call read_2e_integrals(working_dir, nBas, ERI_phys)
  !call wall_time(t1)
  !call Hartree_matrix_AO_basis(nBas, P, ERI_phys, J_deb)
  !call wall_time(t2)
  !print*, " J_deb built in (sec):", (t2-t1)
  !call wall_time(t1)
  !call exchange_matrix_AO_basis(nBas, P, ERI_phys, K_deb)
  !call wall_time(t2)
  !print*, " K_deb built in (sec):", (t2-t1)
  !print*, "max error on J = ", maxval(dabs(J - J_deb))
  !diff = 0.d0
  !do ii = 1, nBas
  !  do jj = 1, nBas
  !    diff_loc = dabs(J(jj,ii) - J_deb(jj,ii))
  !    if(diff_loc .gt. 1d-10) then
  !      print*, 'error in J on: ', jj, ii
  !      print*, J(jj,ii), J_deb(jj,ii)
  !      stop
  !    endif
  !    diff = diff + diff_loc
  !  enddo
  !enddo
  !print*, 'total diff on J = ', diff
  !print*, "max error on K = ", maxval(dabs(K - K_deb))
  !diff = 0.d0
  !do ii = 1, nBas
  !  do jj = 1, nBas
  !    diff_loc = dabs(K(jj,ii) - K_deb(jj,ii))
  !    if(diff_loc .gt. 1d-10) then
  !      print*, 'error in K on: ', jj, ii
  !      print*, K(jj,ii), K_deb(jj,ii)
  !      stop
  !    endif
  !    diff = diff + diff_loc
  !  enddo
  !enddo
  !print*, 'total diff on K = ', diff
  !stop

! Initialization

  n_diis        = 0
  F_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0
  rcond         = 0d0

  Conv = 1d0
  nSCF = 0

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
    call Hartree_matrix_AO_basis_hpc (nBas, ERI_size, P(1,1), ERI_chem(1), J(1,1))
    call exchange_matrix_AO_basis_hpc(nBas, ERI_size, P(1,1), ERI_chem(1), K(1,1))
    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)

    ! Check convergence 
    call dgemm("N", "N", nBas, nBas, nBas, 1.d0, &
               S(1,1), nBas, P(1,1), nBas,       &
               0.d0, tmp1(1,1), nBas)
    call dgemm("N", "N", nBas, nBas, nBas, 1.d0, &
               tmp1(1,1), nBas, F(1,1), nBas,    &
               0.d0, err(1,1), nBas)
    call dgemm("N", "T", nBas, nBas, nBas, 1.d0, &
               F(1,1), nBas, tmp1(1,1), nBas,    &
               -1.d0, err(1,1), nBas)
    !err = matmul(F, matmul(P, S)) - matmul(matmul(S, P), F)
    if(nSCF > 1) Conv = maxval(abs(err))

    ! Kinetic energy
    call dgemm("N", "N", nBas, nBas, nBas, 1.d0, &
               P(1,1), nBas, T(1,1), nBas,       &
               0.d0, tmp1(1,1), nBas)
    ET = trace_matrix(nBas, tmp1(1,1))

    ! Potential energy
    call dgemm("N", "N", nBas, nBas, nBas, 1.d0, &
               P(1,1), nBas, V(1,1), nBas,       &
               0.d0, tmp1(1,1), nBas)
    EV = trace_matrix(nBas, tmp1(1,1))

    ! Hartree energy
    call dgemm("N", "N", nBas, nBas, nBas, 1.d0, &
               P(1,1), nBas, J(1,1), nBas,       &
               0.d0, tmp1(1,1), nBas)
    EJ = 0.5d0*trace_matrix(nBas, tmp1(1,1))

    ! Exchange energy
    call dgemm("N", "N", nBas, nBas, nBas, 1.d0, &
               P(1,1), nBas, K(1,1), nBas,       &
               0.d0, tmp1(1,1), nBas)
    EK = 0.25d0*trace_matrix(nBas, tmp1(1,1))

    ! Total energy
    ERHF = ET + EV + EJ + EK

    ! DIIS extrapolation
    if(max_diis > 1) then
      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBas_Sq,nBas_Sq,n_diis,err_diis,F_diis,err,F)
    endif

    ! Level shift
    if(level_shift > 0d0 .and. Conv > thresh) then
      call level_shifting(level_shift,nBas,nOrb,nO,S,c,F)
    endif

    ! Diagonalize Fock matrix
    call dgemm("N", "N", nBas, nOrb, nBas, 1.d0, &
               F(1,1), nBas, X(1,1), nBas,       &
               0.d0, FX(1,1), nBas)
    call dgemm("T", "N", nOrb, nOrb, nBas, 1.d0, &
               X(1,1), nBas, FX(1,1), nBas,      &
               0.d0, Fp(1,1), nOrb)
    !Fp = matmul(transpose(X), matmul(F, X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nOrb,cp,eHF)
    !c = matmul(X, cp)
    call dgemm("N", "N", nBas, nOrb, nOrb, 1.d0, &
               X(1,1), nBas, cp(1,1), nOrb,      &
               0.d0, c(1,1), nBas)
    

    ! Density matrix
    call dgemm('N', 'T', nBas, nBas, nO, 2.d0, &
               c(1,1), nBas, c(1,1), nBas,     &
               0.d0, P(1,1), nBas)
 
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

    deallocate(J,K,err,cp,Fp,err_diis,F_diis)
    deallocate(tmp1, FX, ERI_chem)

    stop

  end if

! Compute dipole moments

  call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_RHF(nBas,nOrb,nO,eHF,c,ENuc,ET,EV,EJ,EK,ERHF,dipole)

! Testing zone

  if(dotest) then
 
    call dump_test_value('R','RHF energy',ERHF)
    call dump_test_value('R','RHF HOMO energy',eHF(nO))
    call dump_test_value('R','RHF LUMO energy',eHF(nO+1))
    call dump_test_value('R','RHF dipole moment',norm2(dipole))

  end if

  deallocate(J,K,err,cp,Fp,err_diis,F_diis)
  deallocate(tmp1, FX, ERI_chem)

end subroutine 
