subroutine MOM_cRHF(dotest,maxSCF,thresh,max_diis,guess_type,level_shift,writeMOs,ENuc, & 
                nBas,nO,S,T,V,ERI,CAP,X,ERHF,eHF,c,P,F,occupationsGuess)

! Perform complex restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest,writeMOs

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: occupationsGuess(nO)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: CAP(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  
  ! Local variables

  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: n_diis
  integer,allocatable           :: occupations(:)
  
  complex*16                    :: ET
  complex*16                    :: EV
  complex*16                    :: EJ
  complex*16                    :: EK
  complex*16                    :: EW

  double precision              :: Conv
  double precision              :: rcond
  double precision,external     :: trace_matrix
  complex*16,external           :: complex_trace_matrix

  complex*16,allocatable        :: J(:,:)
  complex*16,allocatable        :: K(:,:)
  complex*16,allocatable        :: cp(:,:)
  complex*16,allocatable        :: cGuess(:,:)
  complex*16,allocatable        :: Fp(:,:)
  complex*16,allocatable        :: err(:,:)
  complex*16,allocatable        :: err_diis(:,:)
  complex*16,allocatable        :: F_diis(:,:)
  complex*16,allocatable        :: Hc(:,:)
  double precision,allocatable  :: tmp(:,:)


! Output variables

  complex*16,intent(out)        :: ERHF
  complex*16,intent(out)        :: eHF(nBas)
  complex*16,intent(inout)      :: c(nBas,nBas)
  complex*16,intent(out)        :: P(nBas,nBas)
  complex*16,intent(inout)      :: F(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'*************************************'
  write(*,*)'* Complex Restricted HF Calculation *'
  write(*,*)'*************************************'
  write(*,*)

! Useful quantities

  nBasSq = nBas*nBas

! Memory allocation

 allocate(err_diis(nBasSq,max_diis))
 allocate(F_diis(nBasSq,max_diis))
 allocate(Hc(nBas,nBas))
 allocate(J(nBas,nBas))
 allocate(K(nBas,nBas))
 allocate(err(nBas,nBas))
 allocate(cp(nBas,nBas))
 allocate(Fp(nBas,nBas))
 allocate(cGuess(nBas,nBas),occupations(nO))

! Define core Hamiltonian with CAP part
  Hc(:,:) = cmplx(T+V,CAP,kind=8)
 
! Guess coefficients and density matrix
  if(guess_type ==6) then
    allocate(tmp(nBas,nBas))
    call read_matin(nBas,nBas,tmp,"real_MOs_alpha.dat")
    c(:,:) = cmplx(tmp, 0d0,kind=8)
    call read_matin(nBas,nBas,tmp,"imag_MOs_alpha.dat")
    c(:,:) = c(:,:) + cmplx(0d0, tmp, kind=8)
    deallocate(tmp)
  end if

! Guess coefficients and density matrix
  cGuess = c 
  occupations = occupationsGuess
  
  print *, "Ground state orbital occupations for MOM-guess:"
  print *, "Alpha:"
  print *, occupations(1:nO)
  print *, "Beta:"
  print *, occupations(1:nO)
 
  P(:,:) = 2d0*matmul(c(:,occupations(1:nO)),transpose(c(:,occupations(1:nO))))

  ! Initialization
  n_diis          = 0
  F_diis(:,:)   = cmplx(0d0,0d0,kind=8)
  err_diis(:,:) = cmplx(0d0,0d0,kind=8)
  rcond = 0d0 
  Conv   = 1d0
  nSCF   = 0
!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A36,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(RHF)','|','Re(EJ(RHF))','|','Re(EK(RHF))','|','Conv','|'
  write(*,*)'-------------------------------------------------------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Increment 

    nSCF = nSCF + 1
    ! Build Fock matrix
     
    call complex_Hartree_matrix_AO_basis(nBas,P,ERI,J)
    call complex_exchange_matrix_AO_basis(nBas,P,ERI,K)
    
    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)

    ! Check convergence 

    err = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
    if(nSCF > 1) Conv = maxval(abs(err))

    ! Kinetic energy

    ET = cmplx(trace_matrix(nBas,real(matmul(P,T))),trace_matrix(nBas,aimag(matmul(P,T))),kind=8)
    
    ! Potential energy

    EV = cmplx(trace_matrix(nBas,real(matmul(P,V))),trace_matrix(nBas,aimag(matmul(P,V))),kind=8)

    ! CAP energy

    EW = complex_trace_matrix(nBas,matmul(P,(0d0,1d0)*CAP))
    
    ! Hartree energy

    EJ = 0.5d0*cmplx(trace_matrix(nBas,real(matmul(P,J))),trace_matrix(nBas,aimag(matmul(P,J))),kind=8)
    

    ! Exchange energy

    EK = 0.25d0*cmplx(trace_matrix(nBas,real(matmul(P,K))),trace_matrix(nBas,aimag(matmul(P,K))),kind=8)

    ! Total energy

    ERHF = ET + EV + EW + EJ + EK

    ! DIIS extrapolation  !

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call C_complex_DIIS_extrapolation(rcond,nBasSq,nBasSq,n_diis,err_diis,F_diis,err,F)
    end if
    ! Level shift
    if(level_shift > 0d0 .and. Conv > thresh) call complex_level_shifting(level_shift,nBas,nBas,nO,S,c,F)
    
    ! Diagonalize Fock matrix

    Fp = matmul(transpose(X(:,:)),matmul(F(:,:),X(:,:)))
    cp(:,:) = Fp(:,:)
    call complex_diagonalize_matrix(nBas,cp,eHF)
    call complex_orthogonalize_matrix(nBas,cp)
    c = matmul(X,cp)

    ! Density matrix with MOM
    call complex_MOM_density_matrix(nBas, nBas, nO, S, c, cGuess, occupations, occupationsGuess, P)
    P = 2*P
    ! Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,A1,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,E10.2,1X,A1,1X)') &
     '|',nSCF,'|',real(ERHF + ENuc),'+',aimag(ERHF),'i','|',real(EJ),'|',real(EK),'|',Conv,'|'
  end do
  write(*,*)'-------------------------------------------------------------------------------------------------'
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

    stop

  end if

  call print_cRHF(nBas,nBas,nO,eHF,C,ENuc,ET,EV,EW,EJ,EK,ERHF)

! Write MOs
  if(writeMOs) then
        call write_matout(nBas,nBas,real(c),'real_MOs_alpha.dat')
        call write_matout(nBas,nBas,aimag(c),'imag_MOs_alpha.dat')
        call write_matout(nBas,nBas,real(c),'real_MOs_beta.dat')
        call write_matout(nBas,nBas,aimag(c),'imag_MOs_beta.dat')
  endif

  ! Testing zone

  if(dotest) then
 
    print *, "Test for cRHF not implemented"
!   call dump_test_value('R','RHF energy',ERHF)
!   call dump_test_value('R','RHF HOMO energy',eHF(nO))
!   call dump_test_value('R','RHF LUMO energy',eHF(nO+1))
!   call dump_test_value('R','RHF dipole moment',norm2(dipole))

  end if
  deallocate(J,K,err,cp,Fp,err_diis,F_diis,Hc)
end subroutine 
