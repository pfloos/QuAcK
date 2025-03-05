subroutine cRHF(dotest,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
                nBas,nO,S,T,V,ERI,dipole_int,X,ERHF,eHF,c,P)

! Perform complex restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: n_diis
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: EK
  double precision              :: dipole(ncart)

  double precision              :: Conv
  double precision              :: rcond
  double precision,external     :: trace_matrix

  double precision              :: eta
  double precision,allocatable  :: W(:,:)
  complex*16,allocatable        :: Hc(:,:)
  complex*16,allocatable        :: J(:,:)
  complex*16,allocatable        :: K(:,:)
  complex*16,allocatable        :: cp(:,:)
  complex*16,allocatable        :: F(:,:)
  complex*16,allocatable        :: Fp(:,:)
  complex*16,allocatable        :: err(:,:)
  complex*16,allocatable        :: err_diis(:,:)
  complex*16,allocatable        :: F_diis(:,:)

! Output variables

  complex*16,intent(out)        :: ERHF
  complex*16,intent(out)        :: eHF(nBas)
  complex*16,intent(inout)      :: c(nBas,nBas)
  complex*16,intent(out)        :: P(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'*************************************'
  write(*,*)'* Complex Restricted HF Calculation *'
  write(*,*)'*************************************'
  write(*,*)

! Useful quantities

  nBasSq = nBas*nBas
  eta = 0.01d0

! Memory allocation

  allocate(J(nBas,nBas),K(nBas,nBas),err(nBas,nBas),cp(nBas,nBas),F(nBas,nBas), &
           Fp(nBas,nBas),err_diis(nBasSq,max_diis),F_diis(nBasSq,max_diis),     & 
           Hc(nBas,nBas),W(nBas,nBas))

! Read CAP integrals from file
  call read_CAP_integrals(nBas,W)
  W(:,:) = -eta*W(:,:)

! Define core Hamiltonian with CAP part

  Hc(:,:) = cmplx(T+V,W,kind=8)

! Guess coefficients and density matrix

  call complex_mo_guess(nBas,nBas,guess_type,S,Hc,X,c)
  P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

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
  write(*,*)'--------------------------------------------------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A36,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(RHF)','|','EJ(RHF)','|','EK(RHF)','|','Conv','|'
  write(*,*)'--------------------------------------------------------------------------------------------------'

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

    ET = trace_matrix(nBas,matmul(P,T))
    ! Potential energy

    EV = trace_matrix(nBas,matmul(P,V))

    ! Hartree energy

    EJ = 0.5d0*trace_matrix(nBas,matmul(P,J))

    ! Exchange energy

    EK = 0.25d0*trace_matrix(nBas,matmul(P,K))

    ! Total energy

    ERHF = ET + EV + EJ + EK

!    ! DIIS extrapolation (fix later) !
!
!    if(max_diis > 1) then
!
!      n_diis = min(n_diis+1,max_diis)
!      call complex_DIIS_extrapolation(rcond,nBasSq,nBasSq,n_diis,err_diis,F_diis,err,F)
!
!    end if
!
    ! Level shift
    if(level_shift > 0d0 .and. Conv > thresh) call level_shifting(level_shift,nBas,nBas,nO,S,c,F)
    
    ! Diagonalize Fock matrix

    Fp = matmul(transpose(X),matmul(F,X))
    cp(:,:) = Fp(:,:)
    write(*,*) nBas
    call complex_diagonalize_matrix(nBas,cp,eHF)
    c = matmul(X,cp)
    ! Density matrix

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

    ! Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F32.10,1X,A1,1X,F32.10,A1,1X,A1,1X,F32.10,1X,A1,1X,F32.10,1X,A1,1X,E10.2,1X,A1,1X)') &
      '|',nSCF,'|',real(ERHF),'+',aimag(ERHF),'i','|',EJ,'|',EK,'|',Conv,'|'
    write(*,*) real(ERHF),'+',aimag(ERHF),'i'
  end do
  write(*,*)'--------------------------------------------------------------------------------------------------'
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

! Compute dipole moments

  call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_cRHF(nBas,nBas,nO,eHF,C,ENuc,ET,EV,EJ,EK,ERHF,dipole)
! Testing zone

  if(dotest) then
 
!   call dump_test_value('R','RHF energy',ERHF)
!   call dump_test_value('R','RHF HOMO energy',eHF(nO))
!   call dump_test_value('R','RHF LUMO energy',eHF(nO+1))
!   call dump_test_value('R','RHF dipole moment',norm2(dipole))

  end if

end subroutine 
