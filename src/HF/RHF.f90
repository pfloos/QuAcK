subroutine RHF(maxSCF,thresh,max_diis,guess_type,nNuc,ZNuc,rNuc,ENuc,nBas,nO,S,T,V,Hc,F,ERI,dipole_int,X,ERHF,e,c,P,Vx)

! Perform restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF,max_diis,guess_type
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
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
  double precision              :: Gap 
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision,allocatable  :: error(:,:)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: ON(:)

! Output variables

  double precision,intent(out)  :: ERHF
  double precision,intent(out)  :: e(nBas)
  double precision,intent(out)  :: c(nBas,nBas)
  double precision,intent(out)  :: P(nBas,nBas)
  double precision,intent(out)  :: Vx(nBas)
  double precision,intent(out)  :: F(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|    Restricted Hartree-Fock calculation       |'
  write(*,*)'************************************************'
  write(*,*)

! Useful quantities

  nBasSq = nBas*nBas

! Memory allocation

  allocate(J(nBas,nBas),K(nBas,nBas),error(nBas,nBas),        &
           cp(nBas,nBas),Fp(nBas,nBas),ON(nBas),              &
           error_diis(nBasSq,max_diis),F_diis(nBasSq,max_diis))

! Guess coefficients and eigenvalues

  call mo_guess(nBas,nO,guess_type,S,Hc,ERI,J,K,X,cp,F,Fp,e,c,P)

! ON(:) = 0d0
! do i=1,nO
!    ON(i) = 1d0
!    ON(i) = dble(2*i-1)
! end do

! call density_matrix(nBas,ON,c,P)
  
! Initialization

  n_diis = 0
  F_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  Conv = 1d0
  nSCF = 0
  rcond = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| RHF calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') & 
            '|','#','|','HF energy','|','Conv','|','HL Gap','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Build Fock matrix
    
    call Coulomb_matrix_AO_basis(nBas,P,ERI,J)
    call exchange_matrix_AO_basis(nBas,P,ERI,K)
    
    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)

!   Check convergence 

    error = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
    Conv = maxval(abs(error))

!   DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    if(abs(rcond) > 1d-7) then
      call DIIS_extrapolation(rcond,nBasSq,nBasSq,n_diis,error_diis,F_diis,error,F)
    else
      n_diis = 0
    end if

!  Diagonalize Fock matrix

    Fp = matmul(transpose(X),matmul(F,X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp,e)
    c = matmul(X,cp)

!   Density matrix

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

!   call density_matrix(nBas,ON,c,P)
  
!   Compute HF energy

    ERHF = trace_matrix(nBas,matmul(P,Hc)) &
        + 0.5d0*trace_matrix(nBas,matmul(P,J)) &
        + 0.25d0*trace_matrix(nBas,matmul(P,K))

!   Compute HOMO-LUMO gap

    if(nBas > nO) then 

      Gap = e(nO+1) - e(nO)

    else

      Gap = 0d0

    endif

!  Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') & 
      '|',nSCF,'|',ERHF+ENuc,'|',Conv,'|',Gap,'|'

  enddo
  write(*,*)'----------------------------------------------------'
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

  endif

! Compute HF energy

  ET = trace_matrix(nBas,matmul(P,T))
  EV = trace_matrix(nBas,matmul(P,V))
  EJ = 0.5d0*trace_matrix(nBas,matmul(P,J))
  EK = 0.25d0*trace_matrix(nBas,matmul(P,K))
  ERHF = ET + EV + EJ + EK

! Compute dipole moments

  call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_RHF(nBas,nO,e,C,ENuc,ET,EV,EJ,EK,ERHF,dipole)

! Compute Vx for post-HF calculations

  call mo_fock_exchange_potential(nBas,c,K,Vx)

end subroutine RHF
