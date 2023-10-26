subroutine RMOM(maxSCF,thresh,max_diis,nBas,nO,S,T,V,Hc,ERI,X,ENuc,ERHF,c,e,P)

! Maximum overlap method

  implicit none

! Input variables

  integer,intent(in)            :: maxSCF,max_diis
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBas,nO
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas) 
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas),X(nBas,nBas)

! Local variables

  integer                       :: iBas,jBas
  integer                       :: nSCF,nBasSq,n_diis
  double precision              :: ET,EV,EJ,EK,Conv,Gap 
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision,allocatable  :: error(:,:),error_diis(:,:),F_diis(:,:)
  double precision,allocatable  :: J(:,:),K(:,:),cp(:,:),F(:,:),Fp(:,:)
  double precision,allocatable  :: cG(:,:),ON(:)

! Output variables

  double precision,intent(inout):: ERHF,c(nBas,nBas),e(nBas),P(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|           Maximum overlap method             |'
  write(*,*)'************************************************'
  write(*,*)

! Useful quantities

  nBasSq = nBas*nBas

! Memory allocation

  allocate(J(nBas,nBas),K(nBas,nBas),error(nBas,nBas), &
           cp(nBas,nBas),Fp(nBas,nBas),F(nBas,nBas),   &
           cG(nBas,nBas),ON(nBas), &
           error_diis(nBasSq,max_diis),F_diis(nBasSq,max_diis))

! Set up guess orbitals
  
  cG(:,:) = c(:,:)

! Set up occupation numbers

  ON(1:nO) = 1d0
  ON(nO+1:nBas) = 0d0

! HOMO-LUMO transition

  ON(nO) = 0d0
  ON(nO+1) = 1d0

  write(*,*)
  write(*,*) ' --- Initial MO occupations --- '
  write(*,*)
  call matout(nBas,1,ON)
  write(*,*)

! Compute density matrix
  
  call density_matrix(nBas,ON,c,P)

! Initialization

  n_diis = 0
  F_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  Conv = 1d0
  nSCF = 0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| MOM calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') & 
            '|','#','|','HF energy','|','Conv','|','HL Gap','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Build Fock matrix
    
    call Hartree_matrix_AO_basis(nBas,P,ERI,J)
    call exchange_matrix_AO_basis(nBas,P,ERI,K)
    
    F(:,:) = Hc(:,:) + J(:,:) + 0.5*K(:,:)

!   Check convergence 

    error = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
    Conv = maxval(abs(error))

!   DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nBasSq,nBasSq,n_diis,error_diis,F_diis,error,F)

!   Reset DIIS if required

    if(abs(rcond) < 1d-15) n_diis = 0

!  Diagonalize Fock matrix

    Fp = matmul(transpose(X),matmul(F,X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp,e)
    c = matmul(X,cp)

!   MOM overlap

    call MOM_overlap(nBas,nO,S,cG,c,ON)

!   Density matrix

    call density_matrix(nBas,ON,c,P)

!   Compute HF energy

    ERHF = trace_matrix(nBas,matmul(P,Hc))      &
         + 0.5d0*trace_matrix(nBas,matmul(P,J)) &
         + 0.5d0*trace_matrix(nBas,matmul(P,K))

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

  write(*,*) 
  write(*,*) ' --- Final MO occupations --- '
  write(*,*)
  call matout(nBas,1,ON)
  write(*,*)

! Compute HF energy

  ET   = trace_matrix(nBas,matmul(P,T))
  EV   = trace_matrix(nBas,matmul(P,V))
  EJ   = 0.5d0*trace_matrix(nBas,matmul(P,J))
  EK   = 0.5d0*trace_matrix(nBas,matmul(P,K))
  ERHF = ET + EV + EJ + EK

  call print_RHF(nBas,nO,e,c,ENuc,ET,EV,EJ,EK,ERHF)

end subroutine 
