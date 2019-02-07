subroutine UHF(nEl,nBas,S,T,V,Hc,G,X,ENuc,EHF,c,e,P,F)

! Perform unrestricted Hartree-Fock calculation

  implicit none

! Input variables

  integer,intent(in)            :: nEl(nspin),nBas
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: S(nBas,nBas),T(nBas,nBas),V(nBas,nBas),Hc(nBas,nBas) 
  double precision,intent(in)   :: G(nBas,nBas,nBas,nBas),X(nBas,nBas)

! Local variables

  logical                       :: random_guess,core_guess,DIIS
  integer,parameter             :: maxSCF = 64
  double precision,parameter    :: thresh = 1d-6
  integer                       :: nO,nSCF,nBasSq,n_diis,ispin,jspin
  double precision              :: ET(spin),EV(nspin),EJ(nspin,nspin),EK(nspin,nspin),Conv(nspin),Gap(spin)
  double precision              :: trace_matrix
  double precision,allocatable  :: FPS_SPF(:,:,:),error_diis(:,:,:),F_diis(:,:,:)
  double precision,allocatable  :: J(:,:,:),K(:,:,:),cp(:,:,:),cO(:,:,:)

! Output variables

  double precision,intent(out)  :: EHF,c(nBas,nBas,nspin),e(nBas,nspin),P(nBas,nBas,nspin),F(nBas,nBas,nspin)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|  Unrestricted Hartree-Fock calculation       |'
  write(*,*)'************************************************'
  write(*,*)

! Useful quantities

  nBasSq = nBas*nBas

! Number of occupied orbitals

! Initialize DIIS variables 

  DIIS = .false.
  n_diis = 5
  if(.not.DIIS) n_diis = 1

! Type of guess 
  
  random_guess = .false.
  core_guess   = .true.

! Memory allocation

  allocate(J(nBas,nBas,nspin),K(nBas,nBas,nspin),      & 
           cp(nBas,nBas,nspin),cO(nBas,nO,nspin),FPS_SPF(nBas,nBas,nspin), &
           error_diis(nBasSq,n_diis,nspin),F_diis(nBasSq,n_diis,spin))

! Guess coefficients and eigenvalues

  if(random_guess) then

    call random_number(c)

  elseif(core_guess) then

    cp(:,:,ispin) = matmul(transpose(X(:,:)),matmul(Hc(:,:),X(:,:)))
    call diagonalize_matrix(nBas,cp(:,:,ispin),e(:,:,ispin))
    c(:,:,ispin) = matmul(X,cp(:,:,ispin)

  endif
  
! Occupied orbitals

  cO(1:nBas,1:nO,1:nspin) = c(1:nBas,1:nO,1:nspin)

! Initialization

  Conv(:) = 1d0
  nSCF = 0
  F_diis(:,:,:)     = 0d0
  error_diis(:,:,:) = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| UHF calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') & 
            '|','#','|','HF energy','|','Conv','|','HL Gap','|'
  write(*,*)'----------------------------------------------------'

  do while(maxval(Conv) > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Density matrix

    P(:,:,ispin) = 2d0*matmul(cO(:,:,ispin),transpose(cO(:,:,ispin)))

!   Build Fock matrix
    
    call Coulomb_matrix_AO_basis(nBas,P(:,:,ispin),G,J(:,:,ispin))
    call exchange_matrix_AO_basis(nBas,P(:,:,ispin),G,K(:,:,ispin))
    
    F(:,:,ispin) = Hc(:,:) + J(:,:,ispin) + J(:,:,mod(ispin,2)+1) + K(:,:,ispin)

!   Check convergence 

    FPS_SPF(:,:,ispin) = matmul(F(:,:,ispin),matmul(P(:,:,ispin),S)) & 
                       - matmul(matmul(S,P(:,:,ispin)),F(:,:,ispin))
    Conv(ispin) = maxval(abs(FPS_SPF(:,:,ispin)))

!   DIIS extrapolation

    call prepend(nBasSq,n_diis,error_diis(:,:,ispin),FPS_SPF(:,:,ispin))
    call prepend(nBasSq,n_diis,F_diis(:,:,ispin),F(:,:,ispin))
    call diis(nBasSq,min(n_diis,nSCF),error_diis(:,:,ispin),F_diis(:,:,ispin),F(:,:,ispin))

!  Diagonalize Fock matrix

    cp(:,:,ispin) = matmul(transpose(X),matmul(F(:,:,ispin),X))
    call diagonalize_matrix(nBas,cp(:,:,ispin),e(:,:,ispin))

    c(:,:,ispin) = matmul(X,cp(:,:,ispin))
    cO(1:nBas,1:nO,1:nspin) = c(1:nBas,1:nO,1:nspin)

!   Compute HF energy

!   EHF = 0.5d0*trace_matrix(nBas,matmul(P,Hc+F))

!   Compute HOMO-LUMO gap

    if(nBas > nO) then 

      Gap(:) = e(nO+1,:) - e(nO,:)

    else

      Gap(:) = 0d0

    endif

!  Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') & 
      '|',nSCF,'|',EHF+ENuc,'|',maxval(Conv),'|',minval(Gap),'|'

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

  EHF = 0d0
  do ispin=1,nspin

    ET(ispin) = trace_matrix(nBas,matmul(P(:,:,ispin),T(:,:,ispin))
    EV(ispin) = trace_matrix(nBas,matmul(P(:,:,ispin),V(:,:,ispin))

    EHF = EHF + ET(ispin) + EV(ispin)

    do jspin=1,nspin

      EJ(ispin,jspin) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,ispin),J(:,:,jspin))
      EK(ispin,jspin) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,ispin,K(:,:,jspin)))

      EHF = EHF + EJ(ispin,jspin) + EK(ispin,jspin)

    enddo

  enddo


! call print_UHF(nBas,nO,e,C,ENuc,ET,EV,EJ,EK,EHF)

end subroutine UHF
