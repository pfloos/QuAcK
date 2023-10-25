subroutine ROHF(maxSCF,thresh,max_diis,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
                nBas,nO,S,T,V,Hc,ERI,dipole_int,X,EHF,e,c,Ptot)

! Perform restricted open-shell Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: mix 
  double precision,intent(in)   :: level_shift
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: nBas

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nO(nspin)
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
  double precision              :: conv
  double precision              :: rcond
  double precision              :: ET(nspin)
  double precision              :: EV(nspin)
  double precision              :: EJ(nsp)
  double precision              :: Ex(nspin)
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

  double precision,intent(out)  :: EHF
  double precision,intent(out)  :: e(nBas)
  double precision,intent(out)  :: c(nBas,nBas)
  double precision,intent(out)  :: Ptot(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'*      Restricted Open-Shell Hartree-Fock      *'
  write(*,*)'************************************************'
  write(*,*)

! Useful stuff

  nBasSq = nBas*nBas

! Memory allocation

  allocate(J(nBas,nBas,nspin),F(nBas,nBas,nspin),Fp(nBas,nBas),Ftot(nBas,nBas),   & 
           P(nBas,nBas,nspin),K(nBas,nBas,nspin),err(nBas,nBas),cp(nBas,nBas), &
           err_diis(nBasSq,max_diis),F_diis(nBasSq,max_diis))

! Guess coefficients and demsity matrices

  call mo_guess(nBas,guess_type,S,Hc,X,c)
  do ispin=1,nspin
    P(:,:,ispin) = matmul(c(:,1:nO(ispin)),transpose(c(:,1:nO(ispin))))
  end do

! Initialization

  n_diis          = 0
  F_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0
  rcond         = 0d0

  nSCF = 0
  conv = 1d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'----------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X)') & 
            '|','#','|','E(ROHF)','|','Ex(ROHF)','|','Conv','|'
  write(*,*)'----------------------------------------------------------'
  
  do while(conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Build Coulomb repulsion

    do ispin=1,nspin
      call Coulomb_matrix_AO_basis(nBas,P(:,:,ispin),ERI(:,:,:,:),J(:,:,ispin))
    end do

!   Compute exchange potential

    do ispin=1,nspin
      call exchange_matrix_AO_basis(nBas,P(:,:,ispin),ERI(:,:,:,:),K(:,:,ispin))
    end do
 
!   Build Fock operator

    do ispin=1,nspin
      F(:,:,ispin) = Hc(:,:) + J(:,:,ispin) + J(:,:,mod(ispin,2)+1) + K(:,:,ispin)
    end do

    call ROHF_fock_matrix(nBas,nO(1),nO(2),S,c,F(:,:,1),F(:,:,2),Ftot)

!   Check convergence 

    err(:,:) = matmul(Ftot(:,:),matmul(Ptot(:,:),S(:,:))) - matmul(matmul(S(:,:),Ptot(:,:)),Ftot(:,:))
    if(nSCF > 1) conv = maxval(abs(err(:,:)))
    
!   DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBasSq,nBasSq,n_diis,err_diis(:,1:n_diis),F_diis(:,1:n_diis),err,Ftot)

    end if

!   Level-shifting

    if(level_shift > 0d0 .and. Conv > thresh) then

      do ispin=1,nspin
        call level_shifting(level_shift,nBas,maxval(nO),S,c,Ftot)
      end do

    end if

!  Transform Fock matrix in orthogonal basis

    Fp(:,:) = matmul(transpose(X(:,:)),matmul(Ftot(:,:),X(:,:)))

!  Diagonalize Fock matrix to get eigenvectors and eigenvalues

    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp,e)
    
!   Back-transform eigenvectors in non-orthogonal basis

    c(:,:) = matmul(X(:,:),cp(:,:))

!   Compute density matrix 

    do ispin=1,nspin
      P(:,:,ispin) = matmul(c(:,1:nO(ispin)),transpose(c(:,1:nO(ispin))))
    end do
    Ptot(:,:) = P(:,:,1) + P(:,:,2) 

!------------------------------------------------------------------------
!   Compute ROHF energy
!------------------------------------------------------------------------

!  Kinetic energy

    do ispin=1,nspin
      ET(ispin) = trace_matrix(nBas,matmul(P(:,:,ispin),T(:,:)))
    end do

!  Potential energy

    do ispin=1,nspin
      EV(ispin) = trace_matrix(nBas,matmul(P(:,:,ispin),V(:,:)))
    end do

!  Coulomb energy

    EJ(1) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,1)))
    EJ(2) = trace_matrix(nBas,matmul(P(:,:,1),J(:,:,2)))
    EJ(3) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,2)))

!   Exchange energy

    do ispin=1,nspin
      Ex(ispin) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,ispin),K(:,:,ispin)))
    end do

!   Total energy

    EHF = sum(ET(:)) + sum(EV(:)) + sum(EJ(:)) + sum(Ex(:))

!   Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X)') & 
      '|',nSCF,'|',EHF + ENuc,'|',sum(Ex(:)),'|',conv,'|'
 
  end do
  write(*,*)'----------------------------------------------------------'
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

! Compute final UHF energy

  call dipole_moment(nBas,Ptot,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_ROHF(nBas,nO,e,c,ENuc,ET,EV,EJ,Ex,EHF,dipole)

end subroutine 
