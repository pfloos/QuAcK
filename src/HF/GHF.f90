subroutine GHF(maxSCF,thresh,max_diis,guess_type,mix,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
               nBas,nBas2,nO,S,T,V,Hc,ERI,dipole_int,X,EHF,e,c,P)

! Perform unrestricted Hartree-Fock calculation

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
  integer,intent(in)            :: nBas2

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
  integer                       :: nBas2Sq
  integer                       :: nOa
  integer                       :: nOb
  integer                       :: n_diis
  double precision              :: Conv
  double precision              :: rcond
  double precision              :: ET(nspin)
  double precision              :: EV(nspin)
  double precision              :: EJ(nsp)
  double precision              :: Ex(nspin)
  double precision              :: dipole(ncart)

  double precision,allocatable  :: Caa(:,:),Cab(:,:),Cba(:,:),Cbb(:,:)
  double precision,allocatable  :: Jaa(:,:),Jab(:,:),Jba(:,:),Jbb(:,:)
  double precision,allocatable  :: Kaa(:,:),Kab(:,:),Kba(:,:),Kbb(:,:)
  double precision,allocatable  :: Faa(:,:),Fab(:,:),Fba(:,:),Fbb(:,:)
  double precision,allocatable  :: Paa(:,:),Pab(:,:),Pba(:,:),Pbb(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: Cp(:,:)
  double precision,allocatable  :: err(:,:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,external     :: trace_matrix

  integer                       :: ispin

! Output variables

  double precision,intent(out)  :: EHF
  double precision,intent(out)  :: e(nBas2)
  double precision,intent(out)  :: C(nBas2,nBas2)
  double precision,intent(out)  :: P(nBas2,nBas2)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'*    Unrestricted Hartree-Fock calculation     *'
  write(*,*)'************************************************'
  write(*,*)

! Useful stuff

  nBasSq = nBas*nBas
  nBas2Sq = nBas2*nBas2

  nOa = nO(1)
  nOb = nO(2)

! Memory allocation

  allocate(Caa(nBas,nBas),Cab(nBas,nBas),Cba(nBas,nBas),Cbb(nBas,nBas),     &
           Jaa(nBas,nBas),Jab(nBas,nBas),Jba(nBas,nBas),Jbb(nBas,nBas),     &
           Kaa(nBas,nBas),Kab(nBas,nBas),Kba(nBas,nBas),Kbb(nBas,nBas),     &
           Faa(nBas,nBas),Fab(nBas,nBas),Fba(nBas,nBas),Fbb(nBas,nBas),     &
           Paa(nBas,nBas),Pab(nBas,nBas),Pba(nBas,nBas),Pbb(nBas,nBas),     &
           F(nBas2,nBas2),Fp(nBas2,nBas2),Cp(nBas2,nBas2),err(nBas2,nBas2), &
           err_diis(nBas2Sq,max_diis),F_diis(nBas2Sq,max_diis))

! Guess coefficients and demsity matrices

! do ispin=1,nspin
!   call mo_guess(nBas,guess_type,S,Hc,X,c(:,:,ispin))
!   P(:,:,ispin) = matmul(c(:,1:nO(ispin),ispin),transpose(c(:,1:nO(ispin),ispin)))
! end do

! Initialization

  nSCF = 0
  Conv = 1d0

  n_diis        = 0
  F_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0

! Construct super overlap matrix

! TO DO

! Construct super orthogonalization matrix

! TO DO

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  write(*,*)
  write(*,*)'----------------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X)') & 
            '|','#','|','E(UHF)','|','Ex(UHF)','|','Conv','|'
  write(*,*)'----------------------------------------------------------'
  
  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Build individual Hartree matrices

    call Hartree_matrix_AO_basis(nBas,Paa,ERI,Jaa)
    call Hartree_matrix_AO_basis(nBas,Pab,ERI,Jab)
    call Hartree_matrix_AO_basis(nBas,Pba,ERI,Jba)
    call Hartree_matrix_AO_basis(nBas,Pbb,ERI,Jbb)

!   Compute individual exchange matrices

    call exchange_matrix_AO_basis(nBas,Paa,ERI,Kaa)
    call exchange_matrix_AO_basis(nBas,Pab,ERI,Kab)
    call exchange_matrix_AO_basis(nBas,Pba,ERI,Kba)
    call exchange_matrix_AO_basis(nBas,Pbb,ERI,Kbb)
 
!   Build individual Fock matrices

    Faa(:,:) = Hc(:,:) + Jaa(:,:) + Jab(:,:) + Kaa(:,:)
    Fab(:,:) = Hc(:,:) + Jab(:,:) + Jba(:,:) + Kab(:,:)
    Fba(:,:) = Hc(:,:) + Jba(:,:) + Jab(:,:) + Kba(:,:)
    Fbb(:,:) = Hc(:,:) + Jbb(:,:) + Jba(:,:) + Kbb(:,:)

!  Build super Fock matrix

    F(     1:nBas ,     1:nBas ) = Faa(1:nBas,1:nBas)
    F(nBas+1:nBas2,     1:nBas ) = Fab(1:nBas,1:nBas)
    F(     1:nBas ,nBas+1:nBas2) = Fba(1:nBas,1:nBas)
    F(nBas+1:nBas2,nBas+1:nBas2) = Fbb(1:nBas,1:nBas)

!   Check convergence 

!   err(:,:) = matmul(F(:,:),matmul(P(:,:),S(:,:))) - matmul(matmul(S(:,:),P(:,:)),F(:,:))

!   if(nSCF > 1) conv = maxval(abs(err(:,:)))
    
!   DIIS extrapolation

!   if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBas2Sq,nBas2Sq,n_diis,err_diis(:,1:n_diis),F_diis(:,1:n_diis),err,F)

!   end if

!   Level-shifting

    if(level_shift > 0d0 .and. Conv > thresh) then

      call level_shifting(level_shift,nBas,nOa+nOb,S,C,F)

    end if

!  Transform Fock matrix in orthogonal basis

    Fp(:,:) = matmul(transpose(X),matmul(F,X))

!  Diagonalize Fock matrix to get eigenvectors and eigenvalues

    Cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas2,Cp,e)
    
!   Back-transform eigenvectors in non-orthogonal basis

    C(:,:) = matmul(X,Cp)

!   Form individual coefficient matrices

!   TO DO

!   Mix guess for UHF solution in singlet states

!   if(nSCF == 1) call mix_guess(nBas,nO,mix,c)

!   Compute individual density matrices

    Paa(:,:) = matmul(Caa(:,1:nOa),transpose(Caa(:,1:nOa)))
    Pab(:,:) = matmul(Cab(:,1:nOb),transpose(Cab(:,1:nOb)))
    Pba(:,:) = matmul(Cba(:,1:nOa),transpose(Cba(:,1:nOa)))
    Pbb(:,:) = matmul(Cbb(:,1:nOb),transpose(Cbb(:,1:nOb)))
 
!------------------------------------------------------------------------
!   Compute UHF energy
!------------------------------------------------------------------------

!  Kinetic energy

!   do ispin=1,nspin
!     ET(ispin) = trace_matrix(nBas,matmul(P(:,:,ispin),T(:,:)))
!   end do

!  Potential energy

!   do ispin=1,nspin
!     EV(ispin) = trace_matrix(nBas,matmul(P(:,:,ispin),V(:,:)))
!   end do

!  Hartree energy

!   EJ(1) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,1),J(:,:,1)))
!   EJ(2) = trace_matrix(nBas,matmul(P(:,:,1),J(:,:,2)))
!   EJ(3) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,2),J(:,:,2)))

!   Exchange energy

!   do ispin=1,nspin
!     Ex(ispin) = 0.5d0*trace_matrix(nBas,matmul(P(:,:,ispin),K(:,:,ispin)))
!   end do

!   Total energy

!   EHF = sum(ET(:)) + sum(EV(:)) + sum(EJ(:)) + sum(Ex(:))

!   Dump results

!   write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X)') & 
!     '|',nSCF,'|',EHF + ENuc,'|',sum(Ex(:)),'|',conv,'|'
 
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

! call dipole_moment(nBas,P(:,:,1)+P(:,:,2),nNuc,ZNuc,rNuc,dipole_int,dipole)
! call print_GHF(nBas,nO,S,e,c,ENuc,ET,EV,EJ,Ex,EHF,dipole)

end subroutine 
