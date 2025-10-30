subroutine ensembleRHF(dotest,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
                       nBas,nOrb,nO,S,T,V,Hc,ERI,dipole_int,X,eERHF,eweight,eforward,eHF,c,P_tot,F)

! P_toterform restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest
  logical,intent(in)            :: eforward

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: eweight
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
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: file_exists
  integer                       :: iunit,iunit2
  integer                       :: iorb
  integer                       :: ibas,jbas,kbas,lbas
  integer                       :: nSCF
  integer                       :: nBas_Sq
  integer                       :: n_diis
  double precision              :: trace_1rdm
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: EK
  double precision              :: Ecore
  double precision              :: Eee
  double precision              :: dipole(ncart)

  double precision              :: Val
  double precision              :: Conv
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision,allocatable  :: Occ(:)
  double precision,allocatable  :: err(:,:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: P_delta(:,:)
  double precision,allocatable  :: P_mo(:,:)

! Output variables

  double precision,intent(out)  :: eERHF
  double precision,intent(out)  :: eHF(nOrb)
  double precision,intent(inout):: c(nBas,nOrb)
  double precision,intent(out)  :: P_tot(nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'****************************'
  write(*,*)'* Ensemble RHF Calculation *'
  write(*,*)'****************************'
  write(*,*)

! Useful quantities
  nBas_Sq = nBas*nBas

! Memory allocation

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))

  allocate(err(nBas,nBas))
  allocate(P_delta(nBas,nBas))

  allocate(Occ(nOrb))
  allocate(P_mo(nOrb,nOrb))
  allocate(cp(nOrb,nOrb))
  allocate(Fp(nOrb,nOrb))

  allocate(err_diis(nBas_Sq,max_diis))
  allocate(F_diis(nBas_Sq,max_diis))

! Guess coefficients and density matrix

  call mo_guess(nBas,nOrb,guess_type,S,Hc,X,c)
  P_tot(:,:) = 2d0 * matmul(c(:,1:nO), transpose(c(:,1:nO)))
  if(eweight>1d-8) then
   if(eforward) then
    P_delta(:,:) = 0d0
    if(nO==nOrb) then
     write(*,*) ' There are no virtual orbitals to fill in forward mode. Hint: Increase the size of the basis.'
     stop
    endif
    if(nO+1<=nOrb) then
     P_delta(:,:) = 2d0*matmul(c(:,1:nO+1),transpose(c(:,1:nO+1)))
    endif
    P_tot(:,:) = (1d0-eweight)*P_tot(:,:) + eweight*P_delta(:,:)
   else
    P_delta(:,:) = 0d0
    if(nO>=2) then
     P_delta(:,:) = 2d0*matmul(c(:,1:nO-1),transpose(c(:,1:nO-1)))
    endif
    P_tot(:,:) = (1d0-eweight)*P_tot(:,:) + eweight*P_delta(:,:)
   endif
  endif
  if(guess_type == 5) then
   inquire(file='P_tot_ao_bin', exist=file_exists)
   if(file_exists) then
    write(*,*) 'Reading P_tot_ao_bin matrix...'
    open(unit=314, form='unformatted', file='P_tot_ao_bin', status='old')
    do
     read(314) ibas,jbas,Val
     if(ibas==0 .and. jbas==0) exit
     P_tot(ibas,jbas)=Val 
    enddo
    close(314)
   endif
  endif

! call dgemm('N', 'T', nBas, nBas, nO, 2.d0, &
!            c(1,1), nBas, c(1,1), nBas,     &
!            0.d0, P_tot(1,1), nBas)

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
            '|','#','|','E(eRHF)','|','EJ(eRHF)','|','EK(eRHF)','|','Conv','|'
  write(*,*)'-----------------------------------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Increment 

    nSCF = nSCF + 1

    ! Build Fock matrix
    
    call Hartree_matrix_AO_basis(nBas,P_tot,ERI,J)
    call exchange_matrix_AO_basis(nBas,P_tot,ERI,K)
    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)

    ! Check convergence 

    err = matmul(F,matmul(P_tot,S)) - matmul(matmul(S,P_tot),F)
    if(nSCF > 2) Conv = maxval(abs(err))

    ! Kinetic energy

    ET = trace_matrix(nBas,matmul(P_tot,T))

    ! P_tototential energy

    EV = trace_matrix(nBas,matmul(P_tot,V))

    ! Hartree energy

    EJ = 0.5d0*trace_matrix(nBas,matmul(P_tot,J))

    ! Exchange energy

    EK = 0.25d0*trace_matrix(nBas,matmul(P_tot,K))

    ! Total energy

    eERHF = ET + EV + EJ + EK

    ! DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBas_Sq,nBas_Sq,n_diis,err_diis,F_diis,err,F)

    end if

    ! Level shift

    if(level_shift > 0d0 .and. Conv > thresh) then
      call level_shifting(level_shift,nBas,nOrb,nO,S,c,F)
    endif

    ! Diagonalize Fock matrix

      Fp = matmul(transpose(X),matmul(F,X))
      cp(:,:) = Fp(:,:)
      call diagonalize_matrix(nOrb,cp,eHF)
      c = matmul(X,cp)

    ! Density matrix

    P_tot(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))
    if(eweight>1d-8) then
     if(eforward) then
      P_delta(:,:) = 0d0
      if(nO+1<=nOrb) then
       P_delta(:,:) = 2d0*matmul(c(:,1:nO+1),transpose(c(:,1:nO+1)))
      endif
      P_tot(:,:) = (1d0-eweight)*P_tot(:,:) + eweight*P_delta(:,:)
     else
      P_delta(:,:) = 0d0
      if(nO>=2) then
       P_delta(:,:) = 2d0*matmul(c(:,1:nO-1),transpose(c(:,1:nO-1)))
      endif
      P_tot(:,:) = (1d0-eweight)*P_tot(:,:) + eweight*P_delta(:,:)
     endif
    endif
!   call dgemm('N', 'T', nBas, nBas, nO, 2.d0, &
!              c(1,1), nBas, c(1,1), nBas,     &
!              0.d0, P_tot(1,1), nBas)

    ! Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,E10.2,1X,A1,1X)') &
      '|',nSCF,'|',eERHF + ENuc,'|',EJ,'|',EK,'|',Conv,'|'

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

    !deallocate(J,K,err,cp,Fp,err_diis,F_diis)

    write(*,*) ' Warning! Convergence failed at Hartree-Fock level.'

  end if

! Compute dipole moments

  call Hartree_matrix_AO_basis(nBas,P_tot,ERI,J)
  call exchange_matrix_AO_basis(nBas,P_tot,ERI,K)
  F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)
  Fp = matmul(transpose(X),matmul(F,X))
  cp(:,:) = Fp(:,:)
  call diagonalize_matrix(nOrb,cp,eHF)
  c = matmul(X,cp)

  call dipole_moment(nBas,P_tot,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_eRHF(nBas, nOrb, nO, eHF, c, ENuc, ET, EV, EJ, EK, eERHF, dipole)
  
! Build NOs and occ numbers
  P_mo = -matmul(transpose(X),matmul(P_tot,X))
  call diagonalize_matrix(nOrb,P_mo,Occ)
  write(*,*)
  write(*,*) ' Occupation numbers'
  Occ=-Occ
  trace_1rdm=0d0
  do ibas=1,nOrb
   write(*,'(I7,F15.8)') ibas,Occ(ibas)
   trace_1rdm=trace_1rdm+Occ(ibas)
  enddo
  write(*,'(A33,1X,F16.10,A3)') ' Trace [ 1D^NO ]     = ',trace_1rdm,'   '

! Testing zone

  if(dotest) then

! TODO
 
  end if

! Memory deallocation

  deallocate(J,K,err,P_delta,cp,Fp,err_diis,F_diis)
  deallocate(P_mo,Occ)


end subroutine 

