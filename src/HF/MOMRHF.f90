subroutine MOMRHF(dotest,doaordm,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
               nBas,nOrb,nO,S,T,V,Hc,ERI,dipole_int,X,ERHF,eHF,c,P,F,occupationsGuess)

! Perform restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest
  logical,intent(in)            :: doaordm

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  integer,intent(in)            :: guess_type
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nNuc
  integer,intent(in)            :: occupationsGuess(nO,nspin)
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
  double precision,allocatable  :: err(:,:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: Fp(:,:)

  double precision,allocatable  :: cGuess(:,:)
  integer,allocatable           :: occupations(:)
  integer                       :: i
  logical                       :: alphaEqualsBeta

! Output variables

  double precision,intent(out)  :: ERHF
  double precision,intent(inout):: eHF(nOrb)
  double precision,intent(inout):: c(nBas,nOrb)
  double precision,intent(out)  :: P(nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'****************************************************'
  write(*,*)'* Maximum Overlap Method Restricted HF Calculation *'
  write(*,*)'****************************************************'
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
  
  allocate(cGuess(nBas,nOrb),occupations(nO)) 

! Read input occupations (ground state) and prepare guess
  ! Check if alpha equals beta occupation
  alphaEqualsBeta = .true. 
  do i=1,nO
    if(occupationsGuess(i,1)/= occupationsGuess(i,2)) then
      alphaEqualsBeta = .false.
    end if
  end do
  if(.not. alphaEqualsBeta) then
      print *, ""
      print *, "Warning: Alpha occupation does not match Beta Occupation."
      print *, "Proceeding with alpha occupations..."
      print *, ""
  end if

  print *, ""
  print *, "Ground state orbital occupations for MOM-guess:"
  print *, occupationsGuess(:,1)
  print *, ""
  
  cGuess = c
  occupations(:) = occupationsGuess(:,1)
  P = 2 * matmul(c(:,occupations(1:nO)),transpose(c(:,occupations(1:nO))))

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
    
    call Hartree_matrix_AO_basis(nBas,P,ERI,J)
    call exchange_matrix_AO_basis(nBas,P,ERI,K)
    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)

    ! Check convergence 

    err = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)
    if(nSCF > 2) Conv = maxval(abs(err))

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

    call MOM_density_matrix(nBas, nOrb, nO, S, c, cGuess, &
                              occupations, occupationsGuess(:,1), P)
    P  = 2 * P
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

    write(*,*) ' Warning! Convergence failed at Hartree-Fock level.'

  end if

  call print_RHF(nBas,nOrb,nO,eHF,c,ENuc,ET,EV,EJ,EK,ERHF,dipole)

! Print the 1-RDM and 2-RDM in AO basis
  if(doaordm) then
   call print_RHF_AO_rdms(nBas,ENuc,S,T,V,P,ERI)
  endif

! Testing zone

  if(dotest) then
 
    call dump_test_value('R','RHF energy',ERHF)
    call dump_test_value('R','RHF HOMO energy',eHF(nO))
    call dump_test_value('R','RHF LUMO energy',eHF(nO+1))
    call dump_test_value('R','RHF dipole moment',norm2(dipole))

  end if

! Memory deallocation

  deallocate(J,K,err,cp,Fp,err_diis,F_diis)

end subroutine 

