subroutine RHF(dotest,doaordm,maxSCF,thresh,max_diis,guess_type,level_shift,nNuc,ZNuc,rNuc,ENuc, & 
               nBas,nOrb,nO,S,T,V,Hc,ERI,dipole_int,X,ERHF,eHF,c,P,F)

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
  double precision              :: trace_1rdm
  double precision              :: trace_2rdm
  double precision,external     :: trace_matrix
  double precision,allocatable  :: err(:,:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: AO_1rdm(:,:)
  double precision,allocatable  :: AO_2rdm(:,:,:,:)

! Output variables

  double precision,intent(out)  :: ERHF
  double precision,intent(out)  :: eHF(nOrb)
  double precision,intent(inout):: c(nBas,nOrb)
  double precision,intent(out)  :: P(nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'*****************************'
  write(*,*)'* Restricted HF Calculation *'
  write(*,*)'*****************************'
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

! Guess coefficients and density matrix

  call mo_guess(nBas,nOrb,guess_type,S,Hc,X,c)
  P(:,:) = 2d0 * matmul(c(:,1:nO), transpose(c(:,1:nO)))
  if(guess_type == 5) then
   inquire(file='P_ao_bin', exist=file_exists)
   if(file_exists) then
    write(*,*) 'Reading P_ao_bin matrix...'
    open(unit=314, form='unformatted', file='P_ao_bin', status='old')
    do
     read(314) ibas,jbas,Val
     if(ibas==0 .and. jbas==0) exit
     P(ibas,jbas)=Val 
    enddo
    close(314)
   endif
  endif

! call dgemm('N', 'T', nBas, nBas, nO, 2.d0, &
!            c(1,1), nBas, c(1,1), nBas,     &
!            0.d0, P(1,1), nBas)

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

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))
!   call dgemm('N', 'T', nBas, nBas, nO, 2.d0, &
!              c(1,1), nBas, c(1,1), nBas,     &
!              0.d0, P(1,1), nBas)

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

    !deallocate(J,K,err,cp,Fp,err_diis,F_diis)

    write(*,*) ' Warning! Convergence failed at Hartree-Fock level.'

  end if

! Compute dipole moments

  call Hartree_matrix_AO_basis(nBas,P,ERI,J)
  call exchange_matrix_AO_basis(nBas,P,ERI,K)
  F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:)
  Fp = matmul(transpose(X),matmul(F,X))
  cp(:,:) = Fp(:,:)
  call diagonalize_matrix(nOrb,cp,eHF)
  c = matmul(X,cp)

  call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_RHF(nBas,nOrb,nO,eHF,c,ENuc,ET,EV,EJ,EK,ERHF,dipole)

! Print the 1-RDM and 2-RDM in AO basis
  if(doaordm) then
   write(*,*)
   write(*,'(a)') ' -------------------------------------------'
   write(*,'(a)') ' Computing and printing RDMs in the AO basis'
   write(*,'(a)') ' -------------------------------------------'
   allocate(AO_1rdm(nBas,nBas),AO_2rdm(nBas,nBas,nBas,nBas))
   AO_1rdm(:,:)=P(:,:)
   AO_2rdm(:,:,:,:)=0d0
   do ibas=1,nBas
    do jbas=1,nBas
     do kbas=1,nBas
      do lbas=1,nBas
       ! Hartree
       AO_2rdm(ibas,jbas,kbas,lbas)=AO_2rdm(ibas,jbas,kbas,lbas)+0.5d0*P(ibas,kbas)*P(jbas,lbas)
       ! Exchange
       AO_2rdm(ibas,jbas,kbas,lbas)=AO_2rdm(ibas,jbas,kbas,lbas)-0.25d0*P(ibas,lbas)*P(jbas,kbas)
      enddo
     enddo
    enddo
   enddo
   trace_1rdm=0d0
   trace_2rdm=0d0
   iunit=312
   iunit2=313
   open(unit=iunit,form='unformatted',file='rhf_ao_1rdm')
   open(unit=iunit2,form='unformatted',file='rhf_ao_2rdm')
   Ecore=0d0; Eee=0d0;
   do ibas=1,nBas
    do jbas=1,nBas
     trace_1rdm=trace_1rdm+AO_1rdm(ibas,jbas)*S(ibas,jbas)
     write(iunit) ibas,jbas,AO_1rdm(ibas,jbas)
     Ecore=Ecore+AO_1rdm(ibas,jbas)*(T(ibas,jbas)+V(ibas,jbas))
     do kbas=1,nBas
      do lbas=1,nBas
       trace_2rdm=trace_2rdm+AO_2rdm(ibas,jbas,kbas,lbas)*S(ibas,kbas)*S(jbas,lbas)
       write(iunit2) ibas,jbas,kbas,lbas,AO_2rdm(ibas,jbas,kbas,lbas)
       Eee=Eee+AO_2rdm(ibas,jbas,kbas,lbas)*ERI(ibas,jbas,kbas,lbas)
      enddo
     enddo
    enddo
   enddo
   write(iunit) 0,0,0d0
   write(iunit2) 0,0,0,0,0d0
   close(iunit)
   close(iunit2)
   deallocate(AO_1rdm,AO_2rdm)
   write(*,'(a)') ' Energies computed using the 1-RDM and the 2-RDM in the AO basis'
   write(*,*)
   write(*,'(a,f17.8)') '   Hcore (T+V) ',Ecore
   write(*,'(a,f17.8)') '     Vee (Hxc) ',Eee
   write(*,'(a,f17.8)') '   Eelectronic ',Ecore+Eee
   write(*,*)           ' --------------'
   write(*,'(a,f17.8)') '        Etotal ',Ecore+Eee+ENuc
   write(*,*)           ' --------------'
   write(*,'(a,f17.8)') '   Tr[ 1D^AO ] ',trace_1rdm
   write(*,'(a,f17.8)') '   Tr[ 2D^AO ] ',trace_2rdm
   write(*,*)
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
