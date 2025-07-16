subroutine qsRGWim(dotest,maxSCF,thresh,max_diis,level_shift,eta,shift,nNuc,ZNuc,rNuc,ENuc, & 
                  nBas,nOrb,nO,verbose,S,T,V,Hc,ERI,dipole_int,X,EqsGW,eqsGW_state,c,P,F,   &
                  nfreqs,ntimes,wcoord,wweight)

! Perform restricted Hartree-Fock calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: verbose
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: level_shift

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: eta
  double precision,intent(in)   :: shift
  integer,intent(in)            :: nfreqs
  integer,intent(in)            :: ntimes
  double precision,intent(in)   :: wcoord(nfreqs),wweight(nfreqs)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas) 
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: offdiag0
  logical                       :: file_exists
  integer                       :: iorb,jorb,korb,lorb
  integer                       :: nSCF
  integer                       :: nBas_Sq
  integer                       :: n_diis
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: EK
  double precision              :: dipole(ncart)

  double precision              :: Conv
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision,allocatable  :: Sigc_mo(:,:)
  double precision,allocatable  :: Sigc(:,:)
  double precision,allocatable  :: err(:,:)
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: vMAT(:,:)
  double precision,allocatable  :: ERI_MO(:,:,:,:)

! Output variables

  double precision,intent(out)  :: EqsGW
  double precision,intent(out)  :: eqsGW_state(nOrb)
  double precision,intent(inout):: c(nBas,nOrb)
  double precision,intent(out)  :: P(nBas,nBas)
  double precision,intent(out)  :: F(nBas,nBas)

! Hello world

  write(*,*)
  write(*,*)'****************************************************'
  write(*,*)'* Restricted qsGW Calculation (using imag. freqs.) *'
  write(*,*)'****************************************************'
  write(*,*)

! Useful quantities
  nBas_Sq = nBas*nBas
  offdiag0=.false. ! Set it to true if you want to try qsGW version 2 where all the off-diag. elements of Sigc_mo are eval at Fermi level

! Memory allocation

  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(Sigc(nBas,nBas))

  allocate(err(nBas,nBas))

  allocate(Sigc_mo(nOrb,nOrb))
  allocate(cp(nOrb,nOrb))
  allocate(Fp(nOrb,nOrb))

  allocate(vMAT(nOrb*nOrb,nOrb*nOrb))

  allocate(err_diis(nBas_Sq,max_diis))
  allocate(F_diis(nBas_Sq,max_diis))

! Guess coefficients and density matrix

  P(:,:) = 2d0 * matmul(c(:,1:nO), transpose(c(:,1:nO)))
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
            '|','#','|','E(qsGW)','|','EJ(qsGW)','|','EK(qsGW)','|','Conv','|'
  write(*,*)'-----------------------------------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Increment 

    nSCF = nSCF + 1

    ! Build Fock matrix
    
    call Hartree_matrix_AO_basis(nBas,P,ERI,J)
    call exchange_matrix_AO_basis(nBas,P,ERI,K)
    allocate(ERI_MO(nOrb,nOrb,nOrb,nOrb))
    call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI,ERI_MO)
    do iorb=1,nOrb
     do jorb=1,nOrb
      do korb=1,nOrb
       do lorb=1,nOrb
        vMAT(1+(korb-1)+(iorb-1)*nOrb,1+(lorb-1)+(jorb-1)*nOrb)=ERI_MO(iorb,jorb,korb,lorb)
       enddo
      enddo
     enddo
    enddo
    deallocate(ERI_MO)
    call sigc_MO_basis_RHF(nOrb,nO,offdiag0,eta,shift,eqsGW_state,vMAT,nfreqs, &
                           ntimes,wcoord,wweight,Sigc_mo)
    Sigc_mo = 0.5d0 * (Sigc_mo + transpose(Sigc_mo))
    if(verbose/=0) then    
     write(*,*) 'Sigma_c MO'
     do iorb=1,nOrb
      write(*,'(*(f10.5))') Sigc_mo(iorb,:)
     enddo
    endif
    call MOtoAO(nBas,nOrb,S,c,Sigc_mo,Sigc)
    if(verbose/=0) then
     write(*,*) 'Sigma_c AO'
     do iorb=1,nBas
      write(*,'(*(f10.5))') Sigc(iorb,:)
     enddo
    endif


    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + Sigc(:,:)

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

    EqsGW = ET + EV + EJ + EK

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
      call diagonalize_matrix(nOrb,cp,eqsGW_state)
      c = matmul(X,cp)

    ! Density matrix

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))
!   call dgemm('N', 'T', nBas, nBas, nO, 2.d0, &
!              c(1,1), nBas, c(1,1), nBas,     &
!              0.d0, P(1,1), nBas)

    ! Dump results

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,F16.10,1X,A1,1X,E10.2,1X,A1,1X)') &
      '|',nSCF,'|',EqsGW + ENuc,'|',EJ,'|',EK,'|',Conv,'|'

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

    deallocate(J,K,Sigc,Sigc_mo,err,cp,Fp,err_diis,F_diis)
    deallocate(vMAT)

    write(*,*) ' Warning! Convergence failed at Hartree-Fock level.'

  end if

! Compute dipole moments

  call Hartree_matrix_AO_basis(nBas,P,ERI,J)
  call exchange_matrix_AO_basis(nBas,P,ERI,K)
  F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + Sigc(:,:)
  Fp = matmul(transpose(X),matmul(F,X))
  cp(:,:) = Fp(:,:)
  call diagonalize_matrix(nOrb,cp,eqsGW_state)
  c = matmul(X,cp)

  call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int,dipole)
  call print_RqsGWi(nBas,nOrb,nO,eqsGW_state,c,ENuc,ET,EV,EJ,EK,EqsGW,dipole)

! Testing zone

  if(dotest) then
 
    call dump_test_value('R','qsGW energy',EqsGW)
    call dump_test_value('R','qsGW HOMO energy',eqsGW_state(nO))
    call dump_test_value('R','qsGW LUMO energy',eqsGW_state(nO+1))
    call dump_test_value('R','qsGW dipole moment',norm2(dipole))

  end if

! Memory deallocation

  deallocate(J,K,Sigc,Sigc_mo,err,cp,Fp,err_diis,F_diis)
  deallocate(vMAT)

end subroutine 
