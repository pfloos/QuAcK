subroutine qsRGF2(dotest,maxSCF,thresh,max_diis,dophBSE,doppBSE,TDA,  &
                  dBSE,dTDA,singlet,triplet,eta,doSRG,nNuc,ZNuc,      &
                  rNuc,ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,S,X,T,V,Hc, & 
                  ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)

! Perform a quasiparticle self-consistent GF2 calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas,nOrb,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: cHF(nBas,nOrb)
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nOrb)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_MO(nOrb,nOrb,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: nBas_Sq
  integer                       :: ispin
  integer                       :: n_diis
  double precision              :: EqsGF2
  double precision              :: Conv
  double precision              :: rcond
  double precision              :: flow
  double precision,external     :: trace_matrix
  double precision              :: dipole(ncart)
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: Ex
  double precision              :: Ec
  double precision              :: EcBSE(nspin)

  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: eGF(:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: SigCp(:,:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: err(:,:)

! Hello world


  write(*,*)
  write(*,*)'********************************'
  write(*,*)'* Restricted qsGF2 Calculation *'
  write(*,*)'********************************'
  write(*,*)

! SRG regularization

  flow = 500d0

  if(doSRG) then

    write(*,*) '*** SRG regularized qsGF2 scheme ***'
    write(*,*)

  end if

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsGF2 !!'
  write(*,*)

! Stuff 

  nBas_Sq = nBas*nBas

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  allocate(eGF(nOrb))

  allocate(c(nBas,nOrb))

  allocate(cp(nOrb,nOrb))
  allocate(Fp(nOrb,nOrb))

  allocate(P(nBas,nBas))
  allocate(F(nBas,nBas))
  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(err(nBas,nBas))

  allocate(Z(nOrb))
  allocate(SigC(nOrb,nOrb))

  allocate(SigCp(nBas,nBas))

  allocate(err_diis(nBas_Sq,max_diis))
  allocate(F_diis(nBas_Sq,max_diis))

! Initialization
  
  nSCF          = -1
  n_diis        = 0
  ispin         = 1
  Conv          = 1d0
  P(:,:)        = PHF(:,:)
  eGF(:)        = eHF(:)
  c(:,:)        = cHF(:,:)
  F_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0
  rcond         = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Increment

    nSCF = nSCF + 1

    ! Buid Hartree matrix

    call Hartree_matrix_AO_basis(nBas,P,ERI_AO,J)

    ! Compute exchange part of the self-energy 

    call exchange_matrix_AO_basis(nBas,P,ERI_AO,K)

    ! AO to MO transformation of two-electron integrals

    call AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI_MO)

    ! Compute self-energy and renormalization factor

    if(doSRG) then

      call RGF2_SRG_self_energy(flow,nOrb,nC,nO,nV,nR,eGF,ERI_MO,Ec,SigC,Z)

    else

      call RGF2_self_energy(eta,nOrb,nC,nO,nV,nR,eGF,ERI_MO,Ec,SigC,Z)

    end if

    print*,Ec

    ! Make correlation self-energy Hermitian and transform it back to AO basis
   
    SigC = 0.5d0*(SigC + transpose(SigC))

    call MOtoAO(nBas,nOrb,S,c,SigC,SigCp)
 
    ! Solve the quasi-particle equation

    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + SigCp(:,:)
    if(nBas .ne. nOrb) then
      call AOtoMO(nBas,nOrb,c,F,Fp)
      call MOtoAO(nBas,nOrb,S,c,Fp,F)
    endif

    ! Compute commutator and convergence criteria

    err = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)

    Conv = maxval(abs(err))

    ! Kinetic energy

    ET = trace_matrix(nBas,matmul(P,T))

    ! Potential energy

    EV = trace_matrix(nBas,matmul(P,V))

    ! Hartree energy

    EJ = 0.5d0*trace_matrix(nBas,matmul(P,J))

    ! Exchange energy

    Ex = 0.25d0*trace_matrix(nBas,matmul(P,K))

    ! Total energy

    EqsGF2 = ET + EV + EJ + Ex + Ec

    ! DIIS extrapolation 

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBas_Sq,nBas_Sq,n_diis,err_diis,F_diis,err,F)

    end if

    ! Diagonalize Hamiltonian in AO basis

    if(nBas == nOrb) then
      Fp = matmul(transpose(X),matmul(F,X))
      cp(:,:) = Fp(:,:)
      call diagonalize_matrix(nOrb,cp,eGF)
      c = matmul(X,cp)
    else
      Fp = matmul(transpose(c),matmul(F,c))
      cp(:,:) = Fp(:,:)
      call diagonalize_matrix(nOrb,cp,eGF)
      c = matmul(c,cp)
    endif

    call AOtoMO(nBas,nOrb,c,SigCp,SigC)

    ! Compute new density matrix in the AO basis

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

    ! Print results

    call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsRGF2(nBas,nOrb,nO,nSCF,Conv,thresh,eHF,eGF,&
                      c,SigC,Z,ENuc,ET,EV,EJ,Ex,Ec,EqsGF2,dipole)

  end do
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Did it actually converge?

  if(nSCF == maxSCF+1) then

    write(*,*)
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)'                 Convergence failed                 '
    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    write(*,*)

    deallocate(c,cp,P,F,Fp,J,K,SigC,SigCp,Z,err,err_diis,F_diis)
    stop

  end if

! Deallocate memory

  deallocate(c,cp,P,F,Fp,J,K,SigC,SigCp,Z,err,err_diis,F_diis)

! Perform phBSE@GF2 calculation

  if(dophBSE) then

    call RGF2_phBSE(TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,&
                    nV,nR,nS,ERI_MO,dipole_int_MO,eGF,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@qsGF2 correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@qsGF2 correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@qsGF2 correlation energy           =',sum(EcBSE(:))
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@qsGF2 total energy                 =',ENuc + EqsGF2 + sum(EcBSE(:))
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if


! Perform ppBSE@GF2 calculation

  if(doppBSE) then

    call RGF2_ppBSE(TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,&
                    nC,nO,nV,nR,ERI_MO,dipole_int_MO,eGF,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGF2 correlation energy (singlet) =',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGF2 correlation energy (triplet) =',3d0*EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGF2 correlation energy           =',EcBSE(1) + 3d0*EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGF2 total energy                 =',ENuc + EqsGF2 + EcBSE(1) + 3d0*EcBSE(2),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Testing zone

  if(dotest) then

    call dump_test_value('R','qsGF2 correlation energy',Ec)
    call dump_test_value('R','qsGF2 HOMO energy',eGF(nO))
    call dump_test_value('R','qsGF2 LUMO energy',eGF(nO+1))

  end if

end subroutine 
