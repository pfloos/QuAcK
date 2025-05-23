subroutine complex_qsRGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2, &
                 TDA_W,TDA,dBSE,dTDA,doppBSE,singlet,triplet,eta,doSRG,nNuc,ZNuc,rNuc,         &
                 ENuc,nBas,nOrb,nC,nO,nV,nR,nS,ERHF,S,X,T,V,Hc,ERI_AO,                         &
                 ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF,                               &
                 CAP_AO,CAP_MO)

! Perform a quasiparticle self-consistent GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  complex*16,intent(in)         :: ERHF
  complex*16,intent(in)         :: eHF(nOrb)
  complex*16,intent(in)         :: cHF(nBas,nOrb)
  complex*16,intent(in)         :: PHF(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: CAP_AO(nBas,nBas)
  complex*16,intent(inout)      :: CAP_MO(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  complex*16,intent(inout)      :: ERI_MO(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  complex*16,intent(inout)      :: dipole_int_MO(nOrb,nOrb,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: nBas_Sq
  integer                       :: ispin
  integer                       :: ixyz
  integer                       :: n_diis
  complex*16                    :: ET
  complex*16                    :: EV
  complex*16                    :: EJ
  complex*16                    :: EK
  complex*16                    :: EqsGW
  complex*16                    :: EW
  complex*16                    :: EcRPA
  complex*16                    :: EcBSE(nspin)
  complex*16                    :: EcGM
  double precision              :: Conv
  double precision              :: rcond
  complex*16,external           :: complex_trace_matrix
  complex*16                    :: dipole(ncart)

  double precision              :: flow

  logical                       :: plot_self = .false.
  logical                       :: dRPA_W  = .true.
  logical                       :: print_W = .false.
  complex*16,allocatable        :: err_diis(:,:)
  complex*16,allocatable        :: F_diis(:,:)
  complex*16,allocatable        :: Aph(:,:)
  complex*16,allocatable        :: Bph(:,:)
  complex*16,allocatable        :: Om(:)
  complex*16,allocatable        :: XpY(:,:)
  complex*16,allocatable        :: XmY(:,:)
  complex*16,allocatable        :: rho(:,:,:)
  complex*16,allocatable        :: c(:,:)
  complex*16,allocatable        :: cp(:,:)
  complex*16,allocatable        :: eGW(:)
  complex*16,allocatable        :: P(:,:)
  complex*16,allocatable        :: F(:,:)
  complex*16,allocatable        :: Fp(:,:)
  complex*16,allocatable        :: J(:,:)
  complex*16,allocatable        :: K(:,:)
  complex*16,allocatable        :: SigC(:,:)
  complex*16,allocatable        :: SigCp(:,:)
  complex*16,allocatable        :: Z(:)
  complex*16,allocatable        :: err(:,:)

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted qsGW Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs and CAP in MO basis will be overwritten in qsGW !!'
  write(*,*)

! Stuff 

  nBas_Sq = nBas*nBas

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamical screening!'
    write(*,*)
  end if

! SRG regularization

  flow = 500d0

  if(doSRG) then

    write(*,*) '*** SRG regularized qsGW scheme ***'
    write(*,*)

  end if

! Memory allocation

  allocate(eGW(nOrb))
  allocate(Z(nOrb))

  allocate(c(nBas,nOrb))

  allocate(cp(nOrb,nOrb))
  allocate(Fp(nOrb,nOrb))
  allocate(SigC(nOrb,nOrb))

  allocate(P(nBas,nBas))
  allocate(F(nBas,nBas))
  allocate(J(nBas,nBas))
  allocate(K(nBas,nBas))
  allocate(err(nBas,nBas))
  allocate(SigCp(nBas,nBas))

  allocate(Aph(nS,nS))
  allocate(Bph(nS,nS))
  allocate(Om(nS))
  allocate(XpY(nS,nS))
  allocate(XmY(nS,nS))
  allocate(rho(nOrb,nOrb,nS))

  allocate(err_diis(nBas_Sq,max_diis))
  allocate(F_diis(nBas_Sq,max_diis))

! Initialization
  
  nSCF          = -1
  n_diis        = 0
  ispin         = 1
  Conv          = 1d0
  P(:,:)        = PHF(:,:)
  eGW(:)        = eHF(:)
  c(:,:)        = cHF(:,:)
  F_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0
  rcond          = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

    ! Increment

    nSCF = nSCF + 1

    ! Build Hartree-exchange matrix

    call complex_Hartree_matrix_AO_basis(nBas,P,ERI_AO,J)
    call complex_exchange_matrix_AO_basis(nBas,P,ERI_AO,K)

    ! AO to MO transformation of two-electron integrals

    do ixyz=1,ncart
      call complex_AOtoMO(nBas,nOrb,c,dipole_int_AO(1,1,ixyz),dipole_int_MO(1,1,ixyz))
    end do

    call complex_AOtoMO_ERI_RHF(nBas,nOrb,c,ERI_AO,ERI_MO)

    ! Compute linear response

                   call complex_phRLR_A(ispin,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,eGW,ERI_MO,Aph)
    if(.not.TDA_W) call complex_phRLR_B(ispin,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)

    call complex_phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
    if(print_W) call print_excitation_energies('phRPA@GW@RHF','singlet',nS,Om)

    call complex_RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI_MO,XpY,rho)
    
    if(doSRG) then 
      call complex_RGW_SRG_self_energy(flow,eta,nBas,nOrb,nC,nO,nV,nR,nS,eGW,Om,&
                                   rho,EcGM,SigC,Z)
    else
      call complex_RGW_self_energy(eta,nBas,nOrb,nC,nO,nV,nR,nS,eGW,Om,rho,&
                                   EcGM,SigC,Z)
    end if
    ! Make correlation self-energy Hermitian and transform it back to AO basis
   
    SigC = 0.5d0*(SigC + transpose(SigC))

    call complex_MOtoAO(nBas,nOrb,S,c,SigC,SigCp)
 
    ! Solve the quasi-particle equation

    F(:,:) = cmplx(Hc(:,:),CAP_AO(:,:),kind=8) + J(:,:) + 0.5d0*K(:,:) + SigCp(:,:)
    if(nBas .ne. nOrb) then
      call complex_complex_AOtoMO(nBas,nOrb,c(1,1),F(1,1),Fp(1,1))
      call complex_MOtoAO(nBas,nOrb,S(1,1),c(1,1),Fp(1,1),F(1,1))
    endif
    
    ! Compute commutator and convergence criteria

    err = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)

    if(nSCF > 1) Conv = maxval(abs(err))

    ! Kinetic energy
    
    ET = complex_trace_matrix(nBas,matmul(P,T))

    ! Potential energy

    EV = complex_trace_matrix(nBas,matmul(P,V))

    ! Hartree energy

    EJ = 0.5d0*complex_trace_matrix(nBas,matmul(P,J))

    ! Exchange energy

    EK = 0.25d0*complex_trace_matrix(nBas,matmul(P,K))
    
    ! CAP energy

    EW = complex_trace_matrix(nBas,matmul(P,(0d0,1d0)*CAP_AO))
   
    ! Total energy

    EqsGW = ET + EV + EJ + EK + EW 

    ! DIIS extrapolation 

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call complex_DIIS_extrapolation(rcond,nBas_Sq,nBas_Sq,n_diis,err_diis,F_diis,err,F)

    end if

    ! Diagonalize Hamiltonian in AO basis

    if(nBas .eq. nOrb) then
      Fp = matmul(transpose(X),matmul(F,X))
      cp(:,:) = Fp(:,:)
      call complex_diagonalize_matrix(nOrb,cp,eGW)
      call complex_orthogonalize_matrix(nOrb,cp)
      c = matmul(X,cp)
    else
      Fp = matmul(transpose(c),matmul(F,c))
      cp(:,:) = Fp(:,:)
      call complex_diagonalize_matrix(nOrb,cp,eGW)
      call complex_orthogonalize_matrix(nOrb,cp)
      c = matmul(c,cp)
    endif

    call complex_complex_AOtoMO(nBas,nOrb,c,SigCp,SigC)
    
    ! Density matrix

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

    ! Print results

    !call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_complex_qsRGW(nBas,nOrb,nO,nSCF,Conv,thresh,eHF,eGW,c,SigC,Z, &
                     ENuc,ET,EV,EW,EJ,EK,EcGM,EcRPA,EqsGW,dipole)

  end do
  
! Plot self-energy, renormalization factor, and spectral function
!
  if(plot_self) call complex_RGW_plot_self_energy(nOrb,eta,nC,nO,nV,nR,nS,real(eHF),aimag(eHF),real(eGW),aimag(eGW),Om,rho)
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Did it actually converge?

!  if(nSCF == maxSCF+1) then
!
!    write(*,*)
!    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(*,*)'                 Convergence failed                 '
!    write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!    write(*,*)
!
!    deallocate(c,cp,P,F,Fp,J,K,SigC,SigCp,Z,Om,XpY,XmY,rho,err,err_diis,F_diis)
!    stop
!
!  end if
!
!! Deallocate memory
!
!  deallocate(c,cp,P,F,Fp,J,K,SigC,SigCp,Z,Om,XpY,XmY,rho,err,err_diis,F_diis)
!  
!! Perform BSE calculation
!
!  if(dophBSE) then
!
!    call RGW_phBSE(dophBSE2,exchange_kernel,TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta, &
!                   nOrb,nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO,eGW,eGW,EcBSE)
!
!    write(*,*)
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@RHF correlation energy           = ',sum(EcBSE),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@RHF total       energy           = ',ENuc + EqsGW + sum(EcBSE),' au'
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,*)
!
!!   Compute the BSE correlation energy via the adiabatic connection 
!
!    if(doACFDT) then
!
!      call RGW_phACFDT(exchange_kernel,doXBS,TDA_W,TDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI_MO,eGW,eGW,EcBSE)
!
!      write(*,*)
!      write(*,*)'-------------------------------------------------------------------------------'
!      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
!      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
!      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@RHF correlation energy           = ',sum(EcBSE),' au'
!      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@RHF total       energy           = ',ENuc + EqsGW + sum(EcBSE),' au'
!      write(*,*)'-------------------------------------------------------------------------------'
!      write(*,*)
!
!    end if
!
!  end if
!
!  if(doppBSE) then
!      
!    call RGW_ppBSE(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO,eHF,eGW,EcBSE)
!   
!    write(*,*)
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGW@RHF correlation energy           = ',sum(EcBSE),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGW@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,*)
!
!  end if
!
!! Testing zone
!
!  if(dotest) then
!
!    call dump_test_value('R','qsGW correlation energy',EcRPA)
!    call dump_test_value('R','qsGW HOMO energy',eGW(nO))
!    call dump_test_value('R','qsGW LUMO energy',eGW(nO+1))
!
!  end if
!
end subroutine 
