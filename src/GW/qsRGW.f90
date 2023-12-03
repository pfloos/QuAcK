subroutine qsRGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE, & 
                 singlet,triplet,eta,regularize,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,ERHF,S,X,T,V,Hc,ERI_AO,     & 
                 ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)

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
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: X(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_MO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(inout):: dipole_int_MO(nBas,nBas,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: ispin
  integer                       :: ixyz
  integer                       :: n_diis
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: EK
  double precision              :: EqsGW
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcGM
  double precision              :: Conv
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision              :: dipole(ncart)

  logical                       :: dRPA = .true.
  logical                       :: print_W = .false.
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: eGW(:)
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
  write(*,*)'*******************************'
  write(*,*)'* Restricted qsGW Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsGW !!'
  write(*,*)

! Stuff 

  nBasSq = nBas*nBas

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamical screening!'
    write(*,*)
  end if

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  allocate(eGW(nBas),c(nBas,nBas),cp(nBas,nBas),P(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas), &
           J(nBas,nBas),K(nBas,nBas),SigC(nBas,nBas),SigCp(nBas,nBas),Z(nBas),           & 
           Aph(nS,nS),Bph(nS,nS),Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nBas,nBas,nS),         &
           err(nBas,nBas),err_diis(nBasSq,max_diis),F_diis(nBasSq,max_diis))

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

    call Hartree_matrix_AO_basis(nBas,P,ERI_AO,J)
    call exchange_matrix_AO_basis(nBas,P,ERI_AO,K)

    ! AO to MO transformation of two-electron integrals

    do ixyz=1,ncart
      call AOtoMO(nBas,c,dipole_int_AO(:,:,ixyz),dipole_int_MO(:,:,ixyz))
    end do

    call AOtoMO_ERI_RHF(nBas,c,ERI_AO,ERI_MO)

    ! Compute linear response

    call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI_MO,Aph)
    if(.not.TDA_W) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)

    call phLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
    if(print_W) call print_excitation_energies('phRPA@GW@RHF','singlet',nS,Om)

    ! Compute correlation part of the self-energy 

    call GW_excitation_density(nBas,nC,nO,nR,nS,ERI_MO,XpY,rho)

    if(regularize) call GW_regularization(nBas,nC,nO,nV,nR,nS,eGW,Om,rho)

    call GW_self_energy(eta,nBas,nC,nO,nV,nR,nS,eGW,Om,rho,EcGM,SigC,Z)

    ! Make correlation self-energy Hermitian and transform it back to AO basis
   
    SigC = 0.5d0*(SigC + transpose(SigC))

    call MOtoAO(nBas,S,c,SigC,SigCp)
 
    ! Solve the quasi-particle equation

    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + SigCp(:,:)

    ! Compute commutator and convergence criteria

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
      call DIIS_extrapolation(rcond,nBasSq,nBasSq,n_diis,err_diis,F_diis,err,F)

    end if

    ! Diagonalize Hamiltonian in AO basis

    Fp = matmul(transpose(X),matmul(F,X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp,eGW)
    c = matmul(X,cp)
    call AOtoMO(nBas,c,SigCp,SigC)

    ! Density matrix

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

    ! Print results

    call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsRGW(nBas,nO,nSCF,Conv,thresh,eHF,eGW,c,SigCp,Z,ENuc,ET,EV,EJ,EK,EcGM,EcRPA,EqsGW,dipole)

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

    stop

  end if

! Deallocate memory

  deallocate(c,cp,P,F,Fp,J,K,SigC,SigCp,Z,Om,XpY,XmY,rho,err,err_diis,F_diis)

! Perform BSE calculation

  if(dophBSE) then

    call GW_phBSE(dophBSE2,TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO,eGW,eGW,EcBSE)

    if(exchange_kernel) then

      EcBSE(1) = 0.5d0*EcBSE(1)
      EcBSE(2) = 1.5d0*EcBSE(2)

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@RHF correlation energy           = ',sum(EcBSE),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW@RHF total       energy           = ',ENuc + EqsGW + sum(EcBSE),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '------------------------------------------------------'
      write(*,*) 'Adiabatic connection version of BSE correlation energy'
      write(*,*) '------------------------------------------------------'
      write(*,*)

      if(doXBS) then

        write(*,*) '*** scaled screening version (XBS) ***'
        write(*,*)

      end if

      call GW_phACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,dophBSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,eGW,eGW,EcBSE)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@RHF correlation energy           = ',sum(EcBSE),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@qsGW@RHF total       energy           = ',ENuc + EqsGW + sum(EcBSE),' au'
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

  if(doppBSE) then
      
    call GW_ppBSE(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO,eHF,eGW,EcBSE)
 
    EcBSE(2) = 3d0*EcBSE(2)
   
    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGW@RHF correlation energy           = ',sum(EcBSE),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@qsGW@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Testing zone

  if(dotest) then

    call dump_test_value('R','qsGW correlation energy',EcRPA)
    call dump_test_value('R','qsGW HOMO energy',eGW(nO))
    call dump_test_value('R','qsGW LUMO energy',eGW(nO+1))

  end if

end subroutine 
