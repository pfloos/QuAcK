subroutine qsGGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE, & 
                 eta,regularize,nNuc,ZNuc,rNuc,ENuc,nBas,nBas2,nC,nO,nV,nR,nS,EGHF,Ov,Or,T,V,Hc,ERI_AO,  & 
                 ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)

! Generalized version of quasiparticle self-consistent GW 

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
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nNuc
  double precision,intent(in)   :: ZNuc(nNuc)
  double precision,intent(in)   :: rNuc(nNuc,ncart)
  double precision,intent(in)   :: ENuc

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nBas2
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: eHF(nBas2)
  double precision,intent(in)   :: cHF(nBas2,nBas2)
  double precision,intent(in)   :: PHF(nBas2,nBas2)
  double precision,intent(in)   :: Ov(nBas,nBas)
  double precision,intent(in)   :: T(nBas,nBas)
  double precision,intent(in)   :: V(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: Or(nBas,nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(inout):: ERI_MO(nBas2,nBas2,nBas2,nBas2)
  double precision,intent(in)   :: dipole_int_AO(nBas,nBas,ncart)
  double precision,intent(inout):: dipole_int_MO(nBas2,nBas2,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: nBas2Sq
  integer                       :: ispin
  integer                       :: ixyz
  integer                       :: n_diis
  double precision              :: ET,ETaa,ETbb
  double precision              :: EV,EVaa,EVbb
  double precision              :: EJ,EJaaaa,EJaabb,EJbbaa,EJbbbb
  double precision              :: EK,EKaaaa,EKabba,EKbaab,EKbbbb
  double precision              :: EqsGW
  double precision              :: EcRPA
  double precision              :: EcBSE
  double precision              :: EcGM
  double precision              :: Conv
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision              :: dipole(ncart)

  logical                       :: dRPA = .true.
  logical                       :: print_W = .true.
  double precision,allocatable  :: err_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)

  double precision,allocatable  :: Ca(:,:),Cb(:,:)
  double precision,allocatable  :: ERI_tmp(:,:,:,:)
  double precision,allocatable  :: Jaa(:,:),Jbb(:,:)
  double precision,allocatable  :: Kaa(:,:),Kab(:,:),Kba(:,:),Kbb(:,:)
  double precision,allocatable  :: Faa(:,:),Fab(:,:),Fba(:,:),Fbb(:,:)
  double precision,allocatable  :: Paa(:,:),Pab(:,:),Pba(:,:),Pbb(:,:)

  double precision,allocatable  :: C(:,:)
  double precision,allocatable  :: Cp(:,:)
  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: H(:,:)
  double precision,allocatable  :: S(:,:)
  double precision,allocatable  :: X(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: SigCp(:,:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: err(:,:)

! Hello world

  write(*,*)
  write(*,*)'********************************'
  write(*,*)'* Generalized qsGW Calculation *'
  write(*,*)'********************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsGW !!'
  write(*,*)

! Stuff 

  nBasSq = nBas*nBas
  nBas2Sq = nBas2*nBas2

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
    write(*,*)
  end if

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  allocate(Ca(nBas,nBas2),Cb(nBas,nBas2),P(nBas2,nBas2),Jaa(nBas,nBas),Jbb(nBas,nBas), &
           Kaa(nBas,nBas),Kab(nBas,nBas),Kba(nBas,nBas),Kbb(nBas,nBas),   &
           Faa(nBas,nBas),Fab(nBas,nBas),Fba(nBas,nBas),Fbb(nBas,nBas),   &
           Paa(nBas,nBas),Pab(nBas,nBas),Pba(nBas,nBas),Pbb(nBas,nBas),   &
           F(nBas2,nBas2),Fp(nBas2,nBas2),C(nBas2,nBas2),Cp(nBas2,nBas2), &
           H(nBas2,nBas2),S(nBas2,nBas2),X(nBas2,nBas2),err(nBas2,nBas2), &
           err_diis(nBas2Sq,max_diis),F_diis(nBas2Sq,max_diis),           &
           eGW(nBas2),SigC(nBas2,nBas2),SigCp(nBas,nBas),Z(nBas2),Aph(nS,nS),Bph(nS,nS),   & 
           Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nBas2,nBas2,nS))

! Initialization
  
  nSCF          = -1
  n_diis        = 0
  ispin         = 3
  Conv          = 1d0
  P(:,:)        = PHF(:,:)
  eGW(:)        = eHF(:)
  c(:,:)        = cHF(:,:)
  F_diis(:,:)   = 0d0
  err_diis(:,:) = 0d0
  rcond         = 0d0

! Construct super overlap matrix

  S(      :     ,       :    ) = 0d0
  S(     1:nBas ,     1:nBas ) = Ov(1:nBas,1:nBas)
  S(nBas+1:nBas2,nBas+1:nBas2) = Ov(1:nBas,1:nBas)

! Construct super orthogonalization matrix

  X(      :     ,       :    ) = 0d0
  X(     1:nBas ,     1:nBas ) = Or(1:nBas,1:nBas)
  X(nBas+1:nBas2,nBas+1:nBas2) = Or(1:nBas,1:nBas)

! Construct super orthogonalization matrix

  H(      :     ,       :    ) = 0d0
  H(     1:nBas ,     1:nBas ) = Hc(1:nBas,1:nBas)
  H(nBas+1:nBas2,nBas+1:nBas2) = Hc(1:nBas,1:nBas)

! Construct super density matrix

  P(:,:) = matmul(C(:,1:nO),transpose(C(:,1:nO)))

  Paa(:,:) = P(     1:nBas ,     1:nBas )
  Pab(:,:) = P(     1:nBas ,nBas+1:nBas2)
  Pba(:,:) = P(nBas+1:nBas2,     1:nBas )
  Pbb(:,:) = P(nBas+1:nBas2,nBas+1:nBas2)

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

!   Increment

    nSCF = nSCF + 1

!   Buid Hartree matrix

    call Hartree_matrix_AO_basis(nBas,Paa,ERI_AO,Jaa)
    call Hartree_matrix_AO_basis(nBas,Pbb,ERI_AO,Jbb)

!   Compute exchange part of the self-energy 

    call exchange_matrix_AO_basis(nBas,Paa,ERI_AO,Kaa)
    call exchange_matrix_AO_basis(nBas,Pba,ERI_AO,Kab)
    call exchange_matrix_AO_basis(nBas,Pab,ERI_AO,Kba)
    call exchange_matrix_AO_basis(nBas,Pbb,ERI_AO,Kbb)

!   Build individual Fock matrices

    Faa(:,:) = Hc(:,:) + Jaa(:,:) + Jbb(:,:) + Kaa(:,:) 
    Fab(:,:) =                               + Kab(:,:)
    Fba(:,:) =                               + Kba(:,:)
    Fbb(:,:) = Hc(:,:) + Jbb(:,:) + Jaa(:,:) + Kbb(:,:)

!   Build super Fock matrix

    F(     1:nBas ,     1:nBas ) = Faa(1:nBas,1:nBas)
    F(     1:nBas ,nBas+1:nBas2) = Fab(1:nBas,1:nBas)
    F(nBas+1:nBas2,     1:nBas ) = Fba(1:nBas,1:nBas)
    F(nBas+1:nBas2,nBas+1:nBas2) = Fbb(1:nBas,1:nBas)

!   AO to MO transformation of two-electron integrals

    allocate(ERI_tmp(nBas2,nBas2,nBas2,nBas2))
 
    Ca(:,:) = C(1:nBas,1:nBas2)
    Cb(:,:) = C(nBas+1:nBas2,1:nBas2)

    do ixyz=1,ncart
        call AOtoMO_GHF(nBas,nBas2,Ca,Cb,dipole_int_AO(:,:,ixyz),dipole_int_MO(:,:,ixyz))
    end do

    call AOtoMO_ERI_GHF(nBas,nBas2,Ca,Ca,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_tmp(:,:,:,:)

    call AOtoMO_ERI_GHF(nBas,nBas2,Ca,Cb,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

    call AOtoMO_ERI_GHF(nBas,nBas2,Cb,Ca,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

    call AOtoMO_ERI_GHF(nBas,nBas2,Cb,Cb,ERI_AO,ERI_tmp)
    ERI_MO(:,:,:,:) = ERI_MO(:,:,:,:) + ERI_tmp(:,:,:,:)

    deallocate(ERI_tmp)

!   Compute linear response

    call phLR_A(ispin,dRPA,nBas2,nC,nO,nV,nR,nS,1d0,eGW,ERI_MO,Aph)
    if(.not.TDA_W) call phLR_B(ispin,dRPA,nBas2,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)

    call phLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)
    if(print_W) call print_excitation_energies('phRPA@GGW','spinorbital',nS,Om)

!   Compute correlation part of the self-energy 

    call GW_excitation_density(nBas2,nC,nO,nR,nS,ERI_MO,XpY,rho)

    if(regularize) call GW_regularization(nBas2,nC,nO,nV,nR,nS,eGW,Om,rho)

    call GGW_self_energy(eta,nBas2,nC,nO,nV,nR,nS,eGW,Om,rho,EcGM,SigC,Z)

!   Make correlation self-energy Hermitian and transform it back to AO basis
   
    SigC = 0.5d0*(SigC + transpose(SigC))

    call MOtoAO_GHF(nBas2,nBas,S,Ca,Cb,SigC,SigCp)

!   ... and add self-energy 

    F(:,:) = F(:,:) + SigCp(:,:)

!   Compute commutator and convergence criteria

    err = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)

    if(nSCF > 1) Conv = maxval(abs(err))

!   DIIS extrapolation 

    if(max_diis > 1) then
      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBas2Sq,nBas2Sq,n_diis,err_diis,F_diis,err,F)
    end if

!  Transform Fock matrix in orthogonal basis

    Fp(:,:) = matmul(transpose(X),matmul(F,X))

!  Diagonalize Fock matrix to get eigenvectors and eigenvalues

    Cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas2,Cp,eGW)

!   Back-transform eigenvectors in non-orthogonal basis

    C(:,:) = matmul(X,Cp)

    call AOtoMO_GHF(nBas,nBas2,Ca,Cb,SigCp,SigC)

!   Form super density matrix

    P(:,:) = matmul(C(:,1:nO),transpose(C(:,1:nO)))

!   Compute individual density matrices

    Paa(:,:) = P(     1:nBas ,     1:nBas )
    Pab(:,:) = P(     1:nBas ,nBas+1:nBas2)
    Pba(:,:) = P(nBas+1:nBas2,     1:nBas )
    Pbb(:,:) = P(nBas+1:nBas2,nBas+1:nBas2)

    !------------------------------------------------------------------------
    !   Compute total energy
    !------------------------------------------------------------------------

!  Kinetic energy

    ETaa = trace_matrix(nBas,matmul(Paa,T))
    ETbb = trace_matrix(nBas,matmul(Pbb,T))

    ET = ETaa + ETbb

!  Potential energy

   EVaa = trace_matrix(nBas,matmul(Paa,V))
   EVbb = trace_matrix(nBas,matmul(Pbb,V))

   EV = EVaa + EVbb

!  Hartree energy

    EJaaaa = 0.5d0*trace_matrix(nBas,matmul(Paa,Jaa))
    EJaabb = 0.5d0*trace_matrix(nBas,matmul(Paa,Jbb))
    EJbbaa = 0.5d0*trace_matrix(nBas,matmul(Pbb,Jaa))
    EJbbbb = 0.5d0*trace_matrix(nBas,matmul(Pbb,Jbb))

    EJ = EJaaaa + EJaabb + EJbbaa + EJbbbb

!   Exchange energy

    EKaaaa = 0.5d0*trace_matrix(nBas,matmul(Paa,Kaa))
    EKabba = 0.5d0*trace_matrix(nBas,matmul(Pab,Kba))
    EKbaab = 0.5d0*trace_matrix(nBas,matmul(Pba,Kab))
    EKbbbb = 0.5d0*trace_matrix(nBas,matmul(Pbb,Kbb))

    EK = EKaaaa + EKabba + EKbaab + EKbbbb

!   Total energy

    EqsGW = ET + EV + EJ + EK 

    ! Print results

    call dipole_moment(nBas2,P,nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsGGW(nBas,nBas2,nO,nSCF,Conv,thresh,eHF,eGW,c,Ov,SigC,Z,ENuc,ET,EV,EJ,EK,EcGM,EcRPA,EqsGW,dipole)

  enddo
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

  endif

! Perform BSE calculation

  if(dophBSE) then

    call GGW_phBSE(dophBSE2,TDA_W,TDA,dBSE,dTDA,eta,nBas2,nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO,eGW,eGW,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@qsGW correlation energy           =',EcBSE,' au'
    write(*,'(2X,A50,F20.10,A4)') 'Tr@BSE@qsGW total energy                 =',ENuc + EqsGW + EcBSE,' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

!   if(doACFDT) then

!     write(*,*) '------------------------------------------------------'
!     write(*,*) 'Adiabatic connection version of BSE correlation energy'
!     write(*,*) '------------------------------------------------------'
!     write(*,*)

!     if(doXBS) then

!       write(*,*) '*** scaled screening version (XBS) ***'
!       write(*,*)

!     end if

!     call GW_phACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,dophBSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,eGW,eGW,EcBSE)

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy (singlet) =',EcBSE(1)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy (triplet) =',EcBSE(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy           =',EcBSE(1) + EcBSE(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW total energy                 =',ENuc + EqsGW + EcBSE(1) + EcBSE(2)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

  end if

! if(doppBSE) then
!     
!   call GW_ppBSE(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO,eHF,eGW,EcBSE)
!   
!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@qsGW correlation energy (singlet) =',EcBSE(1)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@qsGW correlation energy (triplet) =',3d0*EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@qsGW correlation energy =',EcBSE(1) + 3d0*EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@qsGW total energy =',ENuc + EGHF + EcBSE(1) + 3d0*EcBSE(2)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

! end if

! Testing zone

  if(dotest) then

    call dump_test_value('G','qsGW correlation energy',EcRPA)
    call dump_test_value('G','qsGW HOMO energy',eGW(nO))
    call dump_test_value('G','qsGW LUMO energy',eGW(nO+1))

  end if

end subroutine 
