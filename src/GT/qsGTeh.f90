subroutine qsGTeh(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,BSE,BSE2,TDA_T,TDA,          & 
                  dBSE,dTDA,evDyn,singlet,triplet,eta,regularize,nNuc,ZNuc,rNuc,ENuc,nBas,nC,nO,nV,nR,nS,ERHF, &
                  S,X,T,V,Hc,ERI_AO,ERI_MO,dipole_int_AO,dipole_int_MO,PHF,cHF,eHF)

! Perform a quasiparticle self-consistent GTeh calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: BSE
  logical,intent(in)            :: BSE2
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
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
  double precision,intent(in)   :: dipole_int_MO(nBas,nBas,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: nBasSq
  integer                       :: ispin
  integer                       :: n_diis
  double precision              :: ET
  double precision              :: EV
  double precision              :: EJ
  double precision              :: Ex
  double precision              :: EqsGT
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM
  double precision              :: Conv
  double precision              :: rcond
  double precision,external     :: trace_matrix
  double precision              :: dipole(ncart)

  logical                       :: print_T = .true.
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: F_diis(:,:)
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rhoL_RPA(:,:,:)
  double precision,allocatable  :: rhoR_RPA(:,:,:)
  double precision,allocatable  :: c(:,:)
  double precision,allocatable  :: cp(:,:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: P(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: Fp(:,:)
  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: SigCp(:,:)
  double precision,allocatable  :: SigCm(:,:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: error(:,:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|     Self-consistent qsGTeh calculation       |'
  write(*,*)'************************************************'
  write(*,*)

! Warning 

  write(*,*) '!! ERIs in MO basis will be overwritten in qsGTeh !!'
  write(*,*)

! Stuff 

  nBasSq = nBas*nBas

! TDA for T

  if(TDA_T) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
    write(*,*)
  end if

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Memory allocation

  allocate(eGT(nBas),eOld(nBas),c(nBas,nBas),cp(nBas,nBas),P(nBas,nBas),F(nBas,nBas),Fp(nBas,nBas), &
           J(nBas,nBas),K(nBas,nBas),SigC(nBas,nBas),SigCp(nBas,nBas),SigCm(nBas,nBas),Z(nBas),     & 
           OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rhoL_RPA(nBas,nBas,nS),rhoR_RPA(nBas,nBas,nS),   &
           error(nBas,nBas),error_diis(nBasSq,max_diis),F_diis(nBasSq,max_diis))

! Initialization
  
  nSCF            = -1
  n_diis          = 0
  ispin           = 2
  Conv            = 1d0
  P(:,:)          = PHF(:,:)
  eGT(:)          = eHF(:)
  eOld(:)         = eHF(:)
  c(:,:)          = cHF(:,:)
  F_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  rcond           = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

    ! Increment

    nSCF = nSCF + 1

    ! Buid Coulomb matrix

    call Coulomb_matrix_AO_basis(nBas,P,ERI_AO,J)

    ! Compute exchange part of the self-energy 

    call exchange_matrix_AO_basis(nBas,P,ERI_AO,K)

    ! AO to MO transformation of two-electron integrals

    call AOtoMO_integral_transform(1,1,1,1,nBas,c,ERI_AO,ERI_MO)

    ! Compute linear response

    call phLR(ispin,.false.,TDA_T,eta,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI_MO,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
    if(print_T) call print_excitation('RPA@qsGTeh  ',ispin,nS,OmRPA)

    ! Compute correlation part of the self-energy 

    call GTeh_excitation_density(nBas,nC,nO,nR,nS,ERI_MO,XpY_RPA,XmY_RPA,rhoL_RPA,rhoR_RPA)

     if(regularize) then

!      call regularized_self_energy_correlation(eta,nBas,nC,nO,nV,nR,nS,eGT,OmRPA,rhoL_RPA,rhoR_RPA,EcGM,SigC)
!      call regularized_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS,eGT,OmRPA,rhoL_RPA,rhoR_RPA,Z)

     else

       call GTeh_self_energy(eta,nBas,nC,nO,nV,nR,nS,eGT,OmRPA,rhoL_RPA,rhoR_RPA,EcGM,SigC)
       call GTeh_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS,eGT,OmRPA,rhoL_RPA,rhoR_RPA,Z)

     endif

    ! Make correlation self-energy Hermitian and transform it back to AO basis
   
    SigCp = 0.5d0*(SigC + transpose(SigC))
    SigCm = 0.5d0*(SigC - transpose(SigC))

    call MOtoAO_transform(nBas,S,c,SigCp)
 
    ! Solve the quasi-particle equation

    F(:,:) = Hc(:,:) + J(:,:) + 0.5d0*K(:,:) + SigCp(:,:)

    ! Compute commutator and convergence criteria

    error = matmul(F,matmul(P,S)) - matmul(matmul(S,P),F)

    ! DIIS extrapolation 

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nBasSq,nBasSq,n_diis,error_diis,F_diis,error,F)

    end if

    ! Diagonalize Hamiltonian in AO basis

    Fp = matmul(transpose(X),matmul(F,X))
    cp(:,:) = Fp(:,:)
    call diagonalize_matrix(nBas,cp,eGT)
    c = matmul(X,cp)
    SigCp = matmul(transpose(c),matmul(SigCp,c))

    ! Compute new density matrix in the AO basis

    P(:,:) = 2d0*matmul(c(:,1:nO),transpose(c(:,1:nO)))

    ! Save quasiparticles energy for next cycle

    Conv = maxval(abs(error))
    eOld(:) = eGT(:)

    !------------------------------------------------------------------------
    !   Compute total energy
    !------------------------------------------------------------------------

    ! Kinetic energy

    ET = trace_matrix(nBas,matmul(P,T))

    ! Potential energy

    EV = trace_matrix(nBas,matmul(P,V))

    ! Coulomb energy

    EJ = 0.5d0*trace_matrix(nBas,matmul(P,J))

    ! Exchange energy

    Ex = 0.25d0*trace_matrix(nBas,matmul(P,K))

    ! Total energy

    EqsGT = ET + EV + EJ + Ex 

    ! Print results

    call dipole_moment(nBas,P,nNuc,ZNuc,rNuc,dipole_int_AO,dipole)
    call print_qsGTeh(nBas,nO,nSCF,Conv,thresh,eHF,eGT,c,SigCp,Z,ENuc,ET,EV,EJ,Ex,EcGM,EcRPA,EqsGT,dipole)

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

! Deallocate memory

  deallocate(c,cp,P,F,Fp,J,K,SigC,SigCp,SigCm,Z,OmRPA,XpY_RPA,XmY_RPA,rhoL_RPA,rhoR_RPA,error,error_diis,F_diis)

! Perform BSE calculation

! if(BSE) then

!   call Bethe_Salpeter(BSE2,TDA_T,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int_MO, & 
!                       eGW,eGW,EcBSE)

!   if(exchange_kernel) then

!     EcBSE(1) = 0.5d0*EcBSE(1)
!     EcBSE(2) = 1.5d0*EcBSE(2)

!   end if

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW correlation energy (singlet) =',EcBSE(1)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW correlation energy (triplet) =',EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW correlation energy           =',EcBSE(1) + EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@qsGW total energy                 =',ENuc + EqsGW + EcBSE(1) + EcBSE(2)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

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

!     call ACFDT(exchange_kernel,doXBS,.true.,TDA_T,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,eGW,eGW,EcAC)

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy (singlet) =',EcAC(1)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy (triplet) =',EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW correlation energy           =',EcAC(1) + EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@qsGW total energy                 =',ENuc + EqsGW + EcAC(1) + EcAC(2)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

! end if

end subroutine
