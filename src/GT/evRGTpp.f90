subroutine evRGTpp(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,BSE,TDA_T,TDA,dBSE,dTDA,singlet,triplet, & 
                   linearize,eta,regularize,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform eigenvalue self-consistent calculation with a T-matrix self-energy (evGT)

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
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)


! Local variables

  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: rcond
  double precision              :: Conv
  integer                       :: ispin
  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  double precision              :: EcGM
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)
  double precision,allocatable  :: Om1s(:),Om1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Om2s(:),Om2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)
  double precision,allocatable  :: Sig(:)
  double precision,allocatable  :: Z(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Restricted evGTpp Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! TDA for T

  if(TDA_T) then
    write(*,*) 'Tamm-Dancoff approximation for pp T-matrix!'
    write(*,*)
  end if

! Dimensions of the pp-RPA linear reponse matrices

  nOOs = nO*(nO+1)/2
  nVVs = nV*(nV+1)/2

  nOOt = nO*(nO-1)/2
  nVVt = nV*(nV-1)/2

! Memory allocation

  allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),        & 
           Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs),        & 
           rho1s(nOrb,nOrb,nVVs),rho2s(nOrb,nOrb,nOOs),     & 
           Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),        & 
           Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt),        & 
           rho1t(nOrb,nOrb,nVVt),rho2t(nOrb,nOrb,nOOt),     &
           eGT(nOrb),eOld(nOrb),Z(nOrb),Sig(nOrb),          &
           error_diis(nOrb,max_diis),e_diis(nOrb,max_diis))

! Initialization

  nSCF            = 0
  n_diis          = 0
  Conv            = 1d0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGT(:)          = eHF(:)
  eOld(:)         = eGT(:)
  Z(:)            = 1d0
  rcond           = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

  !----------------------------------------------
  ! alpha-beta block
  !----------------------------------------------

    ispin  = 1
 
    ! Compute linear response
 
    allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

    if(.not.TDA_T) call ppRLR_B(ispin,nOrb,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                   call ppRLR_C(ispin,nOrb,nC,nO,nV,nR,nVVs,1d0,eGT,ERI,Cpp)
                   call ppRLR_D(ispin,nOrb,nC,nO,nV,nR,nOOs,1d0,eGT,ERI,Dpp)

    call ppRLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(ispin))

    deallocate(Bpp,Cpp,Dpp)
 
  !----------------------------------------------
  ! alpha-alpha block
  !----------------------------------------------

    ispin  = 2

  ! Compute linear response

    allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))

    if(.not.TDA_T) call ppRLR_B(ispin,nOrb,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                   call ppRLR_C(ispin,nOrb,nC,nO,nV,nR,nVVt,1d0,eGT,ERI,Cpp)
                   call ppRLR_D(ispin,nOrb,nC,nO,nV,nR,nOOt,1d0,eGT,ERI,Dpp)

    call ppRLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(ispin))

    deallocate(Bpp,Cpp,Dpp)

    EcRPA(1) = 1d0*EcRPA(1)
    EcRPA(2) = 3d0*EcRPA(2)

  !----------------------------------------------
  ! Compute T-matrix version of the self-energy 
  !----------------------------------------------
 
    ispin = 3
    call RGTpp_excitation_density(ispin,nOrb,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)
 
    ispin = 4
    call RGTpp_excitation_density(ispin,nOrb,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)
 
  !----------------------------------------------
  ! Compute T-matrix version of the self-energy 
  !----------------------------------------------

    if(regularize) then 
      call GTpp_regularization(eta,nOrb,nC,nO,nV,nR,nOOs,nVVs,eGT,Om1s,rho1s,Om2s,rho2s)
      call GTpp_regularization(eta,nOrb,nC,nO,nV,nR,nOOt,nVVt,eGT,Om1t,rho1t,Om2t,rho2t)
    end if 

    call RGTpp_self_energy_diag(eta,nOrb,nC,nO,nV,nR,nooS,nVVt,nOOt,nVVt,eGT,Om1s,rho1s,Om2s,rho2s, & 
                               Om1t,rho1t,Om2t,rho2t,EcGM,Sig,Z)

  !----------------------------------------------
  ! Solve the quasi-particle equation
  !----------------------------------------------

    if(linearize) then
 
      write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
      write(*,*)
 
      eGT(:) = eHF(:) + Sig(:)
 
    else
 
      write(*,*) ' *** Quasiparticle energies obtained by root search *** '
      write(*,*)
 
     call RGTpp_QP_graph(eta,nOrb,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF,Om1s,rho1s,Om2s,rho2s, & 
                        Om1t,rho1t,Om2t,rho2t,eOld,eOld,eGT,Z)
 
    end if

    ! Convergence criteria

    Conv = maxval(abs(eGT(:) - eOld(:)))

  !----------------------------------------------
  ! Dump results
  !----------------------------------------------

    call print_evRGTpp(nOrb,nO,nSCF,Conv,eHF,ENuc,ERHF,Sig,Z,eGT,EcGM,EcRPA)

    ! DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call DIIS_extrapolation(rcond,nOrb,nOrb,n_diis,error_diis,e_diis,eGT(:)-eOld(:),eGT(:))

    end if

    ! Save quasiparticles energy for next cycle

    eOld(:) = eGT(:)

    ! Increment

    nSCF = nSCF + 1

  end do
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Perform BSE calculation

  if(BSE) then

    call RGTpp_phBSE(exchange_kernel,TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, &
                    Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,Om2t,X2t,Y2t,rho1t,rho2t,                      &
                    ERI,dipole_int,eGT,eGT,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@evGTpp correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@evGTpp correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@evGTpp correlation energy           =',sum(EcBSE)
    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@evGTpp total energy                 =',ENuc + ERHF + sum(EcBSE)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      call RGTpp_phACFDT(exchange_kernel,doXBS,.false.,TDA_T,TDA,BSE,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS, &
                        nOOs,nVVs,nOOt,nVVt,Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,               &
                        Om2t,X2t,Y2t,rho1t,rho2t,ERI,eGT,eGT,EcBSE)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10)') 'AC@phBSE@evGTpp correlation energy (singlet) =',EcBSE(1)
      write(*,'(2X,A50,F20.10)') 'AC@phBSE@evGTpp correlation energy (triplet) =',EcBSE(2)
      write(*,'(2X,A50,F20.10)') 'AC@phBSE@evGTpp correlation energy           =',sum(EcBSE)
      write(*,'(2X,A50,F20.10)') 'AC@phBSE@evGTpp total energy                 =',ENuc + ERHF + sum(EcBSE)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

! Testing zone

  if(dotest) then

    call dump_test_value('R','evGTpp correlation energy',sum(EcRPA))
    call dump_test_value('R','evGTpp HOMO energy',eGT(nO))
    call dump_test_value('R','evGTpp LUMO energy',eGT(nO+1))

  end if

end subroutine 
