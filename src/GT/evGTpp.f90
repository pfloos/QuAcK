subroutine evGTpp(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS, & 
                  BSE,TDA_T,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,regularize,nBas, & 
                  nC,nO,nV,nR,nS,ENuc,ERHF,ERI_AO,ERI_MO,dipole_int,PHF,cHF,eHF,Vxc,eG0T0)

! Perform eigenvalue self-consistent calculation with a T-matrix self-energy (evGT)

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
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: Vxc(nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_MO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: eG0T0(nBas)


! Local variables

  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: rcond
  double precision              :: Conv
  integer                       :: ispin
  integer                       :: iblock
  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  double precision              :: EcGM
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: Om1s(:),Om1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Om2s(:),Om2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)
  double precision,allocatable  :: SigX(:)
  double precision,allocatable  :: SigT(:)
  double precision,allocatable  :: Z(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|      Self-consistent evGT calculation        |'
  write(*,*)'************************************************'
  write(*,*)

! Dimensions of the pp-RPA linear reponse matrices

  nOOs = nO*nO
  nVVs = nV*nV

  nOOt = nO*(nO - 1)/2
  nVVt = nV*(nV - 1)/2

! Memory allocation

  allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),        & 
           Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs),        & 
           rho1s(nBas,nBas,nVVs),rho2s(nBas,nBas,nOOs),        & 
           Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),        & 
           Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt),        & 
           rho1t(nBas,nBas,nVVt),rho2t(nBas,nBas,nOOt),        &
           eGT(nBas),eOld(nBas),Z(nBas),SigX(nBas),SigT(nBas), &
           error_diis(nBas,max_diis),e_diis(nBas,max_diis))

! Compute the exchange part of the self-energy

  call self_energy_exchange_diag(nBas,cHF,PHF,ERI_AO,SigX)

! Initialization

  nSCF            = 0
  n_diis          = 0
  Conv            = 1d0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGT(:)          = eG0T0(:)
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
    iblock = 3
 
    ! Compute linear response
 
    call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,eGT,ERI_MO,  & 
                            Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(ispin))
 
  !----------------------------------------------
  ! alpha-alpha block
  !----------------------------------------------

    ispin  = 2
    iblock = 4

  ! Compute linear response

    call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,eGT,ERI_MO,  & 
                          Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(ispin))

    EcRPA(1) = EcRPA(1) - EcRPA(2)
    EcRPA(2) = 3d0*EcRPA(2)

  !----------------------------------------------
  ! Compute T-matrix version of the self-energy 
  !----------------------------------------------

    EcGM    = 0d0
    SigT(:) = 0d0
    Z(:)    = 0d0
 
    iblock =  3
 
    call GTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI_MO, &
                                 X1s,Y1s,rho1s,X2s,Y2s,rho2s)
 
    call GTpp_self_energy_diag(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,eGT, & 
                               Om1s,rho1s,Om2s,rho2s,EcGM,SigT)
 
    call GTpp_renormalization_factor(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,eGT, & 
                                     Om1s,rho1s,Om2s,rho2s,Z)
 
    iblock =  4
 
    call GTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI_MO, &
                                 X1t,Y1t,rho1t,X2t,Y2t,rho2t)
 
    call GTpp_self_energy_diag(eta,nBas,nC,nO,nV,nR,nOOt,nVVt,eGT, & 
                               Om1t,rho1t,Om2t,rho2t,EcGM,SigT)
 
    call GTpp_renormalization_factor(eta,nBas,nC,nO,nV,nR,nOOt,nVVt,eGT, & 
                                     Om1t,rho1t,Om2t,rho2t,Z)
 
    Z(:) = 1d0/(1d0 - Z(:))

    ! Solve the quasi-particle equation

  !----------------------------------------------
  ! Solve the quasi-particle equation
  !----------------------------------------------

    eGT(:) = eHF(:) + SigX(:) + SigT(:) - Vxc(:)

    ! Convergence criteria

    Conv = maxval(abs(eGT(:) - eOld(:)))

  !----------------------------------------------
  ! Dump results
  !----------------------------------------------

    call print_evGTpp(nBas,nO,nSCF,Conv,eHF,ENuc,ERHF,SigT,Z,eGT,EcGM,EcRPA)

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nBas,nBas,n_diis,error_diis,e_diis,eGT(:)-eOld(:),eGT(:))

    ! Reset DIIS if required

    if(abs(rcond) < 1d-15) n_diis = 0

    ! Save quasiparticles energy for next cycle

    eOld(:) = eGT(:)

    ! Increment

    nSCF = nSCF + 1

  enddo
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Perform BSE calculation

  if(BSE) then

    call GTpp_phBSE(TDA_T,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,   &
                    Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,Om2t,X2t,Y2t,rho1t,rho2t, &
                    ERI_MO,dipole_int,eGT,eGT,EcBSE)

    if(exchange_kernel) then

      EcRPA(1) = 0.5d0*EcRPA(1)
      EcRPA(2) = 1.5d0*EcRPA(1)

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGT correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGT correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGT correlation energy           =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGT total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
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

      call ACFDT_Tmatrix(exchange_kernel,doXBS,.false.,TDA_T,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS, &
                         nOOs,nVVs,nOOt,nVVt,Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,     &
                         Om2t,X2t,Y2t,rho1t,rho2t,ERI_MO,eGT,eGT,EcAC)

      if(exchange_kernel) then

        EcAC(1) = 0.5d0*EcAC(1)
        EcAC(2) = 1.5d0*EcAC(2)

      end if

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evGT correlation energy (singlet) =',EcAC(1)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evGT correlation energy (triplet) =',EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evGT correlation energy           =',EcAC(1) + EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evGT total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

end subroutine 
