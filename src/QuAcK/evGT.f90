subroutine evGT(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS, & 
                BSE,TDA_W,TDA,dBSE,dTDA,evDyn,singlet,triplet,      & 
                eta,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF,eG0T0)

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
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: eG0T0(nBas)


! Local variables

  logical                       :: linear_mixing
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: rcond
  double precision              :: Conv
  integer                       :: ispin
  integer                       :: iblock
  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  double precision              :: dERI
  double precision              :: xERI
  double precision              :: alpha
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: Omega1s(:),Omega1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Omega2s(:),Omega2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)
  double precision,allocatable  :: SigT(:)
  double precision,allocatable  :: Z(:)

  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)
  double precision,allocatable  :: XmY(:,:,:)
  double precision,allocatable  :: rho(:,:,:,:)

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

  allocate(Omega1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs), & 
           Omega2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs), & 
           rho1s(nBas,nO,nVVs),rho2s(nBas,nV,nOOs), & 
           Omega1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt), & 
           Omega2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt), & 
           rho1t(nBas,nO,nVVt),rho2t(nBas,nV,nOOt),     &
           eGT(nBas),eOld(nBas),Z(nBas),SigT(nBas),     &
           error_diis(nBas,max_diis),e_diis(nBas,max_diis))

! Initialization

  nSCF            = 0
  n_diis          = 0
  Conv            = 1d0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGT(:)          = eG0T0(:)
  eOld(:)         = eGT(:)
  Z(:)            = 1d0

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
 
    call linear_response_pp(iblock,.true.,.false.,nBas,nC,nO,nV,nR,nOOs,nVVs,eHF(:),ERI(:,:,:,:),  & 
                            Omega1s(:),X1s(:,:),Y1s(:,:),Omega2s(:),X2s(:,:),Y2s(:,:),EcRPA(ispin))
 
!   EcRPA(ispin) = 1d0*EcRPA(ispin)
 
!   call print_excitation('pp-RPA (N+2)',iblock,nVVs,Omega1s(:))
!   call print_excitation('pp-RPA (N-2)',iblock,nOOs,Omega2s(:))

  !----------------------------------------------
  ! alpha-alpha block
  !----------------------------------------------

    ispin  = 2
    iblock = 4

  ! Compute linear response

    call linear_response_pp(iblock,.true.,.false.,nBas,nC,nO,nV,nR,nOOt,nVVt,eHF(:),ERI(:,:,:,:),  & 
                          Omega1t(:),X1t(:,:),Y1t(:,:),Omega2t(:),X2t(:,:),Y2t(:,:),EcRPA(ispin))

!   EcRPA(ispin) = 2d0*EcRPA(ispin)
!   EcRPA(ispin) = 3d0*EcRPA(ispin)

!   call print_excitation('pp-RPA (N+2)',iblock,nVVt,Omega1t(:))
!   call print_excitation('pp-RPA (N-2)',iblock,nOOt,Omega2t(:))

  !----------------------------------------------
  ! Compute T-matrix version of the self-energy 
  !----------------------------------------------

    SigT(:) = 0d0
    Z(:)    = 0d0
 
    iblock =  3
    dERI   = +1d0
    xERI   = +0d0
    alpha  = +1d0
 
    call excitation_density_Tmatrix(iblock,dERI,xERI,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI(:,:,:,:), &
                                    X1s(:,:),Y1s(:,:),rho1s(:,:,:),X2s(:,:),Y2s(:,:),rho2s(:,:,:))
 
    call self_energy_Tmatrix_diag(alpha,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,eGT(:), & 
                                  Omega1s(:),rho1s(:,:,:),Omega2s(:),rho2s(:,:,:),SigT(:))
 
    call renormalization_factor_Tmatrix(alpha,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,eGT(:), & 
                                        Omega1s(:),rho1s(:,:,:),Omega2s(:),rho2s(:,:,:),Z(:))
 
    iblock =  4
    dERI   = +1d0
    xERI   = -1d0
    alpha  = +1d0
 
    call excitation_density_Tmatrix(iblock,dERI,xERI,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI(:,:,:,:), &
                                    X1t(:,:),Y1t(:,:),rho1t(:,:,:),X2t(:,:),Y2t(:,:),rho2t(:,:,:))
 
    call self_energy_Tmatrix_diag(alpha,eta,nBas,nC,nO,nV,nR,nOOt,nVVt,eGT(:), & 
                                  Omega1t(:),rho1t(:,:,:),Omega2t(:),rho2t(:,:,:),SigT(:))
 
    call renormalization_factor_Tmatrix(alpha,eta,nBas,nC,nO,nV,nR,nOOt,nVVt,eGT(:), & 
                                        Omega1t(:),rho1t(:,:,:),Omega2t(:),rho2t(:,:,:),Z(:))
 
    Z(:) = 1d0/(1d0 - Z(:))

    ! Solve the quasi-particle equation

  !----------------------------------------------
  ! Solve the quasi-particle equation
  !----------------------------------------------

    eGT(:) = eHF(:) + SigT(:)

    ! Convergence criteria

    Conv = maxval(abs(eGT(:) - eOld(:)))

  !----------------------------------------------
  ! Dump results
  !----------------------------------------------

    call print_evGT(nBas,nO,nSCF,Conv,eHF(:),SigT(:),Z(:),eGT(:))

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

! Compute the ppRPA correlation energy

  ispin  = 1
  iblock = 3
  call linear_response_pp(iblock,.false.,.false.,nBas,nC,nO,nV,nR,nOOs,nVVs,eGT(:),ERI(:,:,:,:),  &
                          Omega1s(:),X1s(:,:),Y1s(:,:),Omega2s(:),X2s(:,:),Y2s(:,:),EcRPA(ispin))
  ispin  = 2
  iblock = 4
  call linear_response_pp(iblock,.false.,.false.,nBas,nC,nO,nV,nR,nOOt,nVVt,eGT(:),ERI(:,:,:,:),  &
                          Omega1t(:),X1t(:,:),Y1t(:,:),Omega2t(:),X2t(:,:),Y2t(:,:),EcRPA(ispin))
  EcRPA(1) = EcRPA(1) - EcRPA(2)
  EcRPA(2) = 3d0*EcRPA(2)

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA@evGT correlation energy (singlet) =',EcRPA(1)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA@evGT correlation energy (triplet) =',EcRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA@evGT correlation energy           =',EcRPA(1) + EcRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@ppRPA@evGT total energy                 =',ENuc + ERHF + EcRPA(1) + EcRPA(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)


! Perform BSE calculation

  if(BSE) then

    call Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta, &
                        nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGT,eGT,EcRPA,EcBSE)

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

      call ACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,singlet,triplet,eta, &
                 nBas,nC,nO,nV,nR,nS,ERI,eGT,eGT,EcAC)

      if(exchange_kernel) then

        EcAC(1) = 0.5d0*EcAC(1)
        EcAC(2) = 1.5d0*EcAC(1)

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

end subroutine evGT
