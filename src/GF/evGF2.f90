subroutine evGF2(BSE,TDA,dBSE,dTDA,evDyn,maxSCF,thresh,max_diis,singlet,triplet, &
                 linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform eigenvalue self-consistent second-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: max_diis
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nO
  integer,intent(in)            :: nC
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: Ec
  double precision              :: EcBSE(nspin)
  double precision              :: Conv
  double precision              :: rcond
  double precision,allocatable  :: eGF2(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|  Second-order Green function calculation     |'
  write(*,*)'************************************************'
  write(*,*)

! Memory allocation

  allocate(SigC(nBas),Z(nBas),eGF2(nBas),eOld(nBas),error_diis(nBas,max_diis),e_diis(nBas,max_diis))

! Initialization

  Conv            = 1d0
  nSCF            = 0
  n_diis          = 0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGF2(:)         = eHF(:)
  eOld(:)         = eHF(:)
  rcond           = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Frequency-dependent second-order contribution

    if(regularize) then 

      call regularized_self_energy_GF2_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,eGF2,ERI,SigC,Z)

    else

      call self_energy_GF2_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,eGF2,ERI,SigC,Z)

    end if

    if(linearize) then

      eGF2(:) = eHF(:) + Z(:)*SigC(:)

    else

      eGF2(:) = eHF(:) + SigC(:)

    end if

    Conv = maxval(abs(eGF2 - eOld))

    ! Print results

    call MP2(regularize,nBas,nC,nO,nV,nR,ERI,ENuc,EHF,eGF2,Ec)
    call print_evGF2(nBas,nO,nSCF,Conv,eHF,SigC,Z,eGF2,ENuc,ERHF,Ec)

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nBas,nBas,n_diis,error_diis,e_diis,eGF2-eOld,eGF2)

    if(abs(rcond) < 1d-15) n_diis = 0

    eOld(:) = eGF2(:)

    ! Increment

    nSCF = nSCF + 1

  end do
!------------------------------------------------------------------------
! End main SCF loop
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

! Perform BSE2 calculation

  if(BSE) then

    call BSE2(TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGF2,EcBSE)

  end if

end subroutine evGF2
