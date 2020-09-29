subroutine lCCD(maxSCF,thresh,max_diis,nBasin,nCin,nOin,nVin,nRin,ERI,ENuc,ERHF,eHF)

! Ladder CCD module

  implicit none

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh

  integer,intent(in)            :: nBasin
  integer,intent(in)            :: nCin
  integer,intent(in)            :: nOin
  integer,intent(in)            :: nVin
  integer,intent(in)            :: nRin
  double precision,intent(in)   :: ENuc,ERHF
  double precision,intent(in)   :: eHF(nBasin)
  double precision,intent(in)   :: ERI(nBasin,nBasin,nBasin,nBasin)

! Local variables

  integer                       :: nBas
  integer                       :: nC
  integer                       :: nO
  integer                       :: nV
  integer                       :: nR
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: EcMP2
  double precision              :: ECCD,EcCCD
  double precision,allocatable  :: seHF(:)
  double precision,allocatable  :: sERI(:,:,:,:)
  double precision,allocatable  :: dbERI(:,:,:,:)

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta_OOVV(:,:,:,:)

  double precision,allocatable  :: OOOO(:,:,:,:)
  double precision,allocatable  :: OOVV(:,:,:,:)
  double precision,allocatable  :: OVOV(:,:,:,:)
  double precision,allocatable  :: VVVV(:,:,:,:)

  double precision,allocatable  :: X1(:,:,:,:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: X3(:,:)
  double precision,allocatable  :: X4(:,:,:,:)

  double precision,allocatable  :: u(:,:,:,:)
  double precision,allocatable  :: v(:,:,:,:)

  double precision,allocatable  :: r2(:,:,:,:)
  double precision,allocatable  :: t2(:,:,:,:)

  integer                       :: n_diis
  double precision              :: rcond
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: t_diis(:,:)

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|   ladder CCD calculation           |'
  write(*,*)'**************************************'
  write(*,*)

! Spatial to spin orbitals

  nBas = 2*nBasin
  nC   = 2*nCin
  nO   = 2*nOin
  nV   = 2*nVin
  nR   = 2*nRin

  allocate(seHF(nBas),sERI(nBas,nBas,nBas,nBas))

  call spatial_to_spin_MO_energy(nBasin,eHF,nBas,seHF)
  call spatial_to_spin_ERI(nBasin,ERI,nBas,sERI)

! Antysymmetrize ERIs

  allocate(dbERI(nBas,nBas,nBas,nBas))

  call antisymmetrize_ERI(2,nBas,sERI,dbERI)

  deallocate(sERI)

! Form energy denominator

  allocate(eO(nO),eV(nV))
  allocate(delta_OOVV(nO,nO,nV,nV))

  eO(:) = seHF(1:nO)
  eV(:) = seHF(nO+1:nBas)

  call form_delta_OOVV(nC,nO,nV,nR,eO,eV,delta_OOVV)

  deallocate(seHF)

! Create integral batches

  allocate(OOOO(nO,nO,nO,nO),OOVV(nO,nO,nV,nV),OVOV(nO,nV,nO,nV),VVVV(nV,nV,nV,nV))

  OOOO(:,:,:,:) = dbERI(   1:nO  ,   1:nO  ,   1:nO  ,   1:nO  )
  OOVV(:,:,:,:) = dbERI(   1:nO  ,   1:nO  ,nO+1:nBas,nO+1:nBas)
  OVOV(:,:,:,:) = dbERI(   1:nO  ,nO+1:nBas,   1:nO  ,nO+1:nBas)
  VVVV(:,:,:,:) = dbERI(nO+1:nBas,nO+1:nBas,nO+1:nBas,nO+1:nBas)

  deallocate(dbERI)
 
! MP2 guess amplitudes

  allocate(t2(nO,nO,nV,nV))

  t2(:,:,:,:) = -OOVV(:,:,:,:)/delta_OOVV(:,:,:,:)

  call CCD_correlation_energy(nC,nO,nV,nR,OOVV,t2,EcMP2)

! Memory allocation for DIIS

  allocate(error_diis(nO*nO*nV*nV,max_diis),t_diis(nO*nO*nV*nV,max_diis))

! Initialization

  allocate(r2(nO,nO,nV,nV),u(nO,nO,nV,nV),v(nO,nO,nV,nV))
  allocate(X1(nO,nO,nO,nO),X2(nV,nV),X3(nO,nO),X4(nO,nO,nV,nV))

  Conv = 1d0
  nSCF = 0

  n_diis          = 0
  t_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| ladder CCD calculation                           |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(CCD)','|','Ec(CCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Compute residual

    call form_ladder_r(nC,nO,nV,nR,OOOO,OOVV,VVVV,t2,r2)

    r2(:,:,:,:) = OOVV(:,:,:,:) + delta_OOVV(:,:,:,:)*t2(:,:,:,:) + r2(:,:,:,:) 

!   Check convergence 

    Conv = maxval(abs(r2(nC+1:nO,nC+1:nO,1:nV-nR,1:nV-nR)))
  
!   Update amplitudes

    t2(:,:,:,:) = t2(:,:,:,:) - r2(:,:,:,:)/delta_OOVV(:,:,:,:)

!   Compute correlation energy

    call CCD_correlation_energy(nC,nO,nV,nR,OOVV,t2,EcCCD)

!   Dump results

    ECCD = ERHF + EcCCD

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nO*nO*nV*nV,nO*nO*nV*nV,n_diis,error_diis,t_diis,-r2/delta_OOVV,t2)

    !  Reset DIIS if required

    if(abs(rcond) < 1d-15) n_diis = 0

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',ECCD+ENuc,'|',EcCCD,'|',Conv,'|'

  enddo
  write(*,*)'----------------------------------------------------'
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

    stop

  endif

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'              ladder CCD energy                     '
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A30,1X,F15.10)')' E(lCCD)  = ',ECCD
  write(*,'(1X,A30,1X,F15.10)')' Ec(lCCD) = ',EcCCD
  write(*,*)'----------------------------------------------------'
  write(*,*)

end subroutine lCCD
