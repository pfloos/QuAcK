subroutine drCCD(dotest,maxSCF,thresh,max_diis,nBasin,nCin,nOin,nVin,nRin,ERI,ENuc,ERHF,eHF)

! Direct ring CCD module

  implicit none

! Input variables

  logical,intent(in)            :: dotest

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
  double precision              :: ECC,EcCC
  double precision,allocatable  :: seHF(:)
  double precision,allocatable  :: sERI(:,:,:,:)

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta_OOVV(:,:,:,:)

  double precision,allocatable  :: OOVV(:,:,:,:)
  double precision,allocatable  :: OVVO(:,:,:,:)

  double precision,allocatable  :: r2(:,:,:,:)
  double precision,allocatable  :: t2(:,:,:,:)


  integer                       :: n_diis
  double precision              :: rcond
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: t_diis(:,:)

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|     direct ring CCD calculation    |'
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

! Form energy denominator

  allocate(eO(nO),eV(nV))
  allocate(delta_OOVV(nO,nO,nV,nV))

  eO(:) = seHF(1:nO)
  eV(:) = seHF(nO+1:nBas)

  call form_delta_OOVV(nC,nO,nV,nR,eO,eV,delta_OOVV)

  deallocate(seHF)

! Create integral batches

  allocate(OOVV(nO,nO,nV,nV),OVVO(nO,nV,nV,nO))

  OOVV(:,:,:,:) = sERI(   1:nO  ,   1:nO  ,nO+1:nBas,nO+1:nBas)
  OVVO(:,:,:,:) = sERI(   1:nO  ,nO+1:nBas,nO+1:nBas,   1:nO  )
 
  deallocate(sERI)

! MP2 guess amplitudes

  allocate(t2(nO,nO,nV,nV))

  t2(:,:,:,:) = -OOVV(:,:,:,:)/delta_OOVV(:,:,:,:)

  call CCD_correlation_energy(nC,nO,nV,nR,OOVV,t2,EcMP2)

! Memory allocation for DIIS

  allocate(error_diis(nO*nO*nV*nV,max_diis),t_diis(nO*nO*nV*nV,max_diis))

! Initialization

  allocate(r2(nO,nO,nV,nV))

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
  write(*,*)'| direct ring CCD calculation                      |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(CCD)','|','Ec(CCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Compute residual

    call form_ring_r(nC,nO,nV,nR,OVVO,OOVV,t2,r2)

    r2(:,:,:,:) = OOVV(:,:,:,:) + delta_OOVV(:,:,:,:)*t2(:,:,:,:) + r2(:,:,:,:) 

!   Check convergence 

    Conv = maxval(abs(r2(nC+1:nO,nC+1:nO,1:nV-nR,1:nV-nR)))

!   Update amplitudes

    t2(:,:,:,:) = t2(:,:,:,:) - r2(:,:,:,:)/delta_OOVV(:,:,:,:)

!   Compute correlation energy

    call CCD_correlation_energy(nC,nO,nV,nR,OOVV,t2,EcCC)
    EcCC = 2d0*EcCC

!   Dump results

    ECC = ERHF + EcCC

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nO*nO*nV*nV,nO*nO*nV*nV,n_diis,error_diis,t_diis,-r2/delta_OOVV,t2)

    !  Reset DIIS if required

    if(abs(rcond) < 1d-15) n_diis = 0

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',ECC+ENuc,'|',EcCC,'|',Conv,'|'

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
  write(*,*)'              direct ring CCD energy                '
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A30,1X,F15.10)')' E(drCCD) = ',ECC  
  write(*,'(1X,A30,1X,F15.10)')' Ec(drCCD) = ',EcCC 
  write(*,*)'----------------------------------------------------'
  write(*,*)

  if(dotest) then

    call dump_test_value('R','drCCD correlation energy',EcCC)

  end if

end subroutine 
