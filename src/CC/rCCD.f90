subroutine rCCD(dotest,maxSCF,thresh,max_diis,nBasin,nCin,nOin,nVin,nRin,ERI,ENuc,ERHF,eHF,eGW)

! Ring CCD module

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
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBasin)
  double precision,intent(in)   :: eGW(nBasin)
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
  double precision,allocatable  :: seGW(:)
  double precision,allocatable  :: sERI(:,:,:,:)
  double precision,allocatable  :: dbERI(:,:,:,:)

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: delta_OOVV(:,:,:,:)

  double precision,allocatable  :: OOVV(:,:,:,:)
  double precision,allocatable  :: OVVO(:,:,:,:)

  double precision,allocatable  :: r(:,:,:,:)
  double precision,allocatable  :: t(:,:,:,:)

  integer                       :: n_diis
  double precision              :: rcond
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: t_diis(:,:)

  logical                       :: do_EE_EOM_CC_1h1p = .true.

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|     ring CCD calculation           |'
  write(*,*)'**************************************'
  write(*,*)

! Spatial to spin orbitals

  nBas = 2*nBasin
  nC   = 2*nCin
  nO   = 2*nOin
  nV   = 2*nVin
  nR   = 2*nRin

  allocate(seHF(nBas),seGW(nBas),sERI(nBas,nBas,nBas,nBas))

  call spatial_to_spin_MO_energy(nBasin,eHF,nBas,seHF)
  call spatial_to_spin_MO_energy(nBasin,eGW,nBas,seGW)
  call spatial_to_spin_ERI(nBasin,ERI,nBas,sERI)

! Antysymmetrize ERIs

  allocate(dbERI(nBas,nBas,nBas,nBas))

  call antisymmetrize_ERI(2,nBas,sERI,dbERI)

  deallocate(sERI)

! Form energy denominator

  allocate(eO(nO),eV(nV))
  allocate(delta_OOVV(nO,nO,nV,nV))

  eO(:) = seGW(1:nO)
  eV(:) = seGW(nO+1:nBas)

  call form_delta_OOVV(nC,nO,nV,nR,eO,eV,delta_OOVV)

! deallocate(seHF,seGW)

! Create integral batches

  allocate(OOVV(nO,nO,nV,nV),OVVO(nO,nV,nV,nO))

  OOVV(:,:,:,:) = dbERI(   1:nO   ,   1:nO  ,nO+1:nBas,nO+1:nBas)
  OVVO(:,:,:,:) = dbERI(   1:nO   ,nO+1:nBas,nO+1:nBas,   1:nO  )

  deallocate(dbERI)
 
! MP2 guess amplitudes

  allocate(t(nO,nO,nV,nV))

  t(:,:,:,:) = -OOVV(:,:,:,:)/delta_OOVV(:,:,:,:)

  call CCD_correlation_energy(nC,nO,nV,nR,OOVV,t,EcMP2)

! Memory allocation for DIIS

  allocate(error_diis(nO*nO*nV*nV,max_diis),t_diis(nO*nO*nV*nV,max_diis))

! Initialization

  allocate(r(nO,nO,nV,nV))

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
  write(*,*)'| ring CCD calculation                             |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(CCD)','|','Ec(CCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Compute residual

    call form_ring_r(nC,nO,nV,nR,OVVO,OOVV,t,r)

    r(:,:,:,:) = OOVV(:,:,:,:) + delta_OOVV(:,:,:,:)*t(:,:,:,:) + r(:,:,:,:) 

!   Check convergence 

    Conv = maxval(abs(r(nC+1:nO,nC+1:nO,1:nV-nR,1:nV-nR)))
  
!   Update amplitudes

    t(:,:,:,:) = t(:,:,:,:) - r(:,:,:,:)/delta_OOVV(:,:,:,:)

!   Compute correlation energy

    call CCD_correlation_energy(nC,nO,nV,nR,OOVV,t,EcCC)

!   Dump results

    ECC = ERHF + EcCC

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nO*nO*nV*nV,nO*nO*nV*nV,n_diis,error_diis,t_diis,-r/delta_OOVV,t)

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
  write(*,*)'              ring CCD energy                       '
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A30,1X,F15.10)')' E(rCCD) = ',ECC  
  write(*,'(1X,A30,1X,F15.10)')' Ec(rCCD) = ',EcCC 
  write(*,*)'----------------------------------------------------'
  write(*,*)

! write(*,*)
! write(*,*)'----------------------------------------------------'
! write(*,*)'              ring CCD amplitudes                   '
! write(*,*)'----------------------------------------------------'
! call matout(nO*nO,nV*nV,t)
! write(*,*)

!------------------------------------------------------------------------
! EOM section
!------------------------------------------------------------------------

! EE-EOM-CCD (1h1p)

  if(do_EE_EOM_CC_1h1p) call EE_EOM_CCD_1h1p(nC,nO,nV,nR,eO,eV,OOVV,OVVO,t)


  if(dotest) then

    call dump_test_value('R','rCCD correlation energy',EcCC)

  end if

end subroutine 
