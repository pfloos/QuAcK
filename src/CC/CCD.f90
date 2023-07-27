subroutine CCD(maxSCF,thresh,max_diis,nBasin,nCin,nOin,nVin,nRin,ERI,ENuc,ERHF,eHF)

! CCD module

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
  double precision              :: EcMP2,EcMP3,EcMP4
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

  double precision,allocatable  :: OVVO(:,:,:,:)

  double precision,allocatable  :: X1(:,:,:,:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: X3(:,:)
  double precision,allocatable  :: X4(:,:,:,:)

  double precision,allocatable  :: u(:,:,:,:)
  double precision,allocatable  :: v(:,:,:,:)

  double precision,allocatable  :: r(:,:,:,:)
  double precision,allocatable  :: t(:,:,:,:)

  integer                       :: n_diis,i,j,a,b
  double precision              :: rcond
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: t_diis(:,:)

  logical                       :: do_EE_EOM_CC_1h1p = .false.
  logical                       :: do_EA_EOM_CC_1p   = .false.
  logical                       :: do_IP_EOM_CC_1h   = .false.
  logical                       :: do_DEA_EOM_CC_2p  = .false.
  logical                       :: do_DIP_EOM_CC_2h  = .false.

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|          CCD calculation           |'
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

  allocate(eO(nO-nC),eV(nV-nR))
  allocate(delta_OOVV(nO-nC,nO-nC,nV-nR,nV-nR))

  eO(:) = seHF(nC+1:nO)
  eV(:) = seHF(nO+1:nBas-nR)

  call form_delta_OOVV(nC,nO,nV,nR,eO,eV,delta_OOVV)

  deallocate(seHF)

! Create integral batches

  allocate(OOOO(nO-nC,nO-nC,nO-nC,nO-nC),OOVV(nO-nC,nO-nC,nV-nR,nV-nR), & 
           OVOV(nO-nC,nV-nR,nO-nC,nV-nR),VVVV(nV-nR,nV-nR,nV-nR,nV-nR))

  OOOO(:,:,:,:) = dbERI(nC+1:nO     ,nC+1:nO     ,nC+1:nO     ,nC+1:nO     )
  OOVV(:,:,:,:) = dbERI(nC+1:nO     ,nC+1:nO     ,nO+1:nBas-nR,nO+1:nBas-nR)
  OVOV(:,:,:,:) = dbERI(nC+1:nO     ,nO+1:nBas-nR,nC+1:nO     ,nO+1:nBas-nR)
  VVVV(:,:,:,:) = dbERI(nO+1:nBas-nR,nO+1:nBas-nR,nO+1:nBas-nR,nO+1:nBas-nR)

  allocate(OVVO(nO-nC,nV-nR,nV-nR,nO-nC))
  OVVO(:,:,:,:) = dbERI(nC+1:nO,nO+1:nBas-nR,nO+1:nBas-nR,nC+1:nO)

  deallocate(dbERI)
 
! MP2 guess amplitudes

  allocate(t(nO-nC,nO-nC,nV-nR,nV-nR))

  t(:,:,:,:) = -OOVV(:,:,:,:)/delta_OOVV(:,:,:,:)

  call CCD_correlation_energy(nC,nO,nV,nR,OOVV,t,EcMP2)
  EcMP4 = 0d0

! Memory allocation for DIIS

  allocate(error_diis((nO-nR)**2*(nV-nR)**2,max_diis),t_diis((nO-nR)**2*(nV-nR)**2,max_diis))

! Initialization

  allocate(r(nO-nC,nO-nC,nV-nR,nV-nR),u(nO-nC,nO-nC,nV-nR,nV-nR),v(nO-nC,nO-nC,nV-nR,nV-nR))
  allocate(X1(nO-nC,nO-nC,nO-nC,nO-nC),X2(nV-nR,nV-nR),X3(nO-nC,nO-nC),X4(nO-nC,nO-nC,nV-nR,nV-nR))

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
  write(*,*)'| CCD calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(CCD)','|','Ec(CCD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Form linear array

    call form_u(nC,nO,nV,nR,OOOO,VVVV,OVOV,t,u)

!   Form interemediate arrays

    call form_X(nC,nO,nV,nR,OOVV,t,X1,X2,X3,X4)

!   Form quadratic array

    call form_v(nC,nO,nV,nR,X1,X2,X3,X4,t,v)

!   Compute residual

    r(:,:,:,:) = OOVV(:,:,:,:) + delta_OOVV(:,:,:,:)*t(:,:,:,:) + u(:,:,:,:) + v(:,:,:,:)

!   Check convergence 

    Conv = maxval(abs(r(:,:,:,:)))
  
!   Update amplitudes

    t(:,:,:,:) = t(:,:,:,:) - r(:,:,:,:)/delta_OOVV(:,:,:,:)

!   Compute correlation energy

    call CCD_correlation_energy(nC,nO,nV,nR,OOVV,t,EcCCD)

    if(nSCF == 1) call MP3_correlation_energy(nC,nO,nV,nR,OOVV,t,v,delta_OOVV,EcMP3)

!   Dump results

    ECCD = ERHF + EcCCD

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,(nO-nC)**2*(nV-nR)**2,(nO-nC)**2*(nV-nR)**2,n_diis,error_diis,t_diis,-r/delta_OOVV,t)

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
  write(*,*)'              CCD energy                            '
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A30,1X,F15.10)')' E(CCD) = ',ECCD
  write(*,'(1X,A30,1X,F15.10)')' Ec(CCD) = ',EcCCD
  write(*,*)'----------------------------------------------------'
  write(*,*)

! Moller-Plesset energies

  write(*,*)
  write(*,'(1X,A15,1X,F10.6)') 'Ec(MP2)     = ',EcMP2
  write(*,'(1X,A15,1X,F10.6)') 'Ec(MP3)     = ',EcMP3
  write(*,'(1X,A15,1X,F10.6)') 'Ec(MP4-SDQ) = ',EcMP4
  write(*,*)

!------------------------------------------------------------------------
! EOM section
!------------------------------------------------------------------------

! EE-EOM-CCD (1h1p)

  if(do_EE_EOM_CC_1h1p)  call EE_EOM_CCD_1h1p(nC,nO,nV,nR,eO,eV,OOVV,OVVO,t)

! EA-EOM (1p)

! if(do_EA-EOM-CC_1p) call EA-EOM-CCD_1p()

! IP-EOM-CCD(1h)

! if(do_IP-EOM-CC_1h) call IP-EOM-CCD_1h()

! DEA-EOM (2p)

  if(do_DEA_EOM_CC_2p) call DEA_EOM_CCD_2p(nC,nO,nV,nR,eV,OOVV,VVVV,t)

! DIP-EOM-CCD(2h)

  if(do_DIP_EOM_CC_2h) call DIP_EOM_CCD_2h(nC,nO,nV,nR,eO,OOVV,OOOO,t)

end subroutine 
