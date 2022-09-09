subroutine CCSD(BSE,maxSCF,thresh,max_diis,doCCSDT,nBasin,nCin,nOin,nVin,nRin,ERI,ENuc,ERHF,eHF)

! CCSD module

  implicit none

! Input variables

  logical,intent(in)            :: BSE
  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh

  logical,intent(in)            :: doCCSDT
  integer,intent(in)            :: nBasin
  integer,intent(in)            :: nCin
  integer,intent(in)            :: nOin
  integer,intent(in)            :: nVin
  integer,intent(in)            :: nRin
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBasin)
  double precision,intent(in)   :: ERI(nBasin,nBasin,nBasin,nBasin)

! Local variables

  double precision              :: start_CCSDT,end_CCSDT,t_CCSDT
  integer                       :: nBas
  integer                       :: nC
  integer                       :: nO
  integer                       :: nV
  integer                       :: nR
  integer                       :: nSCF
  double precision              :: Conv
  double precision              :: EcMP2
  double precision              :: ECCSD
  double precision              :: EcCCSD
  double precision              :: EcCCT

  double precision,allocatable  :: seHF(:)
  double precision,allocatable  :: sERI(:,:,:,:)
  double precision,allocatable  :: dbERI(:,:,:,:)
  double precision,allocatable  :: delta_OV(:,:)
  double precision,allocatable  :: delta_OOVV(:,:,:,:)

  double precision,allocatable  :: OOOO(:,:,:,:)
  double precision,allocatable  :: OOOV(:,:,:,:)
  double precision,allocatable  :: OVOO(:,:,:,:)
  double precision,allocatable  :: VOOO(:,:,:,:)
  double precision,allocatable  :: OOVV(:,:,:,:)
  double precision,allocatable  :: OVVO(:,:,:,:)
  double precision,allocatable  :: OVVV(:,:,:,:)
  double precision,allocatable  :: VOVV(:,:,:,:)
  double precision,allocatable  :: VVVO(:,:,:,:)
  double precision,allocatable  :: VVVV(:,:,:,:)

  double precision,allocatable  :: eO(:)
  double precision,allocatable  :: eV(:)
  double precision,allocatable  :: hvv(:,:)
  double precision,allocatable  :: hoo(:,:)
  double precision,allocatable  :: hvo(:,:)
  double precision,allocatable  :: gvv(:,:)
  double precision,allocatable  :: goo(:,:)
  double precision,allocatable  :: aoooo(:,:,:,:)
  double precision,allocatable  :: bvvvv(:,:,:,:)
  double precision,allocatable  :: hovvo(:,:,:,:)

  double precision,allocatable  :: r1(:,:)
  double precision,allocatable  :: r2(:,:,:,:)

  double precision,allocatable  :: t1(:,:)
  double precision,allocatable  :: t2(:,:,:,:)
  double precision,allocatable  :: tau(:,:,:,:)

  integer                       :: n_diis
  double precision              :: rcond1
  double precision              :: rcond2
  double precision,allocatable  :: err1_diis(:,:)
  double precision,allocatable  :: err2_diis(:,:)
  double precision,allocatable  :: t1_diis(:,:)
  double precision,allocatable  :: t2_diis(:,:)

! Hello world

  write(*,*)
  write(*,*)'**************************************'
  write(*,*)'|         CCSD calculation           |'
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

  if(BSE) then

    call static_screening(nBas,nC,nO,nV,nR,seHF,sERI,dbERI)

  else

    call antisymmetrize_ERI(2,nBas,sERI,dbERI)

  end if

  deallocate(sERI)

! Form energy denominator

  allocate(eO(nO),eV(nV))
  allocate(delta_OV(nO,nV),delta_OOVV(nO,nO,nV,nV))

  eO(:) = seHF(1:nO)
  eV(:) = seHF(nO+1:nBas)

  call form_delta_OV(nC,nO,nV,nR,eO,eV,delta_OV)
  call form_delta_OOVV(nC,nO,nV,nR,eO,eV,delta_OOVV)

  deallocate(seHF)

! Create integral batches

  allocate(OOOO(nO,nO,nO,nO),                                     & 
           OOOV(nO,nO,nO,nV),OVOO(nO,nV,nO,nO),VOOO(nV,nO,nO,nO), &
           OOVV(nO,nO,nV,nV),OVVO(nO,nV,nV,nO),                   & 
           OVVV(nO,nV,nV,nV),VOVV(nV,nO,nV,nV),VVVO(nV,nV,nV,nO), & 
           VVVV(nV,nV,nV,nV))

  OOOO(:,:,:,:) = dbERI(   1:nO  ,   1:nO  ,   1:nO  ,   1:nO  )
  OOOV(:,:,:,:) = dbERI(   1:nO  ,   1:nO  ,   1:nO  ,nO+1:nBas)
  OVOO(:,:,:,:) = dbERI(   1:nO  ,nO+1:nBas,   1:nO  ,   1:nO  )
  VOOO(:,:,:,:) = dbERI(nO+1:nBas,   1:nO  ,   1:nO  ,   1:nO  )
  OOVV(:,:,:,:) = dbERI(   1:nO  ,   1:nO  ,nO+1:nBas,nO+1:nBas)
  OVVO(:,:,:,:) = dbERI(   1:nO  ,nO+1:nBas,nO+1:nBas,   1:nO  )
  OVVV(:,:,:,:) = dbERI(   1:nO  ,nO+1:nBas,nO+1:nBas,nO+1:nBas)
  VOVV(:,:,:,:) = dbERI(nO+1:nBas,   1:nO  ,nO+1:nBas,nO+1:nBas)
  VVVO(:,:,:,:) = dbERI(nO+1:nBas,nO+1:nBas,nO+1:nBas,   1:nO  )
  VVVV(:,:,:,:) = dbERI(nO+1:nBas,nO+1:nBas,nO+1:nBas,nO+1:nBas)

  deallocate(dbERI)
 
! MP2 guess amplitudes

  allocate(t1(nO,nV),t2(nO,nO,nV,nV),tau(nO,nO,nV,nV))

  t1(:,:)     = 0d0
  t2(:,:,:,:) = -OOVV(:,:,:,:)/delta_OOVV(:,:,:,:)
  call form_tau(nC,nO,nV,nR,t1,t2,tau)

  EcMP2 = 0.5d0*dot_product(pack(OOVV,.true.),pack(tau,.true.))
  write(*,'(1X,A20,1X,F10.6)') 'Ec(MP2) = ',EcMP2

! Initialization

  allocate(hvv(nV,nV),hoo(nO,nO),hvo(nV,nO), &
           gvv(nV,nV),goo(nO,nO), & 
           aoooo(nO,nO,nO,nO),bvvvv(nV,nV,nV,nV),hovvo(nO,nV,nV,nO), &
           r1(nO,nV),r2(nO,nO,nV,nV))

! Memory allocation for DIIS

  allocate(err1_diis(nO*nV      ,max_diis),t1_diis(nO*nV      ,max_diis), & 
           err2_diis(nO*nO*nV*nV,max_diis),t2_diis(nO*nO*nV*nV,max_diis))

  Conv = 1d0
  nSCF = 0
  ECCSD = ERHF

  n_diis         = 0
  t1_diis(:,:)   = 0d0
  t2_diis(:,:)   = 0d0
  err1_diis(:,:) = 0d0
  err2_diis(:,:) = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------
  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'| CCSD calculation                                  |'
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A1,1X,A3,1X,A1,1X,A16,1X,A1,1X,A10,1X,A1,1X,A10,1X,A1,1X)') &
            '|','#','|','E(CCSD)','|','Ec(CCSD)','|','Conv','|'
  write(*,*)'----------------------------------------------------'

  do while(Conv > thresh .and. nSCF < maxSCF)

!   Increment 

    nSCF = nSCF + 1

!   Scuseria Eqs. (5), (6) and (7)

    call form_h(nC,nO,nV,nR,eO,eV,OOVV,t1,tau,hvv,hoo,hvo)

!   Scuseria Eqs. (9), (10), (11), (12) and (13)

    call form_g(nC,nO,nV,nR,hvv,hoo,VOVV,OOOV,t1,gvv,goo)

    call form_abh(nC,nO,nV,nR,OOOO,OVOO,OOVV,VVVV,VOVV,OVVO,OVVV,t1,tau,aoooo,bvvvv,hovvo)

!   Compute residuals

    call form_r1(nC,nO,nV,nR,OVVO,OVVV,OOOV,hvv,hoo,hvo,t1,t2,tau,r1)

    call form_r2(nC,nO,nV,nR,OOVV,OVOO,OVVV,OVVO,gvv,goo,aoooo,bvvvv,hovvo,t1,t2,tau,r2)

!   Check convergence 

    Conv = max(maxval(abs(r1(nC+1:nO,1:nV-nR))),maxval(abs(r2(nC+1:nO,nC+1:nO,1:nV-nR,1:nV-nR))))

!   Update 

    t1(:,:)      = t1(:,:)     - r1(:,:)    /delta_OV  (:,:)
    t2(:,:,:,:)  = t2(:,:,:,:) - r2(:,:,:,:)/delta_OOVV(:,:,:,:)

    call form_tau(nC,nO,nV,nR,t1,t2,tau)
 
!   Compute correlation energy

    call CCSD_correlation_energy(nC,nO,nV,nR,OOVV,tau,EcCCSD) 

!   Dump results

    ECCSD = ERHF + EcCCSD

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
!   call DIIS_extrapolation(rcond1,nO*nV      ,nO*nV      ,n_diis,err1_diis,t1_diis,-r1/delta_OV  ,t1)
!   call DIIS_extrapolation(rcond2,nO*nO*nV*nV,nO*nO*nV*nV,n_diis,err2_diis,t2_diis,-r2/delta_OOVV,t2)

    !  Reset DIIS if required

!   if(min(abs(rcond1),abs(rcond2)) < 1d-15) n_diis = 0

    write(*,'(1X,A1,1X,I3,1X,A1,1X,F16.10,1X,A1,1X,F10.6,1X,A1,1X,F10.6,1X,A1,1X)') &
      '|',nSCF,'|',ECCSD+ENuc,'|',EcCCSD,'|',Conv,'|'

  end do
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

  end if

  write(*,*)
  write(*,*)'----------------------------------------------------'
  write(*,*)'                 CCSD energy                        '
  write(*,*)'----------------------------------------------------'
  write(*,'(1X,A20,1X,F15.10)')' E(CCSD)  = ',ENuc+ECCSD  
  write(*,'(1X,A20,1X,F10.6)')' Ec(CCSD) = ',EcCCSD 
  write(*,*)'----------------------------------------------------'
  write(*,*)

! Deallocate memory

  deallocate(hvv,hoo,hvo,         &
             delta_OV,delta_OOVV, &
             gvv,goo,             &
             aoooo,bvvvv,hovvo,   &
             tau,                 &
             r1,r2)

!------------------------------------------------------------------------
! (T) correction
!------------------------------------------------------------------------
  if(doCCSDT) then

    call cpu_time(start_CCSDT)
    call CCSDT(nC,nO,nV,nR,eO,eV,OOVV,VVVO,VOOO,t1,t2,EcCCT)
    call cpu_time(end_CCSDT)

     t_CCSDT = end_CCSDT - start_CCSDT
     write(*,'(A65,1X,F9.3,A8)') 'Total CPU time for (T) = ',t_CCSDT,' seconds'
     write(*,*)

    write(*,*)
    write(*,*)'----------------------------------------------------'
    write(*,*)'                 CCSD(T) energy                     '
    write(*,*)'----------------------------------------------------'
    write(*,'(1X,A20,1X,F15.10)')' E(CCSD(T))  = ',ENuc + ECCSD + EcCCT
    write(*,'(1X,A20,1X,F10.6)')' Ec(CCSD(T)) = ',EcCCSD + EcCCT
    write(*,*)'----------------------------------------------------'
    write(*,*)

  end if

end subroutine CCSD
