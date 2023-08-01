subroutine evUGF2(maxSCF,thresh,max_diis,BSE,TDA,dBSE,dTDA,spin_conserved,spin_flip, &
                  eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EUHF,ERI_aaaa,ERI_aabb,ERI_bbbb,   & 
                  dipole_int_aa,dipole_int_bb,cHF,eHF)

! Perform self-consistent eigenvalue-only GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)

  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  logical                       :: linear_mixing
  integer                       :: is
  integer                       :: ispin
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: rcond(nspin)
  double precision              :: Conv
  double precision              :: Ec(nsp)
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: alpha
  double precision,allocatable  :: error_diis(:,:,:)
  double precision,allocatable  :: e_diis(:,:,:)
  double precision,allocatable  :: eGF2(:,:)
  double precision,allocatable  :: eOld(:,:)
  double precision,allocatable  :: Z(:,:)
  integer                       :: nS_aa,nS_bb,nS_sc
  double precision,allocatable  :: SigC(:,:)

! Hello world

  write(*,*)
  write(*,*)'**************************************************'
  write(*,*)'| Self-consistent unrestricted evGF2 calculation |'
  write(*,*)'**************************************************'
  write(*,*)

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Linear mixing

  linear_mixing = .false.
  alpha = 0.2d0

! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  allocate(eGF2(nBas,nspin),eOld(nBas,nspin),Z(nBas,nspin),SigC(nBas,nspin), &
           error_diis(nBas,max_diis,nspin),e_diis(nBas,max_diis,nspin))

! Initialization

  nSCF              = 0
  ispin             = 1
  n_diis            = 0
  Conv              = 1d0
  e_diis(:,:,:)     = 0d0
  error_diis(:,:,:) = 0d0
  eGF2(:,:)         = eHF(:,:)
  eOld(:,:)         = eHF(:,:)
  Z(:,:)            = 1d0
  rcond(:)          = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

    !------------------------------------------------!
    ! Compute self-energy and renormalization factor !
    !------------------------------------------------!

    if(regularize) then

      call UGF2_reg_self_energy_diag(nBas,nC,nO,nV,nR,eta,ERI_aaaa,ERI_aabb,ERI_bbbb,eHF,eGF2,SigC,Z)

    else

      call UGF2_self_energy_diag(nBas,nC,nO,nV,nR,eta,ERI_aaaa,ERI_aabb,ERI_bbbb,eHF,eGF2,SigC,Z)

    end if

    !-----------------------------------!
    ! Solve the quasi-particle equation !
    !-----------------------------------!

    eGF2(:,:) = eHF(:,:) + SigC(:,:)

    ! Convergence criteria

    Conv = maxval(abs(eGF2(:,:) - eOld(:,:)))

    ! Compute MP2 correlation energy

    call UMP2(nBas,nC,nO,nV,nR,ERI_aaaa,ERI_aabb,ERI_bbbb,ENuc,EUHF,eGF2,Ec)

    ! Print results
    
    call print_evUGF2(nBas,nO,nSCF,Conv,eHF,ENuc,EUHF,SigC,Z,eGF2,Ec)

    ! Linear mixing or DIIS extrapolation

    if(linear_mixing) then
 
      eGF2(:,:) = alpha*eGF2(:,:) + (1d0 - alpha)*eOld(:,:)
 
    else

      n_diis = min(n_diis+1,max_diis)
      do is=1,nspin
        call DIIS_extrapolation(rcond(ispin),nBas,nBas,n_diis,error_diis(:,1:n_diis,is), & 
                                e_diis(:,1:n_diis,is),eGF2(:,is)-eOld(:,is),eGF2(:,is))
      end do

!    Reset DIIS if required

      if(minval(rcond(:)) < 1d-15) n_diis = 0

    endif

    ! Save quasiparticles energy for next cycle

    eOld(:,:) = eGF2(:,:)

    ! Increment

    nSCF = nSCF + 1

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

  deallocate(eOld,Z,SigC,error_diis,e_diis)

! Perform BSE calculation

  if(BSE) then

    print*,'!!! BSE2 NYI for evUGF2 !!!'

  endif

end subroutine 
