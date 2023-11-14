subroutine evGGF2(dotest,dophBSE,doppBSE,TDA,dBSE,dTDA,maxSCF,thresh,max_diis, &
                  linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform eigenvalue self-consistent second-order Green function calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  integer,intent(in)            :: maxSCF
  double precision,intent(in)   :: thresh
  integer,intent(in)            :: max_diis
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
  double precision              :: EcBSE
  double precision              :: Conv
  double precision              :: rcond
  double precision,allocatable  :: eGF(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Generalized evGF2 Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! Memory allocation

  allocate(SigC(nBas),Z(nBas),eGF(nBas),eOld(nBas),error_diis(nBas,max_diis),e_diis(nBas,max_diis))

! Initialization

  Conv            = 1d0
  nSCF            = 0
  n_diis          = 0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGF(:)         = eHF(:)
  eOld(:)         = eHF(:)
  rcond           = 0d0

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Frequency-dependent second-order contribution

    if(regularize) then 

!     call GGF2_reg_self_energy_diag(eta,nBas,nC,nO,nV,nR,eGF,ERI,SigC,Z)

    else

      call GGF2_self_energy_diag(eta,nBas,nC,nO,nV,nR,eGF,ERI,SigC,Z)

    end if

    ! Solve the quasi-particle equation
 
    if(linearize) then
 
    else
 
      write(*,*) ' *** Quasiparticle energies obtained by root search (experimental) *** '
      write(*,*)
 
      call GGF2_QP_graph(eta,nBas,nC,nO,nV,nR,eHF,ERI,eOld,eOld,eGF,Z)
 
    end if

    Conv = maxval(abs(eGF - eOld))

    ! Print results

    call GMP2(.false.,regularize,nBas,nC,nO,nV,nR,ERI,ENuc,ERHF,eGF,Ec)
    call print_evGF2(nBas,nO,nSCF,Conv,eHF,SigC,Z,eGF,ENuc,ERHF,Ec)

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call DIIS_extrapolation(rcond,nBas,nBas,n_diis,error_diis,e_diis,eGF-eOld,eGF)

    if(abs(rcond) < 1d-15) n_diis = 0

    eOld(:) = eGF(:)

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

  if(dophBSE) then 
  
    call GGF2_phBSE2(TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@evGGF2 correlation energy          =',EcBSE,' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@evGGF2 total energy                =',ENuc + ERHF + EcBSE,' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Perform ppBSE2 calculation

! if(doppBSE) then

!   call GGF2_ppBSE2(TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,ERI,dipole_int,eGF,EcBSE)

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGGF2 correlation energy          =',EcBSE),' au'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGGF2 total energy                =',ENuc + ERHF + EcBSE,' au'
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

! end if

! Testing zone

  if(dotest) then

    call dump_test_value('G','evGF2 correlation energy',Ec)
    call dump_test_value('G','evGF2 HOMO energy',eGF(nO))
    call dump_test_value('G','evGF2 LUMO energy',eGF(nO+1))

  end if

end subroutine 
