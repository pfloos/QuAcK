subroutine complex_evRGF2(dotest,dophBSE,doppBSE,TDA,dBSE,dTDA,maxSCF,thresh,max_diis,singlet,triplet, &
                 linearize,eta,regularize,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

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
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nO
  integer,intent(in)            :: nC
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  complex*16,intent(in)         :: ERHF
  complex*16,intent(in)         :: eHF(nOrb)
  complex*16,intent(in)         :: ERI(nOrb,nOrb,nOrb,nOrb)
  complex*16,intent(in)         :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: Ec
  double precision              :: flow
  double precision              :: EcBSE(nspin)
  double precision              :: Conv
  double precision              :: rcond
  double precision,allocatable  :: Re_eHF(:)
  double precision,allocatable  :: Im_eHF(:)
  complex*16,allocatable        :: eGF(:)
  double precision,allocatable  :: Re_eGF(:)
  double precision,allocatable  :: Im_eGF(:)
  complex*16,allocatable        :: eOld(:)
  double precision,allocatable  :: Re_eOld(:)
  double precision,allocatable  :: Im_eOld(:)
  double precision,allocatable  :: Re_SigC(:)
  double precision,allocatable  :: Im_SigC(:)
  double precision,allocatable  :: Re_Z(:)
  double precision,allocatable  :: Im_Z(:)
  complex*16,allocatable        :: error_diis(:,:)
  complex*16,allocatable        :: e_diis(:,:)

! Hello world


  write(*,*)
  write(*,*)'********************************'
  write(*,*)'* Restricted evGF2 Calculation *'
  write(*,*)'********************************'
  write(*,*)

! Memory allocation

  allocate(Re_SigC(nOrb),Im_SigC(nOrb), Re_Z(nOrb),Im_Z(nOrb), eGF(nOrb),&
          Re_eGF(nOrb),Im_eGF(nOrb), eOld(nOrb),Re_eOld(nOrb),Im_eOld(nOrb),&
          error_diis(nOrb,max_diis), e_diis(nOrb,max_diis),Re_eHF(nOrb),Im_eHF(nOrb))

! Initialization

  Conv            = 1d0
  nSCF            = 0
  n_diis          = 0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  Re_eHF(:)       = real(eHF(:))
  Im_eHF(:)       = aimag(eHF(:))
  eGF(:)         = eHF(:)
  Re_eGF(:)       = Re_eHF(:)
  Im_eGF(:)       = Im_eHF(:)
  eOld(:)         = eHF(:)
  Re_eOld(:)    = Re_eHF(:)
  Im_eOld(:)    = Im_eHF(:)
  rcond           = 0d0
  Re_Z(:)       = 0d0
  Im_Z(:)       = 0d0
  flow = 500d0 

!------------------------------------------------------------------------
! Main SCF loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF < maxSCF)

    ! Frequency-dependent second-order contribution

    if(regularize) then 
      
      call complex_cRGF2_reg_self_energy_diag(flow,eta,nOrb,nC,nO,nV,nR,Re_eGF,Im_eGF,ERI,Re_SigC,Im_SigC,Re_Z,Im_Z)
      
    else

      call complex_cRGF2_self_energy_diag(eta,nOrb,nC,nO,nV,nR,Re_eGF,Im_eGF,ERI,Re_SigC,Im_SigC,Re_Z,Im_Z)

    end if

    ! Solve the quasi-particle equation
 
    if(linearize) then
      Re_eGF(:) = Re_eHF(:) + Re_SigC(:)
      Im_eGF(:) = Im_eHF(:) + Im_SigC(:)
      eGF = cmplx(Re_eGF,Im_eGF,kind=8)
    else
 
      write(*,*) ' *** Quasiparticle energies obtained by root search *** '
      write(*,*)
 
      call complex_cRGF2_QP_graph(flow,regularize,eta,nOrb,nC,nO,nV,nR,Re_eHF,Im_eHF,&
              ERI,Re_eOld,Im_eOld,Re_eOld,Im_eOld,Re_eGF,Im_eGF,Re_Z,Im_Z)
      eGF = cmplx(Re_eGF,Im_eGF,kind=8)
    end if

    Conv = maxval(abs(eGF - eOld))

    ! Print results

    !call RMP2(.false.,regularize,nOrb,nC,nO,nV,nR,ERI,ENuc,ERHF,eGF,Ec)
    call print_complex_evRGF2(nOrb,nO,nSCF,Conv,Re_eHF,Im_eHF,ENuc,ERHF,Re_SigC,Im_SigC,Re_Z,Im_Z,Re_eGF,Im_eGF)

    ! DIIS extrapolation

    n_diis = min(n_diis+1,max_diis)
    call complex_DIIS_extrapolation(rcond,nOrb,nOrb,n_diis,error_diis,e_diis,eGF-eOld,eGF)

    if(abs(rcond) < 1d-15) n_diis = 0

    Re_eGF(:) = real(eGF(:))
    Im_eGF(:) = aimag(eGF(:))
    eOld(:) = eGF(:)
    Re_eOld(:) = real(eOld(:))
    Im_eOld(:) = aimag(eOld(:))

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

!! Perform BSE@GF2 calculation
!
!  if(dophBSE) then 
!  
!    call RGF2_phBSE(TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eGF,EcBSE)
!
!    write(*,*)
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@evGF2 correlation energy (singlet) =',EcBSE(1)
!    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@evGF2 correlation energy (triplet) =',EcBSE(2)
!    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@evGF2 correlation energy           =',sum(EcBSE(:))
!    write(*,'(2X,A50,F20.10)') 'Tr@phBSE@evGF2 total energy                 =',ENuc + ERHF + sum(EcBSE(:))
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,*)
!
!  end if
!
!! Perform ppBSE@GF2 calculation
!
!  if(doppBSE) then
!
!    call RGF2_ppBSE(TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,ERI,dipole_int,eGF,EcBSE)
!
!    write(*,*)
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGF2 correlation energy (singlet) =',EcBSE(1),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGF2 correlation energy (triplet) =',3d0*EcBSE(2),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGF2 correlation energy           =',EcBSE(1) + 3d0*EcBSE(2),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGF2 total energy                 =',ENuc + ERHF + EcBSE(1) + 3d0*EcBSE(2),' au'
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,*)
!
!  end if
!
!! Testing zone
!
!  if(dotest) then
!
!    call dump_test_value('R','evGF2 correlation energy',Ec)
!    call dump_test_value('R','evGF2 HOMO energy',eGF(nO))
!    call dump_test_value('R','evGF2 LUMO energy',eGF(nO+1))
!
!  end if

  deallocate(Re_SigC,Im_SigC, Re_Z,Im_Z, eGF,Re_eGF,Im_eGF, eOld,Re_eOld,Im_eOld, error_diis, e_diis)

end subroutine 
