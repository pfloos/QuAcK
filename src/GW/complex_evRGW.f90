subroutine complex_evRGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE, & 
                 singlet,triplet,linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform self-consistent eigenvalue-only GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: ENuc
  complex*16,intent(in)         :: ERHF
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  complex*16,intent(in)         :: eHF(nOrb)
  complex*16,intent(in)         :: ERI(nOrb,nOrb,nOrb,nOrb)
  complex*16,intent(in)         :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  logical                       :: dRPA = .true.
  integer                       :: ispin
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: flow
  double precision              :: rcond
  double precision              :: Conv
  complex*16                    :: EcRPA
  complex*16                    :: EcBSE(nspin)
  complex*16                    :: EcGM
  double precision              :: alpha
  complex*16,allocatable        :: Aph(:,:)
  complex*16,allocatable        :: Bph(:,:)
  complex*16,allocatable        :: error_diis(:,:)
  complex*16,allocatable        :: e_diis(:,:)
  complex*16,allocatable        :: eGW(:)
  complex*16,allocatable        :: eOld(:)
  double precision,allocatable  :: Re_eGW(:)
  double precision,allocatable  :: Im_eGW(:)
  double precision,allocatable  :: Re_eOld(:)
  double precision,allocatable  :: Im_eOld(:)
  double precision,allocatable  :: Re_eHF(:)
  double precision,allocatable  :: Im_eHF(:)
  double precision,allocatable  :: Re_Z(:)
  double precision,allocatable  :: Im_Z(:)
  double precision,allocatable  :: Re_SigC(:)
  double precision,allocatable  :: Im_SigC(:)
  complex*16,allocatable        :: Om(:)
  complex*16,allocatable        :: XpY(:,:)
  complex*16,allocatable        :: XmY(:,:)
  complex*16,allocatable        :: rho(:,:,:)

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted evGW Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
    write(*,*)
  end if

! SRG regularization

  flow = 100d0

  if(doSRG) then

    write(*,*) '*** SRG regularized evGW scheme ***'
    write(*,*)

  end if

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),Re_eGW(nOrb),Im_eGW(nOrb),eGW(nOrb),Re_eOld(nOrb),Im_eOld(nOrb),&
           eOld(nOrb),Re_Z(nOrb),Im_Z(nOrb),Re_SigC(nOrb),Im_SigC(nOrb), & 
           Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nOrb,nOrb,nS),error_diis(nOrb,max_diis),e_diis(nOrb,max_diis),&
           Re_eHF(nOrb),Im_eHF(nOrb))

! Initialization

  nSCF            = 0
  ispin           = 1
  n_diis          = 0
  Conv            = 1d0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  Re_eHF(:)       = real(eHF(:))
  Im_eHF(:)       = aimag(eHF(:))
  eGW(:)          = eHF(:)
  Re_eGW(:)       = Re_eHF(:)
  Im_eGW(:)       = Im_eHF(:)
  eOld(:)         = eGW(:)
  Re_eOld(:)      = Re_eGW(:)
  Im_eOld(:)      = Im_eGW(:)
  Re_Z(:)         = 1d0
  Im_Z(:)         = 0d0
  rcond           = 0d0
  
!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

   ! Compute screening

                   call complex_phRLR_A(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
    if(.not.TDA_W) call complex_phRLR_B(ispin,dRPA,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

    call complex_phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

    ! Compute spectral weights

    call complex_RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

    ! Compute correlation part of the self-energy 
    ! Implement here the srg version if(doSRG) .. complex_RGW_SRG_self_energy_diag
    call complex_RGW_self_energy_diag(eta,nBas,nOrb,nC,nO,nV,nR,nS,Re_eGW,Im_eGW,Om,rho,&
            EcGM,Re_SigC,Im_SigC,Re_Z,Im_Z)
      
    ! Solve the quasi-particle equation

    if(linearize) then 
 
      write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
      write(*,*)

      Re_eGW(:) = Re_eHF(:) + Re_SigC(:)
      Im_eGW(:) = Im_eHF(:) + Im_SigC(:)
      eGW = cmplx(Re_eGW,Im_eGW,kind=8)
    else 

      write(*,*) ' *** Quasiparticle energies obtained by root search *** '
      write(*,*)
  
      call complex_RGW_QP_graph(doSRG,eta,flow,nOrb,nC,nO,nV,nR,nS,Re_eHF,Im_eHF,Om,&
              rho,Re_eOld,Im_eOld,Re_eOld,Im_eOld,Re_eGW,Im_eGW,Re_Z,Im_Z)
 
      eGW = cmplx(Re_eGW,Im_eGW,kind=8)
    end if

    ! Convergence criteria

    Conv = maxval(abs(eGW - eOld))

    ! Print results

    call print_complex_evRGW(nBas,nO,nSCF,Conv,Re_eHF,Im_eHF,ENuc,ERHF,Re_SigC,Im_SigC,&
            Re_Z,Im_Z,Re_eGW,Im_eGW,EcRPA,EcGM)
    ! DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      call complex_DIIS_extrapolation(rcond,nOrb,nOrb,n_diis,error_diis,e_diis,eGW-eOld,eGW)

    end if

    ! Save quasiparticles energy for next cycle
    Re_eGW(:) = real(eGW(:))
    Im_eGW(:) = aimag(eGW(:))
    eOld(:) = eGW(:)
    Re_eOld(:) = real(eOld(:))
    Im_eOld(:) = aimag(eOld(:))

    ! Increment

    nSCF = nSCF + 1

  end do
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

  end if

!--------------------!
! Cumulant expansion !
!--------------------!
!
!! call RGWC(dotest,eta,nOrb,nC,nO,nV,nR,nS,Om,rho,eHF,eGW,eGW,Z)
!
!! Deallocate memory
!
!  deallocate(Aph,Bph,eOld,Z,SigC,Om,XpY,XmY,rho,error_diis,e_diis)
!
!! Perform BSE calculation
!
!  if(dophBSE) then
!
!    call RGW_phBSE(dophBSE2,exchange_kernel,TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta, & 
!                   nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eGW,eGW,EcBSE)
!
!    write(*,*)
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@RHF correlation energy           = ',sum(EcBSE),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,*)
!
!!   Compute the BSE correlation energy via the adiabatic connection 
!
!    if(doACFDT) then
!
!      call RGW_phACFDT(exchange_kernel,doXBS,TDA_W,TDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,eGW,eGW,EcBSE)
!
!      write(*,*)
!      write(*,*)'-------------------------------------------------------------------------------'
!      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
!      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
!      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@RHF correlation energy           = ',sum(EcBSE),' au'
!      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
!      write(*,*)'-------------------------------------------------------------------------------'
!      write(*,*)
!
!    end if
!
!  end if
!
!  if(doppBSE) then
!
!    call RGW_ppBSE(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE)
!
!    write(*,*)
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGW@RHF correlation energy (singlet) = ',EcBSE(1),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGW@RHF correlation energy (triplet) = ',EcBSE(2),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGW@RHF correlation energy           = ',sum(EcBSE),' au'
!    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@evGW@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
!    write(*,*)'-------------------------------------------------------------------------------'
!    write(*,*)
!
!  end if
!
!! Testing zone
!
!  if(dotest) then
!
!    call dump_test_value('R','evGW correlation energy',EcRPA)
!    call dump_test_value('R','evGW HOMO energy',eGW(nO))
!    call dump_test_value('R','evGW LUMO energy',eGW(nO+1))
!
!  end if
!
end subroutine 
