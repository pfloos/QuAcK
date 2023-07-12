subroutine evGTeh(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,BSE,BSE2,TDA_T,TDA,dBSE,dTDA,evDyn,ppBSE, & 
                  singlet,triplet,linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_AO,ERI_MO,dipole_int,PHF,  & 
                  cHF,eHF,Vxc,eG0T0)

! Perform self-consistent eigenvalue-only ehGT calculation

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: BSE
  logical,intent(in)            :: BSE2
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: ppBSE
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: Vxc(nBas)
  double precision,intent(in)   :: eG0T0(nBas)
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_MO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  logical                       :: linear_mixing
  integer                       :: ispin
  integer                       :: nSCF
  integer                       :: n_diis
  integer                       :: i,a,jb,p
  double precision              :: rcond
  double precision              :: Conv
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcppBSE(nspin)
  double precision              :: EcGM
  double precision              :: alpha
  double precision              :: Dpijb,Dpajb
  double precision,allocatable  :: error_diis(:,:)
  double precision,allocatable  :: e_diis(:,:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eOld(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: SigX(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rhoL_RPA(:,:,:)
  double precision,allocatable  :: rhoR_RPA(:,:,:)
  
  double precision,allocatable  :: eGTlin(:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|    Self-consistent evGTeh calculation        |'
  write(*,*)'************************************************'
  write(*,*)

! TDA for T

  if(TDA_T) then 
    write(*,*) 'Tamm-Dancoff approximation for eh T-matrix!'
    write(*,*)
  end if

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Linear mixing

  linear_mixing = .false.
  alpha = 0.2d0

! Memory allocation

  allocate(eGT(nBas),eOld(nBas),Z(nBas),SigX(nBas),SigC(nBas),OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS), & 
           rhoL_RPA(nBas,nBas,nS),rhoR_RPA(nBas,nBas,nS),error_diis(nBas,max_diis),e_diis(nBas,max_diis),eGTlin(nBas))

! Compute the exchange part of the self-energy

  call self_energy_exchange_diag(nBas,cHF,PHF,ERI_AO,SigX)

! Initialization

  nSCF            = 0
  ispin           = 2
  n_diis          = 0
  Conv            = 1d0
  e_diis(:,:)     = 0d0
  error_diis(:,:) = 0d0
  eGT(:)          = eG0T0(:)
  eOld(:)         = eGT(:)
  Z(:)            = 1d0
  rcond           = 0d0


  
!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

   ! Compute screening

    call phLR(ispin,.false.,TDA_T,eta,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI_MO,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

   ! Compute spectral weights

    call GTeh_excitation_density(nBas,nC,nO,nR,nS,ERI_MO,XpY_RPA,XmY_RPA,rhoL_RPA,rhoR_RPA)

    ! Compute correlation part of the self-energy 

    if(regularize) then 

!     call regularized_self_energy_correlation_diag(COHSEX,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA,rho_RPA,EcGM,SigC)
!     call renormalization_factor_SRG(eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA,rho_RPA,Z)

    else

      call GTeh_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS,eGT,OmRPA,rhoL_RPA,rhoR_RPA,EcGM,SigC)
      call GTeh_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS,eGT,OmRPA,rhoL_RPA,rhoR_RPA,Z)

    end if

    ! Solve the quasi-particle equation

    eGTlin(:) = eHF(:) + SigX(:) + SigC(:) - Vxc(:)

    ! Linearized or graphical solution?

    if(linearize) then 
 
       write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
       write(*,*)

       eGT(:) = eGTlin(:)

    else 

!      write(*,*) ' *** Quasiparticle energies obtained by root search (experimental) *** '
!      write(*,*)
! 
!      call QP_graph(nBas,nC,nO,nV,nR,nS,eta,eHF,SigX,Vxc,OmRPA,rho_RPA,eGWlin,eGW,regularize)
 
    end if

    ! Convergence criteria

    Conv = maxval(abs(eGT - eOld))

    ! Print results

    call print_evGTeh(nBas,nO,nSCF,Conv,eHF,ENuc,ERHF,SigC,Z,eGT,EcRPA,EcGM)

    ! Linear mixing or DIIS extrapolation

    if(linear_mixing) then
 
      eGT(:) = alpha*eGT(:) + (1d0 - alpha)*eOld(:)
 
    else

      n_diis = min(n_diis+1,max_diis)
      if(abs(rcond) > 1d-7) then
        call DIIS_extrapolation(rcond,nBas,nBas,n_diis,error_diis,e_diis,eGT-eOld,eGT)
      else
        n_diis = 0
      end if

    end if

    ! Save quasiparticles energy for next cycle

    eOld(:) = eGT(:)

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

! Deallocate memory

  deallocate(eOld,Z,SigC,OmRPA,XpY_RPA,XmY_RPA,rhoL_RPA,rhoR_RPA,error_diis,e_diis)

! Perform BSE calculation

! if(BSE) then

!   call Bethe_Salpeter(BSE2,TDA_T,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eGW,eGW,EcBSE)

!   if(exchange_kernel) then

!     EcBSE(1) = 0.5d0*EcBSE(1)
!     EcBSE(2) = 1.5d0*EcBSE(2)

!   end if

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW correlation energy (singlet) =',EcBSE(1)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW correlation energy (triplet) =',EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW correlation energy           =',EcBSE(1) + EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@evGW total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

!   if(doACFDT) then

!     write(*,*) '------------------------------------------------------'
!     write(*,*) 'Adiabatic connection version of BSE correlation energy'
!     write(*,*) '------------------------------------------------------'
!     write(*,*)

!     if(doXBS) then

!       write(*,*) '*** scaled screening version (XBS) ***'
!       write(*,*)

!     end if

!     call ACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,eGW,eGW,EcAC)

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy (singlet) =',EcAC(1)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy (triplet) =',EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy           =',EcAC(1) + EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

! end if

! if(ppBSE) then

!   call Bethe_Salpeter_pp(TDA_W,TDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eHF,eGW,EcppBSE)

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 correlation energy (singlet) =',EcppBSE(1)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 correlation energy (triplet) =',3d0*EcppBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 correlation energy =',EcppBSE(1) + 3d0*EcppBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@ppBSE@G0W0 total energy =',ENuc + ERHF + EcppBSE(1) + 3d0*EcppBSE(2)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

!   nBas2 = 2*nBas
!   nO2   = 2*nO
!   nV2   = 2*nV
!   nC2   = 2*nC
!   nR2   = 2*nR
!   nS2   = nO2*nV2
!
!   allocate(seHF(nBas2),seGW(nBas2),sERI(nBas2,nBas2,nBas2,nBas2))
!
!   call spatial_to_spin_MO_energy(nBas,eHF,nBas2,seHF)
!   call spatial_to_spin_MO_energy(nBas,eGW,nBas2,seGW)
!   call spatial_to_spin_ERI(nBas,ERI_MO,nBas2,sERI)
!
!   call Bethe_Salpeter_pp_so(TDA_W,TDA,singlet,triplet,eta,nBas2,nC2,nO2,nV2,nR2,nS2,sERI,dipole_int,seHF,seGW,EcppBSE)

! end if

end subroutine 
