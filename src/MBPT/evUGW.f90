subroutine evUGW(maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,COHSEX,BSE,TDA_W,TDA,    & 
                G0W,GW0,dBSE,dTDA,evDyn,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS,ENuc, &
                ERHF,Hc,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,dipole_int_aa,dipole_int_bb,PHF,cHF,eHF,eG0W0)

! Perform self-consistent eigenvalue-only GW calculation

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
  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: G0W
  logical,intent(in)            :: GW0
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  double precision,intent(in)   :: eta

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)

  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: PHF(nBas,nBas,nspin) 
  double precision,intent(in)   :: eG0W0(nBas,nspin)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_abab(nBas,nBas,nBas,nBas)
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
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM
  double precision              :: alpha
  double precision,allocatable  :: error_diis(:,:,:)
  double precision,allocatable  :: e_diis(:,:,:)
  double precision,allocatable  :: eGW(:,:)
  double precision,allocatable  :: eOld(:,:)
  double precision,allocatable  :: Z(:,:)
  integer                       :: nS_aa,nS_bb,nS_sc
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:,:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|       Self-consistent evGW calculation       |'
  write(*,*)'************************************************'
  write(*,*)

! COHSEX approximation

  if(COHSEX) then 
    write(*,*) 'COHSEX approximation activated!'
    write(*,*)
  end if

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
    write(*,*)
  end if

! TDA 

  if(TDA) then 
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! GW0

  if(GW0) then 
    write(*,*) 'GW0 scheme activated!'
    write(*,*)
  end if

! G0W

  if(G0W) then 
    write(*,*) 'G0W scheme activated!'
    write(*,*)
  end if

! Linear mixing

  linear_mixing = .false.
  alpha = 0.2d0

! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  allocate(eGW(nBas,nspin),eOld(nBas,nspin),Z(nBas,nspin),SigC(nBas,nspin),OmRPA(nS_sc),XpY_RPA(nS_sc,nS_sc), &
           XmY_RPA(nS_sc,nS_sc),rho_RPA(nBas,nBas,nS_sc,nspin),error_diis(nBas,max_diis,nspin),e_diis(nBas,max_diis,nspin))

! Initialization

  nSCF              = 0
  ispin             = 1
  n_diis            = 0
  Conv              = 1d0
  e_diis(:,:,:)     = 0d0
  error_diis(:,:,:) = 0d0
  eGW(:,:)          = eG0W0(:,:)
  eOld(:,:)         = eGW(:,:)
  Z(:,:)            = 1d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

   ! Compute screening

    if(.not. GW0 .or. nSCF == 0) then

      call unrestricted_linear_response(ispin,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
                                        eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)


    endif

    !----------------------!
    ! Excitation densities !
    !----------------------!

    call unrestricted_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

    !------------------------------------------------!
    ! Compute self-energy and renormalization factor !
    !------------------------------------------------!

    if(G0W) then

      call unrestricted_self_energy_correlation_diag(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,OmRPA,rho_RPA,SigC)
      call unrestricted_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,OmRPA,rho_RPA,Z)

    else 

      call unrestricted_self_energy_correlation_diag(eta,nBas,nC,nO,nV,nR,nS_sc,eGW,OmRPA,rho_RPA,SigC)
      call unrestricted_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS_sc,eGW,OmRPA,rho_RPA,Z)

    endif

    !-----------------------------------!
    ! Solve the quasi-particle equation !
    !-----------------------------------!

    eGW(:,:) = eHF(:,:) + SigC(:,:)

    ! Convergence criteria

    Conv = maxval(abs(eGW(:,:) - eOld(:,:)))

    ! Print results

    call print_evUGW(nBas,nO,nSCF,Conv,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA)

    ! Linear mixing or DIIS extrapolation

    if(linear_mixing) then
 
      eGW(:,:) = alpha*eGW(:,:) + (1d0 - alpha)*eOld(:,:)
 
    else

      n_diis = min(n_diis+1,max_diis)
      do is=1,nspin
        call DIIS_extrapolation(rcond(ispin),nBas,nBas,n_diis,error_diis(:,:,is), & 
                                e_diis(:,:,is),eGW(:,is)-eOld(:,is),eGW(:,is))
      end do

!    Reset DIIS if required

      if(minval(rcond(:)) < 1d-15) n_diis = 0

    endif

    ! Save quasiparticles energy for next cycle

    eOld(:,:) = eGW(:,:)

    ! Increment

    nSCF = nSCF + 1

  enddo
!------------------------------------------------------------------------
! End main loop
!------------------------------------------------------------------------

! Plot stuff

! call plot_GW(nBas,nC,nO,nV,nR,nS,eHF,eGW,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin))

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

  deallocate(eOld,Z,SigC,OmRPA,XpY_RPA,XmY_RPA,rho_RPA,error_diis,e_diis)

! Perform BSE calculation

  if(BSE) then

    call unrestricted_Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS, &
                                     ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,dipole_int_aa,dipole_int_bb,eGW,eGW,EcBSE)

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

!     call ACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,eGW,eGW,EcAC)

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy (singlet) =',EcAC(1)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy (triplet) =',EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW correlation energy           =',EcAC(1) + EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@evGW total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

  endif

end subroutine evUGW