subroutine evUGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,BSE,TDA_W,TDA,dBSE,dTDA, & 
                 spin_conserved,spin_flip,linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,   &
                 EUHF,S,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,cHF,eHF)

! Perform self-consistent eigenvalue-only GW calculation

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  integer,intent(in)            :: maxSCF
  integer,intent(in)            :: max_diis
  double precision,intent(in)   :: thresh
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
  logical,intent(in)            :: linearize
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
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  logical                       :: dRPA
  logical                       :: linear_mixing
  integer                       :: is
  integer                       :: ispin
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: rcond(nspin)
  double precision              :: Conv
  double precision              :: EcRPA
  double precision              :: EcGM(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: alpha
  double precision,allocatable  :: error_diis(:,:,:)
  double precision,allocatable  :: e_diis(:,:,:)
  double precision,allocatable  :: eGW(:,:)
  double precision,allocatable  :: eOld(:,:)
  double precision,allocatable  :: Z(:,:)
  integer                       :: nSa,nSb,nSt
  double precision,allocatable  :: SigC(:,:)

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:,:)

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'| Unrestricted evGW Calculation *'
  write(*,*)'*********************************'
  write(*,*)

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

! Initialization

  EcRPA = 0d0
  dRPA = .true.

! Linear mixing

  linear_mixing = .false.
  alpha = 0.2d0

! Memory allocation

  nSa = nS(1)
  nSb = nS(2)
  nSt = nSa + nSb

  allocate(eGW(nBas,nspin),eOld(nBas,nspin),Z(nBas,nspin),SigC(nBas,nspin), &
           Aph(nSt,nSt),Bph(nSt,nSt),Om(nSt),XpY(nSt,nSt),XmY(nSt,nSt),     &
           rho(nBas,nBas,nSt,nspin),error_diis(nBas,max_diis,nspin),e_diis(nBas,max_diis,nspin))

! Initialization

  nSCF              = 0
  ispin             = 1
  n_diis            = 0
  Conv              = 1d0
  e_diis(:,:,:)     = 0d0
  error_diis(:,:,:) = 0d0
  eGW(:,:)          = eHF(:,:)
  eOld(:,:)         = eGW(:,:)
  Z(:,:)            = 1d0
  rcond(:)          = 0d0

!------------------------------------------------------------------------
! Main loop
!------------------------------------------------------------------------

  do while(Conv > thresh .and. nSCF <= maxSCF)

    ! Compute screening

    call phULR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,1d0,eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
    if(.not.TDA) call phULR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nSa,nSb,nSt,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)
    
    call phULR(TDA_W,nSa,nSb,nSt,Aph,Bph,EcRPA,Om,XpY,XmY)

    !----------------------!
    ! Excitation densities !
    !----------------------!

    call UGW_excitation_density(nBas,nC,nO,nR,nSa,nSb,nSt,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY,rho)

    !------------------------------------------------!
    ! Compute self-energy and renormalization factor !
    !------------------------------------------------!

    if(regularize) then
      do is=1,nspin
        call GW_regularization(nBas,nC(is),nO(is),nV(is),nR(is),nSt,eGW(:,is),Om,rho(:,:,:,is))
      end do
    end if

    call UGW_self_energy_diag(eta,nBas,nC,nO,nV,nR,nSt,eGW,Om,rho,SigC,Z,EcGM)

    !-----------------------------------!
    ! Solve the quasi-particle equation !
    !-----------------------------------!

    if(linearize) then

      write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
      write(*,*)

      eGW(:,:) = eHF(:,:) + SigC(:,:)

    else

      write(*,*) ' *** Quasiparticle energies obtained by root search (experimental) *** '
      write(*,*)

      do is=1,nspin

        write(*,*)'-----------------------------------------------------'
        if(is==1) write(*,*)'    Spin-up   orbitals    '
        if(is==2) write(*,*)'    Spin-down orbitals    '

        call UGW_QP_graph(eta,nBas,nC(is),nO(is),nV(is),nR(is),nSt,eHF(:,is), &
                          Om,rho(:,:,:,is),eOld(:,is),eOld(:,is),eGW(:,is),Z(:,is))
      end do

    end if

    ! Convergence criteria

    Conv = maxval(abs(eGW(:,:) - eOld(:,:)))

    ! Print results

    call print_evUGW(nBas,nO,nSCF,Conv,eHF,ENuc,EUHF,SigC,Z,eGW,EcRPA,EcGM)

    ! Linear mixing or DIIS extrapolation

    if(linear_mixing) then
 
      eGW(:,:) = alpha*eGW(:,:) + (1d0 - alpha)*eOld(:,:)
 
    else

      n_diis = min(n_diis+1,max_diis)
      do is=1,nspin
        call DIIS_extrapolation(rcond(ispin),nBas,nBas,n_diis,error_diis(:,1:n_diis,is), & 
                                e_diis(:,1:n_diis,is),eGW(:,is)-eOld(:,is),eGW(:,is))
      end do

!    Reset DIIS if required

      if(minval(rcond(:)) < 1d-15) n_diis = 0

    end if

    ! Save quasiparticles energy for next cycle

    eOld(:,:) = eGW(:,:)

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

  end if

! Deallocate memory

  deallocate(eOld,Z,SigC,Om,XpY,XmY,rho,error_diis,e_diis)

! Perform BSE calculation

  if(BSE) then

    call UGW_phBSE(TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS, &
                   S,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,cHF,eGW,eGW,EcBSE)

    if(exchange_kernel) then

      EcBSE(1) = 0.5d0*EcBSE(1)
      EcBSE(2) = 0.5d0*EcBSE(2)

    else

      EcBSE(2) = 0.0d0

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evUGW correlation energy (spin-conserved) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evUGW correlation energy (spin-flip)      =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evUGW correlation energy                  =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@evUGW total energy                        =',ENuc + EUHF + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '--------------------------------------------------------------'
      write(*,*) ' Adiabatic connection version of BSE@evUGW correlation energy '
      write(*,*) '--------------------------------------------------------------'
      write(*,*)

      if(doXBS) then

        write(*,*) '*** scaled screening version (XBS) ***'
        write(*,*)

      end if

      call UGW_phACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,spin_conserved,spin_flip, &
                       eta,nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,eGW,eGW,EcAC)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evUGW correlation energy (spin-conserved) =',EcAC(1)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evUGW correlation energy (spin-flip)      =',EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evUGW correlation energy                  =',EcAC(1) + EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@evUGW total energy                        =',ENuc + EUHF + EcAC(1) + EcAC(2)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

! Testing zone
  
  if(dotest) then
  
    call dump_test_value('U','evUGW correlation energy',EcRPA)
    call dump_test_value('U','evUGW HOMOa energy',eGW(nO(1),1))
    call dump_test_value('U','evUGW LUMOa energy',eGW(nO(1)+1,1))
    call dump_test_value('U','evUGW HOMOa energy',eGW(nO(2),2))
    call dump_test_value('U','evUGW LUMOa energy',eGW(nO(2)+1,2))

  end if

end subroutine 
