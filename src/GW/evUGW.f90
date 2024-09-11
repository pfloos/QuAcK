subroutine evUGW(dotest,maxSCF,thresh,max_diis,doACFDT,exchange_kernel,doXBS,BSE,TDA_W,TDA,dBSE,dTDA, & 
                 spin_conserved,spin_flip,linearize,eta,doSRG,nBas,nC,nO,nV,nR,nS,ENuc,   &
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
  logical,intent(in)            :: doSRG

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
  integer                       :: is
  integer                       :: ispin
  integer                       :: nSCF
  integer                       :: n_diis
  double precision              :: flow
  double precision              :: rcond(nspin)
  double precision              :: Conv
  double precision              :: EcRPA(nspin)
  double precision              :: EcGM(nspin)
  double precision              :: EcBSE(nspin)
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

! SRG regularization

  flow = 500d0

  if(doSRG) then

    write(*,*) '*** SRG regularized evGW scheme ***'
    write(*,*)

  end if

! Initialization

  EcRPA(:) = 0d0
  dRPA = .true.

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
    
    call phULR(TDA_W,nSa,nSb,nSt,Aph,Bph,EcRPA(ispin),Om,XpY,XmY)

    !----------------------!
    ! Excitation densities !
    !----------------------!

    call UGW_excitation_density(nBas,nC,nO,nR,nSa,nSb,nSt,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY,rho)

    !------------------------------------------------!
    ! Compute self-energy and renormalization factor !
    !------------------------------------------------!

    if(doSRG) then
      call UGW_SRG_self_energy_diag(flow,nBas,nC,nO,nV,nR,nSt,eGW,Om,rho,SigC,Z,EcGM)
    else
      call UGW_self_energy_diag(eta,nBas,nC,nO,nV,nR,nSt,eGW,Om,rho,SigC,Z,EcGM)
    end if


    !-----------------------------------!
    ! Solve the quasi-particle equation !
    !-----------------------------------!

    if(linearize) then

      write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
      write(*,*)

      eGW(:,:) = eHF(:,:) + SigC(:,:)

    else

      write(*,*) ' *** Quasiparticle energies obtained by root search *** '
      write(*,*)

      do is=1,nspin

        write(*,*)'-----------------------------------------------------'
        if(is==1) write(*,*)'    Spin-up   orbitals    '
        if(is==2) write(*,*)'    Spin-down orbitals    '

        call UGW_QP_graph(doSRG,eta,flow,nBas,nC(is),nO(is),nV(is),nR(is),nSt,eHF(:,is), &
                          Om,rho(:,:,:,is),eOld(:,is),eOld(:,is),eGW(:,is),Z(:,is))
      end do

    end if

    ! Convergence criteria

    Conv = maxval(abs(eGW(:,:) - eOld(:,:)))

    ! Print results

    call print_evUGW(nBas,nO,nSCF,Conv,eHF,ENuc,EUHF,SigC,Z,eGW,EcRPA(ispin),EcGM)

    ! Linear mixing or DIIS extrapolation

    if(max_diis > 1) then

      n_diis = min(n_diis+1,max_diis)
      do is=1,nspin
        call DIIS_extrapolation(rcond(ispin),nBas,nBas,n_diis,error_diis(:,1:n_diis,is), & 
                                e_diis(:,1:n_diis,is),eGW(:,is)-eOld(:,is),eGW(:,is))
      end do

    end if

    ! Save quasiparticles energy for next cycle

    eOld(:,:) = eGW(:,:)

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

  deallocate(eOld,Z,SigC,Om,XpY,XmY,rho,error_diis,e_diis)

! Perform BSE calculation

  if(BSE) then

    call UGW_phBSE(exchange_kernel,TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS, &
                   S,ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,cHF,eGW,eGW,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@UHF correlation energy (spin-conserved) =',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@UHF correlation energy (spin-flip)      =',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@UHF correlation energy                  =',sum(EcBSE),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@evGW@UHF total energy                        =',ENuc + EUHF + sum(EcBSE),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      call UGW_phACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,spin_conserved,spin_flip, &
                       eta,nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,eGW,eGW,EcRPA)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@UHF correlation energy (spin-conserved) =',EcRPA(1),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@UHF correlation energy (spin-flip)      =',EcRPA(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@UHF correlation energy                  =',sum(EcRPA),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@evGW@UHF total energy                        =',ENuc + EUHF + sum(EcRPA),' au'
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

! Testing zone
  
  if(dotest) then
  
    call dump_test_value('U','evGW correlation energy',EcRPA)
    call dump_test_value('U','evGW HOMOa energy',eGW(nO(1),1))
    call dump_test_value('U','evGW LUMOa energy',eGW(nO(1)+1,1))
    call dump_test_value('U','evGW HOMOa energy',eGW(nO(2),2))
    call dump_test_value('U','evGW LUMOa energy',eGW(nO(2)+1,2))

  end if

end subroutine 
