subroutine UG0W0(doACFDT,exchange_kernel,doXBS,BSE,TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip,      &
                 linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EUHF,S,ERI_aaaa,ERI_aabb,ERI_bbbb, & 
                 dipole_int_aa,dipole_int_bb,cHF,eHF)

! Perform unrestricted G0W0 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

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
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: S(nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)

! Local variables

  logical                       :: print_W = .true.
  logical                       :: dRPA
  integer                       :: is
  integer                       :: ispin
  double precision              :: EcRPA
  double precision              :: EcGM(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: Z(:,:)
  integer                       :: nS_aa,nS_bb,nS_sc
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:,:)

  double precision,allocatable  :: eGWlin(:,:)

! Output variables

  double precision              :: eGW(nBas,nspin)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0W0 calculation           |'
  write(*,*)'|         *** Unrestricted version ***         |'
  write(*,*)'************************************************'
  write(*,*)

! Initialization

  EcRPA = 0d0
  dRPA = .true.

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

! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  allocate(SigC(nBas,nspin),Z(nBas,nspin),eGWlin(nBas,nspin), &
           Om(nS_sc),XpY(nS_sc,nS_sc),XmY(nS_sc,nS_sc),rho(nBas,nBas,nS_sc,nspin))

!-------------------!
! Compute screening !
!-------------------!

! Spin-conserving transitions

  ispin = 1

  call phULR(ispin,dRPA,TDA_W,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
             eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,Om,rho,EcRPA,Om,XpY,XmY)

  if(print_W) call print_excitation_energies('phRPA@UHF',5,nS_sc,Om)

!----------------------!
! Excitation densities !
!----------------------!
 
  call UGW_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY,rho)

!------------------------------------------------!
! Compute self-energy and renormalization factor !
!------------------------------------------------!

  if(regularize) then
    do is=1,nspin
      call GW_regularization(nBas,nC(is),nO(is),nV(is),nR(is),nS_sc,eHF(:,is),Om,rho(:,:,:,is))
    end do
  end if

  call UGW_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,Om,rho,SigC,Z,EcGM)

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  eGWlin(:,:) = eHF(:,:) + Z(:,:)*SigC(:,:) 

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGW(:,:) = eGWlin(:,:)

  else 
  
  ! Find graphical solution of the QP equation

    do is=1,nspin
      call UGW_QP_graph(eta,nBas,nC(is),nO(is),nV(is),nR(is),nS_sc,eHF(:,is), & 
                        Om,rho(:,:,:,is),eHF(:,is),eGW(:,is),Z(:,is))
    end do
 
  end if

! Compute RPA correlation energy

  call phULR(ispin,dRPA,TDA_W,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
             eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,Om,rho,EcRPA,Om,XpY,XmY)

! Dump results

  call print_UG0W0(nBas,nO,eHF,ENuc,EUHF,SigC,Z,eGW,EcRPA,EcGM)

! Free memory

  deallocate(Om,XpY,XmY,rho)

! Perform BSE calculation

  if(BSE) then

    call UGW_phBSE(TDA_W,TDA,dBSE,dTDA,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS,S, &
                   ERI_aaaa,ERI_aabb,ERI_bbbb,dipole_int_aa,dipole_int_bb,cHF,eHF,eGW,EcBSE)

    if(exchange_kernel) then
 
      EcBSE(1) = 0.5d0*EcBSE(1)
      EcBSE(2) = 0.5d0*EcBSE(2)
 
    else
 
      EcBSE(2) = 0.0d0

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@UG0W0 correlation energy (spin-conserved) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@UG0W0 correlation energy (spin-flip)      =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@UG0W0 correlation energy                  =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@UG0W0 total energy                        =',ENuc + EUHF + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '------------------------------------------------------------'
      write(*,*) 'Adiabatic connection version of BSE@UG0W0 correlation energy'
      write(*,*) '------------------------------------------------------------'
      write(*,*) 

      if(doXBS) then 

        write(*,*) '*** scaled screening version (XBS) ***'
        write(*,*)

      end if

      call UGW_phACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,spin_conserved,spin_flip,eta, & 
                       nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,eHF,eGW,EcAC)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10)') 'AC@BSE@UG0W0 correlation energy (spin-conserved) =',EcAC(1)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@UG0W0 correlation energy (spin-flip)      =',EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@UG0W0 correlation energy                  =',EcAC(1) + EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@UG0W0 total energy                        =',ENuc + EUHF + EcAC(1) + EcAC(2)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

end subroutine 
