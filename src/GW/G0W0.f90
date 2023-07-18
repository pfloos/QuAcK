subroutine G0W0(doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,evDyn,doppBSE, & 
                singlet,triplet,linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_AO,ERI_MO,    & 
                dipole_int,PHF,cHF,eHF,Vxc)

! Perform G0W0 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: evDyn
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
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI_AO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_MO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: Vxc(nBas)
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: PHF(nBas,nBas)

! Local variables

  logical                       :: print_W = .true.
  logical                       :: dRPA
  integer                       :: ispin
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: SigX(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)

  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: eGWlin(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0W0 calculation           |'
  write(*,*)'************************************************'
  write(*,*)

! Initialization

  dRPA = .true.
  EcRPA = 0d0

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

! Spin manifold 

  ispin = 1

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),SigC(nBas),SigX(nBas),Z(nBas),Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nBas,nBas,nS), & 
           eGW(nBas),eGWlin(nBas))

!-------------------!
! Compute screening !
!-------------------!

  call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI_MO,Aph)
  if(.not.TDA_W) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)

  call phLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

  if(print_W) call print_excitation('RPA@HF      ',ispin,nS,Om)

!--------------------------!
! Compute spectral weights !
!--------------------------!

  call GW_excitation_density(nBas,nC,nO,nR,nS,ERI_MO,XpY,rho)

!------------------------!
! Compute GW self-energy !
!------------------------!

  call self_energy_exchange_diag(nBas,cHF,PHF,ERI_AO,SigX)

  if(regularize) then 

    call regularized_self_energy_correlation_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigC)
    call regularized_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,Z)

  else

    call GW_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigC,Z)

  end if

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  eGWlin(:) = eHF(:) + Z(:)*(SigX(:) + SigC(:) - Vxc(:))

  ! Linearized or graphical solution?

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGW(:) = eGWlin(:)

  else 

    write(*,*) ' *** Quasiparticle energies obtained by root search (experimental) *** '
    write(*,*)
  
    call QP_graph(nBas,nC,nO,nV,nR,nS,eta,eHF,SigX,Vxc,Om,rho,eGWlin,eGW,regularize)

    ! Find all the roots of the QP equation if necessary

    ! call QP_roots(nBas,nC,nO,nV,nR,nS,eta,eHF,Om,rho,eGWlin) 
 
  end if

! Compute the RPA correlation energy

  call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI_MO,Aph)
  if(.not.TDA_W) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)

  call phLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

!--------------!
! Dump results !
!--------------!

  call print_G0W0(nBas,nO,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA,EcGM)

! Deallocate memory

  deallocate(SigC,Z,Om,XpY,XmY,rho,eGWlin)

! Perform BSE calculation

  if(dophBSE) then

    call GW_phBSE(dophBSE2,TDA_W,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eHF,eGW,EcBSE)

    if(exchange_kernel) then
 
      EcBSE(1) = 0.5d0*EcBSE(1)
      EcBSE(2) = 1.5d0*EcBSE(2)
 
    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0 correlation energy (singlet) =',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0 correlation energy (triplet) =',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0 correlation energy           =',EcBSE(1) + EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0 total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '-------------------------------------------------------------'
      write(*,*) ' Adiabatic connection version of BSE@G0W0 correlation energy '
      write(*,*) '-------------------------------------------------------------'
      write(*,*) 

      if(doXBS) then 

        write(*,*) '*** scaled screening version (XBS) ***'
        write(*,*)

      end if

      call GW_phACFDT(exchange_kernel,doXBS,dRPA,TDA_W,TDA,dophBSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,eHF,eGW,EcAC)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0W0 correlation energy (singlet) =',EcAC(1),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0W0 correlation energy (triplet) =',EcAC(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0W0 correlation energy           =',EcAC(1) + EcAC(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0W0 total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2),' au'
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

  if(doppBSE) then

    call GW_ppBSE(TDA_W,TDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eHF,eGW,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 correlation energy (singlet) =',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 correlation energy (triplet) =',3d0*EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 correlation energy           =',EcBSE(1) + 3d0*EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 total energy                 =',ENuc + ERHF + EcBSE(1) + 3d0*EcBSE(2),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine 
