subroutine G0W0(doACFDT,exchange_kernel,doXBS,COHSEX,BSE,BSE2,TDA_W,TDA,dBSE,dTDA,evDyn,ppBSE, & 
                singlet,triplet,linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,   & 
                ERI_AO,ERI_MO,dipole_int,PHF,cHF,eHF,Vxc)

! Perform G0W0 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: BSE
  logical,intent(in)            :: BSE2
  logical,intent(in)            :: ppBSE
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
  integer                       :: ispin
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcppBSE(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: SigX(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

  double precision,allocatable  :: eGW(:)
  double precision,allocatable  :: eGWlin(:)

  integer                       :: nBas2
  integer                       :: nC2
  integer                       :: nO2
  integer                       :: nV2
  integer                       :: nR2
  integer                       :: nS2

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0W0 calculation           |'
  write(*,*)'************************************************'
  write(*,*)

! Initialization

  EcRPA = 0d0

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

! Spin manifold 

  ispin = 1

! Memory allocation

  allocate(SigC(nBas),SigX(nBas),Z(nBas),OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS),rho_RPA(nBas,nBas,nS),eGW(nBas),eGWlin(nBas))

!-------------------!
! Compute screening !
!-------------------!

  call phLR(ispin,.true.,TDA_W,eta,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI_MO,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

  if(print_W) call print_excitation('RPA@HF      ',ispin,nS,OmRPA)

!--------------------------!
! Compute spectral weights !
!--------------------------!

  call GW_excitation_density(nBas,nC,nO,nR,nS,ERI_MO,XpY_RPA,rho_RPA)

!------------------------!
! Compute GW self-energy !
!------------------------!

  call self_energy_exchange_diag(nBas,cHF,PHF,ERI_AO,SigX)

  if(regularize) then 

    call regularized_self_energy_correlation_diag(COHSEX,eta,nBas,nC,nO,nV,nR,nS,eHF,OmRPA,rho_RPA,EcGM,SigC)
    call regularized_renormalization_factor(COHSEX,eta,nBas,nC,nO,nV,nR,nS,eHF,OmRPA,rho_RPA,Z)

  else

    call GW_self_energy_diag(COHSEX,eta,nBas,nC,nO,nV,nR,nS,eHF,OmRPA,rho_RPA,EcGM,SigC,Z)

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
  
    call QP_graph(nBas,nC,nO,nV,nR,nS,eta,eHF,SigX,Vxc,OmRPA,rho_RPA,eGWlin,eGW,regularize)

    ! Find all the roots of the QP equation if necessary

    ! call QP_roots(nBas,nC,nO,nV,nR,nS,eta,eHF,OmRPA,rho_RPA,eGWlin) 
 
  end if

! Compute the RPA correlation energy

  call phLR(ispin,.true.,TDA_W,eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI_MO,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

!--------------!
! Dump results !
!--------------!

  call print_G0W0(nBas,nO,eHF,ENuc,ERHF,SigC,Z,eGW,EcRPA,EcGM)

! Deallocate memory

! deallocate(SigC,Z,OmRPA,XpY_RPA,XmY_RPA,rho_RPA,eGWlin)

! Plot stuff

!  call plot_GW(nBas,nC,nO,nV,nR,nS,eHF,eGW,OmRPA,rho_RPA)

! Perform BSE calculation

  if(BSE) then

    call GW_phBSE(BSE2,TDA_W,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eHF,eGW,EcBSE)

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

      call GW_phACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,eHF,eGW,EcAC)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@G0W0 correlation energy (singlet) =',EcAC(1),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@G0W0 correlation energy (triplet) =',EcAC(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@G0W0 correlation energy           =',EcAC(1) + EcAC(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@G0W0 total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2),' au'
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

  if(ppBSE) then

    call GW_ppBSE(TDA_W,TDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eHF,eGW,EcppBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 correlation energy (singlet) =',EcppBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 correlation energy (triplet) =',3d0*EcppBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 correlation energy           =',EcppBSE(1) + 3d0*EcppBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 total energy                 =',ENuc + ERHF + EcppBSE(1) + 3d0*EcppBSE(2),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

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
!   call  GW_ppBSE_so(TDA_W,TDA,singlet,triplet,eta,nBas2,nC2,nO2,nV2,nR2,nS2,sERI,dipole_int,seHF,seGW,EcppBSE)

  end if

! if(BSE) call ufBSE(nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF,eGW)
! if(BSE) call ufXBSE(nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF,OmRPA,rho_RPA)

  if(BSE) call XBSE(TDA_W,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eHF,eGW,EcBSE)


end subroutine 
