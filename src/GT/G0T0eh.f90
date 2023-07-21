subroutine G0T0eh(doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_T,TDA,dBSE,dTDA,doppBSE, & 
                  singlet,triplet,linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,   & 
                  ERI_AO,ERI_MO,dipole_int,PHF,cHF,eHF,Vxc)

! Perform ehG0T0 calculation

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
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
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
  logical                       :: dRPA = .false.
  integer                       :: ispin
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: SigX(:)
  double precision,allocatable  :: Sig(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rhoL(:,:,:,:)
  double precision,allocatable  :: rhoR(:,:,:,:)

  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eGTlin(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|        One-shot eh G0T0 calculation          |'
  write(*,*)'************************************************'
  write(*,*)

! Initialization

  EcRPA = 0d0

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

! Spin manifold 

  ispin = 2

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),Sig(nBas),SigX(nBas),Z(nBas),Om(nS),XpY(nS,nS),XmY(nS,nS), & 
           rhoL(nBas,nBas,nS,2),rhoR(nBas,nBas,nS,2),eGT(nBas),eGTlin(nBas))

!-------------------!
! Compute screening !
!-------------------!

  call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI_MO,Aph)
  if(.not.TDA_T) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)

  call phLR(TDA_T,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

  if(print_W) call print_excitation('RPA@HF      ',ispin,nS,Om)

!--------------------------!
! Compute spectral weights !
!--------------------------!

  call GTeh_excitation_density(nBas,nC,nO,nR,nS,ERI_MO,XpY,XmY,rhoL,rhoR)

!------------------------!
! Compute GW self-energy !
!------------------------!

  call self_energy_exchange_diag(nBas,cHF,PHF,ERI_AO,SigX)

  if(regularize) then 

!   call regularized_self_energy_correlation_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,Sig)
!   call regularized_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,Z)

  else

    call GTeh_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rhoL,rhoR,EcGM,Sig,Z)

  end if

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  eGTlin(:) = eHF(:) + Z(:)*(SigX(:) + Sig(:) - Vxc(:))

  ! Linearized or graphical solution?

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGT(:) = eGTlin(:)

  else 

    write(*,*) ' *** Quasiparticle energies obtained by root search (experimental) *** '
    write(*,*)
  
!   call QP_graph(nBas,nC,nO,nV,nR,nS,eta,eHF,SigX,Vxc,Om,rho,eGTlin,eGT)

    ! Find all the roots of the QP equation if necessary

    ! call QP_roots(nBas,nC,nO,nV,nR,nS,eta,eHF,Om,rho,eGTlin)
 
  end if

! Compute the RPA correlation energy

  call phLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGT,ERI_MO,Aph)
  if(.not.TDA_T) call phLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI_MO,Bph)

  call phLR(TDA_T,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

!--------------!
! Dump results !
!--------------!

  call print_G0T0eh(nBas,nO,eHF,ENuc,ERHF,Sig,Z,eGT,EcRPA,EcGM)

! Deallocate memory

! deallocate(Sig,Z,Om,XpY,XmY,rho,eGTlin)

! Perform BSE calculation

! if(BSE) then

!   call Bethe_Salpeter(BSE2,TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eHF,eGW,EcBSE)

!   if(exchange_kernel) then
!
!     EcBSE(1) = 0.5d0*EcBSE(1)
!     EcBSE(2) = 1.5d0*EcBSE(2)
!
!   end if

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0 correlation energy (singlet) =',EcBSE(1),' au'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0 correlation energy (triplet) =',EcBSE(2),' au'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0 correlation energy           =',EcBSE(1) + EcBSE(2),' au'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0 total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2),' au'
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

!   if(doACFDT) then

!     write(*,*) '-------------------------------------------------------------'
!     write(*,*) ' Adiabatic connection version of BSE@G0W0 correlation energy '
!     write(*,*) '-------------------------------------------------------------'
!     write(*,*) 

!     if(doXBS) then 

!       write(*,*) '*** scaled screening version (XBS) ***'
!       write(*,*)

!     end if

!     call ACFDT(exchange_kernel,doXBS,.true.,TDA_T,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,eHF,eGW,EcAC)

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@G0W0 correlation energy (singlet) =',EcAC(1),' au'
!     write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@G0W0 correlation energy (triplet) =',EcAC(2),' au'
!     write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@G0W0 correlation energy           =',EcAC(1) + EcAC(2),' au'
!     write(*,'(2X,A50,F20.10,A3)') 'AC@BSE@G0W0 total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2),' au'
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

! end if

! if(ppBSE) then

!   call Bethe_Salpeter_pp(TDA_T,TDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eHF,eGW,EcBSE)

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 correlation energy (singlet) =',EcBSE(1),' au'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 correlation energy (triplet) =',3d0*EcBSE(2),' au'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 correlation energy           =',EcBSE(1) + 3d0*EcBSE(2),' au'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0 total energy                 =',ENuc + ERHF + EcBSE(1) + 3d0*EcBSE(2),' au'
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

! end if

! if(BSE) call ufBSE(nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF,eGW)
! if(BSE) call ufXBSE(nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,eHF,Om,rho)

! if(BSE) call XBSE(TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eHF,eGW,EcBSE)


end subroutine 
