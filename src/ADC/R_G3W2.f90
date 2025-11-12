subroutine R_G3W2(dotest,TDA_W,singlet,triplet,linearize,eta,doSRG,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform single-shot G3W2 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: TDA_W
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
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)
  double precision,intent(in)   :: eHF(nOrb)

! Local variables

  logical                       :: print_W   = .false.
  logical                       :: plot_self = .false.
  logical                       :: dRPA_W
  integer                       :: isp_W
  double precision              :: flow
  double precision              :: EcRPA
  double precision              :: EcGM
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)

  double precision,allocatable  :: eQPlin(:)
  double precision,allocatable  :: eQP(:)

! Output variables

  ! None

! Hello world

  write(*,*)
  write(*,*)'*******************************'
  write(*,*)'* Restricted G3W2 Calculation *'
  write(*,*)'*******************************'
  write(*,*)

! Spin manifold and TDA for dynamical screening

  isp_W = 1
  dRPA_W = .true.

! SRG regularization

  flow = 500d0

  if(doSRG) then

    write(*,*) '*** SRG regularized G3W2 scheme ***'
    write(*,*)

  end if

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),SigC(nOrb),Z(nOrb),Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nOrb,nOrb,nS), & 
           eQP(nOrb),eQPlin(nOrb))

!-------------------!
! Compute screening !
!-------------------!

                 call phRLR_A(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  if(.not.TDA_W) call phRLR_B(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

  if(print_W) call print_excitation_energies('phRPA@RHF','singlet',nS,Om)

!--------------------------!
! Compute spectral weights !
!--------------------------!

  call RGW_excitation_density(nOrb,nC,nO,nR,nS,ERI,XpY,rho)

!------------------------!
! Compute GW self-energy !
!------------------------!

  if(doSRG) then 

!   call R_G3W2_SRG_self_energy_diag(flow,nBas,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigC,Z)

  else

    call R_G3W2_self_energy_diag(eta,nBas,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,ERI,EcGM,SigC,Z)

  end if
  
!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  ! Linearized or graphical solution?

  eQPlin(:) = eHF(:) + Z(:)*SigC(:)
  
  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eQP(:) = eQPlin(:)

  else 

     write(*,*) ' *** Quasiparticle energies obtained by root search *** '
     write(*,*)

     call R_G3W2_QP_graph(doSRG,eta,flow,nBas,nOrb,nC,nO,nV,nR,nS,eHF,Om,rho,ERI,eQPlin,eHF,eQP,Z)


  end if

! Plot self-energy, renormalization factor, and spectral function

! if(plot_self) call RGW_plot_self_energy(nOrb,eta,nC,nO,nV,nR,nS,eHF,eQP,Om,rho)
  
! Compute the RPA correlation energy

                 call phRLR_A(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,eQP,ERI,Aph)
  if(.not.TDA_W) call phRLR_B(isp_W,dRPA_W,nOrb,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phRLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

!--------------!
! Dump results !
!--------------!

  call print_R_G3W2(nOrb,nC,nO,nV,nR,eHF,ENuc,ERHF,SigC,Z,eQP,EcRPA,EcGM)
  
! Testing zone

  if(dotest) then

    call dump_test_value('R','G3W2 correlation energy',EcRPA)
    call dump_test_value('R','G3W2 HOMO energy',eQP(nO))
    call dump_test_value('R','G3W2 LUMO energy',eQP(nO+1))

  end if

end subroutine 
