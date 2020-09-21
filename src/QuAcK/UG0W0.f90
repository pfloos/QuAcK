subroutine UG0W0(doACFDT,exchange_kernel,doXBS,COHSEX,BSE,TDA_W,TDA,dBSE,dTDA,evDyn,  &
                 singlet_manifold,triplet_manifold,linearize,eta,nBas,nC,nO,nV,nR,nS, &
                 ENuc,EUHF,Hc,ERI_aa,ERI_ab,ERI_bb,PHF,cHF,eHF,eGW)

! Perform unrestricted G0W0 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

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
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EUHF
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: PHF(nBas,nBas,nspin)
  double precision,intent(in)   :: Hc(nBas,nBas,nspin)
  double precision,intent(in)   :: ERI_aa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_ab(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bb(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: print_W = .true.
  integer                       :: ispin
  integer                       :: iblock
  integer                       :: bra
  integer                       :: ket
  integer                       :: nSa
  integer                       :: nSb
  integer                       :: nSt
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: Z(:,:)
  double precision,allocatable  :: Omega(:)
  double precision,allocatable  :: XpY_a(:,:)
  double precision,allocatable  :: XpY_b(:,:)
  double precision,allocatable  :: XmY_a(:,:)
  double precision,allocatable  :: XmY_b(:,:)
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

  EcRPA(:) = 0d0

! COHSEX approximation

  if(COHSEX) write(*,*) 'COHSEX approximation activated!'
  write(*,*)

! TDA for W

  if(TDA_W) write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
  write(*,*)

! TDA 

  if(TDA) write(*,*) 'Tamm-Dancoff approximation activated!'
  write(*,*)

! Memory allocation

  nSa = nS(1)
  nSb = nS(2)
  nSt = nSa + nSb

  allocate(SigC(nBas,nspin),Z(nBas,nspin),Omega(nSt),XpY_a(nSa,nSa),XpY_b(nSb,nSb),XmY_a(nSa,nSa),XmY_b(nSb,nSb), & 
           rho(nBas,nBas,nSt,nspin),eGWlin(nBas,nspin))

! Compute linear response

!----------------------------------------------
! alpha-alpha block
!----------------------------------------------

  ispin  = 1
  iblock = 3

  call linear_response(iblock,.true.,TDA_W,.false.,eta,nBas,nC(ispin),nO(ispin),nV(ispin),nR(ispin),nSa,1d0, &
                       eHF(:,ispin),ERI_aa,rho(:,:,1:nSa,ispin),EcRPA(ispin),Omega(1:nSa),XpY_a,XmY_a)

  if(print_W) call print_excitation('RPA@HF alpha',iblock,nSa,Omega(1:nSa))

!----------------------------------------------
! alpha-beta block
!----------------------------------------------

  ispin  = 2
  iblock = 3

  call linear_response(iblock,.true.,TDA_W,.false.,eta,nBas,nC(ispin),nO(ispin),nV(ispin),nR(ispin),nSb,1d0, &
                       eHF(:,ispin),ERI_bb,rho(:,:,nSa+1:nSt,ispin),EcRPA(ispin),Omega(nSa+1:nSt),XpY_b,XmY_b)

  if(print_W) call print_excitation('RPA@HF beta ',iblock,nSb,Omega(nSa+1:nSt))

!----------------------------------------------
! Excitation densities for alpha self-energy 
!----------------------------------------------
 
  call unrestricted_excitation_density(nBas,nC,nO,nR,nSa,nSb,nSt,ERI_aa,ERI_ab,ERI_bb,XpY_a,XpY_b,rho)

!----------------------
! Compute self-energy
!----------------------

  call unrestricted_self_energy_correlation_diag(eta,nBas,nC,nO,nV,nR,nSa,nSb,nSt,eHF,Omega,rho,SigC)

! Compute renormalization factor

! call renormalization_factor(COHSEX,SOSEX,eta,nBas,nC,nO,nV,nR,nS,eHF, & 
!                             Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),Z(:))

! Solve the quasi-particle equation
  Z(:,:) = 1d0
  eGWlin(:,:) = eHF(:,:) + Z(:,:)*SigC(:,:)

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGW(:,:) = eGWlin(:,:)

  else 
  
  ! Find graphical solution of the QP equation

! do is=1,nspin
!   call QP_graph(nBas,nC(:,is),nO(:,is),nV(:,is),nR(:,is),nS(:,is),eta,eHF(:,is),Omega(:,is), & 
!                 rho(:,:,:,ispin),eGWlin(:,is),eGW(:,is))
! end do
 
  end if

! Dump results

  do ispin=1,nspin
    call print_G0W0(nBas,nO(ispin),eHF(:,ispin),ENuc,EUHF,SigC(:,ispin),Z(:,ispin),eGW(:,ispin),EcRPA(ispin))
  end do

! Compute the RPA correlation energy

! call linear_response(ispin,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI, & 
!                      rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A50,F20.10)') 'Tr@RPA@G0W0 correlation energy (singlet) =',EcRPA(1)
  write(*,'(2X,A50,F20.10)') 'Tr@RPA@G0W0 correlation energy (triplet) =',EcRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@RPA@G0W0 correlation energy           =',EcRPA(1) + EcRPA(2)
  write(*,'(2X,A50,F20.10)') 'Tr@RPA@G0W0 total energy                 =',ENuc + EUHF + EcRPA(1) + EcRPA(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

! Perform BSE calculation

! if(BSE) then

!   call Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,singlet_manifold,triplet_manifold,eta, &
!                       nBas,nC,nO,nV,nR,nS,ERI,eHF,eGW,Omega,XpY,XmY,rho,EcRPA,EcBSE)

!   if(exchange_kernel) then
!
!     EcRPA(1) = 0.5d0*EcRPA(1)
!     EcRPA(2) = 1.5d0*EcRPA(1)
!
!   end if

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@G0W0 correlation energy (singlet) =',EcBSE(1)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@G0W0 correlation energy (triplet) =',EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@G0W0 correlation energy           =',EcBSE(1) + EcBSE(2)
!   write(*,'(2X,A50,F20.10)') 'Tr@BSE@G0W0 total energy                 =',ENuc + EUHF + EcBSE(1) + EcBSE(2)
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

!     call ACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,singlet_manifold,triplet_manifold,eta, & 
!                nBas,nC,nO,nV,nR,nS,ERI,eHF,eGW,Omega,XpY,XmY,rho,EcAC)

!     if(exchange_kernel) then
! 
!       EcAC(1) = 0.5d0*EcAC(1)
!       EcAC(2) = 1.5d0*EcAC(1)
! 
!     end if

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@G0W0 correlation energy (singlet) =',EcAC(1)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@G0W0 correlation energy (triplet) =',EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@G0W0 correlation energy           =',EcAC(1) + EcAC(2)
!     write(*,'(2X,A50,F20.10)') 'AC@BSE@G0W0 total energy                 =',ENuc + EUHF + EcAC(1) + EcAC(2)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

! end if

end subroutine UG0W0
