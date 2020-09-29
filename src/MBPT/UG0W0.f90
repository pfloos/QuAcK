subroutine UG0W0(doACFDT,exchange_kernel,doXBS,COHSEX,BSE,TDA_W,TDA,dBSE,dTDA,evDyn,spin_conserved,spin_flip, &
                 linearize,eta,nBas,nC,nO,nV,nR,nS,ENuc,EUHF,Hc,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,          & 
                 dipole_int_aa,dipole_int_bb,PHF,cHF,eHF,eGW)

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
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip
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
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_abab(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)

! Local variables

  logical                       :: print_W = .true.
  integer                       :: is
  integer                       :: ispin
  double precision              :: EcRPA
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision,allocatable  :: SigC(:,:)
  double precision,allocatable  :: Z(:,:)
  integer                       :: nS_aa,nS_bb,nS_sc
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:,:)

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

! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  allocate(SigC(nBas,nspin),Z(nBas,nspin),OmRPA(nS_sc),XpY_RPA(nS_sc,nS_sc),XmY_RPA(nS_sc,nS_sc), &
           rho_RPA(nBas,nBas,nS_sc,nspin),eGWlin(nBas,nspin))

!-------------------!
! Compute screening !
!-------------------!

! Spin-conserving transition

  ispin = 1

  call unrestricted_linear_response(ispin,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
                                    eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

  if(print_W) call print_excitation('RPA@UHF     ',5,nS_sc,OmRPA)

!----------------------!
! Excitation densities !
!----------------------!
 
  call unrestricted_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

!---------------------!
! Compute self-energy !
!---------------------!

  call unrestricted_self_energy_correlation_diag(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,OmRPA,rho_RPA,SigC)

!--------------------------------!
! Compute renormalization factor !
!--------------------------------!

  call unrestricted_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,OmRPA,rho_RPA,Z)

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
    call unrestricted_QP_graph(nBas,nC(is),nO(is),nV(is),nR(is),nS_sc,eta,eHF(:,is),OmRPA, & 
                               rho_RPA(:,:,:,is),eGWlin(:,is),eGW(:,is))
  end do
 
  end if

! Compute RPA correlation energy

  call unrestricted_linear_response(ispin,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
                                    eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

! Dump results

  call print_UG0W0(nBas,nO,eHF,ENuc,EUHF,SigC,Z,eGW,EcRPA)

! Free memory

  deallocate(OmRPA,XpY_RPA,XmY_RPA,rho_RPA)

! Perform BSE calculation

  if(BSE) then

    call unrestricted_Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS, &
                                     ERI_aaaa,ERI_aabb,ERI_bbbb,ERI_abab,dipole_int_aa,dipole_int_bb,eHF,eGW,EcBSE)

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

  end if

end subroutine UG0W0