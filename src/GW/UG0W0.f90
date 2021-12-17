subroutine UG0W0(doACFDT,exchange_kernel,doXBS,COHSEX,BSE,TDA_W,TDA,dBSE,dTDA,evDyn,spin_conserved,spin_flip, &
                 linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,EUHF,S,ERI,ERI_aaaa,ERI_aabb,ERI_bbbb,     & 
                 dipole_int_aa,dipole_int_bb,PHF,cHF,eHF,Vxc,eGW)

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
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int_aa(nBas,nBas,ncart)
  double precision,intent(in)   :: dipole_int_bb(nBas,nBas,ncart)
  double precision,intent(in)   :: eHF(nBas,nspin)
  double precision,intent(in)   :: cHF(nBas,nBas,nspin)
  double precision,intent(in)   :: PHF(nBas,nBas,nspin)
  double precision,intent(in)   :: Vxc(nBas,nspin)

! Local variables

  logical                       :: print_W = .true.
  integer                       :: is
  integer                       :: ispin
  double precision              :: EcRPA
  double precision              :: EcGM(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision,allocatable  :: SigX(:,:)
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

  allocate(SigX(nBas,nspin),SigC(nBas,nspin),Z(nBas,nspin),eGWlin(nBas,nspin), &
           OmRPA(nS_sc),XpY_RPA(nS_sc,nS_sc),XmY_RPA(nS_sc,nS_sc),rho_RPA(nBas,nBas,nS_sc,nspin))

!-------------------!
! Compute screening !
!-------------------!

! Spin-conserving transition

  ispin = 1

  call unrestricted_linear_response(ispin,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
                                    eHF,ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

  if(print_W) call print_excitation('RPA@UHF     ',5,nS_sc,OmRPA)

!----------------------!
! Excitation densities !
!----------------------!
 
  call unrestricted_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

!---------------------!
! Compute self-energy !
!---------------------!

  do is=1,nspin
    call self_energy_exchange_diag(nBas,cHF(:,:,is),PHF(:,:,is),ERI,SigX(:,is))
  end do

  call unrestricted_self_energy_correlation_diag(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,OmRPA,rho_RPA,SigC,EcGM)

!--------------------------------!
! Compute renormalization factor !
!--------------------------------!

  call unrestricted_renormalization_factor(eta,nBas,nC,nO,nV,nR,nS_sc,eHF,OmRPA,rho_RPA,Z)

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  eGWlin(:,:) = eHF(:,:) + Z(:,:)*(SigX(:,:) + SigC(:,:) - Vxc(:,:))

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGW(:,:) = eGWlin(:,:)

  else 
  
  ! Find graphical solution of the QP equation

    do is=1,nspin
      call unrestricted_QP_graph(nBas,nC(is),nO(is),nV(is),nR(is),nS_sc,eta,eHF(:,is),SigX(:,is),Vxc(:,is), & 
                                 OmRPA,rho_RPA(:,:,:,is),eGWlin(:,is),eGW(:,is))
    end do
 
  end if

! Compute RPA correlation energy

  call unrestricted_linear_response(ispin,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, &
                                    eGW,ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)

! Dump results

  call print_UG0W0(nBas,nO,eHF,ENuc,EUHF,SigC,Z,eGW,EcRPA,EcGM)

! Free memory

  deallocate(OmRPA,XpY_RPA,XmY_RPA,rho_RPA)

! Perform BSE calculation

  if(BSE) then

    call unrestricted_Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,spin_conserved,spin_flip,eta,nBas,nC,nO,nV,nR,nS,S, &
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

      call unrestricted_ACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,spin_conserved,spin_flip,eta, & 
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

end subroutine UG0W0
