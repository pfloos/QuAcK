subroutine G0W0_SOSEX(doACFDT,exchange_kernel,doXBS,BSE,TDA_W,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta, & 
                      nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_AO,ERI_MO,dipole_int,PHF,cHF,eHF,Vxc,eSOSEX)

! Perform the SOSEX extension of G0W0

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
  logical,intent(in)            :: evDyn
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  double precision,intent(in)   :: eta

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
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: SigX(:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: OmRPA(:,:)
  double precision,allocatable  :: XpY_RPA(:,:,:)
  double precision,allocatable  :: XmY_RPA(:,:,:)
  double precision,allocatable  :: rho_RPA(:,:,:,:)

! Output variables

  double precision              :: eSOSEX(nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot SOSEX calculation          |'
  write(*,*)'************************************************'
  write(*,*)

! Initialization

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

! Memory allocation

  allocate(SigC(nBas),SigX(nBas),Z(nBas),OmRPA(nS,nspin),XpY_RPA(nS,nS,nspin),XmY_RPA(nS,nS,nspin),rho_RPA(nBas,nBas,nS,nspin))

!-------------------!
! Compute screening !
!-------------------!

  do ispin=1,nspin
    call linear_response(ispin,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI_MO, & 
                         OmRPA(:,ispin),rho_RPA(:,:,:,ispin),EcRPA(ispin),OmRPA(:,ispin),XpY_RPA(:,:,ispin),XmY_RPA(:,:,ispin))
    if(print_W) call print_excitation('RPA@HF      ',ispin,nS,OmRPA)
  end do

!--------------------------!
! Compute spectral weights !
!--------------------------!

  call excitation_density_SOSEX(nBas,nC,nO,nR,nS,ERI_MO,XpY_RPA,rho_RPA)

!------------------------!
! Compute GW self-energy !
!------------------------!

  call self_energy_exchange_diag(nBas,cHF,PHF,ERI_AO,SigX)
  call self_energy_correlation_SOSEX_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,OmRPA,rho_RPA,EcGM,SigC)

!--------------------------------!
! Compute renormalization factor !
!--------------------------------!

  call renormalization_factor_SOSEX(eta,nBas,nC,nO,nV,nR,nS,eHF,OmRPA,rho_RPA,Z)

!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  eSOSEX(:) = eHF(:) + Z(:)*(SigX(:) + SigC(:) - Vxc(:))

! Compute the RPA correlation energy

  do ispin=1,nspin
    call linear_response(ispin,.true.,TDA_W,.false.,eta,nBas,nC,nO,nV,nR,nS,1d0,eSOSEX,ERI_MO,OmRPA(:,ispin), & 
                         rho_RPA(:,:,:,ispin),EcRPA(ispin),OmRPA(:,ispin),XpY_RPA(:,:,ispin),XmY_RPA(:,:,ispin))
  end do

!--------------!
! Dump results !
!--------------!

  call print_SOSEX(nBas,nO,eHF,ENuc,ERHF,SigC,Z,eSOSEX,EcRPA,EcGM)

! Deallocate memory

  deallocate(SigC,Z,OmRPA,XpY_RPA,XmY_RPA,rho_RPA)

! Perform BSE calculation

  if(BSE) then

    call Bethe_Salpeter(TDA_W,TDA,dBSE,dTDA,evDyn,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,dipole_int,eHF,eSOSEX,EcBSE)

    if(exchange_kernel) then
 
      EcBSE(1) = 0.5d0*EcBSE(1)
      EcBSE(2) = 1.5d0*EcBSE(2)
 
    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@SOSEX correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@SOSEX correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@SOSEX correlation energy           =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@SOSEX total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '--------------------------------------------------------------'
      write(*,*) ' Adiabatic connection version of BSE@SOSEX correlation energy '
      write(*,*) '--------------------------------------------------------------'
      write(*,*) 

      if(doXBS) then 

        write(*,*) '*** scaled screening version (XBS) ***'
        write(*,*)

      end if

      call ACFDT(exchange_kernel,doXBS,.true.,TDA_W,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI_MO,eHF,eSOSEX,EcAC)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10)') 'AC@BSE@SOSEX correlation energy (singlet) =',EcAC(1)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@SOSEX correlation energy (triplet) =',EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@SOSEX correlation energy           =',EcAC(1) + EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@SOSEX total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

end subroutine G0W0_SOSEX
