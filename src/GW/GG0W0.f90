subroutine GG0W0(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,dophBSE2,TDA_W,TDA,dBSE,dTDA,doppBSE, & 
                 linearize,eta,doSRG,do_linDM_GW,nBas,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF,eGW_out)
  
! Perform G0W0 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: dotest

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
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: doSRG
  logical,intent(in)            :: do_linDM_GW

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: eHF(nBas)

! Local variables

  logical                       :: print_W = .true.
  logical                       :: dRPA
  double precision              :: flow
  double precision              :: EcRPA
  double precision              :: EcBSE
  double precision              :: EcGM
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)
  double precision,allocatable  :: rho(:,:,:)

  double precision,allocatable  :: J(:,:)
  double precision,allocatable  :: K(:,:)
  double precision,allocatable  :: F(:,:)
  double precision,allocatable  :: linDM(:,:)
  integer                       :: p
  double precision,allocatable  :: eHFlinDM(:)
  double precision,allocatable  :: occ_nb(:)
  
  double precision,allocatable  :: eGWlin(:)
  double precision,allocatable  :: eGW(:)

! Output variables

  double precision,intent(out)  :: eGW_out(nBas)

! Hello world

  write(*,*)
  write(*,*)'********************************'
  write(*,*)'* Generalized G0W0 Calculation *'
  write(*,*)'********************************'
  write(*,*)

! Initialization

  dRPA = .true.

! TDA for W

  if(TDA_W) then 
    write(*,*) 'Tamm-Dancoff approximation for dynamic screening!'
    write(*,*)
  end if

! SRG regularization

  flow = 500d0

  if(doSRG) then

    write(*,*) '*** SRG regularized G0W0 scheme ***'
    write(*,*)

  end if

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),SigC(nBas),Z(nBas),Om(nS),XpY(nS,nS),XmY(nS,nS),rho(nBas,nBas,nS), & 
           eGW(nBas),eGWlin(nBas))

!-------------------!
! Compute screening !
!-------------------!

                 call phGLR_A(dRPA,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI,Aph)
  if(.not.TDA_W) call phGLR_B(dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)
    call matout(nS,nS, Aph)
    call matout(nS,nS, Bph)

  call phGLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

  if(print_W) call print_excitation_energies('phRPA@GHF','generalized',nS,Om)

!--------------------------!
! Compute spectral weights !
!--------------------------!

  call GGW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY,rho)

!------------------------!
! Compute GW self-energy !
!------------------------!

  if(doSRG) then 
    call GGW_SRG_self_energy_diag(flow,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigC,Z)
  else
    call GGW_self_energy_diag(eta,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,EcGM,SigC,Z,ERI)
  end if

!--------------------------------------!
! Linearized density matrix correction !
!--------------------------------------!

  allocate(linDM(nBas,nBas),eHFlinDM(nBas))
  linDM(:,:) = 0d0
  if(do_linDM_GW) then
     call G_linDM_GW(nBas,nC,nO,nV,nR,nS,eHF,Om,rho,eta,linDM)

     allocate(J(nBas,nBas))
     allocate(K(nBas,nBas))
     allocate(F(nBas,nBas))
     call Hartree_matrix_AO_basis(nBas,linDM,ERI,J)
     call exchange_matrix_AO_basis(nBas,linDM,ERI,K)
     F(:,:) = J(:,:) + K(:,:)

     do p=nC+1,nBas-nR
        eHFlinDM(p) = eHF(p) + F(p,p)
     end do

     do p=nC+1,nO
        linDM(p,p) = linDM(p,p) + 1d0
     end do

     ! allocate(occ_nb(nBas))
     ! occ_nb(:) = 0d0
     ! call diagonalize_matrix(nBas,linDM,occ_nb)
     ! call vecout(nBas,occ_nb)
     ! deallocate(occ_nb)
     
     deallocate(J,K,F,linDM)
  end if
  
!-----------------------------------!
! Solve the quasi-particle equation !
!-----------------------------------!

  ! Linearized or graphical solution?

  if(do_linDM_GW) then
     eGWlin(:) = eHFlinDM(:) + Z(:)*SigC(:)
  else
     eGWlin(:) = eHF(:) + Z(:)*SigC(:)
  end if

  if(linearize) then 
 
    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGW(:) = eGWlin(:)

  else 

    write(*,*) ' *** Quasiparticle energies obtained by root search *** '
    write(*,*)
  
     if(do_linDM_GW) then
        call GGW_QP_graph(doSRG,eta,flow,nBas,nC,nO,nV,nR,nS,eHFlinDM,Om,rho,eGWlin,eHF,eGW,Z)
     else
        call GGW_QP_graph(doSRG,eta,flow,nBas,nC,nO,nV,nR,nS,eHF,Om,rho,eGWlin,eHF,eGW,Z)
     end if
        
  end if

! Compute the RPA correlation energy

                 call phGLR_A(dRPA,nBas,nC,nO,nV,nR,nS,1d0,eGW,ERI,Aph)
  if(.not.TDA_W) call phGLR_B(dRPA,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phGLR(TDA_W,nS,Aph,Bph,EcRPA,Om,XpY,XmY)

!--------------!
! Dump results !
!--------------!

  call print_GG0W0(nBas,nC,nO,nV,nR,eHF,ENuc,EGHF,SigC,Z,eGW,EcRPA,EcGM)

  eGW_out(:) = eGW(:)

! Deallocate memory

  deallocate(SigC,Z,Om,XpY,XmY,rho)

! Perform BSE calculation

  if(dophBSE) then

    call GGW_phBSE(dophBSE2,TDA_W,TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE)
    call GGW_phBSE_omegazero(dophBSE2,TDA_W,TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE) 
    call GGW_phBSE_qs(dophBSE2,TDA_W,TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE) 

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0@GHF correlation energy = ',EcBSE,' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@BSE@G0W0@GHF total energy       = ',ENuc + EGHF + EcBSE,' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

  if(doppBSE) then

    call GGW_ppBSE(TDA_W,TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE)
    !call GGW_ppBSE_omegazero(TDA_W,TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE)
    !call GGW_ppBSE_qs(TDA_W,TDA,dBSE,dTDA,eta,nBas,nC,nO,nV,nR,nS,ERI,dipole_int,eHF,eGW,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0@GHF correlation energy         = ',EcBSE,' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0W0@GHF total energy               = ',ENuc + EGHF + EcBSE,' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Testing zone

  if(dotest) then

    call dump_test_value('G','RPA@G0W0 correlation energy',EcRPA)
    call dump_test_value('G','Tr@ppBSE@G0W0 correlation energy',EcBSE)
    call dump_test_value('G','G0W0 HOMO energy',eGW(nO))
    call dump_test_value('G','G0W0 LUMO energy',eGW(nO+1))

  end if

end subroutine 
