subroutine GG0T0pp(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,doppBSE, & 
                   linearize,eta,regularize,nOrb,nC,nO,nV,nR,nS,ENuc,EGHF,ERI,dipole_int,eHF)

! Perform one-shot calculation with a T-matrix self-energy (G0T0)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dotest

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dophBSE
  logical,intent(in)            :: doppBSE
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: dBSE
  logical,intent(in)            :: dTDA
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: EGHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  logical                       :: print_T = .false.

  integer                       :: nOO
  integer                       :: nVV
  double precision              :: EcRPA
  double precision              :: EcBSE
  double precision              :: EcGM
  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)
  double precision,allocatable  :: Om1(:)
  double precision,allocatable  :: X1(:,:)
  double precision,allocatable  :: Y1(:,:)
  double precision,allocatable  :: rho1(:,:,:)
  double precision,allocatable  :: Om2(:)
  double precision,allocatable  :: X2(:,:)
  double precision,allocatable  :: Y2(:,:)
  double precision,allocatable  :: rho2(:,:,:)
  double precision,allocatable  :: Sig(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eGTlin(:)
  double precision, allocatable :: Om(:), R(:,:)


! Output variables

! Hello world

  write(*,*)
  write(*,*)'**********************************'
  write(*,*)'* Generalized G0T0pp Calculation *'
  write(*,*)'**********************************'
  write(*,*)


! TDA for T

  if(TDA_T) then
    write(*,*) 'Tamm-Dancoff approximation activated for pp T-matrix!'
    write(*,*)
  end if

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Dimensions of the pp-RPA linear reponse matrices

  nOO = nO*(nO - 1)/2
  nVV = nV*(nV - 1)/2

! Memory allocation

  allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO), & 
           rho1(nOrb,nOrb,nVV),rho2(nOrb,nOrb,nOO),Sig(nOrb),Z(nOrb),eGT(nOrb),eGTlin(nOrb))

!----------------------------------------------
! Compute T-matrix
!----------------------------------------------

! Compute linear response

  allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))

                 call ppGLR_C(nOrb,nC,nO,nV,nR,nVV,1d0,eHF,ERI,Cpp)
                 call ppGLR_D(nOrb,nC,nO,nV,nR,nOO,1d0,eHF,ERI,Dpp)
  if(.not.TDA_T) call ppGLR_B(nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)

  call ppGLR(TDA_T,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA)

  deallocate(Bpp,Cpp,Dpp)

  if(print_T) call print_excitation_energies('ppRPA@GHF','2p',nVV,Om1)
  if(print_T) call print_excitation_energies('ppRPA@FHF','2h',nOO,Om2)

!----------------------------------------------
! Compute excitation densities
!----------------------------------------------

  call GGTpp_excitation_density(nOrb,nC,nO,nV,nR,nOO,nVV,ERI,X1,Y1,rho1,X2,Y2,rho2)

!----------------------------------------------
! Compute T-matrix version of the self-energy 
!----------------------------------------------

  if(regularize) call GTpp_regularization(nOrb,nC,nO,nV,nR,nOO,nVV,eHF,Om1,rho1,Om2,rho2)

  call GGTpp_self_energy_diag(eta,nOrb,nC,nO,nV,nR,nOO,nVV,eHF,Om1,rho1,Om2,rho2,EcGM,Sig,Z)

!----------------------------------------------
! Solve the quasi-particle equation
!----------------------------------------------

  eGTlin(:) = eHF(:) + Z(:)*Sig(:)

  if(linearize) then

    write(*,*) ' *** Quasiparticle energies obtained by linearization *** '
    write(*,*)

    eGT(:) = eGTlin(:)

  else

    write(*,*) ' *** Quasiparticle energies obtained by root search *** '
    write(*,*)
     
   call GGTpp_QP_graph(eta,nOrb,nC,nO,nV,nR,nOO,nVV,eHF,Om1,rho1,Om2,rho2,eGTlin,eHF,eGT,Z)

  end if

! call GGTpp_plot_self_energy(nOrb,nC,nO,nV,nR,nOO,nVV,eHF,eGT,Om1,rho1,Om2,rho2)

!----------------------------------------------
! Dump results
!----------------------------------------------

! Compute the ppRPA correlation energy

  allocate(Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))

                 call ppGLR_C(nOrb,nC,nO,nV,nR,nVV,1d0,eGT,ERI,Cpp)
                 call ppGLR_D(nOrb,nC,nO,nV,nR,nOO,1d0,eGT,ERI,Dpp)
  if(.not.TDA_T) call ppGLR_B(nOrb,nC,nO,nV,nR,nOO,nVV,1d0,ERI,Bpp)

  call ppGLR(TDA_T,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcRPA)

  deallocate(Bpp,Cpp,Dpp)

  call print_GG0T0pp(nOrb,nO,eHF,ENuc,EGHF,Sig,Z,eGT,EcGM,EcRPA)

! Perform BSE calculation

! if(dophBSE) then

!   call GGTpp_phBSE(TDA_T,TDA,dBSE,dTDA,eta,nOrb,nC,nO,nV,nR,nS,nOO,nVV, &
!                   Om1,X1,Y1,Om2,X2,Y2,rho1,rho2,ERI,dipole_int,eHF,eGT,EcBSE)

!   if(exchange_kernel) then

!     EcBSE = 0.5d0*EcBSE

!   end if

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@G0T0pp@GHF correlation energy       = ',EcBSE,' au'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@G0T0pp@GHF total energy             = ',ENuc + EGHF + EcBSE,' au'
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

!   if(doACFDT) then

!     write(*,*) '--------------------------------------------------------'
!     write(*,*) 'Adiabatic connection version of phBSE correlation energy'
!     write(*,*) '--------------------------------------------------------'
!     write(*,*)

!     if(doXBS) then

!       write(*,*) '*** scaled screening version (XBS) ***'
!       write(*,*)

!     end if

!     call RGTpp_phACFDT(exchange_kernel,doXBS,.false.,TDA_T,TDA,dophBSE,eta,nOrb,nC,nO,nV,nR,nS, &
!                       nOO,nVV,Om1,X1,Y1,Om2,X2,Y2,rho1,rho2,ERI,eHF,eGT,EcBSE)

!     if(exchange_kernel) then

!       EcBSE = 0.5d0*EcBSE

!     end if

!     write(*,*)
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0T0pp@GHF correlation energy           = ',EcBSE,' au'
!     write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0T0pp@GHF total       energy           = ',ENuc + EGHF + EcBSE,' au'
!     write(*,*)'-------------------------------------------------------------------------------'
!     write(*,*)

!   end if

! end if

! if(doppBSE) then

!   call GGTpp_ppBSE(TDA_T,TDA,dBSE,dTDA,eta,nOrb,nC,nO,nV,nR,nOO,nVV, &
!                   Om1,X1,Y1,Om2,X2,Y2,rho1,rho2,ERI,dipole_int,eHF,eGT,EcBSE)

!   write(*,*)
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0T0pp@GHF correlation energy           = ',EcBSE,' au'
!   write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0T0pp@GHF total       energy           = ',ENuc + EGHF + EcBSE,' au'
!   write(*,*)'-------------------------------------------------------------------------------'
!   write(*,*)

! end if

! Testing zone

  if(dotest) then
  
    call dump_test_value('G','G0T0pp correlation energy',EcRPA)
    call dump_test_value('G','G0T0pp HOMO energy',eGT(nO))
    call dump_test_value('G','G0T0pp LUMO energy',eGT(nO+1))

  end if

end subroutine 
