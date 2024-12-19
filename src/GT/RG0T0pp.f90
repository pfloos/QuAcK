subroutine RG0T0pp(dotest,doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,doppBSE,singlet,triplet, & 
                   linearize,eta,regularize,nBas,nOrb,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

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
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta
  logical,intent(in)            :: regularize

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nOrb
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nOrb)
  double precision,intent(in)   :: ERI(nOrb,nOrb,nOrb,nOrb)
  double precision,intent(in)   :: dipole_int(nOrb,nOrb,ncart)

! Local variables

  logical                       :: print_T   = .true.
  logical                       :: plot_self = .false.

  integer                       :: isp_T
  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  integer                       :: n_states, n_states_diag
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)
  double precision,allocatable  :: Om1s(:),Om1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Om2s(:),Om2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)
  double precision,allocatable  :: Sig(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eGTlin(:)
  integer,          allocatable :: supp_data_int(:)
  double precision, allocatable :: supp_data_dbl(:)
  double precision, allocatable :: Om(:), R(:,:)


! Output variables

! Hello world

  write(*,*)
  write(*,*)'*********************************'
  write(*,*)'* Restricted G0T0pp Calculation *'
  write(*,*)'*********************************'
  write(*,*)

! TDA for T

  if(TDA_T) then
    write(*,*) 'Tamm-Dancoff approximation activated for pp T-matrix!'
    write(*,*)
  end if

! Dimensions of the pp-RPA linear reponse matrices

  nOOs = nO*(nO + 1)/2
  nVVs = nV*(nV + 1)/2

  nOOt = nO*(nO - 1)/2
  nVVt = nV*(nV - 1)/2

! Memory allocation

  allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),    & 
           Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs),    & 
           rho1s(nOrb,nOrb,nVVs),rho2s(nOrb,nOrb,nOOs), & 
           Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),    & 
           Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt),    & 
           rho1t(nOrb,nOrb,nVVt),rho2t(nOrb,nOrb,nOOt), & 
           Sig(nOrb),Z(nOrb),eGT(nOrb),eGTlin(nOrb))

!----------------------------------------------
! alpha-beta block
!----------------------------------------------

  isp_T  = 1

! Compute linear response

  allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

  if(.not.TDA_T) call ppRLR_B(isp_T,nOrb,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                 call ppRLR_C(isp_T,nOrb,nC,nO,nV,nR,nVVs,1d0,eHF,ERI,Cpp)
                 call ppRLR_D(isp_T,nOrb,nC,nO,nV,nR,nOOs,1d0,eHF,ERI,Dpp)

  call ppRLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(isp_T))

  deallocate(Bpp,Cpp,Dpp)

  if(print_T) call print_excitation_energies('ppRPA@RHF','2p (singlets)',nVVs,Om1s)
  if(print_T) call print_excitation_energies('ppRPA@RHF','2h (singlets)',nOOs,Om2s)

!----------------------------------------------
! alpha-alpha block
!----------------------------------------------

  isp_T  = 2

! Compute linear response

  allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))

  if(.not.TDA_T) call ppRLR_B(isp_T,nOrb,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                 call ppRLR_C(isp_T,nOrb,nC,nO,nV,nR,nVVt,1d0,eHF,ERI,Cpp)
                 call ppRLR_D(isp_T,nOrb,nC,nO,nV,nR,nOOt,1d0,eHF,ERI,Dpp)

  call ppRLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(isp_T))

  deallocate(Bpp,Cpp,Dpp)

  if(print_T) call print_excitation_energies('ppRPA@RHF','2p (triplets)',nVVt,Om1t)
  if(print_T) call print_excitation_energies('ppRPA@RHF','2h (triplets)',nOOt,Om2t)

!----------------------------------------------
! Compute excitation densities
!----------------------------------------------

  isp_T = 1

  call RGTpp_excitation_density(isp_T,nOrb,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

  isp_T = 2

  call RGTpp_excitation_density(isp_T,nOrb,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

!----------------------------------------------
! Compute T-matrix version of the self-energy 
!----------------------------------------------

  if(regularize) then 
    call GTpp_regularization(nOrb,nC,nO,nV,nR,nOOs,nVVs,eHF,Om1s,rho1s,Om2s,rho2s)
    call GTpp_regularization(nOrb,nC,nO,nV,nR,nOOt,nVVt,eHF,Om1t,rho1t,Om2t,rho2t)
  end if

  call RGTpp_self_energy_diag(eta,nOrb,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF,Om1s,rho1s,Om2s,rho2s, & 
                              Om1t,rho1t,Om2t,rho2t,EcGM,Sig,Z)

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
     
   call RGTpp_QP_graph(eta,nOrb,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF,Om1s,rho1s,Om2s,rho2s, & 
                       Om1t,rho1t,Om2t,rho2t,eGTlin,eHF,eGT,Z)

  end if

  if(plot_self) call RGTpp_plot_self_energy(nOrb,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF,eGT,Om1s,rho1s,Om2s,rho2s, &
                                            Om1t,rho1t,Om2t,rho2t)

!----------------------------------------------
! Dump results
!----------------------------------------------

! Compute the ppRPA correlation energy

  isp_T  = 1

  allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

  if(.not.TDA_T) call ppRLR_B(isp_T,nOrb,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                 call ppRLR_C(isp_T,nOrb,nC,nO,nV,nR,nVVs,1d0,eGT,ERI,Cpp)
                 call ppRLR_D(isp_T,nOrb,nC,nO,nV,nR,nOOs,1d0,eGT,ERI,Dpp)

  call ppRLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(isp_T))
!  call ppLR_transition_vectors(.true.,nBas,nC,nO,nV,nR,nOOs,nVVs,dipole_int,Om1s,X1s,Y1s,Om2s,X2s,Y2s)

  deallocate(Bpp,Cpp,Dpp)

  isp_T  = 2

  allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))

  if(.not.TDA_T) call ppRLR_B(isp_T,nOrb,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                 call ppRLR_C(isp_T,nOrb,nC,nO,nV,nR,nVVt,1d0,eGT,ERI,Cpp)
                 call ppRLR_D(isp_T,nOrb,nC,nO,nV,nR,nOOt,1d0,eGT,ERI,Dpp)

  call ppRLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(isp_T))
!  call ppLR_transition_vectors(.false.,nBas,nC,nO,nV,nR,nOOt,nVVt,dipole_int,Om1t,X1t,Y1t,Om2t,X2t,Y2t)

  deallocate(Bpp,Cpp,Dpp)

  EcRPA(1) = 1d0*EcRPA(1)
  EcRPA(2) = 3d0*EcRPA(2)

  call print_RG0T0pp(nOrb,nO,eHF,ENuc,ERHF,Sig,Z,eGT,EcGM,EcRPA)

! Perform BSE calculation

  if(dophBSE) then

    call RGTpp_phBSE(exchange_kernel,TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt, & 
                    Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,Om2t,X2t,Y2t,rho1t,rho2t,      &
                    ERI,dipole_int,eHF,eGT,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@G0T0pp@RHF correlation energy (singlet) = ',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@G0T0pp@RHF correlation energy (triplet) = ',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@G0T0pp@RHF correlation energy           = ',sum(EcBSE),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@G0T0pp@RHF total energy                 = ',ENuc + ERHF + sum(EcBSE),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      call RGTpp_phACFDT(exchange_kernel,doXBS,.false.,TDA_T,TDA,dophBSE,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nS, &
                        nOOs,nVVs,nOOt,nVVt,Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,     & 
                        Om2t,X2t,Y2t,rho1t,rho2t,ERI,eHF,eGT,EcBSE)

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0T0pp@RHF correlation energy (singlet) = ',EcBSE(1),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0T0pp@RHF correlation energy (triplet) = ',EcBSE(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0T0pp@RHF correlation energy           = ',sum(EcBSE),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0T0pp@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

  if(doppBSE) then

    call RGTpp_ppBSE(TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nOrb,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt, &
                     ERI,dipole_int,eHF,eGT,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0T0pp@RHF correlation energy (singlet) = ',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0T0pp@RHF correlation energy (triplet) = ',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0T0pp@RHF correlation energy           = ',sum(EcBSE),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0T0pp@RHF total       energy           = ',ENuc + ERHF + sum(EcBSE),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

! Testing zone

  if(dotest) then
  
    call dump_test_value('R','G0T0pp correlation energy',sum(EcRPA))
    call dump_test_value('R','G0T0pp HOMO energy',eGT(nO))
    call dump_test_value('R','G0T0pp LUMO energy',eGT(nO+1))

  end if

end subroutine 
