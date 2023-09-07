subroutine G0T0pp(doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,doppBSE,singlet,triplet, & 
                  linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,dipole_int,eHF)

! Perform one-shot calculation with a T-matrix self-energy (G0T0)

  implicit none
  include 'parameters.h'

! Input variables

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
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin
  integer                       :: iblock
  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
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

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0T0pp calculation           |'
  write(*,*)'************************************************'
  write(*,*)


! TDA for T

  if(TDA_T) then
    write(*,*) 'Tamm-Dancoff approximation for pp T-matrix!'
    write(*,*)
  end if

! TDA 

  if(TDA) then
    write(*,*) 'Tamm-Dancoff approximation activated!'
    write(*,*)
  end if

! Dimensions of the pp-RPA linear reponse matrices

! nOOs = nO*(nO + 1)/2
! nVVs = nV*(nV + 1)/2

  nOOs = nO*nO
  nVVs = nV*nV

  nOOt = nO*(nO - 1)/2
  nVVt = nV*(nV - 1)/2

! Memory allocation

  allocate(Om1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs),    & 
           Om2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs),    & 
           rho1s(nBas,nBas,nVVs),rho2s(nBas,nBas,nOOs), & 
           Om1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt),    & 
           Om2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt),    & 
           rho1t(nBas,nBas,nVVt),rho2t(nBas,nBas,nOOt), & 
           Sig(nBas),Z(nBas),eGT(nBas),eGTlin(nBas))

!----------------------------------------------
! alpha-beta block
!----------------------------------------------

  ispin  = 1
! iblock = 1
  iblock = 3

! Compute linear response

  allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

  if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                 call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVs,1d0,eHF,ERI,Cpp)
                 call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOs,1d0,eHF,ERI,Dpp)

  call ppLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

  call print_excitation_energies('ppRPA@HF (N+2)',iblock,nVVs,Om1s(:))
  call print_excitation_energies('ppRPA@HF (N-2)',iblock,nOOs,Om2s(:))

!----------------------------------------------
! alpha-alpha block
!----------------------------------------------

  ispin  = 2
! iblock = 2
  iblock = 4

! Compute linear response

  allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))

  if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                 call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVt,1d0,eHF,ERI,Cpp)
                 call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOt,1d0,eHF,ERI,Dpp)

  call ppLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

  call print_excitation_energies('ppRPA (N+2)',iblock,nVVt,Om1t)
  call print_excitation_energies('ppRPA (N-2)',iblock,nOOt,Om2t)

!----------------------------------------------
! Compute excitation densities
!----------------------------------------------

! iblock = 1
  iblock = 3
  call GTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

! iblock = 2
  iblock = 4
  call GTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

!----------------------------------------------
! Compute T-matrix version of the self-energy 
!----------------------------------------------

  if(regularize) then 
    call GTpp_regularization(nBas,nC,nO,nV,nR,nOOs,nVVs,eHF,Om1s,rho1s,Om2s,rho2s)
    call GTpp_regularization(nBas,nC,nO,nV,nR,nOOt,nVVt,eHF,Om1t,rho1t,Om2t,rho2t)
  end if

  call GTpp_self_energy_diag(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF,Om1s,rho1s,Om2s,rho2s, & 
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

    write(*,*) ' *** Quasiparticle energies obtained by root search (experimental) *** '
    write(*,*)
     
   call GTpp_QP_graph(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF,Om1s,rho1s,Om2s,rho2s, & 
                      Om1t,rho1t,Om2t,rho2t,eGTlin,eGT,Z)

  end if

! call GTpp_plot_self_energy(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF,eGT,Om1s,rho1s,Om2s,rho2s, &
!                            Om1t,rho1t,Om2t,rho2t)

!----------------------------------------------
! Dump results
!----------------------------------------------

! Compute the ppRPA correlation energy

  ispin  = 1
! iblock = 1
  iblock = 3

  allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

  if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,1d0,ERI,Bpp)
                 call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVs,1d0,eGT,ERI,Cpp)
                 call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOs,1d0,eGT,ERI,Dpp)

  call ppLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

  ispin  = 2
! iblock = 2
  iblock = 4

  allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))

  if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,1d0,ERI,Bpp)
                 call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVt,1d0,eGT,ERI,Cpp)
                 call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOt,1d0,eGT,ERI,Dpp)

  call ppLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

  EcRPA(1) = EcRPA(1) - EcRPA(2)
  EcRPA(2) = 3d0*EcRPA(2)

  call print_G0T0pp(nBas,nO,eHF,ENuc,ERHF,Sig,Z,eGT,EcGM,EcRPA)

! Perform BSE calculation

  if(dophBSE) then

    call GTpp_phBSE(TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,nOOt,nVVt,   & 
                    Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,Om2t,X2t,Y2t,rho1t,rho2t, &
                    ERI,dipole_int,eHF,eGT,EcBSE)

    if(exchange_kernel) then

      EcBSE(1) = 0.5d0*EcBSE(1)
      EcBSE(2) = 1.5d0*EcBSE(1)

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@G0T0pp correlation energy (singlet) =',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@G0T0pp correlation energy (triplet) =',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@G0T0pp correlation energy           =',EcBSE(1) + EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@phBSE@G0T0pp total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '--------------------------------------------------------'
      write(*,*) 'Adiabatic connection version of phBSE correlation energy'
      write(*,*) '--------------------------------------------------------'
      write(*,*)

      if(doXBS) then

        write(*,*) '*** scaled screening version (XBS) ***'
        write(*,*)

      end if

      call GTpp_phACFDT(exchange_kernel,doXBS,.false.,TDA_T,TDA,dophBSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS, &
                        nOOs,nVVs,nOOt,nVVt,Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,     & 
                        Om2t,X2t,Y2t,rho1t,rho2t,ERI,eHF,eGT,EcBSE)

      if(exchange_kernel) then

        EcBSE(1) = 0.5d0*EcBSE(1)
        EcBSE(2) = 1.5d0*EcBSE(2)

      end if

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0T0pp correlation energy (singlet) =',EcBSE(1),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0T0pp correlation energy (triplet) =',EcBSE(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0T0pp correlation energy           =',EcBSE(1) + EcBSE(2),' au'
      write(*,'(2X,A50,F20.10,A3)') 'AC@phBSE@G0T0pp total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2),' au'
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

  if(doppBSE) then

    call GTpp_ppBSE(TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,      &
                    Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,Om2t,X2t,Y2t,rho1t,rho2t, &
                    ERI,dipole_int,eHF,eGT,EcBSE)

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0T0pp correlation energy (singlet) =',EcBSE(1),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0T0pp correlation energy (triplet) =',EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0T0pp correlation energy =',EcBSE(1) + EcBSE(2),' au'
    write(*,'(2X,A50,F20.10,A3)') 'Tr@ppBSE@G0T0pp total energy =',ENuc + ERHF + EcBSE(1) + EcBSE(2),' au'
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine 
