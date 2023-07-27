subroutine G0T0pp(doACFDT,exchange_kernel,doXBS,dophBSE,TDA_T,TDA,dBSE,dTDA,doppBSE,singlet,triplet, & 
                  linearize,eta,regularize,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI_MO,dipole_int,PHF,cHF,eHF)

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
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: ERI_MO(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)

! Local variables

  integer                       :: ispin
  integer                       :: iblock
  integer                       :: nOOab,nOOaa
  integer                       :: nVVab,nVVaa
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)
  double precision,allocatable  :: Om1ab(:),Om1aa(:)
  double precision,allocatable  :: X1ab(:,:),X1aa(:,:)
  double precision,allocatable  :: Y1ab(:,:),Y1aa(:,:)
  double precision,allocatable  :: rho1ab(:,:,:),rho1aa(:,:,:)
  double precision,allocatable  :: Om2ab(:),Om2aa(:)
  double precision,allocatable  :: X2ab(:,:),X2aa(:,:)
  double precision,allocatable  :: Y2ab(:,:),Y2aa(:,:)
  double precision,allocatable  :: rho2ab(:,:,:),rho2aa(:,:,:)
  double precision,allocatable  :: Sig(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: eGT(:)
  double precision,allocatable  :: eGTlin(:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0T0 calculation           |'
  write(*,*)'************************************************'
  write(*,*)

! Dimensions of the pp-RPA linear reponse matrices

  nOOab = nO*nO
  nVVab = nV*nV

  nOOaa = nO*(nO - 1)/2
  nVVaa = nV*(nV - 1)/2

! Memory allocation

  allocate(Om1ab(nVVab),X1ab(nVVab,nVVab),Y1ab(nOOab,nVVab), & 
           Om2ab(nOOab),X2ab(nVVab,nOOab),Y2ab(nOOab,nOOab), & 
           rho1ab(nBas,nBas,nVVab),rho2ab(nBas,nBas,nOOab),  & 
           Om1aa(nVVaa),X1aa(nVVaa,nVVaa),Y1aa(nOOaa,nVVaa), & 
           Om2aa(nOOaa),X2aa(nVVaa,nOOaa),Y2aa(nOOaa,nOOaa), & 
           rho1aa(nBas,nBas,nVVaa),rho2aa(nBas,nBas,nOOaa),  & 
           Sig(nBas),Z(nBas),eGT(nBas),eGTlin(nBas))

!----------------------------------------------
! alpha-beta block
!----------------------------------------------

  ispin  = 1
  iblock = 3

! Compute linear response

  allocate(Bpp(nVVab,nOOab),Cpp(nVVab,nVVab),Dpp(nOOab,nOOab))

  if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOab,nVVab,1d0,ERI_MO,Bpp)
                 call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVab,1d0,eHF,ERI_MO,Cpp)
                 call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOab,1d0,eHF,ERI_MO,Dpp)

  call ppLR(TDA_T,nOOab,nVVab,Bpp,Cpp,Dpp,Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

  call print_excitation('pp-RPA (N+2)',iblock,nVVab,Om1ab(:))
  call print_excitation('pp-RPA (N-2)',iblock,nOOab,Om2ab(:))

!----------------------------------------------
! alpha-alpha block
!----------------------------------------------

  ispin  = 2
  iblock = 4

! Compute linear response

  allocate(Bpp(nVVaa,nOOaa),Cpp(nVVaa,nVVaa),Dpp(nOOaa,nOOaa))

  if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOaa,nVVaa,1d0,ERI_MO,Bpp)
                 call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVaa,1d0,eHF,ERI_MO,Cpp)
                 call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOaa,1d0,eHF,ERI_MO,Dpp)

  call ppLR(TDA_T,nOOaa,nVVaa,Bpp,Cpp,Dpp,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

  call print_excitation('pp-RPA (N+2)',iblock,nVVaa,Om1aa(:))
  call print_excitation('pp-RPA (N-2)',iblock,nOOaa,Om2aa(:))

!----------------------------------------------
! Compute T-matrix version of the self-energy 
!----------------------------------------------

  EcGM   = 0d0
  Sig(:) = 0d0
  Z(:)   = 0d0

  iblock = 3

  call GTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOab,nVVab,ERI_MO,X1ab,Y1ab,rho1ab,X2ab,Y2ab,rho2ab)

  if(regularize) then

    call regularized_self_energy_Tmatrix_diag(eta,nBas,nC,nO,nV,nR,nOOab,nVVab,eHF,Om1ab,rho1ab,Om2ab,rho2ab,EcGM,Sig)
    call regularized_renormalization_factor_Tmatrix(eta,nBas,nC,nO,nV,nR,nOOab,nVVab,eHF,Om1ab,rho1ab,Om2ab,rho2ab,Z)

  else

    call GTpp_self_energy_diag(eta,nBas,nC,nO,nV,nR,nOOab,nVVab,eHF,Om1ab,rho1ab,Om2ab,rho2ab,EcGM,Sig,Z)

  end if

  iblock = 4

  call GTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOaa,nVVaa,ERI_MO,X1aa,Y1aa,rho1aa,X2aa,Y2aa,rho2aa)

  if(regularize) then

    call regularized_self_energy_Tmatrix_diag(eta,nBas,nC,nO,nV,nR,nOOaa,nVVaa,eHF,Om1aa,rho1aa,Om2aa,rho2aa,EcGM,Sig)
    call regularized_renormalization_factor_Tmatrix(eta,nBas,nC,nO,nV,nR,nOOaa,nVVaa,eHF,Om1aa,rho1aa,Om2aa,rho2aa,Z)

  else

    call GTpp_self_energy_diag(eta,nBas,nC,nO,nV,nR,nOOaa,nVVaa,eHF,Om1aa,rho1aa,Om2aa,rho2aa,EcGM,Sig,Z)

  end if

  Z(:) = 1d0/(1d0 - Z(:))

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
     
   call GTpp_QP_graph(eta,nBas,nC,nO,nV,nR,nOOaa,nVVaa,eHF,Om1aa,rho1aa,Om2aa,rho2aa,eGTlin,eGT)

  end if

!----------------------------------------------
! Dump results
!----------------------------------------------

! Compute the ppRPA correlation energy

  ispin  = 1
  iblock = 3

  allocate(Bpp(nVVab,nOOab),Cpp(nVVab,nVVab),Dpp(nOOab,nOOab))

  if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOab,nVVab,1d0,ERI_MO,Bpp)
                 call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVab,1d0,eGT,ERI_MO,Cpp)
                 call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOab,1d0,eGT,ERI_MO,Dpp)

  call ppLR(TDA_T,nOOab,nVVab,Bpp,Cpp,Dpp,Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

  ispin  = 2
  iblock = 4

  allocate(Bpp(nVVaa,nOOaa),Cpp(nVVaa,nVVaa),Dpp(nOOaa,nOOaa))

  if(.not.TDA_T) call ppLR_B(iblock,nBas,nC,nO,nV,nR,nOOaa,nVVaa,1d0,ERI_MO,Bpp)
                 call ppLR_C(iblock,nBas,nC,nO,nV,nR,nVVaa,1d0,eGT,ERI_MO,Cpp)
                 call ppLR_D(iblock,nBas,nC,nO,nV,nR,nOOaa,1d0,eGT,ERI_MO,Dpp)

  call ppLR(TDA_T,nOOaa,nVVaa,Bpp,Cpp,Dpp,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,EcRPA(ispin))

  deallocate(Bpp,Cpp,Dpp)

  EcRPA(1) = EcRPA(1) - EcRPA(2)
  EcRPA(2) = 3d0*EcRPA(2)

  call print_G0T0pp(nBas,nO,eHF,ENuc,ERHF,Sig,Z,eGT,EcGM,EcRPA)

! Perform BSE calculation

  if(dophBSE) then

    call GTpp_phBSE(TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,nOOab,nVVab,nOOaa,nVVaa,   & 
                    Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,rho1ab,rho2ab,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,rho1aa,rho2aa, &
                    ERI_MO,dipole_int,eHF,eGT,EcBSE)

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
                        nOOab,nVVab,nOOaa,nVVaa,Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,rho1ab,rho2ab,Om1aa,X1aa,Y1aa,     & 
                        Om2aa,X2aa,Y2aa,rho1aa,rho2aa,ERI_MO,eHF,eGT,EcBSE)

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

    call GTpp_ppBSE(TDA_T,TDA,dBSE,dTDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nOOab,nVVab,nOOaa,nVVaa,      &
                    Om1ab,X1ab,Y1ab,Om2ab,X2ab,Y2ab,rho1ab,rho2ab,Om1aa,X1aa,Y1aa,Om2aa,X2aa,Y2aa,rho1aa,rho2aa, &
                    ERI_MO,dipole_int,eHF,eGT,EcBSE)

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
