subroutine G0T0(doACFDT,exchange_kernel,doXBS,BSE,TDA,singlet_manifold,triplet_manifold, & 
                linearize,eta,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,eHF)

! Perform one-shot calculation with a T-matrix self-energy (G0T0)

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: doACFDT
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: doXBS
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  logical,intent(in)            :: linearize
  double precision,intent(in)   :: eta

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

! Local variables

  integer                       :: ispin
  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcAC(nspin)
  double precision,allocatable  :: Omega1s(:),Omega1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Omega2s(:),Omega2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)
  double precision,allocatable  :: rho1st(:,:,:),rho2st(:,:,:)
  double precision,allocatable  :: SigT(:)
  double precision,allocatable  :: Z(:)

  double precision,allocatable  :: eG0T0(:)

  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)
  double precision,allocatable  :: XmY(:,:,:)
  double precision,allocatable  :: rho(:,:,:,:)

! Output variables

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0T0 calculation           |'
  write(*,*)'************************************************'
  write(*,*)

! Dimensions of the rr-RPA linear reponse matrices

  nOOs = nO*(nO + 1)/2
  nVVs = nV*(nV + 1)/2

  nOOt = nO*(nO - 1)/2
  nVVt = nV*(nV - 1)/2

! Memory allocation

  allocate(Omega1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs), & 
           Omega2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs), & 
           rho1s(nBas,nO,nVVs),rho2s(nBas,nV,nOOs), & 
           Omega1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt), & 
           Omega2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt), & 
           rho1t(nBas,nO,nVVt),rho2t(nBas,nV,nOOt), & 
           rho1st(nBas,nO,nVVt),rho2st(nBas,nV,nOOt), & 
           SigT(nBas),Z(nBas),eG0T0(nBas))

!----------------------------------------------
! Singlet manifold
!----------------------------------------------

 ispin = 1

! Compute linear response

  call linear_response_pp(ispin,.true.,.false.,nBas,nC,nO,nV,nR,nOOs,nVVs,eHF(:),ERI(:,:,:,:),  & 
                          Omega1s(:),X1s(:,:),Y1s(:,:),Omega2s(:),X2s(:,:),Y2s(:,:),   & 
                          EcRPA(ispin))

  EcRPA(ispin) = 1d0*EcRPA(ispin)

  call print_excitation('pp-RPA (N+2)',ispin,nVVs,Omega1s(:))
  call print_excitation('pp-RPA (N-2)',ispin,nOOs,Omega2s(:))

!----------------------------------------------
! Triplet manifold
!----------------------------------------------

 ispin = 2

! Compute linear response

  call linear_response_pp(ispin,.true.,.false.,nBas,nC,nO,nV,nR,nOOt,nVVt,eHF(:),ERI(:,:,:,:),  & 
                          Omega1t(:),X1t(:,:),Y1t(:,:),Omega2t(:),X2t(:,:),Y2t(:,:),   & 
                          EcRPA(ispin))

  EcRPA(ispin) = 3d0*EcRPA(ispin)

  call print_excitation('pp-RPA (N+2)',ispin,nVVt,Omega1t(:))
  call print_excitation('pp-RPA (N-2)',ispin,nOOt,Omega2t(:))

!-----------------------------------------------------------
! Compute excitation densities for the T-matrix self-energy
!-----------------------------------------------------------

  call excitation_density_Tmatrix(nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,ERI(:,:,:,:),rho1st(:,:,:),rho2st(:,:,:), &
                                  X1s(:,:),Y1s(:,:),rho1s(:,:,:),X2s(:,:),Y2s(:,:),rho2s(:,:,:), &
                                  X1t(:,:),Y1t(:,:),rho1t(:,:,:),X2t(:,:),Y2t(:,:),rho2t(:,:,:))

!----------------------------------------------
! Compute T-matrix version of the self-energy 
!----------------------------------------------

! rho2s(:,:,:) = 0d0
! rho2t(:,:,:) = 0d0
! rho2st(:,:,:) = 0d0

  call self_energy_Tmatrix_diag(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF(:), & 
                                Omega1s(:),rho1s(:,:,:),Omega2s(:),rho2s(:,:,:), & 
                                Omega1t(:),rho1t(:,:,:),Omega2t(:),rho2t(:,:,:), & 
                                rho1st(:,:,:),rho2st(:,:,:),SigT(:))

! Compute renormalization factor for T-matrix self-energy

  call renormalization_factor_Tmatrix(eta,nBas,nC,nO,nV,nR,nOOs,nVVs,nOOt,nVVt,eHF(:), & 
                                      Omega1s(:),rho1s(:,:,:),Omega2s(:),rho2s(:,:,:), & 
                                      Omega1t(:),rho1t(:,:,:),Omega2t(:),rho2t(:,:,:), & 
                                      rho1st(:,:,:),rho2st(:,:,:),Z(:))

!----------------------------------------------
! Solve the quasi-particle equation
!----------------------------------------------

  if(linearize) then

    eG0T0(:) = eHF(:) + Z(:)*SigT(:)

  else
  
    eG0T0(:) = eHF(:) + SigT(:)

  end if

!----------------------------------------------
! Dump results
!----------------------------------------------

  call print_G0T0(nBas,nO,eHF(:),ENuc,ERHF,SigT(:),Z(:),eG0T0(:),EcRPA(:))

! Perform BSE calculation

  if(BSE) then

     allocate(Omega(nS,nspin),XpY(nS,nS,nspin),XmY(nS,nS,nspin),rho(nBas,nBas,nS,nspin))

    call Bethe_Salpeter(TDA,singlet_manifold,triplet_manifold,eta, &
                        nBas,nC,nO,nV,nR,nS,ERI,eHF,eG0T0,Omega,XpY,XmY,rho,EcRPA,EcBSE)

    if(exchange_kernel) then

      EcRPA(1) = 0.5d0*EcRPA(1)
      EcRPA(2) = 1.5d0*EcRPA(1)

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@G0T0 correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@G0T0 correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@G0T0 correlation energy           =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A50,F20.10)') 'Tr@BSE@G0T0 total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(doACFDT) then

      write(*,*) '------------------------------------------------------'
      write(*,*) 'Adiabatic connection version of BSE correlation energy'
      write(*,*) '------------------------------------------------------'
      write(*,*)

      if(doXBS) then

        write(*,*) '*** scaled screening version (XBS) ***'
        write(*,*)

      end if

      call ACFDT(exchange_kernel,doXBS,.true.,TDA,BSE,singlet_manifold,triplet_manifold,eta, &
                 nBas,nC,nO,nV,nR,nS,ERI,eHF,eG0T0,Omega,XpY,XmY,rho,EcAC)

      if(exchange_kernel) then

        EcAC(1) = 0.5d0*EcAC(1)
        EcAC(2) = 1.5d0*EcAC(1)

      end if

      write(*,*)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,'(2X,A50,F20.10)') 'AC@BSE@G0T0 correlation energy (singlet) =',EcAC(1)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@G0T0 correlation energy (triplet) =',EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@G0T0 correlation energy           =',EcAC(1) + EcAC(2)
      write(*,'(2X,A50,F20.10)') 'AC@BSE@G0T0 total energy                 =',ENuc + ERHF + EcAC(1) + EcAC(2)
      write(*,*)'-------------------------------------------------------------------------------'
      write(*,*)

    end if

  end if

end subroutine G0T0
