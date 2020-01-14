subroutine G0W0(COHSEX,SOSEX,BSE,TDA,singlet_manifold,triplet_manifold,eta, & 
                nBas,nC,nO,nV,nR,nS,ENuc,ERHF,Hc,H,ERI,PHF,cHF,eHF,eG0W0)

! Perform G0W0 calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: COHSEX
  logical,intent(in)            :: SOSEX
  logical,intent(in)            :: BSE
  logical,intent(in)            :: TDA
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  double precision,intent(in)   :: eta

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: cHF(nBas,nBas)
  double precision,intent(in)   :: PHF(nBas,nBas)
  double precision,intent(in)   :: Hc(nBas,nBas)
  double precision,intent(in)   :: H(nBas,nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dRPA
  integer                       :: ispin,jspin
  double precision              :: EcRPA(nspin)
  double precision              :: EcBSE(nspin)
  double precision              :: EcGM
  double precision,allocatable  :: SigC(:)
  double precision,allocatable  :: Z(:)
  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)
  double precision,allocatable  :: XmY(:,:,:)
  double precision,allocatable  :: rho(:,:,:,:)
  double precision,allocatable  :: rhox(:,:,:,:)

  logical                       :: AC
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: EcACBSE(:,:)
  double precision,allocatable  :: EcAC(:,:)

! Output variables

  double precision              :: eG0W0(nBas)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|          One-shot G0W0 calculation           |'
  write(*,*)'************************************************'
  write(*,*)

! SOSEX correction

  if(SOSEX) write(*,*) 'SOSEX correction activated!'
  write(*,*)

! COHSEX approximation

  if(COHSEX) write(*,*) 'COHSEX approximation activated!'
  write(*,*)

! Switch off exchange for G0W0

  dRPA = .true.

! Spin manifold 

  ispin = 1

! Memory allocation

  allocate(SigC(nBas),Z(nBas),Omega(nS,nspin),XpY(nS,nS,nspin),XmY(nS,nS,nspin),  & 
           rho(nBas,nBas,nS,nspin),rhox(nBas,nBas,nS,nspin))

  AC = .true.
  allocate(EcACBSE(nAC,nspin),EcAC(nAC,nspin))

! Compute linear response

  call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI, & 
                       rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

! Compute correlation part of the self-energy 

  call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

  if(SOSEX) call excitation_density_SOSEX(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rhox(:,:,:,ispin))

  call self_energy_correlation_diag(COHSEX,SOSEX,eta,nBas,nC,nO,nV,nR,nS,eHF, & 
                                    Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),EcGM,SigC)

! Compute renormalization factor

  call renormalization_factor(COHSEX,SOSEX,eta,nBas,nC,nO,nV,nR,nS,eHF, & 
                              Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin),Z(:))

! Solve the quasi-particle equation

  eG0W0(:) = eHF(:) + Z(:)*SigC(:)

! Dump results

  call print_excitation('RPA         ',ispin,nS,Omega(:,ispin))
  call print_G0W0(nBas,nO,eHF,ENuc,ERHF,SigC,Z,eG0W0,EcRPA(ispin),EcGM)

! Plot stuff

  call plot_GW(nBas,nC,nO,nV,nR,nS,eHF,eG0W0,Omega(:,ispin),rho(:,:,:,ispin),rhox(:,:,:,ispin))
 
! Perform BSE calculation

  if(BSE) then

   ! Singlet manifold

   if(singlet_manifold) then

      ispin = 1
      EcBSE(ispin) = 0d0

      call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI, &
                           rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
      call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,1d0,eG0W0,ERI, &
                           rho(:,:,:,ispin),EcBSE(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    end if

   ! Triplet manifold

   if(triplet_manifold) then

      ispin = 2
      EcBSE(ispin) = 0d0

      call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,1d0,eHF,ERI, &
                           rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
      call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,1d0,eG0W0,ERI, &
                           rho(:,:,:,ispin),EcBSE(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
      call print_excitation('BSE  ',ispin,nS,Omega(:,ispin))

    end if

    write(*,*)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,'(2X,A40,F15.6)') 'BSE@G0W0 correlation energy (singlet) =',EcBSE(1)
    write(*,'(2X,A40,F15.6)') 'BSE@G0W0 correlation energy (triplet) =',EcBSE(2)
    write(*,'(2X,A40,F15.6)') 'BSE@G0W0 correlation energy           =',EcBSE(1) + EcBSE(2)
    write(*,'(2X,A40,F15.6)') 'BSE@G0W0 total energy                 =',ENuc + ERHF + EcBSE(1) + EcBSE(2)
    write(*,*)'-------------------------------------------------------------------------------'
    write(*,*)

!   Compute the BSE correlation energy via the adiabatic connection 

    if(AC) then

      write(*,*) '------------------------------------------------------'
      write(*,*) 'Adiabatic connection version of BSE correlation energy'
      write(*,*) '------------------------------------------------------'
      write(*,*) 
 
      if(singlet_manifold) then

        ispin = 1
        EcACBSE(:,ispin) = 0d0

        write(*,*) '--------------'
        write(*,*) 'Singlet states'
        write(*,*) '--------------'
        write(*,*) 

        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','EcBSE(lambda)','Tr(V x P_lambda)'
        write(*,*) '-----------------------------------------------------------------------------------'

        do iAC=1,nAC
 
          lambda = rAC(iAC)

          call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,lambda,eHF,ERI, &
                               rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
          call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))
  
          call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,lambda,eG0W0,ERI, &
                           rho(:,:,:,ispin),EcACBSE(iAC,ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

          call Ec_AC(ispin,dRPA,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),EcAC(iAC,ispin))

          write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcACBSE(iAC,ispin),EcAC(iAC,ispin)

        end do

        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,'(2X,A50,1X,F15.6)') ' Ec(BSE) via Gauss-Legendre quadrature:',0.5d0*dot_product(wAC,EcAC(:,ispin))
        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,*)

      end if
 
      if(triplet_manifold) then

        ispin = 2
        EcACBSE(:,ispin) = 0d0

        write(*,*) '--------------'
        write(*,*) 'Triplet states'
        write(*,*) '--------------'
        write(*,*) 

        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','EcBSE(lambda)','Tr(V x P_lambda)'
        write(*,*) '-----------------------------------------------------------------------------------'

        do iAC=1,nAC
 
          lambda = rAC(iAC)

          call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,lambda,eHF,ERI, &
                               rho(:,:,:,ispin),EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
          call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))
         
          call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,lambda,eG0W0,ERI, &
                               rho(:,:,:,ispin),EcACBSE(iAC,ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

          call Ec_AC(ispin,dRPA,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),EcAC(iAC,ispin))

          write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcACBSE(iAC,ispin),EcAC(iAC,ispin)

        end do

        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,'(2X,A50,1X,F15.6)') ' Ec(BSE) via Gauss-Legendre quadrature:',0.5d0*dot_product(wAC,EcAC(:,ispin))
        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,*)

      end if
 
    end if

  end if

end subroutine G0W0
