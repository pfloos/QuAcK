subroutine TDHF(singlet_manifold,triplet_manifold,nBas,nC,nO,nV,nR,nS,ENuc,ERHF,ERI,e)

! Perform random phase approximation calculation

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: ENuc
  double precision,intent(in)   :: ERHF
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dRPA
  logical                       :: TDA
  logical                       :: BSE
  integer                       :: ispin
  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)
  double precision,allocatable  :: XmY(:,:,:)

  double precision              :: rho
  double precision              :: EcRPA(nspin)

  logical                       :: AC
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: EcACRPA(:,:)
  double precision,allocatable  :: EcAC(:,:)

! Hello world

  write(*,*)
  write(*,*)'************************************************'
  write(*,*)'|  Time-dependent Hartree-Fock calculation     |'
  write(*,*)'************************************************'
  write(*,*)

! Initialization

  EcRPA(:) = 0d0

! Switch on exchange for TDHF

  dRPA = .false.
 
! Switch off Tamm-Dancoff approximation for TDHF

  TDA = .false.
 
! Switch off Bethe-Salpeter equation for TDHF

  BSE = .false. 

! Memory allocation

  allocate(Omega(nS,nspin),XpY(nS,nS,nspin),XmY(nS,nS,nspin))

  AC = .true.
  allocate(EcACRPA(nAC,nspin),EcAC(nAC,nspin))

! Singlet manifold

  if(singlet_manifold) then 

    ispin = 1

    call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,1d0,e,ERI,rho, &
                         EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('TDHF ',ispin,nS,Omega(:,ispin))

  endif

! Triplet manifold 

  if(triplet_manifold) then 

    ispin = 2

    call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,1d0,e,ERI,rho, &
                         EcRPA(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
    call print_excitation('TDHF ',ispin,nS,Omega(:,ispin))

  endif

  write(*,*)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,'(2X,A40,F15.6)') 'RPA@TDHF correlation energy (singlet) =',EcRPA(1)
  write(*,'(2X,A40,F15.6)') 'RPA@TDHF correlation energy (triplet) =',EcRPA(2)
  write(*,'(2X,A40,F15.6)') 'RPA@TDHF correlation energy           =',EcRPA(1) + EcRPA(2)
  write(*,'(2X,A40,F15.6)') 'RPA@TDHF total energy                 =',ENuc + ERHF + EcRPA(1) + EcRPA(2)
  write(*,*)'-------------------------------------------------------------------------------'
  write(*,*)

!   Compute the correlation energy via the adiabatic connection 

    if(AC) then

      write(*,*) '------------------------------------------------------'
      write(*,*) 'Adiabatic connection version of RPA correlation energy'
      write(*,*) '------------------------------------------------------'
      write(*,*) 
 
      if(singlet_manifold) then

        ispin = 1
        EcACRPA(:,ispin) = 0d0

        write(*,*) '--------------'
        write(*,*) 'Singlet states'
        write(*,*) '--------------'
        write(*,*) 

        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','EcRPA(lambda)','Tr(V x P_lambda)'
        write(*,*) '-----------------------------------------------------------------------------------'

        do iAC=1,nAC
 
          lambda = rAC(iAC)
  
          call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,lambda,e,ERI, &
                           rho,EcACRPA(iAC,ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

          call Ec_AC(ispin,dRPA,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),EcAC(iAC,ispin))

          write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcACRPA(iAC,ispin),EcAC(iAC,ispin)

        end do

        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,'(2X,A50,1X,F15.6)') ' Ec(RPA) via Gauss-Legendre quadrature:',0.5d0*dot_product(wAC,EcAC(:,ispin))
        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,*)

      end if
 
      if(triplet_manifold) then

        ispin = 2
        EcACRPA(:,ispin) = 0d0

        write(*,*) '--------------'
        write(*,*) 'Triplet states'
        write(*,*) '--------------'
        write(*,*) 

        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','EcRPA(lambda)','Tr(V x P_lambda)'
        write(*,*) '-----------------------------------------------------------------------------------'

        do iAC=1,nAC
 
          lambda = rAC(iAC)
         
          call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,lambda,e,ERI, &
                               rho,EcACRPA(iAC,ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

          call Ec_AC(ispin,dRPA,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),EcAC(iAC,ispin))

          write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcACRPA(iAC,ispin),EcAC(iAC,ispin)

        end do

        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,'(2X,A50,1X,F15.6)') ' Ec(RPA) via Gauss-Legendre quadrature:',0.5d0*dot_product(wAC,EcAC(:,ispin))
        write(*,*) '-----------------------------------------------------------------------------------'
        write(*,*)

      end if
 
    end if

end subroutine TDHF
