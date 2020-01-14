subroutine ACDFT(scaled_screening,dRPA,TDA,BSE,singlet_manifold,triplet_manifold, & 
                 nBas,nC,nO,nV,nR,nS,ERI,e,Omega,XpY,XmY,rho)

! Compute the correlation energy via the adiabatic connection dissipation fluctuation theorem

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: scaled_screening
  logical,intent(in)            :: dRPA
  logical,intent(in)            :: TDA
  logical,intent(in)            :: BSE
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold

  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

  double precision              :: Omega(nS,nspin)
  double precision              :: XpY(nS,nS,nspin)
  double precision              :: XmY(nS,nS,nspin)
  double precision              :: rho(nBas,nBas,nS,nspin)

! Local variables

  integer                       :: ispin
  logical                       :: adiabatic_connection
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: Ec(:,:)
  double precision,allocatable  :: EcAC(:,:)

! Memory allocation

  allocate(Ec(nAC,nspin),EcAC(nAC,nspin))

  if(singlet_manifold) then

    ispin = 1
    Ec(:,ispin)   = 0d0
    EcAC(:,ispin) = 0d0

    write(*,*) '--------------'
    write(*,*) 'Singlet states'
    write(*,*) '--------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(V x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      if(scaled_screening) then

        call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,lambda,e,ERI, &
                             rho(:,:,:,ispin),Ec(iAC,ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
        call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

      end if
  
      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,lambda,e,ERI, &
                       rho(:,:,:,ispin),Ec(iAC,ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

      call Ec_AC(ispin,dRPA,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),EcAC(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,Ec(iAC,ispin),EcAC(iAC,ispin)

    end do

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',0.5d0*dot_product(wAC,EcAC(:,ispin))
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

  end if
 
  if(triplet_manifold) then

    ispin = 2
    Ec(:,ispin)   = 0d0
    EcAC(:,ispin) = 0d0

    write(*,*) '--------------'
    write(*,*) 'Triplet states'
    write(*,*) '--------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(V x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      if(scaled_screening) then

        call linear_response(ispin,dRPA,TDA,.false.,nBas,nC,nO,nV,nR,nS,lambda,e,ERI, &
                             rho(:,:,:,ispin),Ec(iAC,ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
        call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

      end if  

      call linear_response(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS,lambda,e,ERI, &
                           rho(:,:,:,ispin),Ec(iAC,ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

      call Ec_AC(ispin,dRPA,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),EcAC(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,Ec(iAC,ispin),EcAC(iAC,ispin)

    end do

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',0.5d0*dot_product(wAC,EcAC(:,ispin))
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine ACDFT
