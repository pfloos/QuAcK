subroutine ACFDT(exchange_kernel,doXBS,dRPA,TDA,BSE,singlet_manifold,triplet_manifold,eta, & 
                 nBas,nC,nO,nV,nR,nS,ERI,e,Omega,XpY,XmY,rho,EcAC)

! Compute the correlation energy via the adiabatic connection fluctuation dissipation theorem

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doXBS
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: dRPA
  logical,intent(in)            :: TDA
  logical,intent(in)            :: BSE
  logical,intent(in)            :: singlet_manifold
  logical,intent(in)            :: triplet_manifold

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS
  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

  double precision              :: Omega(nS,nspin)
  double precision              :: XpY(nS,nS,nspin)
  double precision              :: XmY(nS,nS,nspin)
  double precision              :: rho(nBas,nBas,nS,nspin)

! Local variables

  integer                       :: ispin
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: Ec(:,:)

! Output variables

  double precision,intent(out)  :: EcAC(nspin)

! Memory allocation

  allocate(Ec(nAC,nspin))

! Antisymmetrized kernel version

  if(exchange_kernel) then

    write(*,*)
    write(*,*) '*** Exchange kernel version ***'
    write(*,*)

  end if

! Singlet manifold

  if(singlet_manifold) then

    ispin = 1
    EcAC(ispin) = 0d0
    Ec(:,ispin) = 0d0

    write(*,*) '--------------'
    write(*,*) 'Singlet states'
    write(*,*) '--------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(K x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      if(doXBS) then

        call linear_response(ispin,dRPA,TDA,.false.,eta,nBas,nC,nO,nV,nR,nS,lambda,e,ERI, &
                             rho(:,:,:,ispin),EcAC(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
        call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

      end if

      call linear_response(ispin,dRPA,TDA,BSE,eta,nBas,nC,nO,nV,nR,nS,lambda,e,ERI, &
                       rho(:,:,:,ispin),EcAC(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

      call ACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

  end if
 
! Triplet manifold

  if(triplet_manifold) then

    ispin = 2
    EcAC(ispin) = 0d0
    Ec(:,ispin) = 0d0

    write(*,*) '--------------'
    write(*,*) 'Triplet states'
    write(*,*) '--------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(K x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

!     if(doXBS) then

!       call linear_response(ispin,dRPA,TDA,.false.,eta,nBas,nC,nO,nV,nR,nS,lambda,e,ERI, &
!                            rho(:,:,:,ispin),EcAC(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))
!       call excitation_density(nBas,nC,nO,nR,nS,ERI,XpY(:,:,ispin),rho(:,:,:,ispin))

!     end if  

!     call linear_response(ispin,dRPA,TDA,BSE,eta,nBas,nC,nO,nV,nR,nS,lambda,e,ERI, &
!                          rho(:,:,:,ispin),EcAC(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

!     call ACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),Ec(iAC,ispin))

!     write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine ACFDT
