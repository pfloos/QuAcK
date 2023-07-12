subroutine ACFDT_pp(TDA,singlet,triplet,nBas,nC,nO,nV,nR,nS,ERI,e,EcAC)

! Compute the correlation energy via the adiabatic connection fluctuation dissipation theorem for pp sector

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: TDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: Ec(:,:)

  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt

  double precision,allocatable  :: Omega1s(:),Omega1t(:)
  double precision,allocatable  :: X1s(:,:),X1t(:,:)
  double precision,allocatable  :: Y1s(:,:),Y1t(:,:)
  double precision,allocatable  :: rho1s(:,:,:),rho1t(:,:,:)
  double precision,allocatable  :: Omega2s(:),Omega2t(:)
  double precision,allocatable  :: X2s(:,:),X2t(:,:)
  double precision,allocatable  :: Y2s(:,:),Y2t(:,:)
  double precision,allocatable  :: rho2s(:,:,:),rho2t(:,:,:)

! Output variables

  double precision,intent(out)  :: EcAC(nspin)

! Useful quantities

  nOOs = nO*(nO+1)/2
  nVVs = nV*(nV+1)/2

  nOOt = nO*(nO-1)/2
  nVVt = nV*(nV-1)/2

! Memory allocation

  allocate(Omega1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs), &
           Omega2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs), &
           rho1s(nBas,nBas,nVVs),rho2s(nBas,nBas,nOOs), &
           Omega1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt), &
           Omega2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt), &
           rho1t(nBas,nBas,nVVt),rho2t(nBas,nBas,nOOt))
  allocate(Ec(nAC,nspin))

! Antisymmetrized kernel version

  EcAC(:) = 0d0
  Ec(:,:) = 0d0

! Singlet manifold

  if(singlet) then

    ispin = 1

    write(*,*) '--------------'
    write(*,*) 'Singlet states'
    write(*,*) '--------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(K x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      call ppLR(ispin,TDA,nBas,nC,nO,nV,nR,nOOs,nVVs,lambda,e,ERI,Omega1s,X1s,Y1s,Omega2s,X2s,Y2s,EcAC(ispin))

      call ACFDT_pp_correlation_energy(ispin,nBas,nC,nO,nV,nR,nS,ERI,nOOs,nVVs,X1s,Y1s,X2s,Y2s,Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

  end if
 
! Triplet manifold

  if(triplet) then

    ispin = 2

    write(*,*) '--------------'
    write(*,*) 'Triplet states'
    write(*,*) '--------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(K x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      ! Initialize T matrix

      call ppLR(ispin,TDA,nBas,nC,nO,nV,nR,nOOt,nVVt,lambda,e,ERI,Omega1t,X1t,Y1t,Omega2t,X2t,Y2t,EcAC(ispin))

      call ACFDT_pp_correlation_energy(ispin,nBas,nC,nO,nV,nR,nS,ERI,nOOt,nVVt,X1t,Y1t,X2t,Y2t,Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 1.5d0*dot_product(wAC,Ec(:,ispin))

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine ACFDT_pp
