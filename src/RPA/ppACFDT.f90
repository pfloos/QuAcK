subroutine ppACFDT(TDA,singlet,triplet,nBas,nC,nO,nV,nR,ERI,e,EcAC)

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

  double precision,intent(in)   :: e(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: Ec(:,:)

  integer                       :: nOO
  integer                       :: nVV

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

! Output variables

  double precision,intent(out)  :: EcAC(nspin)

! Hello world

  write(*,*) '--------------------------------------------------------'
  write(*,*) 'Adiabatic connection version of ppRPA correlation energy'
  write(*,*) '--------------------------------------------------------'
  write(*,*)

! Memory allocation

  allocate(Ec(nAC,nspin))

! Antisymmetrized kernel version

  EcAC(:) = 0d0
  Ec(:,:) = 0d0

! Singlet manifold

  if(singlet) then

    ispin = 1

    nOO = nO*(nO+1)/2
    nVV = nV*(nV+1)/2

    allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),rho1(nBas,nBas,nVV), & 
             Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO),rho2(nBas,nBas,nOO), &
             Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))

    write(*,*) '--------------'
    write(*,*) 'Singlet states'
    write(*,*) '--------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(K x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      if(.not.TDA) call ppLR_B(ispin,nBas,nC,nO,nV,nR,nOO,nVV,lambda,ERI,Bpp)
                   call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVV,lambda,e,ERI,Cpp)
                   call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOO,lambda,e,ERI,Dpp)

      call ppLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcAC(ispin))

      call ppACFDT_correlation_energy(ispin,nBas,nC,nO,nV,nR,ERI,nOO,nVV,X1,Y1,X2,Y2,Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

    deallocate(Om1,X1,Y1,rho1,Om2,X2,Y2,rho2,Bpp,Cpp,Dpp)

  end if
 
! Triplet manifold

  if(triplet) then

    ispin = 2

    nOO = nO*(nO-1)/2
    nVV = nV*(nV-1)/2

    allocate(Om1(nVV),X1(nVV,nVV),Y1(nOO,nVV),rho1(nBas,nBas,nVV), & 
             Om2(nOO),X2(nVV,nOO),Y2(nOO,nOO),rho2(nBas,nBas,nOO), &
             Bpp(nVV,nOO),Cpp(nVV,nVV),Dpp(nOO,nOO))

    write(*,*) '--------------'
    write(*,*) 'Triplet states'
    write(*,*) '--------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(K x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      if(.not.TDA) call ppLR_B(ispin,nBas,nC,nO,nV,nR,nOO,nVV,lambda,ERI,Bpp)
                   call ppLR_C(ispin,nBas,nC,nO,nV,nR,nVV,lambda,e,ERI,Cpp)
                   call ppLR_D(ispin,nBas,nC,nO,nV,nR,nOO,lambda,e,ERI,Dpp)

      call ppLR(TDA,nOO,nVV,Bpp,Cpp,Dpp,Om1,X1,Y1,Om2,X2,Y2,EcAC(ispin))

      call ppACFDT_correlation_energy(ispin,nBas,nC,nO,nV,nR,ERI,nOO,nVV,X1,Y1,X2,Y2,Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 1.5d0*dot_product(wAC,Ec(:,ispin))

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

    deallocate(Om1,X1,Y1,rho1,Om2,X2,Y2,rho2,Bpp,Cpp,Dpp)

  end if

end subroutine 
