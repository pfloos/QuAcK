subroutine ACFDT_Tmatrix(exchange_kernel,doXBS,dRPA,TDA_T,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS, &
                         ERI,eT,eGT,EcAC)

! Compute the correlation energy via the adiabatic connection fluctuation dissipation theorem for the T-matrix

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doXBS
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: dRPA
  logical,intent(in)            :: TDA_T
  logical,intent(in)            :: TDA
  logical,intent(in)            :: BSE
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: eT(nBas)
  double precision,intent(in)   :: eGT(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin
  integer                       :: isp_T
  integer                       :: iblock
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: Ec(:,:)

  integer                       :: nOOs,nOOt
  integer                       :: nVVs,nVVt

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: TA(:,:)
  double precision,allocatable  :: TB(:,:)
  double precision,allocatable  :: Omega(:,:)
  double precision,allocatable  :: XpY(:,:,:)
  double precision,allocatable  :: XmY(:,:,:)

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

  nOOs = nO*nO
  nVVs = nV*nV

  nOOt = nO*(nO-1)/2
  nVVt = nV*(nV-1)/2

! Memory allocation

  allocate(Omega1s(nVVs),X1s(nVVs,nVVs),Y1s(nOOs,nVVs), &
           Omega2s(nOOs),X2s(nVVs,nOOs),Y2s(nOOs,nOOs), &
           rho1s(nBas,nBas,nVVs),rho2s(nBas,nBas,nOOs), &
           Omega1t(nVVt),X1t(nVVt,nVVt),Y1t(nOOt,nVVt), &
           Omega2t(nOOt),X2t(nVVt,nOOt),Y2t(nOOt,nOOt), &
           rho1t(nBas,nBas,nVVt),rho2t(nBas,nBas,nOOt))
  allocate(TA(nS,nS),TB(nS,nS),Omega(nS,nspin),XpY(nS,nS,nspin),XmY(nS,nS,nspin))
  allocate(Ec(nAC,nspin))

! Antisymmetrized kernel version

  if(exchange_kernel) then

    write(*,*)
    write(*,*) '*** Exchange kernel version ***'
    write(*,*)

  end if

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

      ! Initialize T matrix

      TA(:,:) = 0d0
      TB(:,:) = 0d0

      if(doXBS) then

        isp_T  = 1
        iblock = 3

        call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOs,nVVs,lambda,eT,ERI,  &
                                Omega1s,X1s,Y1s,Omega2s,X2s,Y2s,EcRPA(isp_T))

        call excitation_density_Tmatrix(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

                     call static_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,lambda,ERI,Omega1s,rho1s,Omega2s,rho2s,TA)
        if(.not.TDA) call static_Tmatrix_B(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,lambda,ERI,Omega1s,rho1s,Omega2s,rho2s,TB)

        isp_T  = 2
        iblock = 4

        call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOt,nVVt,lambda,eT,ERI,  &
                                Omega1t,X1t,Y1t,Omega2t,X2t,Y2t,EcRPA(isp_T))

        call excitation_density_Tmatrix(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

                     call static_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,lambda,ERI,Omega1t,rho1t,Omega2t,rho2t,TA)
        if(.not.TDA) call static_Tmatrix_B(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,lambda,ERI,Omega1t,rho1t,Omega2t,rho2t,TB)

      end if

      call linear_response_BSE(ispin,.false.,TDA,eta,nBas,nC,nO,nV,nR,nS,lambda,eGT,ERI,TA,TB, &
                               EcAC(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

      call ACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    if(exchange_kernel) EcAC(ispin) = 0.5d0*EcAC(ispin)

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

      TA(:,:) = 0d0
      TB(:,:) = 0d0

      if(doXBS) then

        isp_T  = 1
        iblock = 3

        call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOs,nVVs,lambda,eT,ERI,  &
                                Omega1s,X1s,Y1s,Omega2s,X2s,Y2s,EcRPA(isp_T))

        call excitation_density_Tmatrix(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

                     call static_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,lambda,ERI,Omega1s,rho1s,Omega2s,rho2s,TA)
        if(.not.TDA) call static_Tmatrix_B(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,lambda,ERI,Omega1s,rho1s,Omega2s,rho2s,TB)

        isp_T  = 2
        iblock = 4

        call linear_response_pp(iblock,TDA_T,nBas,nC,nO,nV,nR,nOOt,nVVt,lambda,eT,ERI,  &
                                Omega1t,X1t,Y1t,Omega2t,X2t,Y2t,EcRPA(isp_T))

        call excitation_density_Tmatrix(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

                     call static_Tmatrix_A(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,lambda,ERI,Omega1t,rho1t,Omega2t,rho2t,TA)
        if(.not.TDA) call static_Tmatrix_B(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,lambda,ERI,Omega1t,rho1t,Omega2t,rho2t,TB)

      end if  

      call linear_response_BSE(ispin,.false.,TDA,eta,nBas,nC,nO,nV,nR,nS,lambda,eGT,ERI,TA,TB, &
                               EcAC(ispin),Omega(:,ispin),XpY(:,:,ispin),XmY(:,:,ispin))

      call ACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,ERI,XpY(:,:,ispin),XmY(:,:,ispin),Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    if(exchange_kernel) EcAC(ispin) = 1.5d0*EcAC(ispin)

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine ACFDT_Tmatrix
