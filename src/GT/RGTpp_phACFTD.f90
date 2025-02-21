subroutine RGTpp_phACFDT(exchange_kernel,doXBS,dRPA,TDA_T,TDA,BSE,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS, &
                        nOOs,nVVs,nOOt,nVVt,Om1s,X1s,Y1s,Om2s,X2s,Y2s,rho1s,rho2s,Om1t,X1t,Y1t,           & 
                        Om2t,X2t,Y2t,rho1t,rho2t,ERI,eT,eGT,EcAC)

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

  integer,intent(in)            :: nOOs
  integer,intent(in)            :: nOOt
  integer,intent(in)            :: nVVs
  integer,intent(in)            :: nVVt

  double precision,intent(in)   :: Om1s(nVVs)
  double precision,intent(in)   :: X1s(nVVs,nVVs)
  double precision,intent(in)   :: Y1s(nOOs,nVVs)
  double precision,intent(in)   :: Om2s(nOOs)
  double precision,intent(in)   :: X2s(nVVs,nOOs)
  double precision,intent(in)   :: Y2s(nOOs,nOOs)
  double precision,intent(in)   :: rho1s(nBas,nBas,nVVs)
  double precision,intent(in)   :: rho2s(nBas,nBas,nOOs)
  double precision,intent(in)   :: Om1t(nVVt)
  double precision,intent(in)   :: X1t(nVVt,nVVt)
  double precision,intent(in)   :: Y1t(nOOt,nVVt)
  double precision,intent(in)   :: Om2t(nOOt)
  double precision,intent(in)   :: X2t(nVVt,nOOt)
  double precision,intent(in)   :: Y2t(nOOt,nOOt)
  double precision,intent(in)   :: rho1t(nBas,nBas,nVVt)
  double precision,intent(in)   :: rho2t(nBas,nBas,nOOt)

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

  double precision              :: EcRPA(nspin)
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: Bpp(:,:)
  double precision,allocatable  :: Cpp(:,:)
  double precision,allocatable  :: Dpp(:,:)
  double precision,allocatable  :: TAs(:,:)
  double precision,allocatable  :: TBs(:,:)
  double precision,allocatable  :: TAt(:,:)
  double precision,allocatable  :: TBt(:,:)
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)

! Output variables

  double precision,intent(out)  :: EcAC(nspin)

! Memory allocation

  allocate(Aph(nS,nS),Bph(nS,nS),TAs(nS,nS),TBs(nS,nS),TAt(nS,nS),TBt(nS,nS),Om(nS),XpY(nS,nS),XmY(nS,nS))
  allocate(Ec(nAC,nspin))

! Hello World

  write(*,*) '-------------------------------------------------------------'
  write(*,*) ' Adiabatic connection version of BSE@GTpp correlation energy '
  write(*,*) '-------------------------------------------------------------'
  write(*,*)

! eXtended BSE

  if(doXBS) then

    write(*,*) '*** scaled screening version (XBS) ***'
    write(*,*)

  end if

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

      if(doXBS) then

        isp_T  = 1
        iblock = 3

        allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

        if(.not.TDA_T) call ppRLR_B(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,lambda,ERI,Bpp)
                       call ppRLR_C(iblock,nBas,nC,nO,nV,nR,nVVs,lambda,eT,ERI,Cpp)
                       call ppRLR_D(iblock,nBas,nC,nO,nV,nR,nOOs,lambda,eT,ERI,Dpp)

        call ppRLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(isp_T))

        deallocate(Bpp,Cpp,Dpp)

        call RGTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

                     call RGTpp_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,lambda,Om1s,rho1s,Om2s,rho2s,TAs)
        if(.not.TDA) call RGTpp_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,lambda,Om1s,rho1s,Om2s,rho2s,TBs)

        isp_T  = 2
        iblock = 4

        allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))

        if(.not.TDA_T) call ppRLR_B(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,lambda,ERI,Bpp)
                       call ppRLR_C(iblock,nBas,nC,nO,nV,nR,nVVt,lambda,eT,ERI,Cpp)
                       call ppRLR_D(iblock,nBas,nC,nO,nV,nR,nOOt,lambda,eT,ERI,Dpp)

        call ppRLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(isp_T))

        deallocate(Bpp,Cpp,Dpp)

        call RGTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

                     call RGTpp_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,lambda,Om1t,rho1t,Om2t,rho2t,TAt)
        if(.not.TDA) call RGtpp_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,lambda,Om1t,rho1t,Om2t,rho2t,TBt)

      end if

                   call phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,eGT,ERI,Aph)
      if(.not.TDA) call phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,Bph)
 
                   Aph(:,:) = Aph(:,:) - TAt(:,:) - TAs(:,:)
      if(.not.TDA) Bph(:,:) = Bph(:,:) - TBt(:,:) - TBs(:,:)
 
      call phRLR(TDA,nS,Aph,Bph,EcAc(ispin),Om,XpY,XmY)

      call phACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,ERI,XpY,XmY,Ec(iAC,ispin))

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

      if(doXBS) then

        isp_T  = 1
        iblock = 3

        allocate(Bpp(nVVs,nOOs),Cpp(nVVs,nVVs),Dpp(nOOs,nOOs))

        if(.not.TDA_T) call ppRLR_B(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,lambda,ERI,Bpp)
                       call ppRLR_C(iblock,nBas,nC,nO,nV,nR,nVVs,lambda,eT,ERI,Cpp)
                       call ppRLR_D(iblock,nBas,nC,nO,nV,nR,nOOs,lambda,eT,ERI,Dpp)

        call ppRLR(TDA_T,nOOs,nVVs,Bpp,Cpp,Dpp,Om1s,X1s,Y1s,Om2s,X2s,Y2s,EcRPA(isp_T))

        deallocate(Bpp,Cpp,Dpp)

        call RGTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOs,nVVs,ERI,X1s,Y1s,rho1s,X2s,Y2s,rho2s)

                     call RGTpp_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,lambda,Om1s,rho1s,Om2s,rho2s,TAs)
        if(.not.TDA) call RGTpp_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOOs,nVVs,lambda,Om1s,rho1s,Om2s,rho2s,TBs)

        isp_T  = 2
        iblock = 4

        allocate(Bpp(nVVt,nOOt),Cpp(nVVt,nVVt),Dpp(nOOt,nOOt))

        if(.not.TDA_T) call ppRLR_B(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,lambda,ERI,Bpp)
                       call ppRLR_C(iblock,nBas,nC,nO,nV,nR,nVVt,lambda,eT,ERI,Cpp)
                       call ppRLR_D(iblock,nBas,nC,nO,nV,nR,nOOt,lambda,eT,ERI,Dpp)

        call ppRLR(TDA_T,nOOt,nVVt,Bpp,Cpp,Dpp,Om1t,X1t,Y1t,Om2t,X2t,Y2t,EcRPA(isp_T))

        deallocate(Bpp,Cpp,Dpp)

        call RGTpp_excitation_density(iblock,nBas,nC,nO,nV,nR,nOOt,nVVt,ERI,X1t,Y1t,rho1t,X2t,Y2t,rho2t)

                     call RGTpp_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,lambda,Om1t,rho1t,Om2t,rho2t,TAt)
        if(.not.TDA) call RGTpp_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,nOOt,nVVt,lambda,Om1t,rho1t,Om2t,rho2t,TBt)

      end if

                   call phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,eGT,ERI,Aph)
      if(.not.TDA) call phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,Bph)

                   Aph(:,:) = Aph(:,:) + TAs(:,:) - TAt(:,:)
      if(.not.TDA) Bph(:,:) = Bph(:,:) + TBs(:,:) - TBt(:,:)

      call phRLR(TDA,nS,Aph,Bph,EcAc(ispin),Om,XpY,XmY)

      call phACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,ERI,XpY,XmY,Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    if(exchange_kernel) EcAC(ispin) = 1.5d0*EcAC(ispin)

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

  end if

end subroutine 
