subroutine RGW_phACFDT(exchange_kernel,doXBS,TDA_W,TDA,singlet,triplet,eta,nBas,nC,nO,nV,nR,nS,ERI,eW,eGW,EcAC)

! Compute the correlation energy via the adiabatic connection fluctuation dissipation theorem

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doXBS
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: singlet
  logical,intent(in)            :: triplet

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS
  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dRPA   = .false.
  logical                       :: dRPA_W = .true.
  integer                       :: ispin
  integer                       :: isp_W
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: Ec(:,:)

  double precision              :: EcRPA
  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: KA(:,:)
  double precision,allocatable  :: KB(:,:)
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:)

  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)

! Output variables

  double precision,intent(out)  :: EcAC(nspin)

! Memory allocation

  allocate(Ec(nAC,nspin))
  allocate(Aph(nS,nS),Bph(nS,nS),KA(nS,nS),KB(nS,nS),OmRPA(nS),XpY_RPA(nS,nS),XmY_RPA(nS,nS), & 
           rho_RPA(nBas,nBas,nS),Om(nS),XpY(nS,nS),XmY(nS,nS))

! Hello World

  write(*,*) '-----------------------------------------------------------'
  write(*,*) ' Adiabatic connection version of BSE@GW correlation energy '
  write(*,*) '-----------------------------------------------------------'
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

! Compute (singlet) RPA screening 

  isp_W = 1

                 call phRLR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,1d0,eW,ERI,Aph)
  if(.not.TDA_W) call phRLR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,1d0,ERI,Bph)

  call phRLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call RGW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)

               call RGW_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KA)
  if(.not.TDA) call RGW_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,KB)

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

                       call phRLR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,lambda,eW,ERI,Aph)
        if(.not.TDA_W) call phRLR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,lambda,ERI,Bph)
   
        call phRLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
        call RGW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)
   
                     call RGW_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,OmRPA,rho_RPA,KA)
        if(.not.TDA) call RGW_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,OmRPA,rho_RPA,KB)

      end if

                   call phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,eGW,ERI,Aph)
      if(.not.TDA) call phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,Bph)

                   Aph(:,:) = Aph(:,:) + KA(:,:)
      if(.not.TDA) Bph(:,:) = Bph(:,:) + KB(:,:)

      call phRLR(TDA,nS,Aph,Bph,EcAC(ispin),Om,XpY,XmY)

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

                       call phRLR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,lambda,eW,ERI,Aph)
        if(.not.TDA_W) call phRLR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS,lambda,ERI,Bph)
  
        call phRLR(TDA_W,nS,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
        call RGW_excitation_density(nBas,nC,nO,nR,nS,ERI,XpY_RPA,rho_RPA)
  
                     call RGW_phBSE_static_kernel_A(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,OmRPA,rho_RPA,KA)
        if(.not.TDA) call RGW_phBSE_static_kernel_B(eta,nBas,nC,nO,nV,nR,nS,lambda,ERI,OmRPA,rho_RPA,KB)

      end if  

                   call phRLR_A(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,eGW,ERI,Aph)
      if(.not.TDA) call phRLR_B(ispin,dRPA,nBas,nC,nO,nV,nR,nS,lambda,ERI,Bph)

                   Aph(:,:) = Aph(:,:) + KA(:,:)
      if(.not.TDA) Bph(:,:) = Bph(:,:) + KB(:,:)

      call phRLR(TDA,nS,Aph,Bph,EcAC(ispin),Om,XpY,XmY)

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
