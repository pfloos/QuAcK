subroutine UGW_phACFDT(exchange_kernel,doXBS,TDA_W,TDA,spin_conserved,spin_flip,eta, &
                       nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,eW,eGW,EcAC)

! Compute the correlation energy via the adiabatic connection fluctuation dissipation theorem

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doXBS
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: exchange_kernel
  logical,intent(in)            :: spin_conserved
  logical,intent(in)            :: spin_flip

  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC(nspin)
  integer,intent(in)            :: nO(nspin)
  integer,intent(in)            :: nV(nspin)
  integer,intent(in)            :: nR(nspin)
  integer,intent(in)            :: nS(nspin)
  double precision,intent(in)   :: eW(nBas,nspin)
  double precision,intent(in)   :: eGW(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)

! Local variables

  logical                       :: dRPA   = .false.
  logical                       :: dRPA_W = .true.

  integer                       :: ispin
  integer                       :: isp_W
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: Ec(:,:)

  double precision,allocatable  :: Aph(:,:)
  double precision,allocatable  :: Bph(:,:)
  double precision,allocatable  :: KA(:,:)
  double precision,allocatable  :: KB(:,:)

  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:,:)

  integer                       :: nS_aa,nS_bb,nS_sc
  integer                       :: nS_ab,nS_ba,nS_sf
  double precision,allocatable  :: Om(:)
  double precision,allocatable  :: XpY(:,:)
  double precision,allocatable  :: XmY(:,:)

! Output variables

  double precision,intent(out)  :: EcAC(nspin)

! Memory allocation

  allocate(Ec(nAC,nspin))

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
  EcRPA = 0d0

  ! Memory allocation

  nS_aa = nS(1)
  nS_bb = nS(2)
  nS_sc = nS_aa + nS_bb

  nS_ab = (nO(1) - nC(1))*(nV(2) - nR(2))
  nS_ba = (nO(2) - nC(2))*(nV(1) - nR(1))
  nS_sf = nS_ab + nS_ba

  allocate(OmRPA(nS_sc),XpY_RPA(nS_sc,nS_sc),XmY_RPA(nS_sc,nS_sc),rho_RPA(nBas,nBas,nS_sc,nspin))
  allocate(Aph(nS_sc,nS_sc),Bph(nS_sc,nS_sc),KA(nS_sc,nS_sc),KB(nS_sc,nS_sc))

               call phULR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,1d0,eW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
  if(.not.TDA) call phULR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,1d0,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

  call phULR(TDA_W,nS_aa,nS_bb,nS_sc,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call UGW_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

  call UGW_phBSE_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, & 
                                 ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,KA)
  call UGW_phBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0, & 
                                 ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,KB)

! Spin-conserved manifold

  if(spin_conserved) then

    ispin = 1

    allocate(Om(nS_sc),XpY(nS_sc,nS_sc),XmY(nS_sc,nS_sc))

    write(*,*) '------------------------'
    write(*,*) 'Spin-conserved manifold '
    write(*,*) '------------------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(K x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      if(doXBS) then

                     call phULR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,lambda,eW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
        if(.not.TDA) call phULR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)
   
        call phULR(TDA_W,nS_aa,nS_bb,nS_sc,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
        call UGW_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

        call UGW_phBSE_static_kernel_A(ispin,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,lambda, & 
                                       ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,KA)
        call UGW_phBSE_static_kernel_B(ispin,eta,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,lambda, & 
                                       ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,KB)

      end if

                   call phULR_A(isp_W,dRPA,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,lambda,eW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
      if(.not.TDA) call phULR_B(isp_W,dRPA,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

                   Aph(:,:) = Aph(:,:) + KA(:,:)
      if(.not.TDA) Bph(:,:) = Bph(:,:) + KB(:,:)

      call phULR(TDA_W,nS_aa,nS_bb,nS_sc,Aph,Bph,EcAC(ispin),OmRPA,XpY_RPA,XmY_RPA)

      call phUACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,nS_aa,nS_bb,nS_sc, &
                                       ERI_aaaa,ERI_aabb,ERI_bbbb,XpY,XmY,Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    if(exchange_kernel) EcAC(ispin) = 0.5d0*EcAC(ispin)

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

    deallocate(Om,XpY,XmY)

  end if
 
! spin-flip manifold

  if(spin_flip) then

    ispin = 2

    ! Memory allocation

    allocate(Om(nS_sf),XpY(nS_sf,nS_sf),XmY(nS_sf,nS_sf))

    write(*,*) '--------------------'
    write(*,*) ' Spin-flip manifold '
    write(*,*) '--------------------'
    write(*,*) 

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A15,1X,A30,1X,A30)') 'lambda','Ec(lambda)','Tr(K x P_lambda)'
    write(*,*) '-----------------------------------------------------------------------------------'

    do iAC=1,nAC
 
      lambda = rAC(iAC)

      if(doXBS) then

                     call phULR_A(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,lambda,eW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
        if(.not.TDA) call phULR_B(isp_W,dRPA_W,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

        call phULR(TDA_W,nS_aa,nS_bb,nS_sc,Aph,Bph,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
        call UGW_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

      end if

                   call phULR_A(isp_W,dRPA,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sf,lambda,eW,ERI_aaaa,ERI_aabb,ERI_bbbb,Aph)
      if(.not.TDA) call phULR_B(isp_W,dRPA,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sf,lambda,ERI_aaaa,ERI_aabb,ERI_bbbb,Bph)

                   Aph(:,:) = Aph(:,:) + KA(:,:)
      if(.not.TDA) Bph(:,:) = Bph(:,:) + KB(:,:)

      call phULR(TDA_W,nS_aa,nS_bb,nS_sf,Aph,Bph,EcAC(ispin),OmRPA,XpY_RPA,XmY_RPA)

      call phUACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,nS_ab,nS_ba,nS_sf, &
                                       ERI_aaaa,ERI_aabb,ERI_bbbb,XpY,XmY,Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    if(exchange_kernel) EcAC(ispin) = 0.5d0*EcAC(ispin)

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

    deallocate(Om,XpY,XmY)

  end if

end subroutine 
