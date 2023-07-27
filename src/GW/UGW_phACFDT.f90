subroutine UGW_phACFDT(exchange_kernel,doXBS,dRPA,TDA_W,TDA,BSE,spin_conserved,spin_flip,eta, &
                       nBas,nC,nO,nV,nR,nS,ERI_aaaa,ERI_aabb,ERI_bbbb,eW,e,EcAC)

! Compute the correlation energy via the adiabatic connection fluctuation dissipation theorem

  implicit none
  include 'parameters.h'
  include 'quadrature.h'

! Input variables

  logical,intent(in)            :: doXBS
  logical,intent(in)            :: dRPA
  logical,intent(in)            :: TDA_W
  logical,intent(in)            :: TDA
  logical,intent(in)            :: BSE
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
  double precision,intent(in)   :: e(nBas,nspin)
  double precision,intent(in)   :: ERI_aaaa(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_aabb(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: ERI_bbbb(nBas,nBas,nBas,nBas)

! Local variables

  integer                       :: ispin
  integer                       :: isp_W
  integer                       :: iAC
  double precision              :: lambda
  double precision,allocatable  :: Ec(:,:)

  double precision              :: EcRPA
  double precision,allocatable  :: OmRPA(:)
  double precision,allocatable  :: XpY_RPA(:,:)
  double precision,allocatable  :: XmY_RPA(:,:)
  double precision,allocatable  :: rho_RPA(:,:,:,:)

  integer                       :: nS_aa,nS_bb,nS_sc
  double precision,allocatable  :: Om_sc(:)
  double precision,allocatable  :: XpY_sc(:,:)
  double precision,allocatable  :: XmY_sc(:,:)

  integer                       :: nS_ab,nS_ba,nS_sf
  double precision,allocatable  :: Om_sf(:)
  double precision,allocatable  :: XpY_sf(:,:)
  double precision,allocatable  :: XmY_sf(:,:)

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

  call phULR(isp_W,.true.,TDA_W,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,1d0,eW, &
             ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
  call UGW_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

! Spin-conserved manifold

  if(spin_conserved) then

    ispin = 1

    allocate(Om_sc(nS_sc),XpY_sc(nS_sc,nS_sc),XmY_sc(nS_sc,nS_sc))

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

        call phULR(isp_W,.true.,TDA_W,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,lambda,eW, &
                   ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
        call UGW_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

      end if

      call phULR(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,lambda,e, &
                 ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcAC(ispin),Om_sc,XpY_sc,XmY_sc)

      call phUACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,nS_aa,nS_bb,nS_sc, &
                                       ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_sc,XmY_sc,Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    if(exchange_kernel) EcAC(ispin) = 0.5d0*EcAC(ispin)

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

    deallocate(Om_sc,XpY_sc,XmY_sc)

  end if
 
! spin-flip manifold

  if(spin_flip) then

    ispin = 2

    ! Memory allocation

    allocate(Om_sf(nS_sf),XpY_sf(nS_sf,nS_sf),XmY_sf(nS_sf,nS_sf))

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

        call phULR(isp_W,.true.,TDA_W,.false.,nBas,nC,nO,nV,nR,nS_aa,nS_bb,nS_sc,nS_sc,lambda,eW, &
                   ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcRPA,OmRPA,XpY_RPA,XmY_RPA)
        call UGW_excitation_density(nBas,nC,nO,nR,nS_aa,nS_bb,nS_sc,ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_RPA,rho_RPA)

      end if

      call phULR(ispin,dRPA,TDA,BSE,nBas,nC,nO,nV,nR,nS_ab,nS_ba,nS_sf,nS_sc,lambda,e, &
                 ERI_aaaa,ERI_aabb,ERI_bbbb,OmRPA,rho_RPA,EcAC(ispin),Om_sf,XpY_sf,XmY_sf)

      call phUACFDT_correlation_energy(ispin,exchange_kernel,nBas,nC,nO,nV,nR,nS,nS_ab,nS_ba,nS_sf, &
                                       ERI_aaaa,ERI_aabb,ERI_bbbb,XpY_sf,XmY_sf,Ec(iAC,ispin))

      write(*,'(2X,F15.6,1X,F30.15,1X,F30.15)') lambda,EcAC(ispin),Ec(iAC,ispin)

    end do

    EcAC(ispin) = 0.5d0*dot_product(wAC,Ec(:,ispin))

    if(exchange_kernel) EcAC(ispin) = 0.5d0*EcAC(ispin)

    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,'(2X,A50,1X,F15.6)') ' Ec(AC) via Gauss-Legendre quadrature:',EcAC(ispin)
    write(*,*) '-----------------------------------------------------------------------------------'
    write(*,*)

    deallocate(Om_sf,XpY_sf,XmY_sf)

  end if

end subroutine 
