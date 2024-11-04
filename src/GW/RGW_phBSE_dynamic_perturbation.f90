subroutine RGW_phBSE_dynamic_perturbation(dophBSE2,dTDA,eta,nBas,nC,nO,nV,nR,nS,eW,eGW,ERI,dipole_int, & 
                                         OmRPA,rho_RPA,OmBSE,XpY,XmY,KA_sta,KB_sta)

! Compute dynamical effects via perturbation theory for BSE

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: dophBSE2
  logical,intent(in)            :: dTDA 
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eW(nBas)
  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: dipole_int(nBas,nBas,ncart)
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: rho_RPA(nBas,nBas,nS)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)
  double precision,intent(in)   :: KA_sta(nS,nS)
  double precision,intent(in)   :: KB_sta(nS,nS)

! Local variables

  integer                       :: ia

  integer                       :: maxS = 10

  double precision,allocatable  :: Om_dyn(:)
  double precision,allocatable  :: Z_dyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  :: KAp_dyn(:,:)
  double precision,allocatable  :: KAm_dyn(:,:)
  double precision,allocatable  :: ZAp_dyn(:,:)
  double precision,allocatable  :: ZAm_dyn(:,:)

  double precision,allocatable  :: KB_dyn(:,:)

  double precision,allocatable  :: W(:,:,:,:)

! Memory allocation

  allocate(Om_dyn(maxS),Z_dyn(maxS),X(nS),Y(nS),KAp_dyn(nS,nS),ZAp_dyn(nS,nS), &
           KAm_dyn(nS,nS),ZAm_dyn(nS,nS),KB_dyn(nS,nS))

  if(dTDA) then 
    write(*,*)
    write(*,*) '*** dynamical TDA activated ***'
    write(*,*)
  end if

  if(dophBSE2) then 

    write(*,*)
    write(*,*) '*** Second-order BSE dynamic kernel activated! ***'
    write(*,*)

    allocate(W(nBas,nBas,nBas,nBas))
    call RGW_phBSE_static_kernel(eta,nBas,nC,nO,nV,nR,nS,1d0,ERI,OmRPA,rho_RPA,W)

  end if

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static Bethe-Salpeter excitation energies                     '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do ia=1,min(nS,maxS)

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    ! Resonant part of the BSE correction for dynamical TDA

    call RGW_phBSE_dynamic_kernel_A(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,OmRPA,rho_RPA,+OmBSE(ia),KAp_dyn,ZAp_dyn)

    if(dophBSE2) call RGW_phBSE2_dynamic_kernel_A(eta,nBas,nC,nO,nV,nR,nS,eGW,W,OmBSE(ia),KAp_dyn,ZAp_dyn)

    if(dTDA) then 

      Z_dyn(ia)  = dot_product(X,matmul(ZAp_dyn,X))
      Om_dyn(ia) = dot_product(X,matmul(KAp_dyn - KA_sta,X))

    else

      ! Resonant and anti-resonant part of the BSE correction

      call RGW_phBSE_dynamic_kernel_A(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,OmRPA,rho_RPA,-OmBSE(ia),KAm_dyn,ZAm_dyn)
      call RGW_phBSE_dynamic_kernel_B(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW,OmRPA,rho_RPA,KB_dyn)

      if(dophBSE2) call RGW_phBSE2_dynamic_kernel_B(eta,nBas,nC,nO,nV,nR,nS,eGW,W,KB_dyn)

      ! Renormalization factor of the resonant and anti-resonant parts

      Z_dyn(ia)  = dot_product(X,matmul(ZAp_dyn,X)) &
                 + dot_product(Y,matmul(ZAm_dyn,Y))

      Om_dyn(ia) = dot_product(X,matmul(KAp_dyn - KA_sta,X)) &
                 - dot_product(Y,matmul(KAm_dyn - KA_sta,Y)) &
                 + dot_product(X,matmul(KB_dyn  - KB_sta,Y)) &
                 - dot_product(Y,matmul(KB_dyn  - KB_sta,X))

    end if

    Z_dyn(ia)  = 1d0/(1d0 - Z_dyn(ia))
    Om_dyn(ia) = Z_dyn(ia)*Om_dyn(ia)

    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
      ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+Om_dyn(ia))*HaToeV,Om_dyn(ia)*HaToeV,Z_dyn(ia)

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine 
