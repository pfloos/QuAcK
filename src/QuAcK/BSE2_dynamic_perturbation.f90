subroutine BSE2_dynamic_perturbation(TDA,eta,nBas,nC,nO,nV,nR,nS,eHF,eGF,OmBSE,XpY,XmY)

! Compute dynamical effects via perturbation theory for BSE

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  logical                       :: dTDA = .false.
  integer                       :: ia
  integer,parameter             :: maxS = 10
  double precision              :: gapGF

  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: ZDyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  ::  Ap_dyn(:,:)
  double precision,allocatable  ::  Am_dyn(:,:)
  double precision,allocatable  :: ZAp_dyn(:,:)
  double precision,allocatable  :: ZAm_dyn(:,:)

  double precision,allocatable  ::  Bp_dyn(:,:)
  double precision,allocatable  ::  Bm_dyn(:,:)
  double precision,allocatable  :: ZBp_dyn(:,:)
  double precision,allocatable  :: ZBm_dyn(:,:)

! Memory allocation

  allocate(OmDyn(nS),ZDyn(nS),X(nS),Y(nS),Ap_dyn(nS,nS),ZAp_dyn(nS,nS))

  if(.not.dTDA) allocate(Am_dyn(nS,nS),ZAm_dyn(nS,nS),Bp_dyn(nS,nS),Bm_dyn(nS,nS),ZBp_dyn(nS,nS),ZBm_dyn(nS,nS))

  gapGF = eGF(nO+1) - eGF(nO) 

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static 2nd-order Bethe-Salpeter excitation energies '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do ia=1,min(nS,maxS)

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    ! First-order correction 

    if(dTDA) then 

      ! Resonant part of the BSE correction for dynamical TDA

      call BSE2_A_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eHF(:),eGF(:),OmBSE(ia),Ap_dyn(:,:))

      ! Renormalization factor of the resonant parts for dynamical TDA

      call BSE2_ZA_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eHF(:),eGF(:),OmBSE(ia),ZAp_dyn(:,:))

      ZDyn(ia)  = dot_product(X(:),matmul(ZAp_dyn(:,:),X(:)))
      OmDyn(ia) = dot_product(X(:),matmul(Ap_dyn(:,:),X(:)))

    else

      ! Resonant and anti-resonant part of the BSE correction

!     call Bethe_Salpeter_AB_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eHF(:),eGF(:),OmBSE(ia), & 
!                                           Ap_dyn(:,:),Am_dyn(:,:),Bp_dyn(:,:),Bm_dyn(:,:))

      ! Renormalization factor of the resonant and anti-resonant parts

!     call Bethe_Salpeter_ZAB_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eHF(:),eGF(:),OmBSE(ia), &
!                                            ZAp_dyn(:,:),ZAm_dyn(:,:),ZBp_dyn(:,:),ZBm_dyn(:,:))

      ZDyn(ia)  = dot_product(X(:),matmul(ZAp_dyn(:,:),X(:))) &
                - dot_product(Y(:),matmul(ZAm_dyn(:,:),Y(:))) &
                + dot_product(X(:),matmul(ZBp_dyn(:,:),Y(:))) & 
                - dot_product(Y(:),matmul(ZBm_dyn(:,:),X(:)))  

      OmDyn(ia) = dot_product(X(:),matmul(Ap_dyn(:,:),X(:))) &
                - dot_product(Y(:),matmul(Am_dyn(:,:),Y(:))) &
                + dot_product(X(:),matmul(Bp_dyn(:,:),Y(:))) & 
                - dot_product(Y(:),matmul(Bm_dyn(:,:),X(:)))  

    end if

    ZDyn(ia) = 1d0/(1d0 - ZDyn(ia))
    OmDyn(ia) = ZDyn(ia)*OmDyn(ia)

    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
      ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV,OmDyn(ia)*HaToeV,ZDyn(ia)

    if(OmBSE(ia) > gapGF) write(*,*) ' !!! BSE2 neutral excitation larger than the GF2 gap !!! '

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine BSE2_dynamic_perturbation
