subroutine Bethe_Salpeter_dynamic_perturbation(TDA,eta,nBas,nC,nO,nV,nR,nS,eGW,OmRPA,OmBSE,XpY,XmY,rho)

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

  double precision,intent(in)   :: eGW(nBas)
  double precision,intent(in)   :: OmRPA(nS)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)
  double precision,intent(in)   :: rho(nBas,nBas,nS)

! Local variables

  logical                       :: TDA_dyn = .false.
  integer                       :: ia
  integer,parameter             :: maxS = 10
  double precision              :: gapGW

  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: ZDyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)
  double precision,allocatable  :: A_dyn(:,:)
  double precision,allocatable  :: B_dyn(:,:)
  double precision,allocatable  :: ZA_dyn(:,:)
  double precision,allocatable  :: ZB_dyn(:,:)

! Memory allocation

  allocate(OmDyn(nS),ZDyn(nS),X(nS),Y(nS),A_dyn(nS,nS),ZA_dyn(nS,nS))

  if(.not.TDA_dyn) allocate(B_dyn(nS,nS),ZB_dyn(nS,nS))

  gapGW = eGW(nO+1) - eGW(nO) 

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static Bethe-Salpeter excitation energies                                 '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do ia=1,min(nS,maxS)

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    ! Resonant part of the BSE correction

    call Bethe_Salpeter_A_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW(:),OmRPA(:),OmBSE(ia),rho(:,:,:),A_dyn(:,:))

    ! Renormalization factor of the resonant part

    call Bethe_Salpeter_ZA_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW(:),OmRPA(:),OmBSE(ia),rho(:,:,:),ZA_dyn(:,:))

    ! First-order correction 

    if(TDA_dyn) then 

      ZDyn(ia)  = dot_product(X(:),matmul(ZA_dyn(:,:),X(:)))
      OmDyn(ia) = dot_product(X(:),matmul(A_dyn(:,:),X(:)))

    else

      ! Anti-resonant part of the BSE correction

      call Bethe_Salpeter_B_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW(:),OmRPA(:),OmBSE(ia),rho(:,:,:),B_dyn(:,:))

      ! Renormalization factor of the anti-resonant part

      call Bethe_Salpeter_ZB_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW(:),OmRPA(:),OmBSE(ia),rho(:,:,:),ZB_dyn(:,:))

      ZDyn(ia)  = dot_product(X(:),matmul(ZA_dyn(:,:),X(:))) &
                - dot_product(Y(:),matmul(ZA_dyn(:,:),Y(:))) &
                + dot_product(X(:),matmul(ZB_dyn(:,:),Y(:))) & 
                - dot_product(Y(:),matmul(ZB_dyn(:,:),X(:)))  

      OmDyn(ia) = dot_product(X(:),matmul(A_dyn(:,:),X(:))) &
                - dot_product(Y(:),matmul(A_dyn(:,:),Y(:))) &
                + dot_product(X(:),matmul(B_dyn(:,:),Y(:))) & 
                - dot_product(Y(:),matmul(B_dyn(:,:),X(:)))  

    end if

    ZDyn(ia) = 1d0/(1d0 - ZDyn(ia))
    OmDyn(ia) = ZDyn(ia)*OmDyn(ia)

    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
      ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV,OmDyn(ia)*HaToeV,ZDyn(ia)

    if(OmBSE(ia) > gapGW) write(*,*) ' !!! BSE neutral excitation larger than the GW gap !!! '

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine Bethe_Salpeter_dynamic_perturbation
