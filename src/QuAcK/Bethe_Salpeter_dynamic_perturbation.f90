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

  integer                       :: ia
  integer,parameter             :: maxS = 10

  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: ZDyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)
  double precision,allocatable  :: A_dyn(:,:)
  double precision,allocatable  :: B_dyn(:,:)
  double precision,allocatable  :: Z_dyn(:,:)

! Memory allocation

  allocate(OmDyn(nS),ZDyn(nS),X(nS),Y(nS),A_dyn(nS,nS),Z_dyn(nS,nS))

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static Bethe-Salpeter excitation energies                                 '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'
  do ia=1,min(nS,maxS)

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    call Bethe_Salpeter_A_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW(:),OmRPA(:),OmBSE(ia),rho(:,:,:),A_dyn(:,:))
    call Bethe_Salpeter_Z_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW(:),OmRPA(:),OmBSE(ia),rho(:,:,:),Z_dyn(:,:))

    ! First-order correction 

    if(.true.) then 

      ZDyn(ia) = dot_product(X(:),matmul(Z_dyn(:,:),X(:)))

    else

      allocate(B_dyn(nS,nS))
      call Bethe_Salpeter_B_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,eGW(:),OmRPA(:),OmBSE(ia),rho(:,:,:),B_dyn(:,:))

      OmDyn(ia) = dot_product(X(:),matmul(A_dyn(:,:),X(:))) &
                - dot_product(Y(:),matmul(A_dyn(:,:),Y(:))) &
                + dot_product(X(:),matmul(B_dyn(:,:),Y(:))) & 
                - dot_product(Y(:),matmul(B_dyn(:,:),X(:)))  

    end if

    ! Renormalization factor

    ZDyn(ia) = 1d0/(1d0 - ZDyn(ia))
    OmDyn(ia) = ZDyn(ia)*dot_product(X(:),matmul(A_dyn(:,:),X(:)))

    write(*,'(2X,I5,5X,F15.6,5X,F15.6,5X,F15.6,5X,F15.6)') & 
      ia,OmBSE(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV,OmDyn(ia)*HaToeV,ZDyn(ia)

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine Bethe_Salpeter_dynamic_perturbation
