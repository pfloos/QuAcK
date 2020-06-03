subroutine BSE2_dynamic_perturbation(ispin,eta,nBas,nC,nO,nV,nR,nS,ERI,eHF,eGF,OmBSE,XpY,XmY)

! Compute dynamical effects via perturbation theory for BSE

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: ispin
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas
  integer,intent(in)            :: nC
  integer,intent(in)            :: nO
  integer,intent(in)            :: nV
  integer,intent(in)            :: nR
  integer,intent(in)            :: nS

  double precision,intent(in)   :: ERI(nBas,nBas,nBas,nBas)
  double precision,intent(in)   :: eHF(nBas)
  double precision,intent(in)   :: eGF(nBas)
  double precision,intent(in)   :: OmBSE(nS)
  double precision,intent(in)   :: XpY(nS,nS)
  double precision,intent(in)   :: XmY(nS,nS)

! Local variables

  logical                       :: dTDA = .true.
  integer                       :: ia
  integer,parameter             :: maxS = 10
  double precision              :: gapGF

  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: ZDyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)

  double precision,allocatable  ::  A_dyn(:,:)
  double precision,allocatable  :: ZA_dyn(:,:)

  double precision,allocatable  ::  B_dyn(:,:)
  double precision,allocatable  :: ZB_dyn(:,:)

! Memory allocation

  allocate(OmDyn(nS),ZDyn(nS),X(nS),Y(nS),A_dyn(nS,nS),ZA_dyn(nS,nS))

  if(.not.dTDA) allocate(B_dyn(nS,nS),ZB_dyn(nS,nS))

  ! Print main components of transition vectors

  call print_transition_vectors(nBas,nC,nO,nV,nR,nS,OmBSE,XpY,XmY)

  gapGF = eGF(nO+1) - eGF(nO) 

  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) ' First-order dynamical correction to static 2nd-order Bethe-Salpeter excitation energies '
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(A58,F10.6,A3)') ' BSE neutral excitation must be lower than the GF2 gap = ',gapGF*HaToeV,' eV'
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A20,1X,A20,1X,A20,1X,A20)') '#','Static (eV)','Dynamic (eV)','Correction (eV)','Renorm. (eV)'
  write(*,*) '---------------------------------------------------------------------------------------------------'

  do ia=1,min(nS,maxS)

    X(:) = 0.5d0*(XpY(ia,:) + XmY(ia,:))
    Y(:) = 0.5d0*(XpY(ia,:) - XmY(ia,:))

    ! Resonant part of the BSE correction for dynamical TDA

    call BSE2_A_matrix_dynamic(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0, & 
                               ERI(:,:,:,:),eHF(:),eGF(:),OmBSE(ia),A_dyn(:,:),ZA_dyn(:,:))

    if(dTDA) then 

      ZDyn(ia)  = dot_product(X(:),matmul(ZA_dyn(:,:),X(:)))
      OmDyn(ia) = dot_product(X(:),matmul(A_dyn(:,:),X(:)))

    else

      ! Anti-resonant part of the BSE correction

      call BSE2_B_matrix_dynamic(ispin,eta,nBas,nC,nO,nV,nR,nS,1d0, & 
                                 ERI(:,:,:,:),eHF(:),eGF(:),OmBSE(ia),B_dyn(:,:),ZB_dyn(:,:))

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

  end do
  write(*,*) '---------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine BSE2_dynamic_perturbation
