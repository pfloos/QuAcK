subroutine BSE_dynamic_perturbation(TDA,eta,nBas,nC,nO,nV,nR,nS,OmRPA,OmBSE,XpY,XmY,rho)

! Compute dynamical effects via perturbation theory for BSE

  implicit none
  include 'parameters.h'

! Input variables

  logical,intent(in)            :: TDA
  double precision,intent(in)   :: eta
  integer,intent(in)            :: nBas,nC,nO,nV,nR,nS

  double precision              :: OmRPA(nS)
  double precision              :: OmBSE(nS)
  double precision              :: XpY(nS,nS)
  double precision              :: XmY(nS,nS)
  double precision              :: rho(nBas,nBas,nS)

! Local variables

  integer                       :: ia
  integer,parameter             :: maxS = 10

  double precision,allocatable  :: OmDyn(:)
  double precision,allocatable  :: X(:)
  double precision,allocatable  :: Y(:)
  double precision,allocatable  :: A_dyn(:,:)
  double precision,allocatable  :: B_dyn(:,:)

! Memory allocation

  allocate(OmDyn(nS),X(nS),Y(nS),A_dyn(nS,nS),B_dyn(nS,nS))

  write(*,*) '-----------------------------------------------------------------------------------------------------------'
  write(*,'(2X,A5,1X,A30,1X,A30,1X,A30)') '#','Static excitation (eV)','Dynamic correction (eV)','Dynamic excitation (eV)'
  write(*,*) '-----------------------------------------------------------------------------------------------------------'

  do ia=1,maxS

    X(:) = 0.5d0*(XpY(:,ia) + XmY(:,ia))
    Y(:) = 0.5d0*(XpY(:,ia) - XmY(:,ia))

    call Bethe_Salpeter_A_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,OmRPA(:),OmBSE(:),rho(:,:,:),A_dyn(:,:))

    if(TDA) then 
      B_dyn(:,:) = 0d0
    else
      call Bethe_Salpeter_B_matrix_dynamic(eta,nBas,nC,nO,nV,nR,nS,1d0,OmRPA(:),OmBSE(:),rho(:,:,:),B_dyn(:,:))
    end if

    OmDyn(ia) = dot_product(X(:),matmul(A_dyn(:,:),X(:))) - dot_product(Y(:),matmul(A_dyn(:,:),Y(:))) &
              + dot_product(X(:),matmul(B_dyn(:,:),Y(:))) - dot_product(Y(:),matmul(B_dyn(:,:),X(:)))  

    write(*,'(2X,I5,15X,F15.6,15X,F15.6,15X,F15.6)') ia,OmBSE(ia)*HaToeV,OmDyn(ia)*HaToeV,(OmBSE(ia)+OmDyn(ia))*HaToeV

  end do
  write(*,*) '-----------------------------------------------------------------------------------------------------------'
  write(*,*) 

end subroutine BSE_dynamic_perturbation
