subroutine MFL20_lda_correlation_Levy_Zahariev_shift(nEns,wEns,nGrid,weight,rho,EcLZ)

! Compute Marut-Fromager-Loos's LDA contribution to Levy-Zahariev shift

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid,nspin)

! Local variables

  logical                       :: LDA_centered = .false.
  integer                       :: iEns
  double precision              :: EcLZLDA(nsp)
  double precision,allocatable  :: aMFL(:,:)
  double precision,allocatable  :: EcLZeLDA(:,:)

! Output variables

  double precision,intent(out)  :: EcLZ(nsp)


! Allocation

  allocate(aMFL(3,nEns),EcLZeLDA(nsp,nEns))

! Parameters for weight-dependent LDA correlation functional

  aMFL(1,1) = -0.0238184d0
  aMFL(2,1) = +0.00575719d0
  aMFL(3,1) = +0.0830576d0

  aMFL(1,2) = -0.0282814d0
  aMFL(2,2) = +0.00340758d0
  aMFL(3,2) = +0.0663967d0

  aMFL(1,3) = -0.0144633d0
  aMFL(2,3) = -0.0504501d0
  aMFL(3,3) = +0.0331287d0

! Compute correlation energy for ground, singly-excited and doubly-excited states

  do iEns=1,nEns

    call elda_correlation_Levy_Zahariev_shift(nEns,aMFL(:,iEns),nGrid,weight(:),rho(:,:),EcLZeLDA(:,iEns))

  end do

! LDA-centered functional

  EcLZLDA(:) = 0d0

  if(LDA_centered) then

    call VWN5_lda_correlation_Levy_Zahariev_shift(nGrid,weight(:),rho(:,:),EcLZLDA(:))

    do iEns=1,nEns
 
      EcLZeLDA(:,iEns) = EcLZeLDA(:,iEns) + EcLZLDA(:)- EcLZeLDA(:,1)
 
    end do

  end if

! Weight-denpendent functional for ensembles

  EcLZ(:) = 0d0

  do iEns=1,nEns

    EcLZ(:) = EcLZ(:) + wEns(iEns)*EcLZeLDA(:,iEns)

  enddo

end subroutine MFL20_lda_correlation_Levy_Zahariev_shift
