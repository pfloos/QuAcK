subroutine RMFL20_lda_correlation_Levy_Zahariev_shift(nEns,wEns,nGrid,weight,rho,EcLZ)

! Compute the restricted Marut-Fromager-Loos LDA correlation contribution to Levy-Zahariev shift

  implicit none

  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rho(nGrid)

! Local variables

  logical                       :: LDA_centered = .true.
  integer                       :: iEns
  double precision              :: EcLZLDA
  double precision,allocatable  :: aMFL(:,:)
  double precision,allocatable  :: EcLZeLDA(:)

! Output variables

  double precision,intent(out)  :: EcLZ


! Allocation

  allocate(aMFL(3,nEns),EcLZeLDA(nEns))

! Parameters for weight-dependent LDA correlation functional

  aMFL(1,1) = -0.0238184d0
  aMFL(2,1) = +0.00540994d0
  aMFL(3,1) = +0.0830766d0

  aMFL(1,2) = -0.0144633d0
  aMFL(2,2) = -0.0506019d0
  aMFL(3,2) = +0.0331417d0

! Compute correlation energy for ground, singly-excited and doubly-excited states

  do iEns=1,nEns

    call restricted_elda_correlation_Levy_Zahariev_shift(nEns,aMFL(:,iEns),nGrid,weight(:),rho(:),EcLZeLDA(iEns))

  end do

! LDA-centered functional

  EcLZLDA = 0d0

  if(LDA_centered) then

    call RVWN5_lda_correlation_Levy_Zahariev_shift(nGrid,weight(:),rho(:),EcLZLDA)

    do iEns=1,nEns
 
      EcLZeLDA(iEns) = EcLZeLDA(iEns) + EcLZLDA - EcLZeLDA(1)
 
    end do

  end if

! Weight-denpendent functional for ensembles

  EcLZ = 0d0

  do iEns=1,nEns

    EcLZ = EcLZ + wEns(iEns)*EcLZeLDA(iEns)

  enddo

end subroutine RMFL20_lda_correlation_Levy_Zahariev_shift
