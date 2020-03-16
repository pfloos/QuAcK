subroutine MFL20_lda_correlation_derivative_discontinuity(nEns,wEns,nGrid,weight,rhow,Ec)

! Compute eLDA correlation part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid,nspin)

! Local variables

  integer                       :: iEns,jEns
  double precision,allocatable  :: aMFL(:,:)
  double precision              :: dEc(nsp,nEns)
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: Ec(nsp,nEns)

! Allocation

  allocate(aMFL(3,nEns))

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

    call elda_correlation_energy(nEns,aMFL(:,iEns),nGrid,weight(:),rhow(:,:),dEc(:,iEns))

  end do

  Ec(:,:) = 0d0

  do iEns=1,nEns
    do jEns=1,nEns

      Ec(:,iEns) = Ec(:,iEns) + (Kronecker_delta(iEns,jEns) - wEns(jEns))*(dEc(:,jEns) - dEc(:,1))

    end do
  end do

end subroutine MFL20_lda_correlation_derivative_discontinuity
