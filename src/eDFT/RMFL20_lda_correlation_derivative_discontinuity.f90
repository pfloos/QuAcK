subroutine RMFL20_lda_correlation_derivative_discontinuity(nEns,wEns,nGrid,weight,rhow,EcDD)

! Compute the restricted version of the eLDA correlation part of the derivative discontinuity

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nEns
  double precision,intent(in)   :: wEns(nEns)
  integer,intent(in)            :: nGrid
  double precision,intent(in)   :: weight(nGrid)
  double precision,intent(in)   :: rhow(nGrid)

! Local variables

  integer                       :: iEns,jEns
  double precision,allocatable  :: aMFL(:,:)
  double precision              :: dEcdw(nEns)
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: EcDD(nEns)

! Allocation

  allocate(aMFL(3,nEns))

! Parameters for weight-dependent LDA correlation functional

  aMFL(1,1) = -0.0238184d0
  aMFL(2,1) = +0.00540994d0
  aMFL(3,1) = +0.0830766d0

  aMFL(1,2) = -0.0144633d0
  aMFL(2,2) = -0.0506019d0
  aMFL(3,2) = +0.0331417d0

! Compute correlation energy for ground, singly-excited and doubly-excited states

  do iEns=1,nEns

    call restricted_elda_correlation_energy(aMFL(:,iEns),nGrid,weight(:),rhow(:),dEcdw(iEns))

  end do

  EcDD(:) = 0d0

  do iEns=1,nEns
    do jEns=1,nEns

      EcDD(iEns) = EcDD(iEns) + (Kronecker_delta(iEns,jEns) - wEns(jEns))*(dEcdw(jEns) - dEcdw(1))

    end do
  end do

end subroutine RMFL20_lda_correlation_derivative_discontinuity
