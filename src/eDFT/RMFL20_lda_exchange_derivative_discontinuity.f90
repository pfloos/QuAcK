subroutine RMFL20_lda_exchange_derivative_discontinuity(nEns,wEns,nGrid,weight,rhow,ExDD)

! Compute the restricted version of the eLDA exchange part of the derivative discontinuity

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
  double precision              :: Cx(nEns)
  double precision              :: dExdw(nEns)
  double precision,external     :: Kronecker_delta

! Output variables

  double precision,intent(out)  :: ExDD(nEns)

! Weight-dependent Cx coefficient for RMFL20 exchange functional

  Cx(1) = Cx0
  Cx(2) = Cx1

! Compute correlation energy for ground, singly-excited and doubly-excited states

  do iEns=1,nEns

    call restricted_elda_exchange_energy(nEns,Cx(iEns),nGrid,weight(:),rhow(:),dExdw(iEns))

  end do

  ExDD(:) = 0d0

  do iEns=1,nEns
    do jEns=1,nEns

      ExDD(iEns) = ExDD(iEns) + (Kronecker_delta(iEns,jEns) - wEns(jEns))*(dExdw(jEns) - dExdw(1))

    end do
  end do

end subroutine RMFL20_lda_exchange_derivative_discontinuity
