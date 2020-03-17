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

  integer                       :: iG
  double precision              :: Cx0
  double precision              :: Cx1
  double precision              :: rw
  double precision              :: dExdw

! Output variables

  double precision,intent(out)  :: ExDD(nEns)

! Weight-dependent Cx coefficient for RMFL20 exchange functional

  Cx0   = -(4d0/3d0)*(1d0/pi)**(1d0/3d0)
  Cx1   = -(176d0/105d0)*(1d0/pi)**(1d0/3d0)

! Compute correlation energy for ground, singly-excited and doubly-excited states

  dExdw = 0d0

  do iG=1,nGrid

    rw = max(0d0,rhow(iG))

    if(rw > threshold) then

      dExdw = dExdw + weight(iG)*(Cx1 - Cx0)*rw**(4d0/3d0)

    end if

  end do

  ExDD(1) =      - wEns(2) *dExdw
  ExDD(2) = (1d0 - wEns(2))*dExdw

end subroutine RMFL20_lda_exchange_derivative_discontinuity
