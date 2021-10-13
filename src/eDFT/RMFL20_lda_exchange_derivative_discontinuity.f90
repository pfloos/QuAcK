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
  integer                       :: iG
  double precision              :: r
  double precision              :: dExdw(nEns)
  double precision,external     :: Kronecker_delta

  double precision,parameter    :: Cx0   = - (4d0/3d0)*(1d0/pi)**(1d0/3d0)
  double precision,parameter    :: Cx1   = - (176d0/105d0)*(1d0/pi)**(1d0/3d0)

! Output variables

  double precision,intent(out)  :: ExDD(nEns)

! Compute correlation energy for ground- and doubly-excited states

  dExdw(:) = 0d0

  do iG=1,nGrid
    
    r = max(0d0,rhow(iG))
    
    if(r > threshold) then
 
      dExdw(1) = 0d0
      dExdw(2) = dExdw(2) + weight(iG)*(Cx1 - Cx0)*r**(4d0/3d0)

    end if
     
  end do 

  ExDD(:) = 0d0

  do iEns=1,nEns
    do jEns=2,nEns

      ExDD(iEns) = ExDD(iEns) + (Kronecker_delta(iEns,jEns) - wEns(jEns))*dExdw(jEns)

    end do
  end do

end subroutine RMFL20_lda_exchange_derivative_discontinuity
