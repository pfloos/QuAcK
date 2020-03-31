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

! Output variables

  double precision,intent(out)  :: ExDD(nEns)

! Compute correlation energy for ground- and doubly-excited states

  dExdw(:) = 0d0

  do iG=1,nGrid
    
    r = max(0d0,rhow(iG))
    
    if(r > threshold) then
      dExdw(1) = dExdw(1) + weight(iG)*Cx0*r**(4d0/3d0)
      dExdw(2) = dExdw(2) + weight(iG)*Cx1*r**(4d0/3d0)
    end if
     
  end do 

  ExDD(:) = 0d0

  do iEns=1,nEns
    do jEns=1,nEns

      ExDD(iEns) = ExDD(iEns) + (Kronecker_delta(iEns,jEns) - wEns(jEns))*(dExdw(jEns) - dExdw(1))    

    end do
  end do

end subroutine RMFL20_lda_exchange_derivative_discontinuity
