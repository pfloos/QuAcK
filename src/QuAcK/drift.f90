subroutine drift(nWalk,r,r12,g,dg,F)

! Compute quantum force

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: nWalk
  double precision,intent(in)   :: r(nWalk,2,3),r12(nWalk),g(nWalk,2),dg(nWalk,2,3)

! Local variables

  logical                       :: smoothDrift
  double precision              :: rij,rijSq,w,wSq,damp
  integer                       :: iW

! Output variables

  double precision,intent(out)  :: F(nWalk,2,3)

! Compute

  smoothDrift = .false.
  w = 0.1d0
  wSq = w*w
  
  do iW=1,nWalk

    rij = r12(iW)
    rijSq = rij*rij

    F(iW,1,1:3) = dg(iW,1,1:3)/g(iW,1)
    F(iW,2,1:3) = dg(iW,2,1:3)/g(iW,2)

    if(smoothDrift) then 
      damp = 1d0 + 2d0*w/sqrt(pi)*rij*exp(-wSq*rijSq)/erfc(w*rij)
    else
      damp = 1d0
    endif

    F(iW,1,1:3) = F(iW,1,1:3) - damp*(r(iW,2,1:3) - r(iW,1,1:3))/rijSq
    F(iW,2,1:3) = F(iW,2,1:3) - damp*(r(iW,2,1:3) - r(iW,1,1:3))/rijSq

  enddo

! print*,' F',F

end subroutine drift
