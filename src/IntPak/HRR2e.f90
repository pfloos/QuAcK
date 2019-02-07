recursive function HRR2e(AngMomBra,AngMomKet,       &
                         maxm,Om,ExpZi,ExpY,        &
                         CenterAB,CenterZA,CenterY) &
                   result(a1a2b1b2)

! Horintal recurrence relations for two-electron integrals

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: AngMomBra(2,3),AngMomKet(2,3)
  integer,intent(in)            :: maxm
  double precision,intent(in)   :: Om(0:maxm),ExpZi(2),ExpY(2,2)
  double precision,intent(in)   :: CenterAB(2,3),CenterZA(2,3),CenterY(2,2,3)

! Local variables

  logical                       :: NegAngMomKet(2)
  integer                       :: TotAngMomBra(2),TotAngMomKet(2)
  integer                       :: a1p(2,3),b1m(2,3),a2p(2,3),b2m(2,3)
  integer                       :: i,j,xyz
  double precision              :: VRR2e

! Output variables

  double precision              :: a1a2b1b2

  do i=1,2
    NegAngMomKet(i) = AngMomKet(i,1) < 0 .or. AngMomKet(i,2) < 0 .or. AngMomKet(i,3) < 0
    TotAngMomBra(i) = AngMomBra(i,1) + AngMomBra(i,2) + AngMomBra(i,3)
    TotAngMomKet(i) = AngMomKet(i,1) + AngMomKet(i,2) + AngMomKet(i,3)
  enddo

!------------------------------------------------------------------------
! Termination condition
!------------------------------------------------------------------------
!  if(NegAngMomKet(1) .or. NegAngMomKet(2)) then
!    a1a2b1b2 = 0d0
!------------------------------------------------------------------------
! 1st and 2nd vertical recurrence relations: <a1a2|00>
!------------------------------------------------------------------------
!  elseif(TotAngMomKet(1) == 0 .and. TotAngMomKet(2) == 0) then
  if(TotAngMomKet(1) == 0 .and. TotAngMomKet(2) == 0) then
    a1a2b1b2 = VRR2e(0,AngMomBra,maxm,Om,ExpZi,ExpY,CenterZA,CenterY)
!------------------------------------------------------------------------
! 1st horizontal recurrence relation (2 terms): <a1a2|b1+0>
!------------------------------------------------------------------------
  elseif(TotAngMomKet(2) == 0) then
    do i=1,2
      do j=1,3
        a1p(i,j) = AngMomBra(i,j)
        b1m(i,j) = AngMomKet(i,j)
      enddo
    enddo
! Loop over cartesian directions
    xyz = 0
    if    (AngMomKet(1,1) > 0) then
      xyz = 1
    elseif(AngMomKet(1,2) > 0) then
      xyz = 2
    elseif(AngMomKet(1,3) > 0) then
      xyz = 3
    else
      write(*,*) 'xyz = 0 in HRR2e!'
    endif
! End loop over cartesian directions
    a1p(1,xyz) = a1p(1,xyz) + 1
    b1m(1,xyz) = b1m(1,xyz) - 1
    a1a2b1b2 = HRR2e(a1p,b1m,maxm,Om,ExpZi,ExpY,CenterAB,CenterZA,CenterY) &
             + CenterAB(1,xyz)*HRR2e(AngMomBra,b1m,maxm,Om,ExpZi,ExpY,CenterAB,CenterZA,CenterY)
!------------------------------------------------------------------------
! 2nd horizontal recurrence relation (2 terms): <a1a2|b1b2+>
!------------------------------------------------------------------------
  else
    do i=1,2
      do j=1,3
        a2p(i,j) = AngMomBra(i,j)
        b2m(i,j) = AngMomKet(i,j)
      enddo
    enddo
! Loop over cartesian directions
    xyz = 0
    if    (AngMomKet(2,1) > 0) then
      xyz = 1
    elseif(AngMomKet(2,2) > 0) then
      xyz = 2
    elseif(AngMomKet(2,3) > 0) then
      xyz = 3
    else
      write(*,*) 'xyz = 0 in HRR2e!'
    endif
! End loop over cartesian directions
    a2p(2,xyz) = a2p(2,xyz) + 1
    b2m(2,xyz) = b2m(2,xyz) - 1
    a1a2b1b2 = HRR2e(a2p,b2m,maxm,Om,ExpZi,ExpY,CenterAB,CenterZA,CenterY) &
             + CenterAB(2,xyz)*HRR2e(AngMomBra,b2m,maxm,Om,ExpZi,ExpY,CenterAB,CenterZA,CenterY)
  endif

end function HRR2e
