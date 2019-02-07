recursive function HRR3e(AngMomBra,AngMomKet,maxm,Om,ExpZ,CenterAB,CenterZA,DY0,DY1,D2Y0,D2Y1) &
                   result(a1a2a3b1b2b3)

! Horizontal recurrence relations for three-electron integrals

  implicit none
  include 'parameters.h'


! Input variables

  integer,intent(in)            :: AngMomBra(3,3),AngMomKet(3,3)
  integer,intent(in)            :: maxm
  double precision,intent(in)   :: Om(0:maxm),ExpZ(3),CenterAB(3,3),CenterZA(3,3)
  double precision,intent(in)   :: DY0(3),DY1(3),D2Y0(3,3),D2Y1(3,3)

! Local variables

  logical                       :: NegAngMomKet(3)
  integer                       :: TotAngMomBra(3),TotAngMomKet(3)
  integer                       :: a1p(3,3),b1m(3,3),a2p(3,3),b2m(3,3),a3p(3,3),b3m(3,3)
  integer                       :: i,j,xyz
  double precision              :: VRR3e

! Output variables

  double precision              :: a1a2a3b1b2b3

  do i=1,3
    NegAngMomKet(i) = AngMomKet(i,1) < 0 .or. AngMomKet(i,2) < 0 .or. AngMomKet(i,3) < 0
    TotAngMomBra(i) = AngMomBra(i,1) + AngMomBra(i,2) + AngMomBra(i,3)
    TotAngMomKet(i) = AngMomKet(i,1) + AngMomKet(i,2) + AngMomKet(i,3)
  enddo

!------------------------------------------------------------------------
! Termination condition
!------------------------------------------------------------------------
  if(NegAngMomKet(1) .or. NegAngMomKet(2) .or. NegAngMomKet(3)) then
    a1a2a3b1b2b3 = 0d0
!------------------------------------------------------------------------
! 1st and 2nd vertical recurrence relations: <a1a2a3|000>
!------------------------------------------------------------------------
  elseif(TotAngMomKet(1) == 0 .and. TotAngMomKet(2) == 0 .and. TotAngMomKet(3) == 0) then
    a1a2a3b1b2b3 = VRR3e(0,AngMomBra,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1)
!------------------------------------------------------------------------
! 1st horizontal recurrence relation (2 terms): <a1a2a3|b1+00>
!------------------------------------------------------------------------
  elseif(TotAngMomKet(2) == 0 .and. TotAngMomKet(3) == 0) then
    do i=1,3
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
      write(*,*) 'xyz = 0 in HRR3e!'
    endif
! End loop over cartesian directions
    a1p(1,xyz) = a1p(1,xyz) + 1
    b1m(1,xyz) = b1m(1,xyz) - 1
    a1a2a3b1b2b3 = HRR3e(a1p,b1m,maxm,Om,ExpZ,CenterAB,CenterZA,DY0,DY1,D2Y0,D2Y1) &
                 + CenterAB(1,xyz)*                                                    &
                   HRR3e(AngMomBra,b1m,maxm,Om,ExpZ,CenterAB,CenterZA,DY0,DY1,D2Y0,D2Y1)
!------------------------------------------------------------------------
! 2nd horizontal recurrence relation (2 terms):  <a1a2a3|b1b2+0>
!------------------------------------------------------------------------
  elseif(TotAngMomKet(3) == 0) then
    do i=1,3
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
      write(*,*) 'xyz = 0 in HRR3e!'
    endif
! End loop over cartesian directions
    a2p(2,xyz) = a2p(2,xyz) + 1
    b2m(2,xyz) = b2m(2,xyz) - 1
    a1a2a3b1b2b3 = HRR3e(a2p,b2m,maxm,Om,ExpZ,CenterAB,CenterZA,DY0,DY1,D2Y0,D2Y1) &
                 + CenterAB(2,xyz)*                                                    &
                   HRR3e(AngMomBra,b2m,maxm,Om,ExpZ,CenterAB,CenterZA,DY0,DY1,D2Y0,D2Y1)
!------------------------------------------------------------------------
! 3rd horizontal recurrence relation (2 terms):  <a1a2a3|b1b2b3+>
!------------------------------------------------------------------------
  else
    do i=1,3
      do j=1,3
        a3p(i,j) = AngMomBra(i,j)
        b3m(i,j) = AngMomKet(i,j)
      enddo
    enddo
! Loop over cartesian directions
    xyz = 0
    if    (AngMomKet(3,1) > 0) then
      xyz = 1
    elseif(AngMomKet(3,2) > 0) then
      xyz = 2
    elseif(AngMomKet(3,3) > 0) then
      xyz = 3
    else
      write(*,*) 'xyz = 0 in HRR3e!'
    endif
! End loop over cartesian directions
    a3p(3,xyz) = a3p(3,xyz) + 1
    b3m(3,xyz) = b3m(3,xyz) - 1
    a1a2a3b1b2b3 = HRR3e(a3p,b3m,maxm,Om,ExpZ,CenterAB,CenterZA,DY0,DY1,D2Y0,D2Y1) &
                 + CenterAB(3,xyz)*                                                    &
                   HRR3e(AngMomBra,b3m,maxm,Om,ExpZ,CenterAB,CenterZA,DY0,DY1,D2Y0,D2Y1)
  endif

end function HRR3e
