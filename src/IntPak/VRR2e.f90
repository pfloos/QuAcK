recursive function VRR2e(m,AngMomBra,maxm,Om,ExpZi,ExpY,CenterZA,CenterY) &
                   result(a1a2)

! Compute two-electron integrals over Gaussian geminals

  implicit none
  include 'parameters.h'

! Input variables

  integer,intent(in)            :: m
  integer,intent(in)            :: AngMomBra(2,3)
  integer,intent(in)            :: maxm
  double precision,intent(in)   :: Om(0:maxm),ExpZi(2),ExpY(2,2)
  double precision,intent(in)   :: CenterZA(2,3),CenterY(2,2,3)

! Local variables

  logical                       :: NegAngMomBra(2)
  integer                       :: TotAngMomBra(2)
  integer                       :: a1m(2,3),a2m(2,3)
  integer                       :: a1mm(2,3),a2mm(2,3)
  integer                       :: a1m2m(2,3)
  double precision              :: fZ(2)
  integer                       :: i,j,xyz

! Output variables

  double precision              :: a1a2

  do i=1,2
    NegAngMomBra(i) = AngMomBra(i,1) < 0 .or. AngMomBra(i,2) < 0 .or. AngMomBra(i,3) < 0
    TotAngMomBra(i) = AngMomBra(i,1) + AngMomBra(i,2) + AngMomBra(i,3)
  end do

  fZ(1) = ExpY(1,2)*ExpZi(1)
  fZ(2) = ExpY(1,2)*ExpZi(2)

!------------------------------------------------------------------------
! Termination condition
!------------------------------------------------------------------------
!  if(NegAngMomBra(1) .or. NegAngMomBra(2)) then
!    a1a2 = 0d0
!------------------------------------------------------------------------
! Fundamental integral: (00|00)^m
!------------------------------------------------------------------------
!  elseif(TotAngMomBra(1) == 0 .and. TotAngMomBra(2) == 0) then
  if(TotAngMomBra(1) == 0 .and. TotAngMomBra(2) == 0) then
    a1a2 = Om(m)
!------------------------------------------------------------------------
! 1st vertical recurrence relation (4 terms): (a+0|00)^m
!------------------------------------------------------------------------
  elseif(TotAngMomBra(2) == 0) then
    do i=1,2
      do j=1,3
        a1m(i,j)  = AngMomBra(i,j)
        a1mm(i,j) = AngMomBra(i,j)
      end do
    end do
! Loop over cartesian directions
    xyz = 0
    if    (AngMomBra(1,1) > 0) then
      xyz = 1
    elseif(AngMomBra(1,2) > 0) then
      xyz = 2
    elseif(AngMomBra(1,3) > 0) then
      xyz = 3
    else
      write(*,*) 'xyz = 0 in VRR2e!'
    end if
! End loop over cartesian directions
    a1m(1,xyz)  = a1m(1,xyz)  - 1
    a1mm(1,xyz) = a1mm(1,xyz) - 2
    if(AngMomBra(1,xyz) <= 0) then
      a1a2 = 0d0
    elseif(AngMomBra(1,xyz) == 1) then
      a1a2 = CenterZA(1,xyz)*VRR2e(m,a1m,maxm,Om,ExpZi,ExpY,CenterZA,CenterY)          &
           - fZ(1)*CenterY(1,2,xyz)*VRR2e(m+1,a1m,maxm,Om,ExpZi,ExpY,CenterZA,CenterY) 
    else
      a1a2 = CenterZA(1,xyz)*VRR2e(m,a1m,maxm,Om,ExpZi,ExpY,CenterZA,CenterY)          &
           - fZ(1)*CenterY(1,2,xyz)*VRR2e(m+1,a1m,maxm,Om,ExpZi,ExpY,CenterZA,CenterY) &
           + 0.5d0*dble(AngMomBra(1,xyz)-1)*ExpZi(1)*(                                 &
               VRR2e(m,a1mm,maxm,Om,ExpZi,ExpY,CenterZA,CenterY)                       &
               - fZ(1)*VRR2e(m+1,a1mm,maxm,Om,ExpZi,ExpY,CenterZA,CenterY))
    end if
!------------------------------------------------------------------------
! 2nd vertical recurrence relation (5 terms): (a0|c+0)^m
!------------------------------------------------------------------------
  else
    do i=1,2
      do j=1,3
        a2m(i,j)   = AngMomBra(i,j)
        a2mm(i,j)  = AngMomBra(i,j)
        a1m2m(i,j) = AngMomBra(i,j)
      end do
    end do
! Loop over cartesian directions
    xyz = 0
    if    (AngMomBra(2,1) > 0) then
      xyz = 1
    elseif(AngMomBra(2,2) > 0) then
      xyz = 2
    elseif(AngMomBra(2,3) > 0) then
      xyz = 3
    else
      write(*,*) 'xyz = 0 in VRR2e!'
    end if
! End loop over cartesian directions
    a2m(2,xyz)   = a2m(2,xyz)   - 1
    a2mm(2,xyz)  = a2mm(2,xyz)  - 2
    a1m2m(1,xyz) = a1m2m(1,xyz) - 1
    a1m2m(2,xyz) = a1m2m(2,xyz) - 1
    if(AngMomBra(2,xyz) <= 0) then
      a1a2 = 0d0
    elseif(AngMomBra(2,xyz) == 1) then
        a1a2 = CenterZA(2,xyz)*VRR2e(m,a2m,maxm,Om,ExpZi,ExpY,CenterZA,CenterY)          &
             + fZ(2)*CenterY(1,2,xyz)*VRR2e(m+1,a2m,maxm,Om,ExpZi,ExpY,CenterZA,CenterY)  
    else
      a1a2 = CenterZA(2,xyz)*VRR2e(m,a2m,maxm,Om,ExpZi,ExpY,CenterZA,CenterY)          &
           + fZ(2)*CenterY(1,2,xyz)*VRR2e(m+1,a2m,maxm,Om,ExpZi,ExpY,CenterZA,CenterY) &
           + 0.5d0*dble(AngMomBra(2,xyz)-1)*ExpZi(2)*(                                 &
               VRR2e(m,a2mm,maxm,Om,ExpZi,ExpY,CenterZA,CenterY)                       &
             - fZ(2)*VRR2e(m+1,a2mm,maxm,Om,ExpZi,ExpY,CenterZA,CenterY))               
    end if
    if(AngMomBra(1,xyz) > 0) &
      a1a2 = a1a2 &
           + 0.5d0*dble(AngMomBra(1,xyz))*fZ(2)*ExpZi(1)*VRR2e(m+1,a1m2m,maxm,Om,ExpZi,ExpY,CenterZA,CenterY)
  end if

end function VRR2e
