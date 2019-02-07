recursive function VRR3e(m,AngMomBra,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
                   result(a1a2a3)

! Vertical recurrence relations for three-electron integrals

  implicit none
  include 'parameters.h'


! Input variables

  integer,intent(in)            :: m
  integer,intent(in)            :: AngMomBra(3,3)
  integer,intent(in)            :: maxm
  double precision,intent(in)   :: Om(0:maxm),ExpZ(3),CenterZA(3,3)
  double precision,intent(in)   :: DY0(3),DY1(3),D2Y0(3,3),D2Y1(3,3)

! Local variables

  logical                       :: NegAngMomBra(3)
  integer                       :: TotAngMomBra(3)
  integer                       :: a1m(3,3),a2m(3,3),a3m(3,3)
  integer                       :: a1mm(3,3),a2mm(3,3),a3mm(3,3)
  integer                       :: a1m2m(3,3),a1m3m(3,3),a2m3m(3,3)
  integer                       :: i,j,xyz

! Output variables

  double precision              :: a1a2a3

  do i=1,3
    NegAngMomBra(i) = AngMomBra(i,1) < 0 .or. AngMomBra(i,2) < 0 .or. AngMomBra(i,3) < 0
    TotAngMomBra(i) = AngMomBra(i,1) + AngMomBra(i,2) + AngMomBra(i,3)
  enddo

!------------------------------------------------------------------------
! Termination condition
!------------------------------------------------------------------------
  if(NegAngMomBra(1) .or. NegAngMomBra(2) .or. NegAngMomBra(3)) then
    a1a2a3 = 0d0
!------------------------------------------------------------------------
! Fundamental integral: (000|000)^m
!------------------------------------------------------------------------
  elseif(TotAngMomBra(1) == 0 .and. TotAngMomBra(2) == 0 .and. TotAngMomBra(3) == 0) then
    a1a2a3 = Om(m)
!------------------------------------------------------------------------
! 1st vertical recurrence relation (4 terms): (a1+00|000)^m
!------------------------------------------------------------------------
  elseif(TotAngMomBra(2) == 0 .and. TotAngMomBra(3) == 0) then
    do i=1,3
      do j=1,3
        a1m(i,j)  = AngMomBra(i,j)
        a1mm(i,j) = AngMomBra(i,j)
      enddo
    enddo
! Loop over cartesian directions
    xyz = 0
    if    (AngMomBra(1,1) > 0) then
      xyz = 1
    elseif(AngMomBra(1,2) > 0) then
      xyz = 2
    elseif(AngMomBra(1,3) > 0) then
      xyz = 3
    else
      write(*,*) 'xyz = 0 in VRR3e!'
    endif
! End loop over cartesian directions
    a1m(1,xyz)  = a1m(1,xyz)  - 1
    a1mm(1,xyz) = a1mm(1,xyz) - 2
    if(AngMomBra(1,xyz) == 0) then
      a1a2a3 = 0d0
    elseif(AngMomBra(1,xyz) == 1) then
      a1a2a3 =                           (CenterZA(1,xyz) - DY0(1))*VRR3e(m,  a1m, maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -                                    (DY1(1) - DY0(1))*VRR3e(m+1,a1m, maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1)  
    else
      a1a2a3 =                           (CenterZA(1,xyz) - DY0(1))*VRR3e(m,  a1m, maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -                                    (DY1(1) - DY0(1))*VRR3e(m+1,a1m, maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             + dble(AngMomBra(1,xyz)-1)*(0.5d0/ExpZ(1) - D2Y0(1,1))*VRR3e(m,  a1mm,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -     dble(AngMomBra(1,xyz)-1)*(D2Y1(1,1) - D2Y0(1,1))*VRR3e(m+1,a1mm,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1)
    endif
!------------------------------------------------------------------------
! 2nd vertical recurrence relation (6 terms): (a1a2+0|000)^m
!------------------------------------------------------------------------
  elseif(TotAngMomBra(3) == 0) then
    do i=1,3
      do j=1,3
        a2m(i,j)   = AngMomBra(i,j)
        a2mm(i,j)  = AngMomBra(i,j)
        a1m2m(i,j) = AngMomBra(i,j)
      enddo
    enddo
! Loop over cartesian directions
    xyz = 0
    if    (AngMomBra(2,1) > 0) then
      xyz = 1
    elseif(AngMomBra(2,2) > 0) then
      xyz = 2
    elseif(AngMomBra(2,3) > 0) then
      xyz = 3
    else
      write(*,*) 'xyz = 0 in VRR3e!'
    endif
! End loop over cartesian directions
    a2m(2,xyz)   = a2m(2,xyz)   - 1
    a2mm(2,xyz)  = a2mm(2,xyz)  - 2
    a1m2m(1,xyz) = a1m2m(1,xyz) - 1
    a1m2m(2,xyz) = a1m2m(2,xyz) - 1
    if(AngMomBra(2,xyz) == 0) then
      a1a2a3 = 0d0
    elseif(AngMomBra(2,xyz) == 1) then
      a1a2a3 =                           (CenterZA(2,xyz) - DY0(2))*VRR3e(m,  a2m,  maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -                                    (DY1(2) - DY0(2))*VRR3e(m+1,a2m,  maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1)  
    else
      a1a2a3 =                           (CenterZA(2,xyz) - DY0(2))*VRR3e(m,  a2m,  maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -                                    (DY1(2) - DY0(2))*VRR3e(m+1,a2m,  maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             + dble(AngMomBra(2,xyz)-1)*(0.5d0/ExpZ(2) - D2Y0(2,2))*VRR3e(m,  a2mm, maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -     dble(AngMomBra(2,xyz)-1)*(D2Y1(2,2) - D2Y0(2,2))*VRR3e(m+1,a2mm, maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) 
    endif
    if(AngMomBra(1,xyz) > 0) &
      a1a2a3 = a1a2a3 &
             +                  dble(AngMomBra(1,xyz))*(-D2Y0(2,1))*VRR3e(m,  a1m2m,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -       dble(AngMomBra(1,xyz))*(D2Y1(2,1) - D2Y0(2,1))*VRR3e(m+1,a1m2m,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1)  
!------------------------------------------------------------------------
! 3rd vertical recurrence relation (8 terms): (a1a2a3+|000)^m
!------------------------------------------------------------------------
  else
    do i=1,3
      do j=1,3
        a3m(i,j)   = AngMomBra(i,j)
        a3mm(i,j)  = AngMomBra(i,j)
        a1m3m(i,j) = AngMomBra(i,j)
        a2m3m(i,j) = AngMomBra(i,j)
      enddo
    enddo
! Loop over cartesian directions
    xyz = 0
    if    (AngMomBra(3,1) > 0) then
      xyz = 1
    elseif(AngMomBra(3,2) > 0) then
      xyz = 2
    elseif(AngMomBra(3,3) > 0) then
      xyz = 3
    else
      write(*,*) 'xyz = 0 in VRR3e!'
    endif
! End loop over cartesian directions
    a3m(3,xyz)   = a3m(3,xyz)   - 1
    a3mm(3,xyz)  = a3mm(3,xyz)  - 2
    a1m3m(1,xyz) = a1m3m(1,xyz) - 1
    a1m3m(3,xyz) = a1m3m(3,xyz) - 1
    a2m3m(2,xyz) = a2m3m(2,xyz) - 1
    a2m3m(3,xyz) = a2m3m(3,xyz) - 1
    if(AngMomBra(3,xyz) == 0) then
      a1a2a3 = 0d0
    elseif(AngMomBra(3,xyz) == 1) then
      a1a2a3 =                           (CenterZA(3,xyz) - DY0(3))*VRR3e(m,  a3m,  maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -                                    (DY1(3) - DY0(3))*VRR3e(m+1,a3m,  maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1)  
    else
      a1a2a3 =                           (CenterZA(3,xyz) - DY0(3))*VRR3e(m,  a3m,  maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -                                    (DY1(3) - DY0(3))*VRR3e(m+1,a3m,  maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             + dble(AngMomBra(3,xyz)-1)*(0.5d0/ExpZ(3) - D2Y0(3,3))*VRR3e(m,  a3mm, maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -     dble(AngMomBra(3,xyz)-1)*(D2Y1(3,3) - D2Y0(3,3))*VRR3e(m+1,a3mm, maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1)
    endif
    if(dble(AngMomBra(1,xyz)) > 0) &
      a1a2a3 = a1a2a3 &
             +                  dble(AngMomBra(1,xyz))*(-D2Y0(3,1))*VRR3e(m,  a1m3m,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -       dble(AngMomBra(1,xyz))*(D2Y1(3,1) - D2Y0(3,1))*VRR3e(m+1,a1m3m,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1)  
    if(dble(AngMomBra(2,xyz)) > 0) &
      a1a2a3 = a1a2a3 &
             +                  dble(AngMomBra(2,xyz))*(-D2Y0(3,2))*VRR3e(m,  a2m3m,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1) &
             -       dble(AngMomBra(2,xyz))*(D2Y1(3,2) - D2Y0(3,2))*VRR3e(m+1,a2m3m,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1)  
  endif

end function VRR3e
