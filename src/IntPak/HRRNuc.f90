recursive function HRRNuc(AngMomA,AngMomB,maxm,Om,ExpPi,CenterAB,CenterPA,CenterPC) &
                   result(Gab)

! Horizontal recurrence relation for one-electron nuclear attraction integrals

  implicit none

! Input variables

  integer,intent(in)            :: AngMomA(3),AngMomB(3)
  integer,intent(in)            :: maxm
  double precision,intent(in)   :: Om(0:maxm)
  double precision,intent(in)   :: ExpPi
  double precision,intent(in)   :: CenterAB(3),CenterPA(3),CenterPC(3)

! Local variables

  logical                       :: NegAngMomB
  integer                       :: TotAngMomA,TotAngMomB
  integer                       :: xyz,ap(3),bm(3)
  integer                       :: i
  double precision              :: VRRNuc

! Output variables

  double precision              :: Gab

  NegAngMomB = AngMomB(1) < 0 .or. AngMomB(2) < 0 .or. AngMomB(3) < 0

  TotAngMomA = AngMomA(1) + AngMomA(2) + AngMomA(3)
  TotAngMomB = AngMomB(1) + AngMomB(2) + AngMomB(3)

!------------------------------------------------------------------------
! Termination condition
!------------------------------------------------------------------------
  if(NegAngMomB) then
    Gab = 0d0
  else
!------------------------------------------------------------------------
! Vertical recurrence relations: (a|0)
!------------------------------------------------------------------------
    if(TotAngMomB == 0) then
      Gab = VRRNuc(0,AngMomA,maxm,Om,ExpPi,CenterAB,CenterPA,CenterPC)
    else
!------------------------------------------------------------------------
! 1st horizontal recurrence relation (2 terms): (a|b+)
!------------------------------------------------------------------------
      do i=1,3
        ap(i) = AngMomA(i)
        bm(i) = AngMomB(i)
      enddo
! Loop over cartesian directions
      xyz = 0
      if    (AngMomB(1) > 0) then
        xyz = 1
      elseif(AngMomB(2) > 0) then
        xyz = 2
      elseif(AngMomB(3) > 0) then
        xyz = 3
      else
        write(*,*) 'xyz = 0 in HRRNuc!'
      endif
! End loop over cartesian directions
      ap(xyz) = ap(xyz) + 1
      bm(xyz) = bm(xyz) - 1
      Gab = HRRNuc(ap,bm,maxm,Om,ExpPi,CenterAB,CenterPA,CenterPC)                    &
          + CenterAB(xyz)*HRRNuc(AngMomA,bm,maxm,Om,ExpPi,CenterAB,CenterPA,CenterPC)
    endif
  endif

end function HRRNuc
