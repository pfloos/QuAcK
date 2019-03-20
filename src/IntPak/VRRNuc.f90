recursive function VRRNuc(m,AngMomA,maxm,Om,ExpPi,CenterAB,CenterPA,CenterPC) &
                   result(Ga)

! Compute two-electron integrals over Gaussian geminals

  implicit none

! Input variables

  integer,intent(in)            :: m
  integer,intent(in)            :: AngMomA(3)
  integer,intent(in)            :: maxm
  double precision,intent(in)   :: Om(0:maxm)
  double precision,intent(in)   :: ExpPi
  double precision,intent(in)   :: CenterAB(3),CenterPA(3),CenterPC(3)

! Local variables

  logical                       :: NegAngMomA
  integer                       :: TotAngMomA
  integer                       :: xyz,am(3),amm(3)
  integer                       :: i

! Output variables

  double precision              :: Ga

  NegAngMomA = AngMomA(1) < 0 .or. AngMomA(2) < 0 .or. AngMomA(3) < 0
  TotAngMomA = AngMomA(1) + AngMomA(2) + AngMomA(3)

!------------------------------------------------------------------------
! Termination condition
!------------------------------------------------------------------------

  if(NegAngMomA) then

    Ga = 0d0

  else
!------------------------------------------------------------------------
! Fundamental integral: (0|0)^m
!------------------------------------------------------------------------
    if(TotAngMomA == 0) then

      Ga = Om(m)

    else
!------------------------------------------------------------------------
! Vertical recurrence relation (4 terms): (a+|0)^m
!------------------------------------------------------------------------
      do i=1,3
        am(i)  = AngMomA(i)
        amm(i) = AngMomA(i)
      end do
! Loop over cartesian directions
      xyz = 0
      if    (AngMomA(1) > 0) then
        xyz = 1
      elseif(AngMomA(2) > 0) then
        xyz = 2
      elseif(AngMomA(3) > 0) then
        xyz = 3
      else
        write(*,*) 'xyz = 0 in VRRNuc!'
      end if
! End loop over cartesian directions
      am(xyz)  = am(xyz)  - 1
      amm(xyz) = amm(xyz) - 2
      Ga = CenterPA(xyz)*VRRNuc(m,am,maxm,Om,ExpPi,CenterAB,CenterPA,CenterPC)                   &
         + 0.5d0*dble(am(xyz))*ExpPi*VRRNuc(m,amm,maxm,Om,ExpPi,CenterAB,CenterPA,CenterPC)      &
         - CenterPC(xyz)*ExpPi*VRRNuc(m+1,am,maxm,Om,ExpPi,CenterAB,CenterPA,CenterPC)           &
         - 0.5d0*dble(am(xyz))*ExpPi**2*VRRNuc(m+1,amm,maxm,Om,ExpPi,CenterAB,CenterPA,CenterPC)
    end if
  end if

end function VRRNuc
