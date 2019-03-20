recursive function HRROv(AngMomA,AngMomB,ExpPi,CenterAB,CenterPA) &
                   result(Gab)

! Horizontal recurrence relations for one-electron overlap integrals

  implicit none

! Input variables
  integer,intent(in)            :: AngMomA,AngMomB
  double precision,intent(in)   :: ExpPi
  double precision,intent(in)   :: CenterAB,CenterPA

! Local variables
  double precision              :: VRROv
  double precision              :: Gab

  if(AngMomB < 0) then
    Gab = 0d0
  else
    if(AngMomB == 0) then
      Gab = VRROv(AngMomA,ExpPi,CenterPA)
    else
      Gab = HRROv(AngMomA+1,AngMomB-1,ExpPi,CenterAB,CenterPA)        &
          + CenterAB*HRROv(AngMomA,AngMomB-1,ExpPi,CenterAB,CenterPA)
    end if
  end if

end function HRROv
