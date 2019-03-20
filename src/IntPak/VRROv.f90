recursive function VRROv(AngMomA,ExpPi,CenterPA) &
                   result(Ga)

! Compute two-electron integrals over Gaussian geminals

  implicit none

! Input variables

  integer,intent(in)            :: AngMomA
  double precision,intent(in)   :: ExpPi
  double precision,intent(in)   :: CenterPA

! Output variables

  double precision              :: Ga

  if(AngMomA < 0) then
    Ga = 0d0
  else
    if(AngMomA == 0) then
      Ga = 1d0
    else
      Ga = CenterPA*VRROv(AngMomA-1,ExpPi,CenterPA) + 0.5d0*dble(AngMomA-1)*ExpPi*VRROv(AngMomA-2,ExpPi,CenterPA)
    end if
  end if

end function VRROv
