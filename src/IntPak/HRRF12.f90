recursive function HRRF12(AngMomA,AngMomB,AngMomC,AngMomD,fG,gP,gG,gQ,ExpPGQi, &
                        CenterPQSq,CenterRA,CenterRC,CenterAB,CenterCD)      &
                   result(Gabcd)

! Compute two-electron integrals over Gaussian geminals

  implicit none

! Input variables
  integer,intent(in)            :: AngMomA,AngMomB,AngMomC,AngMomD
  double precision,intent(in)   :: ExpPGQi
  double precision,intent(in)   :: fG,gP,gG,gQ
  double precision,intent(in)   :: CenterPQSq,CenterRA,CenterRC
  double precision,intent(in)   :: CenterAB,CenterCD

! Local variables
  double precision              :: VRRF12
  double precision              :: Gabcd

  If(AngMomB < 0 .or. AngMomD < 0) then
    Gabcd = 0d0
  Else
    If(AngMomB == 0 .and. AngMomD == 0) then
      Gabcd = VRRF12(AngMomA,AngMomC,fG,gP,gG,gQ,ExpPGQi,CenterPQSq,CenterRA,CenterRC)
    Else
      If(AngMomD == 0) then
        Gabcd = HRRF12(AngMomA+1,AngMomB-1,AngMomC,AngMomD,fG,gP,gG,gQ,ExpPGQi, &
                  CenterPQSq,CenterRA,CenterRC,CenterAB,CenterCD)             &
              + CenterAB*HRRF12(AngMomA,AngMomB-1,AngMomC,AngMomD,fG,gP,gG,gQ,  &
                  ExpPGQi,CenterPQSq,CenterRA,CenterRC,CenterAB,CenterCD)
      Else
        Gabcd = HRRF12(AngMomA,AngMomB,AngMomC+1,AngMomD-1,fG,gP,gG,gQ,ExpPGQi, &
                  CenterPQSq,CenterRA,CenterRC,CenterAB,CenterCD)             &
              + CenterCD*HRRF12(AngMomA,AngMomB,AngMomC,AngMomD-1,fG,gP,gG,gQ,  &
                  ExpPGQi,CenterPQSq,CenterRA,CenterRC,CenterAB,CenterCD)
      EndIf
    EndIf
  EndIf

end function HRRF12
