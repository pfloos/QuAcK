recursive function VRRF12(AngMomA,AngMomC,fG,gP,gG,gQ,ExpPGQi,CenterPQSq,CenterRA,CenterRC) &
                   result(Gac)

! Compute two-electron integrals over Gaussian geminals

  implicit none

! Input variables

  integer,intent(in)            :: AngMomA,AngMomC
  double precision,intent(in)   :: ExpPGQi
  double precision,intent(in)   :: fG,gP,gG,gQ
  double precision,intent(in)   :: CenterPQSq,CenterRA,CenterRC

! Output variables

  double precision              :: Gac

  if(AngMomA < 0 .or. AngMomC < 0) then
    Gac = 0d0
  else
    if(AngMomA == 0 .and. AngMomC == 0) then
      Gac = sqrt(fG)*exp(-CenterPQSq/ExpPGQi)
    else
      If(AngMomC == 0) then
        Gac = CenterRA*VRRF12(AngMomA-1,AngMomC,fG,gP,gG,gQ,ExpPGQi,CenterPQSq,CenterRA,CenterRC)           &
            + dble(AngMomA-1)*gP*VRRF12(AngMomA-2,AngMomC,fG,gP,gG,gQ,ExpPGQi,CenterPQSq,CenterRA,CenterRC)
      else
        Gac = CenterRC*VRRF12(AngMomA,AngMomC-1,fG,gP,gG,gQ,ExpPGQi,CenterPQSq,CenterRA,CenterRC)           &
            + dble(AngMomA)*gG*VRRF12(AngMomA-1,AngMomC-1,fG,gP,gG,gQ,ExpPGQi,CenterPQSq,CenterRA,CenterRC) &
            + dble(AngMomC-1)*gQ*VRRF12(AngMomA,AngMomC-2,fG,gP,gG,gQ,ExpPGQi,CenterPQSq,CenterRA,CenterRC)
      endIf
    endIf
  endIf

end function VRRF12
