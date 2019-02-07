function RRKin(AngMomA,AngMomB,ExpA,ExpB,ExpPi,CenterAB,CenterPA) &
         result(Gab)

! Recurrence relation for one-electron kinetic integrals

  implicit none

! Input variables
  integer,intent(in)            :: AngMomA,AngMomB
  double precision,intent(in)   :: ExpA,ExpB,ExpPi
  double precision,intent(in)   :: CenterAB,CenterPA

! Local variables
  double precision              :: HRROv
  double precision              :: a,b,s1,s2,s3,s4
  double precision              :: Gab

   a = dble(AngMomA)
   b = dble(AngMomB)

   s1 = HRROv(AngMomA-1,AngMomB-1,ExpPi,CenterAB,CenterPA)
   s2 = HRROv(AngMomA+1,AngMomB-1,ExpPi,CenterAB,CenterPA)
   s3 = HRROv(AngMomA-1,AngMomB+1,ExpPi,CenterAB,CenterPA)
   s4 = HRROv(AngMomA+1,AngMomB+1,ExpPi,CenterAB,CenterPA)

   Gab = 0.5d0*a*b*s1 - ExpA*b*s2 - a*ExpB*s3 + 2d0*ExpA*ExpB*s4
          

end function RRKin
