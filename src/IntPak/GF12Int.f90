function GF12Int(ExpG,ExpA,CenterA,AngMomA,ExpB,CenterB,AngMomB,ExpC,CenterC,AngMomC,ExpD,CenterD,AngMomD)

! Compute two-electron integrals over Gaussian geminals

  implicit none

! Input variables

  double precision,intent(in)   :: ExpG
  double precision,intent(in)   :: ExpA,ExpB,ExpC,ExpD
  double precision,intent(in)   :: CenterA(3),CenterB(3),CenterC(3),CenterD(3)
  integer,intent(in)            :: AngMomA(3),AngMomB(3),AngMomC(3),AngMomD(3)


! Local variables

  double precision              :: ExpAi,ExpBi,ExpCi,ExpDi,ExpGi
  double precision              :: ExpP,ExpQ,ExpPi,ExpQi,ExpPGQi
  double precision              :: CenterP(3),CenterQ(3),CenterAB(3),CenterCD(3),CenterPQSq(3),CenterRA(3),CenterRC(3)
  double precision              :: NormABSq,NormCDSq
  double precision              :: GAB,GCD
  double precision              :: fP,fG,fQ,gP,gG,gQ
  double precision              :: HRRF12

  integer                       :: i
  double precision              :: pi
  double precision              :: start_RR,finish_RR,t_RR
  double precision              :: Gabcd(3)

! Output variables
  double precision              :: GF12Int

  pi = 4d0*atan(1d0)

! Pre-computed shell quantities

  ExpAi = 1d0/ExpA
  ExpBi = 1d0/ExpB
  ExpCi = 1d0/ExpC
  ExpDi = 1d0/ExpD
  ExpGi = 1d0/ExpG

! Pre-computed quantities for shell-pair AB

  ExpP  = ExpA + ExpB
  ExpPi = 1d0/ExpP

  NormABSq = 0d0
  Do i=1,3
    CenterP(i) = (ExpA*CenterA(i) + ExpB*CenterB(i))*ExpPi
    CenterAB(i) = CenterA(i) - CenterB(i)
    NormABSq = NormABSq + CenterAB(i)**2
  Enddo

  GAB = (pi*ExpPi)**(1.5d0)*exp(-NormABSq/(ExpAi+ExpBi))

! Pre-computed quantities for shell-pair CD

  ExpQ  = ExpC + ExpD
  ExpQi = 1d0/ExpQ

  NormCDSq = 0d0
  Do i=1,3
    CenterQ(i) = (ExpC*CenterC(i) + ExpD*CenterD(i))*ExpQi
    CenterCD(i) = CenterC(i) - CenterD(i)
    NormCDSq = NormCDSq + CenterCD(i)**2
  Enddo

  GCD = (pi*ExpQi)**(1.5d0)*exp(-NormCDSq/(ExpCi+ExpDi))

! Pre-computed shell-quartet quantities

  ExpPGQi = ExpPi + ExpGi + ExpQi

  Do i=1,3
    CenterPQSq(i) = (CenterP(i) - CenterQ(i))**2
  Enddo

  fP = ExpPi/ExpPGQi
  fG = ExpGi/ExpPGQi
  fQ = ExpQi/ExpPGQi

  gP = (1d0 - fP)*0.5d0*ExpPi
  gG = fP*0.5d0*expQi
  gQ = (1d0 - fQ)*0.5d0*ExpQi

  do i=1,3
    CenterRA(i) = CenterP(i) - CenterA(i) + fP*(CenterQ(i) - CenterP(i))
    CenterRC(i) = CenterQ(i) - CenterC(i) + fQ*(CenterP(i) - CenterQ(i))
  enddo
!------------------------------------------------------------------------
! Launch reccurence relations!
!------------------------------------------------------------------------
   call cpu_time(start_RR)
! Loop over cartesian directions
  Do i=1,3
    Gabcd(i) = HRRF12(AngMomA(i),AngMomB(i),AngMomC(i),AngMomD(i),fG,gP,gG,gQ,ExpPGQi, &
                 CenterPQSq(i),CenterRA(i),CenterRC(i),CenterAB(i),CenterCD(i))
  Enddo
  call cpu_time(finish_RR)

! Print result

  GF12Int = GAB*GCD*Gabcd(1)*Gabcd(2)*Gabcd(3)
  t_RR = finish_RR - start_RR

end function GF12Int
