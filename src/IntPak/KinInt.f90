subroutine KinInt(npKin,nSigpKin,ExpA,CenterA,AngMomA,ExpB,CenterB,AngMomB,pKin)

! Compute one-electron kinetic integrals

  implicit none

! Input variables

  double precision,intent(in)   :: ExpA,ExpB
  double precision,intent(in)   :: CenterA(3),CenterB(3)
  integer,intent(in)            :: AngMomA(3),AngMomB(3)


! Local variables

  double precision              :: ExpAi,ExpBi
  double precision              :: ExpP,ExpPi
  double precision              :: CenterP(3),CenterAB(3),CenterPA(3)
  double precision              :: NormABSq
  double precision              :: GAB
  double precision              :: HRROv,RRKin

  integer                       :: i
  double precision              :: pi
  double precision              :: start_RR,finish_RR,t_RR
  double precision              :: s(3),k(3)

! Output variables

  integer,intent(inout)         :: npKin,nSigpKin
  double precision,intent(out)  :: pKin

  pi = 4d0*atan(1d0)

! Pre-computed shell quantities

  ExpAi = 1d0/ExpA
  ExpBi = 1d0/ExpB

! Pre-computed quantities for shell-pair AB

  ExpP  = ExpA + ExpB
  ExpPi = 1d0/ExpP

  NormABSq = 0d0
  Do i=1,3
    CenterP(i) = (ExpA*CenterA(i) + ExpB*CenterB(i))*ExpPi
    CenterPA(i) = CenterP(i) - CenterA(i)
    CenterAB(i) = CenterA(i) - CenterB(i)
    NormABSq = NormABSq + CenterAB(i)**2
  Enddo

  GAB = (pi*ExpPi)**(1.5d0)*exp(-NormABSq/(ExpAi+ExpBi))

!------------------------------------------------------------------------
! Launch reccurence relations!
!------------------------------------------------------------------------
   call cpu_time(start_RR)
! Loop over cartesian directions
  Do i=1,3
    s(i) = HRROv(AngMomA(i),AngMomB(i),ExpPi,CenterAB(i),CenterPA(i))
    k(i) = RRKin(AngMomA(i),AngMomB(i),ExpA,ExpB,ExpPi,CenterAB(i),CenterPA(i))
  Enddo
  call cpu_time(finish_RR)

  pKin = k(1)*s(2)*s(3) + s(1)*k(2)*s(3) + s(1)*s(2)*k(3)
  pKin = GAB*pKin
  t_RR = finish_RR - start_RR

! Print result
  npKin = npKin + 1
  if(abs(pKin) > 1d-15) then
    nSigpKin = nSigpKin + 1
  end if

end subroutine KinInt
