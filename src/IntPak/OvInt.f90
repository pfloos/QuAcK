subroutine OvInt(npOv,nSigpOv,ExpA,CenterA,AngMomA,ExpB,CenterB,AngMomB,pOv)

! Compute one-electron overlap integrals

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
  double precision              :: G
  double precision              :: HRROv

  integer                       :: i
  double precision              :: pi
  double precision              :: start_RR,finish_RR,t_RR
  double precision              :: Gab(3)

! Output variables

  integer,intent(inout)         :: npOv,nSigpOv
  double precision,intent(out)  :: pOv

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

  G = (pi*ExpPi)**(1.5d0)*exp(-NormABSq/(ExpAi+ExpBi))

!------------------------------------------------------------------------
! Launch reccurence relations!
!------------------------------------------------------------------------
   call cpu_time(start_RR)
! Loop over cartesian directions
  Do i=1,3
    Gab(i) = HRROv(AngMomA(i),AngMomB(i),ExpPi,CenterAB(i),CenterPA(i))
  Enddo
  call cpu_time(finish_RR)

  pOv = G*Gab(1)*Gab(2)*Gab(3)
  t_RR = finish_RR - start_RR

! Print result
  npOv = npOv + 1
  if(abs(pOv) > 1d-15) then
    nSigpOv = nSigpOv + 1
  endif

end subroutine OvInt
