subroutine NucInt(debug,npNuc,nSigpNuc, &
                  ExpA,CenterA,AngMomA, &
                  ExpB,CenterB,AngMomB, & 
                  CenterC,              &
                  pNuc)

! Compute recursively the primitive one-electron nuclear attraction integrals

  implicit none

! Input variables

  logical,intent(in)            :: debug
  double precision,intent(in)   :: ExpA,ExpB
  double precision,intent(in)   :: CenterA(3),CenterB(3),CenterC(3)
  integer,intent(in)            :: AngMomA(3),AngMomB(3)

! Local variables

  double precision              :: ExpAi,ExpBi
  integer                       :: TotAngMomA,TotAngMomB
  double precision              :: ExpP,ExpPi
  double precision              :: CenterP(3),CenterAB(3),CenterPA(3),CenterPC(3)
  double precision              :: NormABSq,NormPCSq
  double precision              :: G
  double precision,allocatable  :: Om(:)
  double precision              :: HRRNuc
  double precision              :: Gab

  double precision              :: pi
  integer                       :: i,maxm
  double precision              :: start_Om,finish_Om,start_RR,finish_RR,t_Om,t_RR

! Output variables

  integer,intent(inout)           :: npNuc,nSigpNuc
  double precision,intent(out)    :: pNuc

  pi = 4d0*atan(1d0)

! Pre-computed shell quantities

  ExpAi = 1d0/ExpA
  ExpBi = 1d0/ExpB

! Pre-computed quantities for shell-pair AB

  ExpP  = ExpA + ExpB
  ExpPi = 1d0/ExpP

  NormABSq = 0d0
  NormPCSq = 0d0
  do i=1,3
    CenterP(i) = (ExpA*CenterA(i) + ExpB*CenterB(i))*ExpPi
    CenterAB(i) = CenterA(i) - CenterB(i)
    CenterPA(i) = CenterP(i) - CenterA(i)
    CenterPC(i) = CenterP(i) - CenterC(i)
    NormABSq = NormABSq + CenterAB(i)**2
    NormPCSq = NormPCSq + CenterPC(i)**2
  end do

  G = (pi*ExpPi)**(1.5d0)*exp(-NormABSq/(ExpAi+ExpBi))

! Total angular momemtum

  TotAngMomA = AngMomA(1) + AngMomA(2) + AngMomA(3)
  TotAngMomB = AngMomB(1) + AngMomB(2) + AngMomB(3)

  maxm = TotAngMomA + TotAngMomB

! Pre-compute (0|V|0)^m

  allocate(Om(0:maxm))
  call cpu_time(start_Om)
  call CalcOmNuc(maxm,ExpPi,NormPCSq,Om)
  call cpu_time(finish_Om)

! Print (0|V|0)^m

  if(debug) then
    write(*,*) '(0|V|0)^m'
    do i=0,maxm
      write(*,*) i,Om(i)
    end do
    write(*,*)
  end if

!------------------------------------------------------------------------
! Launch reccurence relations!
!------------------------------------------------------------------------
   call cpu_time(start_RR)
   Gab = HRRNuc(AngMomA,AngMomB,maxm,Om,ExpPi,CenterAB,CenterPA,CenterPC)
   call cpu_time(finish_RR)

! Timings

  t_Om = finish_Om - start_Om
  t_RR = finish_RR - start_RR

! Print result

  pNuc = G*Gab

  npNuc = npNuc + 1
  if(abs(pNuc) > 1d-15) then
    nSigpNuc = nSigpNuc + 1
!    write(*,'(A10,1X,F16.10,1X,I6,1X,I6)') '[a|V|b] = ',pNuc
  end if

! Deallocate arrays

  deallocate(Om)

end subroutine NucInt
