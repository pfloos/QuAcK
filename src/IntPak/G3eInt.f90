function G3eInt(debug,iType,                &
                ExpG13,ExpG23,              &
                ExpBra,CenterBra,AngMomBra, &
                ExpKet,CenterKet,AngMomKet)

! Compute two-electron integrals over the Yukawa operator

  implicit none
  include 'parameters.h'


! Input variables

  logical,intent(in)            :: debug
  integer,intent(in)            :: iType
  double precision,intent(in)   :: ExpG13,ExpG23
  double precision,intent(in)   :: ExpBra(3),ExpKet(3)
  double precision,intent(in)   :: CenterBra(3,3),CenterKet(3,3)
  integer,intent(in)            :: AngMomBra(3,3),AngMomKet(3,3)

! Local variables

  double precision              :: ExpG(3,3)
  integer                       :: TotAngMomBra(3),TotAngMomKet(3)
  double precision              :: ExpZ(3)
  double precision              :: CenterZ(3,3),CenterAB(3,3),CenterZA(3,3)
  double precision              :: NormABSq(3)
  double precision              :: GAB(3)
  double precision,allocatable  :: Om(:)
  double precision              :: HRR3e,VRR3e

  double precision              :: DY0(3),DY1(3),D2Y0(3,3),D2Y1(3,3)
  double precision              :: delta0,delta1,Y0,Y1

  integer                       :: i,j,maxm
  double precision              :: start_Om,finish_Om,t_Om,start_RR,finish_RR,t_RR
  double precision              :: a1a2a3b1b2b3

! Output variables
  double precision              :: G3eInt

! Gaussian geminal exponents

  ExpG = 0d0
  ExpG(1,3) = ExpG13
  ExpG(2,3) = ExpG23

! Pre-computed quantities for shell-pair

  do i=1,3
    ExpZ(i)  = ExpBra(i) + ExpKet(i)
  end do 

  NormABSq = 0d0
  do i=1,3
    do j=1,3
      CenterZ(i,j) = (ExpBra(i)*CenterBra(i,j) + ExpKet(i)*CenterKet(i,j))/ExpZ(i)
      CenterAB(i,j) = CenterBra(i,j) - CenterKet(i,j)
      CenterZA(i,j) = CenterZ(i,j) - CenterBra(i,j)
      NormABSq(i) = NormABSq(i) + CenterAB(i,j)**2
    end do
  end do

  do i=1,3
    GAB(i) = (pi/ExpZ(i))**(1.5d0)*exp(-ExpBra(i)*ExpKet(i)*NormABSq(i)/ExpZ(i))
  end do

! Pre-computed shell-sextet quantities

  call FormVRR3e(ExpZ,ExpG,CenterZ,DY0,DY1,D2Y0,D2Y1,delta0,delta1,Y0,Y1)

! Total angular momemtum

  maxm = 0
  do i=1,3
    TotAngMomBra(i) = AngMomBra(i,1) + AngMomBra(i,2) + AngMomBra(i,3)
    TotAngMomKet(i) = AngMomKet(i,1) + AngMomKet(i,2) + AngMomKet(i,3)
    maxm = maxm + TotAngMomBra(i) + TotAngMomKet(i)
  end do

! Pre-compute (000|000)^m

  allocate(Om(0:maxm))
  call cpu_time(start_Om)
  call CalcOm3e(maxm,delta0,delta1,Y0,Y1,Om)
  call cpu_time(finish_Om)

! Print (000|000)^m

  if(.false.) then
    write(*,*) '(000|000)^m'
    do i=0,maxm
      write(*,*) i,Om(i)
    end do
    write(*,*)
  end if

!------------------------------------------------------------------------
! Launch reccurence relations!
!------------------------------------------------------------------------
  call cpu_time(start_RR)
  if(TotAngMomKet(1) == 0 .and. TotAngMomKet(2) == 0 .and. TotAngMomKet(3) == 0) then
    if(TotAngMomBra(1) == 0 .and. TotAngMomBra(2) == 0 .and.  TotAngMomBra(3) == 0) then
      a1a2a3b1b2b3 = Om(0)
    else
    a1a2a3b1b2b3 = VRR3e(0,AngMomBra,maxm,Om,ExpZ,CenterZA,DY0,DY1,D2Y0,D2Y1)
    end if
  else
    a1a2a3b1b2b3 = HRR3e(AngMomBra,AngMomKet,maxm,Om,ExpZ,CenterAB,CenterZA,DY0,DY1,D2Y0,D2Y1)
  end if

  
  call cpu_time(finish_RR)

! Timings

  t_Om = finish_Om - start_Om
  t_RR = finish_RR - start_RR

! Print result

  G3eInt = GAB(1)*GAB(2)*GAB(3)*a1a2a3b1b2b3

end function G3eInt
