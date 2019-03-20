function G2eInt(debug,iType,                &
                ExpG,                       &
                ExpBra,CenterBra,AngMomBra, &
                ExpKet,CenterKet,AngMomKet)

! Compute recursively the primitive two-electron integral [ab|cd]

  implicit none
  include 'parameters.h'


! Input variables

  logical,intent(in)            :: debug
  integer,intent(in)            :: iType
  double precision,intent(in)   :: ExpBra(2),ExpKet(2)
  double precision,intent(in)   :: ExpG
  double precision,intent(in)   :: CenterBra(2,3),CenterKet(2,3)
  integer,intent(in)            :: AngMomBra(2,3),AngMomKet(2,3)

! Local variables

  integer                       :: TotAngMomBra(3),TotAngMomKet(3)
  double precision              :: ExpZi(2),ExpY(2,2)
  double precision              :: CenterZ(2,3),CenterAB(2,3),CenterZA(2,3),CenterY(2,2,3)
  double precision              :: NormABSq(2),NormYSq(2,2)
  double precision              :: GAB(2)
  double precision,allocatable  :: Om(:)
  double precision              :: fG
  double precision              :: HRR2e,VRR2e
  double precision              :: a1a2b1b2

  integer                       :: i,j,k,maxm
  double precision              :: start_Om,finish_Om,start_RR,finish_RR,t_Om,t_RR

! Output variables
  double precision              :: G2eInt

! Pre-computed shell-pair quantities

  do i=1,2
    ExpZi(i)  = 1d0/(ExpBra(i) + ExpKet(i))
  end do

  NormABSq = 0d0
  do j=1,3
    do i=1,2
      CenterZ(i,j) = (ExpBra(i)*CenterBra(i,j) + ExpKet(i)*CenterKet(i,j))*ExpZi(i)
      CenterAB(i,j) = CenterBra(i,j) - CenterKet(i,j)
      CenterZA(i,j) = CenterZ(i,j) - CenterBra(i,j)
      NormABSq(i) = NormABSq(i) + CenterAB(i,j)**2
    end do
  end do

  do i=1,2
    GAB(i) = (pi*ExpZi(i))**(1.5d0)*exp(-ExpBra(i)*ExpKet(i)*NormABSq(i)*ExpZi(i))
  end do

! Pre-computed shell-quartet quantities

  do i=1,2
    do j=1,2
      ExpY(i,j) = 1d0/(ExpZi(i) + ExpZi(j))
    end do
  end do

  do i=1,2
    do j=1,2
      NormYSq(i,j) = 0d0
      do k=1,3   
        CenterY(i,j,k) = CenterZ(i,k) - CenterZ(j,k)
        NormYSq(i,j) = NormYSq(i,j) + CenterY(i,j,k)**2
      end do
    end do
  end do

!  fG = (ExpZ(1)*ExpZ(2)*ExpG)/(ExpZ(1)*ExpZ(2) + ExpZ(1)*ExpG + ExpZ(2)*ExpG)
  fG = 1d0/(ExpZi(1) + 1d0/ExpG + ExpZi(2))

! Total angular momemtum

  maxm = 0
  do i=1,2
    TotAngMomBra(i) = AngMomBra(i,1) + AngMomBra(i,2) + AngMomBra(i,3)
    TotAngMomKet(i) = AngMomKet(i,1) + AngMomKet(i,2) + AngMomKet(i,3)
    maxm = maxm + TotAngMomBra(i) + TotAngMomKet(i)
  end do

! Pre-compute (00|00)^m

  allocate(Om(0:maxm))
  call cpu_time(start_Om)

  if(iType == 1) then
    call CalcOmERI(maxm,ExpY(1,2),NormYSq(1,2),Om)
  elseif(iType == 3) then
    call CalcOmYuk(maxm,ExpG,ExpY(1,2),fG,NormYSq(1,2),Om)
  elseif(iType == 4) then
    call CalcOmErf(maxm,ExpY(1,2),fG,NormYSq(1,2),Om)
  end if

  call cpu_time(finish_Om)

! Print (00|00)^m

  if(debug) then
    write(*,*) '(00|00)^m'
    do i=0,maxm
      write(*,*) i,Om(i)
    end do
    write(*,*)
  end if

!------------------------------------------------------------------------
! Launch reccurence relations!
!------------------------------------------------------------------------
  call cpu_time(start_RR)

  if(TotAngMomKet(1) == 0 .and. TotAngMomKet(2) == 0) then
    if(TotAngMomBra(1) == 0 .and. TotAngMomBra(2) == 0) then
      a1a2b1b2 = Om(0)
    else
      a1a2b1b2 = VRR2e(0,AngMomBra,maxm,Om,ExpZi,ExpY,CenterZA,CenterY)
    end if
  else
    a1a2b1b2 = HRR2e(AngMomBra,AngMomKet,maxm,Om,ExpZi,ExpY,CenterAB,CenterZA,CenterY)
  end if

  call cpu_time(finish_RR)

! Timings

  t_Om = finish_Om - start_Om
  t_RR = finish_RR - start_RR

! Print result

  G2eInt = GAB(1)*GAB(2)*a1a2b1b2

end function G2eInt
