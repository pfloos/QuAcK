subroutine FormVRR3e(ExpZ,ExpG,CenterZ,DY0,DY1,D2Y0,D2Y1,delta0,delta1,Y0,Y1)

! Form stuff we need...

  implicit none
  include 'parameters.h'


! Input variables

  double precision,intent(in)   :: ExpZ(3),ExpG(3,3)
  double precision,intent(in)   :: CenterZ(3,3)

! Local variables

  integer                       :: i,j,k,l
  double precision              :: ZetaMat(3,3)
  double precision              :: CMat(3,3),GMat(3,3)
  double precision              :: Delta0Mat(3,3),Delta1Mat(3,3)
  double precision              :: InvDelta0Mat(3,3),InvDelta1Mat(3,3)
  double precision              :: CenterY(3,3,3)
  double precision              :: YMat(3,3),Y2Mat(3,3)
  double precision              :: DYMat(3,3,3),D2YMat(3,3,3,3)
  double precision              :: D0Mat(3,3),D1Mat(3,3)

  double precision              :: KappaCross

! Output variables

  double precision,intent(out)  :: DY0(3),DY1(3),D2Y0(3,3),D2Y1(3,3)
  double precision,intent(out)  :: delta0,delta1,Y0,Y1

! Initalize arrays

  ZetaMat = 0d0
  CMat    = 0d0
  GMat    = 0d0
  YMat    = 0d0
  Y2Mat   = 0d0
  D0Mat   = 0d0
  D1Mat   = 0d0

! Form the zeta matrix Eq. (15a)

  do i=1,3
    ZetaMat(i,i) = ExpZ(i)
  enddo

!  print*,'Zeta'
!  call matout(3,3,ZetaMat)

! Form the C matrix Eq. (15a)

  CMat(1,1) = 1d0
  CMat(2,2) = 1d0
  CMat(1,2) = -1d0
  CMat(2,1) = -1d0

!  print*,'C'
!  call matout(3,3,CMat)

! Form the G matrix Eq. (15b)

  do i=1,3
    do j=1,i-1
      GMat(i,j) = - ExpG(j,i)
    enddo
    do j=i+1,3
      GMat(i,j) = - ExpG(i,j)
    enddo
  enddo

  do i=1,3
    do j=1,i-1
      GMat(i,i) = GMat(i,i) + ExpG(j,i)
    enddo
    do j=i+1,3
      GMat(i,i) = GMat(i,i) + ExpG(i,j)
    enddo
  enddo

!  print*,'G'
!  call matout(3,3,GMat)

! Form the Y and Y^2 matrices Eq. (16b)

  do i=1,3
    do j=i+1,3
      do k=1,3
        CenterY(i,j,k) = CenterZ(i,k) - CenterZ(j,k)
        Y2Mat(i,j) = Y2Mat(i,j) + CenterY(i,j,k)**2
      enddo
      YMat(i,j) = sqrt(Y2Mat(i,j))
    enddo
  enddo

!  print*,'Y'
!  call matout(3,3,YMat)

!  print*,'Y2'
!  call matout(3,3,Y2Mat)

! Form the delta0 and delta1 matrices Eq. (14)

  do i=1,3
    do j=1,3
      Delta0Mat(i,j) = ZetaMat(i,j)   + GMat(i,j)
      Delta1Mat(i,j) = Delta0Mat(i,j) + CMat(i,j)
    enddo
  enddo

! Form the DY and D2Y matrices

  do i=1,3
    do j=1,3
      do k=1,3
        DYMat(i,j,k) = KappaCross(i,j,k)*YMat(j,k)/ExpZ(i)
        do l=1,3
          D2YMat(i,j,k,l) = 0.5d0*KappaCross(i,k,l)*KappaCross(j,k,l)/(ExpZ(i)*ExpZ(j))
        enddo
      enddo
    enddo
  enddo

! Compute the inverse of the Delta0 and Delta1 matrices

!  InvDelta0Mat = Delta0Mat
!  InvDelta1Mat = Delta1Mat
  do i=1,3
    do j=1,3
      InvDelta0Mat(i,j) = Delta0Mat(i,j)
      InvDelta1Mat(i,j) = Delta1Mat(i,j)
    enddo
  enddo
!  call amove(3,3,Delta0Mat,InvDelta0Mat)
!  call amove(3,3,Delta1Mat,InvDelta1Mat)

  call CalcInv3(InvDelta0Mat,delta0)
  call CalcInv3(InvDelta1Mat,delta1)

!  call matout(3,3,InvDelta0Mat)
!  call matout(3,3,InvDelta1Mat)
!  print*, 'delta0,delta1 = ',delta0,delta1

! Form the Delta matrix Eq. (16a)

  do i=1,3
    do j=1,3
      do k=1,3
        do l=1,3
          D0Mat(i,j) = D0Mat(i,k) + ZetaMat(i,k)*InvDelta0Mat(k,l)*ZetaMat(l,j)
          D1Mat(i,j) = D1Mat(i,k) + ZetaMat(i,k)*InvDelta1Mat(k,l)*ZetaMat(l,j)
        enddo
      enddo
    enddo
  enddo

! Form the derivative matrices

  do i=1,3
    call CalcTrAB(3,D0Mat,D2YMat,DY0(i))
    call CalcTrAB(3,D1Mat,D2YMat,DY1(i))
    do j=1,3
      call CalcTrAB(3,D0Mat,D2YMat,D2Y0(i,j))
      call CalcTrAB(3,D1Mat,D2YMat,D2Y1(i,j))
    enddo
  enddo

! Compute Y0 and Y1

  call CalcTrAB(3,D0Mat,Y2Mat,Y0)
  call CalcTrAB(3,D1Mat,Y2Mat,Y1)

end subroutine FormVRR3e
