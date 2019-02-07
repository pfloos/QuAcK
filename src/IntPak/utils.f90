!------------------------------------------------------------------------

function KroneckerDelta(i,j) result(delta)

! Kronecker Delta

  implicit none

! Input variables

  integer,intent(in)            :: i,j


! Output variables

  integer                       :: delta

  if(i == j) then
    delta = 1
  else
    delta = 0
  endif

end function KroneckerDelta

!------------------------------------------------------------------------

subroutine matout(m,n,A)

! Print the MxN array A

  implicit none

  integer,parameter             :: ncol = 5
  double precision,parameter    :: small = 1d-10
  integer,intent(in)            :: m,n
  double precision,intent(in)   :: A(m,n)
  double precision              :: B(ncol)
  integer                       :: ilower,iupper,num,i,j
  
  do ilower=1,n,ncol
    iupper = min(ilower + ncol - 1,n)
    num = iupper - ilower + 1
    write(*,'(3X,10(7X,I6))') (j,j=ilower,iupper)
    do i=1,m
      do j=ilower,iupper
        B(j-ilower+1) = A(i,j)
      enddo
      do j=1,num
        if(abs(B(j)) < small) B(j) = 0d0
      enddo
      write(*,'(I7,10F15.8)') i,(B(j),j=1,num)
    enddo
  enddo

end subroutine matout

!------------------------------------------------------------------------

subroutine CalcTrAB(n,A,B,Tr)

! Calculate the trace of the square matrix A.B

  implicit none

! Input variables

  integer,intent(in)            :: n
  double precision,intent(in)   :: A(n,n),B(n,n)

! Local variables

  integer                       :: i,j

! Output variables

  double precision,intent(out)  :: Tr

  Tr = 0d0
  do i=1,n
    do j=1,n
      Tr = Tr + A(i,j)*B(j,i)
    enddo
  enddo

end subroutine CalcTrAB

!------------------------------------------------------------------------

function EpsilonSwitch(i,j) result(delta)

! Epsilon function 

  implicit none

! Input variables

  integer,intent(in)            :: i,j
  integer                       :: delta

  if(i <= j) then
    delta = 1
  else
    delta = -1
  endif

end function EpsilonSwitch

!------------------------------------------------------------------------

function KappaCross(i,j,k) result(kappa)

! kappa(i,j,k) = eps(i,j) delta(i,k) - eps(k,i) delta(i,j)

  implicit none

! Input variables

  integer,intent(in)            :: i,j,k

! Local variables 

  integer                       :: EpsilonSwitch,KroneckerDelta
  double precision              :: kappa

  kappa = dble(EpsilonSwitch(i,j)*KroneckerDelta(i,k) - EpsilonSwitch(k,i)*KroneckerDelta(i,j))

end function KappaCross

!------------------------------------------------------------------------

subroutine CalcInv3(a,det)

! Calculate the inverse and the determinant of a 3x3 matrix

  implicit none

  double precision,intent(inout)  :: a(3,3)
  double precision, intent(inout) :: det
  double precision                :: b(3,3)
  integer                         :: i,j

  det = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
      - a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
      + a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))

  do i=1,3
    b(i,1) = a(i,1)
    b(i,2) = a(i,2)
    b(i,3) = a(i,3)
  enddo

  a(1,1) = b(2,2)*b(3,3) - b(2,3)*b(3,2)
  a(2,1) = b(2,3)*b(3,1) - b(2,1)*b(3,3)
  a(3,1) = b(2,1)*b(3,2) - b(2,2)*b(3,1)

  a(1,2) = b(1,3)*b(3,2) - b(1,2)*b(3,3)
  a(2,2) = b(1,1)*b(3,3) - b(1,3)*b(3,1)
  a(3,2) = b(1,2)*b(3,1) - b(1,1)*b(3,2)

  a(1,3) = b(1,2)*b(2,3) - b(1,3)*b(2,2)
  a(2,3) = b(1,3)*b(2,1) - b(1,1)*b(2,3)
  a(3,3) = b(1,1)*b(2,2) - b(1,2)*b(2,1)

  do i=1,3
    do j=1,3
      a(i,j) = a(i,j)/det
    enddo
  enddo

end subroutine CalcInv3

!------------------------------------------------------------------------

subroutine CalcInv4(a,det)

  implicit none

  double precision,intent(inout) :: a(4,4)
  double precision,intent(inout) :: det
  double precision               :: b(4,4)
  integer                        :: i,j

  det = a(1,1)*(a(2,2)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))  &
               -a(2,3)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))  &
               +a(2,4)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))) &
      - a(1,2)*(a(2,1)*(a(3,3)*a(4,4)-a(3,4)*a(4,3))  &
               -a(2,3)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))  &
               +a(2,4)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))) &
      + a(1,3)*(a(2,1)*(a(3,2)*a(4,4)-a(3,4)*a(4,2))  &
               -a(2,2)*(a(3,1)*a(4,4)-a(3,4)*a(4,1))  &
               +a(2,4)*(a(3,1)*a(4,2)-a(3,2)*a(4,1))) &
      - a(1,4)*(a(2,1)*(a(3,2)*a(4,3)-a(3,3)*a(4,2))  &
               -a(2,2)*(a(3,1)*a(4,3)-a(3,3)*a(4,1))  &
               +a(2,3)*(a(3,1)*a(4,2)-a(3,2)*a(4,1)))
  do i=1,4
    b(1,i) = a(1,i)
    b(2,i) = a(2,i)
    b(3,i) = a(3,i)
    b(4,i) = a(4,i)
  enddo

  a(1,1) =  b(2,2)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))-b(2,3)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))+b(2,4)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))
  a(2,1) = -b(2,1)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))+b(2,3)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))-b(2,4)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))
  a(3,1) =  b(2,1)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))-b(2,2)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))+b(2,4)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))
  a(4,1) = -b(2,1)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))+b(2,2)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))-b(2,3)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))

  a(1,2) = -b(1,2)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))+b(1,3)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))-b(1,4)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))
  a(2,2) =  b(1,1)*(b(3,3)*b(4,4)-b(3,4)*b(4,3))-b(1,3)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))+b(1,4)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))
  a(3,2) = -b(1,1)*(b(3,2)*b(4,4)-b(3,4)*b(4,2))+b(1,2)*(b(3,1)*b(4,4)-b(3,4)*b(4,1))-b(1,4)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))
  a(4,2) =  b(1,1)*(b(3,2)*b(4,3)-b(3,3)*b(4,2))-b(1,2)*(b(3,1)*b(4,3)-b(3,3)*b(4,1))+b(1,3)*(b(3,1)*b(4,2)-b(3,2)*b(4,1))

  a(1,3) =  b(1,2)*(b(2,3)*b(4,4)-b(2,4)*b(4,3))-b(1,3)*(b(2,2)*b(4,4)-b(2,4)*b(4,2))+b(1,4)*(b(2,2)*b(4,3)-b(2,3)*b(4,2))
  a(2,3) = -b(1,1)*(b(2,3)*b(4,4)-b(2,4)*b(4,3))+b(1,3)*(b(2,1)*b(4,4)-b(2,4)*b(4,1))-b(1,4)*(b(2,1)*b(4,3)-b(2,3)*b(4,1))
  a(3,3) =  b(1,1)*(b(2,2)*b(4,4)-b(2,4)*b(4,2))-b(1,2)*(b(2,1)*b(4,4)-b(2,4)*b(4,1))+b(1,4)*(b(2,1)*b(4,2)-b(2,2)*b(4,1))
  a(4,3) = -b(1,1)*(b(2,2)*b(4,3)-b(2,3)*b(4,2))+b(1,2)*(b(2,1)*b(4,3)-b(2,3)*b(4,1))-b(1,3)*(b(2,1)*b(4,2)-b(2,2)*b(4,1))

  a(1,4) = -b(1,2)*(b(2,3)*b(3,4)-b(2,4)*b(3,3))+b(1,3)*(b(2,2)*b(3,4)-b(2,4)*b(3,2))-b(1,4)*(b(2,2)*b(3,3)-b(2,3)*b(3,2))
  a(2,4) =  b(1,1)*(b(2,3)*b(3,4)-b(2,4)*b(3,3))-b(1,3)*(b(2,1)*b(3,4)-b(2,4)*b(3,1))+b(1,4)*(b(2,1)*b(3,3)-b(2,3)*b(3,1))
  a(3,4) = -b(1,1)*(b(2,2)*b(3,4)-b(2,4)*b(3,2))+b(1,2)*(b(2,1)*b(3,4)-b(2,4)*b(3,1))-b(1,4)*(b(2,1)*b(3,2)-b(2,2)*b(3,1))
  a(4,4) =  b(1,1)*(b(2,2)*b(3,3)-b(2,3)*b(3,2))-b(1,2)*(b(2,1)*b(3,3)-b(2,3)*b(3,1))+b(1,3)*(b(2,1)*b(3,2)-b(2,2)*b(3,1))

  do i=1,4
    do j=1,4
      a(i,j) = a(i,j)/det
    enddo
  enddo

end subroutine CalcInv4


!double precision function binom(i,j)
!  implicit none
!  integer,intent(in)             :: i,j
!  double precision               :: logfact
!  integer, save                  :: ifirst
!  double precision, save         :: memo(0:15,0:15)
!  integer                        :: k,l
!  if (ifirst == 0) then
!    ifirst = 1
!    do k=0,15
!      do l=0,15
!        memo(k,l) = dexp( logfact(k)-logfact(l)-logfact(k-l) )
!      enddo
!    enddo
!  endif
!  if ( (i<=15).and.(j<=15) ) then
!    binom = memo(i,j)
!  else
!    binom = dexp( logfact(i)-logfact(j)-logfact(i-j) )
!  endif
!end
!
!double precision function fact(n)
!  implicit none
!  integer                        :: n
!  double precision, save         :: memo(1:100)
!  integer, save                  :: memomax = 1
!
!  if (n<=memomax) then
!    if (n<2) then
!      fact = 1.d0
!    else
!      fact = memo(n)
!    endif
!    return
!  endif
!
!  integer                        :: i
!  memo(1) = 1.d0
!  do i=memomax+1,min(n,100)
!    memo(i) = memo(i-1)*dble(i)
!  enddo
!  memomax = min(n,100)
!  double precision :: logfact
!  fact = dexp(logfact(n))
!end function
!
!double precision function logfact(n)
!  implicit none
!  integer                        :: n
!  double precision, save         :: memo(1:100)
!  integer, save                  :: memomax = 1
!
!  if (n<=memomax) then
!    if (n<2) then
!      logfact = 0.d0
!    else
!      logfact = memo(n)
!    endif
!    return
!  endif
!
!  integer                        :: i
!  memo(1) = 0.d0
!  do i=memomax+1,min(n,100)
!    memo(i) = memo(i-1)+dlog(dble(i))
!  enddo
!  memomax = min(n,100)
!  logfact = memo(memomax)
!  do i=101,n
!    logfact += dlog(dble(i))
!  enddo
!end function
!
!double precision function dble_fact(n)
!  implicit none
!  integer :: n
!  double precision :: dble_fact_even, dble_fact_odd
!
!  dble_fact = 1.d0
!
!  if(n.lt.0) return
!
!  if(iand(n,1).eq.0)then
!    dble_fact = dble_fact_even(n)
!  else
!    dble_fact= dble_fact_odd(n)
!  endif
!
!end function
!
!double precision function dble_fact_even(n) result(fact2)
!  implicit none
!  integer                        :: n,k
!  double precision, save         :: memo(0:100)
!  integer, save                  :: memomax = 0
!  double precision               :: prod
!
!
!  if (n <= memomax) then
!    if (n < 2) then
!      fact2 = 1.d0
!    else
!      fact2 = memo(n)
!    endif
!    return
!  endif
!
!  integer                        :: i
!  memo(0)=1.d0
!  memo(1)=1.d0
!  do i=memomax+2,min(n,100),2
!    memo(i) = memo(i-2)* dble(i)
!  enddo
!  memomax = min(n,100)
!  fact2 = memo(memomax)
!
!  if (n > 100) then
!    double precision :: dble_logfact
!    fact2 = dexp(dble_logfact(n))
!  endif
!
!end function
!
!double precision function dble_fact_odd(n) result(fact2)
!  implicit none
!  integer                        :: n
!  double precision, save         :: memo(1:100)
!  integer, save                  :: memomax = 1
!
!  if (n<=memomax) then
!    if (n<3) then
!      fact2 = 1.d0
!    else
!      fact2 = memo(n)
!    endif
!    return
!  endif
!
!  integer                        :: i
!  memo(1) = 1.d0
!  do i=memomax+2,min(n,99),2
!    memo(i) = memo(i-2)* dble(i)
!  enddo
!  memomax = min(n,99)
!  fact2 = memo(memomax)
!
!  if (n > 99) then
!    double precision :: dble_logfact
!    fact2 = dexp(dble_logfact(n))
!  endif
!
!end function

