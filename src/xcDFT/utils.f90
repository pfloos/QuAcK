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
    write(*,'(3X,10(9X,I6))') (j,j=ilower,iupper)
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
function trace_matrix(n,A) result(Tr)

! Calculate the trace of the square matrix A

  implicit none

! Input variables

  integer,intent(in)            :: n
  double precision,intent(in)   :: A(n,n)

! Local variables

  integer                       :: i

! Output variables

  double precision              :: Tr

  Tr = 0d0
  do i=1,n
    Tr = Tr + A(i,i)
  enddo

end function trace_matrix
!------------------------------------------------------------------------
subroutine prepend(N,M,A,b)

! Prepend the vector b of size N into the matrix A of size NxM

  implicit none

! Input variables

  integer,intent(in)            :: N,M
  double precision,intent(in)   :: b(N)

! Local viaruabkes

  integer                       :: i,j

! Output variables

  double precision,intent(out)  :: A(N,M)


! print*,'b in append'
! call matout(N,1,b)

  do i=1,N
    do j=M-1,1,-1
      A(i,j+1) = A(i,j)
    enddo
    A(i,1) = b(i)
  enddo

end subroutine prepend
!------------------------------------------------------------------------
subroutine append(N,M,A,b)

! Append the vector b of size N into the matrix A of size NxM

  implicit none

! Input variables

  integer,intent(in)            :: N,M
  double precision,intent(in)   :: b(N)

! Local viaruabkes

  integer                       :: i,j

! Output variables

  double precision,intent(out)  :: A(N,M)

  do i=1,N
    do j=2,M
      A(i,j-1) = A(i,j)
    enddo
    A(i,M) = b(i)
  enddo

end subroutine append
!------------------------------------------------------------------------
subroutine AtDA(N,A,D,B)

! Perform B = At.D.A where A is a NxN matrix and D is a diagonal matrix given 
! as a vector of length N

  implicit none

! Input variables

  integer,intent(in)            :: N
  double precision,intent(in)   :: A(N,N),D(N)

! Local viaruabkes

  integer                       :: i,j,k

! Output variables

  double precision,intent(out)  :: B(N,N)

  B = 0d0

  do i=1,N
    do j=1,N
      do k=1,N
        B(i,k) = B(i,k) + A(j,i)*D(j)*A(j,k)
      enddo
    enddo
  enddo

end subroutine AtDA
!------------------------------------------------------------------------
subroutine ADAt(N,A,D,B)

! Perform B = A.D.At where A is a NxN matrix and D is a diagonal matrix given 
! as a vector of length N

  implicit none

! Input variables

  integer,intent(in)            :: N
  double precision,intent(in)   :: A(N,N),D(N)

! Local viaruabkes

  integer                       :: i,j,k

! Output variables

  double precision,intent(out)  :: B(N,N)

  B = 0d0

  do i=1,N
    do j=1,N
      do k=1,N
        B(i,k) = B(i,k) + A(i,j)*D(j)*A(k,j)
      enddo
    enddo
  enddo

end subroutine ADAt
!------------------------------------------------------------------------
subroutine DA(N,D,A)

! Perform A <- D.A where A is a NxN matrix and D is a diagonal matrix given 
! as a vector of length N

  implicit none

  integer,intent(in)            :: N
  integer                       :: i,j
  double precision,intent(in)   :: D(N)
  double precision,intent(inout):: A(N,N)

  do i=1,N
    do j=1,N
      A(i,j) = D(i)*A(i,j)
    enddo
  enddo

end subroutine DA

!------------------------------------------------------------------------
subroutine AD(N,A,D)

! Perform A <- A.D where A is a NxN matrix and D is a diagonal matrix given 
! as a vector of length N

  implicit none

  integer,intent(in)            :: N
  integer                       :: i,j
  double precision,intent(in)   :: D(N)
  double precision,intent(inout):: A(N,N)

  do i=1,N
    do j=1,N
      A(i,j) = A(i,j)*D(j)
    enddo
  enddo

end subroutine AD
!------------------------------------------------------------------------
recursive function fac(n) result(fact)

  implicit none
  integer :: fact
  integer, intent(in) :: n

  if (n == 0) then
     fact = 1
  else
     fact = n * fac(n-1)
  end if

end function fac

function dfac(n) result(fact)

  implicit none
  double precision :: fact
  integer, intent(in) :: n
  integer             :: fac

  fact = dble(fac(n))

end function dfac
