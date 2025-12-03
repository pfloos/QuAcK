
!------------------------------------------------------------------------
function SRG_reg(s,x) result(val)

! SRG regulization scheme

  implicit none

! Input variables

  double precision,intent(in)   :: s
  double precision,intent(in)   :: x

! Local variables

  double precision,parameter    :: thres = 1d-15

! Output variables

  double precision              :: val

  if(abs(x) < thres) then
    val = 0d0
  else
    val = (1d0 - exp(-2d0*s*x*x))/x
  end if

end function 

!------------------------------------------------------------------------
function Heaviside_step(x) result(val)

! Kronecker Delta

  implicit none

! Input variables

  double precision,intent(in)   :: x

! Output variables

  double precision              :: val

  if(x >=0d0 ) then
    val = 1d0
  else
    val = 0d0
  end if

end function 

function Kronecker_delta(i,j) result(delta)

! Kronecker Delta

  implicit none

! Input variables

  integer,intent(in)            :: i,j

! Output variables

  double precision              :: delta

  if(i == j) then
    delta = 1d0
  else
    delta = 0d0
  end if

end function 

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
  end if

end function 

!------------------------------------------------------------------------
subroutine diagonal_matrix(N,D,A)

! Construct diagonal matrix A from vector D

  implicit none

  integer,intent(in)            :: N
  double precision,intent(in)   :: D(N)
  double precision,intent(out)  :: A(N,N)

  integer                       :: i

  A(:,:) = 0d0
  do i=1,N
    A(i,i) = D(i)
  end do

end subroutine

!------------------------------------------------------------------------
subroutine matrix_exponential(N, A, ExpA)

! Compute Exp(A), where A is skew-symmetric

  implicit none

  integer,intent(in)            :: N
  integer                       :: i
  double precision,intent(in)   :: A(N,N)
  double precision,allocatable  :: W(:,:)
  double precision,allocatable  :: tau(:)
  double precision,allocatable  :: t(:,:)
  double precision,intent(out)  :: ExpA(N,N)

! Memory allocation

  allocate(W(N,N), tau(N), t(N,N))

! Initialize

  ExpA(:,:) = 0d0

! Diagonalize

  W(:,:) = - matmul(A, A)
  call diagonalize_matrix(N, W, tau)

  do i=1,N
    tau(i) = max(abs(tau(i)),1d-14)
  end do
  tau(:) = sqrt(abs(tau(:)))

! Construct cos part

  call diagonal_matrix(N, cos(tau), t)
  t(:,:) = matmul(t, transpose(W))
  ExpA(:,:) = ExpA(:,:) + matmul(W, t) 

! Construct sin part

  call diagonal_matrix(N, sin(tau)/tau, t)
  t(:,:) = matmul(t, transpose(W))
  t(:,:) = matmul(t, A)
  ExpA(:,:) = ExpA(:,:) + matmul(W, t)

  deallocate(W, tau, t)

end subroutine

!------------------------------------------------------------------------
subroutine matrix_exp(N, A, ExpA)

! Compute Exp(A), where A is skew symmetric, by using iA = U lambda U^H  and hermicity of iA

  implicit none

  integer,intent(in)            :: N
  double precision,intent(in)   :: A(N,N)
  integer                       :: j
  double precision,allocatable  :: lambda(:)
  complex*16                    :: alpha
  complex*16,allocatable        :: U(:,:)
  complex*16,allocatable        :: B(:,:)
  complex*16,allocatable        :: cExpA(:,:)
  double precision,intent(out)  :: ExpA(N,N)

! Memory allocation

  allocate(U(N,N),cExpA(N,N),lambda(N),B(N,N))

! Initialize

  ExpA(:,:) = 0d0

! Diagonalize

  U = cmplx(0d0,1d0,kind=8)*A
  call diagonalize_hmatrix(N,U,lambda)
  B = U
  do j = 1, N
    alpha = cmplx(cos(-lambda(j)), sin(-lambda(j)), kind=8)  ! exp(-i*Lambda_j)
    B(:,j) = B(:,j)*alpha
  end do
  ! ExpA = B * U^H
  call zgemm('N','C', N, N, N, cmplx(1.0d0,0.0d0,kind=8), B, N, U, N, cmplx(0.0d0,0.0d0,kind=8), cExpA, N)
  if(maxval(aimag(cExpA))>1e-14) write(*,*) "Problem in exp(A), exp(A) is not real !!!"
  ExpA = real(cExpA)

  deallocate(U,cExpA,B,lambda)

end subroutine

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
      end do
      do j=1,num
        if(abs(B(j)) < small) B(j) = 0d0
      end do
      write(*,'(I7,10F15.8)') i,(B(j),j=1,num)
    end do
  end do

end subroutine 

!------------------------------------------------------------------------
subroutine vecout(m,A)

! Print the N vector A

  implicit none

  integer,intent(in)            :: m
  double precision,intent(in)   :: A(m)
 
  call matout(m,1,A) 

end subroutine 

!------------------------------------------------------------------------
subroutine trace_vector(n,v,Tr)

! Calculate the trace of the vector v of length n
!!! Please use the intrinsic fortran sum()  !!!

  implicit none

! Input variables

  integer,intent(in)            :: n
  double precision,intent(in)   :: v(n)

! Local variables

  integer                       :: i

! Output variables

  double precision,intent(out)  :: Tr

  Tr = 0d0
  do i=1,n
    Tr = Tr + v(i)
  end do

end subroutine 

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
  end do

end function 

function complex_trace_matrix(n,A) result(Tr)

! Calculate the trace of the square matrix A

  implicit none

! Input variables

  integer,intent(in)            :: n
  complex*16,intent(in)         :: A(n,n)

! Local variables

  integer                       :: i

! Output variables

  complex*16                    :: Tr

  Tr = (0d0,0d0)
  do i=1,n
    Tr = Tr + A(i,i)
  end do

end function 
!------------------------------------------------------------------------
subroutine compute_error(nData,Mean,Var,Error)

! Calculate the statistical error

  implicit none

! Input variables

  double precision,intent(in)   :: nData,Mean(3)

! Output variables

  double precision,intent(out)  :: Error(3)
  double precision,intent(inout):: Var(3)
  
  Error = sqrt((Var-Mean**2/nData)/nData/(nData-1d0))

end subroutine 

!------------------------------------------------------------------------
subroutine identity_matrix(N,A)

! Set the matrix A to the identity matrix

  implicit none

! Input variables

  integer,intent(in)            :: N

! Local viaruabkes

  integer                       :: i

! Output variables

  double precision,intent(out)  :: A(N,N)
  
  A = 0d0

  do i=1,N
    A(i,i) = 1d0
  end do
     
end subroutine 

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
    end do
    A(i,1) = b(i)
  end do

end subroutine 

subroutine complex_prepend(N,M,A,b)

! Prepend the vector b of size N into the matrix A of size NxM

  implicit none

! Input variables

  integer,intent(in)            :: N,M
  complex*16,intent(in)         :: b(N)

! Local variables

  integer                       :: i,j

! Output variables

  complex*16,intent(out)        :: A(N,M)


! print*,'b in append'
! call complex_matout(N,1,b)

  do i=1,N
    do j=M-1,1,-1
      A(i,j+1) = A(i,j)
    end do
    A(i,1) = b(i)
  end do

end subroutine


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
    end do
    A(i,M) = b(i)
  end do

end subroutine 

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
      end do
    end do
  end do

end subroutine 

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

  double precision, allocatable :: tmp(:,:)

  allocate(tmp(N,N))
  !$OMP PARALLEL DEFAULT(NONE) PRIVATE(i, j) SHARED(N, A, D, tmp)
  !$OMP DO
  do i = 1, N
    do j = 1, N
      tmp(i,j) = D(i) * A(j,i)
    enddo
  enddo
  !$OMP END DO
  !$OMP END PARALLEL
  call dgemm("N", "N", N, N, N, 1.d0, A(1,1), N, tmp(1,1), N, 0.d0, B(1,1), N)
  deallocate(tmp)

!  B = 0d0
!  do i=1,N
!    do j=1,N
!      do k=1,N
!        B(i,k) = B(i,k) + A(i,j)*D(j)*A(k,j)
!      end do
!    end do
!  end do

end subroutine 
!------------------------------------------------------------------------
subroutine DA(N,D,A)

! Perform A <- D.A where A is a NxN matrix and D is a diagonal matrix given 
! as a vector of length N

  implicit none

  integer,intent(in)            :: N
  integer                       :: i,j,k
  double precision,intent(in)   :: D(N)
  double precision,intent(inout):: A(N,N)

  do i=1,N
    do j=1,N
      A(i,j) = D(i)*A(i,j)
    end do
  end do

end subroutine 

!------------------------------------------------------------------------
subroutine AD(N,A,D)

! Perform A <- A.D where A is a NxN matrix and D is a diagonal matrix given 
! as a vector of length N

  implicit none

  integer,intent(in)            :: N
  integer                       :: i,j,k
  double precision,intent(in)   :: D(N)
  double precision,intent(inout):: A(N,N)

  do i=1,N
    do j=1,N
      A(i,j) = A(i,j)*D(j)
    end do
  end do

end subroutine 

!------------------------------------------------------------------------
subroutine print_warning(message)

! Print warning

  implicit none

  character(len=*),intent(in)            :: message

  write(*,*) message

end subroutine 

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
    end do
  end do

end subroutine 

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
  end if

end function 

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

end function 

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
  end do

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
    end do
  end do

end subroutine 

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
  end do

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
    end do
  end do

end subroutine 

subroutine wall_time(t)
  implicit none
  double precision, intent(out)  :: t
  integer*8                        :: c
  integer*8, save                  :: rate = 0
  if (rate == 0) then
    CALL SYSTEM_CLOCK(count_rate=rate)
  end if
  CALL SYSTEM_CLOCK(count=c)
  t = dble(c)/dble(rate)
end subroutine 

! ---
! Compute <u|u>
double precision function u_dot_u(u, sze)

  implicit none

  integer,          intent(in) :: sze
  double precision, intent(in) :: u(sze)

  double precision, external   :: ddot

  !DIR$ FORCEINLINE
  u_dot_u = ddot(sze,u,1,u,1)

end

! ---

! Orthogonalization using Q.R factorization
subroutine ortho_qr(A, LDA, m, n)

  ! A : matrix to orthogonalize
  ! LDA : leftmost dimension of A
  ! m : Number of rows of A
  ! n : Number of columns of A
  ! /!\ int(WORK(1)) becomes negative when WORK(1) > 2147483648

  integer,          intent(in)    :: m, n, LDA
  double precision, intent(inout) :: A(LDA,n)

  integer                        :: LWORK, INFO
  double precision, allocatable  :: TAU(:), WORK(:)

  allocate(TAU(min(m, n)))

  allocate(WORK(1))

  LWORK = -1

  call dgeqrf(m, n, A(1,1), LDA, TAU, WORK, LWORK, INFO)

  if(INFO .lt. 0) then
    print*, ' dgeqrf failed !! ', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))

  deallocate(WORK)
  allocate(WORK(LWORK))

  call dgeqrf(m, n, A(1,1), LDA, TAU, WORK, LWORK, INFO)

  if(INFO .lt. 0) then
    print*, ' dgeqrf failed !! ', INFO
    stop
  endif

  LWORK = -1
  call dorgqr(m, n, n, A(1,1), LDA, TAU, WORK, LWORK, INFO)

  if(INFO .lt. 0) then
    print*, ' dorgqr failed !! ', INFO
    stop
  endif

  LWORK = max(n, int(WORK(1)))
 
  deallocate(WORK)
  allocate(WORK(LWORK))

  call dorgqr(m, n, n, A(1,1), LDA, TAU, WORK, LWORK, INFO)

  if(INFO .lt. 0) then
    print*, ' dorgqr failed !! ', INFO
    stop
  endif

  deallocate(WORK, TAU)

end

! ---

! Normalizes vector u
subroutine normalize(u, sze)

  implicit none

  integer,          intent(in)    :: sze
  double precision, intent(inout) :: u(sze)

  integer                         :: i
  double precision                :: d
  double precision, external      :: dnrm2

  !DIR$ FORCEINLINE
  d = dnrm2(sze, u, 1)
  if (d /= 0.d0) then
    d = 1.d0/d 
  endif
  if (d /= 1.d0) then
    !DIR$ FORCEINLINE
    call dscal(sze, d, u, 1)
  endif

end

! ---

subroutine diag_nonsym_right(n, A, A_ldim, V, V_ldim, energy, E_ldim)

  implicit none

  integer,          intent(in)  :: n, A_ldim, V_ldim, E_ldim
  double precision, intent(in)  :: A(A_ldim,n)
  double precision, intent(out) :: energy(E_ldim), V(V_ldim,n)

  character*1                   :: JOBVL, JOBVR, BALANC, SENSE
  integer                       :: i, j
  integer                       :: ILO, IHI, lda, ldvl, ldvr, LWORK, INFO
  double precision              :: ABNRM
  integer,          allocatable :: iorder(:), IWORK(:)
  double precision, allocatable :: WORK(:), SCALE_array(:), RCONDE(:), RCONDV(:)
  double precision, allocatable :: Atmp(:,:), WR(:), WI(:), VL(:,:), VR(:,:), Vtmp(:)
  double precision, allocatable :: energy_loc(:), V_loc(:,:)

  !print*, " n = ", n
  allocate(Atmp(n,n))
  allocate(WI(n))
  allocate(WR(n))
  allocate(VR(n,n))
  allocate(VL(1,1))

  do i = 1, n
    do j = 1, n
      Atmp(j,i) = A(j,i)
    enddo
  enddo

  JOBVL  = "N" ! computes the left  eigenvectors
  JOBVR  = "V" ! computes the right eigenvectors
  BALANC = "B" ! Diagonal scaling and Permutation for optimization
  SENSE  = "V" ! Determines which reciprocal condition numbers are computed
  lda  = n
  ldvr = n
  ldvl = 1

  allocate(WORK(1), SCALE_array(n), RCONDE(n), RCONDV(n), IWORK(2*n-2))

  LWORK = -1 ! to ask for the optimal size of WORK
  call dgeevx(BALANC, JOBVL, JOBVR, SENSE,                  & ! CHARACTERS
              n, Atmp, lda,                                 & ! MATRIX TO DIAGONALIZE
              WR, WI,                                       & ! REAL AND IMAGINARY PART OF EIGENVALUES
              VL, ldvl, VR, ldvr,                           & ! LEFT AND RIGHT EIGENVECTORS
              ILO, IHI, SCALE_array, ABNRM, RCONDE, RCONDV, & ! OUTPUTS OF OPTIMIZATION
              WORK, LWORK, IWORK, INFO)

  if(INFO .ne. 0) then
    print*, 'dgeevx failed !!', INFO
    stop
  endif

  LWORK = max(int(work(1)), 1) ! this is the optimal size of WORK
  deallocate(WORK)
  allocate(WORK(LWORK))
  call dgeevx(BALANC, JOBVL, JOBVR, SENSE,                  &
              n, Atmp, lda,                                 &
              WR, WI,                                       &
              VL, ldvl, VR, ldvr,                           &
              ILO, IHI, SCALE_array, ABNRM, RCONDE, RCONDV, &
              WORK, LWORK, IWORK, INFO)

  if(INFO .ne. 0) then
    print*, 'dgeevx failed !!', INFO
    stop
  endif

  deallocate(WORK, SCALE_array, RCONDE, RCONDV, IWORK)
  deallocate(VL, Atmp)


  allocate(energy_loc(n), V_loc(n,n))
  energy_loc = 0.d0
  V_loc = 0.d0

  i = 1
  do while(i .le. n)

    !print*, i, WR(i), WI(i)

    if(dabs(WI(i)) .gt. 1d-7) then

      !print*, ' Found an imaginary component to eigenvalue'
      !print*, ' Re(i) + Im(i)', i, WR(i), WI(i)

      energy_loc(i) = WR(i)
      do j = 1, n
        V_loc(j,i) = WR(i) * VR(j,i) - WI(i) * VR(j,i+1)
      enddo
      energy_loc(i+1) = WI(i)
      do j = 1, n
        V_loc(j,i+1) = WR(i) * VR(j,i+1) + WI(i) * VR(j,i)
      enddo
      i = i + 2

    else

      energy_loc(i) = WR(i)
      do j = 1, n
        V_loc(j,i) = VR(j,i)
      enddo
      i = i + 1

    endif

  enddo

  deallocate(WR, WI, VR)


  ! ordering
!  do j = 1, n
!    write(444, '(100(1X, F16.10))') (V_loc(j,i), i=1,5)
!  enddo
  allocate(iorder(n))
  do i = 1, n
    iorder(i) = i
  enddo
  call quick_sort(energy_loc, iorder, n)
  do i = 1, n
    energy(i) = energy_loc(i)
    do j = 1, n
      V(j,i) = V_loc(j,iorder(i))
    enddo
  enddo
  deallocate(iorder)
!  do j = 1, n
!    write(445, '(100(1X, F16.10))') (V_loc(j,i), i=1,5)
!  enddo
  deallocate(V_loc, energy_loc)

end

! ---

subroutine lower_case(inp_str, out_str)

  implicit none

  character(len=*),            intent(in)  :: inp_str
  character(len=len(inp_str)), intent(out) :: out_str

  integer                                  :: i
        
  out_str = inp_str

  do i = 1, len(inp_str)
    if(ichar(inp_str(i:i)) >= ichar('A') .and. ichar(inp_str(i:i)) <= ichar('Z')) then
      out_str(i:i) = char(ichar(inp_str(i:i)) + 32)
    endif
  enddo

end

! ---

integer*8 function Yoshimine_4ind(a, b, c, d)

  implicit none
  integer*8, intent(in) :: a, b, c, d
  integer*8, external :: Yoshimine_2ind

  Yoshimine_4ind = Yoshimine_2ind(Yoshimine_2ind(a, b), &
                                  Yoshimine_2ind(c, d))

  return
end

! ---

integer*8 function Yoshimine_2ind(a, b)

  implicit none
  integer*8, intent(in) :: a, b

  if(a > b) then
    !Yoshimine_2ind = (a * (a - 1)) / 2 + b
    Yoshimine_2ind = shiftr(a * (a - 1), 1) + b
  else
    !Yoshimine_2ind = (b * (b - 1)) / 2 + a
    Yoshimine_2ind = shiftr(b * (b - 1), 1) + a
  endif

  return
end

! ---

