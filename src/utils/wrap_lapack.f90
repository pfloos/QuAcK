subroutine diagonalize_general_matrix_LR(N,A,WR,VL,VR)

! Diagonalize a non-symmetric square matrix

  implicit none

! Input variables

  integer,intent(in)            :: N
  double precision,intent(inout):: A(N,N)
  double precision,intent(out)  :: VL(N,N)
  double precision,intent(out)  :: VR(N,N)
  double precision,intent(out)  :: WR(N)

! Local variables

  integer                       :: i
  double precision              :: tmp
  integer                       :: lwork,info
  double precision,allocatable  :: work(:),WI(:)

! Memory allocation

  allocate(work(1),WI(N))
  
  lwork = -1
  call dgeev('V','V',N,A,N,WR,WI,VL,N,VR,N,work,lwork,info)
  lwork = int(work(1))

  deallocate(work)

  allocate(work(lwork))

  call dgeev('V','V',N,A,N,WR,WI,VL,N,VR,N,work,lwork,info)

  do i=1,N
    tmp = dot_product(vl(:,i),vr(:,i))
    vl(:,i) = vl(:,i)/tmp
  end do

  call matout(N,N,matmul(transpose(VL),VR))

  deallocate(work,WI)

  if(info /= 0) then 
    print*,'Problem in diagonalize_general_matrix (dgeev)!!'
  end if

end subroutine 

subroutine diagonalize_general_matrix(N,A,WR,VR)

! Diagonalize a non-symmetric square matrix

  implicit none

! Input variables

  integer :: i,j,k
  integer,intent(in)            :: N
  double precision,intent(inout):: A(N,N)
  double precision,intent(out)  :: VR(N,N)
  double precision,intent(out)  :: WR(N)

! Local variables

  integer                       :: lwork,info
  double precision,allocatable  :: work(:),WI(:),VL(:,:)

! Memory allocation

  allocate(work(1),WI(N),VL(N,N))
  
  lwork = -1
  call dgeev('V','V',N,A,N,WR,WI,VL,N,VR,N,work,lwork,info)
  lwork = int(work(1))

  deallocate(work)

  allocate(work(lwork))

  call dgeev('V','V',N,A,N,WR,WI,VL,N,VR,N,work,lwork,info)

  deallocate(work,WI,VL)

  if(info /= 0) then 
    print*,'Problem in diagonalize_general_matrix (dgeev)!!'
  end if

end subroutine 
 
subroutine diagonalize_matrix(N,A,e)

! Diagonalize a square matrix

  implicit none

! Input variables

  integer,intent(in)            :: N
  double precision,intent(inout):: A(N,N)
  double precision,intent(out)  :: e(N)

! Local variables

  integer                       :: lwork,info
  double precision,allocatable  :: work(:)

! Memory allocation

  allocate(work(1))
  
  lwork = -1
  call dsyev('V','U',N,A,N,e,work,lwork,info)
  lwork = int(work(1))

  deallocate(work)

  allocate(work(lwork))

  call dsyev('V','U',N,A,N,e,work,lwork,info)

  deallocate(work)
 
  if(info /= 0) then 
    print*,'Problem in diagonalize_matrix (dsyev)!!'
  end if
  
end subroutine 

subroutine diagonalize_hmatrix(N,A,e)

! Diagonalize a complex square hermitian matrix

  implicit none

! Input variables

  integer,intent(in)            :: N
  complex*16,intent(inout)      :: A(N,N)
  double precision,intent(out)  :: e(N)

! Local variables

  integer                       :: lwork,info
  complex*16,allocatable        :: work(:)

! Memory allocation

  allocate(work(1))
  
  lwork = -1
  call zheev('V','U',N,A,N,e,work,lwork,info)
  lwork = int(work(1))

  deallocate(work)

  allocate(work(lwork))

  call zheev('V','U',N,A,N,e,work,lwork,info)

  deallocate(work)
 
  if(info /= 0) then 
    print*,'Problem in diagonalize_matrix (dsyev)!!'
  end if
  
end subroutine 

subroutine complex_diagonalize_matrix(N,A,e)

! Diagonalize a general complex matrix

  implicit none

! Input variables

  integer,intent(in)            :: N
  complex*16,intent(inout)      :: A(N,N)
  complex*16,intent(out)        :: e(N)

! Local variables

  integer                       :: lwork,info
  double precision,allocatable  :: rwork(:)
  complex*16,allocatable        :: work(:)
  complex*16,allocatable        :: VL(:,:)
  complex*16,allocatable        :: VR(:,:)

! Memory allocation
  allocate(work(1),rwork(2*N),VL(1,1),VR(N,N))
  lwork = -1
  call zgeev('N','V',N,A,N,e,VL,1,VR,N,work,lwork,rwork,info)
  lwork = max(1,int(real(work(1))))
  
  deallocate(work)
  allocate(work(lwork))

  call zgeev('N','V',N,A,N,e,VL,N,VR,N,work,lwork,rwork,info)
  call complex_sort_eigenvalues(N,e,VR)
  

  deallocate(work)
  A = VR

  if(info /= 0) then 
    print*,'Problem in diagonalize_matrix (zgeev)!!'
  end if
  
end subroutine 

subroutine complex_diagonalize_matrix_without_sort(N,A,e)

! Diagonalize a general complex matrix

  implicit none

! Input variables

  integer,intent(in)            :: N
  complex*16,intent(inout)      :: A(N,N)
  complex*16,intent(out)        :: e(N)

! Local variables

  integer                       :: lwork,info
  double precision,allocatable  :: rwork(:)
  complex*16,allocatable        :: work(:)
  complex*16,allocatable        :: VL(:,:)
  complex*16,allocatable        :: VR(:,:)

! Memory allocation
  allocate(work(1),rwork(2*N),VL(1,1),VR(N,N))
  lwork = -1
  call zgeev('N','V',N,A,N,e,VL,1,VR,N,work,lwork,rwork,info)
  lwork = max(1,int(real(work(1))))
  
  deallocate(work)
  allocate(work(lwork))

  call zgeev('N','V',N,A,N,e,VL,N,VR,N,work,lwork,rwork,info)
  

  deallocate(work)
  A = VR

  if(info /= 0) then 
    print*,'Problem in diagonalize_matrix (zgeev)!!'
  end if
  
end subroutine 


subroutine svd(N,A,U,D,Vt)

  ! Compute A = U.D.Vt
  ! Dimension of A is NxN

  implicit none

  integer, intent(in)             :: N
  double precision,intent(in)     :: A(N,N)
  double precision,intent(out)    :: U(N,N)
  double precision,intent(out)    :: Vt(N,N)
  double precision,intent(out)    :: D(N)
  double precision,allocatable    :: work(:)
  integer                         :: info,lwork

  double precision,allocatable    :: scr(:,:)

  allocate (scr(N,N))

  scr(:,:) = A(:,:)

  ! Find optimal size for temporary arrays

  allocate(work(1))

  lwork = -1
  call dgesvd('A','A',N,N,scr,N,D,U,N,Vt,N,work,lwork,info)
  lwork = int(work(1))

  deallocate(work)

  allocate(work(lwork))

  call dgesvd('A','A',N,N,scr,N,D,U,N,Vt,N,work,lwork,info)

  deallocate(work,scr)

  if (info /= 0) then
    print *,  info, ': SVD failed'
    stop
  end if

end

subroutine inverse_matrix(N,A,B)

! Returns the inverse of the square matrix A in B

  implicit none

  integer,intent(in)             :: N
  double precision, intent(in)   :: A(N,N)
  double precision, intent(out)  :: B(N,N)

  integer                        :: info,lwork
  integer, allocatable           :: ipiv(:)
  double precision,allocatable   :: work(:)

  allocate (ipiv(N),work(N*N))
  lwork = size(work)

  B(1:N,1:N) = A(1:N,1:N)

  call dgetrf(N,N,B,N,ipiv,info)

  if (info /= 0) then

    print*,info
    stop 'error in inverse (dgetrf)!!'

  end if

  call dgetri(N,B,N,ipiv,work,lwork,info)

  if (info /= 0) then

    print *,  info
    stop 'error in inverse (dgetri)!!'

  end if

  deallocate(ipiv,work)

end subroutine 

subroutine pseudo_inverse_matrix(N,A,B)

! Compute Pseudo-inverse B for real symmetric Matrix A, i.e. ABA = A.
! If A = U D U^T --> B = U D^+ U^T where D^+_ii = 1/D_ii if D_ii!=0 and 0 otherwise

  implicit none

  !Input
  integer,intent(in)            :: N
  double precision,intent(in)   :: A(N,N)
  
  ! local
  integer                       :: i
  double precision,allocatable  :: U(:,:)
  double precision,allocatable  :: DPUT(:,:)
  double precision,allocatable  :: e(:)

  !Output
  double precision,intent(out)  :: B(N,N)

  allocate(U(N,N),e(N),DPUT(N,N))
  U(:,:) = A(:,:)
  call diagonalize_matrix(N,U,e)
  DPUT(:,:) = 0d0
  ! build D^+ U^T
  do i=1,N
    if(abs(e(i))>1e-6) then
      DPUT(i,:) = U(:,i)/e(i)     
    endif
  end do

  ! B = U D^+ U^T 
  call dgemm("N", "N", N, N, N, 1.d0, &
             U(1,1), N, DPUT(1,1), N,   &
             0.d0, B(1,1), N)

  deallocate(U,e,DPUT)

end subroutine

subroutine complex_inverse_matrix(N,A,B)

! Returns the inverse of the complex square matrix A in B

  implicit none

  integer,intent(in)             :: N
  complex*16, intent(in)         :: A(N,N)
  complex*16, intent(out)        :: B(N,N)

  integer                        :: info,lwork
  integer, allocatable           :: ipiv(:)
  complex*16,allocatable         :: work(:)
  
  allocate (ipiv(N),work(N*N))
  lwork = size(work)

  B(1:N,1:N) = A(1:N,1:N)

  call zgetrf(N,N,B,N,ipiv,info)

  if (info /= 0) then

    print*,info
    stop 'error in inverse (zgetri)!!'

  end if

  call zgetri(N,B,N,ipiv,work,lwork,info)

  if (info /= 0) then

    print *,  info
    stop 'error in inverse (zgetri)!!'

  end if

  deallocate(ipiv,work)

end subroutine 
subroutine linear_solve(N,A,b,x,rcond)

! Solve the linear system A.x = b where A is a NxN matrix
! and x and x are vectors of size N

  implicit none

  integer,intent(in)             :: N
  double precision,intent(out)   :: A(N,N),b(N),rcond
  double precision,intent(out)   :: x(N)

  integer                        :: info,lwork
  double precision               :: ferr,berr
  integer,allocatable            :: ipiv(:),iwork(:)
  double precision,allocatable   :: AF(:,:),work(:)

  ! Find optimal size for temporary arrays

  allocate(work(1))
  allocate(AF(N,N),ipiv(N),iwork(N))

  lwork = -1
  call dsysvx('N','U',N,1,A,N,AF,N,ipiv,b,N,x,N,rcond,ferr,berr,work,lwork,iwork,info)
  lwork = int(work(1))

  deallocate(work)

  allocate(work(lwork))

  call dsysvx('N','U',N,1,A,N,AF,N,ipiv,b,N,x,N,rcond,ferr,berr,work,lwork,iwork,info)

! if (info /= 0) then

!   print *,  info
!   stop 'error in linear_solve (dsysvx)!!'

! end if
deallocate(work,ipiv,iwork,AF)
end subroutine 

subroutine complex_linear_solve(N,A,b,x,rcond)

! Solve the linear system A.x = b where A is a NxN matrix
! and x and x are vectors of size N

  implicit none

  integer,intent(in)             :: N
  complex*16,intent(out)         :: A(N,N),b(N)
  complex*16,intent(out)         :: x(N)
  double precision,intent(out)   :: rcond

  integer                        :: info,lwork
  double precision               :: ferr,berr
  integer,allocatable            :: ipiv(:)
  double precision,allocatable   :: rwork(:)
  complex*16,allocatable         :: AF(:,:),work(:)

  ! Find optimal size for temporary arrays

  allocate(work(1))
  allocate(AF(N,N),ipiv(N),rwork(N))

  lwork = -1
  call zsysvx('N','U',N,1,A,N,AF,N,ipiv,b,N,x,N,rcond,ferr,berr,work,lwork,rwork,info)
  lwork = max(1,int(real(work(1))))

  deallocate(work)

  allocate(work(lwork))

  call zsysvx('N','U',N,1,A,N,AF,N,ipiv,b,N,x,N,rcond,ferr,berr,work,lwork,rwork,info)

! if (info /= 0) then
!
!   print *,  info
!   stop 'error in linear_solve (zsysv)!!'
!
! end if
deallocate(work,ipiv,rwork,AF)
end subroutine 

subroutine easy_linear_solve(N,A,b,x)

! Solve the linear system A.x = b where A is a NxN matrix
! and x and x are vectors of size N

  implicit none

  integer,intent(in)             :: N
  double precision,intent(in)    :: A(N,N),b(N)
  double precision,intent(out)   :: x(N)

  integer                        :: info,lwork
  integer,allocatable            :: ipiv(:)
  double precision,allocatable   :: work(:)

  allocate(ipiv(N),work(N*N))
  lwork = size(work)

  x = b

  call dsysv('U',N,1,A,N,ipiv,x,N,work,lwork,info)

  if (info /= 0) then

    print *,  info
    stop 'error in linear_solve (dsysv)!!'

  end if

end subroutine 

