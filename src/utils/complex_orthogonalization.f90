subroutine complex_orthogonalize_matrix(N, vectors)
  
  ! Input variables
  implicit none
  integer, intent(in) :: N
  complex*16, intent(inout) :: vectors(N, N)

  ! Local variables
  integer :: i, j
  complex*16,allocatable :: L(:,:),Linv(:,:)
  complex*16 :: proj
  complex*16 :: norm

  ! Copy the input matrix to a temporary matrix
  allocate(L(N,N),Linv(N,N))
  L = matmul(transpose(vectors) ,vectors)
  call complex_cholesky_decomp(N,L)
  call complex_inverse_matrix(N,L,Linv)
  vectors = matmul(vectors,transpose(Linv))
  deallocate(L,Linv)
end subroutine
subroutine complex_orthonormalize(N,vectors,A)
 
  ! Orthonormalize vectors matrix, such that vectors^T A vectors = Identity
  ! A and vectors are assumed quadratic NxN matrices

  ! Input variables
  implicit none
  integer, intent(in) :: N
  complex*16, intent(inout) :: vectors(N, N)
  complex*16, intent(inout) :: A(N, N)

  ! Local variables
  integer :: i, j
  complex*16,allocatable :: L(:,:),Linv(:,:)
  complex*16 :: proj
  complex*16 :: norm

  ! Copy the input matrix to a temporary matrix
  allocate(L(N,N),Linv(N,N))
  L = matmul(matmul(transpose(vectors),A),vectors)
  call complex_cholesky_decomp(N,L)
  call complex_inverse_matrix(N,L,Linv)
  vectors = matmul(vectors,transpose(Linv))
  deallocate(L,Linv)
end subroutine
subroutine complex_normalize_RPA(nS,XYYX)
 
  ! Orthonormalize vectors matrix, such that RPA^T (1 0; 0 -1) RPA = Identity

  ! Input variables
  implicit none
  integer, intent(in) :: nS
  complex*16, intent(inout) :: XYYX(2*nS, 2*nS)

  ! Local variables
  integer :: i
  complex*16,allocatable :: A(:,:)
  
  allocate(A(2*nS,2*nS))
  A(:,:) = (0d0,0d0)
  do i=1,nS
    A(i,i) = 1
    A(i+nS,i+nS) = -1
  end do
  call complex_orthonormalize(2*nS,XYYX,A)
  deallocate(A)
end subroutine
subroutine complex_gram_schmidt(N, vectors)
  
  ! Input variables
  implicit none
  integer, intent(in) :: N
  complex*16, intent(inout) :: vectors(N, N)

  ! Local variables
  integer :: i, j
  complex*16,allocatable :: Mtmp(:,:)
  complex*16 :: proj
  complex*16 :: norm

  ! Copy the input matrix to a temporary matrix
  allocate(Mtmp(N,N))
  Mtmp(:, :) = vectors(:,:)
  
  ! Orthonormalize the vectors
  do i = 1, N

      ! Orthonormalize with respect to the previous vectors in Mtmp
      do j = 1, i-1
          ! Compute the dot product (scalar product) of vectors j and i
          call complex_dot_product(N,Mtmp(:, j), Mtmp(:, i),proj)

          ! Subtract the projection onto the j-th vector
          Mtmp(:, i) = Mtmp(:, i) - proj * Mtmp(:, j)
      end do

      ! Normalize the vector
      call complex_dot_product(N,Mtmp(:, i), Mtmp(:, i),proj)
      norm = sqrt(proj)

      if (abs(norm) > 1.0d-10) then
          ! Normalize the i-th vector and store it back in vectors
          vectors(:, i) = Mtmp(:, i) / norm
      else
          print*, "Error: Norm of eigenvector < 1e-10 !!!"
          stop
      end if
  end do
  deallocate(Mtmp)
end subroutine

subroutine complex_dot_product(N,vec1,vec2,proj)
  ! Input
  complex*16,intent(in) :: vec1(N),vec2(N)
  !Output
  complex*16, intent(out) :: proj
  ! Local variables
  integer :: i
  
  proj = 0d0
  do i=1,N
     proj = proj +vec1(i)*vec2(i)
  end do
end subroutine

subroutine complex_cholesky_decomp(n,A)
    implicit none
    integer, intent(in) :: n               ! Matrix size
    complex*16, intent(inout) :: A(n, n)   ! Output: Lower triangular Cholesky factor L
    integer :: i, j, k
    complex*16 :: s

    ! Perform Cholesky decomposition
    do i = 1, n
        do j = 1, i
            s = A(i, j)
            do k = 1, j - 1
                s = s - A(i, k) * A(j, k)
            end do

            if (i > j) then
                A(i, j) = s / A(j, j)  ! Compute lower triangular elements
            else
                A(i, i) = sqrt(s)
            end if
        end do
    end do

    ! Zero out upper triangular part
    do i = 1, n
        do j = i + 1, n
            A(i, j) = cmplx(0.0d0, 0.0d0, kind=8)
        end do
    end do
end subroutine
