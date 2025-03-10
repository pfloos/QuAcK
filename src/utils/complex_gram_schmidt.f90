subroutine complex_gram_schmidt(N, vectors)
  
  ! Input variables
  implicit none
  integer, intent(in) :: N
  complex*16, intent(inout) :: vectors(N, N)

  ! Local variables
  integer :: i, j
  complex*16 :: Mtmp(N,N)
  complex*16 :: proj
  complex*16 :: norm

  ! Copy the input matrix to a temporary matrix
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
