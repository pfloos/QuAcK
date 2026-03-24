subroutine complex_sort_LR_eigenvalues(N, eigvals, eigvecsL,eigvecsR)
  ! Sort eigenvalues and corresponding eigenvectors wrt the real part of the eigenvalue
  implicit none
  integer, intent(in) :: N
  complex*16, intent(inout) :: eigvals(N)
  complex*16, intent(inout) :: eigvecsL(N, N)
  complex*16, intent(inout) :: eigvecsR(N, N)

  integer :: i, j
  complex*16 :: temp_eigval
  complex*16 :: temp_vec(N)
  logical :: swapped

  do i = 1, N-1
    swapped = .FALSE.
    do j = 1, N-i
      if (REAL(eigvals(j)) > REAL(eigvals(j+1))) then
        ! Swap eigenvalues
        temp_eigval = eigvals(j)
        eigvals(j) = eigvals(j+1)
        eigvals(j+1) = temp_eigval

        ! Swap corresponding eigenvectors
!        temp_vec = eigvecsL(j,:)
!        eigvecsL(j,:) = eigvecsL(j+1,:)
!        eigvecsL(j+1,:) = temp_vec
        
        temp_vec = eigvecsR(:,j)
        eigvecsR(:,j) = eigvecsR(:,j+1)
        eigvecsR(:,j+1) = temp_vec
        
        temp_vec = eigvecsL(:,j)
        eigvecsL(:,j) = eigvecsL(:,j+1)
        eigvecsL(:,j+1) = temp_vec

        swapped = .TRUE.
      end if
    end do
    ! If no swaps were made, the array is already sorted
    if (.not. swapped) exit
  end do

end subroutine 

