subroutine complex_sort_eigenvalues_RPA(n, evals, evecs)
  
! Sorts the eigenvecs and eigenvals like in the rpa problem, i.e. first the positive ones in ascending order than the negative ones
! in descending order

  implicit none
  integer, intent(in)           :: n
  complex*16, intent(inout)     :: evals(n)
  complex*16, intent(inout)     :: evecs(n, n)
  integer                       :: i, j, k
  complex*16                    :: temp_val
  complex*16                    :: temp_vec(n)
  
  ! bubble sort (or any sorting algorithm can be applied)
  do i = 1, n-1
     do j = i+1, n
        ! sorting condition: first ascending for positives, then by absolute value for negatives
        if ((real(evals(i)) > 0 .and. real(evals(j)) > 0 .and. real(evals(i)) > real(evals(j))) .or. &
            (real(evals(i)) < 0 .and. real(evals(j)) < 0 .and. abs(real(evals(i))) > abs(real(evals(j)))) .or. &
            (real(evals(i)) < 0 .and. real(evals(j)) > 0)) then
           
           ! swap eigenvalues
           temp_val = evals(i)
           evals(i) = evals(j)
           evals(j) = temp_val
           
           ! swap corresponding eigenvectors
           temp_vec = evecs(:, i)
           evecs(:, i) = evecs(:, j)
           evecs(:, j) = temp_vec
        end if
     end do
  end do
  
end subroutine

