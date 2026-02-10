subroutine complex_sort_eigenvalues_RPA(n, evals, evecs)
  
! Sorts the eigenvecs and eigenvals like in the rpa problem, i.e. first the positive ones in ascending order than the negative ones
! in descending order. Then it swaps the negative eigenvalues which correspond to excitations with their positive counterparts. They
! are identified by computing the norm of the corresponding X and Y, if the X contribution is dominating it is considered as an
! exciation.

  implicit none
  integer, intent(in)           :: n
  complex*16, intent(inout)     :: evals(n)
  complex*16, intent(inout)     :: evecs(n, n)
  integer                       :: i, j, k,nhalf,counter
  complex*16                    :: temp_val
  complex*16                    :: temp_vec(n)
  
  ! bubble sort 
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
 
  ! Swap the right negative excitation energies in
  nhalf = n/2
  counter = 0
  do i=1,nhalf
    if(sum(abs(evecs(1:nhalf,i))**2) - sum(abs(evecs(nhalf+1:n,i))**2)<0d0) then
    
      j = nhalf + i
      ! swap eigenvalues
      temp_val = evals(i)
      evals(i) = evals(j)
      evals(j) = temp_val
      ! swap corresponding eigenvectors
      temp_vec = evecs(:, i)
      evecs(:, i) = evecs(:, j)
      evecs(:, j) = temp_vec
      counter = counter + 1  
    end if
  end do
  if(counter/=0) then
    print *, "Attention:", counter, "excitations with negative excitation energies have been found!" 
  end if
end subroutine

