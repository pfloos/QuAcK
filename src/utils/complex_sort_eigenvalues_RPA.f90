subroutine complex_sort_eigenvalues_RPA(n, evals, evecs)
  implicit none
  integer, intent(in)           :: n
  complex*16, intent(inout)     :: evals(n)
  complex*16, intent(inout)     :: evecs(n, n)
  integer                       :: i, j, k, nhalf, counter
  complex*16                    :: temp_val
  complex*16                    :: temp_vec(n)

  ! -------------------------------
  ! Step 1: initial bubble sort
  ! Positives ascending, negatives by abs descending
  ! -------------------------------
  do i = 1, n-1
     do j = i+1, n
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

  ! -------------------------------
  ! Step 2: swap negative excitations if X dominates Y
  ! -------------------------------
  nhalf = n/2
  counter = 0
  do i=1,nhalf
    if(sum(abs(evecs(1:nhalf,i))**2) - sum(abs(evecs(nhalf+1:n,i))**2) < 0d0) then
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
  if(counter /= 0) then
    print *, "Attention:", counter, "excitations with negative excitation energies have been found!" 
  end if

  ! -------------------------------
  ! Step 3: separate positives and negatives
  ! -------------------------------
  ! Assuming after step 2, first nhalf are excitations (positive) and last nhalf are de-excitations (negative)
  
  ! Sort positive excitations ascending
  do i = 1, nhalf-1
     do j = i+1, nhalf
        if (real(evals(i)) > real(evals(j))) then
          temp_val = evals(i)
          evals(i) = evals(j)
          evals(j) = temp_val
          temp_vec = evecs(:, i)
          evecs(:, i) = evecs(:, j)
          evecs(:, j) = temp_vec
        end if
     end do
  end do

  ! Sort negative de-excitations descending
  do i = nhalf+1, n-1
     do j = i+1, n
        if (real(evals(i)) < real(evals(j))) then
          temp_val = evals(i)
          evals(i) = evals(j)
          evals(j) = temp_val
          temp_vec = evecs(:, i)
          evecs(:, i) = evecs(:, j)
          evecs(:, j) = temp_vec
        end if
     end do
  end do

end subroutine
