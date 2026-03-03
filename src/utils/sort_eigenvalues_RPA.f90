subroutine sort_eigenvalues_RPA(n, evals, evecs)
  implicit none
  integer, intent(in)              :: n
  double precision, intent(inout)  :: evals(n)
  double precision, intent(inout)  :: evecs(n, n)
  integer                          :: i, j, k, nhalf, counter
  double precision                 :: temp_val
  double precision                 :: temp_vec(n)
  double precision                 :: etanorm
  double precision                 :: threshold = 1d-8

  ! -------------------------------
  ! Step 1: initial bubble sort
  ! Sort in descending order by real part
  ! -------------------------------
  do i = 1, n-1
   do j = i+1, n

     if (evals(i) < evals(j)) then
        ! swap
        temp_val = evals(i)
        evals(i) = evals(j)
        evals(j) = temp_val

        temp_vec = evecs(:, i)
        evecs(:, i) = evecs(:, j)
        evecs(:, j) = temp_vec
     end if

   end do
end do

  ! -------------------------------
  ! Step 2: classify first half
  !  - eta-norm > 0  -> excitation
  !  - eta-norm < 0  -> swap
  ! -------------------------------
  nhalf = n/2
  counter = 0

  do i = 1, nhalf

     etanorm = sum(abs(evecs(1:nhalf,i))**2) - &
               sum(abs(evecs(nhalf+1:n,i))**2)

     if(abs(etanorm) < threshold) then
       print *, 'Mode with vanishing eta-norm has been found !'
       print *, 'Omega:', evals(i)
       print *, 'eta-norm', etanorm
     end if

     ! normal case: classify by eta-norm
     if (etanorm < -threshold) then

        j = n + 1 - i

        temp_val = evals(i)
        evals(i) = evals(j)
        evals(j) = temp_val

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
  ! Assuming after step 2, first nhalf are excitations (positive eta norm) and last nhalf are de-excitations (negative eta norm)
  
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

  ! Sort negative de-excitations descending (real value), and descding in imaginary value
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
